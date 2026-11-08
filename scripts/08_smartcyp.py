#!/usr/bin/env python
"""Run SMARTCyp over the curated molecules and record its per-atom scores.

SMARTCyp (Rydberg & Olsen) is the standard published method for this task: a
lookup of DFT-derived activation energies by SMARTS atom type, combined with a
2D accessibility term. It is the external yardstick this project otherwise
lacks -- without it, the models here are only compared against each other.

It is a Java program and is not redistributed in this repository. Build it
yourself from the CDK port (LGPL-3.0) and point this script at the jar:

    git clone https://github.com/cdk/smartcyp
    cd smartcyp && mvn -B -DskipTests package
    python scripts/08_smartcyp.py --jar smartcyp/target/smartcyp.jar

Atom indexing: SMARTCyp labels atoms `Element.N`, 1-based over the input order.
CDK parses SMILES preserving atom order, and the curated SMILES are already
RDKit-canonical, so N-1 is the index used everywhere else in this project. The
element in the label is checked against the molecule for every single row
rather than assumed -- a silent off-by-one here would make SMARTCyp look much
worse than it is, which is the most misleading way for a baseline to fail.

Only carbon, nitrogen and sulfur are scored by SMARTCyp; atoms it does not
report are written with `UNREACTIVE_SCORE`.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rdkit import Chem, RDLogger  # noqa: E402

from som.paths import INTERIM, PROCESSED  # noqa: E402

RDLogger.DisableLog("rdApp.*")

OUT = INTERIM / "smartcyp.csv"

#: SMARTCyp's own "will not react" sentinel is an energy of 999; scores for
#: such atoms land just under 1000. Unreported atoms get the same treatment.
UNREACTIVE_SCORE = 999.0

ATOM_LABEL = re.compile(r"^([A-Za-z]+)\.(\d+)$")

FIELDS = ["smiles", "atom_idx", "score", "energy", "score_2d6", "score_2c"]


def load_smiles() -> list[str]:
    path = PROCESSED / "molecules.csv"
    if not path.exists():
        sys.exit(f"{path} not found - run `make curate` first")
    with path.open(newline="") as fh:
        return [r["smiles"] for r in csv.DictReader(fh)]


def run_smartcyp(
    jar: Path, smiles: list[str], workdir: Path, stack: str
) -> Path | None:
    """Run the jar on one batch of SMILES; None if the JVM died.

    CDK's `DeduceBondSystemTool.fixAromaticBondOrders` recurses deeply on some
    fused aromatics and overflows the default JVM stack, which is why `-Xss` is
    raised here. The CDK port of SMARTCyp carries a "not tested yet" note from
    its authors, so the caller also degrades to one molecule at a time rather
    than assuming a clean run.
    """
    for f in glob.glob(str(workdir / "SMARTCyp_Results_*.csv")):
        os.remove(f)
    smi_file = workdir / "input.smi"
    smi_file.write_text("\n".join(smiles) + "\n")

    proc = subprocess.run(
        ["java", f"-Xss{stack}", "-jar", str(jar.resolve()), smi_file.name],
        cwd=workdir, capture_output=True, text=True,
    )
    hits = sorted(glob.glob(str(workdir / "SMARTCyp_Results_*.csv")))
    if proc.returncode != 0 or not hits:
        return None
    return Path(hits[-1])


def run_chunked(
    jar: Path, smiles: list[str], workdir: Path, stack: str, chunk: int
) -> tuple[list[dict], list[str]]:
    """Run in chunks, falling back to single molecules when a chunk fails."""
    rows: list[dict] = []
    failed: list[str] = []
    total_stats = {"element_mismatch": 0}

    for start in range(0, len(smiles), chunk):
        batch = smiles[start : start + chunk]
        out = run_smartcyp(jar, batch, workdir, stack)
        if out is not None:
            got, st = parse(out, batch)
            rows.extend(got)
            total_stats["element_mismatch"] += st["element_mismatch"]
            print(f"  [{start + len(batch):>4}/{len(smiles)}] ok", flush=True)
            continue

        print(f"  [{start + len(batch):>4}/{len(smiles)}] chunk failed, "
              "retrying one by one", flush=True)
        for smi in batch:
            one = run_smartcyp(jar, [smi], workdir, stack)
            if one is None:
                failed.append(smi)
                continue
            got, st = parse(one, [smi])
            rows.extend(got)
            total_stats["element_mismatch"] += st["element_mismatch"]

    if total_stats["element_mismatch"]:
        sys.exit(
            f"{total_stats['element_mismatch']} element mismatches - "
            "atom mapping is wrong, refusing to write"
        )
    return rows, failed


def _f(value: str, default: float = UNREACTIVE_SCORE) -> float:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return default
    return v if v == v else default


def parse(result_csv: Path, smiles: list[str]) -> tuple[list[dict], dict]:
    rows: list[dict] = []
    stats = {"parsed": 0, "bad_label": 0, "element_mismatch": 0,
             "out_of_range": 0, "molecules": set()}

    mols = [Chem.MolFromSmiles(s) for s in smiles]

    with result_csv.open(newline="") as fh:
        for r in csv.DictReader(fh):
            try:
                mol_i = int(r["Molecule"]) - 1
            except (TypeError, ValueError):
                continue
            if not (0 <= mol_i < len(smiles)):
                continue
            mol = mols[mol_i]
            if mol is None:
                continue

            m = ATOM_LABEL.match((r.get("Atom") or "").strip())
            if not m:
                stats["bad_label"] += 1
                continue
            element, num = m.group(1), int(m.group(2))
            idx = num - 1
            if not (0 <= idx < mol.GetNumAtoms()):
                stats["out_of_range"] += 1
                continue
            if mol.GetAtomWithIdx(idx).GetSymbol().upper() != element.upper():
                stats["element_mismatch"] += 1
                continue

            stats["parsed"] += 1
            stats["molecules"].add(mol_i)
            rows.append({
                "smiles": smiles[mol_i],
                "atom_idx": idx,
                "score": round(_f(r.get("Score")), 4),
                "energy": round(_f(r.get("Energy")), 4),
                "score_2d6": round(_f(r.get("2D6score")), 4),
                "score_2c": round(_f(r.get("2Cscore")), 4),
            })
    return rows, stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--jar", default=os.environ.get("SMARTCYP_JAR", ""),
        help="path to smartcyp.jar (or set SMARTCYP_JAR)",
    )
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--keep-workdir", action="store_true")
    ap.add_argument("--stack", default="512m", help="JVM -Xss stack size")
    ap.add_argument("--chunk", type=int, default=100)
    args = ap.parse_args()

    if not args.jar:
        sys.exit("need --jar /path/to/smartcyp.jar (see the module docstring)")
    jar = Path(args.jar)
    if not jar.exists():
        sys.exit(f"jar not found: {jar}")
    if shutil.which("java") is None:
        sys.exit("java not found on PATH")

    smiles = load_smiles()
    if args.limit:
        smiles = smiles[: args.limit]
    print(f"running SMARTCyp on {len(smiles)} molecules", flush=True)

    workdir = Path(tempfile.mkdtemp(prefix="smartcyp_"))
    try:
        rows, failed = run_chunked(
            jar, smiles, workdir, args.stack, args.chunk
        )
    finally:
        if not args.keep_workdir:
            shutil.rmtree(workdir, ignore_errors=True)
        else:
            print(f"workdir kept at {workdir}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    covered = len({r["smiles"] for r in rows})
    print(f"\nwrote {OUT}")
    print(f"  scored atoms        : {len(rows)}")
    print(f"  molecules covered   : {covered} / {len(smiles)}")
    if failed:
        print(f"  molecules SMARTCyp could not process : {len(failed)}")
        for s in failed[:10]:
            print(f"      {s}")
        if len(failed) > 10:
            print(f"      ... and {len(failed) - 10} more")
        (INTERIM / "smartcyp_failed.txt").write_text("\n".join(failed) + "\n")
        print(f"  (listed in {INTERIM / 'smartcyp_failed.txt'})")


if __name__ == "__main__":
    main()
