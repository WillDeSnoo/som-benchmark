#!/usr/bin/env python
"""Precompute C-H bond dissociation energies with ALFABET.

This is the one stage that does NOT run in the `som-bench` environment.
ALFABET depends on TensorFlow and hard-pins `scikit-learn==0.24.2`, neither of
which belongs in the main environment. It therefore runs in its own
`alfabet-env` interpreter and communicates through a CSV, which the rest of the
pipeline treats as just another precomputed input:

    conda run -n alfabet-env python scripts/07_bde.py

(The sklearn pin is only needed by `alfabet.neighbors`, the applicability-domain
check, which the prediction path never imports -- so `alfabet` is installed with
`--no-deps` against modern dependencies. See README.)

Why this is worth doing at all: ALFABET is itself a graph network, so its output
is in principle a function of the molecular graph, which the GNN ablation showed
is not enough on its own. The point is that ALFABET was trained on ~290k DFT
bond dissociations. It transfers knowledge that 679 molecules cannot teach. This
is knowledge distillation, not feature engineering.

Atom indexing: ALFABET canonicalises with `Chem.MolToSmiles` and fragments
`Chem.AddHs(mol)`. Since the curated SMILES are already RDKit-canonical and
`AddHs` appends hydrogens after the heavy atoms, heavy-atom indices are
preserved. That assumption is verified for every molecule rather than trusted.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MOLECULES = ROOT / "data" / "processed" / "molecules.csv"
OUT = ROOT / "data" / "interim" / "alfabet_bde.csv"

#: stand-in BDE for an atom with no abstractable hydrogen. Chosen above any
#: real C-H bond so that tree splits and the GNN both read it as "unreactive"
#: rather than as a very weak bond, which a 0 fill would imply.
NO_H_BDE = 140.0


def load_smiles() -> list[str]:
    if not MOLECULES.exists():
        sys.exit(f"{MOLECULES} not found - run `make curate` first")
    with MOLECULES.open(newline="") as fh:
        return [r["smiles"] for r in csv.DictReader(fh)]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0, help="only the first N (smoke test)")
    args = ap.parse_args()

    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from alfabet import model as alfabet_model

    smiles = load_smiles()
    if args.limit:
        smiles = smiles[: args.limit]
    print(f"predicting BDEs for {len(smiles)} molecules", flush=True)

    # drop_duplicates=False: symmetry-equivalent bonds must each keep a row,
    # otherwise most atoms of a symmetric molecule would come back empty.
    pred = alfabet_model.predict(
        smiles, drop_duplicates=False, batch_size=args.batch_size, verbose=True
    )
    print(f"got {len(pred)} bond predictions", flush=True)

    # ---- map bond_index (on the AddHs molecule) back to the heavy atom ----
    rows: list[dict] = []
    skipped_mismatch = 0
    for smi, grp in pred.groupby("molecule"):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        molh = Chem.AddHs(mol)
        n_heavy = mol.GetNumAtoms()

        # verify the index assumption rather than trusting it
        ok = all(
            molh.GetAtomWithIdx(i).GetSymbol() == mol.GetAtomWithIdx(i).GetSymbol()
            for i in range(n_heavy)
        )
        if not ok:
            skipped_mismatch += 1
            continue

        best: dict[int, tuple[float, float]] = {}
        for _, r in grp.iterrows():
            bi = int(r["bond_index"])
            if bi >= molh.GetNumBonds():
                continue
            bond = molh.GetBondWithIdx(bi)
            a, b = bond.GetBeginAtom(), bond.GetEndAtom()
            # keep only X-H bonds, and record against the heavy partner
            if a.GetAtomicNum() == 1 and b.GetAtomicNum() != 1:
                heavy = b.GetIdx()
            elif b.GetAtomicNum() == 1 and a.GetAtomicNum() != 1:
                heavy = a.GetIdx()
            else:
                continue
            if heavy >= n_heavy:
                continue
            bde, bdfe = float(r["bde_pred"]), float(r["bdfe_pred"])
            if bde != bde:  # NaN
                continue
            prev = best.get(heavy)
            if prev is None or bde < prev[0]:
                best[heavy] = (bde, bdfe)

        valid = bool(grp["is_valid"].all()) if "is_valid" in grp else True
        for heavy, (bde, bdfe) in sorted(best.items()):
            rows.append({
                "smiles": smi,
                "atom_idx": heavy,
                "bde": round(bde, 4),
                "bdfe": round(bdfe, 4),
                "is_valid": int(valid),
            })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["smiles", "atom_idx", "bde", "bdfe", "is_valid"]
        )
        w.writeheader()
        for r in rows:
            w.writerow(r)

    covered = len({r["smiles"] for r in rows})
    invalid = len({r["smiles"] for r in rows if not r["is_valid"]})
    print(f"\nwrote {OUT}")
    print(f"  atoms with a predicted X-H BDE : {len(rows)}")
    print(f"  molecules covered              : {covered} / {len(smiles)}")
    print(f"  molecules outside ALFABET's domain: {invalid}")
    if skipped_mismatch:
        print(f"  skipped (atom-order mismatch)  : {skipped_mismatch}")


if __name__ == "__main__":
    main()
