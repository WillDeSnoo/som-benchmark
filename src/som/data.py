"""Download and curate the CYP site-of-metabolism data.

The upstream SDFs label sites of metabolism as 1-indexed heavy-atom numbers into
the *as-distributed* atom ordering. Standardization can delete atoms (salt
stripping) and reorder them (canonicalization), so the single most important job
in this module is to carry those labels through every transformation without
silently corrupting them.

The strategy:

1. Tag every atom with its original 0-based index in an atom property before
   touching the molecule.
2. Standardize.
3. Read the tags back off the surviving atoms to learn the old -> new mapping.
4. Re-express the labels as positions in the *canonical SMILES output order*, so
   that `Chem.MolFromSmiles(canonical_smiles)` reproduces exactly the atom
   ordering the labels refer to. The curated table therefore needs nothing but a
   SMILES string and a list of integers, and the SDFs are never read again.

A round-trip assertion guards step 4 for every record; a record that fails it is
excluded rather than used.
"""

from __future__ import annotations

import csv
import urllib.request
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

from . import ISOFORMS, UPSTREAM_COMMIT, UPSTREAM_PATH, UPSTREAM_REPO
from .paths import PROCESSED, RAW

RDLogger.DisableLog("rdApp.*")

ORIG_IDX_PROP = "_somOrigIdx"

#: SOM annotation tiers present upstream. `PRIMARY_SOM` is on every record;
#: the other two are sparser. The headline task treats all three as positives
#: (the usual convention for this dataset); `primary_only` is reported as a
#: sensitivity analysis.
SOM_FIELDS = ("PRIMARY_SOM", "SECONDARY_SOM", "TERTIARY_SOM")


# --------------------------------------------------------------------------
# download
# --------------------------------------------------------------------------

def raw_url(isoform: str) -> str:
    return (
        f"https://raw.githubusercontent.com/{UPSTREAM_REPO}/"
        f"{UPSTREAM_COMMIT}/{UPSTREAM_PATH}/{isoform}.sdf"
    )


def download(force: bool = False) -> list[Path]:
    """Fetch the nine per-isoform SDFs to data/raw/ at the pinned commit."""
    RAW.mkdir(parents=True, exist_ok=True)
    out = []
    for iso in ISOFORMS:
        dest = RAW / f"{iso}.sdf"
        if dest.exists() and not force:
            print(f"  {iso}.sdf  already present ({dest.stat().st_size:,} bytes)")
        else:
            url = raw_url(iso)
            with urllib.request.urlopen(url, timeout=120) as r:
                dest.write_bytes(r.read())
            print(f"  {iso}.sdf  downloaded ({dest.stat().st_size:,} bytes)")
        out.append(dest)
    return out


# --------------------------------------------------------------------------
# property access (the upstream Citation fields are not valid UTF-8)
# --------------------------------------------------------------------------

def safe_prop(mol: Chem.Mol, key: str) -> str | None:
    """Read a mol property, tolerating the non-UTF-8 bytes in some Citations."""
    if not mol.HasProp(key):
        return None
    try:
        return mol.GetProp(key)
    except (UnicodeDecodeError, RuntimeError):
        return None


def parse_som_field(value: str | None) -> list[int]:
    """Parse a SOM property into 0-based atom indices.

    Upstream values are whitespace- or comma-separated 1-indexed atom numbers.
    Non-integer tokens are reported by returning them separately so the caller
    can log an exclusion rather than guess.
    """
    if not value:
        return []
    out = []
    for tok in value.replace(",", " ").split():
        out.append(int(tok) - 1)  # may raise ValueError; caller handles
    return out


# --------------------------------------------------------------------------
# standardization
# --------------------------------------------------------------------------

def build_standardizer():
    params = rdMolStandardize.CleanupParameters()
    largest = rdMolStandardize.LargestFragmentChooser(params)
    uncharger = rdMolStandardize.Uncharger()

    def run(mol: Chem.Mol) -> Chem.Mol:
        m = rdMolStandardize.Cleanup(mol, params)
        m = largest.choose(m)
        m = uncharger.uncharge(m)
        Chem.SanitizeMol(m)
        # Drop inherited mol-level properties (the upstream `Citation` fields
        # contain non-UTF-8 bytes that break RDKit's bulk property accessors
        # downstream). Atom-level properties, including our index tags, stay.
        for name in list(m.GetPropNames()):
            m.ClearProp(name)
        return m

    return run


def canonical_order_map(mol: Chem.Mol) -> tuple[str, dict[int, int]]:
    """Canonical SMILES plus a map from current atom index -> SMILES position.

    `_smilesAtomOutputOrder` lists, for each position in the SMILES, the atom
    index it came from. Inverting it gives what we need.
    """
    smiles = Chem.MolToSmiles(mol)
    # Read the property directly rather than via GetPropsAsDict: some upstream
    # records carry non-UTF-8 bytes in `Citation`, which makes the bulk
    # accessor raise before it ever reaches the key we want.
    raw = mol.GetProp("_smilesAtomOutputOrder")
    order = [int(x) for x in raw.strip("[] ,").replace(",", " ").split()]
    return smiles, {old: new for new, old in enumerate(order)}


# --------------------------------------------------------------------------
# curation
# --------------------------------------------------------------------------

@dataclass
class Exclusion:
    isoform: str
    record: int
    mol_id: str
    reason: str
    detail: str = ""


@dataclass
class CurationReport:
    records_read: int = 0
    exclusions: list[Exclusion] = field(default_factory=list)
    duplicate_instances_merged: int = 0
    duplicate_label_conflicts: int = 0
    labels_dropped_by_standardization: int = 0

    def count(self, reason: str) -> int:
        return sum(1 for e in self.exclusions if e.reason == reason)


def curate(
    raw_dir: Path | None = None,
    som_fields: tuple[str, ...] = SOM_FIELDS,
) -> tuple[list[dict], CurationReport]:
    """Standardize, relabel and deduplicate. Returns rows and a report.

    One row per surviving (molecule, isoform) instance:
        isoform, mol_id, smiles, n_heavy_atoms, som_primary, som_all, n_sources
    where the som_* fields are sorted 0-based atom indices into the canonical
    SMILES atom order.
    """
    raw_dir = raw_dir or RAW
    std = build_standardizer()
    report = CurationReport()

    # (smiles, isoform) -> accumulated labels
    acc: dict[tuple[str, str], dict] = {}

    for iso in ISOFORMS:
        path = raw_dir / f"{iso}.sdf"
        supplier = Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True)
        for rec_i, mol in enumerate(supplier):
            report.records_read += 1
            mol_id = "?"

            if mol is None:
                report.exclusions.append(
                    Exclusion(iso, rec_i, mol_id, "rdkit_parse_failed")
                )
                continue

            mol_id = safe_prop(mol, "ID") or mol.GetProp("_Name") or f"{iso}:{rec_i}"
            n_in = mol.GetNumAtoms()

            # --- read labels in the as-distributed ordering -----------------
            tiers: dict[str, list[int]] = {}
            bad_token = None
            for fld in som_fields:
                try:
                    tiers[fld] = parse_som_field(safe_prop(mol, fld))
                except ValueError as exc:
                    bad_token = f"{fld}: {exc}"
                    break
            if bad_token is not None:
                report.exclusions.append(
                    Exclusion(iso, rec_i, mol_id, "unparsable_som_field", bad_token)
                )
                continue

            primary_in = tiers.get("PRIMARY_SOM", [])
            all_in = sorted({i for v in tiers.values() for i in v})

            if not primary_in:
                report.exclusions.append(
                    Exclusion(iso, rec_i, mol_id, "no_primary_som")
                )
                continue

            out_of_range = [i for i in all_in if not (0 <= i < n_in)]
            if out_of_range:
                report.exclusions.append(
                    Exclusion(
                        iso, rec_i, mol_id, "som_index_out_of_range",
                        f"n_atoms={n_in} bad={[i + 1 for i in out_of_range]}",
                    )
                )
                continue

            # --- tag, standardize, recover the mapping ----------------------
            for atom in mol.GetAtoms():
                atom.SetIntProp(ORIG_IDX_PROP, atom.GetIdx())
            try:
                clean = std(mol)
            except Exception as exc:  # noqa: BLE001 - report, don't crash
                report.exclusions.append(
                    Exclusion(iso, rec_i, mol_id, "standardization_failed", str(exc))
                )
                continue

            if clean.GetNumAtoms() == 0:
                report.exclusions.append(
                    Exclusion(iso, rec_i, mol_id, "empty_after_standardization")
                )
                continue

            surviving = {}
            for atom in clean.GetAtoms():
                if atom.HasProp(ORIG_IDX_PROP):
                    surviving[atom.GetIntProp(ORIG_IDX_PROP)] = atom.GetIdx()

            if len(surviving) != clean.GetNumAtoms():
                report.exclusions.append(
                    Exclusion(
                        iso, rec_i, mol_id, "lost_atom_provenance",
                        f"tagged={len(surviving)} atoms={clean.GetNumAtoms()}",
                    )
                )
                continue

            lost = [i for i in all_in if i not in surviving]
            if lost:
                report.labels_dropped_by_standardization += 1
                report.exclusions.append(
                    Exclusion(
                        iso, rec_i, mol_id, "som_atom_removed_by_standardization",
                        f"lost={[i + 1 for i in lost]}",
                    )
                )
                continue

            # --- basic sanity filters ---------------------------------------
            if clean.GetNumAtoms() < 2:
                report.exclusions.append(
                    Exclusion(iso, rec_i, mol_id, "fewer_than_2_heavy_atoms")
                )
                continue
            if not any(a.GetSymbol() == "C" for a in clean.GetAtoms()):
                report.exclusions.append(
                    Exclusion(iso, rec_i, mol_id, "no_carbon")
                )
                continue
            if len(Chem.GetMolFrags(clean)) != 1:
                report.exclusions.append(
                    Exclusion(iso, rec_i, mol_id, "disconnected_after_standardization")
                )
                continue

            # --- re-express labels in canonical SMILES order ----------------
            smiles, to_canon = canonical_order_map(clean)
            try:
                primary_c = sorted({to_canon[surviving[i]] for i in primary_in})
                all_c = sorted({to_canon[surviving[i]] for i in all_in})
            except KeyError as exc:
                report.exclusions.append(
                    Exclusion(iso, rec_i, mol_id, "canonical_remap_failed", str(exc))
                )
                continue

            # round-trip guard: reloading the SMILES must reproduce the
            # element at every labelled position
            probe = Chem.MolFromSmiles(smiles)
            if probe is None or probe.GetNumAtoms() != clean.GetNumAtoms():
                report.exclusions.append(
                    Exclusion(iso, rec_i, mol_id, "smiles_roundtrip_failed")
                )
                continue
            inv = {v: k for k, v in to_canon.items()}
            mismatch = [
                j for j in all_c
                if probe.GetAtomWithIdx(j).GetSymbol()
                != clean.GetAtomWithIdx(inv[j]).GetSymbol()
            ]
            if mismatch:
                report.exclusions.append(
                    Exclusion(
                        iso, rec_i, mol_id, "label_roundtrip_mismatch",
                        f"positions={mismatch}",
                    )
                )
                continue

            # --- accumulate, merging duplicates -----------------------------
            key = (smiles, iso)
            if key in acc:
                prev = acc[key]
                report.duplicate_instances_merged += 1
                if set(prev["som_all"]) != set(all_c):
                    report.duplicate_label_conflicts += 1
                # union: both are reported experimental observations
                prev["som_primary"] = sorted(set(prev["som_primary"]) | set(primary_c))
                prev["som_all"] = sorted(set(prev["som_all"]) | set(all_c))
                prev["n_sources"] += 1
                prev["mol_id"] = min(prev["mol_id"], mol_id)
            else:
                acc[key] = {
                    "isoform": iso,
                    "mol_id": mol_id,
                    "smiles": smiles,
                    "n_heavy_atoms": clean.GetNumAtoms(),
                    "som_primary": primary_c,
                    "som_all": all_c,
                    "n_sources": 1,
                }

    rows = sorted(acc.values(), key=lambda r: (r["smiles"], r["isoform"]))
    return rows, report


def murcko_scaffold(smiles: str) -> str:
    """Bemis-Murcko scaffold SMILES; empty string for acyclic molecules."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception:  # noqa: BLE001
        return ""
    return Chem.MolToSmiles(scaf) if scaf.GetNumAtoms() else ""


def molecule_table(rows: list[dict]) -> list[dict]:
    """One row per unique molecule, with scaffold and the isoforms it appears in."""
    by_smiles: dict[str, dict] = {}
    isos = defaultdict(set)
    for r in rows:
        isos[r["smiles"]].add(r["isoform"])
        if r["smiles"] not in by_smiles:
            mol = Chem.MolFromSmiles(r["smiles"])
            by_smiles[r["smiles"]] = {
                "smiles": r["smiles"],
                "mol_id": r["mol_id"],
                "n_heavy_atoms": r["n_heavy_atoms"],
                "scaffold": murcko_scaffold(r["smiles"]),
                "n_rings": rdMolDescriptors.CalcNumRings(mol),
                "mw": round(rdMolDescriptors.CalcExactMolWt(mol), 3),
            }
    out = []
    for smi, rec in by_smiles.items():
        rec = dict(rec)
        rec["isoforms"] = " ".join(sorted(isos[smi]))
        rec["n_isoforms"] = len(isos[smi])
        out.append(rec)
    return sorted(out, key=lambda r: r["smiles"])


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------

def _write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(r[k]) for k in fields})


def _fmt(v):
    if isinstance(v, (list, tuple)):
        return " ".join(str(x) for x in v)
    return v


def write_curated(rows: list[dict], report: CurationReport) -> dict[str, Path]:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    inst = PROCESSED / "instances.csv"
    mols = PROCESSED / "molecules.csv"
    excl = PROCESSED / "exclusions.csv"

    _write_csv(inst, rows, [
        "isoform", "mol_id", "smiles", "n_heavy_atoms",
        "som_primary", "som_all", "n_sources",
    ])
    _write_csv(mols, molecule_table(rows), [
        "smiles", "mol_id", "scaffold", "n_heavy_atoms", "n_rings", "mw",
        "isoforms", "n_isoforms",
    ])
    _write_csv(
        excl,
        [vars(e) for e in report.exclusions],
        ["isoform", "record", "mol_id", "reason", "detail"],
    )
    return {"instances": inst, "molecules": mols, "exclusions": excl}


def load_instances(path: Path | None = None) -> list[dict]:
    """Read instances.csv back, parsing the label columns into lists of ints."""
    path = path or PROCESSED / "instances.csv"
    out = []
    with path.open(newline="") as fh:
        for r in csv.DictReader(fh):
            out.append({
                "isoform": r["isoform"],
                "mol_id": r["mol_id"],
                "smiles": r["smiles"],
                "n_heavy_atoms": int(r["n_heavy_atoms"]),
                "som_primary": [int(x) for x in r["som_primary"].split()],
                "som_all": [int(x) for x in r["som_all"].split()],
                "n_sources": int(r["n_sources"]),
            })
    return out


def load_molecules(path: Path | None = None) -> dict[str, dict]:
    path = path or PROCESSED / "molecules.csv"
    out = {}
    with path.open(newline="") as fh:
        for r in csv.DictReader(fh):
            out[r["smiles"]] = {
                "smiles": r["smiles"],
                "mol_id": r["mol_id"],
                "scaffold": r["scaffold"],
                "n_heavy_atoms": int(r["n_heavy_atoms"]),
                "n_rings": int(r["n_rings"]),
                "mw": float(r["mw"]),
                "isoforms": r["isoforms"].split(),
                "n_isoforms": int(r["n_isoforms"]),
            }
    return out
