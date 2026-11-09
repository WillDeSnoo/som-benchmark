"""Atom featurization for the three modelling approaches.

Three deliberately different representations, so that the comparison says
something:

``fp``
    Hashed circular atom environments (Morgan identifiers centred on the atom,
    radii 0-2) plus the isoform. Pure topology, no chemistry put in by hand.
    This is the "just throw a fingerprint at it" baseline.

``phys``
    A small, interpretable set of cheap physics-flavoured descriptors: simple
    Hückel frontier densities and free valence, Gasteiger charges, E-state
    indices, Crippen atomic contributions, and basic valence/ring context. No
    hashed environment information at all, so it cannot smuggle in the
    baseline's representation. 55 features instead of a thousand.

``graph``
    Basic atom and bond features for the GNN to build its own representation
    from. Deliberately excludes the Hückel descriptors, so that "learned
    representation" and "hand-built physics" are compared fairly.

Features are computed once per unique molecule and reused across the isoforms
that molecule appears in.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.EState import EStateIndices

from . import ISOFORMS
from .huckel import HUCKEL_FIELDS, as_vector, huckel_atoms

RDLogger.DisableLog("rdApp.*")

FP_BITS = 1024
FP_RADIUS = 2

_ELEMENTS = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]
_HYBRIDS = [
    Chem.HybridizationType.SP,
    Chem.HybridizationType.SP2,
    Chem.HybridizationType.SP3,
]


# --------------------------------------------------------------------------
# containers
# --------------------------------------------------------------------------

@dataclass
class AtomTable:
    """Flat atom-level design matrix with the bookkeeping to regroup it.

    `inst_id` indexes into the instance list, so all atoms of one
    (molecule, isoform) instance share a value. That is what the top-k metrics
    group on.
    """
    X: np.ndarray
    y: np.ndarray
    inst_id: np.ndarray
    atom_idx: np.ndarray
    feature_names: list[str]

    def __len__(self) -> int:
        return int(self.X.shape[0])


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _onehot(value, options) -> list[float]:
    vec = [0.0] * (len(options) + 1)
    try:
        vec[options.index(value)] = 1.0
    except ValueError:
        vec[-1] = 1.0
    return vec


def _safe(x: float) -> float:
    return float(x) if np.isfinite(x) else 0.0


def isoform_onehot(isoform: str) -> list[float]:
    return [1.0 if isoform == i else 0.0 for i in ISOFORMS]


def _ring_info(mol: Chem.Mol):
    ri = mol.GetRingInfo()
    min_ring, n_rings = [], []
    for a in mol.GetAtoms():
        sizes = [len(r) for r in ri.AtomRings() if a.GetIdx() in r]
        min_ring.append(min(sizes) if sizes else 0)
        n_rings.append(len(sizes))
    return min_ring, n_rings


# --------------------------------------------------------------------------
# per-molecule feature blocks
# --------------------------------------------------------------------------

def basic_atom_block(mol: Chem.Mol) -> tuple[np.ndarray, list[str]]:
    """Valence, element, ring and hybridization context. Shared by phys/graph."""
    min_ring, n_rings = _ring_info(mol)
    rows = []
    for a in mol.GetAtoms():
        i = a.GetIdx()
        nbrs = list(a.GetNeighbors())
        aromatic_nbrs = sum(1 for n in nbrs if n.GetIsAromatic())
        hetero_nbrs = sum(1 for n in nbrs if n.GetSymbol() not in ("C", "H"))
        has_double = any(
            b.GetBondType() == Chem.BondType.DOUBLE for b in a.GetBonds()
        )
        has_triple = any(
            b.GetBondType() == Chem.BondType.TRIPLE for b in a.GetBonds()
        )
        benzylic = (
            a.GetSymbol() == "C"
            and not a.GetIsAromatic()
            and a.GetTotalNumHs() > 0
            and aromatic_nbrs > 0
        )
        rows.append(
            _onehot(a.GetSymbol(), _ELEMENTS)
            + _onehot(a.GetHybridization(), _HYBRIDS)
            + [
                float(a.GetTotalNumHs()),
                float(a.GetDegree()),
                float(a.GetFormalCharge()),
                float(a.GetIsAromatic()),
                float(a.IsInRing()),
                float(min_ring[i]),
                float(n_rings[i]),
                float(aromatic_nbrs),
                float(hetero_nbrs),
                float(has_double),
                float(has_triple),
                float(benzylic),
                float(a.GetDegree() == 1),
            ]
        )
    names = (
        [f"el_{e}" for e in _ELEMENTS] + ["el_other"]
        + ["hyb_sp", "hyb_sp2", "hyb_sp3", "hyb_other"]
        + [
            "n_hydrogens", "heavy_degree", "formal_charge", "is_aromatic",
            "in_ring", "min_ring_size", "n_rings_at_atom", "aromatic_neighbours",
            "hetero_neighbours", "has_double_bond", "has_triple_bond",
            "is_benzylic_or_allylic_CH", "is_terminal",
        ]
    )
    return np.asarray(rows, dtype=np.float64), names


def charge_block(mol: Chem.Mol) -> tuple[np.ndarray, list[str]]:
    """Gasteiger charges, E-state indices and Crippen atomic contributions."""
    work = Chem.Mol(mol)
    try:
        AllChem.ComputeGasteigerCharges(work)
        q = [
            _safe(float(a.GetProp("_GasteigerCharge")))
            for a in work.GetAtoms()
        ]
    except Exception:  # noqa: BLE001
        q = [0.0] * mol.GetNumAtoms()
    try:
        estate = [_safe(v) for v in EStateIndices(mol)]
    except Exception:  # noqa: BLE001
        estate = [0.0] * mol.GetNumAtoms()
    try:
        crippen = rdMolDescriptors._CalcCrippenContribs(mol)
    except Exception:  # noqa: BLE001
        crippen = [(0.0, 0.0)] * mol.GetNumAtoms()

    rows = []
    for a in mol.GetAtoms():
        i = a.GetIdx()
        nq = [q[n.GetIdx()] for n in a.GetNeighbors()] or [0.0]
        logp, mr = crippen[i]
        rows.append([
            q[i], estate[i], _safe(logp), _safe(mr),
            max(nq), min(nq), float(np.mean(nq)),
        ])
    names = [
        "gasteiger_charge", "estate_index", "crippen_logp", "crippen_mr",
        "neighbour_charge_max", "neighbour_charge_min", "neighbour_charge_mean",
    ]
    return np.asarray(rows, dtype=np.float64), names


def huckel_block(mol: Chem.Mol) -> tuple[np.ndarray, list[str]]:
    """Hückel descriptors, plus their maximum over heavy neighbours.

    The neighbour aggregation is the point of including it: an sp3 benzylic
    carbon has no pi system of its own, but it sits next to one, and that is
    exactly the situation where hydrogen abstraction is easy. Propagating the
    neighbour's frontier density gives the model a handle on that.
    """
    recs = huckel_atoms(mol)
    base = np.asarray([as_vector(r) for r in recs], dtype=np.float64)

    prop_fields = ["f_electrophilic", "f_radical", "free_valence", "q_pi"]
    cols = [HUCKEL_FIELDS.index(f) for f in prop_fields]
    extra = np.zeros((mol.GetNumAtoms(), 2 * len(cols) + 1))
    for a in mol.GetAtoms():
        i = a.GetIdx()
        nb = [n.GetIdx() for n in a.GetNeighbors()]
        if nb:
            sub = base[np.asarray(nb)][:, cols]
            extra[i, : len(cols)] = sub.max(axis=0)
            extra[i, len(cols) : 2 * len(cols)] = sub.mean(axis=0)
            extra[i, -1] = float(sum(1 for j in nb if recs[j].in_pi))

    names = (
        [f"huckel_{f}" for f in HUCKEL_FIELDS]
        + [f"nbr_max_huckel_{f}" for f in prop_fields]
        + [f"nbr_mean_huckel_{f}" for f in prop_fields]
        + ["n_pi_neighbours"]
    )
    return np.hstack([base, extra]), names


def fp_block(mol: Chem.Mol, n_bits: int = FP_BITS, radius: int = FP_RADIUS):
    """Hashed Morgan atom environments centred on each atom, radii 0..radius."""
    n = mol.GetNumAtoms()
    out = np.zeros((n, n_bits), dtype=np.float64)
    for i in range(n):
        fp = AllChem.GetMorganFingerprint(mol, radius, fromAtoms=[i])
        for ident, count in fp.GetNonzeroElements().items():
            out[i, ident % n_bits] += count
    names = [f"fp_{b}" for b in range(n_bits)]
    return out, names


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------

#: cached BDE table, loaded on first use by the `*_bde` feature kinds
_BDE_TABLE: dict | None = None


def _bde_table() -> dict:
    global _BDE_TABLE
    if _BDE_TABLE is None:
        from .bde import bde_path, load_bde

        _BDE_TABLE = load_bde()
        if not _BDE_TABLE:
            raise FileNotFoundError(
                f"{bde_path()} not found. Run `make bde` first "
                "(it needs the separate alfabet-env; see README)."
            )
    return _BDE_TABLE


def molecule_features(
    mol: Chem.Mol, kind: str, smiles: str | None = None
) -> tuple[np.ndarray, list[str]]:
    """Feature block for every atom of one molecule, excluding the isoform.

    The `*_bde` kinds append the ALFABET bond-dissociation-energy block to the
    corresponding base kind, so that `phys_bde` is exactly `phys` plus BDE and
    `graph_bde` is exactly `graph` plus BDE. Keeping the extension strictly
    additive is what makes the paired ablations interpretable.
    """
    if kind == "fp":
        return fp_block(mol)

    base = {
        "phys": [basic_atom_block, charge_block, huckel_block],
        "graph": [basic_atom_block, charge_block],
        "phys_bde": [basic_atom_block, charge_block, huckel_block],
        "graph_bde": [basic_atom_block, charge_block],
    }.get(kind)
    if base is None:
        raise ValueError(f"unknown feature kind: {kind!r}")

    blocks = [fn(mol) for fn in base]
    if kind.endswith("_bde"):
        from .bde import bde_block

        smi = smiles if smiles is not None else Chem.MolToSmiles(mol)
        blocks.append(bde_block(mol, smi, _bde_table()))

    X = np.hstack([b[0] for b in blocks])
    names = [n for b in blocks for n in b[1]]
    return X, names


def build_atom_table(
    instances: list[dict],
    kind: str,
    label_field: str = "som_all",
    with_isoform: bool = True,
) -> AtomTable:
    """Featurize every atom of every instance.

    Molecules are featurized once and reused across isoforms, which matters
    because 679 molecules expand to 2003 instances.
    """
    cache: dict[str, np.ndarray] = {}
    names: list[str] | None = None

    Xs, ys, inst_ids, atom_idxs = [], [], [], []
    for inst_id, row in enumerate(instances):
        smi = row["smiles"]
        if smi not in cache:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                raise ValueError(f"curated SMILES failed to parse: {smi}")
            block, blk_names = molecule_features(mol, kind, smiles=smi)
            cache[smi] = block
            if names is None:
                names = list(blk_names)
        block = cache[smi]
        n = block.shape[0]
        if n != row["n_heavy_atoms"]:
            raise ValueError(
                f"atom count mismatch for {smi}: {n} vs {row['n_heavy_atoms']}"
            )

        if with_isoform:
            iso = np.tile(isoform_onehot(row["isoform"]), (n, 1))
            Xs.append(np.hstack([block, iso]))
        else:
            Xs.append(block)

        labels = np.zeros(n, dtype=np.int8)
        labels[np.asarray(row[label_field], dtype=int)] = 1
        ys.append(labels)
        inst_ids.append(np.full(n, inst_id, dtype=np.int32))
        atom_idxs.append(np.arange(n, dtype=np.int32))

    assert names is not None
    if with_isoform:
        names = names + [f"isoform_{i}" for i in ISOFORMS]

    return AtomTable(
        X=np.vstack(Xs),
        y=np.concatenate(ys),
        inst_id=np.concatenate(inst_ids),
        atom_idx=np.concatenate(atom_idxs),
        feature_names=names,
    )
