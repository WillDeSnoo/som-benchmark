"""Bond dissociation energy features, precomputed by ALFABET.

`scripts/07_bde.py` writes `data/interim/alfabet_bde.csv` from a separate
environment; this module turns it into per-atom features.

The absolute BDE is the least interesting of these. Hydrogen abstraction
competes *within* a molecule, and the quantity a model needs is "how much
weaker is this bond than the weakest one in the same molecule" -- so the
within-molecule delta and rank are the features that matter for a top-k metric.
Absolute BDE largely just restates atom type (benzylic ~90, aromatic ~112),
which the model already knows.

Atoms with no abstractable hydrogen get `NO_H_BDE`, a value above any real C-H
bond, plus an explicit indicator. Filling with zero would make a tree read them
as the weakest bonds in the molecule, which is the opposite of the truth.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .paths import INTERIM

#: stand-in for "no abstractable hydrogen"; must match scripts/07_bde.py
NO_H_BDE = 140.0

BDE_FIELDS = [
    "bde_has",
    "bde",
    "bde_minus_mol_min",
    "bde_rank_in_mol",
    "bde_is_mol_min",
    "bdfe_minus_mol_min",
    "bde_nbr_min",
    "bde_outside_domain",
]


def bde_path() -> Path:
    return INTERIM / "alfabet_bde.csv"


def load_bde(path: Path | None = None) -> dict[str, dict[int, tuple[float, float]]]:
    """{smiles: {atom_idx: (bde, bdfe)}}, plus a domain flag under key -1."""
    path = path or bde_path()
    out: dict[str, dict[int, tuple[float, float]]] = {}
    if not path.exists():
        return out
    with path.open(newline="") as fh:
        for r in csv.DictReader(fh):
            d = out.setdefault(r["smiles"], {})
            d[int(r["atom_idx"])] = (float(r["bde"]), float(r["bdfe"]))
            if not int(r["is_valid"]):
                d[-1] = (1.0, 1.0)
    return out


def bde_block(mol, smiles: str, table: dict) -> tuple[np.ndarray, list[str]]:
    """Per-atom BDE features for one molecule."""
    n = mol.GetNumAtoms()
    rec = table.get(smiles, {})
    outside = 1.0 if -1 in rec else 0.0

    bde = np.full(n, NO_H_BDE, dtype=np.float64)
    bdfe = np.full(n, NO_H_BDE, dtype=np.float64)
    has = np.zeros(n, dtype=np.float64)
    for i, (b, bf) in rec.items():
        if i < 0 or i >= n:
            continue
        bde[i] = b
        bdfe[i] = bf
        has[i] = 1.0

    real = bde[has > 0]
    if real.size:
        mol_min = float(real.min())
        fmin = float(bdfe[has > 0].min())
        # rank among atoms that actually have an abstractable hydrogen
        order = np.argsort(np.argsort(bde + (1 - has) * 1e6))
        rank = order / max(1, n - 1)
    else:
        mol_min = NO_H_BDE
        fmin = NO_H_BDE
        rank = np.zeros(n)

    delta = bde - mol_min
    fdelta = bdfe - fmin
    is_min = ((has > 0) & (np.abs(bde - mol_min) < 1e-9)).astype(np.float64)

    nbr_min = np.full(n, NO_H_BDE, dtype=np.float64)
    for a in mol.GetAtoms():
        vals = [bde[x.GetIdx()] for x in a.GetNeighbors()]
        if vals:
            nbr_min[a.GetIdx()] = min(vals)

    X = np.column_stack([
        has, bde, delta, rank, is_min, fdelta, nbr_min,
        np.full(n, outside),
    ])
    return X, list(BDE_FIELDS)
