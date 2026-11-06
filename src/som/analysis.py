"""Failure analysis: where do these models actually go wrong.

Three cuts, matching the three ways this task can mislead you:

**By CYP isoform.** The nine isoforms are not equally represented (3A4 has 475
instances, 2A6 has 105) and they have different chemistry. A single pooled
number hides the fact that the rare isoforms may be carried by the common ones.

**By site type.** A positive atom is classified by the transformation it
implies -- aromatic hydroxylation, N-dealkylation, S-oxidation, and so on --
from its own element and its immediate environment. Per-site-type recall says
which chemistry the model has actually learned, as opposed to which it is
getting for free because the dataset is mostly aromatic and benzylic carbons.

**By how far the test molecule sits from training.** The nearest-neighbour
Tanimoto similarity of each test molecule to the training set turns the
random-versus-scaffold comparison from two numbers into a curve, which is the
more useful form: it says how performance decays with novelty rather than just
that it does.
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

_HALOGENS = {"F", "Cl", "Br", "I"}

SITE_TYPES = [
    "S-oxidation",
    "P-oxidation / desulfuration",
    "N-oxidation",
    "O-site (direct)",
    "oxidative dehalogenation",
    "aromatic C-H hydroxylation",
    "aromatic C (substituted)",
    "N-dealkylation (alpha to N)",
    "O-dealkylation (alpha to O)",
    "benzylic / allylic hydroxylation",
    "methyl hydroxylation",
    "aliphatic C-H hydroxylation",
    "other",
]


def classify_site(mol: Chem.Mol, idx: int) -> str:
    """Name the transformation implied by metabolism at this atom.

    Order matters: a carbon alpha to nitrogen is called N-dealkylation even if
    it is also a methyl, because that is the transformation a chemist would
    name, and lumping it in with methyl hydroxylation would hide a mechanism.
    """
    a = mol.GetAtomWithIdx(idx)
    sym = a.GetSymbol()

    if sym == "S":
        return "S-oxidation"
    if sym == "P":
        return "P-oxidation / desulfuration"
    if sym == "N":
        return "N-oxidation"
    if sym == "O":
        return "O-site (direct)"
    if sym in _HALOGENS:
        return "oxidative dehalogenation"
    if sym != "C":
        return "other"

    if a.GetIsAromatic():
        return (
            "aromatic C-H hydroxylation" if a.GetTotalNumHs() > 0
            else "aromatic C (substituted)"
        )

    nbrs = list(a.GetNeighbors())
    if any(n.GetSymbol() == "N" for n in nbrs):
        return "N-dealkylation (alpha to N)"
    if any(n.GetSymbol() == "O" for n in nbrs):
        return "O-dealkylation (alpha to O)"
    if any(n.GetIsAromatic() for n in nbrs):
        return "benzylic / allylic hydroxylation"
    if any(
        b.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE)
        for n in nbrs for b in n.GetBonds()
    ):
        return "benzylic / allylic hydroxylation"
    if a.GetTotalNumHs() == 3:
        return "methyl hydroxylation"
    return "aliphatic C-H hydroxylation"


def site_types_for(smiles: str, positives: list[int]) -> dict[int, str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    return {i: classify_site(mol, i) for i in positives}


# --------------------------------------------------------------------------
# per-atom rank within its instance
# --------------------------------------------------------------------------

def atom_ranks(
    p: np.ndarray, inst: np.ndarray, seed: int = 0
) -> np.ndarray:
    """Rank of each atom within its own instance, 0 = highest scoring.

    Ties are broken by a seeded permutation, matching `metrics.topk_hits`, so
    that the per-atom and per-instance views of the same predictions agree.
    """
    rng = np.random.default_rng(seed)
    ranks = np.empty(p.size, dtype=np.int32)
    order = np.argsort(inst, kind="stable")
    sorted_inst = inst[order]
    bounds = np.flatnonzero(np.diff(sorted_inst)) + 1
    for chunk in np.split(order, bounds):
        if chunk.size == 0:
            continue
        jitter = rng.permutation(chunk.size)
        local = np.lexsort((jitter, -p[chunk]))
        ranks[chunk[local]] = np.arange(chunk.size, dtype=np.int32)
    return ranks


def recall_at_k_by_group(
    y: np.ndarray,
    ranks: np.ndarray,
    groups: np.ndarray,
    k: int = 2,
) -> dict[str, dict]:
    """Fraction of positive atoms in each group that land in the top k."""
    out: dict[str, dict] = {}
    pos = y == 1
    for g in np.unique(groups[pos]):
        sel = pos & (groups == g)
        n = int(sel.sum())
        hit = int((ranks[sel] < k).sum())
        out[str(g)] = {
            "n_positive_atoms": n,
            f"recall_at_{k}": hit / n if n else float("nan"),
        }
    return out


def instance_topk_table(
    y: np.ndarray, p: np.ndarray, inst: np.ndarray, ks=(1, 2), seed: int = 0
) -> dict[int, dict]:
    """Per-instance top-k outcome, keyed by instance id."""
    ranks = atom_ranks(p, inst, seed=seed)
    out: dict[int, dict] = {}
    for iid in np.unique(inst):
        sel = inst == iid
        yy, rr = y[sel], ranks[sel]
        if yy.sum() == 0:
            continue
        best = int(rr[yy == 1].min())
        rec = {"best_rank": best, "n_atoms": int(sel.sum()),
               "n_positive": int(yy.sum())}
        for k in ks:
            rec[f"top{k}"] = bool(best < k)
        out[int(iid)] = rec
    return out


def bucket(value: float, edges: list[float], labels: list[str]) -> str:
    for e, lab in zip(edges, labels):
        if value <= e:
            return lab
    return labels[-1]
