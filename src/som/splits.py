"""Train / validation / test splits.

Two things here are easy to get wrong and are the main reason SOM numbers in
the literature are hard to compare:

**Grouping by molecule, not by instance.** A molecule appears in up to nine
isoform-specific instances. Splitting on instances would put the same structure
in both train and test with a different isoform label attached, which leaks
badly: the atom environments are identical and the labels are strongly
correlated across isoforms. Every split here partitions *molecules*.

**Random scaffold splits, not one fixed scaffold split.** The common
"Bemis-Murcko split" assigns the largest scaffold groups to train
deterministically, which gives exactly one arbitrary partition and no way to
put an error bar on it. Here the scaffold *groups* are shuffled per seed and
drawn into the test set until it reaches the target size, so the scaffold
constraint is respected and the split still varies across seeds.

The test set is carved out first and the validation set comes out of the
remaining training molecules, so nothing used for tuning or early stopping ever
touches test.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from .features import AtomTable

RDLogger.DisableLog("rdApp.*")

SPLIT_MODES = ("random", "scaffold")


@dataclass
class Split:
    mode: str
    seed: int
    train_inst: np.ndarray
    val_inst: np.ndarray
    test_inst: np.ndarray
    train_mols: list[str]
    val_mols: list[str]
    test_mols: list[str]

    def summary(self) -> dict:
        return {
            "mode": self.mode,
            "seed": self.seed,
            "train_instances": int(self.train_inst.size),
            "val_instances": int(self.val_inst.size),
            "test_instances": int(self.test_inst.size),
            "train_molecules": len(self.train_mols),
            "val_molecules": len(self.val_mols),
            "test_molecules": len(self.test_mols),
        }


def group_key(smiles: str, molecules: dict[str, dict], mode: str) -> str:
    """The unit that must not straddle the split boundary."""
    if mode == "random":
        return smiles
    if mode == "scaffold":
        scaf = molecules[smiles]["scaffold"]
        # Acyclic molecules have no Bemis-Murcko scaffold. Lumping all of them
        # into one group would force them all onto the same side of every
        # split; each is instead treated as its own group.
        return scaf if scaf else f"__acyclic__{smiles}"
    raise ValueError(f"unknown split mode: {mode!r}")


def _draw_groups(
    groups: dict[str, list[str]], target: int, rng: np.random.Generator
) -> tuple[list[str], list[str]]:
    """Shuffle whole groups into a holdout until it reaches `target` molecules."""
    keys = sorted(groups)
    rng.shuffle(keys)
    held: list[str] = []
    rest: list[str] = []
    for k in keys:
        if len(held) < target:
            held.extend(groups[k])
        else:
            rest.extend(groups[k])
    return held, rest


def make_split(
    instances: list[dict],
    molecules: dict[str, dict],
    mode: str,
    seed: int,
    test_frac: float = 0.2,
    val_frac: float = 0.1,
) -> Split:
    """Partition molecules, then map the partition onto instance indices."""
    rng = np.random.default_rng(seed)
    all_mols = sorted({r["smiles"] for r in instances})

    groups: dict[str, list[str]] = {}
    for smi in all_mols:
        groups.setdefault(group_key(smi, molecules, mode), []).append(smi)

    n_total = len(all_mols)
    test_mols, remaining = _draw_groups(
        groups, int(round(test_frac * n_total)), rng
    )

    # validation comes out of the remaining training molecules, under the same
    # grouping rule, so tuning never sees a test scaffold
    rem_groups: dict[str, list[str]] = {}
    for smi in remaining:
        rem_groups.setdefault(group_key(smi, molecules, mode), []).append(smi)
    val_mols, train_mols = _draw_groups(
        rem_groups, int(round(val_frac * n_total)), rng
    )

    test_set, val_set = set(test_mols), set(val_mols)
    train_i, val_i, test_i = [], [], []
    for i, row in enumerate(instances):
        smi = row["smiles"]
        if smi in test_set:
            test_i.append(i)
        elif smi in val_set:
            val_i.append(i)
        else:
            train_i.append(i)

    return Split(
        mode=mode,
        seed=seed,
        train_inst=np.asarray(train_i, dtype=np.int64),
        val_inst=np.asarray(val_i, dtype=np.int64),
        test_inst=np.asarray(test_i, dtype=np.int64),
        train_mols=sorted(train_mols),
        val_mols=sorted(val_mols),
        test_mols=sorted(test_mols),
    )


def atom_mask(table: AtomTable, inst_ids: np.ndarray) -> np.ndarray:
    """Boolean mask over atom rows belonging to the given instances."""
    return np.isin(table.inst_id, inst_ids)


# --------------------------------------------------------------------------
# split diagnostics
# --------------------------------------------------------------------------

def _morgan(mol: Chem.Mol):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return gen.GetFingerprint(mol)


def nearest_neighbour_similarity(
    query_smiles: list[str], reference_smiles: list[str]
) -> np.ndarray:
    """For each query molecule, its maximum Tanimoto similarity to the reference.

    This is what quantifies *how much harder* the scaffold split actually is:
    if the nearest training neighbour of a test molecule is still at Tanimoto
    0.8, the split has not really separated anything.
    """
    ref = []
    for smi in reference_smiles:
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            ref.append(_morgan(m))
    out = []
    for smi in query_smiles:
        m = Chem.MolFromSmiles(smi)
        if m is None or not ref:
            out.append(float("nan"))
            continue
        sims = DataStructs.BulkTanimotoSimilarity(_morgan(m), ref)
        out.append(max(sims) if sims else float("nan"))
    return np.asarray(out, dtype=np.float64)


def split_diagnostics(split: Split, molecules: dict[str, dict]) -> dict:
    """Scaffold overlap and train-test similarity for one split."""
    tr_scaf = {molecules[s]["scaffold"] for s in split.train_mols}
    te_scaf = {molecules[s]["scaffold"] for s in split.test_mols}
    tr_scaf.discard("")
    te_scaf.discard("")
    shared = tr_scaf & te_scaf

    sims = nearest_neighbour_similarity(split.test_mols, split.train_mols)
    finite = sims[np.isfinite(sims)]
    return {
        "mode": split.mode,
        "seed": split.seed,
        "test_molecules": len(split.test_mols),
        "train_scaffolds": len(tr_scaf),
        "test_scaffolds": len(te_scaf),
        "shared_scaffolds": len(shared),
        "nn_tanimoto_mean": float(finite.mean()) if finite.size else float("nan"),
        "nn_tanimoto_median": float(np.median(finite)) if finite.size else float("nan"),
        "nn_tanimoto_p90": float(np.percentile(finite, 90)) if finite.size else float("nan"),
        "frac_test_with_nn_above_0.6": (
            float((finite > 0.6).mean()) if finite.size else float("nan")
        ),
    }
