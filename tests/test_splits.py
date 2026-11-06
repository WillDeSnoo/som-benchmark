"""The split guarantees. If any of these break, every number is wrong."""

import numpy as np
import pytest

from som.data import load_instances, load_molecules
from som.paths import PROCESSED
from som.splits import SPLIT_MODES, make_split


@pytest.fixture(scope="module")
def data():
    if not (PROCESSED / "instances.csv").exists():
        pytest.skip("run `make curate` first")
    return load_instances(), load_molecules()


@pytest.mark.parametrize("mode", SPLIT_MODES)
@pytest.mark.parametrize("seed", [0, 1, 7])
def test_molecule_partitions_are_disjoint(data, mode, seed):
    instances, molecules = data
    sp = make_split(instances, molecules, mode, seed)
    tr, va, te = set(sp.train_mols), set(sp.val_mols), set(sp.test_mols)
    assert not tr & te
    assert not tr & va
    assert not va & te
    assert len(tr | va | te) == len({r["smiles"] for r in instances})


@pytest.mark.parametrize("mode", SPLIT_MODES)
@pytest.mark.parametrize("seed", [0, 1, 7])
def test_instances_of_one_molecule_never_straddle_the_split(data, mode, seed):
    """The leak this project most needs to avoid.

    One molecule appears in up to nine isoform instances. If any of them land
    on different sides, the test set contains structures the model trained on.
    """
    instances, molecules = data
    sp = make_split(instances, molecules, mode, seed)
    side = {}
    for name, ids in (("train", sp.train_inst), ("val", sp.val_inst),
                      ("test", sp.test_inst)):
        for i in ids:
            smi = instances[int(i)]["smiles"]
            if smi in side:
                assert side[smi] == name, (
                    f"{smi} appears in both {side[smi]} and {name}"
                )
            side[smi] = name
    assert len(side) == len({r["smiles"] for r in instances})


@pytest.mark.parametrize("seed", [0, 1, 7])
def test_scaffold_split_shares_no_scaffold(data, seed):
    instances, molecules = data
    sp = make_split(instances, molecules, "scaffold", seed)
    tr = {molecules[s]["scaffold"] for s in sp.train_mols}
    te = {molecules[s]["scaffold"] for s in sp.test_mols}
    va = {molecules[s]["scaffold"] for s in sp.val_mols}
    tr.discard("")
    te.discard("")
    va.discard("")
    assert not tr & te, "a scaffold appears in both train and test"
    assert not va & te, "a scaffold appears in both validation and test"


@pytest.mark.parametrize("mode", SPLIT_MODES)
def test_every_instance_is_assigned_exactly_once(data, mode):
    instances, molecules = data
    sp = make_split(instances, molecules, mode, 0)
    allids = np.concatenate([sp.train_inst, sp.val_inst, sp.test_inst])
    assert allids.size == len(instances)
    assert np.unique(allids).size == len(instances)


@pytest.mark.parametrize("mode", SPLIT_MODES)
def test_splits_are_deterministic_given_the_seed(data, mode):
    instances, molecules = data
    a = make_split(instances, molecules, mode, 3)
    b = make_split(instances, molecules, mode, 3)
    assert a.test_mols == b.test_mols
    assert a.train_mols == b.train_mols


@pytest.mark.parametrize("mode", SPLIT_MODES)
def test_different_seeds_give_different_splits(data, mode):
    instances, molecules = data
    a = make_split(instances, molecules, mode, 0)
    b = make_split(instances, molecules, mode, 1)
    assert a.test_mols != b.test_mols


@pytest.mark.parametrize("mode", SPLIT_MODES)
def test_test_fraction_is_roughly_respected(data, mode):
    instances, molecules = data
    n_mols = len({r["smiles"] for r in instances})
    sp = make_split(instances, molecules, mode, 0, test_frac=0.2)
    frac = len(sp.test_mols) / n_mols
    # whole scaffold groups are indivisible, so the realised fraction overshoots
    assert 0.15 < frac < 0.30, frac


def test_scaffold_split_is_harder_than_random(data):
    """Sanity check that the scaffold split actually separates chemistry."""
    from som.splits import split_diagnostics

    instances, molecules = data
    rnd = [
        split_diagnostics(make_split(instances, molecules, "random", s), molecules)
        for s in range(3)
    ]
    scf = [
        split_diagnostics(make_split(instances, molecules, "scaffold", s), molecules)
        for s in range(3)
    ]
    rnd_sim = np.mean([d["nn_tanimoto_median"] for d in rnd])
    scf_sim = np.mean([d["nn_tanimoto_median"] for d in scf])
    assert scf_sim < rnd_sim, (scf_sim, rnd_sim)
    assert all(d["shared_scaffolds"] == 0 for d in scf)
