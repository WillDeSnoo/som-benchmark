"""BDE feature block.

The chemistry checks here double as a validation of the whole ALFABET
pipeline: if the bond-index to heavy-atom mapping in `scripts/07_bde.py` were
wrong, toluene's methyl would not come out ~23 kcal/mol below its ring carbons.
"""

import numpy as np
import pytest
from rdkit import Chem

from som.bde import BDE_FIELDS, NO_H_BDE, bde_block, bde_path, load_bde
from som.features import molecule_features


@pytest.fixture(scope="module")
def table():
    if not bde_path().exists():
        pytest.skip("run `make bde` first")
    t = load_bde()
    if not t:
        pytest.skip("empty BDE table")
    return t


def test_block_shape_and_fields():
    mol = Chem.MolFromSmiles("CCO")
    X, names = bde_block(mol, "CCO", {})
    assert X.shape == (3, len(BDE_FIELDS))
    assert names == BDE_FIELDS
    assert np.all(np.isfinite(X))


def test_missing_molecule_falls_back_to_no_hydrogen_value():
    mol = Chem.MolFromSmiles("CCO")
    X, names = bde_block(mol, "not-in-table", {})
    bde = X[:, names.index("bde")]
    has = X[:, names.index("bde_has")]
    assert np.allclose(bde, NO_H_BDE)
    assert np.allclose(has, 0.0)


def test_no_hydrogen_atoms_get_the_high_sentinel_not_zero():
    """A zero fill would make a tree read these as the weakest bonds."""
    mol = Chem.MolFromSmiles("Cc1ccccc1")
    fake = {"Cc1ccccc1": {0: (89.6, 82.0)}}  # only the methyl has a BDE
    X, names = bde_block(mol, "Cc1ccccc1", fake)
    bde = X[:, names.index("bde")]
    has = X[:, names.index("bde_has")]
    assert bde[0] == pytest.approx(89.6)
    assert has[0] == 1.0
    # every other atom is unmeasured and must sit high, not at zero
    assert np.all(bde[1:] == NO_H_BDE)
    assert np.all(has[1:] == 0.0)
    assert bde.min() == pytest.approx(89.6)


def test_delta_is_zero_at_the_weakest_bond():
    mol = Chem.MolFromSmiles("Cc1ccccc1")
    fake = {"x": {0: (89.6, 82.0), 2: (110.9, 103.0), 3: (110.7, 102.8)}}
    X, names = bde_block(mol, "x", fake)
    delta = X[:, names.index("bde_minus_mol_min")]
    is_min = X[:, names.index("bde_is_mol_min")]
    assert delta[0] == pytest.approx(0.0)
    assert is_min[0] == 1.0
    assert delta[2] == pytest.approx(21.3, abs=0.01)
    assert is_min[2] == 0.0


def test_phys_bde_is_phys_plus_bde(table):
    mol = Chem.MolFromSmiles("Cc1ccccc1")
    smi = Chem.MolToSmiles(mol)
    Xa, na = molecule_features(mol, "phys", smiles=smi)
    Xb, nb = molecule_features(mol, "phys_bde", smiles=smi)
    assert nb[: len(na)] == na
    assert np.allclose(Xb[:, : len(na)], Xa)
    assert nb[len(na):] == BDE_FIELDS


def test_graph_bde_is_graph_plus_bde(table):
    mol = Chem.MolFromSmiles("Cc1ccccc1")
    smi = Chem.MolToSmiles(mol)
    Xa, na = molecule_features(mol, "graph", smiles=smi)
    Xb, nb = molecule_features(mol, "graph_bde", smiles=smi)
    assert nb[: len(na)] == na
    assert np.allclose(Xb[:, : len(na)], Xa)


def test_toluene_benzylic_is_far_weaker_than_its_ring(table):
    """Literature: benzylic C-H ~89.7, aromatic C-H ~113 kcal/mol."""
    smi = Chem.MolToSmiles(Chem.MolFromSmiles("Cc1ccccc1"))
    if smi not in table:
        pytest.skip("toluene not in this dataset")
    mol = Chem.MolFromSmiles(smi)
    X, names = bde_block(mol, smi, table)
    bde = X[:, names.index("bde")]
    has = X[:, names.index("bde_has")]

    methyl = [a.GetIdx() for a in mol.GetAtoms()
              if a.GetSymbol() == "C" and not a.GetIsAromatic()]
    ring_ch = [a.GetIdx() for a in mol.GetAtoms()
               if a.GetIsAromatic() and a.GetTotalNumHs() == 1]
    assert len(methyl) == 1 and len(ring_ch) == 5

    assert bde[methyl[0]] == pytest.approx(89.7, abs=3.0)
    for i in ring_ch:
        assert has[i] == 1.0
        assert bde[i] == pytest.approx(112.0, abs=4.0)
    assert min(bde[i] for i in ring_ch) - bde[methyl[0]] > 15.0


def test_benzene_ring_carbons_are_equivalent(table):
    smi = Chem.MolToSmiles(Chem.MolFromSmiles("c1ccccc1"))
    if smi not in table:
        pytest.skip("benzene not in this dataset")
    mol = Chem.MolFromSmiles(smi)
    X, names = bde_block(mol, smi, table)
    bde = X[:, names.index("bde")]
    assert bde.max() - bde.min() < 0.5, "symmetry-equivalent C-H must agree"


def test_dataset_bde_values_are_physically_plausible(table):
    vals = np.array([b for rec in table.values()
                     for k, (b, _) in rec.items() if k >= 0])
    assert vals.size > 1000
    # no real X-H bond in a drug-like molecule sits outside this window
    assert vals.min() > 40.0, vals.min()
    assert vals.max() < NO_H_BDE, vals.max()
    assert 90.0 < np.median(vals) < 115.0, np.median(vals)
