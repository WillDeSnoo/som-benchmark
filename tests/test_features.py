"""Feature-set properties, including the one the GNN ablation depends on.

The `gnn` vs `gnn_huckel` comparison is only meaningful if the plain GNN
genuinely cannot see the Hückel descriptors and the ablated one sees exactly
the set `phys_gbm` uses. That is a property of the feature blocks, so it is
asserted here rather than assumed.
"""

import numpy as np
import pytest
from rdkit import Chem

from som import ISOFORMS
from som.features import build_atom_table, molecule_features

SMILES = "COc1ccccc1CN(C)CCSc1ccccc1F"


@pytest.fixture(scope="module")
def mol():
    m = Chem.MolFromSmiles(SMILES)
    assert m is not None
    return m


def test_graph_features_contain_no_huckel_terms(mol):
    _, names = molecule_features(mol, "graph")
    offenders = [n for n in names if "huckel" in n.lower() or "pi_neighbours" in n]
    assert offenders == [], offenders


def test_phys_features_do_contain_huckel_terms(mol):
    _, names = molecule_features(mol, "phys")
    assert any(n.startswith("huckel_") for n in names)
    assert any(n.startswith("nbr_max_huckel_") for n in names)


def test_phys_is_graph_plus_huckel(mol):
    """The ablation is additive: phys = graph features, then Hückel.

    If this ever stopped holding, `gnn` and `gnn_huckel` would differ by more
    than the physics and the comparison would no longer isolate it.
    """
    Xg, ng = molecule_features(mol, "graph")
    Xp, np_names = molecule_features(mol, "phys")
    assert np_names[: len(ng)] == ng
    assert np.allclose(Xp[:, : len(ng)], Xg)
    extra = np_names[len(ng):]
    assert extra, "phys must add something"
    assert all(
        n.startswith(("huckel_", "nbr_max_huckel_", "nbr_mean_huckel_"))
        or n == "n_pi_neighbours"
        for n in extra
    ), extra


def test_fp_features_are_disjoint_from_physics(mol):
    """The fingerprint baseline must not contain hand-built descriptors."""
    _, names = molecule_features(mol, "fp")
    assert all(n.startswith("fp_") for n in names)


def test_feature_blocks_are_finite(mol):
    for kind in ("fp", "phys", "graph"):
        X, names = molecule_features(mol, kind)
        assert X.shape[0] == mol.GetNumAtoms()
        assert X.shape[1] == len(names)
        assert np.all(np.isfinite(X)), kind


def test_unknown_feature_kind_raises(mol):
    with pytest.raises(ValueError):
        molecule_features(mol, "nonsense")


def test_build_atom_table_puts_labels_on_the_right_atoms():
    instances = [
        {"isoform": "3A4", "mol_id": "x", "smiles": "CCOc1ccccc1",
         "n_heavy_atoms": 9, "som_primary": [1], "som_all": [1, 3],
         "n_sources": 1},
        {"isoform": "2D6", "mol_id": "x", "smiles": "CCOc1ccccc1",
         "n_heavy_atoms": 9, "som_primary": [0], "som_all": [0],
         "n_sources": 1},
    ]
    table = build_atom_table(instances, "phys", label_field="som_all")
    assert len(table) == 18
    first = table.inst_id == 0
    assert table.y[first].tolist() == [0, 1, 0, 1, 0, 0, 0, 0, 0]
    second = table.inst_id == 1
    assert table.y[second].tolist() == [1, 0, 0, 0, 0, 0, 0, 0, 0]


def test_isoform_one_hot_is_appended_and_differs_between_instances():
    instances = [
        {"isoform": "1A2", "mol_id": "x", "smiles": "CCO", "n_heavy_atoms": 3,
         "som_primary": [0], "som_all": [0], "n_sources": 1},
        {"isoform": "3A4", "mol_id": "x", "smiles": "CCO", "n_heavy_atoms": 3,
         "som_primary": [0], "som_all": [0], "n_sources": 1},
    ]
    table = build_atom_table(instances, "phys")
    assert table.feature_names[-len(ISOFORMS):] == [f"isoform_{i}" for i in ISOFORMS]
    iso_block = table.X[:, -len(ISOFORMS):]
    a = iso_block[table.inst_id == 0][0]
    b = iso_block[table.inst_id == 1][0]
    assert a.sum() == 1 and b.sum() == 1
    assert not np.array_equal(a, b)
    # the chemistry half of the row must be identical for the same molecule
    chem = table.X[:, : -len(ISOFORMS)]
    assert np.allclose(chem[table.inst_id == 0], chem[table.inst_id == 1])


def test_atom_count_mismatch_is_caught():
    bad = [{"isoform": "3A4", "mol_id": "x", "smiles": "CCO",
            "n_heavy_atoms": 4, "som_primary": [0], "som_all": [0],
            "n_sources": 1}]
    with pytest.raises(ValueError, match="atom count mismatch"):
        build_atom_table(bad, "phys")


def test_the_huckel_ablation_model_uses_the_graph_feature_set(monkeypatch):
    """`basic_gbm` must be `phys_gbm` minus the Hückel block, nothing else.

    The phys_gbm - basic_gbm comparison is only an ablation of Hückel theory
    if basic_gbm is fitted on exactly the `phys` columns with huckel_block
    dropped. `graph` is that set (asserted above to be a strict prefix of
    `phys`), so this pins the dispatch to it by intercepting the fit rather
    than by matching source text.
    """
    from som import models

    seen = {}

    def spy(name, kind, *a, **kw):
        seen["kind"] = kind
        raise _Stop

    class _Stop(Exception):
        pass

    monkeypatch.setattr(models, "_fit_gbm", spy)
    try:
        models.fit_predict("basic_gbm", [], None, "som_all", 0)
    except Exception:
        pass
    assert seen.get("kind") == "graph", (
        "basic_gbm must be fitted on the 'graph' feature kind; if the "
        f"dispatch changed (saw {seen.get('kind')!r}), the Hückel ablation "
        "no longer isolates Hückel"
    )
