"""Curation correctness, above all that the SOM labels survive standardization.

The labels are atom *indices*. Standardization strips salts and canonical
SMILES output reorders atoms, so an off-by-one or a lost mapping would quietly
move every label onto a neighbouring atom and still produce a plausible-looking
model. These tests check the mapping against chemistry that is not in doubt.
"""

import pytest
from rdkit import Chem

from som.data import (
    canonical_order_map,
    build_standardizer,
    load_instances,
    load_molecules,
    murcko_scaffold,
    parse_som_field,
)
from som.paths import PROCESSED

_HALOGENS = {"F", "Cl", "Br", "I"}


# --------------------------------------------------------------------------
# unit level
# --------------------------------------------------------------------------

def test_parse_som_field_is_one_indexed():
    assert parse_som_field("7") == [6]
    assert parse_som_field("1 2 3") == [0, 1, 2]
    assert parse_som_field("4,9") == [3, 8]
    assert parse_som_field("") == []
    assert parse_som_field(None) == []
    with pytest.raises(ValueError):
        parse_som_field("not_an_atom")


def test_canonical_order_map_round_trips():
    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    smiles, mapping = canonical_order_map(mol)
    reloaded = Chem.MolFromSmiles(smiles)
    assert reloaded.GetNumAtoms() == mol.GetNumAtoms()
    for old, new in mapping.items():
        assert (
            mol.GetAtomWithIdx(old).GetSymbol()
            == reloaded.GetAtomWithIdx(new).GetSymbol()
        )


def test_standardizer_keeps_atom_tags_through_salt_stripping():
    """A labelled atom on the parent must keep its tag when the salt goes."""
    std = build_standardizer()
    mol = Chem.MolFromSmiles("CCN.Cl")
    for a in mol.GetAtoms():
        a.SetIntProp("_somOrigIdx", a.GetIdx())
    clean = std(mol)
    assert clean.GetNumAtoms() == 3  # chloride removed
    tags = sorted(a.GetIntProp("_somOrigIdx") for a in clean.GetAtoms())
    assert tags == [0, 1, 2]


def test_standardizer_clears_mol_properties():
    """Upstream Citation fields hold non-UTF-8 bytes; they must not survive."""
    std = build_standardizer()
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("Citation", "whatever")
    clean = std(mol)
    assert "Citation" not in list(clean.GetPropNames())


def test_murcko_scaffold_of_acyclic_is_empty():
    assert murcko_scaffold("CCCCO") == ""
    assert murcko_scaffold("c1ccccc1CC") == "c1ccccc1"


# --------------------------------------------------------------------------
# curated dataset level
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def curated():
    if not (PROCESSED / "instances.csv").exists():
        pytest.skip("run `make curate` first")
    return load_instances(), load_molecules()


def test_every_instance_has_at_least_one_label(curated):
    instances, _ = curated
    for row in instances:
        assert row["som_primary"], row["mol_id"]
        assert set(row["som_primary"]) <= set(row["som_all"])


def test_labels_are_in_range_and_smiles_reparses(curated):
    instances, _ = curated
    for row in instances:
        mol = Chem.MolFromSmiles(row["smiles"])
        assert mol is not None, row["smiles"]
        assert mol.GetNumAtoms() == row["n_heavy_atoms"]
        for i in row["som_all"]:
            assert 0 <= i < row["n_heavy_atoms"], (row["mol_id"], i)


def test_no_duplicate_molecule_isoform_pairs(curated):
    instances, _ = curated
    seen = set()
    for row in instances:
        key = (row["smiles"], row["isoform"])
        assert key not in seen, key
        seen.add(key)


def test_fluorine_is_essentially_never_a_site_of_metabolism(curated):
    """C-F is not cleaved by CYPs.

    This is the strongest available check that the index mapping is right: if
    the labels had drifted by one atom, fluorine would start appearing as a
    site at roughly its background frequency instead of almost never.
    """
    instances, _ = curated
    n_f_atoms = n_f_labelled = 0
    for row in instances:
        mol = Chem.MolFromSmiles(row["smiles"])
        f_idx = {a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "F"}
        n_f_atoms += len(f_idx)
        n_f_labelled += len(f_idx & set(row["som_all"]))
    assert n_f_atoms > 200, "expected plenty of fluorines in a drug-like set"
    assert n_f_labelled / n_f_atoms < 0.02, (n_f_labelled, n_f_atoms)


def test_labelled_carbons_overwhelmingly_bear_hydrogens(curated):
    """CYP carbon oxidation proceeds by hydrogen abstraction."""
    instances, _ = curated
    with_h = total = 0
    for row in instances:
        mol = Chem.MolFromSmiles(row["smiles"])
        for i in row["som_all"]:
            a = mol.GetAtomWithIdx(i)
            if a.GetSymbol() != "C":
                continue
            total += 1
            with_h += a.GetTotalNumHs() > 0
    assert total > 2000
    assert with_h / total > 0.95, with_h / total


def test_sulfur_is_enriched_among_labelled_atoms(curated):
    """S-oxidation is a classic CYP reaction, so sulfur should be over-represented."""
    instances, _ = curated
    s_all = s_lab = n_all = n_lab = 0
    for row in instances:
        mol = Chem.MolFromSmiles(row["smiles"])
        labels = set(row["som_all"])
        for a in mol.GetAtoms():
            n_all += 1
            is_s = a.GetSymbol() == "S"
            s_all += is_s
            if a.GetIdx() in labels:
                n_lab += 1
                s_lab += is_s
    enrichment = (s_lab / n_lab) / (s_all / n_all)
    assert enrichment > 2.0, enrichment


def test_oxygen_is_depleted_among_labelled_atoms(curated):
    """O-dealkylation is annotated at the carbon, not at the oxygen."""
    instances, _ = curated
    o_all = o_lab = n_all = n_lab = 0
    for row in instances:
        mol = Chem.MolFromSmiles(row["smiles"])
        labels = set(row["som_all"])
        for a in mol.GetAtoms():
            n_all += 1
            is_o = a.GetSymbol() == "O"
            o_all += is_o
            if a.GetIdx() in labels:
                n_lab += 1
                o_lab += is_o
    enrichment = (o_lab / n_lab) / (o_all / n_all)
    assert enrichment < 0.5, enrichment


def test_molecules_table_covers_every_instance(curated):
    instances, molecules = curated
    for row in instances:
        assert row["smiles"] in molecules
        assert row["isoform"] in molecules[row["smiles"]]["isoforms"]
