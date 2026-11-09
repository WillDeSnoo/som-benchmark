"""Keep docs/HUCKEL_DESCRIPTORS.md in step with the code.

A reference table of feature names is the kind of document that is correct on
the day it is written and quietly wrong six months later. These tests fail if
a descriptor is added, removed or renamed without the document following, and
if a Streitwieser parameter is changed without its row being updated.
"""

import re
from pathlib import Path

import pytest
from rdkit import Chem

from som.features import molecule_features
from som.huckel import _PARAMS

DOC = Path(__file__).resolve().parent.parent / "docs" / "HUCKEL_DESCRIPTORS.md"

#: a molecule with aromatic, aliphatic, heteroatom and halogen environments,
#: so that every feature column is exercised
PROBE = "COc1ccccc1CN(C)CCSc1ccccc1F"


def documented_names(text):
    """Feature names in the doc's tables: the `code`-quoted first cells."""
    found = set()
    for line in text.splitlines():
        m = re.match(r"^\|\s*`([a-z0-9_]+)`\s*\|", line)
        if m:
            found.add(m.group(1))
    return found


@pytest.fixture(scope="module")
def doc_text():
    assert DOC.exists(), f"{DOC} is missing"
    return DOC.read_text()


@pytest.fixture(scope="module")
def huckel_feature_names():
    mol = Chem.MolFromSmiles(PROBE)
    _, names = molecule_features(mol, "phys")
    return {
        n for n in names
        if "huckel" in n.lower() or n == "n_pi_neighbours"
    }


def test_every_huckel_feature_is_documented(doc_text, huckel_feature_names):
    missing = sorted(huckel_feature_names - documented_names(doc_text))
    assert not missing, (
        f"these feature columns exist in code but are not in {DOC.name}: "
        f"{missing}"
    )


def test_no_documented_feature_has_been_removed(doc_text, huckel_feature_names):
    documented = documented_names(doc_text)
    # the parameter table's rows are elements, not features; only names that
    # look like feature columns are compared
    feature_like = {d for d in documented
                    if "huckel" in d or d == "n_pi_neighbours"}
    stale = sorted(feature_like - huckel_feature_names)
    assert not stale, (
        f"{DOC.name} documents feature columns that no longer exist: {stale}"
    )


def test_the_documented_count_matches(doc_text, huckel_feature_names):
    """The prose says how many there are; that number has to be right too."""
    m = re.search(r"\*\*(\d+) of the (\d+) features in `phys`", doc_text)
    assert m, "the doc should state how many features are Huckel-derived"
    claimed_huckel, claimed_total = int(m.group(1)), int(m.group(2))

    mol = Chem.MolFromSmiles(PROBE)
    _, names = molecule_features(mol, "phys")
    assert claimed_huckel == len(huckel_feature_names)
    assert claimed_total == len(names)


def test_streitwieser_table_matches_the_code(doc_text):
    """Every (atom, mode) in the code appears in the doc with the same h, k, n."""
    rows = {}
    for line in doc_text.splitlines():
        m = re.match(
            r"^\|\s*([A-Z][a-z]?)\s*\|\s*(pi|donor)\s*\|\s*(−?[\d.]+)\s*\|"
            r"\s*([\d.]+)\s*\|\s*(\d+)\s*\|",
            line,
        )
        if m:
            sym, mode, h, k, n = m.groups()
            rows[(sym, mode)] = (float(h.replace("−", "-")),
                                 float(k), int(n))

    assert set(rows) == set(_PARAMS), (
        f"doc/code mismatch in which atoms are parameterised: "
        f"only in doc {sorted(set(rows) - set(_PARAMS))}, "
        f"only in code {sorted(set(_PARAMS) - set(rows))}"
    )
    for key, coded in _PARAMS.items():
        assert rows[key] == pytest.approx(coded), (
            f"{key} is {coded} in code but {rows[key]} in {DOC.name}"
        )
