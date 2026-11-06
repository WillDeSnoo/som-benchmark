"""Check the Hückel solver against textbook results.

These are the numbers any physical-organic text gives for the standard small
pi systems. They pin down the sign convention, the degenerate-orbital
averaging, and the free-valence definition.
"""

import math

import pytest
from rdkit import Chem

from som.huckel import huckel_atoms, pi_system

SQRT3 = math.sqrt(3.0)


def mol(smiles):
    m = Chem.MolFromSmiles(smiles)
    assert m is not None, smiles
    return m


def test_ethylene():
    recs = huckel_atoms(mol("C=C"))
    assert all(r.in_pi for r in recs)
    for r in recs:
        assert r.frag_size == 2
        assert r.q_pi == pytest.approx(1.0, abs=1e-9)
        assert r.bond_order_sum == pytest.approx(1.0, abs=1e-9)
        assert r.free_valence == pytest.approx(SQRT3 - 1.0, abs=1e-9)
        assert r.f_electrophilic == pytest.approx(0.5, abs=1e-9)
        assert r.f_nucleophilic == pytest.approx(0.5, abs=1e-9)
    # E = alpha + x*beta with beta < 0, so the bonding orbital has x = +1
    assert recs[0].frag_x_homo == pytest.approx(1.0, abs=1e-9)
    assert recs[0].frag_x_lumo == pytest.approx(-1.0, abs=1e-9)
    assert recs[0].frag_gap == pytest.approx(2.0, abs=1e-9)


def test_benzene_eigenvalues_and_symmetry():
    recs = huckel_atoms(mol("c1ccccc1"))
    assert all(r.in_pi for r in recs)
    # +2, +1, +1, -1, -1, -2  =>  HOMO x = +1, LUMO x = -1
    assert recs[0].frag_x_homo == pytest.approx(1.0, abs=1e-8)
    assert recs[0].frag_x_lumo == pytest.approx(-1.0, abs=1e-8)
    for r in recs:
        assert r.q_pi == pytest.approx(1.0, abs=1e-8)
        assert r.bond_order_sum == pytest.approx(4.0 / 3.0, abs=1e-8)
        assert r.free_valence == pytest.approx(SQRT3 - 4.0 / 3.0, abs=1e-8)


def test_benzene_degenerate_homo_is_averaged():
    """The whole point of the degeneracy handling.

    Either HOMO of benzene on its own puts unequal density on the six carbons;
    only the average over the degenerate pair restores the six-fold symmetry.
    """
    recs = huckel_atoms(mol("c1ccccc1"))
    for r in recs:
        assert r.f_electrophilic == pytest.approx(1.0 / 6.0, abs=1e-8)
        assert r.f_nucleophilic == pytest.approx(1.0 / 6.0, abs=1e-8)
    spread = max(r.f_electrophilic for r in recs) - min(
        r.f_electrophilic for r in recs
    )
    assert spread < 1e-8


def test_butadiene_free_valence_terminal_gt_central():
    m = mol("C=CC=C")
    recs = huckel_atoms(m)
    assert all(r.in_pi for r in recs)
    assert recs[0].frag_x_homo == pytest.approx(0.618034, abs=1e-5)
    assert recs[0].frag_x_lumo == pytest.approx(-0.618034, abs=1e-5)
    # canonical SMILES order for C=CC=C is terminal, central, central, terminal
    term = [recs[0], recs[3]]
    cent = [recs[1], recs[2]]
    # bond orders are p12 = p34 = 0.894427, p23 = 0.447214, so the Coulson free
    # valences are sqrt(3) - 0.894427 and sqrt(3) - 1.341641
    for t in term:
        assert t.free_valence == pytest.approx(SQRT3 - 0.894427, abs=1e-5)
    for c in cent:
        assert c.free_valence == pytest.approx(SQRT3 - 1.341641, abs=1e-5)
    assert min(t.free_valence for t in term) > max(c.free_valence for c in cent)
    # terminal positions also carry the larger HOMO density
    assert min(t.f_electrophilic for t in term) > max(
        c.f_electrophilic for c in cent
    )


def test_allyl_cation_and_anion_electron_count():
    """A charged carbon still contributes a p orbital to the pi system."""
    cation = huckel_atoms(mol("[CH2+]C=C"))
    anion = huckel_atoms(mol("[CH2-]C=C"))
    for recs, n_elec in ((cation, 2.0), (anion, 4.0)):
        inpi = [r for r in recs if r.in_pi]
        assert len(inpi) == 3, "all three allyl carbons belong to the pi system"
        assert all(r.frag_size == 3 for r in inpi)
        assert sum(r.q_pi for r in recs) == pytest.approx(n_elec, abs=1e-8)
    # allyl eigenvalues are +sqrt2, 0, -sqrt2; the cation holds 2 electrons so
    # its HOMO is the bonding orbital, the anion holds 4 so its HOMO is the
    # non-bonding one
    assert cation[0].frag_x_homo == pytest.approx(math.sqrt(2.0), abs=1e-6)
    assert anion[0].frag_x_homo == pytest.approx(0.0, abs=1e-6)
    # the allyl frontier orbital has a node at the central carbon
    cent = [
        r for r, a in zip(anion, mol("[CH2-]C=C").GetAtoms())
        if a.GetDegree() == 2
    ]
    assert len(cent) == 1
    assert cent[0].f_electrophilic == pytest.approx(0.0, abs=1e-8)


def test_naphthalene_alpha_more_reactive_than_beta():
    """Hückel reproduces the known alpha-preference of naphthalene."""
    m = mol("c1ccc2ccccc2c1")
    recs = huckel_atoms(m)
    # alpha carbons (positions 1,4,5,8) are the CH carbons adjacent to a
    # fusion carbon; beta carbons are the other CH carbons
    fused = {
        a.GetIdx() for a in m.GetAtoms()
        if a.GetTotalNumHs() == 0 and a.GetIsAromatic()
    }
    alpha, beta = [], []
    for a in m.GetAtoms():
        if a.GetTotalNumHs() != 1:
            continue
        nbrs = {n.GetIdx() for n in a.GetNeighbors()}
        (alpha if nbrs & fused else beta).append(recs[a.GetIdx()])
    assert len(alpha) == 4 and len(beta) == 4
    assert min(r.free_valence for r in alpha) > max(r.free_valence for r in beta)
    assert min(r.f_electrophilic for r in alpha) > max(
        r.f_electrophilic for r in beta
    )


def test_disconnected_pi_systems_solved_separately():
    """Diphenylmethane: both rings must get full frontier density."""
    recs = huckel_atoms(mol("c1ccccc1Cc1ccccc1"))
    aromatic = [r for r in recs if r.in_pi]
    assert len(aromatic) == 12
    assert all(r.frag_size == 6 for r in aromatic)
    # if the fragments had been solved as one matrix, one ring would be ~zero
    assert all(r.f_electrophilic > 0.1 for r in aromatic)


def test_sp3_atoms_are_not_in_the_pi_system():
    m = mol("CCCC")
    members, _ = pi_system(m)
    assert members == set()
    recs = huckel_atoms(m)
    assert not any(r.in_pi for r in recs)
    assert all(r.free_valence == 0.0 for r in recs)


def test_aliphatic_ether_not_swept_into_pi_system():
    """An isolated ether oxygen has no pi system to conjugate with."""
    m = mol("CCOCC")
    members, _ = pi_system(m)
    assert members == set()


def test_anisole_oxygen_is_a_donor_in_the_ring_system():
    m = mol("COc1ccccc1")
    members, kinds = pi_system(m)
    o_idx = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O"][0]
    assert o_idx in members
    assert kinds[o_idx] == "donor"
    recs = huckel_atoms(m)
    assert recs[o_idx].frag_size == 7  # six ring carbons + the oxygen
    # donating oxygen raises pi density at the ortho/para carbons
    ring = [r for i, r in enumerate(recs) if r.in_pi and i != o_idx]
    assert max(r.q_pi for r in ring) > 1.0


def test_pyridine_nitrogen_takes_one_electron():
    m = mol("c1ccncc1")
    recs = huckel_atoms(m)
    n_idx = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "N"][0]
    assert recs[n_idx].in_pi
    assert sum(r.q_pi for r in recs if r.in_pi) == pytest.approx(6.0, abs=1e-8)
    # electronegative N withdraws density from the ring
    assert recs[n_idx].q_pi > 1.0


def test_pyrrole_nitrogen_donates_a_pair():
    m = mol("c1cc[nH]c1")
    recs = huckel_atoms(m)
    assert sum(r.q_pi for r in recs if r.in_pi) == pytest.approx(6.0, abs=1e-8)


def test_every_curated_molecule_solves():
    """Smoke test over the real dataset, if it has been curated."""
    from som.data import load_instances
    from som.paths import PROCESSED

    if not (PROCESSED / "instances.csv").exists():
        pytest.skip("run `make curate` first")
    seen = set()
    for row in load_instances():
        if row["smiles"] in seen:
            continue
        seen.add(row["smiles"])
        m = Chem.MolFromSmiles(row["smiles"])
        recs = huckel_atoms(m)
        assert len(recs) == m.GetNumAtoms()
        for r in recs:
            for name in ("q_pi", "free_valence", "f_electrophilic",
                         "self_polarizability"):
                v = getattr(r, name)
                assert math.isfinite(v), (row["smiles"], name, v)
