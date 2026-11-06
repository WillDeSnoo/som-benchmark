"""Simple Hückel theory, used as a source of cheap frontier-orbital atom descriptors.

This is deliberately the 1950s model: one p-orbital per conjugated atom,
heteroatoms entered through the usual Coulomb/resonance parameters
h_X and k_CX (Streitwieser's values), everything in units of beta.

Energies are written E = alpha + x*beta. Because beta is negative, a *larger*
x means a *more stable* orbital, so the occupied set is the one with the largest
x values and the HOMO is the occupied orbital with the smallest x. Getting this
backwards is the classic way to produce plausible-looking nonsense, so the sign
convention is asserted in the tests.

Two details that matter for correctness:

* **Degenerate frontier orbitals.** In benzene the HOMO is doubly degenerate and
  each individual eigenvector has an arbitrary, uneven distribution; only the
  average over the degenerate pair has the correct six-fold symmetry. Frontier
  densities are therefore averaged over all orbitals within `DEGEN_TOL` of the
  frontier energy.
* **Disconnected pi systems.** A molecule with two aromatic rings joined by a
  CH2 has a block-diagonal Hückel matrix, and a global HOMO localised on just
  one ring would leave the other ring with zero frontier density. Each connected
  pi fragment is therefore solved separately, and the fragment's own frontier
  eigenvalues are reported alongside the per-atom densities so that a model can
  still learn which fragment is the reactive one.

Quantities returned per conjugated atom: pi electron density, Coulson free
valence, electrophilic / nucleophilic / radical frontier densities, atom
self-polarizability, and total pi bond order.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from rdkit import Chem

#: orbitals within this distance (in units of beta) count as degenerate
DEGEN_TOL = 1e-6

SQRT3 = math.sqrt(3.0)

# (h_X, k_CX, n_pi_electrons) for the two common bonding situations of each
# element: "donor" = contributes a lone pair (2 electrons), "pi" = contributes
# one electron through a formal multiple bond.
# Values follow Streitwieser's standard table; P and Se are rough analogues and
# are flagged as approximate in the README.
_PARAMS: dict[tuple[str, str], tuple[float, float, int]] = {
    ("C", "pi"):     (0.0, 1.00, 1),
    ("C", "donor"):  (0.0, 1.00, 2),   # carbanion
    ("N", "pi"):     (0.5, 1.00, 1),   # pyridine-type / imine
    ("N", "donor"):  (1.5, 0.80, 2),   # pyrrole-type / amine / amide
    ("O", "pi"):     (1.0, 1.00, 1),   # carbonyl
    ("O", "donor"):  (2.0, 0.80, 2),   # ether / hydroxyl / furan
    ("S", "pi"):     (1.0, 1.00, 1),   # thione
    ("S", "donor"):  (1.11, 0.69, 2),  # thioether / thiophene
    ("P", "pi"):     (0.4, 0.80, 1),   # approximate
    ("P", "donor"):  (0.8, 0.60, 2),   # approximate
    ("Se", "pi"):    (0.9, 0.90, 1),   # approximate
    ("Se", "donor"): (1.0, 0.60, 2),   # approximate
    ("F", "donor"):  (3.0, 0.70, 2),
    ("Cl", "donor"): (2.0, 0.40, 2),
    ("Br", "donor"): (1.5, 0.30, 2),
    ("I", "donor"):  (1.2, 0.25, 2),
    ("B", "pi"):     (-1.0, 0.70, 0),
}

_HALOGENS = {"F", "Cl", "Br", "I"}


@dataclass
class HuckelAtom:
    """Per-atom Hückel descriptors. All zero/NaN-free; `in_pi` says if it's real."""
    in_pi: bool = False
    q_pi: float = 0.0            # pi electron density
    free_valence: float = 0.0    # Coulson free valence
    bond_order_sum: float = 0.0  # sum of pi bond orders at this atom
    f_electrophilic: float = 0.0  # HOMO density (electrophilic attack)
    f_nucleophilic: float = 0.0   # LUMO density (nucleophilic attack)
    f_radical: float = 0.0        # mean of the two
    self_polarizability: float = 0.0
    frag_x_homo: float = 0.0     # fragment frontier eigenvalues, units of beta
    frag_x_lumo: float = 0.0
    frag_gap: float = 0.0
    frag_size: int = 0


def _classify(atom: Chem.Atom) -> str | None:
    """Return 'pi', 'donor', or None if the atom has no p-orbital in the system."""
    sym = atom.GetSymbol()
    if sym in _HALOGENS:
        # only conjugated if attached to something unsaturated
        return "donor"
    has_multiple = any(
        b.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE)
        and b.GetOtherAtom(atom).GetAtomicNum() > 1
        for b in atom.GetBonds()
    )
    if atom.GetIsAromatic():
        if sym == "C":
            return "pi"
        # aromatic heteroatom: 3 connections (pyrrole/furan) donates a pair,
        # 2 connections (pyridine) contributes one electron
        heavy_deg = atom.GetDegree()
        if sym in ("N", "P"):
            return "donor" if (heavy_deg + atom.GetTotalNumHs()) >= 3 else "pi"
        if sym in ("O", "S", "Se"):
            return "donor"
        return "pi"
    if has_multiple:
        return "pi"
    if sym in ("N", "O", "S", "P", "Se"):
        return "donor"
    if sym == "C":
        # A charged carbon still carries a p orbital: an empty one in a
        # carbocation (the allyl cation is the textbook case) or a filled one
        # in a carbanion. Neutral sp3 carbon has none.
        if atom.GetFormalCharge() > 0:
            return "pi"
        if atom.GetFormalCharge() < 0:
            return "donor"
        return None
    if sym == "B":
        return "pi"
    return None


def pi_system(mol: Chem.Mol) -> tuple[set[int], dict[int, str]]:
    """Identify conjugated atoms.

    Seeds on aromatic atoms and atoms in heavy-atom multiple bonds, then grows
    by one shell to pick up lone-pair donors (ethers, amines, halogens) attached
    to the seed. A donor only joins if it is bonded to a seed atom, which stops
    every aliphatic ether oxygen in the molecule from being swept in.
    """
    seeds: set[int] = set()
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            seeds.add(atom.GetIdx())
            continue
        if any(
            b.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE)
            and b.GetOtherAtom(atom).GetAtomicNum() > 1
            for b in atom.GetBonds()
        ):
            seeds.add(atom.GetIdx())

    members = set(seeds)
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in members:
            continue
        # Atoms that are not seeds but still carry a p orbital: lone-pair
        # donors (ethers, amines, halogens) and empty-p centres (carbocations,
        # boron). They join only if bonded to a seed, which stops every
        # aliphatic heteroatom in the molecule from being swept in.
        if _classify(atom) is None:
            continue
        if any(nb.GetIdx() in seeds for nb in atom.GetNeighbors()):
            members.add(idx)

    kinds: dict[int, str] = {}
    for idx in list(members):
        kind = _classify(mol.GetAtomWithIdx(idx))
        if kind is None:
            members.discard(idx)
        else:
            kinds[idx] = kind
    return members, kinds


def _params(atom: Chem.Atom, kind: str) -> tuple[float, float, int]:
    sym = atom.GetSymbol()
    key = (sym, kind)
    if key not in _PARAMS:
        key = (sym, "donor" if ("donor") in kind else "pi")
    if key not in _PARAMS:
        return (0.0, 1.0, 1)  # fall back to a carbon-like centre
    h, k, n_e = _PARAMS[key]
    # formal charge shifts the electron count: a cation has one fewer p electron
    n_e = int(min(2, max(0, n_e - atom.GetFormalCharge())))
    return h, k, n_e


def _fragments(mol: Chem.Mol, members: set[int]) -> list[list[int]]:
    """Connected components of the pi system (bonds between two members)."""
    adj: dict[int, set[int]] = {i: set() for i in members}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in members and b in members:
            adj[a].add(b)
            adj[b].add(a)
    seen: set[int] = set()
    out = []
    for start in sorted(members):
        if start in seen:
            continue
        comp, stack = [], [start]
        seen.add(start)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nb in adj[cur]:
                if nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        out.append(sorted(comp))
    return out


def _frontier_density(
    coef: np.ndarray, x: np.ndarray, level: int, n: int
) -> np.ndarray:
    """Mean |c|^2 over orbitals degenerate with orbital `level`.

    `coef` columns are eigenvectors ordered by descending x. Returns zeros if
    the requested level does not exist (e.g. LUMO of a fully filled fragment of
    size 1).
    """
    if level < 0 or level >= n:
        return np.zeros(n)
    target = x[level]
    cols = np.flatnonzero(np.abs(x - target) <= DEGEN_TOL)
    return (coef[:, cols] ** 2).mean(axis=1)


def huckel_atoms(mol: Chem.Mol) -> list[HuckelAtom]:
    """Solve simple Hückel per pi fragment; return one record per atom."""
    out = [HuckelAtom() for _ in range(mol.GetNumAtoms())]
    members, kinds = pi_system(mol)
    if not members:
        return out

    for frag in _fragments(mol, members):
        n = len(frag)
        pos = {idx: i for i, idx in enumerate(frag)}
        h = np.zeros(n)
        kself = np.zeros(n)
        n_elec = 0
        for idx in frag:
            atom = mol.GetAtomWithIdx(idx)
            hi, ki, ne = _params(atom, kinds[idx])
            h[pos[idx]] = hi
            kself[pos[idx]] = ki
            n_elec += ne

        A = np.diag(h)
        bonds = []
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a in pos and b in pos:
                ia, ib = pos[a], pos[b]
                sa = mol.GetAtomWithIdx(a).GetSymbol()
                sb = mol.GetAtomWithIdx(b).GetSymbol()
                if sa == "C" and sb == "C":
                    k = 1.0
                elif sa == "C":
                    k = kself[ib]
                elif sb == "C":
                    k = kself[ia]
                else:
                    # two heteroatoms bonded directly: geometric mean of their
                    # C-X resonance integrals, an approximation
                    k = math.sqrt(max(kself[ia] * kself[ib], 1e-12))
                A[ia, ib] = A[ib, ia] = k
                bonds.append((ia, ib))

        # eigh returns ascending; flip so column 0 is the most stable orbital
        x, coef = np.linalg.eigh(A)
        order = np.argsort(x)[::-1]
        x, coef = x[order], coef[:, order]

        # fill from the most stable orbital
        occ = np.zeros(n)
        remaining = min(n_elec, 2 * n)
        for i in range(n):
            if remaining <= 0:
                break
            occ[i] = min(2, remaining)
            remaining -= occ[i]
        n_occ = int(np.count_nonzero(occ))
        homo = n_occ - 1
        lumo = n_occ if n_occ < n else -1

        q = (coef**2 * occ).sum(axis=1)
        f_e = _frontier_density(coef, x, homo, n)
        f_n = _frontier_density(coef, x, lumo, n) if lumo >= 0 else np.zeros(n)

        # pi bond orders and Coulson free valence
        bo_sum = np.zeros(n)
        for ia, ib in bonds:
            p = float((occ * coef[ia, :] * coef[ib, :]).sum())
            bo_sum[ia] += p
            bo_sum[ib] += p

        # atom self-polarizability, in units of 1/beta
        occ_i = np.flatnonzero(occ > 0)
        unocc_i = np.flatnonzero(occ < 2)
        pol = np.zeros(n)
        for r in occ_i:
            for s in unocc_i:
                if r == s:
                    continue
                dx = x[r] - x[s]
                if abs(dx) < 1e-9:
                    continue
                pol += 4.0 * (coef[:, r] ** 2) * (coef[:, s] ** 2) / dx

        x_homo = float(x[homo]) if homo >= 0 else 0.0
        x_lumo = float(x[lumo]) if lumo >= 0 else 0.0
        for idx in frag:
            i = pos[idx]
            rec = out[idx]
            rec.in_pi = True
            rec.q_pi = float(q[i])
            rec.bond_order_sum = float(bo_sum[i])
            rec.free_valence = float(SQRT3 - bo_sum[i])
            rec.f_electrophilic = float(f_e[i])
            rec.f_nucleophilic = float(f_n[i])
            rec.f_radical = float(0.5 * (f_e[i] + f_n[i]))
            rec.self_polarizability = float(pol[i])
            rec.frag_x_homo = x_homo
            rec.frag_x_lumo = x_lumo
            rec.frag_gap = float(x_homo - x_lumo) if lumo >= 0 else 0.0
            rec.frag_size = n

    return out


#: field order used when flattening into a feature matrix
HUCKEL_FIELDS = [
    "in_pi", "q_pi", "free_valence", "bond_order_sum",
    "f_electrophilic", "f_nucleophilic", "f_radical",
    "self_polarizability", "frag_x_homo", "frag_x_lumo", "frag_gap", "frag_size",
]


def as_vector(rec: HuckelAtom) -> list[float]:
    return [float(getattr(rec, f)) for f in HUCKEL_FIELDS]
