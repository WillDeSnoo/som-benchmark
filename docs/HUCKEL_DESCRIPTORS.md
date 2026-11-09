# The Hückel descriptor set

Reference for the Hückel-derived features in the `phys` feature set, produced
by [`src/som/huckel.py`](../src/som/huckel.py) and assembled in
[`src/som/features.py`](../src/som/features.py).

**21 of the 55 features in `phys` are Hückel-derived**: 12 computed on the atom
itself, 8 aggregated over its heavy neighbours, plus a neighbour count. The
other 34 are the basic atom and bond features that the plain GNN also receives,
so `phys` = `graph` + this file, strictly additively. That is asserted in
`tests/test_features.py`, because the `gnn` against `gnn_huckel` ablation is
only meaningful if it holds.

`tests/test_huckel_docs.py` asserts that the tables below list exactly the
Hückel feature columns the code emits, so this document cannot drift from the
implementation.

## The method

One p orbital per conjugated atom, energies as **E = α + xβ**, where α is the
Coulomb integral for carbon and β the carbon–carbon resonance integral. β is
negative, so a *larger* x means a *more stable* orbital. Heteroatoms enter
through the standard Streitwieser substitutions

    α_X = α + h_X·β        β_CX = k_CX·β

Two implementation details that are easy to get wrong, and are tested:

- **Degenerate frontier orbitals.** Benzene's HOMO is doubly degenerate, and
  each individual eigenvector has an arbitrary, lopsided distribution. Only the
  average over the degenerate pair restores the six-fold symmetry. Frontier
  densities are averaged over all orbitals within `DEGEN_TOL` (1e-6) of the
  frontier energy. Without this, benzene's six equivalent carbons get six
  different reactivities.
- **Disconnected π systems.** Diphenylmethane's Hückel matrix is block
  diagonal. A single global HOMO would localise on one ring and leave the other
  with zero frontier density, which is chemically wrong. Each connected π
  fragment is solved separately, and the fragment's own frontier eigenvalues
  are carried as features so a model can still learn which fragment is the
  reactive one.

## Computed on the atom

| Feature | Meaning |
| --- | --- |
| `huckel_in_pi` | 1 if the atom carries a p orbital in some π system, else 0. Zero for every sp³ carbon, which is why the neighbour block exists. |
| `huckel_q_pi` | π electron density on the atom, Σ over occupied orbitals of n·c². Electron-rich positions are the ones an electrophilic oxidant prefers. |
| `huckel_free_valence` | Coulson free valence, √3 − Σ(π bond orders at the atom). How much π bonding capacity is unused; the classical index for radical reactivity. |
| `huckel_bond_order_sum` | Σ of π bond orders to the atom's neighbours, the quantity free valence is derived from. |
| `huckel_f_electrophilic` | Frontier density from the HOMO. Susceptibility to electrophilic attack — the relevant one here, since Compound I is electrophilic. |
| `huckel_f_nucleophilic` | Frontier density from the LUMO. Kept for completeness and as a control. |
| `huckel_f_radical` | Mean of the electrophilic and nucleophilic frontier densities; susceptibility to radical attack. |
| `huckel_self_polarizability` | π_ii, how much the atom's own charge responds to a perturbation at itself. A second-order reactivity measure that does not vanish where frontier density does. |
| `huckel_frag_x_homo` | HOMO eigenvalue of this atom's π fragment, as x in E = α + xβ. Smaller x means a higher-lying HOMO, i.e. more easily oxidised. |
| `huckel_frag_x_lumo` | LUMO eigenvalue of the same fragment. |
| `huckel_frag_gap` | HOMO–LUMO gap of the fragment. |
| `huckel_frag_size` | Number of atoms in the fragment; distinguishes an isolated alkene from an extended aromatic. |

## Aggregated over heavy neighbours

| Feature | Meaning |
| --- | --- |
| `nbr_max_huckel_f_electrophilic` | Largest HOMO frontier density among the atom's heavy neighbours. |
| `nbr_mean_huckel_f_electrophilic` | Mean of the same. |
| `nbr_max_huckel_f_radical` | Largest neighbour radical frontier density. |
| `nbr_mean_huckel_f_radical` | Mean of the same. |
| `nbr_max_huckel_free_valence` | Largest neighbour free valence. |
| `nbr_mean_huckel_free_valence` | Mean of the same. |
| `nbr_max_huckel_q_pi` | Largest neighbour π electron density. |
| `nbr_mean_huckel_q_pi` | Mean of the same. |
| `n_pi_neighbours` | How many heavy neighbours are part of a π system. |

**Why this block exists**, and it is the most consequential choice in the
feature set. An sp³ carbon has no π system of its own, so every per-atom
Hückel quantity above is identically zero for it. But benzylic and allylic
positions are among the most metabolically vulnerable sites there are, exactly
because the radical formed on hydrogen abstraction is stabilised by the
*adjacent* π system. Propagating the neighbours' values onto the atom is what
lets a descriptor set built on π theory say anything at all about saturated
carbons.

It works. In the readable linear model
(`results/phys_linear_coefficients_scaffold.csv`, fitted only to be read and
never used for predictions), `huckel_free_valence` and
`nbr_max_huckel_free_valence` are both among the largest positive
coefficients.

Message passing over the molecular graph could in principle derive these
neighbour aggregates for itself, so including them means the GNN and the GBM
see an identical feature set and the comparison is about the estimator rather
than the inputs.

## Streitwieser parameters

`h` shifts the atom's Coulomb integral, `k` scales the resonance integral to
carbon, `n` is the number of π electrons the atom contributes. "pi" is the
one-electron case (pyridine-type N, carbonyl O); "donor" is the lone-pair case
(pyrrole-type N, ether O, halogen).

| Atom | Mode | h | k | n | Typical case |
| --- | --- | --- | --- | --- | --- |
| C | pi | 0.0 | 1.00 | 1 | ordinary sp² carbon |
| C | donor | 0.0 | 1.00 | 2 | carbanion |
| N | pi | 0.5 | 1.00 | 1 | pyridine-type, imine |
| N | donor | 1.5 | 0.80 | 2 | pyrrole-type, amine, amide |
| O | pi | 1.0 | 1.00 | 1 | carbonyl |
| O | donor | 2.0 | 0.80 | 2 | ether, hydroxyl, furan |
| S | pi | 1.0 | 1.00 | 1 | thione |
| S | donor | 1.11 | 0.69 | 2 | thioether, thiophene |
| P | pi | 0.4 | 0.80 | 1 | — |
| P | donor | 0.8 | 0.60 | 2 | — |
| Se | pi | 0.9 | 0.90 | 1 | — |
| Se | donor | 1.0 | 0.60 | 2 | — |
| B | pi | −1.0 | 0.70 | 0 | empty p orbital |
| F | donor | 3.0 | 0.70 | 2 | — |
| Cl | donor | 2.0 | 0.40 | 2 | — |
| Br | donor | 1.5 | 0.30 | 2 | — |
| I | donor | 1.2 | 0.25 | 2 | — |

A charged carbon still carries a p orbital: C⁺ is treated as `pi` with an empty
orbital (the allyl cation is the textbook case) and C⁻ as `donor`.

**The P and Se values are rough analogues rather than fitted parameters.** They
are the least trustworthy numbers here and are flagged as such in the code.

## What is tested

`tests/test_huckel.py` checks the implementation against results that can be
derived by hand: ethylene and benzene orbital energies, butadiene bond orders
and free valences, the allyl cation and anion electron counts and the node at
the central carbon, pyridine versus pyrrole nitrogen, and naphthalene's known
α-over-β preference.
