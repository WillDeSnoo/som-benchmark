# Generated results

## Curation

- records read from the upstream SDFs: **2003**
- excluded: **0**
- duplicate (molecule, isoform) instances merged: **0** (label disagreements: 0)
- surviving instances: **2003** over **679** unique molecules
- Bemis-Murcko scaffolds: **390** (+45 acyclic molecules)
- candidate atoms: **41878**, of which **4163** are labelled sites (9.94%)
- primary-only positives: **2853** (6.81%)

## Headline results

Mean over repeated splits, with the 95% CI of the mean.

### random split

| approach | top-1 | top-2 | ROC-AUC | avg precision | ECE |
|---|---|---|---|---|---|
| atom-type frequency | 0.419 [0.392, 0.446] | 0.597 [0.582, 0.612] | 0.799 [0.786, 0.812] | 0.278 [0.263, 0.292] | 0.023 |
| SMARTCyp (published method) | 0.629 [0.598, 0.660] | 0.760 [0.729, 0.792] | 0.822 [0.811, 0.833] | 0.340 [0.321, 0.359] | 0.016 |
| fingerprint + GBM | 0.622 [0.586, 0.657] | 0.745 [0.716, 0.775] | 0.879 [0.868, 0.889] | 0.523 [0.495, 0.550] | 0.023 |
| Hückel physics + GBM | 0.684 [0.664, 0.703] | 0.822 [0.794, 0.849] | 0.901 [0.891, 0.911] | 0.589 [0.558, 0.620] | 0.035 |
| Hückel + BDE + GBM | 0.706 [0.680, 0.732] | 0.828 [0.806, 0.849] | 0.916 [0.907, 0.925] | 0.631 [0.602, 0.661] | 0.024 |
| message-passing GNN | 0.667 [0.638, 0.696] | 0.794 [0.759, 0.830] | 0.882 [0.868, 0.897] | 0.543 [0.510, 0.577] | 0.023 |
| GNN + Hückel features | 0.667 [0.637, 0.698] | 0.779 [0.746, 0.813] | 0.884 [0.870, 0.898] | 0.545 [0.512, 0.579] | 0.031 |
| GNN + BDE | 0.658 [0.635, 0.681] | 0.787 [0.763, 0.812] | 0.883 [0.867, 0.898] | 0.550 [0.510, 0.590] | 0.029 |

### scaffold split

| approach | top-1 | top-2 | ROC-AUC | avg precision | ECE |
|---|---|---|---|---|---|
| atom-type frequency | 0.425 [0.402, 0.448] | 0.607 [0.583, 0.631] | 0.802 [0.790, 0.814] | 0.274 [0.261, 0.287] | 0.025 |
| SMARTCyp (published method) | 0.635 [0.607, 0.663] | 0.770 [0.744, 0.795] | 0.826 [0.812, 0.840] | 0.341 [0.311, 0.370] | 0.018 |
| fingerprint + GBM | 0.601 [0.568, 0.635] | 0.741 [0.721, 0.760] | 0.871 [0.865, 0.878] | 0.499 [0.477, 0.522] | 0.018 |
| Hückel physics + GBM | 0.660 [0.626, 0.695] | 0.811 [0.787, 0.834] | 0.887 [0.876, 0.899] | 0.531 [0.499, 0.562] | 0.019 |
| Hückel + BDE + GBM | 0.674 [0.650, 0.698] | 0.812 [0.789, 0.834] | 0.907 [0.899, 0.915] | 0.584 [0.552, 0.616] | 0.019 |
| message-passing GNN | 0.642 [0.611, 0.673] | 0.785 [0.762, 0.808] | 0.873 [0.862, 0.884] | 0.501 [0.453, 0.548] | 0.032 |
| GNN + Hückel features | 0.644 [0.618, 0.671] | 0.768 [0.747, 0.789] | 0.872 [0.863, 0.880] | 0.495 [0.457, 0.532] | 0.030 |
| GNN + BDE | 0.623 [0.583, 0.663] | 0.770 [0.747, 0.794] | 0.877 [0.867, 0.886] | 0.513 [0.465, 0.561] | 0.027 |

## Random minus scaffold

| approach | metric | random | scaffold | gap [95% CI] |
|---|---|---|---|---|
| atom-type frequency | top1 | 0.419 | 0.425 | -0.006 [-0.037, 0.025] |
| atom-type frequency | top2 | 0.597 | 0.607 | -0.009 [-0.034, 0.015] |
| SMARTCyp (published method) | top1 | 0.629 | 0.635 | -0.006 [-0.042, 0.030] |
| SMARTCyp (published method) | top2 | 0.760 | 0.770 | -0.009 [-0.044, 0.026] |
| fingerprint + GBM | top1 | 0.622 | 0.601 | 0.021 [-0.022, 0.063] |
| fingerprint + GBM | top2 | 0.745 | 0.741 | 0.005 [-0.026, 0.035] |
| Hückel physics + GBM | top1 | 0.684 | 0.660 | 0.024 [-0.011, 0.058] |
| Hückel physics + GBM | top2 | 0.822 | 0.811 | 0.011 [-0.020, 0.043] |
| Hückel + BDE + GBM | top1 | 0.706 | 0.674 | 0.032 [0.001, 0.063] |
| Hückel + BDE + GBM | top2 | 0.828 | 0.812 | 0.016 [-0.011, 0.043] |
| message-passing GNN | top1 | 0.667 | 0.642 | 0.025 [-0.011, 0.062] |
| message-passing GNN | top2 | 0.794 | 0.785 | 0.009 [-0.027, 0.046] |
| GNN + Hückel features | top1 | 0.667 | 0.644 | 0.023 [-0.012, 0.058] |
| GNN + Hückel features | top2 | 0.779 | 0.768 | 0.011 [-0.023, 0.045] |
| GNN + BDE | top1 | 0.658 | 0.623 | 0.035 [-0.006, 0.075] |
| GNN + BDE | top2 | 0.787 | 0.770 | 0.017 [-0.012, 0.047] |

## Recall at top-2 by site type (scaffold split)

| site type | n atoms | atom-type frequency | SMARTCyp (published method) | fingerprint + GBM | Hückel physics + GBM | Hückel + BDE + GBM | message-passing GNN | GNN + Hückel features | GNN + BDE |
|---|---|---|---|---|---|---|---|---|---|
| aromatic C-H hydroxylation | 2412 | 0.103 | 0.256 | 0.270 | 0.358 | 0.329 | 0.312 | 0.308 | 0.269 |
| N-dealkylation (alpha to N) | 2221 | 0.487 | 0.751 | 0.692 | 0.752 | 0.773 | 0.717 | 0.718 | 0.732 |
| O-dealkylation (alpha to O) | 1254 | 0.618 | 0.676 | 0.727 | 0.748 | 0.754 | 0.778 | 0.731 | 0.772 |
| benzylic / allylic hydroxylation | 1082 | 0.378 | 0.539 | 0.438 | 0.435 | 0.493 | 0.449 | 0.444 | 0.433 |
| S-oxidation | 435 | 0.910 | 0.662 | 0.749 | 0.703 | 0.736 | 0.690 | 0.634 | 0.703 |
| aliphatic C-H hydroxylation | 365 | 0.184 | 0.321 | 0.340 | 0.419 | 0.403 | 0.427 | 0.408 | 0.395 |
| N-oxidation | 319 | 0.210 | 0.407 | 0.389 | 0.342 | 0.345 | 0.320 | 0.298 | 0.364 |
| methyl hydroxylation | 290 | 0.665 | 0.083 | 0.241 | 0.190 | 0.217 | 0.214 | 0.214 | 0.193 |
| O-site (direct) | 91 | 0.000 | 0.000 | 0.077 | 0.143 | 0.165 | 0.055 | 0.143 | 0.099 |
| P-oxidation / desulfuration | 36 | 1.000 | 0.000 | 0.667 | 0.667 | 0.667 | 0.917 | 0.889 | 1.000 |
| oxidative dehalogenation | 31 | 0.903 | 0.290 | 0.903 | 0.903 | 0.903 | 0.903 | 0.871 | 0.903 |
| aromatic C (substituted) | 22 | 0.000 | 0.000 | 0.000 | 0.000 | 0.045 | 0.000 | 0.000 | 0.000 |

## Top-2 by isoform (scaffold split)

| isoform | n test inst. | atom-type frequency | SMARTCyp (published method) | fingerprint + GBM | Hückel physics + GBM | Hückel + BDE + GBM | message-passing GNN | GNN + Hückel features | GNN + BDE |
|---|---|---|---|---|---|---|---|---|---|
| CYP1A2 | 59 | 0.663 | 0.751 | 0.730 | 0.796 | 0.806 | 0.787 | 0.769 | 0.783 |
| CYP2A6 | 24 | 0.658 | 0.811 | 0.673 | 0.796 | 0.818 | 0.766 | 0.768 | 0.734 |
| CYP2B6 | 34 | 0.617 | 0.734 | 0.735 | 0.801 | 0.772 | 0.770 | 0.724 | 0.722 |
| CYP2C19 | 46 | 0.601 | 0.794 | 0.777 | 0.841 | 0.829 | 0.771 | 0.744 | 0.776 |
| CYP2C8 | 29 | 0.536 | 0.790 | 0.701 | 0.819 | 0.831 | 0.760 | 0.762 | 0.754 |
| CYP2C9 | 48 | 0.625 | 0.808 | 0.775 | 0.842 | 0.845 | 0.764 | 0.758 | 0.761 |
| CYP2D6 | 60 | 0.590 | 0.806 | 0.757 | 0.848 | 0.842 | 0.845 | 0.826 | 0.797 |
| CYP2E1 | 33 | 0.686 | 0.775 | 0.735 | 0.798 | 0.775 | 0.754 | 0.780 | 0.768 |
| CYP3A4 | 98 | 0.554 | 0.726 | 0.734 | 0.780 | 0.790 | 0.793 | 0.762 | 0.776 |

## Top-1 by similarity to the nearest training molecule

| NN Tanimoto | n test inst. | atom-type frequency | SMARTCyp (published method) | fingerprint + GBM | Hückel physics + GBM | Hückel + BDE + GBM | message-passing GNN | GNN + Hückel features | GNN + BDE |
|---|---|---|---|---|---|---|---|---|---|
| <=0.30 | 247 | 0.371 | 0.585 | 0.523 | 0.606 | 0.626 | 0.585 | 0.581 | 0.573 |
| 0.30-0.45 | 239 | 0.482 | 0.666 | 0.645 | 0.699 | 0.729 | 0.687 | 0.694 | 0.695 |
| 0.45-0.60 | 160 | 0.443 | 0.663 | 0.644 | 0.672 | 0.680 | 0.679 | 0.674 | 0.647 |
| 0.60-0.75 | 145 | 0.417 | 0.637 | 0.675 | 0.713 | 0.734 | 0.691 | 0.688 | 0.654 |
| >0.75 | 49 | 0.318 | 0.546 | 0.573 | 0.763 | 0.719 | 0.623 | 0.674 | 0.635 |

## Paired model comparisons

Paired t-tests on the per-split differences. Four metrics are tested per pair, so a single p just under 0.05 should be read with that multiplicity in mind.

**phys_bde_gbm − smartcyp** — best model vs the published external method (paired over 20 splits)

| metric | difference [95% CI] | p |
|---|---|---|
| top1 | +0.0578 [+0.0361, +0.0796] | 0.0000 |
| top2 | +0.0549 [+0.0398, +0.0699] | 0.0000 |
| roc_auc | +0.0875 [+0.0798, +0.0952] | 0.0000 |
| avg_precision | +0.2671 [+0.2468, +0.2875] | 0.0000 |

**smartcyp − prior** — does the published method beat the trivial reference? (paired over 20 splits)

| metric | difference [95% CI] | p |
|---|---|---|
| top1 | +0.2100 [+0.1852, +0.2348] | 0.0000 |
| top2 | +0.1628 [+0.1418, +0.1837] | 0.0000 |
| roc_auc | +0.0236 [+0.0152, +0.0319] | 0.0000 |
| avg_precision | +0.0646 [+0.0465, +0.0827] | 0.0000 |

**phys_bde_gbm − phys_gbm** — does adding transferred DFT bond dissociation energies help? (paired over 20 splits)

| metric | difference [95% CI] | p |
|---|---|---|
| top1 | +0.0179 [+0.0061, +0.0297] | 0.0050 |
| top2 | +0.0036 [-0.0051, +0.0123] | 0.4027 |
| roc_auc | +0.0174 [+0.0125, +0.0222] | 0.0000 |
| avg_precision | +0.0475 [+0.0375, +0.0575] | 0.0000 |

**gnn_bde − gnn** — same question, for the GNN (paired over 20 splits)

| metric | difference [95% CI] | p |
|---|---|---|
| top1 | -0.0136 [-0.0317, +0.0044] | 0.1310 |
| top2 | -0.0108 [-0.0295, +0.0079] | 0.2430 |
| roc_auc | +0.0022 [-0.0053, +0.0097] | 0.5428 |
| avg_precision | +0.0096 [-0.0079, +0.0271] | 0.2661 |

**gnn_huckel − gnn** — does giving the GNN the Hückel descriptors help? (paired over 20 splits)

| metric | difference [95% CI] | p |
|---|---|---|
| top1 | +0.0017 [-0.0121, +0.0155] | 0.8012 |
| top2 | -0.0157 [-0.0305, -0.0009] | 0.0393 |
| roc_auc | +0.0003 [-0.0070, +0.0075] | 0.9390 |
| avg_precision | -0.0020 [-0.0189, +0.0149] | 0.8083 |

**phys_gbm − gnn_huckel** — same features, different model class (paired over 20 splits)

| metric | difference [95% CI] | p |
|---|---|---|
| top1 | +0.0160 [-0.0050, +0.0370] | 0.1267 |
| top2 | +0.0423 [+0.0296, +0.0550] | 0.0000 |
| roc_auc | +0.0164 [+0.0112, +0.0216] | 0.0000 |
| avg_precision | +0.0399 [+0.0291, +0.0507] | 0.0000 |

**phys_gbm − fp_gbm** — physics descriptors vs hashed environments (paired over 20 splits)

| metric | difference [95% CI] | p |
|---|---|---|
| top1 | +0.0606 [+0.0448, +0.0764] | 0.0000 |
| top2 | +0.0732 [+0.0597, +0.0866] | 0.0000 |
| roc_auc | +0.0192 [+0.0137, +0.0248] | 0.0000 |
| avg_precision | +0.0489 [+0.0334, +0.0644] | 0.0000 |

**phys_gbm − prior** — best model vs the frequency reference (paired over 20 splits)

| metric | difference [95% CI] | p |
|---|---|---|
| top1 | +0.2500 [+0.2262, +0.2738] | 0.0000 |
| top2 | +0.2141 [+0.1939, +0.2342] | 0.0000 |
| roc_auc | +0.0937 [+0.0864, +0.1009] | 0.0000 |
| avg_precision | +0.2842 [+0.2622, +0.3063] | 0.0000 |

## Split diagnostics (seed 0)

| mode | test molecules | shared scaffolds | median NN Tanimoto | frac. test with NN > 0.6 |
|---|---|---|---|---|
| random | 136 | 34 | 0.491 | 0.368 |
| scaffold | 136 | 0 | 0.343 | 0.103 |

