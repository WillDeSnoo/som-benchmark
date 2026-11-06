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
| fingerprint + GBM | 0.622 [0.586, 0.657] | 0.745 [0.716, 0.775] | 0.879 [0.868, 0.889] | 0.523 [0.495, 0.550] | 0.023 |
| Hückel physics + GBM | 0.684 [0.664, 0.703] | 0.822 [0.794, 0.849] | 0.901 [0.891, 0.911] | 0.589 [0.558, 0.620] | 0.035 |
| message-passing GNN | 0.667 [0.638, 0.696] | 0.794 [0.759, 0.830] | 0.882 [0.868, 0.897] | 0.543 [0.510, 0.577] | 0.023 |
| GNN + Hückel features | 0.667 [0.637, 0.698] | 0.779 [0.746, 0.813] | 0.884 [0.870, 0.898] | 0.545 [0.512, 0.579] | 0.031 |

### scaffold split

| approach | top-1 | top-2 | ROC-AUC | avg precision | ECE |
|---|---|---|---|---|---|
| atom-type frequency | 0.425 [0.402, 0.448] | 0.607 [0.583, 0.631] | 0.802 [0.790, 0.814] | 0.274 [0.261, 0.287] | 0.025 |
| fingerprint + GBM | 0.601 [0.568, 0.635] | 0.741 [0.721, 0.760] | 0.871 [0.865, 0.878] | 0.499 [0.477, 0.522] | 0.018 |
| Hückel physics + GBM | 0.660 [0.626, 0.695] | 0.811 [0.787, 0.834] | 0.887 [0.876, 0.899] | 0.531 [0.499, 0.562] | 0.019 |
| message-passing GNN | 0.642 [0.611, 0.673] | 0.785 [0.762, 0.808] | 0.873 [0.862, 0.884] | 0.501 [0.453, 0.548] | 0.032 |
| GNN + Hückel features | 0.644 [0.618, 0.671] | 0.768 [0.747, 0.789] | 0.872 [0.863, 0.880] | 0.495 [0.457, 0.532] | 0.030 |

## Random minus scaffold

| approach | metric | random | scaffold | gap [95% CI] |
|---|---|---|---|---|
| atom-type frequency | top1 | 0.419 | 0.425 | -0.006 [-0.037, 0.025] |
| atom-type frequency | top2 | 0.597 | 0.607 | -0.009 [-0.034, 0.015] |
| fingerprint + GBM | top1 | 0.622 | 0.601 | 0.021 [-0.022, 0.063] |
| fingerprint + GBM | top2 | 0.745 | 0.741 | 0.005 [-0.026, 0.035] |
| Hückel physics + GBM | top1 | 0.684 | 0.660 | 0.024 [-0.011, 0.058] |
| Hückel physics + GBM | top2 | 0.822 | 0.811 | 0.011 [-0.020, 0.043] |
| message-passing GNN | top1 | 0.667 | 0.642 | 0.025 [-0.011, 0.062] |
| message-passing GNN | top2 | 0.794 | 0.785 | 0.009 [-0.027, 0.046] |
| GNN + Hückel features | top1 | 0.667 | 0.644 | 0.023 [-0.012, 0.058] |
| GNN + Hückel features | top2 | 0.779 | 0.768 | 0.011 [-0.023, 0.045] |

## Recall at top-2 by site type (scaffold split)

| site type | n atoms | atom-type frequency | fingerprint + GBM | Hückel physics + GBM | message-passing GNN | GNN + Hückel features |
|---|---|---|---|---|---|---|
| aromatic C-H hydroxylation | 2412 | 0.103 | 0.270 | 0.358 | 0.312 | 0.308 |
| N-dealkylation (alpha to N) | 2221 | 0.487 | 0.692 | 0.752 | 0.717 | 0.718 |
| O-dealkylation (alpha to O) | 1254 | 0.618 | 0.727 | 0.748 | 0.778 | 0.731 |
| benzylic / allylic hydroxylation | 1082 | 0.378 | 0.438 | 0.435 | 0.449 | 0.444 |
| S-oxidation | 435 | 0.910 | 0.749 | 0.703 | 0.690 | 0.634 |
| aliphatic C-H hydroxylation | 365 | 0.184 | 0.340 | 0.419 | 0.427 | 0.408 |
| N-oxidation | 319 | 0.210 | 0.389 | 0.342 | 0.320 | 0.298 |
| methyl hydroxylation | 290 | 0.665 | 0.241 | 0.190 | 0.214 | 0.214 |
| O-site (direct) | 91 | 0.000 | 0.077 | 0.143 | 0.055 | 0.143 |
| P-oxidation / desulfuration | 36 | 1.000 | 0.667 | 0.667 | 0.917 | 0.889 |
| oxidative dehalogenation | 31 | 0.903 | 0.903 | 0.903 | 0.903 | 0.871 |
| aromatic C (substituted) | 22 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |

## Top-2 by isoform (scaffold split)

| isoform | n test inst. | atom-type frequency | fingerprint + GBM | Hückel physics + GBM | message-passing GNN | GNN + Hückel features |
|---|---|---|---|---|---|---|
| CYP1A2 | 59 | 0.663 | 0.730 | 0.796 | 0.787 | 0.769 |
| CYP2A6 | 24 | 0.658 | 0.673 | 0.796 | 0.766 | 0.768 |
| CYP2B6 | 34 | 0.617 | 0.735 | 0.801 | 0.770 | 0.724 |
| CYP2C19 | 46 | 0.601 | 0.777 | 0.841 | 0.771 | 0.744 |
| CYP2C8 | 29 | 0.536 | 0.701 | 0.819 | 0.760 | 0.762 |
| CYP2C9 | 48 | 0.625 | 0.775 | 0.842 | 0.764 | 0.758 |
| CYP2D6 | 60 | 0.590 | 0.757 | 0.848 | 0.845 | 0.826 |
| CYP2E1 | 33 | 0.686 | 0.735 | 0.798 | 0.754 | 0.780 |
| CYP3A4 | 98 | 0.554 | 0.734 | 0.780 | 0.793 | 0.762 |

## Top-1 by similarity to the nearest training molecule

| NN Tanimoto | n test inst. | atom-type frequency | fingerprint + GBM | Hückel physics + GBM | message-passing GNN | GNN + Hückel features |
|---|---|---|---|---|---|---|
| <=0.30 | 247 | 0.371 | 0.523 | 0.606 | 0.585 | 0.581 |
| 0.30-0.45 | 239 | 0.482 | 0.645 | 0.699 | 0.687 | 0.694 |
| 0.45-0.60 | 160 | 0.443 | 0.644 | 0.672 | 0.679 | 0.674 |
| 0.60-0.75 | 145 | 0.417 | 0.675 | 0.713 | 0.691 | 0.688 |
| >0.75 | 49 | 0.318 | 0.573 | 0.763 | 0.623 | 0.674 |

## Paired model comparisons

Paired t-tests on the per-split differences. Four metrics are tested per pair, so a single p just under 0.05 should be read with that multiplicity in mind.

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

