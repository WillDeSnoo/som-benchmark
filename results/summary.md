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
| message-passing GNN | 0.682 [0.661, 0.704] | 0.797 [0.774, 0.820] | 0.892 [0.880, 0.904] | 0.562 [0.535, 0.589] | 0.021 |
| GNN + Hückel features | 0.673 [0.648, 0.697] | 0.786 [0.756, 0.817] | 0.883 [0.871, 0.895] | 0.548 [0.510, 0.586] | 0.025 |
| GNN + BDE | 0.652 [0.623, 0.682] | 0.782 [0.759, 0.804] | 0.886 [0.872, 0.899] | 0.556 [0.519, 0.593] | 0.068 |

### scaffold split

| approach | top-1 | top-2 | ROC-AUC | avg precision | ECE |
|---|---|---|---|---|---|
| atom-type frequency | 0.425 [0.402, 0.448] | 0.607 [0.583, 0.631] | 0.802 [0.790, 0.814] | 0.274 [0.261, 0.287] | 0.025 |
| SMARTCyp (published method) | 0.635 [0.607, 0.663] | 0.770 [0.744, 0.795] | 0.826 [0.812, 0.840] | 0.341 [0.311, 0.370] | 0.018 |
| fingerprint + GBM | 0.601 [0.568, 0.635] | 0.741 [0.721, 0.760] | 0.871 [0.865, 0.878] | 0.499 [0.477, 0.522] | 0.018 |
| Hückel physics + GBM | 0.660 [0.626, 0.695] | 0.811 [0.787, 0.834] | 0.887 [0.876, 0.899] | 0.531 [0.499, 0.562] | 0.019 |
| Hückel + BDE + GBM | 0.674 [0.650, 0.698] | 0.812 [0.789, 0.834] | 0.907 [0.899, 0.915] | 0.584 [0.552, 0.616] | 0.019 |
| message-passing GNN | 0.643 [0.615, 0.670] | 0.788 [0.771, 0.804] | 0.880 [0.869, 0.891] | 0.504 [0.459, 0.549] | 0.026 |
| GNN + Hückel features | 0.640 [0.617, 0.663] | 0.770 [0.746, 0.794] | 0.868 [0.861, 0.875] | 0.487 [0.444, 0.530] | 0.028 |
| GNN + BDE | 0.623 [0.587, 0.658] | 0.778 [0.751, 0.805] | 0.882 [0.872, 0.891] | 0.521 [0.475, 0.568] | 0.026 |

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
| message-passing GNN | top1 | 0.682 | 0.643 | 0.040 [0.009, 0.070] |
| message-passing GNN | top2 | 0.797 | 0.788 | 0.009 [-0.016, 0.033] |
| GNN + Hückel features | top1 | 0.673 | 0.640 | 0.033 [0.004, 0.062] |
| GNN + Hückel features | top2 | 0.786 | 0.770 | 0.016 [-0.017, 0.050] |
| GNN + BDE | top1 | 0.652 | 0.623 | 0.030 [-0.010, 0.070] |
| GNN + BDE | top2 | 0.782 | 0.778 | 0.004 [-0.027, 0.034] |

## Recall at top-2 by site type (scaffold split)

| site type | n atoms | atom-type frequency | SMARTCyp (published method) | fingerprint + GBM | Hückel physics + GBM | Hückel + BDE + GBM | message-passing GNN | GNN + Hückel features | GNN + BDE |
|---|---|---|---|---|---|---|---|---|---|
| aromatic C-H hydroxylation | 2412 | 0.103 | 0.256 | 0.270 | 0.358 | 0.329 | 0.301 | 0.295 | 0.255 |
| N-dealkylation (alpha to N) | 2221 | 0.487 | 0.751 | 0.692 | 0.752 | 0.773 | 0.726 | 0.734 | 0.737 |
| O-dealkylation (alpha to O) | 1254 | 0.618 | 0.676 | 0.727 | 0.748 | 0.754 | 0.752 | 0.745 | 0.787 |
| benzylic / allylic hydroxylation | 1082 | 0.378 | 0.539 | 0.438 | 0.435 | 0.493 | 0.445 | 0.435 | 0.449 |
| S-oxidation | 435 | 0.910 | 0.662 | 0.749 | 0.703 | 0.736 | 0.814 | 0.651 | 0.724 |
| aliphatic C-H hydroxylation | 365 | 0.184 | 0.321 | 0.340 | 0.419 | 0.403 | 0.447 | 0.389 | 0.400 |
| N-oxidation | 319 | 0.210 | 0.407 | 0.389 | 0.342 | 0.345 | 0.326 | 0.291 | 0.376 |
| methyl hydroxylation | 290 | 0.665 | 0.083 | 0.241 | 0.190 | 0.217 | 0.231 | 0.155 | 0.190 |
| O-site (direct) | 91 | 0.000 | 0.000 | 0.077 | 0.143 | 0.165 | 0.055 | 0.132 | 0.099 |
| P-oxidation / desulfuration | 36 | 1.000 | 0.000 | 0.667 | 0.667 | 0.667 | 0.972 | 0.917 | 0.917 |
| oxidative dehalogenation | 31 | 0.903 | 0.290 | 0.903 | 0.903 | 0.903 | 0.903 | 0.903 | 0.903 |
| aromatic C (substituted) | 22 | 0.000 | 0.000 | 0.000 | 0.000 | 0.045 | 0.000 | 0.000 | 0.000 |

## Top-2 by isoform (scaffold split)

| isoform | n test inst. | atom-type frequency | SMARTCyp (published method) | fingerprint + GBM | Hückel physics + GBM | Hückel + BDE + GBM | message-passing GNN | GNN + Hückel features | GNN + BDE |
|---|---|---|---|---|---|---|---|---|---|
| CYP1A2 | 59 | 0.663 | 0.751 | 0.730 | 0.796 | 0.806 | 0.758 | 0.777 | 0.774 |
| CYP2A6 | 24 | 0.658 | 0.811 | 0.673 | 0.796 | 0.818 | 0.775 | 0.741 | 0.771 |
| CYP2B6 | 34 | 0.617 | 0.734 | 0.735 | 0.801 | 0.772 | 0.761 | 0.724 | 0.730 |
| CYP2C19 | 46 | 0.601 | 0.794 | 0.777 | 0.841 | 0.829 | 0.776 | 0.791 | 0.797 |
| CYP2C8 | 29 | 0.536 | 0.790 | 0.701 | 0.819 | 0.831 | 0.817 | 0.739 | 0.777 |
| CYP2C9 | 48 | 0.625 | 0.808 | 0.775 | 0.842 | 0.845 | 0.786 | 0.764 | 0.786 |
| CYP2D6 | 60 | 0.590 | 0.806 | 0.757 | 0.848 | 0.842 | 0.825 | 0.833 | 0.810 |
| CYP2E1 | 33 | 0.686 | 0.775 | 0.735 | 0.798 | 0.775 | 0.767 | 0.779 | 0.742 |
| CYP3A4 | 98 | 0.554 | 0.726 | 0.734 | 0.780 | 0.790 | 0.799 | 0.750 | 0.776 |

## Top-1 by similarity to the nearest training molecule

| NN Tanimoto | n test inst. | atom-type frequency | SMARTCyp (published method) | fingerprint + GBM | Hückel physics + GBM | Hückel + BDE + GBM | message-passing GNN | GNN + Hückel features | GNN + BDE |
|---|---|---|---|---|---|---|---|---|---|
| <=0.30 | 247 | 0.371 | 0.585 | 0.523 | 0.606 | 0.626 | 0.576 | 0.593 | 0.572 |
| 0.30-0.45 | 239 | 0.482 | 0.666 | 0.645 | 0.699 | 0.729 | 0.695 | 0.698 | 0.694 |
| 0.45-0.60 | 160 | 0.443 | 0.663 | 0.644 | 0.672 | 0.680 | 0.704 | 0.675 | 0.635 |
| 0.60-0.75 | 145 | 0.417 | 0.637 | 0.675 | 0.713 | 0.734 | 0.712 | 0.679 | 0.665 |
| >0.75 | 49 | 0.318 | 0.546 | 0.573 | 0.763 | 0.719 | 0.627 | 0.630 | 0.599 |

## Paired model comparisons

Paired t-tests on the per-split differences. The pooled column averages the random and scaffold splits together; the scaffold column is the harder split on its own, and is the one to quote when the claim is about generalisation to new scaffolds. Four metrics are tested per pair, so a single p just under 0.05 should be read with that multiplicity in mind.

**phys_bde_gbm − smartcyp** — best model vs the published external method (pooled over 20 splits; scaffold column is the 10 scaffold splits alone)

| metric | pooled difference [95% CI] | p | scaffold only [95% CI] | p |
|---|---|---|---|---|
| top1 | +0.0578 [+0.0361, +0.0796] | 0.0000 | +0.0389 [+0.0085, +0.0692] | 0.0177 |
| top2 | +0.0549 [+0.0398, +0.0699] | 0.0000 | +0.0421 [+0.0245, +0.0596] | 0.0004 |
| roc_auc | +0.0875 [+0.0798, +0.0952] | 0.0000 | +0.0810 [+0.0682, +0.0938] | 0.0000 |
| avg_precision | +0.2671 [+0.2468, +0.2875] | 0.0000 | +0.2431 [+0.2149, +0.2712] | 0.0000 |

**phys_bde_gbm − gnn** — best model vs the tuned GNN (pooled over 20 splits; scaffold column is the 10 scaffold splits alone)

| metric | pooled difference [95% CI] | p | scaffold only [95% CI] | p |
|---|---|---|---|---|
| top1 | +0.0274 [+0.0123, +0.0425] | 0.0012 | +0.0313 [+0.0052, +0.0574] | 0.0240 |
| top2 | +0.0273 [+0.0165, +0.0381] | 0.0000 | +0.0236 [+0.0074, +0.0398] | 0.0093 |
| roc_auc | +0.0255 [+0.0197, +0.0313] | 0.0000 | +0.0268 [+0.0155, +0.0382] | 0.0005 |
| avg_precision | +0.0744 [+0.0587, +0.0901] | 0.0000 | +0.0800 [+0.0500, +0.1099] | 0.0002 |

**phys_bde_gbm − gnn_huckel** — best model vs the GNN given the Hückel descriptors (pooled over 20 splits; scaffold column is the 10 scaffold splits alone)

| metric | pooled difference [95% CI] | p | scaffold only [95% CI] | p |
|---|---|---|---|---|
| top1 | +0.0336 [+0.0197, +0.0475] | 0.0001 | +0.0339 [+0.0123, +0.0555] | 0.0062 |
| top2 | +0.0414 [+0.0265, +0.0564] | 0.0000 | +0.0415 [+0.0211, +0.0619] | 0.0013 |
| roc_auc | +0.0358 [+0.0312, +0.0405] | 0.0000 | +0.0388 [+0.0318, +0.0458] | 0.0000 |
| avg_precision | +0.0902 [+0.0773, +0.1031] | 0.0000 | +0.0969 [+0.0734, +0.1204] | 0.0000 |

**phys_bde_gbm − gnn_bde** — best model vs the GNN given the bond dissociation energies (pooled over 20 splits; scaffold column is the 10 scaffold splits alone)

| metric | pooled difference [95% CI] | p | scaffold only [95% CI] | p |
|---|---|---|---|---|
| top1 | +0.0523 [+0.0344, +0.0703] | 0.0000 | +0.0512 [+0.0174, +0.0850] | 0.0075 |
| top2 | +0.0399 [+0.0263, +0.0534] | 0.0000 | +0.0336 [+0.0097, +0.0576] | 0.0113 |
| roc_auc | +0.0279 [+0.0227, +0.0331] | 0.0000 | +0.0254 [+0.0176, +0.0332] | 0.0000 |
| avg_precision | +0.0688 [+0.0536, +0.0839] | 0.0000 | +0.0625 [+0.0385, +0.0865] | 0.0002 |

**smartcyp − prior** — does the published method beat the trivial reference? (pooled over 20 splits; scaffold column is the 10 scaffold splits alone)

| metric | pooled difference [95% CI] | p | scaffold only [95% CI] | p |
|---|---|---|---|---|
| top1 | +0.2100 [+0.1852, +0.2348] | 0.0000 | +0.2102 [+0.1654, +0.2549] | 0.0000 |
| top2 | +0.1628 [+0.1418, +0.1837] | 0.0000 | +0.1627 [+0.1236, +0.2018] | 0.0000 |
| roc_auc | +0.0236 [+0.0152, +0.0319] | 0.0000 | +0.0242 [+0.0078, +0.0406] | 0.0087 |
| avg_precision | +0.0646 [+0.0465, +0.0827] | 0.0000 | +0.0667 [+0.0316, +0.1018] | 0.0020 |

**phys_bde_gbm − phys_gbm** — does adding transferred DFT bond dissociation energies help? (pooled over 20 splits; scaffold column is the 10 scaffold splits alone)

| metric | pooled difference [95% CI] | p | scaffold only [95% CI] | p |
|---|---|---|---|---|
| top1 | +0.0179 [+0.0061, +0.0297] | 0.0050 | +0.0137 [-0.0065, +0.0339] | 0.1590 |
| top2 | +0.0036 [-0.0051, +0.0123] | 0.4027 | +0.0010 [-0.0142, +0.0163] | 0.8814 |
| roc_auc | +0.0174 [+0.0125, +0.0222] | 0.0000 | +0.0198 [+0.0112, +0.0285] | 0.0006 |
| avg_precision | +0.0475 [+0.0375, +0.0575] | 0.0000 | +0.0529 [+0.0344, +0.0713] | 0.0001 |

**gnn_bde − gnn** — same question, for the GNN (pooled over 20 splits; scaffold column is the 10 scaffold splits alone)

| metric | pooled difference [95% CI] | p | scaffold only [95% CI] | p |
|---|---|---|---|---|
| top1 | -0.0249 [-0.0387, -0.0111] | 0.0013 | -0.0200 [-0.0425, +0.0026] | 0.0761 |
| top2 | -0.0126 [-0.0261, +0.0009] | 0.0665 | -0.0100 [-0.0343, +0.0143] | 0.3756 |
| roc_auc | -0.0024 [-0.0085, +0.0037] | 0.4149 | +0.0014 [-0.0096, +0.0124] | 0.7754 |
| avg_precision | +0.0056 [-0.0112, +0.0225] | 0.4909 | +0.0175 [-0.0005, +0.0354] | 0.0554 |

**gnn_huckel − gnn** — does giving the GNN the Hückel descriptors help? (pooled over 20 splits; scaffold column is the 10 scaffold splits alone)

| metric | pooled difference [95% CI] | p | scaffold only [95% CI] | p |
|---|---|---|---|---|
| top1 | -0.0062 [-0.0191, +0.0066] | 0.3244 | -0.0027 [-0.0218, +0.0165] | 0.7618 |
| top2 | -0.0141 [-0.0268, -0.0015] | 0.0299 | -0.0179 [-0.0380, +0.0022] | 0.0745 |
| roc_auc | -0.0103 [-0.0158, -0.0048] | 0.0009 | -0.0119 [-0.0229, -0.0009] | 0.0364 |
| avg_precision | -0.0158 [-0.0282, -0.0033] | 0.0155 | -0.0169 [-0.0380, +0.0041] | 0.1026 |

**phys_gbm − gnn_huckel** — same features, different model class (pooled over 20 splits; scaffold column is the 10 scaffold splits alone)

| metric | pooled difference [95% CI] | p | scaffold only [95% CI] | p |
|---|---|---|---|---|
| top1 | +0.0157 [-0.0005, +0.0319] | 0.0567 | +0.0202 [-0.0108, +0.0512] | 0.1746 |
| top2 | +0.0379 [+0.0247, +0.0510] | 0.0000 | +0.0405 [+0.0191, +0.0619] | 0.0021 |
| roc_auc | +0.0185 [+0.0138, +0.0231] | 0.0000 | +0.0189 [+0.0126, +0.0253] | 0.0001 |
| avg_precision | +0.0427 [+0.0317, +0.0536] | 0.0000 | +0.0440 [+0.0226, +0.0655] | 0.0012 |

**phys_gbm − fp_gbm** — physics descriptors vs hashed environments (pooled over 20 splits; scaffold column is the 10 scaffold splits alone)

| metric | pooled difference [95% CI] | p | scaffold only [95% CI] | p |
|---|---|---|---|---|
| top1 | +0.0606 [+0.0448, +0.0764] | 0.0000 | +0.0592 [+0.0388, +0.0796] | 0.0001 |
| top2 | +0.0732 [+0.0597, +0.0866] | 0.0000 | +0.0699 [+0.0466, +0.0932] | 0.0001 |
| roc_auc | +0.0192 [+0.0137, +0.0248] | 0.0000 | +0.0159 [+0.0074, +0.0243] | 0.0021 |
| avg_precision | +0.0489 [+0.0334, +0.0644] | 0.0000 | +0.0314 [+0.0175, +0.0452] | 0.0006 |

**phys_gbm − prior** — best model vs the frequency reference (pooled over 20 splits; scaffold column is the 10 scaffold splits alone)

| metric | pooled difference [95% CI] | p | scaffold only [95% CI] | p |
|---|---|---|---|---|
| top1 | +0.2500 [+0.2262, +0.2738] | 0.0000 | +0.2353 [+0.1979, +0.2728] | 0.0000 |
| top2 | +0.2141 [+0.1939, +0.2342] | 0.0000 | +0.2037 [+0.1677, +0.2397] | 0.0000 |
| roc_auc | +0.0937 [+0.0864, +0.1009] | 0.0000 | +0.0854 [+0.0761, +0.0946] | 0.0000 |
| avg_precision | +0.2842 [+0.2622, +0.3063] | 0.0000 | +0.2569 [+0.2288, +0.2849] | 0.0000 |

## The GNN hyperparameter search

Selected on validation average precision. The test set plays no part in any row below. Trial 0 is the configuration used before tuning, so the search can only improve on it.

**message-passing GNN, scaffold split** (full search, 12 configurations)

| configuration | epochs | validation AP |
|---|---|---|
| lr=0.0003 | 58 | 0.5663 **<-** |
| pos_weight=3.0 | 38 | 0.5649 |
| dropout=0.1 | 54 | 0.5644 |
| layers=2 | 36 | 0.5622 |
| base (pre-tuning) | 45 | 0.5595 |
| dropout=0.25 | 62 | 0.5591 |
| hidden=128, layers=3, dropout=0.1, pos_weight=3.0 | 39 | 0.5582 |
| weight_decay=0.001 | 45 | 0.5538 |
| layers=3 | 40 | 0.5494 |
| hidden=64 | 46 | 0.5453 |
| lr=0.003 | 35 | 0.5435 |
| hidden=64, layers=3, dropout=0.1, weight_decay=0.001 | 54 | 0.5416 |

**message-passing GNN, random split** (shortlist search, 3 configurations)

| configuration | epochs | validation AP |
|---|---|---|
| lr=0.0003 | 85 | 0.6241 **<-** |
| pos_weight=3.0 | 47 | 0.6192 |
| dropout=0.1 | 58 | 0.6097 |

**GNN + Hückel features, scaffold split** (full search, 12 configurations)

| configuration | epochs | validation AP |
|---|---|---|
| layers=2 | 34 | 0.5513 **<-** |
| pos_weight=3.0 | 31 | 0.5495 |
| hidden=128, layers=3, dropout=0.1, pos_weight=3.0 | 31 | 0.5476 |
| lr=0.0003 | 38 | 0.5475 |
| hidden=64, layers=3, dropout=0.1, weight_decay=0.001 | 29 | 0.5462 |
| layers=3 | 32 | 0.5459 |
| hidden=64 | 33 | 0.5449 |
| weight_decay=0.001 | 35 | 0.5427 |
| dropout=0.25 | 35 | 0.5421 |
| base (pre-tuning) | 42 | 0.5420 |
| dropout=0.1 | 31 | 0.5409 |
| lr=0.003 | 33 | 0.5381 |

**GNN + Hückel features, random split** (shortlist search, 3 configurations)

| configuration | epochs | validation AP |
|---|---|---|
| layers=2 | 63 | 0.6079 **<-** |
| hidden=128, layers=3, dropout=0.1, pos_weight=3.0 | 60 | 0.6061 |
| pos_weight=3.0 | 48 | 0.5975 |

**GNN + BDE, scaffold split** (full search, 12 configurations)

| configuration | epochs | validation AP |
|---|---|---|
| hidden=64 | 42 | 0.6093 **<-** |
| layers=2 | 41 | 0.6077 |
| hidden=128, layers=3, dropout=0.1, pos_weight=3.0 | 41 | 0.6039 |
| lr=0.0003 | 41 | 0.6007 |
| dropout=0.25 | 38 | 0.5996 |
| hidden=64, layers=3, dropout=0.1, weight_decay=0.001 | 39 | 0.5995 |
| lr=0.003 | 34 | 0.5990 |
| pos_weight=3.0 | 33 | 0.5983 |
| dropout=0.1 | 34 | 0.5978 |
| layers=3 | 38 | 0.5963 |
| weight_decay=0.001 | 33 | 0.5955 |
| base (pre-tuning) | 33 | 0.5913 |

**GNN + BDE, random split** (shortlist search, 3 configurations)

| configuration | epochs | validation AP |
|---|---|---|
| hidden=128, layers=3, dropout=0.1, pos_weight=3.0 | 58 | 0.6113 **<-** |
| hidden=64 | 80 | 0.6038 |
| layers=2 | 49 | 0.5898 |

## Split diagnostics (seed 0)

| mode | test molecules | shared scaffolds | median NN Tanimoto | frac. test with NN > 0.6 |
|---|---|---|---|---|
| scaffold | 136 | 0 | 0.343 | 0.103 |
| random | 136 | 34 | 0.491 | 0.368 |

