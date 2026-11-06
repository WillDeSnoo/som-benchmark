#!/usr/bin/env python
"""Stage 3: train and evaluate every model on every split, repeated over seeds.

For each (model, split mode, seed) this writes the raw per-atom test
predictions to data/interim/predictions/ so that stage 4 can do the failure
analysis without retraining anything, and appends a row of metrics to
results/metrics_per_run.csv. Aggregation across seeds goes to
results/metrics_summary.csv.

The test set is used exactly once per run, to produce those predictions.
Hyperparameters come from the validation split of the first seed only.
"""

import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np

from som.data import load_instances, load_molecules
from som.metrics import TOP_KS, bootstrap_ci, evaluate, mean_ci_across_seeds
from som.models import MODELS, clear_caches, fit_predict
from som.paths import INTERIM, RESULTS, ensure_dirs
from som.splits import SPLIT_MODES, make_split, split_diagnostics

METRIC_KEYS = (
    [f"top{k}" for k in TOP_KS]
    + [f"top{k}_pessimistic" for k in TOP_KS]
    + ["roc_auc", "avg_precision", "brier", "ece"]
)
CI_METRICS = ("top1", "top2", "roc_auc", "avg_precision")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--models", nargs="+", default=list(MODELS), choices=MODELS)
    ap.add_argument("--modes", nargs="+", default=list(SPLIT_MODES),
                    choices=SPLIT_MODES)
    ap.add_argument("--label-field", default="som_all",
                    choices=["som_all", "som_primary"])
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--test-frac", type=float, default=0.2)
    ap.add_argument("--val-frac", type=float, default=0.1)
    args = ap.parse_args()

    ensure_dirs()
    pred_dir = INTERIM / "predictions"
    pred_dir.mkdir(parents=True, exist_ok=True)

    instances = load_instances()
    molecules = load_molecules()
    mol_of_inst = {i: r["smiles"] for i, r in enumerate(instances)}
    print(f"{len(instances)} instances over {len(molecules)} molecules, "
          f"label field = {args.label_field}\n")

    clear_caches()
    rows = []
    diag_rows = []

    for mode in args.modes:
        for seed in range(args.seeds):
            split = make_split(
                instances, molecules, mode, seed,
                test_frac=args.test_frac, val_frac=args.val_frac,
            )
            if seed == 0:
                diag = split_diagnostics(split, molecules)
                diag_rows.append(diag)
                print(f"[{mode} seed {seed}] "
                      f"train/val/test molecules = {len(split.train_mols)}/"
                      f"{len(split.val_mols)}/{len(split.test_mols)}, "
                      f"shared scaffolds = {diag['shared_scaffolds']}, "
                      f"median NN Tanimoto = {diag['nn_tanimoto_median']:.3f}")

            for name in args.models:
                t0 = time.time()
                pred = fit_predict(
                    name, instances, split,
                    label_field=args.label_field, seed=seed, tune=True,
                )
                raw = evaluate(pred.y_test, pred.p_test, pred.inst_test, seed=seed)
                cal = (
                    evaluate(pred.y_test, pred.p_test_calibrated,
                             pred.inst_test, seed=seed)
                    if pred.p_test_calibrated is not None else {}
                )
                ci = bootstrap_ci(
                    pred.y_test, pred.p_test, pred.inst_test, mol_of_inst,
                    metrics=CI_METRICS, n_boot=args.n_boot, seed=seed,
                )
                dt = time.time() - t0

                np.savez_compressed(
                    pred_dir / f"{name}_{mode}_seed{seed}.npz",
                    y=pred.y_test, p=pred.p_test,
                    p_cal=(pred.p_test_calibrated
                           if pred.p_test_calibrated is not None
                           else np.full_like(pred.p_test, np.nan)),
                    inst=pred.inst_test, atom=pred.atom_test,
                )

                row = {
                    "model": name, "mode": mode, "seed": seed,
                    "seconds": round(dt, 1),
                    "n_test_instances": raw["n_instances"],
                    "n_test_atoms": raw["n_atoms"],
                    "positive_rate": round(raw["positive_rate"], 5),
                }
                for k in METRIC_KEYS:
                    row[k] = raw.get(k, float("nan"))
                row["ece_calibrated"] = cal.get("ece", float("nan"))
                row["brier_calibrated"] = cal.get("brier", float("nan"))
                for m, (lo, hi) in ci.items():
                    row[f"{m}_boot_lo"] = lo
                    row[f"{m}_boot_hi"] = hi
                row["params"] = json.dumps(pred.chosen_params)
                rows.append(row)

                print(f"  {name:9} {mode:8} seed {seed:2}  "
                      f"top1 {raw['top1']:.3f}  top2 {raw['top2']:.3f}  "
                      f"auc {raw['roc_auc']:.3f}  ap {raw['avg_precision']:.3f}  "
                      f"ece {raw['ece']:.3f}  ({dt:.0f}s)")

                if name == "phys_gbm" and pred.feature_importance and seed == 0:
                    out = RESULTS / f"phys_linear_coefficients_{mode}.csv"
                    with out.open("w", newline="") as fh:
                        w = csv.writer(fh)
                        w.writerow(["feature", "coefficient"])
                        for f, c in sorted(
                            pred.feature_importance.items(),
                            key=lambda kv: -abs(kv[1]),
                        ):
                            w.writerow([f, f"{c:.6f}"])

    # ---------------- write per-run table ----------------
    # Merge with whatever is already there, keyed on (model, mode, seed), so
    # that a later run covering only one model adds to the table instead of
    # replacing it. Rows recomputed in this run win.
    per_run = RESULTS / "metrics_per_run.csv"
    fields = list(rows[0].keys())
    merged: dict[tuple[str, str, str], dict] = {}
    numeric = set(METRIC_KEYS) | {"ece_calibrated", "brier_calibrated"}
    if per_run.exists():
        with per_run.open(newline="") as fh:
            for old in csv.DictReader(fh):
                # values come back as strings; the aggregation needs floats
                for k in numeric:
                    if k in old:
                        try:
                            old[k] = float(old[k])
                        except (TypeError, ValueError):
                            old[k] = float("nan")
                merged[(old["model"], old["mode"], str(old["seed"]))] = old
    for r in rows:
        merged[(r["model"], r["mode"], str(r["seed"]))] = r

    def sort_key(item):
        (model, mode, seed), _ = item
        order = list(MODELS)
        return (mode, order.index(model) if model in order else 99, int(seed))

    with per_run.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for _, r in sorted(merged.items(), key=sort_key):
            w.writerow(r)
    print(f"\nwrote {per_run} ({len(merged)} rows, {len(rows)} from this run)")
    rows = [dict(r) for _, r in sorted(merged.items(), key=sort_key)]

    diag_path = RESULTS / "split_diagnostics.csv"
    if diag_rows:
        with diag_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(diag_rows[0].keys()))
            w.writeheader()
            for r in diag_rows:
                w.writerow(r)
        print(f"wrote {diag_path}")

    # ---------------- aggregate across seeds ----------------
    summary = []
    all_modes = [m for m in SPLIT_MODES if any(r["mode"] == m for r in rows)]
    all_models = [m for m in MODELS if any(r["model"] == m for r in rows)]
    for mode in all_modes:
        for name in all_models:
            sub = [r for r in rows if r["model"] == name and r["mode"] == mode]
            if not sub:
                continue
            agg = {"model": name, "mode": mode, "n_seeds": len(sub)}
            for k in METRIC_KEYS + ["ece_calibrated", "brier_calibrated"]:
                s = mean_ci_across_seeds([r[k] for r in sub])
                agg[f"{k}_mean"] = round(s["mean"], 4)
                agg[f"{k}_lo"] = round(s["lo"], 4)
                agg[f"{k}_hi"] = round(s["hi"], 4)
            summary.append(agg)

    summ_path = RESULTS / "metrics_summary.csv"
    with summ_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
        w.writeheader()
        for r in summary:
            w.writerow(r)
    print(f"wrote {summ_path}")

    print("\n" + "=" * 78)
    print("across-seed means with 95% CI of the mean")
    print("=" * 78)
    print(f"{'mode':<10}{'model':<11}{'top-1':>18}{'top-2':>18}{'ROC-AUC':>16}")
    for r in summary:
        def f(k):
            return f"{r[k + '_mean']:.3f} [{r[k + '_lo']:.3f},{r[k + '_hi']:.3f}]"
        print(f"{r['mode']:<10}{r['model']:<11}{f('top1'):>18}"
              f"{f('top2'):>18}{f('roc_auc'):>16}")
    print("=" * 78)


if __name__ == "__main__":
    main()
