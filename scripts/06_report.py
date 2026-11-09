#!/usr/bin/env python
"""Stage 6: render the result CSVs as markdown tables into results/summary.md.

The README quotes numbers from this file, so they are generated rather than
transcribed by hand.
"""

import argparse
import csv
import json

from som.paths import RESULTS
from som.viz import MODEL_LABELS

ORDER = ["prior", "smartcyp", "fp_gbm", "basic_gbm", "phys_gbm",
         "phys_bde_gbm", "gnn", "gnn_huckel", "gnn_bde"]


def read(name):
    path = RESULTS / name
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def fv(row, key, nd=3):
    try:
        return f"{float(row[key]):.{nd}f}"
    except (KeyError, ValueError, TypeError):
        return "-"


def ci(row, key, nd=3):
    try:
        return (f"{float(row[key + '_mean']):.{nd}f} "
                f"[{float(row[key + '_lo']):.{nd}f}, "
                f"{float(row[key + '_hi']):.{nd}f}]")
    except (KeyError, ValueError, TypeError):
        return "-"


def sort_models(rows):
    idx = {m: i for i, m in enumerate(ORDER)}
    return sorted(rows, key=lambda r: idx.get(r["model"], 99))


def main() -> None:
    argparse.ArgumentParser(description=__doc__).parse_args()
    out = ["# Generated results", ""]

    cur_path = RESULTS / "curation_summary.json"
    if cur_path.exists():
        c = json.loads(cur_path.read_text())
        out += [
            "## Curation", "",
            f"- records read from the upstream SDFs: **{c['records_read']}**",
            f"- excluded: **{c['excluded_total']}**"
            + (f" ({c['excluded_by_reason']})" if c["excluded_by_reason"] else ""),
            f"- duplicate (molecule, isoform) instances merged: "
            f"**{c['duplicate_instances_merged']}** "
            f"(label disagreements: {c['duplicate_label_conflicts']})",
            f"- surviving instances: **{c['instances']}** over "
            f"**{c['molecules']}** unique molecules",
            f"- Bemis-Murcko scaffolds: **{c['scaffolds']}** "
            f"(+{c['acyclic_molecules']} acyclic molecules)",
            f"- candidate atoms: **{c['candidate_atoms']}**, of which "
            f"**{c['positive_atoms_all_tiers']}** are labelled sites "
            f"({100 * c['positive_atoms_all_tiers'] / c['candidate_atoms']:.2f}%)",
            f"- primary-only positives: **{c['positive_atoms_primary_only']}** "
            f"({100 * c['positive_atoms_primary_only'] / c['candidate_atoms']:.2f}%)",
            "",
        ]

    summ = read("metrics_summary.csv")
    if summ:
        out += ["## Headline results", "",
                "Mean over repeated splits, with the 95% CI of the mean.", ""]
        for mode in ("random", "scaffold"):
            sub = sort_models([r for r in summ if r["mode"] == mode])
            if not sub:
                continue
            out += [
                f"### {mode} split", "",
                "| approach | top-1 | top-2 | ROC-AUC | avg precision | ECE |",
                "|---|---|---|---|---|---|",
            ]
            for r in sub:
                out.append(
                    f"| {MODEL_LABELS.get(r['model'], r['model'])} "
                    f"| {ci(r, 'top1')} | {ci(r, 'top2')} "
                    f"| {ci(r, 'roc_auc')} | {ci(r, 'avg_precision')} "
                    f"| {fv(r, 'ece_mean')} |"
                )
            out.append("")

    gap = read("split_gap.csv")
    if gap:
        out += ["## Random minus scaffold", "",
                "| approach | metric | random | scaffold | gap [95% CI] |",
                "|---|---|---|---|---|"]
        for r in sort_models(gap):
            out.append(
                f"| {MODEL_LABELS.get(r['model'], r['model'])} | {r['metric']} "
                f"| {fv(r, 'random_mean')} | {fv(r, 'scaffold_mean')} "
                f"| {fv(r, 'gap')} [{fv(r, 'gap_lo')}, {fv(r, 'gap_hi')}] |"
            )
        out.append("")

    site = read("per_site_type.csv")
    if site:
        models = [m for m in ORDER if m in {r["model"] for r in site}]
        by = {(r["model"], r["site_type"]): r
              for r in site if r["mode"] == "scaffold"}
        types = []
        counts = {}
        for r in site:
            if r["mode"] != "scaffold":
                continue
            t = r["site_type"]
            counts[t] = max(counts.get(t, 0),
                            int(float(r["n_positive_atoms_pooled"])))
            if t not in types:
                types.append(t)
        types.sort(key=lambda t: -counts[t])
        out += ["## Recall at top-2 by site type (scaffold split)", "",
                "| site type | n atoms | "
                + " | ".join(MODEL_LABELS.get(m, m) for m in models) + " |",
                "|---" * (2 + len(models)) + "|"]
        for t in types:
            cells = [fv(by.get((m, t), {}), "recall_at_2_pooled")
                     for m in models]
            out.append(f"| {t} | {counts[t]} | " + " | ".join(cells) + " |")
        out.append("")

    iso = read("per_isoform.csv")
    if iso:
        models = [m for m in ORDER if m in {r["model"] for r in iso}]
        sub = [r for r in iso if r["mode"] == "scaffold"]
        isoforms = sorted({r["isoform"] for r in sub})
        by = {(r["model"], r["isoform"]): r for r in sub}
        out += ["## Top-2 by isoform (scaffold split)", "",
                "| isoform | n test inst. | "
                + " | ".join(MODEL_LABELS.get(m, m) for m in models) + " |",
                "|---" * (2 + len(models)) + "|"]
        for i in isoforms:
            n = fv(by.get((models[0], i), {}), "mean_test_instances", 0)
            cells = [fv(by.get((m, i), {}), "top2_mean") for m in models]
            out.append(f"| CYP{i} | {n} | " + " | ".join(cells) + " |")
        out.append("")

    nov = read("novelty_curve.csv")
    if nov:
        models = [m for m in ORDER if m in {r["model"] for r in nov}]
        sub = [r for r in nov if r["mode"] == "pooled"]
        order = ["<=0.30", "0.30-0.45", "0.45-0.60", "0.60-0.75", ">0.75"]
        buckets = [b for b in order
                   if b in {r["nn_tanimoto_bucket"] for r in sub}]
        by = {(r["model"], r["nn_tanimoto_bucket"]): r for r in sub}
        out += ["## Top-1 by similarity to the nearest training molecule", "",
                "| NN Tanimoto | n test inst. | "
                + " | ".join(MODEL_LABELS.get(m, m) for m in models) + " |",
                "|---" * (2 + len(models)) + "|"]
        for b in buckets:
            n = fv(by.get((models[0], b), {}), "mean_test_instances", 0)
            cells = [fv(by.get((m, b), {}), "top1_mean") for m in models]
            out.append(f"| {b} | {n} | " + " | ".join(cells) + " |")
        out.append("")

    # ---- paired comparisons --------------------------------------------
    # Models are run on the same (mode, seed) splits, so differences can be
    # paired, which is far more sensitive than comparing two overlapping CIs.
    per_run = read("metrics_per_run.csv")
    if per_run:
        import numpy as np
        from scipy import stats

        by_split: dict[tuple[str, str], dict[str, dict]] = {}
        for r in per_run:
            by_split.setdefault((r["mode"], r["seed"]), {})[r["model"]] = r

        pairs = [
            ("phys_bde_gbm", "smartcyp",
             "best model vs the published external method"),
            ("phys_bde_gbm", "gnn",
             "best model vs the tuned GNN"),
            ("phys_bde_gbm", "gnn_huckel",
             "best model vs the GNN given the Hückel descriptors"),
            ("phys_bde_gbm", "gnn_bde",
             "best model vs the GNN given the bond dissociation energies"),
            ("smartcyp", "prior",
             "does the published method beat the trivial reference?"),
            ("phys_bde_gbm", "phys_gbm",
             "does adding transferred DFT bond dissociation energies help?"),
            ("gnn_bde", "gnn",
             "same question, for the GNN"),
            ("gnn_huckel", "gnn",
             "does giving the GNN the Hückel descriptors help?"),
            ("phys_gbm", "gnn_huckel",
             "same features, different model class"),
            ("phys_gbm", "basic_gbm",
             "what the 21 Hückel columns contribute, estimator held fixed"),
            ("basic_gbm", "gnn",
             "identical features, different estimator: the cleanest "
             "model-class comparison here"),
            ("basic_gbm", "fp_gbm",
             "hand-built atom/bond features vs hashed environments, "
             "no Hückel theory involved"),
            ("phys_gbm", "fp_gbm",
             "physics descriptors vs hashed environments"),
            ("phys_gbm", "prior",
             "best model vs the frequency reference"),
        ]
        def paired(pairs_a, pairs_b, metric, rows):
            """Paired t-test of model a against b over the given splits."""
            xa = np.array([float(v[pairs_a][metric]) for v in rows])
            xb = np.array([float(v[pairs_b][metric]) for v in rows])
            d = xa - xb
            _, p = stats.ttest_rel(xa, xb)
            half = (stats.t.ppf(0.975, d.size - 1)
                    * d.std(ddof=1) / np.sqrt(d.size))
            return (f"{d.mean():+.4f} "
                    f"[{d.mean() - half:+.4f}, {d.mean() + half:+.4f}]",
                    f"{p:.4f}")

        lines = []
        for a, b, why in pairs:
            avail = [v for k, v in by_split.items() if a in v and b in v]
            scaf = [v for k, v in by_split.items()
                    if a in v and b in v and k[0] == "scaffold"]
            if len(avail) < 3:
                continue
            lines.append(f"**{a} − {b}** — {why} (pooled over "
                         f"{len(avail)} splits; scaffold column is the "
                         f"{len(scaf)} scaffold splits alone)")
            lines.append("")
            lines.append("| metric | pooled difference [95% CI] | p "
                         "| scaffold only [95% CI] | p |")
            lines.append("|---|---|---|---|---|")
            for metric in ("top1", "top2", "roc_auc", "avg_precision"):
                pd, pp = paired(a, b, metric, avail)
                if len(scaf) >= 3:
                    sd, sp = paired(a, b, metric, scaf)
                else:
                    sd = sp = "-"
                lines.append(f"| {metric} | {pd} | {pp} | {sd} | {sp} |")
            lines.append("")
        if lines:
            out += ["## Paired model comparisons", "",
                    "Paired t-tests on the per-split differences. The pooled "
                    "column averages the random and scaffold splits together; "
                    "the scaffold column is the harder split on its own, and "
                    "is the one to quote when the claim is about "
                    "generalisation to new scaffolds. Four metrics are tested "
                    "per pair, so a single p just under 0.05 should be read "
                    "with that multiplicity in mind.", ""] + lines

    tun = read("gnn_tuning.csv")
    if tun:
        out += ["## The GNN hyperparameter search", "",
                "Selected on validation average precision. The test set plays "
                "no part in any row below. Trial 0 is the configuration used "
                "before tuning, so the search can only improve on it.", ""]
        for model in [m for m in ORDER if m in {r["model"] for r in tun}]:
            for mode in ("scaffold", "random"):
                sub = [r for r in tun
                       if r["model"] == model and r["mode"] == mode]
                if not sub:
                    continue
                best = max(sub, key=lambda r: float(r["val_avg_precision"]))
                out += [f"**{MODEL_LABELS.get(model, model)}, {mode} split** "
                        f"({sub[0]['scope']} search, {len(sub)} "
                        f"configurations)", "",
                        "| configuration | epochs | validation AP |",
                        "|---|---|---|"]
                for r in sorted(sub,
                                key=lambda r: -float(r["val_avg_precision"])):
                    params = json.loads(r["params"])
                    label = ", ".join(f"{k}={v}" for k, v in params.items()) \
                        or "base (pre-tuning)"
                    mark = " **<-**" if r is best else ""
                    out.append(f"| {label} | {r['epochs_run']} "
                               f"| {fv(r, 'val_avg_precision', 4)}{mark} |")
                out.append("")

    diag = read("split_diagnostics.csv")
    if diag:
        out += ["## Split diagnostics (seed 0)", "",
                "| mode | test molecules | shared scaffolds | "
                "median NN Tanimoto | frac. test with NN > 0.6 |",
                "|---|---|---|---|---|"]
        for r in diag:
            out.append(
                f"| {r['mode']} | {r['test_molecules']} "
                f"| {r['shared_scaffolds']} "
                f"| {fv(r, 'nn_tanimoto_median')} "
                f"| {fv(r, 'frac_test_with_nn_above_0.6')} |"
            )
        out.append("")

    path = RESULTS / "summary.md"
    path.write_text("\n".join(out) + "\n")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
