#!/usr/bin/env python
"""Stage 5: figures, from the CSVs written by stages 3 and 4.

Every figure with more than two series has the matching CSV in results/ as its
table view, so identity never rests on colour alone.
"""

import argparse
import csv
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from som.models import MODELS
from som.paths import FIGURES, RESULTS
from som.viz import (
    AXIS,
    INK_2,
    MODEL_LABELS,
    MODE_LABELS,
    MUTED,
    ON_FILL,
    SERIES,
    footnote,
    label_bars,
    save,
    strip_frame,
    use_style,
)

ORDER = ["prior", "fp_gbm", "phys_gbm", "gnn"]


def read(path):
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def f(row, key, default=float("nan")):
    try:
        return float(row[key])
    except (KeyError, ValueError, TypeError):
        return default


def models_present(rows):
    seen = {r["model"] for r in rows}
    return [m for m in ORDER if m in seen]


# --------------------------------------------------------------------------

def fig_main(summary):
    """Top-1 and top-2 per model, random vs scaffold split. Two series only."""
    present = models_present(summary)
    lookup = {(r["model"], r["mode"]): r for r in summary}
    modes = ["random", "scaffold"]

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))
    x = np.arange(len(present))
    width = 0.36

    for ax, metric, title in zip(
        axes, ["top1", "top2"],
        ["Top-1 accuracy", "Top-2 accuracy"],
    ):
        for s, mode in enumerate(modes):
            means, los, his = [], [], []
            for m in present:
                r = lookup.get((m, mode))
                means.append(f(r, f"{metric}_mean") if r else float("nan"))
                los.append(f(r, f"{metric}_lo") if r else float("nan"))
                his.append(f(r, f"{metric}_hi") if r else float("nan"))
            means = np.asarray(means)
            err = np.vstack([means - np.asarray(los), np.asarray(his) - means])
            err[~np.isfinite(err)] = 0.0
            # 2px surface gap between adjacent bars
            bars = ax.bar(
                x + (s - 0.5) * (width + 0.018), means, width,
                color=SERIES[s], label=MODE_LABELS[mode],
                edgecolor="#fcfcfb", linewidth=1.0, zorder=3,
            )
            ax.errorbar(
                x + (s - 0.5) * (width + 0.018), means, yerr=err,
                fmt="none", ecolor=INK_2, elinewidth=1.0, capsize=2.5, zorder=4,
            )
            label_bars(ax, bars, means, fmt="{:.2f}")

        ax.set_xticks(x)
        ax.set_xticklabels(
            [MODEL_LABELS[m].replace(" + ", "\n+ ") for m in present], fontsize=8
        )
        ax.set_ylim(0, 1.0)
        ax.set_ylabel("per-molecule accuracy")
        ax.set_title(title, loc="left")
        ax.grid(axis="x", visible=False)
        strip_frame(ax)

    axes[0].legend(loc="upper left", ncols=2, bbox_to_anchor=(0, 1.0))
    fig.suptitle(
        "Finding the site of metabolism: four approaches, two splits",
        x=0.008, ha="left", fontsize=12, fontweight="semibold",
    )
    footnote(
        fig,
        "Mean over 10 repeated splits; bars show the 95% CI of the mean. "
        "Error bars fold in both split and training variability.",
    )
    fig.tight_layout(rect=(0, 0.03, 1, 0.94))
    save(fig, FIGURES / "fig1_model_comparison.png")


def fig_per_isoform(rows):
    """Dot plot: top-2 by CYP isoform, scaffold split."""
    rows = [r for r in rows if r["mode"] == "scaffold"]
    if not rows:
        return
    present = models_present(rows)
    isoforms = sorted({r["isoform"] for r in rows})
    by = {(r["model"], r["isoform"]): r for r in rows}

    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    y = np.arange(len(isoforms))
    offs = np.linspace(-0.22, 0.22, len(present))

    for i, m in enumerate(present):
        xs = [f(by.get((m, iso)), "top2_mean") for iso in isoforms]
        los = [f(by.get((m, iso)), "top2_lo") for iso in isoforms]
        his = [f(by.get((m, iso)), "top2_hi") for iso in isoforms]
        ax.hlines(
            y + offs[i], los, his, color=SERIES[i], linewidth=1.6,
            alpha=0.55, zorder=2,
        )
        ax.scatter(
            xs, y + offs[i], s=58, color=SERIES[i], label=MODEL_LABELS[m],
            edgecolor="#fcfcfb", linewidth=1.4, zorder=3,
        )

    n_lab = {
        iso: f(by.get((present[0], iso)), "mean_test_instances")
        for iso in isoforms
    }
    ax.set_yticks(y)
    ax.set_yticklabels([f"CYP{iso}" for iso in isoforms])
    ax.set_ylim(-0.6, len(isoforms) - 0.4)
    ax.invert_yaxis()
    ax.set_xlabel("top-2 accuracy (scaffold split)")
    ax.set_xlim(0, 1.0)
    ax.grid(axis="y", visible=False)
    strip_frame(ax)

    for yi, iso in zip(y, isoforms):
        ax.text(
            1.005, yi, f"n≈{n_lab[iso]:.0f}", transform=ax.get_yaxis_transform(),
            va="center", fontsize=7, color=MUTED,
        )

    ax.legend(loc="lower left", ncols=2, bbox_to_anchor=(0, -0.02))
    ax.set_title(
        "Per-isoform performance, and how uneven the support is",
        loc="left",
    )
    footnote(
        fig,
        "Mean over 10 scaffold splits; lines are 95% CIs. n is the mean number "
        "of test instances per split. Table view: results/per_isoform.csv",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, FIGURES / "fig2_per_isoform.png")


def fig_site_type(rows, min_atoms=30):
    """Recall at top-2 by the transformation the labelled atom implies."""
    rows = [r for r in rows if r["mode"] == "scaffold"]
    if not rows:
        return
    show = [m for m in ["fp_gbm", "phys_gbm", "gnn"]
            if m in {r["model"] for r in rows}]
    by = {(r["model"], r["site_type"]): r for r in rows}

    counts = {}
    for r in rows:
        counts[r["site_type"]] = max(
            counts.get(r["site_type"], 0), int(f(r, "n_positive_atoms_pooled", 0))
        )
    types = [t for t, n in sorted(counts.items(), key=lambda kv: -kv[1])
             if n >= min_atoms]
    if not types:
        return

    fig, ax = plt.subplots(figsize=(8.2, 0.52 * len(types) + 2.0))
    y = np.arange(len(types))
    h = 0.78 / len(show)

    for i, m in enumerate(show):
        vals = [f(by.get((m, t)), "recall_at_2_pooled") for t in types]
        ax.barh(
            y + (i - (len(show) - 1) / 2) * (h + 0.012), vals, h,
            color=SERIES[i], label=MODEL_LABELS[m],
            edgecolor="#fcfcfb", linewidth=1.0, zorder=3,
        )

    ax.set_yticks(y)
    ax.set_yticklabels(types, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("fraction of true sites ranked in the molecule's top 2")
    ax.grid(axis="y", visible=False)
    strip_frame(ax)

    for yi, t in zip(y, types):
        ax.text(
            1.005, yi, f"n={counts[t]}", transform=ax.get_yaxis_transform(),
            va="center", fontsize=7, color=MUTED,
        )

    ax.legend(loc="lower right", ncols=1)
    ax.set_title(
        "Which metabolic chemistry each model actually finds", loc="left"
    )
    footnote(
        fig,
        f"Scaffold split, positive atoms pooled over 10 seeds; site types with "
        f"at least {min_atoms} labelled atoms. Table view: "
        "results/per_site_type.csv",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, FIGURES / "fig3_site_type.png")


def fig_novelty(rows):
    """Accuracy against similarity of the test molecule to the training set."""
    pooled = [r for r in rows if r["mode"] == "pooled"]
    rows = pooled or [r for r in rows if r["mode"] == "scaffold"]
    if not rows:
        return
    present = models_present(rows)
    buckets = []
    for r in rows:
        if r["nn_tanimoto_bucket"] not in buckets:
            buckets.append(r["nn_tanimoto_bucket"])
    order = ["<=0.30", "0.30-0.45", "0.45-0.60", "0.60-0.75", ">0.75"]
    buckets = [b for b in order if b in buckets]
    by = {(r["model"], r["nn_tanimoto_bucket"]): r for r in rows}

    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    x = np.arange(len(buckets))
    for i, m in enumerate(present):
        ys = [f(by.get((m, b)), "top1_mean") for b in buckets]
        lo = [f(by.get((m, b)), "top1_lo") for b in buckets]
        hi = [f(by.get((m, b)), "top1_hi") for b in buckets]
        ax.fill_between(x, lo, hi, color=SERIES[i], alpha=0.13, linewidth=0)
        ax.plot(x, ys, color=SERIES[i], marker="o", markersize=7,
                markeredgecolor="#fcfcfb", markeredgewidth=1.4,
                label=MODEL_LABELS[m], zorder=3)

    ns = [f(by.get((present[0], b)), "mean_test_instances") for b in buckets]
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{b}\nn≈{n:.0f}" for b, n in zip(buckets, ns)], fontsize=8
    )
    ax.set_xlabel(
        "Tanimoto similarity of the test molecule to its nearest training neighbour"
    )
    ax.set_ylabel("top-1 accuracy")
    ax.set_ylim(0, 1.0)
    ax.grid(axis="x", visible=False)
    strip_frame(ax)
    ax.legend(loc="upper left", ncols=2)
    ax.set_title(
        "Performance decays with distance from the training set", loc="left"
    )
    footnote(
        fig,
        "Test instances from both split types pooled per seed, then averaged "
        "over 10 seeds; bands are 95% CIs. Table view: "
        "results/novelty_curve.csv",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, FIGURES / "fig4_novelty.png")


def fig_calibration(rows):
    """Reliability curves on the scaffold split."""
    rows = [r for r in rows if r["mode"] == "scaffold"]
    if not rows:
        return
    present = models_present(rows)
    by = defaultdict(list)
    for r in rows:
        by[r["model"]].append(r)

    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    ax.plot([0, 0.75], [0, 0.75], color=AXIS, linewidth=1.2,
            linestyle=(0, (4, 3)), zorder=1)
    ax.text(0.60, 0.645, "perfect calibration", fontsize=7.5, color=MUTED,
            rotation=38, ha="center", va="center")

    for i, m in enumerate(present):
        pts = sorted(by[m], key=lambda r: int(r["bin"]))
        xs = [f(r, "mean_predicted") for r in pts]
        ys = [f(r, "observed_frequency") for r in pts]
        ax.plot(xs, ys, color=SERIES[i], marker="o", markersize=6.5,
                markeredgecolor="#fcfcfb", markeredgewidth=1.3,
                label=MODEL_LABELS[m], zorder=3)

    ax.set_xlabel("mean predicted probability")
    ax.set_ylabel("observed fraction of true sites")
    ax.set_xlim(0, 0.78)
    ax.set_ylim(0, 0.78)
    ax.set_aspect("equal")
    strip_frame(ax)
    ax.legend(loc="upper left")
    ax.set_title("Are the probabilities believable?", loc="left")
    footnote(
        fig,
        "Equal-count bins, scaffold split, averaged over 10 seeds. "
        "Table view: results/calibration.csv",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    save(fig, FIGURES / "fig5_calibration.png")


def main() -> None:
    argparse.ArgumentParser(description=__doc__).parse_args()
    use_style()
    FIGURES.mkdir(parents=True, exist_ok=True)
    print("writing figures:")
    fig_main(read(RESULTS / "metrics_summary.csv"))
    fig_per_isoform(read(RESULTS / "per_isoform.csv"))
    fig_site_type(read(RESULTS / "per_site_type.csv"))
    fig_novelty(read(RESULTS / "novelty_curve.csv"))
    fig_calibration(read(RESULTS / "calibration.csv"))


if __name__ == "__main__":
    main()
