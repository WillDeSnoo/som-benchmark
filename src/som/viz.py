"""Shared plotting tokens and helpers.

Colours, surfaces and ink come from a validated categorical palette, used in
fixed slot order and never cycled. Three of the light-mode slots sit below 3:1
contrast against the chart surface, so every figure that uses more than two
series is accompanied by the corresponding CSV in results/ as its table view.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

#: categorical slots, in fixed order - never cycled, never reordered
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]

#: single hue for sequential/ordinal magnitude
SEQ = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

#: label colour when drawn inside a fill
ON_FILL = {"#eda100": "#0b0b0b"}

MODEL_LABELS = {
    "prior": "atom-type frequency",
    "fp_gbm": "fingerprint + GBM",
    "phys_gbm": "Hückel physics + GBM",
    "gnn": "message-passing GNN",
}

#: shorter forms for crowded category axes
MODEL_LABELS_SHORT = {
    "prior": "atom-type\nfrequency",
    "fp_gbm": "fingerprint\n+ GBM",
    "phys_gbm": "Hückel\n+ GBM",
    "gnn": "GNN",
}

MODE_LABELS = {"random": "random split", "scaffold": "scaffold split"}


def use_style() -> None:
    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.titlecolor": INK,
        "axes.labelsize": 9,
        "axes.labelcolor": INK_2,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK_2,
        "ytick.labelcolor": INK_2,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.frameon": False,
        "legend.fontsize": 8,
        "lines.linewidth": 2.0,
        "lines.markersize": 7,
        "figure.dpi": 150,
    })


def strip_frame(ax, keep=("left", "bottom")) -> None:
    for side, spine in ax.spines.items():
        spine.set_visible(side in keep)


def label_bars(
    ax, bars, values, fmt="{:.2f}", inside=False, color=None, tops=None
) -> None:
    """Direct value labels. Required wherever a low-contrast slot is used.

    `tops` is the height the label must clear - pass the upper CI ends so the
    text sits above the error-bar cap rather than on top of it.
    """
    for i, (b, v) in enumerate(zip(bars, values)):
        if v is None or v != v:  # NaN
            continue
        if inside:
            ax.text(
                b.get_x() + b.get_width() / 2, b.get_height() / 2,
                fmt.format(v), ha="center", va="center", fontsize=7.5,
                color=color or "#ffffff",
            )
            continue
        y = b.get_height()
        if tops is not None:
            t = tops[i]
            if t == t:  # not NaN
                y = max(y, t)
        ax.text(
            b.get_x() + b.get_width() / 2, y + 0.022,
            fmt.format(v), ha="center", va="bottom", fontsize=7.5,
            color=INK_2,
        )


def footnote(fig, text: str) -> None:
    fig.text(0.01, 0.005, text, fontsize=7, color=MUTED, ha="left", va="bottom")


def save(fig, path, pad=0.35) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", pad_inches=pad)
    plt.close(fig)
    print(f"  wrote {path}")
