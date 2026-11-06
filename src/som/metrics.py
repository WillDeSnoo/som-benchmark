"""Evaluation metrics for atom-level site-of-metabolism prediction.

Reported quantities:

* **top-1 / top-2 per-instance accuracy** -- the fraction of (molecule, isoform)
  instances for which at least one true site appears among the k
  highest-scoring atoms. This is the number the SOM literature quotes, and it
  is the one a medicinal chemist actually cares about: "did it point me at the
  right position".
* **atom-level ROC-AUC and average precision** -- ranking quality over all
  candidate atoms. With roughly 10% positives, average precision is the more
  honest of the two, and ROC-AUC is included because everyone else reports it.
* **calibration** -- expected calibration error, Brier score, and a reliability
  curve, because a model whose probabilities mean nothing is much less useful
  for prioritising experiments even if its ranking is fine.

Two choices worth stating plainly:

**Ties.** Several atoms frequently receive identical scores, especially
symmetry-equivalent ones. Breaking ties in index order silently flatters a
model whose first atoms happen to be the labelled ones. Ties are therefore
broken by a seeded random permutation, which is unbiased in expectation; the
pessimistic variant (a tie counts as a miss) is also computed so the gap is
visible rather than hidden.

**Bootstrap unit.** Resampling is over *molecules*, not atoms and not
instances. Atoms within a molecule are strongly dependent and instances of one
molecule share a structure, so resampling either of those would give
intervals that are far too narrow.
"""

from __future__ import annotations

import numpy as np
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

TOP_KS = (1, 2, 3)


# --------------------------------------------------------------------------
# per-instance top-k
# --------------------------------------------------------------------------

def _instance_slices(inst_id: np.ndarray) -> dict[int, np.ndarray]:
    order = np.argsort(inst_id, kind="stable")
    sorted_ids = inst_id[order]
    bounds = np.flatnonzero(np.diff(sorted_ids)) + 1
    out = {}
    for chunk in np.split(order, bounds):
        if chunk.size:
            out[int(inst_id[chunk[0]])] = chunk
    return out


def topk_hits(
    y: np.ndarray,
    p: np.ndarray,
    inst_id: np.ndarray,
    ks=TOP_KS,
    rng: np.random.Generator | None = None,
    tie_policy: str = "random",
) -> dict[int, np.ndarray]:
    """Per-instance top-k hit indicators.

    Returns {k: bool array over instances}, ordered by sorted instance id.
    Instances with no positive atom are skipped (they carry no information
    about ranking); in this dataset there are none, since every record has at
    least one PRIMARY_SOM.
    """
    rng = rng or np.random.default_rng(0)
    slices = _instance_slices(inst_id)
    out = {k: [] for k in ks}
    for iid in sorted(slices):
        idx = slices[iid]
        yy, pp = y[idx], p[idx]
        n = pp.size
        # The draw happens before the skip, so that the random stream depends
        # only on the instance sizes and not on which instances have labels.
        # `analysis.atom_ranks` draws the same way, which is what keeps the
        # per-atom and per-instance views of the same predictions consistent.
        jitter = rng.permutation(n)
        if yy.sum() == 0:
            continue
        if tie_policy == "random":
            order = np.lexsort((jitter, -pp))
        elif tie_policy == "pessimistic":
            # among equal scores, put the positives last
            order = np.lexsort((yy, -pp))
        else:
            raise ValueError(f"unknown tie_policy {tie_policy!r}")
        ranked = yy[order]
        for k in ks:
            out[k].append(bool(ranked[: min(k, n)].sum() > 0))
    return {k: np.asarray(v, dtype=bool) for k, v in out.items()}


# --------------------------------------------------------------------------
# calibration
# --------------------------------------------------------------------------

def expected_calibration_error(
    y: np.ndarray, p: np.ndarray, n_bins: int = 10
) -> tuple[float, list[dict]]:
    """ECE with equal-count bins, plus the reliability curve.

    Equal-count rather than equal-width bins, because predicted probabilities
    for this task pile up near zero and equal-width bins would leave most bins
    nearly empty and the ECE dominated by one of them.
    """
    n = p.size
    if n == 0:
        return float("nan"), []
    order = np.argsort(p, kind="stable")
    bins = np.array_split(order, min(n_bins, n))
    ece = 0.0
    curve = []
    for b in bins:
        if b.size == 0:
            continue
        conf = float(p[b].mean())
        acc = float(y[b].mean())
        ece += (b.size / n) * abs(acc - conf)
        curve.append({
            "n": int(b.size),
            "mean_predicted": conf,
            "observed_frequency": acc,
        })
    return float(ece), curve


# --------------------------------------------------------------------------
# full evaluation
# --------------------------------------------------------------------------

def evaluate(
    y: np.ndarray,
    p: np.ndarray,
    inst_id: np.ndarray,
    ks=TOP_KS,
    seed: int = 0,
    with_curve: bool = False,
) -> dict:
    """All metrics for one set of predictions."""
    rng = np.random.default_rng(seed)
    res: dict = {}

    hits = topk_hits(y, p, inst_id, ks=ks, rng=rng, tie_policy="random")
    for k, v in hits.items():
        res[f"top{k}"] = float(v.mean()) if v.size else float("nan")
    pess = topk_hits(y, p, inst_id, ks=ks, tie_policy="pessimistic")
    for k, v in pess.items():
        res[f"top{k}_pessimistic"] = float(v.mean()) if v.size else float("nan")

    res["n_instances"] = int(hits[ks[0]].size)
    res["n_atoms"] = int(y.size)
    res["positive_rate"] = float(y.mean()) if y.size else float("nan")

    if y.size and 0 < y.sum() < y.size:
        res["roc_auc"] = float(roc_auc_score(y, p))
        res["avg_precision"] = float(average_precision_score(y, p))
    else:
        res["roc_auc"] = float("nan")
        res["avg_precision"] = float("nan")

    pc = np.clip(p, 0.0, 1.0)
    res["brier"] = float(brier_score_loss(y, pc)) if y.size else float("nan")
    ece, curve = expected_calibration_error(y, pc)
    res["ece"] = ece
    if with_curve:
        res["reliability_curve"] = curve
    return res


# --------------------------------------------------------------------------
# bootstrap over molecules
# --------------------------------------------------------------------------

def bootstrap_ci(
    y: np.ndarray,
    p: np.ndarray,
    inst_id: np.ndarray,
    mol_of_inst: dict[int, str],
    metrics=("top1", "top2", "roc_auc", "avg_precision"),
    n_boot: int = 500,
    alpha: float = 0.05,
    seed: int = 0,
) -> dict[str, tuple[float, float]]:
    """Percentile bootstrap CIs, resampling molecules with replacement."""
    rng = np.random.default_rng(seed)
    by_mol: dict[str, list[int]] = {}
    for iid in np.unique(inst_id):
        by_mol.setdefault(mol_of_inst[int(iid)], []).append(int(iid))
    mols = sorted(by_mol)
    if len(mols) < 2:
        return {m: (float("nan"), float("nan")) for m in metrics}

    rows_of_inst = _instance_slices(inst_id)
    draws: dict[str, list[float]] = {m: [] for m in metrics}

    for _ in range(n_boot):
        picked = rng.choice(len(mols), size=len(mols), replace=True)
        rows, ids = [], []
        # resampled molecules need distinct instance ids, or repeated draws of
        # the same molecule would be merged into one instance
        next_id = 0
        for j in picked:
            for iid in by_mol[mols[j]]:
                r = rows_of_inst[iid]
                rows.append(r)
                ids.append(np.full(r.size, next_id, dtype=np.int64))
                next_id += 1
        if not rows:
            continue
        ridx = np.concatenate(rows)
        rid = np.concatenate(ids)
        m = evaluate(y[ridx], p[ridx], rid, seed=int(rng.integers(1 << 30)))
        for name in metrics:
            draws[name].append(m.get(name, float("nan")))

    out = {}
    for name, vals in draws.items():
        arr = np.asarray(vals, dtype=np.float64)
        arr = arr[np.isfinite(arr)]
        if arr.size < 10:
            out[name] = (float("nan"), float("nan"))
        else:
            out[name] = (
                float(np.percentile(arr, 100 * alpha / 2)),
                float(np.percentile(arr, 100 * (1 - alpha / 2))),
            )
    return out


# --------------------------------------------------------------------------
# aggregation across seeds
# --------------------------------------------------------------------------

def mean_ci_across_seeds(values: list[float], alpha: float = 0.05) -> dict:
    """Mean and t-based CI of the mean across repeated splits.

    This is the headline interval: it folds in both the variability of the
    split and the variability of training, which is what someone asking "how
    well does this method work" wants to know. A bootstrap CI computed inside a
    single split answers a narrower question and is reported separately.
    """
    from scipy import stats

    arr = np.asarray([v for v in values if np.isfinite(v)], dtype=np.float64)
    if arr.size == 0:
        return {"mean": float("nan"), "lo": float("nan"), "hi": float("nan"), "n": 0}
    if arr.size == 1:
        return {
            "mean": float(arr[0]), "sd": 0.0,
            "lo": float(arr[0]), "hi": float(arr[0]), "n": 1,
        }
    mean = float(arr.mean())
    sd = float(arr.std(ddof=1))
    half = float(stats.t.ppf(1 - alpha / 2, arr.size - 1) * sd / np.sqrt(arr.size))
    return {
        "mean": mean, "sd": sd, "lo": mean - half, "hi": mean + half,
        "n": int(arr.size),
    }
