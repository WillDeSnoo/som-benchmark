"""Metric behaviour on hand-built cases, where the right answer is obvious."""

import numpy as np
import pytest

from som.metrics import (
    bootstrap_ci,
    evaluate,
    expected_calibration_error,
    mean_ci_across_seeds,
    topk_hits,
)


def test_topk_perfect_and_worst():
    # two instances of four atoms each, one positive each
    y = np.array([1, 0, 0, 0, 0, 0, 0, 1])
    inst = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    perfect = np.array([0.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.9])
    hits = topk_hits(y, perfect, inst, ks=(1, 2))
    assert hits[1].tolist() == [True, True]

    worst = 1.0 - perfect
    hits = topk_hits(y, worst, inst, ks=(1, 2))
    assert hits[1].tolist() == [False, False]
    # the positive is ranked last of four, so even top-3 misses it
    assert topk_hits(y, worst, inst, ks=(3,))[3].tolist() == [False, False]


def test_topk_counts_any_positive_in_the_top_k():
    y = np.array([0, 1, 0, 0])
    inst = np.zeros(4, dtype=int)
    p = np.array([0.9, 0.8, 0.2, 0.1])
    assert topk_hits(y, p, inst, ks=(1,))[1].tolist() == [False]
    assert topk_hits(y, p, inst, ks=(2,))[2].tolist() == [True]


def test_instances_without_positives_are_skipped():
    y = np.array([1, 0, 0, 0])
    inst = np.array([0, 0, 1, 1])
    p = np.array([0.9, 0.1, 0.5, 0.4])
    hits = topk_hits(y, p, inst, ks=(1,))
    assert hits[1].size == 1


def test_pessimistic_tie_policy_is_a_lower_bound():
    """All atoms tied: optimistic ordering could call this a hit, the
    pessimistic one must not."""
    y = np.array([1, 0, 0, 0])
    inst = np.zeros(4, dtype=int)
    p = np.full(4, 0.25)
    pess = topk_hits(y, p, inst, ks=(1,), tie_policy="pessimistic")
    assert pess[1].tolist() == [False]
    # with a seeded random tie-break the hit rate should be about 1/4 over
    # many draws, i.e. strictly between the pessimistic and optimistic bounds
    rate = np.mean([
        topk_hits(y, p, inst, ks=(1,), rng=np.random.default_rng(s))[1][0]
        for s in range(400)
    ])
    assert 0.15 < rate < 0.35, rate


def test_random_tie_break_is_reproducible_given_the_seed():
    y = np.array([1, 0, 0, 0])
    inst = np.zeros(4, dtype=int)
    p = np.full(4, 0.25)
    a = topk_hits(y, p, inst, ks=(1,), rng=np.random.default_rng(11))
    b = topk_hits(y, p, inst, ks=(1,), rng=np.random.default_rng(11))
    assert a[1].tolist() == b[1].tolist()


def test_ece_is_zero_for_a_perfectly_calibrated_predictor():
    rng = np.random.default_rng(0)
    p = rng.uniform(0, 1, 40000)
    y = (rng.uniform(0, 1, 40000) < p).astype(int)
    ece, curve = expected_calibration_error(y, p, n_bins=10)
    assert ece < 0.02, ece
    assert len(curve) == 10
    for b in curve:
        assert abs(b["mean_predicted"] - b["observed_frequency"]) < 0.05


def test_ece_detects_systematic_overconfidence():
    rng = np.random.default_rng(1)
    truth = rng.uniform(0, 1, 20000)
    y = (rng.uniform(0, 1, 20000) < truth).astype(int)
    # claim far more confidence than warranted
    p = np.clip(truth * 1.8, 0, 1)
    ece, _ = expected_calibration_error(y, p)
    assert ece > 0.1, ece


def test_equal_count_bins_are_actually_equal_count():
    p = np.concatenate([np.full(900, 0.01), np.linspace(0.2, 0.9, 100)])
    y = (np.random.default_rng(0).uniform(size=p.size) < p).astype(int)
    _, curve = expected_calibration_error(y, p, n_bins=10)
    sizes = [b["n"] for b in curve]
    assert max(sizes) - min(sizes) <= 1, sizes


def test_evaluate_reports_expected_keys():
    y = np.array([1, 0, 0, 0, 0, 1, 0, 0])
    inst = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    p = np.array([0.9, 0.2, 0.1, 0.05, 0.1, 0.8, 0.2, 0.1])
    res = evaluate(y, p, inst, with_curve=True)
    for key in ("top1", "top2", "top1_pessimistic", "roc_auc",
                "avg_precision", "brier", "ece", "reliability_curve"):
        assert key in res
    assert res["top1"] == pytest.approx(1.0)
    assert res["n_instances"] == 2
    assert res["positive_rate"] == pytest.approx(0.25)


def test_bootstrap_interval_brackets_the_point_estimate():
    rng = np.random.default_rng(3)
    n_mol = 120
    y, p, inst = [], [], []
    mol_of_inst = {}
    for m in range(n_mol):
        k = 8
        lab = np.zeros(k, dtype=int)
        lab[0] = 1
        score = rng.uniform(0, 1, k)
        score[0] += 0.8  # usually, but not always, ranked first
        y.append(lab)
        p.append(score)
        inst.append(np.full(k, m))
        mol_of_inst[m] = f"mol{m}"
    y = np.concatenate(y)
    p = np.concatenate(p)
    inst = np.concatenate(inst)

    point = evaluate(y, p, inst)["top1"]
    ci = bootstrap_ci(y, p, inst, mol_of_inst, metrics=("top1",), n_boot=200)
    lo, hi = ci["top1"]
    assert lo <= point <= hi, (lo, point, hi)
    assert hi - lo < 0.3


def test_mean_ci_across_seeds():
    s = mean_ci_across_seeds([0.5, 0.52, 0.48, 0.51, 0.49])
    assert s["mean"] == pytest.approx(0.5, abs=1e-9)
    assert s["lo"] < s["mean"] < s["hi"]
    assert s["n"] == 5
    single = mean_ci_across_seeds([0.7])
    assert single["lo"] == single["hi"] == 0.7
    assert mean_ci_across_seeds([])["n"] == 0


def test_per_atom_and_per_instance_views_agree():
    """`analysis.atom_ranks` must break ties the same way `topk_hits` does.

    Stage 4 reads the per-atom ranks while stage 3 reports the per-instance
    top-k, both from the same saved probabilities. If the two tie-break
    streams diverged, the failure analysis would contradict the headline
    numbers for no visible reason.
    """
    from som.analysis import instance_topk_table

    rng = np.random.default_rng(5)
    y, p, inst = [], [], []
    for m in range(60):
        k = int(rng.integers(4, 12))
        lab = np.zeros(k, dtype=int)
        lab[rng.integers(0, k)] = 1
        # deliberately coarse scores, so ties are common
        y.append(lab)
        p.append(np.round(rng.uniform(0, 1, k), 1))
        inst.append(np.full(k, m))
    y = np.concatenate(y)
    p = np.concatenate(p)
    inst = np.concatenate(inst)

    for seed in (0, 1, 4):
        hits = topk_hits(y, p, inst, ks=(1, 2), rng=np.random.default_rng(seed))
        table = instance_topk_table(y, p, inst, ks=(1, 2), seed=seed)
        for k in (1, 2):
            from_table = [table[i][f"top{k}"] for i in sorted(table)]
            assert from_table == hits[k].tolist(), (seed, k)


def test_nan_values_are_dropped_not_propagated():
    s = mean_ci_across_seeds([0.4, float("nan"), 0.6])
    assert s["n"] == 2
    assert s["mean"] == pytest.approx(0.5)
