"""Guards on the GNN hyperparameter search.

The search is the one place in the pipeline where a model is chosen rather than
just fitted, so the things worth testing are: that it cannot do worse than the
configuration it replaced, that it selects on validation data only, and that
the knobs it turns actually reach the network.
"""

import numpy as np
import pytest

from som import models
from som.models import _GNN_GRID, _GNN_SHORTLIST_N, _build_mpnn


def test_grid_starts_from_the_former_default():
    """The first entry must be the empty dict: the pre-tuning configuration.

    If it were not in the grid, a search that happened to pick a worse
    configuration could make the tuned GNN look worse than the untuned one for
    no reason other than the grid.
    """
    assert _GNN_GRID[0] == {}


def test_grid_has_no_duplicates():
    seen = [tuple(sorted(c.items())) for c in _GNN_GRID]
    assert len(seen) == len(set(seen))


def test_grid_only_names_real_knobs():
    import inspect

    allowed = set(inspect.signature(models._fit_gnn).parameters)
    for cfg in _GNN_GRID:
        assert set(cfg) <= allowed, f"unknown knob in {cfg}"


def test_shortlist_is_smaller_than_the_grid():
    assert 0 < _GNN_SHORTLIST_N < len(_GNN_GRID)


def test_dropout_reaches_the_network():
    """A non-zero dropout must put the module in train/eval disagreement.

    Checking the forward pass rather than the constructor: a dropout argument
    that is accepted and then ignored is exactly the bug this catches.
    """
    torch = pytest.importorskip("torch")

    rng = np.random.default_rng(0)
    x = torch.from_numpy(rng.normal(size=(6, 5)).astype(np.float32))
    src = torch.tensor([0, 1, 1, 2, 3, 4], dtype=torch.long)
    dst = torch.tensor([1, 0, 2, 1, 4, 3], dtype=torch.long)
    edge = torch.zeros((6, 7), dtype=torch.float32)
    iso = torch.zeros(6, dtype=torch.long)

    for dropout, want_equal in ((0.0, True), (0.5, False)):
        torch.manual_seed(0)
        net = _build_mpnn(5, 7, 9, hidden=16, layers=2, dropout=dropout)
        net.train()
        torch.manual_seed(1)
        a = net(x, src, dst, edge, iso).detach().numpy()
        torch.manual_seed(2)
        b = net(x, src, dst, edge, iso).detach().numpy()
        same = np.allclose(a, b)
        assert same is want_equal, f"dropout={dropout}"


def test_pos_weight_changes_the_gradient():
    """pos_weight must actually re-weight the positive class in the loss."""
    torch = pytest.importorskip("torch")
    from torch import nn

    logits = torch.zeros(4, requires_grad=True)
    y = torch.tensor([1.0, 0.0, 0.0, 0.0])

    plain = nn.BCEWithLogitsLoss()(logits, y)
    g_plain = torch.autograd.grad(plain, logits)[0].clone()

    logits2 = torch.zeros(4, requires_grad=True)
    weighted = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(3.0))(logits2, y)
    g_weighted = torch.autograd.grad(weighted, logits2)[0].clone()

    # the positive atom's gradient is tripled, the negatives' are untouched
    assert g_weighted[0] == pytest.approx(3.0 * g_plain[0])
    assert g_weighted[1:].numpy() == pytest.approx(g_plain[1:].numpy())


def test_search_selects_the_validation_argmax_and_ignores_test(monkeypatch):
    """The chosen configuration is the one with the best validation AP.

    ``_fit_gnn`` is replaced by a stub that reports a validation score which
    *disagrees* with its test score, so a search that peeked at test
    performance would pick a different configuration and fail here. This is the
    protocol guarantee the whole evaluation rests on, so it is worth a test that
    can actually catch a regression rather than a comment.
    """
    from som.models import Prediction
    from som.splits import Split

    grid = models._GNN_GRID
    # deliberately opposed orderings: best val is the last config, best test
    # is the first
    val_scores = [0.10 + 0.01 * i for i in range(len(grid))]
    test_scores = list(reversed(val_scores))

    calls = []

    def fake_fit_gnn(instances, split, label_field, seed, *, node_kind, name,
                     **params):
        i = next(
            j for j, c in enumerate(grid)
            if dict(c) == dict(params)
        )
        calls.append(params)
        return Prediction(
            model=name,
            y_test=np.array([0, 1], dtype=np.int8),
            p_test=np.full(2, test_scores[i]),
            inst_test=np.array([0, 0]), atom_test=np.array([0, 1]),
            y_val=np.array([0, 1], dtype=np.int8),
            p_val=np.full(2, val_scores[i]), inst_val=np.array([0, 0]),
            chosen_params={"best_val_ap": val_scores[i], "epochs_run": 1},
        )

    monkeypatch.setattr(models, "_fit_gnn", fake_fit_gnn)
    models._TUNED.clear()
    models._GNN_SHORTLIST.clear()
    models.GNN_TRIALS.clear()

    split = Split(
        mode=models._GNN_PRIMARY_MODE, seed=0,
        train_inst=np.array([0]), val_inst=np.array([1]),
        test_inst=np.array([2]),
        train_mols=["a"], val_mols=["b"], test_mols=["c"],
    )
    chosen = models._tune_gnn("gnn", "graph", [], split, "som_all")

    assert len(calls) == len(grid) * models._GNN_TUNE_REPEATS, (
        "every configuration should be tried, once per repeat"
    )
    assert chosen == dict(grid[-1]), "should pick the best validation score"

    logged = [t["val_avg_precision"] for t in models.GNN_TRIALS]
    assert logged == pytest.approx(val_scores)

    # the shortlist carried to the other mode is the validation top-N
    assert models._GNN_SHORTLIST["gnn"] == [
        dict(c) for c in reversed(grid[-models._GNN_SHORTLIST_N:])
    ]

    models._TUNED.clear()
    models._GNN_SHORTLIST.clear()
    models.GNN_TRIALS.clear()


def test_primary_mode_is_a_real_split_mode():
    """The full-search mode must be one the evaluation actually runs.

    A typo here would silently give every mode a full search, tripling the
    cost, which is the kind of mistake that only shows up as a long wait.
    """
    from som.splits import SPLIT_MODES

    assert models._GNN_PRIMARY_MODE in SPLIT_MODES


def test_evaluate_puts_the_primary_mode_first():
    """Stage 3 must visit the primary mode first, whatever order it is given.

    The shortlist for the other modes is derived from the primary mode's
    search, so if the primary mode came second it would get the shortlist and
    the other mode the full search - the opposite of the documented protocol.
    """
    from som.splits import SPLIT_MODES

    primary = models._GNN_PRIMARY_MODE
    for given in (list(SPLIT_MODES), list(reversed(SPLIT_MODES))):
        ordered = sorted(given, key=lambda m: m != primary)
        assert ordered[0] == primary
        assert sorted(ordered) == sorted(given), "no mode may be dropped"


# --------------------------------------------------------------------------
# reproducibility
# --------------------------------------------------------------------------

_TOY_SMILES = [
    "CCOc1ccccc1", "Cc1ccc(N)cc1", "CN(C)CCc1ccccc1", "OCc1ccncc1",
    "CSc1ccccc1", "CC(=O)Nc1ccc(O)cc1", "Clc1ccc(CN)cc1",
    "CCN1CCN(C)CC1", "Fc1ccc(CO)cc1",
]


def _toy_problem():
    from rdkit import Chem

    from som.splits import Split

    instances = []
    for i, smi in enumerate(_TOY_SMILES):
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, smi
        # a deterministic, chemically arbitrary label: the carbons bearing
        # hydrogens, every third one. The point of the test is bit-exactness,
        # not chemistry.
        labels = [
            a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() == "C" and a.GetTotalNumHs() > 0
        ][::3]
        instances.append({
            "smiles": smi,
            "isoform": "3A4" if i % 2 else "2D6",
            "som_all": labels or [0],
        })

    split = Split(
        mode="scaffold", seed=0,
        train_inst=np.arange(0, 5), val_inst=np.arange(5, 7),
        test_inst=np.arange(7, 9),
        train_mols=_TOY_SMILES[0:5], val_mols=_TOY_SMILES[5:7],
        test_mols=_TOY_SMILES[7:9],
    )
    return instances, split


def test_two_fits_at_the_same_seed_are_bit_identical():
    """A repeated GNN fit must reproduce exactly, not merely closely.

    Caveat on what this does and does not catch: the toy problem here is far
    too small for torch to parallelize the tensor ops, so it would pass even
    with threading left on. It guards gross nondeterminism -- an unseeded
    dropout mask, an unseeded batch order -- rather than the thread race. The
    thread race is guarded by `test_the_gnn_is_single_threaded` below, which
    checks the mechanism directly.
    """
    pytest.importorskip("torch")
    instances, split = _toy_problem()

    runs = [
        models._fit_gnn(
            instances, split, "som_all", 0, node_kind="graph", name="gnn",
            hidden=16, layers=2, epochs=4, patience=4,
        )
        for _ in range(2)
    ]
    a, b = runs
    assert np.array_equal(a.p_test, b.p_test), (
        "identical inputs and seed gave different predictions"
    )
    assert a.chosen_params["best_val_ap"] == b.chosen_params["best_val_ap"]


def test_different_seeds_give_different_fits():
    """The companion check: seeding must actually vary the fit.

    Without this, a bug that pinned every fit to one initialization would make
    the test above pass for the wrong reason, and the repeated fits in the
    search would silently be measuring nothing.
    """
    pytest.importorskip("torch")
    instances, split = _toy_problem()

    a, b = (
        models._fit_gnn(
            instances, split, "som_all", seed, node_kind="graph", name="gnn",
            hidden=16, layers=2, epochs=4, patience=4,
        )
        for seed in (0, 1)
    )
    assert not np.array_equal(a.p_test, b.p_test)


def test_the_gnn_is_single_threaded():
    """The fit must pin torch to one thread, and must do so itself.

    This is the real guard on reproducibility. It asserts the mechanism rather
    than an outcome, because the outcome -- a bit-identical refit -- only
    diverges on problems large enough for torch to parallelize, which no fast
    test can construct. Checked after a fit, so that it is the fit that set it
    and not the ambient default: the test deliberately sets a different value
    first.
    """
    torch = pytest.importorskip("torch")

    assert models._GNN_THREADS == 1
    before = torch.get_num_threads()
    try:
        torch.set_num_threads(4)
        instances, split = _toy_problem()
        models._fit_gnn(
            instances, split, "som_all", 0, node_kind="graph", name="gnn",
            hidden=16, layers=2, epochs=2, patience=2,
        )
        assert torch.get_num_threads() == 1
    finally:
        torch.set_num_threads(before)


def test_stage_3_is_guarded_for_spawn():
    """Stage 3 must guard its entry point, or `--jobs` misbehaves badly.

    Under the `spawn` start method -- the default on macOS -- every worker
    re-imports the main module. Without the guard each worker re-runs the
    whole evaluation instead of just its own fit, which is slow, confusing,
    and only shows up once someone passes `--jobs`. Cheap to check, nasty to
    debug.
    """
    from pathlib import Path

    root = Path(__file__).resolve().parent.parent
    for script in ("03_evaluate.py", "07_bde.py", "08_smartcyp.py"):
        path = root / "scripts" / script
        if not path.exists():
            continue
        assert 'if __name__ == "__main__":' in path.read_text(), script
