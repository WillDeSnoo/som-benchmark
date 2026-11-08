"""The four modelling approaches.

``prior``
    Not a model so much as a reference point: the training-set frequency of
    being a site of metabolism for each coarse atom type
    (element, number of hydrogens, aromatic, in-ring), with Laplace smoothing.
    It learns nothing about molecular context. Including it is the difference
    between "our model reaches top-1 of X" and knowing how much of X is just
    "CYPs attack methyl groups and benzylic positions". Any approach that does
    not clear this by a wide margin is not doing useful work.

``fp_gbm``
    Gradient boosting on hashed circular atom environments. The standard
    cheap-and-strong cheminformatics baseline.

``phys_gbm``
    Gradient boosting on the small interpretable physics descriptor set
    (Hückel frontier densities and free valence, Gasteiger charges, E-state,
    Crippen contributions, valence/ring context). A linear model on the same
    features is also fitted, purely so its coefficients can be read.

``gnn``
    A message-passing network over the molecular graph, written in plain
    PyTorch to avoid a PyTorch Geometric dependency. It gets only basic atom
    and bond features, *not* the Hückel descriptors, so that "learned
    representation" is compared against "hand-built physics" rather than
    quietly containing it. Width, depth, learning rate, dropout, weight decay
    and positive-class weight are searched; each GNN variant gets its own
    search, so the ablations below compare equally-tuned models.

``gnn_huckel``
    The same network given exactly the feature set ``phys_gbm`` uses, Hückel
    descriptors included. This is the ablation that says whether the physics
    and the learned representation carry the *same* information or different
    information: if the GNN can already infer frontier-orbital structure by
    passing messages, adding it explicitly should change nothing, and if it
    cannot, this should close the gap to ``phys_gbm``. Note that the physics
    block includes neighbour-aggregated Hückel terms, which message passing
    could in principle derive for itself; they are kept so that the GNN and the
    GBM see an identical feature set and the comparison is about the model
    rather than the inputs.

Hyperparameters are chosen on the validation split of the first seed only and
then reused for every later seed. The test set is never touched by tuning.

The GBM and the GNN get the same protocol but not the same budget. The GBM has
four configurations and one fit each; the GNN has twelve and three fits each,
because a single GNN fit is a noisy measurement of a configuration's quality
(see ``_GNN_TUNE_REPEATS``) while a GBM fit at a fixed seed is not. Both search
on validation average precision. The search is written out to
results/gnn_tuning.csv so that it can be read rather than taken on trust.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .features import AtomTable, build_atom_table, isoform_onehot
from .splits import Split, atom_mask

MODELS = (
    "prior", "smartcyp", "fp_gbm", "phys_gbm", "gnn", "gnn_huckel",
    "phys_bde_gbm", "gnn_bde",
)

#: SMARTCyp ships isoform-specific variants for 2D6 and for the 2C subfamily;
#: everything else uses its general score.
_SMARTCYP_COLUMN = {
    "2D6": "score_2d6",
    "2C8": "score_2c", "2C9": "score_2c", "2C19": "score_2c",
}

#: cache of featurized tables, keyed by (feature kind, label field)
_TABLE_CACHE: dict[tuple[str, str], AtomTable] = {}

#: hyperparameters chosen on seed 0, keyed by (model, split mode)
_TUNED: dict[tuple[str, str], dict] = {}

_GBM_GRID = [
    {"learning_rate": 0.05, "max_leaf_nodes": 31},
    {"learning_rate": 0.05, "max_leaf_nodes": 63},
    {"learning_rate": 0.10, "max_leaf_nodes": 31},
    {"learning_rate": 0.10, "max_leaf_nodes": 63},
]

#: Coarse search space for the GNN. The first entry is the configuration used
#: before any tuning was done, so the search can only improve on it (on the
#: validation split, which is what it selects on).
#:
#: It is a hand-chosen path through the space rather than a full factorial:
#: one-at-a-time moves away from the base along each axis that plausibly
#: matters at 679 molecules, plus two combinations of the moves that small-data
#: graph models usually want together (smaller, shallower, regularized). A full
#: factorial over these six axes would be a few hundred fits at ~90 s each,
#: times ``_GNN_TUNE_REPEATS``, times three variants and two split modes.
#:
#: Capacity (``hidden``, ``layers``) and regularization (``dropout``,
#: ``weight_decay``) are the axes expected to matter; ``pos_weight`` is here
#: because ~10% of atoms are positive and the original loss ignored that.
#: Moves *upward* in capacity (``hidden`` 256, ``layers`` 6) were left out:
#: they are the most expensive fits and the least likely to help a model that
#: is already the weaker estimator on this sample size, so the budget went to
#: repeats instead. That is a stated limit of the search, not a finding.
_GNN_GRID = [
    {},                                                    # the former default
    {"hidden": 64},
    {"layers": 2},
    {"layers": 3},
    {"lr": 3e-4},
    {"lr": 3e-3},
    {"dropout": 0.1},
    {"dropout": 0.25},
    {"weight_decay": 1e-3},
    {"pos_weight": 3.0},
    {"hidden": 64, "layers": 3, "dropout": 0.1, "weight_decay": 1e-3},
    {"hidden": 128, "layers": 3, "dropout": 0.1, "pos_weight": 3.0},
]

#: Threads for a GNN fit. One, deliberately.
#:
#: `index_add_` over several CPU threads accumulates in a nondeterministic
#: order, and those rounding differences grow through training into validation
#: scores that differ by ~0.013 AP between repeats of the *same* configuration
#: at the *same* seed -- larger than the genuine seed-to-seed spread (~0.009)
#: and comparable to the gaps between neighbouring grid points. A search run
#: that way is partly a search over thread scheduling.
#:
#: Single-threaded is bit-exact, and on this workload it costs almost nothing:
#: the graphs are small enough that threading barely helps, measured at
#: 2.00 s/epoch on one thread against 1.86 s on eight and 1.53 s on two (eight
#: is slower than two -- thread overhead dominates the small matmuls). Paying
#: ~16% wall time for an exactly reproducible model is worth it, and a
#: one-thread fit also leaves the rest of the machine free.
_GNN_THREADS = 1

#: Each configuration is fitted this many times, with different torch seeds,
#: and scored on the mean validation AP. One fit is not a reliable measurement:
#: `index_add_` over several CPU threads accumulates in a nondeterministic
#: order, and the initialization differs too, so repeated fits of the *same*
#: configuration vary by an amount comparable to the gaps between neighbouring
#: configurations. Selecting the argmax over single fits would partly be
#: selecting that noise. See results/gnn_tuning.csv for the per-trial spread.
_GNN_TUNE_REPEATS = 3

#: The split mode that gets the full search. The other mode is searched over a
#: shortlist of the configurations that did best here, which keeps the cost to
#: roughly a third of a full search while still letting the two modes disagree.
_GNN_PRIMARY_MODE = "scaffold"
_GNN_SHORTLIST_N = 3

#: Worker processes for the search. Because every fit is single-threaded and
#: self-contained, running them concurrently cannot change any result -- it is
#: purely a wall-clock choice, and the chosen configuration is identical at
#: any value. 1 keeps everything in this process; stage 3 raises it.
_GNN_TUNE_JOBS = 1

#: every (model, mode, config) trial and its validation AP, for the record
GNN_TRIALS: list[dict] = []

#: shortlists carried from the primary mode, keyed by model name
_GNN_SHORTLIST: dict[str, list[dict]] = {}


def get_table(kind: str, instances: list[dict], label_field: str) -> AtomTable:
    key = (kind, label_field)
    if key not in _TABLE_CACHE:
        _TABLE_CACHE[key] = build_atom_table(instances, kind, label_field=label_field)
    return _TABLE_CACHE[key]


def clear_caches() -> None:
    _TABLE_CACHE.clear()
    _TUNED.clear()
    _GNN_SHORTLIST.clear()
    GNN_TRIALS.clear()


@dataclass
class Prediction:
    """Predictions on the held-out test atoms, plus what is needed to audit them."""
    model: str
    y_test: np.ndarray
    p_test: np.ndarray
    inst_test: np.ndarray
    atom_test: np.ndarray
    y_val: np.ndarray
    p_val: np.ndarray
    inst_val: np.ndarray
    p_test_calibrated: np.ndarray | None = None
    feature_importance: dict[str, float] = field(default_factory=dict)
    chosen_params: dict = field(default_factory=dict)
    n_train_atoms: int = 0


# --------------------------------------------------------------------------
# trivial frequency reference
# --------------------------------------------------------------------------

def _atom_type_keys(instances: list[dict], table: AtomTable) -> np.ndarray:
    """Coarse atom-type key per row, for the frequency reference."""
    cache: dict[str, list[str]] = {}
    keys = []
    for iid, aidx in zip(table.inst_id, table.atom_idx):
        smi = instances[int(iid)]["smiles"]
        if smi not in cache:
            mol = Chem.MolFromSmiles(smi)
            cache[smi] = [
                f"{a.GetSymbol()}|H{a.GetTotalNumHs()}|"
                f"{int(a.GetIsAromatic())}|{int(a.IsInRing())}"
                for a in mol.GetAtoms()
            ]
        keys.append(cache[smi][int(aidx)])
    return np.asarray(keys, dtype=object)


def _fit_prior(
    instances: list[dict], table: AtomTable, split: Split, label_field: str
) -> Prediction:
    keys = _atom_type_keys(instances, table)
    tr = atom_mask(table, split.train_inst)
    va = atom_mask(table, split.val_inst)
    te = atom_mask(table, split.test_inst)

    base = float(table.y[tr].mean())
    pos: dict[str, float] = {}
    tot: dict[str, float] = {}
    for k, yy in zip(keys[tr], table.y[tr]):
        tot[k] = tot.get(k, 0.0) + 1.0
        pos[k] = pos.get(k, 0.0) + float(yy)

    def lookup(arr):
        out = np.empty(arr.size, dtype=np.float64)
        for i, k in enumerate(arr):
            n = tot.get(k, 0.0)
            # Laplace smoothing towards the global rate; unseen types get it
            out[i] = (pos.get(k, 0.0) + 1.0 * base) / (n + 1.0)
        return out

    return Prediction(
        model="prior",
        y_test=table.y[te], p_test=lookup(keys[te]),
        inst_test=table.inst_id[te], atom_test=table.atom_idx[te],
        y_val=table.y[va], p_val=lookup(keys[va]), inst_val=table.inst_id[va],
        n_train_atoms=int(tr.sum()),
        chosen_params={"n_atom_types": len(tot)},
    )


# --------------------------------------------------------------------------
# gradient boosting on tabular features
# --------------------------------------------------------------------------

def _load_smartcyp() -> dict[tuple[str, int], dict[str, float]]:
    import csv as _csv

    from .paths import INTERIM

    path = INTERIM / "smartcyp.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Build the SMARTCyp jar and run "
            "`scripts/08_smartcyp.py --jar ...` first; see its docstring."
        )
    out: dict[tuple[str, int], dict[str, float]] = {}
    with path.open(newline="") as fh:
        for r in _csv.DictReader(fh):
            out[(r["smiles"], int(r["atom_idx"]))] = {
                "score": float(r["score"]),
                "score_2d6": float(r["score_2d6"]),
                "score_2c": float(r["score_2c"]),
            }
    return out


def _fit_smartcyp(
    instances: list[dict], table: AtomTable, split: Split, label_field: str
) -> Prediction:
    """SMARTCyp as an external baseline.

    Nothing is trained. The only thing fitted is a monotone isotonic map from
    SMARTCyp's score to a probability, fitted on the *training* atoms purely so
    that Brier score and ECE are defined. Because the map is monotone and
    decreasing, it cannot change the ranking, so top-k, ROC-AUC and average
    precision are exactly SMARTCyp's own out-of-the-box numbers.

    Atoms SMARTCyp does not score (oxygens, halogens) get its own "will not
    react" sentinel.
    """
    from .paths import INTERIM  # noqa: F401  (import error surfaces earlier)

    scores_tab = _load_smartcyp()
    raw = np.empty(table.X.shape[0], dtype=np.float64)
    for row, (iid, aidx) in enumerate(zip(table.inst_id, table.atom_idx)):
        inst = instances[int(iid)]
        rec = scores_tab.get((inst["smiles"], int(aidx)))
        if rec is None:
            raw[row] = 999.0
            continue
        col = _SMARTCYP_COLUMN.get(inst["isoform"], "score")
        raw[row] = rec[col]

    tr = atom_mask(table, split.train_inst)
    va = atom_mask(table, split.val_inst)
    te = atom_mask(table, split.test_inst)

    iso = IsotonicRegression(
        increasing=False, out_of_bounds="clip", y_min=0.0, y_max=1.0
    )
    iso.fit(raw[tr], table.y[tr])

    return Prediction(
        model="smartcyp",
        y_test=table.y[te], p_test=iso.predict(raw[te]),
        inst_test=table.inst_id[te], atom_test=table.atom_idx[te],
        y_val=table.y[va], p_val=iso.predict(raw[va]),
        inst_val=table.inst_id[va],
        n_train_atoms=int(tr.sum()),
        chosen_params={"note": "external method; only the score->probability "
                               "map is fitted, on train atoms, and it is "
                               "monotone so rankings are unchanged"},
    )


def _fit_gbm(
    name: str,
    kind: str,
    instances: list[dict],
    split: Split,
    label_field: str,
    seed: int,
    tune: bool,
) -> Prediction:
    table = get_table(kind, instances, label_field)
    tr = atom_mask(table, split.train_inst)
    va = atom_mask(table, split.val_inst)
    te = atom_mask(table, split.test_inst)

    Xtr, ytr = table.X[tr], table.y[tr]
    Xva, yva = table.X[va], table.y[va]

    key = (name, split.mode)
    if tune and key not in _TUNED:
        best, best_score = None, -np.inf
        for params in _GBM_GRID:
            clf = HistGradientBoostingClassifier(
                max_iter=300, early_stopping=False, l2_regularization=1.0,
                min_samples_leaf=20, random_state=0, **params,
            )
            clf.fit(Xtr, ytr)
            score = average_precision_score(yva, clf.predict_proba(Xva)[:, 1])
            if score > best_score:
                best, best_score = params, score
        _TUNED[key] = dict(best)
    params = _TUNED.get(key, _GBM_GRID[0])

    clf = HistGradientBoostingClassifier(
        max_iter=300, early_stopping=False, l2_regularization=1.0,
        min_samples_leaf=20, random_state=seed, **params,
    )
    clf.fit(Xtr, ytr)

    p_te = clf.predict_proba(table.X[te])[:, 1]
    p_va = clf.predict_proba(Xva)[:, 1]

    importance: dict[str, float] = {}
    if kind == "phys":
        # a linear model on the same features, fitted only so that signed
        # coefficients can be reported; it is not used for predictions
        lin = make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=2000, C=1.0, random_state=seed),
        )
        try:
            lin.fit(Xtr, ytr)
            coefs = lin[-1].coef_[0]
            importance = {
                n: float(c) for n, c in zip(table.feature_names, coefs)
            }
        except Exception:  # noqa: BLE001
            importance = {}

    return Prediction(
        model=name,
        y_test=table.y[te], p_test=p_te,
        inst_test=table.inst_id[te], atom_test=table.atom_idx[te],
        y_val=yva, p_val=p_va, inst_val=table.inst_id[va],
        feature_importance=importance,
        chosen_params=dict(params),
        n_train_atoms=int(tr.sum()),
    )


# --------------------------------------------------------------------------
# graph neural network
# --------------------------------------------------------------------------

_BOND_TYPES = [
    Chem.BondType.SINGLE,
    Chem.BondType.DOUBLE,
    Chem.BondType.TRIPLE,
    Chem.BondType.AROMATIC,
]


def _bond_features(bond: Chem.Bond) -> list[float]:
    vec = [0.0] * (len(_BOND_TYPES) + 1)
    try:
        vec[_BOND_TYPES.index(bond.GetBondType())] = 1.0
    except ValueError:
        vec[-1] = 1.0
    return vec + [float(bond.IsInRing()), float(bond.GetIsConjugated())]


def _graph_cache(instances: list[dict], kind: str = "graph"):
    """Per-molecule node features and edge lists, built once."""
    from .features import molecule_features

    cache: dict[str, dict] = {}
    for row in instances:
        smi = row["smiles"]
        if smi in cache:
            continue
        mol = Chem.MolFromSmiles(smi)
        nf, _ = molecule_features(mol, kind, smiles=smi)
        src, dst, ef = [], [], []
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            feat = _bond_features(b)
            src += [i, j]
            dst += [j, i]
            ef += [feat, feat]
        cache[smi] = {
            "x": np.asarray(nf, dtype=np.float32),
            "src": np.asarray(src, dtype=np.int64),
            "dst": np.asarray(dst, dtype=np.int64),
            "edge": np.asarray(ef, dtype=np.float32).reshape(len(ef), -1),
        }
    return cache


def _build_mpnn(
    n_node_feat: int, n_edge_feat: int, n_iso: int, hidden: int, layers: int,
    dropout: float = 0.0,
):
    """Residual message-passing network.

    ``dropout`` is applied to the residual update and inside the read-out head,
    not to the aggregated messages themselves: dropping messages would make the
    effective neighbourhood of an atom random, which is a different (and for
    this task worse) kind of regularization than dropping hidden units.
    """
    import torch
    from torch import nn

    class Layer(nn.Module):
        def __init__(self):
            super().__init__()
            self.msg = nn.Sequential(
                nn.Linear(2 * hidden + n_edge_feat, hidden), nn.ReLU(),
                nn.Linear(hidden, hidden),
            )
            self.upd = nn.Sequential(
                nn.Linear(2 * hidden, hidden), nn.ReLU(),
                nn.Dropout(dropout), nn.Linear(hidden, hidden),
            )
            self.norm = nn.LayerNorm(hidden)

        def forward(self, h, src, dst, edge):
            m = self.msg(torch.cat([h[src], h[dst], edge], dim=1))
            agg = torch.zeros_like(h).index_add_(0, dst, m)
            return self.norm(h + self.upd(torch.cat([h, agg], dim=1)))

    class MPNN(nn.Module):
        def __init__(self):
            super().__init__()
            self.embed = nn.Sequential(nn.Linear(n_node_feat, hidden), nn.ReLU())
            self.iso = nn.Embedding(n_iso, hidden)
            self.layers = nn.ModuleList(Layer() for _ in range(layers))
            self.head = nn.Sequential(
                nn.Linear(hidden, hidden), nn.ReLU(),
                nn.Dropout(dropout), nn.Linear(hidden, 1),
            )

        def forward(self, x, src, dst, edge, iso_idx):
            h = self.embed(x) + self.iso(iso_idx)
            for layer in self.layers:
                h = layer(h, src, dst, edge)
            return self.head(h).squeeze(-1)

    return MPNN()


def _collate(batch, cache, instances, iso_index):
    """Concatenate molecules into one disconnected graph."""
    import torch

    xs, srcs, dsts, efs, isos, ys, iids, aidx = [], [], [], [], [], [], [], []
    offset = 0
    for iid in batch:
        row = instances[iid]
        g = cache[row["smiles"]]
        n = g["x"].shape[0]
        xs.append(g["x"])
        srcs.append(g["src"] + offset)
        dsts.append(g["dst"] + offset)
        efs.append(g["edge"])
        isos.append(np.full(n, iso_index[row["isoform"]], dtype=np.int64))
        lab = np.zeros(n, dtype=np.float32)
        lab[np.asarray(row["_labels"], dtype=int)] = 1.0
        ys.append(lab)
        iids.append(np.full(n, iid, dtype=np.int64))
        aidx.append(np.arange(n, dtype=np.int64))
        offset += n

    edge = np.vstack(efs) if efs else np.zeros((0, 1), dtype=np.float32)
    return {
        "x": torch.from_numpy(np.vstack(xs)),
        "src": torch.from_numpy(np.concatenate(srcs)),
        "dst": torch.from_numpy(np.concatenate(dsts)),
        "edge": torch.from_numpy(edge),
        "iso": torch.from_numpy(np.concatenate(isos)),
        "y": torch.from_numpy(np.concatenate(ys)),
        "inst": np.concatenate(iids),
        "atom": np.concatenate(aidx),
    }


def _fit_gnn(
    instances: list[dict],
    split: Split,
    label_field: str,
    seed: int,
    node_kind: str = "graph",
    name: str = "gnn",
    hidden: int = 128,
    layers: int = 4,
    epochs: int = 120,
    batch_size: int = 32,
    lr: float = 1e-3,
    patience: int = 20,
    dropout: float = 0.0,
    weight_decay: float = 1e-5,
    pos_weight: float = 1.0,
) -> Prediction:
    import torch
    from torch import nn

    from . import ISOFORMS

    torch.manual_seed(seed)
    np.random.seed(seed)
    torch.set_num_threads(_GNN_THREADS)

    work = [dict(r) for r in instances]
    for r in work:
        r["_labels"] = r[label_field]

    cache = _graph_cache(work, kind=node_kind)
    iso_index = {iso: i for i, iso in enumerate(ISOFORMS)}
    n_node = next(iter(cache.values()))["x"].shape[1]
    n_edge = next(iter(cache.values()))["edge"].shape[1]

    model = _build_mpnn(n_node, n_edge, len(ISOFORMS), hidden, layers, dropout)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=8)
    # pos_weight > 1 up-weights the ~10% of atoms that are sites. It changes the
    # probability scale (and so the Brier score and ECE), but recalibration on
    # the validation split is monotone and absorbs most of that, so it mainly
    # affects where the model spends its capacity.
    lossf = nn.BCEWithLogitsLoss(
        pos_weight=torch.tensor(float(pos_weight)) if pos_weight != 1.0 else None
    )

    # standardize node features using TRAIN molecules only
    train_mols = {work[i]["smiles"] for i in split.train_inst}
    stack = np.vstack([cache[s]["x"] for s in sorted(train_mols)])
    mu = stack.mean(axis=0, keepdims=True)
    sd = stack.std(axis=0, keepdims=True)
    sd[sd < 1e-8] = 1.0
    for g in cache.values():
        g["x"] = ((g["x"] - mu) / sd).astype(np.float32)

    train_ids = list(map(int, split.train_inst))
    val_batch = _collate(list(map(int, split.val_inst)), cache, work, iso_index)
    test_batch = _collate(list(map(int, split.test_inst)), cache, work, iso_index)

    rng = np.random.default_rng(seed)
    best_score, best_state, since = -np.inf, None, 0

    epochs_run = 0
    for _ in range(epochs):
        epochs_run += 1
        model.train()
        order = rng.permutation(len(train_ids))
        for s in range(0, len(order), batch_size):
            ids = [train_ids[j] for j in order[s : s + batch_size]]
            b = _collate(ids, cache, work, iso_index)
            opt.zero_grad()
            logits = model(b["x"], b["src"], b["dst"], b["edge"], b["iso"])
            loss = lossf(logits, b["y"])
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()

        model.eval()
        with torch.no_grad():
            pv = torch.sigmoid(
                model(val_batch["x"], val_batch["src"], val_batch["dst"],
                      val_batch["edge"], val_batch["iso"])
            ).numpy()
        yv = val_batch["y"].numpy()
        score = (
            average_precision_score(yv, pv) if 0 < yv.sum() < yv.size else -np.inf
        )
        sched.step(score)
        if score > best_score:
            best_score = score
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            since = 0
        else:
            since += 1
            if since >= patience:
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        p_te = torch.sigmoid(
            model(test_batch["x"], test_batch["src"], test_batch["dst"],
                  test_batch["edge"], test_batch["iso"])
        ).numpy()
        p_va = torch.sigmoid(
            model(val_batch["x"], val_batch["src"], val_batch["dst"],
                  val_batch["edge"], val_batch["iso"])
        ).numpy()

    return Prediction(
        model=name,
        y_test=test_batch["y"].numpy().astype(np.int8),
        p_test=p_te,
        inst_test=test_batch["inst"],
        atom_test=test_batch["atom"],
        y_val=val_batch["y"].numpy().astype(np.int8),
        p_val=p_va,
        inst_val=val_batch["inst"],
        chosen_params={
            "hidden": hidden, "layers": layers, "lr": lr,
            "batch_size": batch_size, "dropout": dropout,
            "weight_decay": weight_decay, "pos_weight": pos_weight,
            "node_kind": node_kind,
            "n_node_features": int(n_node),
            "epochs_run": int(epochs_run),
            "best_val_ap": float(best_score),
        },
        n_train_atoms=int(sum(cache[work[i]["smiles"]]["x"].shape[0]
                              for i in split.train_inst)),
    )


def _worker_init() -> None:
    """Keep a worker to one thread, including inside the BLAS libraries.

    `_fit_gnn` sets the torch thread count itself, but BLAS reads its
    environment at import, which has already happened by then.
    """
    import os

    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[var] = "1"


def _gnn_trial(task):
    """One fit, as a picklable top-level function for the process pool."""
    instances, split, label_field, seed, node_kind, name, params = task
    pred = _fit_gnn(
        instances, split, label_field, seed,
        node_kind=node_kind, name=name, **params,
    )
    return (
        float(pred.chosen_params["best_val_ap"]),
        int(pred.chosen_params["epochs_run"]),
    )


def _run_trials(tasks: list) -> list:
    """Fit every (configuration, repeat) task, in parallel if asked.

    The caller's entry point must be guarded by ``if __name__ ==
    "__main__":``. That is the documented requirement of the ``spawn`` start
    method, which is the default on macOS: each worker re-imports the main
    module, so an unguarded script re-runs itself in every worker. Stage 3 is
    guarded (and ``test_stage_3_is_guarded_for_spawn`` keeps it that way).

    Falls back to running in this process if the pool cannot be created, since
    a search that is merely slow is much better than one that does not run.
    The results are identical either way -- each fit is single-threaded and
    independent of the others.
    """
    if _GNN_TUNE_JOBS <= 1 or len(tasks) < 2:
        return [_gnn_trial(t) for t in tasks]
    from concurrent.futures import ProcessPoolExecutor
    try:
        with ProcessPoolExecutor(
            max_workers=_GNN_TUNE_JOBS, initializer=_worker_init
        ) as ex:
            return list(ex.map(_gnn_trial, tasks))
    except Exception as exc:  # noqa: BLE001
        first = str(exc).strip().splitlines()[0] if str(exc).strip() else exc
        print(f"    parallel search unavailable, running serially: {first}")
        print("    (if this mentions bootstrapping, the calling script is "
              "missing its `if __name__ == \"__main__\":` guard)")
        return [_gnn_trial(t) for t in tasks]


def _tune_gnn(
    name: str,
    node_kind: str,
    instances: list[dict],
    split: Split,
    label_field: str,
) -> dict:
    """Pick GNN hyperparameters on the validation split, once per (model, mode).

    Selection is on validation average precision only. Every trial does compute
    test predictions, because the trial reuses ``_fit_gnn`` unchanged rather
    than a second near-copy of the training loop, but those predictions are
    discarded here and never enter the choice. The returned configuration is
    then refitted per seed in the normal way.

    Each GNN variant gets its own full search in the primary mode, so the
    ablations (``gnn_huckel`` against ``gnn``, ``gnn_bde`` against ``gnn``) are
    comparisons between equally-tuned models rather than between a tuned model
    and a default one.
    """
    key = (name, split.mode)
    if key in _TUNED:
        return _TUNED[key]

    if split.mode == _GNN_PRIMARY_MODE or name not in _GNN_SHORTLIST:
        grid = _GNN_GRID
        scope = "full"
    else:
        grid = _GNN_SHORTLIST[name]
        scope = "shortlist"

    # the same set of torch seeds for every configuration, so that the
    # comparison is between configurations and not between seed draws
    tasks = [
        (instances, split, label_field, rep, node_kind, name, params)
        for params in grid
        for rep in range(_GNN_TUNE_REPEATS)
    ]
    print(f"    tune {name} {split.mode} [{scope}]: {len(grid)} "
          f"configurations x {_GNN_TUNE_REPEATS} fits"
          + (f", {_GNN_TUNE_JOBS} workers" if _GNN_TUNE_JOBS > 1 else ""))
    results = _run_trials(tasks)

    best, best_score = None, -np.inf
    for i, params in enumerate(grid):
        chunk = results[i * _GNN_TUNE_REPEATS : (i + 1) * _GNN_TUNE_REPEATS]
        scores = [c[0] for c in chunk]
        epochs = [c[1] for c in chunk]
        score = float(np.mean(scores))
        sd = float(np.std(scores, ddof=1)) if len(scores) > 1 else 0.0
        GNN_TRIALS.append({
            "model": name, "mode": split.mode, "scope": scope,
            "trial": i, "params": dict(params),
            "repeats": len(scores),
            "epochs_run": int(np.mean(epochs)),
            "val_avg_precision": score,
            "val_avg_precision_sd": sd,
            "val_avg_precision_runs": scores,
        })
        print(f"      {i + 1:2}/{len(grid)} {str(params or 'base'):58} "
              f"val AP {score:.4f} +- {sd:.4f}")
        if score > best_score:
            best, best_score = dict(params), score

    if split.mode == _GNN_PRIMARY_MODE:
        ranked = sorted(
            (t for t in GNN_TRIALS if t["model"] == name
             and t["mode"] == split.mode),
            key=lambda t: -t["val_avg_precision"],
        )
        _GNN_SHORTLIST[name] = [
            dict(t["params"]) for t in ranked[:_GNN_SHORTLIST_N]
        ]

    _TUNED[key] = best
    print(f"    tune {name} {split.mode}: chose {best or 'base'} "
          f"(val AP {best_score:.4f})")
    return best


# --------------------------------------------------------------------------
# dispatch and calibration
# --------------------------------------------------------------------------

def calibrate(pred: Prediction) -> Prediction:
    """Fit isotonic regression on the validation split and apply it to test.

    Reported alongside the raw probabilities rather than instead of them: the
    raw numbers say how well calibrated the model is out of the box, which is
    the honest thing to show, and the calibrated numbers say how much of the
    miscalibration a held-out recalibration step can remove. Monotone
    recalibration cannot change the ranking, so top-k and AUC are unaffected.
    """
    if pred.y_val.size == 0 or not (0 < pred.y_val.sum() < pred.y_val.size):
        return pred
    iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
    iso.fit(pred.p_val, pred.y_val)
    pred.p_test_calibrated = iso.predict(pred.p_test)
    return pred


def fit_predict(
    name: str,
    instances: list[dict],
    split: Split,
    label_field: str = "som_all",
    seed: int = 0,
    tune: bool = True,
) -> Prediction:
    if name == "prior":
        table = get_table("phys", instances, label_field)
        pred = _fit_prior(instances, table, split, label_field)
    elif name == "smartcyp":
        table = get_table("phys", instances, label_field)
        pred = _fit_smartcyp(instances, table, split, label_field)
    elif name == "fp_gbm":
        pred = _fit_gbm("fp_gbm", "fp", instances, split, label_field, seed, tune)
    elif name == "phys_gbm":
        pred = _fit_gbm("phys_gbm", "phys", instances, split, label_field, seed, tune)
    elif name == "phys_bde_gbm":
        pred = _fit_gbm("phys_bde_gbm", "phys_bde", instances, split,
                        label_field, seed, tune)
    elif name in ("gnn", "gnn_huckel", "gnn_bde"):
        node_kind = {"gnn": "graph", "gnn_huckel": "phys",
                     "gnn_bde": "graph_bde"}[name]
        params = (
            _tune_gnn(name, node_kind, instances, split, label_field)
            if tune else _TUNED.get((name, split.mode), {})
        )
        pred = _fit_gnn(instances, split, label_field, seed,
                        node_kind=node_kind, name=name, **(params or {}))
    else:
        raise ValueError(f"unknown model: {name!r}")
    return calibrate(pred)
