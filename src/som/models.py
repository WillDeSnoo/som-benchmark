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
    quietly containing it.

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

MODELS = ("prior", "fp_gbm", "phys_gbm", "gnn", "gnn_huckel")

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


def get_table(kind: str, instances: list[dict], label_field: str) -> AtomTable:
    key = (kind, label_field)
    if key not in _TABLE_CACHE:
        _TABLE_CACHE[key] = build_atom_table(instances, kind, label_field=label_field)
    return _TABLE_CACHE[key]


def clear_caches() -> None:
    _TABLE_CACHE.clear()
    _TUNED.clear()


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
        nf, _ = molecule_features(mol, kind)
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


def _build_mpnn(n_node_feat: int, n_edge_feat: int, n_iso: int, hidden: int, layers: int):
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
                nn.Linear(2 * hidden, hidden), nn.ReLU(), nn.Linear(hidden, hidden)
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
                nn.Linear(hidden, hidden), nn.ReLU(), nn.Linear(hidden, 1)
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
) -> Prediction:
    import torch
    from torch import nn

    from . import ISOFORMS

    torch.manual_seed(seed)
    np.random.seed(seed)
    torch.set_num_threads(max(1, min(8, torch.get_num_threads())))

    work = [dict(r) for r in instances]
    for r in work:
        r["_labels"] = r[label_field]

    cache = _graph_cache(work, kind=node_kind)
    iso_index = {iso: i for i, iso in enumerate(ISOFORMS)}
    n_node = next(iter(cache.values()))["x"].shape[1]
    n_edge = next(iter(cache.values()))["edge"].shape[1]

    model = _build_mpnn(n_node, n_edge, len(ISOFORMS), hidden, layers)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-5)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=8)
    lossf = nn.BCEWithLogitsLoss()

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

    for _ in range(epochs):
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
            "batch_size": batch_size, "node_kind": node_kind,
            "n_node_features": int(n_node),
            "best_val_ap": float(best_score),
        },
        n_train_atoms=int(sum(cache[work[i]["smiles"]]["x"].shape[0]
                              for i in split.train_inst)),
    )


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
    elif name == "fp_gbm":
        pred = _fit_gbm("fp_gbm", "fp", instances, split, label_field, seed, tune)
    elif name == "phys_gbm":
        pred = _fit_gbm("phys_gbm", "phys", instances, split, label_field, seed, tune)
    elif name == "gnn":
        pred = _fit_gnn(instances, split, label_field, seed,
                        node_kind="graph", name="gnn")
    elif name == "gnn_huckel":
        pred = _fit_gnn(instances, split, label_field, seed,
                        node_kind="phys", name="gnn_huckel")
    else:
        raise ValueError(f"unknown model: {name!r}")
    return calibrate(pred)
