#!/usr/bin/env python
"""Stage 4: failure analysis from the saved predictions. Trains nothing.

Writes:
  results/per_isoform.csv     top-1/top-2 per CYP isoform
  results/per_site_type.csv   recall at top-2 per implied transformation
  results/novelty_curve.csv   accuracy vs similarity of the test molecule to train
  results/split_gap.csv       random minus scaffold, per model, with an interval
  results/calibration.csv     reliability curves
"""

import argparse
import csv
from collections import defaultdict

import numpy as np

from som.analysis import (
    atom_ranks,
    bucket,
    instance_topk_table,
    site_types_for,
)
from som.data import load_instances, load_molecules
from som.metrics import expected_calibration_error, mean_ci_across_seeds
from som.models import MODELS
from som.paths import INTERIM, RESULTS
from som.splits import SPLIT_MODES, make_split, nearest_neighbour_similarity

SIM_EDGES = [0.3, 0.45, 0.6, 0.75]
SIM_LABELS = ["<=0.30", "0.30-0.45", "0.45-0.60", "0.60-0.75", ">0.75"]

ATOM_EDGES = [20.0, 30.0, 40.0]
ATOM_LABELS = ["<=20", "21-30", "31-40", ">40"]


def write_csv(path, rows, fields=None):
    if not rows:
        print(f"  (nothing to write for {path.name})")
        return
    fields = fields or list(rows[0].keys())
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"  wrote {path}")


def load_pred(model, mode, seed):
    path = INTERIM / "predictions" / f"{model}_{mode}_seed{seed}.npz"
    if not path.exists():
        return None
    d = np.load(path)
    return {k: d[k] for k in d.files}


def agg(values):
    s = mean_ci_across_seeds(values)
    return s["mean"], s["lo"], s["hi"], s["n"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--models", nargs="+", default=list(MODELS))
    ap.add_argument("--modes", nargs="+", default=list(SPLIT_MODES))
    ap.add_argument("--label-field", default="som_all",
                    choices=["som_all", "som_primary"])
    args = ap.parse_args()

    instances = load_instances()
    molecules = load_molecules()
    RESULTS.mkdir(parents=True, exist_ok=True)

    # site-type lookup per (instance, atom)
    site_type_of: dict[tuple[int, int], str] = {}
    for iid, row in enumerate(instances):
        for a, t in site_types_for(row["smiles"], row[args.label_field]).items():
            site_type_of[(iid, a)] = t

    # nearest-neighbour similarity of each test molecule, per split
    nn_sim: dict[tuple[str, int], dict[str, float]] = {}
    for mode in args.modes:
        for seed in range(args.seeds):
            split = make_split(instances, molecules, mode, seed)
            sims = nearest_neighbour_similarity(split.test_mols, split.train_mols)
            nn_sim[(mode, seed)] = dict(zip(split.test_mols, sims))

    iso_rows, site_rows, nov_rows, cal_rows, atom_rows = [], [], [], [], []
    gap_source: dict[tuple[str, str, str], list[float]] = defaultdict(list)

    for model in args.models:
        for mode in args.modes:
            per_iso: dict[str, dict[str, list[float]]] = defaultdict(
                lambda: defaultdict(list))
            per_site: dict[str, list[tuple[int, int]]] = defaultdict(list)
            per_nov: dict[str, dict[str, list[float]]] = defaultdict(
                lambda: defaultdict(list))
            per_atoms: dict[str, dict[str, list[float]]] = defaultdict(
                lambda: defaultdict(list))
            curves: list[list[dict]] = []

            found = 0
            for seed in range(args.seeds):
                d = load_pred(model, mode, seed)
                if d is None:
                    continue
                found += 1
                y, p, inst, atom = d["y"], d["p"], d["inst"], d["atom"]

                # ---- per-instance outcomes ----
                table = instance_topk_table(y, p, inst, ks=(1, 2), seed=seed)
                iso_seed: dict[str, list[bool]] = defaultdict(list)
                iso_seed2: dict[str, list[bool]] = defaultdict(list)
                nov_seed: dict[str, list[bool]] = defaultdict(list)
                nov_seed2: dict[str, list[bool]] = defaultdict(list)
                at_seed: dict[str, list[bool]] = defaultdict(list)
                at_seed2: dict[str, list[bool]] = defaultdict(list)

                for iid, rec in table.items():
                    row = instances[iid]
                    iso_seed[row["isoform"]].append(rec["top1"])
                    iso_seed2[row["isoform"]].append(rec["top2"])
                    sim = nn_sim[(mode, seed)].get(row["smiles"], float("nan"))
                    if np.isfinite(sim):
                        b = bucket(sim, SIM_EDGES, SIM_LABELS)
                        nov_seed[b].append(rec["top1"])
                        nov_seed2[b].append(rec["top2"])
                    ab = bucket(float(row["n_heavy_atoms"]), ATOM_EDGES, ATOM_LABELS)
                    at_seed[ab].append(rec["top1"])
                    at_seed2[ab].append(rec["top2"])

                for k, v in iso_seed.items():
                    per_iso[k]["top1"].append(float(np.mean(v)))
                    per_iso[k]["n"].append(len(v))
                for k, v in iso_seed2.items():
                    per_iso[k]["top2"].append(float(np.mean(v)))
                for k, v in nov_seed.items():
                    per_nov[k]["top1"].append(float(np.mean(v)))
                    per_nov[k]["n"].append(len(v))
                for k, v in nov_seed2.items():
                    per_nov[k]["top2"].append(float(np.mean(v)))
                for k, v in at_seed.items():
                    per_atoms[k]["top1"].append(float(np.mean(v)))
                    per_atoms[k]["n"].append(len(v))
                for k, v in at_seed2.items():
                    per_atoms[k]["top2"].append(float(np.mean(v)))

                # ---- per-positive-atom recall at top-2, by site type ----
                ranks = atom_ranks(p, inst, seed=seed)
                for j in np.flatnonzero(y == 1):
                    t = site_type_of.get((int(inst[j]), int(atom[j])), "other")
                    per_site[t].append((1 if ranks[j] < 2 else 0, 1))

                # ---- calibration curve ----
                _, curve = expected_calibration_error(y, np.clip(p, 0, 1))
                curves.append(curve)

                for metric, val in (
                    ("top1", float(np.mean([r["top1"] for r in table.values()]))),
                    ("top2", float(np.mean([r["top2"] for r in table.values()]))),
                ):
                    gap_source[(model, metric, mode)].append(val)

            if found == 0:
                continue

            for iso, vals in sorted(per_iso.items()):
                m1, l1, h1, n = agg(vals["top1"])
                m2, l2, h2, _ = agg(vals["top2"])
                iso_rows.append({
                    "model": model, "mode": mode, "isoform": iso,
                    "mean_test_instances": round(float(np.mean(vals["n"])), 1),
                    "n_seeds": n,
                    "top1_mean": round(m1, 4), "top1_lo": round(l1, 4),
                    "top1_hi": round(h1, 4),
                    "top2_mean": round(m2, 4), "top2_lo": round(l2, 4),
                    "top2_hi": round(h2, 4),
                })

            for t, pairs in sorted(per_site.items(), key=lambda kv: -len(kv[1])):
                hits = sum(h for h, _ in pairs)
                n = len(pairs)
                site_rows.append({
                    "model": model, "mode": mode, "site_type": t,
                    "n_positive_atoms_pooled": n,
                    "recall_at_2_pooled": round(hits / n, 4) if n else "",
                })

            for b in SIM_LABELS:
                if b not in per_nov:
                    continue
                m1, l1, h1, n = agg(per_nov[b]["top1"])
                m2, l2, h2, _ = agg(per_nov[b]["top2"])
                nov_rows.append({
                    "model": model, "mode": mode, "nn_tanimoto_bucket": b,
                    "mean_test_instances": round(float(np.mean(per_nov[b]["n"])), 1),
                    "n_seeds": n,
                    "top1_mean": round(m1, 4), "top1_lo": round(l1, 4),
                    "top1_hi": round(h1, 4),
                    "top2_mean": round(m2, 4), "top2_lo": round(l2, 4),
                    "top2_hi": round(h2, 4),
                })

            for b in ATOM_LABELS:
                if b not in per_atoms:
                    continue
                m1, l1, h1, n = agg(per_atoms[b]["top1"])
                m2, _, _, _ = agg(per_atoms[b]["top2"])
                atom_rows.append({
                    "model": model, "mode": mode, "heavy_atom_bucket": b,
                    "mean_test_instances": round(float(np.mean(per_atoms[b]["n"])), 1),
                    "top1_mean": round(m1, 4), "top1_lo": round(l1, 4),
                    "top1_hi": round(h1, 4), "top2_mean": round(m2, 4),
                })

            # average the reliability curves bin-by-bin across seeds
            if curves:
                n_bins = min(len(c) for c in curves)
                for b in range(n_bins):
                    cal_rows.append({
                        "model": model, "mode": mode, "bin": b,
                        "mean_predicted": round(
                            float(np.mean([c[b]["mean_predicted"] for c in curves])), 5),
                        "observed_frequency": round(
                            float(np.mean([c[b]["observed_frequency"] for c in curves])), 5),
                        "mean_n": round(
                            float(np.mean([c[b]["n"] for c in curves])), 1),
                    })

    # ---- pooled novelty curve ----
    # Both split families are pooled per seed here. Each on its own covers only
    # part of the similarity range (the random split rarely produces a test
    # molecule with no close training neighbour, the scaffold split rarely
    # produces one with a very close neighbour), so pooling them is what turns
    # two disconnected sets of points into a single decay curve.
    for model in args.models:
        per_nov: dict[str, dict[str, list[float]]] = defaultdict(
            lambda: defaultdict(list))
        for seed in range(args.seeds):
            h1: dict[str, list[bool]] = defaultdict(list)
            h2: dict[str, list[bool]] = defaultdict(list)
            for mode in args.modes:
                d = load_pred(model, mode, seed)
                if d is None:
                    continue
                table = instance_topk_table(
                    d["y"], d["p"], d["inst"], ks=(1, 2), seed=seed
                )
                for iid, rec in table.items():
                    row = instances[iid]
                    sim = nn_sim[(mode, seed)].get(row["smiles"], float("nan"))
                    if not np.isfinite(sim):
                        continue
                    b = bucket(sim, SIM_EDGES, SIM_LABELS)
                    h1[b].append(rec["top1"])
                    h2[b].append(rec["top2"])
            for k, v in h1.items():
                per_nov[k]["top1"].append(float(np.mean(v)))
                per_nov[k]["n"].append(len(v))
            for k, v in h2.items():
                per_nov[k]["top2"].append(float(np.mean(v)))

        for b in SIM_LABELS:
            if b not in per_nov:
                continue
            m1, l1, h1v, n = agg(per_nov[b]["top1"])
            m2, l2, h2v, _ = agg(per_nov[b]["top2"])
            nov_rows.append({
                "model": model, "mode": "pooled", "nn_tanimoto_bucket": b,
                "mean_test_instances": round(float(np.mean(per_nov[b]["n"])), 1),
                "n_seeds": n,
                "top1_mean": round(m1, 4), "top1_lo": round(l1, 4),
                "top1_hi": round(h1v, 4),
                "top2_mean": round(m2, 4), "top2_lo": round(l2, 4),
                "top2_hi": round(h2v, 4),
            })

    # ---- random vs scaffold gap ----
    gap_rows = []
    for model in args.models:
        for metric in ("top1", "top2"):
            r = gap_source.get((model, metric, "random"), [])
            s = gap_source.get((model, metric, "scaffold"), [])
            if not r or not s:
                continue
            ra, rl, rh, _ = agg(r)
            sa, sl, sh, _ = agg(s)
            # unpaired difference of means: the two split families are drawn
            # independently, so the seed index does not pair them
            na, nb = len(r), len(s)
            va = np.var(r, ddof=1) if na > 1 else 0.0
            vb = np.var(s, ddof=1) if nb > 1 else 0.0
            se = float(np.sqrt(va / na + vb / nb))
            gap_rows.append({
                "model": model, "metric": metric,
                "random_mean": round(ra, 4), "random_lo": round(rl, 4),
                "random_hi": round(rh, 4),
                "scaffold_mean": round(sa, 4), "scaffold_lo": round(sl, 4),
                "scaffold_hi": round(sh, 4),
                "gap": round(ra - sa, 4),
                "gap_lo": round(ra - sa - 1.96 * se, 4),
                "gap_hi": round(ra - sa + 1.96 * se, 4),
            })

    print("writing analysis tables:")
    write_csv(RESULTS / "per_isoform.csv", iso_rows)
    write_csv(RESULTS / "per_site_type.csv", site_rows)
    write_csv(RESULTS / "novelty_curve.csv", nov_rows)
    write_csv(RESULTS / "per_molecule_size.csv", atom_rows)
    write_csv(RESULTS / "split_gap.csv", gap_rows)
    write_csv(RESULTS / "calibration.csv", cal_rows)

    if gap_rows:
        print("\nrandom minus scaffold (top-1):")
        for r in gap_rows:
            if r["metric"] == "top1":
                print(f"  {r['model']:<10} random {r['random_mean']:.3f}  "
                      f"scaffold {r['scaffold_mean']:.3f}  "
                      f"gap {r['gap']:+.3f} "
                      f"[{r['gap_lo']:+.3f},{r['gap_hi']:+.3f}]")


if __name__ == "__main__":
    main()
