#!/usr/bin/env python
"""Stage 2: standardize, relabel, deduplicate; write the curated tables.

Prints a curation report. Every excluded record is written to
data/processed/exclusions.csv with the reason, so nothing disappears silently.
"""

import argparse
import json
from collections import Counter

from som import ISOFORMS
from som.data import curate, molecule_table, write_curated
from som.paths import RESULTS


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--primary-only", action="store_true",
        help="treat only PRIMARY_SOM as positive (default: primary+secondary+tertiary)",
    )
    args = ap.parse_args()

    fields = ("PRIMARY_SOM",) if args.primary_only else (
        "PRIMARY_SOM", "SECONDARY_SOM", "TERTIARY_SOM"
    )

    rows, report = curate(som_fields=fields)
    paths = write_curated(rows, report)
    mols = molecule_table(rows)

    n_atoms = sum(r["n_heavy_atoms"] for r in rows)
    n_pos = sum(len(r["som_all"]) for r in rows)
    n_pos_primary = sum(len(r["som_primary"]) for r in rows)
    scaffolds = Counter(m["scaffold"] for m in mols)
    acyclic = scaffolds.get("", 0)

    print("=" * 66)
    print("CURATION REPORT")
    print("=" * 66)
    print(f"records read from SDFs           {report.records_read:>8}")
    print(f"excluded                         {len(report.exclusions):>8}")
    for reason, n in sorted(Counter(e.reason for e in report.exclusions).items()):
        print(f"    {reason:<32} {n:>6}")
    print(f"duplicate instances merged       {report.duplicate_instances_merged:>8}")
    print(f"  of which labels disagreed      {report.duplicate_label_conflicts:>8}")
    print("-" * 66)
    print(f"surviving (molecule, isoform)    {len(rows):>8}")
    print(f"unique molecules                 {len(mols):>8}")
    print(f"unique Bemis-Murcko scaffolds    {len(scaffolds) - (1 if acyclic else 0):>8}"
          f"   (+{acyclic} acyclic)")
    print(f"candidate atoms (all heavy)      {n_atoms:>8}")
    print(f"positive atoms (all tiers)       {n_pos:>8}"
          f"   = {100 * n_pos / n_atoms:.2f}% of candidates")
    print(f"positive atoms (primary only)    {n_pos_primary:>8}"
          f"   = {100 * n_pos_primary / n_atoms:.2f}% of candidates")
    print("-" * 66)
    print(f"{'isoform':<9}{'instances':>10}{'atoms':>9}{'positives':>11}{'% pos':>8}")
    per_iso = {}
    for iso in ISOFORMS:
        sub = [r for r in rows if r["isoform"] == iso]
        a = sum(r["n_heavy_atoms"] for r in sub)
        p = sum(len(r["som_all"]) for r in sub)
        per_iso[iso] = {"instances": len(sub), "atoms": a, "positives": p}
        print(f"{iso:<9}{len(sub):>10}{a:>9}{p:>11}{100 * p / a:>7.2f}%")
    print("=" * 66)
    for name, p in paths.items():
        print(f"wrote {name:<11} {p}")

    RESULTS.mkdir(parents=True, exist_ok=True)
    summary = {
        "som_fields": list(fields),
        "records_read": report.records_read,
        "excluded_total": len(report.exclusions),
        "excluded_by_reason": dict(Counter(e.reason for e in report.exclusions)),
        "duplicate_instances_merged": report.duplicate_instances_merged,
        "duplicate_label_conflicts": report.duplicate_label_conflicts,
        "instances": len(rows),
        "molecules": len(mols),
        "scaffolds": len(scaffolds) - (1 if acyclic else 0),
        "acyclic_molecules": acyclic,
        "candidate_atoms": n_atoms,
        "positive_atoms_all_tiers": n_pos,
        "positive_atoms_primary_only": n_pos_primary,
        "per_isoform": per_iso,
    }
    out = RESULTS / "curation_summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"wrote summary     {out}")


if __name__ == "__main__":
    main()
