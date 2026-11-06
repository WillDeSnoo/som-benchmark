"""Canonical locations for data, results and figures.

Everything is resolved relative to the repository root so that the stages can be
run from any working directory.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

DATA = ROOT / "data"
RAW = DATA / "raw"
INTERIM = DATA / "interim"
PROCESSED = DATA / "processed"

RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"


def ensure_dirs() -> None:
    for d in (RAW, INTERIM, PROCESSED, RESULTS, FIGURES):
        d.mkdir(parents=True, exist_ok=True)
