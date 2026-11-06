#!/usr/bin/env python
"""Stage 1: fetch the raw per-isoform SDFs from the pinned upstream commit.

Nothing is redistributed by this repository; see DATA.md.
"""

import argparse

from som import UPSTREAM_COMMIT, UPSTREAM_REPO
from som.data import download
from som.paths import RAW


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true", help="re-download existing files")
    args = ap.parse_args()

    print(f"source : github.com/{UPSTREAM_REPO}")
    print(f"commit : {UPSTREAM_COMMIT}")
    print(f"dest   : {RAW}\n")
    download(force=args.force)
    print("\nSee DATA.md for provenance and licence. Cite the primary sources.")


if __name__ == "__main__":
    main()
