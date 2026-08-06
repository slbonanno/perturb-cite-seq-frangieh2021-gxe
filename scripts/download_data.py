#!/usr/bin/env python
"""Fetch the raw h5ads into data/raw/.

    python scripts/download_data.py [--overwrite]

Roughly 3-5 GB total. Idempotent: existing files are skipped.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.config import load_config           # noqa: E402
from src.data_io import fetch_raw            # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    cfg = load_config()
    files = fetch_raw(cfg, overwrite=args.overwrite)
    print("\nReady:")
    for k, v in files.items():
        print(f"  {k:8s} {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
