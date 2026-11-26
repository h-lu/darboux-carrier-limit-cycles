#!/usr/bin/env python3
"""Recompute exact ranks and integer tables, then compare the reference files."""
import argparse
import json
from pathlib import Path
import sys

import certificates
import ranks


def compare(path, actual):
    expected = path.read_text(encoding="utf-8")
    if actual != expected:
        raise ValueError(f"{path.name} differs from the recomputed result")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--part", choices=("all", "ranks", "certificates"), default="all")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--reference", type=Path, default=Path(__file__).with_name("results"),
                      help="compare with files in this directory (default: supplied results)")
    mode.add_argument("--write", type=Path, help="write recomputed files in this directory")
    args = parser.parse_args()
    products = {}
    def json_text(value):
        return json.dumps(value, indent=2, sort_keys=True) + "\n"
    if args.part in ("all", "ranks"):
        products["ranks.json"] = json_text(ranks.run())
    if args.part in ("all", "certificates"):
        summary, table = certificates.run()
        products["certificates.json"] = json_text(summary)
        products["degrees.csv"] = table
    for name, value in products.items():
        if args.write:
            args.write.mkdir(exist_ok=True, parents=True)
            (args.write/name).write_text(value, encoding="utf-8")
        else:
            compare(args.reference/name, value)
        print(f"{name} OK")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, ZeroDivisionError) as exc:
        sys.exit(f"Verification failed: {exc}")
