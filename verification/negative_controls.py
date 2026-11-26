#!/usr/bin/env python3
"""Run genuine rank and degree-certificate failures under the invoking Python mode."""
import csv
import io
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def expect_failure(args, text):
    flags = ["-B"] + ["-O"] * sys.flags.optimize
    proc = subprocess.run([sys.executable, *flags, *args], cwd=HERE,
                          capture_output=True, text=True, check=False)
    if proc.returncode == 0 or text not in proc.stderr:
        raise RuntimeError(f"expected a nonzero exit containing {text!r}, got "
                           f"{proc.returncode}\n{proc.stdout}\n{proc.stderr}")
    print(f"rejected ({proc.returncode}): {text}")


def main():
    # Ten Taylor rows cannot witness the unweighted rank 22. The rank checker must fail.
    expect_failure(["-c", "from ranks import test_block; test_block(4,4,2,2,12,10,0)"],
                   "Taylor rank (4, 4, 2, 2, 12, 0, 0): 10 vs 22")
    # The full normal entry point reads the tampered D=5 entry after recomputing the ledger.
    with tempfile.TemporaryDirectory() as directory:
        dst = Path(directory)
        original = HERE/"results"
        (dst/"certificates.json").write_bytes((original/"certificates.json").read_bytes())
        rows = list(csv.DictReader(io.StringIO((original/"degrees.csv").read_text())))
        row = next(r for r in rows if r["D"] == "5")
        if row["mechanical"] != "4":
            raise ValueError("expected clean mechanical degree-five fixture")
        row["mechanical"] = "0"
        with (dst/"degrees.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        expect_failure(["verify.py", "--part", "certificates", "--reference", str(dst)],
                       "degrees.csv differs from the recomputed result")
    print(f"negative controls passed (optimization level {sys.flags.optimize})")


if __name__ == "__main__":
    main()
