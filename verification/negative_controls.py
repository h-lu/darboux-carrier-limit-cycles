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
    # Alter a published comparator, then a losing grid entry with the maximum intact.
    with tempfile.TemporaryDirectory() as directory:
        dst = Path(directory)
        original = HERE/"results"
        (dst/"certificates.json").write_bytes((original/"certificates.json").read_bytes())
        rows = list(csv.DictReader(io.StringIO((original/"degrees.csv").read_text())))
        for degree, field, clean, replacement in (("5", "mechanical", "4", "0"),
                                                  ("95", "weighted", "47783", "48071")):
            row = next(r for r in rows if r["D"] == degree)
            if row[field] != clean:
                raise ValueError(f"unexpected clean {field} entry at degree {degree}")
            row[field] = replacement
            with (dst/"degrees.csv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(row), lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            print(f"altered {field} at D={degree} from {clean} to {replacement}")
            expect_failure(["verify.py", "--part", "certificates", "--reference", str(dst)],
                           "degrees.csv differs from the recomputed result")
            row[field] = clean
    print(f"negative controls passed (optimization level {sys.flags.optimize})")


if __name__ == "__main__":
    main()
