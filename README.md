# First-return periods on Darboux grids and hyperbolic limit cycles of polynomial fields

Haibo Lu, Shanghai Institute of Technology. English manuscript R2.

[Paper](manuscript/main.pdf) · [TeX entry point](manuscript/main.tex) · [Bibliography](manuscript/references.bib) · [Exact reproduction](verification/README.md) · [Degree-by-degree bounds and grid witnesses](verification/results/degrees.csv).

The paper proves the rank and kernel formulas on dyadic blocks, realizes their simple period zeros in one field and gives the two all-degree estimates. The finite ledger records separate maxima for the unweighted and weighted grids.

Run from this directory using Python 3.9 or later. All exact calculations use the standard library.

```sh
sha256sum -c SHA256SUMS
python3 -B verification/verify.py
python3 -B -O verification/verify.py
python3 -B verification/negative_controls.py
python3 -B -O verification/negative_controls.py
```

For a single total degree budget of at least 15, compute the square-grid witnesses used in the all-degree proof.

```sh
python3 -B verification/certificates.py 1000000 --unweighted
python3 -B verification/certificates.py 1000000
```

To build the paper with the supplied figures, run the following from this directory.

```sh
mkdir -p build
cd manuscript
SOURCE_DATE_EPOCH=1795737600 latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../build main.tex
```

The result is `build/main.pdf`. Both figures can also be redrawn with `python3 manuscript/figures/draw.py` from the bundle root using NumPy and Matplotlib. The manifest covers the supplied PDF, source, figure files, verification programs, results and README files.
