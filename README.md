# Carrier-weighted Darboux polycycle grids

Haibo Lu, Shanghai Institute of Technology. English manuscript R1.

[Read the paper](manuscript/main.pdf). The complete proof, both figures and the finite comparison are in the paper. [TeX entry point](manuscript/main.tex) · [Bibliography](manuscript/references.bib) · [Verification methods and results](verification/README.md).

From the repository root, check the scientific files and reproduce the exact calculations with Python 3.9 or later. These commands use only the standard library.

```sh
sha256sum -c SHA256SUMS
python3 -B verification/verify.py
python3 -B -O verification/verify.py
python3 -B verification/negative_controls.py
python3 -B -O verification/negative_controls.py
```

For the explicit all-degree choice, supply an integer total degree budget of at least 15.

```sh
python3 -B verification/certificates.py 1000000
```

Build the manuscript with pdfLaTeX and BibTeX. The vector figure PDFs are included. To redraw them, run `python3 manuscript/figures/draw.py` with NumPy and Matplotlib installed. To build the paper using the included figures, run the following from the repository root.

```sh
mkdir -p build
cd manuscript
SOURCE_DATE_EPOCH=1795651200 latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../build main.tex
```

The rebuilt PDF is `build/main.pdf`. The manifest identifies the distributed bytes. Font and TeX installations can affect newly built PDFs.
