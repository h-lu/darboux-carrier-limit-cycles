# First-return periods on Darboux grids and hyperbolic limit cycles of polynomial fields

Haibo Lu · Shanghai Institute of Technology

[Read the paper](manuscript/main.pdf) · [TeX source](manuscript/main.tex) · [Exact reproduction guide](verification/README.md) · [Finite-degree data](verification/results/degrees.csv) · [Frozen manuscript: v1.0.0](https://github.com/h-lu/darboux-carrier-limit-cycles/tree/v1.0.0)

## What the paper proves

How many hyperbolic limit cycles can coexist in **one** real planar polynomial vector field of total degree at most $D$? Here $D$ bounds the larger of the two component degrees. It is not a monomial budget, nor is the field required to have degree exactly $D$. The paper gives constructive lower bounds, not upper bounds or a determination of the Hilbert number.

For every sufficiently large integer $D$, Theorem 1.1 constructs fields with at least

$$
\frac{c_p D^2}{(\ln 2)^2}\bigl((\ln D)^2-2\ln D\,\ln\ln D\bigr)-K D^2\ln D
$$

distinct hyperbolic limit cycles, with $c_0=3/16$ for the unweighted grid and $c_1=5/24$ for the weighted grid. Each construction uses one perturbation to realize the cycles simultaneously across selected annuli. The constants $c_p$ are leading coefficients in units of $(\ln 2)^{-2}$; they are not exact finite-degree counts.

## Why the multiplier matters

1. **Keep the closed curves.** The zero lines of the Chebyshev polynomials $T_N(x)$ and $T_M(y)$ bound rectangular cells with nested periodic ovals. Multiplying the unperturbed field by $C=1-\kappa T_M(y)$, with $0<\kappa<1/2$, preserves the ovals inside these cells and changes the speed along them.
2. **Change the first variation.** When the polynomial perturbation is added independently of $C$, the first variation of the return map uses the reciprocal weight $1/C$. A relatively exact form that integrates to zero without the weight can give a nonzero weighted period. The paper computes the resulting ranks under the triangular total-degree cutoff. For example, at $N=M=4$, $D=11$ and $\lambda=\pi$, the rank of the 42-form odd block is 22 without the multiplier and 33 with generic multiplier strength. These are **ranks**, not direct counts of cycles.
3. **Turn rank into cycles.** A rank-$r$ block supplies $r-1$ simple return zeros. A hierarchy of coefficients makes the selected simple zeros persist in one field, where they yield hyperbolic limit cycles. The proof fixes a transcendental anisotropy to separate the two frequency indices in the Laplace representation. Finite rational and modular rank checks test examples; they do not substitute for that proof.

The multiplier costs $M$ degrees: the unweighted carrier has degree $N+M-1$ and the weighted one has degree $N+2M-1$. A higher leading coefficient does not imply a better certificate at each finite budget.

## Equal-degree comparisons

The [degree ledger](verification/results/degrees.csv) separately maximizes the two specified grid certificates at every integer budget from 3 to 350 and records each maximizing pair $(N,M)$. The weighted grid first fits at $D=5$. There are no ties between the two grid certificates on their common domain:

| Larger certificate | Integer degree budgets through 350 |
| --- | --- |
| Weighted grid | 79–94, 159–205, 319–350 |
| Unweighted grid | 5–78, 95–158, 206–318 |

For example, at $D=79$ the weighted and unweighted counts are **25,947** and **24,879**, respectively. At $D=255$ the order reverses: the unweighted $128\times128$ square gives **594,111**, while the best weighted rectangle, $128\times64$, gives **561,279**. Weighting that same square would require degree 383.

The leading coefficients of the previous mechanical renewal, the unweighted grid, and the weighted grid are respectively $1/8$, $3/16$, and $5/24$ in the same units. The grid construction already improves the leading coefficient without weighting; the multiplier gives a further increase at a degree cost. Section 6 explains the comparison with Qin–Sun's ternary seeds, distinguishing the original subsequence estimate from its extension to all degrees.

The [comparison ledger](verification/results/comparisons.csv) also compares the two grids and their maximum with four specified baselines. The classical column counts isolated cycles, not necessarily hyperbolic ones; it still exceeds the combined grid certificate at budgets 5–11, 13 and 14. The two-source column uses the **empty-seed** version, and the ternary column is unavailable before degree 31. These are comparisons of particular constructions, not a ranking of all possible fields or all variants of those methods.

## Reproduce the calculations

Run from the repository root with Python 3.9 or later. The exact checks use the standard library and need no network or computer algebra system.

```sh
sha256sum -c SHA256SUMS
python3 -B verification/verify.py
python3 -B -O verification/verify.py
python3 -B verification/negative_controls.py
python3 -B -O verification/negative_controls.py
```

- `verify.py` recomputes and compares four archived files: rank checks, cycle certificates, the finite-degree ledger, and the comparison intervals. For a smaller check use `--part ranks` or `--part certificates`. To write new output separately, use `--write <directory>`.
- `negative_controls.py` tests that the verifier rejects deliberately insufficient rank evidence and altered data, including an error hidden by an unchanged combined maximum.
- For a single degree budget of at least 15, the commands below compute the **square-grid witnesses used in the all-degree proof**, not the optimized finite table.

```sh
python3 -B verification/certificates.py 1000000 --unweighted
python3 -B verification/certificates.py 1000000
```

## Rebuild the paper

With `latexmk` and a TeX installation that includes the packages listed in [main.tex](manuscript/main.tex), use the supplied figure PDFs:

```sh
mkdir -p build
cd manuscript
SOURCE_DATE_EPOCH=1795737600 latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../build main.tex
```

The result is `build/main.pdf` relative to the repository root. Redrawing the two figures is optional and requires NumPy and Matplotlib: run `python3 manuscript/figures/draw.py` from the root in a working copy. See the [technical reproduction guide](verification/README.md) for the exact arithmetic, row conventions, ledger columns, and negative tests.

## Repository map

| Location | Contents |
| --- | --- |
| [`manuscript/`](manuscript/) | Final PDF, TeX, bibliography and figures |
| [`verification/`](verification/) | Exact rank and cycle checks with a technical guide |
| [`verification/results/`](verification/results/) | Expected results and finite comparison data |
| [`SHA256SUMS`](SHA256SUMS) | Integrity hashes for the manuscript and reproduction files |

For citation, use the author and title above and identify the frozen manuscript by tag `v1.0.0`. The tag fixes the manuscript and scientific artifacts even if this landing page is clarified later.
