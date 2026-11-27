# An Improved Asymptotic Lower Bound for Hyperbolic Limit Cycles of Planar Polynomial Vector Fields via Darboux Grids

Haibo Lu · Shanghai Institute of Technology

[Paper PDF](manuscript/main.pdf) · [TeX source](manuscript/main.tex) · [Reproduction guide](verification/README.md) · [Degree-by-degree results](verification/results/degrees.csv)

## Main result

For every sufficiently large integer $D$, a single real planar polynomial vector field of total degree at most $D$ can have at least

$$
\frac{5D^2}{24(\ln 2)^2}
\left((\ln D)^2-2\ln D\,\ln\ln D\right)-K D^2\ln D
$$

distinct hyperbolic limit cycles, for a constant $K>0$. The total degree of a vector field is the larger of the two component degrees. The bound concerns cycles coexisting in one field under this degree budget; there is no restriction on the number of nonzero monomials.

The three leading coefficients below multiply the same scale, $D^2(\ln D)^2/(\ln 2)^2$:

| Construction | Leading coefficient |
| --- | ---: |
| Earlier mechanical renewal | $1/8$ |
| Unweighted Darboux grid | $3/16$ |
| Weighted Darboux grid | $5/24$ |

The unweighted grid raises the earlier coefficient by a factor of $3/2$. Adding the multiplier raises it by a further factor of $10/9$ at the same degree budget. These are improvements in the leading coefficients of constructive lower bounds, rather than comparisons at every finite degree.

## How the construction produces more cycles

The zero lines of two Chebyshev polynomials form rectangular cells. A Darboux first integral gives nested periodic ovals in each cell, and small polynomial perturbations can turn selected ovals into limit cycles.

Multiplying the unperturbed field by a positive polynomial preserves those ovals but changes the speed along them. When the perturbation is added separately, its first-return integral acquires the reciprocal of the multiplier as a weight. Some perturbations whose unweighted integrals vanish then contribute to the return displacement. This increases the number of independent return integrals. The proof uses the two Chebyshev frequency indices to determine their ranks under the total-degree cutoff, then chooses a single perturbation whose simple return zeros persist simultaneously on many annuli.

For example, at grid orders $N=M=4$ and budget $D=11$, the same 42-dimensional odd block has rank 22 without weighting and rank 33 with generic weighting. These are ranks, not counts of limit cycles. The cycle certificate also accounts for which cell layers each block is used on and how many simple zeros are realized there.

## What changes at finite degree

The unweighted carrier has degree $N+M-1$; the weighted carrier has degree $N+2M-1$. At a fixed budget, the extra $M$ degrees can force the weighted construction to use a smaller grid. This explains why the advantage reverses even though the weighted leading coefficient is larger:

| Total degree budget | Unweighted certificate | Weighted certificate |
| ---: | ---: | ---: |
| 78 | 24,399 | 20,799 |
| 79 | 24,879 | 25,947 |
| 95 | 48,071 | 47,783 |
| 255 | 594,111 | 561,279 |
| 319 | 793,839 | 875,631 |

The [finite ledger](verification/results/degrees.csv) separately optimizes the two grid certificates at every integer budget from 3 through 350, recording a witnessing pair $(N,M)$ for each maximum. The weighted grid is available from $D=5$ and gives the larger certificate on 79–94, 159–205, and 319–350. The unweighted grid gives the larger certificate at every other budget in their common range through 350.

[Section 6 of the paper](manuscript/main.pdf) aligns the degree budgets and cycle conventions of the literature comparisons. It treats Qin–Sun's ternary seeds separately on their original degree subsequence and after extension to all sufficiently large budgets. The [comparison ledger](verification/results/comparisons.csv) compares specified constructions; its classical column counts isolated cycles, which need not be hyperbolic.

## Reproduce the results

From the repository root, use Python 3.9 or later. The numerical checks below are exact, use the standard library, and require no network, plotting package, or computer algebra system.

~~~sh
sha256sum -c SHA256SUMS
python3 -B verification/verify.py
python3 -B -O verification/verify.py
python3 -B verification/negative_controls.py
python3 -B -O verification/negative_controls.py
~~~

- **Full verification:** recomputes the rank checks, cycle certificates, and both finite ledgers and compares their contents with the archived results.
- **Focused verification:** add `--part ranks` or `--part certificates` to `verify.py`. Use `--write <directory>` to save recomputed results separately.
- **Negative controls:** check that missing rank evidence, changed counts, and an incorrectly filled unavailable comparison are rejected. Running under `-O` confirms that checks are not silently removed by Python optimization.

To check one of the square-grid witnesses used in the asymptotic proof, supply the desired integer budget (at least 15). These commands evaluate the proof's witnesses, not the optimized finite ledger:

~~~sh
python3 -B verification/certificates.py 1000000 --unweighted
python3 -B verification/certificates.py 1000000
~~~

See the [technical reproduction guide](verification/README.md) for the rank arithmetic, certificate definitions, and ledger columns. The finite exact checks complement the proofs in the paper.

## Rebuild the PDF

With `latexmk` and the TeX packages listed in [the source](manuscript/main.tex), run the following from the repository root. The figure PDFs are included.

~~~sh
mkdir -p build
cd manuscript
SOURCE_DATE_EPOCH=1795737600 latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../build main.tex
~~~

The rebuilt PDF is `build/main.pdf` relative to the repository root. Recreating the figure PDFs is optional and requires NumPy and Matplotlib; run `python3 manuscript/figures/draw.py` from the root of a working copy.

## Versions and files

[`v1.0.1`](https://github.com/h-lu/darboux-carrier-limit-cycles/tree/v1.0.1) is the edition with the revised title and abstract. The theorems, proofs, figures, and reproduction code are unchanged from [`v1.0.0`](https://github.com/h-lu/darboux-carrier-limit-cycles/tree/v1.0.0). [`v1.0.2`](https://github.com/h-lu/darboux-carrier-limit-cycles/tree/v1.0.2) adds the author ORCID, manuscript date, and AI-use declaration. When citing or reproducing the paper, identify the tag you used.

| Path | Contents |
| --- | --- |
| [`manuscript/`](manuscript/) | PDF, LaTeX, bibliography, and figures |
| [`verification/`](verification/) | Rank and cycle checks and the technical guide |
| [`verification/results/`](verification/results/) | Archived exact results and comparison tables |
| [`SHA256SUMS`](SHA256SUMS) | File integrity hashes |
