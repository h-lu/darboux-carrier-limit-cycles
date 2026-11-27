# Exact reproduction

Run `python3 -B verification/verify.py` from the package root. The entry point recomputes all four result files and compares their complete contents. Use `--part ranks` or `--part certificates` for a single group, and `--write DIR` for a separate copy of the computed results. All scientific calculations use integers, rational arithmetic or prime fields. The programs need no network, CAS, source archive or plotting package.

## Periods and kernels

[ranks.py](ranks.py) computes Taylor coefficients from the angle chart in the proof of the block rank theorem. Write $t=u^2$ in one Morse coordinate. The inverse angle derivative is $J(t)=\sqrt{t/(e^t-1)}$. For normalized frequency $c$, write $\cos(c\theta)=K(t)$ and $\sin(c\theta)=uS(t)$. The two recursions follow from $2K'=-cJS$ and $S+2tS'=cJK$. With signed carrier strength $d$, the weight is $W=e^{t/2}/(1-de^{-t/2})$. Green's formula uses the density $WJK$ and the derivative of $uWS$ with respect to $u$. The program differentiates the weighted vertical primitive before integrating the horizontal-form columns.

For the term $u^{2i}v^{2j}$ in the area density, put $n=i+j$. Integration on the energy ellipse contributes the factor

```text
(pi/sqrt(lambda)) * 2^(1-n)/(n+1)! * lambda^(-i) * (2i)!/i! * (2j)!/j!
```

to the coefficient of $h^{n+1}$. The first three factors up to and including $1/(n+1)!$ are common to a Taylor row. The program removes this nonzero row factor and keeps both orientation signs and the $1/M,1/N$ scales. [relations.py](relations.py) expands the displayed polynomial kernel generators in these same columns. At the rational sample points it checks every coefficient through the stated jet order and the independence of the generator vectors as polynomial forms.

[The rank results](results/ranks.json) use $\lambda=101/37$, $\kappa=1/3$ for the weighted field, and $\kappa=0$ for the unweighted field. `rows` starts at the coefficient of $h$. `prime: 0` denotes rational elimination. Eight small block tests include both signs of the carrier and nonsquare grids. The odd block at $N=M=4,D=11$ has 42 input forms and ranks 22 and 33. Figure 2 shows ten frequency pairs on the smaller unweighted border. Six contribute one form and four contribute two. The 14 border columns and four added forms give an augmented rank of 18. In the full weighted block the kernel has six inherited relations and three multiples of the foliation form. The simple example of a visible weighted direction in the paper has first two normalized rows 0 and 3/16 when the displayed form is doubled.

For each of two primes the program tests the nine blocks at $N=M=8,D=23$, two nonsquare examples and the odd block at $N=M=16,D=47$. It also computes the folded zeros and poles as exact fractions. `shared_upper_x` counts high pairs whose upper horizontal index occurs on the unweighted border in another vertical class. At $D=47$, 28 of the 64 high pairs do so. The augmented rank is 264 and the weighted rank is 492. Nonzero modular minors witness independence at the indicated finite rational specialization. The all-size theorem uses the analytic lifting and joint residue proof in the manuscript.

## Cycle certificates

[certificates.py](certificates.py) reconstructs three grid certificates from the frequency pairs.

1. For the general rectangle it counts odd-step lattice pairs through the total-degree cutoff, applies the block rank, subtracts one per nonempty layer cell and sums with the layer multiplicities. The scan exhausts admissible even grid orders and both carriers. Independent literal enumeration uses the two unweighted border conditions and the two weighted kernel summands to check each family separately through degree 24 and selected larger witnesses.
2. For the unweighted square at $D=2m2^r-1$ it compares the block sum with the first closed polynomial in the square-grid corollary.
3. For the weighted square at $D=3m2^r-1$ it makes the same comparison with the second polynomial. Each symmetric test also checks the exact triangle count, including equality at $u/2^k+v/2^l=j$.

[The certificate results](results/certificates.json) give the 72 symmetric counts for $m=1,2,3,5$ and $r=1,\ldots,9$. They also include the integer-rounded witnesses for both multiplier choices in Theorem 1.1. For every budget from 15 through 20,000 the program checks the binary scale and the degree of each square witness against its budget. At a budget of 31, for example, the rounded weighted square has degree 23. Through degree 350 it compares each rounded certificate with the separate maximum.

[The finite ledger](results/degrees.csv) gives the unweighted and weighted maxima at each total degree from 3 through 350. The columns headed `unweighted` and `weighted` each have their own witness `N,M`. Ties use increasing grid orders. The leading `grid` column takes the larger of the two maxima with precedence given to the unweighted family in a tie. An empty weighted entry means that no even grid fits. The program also checks the replication envelope of each grid family.

[Comparison intervals](results/comparisons.csv) give all strict inequalities, ties and unavailable degrees from 5 through 350. The calculation reads the named columns of the freshly recomputed degree table. The first pair compares the weighted grid to the unweighted grid. The next twelve pairs compare each grid and their maximum separately to the four published protocols. A blank ternary entry at degrees 5–30 is classified as unavailable. The two-source protocol uses the empty initial family, whose zero count at small degrees is available for comparison. Interval endpoints are inclusive. At degree 255 the unweighted witness is the 128 by 128 square with 594,111 cycles, and the weighted witness is the 128 by 64 rectangle with 561,279 cycles.

For the comparison columns, Section 6 of the paper specifies the Han–Li tree and direct seeds, empty-seed two-source recurrence, ternary seeds and binary mechanical recurrence. The classical column is a bound for the number of isolated cycles. The other columns count hyperbolic cycles. A blank ternary entry precedes its first seed. Grid maxima refer to the specified projections of the perturbation space and may correspond to different fields at different budgets.

## Negative controls

`negative_controls.py` starts fresh interpreters with the invoking optimization level. One test supplies only ten Taylor rows to a block of rank 22. Two further tests modify a temporary copy of the ledger. One replaces the mechanical entry 4 at degree five by zero. The other replaces the losing weighted entry at degree 95 by the combined maximum, leaving that maximum intact. A fourth test changes an unavailable ternary comparison to a tie while leaving all counts intact. Each child must exit unsuccessfully with the identified rank or data mismatch. An unexpected successful child exit makes the enclosing test fail. Run the enclosing script once normally and once with `-O`.
