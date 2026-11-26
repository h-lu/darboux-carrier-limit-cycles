# Exact reproduction

Run `python3 -B verification/verify.py` from the repository root. The entry point recomputes all three result files and compares their complete contents. Use `--part ranks` or `--part certificates` for a single group, and `--write DIR` for a separate copy of the computed results. All scientific calculations use integers, rational arithmetic or prime fields. The programs need no network, CAS, source archive or plotting package.

## Periods and kernels

[ranks.py](ranks.py) computes Taylor coefficients from the angle chart in the proof of the block rank theorem. Write $t=u^2$ in one Morse coordinate. The inverse angle derivative is $J(t)=\sqrt{t/(e^t-1)}$. For normalized frequency $c$, write $\cos(c\theta)=K(t)$ and $\sin(c\theta)=uS(t)$. The two recursions follow from $2K'=-cJS$ and $S+2tS'=cJK$. With signed carrier strength $d$, the weight is $W=e^{t/2}/(1-de^{-t/2})$. Green's formula uses the density $WJK$ and the derivative of $uWS$ with respect to $u$. The program differentiates the weighted vertical primitive before integrating the horizontal-form columns.

For the term $u^{2i}v^{2j}$ in the area density, put $n=i+j$. Integration on the energy ellipse contributes the factor

```text
(pi/sqrt(lambda)) * 2^(1-n)/(n+1)! * lambda^(-i) * (2i)!/i! * (2j)!/j!
```

to the coefficient of $h^{n+1}$. The first three factors up to and including $1/(n+1)!$ are common to a Taylor row. The program removes this nonzero row factor and keeps both orientation signs and the $1/M,1/N$ scales. [relations.py](relations.py) expands the displayed polynomial kernel generators in these same columns. At the rational sample points it checks every coefficient through the stated jet order and the independence of the generator vectors as polynomial forms.

[The rank results](results/ranks.json) use $\lambda=101/37$, $\kappa=1/3$ for the weighted field, and $\kappa=0$ for the unweighted field. `rows` starts at the coefficient of $h$. `prime: 0` denotes rational elimination. Eight small block tests include both signs of the carrier and nonsquare grids. The odd block at $N=M=4,D=11$ has ranks 22 and 33. Its augmented matrix (unweighted border plus the high terms of the division remainder) has 18 independent columns. The simple example of a visible weighted direction in the paper has first two normalized rows 0 and 3/16 when the displayed form is doubled.

For each of two primes the program tests the nine blocks at $N=M=8,D=23$, two nonsquare examples and the odd block at $N=M=16,D=47$. It also computes the folded zeros and poles as exact fractions. `shared_upper_x` counts high pairs whose upper horizontal index occurs on the unweighted border in another vertical class. At $D=47$, 28 of the 64 high pairs do so. The augmented rank is 264 and the weighted rank is 492. Nonzero modular minors witness independence at the indicated finite rational specialization. The all-size theorem uses the analytic lifting and joint residue proof in the manuscript.

## Cycle certificates

[certificates.py](certificates.py) reconstructs three grid certificates from the frequency pairs.

1. For the general rectangle it counts odd-step lattice pairs through the total-degree cutoff, applies the block rank, subtracts one per nonempty layer cell and sums with the layer multiplicities. The scan exhausts admissible even grid orders and both carriers. Independent literal enumeration checks all grid maxima through degree 24 and selected larger witnesses.
2. For the unweighted square at $D=2m2^r-1$ it compares the block sum with the first closed polynomial in the square-grid corollary.
3. For the weighted square at $D=3m2^r-1$ it makes the same comparison with the second polynomial. Each symmetric test also checks the exact triangle count, including equality on its diagonal.

[The certificate results](results/certificates.json) give the 72 symmetric counts for $m=1,2,3,5$ and $r=1,\ldots,9$. They also include the integer-rounded witnesses from the proof of the all-degree theorem at selected budgets. For all budgets from 15 to 20,000 the program checks both binary rounding inequalities. At overlapping scan budgets it compares the rounded certificate with the full grid search.

[The finite ledger](results/degrees.csv) gives the maximizing grid witness at every $D=3,\ldots,350$. A tie is resolved in increasing order of ($p,N,M$). Each published comparison family uses the replication and unused-degree rule of the finite supplement. Both the prefix-residue and triangular-mode formulas check the mechanical births. A blank ternary entry means that its smallest seed exceeds the budget. Counts on different rows may belong to different fields.

## Negative controls

`negative_controls.py` starts fresh interpreters with the invoking optimization level. One test supplies only ten Taylor rows to a block of rank 22. A second test changes the mechanical entry at degree five to zero in a temporary copy and runs the ordinary certificate entry point. Each child must exit unsuccessfully with the identified mathematical or data mismatch. Run the enclosing script once normally and once with `-O`. An unexpected successful child exit makes the enclosing test fail.
