"""Exact Taylor coefficients from the signed Morse chart and joint-pole tests."""
from fractions import Fraction
from functools import lru_cache
from math import gcd, lcm, prod

from certificates import block_rank, pairs, require, steps
from relations import exact_x, verify as verify_relations


class Field:
    """Rational numbers, or reductions in a specified prime field."""
    def __init__(self, prime):
        self.prime = prime

    def reduce(self, x):
        return x % self.prime if self.prime else x

    def ratio(self, n, d=1):
        if self.prime:
            return n * pow(d, -1, self.prime) % self.prime
        return Fraction(n, d)


def product(a, b, field):
    return [field.reduce(sum(a[i]*b[n-i] for i in range(n+1)))
            for n in range(len(a))]


@lru_cache(maxsize=None)
def geometry(rows, prime):
    """Coefficients of J=theta'(u), in t=u^2.  2 H J'+H'J=0, H=(exp(t)-1)/t."""
    field = Field(prime)
    h, j, moment = [1], [1], [1]
    for n in range(1, rows):
        h.append(field.reduce(h[-1] * field.ratio(1, n+1)))
        coeff = sum((2*n-i)*h[i]*j[n-i] for i in range(1, n+1))
        j.append(field.reduce(-coeff*field.ratio(1, 2*n)))
        moment.append(field.reduce(2*(2*n-1)*moment[-1]))
    return j, moment


@lru_cache(maxsize=None)
def factors(freq, signed_kappa, rows, prime):
    """Even density and differentiated odd primitive in one coordinate."""
    f = Field(prime)
    gamma = f.ratio(freq.numerator, freq.denominator)
    kappa = f.ratio(signed_kappa.numerator, signed_kappa.denominator)
    jac, moment = geometry(rows, prime)
    # cos(gamma*theta(u))=K(t), sin(gamma*theta(u))=u*S(t).
    cosines, sines = [1], []
    growing, falling, weight = [1], [1], []
    inv_constant = f.ratio(1, 1-kappa)
    for n in range(rows):
        if n:
            growing.append(f.reduce(growing[-1]*f.ratio(1, 2*n)))
            falling.append(f.reduce(falling[-1]*f.ratio(-1, 2*n)))
        value = growing[n] + kappa*sum(falling[i]*weight[n-i] for i in range(1, n+1))
        weight.append(f.reduce(value*inv_constant))
        jc = sum(jac[i]*cosines[n-i] for i in range(n+1))
        sines.append(f.reduce(gamma*f.ratio(1, 2*n+1)*jc))
        js = sum(jac[i]*sines[n-i] for i in range(n+1))
        cosines.append(f.reduce(-gamma*f.ratio(1, 2*n+2)*js))
    primitive = product(weight, sines, f)
    density = product(weight, product(jac, cosines[:rows], f), f)
    return ([f.reduce(moment[n]*density[n]) for n in range(rows)],
            [f.reduce(moment[n]*(2*n+1)*primitive[n]) for n in range(rows)])


def jets(n, m, k, ell, cutoff, kappa, rows, prime, lam=Fraction(101, 37)):
    """Columns with the common nonzero factor in each Taylor row divided out."""
    f = Field(prime)
    x_step, y_step = n >> k, m >> ell
    signed = kappa if ell > 1 else -kappa
    inv_lam = f.ratio(lam.denominator, lam.numerator)
    inv_powers = [1]
    for i in range(1, rows):
        inv_powers.append(f.reduce(inv_powers[-1]*inv_lam))
    x_factors = {}
    columns = {}
    for a, b in pairs(x_step, y_step, cutoff):
        if a not in x_factors:
            x_factors[a] = [[f.reduce(w*p) for w, p in zip(line, inv_powers)]
                            for line in factors(Fraction(a, n), 0, rows, prime)]
        x_density, x_deriv = x_factors[a]
        y_density, y_deriv = factors(Fraction(b, m), signed, rows, prime)
        phase = (-1)**((a//x_step-1)//2 + (b//y_step-1)//2)
        for form, left, right, scale in (("E", x_deriv, y_density, f.ratio(phase, m)),
                                         ("F", x_density, y_deriv, f.ratio(-phase, n))):
            columns[form, a, b] = [f.reduce(scale*v) for v in product(left, right, f)]
    return columns


def rank(vectors, prime):
    """Incremental exact elimination. In characteristic zero, use primitive integer rows."""
    basis = {}
    for vector in vectors:
        if prime:
            row = list(vector)
        else:
            denominator = lcm(*(v.denominator for v in vector))
            row = [int(v*denominator) for v in vector]
        for i in range(len(row)):
            if not row[i]:
                continue
            if i not in basis:
                if prime:
                    inv = pow(row[i], -1, prime)
                    row[i:] = [v*inv % prime for v in row[i:]]
                basis[i] = row
                break
            other = basis[i]
            if prime:
                scale = row[i]
                tail = [(a-scale*b) % prime for a, b in zip(row[i+1:], other[i+1:])]
            else:
                x, y = other[i], row[i]
                common = gcd(x, y)
                x, y = x//common, y//common
                tail = [x*a-y*b for a, b in zip(row[i+1:], other[i+1:])]
                content = gcd(*tail)
                if content > 1:
                    tail = [v//content for v in tail]
            row[i] = 0
            row[i+1:] = tail
    return len(basis)


def base_and_poles(c):
    base = c % 2
    if base > 1:
        base = 2 - base
    require(0 < base < 1, "non-dyadic base frequency")
    signed = base if (c-base) % 2 == 0 else -base
    m = (c-signed)/2
    require(m.denominator == 1, "frequency shift")
    zeros = [2*j-1+signed for j in range(1, int(m)+1)]
    poles = [1-2*j-signed for j in range(1, int(m)+1)]
    return base, zeros, poles


def collision_data(n, m, k, ell, s):
    step_x, step_y = n >> k, m >> ell
    pure = [(a,b) for a,b in pairs(step_x, step_y, s-m) if a<n or b<m]
    high = [(a,b) for a,b in pairs(step_x, step_y, s-n) if b<m and a+b>s-n-m]
    clashes = sum(any(a+n == x and b != y for x,y in pure) for a,b in high)
    for a,b in high:
        alpha = 1 + Fraction(a, n)
        beta = Fraction(b, m)
        base, zeros, poles = base_and_poles(alpha)
        pole = min(poles)
        require(pole == -Fraction(a,n) and len(set(poles)) == len(poles), "top pole location")
        residue = prod(pole-z for z in zeros)/prod(pole-z for z in poles if z != pole)
        require(residue != 0, "cancelled residue")
        competitors = list(pure)
        for c,d in high:
            if d != b or c > a:
                continue
            competitors.append((abs(c-n), d))
            if c < a:
                competitors.append((c+n, d))
        for x,y in competitors:
            xb, _, other_poles = base_and_poles(Fraction(x,n))
            yb, _, _ = base_and_poles(Fraction(y,m))
            if (xb, yb) == (base, beta):
                require(pole not in other_poles, "collision within a separated frequency pair")
        for w in [Fraction(3,2), Fraction(2)]:
            actual = pole/(2*w)*(alpha*(w-1)-beta*(pole-1))*residue
            require(actual == pole/(2*w)*alpha*(w-1+beta)*residue != 0,
                    "joint residue cofactor")
    return pure, high, clashes


def test_block(n, m, k, ell, s, rows, prime, augmented=True):
    lam, kap = Fraction(101,37), Fraction(1,3)
    f = Field(prime)
    unweighted = jets(n, m, k, ell, s, 0, rows, prime)
    weighted = jets(n, m, k, ell, s, kap, rows, prime)
    data = {}
    def count(x):
        return len(pairs(n >> k, m >> ell, x))
    for p, columns in enumerate([unweighted, weighted]):
        value = rank(columns.values(), prime)
        expected = block_rank(count, s, n, m, p)
        require(value == expected, f"Taylor rank {(n,m,k,ell,s,p,prime)}: {value} vs {expected}")
        name = ("unweighted", "weighted")[p]
        data[name] = value
        if not prime:
            data[name + "_kernel"] = verify_relations(n, m, k, ell, s, lam, kap, columns,
                                                      p, value, rank)
    if augmented:
        basis_pairs, high, clashes = collision_data(n, m, k, ell, s)
        border = []
        for a,b in basis_pairs:
            if b < m:
                border.append(unweighted["E",a,b])
            if a < n:
                border.append(unweighted["F",a,b])
        coef = f.ratio(n*lam.numerator, lam.denominator)
        for a,b in high:
            sign = 1 if n>a else -1
            e1, e2 = (unweighted["E",x,b] for x in (a+n, abs(a-n)))
            f1, f2 = (unweighted["F",x,b] for x in (a+n, abs(a-n)))
            # Twice D_lambda(A); the sign is required in the folded sine index.
            border.append([f.reduce(m*(u+v) + coef*(w+sign*z))
                           for u,v,w,z in zip(e1,e2,f1,f2)])
        measured = rank(border, prime)
        require(measured == len(border), "augmented Taylor matrix")
        data.update(augmented=measured, high_pairs=len(high), shared_upper_x=clashes)
    return {"N": n, "M": m, "k": k, "l": ell, "S": s, "rows": rows, "prime": prime, **data}


def run():
    # Both carrier signs, unequal grid orders, and a grid with odd multiplier three.
    rational = [test_block(n,m,k,l,s,rows,0)
                for n,m,s,rows in [(2,2,6,12), (4,4,12,38), (4,2,8,20), (6,2,10,20)]
                for k in range(1,len(steps(n))+1) for l in range(1,len(steps(m))+1)]
    example = next(x for x in rational if (x["N"], x["M"], x["k"], x["l"]) == (4,4,2,2))
    require((example["unweighted"],example["weighted"],example["augmented"],
             example["high_pairs"],example["shared_upper_x"]) == (22,33,18,4,1),
            "worked collision block")
    lam = Fraction(101,37)
    visible = jets(2,2,1,1,6,Fraction(1,3),3,0)
    direction = exact_x(1,1,2,lam)
    initial = [sum(c*visible[name,a,b][i] for name,a,b,c in direction) for i in (0,1)]
    require(initial == [0, Fraction(3,16)], "visible weighted direction at the center")
    modulus = [1000000007, 998244353]
    modular = []
    for p in modulus:
        for k in range(1,4):
            for l in range(1,4):
                modular.append(test_block(8,8,k,l,24,135,p))
        for n,m,k,l,s,rows in [(16,8,4,3,32,225), (8,16,3,4,40,365),
                                (16,16,4,4,48,505)]:
            modular.append(test_block(n,m,k,l,s,rows,p))
    high = modular[-1]
    require((high["weighted"],high["augmented"],high["shared_upper_x"]) == (492,264,28),
            "large collision block")
    symmetric = jets(2,2,1,1,6,0,12,0,lam=Fraction(1))
    require(rank(symmetric.values(), 0) < block_rank(lambda x: len(pairs(1,1,x)),6,2,2,0),
            "symmetric anisotropy should lose rank")
    return {"lambda": "101/37", "carrier_kappa": "1/3", "rational_blocks": rational,
            "modular_blocks": modular, "twice_visible_direction_initial_rows": list(map(str, initial)),
            "lambda_one_unweighted_example": rank(symmetric.values(),0)}
