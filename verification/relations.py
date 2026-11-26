"""Check the polynomial kernel generators against independently computed period jets."""
from fractions import Fraction

from certificates import pairs, require


def exact_x(a, b, n, lam):
    """Twice K_x(T_a T_b), written as (mode, x-index, y-index, coefficient)."""
    sign = 1 if n > a else -1
    return [("E", a+n, b, b), ("E", abs(a-n), b, b),
            ("F", a+n, b, lam*n+a), ("F", abs(a-n), b, sign*(lam*n-a))]


def exact_y(a, b, m, lam):
    """Twice K_y(T_a T_b)."""
    sign = 1 if m > b else -1
    return [("F", a, b+m, lam*a), ("F", a, abs(b-m), lam*a),
            ("E", a, b+m, m+lam*b), ("E", a, abs(b-m), sign*(m-lam*b))]


def multiply_carrier(terms, m, kap):
    result = list(terms)
    for mode, a, b, c in terms:
        sign = -1 if mode == "E" and b < m else 1
        result.extend([(mode, a, b+m, -kap*c/2),
                       (mode, a, abs(b-m), -kap*c*sign/2)])
    return result


def multiple_omega(a, b, n, m, lam):
    """Four times the scalar T_a T_b times the foliation form Omega."""
    result = []
    for x, sign_x in [(a+n, 1), (abs(a-n), 1 if n>a else -1)]:
        for y, sign_y in [(b+m, 1), (abs(b-m), 1 if m>b else -1)]:
            result.extend([("E", x, y, m*sign_y), ("F", x, y, lam*n*sign_x)])
    return result


def verify(n, m, k, ell, s, lam, kap, columns, p, expected_rank, rank_fn):
    universe = list(columns)
    generators = []
    def take(terms):
        vector = dict.fromkeys(universe, Fraction(0))
        for mode, a, b, coeff in terms:
            require((mode,a,b) in columns and a+b-1 <= s-1, "kernel index outside budget")
            vector[mode,a,b] += coeff
        entries = [vector[key] for key in universe]
        jets = [sum(c*v for c,v in zip(entries, row)) for row in zip(*columns.values())]
        require(all(x == 0 for x in jets), "nonzero jet of a kernel generator")
        generators.append(entries)

    x, y = n >> k, m >> ell
    lower_cutoff = s - (m if p else 0)
    for a,b in pairs(x, y, lower_cutoff-n):
        v = exact_x(a,b,n,lam)
        take(multiply_carrier(v,m,kap) if p else v)
    for a,b in pairs(x, y, lower_cutoff-m):
        v = exact_y(a,b,m,lam)
        take(multiply_carrier(v,m,kap) if p else v)
    if p:
        for a,b in pairs(x, y, s-n-m):
            if b < m:
                take(multiple_omega(a,b,n,m,lam))
    nullity = len(universe)-expected_rank
    require(len(generators) == rank_fn(generators, 0) == nullity,
            "kernel generators fail to span the expected polynomial dimension")
    return nullity
