"""Integer certificates for the Darboux grid and the specified comparison ledgers."""
import argparse
import csv
import json
from functools import lru_cache
from io import StringIO

SCAN_LIMIT = 350


def require(condition, message):
    if not condition:
        raise ValueError(message)


def steps(n):
    """One elementary frequency step per layer, in layer order."""
    result = []
    while n % 2 == 0 and n > 0:
        n //= 2
        result.append(n)
    return result


@lru_cache(maxsize=None)
def lattice_table(x_step, y_step, budget):
    """P(X) for 0 <= X <= budget, removing the first vertical column."""
    counts = [0] * (budget + 1)
    for x in range(x_step + y_step, budget + 1):
        first_column = (x - x_step + y_step) // (2 * y_step)
        rest = counts[x - 2 * x_step] if x >= 2 * x_step else 0
        counts[x] = first_column + rest
    return counts


def pairs(x_step, y_step, cutoff):
    return [(a, b) for a in range(x_step, cutoff + 1, 2 * x_step)
            for b in range(y_step, cutoff - a + 1, 2 * y_step)]


def block_rank(count, budget, n, m, carrier):
    """count(X) is the number of admissible frequency pairs below X."""
    if carrier:
        return (2 * count(budget) - 2 * count(budget - n - m)
                - count(budget - 2 * m) + count(budget - n - 2 * m))
    return 2 * count(budget) - count(budget - n) - count(budget - m)


def closed_counts(q, depth, multiplier):
    pure = multiplier**2 * (3*q*q*depth**2 - 2*q*q + 6*q - 4) // 4
    weighted = multiplier**2 * (15*q*q*depth**2 - 2*q*q + 10*q - 8) // 8
    return pure, weighted


def integer_witness(degree, carrier=1):
    """The all-degree square construction with multiplier 1 or 1-kappa*T_M."""
    require(degree >= 15, "the binary-length choice is specified for D >= 15")
    require(carrier in (0, 1), "choose multiplier 0 or 1")
    budget = degree + 1
    bits = budget.bit_length()
    r = (budget // bits).bit_length() - 1
    q = 1 << r
    cost = 2 + carrier
    m = budget // (cost*q)
    n = q*m
    require(q*bits <= budget < 2*q*bits and m >= 1 and r >= 1, "binary scale selection")
    require(cost*n <= budget < cost*(n+q), "rounded grid degree")
    return {"D": degree, "carrier": carrier, "grid_degree": cost*n-1,
            "N": n, "M": n, "q": q, "r": r, "m": m, "cycles": closed_counts(q, r, m)[carrier]}


def check_closed():
    fixtures = ([3, 45, 411, 2967, 18735, 108639, 594111, 3113343, 15794943],
                [8, 120, 1073, 7635, 47783, 275535, 1501343, 7848255, 39748223])
    records = []
    for depth in range(1, 10):
        q = 2**depth
        for multiplier in [1, 2, 3, 5]:
            n = multiplier * q
            for p, closed in enumerate(closed_counts(q, depth, multiplier)):
                budget = (2 + p) * n
                total = 0
                for k in range(1, depth + 1):
                    for ell in range(1, depth + 1):
                        x_step, y_step = n >> k, n >> ell

                        def count(x):
                            return sum((x - a + y_step) // (2 * y_step)
                                       for a in range(x_step, x - y_step + 1, 2*x_step))

                        for h in range(4):
                            lhs = 8 * count(h * n)
                            rhs = h*h * (1 << (k + ell)) + (2*h * (1 << k) if k == ell else 0)
                            require(lhs == rhs, f"triangle boundary at {(k, ell, h)}")
                        rank = block_rank(count, budget, n, n, p)
                        total += x_step * y_step * max(0, rank - 1)
                require(total == closed == multiplier**2 * fixtures[p][depth-1],
                        f"symmetric certificate at {(depth, multiplier, p)}")
                records.append({"carrier": p, "m": multiplier, "r": depth,
                                "D": budget-1, "cycles": total})
    return records


def best_grids(limit):
    """A separate maximum for each multiplier. None means no admissible grid."""
    best = [[None] * (limit + 1) for _ in (0, 1)]
    last = limit + 1
    for p in (0, 1):
        for n in range(2, last - 1, 2):
            for m in range(2, (last - n) // (1 + p) + 1, 2):
                first = n + (1 + p)*m
                totals = [0] * (last + 1)
                for x_step in steps(n):
                    for y_step in steps(m):
                        table = lattice_table(min(x_step, y_step), max(x_step, y_step), last)
                        weight = x_step * y_step
                        for budget in range(first, last + 1):
                            if p:
                                rank = (2*table[budget] - 2*table[budget-n-m]
                                        - table[budget-2*m] + table[budget-n-2*m])
                            else:
                                rank = 2*table[budget] - table[budget-n] - table[budget-m]
                            totals[budget] += weight * max(0, rank - 1)
                for budget in range(first, last + 1):
                    incumbent = best[p][budget-1]
                    if incumbent is None or totals[budget] > incumbent[0]:
                        best[p][budget-1] = (totals[budget], n, m)
    return best


def direct_grid(degree, n, m, p):
    total = 0
    for u in steps(n):
        for v in steps(m):
            # Direct enumeration, independent of the cumulative-table algorithm.
            def border(modes):
                return sum(int(a < n) + int(b < m) for a, b in modes)

            budget = degree+1
            modes = pairs(u, v, budget)
            if p:
                smaller = pairs(u, v, budget-m)
                unweighted_nullity = 2*len(smaller) - border(smaller)
                tangent_nullity = sum(b < m for a, b in pairs(u, v, budget-n-m))
                rank = 2*len(modes) - unweighted_nullity - tangent_nullity
            else:
                rank = border(modes)
            total += u*v * max(0, rank - 1)
    return total


def mechanical(degree):
    t = (degree + 1) // 2
    total = 0
    for stage in range(1, t.bit_length()):
        u = t >> (t.bit_length() - 1 - stage)
        eta = u % 2
        parent = u // 2
        modes = residues = 0
        for k in range(stage):
            a = 2**k
            h, e = divmod(parent, a)
            modes += 4*a*(h*(u - a*h) - 1)
            residues += e*(e + eta)
        other = stage*(u*u - eta) - 4*residues - 4*(2**stage - 1)
        require(modes == other, f"mechanical births at prefix {u}")
        total = 4*total + modes
    return total


def comparison_ledgers(limit):
    """Native states and the shared copy/padding operation in Lu's finite supplement."""
    seeds = {2: 4, 4: 28, 5: 37, 6: 53, 7: 74, 8: 96, 9: 120, 10: 142,
             13: 212, 17: 384, 21: 568, 31: 1184, 35: 1536, 39: 1920, 43: 2272}
    tree = {3: (1, 13), 4: (1, 20), 5: (2, 28), 6: (2, 35)}
    classical, two_source, ternary, mech = ([0]*(limit+1) for _ in range(4))
    for d in range(2, limit+1):
        if d >= 7:
            parent = (d - 1) // 2
            r, old = tree[parent]
            births = (4*r + 1)*parent - (2*r if d % 2 else 0)
            tree[d] = (2*r + 1, 4*old + births)
        classical[d] = max(tree.get(d, (0, 0))[1], seeds.get(d, 0))
        if d >= 11:
            parent = (d - 5) // 2
            two_source[d] = 4*two_source[parent] + 4*(parent+2)*(parent+1)
        mech[d] = mechanical(d)
    r = 2
    while 4*3**r - 5 <= limit:
        q = 3**r
        ternary[4*q-5] = 4*r*r*(q//3)**2 - (q-1)**2
        r += 1

    def envelope(ledger, d):
        return max(ledger[j]*((d+1)//(j+1))**2 for j in range(2, d+1))

    return [tuple(envelope(ledger, d) for ledger in (classical, two_source, ternary, mech))
            if d >= 2 else (0, 0, 0, 0) for d in range(limit+1)]


def run():
    symmetric = check_closed()
    limit = SCAN_LIMIT
    for u in range(1, 9):
        for v in range(1, 9):
            table = lattice_table(u, v, 40)
            require(all(table[x] == len(pairs(u, v, x)) for x in range(41)),
                    f"pair counter with steps {(u, v)}")
    by_family, other = best_grids(limit), comparison_ledgers(limit)
    grid = [None] * (limit+1)
    for d in range(3, limit+1):
        unweighted, weighted = by_family[0][d], by_family[1][d]
        p = int(weighted is not None and weighted[0] > unweighted[0])
        value, n, m = by_family[p][d]
        grid[d] = (value, p, n, m)
    for d in range(3, 25):
        for p in (0, 1):
            candidates = ((direct_grid(d, n, m, p), n, m) for n in range(2, d, 2)
                          for m in range(2, (d+1-n)//(1+p)+1, 2))
            expected = max(candidates, key=lambda entry: entry[0], default=None)
            require(expected == by_family[p][d], f"exhaustive search at {(d, p)}")
    historical_rows = {11: (153,80,0,72), 23: (833,440,0,676),
                       31: (1801,1320,80,1940), 63: (9065,7912,320,12756),
                       103: (23029,24792,2240,34948), 319: (280073,394472,40256,652116)}
    for d, reference in historical_rows.items():
        require(other[d] == reference, f"supplement table at {d}")
    grid_fixtures = {3: 3, 5: 11, 7: 45, 11: 140, 23: 1161, 47: 7907, 63: 18735,
                     95: 48071, 191: 275535, 255: 594111}
    for d, reference in grid_fixtures.items():
        val, p, n, m = grid[d]
        require(val == reference == direct_grid(d, n, m, p), f"grid witness at {d}")
    require(mechanical(5) == other[5][3] == 4, "mechanical degree-five certificate")
    require(all(other[d][3] == mechanical(d) for d in range(3, limit+1)),
            "mechanical envelope differs from native output")
    pure_wins = sum(by_family[0][d][0] > other[d][3] > 0 for d in range(5, limit+1))
    all_pure_wins = sum(by_family[0][d][0] > other[d][3] for d in range(3, limit+1))
    first_weighted = next(d for d in range(5, limit+1) if by_family[1][d][0] > by_family[0][d][0])
    require(pure_wins == 346 and all_pure_wins == 348 and first_weighted == 79, "separate finite degree comparisons")
    for p, first in ((0, 3), (1, 5)):
        for d in range(first, limit+1):
            copied = max(by_family[p][j][0]*((d+1)//(j+1))**2 for j in range(first, d+1))
            require(copied == by_family[p][d][0], f"grid replication bound at {(d, p)}")
    for p, witness in enumerate(((24879, 32, 32), (25947, 16, 32))):
        value, n, m = witness
        require(witness == by_family[p][79] and value == direct_grid(79, n, m, p),
                f"degree 79 witness, multiplier {p}")
    out = StringIO()
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(["D", "grid", "carrier", "N", "M", "classical", "two_source",
                     "Qin_Sun", "mechanical", "unweighted", "unweighted_N", "unweighted_M",
                     "weighted", "weighted_N", "weighted_M"])
    for d in range(3, limit+1):
        val, p, n, m = grid[d]
        c, t, q, mech = other[d]
        writer.writerow([d, val, p, n, m, c, t, q if q else "", mech,
                         *by_family[0][d], *(by_family[1][d] or ("", "", ""))])
    for p in (0, 1):
        for d in range(15, 20001):
            witness = integer_witness(d, p)
            if d <= limit:
                require(witness["cycles"] <= by_family[p][d][0], f"rounded witness at {(d, p)}")
    stats = {"symmetric_certificates": symmetric, "scan_degrees": [3, limit],
             "rounded_witness_budgets_checked": [15, 20000],
             "rounded_witness_examples": [integer_witness(d,p) for d in (31,103,1000,10000,10**6,10**9)
                                          for p in (0,1)],
             "strict_unweighted_over_mechanical_5_to_350": pure_wins,
             "strict_unweighted_over_mechanical_3_to_350": all_pure_wins,
             "first_strict_weighted_over_unweighted": first_weighted, "mechanical_D5": mechanical(5)}
    return stats, out.getvalue()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute a square-grid witness at fixed total degree.")
    parser.add_argument("degree", type=int, help="total degree budget, at least 15")
    parser.add_argument("--unweighted", action="store_true", help="use the constant multiplier")
    args = parser.parse_args()
    print(json.dumps(integer_witness(args.degree, 0 if args.unweighted else 1), indent=2))
