"""Screen ordinary Hilbert class fields of real quadratic fields at a fixed budget."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from heapq import heapify, heappop, heappush
from math import gcd, lcm, prod

from noninterlacing.pari import pari
from unramified_power_generators import collision_valuation


@lru_cache(None)
def residue_factors(p, f):
    gp = pari()
    return [(t, sum(int(gp.moebius(u))*p**(t//u) for u in range(1, t+1) if t%u == 0)//t)
            for t in range(1, f+1) if f%t == 0]


def cross_resultants(p, e, f, places):
    # Each bin is a residue minimal polynomial. Placing the next completion
    # in a bin already containing k completions costs k*2*e*f^2/t.
    slots = [(0, 2*e*f*f//t, min(count, places)) for t, count in residue_factors(p, f)]
    heapify(slots)
    total = 0
    while places:
        cost, step, count = heappop(slots)
        used = min(count, places)
        total += used*cost
        places -= used
        heappush(slots, (cost+step, step, count))
    return total


def screen_field(D, b, n):
    gp = pari()
    group = list(map(int, b.bnf_get_cyc()))
    forced_primes = []
    bound = D**(n//2)
    for prime in map(int, gp.primes([2, n])):
        symbol = int(gp.kronecker(D, prime))
        if symbol == -1:
            e, f = 1, 2
        else:
            e = 1 if symbol else 2
            ideal = gp.idealprimedec(b, prime)[0]
            coords = gp.bnfisprincipal(b, ideal, 0)
            f = lcm(*(a//gcd(a, int(c)) for a, c in zip(group, coords)))
        places = n//(e*f)
        if places <= residue_factors(prime, f)[-1][1]:
            continue
        forced_primes.append(prime)
        increment = cross_resultants(prime, e, f, places)
        if e == 1:
            increment = max(increment, collision_valuation(n, prime**f))
        # Each common index divisor contributes at least two to v(disc polynomial / disc field).
        bound *= prime**max(2, increment)
    return bound, forced_primes


def census(budget, min_degree):
    gp = pari()
    sizes = Counter()
    survivors = []
    locally_monogenic = []
    for D in range(5, budget**2+1):
        if not gp.isfundamental(D):
            continue
        b = gp.bnfinit(gp.quadpoly(D, "y"), 1)
        assert gp.bnfcertify(b)
        n = 2*prod(map(int, b.bnf_get_cyc()))
        if n < min_degree:
            continue
        sizes[n] += 1
        lower, forced = screen_field(D, b, n)
        if lower <= budget**n:
            survivors.append((n, D, forced))
        if not forced:
            locally_monogenic.append((n, D))
    print("degrees and field counts", sorted(sizes.items()))
    print("screened", sum(sizes.values()), "surviving", len(survivors))
    print("unexcluded (degree, base discriminant, common index primes)", sorted(survivors))
    print("without common index divisors", len(locally_monogenic))


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--budget", type=int, default=200)
    parser.add_argument("--min-degree", type=int, default=32)
    args = parser.parse_args()
    census(args.budget, args.min_degree)
