"""Exhaust normalized real cyclotomic fibres at the endpoint x=2."""

import argparse
from functools import lru_cache
from itertools import product
from math import comb

from flint import arb, fmpq, fmpq_poly
from noninterlacing.kronecker import cyclotomic_products, product_poly
from noninterlacing.pari import pari, to_pari
from noninterlacing.poly import X, is_irreducible, real_cyclotomic
from noninterlacing.realroots import RootData, precision


@lru_cache(None)
def outer_polynomials(shift):
    out = []
    for m in range(1, 43):
        g = real_cyclotomic(m)(X-shift)
        if g.degree() < 2:
            continue
        if pari().polsturm(to_pari(g), [0, shift + 2]) == 1:
            roots = RootData(g, check=False).balls(128)
            out.append((m, g, roots[0], roots[-1]))
    return out


@lru_cache(None)
def small_outer_polynomials():
    """One root in (0,7], all others in [-2,0), degree two or three."""
    from flint import fmpz_poly

    out = []
    for q in (2, 3):
        # Bound each elementary symmetric function, separating the positive root.
        bounds = [comb(q-1, j)*2**j + 7*comb(q-1, j-1)*2**(j-1)
                  for j in range(1, q+1)]
        for descending in product(*(range(-v, v+1) for v in bounds)):
            g = fmpz_poly([*reversed(descending), 1])
            if g(0) >= 0:
                continue
            gp = to_pari(g)
            if (pari().polsturm(gp, [-2, 0]) != q-1
                    or pari().polsturm(gp, [0, 7]) != 1 or not is_irreducible(g)):
                continue
            roots = RootData(g, check=False).balls(128)
            out.append((str(g), g, roots[0], roots[-1]))
    return out


def fibre_bounds(t, endpoint):
    # Read extrema from right to left, starting with a local minimum.
    critical = RootData(t.derivative(), check=False).balls(128)
    values = [t(r) for r in reversed(critical)]
    lows, highs = values[::2], values[1::2] + [t(endpoint)]
    (lows if t.degree() % 2 else highs).append(arb(t(-2)))
    # Selecting any one lower and upper constraint gives a safe enclosure.
    return max(lows, key=float), min(highs, key=float)


def census(b, cap):
    endpoint = arb(cap.p)/arb(cap.q)
    count = 0
    hits = []
    for indices in cyclotomic_products(b-1, up_to_reflection=False):
        if min(indices) <= 2:
            continue
        count += 1
        t = (X-2)*product_poly(indices)
        # A level strictly between the smallest and largest outer root has
        # simple real fibres. If both tails stop by +/-1, its nonzero integer
        # norm excludes the outer polynomial (degree >= 2).
        tp = to_pari(t)
        if (pari().polsturm(tp+1, [-2, 2]) < b
                and (fmpq_poly(t)(cap) < 1 or pari().polsturm(tp-1) < b)):
            continue
        low, high = fibre_bounds(t, endpoint)
        if low >= 0 or high <= 0:
            continue
        shift = next((c for c in (-1, 0, 1) if low >= c-2 and high <= c+2), None)
        if shift is not None:
            candidates = outer_polynomials(shift)
        elif ((low >= -arb(17)/12 and high <= arb(10)/3)
              or (low >= -arb(1)/2 and high <= 7)):
            # The norms of y(y+1)^2 and y, respectively, force degree <= 3.
            candidates = small_outer_polynomials()
        else:
            raise ArithmeticError(f"interval requires another outer search: {indices}, {low}, {high}")
        for m, outer, least, greatest in candidates:
            if least < low or greatest > high:
                continue
            f = outer(t)
            if f(2) >= 0 or f(-2) == 0 or not is_irreducible(f):
                continue
            if pari().polsturm(to_pari(f), [-2, 2]) != f.degree()-1:
                continue
            if fmpq_poly(f)(cap) < 0:
                continue
            r = RootData(f, check=False).balls(128)[-1]
            hits.append((float((r+(r*r-4).sqrt())/2), f.degree(), indices, m, shift))
    return count, sorted(hits)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-degree", type=int, default=18)
    parser.add_argument("--all-measures", action="store_true",
                        help="classify inner degrees 3 and above without a measure cutoff")
    args = parser.parse_args()
    with precision(160):
        for b in range(3 if args.all_measures else 2, args.max_degree+1):
            cap = fmpq(4*b) if args.all_measures else fmpq(41,20)
            count, hits = census(b, cap)
            print("inner degree", b, "fibres", count, "trace cap", cap, "hits", hits, flush=True)
        count, hits = census(3, fmpq(23,10))
        print("cubic inner degree, trace root <= 23/10", hits)
