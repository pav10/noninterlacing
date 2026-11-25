"""Enumerate Salem times cyclotomic witnesses with <=5 terms and measure <=5/4."""

from collections import Counter
from functools import lru_cache
from math import gcd

from flint import fmpq, fmpq_poly, fmpz_poly
from noninterlacing.pari import pari, to_pari
from noninterlacing.poly import X, chebyshev_C, real_cyclotomic
from noninterlacing.realroots import RootData, precision
from salem_five_terms import salem_trace


@lru_cache(None)
def C(n):
    return chebyshev_C(n)


def trace_of_reciprocal(p):
    d = p.degree()//2
    assert p.degree() == 2*d and p.coeffs() == p.coeffs()[::-1] and p[2*d] == 1
    return sum((int(p[d+j])*C(j) for j in range(1, d+1) if p[d+j]), fmpz_poly([p[d]]))


@lru_cache(None)
def cyclotomic_mask(scale):
    return (X**(12*scale)-1)**4


def check_candidate(p, b, scale, counts):
    # All cyclotomic roots have order dividing 12*scale, with multiplicity <=4.
    cofactor = p.gcd(cyclotomic_mask(scale))
    residual = p // cofactor
    if residual.degree() < 4 or residual.degree()%2 or residual[0] != 1:
        return None
    d = residual.degree()//2
    m = d-1
    if 64*m*5**(2*b) < 2**(scale*(m-1))*4**(2*b):
        return None
    tr = trace_of_reciprocal(residual)
    # At t=2 the trace must be negative. The trace of 5/4 is 41/20.
    if tr(2) >= 0 or fmpq_poly(tr)(fmpq(41, 20)) < 0:
        return None
    # A real-rooted polynomial satisfies Laguerre's inequality everywhere.
    der = tr.derivative()
    second = der.derivative()
    for x in (0, 1, -1, 2, -2):
        if der(x)**2 < tr(x)*second(x):
            return None
    counts["Sturm tests"] += 1
    gp = to_pari(tr)
    if pari().polsturm(gp, [-2, 2]) != d-1 or not pari().polisirreducible(gp):
        return None
    return tuple(map(int, tr.coeffs()))


def census():
    counts = Counter()
    hits = {}
    # Degree >=5 cannot occur in an interval of length 2.
    assert 5**5 * 4**4 < 2**20
    for k in range(1, 129):
        h = real_cyclotomic(k)
        if 2 <= h.degree() <= 4:
            roots = RootData(h).balls(128)
            assert roots[-1]-roots[0] > 2
    for scale in (1, 2):
        assert 64*88*5**(2*scale*113) < 2**(scale*87)*4**(2*scale*113)
        assert 89*5**(2*scale) < 88*2**scale*4**(2*scale)

    def record(params, result):
        if result:
            hits.setdefault(result, []).append(params)

    for a in range(2, 113):
        for b in range(1, a):
            if gcd(a, b) != 1:
                continue
            for u in range(-5, 6):
                if not u:
                    continue
                for v in range(-7, 8):
                    if not v or 2+2*u+v > 0:
                        continue
                    if abs(u) == abs(v) == 1:
                        if (u, v) == (-1, -1) and a <= 12:
                            tr = salem_trace(a, b)
                            if tr is not None and fmpq_poly(tr)(fmpq(41, 20)) >= 0:
                                record((a, b, u, v), tuple(map(int, tr.coeffs())))
                        continue
                    p = X**(2*a)+u*(X**(a+b)+X**(a-b))+v*X**a+1
                    counts["five-term inputs"] += 1
                    record((a, b, u, v), check_candidate(p, b, 1, counts))
    for n in range(4, 225):
        for i in range(1, (n+1)//2):
            if gcd(i, n) != 1:
                continue
            for u in (-3, -2, 2, 3):
                for sign in (-1, 1):
                    p = X**n+u*X**(n-i)+sign*u*X**i+sign
                    counts["four-term inputs"] += 1
                    record((n, i, u, sign, "four terms"), check_candidate(p, n-2*i, 2, counts))

    expected_pairs = {(5, 1), (5, 2), (6, 1), (7, 4), (7, 5), (8, 7), (9, 8), (10, 9)}
    assert {w for group in hits.values() for w in group} == {
        (a, b, -1, -1) for a, b in expected_pairs}
    assert len(hits) == 6
    for coeffs, witnesses in hits.items():
        tr = fmpz_poly(list(coeffs))
        root = RootData(tr).balls(128)[-1]
        print("Salem degree", 2*tr.degree(), "measure", (root+(root*root-4).sqrt())/2,
              "witnesses", witnesses)
    print(dict(counts))
    return hits


if __name__ == "__main__":
    with precision(160):
        census()
