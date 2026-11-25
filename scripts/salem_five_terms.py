"""Classify height-one Salem times cyclotomic polynomials with at most five terms."""

from math import gcd

from flint import arb, fmpz_poly
from noninterlacing.pari import pari, to_pari
from noninterlacing.poly import X, chebyshev_C
from noninterlacing.realroots import RootData, precision


def salem_trace(a, b):
    trace = chebyshev_C(a) - chebyshev_C(b) - 1
    outlier = None
    for factor, multiplicity in trace.factor()[1]:
        interior = int(pari().polsturm(to_pari(factor), [-2, 2]))
        if interior == factor.degree():
            continue
        if (outlier is not None or multiplicity != 1 or factor.degree() < 2
                or interior != factor.degree()-1 or factor(2) >= 0):
            return None
        outlier = factor
    return outlier


def census():
    lehmer_trace = X**5 + X**4 - 5*X**3 - 5*X**2 + 4*X + 3
    least = RootData(lehmer_trace).balls(128)[-1]
    lehmer = fmpz_poly([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1])
    assert X**12-X**7-X**6-X**5+1 == (X**2-X+1)*lehmer
    seed = X**3-X-1
    reverse = 1-X**2-X**3
    assert X**8*seed-reverse == (X-1)*lehmer
    assert 4*arb.pi() < 13
    pairs = []
    ties = []
    factors = set()
    for a in range(2, 13):
        for b in range(1, a):
            if gcd(a, b) != 1:
                continue
            trace = salem_trace(a, b)
            if trace is None:
                continue
            root = RootData(trace).balls(128)[-1]
            if trace == lehmer_trace:
                ties.append((a, b))
            else:
                assert root > least
            pairs.append((a, b))
            factors.add(tuple(map(int, trace.coeffs())))
    assert ties == [(6, 1), (7, 4)]
    print("admissible pairs", pairs)
    print("distinct Salem numbers", len(factors))
    print("Lehmer pairs", ties, "minimum", (least+(least*least-4).sqrt())/2)


if __name__ == "__main__":
    with precision(160):
        census()
