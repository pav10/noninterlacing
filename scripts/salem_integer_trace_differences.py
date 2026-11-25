"""Classify Salem factors with two positive trace powers differing by an integer."""

from collections import defaultdict
from math import comb

from flint import arb, fmpz_poly

from noninterlacing.pari import pari, to_pari
from noninterlacing.poly import X, chebyshev_C
from noninterlacing.realroots import RootData, precision
from salem_dobrowolski import PRIMES, certificate_gap, gram_determinants


def check_log_bound(max_degree, denominator, base_weight, prime_count):
    weights = {1: base_weight, **dict.fromkeys(PRIMES[:prime_count], 1)}
    for d in range(2, max_degree+1):
        n = 2*d
        grams = gram_determinants(n*sum(weights.values()), weights.values())
        assert certificate_gap(n, weights, arb(1)/denominator, grams, False) > 0


def check_span_cutoff(d, width, twice_power_log):
    # This left hand side is increasing for integer degrees from d onward.
    growth = arb(4)/width
    assert growth > arb(d)/(d-1)
    assert (d-2)*growth.log() - arb(64*(d-1)).log() > twice_power_log


def classify():
    for c in (1, 2, 3):
        width = 4 - arb(c*c)/4
        bound = 1 if c == 1 else 2
        assert 2*(arb(bound).cosh()-1) > c
        check_span_cutoff(200, width, 2*bound)
    check_log_bound(199, 100, 8, 15)

    for c in (1, 2, 3):
        bound = 5 if c == 1 else 6
        assert 2*arb(bound).sinh()/100 > c
        check_span_cutoff(70, 4-c, 2*bound)
    check_log_bound(69, 55, 6, 10)

    powers = [chebyshev_C(i) for i in range(301)]
    levels = [2*(arb(i)/55).cosh() for i in range(301)]
    assert levels[300]-levels[299] > 3

    hits = defaultdict(list)
    counts = [0]*4
    for a in range(2, 300):
        for b in range(1, a):
            delta = levels[a]-levels[b]
            for c in (1, 2, 3):
                if delta > c:
                    continue
                counts[c] += 1
                trace = powers[a]-powers[b]-c
                # The unique zero above 2 is simple. Only its factor is relevant.
                exceptional, = [f for f, e in trace.factor()[1] if f(2) < 0]
                degree = exceptional.degree()
                if 2 <= degree < 70 and pari().polsturm(to_pari(exceptional), [-2, 2]) == degree-1:
                    hits[tuple(map(int, exceptional.coeffs()))].append((a, b, c))

    lehmer = (3, 4, -5, -5, 1, 1)
    greedy = 1-X-X**12-X**31-X**44-X**63-X**74+X**75
    other = X**32-X**31-2*X**16-X+1
    lehmer_poly = fmpz_poly([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1])
    assert greedy.gcd(other) == (X**2+1)*lehmer_poly
    k = 1024
    assert 63**k * 4**(3*k) < 64**k * comb(2*k, k)**3
    assert hits[lehmer] == [(6, 1, 1), (7, 4, 1), (16, 15, 2)]
    least = RootData(fmpz_poly(list(lehmer))).balls(128)[-1]
    for coeffs, parameters in hits.items():
        t = fmpz_poly(list(coeffs))
        r = RootData(t).balls(128)[-1]
        assert coeffs == lehmer or r > least
        print(t, parameters)
    assert counts[1:] == [2675, 4490, 5960]
    assert len(hits) == 19 and sum(map(len, hits.values())) == 29
    print("pairs tested by difference", counts[1:], "distinct Salem numbers", len(hits))


if __name__ == "__main__":
    with precision(160):
        classify()
