"""Certify nonreal disk zeros of greedy digit series and Lehmer's Parry cofactor."""

import cmath
from math import pi

from flint import acb, arb, fmpz_poly

from noninterlacing.pari import pari, to_pari
from noninterlacing.poly import X, chebyshev_C
from noninterlacing.realroots import precision


def disk_parameters(n, ones):
    # Newton only chooses a disk; its rounded center is treated as exact.
    v = 2j*pi - cmath.log(1-cmath.exp(-2j*pi/n))
    for _ in range(12):
        terms = [cmath.exp(-j*v/n) for j in ones]
        v -= (1-sum(terms))/sum(j*t/n for j, t in zip(ones, terms))
    w = acb(arb(round(10**8*v.real))/10**8, arb(round(10**8*v.imag))/10**8)
    terms = [(-j*w/n).exp() for j in ones]
    derivative = abs(sum(j*t/n for j, t in zip(ones, terms)))
    return w, terms, derivative


def tail_allowance(n, ones, params, radius):
    w, terms, derivative = params
    s = arb(radius)/100
    assert w.real > s and w.imag > s and w.imag+s < n*arb.pi()
    error = abs(1-sum(terms)) + sum(
        abs(t)*((j*s/n).exp()-1-j*s/n) for j, t in zip(ones, terms))
    lo = (w.real-s)/n
    # If exp(-lo*next_position) is below this, the geometric tail fits.
    return lo, (derivative*s-error)*(1-(-(n-1)*lo).exp())


def check_first_gap(n):
    ones = (1, n)
    params = disk_parameters(n, ones)
    lo, allowance = tail_allowance(n, ones, params, 35)
    assert (-(2*n-1)*lo).exp() < allowance


def check_gap_tree(n):
    pending = [(1, n)]
    checked = 0
    while pending:
        ones = pending.pop()
        params = disk_parameters(n, ones)
        first = ones[-1]+n-1
        thresholds = []
        for radius in (10, 20, 30, 40, 50, 60, 75):
            lo, allowance = tail_allowance(n, ones, params, radius)
            if allowance > 0:
                j = max(first, int(float(-allowance.log()/lo))+1)
                while not (-j*lo).exp() < allowance:
                    j += 1
                thresholds.append(j)
        if not thresholds:
            raise ArithmeticError(f"uncertified digit prefix {ones}")
        pending.extend(ones+(j,) for j in range(first, min(thresholds)))
        checked += 1
    return checked


def lehmer_greedy():
    P = fmpz_poly([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1])
    beta, = [z.real for z, _ in P.complex_roots() if z.imag.is_zero() and z.real > 1]
    state = X**0
    seen = {}
    digits = []
    while tuple(state.coeffs()) not in seen:
        assert len(digits) < 76
        seen[tuple(state.coeffs())] = len(digits)
        value = (X*state) % P
        if value == 1 or value(beta) > 1:
            digit = 1
        else:
            assert value(beta) < 1
            digit = 0
        digits.append(digit)
        state = value-digit
    pre = seen[tuple(state.coeffs())]
    period = len(digits)-pre
    assert (pre, period) == (1, 74)
    assert [i+1 for i, bit in enumerate(digits) if bit] == [1, 12, 31, 44, 63]
    rel = fmpz_poly([1]+[-i for i in digits])-X**period*fmpz_poly([1,-1])
    cofactor, remainder = divmod(rel, P)
    assert not remainder
    cyclo = X**0
    for order in (2, 4, 12, 22):
        cyclo *= fmpz_poly.cyclotomic(order)
    extra, remainder = divmod(cofactor, cyclo)
    assert not remainder and extra.degree() == 48 and pari().polisirreducible(to_pari(extra))
    assert extra.coeffs() == extra.coeffs()[::-1]
    trace = sum((int(extra[24+j])*chebyshev_C(j) for j in range(1, 25)),
                fmpz_poly([extra[24]]))
    t = to_pari(trace)
    assert pari().polsturm(t) == pari().polsturm(t, [-2, 2]) == 20
    assert trace(2)*trace(-2) != 0
    print("Lehmer preperiod", pre, "period", period, "noncyclotomic cofactor degree", extra.degree())


if __name__ == "__main__":
    with precision(160):
        lehmer_greedy()
        count = sum(check_gap_tree(n) for n in range(18, 46))
        for n in range(46, 4096):
            check_first_gap(n)
        print("certified positions 18 through 4095;", count, "prefixes for 18 through 45")
