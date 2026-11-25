"""Bound the optimum autocorrelation cost from both sides for explicit Salem fields."""

from fractions import Fraction
from math import gcd

from flint import fmpz_poly

from noninterlacing.pari import pari, to_pari
from noninterlacing.poly import X, real_cyclotomic


def reciprocal_lift(trace):
    coefficients = trace.coeffs()[::-1]
    lift = fmpz_poly([coefficients[0]])
    for j, c in enumerate(coefficients[1:], start=1):
        lift = (X**2+1)*lift+c*X**j
    return lift


def cost(witness):
    degree = witness.degree()
    correlations = witness*fmpz_poly(witness.coeffs()[::-1])
    energy = int(correlations[degree])
    tail = [abs(int(correlations[degree+j])) for j in range(1, degree+1)]
    return Fraction(4*sum(tail)*sum(j*v for j, v in enumerate(tail, start=1)), energy**2)


def certify(conductor, t):
    f = real_cyclotomic(conductor)
    h = (X-t)*f-f.derivative()
    # The secular equation is increasing between its consecutive poles.
    assert Fraction(int(h(-2)), int(f(-2))) < 0 and t > 2
    assert pari().polisirreducible(to_pari(h))
    n = 2*h.degree()
    M = conductor//gcd(2, conductor)
    binomial = X**M + (1 if conductor%2 == 0 else -1)
    multiplier, remainder = divmod(binomial, fmpz_poly.cyclotomic(conductor))
    assert not remainder
    upper = cost(reciprocal_lift(h)*multiplier)
    lower = M/(1+Fraction(M*M, 2*(t-2)))**2
    assert lower <= upper
    print("conductor", conductor, "Salem degree", n, "optimal cost / n between",
          (lower*10**6//n)/10**6, "and", -((-upper*10**6//n)/10**6))


if __name__ == "__main__":
    for conductor in (30, 210, 2310):
        certify(conductor, 10**9)
