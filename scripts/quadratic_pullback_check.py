"""Check the homogeneous quadratic pullback identity, including nonunit resultants."""

from itertools import product
from flint import fmpz_poly
from noninterlacing.poly import X
from noncyclotomic_certificates import homogeneous_value


def check():
    seed = X**3 - 4*X - 2
    checked = 0
    for a, b, c, e, f, g in product((-1, 0, 1), repeat=6):
        first = fmpz_poly([c, b, a])
        second = fmpz_poly([g, f, e])
        rho = (a*g-c*e)**2 - (a*f-b*e)*(b*g-c*f)
        if not rho:
            continue
        delta = (b-f*X)**2 - 4*(a-e*X)*(c-g*X)
        pullback = homogeneous_value(seed, first, second)
        expected = seed.discriminant()**2 * rho**6 * seed.resultant(delta)
        assert abs(pullback.discriminant()) == abs(expected)
        checked += 1
    print(f"{checked} degree-two maps matched the binary discriminant formula")


if __name__ == "__main__":
    check()
