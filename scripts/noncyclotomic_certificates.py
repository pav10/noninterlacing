"""Reproduce the discriminant and Salem certificates in notes/noncyclotomic.md."""

from itertools import product
from math import exp, gcd, isqrt, lcm, log, prod

from flint import fmpz, fmpz_mat, fmpz_poly

from noninterlacing.cyclotomic import factorint, real_cyclotomic_degree
from noninterlacing.kronecker import cyclotomic_products, product_poly
from noninterlacing.pari import pari, to_pari
from noninterlacing.poly import X, chebyshev_C, real_cyclotomic
from noninterlacing.realroots import interlaces, sturm_count


def homogeneous_value(f, a, b):
    """Evaluate b^deg(f) f(a/b) over the integer polynomial ring."""
    value = fmpz_poly([f.leading_coefficient()])
    power = fmpz_poly([1])
    for i in range(f.degree() - 1, -1, -1):
        power *= b
        value = value * a + f[i] * power
    return value


def discriminant_certificate():
    n, q, w = 765765, 9189181, 5221919
    primes = list(factorint(n))
    assert fmpz(q).is_prime() and pow(w, n, q) == 1
    assert all(pow(w, n // p, q) != 1 for p in primes)
    for j, expected in [(1, 1), (2, q - 1)]:
        u = pow(w, j, q)
        beta = (u + pow(u, -1, q)) % q
        assert gcd(j, n) == 1
        assert pow(beta**2 + 4, (q - 1) // 2, q) == expected

    # Compare 200^(1/6) prod p^(1-1/(p-1))/(p-1) to 1 in integers.
    exponent = lcm(6, *(p - 1 for p in primes))
    numerator = 200 ** (exponent // 6) * prod(
        p ** (exponent - exponent // (p - 1)) for p in primes
    )
    denominator = prod((p - 1) ** exponent for p in primes)
    assert numerator < denominator
    print("first degree", 2 * real_cyclotomic_degree(n),
          "rd/degree upper bound", (numerator / denominator) ** (1 / exponent))

    for m in (9, 21, 45, 105):
        f = real_cyclotomic(m)
        lift = homogeneous_value(f, X**2 - 1, X)
        norm = f.resultant(X**2 + 4)
        assert lift.discriminant() == f.discriminant() ** 2 * norm


def field_search(base):
    """Exhaust all power-basis generators in the Salem trace region up to 4."""
    d = base.degree()
    basis = [X**0] + [chebyshev_C(i) for i in range(1, d)]
    matrices = []
    for h in basis:
        columns = [(h * X**i) % base for i in range(d)]
        matrices.append(fmpz_mat([[v[j] for v in columns] for j in range(d)]))
    def trace(matrix):
        return sum(matrix[i, i] for i in range(d))

    gram = fmpz_mat([[trace(a * b) for b in matrices]
                     for a in matrices])
    dual = gram.inv()
    bound = 4 * (d - 1) + 16
    radii = [isqrt(int(bound * dual[j, j])) for j in range(d)]
    found = set()
    for c in product(*(range(-r, r + 1) for r in radii)):
        row = fmpz_mat([c])
        if (row * gram * row.transpose())[0, 0] >= bound:
            continue
        matrix = sum((a * m for a, m in zip(c, matrices)), fmpz_mat(d, d))
        h = matrix.charpoly()
        if h.discriminant() != base.discriminant() or h(2) >= 0 or h(4) < 0:
            continue
        if sturm_count(h, -2, 2) == d - 1:
            found.add(tuple(map(int, h.coeffs())))
    print("Salem trace generators up to trace root 4, coefficients ascending", sorted(found))


def salem_certificates():
    base = real_cyclotomic(21)
    candidates = list(cyclotomic_products(5, up_to_reflection=False))
    assert all(not interlaces(base, product_poly(idx)) for idx in candidates)
    print("degree-five Kronecker polynomials excluded", len(candidates))
    field_search(base)

    # beta = 1/(A+1)-2, and 2+beta is a unit.
    assert abs(base(-2)) == 1
    pisot = homogeneous_value(base, -2 * X - 1, X + 1)
    assert pisot.leading_coefficient() == 1 and sturm_count(pisot, -1, 1) == 5
    assert pisot(43) < 0 < pisot(44)
    previous, current = X**0, X - 1
    for k in range(1, 65):
        if k in (1, 4, 12, 32, 64):
            norm = homogeneous_value(pisot, current, previous)
            factor, = [h for h, _ in norm.factor()[1] if h(2) < 0]
            assert factor(-2) != 0
            assert pari().polsturm(to_pari(factor), [-2, 2]) == factor.degree() - 1
            d = factor.degree()
            ratio = exp(log(int(factor.discriminant())) / d) / d
            print("k", k, "Salem trace degree", d, "Salem degree", 2 * d, "rd/degree", ratio)
        previous, current = current, X * current - previous


if __name__ == "__main__":
    discriminant_certificate()
    salem_certificates()
