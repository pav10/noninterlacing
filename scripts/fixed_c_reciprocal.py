"""Exact fixed-root-discriminant census for x^m Psi_n(x-1/x), C=60 and 200."""

from math import exp, gcd, log, prod

from flint import fmpz
from noninterlacing.cyclotomic import factorint, real_cyclotomic_degree, real_cyclotomic_discriminant_factors
from noninterlacing.poly import X, real_cyclotomic
from noncyclotomic_certificates import homogeneous_value


# n, split prime, primitive nth root modulo that prime, two exponents
WITNESSES = {51: (103, 4, 1, 4), 120: (241, 49, 1, 23), 255: (1021, 625, 1, 2)}


def certify(n, discriminant):
    q, w, i, j = WITNESSES[n]
    assert fmpz(q).is_prime() and pow(w, n, q) == 1
    assert all(pow(w, n//p, q) != 1 for p in factorint(n))
    signs = set()
    for k in (i, j):
        assert gcd(k, n) == 1
        u = pow(w, k, q)
        beta = (u + pow(u, -1, q)) % q
        signs.add(pow(beta*beta + 4, (q-1)//2, q))
    assert signs == {1, q-1}
    lifted = homogeneous_value(real_cyclotomic(n), X**2-1, X)
    assert int(lifted.discriminant()) == discriminant


def census(cap):
    hits = []
    for n in range(3, cap**2 + 1):
        if n % 4 == 2 or (4*n > cap**2 and len(factorint(n)) > 1):
            continue
        m = real_cyclotomic_degree(n)
        disc = prod(p**e for p, e in real_cyclotomic_discriminant_factors(n).items())
        if 2**m * disc > cap**m:
            continue
        f = real_cyclotomic(n)
        total = disc**2 * int(f.resultant(X**2+4))
        if total <= cap**(2*m):
            hits.append((2*m, n, total))
    degree = max(d for d, _, _ in hits)
    for d, n, total in hits:
        if d == degree:
            certify(n, total)
            print(f"C={cap}, conductor {n}, degree {d}, root discriminant {exp(log(total)/d):.8f}")


if __name__ == "__main__":
    for cap in (60, 200):
        census(cap)
