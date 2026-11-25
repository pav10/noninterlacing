"""Certified finite-degree Mahler measure bounds for Salem numbers.

Run with the cloud Python. The parameter choices affect strength only.
"""

from math import comb
from flint import arb, ctx, fmpz, fmpz_mat

PRIMES = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
          47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109,
          113, 127, 131, 137, 139, 149, 151, 157, 163, 167)

# Salem degree, base multiplicity, prime multiplicities, lower bounds for
# log(lambda) in millionths (coarse row estimate, reciprocal-pair estimate).
CASES = [
    (32, 45, [3, 5, 8, 8, 8, 8, 6, 4, 3], 32920, 34785),
    (128, 37, [2, 4, 5, 6, 6, 6, 6, 5, 5, 4, 3, 2, 2], 21561, 21927),
    (512, 36, [1, 3, 4, 5, 5, 6, 6, 6, 5, 5, 4, 4, 4, 3, 2, 2, 1, 1],
     14949, 15019),
    (1024, 35, [1, 2, 4, 4, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 3, 3, 2, 2, 1, 1],
     12738, 12769),
    (4096, 39, [1, 2, 4, 4, 5, 5, 5, 6, 6, 6, 5, 5, 5, 5, 5, 4, 4, 3, 3, 3, 2, 2, 1],
     9550, 9557),
    (1048576, 45, [1, 2, 3, 3, 4, 4, 4] + [5]*19 + [4]*6 + [3]*5, 4026, 4026),
]


def gram_determinants(size, orders):
    """Squared volumes for the rows binom(j,r), 0 <= j < size."""
    k = max(orders)
    choices = [comb(size, h + 1) for h in range(2*k - 1)]
    entries = [[sum(comb(r + s - t, r) * comb(r, t) * choices[r + s - t]
                    for t in range(min(r, s) + 1)) for s in range(k)] for r in range(k)]
    return {u: arb(fmpz_mat([row[:u] for row in entries[:u]]).det()).log()
            for u in set(orders)}


def certificate_gap(n, powers, t, log_grams, refined):
    """Positive lower minus upper log determinant excludes log(lambda) <= t."""
    size = n * sum(powers.values())
    d = n // 2
    log_disc = 4 * (d * arb(d).log() - arb(d + 1).lgamma())
    lower = sum(u*u for u in powers.values()) * log_disc / 2
    lower += n * powers[1] * sum(u * arb(p).log() for p, u in powers.items() if p > 1)
    upper = n * sum(log_grams[u] for u in powers.values()) / 2
    for m, u in powers.items():
        w = m * t
        z = w * u * (size - u)
        if refined:
            # Both bounds majorize the same nondecreasing finite cosh sum.
            h1 = z.cosh().log()
            h2 = z - u*u * (-(-2*w).expm1()).log() - log_grams[u]
            pair = (size*w).sinh() / (size*w.sinh())
            upper += (h1 if h1 < h2 else h2) + (1 - 1/pair**2).log()/2
        else:
            upper += z
    return lower - upper


def abelian_trace_gap(n, t):
    """Prime product forced into the trace discriminant if log(lambda) <= t."""
    log_product = arb(0)
    for p in range(5, 10001):
        if fmpz(p).is_prime() and (n // 2) * (arb(p)/4).log() > (p + 1) * t:
            log_product += arb(p).log()
    return log_product - 2 * arb(n).log() - 4*t


def main():
    with ctx.workprec(192):
        for n, base, weights, coarse, refined in CASES:
            powers = {1: base, **dict(zip(PRIMES, weights))}
            grams = gram_determinants(n * sum(powers.values()), powers.values())
            for value, use_pairs in ((coarse, False), (refined, True)):
                gap = certificate_gap(n, powers, arb(value)/10**6, grams, use_pairs)
                if not gap > 0:
                    raise ArithmeticError(f"degree {n}, uncertified log determinant gap {gap}")
            print(f"degree {n:4d}   coarse {coarse/1e6:.6f}   paired {refined/1e6:.6f}")
        for n, thousandths in [(32, 1166), (128, 3746), (512, 10361),
                               (1024, 16657), (4096, 40307), (1048576, 935252)]:
            if not abelian_trace_gap(n, arb(thousandths)/1000) > 0:
                raise ArithmeticError(f"uncertified abelian trace bound for degree {n}")
            print(f"degree {n:4d}   abelian trace {thousandths/1000:.3f}")


if __name__ == "__main__":
    main()
