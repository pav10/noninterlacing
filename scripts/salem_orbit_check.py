"""Check the orbit-block comparison and the full determinant on Lehmer's polynomial."""

from math import comb
from flint import acb_mat, arb, ctx, fmpz_poly

from noninterlacing.pari import pari, to_pari, from_pari
from noninterlacing.poly import X, chebyshev_C
from salem_dobrowolski import gram_determinants


def log_determinant(p, powers):
    gen = pari().Mod(pari()("x"), to_pari(p))
    polys = {m: from_pari(pari().minpoly(gen**m)).numer() for m in powers}
    result = arb(0)
    for m, u in powers.items():
        q = polys[m]
        result += u*u*arb(abs(q.discriminant())).log()/2
        for k, v in powers.items():
            if k < m:
                result += u*v*arb(abs(q.resultant(polys[k]))).log()
    return result


def check(p, powers, stride):
    n = p.degree()
    d = n//2
    assert n == 2*d and p.coeffs() == list(reversed(p.coeffs()))
    trace = fmpz_poly([p[d]])
    for j in range(1, d+1):
        if p[d+j]:
            trace += p[d+j]*chebyshev_C(j)
    assert trace(-2) != 0 and trace(2) < 0
    assert pari().polsturm(to_pari(trace), [-2, 2]) == d-1
    size_per_root = sum(powers.values())
    ncols = n * size_per_root
    square_sum = sum(u*u for u in powers.values())
    weighted_sum = sum(m*u for m, u in powers.items())
    slope = sum(m*u*(ncols-u) for m, u in powers.items())
    roots = [z for z, _ in p.complex_roots()]
    salem, = [z.real for z in roots if z.imag.is_zero() and z.real > 1]
    log_measure = salem.log()
    representatives = [z for z in roots if z.imag > 0] + [salem]
    assert 2*len(representatives) == n and 1 <= stride <= n//2

    log_blocks = arb(0)
    for z in representatives:
        rows = []
        for m, u in powers.items():
            for exponent in (m, -m):
                node = z**exponent
                values = [node**0]
                for _ in range(1, ncols):
                    values.append(values[-1] * node)
                for r in range(u):
                    rows.append([comb(j, r)*x for j, x in enumerate(values)])
        a = acb_mat(rows)
        gram = a*a.conjugate().transpose()
        log_blocks += gram.det().real.log()/2

    log_full = log_determinant(p, powers)
    grams = gram_determinants(ncols, powers.values())
    log_coarse = slope*log_measure + n*sum(grams[u] for u in powers.values())/2
    bound = ((n-2)*square_sum*(arb(ncols)/stride).log()/2 + square_sum*arb(ncols).log()
             - 2*size_per_root**2*(1-1/salem).log()
             + (2*stride+2)*size_per_root*weighted_sum*log_measure)
    savings = log_coarse-log_blocks
    assert log_full < log_blocks and 0 < savings and savings < bound
    print("n", n, "columns", ncols, "orbit saving / slope", float(savings/slope),
          "upper bound", float(bound/slope), "cross-pair saving / slope",
          float((log_blocks-log_full)/slope))


if __name__ == "__main__":
    with ctx.workprec(512):
        lehmer = fmpz_poly([1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1])
        check(lehmer, {1: 2, 2: 1, 3: 1}, 2)
        check(lehmer, {1: 3, 2: 1, 3: 1, 5: 1}, 3)
        seed = X**3-X-1
        reverse = fmpz_poly(list(reversed(seed.coeffs())))
        for k in (63, 127):
            p, = [q for q, _ in (X**k*seed-reverse).factor()[1] if q(1) < 0]
            check(p, {1: 2, 2: 1, 3: 1}, 4)
