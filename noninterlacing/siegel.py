"""Siegel's formula for ``zeta_K(-1)`` and the "Siegel ceiling" for ``m(O_K)``.

For a totally real field ``K`` of degree ``d`` the restriction to the diagonal
of the Hecke Eisenstein series of parallel weight 2 is a modular form of
weight ``2d`` for ``SL_2(Z)``:

    zeta_K(-1)/2^d + sum_{k>=1} s_K(k) q^k,
    s_K(k) = sum_{nu in d^{-1}, nu >> 0, Tr nu = k} sigma_1(nu d),

where ``d`` is the different and ``sigma_1(a) = sum_{b | a} N(b)``.  Since a
modular form of weight ``w`` is determined by its first ``r = dim M_w``
coefficients ``a_1..a_r``, Siegel obtained rational ``b_k(w)`` with

    zeta_K(-1) = 2^d * sum_{k=1}^{r} b_k(2d) s_K(k).

As ``zeta_K(-1) != 0`` some ``s_K(k)``, ``k <= r``, is non-zero, i.e.
``m(O_K) <= r = dim M_{2d}`` (the Siegel ceiling); in particular
``m(O_K) = 1`` when ``r = 1`` (``d <= 5`` and ``d = 7``), which is the content
of Yatsyna's Theorem 3.  Since ``Z[alpha] subset O_K``, every monic totally
real ``f`` of degree ``d`` satisfies ``m(f) <= dim M_{2d}``.

The functions here compute ``b_k(w)`` exactly, tabulate signs, and evaluate
the right-hand side of Siegel's formula from our enumeration of totally
positive codifferent elements (ideal factorisation through PARI), which can
be compared with PARI's ``zeta_K(-1)``: an independent end-to-end check
that the enumeration is complete.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import comb


@lru_cache(maxsize=None)
def bernoulli(n: int) -> Fraction:
    """Bernoulli numbers with ``B_1 = -1/2``."""
    B = [Fraction(1)]
    for m in range(1, n + 1):
        B.append(-sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1))
    return B[n]


def sigma(k: int, n: int) -> int:
    s, i = 0, 1
    while i * i <= n:
        if n % i == 0:
            s += i ** k
            j = n // i
            if j != i:
                s += j ** k
        i += 1
    return s


def eisenstein(k: int, prec: int) -> list:
    """q-expansion of ``E_k = 1 - (2k/B_k) sum sigma_{k-1}(n) q^n`` (``k >= 4`` even)."""
    c = Fraction(-2 * k) / bernoulli(k)
    return [Fraction(1)] + [c * sigma(k - 1, n) for n in range(1, prec)]


def _mul(a: list, b: list, prec: int) -> list:
    out = [Fraction(0)] * prec
    for i, x in enumerate(a[:prec]):
        if x:
            for j in range(min(len(b), prec - i)):
                out[i + j] += x * b[j]
    return out


def dim_M(w: int) -> int:
    """Dimension of ``M_w(SL_2(Z))``."""
    if w < 0 or w % 2:
        return 0
    if w == 2:
        return 0
    return w // 12 + (0 if w % 12 == 2 else 1)


def modular_basis(w: int, prec: int) -> list:
    """The basis ``E_4^a E_6^b`` (``4a + 6b = w``) of ``M_w``."""
    E4, E6 = eisenstein(4, prec), eisenstein(6, prec)
    out = []
    for b in range(w // 6 + 1):
        rest = w - 6 * b
        if rest % 4:
            continue
        a = rest // 4
        f = [Fraction(1)] + [Fraction(0)] * (prec - 1)
        for _ in range(a):
            f = _mul(f, E4, prec)
        for _ in range(b):
            f = _mul(f, E6, prec)
        out.append(f)
    return out


def _solve(A: list, y: list) -> list:
    """Exact Gaussian elimination for a square rational system."""
    n = len(A)
    M = [list(A[i]) + [y[i]] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        M[c] = [v / piv for v in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                fac = M[r][c]
                M[r] = [u - fac * v for u, v in zip(M[r], M[c])]
    return [M[i][n] for i in range(n)]


@lru_cache(maxsize=None)
def siegel_coefficients(w: int) -> tuple:
    """Rationals ``b_1..b_r`` (``r = dim M_w``) with ``a_0 = sum b_k a_k`` on ``M_w``."""
    r = dim_M(w)
    if r == 0:
        return ()
    basis = modular_basis(w, r + 1)
    A = [[f[k] for k in range(1, r + 1)] for f in basis]
    y = [f[0] for f in basis]
    return tuple(_solve(A, y))


def siegel_ceiling(d: int) -> int:
    """``dim M_{2d}``: upper bound for ``m(O_K)`` (hence for every ``m(f)``) in degree ``d``."""
    return dim_M(2 * d)


def coefficient_table(degrees) -> list:
    """Rows ``(d, dim M_{2d}, b_1(2d), ..., sign pattern)``."""
    rows = []
    for d in degrees:
        b = siegel_coefficients(2 * d)
        rows.append({
            "d": d,
            "dim M_2d": len(b),
            "b": [str(x) for x in b],
            "signs": "".join("+" if x > 0 else "-" if x < 0 else "0" for x in b),
            # s_1 = 0 (non-interlacing O_K) needs sum_{k>=2} b_k s_k to have sign (-1)^d
            "sign zeta_K(-1)": "+" if d % 2 == 0 else "-",
        })
    return rows


def noninterlacing_sign_test(d: int) -> str:
    """What Siegel's formula says about ``s_1 = 0`` (``m(O_K) >= 2``) in degree ``d``.

    If ``dim M_{2d} = 1`` it is impossible.  If ``dim M_{2d} = 2`` it forces
    ``zeta_K(-1) = 2^d b_2 s_2`` with ``s_2 > 0``; this is consistent with the
    sign ``(-1)^d`` of ``zeta_K(-1)`` iff ``(-1)^d b_2 > 0``.
    """
    b = siegel_coefficients(2 * d)
    if len(b) <= 1:
        return "impossible (dim M_2d <= 1)"
    if len(b) == 2:
        ok = ((-1) ** d) * b[1] > 0
        return "allowed by signs" if ok else "excluded by signs"
    return "not decided by signs (dim M_2d >= 3)"


# ---------------------------------------------------------------------------
# Evaluating Siegel's formula from enumerated codifferent elements (PARI)
# ---------------------------------------------------------------------------

def sigma1_ideal(nf, ideal) -> int:
    """``sigma_1`` of an integral ideal (PARI ideal in ``nf``)."""
    from .pari import pari

    P = pari()
    fa = P.idealfactor(nf, ideal)
    s = 1
    for i in range(int(P.matsize(fa)[0])):
        pr, e = fa[i, 0], int(fa[i, 1])
        q = int(P.idealnorm(nf, pr))
        s *= (q ** (e + 1) - 1) // (q - 1)
    return s


def s_K(f, k: int, order=None) -> tuple:
    """``(s_K(k), number of nu)`` for ``K = Q[x]/(f)`` via enumeration + PARI."""
    from .orders import Order
    from .pari import nfinit, pari, to_pari

    P = pari()
    O = order or Order.maximal(f)
    nf = nfinit(O.f)
    diff = nf.nf_get_diff()
    total, count = 0, 0
    for nu in O.totally_positive_codifferent(k):
        I = P.idealmul(nf, to_pari(nu), diff)
        total += sigma1_ideal(nf, I)
        count += 1
    return total, count


def siegel_zeta(f, max_k: int | None = None) -> dict:
    """``zeta_K(-1)`` from Siegel's formula with enumerated ``s_K(k)``; compare to PARI."""
    from .orders import Order
    from .pari import zeta_at_minus_one
    from .poly import poly

    f = poly(f)
    d = f.degree()
    b = siegel_coefficients(2 * d)
    O = Order.maximal(f)
    s, counts = [], []
    for k in range(1, len(b) + 1):
        sk, ck = s_K(f, k, O)
        s.append(sk)
        counts.append(ck)
    val = Fraction(2 ** d) * sum((bk * sk for bk, sk in zip(b, s)), Fraction(0))
    z, x = zeta_at_minus_one(f)
    from .pari import pari

    err = abs(float(pari()(f"{val.numerator}/{val.denominator}") - x))
    return {"degree": d, "b": [str(x) for x in b], "s_K": s, "counts": counts,
            "siegel": str(val), "pari": str(x), "pari_rational_guess": str(z),
            "abs_error": err, "agree": err < 1e-9 * max(1.0, abs(float(val)))}
