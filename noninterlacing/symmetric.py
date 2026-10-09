"""Integer symmetric matrices and interlacing (Lemma 4.1 of the paper).

If ``f`` is irreducible and the minimal polynomial of an integer symmetric
``n x n`` matrix ``A`` (so ``charpoly(A) = f^m``, ``n = dm``), then for the
principal submatrix ``A'`` (last row/column deleted)

    h = charpoly(A') / f^(m-1)

is a monic integer interlacer of ``f``.  Hence non-interlacing polynomials
are never minimal polynomials of integer symmetric matrices (Estes--Guralnick
counterexamples).  :func:`cauchy_interlacer` makes this constructive;
:func:`find_symmetric_realization` is a (non-exhaustive) search for an
integer symmetric matrix with characteristic polynomial ``f``, built by
the same interlacing recursion.
"""

from __future__ import annotations

import itertools
import random
from typing import Optional, Sequence

from flint import fmpz_mat, fmpz_poly

from .poly import PolyLike, poly, to_str
from .realroots import interlaces


def charpoly(A: Sequence[Sequence[int]]) -> fmpz_poly:
    return fmpz_mat([list(map(int, r)) for r in A]).charpoly()


def minimal_polynomial(A: Sequence[Sequence[int]]) -> fmpz_poly:
    return fmpz_mat([list(map(int, r)) for r in A]).minpoly()


def cauchy_interlacer(A: Sequence[Sequence[int]], delete: Optional[int] = None) -> dict:
    """The interlacer of Lemma 4.1 obtained from an integer symmetric matrix."""
    A = [list(map(int, r)) for r in A]
    n = len(A)
    if any(A[i][j] != A[j][i] for i in range(n) for j in range(n)):
        raise ValueError("matrix is not symmetric")
    chi = charpoly(A)
    _, fac = chi.factor()
    if len(fac) != 1:
        raise ValueError("characteristic polynomial is not a power of an irreducible polynomial")
    f, m = fac[0]
    i = n - 1 if delete is None else delete
    Ap = [[A[r][c] for c in range(n) if c != i] for r in range(n) if r != i]
    chi_p = charpoly(Ap) if Ap else fmpz_poly([1])
    q, rem = divmod(chi_p, f ** (m - 1))
    if not rem.is_zero():  # pragma: no cover - contradicts Cauchy interlacing
        raise AssertionError("f^(m-1) does not divide charpoly(A')")
    return {"f": f, "multiplicity": m, "h": q, "interlaces": interlaces(f, q) if f.degree() > 1 else True}


def random_symmetric(n: int, bound: int = 2, rng: Optional[random.Random] = None) -> list:
    rng = rng or random.Random()
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            A[i][j] = A[j][i] = rng.randint(-bound, bound)
    return A


def find_symmetric_realization(f: PolyLike, max_nodes: int = 20000,
                               entry_bound: Optional[int] = None) -> Optional[list]:
    """Search an integer symmetric ``d x d`` matrix with characteristic polynomial ``f``.

    Recursion: such a matrix is ``[[A', b], [b^T, a]]`` with ``A'`` a
    realisation of a monic interlacer ``g`` of ``f``, ``a = tr f - tr g``,
    and ``b`` an integer vector with ``|b|^2 = (sum r_i^2 - sum s_i^2 - a^2)/2``.
    The search tries interlacers and vectors ``b`` of the right norm in turn;
    it is *not* exhaustive (realisations of ``g`` other than the first found
    are not tried), so ``None`` does not prove non-realisability -- for that
    use the interlacing test (Lemma 4.1).
    """
    from .interlace import Interlacing

    f = poly(f)
    budget = [max_nodes]

    def p2(h: fmpz_poly) -> int:  # sum of squares of roots
        e = h.degree()
        s1 = -int(h[e - 1]) if e >= 1 else 0
        s2 = int(h[e - 2]) if e >= 2 else 0
        return s1 * s1 - 2 * s2

    def rec(h: fmpz_poly) -> Optional[list]:
        budget[0] -= 1
        if budget[0] < 0:
            return None
        e = h.degree()
        if e == 1:
            return [[-int(h[0])]]
        try:
            gs = Interlacing(h, check=False).interlacers(1, limit=200)
        except ValueError:
            return None
        for g in gs:
            a = -int(h[e - 1]) + int(g[e - 2])
            nb2 = p2(h) - p2(g) - a * a
            if nb2 < 0 or nb2 % 2:
                continue
            nb2 //= 2
            Ap = rec(g)
            if Ap is None:
                continue
            for b in _vectors_of_norm(e - 1, nb2):
                M = [row + [b[i]] for i, row in enumerate(Ap)] + [list(b) + [a]]
                if charpoly(M) == h:
                    return M
        return None

    return rec(f)


def _vectors_of_norm(n: int, N: int):
    """Integer vectors in ``Z^n`` with ``|v|^2 = N`` (small cases)."""
    if n == 0:
        if N == 0:
            yield []
        return
    from math import isqrt

    r = isqrt(N)
    for x in range(-r, r + 1):
        for rest in _vectors_of_norm(n - 1, N - x * x):
            yield [x] + rest
