"""Small exact lattice helpers (unimodular completions, LLL on real bases)."""

from __future__ import annotations

from typing import Sequence

from flint import arb, fmpz_mat


def egcd(a: int, b: int) -> tuple[int, int, int]:
    """Return ``(g, x, y)`` with ``g = gcd(a, b) >= 0`` and ``a x + b y = g``."""
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q, a, b = a // b, b, a - (a // b) * b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    if a < 0:
        a, x0, y0 = -a, -x0, -y0
    return a, x0, y0


def unimodular_completion(t: Sequence[int]) -> tuple[int, list[list[int]]]:
    """Unimodular ``W`` (as a list of columns) with ``t . W = (g, 0, ..., 0)``.

    ``g = gcd(t)``.  Column ``0`` of ``W`` is a vector of ``t``-value ``g``;
    the remaining columns form a basis of the kernel of ``t`` on ``Z^d``.
    """
    d = len(t)
    v = [int(a) for a in t]
    cols = [[1 if r == c else 0 for r in range(d)] for c in range(d)]
    # move a nonzero entry to position 0
    for c in range(d):
        if v[c] != 0:
            v[0], v[c] = v[c], v[0]
            cols[0], cols[c] = cols[c], cols[0]
            break
    for c in range(1, d):
        if v[c] == 0:
            continue
        g, x, y = egcd(v[0], v[c])
        a, b = v[0] // g, v[c] // g
        c0 = [x * p + y * q for p, q in zip(cols[0], cols[c])]
        cc = [-b * p + a * q for p, q in zip(cols[0], cols[c])]
        cols[0], cols[c] = c0, cc
        v[0], v[c] = g, 0
    if v[0] < 0:
        v[0] = -v[0]
        cols[0] = [-p for p in cols[0]]
    return v[0], cols


def arb_to_scaled_int(x: arb, shift: int) -> int:
    """``round(mid(x) * 2^shift)`` as a Python integer."""
    m, e = x.mid().man_exp()
    m, e = int(m), int(e) + shift
    if e >= 0:
        return m << e
    k = -e
    return (m + (1 << (k - 1))) >> k


def lll_columns(cols: Sequence[Sequence[arb]], shift: int) -> list[list[int]]:
    """LLL-reduce real column vectors; return the unimodular transform ``U``.

    ``cols[c][r]`` is coordinate ``r`` of basis vector ``c``.  The vectors are
    scaled by ``2^shift`` and rounded to integers before calling FLINT's LLL,
    so ``shift`` must be large enough that rounding errors are negligible
    relative to the shortest vectors.  The returned ``U`` (list of rows,
    ``n x n``) satisfies: new basis vector ``k`` = ``sum_c U[c][k] * cols[c]``.
    Correctness of downstream computations never depends on the quality of
    the reduction, only on ``U`` being unimodular, which it always is.
    """
    n = len(cols)
    if n == 0:
        return []
    rows = [[arb_to_scaled_int(x, shift) for x in col] for col in cols]
    M = fmpz_mat(rows)
    _, T = M.lll(transform=True)  # T * M = reduced, rows are new vectors
    T = [[int(T[i, j]) for j in range(n)] for i in range(n)]
    # new vector i = sum_c T[i][c] * old vector c  ->  U[c][i] = T[i][c]
    return [[T[i][c] for i in range(n)] for c in range(n)]


def mat_mul_int(A: Sequence[Sequence[int]], B: Sequence[Sequence[int]]) -> list[list[int]]:
    """Product of integer matrices given as lists of rows."""
    Bt = list(zip(*B))
    return [[sum(a * b for a, b in zip(row, col)) for col in Bt] for row in A]
