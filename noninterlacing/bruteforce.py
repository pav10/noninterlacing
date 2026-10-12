"""Naive reference method: the coefficient bounding box of the interlacing polytope.

For a level ``j`` the interlacers ``g = j x^{d-1} + c_{d-2} x^{d-2} + ... + c_0``
form the integer points of the polytope

    P_j(f) = { c in R^{d-1} : g(r_i) / f'(r_i) > 0  for every root r_i }.

This module computes the exact (LP) bounding box of ``P_j(f)`` and checks
every integer point of the box.  It is exponentially slower than
:mod:`noninterlacing.simplex` but completely independent of it, and is used
for cross-validation (and to reproduce the "bounding box" count quoted in
the paper).
"""

from __future__ import annotations

import itertools
import math

import numpy as np

from .poly import PolyLike, poly
from .realroots import RootData, interlaces


def interlacing_polytope(f: PolyLike, level: int = 1):
    """Return ``(A, b)`` with ``P_j(f) = { c : A c < b }`` (floating point)."""
    rd = RootData(poly(f))
    r = rd.floats()
    d = rd.d
    fp = np.array([float(rd.fprime(float(x))) for x in r])
    w = 1.0 / fp
    V = np.vander(r, d, increasing=True)  # columns r^0 .. r^{d-1}
    # mu_i = w_i (j r_i^{d-1} + sum_k c_k r_i^k) > 0  <=>  -w_i V_i[:d-1] c < j w_i r_i^{d-1}
    A = -(w[:, None] * V[:, : d - 1])
    b = level * w * V[:, d - 1]
    return A, b


def bounding_box(f: PolyLike, level: int = 1):
    """Exact-LP bounding box ``[(lo_k, hi_k)]`` of ``P_j(f)``, ``k = 0..d-2``."""
    from scipy.optimize import linprog

    A, b = interlacing_polytope(f, level)
    n = A.shape[1]
    box = []
    for k in range(n):
        c = np.zeros(n)
        c[k] = 1.0
        lo = linprog(c, A_ub=A, b_ub=b, bounds=[(None, None)] * n, method="highs")
        hi = linprog(-c, A_ub=A, b_ub=b, bounds=[(None, None)] * n, method="highs")
        if lo.status != 0 or hi.status != 0:
            raise RuntimeError("LP failed")
        box.append((lo.fun, -hi.fun))
    return box


def integer_box(f: PolyLike, level: int = 1):
    """Integer ranges of the bounding box and the number of lattice points in it."""
    box = bounding_box(f, level)
    ranges = [(math.ceil(lo - 1e-9), math.floor(hi + 1e-9)) for lo, hi in box]
    size = 1
    for lo, hi in ranges:
        size *= max(0, hi - lo + 1)
    return ranges, size


def brute_force_interlacers(f: PolyLike, level: int = 1, max_box: int = 5_000_000):
    """All level-``j`` interlacers by scanning the integer bounding box.

    A vectorised floating-point filter (with a generous margin) discards
    points that are clearly outside; every survivor is checked exactly with
    the Sturm/Cauchy-index test :func:`noninterlacing.realroots.interlaces`.
    """
    f = poly(f)
    d = f.degree()
    ranges, size = integer_box(f, level)
    if size > max_box:
        raise ValueError(f"bounding box has {size} points (> max_box={max_box})")
    A, b = interlacing_polytope(f, level)
    if size == 0:
        return [], size
    axes = [np.arange(lo, hi + 1) for lo, hi in ranges]
    found = []
    # iterate over all but the first coordinate in Python, vectorise the first
    for rest in itertools.product(*axes[1:]):
        rest = np.array(rest, dtype=float)
        base = A[:, 1:] @ rest if len(rest) else np.zeros(A.shape[0])
        c0 = axes[0].astype(float)
        vals = base[None, :] + c0[:, None] * A[:, 0][None, :]
        ok = np.all(vals < b[None, :] + 1e-7 * (1 + np.abs(b[None, :])), axis=1)
        for idx in np.nonzero(ok)[0]:
            coeffs = [int(axes[0][idx])] + [int(x) for x in rest] + [level]
            from flint import fmpz_poly

            g = fmpz_poly(coeffs)
            if interlaces(f, g):
                found.append(g)
    found.sort(key=lambda g: [int(c) for c in reversed(g.coeffs())])
    return found, size


def lp_search_interlacers(f: PolyLike, level: int = 1, max_nodes: int = 10_000_000):
    """Second reference method: depth-first search over the raw coefficients.

    Coordinates ``c_{d-2}, c_{d-3}, ..., c_0`` are fixed one at a time; the
    admissible range of the next one is the exact LP range over the slice of
    ``P_j(f)`` (HiGHS), widened by a margin.  Leaves are decided by the exact
    Sturm test.  No lattice reduction, convex hulls or embeddings are used,
    so this is independent of :mod:`noninterlacing.simplex`.
    Returns ``(interlacers, nodes)``.
    """
    from flint import fmpz_poly
    from scipy.optimize import linprog

    f = poly(f)
    d = f.degree()
    A, b = interlacing_polytope(f, level)  # A c < b, c = (c_0, ..., c_{d-2})
    n = d - 1
    order = list(range(n - 1, -1, -1))  # fix c_{d-2} first
    found = []
    nodes = [0]
    fixed = {}

    def bounds_for(k):
        free = [i for i in range(n) if i not in fixed]
        rhs = b - sum(A[:, i] * v for i, v in fixed.items()) if fixed else b.copy()
        Af = A[:, free]
        pos = free.index(k)
        c = np.zeros(len(free))
        c[pos] = 1.0
        lo = linprog(c, A_ub=Af, b_ub=rhs, bounds=[(None, None)] * len(free), method="highs")
        if lo.status != 0:
            return None
        hi = linprog(-c, A_ub=Af, b_ub=rhs, bounds=[(None, None)] * len(free), method="highs")
        return lo.fun, -hi.fun

    def rec(depth):
        nodes[0] += 1
        if nodes[0] > max_nodes:
            raise RuntimeError("node limit exceeded")
        k = order[depth]
        if depth == n - 1:
            rhs = b - sum(A[:, i] * v for i, v in fixed.items())
            a = A[:, k]
            lo, hi = -np.inf, np.inf
            for ai, ri in zip(a, rhs):
                if ai > 0:
                    hi = min(hi, ri / ai)
                elif ai < 0:
                    lo = max(lo, ri / ai)
                elif ri <= -1e-9:
                    return
            for v in range(math.ceil(lo - 1e-6), math.floor(hi + 1e-6) + 1):
                coeffs = [0] * (d + 1)
                coeffs[d - 1] = level
                for i, w in fixed.items():
                    coeffs[i] = w
                coeffs[k] = v
                g = fmpz_poly(coeffs)
                if interlaces(f, g):
                    found.append(g)
            return
        r = bounds_for(k)
        if r is None:
            return
        lo, hi = r
        for v in range(math.ceil(lo - 1e-6), math.floor(hi + 1e-6) + 1):
            fixed[k] = v
            rec(depth + 1)
            del fixed[k]

    rec(0)
    found.sort(key=lambda g: [int(c) for c in reversed(g.coeffs())])
    return found, nodes[0]
