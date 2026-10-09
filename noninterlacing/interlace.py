"""Interlacing polynomials of every level.

Definitions
-----------
Let ``f in Z[x]`` be monic, squarefree, of degree ``d`` with all roots real,
``r_1 > r_2 > ... > r_d``.  For an integer ``j >= 1`` call ``g in Z[x]`` a
**level-j interlacer** (or *j-interlacer*) of ``f`` if ``deg g = d - 1``, the
leading coefficient of ``g`` is ``j`` and the roots ``s_1 > ... > s_{d-1}`` of
``g`` satisfy ``r_i > s_i > r_{i+1}``.  ``f`` is **j-interlaced** if a level-j
interlacer exists; *interlacing* means 1-interlaced.

Codifferent dictionary (Yatsyna)
--------------------------------
With ``alpha`` a root of ``f`` (in the etale algebra ``Q[x]/(f)``),
``Z[alpha]^vee = f'(alpha)^{-1} Z[alpha]`` (Euler), and

    g  |->  mu = g(alpha) / f'(alpha)

is a bijection from level-j interlacers onto totally positive elements of
``Z[alpha]^vee`` of trace ``j``: the trace is the ``x^{d-1}`` coefficient of
``g`` (Euler's identity ``sum g(r)/f'(r) = [x^{d-1}] g``), and positivity of
all ``g(r)/f'(r)`` is equivalent to ``g`` alternating in sign along the roots,
i.e. to strict interlacing.  Consequently

* the set of levels ``S(f) = { j : f is j-interlaced }`` is closed under
  addition (sum of interlacers is an interlacer), contains ``d`` (``g = f'``)
  and ``1 + d*N`` for large ``N``; hence it is a **numerical semigroup**:
  ``f`` fails to be ``j``-interlaced for only finitely many ``j`` (the gaps);
* ``m(f) := min S(f) = m(Z[alpha])`` is Yatsyna's minimal trace, and ``f``
  is non-interlacing iff ``m(f) >= 2``;
* ``m(f) <= m(O_K) <= dim M_{2d}(SL_2(Z))`` (Siegel), see
  :mod:`noninterlacing.siegel`.

The counting problem is a lattice-point problem in a regular simplex,
solved by :class:`noninterlacing.simplex.SimplexLattice`.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Iterable, Optional

from flint import arb, fmpz_poly

from .poly import PolyLike, chebyshev_C, poly, to_str
from .realroots import RootData, interlaces, precision
from .semigroup import NumericalSemigroup
from .simplex import EnumerationResult, SimplexLattice


def siegel_ceiling(d: int) -> int:
    """``dim M_{2d}(SL_2(Z))``: an upper bound for ``m(O_K)`` in degree ``d``."""
    w = 2 * d
    if w == 2:
        return 0
    return w // 12 + (0 if w % 12 == 2 else 1)


class Interlacing:
    """All interlacing data of a monic totally real squarefree ``f``.

    >>> I = Interlacing("x^6 - x^5 - 6*x^4 + 6*x^3 + 8*x^2 - 8*x + 1")
    >>> I.count(1), I.min_level()
    (0, 2)
    """

    def __init__(self, f: PolyLike, check: bool = True):
        self.f = poly(f)
        if int(self.f.leading_coefficient()) != 1:
            raise ValueError("f must be monic")
        self.d = self.f.degree()
        self.roots = RootData(self.f, check=check)
        fl = self.roots.floats()
        # shift of the power basis (x - c)^k improves conditioning; Z[x-c] = Z[x]
        self.c = int(round(float(fl.mean()))) if self.d else 0
        shift = fmpz_poly([-self.c, 1])
        self._basis = [shift ** k for k in range(self.d)]
        self.lattice = SimplexLattice(
            self._embedding,
            [0] * (self.d - 1) + [1],
            verify=self._verify,
        )

    # -- lattice plumbing ---------------------------------------------------
    def _embedding(self, prec: int):
        balls = self.roots.balls(prec + 32)
        fp = self.roots.fprime
        with precision(prec + 32):
            out = []
            for r in balls:
                inv = 1 / fp(r)
                z = r - self.c
                row, p = [], arb(1)
                for _ in range(self.d):
                    row.append(p * inv)
                    p = p * z
                out.append(row)
        return out

    def coords_to_poly(self, a: Iterable[int]) -> fmpz_poly:
        g = fmpz_poly([0])
        for ak, b in zip(a, self._basis):
            if ak:
                g += int(ak) * b
        return g

    def _verify(self, a) -> bool:
        return self.roots.is_positive_interlacer(self.coords_to_poly(a))

    # -- public API ---------------------------------------------------------
    def enumerate(self, level: int = 1, mode: str = "list", **kw) -> EnumerationResult:
        """Raw enumeration result (points are coordinate vectors)."""
        return self.lattice.enumerate(level, mode=mode, **kw)

    def interlacers(self, level: int = 1, limit: Optional[int] = None, **kw) -> list:
        """All level-``level`` interlacers as ``fmpz_poly`` (sorted)."""
        res = self.lattice.enumerate(level, mode="list", limit=limit, **kw)
        gs = [self.coords_to_poly(a) for a in res.points]
        gs.sort(key=lambda g: [int(c) for c in reversed(g.coeffs())])
        return gs

    def count(self, level: int = 1, **kw) -> int:
        """Number ``I_level(f)`` of level-``level`` interlacers."""
        res = self.lattice.enumerate(level, mode="count", **kw)
        if not res.complete:
            raise RuntimeError("enumeration aborted before completion")
        return res.count

    def margin(self, level: int = 1, **kw):
        """Interlacing margin ``eps_j(f) = max_g min_i g(r_i)/f'(r_i)`` (MILP).

        The maximum runs over all ``g in Z[x]`` of degree ``<= d-1`` with
        ``x^{d-1}``-coefficient ``j``; ``f`` is ``j``-interlaced iff the margin
        is positive.  Returns a :class:`~noninterlacing.simplex.MaxMinResult`
        whose ``point`` converts to the optimal ``g`` via :meth:`coords_to_poly`.
        """
        return self.lattice.max_min(level, **kw)

    def find(self, level: int = 1, method: str = "auto", **kw) -> Optional[fmpz_poly]:
        """Some level-``level`` interlacer, or ``None`` if there is none.

        ``method``: ``"dfs"`` (complete enumeration, stops at the first
        point), ``"milp"`` (HiGHS branch and bound on the margin; any
        interlacer it returns is verified exactly, a negative margin is a
        floating-point proof of non-existence) or ``"auto"`` (DFS up to degree
        20, MILP above, DFS again if the MILP is inconclusive).
        """
        if method == "auto":
            method = "dfs" if self.d <= 20 else "milp"
        if method == "milp":
            mm = self.lattice.max_min(level)
            if mm.verified:
                return self.coords_to_poly(mm.point)
            if mm.exists is False:
                return None
            method = "dfs"
        res = self.lattice.enumerate(level, mode="first", **kw)
        if res.points:
            return self.coords_to_poly(res.points[0])
        if not res.complete:
            raise RuntimeError("enumeration aborted before completion")
        return None

    def is_interlaced(self, level: int = 1, **kw) -> bool:
        return self.find(level, **kw) is not None

    def min_level(self, method: str = "auto") -> int:
        """``m(f) = m(Z[alpha])``: the least level ``j`` with a ``j``-interlacer."""
        j = 1
        while True:
            if self.find(j, method=method) is not None:
                return j
            j += 1
            if j > self.d:  # f' is always a level-d interlacer
                raise AssertionError("f' should interlace f")  # pragma: no cover

    def levels(self, method: str = "auto") -> NumericalSemigroup:
        """The numerical semigroup ``S(f)`` of interlacing levels."""
        m = self.min_level(method)
        return NumericalSemigroup.from_oracle(lambda j: self.find(j, method=method) is not None, m)

    def certificate(self, level: int = 1) -> dict:
        """A checkable witness: an interlacer, or the fact that none exists."""
        g = self.find(level)
        out = {"f": to_str(self.f), "level": level, "interlaced": g is not None}
        if g is not None:
            out["interlacer"] = to_str(g)
            out["exact_check"] = interlaces(self.f, g) and int(g.leading_coefficient()) == level
        return out


@lru_cache(maxsize=128)
def _cached(coeffs: tuple) -> Interlacing:
    return Interlacing(fmpz_poly(list(coeffs)))


def interlacing(f: PolyLike) -> Interlacing:
    f = poly(f)
    return _cached(tuple(int(c) for c in f.coeffs()))


def interlacers(f: PolyLike, level: int = 1, limit: Optional[int] = None) -> list:
    """All level-``level`` interlacers of ``f``."""
    return interlacing(f).interlacers(level, limit=limit)


def count_interlacers(f: PolyLike, level: int = 1) -> int:
    return interlacing(f).count(level)


def find_interlacer(f: PolyLike, level: int = 1) -> Optional[fmpz_poly]:
    return interlacing(f).find(level)


def is_interlacing(f: PolyLike, level: int = 1) -> bool:
    return interlacing(f).is_interlaced(level)


def min_level(f: PolyLike) -> int:
    return interlacing(f).min_level()


def level_semigroup(f: PolyLike) -> NumericalSemigroup:
    return interlacing(f).levels()
