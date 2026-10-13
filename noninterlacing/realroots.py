"""Certified real roots and exact sign/interlacing tests.

Two independent rigorous tools are provided:

* :class:`RootData` -- isolating balls (``arb``) for the real roots of a
  squarefree, real-rooted integer polynomial, computed with FLINT/Arb and
  refined on demand.  Signs of other polynomials at the roots are decided
  by ball arithmetic (with a gcd test for exact zeros).
* Sturm / Cauchy index machinery in exact rational arithmetic
  (:func:`cauchy_index`, :func:`tarski_query`, :func:`interlaces`), which
  never touches floating point at all.
"""

from __future__ import annotations

from contextlib import contextmanager
from functools import lru_cache
from typing import Iterable, Sequence

import numpy as np
from flint import arb, arb_poly, ctx, fmpq, fmpq_poly, fmpz_poly

from .poly import PolyLike, poly


class NotTotallyReal(ValueError):
    """Raised when a polynomial is not squarefree with only real roots."""


@contextmanager
def precision(bits: int):
    old = ctx.prec
    ctx.prec = bits
    try:
        yield
    finally:
        ctx.prec = old


def _as_fmpq_poly(p) -> fmpq_poly:
    if isinstance(p, fmpq_poly):
        return p
    return fmpq_poly(poly(p))


# ---------------------------------------------------------------------------
# Exact Sturm machinery
# ---------------------------------------------------------------------------

def _primitive_positive(p: fmpq_poly) -> fmpq_poly:
    """Positive rational multiple of ``p`` with coprime integer coefficients."""
    if p.is_zero():
        return p
    num = p.numer()  # fmpz_poly, p = num / den with den > 0
    c = int(num.content())
    return fmpq_poly([int(a) // c for a in num.coeffs()])


def signed_remainder_sequence(p: PolyLike, q: PolyLike) -> list[fmpq_poly]:
    """``[p, q, -rem(p, q), ...]`` up to positive scalar factors."""
    a, b = _as_fmpq_poly(p), _as_fmpq_poly(q)
    seq = [_primitive_positive(a)]
    if b.is_zero():
        return seq
    seq.append(_primitive_positive(b))
    while True:
        r = seq[-2] % seq[-1]
        if r.is_zero():
            return seq
        seq.append(_primitive_positive(-r))


def _sign(x) -> int:
    return (x > 0) - (x < 0)


def _variations(signs: Iterable[int]) -> int:
    s = [x for x in signs if x != 0]
    return sum(1 for u, v in zip(s, s[1:]) if u != v)


def _var_at_infinity(seq: Sequence[fmpq_poly], positive: bool) -> int:
    signs = []
    for p in seq:
        lc = _sign(p[p.degree()])
        if not positive and p.degree() % 2:
            lc = -lc
        signs.append(lc)
    return _variations(signs)


def _var_at(seq: Sequence[fmpq_poly], x: fmpq) -> int:
    return _variations(_sign(p(x)) for p in seq)


def cauchy_index(q: PolyLike, p: PolyLike) -> int:
    """Cauchy index of ``q/p`` over the real line (exact).

    Each real pole where ``q/p`` jumps from ``-oo`` to ``+oo`` counts ``+1``,
    the opposite jump ``-1``.  For ``p`` squarefree with simple real roots
    ``r`` and ``deg q < deg p`` this is ``sum_r sign(q(r)/p'(r))``.
    """
    seq = signed_remainder_sequence(p, q)
    return _var_at_infinity(seq, False) - _var_at_infinity(seq, True)


def sturm_count(f: PolyLike, a=None, b=None) -> int:
    """Number of distinct real roots of ``f`` (in ``(a, b]`` if given)."""
    f = _as_fmpq_poly(f)
    seq = signed_remainder_sequence(f, f.derivative())
    va = _var_at_infinity(seq, False) if a is None else _var_at(seq, fmpq(a))
    vb = _var_at_infinity(seq, True) if b is None else _var_at(seq, fmpq(b))
    return va - vb


def tarski_query(q: PolyLike, f: PolyLike) -> int:
    """``sum_{f(r)=0, r real} sign q(r)`` (exact, Sturm--Tarski)."""
    f = _as_fmpq_poly(f)
    q = _as_fmpq_poly(q)
    return cauchy_index((f.derivative() * q) % f, f)


def is_totally_real(f: PolyLike) -> bool:
    """``f`` squarefree with all roots real (exact)."""
    f = poly(f)
    d = f.degree()
    if d < 1:
        return False
    if f.gcd(f.derivative()).degree() > 0:
        return False
    return sturm_count(f) == d


def is_totally_positive_exact(q: PolyLike, f: PolyLike) -> bool:
    """Whether ``q(r) > 0`` at every root ``r`` of ``f`` (exact)."""
    f = poly(f)
    return tarski_query(q, f) == f.degree()


def interlaces(f: PolyLike, g: PolyLike) -> bool:
    """Exact test that the roots of ``g`` strictly interlace those of ``f``.

    True iff ``f`` has ``d`` distinct real roots, ``deg g = d - 1`` and every
    open gap between consecutive roots of ``f`` contains exactly one root of
    ``g`` (so ``g`` has no root in common with ``f``).  The sign of the
    leading coefficient of ``g`` is irrelevant here; with a positive
    leading coefficient ``j`` the residues ``g(r)/f'(r)`` are all positive
    and ``g(alpha)/f'(alpha)`` is a totally positive codifferent element of
    trace ``j``.
    """
    f, g = poly(f), poly(g)
    d = f.degree()
    if d < 1 or g.degree() != d - 1:
        return False
    return abs(cauchy_index(g, f)) == d


# ---------------------------------------------------------------------------
# Certified numerical roots
# ---------------------------------------------------------------------------

class RootData:
    """Certified real roots of a squarefree real-rooted integer polynomial.

    ``balls(prec)`` returns isolating ``arb`` balls in increasing order with
    relative accuracy at least ``prec`` bits.
    """

    def __init__(self, f: PolyLike, check: bool = True):
        self.f = poly(f)
        self.d = self.f.degree()
        if self.d < 1:
            raise NotTotallyReal("polynomial must have positive degree")
        if check and not is_totally_real(self.f):
            raise NotTotallyReal(f"{self.f} is not squarefree with all roots real")
        self.fprime = self.f.derivative()
        self._balls: dict[int, list[arb]] = {}

    def balls(self, prec: int = 128) -> list[arb]:
        if prec not in self._balls:
            with precision(prec):
                roots = self.f.complex_roots()
            out = []
            for z, mult in roots:
                if mult != 1 or not z.imag.is_zero():
                    raise NotTotallyReal(f"{self.f} has a non-real or repeated root")
                out.append(z.real)
            out.sort(key=lambda b: float(b.mid()))
            self._balls[prec] = out
        return self._balls[prec]

    def floats(self) -> np.ndarray:
        return np.array([float(b.mid()) for b in self.balls(64)])

    @property
    def span(self) -> float:
        r = self.floats()
        return float(r[-1] - r[0])

    def _zero_indices(self, q: fmpz_poly) -> set[int]:
        g = self.f.gcd(q)
        if g.degree() <= 0:
            return set()
        # The roots of g are roots of f; locate them by ball overlap.
        prec = 128
        while True:
            fb = self.balls(prec)
            with precision(prec):
                gb = [z.real for z, _ in g.complex_roots()]
            hits = set()
            ok = True
            for b in gb:
                idx = [i for i, a in enumerate(fb) if a.overlaps(b)]
                if len(idx) != 1:
                    ok = False
                    break
                hits.add(idx[0])
            if ok and len(hits) == g.degree():
                return hits
            prec *= 2

    def signs(self, q: PolyLike) -> list[int]:
        """Exact signs of ``q`` at the roots (``0`` exactly where ``q`` vanishes)."""
        if isinstance(q, fmpq_poly):
            q = q.numer()  # positive denominator does not change signs
        q = poly(q)
        zeros = self._zero_indices(q)
        signs = [0] * self.d
        todo = [i for i in range(self.d) if i not in zeros]
        prec = 128
        coeffs = [int(c) for c in q.coeffs()]
        while todo:
            fb = self.balls(prec)
            with precision(prec):
                qa = arb_poly(coeffs)
                rest = []
                for i in todo:
                    v = qa(fb[i])
                    if v > 0:
                        signs[i] = 1
                    elif v < 0:
                        signs[i] = -1
                    else:
                        rest.append(i)
            todo = rest
            prec *= 2
            if prec > 1 << 20:  # pragma: no cover - would indicate a bug
                raise RuntimeError("sign determination did not converge")
        return signs

    def residue_signs(self, g: PolyLike) -> list[int]:
        """Signs of ``g(r) / f'(r)``; all ``+1`` iff ``g`` is a positive interlacer."""
        s = self.signs(g)
        n = self.d
        # sign f'(r_i) = (-1)^(d-1-i) for increasing roots r_0 < ... < r_{d-1}
        return [s[i] * (1 if (n - 1 - i) % 2 == 0 else -1) for i in range(n)]

    def is_totally_positive(self, q: PolyLike) -> bool:
        return all(s > 0 for s in self.signs(q))

    def is_positive_interlacer(self, g: PolyLike) -> bool:
        """``g`` strictly interlaces ``f`` and has positive leading coefficient."""
        g = poly(g)
        if g.degree() != self.d - 1 or int(g.leading_coefficient()) <= 0:
            return False
        return all(s > 0 for s in self.residue_signs(g))


@lru_cache(maxsize=256)
def _rootdata_cached(coeffs: tuple) -> RootData:
    return RootData(fmpz_poly(list(coeffs)))


def root_data(f: PolyLike) -> RootData:
    """Cached :class:`RootData` for ``f``."""
    f = poly(f)
    return _rootdata_cached(tuple(int(c) for c in f.coeffs()))
