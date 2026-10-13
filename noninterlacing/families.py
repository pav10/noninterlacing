"""The families of Theorems 3.1 and 3.3: ``f_N = m_beta(x^2 - x - N)`` and towers."""

from __future__ import annotations

import math
from typing import Optional

from flint import arb, fmpq, fmpq_poly, fmpz_poly

from .poly import (PolyLike, compose, discriminant, is_irreducible, mahler_measure,
                   poly, real_cyclotomic, to_str)
from .realroots import RootData, is_totally_real


def f_N(N: int, base: Optional[PolyLike] = None) -> fmpz_poly:
    """``m_beta(x^2 - x - N)`` with ``m_beta`` the base (default ``F_21``)."""
    base = real_cyclotomic(21) if base is None else poly(base)
    return compose(base, fmpz_poly([-int(N), -1, 1]))


def tower(base: PolyLike, shifts) -> list:
    """Iterate ``f_{i+1}(x) = f_i(x^2 - x - M_i)``; returns all levels."""
    f = poly(base)
    out = [f]
    for M in shifts:
        f = compose(f, fmpz_poly([-int(M), -1, 1]))
        out.append(f)
    return out


def min_shift_for_total_reality(base: PolyLike) -> int:
    """Least ``N`` with ``1 + 4N + 4 tau(beta) > 0`` for all conjugates."""
    r = RootData(poly(base)).floats()
    N = math.floor(-(1 + 4 * r.min()) / 4) + 1
    while 1 + 4 * N + 4 * r.min() <= 0:
        N += 1
    return N


def min_shift_for_mahler(base: PolyLike) -> int:
    """Least ``N`` making every root of ``f_N`` exceed 1 in absolute value."""
    r = RootData(poly(base)).floats()
    N = min_shift_for_total_reality(base)
    while True:
        disc = 1 + 4 * N + 4 * r
        if disc.min() > 9:  # sqrt > 3  <=>  negative root < -1; positive root > 1 then too
            return N
        N += 1


def norm_value(base: fmpz_poly, t: fmpq) -> fmpq:
    """``N_{F/Q}(t + beta) = (-1)^e m_beta(-t)``."""
    e = base.degree()
    v = fmpq_poly(base)(-t)
    return v if e % 2 == 0 else -v


def check_main_theorem(N: int, base: Optional[PolyLike] = None, interlacing: bool = True) -> dict:
    """Verify items (i)-(v) of Theorem 3.1 for one ``N`` (exact where possible)."""
    base = real_cyclotomic(21) if base is None else poly(base)
    e = base.degree()
    f = f_N(N, base)
    out = {"N": N, "degree": f.degree(), "f": to_str(f)}
    tr = is_totally_real(f)
    out["totally_real"] = tr
    out["irreducible"] = is_irreducible(f)
    if not tr:
        return out
    rd = RootData(f)
    rb = RootData(base).floats()
    span = rd.span
    out["span"] = span
    out["span_formula"] = math.sqrt(1 + 4 * N + 4 * rb.max())
    out["span_ok"] = abs(span - out["span_formula"]) < 1e-9 * max(1.0, span)
    # discriminant: disc(m_beta)^2 * N(1 + 4N + 4 beta) = disc^2 * 4^e N(N + 1/4 + beta)
    D = discriminant(f)
    pred = discriminant(base) ** 2 * (4 ** e) * norm_value(base, fmpq(4 * N + 1, 4))
    out["disc"] = D
    out["disc_formula_ok"] = fmpq(D) == pred
    # Mahler measure
    M = mahler_measure(f, prec=256)
    out["mahler"] = float(M)
    if all(abs(x) > 1 for x in rd.floats()):
        pm = norm_value(base, fmpq(N))
        out["mahler_formula"] = int(pm.p) if pm.q == 1 else pm
        out["mahler_formula_ok"] = bool(M.contains(arb(int(pm.p)) / int(pm.q)))
    else:
        out["mahler_formula"] = None
    if interlacing and out["irreducible"]:
        from .interlace import Interlacing

        I = Interlacing(f)
        res = I.enumerate(1, mode="count")
        out["I(f)"] = res.count
        out["search_nodes"] = res.nodes
        out["noninterlacing"] = res.count == 0 and res.complete
        out["margin"] = I.margin(1).value
    return out
