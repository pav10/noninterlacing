"""Optional PARI/GP access through :mod:`cypari2`.

PARI is used for maximal orders, subfields, Dedekind zeta values and ideal
factorisation.  Everything else in the package works without it.
"""

from __future__ import annotations

from functools import lru_cache

from flint import fmpq, fmpq_poly, fmpz_poly

_PARI = None


class PariUnavailable(RuntimeError):
    pass


def pari():
    """The shared :class:`cypari2.Pari` instance (raises if not installed)."""
    global _PARI
    if _PARI is None:
        try:
            import cypari2
        except ImportError as exc:  # pragma: no cover
            raise PariUnavailable(
                "this feature needs PARI/GP: pip install cypari2"
            ) from exc
        _PARI = cypari2.Pari()
        try:
            _PARI.allocatemem(1 << 30, silent=True)
        except Exception:  # pragma: no cover
            _PARI.allocatemem(1 << 28)
        _PARI.set_real_precision(60)
    return _PARI


def have_pari() -> bool:
    try:
        pari()
        return True
    except PariUnavailable:
        return False


def to_pari(p) -> "object":
    """``fmpz_poly``/``fmpq_poly`` -> PARI polynomial in ``x``."""
    P = pari()
    if isinstance(p, fmpz_poly):
        coeffs = [int(c) for c in p.coeffs()]
        return P.Pol(list(reversed(coeffs)) or [0])
    if isinstance(p, fmpq_poly):
        coeffs = [P(f"{int(c.p)}/{int(c.q)}") for c in p.coeffs()]
        return P.Pol(list(reversed(coeffs)) or [0])
    raise TypeError(type(p))


def from_pari(p) -> fmpq_poly:
    """PARI polynomial (or scalar, or polmod) with rational coefficients -> ``fmpq_poly``."""
    P = pari()
    t = str(P.type(p))
    if t == "t_POLMOD":
        p = P.lift(p)
        t = str(P.type(p))
    if t == "t_VEC" and len(p) == 1:
        p = p[0]
        t = str(P.type(p))
    if t in ("t_INT", "t_FRAC"):
        return fmpq_poly([fmpq(int(P.numerator(p)), int(P.denominator(p)))])
    out = []
    for c in P.Vecrev(p):
        out.append(fmpq(int(P.numerator(c)), int(P.denominator(c))))
    return fmpq_poly(out)


def nfbasis(f: fmpz_poly) -> list:
    """Integral basis of ``Q[x]/(f)`` as ``fmpq_poly`` (PARI ``nfbasis``)."""
    return [from_pari(b) for b in pari().nfbasis(to_pari(f))]


def nfdisc(f: fmpz_poly) -> int:
    return int(pari().nfdisc(to_pari(f)))


def subfields(f: fmpz_poly, degree: int | None = None) -> list:
    """Subfields of ``Q[x]/(f)``: list of ``(g, h)`` with ``g(h(alpha)) = 0``.

    ``g`` is the defining polynomial (``fmpz_poly``) of the subfield and ``h``
    (``fmpq_poly``) expresses a root of ``g`` in terms of ``alpha``.
    """
    P = pari()
    res = P.nfsubfields(to_pari(f), degree) if degree else P.nfsubfields(to_pari(f))
    out = []
    for pair in res:
        g = from_pari(pair[0])
        out.append((fmpz_poly([int(c) for c in g.numer().coeffs()]), from_pari(pair[1])))
    return out


def polredabs(f: fmpz_poly) -> fmpz_poly:
    g = from_pari(pari().polredabs(to_pari(f)))
    return fmpz_poly([int(c) for c in g.numer().coeffs()])


@lru_cache(maxsize=64)
def _nfinit_cached(coeffs: tuple):
    return pari().nfinit(to_pari(fmpz_poly(list(coeffs))))


def nfinit(f: fmpz_poly):
    return _nfinit_cached(tuple(int(c) for c in f.coeffs()))


def zeta_at_minus_one(f: fmpz_poly, bitprec: int = 256):
    """``zeta_K(-1)`` for ``K = Q[x]/(f)`` via PARI ``lfun``.

    Returns ``(q, x)``: ``x`` the numerical value (``bitprec`` bits) and ``q``
    its best rational approximation with denominator below ``10^12``
    (``zeta_K(-1)`` is rational; the guess is only used for display, the
    comparisons in :mod:`noninterlacing.siegel` are done numerically).
    """
    P = pari()
    old = P.get_real_precision_bits()
    P.set_real_precision_bits(bitprec)
    try:
        val = P.lfun(P.nfinit(to_pari(f)), -1)
    finally:
        P.set_real_precision_bits(old)
    r = P.bestappr(val, 10 ** 12)
    return fmpq(int(P.numerator(r)), int(P.denominator(r))), val
