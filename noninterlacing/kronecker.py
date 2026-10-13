"""Interlacing levels of "cyclotomic" polynomials: all roots in [-2, 2].

By Kronecker, a monic ``f in Z[x]`` with all roots in ``[-2, 2]`` is a product
of the polynomials ``Psi_n`` = minimal polynomial of ``2cos(2pi/n)``
(``Psi_1 = x - 2``, ``Psi_2 = x + 2``, ``Psi_3 = x + 1``, ``Psi_4 = x``,
``Psi_6 = x - 1``, ``deg Psi_n = phi(n)/2`` for ``n >= 3``).  We consider the
squarefree ones (distinct factors).

Level 1.  A *monic* interlacer ``g`` of such an ``f`` has its roots strictly
inside ``(-2, 2)``, hence is itself a product of distinct ``Psi_n``.  With
``x = z + 1/z``, ``F(z) = z^d f(x)`` and ``G(z) = (z^2 - 1) z^{d-1} g(x)``
are products of cyclotomic polynomials ``Phi_n`` of the same degree ``2d``
whose zeros interlace on the unit circle -- the situation classified by
Beukers--Heckman (finite hypergeometric groups, Invent. Math. 95 (1989))
and made explicit by McKee--Smyth (Cent. Eur. J. Math. 11 (2013)).
:func:`circle_pair` returns this pair.

Higher levels.  For ``j >= 2`` the interlacer ``g`` has leading coefficient
``j`` and ``g/j`` is not integral in general, so its roots are not Kronecker
numbers and the circle classification does not apply; the levels are
computed with the codifferent machinery of :mod:`noninterlacing.interlace`.
Useful general facts: ``S(f_1 f_2) supseteq S(f_1) + S(f_2)`` for coprime
``f_1, f_2`` (``Z[x]/(f_1 f_2)`` embeds in ``Z[x]/(f_1) x Z[x]/(f_2)``, whose
codifferent is the product), hence ``m(f_1 f_2) <= m(f_1) + m(f_2)``; and
``m(Psi_n) = m(O_{F_n})``.
"""

from __future__ import annotations

import itertools
import time
from dataclasses import dataclass
from typing import Iterator, Optional

from flint import fmpz_poly

from .cyclotomic import euler_phi
from .poly import canonical_form, reflect, to_str


def psi(n: int) -> fmpz_poly:
    """``Psi_n``: minimal polynomial of ``2cos(2pi/n)``."""
    return fmpz_poly.cos_minpoly(n)


def psi_degree(n: int) -> int:
    return 1 if n <= 2 else euler_phi(n) // 2


def psi_indices(max_degree: int) -> list:
    """All ``n`` with ``deg Psi_n <= max_degree`` (``phi(n) >= sqrt(n/2)``)."""
    out = []
    for n in range(1, 8 * max_degree * max_degree + 10):
        if psi_degree(n) <= max_degree:
            out.append(n)
    return out


def reflected_index(n: int) -> int:
    """``Psi_n(-x) = +-Psi_{n'}(x)``: ``n' = 2n`` (n odd), ``n/2`` (n = 2 mod 4), ``n`` (4 | n)."""
    if n % 2:
        return 2 * n
    if n % 4 == 2:
        return n // 2
    return n


def cyclotomic_products(d: int, up_to_reflection: bool = True) -> Iterator[tuple]:
    """Squarefree products of distinct ``Psi_n`` of total degree ``d``.

    Yields sorted index tuples; with ``up_to_reflection`` only one of
    ``f(x)`` and ``(-1)^d f(-x)`` is produced.
    """
    idx = psi_indices(d)
    idx.sort(key=lambda n: (-psi_degree(n), n))

    def rec(start: int, remaining: int, chosen: list):
        if remaining == 0:
            yield tuple(sorted(chosen))
            return
        for i in range(start, len(idx)):
            n = idx[i]
            e = psi_degree(n)
            if e <= remaining:
                chosen.append(n)
                yield from rec(i + 1, remaining - e, chosen)
                chosen.pop()

    seen = set()
    for t in rec(0, d, []):
        if up_to_reflection:
            r = tuple(sorted(reflected_index(n) for n in t))
            key = min(t, r)
            if key in seen:
                continue
            seen.add(key)
            yield key
        else:
            yield t


def product_poly(indices) -> fmpz_poly:
    f = fmpz_poly([1])
    for n in indices:
        f *= psi(n)
    return f


def label(indices) -> str:
    return "*".join(f"Psi_{n}" for n in indices)


def circle_pair(f: fmpz_poly, g: fmpz_poly) -> tuple:
    """``(F, G)`` with ``F(z) = z^d f(z + 1/z)``, ``G(z) = (z^2-1) z^{d-1} g(z + 1/z)``."""
    d = f.degree()

    def lift(p: fmpz_poly, k: int) -> fmpz_poly:
        # z^k p(z + 1/z) as a polynomial in z
        out = fmpz_poly([0])
        zz = fmpz_poly([1, 0, 1])  # z^2 + 1 = z (z + 1/z)
        for i in range(p.degree() + 1):
            c = int(p[i])
            if c:
                out += c * zz ** i * fmpz_poly([0] * (k - i) + [1])
        return out

    return lift(f, d), fmpz_poly([-1, 0, 1]) * lift(g, d - 1)


@dataclass
class KroneckerLevels:
    indices: tuple
    f: fmpz_poly
    degree: int
    m: int
    semigroup: str
    gaps: list
    counts: dict
    margins: dict
    interlacer_factors: Optional[list]

    def as_dict(self) -> dict:
        return {"factors": label(self.indices), "indices": list(self.indices), "f": to_str(self.f),
                "degree": self.degree, "m": self.m, "semigroup": self.semigroup, "gaps": self.gaps,
                "counts": self.counts, "margins": self.margins,
                "monic_interlacer_factors": self.interlacer_factors}


def factor_as_psi(g: fmpz_poly) -> Optional[list]:
    """Indices ``n`` with ``g = prod Psi_n`` (``None`` if ``g`` is not of that form)."""
    _, fac = g.factor()
    out = []
    for p, e in fac:
        if e != 1 or not p.is_cyclotomic() and not _is_psi(p):
            return None
        n = _psi_index(p)
        if n is None:
            return None
        out.append(n)
    return sorted(out)


def _psi_index(p: fmpz_poly) -> Optional[int]:
    e = p.degree()
    for n in psi_indices(e):
        if psi_degree(n) == e and psi(n) == p:
            return n
    return None


def _is_psi(p: fmpz_poly) -> bool:
    return _psi_index(p) is not None


def kronecker_levels(indices, count_levels=(1, 2)) -> KroneckerLevels:
    from .interlace import Interlacing

    f = product_poly(indices)
    I = Interlacing(f, check=False)
    S = I.levels()
    counts = {j: I.count(j) for j in count_levels}
    margins = {j: round(I.margin(j).value, 8) for j in (1, 2)}
    fac = None
    if counts.get(1):
        g = I.find(1)
        fac = factor_as_psi(g)
    return KroneckerLevels(tuple(indices), f, f.degree(), S.multiplicity, str(S), S.gaps, counts,
                           margins, fac)


def kronecker_census(dmax: int, dmin: int = 1, irreducible_only: bool = False,
                     verbose: bool = False) -> list:
    """Level data of every squarefree cyclotomic polynomial of degree ``dmin..dmax``."""
    out = []
    for d in range(dmin, dmax + 1):
        for t in cyclotomic_products(d):
            if irreducible_only and len(t) != 1:
                continue
            r = kronecker_levels(t)
            out.append(r)
            if verbose:
                print(r.degree, label(t), r.semigroup, r.counts, flush=True)
    return out


# ---------------------------------------------------------------------------
# Level 1: the McKee--Smyth / Beukers--Heckman classification
# ---------------------------------------------------------------------------
#
# For f = Psi_1^a Psi_2^b prod_{n in N} Psi_n (squarefree, roots in [-2, 2])
# and a monic g of degree d - 1, put
#     F = Phi_1^a Phi_2^b prod_{n in N} Phi_n,
#     G = Phi_1^(1-a) Phi_2^(1-b) z^(d-1) g(z + 1/z)
# (z^(phi(n)/2) Psi_n(z + 1/z) = Phi_n(z) for n >= 3).  Then g interlaces f on
# [-2, 2] iff {F, G} is an interlacing cyclotomic pair on the unit circle.
# McKee--Smyth, Theorem 1.1: every interlacing pair is {P(z^l), Q(z^l)} with
# {P, Q} or {P(-z), Q(-z)} in their Table 1 (two families) or Table 2 (26
# sporadic pairs).

# Table 2 of McKee--Smyth: label -> (P, Q, P_-, Q_-) as products of Phi_n
MS_SPORADIC = {
    "A": ((12,), (1, 2, 3), (12,), (1, 2, 6)),
    "B": ((3, 12), (1, 2, 8), (6, 12), (1, 2, 8)),
    "C": ((3, 12), (1, 2, 5), (6, 12), (1, 2, 10)),
    "D": ((9,), (1, 2, 4, 6), (18,), (1, 2, 3, 4)),
    "E": ((9,), (1, 2, 8), (18,), (1, 2, 8)),
    "F": ((9,), (1, 2, 5), (18,), (1, 2, 10)),
    "G": ((2, 18), (1, 3, 12), (2, 6, 12), (1, 9)),
    "H": ((2, 18), (1, 3, 5), (2, 6, 10), (1, 9)),
    "I": ((2, 18), (1, 7), (2, 14), (1, 9)),
    "J": ((2, 14), (1, 3, 12), (2, 6, 12), (1, 7)),
    "K": ((2, 14), (1, 3, 5), (2, 6, 10), (1, 7)),
    "L": ((30,), (1, 2, 18), (15,), (1, 2, 9)),
    "M": ((30,), (1, 2, 3, 12), (15,), (1, 2, 6, 12)),
    "N": ((30,), (1, 2, 3, 8), (15,), (1, 2, 6, 8)),
    "O": ((30,), (1, 2, 4, 5), (15,), (1, 2, 4, 10)),
    "P": ((30,), (1, 2, 3, 5), (15,), (1, 2, 6, 10)),
    "Q": ((30,), (1, 2, 7), (15,), (1, 2, 14)),
    "R": ((30,), (1, 2, 4, 12), (15,), (1, 2, 4, 12)),
    "S": ((30,), (1, 2, 4, 8), (15,), (1, 2, 4, 8)),
    "T": ((20,), (1, 2, 3, 12), (20,), (1, 2, 6, 12)),
    "U": ((20,), (1, 2, 3, 8), (20,), (1, 2, 6, 8)),
    "V": ((20,), (1, 2, 7), (20,), (1, 2, 14)),
    "W": ((20,), (1, 2, 9), (20,), (1, 2, 18)),
    "X": ((24,), (1, 2, 4, 5), (24,), (1, 2, 4, 10)),
    "Y": ((24,), (1, 2, 7), (24,), (1, 2, 14)),
    "Z": ((24,), (1, 2, 9), (24,), (1, 2, 18)),
}


def _phi_set(p: fmpz_poly) -> frozenset:
    _, fac = p.factor()
    out = []
    for q, e in fac:
        n = int(q.is_cyclotomic())
        if e != 1 or n == 0:
            raise ValueError("not a squarefree product of cyclotomic polynomials")
        out.append(n)
    return frozenset(out)


def _phi_prod(ns) -> fmpz_poly:
    p = fmpz_poly([1])
    for n in ns:
        p *= fmpz_poly.cyclotomic(n)
    return p


def _minus(ns: frozenset) -> frozenset:
    return frozenset(reflected_index(n) for n in ns)


def _dilate(ns: frozenset, l: int) -> frozenset:
    p = fmpz_poly([1])
    for n in ns:
        p *= fmpz_poly.cyclotomic(n)(fmpz_poly([0] * l + [1]))
    return _phi_set(p)


def ms_primitive_pairs(max_degree: int) -> set:
    """Primitive interlacing cyclotomic pairs of degree <= ``max_degree``
    (as frozensets of two frozensets of cyclotomic indices), closed under
    ``z -> -z``."""
    from math import gcd

    z = fmpz_poly([0, 1])
    one = fmpz_poly([1])
    pairs = set()

    def add(P: fmpz_poly, Q: fmpz_poly):
        a, b = _phi_set(P), _phi_set(Q)
        pairs.add(frozenset((a, b)))
        pairs.add(frozenset((_minus(a), _minus(b))))

    for n in range(1, max_degree + 1):  # family 1
        N = n + 1
        P = (z ** N - one) // (z - one)
        for j in range(1, N // 2 + 1):
            if gcd(j, N) == 1:
                add(P, (z ** j - one) * (z ** (N - j) - one) // (z - one))
    for n in range(3, max_degree + 1):  # family 2
        P = z ** n + one
        for j in range(1, n):
            if gcd(j, 2 * n) == 1:
                add(P, (z ** j - one) * (z ** (n - j) + one))
    for P, Q, Pm, Qm in MS_SPORADIC.values():
        if _phi_prod(P).degree() <= max_degree:
            add(_phi_prod(P), _phi_prod(Q))
            add(_phi_prod(Pm), _phi_prod(Qm))
    return pairs


def ms_all_pairs(max_degree: int) -> set:
    """All interlacing cyclotomic pairs of degree <= ``max_degree``."""
    out = set()
    for pair in ms_primitive_pairs(max_degree):
        a, b = tuple(pair)
        e = _phi_prod(a).degree()
        for l in range(1, max_degree // e + 1):
            out.add(frozenset((_dilate(a, l), _dilate(b, l))))
    return out


def circle_index_set(indices) -> frozenset:
    """``f = prod Psi_n  ->  F = prod Phi_n`` (the same index set)."""
    return frozenset(indices)


def ms_monic_interlacers(indices, pairs: set) -> list:
    """Monic interlacers of ``prod Psi_n`` predicted by McKee--Smyth, as sorted index tuples."""
    N = circle_index_set(indices)
    out = []
    for pair in pairs:
        if N in pair:
            other = next(iter(pair - {N})) if len(pair) == 2 else N
            g = tuple(sorted(n for n in other if n >= 3 or n not in (1, 2)))
            # remove the Phi_1/Phi_2 factors that belong to G for missing roots +-2
            g = tuple(sorted(n for n in other if not (n in (1, 2) and n not in N)))
            out.append(g)
    return sorted(out)


# ---------------------------------------------------------------------------
# Higher levels: a lower bound from the arithmetic--geometric mean inequality
# ---------------------------------------------------------------------------

def amgm_lower_bound(f: fmpz_poly) -> tuple:
    """``m(f) >= ceil( sum_i e_i |disc(f_i) Res(f_i, f/f_i)|^(-1/e_i) )``.

    For a level-``j`` interlacer ``g`` and an irreducible factor ``f_i`` of
    degree ``e_i`` with root ``alpha_i``, ``mu_i = g(alpha_i)/f'(alpha_i)`` is
    totally positive in ``Q(alpha_i)`` and ``N(g(alpha_i))`` is a non-zero
    integer, so ``Tr mu_i >= e_i N(mu_i)^(1/e_i) >= e_i |N f'(alpha_i)|^(-1/e_i)``
    with ``N f'(alpha_i) = +- disc(f_i) Res(f_i, f/f_i)``; and ``j = sum_i Tr mu_i``.
    For irreducible ``f`` this is Yatsyna's ``|disc f| >= d^d`` for interlacing ``f``.
    Returns ``(bound, real_sum)``.
    """
    import math

    _, fac = f.factor()
    total = 0.0
    for fi, _e in fac:
        e = fi.degree()
        rest = f // fi
        D = abs(int(fi.discriminant())) * (abs(int(fi.resultant(rest))) if rest.degree() > 0 else 1)
        total += e * D ** (-1.0 / e) if D else math.inf
    return max(1, math.ceil(total - 1e-12)), total


def circle_lift(indices, g: fmpz_poly) -> tuple:
    """``(F, G)`` on the unit circle for ``f = prod Psi_n`` and an interlacer ``g``.

    ``F = prod Phi_n`` and ``G = Phi_1^(1-a) Phi_2^(1-b) z^(d-1) g(z + 1/z)``
    where ``a, b`` record the roots ``2, -2`` of ``f``.  For a level-``j``
    interlacer, ``F`` and ``G`` interlace on the circle and ``H = F + G`` has
    leading coefficient ``j + 1`` and all zeros in ``|z| < 1``; ``H`` is
    2-reciprocal and ``F = (H + H*)/2`` (``(H - H*)/2`` if ``f(2) = 0``).  For ``j = 1`` this is the
    "disc-bionic" polynomial of McKee--Smyth (Theorem 2.1).
    """
    N = set(indices)
    d = sum(psi_degree(n) for n in N)
    F = _phi_prod(sorted(N))
    zz = fmpz_poly([1, 0, 1])
    k = d - 1
    lift = fmpz_poly([0])
    for i in range(g.degree() + 1):
        c = int(g[i])
        if c:
            lift += c * zz ** i * fmpz_poly([0] * (k - i) + [1])
    G = lift
    for m in (1, 2):
        if m not in N:
            G *= fmpz_poly.cyclotomic(m)
    return F, G


def check_disc_polynomial(indices, g: fmpz_poly) -> dict:
    import numpy as np

    F, G = circle_lift(indices, g)
    H = F + G
    r = np.roots([float(c) for c in reversed(H.coeffs())])
    Hstar = fmpz_poly(list(reversed(H.coeffs())) + [0] * 0)
    two_reciprocal = all((int(a) - int(b)) % 2 == 0 for a, b in
                         zip(H.coeffs(), list(reversed(H.coeffs()))))
    Hs = _reverse(H, H.degree())
    # Phi_1 = z - 1 is anti-reciprocal: F* = (-1)^a F, G* = -(-1)^a G
    a = 1 if 1 in set(indices) else 0
    rec = (H + Hs) if a == 0 else (H - Hs)
    return {"deg": H.degree(), "lead": int(H.leading_coefficient()), "max|root|": float(np.max(np.abs(r))),
            "two_reciprocal": two_reciprocal, "F_recovered": rec == 2 * F}


def _reverse(p: fmpz_poly, n: int) -> fmpz_poly:
    c = [int(x) for x in p.coeffs()] + [0] * (n + 1 - len(p.coeffs()))
    return fmpz_poly(list(reversed(c)))


# ---------------------------------------------------------------------------
# Decomposition along unit resultants
# ---------------------------------------------------------------------------

def is_prime_power(n: int) -> bool:
    from .cyclotomic import factorint

    return n > 1 and len(factorint(n)) == 1


def linked(m: int, n: int) -> bool:
    """``Res(Psi_m, Psi_n) != +-1``  iff  ``m/n`` or ``n/m`` is a prime power."""
    a, b = max(m, n), min(m, n)
    return a % b == 0 and is_prime_power(a // b)


def components(indices) -> list:
    """Connected components of the graph ``m ~ n  iff  m/n is a prime power``.

    If ``f = f_1 f_2`` with ``Res(f_1, f_2) = +-1`` then
    ``Z[x]/(f) = Z[x]/(f_1) x Z[x]/(f_2)`` (CRT), so the totally positive
    codifferent elements -- i.e. the interlacers -- split:
    every level-``j`` interlacer of ``f`` is uniquely ``g_1 f_2 + g_2 f_1`` with
    ``g_i`` interlacers of ``f_i`` of levels adding to ``j``.  Hence
    ``S(f) = S(f_1) + S(f_2)``, ``m(f) = m(f_1) + m(f_2)`` and
    ``I_j(f) = sum_{j_1 + j_2 = j} I_{j_1}(f_1) I_{j_2}(f_2)``.
    """
    idx = list(indices)
    comp, seen = [], set()
    for n in idx:
        if n in seen:
            continue
        stack, cur = [n], []
        seen.add(n)
        while stack:
            a = stack.pop()
            cur.append(a)
            for b in idx:
                if b not in seen and linked(a, b):
                    seen.add(b)
                    stack.append(b)
        comp.append(tuple(sorted(cur)))
    return sorted(comp)


def convolve_counts(*count_lists) -> list:
    """Convolution of interlacer count sequences ``[I_0 = 0, I_1, I_2, ...]``."""
    out = [1]
    for c in count_lists:
        new = [0] * (len(out) + len(c) - 1)
        for i, a in enumerate(out):
            for k, b in enumerate(c):
                new[i + k] += a * b
        out = new
    return out
