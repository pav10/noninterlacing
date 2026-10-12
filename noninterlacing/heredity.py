"""Subfields, the trace-down lemma and heredity certificates.

Trace-down (Lemma 2.1 of the paper): for an order ``O`` of ``K``, a subfield
``F`` and ``mu in O^vee`` totally positive, ``nu = Tr_{K/F}(mu)`` lies in
``(O cap F)^vee``, is totally positive and has the same absolute trace.
Hence ``m(O) >= m(O cap F)``.  With ``O = Z[alpha]`` this bounds the
interlacing level of ``f`` from below by data of a smaller field; if
``O_F subset Z[alpha]`` it gives the paper's heredity corollary
``m(O_F) >= 2  =>  f non-interlacing``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from flint import fmpq, fmpq_mat, fmpq_poly, fmpz_poly

from .orders import Order, integer_kernel
from .poly import PolyLike, poly, to_str


def _powers(h: fmpq_poly, f: fmpz_poly, e: int) -> list:
    F = fmpq_poly(f)
    out, p = [], fmpq_poly([1])
    for _ in range(e):
        out.append(p)
        p = (p * h) % F
    return out


def subfield_minpoly(f: PolyLike, h) -> fmpz_poly:
    """Minimal polynomial of ``beta = h(alpha)`` (exact, via a resultant)."""
    from .pari import from_pari, pari, to_pari

    P = pari()
    f = poly(f)
    h = fmpq_poly(h)
    # char poly of multiplication by h, then its squarefree part
    cp = from_pari(P.charpoly(P.Mod(to_pari(h), to_pari(f))))
    g = cp.numer()
    sq = g.factor()[1]
    if len(sq) != 1:
        raise ValueError("h(alpha) does not generate a field (algebra not a field?)")
    m = sq[0][0]
    if int(m.leading_coefficient()) < 0:
        m = -m
    return fmpz_poly([int(c) for c in m.coeffs()])


def intersect_with_subfield(O: Order, h, g: Optional[PolyLike] = None) -> Order:
    """The order ``O cap F`` of ``F = Q(beta)``, ``beta = h(alpha)``, as an
    :class:`Order` of ``Q[y]/(g)`` (``g`` = minimal polynomial of ``beta``)."""
    h = fmpq_poly(h)
    g = poly(g) if g is not None else subfield_minpoly(O.f, h)
    e = g.degree()
    bpows = _powers(h, O.f, e)
    P = [O.coordinates(b) for b in bpows]  # e vectors of length d (columns of P)
    d = O.d
    Pm = fmpq_mat([[P[i][r] for i in range(e)] for r in range(d)])  # d x e
    # left null space of P: rows u with u^T P = 0
    rows = []
    if e < d:
        from flint import fmpz_mat

        den = 1
        for vec in P:
            for x in vec:
                den = den * int(x.q) // _g(den, int(x.q))
        PT = fmpz_mat([[int(x * den) for x in vec] for vec in P])  # e x d
        N, nullity = PT.nullspace()
        for c in range(nullity):
            rows.append([int(N[r, c]) for r in range(d)])
    Lam = integer_kernel(rows, d) if rows else [[int(i == j) for j in range(d)] for i in range(d)]
    # back to beta-coordinates: c = (P^T P)^{-1} P^T v
    PtP_inv = (Pm.transpose() * Pm).inv()
    basis = []
    for v in Lam:
        vv = fmpq_mat([[x] for x in v])
        c = PtP_inv * (Pm.transpose() * vv)
        basis.append(fmpq_poly([c[i, 0] for i in range(e)]))
    return Order.from_lattice(g, basis)


def relative_trace(O: Order, mu, h, g: Optional[PolyLike] = None) -> fmpq_poly:
    """``Tr_{K/F}(mu)`` as a polynomial in ``beta = h(alpha)`` (exact)."""
    h = fmpq_poly(h)
    g = poly(g) if g is not None else subfield_minpoly(O.f, h)
    e = g.degree()
    F = fmpq_poly(O.f)
    mu = fmpq_poly(mu) % F
    bpows = _powers(h, O.f, 2 * e - 1)
    r = fmpq(O.d, e)
    G = fmpq_mat([[O.trace(bpows[i + k]) / r for i in range(e)] for k in range(e)])
    rhs = fmpq_mat([[O.trace((bpows[k] * mu) % F)] for k in range(e)])
    c = G.solve(rhs)
    return fmpq_poly([c[i, 0] for i in range(e)])


def _g(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return abs(a)


# ---------------------------------------------------------------------------
# Certificates
# ---------------------------------------------------------------------------

@dataclass
class SubfieldBound:
    degree: int
    g: fmpz_poly
    h: fmpq_poly
    intersection_disc: int
    field_disc: Optional[int]
    maximal: Optional[bool]
    m: int

    def as_dict(self) -> dict:
        return {
            "degree": self.degree,
            "subfield_poly": to_str(self.g, "y"),
            "embedding": to_str(self.h),
            "disc(O cap F)": self.intersection_disc,
            "disc(F)": self.field_disc,
            "O cap F maximal": self.maximal,
            "m(O cap F)": self.m,
        }


@dataclass
class HeredityReport:
    f: fmpz_poly
    bounds: list = field(default_factory=list)

    @property
    def lower_bound(self) -> int:
        return max([b.m for b in self.bounds] + [1])

    @property
    def certifies_noninterlacing(self) -> bool:
        return self.lower_bound >= 2

    def as_dict(self) -> dict:
        return {
            "f": to_str(self.f),
            "lower bound for m(Z[alpha])": self.lower_bound,
            "subfields": [b.as_dict() for b in self.bounds],
        }


def subfield_bound(O: Order, g: PolyLike, h) -> SubfieldBound:
    g = poly(g)
    Oi = intersect_with_subfield(O, h, g)
    try:
        from .pari import nfdisc

        dF = nfdisc(g)
    except Exception:  # pragma: no cover - PARI missing
        dF = None
    maximal = None if dF is None else (Oi.discriminant == dF)
    return SubfieldBound(g.degree(), g, fmpq_poly(h), Oi.discriminant, dF, maximal, Oi.min_trace())


def heredity_report(f: PolyLike, order: Optional[Order] = None, min_degree: int = 2) -> HeredityReport:
    """Lower bounds ``m(O) >= m(O cap F)`` over all proper subfields (PARI).

    ``order`` defaults to ``Z[alpha]``.  Subfields of degree 2-5 and 7 always
    give ``m = 1`` (Siegel/Yatsyna) and are skipped unless ``min_degree`` is
    lowered *and* ``skip_trivial`` logic is not wanted.
    """
    from .pari import subfields

    f = poly(f)
    O = order or Order.monogenic(f)
    rep = HeredityReport(f)
    d = f.degree()
    for g, h in subfields(f):
        e = g.degree()
        if e < min_degree or e >= d:
            continue
        if e in (2, 3, 4, 5, 7):
            continue  # m(O_F) = 1 for these degrees, hence m(O cap F) = 1
        rep.bounds.append(subfield_bound(O, g, h))
    return rep


def composition_bound(base: PolyLike, inner: PolyLike) -> SubfieldBound:
    """Bound for ``f = base(inner(x))`` from the subfield ``Q(inner(alpha))``.

    No PARI needed for the containment ``Z[beta] subset Z[alpha]``; if
    ``Z[beta]`` is the maximal order of ``F`` (e.g. real cyclotomic bases),
    this is the paper's heredity certificate.
    """
    base, inner = poly(base), poly(inner)
    f = base(inner)
    O = Order.monogenic(f)
    return subfield_bound(O, base, fmpq_poly(inner))


# ---------------------------------------------------------------------------
# Upper bounds by lifting from subfields
# ---------------------------------------------------------------------------

@dataclass
class LiftedWitness:
    """A totally positive element of ``O^vee`` obtained from a subfield.

    If ``nu in O_F^vee`` is totally positive and ``Tr_{K/F}(O) subset O_F``
    (always true for ``O = O_K``), then ``nu in O^vee`` with
    ``Tr_K(nu) = [K:F] Tr_F(nu)``; hence ``m(O) <= [K:F] m(O_F)``.  The
    membership ``nu in O^vee`` is re-checked exactly for the given order.
    """

    subfield_degree: int
    g: fmpz_poly
    m_F: int
    trace: int
    element: fmpq_poly

    def as_dict(self) -> dict:
        return {"subfield_degree": self.subfield_degree, "subfield_poly": to_str(self.g, "y"),
                "m(O_F)": self.m_F, "trace in K": self.trace, "element": to_str(self.element)}


def in_codifferent(O: Order, mu) -> bool:
    F = fmpq_poly(O.f)
    mu = fmpq_poly(mu) % F
    return all(O.trace((mu * w) % F).q == 1 for w in O.basis)


def lifted_witnesses(O: Order, max_index: int = 4) -> list:
    """Witnesses for ``m(O) <= [K:F] m(O_F)`` from proper subfields of index
    ``<= max_index`` (PARI ``nfsubfields`` restricted to those degrees)."""
    from .pari import subfields

    F = fmpq_poly(O.f)
    out = []
    cands = []
    for r in range(2, max_index + 1):
        if O.d % r == 0 and O.d // r >= 1:
            cands += subfields(O.f, O.d // r)
    for g, h in cands:
        e = g.degree()
        if e in (1, O.d):
            continue
        OF = Order.maximal(g)
        mF = OF.min_trace()
        nu = OF.find(mF)
        mu = nu(fmpq_poly(h)) % F
        if not in_codifferent(O, mu) or not O.roots.is_totally_positive(mu):
            continue
        tr = O.trace(mu)
        out.append(LiftedWitness(e, g, mF, int(tr.p), mu))
    out.sort(key=lambda w: w.trace)
    return out


def lift_witness_at(O: Order, j: int, time_limit: float = 30.0) -> Optional[LiftedWitness]:
    """A totally positive element of ``O^vee`` of trace ``j`` lifted from a
    subfield ``F`` of index ``r | j`` (``r >= 2``) from an element of
    ``O_F^vee`` of trace ``j / r``; ``None`` if none is found within budget."""
    from .pari import subfields

    F = fmpq_poly(O.f)
    for r in range(2, j + 1):
        if j % r or O.d % r:
            continue
        for g, h in subfields(O.f, O.d // r):
            OF = Order.maximal(g)
            nu = OF.find(j // r, time_limit=time_limit)
            if nu is None:
                continue
            mu = nu(fmpq_poly(h)) % F
            if in_codifferent(O, mu) and O.roots.is_totally_positive(mu):
                return LiftedWitness(g.degree(), g, j // r, int(O.trace(mu).p), mu)
    return None
