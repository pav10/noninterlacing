"""Orders of totally real etale algebras and their codifferents.

An :class:`Order` is a full-rank subring ``O`` of ``A = Q[x]/(f)`` with ``f``
monic, squarefree and real-rooted (``A`` is a product of totally real
fields; for irreducible ``f`` it is the field ``K = Q(alpha)``).  It is given
by a ``Z``-basis of rational polynomials reduced modulo ``f``.

For each order we can compute its codifferent ``O^vee = {x : Tr(xO) in Z}``
and Yatsyna's invariant

    m(O) = min { Tr(mu) : mu in O^vee totally positive },

together with the numerical semigroup of all traces of totally positive
codifferent elements, the elements of a given trace, and the *margin*
``max { min_i sigma_i(mu) : mu in O^vee, Tr mu = j }``.

``Order.monogenic(f)`` is ``Z[alpha]`` (whose totally positive codifferent
elements of trace ``j`` are exactly the level-``j`` interlacers of ``f``,
:mod:`noninterlacing.interlace`); ``Order.maximal(f)`` is ``O_K`` (PARI).
Monotonicity ``O_1 subset O_2  =>  m(O_1) <= m(O_2)`` and the trace-down
inequality ``m(O) >= m(O cap F)`` (subfield ``F``) are implemented in
:mod:`noninterlacing.heredity`.
"""

from __future__ import annotations

from fractions import Fraction
from functools import cached_property
from typing import Iterable, Optional, Sequence

from flint import arb, arb_poly, fmpq, fmpq_mat, fmpq_poly, fmpz_mat, fmpz_poly

from .lattice import lll_columns
from .poly import PolyLike, poly, to_str
from .realroots import RootData, precision
from .semigroup import NumericalSemigroup
from .simplex import SimplexLattice


def power_sums(f: fmpz_poly, count: int) -> list[int]:
    """Exact power sums ``p_k = sum r^k`` (``k = 0..count-1``) of the roots of monic ``f``."""
    d = f.degree()
    a = [int(f[d - i]) for i in range(d + 1)]  # a[0] = 1, f = sum a[i] x^(d-i)
    p = [d]
    for k in range(1, count):
        s = k * a[k] if k <= d else 0
        for i in range(1, min(k, d + 1)):
            s += a[i] * p[k - i]
        p.append(-s)
    return p


def integer_kernel(rows: Sequence[Sequence[int]], ncols: int) -> list[list[int]]:
    """Z-basis (list of vectors) of ``{v in Z^n : A v = 0}`` for integer ``A``.

    Unimodular column reduction (extended gcd), so the result is a basis of
    the full kernel lattice, not of a finite-index sublattice.
    """
    n = ncols
    U = [[int(i == j) for j in range(n)] for i in range(n)]  # columns of U are U[c]
    A = [list(map(int, r)) for r in rows]
    piv = 0
    for r in range(len(A)):
        if piv >= n:
            break
        row = A[r]
        for c in range(piv + 1, n):
            if row[c] == 0:
                continue
            a, b = row[piv], row[c]
            if a == 0:
                for R in A:
                    R[piv], R[c] = R[c], R[piv]
                U[piv], U[c] = U[c], U[piv]
                continue
            from .lattice import egcd

            g, x, y = egcd(a, b)
            p, q = a // g, b // g
            for R in A:
                u, v = R[piv], R[c]
                R[piv], R[c] = x * u + y * v, -q * u + p * v
            cu, cv = U[piv], U[c]
            U[piv] = [x * s + y * t for s, t in zip(cu, cv)]
            U[c] = [-q * s + p * t for s, t in zip(cu, cv)]
        if A[r][piv] != 0:
            piv += 1
    return [list(U[c]) for c in range(piv, n)]


def _lll_reduce_int_vectors(vecs: list[list[int]]) -> list[list[int]]:
    if not vecs:
        return vecs
    M = fmpz_mat(vecs)
    L = M.lll()
    return [[int(L[i, j]) for j in range(L.ncols())] for i in range(L.nrows())]


class Order:
    """An order of ``Q[x]/(f)`` with an explicit Z-basis."""

    def __init__(self, f: PolyLike, basis: Iterable, name: str = ""):
        self.f = poly(f)
        if int(self.f.leading_coefficient()) != 1:
            raise ValueError("defining polynomial must be monic")
        self.d = self.f.degree()
        F = fmpq_poly(self.f)
        self.basis = [fmpq_poly(b) % F for b in basis]
        if len(self.basis) != self.d:
            raise ValueError("an order needs exactly deg f basis elements")
        self.name = name
        self.roots = RootData(self.f)
        if self.discriminant == 0:
            raise ValueError("basis elements are linearly dependent")

    # -------------------------------------------------------------- creation
    @classmethod
    def monogenic(cls, f: PolyLike) -> "Order":
        f = poly(f)
        return cls(f, [fmpq_poly([0] * k + [1]) for k in range(f.degree())], name="Z[a]")

    @classmethod
    def maximal(cls, f: PolyLike) -> "Order":
        """The maximal order (PARI ``nfbasis``; ``f`` irreducible)."""
        from .pari import nfbasis

        f = poly(f)
        return cls(f, nfbasis(f), name="O_K")

    @classmethod
    def from_lattice(cls, f: PolyLike, elements: Iterable) -> "Order":
        """The Z-module spanned by ``elements`` (must have full rank)."""
        f = poly(f)
        F = fmpq_poly(f)
        els = [fmpq_poly(e) % F for e in elements]
        return cls(f, _zbasis(els, f.degree()))

    @classmethod
    def generated_by(cls, f: PolyLike, generators: Iterable) -> "Order":
        """The ring ``Z[generators]`` (closure under multiplication)."""
        f = poly(f)
        F = fmpq_poly(f)
        d = f.degree()
        gens = [fmpq_poly(g) % F for g in generators]
        span = [fmpq_poly([1])]
        while True:
            new = span + [(a * g) % F for a in span for g in gens]
            basis = _zbasis(new, d, allow_deficient=True)
            if len(basis) == len(span) and _same_lattice(basis, span, d):
                break
            span = basis
        if len(span) != d:
            raise ValueError("generators do not generate an order of full rank")
        return cls(f, span)

    # -------------------------------------------------------------- algebra
    def coords_matrix(self) -> fmpq_mat:
        """Columns: power-basis coordinates of the basis elements."""
        d = self.d
        rows = [[self.basis[c][r] for c in range(d)] for r in range(d)]
        return fmpq_mat(rows)

    @cached_property
    def _coords_inv(self) -> fmpq_mat:
        return self.coords_matrix().inv()

    def coordinates(self, x: PolyLike | fmpq_poly) -> list:
        """Coordinates (rationals) of ``x`` in the basis of the order."""
        x = fmpq_poly(x) % fmpq_poly(self.f)
        v = fmpq_mat([[x[r]] for r in range(self.d)])
        w = self._coords_inv * v
        return [w[i, 0] for i in range(self.d)]

    def contains(self, x) -> bool:
        return all(c.q == 1 for c in self.coordinates(x))

    def is_ring(self) -> bool:
        F = fmpq_poly(self.f)
        if not self.contains(fmpq_poly([1])):
            return False
        return all(self.contains((a * b) % F) for a in self.basis for b in self.basis)

    def contains_order(self, other: "Order") -> bool:
        return all(self.contains(b) for b in other.basis)

    def index_in(self, other: "Order") -> Fraction:
        """``[other : self]`` (for ``self subset other``)."""
        D1, D2 = self.discriminant, other.discriminant
        r = Fraction(D1, D2)
        from math import isqrt

        num, den = r.numerator, r.denominator
        sn, sd = isqrt(num), isqrt(den)
        if sn * sn != num or sd * sd != den:
            raise ValueError("discriminant ratio is not a square")
        return Fraction(sn, sd)

    @cached_property
    def _power_sums(self) -> list[int]:
        return power_sums(self.f, self.d)

    def trace(self, x) -> fmpq:
        x = fmpq_poly(x) % fmpq_poly(self.f)
        p = self._power_sums
        return sum((x[k] * p[k] for k in range(self.d)), fmpq(0))

    def trace_form(self) -> fmpq_mat:
        F = fmpq_poly(self.f)
        return fmpq_mat([[self.trace((a * b) % F) for b in self.basis] for a in self.basis])

    @cached_property
    def discriminant(self) -> int:
        D = self.trace_form().det()
        if D.q != 1:
            raise ValueError("basis does not consist of algebraic integers")
        return int(D.p)

    @cached_property
    def codifferent(self) -> list:
        """Trace-dual basis ``w_k^vee`` (``Tr(w_a^vee w_b) = delta_ab``)."""
        Ginv = self.trace_form().inv()
        d = self.d
        out = []
        for a in range(d):
            e = fmpq_poly([0])
            for b in range(d):
                if Ginv[a, b] != 0:
                    e += Ginv[a, b] * self.basis[b]
            out.append(e)
        return out

    @cached_property
    def codifferent_traces(self) -> list[int]:
        tr = [self.trace(w) for w in self.codifferent]
        if any(t.q != 1 for t in tr):
            raise ValueError("1 is not in the order: codifferent traces are not integral")
        return [int(t.p) for t in tr]

    def element(self, coords: Sequence[int], dual: bool = True) -> fmpq_poly:
        basis = self.codifferent if dual else self.basis
        e = fmpq_poly([0])
        for c, w in zip(coords, basis):
            if c:
                e += int(c) * w
        return e

    # ------------------------------------------------------ codifferent lattice
    def _embedding(self, prec: int):
        balls = self.roots.balls(prec + 32)
        out = []
        with precision(prec + 32):
            polys = []
            for w in self.codifferent:
                num = [int(c) for c in w.numer().coeffs()]
                polys.append((arb_poly(num), int(w.denom())))
            for r in balls:
                out.append([P(r) / q for P, q in polys])
        return out

    def _verify(self, coords) -> bool:
        return self.roots.is_totally_positive(self.element(coords))

    @cached_property
    def lattice(self) -> SimplexLattice:
        return SimplexLattice(self._embedding, self.codifferent_traces, verify=self._verify)

    def totally_positive_codifferent(self, trace: int, limit: Optional[int] = None) -> list:
        """All ``mu in O^vee``, ``mu >> 0`` with ``Tr(mu) = trace`` (``fmpq_poly`` in ``alpha``)."""
        res = self.lattice.enumerate(trace, mode="list", limit=limit)
        return [self.element(a) for a in res.points]

    def count(self, trace: int) -> int:
        res = self.lattice.enumerate(trace, mode="count")
        return res.count

    def find(self, trace: int, method: str = "auto"):
        if method == "auto":
            method = "dfs" if self.d <= 20 else "milp"
        if method == "milp":
            mm = self.lattice.max_min(trace)
            if mm.verified:
                return self.element(mm.point)
            if mm.exists is False:
                return None
        res = self.lattice.enumerate(trace, mode="first")
        return self.element(res.points[0]) if res.points else None

    def margin(self, trace: int = 1):
        return self.lattice.max_min(trace)

    def min_trace(self, method: str = "auto") -> int:
        """Yatsyna's ``m(O)``."""
        j = 1
        while self.find(j, method) is None:
            j += 1
        return j

    def trace_semigroup(self, method: str = "auto") -> NumericalSemigroup:
        m = self.min_trace(method)
        return NumericalSemigroup.from_oracle(lambda j: self.find(j, method) is not None, m)

    def __repr__(self) -> str:
        nm = f" {self.name}" if self.name else ""
        return f"<Order{nm} of Q[x]/({to_str(self.f)}), disc {self.discriminant}>"


def _zbasis(elements: list, d: int, allow_deficient: bool = False) -> list:
    """Z-basis of the lattice spanned by rational polynomials (degree < d)."""
    if not elements:
        return []
    den = 1
    for e in elements:
        den = den * int(e.denom()) // _gcd(den, int(e.denom()))
    vecs = []
    for e in elements:
        v = [int(e[k] * den) for k in range(d)]
        vecs.append(v)
    M = fmpz_mat(vecs)
    H = M.hnf()
    rows = [[int(H[i, j]) for j in range(d)] for i in range(H.nrows())]
    rows = [r for r in rows if any(r)]
    if len(rows) != d and not allow_deficient:
        raise ValueError("elements do not span a full-rank lattice")
    rows = _lll_reduce_int_vectors(rows)
    return [fmpq_poly([fmpq(c, den) for c in r]) for r in rows]


def _same_lattice(a: list, b: list, d: int) -> bool:
    if len(a) != len(b):
        return False
    den = 1
    for e in a + b:
        den = den * int(e.denom()) // _gcd(den, int(e.denom()))
    A = fmpz_mat([[int(e[k] * den) for k in range(d)] for e in a]).hnf()
    B = fmpz_mat([[int(e[k] * den) for k in range(d)] for e in b]).hnf()
    return A == B


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return abs(a)
