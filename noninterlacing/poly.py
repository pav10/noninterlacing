"""Integer polynomial utilities.

Polynomials are represented internally as :class:`flint.fmpz_poly`.  Every
public function accepts any *polynomial-like* object:

* an ``fmpz_poly``;
* a string such as ``"x^6 - x^5 - 6*x^4 + 6*x^3 + 8*x^2 - 8*x + 1"``
  (``^`` or ``**`` for powers, implicit multiplication ``2x`` allowed,
  parentheses and composition such as ``"(x^2-x-3)^2 + 1"`` allowed);
* a list/tuple of integers in **descending** degree order, e.g.
  ``[1, -1, -6, 6, 8, -8, 1]`` (the convention used in print and in PARI's
  ``Pol``).
"""

from __future__ import annotations

import math
import re
from typing import Sequence, Union

from flint import arb, ctx, fmpq_poly, fmpz, fmpz_poly

PolyLike = Union[fmpz_poly, fmpq_poly, str, Sequence[int], int]

X = fmpz_poly([0, 1])


# ---------------------------------------------------------------------------
# Parsing and printing
# ---------------------------------------------------------------------------

_TOKEN = re.compile(r"\s*(?:(\d+)|([A-Za-z_]\w*)|(\*\*|[-+*/^()]))")


class _Parser:
    """Tiny recursive-descent parser for integer polynomial expressions."""

    def __init__(self, text: str, var: str | None):
        self.tokens = []
        pos = 0
        text = text.strip()
        while pos < len(text):
            m = _TOKEN.match(text, pos)
            if not m or m.end() == pos:
                raise ValueError(f"cannot parse polynomial near {text[pos:]!r}")
            num, name, op = m.groups()
            if num is not None:
                self.tokens.append(("num", int(num)))
            elif name is not None:
                self.tokens.append(("var", name))
            else:
                self.tokens.append(("op", "^" if op == "**" else op))
            pos = m.end()
        self.i = 0
        self.var = var

    def peek(self):
        return self.tokens[self.i] if self.i < len(self.tokens) else (None, None)

    def take(self):
        tok = self.peek()
        self.i += 1
        return tok

    def parse(self) -> fmpq_poly:
        value = self.expr()
        if self.i != len(self.tokens):
            raise ValueError(f"unexpected token {self.peek()[1]!r}")
        return value

    def expr(self) -> fmpq_poly:
        value = self.term()
        while self.peek() in (("op", "+"), ("op", "-")):
            op = self.take()[1]
            rhs = self.term()
            value = value + rhs if op == "+" else value - rhs
        return value

    def term(self) -> fmpq_poly:
        value = self.unary()
        while True:
            kind, tok = self.peek()
            if (kind, tok) == ("op", "*"):
                self.take()
                value = value * self.unary()
            elif (kind, tok) == ("op", "/"):
                self.take()
                den = self.unary()
                if den.degree() != 0:
                    raise ValueError("division only by constants is supported")
                value = value / den[0]
            elif kind in ("num", "var") or (kind, tok) == ("op", "("):
                value = value * self.power()  # implicit multiplication
            else:
                return value

    def unary(self) -> fmpq_poly:
        if self.peek() == ("op", "-"):
            self.take()
            return -self.unary()
        if self.peek() == ("op", "+"):
            self.take()
            return self.unary()
        return self.power()

    def power(self) -> fmpq_poly:
        base = self.atom()
        if self.peek() == ("op", "^"):
            self.take()
            exp = self.unary()
            if exp.degree() > 0 or exp[0].denom() != 1 or exp[0] < 0:
                raise ValueError("exponents must be non-negative integers")
            return base ** int(exp[0])
        return base

    def atom(self) -> fmpq_poly:
        kind, tok = self.take()
        if kind == "num":
            return fmpq_poly([tok])
        if kind == "var":
            if self.var is None:
                self.var = tok
            elif tok != self.var:
                raise ValueError(f"polynomial uses two variables: {self.var!r}, {tok!r}")
            return fmpq_poly([0, 1])
        if (kind, tok) == ("op", "("):
            value = self.expr()
            if self.take() != ("op", ")"):
                raise ValueError("unbalanced parentheses")
            return value
        raise ValueError(f"unexpected token {tok!r}")


def parse_poly(text: str, var: str | None = None) -> fmpz_poly:
    """Parse an integer polynomial from a string."""
    q = _Parser(text, var).parse()
    if q.denom() != 1:
        raise ValueError(f"polynomial {text!r} does not have integer coefficients")
    return fmpz_poly([int(c) for c in q.numer().coeffs()])


def poly(p: PolyLike) -> fmpz_poly:
    """Coerce a polynomial-like object to :class:`fmpz_poly`."""
    if isinstance(p, fmpz_poly):
        return p
    if isinstance(p, fmpq_poly):
        if p.denom() != 1:
            raise ValueError("polynomial does not have integer coefficients")
        return fmpz_poly([int(c) for c in p.numer().coeffs()])
    if isinstance(p, str):
        return parse_poly(p)
    if isinstance(p, (int, fmpz)):
        return fmpz_poly([int(p)])
    if isinstance(p, (list, tuple)):
        return fmpz_poly([int(c) for c in reversed(p)])
    try:  # numpy arrays and other sequences
        return fmpz_poly([int(c) for c in reversed(list(p))])
    except TypeError as exc:  # pragma: no cover - defensive
        raise TypeError(f"cannot interpret {p!r} as a polynomial") from exc


def coeffs_desc(p: PolyLike) -> list[int]:
    """Coefficients in descending degree order."""
    return [int(c) for c in reversed(poly(p).coeffs())]


def coeffs_asc(p: PolyLike) -> list[int]:
    """Coefficients in ascending degree order."""
    return [int(c) for c in poly(p).coeffs()]


def to_str(p: PolyLike, var: str = "x") -> str:
    """Human readable form, e.g. ``x^6 - x^5 - 6*x^4 + 1``."""
    if isinstance(p, fmpq_poly) and p.denom() != 1:
        den = int(p.denom())
        return f"({to_str(p.numer(), var)})/{den}"
    p = poly(p)
    if p.is_zero():
        return "0"
    parts = []
    for k in range(p.degree(), -1, -1):
        c = int(p[k])
        if c == 0:
            continue
        sign = "-" if c < 0 else "+"
        a = abs(c)
        if k == 0:
            body = str(a)
        else:
            mon = var if k == 1 else f"{var}^{k}"
            body = mon if a == 1 else f"{a}*{mon}"
        parts.append((sign, body))
    first_sign, first_body = parts[0]
    out = ("-" if first_sign == "-" else "") + first_body
    for sign, body in parts[1:]:
        out += f" {sign} {body}"
    return out


# ---------------------------------------------------------------------------
# Basic operations
# ---------------------------------------------------------------------------

def degree(p: PolyLike) -> int:
    return poly(p).degree()


def is_monic(p: PolyLike) -> bool:
    p = poly(p)
    return p.degree() >= 0 and int(p.leading_coefficient()) == 1


def compose(f: PolyLike, g: PolyLike) -> fmpz_poly:
    """Return ``f(g(x))``."""
    return poly(f)(poly(g))


def translate(f: PolyLike, n: int) -> fmpz_poly:
    """Return ``f(x + n)``."""
    return poly(f)(fmpz_poly([int(n), 1]))


def reflect(f: PolyLike) -> fmpz_poly:
    """Return ``(-1)^d f(-x)`` (monic if ``f`` is): roots ``r -> -r``."""
    f = poly(f)
    d = f.degree()
    g = f(fmpz_poly([0, -1]))
    return -g if d % 2 else g


def chebyshev_C(k: int) -> fmpz_poly:
    """Monic Chebyshev polynomial ``C_k(x) = 2 T_k(x/2)``, ``C_k(2cos t) = 2cos(kt)``."""
    if k < 0:
        raise ValueError("k must be non-negative")
    a, b = fmpz_poly([2]), X
    if k == 0:
        return a
    for _ in range(k - 1):
        a, b = b, X * b - a
    return b


def real_cyclotomic(n: int) -> fmpz_poly:
    """Minimal polynomial of ``beta_n = 2cos(2 pi / n) = zeta_n + zeta_n^{-1}``."""
    if n < 1:
        raise ValueError("n must be positive")
    return fmpz_poly.cos_minpoly(n)


def is_squarefree(f: PolyLike) -> bool:
    f = poly(f)
    return f.gcd(f.derivative()).degree() == 0


def is_irreducible(f: PolyLike) -> bool:
    f = poly(f)
    if f.degree() < 1:
        return False
    content, factors = f.factor()
    return abs(int(content)) == 1 and len(factors) == 1 and factors[0][1] == 1


def factor(f: PolyLike) -> list[tuple[fmpz_poly, int]]:
    return poly(f).factor()[1]


def discriminant(f: PolyLike) -> int:
    return int(poly(f).discriminant())


def resultant(f: PolyLike, g: PolyLike) -> int:
    return int(poly(f).resultant(poly(g)))


def content_free(p: fmpz_poly) -> fmpz_poly:
    """Divide by the (positive) content."""
    c = int(p.content())
    if c in (0, 1):
        return p
    return fmpz_poly([int(a) // c for a in p.coeffs()])


# ---------------------------------------------------------------------------
# Translation classes
# ---------------------------------------------------------------------------

def translation_normal_form(f: PolyLike) -> tuple[fmpz_poly, int]:
    """Translate ``f`` so that the root mean lies in ``[0, 1)``.

    Returns ``(f(x + n), n)``.  Two monic polynomials are equivalent under
    ``x -> x + n`` iff their normal forms coincide.
    """
    f = poly(f)
    d = f.degree()
    if d < 1:
        return f, 0
    a = int(f[d - 1])  # = -(sum of roots) for monic f
    n = (-a) // d
    return translate(f, n), n


def canonical_form(f: PolyLike, reflections: bool = True) -> fmpz_poly:
    """Canonical representative modulo ``x -> x + n`` (and ``x -> -x``)."""
    g, _ = translation_normal_form(f)
    if not reflections:
        return g
    h, _ = translation_normal_form(reflect(f))
    return min(g, h, key=coeffs_desc)


def z_equivalent(f: PolyLike, g: PolyLike, reflections: bool = False) -> bool:
    """Whether ``g(x) = f(x + n)`` (or ``(-1)^d f(-x + n)``) for some integer ``n``."""
    return canonical_form(f, reflections) == canonical_form(g, reflections)


# ---------------------------------------------------------------------------
# Numerical invariants
# ---------------------------------------------------------------------------

def mahler_measure(f: PolyLike, prec: int = 128) -> arb:
    """Mahler measure ``|lc| * prod max(1, |root|)`` as a certified ``arb`` ball."""
    f = poly(f)
    old = ctx.prec
    ctx.prec = prec
    try:
        m = arb(abs(int(f.leading_coefficient())))
        for root, mult in f.complex_roots():
            a = abs(root)
            if a > 1:
                m *= a ** mult
            elif not a < 1:
                # ball straddles 1: max(1, |r|) = max(1, upper) safe enclosure
                m *= arb(1).union(a) ** mult if a.upper() > 1 else 1
        return m
    finally:
        ctx.prec = old


def house(f: PolyLike) -> float:
    """Largest absolute value of a root."""
    return max(float(abs(r)) for r, _ in poly(f).complex_roots())


def bits(p: PolyLike) -> int:
    return int(poly(p).height_bits())


def log_height(f: PolyLike) -> float:
    return math.log(max(abs(int(c)) for c in poly(f).coeffs()))
