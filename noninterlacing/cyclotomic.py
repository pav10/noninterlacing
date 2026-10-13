"""Real cyclotomic bases, Yatsyna's admissible conductors, Chebyshev maps and
the set of degrees reachable by heredity.

``beta_n = 2 cos(2 pi / n)`` generates ``F_n = Q(zeta_n)^+`` and
``O_{F_n} = Z[beta_n]`` for every ``n`` (Washington, Prop. 2.16), so for real
cyclotomic fields ``m(Z[beta_n]) = m(O_{F_n})`` and the interlacing level of
``m_{beta_n}`` is an invariant of the field.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math
from math import gcd
from typing import Iterable, Optional

from flint import fmpz_poly

from .poly import chebyshev_C, is_irreducible, real_cyclotomic


def factorint(n: int) -> dict:
    out, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def is_prime(n: int) -> bool:
    return n >= 2 and factorint(n) == {n: 1}


def euler_phi(n: int) -> int:
    r = n
    for p in factorint(n):
        r = r // p * (p - 1)
    return r


def is_squarefree_int(n: int) -> bool:
    return all(e == 1 for e in factorint(n).values())


def real_cyclotomic_degree(n: int) -> int:
    return 1 if n <= 2 else euler_phi(n) // 2


def is_admissible(n: int) -> bool:
    """Conductors of Yatsyna's theorem as quoted in the paper (Thm 2.4):
    ``n`` squarefree, composite, not twice a prime, ``n > 20``, ``n != 30``."""
    if n <= 20 or n == 30 or not is_squarefree_int(n) or is_prime(n):
        return False
    if n % 2 == 0 and is_prime(n // 2):
        return False
    return True


def admissible_divisors(n: int) -> list:
    return [m for m in range(21, n + 1) if n % m == 0 and is_admissible(m)]


def has_admissible_divisor(n: int) -> bool:
    return bool(admissible_divisors(n))


def canonical_conductor(n: int) -> int:
    """``F_n = F_{n/2}`` for ``n = 2 mod 4``; return the conductor of ``F_n``."""
    return n // 2 if n % 4 == 2 else n


def chebyshev_composition(n: int, k: int) -> fmpz_poly:
    """``m_{beta_n}(C_k(x))``; equals ``m_{beta_{nk}}`` when irreducible."""
    return real_cyclotomic(n)(chebyshev_C(k))


def chebyshev_irreducible_predicted(n: int, k: int) -> bool:
    """``phi(nk) = k phi(n)``, i.e. every prime factor of ``k`` divides ``n``."""
    return all(n % p == 0 for p in factorint(k))


def odd_degree_conductors(limit: int) -> list:
    """All ``n <= limit`` (``n != 2 mod 4``) with ``phi(n)/2`` odd.

    These are ``n = 4`` and ``n = p^k`` with ``p = 3 mod 4`` (``k >= 1``); the
    prime powers with ``k >= 2`` (27, 49, 81, 121, ...) are the ones
    Proposition 5.7 of the paper does not mention.
    """
    out = []
    for n in range(3, limit + 1):
        if n % 4 == 2:
            continue
        if real_cyclotomic_degree(n) % 2 == 1:
            out.append(n)
    return out


def phi_table(limit: int) -> list:
    phi = list(range(limit + 1))
    for p in range(2, limit + 1):
        if phi[p] == p:
            for m in range(p, limit + 1, p):
                phi[m] -= phi[m] // p
    return phi


def _phi_inverse_bound(v: int) -> int:
    """An ``N`` with ``phi(n) > v`` for all ``n > N`` (Rosser--Schoenfeld:
    ``phi(n) > n / (e^gamma log log n + 3 / log log n)`` for ``n >= 3``)."""
    eg = 1.7810724179901979
    N = max(30, 2 * v)
    while True:
        ll = math.log(math.log(N))
        if N / (eg * ll + 3 / ll) > v:
            # the bound n / (...) is increasing for n >= N here, so N works
            return N
        N = int(N * 1.2) + 1


def reachable_degrees(bound: int) -> set:
    """The set ``D cap [1, bound]`` of Prop. 5.4:
    ``{ 2^j phi(n)/2 : n has an admissible divisor, j >= 0 }``."""
    nmax = _phi_inverse_bound(2 * bound)
    phi = phi_table(nmax)
    adm = [m for m in range(21, nmax + 1) if is_admissible(m)]
    has = bytearray(nmax + 1)
    for m in adm:
        for k in range(m, nmax + 1, m):
            has[k] = 1
    out = set()
    for n in range(21, nmax + 1):
        if has[n]:
            e = phi[n] // 2
            while e <= bound:
                out.add(e)
                e *= 2
    return out


@dataclass
class CyclotomicLevel:
    n: int
    degree: int
    admissible: bool
    admissible_divisors: list
    m: int
    margin: float
    method: str

    def as_dict(self) -> dict:
        return dict(self.__dict__)


def cyclotomic_level(n: int, method: str = "auto") -> CyclotomicLevel:
    """``m(O_{F_n})`` computed directly from ``m_{beta_n}``."""
    from .interlace import Interlacing

    f = real_cyclotomic(n)
    I = Interlacing(f)
    mm = I.margin(1)
    m = I.min_level(method)
    used = method if method != "auto" else ("dfs" if f.degree() <= 20 else "milp")
    return CyclotomicLevel(n, f.degree(), is_admissible(n), admissible_divisors(n), m, mm.value, used)


def cyclotomic_table(nmin: int = 3, nmax: int = 100, max_degree: int = 60,
                     min_degree: int = 6) -> list:
    """Levels of all real cyclotomic fields with conductor in range (one per field)."""
    out = []
    for n in range(nmin, nmax + 1):
        if n % 4 == 2:
            continue
        e = real_cyclotomic_degree(n)
        if e < min_degree or e > max_degree:
            continue
        out.append(cyclotomic_level(n))
    return out


# ---------------------------------------------------------------------------
# Upper bounds from composita: m(Psi_n) <= 2 m(Psi_{n1}) m(Psi_{n2})
# ---------------------------------------------------------------------------

def psi_interlacing_predicted(n: int) -> bool:
    """McKee--Smyth: ``Psi_n`` is interlacing iff rad(odd part of n) is 1, a prime or 15."""
    odd = [p for p in factorint(n) if p != 2]
    return len(odd) <= 1 or sorted(odd) == [3, 5]


def compositum_witness(n: int, n1: int, n2: int, g1=None, g2=None) -> fmpz_poly:
    """A level-``2 j1 j2`` interlacer of ``Psi_n`` built from interlacers of ``Psi_{n1}``, ``Psi_{n2}``.

    ``n = n1 n2`` with ``gcd(n1, n2) = 1`` and ``n1, n2 > 2``; ``g_i`` is a level-``j_i``
    interlacer of ``Psi_{n_i}`` (found if not given).  ``F = Q(zeta_{n1})^+ Q(zeta_{n2})^+``
    has index 2 in ``K = Q(zeta_n)^+`` and, the discriminants being coprime,
    ``O_F = O_1 (x) O_2`` and ``O_F^vee = O_1^vee O_2^vee``.  So ``mu = mu_1 mu_2`` with
    ``mu_i = g_i(beta_{n_i}) / Psi_{n_i}'(beta_{n_i})`` is totally positive in ``O_F^vee``
    of trace ``j1 j2``, hence in ``O_K^vee`` of trace ``2 j1 j2``.  The result
    ``g = mu Psi_n'(alpha) mod Psi_n`` is checked exactly (integrality and interlacing).
    """
    from flint import fmpq_poly

    from .interlace import Interlacing
    from .realroots import interlaces

    if n1 * n2 != n or gcd(n1, n2) != 1 or min(n1, n2) <= 2:
        raise ValueError("need n = n1 n2 with coprime n1, n2 > 2")
    P = fmpq_poly(real_cyclotomic(n))
    mu = fmpq_poly([1])
    for ni, other, gi in ((n1, n2, g1), (n2, n1, g2)):
        fi = real_cyclotomic(ni)
        if gi is None:
            gi = Interlacing(fi).find(1)
            if gi is None:
                raise ValueError(f"Psi_{ni} has no level-1 interlacer")
        # beta_{n_i} = 2cos(2 pi other / n) = C_other(alpha)
        b = fmpq_poly(chebyshev_C(other)) % P
        num = fmpq_poly(gi)(b) % P
        den = fmpq_poly(fi.derivative())(b) % P
        gg, inv, _ = den.xgcd(P)
        if gg != 1:
            raise ArithmeticError("Psi_ni'(beta) not invertible")
        mu = mu * num * inv % P
    g = mu * fmpq_poly(P.derivative()) % P
    if g.denom() != 1:
        raise ArithmeticError("compositum witness is not integral")
    g = fmpz_poly(g.numer())
    if not interlaces(real_cyclotomic(n), g):
        raise ArithmeticError("compositum witness does not interlace")
    return g


def compositum_splits(n: int) -> list:
    """Coprime splits ``n = n1 n2`` (``n1 < n2``, both > 2) with ``Psi_{n1}``, ``Psi_{n2}`` interlacing."""
    out = []
    for n1 in range(3, n):
        n2, r = divmod(n, n1)
        if r == 0 and n1 < n2 and gcd(n1, n2) == 1 and n2 > 2 \
                and psi_interlacing_predicted(n1) and psi_interlacing_predicted(n2):
            out.append((n1, n2))
    return out
