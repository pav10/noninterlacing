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
