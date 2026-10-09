"""Numerical semigroups (used for the set of interlacing levels)."""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from typing import Callable, Iterable, Optional


@dataclass(frozen=True)
class NumericalSemigroup:
    """A numerical semigroup given by its Apery set w.r.t. its multiplicity.

    ``apery[r]`` is the least element congruent to ``r`` modulo
    ``multiplicity``.  Then ``S = { apery[r] + k * multiplicity : k >= 0 }``.
    """

    multiplicity: int
    apery: tuple

    @classmethod
    def from_generators(cls, gens: Iterable[int]) -> "NumericalSemigroup":
        gens = sorted({int(g) for g in gens if g > 0})
        if not gens:
            raise ValueError("need at least one positive generator")
        g = 0
        for x in gens:
            g = gcd(g, x)
        if g != 1:
            raise ValueError("generators must be coprime")
        m = gens[0]
        # shortest paths in Z/m with edge weights gens (Dijkstra-free: Bellman style)
        INF = float("inf")
        best = [INF] * m
        best[0] = 0
        changed = True
        while changed:
            changed = False
            for r in range(m):
                if best[r] == INF:
                    continue
                for x in gens[1:]:
                    v = best[r] + x
                    rr = v % m
                    if v < best[rr]:
                        best[rr] = v
                        changed = True
        return cls(m, tuple(int(b) for b in best))

    @classmethod
    def from_oracle(
        cls,
        member: Callable[[int], bool],
        multiplicity: int,
        upper: Optional[int] = None,
    ) -> "NumericalSemigroup":
        """Build from a membership oracle, knowing the multiplicity ``m``.

        For each residue ``r`` mod ``m`` the oracle is queried on
        ``r, r + m, r + 2m, ...`` until a member is found (closure under
        addition guarantees everything above it in that class is a member).
        ``upper`` caps the search (raise if exceeded).
        """
        m = multiplicity
        apery = [0] * m
        for r in range(1, m):
            x = r
            while not member(x):
                x += m
                if upper is not None and x > upper:
                    raise RuntimeError(f"no element found in class {r} mod {m} below {upper}")
            apery[r] = x
        return cls(m, tuple(apery))

    def __contains__(self, x: int) -> bool:
        if x < 0:
            return False
        return x >= self.apery[x % self.multiplicity]

    @property
    def frobenius(self) -> int:
        """Largest integer not in ``S`` (``-1`` if ``S = N``)."""
        return max(self.apery) - self.multiplicity

    @property
    def gaps(self) -> list:
        return [x for x in range(1, self.frobenius + 1) if x not in self]

    @property
    def genus(self) -> int:
        return len(self.gaps)

    @property
    def conductor(self) -> int:
        return self.frobenius + 1

    @property
    def generators(self) -> list:
        """Minimal generating set."""
        m = self.multiplicity
        gens = [m]
        for w in sorted(self.apery):
            if w == 0:
                continue
            # w is a minimal generator iff it is not w' + s with w' in Apery, s in S\{0}
            if not any(w - w2 in self and w2 != 0 and w - w2 > 0 for w2 in self.apery):
                gens.append(w)
        return sorted(gens)

    def elements(self, upto: int) -> list:
        return [x for x in range(0, upto + 1) if x in self]

    def __str__(self) -> str:
        gens = ", ".join(map(str, self.generators))
        return f"<{gens}>"
