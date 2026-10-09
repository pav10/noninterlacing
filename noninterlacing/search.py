"""Searches for non-interlacing polynomials.

* :func:`totally_real_census` -- every monic, squarefree, totally real
  ``f in Z[x]`` of degree ``d`` and span ``<= S``, one per class under
  ``x -> x + n`` (and optionally ``x -> -x``).  Classical method of
  Robinson: the normalised derivatives ``p_k = f^{(k)}/k!`` are all totally
  real; ``p_k = (known part) + a_k`` and its critical points are the
  (already known) roots of ``p_{k+1}``, so real-rootedness of ``p_k`` and the
  location of its roots confine ``a_k`` to an explicit interval.
* :func:`levels_census` -- the census annotated with the interlacing level.
* :func:`random_search` -- random totally real polynomials from sampled root
  configurations (the experiments (i)-(ii) of Section 5.4 of the paper).
"""

from __future__ import annotations

import math
import random
import time
from dataclasses import dataclass
from math import comb
from typing import Iterator, Optional

import numpy as np
from flint import fmpz_poly

from .poly import canonical_form, coeffs_desc, is_irreducible, to_str
from .realroots import is_totally_real


def _real_roots_float(c_desc) -> Optional[np.ndarray]:
    """Sorted real roots (floats) or ``None`` if some root is clearly non-real."""
    c = np.asarray(c_desc, dtype=float)
    n = len(c) - 1
    if n <= 0:
        return np.array([])
    if n == 1:
        return np.array([-c[1] / c[0]])
    C = np.zeros((n, n))
    C[0, :] = -c[1:] / c[0]
    C[np.arange(1, n), np.arange(n - 1)] = 1.0
    r = np.linalg.eigvals(C)
    if np.max(np.abs(r.imag)) > 1e-6 * max(1.0, float(np.max(np.abs(r.real)))):
        return None
    return np.sort(r.real)


def totally_real_census(d: int, max_span: float, reflections: bool = True,
                        irreducible_only: bool = False,
                        max_count: Optional[int] = None) -> Iterator[fmpz_poly]:
    """Yield the census polynomials (canonical representatives)."""
    S = float(max_span)
    seen = set()
    eps = 1e-9
    count = 0
    binom = [[comb(i, k) for k in range(d + 1)] for i in range(d + 1)]
    for top in range(-(d - 1), 1):  # a_{d-1}: root mean in [0, 1)
        a = [0] * (d + 1)
        a[d] = 1
        a[d - 1] = top
        m = -top / d
        A0 = m - (d - 1) * S / d - eps
        B0 = m + (d - 1) * S / d + eps
        # node: (k, coefficients, roots of p_{k+1}, window [A, B] for the roots of f)
        stack = [(d - 2, a, np.array([m]), A0, B0)]
        while stack:
            k, a, xi, A, B = stack.pop()
            if k < 0:
                f = fmpz_poly(a)
                if xi[-1] - xi[0] > S + 1e-9 or not is_totally_real(f):
                    continue
                if irreducible_only and not is_irreducible(f):
                    continue
                key = canonical_form(f, reflections)
                kt = tuple(int(c) for c in key.coeffs())
                if kt in seen:
                    continue
                seen.add(kt)
                yield key
                count += 1
                if max_count is not None and count >= max_count:
                    return
                continue
            D = d - k  # degree of p_k = sum_{i>=k} C(i,k) a_i x^{i-k}
            pk_desc = [binom[i][k] * a[i] for i in range(d, k, -1)] + [0]
            Pc = np.array(pk_desc, dtype=float)
            vals = -np.polyval(Pc, xi)
            lo, hi = -math.inf, math.inf
            for i in range(1, D):
                v = vals[i - 1]
                if (D - i) % 2 == 0:
                    lo = max(lo, v)
                else:
                    hi = min(hi, v)
            lo = max(lo, -np.polyval(Pc, B))
            if D % 2 == 0:
                lo = max(lo, -np.polyval(Pc, A))
            else:
                hi = min(hi, -np.polyval(Pc, A))
            klo = math.ceil(lo - 1e-7 * (1 + abs(lo)))
            khi = math.floor(hi + 1e-7 * (1 + abs(hi)))
            for ak in range(klo, khi + 1):
                b = list(a)
                b[k] = ak
                Pc[-1] = ak
                r = _real_roots_float(Pc)
                if r is None or r[-1] - r[0] > S + 1e-7:
                    continue
                # roots of f lie in [r_max - S, r_min + S]
                stack.append((k - 1, b, r, max(A, r[-1] - S - eps), min(B, r[0] + S + eps)))


@dataclass
class CensusEntry:
    f: fmpz_poly
    span: float
    irreducible: bool
    level: int
    margin: float

    def as_dict(self) -> dict:
        return {"f": to_str(self.f), "span": round(self.span, 6),
                "irreducible": self.irreducible, "level": self.level,
                "margin": round(self.margin, 6)}


def levels_census(d: int, max_span: float, irreducible_only: bool = True,
                  reflections: bool = True) -> list:
    """Census with ``m(f)`` (and the level-1 margin) for each polynomial."""
    from .interlace import Interlacing

    out = []
    for f in totally_real_census(d, max_span, reflections, irreducible_only):
        I = Interlacing(f, check=False)
        lvl = I.min_level()
        mm = I.margin(1).value
        out.append(CensusEntry(f, I.roots.span, is_irreducible(f), lvl, mm))
    return out


def random_totally_real(d: int, rng: random.Random, low: float = -8, high: float = 8,
                        cluster: Optional[float] = None, tries: int = 1000) -> Optional[fmpz_poly]:
    """A random monic totally real squarefree polynomial (rounded root sample)."""
    for _ in range(tries):
        if cluster is None:
            roots = [rng.uniform(low, high) for _ in range(d)]
        else:
            # clustered: groups of nearby roots (large discriminant regime of the paper)
            centres = [rng.uniform(low, high) for _ in range(max(1, d // 3))]
            roots = [rng.choice(centres) + rng.uniform(-cluster, cluster) for _ in range(d)]
        c = np.poly(roots)
        coeffs = [int(round(x)) for x in c]
        f = fmpz_poly(list(reversed(coeffs)))
        if is_totally_real(f) and is_irreducible(f):
            return f
    return None


def random_search(d: int, samples: int, seed: int = 0, cluster: Optional[float] = None,
                  low: float = -8, high: float = 8, time_limit: Optional[float] = None,
                  stop_at_first: bool = False) -> dict:
    """Random totally real irreducible polynomials, tested for interlacing."""
    from .interlace import Interlacing

    rng = random.Random(seed)
    start = time.time()
    tested, noninterlacing, margins = 0, [], []
    while tested < samples:
        if time_limit is not None and time.time() - start > time_limit:
            break
        f = random_totally_real(d, rng, low, high, cluster)
        if f is None:
            continue
        I = Interlacing(f, check=False)
        g = I.find(1)
        tested += 1
        if g is None:
            noninterlacing.append(to_str(f))
            if stop_at_first:
                break
    return {"degree": d, "tested": tested, "noninterlacing": noninterlacing,
            "seconds": round(time.time() - start, 2)}
