"""Fibres of the relative trace between subfields of a real cyclotomic field.

Inside ``K = Q(zeta_n)^+ = Q(alpha)`` (``alpha = 2cos(2 pi/n)``, ``O_K = Z[alpha]``)
every subfield is ``K^H`` for a subgroup ``H`` of ``G = (Z/n)^* / {+-1}``, which
acts by ``sigma_a(alpha) = C_a(alpha)``.  Elements are ``fmpq_poly`` in ``alpha``
reduced mod ``Psi_n``; ``sigma_k`` (``k`` in ``G``) is the real embedding
``alpha -> 2cos(2 pi k/n)``.

For subfields ``F'' = K^{H''} subset F' = K^{H'}`` (``H'' < H'``) the relative trace
maps the totally positive elements of ``O_{F''}^vee`` of trace ``j`` onto totally
positive elements of ``O_{F'}^vee`` of trace ``j`` (``Tr_{F'}(Tr(mu) x) = Tr_{F''}(mu x)``
for ``x in O_{F'}``).  Hence the level-``j`` elements of ``F''`` are the union, over
the level-``j`` elements ``nu`` of ``F'``, of the fibres

    { mu in O_{F''}^vee : Tr_{F''/F'}(mu) = nu,  mu >> 0 },

and Galois-conjugate ``nu`` have conjugate fibres, so one ``nu`` per orbit is
enough.  A fibre is the set of points of a coset of the lattice
``ker Tr_{F''/F'}`` lying in a product of ``[F':Q]`` open simplices of dimension
``t - 1`` (``t = [F'':F']``).  :meth:`CyclotomicField.fibre` enumerates it with a
Fincke--Pohst search over the circumscribed ball (in coordinates scaled by
``1/tau(nu)``, an LLL-reduced basis, floating point with a safety margin) and
checks every candidate exactly with Sturm sequences.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from functools import cached_property
from math import gcd
from typing import Iterable, Optional

import numpy as np
from flint import arb, arb_poly, ctx, fmpq, fmpq_mat, fmpq_poly, fmpz_mat, fmpz_poly

from .lattice import lll_columns
from .poly import chebyshev_C, real_cyclotomic
from .realroots import is_totally_positive_exact


class precision:
    def __init__(self, bits: int):
        self.bits = bits

    def __enter__(self):
        self.old = ctx.prec
        ctx.prec = self.bits

    def __exit__(self, *exc):
        ctx.prec = self.old


@dataclass
class FibreResult:
    points: list = field(default_factory=list)  # fmpq_poly in alpha
    candidates: int = 0  # lattice points inside the ball (checked exactly)
    nodes: int = 0
    rank: int = 0
    seconds: float = 0.0


class CyclotomicField:
    """``K = Q(zeta_n)^+`` with its Galois action and subfield codifferents."""

    def __init__(self, n: int):
        if n < 3 or n % 4 == 2:
            raise ValueError("use a conductor n >= 3, n not 2 mod 4")
        self.n = n
        self.f = real_cyclotomic(n)
        self.d = self.f.degree()
        self.F = fmpq_poly(self.f)
        self.G = [a for a in range(1, n // 2 + 1) if gcd(a, n) == 1]
        if len(self.G) != self.d:
            raise ArithmeticError("Galois group size mismatch")
        self._cheb = {}
        self._mats = {}
        self._fibre_cache = {}

    # ------------------------------------------------------------ the group
    def rep(self, a: int) -> int:
        a %= self.n
        return min(a, self.n - a)

    def mul(self, a: int, b: int) -> int:
        return self.rep(a * b)

    def subgroup(self, gens: Iterable[int]) -> list:
        H = {1}
        frontier = [1]
        gens = [self.rep(g) for g in gens]
        while frontier:
            x = frontier.pop()
            for g in gens:
                y = self.mul(x, g)
                if y not in H:
                    H.add(y)
                    frontier.append(y)
        return sorted(H)

    def fixer(self, m: int) -> list:
        """The subgroup fixing ``F_m = Q(zeta_m)^+`` (``m | n``): ``a = +-1 mod m``."""
        if self.n % m:
            raise ValueError("m must divide n")
        return [a for a in self.G if a % m in (1, m - 1) or m <= 2]

    def cosets(self, H: list, within: Optional[list] = None) -> list:
        """Representatives of ``within / H`` (``within`` defaults to ``G``)."""
        Hs = set(H)
        seen, reps = set(), []
        for a in (within if within is not None else self.G):
            if a in seen:
                continue
            reps.append(a)
            seen.update(self.mul(a, h) for h in Hs)
        return reps

    # ---------------------------------------------------- exact arithmetic
    def cheb(self, a: int) -> fmpq_poly:
        a = self.rep(a)
        if a not in self._cheb:
            self._cheb[a] = fmpq_poly(chebyshev_C(a)) % self.F
        return self._cheb[a]

    def act_matrix(self, a: int) -> fmpz_mat:
        """Matrix of ``sigma_a`` on the power basis: column ``i`` = ``C_a(alpha)^i mod Psi_n``."""
        a = self.rep(a)
        if a not in self._mats:
            C = fmpz_poly([int(c) for c in self.cheb(a).numer().coeffs()])
            col, cols = fmpz_poly([1]), []
            for _ in range(self.d):
                cols.append([int(col[r]) for r in range(self.d)])
                col = (col * C) % self.f
            self._mats[a] = fmpz_mat([[cols[c][r] for c in range(self.d)] for r in range(self.d)])
        return self._mats[a]

    def act(self, a: int, x: fmpq_poly) -> fmpq_poly:
        """``sigma_a(x)``: substitute ``alpha -> C_a(alpha)`` (mod ``Psi_n``)."""
        x = fmpq_poly(x) % self.F
        v = fmpq_mat(self.d, 1, [x[i] for i in range(self.d)])
        w = fmpq_mat(self.act_matrix(a)) * v
        return fmpq_poly([w[i, 0] for i in range(self.d)])

    def trace(self, H: list, x) -> fmpq_poly:
        out = fmpq_poly([0])
        for h in H:
            out += self.act(h, x)
        return out

    @cached_property
    def codifferent_K(self) -> list:
        """``alpha^i / Psi_n'(alpha)``: a basis of ``O_K^vee``."""
        _, inv, _ = fmpq_poly(self.f.derivative()).xgcd(self.F)
        return [(fmpq_poly([0] * i + [1]) * inv) % self.F for i in range(self.d)]

    def codifferent(self, H: list) -> list:
        """A Z-basis of ``O_{K^H}^vee = Tr_H(O_K^vee)`` (rank ``d/|H|``)."""
        from .orders import _zbasis

        if len(H) == 1:
            return list(self.codifferent_K)
        imgs = [self.trace(H, b) for b in self.codifferent_K]
        basis = _zbasis(imgs, self.d, allow_deficient=True)
        if len(basis) * len(H) != self.d:
            raise ArithmeticError("codifferent rank mismatch")
        return basis

    # ------------------------------------------------------- numerics
    def root(self, k: int, prec: int) -> arb:
        with precision(prec):
            return 2 * (2 * arb.pi() * k / self.n).cos()

    def values(self, x: fmpq_poly, ks: list, prec: int) -> list:
        num = arb_poly([int(c) for c in x.numer().coeffs()])
        q = int(x.denom())
        with precision(prec):
            return [num(self.root(k, prec)) / q for k in ks]

    def is_totally_positive(self, x: fmpq_poly, H: Optional[list] = None, prec: int = 256) -> bool:
        """Exact: certified ``arb`` signs at the embeddings of ``K^H``, Sturm only if a ball meets 0."""
        ks = self.cosets(H) if H else self.G
        vals = self.values(x, ks, prec)
        if any(v < 0 for v in vals):
            return False
        if all(v > 0 for v in vals):
            return True
        return is_totally_positive_exact(fmpz_poly([int(c) for c in x.numer().coeffs()]), self.f)

    def _fibre_data(self, H_small: list, H_big: list):
        key = (tuple(H_small), tuple(H_big))
        if key not in self._fibre_cache:
            from .pari import pari

            P = pari()
            W = self.codifferent(H_small)
            R = self.cosets(H_small, within=H_big)
            taus = self.cosets(H_big)
            if len(taus) * len(R) != len(W):
                raise ArithmeticError("index mismatch")
            TW = [self.trace(R, w) for w in W]
            den = 1
            for x in TW:
                den = den * int(x.denom()) // gcd(den, int(x.denom()))
            M = [[int(TW[c][r] * den) for c in range(len(W))] for r in range(self.d)]
            Mp = P.matrix(self.d, len(W), [x for row in M for x in row])
            Kb = P.matkerint(Mp)
            m = int(P.matsize(Kb)[1])
            if m != len(W) - len(taus):
                raise ArithmeticError("kernel rank mismatch")
            kern = [[int(Kb[r, c]) for r in range(len(W))] for c in range(m)]
            lam = [_comb(col, W) for col in kern]
            ks = [self.mul(k, r) for k in taus for r in R]
            lam_vals = [self.values(x, ks, 384) for x in lam]
            self._fibre_cache[key] = (W, R, taus, den, Mp, lam, lam_vals)
        return self._fibre_cache[key]

    def orbit_key(self, x: fmpq_poly, H: list) -> tuple:
        """A key constant on ``G``-orbits of elements of ``K^H`` (rounded embedding values)."""
        ks = self.cosets(H)
        v = {k: float(arb(z).mid()) for k, z in zip(ks, self.values(x, ks, 128))}
        idx = {}
        for k in ks:  # map each element of G to its coset representative
            for h in H:
                idx[self.mul(k, h)] = k
        best = None
        for a in ks:
            t = tuple(round(v[idx[self.mul(k, a)]], 8) for k in ks)
            if best is None or t < best:
                best = t
        return best

    def orbit_size(self, x: fmpq_poly, H: list) -> int:
        ks = self.cosets(H)
        vals = [float(arb(z).mid()) for z in self.values(x, ks, 128)]
        pos = {k: i for i, k in enumerate(ks)}
        idx = {}
        for k in ks:
            for h in H:
                idx[self.mul(k, h)] = k
        images = set()
        for a in ks:
            images.add(tuple(round(vals[pos[idx[self.mul(k, a)]]], 8) for k in ks))
        return len(images)

    # --------------------------------------------------------- fibres
    def fibre(self, H_small: list, H_big: list, nu: fmpq_poly, prec: int = 384,  # = cached precision
              margin: float = 1e-6) -> FibreResult:
        """All ``mu in O_{K^{H_small}}^vee``, ``mu >> 0`` with ``Tr(mu) = nu`` in ``K^{H_big}``."""
        from .pari import pari

        t0 = time.time()
        res = FibreResult()
        P = pari()
        W, R, taus, den, Mp, lam, lam_vals = self._fibre_data(H_small, H_big)
        t, e1, m = len(R), len(taus), len(lam)
        res.rank = m
        nu = fmpq_poly(nu) % self.F
        nu_vals = self.values(nu, taus, prec)
        if not all(v > 0 for v in nu_vals):
            return res  # nu not totally positive: empty fibre
        if (nu * den).denom() != 1:
            return res  # nu is not in Tr(O_{F''}^vee): empty fibre
        rhs = [int(nu[r] * den) for r in range(self.d)]
        sol = P.matsolvemod(Mp, 0, P.Col(rhs))
        if str(P.type(sol)) == "t_INT":
            return res  # nu is not in Tr(O_{F''}^vee): empty fibre
        mu0 = _comb([int(x) for x in sol], W)
        # weighted coordinates u_(tau, r) = (sigma_{tau r}(mu) - tau(nu)/t) / tau(nu)
        rows = [(k, r, i) for i, k in enumerate(taus) for r in R]
        ks = [self.mul(k, r) for (k, r, i) in rows]
        with precision(prec):
            inv_nu = [1 / v for v in nu_vals]
            cols = []
            for vals in lam_vals:
                cols.append([v * inv_nu[i] for v, (k, r, i) in zip(vals, rows)])
            v0 = self.values(mu0, ks, prec)
            c = [(v - nu_vals[i] / t) * inv_nu[i] for v, (k, r, i) in zip(v0, rows)]
        bound = e1 * (t - 1) / t
        if m == 0:
            pts = [[]]
            U = []
        else:
            U = lll_columns(cols, shift=prec // 2)
            with precision(prec):
                newcols = []
                for kk in range(m):
                    newcols.append([sum((cols[cc][r] * U[cc][kk] for cc in range(m) if U[cc][kk]),
                                        arb(0)) for r in range(len(rows))])
            B = np.array([[float(newcols[kk][r].mid()) for kk in range(m)] for r in range(len(rows))])
            cv = np.array([float(x.mid()) for x in c])
            pts = _fincke_pohst(B, cv, bound * (1 + margin) + margin, res)
        for zp in pts:
            if m:
                # float pre-test of positivity (u > -1/t); only clear violations are dropped
                u = B @ np.array(zp, dtype=float) + cv
                if float(u.min()) < -1.0 / t - margin:
                    continue
            z = [sum(U[cc][kk] * zp[kk] for kk in range(m)) for cc in range(m)] if m else []
            mu = mu0
            for cz, x in zip(z, lam):
                if cz:
                    mu = mu + cz * x
            res.candidates += 1
            if self.is_totally_positive(mu, H_small):
                res.points.append(mu)
        res.seconds = time.time() - t0
        return res


def _comb(z, W) -> fmpq_poly:
    out = fmpq_poly([0])
    for c, w in zip(z, W):
        if c:
            out += int(c) * w
    return out


def _fincke_pohst(B: np.ndarray, c: np.ndarray, r2: float, res: FibreResult):
    """Yield all integer ``z`` with ``|| B z + c ||^2 <= r2`` (Schnorr--Euchner)."""
    Q, Rm = np.linalg.qr(B)  # B = Q Rm, Rm upper triangular m x m
    m = Rm.shape[0]
    qc = Q.T @ c
    rest = float(c @ c - qc @ qc)
    r2 = r2 - max(rest, 0.0)
    if r2 < 0:
        return
    # || Rm z + qc ||^2 : coordinate i depends on z_i..z_{m-1}
    z = [0] * m

    def rec(i: int, partial: float):
        res.nodes += 1
        s = qc[i] + sum(Rm[i, j] * z[j] for j in range(i + 1, m))
        rii = Rm[i, i]
        centre = -s / rii
        room = r2 - partial
        if room < 0:
            return
        half = math.sqrt(room) / abs(rii)
        lo, hi = math.ceil(centre - half), math.floor(centre + half)
        for k in range(lo, hi + 1):
            v = rii * k + s
            p = partial + v * v
            if p > r2:
                continue
            z[i] = k
            if i == 0:
                yield list(z)
            else:
                yield from rec(i - 1, p)
        z[i] = 0

    yield from rec(m - 1, 0.0)


@dataclass
class Stage:
    """Level-``j`` elements of one subfield ``K^H`` of a tower, up to Galois."""

    H: list
    degree: int
    reps: list  # one element per G-orbit
    count: int  # number of elements (sum of orbit sizes)
    candidates: int = 0
    seconds: float = 0.0


def tower_levels(K: CyclotomicField, chain: list, j: int, log=None) -> list:
    """Level-``j`` elements up the tower ``K^{chain[0]} subset K^{chain[1]} subset ...``.

    ``chain`` is a decreasing list of subgroups starting with ``G`` (the field ``Q``);
    each stage is obtained from the previous one by fibres, so the count at every
    stage is complete.  Returns one :class:`Stage` per subgroup.
    """
    if sorted(chain[0]) != K.G:
        raise ValueError("the chain must start at G (the field Q)")
    stages = [Stage(chain[0], 1, [fmpq_poly([j])], 1)]
    for H_big, H_small in zip(chain, chain[1:]):
        t0 = time.time()
        found, cand = {}, 0
        for nu in stages[-1].reps:
            fr = K.fibre(H_small, H_big, nu)
            cand += fr.candidates
            for mu in fr.points:
                found.setdefault(K.orbit_key(mu, H_small), mu)
        reps = list(found.values())
        count = sum(K.orbit_size(mu, H_small) for mu in reps)
        st = Stage(H_small, K.d // len(H_small), reps, count, cand, round(time.time() - t0, 2))
        stages.append(st)
        if log:
            log(f"degree {st.degree}: {count} elements of trace {j} in {len(reps)} orbits "
                f"({cand} candidates, {st.seconds} s)")
    return stages

