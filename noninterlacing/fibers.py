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

from .lattice import arb_to_scaled_int, lll_columns
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
        self._codiff = {}

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

    def embed(self, x: fmpq_poly, m: int) -> fmpq_poly:
        """An element ``x(alpha_m)`` of ``F_m = Q(zeta_m)^+`` (``m | n``) in the coordinates of ``K``.

        ``alpha_m = 2cos(2 pi/m) = C_{n/m}(alpha)``; substitution is exact.
        """
        if self.n % m:
            raise ValueError("m must divide n")
        C = self.cheb(self.n // m)
        out = fmpq_poly([0])
        for c in reversed(fmpq_poly(x).coeffs()):
            out = (out * C + c) % self.F
        return out

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

        key = tuple(H)
        if key not in self._codiff:
            if len(H) == 1:  # an HNF/LLL-reduced basis: alpha^i/Psi'(alpha) has huge coefficients
                basis = _zbasis(list(self.codifferent_K), self.d)
            else:
                imgs = [self.trace(H, b) for b in self.codifferent_K]
                basis = _zbasis(imgs, self.d, allow_deficient=True)
            if len(basis) * len(H) != self.d:
                raise ArithmeticError("codifferent rank mismatch")
            self._codiff[key] = basis
        return self._codiff[key]

    # ------------------------------------------------------- numerics
    def root(self, k: int, prec: int) -> arb:
        with precision(prec):
            return 2 * (2 * arb.pi() * k / self.n).cos()

    def values(self, x: fmpq_poly, ks: list, prec: int) -> list:
        """``sigma_k(x)`` for ``k`` in ``ks`` as ``arb`` balls of relative accuracy ``~2^-(prec-16)``.

        Elements of subfields written in the power basis of ``alpha`` can have huge
        coefficients and cancel heavily, so the working precision is raised until every
        ball is that accurate.
        """
        num = arb_poly([int(c) for c in x.numer().coeffs()])
        q = int(x.denom())
        bits = max([abs(int(c)).bit_length() for c in x.numer().coeffs()] + [1])
        work = max(prec, bits + 2 * self.d) + 32
        while True:
            with precision(work):
                vals = [num(self.root(k, work)) / q for k in ks]
            if all(_accurate(v, prec - 16) for v in vals):
                return vals
            work *= 2

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
        """Lattice data of ``Tr: O_{F''}^vee -> O_{F'}^vee`` in the bases ``W``, ``W1`` (cached)."""
        key = (tuple(H_small), tuple(H_big))
        if key not in self._fibre_cache:
            from .pari import pari

            P = pari()
            W = self.codifferent(H_small)
            W1 = self.codifferent(H_big)
            R = self.cosets(H_small, within=H_big)
            taus = self.cosets(H_big)
            e1, e2 = len(W1), len(W)
            if len(taus) != e1 or e1 * len(R) != e2:
                raise ArithmeticError("index mismatch")
            # coordinates with respect to W1: x = sum y_i W1_i  <=>  y = Cinv x[rows1]
            C = fmpq_mat([[W1[c][r] for c in range(e1)] for r in range(self.d)])
            rows1 = _pivot_rows(C)
            Cinv = fmpq_mat([[C[r, c] for c in range(e1)] for r in rows1]).inv()
            coords = lambda x: _coords(x, C, rows1, Cinv)
            Tm = []
            for w in W:
                y = coords(self.trace(R, w))
                if y is None:
                    raise ArithmeticError("Tr(O_F''^vee) is not inside O_F'^vee")
                Tm.append(y)
            Tp = P.matrix(e1, e2, [Tm[c][r] for r in range(e1) for c in range(e2)])
            HU = P.mathnf(Tp, 1)  # Tp * U = [0 | H]; Tr is onto O_F'^vee, so H = 1
            Hm, Um = HU[0], HU[1]
            if [[int(Hm[r, c]) for c in range(e1)] for r in range(e1)] != [[int(r == c) for c in range(e1)] for r in range(e1)]:
                raise ArithmeticError("the relative trace is not onto O_F'^vee")
            m = e2 - e1
            Ulast = fmpz_mat([[int(Um[r, m + c]) for c in range(e1)] for r in range(e2)])
            kz = fmpz_mat([[int(Um[r, c]) for r in range(e2)] for c in range(m)]).lll() if m else None
            kern = [[int(kz[i, r]) for r in range(e2)] for i in range(m)]
            Kz = fmpz_mat(kern) if m else None
            kred = (Kz, (fmpq_mat(Kz) * fmpq_mat(Kz.transpose())).inv()) if m else None
            solver = (C, rows1, Cinv, Ulast, kred)
            lam = [_comb(col, W) for col in kern]
            ks = [self.mul(k, r) for k in taus for r in R]
            lam_vals = [self.values(x, ks, 384) for x in lam]
            self._fibre_cache[key] = (W, R, taus, solver, lam, lam_vals)
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
              margin: float = 1e-6, scale: float = 1.0, limit: Optional[int] = None,
              method: str = "auto", time_limit: float = 600.0) -> FibreResult:
        """All ``mu in O_{K^{H_small}}^vee``, ``mu >> 0`` with ``Tr(mu) = nu`` in ``K^{H_big}``.

        ``scale < 1`` shrinks the search ball around the centre of the fibre (the points
        with all conjugates over each ``tau`` close to ``tau(nu)/t``): the result is then a
        subset, useful for finding witnesses, not for proving emptiness.  ``limit`` stops
        after that many points.  ``method="milp"`` looks for one point with HiGHS (for
        fibres with many points, where the circumscribed ball is far too large); it finds
        witnesses but proves nothing when it finds none.
        """
        from .pari import pari

        t0 = time.time()
        res = FibreResult()
        P = pari()
        W, R, taus, solver, lam, lam_vals = self._fibre_data(H_small, H_big)
        t, e1, m = len(R), len(taus), len(lam)
        res.rank = m
        nu = fmpq_poly(nu) % self.F
        nu_vals = self.values(nu, taus, prec)
        if not all(v > 0 for v in nu_vals):
            return res  # nu not totally positive: empty fibre
        C, rows1, Cinv, Ulast, kred = solver
        y = _coords(nu, C, rows1, Cinv)
        if y is None:
            return res  # nu is not in O_{F'}^vee = Tr(O_{F''}^vee): empty fibre
        z0 = Ulast * fmpz_mat(e1, 1, y)
        if m:  # reduce z0 modulo the kernel (exact Babai in coefficient space) to keep mu0 small
            Kz, KKinv = kred
            w = KKinv * (fmpq_mat(Kz) * fmpq_mat(z0))
            kk = fmpz_mat(m, 1, [_round_q(w[i, 0]) for i in range(m)])
            z0 = z0 - Kz.transpose() * kk
        mu0 = _comb([int(z0[i, 0]) for i in range(len(W))], W)
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
        bound = scale * e1 * (t - 1) / t
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
                # Babai: move mu0 next to the centre so that the float search sees a small target
                lam_red = [_comb([U[cc][kk] for cc in range(m)], lam) for kk in range(m)]
                for _ in range(4):
                    kb = _babai(newcols, c, prec)
                    if not any(kb):
                        break
                    c = [c[r] + sum((newcols[kk][r] * kb[kk] for kk in range(m) if kb[kk]), arb(0))
                         for r in range(len(rows))]
                    mu0 = mu0 + _comb(kb, lam_red)
            if not all(_accurate(x, 64, absolute=True) for col in newcols for x in col) or \
                    not all(_accurate(x, 64, absolute=True) for x in c):
                raise ArithmeticError("fibre data not accurate enough for the float search")
            B = np.array([[float(newcols[kk][r].mid()) for kk in range(m)] for r in range(len(rows))])
            cv = np.array([float(x.mid()) for x in c])
            if float(np.max(np.abs(cv))) > 1e6 or float(np.max(np.abs(B))) > 1e9:
                raise ArithmeticError(f"fibre lattice badly scaled for the float search: "
                                      f"target {float(np.max(np.abs(cv))):.3g}, basis {float(np.max(np.abs(B))):.3g}")
            if method == "milp":
                pts = _milp_point(B, cv, -1.0 / t + 1e-7, time_limit)
            elif method == "lp":
                pts = _lp_enumerate(B, cv, -1.0 / t - margin, res)
            elif method == "fp":
                pts = _fincke_pohst(B, cv, bound * (1 + margin) + margin, res)
            else:  # Fincke-Pohst while cheap, LP-pruned search if it runs past its node budget
                try:
                    pts = list(_fincke_pohst(B, cv, bound * (1 + margin) + margin, res, budget=50000))
                except _Budget:
                    pts = _lp_enumerate(B, cv, -1.0 / t - margin, res)
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
                if limit is not None and len(res.points) >= limit:
                    break
        res.seconds = time.time() - t0
        return res


def _accurate(v: arb, bits: int, absolute: bool = False) -> bool:
    """Radius of ``v`` below ``2^-bits`` (times ``max(1, |v|)`` unless ``absolute``)."""
    r = float(v.rad())
    if not math.isfinite(r) or not math.isfinite(float(v.mid())):
        return False
    scale = 1.0 if absolute else max(1.0, abs(float(v.mid())))
    return r <= scale * 2.0 ** -bits


def _coords(x: fmpq_poly, C: fmpq_mat, rows: list, Cinv: fmpq_mat):
    """Integer coordinates of ``x`` in the basis whose coordinate columns are ``C``, or ``None``."""
    e = C.ncols()
    y = Cinv * fmpq_mat(e, 1, [x[r] for r in rows])
    if any(y[i, 0].q != 1 for i in range(e)):
        return None
    yz = [int(y[i, 0].p) for i in range(e)]
    for r in range(C.nrows()):
        if sum((C[r, i] * yz[i] for i in range(e) if yz[i]), fmpq(0)) != x[r]:
            return None
    return yz


def _round_q(q: fmpq) -> int:
    p, d = int(q.p), int(q.q)
    return (2 * p + d) // (2 * d)


def _pivot_rows(H) -> list:
    """Indices of rows of a full-column-rank matrix forming an invertible square block."""
    rows = []
    for r in range(H.nrows()):
        cand = rows + [r]
        M = fmpq_mat([[H[i, c] for c in range(H.ncols())] for i in cand])
        if M.rank() == len(cand):
            rows = cand
            if len(rows) == H.ncols():
                return rows
    raise ArithmeticError("matrix does not have full column rank")


def _babai(cols: list, c: list, prec: int) -> list:
    """Integer ``k`` rounding the least-squares solution of ``sum_k cols[k] k ~ -c`` (``arb``)."""
    from flint import arb_mat

    m, N = len(cols), len(c)
    with precision(prec):
        B = arb_mat(N, m, [cols[k][r] for r in range(N) for k in range(m)])
        cv = arb_mat(N, 1, [c[r] for r in range(N)])
        G = B.transpose() * B
        w = G.solve(-(B.transpose() * cv))
        return [arb_to_scaled_int(w[k, 0], 0) for k in range(m)]  # exact rounding of the midpoint


def _comb(z, W) -> fmpq_poly:
    out = fmpq_poly([0])
    for c, w in zip(z, W):
        if c:
            out += int(c) * w
    return out


def _milp_point(B: np.ndarray, c: np.ndarray, lower: float, time_limit: float) -> list:
    """At most one integer ``z`` with ``B z + c >= lower`` (HiGHS feasibility; floats)."""
    from scipy.optimize import Bounds, LinearConstraint, milp

    from .simplex import _Silence

    m = B.shape[1]
    with _Silence():
        r = milp(np.zeros(m), constraints=[LinearConstraint(B, lower - c, np.inf)],
                 integrality=np.ones(m), bounds=Bounds(-np.inf, np.inf),
                 options={"disp": False, "time_limit": time_limit})
    return [[int(round(v)) for v in r.x]] if r.x is not None else []


def _lp_enumerate(B: np.ndarray, c: np.ndarray, lower: float, res: FibreResult):
    """Yield all integer ``z`` with ``B z + c >= lower``: depth-first, each range from two LPs.

    Every coordinate range is the exact projection of the (relaxed) polytope, so no node is
    wasted on parts of the ball outside it; ``lower`` already carries the safety margin.
    """
    from scipy.optimize import linprog

    m = B.shape[1]
    z = [0] * m

    def rng(i):
        # variables z_0..z_i free, z_{i+1}.. fixed: A x <= b with A = -B[:, :i+1]
        A = -B[:, :i + 1]
        b = c - lower + B[:, i + 1:] @ np.array(z[i + 1:], dtype=float)
        e = np.zeros(i + 1)
        e[i] = 1.0
        lo = linprog(e, A_ub=A, b_ub=b, bounds=[(None, None)] * (i + 1), method="highs")
        if lo.status != 0:
            return None
        hi = linprog(-e, A_ub=A, b_ub=b, bounds=[(None, None)] * (i + 1), method="highs")
        if hi.status != 0:
            return None
        return math.ceil(lo.fun - 1e-6), math.floor(-hi.fun + 1e-6)

    def rec(i):
        res.nodes += 1
        r = rng(i)
        if r is None:
            return
        for k in range(r[0], r[1] + 1):
            z[i] = k
            if i == 0:
                yield list(z)
            else:
                yield from rec(i - 1)
        z[i] = 0

    yield from rec(m - 1)


class _Budget(Exception):
    pass


def _fincke_pohst(B: np.ndarray, c: np.ndarray, r2: float, res: FibreResult, budget: Optional[int] = None):
    """Yield all integer ``z`` with ``|| B z + c ||^2 <= r2`` (Schnorr--Euchner).

    With ``budget``, raise :class:`_Budget` once that many nodes have been visited."""
    start = res.nodes
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
        if budget is not None and res.nodes - start > budget:
            raise _Budget
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


def tower_levels(K: CyclotomicField, chain: list, j: int, log=None, start: Optional[list] = None,
                 on_stage=None) -> list:
    """Level-``j`` elements up the tower ``K^{chain[0]} subset K^{chain[1]} subset ...``.

    ``chain`` is a decreasing list of subgroups starting with ``G`` (the field ``Q``);
    each stage is obtained from the previous one by fibres, so the count at every
    stage is complete.  Returns one :class:`Stage` per subgroup.
    """
    if sorted(chain[0]) != K.G:
        raise ValueError("the chain must start at G (the field Q)")
    stages = list(start) if start else [Stage(chain[0], 1, [fmpq_poly([j])], 1)]
    if [s.H for s in stages] != chain[:len(stages)]:
        raise ValueError("resumed stages do not match the chain")
    for H_big, H_small in list(zip(chain, chain[1:]))[len(stages) - 1:]:
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
        if on_stage:
            on_stage(stages)
        if log:
            log(f"degree {st.degree}: {count} elements of trace {j} in {len(reps)} orbits "
                f"({cand} candidates, {st.seconds} s)")
    return stages


def interlacer_from_element(K: CyclotomicField, mu: fmpq_poly) -> fmpz_poly:
    """``g = mu Psi_n'(alpha) mod Psi_n`` for ``mu in O_K^vee``, checked to be an interlacer of ``Psi_n``."""
    from .realroots import interlaces

    g = fmpq_poly(mu) * fmpq_poly(K.f.derivative()) % K.F
    if g.denom() != 1:
        raise ArithmeticError("mu is not in O_K^vee")
    gz = fmpz_poly([int(c) for c in g.numer().coeffs()])
    if int(gz.leading_coefficient()) <= 0 or not interlaces(K.f, gz):
        raise ArithmeticError("mu is not totally positive")
    return gz


def central_search(K: CyclotomicField, H: list, j: int, scale: float, time_limit: float = 300.0,
                   max_orbits: Optional[int] = None, log=None) -> Optional[fmpz_poly]:
    """Look for a level-``j`` interlacer of ``Psi_n`` in the fibres over central elements of ``K^H``.

    The fibre over ``nu`` has volume proportional to ``prod tau(nu)``, largest for ``nu``
    near the centre of the trace-``j`` simplex of ``K^H``.  This takes the elements of
    ``K^H`` within ``scale`` of the centre (:meth:`CyclotomicField.fibre` over ``Q``), one per
    Galois orbit, most central first, and searches each fibre with HiGHS.  A result is a
    checked interlacer; ``None`` proves nothing.
    """
    import math

    cen = K.fibre(H, K.G, fmpq_poly([j]), scale=scale, limit=5000)
    taus = K.cosets(H)
    reps = {}
    for nu in cen.points:
        reps.setdefault(K.orbit_key(nu, H), nu)
    order = sorted(reps.values(), key=lambda nu: -sum(math.log(float(v.mid())) for v in K.values(nu, taus, 128)))
    if log:
        log(f"{len(cen.points)} central elements of trace {j} in {len(order)} orbits")
    for i, nu in enumerate(order[:max_orbits]):
        fr = K.fibre([1], H, nu, limit=1, method="milp", time_limit=time_limit)
        if log:
            log(f"orbit {i}: {len(fr.points)} point(s), {fr.seconds:.1f} s")
        if fr.points:
            return interlacer_from_element(K, fr.points[0])
    return None

