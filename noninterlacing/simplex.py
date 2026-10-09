"""Enumeration of lattice points in the trace-``j`` simplex.

Setting
-------
Let ``L`` be a full-rank lattice in ``R^d`` given by a basis ``w_0..w_{d-1}``
(columns of an embedding matrix ``E``) whose coordinate sums ``tr(w_k)`` are
integers.  For ``j >= 1`` put

    Delta_j = { x in R^d : x_i > 0 for all i,  x_1 + ... + x_d = j }.

For ``L`` the codifferent of an order ``O`` in a totally real (etale) algebra,
embedded by its ``d`` real embeddings, ``L cap Delta_j`` is exactly the set
of totally positive codifferent elements of trace ``j``; for ``O = Z[alpha]``
these are the interlacers of level ``j`` (see :mod:`noninterlacing.interlace`).

Algorithm
---------
1. ``L_0 = L cap {sum = 0}`` is computed exactly (unimodular completion of the
   integer trace vector), giving the affine lattice
   ``L cap {sum = j} = t_j + L_0``.
2. The basis of ``L_0`` is LLL-reduced for the Euclidean (= trace) form, in
   which ``Delta_j`` is a *regular* simplex, and orthogonalised:
   ``B = Q R`` with ``R`` upper triangular.  In the coordinates
   ``y = Q^T (x - j/d)`` the point with integer coordinates ``k`` is
   ``y = s_j + R k``, so ``y_i`` depends only on ``k_i, ..., k_{n-1}``.
3. Depth-first search fixes ``k_{n-1}, k_{n-2}, ..., k_0``.  At depth ``i``
   the admissible values of ``y_i`` given ``y_{i+1..}`` are read off the
   *facets of the projection* of ``Delta_j`` onto the coordinates
   ``y_i..y_{n-1}`` (precomputed once per lattice with Qhull; they scale
   with ``j``).  This is the exact slice, so every node of the search tree is
   a projection of a genuine point of the polytope: the tree is as small as
   it can be for the chosen basis order.  For very high dimension, levels
   whose projection has too many facets fall back to a cheaper valid
   relaxation (ball + individual half-spaces).
4. Leaves are certified: points farther than a margin from the boundary are
   accepted on floating-point evidence (errors are ~1e-15 relative), the
   remaining few are decided by an exact callback.

All search bounds are relaxed by a tolerance, so floating point error can
only add candidate nodes, never lose genuine points.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np
from flint import arb, arb_mat

from .lattice import lll_columns, unimodular_completion
from .realroots import precision

try:  # scipy is optional; without it every level uses the relaxation
    from scipy.spatial import ConvexHull
    from scipy.spatial import QhullError  # type: ignore[attr-defined]
except Exception:  # pragma: no cover
    ConvexHull = None

    class QhullError(Exception):
        pass


class _Abort(Exception):
    pass


@dataclass
class EnumerationResult:
    """Outcome of an enumeration at one trace level."""

    level: int
    count: int
    points: list = field(default_factory=list)
    nodes: int = 0
    complete: bool = True
    exact_checks: int = 0
    rejected: int = 0
    seconds: float = 0.0

    @property
    def exists(self) -> bool:
        return self.count > 0


def _cyclic_facet_bound(nverts: int, dim: int) -> float:
    """Upper bound theorem: facets of a ``dim``-polytope with ``nverts`` vertices."""
    if dim <= 1:
        return 2
    if nverts <= dim + 1:
        return nverts
    if dim % 2 == 0:
        m = dim // 2
        return nverts / (nverts - m) * math.comb(nverts - m, m)
    m = (dim - 1) // 2
    return 2 * math.comb(nverts - m - 1, m)


class SimplexLattice:
    """Lattice points of ``L`` in the open simplices ``Delta_j``.

    Parameters
    ----------
    embedding:
        Callable ``prec -> E`` returning the ``d x d`` embedding matrix (list
        of rows of ``arb``; ``E[i][k]`` = coordinate ``i`` of basis vector
        ``k``) with at least ``prec`` bits of relative accuracy.
    traces:
        Exact integer coordinate sums of the basis vectors.
    verify:
        Callable taking an integer coordinate vector (w.r.t. the input
        basis) and returning whether the point is *strictly* inside the
        positive orthant; it must be exact.  Used for points near the
        boundary (and for every point when ``rigorous=True``).
    """

    def __init__(
        self,
        embedding: Callable[[int], Sequence[Sequence[arb]]],
        traces: Sequence[int],
        verify: Optional[Callable[[list], bool]] = None,
        prec: int = 256,
        facet_cap: int = 60000,
        hull_bound_cap: float = 100000,
    ):
        self._embedding = embedding
        self.traces = [int(t) for t in traces]
        self.d = len(self.traces)
        self.n = self.d - 1
        self.verify = verify
        self.facet_cap = facet_cap
        self.hull_bound_cap = hull_bound_cap
        self.gcd, W = unimodular_completion(self.traces)
        if self.gcd == 0:
            raise ValueError("all traces vanish: not a lattice with a trace functional")
        self.p1 = W[0]
        self.Kcols = W[1:]
        self._setup(prec)

    # ------------------------------------------------------------------ setup
    def _try_setup(self, prec: int) -> bool:
        """Build the reduced basis at working precision ``prec``; report success."""
        d, n = self.d, self.n
        with precision(prec):
            E = arb_mat([list(row) for row in self._embedding(prec)])
            self.e1 = [sum((E[r, k] * self.p1[k] for k in range(d)), arb(0)) for r in range(d)]
            if n == 0:
                self.E = E
                return True
            K = arb_mat([[self.Kcols[c][r] for c in range(n)] for r in range(d)])
            B0 = E * K
            cols = [[B0[r, c] for r in range(d)] for c in range(n)]
            norms = [math.sqrt(sum(float(x) ** 2 for x in col)) for col in cols]
            shift = 64 - math.floor(math.log2(max(min(norms), 1e-300)))
            U = lll_columns(cols, shift)
            M = [[sum(self.Kcols[c][r] * U[c][k] for c in range(n)) for k in range(n)] for r in range(d)]
            B = E * arb_mat(M)
        for c in range(n):
            cn = math.sqrt(sum(float(B[r, c]) ** 2 for r in range(d)))
            if any(float(B[r, c].rad()) > cn * 2.0 ** -80 for r in range(d)):
                return False
        self.E, self.M, self.B_arb = E, M, B
        return True

    def _setup(self, prec: int) -> None:
        d, n = self.d, self.n
        for _ in range(8):
            if self._try_setup(prec):
                break
            prec *= 2
        else:  # pragma: no cover
            raise RuntimeError("could not reach working precision")
        self.prec = prec
        if n == 0:
            return
        Bf = np.array([[float(self.B_arb[r, c]) for c in range(n)] for r in range(d)])
        Q, R = np.linalg.qr(Bf)
        sgn = np.sign(np.diag(R))
        sgn[sgn == 0] = 1
        self.Q = Q * sgn
        self.R = R * sgn[:, None]
        self.Bf = Bf
        self.V = self.Q - self.Q.sum(axis=0) / d  # vertices of Delta_1 in y-coordinates
        self.rho1 = float(np.max(np.linalg.norm(self.V, axis=1)))
        self._levels = [None] * n  # facet data, computed lazily per depth
        self._level_cache: dict[int, tuple] = {}
        self._tf_cache: dict[int, np.ndarray] = {}

    def _facets(self, i: int):
        """Facet data for the slice at coordinate ``i`` (see module docstring)."""
        n, d = self.n, self.d
        if i == 0:
            A = -self.Q  # x_r = j/d + Q[r] . y >= 0  <=>  -Q[r] . y - j/d <= 0
            b = np.full(d, -1.0 / d)
            return ("facets", A[:, 0].copy(), A[:, 1:].copy(), b)
        dim = n - i
        pts = self.V[:, i:]
        if dim == 1:
            lo, hi = float(pts.min()), float(pts.max())
            A = np.array([[1.0], [-1.0]])
            b = np.array([-hi, lo])
            return ("facets", A[:, 0].copy(), A[:, 1:].copy(), b)
        if ConvexHull is not None and _cyclic_facet_bound(d, dim) <= self.hull_bound_cap:
            for opts in (None, "QJ"):
                try:
                    hull = ConvexHull(pts, qhull_options=opts) if opts else ConvexHull(pts)
                except (QhullError, ValueError):
                    continue
                eq = np.unique(np.round(hull.equations, 12), axis=0)
                if len(eq) > self.facet_cap:
                    break
                A, b = eq[:, :-1], eq[:, -1]
                if np.max(pts @ A.T + b) > 1e-8:
                    continue
                return ("facets", A[:, 0].copy(), A[:, 1:].copy(), b)
        Qi = self.Q[:, i]
        Qk = self.Q[:, i + 1:]
        w = np.linalg.norm(self.Q[:, :i], axis=1)
        return ("ball", Qi, Qk, w)

    def _level_data(self, i: int):
        L = self._levels[i]
        if L is None:
            L = self._levels[i] = self._facets(i)
        return L

    @property
    def facet_counts(self) -> list:
        """Facet counts per depth (``None``: relaxation used; forces all levels)."""
        return [len(L[1]) if L[0] == "facets" else None
                for L in (self._level_data(i) for i in range(self.n))]

    # ------------------------------------------------------------ level data
    def _level(self, j: int):
        if j in self._level_cache:
            return self._level_cache[j]
        d, n = self.d, self.n
        q = j // self.gcd
        a0 = [p * q for p in self.p1]
        c0 = j / d
        with precision(self.prec):
            t = [x * q for x in self.e1]
            s, a0 = self._centre(t, a0, c0)
        self._level_cache[j] = (s, a0)
        return s, a0

    def _centre(self, t, a0, c0):
        d, n = self.d, self.n
        for _ in range(60):
            v = np.array([c0 - float(x) for x in t])
            z = np.linalg.solve(self.R, self.Q.T @ v)
            kk = np.round(z).astype(object)
            if not any(kk):
                break
            kk = [int(x) for x in kk]
            for r in range(d):
                t[r] = t[r] + sum((self.B_arb[r, c] * kk[c] for c in range(n) if kk[c]), arb(0))
                a0[r] += sum(self.M[r][c] * kk[c] for c in range(n))
        tf = np.array([float(x) for x in t])
        self._tf_cache[round(c0 * d)] = tf
        return self.Q.T @ (tf - c0), a0

    def point_coords(self, a0: list, k: Sequence[int]) -> list:
        return [a0[r] + sum(self.M[r][c] * k[c] for c in range(self.n)) for r in range(self.d)]

    # ------------------------------------------------------------- the search
    def _range(self, i: int, y: np.ndarray, j: int, tol: float):
        kind, a0, ak, b = self._level_data(i)
        if kind == "facets":
            rhs = tol - j * b
            if ak.shape[1]:
                rhs = rhs - ak @ y[i + 1:]
            lo, hi = -math.inf, math.inf
            pos = a0 > 1e-12
            neg = a0 < -1e-12
            if pos.any():
                hi = float(np.min(rhs[pos] / a0[pos]))
            if neg.any():
                lo = float(np.max(rhs[neg] / a0[neg]))
            zero = ~(pos | neg)
            if zero.any() and np.min(rhs[zero]) < -tol:
                return 1.0, -1.0
            return lo, hi
        # relaxation: ball of radius j*rho1 and each half-space separately
        Qi, Qk, w = a0, ak, b
        yk = y[i + 1:]
        R2 = (j * self.rho1) ** 2 - float(yk @ yk) + tol
        if R2 < 0:
            return 1.0, -1.0
        Rr = math.sqrt(R2)
        c = j / self.d + Qk @ yk + tol
        A = np.hypot(Qi, w)
        lo, hi = -Rr, Rr
        with np.errstate(divide="ignore", invalid="ignore"):
            kappa = np.where(A > 0, -c / (Rr * A), np.where(c >= 0, -2.0, 2.0))
        if np.any(kappa > 1):
            return 1.0, -1.0
        phi = np.arctan2(w, Qi)
        act = kappa > -1
        if act.any():
            ac = np.arccos(kappa[act])
            th_hi = np.minimum(math.pi, phi[act] + ac)
            th_lo = np.maximum(0.0, phi[act] - ac)
            lo = max(lo, float(np.max(Rr * np.cos(th_hi))))
            hi = min(hi, float(np.min(Rr * np.cos(th_lo))))
        return lo, hi

    def enumerate(
        self,
        j: int,
        mode: str = "list",
        limit: Optional[int] = None,
        node_limit: Optional[int] = None,
        time_limit: Optional[float] = None,
        rigorous: bool = False,
    ) -> EnumerationResult:
        """Enumerate ``L cap Delta_j``.

        ``mode``: ``"list"`` (return integer coordinate vectors), ``"count"``
        (count only) or ``"first"`` (stop at the first point).  ``limit``
        stops after that many points.  ``rigorous`` verifies every point with
        the exact callback.
        """
        if mode not in ("list", "count", "first"):
            raise ValueError("mode must be 'list', 'count' or 'first'")
        start = time.time()
        res = EnumerationResult(level=j, count=0)
        if j <= 0 or j % self.gcd:
            res.seconds = time.time() - start
            return res
        if mode == "first":
            limit = 1
        want_points = mode != "count"
        if self.n == 0:
            a = [p * (j // self.gcd) for p in self.p1]
            res.count = 1
            if want_points:
                res.points.append(a)
            return res

        s, a0 = self._level(j)
        n, d, R, Q = self.n, self.d, self.R, self.Q
        tol = 1e-9 * (j + float(np.max(np.abs(s))))
        delta = 1e-7 * j
        verify = self.verify
        k = [0] * n
        y = np.zeros(n)
        deadline = None if time_limit is None else start + time_limit

        def check_exact(kvec) -> bool:
            res.exact_checks += 1
            if verify is None:
                return True
            return bool(verify(self.point_coords(a0, kvec)))

        def emit(kvec):
            if want_points:
                res.points.append(self.point_coords(a0, kvec))
            res.count += 1
            if limit is not None and res.count >= limit:
                raise _Abort

        def order(lo_k: int, hi_k: int, centre: float):
            if mode != "first":
                return range(lo_k, hi_k + 1)
            kc = min(max(int(round(centre)), lo_k), hi_k)
            out = [kc]
            step = 1
            while kc - step >= lo_k or kc + step <= hi_k:
                if kc + step <= hi_k:
                    out.append(kc + step)
                if kc - step >= lo_k:
                    out.append(kc - step)
                step += 1
            return out

        def leaf(c: float, klo: int, khi: int):
            # x_r(kk) = j/d + Q[r,0] (c + R00 kk) + Q[r,1:] . y[1:]
            base = j / d + Q[:, 0] * c + Q[:, 1:] @ y[1:]
            slope = Q[:, 0] * R[0, 0]
            slo, shi = klo, khi
            with np.errstate(divide="ignore"):
                bnd = (delta - base) / slope
            pos, neg = slope > 0, slope < 0
            if pos.any():
                slo = max(slo, math.ceil(float(np.max(bnd[pos]))))
            if neg.any():
                shi = min(shi, math.floor(float(np.min(bnd[neg]))))
            zero = ~(pos | neg)
            if zero.any() and np.min(base[zero]) < delta:
                slo, shi = khi + 1, khi
            bulk = not want_points and limit is None and not rigorous
            if bulk:
                if shi >= slo:
                    unsure = list(range(klo, slo)) + list(range(shi + 1, khi + 1))
                else:
                    unsure = list(range(klo, khi + 1))
                for kk in unsure:
                    k[0] = kk
                    if check_exact(k):
                        res.count += 1
                    else:
                        res.rejected += 1
                if shi >= slo:
                    res.count += shi - slo + 1
                return
            for kk in order(klo, khi, (klo + khi) / 2):
                k[0] = kk
                if rigorous or not (slo <= kk <= shi):
                    if not check_exact(k):
                        res.rejected += 1
                        continue
                emit(k)

        def rec(i: int, partial: np.ndarray):
            res.nodes += 1
            if node_limit is not None and res.nodes > node_limit:
                res.complete = False
                raise _Abort
            if deadline is not None and (res.nodes & 1023) == 0 and time.time() > deadline:
                res.complete = False
                raise _Abort
            c = float(partial[i])
            lo, hi = self._range(i, y, j, tol)
            if lo > hi:
                return
            rii = R[i, i]
            klo = math.ceil((lo - tol - c) / rii)
            khi = math.floor((hi + tol - c) / rii)
            if klo > khi:
                return
            if i == 0:
                leaf(c, klo, khi)
                return
            col = R[:i, i]
            for kk in order(klo, khi, ((lo + hi) / 2 - c) / rii):
                k[i] = kk
                y[i] = c + rii * kk
                rec(i - 1, partial[:i] + col * kk)
            k[i] = 0

        try:
            rec(n - 1, s.copy())
        except _Abort:
            pass
        res.seconds = time.time() - start
        return res


@dataclass
class MaxMinResult:
    """Result of :meth:`SimplexLattice.max_min`.

    ``value`` approximates ``max { min_i x_i : x in L, sum x = j }`` (the
    *margin*); the level is attained (``L cap Delta_j`` nonempty) iff it is
    positive.  ``point`` is a maximiser (integer coordinates) and
    ``verified`` tells whether its total positivity was confirmed exactly.
    """

    level: int
    value: float
    point: Optional[list]
    verified: Optional[bool]
    status: str
    seconds: float

    @property
    def exists(self) -> Optional[bool]:
        if self.verified:
            return True
        if self.value < -1e-6 * max(1, self.level):
            return False
        return None  # undecided by the MILP (|value| within tolerance)


class _Silence:
    """Silence C-level stdout (HiGHS prints debugging lines)."""

    def __enter__(self):
        import os
        import sys

        sys.stdout.flush()
        self._fd = os.dup(1)
        self._null = os.open(os.devnull, os.O_WRONLY)
        os.dup2(self._null, 1)
        return self

    def __exit__(self, *exc):
        import os

        os.dup2(self._fd, 1)
        os.close(self._fd)
        os.close(self._null)
        return False


def _max_min(self: SimplexLattice, j: int, time_limit: Optional[float] = None) -> MaxMinResult:
    """Maximise ``min_i x_i`` over the affine lattice ``L cap {sum = j}`` (MILP).

    Solved with HiGHS (``scipy.optimize.milp``) in the LLL-reduced
    coordinates.  A positive optimum is turned into a rigorous witness by
    exact verification of the maximiser; a clearly negative optimum is a
    (floating point, branch-and-bound) proof that level ``j`` is empty.
    """
    from scipy.optimize import Bounds, LinearConstraint, milp

    start = time.time()
    if j <= 0 or j % self.gcd:
        return MaxMinResult(j, -math.inf, None, None, "empty level", 0.0)
    if self.n == 0:
        a = [p * (j // self.gcd) for p in self.p1]
        return MaxMinResult(j, float(j), a, True, "dimension 0", 0.0)
    s, a0 = self._level(j)
    tf = self._tf_cache[j] if j in self._tf_cache else self.Q @ s + j / self.d
    d, n = self.d, self.n
    c = np.zeros(n + 1)
    c[-1] = -1.0
    A = np.hstack([-self.Bf, np.ones((d, 1))])  # t - B k <= tf
    cons = LinearConstraint(A, -np.inf, tf)
    integrality = np.array([1] * n + [0])
    # finite box for k: x_r >= -j and sum x = j give |x - tf| <= |x| + |tf|
    xmax = (d + 1) * j * math.sqrt(d) + float(np.linalg.norm(tf))
    Rinv = np.linalg.inv(self.R)
    kb = np.ceil(np.linalg.norm(Rinv, axis=1) * xmax) + 1
    lb = np.concatenate([-kb, [-float(j)]])
    ub = np.concatenate([kb, [float(j)]])
    r = None
    for presolve in (True, False):
        opts = {"disp": False, "presolve": presolve}
        if time_limit is not None:
            opts["time_limit"] = time_limit
        with _Silence():
            r = milp(c, constraints=[cons], integrality=integrality, bounds=Bounds(lb, ub), options=opts)
        if r.x is not None:
            break
    if r.x is None:
        return MaxMinResult(j, math.nan, None, None, r.message, time.time() - start)
    k = [int(round(v)) for v in r.x[:n]]
    x = tf + self.Bf @ np.array(k, dtype=float)
    value = float(np.min(x))
    point = self.point_coords(a0, k)
    verified = None
    if value > -1e-6 * j and self.verify is not None:
        verified = bool(self.verify(point))
    status = "optimal" if r.status == 0 else r.message
    return MaxMinResult(j, value, point, verified, status, time.time() - start)


SimplexLattice.max_min = _max_min
