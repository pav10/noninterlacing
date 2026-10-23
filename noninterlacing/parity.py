"""Arrays of elements of ``F_n = Q(zeta_n)^+`` and the three-block parity identity.

For ``delta`` in ``F_n`` the array is ``c(x) = Tr_{Q(zeta_n)/Q}(delta zeta_n^x)`` for ``x`` in
``Z/n``; it equals ``Tr_{F_n/Q}(delta C_x(alpha))`` with ``C_x(2cos t) = 2cos(xt)``.  For
pairwise coprime blocks ``n = n_1 n_2 n_3`` the coordinates are ``x mod n_i`` (CRT), the
trace to ``F_{n_1 n_2}`` is the restriction to the plane ``x_3 = 0``, and products of elements
of coprime conductors have tensor-product arrays.  See ``notes/parity.tex``.
"""

from __future__ import annotations

from collections import Counter

from flint import fmpq_poly

from .fibers import CyclotomicField, tower_levels


def power_sums(K: CyclotomicField) -> list:
    """``Tr_{F_n/Q}(alpha^i)`` for ``0 <= i < d`` (Newton's identities)."""
    c, d = [int(a) for a in K.f.coeffs()], K.d
    s = [d]
    for k in range(1, d):
        v = -k * c[d - k]
        for i in range(1, k):
            v -= c[d - i] * s[k - i]
        s.append(v)
    return s


def trace_array(K: CyclotomicField, x, sums: list = None) -> list:
    """The array ``c(m) = Tr_{F_n/Q}(x C_m(alpha))``, ``m = 0, ..., n-1``; integral iff ``x`` is in the codifferent."""
    s = sums or power_sums(K)
    out = []
    for m in range(K.n):
        y = (fmpq_poly(x) * K.cheb(m)) % K.F
        t = sum(y[i] * s[i] for i in range(K.d))
        if t.q != 1:
            raise ArithmeticError("x is not in the codifferent of F_n")
        out.append(int(t.p))
    return out


def block_array(c: list, blocks: tuple) -> dict:
    """``c`` keyed by the CRT coordinates ``(m mod n_1, ..., m mod n_k)``."""
    return {tuple(m % q for q in blocks): v for m, v in enumerate(c)}


def reps(q: int) -> list:
    """``T = {1, ..., (q-1)/2}``: ``Z/q minus 0`` is ``T`` together with ``-T``."""
    return list(range(1, (q - 1) // 2 + 1))


def parity_terms(C: dict, blocks: tuple) -> tuple:
    """``(X, Y, Z)`` of the parity identity; ``X + Y + Z`` is even for every array of the codifferent.

    ``X``, ``Y``, ``Z`` are the sums of ``C`` over ``0 x T_2 x T_3``, ``T_1 x 0 x (-T_3)`` and
    ``T_1 x T_2 x 0``.
    """
    n1, n2, n3 = blocks
    T1, T2, T3 = reps(n1), reps(n2), reps(n3)
    X = sum(C[(0, y, z)] for y in T2 for z in T3)
    Y = sum(C[(x, 0, (-z) % n3)] for x in T1 for z in T3)
    Z = sum(C[(x, y, 0)] for x in T1 for y in T2)
    return X, Y, Z


def corner_parity(C: dict, m1: int, m2: int) -> int:
    """``tau(C) = sum of C over T_1 x T_2`` mod 2 (independent of the ``T_i`` when the axes are even)."""
    return sum(C[(x, y)] for x in reps(m1) for y in reps(m2)) % 2


def _orbit_elements(K: CyclotomicField, stage) -> list:
    elems = {}
    for mu in stage.reps:
        for a in K.G:
            x = K.act(a, mu)
            elems[str(x)] = x
    if len(elems) != stage.count:
        raise ArithmeticError("orbit expansion does not match the count")
    return list(elems.values())


def level_two_census(p: int, q: int) -> dict:
    """Level-2 elements of ``F_pq`` sorted by their traces to ``F_p``, ``F_q`` and by ``tau``.

    An element has *doubled* traces if ``Tr_{F_pq/F_p} beta = 2a`` and ``Tr_{F_pq/F_q} beta = 2b``
    with ``a``, ``b`` level-1 elements; it is then a product ``ab`` iff its array is
    ``(2c_a) (x) (2c_b) / 4``.  Returns the total and the counts of the classes
    ``"product"``, ``"tau=0"``, ``"tau=1"`` (non-products with doubled traces) and ``"other"``.
    """
    K = CyclotomicField(p * q)
    s = power_sums(K)

    def axis_set(r: int, other: int) -> set:
        # level-1 elements a of F_r; in F_pq their array is c_a(x) Tr(zeta_other^y), so c_a = axis / (other - 1)
        H = K.fixer(r)
        elems = [fmpq_poly([1])] if len(H) == K.d else _orbit_elements(K, tower_levels(K, [K.G, H], 1)[-1])
        out = set()
        for a in elems:
            C = block_array(trace_array(K, a, s), (p, q))
            ax = [C[(x, 0)] if r == p else C[(0, x)] for x in range(r)]
            out.add(tuple(v // (other - 1) for v in ax))
        return out

    L1p, L1q = axis_set(p, q), axis_set(q, p)
    st = tower_levels(K, [K.G, K.fixer(q), [1]], 2)[-1]
    stats = Counter()
    for beta in _orbit_elements(K, st):
        C = block_array(trace_array(K, beta, s), (p, q))
        A = [C[(x, 0)] for x in range(p)]
        B = [C[(0, y)] for y in range(q)]
        doubled = (all(v % 2 == 0 for v in A + B) and tuple(v // 2 for v in A) in L1p
                   and tuple(v // 2 for v in B) in L1q)
        if not doubled:
            stats["other"] += 1
        elif all(4 * C[(x, y)] == A[x] * B[y] for x in range(p) for y in range(q)):
            stats["product"] += 1
        else:
            stats[f"tau={corner_parity(C, p, q)}"] += 1
    return {"field": f"F_{p * q}", "degree": K.d, "level2": st.count,
            **{k: stats[k] for k in ("product", "tau=0", "tau=1", "other")}}


def three_q_family(q: int) -> list:
    """``a(1 +- 1/theta)`` in ``F_3q``, ``a = (2 - alpha_q)/q``, ``theta = (z_3 - 1/z_3)(z_q - 1/z_q)``.

    Returns the two elements in the coordinates of ``F_3q``; they lie in the codifferent and
    have trace 2, and are totally positive iff ``12 sin^2(pi/q) > 1``, i.e. ``q <= 7``.
    """
    n = 3 * q
    K = CyclotomicField(n)
    al = fmpq_poly([0, 1])
    h = next(a for a in range(1, n) if a % 3 == 2 and a % q == 1)  # fixes F_q, inverts zeta_3
    theta = al - K.act(h, al)
    _, inv, _ = theta.xgcd(K.F)
    a = (2 - K.embed(al, q)) / q
    return [(a * (1 + s * inv)) % K.F for s in (1, -1)]
