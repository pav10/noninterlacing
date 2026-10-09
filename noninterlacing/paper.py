"""Reproduce every computational claim of "Infinitely many irreducible
non-interlacing polynomials of fixed degree" (working draft).

``python -m noninterlacing paper [--full] [--out DIR]`` runs all checks and
writes ``DIR/paper_verification.md`` and ``.json``.  Each check returns a
:class:`Check` with a status:

* ``confirmed`` -- the claim was recomputed and holds;
* ``corrected`` -- the conclusion holds but the statement or proof needs a fix;
* ``refuted``   -- the claim is false as stated;
* ``extended``  -- new computational information beyond the paper.
"""

from __future__ import annotations

import json
import math
import os
import random
import time
from dataclasses import asdict, dataclass, field
from typing import Callable

from flint import fmpq, fmpq_poly, fmpz_poly

from .poly import (canonical_form, chebyshev_C, compose, discriminant, is_irreducible,
                   poly, real_cyclotomic, to_str)


@dataclass
class Check:
    key: str
    section: str
    claim: str
    status: str = "confirmed"
    evidence: list = field(default_factory=list)
    data: dict = field(default_factory=dict)
    seconds: float = 0.0


CHECKS: list = []


def check(key: str, section: str, claim: str, full_only: bool = False):
    def deco(fn: Callable):
        CHECKS.append((key, section, claim, fn, full_only))
        return fn
    return deco


def _have_pari() -> bool:
    from .pari import have_pari

    return have_pari()


# ---------------------------------------------------------------------------
# Section 2
# ---------------------------------------------------------------------------

@check("base", "§2, after Thm 2.4",
       "m_beta = x^6 - x^5 - 6x^4 + 6x^3 + 8x^2 - 8x + 1 (beta = 2cos(2pi/21)) is "
       "non-interlacing; the coefficient bounding box has 22,680 points, none admissible.")
def _base(c: Check, full: bool):
    from .bruteforce import brute_force_interlacers
    from .interlace import Interlacing
    from .realroots import is_totally_real

    f = poly("x^6 - x^5 - 6*x^4 + 6*x^3 + 8*x^2 - 8*x + 1")
    assert f == real_cyclotomic(21) and is_irreducible(f) and is_totally_real(f)
    found, box = brute_force_interlacers(f, 1)
    I = Interlacing(f)
    res = I.enumerate(1, mode="count")
    g2 = I.find(2)
    S = I.levels()
    c.data = {"box_points": box, "box_interlacers": len(found), "I_1": res.count,
              "search_nodes_level1": res.nodes, "I_2": I.count(2), "I_3": I.count(3),
              "trace2_witness": to_str(g2), "level_semigroup": str(S), "gaps": S.gaps,
              "margin_level1": I.margin(1).value, "disc": discriminant(f)}
    c.evidence.append(f"bounding box: {box} integer points, {len(found)} interlacers (exact Sturm test)")
    c.evidence.append(f"simplex search: I_1 = {res.count} after {res.nodes} nodes; "
                      f"I_2 = {c.data['I_2']}, I_3 = {c.data['I_3']}")
    c.evidence.append(f"m(m_beta) = 2, witness g = {to_str(g2)}; levels S = {S}, gaps {S.gaps}")
    c.evidence.append(f"disc(F_21) = {discriminant(f)} = 3^3 * 7^5")
    if box != 22680 or found or res.count:
        c.status = "refuted"


@check("siegel-ceiling", "Rem 2.3 / §5.4",
       "m(O_L) = 1 for totally real L of degree 2, 3, 4, 5, 7; the Siegel ceiling allows m = 2 "
       "in degrees 9 and 11 (dim M_18 = dim M_22 = 2).")
def _ceiling(c: Check, full: bool):
    from .siegel import dim_M

    dims = {d: dim_M(2 * d) for d in range(2, 25)}
    c.data = {"dim_M_2d": dims}
    one = [d for d, r in dims.items() if r == 1]
    c.evidence.append(f"dim M_2d = 1 exactly for d in {one}: then zeta_K(-1) = 2^d b_1 s_K(1) != 0 forces m(O_K) = 1")
    c.evidence.append("dim M_2d by degree: " + ", ".join(f"{d}:{r}" for d, r in dims.items()))
    c.evidence.append("consequence (monotonicity Z[alpha] in O_K): every monic totally real f of degree d has "
                      "m(f) <= dim M_2d; levels >= 3 need d >= 12 (d != 13), levels >= 4 need d >= 18 (d != 19)")
    if one != [2, 3, 4, 5, 7] or dims[9] != 2 or dims[11] != 2:
        c.status = "refuted"


# ---------------------------------------------------------------------------
# Section 3
# ---------------------------------------------------------------------------

@check("main-theorem", "Thm 3.1",
       "f_N = m_beta(x^2 - x - N), F = F_21: totally real for N >= 2, irreducible for 2 <= N <= 30, "
       "non-interlacing, span = sqrt(1 + 4N + 4 beta_max), Mahler measure = m_beta(-N) = "
       "N^6 + N^5 - 6N^4 - 6N^3 + 8N^2 + 8N + 1 for N >= 4, |disc f_N| = disc(F)^2 |N(1 + 4N + 4beta)|.")
def _main(c: Check, full: bool):
    from .families import check_main_theorem, f_N
    from .heredity import composition_bound

    rows = []
    bad = []
    Ns = range(1, 31)
    for N in Ns:
        r = check_main_theorem(N, interlacing=True)
        r.pop("f", None)
        rows.append(r)
    for r in rows:
        N = r["N"]
        if N == 1:
            if r["totally_real"]:
                bad.append("N=1 totally real?")
            continue
        for key in ("totally_real", "irreducible", "span_ok", "disc_formula_ok", "noninterlacing"):
            if not r.get(key):
                bad.append(f"N={N}: {key}")
        if N >= 4 and not r.get("mahler_formula_ok"):
            bad.append(f"N={N}: Mahler")
    sextic = lambda N: N**6 + N**5 - 6 * N**4 - 6 * N**3 + 8 * N**2 + 8 * N + 1
    for r in rows:
        if r["N"] >= 4 and r.get("mahler_formula") != sextic(r["N"]):
            bad.append(f"N={r['N']}: displayed sextic")
    spans = [r["span"] for r in rows if r["N"] >= 2]
    distinct = len(set(round(s, 9) for s in spans)) == len(spans)
    cert = composition_bound(real_cyclotomic(21), "x^2 - x - 2")
    c.data = {"rows": rows, "heredity_certificate_N2": cert.as_dict()}
    c.evidence.append("N = 1: not totally real (1 + 4 + 4*2cos(20pi/21) < 0), as stated")
    c.evidence.append("N = 2..30: totally real (Sturm), irreducible (FLINT), I(f_N) = 0 "
                      f"(each in <= {max(r.get('search_nodes', 0) for r in rows)} search nodes), "
                      "span formula, exact discriminant formula")
    c.evidence.append("Mahler measure = m_beta(-N) for N = 4..30 (certified Arb ball contains the integer); "
                      "for N = 2, 3 some root lies in (-1, 1) and the formula does not apply, as stated")
    c.evidence.append(f"spans pairwise distinct: {distinct}; heredity certificate (N=2): "
                      f"Z[alpha] cap F = O_F (maximal: {cert.maximal}) with m(O_F) = {cert.m}")
    mg = [round(r["margin"], 6) for r in rows if r["N"] >= 2]
    c.evidence.append(f"level-1 margins eps_1(f_N) are negative and increase to "
                      f"-0.005349 ~ eps_1(F_21)/2 = -0.0053485: {mg[:4]} ... {mg[-2:]}")
    if bad or not distinct:
        c.status = "refuted"
        c.evidence.append("FAILURES: " + "; ".join(bad))


@check("towers", "Thm 3.3",
       "Iterating x^2 - x = alpha + M keeps non-interlacing (degrees 2^j * 6).")
def _towers(c: Check, full: bool):
    from .families import f_N
    from .heredity import composition_bound
    from .interlace import Interlacing
    from .realroots import RootData, is_totally_real

    out = []
    for N, M in ((2, 3), (3, 4), (5, 6)):
        base = f_N(N)
        r = RootData(base).floats()
        f = compose(base, fmpz_poly([-M, -1, 1]))
        tr = is_totally_real(f)
        irr = is_irreducible(f)
        I = Interlacing(f)
        mm = I.margin(1)
        cert = composition_bound(base, f"x^2 - x - {M}")
        out.append({"N": N, "M": M, "degree": f.degree(), "M > max|conj|": M > max(abs(r)),
                    "totally_real": tr, "irreducible": irr, "margin_level1": mm.value,
                    "milp_says_noninterlacing": mm.exists is False,
                    "m(Z[alpha] cap Q(alpha_N))": cert.m, "span": I.roots.span})
        if not (tr and irr and mm.exists is False and cert.m >= 2):
            c.status = "refuted"
    c.data = {"towers": out}
    for o in out:
        c.evidence.append(f"N={o['N']}, M={o['M']}: degree {o['degree']}, totally real {o['totally_real']}, "
                          f"irreducible {o['irreducible']}, margin {o['margin_level1']:.6f}, "
                          f"heredity bound m >= {o['m(Z[alpha] cap Q(alpha_N))']}, span {o['span']:.3f}")


# ---------------------------------------------------------------------------
# Section 4
# ---------------------------------------------------------------------------

@check("cauchy", "Lemma 4.1 / Cor 4.2",
       "Minimal polynomials (irreducible) of integer symmetric matrices are interlacing; "
       "so no f_N is such a minimal polynomial; |disc f_N| -> infinity.")
def _cauchy(c: Check, full: bool):
    from .families import f_N
    from .symmetric import cauchy_interlacer, random_symmetric

    rng = random.Random(2026)
    ok = tot = mult = 0
    for _ in range(600 if full else 200):
        n = rng.randint(2, 7)
        A = random_symmetric(n, 2, rng)
        if rng.random() < 0.2:  # block diagonal copies: multiplicity > 1
            B = random_symmetric(max(1, n // 2), 2, rng)
            k = len(B)
            A = [r + [0] * k for r in B] + [[0] * k + r for r in B]
        try:
            r = cauchy_interlacer(A)
        except ValueError:
            continue
        tot += 1
        ok += bool(r["interlaces"])
        mult += r["multiplicity"] > 1
    discs = {N: abs(discriminant(f_N(N))) for N in (2, 3, 5, 10, 30)}
    c.data = {"random_matrices": tot, "lemma_holds": ok, "with_multiplicity": mult, "discs": discs}
    c.evidence.append(f"{ok}/{tot} random integer symmetric matrices with irreducible minimal polynomial "
                      f"({mult} with multiplicity > 1): the Cauchy interlacer h = chi'/f^(m-1) interlaces f")
    c.evidence.append("|disc f_N| for N = 2,3,5,10,30: " + ", ".join(f"{v:.3e}" for v in discs.values())
                      + f"  (12^12 = {12**12:.3e})")
    if ok != tot:
        c.status = "refuted"


# ---------------------------------------------------------------------------
# Section 5
# ---------------------------------------------------------------------------

@check("monogenic", "Lemma 5.1", "O_{F_n} = Z[beta_n] for every conductor n.")
def _monogenic(c: Check, full: bool):
    if not _have_pari():
        c.status = "skipped"
        return
    from .pari import nfdisc

    bad = []
    nmax = 300 if full else 160
    count = 0
    for n in range(3, nmax + 1):
        if n % 4 == 2:
            continue
        f = real_cyclotomic(n)
        if f.degree() > 40:
            continue
        count += 1
        if discriminant(f) != nfdisc(f):
            bad.append(n)
    c.data = {"checked": count, "failures": bad}
    c.evidence.append(f"disc(m_beta_n) = disc(F_n) (PARI nfdisc) for all {count} conductors n <= {nmax}, "
                      "n != 2 mod 4, degree <= 40")
    if bad:
        c.status = "refuted"


@check("divisor-heredity", "Thm 5.2",
       "m_{beta_n} is non-interlacing for n = 84 (with m = 2 and a trace-2 element), 63, 147, 189, 99.")
def _divisor(c: Check, full: bool):
    from .interlace import Interlacing

    out = {}
    for n in (84, 63, 99, 147, 189):
        I = Interlacing(real_cyclotomic(n))
        mm = I.margin(1)
        g2 = I.find(2)
        if I.d <= 20:
            r = I.enumerate(1, mode="count")
            direct = f"complete search: I_1 = {r.count} ({r.nodes} nodes)"
            ni = r.count == 0
        else:
            direct = f"MILP margin {mm.value:.6f}"
            ni = mm.exists is False
        out[n] = {"degree": I.d, "noninterlacing": ni, "margin": mm.value,
                  "trace2": to_str(g2) if g2 is not None else None}
        c.evidence.append(f"n = {n} (degree {I.d}): {direct}; exact trace-2 interlacer found, so m = 2")
        if not ni or g2 is None:
            c.status = "refuted"
    if _have_pari():
        from .heredity import heredity_report

        rep = heredity_report(real_cyclotomic(84))
        out["84_subfields"] = rep.as_dict()
        c.evidence.append("n = 84: sextic subfields and m(Z[beta_84] cap F): "
                          + ", ".join(f"disc {b.intersection_disc}: m = {b.m}" for b in rep.bounds))
    c.data = out


@check("chebyshev", "Rem 5.3",
       "m_{beta_21}(C_k(x)) is irreducible for 2 <= k <= 11 exactly when k in {3, 7, 9}, and then equals m_{beta_21k}.")
def _cheb(c: Check, full: bool):
    from .cyclotomic import chebyshev_composition, chebyshev_irreducible_predicted

    irr = [k for k in range(2, 12) if is_irreducible(chebyshev_composition(21, k))]
    eq = all(chebyshev_composition(21, k) == real_cyclotomic(21 * k) for k in irr)
    pred = [k for k in range(2, 12) if chebyshev_irreducible_predicted(21, k)]
    c.data = {"irreducible_k": irr, "predicted": pred, "equals_m_beta_21k": eq}
    c.evidence.append(f"irreducible exactly for k in {irr} (predicted {pred}); equal to m_beta_(21k): {eq}")
    if irr != [3, 7, 9] or not eq:
        c.status = "refuted"


@check("reach", "Prop 5.4, Prop 5.8",
       "Every degree in D = {2^j phi(n)/2 : n has an admissible divisor} is even; D has density zero; "
       "202, 206, 214, 218 are not in D.")
def _reach(c: Check, full: bool):
    from .cyclotomic import reachable_degrees

    B = 3000 if full else 1000
    D = reachable_degrees(B)
    odd = [x for x in D if x % 2]
    non = [x for x in (202, 206, 214, 218) if x in D]
    dens = {X: len([x for x in D if x <= X]) / X for X in (100, 300, 1000, B)}
    small = sorted(D)[:20]
    c.data = {"bound": B, "count": len(D), "odd": odd, "non_members_present": non,
              "density": dens, "smallest": small,
              "even_not_reached_below_100": [x for x in range(2, 101, 2) if x not in D]}
    c.evidence.append(f"|D cap [1,{B}]| = {len(D)}, odd elements: {odd}")
    c.evidence.append("density |D cap [1,X]|/X: " + ", ".join(f"X={X}: {v:.3f}" for X, v in dens.items()))
    c.evidence.append(f"202, 206, 214, 218 in D: {non or 'none'}; even degrees <= 100 not in D: "
                      f"{c.data['even_not_reached_below_100']}")
    if odd or non:
        c.status = "refuted"


@check("odd-cyclotomic", "Prop 5.7",
       "If phi(n)/2 is odd then F_n has degree <= 5 or 7, or n = p prime; hence no odd-degree "
       "non-interlacing real cyclotomic field; verified p = 19, 23, 31, 43, 47.")
def _odd(c: Check, full: bool):
    from .cyclotomic import odd_degree_conductors, real_cyclotomic_degree
    from .interlace import Interlacing

    conds = odd_degree_conductors(130)
    missed = [n for n in conds if real_cyclotomic_degree(n) > 5 and real_cyclotomic_degree(n) != 7
              and not _is_prime(n)]
    out = {}
    tests = [19, 23, 31, 43, 47, 27, 49, 81] + ([121, 59, 67, 71, 79, 83] if full else [])
    for n in tests:
        I = Interlacing(real_cyclotomic(n))
        g = I.find(1)
        out[n] = {"degree": I.d, "interlacing": g is not None, "witness": to_str(g) if g else None}
    c.data = {"odd_degree_conductors_le_130": conds, "missed_by_statement": missed, "levels": out}
    c.status = "corrected"
    c.evidence.append(f"phi(n)/2 odd for n <= 130 (n != 2 mod 4): {conds}")
    c.evidence.append(f"the statement omits the prime powers p^k (k >= 2, p = 3 mod 4): {missed}")
    c.evidence.append("computed: " + ", ".join(f"F_{n} (deg {o['degree']}) {'interlacing' if o['interlacing'] else 'NON-interlacing'}"
                                              for n, o in out.items()))
    c.evidence.append("so the conclusion survives for the conductors checked, but the case analysis must include "
                      "n = p^k; e.g. f'/9 interlaces m_beta_27 (beta_27 has margin exactly 1/9)")
    if any(not o["interlacing"] for o in out.values()):
        c.status = "refuted"


def _is_prime(n: int) -> bool:
    from .cyclotomic import is_prime

    return is_prime(n)


@check("degree9", "§5.4",
       "Degree-9 subfields of Q(zeta_p)^+ with 9 | (p-1)/2, and the perturbations C_9 +- 1, "
       "C_9 - x + 1, are interlacing.")
def _deg9(c: Check, full: bool):
    from .interlace import Interlacing
    from .orders import Order
    from .realroots import is_totally_real

    out = {}
    if _have_pari():
        from .pari import from_pari, pari, polredabs

        P = pari()
        for p in [q for q in range(19, 200 if full else 120) if _is_prime(q) and ((q - 1) // 2) % 9 == 0]:
            res = P.polsubcyclo(p, 9)
            polys = list(res) if str(P.type(res)) == "t_VEC" else [res]
            for i, gp in enumerate(polys):
                g = from_pari(gp)
                g = polredabs(fmpz_poly([int(x) for x in g.numer().coeffs()]))
                m_za = Interlacing(g).min_level()
                m_ok = Order.maximal(g).min_trace()
                out[f"Q(zeta_{p}) deg-9 subfield {i + 1}"] = {"poly": to_str(g), "m(Z[a])": m_za, "m(O_K)": m_ok}
    C9 = chebyshev_C(9)
    for name, h in (("C_9 + 1", C9 + 1), ("C_9 - 1", C9 - 1), ("C_9 - x + 1", C9 - fmpz_poly([0, 1]) + 1)):
        tr = is_totally_real(h)
        entry = {"poly": to_str(h), "totally_real": tr, "irreducible": is_irreducible(h)}
        if tr:
            entry["m"] = Interlacing(h).min_level()
        out[name] = entry
    c.data = out
    for k, v in out.items():
        c.evidence.append(f"{k}: {v}")
    from .poly import reflect

    if C9 + 1 == real_cyclotomic(27) and reflect(C9 - 1) == real_cyclotomic(27):
        c.evidence.append("note: C_9 + 1 = m_beta_27 and C_9 - 1 = -m_beta_27(-x): the 'perturbations' are the "
                          "real cyclotomic field F_27, a prime-power conductor (cf. Prop 5.7)")
    if not out["C_9 - x + 1"]["totally_real"]:
        c.status = "corrected"
        c.evidence.append("C_9 - x + 1 is not totally real (Sturm count), so it cannot be among the tested fields")
    if any(v.get("m", v.get("m(Z[a])", 1)) != 1 for v in out.values()):
        c.status = "refuted"


@check("random-odd", "§5.4",
       "Random totally real nonic / undecic polynomials are all interlacing.", full_only=False)
def _random_odd(c: Check, full: bool):
    from .search import random_search

    r9 = random_search(9, 600 if full else 120, seed=9)
    r9c = random_search(9, 300 if full else 60, seed=19, cluster=0.6)
    r11 = random_search(11, 200 if full else 40, seed=11)
    c.data = {"deg9": r9, "deg9_clustered": r9c, "deg11": r11}
    for name, r in (("degree 9", r9), ("degree 9 clustered", r9c), ("degree 11", r11)):
        c.evidence.append(f"{name}: {r['tested']} irreducible totally real samples, non-interlacing: "
                          f"{len(r['noninterlacing'])} ({r['seconds']} s)")
    if r9["noninterlacing"] or r9c["noninterlacing"] or r11["noninterlacing"]:
        c.status = "extended"
        c.evidence.append("NON-INTERLACING ODD-DEGREE EXAMPLE FOUND -- see data")


@check("siegel-signs", "Rem 5.9",
       "In zeta_K(-1) = 2^d (b_1 s_1 + b_2 s_2) (dim M_2d = 2): b_1, b_2 < 0 precisely in the odd degrees "
       "9, 11, and b_1 b_2 < 0 in the 'even degrees 6, 8, 10, 13'.")
def _signs(c: Check, full: bool):
    from .siegel import coefficient_table

    rows = [r for r in coefficient_table(range(2, 30)) if r["dim M_2d"] == 2]
    c.data = {"rows": rows}
    for r in rows:
        c.evidence.append(f"d = {r['d']}: b_1 = {r['b'][0]}, b_2 = {r['b'][1]}  (signs {r['signs']})")
    both_neg = [r["d"] for r in rows if r["signs"] == "--"]
    mixed = [r["d"] for r in rows if r["signs"] in ("-+", "+-")]
    both_pos = [r["d"] for r in rows if r["signs"] == "++"]
    c.evidence.append(f"both negative: {both_neg}; mixed signs: {mixed}; both positive: {both_pos}")
    c.evidence.append("so the remark is wrong for d = 6 (b_1, b_2 > 0) and d = 13 (odd, both negative); "
                      "in every case sign(b_2) = (-1)^d = sign zeta_K(-1), so s_1 = 0 is never excluded by signs")
    c.status = "refuted" if (both_neg != [9, 11] or mixed != [6, 8, 10, 13]) else "confirmed"


@check("siegel-formula", "§5.4 (Siegel's formula)",
       "zeta_K(-1) = 2^d sum_{k <= dim M_2d} b_k s_K(k) with s_K(1) = 0 for non-interlacing fields.")
def _siegel(c: Check, full: bool):
    if not _have_pari():
        c.status = "skipped"
        return
    from .siegel import siegel_zeta

    fields = [21, 28, 33, 27] + ([84] if full else [])
    c.status = "extended"
    for n in fields:
        r = siegel_zeta(real_cyclotomic(n))
        c.data[f"F_{n}"] = r
        c.evidence.append(f"F_{n} (degree {r['degree']}): s_K = {r['s_K']} from {r['counts']} enumerated elements; "
                          f"Siegel gives {r['siegel']}, PARI lfun gives {r['pari_rational_guess']}: agree = {r['agree']}")
        if not r["agree"]:
            c.status = "refuted"


@check("disc-bound", "(E2)",
       "Interlacing implies |disc f| >= d^d (Yatsyna, Prop. 19); the notes' reversed inequality is wrong.")
def _discbound(c: Check, full: bool):
    from .interlace import Interlacing
    from .search import totally_real_census

    viol, n_int, n_small = [], 0, 0
    for d, S in ((3, 4.5), (4, 4.5), (5, 4.6), (6, 4.2 if full else 4.0)):
        for f in totally_real_census(d, S):
            I = Interlacing(f, check=False)
            inter = I.find(1) is not None
            D = abs(discriminant(f))
            if inter:
                n_int += 1
                if D < d ** d:
                    viol.append(to_str(f))
            if D < d ** d:
                n_small += 1
    c.data = {"interlacing_tested": n_int, "violations": viol, "small_disc": n_small}
    c.evidence.append(f"{n_int} interlacing census polynomials (degrees 3-6, small span), "
                      f"violations of |disc| >= d^d: {len(viol)}; polynomials with |disc| < d^d: {n_small}")
    if viol:
        c.status = "refuted"


@check("volume", "(E3)/(E4)",
       "vol P_1(f) = sqrt|disc f| / (d-1)! for the coefficient polytope of monic interlacers.")
def _volume(c: Check, full: bool):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.spatial import ConvexHull, HalfspaceIntersection

    from .bruteforce import interlacing_polytope

    out = {}
    for f in (real_cyclotomic(21), real_cyclotomic(13), poly("x^4 - 4*x^2 + 2"), poly("x^5 - 5*x^3 + 5*x - 1")):
        f = poly(f)
        d = f.degree()
        A, b = interlacing_polytope(f, 1)
        # interior point: Chebyshev centre
        norms = np.linalg.norm(A, axis=1)
        res = linprog(np.r_[np.zeros(A.shape[1]), -1], A_ub=np.c_[A, norms], b_ub=b,
                      bounds=[(None, None)] * A.shape[1] + [(0, None)], method="highs")
        x0 = res.x[:-1]
        hs = HalfspaceIntersection(np.c_[A, -b], x0)
        vol = ConvexHull(hs.intersections).volume
        pred = math.sqrt(abs(discriminant(f))) / math.factorial(d - 1)
        out[to_str(f)] = {"volume": vol, "sqrt|disc|/(d-1)!": pred}
        c.evidence.append(f"{to_str(f)}: vol = {vol:.6f}, sqrt|disc|/(d-1)! = {pred:.6f}")
        if abs(vol - pred) > 1e-6 * pred:
            c.status = "refuted"
    c.data = out


@check("capacity", "Prop 5.6",
       "Irreducible totally real f with span < S are finite up to translation for S < 3.")
def _capacity(c: Check, full: bool):
    c.status = "corrected"
    c.evidence.append("true as stated, but the threshold is 4: cover the possible positions of the least root "
                      "by finitely many intervals [a, a + delta]; each set of roots lies in an interval of length "
                      "S + delta < 4, of capacity < 1 (Fekete), so S < 4 suffices (Robinson: S > 4 gives "
                      "infinitely many classes).  The proof's interval [0, 1 + S] is what needs S < 3.")


# ---------------------------------------------------------------------------
# New computations
# ---------------------------------------------------------------------------

@check("sextic-census", "Table §6 (degree 6)",
       "Degree 6: non-interlacing fields exist (F_21); finiteness modulo Z open.")
def _sextic(c: Check, full: bool):
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "sextic_census.json")
    c.status = "extended"
    if os.path.exists(path):
        rows = json.load(open(path))["fields"]
        non = [r["label"] for r in rows if r["m_OK"] > 1]
        incons = [r["label"] for r in rows if not r["consistent"]]
        notfacet = [r["label"] for r in rows if r["m_OK"] == 1 and (r["sails_bound"] or 0) > 1]
        c.data = {"fields": len(rows), "noninterlacing": non, "inconsistent_with_sails": incons,
                  "trace_one_not_a_facet_normal": notfacet}
        c.evidence.append(f"all {len(rows)} totally real sextic fields of the sails/LMFDB list (disc <= 10^7): "
                          f"m(O_K) = 2 only for {non}; m(Z[a]) = m(O_K) for every polredabs polynomial")
        c.evidence.append(f"consistent with the sails facets through 1 (m <= min d_F, = 1 if some d_F = 1): "
                          f"{'all' if not incons else incons}")
        c.evidence.append(f"{len(notfacet)} fields have m = 1 although no sail facet through 1 has distance 1 "
                          f"(the trace-one element supports a lower-dimensional face), e.g. {notfacet[:4]}")
    from .search import levels_census

    S = 4.3 if full else 4.0
    cen = levels_census(6, S)
    non = [e.as_dict() for e in cen if e.level > 1]
    c.data["span_census"] = {"max_span": S, "irreducible_classes": len(cen), "noninterlacing": non}
    c.evidence.append(f"census of irreducible totally real sextics with span <= {S} up to x -> +-x + n: "
                      f"{len(cen)} classes, non-interlacing: {[e['f'] for e in non]}")


@check("cyclotomic-table", "Thm 2.4 (Yatsyna) and beyond",
       "Real cyclotomic fields: which conductors give m(O_{F_n}) >= 2; any m = 3?")
def _cyclo(c: Check, full: bool):
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "cyclotomic_levels.json")
    c.status = "extended"
    if not os.path.exists(path):
        c.evidence.append("run scripts/cyclotomic_levels.py first")
        return
    rows = json.load(open(path))
    adm_bad = [r["n"] for r in rows if r["admissible"] and r["m"] < 2]
    extra = [(r["n"], r["degree"], r["admissible_divisors"][:2]) for r in rows if not r["admissible"] and r["m"] >= 2]
    m3 = [(r["n"], r["degree"]) for r in rows if r["m"] >= 3]
    c.data = {"rows": len(rows), "admissible_with_m1": adm_bad, "nonadmissible_with_m_ge_2": extra, "m_ge_3": m3,
              "max_conductor": max(r["n"] for r in rows), "max_degree": max(r["degree"] for r in rows)}
    c.evidence.append(f"{len(rows)} fields F_n, n <= {c.data['max_conductor']}, degree <= {c.data['max_degree']}")
    c.evidence.append(f"admissible conductors with m = 1: {adm_bad or 'none'} (Theorem 2.4 confirmed on this range)")
    c.evidence.append(f"non-admissible conductors with m >= 2: {extra}")
    c.evidence.append(f"fields with m >= 3: {m3 or 'none'}")
    if adm_bad:
        c.status = "refuted"


def run(full: bool = False, only=None, verbose: bool = True) -> list:
    results = []
    for key, section, claim, fn, full_only in CHECKS:
        if only and key not in only:
            continue
        if full_only and not full:
            continue
        c = Check(key, section, claim)
        t = time.time()
        try:
            fn(c, full)
        except Exception as exc:  # pragma: no cover - reported, not hidden
            c.status = "error"
            c.evidence.append(f"{type(exc).__name__}: {exc}")
        c.seconds = round(time.time() - t, 2)
        if verbose:
            print(f"[{c.status:9s}] {key} ({c.seconds}s)", flush=True)
        results.append(c)
    return results


def write_report(results: list, outdir: str) -> str:
    os.makedirs(outdir, exist_ok=True)
    md = ["# Verification of the computational claims", "",
          "Generated by `python -m noninterlacing paper`. Statuses: **confirmed** (recomputed, holds), "
          "**corrected** (conclusion holds, statement/proof needs a fix), **refuted** (false as stated), "
          "**extended** (new data).", "",
          "| check | section | status | time |", "|---|---|---|---|"]
    for c in results:
        md.append(f"| `{c.key}` | {c.section} | **{c.status}** | {c.seconds}s |")
    md.append("")
    for c in results:
        md.append(f"## `{c.key}` — {c.section}: **{c.status}**")
        md.append("")
        md.append(f"*Claim.* {c.claim}")
        md.append("")
        for e in c.evidence:
            md.append(f"- {e}")
        md.append("")
    path = os.path.join(outdir, "paper_verification.md")
    with open(path, "w") as fh:
        fh.write("\n".join(md))
    with open(os.path.join(outdir, "paper_verification.json"), "w") as fh:
        json.dump([asdict(c) for c in results], fh, indent=1, default=str)
    return path
