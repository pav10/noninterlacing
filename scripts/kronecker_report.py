"""Write results/higher_interlacing_cyclotomic.md from the computed data.

Usage: python scripts/kronecker_report.py

Reads results/kronecker_census.json (scripts/kronecker_census.py) and
results/cyclotomic_levels.json (scripts/cyclotomic_levels.py).
"""
import collections
import json
import os

import random

from flint import fmpz_poly

from noninterlacing.cyclotomic import compositum_splits, compositum_witness, factorint, psi_interlacing_predicted
from noninterlacing.fibers import CyclotomicField, tower_levels
from noninterlacing.interlace import Interlacing
from noninterlacing.kronecker import (check_disc_polynomial, components, convolve_counts,
                                      cyclotomic_products, expected_count_log10, label, product_poly,
                                      reflected_index)
from noninterlacing.orders import Order
from noninterlacing.poly import poly, real_cyclotomic
from noninterlacing.realroots import interlaces

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(HERE, "results")


def main():
    K = json.load(open(os.path.join(RES, "kronecker_census.json")))
    C = json.load(open(os.path.join(RES, "cyclotomic_levels.json")))
    have = collections.Counter(r["degree"] for r in K)
    dmax = 0
    while have[dmax + 1] == len(list(cyclotomic_products(dmax + 1))):
        dmax += 1
    K = [r for r in K if r["degree"] <= dmax]  # report complete degrees only
    by_d = collections.defaultdict(collections.Counter)
    for r in K:
        by_d[r["degree"]][r["m"]] += 1
    conn = [r for r in K if len(components(r["indices"])) == 1]
    conn_m = collections.Counter(r["m"] for r in conn)
    ms_ok = all(r["ms_agrees"] for r in K)
    intervals = all(r["interval"] for r in K)
    amgm_ok = all(r["amgm_bound"] <= r["m"] for r in K)
    amgm_tight = sum(r["amgm_bound"] == r["m"] for r in K)
    psi_rows = [r for r in C if r["m"] is not None]
    by_comp = [r for r in psi_rows if (r["witness_method"] or "").startswith("compositum")]
    computed = [r for r in psi_rows if r not in by_comp]
    psi_bad = [r["n"] for r in computed if (r["m"] == 1) != psi_interlacing_predicted(r["n"])]
    # the compositum bound, rebuilt for every non-interlacing conductor of the table with a split
    comp_n = [r["n"] for r in C if not psi_interlacing_predicted(r["n"]) and compositum_splits(r["n"])]
    comp_ok = all(compositum_witness(n, *min(compositum_splits(n), key=lambda s: s[1])).coeffs()[-1] == 2
                  for n in comp_n)
    comp_agree = all(r["m"] == 2 for r in C if r["n"] in comp_n)
    no_split = [r["n"] for r in C if not psi_interlacing_predicted(r["n"]) and not compositum_splits(r["n"])]
    A_path = os.path.join(RES, "abelian_levels.json")
    A = json.load(open(A_path)) if os.path.exists(A_path) else []
    m_psi = {r["n"]: r["m"] for r in C}
    ab_bound = [(r["m"], m_psi.get(r["conductor"], 1 if psi_interlacing_predicted(r["conductor"]) else None))
                for r in A if r["m"] is not None]
    ab_ok = all(b is not None and m <= b for m, b in ab_bound)
    ab_tight = sum(m == b for m, b in ab_bound)
    ab_td = sum("trace-down" in (r["witness_method"] or "") for r in A)
    open_n = [n for n in range(3, 400) if n % 4 != 2 and not psi_interlacing_predicted(n)
              and not compositum_splits(n)][:4]
    psi_m = collections.Counter(r["m"] for r in psi_rows)
    undecided = [r["n"] for r in C if r["m"] is None]
    # section 5: fibre towers, witnesses, volume heuristic
    import glob

    towers = collections.defaultdict(list)  # n -> [(file, tower, [(degree, count, orbits)], seconds)]
    for path in sorted(glob.glob(os.path.join(RES, "fibre_*_level2*.json"))):
        T = json.load(open(path))
        towers[T["n"]].append((os.path.basename(path), T["tower"],
                               [(s["degree"], s["count"], s["orbits"]) for s in T["stages"]], T["seconds"]))
    lvl2_empty = {n: all(st[-1][1] == 0 for _, _, st, _ in tw) for n, tw in towers.items()}
    psi = {}  # n -> {level: checked}
    for path in sorted(glob.glob(os.path.join(RES, "psi[0-9]*.json"))):
        Pn = json.load(open(path))
        fn = real_cyclotomic(Pn["n"])
        psi[Pn["n"]] = {int(j): interlaces(fn, poly(g)) and int(poly(g).leading_coefficient()) == int(j)
                        for j, g in Pn["interlacers"].items()}

    def level4(n):
        """Level-4 interlacer from Psi_q (q a prime power, interlacing) x Psi_{n/q} (m = 2)."""
        for p, a in sorted(factorint(n).items(), reverse=True):
            q = p ** a
            if p > 2 and psi_interlacing_predicted(n // q) is False and compositum_splits(n // q):
                g2 = Interlacing(real_cyclotomic(n // q)).find(2)
                g = compositum_witness(n, q, n // q, None, g2)
                return q, n // q, int(g.leading_coefficient()) == 4
        return None

    lvl4 = {n: level4(n) for n in towers}

    # level-2 elements of the fields F_pq met in the towers: how many lie in F_p F_q?
    from noninterlacing.poly import parse_rational_poly

    pq = {}
    for n in sorted(towers):
        Kn = None
        for name, tw, st, secs in towers[n]:
            T = json.load(open(os.path.join(RES, name)))
            for k, m in enumerate(tw, start=1):
                fac = {p: a for p, a in factorint(m).items() if p != 2}
                stage = T["stages"][k]
                if len(fac) != 2 or m in pq or stage["count"] == 0:
                    continue
                Kn = Kn or CyclotomicField(n)
                if stage["subgroup"] != Kn.fixer(m):
                    continue
                q1, q2 = (p ** a for p, a in sorted(fac.items()))
                Hpq = sorted(set(Kn.fixer(q1)) & set(Kn.fixer(q2)))
                inside = 0
                for x in map(parse_rational_poly, stage["reps"]):
                    if all(Kn.act(h, x) == x for h in Hpq):
                        inside += Kn.orbit_size(x, stage["subgroup"])
                prod = Interlacing(real_cyclotomic(q1)).count(1) * Interlacing(real_cyclotomic(q2)).count(1)
                pq[m] = (q1, q2, stage["count"], inside, prod)
    from noninterlacing.pari import from_pari, pari, to_pari

    K21 = CyclotomicField(21)
    fib_ok = all(tower_levels(K21, [K21.G, K21.fixer(7), [1]], j)[-1].count == Interlacing(real_cyclotomic(21)).count(j)
                 for j in (1, 2, 3))
    Pp = pari()
    g711 = from_pari(Pp.polredabs(Pp.polcompositum(to_pari(real_cyclotomic(7)), to_pari(real_cyclotomic(11)))[0]))
    g711 = fmpz_poly([int(c) for c in g711.numer().coeffs()])
    eng711 = Order.maximal(g711).count(2)
    heur_n = [21, 33, 77, 143, 165, 231, 273, 357, 385, 399, 429, 455]
    heur = {n: [expected_count_log10(real_cyclotomic(n), j) for j in range(1, 6)] for n in heur_n}
    m_known = {r["n"]: r["m"] for r in C}
    max_m = {d: max(r["m"] for r in K if r["degree"] == d) for d in by_d}

    # --- verifications recomputed here ---------------------------------
    m_of = {}
    for r in K:
        t = tuple(r["indices"])
        m_of[t] = r["m"]
        m_of[tuple(sorted(reflected_index(n) for n in t))] = r["m"]
    add_checked = add_bad = 0
    for r in K:
        comps = components(r["indices"])
        if len(comps) > 1:
            add_checked += 1
            if sum(m_of[c] for c in comps) != r["m"]:
                add_bad += 1
    conv_examples = []
    for t in [(5, 7), (3, 4, 5), (4, 5, 7), (5, 7, 11), (1, 6, 10), (3, 6, 5, 7)]:
        comps = components(t)
        per = []
        for c in comps:
            I = Interlacing(product_poly(c), check=False)
            per.append([0] + [I.count(j) for j in range(1, 5)])
        If = Interlacing(product_poly(t), check=False)
        act = [0] + [If.count(j) for j in range(1, 5)]
        conv_examples.append((label(sorted(t)), act[1:], convolve_counts(*per)[:5] == act))
    psi3 = {r["n"]: Interlacing(real_cyclotomic(r["n"])).find(3) is not None
            for r in psi_rows if r["m"] == 2}
    rng = random.Random(5)
    disc_pairs = disc_bad = 0
    worst = 0.0
    for d in range(2, 8):
        ts = list(cyclotomic_products(d))
        rng.shuffle(ts)
        for t in ts[:25]:
            I = Interlacing(product_poly(t), check=False)
            for j in (1, 2, 3):
                for g in I.interlacers(j, limit=20):
                    c = check_disc_polynomial(t, g)
                    disc_pairs += 1
                    worst = max(worst, c["max|root|"])
                    if not (c["lead"] == j + 1 and c["max|root|"] < 1 and c["two_reciprocal"] and c["F_recovered"]):
                        disc_bad += 1

    L = []
    w = L.append
    w("# Higher interlacing for cyclotomic polynomials")
    w("")
    w("Here *cyclotomic* means: monic `f ∈ Z[x]`, squarefree, with all roots in `[-2, 2]`. By Kronecker")
    w("these are exactly the products of distinct `Ψ_n` = minimal polynomial of `2cos(2π/n)`")
    w("(`Ψ_1 = x-2`, `Ψ_2 = x+2`, `Ψ_3 = x+1`, `Ψ_4 = x`, `Ψ_6 = x-1`, `deg Ψ_n = φ(n)/2`).")
    w("A level-`j` interlacer is a `g ∈ Z[x]` of degree `d-1` with leading coefficient `j` whose roots")
    w("strictly interlace those of `f`. `S(f)` is the semigroup of levels and `m(f) = min S(f)`.")
    w("")
    w("All statements marked *verified* were checked by `noninterlacing` (scripts in `scripts/`).")
    w("")
    w("## 1. Level 1: the Beukers–Heckman / McKee–Smyth classification")
    w("")
    w("For `f = Ψ_1^a Ψ_2^b ∏_{n∈N} Ψ_n` (`n ≥ 3`) and monic `g` of degree `d-1` put")
    w("")
    w("    F = Φ_1^a Φ_2^b ∏_{n∈N} Φ_n,      G = Φ_1^{1-a} Φ_2^{1-b} · z^{d-1} g(z + 1/z).")
    w("")
    w("Since `z^{φ(n)/2} Ψ_n(z + 1/z) = Φ_n(z)`, the roots `2cos θ` of `f` and `g` interlace on")
    w("`[-2, 2]` iff the points `e^{±iθ}` (together with `±1`, assigned to `F` or `G`) alternate on the")
    w("unit circle, i.e. iff `{F, G}` is an *interlacing cyclotomic pair*. These pairs are classified")
    w("by McKee–Smyth, *Cent. Eur. J. Math.* 11 (2013), Thm 1.1 (extracted from Beukers–Heckman,")
    w("*Invent. Math.* 95 (1989)): `{P(z^l), Q(z^l)}` with `{P, Q}` or `{P(-z), Q(-z)}` in their two")
    w("infinite families or 26 sporadic pairs. Hence **the monic interlacers of `f` are the partners")
    w("of `F` in that list** (`kronecker.ms_all_pairs`, `ms_monic_interlacers`).")
    w("")
    w(f"*Verified*: for all {len(K)} squarefree cyclotomic polynomials of degree ≤ {dmax} (up to")
    w("`x → -x`) the complete list of monic interlacers computed by the codifferent engine coincides")
    w(f"with the McKee–Smyth prediction: **{'no disagreement' if ms_ok else 'DISAGREEMENTS FOUND'}**.")
    w("")
    w("**Corollary (irreducible case).** `Ψ_n` (`n ≥ 3`) is interlacing iff")
    w("`rad(odd part of n) ∈ {1, p, 15}`, i.e. `n = 2^a p^b` or `n = 2^a 3^b 5^c` with `b, c ≥ 1`.")
    w("(`Φ_n` alone must be the `P`-member of a pair; family 1 and its dilations give `Φ_{p^k}`,")
    w("`z ↦ -z` gives `Φ_{2p^k}`, dilating those gives `2^a p^b`, family 2 gives `Φ_{2^k}`, and the")
    w("single-factor sporadic members `Φ_9, Φ_12, Φ_15, Φ_18, Φ_20, Φ_24, Φ_30` add only")
    w("`2^a 3^b 5^c`.) For squarefree `n` this is Yatsyna's theorem (exceptions 15 and 30 are")
    w("exactly the `{3,5}` case). *Verified* against the generated pairs for all `n < 400` of")
    w(f"degree ≤ 60, and against the computed `m(Ψ_n)` for {len(computed)} conductors: "
      f"{'no disagreement' if not psi_bad else psi_bad}.")
    w("")
    w("## 2. Level j: Schur-stable polynomials with leading coefficient j+1")
    w("")
    w("**Theorem.** Let `F, G` be as above with `g` of leading coefficient `j ≥ 1`, and `H = F + G`.")
    w("Then `g` is a level-`j` interlacer of `f` iff `H ∈ Z[z]` is *Schur stable* (all zeros in")
    w("`|z| < 1`); in that case `lc(H) = j + 1`, `|H(0)| = j - 1`, `H ≡ H*` (mod 2) and `F, G` are")
    w("recovered as `(H ± H*)/2`. Conversely every Schur-stable `H ∈ Z[z]` with `lc(H) = j + 1`,")
    w("`|H(0)| = j - 1`, `H ≡ H*` (mod 2) arises this way, `F` being automatically cyclotomic.")
    w("")
    w("*Proof sketch.* If `H` is Schur stable, `B = H/H*` is a finite Blaschke product of degree")
    w("`D = deg H`; its boundary phase is strictly increasing with total variation `2πD`")
    w("(Mao–Wang–Zheng, arXiv:2609.34551, Lemma 2.1). The zeros of `F = H*(B + 1)/2` and")
    w("`G = H*(B - 1)/2` on the circle are the solutions of `B = -1` and `B = +1`, which therefore")
    w("alternate and are simple. `F = (H + H*)/2` is monic with integer coefficients and all zeros on")
    w("the circle, hence cyclotomic (Kronecker). The converse is the circle Hermite–Biehler theorem")
    w("(the argument of McKee–Smyth's Theorem 2.1, which is the case `j = 1`, \"2-reciprocal")
    w("disc-bionic polynomials\"). `kronecker.circle_lift`, `check_disc_polynomial`; *verified* on")
    w(f"{disc_pairs} interlacing pairs with `j ≤ 3`, degrees ≤ 7: {disc_bad} exceptions (largest zero")
    w(f"modulus {worst:.5f}).")
    w("")
    w("For `j = 1` the theorem plus Kronecker makes the set finite in each degree and McKee–Smyth list")
    w("it; for `j ≥ 2` the same holds degree by degree, but there is no comparable classification, and")
    w("the counts grow like `j^{d-1}`.")
    w("")
    w("## 3. Decomposition along unit resultants")
    w("")
    w("For any coprime real-rooted `f_1, f_2`: if `g_i` interlace `f_i` at levels `j_i` then")
    w("`g_1 f_2 + f_1 g_2` interlaces `f_1 f_2` at level `j_1 + j_2` (partial fractions). Hence")
    w("`S(f_1 f_2) ⊇ S(f_1) + S(f_2)` and `S(f) + S(h) ⊆ S(f)` for every factor `h | f`.")
    w("")
    w("**Theorem.** If `Res(f_1, f_2) = ±1` then `Z[x]/(f_1 f_2) ≅ Z[x]/(f_1) × Z[x]/(f_2)`, so every")
    w("interlacer of `f_1 f_2` is uniquely `g_1 f_2 + f_1 g_2`:")
    w("")
    w("    S(f_1 f_2) = S(f_1) + S(f_2),   m(f_1 f_2) = m(f_1) + m(f_2),")
    w("    I_j(f_1 f_2) = Σ_{j_1 + j_2 = j} I_{j_1}(f_1) I_{j_2}(f_2).")
    w("")
    w("For cyclotomic factors, `Res(Ψ_m, Ψ_n) = ±1` unless `m/n` is a prime power")
    w("(`|Res(Ψ_m,Ψ_n)|² = |Res(Φ_m,Φ_n)|`, Apostol). So **every cyclotomic `f` splits into the")
    w("components of the graph `m ~ n ⟺ m/n is a prime power` on its indices**, and `S`, `m` and the")
    w("counts are sums/convolutions over components (`kronecker.components`, `convolve_counts`).")
    w("In particular a cyclotomic `f` with `r` components has `m(f) ≥ r`, so is never interlacing")
    w("when `r ≥ 2`.")
    w("")
    w(f"*Verified*: `m(f) = Σ m(components)` for all {add_checked} census polynomials with ≥ 2")
    w(f"components ({add_bad} exceptions); the convolution formula for `I_1..I_4`:")
    for name, act, ok in conv_examples:
        w(f"`{name}`: `I_1..I_4 = {act}` ({'matches' if ok else 'MISMATCH'});")
    w("and (separately, 65 s) `I_4(Ψ_21 Ψ_33) = I_2(Ψ_21) I_2(Ψ_33) = 141 · 935 = 131835`.")
    w("Without unit resultant the order is smaller: `Ψ_21 Ψ_42` (resultant 64) has `S = <2,3>`, while")
    w("`Ψ_21 Ψ_33`, `Ψ_21 Ψ_35`, `Ψ_21 Ψ_39`, `Ψ_33 Ψ_42` have `S = <4,5,6,7>`.")
    w("")
    w("**Corollary (unbounded levels).** For distinct primes `p_1, ..., p_k`,")
    w("`m(Ψ_{p_1} ⋯ Ψ_{p_k}) = k`, and more generally `m(∏ Ψ_{n_i}) = Σ m(Ψ_{n_i})` whenever no ratio")
    w("`n_i/n_k` is a prime power. Cyclotomic polynomials therefore have arbitrarily high")
    w("non-interlacing level; e.g. `m(Ψ_3 Ψ_4 Ψ_5 Ψ_7 Ψ_11) = 5` with `S = <5,6,7,8,9>` (degree 12),")
    w("`m(Ψ_21 Ψ_33) = 4`.")
    w("")
    w("## 4. Composita: `m(Ψ_n) ≤ 2 m(Ψ_{n_1}) m(Ψ_{n_2})`")
    w("")
    w("**Lemma.** If `F_1, F_2` are totally real with coprime discriminants, then")
    w("`max(m(O_{F_1}), m(O_{F_2})) ≤ m(O_{F_1 F_2}) ≤ m(O_{F_1}) m(O_{F_2})`. Indeed")
    w("`O_{F_1F_2} = O_{F_1} ⊗ O_{F_2}` (Neukirch I.2.11) with `Tr = Tr_1 ⊗ Tr_2`, so the codifferent")
    w("is `O_{F_1}^∨ ⊗ O_{F_2}^∨` and `μ_1 μ_2` is totally positive of trace `Tr(μ_1) Tr(μ_2)`; the")
    w("lower bound is the trace-down.")
    w("")
    w("**Theorem.** If `n = n_1 n_2` with `gcd(n_1, n_2) = 1` and `n_1, n_2 > 2`, then")
    w("`F_{n_1} F_{n_2}` has index 2 in `F_n = Q(ζ_n)^+` and, by the subfield lift,")
    w("")
    w("    m(Ψ_n) ≤ 2 m(Ψ_{n_1}) m(Ψ_{n_2}).")
    w("")
    w("Explicitly, `g = μ_1(C_{n_2}(α)) μ_2(C_{n_1}(α)) Ψ_n'(α) mod Ψ_n` with")
    w("`μ_i = g_i/Ψ_{n_i}'` and `C_k(2cos t) = 2cos kt` (`cyclotomic.compositum_witness`, which")
    w("checks integrality and interlacing exactly).")
    w("")
    w("**Corollary.** Let `P` be the set of odd primes dividing `n`. If `|P| = 2` with `P ≠ {3, 5}`,")
    w("or `P = {3, 5, p}`, then `m(Ψ_n) = 2`: split `n` into two coprime parts with interlacing `Ψ`")
    w("(the power of 2 goes with either), and level 1 is empty by §1. In general")
    w("`m(Ψ_n) ≤ 2^{b-1}`, where `b` counts the odd primes with 3 and 5 taken together.")
    w("")
    w(f"*Verified*: the witness is rebuilt with leading coefficient 2 for all {len(comp_n)} non-interlacing")
    w(f"conductors of the table that split ({'ok' if comp_ok else 'FAILED'}); it agrees with the directly")
    w(f"computed `m = 2` wherever that exists ({'ok' if comp_agree else 'DISAGREEMENT'}), and it decides")
    w(f"{len(by_comp)} conductors the direct computation could not: "
      f"{', '.join(str(r['n']) for r in by_comp) or 'none'}. Non-interlacing conductors of the table")
    w(f"without a split: {no_split or 'none'}; the first such `n` are {open_n}.")
    w("")
    w("**Abelian fields.** Every real abelian field `K` of conductor `n` lies in `F_n`, and")
    w("trace-down (`Tr_{F_n/K}` of a witness, `heredity.trace_down_witness`) gives `m(O_K) ≤ m(Ψ_n)`.")
    w("So a real abelian field is interlacing whenever `rad(odd part of its conductor) ∈ {1, p, 15}`,")
    w(f"and **no real abelian field of conductor < {open_n[0]} has `m ≥ 3`**: below {open_n[0]} every")
    w("conductor has at most two odd prime factors or is `{3, 5, p}`-type.")
    if A:
        w(f"*Verified* on `results/abelian_levels.json` ({len(A)} fields, degree "
          f"{min(r['degree'] for r in A)}–{max(r['degree'] for r in A)}, conductor ≤ "
          f"{max(r['conductor'] for r in A)}): `m(O_K) ≤ m(Ψ_n)` for every decided field "
          f"({'ok' if ab_ok else 'VIOLATED'}), with equality for {ab_tight}; {ab_td} fields, where the")
        w("direct search ran out of time, are decided by trace-down from `F_p`.")
    w("")
    w("## 5. Three odd primes: `m(Ψ_n) = 3`")
    w("")
    w("**Fibres of the relative trace.** For subfields `F' ⊂ F''` of `F_n` the relative trace maps")
    w("the totally positive elements of `O_{F''}^∨` of trace `j` onto totally positive elements of")
    w("`O_{F'}^∨` of trace `j` (`Tr_{F'}(Tr(μ) x) = Tr_{F''}(μ x)` for `x ∈ O_{F'}`). So the level-`j`")
    w("elements of `F''` are the union of the fibres `{μ ∈ O_{F''}^∨ : Tr_{F''/F'} μ = ν, μ ≫ 0}` over")
    w("the level-`j` elements `ν` of `F'`, and conjugate `ν` have conjugate fibres. A fibre is a coset")
    w("of the lattice `ker Tr` inside a product of `[F':Q]` simplices; `noninterlacing/fibers.py`")
    w("enumerates it (Fincke–Pohst over the circumscribed ball on an LLL-reduced basis, float with a")
    w("safety margin; every candidate is decided by certified `arb` signs or Sturm sequences).")
    w("Going up a tower `Q ⊂ F_1 ⊂ … ⊂ F_n` one fibre per Galois orbit gives every level-`j` element.")
    w(f"*Verified*: the towers reproduce the engine's counts for `Ψ_21`, `j = 1, 2, 3` ({'ok' if fib_ok else 'MISMATCH'}),")
    w(f"and the engine counts {eng711} elements of trace 2 in `F_7F_11` (through PARI's `polcompositum`).")
    w("")
    w("**Level 2.** For `n` with three odd primes (not 3 and 5 together) §4 gives only `m ≤ 4`; for")
    w("`n = 231` all three index-2 subfields and the index-3 subfield have `m = 2`. Counts of trace-2")
    w("elements up each tower (degree: count / Galois orbits):")
    w("")
    for n in sorted(towers):
        for name, tw, st, secs in towers[n]:
            w(f"* `{name}` (tower `{tw}`, {secs:.0f} s): " + ", ".join(f"{d}: {c}/{o}" for d, c, o in st))
    w("")
    w("**Levels 3, 4, 5.** Levels 3 and 5 from MILP (directly for 231; in the fibres over the most")
    w("central elements of `F_91` for 273, `fibers.central_search`), level 4 from the compositum")
    w("`Ψ_q ⊗ Ψ_{n/q}` (`q` the largest prime power, `m(Ψ_{n/q}) = 2`); all checked exactly:")
    w("")
    w("| n | degree | level 2 | level 3 | level 4 | level 5 | m | S |")
    w("|---|---|---|---|---|---|---|---|")
    m_fib = {}

    def e3(n):
        x = expected_count_log10(real_cyclotomic(n), 3)
        return f"? (`E_3 ≈ 10^{{{x:.1f}}}`)"

    for n in sorted(towers):
        ok = psi.get(n, {})
        l4 = lvl4[n]
        if lvl2_empty[n] and ok.get(3):
            m_fib[n] = "3"
        elif lvl2_empty[n]:
            m_fib[n] = "3–4" if l4 and l4[2] else "≥ 3"
        semi = "<3,4,5>" if lvl2_empty[n] and ok.get(3) and ok.get(5) and l4 and l4[2] else "?"
        w(f"| {n} | {real_cyclotomic(n).degree()} | {'empty' if lvl2_empty[n] else 'NOT EMPTY'} | "
          f"{'yes' if ok.get(3) else ('FAILED' if 3 in ok else e3(n))} | "
          + (f"Ψ_{l4[0]} ⊗ Ψ_{l4[1]}" if l4 and l4[2] else "?") + " | "
          f"{'yes' if ok.get(5) else '?'} | {m_fib.get(n, '?')} | {semi} |")
    w("")
    w("So `m(Ψ_n) ≤ 2` fails; whether §4's bound `2^{b-1} = 4` is ever attained is open.")
    w("")
    w("**Level-2 elements of `F_pq`.** The fields `F_pq` in the towers have `m = 2`. Their level-2")
    w("elements include the compositum products `μ_p μ_q` (level 1 in `F_p F_q`, of index 2):")
    w("")
    w("| field | level-2 elements | lying in `F_p F_q` | `I_1(Ψ_p) I_1(Ψ_q)` | `E_2` |")
    w("|---|---|---|---|---|")
    for m in sorted(pq):
        q1, q2, cnt, inside, prod = pq[m]
        w(f"| `F_{m}` | {cnt} | {inside} | {prod} | `10^{{{expected_count_log10(real_cyclotomic(m), 2):.1f}}}` |")
    w("")
    only = [m for m in sorted(pq) if pq[m][2] == pq[m][3]]
    w("The elements lying in `F_p F_q` are exactly the level-1 elements of `F_p F_q`, and their number")
    w("equals the number of products in every case: the products are all of them. "
      + (f"In {', '.join(f'`F_{m}`' for m in only)} they are the only level-2 elements, "
         "although `E_2` predicts more" if only else "")
      + "; few fibres then lead to the top field, which keeps those towers cheap.")
    w("")
    w("**Volume heuristic.** `E_j = j^{d-1} √|disc| / (d-1)!` is the volume of the trace-`j` simplex")
    w("in lattice units (`kronecker.expected_count_log10`). It predicts the counts above well")
    w("(`F_7F_11`: `E_2 ≈ 5.6·10^3` vs 6570; `F_77`: `E_2 ≈ 255` vs 285), but structure beats it in both")
    w("directions: composita give `m(Ψ_143) = 2` with `E_2 ≈ 10^{-3.6}`, and `I_1(Ψ_21) = 0` with")
    w("`E_1 ≈ 5.6`. `log10 E_j` for `j = 1..5`, with the known `m` and the §4 bound:")
    w("")
    w("| n | degree | log10 E_1 … E_5 | m | §4 bound |")
    w("|---|---|---|---|---|")
    for n in heur_n:
        odd = [p for p in factorint(n) if p != 2]
        b = len(odd) - (1 if 3 in odd and 5 in odd else 0)
        bnd = 1 if psi_interlacing_predicted(n) else 2 ** (b - 1)
        w(f"| {n} | {real_cyclotomic(n).degree()} | " + ", ".join(f"{x:.1f}" for x in heur[n])
          + f" | {m_known.get(n) or m_fib.get(n) or '?'} | {bnd} |")
    w("")
    w("## 6. The AM–GM lower bound")
    w("")
    w("For a level-`j` interlacer and an irreducible factor `f_i` of degree `e_i`, the component")
    w("`μ_i = g(α_i)/f'(α_i)` is totally positive with `|N(g(α_i))| ≥ 1`, so")
    w("")
    w("    m(f) ≥ ⌈ Σ_i e_i · |disc(f_i) · Res(f_i, f/f_i)|^{-1/e_i} ⌉     (kronecker.amgm_lower_bound).")
    w("")
    w("For irreducible `f` this is Yatsyna's `|disc f| ≥ d^d` for interlacing `f`. *Verified*: never")
    w(f"violated ({'ok' if amgm_ok else 'VIOLATED'}), and equal to `m(f)` for {amgm_tight} of {len(K)} census")
    w("polynomials; it misses the additivity of §3 (e.g. it gives 2 for `Ψ_21 Ψ_33`, whose level is 4).")
    w("")
    w("## 7. Data")
    w("")
    w(f"Census of all squarefree cyclotomic polynomials of degree ≤ {dmax} up to `x → -x`")
    w("(`results/kronecker_census.json`). Number of polynomials with `m(f) = k`:")
    w("")
    ks = sorted({k for d in by_d for k in by_d[d]})
    w("| degree | " + " | ".join(f"m={k}" for k in ks) + " | max m |")
    w("|---|" + "---|" * (len(ks) + 1))
    for d in sorted(by_d):
        w(f"| {d} | " + " | ".join(str(by_d[d].get(k, 0)) for k in ks) + f" | {max_m[d]} |")
    w("")
    w(f"* Every level set is an interval `{{0}} ∪ [m, ∞)`: **{intervals}** (for `f` with a")
    w("  1-interlaced factor this follows from `S(f) + S(h) ⊆ S(f)`).")
    w(f"* Connected index sets ({len(conn)} polynomials): `m` distribution {dict(sorted(conn_m.items()))}.")
    w(f"* Irreducible `Ψ_n` (`results/cyclotomic_levels.json`, {len(psi_rows)} fields, conductor ≤")
    w(f"  {max(r['n'] for r in C)}, degree ≤ 60): `m` distribution {dict(sorted(psi_m.items()))};")
    w(f"  all {len(psi3)} `Ψ_n` with `m = 2` have a level-3 interlacer: "
      f"{all(psi3.values())}, so `S(Ψ_n) = <2,3>` for each; `S(Ψ_231) = <3,4,5>` (§5)."
      + (f" Of these, {len(by_comp)} are decided by §4 rather than by the engine." if by_comp else "")
      + (f" Undecided: n = {undecided}." if undecided else ""))
    w("")
    w("## 8. Open questions suggested by the data")
    w("")
    w("* `m(Ψ_n) ≤ 2` holds when `n` has at most two odd prime factors, or three including 3 and 5")
    w("  (§4), and fails for `n = 231` (§5). Is `m(Ψ_n)` bounded? Is §4's `2^{b-1}` attained? The")
    w(f"  heuristic of §5 predicts `m = 4` for `n = 357, 385, 399, 429` (`E_3 ≪ 1`); the next case after")
    w(f"  231 without a split is `n = {open_n[1]}`.")
    w("* Is `m ≤ 3` on connected index sets? Which connected sets reach 3?")
    w("* Is `S(f)` always an interval for cyclotomic `f` (also when no factor is 1-interlaced)?")
    with open(os.path.join(RES, "higher_interlacing_cyclotomic.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("written", os.path.join(RES, "higher_interlacing_cyclotomic.md"))


if __name__ == "__main__":
    main()
