"""Write results/higher_interlacing_cyclotomic.md from the computed data.

Usage: python scripts/kronecker_report.py

Reads results/kronecker_census.json (scripts/kronecker_census.py) and
results/cyclotomic_levels.json (scripts/cyclotomic_levels.py).
"""
import collections
import json
import os

import random

from noninterlacing.interlace import Interlacing
from noninterlacing.kronecker import (check_disc_polynomial, components, convolve_counts,
                                      cyclotomic_products, label, product_poly, reflected_index)
from noninterlacing.poly import real_cyclotomic

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(HERE, "results")


def predicted_psi_interlacing(n):
    from noninterlacing.cyclotomic import factorint

    odd = [p for p in factorint(n) if p != 2]
    r = 1
    for p in odd:
        r *= p
    return r == 1 or r == 15 or len(odd) == 1


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
    psi_bad = [r["n"] for r in psi_rows if (r["m"] == 1) != predicted_psi_interlacing(r["n"])]
    psi_m = collections.Counter(r["m"] for r in psi_rows)
    undecided = [r["n"] for r in C if r["m"] is None]
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
    w(f"degree ≤ 60, and against the computed `m(Ψ_n)` for {len(psi_rows)} conductors: "
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
    w("## 4. The AM–GM lower bound")
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
    w("## 5. Data")
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
    w(f"  all {len(psi3)} non-interlacing `Ψ_n` have a level-3 interlacer: "
      f"{all(psi3.values())}, so `S(Ψ_n) = <2,3>` for each."
      + (f" Undecided (time budget): n = {undecided}." if undecided else ""))
    w("")
    w("## 6. Open questions suggested by the data")
    w("")
    w("* Is `m(Ψ_n) ≤ 2` for every `n` (equivalently `m(O_{F_n}) ≤ 2`)? The subfield lift")
    w("  `m(F_n) ≤ [F_n : F] m(F)` proves it whenever `F_n` has an interlacing subfield of index 2")
    w("  (e.g. `3 ‖ n` or `4 | n` with `Ψ_{n/3}` resp. `Ψ_{n/2}` interlacing).")
    w("* Is `m ≤ 3` on connected index sets? Which connected sets reach 3?")
    w("* Is `S(f)` always an interval for cyclotomic `f` (also when no factor is 1-interlaced)?")
    with open(os.path.join(RES, "higher_interlacing_cyclotomic.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("written", os.path.join(RES, "higher_interlacing_cyclotomic.md"))


if __name__ == "__main__":
    main()
