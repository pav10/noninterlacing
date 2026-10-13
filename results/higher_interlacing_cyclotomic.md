# Higher interlacing for cyclotomic polynomials

Here *cyclotomic* means: monic `f ∈ Z[x]`, squarefree, with all roots in `[-2, 2]`. By Kronecker
these are exactly the products of distinct `Ψ_n` = minimal polynomial of `2cos(2π/n)`
(`Ψ_1 = x-2`, `Ψ_2 = x+2`, `Ψ_3 = x+1`, `Ψ_4 = x`, `Ψ_6 = x-1`, `deg Ψ_n = φ(n)/2`).
A level-`j` interlacer is a `g ∈ Z[x]` of degree `d-1` with leading coefficient `j` whose roots
strictly interlace those of `f`. `S(f)` is the semigroup of levels and `m(f) = min S(f)`.

All statements marked *verified* were checked by `noninterlacing` (scripts in `scripts/`).

## 1. Level 1: the Beukers–Heckman / McKee–Smyth classification

For `f = Ψ_1^a Ψ_2^b ∏_{n∈N} Ψ_n` (`n ≥ 3`) and monic `g` of degree `d-1` put

    F = Φ_1^a Φ_2^b ∏_{n∈N} Φ_n,      G = Φ_1^{1-a} Φ_2^{1-b} · z^{d-1} g(z + 1/z).

Since `z^{φ(n)/2} Ψ_n(z + 1/z) = Φ_n(z)`, the roots `2cos θ` of `f` and `g` interlace on
`[-2, 2]` iff the points `e^{±iθ}` (together with `±1`, assigned to `F` or `G`) alternate on the
unit circle, i.e. iff `{F, G}` is an *interlacing cyclotomic pair*. These pairs are classified
by McKee–Smyth, *Cent. Eur. J. Math.* 11 (2013), Thm 1.1 (extracted from Beukers–Heckman,
*Invent. Math.* 95 (1989)): `{P(z^l), Q(z^l)}` with `{P, Q}` or `{P(-z), Q(-z)}` in their two
infinite families or 26 sporadic pairs. Hence **the monic interlacers of `f` are the partners
of `F` in that list** (`kronecker.ms_all_pairs`, `ms_monic_interlacers`).

*Verified*: for all 5678 squarefree cyclotomic polynomials of degree ≤ 11 (up to
`x → -x`) the complete list of monic interlacers computed by the codifferent engine coincides
with the McKee–Smyth prediction: **no disagreement**.

**Corollary (irreducible case).** `Ψ_n` (`n ≥ 3`) is interlacing iff
`rad(odd part of n) ∈ {1, p, 15}`, i.e. `n = 2^a p^b` or `n = 2^a 3^b 5^c` with `b, c ≥ 1`.
(`Φ_n` alone must be the `P`-member of a pair; family 1 and its dilations give `Φ_{p^k}`,
`z ↦ -z` gives `Φ_{2p^k}`, dilating those gives `2^a p^b`, family 2 gives `Φ_{2^k}`, and the
single-factor sporadic members `Φ_9, Φ_12, Φ_15, Φ_18, Φ_20, Φ_24, Φ_30` add only
`2^a 3^b 5^c`.) For squarefree `n` this is Yatsyna's theorem (exceptions 15 and 30 are
exactly the `{3,5}` case). *Verified* against the generated pairs for all `n < 400` of
degree ≤ 60, and against the computed `m(Ψ_n)` for 146 conductors: no disagreement.

## 2. Level j: Schur-stable polynomials with leading coefficient j+1

**Theorem.** Let `F, G` be as above with `g` of leading coefficient `j ≥ 1`, and `H = F + G`.
Then `g` is a level-`j` interlacer of `f` iff `H ∈ Z[z]` is *Schur stable* (all zeros in
`|z| < 1`); in that case `lc(H) = j + 1`, `|H(0)| = j - 1`, `H ≡ H*` (mod 2) and `F, G` are
recovered as `(H ± H*)/2`. Conversely every Schur-stable `H ∈ Z[z]` with `lc(H) = j + 1`,
`|H(0)| = j - 1`, `H ≡ H*` (mod 2) arises this way, `F` being automatically cyclotomic.

*Proof sketch.* If `H` is Schur stable, `B = H/H*` is a finite Blaschke product of degree
`D = deg H`; its boundary phase is strictly increasing with total variation `2πD`
(Mao–Wang–Zheng, arXiv:2609.34551, Lemma 2.1). The zeros of `F = H*(B + 1)/2` and
`G = H*(B - 1)/2` on the circle are the solutions of `B = -1` and `B = +1`, which therefore
alternate and are simple. `F = (H + H*)/2` is monic with integer coefficients and all zeros on
the circle, hence cyclotomic (Kronecker). The converse is the circle Hermite–Biehler theorem
(the argument of McKee–Smyth's Theorem 2.1, which is the case `j = 1`, "2-reciprocal
disc-bionic polynomials"). `kronecker.circle_lift`, `check_disc_polynomial`; *verified* on
3357 interlacing pairs with `j ≤ 3`, degrees ≤ 7: 0 exceptions (largest zero
modulus 0.99996).

For `j = 1` the theorem plus Kronecker makes the set finite in each degree and McKee–Smyth list
it; for `j ≥ 2` the same holds degree by degree, but there is no comparable classification, and
the counts grow like `j^{d-1}`.

## 3. Decomposition along unit resultants

For any coprime real-rooted `f_1, f_2`: if `g_i` interlace `f_i` at levels `j_i` then
`g_1 f_2 + f_1 g_2` interlaces `f_1 f_2` at level `j_1 + j_2` (partial fractions). Hence
`S(f_1 f_2) ⊇ S(f_1) + S(f_2)` and `S(f) + S(h) ⊆ S(f)` for every factor `h | f`.

**Theorem.** If `Res(f_1, f_2) = ±1` then `Z[x]/(f_1 f_2) ≅ Z[x]/(f_1) × Z[x]/(f_2)`, so every
interlacer of `f_1 f_2` is uniquely `g_1 f_2 + f_1 g_2`:

    S(f_1 f_2) = S(f_1) + S(f_2),   m(f_1 f_2) = m(f_1) + m(f_2),
    I_j(f_1 f_2) = Σ_{j_1 + j_2 = j} I_{j_1}(f_1) I_{j_2}(f_2).

For cyclotomic factors, `Res(Ψ_m, Ψ_n) = ±1` unless `m/n` is a prime power
(`|Res(Ψ_m,Ψ_n)|² = |Res(Φ_m,Φ_n)|`, Apostol). So **every cyclotomic `f` splits into the
components of the graph `m ~ n ⟺ m/n is a prime power` on its indices**, and `S`, `m` and the
counts are sums/convolutions over components (`kronecker.components`, `convolve_counts`).
In particular a cyclotomic `f` with `r` components has `m(f) ≥ r`, so is never interlacing
when `r ≥ 2`.

*Verified*: `m(f) = Σ m(components)` for all 4047 census polynomials with ≥ 2
components (0 exceptions); the convolution formula for `I_1..I_4`:
`Psi_5*Psi_7`: `I_1..I_4 = [0, 6, 45, 155]` (matches);
`Psi_3*Psi_4*Psi_5`: `I_1..I_4 = [0, 0, 2, 9]` (matches);
`Psi_4*Psi_5*Psi_7`: `I_1..I_4 = [0, 0, 6, 51]` (matches);
`Psi_5*Psi_7*Psi_11`: `I_1..I_4 = [0, 0, 30, 705]` (matches);
`Psi_1*Psi_6*Psi_10`: `I_1..I_4 = [0, 0, 2, 9]` (matches);
`Psi_3*Psi_5*Psi_6*Psi_7`: `I_1..I_4 = [0, 0, 6, 63]` (matches);
and (separately, 65 s) `I_4(Ψ_21 Ψ_33) = I_2(Ψ_21) I_2(Ψ_33) = 141 · 935 = 131835`.
Without unit resultant the order is smaller: `Ψ_21 Ψ_42` (resultant 64) has `S = <2,3>`, while
`Ψ_21 Ψ_33`, `Ψ_21 Ψ_35`, `Ψ_21 Ψ_39`, `Ψ_33 Ψ_42` have `S = <4,5,6,7>`.

**Corollary (unbounded levels).** For distinct primes `p_1, ..., p_k`,
`m(Ψ_{p_1} ⋯ Ψ_{p_k}) = k`, and more generally `m(∏ Ψ_{n_i}) = Σ m(Ψ_{n_i})` whenever no ratio
`n_i/n_k` is a prime power. Cyclotomic polynomials therefore have arbitrarily high
non-interlacing level; e.g. `m(Ψ_3 Ψ_4 Ψ_5 Ψ_7 Ψ_11) = 5` with `S = <5,6,7,8,9>` (degree 12),
`m(Ψ_21 Ψ_33) = 4`.

## 4. Composita: `m(Ψ_n) ≤ 2 m(Ψ_{n_1}) m(Ψ_{n_2})`

**Lemma.** If `F_1, F_2` are totally real with coprime discriminants, then
`max(m(O_{F_1}), m(O_{F_2})) ≤ m(O_{F_1 F_2}) ≤ m(O_{F_1}) m(O_{F_2})`. Indeed
`O_{F_1F_2} = O_{F_1} ⊗ O_{F_2}` (Neukirch I.2.11) with `Tr = Tr_1 ⊗ Tr_2`, so the codifferent
is `O_{F_1}^∨ ⊗ O_{F_2}^∨` and `μ_1 μ_2` is totally positive of trace `Tr(μ_1) Tr(μ_2)`; the
lower bound is the trace-down.

**Theorem.** If `n = n_1 n_2` with `gcd(n_1, n_2) = 1` and `n_1, n_2 > 2`, then
`F_{n_1} F_{n_2}` has index 2 in `F_n = Q(ζ_n)^+` and, by the subfield lift,

    m(Ψ_n) ≤ 2 m(Ψ_{n_1}) m(Ψ_{n_2}).

Explicitly, `g = μ_1(C_{n_2}(α)) μ_2(C_{n_1}(α)) Ψ_n'(α) mod Ψ_n` with
`μ_i = g_i/Ψ_{n_i}'` and `C_k(2cos t) = 2cos kt` (`cyclotomic.compositum_witness`, which
checks integrality and interlacing exactly).

**Corollary.** Let `P` be the set of odd primes dividing `n`. If `|P| = 2` with `P ≠ {3, 5}`,
or `P = {3, 5, p}`, then `m(Ψ_n) = 2`: split `n` into two coprime parts with interlacing `Ψ`
(the power of 2 goes with either), and level 1 is empty by §1. In general
`m(Ψ_n) ≤ 2^{b-1}`, where `b` counts the odd primes with 3 and 5 taken together.

*Verified*: the witness is rebuilt with leading coefficient 2 for all 58 non-interlacing
conductors of the table that split (ok); it agrees with the directly
computed `m = 2` wherever that exists (ok), and it decides
9 conductors the direct computation could not: 115, 119, 133, 143, 145, 155, 165, 195, 396. Non-interlacing conductors of the table
without a split: [231]; the first such `n` are [231, 273, 357, 385].

## 5. The AM–GM lower bound

For a level-`j` interlacer and an irreducible factor `f_i` of degree `e_i`, the component
`μ_i = g(α_i)/f'(α_i)` is totally positive with `|N(g(α_i))| ≥ 1`, so

    m(f) ≥ ⌈ Σ_i e_i · |disc(f_i) · Res(f_i, f/f_i)|^{-1/e_i} ⌉     (kronecker.amgm_lower_bound).

For irreducible `f` this is Yatsyna's `|disc f| ≥ d^d` for interlacing `f`. *Verified*: never
violated (ok), and equal to `m(f)` for 4209 of 5678 census
polynomials; it misses the additivity of §3 (e.g. it gives 2 for `Ψ_21 Ψ_33`, whose level is 4).

## 6. Data

Census of all squarefree cyclotomic polynomials of degree ≤ 11 up to `x → -x`
(`results/kronecker_census.json`). Number of polynomials with `m(f) = k`:

| degree | m=1 | m=2 | m=3 | m=4 | m=5 | max m |
|---|---|---|---|---|---|---|
| 1 | 3 | 0 | 0 | 0 | 0 | 1 |
| 2 | 7 | 2 | 0 | 0 | 0 | 2 |
| 3 | 11 | 8 | 0 | 0 | 0 | 2 |
| 4 | 19 | 21 | 3 | 0 | 0 | 3 |
| 5 | 23 | 50 | 9 | 0 | 0 | 3 |
| 6 | 27 | 97 | 36 | 0 | 0 | 3 |
| 7 | 33 | 160 | 85 | 7 | 0 | 4 |
| 8 | 36 | 260 | 187 | 24 | 0 | 4 |
| 9 | 50 | 373 | 350 | 76 | 0 | 4 |
| 10 | 43 | 569 | 624 | 184 | 7 | 5 |
| 11 | 53 | 769 | 1050 | 392 | 30 | 5 |

* Every level set is an interval `{0} ∪ [m, ∞)`: **True** (for `f` with a
  1-interlaced factor this follows from `S(f) + S(h) ⊆ S(f)`).
* Connected index sets (1631 polynomials): `m` distribution {1: 305, 2: 1200, 3: 126}.
* Irreducible `Ψ_n` (`results/cyclotomic_levels.json`, 155 fields, conductor ≤
  396, degree ≤ 60): `m` distribution {1: 97, 2: 58};
  all 58 non-interlacing `Ψ_n` have a level-3 interlacer: True, so `S(Ψ_n) = <2,3>` for each. Of these, 9 are decided by §4 rather than by the engine. Undecided: n = [231].

## 7. Open questions suggested by the data

* Is `m(Ψ_n) ≤ 2` for every `n` (equivalently `m(O_{F_n}) ≤ 2`)? §4 proves it when `n` has
  at most two odd prime factors, or exactly three including 3 and 5; open from `n = 231`
  (degree 60, `2 ≤ m ≤ 4`) on, where §4 gives only `m ≤ 2^{b-1}`.
* Is `m ≤ 3` on connected index sets? Which connected sets reach 3?
* Is `S(f)` always an interval for cyclotomic `f` (also when no factor is 1-interlaced)?
