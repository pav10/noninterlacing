# noninterlacing — interlacing polynomials of every level

A computational companion to the draft *"Infinitely many irreducible
non-interlacing polynomials of fixed degree, with unbounded span and Mahler
measure"*: it computes interlacing polynomials of every level efficiently and
exactly, finds and certifies non-interlacing polynomials, works with
codifferents of arbitrary totally real orders, and re-runs every computational
claim of the draft ([`results/paper_verification.md`](results/paper_verification.md)).

## Definitions

Let `f ∈ Z[x]` be monic, squarefree, of degree `d`, with real roots
`r_1 > r_2 > … > r_d`. For `j ≥ 1` a **level-j interlacer** is a `g ∈ Z[x]` of
degree `d − 1` with leading coefficient `j` whose roots satisfy
`r_i > s_i > r_{i+1}`. `f` is **j-interlaced** if one exists; *interlacing* means
1-interlaced.

With `α` a root of `f`, `Z[α]^∨ = f'(α)^{-1} Z[α]` (Euler), and

    g  ↦  μ = g(α)/f'(α)

is a bijection between level-`j` interlacers and totally positive elements of
the codifferent `Z[α]^∨` of trace `j` (the trace is the `x^{d-1}` coefficient
of `g`; positivity of every `g(r_i)/f'(r_i)` is exactly strict interlacing).
Consequences that the program uses and reports:

* **levels form a numerical semigroup** `S(f) = {j : f is j-interlaced}`: it is
  closed under addition, contains `d` (`g = f'`) and `1 + dN` for large `N`.
  Its gaps are the levels at which `f` is *not* interlaced; there are finitely
  many;
* `m(f) := min S(f) = m(Z[α])` (Yatsyna's minimal trace); `f` is
  non-interlacing iff `m(f) ≥ 2`;
* **Siegel ceiling**: `m(f) ≤ m(O_K) ≤ dim M_{2d}(SL_2(Z))`, so `m(f) = 1` in
  degrees 2–5 and 7, `m(f) ≤ 2` for `d ≤ 11` and `d = 13`, and level 3 needs
  `d ≥ 12`;
* the **margin** `ε_j(f) = max_g min_i g(r_i)/f'(r_i)` (over `g` with
  `x^{d-1}`-coefficient `j`) is positive iff `f` is `j`-interlaced. It measures
  how far from (non-)interlacing `f` is.

## Install

```bash
pip install -e .            # numpy, scipy, python-flint
pip install cypari2         # optional: maximal orders, subfields, zeta values (PARI/GP)
python -m pytest tests      # 51 tests, ~1.5 min
```

## Command line

```bash
python -m noninterlacing level "x^6 - x^5 - 6*x^4 + 6*x^3 + 8*x^2 - 8*x + 1"
#   min level m(f): 2, level semigroup <2, 3>, gaps [1], margins, a level-2 interlacer
python -m noninterlacing interlacers "x^6 - 7*x^4 + 14*x^2 - 7"          # all monic interlacers
python -m noninterlacing interlacers "<f>" --level 3 --count             # I_3(f)
python -m noninterlacing check "x^6 - 7*x^4 + 14*x^2 - 7" "x^5 - 5*x^3 + 5*x"   # exact test
python -m noninterlacing order "<f>" --maximal --trace 2                 # m(O_K), elements of O_K^∨
python -m noninterlacing heredity "<f>"                                  # m(Z[α]) ≥ m(Z[α] ∩ F)
python -m noninterlacing cyclotomic 20..100                              # levels of Q(ζ_n)^+
python -m noninterlacing family --N 2..30                                # Theorem 3.1 for f_N
python -m noninterlacing siegel --degrees 2..20                          # b_k(2d), sign tests
python -m noninterlacing siegel --zeta "<f>"                             # Siegel's formula vs PARI
python -m noninterlacing census -d 6 -S 4.3                              # all sextics of span ≤ 4.3
python -m noninterlacing random -d 9 -n 500                              # random odd-degree search
python -m noninterlacing paper [--full]                                  # verify the draft
```
`--json` (before the subcommand) gives machine-readable output.

## Python

```python
from noninterlacing import Interlacing, Order, real_cyclotomic, to_str

I = Interlacing(real_cyclotomic(21))      # beta = 2cos(2pi/21)
I.count(1), I.count(2), I.count(3)        # (0, 141, 1324)
I.levels()                                # <2, 3>   (gaps: [1])
I.margin(1).value                         # -0.0106972...
I.level_certificate()                     # m = 2 with proofs: level 1 empty, exact level-2 witness
[to_str(g) for g in I.interlacers(2)][:3]

O = Order.maximal("x^6 - x^5 - 7*x^4 + 2*x^3 + 7*x^2 - 2*x - 1")   # needs cypari2
O.min_trace(), O.trace_semigroup(), O.totally_positive_codifferent(1)
```

## How it works

**Enumeration** (`simplex.py`). The trace-`j` totally positive codifferent
elements are the points of the affine lattice `L ∩ {Σx = j}` inside the open
simplex `Δ_j = {x > 0, Σx = j}`, a *regular* simplex in the trace form. The
trace-0 sublattice is computed exactly (unimodular gcd completion), LLL-reduced
(FLINT, on a high-precision Arb embedding), and orthogonalised. A depth-first
search fixes one integer coordinate at a time. At each depth the admissible
interval is read off the **facets of the projection of `Δ_j`** onto the
remaining Gram–Schmidt coordinates (Qhull, computed once per lattice, lazily,
scaled by `j`). This is the exact slice, so every node of the tree is the
projection of a real point of the polytope. For `m_β` (`n = 21`), `I = 0` is
proved with 26 nodes, against 22,680 points in the coefficient bounding box.
For every `f_N` of Theorem 3.1 it also takes 26 nodes, because the reduced
basis "sees" the subfield. Bounds carry tolerances, so rounding can only add
nodes; points within `10⁻⁷ j` of the boundary are decided exactly.

**High degree.** Projections of `Δ` can have too many facets for `d ≳ 24`; those
depths fall back to a valid ball/half-space relaxation. Existence questions are
then decided by **HiGHS MILP** (`scipy.optimize.milp`): infeasibility of
`{k : x(k) ≥ −tol}` proves a level empty (floating point, branch and bound), and
any point found is verified exactly. `level_certificate()` also **lifts
witnesses from subfields**: a totally positive `ν ∈ O_F^∨` lies in `O_K^∨` with
trace `[K:F]·Tr ν`, so `m(O_K) ≤ [K:F]·m(O_F)`. For `F_77` (degree 30) the
degree-15 subfield gives the trace-2 witness in 0.1 s.

**Exactness.** Real roots are certified Arb balls (FLINT), and signs of
`g(r_i)` are decided by ball arithmetic with a gcd test for exact zeros. The
independent check `realroots.interlaces(f, g)` uses Sturm/Cauchy index
`Ind(g/f) = ±d` in exact rational arithmetic only. `bruteforce.py` (LP bounding
box + exact tests) is the reference implementation used in the tests.

**Orders** (`orders.py`). `Order(f, basis)` is any order of `Q[x]/(f)`
(`monogenic`, `maximal` via PARI, `generated_by`, `from_lattice`), with an exact
trace form, codifferent, `m(O)`, the trace semigroup, elements of a given
trace, and margins. `heredity.py` provides `O ∩ F` for a subfield (exact
saturation), the relative trace `Tr_{K/F}`, the trace-down bound
`m(O) ≥ m(O ∩ F)`, the paper's heredity certificate, and subfield lifts.

**Siegel's formula** (`siegel.py`). Exact `b_k(2d)` with
`ζ_K(−1) = 2^d Σ_{k ≤ dim M_{2d}} b_k s_K(k)`,
`s_K(k) = Σ_{ν ∈ 𝔡^{-1}, ν ≫ 0, Tr ν = k} σ_1(ν𝔡)`. `siegel_zeta(f)` evaluates the
right-hand side from the enumeration (ideal factorisation by PARI) and compares
it with PARI's `ζ_K(−1)`. This is an end-to-end completeness test of the
enumeration: F₂₁ gives `s = (0, 16380)` from its 141 trace-2 elements, i.e.
`16/3`.

**Other tools.** `search.py` (Robinson-style census of all totally real
polynomials of given degree and span, random searches), `cyclotomic.py`
(admissible conductors, Chebyshev compositions, the reachable degrees
`D`), `families.py` (`f_N`, towers, Theorem 3.1 checks), `symmetric.py`
(constructive Lemma 4.1: interlacer from an integer symmetric matrix; search
for symmetric realisations), `sailsdata.py` (bridge to
[pav10/sails](https://github.com/pav10/sails)).

## Findings

Full details are in [`results/paper_verification.md`](results/paper_verification.md) and
[`results/higher_interlacing_cyclotomic.md`](results/higher_interlacing_cyclotomic.md).

* **A non-interlacing octic** (the draft leaves degree 8 open):

      f = x^8 - 4x^7 - x^6 + 17x^5 - 5x^4 - 23x^3 + 6x^2 + 9x - 1 = h(x^2 - x),
      h = x^4 - 7x^3 + 15x^2 - 9x - 1   (the totally real quartic of discriminant 725).

  `f` is irreducible and totally real with span 3.798, `O_K = Z[α]`, disc
  `309593125 = 5^4·19·29^2·31`, and `m(f) = m(O_K) = 2`, with `S = <2,3>`. Evidence:
  * the simplex search;
  * an independent LP coefficient search, which produces the same 406 level-2 interlacers;
  * the MILP margin −0.00695;
  * Siegel's formula: `s_K(1) = 0` and `s_K(2) = 137088` give `ζ_K(−1) = 3584/15`, matching PARI.

  Its proper subfields `Q(√5)` and the quartic field have `m = 1`. So this is the
  "genuinely new obstruction" the draft says a degree-8 example would need. It is
  sporadic: `h(x^2 − x − N)` is interlacing for `N = 1..12`. It is the only
  non-interlacing class among the 26 octics of span ≤ 4.
* **Higher interlacing of cyclotomic polynomials** (roots in `[−2,2]`, i.e. products of
  distinct `Ψ_n = minpoly(2cos(2π/n))`):
  * Level 1 is the Beukers–Heckman/McKee–Smyth classification of interlacing
    cyclotomic pairs; the engine agrees with it on every squarefree cyclotomic
    polynomial of degree ≤ 12 (9363 polynomials up to `x → −x`).
  * Consequence: `Ψ_n` is interlacing iff `rad(odd part of n) ∈ {1, p, 15}`.
  * Level `j` ↔ Schur-stable `H ∈ Z[z]` with `lc(H) = j+1` and `|H(0)| = j−1`, via
    the Blaschke-phase argument of Mao–Wang–Zheng (arXiv:2609.34551).
  * Unit resultants split everything: `S`, `m` and the counts `I_j` are
    sums/convolutions over the components of the graph `m ~ n ⟺ m/n` is a prime
    power. So `m(Ψ_{p_1}⋯Ψ_{p_k}) = k`, and levels are unbounded.
  * Composita: for coprime `n = n_1 n_2` with `n_1, n_2 > 2`,
    `m(Ψ_n) ≤ 2 m(Ψ_{n_1}) m(Ψ_{n_2})`. `F_{n_1}F_{n_2}` has index 2 in `F_n`, its
    ring of integers is `O_{F_{n_1}} ⊗ O_{F_{n_2}}`, and the product of two witnesses is a
    witness. Hence `m(Ψ_n) = 2` whenever `n` has two odd prime factors (other than
    `{3, 5}`) or three including 3 and 5.
  * **`m(Ψ_231) = 3`** (`231 = 3·7·11`, degree 60), the first real cyclotomic field
    that is not 2-interlacing. Level 2 is empty, which takes a new method, below. Levels 3
    and 5 come from MILP witnesses, level 4 from a compositum, so `S(Ψ_231) = <3,4,5>`.
    The method: the relative trace maps level-`j` elements of a field onto level-`j`
    elements of each subfield. So level 2 can be enumerated up a tower
    `F_11 ⊂ F_7F_11 ⊂ F_77 ⊂ F_231` (80, 6570, 285 and 0 elements), one fibre per
    Galois orbit (`noninterlacing/fibers.py`). The engine confirms the 6570.
  * Abelian fields: trace-down gives `m(O_K) ≤ m(Ψ_n)` for a real abelian field `K` of
    conductor `n`. So no real abelian field of conductor < 231 has `m ≥ 3`, and `K` is
    interlacing whenever `Ψ_n` is.
  * Every `Ψ_n` with `m = 2` in the table (58 conductors) has `S = <2,3>`.

* **Confirmed:** `I(m_β) = 0` (box of 22,680 points); Theorem 3.1 (i)–(v) for
  `N = 2..30` (total reality from `N = 2`, irreducibility, `I(f_N) = 0`, span,
  `M(f_N) = m_β(−N)` from `N = 4`, the exact discriminant formula); towers in
  degree 24; Lemma 4.1; monogenicity of `F_n`; divisor heredity for
  `n = 84, 63, 99, 147, 189` (all `m = 2`, up to degree 54); the Chebyshev remark
  (`k ∈ {3, 7, 9}`); evenness and sparsity of `D`; the volume formula
  `vol P_1(f) = √|disc f|/(d−1)!`; `|disc f| ≥ d^d` for interlacing `f`.
* **Corrected:** Prop 5.7 omits the conductors `p^k`, `k ≥ 2`, `p ≡ 3 (mod 4)`
  (27, 49, 81, 121, …). The fields checked (`F_27`, `F_49`, `F_81`) are all
  interlacing, so the conclusion survives for them. The tested "perturbations"
  `C_9 ± 1` are `±m_{β_27}(±x)`, i.e. exactly such a field, and `C_9 − x + 1` is
  not totally real. Prop 5.6 holds with threshold 4 instead of 3.
* **Refuted:** Remark 5.9's sign pattern. For `d = 6` both `b_1, b_2 > 0`; for
  `d = 13` (odd) both are negative. In all degrees with `dim M_{2d} = 2`,
  `sign b_2 = (−1)^d = sign ζ_K(−1)`, so the signs never exclude `s_1 = 0`.
* **New data.**
  * Among all 398 totally real sextic fields in the sails/LMFDB list
    (disc ≤ 10⁷), **F₂₁ is the only one with `m(O_K) = 2`**, and
    `m(Z[a]) = m(O_K)` for every polredabs polynomial.
  * The census of irreducible totally real sextics of span ≤ 4.3 (226 classes up to
    `x → ±x + n`) contains exactly one non-interlacing class, `F_21`.
  * A table of `m(F_n)` for the real cyclotomic fields with `n ≤ 400` and degree ≤ 60
    ([`results/cyclotomic_levels.json`](results/cyclotomic_levels.json)). 155 of the
    156 fields are decided, 9 of them by the compositum bound and `n = 231` (`m = 3`)
    by fibres.
  * All 141 real abelian fields of degree 6–24 and conductor ≤ 113
    ([`results/abelian_levels.json`](results/abelian_levels.json)): 126 have `m = 1` and
    15 have `m = 2`. Four of them are decided by trace-down from `Q(ζ_p)^+`
    (`p = 73, 89, 97, 109`) after the direct search ran out of time.
  * Margins: whenever the obstruction is inherited from `F`, the computed
    level-1 margin is `ε_1(K) = ε_1(F)/[K:F]` (`F_84, F_63, F_147, F_189` over
    `F_21`, `F_99` over `F_33`), and `ε_1(f_N) → ε_1(F_21)/2`. Trace-down
    proves `ε_1(K) ≤ ε_1(F)/[K:F]`; the equality is empirical.
  * 22 sextic fields have `m = 1` although no sail facet through `1` has
    integer distance 1.

## Relation to sails

[pav10/sails](https://github.com/pav10/sails) computes sails (convex hulls of
totally positive integers) of totally real orders. Sail facets are given by
totally positive *codifferent* normals `δ_F` with integer distance
`d_F = min Tr(δ_F x)`. For a facet through the vertex `1`, `d_F = Tr(δ_F)`,
hence `m(O) ≤ min{d_F : 1 ∈ F}`, and a facet with `d_F = 1` *is* an interlacer.
`scripts/sextic_census.py` reads sails' published sextic records (read-only).
All 398 fields agree with this bound. The converse fails in 22 fields: there a
trace-one `δ` supports the sail only along a lower-dimensional face through 1.

## Layout

| module | contents |
|---|---|
| `poly.py` | parsing/printing, compositions, Chebyshev `C_k`, `m_{β_n}`, Mahler measure, translation classes |
| `realroots.py` | Sturm sequences, Cauchy index, exact interlacing test, certified root balls |
| `simplex.py` | the enumeration engine, MILP margin and feasibility |
| `interlace.py` | level-`j` interlacers, counts, `m(f)`, level semigroup, certificates |
| `orders.py`, `heredity.py` | orders, codifferents, `m(O)`, subfields, trace-down (bounds and witnesses), lifts |
| `cyclotomic.py`, `families.py` | real cyclotomic bases, admissibility, compositum witnesses, `D`, `f_N`, towers |
| `siegel.py` | modular forms, Siegel coefficients, `ζ_K(−1)` check |
| `search.py`, `symmetric.py`, `sailsdata.py` | census/random search, Lemma 4.1, sails bridge |
| `kronecker.py` | cyclotomic polynomials: McKee–Smyth pairs, components, AM–GM bound, volume heuristic, level-`j` disc polynomials |
| `bruteforce.py` | reference methods: bounding-box scan, LP-bounded coefficient search |
| `paper.py`, `cli.py` | verification report, command line |
| `jobs.py` | one computation in a child process with memory and wall-clock limits (used by the table scripts) |
| `fibers.py` | real cyclotomic fields: Galois action, subfield codifferents, fibres of the relative trace, level-`j` elements up a tower |
| `scripts/` | `sextic_census.py`, `cyclotomic_levels.py`, `abelian_levels.py`, `span_census.py`, `kronecker_census.py`, `fibre_levels.py`, `kronecker_report.py` (produce `results/`) |
