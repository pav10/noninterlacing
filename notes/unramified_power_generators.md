# Power generators in totally real unramified extensions

For a primitive integral element $a$ of a degree-$n$ field $K$, the
polynomial discriminant is $D_K[\mathcal O_K:\mathbf Z[a]]^2$ in absolute
value. Both factors matter at a fixed budget $C^n$.

For perspective, the asymptotic Odlyzko lower bounds for **totally real**
field root discriminants are $4\pi e^{1+\gamma}=60.83\ldots$
unconditionally and $8\pi e^{\gamma+\pi/2}=215.33\ldots$ under GRH
(A. M. Odlyzko, *Bounds for discriminants and related estimates for class
numbers, regulators and zeros of zeta functions*, 1990). Thus the literal
choice $C=200$ serves as a finite-degree benchmark under GRH. An
infinite-degree construction could target a larger absolute constant.

## A local obstruction that applies to every generator

Suppose $s$ embeddings of $K$ into $\overline{\mathbf Q}_p$ take values in
one unramified extension of $\mathbf Q_p$ with residue size $Q=p^f$. At
level $j$ there are $Q^j$ residue classes modulo $p^j$. Put
$(u_j,r_j)=\operatorname{divmod}(s,Q^j)$. Distributing the conjugates as
evenly as possible minimizes pairs in the same class. Counting ordered
pairs gives
$$
 v_p(\operatorname{disc}(a))\ \ge\
 V(s,Q)=\sum_{Q^j<s,\ j\ge1}
       [Q^j u_j(u_j-1)+2r_j u_j].
$$
Each pair in this local field contributes once for each level up to the
valuation of its difference. All other conjugate differences have
nonnegative valuation. If $p$ is unramified in $K$, then $p^{V(s,Q)}$
divides the square of the index, so it multiplies the field discriminant
in the lower bound. Already level one gives $V\ge s^2/Q-s$ for $s\ge Q$.
In any sequence with $|\operatorname{disc}(a)|\le C^n$, this forces
$s=O_{p,f,C}(\sqrt n)$. In particular, a positive proportion of embeddings
in a fixed unramified local field excludes every fixed polynomial budget
as $n$ grows. This also applies when only a subset of places has bounded
residue degree, by taking a common multiple of those degrees for $f$.

**A uniform bound even with unrestricted ramification.** Fix a rational
prime $p$ and let $N_f$ be the number of prime ideals of $K$ of norm
$Q=p^f$. For any primitive integral generator $a$ one has
$$
 v_p(\operatorname{disc}(a))\ \ge\ (fN_f)^2/Q-fN_f.
$$
Take the embeddings inducing precisely those places. Give each embedding
$i$ above a place of ramification index $e_i$ weight $w_i=1/e_i$. Their
total weight is $S=fN_f$, since a place contributes $e_i f$ embeddings.
Their residues lie in the same finite field $\mathbf F_Q$. For a residue
$t$, its Teichmüller lift is the same in every image field containing that
residue (the unique solution of $X^Q=X$ reducing to $t$). Thus a conjugate
of residue $t$ differs from this common lift by an element of valuation at
least $1/e_i$, normalizing $v_p(p)=1$. Two conjugates of the same residue
contribute at least $1/\max(e_i,e_j)\ge w_i w_j$ to the valuation of their
difference. If $S_t$ is the weight in residue $t$, counting ordered
distinct pairs gives at least
$\sum_t S_t^2-\sum_i w_i^2\ge S^2/Q-S$. This proves the inequality.

In particular, if $|\operatorname{disc}(a)|\le C^n$, then
$$
 N_f\le \frac{Q+\sqrt{Q^2+4Qn\log_p C}}{2f}.
$$
This $O_{p,f,C}(\sqrt n)$ bound uses no assumption on splitting or
ramification. It rules out any tower with more than square-root growth
in the number of prime ideals of any one fixed norm, even if their
ramification indices grow. Conversely, small values of these counts
alone give no upper bound on the index of a generator.

In the ordinary Hilbert class field $H$ of a real quadratic field of
fundamental discriminant $D$, the real places split and $n=[H:\mathbf Q]=2h$.
The field discriminant equals $D^{n/2}$. An inert rational prime $p$ gives
a principal prime $(p)$ of the quadratic field. Its Artin symbol is
trivial over the base, so all $n$ embeddings land in the same unramified
quadratic local field. This gives the following lower bounds for every
primitive integral generator.

| $D$ | $n$ | Field root discriminant | $p$ | $V(n,p^2)$ | Polynomial root discriminant at least (rounded down) |
|---:|---:|---:|---:|---:|---:|
| 7057 | 42 | 84.01 (approx.) | 5 | 34 | 309.13 |
| 32401 | 90 | 180.00 (approx.) | 7 | 82 | 1059.88 |
| 36865 | 104 | 192.00 (approx.) | 7 | 122 | 1882.22 |

Each field meets the field budget 200 and each fails the polynomial
budget 200. The class numbers are certified unconditionally with PARI.

## An absolute generator in a cubic-base class field

The totally real cubic field defined by $y^3+y^2-17y+16$ has discriminant
8069 and cyclic class group of order four. Its ordinary Hilbert class field
has degree 12 over the rationals. The intermediate degree-six field and
the degree-twelve field have the same field root discriminant
$8069^{1/3}$. Absolute generator reduction produces the following indices
and polynomial root discriminants.

| Absolute degree | Index | Polynomial root discriminant |
|---:|---:|---:|
| 6 | 7 | 38.36830 |
| 12 | 17 | 32.16237 |

The degree-twelve polynomial is
$$
x^{12}-5x^{11}-3x^{10}+42x^9-25x^8-102x^7+97x^6+68x^5
-70x^4-21x^3+17x^2+3x-1.
$$
Both fields contain the non-Galois cubic base, so neither is contained in
an abelian cyclotomic extension of the rationals. Choosing a fresh
absolute generator can lower the polynomial root discriminant when passing
up a tower. These values are attained upper bounds. The construction gives
no infinite sequence of bounded polynomial root discriminant.

Run `python scripts/unramified_power_generators.py` for the exact budget
checks, certified class groups, class field construction, irreducibility,
real-root counts and index identities.

## Census of quadratic-base ordinary Hilbert class fields at budget 200

Among all real quadratic bases whose ordinary Hilbert class field has
absolute degree $n\ge32$ and field root discriminant at most 200, there
are 181 bases. Local divisibility of the power-basis index rules out
**160** at polynomial discriminant $200^n$. Of the 21 surviving fields,
16 have no common index divisor, and five have one or more forced primes
whose proved contribution still fits the budget. The highest-degree
survivor comes from $D=14401$, with $n=86$. Its quadratic base has cyclic class group of order 43. Conjugation acts
by inversion on this group, so the absolute field is dihedral and hence
noncyclotomic. It has no common index divisor. For this survivor any single prime can be avoided in a suitably
chosen generator's index. The possible sizes of the global indices are
still unknown here.

The search runs over all positive fundamental $D$ from 5 through 40000,
since the absolute field root discriminant is $\sqrt D$. It certifies
each ordinary class group using PARI. At an inert rational prime the local
data are $e=1,f=2$. At a split or ramified rational prime take $f$ to be
the order of the ideal class of a prime of the quadratic base above it,
and $e$ the ramification index in that base. These follow by the Artin
map and unramifiedness over the base. All $s=n/(ef)$ completions have the
same local degrees.

Let $M_t(p)$ be the number of monic irreducibles of degree $t$ over
$\mathbf F_p$, computed by Möbius inversion. The common-index-divisor
criterion gives $p$ as a forced prime exactly when $s>M_f(p)$. Indeed,
index prime to $p$ requires distinct degree-$f$ residue generators at
the $s$ places. Conversely those distinct irreducibles can be used to
generate the product of local integer rings: in each completion add a
uniformizer to a lift of a residue generator from its maximal unramified
subextension. The resulting element generates the local integer ring
(the unramified lift is recovered by Hensel lifting in the complete
generated order, then so is the uniformizer). Distinct residue minimal
polynomials separate the factors by the Chinese remainder theorem.
Approximation modulo $p$ then gives a generator of the global order
locally at $p$. Primes $p>n$ cannot be forced: translating a fixed
degree-$f$ residue generator by $\mathbf F_p$ produces $p$ distinct
degree-$f$ elements, hence $fM_f(p)\ge p>n$.

For the size of the forced index, assign each completion to the residue
minimal polynomial of the chosen generator there. It may have any degree
$t\mid f$. Two completions assigned the same irreducible of degree $t$
contribute at least $2ef^2/t$ to the valuation of the *square* of the
index through the cross resultant of their local minimal polynomials.
Indeed at each of its $t$ roots there are $ef/t$ embeddings per completion;
a pair reducing to the same root differs by an element of valuation at
least $1/e$, by its common Teichmüller lift. Count both orders of each
pair. Discriminants internal to the factors already contain their local
field discriminants. Minimizing this additional cost over all assignments
gives a lower bound for every generator. For each of the $M_t(p)$ bins
use the successive marginal costs $0,2ef^2/t,4ef^2/t,\ldots$ and take
the $s$ cheapest slots. At unramified primes the earlier multilevel
collision bound may be larger, so use the maximum. Multiply the resulting
prime powers into the field discriminant and compare integers.

Run `python scripts/quadratic_hilbert_screen.py` for the complete list of
surviving base discriminants and their forced primes. This is a census
of one class-field construction and of necessary local conditions.
