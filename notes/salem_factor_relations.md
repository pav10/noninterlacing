# Detecting a Salem factor when the other factors are unrestricted

Let $\lambda$ be Salem of degree $n=2d$, $L=\log\lambda$, and write
$T_j=\lambda^j+\lambda^{-j}$. The following statements concern
the irreducible factor containing $\lambda$. The other factors of a
polynomial relation may have any root signature.

1. There are exactly **19** Salem numbers for which $T_a-T_b$ is rational
   for some integers $a>b\ge1$. The least is Lehmer's number. There are
   29 triples $(a,b,c)$ with $T_a-T_b=c$ across those numbers.
2. Only finitely many Salem numbers can be roots of *any* integer
   polynomial with at most five nonzero coefficients, all $\pm1$.
   Neither reciprocity of the polynomial nor cyclotomicity of the cofactor
   is required. We give a finiteness proof for this larger class. Its
   minimum is not determined by the census in statement 1.

## Integer differences at unrestricted exponents

A rational difference is an algebraic integer. It is positive at the
exceptional root and strictly below 4 at any unit-circle conjugate.
Thus $c=1,2,3$. Put $h=a-b$. The identity
$$ (z^h-1)(z^{a+b}-1)=c z^a $$
holds at every conjugate $z$ of $\lambda$. For a unit-circle conjugate it
gives $|1-z^h|\ge c/2$, so $z^h+z^{-h}\le2-c^2/4$. At the exceptional
root,
$$2(\cosh(aL)-\cosh(bL))=c\ge2(\cosh(hL)-1).$$
Thus a short interval for the ordinary conjugates and an upper bound
on $hL$ both concern the *same power*. This is what avoids assumptions
on the cofactor. There is a second pair of bounds:
$z^b+z^{-b}\in[-2,2-c]$ and $c\ge 2hL\sinh(bL)$.

For an interval of width $W$ containing all ordinary trace conjugates
of a power $j$, the [discriminant span estimate](salem_five_unrestricted.md)
gives
$$
 (d-2)\log(4/W)-\log(64(d-1))\le 2jL .
$$
For $j=h$ use $W=4-c^2/4$ and $hL<1$ for $c=1$, $hL<2$ otherwise.
The left side exceeds the permitted value at $d=200$ and is increasing
thereafter, giving $d<200$. The [finite-degree determinant](salem_resultant_bounds.md)
certifies $L>1/100$ individually for $2\le d<200$. It uses multiplicity
8 at exponent 1 and multiplicity 1 at the 15 primes starting with 5.
Hence $bL<5$ for $c=1$ and $bL<6$ otherwise. With $j=b$, $W=4-c$, the
span estimate now gives $d<70$. Multiplicity 6 at exponent 1 and the
first ten primes starting with 5 give $L>1/55$ at each remaining degree.

For $a\ge300$, already
$2(\cosh(a/55)-\cosh((a-1)/55))>3$. This adjacent difference increases
with $a$ and with $L$, hence all exponents have been bounded. No
coprimality condition on $a,b$ was imposed: cases involving proper
powers are included.

## Exact factor census

Let $C_j(z+1/z)=z^j+z^{-j}$. Factor $C_a-C_b-c$ over the integers.
On the real interval above 2 it crosses zero just once, simply, by
the hyperbolic cosine formula. Only its monic irreducible factor $f$
with $f(2)<0$ can yield a Salem number. Test $2\le\deg f<70$ and
count exactly $\deg f-1$ roots in $[-2,2]$ by Sturm. Endpoints cannot be
roots of an irreducible nonlinear factor. The reciprocal lift is then
irreducible, since at an interior embedding its defining quadratic has
nonreal roots over a totally real field.

| $c$ | Pairs left by the exponent bounds | Pairs giving a Salem root | Distinct Salem numbers |
|---:|---:|---:|---:|
| 1 | 2675 | 17 | 13 |
| 2 | 4490 | 8 | 8 |
| 3 | 5960 | 4 | 4 |

There are overlaps between rows. The extra $c=2$ pairs are
$(2,1),(3,1),(3,2),(4,1),(4,3),(5,4),(9,8),(16,15)$, and the $c=3$
pairs are $(2,1),(3,1),(3,2),(7,6)$. The $c=1$ pairs are exactly those
in the [height-one classification](salem_minimum.md), even after allowing
arbitrary cofactors. The three triples for Lehmer are $(6,1,1)$,
$(7,4,1)$ and $(16,15,2)$. In the last example $C_{16}-C_{15}-2$
has the linear factor $x$, Lehmer's trace factor, and an irreducible
degree-ten factor with only six roots in $[-2,2]$. Requiring a pure
Salem-times-cyclotomic witness would discard this valid relation.

There is an exact way to separate the two cofactors in this example.
The new relation lifts to $H=x^{32}-x^{31}-2x^{16}-x+1$. For the
[degree-75 greedy relation](salem_greedy_cofactor.md) $G$, one obtains
$$\gcd(G,H)=(x^2+1)P_{\rm Lehmer}.$$
Both input polynomials have noncyclotomic cofactors; the independent
relation eliminates the degree-48 factor exactly.

Run `python scripts/salem_integer_trace_differences.py` for the
integer factorizations, exact Sturm counts, interval-certified cutoffs
and comparison of the exceptional roots.

## A bound for longer products

More generally suppose
$$\prod_{i=1}^s(\lambda^{m_i}-\epsilon_i)=c\lambda^a,$$
where $s\ge1$, $m_i\ge1$, $a$ is integral, $\epsilon_i=\pm1$ and $c$
is a nonzero integer. Positivity and a unit embedding give $1\le c<2^s$.
Evaluation at $1/\lambda$ forces $2a=\sum m_i$ and
$\prod(-\epsilon_i)=1$. Dividing by $\lambda^a$ now expresses $c$ as
a product of $2\sinh(m_i L/2)$ or $2\cosh(m_i L/2)$. Thus for
$h=\min m_i$, $hL<2\operatorname{arsinh}(1)$. At every ordinary
conjugate each factor has modulus at least $c/2^{s-1}$. Using the factor
at $h$, this confines the $h$-power trace to an interval of width
$4(1-c^2/4^s)$, at either end according to its sign. For $d\ge3$,
$$
 4^s\ \ge\ \frac{d-2}{\log(64(d-1))+4\operatorname{arsinh}(1)}.
$$
This uses $-\log(1-u)\ge u$. Fixed $s$ bounds both the Salem degree
and the house, giving finiteness even when the cofactor is unknown.
There is no assertion that every Salem number satisfies such a product
relation.

## Five terms with an arbitrary cofactor and no symmetry assumption

Here is a separate argument that also handles a nonreciprocal polynomial
witness. We use a simultaneous return property of the *ordinary
conjugates of the Salem factor*.

Consider any sequence of Salem numbers with $n\to\infty$. Select at most
three positive integers $h_j$ per member with $h_j L\le\log4$. Eventually
some ordinary conjugate $z$ satisfies $|z^{h_j}-1|<1/4$ for all selected
$j$. To prove this, use the trigonometric polynomial
$$K(w)=|(1+w)/2|^{2k},\quad |w|=1,\qquad k=1024.$$
Its Fourier coefficients are nonnegative, sum to 1, and its constant
coefficient is $p=4^{-k}\binom{2k}{k}$. In the product of $K(z^{h_j})$
the constant coefficient is at least $p^3$, including when some
frequencies coincide. For each nonzero resulting frequency $q$,
$|q|L\le3k\log4$. The average of $z^q$ over the ordinary conjugates
is $o(1)$ uniformly in these frequencies. Indeed, $\lambda^{|q|}$ is
Salem of the same degree, has bounded house, and its first ordinary
power sum is $O_k(\sqrt{n\log n})$ by the
[logarithmic-kernel estimate](salem_decomposable_traces.md) at
smoothing parameter $1-1/n$. The negative frequencies are conjugates
of positive ones. Thus the average product is at least $p^3-o(1)$.
If the simultaneous return failed at every embedding, the product
would everywhere be at most $(63/64)^k<p^3$ because
$|(1+w)/2|^2=1-|1-w|^2/4$. This proves the claim.

Now take a relation with exactly five coefficients $\pm1$, ordered by
increasing exponent, and remove any monomial factor. At the real root
the highest term is bounded by four times the size of the second
highest, so their exponent difference times $L$ is at most $\log4$.
At the reciprocal root the same holds for the lowest pair.

* If both end pairs have opposite signs, choose an ordinary conjugate
  making both gap powers within $1/4$ of 1. Each pair has modulus below
  $1/4$ there and cannot together cancel the middle term of modulus 1.
* If the highest two signs agree, the highest term is at most three
  times the size of the third highest, as the top pair must be cancelled
  by the remaining three terms. Relative to the third highest exponent,
  both higher gaps times $L$ are at most $\log3$. Apply the simultaneous
  return to those two gaps and the gap of the bottom pair. After removing
  a unit phase from each group, the triple sum is within $1/2$ of an
  odd integer and the pair sum within $1/4$ of an even integer.
  Their absolute values cannot agree. If the lowest two signs agree,
  reverse the witness and use the same proof, since a Salem factor is
  reciprocal.

Both cases contradict unbounded Salem degree. With fewer terms there
is no Salem factor at all. Indeed at a nontorsion unit conjugate two
signed monomials cannot cancel; three unit vectors summing to zero form
a rotated equilateral triple. Four unit vectors summing to zero pair
into antipodes: their monic polynomial has vanishing first and third
elementary symmetric sums, hence is even. Each case would give a
root-of-unity quotient between distinct signed powers. Finally every
root of a monic height-one polynomial has modulus less than 2, giving
finiteness at bounded degree.

This finiteness statement controls arbitrarily large cofactors. It
does not classify all nonreciprocal five-term witnesses. The number of
digits in a general low-height multiple remains unbounded.

## A test for witnesses with unboundedly many terms

There is also a cofactor-independent estimate that uses coefficient
correlations in place of a term bound. For any real integer-polynomial
witness $H=\sum a_i x^i$ divisible by the Salem polynomial, form
$$
 E=\sum_i a_i^2,\qquad D_j=\sum_i a_i a_{i+j}\quad(j>0),\qquad
 J=2\sum_{j>0}|D_j|/E,\qquad
 B=\frac{\sum_{j>0} j|D_j|}{\sum_{j>0}|D_j|}.
$$
Missing coefficients are zero. $J$ is nonzero, as one sees by taking
the lag between the two extreme nonzero coefficients. With $m=n-2$ and
$r=1-1/n$, every such witness satisfies
$$
 2(1+1/m) B\log\lambda
 \ \ge\ \frac{rm}{J^2}-\log n+\frac{m-1}{2}\log r-2\log4 .
$$
At each ordinary conjugate $z$, $|H(z)|^2=0$. Average its Fourier
expansion to obtain $m E\le2\sum_{j>0}|D_j||S_j|$, where
$S_j=\sum_z z^j$. The first Fourier mode of the logarithmic-kernel
bound applied to $\lambda^j$ gives
$$
 r|S_j|^2\le m\log n-\binom m2\log r+2(m+1)jL+2m\log4 .
$$
Use Cauchy–Schwarz with weights proportional to $|D_j|$ and average
these squared bounds. This proves the formula. For instance, $B\le Kn$
and $J\le J_0$ along an unbounded-degree sequence imply
$\liminf\log\lambda\ge1/(2KJ_0^2)$. With $n=10^6$, $B\le n$, $J\le4$
the inequality gives $\log\lambda>0.03124$. No restriction on the
number or size of the coefficients is needed. For a witness with $s$
nonzero terms one also has $J\le s-1$ by Cauchy–Schwarz on the
coefficient vector.

Thus a sequence tending to 1 cannot possess witnesses simultaneously
having bounded $J$ and $B/n$. A low-height multiple from lattice
counting need not control either quantity. Searching for cancellations
of the signed sums $D_j$ would test a different property from counting
nonzero coefficients.

[Fourier modes and shifted relations](salem_fourier_correlations.md)
derives a lower bound for multiples supported in two short blocks and
constructs Salem polynomials for which every multiple has large $J^2B/n$.
