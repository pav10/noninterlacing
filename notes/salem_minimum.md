# Attainment of a Salem minimum under a sparsity hypothesis

Write $P_\lambda$ for the monic minimal polynomial of a Salem number.
Allow multiplication by an arbitrary product of cyclotomic polynomials
with multiplicities, including the empty product. Define $s(\lambda)$
as the least number of nonzero coefficients of any resulting polynomial.
This measures the whole product, and the Salem factor occurs once.

**Finiteness.** Let $\theta=1.324717957\ldots$ be the real root of
$x^3-x-1$, the smallest Pisot number by Siegel's theorem. For each fixed
$s$ and each $B<\theta$, only finitely many Salem numbers
$\lambda\le B$ satisfy $s(\lambda)\le s$. In particular, for every
$s\ge5$ the Salem numbers with $s(\lambda)\le s$ have a smallest member:
the set contains Lehmer's number, and its part at or below Lehmer is finite.
There are no restrictions on the degrees, exponents, or coefficient sizes.

**An exact minimum.** With the extra conditions of coefficient height one
and at most five nonzero coefficients in the product, there are precisely
13 Salem numbers. Lehmer's number is the smallest. A finite certificate
and its degree-independent cutoff are given below. The subsequent
[arbitrary-coefficient classification](salem_five_unrestricted.md) proves
the same five-term minimum and identifies all six Salem numbers up to
$5/4$ without a height restriction.

These statements address attainment as well as separation from 1. If the
full Salem set had an unattained infimum, a strictly decreasing sequence
toward it could be chosen below Lehmer's number, unless Lehmer itself
already attained the infimum. Its minimum term count would have to tend
to infinity. The arguments supply no bound on that count in general.

## Accumulation with bounded term count

Here are details of the finiteness argument using the
Boyd–Lawton specialization theorem, the multivariate Kronecker theorem,
and Mahler's derivative inequality. They are included to distinguish the
use of those classical results from the additional role of the Salem
signature.

Take distinct Salem numbers $\lambda_j$ converging to a finite limit
$\beta$, with witnesses $P_j$ of at most $s$ terms. Write $M$ for
multiplicative Mahler measure, so $M(P_j)=\lambda_j$.
For any polynomial $P=\sum_{i=0}^{k-1}c_i x^{e_i}$ with increasing
exponents starting at zero,
$$ |c_i|\le {k-1\choose i}M(P). $$
For completeness the endpoints follow from Jensen's formula. For an
interior $i$, induct on $k$ using $M(P')\le e_{k-1}M(P)$, first on $P$
and then on its reverse. Removing monomial factors leaves Mahler measure
unchanged. With $d=e_{k-1}$ the two bounds on $|c_i|/M(P)$ are
$d\binom{k-2}{i-1}/e_i$ and $d\binom{k-2}{i}/(d-e_i)$. Their weighted
average with weights $e_i/d$ and $1-e_i/d$ is $\binom{k-1}{i}$.
Thus we can pass to a subsequence with fixed coefficients. Degrees go
to infinity after passing to a subsequence, since only finitely many
polynomials of bounded degree occur.

Pass again so each successive exponent gap is constant or tends to
infinity. On compact subsets of the open unit disk $P_j$ converges to
the polynomial $A(x)$ consisting of the terms before the first growing
gap. Its constant term is $\pm1$. If $\beta>1$, the only open-disk zero
of $A$ is $1/\beta$, with multiplicity one. This follows from Rouché's
theorem on small circles, since the only open-disk zero of $P_j$ is
$1/\lambda_j$. Hence the signed reverse of $A$ is monic and its only
root outside the closed unit disk is $\beta$. In particular $M(A)=\beta$.
The minimal polynomial of $\beta$ makes it a Pisot or a Salem number:
if it has a unit-circle conjugate, it shares a root with its reciprocal
and is therefore reciprocal. In this case the other conjugates apart
from the reciprocal pair must all be on the circle.

### Retaining the gaps in a multivariate polynomial

Write the growing gaps as $g_{1j},\ldots,g_{hj}$. Every exponent of
$P_j$ is a fixed integer linear combination of the coordinates of
$v_j=(1,g_{1j},\ldots,g_{hj})$. Pass to a diagonal subsequence on which
for every $c\in\mathbf Z^{h+1}$ the assertion $c\cdot v_j=0$ stabilizes.
The eventual relations form a saturated lattice $R$. Its finite basis
vanishes simultaneously for large $j$. Set $H=\mathbf Z^{h+1}/R$.
This is free abelian and $v_j$ gives a specialization $H\to\mathbf Z$.
The first basis vector has image $u\in H$ sent to 1, so $u$ is primitive.
In a fixed basis of $H$ the shortest nonzero relation for these
specializations tends to infinity, since no fixed nonzero element can
specialize to zero infinitely often. Lifting the exponents gives a single
integer Laurent polynomial $F$ in the group ring of $H$ specializing to
$P_j$. Boyd–Lawton now gives $M(F)=\lim_j M(P_j)=\beta$.

If $\beta=1$, multivariate Kronecker makes $F$ a signed monomial times
a product of cyclotomic polynomials in monomials. Its nonzero
specializations have no nonzero roots off the unit circle, contrary to
$P_j(\lambda_j)=0$. This rules out convergence to 1.

For $\beta>1$ we can present $F$ as a polynomial in one variable $t$
with constant coefficient $A(z)$, where $z$ corresponds to $u$. In
detail, subtract a fixed homomorphism taking $u$ to 1 from a sufficiently
late specialization. The difference kills $u$ and is strictly positive
on every growing-gap vector. Divide by its content to get a surjective
integer grading. Split off its degree-one direction $t$ and extend $u$
to a basis of its kernel. The other kernel variables are denoted $w$.
The terms of degree zero are exactly the initial cluster. Write
$F=\sum_{i=0}^D B_i(z,w)t^i$ with $B_0=A(z)$.

### Why a Salem number cannot be the limit

Jensen in $t$ gives $M_t(F)\ge |A(z)|$ on the $(z,w)$ torus.
Since $M(F)=M(A)$, equality holds almost everywhere by integration of
logarithms. At those points with $A(z)\ne0$ the polynomial in $t$ has
no root in the open unit disk. Factoring it from its constant term gives
$$ |B_i(z,w)|\le {D\choose i}|A(z)|. $$
Continuity gives this inequality everywhere on the torus, including the
zeros of $A$. If $\beta$ were Salem, $A$ would have a root $\xi$ that is
a unit-circle conjugate of $\beta$. Thus $B_i(\xi,w)=0$ identically in
$w$ for every $i$. Every Laurent coefficient in $z$ is divisible by the
minimal polynomial $P_\beta(z)$ (clear monomials first). The divisibility
survives each specialization because $z$ always specializes to $x$.
Thus $P_\beta$ divides every sufficiently late $P_j$. Each has just one
root above one, forcing $\lambda_j=\beta$, contrary to distinctness.

It follows that every finite accumulation point in a fixed term-count
class is Pisot. Compactness and Siegel's theorem give the stated
finiteness below $\theta$. The isolation argument uses the unique
off-circle pair; it does not assert discreteness for general sparse
Mahler measures.

## Consequence for dense coefficient patterns

The sparsity condition can be checked even on polynomials with many
terms. Suppose a Salem polynomial times a cyclotomic product is a sum
of $r$ blocks of the form
$$ c x^a(1+x^h+\cdots+x^{(l-1)h}),\qquad h,l\ge1. $$
The block coefficients, starts, strides, and lengths may all vary.
Multiply by the product of the $r$ binomials $x^h-1$, one for each
block. Each block telescopes against its own binomial and the other
$r-1$ factors produce at most $2^{r-1}$ copies. This gives a monic
witness with at most $r2^r$ nonzero terms. Consequently any decreasing
sequence witnessing an unattained infimum must escape every bounded
number of such blocks, even after cyclotomic multiplication. In
particular this covers arbitrarily long constant runs of coefficients.

## Height one and five terms

A monic Salem polynomial times cyclotomic factors is reciprocal or
antireciprocal. With coefficient height one, the two-term and three-term
possibilities are cyclotomic. A four-term polynomial has the shape
$x^d+u x^{d-i}+v x^i+w$, where $w=\pm1$, $u=\pm1$, and $v=wu$.
It factors as $(x^{d-i}+wu)(x^i+u)$ with constant $w$.
An antireciprocal polynomial cannot have an odd number of terms.
Thus a five-term witness has degree $2a$ and the shape
$$ x^{2a}+u(x^{a+b}+x^{a-b})+v x^a+1,\qquad 0<b<a. $$
Its value at 1 is nonzero and must be negative, since exactly one real
root exceeds 1 and the polynomial is monic. This forces $u=v=-1$.
The trace polynomial is $C_a-C_b-1$, where
$C_m(2\cos t)=2\cos(mt)$. If $g=\gcd(a,b)>1$, the original polynomial
is a polynomial in $x^g$, and rotating its largest root produces other
roots outside the circle. We may therefore take $a,b$ coprime.

**A uniform cutoff.** If $a\ge13$, choose $0<j<a/2$ with
$bj\equiv\pm1\pmod a$ by reflecting a modular inverse if necessary.
On the cell $t=2\pi j/a+\delta$, $|\delta|\le\pi/a$, put
$U(t)=2\cos(at)-2\cos(bt)$. Since
$\cos(a\delta)\le\cos(b\delta)$ there,
$$ U(t)\le 2|\cos(b\delta)-\cos(b\delta\pm2\pi/a)|
          \le4\pi/a<1. $$
At both ends of the cell $U\le0$, while at its center $U>0$. Thus
$U-1$ has a strictly negative interior local maximum. Under $x=2\cos t$
this is also a negative local maximum of the trace polynomial, in
$(-2,2)$. A real-rooted polynomial cannot have such a maximum: away from
its zeros its logarithmic derivative has strictly negative derivative,
so at a negative stationary value its second derivative is positive.
Thus the trace fails total reality. A Salem times cyclotomic polynomial
would have totally real trace.

For $a\le12$ exact factorization and Sturm counts give 17 parameter
pairs and 13 distinct Salem numbers. The least occurs at $(6,1)$ and
$(7,4)$. In the first case the pentanomial is
$(x^2-x+1)P_{\rm Lehmer}(x)$. The certificate compares isolated
exceptional trace roots, using monotonicity of
$(r+\sqrt{r^2-4})/2$. Run `python scripts/salem_five_terms.py`.

## Six height-one terms attain the first possible accumulation point

The threshold $\theta$ in the finiteness result is sharp. The elementary
Salem–Boyd construction with $A=x^3-x-1$ and $A^*=1-x^2-x^3$ gives
$$ P_m=x^m A-A^*,\qquad m\ge8. $$
These polynomials have exactly six nonzero coefficients, all of absolute
value one. Since $P_m(1)=0$, $P_m'(1)=7-m<0$, and $P_m(\theta)>0$,
there is a real root $\lambda_m$ in $(1,\theta)$. For $0<t<1$,
Rouché on the unit circle gives $x^m A-tA^*$ exactly one outside root:
$|A|=|A^*|$ there and $A$ has just one outside root. Taking $t$ to 1
shows that $P_m$ has at most one outside root, counted with multiplicity.
Antireciprocity gives exactly one inside root, $1/\lambda_m$. Because
$P_m$ has unit constant coefficient, the irreducible factor of
$\lambda_m$ must also contain the inside root to have unit norm. Its
degree exceeds two: in degree two the trace $\lambda_m+1/\lambda_m$
would be an integer strictly between 2 and 3. It therefore has
unit-circle conjugates and is Salem. The remaining factors are cyclotomic.

At $m=8$, $P_8=(x-1)P_{\rm Lehmer}$. At $\lambda_m$ the value of
$P_{m+1}$ is $\lambda_m^m(\lambda_m-1)A(\lambda_m)<0$, giving
$\lambda_{m+1}>\lambda_m$. For every fixed $r$ in $(1,\theta)$ one has
$P_m(r)<0$ eventually, so $\lambda_m$ increases to $\theta$.
Thus allowing a sixth height-one term gives infinitely many Salem
numbers *below* the smallest Pisot number, but their first possible
accumulation point is exactly that Pisot number. This is why the
finiteness theorem uses each fixed $B<\theta$.

### References for the classical inputs

- D. W. Boyd, *Kronecker's theorem and Lehmer's problem for polynomials in several variables*, Journal of Number Theory 13 (1981), 116–121.
- W. M. Lawton, *A problem of Boyd concerning geometric means of polynomials*, Journal of Number Theory 16 (1983), 356–362.
- K. Mahler, *On the zeros of the derivative of a polynomial*, Proceedings of the Royal Society A 264 (1961), 145–154.
- C. L. Siegel, *Algebraic integers whose conjugates lie in the unit circle*, Duke Mathematical Journal 11 (1944), 597–602.
