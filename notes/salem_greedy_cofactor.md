# Bounded coefficients and the cofactor in a Salem search

A likely reference for the low-height-multiple question is E. Bombieri
and J. Vaaler, *Polynomials with low height and prescribed vanishing*,
in *Analytic Number Theory and Diophantine Problems* (1987).
The collision argument below is self-contained. The second construction
uses the greedy expansion of 1 in base $\beta$. Its extra nonreal zeros
are related to the beta-conjugates and “lenticular zeros” studied by
J.-L. Verger-Gaugry. We give explicit bounds and certificates here
without making a priority claim.

## What the geometry supplies

A nonzero integer polynomial $P$ with $M(P)<2$ has a nonzero integer
multiple of coefficient height **one**. Remove a monomial factor and
normalize the leading sign first. The leading coefficient has magnitude
one since it is bounded by $M(P)$. Reduce the $2^N$ polynomials with
digits 0 or 1 of degrees below $N$ modulo $P$. In the real vector
space of the quotient use evaluation at the distinct roots and
derivatives through order one less than each multiplicity (real and
imaginary parts for one member of a complex pair). The image of the
integer remainder coefficients is a full lattice, by Hermite
interpolation. The box for an order-$j$ evaluation at $z$ has coordinate
bound $O_P(N^{j+1}\max(1,|z|)^N)$. Counting lattice points by thickening
the box by a fixed fundamental parallelepiped gives
$O_P(N^{O_P(1)}M(P)^N)$ distinct remainders. A collision gives a multiple
with digits in $\{-1,0,1\}$. Restore the removed monomial.

For a Salem polynomial one can require reciprocity as well. In the
totally real trace field use binary sums of $C_j(r)$, $1\le j\le N$,
where $r=\lambda+1/\lambda$ and $C_j(z+1/z)=z^j+z^{-j}$. The $d-1$
ordinary coordinates are $O(N)$, and the exceptional coordinate is
$O_\lambda(\lambda^N)$. The same lattice count gives a collision for
$\lambda<2$. Clear the negative powers of the resulting Laurent
polynomial and normalize the sign to get a monic reciprocal height-one
multiple. Neither counting argument restricts the other factors.

## An exact test on Lehmer's number

Use $r_0=1$, $d_j=\lfloor\beta r_{j-1}\rfloor$,
$r_j=\beta r_{j-1}-d_j$. For Lehmer's number the orbit has preperiod
1 and period 74. Up to position 75 the ones are exactly at
$1,12,31,44,63$. The canonical numerator of
$1-\sum d_j z^j$ is therefore
$$
 H=1-z-z^{12}-z^{31}-z^{44}-z^{63}-z^{74}+z^{75}.
$$
Exact division gives
$$H=P_{\rm Lehmer}\Phi_2\Phi_4\Phi_{12}\Phi_{22} Q_{48}.$$
The degree-48 factor is irreducible and reciprocal. Of the 24 roots of
its trace polynomial, exactly 20 are real, all strictly between $-2$ and
2. Hence it contributes four nonreal roots inside the unit disk and four
outside. Taking a gcd of $H$ with its reverse leaves the entire
degree-48 factor. In contrast, $(x^2-x+1)P_{\rm Lehmer}$ gives a
five-term cyclotomic completion of the same Salem polynomial. The
choice of multiple matters even for the smallest known example.

## An unavoidable cofactor for the greedy construction near 1

For $1<\beta<2$, let $n$ be the position of the second one in the greedy
expansion of 1. It exists, with $n=1+\min\{g\ge1:
\beta^g(\beta-1)\ge1\}$. After every subsequent one the remainder is
at most $\beta-1$. Thus the positions of successive ones differ by at
least $n-1$. The greedy series has the form
$$F(z)=1-z-z^n-\sum_{j\ge 2n-1}d_jz^j.$$
On $|z|=r<1$ the absolute tail is at most $r^{2n-1}/(1-r^{n-1})$.
This applies equally to a finite digit string.

If $n\ge18$, **every** such series with this gap constraint has a
nonreal zero inside the unit disk. The condition holds for
$1<\beta<1.1339024903\ldots$, where the endpoint solves
$x^{16}(x-1)=1$. For $n\ge4096$ there are at least
$$2K_n=2(\lfloor n/128\rfloor-\lceil n/256\rceil+1)$$
distinct nonreal disk zeros of modulus at most
$\exp[-(\log 20-1/2)/n]$.

For a Salem base the only open-disk zero of its minimal polynomial is
the positive real $1/\beta$. For an eventually periodic expansion
(including termination) the zeros just counted must therefore belong to
the **noncyclotomic** part of the canonical cofactor $Q$. Clearing a
period denominator $1-z^p$ does not change open-disk zeros. For
$n\ge4096$ this proves
$$
 \deg Q_{\rm noncyc}\ge 2K_n,\qquad
 \log M(Q)\ge (2K_n/n)(\log20-1/2).
$$
The second bound is Jensen's formula at the unit circle using the unit
constant coefficient. In a hypothetical sequence of such bases tending
to 1, the complementary log measure has liminf at least
$(\log20-1/2)/128$. This conclusion uses periodicity only when speaking
of a polynomial cofactor. Eventual periodicity for arbitrary Salem
bases is not assumed.

### Disks for unbounded $n$

Put $z=\exp(-w/n)$. For each integer $n/256\le k\le n/128$ set
$A=1-\exp(-2\pi i k/n)$ and $w_0=2\pi i k-\log A$ with the principal
logarithm. We have $1/64\le|A|<1/20$ by the elementary sine bounds.
Use $|w-w_0|=1/2$ and compare to $A-e^{-w}$. Its boundary modulus is
at least $(2-e^{1/2})|A|>|A|/3$ by the exponential series. It has just
one zero in this disk. Throughout the disk $0<\Re w<5$, and the error in
replacing $1-e^{-w/n}$ by $A$ is at most $7/n\le|A|/8$, by integration
along a segment in the closed right half-plane. Also $|z|^n\le
e^{1/2}|A|$ and $|z|^{n-1}\le1.7|A|$. The digit tail consequently has
modulus below $|A|/6$. Rouché gives one zero in each disk. The disks
are disjoint (successive imaginary parts of their centers differ by
$2\pi-\pi/n$). They lie in $0<\Im w<\pi n$, giving distinct nonreal
zeros in one half-plane under exponentiation. Conjugating doubles the
count. The stated modulus bound follows from $\Re w\ge\log20-1/2$.

### Finite range of initial gaps

For each $46\le n<4096$ the script selects a rational complex center
near a nonreal zero of $G(w)=1-e^{-w/n}-e^{-w}$. On a disk of radius
$s=0.35$ compare $F(e^{-w/n})$ to $G'(w_0)(w-w_0)$. With
$a=e^{-w_0/n}, b=e^{-w_0}, l=\Re w_0-s$ and $q=e^{-(1-1/n)l}$ the
total error is at most
$$
 |1-a-b|+|a|(e^{s/n}-1-s/n)+|b|(e^s-1-s)+e^{-l}q/(1-q).
$$
Arb verifies that this is strictly below $s|a/n+b|$ and the disk
lies in $\Re w>0$, $0<\Im w<\pi n$. For $18\le n\le45$ the comparison is made adaptive. For a prefix with
ones at the positions $E$ use $G(w)=1-\sum_{j\in E}e^{-jw/n}$ and the
same linear Taylor comparison, bounding the remainder by the sum of the
individual exponential remainders. If the next one occurs at position
$J$ or later, the tail is bounded by $e^{-Jl/n}/(1-q)$. For each
prefix the script certifies an integer threshold $J$ beyond which the
linear term dominates. It recurses on every earlier possible next one
from $\max E+n-1$ through $J-1$. A terminating expansion is covered
by the zero tail. All branches close, with 355 prefixes across the 28
initial gaps. This proof uses just the gap inequality and includes
arbitrary real bases, without an assumption on algebraicity. The
floating-point Newton iteration only selects the centers. Run `python scripts/salem_greedy_cofactor.py`.

## Implication for the minimum problem

Distinct small Salem numbers would have a strong consequence if they
tended to 1: Salem numbers would be dense in $(1,\infty)$. Every
positive integral power of a Salem number is Salem (the exceptional
conjugate is unique, so taking powers preserves degree). For any fixed
$1<a<b$ and $1<\lambda<b/a$, the first power strictly exceeding $a$
lies in $(a,b)$. Thus a single Salem-free open interval away from 1
would already give a uniform gap at 1. Attainment of the infimum above
that gap is a further question.

A bounded-digit relation is insufficient to certify such an interval,
even with reciprocity. For any fixed $t\in(1,2)$ truncate its greedy
series after $N$ digits to $A_N$ and form the monic reciprocal height-one
polynomial $A_N(z)+z^{2N+1}A_N(1/z)$. Its two coefficient blocks do
not overlap. On compact subsets of the disk these converge to the
greedy series, which changes sign at $1/t$ (the first digit is one).
The intermediate value theorem and reciprocity give positive real
roots arbitrarily close to $t$. What matters for a Salem exclusion is
which of the disk zeros belong to the same irreducible factor as $1/t$.
Lehmer's computation shows precisely how that inference can fail.
