# A power-trace bound for Salem numbers with decomposable trace polynomial

Let $\lambda$ be a Salem number of degree $n=2d$. Suppose the degree-$d$ minimal polynomial of $\lambda+\lambda^{-1}$ decomposes over $\mathbf Q$ as $g(T(x))$, with inner degree $b\ge2$ and outer degree $q=d/b\ge2$. Both polynomials can be taken monic integral. Indeed, normalize the inner polynomial to be monic with zero constant term. Its other coefficients are rational symmetric functions of the algebraic-integer roots of any one fibre. The roots of the normalized outer polynomial are then integral as well.

Two consequences are:

1. If $n\ge640 b^2\log n$, then $\log\lambda\ge\log(n/(4b))/(b-1)$. In particular this beats the *order* $(\log\log n/\log n)^3$ along any such sequence with $b=o((\log n)^4/(\log\log n)^3)$. This is a structural restriction. It leaves indecomposable traces untouched and does not supply an unrestricted asymptotic improvement.
2. For quadratic inner degree the sharp minimum is $\lambda=1.547196965626387\ldots$, attained by the degree-eight reciprocal lift of $(x^2-x)^2-3(x^2-x)+1$.

## Arithmetic obstruction in a low power trace

Write $C_j(z+z^{-1})=z^j+z^{-j}$. There is some $1\le k<b$ with
$$|\mathrm{Tr}_{\mathbf Q(\lambda)/\mathbf Q}(\lambda^k)|\ \ge\ q k.$$
To see this, compare the nonconstant coefficients of $T$ to those of the monic $C_b$ in descending degree. If all agree, the exceptional fibre would solve $C_b(x)=C_b(\lambda+\lambda^{-1})>2$. For even $b$ this has at least two roots outside $[-2,2]$, and for odd $b>1$ only one of its $b$ roots is real. Both contradict the Salem signature.

Let the first discrepancy be the coefficient of $x^{b-k}$ and denote its nonzero integer difference by $a$. Newton's identities in a fibre $T(x)-y$ show that its power sums below $k$ agree with those of $C_b(x)$, and the $k$th differs by $-ka$. The fibre sum of $C_k$ over the zeros of $C_b$ is zero, by the cosine formula, or by summing over the zeros of $z^{2b}+1$. Consequently each fibre contributes exactly $-ka$ to the global trace. Summation over the $q$ fibres proves the assertion.

There is also an upper bound on the exceptional trace $r_*=\lambda+\lambda^{-1}$. Every fibre has the same root sum. One entirely interior fibre gives a sum of absolute magnitude strictly below $2b$. In the exceptional fibre the other $b-1$ roots are each larger than $-2$. Thus $\lambda<r_*<4b-2$. This will allow use of $\log\lambda<\log n$ even when no small-measure bound is assumed.

## A discriminant estimate for the unit-circle power sums

Let $m=n-2$ and let $z_i$ range over the unit-circle conjugates. Put $S_k=\sum z_i^k$. For $0<r<1$ the logarithmic kernel gives
$$
\sum_{k\ge1}\frac{r^k}{k}|S_k|^2
 =-\sum_{i,j}\log|1-rz_i\overline{z_j}|.
$$
Off the diagonal, $|1-rz|^2\ge r|1-z|^2$ on the unit circle. The omitted factors in the full polynomial discriminant have absolute product at most $\lambda^2(4\lambda)^{2m}$, since the two exceptional roots are $\lambda,\lambda^{-1}$. The full discriminant is a nonzero integer. Hence, if $\lambda\le B$,
$$
 E=m\log\frac1{1-r}-\binom m2\log r+2(m+1)\log B+2m\log4
 \ \ge\ \sum_{k\ge1}\frac{r^k}{k}|S_k|^2.
$$
Combining with the forced arithmetic trace yields, for at least one $1\le k<b$,
$$
 qk\ \le\ \lambda^k+1+\sqrt{k E/r^k}. \tag{1}
$$

Take $r=1-1/n$. Using $B=n$, $n\ge8$, and $k<b\le n/4$, one gets $E\le5n\log n$ and $r^{-k}<2$. Thus $n\ge640b^2\log n$ makes the square-root term at most $nk/(8b)$, and also $1\le nk/(8b)$. This proves the first stated bound. For fixed $b$ it also proves finiteness: $\lambda<4b-2$ bounds the degree via (1), and bounded degree and house bound the integer coefficients. For fixed $b$ the estimate $\lambda^{b-1}\ge n/(2b)-O_b(\sqrt{n\log n})$ follows from the same calculation.

The sharper numerical test substitutes $B=5/4$ and excludes $b$ if (1) fails for every $1\le k<b$. At Salem degree 4096 it excludes all integer $b$ from 2 through 9. At degree 1048576 it excludes 2 through 59, and at degree $10^8$ it excludes 2 through 80, even before applying $b\mid n/2$. No assertion of existence at the first unexcluded degree is implied.

## Sharp result for inner degree two

Write $T=x^2-a x$ after translating the value of $T$. Each pair of roots has sum $a$. From an interior pair $|a|<4$, and the partner of $r_*>2$ must exceed $-2$, so $a\in\{1,2,3\}$.

Take $g=y^2-3y+1$ and $a=1$. The quartic trace polynomial is irreducible (also checkable by reduction or integer factorization), with three roots in $(-2,2)$ and exceptional root
$$r_0=(1+\sqrt{7+2\sqrt5})/2<11/5.$$
Its reciprocal lift is irreducible because at an interior embedding the defining quadratic has nonreal roots while the trace field is totally real. The resulting Salem number is $(r_0+\sqrt{r_0^2-4})/2$.

To prove minimality it suffices to consider $2<r_*\le r_0$. For $a=1$ put $u=y-1$. Its interior conjugates lie in $[-5/4,1)$, since both roots $x,1-x$ are interior. The exceptional conjugate is $r_*^2-r_*-1\in(1,(1+\sqrt5)/2]$. All conjugates are in $(-2,2)$, so Kronecker gives $u=2\cos(2\pi j/h)$ for a primitive $h$th root. The largest conjugate restricts $h$ to $7,8,9,10$. For $7,8,9$ the least conjugate is below $-5/4$. Conductor 10 yields exactly the displayed example.

For $a=2$ take $u=y$, and for $a=3$ take $u=y+2$. The interior conjugates lie in $[-1,0)$ and $[-1/4,0)$ respectively. The exceptional conjugate is respectively $r_*(r_*-2)$ or $(r_*-1)(r_*-2)$. Both are strictly between zero and $11/25$. A real cyclotomic trace with positive largest conjugate has largest conjugate at least $2\cos(2\pi/5)=(\sqrt5-1)/2$, excluding both cases (the endpoint conductor 1 also satisfies that bound).

Run *python scripts/salem_inner_degree.py* for the exact root count, checked numerical interval for the minimum and finite-degree energy exclusions. This argument uses discriminant integrality and Newton sums. The [orbit-block estimate](salem_orbit_blocks.md) addresses what can be gained directly in the Dobrowolski determinant.

## The fibre through 2 and a finite outer search

Suppose the trace polynomial decomposes as above. Translate the *value* of
the monic integral inner polynomial so that $T(2)=0$. Then
$$
 T(x)=(x-2)\prod_{m\in I}\Psi_m(x),\qquad m\ge3,
$$
where $I$ is a set and $\Psi_m$ is the minimal polynomial of
$2\cos(2\pi/m)$. To prove this, choose a fibre with all $b$ roots in
$(-2,2)$. Rolle's theorem puts all $b-1$ critical points in that interval
and makes them simple. Each monotone branch also has one root of the
exceptional fibre. Its rightmost root is $r_*>2$ and the other roots are
interior. On the increasing rightmost branch the level at 2 lies strictly
between the two fibre levels. On every other branch the root at that
intermediate level lies strictly between the corresponding roots. Thus
$T/(x-2)$ is squarefree with all roots in $(-2,2)$ and Kronecker's theorem
gives the product. The outer roots consist of one positive value $T(r_*)$
and $q-1$ negative values.

For $r_*\le R$ a lower bound on every outer root is the largest local
minimum of $T$, and an upper bound is the smallest local maximum (if
present) or $T(R)$, whichever is smaller. The leftmost branch gives an
additional lower bound $T(-2)$ for odd $b$ and an upper bound for even $b$.
All outer conjugates obey all these inequalities. There is also an exact
shortcut. A negative outer root below $-1$ forces $T+1$ to have $b$ simple
interior roots. A positive outer root above 1 forces $T-1$ to have $b$
simple real roots and $T(R)>1$. These follow by moving between fibre
levels on each branch. If neither tail is possible, all outer roots lie
in $[-1,1]$, contradicting their nonzero integral norm since $q\ge2$ and
none is rational.

When the remaining enclosure fits in $[c-2,c+2]$ for some
$c\in\{-1,0,1\}$, Kronecker gives $g(y)=\Psi_m(y-c)$. In fact $m\le42$.
For $m>42$ let $p$ be the first prime not dividing $m$. If $p<m/6$, both
primitive exponents 1 and $p$ give distinct conjugates strictly greater
than $c+1\ge0$, impossible. Hence $m\le6p$. If $p\ge11$, the product of
all smaller primes divides $m$ and exceeds $6p$. This holds already for 11
and propagates to the next prime by Bertrand's postulate. Hence $p\le7$,
again a contradiction.

Using $R=41/20$ (measure $5/4$), the shortcut and interval tests covered
every endpoint fibre for $2\le b\le18$. No admissible irreducible composite
remained. The numbers of endpoint products before screening, in increasing
order of $b$, were 3, 7, 17, 35, 67, 127, 225, 385, 651, 1065, 1705, 2701,
4191, 6409, 9711, 14505 and 21427. The outer degree is unrestricted.

For $b=3$ increasing $R$ to $23/10$ gives exactly one trace polynomial,
$$
 T=(x-2)(x^2+x-1),\qquad f=T^2-T-1 .
$$
Its reciprocal lift has degree 12 and measure $1.648593855178627\ldots$.
Thus this is the sharp minimum over *all* cubic-inner decompositions.
The trace cap exceeds the exceptional root of the displayed example,
and measure increases with that root. Irreducibility of the lift follows
from an interior trace embedding as in the quadratic case.

Run `python scripts/salem_inner_fibres.py` to reproduce both certificates.
The extrema are enclosed by Arb and the final tests use exact
irreducibility, Sturm counts and rational endpoint evaluation. The script
raises an error for an enclosure it cannot cover with the indicated
translates. The argument gives no bound on indecomposable traces.

There is also a degree descent for any $b$ when $2<r_*<3$. If $T(r_*)>1$,
the rightmost root of $T-1$ is strictly between 2 and $r_*$ and all its
other roots are interior. Its irreducible factor containing that root has
degree at least two since $(2,3)$ contains no integer. Its reciprocal lift
therefore gives a strictly smaller Salem number. Consequently a global
minimum in this range with a decomposable trace would require
$T(r_*)<1$ (equality would give a rational outer root). Its nonzero
integral norm then forces some negative conjugate below $-1$. In this case
both $T$ and $T+1$ are squarefree products of real cyclotomic factors.
This identifies an arithmetic condition on endpoint products for a possible
minimum even beyond the computed degrees.

The same reduction applies to a putative sequence tending to 1 even if
there is no global minimum. Repeat the descent whenever the trace admits
a decomposition with $T(r_*)>1$. The trace degree strictly decreases at
each step, so each member yields a Salem number at most as large with
either an indecomposable trace or an endpoint pair $T,T+1$ as above.
In the latter case every factor $\Psi_m$ of $T+1$ has $m$ divisible by at
least two distinct primes. Indeed $(T+1)(2)=1$, whereas
$\Psi_{p^k}(2)=\Phi_{p^k}(1)=p$ for $p^k\ge3$ and every factor has positive
integral value at 2. The cofactor condition holds at unrestricted inner
degree.

## Classification at every measure for inner degrees 3 through 20

There are four normalized cubic-inner examples, one quartic-inner example
and none with inner degree 5 through 20. Both the outer degree and the
Salem measure are unrestricted. For the four cubics the same inner
polynomial $T=(x-2)\Psi_5$ occurs. The outer polynomials, in increasing
order of Salem measure, are
$$
 y^2-y-1,\quad y^2-y-3,\quad y^2-2y-1,\quad y^2-2y-2.
$$
The unique quartic has $T=(x-2)\Psi_7$ and outer polynomial $y^2-y-1$,
giving Salem degree 16 and measure $1.5042743089473367\ldots$.

The root-sum bound from the first section gives $r_*<4b-2$, so testing
through $4b$ exhausts all measures. For each endpoint fibre use the
translated cyclotomic test above when it applies. Every remaining viable
interval in the calculation satisfies at least one of the following tests.

- All negative outer roots lie in $[-17/12,0]$ and the positive root is
  at most $10/3$. The absolute value of $h(y)=y(y+1)^2$ on the negative
  interval is at most $1/4$ (check the endpoints and the critical points
  $-1,-1/3$). On the positive interval it is at most $1690/27<64$.
  If the outer degree were at least four, the nonzero integral norm of
  $h(y)$ would have absolute value below one.
- All negative outer roots lie in $[-1/2,0]$ and the positive root is
  at most 7. The same norm argument for $y$ gives degree at most three.

The script enumerates all quadratic and cubic integer coefficient vectors
in those boxes by elementary symmetric bounds, then uses exact
irreducibility and Sturm counts. An uncovered interval raises an error.
The last two degrees involved 31405 and 45539 endpoint products.
Run `python scripts/salem_inner_fibres.py --all-measures --max-degree 20`.

### The outer polynomial at unrestricted inner degree

In any such decomposition all negative outer roots are strictly greater
than $-4$. Otherwise the ordinary fibre at the least outer root (an
irrational number below $-4$) and the endpoint fibre at zero force $T+2$
to alternate below $-2$ and above $2$ at the $b-1$ critical points, with
the endpoint at $-2$ strictly on its corresponding side. This follows
from the monotone branches and the interior positions of the ordinary
roots. At the right endpoint $T(2)+2=2$. Subtract the monic
$C_b(x)=2\cos(b\arccos(x/2))$, whose magnitude on $[-2,2]$ is at most 2.
The difference has a zero in each interval between consecutive points of
the list consisting of $-2$ and the critical points, and another at 2.
The resulting $b$ distinct zeros contradict its degree at most $b-1$;
strict inequalities at the critical points exclude the zero polynomial.

Consequently $g(x-2)$ is itself a Salem trace polynomial. Its exceptional
root is $2+T(r_*)$, strictly larger than $r_*$, since
$T(r_*)=(r_*-2)H(r_*)$ with $H(r_*)>H(2)\ge1$. Every decomposition thus
supplies a Salem number of smaller degree and larger measure. This differs
from the earlier descent using the fibre at 1, which decreases measure.
