# The five-term minimum with arbitrary integer coefficients

Let $P$ be the monic minimal polynomial of a Salem number $\lambda$,
multiplied by any product of cyclotomic polynomials. If $P$ has at most
five nonzero terms and $\lambda\le5/4$, there are exactly six possible
Salem numbers. All the witnesses in this range have height one. In
particular **Lehmer's number is the minimum with no coefficient bound
and no degree bound** in this class.

| Salem number (rounded) | Salem degree | Exponent pairs $(a,b)$ |
|---:|---:|---|
| 1.176280818260 | 10 | (6,1), (7,4) |
| 1.200026523987 | 14 | (10,9) |
| 1.202616743689 | 14 | (7,5) |
| 1.216391661138 | 10 | (5,1), (9,8) |
| 1.230391434407 | 10 | (5,2) |
| 1.236317931803 | 16 | (8,7) |

The polynomial for each listed pair is
$x^{2a}-x^{a+b}-x^a-x^{a-b}+1$. The cyclotomic cofactor is essential.
Finding a Salem factor in a polynomial with other noncyclotomic factors
does not put it in this class.

## A span bound at every power

Let $2d$ be the Salem degree, $d\ge3$, and fix $b\ge1$. The ordinary
conjugates of $y=\lambda^b+\lambda^{-b}$ are the $d-1$ conjugates other
than $y$ itself. They are distinct and lie in $[-2,2]$. Indeed the trace
field of $\lambda$ contains $y$, and only one of its embeddings maps
$y$ to the exceptional value above two. Thus $y$ still generates the
full trace field. If its ordinary conjugates lie in an interval of
length $W$, then
$$
 2b\log\lambda\ \ge\ (d-2)\log(4/W)-\log(64(d-1)). \tag{1}
$$

To prove this, put $m=d-1$. In the Vandermonde determinant on the
ordinary roots, use monic Chebyshev polynomials scaled to their
interval. The $j$th polynomial has sup norm $2(W/4)^j$ for $j\ge1$,
and the constant column is 1. By Hadamard the squared determinant is
at most $m^m4^{m-1}(W/4)^{m(m-1)}$. Each difference from $y$ is at
most $\lambda^b+3\le4\lambda^b$. Restoring these differences gives an
upper bound on the nonzero integer discriminant of the minimal
polynomial of $y$. Taking the $m$th root and relaxing $4^{1-1/m}$ to
4 proves (1).

This has a consequence even for **unrestricted term count**. Write a
Salem minimal polynomial on the circle as
$x^d(c_0+\sum_{j=1}^d c_j(x^j+x^{-j}))$, $c_d=1$. If $c_b\ne0$, put
$S_b=\sum_{j\ne b,\ j\ge1}|c_j|$. At every ordinary conjugate the
relation puts $x^b+x^{-b}$ in an interval of width at most
$4S_b/|c_b|$, since the constant term just shifts the interval.
Consequently if $|c_b|\ge(1+\delta)S_b$ for a fixed $\delta>0$, then
along any sequence with $d\to\infty$ one has
$\liminf\lambda\ge\sqrt{1+\delta}$, using $b\le d$ and (1).
For example, a coefficient at least twice the sum of the other
nonconstant trace coefficients forces asymptotic lower bound $\sqrt2$.
Candidates approaching 1 must satisfy
$|c_b|\le(1+o(1))S_b$ uniformly over their nonzero modes.
This necessary condition still allows many coefficient patterns.

## Five and four terms force short intervals

The sparse coefficient bound $|c_i|\le\binom{s-1}{i}M(P)$, with
coefficients indexed in exponent order, was proved in the
[bounded-term note](salem_minimum.md). Here $M(P)=\lambda\le5/4$.
Five terms force reciprocity and the shape
$$
 P=x^{2a}+u(x^{a+b}+x^{a-b})+v x^a+1,\qquad
 0<b<a,\quad |u|\le5,\quad |v|\le7 .
$$
Both coefficients are nonzero. We require $\gcd(a,b)=1$, since
otherwise rotating the exceptional root produces a second outside
root. Also $P(1)\le0$. Height one is covered by the earlier exact
classification. For each unit root $z$ put $Y=z^b+z^{-b}$. The relation
$z^a+z^{-a}=-uY-v$ confines $Y$ to the intersection of $[-2,2]$ with
$|uY+v|\le2$. It has width at most two if $|u|\ge2$. If $|u|=1$ and
$|v|\ge2$, its width is again at most two, by the shift.

With four terms write
$$
 P=x^N+u x^{N-i}+\epsilon u x^i+\epsilon,\qquad
 0<i<N/2,\quad \gcd(N,i)=1,\quad \epsilon=\pm1,\quad |u|\le3 .
$$
For $|u|=1$ this is a product of binomials with unit roots. Otherwise
put $b=N-2i$ and square the equation obtained by dividing $P(z)=0$
by $z^{N/2}$ (either choice of square root). This gives
$$z^N+z^{-N}+2\epsilon=u^2(z^b+z^{-b}+2\epsilon).$$
Thus $Y$ lies in an interval of width at most one. Two terms are
cyclotomic. With three terms the only possible shape is
$x^{2a}+v x^a+1$; rotation excludes $a>1$ and degree two excludes a Salem
number for $a=1$.

## Bounding the contribution of cyclotomic factors

If $z$ is a cyclotomic root, every conjugate of its $Y$ lies in
the same interval of width at most two. A real cyclotomic trace with
that span is an integer. Indeed for a monic irreducible polynomial of
degree $k$ in a width-two interval, the same Vandermonde estimate gives
$|D|\le k^k4^{k-1}2^{-k(k-1)}<1$ for $k\ge5$. Check $k=5$ directly
and use $4k<2^{k-1}$ for $k\ge6$. In degrees 2 through 4, exhaust
$\varphi(h)\le8$ using $\varphi(h)\ge\sqrt{h/2}$; all nonlinear real
cyclotomic traces in this range have certified span greater than two.
Thus $Y\in\{-2,-1,0,1,2\}$ and $(z^b)^{12}=1$.

In the five-term case an interval of width two contains at most three
integers, giving at most six choices of $z^b$. Each fixes the trace
of $z^a$ and therefore gives at most two choices of $z^a$. This latter
trace is also integral, so $(z^a)^{12}=1$. Coprimality determines
$z$ and gives $z^{12}=1$. There are at most twelve distinct cyclotomic
roots. Their multiplicities are at most four: applying $(x\,d/dx)^j$
at a putative nonzero root for $0\le j<s$ would otherwise give a
singular Vandermonde on the $s$ distinct exponents. Hence $2a\le2d+48$.

With four terms the width-one interval gives at most four choices of
$z^b$ and then two of $z^N$. The latter again has integral trace.
Since $\gcd(N,b)\le2$, there are at most sixteen choices of $z$, each
of order dividing 24 and of multiplicity at most three. Hence
$N\le2d+48$. These bounds concern the entire witness, including
repeated cyclotomic factors.

## Exact cutoff and census

Use $s=1$ for five terms and $s=2$ for four in the larger-coefficient
cases. We have $W\le4/2^s$ and $b\le s(m+25)$, $m=d-1$. For a
degree-two trace ($m=1$) the upper bound on the degree of the witness
already suffices. Otherwise (1) implies
$$ 1\le64m\,2^{-s(m-1)}(5/4)^{2s(m+25)}. $$
For both $s$, this upper bound is below 1 at $m=88$, and its ratio at
successive integers is below 1 thereafter. Both comparisons use integer
arithmetic. Thus $d\le88$, and the entire witness has degree at most
224. For height-one pentanomials the earlier analytic cutoff is $a\le12$.

The script removes the gcd with $(x^{12s}-1)^4$, which removes exactly
the cyclotomic cofactor of any admissible candidate. It checks the
residual trace polynomial using exact irreducibility, Sturm counts and
the rational endpoint $41/20=5/4+4/5$. Prior filters use (1), the sign
at 2, and the necessary real-rootedness inequality $T'^2-TT''\ge0$.
The sweep considered 226265 larger-coefficient pentanomials and 61216
quadrinomials, with 1391 residual Sturm tests. None contributed a
witness in range. The independent height-one classification supplies
all eight exponent pairs in the table.

Run `python scripts/salem_five_unrestricted.py` with the project's PARI
dependency installed.
