# Fourier modes, shifted relations and the limits of the correlation bound

Let $n$ be the Salem degree, $m=n-2$, $L=\log\lambda$, and
$\mu_j=m^{-1}\sum_z z^j$, summed over the unit-circle conjugates.
Write $r=1-1/n$ and $A_n=\log n-\tfrac{m-1}{2}\log r+2\log4$.
The [power-sum estimate](salem_factor_relations.md#a-test-for-witnesses-with-unboundedly-many-terms)
applied to $\lambda^j$ gives
$$ rm|\mu_j|^2\le A_n+2(1+1/m)jL. \tag{1} $$

## Using shifts of one relation

Let $H$ be any nonzero polynomial multiple of the Salem polynomial, with
exponent support $S$. Coefficient sizes and other factors are unrestricted.
Choose a nonempty finite set $D$ of nonnegative integer shifts. Put
$K=S+D$, $k=|K|$, $h=|D|$, and let $b$ be the average of $|u-v|$ over
distinct ordered pairs $u,v\in K$. Then
$$
 2(1+1/m)bL\ \ge\ rm\frac{h}{(k-h)(k-1)}-A_n. \tag{2}
$$
This only needs the support of the witness. To prove it, form the Gram
matrix $G=(\mu_{u-v})_{u,v\in K}$ of the monomials evaluated on the
ordinary conjugates with normalized counting measure. Each $x^jH$ for
$j\in D$ supplies a null vector. These polynomials are linearly
independent since multiplication by $H$ is injective. The trace is $k$
and the rank is at most $k-h$. Cauchy–Schwarz on the nonzero eigenvalues
gives $\operatorname{tr}(G^2)\ge k^2/(k-h)$. On removing the $k$ diagonal
entries, the average squared off-diagonal entry is at least
$h/((k-h)(k-1))$. Average (1) over these entries to get (2).

**Two blocks with growing width.** Suppose a multiple has the form
$B(x)+x^M C(x)$, with $\deg B,\deg C\le s$, $s\ge1$. For
$D=\{0,\ldots,s\}$ we have $h=s+1$, $k\le4s+2$ and $b\le M+2s$.
In particular the factor in (2) is at least
$(s+1)/((3s+1)(4s+1))$. If $M+s=O(n)$ and $s=o(n/\log n)$,
then $L\gg1/s$ as $n\to\infty$, with the constant depending on the
degree-ratio bound. This improves the *order* of Dobrowolski when
$s=o((\log n/\log\log n)^3)$. The block coefficients can vary with $n$,
and this proof imposes no reciprocity condition on the witness.

More generally, if $|S+S|\le C|S|$, take $D=S$. For witnesses of degree
$O(n)$ with $C$ fixed this gives the same order $1/|S|$ when
$|S|=o(n/\log n)$. In the other direction, if the minimal polynomial
has no missing coefficients, then $S=\{0,\ldots,n\}$ and any $D$ gives
$k\ge n+h$: in increasing order of the shifts the first interval
contributes $n+1$ integers and each subsequent interval at least one.
Thus $h/((k-h)(k-1))\le1/n$, and (2) gives no positive bound. This
identifies where coefficient information must enter a support-based
approach.

## A witness-free upper bound on the method

For a witness $H$ use the normalized absolute autocorrelation sum $J$
and its mean lag $B$ as in the [preceding note](salem_factor_relations.md#a-test-for-witnesses-with-unboundedly-many-terms).
Define $U=\sup_{j\ge1}|\mu_j|^2/j$. Every witness satisfies
$$ J^2 B\ge 1/U. $$
Indeed $1/J^2\le\sum w_j|\mu_j|^2$ for weights $w_j$ proportional
to the absolute autocorrelations. The right side is at most $UB$.
Likewise the lower bound on $L$ from any one witness in that note, or
from (2), is at most
$$ V=\sup_{j\ge1}\frac{rm|\mu_j|^2-A_n}{2(1+1/m)j}\ \le L. $$
This follows by averaging the inequalities defining $V$. Thus the best
single-mode test dominates these averaged tests for a *known* Salem
polynomial. A witness or support argument is useful when it forces a
large mode across an entire family of unknown polynomials.

In particular, fix $T>1$. A sufficient missing lemma for a uniform gap
would be the existence of constants $c,K>0$ such that every Salem number
below $T$ of sufficiently large degree $n$ has some $1\le j\le Kn$
with $|\mu_j|\ge c$. Equation (1) would then give
$\liminf L\ge c^2/(2K)$ along any such sequence of increasing degrees.
Conversely, if $\lambda\to1$, bounded-degree finiteness forces
$n\to\infty$, and (1) gives $\max_{1\le j\le Kn}|\mu_j|\to0$ for
every fixed $K$. These implications concern a gap at 1. Establishing
attainment of the infimum also requires excluding accumulation from
above at a positive infimum that is not itself a Salem number.

Both quantities have tail bounds independent of the polynomial:
beyond mode $j$ the raw maximum is at most $1/(j+1)$, and the adjusted
maximum is at most $\max(0,rm-A_n)/(2(1+1/m)(j+1))$. The implementation
uses Arb to compare all the modes up to the stopping point and certify
the unique maximizer. If $rm\le A_n$, the adjusted supremum is zero.
Otherwise simultaneous approximation of the finitely many angles
gives a positive adjusted value eventually.

For $x^t(x^3-x-1)-(1-x^2-x^3)$, select the irreducible factor containing
the Salem number. The Rouché argument for the [Salem–Boyd family](salem_minimum.md)
establishes this signature; the script verifies it by Sturm at each
sampled degree. All displayed lower bounds are rounded down.

| $t$ | Salem degree | Raw maximizing mode | Adjusted maximizing mode | $V$ exceeds |
|---:|---:|---:|---:|---:|
| 8 | 10 | 1 | 261 | 0.00180 |
| 31 | 30 | 31 | 61 | 0.01127 |
| 63 | 64 | 63 | 63 | 0.04395 |
| 127 | 122 | 127 | 127 | 0.06887 |
| 255 | 256 | 255 | 255 | 0.08584 |
| 511 | 510 | 511 | 511 | 0.09349 |

The large mode at $j=t$ has an explanation. With $A=x^3-x-1$ and
$R=A^*/A$, the ordinary conjugates satisfy $z^t=R(z)$. The ordinary
conjugates become equidistributed as the degrees tend to infinity
(apply the fixed-mode estimate (1) and the Fourier criterion).
These degrees do tend to infinity since the Salem numbers are distinct
and bounded. $R$ is continuous on the unit circle. Therefore $\mu_t$
tends to its circle mean, which is
$-1-A^*(\theta)/(\theta A'(\theta))=-0.45488440437\ldots$, where
$\theta$ is the real root of $A$. Integrate $R(z)/z$ and subtract its
residue at the unique exterior pole from $R(\infty)=-1$.
Run `python scripts/salem_fourier_modes.py` for the finite certificates.

## Sharp obstruction from cyclotomic orbits

First replace the ordinary nodes by all primitive $q$th roots and allow
any multiple of $\Phi_q$. Put $M=q/\gcd(q,2)$. In this model
$$ \min_H J(H)^2 B(H)=M. $$
Indeed the normalized Ramanujan sum is
$\mu_j=\mu_{\rm Mob}(v)/\varphi(v)$, $v=q/\gcd(j,q)$.
For squarefree $v$, $\varphi(v)^2\ge v/\gcd(2,v)$ by multiplying
$(p-1)^2\ge p$ for odd primes. Non-squarefree $v$ gives a zero sum.
Hence $|\mu_j|^2\le j/M$. The lower bound follows from the witness-free
inequality. Equality is realized by $x^q-1$ for odd $q$ and by
$x^{q/2}+1$ for even $q$: the only nonzero positive-lag autocorrelation
is at $M$, with $J=1$.

### Transfer to actual Salem numbers

Fix $q\ge3$, let $f=\Psi_q$ have degree $d=\varphi(q)/2$, and set
$$ h_t(y)=(y-t)f(y)-f'(y),\qquad t\in\mathbf Z. $$
For sufficiently large $t$, the reciprocal lift $P_t$ is Salem of
degree $n=2d+2$. In fact $h_t/f=y-t-\sum_i1/(y-r_i)$ is strictly
increasing between the roots $r_i$ of $f$. For large $t$ the leftmost
root is above $-2$ (check $h_t(-2)/f(-2)<0$); there is one root
$s_i$ immediately to the left of each pole, and one $y_*>t>2$.
A proper factor away from $y_*$ would be real cyclotomic by Kronecker.
At bounded degree there are only finitely many possible cyclotomic
roots $u$, and $h_t(u)=0$ fixes at most one $t$ since $f$ and $f'$ are
coprime. Thus $h_t$ is irreducible for all sufficiently large $t$.
One interior embedding then gives irreducibility of its reciprocal lift.

The trace identity gives
$\sum_i(r_i-s_i)=y_*-t=\sum_i1/(y_*-r_i)\le d/(t-2)$.
Since $|C_j'|\le j^2$ on $[-2,2]$ for $C_j(2\cos u)=2\cos(ju)$,
the difference of the normalized ordinary and cyclotomic moments
is at most $j^2/(2(t-2))$. For $j\le M$ combine this with the Ramanujan
bound; for $j>M$ use $|\mu_j|\le1$. The result for *every* multiple is
$$
 J^2 B\ge M\big(1+M^2/(2(t-2))\big)^{-2}.
$$
This lower bound tends to the exact optimum as $t\to\infty$ with $q$
fixed. Indeed multiplying $P_t$ by the cyclotomic quotient of the
extremizing binomial by $\Phi_q$ gives, after division by $-t$, a
polynomial converging coefficientwise to $x$ times that binomial.
The cost is continuous there and invariant under scaling and shifting.

For primorial $q$, $M/n=q/(2(\varphi(q)+2))$ grows without bound by
the Euler product. Taking $t$ sufficiently large along these conductors
therefore gives genuine Salem polynomials for which **every** choice of
multiple has $J^2B/n$ arbitrarily large. The Salem numbers in this
construction grow with $t$; it leaves open uniform control close to 1.

| Conductor | Salem degree | Certified interval for the optimal cost $\inf_H J^2B/n$, $t=10^9$ |
|---:|---:|---:|
| 30 | 10 | [1.499999, 1.500001] |
| 210 | 50 | [2.099976, 2.100016] |
| 2310 | 482 | [2.393072, 2.398458] |

Run `python scripts/salem_cyclotomic_correlation.py` for the exact
irreducibility checks, rational lower bounds and explicit integer
witnesses giving the upper bounds.
