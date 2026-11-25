# A splitting condition giving a stronger Salem bound

Let $\lambda$ have Salem degree $n$, put $L=\log\lambda$, and let
$F=\mathbf Q(r)$ with $r=\lambda+\lambda^{-1}$. Write $N_f(p)$ for the
number of prime ideals of the maximal order of $F$ of norm $q=p^f$.
For every integer $k\ge2$,
$$
 fN_f(p)\log p\ \le\ \frac{(n-2)\log k}{k-1}
                     +\frac{(k+1)qL}{3}. \tag{1}
$$
No conditions on the Galois group, ramification, power-basis index, or
decomposition of the trace polynomial enter this estimate.

For fixed $p,f$ and $N_f(p)\ge c\sqrt n$, (1) gives
$L\gg_{p,f,c}1/\log n$. This improves the order of the Dobrowolski bound
on that subclass by a factor of order $(\log n)^2/(\log\log n)^3$.
More generally put $A=fN_f(p)\log p>0$. Taking
$k=1+\lceil4n\log(2n)/A\rceil$ bounds the first term by $A/2$ for all
sufficiently large $n$ at fixed $p,f$ (use $\log k\le2\log(2n)$).
Since $A=O_p(n)$ this gives $L\gg A^2/(qn\log n)$.
For example, any sequence with $L=O((\log\log n/\log n)^3)$ must satisfy
$$
 N_f(p)=O_{p,f}\left(\frac{\sqrt n\,(\log\log n)^{3/2}}{\log n}\right)
$$
at each fixed prime power. For local counts proportional to $n$, a
sufficiently large *fixed* $k$ in (1) gives $L\gg n$.

## Proof

Set $R_m=|N_{\mathbf Q(\lambda)/\mathbf Q}(\lambda^m-1)|$. Divisibility of
$X^{jm}-1$ by $X^m-1$ gives $R_m\mid R_{jm}$. For a unit-circle conjugate
$z$, apply Hadamard to the $k$-row Vandermonde matrix on the nodes
$1,z^m,\ldots,z^{(k-1)m}$. Its determinant has modulus
$\prod_{j=1}^{k-1}|1-z^{jm}|^{k-j}\le k^{k/2}$.
The two exceptional conjugates together contribute
$(\lambda^{jm}-1)(1-\lambda^{-jm})\le\lambda^{jm}$. Multiply over all
embeddings. Since $\sum j(k-j)=k(k^2-1)/6$, dividing logarithms by
$k(k-1)/2$ gives
$$
 \log R_m\le (n-2)\log k/(k-1)+(k+1)mL/3.
$$

Define $C_q(z+z^{-1})=z^q+z^{-q}$. In characteristic $p$ one has
$C_q(X)=X^q$. Thus $C_q(r)-r$ lies in every prime ideal counted by
$N_f(p)$, giving $p^{fN_f(p)}\mid |N_F(C_q(r)-r)|$. The factorization
$$
 C_q(r)-r=\lambda^{-q}(\lambda^{q-1}-1)(\lambda^{q+1}-1)
$$
and the quadratic extension from the trace field to the Salem field give
$|N_F(C_q(r)-r)|^2=R_{q-1}R_{q+1}$. Averaging the two bounds for
$\log R_m$ gives (1). The same argument can count all prime ideals with
residue degree dividing $f$, weighting each by its actual residue degree.

## Finite checks and scope

If $p=2$, $f=1$, and $N_1(2)\ge\sqrt n$, the certified lower bounds for
$L$ are 0.02419 at $n=4096$ ($k=1537$) and 0.01603 at $n=1048576$
($k=40961$), rounded down. These are conditional benchmarks. The script
also tests divisibility of exact norms and the archimedean bound on
Lehmer's example and the quartic-inner degree-16 example, using prime ideals
of the maximal orders. Run `python scripts/salem_local_norms.py`.

For a field that actually occurs in the discriminant census, take the
degree-86 ordinary Hilbert class field of the real quadratic field of
discriminant 14401. The prime 11 is inert in the base and the resulting
principal prime splits over it. The Hilbert class field therefore has 43
primes of norm 121. Any Salem number with this exact trace field has
$n=172$, and $k=2$ in (1) gives
$$
 \lambda \ge \exp((86\log 11-170\log 2)/121)>2.07600.
$$
Thus this field remains eligible for a polynomial root discriminant of
200 but already excludes small Salem numbers. The discriminant budget is
equivalent to an absolute power-basis index at most 3,458,181,874.

There is no general lower bound of order $\sqrt n$ for these local counts
in this argument. Fields with few small-norm prime ideals remain relevant
to the unrestricted Salem problem.
