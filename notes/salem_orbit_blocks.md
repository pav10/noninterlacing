# An asymptotic limit of grouping Salem conjugates

Write $n$ for the Salem degree, $L=\log\lambda$, and $K=\log n/\log\log n$. In the classical Dobrowolski parameter regime, the logarithmic saving from grouping **all chosen powers of each conjugate together with all the powers of its inverse**, after division by the coefficient of $L$ in the determinant bound, is
$$O(K^{-4})=o(K^{-3})$$
uniformly for $1/n\le L\le1$. This estimates even the exact joint Gram volumes. Voutier’s lower bound puts all sufficiently high-degree Salem numbers with $L\le1$ in this range.

Thus the finite-degree gain in [the earlier determinant calculation](salem_resultant_bounds.md) cannot change the asymptotic order by this grouping at these parameters. The estimate leaves open correlations **between different inverse pairs**, and improvements to the arithmetic lower bound. It gives no stronger general lower bound on Salem numbers by itself.

## A finite comparison

Use positive integer exponents $m$ of multiplicities $u_m$, and put
$$
 U=\sum u_m,\quad W=\sum u_m^2,\quad T=\sum m u_m,\quad
 N=nU,\quad A=\sum m u_m(N-u_m).
$$
The confluent evaluation rows have entries $\binom jr x^j$, for $0\le j<N$ and $0\le r<u_m$ at a node $x=\alpha_i^m$. Let $H_m$ be the Gram determinant of the $u_m$ rows at $x=1$. Define $B_0=\exp(AL)\prod H_m^{n/2}$, the coarse archimedean bound. Let $B_1$ be the product of the **exact Euclidean row-block volumes**, one block per inverse pair, including the exceptional real pair. Every block contains $2U$ rows. Block Hadamard gives $|\det V|\le B_1\le B_0$.

For every positive integer $h\le n/2$, the saving satisfies
$$
0\le\log(B_0/B_1)
\le \frac{n-2}{2} W\log(N/h)+W\log N
 +2U^2\log\frac1{1-e^{-L}}+(2h+2)UTL. \tag{1}
$$

**The unit-circle blocks.** Take a representative $z$ of one nonreal inverse pair. Index its nodes by exponents $e\in\{\pm m\}$ of weights $w_e=u_{|e|}$. From the $N$ columns retain the $h$ disjoint sets $j=l+hk$, $0\le k<2U$, $0\le l<h$. On each set the modulus of the confluent determinant is
$$
 h^{W-U}\prod_{a<b}|1-z^{h(b-a)}|^{w_aw_b}.
$$
Indeed $\binom{l+hk}{r}$ has leading coefficient $h^r/r!$ in $k$, so triangular row operations reduce to the usual confluent Vandermonde at the nodes $z^{eh}$. All monomial factors have modulus one. Distinct exponents give distinct nodes since a Salem conjugate has infinite multiplicative order.

The Gram matrix on the retained columns is the sum of $h$ positive definite Gram matrices of size $2U$. Minkowski's determinant inequality supplies an additional factor $h^U$ in the volume. Thus the volume has lower bound $h^W$ times the displayed product of differences. For each positive integer $v$,
$$
 \prod_{\text{nonreal inverse pairs}}|1-z^v|
 =\left(\frac{|\mathrm{Res}(P,X^v-1)|}
 {(\lambda^v-1)(1-\lambda^{-v})}\right)^{1/2}
 \ge e^{-vL/2},
$$
where $P$ is the Salem minimal polynomial. This norm estimate pays for small spacings **across all the blocks at once**. Since $\sum_{a<b}(b-a)w_aw_b\le4UT$, the product of nonreal block volumes is at least $h^{(n-2)W/2}\exp(-2hUTL)$.

**The exceptional block.** Retain the first $U$ and the last $U$ columns. The generalized Vandermonde minor at the positive real nodes $x_e=\exp(eL)$ factors into the ordinary confluent Vandermonde and the Schur polynomial with rectangular partition $((N-2U)^U)$. This follows by dividing the alternant by the ordinary Vandermonde and coalescing the repeated nodes. Row scaling contributes $\prod x_e^{\binom{w_e}{2}}=1$. By positivity of the Schur coefficients, retaining the monomial supported on the largest $U$ variables gives a lower bound $\exp((N-2U)TL)$ for the Schur factor.

In the difference product write $x_b-x_a=\exp((a+b)L/2)\,2\sinh((b-a)L/2)$. The first factors cancel in pairs under $(a,b)\mapsto(-b,-a)$. Each remaining factor exceeds $1-e^{-L}$, and there are at most $2U^2$ factors counted with multiplicity. The exceptional volume is therefore at least $\exp((N-2U)TL)(1-e^{-L})^{2U^2}$.

Finally $H_m\le N^{u_m^2}$ by row Hadamard, since row $r$ has squared norm at most $N^{2r+1}$. Also $A\le NT$. Combining the estimates gives (1).

## At the Dobrowolski scale

Take $u_1\asymp K$ and the first $s\asymp K^2$ primes of multiplicity one. Then $U,W\asymp K^2$, $T\asymp K^4\log K$ by the prime number theorem, and $A\asymp nK^6\log K$. Choose $h=\lfloor n/(\log n)^5\rfloor$ for sufficiently large $n$. Then $\log(N/h)=O(\log\log n)$ and $\log(1/(1-e^{-L}))=O(\log n)$ throughout the stated range. Division of (1) by $A$ gives $O(K^{-4})$. The divisor $A$ matters: a large change in a raw determinant can yield a negligible change in the lower bound for $L$.

For fixed arithmetic input $\log|\det V|\ge D_{\rm arith}$, the coarse inequality reads $L\ge(D_{\rm arith}-\tfrac n2\sum\log H_m)/A$. Using exact orbit blocks adds $\log(B_0/B_1)/A$ on the right. Equation (1) bounds precisely that correction. It does not assume a sequence of Salem numbers approaching one exists.

## Grouping different conjugates in a Galois-invariant partition

The argument also limits a broader use of arithmetic structure. Partition the $n$ base conjugates into blocks permuted by the Galois group, with each block closed under inversion. The block size $b$ is constant by transitivity. In each block take **every** selected power and derivative row of its base conjugates. If $b=n^{o(1)}$, their exact block volumes can save only $o(K^{-3})$ after normalization at the parameters above. For $b$ bounded by a fixed power of $\log n$, the bound is again $O(K^{-4})$. This includes partitions arising by grouping embeddings above a subfield when the inverse-pair condition holds.

Here are details for the extra issue caused by unit and real rows in the same block. Take $1\le h\le n/b$. In every block first handle only its unit-circle rows by the stride argument. Multiplying the resulting minors gives the prefactor $h^{(n-2)W/2}$. For the differences at stride $h$, temporarily restore all the omitted real nodes. The product of **squared** differences within all blocks, with the indicated multiplicities and over distinct nodes, is a nonzero rational algebraic integer by invariance of the partition. Distinct nodes at stride $h$ follow from irreducibility and distinct largest moduli for the power polynomials. Removing the differences incident on real nodes loses at most a factor
$$2^{2bU^2}\exp(hbUTL).$$
Indeed there are at most $2bU^2$ such differences with multiplicity. Only the positive real powers of $\lambda$ have modulus above one, and their weighted log moduli at stride $h$ sum to $hTL$.

There is one block containing the two real conjugates. Write $c=(b-2)U$ and let $H(X)$ be the monic degree-$c$ polynomial vanishing at that block’s unit nodes with their row multiplicities. Column vectors of polynomials $H(X)q(X)$, $\deg q<N-c$, belong to the kernel of the unit-row evaluation. Multiplication by $H$ on coefficient vectors has operator norm at most $\|H\|_1\le2^c$. Thus the real rows after projection off the unit rows have volume at least $2^{-2Uc}$ times their volume as evaluation functionals on $Hq$. To see this, the multiplier operator $T_H$ satisfies $T_HT_H^*\le 2^{2c}\Pi$ for the projection $\Pi$ onto that kernel. Take the determinant of the Gram matrices after applying the $2U$ real rows on both sides.

The jet product rule is triangular with diagonal entries $H(\lambda^{\pm m})$. These contribute a factor of at least $(1-e^{-L})^{2Uc}$, since each unit node is at least $1-e^{-L}$ from every real node. The remaining real determinant has $N-c$ columns, so the Schur minor above supplies at least $\exp((N-c-2U)TL)(1-e^{-L})^{2U^2}$. Multiplying by the unit-row volumes now bounds the whole block product $B_b$ from below. In particular
$$
\log(B_0/B_b)\le
 \frac{n-2}{2}W\log(N/h)+W\log N+
 b(h+1)UTL+4bU^2\log 2+2bU^2\log\frac1{1-e^{-L}} .
$$
Choose $h=\lfloor n/[b(\log n)^5]\rfloor$. For $b=n^{o(1)}$ the first logarithm is $o(\log n)$; for $b$ polylogarithmic it is $O(\log\log n)$. Dividing by $A$ proves both assertions. Partitions selected by angular proximity need not be Galois invariant. The integer-product step is precisely where this distinction enters.

As a finite implementation check, *python scripts/salem_orbit_check.py* computes the Gram blocks with Arb and the full determinant from exact discriminants and resultants. It checks (1) on Lehmer's degree-ten polynomial and on degree-64 and degree-122 Salem–Boyd examples. For multiplicities $(u_1,u_2,u_3)=(2,1,1)$, the normalized within-pair savings at degrees 64 and 122 were 0.02115 and 0.01220, while the remaining cross-pair savings were 0.14120 and 0.15711. These examples have Salem numbers near 1.325. Their observed savings do not establish comparable savings for hypothetical numbers near 1.
