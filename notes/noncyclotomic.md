# Two noncyclotomic constructions

Here $D$ is the degree of a monic irreducible totally real polynomial $P$ and $\mathrm{rd}(P)=|\mathrm{disc}(P)|^{1/D}$. In the Salem computations $P$ is the totally real trace polynomial. Its reciprocal lift has degree $2D$ and Mahler measure equal to the Salem number.

## Non-Galois fields below the polynomial discriminant threshold

Let $n_0=765765=3^2\cdot5\cdot7\cdot11\cdot13\cdot17$, $q=9189181$ and $n=n_0q^j$ for $j\ge0$. Put $\beta=\zeta_n+\zeta_n^{-1}$, let $\Psi_n$ be its minimal polynomial, and write $m=\varphi(n)/2$. Define

$$P_n(X)=X^m\Psi_n(X-X^{-1}).$$

This is an integer polynomial with constant term $\pm1$. Its quadratic factors over the real cyclotomic field $K=\mathbf Q(\beta)$ have the form $X^2-\beta_iX-1$. All roots are real. We will see that $P_n$ is irreducible, each root generates a non-Galois extension of $\mathbf Q$, and

$$\mathrm{rd}(P_n)<0.995207\,D,\qquad D=2m.$$

The first degree is 276480. This strict bound also certifies non-interlacing by the residue product argument in [the discriminant note](high_degree_discriminants.md).

**Irreducibility and nonnormality.** The prime $q$ splits completely at conductor $n_0$. Modulo $q$, the element $w=5221919$ has exact order $n_0$. At the primes corresponding to $w$ and $w^2$, the residues of $\beta$ are 4044811 and 4659871. The two values of $\beta^2+4$ have quadratic characters $+1$ and $-1$, respectively. Raising the conductor to $n_0 q^j$ preserves these residue choices. Locally the $q$-power cyclotomic extension is totally ramified, and the $q$-power roots of unity reduce to 1. Thus the quadratic extension $K(\sqrt{\beta^2+4})/K$ splits at one of the two primes and is inert at the other. In particular its degree is two. It cannot be Galois over $\mathbf Q$, since conjugate primes of $K$ have different relative splitting types. A root $t$ of the quadratic generates this extension because $\beta=t-1/t$.

**Polynomial discriminant.** The absolute resultant of two distinct quadratic factors is $(\beta_i-\beta_j)^2$, by evaluating their difference at the two roots of either one. Hence

$$|\mathrm{disc}(P_n)|=\mathrm{disc}(\Psi_n)^2\,N_{K/\mathbf Q}(\beta^2+4).$$

For $u=\zeta_n^{2a}$ with $a$ a unit modulo $n$, $\beta_a^2+4=6+u+u^{-1}$. Since $9\mid n$, the primitive roots partition into triples $u,\omega u,\omega^2u$ for a primitive cube root $\omega$. The three factors have product $198+u^3+u^{-3}\le200$. Running over primitive roots counts each real embedding twice. Consequently $N(\beta^2+4)^{1/(2m)}\le200^{1/6}$. The real cyclotomic discriminant formula gives

$$
\frac{\mathrm{rd}(P_n)}{D}
\le 200^{1/6}\prod_{p\mid n}\frac{p^{1-1/(p-1)}}{p-1}.
$$

At $n_0$ the bound is approximately 0.995206149. Raising it to the 240th power gives a rational number strictly below 1 by exact integer comparison. For $j>0$ the additional factor for $q$ is less than 1, since $(1+1/(q-1))^{q-1}<e<q$.

**Comparison with Simon’s questions.** Interpreting the displayed targets as $(D/4)^{D/2}$, $D^{D/2}$ and $(D/\log D)^D$, even the last asks for $\mathrm{rd}(P)<D/\log D$. In this family $\beta_i^2+4\ge4$ gives $\mathrm{rd}(P_n)/D\ge\mathrm{rd}(\Psi_n)/m$, which exceeds 0.82 for every $j$. Thus this family isolates the gap between the non-interlacing threshold and Simon’s targets. For a bounded-root-discriminant field tower, the unresolved generator problem is to bound $[\mathcal O_K:\mathbf Z[t]]^{1/[K:\mathbf Q]}$ uniformly. The calculation above works directly with the polynomial discriminant.

## Salem traces with a proof of non-interlacing

Let $F=\mathbf Q(\zeta_{21}+\zeta_{21}^{-1})$, with generator $\beta$ and ring of integers $\mathbf Z[\beta]$. This ring has no totally positive codifferent element of trace one. One can check this in exact arithmetic: every monic interlacer of $\Psi_{21}$ would have degree five and all roots strictly between $-2$ and $2$. Kronecker’s theorem restricts it to a squarefree product of real cyclotomic polynomials. All 154 degree-five products, including candidates at the endpoints, fail the exact interlacing test.

If $F\subset E=\mathbf Q(t)$ and $\mathcal O_F\subseteq\mathbf Z[t]$, a trace-one totally positive element $\mu$ of the codifferent of $\mathbf Z[t]$ would give $\mathrm{Tr}_{E/F}(\mu)$ with the same properties in the codifferent of $\mathcal O_F$. Integrality against each element of $\mathcal O_F$ follows from transitivity of trace. Thus an order inclusion suffices to certify non-interlacing of the minimal polynomial of $t$.

**A small-degree benchmark.** Exactly one trace polynomial with one root in $(2,4]$ and the other five in $(-2,2)$ generates $\mathcal O_F$ as a power order:

$$x^6-3x^5-5x^4+8x^3+4x^2-5x+1.$$

Its Salem number is $3.420749775618197\ldots$, the minimum among Salem traces generating this particular ring. To verify completeness, use the integral basis $1,C_1(\beta),\ldots,C_5(\beta)$ with $C_j(2\cos t)=2\cos(jt)$. For coefficient vector $c$ and trace Gram matrix $G$, the root bounds give $c^tGc<36$. Dual Cauchy–Schwarz bounds $c_i^2$ by $36(G^{-1})_{ii}$. The integer coordinate bounds are $(2,1,2,2,1,2)$. Enumerating this box, testing the characteristic polynomial discriminant for index one, and counting roots by Sturm gives the polynomial above.

**Unbounded degree.** Take $A=1/(2+\beta)-1$. Both $2+\beta$ and $A+1$ are units with integral inverses, so $\mathbf Z[A]=\mathcal O_F$. Its minimal polynomial is

$$h(Y)=Y^6-42Y^5-77Y^4-14Y^3+49Y^2+35Y+7.$$

Five conjugates lie in $(-1,1)$ and the sixth is $A_* =43.76606865\ldots$ in $(43,44)$. These counts follow by Sturm. Define $V_0=1$, $V_1=x-1$ and $V_k=xV_{k-1}-V_{k-2}$. Factor the integer norm polynomial

$$H_k(x)=V_{k-1}(x)^6 h\big(V_k(x)/V_{k-1}(x)\big).$$

For any real conjugate $a$ of $A$, the polynomial $V_k-aV_{k-1}$ has simple real roots. For $k\ge2$ it is the characteristic polynomial of a tridiagonal matrix with unit off-diagonals and diagonal $(1,0,\ldots,0,a)$. The roots interlace the zeros of $V_{k-1}$, which are in $(-2,2)$ by the cosine formula. At $2$ the value is $1-a$, and at $-2$ its sign-adjusted value is $(2k+1)+a(2k-1)$. It follows that the five conjugates with $|a|<1$ give only interior roots, and $A_*$ gives exactly one root above 2 and all others interior. The same conclusion for $k=1$ follows from $V_1-a=x-1-a$.

For the exceptional root $t$, the recurrence gives $(V_k,V_{k-1})=(1)$ in $\mathbf Z[x]$. Substituting $V_k(t)=A_*V_{k-1}(t)$ first shows that $V_{k-1}(t)$ is a unit among the algebraic integers. Its inverse is then in $\mathbf Z[t]$ by its monic minimal polynomial with constant term $\pm1$. Hence $A_*\in\mathbf Z[t]$, and the trace-down argument certifies the irreducible factor containing $t$. The corresponding reciprocal lift is irreducible since the trace field is totally real and at an interior embedding the quadratic $z^2-tz+1$ has nonreal roots.

| $k$ | Certified non-interlacing trace degree | Salem degree |
|---:|---:|---:|
| 1 | 6 | 12 |
| 4 | 24 | 48 |
| 12 | 72 | 144 |
| 32 | 192 | 384 |
| 64 | 384 | 768 |

The Salem numbers satisfy $\lambda^{2k}(\lambda-A_*)=A_*\lambda-1$ and converge to $A_*$ from above. Indeed $\lambda>A_*$ and $\lambda-A_*<A_*/\lambda^{2k-1}$. Distinct $k$ give distinct values, and bounded-degree monic integer polynomials with bounded roots form a finite set. This proves unbounded degree even without irreducibility of every full $H_k$. At $k=64$ the trace polynomial has root discriminant about $9.38$ times its degree. Its non-interlacing certificate comes from trace-down.

## A separate search near 1

Factor $C_d-C_a-C_b+1$ for $0\le b<a<d$, $\gcd(a,b,d)=1$, with $C_0=2$. The sampled degrees were 12, 16, 20, 24, 30, 32, 36, 40, 48, 60, 64, 80 and 96. Among 10845 distinct irreducible factors of degree at least two with negative value at 2, exactly three had Salem measure below $5/4$. Their measures were 1.176280818260, 1.202616743689 and 1.240726423653, with trace degrees 5, 7 and 6. For each potential factor $g$, the tests $g(2)<0<g(41/20)$ and an exact count of $\deg(g)-1$ roots in $(-2,2)$ decide the target condition, since $41/20=5/4+4/5$. This count covers the specified parameter set only.

All three small trace polynomials admit monic interlacers. In order of increasing measure, exact witnesses are $x^4+x^3-3x^2-3x$, $x^6-6x^4+9x^2-1$ and $x^5-x^4-3x^3+3x^2+x-1$.

Run the algebraic certificates with `python scripts/noncyclotomic_certificates.py` and the sparse scan with `python scripts/sparse_salem_traces.py`. Both use the repository’s optional PARI dependency for exact root counts.
