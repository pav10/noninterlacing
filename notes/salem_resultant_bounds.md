# Salem-specific estimates in Dobrowolski’s determinant

Let $\lambda>1$ be a Salem number with minimal polynomial $P$, degree $n=2d$, trace $r=\lambda+\lambda^{-1}$ and trace field $F=\mathbf Q(r)$. All bounds below concern $L=\log\lambda$, with natural logarithms. The trace polynomial has degree $d$.

| Salem degree $n$ | Voutier benchmark | Same determinant, coarse estimate | Reciprocal-pair estimate | With abelian trace field |
|---:|---:|---:|---:|---:|
| 32 | 0.011532 | 0.032920 | 0.034785 | 1.166 |
| 128 | 0.008623 | 0.021561 | 0.021927 | 3.746 |
| 512 | 0.006318 | 0.014949 | 0.015019 | 10.361 |
| 1024 | 0.005448 | 0.012738 | 0.012769 | 16.657 |
| 4096 | 0.004130 | 0.009550 | 0.009557 | 40.307 |
| 1048576 | 0.001706 | 0.004026 | 0.004026 | 935.252 |

The Voutier benchmark is $\tfrac14(\log\log n/\log n)^3$, rounded to nearest. Each of the last three columns gives a **strict lower bound** certified at that exact degree. The coarse and paired columns use identical multiplicities, primes and discriminant input. Their difference isolates the archimedean gain. These are finite-degree comparisons. Optimality among published bounds and an improved asymptotic exponent for general Salem numbers remain unestablished.

## 1. Arithmetic lower bound for the confluent determinant

Choose exponents $m$ consisting of 1 and distinct primes, with multiplicities $u_m$. Write $P_m$ for the monic polynomial on the conjugates $\alpha_i^m$ of $\lambda^m$. No two conjugates of $\lambda$ have a root-of-unity quotient. An automorphism taking the first to $\lambda$ would take the second to a different conjugate of the same modulus. Thus $P_m$ is irreducible of degree $n$ and $\mathbf Q(\lambda^m)=\mathbf Q(\lambda)$. Different $m$ give disjoint root sets because their largest moduli differ.

For each prime $p$, $P(\lambda^p)\equiv P(\lambda)^p=0\bmod p\mathbf Z[\lambda]$. Taking norms gives $p^n\mid\operatorname{Res}(P,P_p)$. Also $|\mathrm{disc}(P_m)|\ge|\mathrm{disc}(P)|$ by the power-order index. Form the $N$-column confluent Vandermonde determinant $V$, with $N=n\sum_m u_m$, on the roots of $P_m$ repeated $u_m$ times. Use derivatives divided by the derivative factorial. Its product formula gives
$$
\log|V|\ \ge\ \tfrac12\sum_m u_m^2\log|\mathrm{disc}(P)|
       +nu_1\sum_{p}u_p\log p .
$$
The reciprocal discriminant identity gives $|\mathrm{disc}(P)|=|f(2)f(-2)|\mathrm{disc}(f)^2$ for the trace polynomial $f$. Here both endpoint values are nonzero integers. Minkowski’s bound in the totally real degree-$d$ field gives $\log|\mathrm{disc}(P)|\ge4\log(d^d/d!)$. These ingredients apply individually to every Salem number of the indicated degree.

## 2. What the two real conjugates save

The rows at a node $x$ can be taken as $\binom jr x^j$, $0\le r<u$, $0\le j<N$. This rescales normalized derivatives and preserves the absolute value of the full determinant, since the conjugate product has modulus one. Let $H(N,u)$ be the Gram determinant of these rows at $x=1$. At every unit-circle node the squared volume is exactly $H(N,u)$.

For the reciprocal nodes $e^w,e^{-w}$, $w=mL$, let $B(w)$ be the product of the two separate row-block volumes. Cauchy–Binet and reversal of column subsets give
$$
 B(w)=\sum_{|J|=u} W_J\exp\!\left(2w\left[\sum_{j\in J}j-u(N-1)/2\right]\right),
 \qquad W_J=\left(\det[\tbinom jr]_{0\le r<u,\ j\in J}\right)^2.
$$
Here the row index runs from 0 to $u-1$. Reversal preserves $W_J$ by the Vandermonde formula for the minor. This is a sum of coshes with nonnegative weights and $B(0)=H(N,u)$. With $z=wu(N-u)$ it follows that $B(w)/H(N,u)\le\cosh z$.

A second upper bound uses the decaying exponential series. Extending the subset sum to all nonnegative integers yields $B(w)\le e^z(1-e^{-2w})^{-u^2}$. To see the normalization, the infinite sum $\sum_J W_J x^{\sum J}$ equals $x^{u(u-1)/2}(1-x)^{-u^2}$ for $0<x<1$. Its Gram matrix has bivariate generating function $[1-x(1+y)(1+v)]^{-1}$. Setting $h=x/(1-x)$, this kernel expands as
$$
\frac{1}{1-x}\sum_{\ell\ge0}\left(\frac{x}{(1-x)^2}\right)^\ell
 \frac{y^\ell}{(1-hy)^{\ell+1}}\frac{v^\ell}{(1-hv)^{\ell+1}} .
$$
The two series factors give unit triangular matrices. The product of the first $u$ diagonal entries gives the claimed determinant.

The two blocks also have a definite angle. Their constant rows have normalized inner product $1/S_N(w)$, where $S_N(w)=\sinh(Nw)/(N\sinh w)$. Projecting the second block off the first reduces the joint volume by at least the corresponding sine factor. This follows by starting Gram–Schmidt on the second block with its constant row. All subsequent projection lengths can only decrease. The joint volume is at most $B(w)\sqrt{1-S_N(w)^{-2}}$.

For a proposed $L\le t$, bound each joint volume at $w=mt$. Both the exact $B(w)$ and the sine factor are nondecreasing on $w>0$ by their symmetric exponential sums. Thus the second upper bound above need only be evaluated at the right endpoint. Block Hadamard now gives the following upper bound on $\log|V|$:
$$
\frac n2\sum_m\log H(N,u_m)+
\sum_m\left[
\min\{\log\cosh z,\ z-u_m^2\log(1-e^{-2w})-\log H(N,u_m)\}
+\tfrac12\log(1-S_N(w)^{-2})
\right].
$$
Replacing each square bracket by $z$ gives the coarse comparison. All choices of nonnegative multiplicities are valid after unused exponents are omitted. The script uses fixed choices found by a numerical search and independently proves the displayed inequalities using integer Gram determinants and 192-bit Arb intervals.

## 3. A Frobenius obstruction, including nonabelian splitting fields

Let $p$ be unramified in the Galois closure of $F$. Its Frobenius conjugacy class acts on the $d$ trace conjugates $r_i$. For one index $i$, let $b$ be the number of distinct images of $i$ under that class. Transitivity makes $b$ independent of $i$. Form all edges $(i,j)$ from these images. Every vertex has $b$ incoming and $b$ outgoing edges. Then
$$
 L\ \ge\ \frac{d}{p+1}\left(\frac{\log p}{b}-\log4\right).
$$
Indeed, multiply $C_p(r_i)-r_j$ over the edges, where $C_p(z+z^{-1})=z^p+z^{-p}$. The product is a rational integer by Galois invariance. None of its factors vanish, since conjugation to the exceptional trace would give a value larger than every $r_j$. At a chosen prime above $p$, the Frobenius permutation supplies $d$ edges vanishing modulo that prime because $C_p(x)\equiv x^p\bmod p$. Unramifiedness gives divisibility of the integer product by $p^d$. Its archimedean absolute value is at most $4^{bd}\lambda^{(p+1)b}$: only the exceptional vertex contributes powers of $\lambda$, with $b$ occurrences in each direction. This proves the estimate. For a completely split prime $b=1$ even if the trace field is non-Galois.

If $F$ is abelian, every unramified prime has $b=1$. Hence a proposed $L\le t$ forces **every** prime satisfying $d\log(p/4)>(p+1)t$ to ramify. In a Galois field of degree $d$, each ramified prime contributes at least $d/2$ to its discriminant exponent by the different bound. Thus their product is at most $\mathrm{rd}(F)^2$. Ordinary Vandermonde–Hadamard on $P$ gives $|\mathrm{disc}(P)|\le n^n\lambda^{2(n-1)}$, and the reciprocal discriminant identity yields $\mathrm{rd}(F)^2\le n^2 e^{4t}$. A product of forced primes exceeding that number excludes $t$. Primes through 10000 suffice for the abelian column.

For fixed $\lambda\le B$ this obstruction eliminates abelian trace fields in sufficiently high degree. It is a computable local criterion also applicable when individual Frobenius orbits in a nonabelian example are known. For context, the relative height theorem of Amoroso and Zannier already implies an asymptotically stronger linear lower bound for $\log\lambda$ in the abelian-trace case: $\lambda$ has degree at most two over $\mathbf Q^{\rm ab}$ and height $(\log\lambda)/n$.

Run the certificates with *python scripts/salem_dobrowolski.py*.

### Literature

- E. Dobrowolski, *On a question of Lehmer and the number of irreducible factors of a polynomial*, Acta Arith. 34 (1979), 391–401. [Original paper](https://doi.org/10.4064/aa-34-4-391-401)
- D. G. Cantor and E. G. Straus, *On a conjecture of D. H. Lehmer*, Acta Arith. 42 (1982/83), 97–100. The confluent determinant method is classical.
- P. Voutier, *An effective lower bound for the height of algebraic numbers*, Acta Arith. 74 (1996), 81–95. [Explicit benchmark](https://doi.org/10.4064/aa-74-1-81-95)
- F. Amoroso and U. Zannier, *A relative Dobrowolski lower bound over abelian extensions*, Ann. Scuola Norm. Sup. Pisa 29 (2000), 711–727.
