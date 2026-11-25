# High-degree totally real polynomials with small discriminant

Write $D(f)=|\operatorname{disc}(f)|$ and $\operatorname{rd}(f)=D(f)^{1/d}$, where $d=\deg f$. We seek monic irreducible polynomials over the integers with distinct real roots. The threshold $D(f)<d^d$ is equivalent to $\operatorname{rd}(f)<d$. The stronger goal $D(f)\le C^d$ asks for bounded root discriminant in unbounded degree.

## 1. An explicit family below the non-interlacing threshold

Let $\Psi_n$ be the minimal polynomial of $2\cos(2\pi/n)$. For odd $n$ with at least two prime divisors, $d=\varphi(n)/2$. The reciprocal lift $z^d\Psi_n(z+z^{-1})$ is the cyclotomic polynomial $\Phi_n(z)$. For any monic trace polynomial $f$, the roots of the lift occur in pairs given by $z^2-rz+1$ for each root $r$ of $f$. The within-pair discriminant is $r^2-4$ and the resultant between two quadratics is $(r_i-r_j)^2$. Consequently

$$
|\operatorname{disc}(z^d f(z+z^{-1}))|
= |f(2)f(-2)|D(f)^2.
$$

For these odd conductors the endpoint values have absolute value one. The standard cyclotomic discriminant formula then gives

$$
\operatorname{rd}(\Psi_n)=\frac{n}{\prod_{p\mid n}p^{1/(p-1)}},
\qquad
\frac{\operatorname{rd}(\Psi_n)}{d}
=2\prod_{p\mid n}\frac{p}{p-1}p^{-1/(p-1)}.
$$

| Conductor | Degree | Root discriminant | Root discriminant / degree |
|---:|---:|---:|---:|
| 105 | 24 | 29.3114 | 1.22131 |
| 1155 | 240 | 253.6822 | 1.05701 |
| 15015 | 2880 | 2663.2053 | 0.92472 |
| 45045 | 8640 | 7989.6160 | 0.92472 |
| 135135 | 25920 | 23968.8480 | 0.92472 |
| 255255 | 46080 | 37927.1770 | 0.82307 |

The exact discriminant at $n=15015$ is $3^{1440}5^{2160}7^{2400}11^{2592}13^{2640}<2880^{2880}$. Multiplying 15015 by arbitrary powers of 3 preserves the ratio and gives examples in arbitrarily high degree. Taking $n$ to be the product of the odd primes up to $y$ makes the ratio tend to zero. The logarithm of each factor in the product is at most $(1-\log p)/(p-1)$, and the sum tends to minus infinity by divergence of the sum of reciprocal primes.

These strict inequalities prove non-interlacing without high-dimensional enumeration. A monic integer interlacer $g$ would give positive residues $g(r)/f^\prime(r)$ summing to 1. Their product is $|\operatorname{Res}(f,g)|/D(f)\ge1/D(f)$, whereas the arithmetic–geometric mean inequality bounds it above by $d^{-d}$. The cyclotomic examples have all roots in $(-2,2)$. Their reciprocal lifts have Mahler measure 1 and furnish a boundary case for the Salem search.

Reproduce the table with:

```sh
python scripts/cyclotomic_discriminants.py 105 1155 15015 45045 135135 255255
```

## 2. Bounded root discriminant and the power-basis index

An infinite tower of totally real fields unramified at finite primes has constant **field** root discriminant. For a generator $\alpha$ with minimal polynomial $f$, the polynomial root discriminant is the field root discriminant times $[\mathcal O_K:\mathbf Z[\alpha]]^{2/d}$. A first computational problem is to choose a tower with explicit equations and try to keep $\log[\mathcal O_K:\mathbf Z[\alpha]]/d$ bounded, using integral bases and exact index calculations.

Screen towers for completely split primes before investing in generator searches. If a fixed rational prime $p$ splits completely, the $d$ conjugates of every integral generator lie in $\mathbf Z_p$. The valuation of the discriminant is twice the sum of the pairwise valuations of their differences. Modulo $p^k$, the pigeonhole principle gives the following lower bound, already sufficient with just $k=1$. Write $d=qp+r$, with $0\le r<p$. Then

$$
v_p(D(f))\ge 2\left[(p-r)\binom q2+r\binom{q+1}2\right].
$$

The right side is $d^2/p-O(d)$. This forces polynomial root discriminants to grow in such a tower regardless of the choice of generator. The same screen applies if $p$ is unramified and a positive fraction of the $d$ embeddings lie above primes of residue degree at most a fixed $B$. These conjugates all lie in the unramified extension of $\mathbf Q_p$ of degree $\operatorname{lcm}(1,\ldots,B)$. Its residue field has fixed finite size, so the same pair count for this subset grows quadratically in $d$.

## 3. Move from cyclotomic roots to one Salem outlier

The target for a Salem trace polynomial is one root $r>2$ and all other roots in $(-2,2)$. Its reciprocal lift has Mahler measure $(r+\sqrt{r^2-4})/2$. Thus for this search the small measure objective belongs to the lift, while the discriminant and interlacing tests can be applied to the totally real trace polynomial. An irreducible trace polynomial of degree at least two in this region gives an irreducible lift. If its defining quadratic split over the totally real trace field, the two roots would remain real under every embedding, contrary to any embedding taking the trace to $(-2,2)$.

A first search used $\Psi_n(x)+a+bx$ for canonical conductors through 300, with $|a|\le4$ and $|b|\le2$. Among the 1545 polynomials of degree at least two with the required endpoint signs, certified root isolation and factorization gave 23 distinct Salem trace factors. The largest trace degree was 7 and none met $D(f)<d^d$. The lowest measure was Lehmer's number, from $\Psi_{36}-x$. The next lowest was approximately 1.20261674, from $\Psi_{40}-x-1$, which has trace polynomial $x^7-8x^5+19x^3-12x-1$. Exact interlacing tests gave monic witnesses for both. In degree 7 the witness was $x^6-6x^4+9x^2-1$.

A larger search can rank seeds by endpoint signs and discriminant margin. Alternating signs at certified critical points of the original polynomial and at $\pm2$ can certify preservation of the interior roots for some perturbations, with direct root isolation available for ambiguous cases.

A second baseline is the family $z^n p(z)\pm p^*(z)$ for $p(z)=z^3-z-1$, removing the linear factors at $\pm1$ and factoring in the trace variable. A scan of exponents 1–35, 47, 63, 95, 127, 191 and 255 recovered Lehmer's number at exponent 8 and reached trace degree 128. The high-degree values cluster around the plastic constant $1.324718$, and the scanned Salem trace factors were above the $d^d$ threshold. This makes the family a control for a search across different sparse seeds or perturbations.
