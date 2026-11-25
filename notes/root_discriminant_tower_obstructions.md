# Bounded polynomial root discriminant and the critical orbits of a tower

For a monic irreducible integer polynomial $f$ write $\mathrm{rd}(f)=|\mathrm{disc}(f)|^{1/\deg f}$. We seek fixed $C$ and unbounded degree with $\mathrm{rd}(f)\le C$. All finite examples below are totally real. The obstructions also apply without a real-root hypothesis.

## Obstruction for all integral quadratic rational substitutions

**Theorem.** Fix an infinite tower of monic irreducible totally real polynomials of unbounded degree, with binary homogenizations $F_j$. Suppose each step has the form
$$F_{j+1}(X,Y)=\pm F_j(A_j(X,Y),B_j(X,Y)),$$
where $A_j,B_j$ are integral binary quadratic forms defining a degree-two map. Then $\mathrm{rd}(f_j)\to\infty$. The maps may vary arbitrarily. Integral unimodular changes of generator on either side of each step are allowed. The theorem concerns full pullbacks. Dividing a pullback by nonunit content, taking an irreducible factor, or choosing a general new power-basis generator leads to a different problem.

### Reduction to one type of quadratic map

Let $d=\deg F$, $\rho=\mathrm{Res}(A,B)$ and $\Delta(U,V)=\mathrm{disc}_{X,Y}(VA-UB)$. The exact binary-form identity in absolute value is
$$
 |\mathrm{disc}(F(A,B))|
 =|\mathrm{disc}(F)|^2 |\rho|^{d(d-1)}|\mathrm{Res}(F,\Delta)|.
$$
Indeed the resultant between $A-r_iB$ and $A-r_jB$ is $\rho(r_i-r_j)^2$, and the discriminant of the individual quadratic is $\Delta(r_i,1)$. Multiplying gives the identity including the leading-coefficient factors, or one may homogenize the monic identity. Each nonzero integer factor is at least one in absolute value. Root discriminant is therefore nondecreasing. A step with $|\rho|>1$ increases its logarithm by at least $(d-1)\log(2)/2$.

For $|\rho|=1$ there are exactly three map types up to pre- and postcomposition by integral unimodular matrices:
$$
 [X^2+Y^2:XY],\qquad [X^2-Y^2:XY],\qquad [X^2:Y^2].
$$
To verify this, write the coefficient rows of $A,B$ in the order $(X^2,XY,Y^2)$ and take their cross product $(s,t,u)$. Then $\rho=t^2-su$. The rows span the full integral plane normal to $(s,t,u)$ since its coordinates have gcd one. A variable change acts by integral congruence on the normal's symmetric matrix with rows $(s,t),(t,u)$ (up to sign); this follows by pairing it with the coefficients of a binary quadratic. Unimodular target changes choose another basis of the plane. Up to sign a definite unimodular binary lattice is the identity lattice: in a reduced basis $|2t|\le s\le u$ gives $1\ge3s^2/4$, hence $s=u=1,t=0$. An indefinite one has a primitive isotropic vector since its determinant is $-1$. Choose a second vector pairing to one and adjust its square to 0 or 1 by adding multiples of the isotropic vector. This gives the hyperbolic plane or $\mathrm{diag}(1,-1)$. Taking the respective orthogonal coefficient planes yields the three maps displayed above.

For the minus-sign map the extra discriminant factor is $|H(2i,1)|^2$ for the transformed input form $H$. Its dehomogenization has integral nonzero leading coefficient and only real roots, so the factor is at least $4^d$. For the squaring map it is $4^d|H(0,1)H(1,0)|\ge4^d$. A right-side unimodular change preserves the binary discriminant. Each such step increases log root discriminant by at least $\log2$. Hence any bounded-root-discriminant tower would eventually consist entirely of the plus-sign type.

### Why the plus-sign type also forces growth

Write $R=[X^2+Y^2:XY]$. Absorb internal coordinate changes into their neighbors. At depth $j$ the composite is $W_j=S_1\circ\cdots\circ S_j$ followed on its input side by a unimodular change, with $S_i=M_i\circ R$. This last change preserves the discriminant and is absorbed into $M_{j+1}$ in the next step. The binary forms $F_0\circ W_j$ may have nonunit leading coefficients, which the homogeneous cost identity allows. Restart with $d_0=\deg F_0\ge3$. Define the primitive rational critical images
$$v_{j,\pm}=W_j(M_{j+1}([\pm2:1])),$$
and write $\|v\|$ for the Euclidean length of a primitive pair representing $v$. By the binary-form estimate (1) below and the exact cost $|F_0(v_{j,+})F_0(v_{j,-})|$, bounded root discriminant would imply
$$
\sum_{j\ge0} b_j<\infty,\qquad b_j=2^{-j}\max_\pm\log\|v_{j,\pm}\| .
$$

**Large matrices cannot remedy the cost.** Set $a_i=\log\|M_i\|_2$, using the Euclidean operator norm. The inverse has the same norm by unimodularity. Primitivity is preserved, and
$$\|M_iR(v)\|\ge e^{-a_i}\|v\|^2.$$
For a matrix with columns $u,v$, the larger squared length of $2u+v,-2u+v$ is at least $4\|u\|^2+\|v\|^2\ge\|M_i\|_2^2$. Thus
$$b_j\ge a_{j+1}-T_j,\qquad T_0=0,\quad
T_j=\sum_{i=1}^j 2^{i-1-j}a_i.$$
Since $2(T_{j+1}-T_j)=a_{j+1}-T_j$, summing shows $T_j\le\frac12\sum_{i<j} b_i$ and $a_{j+1}\le b_j+T_j$. In particular the matrix norms would be bounded.

**Critical images cannot stay bounded for a bounded alphabet.** Read a window of maps from innermost to outermost. Starting with the empty set of rational points, update it by $S'=MR(S\cup\{[1:1],[-1:1]\})$. Thus each point of $S'$ is some critical image propagated through the current window. Let $D(S)$ be the maximum absolute determinant of two primitive representatives. Unimodularity preserves these determinants, and
$$
 D(S')\ge 4,\qquad D(S')\ge \max_{(x,y)\in S}(|x|+|y|)^2 .
$$
The second assertion follows by pairing $R(x,y)$ with $R(1,1)$ and $R(-1,1)$: their absolute determinants are $(x-y)^2$ and $(x+y)^2$. Now $\max_{v\in S}\|v\|\ge\sqrt{D(S)}$. Once $D(S)\ge4$ a maximizer has two nonzero integral coordinates by primitivity, giving $\|v\|_1\ge\|v\|+1/2$. Hence $\sqrt{D(S')}\ge\sqrt{D(S)}+1/2$ at every subsequent step.

If all matrix norms are at most $H$, a window of $r=\lceil4H\rceil+2$ maps therefore produces some critical image of Euclidean length above $2H$. Each additional outer map sends length $V$ to at least $V^2/H$. For a window preceded by $s$ outer maps, at least one $b_j$ with $s\le j<s+r$ is at least $2^{-r}\log2$. Disjoint windows contradict summability. This completes the proof for arbitrary varying matrices. In particular no periodic word in this plus-sign type can be postcritically finite: the determinant growth alone precludes bounded critical images for its repetitions.

There is a slow quantitative lower bound for pure plus-sign towers. At a finite depth $k$, write $H=\max_{i\le k}\|M_i\|_2$ and $S_k=\sum_{j<k}b_j$. The two arguments show
$$ S_k\ge\max\{\log H,\ \lfloor k/r\rfloor\,2^{-r}\log 2\}.$$
If $H\ge(\log k)/(16\log2)$ the first term is at least $\log\log k-O(1)$; otherwise the second term is $\gg k^{3/4}/\log k$ for large $k$. With (1) this gives
$$
\log\mathrm{rd}(f_k)\ge
\frac{d_0-2-\epsilon}{2d_0}\log\log k-O_{F_0,\epsilon}(1).
$$
This estimate addresses boundedness, not Simon's stronger discriminant thresholds. The next section gives a much stronger rate for plain polynomial substitutions.


## Arbitrarily varying monic quadratic substitutions

Suppose $f_{j+1}=f_j\circ T_j$, where $T_j(x)=x^2+b_jx+c_j$ with integer coefficients. Assume the full iterates are irreducible and $\deg f_j=D_j=2^j D_0$. Then
$$
\liminf_{j\to\infty}\frac{\log\mathrm{rd}(f_j)}{\log D_j}\ \ge\ 1 .
$$
This excludes a bounded polynomial root discriminant for the entire class. The coefficients can grow, the substitutions can change at each step, and the conclusion concerns the absolute power order. Passing to a factor of a reducible iterate or changing to a different generator between substitutions is outside the statement.

**Binary-form estimate.** Let $F(A,B)=B^d f(A/B)$ for a fixed irreducible polynomial of degree $d\ge3$. For each $0<\epsilon<d-2$ there is $c>0$ such that, on primitive integer pairs with $H=\max(|A|,|B|)$,
$$
c H^{d-2-\epsilon}\le |F(A,B)|\le c^{-1} H^d. \tag{1}
$$
When $A/B$ is bounded and close to one real root, all the other linear factors stay uniformly separated from zero. Roth's rational approximation theorem bounds the exceptional factor from below by a constant times $|B|^{-2-\epsilon}$. Compactness handles the remaining bounded ratios. For large $|A/B|$ the leading term dominates. The finitely many exceptions can be absorbed in $c$ by irreducibility. This argument is ineffective in $f$ and $\epsilon$.

**Cost of one step.** The discriminant formula for quadratic composition gives
$$
|\mathrm{disc}(f_{j+1})|=|\mathrm{disc}(f_j)|^2 E_j,\qquad
E_j=4^{D_j}|f_j(c_j-b_j^2/4)|.
$$
An even $b_j$ gives $E_j\ge4^{D_j}$, since the evaluation is a nonzero integer. For odd $b_j$, the critical value is an odd numerator over 4. Each preceding monic integer quadratic squares its reduced denominator: the new numerator is the old numerator squared modulo every prime of that denominator. After applying $T_{j-1},\ldots,T_0$ in that order, the result is a primitive $A/B$ with $B=4^{2^j}$. Therefore $E_j=|F_0(A,B)|$. With $D_0\ge3$, (1) gives
$$
\log \mathrm{rd}(f_{j+1})-\log\mathrm{rd}(f_j)
\ge \left(1-\frac{2+\epsilon}{D_0}\right)\log 2+
       \frac{\log c}{2D_j}
$$
in both parity cases after decreasing $c$ to at most one. The last terms have bounded sum. Restart the argument at any fixed later iterate. Its degree can be arbitrarily large, while the ineffective constant for that restart stays fixed as $j\to\infty$. This proves the liminf statement.

## Extension to varying polynomial substitutions of bounded degree

For ordinary monic polynomial substitutions the same asymptotic exclusion holds for **every fixed degree bound**. Suppose $f_{j+1}=f_j(T_j)$ with $T_j\in\mathbf Z[X]$ monic, $2\le b_j=\deg T_j\le B$, and every $f_j$ monic irreducible totally real. Then
$$\liminf_{j\to\infty}\frac{\log\mathrm{rd}(f_j)}{\log\deg f_j}\ge1.$$

Here the input from Diophantine approximation is Schmidt's subspace theorem, specifically its standard corollary $|G(\alpha)|\gg_{\alpha,e,\epsilon}H(G)^{-e-\epsilon}$ for nonzero integral polynomials $G$ of degree at most $e$ and a real algebraic number $\alpha$ of degree greater than $e$. $H$ is the largest absolute coefficient. To recall the corollary, the numbers $1,\alpha,\ldots,\alpha^e$ are linearly independent over $\mathbf Q$. Use the linear-form version of the subspace theorem and induction on the dimension of each exceptional rational subspace; the restricted coefficients remain rationally independent in a rational basis.

For a fixed monic irreducible totally real $f_0$ of degree $D_0>e$ this gives
$$|\mathrm{Res}(f_0,G)|\gg_{f_0,e,\epsilon}
 H(G)^{D_0-e(e+1+\epsilon)}.$$
At most $e$ of the $D_0$ root evaluations of $G$ can have magnitude below $cH(G)$ for a sufficiently small fixed $c>0$. Any $e+1$ evaluations recover all coefficients by a fixed invertible Vandermonde matrix. Apply Schmidt's bound to these at most $e$ small factors, and multiply.

Write $f_j=f_0(W_j)$, $D_j=\deg f_j$, and take
$$\delta(Y)=\mathrm{disc}_X(T_j(X)-Y),\qquad
 G(X)=\mathrm{Res}_Y(\delta(Y),X-W_j(Y)).$$
Then $e=\deg\delta=b_j-1$ and its leading coefficient has absolute value $b_j^{b_j}$. The polynomial $G$ has degree $e$ and leading coefficient of absolute value $b_j^{b_j\deg W_j}$. From the roots of $\delta$, or the composition rule for resultants, the extra absolute discriminant factor is
$$
\frac{|\mathrm{disc}(f_{j+1})|}{|\mathrm{disc}(f_j)|^{b_j}}
=|\mathrm{Res}(f_j,\delta)|=|\mathrm{Res}(f_0,G)|.
$$
Restart with $D_0>(B-1)(B+\epsilon)$ to make the height exponent positive. Since $\deg W_j=D_j/D_0$, the increment of log root discriminant is at least
$$\left(1-\frac{(B-1)(B+\epsilon)}{D_0}\right)\log b_j-O_{B,f_0,\epsilon}(1/D_j).$$
Summing and restarting at arbitrarily large fixed degrees proves the assertion. The constants are ineffective, so this is an asymptotic obstruction and gives no practical first failing degree for $C=200$.

## A critical-orbit criterion for reciprocal substitutions

Set $R([A:B])=[A^2+B^2:AB]$. Let $M_i\in\mathrm{GL}_2(\mathbf Z)$, $S_i=M_i\circ R$, and write $W_j=S_1\circ\cdots\circ S_j$ with the innermost map on the right. The binary polynomial at level $j$ is $\pm F_0\circ W_j$. Assume each level yields a monic irreducible polynomial of doubled degree; in an application to our problem it must also be totally real. Both $R$ and the matrices send primitive integer pairs to primitive pairs. Total reality and monicity still have to be tested separately.

The two rational critical points of $R$ have images $[2:1],[-2:1]$. For step $j+1$ define
$$v_{j,\pm}=W_j\big(M_{j+1}([\pm2:1])\big).$$
Take primitive representatives, including $(1,0)$ at infinity, and let $h(v)=\log\max(|A|,|B|)$. By the binary discriminant law and the reciprocal-lift formula, the exact cost at this step is $E_j=|F_0(v_{j,+})F_0(v_{j,-})|$. Consequently, for $D_0\ge3$,
$$
\sup_j \mathrm{rd}(f_j)<\infty
\quad\Longleftrightarrow\quad
\sum_{j\ge0}2^{-j}\big(h(v_{j,+})+h(v_{j,-})\big)<\infty. \tag{2}
$$
Indeed the increment is $(\log E_j)/(2D_j)$, each integer cost is positive, and (1) compares $\log|F_0(v)|$ from both sides to $h(v)$ up to constants summable with these weights. Monicity requires the additional exact condition $|F_0(W_j([1:0]))|=1$.

**Periodic sequences.** For a period of length $r$, let $\Phi=S_1\circ\cdots\circ S_r$ of degree $b=2^r$. The rational points
$$v_{i,\pm}=S_1\circ\cdots\circ S_{i-1}\big(M_i([\pm2:1])\big),\quad 1\le i\le r,$$
are exactly its set of critical values, allowing repetitions (here $i$ numbers positions in the period). This follows from the chain rule and surjectivity of the suffix maps on the projective line. The sequences in (2) are precisely the forward $\Phi$-orbits of these points, taken one residue class at a time.

For a fixed matrix the height after one step differs from twice the old height by a bounded amount, since $A^2+B^2\ge H^2$ and both the matrix and its inverse have fixed integer entries. It follows that $|h(\Phi(v))-b h(v)|$ is uniformly bounded. Telescoping gives the canonical height $\hat h(v)=\lim b^{-k}h(\Phi^k(v))$ and $|h(v)-\hat h(v)|=O_\Phi(1)$. On rational points $\hat h(v)=0$ exactly for finite forward orbits, since bounded rational height gives a finite set. If one orbit is infinite, its weighted terms in (2) approach a strictly positive limit along its residue class. If all orbits are finite, all the costs are bounded. Thus, subject to the tower hypotheses, bounded polynomial root discriminant is **equivalent to the period map being postcritically finite**, meaning every critical point has finite forward orbit. This also applies after discarding a finite prefix and restarting with the resulting seed.

**Exact search for period maps.** For a word of matrices choose $K=\max_i\|M_i^{-1}\|_\infty$ (maximum absolute row sum). If a primitive pair has height $H>K$, every next ordinary height is at least $H^2/K>H$. Thus a critical orbit either escapes or repeats below the cutoff. The test terminates and uses only integers. Entries are bounded in absolute value and matrices are counted modulo overall sign.

| Period length | Entry bound | Ordered words | Postcritically finite |
|---:|---:|---:|---:|
| 1 | 16 | 2548 | 0 |
| 2 | 3 | 13456 | 0 |
| 3 | 1 | 8000 | 0 |
| 4 | 1 | 160000 | 0 |

All 184004 word presentations fail the necessary dynamical condition, independently of the monicity and reality checks. The algebraic determinant argument above covers all periods and all entry sizes. The table provides independent finite checks. The script also checks the escape-window mechanism exhaustively, including nonperiodic choices. At entry bound 16, the numbers of surviving critical-image sets after one through four steps are 938, 3634, 946 and 0. Run *python scripts/periodic_discriminant_search.py*.

## A fixed-$C$ census in a non-Galois family

For $f_n=\Psi_n$, the minimal polynomial of $2\cos(2\pi/n)$ of degree $m$, consider $P_n=X^m f_n(X-X^{-1})$. Its degree is $2m$, every root is real, and its discriminant is
$$\mathrm{disc}(P_n)=\mathrm{disc}(f_n)^2\,\mathrm{Res}(f_n,X^2+4).$$
The resultant is at least $4^m$. With the exact integer comparison $\mathrm{disc}(P_n)\le C^{2m}$, the largest degrees of irreducible full lifts in the family are as follows.

| Fixed $C$ | Largest degree | Canonical conductor | Polynomial root discriminant |
|---:|---:|---:|---:|
| 60 | 32 | 51 | 59.83604005 |
| 60 | 32 | 120 | 56.02870785 |
| 200 | 128 | 255 | 198.87949898 |

**Completeness.** Use canonical conductors (omit $n\equiv2\bmod4$, which gives the same field and reflected trace polynomial). For conductors with at least two distinct prime divisors the cyclotomic formula gives $\mathrm{rd}(f_n)=n/\prod_{p\mid n}p^{1/(p-1)}\ge\sqrt n$. Here an even canonical conductor is divisible by four. At a prime power the additional endpoint correction divides this expression by at most 2: it is $p^{1/\varphi(n)}$ at an odd prime and $2^{1/m}$ at 2. Thus $\mathrm{rd}(P_n)\ge2\,\mathrm{rd}(f_n)$ bounds the conductor by $C^2/4$ in the first case and $C^2$ in the second. The linear seeds at conductors 1 and 2 cannot affect the maximum. The search includes reducible full lifts in its initial upper-bound calculation, then certifies irreducibility and nonnormality at every maximizing conductor. It makes no claim about individual factors of discarded lifts.

**Field check.** At conductor 255 take the split prime 1021 and the primitive 255th root 625 modulo 1021. The exponents 1 and 2 give the trace residues 411 and 454, for which $\beta^2+4$ has opposite nonzero quadratic characters. Thus the quadratic extension splits at one prime of the trace field and is inert at the other. Its degree is two and it is non-Galois over $\mathbf Q$. Since $t-1/t=\beta$, any root generates the extension. The same checks hold at $(n,q,w)=(51,103,4)$ and $(120,241,49)$, with exponent pairs $(1,4)$ and $(1,23)$ respectively. The certificate script also constructs each extremal polynomial and computes its discriminant directly. Run *python scripts/fixed_c_reciprocal.py*.

The Diophantine inputs are K. F. Roth, *Rational approximations to algebraic numbers*, Mathematika 2 (1955), 1–20, and W. M. Schmidt, *Diophantine Approximation*, Lecture Notes in Mathematics 785 (1980). The script *python scripts/quadratic_pullback_check.py* independently checks the quadratic discriminant formula on 480 maps, including maps of nonunit resultant.
