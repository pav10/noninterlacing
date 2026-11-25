# Testing a fixed bound for polynomial root discriminants

For a fixed $C$, the target in degree $d$ is $0<|\mathrm{disc}(f)|\le C^d$ for monic irreducible totally real $f$, with $d$ unbounded. A tower of fields with bounded root discriminant supplies suitable fields only if the absolute power-basis indices can be kept exponentially bounded. The following experiment computes the polynomial discriminant at each step.

## An exact cost for doubling the degree

Start with a monic irreducible $h$ of degree $d$, with all roots real and outside $[s-2,s+2]$ for an integer $s$. Set
$$g(X)=X^d h(s+X+X^{-1}).$$
Each quadratic factor has two distinct real roots. If $g$ is irreducible, it gives a totally real extension of doubled degree. The within-pair discriminants and the between-pair resultants give
$$
|\mathrm{disc}(g)|=\mathrm{disc}(h)^2\,|h(s-2)h(s+2)|.
$$
One may first replace $h$ by a unimodular fractional-linear transform whose leading coefficient is a unit and normalize its sign. The binary discriminant transformation law preserves its discriminant exactly.

Write $E_j=|h_j(s_j-2)h_j(s_j+2)|$ for each nonzero integer cost after the change of generator. Thus along an irreducible tower in this construction,
$$
\log\mathrm{rd}(g_k)=\log\mathrm{rd}(g_0)+
  \sum_{j=0}^{k-1}\frac{\log E_j}{2^{j+1}\deg(g_0)} .
$$
**Bounded root discriminant is equivalent to convergence of this series.** Cost one is sufficient at any step. A uniform bound on the integer costs would also suffice, as would $\log E_j=O(2^j/(1+j)^2)$. A search can therefore track a fixed budget $\log C-\log\mathrm{rd}(g_0)$ and prune branches as soon as they spend it.

## A degree doubling with exactly preserved root discriminant

For every integer $a\ge4$, take $f_a=(X^2-4)(X-a)+1$. Its roots lie respectively below $-2$, between 2 and 3, and between 3 and $a$. Its endpoint values are both 1. It is irreducible since an integer root would satisfy $X^2-4=\pm1$. The sextic is
$$
g_a=(X^2-1)^2(X^2-aX+1)+X^3.
$$
If the quadratic lift split over the cubic field, its two reciprocal roots would have monic cubic minimal polynomials reciprocal to each other, since the trace recovers the full cubic field. They would be distinct because the trace generates that field. Write one as $h$. The identity $g_a(1)=1$ forces $h(0)=1$ and $h(1)=\pm1$, using $g_a(X)=h(X)X^3 h(1/X)/h(0)$. Writing $h=X^3+uX^2+vX+1$, the coefficient of $X^5$ gives $u+v=-a$, whereas $h(1)=2+u+v=\pm1$ forces $a=1$ or 3. This proves irreducibility for the entire stated family.

For $a=4$ the cubic discriminant is 1957 and the sextic discriminant is $1957^2$. Both power orders have root discriminant $12.5082612\ldots$. In this example the cubic discriminant is squarefree, and the relative quadratic basis has unit discriminant, so both power orders are maximal as well. Varying $a$ gives infinitely many one-step examples with varying constants and fixed degrees.

## Can the step be repeated?

For a sextic $g$ write $G$ for its binary homogenization. PARI's unconditional Thue solver gives all primitive integer vectors $U$ with $G(u,v)=\pm1$. A second cost-one step via a monic unimodular change requires $v_0,w,v_0+4w$ all in $U$ with $|\det(w,v_0)|=1$. Here $w$ is the leading-coefficient vector and the other two vectors give the two endpoints.

For $a=4,\ldots,12,16,24$, the solver returned eight unit vectors for each sextic and zero such triples. For $a=4$ the vectors are $\pm(1,0),\pm(0,1),\pm(1,1),\pm(1,-1)$. Allowing costs above one, the **least** second cost under these transformations and irreducibility is **1801**. For each unit vector the possible second matrix columns differ by its integer multiples, giving integer translations, and changing orientation reflects the variable. The search tests all translates up to the point where both endpoints lie on the same exterior side of the roots and the cost exceeds 1801. Further translation on that side strictly increases the cost.

| Degree | Polynomial root discriminant | Below fixed $C=200$ |
|---:|---:|:---:|
| 3 | 12.5083 | yes |
| 6 | 12.5083 | yes |
| 12 | 23.3609 | yes |
| 24 | 74.5498 | yes |
| 48 | 363.7465 | no |

Only the second cost in this table is proved minimal. The last two moves were selected by a finite search. The degree-48 branch has spent the budget for $C=200$. An infinite branch with summable costs remains unresolved.

Run *python scripts/unit_discriminant_steps.py* for the exact discriminants, Thue checks, minimum-cost certificate and explicit changes of variables.
