# Two bilinear rank defects force a quadratic fan

Let `A=K[x,y,z]/I` have positive length `n`, over any field. Assume

\[
 \dim\langle1,x,y,xy\rangle_A\le3,\qquad
 \dim\langle1,x,z,xz\rangle_A\le3.                       \tag{B1}
\]

Then

\[
                          g(I)=O(n^2).                  \tag{B2}
\]

The quadratic order is attained over every field by the existing
three-pole family. No coordinate is required to generate the algebra,
and no reducedness, separability or generic-evaluation assumption is
made. This gives an inverse statement for two small monomial rank
defects, rather than assuming a curve carrier in advance. The proof
uses elementary bilinear elimination and the existing planar theorem;
its exact priority is unassessed. It does not settle unrestricted `G_3`.

## The algebraic dichotomy

A coordinatewise affine substitution preserves the initial ideal for
each global monomial order: expanding a monomial replaces it by itself
times a nonzero scalar and by proper divisors. Leading monomials of
every polynomial are therefore preserved, as are those of the entire
ideal. This is the established linear-shift principle; see
[Dimitrova--He--Robbiano--Stigler, Proposition 2.7](https://arxiv.org/html/1704.03680v2).

First, any nonzero relation supported on `1,x,y,xy` puts `xy` in
every global monomial initial ideal. If its `xy` coefficient is
nonzero, that term strictly dominates all others. Otherwise its leading
monomial is one of `1,x,y`, which divides `xy`. The same applies to
the other coordinate pairs. Consequently, if all three pairs satisfy
the rank defect, every standard set is an axis allocation and

\[
                       g(I)\le\binom{n+1}{2}.           \tag{B3}
\]

Choose relations `F(x,y),G(x,z)` from (B1). If a coordinate is a
scalar in `A`, eliminating it preserves all global initials and the
all-field planar bound gives `O(n log(n+1))`. Otherwise a relation
with no mixed term has nonzero coefficients at both its variables.
Using it to substitute for `x` in the other relation gives a nonzero
relation on `1,y,z,yz`. Thus (B3) applies in this case.

It remains that both chosen mixed coefficients are nonzero. Normalize
the relations, using scalars `a,b,c,d,e,f` in `K`, as

\[
 F=(x-a)(y-b)-c,\qquad G=(x-d)(z-e)-f.                   \tag{B4}
\]

Their elementary elimination identity is

\[
 (y-b)G-(z-e)F
  =(a-d)(y-b)(z-e)+c(z-e)-f(y-b).                       \tag{B5}
\]

Unless `a=d,c=f=0`, its right side is a nonzero bilinear relation
between `y,z`, again proving (B3). In the remaining case the same
coordinatewise affine shift gives

\[
                             xy,xz\in I.               \tag{B6}
\]

These exhaust the cases over every field. In particular, it is not
valid to omit the common-factor case in (B4) or infer the third rank
defect directly from just (B1).

## The plane-plus-axis exception has planar growth

**Theorem.** Every finite quotient satisfying (B6) obeys

\[
                         g(I)=O(n\log(n+1)).             \tag{B7}
\]

The order is sharp, because taking `x` in `I` includes all planar
quotients. The estimate permits nonreduced quotients and relations
identifying an axis socle with a planar socle.

Put `U=(x)A=span{x,x^2,...}`, `B=K[y,z]` inside `A`, and `r=dim U`.
Then `A=U+B`, and

\[
                            \dim(U\cap B)\le1.          \tag{B8}
\]

Here is a linear proof. If `x=0`, it is immediate. Otherwise
multiplication `T` by `x` on `U` is cyclic with generator `x`.
For `v` in `U cap B`, write `v=c*1+w` with `w` a combination of
positive `y,z` monomials. Equation (B6) gives `Tv=c*x`, so
`U cap B` is contained in `T^{-1}(K*x)`. If `T` is invertible this
preimage has dimension one. If not, its cyclic generator cannot lie
in `im T`: otherwise its entire cyclic span would lie there and `T`
would be surjective. The preimage is then `ker T`, also of dimension
at most one for a cyclic operator. This proves (B8).

Every initial standard set consists of a prefix `x,...,x^a` and a
planar downset containing `1`. Independence gives `a<=r`; the planar
part has at most `dim B=n-r+dim(U cap B)` elements. Therefore

\[
                              a\in\{r-1,r\}.             \tag{B9}
\]

Only feasible nonnegative values are retained. This bound is on one
fixed quotient, including all of its weights.

For either feasible `a`, the literal core
`C_a={1,x,...,x^a}` is independent. Contract its columns once in the
actual vector matroid. Write `k=n-a-1`. The remaining standard
monomials lie in

\[
 E_k=\{y^b z^c:(b,c)\ne(0,0),\ (b+1)(c+1)\le k+1\}.     \tag{B10}
\]

Indeed all their positive divisors also occur among the `k` remaining
terms. This finite contracted matroid has rank exactly `k` whenever
the class is nonempty. Every actual standard set in this class gives
a strictly optimal basis under the two costs `w_y,w_z`: another
contracted basis with smaller or equal cost would contradict the
strict global optimum after adding `C_a`, whose cost is fixed within
the class.

The [general planar hyperbolic-matroid theorem](planar_dyadic_bound.md)
with rank `k` and hyperbolic parameter `H=k+1` bounds these minima by
`O(1+k log(k+2))`. It applies to arbitrary vector matroids, so the
contracted core is not being treated as an algebra ideal. Summing
over the at most two values in (B9) proves (B7). No global support
identity outside these fixed-axis classes is needed.

Together, (B3), (B7) and planar elimination prove (B2).

## What a cyclic coordinate adds

When `1,x,...,x^(n-1)` is a basis and `n>=3`, (B1) has a useful
operator interpretation. Each other coordinate `h` is exactly one of:

* a scalar;
* a nonconstant Möbius function of `x`, with a unit denominator;
* a scalar plus a nonzero rank-one multiplication operator supported
  on the kernel of `x-a`, for some `a` in `K`.

To prove this, a rank relation is
`(alpha+beta*x)+(gamma+delta*x)*h=0`. If its two-by-two coefficient
matrix is nonsingular, the numerator and denominator are coprime
linear polynomials. Their Bézout identity and this relation make the
denominator a unit, giving the Möbius case. If the matrix is singular,
the nonconstant case factors as `(x-a)(h-b)=0`. For
`A=K[t]/F(t)`, its annihilator is one-dimensional, generated by
`F(t)/(t-a)`; multiplication by this element has rank one. Rank at
most two forces `h` scalar, since `1,x,x^2` are independent.

In the common-factor exception both nonconstant residual coordinates
are proportional in that one-dimensional annihilator. They have a
linear relation, so the third pair also satisfies the rank defect.
Thus the sharper exact bound (B3) holds under the cyclic hypothesis.
This interpretation is an elementary consequence of a four-column
dependence, not an asserted new general inverse-product theorem.

## Sharpness, algorithmic consequence and limits

For every field and `n>=2`, take a monic `F(t)` of degree `n` avoiding
`0,1`, for example `t^n-t+1`, and

\[
 I_F=(xy-1,(x-1)z-1,F(x)).                              \tag{B11}
\]

The [existing three-pole proof](weighted_facet_incidence.md#quartic-incidence-is-sharp-for-actual-finite-quotients)
shows that every axis allocation is a basis and is coherent. Hence
`g(I_F)=binom(n+1,2)`, independently of factorization and repeated
roots. Both tests in (B1) hold. This establishes the sharp quadratic
order of the whole class, not a polynomial gain over generic cuts.

Substituting (B2) in the existing
[output-sensitive enumeration theorem](oracle_fan_enumeration.md)
gives

\[
 O(n^3\log(n+1)^2+n^{\omega+2}\log(n+1))                 \tag{B12}
\]

field and exact geometric operations for this class. Coefficient bit
costs and practical superiority remain unestablished. No core CLI or
auditor change is made.

Small-product inverse theorems are established research. In particular,
[Mirandola--Zémor, Theorems 2.5 and 3.3](https://arxiv.org/html/1501.06419)
treat full-support codes and coordinatewise products; their hypotheses
must not be substituted for arbitrary nonreduced quotient ranks.
[Beck--Lecouvey, Theorem 4.1](https://arxiv.org/html/1504.02287)
has explicit unit and finite-subalgebra hypotheses. Those statements
were inspected as scope comparisons, not as a completed priority
assessment or inputs needed for the elementary proof above. Bilinear
elimination, cyclic annihilators and the existing planar bound supply
the argument. The exact theorem-level prior-art comparison remains
open.

The unrestricted target still requires coupled dependencies beyond
these two small rectangles, or a genuinely larger realizable family.
In particular, a superquadratic family must have at least two of its
three coordinate-pair boxes `span{1,x_i,x_j,x_i*x_j}` of full rank
four: any two deficient pairs share an anchor and satisfy (B1).
This is a necessary condition, not a large-fan certificate; fully
generic ideals already satisfy all three full-rank tests.
This result identifies a complete nonamplifying rank-defect regime;
it does not lower the original goal's completion criterion.

Run from the repository root:

    python3 -B scripts/run_bilinear_rank_defect_checks.py
    ./build.sh
    python3 tests/test_smoke.py

The independent controls check the symbolic elimination identity,
finite quotient operator models, complete fans, all-field elementary
cases through finite specializations, and the actual fixed-core
planar contractions. Computational agreement supports these uniform
proofs; it does not establish novelty or the unrestricted endpoint.
