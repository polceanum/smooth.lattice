# Uniform free-prefix filtering and cheap monomial rank obstructions

The unresolved collective wall-cost estimate must use ranks in one
quotient. This note gives a necessary rank defect on an explicit finite
set for each possible wall, a restricted whole-fan consequence, and a
precise obstruction to importing the ordinary three-dimensional
`k`-set bound into arbitrary matroids. None settles unrestricted `G_3`,
produces a polynomial separating family of ideals, or establishes novelty.

## The entire weight-independent cheaper set

Let `I` have colength `n` over a field, and let `q` be a primitive mixed
integer direction. Put `u=q^+`, `v=q^-`, and

\[
 P_q=\operatorname{conv}(\operatorname{Div}u\cup
                         \operatorname{Div}v),\qquad
 L_q=(P_q\cap\mathbb Z^3)\setminus\{u,v\}.                 \tag{R1}
\]

Divisors here mean their exponent vectors. For every positive `w` in
`q^perp`, every monomial in `L_q` has cost strictly below
`t=w.u=w.v`. The maximizing face of `P_q` is exactly the segment
`[u,v]`: positivity makes each divisor box have its unique maximum at
its tip. Primitivity makes `u,v` the segment's only lattice points.

In fact `L_q` is the exact common strictly cheaper set for these weights.
On the support of `q` the polytope has inequalities

\[
 a_i\ge0,\qquad a_i/u_i+a_j/v_j\le1
       \quad(u_i>0,\ v_j>0),                              \tag{R2}
\]

and coordinates off that support are zero. For a nonnegative exponent
outside this polytope, a violated pair inequality gives a nonnegative
balanced weight concentrated on that pair. A small positive balanced
perturbation preserves the strict violation. An off-support coordinate
can instead be given arbitrarily large price. Thus some admissible
weight makes that monomial more expensive than `t`.

**Rank obstruction.** If `q` occurs as an actual positive wall of `I`,
then

\[
 \operatorname{rank}_I L_q\le n-1.                        \tag{R3}
\]

Choose a weight in the relative interior of that wall, avoiding every
extra monomial tie not parallel to `q`. An exchanging pair differs by
`kq`, `k!=0`, so its common cost is at least `t`: one exponent contains
`u` and the other contains `v`. If the strictly cheaper columns `L_q`
already spanned the quotient, the greedy basis would finish before any
such exchange. Nearby weights would have the same basis. This
contradicts the wall. The finite universal ground suffices for the
generic choice of weight; columns outside that ground may still be
used to compute the rank in (R3).

The [compulsory core](wall_coexistence.md) `C_q` is contained in `L_q`
and has independent columns of cardinality `nu(q)-1`. Hence the image
of `L_q` in that *same* contracted vector space satisfies

\[
 \operatorname{rank}_{I/C_q}(L_q\setminus C_q)
       \le n-\nu(q).                                     \tag{R4}
\]

This is a rank condition, not literal containment of `L_q` in either
standard set. Its span is not asserted to be an ideal or a faithful
global fan contraction. In particular, the earlier local-support scope
guards still apply.

For `q=(a,-b,-c)`, with `a,b,c>0`, `P_q` is the pyramid over the
`b` by `c` coordinate rectangle with apex `(a,0,0)`, and

\[
 |L_q|=\sum_{i=0}^{a}
   (\lfloor b(a-i)/a\rfloor+1)
   (\lfloor c(a-i)/a\rfloor+1)-2.                         \tag{R5}
\]

The box with caps `floor(|q_i|/2)` lies inside `L_q`. It therefore
has rank below `n` as well. Its cardinality is at least
`A(q)B(q)/8`, in the notation of the
[free-wall threshold](free_wall_threshold.md). A direction whose free
threshold is much larger than `n` thus needs many actual linear
dependencies among these universally cheaper monomials. Core
independence alone does not supply this information.

The minimal socle realization of a wall shows that (R3) can be sharp
even when `|L_q|` has quadratic order in `n`: every monomial outside
the two divisor boxes is zero there. Conversely, a uniform evaluation
matroid can have the same independent core while `L_q` already has
rank `n`, excluding that wall. These are different quotients, not a
coexistence theorem.

## Exact scale of uniform free-prefix filtering

For a finite quotient in `d` variables, define

\[
 s(I)=\sup\{|B|:B\text{ is an anchored monomial box of rank }<n\}.
                                                               \tag{S1}
\]

Take `s=0` if no such box exists, and permit `s=infinity`. Let `M_*(I)`
be the supremum, over positive generic weights, of the number of
ordinary lattice terms scanned through the last greedy basis selection.
Then

\[
       1+s(I)/d^d\le M_*(I)\le1+(d+1)^d s(I).            \tag{S2}
\]

The constants depend only on dimension. For the upper bound, if the
last selection has cost `tau>0`, its strictly cheaper box with caps
`floor(tau/((d+1)w_i))` cannot span. The enclosing-box inequality (R6b)
below bounds the ordinary rank of that selection by `(d+1)^d s`.
For the lower bound take any nonsaturating box with caps `a` and
prices `w_i=1/(a_i+1)`. Every monomial of cost below one lies inside
that box and cannot yet span. The subbox with caps `floor(a_i/d)`
has all costs strictly below one and at least `|B|/d^d` monomials.
A small generic perturbation preserves its strict separation from all
monomials outside `B`, whose minimum is controlled by the finite axis
border. Greedy completion must follow these terms. Take the supremum.
The separate `n=1` case has `s=0`, `M_*=1`.

These parameters are finite **if and only if each coordinate separately
is cyclic**: its pure powers span the entire quotient. If a coordinate
is not cyclic, its arbitrarily long axis boxes do not span; making its
price much smaller than the other prices gives arbitrarily long
uncompleted free prefixes. Conversely, if every axis spans, its first
`n` powers are a basis by the univariate dependence recurrence. Every
nonsaturating box then has all caps at most `n-2`. In this case
`M_*<=n^d` also follows directly: greedy finishes by cost
`(n-1)*min_i w_i`, and every term through that cost has every exponent
at most `n-1`.

The worst exponent `d` is attained. In the all-field algebra
`K[t]/(t^n)` with every coordinate equal to `t`, take generic weights
very close to each other. Before selecting the cheapest variable's
`(n-1)`-st power, every monomial of total degree at most `n-2`
has been scanned. This gives at least `binom(n+d-2,d)+1` terms.
Also `s` is the maximum of `prod(a_i+1)` subject to `sum a_i<=n-2`,
which has order `n^d`. This algebra has only `d` initials when `n>1`.
Thus high scanning cost alone does not give a large-fan lower bound.

The existing ordered-prefix count and (S2) give, whenever `s` is finite,

\[
 g(I)=O_d((s(I)+1)^{d-1}\log(s(I)+2)^{d-1}).             \tag{S3}
\]

In particular the **whole** class `s(I)=O_d(n)` has the sharp generic
order `O_d(n^{d-1} log(n+1)^{d-1})`. This is a quantitative reduction
to the original free exponent orders, rather than a statement that
its columns are uniform. It leaves quotients with superlinear `s`,
and quotients with a noncyclic coordinate, for a separate joint-rank
argument; the diagonal example prohibits interpreting these as
lower-count witnesses.

## Rectangular rank saturation and a nonuniform family

Say that a quotient has rectangular rank saturation if every anchored
monomial box `B` has rank `min(n,|B|)`. This is weaker than requiring
every set of at most `n` monomials to be independent. It is a finite
condition: boxes with all caps at most `n-1` suffice. Each coordinate
axis box of size `n` is then a basis, so every larger-cap box already
spans.

There is a stronger whole-fan consequence than merely bounding the
products of its wall directions. For fixed dimension `d>=2`, every
positive generic weight selects its complete quotient basis within the
first

\[
                  M=(d+1)^d(n-1)                         \tag{R6a}
\]

ordinary lattice terms, when `n>1`. The case `n=1` needs one term.
To prove this, let `tau>0` be the cost of the last selected monomial.
The box with caps `floor(tau/((d+1)w_i))` has all costs strictly below
`tau`. It cannot span the quotient, since otherwise greedy selection
would have finished earlier. Rectangular saturation makes its
cardinality at most `n-1`. On the other hand the number of lattice
points of cost at most `tau` is at most

\[
 \prod_i(\lfloor\tau/w_i\rfloor+1)
 \le(d+1)^d\prod_i(\lfloor\tau/((d+1)w_i)\rfloor+1)
 \le M.                                                   \tag{R6b}
\]

Thus greedy filtering of one fixed free-prefix length works at *every*
positive generic weight. The complete ordered free prefix determines
the quotient basis. The
[ordered-prefix theorem](../paper/corner_cut_growth/main.tex), with its
explicit aspect-chart proof, now gives

\[
 g(I)=O_d(n^{d-1}\log(n+1)^{d-1}),\qquad
 g(I)=O(n^2\log(n+1)^2)\quad(d=3).                       \tag{R6}
\]

For clarity, that prefix upper bound uses `O(log(M+1)^(d-1))` dyadic
aspect charts. Each chart contains the full free prefix in a box of
volume `O_d(M)`, has `O_d(M)` distinct integer difference hyperplanes,
and hence has `O_d(M^(d-1))` central-arrangement chambers. It counts
ordered prefixes, rather than only their final corner cuts. For a
non-weight global monomial order, a positive generic weight realizes
its finite quotient comparison system, so the same initial ideals are
covered. Generic point ideals belong to this class and attain the
already proved corner-cut scale. Hence (R6) is sharp in order for the
class, with that existing lower-count theorem as input. Its exact
priority remains to be established.

The [existing output-sensitive enumeration theorem](oracle_fan_enumeration.md)
then yields, in three variables,

\[
 O(n^3\log(n+1)^2+n^{\omega+2}\log(n+1)^3)              \tag{R6c}
\]

field and exact geometric operations, from one reduced basis or
commuting multiplication matrices. The algorithm need not be told that
the saturation promise holds; the smaller output bound controls its
running time on this class. This excludes coefficient bit costs and
promises no practical speedup. The oversampling statement itself is a
uniform reduction to the original free exponent orders, including
positive prime-log weights; it is not a prime-distribution theorem.
Neither this class bound nor its enumeration consequence resolves
unrestricted `G_3`.

The class contains nonuniform families. Let `F` be any absolutely
irreducible polynomial over an infinite field with positive coordinate
degrees `d_1,d_2,d_3`, and take `n` sufficiently generic distinct points
on its hypersurface, where

\[
              n\le\prod_i(d_i+1)-1.                     \tag{R7}
\]

For a box with caps `a_i`, the dimension of its polynomial space modulo
`F` is exactly

\[
 \prod_i(a_i+1)-\prod_i\max(a_i-d_i+1,0).                \tag{R8}
\]

Indeed every relation is a multiple of `F`, and coordinate degree is
additive for nonzero polynomial products. If some cap is smaller than
`d_i`, the second term is zero. Otherwise (R8) is at least
`prod_i(d_i+1)-1`: expansion in `a_i-d_i` has nonnegative
coefficients. Thus its value is at least `min(n,|B|)`.
Evaluation at sufficiently generic points has rank the smaller of `n`
and that function-space dimension. Only finitely many boxes need be
considered, so all these open conditions hold simultaneously. This
proves rectangular rank saturation and (R6). Points may be taken over
an algebraically closed infinite field or a suitable extension; no
rational specialization over every ground field is promised.

For example, in characteristic zero use the smooth, absolutely
irreducible affine Fermat hypersurface
`F=x^D+y^D+z^D-1`, `D>=2`, with `n=(D+1)^3-1` generic points.
Its homogenization is smooth by its partial derivatives; a reducible
projective surface would have singular points where its components meet.
Its four displayed monomials have a nontrivial relation and all belong
to the length-`n` universal ground. Thus the residue matroid is
nonuniform for an infinite family in this class. This example proves
that the hypothesis is weaker; it does not give a lower count.

## Commuting multiplication constrains the saturation diagram

In any cyclic finite algebra, let `U(a)` be the span of its anchored
box with caps `a`. If `U(a+e_i)=U(a)`, then `U(a)` is stable under
multiplication by `x_i`. Increasing other caps preserves that
stability, because the coordinate operators commute.

It follows that every coordinatewise minimal spanning box satisfies

\[
                       a_1+a_2+a_3\le n-1.              \tag{R9}
\]

Along a coordinate path from `(0,0,0)` to a minimal spanning cap,
every increment must increase rank. A zero increment would make that
coordinate stable from then on; its final increment would be redundant,
contradicting minimality. Rank starts at one and ends at `n`, proving
(R9). There are at most `binom(n+1,2)` such minimal caps: two with
the same first two caps would be comparable, and those two caps have
sum at most `n-1`. This constrains the full saturation diagram, but
does not yet count the intervening fan walls or prove a cubic
polylogarithmic wall-cost budget.

## Why the spatial k-set theorem cannot replace the rank argument

[Sharir--Smorodinsky--Tardos, Theorem 1.1](https://www.renyi.hu/~tardos/michashakhar.pdf)
proves `O(m*k^(3/2))` for ordinary three-dimensional halfspace `k`-sets.
Its proof uses halving triangles and their antipodality. It is not a
general two-parameter matroid bound. Such an extension is false.

[Eppstein, Section 4, Theorems 5 and 7 and Corollaries 1--2](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf)
supplies rank-`r` transversal matroids with `m=Theta(r)` elements and
`Omega(r^(4/3))` one-parameter optimal bases. Restrict to a finite
interval containing all these strict bases, translate the parameter to
be positive, and perturb coefficients to rational numbers while
preserving one strict witness per basis. Clearing denominators and
translating label coordinates makes the labels nonnegative integers.

Take the direct sum of two such matroids. Embed their cost labels as
`(a_i,b_i,0)` and `(c_j,0,d_j)`. On the chart `w_x=1`, the two
independent parameters are `w_y` and `w_z`. Every pair of strict
factor optima is a strict optimum of the sum. Different factor optima
have different cost slopes, so their resulting spatial profiles are
distinct. This gives

\[
  \Omega(r^{8/3})\quad\hbox{positive profile vertices on }
  m'=\Theta(r)\hbox{ elements and rank }2r.              \tag{R10}
\]

Transversal matroids have linear representations over infinite
extensions: put independent variables in the allowed entries of the
bipartite incidence matrix. A minor is nonzero exactly when its
matching exists, because different matchings give different determinant
monomials. The direct sum preserves representability. Thus (R10)
also excludes the purported general `O(m'*rank^(3/2))` extension for
representable matroids. This product is an elementary consequence of
the published lower bound, not a new lower-bound mechanism.

There is no claimed realization of (R10) by a colength-`2r` ideal.
The integer label sizes need not meet the monomial divisor bound,
and the representation need not have commuting cyclic multiplication.
The planar [dyadic theorem](planar_dyadic_bound.md) already prevents
the required planar factors from being realized at linear colength.
The general-matroid obstruction directs the proof back to the rank
and multiplication constraints above, without changing the goal.

## Verification and remaining obligation

Run from the repository root:

    python3 -B scripts/run_wall_rank_filter_checks.py
    ./build.sh
    python3 tests/test_smoke.py

The checker replays complete saved actual quotient fans, computes the
cheap hull and its contracted ranks from actual columns, and checks
minimal spanning boxes by independent dense elimination. It exhausts
all 512 relevant boxes of an explicit nonuniform eight-point quadric
example, independently filters ordinary free prefixes, and checks the
all-field diagonal scanning obstruction. Separate
finite transversal examples verify the direct-sum profile product by
matching, linear representations and complete exact spatial cones.
These finite checks do not reconstruct Eppstein's asymptotic example
or establish novelty. The required collective bound, the unrestricted
sharp growth and the full priority comparison remain open.
