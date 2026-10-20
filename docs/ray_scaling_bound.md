# Positive ray scaling, tropical counts and the remaining mixed terms

This note identifies a general reason that some component constructions
cannot multiply their fan counts. The proof uses classical matroid greedy
invariance and a tropical polynomial of bounded degree. It gives an exact
count and an enumeration reduction for a structured class, rather than a
new unrestricted bound for $G_3(n)$.

## A rank bound independent of the number of elements

Let $M$ be a finite rank-$r$ matroid. Assign each element a type
$i(e)\in\{1,\ldots,t\}$ and a positive coefficient $\alpha_e$.
For parameters $q_1,\ldots,q_t>0$ give it cost

$$
 c_e(q)=\alpha_e q_{i(e)}.
$$

**Proposition.** Across all positive parameters, the number of bases
that are uniquely optimal is at most

$$
 \binom{r+t-1}{t-1}. \tag{1}
$$

The same bound holds for distinct strictly exposed basis profiles
$\sum_{e\in B}\alpha_e e_{i(e)}$, allowing repeated element weights.
It holds for arbitrary matroids, not just representable ones. Rank zero
has one empty basis. Arbitrary positive real coefficients are permitted.

To prove this, recall that a strictly increasing transformation of
element costs preserves the entire family of minimum matroid bases.
Indeed that family is determined by the ordered cost levels: it consists
of bases taking the full rank of every cumulative cost prefix. Greedy
selection attains all those prefix ranks simultaneously, and the
successive positive cost gaps prove their necessity. The argument also
handles ties.

Apply the transformation $\log$. Put $u_i=\log q_i$ and
$k_i(B)=|\{e\in B:i(e)=i\}|$. The same bases minimize

$$
 \sum_{e\in B}\log c_e(q)
   =\log\prod_{e\in B}\alpha_e+\sum_i k_i(B)u_i. \tag{2}
$$

For a fixed count vector $k$, the product-minimizing basis is independent
of $u$. Thus two distinct uniquely optimal bases cannot have the same
count vector. For profiles, if two representatives have that vector,
both must attain the same minimum coefficient product. At a weight
strictly exposing either profile, (2) and greedy invariance make both
bases optimal, forcing their profiles to coincide. Nonnegative count
vectors of total $r$ are the compositions counted in (1).

Equivalently, define

$$
 \beta_k=\min_{\substack{B\text{ basis}\\k(B)=k}}
               \log\prod_{e\in B}\alpha_e,\qquad
 T(u)=\min_k(\beta_k+k\cdot u), \tag{3}
$$

omitting infeasible $k$. This is a min-tropical polynomial of degree $r$
in $t$ variables. Its full-dimensional domains have different monomial
labels. Formula (1) is the elementary degree/count ceiling, together
with the greedy reduction (2).

The bound is sharp. Use the uniform rank-$r$ matroid on $t$ groups,
each with coefficients $1,\ldots,r$. For any composition $k$ of $r$,
choose $q_i=2/(2k_i+1)$. Precisely the first $k_i$ elements of group
$i$ have cost below one, and every remaining element has cost above
one. Their union is the unique optimal basis.

## Common sums retain the same ceiling

Suppose several matroids, of ranks $r_1,\ldots,r_s$, use the same
$t$ positive ray parameters, with completely arbitrary coefficients
in each component. Their common positive fan, or equivalently the
positive fan of the Minkowski sum of their profile polytopes, has at most

$$
 \binom{r_1+\cdots+r_s+t-1}{t-1} \tag{4}
$$

strict states. Apply (1) to their direct sum. Adding a fixed translation
to each component's profiles changes no fan. The same reasoning permits
elements of component $j$ at positions $v_j+\alpha_e v_{i(e)}$:
its common offset contributes the fixed linear term $r_j w\cdot v_j$.

Here the rays $v_i$ must be common to all components, and only weights
with $w\cdot v_i>0$ are counted. Restricting these parameters to ones
arising from a lower-dimensional weight space can only reduce the count.
Independent common-ray constructions therefore share one rank ceiling;
their separate numbers of states cannot be multiplied.

## Bounded local certificates and traversal

A minimum basis $B$ is certified by its feasible single exchanges:

$$
 \alpha_e q_i\le\alpha_f q_j
 \quad(e\in B,\ f\notin B,\ B-e+f\text{ a basis},
       \ i=i(e),\ j=i(f)). \tag{5}
$$

Necessity is immediate. For sufficiency, if a cost prefix $A$ has
$|B\cap A|<\operatorname{rank}(A)$, choose $f\in A$ outside the
span of $B\cap A$. Its fundamental circuit in $B+f$ contains
$e\in B\setminus A$. Exchanging that $e$ for $f$ strictly improves
cost. Thus (5) forces every prefix rank to be full, which is the
classical greedy optimality criterion. For each ordered pair $i\ne j$
only the smallest feasible ratio $\alpha_f/\alpha_e$ matters:

$$
 q_i/q_j\le
 \min_{\substack{e,f\text{ as in (5)}\\i(e)=i,\ i(f)=j}}
       \alpha_f/\alpha_e. \tag{6}
$$

Same-type comparisons are constant conditions. Consequently a strict
positive basis region has at most $t(t-1)$ interior facets. Its image
under $u_i=\log q_i$ is a polyhedron of difference constraints
$u_i-u_j\le \log(\alpha_f/\alpha_e)$. The logarithmic map is a
diffeomorphism here; it is not a linear transformation of profile space.

For fixed $t$, (1) and (6) give a traversal reduction using at most

$$
 1+t(t-1)\binom{r+t-1}{t-1} \tag{7}
$$

greedy support calls. Start at any strict positive parameter, compute
(6), and cross every interior facet. Positive full-dimensional basis
regions form a connected adjacency graph; a generic path between
witnesses avoids intersections of codimension two. Repeated states
are discarded; a fixed tie order handles persistent duplicate elements
without counting them as different profiles. At most $r(m-r)$ exchange tests construct (6) per new
basis on $m$ elements, in addition to the greedy call's independence
tests. Thus the support-call bound does not hide a claim that inspecting
the ground is free.

For rational coefficients, cones and facet-crossing witnesses can be
handled using their rational linear inequalities in $q$. Numerical
logarithms are unnecessary. The checker implements this traversal for
three rays and independently compares every returned cone with the
inequalities against all bases. This is a research implementation for
the specified cost model, not a general initial-fan enumerator.

## Exact consequence for axis-carrier quotient algebras

Let $I\subset K[x_1,\ldots,x_d]$ have positive finite colength $n$.
Suppose every positive initial standard set contains only $1$ and
pure powers of individual variables. Its actual residue matroid,
restricted to those pure powers and then contracted by $1$, has
rank $n-1$. The remaining costs are $j w_i$, so (1) gives

$$
 g(I)\le\binom{n+d-2}{d-1}. \tag{8}
$$

As usual, every global monomial initial ideal has a strict positive
weight representative for its finite reduced-basis comparisons.
No assumption about arbitrary residue-subset ranks or generic
sampling is required.

A sufficient algebraic condition is that for every $i\ne j$ the ideal
contains a relation

$$
 x_i x_j=a_{ij}x_i+b_{ij}x_j+c_{ij}. \tag{9}
$$

Every mixed product then reduces in total degree to pure powers.
For positive weights each displayed tail is strictly cheaper than
the head. Thus mixed monomials cannot enter a positive greedy basis,
and the pure-axis restriction spans the quotient.

The bound is sharp even for radical rational point ideals. Choose
distinct rational poles $a_2,\ldots,a_d$ and $n$ distinct nodes $b_l$
away from them, and evaluate

$$
 x_1=t,\qquad x_i=(t-a_i)^{-1}\quad(2\le i\le d)
$$

at those nodes. Partial fractions give (9). For any composition
$k_1+\cdots+k_d=n-1$, the $n$ functions
$1,x_i,\ldots,x_i^{k_i}$ are independent: multiply a putative relation
by $\prod_{i\ge2}(t-a_i)^{k_i}$. Its numerator has degree at most
$n-1$, so vanishing at all $n$ nodes makes it zero; uniqueness of
partial fractions makes every coefficient zero. The weights used in
the sharpness argument above realize every such allocation.
In three variables, (8) is exactly $\binom{n+1}{2}$.

This is a general explanation of the earlier three-pole count.
It does not amplify the known quadratic lower family.

## Where published tropical geometry already applies

For a direct sum of rank-one matroids with all $t$ types in each
component, (3) is a sum of tropical linear forms. Its walls form a
tropical hyperplane arrangement. Count vectors are exactly its coarse
types. This is the published injectivity and sharp generic count in
[Dochtermann--Joswig--Sanyal, Theorem 2.9 and Corollary 2.10,
arXiv:1001.0237v2, printed p. 7](https://arxiv.org/pdf/1001.0237v2).
Their Proposition 2.7 gives the monomial-label interpretation used
in that proof. Those statements and proof, with the definitions in
Section 2, were inspected.

The arbitrary-matroid statement above follows by greedy invariance
and the same tropical monomial count; the inspected tropical theorem
is specifically about arrangements. A complete classification of
its further matroid specializations is not asserted.
[Eppstein's primary exposition, Section 1.2 and the opening of
Section 4](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf)
records the greedy framework and exchange viewpoint. These are
classical inputs. This reduction and its elementary enumeration
consequence are not being promoted as a new standalone paper.

## The actual mixed-monomial gap

The ray assumption is essential. In a general exponent lattice the
cost is $a_1w_1+\cdots+a_dw_d$. Its logarithm is
$\log(\sum_i a_i e^{u_i})$, which is not an affine function of $u$.
The bounded-degree polynomial (3) in only $d$ variables does not follow.

There is an actual small quotient showing why retaining just pure-axis
counts loses information. Over $\mathbb F_{101}$ let

$$
 I=(x^2,y^3,z-F),\qquad
 F=7x+7y+10y^2+15xy+23xy^2.
$$

Its length is six. Two strict standard sets are

$$
 \{1,x,y,z,z^2,yz\}\quad\text{at }(20,15,13),\qquad
 \{1,x,y,z,z^2,xz\}\quad\text{at }(31,51,38).
$$

Both have pure-axis counts $(1,1,2)$, but their profiles are
$(1,2,4)$ and $(2,1,4)$. Independent determinant, direct greedy and
strict minimal-border checks revalidate these previously recorded
states. This is a failure of that proposed certificate, not a
superquadratic family or a disproof of the desired extremal bound.

The unrestricted problem therefore needs control of the mixed linear
forms and their interactions, beyond common-ray scaling. Further
enumeration of common-ray components cannot supply the requested
power separation. The main goal, the rational-shadow gap and all
publication criteria remain active.

## The refuted uniform relative-margin question

**Current status:** the uniform hypothesis below is false, even for
actual local quotients over every field. The
[sharp stability theorem](initial_cone_stability.md) gives worst margins
of order $n^{-2}$ in two variables and $n^{-3}$ in three, with constant
initial-fan counts. The conditional packing implication below remains
correct; it cannot be applied to all quotients. The original formulation
and finite diagnostics are retained to identify the failed premise.

There is a more direct candidate for the unrestricted problem. For an
open positive initial cone $C$, let $Z$ be its reduced-basis comparison
vectors, so $C=\{w>0:z\cdot w>0\ (z\in Z)\}$. Define

$$
 \mu(C)=\sup_{w>0}\min_{z\in Z}
                  \frac{z\cdot w}{|z|\cdot w}, \tag{10}
$$

where absolute value is coordinatewise. If there are no comparisons,
the cone is the positive orthant and may be assigned margin one.

**Refuted uniform hypothesis.** There is an absolute $\kappa>0$
such that every positive initial cone of every length-$n$
three-variable quotient over every field has
$\mu(C)\ge\kappa/n$.

This was a proposed proof target, not a theorem inferred from finite checks.
It asks for a robust witness somewhere in the cone, rather than at
its smallest integer witness.

**Conditional consequence.** That hypothesis would imply

$$
 G_3(n)=O(n^2\log^2(n+1)). \tag{11}
$$

Here is a proof of the implication, including the aspect cutoff.
Decrease $\kappa$ to a fixed rational number at most one. Put
$s=\kappa/(2n)$. The hypothesis makes the eroded rational cone

$$
 w>0,\qquad (z-s|z|)\cdot w>0 \quad(z\in Z) \tag{12}
$$

nonempty. Every minimal initial generator has total degree at most
$n$, since all its proper divisors are standard; standard monomials
have degree at most $n-1$. Thus each integer comparison vector has
coordinates bounded by $n$. Clearing the fixed denominator in (12)
gives coefficients $O_\kappa(n^2)$.

The classical cofactor argument now provides an integer interior
witness with all coordinates between one and $O_\kappa(n^4)$:
take the sum of three independent primitive extreme rays of the
closed pointed eroded cone, each bounded by its two-by-two cofactors.
No nonzero defining inequality can vanish on all three independent
rays. The sum is therefore strict, including in every coordinate.
This is the same elementary cone construction used in the earlier
[integer-height comparison](initial_cone_height.md), now applied to
the eroded system.

For that witness the logarithmic aspect center
$(\log(w_x/w_z),\log(w_y/w_z))$ lies in a box of radius
$O_\kappa(\log(n+1))$. Write $P=z^+\cdot w$ and
$N=z^-\cdot w$. Coordinatewise multiplicative perturbations by at
most $e^\rho$ preserve the comparison whenever
$e^{-\rho}P>e^\rho N$, equivalently
$(P-N)/(P+N)>\tanh\rho$. Choose $\rho=s/2$.
Thus every cone contains a logarithmic square of radius
$\kappa/(4n)$ centered in the stated box. Squares from different
initial cones are disjoint. Area packing proves (11).

The known generic lower scale would then match (11). The uniform
hypothesis has now been refuted. The conditional implication is not an
improvement to the established unrestricted bound.

The old hard integer-weight examples do not refute this hypothesis.
For the wholly interior nine-cone family of length $n=2k+1$, its
hard cone is given by
$s>t$, $kt>u$, $u+t>ks$. The witness

$$
 (s,t,u)=(4k+1,4k,4k^2-2k)
$$

has the three relative slacks

$$
 \frac1{8k+1},\qquad\frac1{4k-1},\qquad\frac1{8k+3},
$$

each at least $1/(4n)$. Its large necessary integer height and
projection radius therefore coexist with the proposed robustness
scale. This is an all-$k$ identity, not a numerical rate fit.

The rank and hyperbolic lattice alone are insufficient. For an integer
$M\ge2$, put $n=10(M+5)$ and take four rank-one parallel pairs
with weight vectors

$$
 ((M+1)e_x,Me_y),\quad ((M+2)e_x,(M+1)e_y),\quad
 ((M+4)e_x,(M+3)e_z),\quad ((M+5)e_x,(M+4)e_z).
$$

Add $n-4$ coloops at any other distinct pure-axis points of
$V_n^3$. This is a representable rank-$n$ matroid on distinct
hyperbolic-$n$ weights. One positive cell has

$$
 \frac{M}{M+1}<\frac{w_x}{w_y}<\frac{M+1}{M+2},\qquad
 \frac{M+3}{M+4}<\frac{w_x}{w_z}<\frac{M+4}{M+5}. \tag{13}
$$

In logarithmic aspect coordinates its area is at most
$[M(M+2)(M+3)(M+5)]^{-1}=O(n^{-4})$, using
$\log(1+x)\le x$. Its maximum relative margin is $O(n^{-2})$,
so it violates the proposed scale. However its selected basis
contains $x^{M+2}$ and omits $x^{M+1}$; it is not a downset.
This is not a quotient realization or an extremal lower family.
Any proof of the hypothesis must use additional downset or
multiplication constraints.

The checker also erodes all 98 previously certified nonlinear
rectangle cones by $s=1/(16n)$ and finds strict witnesses. This is
only a feasibility diagnostic on those saved models. The subsequent
actual quotient counterexample shows why these finite successes were
insufficient. Any continuation must use an aggregate count that permits
thin cones, or another global mechanism; further common-ray screens
cannot repair the refuted premise.

## Reproduction

~~~bash
python3 -B scripts/run_ray_scaling_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The report records exact witnesses and relevant source hashes. Existing
research reports, manuscripts, the smooth-number kernels and the paused
ML work are preserved.
