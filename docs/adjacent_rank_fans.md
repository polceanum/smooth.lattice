# Elementary quotient fans, prime digits and lattice gaps

The prime-digit rank recursion preserves actual dependencies. Its block
structure is an iterated uniform-matroid substitution, so its compact
enumerator must be compared with classical tensor-product formulas.
A face lemma for one-dimensional linear quotients, including adjacent
truncations, gives a quantitative consequence
for exact greedy thresholds and the lattice-programming gap fan.
These statements do not resolve the sharp unrestricted $G_3$ problem,
and their exact priority remains unestablished.

## The prime-digit matroid has a classical composition

Let $\mathcal P_{n,M}$ be the column matroid of the truncated Pascal
matrix in [the rank note](oracle_fan_enumeration.md), with $M=p^s$ and
$0<n<M$. Write $h=M/p$, $n=qh+r$, $0\le r<h$. A basis selects $q$
or $q+1$ entries in each fiber $u+h\{0,\ldots,p-1\}$; the fibers
with $q+1$ entries form a basis of $\mathcal P_{r,h}$.

If $r>0$, this identifies $\mathcal P_{n,M}$ with the matroid tensor
product of $\mathcal P_{r,h}$ and the pointed uniform matroid
$U_{q+1,p+1}$. To verify the identification, take a two-sum at each
smaller element $u$ with the distinguished element of its uniform
block. Exactly one of the two component bases contains the gluing
element. A smaller basis containing $u$ therefore takes $q+1$ ordinary
block elements; a smaller basis omitting $u$ takes $q$. These are
precisely all bases in the digit recursion. Equality of basis sets
proves equality of matroids, not an equivalence of the actual codes.

The smaller matroid has no loops, since every Pascal column has first
entry one. Cyclic translation of its labels is an automorphism,
induced by multiplication by $(1+t)^a$. Transitivity and $0<r<h$
exclude coloops. Thus the usual non-loop/non-coloop tensor hypotheses
apply. When $r=0$, use instead the direct sum of $h$ copies of
$U_{q,p}$; no tensor theorem with a looped smaller matroid is assumed.
The cases $n=0,M$ are the zero and free matroids.

For $p=2$ the blocks are the ordinary parallel or series extensions.
Consequently these binary Pascal matroids are direct sums of
series-parallel matroids. For general $p$ the blocks are uniform
matroids, which need not be binary or regular.

[Jaeger--Vertigan--Welsh, Section 4, equations (4.1)--(4.2)](https://gwern.net/doc/math/1990-jaeger.pdf)
states the tensor-product Tutte substitution and identifies stretches
and thickenings. [Moffatt--Noble--Thompson, Corollary 29](https://arxiv.org/pdf/2309.00493)
states and derives the corresponding matroid formula with explicit
non-loop/non-coloop hypotheses. The subset-rank enumerator in (O16)
is a homogeneous form of the Tutte rank enumerator: for rank $n$,

$$
T_{\mathcal P}(x,y)=(x-1)^n
 H_{n,M}(1,y-1,((x-1)(y-1))^{-1}).                       \tag{A1}
$$

This identity is understood after cancellation, directly from the
subset definition. Thus the compact code-weight expression is a
specialized composition of standard matroid invariants. The exact
Pascal/code identification still needs its full coding-literature
comparison; the generic substitution principle is already published.
There is no basis here for claiming a new general weight-enumeration
algorithm or solving a general coding problem.

## A face lemma for two consecutive truncations

Let $\mathcal M$ have rank $r\ge1$ on a finite ground $E$, with arbitrary
labels $a_e\in\mathbb R^d$. For $k=r,r-1$ put

$$
 P_k=\operatorname{conv}\left\{\sum_{e\in B}a_e:
          B\text{ independent in }\mathcal M,\ |B|=k\right\},
 \qquad P_0=\{0\}.
$$

For a cost $w$, let $F_k(w)$ be the minimizing face and
$D_k(w)=\operatorname{span}(F_k(w)-F_k(w))$ its direction space.
Let $t$ be the largest cost level used by greedy selection of an
$r$-basis and put

$$
 s=r-\operatorname{rk}_{\mathcal M}\{e:w\cdot a_e<t\}.
$$

**Theorem 1.** For every $w$, including ties and degenerate projections,

$$
 \begin{cases}
 D_{r-1}(w)\subseteq D_r(w),&s=1,\\
 D_r(w)\subseteq D_{r-1}(w),&s\ge2.
 \end{cases}                                             \tag{A2}
$$

In particular

$$
 \dim(F_r(w)+F_{r-1}(w))=
       \max(\dim F_r(w),\dim F_{r-1}(w)).                 \tag{A3}
$$

**Proof.** First work with indicator-vector polytopes in $\mathbb R^E$.
Order the cost levels through $t$. The standard greedy face description
is a direct sum of matroid minors, one for each level: restrict to the
cumulative prefix and contract the previous prefix. A minimum basis
chooses its prescribed rank increment in every level. In the last
minor the increment is $s$.

If $s=1$, an optimal $(r-1)$-set uses exactly the same lower-level
bases and no last-level element. Its face direction space is therefore
contained in the full-basis face direction space.

If $s\ge2$, the optimal $(r-1)$-face keeps all lower factors and
replaces the last minor by its rank-$(s-1)$ truncation. Every edge
direction of the last full-basis polytope is a difference
$\mathbf1_B-\mathbf1_C$ with $B,C$ bases differing by one element.
They share $s-1\ge1$ elements. Delete a common element from both:
$B',C'$ are independent $(s-1)$-sets with the same difference.
Hence every full-face edge direction occurs in the direction space
of the truncated face. Edge directions span a polytope's direction
space. The unchanged lower factors give the other directions.
This proves the second inclusion. Linear projection by the labels
preserves both inclusions. The minimizing face of a Minkowski sum
is the sum of minimizing faces, so (A3) follows. $\square$

The use of consecutive ranks is necessary. In $U_{3,4}$ take labels
$(1,0,0),(0,1,0),(2,0,0),(0,0,2)$ and $w=(1,1,1)$.
The rank-one and rank-three faces are edges with respective directions
$(1,-1,0)$ and $(2,0,-2)$. Their sum has dimension two. Thus the
analogous statement for ranks one and three is false.

The [flag-matroid Minkowski theorem, Ardila--Benedetti--Doker,
Theorem 5.1 and Definition 5.2](https://fardila.com/Articles/minkowski.pdf),
is an established framework for truncation polytopes. It does not by
itself assert the projected face inclusions (A2); these are proved
above. No novelty is inferred from not finding that exact formulation.

## Every one-dimensional linear quotient has nested face directions

The nesting is more general than truncation. Let columns $c_e$ span
a vector space $V$ over any field, and let $0\ne v\in V$. Compare
their matroid $\mathcal M$ with the matroid $\mathcal N$ of their
images in $V/Kv$, using the same arbitrary real labels $a_e$.
Their ranks differ by one.

**Theorem 2.** For every cost, the direction spaces of the two
minimizing profile faces are nested. Consequently (A3) holds with
$P_{\mathcal M},P_{\mathcal N}$ in place of the consecutive
truncation polytopes.

**Proof.** Let $E_1,\ldots,E_l$ be the ordered cost levels. For every
cumulative prefix $A$,

$$
 \operatorname{rk}_{\mathcal N}(A)
 =\operatorname{rk}_{\mathcal M}(A)
      -\mathbf1_{\{v\in\operatorname{span}(c_e:e\in A)\}}. \tag{A3a}
$$

The indicator changes from zero to one at a single level $j$.
Before that level projection is injective on the entire prefix span.
After it, projection becomes an isomorphism upon contracting the
previous prefix, since the removed line already lies in that span.
Thus the greedy face minors in all levels except $j$ agree.

At level $j$, contract the cheaper columns. The change is precisely
projection of one represented matroid $\mathcal L$ along the nonzero
image $\bar v$. The spans of its connected components form a direct
sum. Decompose $\bar v$ into those component spans.

If only one component has nonzero support, the projection leaves all
other components unchanged and may split that component. Therefore
the component partition of the new matroid refines the old partition.
If at least two components have nonzero support, projection is
injective on each old component span separately. Every old circuit
remains a circuit, so old connected components cannot split. The new
component partition coarsens the old partition.

For any matroid, the direction space of its indicator basis polytope
is exactly the space of vectors whose coordinate sum is zero on
each connected component, including singleton loops and coloops.
One inclusion follows because bases take a fixed rank in each
component. For the reverse inclusion, two elements in a common circuit
give a basis exchange: extend the circuit minus one element to a
basis and exchange either element. These differences span all
zero-sum directions within a component by circuit connectivity.

Partition refinement therefore decreases the direction space, and
coarsening increases it. Adding the unchanged level factors and then
projecting by the labels preserves the corresponding inclusion.
This proves nesting for every cost. $\square$

In particular a new child ray absent from the parent three-coordinate
fan can occur only in the second case: the removed relation must
have support in at least two components of the cost minor where it
first enters the prefix span. This is a necessary condition only.
It identifies a dependency mechanism, not an aggregate bound.

For nested polynomial ideals $J\subset I$ with colengths $n,n-1$,
take the residue columns on the common ground $V_n^3$ and the kernel
line $I/J$. The theorem applies in the original coordinates over
every field. It includes deleting one rational point from a finite
point set. No assumption of generic ranks or diagonal multiplication
matrices is made. A quotient by two dimensions need not have this
property, as the nonadjacent uniform-matroid example above shows.
That example concerns general linear data. A separate actual ideal-chain
counterexample is given below.

## Actual ideal flags: a socle rule, and a limit on local arguments

Every finite algebra admits an actual chain of one-dimensional ideal
quotients after extension to an algebraic closure. Field extension preserves
all ranks of finite monomial column sets and hence all initial staircases.
The commuting coordinate actions have a common eigenvector: take an
eigenspace for the first action, restrict the next action to it, and continue.
Their common eigenline is an ideal. Repeat in the quotient. Thus this
reduction uses polynomial ideals throughout; it does not presume that
arbitrary matroid truncations are algebras.

Let $J\subset I$ have lengths $n,n-1$ and write $I/J=Kv$.
There are scalars $\lambda_i$ with $x_i v=\lambda_i v$.
For a strictly positive weight, filter $A=K[x]/J$ by the spans of
monomials of cost at most $t$. If $v$ first enters at cost $t$, its
nonzero class in $F_tA/F_{<t}A$ is annihilated by every variable:
$x_iv=\lambda_i v$ belongs to $F_tA\subset F_{<t+w_i}A$.
Its associated graded initial form is therefore in the socle. Positivity
is necessary for this argument.

For a fixed global monomial order, let $D$ be the parent standard set and
write the unique normal form

$$
 v=\sum_{a\in D}c_a x^a,\qquad
 a_* = \max_{\prec}\{a:c_a\ne0\}.
$$

The child standard set is exactly $D\setminus\{a_*\}$. Indeed,
$\operatorname{in}(J)\subset\operatorname{in}(I)$, so the child standard
set is a subset of $D$ of size $n-1$. The displayed normal form belongs
to $I$ and excludes $a_*$. This accounts for the whole difference.
In particular $a_*$ is a maximal member of $D$: deleting a nonmaximal
member would fail downward closure. On each open parent weight cone,
the child therefore deletes a most expensive supported monomial, and
only a socle corner of that parent staircase can be deleted. Ties are
resolved by the specified global order. This is an elementary initial-
ideal fact, not a new point-deletion algorithm.

For point ideals the eigenline is the usual separator of the removed
rational point. [Just--Stigler, Section 2 and Lemma 3.4](https://arxiv.org/pdf/0711.3475)
defines separators and gives their unique reduced expression in the
standard basis. That comparison concerns radical rational-point ideals;
it does not provide a count of all weight cones for arbitrary algebras.
[Denham--Schulze--Walther, Definition 2.17(e), Remark 2.18(b) and
Lemma 2.19](https://doi.org/10.1007/s11005-020-01352-3)
supplies the established elementary-quotient realization and the lifting
of a matroid direct sum to a direct sum of its realizing spaces.
Those definitions and the displayed proof of Lemma 2.19 were inspected.
They account for inputs to Theorem 2; its exact projected-fan priority
remains open. No configuration-hypersurface singularity theorem is used.

### Multiplication does not bound the number of merged components

Here is a uniform all-field control on a proposed shortcut. Put

$$
 D_L=\{a\in\mathbb Z_{\ge0}^3:|a|\le L\},\quad
 S_L=\{a:|a|=L\},\quad
 J_L=(x,y,z)^{L+1},\quad
 I_L=J_L+\left(\sum_{a\in S_L}x^a\right),\qquad L\ge1.
$$

All variables kill the added relation. The parent length is
$N=\binom{L+3}{3}$ and the child length is $N-1$, over every field.
At $w=(1,1,1)$, the pivot minor has $m=\binom{L+2}{2}$ singleton
coloop components. The quotient has the unique circuit with support
$S_L$ and merges **all $m$ components**. Thus a single ideal quotient
can merge $\Theta(N^{2/3})$ components, despite commuting actions and
a socle kernel. There is no constant merger-arity bound from those
hypotheses. Nevertheless the child has exactly three initial ideals:
the largest supported monomial is $x^L,y^L$ or $z^L$ according to the
largest variable weight, with the term order resolving ties.
Counting merged components alone is not a count of new cells.

More generally, if $S$ is a subset of the maximal corners of any finite
downset $D$, **every** flag of linear subspaces $R_j\subseteq K^S$
gives a flag of actual ideals $J_D+R_j$. All nonconstant multiples of
the relation rows vanish. The residue matroids split into fixed coloops
$D\setminus S$ and the arbitrary represented quotients $K^S/R_j$.
This elementary construction is already within the repository's
[full-socle moment bound](interleaving_research.md#all-socle-relation-spaces-have-a-small-moment-budget):
its three-variable fans have $O(n^{3/2})$ cells. Consequently local
component behavior alone does not impose the missing global restriction.
A useful argument must also charge the downset size or interactions
between distinct multiplication layers.

### Nonadjacent fans can cross even in an actual ideal composition chain

Let $D$ be the union of the divisor boxes of

$$
 a=x^3,\quad b=y^3,\quad c=z^4,\quad d=x^2y^2.
$$

There are fifteen monomials in $D$ and these four are maximal. With
$C=J_D$ define nested ideals

$$
 J=C+(a-b+c-d),\quad
 H=C+(a-b,c-d),\quad
 I=C+(a-b,a-c,a-d).                                    \tag{A3b}
$$

Their relation spaces have ranks one, two and three over every field,
and $J\subset H\subset I$. Their lengths are $14,13,12$.
All relations lie in the parent monomial socle, so these are actual
polynomial ideals. The three varying residue matroids are respectively
$U_{3,4}$, $U_{1,2}\oplus U_{1,2}$ on the pairs $\{a,b\},\{c,d\}$,
and $U_{1,4}$. The other eleven columns remain fixed coloops.

At $w=(1,1,1)$ the costs of $a,b$ are three and those of $c,d$ are
four. After subtracting the common fixed moment, the minimizing profile
faces are

$$
 F_J=[a+b+c,a+b+d],\quad
 F_H=[a,b]+[c,d],\quad
 F_I=[a,b].                                             \tag{A3c}
$$

Here letters in profile sums denote exponent vectors. The outside
faces are edges of nonparallel directions $(2,2,-4)$ and $(3,-3,0)$;
their sum has dimension two. The middle face is the parallelogram
spanned by these two directions. Both adjacent pairs satisfy Theorem 2,
but the $J$ and $I$ fan edges cross at the strictly positive ray
$(1,1,1)$. This is a new vertex of their common refinement, absent
from both outside diagrams. Complete exact positive fans have respectively
four, four and three cells; the outside overlay has nine cells.

Thus the abstract nonadjacent counterexample above also has an actual,
all-field ideal-chain counterpart. The obstruction survives the eigenline
and socle conditions; pairwise face nesting cannot simply be iterated
over an ideal composition series. It does not refute an aggregate
incidence estimate, nor give a superquadratic lower family. The whole
example belongs to the already bounded full-socle class.

## A linear overlay bound in three cost coordinates

Restrict either pair of profile fans above to $w_i\ge0$, $\sum_iw_i=1$
in dimension three, and retain cells with a strictly positive weight in their
interior. Let $a,b$ be their numbers of maximal cells, and $c$ the
number in the common refinement. Then

$$
 \max(a,b)\le c\le4(a+b)+2.                              \tag{A4}
$$

For the upper bound, two edges of the two planar diagrams cannot
cross transversely in their relative interiors. Such a crossing would
expose an edge in each projected polytope with nonparallel direction
spaces, contrary to Theorems 1--2. More generally the face-dimension
identity says that the refinement
creates no new interior vertex. Every boundary vertex of the
refinement is already a boundary vertex of one input diagram.

A convex subdivision of a triangle into $a$ cells, with unnecessary
collinear degree-two vertices removed, has at most $2a+1$ vertices.
Indeed Euler's formula gives $E=V+a-1$, while all non-corner vertices
have degree at least three and the three corners have degree at least
two. Thus $2E\ge3V-3$. The refinement has at most $2(a+b)+2$ vertices.
A simple planar graph with that many vertices has at most $3V-6$
edges and hence at most $2V-5$ bounded cells. The looser displayed
constant covers the single-cell and flat cases too. The lower bound
holds because a common refinement refines each input partition.

The statement counts profile cells. If different tied bases project
to the same profile, it does not count them separately. The positivity
restriction does not assert a general affine fan on a zero weight face.

There is a necessary condition for the birth of a new interior ray
when rank increases from $r-1$ to $r$. Such a ray exposes a
two-dimensional $F_r(w)$ but an $F_{r-1}(w)$ of dimension at most one.
Equation (A2) forces $s=1$. At this cost the last minor in the greedy
face decomposition has rank one; its nonloop columns are parallel
modulo the span of cheaper columns. Its profile face is the convex
hull of those last-layer labels, added to the lower-level profile
face. Thus only these rank-one last layers can create new rays.
This is a necessary condition, not a bound on how many such events
occur. A collective bound respecting the actual commuting monomial
actions remains a possible route to the unrestricted problem.

One consequence concerns the actual last selected cost. Let
$F_k(w)=\min_{|B|=k}\sum_{e\in B}w\cdot a_e$. Then
$F_r(w)-F_{r-1}(w)$ is the cost of the last greedy element, and is
linear on this common refinement. Its faithful refinement therefore
has $O(a+b)$ cells, rather than a possible quadratic overlay bound.
For a finite polynomial quotient of length at least $k$, the first
$k$ greedy monomials form a downset and lie in $V_k^3$. The
[rational-shadow proof](rational_shadow_bound.md) applies to their
rank-$k$ truncation on this finite ground. Applying it at ranks
$k,k-1$ and then (A4) bounds the refined exact $k$-th-cost function by
$O(k^{23/8}\log^{3/8}(k+1))$ cells, uniformly in the larger quotient
length. The same coarse exponent also follows by applying the shadow
theorem to a direct sum of the two labeled matroids, allowing repeated
labels. The additional content of (A4) is the instance-wise linear
bound in $a+b$, not a new unrestricted power. This is a structural
bound, not a bit-cost or fast unrank implementation claim. It uses
actual ranks, not generic independence.

## A realizable group-algebra interpretation of the gap fan

Let $\Lambda\subseteq\mathbb Z^3$ have full rank and index $m\ge2$,
let $G=\mathbb Z^3/\Lambda$, and work over any field $K$.
The images $g_i$ of the coordinate vectors generate the finite group.
Use the monomial map $x^a\mapsto[g_1]^{a_1}[g_2]^{a_2}[g_3]^{a_3}$
into the group algebra $K[G]$, with kernel $J_\Lambda$.
Set $\eta=\sum_{g\in G}[g]$. Multiplication by every $[g]$ fixes
$\eta$, so $K\eta$ is a one-dimensional ideal in every characteristic.
It is nonzero even when the characteristic divides $m$.
The inverse image $I_\Lambda$ of this ideal has colength $m-1$.

The only dependence among the $m$ group columns modulo $K\eta$ is
that their sum is zero. Any $m-1$ of them are an integral basis,
with determinant $\pm1$ in the reference obtained by eliminating one
class. This proves the rank statement over every field.

For positive $w$, let $b_g(w)$ be the cheapest representative of
class $g$, with a fixed monomial tie rule, and put

$$
 d_g(w)=w\cdot b_g(w),\qquad
 \gamma_\Lambda(w)=\max_{g\in G}d_g(w).
$$

The standard set of $J_\Lambda$ consists of all $m$ representatives.
The standard set of $I_\Lambda$ deletes the most expensive one.
Equivalently the second residue matroid is the rank-$(m-1)$
truncation of the first, on a shared finite ground $V_m^3$.
For their minimum basis sums $F_J,F_I$,

$$
 \gamma_\Lambda(w)=F_J(w)-F_I(w).                         \tag{A5}
$$

On a cone of the group ideal, all representatives are fixed, and
$\gamma$ is the maximum of their linear costs. Its regions of linearity
within that cone are exactly the additional norm-quotient subdivisions.
Thus the fan which retains the group-standard-set label and the gap
label is exactly the common refinement of the two initial fans.
If its cell count is $h_\Lambda$, (A4) gives

$$
 \max(g(J_\Lambda),g(I_\Lambda))\le h_\Lambda
 \le4(g(J_\Lambda)+g(I_\Lambda))+2.                       \tag{A6}
$$

This also holds for noncyclic $G$. In the cyclic case, with the first
coordinate a generator and other residues $u,v$, a concrete quotient is

$$
 I_\Lambda=(1+X+\cdots+X^{m-1},\ Y-X^u,\ Z-X^v),          \tag{A7}
$$

using $0\le u,v<m$. Since $(X-1)(1+\cdots+X^{m-1})=X^m-1$,
this is the asserted quotient of $K[G]$. If $m=p^s$ in characteristic
$p$, the first polynomial is $(X-1)^{m-1}$. This directly connects
the endpoint of the prime-digit family to the same all-fields gap fan.
No prime-specific counting gain follows from that endpoint.

[Hoşten--Sturmfels, Section 5, Theorem 5.1 and Corollary 5.2](https://arxiv.org/pdf/math/0301266)
defines the gap fan by refining the Gröbner fan according to its
piecewise linear gap function. For full-rank $\Lambda$ the continuous
relaxation has value zero at positive cost, giving the function above.
[Aliev, Lemmas 2.1 and 3.1](https://arxiv.org/pdf/1404.6967)
identifies this gap with the directed quotient-lattice diameter and
with $\rho(\Delta_w,\Lambda)-\sum_iw_i$, where
$\Delta_w=\{x\ge0:w\cdot x\le1\}$. Their displayed arguments and
these hypotheses were inspected; no varying-dimension hardness result
is being transferred to three variables.

### The existing corank-one exclusion still applies

The repository already proved a stronger restricted count than a
naive product of group fans and socle corners:
[the one-dimensional binomial-extension theorem, (BE1)--(BE3)](interleaving_research.md#twelfth-feasibility-pass-one-dimensional-quotients-of-binomial-algebras).
For a binomial quotient of length $m$ and an ideal extension with
one-dimensional kernel, it gives

$$
 g(I)=O(N^2\alpha(N)),\qquad
 N=|V_m^3|=O(m\log^2(m+1)),                              \tag{A8}
$$

where $\alpha$ is the inverse Ackermann function. It is an existing
consequence of published normal-form and envelope bounds, not a new
result of this continuation. It covers exactly the norm quotient
above. Combining (A6), (A8) and the existing binomial count gives

$$
 h_\Lambda=O(m^2\log^4(m+1)\,
                   \alpha(m\log^2(m+1))).               \tag{A9}
$$

This is a quantitative transfer to the faithful lattice gap diagram;
the earlier envelope proof already tracks the same representative and
omission labels and yields this bound directly. Equation (A6) supplies
an additional instance-wise comparison, not a new power in (A9).
Its strongest-current priority is unestablished. In particular this
whole corank-one family cannot give a polynomial separation above the
generic corner-cut scale. The gap connection does not reopen a route
already excluded by (A8). Several removed residue directions, growing
codimension, or another genuine dependency mechanism would be needed.
No higher-corank envelope bound or separating family is asserted here.

## A concrete unrestricted proof obligation

The preceding examples prevent deriving a count just from a bounded
number of merged components or from a no-crossing claim for a full ideal
flag. They leave a collective incidence question open. For each strictly
positive primitive interior fan ray $w$, let $d(w)$ be its incident-cell
count and put

$$
 U(I)=\sum_w\left\lceil\sqrt{\|w\|_\infty}\right\rceil d(w).
$$

Changing the high-ray charge in the existing
[weighted-incidence reduction](weighted_facet_incidence.md#the-remaining-weighted-incidence)
gives, for every integer $q\ge1$ and $L_n=\log(n+1)$,

$$
 g(I)=O\left(n^{5/3}L_n+nq^{10/3}L_n+U(I)/q\right).     \tag{A10}
$$

Indeed a cell incident to a ray of height greater than $q^2$ is charged
at least $q$ in $U$. All other nonboundary cells have only rays of
height at most $q^2$ and are covered by the established short-relation
shadow bound with that threshold. Boundary cells use its existing
$O(n^{5/3}L_n)$ bound. Optimizing the integer $q$, with $q=1$ when
the formal optimum is smaller, yields

$$
 g(I)=O\left(n^{5/3}L_n+
      n^{3/13}U(I)^{10/13}L_n^{3/13}\right).             \tag{A11}
$$

Thus a **currently unproved** budget $U(I)=O(n^3L_n^C)$, with fixed
$C\ge0$, would imply the unrestricted bound
$G_3(n)=O(n^{33/13}L_n^{(10C+3)/13})$. This is a precise proof
target, not an established improvement. It still would not establish
sharp growth or, by itself, complete the research goal.

The cubic power in this proposed budget cannot be reduced uniformly:
the actual three-pole curve family in the incidence note has
$\Theta(n^2)$ rays of height $\Omega(n^2)$ along its stated infinite
length subsequence, hence $U=\Omega(n^3)$. Conversely that family's
quadratic fan count and the universal $O(n^2)$ ray-height bound give
$U=O(n^3)$. The earlier first-moment budget is quartic and sharp there;
using the square root avoids mistaking that counterexample for a
refutation of the cubic question. No proof of the question is supplied.
A useful composition-series argument must control this weighted sum,
or give a comparably strong direct cell count, across all layers.

## Comparing incidence, primitive length and local algebra cores

The proposed square-root incidence should not be replaced by a larger
quantity without checking the loss. For the actual undirected walls put
\[
 W(I)=\sum_e\|q_e\|_1,\qquad
 V(I)=\sum_e\nu(q_e).
\]
The primitive mixed direction is counted once per wall piece. Then
\[
                         U(I)\le2W(I)\le2V(I).          \tag{A12}
\]
These are elementary comparisons, not a proof of a cubic budget or a new
unrestricted exponent.

Here is a proof retaining the primitive lattice information. If integer
vectors \(q,r\) each have both signs, and \(a=\|q\|_1,b=\|r\|_1\),
then every coordinate of their cross product has magnitude at most
\[
                         (a-1)(b-1).                   \tag{A13}
\]
For a two-by-two determinant whose two product terms have the same sign,
its magnitude is at most the larger product. Every coordinate of a mixed
integer vector has magnitude at most its one-norm minus one. If the two
products have opposite signs, one of the selected coordinate pairs has
the same sign. That vector has an opposite-sign coordinate outside the
pair, so the sum of the pair's magnitudes is at most its one-norm minus
one. Bound each coordinate of the other vector separately. This proves
(A13), including zero coordinates.

At a compact positive state face with primitive normal \(w\), list its
primitive edge directions cyclically as \(q_1,\ldots,q_d\). Consecutive
directions are independent, and their integer cross product is a nonzero
integer multiple of \(w\). Set \(H=\|w\|_\infty\). By (A13),
\[
 \left\lceil\sqrt H\right\rceil
 \le\left\lceil\sqrt{(\|q_i\|_1-1)(\|q_{i+1}\|_1-1)}\right\rceil
 \le\frac{\|q_i\|_1+\|q_{i+1}\|_1}{2}.
\]
The last inequality uses arithmetic-geometric means and
\(\lceil t\rceil\le t+1\). Summing cyclically bounds this face's
contribution to \(U\) by its total primitive one-norm. A wall has at
most two interior endpoints in the positive weight simplex, hence at
most two such face incidences. This proves the first part of (A12).
Finally \(\nu(q)\ge\|q\|_1\), by expanding the two divisor products,
which proves the second part. The proof also applies when the compact
two-face is the whole lower-dimensional profile polygon.

For comparison, one can form an actual algebra core for each such face.
Take \(J_w=\operatorname{in}_w I\), take the canonical monomial factor
\(f_w\) of its reduced nonmonomial relations, and let
\[
 m_w=\dim_K R/(J_w:f_w),\qquad
                         M(I)=\sum_w m_wd(w).            \tag{A14}
\]
The initial-ideal step is classical:
[Fukuda--Jensen--Thomas, Lemma 2.13 and Corollary 2.14](https://arxiv.org/pdf/math/0509544)
give its reduced basis from the top parts of any adjacent marked basis.
The finite perturbation identity identifies its standard sets with the
original face's incident sets. Applying the
[canonical deflation](canonical_fan_deflation.md) then preserves the
whole local polygon up to fixed profile translation. This is an actual
colon algebra; no arbitrary common vector span is used as a quotient.

The resulting comparisons are
\[
 M(I)\le nW(I),\qquad
 V(I)\le M(I)+O(n^{8/3}\log(n+1)).                       \tag{A15}
\]
For the first, \(m_w\le n\), and summing interior vertex degrees counts
each wall at most twice. Also \(W\ge2|\mathcal E(I)|\). For the
second, every wall touching an interior vertex \(w\) is realized in
that local colon, so its threshold \(\nu(q_e)\le m_w\). Assign each
such wall to one endpoint and sum. Every other wall has both endpoints
on the simplex boundary, so both adjacent cells are boundary-adjacent.
Their simple planar adjacency graph has linearly many edges. The existing
[boundary-cell bound](weighted_facet_incidence.md#three-planes-cover-every-boundary-cell)
is \(O(n^{5/3}\log(n+1))\); multiplying by the maximum cost \(n\)
gives the remainder in (A15).

Consequently a cubic-polylog budget for \(M\) would imply the already
conditional \(12/5\) exponent through \(V\), not only the weaker
\(33/13\) consequence through \(U\). But it is a stronger requirement
than a \(V\) budget, and can overcount badly even when the fan is fixed.

### A sharp linear loss from retaining local core length

Over characteristic zero, for every \(k\ge2\), let
\[
             I_k=(x^2,y^k,z-x-y),\qquad n=2k.             \tag{A16}
\]
The quotient is \(K[x,y]/(x^2,y^k)\). Its entire positive fan is
independent of \(k\) and has five cones. The following table lists the
other minimal heads after eliminating the largest variable.

| Largest variable | Remaining comparison | Other heads |
|---|---|---|
| \(z\) | either | \(x^2,y^k\) |
| \(y\) | \(z>x\) | \(x^2,z^k\) |
| \(y\) | \(x>z\) | \(x^2,xz^{k-1},z^{k+1}\) |
| \(x\) | \(z>y\) | \(y^k,z^2\) |
| \(x\) | \(y>z\) | \(y^2,yz^{k-1},z^{k+1}\) |

For \(y\) largest, use
\((z-x)^k\equiv z^k-kxz^{k-1}\pmod{x^2}\).
For \(x\) largest, use \((z-y)^2\) and
\(kyz^{k-1}-(k-1)z^k\). In all cases \(z^{k+1}\in I_k\), since
\((x+y)^{k+1}\in(x^2,y^k)\). The listed head ideals have colength
\(2k\), so inclusion proves equality with the initial ideals. Once a
variable is eliminated, all remaining relations are homogeneous in two
variables, where comparing the variables determines every comparison
in a fixed degree. This proves completeness for all global orders.

All five cones meet at the unique interior ray \((1,1,1)\). All five
walls have primitive one-norm and threshold two. The essential linear
relation has gcd one; at that ray \(J_w=I_k\) and \(m_w=n\). Thus
\[
       U(I_k)=5,\quad W(I_k)=V(I_k)=10,\quad M(I_k)=5n.   \tag{A17}
\]
The factor \(n\) in (A15) is therefore sharp in order. No uniform
polylogarithmic comparison \(M=O(W\operatorname{polylog}n)\), or
\(M=O(U\operatorname{polylog}n)\), is possible. These examples do not
refute a cubic \(M\) bound; their fan is constant. They show why retaining
the exact local profile polygon can be unnecessarily expensive when the
target counts cones and primitive directions.

The loss also occurs with a full-dimensional original basis-profile
polytope. Replace the last relation in (A16) by \(z-x-y-1\).
Translation preserves all global initial ideals: each extra term of a
translated monomial is a proper divisor. At \(w=(1,1,1)\) its initial
ideal is the original \(I_k\), so (A17) is unchanged. In its ordinary
vector basis, replacing \(1\) by \(z=1+x+y\) is invertible and moves
the profile off the common plane of standard-set moments. The latter
span a polygon, so the enlarged basis-profile polytope is three-dimensional.
This guards against dismissing the example merely because the unshifted
presentation is positively graded.

The next direct target is the **unproved** budget
\(W(I)=O(n^3\log^C(n+1))\). By (A12) it suffices for (A11)'s
\(33/13\) consequence, while a \(V\) budget would give the stronger
\(12/5\) consequence. Neither follows from these comparisons. The known
three-pole \(U=\Omega(n^3)\) family forces the same cubic lower scale
for \(W\). Also the minimum-length single-wall ideals of direction
\((k,k,-k-1)\) have \(V=k^2+3k+1\) but \(W=3k+1\), so replacing
\(V\) by primitive length is a genuine relaxation. The state-polyhedron
dictionary remains the established
[Babson--Onn--Thomas, Theorem 2.4](https://arxiv.org/pdf/math/0207135).
The determinant comparison and five-cone computation are elementary;
exact priority is unassessed and no paper or new extremal exponent is
claimed from this route audit.

## Opposite support faces retain the parent in one elementary quotient

The direction-space nesting in Theorem 2 has a stronger geometric form,
and extends beyond represented matroids. Let \(\mathcal N\) be an
elementary matroid quotient of \(\mathcal M\) on the same ground.
Equivalently \(\delta(A)=r_{\mathcal M}(A)-r_{\mathcal N}(A)\)
is nondecreasing under inclusion, with \(\delta(E)=1\).
This is the standard quotient criterion, for example
[Chen--Fei--Gao--Sun--Zhang, Proposition 2.8](https://alco.centre-mersenne.org/item/10.5802/alco.475.pdf).
Use the same profile notation and fix an arbitrary cost \(w\). Suppose
\(D_{\mathcal M}(w)\) is properly contained in
\(D_{\mathcal N}(w)\). Let \(\lambda\) be a linear functional that
annihilates \(D_{\mathcal M}(w)\) and is nonconstant on
\(F_{\mathcal N}(w)\). Write \(F^-\) and \(F^+\) for the
minimum and maximum \(\lambda\)-faces of that child face.

**Theorem 3.** There exists a translation \(t\) such that
\[
             F_{\mathcal M}(w)+t\ \subseteq\ F^-+F^+.    \tag{A18}
\]
This needs no field representation and allows arbitrary real labels and
linear projections. The two support faces can have smaller dimensions than the
parent face; neither is asserted to contain the parent by itself.

Here is the exact component description behind the assertion. At the
unique cost level where the prefix rank deficit \(\delta\) changes from
zero to one, let \(K\) and \(L\) be the parent and child greedy
minors. The rank criterion shows that \(L\) is an elementary quotient
of \(K\); all other level minors coincide. Decompose \(K\) into
connected components on grounds \(C_i\). Put
\[
 S=\{j:r_L(E'\setminus C_j)=r_K(E'\setminus C_j)\},
 \qquad L_j=(L/(E'\setminus C_j))|C_j,
\]
where \(E'\) is the critical minor's ground. For \(j\in S\), the
matroid \(L_j\) has rank \(r_K(C_j)-1\). Let \(P_i\) be the
projected basis polytope of \(K|C_i\), \(Q_j\) that of \(L_j\), and
\(P_0\) the sum of the unchanged cost-level factors. Then
\[
 \begin{split}
 F_{\mathcal M}(w)&=P_0+\sum_iP_i,\\
 F_{\mathcal N}(w)&=P_0+
    \operatorname{conv}\!\bigcup_{j\in S}
              \left(Q_j+\sum_{i\ne j}P_i\right).        \tag{A19}
 \end{split}
\]
Indeed a child basis is independent in the parent and has one fewer
element. Exactly one old component loses one rank, say \(j\).
Its other elements form a full basis outside \(C_j\), forcing \(j\in S\).
The parent and child restrictions outside \(C_j\) are equal: they are
quotients with the same rank. Contracting that common basis leaves exactly
\(L_j\). Thus all other components choose arbitrary full bases and
component \(j\) chooses a basis of \(L_j\), and conversely these
choices always give a child basis. This proves (A19). In the represented
case, \(S\) is exactly the support of the removed vector in the old
component spans. Neither genericity nor principal truncation is assumed.

This also proves the abstract extension of Theorem 2. If \(S\) has one
element, the child component partition refines the old one. If \(|S|\ge2\),
the child restriction to each \(C_i\) has the same rank as the parent
restriction: otherwise a deficit on \(C_i\) would force \(S\subseteq\{i\}\).
Hence all these restrictions agree, all old circuits remain circuits, and
the component partition coarsens. The direction-space proof therefore
applies without representation as well.

Put \(L=D_{\mathcal M}(w)\). Each connected component's label
differences lie in \(L\), by the connected-component direction-space
description used in Theorem 2. Consequently \(\lambda\) is constant
on every \(P_i\), \(Q_j\) and \(P_0\). The summands indexed by
\(j\) in (A19) occupy parallel affine copies of \(L\). Choose indices
\(a,b\in S\) at the two different extreme \(\lambda\)-levels.
Then \(a\ne b\), and
\[
 \begin{split}
 F^-+F^+
 &\supseteq 2P_0+Q_a+Q_b+
                 \sum_{i\ne a}P_i+\sum_{i\ne b}P_i\\
 &=F_{\mathcal M}(w)+
          \left(P_0+Q_a+Q_b+\sum_{i\ne a,b}P_i\right).
 \end{split}
\]
The last parenthesized polytope is nonempty. Choose any one of its
points as \(t\), proving (A18).

### A parallel edge and a sharp summed lattice-length bound

In three dimensions suppose the parent minimizing face is an edge and
the child minimizing face is two-dimensional. Choose \(\lambda\) on
the child's plane perpendicular to the parent edge. Its opposite
support faces are points or edges parallel to that edge. Hence the child
has at least one parallel edge, and the sum of the lengths of those
parallel edges is at least the parent edge length. For integral labels
and primitive parent direction \(q\), this is the exact integer claim
\[
       \sum_{e\text{ child edge},\ e\parallel q}\ell_q(e)
                     \ \ge\ \ell_q(F_{\mathcal M}(w)),   \tag{A20}
\]
where \(\ell_q\) is the positive multiple of \(q\) in the edge
vector. There are at most two parallel edges in the child polygon.
This strengthens direction-space nesting: the old direction occurs on
the child polygon's boundary, with the stated total length. It gives
no upper bound on the new nonparallel edges or their number of births.

Both one and two inherited edges are necessary, and the bound is sharp
in actual polynomial ideals over every field. Fix \(t\ge1\), use the
monomial ideal complementary to the divisor downset of the displayed
tips, and impose the listed relations in its socle.

* For tips \(a=x^t,b=y^t,c=z^t\), the parent has the relation
  \(a-b\) and the child also has \(a-c\). Their lengths are
  \(3t,3t-1\). At \(w=(1,1,1)\), the parent edge has lattice
  length \(t\) in direction \((1,-1,0)\). The child triangle has
  just one parallel edge, also of length \(t\).
* For tips \(a=x^{2t},b=x^ty^t,c=x^tz^t,d=y^tz^t\), the parent
  relations are \(a-b,c-d\), and the child relations identify all four
  tips. Their lengths are \(3t^2+4t-1,3t^2+4t-2\). The parent
  edge has lattice length \(2t\), while the child parallelogram has
  two parallel edges of length \(t\) each. Thus neither individual
  edge need be as long as the parent edge.

The hypothesis of one removed dimension cannot be dropped. Instead take
\(a=x^tz^t,b=y^tz^t,c=z^{2t},d=x^ty^t\). Impose just \(a-b\)
in the parent and identify all four tips in the child. The lengths are
\(3t^2+4t,3t^2+4t-2\), so this removes two dimensions. At the same
weight the parent edge has direction \((1,-1,0)\), while the child
parallelogram has only directions \((1,0,-1)\) and \((0,1,-1)\).
It has no parallel edge. The functional \((1,1,-2)\) exposes single
points at both extremes, so even (A18) fails. All tips in these examples
are distinct maximal monomials, and the relation coefficients are units;
this proves the assertions in every characteristic. Their complete fans
have respectively \((2,3),(2,4),(2,4)\) cells. None is a large-fan
construction.

There is also a higher-dimensional scope issue in (A18). Two parallel
classes labelled by the segments
\([(0,0,0),(t,0,0)]\) and \([(0,0,h),(0,s,h)]\), with
\(t,s,h>0\), give a planar parent parallelogram. Quotienting their
two column lines by the sum vector gives a child tetrahedron. At cost
zero its two extreme \(z\)-faces are segments, and their Minkowski
sum equals the parent parallelogram. There is no parallel planar facet.
This is projected matroid data, not a three-variable positive-weight
quotient example.

The direct-sum face and flag-Minkowski framework is established; the
inspected [Ardila--Benedetti--Doker, Theorem 5.1](https://fardila.com/Articles/minkowski.pdf)
states the flag-polytope sum. (A18)--(A20) follow from the explicit
one-rank component description above. Principal-truncation base deletion,
for example [Eur, Proposition 3.1.5](https://people.math.harvard.edu/~ceur/pdf/Eur_PhD_Thesis.pdf),
is a related established special case; a general quotient must not
silently be assumed principal or generic. The checker additionally uses a
rank-oracle Vámos matroid with two coloops and its truncation, outside
the represented examples. Its exceptional four-subsets are five of the
six unions of two of four disjoint pairs; the checker
replays every basis exchange and obtains Ingleton values \(15<16\).
The established nonrepresentability of this configuration is also recalled
in [Bernstein--Walsh's discussion of Vámos-like matroids](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v31i1p15/pdf/).
It is not presented as a polynomial-quotient residue matroid. These primary statements were inspected. The
exact projected containment's priority remains unassessed.
The inheritance theorem does not prove an amortized upper bound on \(W\):
new walls can still be born inside a parent cone, and nonparallel edges
can be added at a surviving parent ray. A cubic aggregate needs a bound
on those contributions, rather than another count of merged components.
The unrestricted goal and all paper criteria remain unchanged.

### Primitive perimeter at a ray birth or death

There is also an upper bound on the primitive cost of one such event.
Let the common integer labels have coordinate ranges \(D_i\), and put
\(D=\sum_{i=1}^3D_i\). Consider any matroid quotient that drops rank
by \(c\ge1\), not necessarily one.
At a nonzero cost \(w\), suppose one minimizing profile face is a polygon
and the other is a point or an edge. Their direction spaces then lie in
the common plane \(w^\perp\). For a polygon \(F\) put
\(\pi_1(F)=\sum_{e\text{ edge of }F}\|q_e\|_1\), using primitive
directions and omitting their integer edge multiples. Then
\[
 \begin{array}{ll}
 \pi_1(F)\le2cD,&\text{if the other face is a point},\\
 \pi_1(F)\le(4c+2)D,&\text{if the other face is an edge}.
 \end{array}                                             \tag{A21}
\]
The statement is symmetric in which face is the parent. For the edge
case, the sum of the **actual** one-norm lengths of nonparallel polygon
edges is at most \(4cD\); the two possible parallel edges are counted
primitively in (A21). Their actual lengths need not be bounded by \(D\).

To prove it, use the same strict cost order for parent and child greedy
selection. The rank criterion makes every child greedy rank increment at
most the corresponding parent increment. Thus the child basis is a subset
of the parent basis, with exactly \(c\) elements deleted. This is the
classical greedy-basis containment, also stated in
[Chen--Fei--Gao--Sun--Zhang, Proposition 2.10](https://alco.centre-mersenne.org/item/10.5802/alco.475.pdf).
At each generic secondary objective,
the exposed profiles therefore differ by a sum of \(c\) ground labels.
These label sums have coordinate ranges at most \(cD_i\).

If the smaller face is a point \(p\), all polygon vertices lie in one
translate of those sums or their negatives. A planar convex polygon's
actual one-norm perimeter is twice the sum of its coordinate ranges,
which proves the first bound. Primitive perimeter can only be smaller.

If the smaller face is an edge with direction \(q\), divide the secondary
objective circle in the polygon's plane by the two directions perpendicular
to \(q\). On each open half-circle the smaller face exposes the same
endpoint. The polygon's corresponding boundary chain consequently has
all its vertices in one fixed translate of the label sums. Each coordinate
has at most two monotone phases on this convex chain, so its variation
is at most twice its range. The two chains contain every nonparallel edge
and give the \(4cD\) actual-length bound. There are at most two remaining
edges, both parallel to \(q\). A projected matroid-polytope edge is parallel
to a difference of two ground labels, so its primitive direction satisfies
\(\|q\|_1\le D\). Adding these two primitive contributions proves
(A21). This last standard edge-direction fact follows from basis exchange;
the preimage of a projected edge has an edge projecting nontrivially onto it.

For an actual quotient of lengths \(n,n-1\), take the common universal
monomial ground \(V_n^3\). Then \(D\le3(n-1)\). At every strictly
positive ray which is present in just one of the two fans, the other
minimizing face has dimension at most one. By the local proof of (A12),
\[
 d(w)\left\lceil\sqrt{\|w\|_\infty}\right\rceil
                 \le\pi_1(F(w))\le18(n-1).              \tag{A22}
\]
If the other face is a point the stronger constant \(6(n-1)\) applies.
No radicality, characteristic-zero or generic-rank assumption is used.
The point-face perimeter bound is attained by deleting \((1,1,1)\)
from the two-point scheme \(\{(0,0,0),(1,1,1)\}\): at the equal
weight the parent profile triangle has primitive perimeter six, the child
is a point, and the common ground has \(D=3\). Its total \(W\) falls
from six to zero, so a monotonicity claim for \(W\) would be false.

The linear order in (A22) is sharp in actual radical ideal steps. In the
[three-pole point family](ideal_step_incidence.md), put
\(A=k,B=k+1,C=2k+1\), so \(N=A+B+C=4k+2\).
These three integers are pairwise coprime. Hence
\[
 w=(BC,AC,AB),\qquad H=(k+1)(2k+1)
\]
is primitive. At lengths \(N-1,N-2\), the longer face is a triangle
and the shorter a point; at lengths \(N+1,N\), the longer face is
a point and the shorter a triangle. The common strict prefix has
\(N-2\) monomials, and the three tied candidates are
\(x^A,y^B,z^C\). Selecting respectively one, zero, three or two
of these proves the assertions from the known exact curve fan.
Each triangle has primitive perimeter \(2N\), degree three and
weighted incidence \(3\lceil\sqrt{(k+1)(2k+1)}\rceil\), both
of linear order in the quotient length.

These bounds quantify one creation or disappearance. They do not count
how many events occur over the objective plane, or bound changes at a ray
present on both sides. In particular they do not imply
\(W(A)-W(A/Kv)=O(n^2)\) or the desired cubic total. Establishing an
amortized upper bound still requires controlling those missing contributions.
The [separate ideal-step analysis](ideal_step_incidence.md) now proves
\(\Theta(n^3)\) jumps of both signs for \(U\) along infinite length
subsequences in that fixed curve family. This is consistent with
\(\Theta(n^2)\) new rays at one step, each with only linear cost.
A quadratic per-step recurrence for \(U\) is therefore unavailable.
On the same family \(W\) has the exact quadratic divisor-sum increment
proved in that note; a general quadratic increment for \(W\) is still
an open possibility, not a consequence of (A18)--(A22).

~~~bash
python3 -B scripts/run_adjacent_rank_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The bounded checker compares projected exposed-face spaces by separate
Gaussian calculations, including explicit nongeneric line projections,
exhausts the small positive profile fans and
checks the entire overlay. A nonadjacent-rank counterexample verifies
the abstract scope restriction. A separate actual ideal-chain control
checks multiplication, cyclicity, border relations, complete positive
fans and the crossing in (A3b)--(A3c) in several characteristics.
The simplex socle family checks the unbounded merger arity, including
the unchanged three-cell child fan. Actual group and norm quotients are checked by
integral border relations, determinant-one bases, exhaustive downsets
and independent shortest paths. Uniform block substitutions are checked
against separately enumerated Pascal bases. These controls support the
proofs; they do not establish sharp growth, novelty or publication merit.
