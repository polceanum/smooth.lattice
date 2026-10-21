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
That counterexample concerns general linear data; it has not been
realized here as a pair of polynomial ideals.

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

~~~bash
python3 -B scripts/run_adjacent_rank_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The bounded checker compares projected exposed-face spaces by separate
Gaussian calculations, including explicit nongeneric line projections,
exhausts the small positive profile fans and
checks the entire overlay. A nonadjacent-rank counterexample verifies
the scope restriction. Actual group and norm quotients are checked by
integral border relations, determinant-one bases, exhaustive downsets
and independent shortest paths. Uniform block substitutions are checked
against separately enumerated Pascal bases. These controls support the
proofs; they do not establish sharp growth, novelty or publication merit.
