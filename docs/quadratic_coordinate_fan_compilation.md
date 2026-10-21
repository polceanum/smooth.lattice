# Joint flags for a quadratic coordinate

Let \(A=K[x,y,z]/I\) have length \(n>0\) over any field, and assume
\(I\) contains a monic univariate polynomial of degree two in \(x\).
There is an exact description of its entire positive initial-ideal fan
by two planar flags and their relative partial matching. The matching
is essential: the [equal-marginal construction](planar_layer_rank_obstruction.md)
shows that the two planar systems alone do not determine this fan.

The description below has \(O(n^2\log(n+1))\) occurrences of labelled
matching carriers as the planar weight varies. It also gives an
enumeration reduction using \(O(n^{\omega+2}(1+\log(n+1)))\) field
operations, with exact rational geometry counted separately. Here
\(\omega>2\) is a feasible matrix multiplication exponent. The output
is all initial standard sets and their full cones; a list of all reduced
Gröbner bases is a different output task.

These are structural and arithmetic statements for this class. They
do not improve the existing \(O(n^{7/3})\) bound for the number of its
cones, prove a cubic collective wall-cost bound, establish novelty, or
resolve unrestricted \(G_3(n)\). In particular a carrier can meet many
other carriers; the number of carriers is not a bound on fan cells.

## The two cyclic planar spaces

Write \(B=K[y,z]\), and inside \(A\) put
\[
 L=B\cdot1,\qquad U=B\cdot x,\qquad
 p=\dim L,\quad q=\dim U.
\]
The quadratic relation gives \(L+U=A\), so
\[
                     h=p+q-n=\dim(L\cap U)\le n.       \tag{Q1}
\]
Neither flatness over the quadratic algebra nor monomial planar kernels
is assumed. The polynomial may split, be inseparable, or be irreducible.
The case \(U=0\) is included. Every global initial standard set has
\(x\)-exponents zero or one, because the relation has leading monomial
\(x^2\) for every global order.

Normalize weights by \(w_z=1\), and put \(u=w_y/w_z>0\),
\(t=w_x/w_z>0\). Scan planar monomials in increasing \(u a+b\).
Compress the lower and upper lists independently by greedy linear
independence, retaining ordered bases
\[
 \mathcal A=(a_1,\ldots,a_p),\qquad
 \mathcal B=(b_1,\ldots,b_q),                             \tag{Q2}
\]
where the labels are full exponent triples with first coordinate zero
and one respectively. A rejected column is a combination of earlier
columns of its own list, so it can never be selected in any merged
greedy scan, for any \(t\). Compression is therefore exact.

Both lists can be obtained from the finite planar ground
\[
       V_n=\{(a,b)\ge0:(a+1)(b+1)\le n\}.               \tag{Q3}
\]
Indeed \(L\) and \(U\) are cyclic \(B\)-modules. Their planar initial
bases are downsets of sizes \(p,q\le n\), and hence lie in \(V_n\).
This use of cyclicity is only inside each layer; the whole quotient
can have arbitrary dependencies between them.

## The relative matching and its lines

For fixed internal orders, set
\[
 R(i,j)=\operatorname{rank}(\mathcal A_{\le i}\cup
                                  \mathcal B_{\le j}),
\quad
 d_{ij}=R(i,j-1)+R(i-1,j)-R(i-1,j-1)-R(i,j).              \tag{Q4}
\]
The [two-list identity](rational_shadow_bound.md#two-lists-give-a-partial-matching)
proves that \(d\) is a zero-one partial matching with \(h\) pairs.
Every merged minimum basis consists of all unmatched elements and the
cheaper member of each matched pair. Its cost is
\[
 \sum_{e\text{ unmatched}}c_e+
                         \sum_{d_{ij}=1}\min(c_{a_i},c_{b_j}). \tag{Q5}
\]
This is classical interleaving structure. See
[Gusfield (1984)](https://www.sciencedirect.com/science/article/pii/0166218X84900775);
the publisher abstract identifies the two-ordered-set problem, while
a theorem-level comparison with the full article remains deferred.
The repository note supplies a proof of the identity used here.

A pair contributes the comparison line
\[
       t=(a_{i,y}-b_{j,y})u+a_{i,z}-b_{j,z}.             \tag{Q6}
\]
Within a region of fixed internal orders, its part with \(t>0\) is an
actual wall carrier. Different pairs may give the same line and switch
together. After identifying coincident lines, the arrangement has at
most \(1+h(h+1)/2\) open cells. Formula (Q5) labels each cell without
another independence test.

There is one useful crossing restriction, but it is not a global count.
If \(i<i'\) and \(j>j'\), then
\[
 c_{a_i}-c_{b_j}<c_{a_{i'}}-c_{b_{j'}}
\]
throughout the interval of internal orders. These two lines cannot
cross there. Crossings of distinct matching lines require both endpoint
orders to increase together. The equal-marginal grid demonstrates
that many such crossings really can occur.

## All joint ranks from one rank-profile matrix

Extend \(\mathcal A\), in its given order, to a basis \(C\) of \(A\)
using columns of \(\mathcal B\). Let \(G\) be the \(n\)-by-\(q\)
coordinate matrix of \(\mathcal B\) in \(C\), and reverse its row order
to get \(H\). Write \(\rho(k,j)\) for the rank of the leading
\(k\)-by-\(j\) submatrix of \(H\). Modding out the first \(i\)
coordinate vectors of \(C\) proves
\[
                         R(i,j)=i+\rho(n-i,j).          \tag{Q7}
\]
Let \(P_H\) be the rank-profile matrix, whose leading submatrices
have the same ranks as those of \(H\). Substituting (Q7) into (Q4)
gives the exact index rule
\[
                         d_{ij}=(P_H)_{n-i+1,j}.        \tag{Q8}
\]
Indices here start at one. Thus the matching is obtained from the
last \(p\) rows of one rank-profile matrix, reflected vertically.
The other pivots account for the part of \(U\) outside \(L\).

Rank-profile matrices and their efficient computation are established
linear algebra. [Dumas--Pernet--Sultan, Theorem 3 and Section 4.2](https://arxiv.org/pdf/1601.01798)
give their defining leading-rank property and elimination framework.
Their [earlier paper, Lemma 2, Corollary 2 and Theorem 2](https://arxiv.org/pdf/1301.4438)
give leading-rank recovery and the \(O(mnr^{\omega-2})\) arithmetic
bound. Those statements and the relevant proof passages were inspected.
No new rank-profile algorithm is claimed. The checker uses elementary
column elimination, independently checked against every leading rank.

## A uniform bound on matching-carrier changes

The coordinates of a label in (Q3) are at most \(n-1\). Consequently
every positive equality between two planar monomial costs occurs at
\[
             u=b/a,\qquad 1\le a,b\le n-1.              \tag{Q9}
\]
The resulting \(O(n^2)\) intervals are a universal refinement of both
list orders, independent of the field and coefficients.

There is a smaller aggregate bound than multiplying these intervals
by \(h\). Count an occurrence of a literal pair \((a,b)\) every time it
appears after being absent on the preceding open interval; include
the pairs in the first interval. Equivalently these are the maximal
consecutive interval runs carrying that pair. Then
\[
       N_{\rm carriers}=O(n^2\log(n+1)).                 \tag{Q10}
\]
Restricting a run to the positive part of (Q6) does not increase its
number of components, since a linear function has at most one positive
interval inside a fixed interval.

To prove (Q10), track the first \(k\) accepted planar labels, for every
\(k\le p\) or \(k\le q\). Such a prefix is a downset of size \(k\):
if a proper divisor was rejected, multiplication of its earlier
dependence would reject the purportedly accepted multiple too. Hence
all its labels lie in \(V_k\). The prefix is precisely the greedy
basis of the rank-\(k\) truncation of the planar column matroid on
\(V_k\). The [planar matroid theorem](planar_dyadic_bound.md) therefore
bounds its number of distinct open-interval states by
\(O(k\log(k+1))\), over every field.

A given prefix state cannot disappear and later reappear: its strict
optimality region in the single parameter \(u\) is an interval, being
an intersection of affine cost inequalities. If \(D\) is the total
number of changed prefix sets across adjacent intervals, this gives
\[
 D=O\left(\sum_{k\le p}k\log(k+1)+
                 \sum_{k\le q}k\log(k+1)\right)
   =O(n^2\log(n+1)).                                    \tag{Q11}
\]
At one simultaneous order event let \(S,T\) be the lower and upper
prefix indices whose sets change. The joint rank \(R(i,j)\) is unchanged
when \(i\notin S\) and \(j\notin T\). Thus an entry of \(d\) can
change only in rows of \(S\cup(S+1)\), or columns of
\(T\cup(T+1)\). The literal label at a row or column can change only
at those same indices. Since each matching uses at most one entry in
each row and column, at most \(4(|S|+|T|)\) new labelled pairs appear.
This deliberately loose constant handles simultaneous changes without
expanding a tie into artificial orders. Initially there are at most
\(n\) pairs, proving (Q10).

The bound concerns runs with their literal monomial endpoints. It does
not count just slopes, tangent parameters, marginal flags, or arbitrary
unrealized comparison directions.

In particular the entire support of actual walls with nonzero
\(x\)-component can be covered by \(O(n^2\log(n+1))\) straight
segments or rays in the positive \((u,t)\)-quadrant, allowing endpoints
on its boundary. Each is contained in one actual wall line. Intersections
with other walls do not split a member of this cover. To see this,
take the positive part of (Q6) on each run and its relative closure.
They are actual wall portions by (Q5). Conversely any nonvertical
actual wall has points outside the finite set of vertical cuts (Q9),
so these closures cover it. Gluing collinear touching portions can only
decrease the cover size. This is a collective geometric consequence;
it still does not count the pieces made by crossings.

## A bound on noncentral exposed two-faces

Let the state polyhedron be the convex hull of basis exponent-sum
profiles plus the nonnegative orthant. Among its compact two-dimensional
faces exposed by strictly positive weights, at most
\[
                         O(n^2\log(n+1))                \tag{Q13}
\]
fail to be centrally symmetric. Equivalently this bounds the positive
fan vertices with a noncentrally symmetric dual polygon. This does not
bound the number of centrally symmetric two-faces.

Here is a local proof that also identifies exactly where exceptions
can occur. At a weight inside an open interval of (Q9), formula (Q5)
expresses the minimum-basis cost as a linear function plus a sum of
minima of two linear functions. The exposed face is a translate of
the Minkowski sum of the tied endpoint segments \([a_i,b_j]\).
It is therefore a zonotope and centrally symmetric.

At a planar order event \(u_0\), consider a strictly positive weight
\((t_0,u_0,1)\) which is not an endpoint of a positive literal
matching run. Every matching pair tied at this weight on either side
of \(u_0\) then persists, with the same monomial endpoints, on both
sides. All other pair choices are locally strict and contribute linear
functions. The minimum-basis cost in a neighborhood consequently has
the form
\[
 F(w)=\ell_\pm(w)+\sum_{(a,b)\in T}\min(w\cdot a,w\cdot b),
                    \qquad \pm(u-u_0)>0,                \tag{Q14}
\]
with the same tied-pair set \(T\) on both sides. Continuity makes
\(\ell_+-\ell_-\) a multiple of \(w_y-u_0w_z\).
Concavity of the minimum-basis cost makes the resulting two-piece
function \(\min(\ell_+,\ell_-)\) (or a single linear function).
Thus the exposed face is still a translate of a sum of segments,
possibly with one additional segment from this vertical event.

Every noncentrally symmetric exposed two-face must therefore occur
at an endpoint of a literal carrier run. A positive run has at most
two finite interior endpoints, so (Q10) proves (Q13). This argument
uses literal pairs, not just coincident geometric lines; changing the
monomial endpoints on the same line can change the dual polygon.
At a positive normal a full two-face has a unique normalized weight,
so different faces cannot be charged to the same point.

This also gives a necessary condition directly on the original cone
count. If \(Z(I)\) is the number of centrally symmetric compact
two-faces exposed by strictly positive weights, then
\[
                         g(I)\le2Z(I)+O(n^2\log(n+1)).  \tag{Q15}
\]
For completeness, the number of noncompact facets of the state
polyhedron is \(O(n\log(n+1))\). A facet normal is nonnegative;
one with exactly one zero coordinate is an edge normal of a planar
projection of the basis-profile polytope. On the finite layer ground,
each such projected label lies in the planar hyperbolic cross of
parameter \(2n\): the \(y,z\) projection has parameter \(n\), and
the \(x,y\) and \(x,z\) projections have parameter at most \(2n\).
Repeated projected labels are allowed. The planar matroid theorem
therefore bounds the shadow edges by \(O(n\log(n+1))\). Normals
with two zero coordinates add at most three facets.

Truncate only the unbounded part of the full-dimensional state
polyhedron by a sufficiently distant total-coordinate plane. Its
vertices include all \(g(I)\) initial profiles; its facets are precisely
the retained original facets and at most one new cut facet. The
three-dimensional Euler inequality \(v\le2f-4\), together with
(Q13) and the noncompact facet bound, proves (Q15).
Thus any polynomial separation above the generic scale in this class
would require \(\Omega(g(I))\) centrally symmetric compact two-faces.
More noncentral junctions alone cannot supply it. The matching identity
explains the restriction, but does not bound how many independent
crossings coexist. Its precise priority has not been established.

## Arithmetic compilation of the full fan

Suppose commuting \(n\)-by-\(n\) multiplication matrices and a cyclic
vector are supplied. One marked Gröbner basis can instead supply these
data using the existing [enumeration preparation](oracle_fan_enumeration.md).
Computing a first basis from arbitrary generators is additional work.

Precompute the \(2|V_n|=O(n(1+\log(n+1)))\) layer columns, using
\(O(n^3(1+\log(n+1)))\) field operations with dense matrix-vector
multiplication. For each interval of (Q9):

1. Sort the planar labels at a rational sample \(u\).
2. Extract each ordered independent list by column-rank-profile
   elimination on its \(n\)-by-\(|V_n|\) matrix.
3. Form \(C,G,H\) and use (Q8) to recover the matching.
4. Arrange its lines in the interval and positive \(t\)-half-plane,
   and label each cell by (Q5).

Steps 2--3 cost \(O(n^\omega(1+\log(n+1)))\) field operations per
interval, by the cited rectangular rank-profile bound and linear
solving. There are \(O(n^2)\) intervals, so the total field cost is
\[
                O(n^{\omega+2}(1+\log(n+1))).            \tag{Q12}
\]
Line arrangements, comparisons, labelling and merging patches use no
more operations in \(K\). They require separate exact rational work
and output storage. A direct clipping implementation is polynomial;
(Q12) is not a bound on its combined running time or bit complexity.
It does not assert an implementation of fast matrix multiplication.

Every full positive cone meets an open interval in (Q9). The patches
with a given basis label cover that cone apart from subdivision edges.
Its closure is convex, so taking the convex hull of their vertices
recovers the full closed cone. This also removes artificial vertical
cuts. No selected-weight sampling is used as evidence of fan coverage.

With the existing \(O(n^{7/3})\) class count, direct field elimination
once per cone has envelope \(O(n^{\omega+7/3}(1+\log(n+1)))\).
(Q12) improves that field-only envelope for the stated output.
This comparison is to the repository's general reduction. No
fastest-known or new-publication claim is made.

## What this leaves unresolved

The representation keeps the extension information which marginal
planar data lose, and bounds its total carrier changes. It still leaves
intersections between matching lines with increasing endpoint orders,
as well as the gluing at planar order events. Controlling those
collectively in actual quotients is the missing step. Neither (Q10)
nor the per-interval arrangement bound yields a quadratic full-fan
bound. The existing integer-moment slice ceiling \(O(n^{7/3})\) is
stronger than counting all possible intersections of these carriers.

There is no implication here for prime distribution, smooth-number
random access, unrestricted \(G_3\), or a supergeneric lower family.
The original goal remains active and its success criterion is unchanged.

## Exact checks

Run

~~~bash
python3 -B scripts/run_quadratic_coordinate_fan_checks.py
~~~

The report is `results/research/quadratic_coordinate_fan_checks.json`.
It checks rank-profile leading ranks by independent dense row
elimination; compares (Q7)--(Q8) with every joint prefix rank; checks
the matching basis against direct quotient greedy selection in every
arrangement cell; and compares the reconstructed complete cones with
an independent traversal using border normal forms. It also checks
prefix non-reappearance and the carrier-birth charge, including
simultaneous equalities. Fixtures include nonflat and flat square-zero
extensions, a zero upper layer, distinct two-plane points, and split
and irreducible quadratic algebras in several characteristics.
It recovers every positive exposed two-face from the independently
traversed full fan, tests central symmetry on its exact profile polygon,
and charges every noncentral face to a literal matching-run endpoint.
Finite checks support the proofs; the checker implements elementary
verification algorithms, not the optimized arithmetic bound (Q12).
