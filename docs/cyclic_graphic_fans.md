# Periodic direct products: a forest model and a uniform fan bound

For any field \(K\), integers \(p,q\ge1\), and residues
\(0\le a<p\), \(0\le b<q\), consider the actual coordinate algebra
\[
 A=K[u]/(u^p-1)\times K[v]/(v^q-1),\qquad
 x=(u,v),\quad y=(u^a,0),\quad z=(0,v^b).                 \tag{G1}
\]
Its length is \(n=p+q\), and \(yz=0\). These coordinates generate
the algebra: \(y^p=(1,0)\), \(z^q=(0,1)\), and their products with
powers of \(x\) give both monomial bases. No distinct-root assumption
is needed. The model includes repeated-root characteristics.

**Theorem.** If \(I\) is the kernel of (G1), then, uniformly in all four
parameters and the field,
\[
                    g(I)=O\bigl(n^2\log^4(n+1)\bigr).      \tag{G2}
\]
The entire positive fan is independent of the field. A field-independent
arrangement of \(O(n\log^2(n+1))\) integer comparison planes refines it.
Thus this nongeneric, zero-overlap family cannot give a polynomial
separation above the generic corner-cut scale. The logarithmic exponent
is not asserted sharp.

This is a consequence of the forest reduction below and the published
fixed-dimensional integer-hull theorem. Incidence matroids, Kruskal's
algorithm, integer-hull vertex bounds and hyperplane-arrangement counting
are classical. Exact priority of the combined quantitative statement
remains unassessed. It neither settles unrestricted \(G_3(n)\) nor
completes the research goal.

## Exact incidence columns, over every field

In the bases \(u^r\), \(v^s\), a monomial has column
\[
 x^iy^jz^k\longmapsto
 \begin{cases}
 0,&j,k>0,\\
 e_{i+aj\bmod p},&j>0,k=0,\\
 f_{i+bk\bmod q},&k>0,j=0,\\
 e_{i\bmod p}+f_{i\bmod q},&j=k=0.
 \end{cases}                                             \tag{G3}
\]
Negate the second block of rows and add a ground vertex. The axis
column is the incidence column of an edge from its \(p\)-vertex to
its \(q\)-vertex. Each branch column is an edge from one such vertex
to ground, up to a column sign. This is an oriented graphic matroid,
including in characteristic two. Every subset has rank
\((n+1)-c\), where \(c\) counts graph components including ground.
Leaf elimination gives determinants \(\pm1\) for tree bases, proving
field independence without diagonalizing the group algebras.

Put \(g=\gcd(p,q)\), \(A_0=p+q-g\). The minimal annihilator of
\(x\) is
\[
 \operatorname{lcm}(t^p-1,t^q-1)
                 =(t^p-1)(t^q-1)/(t^g-1).                \tag{G4}
\]
The polynomial gcd identity follows from the Euclidean algorithm over
every field, including inseparable cases. Hence the first \(A_0\) axis
columns are independent, and every later one is a combination of
earlier, strictly cheaper axis columns. Only the axis edges
\[
              (i\bmod p,i\bmod q),\quad0\le i<A_0        \tag{G5}
\]
are needed. They form a forest with \(g\) components. Their weights,
normalized by \(w_x>0\), are the fixed distinct integers \(i\).

For branch minima it suffices to use \(0\le i<p,1\le j\le p\) on
\(Y\), and \(0\le i<q,1\le k\le q\) on \(Z\). Reducing a larger
exponent by a period preserves its column and decreases its positive
cost. These finite grounds contain all actual standard sets after the
dominated axis columns are removed. Greedy is exactly a minimum
spanning tree on the forest plus its labelled ground edges.
Every monomial initial ideal has a generic positive-weight representative
on this finite ground, so counting its open chambers suffices.
Alternatively the ordinary length-\(n\) hyperbolic ground has
\(O(n\log(n+1))\) nonzero columns here, so one initial standard set
can be found by sorting and graph greedy in
\(O(n\log^2(n+1))\) comparisons. This is an application of classical
graphic greedy, not a new general spanning-tree algorithm or bit bound.

## A general counting lemma for a forest with ground edges

Take any forest on \(N\) vertices. Its edge weights are distinct fixed
multiples \(\tau_e w_x\); the ground edges have distinct coefficient
vectors \(v\in\mathbb R^3\) and weights \(v\cdot w\).
Assume also that no ground coefficient vector equals any
\(\tau_e e_x\), excluding persistent ties between the two edge types.
The two independent normalized parameters are \(w_y/w_x,w_z/w_x\).
Process forest edges by increasing \(\tau_e\), recording the singleton
components and each component created by a merge. There are at most
\(2N\) such components, forming a merge hierarchy.

For a hierarchy component \(C\), let \(Q_C\) be the convex hull of
all coefficient vectors of ground edges incident to its vertices, and
write \(h_C\) for its number of vertices. If none are incident, put
\(h_C=0\) and \(m_C=+\infty\). Then the number of uniquely
optimal spanning trees, at generic positive weights, is
\[
                       O\left((N+\sum_C h_C)^2\right).    \tag{G6}
\]
This concerns the specified forest-plus-ground structure; arbitrary
bounded-treewidth graphs are not being assigned this bound.

To prove (G6), the minimum ground cost \(m_C(w)\) is determined by
the vertices of \(Q_C\). Its normal-fan boundaries are supported on
\(O(1+h_C)\) planes: a three-dimensional polytope has linearly many
edges, and polygons, segments and points have the same bound.
Include those planes and, for every vertex \(v\) of \(Q_C\), the
planes \((v-\tau e_x)\cdot w=0\) for the birth and death levels of
\(C\), when finite. There are \(O(N+\sum_C h_C)\) planes in total.

On an open arrangement chamber each component's unique minimum edge
is fixed, as is the comparison of its cost with its birth and death.
A ground edge is selected precisely when it is that minimum for the
component of the forest at its own cost, i.e. when the minimum cost is
strictly between the component's birth and death. Initial singleton
births and final root deaths are infinite in the appropriate direction.
For a forest edge at level \(\tau\), its two child components are
already connected through ground exactly when both minima are below
\(\tau\); the forest edge is rejected in that case and accepted otherwise.
These are ordinary Kruskal connectivity tests, and they determine the
whole tree on a chamber. An arrangement of \(H\) lines in the positive
projective weight triangle has \(O(1+H^2)\) chambers. This proves (G6).
The description also covers roots in different forest components.

More generally, suppose each \(Q_C\) is the convex hull of a constant
number of affine images of integer points of bounded rational polytopes
in at most \(D\) variables, with a fixed number of inequalities, each
of binary encoding length at most \(\varphi\). The established
integer-hull theorem gives
\[
 h_C=O_D(1+\varphi^{D-1}),\qquad
 \#\text{trees}=O_D\bigl(N^2(1+\varphi)^{2D-2}\bigr).     \tag{G7}
\]
A vertex of a projected integer hull comes from a vertex upstairs;
the hull of a fixed number of projected sets has no more vertices
than their total. Constants may depend on \(D\) and the inequality
count. This is a counting consequence, not an implementation claim
for a general integer-hull enumeration routine.

The exact published input is
[Cook--Hartmann--Kannan--McDiarmid, Theorem 2.1, p. 29](https://www.math.uwaterloo.ca/~bico/papers/ippoints.pdf):
for a rational polyhedron in \(D\) variables described by at most \(m\)
inequalities of size at most \(\varphi\), the integer hull has at most
\(2m^D(6D^2\varphi)^{D-1}\) vertices. The complete relevant pages
27--29 were inspected. Equalities can be written as two inequalities;
full dimensionality is not a hypothesis. No new integer-hull theorem
is asserted here.

## Arithmetic descriptions of all axis components

Assume \(p\le q\); exchanging the two branches covers the other case.
Let \(F_\ell\) contain the first \(\ell\) edges of (G5).
It is enough to describe components at their birth: their vertex sets
remain fixed until the next merge.

For \(\ell<q\), every nonisolated component has one \(Y\)-vertex
\(r\) and exactly the \(Z\)-vertices
\[
                   r+tp<\ell,\quad t\ge0.                \tag{G8}
\]
The other \(Z\)-vertices are singletons. This follows because an axis
edge with index below \(q\) is the only such edge at its \(Z\)-vertex.

For \(\ell\ge q\), all \(Z\)-vertices are attached. Contract those
first \(q\) edges. Each remaining edge, indexed \(q+j\), joins the
\(Y\)-vertices \(j\) and \(j+q\bmod p\), where
\(0\le j<\ell-q\le p-g\). The full permutation \(r\mapsto r+q\)
has \(g\) cycles. At least one edge from each cycle is absent, so its
selected components are paths. Thus a component's \(Y\)-vertices are
\[
 U_C=\{r+tq\bmod p:0\le t<L\},                          \tag{G9}
\]
and its \(Z\)-vertices are exactly
\(\{s:0\le s<q,\ s\bmod p\in U_C\}\).
An entire permutation cycle may occur as a vertex set with one edge
removed; its missing edge fixes a valid path start in (G9).

## Fixed-dimensional lifts of the ground-label sets

The following descriptions have binary size \(O(\log(n+1))\) per
inequality. The number of inequalities is fixed.
For a component of type (G9), its \(Y\)-labels \((i,j,0)\) are the
projections of integer solutions
\[
 0\le i<p,\quad1\le j\le p,\quad0\le t<L,\qquad
              i+aj=r+tq+lp.                              \tag{G10}
\]
Eliminating \(i\), whose equality coefficient is one, makes this a
bounded polytope in three integer variables \(j,t,l\), with affine
projection \((r+tq+lp-aj,j,0)\).
The \(Z\)-labels \((i,0,j)\) are projections of
\[
 \begin{gathered}
 0\le i<q,\quad1\le j\le q,\quad0\le t<L,\\
 0\le r+tq+lp<q,\qquad
              i+bj=r+tq+lp+mq,                           \tag{G11}
 \end{gathered}
\]
Its apparent five variables are unnecessary. Put \(h=t+m\). There
exists an integer \(t\) with \(0\le t<L\) and
\(0\le r+tq+lp<q\) precisely when
\(-(L-1)q\le r+lp<q\): the intervals
\([-tq,(1-t)q)\), for \(0\le t<L\), partition that interval.
The equality then gives the equivalent three-variable description
\[
 \begin{gathered}
 1\le j\le q,\qquad -(L-1)q\le r+lp<q,\\
 0\le r+lp+hq-bj<q,\qquad j,l,h\in\mathbb Z,
 \end{gathered}                                             \tag{G12}
\]
with affine projection \((r+lp+hq-bj,0,j)\). Conversely recover
\(t=-\lfloor(r+lp)/q\rfloor\) and \(m=h-t\) to obtain (G11).
The bounds on \(j\), then \(l\), then \(h\) show boundedness.
The largest right-hand side is \(O(n^2)\), still with binary size
\(O(\log(n+1))\).

Strict integer bounds are replaced by non-strict bounds one smaller.
These eliminations preserve integer witnesses; taking a real projection
before the integer hull would be invalid.

For (G8), the \(Y\)-membership condition is a single congruence;
the \(Z\)-condition is \(i+bj=r+tp+mq\), with a bounded interval
for \(t\). Eliminating \(i\) gives at most three integer variables.
Singleton components use one congruence and two integer variables.
Empty branch sets contribute no candidates.
Thus every hierarchy component satisfies (G7) with \(D=3\) and
\(\varphi=O(\log(n+1))\). There are at most \(2n\) components.
Equations (G6)--(G7) prove (G2) and the stated comparison inventory.

This proof bounds the whole family, rather than inferring growth from
finite fans. It is independent of field characteristic, slope size and
the ratio of the periods. It retains original exponent labels and does
not apply a fan-preserving linear change of coordinates.

## Literature boundary, checks and remaining target

[Eppstein's 2023 paper, abstract and Section 6](https://link.springer.com/article/10.1007/s00453-022-01024-9)
gives single-parameter lower bounds even on bounded-treewidth graphs.
[Fernández-Baca--Slutzki, publisher abstract](https://doi.org/10.1006/jagm.1994.1019)
also concerns one parameter. Neither inspected statement supplies a
two-parameter bound for these periodic ground labels. Their scope is
recorded without claiming an exhaustive novelty audit. The earlier
[cographic construction](cographic_order_ideals.md) instead has rank
\((p-1)(q-1)\) and complements of trees; its classical duality is not
being renamed as a new result here.

Run

```bash
python3 -B scripts/run_cyclic_graphic_checks.py
./build.sh
python3 tests/test_smoke.py
```

The report saves commuting cyclic quotient witnesses, independent
incidence ranks and complete graph-versus-quotient cone comparisons,
arithmetic component descriptions, integer-lift membership checks and
the comparison inventory. The checker may enumerate small integer
sets and arrangements; it is not an optimized integer-hull package.
The uniform logarithmic bound relies on (G7), not on measured counts.

The primitive wall spectrum, unrestricted collective budget and sharp
\(G_3\) target remain open. The integer-hull/forest argument is reusable
under its stated hypotheses, but does not cover arbitrary nongeneric
two-plane ranks, general cographic models or forests with varying
internal edge order. Exact priority and a consequential publication
result still require further work. No manuscript or completed-goal
claim follows from this class exclusion alone.
