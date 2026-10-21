# Mixed rank controls shadows and strengthens the fan bound

Let \(A=K[x,y,z]/I\) have positive length \(n\), over any field. Fix
any center \((a,b,c)\in K^3\), and put
\[
 T=((x-a)(y-b),(x-a)(z-c),(y-b)(z-c))A,\qquad
 t=\dim_KT,\quad \tau=1+t,\quad L=1+\log(t+2).
\]
Use \(g(I)\), the reduced-basis union \(\mathcal U(I)\), and total
primitive wall cost \(V(I)\) as in the
[mixed-ideal rank theorem](mixed_ideal_rank_fans.md).
There are uniform bounds
\[
 \boxed{\quad g(I)=O(\tau^{10/3}L^4),\qquad
       |\mathcal U(I)|=O(\tau^4L^4).\quad}                 \tag{R1}
\]
In particular an actual family with \(g(I_n)\ge c n^{2+\epsilon}\)
must satisfy, at every center,
\[
 t=\Omega_{c,\epsilon}\left(
       \frac{n^{3/5+3\epsilon/10}}{\log^{6/5}(n+2)}\right). \tag{R2}
\]
The center of smallest mixed rank can be used. Large mixed rank is
necessary in (R2), not sufficient for a large fan.

For a nonzero primitive integer \(q\), let \(T_I(q)\) count open
positive cones strictly cut by \(q^\perp\), and let \(p_I(q)\) count
actual undirected wall pieces parallel to \(q\). Then
\[
 T_I(q),p_I(q)=O\left(\min\left\{
       \tau^{4/3}L^2,
       \tau(1+\|q\|_1^{1/3})L\right\}\right).             \tag{R3}
\]
These estimates are independent of \(n\). A further bound, retaining
the actual exponent heights, is
\[
 g(I)=O\left(1+n^{5/4}\tau^{13/8}L^{3/8}\right).          \tag{R4}
\]
Take the minimum with (R1) and the existing unrestricted estimates.
At \(t\) proportional to \(n\), (R4) recovers the earlier
\(n^{23/8}\log^{3/8}(n+1)\) exponent. None of these bounds improves
that unrestricted exponent or establishes a polynomial lower
separation. Sharp dependence on \(t\) and exact priority remain open.
The full cost consequence remains \(V(I)\le3n\,g(I)\).
The bound on \(\mathcal U(I)\) counts polynomials. Their degrees,
term counts and coefficient bits can still grow with \(n\).

## A contraction of small rank

Write \(u=x-a,v=y-b,w=z-c\), using exponent labels in the original
coordinates for the polynomial columns \(Q_\alpha=u^{\alpha_1}
v^{\alpha_2}w^{\alpha_3}\). Their positive-prefix and downset
properties are proved in the preceding theorem. No equality of
arbitrary-subset monomial ranks is used.

Greedily compress the four carriers, namely \(1\) and the three
positive pure-power rays, in their fixed positive internal orders.
Let \(C=C_0\sqcup C_1\sqcup C_2\sqcup C_3\) be their union, and
\[
 E=\bigcup_{i<j}(e_i+e_j+H_t),\qquad
 H_t=\{\alpha\ge0:\prod_i(\alpha_i+1)\le t\}.
\]
Discard globally zero polynomial columns if desired. Set
\[
 m_C=|C|,\quad q_0=|E|,\quad b_0=\operatorname{rank}(E)\le t.
\]
The matroid on \(M=C\sqcup E\) has rank \(n\); it contains every
positive optimal basis. The actual-ideal argument gives
\[
       n-t\le m_C\le n+3+3t,\qquad q_0=O(tL^2).          \tag{R5}
\]

The dual restriction to the carrier has rank exactly
\[
 h=\operatorname{rank}_{M^*}(C)=m_C-n+b_0
                         \le3+4t.                      \tag{R6}
\]
Thus descending dual greedy compression of each \(C_i\) retains
at most \(h\) elements, even though the full dual rank may be
larger. Let \(F\subset C\) be the union of these four retained lists.
Every strictly positive optimal dual basis is contained in \(F\cup E\).
Indeed a rejected carrier element is in the dual span of earlier,
strictly more expensive elements of its own group. The usual greedy
rejection rule applies when the groups are interleaved.

It follows that
\[
                   P=C\setminus F                       \tag{R7}
\]
is contained in every strictly positive optimal primal basis.
In particular \(P\) is independent. Contract it in the represented
matroid, retaining the ground \(K=F\cup E\). The resulting matroid
\(N\) has
\[
 \begin{split}
 |K|&\le q_0+4h=O(\tau L^2),\\
 r:=\operatorname{rank}N
     &=n-m_C+|F|=b_0-h+|F|
                       \le b_0+3h\le9+13t.              \tag{R8}
 \end{split}
\]
This improves the rank of the encoding to \(O(\tau)\), with no
extra logarithms. It is one fixed matroid for the whole positive
weight cone.

Its minimum profile function, plus \(w\cdot\sum_{\alpha\in P}\alpha\),
is exactly the original minimum there. Strict optima agree by
construction; continuity gives agreement also on positive walls.
Consequently its positive fan is the original fan, with a fixed
profile translation.

This is a vector-matroid contraction. The span of \(P\) need not be
an algebra ideal. It may contain \(1\) while remaining proper.
The matroid is not asserted to be the residue matroid of a
length-\(r\) algebra; its other-objective bases may fail to be
downsets. Its labels can have size proportional to \(n\), even
at fixed \(t\).

For example
\[
 A=K[x,y,z]/(xy,z,x^{\ell+1},y^{m+1},x^\ell-y^m),
                       \qquad \ell,m\ge2
\]
has \(t=0\), length \(n=\ell+m\), and a rank-one encoding on the
two labels \((\ell,0,0),(0,m,0)\). It has two positive cones and
one wall of cost \(n/\gcd(\ell,m)\), hence cost \(n\) for coprime
axis lengths. The compulsory span contains \(1\) and
is proper, so it is not an algebra ideal. This also rules out
an \(n\)-independent full-cost bound in terms of \(t\).

## Two bounds for every shadow

Let \(\Pi\) be the full profile polytope of \(N\). The published
parametric-matroid estimate is \(O(mr^{1/3})\) for a rank-\(r\)
matroid on \(m\) planar labels. Four objective quadrants and (R8)
therefore give
\[
                         s_q(\Pi)=O(\tau^{4/3}L^2).      \tag{R9}
\]
Rank zero has one profile separately. Degeneracies do not increase
the number of exposed profiles; strict exposing objectives survive
a sufficiently small generic perturbation.

For the height-sensitive bound, the body labels satisfy
\[
       \prod_i(\alpha_i+1)\le4t\qquad(\alpha\in E).       \tag{R10}
\]
The other labels lie on the three coordinate axes. Apply the dyadic
projected parallelogram cover from the
[rational-shadow proof](rational_shadow_bound.md), with \(H=4t\),
to the body labels. Cover each projected axis by a thin open strip
containing exactly the ground points on that line; use a small open
rectangle if the axis collapses to a point. Put these at most three
strip or point covers first in the assignment order.

There are \(O(L)\) ground-avoiding boundary lines in this combined
cover. The \(r\) convex token paths cross each at most twice.
Internal exchanges within an axis label occur at only one parameter
per axis and objective quadrant. For each body parallelogram, the
same two coordinate-variation budgets and lattice-direction estimate
as in the rational-shadow proof apply, including when a token leaves
the label and returns. The summed internal budget is unchanged.
For a rank-two integer projection with primitive kernel \(q\) this
gives
\[
 s_q(\Pi)=O\left(L\left[1+r+
                  (r^2t\|q\|_1)^{1/3}\right]\right)
          =O\left(\tau(1+\|q\|_1^{1/3})L\right).         \tag{R11}
\]
At \(t=0\), there are only axis labels and the finite-line token
argument gives a constant. The strips must be chosen relative to
the entire finite projected ground, including all body points.
An axis point outside the body is not assigned a fictitious
hyperbolic height \(O(t)\).

A strictly cut positive cone injects into a projected profile vertex.
A wall parallel to \(q\) is the segment fiber over one such vertex;
different parallel edges give different vertices. These are the
existing strict-exposure and collapsed-edge fiber arguments.
Combining (R9)--(R11) proves (R3).

## Counting directions rather than all arrangement intersections

Every positive wall direction is parallel to the difference of two
ground labels of \(N\). This is the classical edge-direction property
of matroid base polytopes and their linear images.

There are only \(O(\tau)\) retained carrier labels and \(O(tL^2)\)
body labels. Carrier--carrier and carrier--body pairs give
\(O(\tau^2L^2)\) possible directions.

Body--body pairs have a smaller geometric inventory than their raw
pair count. If \(\alpha,\beta\in E\) and
\(q=\operatorname{primitive}(\alpha-\beta)\) has both signs, then
\[
 \nu(q)=\prod_i(q_i^++1)+\prod_i(q_i^-+1)-2\le8t-2.       \tag{R12}
\]
Indeed \(q^+\le\alpha\), \(q^-\le\beta\) coordinatewise, and use
(R10). Every mixed-sign direction in three coordinates can be
oriented as \((a,-b,-c)\), up to permutation, with \(a\ge1\),
\(b,c\ge0\). The condition \(\nu(q)\le H\) implies
\(a+(b+1)(c+1)\le H+1\). Overcounting signs, permutations and
nonprimitive vectors gives
\[
 \#\{q\text{ up to sign}:\nu(q)\le H\}
        =O((H+1)^2(1+\log(H+2))).                       \tag{R13}
\]
The ordinary divisor sum proves this: there are at most
\((H+1)\sum_{uv\le H+1}1\) triples. Thus all possible positive
wall directions together number \(O(\tau^2L^2)\).

For each direction, (R9) bounds the number of actual wall pieces.
The adjacency graph of a convex subdivision of the positive
triangle is connected, so \(g(I)\le1+\#\{\text{interior walls}\}\).
Summing (R9) over this direction inventory proves (R1)'s fan bound.
The earlier \(O(\tau^{2/3})\) minimal-generator bound proves its
reduced-basis union bound. Rearranging, with \(t\le n\), proves (R2).

## Retaining exponent height

Every carrier label has coordinates at most \(n-1\). Every body
label has coordinates at most \(t+1\le n+1\). Hence \(\Pi\) lies
in a box of side \(O(n\tau)\), by (R8). Add the unit cube to
handle flat polytopes. Its surface area is \(O(n^2\tau^2)\), and
each shadow gains at most six vertices.

The established short-normal lattice-incidence lemma, with (R11),
now gives for \(h_0\ge1\)
\[
 g(I)=O\left(\tau L h_0^{5/3}
                         +n^2\tau^2/h_0\right).
\]
Choose \(h_0=\max(1,(n^2\tau/L)^{3/8})\) to obtain (R4).
This argument retains \(n\). Treating the remaining long axis
labels as a universal ground of length \(O(t)\) would be invalid.

## Verification and exact boundary

Run

~~~bash
python3 -B scripts/run_mixed_rank_shadow_checks.py
~~~

The report is results/research/mixed_rank_shadow_checks.json.
It checks the dual carrier-restriction rank, fixed compulsory columns,
an independent represented contraction, complete original and contracted
fans, every tested shadow quadrant, token boundary and internal
direction budgets, and the three-part direction inventory.
Long pure axes at fixed mixed rank test the height limitation.
Nonideal compulsory spans and non-downset other-objective bases are
explicit guards. Finite checks do not prove asymptotic sharpness.

The exponent \(1/3\) in (R9) is published, not new:
[Dey, Section 4.1 and Theorem 4.1](https://link.springer.com/content/pdf/10.1007/PL00009354.pdf),
together with
[Eppstein, Section 4, Theorem 4](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf).
Their relevant polygon argument and matroid reduction were inspected.
Matroid contraction, the projected-edge dictionary, convex tokens,
divisor counting and lattice incidence are established tools. The
additional inputs here are the small-rank carrier contraction and its
mixed-ideal parameterization. Exact numerical priority of the combined
bounds remains unassessed. No general spatial k-set upper bound for
represented matroids is assumed.

The original unrestricted sharp \(G_3\) or infinite actual polynomial
separation, a consequential theorem at that scope, strongest-current
exact-priority comparison and honest-paper requirements remain open.
