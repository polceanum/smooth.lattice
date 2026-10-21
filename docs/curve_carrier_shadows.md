# Linear shadows for quotients of a fixed affine curve

Let \(C\subset K[x,y,z]\) be a fixed ideal with
\(\dim K[x,y,z]/C\le1\). For every finite-colength ideal \(I\supset C\),
of length \(n\ge1\), the following bounds hold:
\[
 \begin{split}
 T_I(q),\ p_I(q)&=O_C(n) \quad\hbox{for every primitive }q\ne0,\\
 g(I)&=O_C(n^{14/5}),\\
 |\mathcal U(I)|&=O_C(n^{14/5}).
 \end{split}                                                     \tag{C1}
\]
Here \(T_I(q)\) counts open positive initial cones strictly cut by
\(q^\perp\), and \(p_I(q)\) counts actual undirected parallel wall
pieces. The constants are independent of \(n,q\), coefficients of the
additional relations, and their reducedness or locality. They depend
on the carrier data of \(C\), specified below.

This does not require proportional pole divisors, a common pole winner,
or a reduced curve. The full-fan exponent is an upper bound, not a
sharpness assertion. It is weaker than the quadratic bound in the
previous proportional-pole class, but applies to arbitrary finite
quotients of any fixed affine curve. It neither settles unrestricted
\(G_3\) nor rules out superquadratic growth within this larger class.
Exact numerical priority remains unassessed.

## A matroid shadow on finitely many affine lines

Let a rank-\(r\) matroid have a finite set of labels in a union of
\(t\) affine lines in any real vector space. Repeated labels are allowed.
The image of its basis-profile polytope under every linear map to a
plane has
\[
                   s(P)\le4(1+8rt+t)                         \tag{C2}
\]
vertices. Rank zero has one profile separately. The bound is
independent of the lengths, spacing, directions and translations
of the lines. A collapsed image line is a point.

Here is a direct specialization of the established convex-token
argument. Cover each projected line by an open thin strip containing
exactly the ground points on that line. This is possible because the
ground is finite. Cover each collapsed line by a small open rectangle
containing exactly its coincident ground points. The at most \(4t\)
boundary lines avoid all ground points. Assign each point to the
first cover member containing it. Its label still certifies membership
in that affine line or point; overlapping covers cause no problem.

In one objective quadrant, use costs \(\sigma uX+\tau Y\),
\(u>0,\ \sigma,\tau\in\{-1,1\}\). Adjacent greedy exchanges give
\(r\) token paths. Every step has direction
\((-\sigma,\tau u)\), and the successive \(u\)'s are nondecreasing.
For any linear functional, its change has at most one sign reversal
along a token. Each ground-avoiding boundary line is therefore crossed
at most twice by each token.

A label-changing exchange crosses a boundary of the earlier of its
two cover members. There are at most \(8rt\) such exchanges. An exchange
whose endpoints have the same label is parallel to one fixed line.
Nonzero exchanges in that direction can occur at only one objective
parameter in this quadrant. Thus there are at most \(t\) additional
changing parameters. Counting the initial state and the four quadrants
proves (C2), also for segments and points.

Simultaneous equalities can be resolved by adjacent swaps at the same
parameter. The labels and their collinearities must be retained;
perturbing the labels off their lines is unnecessary and would destroy
the hypothesis.

The token construction is classical:
[Eppstein, Section 4, Theorem 4 and its proof](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf)
was inspected. No new general token or parametric-matroid theorem is
claimed. The finite-line specialization above supplies the estimate used
here.

## The finite carrier data of a fixed curve

Let \(J_1,\ldots,J_M\) be the monomial initial ideals of \(C\) for
global monomial orders. There are finitely many. Each has dimension
at most one. Its standard monomials lie on finitely many coordinate
rays and in a fixed finite set.

An elementary decomposition suffices. Choose \(D\ge1\) larger than
every coordinate exponent of every minimal generator of the \(J_j\)'s.
For each \(J_j\), every pair of variables supports some monomial in
\(J_j\), since its dimension is at most one. A standard monomial
therefore has at most one coordinate at least \(D\). The part in
\([0,D-1]^3\) is finite. Outside this box, fix the other two
coordinates: either the entire ray starting at coordinate \(D\)
is standard or none of it is. Indeed every generator exponent in
the varying coordinate is below \(D\).

Consequently each carrier has a cover by at most
\[
                         t=3D^2+D^3                       \tag{C3}
\]
affine lines, allowing one degenerate line for each finite point.
This deliberately loose count is uniform in \(n\). Smaller explicit
ray decompositions may be used.

For each \(j\), restrict the residue-column matroid of \(I\) to
\[
 E_j=\{a\notin J_j:\ \prod_i(a_i+1)\le n\}.                 \tag{C4}
\]
This matroid has rank \(n\): a generic positive weight in the
corresponding \(C\)-cone supplies an \(I\)-standard basis inside it.
Let \(P_j\) be the convex hull of its basis exponent sums.
All its labels are in the fixed line cover, so every planar shadow
of \(P_j\) satisfies (C2) with \(r=n\).

Every positive initial standard set of \(I\) belongs to at least one
\(E_j\). This follows from \(C\subset I\), using a positive generic
weight inside that initial cone and avoiding the walls of \(C\).
Its profile is a vertex of that \(P_j\).
The universal-ground and profile correspondence is the established
[Babson--Onn--Thomas framework, Proposition 2.2 and Theorem 2.4](https://arxiv.org/pdf/math/0207135).

## Counting cuts and parallel pieces

Use the full residue matroid on the ordinary universal ground to
compare support minima. Each \(P_j\) is a subpolytope of its profile
polytope: every one of its bases is an actual monomial vector basis
of \(I\). No algebra quotient or contraction is being substituted.

If an open \(I\)-cone is cut by \(q^\perp\), choose a strictly exposing
positive objective in the intersection. Assign its standard basis to
one carrier \(E_j\) containing it. Its projected profile is a vertex
of the shadow of \(P_j\). Two such strict states cannot have the same
projected profile, since they would tie at every objective in
\(q^\perp\).

For a parallel wall piece, assign one endpoint basis to a carrier.
At a relative interior wall weight the full minimizing face is the
edge parallel to \(q\). The minimizing face of the assigned \(P_j\)
contains that endpoint and is contained in this edge, so it projects
to a vertex. Different parallel edges have different projected
profiles: a fiber over a projected vertex is one segment. This
argument also covers walls that coincide with walls of \(C\);
both endpoint bases need not fit the same carrier.

Thus \(T_I(q)\) and \(p_I(q)\) are bounded by the sum of \(M\) shadow
counts in (C2), proving the first line of (C1). This is a bound on
the positive fan restriction, not a claim that the full universal
profile polytope has only \(O_C(n)\) vertices in every shadow.

## A full-fan bound from lattice incidence

Every coordinate of every \(P_j\) is between zero and \(n(n-1)\).
Add the unit cube to handle flat hulls. The resulting lattice polytope
retains at least as many vertices, has surface area \(O(n^4)\), and
each planar shadow gains at most six vertices.

The short-normal incidence lemma in the
[rational-shadow proof](rational_shadow_bound.md#small-integer-projections-see-all-vertices-with-short-normals)
therefore gives, for \(h\ge1\),
\[
              f_0(P_j)=O(tn h^{3/2}+n^4/h).
\]
Taking \(h=(n^3/t)^{2/5}\) when this is at least one yields
\[
 f_0(P_j)=O(1+tn+t^{2/5}n^{14/5}),\qquad
 g(I)=O\!\left(M(1+tn+t^{2/5}n^{14/5})\right).             \tag{C5}
\]
The incidence lemma and lattice surface estimate are existing inputs;
(C2) removes their previous projection-height loss in this class.
The constants in (C1) are uniform for specified \(M,D\), over every field.

There is a further output-cardinality consequence. An \(I\)-standard
downset in a fixed \(J_j\)-carrier cuts each coordinate ray at most
once. A minimal excluded monomial is either a minimal generator of
\(J_j\), one of those ray frontiers, or a point in the fixed box.
Hence each reduced basis has \(O_C(1)\) polynomials. Taking the union
over (C5) proves the last line of (C1). This counts polynomials,
not coefficient bits or computation time.

## Coordinate tags have only quadratic added order on fixed curves

For the [coordinate tag](projection_tag_fans.md),
\[
 T_{z,1}(I)=(I\cap K[x,y])R+zI,\qquad
 g(T_{z,1}(I))-g(I)=\sum_{q\in Q_{xy}}T_I(q).
\]
Its planar elimination ideal has \(O(n\log(n+1))\) walls.
The first part of (C1) gives
\[
                g(T_{z,1}(I))-g(I)=O_C(n^2\log(n+1)).      \tag{C6}
\]
The same holds in the other coordinates. Thus one coordinate tag
cannot create a polynomial excess above quadratic from a fixed-curve
parent whose own fan has no such excess. A large parent fan is still
possible under (C5); the statement does not bound it quadratically.

## Verification and remaining boundary

Run

~~~bash
python3 -B scripts/run_curve_carrier_shadow_checks.py
~~~

The report is results/research/curve_carrier_shadow_checks.json.
It checks exact strip assignments, adjacent greedy exchanges, every
token boundary crossing, internal collinear events, four full quadrant
sweeps, and independent complete basis-profile hulls. Actual finite
curve quotients are checked by commuting multiplication and complete
border-cone traversal, including nonreduced quotients and finite
characteristics. Their carrier assignments, parallel-edge injections
and coordinate-tag counts are verified separately.
The two fixed curves have parameterizations
\((t,t^{-1},t+t^{-1})\) and \((t,t^2,t^{-1})\).
Their coordinate pole divisors are nonproportional. Integer S-pair
reductions and Laurent identities verify their five and three curve
carriers, and exact cone areas verify full coverage.

The rank/state dictionary, finite Gröbner fans, monomial ray
decomposition, convex tokens, lattice incidence and tag formula are
established ingredients. Precise priority of the finite-line estimate
and its quantitative curve consequences remains unassessed. The
strongest-current unrestricted comparison, middle-cost coexistence
budget, sharp \(G_3\) endpoint and original publication requirements
remain open. Growing carrier complexity is not covered uniformly by
the \(O_C\) notation.
