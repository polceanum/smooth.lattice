# Exact bridge barriers for mixed-hook price sectors

The annihilator assumption can be removed from the low/middle-price
partner bounds. A fixed own-arm prefix interacts with the bridge space
through at most \(t^2+1\) critical flag changes. Between these barriers,
the full original rank test is one opposite-arm span threshold.

For arbitrary mixed constant-axis central hooks over every field, this
bounds the entire positive fan restricted to
\[
                    w_x<t\max(w_y,w_z).                 \tag{BS1}
\]
No unit-prefix regularity is required. The exact barrier and weighted
partner statements below concern zero outer images and the displayed
sector. A stronger coarse full-fan count already follows by combining
the saved bridge and two-colour theorems, including arbitrary outer
images. Thus the sector count is not an asymptotic improvement.
The unrestricted \(G_3\) problem and exact priority remain open.

## Family, partner bounds and sector count

Use
\[
 I=(x^2,y^a,z^a,y^tz^t-xf),\quad a>2t,\quad t\ge2,\quad
 n=4at-2t^2,\quad A=a-t.
\]
The two arm units are \(\alpha+N_Y,\beta+N_Z\), with nonzero constants
\(\alpha,\beta\) and arbitrary \(N_i\in(v)R_A\).
Both outer generator images are zero.

Normalize \(w=(\lambda,u,1)\), and first take the convex wedge
\[
                    \mathcal R=\{u>1,\ 0<\lambda<tu\}.  \tag{BS2}
\]
Fix a pure Y label \(E=(\epsilon,Y,Z)\), with \(Y\ge(3+\epsilon)t\), where
\(t\le Z<2t\) for \(\epsilon=0\), and \(0\le Z<t\) for \(\epsilon=1\).
Its opposite pure Z partners have the other first exponent.
Let \(P\) count distinct partners on genuine open cross-layer wall
pieces. Then
\[
 \begin{array}{ll}
 \epsilon=0,\ 0<\lambda<t:&P\le8t^6,\\
 \epsilon=0,\ t<\lambda<tu:&P\le32at^6/Y,\\
 \epsilon=1,\ 0<\lambda<tu:&P\le32at^6/Y.
 \end{array}                                            \tag{BS3}
\]
Thus in the whole wedge either type of fixed heavy-axis label satisfies
\(YP=O(at^6)\). There is no primal or dual annihilator hypothesis.

The distinct cross-layer directions in \(\mathcal R\) number
\[
                         O(at^7(1+\log(a+1))).           \tag{BS4}
\]
With the established parallel-shadow estimate, their actual edge
contribution is \(O(n^2t^6(1+\log(n+1))^2)\), retaining subdivisions
and recurrent pieces.

Let \(g_{\rm low}(I)\) count positive initial cones whose interiors meet
(BS1). Then
\[
               g_{\rm low}(I)
                 =O(n^2t^6(1+\log(n+1))^2).             \tag{BS5}
\]
This is a whole restricted-fan count. It is weaker than the following
full-family corollary; its purpose here is to verify that the exact
barrier representation also accounts for same-layer boundaries.

## Quantitative comparison with the full hook family

For every compatible flat hook map \(\phi:M\to K[y,z]/M\), where
\(M=(y^a,y^tz^t,z^a)\), the earlier reductions imply
\[
                         g(I_\phi)=O(n^2t^6).           \tag{BS0}
\]
This includes arbitrary compatible images of both outer generators,
all price sectors, zero coefficients and every field. No constant-axis
or unit hypothesis is needed for this comparison.

Indeed, the [bridge reduction](hook_bridge_reduction.md) gives
\(O((at)^2)=O(n^2)\) convex internal-order chambers. On each chamber,
the full original optimum is recovered from one fixed dual matroid
of rank \(k=3t^2\) on at most \(3k=9t^2\) original labels.
Its labels still have first exponent zero or one. For this fixed
matroid the [two-colour theorem](two_colour_threshold_braids.md)
applies with overlap \(h\le k\) and within-colour pair budget
\(P\le\binom{3k}{2}\). Its bound
\(O((h+1)(P+1)+h^2)\) is \(O(k^3)=O(t^6)\).
The dual maximization is covered by negating all costs and replacing
\(\lambda\) by \(-\lambda\); the theorem allows arbitrary real
coefficients and the whole real parameter plane. Restricting a
strict optimal region to a convex chamber cannot split it. Summing
these patch bounds proves (BS0).

This is an explicit corollary of saved results using classical greedy
compression; no new novelty is asserted for the combination. The
conditional dual matrix can change between chambers. The corollary
does not provide a single coefficient-independent sampler of the same
cardinality, or the weighted individual-label estimates in (BS3).
Retain the established quadratic-coordinate \(O(n^{7/3})\) bound and
the box bound \(O(n^{3/2}a)\) whenever smaller. A useful new count must
beat these full-fan ceilings, rather than only extend (BS5) to the
high-price sector.

## A bounded intersection flag

Fix \(E\), and one own-arm order sector \(J\subset(1,\infty)\).
As in the [weighted-partner proof](hook_weighted_partners.md), at most
\(t^2\) sectors suffice: their cuts are \(b/c\), with
\(1\le c\le b\le t-1\). Both own-layer orders are fixed in \(J\);
opposite-arm orders may change.

Contract the earlier same-layer pure Y columns. Their span is fixed
in \(J\). Work in the resulting Y quotient \(V_Y\), keeping Z separately.
All these contractions use genuinely cheaper original columns.

Let \(H\) be the \(t^2\) central bridge columns. Their Y projections,
together with the image of \(E\), span a subspace
\[
                           B\subset V_Y,\qquad
                           d=\dim B\le t^2+1.            \tag{BS6}
\]
Let \(F_j\) be the span of the first \(j\) opposite-layer pure Y
columns. These are upper Y identities when \(\epsilon=0\), and
lower Y sources when \(\epsilon=1\). The flag
\[
                              B\cap F_j
\]
has at most \(d\) strict changes. Its change indices are the critical
barriers. Record their original labels and planar costs \(h_i(u)\).
The \(d+1\) possible intersection states include the initial empty
prefix; the last state need not have dimension \(d\).

This retains full linear information. If \(j_i\le j<j_{i+1}\),
then \(B\cap F_j=B\cap F_{j_i}\). Every linear combination of \(E\),
bridges and opposite-arm columns has Y component in \(B\).
It lies in \(F_j\oplus0\) exactly when it lies in
\(F_{j_i}\oplus0\). Thus contracting the full flag prefix or only
the prefix through its last critical label gives identical linear
relations on all these columns. This remains true after contracting
the bridges. It is not a support-matching argument or an assertion
that the vector quotient is an algebra.

One way to compute the barriers is the exact formula
\[
 \dim(B\cap F_j)=\dim B+\dim F_j-\dim(B+F_j).
\]
Each added flag column increases this intersection dimension by at
most one. No enumeration of possible coefficient specializations is used.

## Exact threshold model inside one barrier state

All central bridges are cheaper than \(E\) in (BS2) under the stated
thresholds: \(Y\ge3t\) for a lower label and \(Y\ge4t\) for an upper label.
Fix a barrier state \(i\), and contract:

- the earlier same-layer pure Y columns;
- the opposite-layer Y prefix through its last critical label;
- every central bridge.

Denote the image of \(E\) by \(e_i\).
The remaining relevant same-layer Z columns form an increasing flag
\(C(u)\): a label with planar coordinates \((r,s)\) becomes cheaper
than \(E\) when \(u>(s-Z)/(Y-r)\).
There are \(t\) short-coordinate values, all below \(Y\).

The opposite-layer Z columns \(D\) have costs \(Ru+S\).
Define \(b_i(u)\) as their least cost threshold spanning \(e_i\)
together with \(C(u)\). Unreachable and already-loop cases give no
exchange. The threshold's selected death column is the possible
opposite partner. For \(u_2>u_1\),
\[
 b_i(u_2)\le b_i(u_1)+\gamma(u_2-u_1),\qquad
 \gamma=\begin{cases}t-1,&\epsilon=0,\\2t-1,&\epsilon=1.
 \end{cases}                                            \tag{BS7}
\]
An earlier spanning representation remains available because \(C\)
only grows and every \(D\) cost has slope at most \(\gamma\).

The barrier state is valid precisely when
\[
                         h_i(u)<b_i(u)<h_{i+1}(u),       \tag{BS8}
\]
omitting a missing endpoint. At a generic exchange this makes all
contracted original columns cheaper, and the additional cheaper
own-arm columns are rank-irrelevant by (BS6).
Consequently the projected nonloop-parallel test equals the full
original cheaper-column test. This exactness is stronger than the
necessary projections used for annihilator exceptions.

Each barrier cost has slope at least \(t\) for \(\epsilon=0\),
and at least \(2t\) for \(\epsilon=1\). Both exceed \(\gamma\).
Thus \(b_i-h_i\) and \(b_i-h_{i+1}\) are strictly decreasing.
The valid state ratios form one interval. Birth of a finite threshold
can occur only once as \(C\) grows, and once \(e_i\) is a loop it
remains one.

## Count the active window

Put \(\Lambda_i(u)=Yu+Z-b_i(u)\). The original price at an exchange is
\[
 \lambda=\Lambda_i(u)\quad(\epsilon=0),\qquad
 \lambda=-\Lambda_i(u)\quad(\epsilon=1).                 \tag{BS9}
\]
For lower \(E\), the two price ranges have the same window bounds as
in the weighted-partner proof:
\[
 \Delta\le\frac{t}{Y-t+1}\quad(0<\lambda<t),\qquad
 \Delta\le\frac{ta}{(Y-2t+1)^2}\quad(t<\lambda<tu).
                                                               \tag{BS10}
\]
For upper \(E\), \(-tu<\Lambda_i<0\). Equation (BS7) makes
\(\Lambda_i\) increase at rate at least \(Y-2t+1\). Since
\(b_i\le(2t-1)u+a-1\), positivity of \(\lambda\) gives
\(u<a/(Y-2t+1)\). Comparing two active ratios again yields the
second bound in (BS10). Intersecting with (BS8) can only shorten
these windows.

In a ratio interval of width \(\Delta\), the lower/same-layer
entry flag changes at most \(t(Y\Delta+1)\) times. Opposite Z
orders change at at most
\[
                         \tfrac12t(t-1)\Delta+t-1
\]
distinct parameters. Between such events the death column is fixed.
The state and price conditions cut at most one subinterval of each
cell. A genuine cross-layer wall cannot be confined to an isolated
order-event ratio.

These counts give at most \(4t^2\) partners per state in the first
range and \(16at^2/Y\) in the other two. There are at most
\(t^2+2\le2t^2\) states and at most \(t^2\) own-order sectors.
This proves (BS3). The count also bounds the unsubdivided eligibility
intervals produced by this model. It does not assert that one
literal pair's entire eligibility is connected.

## Directions and the complete restricted fan

Summing over \(Y\) and its \(t\) possible short coordinates gives
(BS4). The range \(Y<4t\), same-arm exchanges and exchanges involving
a bridge have \(O(at)\) directions by their difference boxes and the
[bridge inventory](hook_bridge_reduction.md).
The [parallel-shadow theorem](cross_layer_direction_budget.md) gives
\(O(n\log(n+1))\) actual pieces per direction, including recurrent
pieces and intersections with other fan walls. Since \(n=\Theta(at)\),
this proves the stated cross-layer edge bound.

It remains to include same-layer walls. In the \((\lambda,u)\) plane
these are vertical. Distinct such lines do not meet in \(\mathcal R\).
After deleting all nonvertical edge pieces, a vertical component that
separates the convex wedge must reach its boundary. Restoring each
nonvertical piece can add at most one region. Hence the number of
cells is at most the number of nonvertical pieces plus the number
of vertical boundary intersections, up to an absolute constant.

Only the boundaries \(\lambda=0\) and \(\lambda=tu\) contribute
transverse vertical intersections. The boundary \(u=1\) is itself
vertical; the infinite-ratio boundary in the normalized weight
triangle meets no finite vertical line. Along the first two boundaries,
the effective planar labels are respectively
\[
                            (y,z),\qquad(y+tx,z).
\]
Every nonzero residual label has product bound \(O(at)=O(n)\)
in these two coordinates. For upper Y labels the shift by \(t\)
is at most a constant factor; for upper Z labels it makes the
short coordinate at most \(2t\). Lower labels are unchanged.
The [planar profile bound](planar_dyadic_bound.md) therefore gives
\(O(n\log(n+1))\) boundary changes.

A vertical wall has nonzero planar profile difference on either
transverse boundary, so its intersection is one of those changes.
There is at most one vertical line at each boundary ratio.
Degenerate subdivisions on the vertical boundary are charged to
their incident nonvertical pieces. This proves the cell bound in
\(\mathcal R\). Transposing Y and Z covers
\(u<1,\lambda<t\), proving (BS5).
An initial cone has a connected intersection with either convex wedge,
so counting cells bounds the number of original cones meeting them.

## Verification and limits

~~~bash
python3 -B scripts/run_hook_barrier_sector_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker verifies the saved bridge and two-colour report manifests
before using their comparison, and records that (BS5) is dominated by
(BS0). It computes the intersection flag by original field ranks,
tests all retained versus omitted prefixes on relevant columns,
constructs the state quotients, and compares candidate death columns
with full original cheaper-column ranks and strict adjacent minima.
Both choices of the heavy-axis label's layer are exercised. The
checker also enumerates a complete original fan, verifies exact
restricted-area coverage, and independently matches vertical boundary
intersections with lower envelopes of all basis profiles. The
source-prefix and bounded-slope ingredients are classical linear and
parametric-matroid operations; compare
[Eppstein, Section 1.2](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf).

No polynomial lower separation or new unrestricted exponent follows.
The exact weighted barrier treatment of the complementary high-price
sector and arbitrary outer images remains open. Those cases already
have the coarse full count (BS0). Exact priority and the original
final-paper requirements remain open.
