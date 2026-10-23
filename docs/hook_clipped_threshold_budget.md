# A clipped threshold amortizes all bridge states

In the constant-axis central-hook family, the bridge states can be
counted together. The crucial observation is that one original label
has just one full opposite-layer death threshold at each generic ratio.
The states are an exact representation of parts of this one curve.
They are not independent curves whose maximum window lengths must be
summed separately.

For
\[
 I=(x^2,y^a,z^a,y^tz^t-xf),\quad a>2t,\quad t\ge2,\quad
 n=4at-2t^2,
\]
assume both outer images vanish and the arm operators are
\(\alpha+N_Y,\beta+N_Z\), with nonzero constants and
\(N_i\in(v)K[q,v]/(q^{a-t},v^t)\) arbitrary. All statements are
over every field. Normalize \(w=(\lambda,u,1)\), \(u>1\), and use
the wedge \(0<\lambda<tu\).

Fix \(E=(\epsilon,Y,Z)\), where \(Y\ge(3+\epsilon)t\), and
\(t\le Z<2t\) for \(\epsilon=0\), \(0\le Z<t\) for
\(\epsilon=1\). Let \(P(E)\) count its distinct genuine pure Z
partners in the other layer, throughout the wedge. Then
\[
                  P(E)\le30at^4/Y+4t^5.                 \tag{CT1}
\]
In particular \(Y\le a/t\) gives \(P(E)\le34at^4/Y\), without
an annihilator hypothesis. The constants are not asserted sharp.

The entire restricted fan consequently satisfies
\[
 g_{\rm low}(I)
   =O\!\left(n^2t^4\bigl(t+\log(n+1)\bigr)\log(n+1)\right),
 \qquad w_x<t\max(w_y,w_z).                              \tag{CT2}
\]
This improves the previous \(n^2t^6\log^2(n+1)\) sector bound.
It can also beat the saved coarse full-fan ceiling \(O(n^2t^6)\)
when thickness grows faster than the logarithm. It is still a
restricted-sector theorem: the complementary high-price sector,
arbitrary outer images, the unrestricted endpoint and exact novelty
remain open. Retain all smaller existing full-fan ceilings.

## One full threshold in each own-arm order sector

Fix an own-arm order sector \(J\subset(1,\infty)\). There are at
most \(t^2\) such sectors by the
[earlier order inventory](hook_weighted_partners.md).
Contract the earlier same-layer pure Y columns and all central
bridges. The first set is fixed on \(J\); the bridges are genuinely
cheaper than \(E\) throughout the wedge at the displayed depth.
This is an exact vector quotient, not an asserted algebra quotient.

The remaining same-layer Z columns form an increasing flag \(C(u)\).
A column \((\epsilon,r,s)\) enters it at
\(u=(s-Z)/(Y-r)\). There are \(t\) possible \(r\)'s, all below
\(Y\). Let \(D\) now contain **both** opposite-layer arms, with
costs \(Ru+S\). Define \(b(u)\) as the least cost prefix of this
whole list that spans the image of \(E\) together with \(C(u)\).
Use \(+\infty\) when unreachable and \(-\infty\) when already a
loop. The unique death column at a generic ratio is the original
opposite-layer partner.

Between entries of \(C\), \(b\) is continuous: it is a minimum of
finitely many maxima of affine costs over spanning subsets. Entries
of \(C\) can only make it jump down. Reachability is born at most
once, and an already-loop image remains a loop. No bound on all
opposite Y slopes is being assumed.

The original exchange price on this curve is
\[
 \lambda_E(u)=Yu+Z-b(u)\quad(\epsilon=0),\qquad
 \lambda_E(u)=b(u)-Yu-Z\quad(\epsilon=1).                 \tag{CT3}
\]
Thus all bridge states at the same ratio must give the same death
column whenever they give an actual cross-arm exchange. Their active
ratio intervals have disjoint interiors. This is equality of the full
original rank test, as proved in the
[critical-bridge representation](hook_barrier_sectors.md).

## Clip the price, then bound its variation

Put
\[
 d=Y-2t+1,\qquad U=a/d.
\]
A Z death column active in the wedge necessarily has \(1<u<U\).
For a lower \(E\), use \(R\le t-1,S<a\) and
\(Yu+Z-(Ru+S)<tu\). For an upper \(E\), use
\(R\le2t-1,S<a\) and \((Ru+S)-Yu-Z>0\).

Define the clipped function
\[
 F(u)=\begin{cases}
 \min(tu,\max(0,\lambda_E(u))),&\epsilon=0,\\
 -\min(tu,\max(0,\lambda_E(u))),&\epsilon=1.
 \end{cases}                                            \tag{CT4}
\]
The extended threshold values have the natural clipped limits.
Every jump of \(F\) is nonnegative, including birth and loop events.
On a Z-death portion inside the wedge, \(F'\ge d\).

On a Y-death portion inside the wedge, \(F'\ge-2t\).
For a lower label, positivity of
\((Y-R)u+Z-S\) and \(S\ge0,Z<2t\) gives \(R-Y<2t\).
For an upper label, the inequality
\((R-Y)u+S-Z<tu\), with \(S\ge t>Z\), gives \(R-Y<t\).
Outside the wedge the clipped slopes are among \(0,t,-t\), so
the same lower bound holds everywhere else. These are bounds on
the clipped curve; the unrestricted Y-column slopes may be large.

Let \(M\) be the total length of its Z-death portions in the wedge
on \(J\cap(1,U)\). Integrating the derivative and including its
nonnegative jumps gives
\[
 dM\le tU+2tU=3tU,
 \qquad M\le\frac{3ta}{d^2}\le\frac{27ta}{Y^2}.          \tag{CT5}
\]
The endpoint difference of \(F\) is at most \(tU\) for either
sign in (CT4). The last inequality uses \(Y\ge3t\).
If \(U\le1\), there is no active portion. This argument counts the
total measure of all states together, without multiplying by their
number.

## Count partner changes on the active union

The critical-bridge representation gives at most \(t^2+2\) states
on \(J\). Each state's barrier and price conditions cut one interval:
its Z threshold has upper slope at most \(t-1\) or \(2t-1\),
whereas each critical Y cost has greater slope. Downward jumps only
preserve this interval property. Fill the finitely many order-event
points when counting interval components. Therefore the active Z
union has at most
\[
                              L\le t^2+2               \tag{CT6}
\]
components, with total length \(M\).

For each same-layer Z short coordinate, its entry parameters are
contained in an arithmetic progression of spacing \(1/(Y-r)\). Across the active
components there are at most \(YM+2L\) relevant entries for that
coordinate, including endpoint conventions. Thus these events
contribute at most \(tYM+2tL\).

Count Z--Z opposite-layer order changes over the **entire** interval
\((1,U)\), rather than separately inside every state. Their distinct
parameters have denominator at most \(t-1\), giving at most
\[
              \tfrac12t(t-1)U+t-1\le t^2U.              \tag{CT7}
\]
Within one active component, a Z partner cannot change elsewhere:
with \(C\) fixed the full threshold is continuous, and two different
Z costs cannot meet without a Z--Z order event. A critical-state
boundary alone cannot change that partner. Each component contributes
one initial partner. Recurrent partners are allowed throughout.

Consequently
\[
 P_J(E)\le tYM+t^2U+(2t+1)L
          \le30at^2/Y+4t^3.                              \tag{CT8}
\]
Here \(U\le3a/Y\) and
\((2t+1)(t^2+2)\le4t^3\) for \(t\ge2\).
Summing the at most \(t^2\) sectors proves (CT1).

## Directions, the restricted fan, and localization of an excess

There are \(t\) short coordinates for each heavy long coordinate.
Summing (CT1) over both layers yields
\[
             D_{\rm low}=O\bigl(at^5(t+\log(a+1))\bigr). \tag{CT9}
\]
Short heavy coordinates, same-arm exchanges and bridge exchanges add
only \(O(at)\) directions by the saved difference inventory.
The established \(O(n\log(n+1))\) parallel-shadow bound retains all
recurrent pieces and subdivisions, giving (CT2) for nonvertical edges.
The same-layer vertical walls are included by the boundary-shadow
sweep in the critical-bridge proof. Transposition covers the other
heavy-axis orientation. No connectedness of a literal pair is used.

For example, take \(a=t^{23}\), so \(t=\Theta(n^{1/24})\).
The earlier full ceilings give \(O(n^{9/4})\) as their minimum,
whereas (CT2) gives \(O(n^{53/24}\log(n+1))\) for this restricted
sector. This is a polynomial saving in the sector alone.

More generally, suppose \(g(I_n)\ge c n^{2+\epsilon}\) and
\(t\le n^{\epsilon/5-\delta}\), for fixed positive
\(\epsilon,\delta\). Then (CT2) is \(o(n^{2+\epsilon})\).
Thus all but a vanishing fraction of the initial cones must avoid
the low wedge. A polynomial excess in this thickness range would
have to be carried by the complementary high-price sector.
This is a localization consequence, not an exclusion of that sector.

## Verification and limits

~~~bash
python3 -B scripts/run_hook_clipped_threshold_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker constructs the full opposite-arm threshold independently
of the bridge-state thresholds, over every exact ratio-order cell in
its finite range. It integrates the clipped curve and all jumps using
rational arithmetic, compares the union measure and partners with
the state representation, checks the event charging and original
cheaper-column ranks, and replays a complete original restricted fan.
An explicitly abstract two-dimensional represented-matroid control
exercises two separate Z phases and an upward clipped jump; it is
not presented as a hook quotient.

The greedy and parametric-matroid setting is classical; compare
[Eppstein, Section 1.2](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf).
The variation argument is supplied explicitly here. These checks and
that comparison do not establish its exact current novelty. The
uniform full-fan advance, polynomial lower separation, strongest-current
priority assessment and original final-paper requirements remain open.
