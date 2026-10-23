# Inverse propagation alone bounds both directed lower-source classes

Long forward propagation alone cannot produce a polynomial excess in
a compatible hook with two unit central arms. Every pure cross-layer
exchange has one lower source. Its eligibility can be bounded using
only the inverse on that source's arm. Combining the two orientations
and the established parallel-wall budget gives a quadratic polylogarithmic
full-fan bound when both inverse long degrees are of order the thickness.
All compatible outer images are allowed, over every field.

This closes the two geometric central-arm candidate even when its
forward expansions are long and contain several independent steps.
It does not close the plus/plus or power/geometric models when an
inverse has large degree. Unrestricted sharp growth, actual polynomial
lower excess, exact strongest-current novelty and the original final-paper
requirements remain open.

## The full bound and its one-arm form

Use the [original compatible hook](hook_bridge_reduction.md), with
\[
 t\ge2,\quad a\ge2t,\quad A=a-t,\quad n=4at-2t^2.
\]
Suppose both central operators are units in their rectangular rings
\(K[q,v]/(q^A,v^t)\). Put
\[
 d_i=\deg_q(W_i^{-1}),\quad D=\max(d_Y,d_Z),\quad
 H=1+\log(a/t),\quad \ell=1+\log(n+1).
\]
Every compatible coefficient of both outer images is arbitrary. Then
\[
                       g(I)=O\bigl(n^2(H+D/t)\ell\bigr).   \tag{IO1}
\]
If both outer images vanish, the stronger bound is
\[
                       g(I)=O\bigl(n^2(1+D/t)\ell\bigr).   \tag{IO2}
\]
No bound on \(\deg_q W_Y\) or \(\deg_q W_Z\) is assumed.
Earlier bounds remain available whenever smaller; neither logarithmic
factor is asserted sharp.

There is also an oriented statement requiring only the chosen Y central
operator to be a unit. Let \(\mathcal D_Y\) be the actual directions
of pure cross-arm exchanges between a lower Y source and an upper Z
column. With \(d=\deg_qW_Y^{-1}\), a coefficient-independent inventory
contains \(\mathcal D_Y\) and has
\[
                O\bigl(a(d+tH)\bigr)\text{ planes},       \tag{IO3}
\]
of original integer height \(O(a)\). If the Y components of both
outer images are zero, its size improves to \(O(a(t+d))\).
The opposite central operator and all its outer components may be
arbitrary, including nonunit. This controls only the specified orientation.
The whole-fan result uses both unit arms.

The inventory in (IO3) covers the relevant nonzero-first-coordinate
directions. It is not a complete comparison arrangement for the entire
fan. No common complete weight list of the size in (IO1) is claimed.

## Zero chosen-arm outer images: use the lower identity chart

Choose Y and write the objective as \((\lambda,u,1)\), with
\(\lambda,u>0\). An original lower central source has label
\[
                       E=(0,M,s),\quad M\ge2t,\quad
                       t\le s<2t,\quad M<a.
\]
Multiply only the Y rows by the unit \(W_Y^{-1}\). Labels and
costs stay unchanged. Pure lower sources become identities on the
Y tail of local inputs \(q\ge t\). Central bridges have head
inputs \(q<t\). Upper inputs with \(q\ge t\) remain in the
tail; those with \(q<t\) have every output of local cost at most
\((t-1+d)u+t-1\).

If \((M-t)u+s-t\) exceeds this bound, contract the span of strictly
cheaper tail identities. Every head-input upper column now has only
head Y support. With zero Y components of the outer images, the whole
remaining matrix splits into its surviving Y tail and its head-Y-plus-Z
part. The proposed exchange crosses those two summands and is impossible.
The chosen identities need not belong to a particular minimum basis:
further contraction of the entire strictly cheaper span preserves the
direct sum, and every minimum basis uses its fixed ranks.
Thus every coherent exchange of this lower source across the arms satisfies
\[
                    (M-t)u+s-t\le(t-1+d)u+t-1.            \tag{IO4}
\]
This is the lower half of the
[unit-window argument](hook_coherent_unit_windows.md); it uses only
the inverse degree, not the earlier maximum of two degrees.

Write the opposite upper label as \(F=(1,r,K)\), \(r<t\).
Equality of their original costs gives
\[
                   K=u(M-r)+s-\lambda\le uM+s.           \tag{IO5}
\]
If \(u\ge1\), (IO4) bounds \(M\le3t+d\).
If \(u\le1\), it bounds \(uM+s\le4t+d\), hence
\(K\le4t+d\). Therefore every such exchange has
\[
                         M\le3t+d\quad\text{or}\quad K\le4t+d.
                                                               \tag{IO6}
\]
One long difference coordinate is in an \(O(t+d)\) strip, while
the other has \(O(a)\) possible values. This proves the sharper
zero-outer form of (IO3). No upper-column price bound or forward
operator degree was used.

## Arbitrary outer images: only the full cheaper prefix must split

Still use the inverse row chart. Put \(B=4t+2d\). For a central
lower source every coherent pure cross-arm exchange satisfies
\[
             uM\le Bu+4t\quad\text{or}\quad
                                  u(a-M)\le Bu+4t.        \tag{IO7}
\]
Here is the prefix argument, retaining every outer coefficient.
Assume both inequalities fail. The common exchanged cost is
\(\tau=uM+s<au\), so no \(y^a\)-source is in the prefix.
Contract all tail lower identities of cost strictly less than \(\tau\).
All tail outputs of head-input upper columns are killed, by the first
failed inequality. Central bridges have only head Y inputs in this chart.

An upper column in the \(\tau\)-prefix has input long label at most
\(M+2t/u\); its inverse output long label is at most
\(M+2t/u+d<a-t\), by the other failed inequality.
The Y projections of \(z^a\)-bridges lie in original long labels
at least \(a-t\), and polynomial multiplication by the inverse
preserves that long-coordinate ideal. Thus the entire projected
\(\tau\)-prefix splits into middle Y rows and head-Y-plus-far-Y-plus-Z
rows. The exchanged Y identity is in the first part; an opposite pure
column is in the second. The full strictly cheaper span is a sum of subspaces in those two
summands. Contracting it preserves the split, so fixed prefix ranks
forbid their exchange, including objectives with additional ties.
Higher-cost columns may still mix these parts and are not discarded.

This is the lower-source proof of the
[outer-window theorem](hook_outer_unit_windows.md). Every occurrence
of propagation in that proof belongs to the inverse chart; the forward
degree is unnecessary for (IO7).

A pure outer Y source has label \((0,a+i,j)\), \(i\ge t\).
It belongs to the original pure central span because \(W_Y\) is a
unit and multiplication by its inverse preserves the tail ideal.
All these central columns cost at most \((a-1)u+2t-1\).
The outer column can therefore be eligible only if
\[
                         (i+1)u+j\le2t-1.                \tag{IO8}
\]
This argument also uses no forward-degree bound.

## The oriented inventory has a translated hyperbola

For any chosen lower Y label \((0,M,s)\) and opposite upper
\((1,r,K)\), put
\[
                   b=M-r,\qquad c=K-s,\qquad S=t+d.
\]
Their tie plane is \(\lambda=bu-c\), so \(c<ub\).
Every actual such direction satisfies
\[
 b\le16S\quad\text{or}\quad c\le16S\quad\text{or}\quad
 |b-a|\le16S\quad\text{or}\quad |b-a|c\le16at.           \tag{IO9}
\]
For a central column the first alternative in (IO7) gives the first
two strips according to \(u\ge1\) or \(u\le1\).
The other alternative gives the translated strip when \(u\ge1\).
When \(u<1\), put \(\Delta=a-M\). Outside that strip,
\(\Delta-B\ge\Delta/2\), \(|b-a|\le2\Delta\), and
\((\Delta-B)u\le4t\). Since \(c<ub<ua\), their product is
at most \(16at\). For an outer source, (IO8) gives \(iu\le2t\)
and \(c<au+2t\), hence \(ic\le4at\); outside the translated
strip \(|b-a|\le i\). Negative \(c\) is already in a short strip.

The three strips contain \(O(aS)\) integer pairs. For the positive
translated product, sum \(O(at/h)\) possible second coordinates over
\(h=|b-a|\), splitting at a constant multiple of \(t\).
The result is \(O(at(1+\log(a/t)))\). This proves (IO3).
Transposition gives the other lower-source orientation with its own
inverse degree.

## Count actual parallel pieces, not just supporting planes

Every nonzero-first-coordinate wall outside the original
\(O(at)\) known-direction set admits a pure cross-arm exchange.
Exactly one of its labels is lower. If both central arms are units,
it is covered by one of the two inventories above. Together with
known directions, there are
\[
                         L=O(a(D+tH))                    \tag{IO10}
\]
possible nonzero-first-coordinate directions. With zero outer images,
\(L=O(a(t+D))\).

The residual represented matroid has rank \(2At\le n\), labels
with first coordinate zero or one and divisor products at most
\(8at=O(n)\). The established
[translated-layer shadow bound](layer_shadow_budget.md) therefore
gives \(O(n\ell)\) actual parallel wall pieces for each primitive
direction with first coordinate one, uniformly in the other entries.
This counts subdivisions at other fan vertices. It applies directly
to the represented matrix; its rank is not misidentified as an algebra
length. The
[finite polytope sweep](cross_layer_direction_budget.md) counts the
whole positive fan by its boundary shadow and these nonzero-first-coordinate
pieces. Hence \(g=O(n\ell(1+L))\). Substituting (IO10) and
\(at=\Theta(n)\) proves (IO1)--(IO2).

For any fixed \(\epsilon>0\), a family in this two-unit class with
\(g\ge c n^{2+\epsilon}\) would consequently require
\[
                      D/t=\Omega(n^\epsilon/\ell).        \tag{IO11}
\]
The logarithmic term in \(H\ell\) is smaller than any positive
power of \(n\). Thus a long forward expansion without a correspondingly
large inverse degree is excluded as a cause of polynomial excess.
Large inverse degree is only a necessary condition, not a lower construction.

## Actual remainder directions can continue through a common vertex

A finite exact example prevents inferring a universal noncrossing rule
from the two-step support restriction. Work over \(\mathbb F_{1009}\),
with \(t=3,a=16\), zero outer images and
\[
\begin{split}
 W_Y&=1+(771q^3+191q^4)v+(985q^6+544q^7+815q^8)v^2,\\
 W_Z&=1+(730q^3+221q^4)v+(430q^6+986q^7+570q^8)v^2.
\end{split}
\]
Both forward and inverse supports lie in the narrow band
\(3j\le i\le4j\). The two primitive remainder directions
\[
             (1,-6,3),\qquad(1,3,-6)
\]
are outside the original known-direction inventory. Both have actual
wall segments continuing on opposite sides of
\(w_0=(3,1,1)/5\). For the first direction the other endpoints are
\((12,5,6)/23\) and \((24,7,6)/37\); transpose the last two
coordinates for the second direction.

At the four segment midpoints the coherent exchanged pairs are,
respectively,
\[
\begin{array}{ll}
 (0,8,4)\leftrightarrow(1,2,7),&
 (0,7,5)\leftrightarrow(1,1,8),\\
 (0,4,8)\leftrightarrow(1,7,2),&
 (0,5,7)\leftrightarrow(1,8,1).
\end{array}
\]
For each pair the entire strictly cheaper span has rank 52; adjoining
either column or both raises that rank by exactly one. Exact adjacent
weights give two distinct minimum bases differing by that single pair,
and their entire normal polygons contain the displayed segment.
This proves continuation of both remainder directions through one
interior vertex. Other walls also meet there; it is not asserted to
be an isolated four-cone crossing. No growing number of intersections,
generic-coefficient assertion or asymptotic lower bound follows from
this finite certificate.

## Long-forward families and verification

For \(M=q^\nu(1+q)v\), take both central arms
\[
                          W=(1-M)^{-1}.
\]
The inverse has long degree at most \(\nu+1\), although \(W\)
can have degree of order \(a\). For \(\nu=O(t)\), (IO1) gives
\(O(n^2\log^2(n+1))\) with every compatible outer image arbitrary,
and (IO2) gives \(O(n^2\log(n+1))\) when they vanish.
Another example is \(W=(1-q+q^\nu v)^{-1}\): its forward long
degree is exactly \(A-1\), while its inverse degree is at most
\(\max(1,\nu)\). Both statements hold in every field.

~~~bash
python3 -B scripts/run_hook_inverse_orientation_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker reconstructs original quotients and both inverse charts,
checks entire middle-prefix splits and independent contraction ranks,
and tests oriented exchanges in complete original fans. It includes
arbitrary nonunit opposite arms in the one-arm statement and arbitrary
outer maps. Long forward degree is checked independently of the small
inverse. No bound on high-cost columns is inferred from a cheaper-prefix split.
The four remainder segments are separately reconstructed in the original
quotient, with whole normal polygons and entire cheaper-span ranks.

The lower identity-chart argument, integer inventories, matroid faces
and the parallel-shadow counting input are established techniques.
The face dictionary is reviewed in
[Feichtner--Sturmfels, Section 2](https://arxiv.org/pdf/math/0411260v1).
Exact strongest-current novelty of this combined inverse-only class
bound is unestablished. The remaining two-step plus and asymmetric
power/geometric models may have large inverses. Their finite remainder
faces do not yet supply an asymptotic lower count. The unrestricted
research endpoint remains open.
