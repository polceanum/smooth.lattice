# Rank departures from one transverse arm control the whole hook fan

A common quotient makes it possible to count many interacting partner
lists together. For every compatible flat hook, a low-rank departure
from a model with one purely transverse arm gives a small exact decision
matrix on each of only \(O(n^2)\) coefficient-independent chambers.
Both outer images are allowed, and no unit or prefix-regularity
hypothesis is needed.

The count is
\[
                         g(I_\phi)=O(n^2(1+\rho)^3),       \tag{RD1}
\]
where \(\rho\) is the least departure rank defined below. This is a
full-fan parameter theorem. It does not improve the unrestricted
worst-case exponent or prove an actual polynomial lower separation.
Exact current novelty and the original final-paper requirements remain open.

## A common-quotient compression lemma

Let \(M,M_0\) be rank-\(R\) matroids on the same labelled ground,
with a common quotient \(Q\) of rank \(R-r\). Partition the ground
into \(k\) groups. On a convex weight region suppose:

1. one fixed basis \(B_0\) is the unique minimum of \(M_0\);
2. the internal cost order of every group is fixed.

Then a fixed independent set \(C\) and a minor of \(M\) recover
every generic optimum throughout the region. The minor has rank at
most \(kr\) and at most \(k^2r\) original labels. Every recovered
basis is \(C\) union the minimum basis of that minor. In particular,
the whole cone patches and support functions agree up to the fixed
linear cost of \(C\).

To prove this, greedy acceptance in a quotient implies acceptance
in the parent: a rank increase in \(Q\) forces a rank increase in
\(M\), or in \(M_0\), on the same prefix. Thus its greedy basis
\(B_Q\) lies in both parent greedy bases, in particular in \(B_0\).
The restriction \(Q|B_0\) has rank \(R-r\), so its dual has rank
\(r\). In each of the \(k\) lists \(B_0\) intersected with a group,
greedily retain a basis of this dual in **descending** cost order.
Let \(D\) be their union. Then \(|D|\le kr\). A dual column
discarded in its list is spanned by more expensive retained columns,
and cannot enter the merged maximum greedy basis. Consequently
\[
              B_0\setminus B_Q\subseteq D,\qquad
              C=B_0\setminus D\subseteq B_Q.             \tag{RD2}
\]
This \(C\) is fixed on the entire region and is compulsory in every
generic optimum of \(M\). Its contraction has rank
\(R-|C|=|D|\le kr\). Increasing greedy compression in each of the
original \(k\) lists now retains at most \(kr\) labels per list.
Their union has at most \(k^2r\) labels and preserves the merged
minimum. This proves the lemma without separately enumerating the
possible \(B_Q\)'s. At \(r=0\) the decision is empty.

For represented matroids of equal full row rank, matrices \(V,V_0\)
with \(\operatorname{rank}(V-V_0)=r\) supply such a common quotient:
project the output by the image of their difference. The common
projected matrix has rank \(R-r\). More generally the lemma applies
to abstract common quotients; representation is not needed for its
greedy argument. The quotient and rank-perturbation background is
classical; see [Geelen--Gerards--Whittle, Section 10 and Proposition 10.1](https://www.math.uwaterloo.ca/~jfgeelen/Publications/structure.pdf#page=22).

For \(k=3\), if the original costs have the binary-coordinate form
\(c_e(\lambda,u)=x_e\lambda+y_eu+z_e\), \(x_e\in\{0,1\}\),
the saved [two-colour theorem](two_colour_threshold_braids.md) counts
one resulting minor by \(O((1+r)^3)\): its rank is at most \(3r\)
and its ground size at most \(9r\). Restriction to a convex region
cannot split a strict optimality cone.

## The departure rank for a compatible hook

Use the [original residual model](hook_bridge_reduction.md) for
\[
 M=(y^a,y^tz^t,z^a),\quad a\ge2t,\quad t\ge2,\quad
 n=4at-2t^2,\qquad A=a-t.
\]
Its output is \(T_Y\oplus T_Z\), of rank \(R=2At\), and its
fixed three ground groups are \(E_Y,E_Z,H\). Upper identity columns
remain at every compatible specialization.

For any \(C(v)\in K[v]/(v^t)\), form a comparison map as follows:
retain every Z component of all three generator images; erase their
Y components, and replace only the central Y image by \(y^tC(z)\).
Call its residual matrix \(V_{Y,C}\), and the actual matrix \(V\).
All these choices remain compatible: the allowed generator-image
coefficients in the original model are independent. Put
\[
 r_{Y,C}=\operatorname{rank}(V-V_{Y,C}).                  \tag{RD3}
\]
Define \(V_{Z,C},r_{Z,C}\) by transposition, and let
\[
            \rho=\min\{r_{Y,C},r_{Z,C}:C\in K[v]/(v^t)\}.\tag{RD4}
\]
The minimum is attained, since attainable ranks form a nonempty finite
set of integers. Any one chosen comparison polynomial gives a valid
upper bound; no procedure for optimizing over all \(C\)'s is required.
This is actual matrix rank, not support size or full bridge-span coupling.

For a central-only hook with arm operators \(W_Y,W_Z\),
\[
 r_{Y,C}=\operatorname{rank}\bigl(m_{W_Y-C(v)}:R_A\to R_A\bigr),
 \qquad R_A=K[q,v]/(q^A,v^t),                             \tag{RD5}
\]
and similarly on Z. Indeed the central bridges supply all inputs
\(q<t\), and the pure central sources supply all inputs \(q\ge t\).
Taking constant \(C\) recovers the ranks of the off-axis operators
when the axis units are constant, but (RD3) also includes outer images
and nonconstant axis polynomials.

## A common quadratic partition for all comparison models

Fix the Y comparison arm. For every local long index \(q\ge t\),
the Y row fiber is an isolated rank-\(t\) component. Its columns are
the \(t\) upper identities and the \(t\) sources obtained by
multiplication by \(C(v)\). No other comparison column reaches this
fiber. Their common long-coordinate cost cancels in comparisons;
all possible positive cross-layer ties are among
\[
                       w_x=jw_z,\quad1\le j<2t.         \tag{RD6}
\]
Within-layer orders are fixed. Thus these \(O(t)\) planes fix the
minimum in every such fiber, including singular or zero \(C\).
No assertion that its source columns are parallel is needed.

After removing the fibers and loops, every remaining comparison label
lies in the union of two boxes:
\[
 \begin{split}
 &\{0,1\}\times[0,2t)\times[0,2a),\\
 &\{0\}\times[a,a+t)\times[0,t).                         \tag{RD7}
 \end{split}
\]
The second is the possible \(y^a\) outer bridge. Its shifts beyond
short Y index \(t-1\) become loops when the Y image is erased.
The other outer bridge already lies in the first box. Differences
within and between these boxes have only \(O(at)\) distinct values:
one transverse difference varies over \(O(t)\) values, even for the
translated second box, and the other over \(O(a)\).

Include their comparison planes, (RD6), and all original three-group
internal-order planes, also \(O(at)\). The resulting arrangement
has \(O(a^2t^2)=O(n^2)\) convex chambers in the positive weight
triangle. It is independent of all coefficients and of \(C\). On each
chamber the comparison minimum \(B_0\) and the three internal orders
are fixed. Transposition gives the Z version.

Apply the lemma to \(V,V_{Y,C}\), or to the Z comparison. Each
chamber has an exact residual decision of rank at most \(3r_{Y,C}\)
on at most \(9r_{Y,C}\) original binary-coordinate labels. The
two-colour count proves (RD1). Combining it with the saved
[actual bridge-rank bound](hook_prefix_bridge_rank.md), with
\(b=\operatorname{rank}H\), gives
\[
                g(I_\phi)=O\bigl(n^2(1+\min(b,\rho))^3\bigr).
                                                               \tag{RD8}
\]
The minor and compulsory set depend on the actual coefficients and
the chamber. A common small weight list of this cardinality is not
inferred, and the compulsory span is not declared to be an algebra ideal.

## Consequences and the unresolved range

Fixed departure rank gives a quadratic full-fan bound independently
of thickness, arbitrary coefficients on the other arm and compatible
outer images. In particular any actual family with
\(g(I_n)\ge c n^{2+\epsilon}\) must satisfy
\[
                             \rho=\Omega(n^{\epsilon/3}).\tag{RD9}
\]
Thus every transverse comparison on either arm must have growing
departure rank. The exact small minor is an additional reconstruction
consequence, rather than just a count of possible bases.

For example, in the central-only constant-axis class let
\(W_Y=1+q^{A-1}v\), and allow the Z arm to be arbitrary. Then
\(r_{Y,1}=t-1\), so the new full bound is \(O(n^2t^3)\), compared
with the coarse \(O(n^2t^6)\) full-hook ceiling. At \(a=t^{23}\)
this is \(O(n^{17/8})\), a polynomial improvement over that
\(O(n^{9/4})\) ceiling. Smaller coefficient-specific regularity
bounds, the quadratic-coordinate bound and the box bound still apply.
This comparison does not assert sharpness or an unrestricted improvement.

The [linear high-price partner family](hook_high_price_partners.md)
has \(r_{Y,1}=2\). Its many partners are compatible with a quadratic
total: the reduction groups them through a common quotient instead
of bounding each partner list. When departure rank and bridge rank
are both large, (RD8) leaves the original collective problem open.

## Verification and primary-source boundary

~~~bash
python3 -B scripts/run_hook_rank_departure_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker reconstructs original quotient columns and commuting
multiplication for both actual and comparison maps, verifies the
common quotient and dual identities, and compares whole original
cone patches with the retained minor. It includes nonzero outer
images, zero and singular transverse bases, both arm orientations,
small characteristics, rational cases, and independent abstract
matrix updates. Complete original fans provide further controls.

Common quotients, represented rank perturbations and greedy compression
are classical. The cited Geelen--Gerards--Whittle Proposition 10.1
identifies perturbation distance through row-space intersections; its
structure-theory claims are not used. The count uses the saved two-colour
theorem and the explicit hook label inventory. This limited comparison
does not establish strongest-current exact novelty of the combined
parameter theorem. The unrestricted endpoint, a polynomial lower
separation, sharpness and final-paper criteria remain unresolved.
