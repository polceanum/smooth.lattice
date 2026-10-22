# Three small bridge blocks isolate the hook counting problem

Let \(B=K[y,z]\), \(a\ge2t\), \(t\ge2\), and
\[
 M=(y^a,y^tz^t,z^a),\qquad C=B/M,\qquad
 m=2at-t^2,\qquad n=2m.
\]
For **every** compatible map \(\phi:M\to C\), take its flat square-zero
extension \(A_\phi\). No genericity or nonvanishing hypothesis is
imposed. The exact residue matrix has a partition into two pure arms
and exactly \(3t^2\) possible bridge columns. On every region of fixed
internal orders this reduces the remaining entire fan patch to a dual
matroid of rank \(3t^2\) on at most \(9t^2\) original labels.

There is also a global counting localization. Let \(\mathcal N\) be
the primitive cross-layer directions obtainable from two labels in one
pure arm or from a pair involving a bridge label. Then
\[
 |\mathcal N|=O(at)=O(n),\qquad
 g(I_\phi)\le O(n^2\log(n+1))+E_{\rm rem}(I_\phi),              \tag{H1}
\]
where \(E_{\rm rem}\) counts actual positive edge pieces with
nonzero first component and direction outside \(\mathcal N\).
Every such direction must be represented by a pure cross-arm exchange.
Parallel pieces and subdivisions at fan vertices are retained in this
count. Thus a polynomial excess for this whole candidate family must
come from \(E_{\rm rem}\); the remaining incidence estimate is not
proved here.

These are exact reductions using classical greedy compression, duality
and the saved shadow/sweep estimates. They do not prove a new full
hook count, survival of support-matching bases, an unrestricted exponent
or current novelty. The family \(a=t^2+1\) has \(n=\Theta(t^3)\),
so the conditional decision rank is \(O(n^{2/3})\), rather than the
linear overlap rank. The reduced matroid can change between regions;
there is no single global rank-\(3t^2\) replacement.

## The original matrix and its three bridge blocks

Write \(D=\{(i,j):0\le i,j<a,\ i<t\text{ or }j<t\}\).
The adjacent monomial syzygies give the complete allowed supports
\[
\begin{array}{c|c}
 f_0=\phi(y^a)&i\ge t\text{ or }j\ge a-t\\
 f_1=\phi(y^tz^t)&i\ge t\text{ or }j\ge t\\
 f_2=\phi(z^a)&i\ge a-t\text{ or }j\ge t,
\end{array}                                                     \tag{H2}
\]
inside \(D\), with independent coefficients. The proof in the
[generic-hook note](generic_hook_fans.md) applies unchanged to
\(a\ge2t\): the two sides of each adjacent syzygy have disjoint
possible surviving supports, and hence each side is zero.

The plain \(D\) columns are compulsory coloops. Every other possibly
nonzero plain column maps into
\[
 T=(y^t,z^t)C=T_Y\oplus T_Z,
 \quad T_Y=y^tC,\quad T_Z=z^tC,
 \quad \dim T_Y=\dim T_Z=l=t(a-t).                              \tag{H3}
\]
The upper monomials \(xy^iz^j\) with \(i,j<t\) are also coloops:
they alone have coordinates in the complementary \(t^2\)-dimensional
box. Contract these two literal fixed sets. The residual rank is
\(2l\), with upper identity columns on \(T_Y,T_Z\) and the plain
columns given by shifted images of \(f_0,f_1,f_2\).
This rank is an encoding rank. Under exceptional specialization the
actual overlap \(B1\cap Bx\) may be smaller than \(T\).

Only the following plain columns can have coordinates in both arms:
\[
\begin{split}
 H_0&=\{(0,a+i,j):0\le i,j<t\},\\
 H_1&=\{(0,t+i,t+j):0\le i,j<t\},\\
 H_2&=\{(0,i,a+j):0\le i,j<t\}.                         \tag{H4}
\end{split}
These are disjoint and have total size \(k=3t^2\). They are possible
bridges in the independent-parameter matrix. A specialization may
make a bridge one-sided or zero; it remains in this third group.

For example, a shifted \(f_0\) has a surviving \(Y\)-component
only for shifts \(0\le i<a-t,0\le j<t\), and a surviving
\(Z\)-component only for \(i,j<t\). Their intersection gives
\(H_0\). The other two calculations give \(H_1,H_2\).
This argument also lists every possibly nonzero column outside \(D\):
the monomial carrier is \(M^2\), and its three disjoint source regions
are the shifts of the three heads just described.

Call the remaining residual groups \(E_Y,E_Z\), according to their
possible coordinate arm, and put \(H=H_0\sqcup H_1\sqcup H_2\).
Their exact sizes and ranks are
\[
 |E_Y|=|E_Z|=3l-2t^2,\quad |H|=3t^2,\qquad
 r(E_Y)=r(E_Z)=l,\quad r(E_Y\cup E_Z)=2l.                       \tag{H5}
\]
The ranks hold at every specialization because the upper identity
columns remain. The two pure spans are disjoint. Their label boxes are
\[
 E_Y\subset[0,2)\times[0,2a)\times[0,2t),\qquad
 E_Z\subset[0,2)\times[0,2t)\times[0,2a).                       \tag{H6}
\]
All labels in this note are original exponent triples. Row operations
do not change them, and the central-generator image \(f_1\) is retained.

## Exact conditional dual compression

Fix the internal cost orders in \(E_Y,E_Z,H\). Greedily compress
\(E_Y,E_Z\) to independent ordered bases \(B_Y,B_Z\), of size \(l\)
each. Any discarded column was a combination of cheaper columns in
its own group and can never enter a merged greedy minimum. Retain all
of \(H\), including specialized loops. The restricted matroid has
\(2l+k\) columns and rank \(2l\), hence its dual has rank \(k\).

A minimum original basis is the complement of a maximum-cost dual
basis. In each of the three fixed ordered lists \(B_Y,B_Z,H\), scan
in descending cost and retain only dual-independent columns. This
leaves at most \(k\) from each list. A discarded dual column is in
the span of more expensive columns in its own list, so cannot enter
a global maximum dual basis. The resulting fixed rank-\(k\) matroid
on at most \(3k\) columns therefore recovers every original minimum
throughout the whole internal-order region, by taking its complement
and restoring the fixed columns.

This proves equality of support functions up to the explicit linear
cost of the restricted ground. It proves entire fan patches, including
their boundaries by continuity, rather than only agreement at sampled
weights. In a representation using \(B_Y\sqcup B_Z\) as reference,
write the bridge columns as a matrix \(F\). A dual representation is
\([ -F^{\mathsf T}\mid I_k]\). No determinant is replaced by matching.

There is a common internal-order inventory of \(O(at)\) planes:
(H6) gives \(O(at)\) difference vectors in each arm, while (H4)
gives \(O(t^2)\) between the three translated bridge boxes. Thus
\(O((at)^2)=O(n^2)\) internal-order chambers suffice. The additional
subdivision made by the small dual fan is still unbounded at the
desired total scale. Counting these chambers alone does not count the
original initials.

## The known directions already have a quadratic cell budget

Define \(\mathcal N\) exactly as in (H1), using possible labels before
specialization and orienting cross-layer differences to first component
one. Such differences are already primitive. Differences within either
box in (H6) give \(O(at)\) vectors. For a bridge box and an arm box,
one transverse difference coordinate ranges over \(O(a)\) values and
the other over \(O(t)\), irrespective of the bridge box's translation.
There are only six such box pairs. This proves \(|\mathcal N|=O(at)\).
It is an exact pair inventory, not the enlarged union of its enclosing
boxes and not a list of input coefficient weights.

All residual labels satisfy
\((x+1)(y+1)(z+1)\le8at=O(n)\) and \(x<2\).
The [translated-layer shadow theorem](layer_shadow_budget.md), applied
directly to this represented matroid of rank \(2l\le n\), bounds
each relevant rational shadow and each actual parallel multiplicity by
\(O(n\log(n+1))\). This use does not pretend that the contracted
matroid is an algebra of length \(2l\).

The [finite polytope sweep](cross_layer_direction_budget.md) bounds the
number of positive vertices by the zero-first-coordinate shadow plus
the number of actual nonzero-first-coordinate edge pieces. The shadow
has \(O(n\log(n+1))\) vertices. Summing the parallel bound over
\(\mathcal N\) gives \(O(n^2\log(n+1))\) known-direction pieces,
proving (H1). If a remaining edge had a same-arm or bridge exchange
direction, it would lie in \(\mathcal N\). Projected matroid-edge
exchange therefore gives the asserted pure cross-arm localization.

The total primitive one-norm length of the known-direction pieces is
also
\[
                O(a n^2\log(n+1))=O(n^3\log(n+1)),             \tag{H7}
\]
since every difference has one-norm \(O(a)\). This covers some
opposite-sign transverse directions as well as same-sign ones; it is
not a full collective wall budget for the family.

The existing box ceiling \(g(I)=O(n^{3/2}a)\) remains valid, since
\(x^2,y^{2a},z^{2a}\) annihilate \(A_\phi\). In particular
\(a=t^2+1\) gives \(O(n^{13/6})\). The present reduction supplies
no stronger unconditional whole-family exponent.

## Verification and boundary

```bash
python3 -B scripts/run_hook_bridge_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker uses actual field ranks, reconstructs the original quotient
columns and compares full cone patches from the small dual with those
from the original matrix. It includes zero maps, small characteristics,
generic fixtures and the existing sparse convolution rank failure.
The 89-cone generic fixture checks one strict internal-order patch per
original cone, as well as full coverage by the original cones. The
uniform compression proof supplies the other internal-order patches.
Finite calibrations are not extrapolated into a growth law.

Greedy compression and duality are classical; the same two-stage
principle is used in the [univariate-cross proof](univariate_cross_fans.md).
The state-polytope dictionary, projected-edge exchange and shadow inputs
retain the sources and priority limits in their linked proofs. The
independent-parameter construction and its finite generic certificates
are unchanged. Exact priority of this combined localization remains
unassessed. A uniform bound for the remaining cross-arm pieces, or an
actual infinite polynomial separation, is still required by the goal.
