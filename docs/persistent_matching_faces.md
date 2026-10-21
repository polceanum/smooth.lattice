# Simultaneous cores at persistent matching faces

The [quadratic-coordinate reduction](quadratic_coordinate_fan_compilation.md)
shows that a polynomial separation in this class would require many
centrally symmetric compact state faces. Here is a further uniform
constraint on those faces, together with a limit of the planar-flag
charging approach. Neither result bounds their total number.

Throughout, \(A=K[x,y,z]/I\) has length \(n\), over any field, and
\(I\) contains a monic univariate quadratic in \(x\). Put
\(u=w_y/w_z\), \(t=w_x/w_z\). A positive exposed compact two-face
is **persistent** if its normalized weight is not an endpoint of any
positive literal matching run in the two-list reduction. At most
\(O(n^2\log(n+1))\) positive fan vertices are excluded by this
condition, including excluded vertices whose dual face happens to be
centrally symmetric. A persistent face is a zonogon by the existing
local support-function proof.

## One common compulsory core for all its wall directions

For a primitive mixed direction \(q\), define
\[
 C_q=(\operatorname{Div}(q^+)\cup\operatorname{Div}(q^-))
                                             \setminus\{q^+,q^-\}.
                                                               \tag{P1}
\]
Let \(\mathcal B(F)\) be the actual initial standard sets whose
profiles are vertices of a persistent face \(F\). Set
\[
 C=\bigcap_{B\in\mathcal B(F)}B,\qquad
 K_F=\bigcup_{q\text{ an edge direction of }F} C_q,
 \qquad r_F=n-|K_F|.
\]
Then
\[
                              K_F\subseteq C.            \tag{P2}
\]
This is literal monomial containment and simultaneous independence in
the original quotient. It improves the individual-wall assertion by
using the same core for the whole exposed polygon. The vector span of
this core is not asserted to be an algebra ideal.

To prove it, near the face weight every monomial's membership in the
selected standard basis has one of three forms. It is constant; it
switches across the line of one persistent matched pair; or it switches
across the one possible vertical planar-order wall. This follows from
the fixed literal tied-pair set and the two fixed unmatched selections
on the sides of a planar order event. A switching indicator is constant
on each of the two open half-planes of its line. Coincident pair lines
may switch several monomials together.

The [compulsory-core theorem](wall_coexistence.md) puts \(C_q\) in
both endpoint standard sets of every actual \(q\)-edge. A zonogon has
two opposite edges in direction \(q\); their dual wall arms point in
opposite directions from the face weight. If a monomial of \(C_q\)
switched on any different line through that point, those two arms would
lie on opposite sides of its switching line. One \(q\)-edge would then
lack it, a contradiction. Nor can it switch on the \(q\)-line itself,
since each \(q\)-edge contains it on both sides. Thus every element of
\(C_q\) is constant across the entire local fan and belongs to \(C\).
This proves (P2), including a possible vertical direction.

Normalize each nonvertical primitive direction to \(q_x=1\), and
let \(m_q\) be the positive integer length of either parallel profile
edge in this primitive direction. If \(d\) distinct nonvertical
directions occur, then
\[
 \sum_{q_x=1}m_q
       =\max_{B\in\mathcal B(F)}\mu_x(B)
         -\min_{B\in\mathcal B(F)}\mu_x(B)
       \le r_F,\qquad d\le r_F.                         \tag{P3}
\]
The equality is the ordinary width formula for a zonogon. After fixing
the common independent monomials \(K_F\), each remaining basis has
\(r_F\) monomials with \(x\)-exponent zero or one. This proves the
inequality without summing residual ranks of different contractions.

## A sharp square-root constraint on concurrence

There is a second constraint independent of core size:
\[
                   d(d+1)\le n.                         \tag{P4}
\]
If \(\epsilon\in\{0,1\}\) records a vertical wall, the number of
incident full cones, equivalently the number of vertices of \(F\), is
therefore
\[
 |\mathcal B(F)|=2(d+\epsilon)
   \le 2\min\left\{r_F,
       \left\lfloor\frac{\sqrt{4n+1}-1}{2}\right\rfloor\right\}+2.
                                                               \tag{P5}
\]
There is at most one vertical direction because it has equation
\(u=u_0\). The additive two in this displayed bound only allows that
direction; no claim of its sharpness is made.

Choose one upper-layer monomial from an active literal pair for each
nonvertical direction. These \(d\) monomials are pairwise incomparable
under divisibility. Otherwise the half-plane where a multiple is
selected would have to be contained in the half-plane where its divisor
is selected, since every actual initial standard set is a downset.
Distinct lines through an interior point have no such local half-plane
containment. Coincident lines were already grouped into one direction.

Decrease \(t\) slightly from the face weight. All of these upper
representatives are then selected simultaneously. If the point is on a
vertical event, a still smaller change in \(u\) selects one side while
preserving all these choices. The upper planar diagram \(D_1\) of that
initial basis contains an antichain of size \(d\).

The following elementary downset estimate supplies the constant in
(P4). Order any planar antichain by increasing first coordinate. Its
coordinates are at least \((i,d-1-i)\), for \(0\le i<d\), because
the second coordinates strictly decrease. Its divisor downset therefore
contains the full triangle of total degree at most \(d-1\), of size
\(d(d+1)/2\). Also \(D_1\subseteq D_0\), since deleting \(x\)
from a standard monomial preserves standardness. Hence
\(n=|D_0|+|D_1|\ge2|D_1|\ge d(d+1)\), as required.

The estimate concerns distinct directions, not the number of literal
pairs. Translated pairs on the same line may have a large total
multiplicity; (P3) retains their integer edge lengths.

## Equality in the concurrence bound over every field

For every integer \(d\ge2\), let
\[
 D_1=\{(a,b)\ge0:a+b\le d-1\},\qquad
 D_0=D_1\cup\{(d-i,i):0\le i<d\}.
\]
Let \(J_D\) be the monomial ideal complementary to the three-variable
downset \(D=(\{0\}\times D_0)\sqcup(\{1\}\times D_1)\), and set
\[
 I_d=J_D+
  (xy^iz^{d-1-i}-y^{d-i}z^i:0\le i<d).                 \tag{P6}
\]
All \(2d\) displayed monomials are distinct maximal monomials of
\(D\), so their classes are socle elements of the monomial quotient.
The relations are disjoint, independent identifications in that socle.
Their coefficients are units over every field. Thus
\[
                 \dim K[x,y,z]/I_d=d(d+1).              \tag{P7}
\]
Every initial standard set is \(D\) with the expensive endpoint of
each pair removed: this gives a downset and already has the correct
colength. There are no further comparisons.

The primitive directions are
\[
                  q_i=(1,2i-d,d-1-2i),\quad0\le i<d.
\]
They give distinct lines
\[
                   t=1+(d-2i)(u-1),                    \tag{P8}
\]
all meeting at \((t,u)=(1,1)\). Their full positive arrangement has
exactly \(2d\) cones and one interior vertex. Its dual polygon is a
zonogon with \(d\) nonvertical directions. The literal matching pairs
persist for all \(u>0\), so the vertex is persistent even though both
planar list orders change at \(u=1\). Equality holds in (P4).

This proves the square-root order cannot be improved. The family has
only \(2d=\Theta(\sqrt n)\) initial ideals and supplies no lower
separation for the original extremal problem.

## Ordinary crossings can lie strictly inside fixed planar flags

The entire face count cannot be charged locally to changes in the
two planar flags. For \(b\ge3\) with \(3\nmid b+1\), use the
two-layer axis hooks
\[
 D_0=\{1,y,y^2,z,z^2,\ldots,z^b\},\qquad D_1=\{1,y,z\},
\]
and identify the disjoint socle pairs
\[
                           xy=z^b,\qquad xz=y^2.        \tag{P9}
\]
As above this is an actual quotient with \(x^2=0\), of length
\(n=b+4\), over every field. Its full fan has four cones, separated
by \(t=b-u\) and \(t=2u-1\), with crossing
\[
                  u=(b+1)/3,\qquad t=(2b-1)/3>0.        \tag{P10}
\]
At this weight both complete ordered planar bases are strict. Indeed
the reduced ratio for \(u\) has numerator \(b+1>b\) and denominator
\(3>2\). No two monomials in \(D_0\), or in \(D_1\), can have
equal planar cost. Their bases are independent monomial hooks, so
neither ordered flag changes anywhere in a neighborhood. The smallest
example has length seven and crossing \((t,u)=(5/3,4/3)\).
Its joint compulsory core has size \(n-2\), and its profile polygon
has \(x\)-width two, attaining (P3).

This refutes assigning every central crossing to a planar-prefix event
at the same \(u\). It does not refute an aggregate quadratic bound
using information from other parameter values.

## Consequence for the open target and prior-art boundary

For the quadratic-coordinate branch, all but
\(O(n^2\log(n+1))\) fan vertices obey (P2)--(P5). High-order
concurrence costs at least quadratic length in its number of distinct
directions, and there is a common rank-sized residual at each such
face. These facts do not bound the number of ordinary crossings with
two directions. Controlling that number, including coexistence in one
quotient, remains the essential counting problem. The unrestricted
three-variable endpoint and the goal's success criteria are unchanged.

Two-list matroid matching, zonogon width and opposite-edge symmetry,
downset antichain estimates, and socle identifications are established
tools. The existing [compulsory-core proof](wall_coexistence.md) and
[two-flag proof](quadratic_coordinate_fan_compilation.md) supply the
precise inputs. The [primitive-wall note](primitive_wall_spectrum.md)
already compares the socle move with Hering--Maclagan's edge-ideal
framework. No new general socle construction or zonotope principle is
claimed. Exact priority of the persistent-face combination, sharp
constant and fixed-flag obstruction remains unassessed.

## Verification

~~~bash
python3 -B scripts/run_persistent_face_checks.py
~~~

The report `results/research/persistent_face_checks.json` records
complete actual fans, exact profile polygons, simultaneous core
containments, antichain representatives, residual widths and valences.
The checker reconstructs arbitrary small coupled and quadratic examples
independently from border normal forms, identifies persistent faces
using the two-flag compiler, and compares the explicit families with
their entire predicted line arrangements. It independently checks
commuting multiplication matrices and all ground monomial successors
in the two new socle families. It tests the strict-flag
crossings and equality in (P4) in several fields. Uniform conclusions
use the proofs, not a growth fit.
