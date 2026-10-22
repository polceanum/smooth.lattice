# Large-product hook walls retain only the central image

The [central-prefix localization](hook_central_prefix.md) removes outer
source exchanges on an explicit nonempty open locus. A further exclusion
removes both outer images from every remaining exchange whose two
transverse displacements have large product. This identifies its geometric
wall locus with that of the central-only extension. It also transfers the
distinct cross-layer direction problem to that smaller class. Parallel
subdivisions and prefix-degenerate strata remain unresolved.

Use the same compatible flat hook extension, over any field:
\[
 M=(y^a,y^tz^t,z^a),\qquad a\ge2t,\quad t\ge2,\qquad
 n=4at-2t^2.
\]
Fix a central image \(f_1\) satisfying the prefix-regular condition (CP3).
The two outer images \(f_0,f_2\) are arbitrary. Write \(I_\phi\) for
the extension ideal, and \(I_\circ\) for the extension with the same
\(f_1\) and both outer images zero. They have the same length. In
particular the latter ideal has the presentation
\[
 I_\circ=(x^2,y^a,z^a,y^tz^t-xf_1).                         \tag{CS1}
\]
The redundant products \(xM,M^2\) follow by multiplying these
relations, using the compatibility of \(f_1\).

## The remaining convolution has a fixed monomial source

For \(a>2t\), prefix regularity implies that the coefficients of
\(y^t,z^t\) in \(f_1\) are nonzero. This subsection needs only those
two nonzero coefficients, including on prefix-degenerate strata. Put
\(A=a-t\), and write
\[
 f_1=y^tU+z^tV,\qquad
 R_Y=K[y,z]/(y^A,z^t),\quad R_Z=K[y,z]/(y^t,z^A).
\]
The two components of \(f_1\) lie on disjoint arms; \(U\in R_Y\)
and \(V\in R_Z\) are units. Divide the two row modules by their
tags \(y^t,z^t\) and multiply by \(U^{-1},V^{-1}\), respectively.
This is an invertible \(K[y,z]\)-linear row transformation. It changes
multiplication by \(f_1\) into the diagonal map
\[
 h\longmapsto(h\bmod(y^A,z^t),\ h\bmod(y^t,z^A)).
\]
Its kernel and cokernel follow from the elementary intersection sequence:
\[
\begin{split}
 J&=(y^A,z^t)\cap(y^t,z^A)=(y^A,y^tz^t,z^A),\\
 0\longrightarrow K[y,z]/J&\longrightarrow R_Y\oplus R_Z
   \xrightarrow{\ (v,w)\mapsto v-w\ }
       K[y,z]/(y^t,z^t)\longrightarrow0.                 \tag{CS1a}
\end{split}
\]
The last map restricts both components to the common \(t\)-by-\(t\)
box. Equality on that box is exactly the condition for gluing the two
rectangular monomial polynomials, which proves exactness directly.
Consequently
\[
 \operatorname{Ann}_{K[y,z]}(f_1)=J,\qquad
 \dim(f_1C)=2t(a-t)-t^2.                                 \tag{CS1b}
\]
In particular every central-source plain column is independent of the
other central-source columns. The upper identity columns become their
actual inverse-unit convolution vectors under this row transformation;
their exponent labels do not change. Thus this normal form preserves
all ranks but does not make the full monomial matroid independent of
\(f_1\). The intersection sequence and unit argument are classical
algebraic inputs, not a new exact-sequence theorem.

## The second outer head is also too expensive

Consider a central pure lower \(Y\)-column and upper \(Z\)-column
\[
 e=(0,Y,Z),\quad 2t\le Y<a,\ t\le Z<2t,\qquad
 f=(1,R,S),\quad 0\le R<t,\ t\le S<a.
\]
Put \(b=Y-R\), \(c=S-Z\), and suppose
\[
                         bc>16at.                         \tag{CS2}
\]
Here an actual positive exchange necessarily has \(b,c>0\).
Normalize its weight to
\[
 w=(\lambda,u,1),\quad\lambda=bu-c>0,\qquad
 \tau=Yu+Z=\lambda+Ru+S.
\]
Since \(b,c<a\), (CS2) implies \(b,c>16t\). The earlier (CP14)
already proves \(\tau<au\).

**Lemma.** Every such actual exchange also satisfies
\[
                              \tau<a.                     \tag{CS3}
\]

Suppose \(\tau\ge a\). The central pure \(Z\)-tail columns span
\(z^{2t}C\), using the nonzero coefficient \([z^t]f_1\). As in
(CP8), they would make \(f\) a loop after contraction of cheaper
columns unless
\[
 \tau\le(2t-1)u+a-1,\qquad
 u\le\frac{a-1-Z}{Y-2t+1}
       <\frac{a}{b-2t}<\frac{4a}{3b}.                     \tag{CS4}
\]
Set \(r=t+\lfloor2tu\rfloor\). Since \(b>8t\), (CS4) gives
\(r<t+a/3<a-t\), while \(r\ge t\). The transposed prefix
basis \(P_Z(r)\cup Q_Z(r)\) from (CP1) is strictly cheaper
than \(f\). Indeed
\[
\begin{split}
 \tau-\max_{Q_Z(r)}w\cdot v
   &\ge r+1-t-(2t-1)u>0,\\
 \tau-\max_{P_Z(r)}w\cdot v
   &=S-(t-1-R)u-t-r+1\\
   &\ge S-3tu-2t+1>0.                                   \tag{CS5}
\end{split}
\]
For the last inequality, \(3tu<4at/b<c/4\), and \(S\ge c>16t\).
These actual cheaper columns span all of \(T_Z\), including \(f\),
which contradicts exchangeability. This proves (CS3). All strict
inequalities also hold at weights with additional ties elsewhere.

Combining (CS3) and (CP14), the common cost satisfies
\[
                       \tau<\min(au,a).                   \tag{CS6}
\]
The transpose proves the corresponding assertion in the other
cross-arm orientation. The constant sixteen is a sufficient constant;
its optimization is not needed for the reduction.

## Literal eligibility and geometric wall loci agree

Every column from the outer source \(y^a\) has cost at least \(au\),
and every column from \(z^a\) has cost at least \(a\). Below both
heads, all vectors at the exchanged cost and all strictly cheaper
vectors therefore depend only on \(f_1\). At a generic point of an
exchange plane, a literal pair is eligible exactly when its two images
are nonzero and parallel modulo the span of cheaper columns. Hence its
eligibility is the same for \(I_\phi\) and \(I_\circ\). If the
common cost reaches either outer head, the pair is ineligible for every
choice of outer images by the two uniform exclusion lemmas. Thus this
comparison works in both directions.

Let \(\mathcal N\) be the exact same-arm/bridge cross-layer inventory
in the [bridge reduction](hook_bridge_reduction.md). Define
\[
\begin{split}
 \mathcal H=\mathcal N\ \cup\ \{
  (1,d,e):\ |d|,|e|<2a,\quad |de|\le16at\}.              \tag{CS7}
\end{split}
\]
All directions here have first coordinate one, so are primitive. An
actual cross-layer direction outside \(\mathcal H\) can only come
from a pure cross-arm pair. If it came from an outer source, (CP6)
would give \(|de|<12at\), a contradiction. It therefore comes from
a central source. Its transverse entries have opposite signs and
magnitudes greater than \(8t\), and it satisfies (CS2).

For each \(q\notin\mathcal H\), define its geometric wall locus as
the closure, within its positive exchange plane, of the relative open
wall intervals. Subdivision points are included in this closure and are
not counted separately. Then
\[
          \mathcal W_q(I_\phi)=\mathcal W_q(I_\circ),
                    \qquad q\notin\mathcal H.             \tag{CS8}
\]
To verify this from literal eligibility, avoid the finitely many
additional, nonparallel projected-label ties. At a fixed cost there
are at most two labels tied along such a plane: their difference is
\(q\), whose first coordinate is one, and labels have first coordinate
zero or one. A nontrivial exchange at any such cost gives a positive
profile-edge length in direction \(q\); maximizing and minimizing the
upper count cannot cancel that difference. The rank tests above thus
give identical presence or absence of a wall on every remaining open
interval. Taking closures proves (CS8), including recurrent intervals.

Equation (CS8) compares the union of geometric wall intervals. It does
not identify their subdivision into polytope edges. More expensive
columns can change the rest of the minimum basis, introducing vertices
on a persistent wall. In particular, neither equality of parallel piece
counts nor an additive comparison of the two full fan counts follows.

The exception inventory cannot simply be omitted. At the saved
[length-102 recurrent-pair fixture](hook_wall_recurrence.md), the central
image is prefix regular over \(\mathbb F_{1000003}\). The pair
\(e=(0,7,3)\), \(f=(1,1,7)\) is eligible at
\(w=(15/14,71/84,1)\) for the full extension, but becomes ineligible
when both outer images are set to zero. Its direction \((1,-6,4)\)
lies in \(\mathcal H\). This is an actual rank comparison at one
field point, not an asymptotic obstruction to (CS11).

## A direction gate for the central-only class

Put \(B_0=2a-1\), \(Q=16at\). The number of signed pairs in the
second set of (CS7) is exactly
\[
 4B_0+1+4\sum_{b=1}^{B_0}\min(B_0,\lfloor Q/b\rfloor)
                  =O(at(1+\log(a+1))).                   \tag{CS9}
\]
The first term accounts for the two axes. The harmonic-sum estimate
proves the bound uniformly, including when the whole square is included.
Together with \(|\mathcal N|=O(at)\) and \(at=\Theta(n)\), this gives
\( |\mathcal H|=O(n(1+\log(n+1)))\).

Let \(\mathcal D_x(I)\) be the set of actual cross-layer directions,
and \(D_x(I)=|\mathcal D_x(I)|\). The stronger set statement is
\[
 \mathcal D_x(I_\phi)\mathbin{\triangle}\mathcal D_x(I_\circ)
                        \subseteq\mathcal H.              \tag{CS10}
\]
Every direction in \(\mathcal H\) has
\(O(n(1+\log(n+1)))\) actual parallel pieces by the translated-layer
shadow theorem. Thus its entire controlled piece count is
\(O(n^2(1+\log(n+1))^2)\). Alternatively the
[cross-layer direction theorem](cross_layer_direction_budget.md) yields
\[
 g(I_\phi)=O\!\left(
 n^2(1+\log(n+1))^2+
 n(1+\log(n+1))D_x(I_\circ)\right).                       \tag{CS11}
\]
All parallel pieces remain included through that theorem.

Consequently a polynomial excess for a prefix-regular full hook family
requires a superlinear direction inventory already in the central-only
class (CS1). For example, \(g(I_\phi)\ge c_0n^{2+\epsilon}\) along
an infinite sequence implies
\(D_x(I_\circ)=\Omega_{c_0,\epsilon}(n^{1+\epsilon}/\log(n+1))\)
for sufficiently large lengths. A near-linear bound for those directions
is not proved here. The earlier \(O(n^{3/2}a)\) box ceiling, the
prefix-degenerate strata, and the unrestricted extremal endpoint remain
unchanged.

## Verification and priority boundary

```bash
python3 -B scripts/run_hook_central_skeleton_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker reconstructs original commuting quotient columns, verifies
the prefix bases, tests the strict opposite-tail and opposite-prefix
certificates, and compares exact eligibility after changing both outer
images. Separate small-characteristic controls use arbitrary central
coefficients with the two leading units to check the monomial annihilator,
inverse-unit normal form and common-box cokernel. An active-pair control
checks the stated exception-inventory guard. Arithmetic controls independently enumerate the signed
hyperbolic inventory. Finite inactive-pair samples do not prove a uniform
central-wall absence theorem or a growth law.

The prefix, greedy-closure, state-polytope and shadow inputs retain the
sources and limitations in their linked proofs. The divisor-sum estimate
is elementary. Exact priority of the combined wall-locus reduction is
unassessed. No new unrestricted exponent, infinite polynomial separation,
algorithmic speedup, or completed extremal goal is asserted.
