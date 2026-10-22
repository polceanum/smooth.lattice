# Profile sections: the geometric relaxation and the fixed-pair boundary

The open collective problem in [the incidence note](ideal_step_incidence.md#a-collective-section-problem-without-an-ideal-chain-recurrence)
is to bound actual positive profile-edge crossings of half-integer
coordinate planes. A near-linear bound at every level would imply a
cubic polylogarithmic actual-length budget. This note audits two tempting
relaxations. Neither proves that bound or resolves unrestricted \(G_3\).

For an actual length-\(n\) quotient, write
\[
 P_I=\operatorname{conv}\{\mu(B): B\text{ a positive initial standard set}\},
 \qquad c_i(k)=\#\{e:\min_e x_i<k+1/2<\max_e x_i\}.
\]
Here \(e\) ranges over actual positive profile edges. Integer endpoints give
\[
             \mathcal L(I)=\sum_{i=1}^3\sum_k c_i(k).
                                                        \tag{PS1}
\]
Crossings are distinct vertices of the genuine rational section of
\(P_I\). This established reduction does not identify the section with
a planar algebra quotient.

## General matroids do not give a rank-only near-linear section bound

[Rörig--Witte--Ziegler, Theorems 1.2 and 4.4, equation (2)](https://arxiv.org/pdf/0710.3116v3)
construct three-dimensional zonotopes with a quadratic number of vertices
in a central planar section. This remains a counterexample to a proposed
*rank-only* near-linear section estimate for arbitrary represented matroids:
a sum of \(r\) segments is a projection of the rank-\(r\) matroid that
chooses one of two parallel columns in each of \(r\) independent classes.
The result is established convex geometry, not a new construction here.

The integer labels and quotient length must also be retained. The following
explicit normalization checks why the published example does not refute
the proposed bound in the actual colength parameter.

Let \(k\ge1\), \(m=4k+1\), \(L=m+1\), and put
\[
 b=(k-i)_{i=0}^{2k},\qquad
 b'=(i-k+1/2)_{i=0}^{2k-1}.
\]
Multiply every row of their matrix (2), for \(d=3\), by the same
\(2L^2\). The resulting \(2m\) integer generators are
\[
 \begin{array}{ll}
 (2L^2b_j,\ 2L^2,\ 2L^2),&
 (2L^2b'_j,\ -2L^2,\ 2L^2),\\
 (2b_j,\ 0,\ 2L),&
 (2b'_j,\ 0,\ -2L).
 \end{array}                                             \tag{PS2}
\]
Take \(Z=\sum_g[0,g]\). This differs from the source's centrally
symmetric convention only by a translation and uniform homothety.
Its central third-coordinate level is
\[
                        t=mL^2+L.
\]
Every vertex has third coordinate a multiple of \(2L\), whereas \(t\)
is an odd multiple of \(L\). Thus no vertex occurs between \(t\) and
\(t+1/2\), and the same edges cross both planes. The source theorem
gives at least \(2m(m+1)\) section vertices. The general quadratic
face bound gives \(\Theta(k^2)\) altogether. In centered coordinates,
fixing the first \(m\) segment choices contributes an odd multiple of
\(L^2\) to the third coordinate; the other \(m\) segments can vary it by
at most \(mL=L^2-L\). Thus no edge from the second family crosses this
level. The first family's edge directions are not coordinate axes, so
their two-dimensional normal cones contain vectors with no zero
coordinate. Every crossing edge is therefore exposed in at least one
open cost orthant. Hence some sign change
of the three coordinates has \(\Omega(k^2)\) *positive* crossings.
This is an abstract matroid statement, not an algebra realization.

In fact **every** such cost orthant has positive profile width at least
\[
                 L^2 k(2k+1)=\Theta(k^4)                 \tag{PS3}
\]
in the first coordinate. To see this, compare costs approaching its
first and third coordinate axes. In the first two blocks the first
coordinates sum to zero, and their total absolute sum is
\(2L^2\{k(k+1)+k^2\}\). The first-axis minimum chooses their negative
half. The third-axis minimum chooses either all or none of those blocks,
since their third coordinates have the same sign. The difference in the
first moment is half their absolute sum. Each remaining block contributes
a nonnegative difference from the first-axis minimum. Sign changes do
not affect this argument. Sufficiently small positive perturbations make
the two cost vectors strict while retaining their first moments.

Every actual length-\(N\) downset has
\(0\le\mu_x\le\binom N2\). Any algebra with these *displayed integer
profiles*, even up to translation, would therefore require
\[
                   \binom N2\ge L^2k(2k+1),
                   \qquad N=\Omega(k^2).                 \tag{PS4}
\]
The abstract rank \(2m=O(k)\) is insufficient. At that necessary length
scale the displayed example's \(\Theta(k^2)\) section complexity is only
\(O(N)\). It supplies no polynomial counterexample to the actual
near-linear section target. This does not exclude another normalization
or another algebra construction; any such construction needs its own
integer-profile and multiplication checks.

## A global fixed pairing cannot give a supergeneric fan

There is a separate uniform limit on the simplest paired-sequence route.
Suppose an actual length-\(n\) quotient has disjoint monomial pairs
\(\{u_j,v_j\}\), \(1\le j\le a\), and a common fixed monomial set
\(C\), such that **every generic positive cost** has standard set
\[
 B_w=C\ \sqcup\ \{\text{the cheaper of }u_j,v_j:1\le j\le a\}.
                                                        \tag{PS5}
\]
The pairs and \(C\) are fixed over the entire positive cost domain.
Then \(a\le n\), and classical line-arrangement geometry gives
\[
       g(I)\le1+\frac{a(a+1)}2,
       \qquad \mathcal L(I)\le2a^2(n-1)\le2n^2(n-1).
                                                        \tag{PS6}
\]
These are direct classical corollaries, not a new general matroid theorem.

Indeed the profile support function in (PS5) is a fixed linear form plus
\(\sum_j\min(w\cdot u_j,w\cdot v_j)\). It is the support function
of a translated zonotope. Pairs with componentwise comparable endpoints
never switch and can be absorbed into \(C\). Every other endpoint is
selected at some positive cost, so belongs to an actual standard downset.
Its total degree is at most \(n-1\), by its divisor chain. Consequently
\[
                       \|u_j-v_j\|_1\le2(n-1).
\]
Let \(p\le a\) be the number of distinct switching planes. In the
normalized positive triangle each is one line, cut into at most \(p\)
segments by the other lines. Coincident pair directions switch together;
their oriented displacement vectors are parallel with the same sign, so
their actual edge lengths add. Thus
\[
 \mathcal L(I)\le p\sum_j\|u_j-v_j\|_1\le2a^2(n-1).
\]
The usual \(p\)-line cell bound proves the first part of (PS6), including
coincidences and missing intersections. Constant pairs only reduce the
bounds; a fan with one cone has no positive edges.

This closes the whole fixed-pair class as a source of a polynomial fan
separation, regardless of the spread of its exponent labels. It includes
the established unaugmented parity-grid construction. It does **not**
apply to the added-point construction, where the doubled pair changes,
or to a representation whose matching changes between planar strips.
The latter changing interactions remain relevant to the unrestricted goal.

## Binary-coordinate sections: integral, but usually not matroids

For a quotient with a monic scalar quadratic in \(x\), retain the finite
ground \(E=\{0,1\}\times V_n^2\), where
\(V_n^2=\{(a,b)\ge0:(a+1)(b+1)\le n\}\). It contains every
positive standard set. Let \(M\) be its actual rank-\(n\) residue
matroid, and define the **full** basis-profile polytope
\[
 Q=\operatorname{conv}\{\mu(S):S\text{ any basis of }M\},\qquad
 Q_k=\operatorname{conv}\{(\mu_y(S),\mu_z(S)):
                                    |S\cap E_1|=k\}.     \tag{PS7}
\]
Here \(E_1\) is the upper layer. These bases need not be downsets.
The distinction between \(Q\) and the positive-profile hull \(P_I\)
in (PS1) is essential. Their positive exposed faces agree; their
unexposed sections need not agree.

The standard matroid-polytope edge theorem gives the exact identities
\[
 \begin{split}
 Q\cap\{x=k\}&=\{k\}\times Q_k,\\
 Q\cap\{x=k+\theta\}
   &=\{k+\theta\}\times\bigl((1-\theta)Q_k+\theta Q_{k+1}\bigr),
                 \qquad 0<\theta<1,                      \tag{PS8}
 \end{split}
\]
whenever the two levels are feasible. These are classical corollaries,
not a new integrality theorem. Apply the colour functional
\(\chi(S)=|S\cap E_1|\) to the unprojected base polytope. Its edges
are single basis exchanges, so their endpoint colour counts differ
by at most one. An integer slice therefore has only original vertices.
A slice at \(k+\theta\) has vertices dividing original edges with
counts \(k,k+1\) in the proportions \(1-\theta,\theta\). This proves
one inclusion in (PS8); convexity proves the reverse inclusion.
Projecting the remaining coordinates proves the displayed identities.
In particular its half-integer sections lie in the half-integer planar
lattice, and
\[
 c_x(k)\le f_0(Q_k)+f_0(Q_{k+1}).                         \tag{PS9}
\]
Every actual positive crossing edge remains an edge of \(Q\). Planar
Minkowski addition has at most the sum of the input vertex counts,
including segments and points. These facts justify (PS9); inclusion
of arbitrary hulls alone would not justify a vertex-count comparison.

Let \(p=\dim K[y,z]1\), \(q=\dim K[y,z]x\), and
\(h=p+q-n\), as in the [four-quotient reduction](quadratic_intersection_fans.md).
The feasible integer levels are precisely
\[
                           n-p\le k\le q.                \tag{PS10}
\]
The two endpoints are attained by extending a basis of one layer with
the other; connectedness of the basis graph fills every intervening
integer. Thus there are \(h+1\) levels. Each \(Q_k\) is an integer
polygon in \([0,n(n-1)]^2\). Applying the classical planar lattice
vertex estimate, and counting every positive vertex in its own level,
rewrites the existing slice ceiling as
\[
             g(I)=O\bigl((h+1)(1+n^{4/3})\bigr).         \tag{PS11}
\]
For completeness, an integer polygon in \([0,R]^2\) has total
\(\ell_1\) perimeter at most \(4R\). Its oriented primitive edge
directions are distinct, and only \(O(T^2)\) integer directions have
one-norm at most \(T\). The smallest possible sum of \(v\) such
norms is therefore \(\Omega(v^{3/2})\). Since every actual edge is
a positive integer multiple of its primitive vector, the perimeter
gives \(v=O(1+R^{2/3})\); points and segments satisfy this separately.
Set \(R=n(n-1)\). Every actual initial profile is a vertex of \(Q\),
and hence also of its own integer section, proving (PS11).
This simply combines the existing integer-slice estimate with the
overlap-rank identity. It is no new unrestricted exponent or novelty
claim. Combined with the earlier fixed-\(h\) estimate, it gives
\[
 g(I)=O\!\left(\min\left\{
 n\log(n+1)(1+h)^{10/3},\ (1+h)(1+n^{4/3})\right\}\right).
\]
Consequently any \(g(I)\ge c n^{2+\epsilon}\) family in this class
requires \(h=\Omega_{c,\epsilon}(n^{2/3+\epsilon})\), for large
\(n\). This strengthens the necessary overlap condition recorded
with the earlier estimate, without supplying such a family.

It is **invalid** to apply the planar matroid theorem directly to each
\(Q_k\) on the original labels. Fixing a colour count need not leave
a matroid. More precisely, the two-dimensional faces of a base polytope
are triangles and squares. An integer colour slice can introduce a new
edge only by taking the diagonal of a square whose two independent
exchanges change the colour counts in opposite directions. The new
direction is a sum of two exchanges. A fixed-count basis family is a
matroid exactly when no such square cuts it diagonally: all remaining
edges are single exchanges, and the usual edge characterization applies.

These structural inputs are in
[Feichtner--Sturmfels, Section 2, Definition 2.1 and Proposition 2.4](https://arxiv.org/pdf/math/0411260v1).
Their edge characterization, statement that faces are matroid polytopes,
and dimension formula give the triangle/square fact: a dimension-two
face has either one three-element connected component or two
two-element components; loops and coloops contribute points.
This is a comparison with established matroid geometry, not an assertion
that the paper proves the desired algebraic section bound.

## An actual all-field section needs two simultaneous exchanges

The failure occurs within the elementary socle construction already
used in this repository. For \(A,B\ge2\), let
\[
 I_{A,B}=(x^2,yz,xy^2,xz^2,y^{B+1},z^{A+1},
                             xy-z^A,xz-y^B).             \tag{PS12}
\]
Its monomial carrier has basis
\(D=\{1,x,y,\ldots,y^B,z,\ldots,z^A,xy,xz\}\).
The four tips \(xy,xz,y^B,z^A\) are maximal and annihilated by every
variable. The two independent socle identifications give an actual
length-\(n=A+B+2\) algebra over every field. Put
\(C=D\setminus\{xy,xz,y^B,z^A\}\). Its residue matroid has exactly
four bases: \(C\), one choice from \(\{xy,z^A\}\), and one from
\(\{xz,y^B\}\). All are standard sets. Their switching planes have
normals \((1,1,-A)\), \((1,-B,1)\) and meet at the positive ray
\((AB-1,A+1,B+1)\); hence all four occur as actual initial sets.

Exactly two bases have upper-layer count two:
\[
 S=C\sqcup\{xy,y^B\},\qquad T=C\sqcup\{z^A,xz\}.       \tag{PS13}
\]
Removing \(xy\) from \(S\), replacing it by \(z^A\) changes the
colour count. Replacing it by \(xz\) leaves the dependent pair
\(xz,y^B\). Thus this fixed-count family fails the basis exchange
axiom on the **actual monomial residues**, not just on a hypothetical
matrix. Its integer section is the segment with profile difference
\[
                     \mu(S)-\mu(T)=(0,B+1,-A-1).         \tag{PS14}
\]
For \(A=r,B=r+1\), \(r\ge2\), this is primitive and is not parallel
to any single difference of two nonzero monomial labels of this algebra:
its \(y\)-component \(B+1\) exceeds their entire \(y\)-range \(B\).
This is an infinite exchange obstruction inside genuine quotients.
The full fan has only four cones; it supplies no extremal lower gain.
The socle construction and the section consequence are not assigned a
novelty or publication claim.

There is now an exact counting formulation that retains the coupling.
For \(u>0\), define
\[
 F_k(u)=\min_{v\in Q_k}(u v_y+v_z),\qquad
                F(u,t)=\min_k\{F_k(u)+kt\}.             \tag{PS15}
\]
The integer sequence \(F_k(u)\) is convex, by (PS8) and convexity of
linear-programming section values. Its adjacent thresholds
\(t_k(u)=F_k(u)-F_{k+1}(u)\) are nonincreasing in \(k\).
An intermediate level \(k\) is strictly active exactly on
\(t_k(u)<t<t_{k-1}(u)\), intersected with \(t>0\).
At the extreme levels omit the absent endpoint inequality.
Counting the positive exposed profiles of (PS15), with each profile's
level retained, recovers the actual whole fan. The unresolved task is
to bound these **active constrained polygons collectively** in an
actual algebra. Integrality or single-exchange planar estimates alone
do not solve it; (PS12)--(PS14) give a reproducible guard against that
shortcut.

## Exact controls and the remaining proof obligation

Run from the repository root:

```bash
python3 -B scripts/run_profile_section_checks.py
./build.sh
python3 tests/test_smoke.py
```

The report is `results/research/profile_section_checks.json`. For
\(k=1,\ldots,5\), one implementation builds every zonotope facet by
planar Minkowski hulls and intersects its boundary with the rational
section. A separate integer-ray sweep extracts edges in all eight cost
orthants. Their section vertices agree, with counts
\(60,180,364,612,924\). Explicit strict integer costs check (PS3) in
each orthant. These models are labelled **abstract matroids**, never
polynomial quotients.

Separate actual controls use point evaluations over \(\mathbb F_{101}\)
for the established four-line parity family at lengths \(4,8,12,16\).
An independent full border-fan traversal agrees with the pair-zone sweep
on every profile edge. Coordinate section sweeps integrate to the known
actual length and check (PS6). Their small maximum section counts support
neither an asymptotic fit nor the unproved unrestricted bound.

The binary-section extension exhausts all 127 nonempty restrictions of
the seven nonzero columns of \(\mathbb F_2^3\), with all 2,186 colourings.
Independent intersections of every profile pair agree with the integer
basis hulls and with Minkowski interpolation at fractions
\(1/3,1/2,2/3\). Signed objective checks verify the value-function
convexity. These are bounded controls for the uniform proof, not a
uniform count experiment. For (PS12), parameters \(r=2,3,5,8\) over
characteristics zero, two, three and 101 independently check the whole
four-cone border fans, commuting multiplication, carrier successors,
ordinary S-pairs, all basis ranks and the failed single exchanges.

The state-profile dictionary is the established
[Babson--Onn--Thomas, Theorem 2.4](https://arxiv.org/pdf/math/0207135).
The prior-art comparison above inspects the displayed construction and
the main theorem and section proof of Rörig--Witte--Ziegler. It establishes
that quadratic zonotope cuts are known, and that their stated integer
normalization fails the required colength budget. It does not establish
the strongest current algebraic section theorem or novelty of (PS3)--(PS6).

The substantive task is unchanged: control crossings or their aggregate
for **one actual quotient with changing dependencies**, retaining both
its rank and its exponent/divisor budget. Summing bounds for independent
dyadic box restrictions is still unjustified. A uniform proof must handle
cross-box exchanges and moving matchings in the same algebra. The sharp
\(G_3\) endpoint, an actual polynomial lower separation, and the required
consequential complexity theorem all remain unresolved.
