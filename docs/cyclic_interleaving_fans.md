# Generic cyclic interleavings: an exact model and two counting gates

The original occurrence/interleaving observations suggest separating monomials
by a residue class and selecting within each stream. This note gives an actual
finite quotient with an exact version of that rule. It is a candidate family
for the original three-variable problem, not a new prime-distribution claim.
Finite-group gradings and their Hilbert schemes are established mathematics.

The main structural conclusions are uniform proofs, not growth fits:

* Every fixed period has \(O_m(n^2\log(n+1))\) initial ideals.
* If the period is at least a fixed positive multiple of the base rank,
  the bound is \(O_c(n^2\log^8(n+1))\).
* For base rank \(r\ge2\), the number of initial ideals is exactly
  independent of the period once \(m\ge r\).

Thus a polynomial separation in this family would require both
\(m\to\infty\) and \(r/m\to\infty\). That intermediate regime remains
open. The unrestricted \(G_3\) bound, novelty assessment and publication
goal are unchanged.

## An actual quotient with uniform ranks in each stream

Fix a base field \(F\), extend it to
\(K=F(U_1,V_1,\ldots,U_r,V_r)\), and let \(B=K^r\). Put
\(U=(U_i)\), \(V=(V_i)\). For \(m\ge1\), define

\[
 A_{r,m}=B\otimes_K K[g]/(g^m-1),\qquad
 x=U g,\quad y=V g,\quad z=g.                           \tag{C1}
\]

These three elements generate the algebra: \(z^{-1}=z^{m-1}\),
and \(xz^{-1}=U\) generates \(B\) by interpolation at the distinct
\(U_i\). Its length is \(n=mr\). More explicitly, if
\(f(T)=\prod_i(T-U_i)\) and \(q(U_i)=V_i\), \(\deg q<r\),
its kernel is

\[
 I_{r,m}=(z^m-1,\ f(xz^{m-1}),\ yz^{m-1}-q(xz^{m-1})).    \tag{C2}
\]

This works in every characteristic, including when it divides \(m\).
No diagonalization or distinct-root assumption on \(g^m-1\) is used.
The function-field extension specifies generic coefficients; it is not
an assertion that a small finite field contains suitable generic points.
The finite conditions determining a whole fan can also be specialized
over suitable extensions. Uniform ranks on *all* monomials are asserted
for the function-field model, not for every finite specialization.

After reducing \(0\le c<m\), a monomial has column

\[
 x^ay^bz^c\longmapsto (U_i^aV_i^b)_{i=1}^r
                 \otimes g^{a+b+c\bmod m}.              \tag{C3}
\]

Any \(r\) distinct bivariate monomials give independent columns:
the evaluation determinant is a nonzero polynomial in the independent
\(U_i,V_i\). Its permutation terms have distinct exponent tuples and
coefficients \(\pm1\), also in characteristic two. Consequently a
finite monomial set \(S\) has rank

\[
 \operatorname{rank}(S)=\sum_{j=0}^{m-1}
 \min\bigl(r,\#\{(a,b):(a,b,c)\in S,\ a+b+c\equiv j\pmod m\}\bigr).
                                                               \tag{C4}
\]

Duplicates with the same \((a,b,j)\) are counted only once.
For every strict positive weight, the initial standard set consists
exactly of the \(r\) cheapest distinct columns in each residue class.
Since \(z^m-1\in I\), only \(c<m\) is needed. The usual universal
ground \(\prod(a_i+1)\le n\) contains every resulting standard set.
This is a commuting cyclic quotient, not merely a chosen abstract matroid.

There is also an elementary converse with a precise scope. Any finite
\(\mathbb Z/m\)-graded quotient, with all three variables of degree one
and \(z^m=1\), has this tensor form for an **arbitrary** planar algebra
\(B=A_0\). Indeed \(xz^{-1},yz^{-1}\) generate \(A_0\), and
\(A_j=A_0z^j\). For such a \(B\), replace each minimum in (C4) by
the actual planar column rank. The uniform minimum rule and the subsequent
specialized counting gates require generic \(B\); they do not apply to
every planar algebra.

## Fixed periods: charge whole facets along one direction

Assume \(m\ge2\). At an interior fan ray write the primitive normal as
\(w=(d p,d q,d s)\), where \(p,q\) are positive coprime integers and
\(s>0\) is rational. Thus \(d=\gcd(w_x,w_y)\).
For costs \((p,q,s)\), let \(t_j\) be the \(r\)-th cost in class
\(j\) and \(T=\max_jt_j\). The \(r\) axis candidates
\((j+km,0,0)\), \(0\le k<r\), and their \(y\)-axis counterparts
give

\[
 T\le(n-1)\min(p,q).                                    \tag{C5}
\]

The exposed two-dimensional profile face has two independent exchange
directions \(\delta_1,\delta_2\). Each exchanges tied monomials in
one class. Their endpoints have nonnegative planar costs at most \(T\),
and \(|\delta_{i,z}|\le m-1\). Eliminating \(s\) gives a nonzero
integer vector perpendicular to \((p,q)\), so

\[
 pq\le2(m-1)T.                                         \tag{C6}
\]

The eliminated vector cannot vanish: if both vertical differences vanish,
the two directions would already be parallel, and otherwise its vanishing
would make them proportional. These are the two excluded alternatives.

Here is a useful lattice count. If \(T\ge8m\min(p,q)\), the rectangle
\(0\le a<T/(2p)\), \(0\le b<T/(2q)\) contains at least
\(T^2/(8mpq)\) points of any class \(a+b\bmod m\).
Count residues along the cheaper coordinate, whose interval has length
at least \(4m\); the other coordinate has at least its interval length
many nonnegative integer choices. All these points have planar cost
strictly below \(T\). For a class attaining \(T\), there are fewer
than \(r\) such points, whence

\[
 T^2\le8npq,\qquad T\le16mn,\qquad pq\le32m^2n.         \tag{C7}
\]

If \(T<8m\min(p,q)\), (C6) gives
\(\max(p,q)<16m^2\) and \(pq<256m^4\). Therefore every interior
ray has

\[
 pq\le D:=32m^2n+256m^4.                                \tag{C8}
\]

There are at most \(D(1+\log D)\) possible positive pairs. Counting
all three ray coordinates separately would waste the structure.

Fix one pair \((p,q)\) and vary \(s>0\). At a two-dimensional
exposed profile face let \(\ell\) be its vertical coordinate range.
It is a positive integer. In the face plane, points of one fixed vertical
coordinate lie on a line, so at most two can be vertices. Thus its number
\(k\) of vertices satisfies

\[
 k\le2(\ell+1)\le4\ell.                                \tag{C9}
\]

As \(s\) increases, the vertical slope of the minimum-cost function is
nonincreasing. The vertical jump intervals of its exposed faces have
disjoint interiors and total length at most \(n(m-1)\). Hence the sum
of the degrees of **all** interior rays having this pair is at most
\(4n(m-1)\). This charges an entire family of facets together.
It is not an individual-cone margin estimate.

Summing over (C8), and adding the existing universal
[boundary-cell bound](weighted_facet_incidence.md), gives

\[
 g(I_{r,m})=O\left(n^{5/3}\log(n+1)
           +n(m-1)D(1+\log D)\right)
          =O_m(n^2\log(n+1)).                            \tag{C10}
\]

Cells avoiding the simplex boundary have an interior ray and are counted
by those incidences. Boundary-adjacent cells are covered by the cited
bound, including closures touching the boundary only at a point.
Lower-dimensional profile polytopes cause no problem: regularize by a
cube; every strictly positive supporting face is merely translated.
For \(m=1\), the known planar corner-cut bound suffices directly.

The same proof bounds normal heights by
\(32m^3n+128m^4\). Indeed one exchange has a nonzero vertical
difference of magnitude at most \(m-1\), so \(d\le m-1\) and
\(s\le T\). Thus the square-root weighted interior incidence is
\(O_m(n^{5/2}\log(n+1))\). This is a family result, not the
unproved unrestricted cubic incidence budget discussed below.

## A shared finite ground when the period is comparable to the rank

There is a second, independent counting gate. Suppose first that \(x\)
is uniquely cheapest. Write

\[
 \alpha=mX,\qquad\beta=Y-X>0,\qquad\gamma=Z-X>0.
\]

In class \(j\), put \(a+b+c=j+mk\). A candidate has original
exponents \((j+mk-b-c,b,c)\) and, up to the common cost \(jX\),
cost \(\alpha k+\beta b+\gamma c\).
Every selected candidate satisfies \(b+c\le r-1\): along a
coordinate path reducing \((b,c)\) while keeping the total degree
fixed, a candidate of greater degree has at least \(r\) strictly
cheaper independent predecessors. This path respects \(c<m\).

Put \(K_0=\lceil(r-1)/m\rceil\). For fixed \(k<K_0\), all
coordinate predecessors in \((b,c)\) are feasible, so a selected
point has \((b+1)(c+1)\le r\). For \(k\ge K_0\), all levels
\(K_0,\ldots,k\) with those predecessors are feasible, giving

\[
 (k-K_0+1)(b+1)(c+1)\le r.                              \tag{C11}
\]

Thus all selected columns, for all colors and all weights in this
sector, belong to one ground of size

\[
 O\bigl(r\log(r+1)(K_0+1+\log(r+1))\bigr).             \tag{C12}
\]

This uses the familiar two- and three-dimensional hyperbolic crosses,
with a bounded number of shifted boundary layers. The \(y\)-cheapest
sector is symmetric.

If \(z\) is uniquely cheapest, set
\(\alpha=mZ\), \(\beta=X-Z>0\), \(\gamma=Y-Z>0\).
The unique allowed lift of \((a,b)\) in class \(j\) has cost

\[
 jZ+\beta a+\gamma b+
 \alpha\left\lceil\frac{a+b-j}{m}\right\rceil.           \tag{C13}
\]

The ceiling is nonnegative. This is increasing in each coordinate, so
selected pairs form a downset of size \(r\). Therefore
\((a+1)(b+1)\le r\), \(a+b\le r-1\), and its lift index is
between zero and \(K_0\). This is another common ground obeying
(C12).

All selections are fixed by the ordering of their common ground.
For a ground of size \(M\), at most \(M^2\) comparison planes give
\(O(M^4)\) projective two-parameter cells. This is the classical
arrangement bound, with no factor for the number of colors. It follows
from (C12) that for every fixed \(c>0\),

\[
 m\ge cr\quad\Longrightarrow\quad
 g(I_{r,m})=O_c(r^4\log^8(r+1))
           =O_c(n^2\log^8(n+1)).                         \tag{C14}
\]

The three sectors cover every open initial cone; an upper bound does
not require that their artificial boundaries be actual walls.

## Exact stabilization for \(m\ge r\)

For \(r\ge2\) and \(m\ge r\), there is more than an upper bound.
The subdivisions inside all three cheapest-variable sectors are independent
of \(m\) after the preceding invertible linear cost changes.

In the \(x\)-cheapest sector, all selected \((k,b,c)\) have
\(k+b+c\le r-1\): first decrease \(b,c\) to zero, then decrease
\(k\). When \(m\ge r\), these are distinct feasible columns along
a decreasing path. Within this finite simplex, the color-\(j\) ground
is simply

\[
 k\ge1\quad\hbox{or}\quad b+c\le j.                     \tag{C15}
\]

All \(j\ge r-1\) have the same unrestricted ground. The only
exceptional grounds are \(j=0,\ldots,r-2\). They do not depend on
\(m\). The multiplicity of the unrestricted ground changes, but repeating
the same selection fan does not refine it. The \(y\) sector is identical.

In the \(z\)-cheapest sector, (C13) reduces to
\(jZ+\beta a+\gamma b+\alpha\mathbf1_{a+b>j}\), with
\(a+b\le r-1\). Again the only exceptional colors are below \(r-1\),
and the same independence follows.

These sector counts really add. If two variables tie for cheapest, the
color-zero constant is cheapest and the next tier contains at least
\(m\) distinct columns of degree \(m\) using those two variables.
Only \(r-1<m\) of them are selected. The two perturbations of the
cheapest tie therefore give different initial sets, independently of
the third, more expensive price. Thus every such tie is an actual wall;
a maximal cone cannot cross from one cheapest-variable sector to another.
For two nonperiodic variables the tier has \(m+1\) columns instead,
which only strengthens the argument.

Consequently

\[
 \boxed{g(I_{r,m})=g(I_{r,r})\quad(r\ge2,\ m\ge r).}     \tag{C16}
\]

The map between sector weights preserves their gaps and sends their
common minimum \(t\) to \((m/m')t\) when comparing periods
\(m,m'\ge r\). It is linear within each sector. The equality concerns
numbers of initial ideals, not equality of standard monomial labels
or quotient lengths. For \(r=1,m\ge2\), there are exactly three
initial sets; \(r=m=1\) has one.

## A coarse uniform bound and the remaining regime

All initial profiles satisfy \(\mu_x,\mu_y\le n(n-1)/2\),
\(0\le\mu_z\le n(m-1)\), and

\[
 \mu_x+\mu_y+\mu_z\equiv r\sum_{j=0}^{m-1}j\pmod m.
\]

The profile differences therefore lie in the index-\(m\) lattice
\(L_m=\{v:\sum_i v_i\equiv0\pmod m\}\). Add a fundamental
parallelepiped with basis \((1,-1,0),(0,1,-1),(0,0,m)\).
It makes a full-dimensional lattice polytope in the same lattice coset,
refines its normal fan, and has normalized volume \(O(n^5)\).
The classical three-dimensional lattice-vertex bound gives

\[
 g(I_{r,m})=O(n^{5/2})                                  \tag{C17}
\]

uniformly in \(m,r\). This is a normalized-volume corollary, not a
new lattice-vertex theorem. Its proof in fact applies to every planar
\(B\) in the preceding converse: multiplication by \(z\) identifies
all graded pieces, so each has dimension \(r\), and the same profile
congruence and coordinate bounds hold. The stronger (C10), (C14) and
stabilization (C16) remain statements about generic \(B\).
Bounds (C10) and (C14) are stronger in their
stated regimes. A possible polynomial amplification must lie between
them: unbounded period with \(m/r\to0\). Exact ranks make that a
defined counting problem; no amplification is currently proved.

The explicit costs sharpen this gate. Put \(L=1+\log(n+1)\). For
\(2\le m\le r\), (C8)--(C12) give, with an absolute constant,

\[
 g(I_{r,m})\le C\min\{n^{5/2},\ m^3n^2L,\ n^8m^{-12}L^8\}.
                                                               \tag{C17a}
\]

For the last term, \(K_0+1\ll r/m\) and (C12) bounds the common
ground by \(O((r^2/m)\log^2(r+1))\); use its fourth-power
arrangement bound and \(r=n/m\). For the middle term use
\(m^2\le n\) in (C8). Consequently a putative family with
\(g\ge c n^{2+\varepsilon}\), for fixed \(c,\varepsilon>0\),
would need

\[
 m=\Omega_{c,\varepsilon}(n^{\varepsilon/3}L^{-1/3}),\qquad
 m=O_{c,\varepsilon}(n^{1/2-\varepsilon/12}L^{2/3}).       \tag{C17b}
\]

These are necessary conditions, not evidence that such a family exists.
For example \(r=m^2\) lies between the two excluded regimes and
retains only exact residue-count ranks. It is a defined next test case,
not an established lower bound or a reason to expand the search framework.

## Why edge length was investigated

Let \(T_{1/2}(I)=\sum_w d(w)\sqrt{\|w\|_\infty}\), over positive
interior rays. For a full-dimensional integer profile polytope, each
facet vertex with adjacent edge vectors \(a,b\) satisfies
\(\sqrt{\|w\|_\infty}\le\sqrt{\|a\|\|b\|}
\le(\|a\|+\|b\|)/2\). Summing shows

\[
 T_{1/2}(I)\le2L_+(I),                                  \tag{C18}
\]

where \(L_+\) sums the lengths of edges incident to a strictly positive
supporting facet; cube regularization handles flat profiles. If the
*unproved* unrestricted budget \(T_{1/2}=O(n^3\log^\sigma n)\) held,
the existing small-normal shadow count would give
\(G_3(n)=O(n^{33/13}\log^{(3+10\sigma)/13}(n+1))\), plus its
smaller boundary term. The three-pole quotient has
\(T_{1/2}=\Theta(n^3)\) along the coprime subsequence from the
[quartic-incidence proof](weighted_facet_incidence.md), so the candidate
power would be sharp for that budget. Neither (C18), those examples,
nor the fixed-period result proves the unrestricted budget. No such
improvement to \(G_3\) is claimed.

## Prior work, verification and limits

[Haiman--Sturmfels, arXiv:math/0201271v1](https://arxiv.org/pdf/math/0201271),
the introduction's grading/Hilbert-function definitions and Section 5,
Remark 5.1, provide the established finite-abelian grading framework.
Their toric/orbit comparison has Hilbert function one in each degree;
here it is \(r\). The construction is not claimed as a new kind of
Hilbert scheme. Only those passages were inspected in this comparison;
no full-paper or strongest-current novelty certificate is asserted.
The matroid/initial dictionary and lattice-vertex tools are those already
used in [the general research note](extremal_initial_ideals.md).
The specialized counting charge and stabilization still need exact
priority assessment.

```bash
python3 -B scripts/run_cyclic_interleaving_checks.py
./build.sh
python3 tests/test_smoke.py
```

The report retains exact complete fans, rank and commuting-multiplication
controls in characteristic two as well as other fields, sector compression,
the fixed-direction facet charge and stabilization witnesses. Computation
supports the uniform proofs; finite counts are not extrapolated into
an extremal claim. Core algorithms, manuscripts and the paused ML work
are preserved. The main goal remains active.
