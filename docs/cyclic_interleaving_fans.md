# Generic cyclic interleavings: exact ranks and a uniform counting bound

The original occurrence/interleaving observations suggest separating monomials
by a residue class and selecting within each stream. This note gives an actual
finite quotient with an exact version of that rule. It is a candidate family
for the original three-variable problem, not a new prime-distribution claim.
Finite-group gradings and their Hilbert schemes are established mathematics.

The main counting conclusion is uniform in both parameters:

\[
 \boxed{g(I_{r,m})=O(n^2\log^2(n+1)),\qquad n=mr.}       \tag{C0}
\]

The implicit constant is absolute, over every base-field characteristic,
with the generic planar columns specified below. A threshold dichotomy
closes the intermediate regime left by the earlier fixed-period and
comparable-period gates. For \(r\ge2\), the number also stabilizes
exactly once \(m\ge r\). An accessible-lattice priority queue selects
the full initial standard set in \(O(n\log(r+1))\) exact operations.

Thus this entire generic cyclic family cannot give a polynomial separation
above the ordinary corner-cut scale. The theorem does not assert a
matching lower bound, apply to arbitrary planar quotients, or settle the
unrestricted \(G_3\) problem. Its precise priority and publication
significance remain unestablished.

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

## Earlier bounds and the formerly open regime

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
them: unbounded period with \(m/r\to0\). This was the remaining
regime before the threshold argument below; (C0) now excludes it as well.

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

These earlier necessary conditions were not evidence that such a family
exists. In particular, \(r=m^2\) lay between the two previous gates.
The following proof replaces this surviving window by the uniform (C0).

## A threshold dichotomy closes every parameter regime

Fix a strict positive exposing weight with a unique cheapest variable
\(h\). The three such sectors cover all open initial cones for an
upper count; their boundaries need not be walls when \(m<r\).
Use the core prices \(\alpha=mw_h>0\) and the two positive gaps
\(\beta,\gamma\) from the earlier sector changes. Let \(\tau_j\)
be the \(r\)-th core cost in color \(j\), after subtracting
\(jw_h\), and put \(t=\tau_0\).

For \(h=x,y\), the feasible cores are

\[
 G_j=\{(k,b,c)\in\mathbb Z_{\ge0}^3:
              c<m,\ b+c\le j+mk\}.
\]

They satisfy \(G_0\subseteq\cdots\subseteq G_{m-1}\) and
\(G_{m-1}+(1,0,0)\subseteq G_0\). For \(h=z\), the feasible
pairs are all \((b,c)\ge0\), with lift
\(k_j=\max(0,\lceil(b+c-j)/m\rceil)\). Their cost decreases
with \(j\), by at most \(\alpha\) in total. In either case,

\[
 \tau_0\ge\tau_1\ge\cdots\ge\tau_{m-1},\qquad
 \tau_0\le\tau_{m-1}+\alpha.                            \tag{C19}
\]

Write \(K,B,C\) for the coordinate maxima over all selected cores.
Every selected core has cost at most \(t\), so
\(\alpha K,\beta B,\gamma C\le t\), and feasibility gives
\(B,C\le m(K+1)-1\).

### Large core threshold: a small axis box

Suppose \(t\ge2\alpha\). In a nonperiodic-cheapest sector,
\(K\ge1\). If \(K\ge4\), all points in the integer box

\[
 \left\lceil K/4\right\rceil\le k\le\lfloor K/2\rfloor,
 \quad 0\le b\le\lfloor B/16\rfloor,
 \quad 0\le c\le\lfloor C/16\rfloor
\]

are feasible in \(G_0\): \(b+c\le m(K+1)/8\le mk\).
Their cost is at most \(5t/8<t\). The number of integer \(k\)
values is at least \((K+1)/8\). If \(1\le K\le3\), use only
\(k=1\); feasibility still holds and its cost is at most
\(\alpha+t/8\le5t/8\). There are fewer than \(r\) points
strictly below the color-zero threshold. In both cases,

\[
 (K+1)(B+1)(C+1)\le2048(r-1).
\]

The original cheapest coordinate is at most \(m(K+1)-1\), while
the other coordinates are bounded by \(B,C\). Hence the original
axis box has volume at most \(2048n\).

For \(h=z\), use the rectangle
\(0\le b\le\lfloor B/16\rfloor\),
\(0\le c\le\lfloor C/16\rfloor\) in color zero. Its lift
obeys

\[
 \alpha k_0+\beta b+\gamma c
 \le t/4+9\alpha/8\le13t/16<t,
\]

since \(k_0\le(B+C)/(16m)+1\le(K+1)/8+1\).
Thus \((B+1)(C+1)\le256(r-1)\); the original periodic
coordinate is below \(m\), so its axis-box volume is at most
\(256n\). We have proved, uniformly in all sectors,

\[
 t\ge2\alpha\quad\Longrightarrow\quad
 \prod_{i=1}^3(1+\max_{a\in D}a_i)\le2048n.             \tag{C20}
\]

### Small core threshold: a thin integer profile

Suppose \(t<2\alpha\). Every selected core has \(k\le1\),
and \(b+c\le j+mk<2m\). There are just two possible core
layers, regardless of the rank. Put
\(\lambda=(\sum k,\sum b,\sum c)\). These are integer
coordinates, related invertibly to the original profile by

\[
 \mu_h=r\sum_{j=0}^{m-1}j+m\lambda_1-\lambda_2-\lambda_3,
 \qquad (\mu_{o_1},\mu_{o_2})=(\lambda_2,\lambda_3).
                                                               \tag{C21}
\]

All these profiles lie in
\([0,n]\times[0,2nm]^2\). The convex hull of any subset of
strictly exposed profiles keeps each of them as a vertex: its strict
exposing functional is still unique on the smaller set. Distinct
initial standard sets cannot share a profile, since that would tie
their costs at every weight. Add a unit cube to handle flat hulls.
Its regularized integer volume is \(O(n^3m^2)\). The classical
three-dimensional lattice-vertex bound therefore gives
\(O(n^{3/2}m)\) such profiles in a sector. If \(m\le r\),
this is \(O(n^2)\).

The cube step preserves the necessary upper bound: its normal fan
refines the original fan, so it has at least as many vertices. No
zero-volume three-dimensional estimate is applied to a flat hull.

### Fixed-volume dyadic boxes count the large-threshold part

For the standard sets in (C20), round each axis side up to a power
of two, then enlarge one side until the product is exactly
\(2^Q\), where \(Q=\lceil\log_2(16384n)\rceil\).
There are only \(\binom{Q+2}{2}=O(\log^2(n+1))\) possible
triples of nonnegative side exponents with that sum. The common
volume satisfies \(2^Q<32768n\).

For a box with sides \(L_1,L_2,L_3\), profiles of its length-\(n\)
standard sets lie in a coordinate box of volume at most
\(n^3L_1L_2L_3=O(n^4)\). Convex hulls of subsets again retain
every strictly exposed profile; cube regularization costs only a
constant factor. The lattice-vertex bound gives \(O(n^2)\)
profiles per box. Together with (C21) and the three sectors this proves
(C0) for \(m\le r\). If \(m>r\ge2\), use (C16) with period
\(r\), whose length \(r^2\) is at most \(n\). The cases
\(r=1\) and \(m=1\) are covered by the explicit three-state
description and the planar bound. This completes the uniform proof.

The geometric input is the established
[Andrews lattice-vertex theorem](https://doi.org/10.1090/S0002-9947-1963-0143105-7),
not a new vertex exponent. An accessible proof retelling is
[Dillon, Theorem 1 and Sections 1--2.2](https://arxiv.org/pdf/2606.30856).
The threshold split and the integer coordinates in (C21) are what
permit its application here. The argument excludes polynomial
amplification in this generic family; it does not improve unrestricted
\(G_3\) or prove that the upper bound is attained.

## Direct selection without a universal monomial ground

The same cores give an exact full-basis selection algorithm. If \(h=x,y\),
start a min-priority queue at \((0,0,0)\) in each color and insert
unseen feasible neighbors obtained by incrementing one core coordinate.
Every feasible point has a path from the root: increase \(k\) first,
then \(b,c\). Every edge strictly increases cost because
\(\alpha,\beta,\gamma>0\). Therefore the next queue minimum
is the next cheapest feasible point. With \(h=z\), use the analogous
two-dimensional queue on \((b,c)\), with its unique lift; its
cost is strictly increasing in either coordinate.

Pop exactly \(r\) times per color. Original-exponent lexicographic
tie breaking agrees with the direct finite-ground oracle, even at
within-sector cost ties. A color uses at most \(3r-2\) inserted
points and at most \(2r-1\) queued points. The resulting cost is

\[
 O(n\log(r+1))\text{ exact arithmetic/comparison operations},
 \quad O(r)\text{ working queue/seen storage},              \tag{C22}
\]

plus \(O(n)\) storage if the whole output is retained. This is
best-first enumeration on an accessible integer ground, using a
standard priority queue. A balanced dictionary for seen membership
gives the deterministic comparison bound; the checker uses Python's
hashed tuple sets. It is not a new general enumeration method,
a bit-complexity bound, or random access to a smooth-number rank.
It requires the generic uniform-rank rule (C4).

## Why a single small axis box is not a uniform theorem

The threshold split is necessary for this proof. For every integer
\(s\ge3\), put \(m=s^3\), \(r=s^4\), \(q=s^2\),
\(R=2m\), \(L=8m^3\), and take weights
\((L,L+1,L+R)\). In each color the relevant monomials are
strictly ordered first by total degree, then by \(c\), then by
\(b\). At most the tiers \(j,j+m\) are needed: the degree-\(m\)
tier already has more than \(r\) points. Within these tiers
\(b<2m=R\), so all boundary comparisons are strict.

The standard set contains \(x^m\) (second in color zero),
\(y^{m-1}\) (rank \(m\) in color \(m-1\)), and \(z^q\)
(its whole degree-\(q\) tier has
\((q+1)(q+2)/2<r\) points). Higher tiers cannot precede these
witnesses, since \(mL>qR\). Consequently

\[
 \frac{\text{axis-box volume}}{n}
 \ge\frac{(m+1)m(q+1)}{mr}>s.                            \tag{C23}
\]

Here both \(m\) and \(r/m=s\) grow. Thus an unconditional
\(O(n)\) box-volume assertion is false even in the former
intermediate window. These are bases of the actual generic quotient,
not abstract downsets. They are precisely in the small-threshold
branch of (C21); they do not contradict (C0).

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
The specialized counting charge, stabilization and threshold-dichotomy
application still need exact priority assessment. The inspected grading
passages do not establish that assessment. A keyword search finding no
equivalent theorem would not establish novelty either. The Andrews
statement was rechecked in Dillon's retelling; that source explicitly
attributes the theorem and proof strategy to earlier work. No new
lattice-vertex theorem or general priority-queue technique is asserted.

```bash
python3 -B scripts/run_cyclic_interleaving_checks.py
./build.sh
python3 tests/test_smoke.py
```

The report retains exact complete fans, rank and commuting-multiplication
controls in characteristic two as well as other fields, sector compression,
the fixed-direction facet charge and stabilization witnesses. It also
compares the new heap oracle with every saved complete-fan basis and
independent finite-ground sorting; explicitly counts the strictly cheaper
integer rectangles in (C20); checks both regularized profile-volume
branches; and independently replays three thin-threshold bases by degree
tiers. The latter includes the growing-period, growing-rank-ratio regime.
Computation supports the uniform proofs; finite counts are not extrapolated
into an extremal claim. Core algorithms, manuscripts and the paused ML
work are preserved. The main goal remains active.
