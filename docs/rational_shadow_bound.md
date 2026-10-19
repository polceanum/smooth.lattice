# Rational shadows of hyperbolic matroid profiles

This note extends the planar exchange argument to rational projections of
integer exponent vectors in any fixed dimension. In three variables it
gives the uniform all-fields bound

$$
G_3(n)=O\!\left(n^{23/8}\log^{3/8}(n+1)\right). \tag{1}
$$

This is an internally proved upper bound, not a sharp asymptotic answer,
a verified novelty claim, or completion of the active goal. The generic
scale remains separated from (1). The arithmetic projection argument is
given explicitly; it does not depend on the earlier continuous-shadow
inequality or on the uninspected part of McMullen's 1975 application.
Its own precise priority comparison remains necessary.

## A height-dependent shadow theorem

Let a rank-$r$ matroid have integer weight vectors $a(e)\in\mathbb N_0^d$
satisfying $\prod_i(a_i(e)+1)\le H$, where $d\ge2$ is fixed and $H\ge1$.
Let $A$ be a rank-two integer $2\times d$ matrix. Write

$$
\delta=\gcd_{i<j}|\det(A_i,A_j)|,\qquad
\Delta(A)=\delta^{-1}\sum_{i<j}|\det(A_i,A_j)|.
$$

Here $A_i$ denotes its $i$th column, zero minors are allowed, and
$\delta>0$. The lattice $\Lambda=A\mathbb Z^d$ has covolume $\delta$;
this is the usual rank-two integer-lattice identity, for example by
Smith normal form. The normalized height $\Delta$ is unchanged by an
invertible rational change of the two image coordinates, whenever the
resulting matrix is integral.

**Theorem 1.** The polygon of projected basis sums has at most

$$
C_d(1+\log(H+1))
  \left[1+r+(r^2H\Delta(A))^{1/3}\right] \tag{2}
$$

vertices. Repeated projected weights and arbitrary matroids are allowed.
We count exposed profiles, not every basis tied at a degenerate objective.
In dimension three, if $q$ is the primitive integer kernel vector of $A$,
then $\Delta(A)=\|q\|_1$. Thus (2) applies to the orthogonal shadow on
$q^\perp$, up to an invertible change of image coordinates.

### A cover by small lattice parallelograms

Put $L=\lfloor\log_2 H\rfloor$. For every nonnegative tuple
$(j_1,\ldots,j_{d-1})$ of sum at most $L$, use an anchored box with sides

$$
b_i=2^{j_i+1}-1\ (i<d),\qquad
b_d=\max(1,\lfloor H/2^{\sum j_i}\rfloor-1).
$$

These $K=\binom{L+d-1}{d-1}$ boxes cover the weight vectors: choose
$j_i=\lfloor\log_2(a_i+1)\rfloor$. All $b_i\ge1$, and

$$
\prod_i(b_i+1)\le2^dH,\qquad b_i b_j\le2^dH. \tag{3}
$$

For one box put $g_i=b_iA_i$, and choose $g_s,g_t$ maximizing
$D=|\det(g_s,g_t)|>0$. Orient this pair so that
$m=\det(A_s,A_t)>0$, and $D=b_s b_t m$.
Cramer's rule writes every $g_i$ in this pair with both coefficients
of absolute value at most one. Define the integer map

$$
T(v)=(\det(v,A_t),\det(A_s,v)),\qquad\det T=m.
$$

The projected box lies in the interior of the parallelogram

$$
P=\{v: |(Tv)_1|<db_s m+\tfrac12,
        |(Tv)_2|<db_t m+\tfrac12\}.
$$

Its four boundary lines contain no point of $\Lambda$, since their
levels are half-integers. Write $W_s=2db_s m+1$ and
$W_t=2db_t m+1$ for its two coordinate widths. We have

$$
\frac{W_sW_t}{\delta m}\le(2d+1)^2\frac D\delta
 \le(2d+1)^2 2^dH\Delta(A). \tag{4}
$$

Order these $K$ parallelograms once and label every projected ground
point by the first one containing it. The cover may overlap; uniqueness
of the label, rather than disjointness of the cover, is all that is used.

### Convex token paths and boundary crossings

Sweep costs $tX+Y$, $t>0$, on the projected ground. The adjacent greedy
exchange proof in the [planar note](planar_dyadic_bound.md) provides
$r$ token paths. Every successful exchange has displacement
$(-u,tu)$ with $u>0$, and its parameter $t$ never decreases along a
token. Simultaneous crossings give consecutive collinear steps.

For any linear functional $\ell(X,Y)=\alpha X+\beta Y$, the sign of
its change along such a step is the sign of $-\alpha+\beta t$.
Consequently $\ell$ is monotone in at most two consecutive phases
along each token. This holds for every linear functional, not just
the two cost coordinates.

A token therefore crosses any line avoiding the ground points at most
twice. An exchange between distinct labels crosses a boundary line of
one of the cover parallelograms. There are only $O_d(L+1)$ distinct
such lines across the entire cover. Indeed each is of the form

$$
\det(v,A_j)=\pm(db_i|\det(A_i,A_j)|+\tfrac12).
$$

For each ordered column pair, $b_i$ takes at most $L+1$ values;
this also holds for the last coordinate, whose side depends only
on $\sum j_i$. Thus $O_d(r(L+1))$ exchanges change labels, even
though the cover itself can have many more boxes.

For exchanges with both endpoint labels equal to a given box, each
$T$ coordinate has total absolute change at most twice its width per token.
Indeed in each monotone phase, the counted endpoint intervals are
disjoint subintervals of the corresponding width interval. Leaving a label
and later returning does not invalidate this budget. Summing over
tokens bounds the two absolute-coordinate budgets by $2rW_s,2rW_t$.

### Integer directions in a parallelogram

We use a lattice version of the elementary direction estimate. Suppose
$k$ vectors in a planar lattice of covolume $c$ have different ray
directions in a pointed cone, and their absolute coordinate sums are
at most $U,V$. Then

$$
k\le6+(432UV/c)^{1/3}. \tag{5}
$$

If $k\ge2$, both budgets are positive; otherwise all vectors would lie
on one line. At least $k/3$ vectors have
$|x|/U+|y|/V\le3/k$. They lie in a centered rectangle of area
$36UV/k^2$. Order them by angle in their pointed cone. The triangles
from the origin between successive vectors have disjoint interiors,
are contained in the rectangle, and each has area at least $c/2$.
Thus $k/3\le1+72UV/(ck^2)$. For $k\ge6$ this gives
$k^3\le432UV/c$; smaller $k$ are covered by (5). No density estimate
for a thin oblique lattice rectangle is being assumed.

For each parameter having exchanges internal to a label, sum their
displacements. The resulting nonzero vector belongs to $\Lambda$
and has direction $(-1,t)$. Distinct parameters give distinct rays
in a pointed cone; their images under $T$ still do so. Their absolute
coordinate sums are bounded by the exchange budgets
$U\le2rW_s,V\le2rW_t$. The transformed lattice has covolume
$\delta m$. Equations (4)--(5), using the sharper first inequality
in (4), give at most $6+C_d(r^2D/\delta)^{1/3}$ parameters per
label. If $r\ge1$, the constant six is absorbed in this bound,
since $D/\delta$ is a positive integer.

### Summing the box budgets

The last step avoids paying the largest budget for every box. For
every fixed coordinate pair $i,j$, the dyadic sides satisfy

$$
\sum_{\text{boxes}}(b_i b_j)^{1/3}
 \le C_d H^{1/3}(L+1). \tag{5a}
$$

For $i,j<d$, put $s=j_i+j_j$ and $t=L-s$. There are $s+1$
ways to choose these two indices, and $O_d((t+1)^{d-3})$ ways
to choose the others with sum at most $t$. Also
$b_i b_j\le4\,2^s$. Summing gives an upper bound proportional to

$$
2^{L/3}(L+1)\sum_{t\ge0}2^{-t/3}(t+1)^{d-3},
$$

which proves (5a). This case only occurs when $d\ge3$.
For $j=d$, put $t=\sum_{k<d,k\ne i}j_k$. Then
$b_i b_d\le2H\,2^{-t}$, there are at most $L-t+1$ choices
for $j_i$, and $O_d((t+1)^{d-3})$ choices for the other indices.
The same convergent sum proves the assertion when $d\ge3$.
When $d=2$ there are no other indices and
$b_1b_2\le2H$ for all $L+1$ boxes, which proves it directly.
For an explicit constant, $C_d=2\cdot5^{d-2}$ suffices in (5a):
use $2^{-1/3}<4/5$ and
$\sum_{t\ge0}\binom{t+d-3}{d-3}(4/5)^t=5^{d-2}$
in the two cases for $d\ge3$.

As $D=\max_{i<j}b_i b_j|\det(A_i,A_j)|$, (5a) implies

$$
\sum_{\text{boxes}}(D/\delta)^{1/3}
 \le C_d(H\Delta)^{1/3}(L+1). \tag{5b}
$$

Here we bound the maximum by the sum and use the elementary
inequality $\sum c_i^{1/3}\le k^{2/3}(\sum c_i)^{1/3}$ for
the fixed number $k=\binom d2$ of minors. Every changing parameter
has a label-changing or an internal exchange. Combine (5b) with
the shared boundary-line budget and the initial state to obtain
(2) in one quadrant. Rank zero has just one state separately.

The other three quadrants have displacements $( -\sigma u,\tau tu)$,
where $\sigma,\tau\in\{-1,1\}$. Every linear functional still changes
sign at most once, and the same proof applies. Four quadrant counts
cover every polygon vertex, including for a segment or singleton.
This proves Theorem 1 without a generic-position assumption.

## Small integer projections see all vertices with short normals

**Lemma 2.** Let $P$ be a full-dimensional lattice polytope in
$\mathbb R^3$, of surface area $S$. For any real $h\ge1$, call a
vertex bad if an incident facet has primitive integer outward normal
of infinity norm greater than $h$. There are at most $6S/h$ bad
vertices. Every other vertex appears as a vertex in some orthogonal
shadow on $q^\perp$, with primitive integer $q\ne0$ and
$\|q\|_\infty\le\lceil3\sqrt h\rceil$.

For a facet with $k$ vertices and primitive normal $z$, a fan
triangulation has $k-2$ integer triangles. Each has area at least
$\|z\|_2/2$, because its integer cross product is a nonzero integer
multiple of $z$. As $k\ge3$, this gives
$k\le6\operatorname{area}(F)/\|z\|_2$. Summing over long-normal
facets bounds the bad vertices by $6S/h$. This facet charging is
classical lattice-polytope geometry, already used in the older note.

At a good vertex choose three independent incident facet normals
$z_1,z_2,z_3$ and put $w=z_1+z_2+z_3$. Their positive sum lies
in the interior of their full-dimensional simplicial cone, hence
strictly exposes the vertex. Also $w\ne0$ and $\|w\|_1\le9h$.
Put $Q=\lceil3\sqrt h\rceil$. Among the $(Q+1)^3$ integer points
of $\{0,\ldots,Q\}^3$, the integer functional $a\mapsto w\cdot a$
takes at most $Q\|w\|_1+1$ values. Since $Q^2\ge\|w\|_1$,
the former number is strictly larger. Two points have the same
value; their difference gives nonzero integer $q$ with
$w\cdot q=0$ and $\|q\|_\infty\le Q$. Dividing by its gcd makes
it primitive. The objective $w$ belongs to $q^\perp$ and still
strictly exposes the vertex, so it exposes its projected vertex.
This is the elementary short-relation pigeonhole argument, not a
new geometry-of-numbers theorem.

Writing $s_q(P)$ for the number of vertices of this shadow, we obtain
the useful general incidence bound

$$
f_0(P)\le6S/h+
 \sum_{\substack{q\text{ primitive, up to sign}\\
                 0<\|q\|_\infty\le\lceil3\sqrt h\rceil}}s_q(P).
\tag{6}
$$

The sum is an upper bound even when one vertex is visible many times.
No averaging over irrational projections is required.

The [five-cone quotient construction](initial_cone_height.md) shows that
the square-root radius cannot be improved uniformly, even for actual
finite quotients over every field. Its exact minimum radius is
$\lfloor k/2\rfloor+1$ while its largest primitive cone ray has height
$k(k+1)$. Yet two planes cover all five cones. Individual radius and
the number of planes in a selected cover therefore require separate
estimates. This does not assert that summing all directions in (6) is
optimal for quotient state polytopes.

## Application to arbitrary three-variable quotients

Let $I\subset K[x,y,z]$ have colength $n\ge1$. Use its rank-$n$
residue matroid on the established universal ground
$V_n^3=\{a:\prod_i(a_i+1)\le n\}$. Let $B$ be the convex hull of
all its basis sums. The universal-ground/greedy/state dictionary of
[Babson--Onn--Thomas, Proposition 2.2 and Theorem 2.4](https://arxiv.org/pdf/math/0207135)
injects initial standard sets into vertices of $B$. The reason for
injectivity is that strict cost orders give a unique greedy basis;
two distinct bases with the same sum cannot both be strict minima.

Every coordinate of $B$ lies in $[0,n(n-1)]$. Regularize by
$P=B+[0,1]^3$, preserving at least as many vertices and making a
full-dimensional lattice polytope. It lies in a box of side
$R=n(n-1)+1$, and $S(P)\le6R^2=O(n^4)$. To see this elementary
surface bound, weight each facet area by the sum of absolute
components of its unit normal. This is at least its area, while
each coordinate projection is counted twice and has area at most
$R^2$. A projected cube has at most six vertices; Minkowski addition
therefore increases each shadow count by at most six.

For primitive $q$ of height at most $Q$, choose any two independent
integer rows orthogonal to $q$ to form $A$. Their cross product is
an integer multiple of $q$, so $\Delta(A)=\|q\|_1\le3Q$.
Theorem 1, with $r=H=n$, gives

$$
s_q(P)=O\!\left(n Q^{1/3}\log(n+1)+1\right). \tag{7}
$$

There are $O(Q^3)$ possible $q$. Lemma 2 and (7) imply

$$
g(I)\le f_0(P)
 =O\!\left(nh^{5/3}\log(n+1)+n^4/h\right). \tag{8}
$$

Taking $h=\max(1,n^{9/8}/\log^{3/8}(n+1))$ proves (1).
All constants are absolute, with no restriction on the characteristic,
radicality or presentation of $I$.

## Algebraic consequence and scope

An $n$-element downset in three coordinates has $O(n^{2/3})$ minimal
outside monomials. For completeness put $t=\lceil n^{1/3}\rceil$.
Those in $[0,t-1]^3$ number at most $t^2$, by antichain projection.
For each coordinate, an outside minimal exponent at least $t$ in
that coordinate supplies $t$ inside points along its predecessor
chain; the chains are disjoint by the antichain property. There are
at most $3n/t$ such generators in all. Each reduced basis has one
polynomial per minimal initial generator. Taking the union over (1)
therefore gives

$$
|\mathcal U(I)|=O\!\left(n^{85/24}\log^{3/8}(n+1)\right). \tag{9}
$$

The packing and union argument are established; (1) is the changed
quantitative input. This is a cardinality bound, not a field-operation,
coefficient-encoding or smooth-number unranking speed claim.

The theorem concerns the same positive weighted exponent lattice as
the original smooth-number interleavings. It bounds how much algebraic
dependencies can enlarge the family of possible standard sets, while
the shadows retain a quantitative arithmetic dependence. It does not
give a theorem about the distribution of primes. The gap to the generic
three-variable scale and the current-priority audit remain open.

## Verification and prior-art boundary

The greedy token/convex-chain viewpoint is classical; see
[Eppstein, Section 4](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf).
The universal ground, state dictionary, historical $O(n^3)$ bound and
universal-basis union argument are from
[Babson--Onn--Thomas, Theorem 2.4 and Corollary 2.5](https://arxiv.org/pdf/math/0207135).
The lattice triangle and short integer relation arguments above are
elementary classical inputs. The product cover, height-dependent
charging theorem (2), and its combination (8) are the exact statements
whose priority must be assessed. An absence of matching search results
does not establish novelty. The unresolved older-paper comparison is
preserved; this proof does not rely on its contents being unknown.

Dyadic families of weights already appear in
[van der Hoeven--Lecerf, *Amortized multi-point evaluation of multivariate
polynomials*](https://www.texmacs.org/joris/genamp/genamp.pdf), Sections
3.2--3.4 and 4.1. Their polylogarithmic family selects polynomials for
heterogeneous quasi-reduction, rather than meeting every initial cone.
Quasi-reduced monomials need not be a quotient basis. Thus this does not
supply a polylogarithmic projection cover or an initial-ideal count.
The inspected PDF is dated December 26, 2021; the
[HTML version](https://www.texmacs.org/joris/genamp/genamp.html) is dated
June 6, 2025, and version equivalence has not been assumed. The dyadic
construction is relevant prior art for the cover argument.

Two further primary comparisons have limited, explicit scope.
[Black, *Small Shadows of Lattice Polytopes*, arXiv v1 (2022),
Theorems 1.1--1.7 and Corollary 2.2](https://arxiv.org/pdf/2204.09129)
concerns selected short coherent paths and monotone diameter. Those
statements do not give the height-dependent bound for every rational
shadow in (2), or the all-vertex incidence count (6).
[Hausbrandt--Ruzika, Corollary 3.4](https://arxiv.org/html/2503.08178v1#S3)
bounds unrestricted multi-parameter matroid cells by the pairwise
comparison arrangement, with $O(m^{2p})$ for $m$ elements and $p$
parameters. Its statement does not give (2). These are comparisons
of specified statements, not certification of the strongest current
bound or of the absence of an equivalent theorem elsewhere.

Reproduce the supporting exact controls from the repository root:

```bash
python3 -B scripts/run_rational_shadow_bound_checks.py
```

The proof, rather than the finite range of these controls, establishes
the uniform bounds. The active goal remains incomplete.

The report checks 340 quadrant sweeps on 17 matroid/quotient models
and five rational kernel directions, replaying every comparison
interval by direct greedy selection. It independently enumerates
40 full basis-profile hulls, tests 1,330 short integer relations,
and certifies 28 vertices of five lattice polytopes by strict
objectives in short rational shadow planes. The controls also
exercise 15 oblique-lattice direction examples, including thin
unimodular shears, and verify 60 weighted dyadic pair sums with
exact rational cube-root enclosures. In dimensions two through
four, 1,464 projected ground points are checked against the
cover and its shared boundary-line bound. None of these finite
numbers is used to infer the uniform growth law.
