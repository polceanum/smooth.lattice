# An arithmetic consequence of the unrestricted counting bound

The rational-shadow theorem can be used to enumerate initial ideals,
rather than merely bound how many exist. This note gives an exact
reduction in the original coordinates, over every field. The polytope
reconstruction, cone traversal and greedy linear algebra are published
tools. They are attributed below; no new general traversal algorithm or
fastest-known implementation is claimed.

Let $I\subset K[x,y,z]$ have positive colength $n$, presented by one
marked reduced Gröbner basis under a known global monomial order,
and let $g=g(I)$. Write
$L=1+\log(n+1)$, and let $\omega$ be any feasible exponent for matrix
multiplication over $K$. There is an algorithm that returns all initial
standard sets, the state polyhedron and all marked reduced Gröbner
bases using

\[
 O(n^3L^2+(g+1)n^\omega L)                               \tag{O1}
\]

operations in $K$. Their union is a universal Gröbner basis; the
output may represent the union as a list of reduced bases with repeated
polynomials. No minimum-size or canonical duplicate-free output is
asserted. Direct traversal of the positive cones uses

\[
 O((g+1)n^{5/3}L)                                       \tag{O2}
\]

exact rational geometric operations and cost comparisons. The input
assumption matters: computing the first Gröbner basis from arbitrary
generators is an additional task. Alternatively one can start with
commuting coordinate-multiplication matrices and a cyclic vector.
Comparisons with the supplied input order are counted as comparisons.

Combining (O1)--(O2) with the
[proved unrestricted rational-shadow bound](rational_shadow_bound.md)
gives

\[
 \boxed{O\!\left(n^{\omega+23/8}L^{11/8}\right)}           \tag{O3}
\]

field operations, and also bounds the combined number of field and exact
geometric arithmetic operations. The latter budget is now
$O(n^{109/24}L^{11/8})$, below the field term for every feasible
$\omega\ge2$. The previous quadratic-in-$g$ hull-rebuilding budget is
retained as a reference below. This is an arithmetic bound, not a
coefficient-bit-complexity, benchmark or
prime-random-access claim. Its current-best status and exact priority
are unestablished. It does not resolve the sharp $G_3$ question.

## A bounded polytope with only relevant vertices

Use the established finite ground

\[
 E=\{a\ge0:(a_x+1)(a_y+1)(a_z+1)\le n\},\qquad
 U=E\cup\bigcup_i(E+e_i).
\]

Its sizes are $O(nL^2)$, and $U$ contains every minimal border head
of an initial standard set. Let $B(I)$ be the convex hull of profiles
of rank-$n$ monomial bases on $E$. The
[Babson--Onn--Thomas state dictionary](https://arxiv.org/pdf/math/0207135)
gives

\[
 \mathcal S=B(I)+\mathbb R_{\ge0}^3
   =\operatorname{conv}\{\mu(B):B\text{ initial}\}
        +\mathbb R_{\ge0}^3.                              \tag{O4}
\]

Each element of $E$ has total degree at most $n-1$, so all basis
profiles have total degree at most $n(n-1)$. Set $R=n^2+1$ and
truncate only the unbounded part:

\[
 Q=\mathcal S\cap\{u_x+u_y+u_z\le R\}.                    \tag{O5}
\]

This is a full-dimensional integer polytope. Its vertices below the
cut are exactly the $g$ initial profiles. Every new vertex is the
intersection of the cut with an unbounded edge of $\mathcal S$.
Such an edge starts at a state vertex and has direction one of the
three coordinate rays. There is at most one such edge in each
direction at each vertex. Hence

\[
 v(Q)\le4g,\qquad f(Q)\le2v(Q)-4.                         \tag{O6}
\]

All vertex coordinates are integers between zero and $R$. This
construction also handles a single state vertex or a flat basis
polytope, without applying a spatial volume bound to a flat hull.

## A vertex oracle with explicit tie handling

Suppose first that $w\in\mathbb Q^3$ has $M=\max_i w_i\le0$.
Maximizing $w\cdot u$ over (O5) uses no positive slack and reduces
to minimizing $-w\cdot\mu$ over the column matroid. To select a
**vertex**, rather than an arbitrary point on a tied face, maximize
lexicographically

\[
 (w\cdot u,-u_x,-u_y,-u_z).
\]

Matroid greedy selection with monomial cost key

\[
 (-w\cdot a,a_x,a_y,a_z)                                  \tag{O7}
\]

does precisely this. Every variable has positive lexicographic cost,
so the selected basis is an actual global initial standard set.

If $M>0$, let $j$ be the **largest index** achieving $M$, with
$x,y,z$ indexed $1,2,3$. All slack goes into coordinate $j$.
The correct greedy key is

\[
 \left(M|a|-w\cdot a,\quad
       a_x-|a|\mathbf1_{j=x},\quad
       a_y-|a|\mathbf1_{j=y},\quad
       a_z-|a|\mathbf1_{j=z},\quad a_j\right).             \tag{O8}
\]

For its selected basis of profile $\mu$, return

\[
 v=\mu+(R-|\mu|)e_j.                                     \tag{O9}
\]

Here $|\mu|$ denotes the sum of its coordinates. To verify both the
key and vertex property, replace $w_i$ symbolically by
$w_i-\epsilon^i$. Its unique largest component is indexed by $j$.
The cost of a monomial after allocating slack is

\[
 \sum_i a_i(w_j-\epsilon^j-w_i+\epsilon^i).
\]

Its coefficients through degree three in $\epsilon$ are the first
four entries of (O8). The additional $\epsilon^4a_j$ breaks the
remaining tie along the ignored coordinate. It leaves (O9) unchanged.
All variable costs are now positive; for $i\ne j$ this follows from
the largest-component rule, and for $j$ from the additional term.
Thus greedy again returns a genuine initial basis. The first four
objectives select a unique profile of $Q$, hence a vertex. A
sufficiently small positive rational $\epsilon$ realizes these
finite comparisons; the implementation compares tuples exactly
instead of choosing a floating tolerance.

These formulas hold for arbitrary, zero and boundary directions.
They use field arithmetic only for independence. No ordering of the
field and no characteristic-zero assumption are used.

## An alternative oracle reconstruction

Find one initial profile $\mu_0$ by positive-cost greedy selection.
The tetrahedron with vertices $\mu_0$ and
$\mu_0+(R-|\mu_0|)e_i$ is a known full-dimensional subset of $Q$.
Reconstruct $Q$ by querying each unverified inner-hull facet normal.
If the support point lies outside that facet, insert the new vertex;
otherwise the plane is a verified supporting facet of $Q$. Keep a
record of verified planes, so that a coplanar triangulation cannot
charge a facet repeatedly.

Every unsuccessful facet test discovers a previously unknown target
vertex; every successful one certifies a target facet. Once all
inner-hull facets are verified, their halfspaces contain $Q$, so the
inner hull equals $Q$. Thus at most

\[
 1+v(Q)+f(Q)\le12g+1                                     \tag{O10}
\]

greedy vertex queries are needed, including initialization. This is
the published oracle reconstruction principle of
[Emiris--Fisikopoulos--Konaxis--Peñaranda, Algorithm 1 and Lemmas 4--6](https://vissarion.github.io/publications/EFKP_IJCGA_13.pdf),
not a new general convex-hull method. Retaining the initial standard
basis label with each returned profile recovers all $g$ labels below
the cut. The cut facet is discarded to recover $\mathcal S$ and its
normal fan. It never meets a state vertex, so its incident normal
cones are unchanged.

## Direct traversal removes the quadratic geometry overhead

The positive weight triangle is subdivided into $g$ convex initial-ideal
cells. Its adjacency graph is simple, planar and connected: two convex
cells share at most one edge, and a generic path in the open triangle
joins their interiors by crossing edges. It has at most $3g$ undirected
adjacencies, including the small cases.

Cone descriptions, flippable facets and neighboring-basis traversal are
established in [Fukuda--Jensen--Thomas, Proposition 2.6, Sections 4.1--4.2
and Algorithm 4.2](https://arxiv.org/pdf/math/0509544). Those statements
and displayed arguments were inspected. Here a neighboring basis is
recomputed with the fast conversion backend, rather than assigning a new
complexity to a general Gröbner flip. The input ideal need not admit a
positive homogeneous grading. This is a composition of existing tools.

### The size of a cone description

Let $D$ be an $n$-point standard set and $H$ its minimal excluded
monomials. For each coordinate $i$, the proper downward $i$-chains
from heads with $a_i>0$ are disjoint subsets of $D$: two heads with
the same other coordinates would be comparable. Thus
$\sum_{a\in H}a_i\le n$. Put $t=\lceil n^{1/3}\rceil$.
Heads with all coordinates below $t$ number at most $t^2$, by antichain
projection; the rest number at most $3n/t$. Hence $|H|=O(n^{2/3})$.

There are at most $n$ terms in each reduced-basis tail, giving
$M=O(n^{5/3})$ head-minus-tail inequalities. Divisor chains bound
standard-monomial degrees by $n-1$ and minimal-head degrees by $n$.
Every inequality normal therefore has integer coordinates of magnitude
at most $n$.

Add nonnegativity and intersect with $\sum w_i=1$. Compute this polygon
in $O(M\log(M+1))$ exact operations as follows. Translate a known strict
interior point $p$ to zero and divide each affine inequality by its positive
slack there, obtaining $v_j\cdot(x-p)\le1$ in two coordinates. Sort the
$v_j$ and compute their hull using the two monotone-chain stacks; after
sorting each point is pushed and popped at most once. Intersect consecutive
hull-vertex supporting lines to obtain the polar polygon. The triangle
constraints ensure boundedness and put zero inside the dual. Redundant and
parallel inequalities and flat profile polytopes are allowed.

### An explicit neighboring weight

For an edge outside the coordinate boundary, let $p$ be its midpoint,
$D_p$ their least common denominator, and $\eta$ an inward
current-cone normal tight along the edge. Put

$$
 u=(\eta_x+\eta_y+\eta_z)(1,1,1)-3\eta,\qquad
 W=2(1+18n^2)D_p p+u.                                  \tag{O10a}
$$

Here $\sum u_i=0$, $\eta\cdot u<0$ and $|u_i|\le6n$.
The weight $W$ is integral and strictly positive. For any integer normal
$\zeta$ with $|\zeta_i|\le n$, a nonzero $\zeta\cdot p$ has
magnitude at least $1/D_p$, while $|\zeta\cdot u|\le18n^2$.
The perturbation therefore preserves every nonzero such comparison.
It does not require avoiding inactive comparison walls through $p$.

Consider the actual other cone along this edge. Its strict inequalities
at $p$ stay strict by this estimate. Its tight inequalities are parallel
to $-\eta$, because $p$ lies in the relative interior of the shared edge.
Those become strict on the other side. Thus $W$ lies inside precisely
that neighboring cone.

Polygon endpoints are normalized primitive rays of height $O(n^2)$,
by the two-row cofactor bound. Hence $D_p=O(n^4)$ and
$\|W\|_\infty=O(n^6)$. Weights and fixed-dimensional geometry have
$O(L)$ bits. A generic initial weight is $(1,n+1,(n+1)^2)$:
base-$(n+1)$ encoding distinguishes every reduced-basis comparison.

### Queries and output

Start at that weight, compute its reduced basis and full cone, and query
across every nonboundary edge by (O10a). Store discovered standard sets
and process each once. Connectivity proves completeness. There is one
initial query and two directed queries per undirected adjacency, at most
$6g+1$ altogether. The Neiger--Schost backend below proves (O1).
Cone computations, inequality scans and deterministic dictionary lookups
give (O2); the existing $O((g+1)nL^2)$ oracle sorting cost is absorbed.
A balanced dictionary on sorted standard sets suffices, using
$\log(g+2)=O(L)$ from the classical cubic count.

The polygons give the positive fan directly. Their nonnegative primitive
vertex rays give all state-polyhedron facet normals, including normals
with zero coordinates. Obtain each supporting constant from its incident
profile. There are $O(g)$ such rays, so no spatial hull rebuilding is
required. Coefficient bit growth over $K=\mathbb Q$ remains separate.

## Arithmetic and geometric costs

First obtain the three multiplication matrices and the cyclic vector.
From one reduced basis, the unconditional column construction in
[Babson--Onn--Thomas, Proposition 4.1](https://arxiv.org/pdf/math/0207135)
costs $O(n^3L^2)$ field operations on $U$. Its columns include those
needed for the multiplication matrices in the input standard basis.
The same budget can generate columns from supplied matrices in divisor
order. This is a one-time cost.

For each query apply
[Neiger--Schost, Theorem 1.7 and Propositions 3.7, 3.14](https://cs.uwaterloo.ca/~eschost/publications/multi-krylov.pdf)
with $D=n$, $r=3$, $m=1$ and generator $1$. It returns the ordered
monomial basis and the reduced relations in $O(n^\omega L)$ field
operations. The theorem permits arbitrary commuting matrices, fields
and module monomial orders. In particular it does not require Jordan
form. Their separate fast construction of multiplication matrices,
Theorem 1.9, has an extra hypothesis on the input leading module; we
do not apply it to arbitrary original-coordinate ideals.

Doubling matrix powers and retaining ordered rank profiles uses
$O(L)$ stages, each with at most $2n$ rows after initial compression.
Retained **monomial labels** must stabilize; unchanged rank is
insufficient. Batch recovery of at most $3n$ border remainders costs
$O(n^\omega)$. The resulting relations belong to $I$; inclusion of
leading ideals and equal colength prove that they form its reduced
Gröbner basis.

This is the published compressed Krylov method, not a new elimination
algorithm. Its Theorem 1.7, Algorithm 1, correctness Lemmas 3.8--3.10
and the displayed complexity and recovery proofs were inspected.
The checker implements the compression with elementary arithmetic;
it does not implement fast matrix multiplication or certify timings.

Equations (O6) and (O10) give an alternative proof of (O1). The older
oracle reconstruction uses the reference geometry budget

\[
 O((g+1)^2\log(g+2)+(g+1)nL^2).                           \tag{O2ref}
\]

For this reference method, rebuild a
three-dimensional hull after each insertion in $O(k\log(k+2))$
geometric operations on $k=O(g+1)$ points, using a standard exact
hull routine. Summing gives $O((g+1)^2\log(g+2))$. The
[Chan 1996 algorithm, Sections 3--4](https://www.cs.jhu.edu/~misha/Spring20/Chan96.pdf)
provides the required bound and explicitly handles degeneracies;
perturbing nonvertices into vertices is unnecessary. The compressed
oracle sorts $O(n)$ labels at $O(L)$ stages, contributing
$O(nL^2)$ comparisons per query. Alternatively, the fixed module
below also has an initial $O(nL)$ labels to sort once per query,
within that comparison budget.
All polyhedral input points and newly queried facet objectives have
polynomially bounded integer height in $n$. These operation counts still do not bound
the bit growth of coefficients in $K=\mathbb Q$.

## What the comparison establishes

| Inputs to the composition | Resulting field-operation bound |
|---|---|
| Babson--Onn--Thomas Theorem 4.2 itself, in three variables | $O(n^7L^{10})$ |
| Classical $g=O(n^3)$, oracle reconstruction and dense conversion | $O(n^6L^2)$ |
| Classical count, same reconstruction and Neiger--Schost conversion | $O(n^{\omega+3}L)$ |
| Rational-shadow count, same reconstruction and Neiger--Schost conversion | $O(n^{\omega+23/8}L^{11/8})$ |

The last two rows include an additional one-time $O(n^3L^2)$ input
conversion cost. The previous version compared only the first, second
and dense shadow bound $O(n^{47/8}L^{19/8})$. That comparison omitted
a stronger applicable published backend. Apply the same backend to
both counts: the shadow theorem supplies only the further $1/8$
power saving, with the displayed logarithmic change. All rows retain
the originally specified coordinates. None asserts fastest-current
status.

The older geometry budget (O2ref) is $O(n^6L)$ with the classical count
and $O(n^{23/4}L^{7/4})$ with the shadow count. Direct traversal (O2)
reduces these to $O(n^{14/3}L)$ and $O(n^{109/24}L^{11/8})$.
Both are below their fast field terms. Thus (O3) now bounds field and
geometric unit-cost arithmetic together. No coefficient-bit bound follows.

[Fukuda--Jensen--Thomas, Section 5 and Theorem 5.1](https://data.math.au.dk/publications/pp/2005/imf-pp-2005-15.pdf)
already analyze output-sensitive Gröbner-fan traversal. That section
assumes a positive homogeneous grading and leaves the flip cost as
a separate parameter. We inspected those hypotheses and the
displayed complexity proof; we do not claim that Gfan lacks other
applicable methods, or that this is a complete later-algorithm audit.
The oracle paper's full resultant analysis and Chan's other
applications were not audited here.

The result is a concrete algorithmic consequence for arbitrary
finite quotients in their original coordinates. The sharp count,
strongest-current comparison, priority of the rational-shadow input,
and the original publication requirements remain unresolved. No
practical speedup or new manuscript is claimed.

## One fixed module preserves every quotient's dependencies

There is an exact univariate formulation of the remaining counting
problem. It uses the actual quotient columns, including their special
dependencies. Set

\[
 J_n=\{(b,c)\ge0:(b+1)(c+1)\le n\},\qquad
 m=|J_n|=\sum_{j=1}^n\lfloor n/j\rfloor=O(nL).
\]

For each $(b,c)\in J_n$ introduce $e_{bc}$ and the single, fixed map

\[
 \Phi:K[T]^m\longrightarrow A,
 \qquad T^a e_{bc}\longmapsto M_x^aM_y^bM_z^c1.            \tag{O11}
\]

Every global standard set lies in $E$ and hence among these chains,
so $\Phi$ is onto. Its kernel has rank $m$ and quotient dimension
$n$. For a global monomial order, order $T^ae_{bc}$ by the original
monomial $x^ay^bz^c$. This is a valid module order: multiplication
by $T$ preserves comparison and increases each chain. Restricting
greedy selection to a ground containing the entire global greedy
basis leaves that basis unchanged. Therefore the module standard
set is exactly the original standard set. Its chain lengths
$d_{bc}$ have sum $n$.

For positive costs $(\alpha,\beta,\gamma)$ these module orders are
weighted degree orders with shifts
$(\beta b+\gamma c)/\alpha$ and the same tie rule. Thus the original
fan is a two-parameter family of degree profiles of **one** module,
not a different interpolation problem at every weight. With rational
costs one can multiply through to use exact integer comparisons;
the statement does not assume integral shifts. The three commuting
matrices and cyclic generator constrain which modules can arise.
Arbitrary interpolation modules are not thereby realized as quotients.

There is also an exact formulation for every cost-prefix rank. At
threshold $t$, put

\[
 \ell_{bc}=\max\left(0,\min\left(
 \left\lfloor\frac{n}{(b+1)(c+1)}\right\rfloor,
 1+\left\lfloor\frac{t-\beta b-\gamma c}{\alpha}\right\rfloor
 \right)\right).                                       \tag{O12}
\]

Use the module order $(a-\ell_{bc},(b,c))$, with any fixed position
tie rule. Every allowed chain prefix has shifted degree below zero,
and every other term has degree at least zero. The first independent
terms in that initial block span its image. Consequently its rank is
the number of module standard terms $T^ae_{bc}$ with $a<\ell_{bc}$.
Positive-cost partial greedy sets are downsets: if a divisor were
dependent on earlier columns, multiplying that relation would make
the monomial dependent as well. Their size is at most $n$, so they
lie in $E$. This proves that truncation in (O12) preserves the rank
of the full positive-cost prefix, including special and nonreduced
quotients.

Applying the same published univariate algorithm to (O11), with
$r=1$, $D=n$ and $m=O(nL)$, again costs $O(n^\omega L)$ per supplied
order. This faithful reduction and ordinary interpolation-basis
machinery do not bound how many degree profiles occur as shifts vary.
A new all-profile count remains necessary for the main goal.

## Exact convex completion, and why cyclic-span ranks lose the fan

For a general surjective map $K[T]^m\to V$, $\dim_KV=n$, let
$v_i$ be its seeds and $T$ its endomorphism. Write $D$ for the
integer allocations $d\ge0$, $\sum_i d_i=n$, for which
$\{T^av_i:0\le a<d_i\}$ is a basis, and put

$$
 q(d)=\sum_i\binom{d_i}{2},\qquad
 f(S)=\dim_K\sum_{i\in S}K[T]v_i.
$$

The nearest-site formulation is exact **on $D$**:

$$
 q(d)+s\cdot d
   =\tfrac12\|d+s\|_2^2-\tfrac12\|s\|_2^2-\tfrac n2.
                                                               \tag{M1}
$$

Its ordinary convex hull is the classical polymatroid base polytope

$$
 \operatorname{conv}D
 =P(f)=\{d\ge0:d(S)\le f(S),\ d([m])=n\}.                       \tag{M2}
$$

Here $f$ is submodular by the subspace dimension formula. For a
proof of the equality, the inclusion from left to right is immediate.
Order the seeds arbitrarily and take the successive rank increments
of their cyclic subspaces. The previously generated space is
$T$-stable. In its quotient, the next seed is cyclic, so its first
rank-increment-many powers form a basis. Thus every greedy vertex
of $P(f)$ belongs to $D$. The classical greedy description proves
(M2). The subspace-rank and base-polytope framework is established;
compare [Brandenburg--Loho--Smith, Sections 1.1 and 2.1](https://arxiv.org/html/2403.07751v1#S1.SS1).

One cannot minimize $q+s\cdot d$ over all integer points of $P(f)$.
The following exact completion retains the missing information. On
the finite chain ground $E'=\{(a,i):0\le a<n,\ 1\le i\le m\}$,
define

$$
 v(d)=\min\left\{\sum_{(a,i)\in B}a:
      B\text{ is a column basis},\ |B\cap E'_i|=d_i\right\},       \tag{M3}
$$

with value $+\infty$ if there is no such basis. This is aggregation
of a linear cost on the bases of one ordinary matroid. Hence $v$
is an M-convex function, by the established aggregation theorem
[Murota, Proposition 4.6(3), equation (4.4)](https://arxiv.org/pdf/2001.03018v4).
Its domain is $P(f)\cap\mathbb Z^m$, and its convex extension is the
lower convex envelope of $q$ on $D$. In particular,

$$
 q(d)\le v(d)\le\binom n2,\qquad
 v(d)=q(d)\ \Longleftrightarrow\ d\in D,                        \tag{M4}
$$

for integer $d\in P(f)$, and for every real shift $s$,

$$
 \min_d(v(d)+s\cdot d)=\min_{d\in D}(q(d)+s\cdot d).             \tag{M5}
$$

To justify the convex-envelope assertion, first minimize the basis
cost $a+s_i$, breaking ties by a fixed position order. Greedy selection
is closed under division by $T$: a dependency expressing an omitted
term by earlier terms, multiplied by $T$, precludes its successor
from being selected. Cayley--Hamilton excludes powers $a\ge n$.
This gives (M5), so both lower convex hulls agree. The aggregated
M-convex function is convex-extensible and has no domain holes;
thus its values at integer points equal that envelope. For a fixed
allocation the smallest possible sum of distinct powers in chain $i$
is $\binom{d_i}{2}$, with equality only for its prefix. This proves
the lower bound and contact criterion in (M4). Every point of
$P(f)=\operatorname{conv}D$ is a convex combination of allocations
with $q(d)\le\binom n2$, proving the upper bound. M-convexity,
aggregation and convex extensibility are classical tools, not new
optimization results. The quadratic on $D$ need not itself satisfy
the unit exchange axiom.

Every $d\in D$ is uniquely selected by the nonnegative integer shift
$s_i=n-d_i\in[0,n]$. Its $n$ prefix terms have shifted cost below
$n$, and every other term has cost at least $n$. Thus $D$ is exactly
the set of full-shift module standard profiles, not merely a larger
feasible set. In (M7), even bounded shifts permit exponentially many
profiles since $m\ge n$. An algorithm for one supplied shift does
not provide a polynomial all-shift enumeration theorem. A counting
argument for the original problem must retain the special
two-parameter plane, rather than enumerate every module shift.

This is consistent with the published shifted-Popov invariance result:
[Jeannerod--Neiger--Schost--Villard, Definitions 1.1--1.2 and Lemmas 3.3--3.4](https://arxiv.org/html/1602.00651v1#S3)
identify pivot degrees for a fixed shift and fixed row module. The
same abstract quotient invariant factors do not imply the same embedded
row module or the same pivot degrees. The displayed proofs of those
lemmas were inspected; no new degree-invariance theorem is claimed.

### Same operator and rank polytope, one initial ideal or every corner cut

The loss in (M2) is already present for actual radical point ideals.
For $n\ge2$, in characteristic zero put
$F_n(x)=\prod_{j=1}^n(x-j)$ and consider

$$
 I_0=(F_n(x),y-1,z-1),\qquad
 I_1=(F_n(x),y-x^n,z-x^{n^2}).                                \tag{M6}
$$

Both algebras have length $n$, with the **same** operator
$T=\operatorname{diag}(1,\ldots,n)$ and cyclic vector $1$ in their
evaluation representations. For the complete seed set $J_n$ of
(O11), every individual seed is cyclic: its nonzero evaluations
multiply an ordinary Vandermonde matrix by an invertible diagonal
matrix. Hence both have
$f(\varnothing)=0$, $f(S)=n$ for every nonempty $S$, the same
polytope $P(f)=n\Delta_{m-1}$, and the same $K[T]$-quotient
invariant factors, with single nonunit factor $F_n$.

For $I_0$ all seeds are identical. Its only prefix allocations are
$ne_i$, while $v_0(d)=\binom n2$ on every integer composition.
Any two occupied prefixes duplicate their constant term; conversely,
any allocation can distribute the independent powers $0,\ldots,n-1$
among its chains. Since $y,z$ reduce to constants, $g(I_0)=1$.

For $I_1$, a monomial with all coordinates below $n$ has column
$(j^{a+nb+n^2c})_{j=1}^n$. The exponent codes are distinct. A
generalized Vandermonde matrix at strictly increasing positive nodes
is nonsingular: a nonzero polynomial with $r$ nonzero terms has at
most $r-1$ distinct positive roots, by Rolle's theorem inductively
after dividing by its lowest power. Thus every prefix allocation
is feasible, $v_1=q$, and

$$
 |D_0|=m,\qquad |D_1|=\binom{n+m-1}{m-1},\qquad
 g(I_1)=C_3(n).                                               \tag{M7}
$$

Indeed all columns on the universal ground $E$ are in uniform
position, so every global order selects exactly its first $n$
monomials. The last equality is with the established corner-cut
class; it is not a new counting theorem. All of (M6)--(M7) also
have a version in every characteristic over $k(t_1,\ldots,t_n)$:
replace $j$ by algebraically independent $t_j$ in $F_n$ and the
evaluation nodes. Distinct exponent codes give distinct determinant
monomials with coefficients $\pm1$, which remain nonzero even in
characteristic two. No positive-root argument is used there.

The incorrect integer relaxation can have an unbounded value gap
even on the original positive two-parameter shift plane. At
$w=(n^2,1,n+1)$, divide costs by $w_x=n^2$. The relaxed first
$n$ terms are $1,y,\ldots,y^{n-1}$ and their value is
$\binom n2/n^2$. The true value for $I_0$ is $\binom n2$, a
factor $n^2$ larger. The exact completion (M3) gives the true value.

Thus knowing cyclic-span ranks, invariant factors, or the allocation
polytope does not preserve the positive fan, even with the operator
held fixed. The useful remaining object is the full valuation $v$,
with its contact condition (M4), restricted to
$s_{bc}=\beta b+\gamma c$ for $\beta,\gamma>0$. A new collective
bound must exploit that valuation and the two additional commuting
multiplications, or return to the rank-preserving hyperbolic incidence
budget. The classical M-convex reformulation alone supplies no such
bound. No unrestricted count, separating lower family, novelty or
publication claim follows from this obstruction.

## Arithmetic seed positions give exact parking valuations

The full valuation can be described without quotient linear algebra in
the uniform arithmetic families. This gives a concrete constraint beyond
the common simplex in (M2), rather than replacing the degree costs by
$q$. The parking process itself is classical. No new unrestricted count
or general parking theorem is asserted.

### Positive monomial-curve samples: a line with release times

Take $n$ distinct positive rational parameters $t_1<\cdots<t_n$ and
the points $(t_j,t_j^A,t_j^B)$, with nonnegative integers $A,B$.
The point ideal has length $n$ since the $x$ coordinates are distinct.
Every set of at most $n$ monomial columns with distinct encoded degrees
$a+Ab+Bc$ is independent, by the generalized Vandermonde argument
above. For the seed $(b,c)$ write $r_{bc}=Ab+Bc$.

More generally let $r_1\le\cdots\le r_m$ be integer seed positions,
and let $d_i\ge0$, $\sum d_i=n$. In (M3), choosing $d_i$ terms
from seed $i$ means assigning its $d_i$ jobs to distinct integer slots
$r_i+a$, $0\le a<n$, with displacement cost $a$. All $n$ slots
must be distinct across seeds. Put

$$
 E_0=r_1,\qquad L_i=\max(0,E_{i-1}-r_i),\qquad
 E_i=\max(E_{i-1},r_i)+d_i.
$$

Then the exact degree valuation is

$$
 \boxed{v(d)=q(d)+\sum_{i=1}^m d_iL_i.}                       \tag{M8}
$$

Equivalently, sort the $n$ releases, with their multiplicities, as
$\rho_1\le\cdots\le\rho_n$; put
$u_1=\rho_1$, $u_j=\max(\rho_j,u_{j-1}+1)$. Then
$v(d)=\sum_j(u_j-\rho_j)$. To prove optimality, uncross any assigned
pair whose release order is reversed: the swap is feasible and
preserves the sum of slots. For increasing releases, the earliest-slot
recursion is componentwise minimal among increasing feasible slot lists.
The possible displacements are at most $n-1$, so the finite chain ground
in (M3) suffices. Grouping equal-seed releases gives (M8), including
zero allocations and coincident seed positions.

Consequently $v=q$ exactly when the occupied intervals
$[r_i,r_i+d_i-1]$ are pairwise disjoint; empty intervals are ignored.
The defect $v-q$ is the explicit inter-seed waiting cost. Both the
identical-seed and separated-seed examples in (M6) are special cases.
These interval constraints encode the coherent arithmetic starts
$Ab+Bc$ of the two additional commuting multiplications.

### Cyclic orbit segments: circular parking

Let $M\ge n$, let $K$ contain a primitive $M$th root $\zeta$, and
assume $\operatorname{char}K$ does not divide $M$. Consider the
$n$ points

$$
 (\zeta^j,\zeta^{Aj},\zeta^{Bj}),\qquad 0\le j<n.              \tag{M9}
$$

Their ideal is
$(\prod_{j=0}^{n-1}(x-\zeta^j),\ y-x^A,\ z-x^B)$.
Columns have residues $a+Ab+Bc\pmod M$; any $n$ distinct residues
give an ordinary Vandermonde matrix in the distinct nodes
$\zeta^{a+Ab+Bc}$. Thus this is an actual commuting quotient with
the uniform distinct-residue rank rule. It includes composite $M$.

Give $d_i$ jobs preferred slot $r_i=Ab_i+Bc_i\pmod M$ in a circle
of $M$ slots. In any order, park each job at its first unoccupied
slot clockwise. Its displacement is the least nonnegative number of
steps from its preference. Then

$$
 \boxed{v(d)=\text{total circular-parking displacement}.}     \tag{M10}
$$

Here is a deterministic optimality proof, also covering $n=M$.
Fix any job order. After earlier jobs have been fixed greedily, let
$c$ be the first available slot for the next job and $t$ its slot in
an optimal completion. If $c$ is unused in the completion, moving
the job there reduces or preserves cost. Otherwise exchange it with
the job at $c$. The first job's saving is the clockwise distance
from $c$ to $t$, while the second job's added displacement is at
most that distance, by the triangle inequality for directed circular
distance. So there is an optimal completion agreeing with this next
greedy choice. Induction proves (M10) for every insertion order.
Before the $j$th insertion only $j-1$ slots are occupied, so its
displacement is at most $j-1<n$. This again justifies (M3)'s finite
ground; it also constructs the bound $v\le\binom n2$ directly.

The contact criterion is now disjointness of the cyclic arcs
$r_i+\{0,\ldots,d_i-1\}\pmod M$. If these are disjoint they
give a basis of cost $q$; conversely equality forces each chain to
use its own first $d_i$ powers. Classical linear probing and its
parking interpretation are described in
[Flajolet--Poblete--Viola, pp. 490--491](https://algo.inria.fr/flajolet/Publications/FlPoVi98.pdf).
Their random-input moment and limit laws are not used here and do not
bound the number of coherent profiles in our deterministic weight plane.

### The remaining counting problem is still two-dimensional

For both models, substitute the exact parking valuation into (M5) with
$s_{bc}=\beta b+\gamma c$, $\beta,\gamma>0$. The original count
is the number of uniquely optimal profiles of this specific integer
parking cost on that plane. Arbitrary release positions or all module
shifts are a larger problem. The formula gives a faithful arithmetic
model for the growing-parameter family; it does not yet bound that
count by the generic scale or produce a separating family.

It is also not the rank rule of the whole repeated-root family below.
For example, $n=5$, characteristic two, $M=8$, the four distinct
residues $\{0,2,4,6\}$ have rank three. Uniform circular parking
would wrongly assign them rank four. The exact prime-digit rank law
must be retained there. Although its binary matroids are series--parallel,
that label alone supplies no new two-parameter count. For perspective,
[Eppstein's WADS 2021 slides, slide 16](https://ics.uci.edu/~eppstein/pubs/Epp-WADS-21-slides.pdf)
list a one-parameter min/max-formula breakpoint question as unresolved.
This scoped historical comparison is not a claim about the strongest
current bound, or an equivalence with the present coherent problem.

## An exact prime-digit rank law in a nonreduced family

Let $K$ have prime characteristic $p$ and consider the actual length-$n$
quotient

\[
 A=K[X,Y,Z]/((X-1)^n,XY-1,Z-Y^k),\qquad k\ge1.            \tag{O13}
\]

With $U=1+t$ in $K[t]/t^n$, the three generators are
$X=U,Y=U^{-1},Z=U^{-k}$. A monomial maps to
$U^{a-b-kc}$. Choose any $M=p^s\ge n$. Since $U^M=1$, its column
depends only on that integer exponent modulo $M$. For a subset
$S\subseteq\{0,\ldots,M-1\}$, define

\[
 R_{n,M}(S)=\operatorname{rank}_K
       \big(\tbinom{\lambda}{j}\big)_{0\le j<n,\lambda\in S}.
\]

For $0<n<M$, write $h=M/p$ and $n=qh+r$, $0\le r<h$. Let
$l_u=|\{v\in\{0,\ldots,p-1\}:u+vh\in S\}|$. Then

\[
 \boxed{R_{n,M}(S)=\sum_{u=0}^{h-1}\min(q,l_u)
       +R_{r,h}(\{u:l_u>q\}).}                          \tag{O14}
\]

The terminal rules are $R_{0,M}=0$ and $R_{M,M}(S)=|S|$.
In particular (O14) is a rank identity, not a rule based only on
which matrix entries are nonzero.

Here is a full proof. Lucas's identity factors the row $j=ah+b$ and
column $\lambda=vh+u$ as
$\binom va\binom ub$ modulo $p$. The first $qh$ rows contain every
$0\le b<h$ for $0\le a<q$. The lower-digit Pascal matrix
$(\binom ub)_{0\le b,u<h}$ is invertible, so a row change separates
these into independent high-digit blocks, one for each $u$. A block
on $l_u$ distinct high digits has rank $\min(q,l_u)$: its rows are
polynomials of consecutive degrees $0,\ldots,q-1$ with invertible
leading coefficients in $K$, giving the Vandermonde rank.

The remaining $r$ rows have $a=q$, $0\le b<r$. If $l_u\le q$,
they add nothing in that fiber. Otherwise the degree-$q$ high-digit
row adds precisely one independent direction there, again by the
Vandermonde determinant. Modulo the first $qh$ rows, the new rows
therefore have coefficients $\binom ub$ on those fibers with
$l_u>q$. Their rank is exactly the last term in (O14). For $q=0$
the same argument simply identifies duplicate lower-digit columns;
for $n=M$ the whole Pascal matrix is invertible. This completes the
induction without any genericity assumption.

The same block argument describes **all** bases when every residue
is available. Each low-digit fiber supplies $q$ or $q+1$ columns;
exactly $r$ fibers supply the extra column, and those $r$ lower
labels must form a basis for $R_{r,h}$. Each fiber needs at least
$q$ columns to span its full block, and more than $q+1$ are dependent.
Thus the number $B(n,M)$ of residue bases satisfies

\[
 B(n,M)=B(r,h)\binom pq^{h-r}\binom p{q+1}^{r},\qquad
 B(0,M)=B(M,M)=1.                                       \tag{O15}
\]

For distinct arbitrary leaf costs, a minimum-cost basis takes the
$q$ cheapest entries in every fiber. It passes the $(q+1)$-st
cheapest cost in each fiber to the smaller rank-$r$ problem and adds
its chosen extra entries. This is exact recursion, not an exhaustive
search. For sparse rank queries, (O14) has just one smaller recursive
call at each digit. Sorting and grouping the supplied $|S|$ labels
at each level gives $O((|S|+1)\log(|S|+2)(s+1))$ integer comparison
and grouping operations and no field elimination. This count is for
an already supplied residue set; constructing cheapest monomial
representatives or enumerating a whole fan is a separate problem.

For example, in characteristic two, $n=5,M=8$ and
$S=\{0,2,4,6\}$ have rank three, although all four residues differ.
This exhibits actual special-sample dependence discarded by a generic
$\min(n,|S|)$ rank rule. For any positive-cost monomial prefix in
(O13), collecting its residues on $E$ and applying (O14) gives its
exact quotient rank.

The dependency code is exactly
$\langle(V-1)^n\rangle\subseteq K[V]/(V^M-1)$: a coefficient
vector lies in the kernel precisely when its polynomial vanishes
modulo $(V-1)^n$. Hence this is repeated-root cyclic-code mathematics.

The recursion has a further exact consequence. Define the homogeneous
subset-rank enumerator

\[
 H_{n,M}(A,B,Y)=\sum_{S\subseteq[M]}
       A^{M-|S|}B^{|S|}Y^{R_{n,M}(S)},
 \qquad [M]=\{0,\ldots,M-1\}.
\]

Grouping each fiber by its number of chosen columns in (O14) gives

\[
 H_{n,M}(A,B,Y)=H_{r,h}(A',B',Y),                         \tag{O16}
\]
\[
 A'=\sum_{l=0}^q\binom pl A^{p-l}B^lY^l,\qquad
 B'=\sum_{l=q+1}^p\binom pl A^{p-l}B^lY^q.
\]

Indeed, fibers with at most $q$ chosen entries contribute $A'$ and
do not enter the smaller rank query. The others contribute $B'$ and
their lower labels form precisely that query. Summing all patterns
is the right side of (O16). The terminal values are
$H_{0,M}=(A+B)^M$ and $H_{M,M}=(A+BY)^M$.

For a finite field of size $Q$ and characteristic $p$, the Hamming
weight enumerator of $C=\langle(V-1)^n\rangle$ is consequently

\[
 W_C(z)=H_{n,M}(1-z,Qz,Q^{-1}).                          \tag{O17}
\]

To see this directly, the number of codewords supported inside $S$
is $Q^{|S|-R(S)}$. Summing this number times
$z^{|S|}(1-z)^{M-|S|}$ over $S$ counts each word of weight $d$
with total contribution $z^d$. This proves (O17) without enumerating
the $Q^{M-n}$ codewords.

There is a division-free integer-polynomial version. Start with
$A=1-z,D=z$. At every digit replace

\[
 (A,D)\longmapsto\left(
 \sum_{l=0}^q\binom pl A^{p-l}D^l,\quad
 \sum_{l=q+1}^p\binom pl A^{p-l}Q^{l-q-1}D^l\right),       \tag{O18}
\]

and replace $(n,M)$ by $(r,h)$. The terminal expression is
$(A+QD)^M$ when $n=0$, or $(A+D)^M$ when $n=M$. This follows
by substituting $B=QD$, $Y=Q^{-1}$ in (O16); the new $B'$ is again
divisible by $Q$. Thus no coefficient divisions are needed.
Precomputing the powers through $p$ at each of at most $s$ steps
gives an arithmetic circuit of size $O(ps+\log(M+1))$ with integer
constants for the complete enumerator. This is a compact expression
or evaluation bound; expanding its coefficients and the bit cost of
large integers are additional tasks. The checker evaluates both
versions independently against explicit polynomial-generated codewords.

[Zhao--Li--Yang--Fu--Shum, arXiv:2304.00762v3](https://arxiv.org/html/2304.00762v3)
was inspected at its code definitions, Theorems 4, 6 and 16 and their
displayed decompositions, and the general-case statement of Theorem 17.
Those results compare related repetition, direct-sum and block structures;
their remaining distribution term must be compared with (O18) before
claiming an improvement. No claim that (O14)--(O18) are new is made.
The general coding literature, including
[Sobhani's matrix-product construction](https://doi.org/10.1016/j.ffa.2016.02.001)
and earlier squaring constructions, remains incompletely audited.
Only that article's publisher abstract and metadata were accessible;
its full theorem comparison remains open. An exponential number of these
bases also gives no coherent-fan lower bound. What remains useful
for this project is the exact prime-digit dependency recursion in
the original monomial coordinates. No reduction of the earlier coupled
three-point rank matrices to these single-prefix matrices is established.

~~~bash
python3 -B scripts/run_oracle_enumeration_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker independently enumerates all monomial bases in small
cyclic commuting models, compares every tested support query with
the resulting explicit polytope, reconstructs the whole truncated
state polytope, and recovers all positive cones and reduced bases.
It includes nonreduced and point quotients, flat state polytopes,
characteristic changes, zero objectives, ties and cut facets. Its
small hull routine enumerates triples for verification; the claimed
geometric bound uses the cited hull routine, not that slow checker.
Compressed bases are compared independently with full-ground
elimination. The fixed-module prefix formulation and prime-digit
recursion are checked against direct ranks, including arbitrary
Pascal subsets and polynomial columns of (O13). Exact witnesses,
all small residue-basis counts and source hashes are in the report.
The direct positive-cone traversal is checked separately against those
exhaustive state polytopes and the older reconstruction. Polar hulls are
compared with independent clipping, every selected neighbor is checked
along its whole shared edge, and all bounded integer comparison normals
are tested for the claimed sign preservation. Inactive midpoint ties are
explicitly retained. The two algorithms recover identical state halfspaces.
The additional module controls enumerate every basis on finite chain
grounds, recover the grouped minimum-degree costs, and compare them with
an independent exact convex-envelope calculation, the full unit exchange
axiom and shifted greedy selection. The paired actual ideals are checked
by integer determinants for every prefix allocation through length five;
symbolic determinant coefficients exercise characteristics two, three,
five and seven through length four. Those symbolic tests concern
rational-function-field realizations, not specialization to those small
finite fields. Four complete paired fans are independently reconstructed.
The parking controls additionally minimize every allocation cost by direct
integer determinants of quotient-column matrices in 88 line and circular
fixtures. They compare these minima with the closed line formula, every
distinct insertion order of each preference multiset, and the exact
prefix-contact criterion. Complete circles, composite periods, coincident
preferences and wraparound are included. A separate rank comparison at
length five distinguishes the nonuniform characteristic-two Pascal columns
from the uniform characteristic-seventeen orbit columns on the same four
residue labels. These checks verify bounded instances of (M8)--(M10);
their uniform proofs are above, and no growth law follows from the tests.
