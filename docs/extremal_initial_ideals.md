# Extremal counts of monomial initial ideals

The active target is a substantial advance on the number of different monomial
initial ideals of a finite polynomial algebra in three variables. The endpoint
is a sharp asymptotic theorem, or an infinite realizable family with a polynomial
separation from the generic scale and a further general consequence. A new
manuscript, a small exponent improvement, or finite examples do not meet that
endpoint. Mathematical importance and publication cannot be guaranteed in
advance; they must be supported by the eventual theorem and its comparison with
existing work.

The [goal record](../results/research/extremal_initial_ideals_goal.json) preserves
the acceptance criteria. Current results include a primary-source comparison,
uniform local bounds and exact distinctions between coefficient classes.
Generic torus-grid graphs reduce to ordinary cuts, while local rectangles
already have coherent dependencies that this reduction misses. Published
conditions for individual neighboring ideals add a further constraint, but
cannot by themselves count one fan. The sharp three-variable count
remains open. A [rational-shadow theorem](rational_shadow_bound.md) now
gives $G_3(n)=O(n^{23/8}\log^{3/8}(n+1))$, with a uniform proof and
an output-cardinality consequence. It does not complete the sharp target
or certify current novelty. The [dyadic planar proof](planar_dyadic_bound.md) establishes
$G_2(n)=\Theta(n\log(n+1))$ over all fields, with a general hyperbolic-weight
matroid bound and a smaller planar universal-basis cardinality bound.
Its exact priority remains under audit. It closes the planar power-gain
route: the proved product transfer cannot yield a spatial power separation
from planar branches. The fixed-degree carrier classification below remains
a sharper restricted comparison.
Finite gradings give a further classical obstruction: repeated exponent
spacings can keep the count linear even when the carrier degree grows.
A sharper Frobenius obstruction classifies the additive prime-power class:
its initial ideals are prime-power boxes, with a sharp polylogarithmic
count. This restricted classification is not the unrestricted target.
The [integer cone-height obstruction](initial_cone_height.md) separates
arithmetic certificate size from fan cardinality: an all-field family has
five initial ideals, a necessary integer weight height $\Theta(n^2)$,
and a necessary projection radius $\Theta(n)$, while exactly two planes
cover its cones. It rules out two individual-cone shortcuts without
excluding a stronger aggregate projection bound.
The [weighted-incidence reduction](weighted_facet_incidence.md) gives
three-plane boundary coverage, but its proposed cubic incidence budget
is now disproved. Finite quotients of a three-pole rational curve have
exactly $\binom{n+1}{2}$ initial ideals and weighted incidence
$\Theta(n^4)$, even for radical rational point ideals and with interior
cells alone counted. The classical incidence order is therefore sharp;
a stronger count requires more than a power saving in that moment.
A nine-cone quotient also keeps the sharp individual certificate and
radius exponents in a wholly interior cell. A selected shadow cover or
a threshold-dependent cell count remains possible; neither has yet
closed the main extremal gap.

The [nested-matroid construction](nested_jet_matroids.md) now gives an
explicit non-binomial family whose entire monomial rank rule is known,
including quotients of embedding dimension three. The ambient truncation
and actual quotient length are computed separately. This supplies a
constructive route to simultaneous dependencies, rather than assuming
that an arbitrary linear matroid admits commuting multiplication.
Classical ordered interpolation proves the rule. Balanced valuations reduce
to slab prefixes, whose fixed-height count is $\Theta_h(n^2\log(n+1))$
for $h\ge2$. A classical relative-box argument keeps their count within
the generic power scale even as height varies. For general valuations,
an infinite length-$M^2$ family has coherent bounding boxes of $M^3$
points and exact layer-rank defects. Thus neither lattice convexity nor
a uniform linear box product extends that slab proof. The large footprint
alone gives no counting lower bound. The continuation below controls the
equal-valuation family without assuming a small box for its standard sets.

The [chain-allocation continuation](nested_chain_allocation.md) describes
the curvilinear subfamily exactly by cumulative integer quotas and a
separable quadratic objective. Feasible unit transfers give all its
positive initial-cone inequalities. Along a valuation-price line,
submodular uncrossing proves that cumulative allocations increase and
each quota can cease to be tight only once. For $n=M^2$, $m=(1,M,M)$,
the stratum with every quota tight has only
$\Theta(M\log(M+1))$ initial sets, even though large coherent boxes
occur there. This excludes that stratum as a power-separation mechanism.
A stronger consecutive-block certificate now bounds the whole family
$I_{n,M,M}$ by $O(n^2\log^4(n+1))$, uniformly in $M$ and without
requiring $n=M^2$. At an active exchange, consecutive tight-quota blocks
form an ordinary lattice cut in a shifted wedge. In the fixed unimodular
coordinates $(a,b,b+c)$ its closed cut has box volume at most $512n$.
All interior facet normals therefore belong to one hyperbolic cross of
size $O(n\log^2(n+1))$; the classical plane-arrangement count gives
the stated bound. The large footprint of the entire standard set is
consistent with these small relative boxes for exchanges. Thus the
entire equal-valuation family is closed as a power-separation route.
The subsequent [cheap-prefix theorem](nested_prefix_bound.md) removes
the unequal-valuation obstacle and bounds the **entire** explicit nested
family by the same order. At a tied feasible exchange, lowering valuation
forces every cheaper monomial below the higher endpoint valuation into
the chosen basis. Its closed cut is an intersection of two positive
halfspaces with at most $8(n-1)$ lattice points. Its box in the original
coordinates has at most $216(n-1)$ points, bounding every active exchange
direction in one hyperbolic cross. This applies with the true length
$n=L-\max_k(k-\nu_m(k))$, including the non-binomial quotients of
embedding dimension three. In fixed dimension $d$, the bound is
$O_d(n^{d-1}\log^{(d-1)^2}(n+1))$ and one common list of that size
samples all the positive fans.

Generic normalized power series with any fixed tuple of leading orders
have the same rank rule on a nonempty Zariski-open coefficient set, so
the bound also covers generic one-parameter jets. Specializing off that
locus can change ranks and even the quotient length; it does not by
itself produce many cones. A two-cone all-field binomial example has a
quadratic primitive direction product and fails every positive valuation
lowering order. This preserves the distinction between a large normal and
a large fan. The whole additive nested family is closed as a power-gain
route. The nonlinear two-parameter rectangles remain a defined next
candidate. Classical matroid geometry, the Hilbert-zonotope sampling
method and arrangement counting are prior art; the exact prefix theorem
and generic-jet consequence have no novelty certificate and do not
complete the unrestricted goal.

## Connection to the original observations

Ordering smooth numbers generated by primes $p_1,\ldots,p_d$ orders their
exponent vectors $a\in\mathbb N^d$ by $a\cdot(\log p_1,\ldots,\log p_d)$.
A finite prefix is closed under taking coordinatewise smaller vectors. In a
polynomial ring, the same vectors index monomials, and such a set is a possible
set of standard monomials. Varying the positive weights therefore connects the
original interleavings to corner cuts and to monomial initial ideals.

Algebraic relations can change which monomial sets are bases of a quotient.
The question here is how much those dependencies can increase the number of
possible standard sets. This remains an exponent-lattice problem linked to the
original observations, but does not by itself give a prime-distribution result
or a faster smooth-number unranking algorithm.
The [five-cone family](initial_cone_height.md) gives a concrete distinction:
$x^2$ and $y^2$ are standard while $xy$ is dependent. A smooth prefix
cannot have this midpoint exclusion. The algebra changes which shapes
are possible; the number of shapes in one fan remains the unresolved
quantitative issue.

## Target and established comparison

For an ideal $I\subset K[x_1,\ldots,x_d]$ with

$$
\dim_K K[x_1,\ldots,x_d]/I=n,
$$

let $g(I)$ count its distinct monomial initial ideals under global monomial
orders. Let $G_d(n)$ be the maximum over such ideals and fields. Point ideals
give a subclass; write their maximum as $F_d(n)\le G_d(n)$.

| Comparison | Status and scope |
|---|---|
| $G_d(n)=O_d(n^{2d(d-1)/(d+1)})$, hence $G_3(n)=O(n^3)$ | Published: Babson--Onn--Thomas, Corollary 2.5, all finite quotients. |
| $G_2(n)=\Theta(n\log(n+1))$; every planar ideal has $g(I)\le1+4n\lceil\log_2n\rceil$ | [Complete dyadic proof](planar_dyadic_bound.md), all fields and all finite planar ideals, with the known generic corner-cut lower bound. The current-priority audit is incomplete. The three-variable endpoint remains open. |
| Generic corner cuts $C_d(n)=\Theta_d(n^{d-1}\log^{d-1}n)$ | The existing repository manuscript's result; its priority assessment is separate from the published upper bound above. |
| $G_3(n)=O(n^{44/15}\log^{4/5}(n+1))$ | Existing repository lattice-shadow candidate, with a remaining classical-priority comparison. This does not complete the present goal. |
| $G_3(n)=O(n^{23/8}\log^{3/8}(n+1))$; $\lvert\mathcal U(I)\rvert=O(n^{85/24}\log^{3/8}(n+1))$ | [Complete rational-shadow proof](rational_shadow_bound.md), all fields and finite quotients, with independently checked arithmetic projections and vertex visibility. Stronger than the older repository candidate; current priority and the sharp endpoint remain open. |
| Binomial ideals: $g(I)=O_d(n^{d-1}\log^{(d-1)^2}(n+1))$ | Existing repository reduction using published binomial normal forms and Minkowski-sum bounds. In three variables this excludes a polynomial separation from the generic scale within this subclass. |
| Intrinsic grading lattice of rank at most two: $g(I)=O_d(n\log(n+1))$; in three variables every nonzero integer grading is covered | [Exchange-lattice proof](profile_lattice_bound.md), all fields, uniform in grading weights. The maximum over this class has the matching planar lower order. Exact priority is unassessed; the full-rank unrestricted target remains open. |
| Additive nested residue matroids and generic fixed-order one-parameter jets: $g(I)=O_d(n^{d-1}\log^{(d-1)^2}(n+1))$ | [Uniform prefix theorem](nested_prefix_bound.md), independent of ambient truncation and valuations when expressed in the actual length. The explicit and generic-jet realizations use characteristic zero; the abstract matroid bound assumes no field. This closes that whole construction route, with exact priority unsettled. |
| Planar ideals containing an irreducible degree-$D$ equation: $g(I)\le5Dn$; sharp order $\Theta_D(n)$ for the maximum over degrees at most fixed $D\ge2$ | Uniform proof and exact wall controls below; the precise priority comparison remains open. This does not bound all ideals by a constant times $n$. |

The published upper bound and state-polyhedron framework are in
[Babson, Onn and Thomas, *The Hilbert zonotope and a polynomial time algorithm
for universal Gröbner bases*](https://arxiv.org/pdf/math/0207135).
The binomial argument and the unfinished lattice-shadow proof are preserved in
the [earlier research note](interleaving_research.md). Neither an arbitrary
representable matroid nor sharp examples for general lattice polytopes suffice:
the data must come from commuting multiplication by the variables on a cyclic
finite quotient.

The precise sharp $G_3(n)$ target is selected here. The inspected sources
establish the broader counting problem and bounds; they do not establish that
this exact formulation is a published open problem or that no later theorem
settles it. A complete current-priority audit remains part of the goal.

## Finite field source comparison

The full later primary version of Stigler and Zhang, *The Number of Gröbner
Bases in Finite Fields*, is available from the
[NSF repository](https://par.nsf.gov/servlets/purl/10177092)
([published chapter](https://doi.org/10.1007/978-3-030-42687-3_9)).
The inspected material includes Theorem 8 and its proof in Section 6,
Sections 2.1--2.2, the discussion and the Boolean table in the appendix.
The exact downloaded PDF has SHA-256
`a052c1b9c4934370b8186b67e671f656d6c811f6baa6072cc7d4d624fd796662`.
This resolves the earlier access gap for this later version.

Use $d$ for the number of variables and $m$ for the number of points, to
avoid confusing the chapter's notation with $G_d(n)$. For points in
$\mathbb F_p^d$, its Theorem 8 gives, when $0<m\le p^d/2$,

$$
O_d\!\left(\left(p^2\lfloor m/p\rfloor+(m\bmod p)^2\right)^{d(d-1)/(d+1)}\right).
$$

The complementary cardinality is used for larger $m$. The coordinate-sum
packing input, writing $m=pq+r$, is

$$
M=\frac{p(p-1)}2q+\frac{r(r-1)}2.
$$

Every standard set lies in the $p$-box. In each fiber parallel to a coordinate
axis it occupies an initial interval; filling intervals to length $p$, with at
most one residual interval, maximizes the coordinate sum. The formula is
attainable by downsets. The checker independently exhausts 1,194 downsets in
five small grids and verifies every coordinate maximum at every cardinality.

Two qualifications matter for the new target. If $p>m$, the bracket in
Theorem 8 is $m^2$, giving the earlier unrestricted exponent. Also these are
radical point ideals over a fixed finite field, not all finite quotients. In
using its volume argument for a state polytope of possibly smaller dimension,
one may replace $P$ by $P+[0,1]^d$ and $M$ by $M+1$. This supplies a
full-dimensional lattice polytope without losing vertices or changing the
asymptotic exponent. A literal zero ambient volume cannot bound the number of
vertices of a smaller-dimensional polytope.

The chapter's printed Theorem 4 states equality of the standard monomial sets
of complementary point sets. That literal assertion fails when their
cardinalities differ. The source it cites supplies the needed, different
statement: [Dimitrova, He, Robbiano and Stigler, *Small Gröbner Fans*, Theorem
4.8 and Corollary 4.11](https://arxiv.org/pdf/1704.03680v2).
For the full grid $B=\{0,\ldots,p-1\}^d$, the correct relation for each order is

$$
D(S^c)=B\setminus\{(p-1)\mathbf1-a:a\in D(S)\}.
$$

This bijection preserves equality of standard sets as the order varies, hence
preserves fans and their counts. It does not identify the two sets themselves.
Thus the complement step in the count has a valid justification. The checker
verifies the reflected formula for every subset of the Boolean three-cube in
every comparison chamber: 3,072 cases. The maxima at cardinalities 0 through 8
are $1,1,3,3,3,3,3,1,1$, agreeing with the later chapter's table.

The complete 1975 McMullen paper's lattice application remains an unresolved
priority comparison for the old lattice-shadow candidate. Its access issue
does not prevent studying realizable lower families for the present target.
No claim of a new classical geometric inequality follows from this iteration.

## A family with nonlinear coordinate relations

Consider

$$
A_{r,s,F}=K[x,y,z]/(x^r,y^s,z-F(x,y)),\qquad r,s\ge2,
$$

where $F$ has no constant term. Substitution identifies this algebra with

$$
B_{r,s}=K[x,y]/(x^r,y^s),
$$

so its length is exactly $rs$. The rectangular monomials $x^iy^j$,
$0\le i<r,0\le j<s$, are a reference basis. Multiplication by $x,y,F$
commutes and the class of 1 is cyclic. Changing $F$ changes the coordinate
embedding and its initial ideals, even though the abstract algebra is the same.

Unlike a binomial family or separate independent arms, a nonlinear $F$ gives
overlapping multi-term dependencies across a growing rectangle. This makes it
a candidate to analyze, not an established route to the target.

The fixed finite probe uses $K=\mathbb F_{101}$ and

$$
F=\sum_{\substack{0\le i<r,\ 0\le j<s\\i+j>0}}
\big((i+1)^2+3(j+1)+5ij\big)x^iy^j.
$$

| Rectangle | Length | Complete initial ideal count | Count for $F=x+y$ |
|---|---:|---:|---:|
| $2\times2$ | 4 | 6 | 5 |
| $2\times3$ | 6 | 12 | 5 |
| $3\times3$ | 9 | 26 | 5 |
| $3\times4$ | 12 | 54 | 5 |

These four counts support no asymptotic fit or separation claim. Every count
includes all positive-weight cones, not a sample of weights. The method is
finite exact fan enumeration within the established state-polyhedron framework.

Every possible standard monomial belongs to

$$
V_n^3=\{a\in\mathbb N^3:(a_1+1)(a_2+1)(a_3+1)\le n\},
$$

since its distinct divisors must also be standard. The implementation computes
all residue columns on this ground set by truncated polynomial products and
independently checks them using powers of the multiplication matrices.
Greedy linear independence in increasing weight order gives a standard basis.
For each visited basis, its minimal excluded monomials are reduced in that
basis. Requiring each excluded monomial to have greater weight than every term
in its reduction gives its entire cone. These relations have the claimed
initial ideal because that monomial ideal and the quotient both have length
$n$. Strict inequalities make the relation set a Gröbner basis there.

The checker walks across all interior facets using rational weights, avoiding
other comparison walls. Distinct cones have disjoint interiors; their exact
areas cover the normalized positive weight simplex. This is also a completeness
check independent of the walk's stopping condition. Any global monomial order
can be represented by a positive weight on the finitely many comparisons of
its Gröbner basis, so no additional monomial initial ideals are omitted.
For the length-4 case a separate permutation-formula determinant calculation
tests all 70 candidate bases and compares minima over every comparison chamber.
Every basis, strict integer weight and rational cone polygon is saved in the
[exact report](../results/research/extremal_initial_ideals_checks.json).

## The homogeneous control has at most five initial ideals

For $F=x+y$, the ideal is homogeneous. Its single linear relation eliminates
the largest variable, so any initial ideal contains exactly one of $x,y,z$.
If it contains $z$, the remaining ideal is $(x^r,y^s)$, giving one initial
ideal. If it contains $x$, substitution gives the homogeneous two-variable
ideal $((z-y)^r,y^s)$. A homogeneous ideal in two variables has at most two
initial ideals: on each homogeneous component the monomial comparisons depend
only on which variable is larger. The case containing $y$ gives at most two
more. Hence the total is at most five over every field.

If the characteristic is zero or greater than $r+s-2$, the coefficient of
$x^{r-1}y^{s-1}$ in $(x+y)^{r+s-2}$ is nonzero, so the nilpotence index of
$z=x+y$ is $r+s-1$. In the case eliminating $x$, the order with $z>y$
has $z^r$ in its initial ideal, whereas the order with $y>z$ has its first
pure $z$ power at $r+s-1$. These are distinct. The case eliminating $y$
is symmetric. The three groups have different degree-one initial monomials,
so the count is exactly five. Small characteristic can merge cones: the
length-4 control in characteristic 2 has three. This elementary control is
not proposed as a new research theorem.

## Support matching fails uniformly

A tempting shortcut would replace exact ranks of residue columns by maximum
matchings in their nonzero supports. It fails even within this realizable
family. The following identity explains a failure at every square size.

In $B_q=K[x,y]/(x^q,y^q)$, $q\ge3$, write

$$
F=ax+by+c_{20}x^2+c_{11}xy+c_{02}y^2+\text{higher terms},\quad ab\ne0.
$$

Set

$$
H=\binom{2q-3}{q-1},\quad
C=b^2c_{20}+ab c_{11}+a^2c_{02},\quad
U=F^{2q-3},\quad W=(xy)^{q-2}F,\quad V=F^{2q-2}.
$$

Then the following polynomial identity holds over every field:

$$
2a^2b^2U-2H(ab)^q W-(q-2)CV=0. \tag{1}
$$

To prove it, put $X=x^{q-1}y^{q-2}$, $Y=x^{q-2}y^{q-1}$, and
$T=x^{q-1}y^{q-1}$. Truncation leaves only these terms in $U$, only
$X,Y,T$ in $W$, and only $T$ in $V$. Direct multinomial expansion gives

$$
\begin{aligned}
W&=aX+bY+c_{11}T,\\
V&=2H(ab)^{q-1}T,\\
U-H(ab)^{q-2}W&=(q-2)H(ab)^{q-3}CT.
\end{aligned}
$$

Multiplication gives (1), with integer coefficients before reduction to $K$.
If the characteristic does not divide $2Hab$, $U$ and $V$ are independent
by degree, while (1) makes the rank of $U,W,V$ exactly two. Nevertheless their
support graph has a matching of size three: assign $U$ to $X$, $W$ to
$Y$, and $V$ to $T$. This includes characteristic zero.

There is a ground-set issue to address: the monomial defining $W$,
$x^{q-2}y^{q-2}z$, has $2(q-1)^2$ divisors, which exceeds the local length
$q^2$ for $q\ge4$. For a counterexample entirely inside the universal
ground set, work over $\mathbb Q$ and put

$$
m=\max(0,2(q-1)^2-q^2),\qquad
I_q=(x^q,y^q,z-F)\cap I(\{(j,0,0):1\le j\le m\}).
$$

The local component is supported at the origin, disjoint from the added points.
The Chinese remainder theorem gives length
$N=q^2+m=\max(q^2,2(q-1)^2)$. All three monomials contain $z$, so their
columns vanish on every added point and retain the same rank and relation.
Their divisor counts are $2q-2$, $2(q-1)^2$, and $2q-1$, respectively;
all are at most $N$.

For an explicit cyclic basis, take the rectangular monomials together with
$x^q,\ldots,x^{q+m-1}$. The extra columns vanish on the local component.
Their evaluation block on the added points is $(j^{q+k})$, a nonsingular
Vandermonde matrix after factoring nonzero $j^q$ from each row. Thus the
quotient has the stated dimension and is a genuine three-variable example.

The checker replays (1) by integer convolution for $q=3,\ldots,8$, checks
rank two independently modulo 101, and verifies the padding blocks. The
resulting lengths are 9, 18, 32, 50, 72, and 98. At $q=3$, for example,
the relation coefficients for $U,W,V$ are $4802,-705894,-1813$.
The report saves all three sparse integer columns.

This rejects a support-only rank characterization even on the admissible
ground set. It does not show these three columns are simultaneously standard,
give a fan lower bound, or resolve the extremal counting target. No novelty
claim is attached to this diagnostic identity.

## A uniform bound from the local Hilbert function

Let $A=K[x_1,\ldots,x_d]/I$ have length $n$, with radical
$(x_1,\ldots,x_d)$. Thus it is local and supported at the origin, over any
field. Write

$$
h_j=\dim_K\mathfrak m^j/\mathfrak m^{j+1},\qquad
S=\sum_{j\ge0}j h_j,\qquad \mathfrak m=(x_1,\ldots,x_d)A.
$$

For **every monomial vector-space basis** $D$ of $A$, not only initial
standard bases,

$$
\sum_{a\in D}|a|\le S. \tag{2}
$$

Choose a basis adapted to the descending filtration by powers of
$\mathfrak m$, assigning grade $j$ to the $h_j$ basis vectors in that
quotient. A monomial of total degree $k$ lies in $\mathfrak m^k$, so its
column has zero entries in rows of grade below $k$. A nonzero determinant of
the columns indexed by $D$ has a nonzero permutation term. Each matched row
therefore has grade at least the degree of its column. Summing proves (2).
No graded homogeneity of $I$ is assumed.

Every initial standard basis gives a distinct exposed moment vertex, so the
state-polyhedron reduction and the classical lattice vertex bound used in
[Babson--Onn--Thomas, Corollary 2.5 and its proof](https://arxiv.org/pdf/math/0207135)
give

$$
g(I)=O_d\big((1+S)^{d(d-1)/(d+1)}\big). \tag{3}
$$

Indeed, the basis moments lie in the nonnegative simplex with total coordinate
sum at most $S$. Add the unit cube to their convex hull if necessary to
handle smaller affine dimension. Its volume is $O_d((1+S)^d)$; the Minkowski
sum preserves at least the original number of vertices.

The proof covers a local algebra supported at any $K$-rational point after
translation. Translation preserves leading monomials under every global
monomial order: each translated monomial expands into itself and its proper
divisors. The leading coefficient is unchanged. This does not extend (2) to
an algebra supported at several different points.

For three variables, a family with $S=O(n^{4/3})$ has $g(I)=O(n^2)$.
Consequently a local family with $g(I)\ge c n^{2+\epsilon}$ needs
$S=\Omega(n^{4/3+2\epsilon/3})$. If $\mathfrak m^L=0$, then
$S\le n(L-1)$, so it also needs
$L=\Omega(n^{1/3+2\epsilon/3})$. These are necessary conditions, not lower
constructions. In particular, short nilpotent filtrations cannot provide the
requested polynomial separation. A long filtration does not ensure one.

The determinant argument is elementary and (3) applies established geometry.
This is a screening bound, with no novelty or sharpness claim.

## The quantitative window for nonlinear rectangles

For $A_{r,s,F}$, the local Hilbert function is that of
$K[x,y]/(x^r,y^s)$, independently of $F$. In particular,

$$
S=\frac{rs(r+s-2)}2.
$$

There are stronger separate coordinate bounds. A residue of $x^a y^b z^c$
in the rectangular reference basis only uses cells $(i,j)$ with
$i\ge a$, $j\ge b$, and $i+j\ge a+b+c$. Applying the same determinant
matching gives, for every monomial basis moment $\mu$,

$$
0\le\mu_x\le A_0=\frac{rs(r-1)}2,\quad
0\le\mu_y\le B_0=\frac{rs(s-1)}2,\quad
\mu_x+\mu_y+\mu_z\le A_0+B_0.
$$

The containing box, enlarged by the unit cube, has volume
$O(A_0 B_0(A_0+B_0))$. The three-dimensional lattice vertex bound yields
the uniform all-field estimate

$$
g(I)=O\!\left(n^{3/2}
       \sqrt{(r-1)(s-1)(r+s-2)}\right)
     =O\!\left(n^2\sqrt{r+s}\right),\qquad n=rs. \tag{4}
$$

For $r=s=q$, (4) is $O(n^{9/4})$. This narrows the possible gain over
the generic $n^2\log^2 n$ scale to at most a quarter in the power of $n$
within square rectangles. It does not rule out such a gain or prove one.
For $F=x+y$, the stronger five-cone argument above applies.

The exact checker now verifies the three moment budgets for every saved
standard basis and for every independent determinant basis in the length-4
exhaustive control, including bases that are not downsets. This tests the
scope of the determinant lemma rather than just the displayed fan counts.

## The fixed probe lies on a special rational surface

The coefficients chosen for the finite probe are not algebraically
independent. Over the integer formal power series ring,

$$
\begin{aligned}
D&=(1-x)^3(1-y)^2,\\
N&=(1+x)(1-y)+3(1-x)^2+5xy(1-x),\\
F+4&=N/D.
\end{aligned}
$$

The last equality follows by summing the three elementary series for
$(i+1)^2$, $j+1$, and $ij$. Reduction modulo $(x^r,y^s)$ gives the
chosen polynomial $F$ at every rectangle size, in every characteristic.
Therefore every fixed-probe ideal contains the same polynomial

$$
h=D(z+4)-N. \tag{5}
$$

Its unique coordinatewise largest monomial is $x^3y^2z$, with coefficient
$-1$. Thus (5) has that leading monomial for every global monomial order.
The fixed probe lies on a single surface with an order-independent initial
monomial, even as its quotient length grows.

This specialization changes actual rank on the admissible ground set. Over
the integers, $D$ has 12 terms and $4D-N$ has 11, the constant term being
zero. Let $E$ be the set of these 23 exponent vectors in $h$. For every
$r\ge8,s\ge6$, all of $E$ belong to $V_{rs}^3$, since the largest
divisor count among them is 24. Equation (5) makes their rank at most 22 in
the fixed probe, over every field.

For contrast, set

$$
F_*=x+y+t x^4y^3,\qquad t\ne0.
$$

Use the same 23 monomial columns. Select the 11 rows in the support of
$4D-N$, followed by the 12 rows $(4+i,3+j)$, $0\le i\le3,0\le j\le2$.
The selected matrix, with constant-in-$z$ columns first, is

$$
\begin{pmatrix}
\operatorname{Id}_{11}&B\\
0&t\operatorname{Id}_{12}
\end{pmatrix}.
$$

The lower shifted rectangle is disjoint from all low rows and from the
linear contributions of $F_*$. Its determinant is $t^{12}$. Thus these
23 columns have full rank over every field when $t\ne0$. The corresponding
minor is consequently a nonzero polynomial in the coefficients of a general
$F$: over their rational function field the rank is 23.

The checker verifies (5) by integer convolution, checks the unit triangular
minor entrywise at $t=1$, and independently computes ranks modulo 101 at
rectangles $8\times6$, $9\times7$, and $12\times8$. They are respectively
22 and 23 in each case. This is a uniform rank distinction, supported by exact
finite witnesses, rather than a rank guess from a random specialization.

The earlier counts remain valid for their stated coefficient choice. This
distinction does not imply that a generic embedding has more cones: enlarging
the set of feasible bases need not increase the number of vertices of their
moment hull. The relation alone also does not exclude a large fan for the
fixed surface. An asymptotic argument must specify which coefficient class it
uses and prove its coherent-basis count.

### A coherent distinction between the embeddings

The difference is visible in actual initial ideals, as well as in arbitrary
column subsets. Let $a,b\ge1$, $t\ne0$, and set

$$
I_*=(x^{2a},y^{2b},z-x-y-tx^ay^b).
$$

Its length is $4ab$. Write $\alpha=z-x-y$. Besides the two pure powers,
the following four polynomials belong to $I_*$:

$$
x^ay^b-t^{-1}\alpha,\qquad
x^a\alpha,\qquad y^b\alpha,\qquad \alpha^2.
$$

Indeed $\alpha=tx^ay^b$, so the last three vanish by the pure powers.
For every positive weight $(X,Y,Z)$ satisfying

$$
\max(X,Y)<Z<aX+bY,
$$

their respective leading monomials are $x^ay^b,x^az,y^bz,z^2$.
Together with $x^{2a},y^{2b}$ they generate a monomial ideal $J$.
Its standard set is

$$
\begin{aligned}
D_*={}&\{x^iy^j:0\le i<2a,\ 0\le j<2b,\ i<a\text{ or }j<b\}\\
 &{}\cup\{x^iy^jz:0\le i<a,\ 0\le j<b\}.
\end{aligned}
$$

It has $3ab+ab=4ab$ elements. Since $J\subseteq\operatorname{in}(I_*)$
and both quotients have the same length, equality follows. This proves a
Gröbner basis throughout the displayed cone, in every field. The integer
weight $(2,2,3)$ lies in it for every $a,b\ge1$.

If $a\ge4,b\ge3$, $D_*$ contains $x^3y^2z$. Equation (5) excludes
that monomial from every initial standard set of the fixed rational probe at
the same rectangle $r=2a,s=2b$. Thus these two embeddings of the same abstract
algebra have different collections of initial ideals. The checker verifies
the six relations, their strict leading inequalities, and the entire selected
basis by independent greedy residue ranks at $(a,b)=(4,3),(5,3),(4,4)$.
It does not enumerate the full fans of these larger examples.

This is a uniform diagnostic distinction, not a comparison of total fan sizes
or a claim that either family achieves polynomial amplification.

## Generic graphs over torus grids reduce to ordinary cuts

There is a broad comparison family where exact rank is tractable. It is
different from the local nilpotent rectangles above.

Let $G=C_r\times C_s$, $n=rs$. Work over a field $K$ of characteristic
not dividing $n$, containing the required roots of unity. Index the torus
grid by $g\in G$, and let $\chi_{a,b}(g)$ be its characters,
$0\le a<r,0\le b<s$. Give each grid point an independent indeterminate
$\lambda_g$. The $n$ points

$$
(\zeta_r^i,\zeta_s^j,\lambda_{i,j})
$$

are rational over $L=K(\lambda_g:g\in G)$. Their algebra has length $n$,
with $x^r=y^s=1$, and $z=F(x,y)$ for the unique interpolating polynomial
of separate degrees below $r,s$.

**Generic rank claim.** Every $n$ distinct columns

$$
v_{a,b,c}(g)=\chi_{a,b}(g)\lambda_g^c,\qquad
0\le a<r,\quad0\le b<s,\quad c\ge0,
$$

are independent over $L$. This holds for composite $r,s$ as well.

Here is a proof using established matroid theory. Group the chosen columns by
their exponent $c$, with $k_c$ columns in a block. Their character columns
are independent because the full character table is invertible. The row
matroid of that block has rank $k_c$ and is transitive under translation
of $G$: translating rows merely rescales character columns. Translating
any one of its row bases uniformly over $G$ gives a distribution in which
each row belongs to the basis with probability $k_c/n$. Intersecting these
bases with a row subset $T$ proves

$$
\operatorname{rank}_c(T)\ge k_c|T|/n.
$$

Thus $\sum_c\operatorname{rank}_c(T)\ge |T|$ for every $T$.
The classical matroid partition theorem supplies disjoint row bases $R_c$
covering all $n$ rows. One primary statement is Theorem 2 in
[Cunningham, *The Coming of the Matroids*](https://www.math.uwaterloo.ca/~whcunnin/coming.of.matroids.pdf);
the original result is Edmonds's matroid partition theorem.
The density condition is established terminology, not a new invariant:
see [Devriendt--Mulas, Theorem 1.1 and its proof, published 2026](https://ajc.maths.uq.edu.au/pdf/94/ajc_v94_p145.pdf).

Expand the selected determinant as a polynomial in the independent row
parameters. The coefficient of $\prod_c\prod_{g\in R_c}\lambda_g^c$
is, up to sign, the product of the nonzero character minors on $R_c$.
Different row partitions give different parameter monomials, so this
coefficient cannot cancel with another partition. The determinant is
nonzero, proving the claim. This argument does not assume that every smaller
minor of a composite Fourier matrix is nonzero.

Every monomial outside

$$
S_{r,s}=\{(a,b,c):0\le a<r,\ 0\le b<s,\ c\ge0\}
$$

has the same residue as a cheaper monomial obtained by reducing $a,b$
modulo $r,s$. The first $n$ monomials of $S_{r,s}$ in a positive
weight order are independent by the claim. Hence they are exactly the
standard set. The entire initial-ideal count is therefore the number of
ordinary $n$-element positive halfspace cuts of $S_{r,s}$.

Only $E=S_{r,s}\cap V_n^3$ is needed: the divisors of any selected monomial
are also selected. Thus $n\le |E|\le n(1+\log n)^2$. For a particular
length, the finitely many nonzero determinant conditions can also be
specialized over an infinite field. The proof does not guarantee a simultaneous
scalar specialization in a given small finite field.

The checker independently expands all 495 tested maximal determinant
polynomials for $C_2\times C_2$ with three $z$-blocks, over the integers,
and all 924 for $C_2\times C_3$ with two blocks in characteristic 7.
Both underlying Fourier matrices have zero $2\times2$ minors. Sample exact
coefficient witnesses and digests of the full witness ledgers are recorded. These are
applications of classical partition theory; priority of the combined
formulation is unestablished.

## Why the generic cut route is a demanding lower bound target

More generally, suppose an algebraic family has its standard sets equal to
ordinary positive-weight $n$-cuts of a configuration $E\subset V_n^3$.
A polynomial gain here would yield a polynomial gain for the classical
three-dimensional $k$-set problem, as follows.

Write $N=|E|$. If $N\ge2n$, add $N-2n$ points with all coordinates
negative. Every positive cut selecting $n$ original points selects these
extra points too. The new configuration has $M=2(N-n)$ points, and the
selected sets have size $M/2$. If $N<2n$, add $2n-N$ points with all
coordinates greater than $n$; they are outside every relevant positive
cut, whose original points have coordinates at most $n-1$. Now $M=2n$
and the selected sets again have size $M/2$.

In both cases $N\le M\le2N$ and distinct cuts remain distinct. Choose one
strict separating certificate for each of the finitely many cuts and perturb
the points sufficiently little to obtain general position while preserving
every certificate. This is the standard strict-separation perturbation,
also given in [Wagner's thesis, Observation 2.4](https://ti.inf.ethz.ch/ew/research/phdtheses/wagner03.pdf).

Consequently, $n^{2+\epsilon}$ such algebraic standard sets, with
$N=O(n\log^2 n)$, would yield $M^{2+\epsilon-o(1)}$ balanced strict
halfspace subsets of an $M$-point configuration in three dimensions.
Wagner's equation (2.7) records the lower construction
$N k^{d-2}\exp(\Omega(\sqrt{\log k}))$. The later
[Kovács--Tóth paper, introduction and Theorem 2](https://arxiv.org/pdf/1704.00229v2)
constructs the corresponding rate for dense halving configurations. These
benchmarks have only a subpolynomial factor above the quadratic scale in
dimension three.

This is a reason to treat a uniform-cut construction as a difficult geometric
research route, not a reason to regard the full algebraic goal as settled or
equivalent to that problem. The broader current-priority audit is still open.
A [September 2026 proposed planar upper-bound improvement](https://math20.org/f/183/)
surfaced during the search, but its full proof was not accessible for inspection.
No improved bound or transfer to general parametric matroid projections is assumed.

## Coherent rank defects persist throughout the local family

Removing zero monomials does not make the local rectangle family a uniform
cut problem. There is an infinite obstruction inside actual cheapest
prefixes, stronger than a circuit on an arbitrary column subset.

Let $q\ge3$, $F=\alpha x+\beta y+\text{higher terms}$, with
$\alpha\beta\ne0$, in $K[x,y]/(x^q,y^q)$. Assume characteristic zero,
or characteristic greater than $2q-2$. The nonzero monomials in the
$y,z$ face are exactly

$$
y^b z^c,\qquad 0\le b<q,\quad0\le c\le2q-2-b.
$$

For the forward direction, the lowest homogeneous part
$y^b(\alpha x+\beta y)^c$ has a surviving nonzero term. In the other
direction, total degree exceeds the last nonzero degree of the rectangle.

At weights $(X,Y,Z)=(4q^2,2q-1,1)$, all of one $y$-row precedes the
next, and the first $q^2$ nonzero monomials occur before any monomial
containing $x$. The first two rows contain $4q-3\le q^2$ terms, so the
prefix includes both $z^{2q-2}$ and $y z^{2q-3}$. Their residues lie
in the one-dimensional socle and satisfy

$$
z^{2q-2}=2\beta\, y z^{2q-3}\ne0.
$$

Indeed their coefficients on $x^{q-1}y^{q-1}$ are
$\binom{2q-2}{q-1}\alpha^{q-1}\beta^{q-1}$ and
$\binom{2q-3}{q-1}\alpha^{q-1}\beta^{q-2}$, respectively. Thus that
first-$n$ nonzero prefix has rank below $n=q^2$.

There is nevertheless a full $y,z$ basis before the first $x$ term.
To see this, put $f(y)=F(0,y)$. Since
$F-f(y)=xU(x,y)$ with $U(0,0)=\alpha\ne0$, the finite implicit
coordinate change gives

$$
K[y,z]/(y^q,(z-f(y))^q)\simeq K[x,y]/(x^q,y^q).
$$

The monomials $y^i z^j$, $0\le i,j<q$, are a basis and have weight
at most $2q(q-1)<X$. The genuine greedy standard set must therefore
reject at least one of the first $n$ nonzero terms. Small generic weight
perturbations preserve this failure and give strict coherent witnesses.

The checker verifies prefixes of lengths 9, 16, 25, 36, 49, and 64; their
ranks are respectively 8, 15, 22, 32, 43, and 55 for the fixed rule modulo 101.
The theorem is uniform in $F$ with the stated hypotheses. It is not a
count amplification theorem: even $F=x+y$ has this defect, while its
entire fan has only five cones.

## Existing fans already contain interior dependent prefixes

For each saved state cone, compare all its polygon vertices with every
inequality saying that each basis monomial is cheaper than every excluded
nonzero monomial. If all hold, the cone equals an ordinary-cut cone.
A violated inequality supplies a strict interior witness where the first
$n$ nonzero terms are dependent. This checks whole cones, rather than just
the originally saved sample weight.

| Length | Existing complete fan count | Cones containing a nonuniform region | Such cones with $x,y,z$ all standard |
|---|---:|---:|---:|
| 4 | 6 | 0 | 0 |
| 6 | 12 | 0 | 0 |
| 9 | 26 | 4 | 2 |
| 12 | 54 | 6 | 4 |

For example, the length-9 counterexamples include an initial standard set
containing all three generators, while a cheaper nonzero monomial is
dependent on the earlier columns. The exact weights, bases and prefix
vectors are in the report. Thus dependent regions cannot be assigned only
to coordinate-face bases. These counts do not say that whole nonuniform
cones have no ordinary-cut weights, or that the total fan is larger than
the ordinary-cut fan. No larger fan enumeration is needed to reject the
face-confinement claim.

## Published neighbor constraints and their quantitative limit

[Hering--Maclagan, *The T-graph of a multigraded Hilbert scheme*](https://arxiv.org/pdf/1110.1861),
Definition 1.2 and Theorem 1.3, give necessary conditions for a pair of
monomial ideals to be endpoints of one torus orbit. For a positive grading
$\mathbb Z^d/\mathbb Zc$, an arrow map is a degree-preserving bijection
between their ideal monomials. Its displacement along $c$ is nonnegative
in one term order and cannot increase on multiples, in either direction.
There must also be such a map for the two colon ideals by a common box
$Q=(x_1^{r_1},\ldots,x_d^{r_d})$. Proposition 3.2 constructs the map from
the intervening ideal. Example 3.7 shows that the first map alone is
insufficient. These statements concern one pair, potentially using a
different intervening ideal for each pair.

They apply to neighboring cones of one finite quotient: take a generic
positive weight in the common wall. Its initial ideal has a positive
$\mathbb Z^d/\mathbb Zc$ grading, where $c$ is the primitive wall normal.
The two small opposite perturbations give the two neighboring monomial
ideals. There are only finitely many relevant Gröbner-basis comparisons;
the wall weight can avoid every nonparallel comparison. This supplies
the required one-dimensional torus orbit after a field extension, if needed.
No converse or bound for the number of simultaneous neighbors follows.

The later [Maclagan--Silversmith paper, published in 2024](https://escholarship.org/content/qt74v009bx/qt74v009bx.pdf)
proves that the spine edges exist over every infinite field (Theorem 1.2).
Its generic uniform-matroid theorem, Theorem 1.3, requires the **standard**
grading in two variables; Section 3.5 gives counterexamples for other
positive gradings. In particular, Example 3.25 has

$$
I_t=(x^7,xy,y^4+t x^6),\qquad \deg(x)=2,\quad\deg(y)=3.
$$

For $t\ne0$, the degree-12 quotient has rank one on
$x^6,x^3y^2,y^4$, with $x^3y^2$ a loop. It is not uniform. Their
Example 3.26 goes further: in length 29, a degree-18 rank-two matroid
on four monomials is nonuniform with no loops or coloops. For parameters
$(1,2,1)$, quotient columns can be chosen as
$(1,0),(0,1),(0,-1),(-1,-1)$, in every characteristic.
Thus loop removal does not repair the general assertion. This matters
when using [Silversmith's earlier preprint](https://arxiv.org/pdf/1911.03569):
its Remark 3.11 announces a forthcoming assertion for arbitrary positive
two-variable gradings. That broader assertion is not the later theorem.
The earlier paper's Theorem 5.11 instead classifies a particular family
of edges by binary necklaces; it is not a count of one three-variable fan.
Neither paper inspected here settles $G_3(n)$.

### A uniform obstruction to counting individually realizable exchanges

The difference between individual and simultaneous realizability can be
exponential. The following elementary construction makes that obstruction
explicit; it is a control from the classical state-polyhedron framework,
not a new asymptotic extremal result.

Fix $q\ge2$ and $k=\lfloor(q+1)/2\rfloor$. Let

$$
L_q=\{(i,j,0):i+j<q\},\qquad
E_q=\{(i,q-i,0):0\le i\le q\}.
$$

For every $k$-subset $S\subset E_q$, put $D_S=L_q\cup S$. These
are downsets of the common cardinality

$$
n_q=\frac{q(q+1)}2+k.
$$

Connect $S$ and $T$ when they differ by replacing one element. This is
the connected Johnson graph $J(q+1,k)$, with

$$
\binom{q+1}{k}\ \ge\ \frac{2^{q+1}}{q+2}
       =\exp(\Omega(\sqrt{n_q}))
$$

vertices. **Every edge is realizable over every field.** To prove this,
write $S\setminus T=\{a\}$, $T\setminus S=\{b\}$, and let $C$ be the
monomial ideal whose standard set is $D_S\cup D_T$. Both $x^a$ and
$x^b$ are in the socle of $K[x,y,z]/C$: all their positive-degree
multiples vanish. Consequently

$$
I_{S,T}=C+(x^a-x^b)
$$

has length $n_q$. If $w_x<w_y$, the cheaper of the two shell monomials
survives; if $w_x>w_y$, the other survives. These are exactly its two
standard sets $D_S,D_T$. Thus every graph edge comes from an actual
binomial ideal, not just a representable matroid or an arrow-map test.

**One fixed ideal can realize at most two of these vertices.** Indeed,
the state vector $s(D)=\sum_{a\in D}a$ of every $D_S$ satisfies

$$
s_z(D_S)=0,\qquad
s_x(D_S)+s_y(D_S)=\sum_{a\in L_q}|a|+kq.
$$

All these state vectors lie on a single affine line. For a fixed ideal,
an initial basis is the unique positive-weight minimum among all its
monomial vector-space bases. Among collinear vectors, only the two
endpoints of the convex hull can be unique minima. Distinct bases with
the same state vector cannot both be unique minima either. This proves
the bound of two even if the fixed ideal is nonhomogeneous and has other
initial bases outside the displayed family.

The checker exhausts all exchanges for $2\le q\le7$, replays the
commuting multiplication maps, and checks both initial bases by independent
modular elimination in characteristics 2, 3 and 101. It also checks the
published degree-12 and degree-18 nonuniform examples. The proof above,
rather than this finite range, gives the uniform statement.

Therefore the size or connectivity of a graph passing all pairwise
realizability tests cannot replace a count of the fan of one fixed ideal.
The example uses parallel walls and does not rule out using arrow maps
*together with* global state geometry. No larger local fan search or new
standalone arrow-map search is justified by the source comparison alone.

## An exact product construction reduces a lower bound to the plane

There is a uniform way to turn a planar counting improvement into the
requested three-variable separation. It constructs one quotient realizing
all combinations, rather than collecting pairs from different ideals.

Let $J_i\subset K[x,t]$, $1\le i\le r$, have finite positive colength
$m_i$. Suppose a common monic $f(x)\in K[x]$ of degree $q$ belongs
to every $J_i$. In $K[x,y_1,\ldots,y_r]$ put

$$
H=\left(f(x),\ y_i y_j\ (i<j),\ y_i J_i(x,y_i)\ (1\le i\le r)\right).
\tag{5}
$$

Then, over **every field**,

$$
\operatorname{length}(H)=q+\sum_i m_i,\qquad
g(H)=\prod_i g(J_i). \tag{6}
$$

More precisely, for every global monomial order $\prec$,

$$
\operatorname{in}_\prec H=
\left(x^q,\ y_i y_j\ (i<j),\
       y_i\operatorname{in}_{\prec_i}J_i(x,y_i)\ (1\le i\le r)\right),
\tag{7}
$$

where $\prec_i$ is the restriction to $K[x,y_i]$.

**Proof of length and (7).** Set $C=K[x]/(f)$. Modulo $f$ and
the cross products, each polynomial has a unique expression as a constant
in $C$ plus terms belonging separately to $y_i C[y_i]$. Because
$f\in J_i$, imposing the remaining relations replaces the $i$-th
summand by a vector-space copy of $y_i(K[x,t]/J_i)$. Hence

$$
K[x,y_1,\ldots,y_r]/H\simeq
C\ \oplus\ \bigoplus_i y_i(K[x,t]/J_i)
\quad\text{as vector spaces}.
\tag{8}
$$

This proves the length, without assuming that $x$ generates any branch.
For a chosen order let $D_i$ be its branch standard set. The right side
of (7) is contained in $\operatorname{in}_\prec H$, since multiplication
by $y_i$ preserves leading monomials. Its standard set is

$$
\{1,x,\ldots,x^{q-1}\}\ \sqcup\
\bigsqcup_i y_i D_i(x,y_i).
\tag{9}
$$

Its cardinality is $q+\sum_i m_i$, equal to the quotient length.
Inclusion of the two finite-colength monomial ideals is therefore equality.

It remains to prove that all branch combinations occur. Choose a positive
rational weight $(a_i,b_i)$ realizing each specified planar initial ideal;
only its finite Gröbner-basis comparisons need to be represented. Normalize
its first coordinate to 1. The single positive weight

$$
(w_x,w_{y_1},\ldots,w_{y_r})=(1,b_1/a_1,\ldots,b_r/a_r)
$$

realizes every branch simultaneously. Different branch tuples give different
sets (9). This proves (6), as well as a product description of the positive
fan in the independent ratios $w_{y_i}/w_x$.

For any length-$m$ planar quotient, the characteristic polynomial of
multiplication by $x$ supplies a monic $f\in J$ of degree exactly $m$.
Take $r=d-1$ identical branches in (5). Thus the unrestricted extremal
functions satisfy

$$
G_d(dm)\ \ge\ G_2(m)^{d-1},\qquad
G_3(3m)\ \ge\ G_2(m)^2. \tag{10}
$$

The maximum can be taken over fields, since the construction uses the same
field as each chosen planar quotient. In particular, a realizable planar
family with $g(J_m)\ge c m^{1+\epsilon}$ would yield
$g(H_m)\ge c'\operatorname{length}(H_m)^{2+2\epsilon}$ in three
variables. Conversely, a bound $G_3(n)\le n^{2+o(1)}$ would force
$G_2(m)\le m^{1+o(1)}$. This is a quantitative link, **not** a proof
of either premise. Applying (10) only to the generic planar scale
$m\log m$ recovers the generic three-variable order, not a power gain.

There is also an explicit universal-basis consequence. If $U_i$ is a
universal Gröbner basis of $J_i$, then

$$
\{f\}\ \cup\ \{y_i y_j:i<j\}\ \cup\ \bigcup_i y_i U_i
\tag{11}
$$

is a universal Gröbner basis of $H$: for each order its leading monomials
generate (7). Its size is at most
$1+\binom r2+\sum_i|U_i|$. Thus the number of initial ideals multiplies
while this explicit universal generating set grows additively. No new
algorithmic operation bound is inferred.

The transfer also holds for arbitrary reduced rational point sets; a
separating first coordinate is unnecessary. Pass to an infinite field
extension if needed. Field extension preserves every monomial initial ideal.
Translate the second coordinate by a constant so that every $b_j\ne0$.
Translation also preserves all global initial ideals: a translated monomial
is itself plus proper divisors, so every polynomial keeps its leading
monomial.

Choose $m$ distinct elements $A'\subset K$ containing every first
coordinate $a_j$, padding with extra elements if some first coordinates
repeat. Now take $f=\prod_{a\in A'}(x-a)$. It is squarefree, has degree
$m$, and belongs to the translated point ideal. The spatial construction
is the ideal of the $3m$ distinct points

$$
(a,0,0)\ (a\in A'),\quad(a_j,b_j,0),\quad(a_j,0,b_j).
$$

Each generator vanishes there, and equality follows from the proved length.
Thus the same all-fields extremal comparison holds for the point subclass:

$$
F_d(dm)\ge F_2(m)^{d-1},\qquad F_3(3m)\ge F_2(m)^2.
\tag{12}
$$

Over any fixed infinite field the point construction stays in that field.
No smoothing of a nonreduced quotient is used for this assertion.

### Verification and the remaining novelty boundary

The checker uses ten planar fixtures over characteristics 2, 3 and 101,
including a shared coordinate that is not cyclic. For every combination
of their complete planar fans, it compares full-ground selection by
commuting quotient matrices with polynomial S-pair reduction for the lifted
generators. The spatial lengths are 9 and 12. Reduced nine- and twelve-point
cases are checked independently by point evaluation. The latter has repeated
first coordinates in its planar input and uses an added axis root; it also
checks translation invariance across every planar comparison interval.
The displayed product counts
are checks of the uniform proof, not a fitted asymptotic family.

Fibre-product presentations and initial ideals have substantial prior art.
[Conca--Römer's preprint, published in 2010](https://arxiv.org/pdf/0707.0678),
Section 1 and Theorem 1.2, studies standard graded fibre products over
$K$ in disjoint variable sets and generic initial ideals for reverse lexicographic
order in characteristic zero. Here the augmented branches are glued over
$K[x]/(f)$, and (7) concerns all orders in fixed coordinates. The inspected
result does not state (10); this comparison does not establish novelty of
the elementary construction or its counting consequence. The precise
fibre-product priority comparison remains open.

This initially suggested focusing on a planar power improvement while
retaining the three-variable endpoint. The later
[dyadic planar theorem](planar_dyadic_bound.md) now rules out such a
power improvement entirely. The transfer remains valid, but planar
branches cannot yield the intended spatial power separation. The
following proposed realization is retained to make the obstruction precise.

There is a specific existing lower construction to test, rather than an
arbitrary search for planar examples. [Eppstein, Theorems 5 and 7 and
Corollaries 1--2](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf),
uses transversal matroids from convex chains to obtain
$\Omega(M^{4/3})$ parametric bases with $\Theta(M)$ elements and rank.
The earlier note already distinguished this from a quotient-algebra
construction. Equations (5)--(10) now supply the spatial product **if**
a planar realization with controlled cost can be found.

Specifically, a realization of those weighted matroids by a length
$L=O(M^\alpha)$ planar quotient, preserving the witnessed exponent-weight
comparisons, would give a spatial family with at least
$\Omega(L^{8/(3\alpha)})$ initial ideals. To beat the quadratic power
requires $\alpha<4/3$; linear overhead would give power $8/3$.
Neither representation by arbitrary vectors nor a large abstract matroid
fan supplies this: multiplication must commute on a cyclic quotient,
and the monomial exponent sizes must fit the same colength budget.
The later planar theorem supplies that quantitative obstruction:
$L=\Omega(M^{4/3}/\log(M+1))$ is necessary. In particular no fixed
$\alpha<4/3$ can work. Arbitrary representability or coefficient tuning
cannot remove this colength cost.

## A degree-sensitive bound for an irreducible plane carrier

Let $I\subset K[x,y]$ have colength $n\ge1$, and suppose
$F\in I$ is irreducible of total degree $D\ge1$. Over every field,

$$
g(I)\le 5Dn. \tag{13}
$$

The hypothesis is scheme-theoretic: the equation $F$ belongs to the ideal.
The argument allows nonreduced finite quotients and singular carriers.
It does not assume an irreducible low-degree equation exists for every $I$.

**The common variation budget.** Scan positive weights $(t,1)$ from
$t=0$ to $t=\infty$. Finite strict weight representations cover all
global monomial initial ideals. Across an actual wall $t=b/a$, with
$a,b$ relatively prime positive integers, the standard-set moment changes by
$(-Ka,Kb)$ for an integer $K>0$. Minimum basis cost is concave, so its
$x$ derivative decreases; equality of the two costs gives the $y$ change.
Each standard set is an $n$-element downset, whose coordinate moments are
between 0 and $\binom n2$: its coordinate fibers have lengths $\ell_j$
summing to $n$, and $\sum_j\binom{\ell_j}{2}\le\binom n2$.
Therefore, summing over **one fixed ideal's** fan,

$$
\sum_{\text{walls}}K(a+b)\le n(n-1). \tag{14}
$$

Also $a,b\le n-1$: a primitive wall direction is the difference direction
of exchanged standard monomials. This can be seen degree by degree in
the weighted initial ideal at the wall. Its two refinements have the same
number of standard monomials in each weighted degree.

At most $D$ wall ratios can be ties for the leading monomial of $F$.
Indeed, its positive-weight upper envelope has at most $D+1$ different
$x$ exponents. Count these walls separately. There are at most
$2(D-1)(n-1)$ remaining possible ratios with $a<D$ or $b<D$.

**A large wall has a short polynomial lift.** Consider another wall with
$a,b\ge D$, and let $x^r y^s$ be the unique leading monomial of $F$
there. Both neighboring standard sets avoid its multiples. A pure head
would put them in a strip of width at most $D$, incompatible with such a
direction. Hence $r,s>0$ and $r+s\le D$.

A line in direction $(a,-b)$ meets the monomials outside
$(x^r y^s)$ in at most two points: at most one in the vertical strip
$x<r$, and at most one in the horizontal strip $y<s$. A changed weighted
degree thus exchanges

$$
u=(i+\ell a,j),\qquad v=(i,j+\ell b),\qquad
0\le i<r,\quad0\le j<s,\quad\ell\ge1.
$$

As $t$ increases, the lower $x$ exponent wins. All changed weighted
degrees move in this direction, so $K$ is the sum of their $\ell$ values.
In particular, any exhibited exchange has $\ell\le K$.

Reduce $x^u$ in the standard basis immediately after the wall, and put
$g=x^u-\operatorname{NF}(x^u)\in I$. All its terms have wall weight
at most $T=u\cdot(b/a,1)$. At least one tail term ties: otherwise
$x^u$ would also be an initial monomial just before the wall. The
two-point property makes that term $x^v$. Every term of $g$ avoids
$x^r y^s$, so $g$ is nonzero and not divisible by $F$.

The endpoint bounds give

$$
\frac{T}{b/a}=i+\ell a+j\frac ab\le(\ell+2)a,
\qquad
T=i\frac ba+j+\ell b\le(\ell+2)b.
$$

Every nonnegative exponent below that weight has total degree at most
$T/\min(b/a,1)$. Consequently

$$
\deg g\le3\ell\max(a,b). \tag{15}
$$

**Bézout charges the wall.** Since $F$ is irreducible and does not divide
$g$, the affine intersection quotient has length at most $D\deg g$.
It surjects onto $K[x,y]/I$, hence

$$
n\le D\deg g\le3D\ell\max(a,b)\le3DK(a+b). \tag{16}
$$

This is the classical Bézout length bound; compare
[Milne, *Elliptic Curves*, I.1.18](https://www.jmilne.org/math/Books/ectext6.pdf).
For completeness, it also follows over $K$ directly: the coprime
homogenizations of $F,g$ form a regular sequence in $K[X,Y,Z]$.
The degree-$h$ quotient has dimension $D\deg g$ for large $h$, by its
two multiplication exact sequences. Evaluation at $Z=1$ surjects that
component onto the affine intersection quotient for large $h$. This
proves the length inequality, including multiplicities and arbitrary
characteristic, without assuming all points are rational.

Equations (14)--(16) allow at most $3D(n-1)$ large walls. Including the
other two types gives

$$
g(I)\le1+D+2(D-1)(n-1)+3D(n-1)\le5Dn,
$$

proving (13). The budget is shared by all walls; independently realizable
pairs from different quotients would not satisfy this argument.

### The fixed-degree order is sharp

Let $B_D(n)$ be the maximum restricted to ideals containing some
irreducible polynomial of degree at most $D$. For every fixed $D\ge2$,

$$
B_D(n)=\Theta_D(n). \tag{17}
$$

The lower family is $I=(xy-1,f(x))$, where $f$ is monic of degree $n$
with nonzero constant coefficient. In $K[x,x^{-1}]/(f)$, any $n$
consecutive Laurent powers form a basis. The cheapest $n$ monomials
outside $(xy)$ always form such an interval. Every allocation
$\{1,x,\ldots,x^a,y,\ldots,y^b\}$, $a+b=n-1$, occurs, for example
at weights $(2b+1,2a+1)$. Thus $g(I)=n$ over every field. Distinct
nonzero roots of $f$ also give reduced point examples over an infinite
field. For $D=1$, the maximum is instead at most two, attained when
$n\ge2$ by $(y-x,x^n)$; a linear equation has only two possible heads.

The argument proves this restricted order. Its novelty and exact priority
remain unestablished. It neither sharpens the unrestricted $G_2$ exponent
nor completes the three-variable goal. It combines established state
geometry and Bézout, and extends the earlier fixed-monomial-curve screening
bound to arbitrary irreducible plane carriers.

### What a finite arithmetic gain does and does not show

The bounded arithmetic test refutes generic maximality already in the
plane. For a primitive eleventh root $\zeta$, take the seven points

$$
(\zeta^j,\zeta^{7j}),\qquad0\le j<7.
$$

Their complete fan has nine initial ideals, versus eight generic planar
corner cuts. The wall ratios $w_x/w_y$ are
$1/6,1/4,1/2,2/3,1,2,4,6$. The rank rule is ordinary Vandermonde:
monomials with the same $i+7j\pmod {11}$ have the same column, and
any seven distinct classes are independent. Primality is not needed
for that consecutive-row argument. Independent elimination on the
seven points $(2^j,2^{7j})$ in $\mathbb F_{23}^2$ verifies the count.

These points satisfy both $y^2=x^3$ and $xy^3=1$. A fixed-curve
construction therefore cannot turn this finite gain into a power gain.
For example, lift the orbit to size $11^k$ with an exponent $a_k$
satisfying $3a_k+1\equiv0\pmod {11^k}$, and select $7^k$ points.
The equation $xy^3=1$ remains, so (13) gives at most $20\,7^k$
initial ideals. A proposed $9^k$ multiplication law cannot persist.
This is a limitation of that particular fixed-curve lift, not a bound on
every arithmetic family whose relations may grow with its parameters.

The checker saves exact partition witnesses for eighteen fixtures. They
include jets on $x^2y+xy^2+1=0$ in characteristics 5, 7 and 101 and
lengths 7, 12 and 18, exercising the Bézout charging step. Its projective
cubic is nonsingular when the characteristic is not 3: the $Z$ derivative
forces $Z=0$, and at $XY(X+Y)=0$ a remaining derivative is nonzero.
It is therefore geometrically irreducible. Hyperbola and line controls
include smaller characteristics, and the cyclic point example is checked
separately. Every complete cone interval, state jump and large-wall
polynomial lift is replayed exactly. A partition row length $\lambda_j$
encodes all $x^iy^j$, $0\le i<\lambda_j$, without losing basis data.

Before the dyadic theorem, this gave a necessary condition: a planar
family with $g(I)\ge c n^{1+\epsilon}$ would need
every irreducible equation in its ideal to have degree
$\Omega(n^\epsilon)$. Applying the product construction to branches
on a bounded-degree irreducible carrier gives only $O_D(N^2)$ initial
ideals at spatial length $N$. The later dyadic theorem rules out the
planar power gain entirely, including families escaping fixed carriers.

## Growing degree can retain a finite grading

The following historical obstruction is sharper for some growing-degree
families than the general dyadic bound. A finite grading puts all basis
moments in a coarser lattice.
The following is an application of classical state geometry and the lattice
polygon vertex bound, not a claimed new lattice theorem. The underlying
framework and vertex estimate are the ones used by
[Babson--Onn--Thomas, Theorem 2.4 and the proof of Corollary 2.5](https://arxiv.org/pdf/math/0207135).

### The intrinsic lattice and its relative area

Take a reduced Gröbner basis of a finite planar ideal $I$. Let
$\Lambda(I)\subseteq\mathbb Z^2$ be generated by the differences between
the exponents of any two nonzero terms in each of its polynomials. Then
$I$ is homogeneous for the group $\mathbb Z^2/\Lambda(I)$. Conversely,
every group grading induced by the variable degrees for which $I$ is
homogeneous kills these differences. Indeed, a reduced basis polynomial
must be homogeneous: otherwise a nonzero homogeneous component not
containing its head would belong to $I$ while consisting entirely of
standard monomials. This is impossible. These two implications also show
that $\Lambda(I)$ is independent of the chosen reduced basis and term order.
Monomial basis polynomials contribute no differences.

For each standard set $B$, write $\mu(B)=\sum_{u\in B}u$. Any homogeneous
monomial basis has the same number of elements in each group degree,
namely the dimension of that component of $K[x,y]/I$. Consequently all
$\mu(B)$ lie in one affine coset of $\Lambda(I)$. This uses a grading,
not a diagonalization of a group action; it remains true when the
characteristic divides the order of the grading group.

Suppose $\Lambda(I)$ has rank two and index $h$. If $\mathcal B$ is
any collection of initial standard sets and
$P=\operatorname{conv}\{\mu(B):B\in\mathcal B\}$ has dimension two,
then its area in the grading lattice is
$\operatorname{area}(P)/h$. Every displayed moment is a vertex: a
strict positive weight realizing its standard set uniquely minimizes the
basis cost. Apply a lattice basis change and the classical planar vertex
bound to obtain, for an absolute constant $C$,

$$
|\mathcal B|\le 2+C\left(\frac{\operatorname{area}(P)}h\right)^{1/3}.
\tag{17}
$$

If $P$ has dimension at most one, it has at most two vertices, so (17)
still holds. If $\Lambda(I)$ itself has rank at most one, this gives
$g(I)\le2$ directly; a monomial ideal has just one state. In the full-rank
case the universal moment box $[0,\binom n2]^2$ also gives

$$
g(I)\le2+C\left(\frac{\binom n2^2}{h}\right)^{1/3}.
\tag{18}
$$

The low-dimensional case cannot be omitted. An ideal homogeneous for a
free grading can also be graded modulo arbitrarily large integers, while
its moment polygon lies on a line. In dimensions above two, a
lower-dimensional polygon or polytope requires its own intersection
lattice; dividing its relative volume by the full ambient index is not
justified.

### A uniform bound for polynomial graphs

Consider, over any field,

$$
I=(x^n,\ y-G(x)),\qquad
G(x)=\sum_{j\in S}c_jx^j,\quad c_j\ne0,\quad
1\le q=\min S<\max S=D<n.
$$

Put $h=\gcd\{j-q:j\in S\}>0$. The displayed generators are a reduced
basis for lexicographic order with $y>x$. Their exponent differences
give

$$
\Lambda(I)=\langle(-q,1),(h,0)\rangle,
\qquad [\mathbb Z^2:\Lambda(I)]=h.
$$

Thus the exact intrinsic grading is $\deg x=1$, $\deg y=q$ modulo $h$.
The quotient has reference basis $1,x,\ldots,x^{n-1}$, so its dimension
in grade $c$ is $\#\{0\le k<n:k\equiv c\pmod h\}$.

For weights $(t,1)$ with $t<1/D$, the initial ideal is $(x^n,y)$,
giving one state. For $t>1/D$, the head of $y-G(x)$ is $x^D$.
Every other initial standard set therefore lies in the strip $x<D$.
The coordinate fibers of a strip downset have lengths at most $D$,
so all its moments simultaneously satisfy

$$
0\le\mu_x\le\frac{n(D-1)}2,\qquad
0\le\mu_y\le\frac{n(n-1)}2.
$$

Apply (17) to these states, whose polygon has area at most
$n^2(D-1)(n-1)/4$, then add the pure-$x$ state. This proves the
uniform estimate

$$
g((x^n,y-G(x)))\le
3+C n\left(\frac{D-1}{h}\right)^{1/3}.
\tag{19}
$$

All constants are independent of the support, coefficients, characteristic,
$D,h$ and $n$. A monomial $G=cx^q$ has a free grading instead, and
at most two initial ideals.

In particular,

$$
g((x^n,y-x-cx^D))=O(n),\qquad c\ne0,\quad2\le D<n,
\tag{20}
$$

uniformly even when $D$ grows with $n$, because $h=D-1$. This includes
$D=p^s$ for prime $p$; the field need not have characteristic $p$.
More generally, $G=x^qH(x^h)$ with $D/h$ bounded has the same linear
upper bound. Neither prime-power degrees nor repeated spacings by
themselves provide the desired growth.

For a polynomial-graph construction with
$g(I)\ge c n^{1+\epsilon}$ along an infinite family, (19) forces

$$
\frac{D-1}{h}=\Omega(n^{3\epsilon}). \tag{21}
$$

In particular, realizing $\Omega(M^{4/3})$ parametric bases at colength
$n=\Theta(M)$ in this class would require $D=\Theta(n)$ and bounded
$h$. A growing degree with a growing proportional spacing does not
satisfy this condition. This was a necessary condition, not a realization
theorem. The later dyadic argument supplies the unrestricted planar
obstruction independently of this lattice calculation.

### Exact controls and what they decide

The public checker exhausts the positive fans of **24 fixtures** over
$\mathbb F_2,\mathbb F_3,\mathbb F_{101}$: six nonmonomial graph
models, one free-graded monomial graph and one monomial ideal per field.
The nonmonomial supports are $(1,4)$, $(1,5)$, $(1,3,5)$,
$(1,3,4)$, $(1,4,7)$ and $(1,4,5)$, at lengths 7, 12 or 18.
The report supplies exact partitions, strict weights and complete cone
intervals. Independent residue and inverse-matrix checks verify the
grading dimensions, support-difference index, affine moment coset,
relative polygon area and both strip budgets. Cases with the
characteristic dividing $h$ are included.

For example, in characteristic 101 at length 18 the commensurate support
$(1,4,7)$ has index 3 and 12 states; $(1,4,5)$ has index 1 and 16 states.
Other characteristics can collapse additional dependencies. These counts
are diagnostics for the proof and do not establish an asymptotic gain
from symmetry breaking. No novelty certificate or unrestricted bound
follows from the fixtures.

## The additive prime-power class has a sharp small count

There is another way a growing-degree construction can remain simple.
Fix a prime power $q=p^s\ge2$ and a field $K$ of characteristic $p$.
Call a polynomial **affine $q$-linearized** here if it has the form

$$
c+\sum_{i=1}^d\sum_{r\ge0}c_{ir}x_i^{q^r},
$$

with only finitely many nonzero coefficients. The constant is allowed.
This is an algebraic condition on equations in characteristic $p$;
it is not a description of every prime-power sequence or of the
original smooth-number order over real logarithmic weights.

The reduction used below is established linearized-polynomial machinery.
[Kuijper--Trautmann, Proposition 7 and Theorems 11--12](https://arxiv.org/pdf/1406.4600)
treat minimal bases, distinct leading positions and Gröbner bases for
modules over the composition ring of univariate linearized polynomials
over finite fields. Their setting differs from ordinary multivariate
polynomial ideals. The elementary bridge here is that $g^{q^r}$ belongs
to the ordinary ideal $(g)$, and distinct pure-variable heads are
relatively prime. No new Frobenius reduction, module algorithm or
coding-theoretic result is claimed. The exact counting formulation's
priority has not been established.
The later published treatment is
[Horlemann-Trautmann--Kuijper (2018), Definition 2.14, Proposition 2.16
and Theorem 2.17](https://dergipark.org.tr/en/download/article-file/388999).
Its Proposition 3.8 analyzes an interpolation-module algorithm in
$q$-degree, reinforcing the need to credit the specialized algorithmic
literature rather than compare only against general ideal algorithms.

### Box classification for every global order

Let $I\subset K[x_1,\ldots,x_d]$ be a proper finite ideal generated by
affine $q$-linearized polynomials. Then its length is necessarily
$n=q^k$ for an integer $k\ge0$. For every global monomial order,

$$
\operatorname{in}(I)=
(x_1^{q^{a_1}},\ldots,x_d^{q^{a_d}}),\qquad
a_i\ge0,\quad\sum_i a_i=k. \tag{22}
$$

Its reduced basis has exactly $d$ affine $q$-linearized polynomials.
In particular,

$$
g(I)\le\binom{k+d-1}{d-1}. \tag{23}
$$

**Proof.** Start from the given generators. If two nonzero polynomials
$f,g$ have heads $c x_i^{q^u}$ and $b x_i^{q^v}$ with $u\ge v$, replace

$$
f\longleftarrow f-c\,b^{-q^{u-v}}g^{q^{u-v}}.
$$

Retain $g$ and delete a zero result. The ordinary ideal is unchanged:
the subtracted power is a multiple of $g$, and the operation is reversible
with $g$ retained. Frobenius preserves the affine $q$-linearized form.
The head of a power is the power of the head for every global order,
so each replacement strictly decreases the head of $f$. The finite
multiset of heads decreases in its well-founded multiset order; the
process terminates. A nonzero constant cannot occur because $I$ is proper.
No perfectness assumption or extraction of coefficient roots is needed.

At termination the remaining heads involve distinct variables. They
are pairwise relatively prime, so the ordinary Buchberger product
criterion makes the remaining generators a Gröbner basis for the
unchanged ideal. Finite colength forces a head in every variable.
The standard set is consequently the box in (22), with length
$q^{\sum_i a_i}$. Frobenius substitutions also reduce all tails while
preserving their affine form and the heads. They yield the ordinary
reduced basis: every tail monomial is standard. This proves (22)--(23).

### Universal weights and sparse output

Every realized allocation $a$ has a coefficient-independent strict weight
witness

$$
w_i=q^{k-a_i}. \tag{24}
$$

Indeed, its reduced head $x_i^{q^{a_i}}$ has weight $q^k=n$.
A nonconstant tail $x_j^{q^r}$ has $r<a_j$, hence weight at most
$n/q$; a constant has weight zero. Thus every global initial ideal
occurs at one of the $\binom{k+d-1}{d-1}$ weights (24). If a composition
is not realizable for a particular ideal, its weight simply gives
another realized state. No coefficient-dependent fan arrangement
must be enumerated to obtain a complete list of representative weights.

There is also a uniform interior region in each realized cone:

$$
\frac{\max_i(w_iq^{a_i})}{\min_i(w_iq^{a_i})}<q
$$

suffices to keep every head above every tail. A reduced polynomial has
at most $k+2$ terms: one head, at most $\sum_j a_j=k$ possible
nonconstant tail terms, and a constant. The universal reduced basis
therefore has at most $d\binom{k+d-1}{d-1}$ polynomials and
$d(k+2)\binom{k+d-1}{d-1}$ coefficient entries. This is a sparse
output bound in field elements, not a bit-cost or running-time claim;
coefficients may themselves require large representations. The
published module methods remain the relevant algorithmic prior art.

There is also an elementary arithmetic-time consequence, without a claim
of algorithmic novelty. Suppose $m$ input generators have maximum
$q$-degree $M$ (largest $r$ in $x_i^{q^r}$). Reduction at equal weights
keeps all $q$-degrees at most $M$ and computes $k$ from its pure heads;
in particular $k\le dM$. At any canonical weight (24), the largest
input term weight is at most $q^{M+k}$, and the smallest variable
weight is at least one. Head elimination never increases that maximum.
Every intermediate pure term therefore has $q$-degree at most $M+k$.
There are at most $H=d(M+k+1)$ such terms. Each input generator's
head decreases at most $H$ times. Reducing the largest reducible tail
first takes at most $H$ further reductions per final polynomial,
since newly introduced terms are strictly smaller.

Thus there are at most $(m+d)H$ reductions per canonical weight.
Each operates on at most $H+1$ terms; coefficient Frobenius powers
can be evaluated by repeated squaring in polynomially many field
operations in $(M+k)\log q$. Enumerating (24) therefore computes
the ordinary universal reduced basis in a number of field operations
polynomial in $m,M,k,\log q$ for fixed $d$. Its exponents are stored
in binary. This consequence uses the standard reduction rather than
expanding a length-$n$ multiplication table; it does not bound field
element bit complexity or improve the specialized published module
algorithms. The checker uses largest-tail-first reduction and verifies
the step budget at the canonical weights.

### The composition bound is attained

For every $q,k,d$, take algebraically independent elements $Z_{ji}$,
$1\le j\le k$, $1\le i\le d$, over $\mathbb F_q$, and set
$K=\mathbb F_q(Z_{ji})$. Let $V\subset K^d$ be the $\mathbb F_q$-span
of the $k$ row vectors $v_j=(Z_{j1},\ldots,Z_{jd})$. It consists of
$q^k$ distinct points. Its point ideal is generated by $q$-linearized
polynomials, as follows.

The first-coordinate subspace $U$ has dimension $k$. Its subspace
polynomial $A(T)=\prod_{u\in U}(T-u)$ is $q$-linearized. One direct
induction uses

$$
A_{U+\mathbb F_q v}(T)=A_U(T)^q-A_U(v)^{q-1}A_U(T).
$$

Each other coordinate is a $q$-linearized polynomial $L_i$ in the
first coordinate: solve the $k$ equations $L_i(Z_{j1})=Z_{ji}$ using
the powers $T,T^q,\ldots,T^{q^{k-1}}$. The resulting interpolation
matrix is invertible. Linearity over $\mathbb F_q$ extends equality
from the $v_j$ to all of $V$. Thus

$$
I(V)=(A(x_1),\ x_2-L_2(x_1),\ldots,x_d-L_d(x_1)).
$$

Here is both the invertibility justification and the sharpness argument.
For any allocation $a$ with $\sum_i a_i=k$, evaluate the $k$ pure terms

$$
S_a=\{x_i^{q^r}:1\le i\le d,\ 0\le r<a_i\}
$$

on $v_1,\ldots,v_k$. The determinant is nonzero: different determinant
permutations give different monomials in the independent $Z_{ji}$,
each with coefficient $1$ or $-1$, so they cannot cancel in any
characteristic. The interpolation matrix above is one such case.

For each $i$, express the column of $x_i^{q^{a_i}}$ in the columns
$S_a$. The resulting $q$-linearized relation vanishes on the basis
vectors and hence on all of $V$. At (24) its head is strictly the
displayed pure power. Therefore $\operatorname{in}(I(V))$ contains
the box ideal (22). Both quotients have length $q^k$, so they are
equal. Every composition is realized, proving equality in (23).
This also holds on a nonempty open set of specializations over
sufficiently large finite extensions of $\mathbb F_q$.

### A multiplicative relation changes the scale

The additive restriction is essential. In the same characteristic,
take $n=q^k$ and

$$
J=(x^n-x-1,\ xy-1).
$$

The first equation is affine $q$-linearized; the second has a mixed
product. This is the earlier hyperbola family, with exactly $n$
initial ideals. The polynomial $x^n-x-1$ has nonzero constant and
derivative $-1$, so over an algebraic closure the example consists
of $n$ distinct points. Its count grows exponentially in $k$, versus
the additive planar maximum $k+1$. It remains linear in $n$ and
therefore does not supply the sought planar power gain or complete
the unrestricted three-variable goal.

This distinguishes two mechanisms: prime-power degree with additive
relations remains a small-box problem, while multiplicative relations
can already destroy that bound. To exceed the generic scale, a new
construction still needs quantitatively stronger interactions than
this one fixed-degree hyperbola. The result does not bound cyclic
root-of-unity prefixes in characteristic coprime to their orbit size.

### Independent verification

`scripts/run_additive_initial_ideals_checks.py` uses only the standard
library and existing exact geometry helpers. It separately implements
Frobenius reduction, ordinary polynomial division and S-pairs, and
Gaussian selection on actual quotient columns. The controls include
nilpotent additive graphs, a noncyclic truncated tensor algebra,
reduced additive point spaces, an affine translate, and base $q=4$
in characteristic 2. The extension field is represented explicitly
as $\mathbb F_2[T]/(T^8+T^4+T^3+T+1)$; every nonzero element's inverse
is checked rather than assuming that modulus is irreducible.

The report saves **14 additive fixtures, 85 distinct initial bases and
2,349 ordinary S-pair checks**, as well as two independent complete
hyperbola fans of lengths 8 and 9. Two- and three-variable fixtures
cover all relevant pure-$q$ comparison chambers; the four-variable
fixture realizes every composition. Exact specialization matrices,
points, generators, reduced bases, boxes and strict witnesses are
saved. None of these finite checks replaces the uniform proof or
establishes novelty. The classification closes an additive lower-bound
route; it does not improve the unrestricted extremal exponent.

## Next proof task and stopping rule

The next task is now directly three-dimensional. The
[dyadic planar bound](planar_dyadic_bound.md) excludes every planar
$n^{1+\epsilon}$ family, including the proposed low-overhead transversal
realization. The common-coordinate product of a fixed number of planar
branches has at most the generic power and logarithmic scale. This is
a uniform obstruction to the route, not a lack of suitable small examples.

The [rational-shadow proof](rational_shadow_bound.md) now extends that
mechanism: every linear functional on a convex token path has at most
two monotone phases. Projected dyadic covers have only logarithmically
many distinct boundary lines, and their internal lattice budgets sum
to a single logarithm. The resulting height-dependent bound controls
all rational shadows. Short integer relations then make every vertex
with bounded incident facet normals visible in a finite set of these
shadows. This yields the stronger uniform three-variable upper bound
above, without assuming a single line counts all spatial cells.

The next proof must close more of the remaining power gap, not merely
prepare another manuscript or enumerate more examples. One concrete
upper route is to improve the charge for vertices incident to large
facet normals by using the common quotient and its commuting
multiplication; the current $O(n^4/h)$ surface charge uses no such
structure. Another is to prove fewer arithmetic projection directions
suffice. A lower route still requires genuinely three-variable
dependencies. The exact planar and rational-shadow priority comparisons
remain necessary. The affine $q$-linearized classification excludes that
additive class as a power-gain construction in dimension three.

The [five-cone family](initial_cone_height.md) now closes the proposed
uniform $O(n)$ integer-weight route: its exact minimum is
$(n+1)^2/4$ along odd lengths. The corresponding square-root dependence
on primitive cone-ray height is also sharp. These are genuine commuting
quotient algebras, rather than unattached rational cones. Neither a hard
weight representation nor a large necessary projection radius implies a
large initial-ideal count: exactly two chosen planes cover all five cones.
The next projection theorem must therefore bound a selected cover or
its incidence multiplicities globally, instead of improving those
individual-cone exponents. The unrestricted upper and lower targets
are unchanged.

The attempted aggregate estimate is now closed by a realizable
counterexample. Write $T(I)=\sum_w\|w\|_\infty d(w)$ over positive
primitive fan rays. The reduction
$g(I)=O(n^{5/3}\log(n+1)+nh^{5/3}\log(n+1)+T(I)/h)$ is valid,
but its classical $T=O(n^4)$ input cannot have a uniform power saving.
For every monic $F$ of degree $n$ with $F(0)F(1)\ne0$,

$$
I_F=(xy-1,\ (x-1)z-1,\ F(x))
$$

has the same positive fan as choosing the $n-1$ cheapest points on
three coordinate axes. Partial fractions prove that every pure-axis
allocation of size $n-1$, together with $1$, is a basis. Thus
$g(I_F)=\binom{n+1}{2}$ and the fan has $(n-1)^2$ interior rays.
A uniform coprime-triple count proves $T(I_F)=\Theta(n^4)$ along
$n\equiv3\pmod6$. The construction works over every field with
$F=t^n-t+1$, and with radical rational points by using a split $F$
over $\mathbb Q$. This refutes $T=O(n^3)$ and every uniform power
saving over the classical quartic incidence exponent. The quotient residue
matroid differs from a simpler uniform-axis matroid with the same
positive fan; confusing those two objects would miss the realization.

At $h=n^2/32$ the family has $\Theta(n^2)$ wholly interior bad cells,
so the order $n^4/h$ is sharp at that threshold. This does not exclude
a better tail estimate at the smaller threshold used in the existing
upper bound. That threshold dependence, or a smaller selected family
of shadows, is the next upper route. Multiplication compatibility
alone cannot make the proposed full moment cubic. A fixed three-pole
curve already supplies the obstruction; no new larger search framework
or standalone paper is warranted by it.

Every plane meets at most $2n-1$ positive cells of this family, by a
decreasing type-rank potential for the three axis streams. Any selected
plane cover consequently needs $\Omega(n)$ members, even with real
coefficients. A polylogarithmic cover for all quotients is false.
For $k$ distinct primitive integer kernels the elementary height count
gives $\sum\|q\|_1^{1/3}=\Omega(k^{10/9})$. Thus the present rational
shadow envelope, summed separately over every plane in any cover of
this family, has budget at least $\Omega(n^{19/9}\log(n+1))$.
That envelope alone cannot reach a sharp $n^2$ power up to logarithms;
the positive traces here have only linear size. The next route must use
quotient-sensitive trace bounds, overlap accounting, a threshold-specific
tail estimate or a different global argument. This is a limitation of
the stated proof budget, not a lower bound on arbitrary enumeration
algorithms or unranking.

The separate nine-cone all-field quotient has a fully interior cell
whose exact minimum integer height is $(n+1)^2/4-2$ and exact minimum
kernel height $(n-1)/2$. Boundary removal therefore cannot restore a
smaller uniform individual radius. The Sabariego--Santos comparison was
expanded to its local-star Hilbert-basis bound; it concerns the star at
the ordinary-degree ray, rather than all positive rays. The strongest
current unrestricted counting comparison and the main endpoint remain
open.

The nonlinear rectangle route remains available. Its filtration bound
identifies a quantitative window, and the rational relation distinguishes
the fixed surface from free-coefficient embeddings. A proposed count must
handle a defined class uniformly. The sparse minor witness $F_*$ is a rank
control, not an extremal family or a reason for larger fan enumeration.

A useful outcome would be a uniform count or an explicit collection of
distinct cones with a provable growth rate. Only a rate polynomially above the
generic scale, with a further general consequence, meets the lower-family
endpoint of this goal. An upper bound excluding that rate would close the
route and justify a different construction.

The new comparisons isolate the missing mechanism: a quantitative count of
coherent basis changes caused by overlapping algebraic dependencies in
interior weight regions. Generic torus-grid graphs lead back to ordinary
cuts. Nonzero-prefix rank defects are real in local algebras, but their
existence or size does not imply a large fan. The next proof must count their
interactions, or bound them strongly enough to close the local route.

Increasing the finite cases, fitting a curve to four counts, or replacing
coefficient rank by support matching cannot serve as that proof. The source
audit must also identify any existing theorem that already covers a proposed
advance before it is presented as new. The existing lattice-shadow work remains
available as a comparison tool, with its priority issue explicit. The
T-graph audit adds a second limitation: individually realizable neighboring
ideals must be controlled simultaneously in one fan. Pairwise arrow-map
conditions do not supply that global count. Any next proposed inequality
must use the common quotient or additional global state geometry, and must
give a quantitative estimate rather than another necessary condition.

## Rectangle ranks and high-order perturbations

The [rectangle rank comparison](rectangle_rank_structure.md) makes two
further distinctions precise. For every $F=\alpha x+\beta y+H$ in a
rectangular algebra, with $\alpha\beta\ne0$, a coordinate-ideal-preserving
automorphism sends $x+y$ to $F$. Hence the image dimensions of all lists
of mixed monomial multiplication operators are fixed, over every field.
They do not determine the selected monomial ranks or the fan: already
$(x^2,y^2,z-x-y)$ and $(x^2,y^2,z-x-y-xy)$ have respectively five and
six initial ideals over $\mathbb Q$.

For any homogeneous plane algebra with top degree $S$, a perturbation
starting at degree $D$ with $2D>S+1$ leaves every component strictly
between degrees $S-D+1$ and $D$ unchanged. The low and high components
together form an exact direct summand of rank $R$. The full monomial
matroid is the direct sum of this part and the unchanged degree pieces.
The homogeneous graph has at most six initial ideals, so the classical
three-dimensional lattice-vertex bound gives
$g(I_F)=O((1+RS)^{3/2})$. In a square rectangle, $S=2q-2$ and
$R=(2q-D)^2$ for $D\ge q$. Any construction with
$g(I_F)\ge c n^{2+\epsilon}$ therefore needs
$2q-D=\Omega(n^{5/12+\epsilon/3})$, where $n=q^2$.

This is a uniform restriction on terminal perturbations, not a new
unrestricted exponent or a completed novelty assessment. General
low-order nonlinear terms remain outside the compression. The next
lower theorem needs an exact selected-column rank or coherent-exchange
description for such terms. Strong Lefschetz and Jordan classifications
are established inputs, and the final lattice estimate is classical;
neither supplies that missing description.

## The exact profile lattice separates free and finite gradings

The [profile-lattice continuation](profile_lattice_bound.md) identifies
the lattice generated by all monomial basis-profile differences with
the support-difference lattice of any reduced Gröbner basis. The ground
must include the one-step border of the universal staircase ground;
otherwise some defining differences disappear. This is a grading
statement over the integers, including torsion in small characteristic.

The rational-shadow charging proof can use the lattice of exchange
directions, rather than the coarser lattice of all projected ground
weights. Its additive constant per box must then be retained. When the
profile lattice has rank two, Cauchy--Binet makes the normalized
projection factor at most one, independently of the plane's height.
This proves $g(I)=O_d(n\log(n+1))$ for the whole rank-at-most-two class,
and an $O_d(n\log(n+1))$ support-oracle call bound for its full profile
polygon. In three variables a nonzero integer grading suffices. The
known planar lower examples show sharp order for this broader class,
without proving a matching lower bound restricted to positive gradings.
The universal reduced-basis cardinality is $O(n^{5/3}\log(n+1))$ in
three variables, by the established generator-packing argument.
Every two-dimensional face of an unrestricted quotient's profile
polytope has the same $O_d(n\log(n+1))$ vertex bound, because a
matroid face is a direct sum of minors. This gives a local bound on
the number of cones incident to a positive ray, while leaving the
total number of rays uncontrolled. The facet bound is sharp over all
such profile polytopes by a planar example lifted with $z-1$; this
does not prove sharpness at wholly interior rays.

For a square nonlinear rectangle, the exact full-rank lattice index is
$h=\gcd\{a+b-1:(a,b)\in\operatorname{supp}H\}$. Classical normalized
volume gives $g(I_F)=O(n^{9/4}/\sqrt h)$. In particular, one homogeneous
nonlinear correction of degree $D$ with $D-1=\Omega(\sqrt n)$ cannot
give a polynomial separation. The terminal-band comparison also gains
this factor, with its full dimension verified before dividing by the
index. Low-order nonlinear terms with mixed degree spacings can have
$h=1$ and remain open.

The available Katoh--Tokuyama concave-surface preprint explicitly
excludes general two-parameter matroids from its proof. The smaller
level exponent cannot be imported without that missing reduction.
The current comparisons do not certify novelty of the exchange-lattice
theorem or complete the strongest-current-bound audit. The unrestricted
endpoint, original success criteria and earlier proofs are preserved.

## Generic multiaffine carriers and the nonuniform-rank boundary

The [generic multiaffine screening proof](generic_multiaffine_bound.md)
gives $g(I)=O(n^2\log(n+1))$ for sufficiently general point configurations
on every geometrically irreducible multiaffine surface in three variables.
The constant is independent of coefficients and characteristic, with the
usual extension-field qualification when a finite field has no rational
points in the required open set. A common list of the same order meets
every initial cone across this entire class.

The mechanism is exact. Generic evaluations truncate the carrier's
monomial rank function. Its possible heads are squarefree, so their
standard grounds lie on the orthant boundary and have only
$O(n\log(n+1))$ relevant monomials. A small perturbation puts that ground
in convex position while preserving its strict prefix sets. The classical
three-dimensional convex-position $k$-set formula then supplies the bound.
The exact primary formula was inspected in Leroux--Rademacher, Section 2.5
and Proposition 2.13. This is a geometric corollary, not a new general
$k$-set theorem or a completed priority assessment.

For $xyz=1$, an explicit geometric-progression point family gives ordinary
Vandermonde columns and realizes the full canonical-prefix fan. Under
Laurent coordinates these are lattice points in triangles with vertices
on three fixed rays. The three-way pattern therefore has a precise
realization, but cannot provide the requested polynomial separation.
The checker covers 20 full fans, 318 independent comparison chambers,
756 symbolic minors (including characteristics two and three), eight
convex-position configurations, eight explicit toric quotients and four
additional carrier point quotients. Two complete common lists are checked
against both fixed-head and changing-head carriers.

The generic condition cannot be omitted: reduced special points on
$xyz=1$ already give a dependent first canonical prefix although no
monomial evaluation is zero. This leaves the preserved local rectangle
family, and special configurations with a prescribed rank rule, logically
open. Further generic interpolation on multiaffine carriers is excluded
as a lower-bound route. The unrestricted $G_3$ target, its success
criteria, the quotient-sensitive upper route and every earlier result
remain in force; this restricted corollary does not complete the goal.

## Bilinear rectangles separate Frobenius torsion from characteristic zero

The [bilinear rectangle proof](bilinear_rectangle_fans.md) supplies an exact
nonuniform rank model and classifies its matching-characteristic prime-power
subclass. Translation takes $(x^r,y^s,z-x-y-xy)$ to
$((X-1)^r,(Y-1)^s,Z-XY)$ and preserves every global initial ideal.
It does not preserve arbitrary selected-column ranks; it preserves ranks
on monomial downsets by a unit triangular column change.

The translated columns are
$\binom{a+c}{i}\binom{b+c}{j}$ for $0\le i<r,0\le j<s$.
In characteristic zero this is rectangular polynomial interpolation at
the integer nodes $(a+c,b+c)$. In characteristic $p$, full prime-power
side lengths collapse these columns to independent periodic residue
classes. This makes the field dependence explicit through established
Frobenius/Lucas arithmetic. The earlier invariance of operator-image
ranks as $F$ varies within one field still does not determine these
selected-column ranks.

For the torsion ideal $(X^r-1,Y^s-1,Z-XY)$ with $2\le r\mid s$,
the exact all-field count is $\binom{r+1}{2}+1$ when $r=s$, and
$\binom{r+2}{2}$ when $r<s$. The proof reduces each shortest
representative to at most two candidates and gives a complete list of
positive weights of height $O(rs)$, improving to $O(r)$ for equal sides.
Thus local bilinear rectangles with prime-power sides in the matching
characteristic have only $O(n)$ initial ideals. Scaling extends this
conclusion to every nonzero choice of the three coefficients of
$\alpha x+\beta y+\gamma xy$.

The checker independently verifies 12 complete torsion fans, 121 initial
ideals, 2,123 ordinary binomial S-pairs and 31,259 normal forms. Direct
local products, translation expansions and Pascal columns agree on
12,054 downset-prefix ranks. Rational border bases and exact full-cone
coverage certify 6, 16 and 32 characteristic-zero states for square
sides 2, 3 and 4, compared with 4, 7 and 11 in their prime-power
characteristics. These finite controls give no asymptotic law in
characteristic zero.

The lattice/shortest-path dictionary and interpolation are classical;
the exact restricted classification has no completed priority audit.
This closes the prime-power bilinear route, while leaving the structured
characteristic-zero interpolation count and general mixed nonlinear
rectangles open. Generic tensor-code rank results do not automatically
apply to these fixed integer nodes. The unrestricted objective and all
earlier criteria remain unchanged; this is not goal completion.

## Exact outer cones and a complete graphic-dual rank model

The [outer-cone theorem](rectangle_boundary_fan.md) gives an exact
prefix rule for $F=y+x(1+y)^d$ in characteristic zero. In the region
$\min(w_x,w_z)>(s-1)w_y$, one iteration of greedy selection separates
into $x$-degree blocks. Within each block the relevant Taylor columns
become disjoint intervals of integer powers of $1+y$; an ordinary
Vandermonde determinant gives their exact prefix ranks. A common list
of at most $2rs=2n$ weights, all of height at most $2n$, serves every
integer $d\ge1$ and nonzero coefficient scaling on this region.
For the bilinear case $d=1$, the exact regional count is $r+1$ when
$r\le s$, and $s+2$ when $r>s$.

This also proves a useful obstruction. At square side $q=2a-1$, the
primitive direction $(a,-(a-1),-a)$ is a coherent global wall and has
product $a(a+1)^2=\Theta(n^{3/2})$. Hence the rank-sized direction
cross that closed the additive nested construction cannot hold for
this bilinear family. Its large walls nevertheless coexist with only
$q+1$ states in the outer region; they are not a fan-count lower bound.

For the whole bilinear fan, the fixed tensor interpolation nodes
$(U,V)=(a+c,b+c)$ have cheapest representative cost
$\alpha U+\beta V-(\alpha+\beta-\gamma)\min(U,V)$ when
$\gamma<\alpha+\beta$. Prefixes are convex integer-node polygons.
This closes the cone $\max(\alpha,\beta)<\gamma<\alpha+\beta$
by classical tensor interpolation, but does not determine ranks in
all the remaining polygons. The outer-cone formula actually fails
outside its hypothesis: at side three and $w=(13,10,5)$, the true
standard set has $z^3$ instead of its predicted $y^2$.

There is a different exact rule for **all subsets**, not just prefixes,
in matching characteristic. Let $q_1=p^a,q_2=p^b$ and take side lengths
$q_1-1,q_2-1$. The [cographic theorem](bilinear_cographic_family.md)
identifies the translated monomial residue matroid with the cographic
matroid of $K_{q_1,q_2}$. Frobenius makes its complete relation space
exactly the row and column sums of the periodic rectangle. For a set
$S$ of residue labels,

$$
\operatorname{rank}(S)=|S|+1-c(K_{q_1,q_2}\setminus S),
$$

where components include isolated vertices. Thus a basis is the
complement of a spanning tree. Using least-cost representatives,
the initial standard set is the complement of a **maximum** spanning
tree. The earlier divisible-period representative formula and a
classical tree algorithm compute it in $O(n\log n)$ comparisons
and $O(n)$ storage. This concerns one weight and arithmetic comparisons,
not all fan cells, bit complexity, or dense reduced-basis output.

The finite controls certify twelve characteristic-zero regional fans,
252 direct selections over two fields, 5,364 block-prefix ranks,
three rational coherent facets, and thirty whole-fan polygon comparisons.
The graphic-dual controls separately check 1,192 residue-subset ranks,
36 tree-selected standard sets in six equal/unequal models, every one
of 4,177 small spanning trees, and 690 independent cycle inequalities.
Their proofs use classical interpolation, graph incidence and matroid
duality. Exact restricted-result priority and the strongest-current
unrestricted comparison remain unresolved.

The next counting task is now concrete: classify the maximizing trees
for the **specified periodic costs**, or exclude this family by a
uniform upper bound. The enormous number of abstract spanning trees
does not show that the three weights realize them. Likewise existing
one-parameter results with freely chosen affine edge weights require
a realization argument before they give any lower bound here.
The separate characteristic-zero polygon-prefix route and the all-field
upper route remain available. The original sharp $G_3$ or polynomial
separation endpoint, its consequence and paper requirements are unchanged.

## Ordered trees, sparse costs and a corank screen

The [ordered cographic continuation](cographic_order_ideals.md) closes
the exact prior-art comparison for the underlying two-circuit duality:
Hidaka--Itoh Section 5.2 and Fu--Ma--Wang Proposition 3.5 explicitly
contain that identity. Martin--Reiner give the cyclotomic-prime version.
The earlier rank model is an application of established mathematics,
rather than a new graph-duality result. These primary passages were
inspected; the full papers and the priority of the additional ordering
statements have not been audited.

The polynomial realization extends to
$J_{r,s}=(S_r(X),S_s(Y),Z-XY)$ over every field and arbitrary periods,
where $S_q(t)=1+\cdots+t^{q-1}$. For equal periods $q$, every basis
that is a monomial downset in the $XY$-head region omits a spanning
tree whose only possible branching vertices are row $q-1$ and column
$q-1$. Its $Z$-axis has length $q-1$ or $q$. The proof compares
neighbor sets of consecutive rows: at most one neighbor disappears,
whereas two shared neighbors would give a four-cycle. It applies to
noncoherent downset bases as well; it does not count them as initial sets.

There is also a weighted algorithmic consequence. The periodic edge
cost is affine on the pieces determined by three diagonal lines. An
edge away from their width-four strips and the square boundary is the
unique minimum of an explicit four-cycle, and hence cannot belong to
a maximum tree. At most $19q$ candidate edges remain. For generic
weights, the omitted tree can be computed and stored in
$O(q\log q)$ arithmetic time and $O(q)$ words. Together with the
on-demand least-representative rule, it represents the standard set
using $O(\sqrt n)$ words, with deterministic logarithmic membership
queries. Listing all $n=(q-1)^2$ monomials still requires linear output.
These bounds do not count bits or dense reduced-basis coefficients.
The candidate graph changes with the prices, so its sparsity alone
does not bound the number of fan cells.

Consecutive binomial nodes strengthen a different supporting result:
the bilinear outer-cone theorem and its exact regional count now hold
in every characteristic. The relevant integer determinant is one,
even when the prefix length exceeds the characteristic. Thus coherent
directions with product $\Theta(n^{3/2})$ occur inside the regular
cographic class too, including radical characteristic-zero geometric-sum
quotients. Compressed basis size and primitive direction product remain
different quantities.

A classical general screen also limits the next construction search.
If every standard set fits in one box of $M$ monomials and the quotient
length is $n$, matroid duality and the established three-dimensional
lattice-vertex bound give $g(I)=O(1+(M-n)^{3/2}\sqrt M)$.
A box with $M=O(n)$ cannot yield a polynomial separation above $n^2$.
In particular, three equal odd-period circuit factors, mixed through
$X=UV,Y=UW,Z=VW$, give a realizable nonregular quotient of length
$(q-1)^3$ with only $O(n^{3/2})$ possible initials. This is a classical
corollary, not a new unrestricted exponent or lattice theorem.

The distinct-prime-period pair-product construction is not closed by
that small-box argument. It realizes
$\mathbb Q(\zeta_{p_1p_2p_3})$ in three generators, with the known
tensor/simplicial rank model and a price-dependent representative
problem. Its elementary period box has size $(p_1p_2p_3)^2$.
Its full coherence count remains to prove; the known number of abstract
simplicial bases is not a lower bound for that count.

The independent controls enumerate all 84 small downset bases, compare
42 sparse/dense tree problems through period 127, check 114,598 explicit
excluded-edge cycles, and compare 54 ordinary polynomial standard sets
over three characteristics. Actual coherent facets and 660 unit
binomial determinants test the small-characteristic argument. The
triple-factor check verifies actual polynomial relations and dual
profile budgets. Exact bounded fans remain finite evidence, not an
asymptotic law: $6,14,26,86$ states at periods $3,4,5,8$ suggest
$2(q-1)(q-2)+2$, which would exclude the equal-period cographic family
too if proved. No larger fan sweep is justified before a classification.

The next task is a uniform coherence count using these constraints,
or a stronger defined prime-period family with a provable count. The
unrestricted sharp $G_3$ or polynomial separation endpoint, additional
consequence, strongest-primary comparison and honest-paper requirements
remain unchanged. None is replaced by the compressed algorithm or the
finite count conjecture.

## Mixed periods: exact box obstruction and a general tower region

The [price-region continuation](cyclotomic_price_regions.md) gives a
structural test of the remaining three-prime construction before full-fan
experiments. It works more generally with pairwise coprime periods
$r,s,t$ over any field, including characteristics dividing the periods.
For $A=K[U,V,W]/(S_r(U),S_s(V),S_t(W))$ with generators
$X=UV,Y=UW,Z=VW$, the quotient length is
$n=(r-1)(s-1)(t-1)$. The three pure-axis minimal-polynomial degrees
are $(r-1)(s-1),(r-1)(t-1),(s-1)(t-1)$. Each pure extent occurs
in a positive initial standard set, so their smallest common coordinate
box has **exactly $n^2$ monomials**. This rules out a common-$O(n)$-box
exclusion, without asserting a lower bound on the number of initials.

There is a coefficient-independent regional theorem beyond this example.
If $F=K[x]/f$ has dimension $D$, $u,v\in F$ are units, and
$A=F[T]/g$ is a monic degree-$R$ extension, use generators
$X=x,Y=uT,Z=vT$. In the strict region
$\min(w_Y,w_Z)>(D-1)w_X$, each $T$-degree block finishes before
the next one starts. Its selected columns are the greedy basis of
$x^a(v/u)^c$ for $a<D,c\le k$, with cost $a+\rho c$ and
$\rho=(w_Z-w_Y)/w_X$. The only possible breakpoints are
$j/h$, $|j|<D,1\le h<R$. Consequently fewer than $2n=2DR$
integer weights of height at most $2n$ cover the whole region when
$D,R\ge2$. The edge cases are stated separately. The result needs
neither a field extension nor separability; a coefficient-dependent
monic $g$ and nilpotents in $F$ do not change its proof.
The linear order is sharp along an explicit family of length $n=2D$:
take $f=\prod_{i=0}^{D-1}(x-2^i)$, $g=T^2$, $u=1$, $v=x^D$.
Ordinary Vandermonde interpolation makes every $D$ of its $2D$
degree-one columns independent. The fan has exactly $D+1$ states,
with breakpoints $\rho=D-2j-1$. This proves sharpness of the regional
order, not a lower bound above the generic scale.

This covers all three sufficiently cheap-coordinate regions of the
mixed circuit product with only $O(n)$ states in total. Also every
periodic label has a representative in $[0,r+s+t]^3$, by rounding
in the explicit index-four kernel sublattice. For comparable periods
and bounded price ratios this bounds all relevant standard sets by a
common $O(n)$-element box. The classical dual-profile vertex estimate
then gives $O(n^2)$ states in that region. The global $n^2$-element
box and the regional $O(n)$-element box have different quantifiers.

A polynomial separation in the comparable-period family would therefore
have to arise with unbounded price ratios **outside** the cheap-coordinate
regions. The range is nonempty: balanced periods of size $q$ and prices
proportional to $(1,q,q)$ lie there for large $q$. The next question is
how to control tower blocks before their previous layers have finished,
or prove many coherent exchanges in this specific intermediate range.
Unbalanced periods remain a separate case. Abstract simplicial-tree
counts still provide no coherent lower bound.

The new controls compare 1,033 algebraic predictions with a separate frontier
greedy algorithm using ordinary monic polynomial multiplication, and
check 2,067 strict minimal-border relations. They include eight explicit
linear-sharpness models with additional controls outside the tower region,
nonreduced coefficient-dependent towers and complete rational-interval coverage
for specified bounded regions, not a full-fan search. They also check
615 constructive label representatives and 36 pure-axis ranks, including
period two and a rank-one tower extension.
The circuit/simplicial dictionary, greedy dictionary and final volume
bound remain established inputs. Priority of the regional common-list
theorem and the strongest-current unrestricted comparison are unverified.
The full objective, further-consequence requirement and honest-paper
criterion are unchanged; none is completed by these exclusions.

## Nonnegative gradings: a global linear class and local fan links

The [Hilbert-budget theorem](nonnegative_graded_fans.md) proves a
global restriction beyond the cheap-coordinate regions. If a finite
three-variable quotient is homogeneous for a nonnegative integer grading
with at least two positive variable degrees, its whole positive fan has
$O(n)$ initial ideals over every field, uniformly in the grading heights.
For a strictly positive grading an explicit bound is $12n+6$; for
$(0,p,q)$, normalized with $p,q$ coprime, it is at most
$2+2(n-h_0)$. The constants are not claimed optimal. Explicit even-
and odd-length ribbon quotients have $\lfloor n/2\rfloor+1$ states,
so the combined class has sharp order $\Theta(n)$. This does not prove
the same lower order in the strictly positive subclass alone.

The key argument is a weighted divisor budget, not a new lattice-volume
bound. Among pairs $(a,b)$ whose full divisor rectangles fit a fixed
graded Hilbert function of total size $n$, there are at most $2n$.
For coprime positive pair degrees $p,q$, split by $pa\ge qb$. Write
$b=r+pk$; the degree $pa+qr$ contains $k+1$ divisors, and the map
$(a,b)\mapsto(pa+qr,k)$ is injective. The other half is symmetric.
For a zero-positive pair the corresponding rectangle budget is at most
$n$. A symmetric optimal matroid exchange between actual adjacent
standard sets supplies an endpoint rectangle. In a graded profile
polygon, its positive coordinate pair determines the whole primitive
edge direction. There are therefore only linearly many possible edge
directions. Intermediate exchanged vector bases need not be downsets;
the argument uses the actual endpoint standard sets.

For grading $(0,p,q)$ the proof also gives a Hilbert-function-dependent
common list of $O(n)$ positive integer weights of height at most $4n$.
Its linear worst-case certificate height is sharp: the same ribbons
have an extreme cone requiring height at least $\lfloor n/2\rfloor+1$.
This list is independent of the coefficients within a fixed Hilbert
function; it is not a single linear list for every Hilbert function
of size $n$. Two zero variable degrees are excluded: $(J(X,Y),Z)$
retains the known planar $\Theta(n\log n)$ maximum. Mixed-sign
integer gradings retain the earlier general rank-two bound.

There is also a consequence for arbitrary finite quotients: at any
rational nonnegative weight with at least two positive coordinates,
the associated graded quotient has one of these gradings. Classical
initial-of-initial refinement therefore gives $O(n)$ distinct local
positive refinements. All non-axial nonnegative rays of the positive
Gröbner fan have linear link complexity. The three coordinate-axis
rays retain the broader planar bound. This does not count how many
rays occur, and hence does not settle or improve the unrestricted
$G_3$ exponent by itself.

The tower mechanism can now be stated precisely. For
$F[T]/T^R$, with $F=K[x]/f$ and generators $x,uT,vT$ for units
$u,v$, the old slot formula is valid for all generic positive weights:
interleaving the independent lower blocks changes no dependence.
For a monic deformation $F[T]/g$, every initial set using only
$T$ degrees below $R$ must still be one of those same slot sets.
The only extra border comparisons for such a set come from at most
$R+1$ pure $Y,Z$ heads of degree $R$. Consequently every additional
initial set contains a monomial reaching the defining relation. The previous focus
on incomplete lower-block overlap was too broad; reaching the monic
relation is necessary. In the mixed-period control $(3,5,7)$, the
head $Y^2Z^4$ has a noncheaper normal-form term $X^7$, and an actual
additional state results. This verifies the mechanism, not a lower
count or an unbounded-price family.

The exact controls test 27,440 arbitrary positive Hilbert budgets,
34,776 rectangle injections and 6,860 zero-positive budgets. They walk
51 complete bounded comparison fans, including nonfree quotients,
even/odd sharpness families and local initial forms; check 174 standard
sets, 954 strict border relations, 136 symmetric exchanges and 1,098
multiplication-closure tests; and replay 24 tower deformations.
The enumeration is independent of the proposed Hilbert candidate list.
These are finite supporting checks for the written uniform arguments.

The matroid/state dictionary, edge-ideal framework, symmetric exchange,
initial refinement and border criterion are established inputs.
Altmann--Sturmfels Proposition 3/Theorem 8, Kreuzer--Robbiano
Proposition 2.2/Theorem 2.4 and Haiman--Sturmfels Corollary 3.4/
Theorem 3.6 were inspected in their primary texts within the stated
scope. A full current-priority comparison of the numerical Hilbert
budget and grading bounds is still necessary. No new paper or
publication claim is made. The original sharp $G_3$ or polynomial-
separation criterion remains active and unachieved.

## Signed gradings and the mixed-period prefix mechanism

The [grading study](nonnegative_graded_fans.md) now proves a linear
bound for every integer grading with all three variable degrees
nonzero, including mixed signs. For opposite pair degrees $p,-q$,
the map

$$
(a,b)\longmapsto
 \left(pa-qb,\min(\lfloor a/q\rfloor,\lfloor b/p\rfloor)\right)
$$

injects each admissible rectangle corner into its Hilbert budget.
There are at most $n$ such corners. Actual positive-weight endpoint
exchanges, with the third degree determining the remaining coordinate,
give $g(I)\le12n+7$. Mixed signs require counting the positive-visible
chain of the state polygon, rather than assuming every polygon vertex
can be realized by a positive weight.

The grading $(0,p,-q)$ remains outside this proof. The explicit
socle-identification family
$(XY,XZ,X^{a+1},Y^{c+1},Z^{c+1},X^a-Y^cZ^c)$ has exactly two
states and shows that a fixed positive pair $(c,c)$ does not determine
the wall direction. It disproves that encoding shortcut, not a linear
bound for the class. The earlier planar $O(n\log n)$ estimate remains.

For the original full-rank candidate, put
$F=K[U,V]/(S_r(U),S_s(V))$, $D=(r-1)(s-1)$ and
$A=F[W]/S_t(W)$, of length $n=D(t-1)$. Group a finite monomial
prefix by $j=(b+c)\bmod t$, and let $V_j$ be its coefficient span
in $F$. The diagonal kernel of $F^t\to A$ gives the exact all-field
identity

$$
\operatorname{rank}_A(E)=\sum_j\dim V_j-\dim\bigcap_jV_j.
$$

In the established cographic two-factor model, the intersection is
controlled by the **sum of the complement graph cycle spaces**.
It is not controlled by the cycle space of their union: even two
spanning trees can have a union with cycles. This is a concrete
obstruction to extending the old connectivity algorithm directly.
The tensor/circuit rank identity and graph duality are classical
inputs, not new general matroid claims.

When $t\min(w_Y,w_Z)>(D-1)w_X$, every standard monomial has
$a<D$ and $b+c<t$. Each standard set is contained in a $tD$-element
auxiliary slot set and omits exactly $D$ elements. There are only
$O(n)$ possible auxiliary sets. Their fixed rank deficit gives a
stronger count: if $m_j$ basis columns occur at level $j$, then
$\delta_j=D-m_j\ge0$ and $\sum_j\delta_j=D$. Thus the total
degree profile varies over an interval of length only $D(t-1)=n$.
In the unimodular coordinates $(a,c,b+c)$, the basis-profile polytope
lies in a box of volume at most $n^4$. The established lattice-vertex
bound gives $O(n^2)$ states, superseding the earlier
$O(n^2\sqrt t)$ corollary and excluding this whole region from a
polynomial separation.

The profile proof only needs containment in $a<D,b+c<t$, so it
counts all initial sets with that property at any price. Permuting
the three coordinate presentations shows that the union of all bases
failing to reach at least one of the pair-degree thresholds $r,s,t$
has only $O(n^2)$ members, even for unbalanced periods. A polynomially
separating family must therefore get its asymptotic growth from bases
reaching all three thresholds; the three witnessing monomials need
not coincide. A small exact basis verifies that simultaneous crossing
can occur, but lies in an already excluded bounded-price region and
supplies no infinite lower family.

The same profile-deficit calculation applies to
$F[W]/g$, for a monic divisor $g$ of $W^t-1$ of degree
$R=t-s$, including coefficients in nonreduced $F$. In a one-period
ground its bound is
$O(1+n(t-1)\sqrt{D(D-1)s})$; fixed $s$ gives $O_s(n^2)$.
This is a rank budget followed by a classical lattice-volume estimate,
not a new matroid-duality or lattice-vertex theorem. The main remaining
candidate mechanism is cost-prefix variation across multiple periods
in every coordinate presentation, outside the existing price-region
exclusions. No larger full-fan sweep is justified by the new result.

Independent checks cover signed Hilbert budgets, actual complete
small signed fans, the zero/opposite-sign boundary, full cost-prefix
ranks by ordinary polynomial arithmetic and graph cycle sums, and
one-period slot omissions with strict border certificates. They include
characteristics dividing periods. These controls support the uniform
arguments and do not prove a fan growth law.

The primary comparison now distinguishes fixed coordinates, generic
coordinates and generating degree from finite length. In particular,
Frels--Schmitz's explicit Notation 2.10 example has infinite colength;
truncating at degree $d+1$ preserves its degree-$d$ mechanism but
gives length $d(d+1)(d+5)/6$. Their generic-family theorem is not
claimed to exclude Artinian ideals. This parameter normalization does
not establish priority of the present bound. Exact inspected theorem
scopes and the now-identified Silversmith primary-PDF comparison are
recorded in the grading note.
The strongest current unrestricted comparison, novelty and the main
sharp-count or separation endpoint remain unresolved.

## Periodic covers and coherent arithmetic

The [periodic quotient study](periodic_quotient_arithmetic.md) separates
two different ordering problems. The full periodic algebra has length
$M=rst$, pure-axis orders $rs,rt,st$ and a profile lattice of index
$M$. Relative volume gives the classical corollary $g(B)=O(M^2)$,
uniformly in periods and characteristic. The circuit quotient has
intrinsic profile lattice exactly $\mathbb Z^3$ in every field.
The monic pair minimal polynomial has nonzero constant and linear
terms $1,-1$, forcing each generator's group degree to be zero.
Thus the cover's profile-index saving is unavailable. An ideal
quotient does not give a general monotonicity
theorem for numbers of initial ideals.

The complementary algebra
$C=B/((U-1)(V-1)(W-1))$, of length $M-n$, gives the known dual
matroid through inverse monomials. It has the fixed degree-six
carrier $(Z-XY)(Y-XZ)(X-YZ)$. However, the actual omitted columns
solve a maximum-cost problem on inverse labels in the finite cover
basis. After reflection the ground generally lacks the constant
monomial, so it cannot contain a global standard set of $C$.
Neither its low-degree carrier nor a count of its ordinary initial
ideals resolves those costs without another argument.

There is an actual coherent characteristic effect, not just possible
abstract torsion. For periods $(5,7,11)$ and the stated positive
near-total-degree weight, a rational standard set has integral index
two and loses exactly one rank modulo two. The characteristic-two
standard set replaces $XY^{11}$ by $Y^{12}$ and has index one.
The full exact bases, cleared integral relation and strict border
checks are saved. This rules out treating all coherent ranks in the
candidate as characteristic independent; it gives no growth law or
comparison of total fan counts.

For every ideal, substitution $X\mapsto X^d$ multiplies length by
$d$ and repeats standard monomials by exponent residues. It preserves
the full fan under the corresponding diagonal change of weights,
so the number of initial ideals is unchanged. An index-$\Delta$
coherent integral lattice basis acquires index $\Delta^d$. Applying
this to the witness yields arbitrarily large coherent torsion with
no increase in fan count. This flat-extension consequence eliminates
a tempting repetition mechanism for a lower-bound construction;
these extensions are not asserted to be the original three-prime
family. More generally every finite monomial self-map of affine
three-space is a permutation of pure-power maps: finiteness of the
origin fiber forces a pure power for each axis. All such pullbacks
therefore preserve fan counts. Maps finite only on a unit-valued
quotient and nonmonomial finite maps are outside this exclusion.
The next mechanism must create genuinely more positive cost
regions across multiple periods, not just larger determinants.

The known cyclotomic/simplicial duality and good-prime framework are
explicitly cited. Arnold's inspected homogeneous lucky-prime statements
are distinguished from the elementary finite-free determinant test
used here. No new duality, flatness theorem, unrestricted exponent,
priority certificate or paper is claimed. The original goal remains
active and unachieved.

## Fixed-cover traces and the remaining dimension gap

The [periodic quotient continuation](periodic_quotient_arithmetic.md#what-can-be-proved-inside-a-fixed-cover-cone)
gives two uniform controls for a cover standard set $E$ of size $M$
and quotient corank $k=M-n$. Every omitted monomial lies in the
coefficient-independent upper shell
$U_k(E)=\{a\in E:|\{b\in E:b\ge a\}|\le k\}$.
If $E$ has $\sigma$ maximal monomials, shifted divisor boxes give
$|U_k(E)|=O(\sigma k\log^2(k+1))$. The factor $\sigma$ is retained;
no assertion that it is bounded or that shell subsets are coherent
is used. Restricting the actual dual optimum to this shell is exact.
For periods greater than two, eight appropriately shifted cover
monomials also form a circuit of the inverse-label dual matroid over
every field. Its uniquely cheapest corner is excluded from a maximum
basis. This gives a linear number of membership tests once the cover is
explicitly given. No uniform $O(k)$ size bound, full-cone count or
sublinear cover construction is claimed; a nonregular matroid need
not have only these cube circuits.

Along a coordinate-price line inside a fixed cover cone, the optimal
omission's integral coordinate sum is strictly monotone at each
generic state change. There are at most $1+k(\max_E a_X-\min_E a_X)$
such states. The earlier rational-shadow theorem can also be applied
with rank $k$ and hyperbolic height $M$, retaining the projection's
arithmetic height. Neither is a full two-dimensional refinement count.
The pruning rule, monotone slope argument and shadow application are
supporting tools, not a new global bound or algorithmic-priority claim.

Complete exact line traces in one fixed cover cone give three states
for $(3,5,7)$ and thirteen for $(5,7,11)$, in both checked fields
$\mathbb Q$ and $\mathbb F_2$. Every latter state reaches all three
period thresholds. All 310 generic comparison intervals on that line
are covered; each distinct state has independent ordinary polynomial
and strict-border checks, and the cover is independently replayed.
Other directions of tie refinement are not counted. The prices here
are bounded in ratio, so this is not a separating sequence or growth
estimate. Equal finite counts in the two fields do not identify their
state sets or their full fans.

Hrubeš--Yehudayoff's inspected Theorem 4 and Lemma 5 give the
established line-breakpoint/shadow and greedy-comparison framework.
Sabariego--Santos's inspected triple-loop results concern coherent
minimum-distance diagrams at total-degree ties; they do not settle
this circuit quotient's positive-cost refinements. The exact current
priority and unrestricted extremal question remain unresolved. The
next proof must control the two-dimensional refinement geometry,
rather than infer it from these one-dimensional traces.


## Boundary geometry and the established linkage limit

The [periodic-cover continuation](periodic_quotient_arithmetic.md#the-exact-prior-art-boundary-for-box-residuals)
now identifies the precise primary source for the old box-residual
(RC1)--(RC3) reduction: Huneke--Ulrich's inspected Lemma 3.1 for
linkage through regular initial monomials. The equal odd-period
$O(n^{3/2})$ exclusion was already proved in (D10); neither fact is a
new count in this continuation.

A general downset boundary argument gives an actual sector-wide cube
candidate bound $3\sum_i|\pi_i(E)|$. Surviving base corners lie in
boundedly many bottom and top coordinate-fiber layers. The proof charges
each first exit of a fixed positive step path to one fiber boundary;
it requires no convexity or bound on the number of cover socles.
Projection sizes cannot simply be replaced by corank. The infinite
pairwise-coprime family $(2q+1,q,q+1)$ has explicit two-arm planar
covers of size $M=\Theta(q^3)$, corank $k=\Theta(q^2)$, and no
local cube exclusions, although every length-$n$ quotient standard set
inside them reaches all three pair-degree thresholds. These examples
are planar and their two-socle upper-shell bound survives. In general,
all states omitting a variable number only $O(n\log(n+1))$: restrict
to the corresponding colength-$n$ planar elimination quotient and apply
the existing all-field planar theorem, then sum over three variables. They refute
the cube-only shortcut, not the joint filter or a genuinely
three-variable count.

For every nonbox full group cover, degeneration already obstructs
transferring global residual fans: the initial annihilator of the
augmentation ideal has dimension one, whereas the annihilator of its
initial ideal is the monomial socle of dimension $\sigma>1$.
The exact defect is $\sigma-1$, over every field. This isolates a
structural failure of the established box linkage argument; it is not
an assertion that every residual or every nonbox cone has the same
fan complexity. Independent norm quotients, corner tests and shortest
paths support these proofs. The unrestricted extremal and priority
requirements remain unchanged and unachieved.


## A uniform upper-shell estimate without a socle factor

The [new shell proof](periodic_quotient_arithmetic.md#a-socle-free-upper-shell-bound)
gives, for every finite downset $E\subset\mathbb N_0^d$ of size $M$,

$$
|\{a\in E:|E\cap(a+\mathbb N_0^d)|\le k\}|
 =O_d\!\left(M^{(d-1)/d}k^{1/d}\log^{(d-1)/d}(M+1)\right).
$$

Maximal occupied dyadic cells supply $O_d(\log^{d-1}(M+1))$ outer
boxes with total volume $O_d(M)$. An adaptive integer exit step bounds
the shell inside a box of volume $V$ by $d k^{1/d}V^{(d-1)/d}$;
Hölder sums these budgets. This removes the unknown number of maximal
monomials from the previous shell bound. Dyadic-cross monomial algebras
have at least half their basis at multiplication rank at most $2^L$
while their length is $\Theta_d(2^L L^{d-1})$. Thus the displayed
logarithmic exponent is necessary for these $M,k$ powers, though the
estimate is not claimed sharp at every pair of parameters.

The same bound controls basis monomials of multiplication rank at most
$k$ in every finite quotient with a downset basis: its standard products
already give at least the upper-interval rank. For a fixed comparable
three-period cover, all omissions therefore lie in one common set of
$O(M^{8/9}\log^{2/3}(M+1))$ monomials throughout the cover cone,
independently of fields and prices. Given that set, this bounds greedy
independence-oracle calls; constructing the cover and shell remains a
separate cost. The logarithmic witness is a monomial quotient with one
initial ideal and is not a fan-count separation.

A staircase-block monomial family also has
$\Theta(M^{8/9})$ surviving **joint** shell/cube candidates at corank
$\Theta(M^{2/3})$, while containing all three variables. Thus a smaller
power cannot follow from arbitrary downset geometry alone. These
one-initial monomial quotients are not asserted periodic group covers;
the next improvement needs their actual tiling or label-cost constraints.

This is a uniform supporting theorem and an algebraic consequence, not
the requested unrestricted extremal advance. Related dyadic weights in
van der Hoeven--Lecerf and the general parametric-matroid arrangement
in Hausbrandt--Ruzika are explicitly compared. Neither supplies the
missing quotient-specific two-dimensional count. Exact priority and all
original success criteria remain open.

## Reproduction

From the repository root:

```bash
python3 -B scripts/run_extremal_initial_ideals_checks.py
python3 -B scripts/run_generic_multiaffine_checks.py
python3 -B scripts/run_bilinear_rectangle_checks.py
python3 -B scripts/run_rectangle_boundary_checks.py
python3 -B scripts/run_bilinear_cographic_checks.py
python3 -B scripts/run_cographic_structure_checks.py
python3 -B scripts/run_cyclotomic_price_checks.py
python3 -B scripts/run_graded_fan_checks.py
python3 -B scripts/run_periodic_quotient_checks.py
python3 -B scripts/run_additive_initial_ideals_checks.py
python3 -B scripts/run_planar_dyadic_bound_checks.py
python3 -B scripts/run_rational_shadow_bound_checks.py
python3 -B scripts/run_initial_cone_height_checks.py
python3 -B scripts/run_weighted_facet_incidence_checks.py
python3 -B scripts/run_algebraic_incidence_controls.py
python3 -B scripts/run_nested_chain_checks.py
python3 -B scripts/run_nested_prefix_checks.py
python3 -B scripts/run_rectangle_rank_checks.py
python3 -B scripts/run_profile_lattice_checks.py
./build.sh
python3 tests/test_smoke.py
```

The new checker uses only Python's standard library and existing research
helpers. Its report records the relevant source hashes, exact witnesses and
scope. These research checks do not extend the CLI auditor's certified range.
The core kernels, old manuscripts, earlier research and paused ML work are
preserved.
