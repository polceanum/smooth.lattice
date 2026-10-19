# Weighted facet incidence and the limit of surface-area charging

The current spatial proof bounds large-normal vertices by surface area.
This note tests a proposed incidence budget for a stronger proof.
An explicit all-field local quotient shows that even facets with
primitive normal height proportional to $n^2$ can retain area proportional
to $n^4$. Their vertex count can remain four. Thus a power improvement
cannot follow merely by reducing that surface-area tail for all quotients.

A three-pole rational-curve construction below now proves that the
classical quartic weighted-incidence order is itself sharp for actual
finite quotients, even radical point ideals. The proposed cubic budget
is false. The sharp $G_3(n)$ endpoint, its priority audit and publication
assessment remain open.

## Three planes cover every boundary cell

Let $I\subset K[x,y,z]$ have colength $n\ge1$, and let $B$ be its
universal residue-matroid basis-profile polytope. Use minimum costs and
intersect its initial cones with the positive weight simplex. Call a cell
boundary-adjacent if its closure meets a coordinate face, including only
at a vertex.

The [classical cofactor bound](initial_cone_height.md) gives a strict
positive integer certificate $c$ for every initial cone, with
$\|c\|_\infty\le6n^2$. Put $Q=12n^2+1$ and

$$
q_1=(Q,-1,-1),\quad q_2=(-1,Q,-1),\quad q_3=(-1,-1,Q). \tag{1}
$$

If a nonzero cone-boundary vector $r$ has $r_i=0$, then
$q_i\cdot r<0$, whereas $q_i\cdot c>0$. A positive combination
of $r$ and $c$ lies on $q_i^\perp$ and remains strictly inside the
initial cone. Each boundary-adjacent cell is therefore visible in one
of the three rational shadows (1). The existing height-dependent shadow
theorem bounds their total by

$$
O(n^{5/3}\log(n+1)). \tag{2}
$$

This covers cells touching the boundary only at a point, not just cells
visible in an unperturbed coordinate projection. It uses the established
integer-certificate and rational-shadow results, without a new walk or
enumeration algorithm.

## The remaining weighted incidence

For each strictly positive primitive interior fan ray $w$, let $d(w)$ be
the number of initial cells incident to it. Define

$$
T(I)=\sum_{w\text{ an interior ray}}\|w\|_\infty\,d(w). \tag{3}
$$

For a full-dimensional $B$, these are its minimum facet normals and
$d(w)$ is the number of vertices of that facet. All incident vertices
are exposed by nearby positive objectives. A lower-dimensional $B$
has no cells whose closure avoids the simplex boundary: its normal
cones contain a nonzero lineality direction, along which a positive
objective can be moved until it reaches a coordinate face. Such cells
are already covered by (2).

At most $T(I)/h$ nonboundary cells have an incident ray of height
greater than $h$. Every remaining cell has three independent rays of
height at most $h$. Their sum is a strict integer objective, and the
previous short-relation argument makes the cell visible in a shadow
with primitive kernel height at most $\lceil3\sqrt h\rceil$.
The rational-shadow theorem therefore gives, uniformly for $h\ge1$,

$$
g(I)=O\!\left(n^{5/3}\log(n+1)
       +nh^{5/3}\log(n+1)+T(I)/h\right). \tag{4}
$$

This is a reduction to an incidence budget, not a proved improvement
to the extremal exponent. The classical lattice-triangle estimate gives
$T(I)\le6S(B+[0,1]^3)=O(n^4)$ and recovers the existing bound.
For strictly positive normals the cube adds a single vertex to each
supporting face, leaving its vertex count and area unchanged.

The earlier candidate $T(I)=O(n^3)$ would have yielded
$G_3(n)=O(n^{9/4}\log^{3/8}(n+1))$. It is **false**, by the
three-pole quotient below. More generally the formal consequence of a
bound $T(I)=O(n^\tau\log^\sigma(n+1))$, $\tau>1$, is the sum of
(2) and $O(n^{(5\tau+3)/8}\log^{(5\sigma+3)/8}(n+1))$.
But no uniform power $\tau<4$ is possible. A stronger count must
control cells more directly, retain information about the threshold,
or select projection planes without bounding this full first moment.

## Quartic area survives at quadratic normal height

For every integer $k\ge2$ and every field, consider

$$
J_k=(xy-x^{k+1},\ x^{2k+1},\ y^{2k},\ z-y^k). \tag{5}
$$

For $t>ks$ and $u>kt$, the displayed generators are a reduced
Gröbner basis with heads $xy,x^{2k+1},y^{2k},z$. The two
noncoprime-head S-pairs reduce to zero using the pure powers; the other
pairs have coprime heads. The standard monomials are

$$
\{x,x^2,\ldots,x^{2k}\}\ \cup\
\{1,y,\ldots,y^{2k-1}\},\qquad n=4k. \tag{6}
$$

The quotient is local: all three variables are nilpotent. Its monomial
residue classes are independent unit vectors. For $1\le a\le k$,
$x^{a+k}=x^ay$; these are exactly the $k$ nontrivial classes on the
$x$ branch. For $0\le b<k$, $y^{b+k}=y^bz$; these are the $k$
nontrivial classes on the $y$ branch. All other classes in (6) have
one representative. Mixed monomials involving $xz$ or $xy^2$ vanish,
as does $z^2$. These assertions also follow by direct division in the
displayed basis and hold in every characteristic.

The residue matroid has $2^{2k}$ bases on its nonzero universal ground.
Its profile polytope is exactly

$$
B_k=a_k+k[0,v_1]+k[0,v_2],\qquad
a_k=(k(2k+1),k(2k-1),0),
$$

$$
v_1=(-k,1,0),\qquad v_2=(0,-k,1). \tag{7}
$$

All possible profile sums are the $(k+1)^2$ grid points
$a_k+i v_1+j v_2$, $0\le i,j\le k$. The four vertices are the
four initial standard sets: each branch simultaneously chooses its
cheaper representative. Their inequalities are the independent choices
$t\gtrless ks$ and $u\gtrless kt$. All four combinations occur at
positive weights.

The cross product

$$
v_1\times v_2=w_k=(1,k,k^2) \tag{8}
$$

is primitive. Thus $B_k$ has area
$k^2\sqrt{1+k^2+k^4}=\Theta(n^4)$ and intrinsic lattice area
$k^2$. Its only interior fan ray is $w_k$, incident to four cells,
so $T(J_k)=4k^2=n^2/4$.

For the exact regularization $P_k=B_k+[0,1]^3$ used by the upper proof,
$B_k$ and $B_k+(1,1,1)$ are two facets, of primitive normals
$-w_k,w_k$. Every other facet normal comes from a coordinate face of
the cube or an edge of $B_k$ and an edge of the cube. Their primitive
heights are at most $k$. At threshold $H=k^2/2$, exactly the two
displayed facets have height greater than $H$. Consequently

$$
S_{>H}(P_k)=2k^2\sqrt{1+k^2+k^4}=\Theta(n^4),\qquad
\#\{\text{vertices incident to these facets}\}=8. \tag{9}
$$

Here $H=n^2/32$. There is therefore no uniform power saving over
$n^4$ for the large-normal surface-area tail, even at this quadratic
threshold and even for local commuting quotients. The usual area charge
at this threshold is quadratic in $n$, while the actual number of bad
vertices is eight. This does not disprove a better *vertex-incidence*
budget. Equation (3) preserves that possibility by counting the four
facet vertices rather than all their lattice area.

The same example rules out a separate proposed shortcut: its strictly
positive facet normal has coordinate ratio $k^2=\Theta(n^2)$.
A universal $O(n)$ aspect-ratio bound on those normals is false.
The classical cofactor argument bounds every such primitive ray height
by $2n^2$, and its minimum positive coordinate is at least one.
Thus the worst positive-facet aspect-ratio exponent is exactly two
along these lengths. This auxiliary sharpness is not the extremal count.

## Scope, attribution and verification

The universal-ground/state-polytope dictionary is established by
[Babson--Onn--Thomas](https://arxiv.org/pdf/math/0207135). The integer
cofactor bound, lattice triangle charge and positive normal-fan arguments
are classical. Equations (2)--(4) explicitly combine them with the
repository's rational-shadow theorem. The binomial quotient (5) is an
elementary construction within an already bounded subclass. Its exact
priority is unassessed, and no new general binomial theorem is claimed.
It excludes a proposed proof mechanism; it neither establishes an
extremal lower family nor supports a publication claim by itself.

Run from the repository root:

```bash
python3 -B scripts/run_weighted_facet_incidence_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker independently compares ordinary polynomial division with
commuting quotient matrices, checks all S-pairs in several fields,
enumerates every profile sum by residue-class convolution, and computes
every regularized facet using integer supporting-plane tests. Exact
cone clipping verifies all four initial cones. The boundary-plane
certificates and weighted incidences are also replayed on the four saved
nonlinear rectangle fans. Commands, witnesses and source hashes are saved
in `results/research/weighted_facet_incidence_checks.json`. These bounded
checks support the uniform proofs, not a fitted growth rate or novelty.

## Downset geometry does not give a cubic incidence budget

The proposed $T(I)=O(n^3)$ cannot be extended to all representable
matroids on the same hyperbolic ground, even if **every** basis selected
by positive weights is an $n$-element downset. The next section goes further: another residue matroid, coming from
an actual quotient, has exactly the same positive fan.

Let $n\ge4$, $r=n-1$, and use the ground

$$
E=\{0\}\cup\{j e_i:1\le j\le r,\ 1\le i\le3\}
   \subset V_n^3.
$$

Take a coloop at $0$ and the uniform matroid $U_{r,3r}$ on the other
points. All points of $V_n^3\setminus E$ are loops. This rank-$n$
matroid is representable over $\mathbb Q$: assign distinct rational
parameters $t_e$ to the axis points and columns
$(0,1,t_e,\ldots,t_e^{r-1})$, and assign $(1,0,\ldots,0)$ to $0$.
The usual Vandermonde determinant proves that every $r$ axis columns
are independent.

For positive weights a minimum basis consists of $0$ and the $r$
cheapest axis points. It therefore has the form

$$
D_{a,b,c}=\{0\}\cup\{j e_1:1\le j\le a\}
 \cup\{j e_2:1\le j\le b\}
 \cup\{j e_3:1\le j\le c\},\qquad a+b+c=r. \tag{10}
$$

These are downsets. Every allocation occurs: the weights
$w_i=1/(a_i+1/2)$ put every chosen axis point strictly below cost one
and every unchosen one above it. There are exactly
$\binom{n+1}{2}$ positive cells.

At a positive fan ray three axis points $A e_1,B e_2,C e_3$ tie at
the selection threshold. There can be only one tied point on each
axis. There are $A+B+C-3$ strictly cheaper axis points, so a
2-dimensional exposed profile face occurs exactly when

$$
A+B+C=n\quad\hbox{or}\quad A+B+C=n+1. \tag{11}
$$

The tied selection chooses two or one, respectively, of the three
points. Its exposed profile face is a triangle and has three incident
cells. The primitive ray is

$$
w_{A,B,C}=\frac{(BC,AC,AB)}{\gcd(BC,AC,AB)}. \tag{12}
$$

Different triples in (11) give different rays. Equality of rays makes
the triples proportional; two triples on the same level must then be
equal. Between the two levels the primitive proportional triple would
have coordinate sum dividing both $n$ and $n+1$, impossible for a
positive triple. There are exactly $(n-1)^2$ interior rays, and

$$
T(M_n)=3\!\sum_{\substack{A,B,C\ge1\\A+B+C\in\{n,n+1\}}}
       \frac{\max(AB,AC,BC)}{\gcd(AB,AC,BC)}. \tag{13}
$$

In particular $T(M_n)=O(n^4)$. It is also $\Omega(n^4)$ along
$n\equiv3\pmod6$. Here is an elementary counting proof, rather than a
fit to computed incidences. Restrict to triples with sum $n$, each
coordinate in $[n/4,n/2]$, and all three congruent to $1$ modulo $6$.
The corresponding triangle in the $(A,B)$ plane has area $n^2/32$;
counting its rows in the congruence lattice gives
$n^2/1152+O(n)$ triples. No pair has a common factor $2$ or $3$.
For each prime $p\ge5$, the triples in this triangle with $p$ dividing
one specified pair number

$$
\frac{n^2}{1152p^2}+O(n/p+1).
$$

This follows by counting the same triangle in a residue class modulo
$6p$; symmetry covers any of the three pairs. Summing the error over
$p\le n/2$ is $O(n\log(n+1))$. The union bound and

$$
3\sum_{p\ge5}p^{-2}\le3\sum_{j\ge5}j^{-2}
 \le3\left(\frac1{25}+\int_5^\infty x^{-2}\,dx\right)
 =\frac{18}{25}<1
$$

leave a positive proportion of pairwise coprime triples for every
sufficiently large $n\equiv3\pmod6$. Their denominator in (12) is
one and their height is at least $n^2/16$. Equation (13) therefore
proves the quartic lower bound.

This matroid cannot be the residue matroid of a length-$n$ quotient.
The monomial $xy$ is a loop. The columns $x,\ldots,x^{n-1}$ span all
axis columns, including $y$. Were these quotient residues, $xy=0$
would imply $y^2=0$ by multiplying an expression for $y$ in those
powers of $x$. But $y^2$ is a nonloop. This contradiction uses
multiplication, while every coherent basis in (10) is still a downset.

Thus the hyperbolic ground, representability, small ground size and
downset property of every positive coherent basis together do not
improve the quartic incidence exponent. The matroid itself is a negative control, rather than a quotient
realization. It is essential to distinguish the full residue matroid
from its positive fan: the following different quotient matroid does
realize that fan. Neither construction gives a polynomial separation
for $G_3(n)$. The initial-cone note gives a separate
[uniform control](initial_cone_height.md#the-hard-cone-can-be-entirely-interior):
removing boundary-adjacent cells does not improve the worst projection
radius. A smaller radius for every remaining cone is therefore not available
either.

## Quartic incidence is sharp for actual finite quotients

Let $F\in K[t]$ be monic of degree $n\ge4$, with $F(0)F(1)\ne0$.
Consider

$$
I_F=(xy-1,\ (x-1)z-1,\ F(x))\subset K[x,y,z]. \tag{14}
$$

Its quotient is $K[t]/(F)$, with $x=t$, $y=t^{-1}$ and
$z=(t-1)^{-1}$. It has length $n$ over any field, including when
$F$ has repeated factors. A uniform choice valid in every
characteristic is $F=t^n-t+1$.

**Theorem.** The positive initial fan of every such $I_F$ is exactly
the fan just computed for $M_n$. In particular

$$
g(I_F)=\binom{n+1}{2},\qquad T(I_F)=T(M_n),
$$

so $T(I_F)=\Theta(n^4)$ along $n\equiv3\pmod6$. These equalities do
not depend on the coefficients, factorization or multiplicities of $F$.

First, the three curve relations

$$
xy=1,\qquad xz=z+1,\qquad yz=z-y \tag{15}
$$

have mixed quadratic heads under every global monomial order. Each
tail monomial strictly divides its head, except for the constant,
which is also smaller. Repeated use of (15) expresses every mixed
monomial as a sum of a constant and pure powers of strictly smaller
positive cost. Thus no mixed monomial can be standard.

Second, every set $D_{a,b,c}$ in (10) is a quotient basis. Multiply
its rational functions by the unit $t^b(t-1)^c$ modulo $F$. The
resulting $n$ polynomials are

$$
t^{b+j}(t-1)^c\ (0\le j\le a),\quad
 t^{b-j}(t-1)^c\ (1\le j\le b),\quad
 t^b(t-1)^{c-j}\ (1\le j\le c). \tag{16}
$$

Their degrees are at most $a+b+c=n-1$. They are linearly independent:
before multiplication, a rational-function identity among a polynomial
and principal parts at the distinct poles $0$ and $1$ forces all
principal parts to vanish, and then the polynomial to vanish. This
argument holds in every characteristic. Hence (16) is a basis of
$K[t]_{<n}$, whose reduction modulo the monic $F$ is injective.
Multiplication by the displayed unit preserves independence.

For any generic positive cost, the $n-1$ cheapest nonconstant axis
points form one such allocation. They and $1$ are independent. A
proper prefix of this set is independent too, and every mixed
monomial is already in the span of cheaper pure powers by (15).
The greedy residue selection therefore chooses exactly those axis
points. No power above $n-1$ can be among the first $n-1$ axis points.
This proves equality of the complete positive fans, not just equality
of their cardinalities. Global monomial orders give no extra initial
ideals by the finite positive-weight representation.

The full residue matroid here is **different** from $M_n$: $xy$ is
nonzero, and $1$ need not be a coloop. For $F=t^n-t+1$, for example,
$y=1-x^{n-1}$. Multiplication is fully consistent. The argument that
excluded $M_n$ cannot exclude this realization of its positive fan.

This sharpness holds already for radical point ideals over
$\mathbb Q$: take $F=\prod_{i=2}^{n+1}(t-i)$. The points are
$(i,1/i,1/(i-1))$. There is no issue of pretending an arbitrary matroid
is realizable.

The obstruction also survives the removal of boundary cells. In the
central coprime-triple family used above, all three cells incident to
each ray have positive counts on every axis. Their closures avoid the
coordinate boundary: if $w_i=0$, the next excluded point on axis $i$
would have zero cost, forcing all the selected axis endpoints to have
zero cost too, and hence $w=0$. Each such cell has at most six facets,
since its threshold comparisons are only between the next excluded
point on one axis and the last included point on another. Therefore
$\Omega(n^2)$ distinct wholly interior cells are incident to rays of
height at least $n^2/16$.

At threshold $h=n^2/32$, the number of such bad cells is consequently
$\Theta(n^2)$, which is the order of $n^4/h$. Thus a bound
$O(n^{4-\epsilon}/h)$ cannot hold uniformly over all thresholds and
quotients for any $\epsilon>0$. This does **not** disprove a sharper
tail estimate at the smaller threshold used to optimize (4). Nor does
it rule out a carefully selected family of shadows. It closes the
proposed power saving in the full incidence moment, even with
multiplication compatibility and with only interior cells counted.

### Every positive shadow trace has only linearly many cells

There is a second obstruction to the current projection strategy. Any
plane through the weight-space origin meets at most $2n-1$ open
positive cells of this same quotient. Consequently every collection of
such planes meeting all its cells has at least

$$
\left\lceil\frac{\binom{n+1}{2}}{2n-1}\right\rceil=\Omega(n)
\tag{17}
$$

members, even if real coefficients and arbitrarily large arithmetic
height are allowed. A polylogarithmic selected-plane cover for every
finite quotient is therefore impossible.

To prove the trace bound, parametrize an open line segment of positive
weights by $w_i(t)=u_i+t v_i$. The order of the three relative rates
$\rho_i(t)=v_i/w_i(t)$ is constant: the sign of each difference is
that of $v_i u_j-v_j u_i$. Rank the types by these rates as $0,1,2$.
At a threshold exchange between points $A e_i$ and $B e_j$, their
costs are equal. The derivative of the difference of costs has the
sign of $\rho_i-\rho_j$. Thus the leaving point always has a strictly
higher type rank than the entering point. The integer potential
$\sum_i\operatorname{rank}(i)a_i$ of the selected allocation strictly
decreases at every changing parameter and lies in $[0,2(n-1)]$.
Simultaneous exchanges can be paired in decreasing rate order and
have the same property. There are at most $2n-2$ changes and hence
at most $2n-1$ cells. Tied rates can be perturbed; every strict cell
previously intersected persists under a sufficiently small line
perturbation.

A plane meets the positive weight simplex in a line. A compact open
subsegment contains witnesses for every strict cell it meets, so the
same bound applies to the entire plane. Equation (17) follows from
the exact cell count. This is a lower bound for methods that enumerate
a complete fan solely by collecting planar shadows. It is not a lower
bound for arbitrary Gröbner-basis enumeration algorithms or for
smooth-number unranking.

The arithmetic shadow envelope has a further explicit limitation.
For a selected set $\mathcal Q$ of distinct primitive integer kernels,
its current summed bound is proportional to

$$
\mathcal E_n(\mathcal Q)
 =n\log(n+1)\sum_{q\in\mathcal Q}\|q\|_1^{1/3}.
$$

At most $O(R^3)$ primitive kernels, up to sign, have height at most
$R$. Any $k$ distinct kernels therefore satisfy
$\sum\|q\|_1^{1/3}=\Omega(k^{10/9})$: at least half have height
$\Omega(k^{1/3})$. Equation (17) implies, for every rational-plane
cover of this family,

$$
\mathcal E_n(\mathcal Q)=\Omega(n^{19/9}\log(n+1)). \tag{18}
$$

Thus choosing fewer or shorter integer planes, while charging **each**
by that same envelope and adding all charges, cannot prove the generic
$n^2$ power up to logarithms for the whole class. This is a limitation
of the stated proof budget, not a lower bound on actual shadow sizes:
the positive traces here have the linear bound just proved. No linear
bound on a full shadow outside positive objectives is asserted. It can be
avoided by quotient-sensitive trace bounds, overlap accounting or a
different global counting argument. It does not preclude improving the
current nonsharp bound; it prevents mistaking that one proof mechanism
for a route all the way to the sharp target.

Partial fractions, the basis argument and rational-curve interpolation
are classical tools. The exact fan and incidence sharpness here are
proved directly; their precise priority has not been established.
Since this family has only quadratic fan growth, it neither closes the
unrestricted extremal problem nor supplies the required polynomial
separation from generic corner cuts. No new paper is prepared on the
strength of this obstruction alone.

### A related source counts a different local star

[Sabariego--Santos, *Triple-loop networks with arbitrarily many minimum
distance diagrams*, arXiv:0705.3631v1](https://arxiv.org/pdf/0705.3631),
Section 3 and Corollary 4.4 with Lemmas 4.2--4.3, relate coherent
minimum-distance diagrams to initial ideals of finite cyclic lattice
ideals. Their diagrams first minimize ordinary path length, then break
ties. They therefore concern cells incident to the ordinary-degree ray,
rather than all positive weight cells. Their Hilbert-basis bound uses
the lattice of homogeneous routing relations. It neither states a
bound on (3) summed over all rays nor covers arbitrary non-binomial
quotients. This scoped comparison identifies relevant published
multiplication structure; it does not establish novelty of the present
controls or the strongest current unrestricted bound.

The additional controls are reproduced by

```bash
python3 -B scripts/run_algebraic_incidence_controls.py
```

Their report records the complete positive matroid fans, independent
Vandermonde-column selections, arithmetic counting fixtures, actual
three-pole quotient and point-ideal fans, complete positive-plane
traces, and the nine-cone controls. The uniform proofs above establish the
obstructions. Computation establishes neither the desired improved
quotient bound nor priority.
