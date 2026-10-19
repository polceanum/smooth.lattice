# Weighted facet incidence and the limit of surface-area charging

The current spatial proof bounds large-normal vertices by surface area.
This note isolates the combinatorial quantity that a stronger proof must
control. An explicit all-field local quotient shows that even facets with
primitive normal height proportional to $n^2$ can retain area proportional
to $n^4$. Their vertex count can remain four. Thus a power improvement
cannot follow merely by reducing that surface-area tail for all quotients.

No stronger incidence inequality is proved here. The sharp $G_3(n)$
endpoint, its priority audit and publication assessment remain open.

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

For example, a new theorem $T(I)=O(n^3)$ would yield
$G_3(n)=O(n^{9/4}\log^{3/8}(n+1))$ by choosing
$h=n^{3/4}/\log^{3/8}(n+1)$. That incidence theorem is **unproved**.
More generally a bound $T(I)=O(n^\tau\log^\sigma(n+1))$, $\tau>1$,
would give the sum of (2) and
$O(n^{(5\tau+3)/8}\log^{(5\sigma+3)/8}(n+1))$.
Obtaining a suitable uniform budget requires a new argument about
simultaneous dependencies in one quotient.

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
