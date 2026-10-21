# Nilpotent carriers and the limit of high-order perturbations

The two-interaction rank theorem is a real distinction, but its
first-order truncation cannot give a polynomial separation above the
generic quadratic scale. A fixed quartic already belongs to every
ideal in that family. The resulting bound holds for arbitrarily many
higher terms, arbitrary coefficients and every field, subject to the
explicit pure-power condition below.

The geometry used here is classical. The bounded-coordinate layer
argument already occurs in the repository's
[relative-socle comparison](interleaving_research.md#a-uniform-relative-socle-bound).
The elementary power-membership identity was also used in the
[earlier hypersurface reduction](hypersurface_truncations.md).
The application to the proposed perturbation class is a feasibility
decision, not a new general lattice-vertex theorem or completion of
the unrestricted impact goal.

## A polynomial that bounds one coordinate in every order

Let $I\subset K[x,y,z]$ have positive colength $n$. Suppose it contains
a polynomial $f$ of total degree $m\ge1$ whose three pure coefficients
$[x^m]f,[y^m]f,[z^m]f$ are all nonzero. Let $S$ be any common bound
on the total degree moment $\sum_{a\in B}|a|$ of every initial
standard set $B$. Set

\[
 b=\min\{S,(m-1)n\}.
\]

Then, with an absolute constant,

\[
 \boxed{g(I)=O\!\left(\min\{(S+1)\sqrt{b+1},\,
                          (b+1)(S+1)^{2/3}\}\right).}       \tag{P1}
\]

This applies over finite fields and to nonreduced ideals. It requires
membership of $f$ in the ideal, not just vanishing on its reduced
zero set.

**Proof.** In any global monomial order, let $x_i$ be the largest
variable. Every monomial of degree at most $m$ is at most $x_i^m$,
by multiplicativity of the order. Hence the leading monomial of $f$
is $x_i^m$, and every standard monomial has $i$th exponent less
than $m$.

Group the initial standard sets according to one such largest
variable, and let $P_i$ be the convex hull of the moment vectors in
one group. Each moment is a distinct exposed vertex of $P_i$: a
strict positive weight exposing it in the full state polyhedron
also exposes it in this subcollection. This is the established
state-polyhedron dictionary of
[Babson--Onn--Thomas, Theorem 2.4](https://arxiv.org/pdf/math/0207135).
It is not necessary that the three groups be matroids, or that one
bounded coordinate work for the whole fan.

After permuting coordinates, $P_i$ lies in
$[0,b]\times[0,S]^2$. Adding the integer unit cube preserves at
least its number of vertices and gives a full-dimensional lattice
polytope of volume at most $(b+1)(S+1)^2$. The classical
three-dimensional lattice-vertex bound gives the first expression
in (P1), including when $P_i$ was flat.

For the second expression, group the **vertices** of $P_i$ by their
integer first coordinate. The vertices at one level are in convex
position in $[0,S]^2$: otherwise a vertex of $P_i$ would be a
convex combination of other vertices at the same level. The planar
lattice-vertex bound gives $O((S+1)^{2/3})$ of them. There are at
most $\lfloor b\rfloor+1$ levels. We use the hull of the vertices
at the level, not a potentially nonintegral slice of $P_i$.
Collinear and singleton level hulls have at most two vertices.
Summing the bounds for the three groups proves (P1). $\square$

Every $n$-element staircase has total degree moment at most
$\binom n2$: in a divisibility-compatible listing its $j$th element
has degree at most $j-1$. Thus, for $m\ge1$, (P1) implies

\[
 g(I)=O\!\left(\min\{\sqrt m\,n^{5/2},\,
                              m n^{7/3}\}\right).           \tag{P2}
\]

One can retain the smaller value of $b$ in (P1) for larger $m$.
For a local algebra supported at the origin the stronger common
budget is

\[
 S=\sum_j j\dim_K\mathfrak m^j/\mathfrak m^{j+1}.             \tag{P3}
\]

Indeed, in a filtration-adapted basis a monomial of degree $d$ has
zero coordinates in rows of grade below $d$. A nonzero determinant
matches its columns to rows of at least those grades. Summing proves
(P3) for every monomial basis, as in the
[existing Hilbert-budget proof](extremal_initial_ideals.md#a-uniform-bound-from-the-local-hilbert-function).

## A uniform bound for truncated hypersurfaces

Let $q$ be homogeneous of fixed degree $a\ge1$, with all three pure
degree-$a$ coefficients nonzero. Let every term of $h$ have ordinary
degree at least $D>a$, and put

\[
 I=(q+h,\mathfrak m^{s+1}),\qquad s\ge a,\qquad
 e=\max\{1,\lceil(s+1)/D\rceil\}.                          \tag{P4}
\]

There is no restriction on the number or coefficients of the terms
of $h$. Since $h^e\in\mathfrak m^{s+1}$ and $q=-h$ in the quotient,

\[
 q^e\in I.                                                \tag{P5}
\]

Its three pure degree-$ae$ coefficients remain nonzero in every
characteristic, even if binomial cross terms disappear. Thus it is
an admissible carrier for (P1), of degree $m=ae$.

The lowest homogeneous part of a product $u(q+h)$ is the lowest
part of $u$ times $q$. This gives the Hilbert function

\[
 h_j=\binom{j+2}{2}-\binom{j-a+2}{2}\quad(0\le j\le s),
                                                               \tag{P6}
\]

where the second term is zero for $j<a$. Its sum is
$n=\binom{s+3}{3}-\binom{s-a+3}{3}=\Theta_a(s^2)$, and
$S\le ns=O_a(n^{3/2})$. Substitution in the first bound of (P1)
proves

\[
 \boxed{g(I)=O_a(\sqrt e\,n^2).}                            \tag{P7}
\]

In particular, whenever $D\ge c s$ for a fixed $c>0$,

\[
 g(I)=O_{a,c}(n^2).                                       \tag{P8}
\]

This closes the whole class as a source of a polynomial gain over
the generic quadratic power. It is not confined to one or two
interaction monomials.

Conversely, a family in (P4) with
$g(I)\ge c_0 n^{2+\epsilon}$, for fixed $a,\epsilon,c_0>0$,
would need

\[
 e=\Omega_a(n^{2\epsilon}),\qquad
 D=O_a(n^{1/2-2\epsilon}).                                \tag{P9}
\]

The number of possible high substitutions must grow; a bounded
interaction depth cannot suffice. These are necessary conditions,
not a construction in the remaining range.

### The two-interaction family

For the independent and repeated equations in
[two_interaction_ranks.md](two_interaction_ranks.md), take

\[
 q=x^2+y^2+z^2,\quad h=Hz^r+Lz^{r+1},\quad
 s=2r+1,\quad D=r+2.
\]

Then $e=2$, $q^2\in I$, and $m=4$. The length is
$n=(s+1)^2$ and the exact local moment is

\[
 S=\frac{s(s+1)(4s+5)}6.
\]

Hence $g(I)=O(n^2)$ uniformly in odd $r\ge3$, all fields and all
binary quadratics $H,L$. The same argument works with arbitrarily
many terms starting at degree $r+2$. The earlier note permits
arbitrary $Q$ in place of $x^2+y^2$; this application requires both
pure coefficients of $Q$ to be nonzero. A missing pure coefficient
cannot silently be supplied by changing coordinates.

The exact rank gap and exceptional-characteristic ranks survive;
neither evades the quartic carrier bound. The checks give complete
fans at specified small parameters, but (P8) is proved without
fitting those counts.

## A coordinate-preprocessing corollary and its limit

Every length-$n$ three-variable ideal contains a nonzero polynomial
of degree at most $t$, where $\binom{t+3}{3}>n$; thus
$t=O(n^{1/3})$. Let $f$ have degree $m\le t$, and let $f_m$ be
its nonzero highest homogeneous part. For an invertible linear
substitution with columns $v_1,v_2,v_3$, the three pure coefficients
are $f_m(v_i)$. The polynomial

\[
 \det(v_1,v_2,v_3)\prod_i f_m(v_i)
\]

is nonzero. Its nonvanishing defines a nonempty open subset of
coordinate changes. Over an infinite field it has rational points;
over a finite field one may need a finite extension. Applying (P2)
after such a change proves

\[
 g(TI)=O(n^{8/3})                                         \tag{P10}
\]

for a nonempty open set of linear changes, over the extension when
needed. This uses only a dimension count, the pure coefficients,
and established lattice geometry. Its exact priority remains
unassessed. It does **not** bound $G_3(n)$ in the originally supplied
coordinates: changing coordinates can change the number of initial
ideals. It is not a justified replacement for the unrestricted
goal. No arithmetic-operation complexity claim is inferred here.

## Why counting normal-form entries is insufficient

A tempting stronger bound would count only hyperplanes comparing a
parent monomial with a nonzero entry of one of the three coordinate
normal forms ($x<2$, $y<2$, or $z<2$). This proposed refinement is
false already at $r=3$, over $\mathbb Q$.

Write $F=q+xyz^3+x^2z^4$ and truncate above degree seven. Direct
multiplication gives the relation

\[
 G=q(x+yz)+x^2yz^3-xz^6
    \equiv F(x+yz-xz^4)\pmod{\mathfrak m^8}.               \tag{P11}
\]

The complete fan has an active wall $X+Y=3Z$ on the positive
segment between $(1,2,1)$ and $(3/2,3/2,1)$. At prices
$(5/4,7/4,1)$ the two displayed high monomials tie strictly above
the other terms of $G$. Exact adjacent standard sets exchange
$xz^6$ with $x^2yz^3$, and the checker verifies both entire
adjacent cones, not just the polynomial's membership.

This wall is absent from all three entry lists. In a coordinate
normal form the other two exponents never decrease, and ordinary
degree never decreases. A parent-minus-entry vector therefore has
at most one positive coordinate and nonpositive coordinate sum.
Neither orientation of $(1,1,-3)$ has both properties. Thus
elimination creates actual walls beyond these entry comparisons.
The single example disproves that universal shortcut; it supplies
no large-fan lower bound.

## Comparison and reproducibility

The vertex exponents in (P1) are classical Andrews bounds, with
[Dillon, Theorem 1 and Sections 1--2.2](https://arxiv.org/pdf/2606.30856)
an inspected proof retelling of earlier work. The original
[Andrews article](https://doi.org/10.1090/S0002-9947-1963-0143105-7)
was not retrieved. The monomial-state dictionary is
[Babson--Onn--Thomas, Theorem 2.4 and Corollary 2.5](https://arxiv.org/pdf/math/0207135).
These sources do not constitute a novelty certificate for the
applications above. The stronger-current unrestricted comparison,
sharp growth and the original publication requirements remain open.

~~~bash
python3 -B scripts/run_power_carrier_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker replays the polynomial power witnesses by independent
convolution, tests reference kernels and moment caps, certifies
complete finite fans by exact cone traversal and total area, and
checks (P11) and both neighboring initial sets. The classical
asymptotic vertex theorem is an input, not something established
by these finite checks. No manuscript or priority claim is added.
