# Generic points on multiaffine surfaces cannot give a power separation

This continuation screens a full realizable class, rather than fitting a
growth law to its first few fans. If $f\in K[x,y,z]$ is geometrically
irreducible and multiaffine (degree at most one in each variable), the
ideal of $n$ sufficiently general points on $V(f)$ satisfies

$$
g(I)=O(n^2\log(n+1)).                                      \tag{M1}
$$

The constant is absolute, independent of the coefficients and
characteristic. Here sufficiently general means a specified nonempty
Zariski-open condition on the ordered point configuration. Over a finite
field that open set need not have rational points; one may pass to an
extension. The assertion is not about arbitrary point sets or arbitrary
quotients containing $f$.

The proof is a corollary of classical convex-position $k$-set geometry,
combined with generic interpolation. It is not a new unrestricted bound
for $G_3$, a novelty certificate, or completion of the active goal. It
rules out polynomial separation from the generic corner-cut scale in
this whole class. The squarefree leading monomials and the generic
evaluation condition are essential to this argument.

## Generic evaluation truncates the carrier's monomial ranks

Let $X$ be a positive-dimensional geometrically integral affine variety,
and let $E$ be any finite collection of monomials. Over a nonempty open
subset of $X^n$, evaluation at $n$ points has, simultaneously for every
$A\subseteq E$, rank

$$
\operatorname{rank}\operatorname{ev}(A)
 =\min\bigl(n,\dim_K\langle A\rangle_{K[X]}\bigr).           \tag{M2}
$$

To see this, choose any $s\le n$ linearly independent functions. They
admit $s$ evaluation points with a nonsingular evaluation matrix: choose
the points successively, since a nonzero function cannot vanish on a
dense set. The corresponding determinant on $X^s$ is not identically
zero. The product of the finitely many required determinants is nonzero
on the integral variety $X^n$. Impose distinctness of the points as
well. Dependent functions give the reverse rank inequality automatically.
This is the usual generic interpolation argument, not a new matroid
truncation principle.

For a fixed global monomial order, take the standard monomials of the
carrier ideal $J=I(X)$. Their residue classes are independent. Each
nonstandard monomial reduces to earlier standard monomials. Thus the
standard set of the generic point ideal is exactly the first $n$
standard monomials of $J$ in that order, provided (M2) is imposed on
the finite relevant ground.

The ground can be chosen uniformly over all orders. Every such prefix
is a downset of size $n$, so an included exponent $a$ satisfies

$$
\prod_i(a_i+1)\le n.                                      \tag{M3}
$$

Consequently the canonical monomials in the hyperbolic cross $V_n$
suffice, for each of the finitely many carrier initial ideals. First
$n$ standard monomials always exist because $K[X]$ is infinite
dimensional. No infinite intersection of open conditions is used.
In particular, for a principal carrier the generic finite-point fan is
determined by its carrier heads and their prefix sets. For multiaffine
carriers with fixed monomial support it is independent of the nonzero
coefficients: the head regions and each $D_S$ are unchanged. Adjacent
carrier regions may merge when their finite prefixes agree.

## A classical bound for points on the orthant boundary

Let $E\subset\mathbb R_{\ge0}^3$ be a finite set of $m$ distinct points,
each with at least one zero coordinate. For $1\le r<m$, the number of
strict $r$-sets of $E$, even allowing all separating directions, is at
most

$$
2r(m-r)-m+2.                                               \tag{M4}
$$

For $r=m$ there is just one set. To justify using the classical convex
position formula, put

$$
\psi_\epsilon(p)=p+\epsilon\|p\|_2^2(1,1,1).
$$

If $p_i=0$, the functional $e_i-2\epsilon p$ evaluates on
$\psi_\epsilon(q)-\psi_\epsilon(p)$ as

$$
q_i+\epsilon\|q-p\|_2^2
 -2\epsilon^2\Bigl(\sum_jp_j\Bigr)
                       (\|q\|_2^2-\|p\|_2^2).             \tag{M5}
$$

For all distinct $p,q\in E$ this is positive when $\epsilon>0$ is
sufficiently small. Every transformed point is therefore an exposed
vertex. There are finitely many strict $r$-sets; choose one separating
functional for each. Their strict gaps persist for small $\epsilon$.
A further sufficiently small perturbation gives general position while
preserving these gaps and exposed vertices.

For $m$ points in convex position and general position in dimension
three, the classical number of oriented $j$-facets is

$$
e_j=2(j+1)(m-j-2),\qquad 0\le j\le m-3.                    \tag{M6}
$$

This is the formula stated in Leroux--Rademacher, Section 2.5 and
Proposition 2.13, PDF pages 6--7, with references to earlier work.
It is an established input, not Theorem 3.12's higher-dimensional
improvement. The complete relevant pages were inspected. Their
Proposition 2.11 also explains why perturbation preserves strict sets.
Source: [author-hosted published text](https://www.math.ucdavis.edu/~lrademac/algebraicksets.pdf),
[DOI](https://doi.org/10.1007/s00454-021-00340-1).

For completeness, the conversion to (M4) is elementary. The convex hull
of sums of $r$ distinct points has one vertex per strict $r$-set. Each
facet comes from a threshold containing exactly three ground points;
there are either $r-1$ or $r-2$ points strictly below it. Each facet is
a triangle, so Euler's formula gives

$$
a_r=2+\tfrac12(e_{r-1}+e_{r-2})=2r(m-r)-m+2,
$$

with nonexistent $e_{-1}$ and $e_{m-2}$ interpreted as zero. The
sum polytope is three dimensional when $m\ge4$ and $1\le r<m$:
its difference span contains all ground differences, by comparing
$r$-subsets sharing $r-1$ elements. The smaller cases are immediate.
This proves (M4) for the original possibly degenerate boundary set.

This step is for uniform subset selection. It cannot be applied to a
dependent quotient matroid merely because its nonzero monomials lie on
the boundary.

## Applying the bound to all multiaffine carriers

Every possible leading monomial of $f$ is squarefree and nonconstant.
Write it as $x_S=\prod_{i\in S}x_i$, with $t=|S|\in\{1,2,3\}$.
The carrier's standard exponents are

$$
D_S=\{a\in\mathbb N_0^3:\ a_i=0\text{ for some }i\in S\}.
$$

They lie on the orthant boundary. The relevant finite ground
$E_{S,n}=D_S\cap V_n$ has size

$$
m_t=t\sum_{j=1}^n\lfloor n/j\rfloor
      -\binom t2 n+\binom t3=O(n\log(n+1)).                \tag{M7}
$$

This follows by inclusion-exclusion on the coordinate planes, axes and
origin. Within the carrier cone with head $x_S$, (M2) makes each
standard set a strict $n$-set of $E_{S,n}$. For $n\ge2$ use (M4);
there are at most seven possible heads. Summing these bounds proves
(M1). Counting some sets outside their carrier cone only enlarges the
upper bound. For $n=1$ the point ideal has a single initial ideal.

There is also one list of $O(n^2\log(n+1))$ positive weights that
meets every initial cone of every generic point ideal in this class.
For each $S$, intersect its uniform-prefix fan with the fixed
arrangement comparing the eight squarefree monomials. That arrangement
has constant size. A prefix region is convex, so these intersections
create only a constant number of pieces per region. Choose one weight
in each nonempty open piece, and take the union over the seven $S$.
Any carrier support has a fixed head in such a piece. This is a
restricted common-list consequence of established state-fan methods,
not a new general enumeration algorithm or a bit-cost claim.

## An explicit realizable toric model

For $f=xyz-1$ the carrier has the same head $xyz$ for every global
order. Hence the complete point-ideal fan is the uniform-prefix fan of
$D_{\{1,2,3\}}$. Its count is independent of sufficiently general
point coordinates, the surface coefficients in the dominant-head
case, and the characteristic.

The required generic rank need not be left existential. For $n\ge2$
put $B=2n+1$ and, over $\mathbb Q$ with $T=2$, take the $n$ points

$$
p_i=(T^i,T^{Bi},T^{-(B+1)i}),\qquad 0\le i<n.               \tag{M8}
$$

For a canonical exponent $a\in D_{\{1,2,3\}}\cap V_n$, its
evaluation column is

$$
(1,T^{e(a)},\ldots,T^{(n-1)e(a)})^\mathsf T,\qquad
e(a)=a_1-a_3+B(a_2-a_3).                                  \tag{M9}
$$

The integer $e(a)$ is injective on this ground. All coordinates of
$a$ are at most $n-1$; the difference of two first Laurent coordinates
has absolute value at most $2(n-1)<B$. Equality in (M9) therefore
forces both Laurent coordinates to agree. Two canonical triples with
the same Laurent coordinates differ by a multiple of $(1,1,1)$,
which must be zero because their minimum coordinate is zero.

Every $n$ columns form an ordinary Vandermonde matrix in the distinct
nodes $T^{e(a)}$. Thus these explicit reduced rational point ideals
realize the entire fan. Over any characteristic use an indeterminate
$T$ over the original field; the Vandermonde determinant is still
nonzero. For numerical checks in a finite field, choose $T$ of order
greater than $2(B+1)(n-1)$.

Writing $(u,v)=(a_1-a_3,a_2-a_3)$ gives the exact alternative model

$$
w_1u+w_2v+(w_1+w_2+w_3)\max(0,-u,-v).
$$

Its sublevel set of height $R$ is the triangle with vertices
$(R/w_1,0),(0,R/w_2),(-R/w_3,-R/w_3)$. This explains the apparent
new three-way interleaving without supplying a new count theorem:
the classical convex-position restriction already prevents a power
separation.

## Scope controls and next decision

Special configurations on the same surface do not satisfy (M2).
For example, four distinct points $(t,t,t^{-2})$ satisfy $x=y$.
At weight $(2,3,5)$ the first four canonical monomials are
$1,x,y,x^2$ and have rank three. These points are reduced and all
monomial evaluations are nonzero. Containing a squarefree carrier
equation alone therefore does not justify uniform-prefix selection.
Fat-point conditions in the preserved local rectangle family are
another source of nonuniform ranks and remain outside (M1).

The checker below compares full small fans with independent comparison
arrangements, verifies the boundary perturbation and classical facet
counts, and constructs the point quotients by polynomial remainders
independently of diagonal evaluation. The Vandermonde controls include
symbolic determinant identities modulo characteristics two and three.
These finite checks support the proofs; no growth fit is used.

```bash
python3 -B scripts/run_generic_multiaffine_checks.py
```

The unrestricted endpoint and strongest-current-bound audit remain
open. The next lower-bound route must retain a prescribed nonuniform
monomial rank rule, or use a carrier beyond this squarefree class.
Generic interpolation on multiaffine surfaces is now excluded as a
source of a polynomial separation. The local low-order mixed rectangle
family and the independent quotient-sensitive upper route remain open.
