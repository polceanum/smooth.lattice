# From arithmetic shadows to vertex counts in fixed dimension

The [hyperbolic matroid shadow theorem](rational_shadow_bound.md) applies
in any fixed ambient dimension. Counting rational **planes**, rather than
their bases or all their normal vectors, extends its vertex consequence.
For every field and fixed \(d\ge2\), it gives

\[
 G_d(n)=O_d\bigl(n^{e_d}(1+\log(n+1))^{b_d}\bigr),\qquad
 b_d=\frac{3(d-1)}{(d+1)(3d-5)},\quad
 e_d=\frac{2d(d-1)}{d+1}-\frac{d-1}{(d+1)(3d-5)}.        \tag{T1}
\]

This reproduces the existing \(23/8\) bound in three variables; it
does **not** improve that target. In four variables it gives
\(G_4(n)=O(n^{33/7}\log^{9/35}(n+1))\), improving the previous
repository exponent \(52/11\) by \(1/77\). In dimensions at least
five the earlier Dey-based bound in
[the lattice-shadow argument](interleaving_research.md) is stronger;
retain the minimum of the two bounds. The new transfer is useful as a
general geometric statement, not a claim that every dimensional
specialization is a new best bound. Exact priority remains unestablished,
and the original sharp-growth or separating-family goal remains open.

## Rational-plane height and two classical lattice facts

For a rational two-dimensional subspace \(E\subset\mathbb R^d\),
let \(\Lambda_E=E\cap\mathbb Z^d\) and
\(H(E)=\operatorname{covol}(\Lambda_E)\). For an integer basis
\(u,v\), its primitive Plücker coordinates are
\(u_iv_j-u_jv_i\), with gcd one, and

\[
 H(E)^2=\sum_{i<j}(u_iv_j-u_jv_i)^2.                     \tag{T2}
\]

Both facts below are classical geometry of numbers. Their short proofs
specify exactly what the application needs.

**Small containing plane.** If \(w\in\mathbb Z^d\) is primitive and
nonzero, \(d\ge3\), there is a rational plane containing \(w\) with

\[
 H(E)\le c_d\|w\|_2^{(d-2)/(d-1)}.                     \tag{T3}
\]

Extend \(w\) to a unimodular integer basis. Orthogonally projecting the
other basis vectors onto \(w^\perp\) gives a lattice of covolume
\(1/\|w\|_2\): the original determinant is its covolume multiplied
by \(\|w\|_2\). Minkowski's first theorem gives a nonzero projected
vector \(z\) of length at most
\(2\sqrt{d-1}\|w\|_2^{-1/(d-1)}\). Lift it to an integer \(v\).
The plane spanned by \(w,v\) has saturated lattice covolume no larger
than \(\|w\wedge v\|_2=\|w\|_2\|z\|_2\), proving (T3).
Saturation can only decrease this covolume. No bound on the size of the
integer lift is needed.

**Number of planes.** For \(T\ge1\), there are \(O_d(T^d)\) rational
planes with \(H(E)\le T\). To see the needed upper bound directly,
choose a Gauss-reduced basis \(u,v\) of \(\Lambda_E\):
\(\|u\|\le\|v\|\) and \(|u\cdot v|\le\|u\|^2/2\).
Their angle has sine at least \(\sqrt3/2\), so
\(\|u\|\ll\sqrt T\) and \(\|v\|\ll T/\|u\|\).
In a dyadic range \(U\le\|u\|<2U\), there are \(O_d(U^d)\)
choices of \(u\). After reduction modulo \(u\), the vector \(v\)
lies in a cylinder of longitudinal length \(O(U)\) and transverse
radius \(O(T/U)\). Both widths are bounded below by positive
dimensional constants. Covering integer points by their unit cubes
bounds their number by \(O_d(U(T/U)^{d-1})\), including oblique
cylinders. Thus that range contributes at most
\(O_d(T^{d-1}U^2)\) basis pairs. Summing the geometric series through
\(U=O(\sqrt T)\) gives \(O_d(T^d)\). Overcounting bases is harmless;
it does not introduce a logarithm.

The height exponent in (T3) cannot be uniformly reduced. For
\(w=(1,t,\ldots,t^{d-1})\), every containing rational plane has
primitive Plücker infinity height at least \(t^{d-2}\), and one has
exactly that height. Here is an elementary verification. Eliminate the
first coordinate of a second integer vector using \(w_1=1\), and
divide the remaining vector \(a=(0,a_2,\ldots,a_d)\) by its gcd.
The minors of \(w,a\) are now primitive, since their first row is
\(a_2,\ldots,a_d\). If \(a_{j+1}-ta_j\ne0\), choose its largest
index \(j\ge2\). Then
\(w_ja_d-w_da_j=t^{d-2}(a_{j+1}-ta_j)\).
If no such index exists, \(|a_d|=t^{d-2}|a_2|\ge t^{d-2}\).
This proves the lower bound. The plane spanned by \(e_1\) and
\((0,1,t,\ldots,t^{d-2})\) proves the upper bound. This sharpness
observation concerns containing one integer objective, not the number
of planes needed by a whole quotient fan.

## A transfer theorem for lattice polytopes

Let \(P\subset\mathbb R^d\), \(d\ge3\), be a full-dimensional
lattice polytope with surface measure \(S\). Suppose every rational
plane \(E\) satisfies

\[
 f_0(\pi_E P)\le A(1+H(E)^\theta),\qquad A\ge1,\ \theta\ge0.
                                                               \tag{T4}
\]

Then, for every \(h\ge1\),

\[
 f_0(P)\le C_{d,\theta}\left(
 A h^{(d-2)(d+\theta)/(d-1)}+S/h\right).                 \tag{T5}
\]

For completeness, triangulate a facet \(F\), using all its \(k_F\)
vertices. It has at least \(k_F-d+1\) top-dimensional simplices, each
of measure at least \(\|z_F\|_2/(d-1)!\), where \(z_F\) is its
primitive integer normal. Therefore

\[
 k_F\le d!\operatorname{area}(F)/\|z_F\|_2.             \tag{T6}
\]

The vertices incident to a facet with \(\|z_F\|_2>h\) number at
most \(d!S/h\). At every other vertex choose \(d\) independent
incident outward facet normals. Their sum is an interior normal
objective of length at most \(dh\), and after division by its gcd
it is a primitive integer objective. By (T3) a rational plane of height
\(O_d(h^{(d-2)/(d-1)})\) contains it. The vertex is strictly exposed
in that plane and hence remains a vertex of the orthogonal projection.
Count those planes by the preceding lattice argument, then sum (T4).
This gives the first term of (T5). The sum may count a vertex many
times; it remains a valid upper bound. This proves (T5) without any
generic-position, simplicity or field assumption.

## Hyperbolic matroid profiles and finite quotients

Let a rank-\(r\ge1\) matroid have weights in \(\mathbb N_0^d\)
satisfying \(\prod_i(a_i+1)\le H\), with \(H\ge1\), and let
\(B\) be its basis-profile polytope. No representability is needed.
Regularize it by \(P=B+[0,1]^d\); the normal fan refines that of
\(B\), so \(f_0(B)\le f_0(P)\).
Its surface measure is \(O_d((1+rH)^{d-1})\). Indeed it lies in a
cube of side \(1+r(H-1)\), and summing the coordinate projections
of facets gives at most twice the sum of the \(d\) projected volumes.

Choose an integer basis of \(\Lambda_E\) as the two rows of a
projection map. Its minors have gcd one, and the normalized minor sum
in the hyperbolic shadow theorem satisfies

\[
 \Delta\le\sqrt{\binom d2}H(E).
\]

The projection differs from orthogonal projection onto \(E\) by an
invertible map. A projected cube has at most \(2d\) vertices, so
the theorem gives (T4) with

\[
 \theta=1/3,\qquad
 A=C_d(1+\log(H+1))\bigl(1+r+(r^2H)^{1/3}\bigr).         \tag{T7}
\]

Optimizing (T5), including the cutoff \(h\ge1\), gives the reusable
matroid bound

\[
 f_0(B)\le C_d\left[
 A+A^{b_d}(1+rH)^{(d-1)(1-b_d)}\right].                  \tag{T8}
\]

To check the cutoff, put \(S_0=C_d(1+rH)^{d-1}\).
If \(S_0\ge A\), choose \(h=(S_0/A)^{b_d}\); otherwise
\(h=1\) gives \(O(A)\). This avoids using a fractional cutoff below
one. Formula (T8) is uniform over arbitrary matroids, including repeated
weights and lower-dimensional profile polytopes.

For an ideal of colength \(n\), use its rank-\(n\) residue matroid
on \(V_n^d=\{a:\prod_i(a_i+1)\le n\}\). Every initial standard
set is a vertex label in its basis-profile polytope. Take \(r=H=n\)
in (T7)--(T8); \(A=O_d(n\log(n+1))\). The exponent is

\[
 2(d-1)-(2d-3)b_d=e_d,
\]

which proves (T1) for \(d\ge3\). In dimension two the existing
unprojected shadow theorem directly gives \(O(n\log(n+1))\), also
the value of (T1). The usual union of reduced bases further gives
cardinality \(O_d(n^{e_d+(d-1)/d}\log^{b_d}(n+1))\);
the union argument itself is established and is not a new algorithm.

## Exact comparison and limits

The first column below is
[Babson--Onn--Thomas, Corollary 2.5 and proof](https://arxiv.org/pdf/math/0207135).
Its state dictionary is also the input used above. The middle column is
the previously proved Dey-based repository bound, not a certified
strongest-current published theorem.

| Dimension | Classical power | Previous Dey-based power | Power in (T1) | Smaller repository power |
|---|---:|---:|---:|---:|
| 3 | 3 | 44/15 | 23/8 | 23/8, already established |
| 4 | 24/5 | 52/11 | 33/7 | 33/7 |
| 5 | 20/3 | 376/57 | 33/5 | 376/57 |
| 6 | 60/7 | 740/87 | 775/91 | 740/87 |

For all \(d\ge5\) the older power is strictly smaller. Subtracting
the powers gives
\[
 E_d^{\rm old}-e_d=
 -\frac{(d-1)(3d^2-19d+23)}
 {3(d+1)(3d-5)(d^2-d-1)},
\]
whose numerator changes sign between four and five. Do not overwrite
the stronger older higher-dimensional bounds with (T1).

Primitive-plane counting is established, much more precisely, in the
classical work discussed by
[Horesh--Karasik (2023), introduction and Theorem 1.1](https://academic.oup.com/qjmath/article/74/4/1253/7191043).
Only that counting statement, its definition of height and the
introductory comparison were inspected; no full-proof audit of their
equidistribution result is claimed. The elementary upper-bound proof
above suffices here. Minkowski's theorem and unimodular completion are
classical as well. The step requiring exact priority assessment is the
combination (T4)--(T8) with the hyperbolic shadow bound.

[Hausbrandt--Ruzika, arXiv:2503.08178v2](https://arxiv.org/html/2503.08178v2),
Theorem 3.2 with proof and Corollary 3.4, still use the full pairwise
comparison arrangement. The latter gives \(O(m^{2p})\) parameter
regions. It is a positive comparison, not an exhaustive audit of all
later bounds. The unavailable concave-surface journal text is not
assumed to strengthen its preprint to general two-parameter matroids.

The main three-variable bound, the absence of a polynomial separating
family, and all publication gates remain unchanged. This transfer
does not repair an individual product-plane hypothesis: its planes
are indexed by covolume and exposed integer objectives. The next main
task remains a collective quotient-dependent count or a realizable
lower construction escaping the already excluded classes.

There is a useful limit on this proof budget. In dimension three,
substituting an envelope \(A=O(n\log^c n)\) and height power
\(\theta\ge0\) into (T5), with the unchanged \(S=O(n^4)\) budget,
gives power \(4-6/(5+\theta)\). Even the hypothetical height-independent
linear envelope would give \(14/5\), still above the generic power two.
This is an algebraic limit of that substitution, not a lower bound on
\(G_3\), a proof of the hypothetical envelope, or an optimality theorem
for all projection methods. Sharper overlap or quotient-dependent
charges remain outside it. Thus improving only this envelope cannot
by itself reach the original sharp-growth target through (T5).

## Reproduction

```bash
python3 -B scripts/run_rational_plane_transfer_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker independently enumerates small rational planes by primitive
Plücker coordinates, checks reduced lattice bases and exact containing
planes, verifies the infinite sharpness identity, and compares projected
matroid profiles against independent determinant enumeration. Its
dimension and exponent calculations are exact. These checks support
the proof, without establishing novelty or completing the active goal.
