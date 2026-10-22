# Four planar quotients isolate quadratic-coordinate coupling

Let \(A=K[x,y,z]/I\) have positive finite length \(n\), over any field,
and suppose \(x^2=\alpha x+\beta\) in \(A\), for scalars
\(\alpha,\beta\in K\). Write \(B=K[y,z]\), and set
\[
 L=B1,\qquad U=Bx,\qquad W=L\cap U,\qquad
 p=\dim L,\quad q=\dim U,\quad h=\dim W=p+q-n.             \tag{1}
\]
The whole positive initial fan has an exact description by
\(O(n\log(n+1))\) fixed represented matroids, each of rank \(h\)
on \(2h\) labelled columns. Each matroid applies on one vertical
strip in \((u,t)=(w_y/w_z,w_x/w_z)\).

Consequently, with absolute constants,
\[
 \boxed{\quad g(I)=O\bigl(n(1+\log(n+1))(1+h)^{10/3}\bigr).\quad}
                                                               \tag{2}
\]
For every fixed \(h\), the worst order is
\(\Theta_h(n\log(n+1))\), with characteristic-zero sharpness even
when \(h>0\). This is a uniform coupling-rank class theorem, not an
improvement of the unrestricted \(G_3\) exponent or a novelty certificate.
When \(h\) grows, retain the
[earlier \(O(n^{7/3})\) class bound](quadratic_coordinate_fan_compilation.md)
if it is smaller. No sharp dependence on growing \(h\) is claimed.

The reduction improves on using all internal orders of the two full
planar bases. Their four **unordered** initial bases suffice to fix
the compulsory columns. For a fixed algebra, constancy of those four
sets keeps the residual matrix fixed as well. The four planar systems
alone do not supply that matrix or determine the fan.

## An actual ideal carrying the overlap

The intersection \(W\) is an ideal of \(A\), not merely an independent
column span. It is stable under \(y,z\) because \(L,U\) are
\(B\)-modules. If \(w\in W\), write \(w=x\ell\) with \(\ell\in L\).
Then
\[
                xw=\alpha w+\beta\ell\in L,
 \qquad xw\in U\quad\hbox{because }w\in L.
\]
Thus \(xW\subseteq W\). The quadratic relation also gives
\(A=L+U\), proving the length identity in (1).

The four spaces
\[
                         L,\quad U,\quad L/W,\quad U/W       \tag{3}
\]
are cyclic \(B\)-modules, of lengths
\(p,q,p-h,q-h\). Each is isomorphic as a \(B\)-module to a planar
polynomial quotient. For \(U\), the cyclic vector is \(x\), not
the unit of \(A\). Zero modules are allowed. Their lengths sum to
exactly \(2n\).

The actual quotient \(\bar A=A/W\) has the direct vector decomposition
\[
                         \bar A=(L/W)\oplus(U/W).           \tag{4}
\]
Its two-layer monomial greedy selection therefore consists of the
independent planar greedy bases of these two summands. Adding \(t\)
to every upper-layer cost does not change either internal selection.
The fan of \(\bar A\) is independent of \(t\).

## Four planar fans fix a small matrix

At a generic planar ratio \(u>0\), denote the initial monomial bases
of (3), in the same order, by
\[
                  B_0(u),\quad B_1(u),\quad C_0(u),\quad C_1(u).
                                                               \tag{5}
\]
Upper labels retain their first exponent one. All four can be chosen
in their usual length-sized planar hyperbolic grounds.
Projection of columns cannot turn a previously dependent column into
an independent one. Hence the two greedy scans give
\[
                         C_i(u)\subseteq B_i(u).            \tag{6}
\]

The positive planar initial cones are intervals. Overlay the four fans,
using changes of basis **sets**, without subdividing at every internal
order equality. If \(K\) is the number of open overlay intervals,
the existing planar theorem gives
\[
 \begin{split}
 K&\le1+\sum_{M\text{ in }(3)}(g(M)-1)\\
  &\le1+4\sum_{r\in\{p,q,p-h,q-h\}}r\lceil\log_2\max(1,r)\rceil\\
  &\le1+8n\lceil\log_2\max(1,n)\rceil.                  \tag{7}
 \end{split}
\]
A zero module has one empty basis. The same is true for length one,
so the endpoint conventions in (7) cause no exception.

Fix one such interval \(J\). Its four sets in (5) are constant.
Put
\[
 C=C_0\sqcup C_1,\qquad
                 E=(B_0\setminus C_0)\sqcup(B_1\setminus C_1).
                                                               \tag{8}
\]
Then \(|C|=n-h\) and \(|E|=2h\). For every generic positive
objective with \(u\in J\), its original standard set satisfies
\[
                              C\subseteq S\subseteq C\sqcup E.
                                                               \tag{9}
\]
For the first inclusion, \(C\) is precisely the standard basis of
the actual quotient \(A/W\), by (4). Initial-ideal inclusion for
nested polynomial ideals puts every quotient standard monomial in
the original standard set. For the second inclusion, a monomial
rejected in its own planar layer is a combination of earlier columns
in that same layer, and is still rejected in any merged greedy scan.

Contract the independent columns \(C\) and retain \(E\). The result
is a fixed rank-\(h\) represented matroid \(N_J\). Each of the two
parts of \(E\) is a basis: a complement to \(C_0\) inside \(L\),
or to \(C_1\) inside \(U\), projects to a basis modulo \(\langle C\rangle\).
Indeed \(\langle C\rangle\) is a direct vector complement to \(W\).

Every basis of \(N_J\), lifted with \(C\), is an actual monomial
vector basis of \(A\). It cannot beat the global minimum. Conversely
(9) supplies an optimal lift at every generic positive objective in
the strip. Thus throughout the strip the minimum support functions
agree, up to the fixed moment of \(C\). Continuity gives agreement
on all walls as well. This proves equality of full cone patches,
not just a finite list of selected weights.

Although \(W\) is an algebra ideal, \(\langle C\rangle\) usually
is not. When \(0<h<n\), it is proper and contains \(1\).
The rank-\(h\) matroid is not asserted to be a length-\(h\) algebra
quotient or to have a length-\(h\) universal monomial ground. Its
labels can remain long. The matrix and the compulsory set can change
at the endpoints of \(J\); no global fixed contraction is inferred.

## Counting the residual choices

For a rank-\(h\) matroid on \(2h\) arbitrarily labelled points in
three-dimensional space, the full profile polytope has
\[
                              O((1+h)^{10/3})               \tag{10}
\]
vertices. Here is the standard-input derivation used in the earlier
[mixed-rank proof](mixed_rank_shadow_fans.md). A projected edge is
parallel to a difference of two labels, giving \(O(h^2)\) directions.
All edges in one direction inject into the vertices of the shadow
along that direction: a fiber over a projected vertex is one segment.
The published one-parameter matroid estimate bounds that shadow by
\(O((2h)h^{1/3})\). Summing gives \(O(h^{10/3})\) edges and hence
vertices. Polygons, segments, points and \(h=0\) are included
separately. Degeneracies do not increase exposed-profile counts.

The one-parameter input is
[Dey, Theorem 4.1 and its matroid consequence in Section 4.1](https://link.springer.com/content/pdf/10.1007/PL00009354.pdf).
The projected-edge step is the established argument in
[Onn--Rozenblit, proof of Theorem 2.3](https://arxiv.org/pdf/1208.5639).
This is not an asserted spatial \(k\)-set theorem for arbitrary matroids.

Every original initial cone meets at least one open strip. Its patch
there is a strict minimum of \(N_J\), and distinct patches with
distinct basis labels give distinct exposed profiles. Equations
(7)--(10) prove (2). More explicitly,
\[
 g(I)\le K\min\left\{\binom{2h}{h},\ C(1+h)^{10/3}\right\}.  \tag{11}
\]
For example \(h=1\) gives at most two patches per strip.
A polynomially separating quadratic-coordinate family would therefore
require
\[
 h=\Omega_{c,\epsilon}\left(
       n^{3(1+\epsilon)/10}/\log^{3/10}(n+1)\right)
       \quad\hbox{if }g(I)\ge c n^{2+\epsilon}.             \tag{12}
\]
This is a necessary condition, not a large-fan criterion.
It holds at every translated center \(x-c\), \(c\in K\): translation
changes monomials only by proper-divisor terms, preserving global
initial ideals and the quadratic-coordinate hypothesis.

## A common small weight list before choosing the coupling

The four labelled planar fans determine a common list of
\[
                   O\bigl(n(1+\log(n+1))(1+h)^4\bigr)       \tag{13}
\]
strict positive integer weights meeting every initial cone of every
compatible ideal having those same four fans. This list can be built
without the entries of the residual matrices. It samples possible
orderings; recovering the basis at a weight still requires its matrix.

On each strip, arrange all pair-comparison lines of the \(2h\) labels
in (8), rather than just the matched pairs for one matrix. There are
\(O(h^2)\) lines and \(O((1+h)^4)\) open chambers. Every ordering of
the residual columns is fixed on a chamber, for every representation,
so greedy choice is fixed there. One weight per chamber suffices.

The weights may be chosen with infinity norm at most \(6n^2\).
Anchor \(y,z\) coordinates are at most \(n-1\), their \(x\) coordinates
are at most one, and the strip boundaries
come from planar head-tail differences with coordinate magnitudes at
most \(n\). All defining integer normals therefore have infinity
norm at most \(n\), including the positivity boundaries. Primitive
extreme rays have height at most \(2n^2\), by the ordinary cross-product
cofactor bound. Sum any three linearly independent feasible primitive
rays of a chamber. The sum is strictly positive and strictly satisfies
every nonzero chamber inequality: equality would force all three rays
onto one plane. Its height is at most \(6n^2\).

This gives an exact geometric enumeration reduction with separate
field arithmetic. After preparing each fixed small matrix, every
sample needs only a rank-\(h\) greedy computation. The list depends
on the four labelled fans, not only their lengths; computing them,
preparing the matrices, coefficient bits and writing the output are
additional costs. No coefficient-independent list for arbitrary
length-\(n\) ideals is asserted.

## Sharp planar order at every fixed overlap rank

Fix \(h\ge1\) and \(m\ge1\), and work over \(\mathbb Q\). Set
\[
 G(T)=\prod_{j=2}^{m+1}(T-j),\qquad F(T)=(T-1)^{2h}G(T),
 \quad n=m+2h,
\]
and take the kernel of
\[
 \mathbb Q[x,y,z]\longrightarrow\mathbb Q[T]/F,
 \qquad x\mapsto(T-1)^hG(T),\quad y\mapsto T,\quad z\mapsto T^n.
                                                               \tag{14}
\]
This map is surjective and has length \(n\). It has \(x^2=0\),
and multiplication by \(x\) has rank exactly \(h\), by the gcd
with \(F\). Since \(L=A\), its intersection rank in (1) is exactly
the prescribed \(h\). Nonconstant multiples \(x(y-1)^j\) survive
for \(0\le j<h\).

On the planar universal ground \(V_n^2\), all exponents \(a+nb\)
are distinct, since \(0\le a<n\). Any \(n\) corresponding columns
in \(\mathbb Q[T]/F\) are independent. Otherwise a polynomial
with at most \(n\) nonzero monomials would vanish at the positive
nodes \(1,2,\ldots,m+1\), with multiplicity \(2h\) at one and
one at the others, giving at least \(n\) positive zeros counted
with multiplicity.

The classical few-term zero bound rules this out. For completeness,
divide by the lowest monomial and differentiate. Rolle's theorem,
counting multiplicities, lowers the zero count by at most one, while
differentiation removes the constant term. Induction bounds positive
zeros of a nonzero \(k\)-term polynomial by \(k-1\).

Thus the planar residue matroid is uniform on this finite ground.
It realizes every length-\(n\) planar corner cut. All these bases
also occur in the full fan: for each planar weight choose \(t\)
large enough to exclude upper columns. Consequently
\[
                       g(I)\ge C_2(n)=\Theta(n\log(n+1)). \tag{15}
\]
The planar corner-cut order is the established
Corteel--Rémond--Schaeffer--Thomas result, restated in
[Elizalde--Galván, Theorem 2.3](https://arxiv.org/html/2312.16353v1).
The proof uses that preprint statement under the existing working
assumption about its journal version. This does not establish priority
of the present coupling formulation. Combining (2) and (15) proves
the fixed-\(h\) order asserted above. For \(h=0\), a planar quotient
with \(x=0\) gives the same order.

## The matrix cannot be discarded

The [existing permutation pair](planar_layer_rank_obstruction.md) has
the same labelled planar column matroids for **all four** modules (3).
The lower and upper bases are its fixed monomial diagrams. The two
quotient diagrams simply delete all lower and upper identified tips,
independently of their matching permutation. Hence there is one strip
for both ideals, with the same compulsory columns and the same anchor
labels. Their rank-\(2k\) anchor matrices encode different identifications,
and their full fan counts are \(2\) and \((k+1)^2\).
This is retained as a correctness guard, not a new lower family.

The reduction therefore keeps compatible joint-rank information in
the fixed matrices. It does not reconstruct the fan from four separate
rank oracles, count carrier segments as cells, or contract a compulsory
span as an algebra. The unrestricted endpoint and the large-\(h\)
collective crossing problem remain open.

## Reproduction and priority boundary

~~~bash
python3 -B scripts/run_quadratic_intersection_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The report is `results/research/quadratic_intersection_checks.json`.
The checker derives all four quotient column maps, builds the fixed
small matrices and compares their reconstructed whole cones with a
separate full border-fan traversal. It checks coordinate commutation,
the quadratic relation, multiplication stability of \(W\), ordinary
S-pairs, finite confluent free-prefix ranks, and the equal-four-system
guard. A separate all-pair arrangement verifies the common sampler,
including its integer heights and its equality for the two different
couplings with identical planar inputs. Bounded controls support the
uniform proofs; they do not
computationally establish an asymptotic or priority.

Cyclic modules, intersections, initial-ideal inclusion, greedy matroid
contraction, planar counting and one-parameter shadow estimates are
established inputs. Two-list interleaving is already the subject of
[Gusfield (1984)](https://www.sciencedirect.com/science/article/pii/0166218X84900775);
only its publisher abstract, rather than its complete theorem text,
has been compared. The four-quotient compression and fixed-overlap
packaging are proved directly here, with exact priority unassessed.
The sharp unrestricted \(G_3\) question, a polynomial lower separation,
the growing-\(h\) collective crossing bound and exact priority remain open.
