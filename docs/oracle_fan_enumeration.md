# An arithmetic consequence of the unrestricted counting bound

The rational-shadow theorem can be used to enumerate initial ideals,
rather than merely bound how many exist. This note gives an exact
reduction in the original coordinates, over every field. The polytope
reconstruction and greedy linear algebra are published tools. They
are attributed below; no new reconstruction algorithm or fastest-known
implementation is claimed.

Let $I\subset K[x,y,z]$ have positive colength $n$, presented by one
marked reduced Gröbner basis under a known global monomial order,
and let $g=g(I)$. Write
$L=1+\log(n+1)$. There is an algorithm that returns all initial
standard sets, the state polyhedron and all marked reduced Gröbner
bases using

\[
 O((g+1)n^3L^2)                                           \tag{O1}
\]

operations in $K$. Their union is a universal Gröbner basis; the
output may represent the union as a list of reduced bases with repeated
polynomials. No minimum-size or canonical duplicate-free output is
asserted. Separately, one can use

\[
 O((g+1)^2\log(g+2)+(g+1)nL^3)                            \tag{O2}
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
 \boxed{O\!\left(n^{47/8}L^{19/8}\right)}                  \tag{O3}
\]

field operations, with a smaller-order rational geometric-operation
budget. This is an arithmetic bound, not a bit-complexity, benchmark or
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

## Reconstruction and recovery of the bases

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

## Arithmetic and geometric costs

Compute all monomial columns on $U$ once. From one reduced basis,
[Babson--Onn--Thomas, Proposition 4.1](https://arxiv.org/pdf/math/0207135)
does this in $O(n^3L^2)$ field operations. From multiplication
matrices, generate columns in divisor order, using one dense
matrix-vector multiplication per monomial, with the same bound.

For each lexicographic query, sort $E$ and perform Gaussian
elimination on its $n$-entry columns, stopping at rank $n$.
This takes $O(n^2|U|)=O(n^3L^2)$ field operations. Column operations
or inversion of the selected $n$-column matrix also give every
border remainder in that bound. There are at most $3n$ minimal
border heads, all in $U$. For each one, output
$x^a-\sum_{b\in B}c_bx^b$. The strict monomial order makes the
head largest. These relations belong to $I$, and their leading
ideal has complement $B$ of size $n$. Inclusion and equal colength
prove that they generate $I$ and form its reduced Gröbner basis.
This is the established greedy/basis-conversion mechanism, with a
self-contained colength justification, not a new elimination method.

Equations (O6) and (O10) now prove (O1). For (O2), rebuild a
three-dimensional hull after each insertion in $O(k\log(k+2))$
geometric operations on $k=O(g+1)$ points, using a standard exact
hull routine. Summing gives $O((g+1)^2\log(g+2))$. The
[Chan 1996 algorithm, Sections 3--4](https://www.cs.jhu.edu/~misha/Spring20/Chan96.pdf)
provides the required bound and explicitly handles degeneracies;
perturbing nonvertices into vertices is unnecessary. Sorting the
$O(nL^2)$ costs at $O(g+1)$ queries contributes the other term.
All polyhedral input points and newly queried facet objectives have
polynomially bounded integer height in $n$. These operation counts still do not bound
the bit growth of coefficients in $K=\mathbb Q$.

## What the comparison establishes

| Inputs to the composition | Resulting field-operation bound |
|---|---|
| Babson--Onn--Thomas Theorem 4.2 itself, in three variables | $O(n^7L^{10})$ |
| Their classical $g=O(n^3)$ count plus established oracle reconstruction and greedy conversion | $O(n^6L^2)$ |
| The repository's unrestricted rational-shadow count | $O(n^{47/8}L^{19/8})$ |

The second row is already a consequence of published tools; the full
improvement from the first row must not be attributed to the shadow
theorem. Its additional saving is the passage from exponent six to
$47/8$. All three rows retain the originally specified coordinates.

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
