# Ordered cographic bases and sparse periodic tree optimization

The cographic model has a broader, characteristic-independent realization.
For equal periods, monomial downsets restrict its trees sharply, and its
specific weighted graph admits a linear-size candidate set. These results
give a compressed initial-standard-set algorithm and a more defined
counting problem. They do not yet bound the whole fan sharply or establish
the active goal's polynomial separation.

## Known duality and the polynomial realization

Put $S_m(t)=1+t+\cdots+t^{m-1}$. For any field and integers $r,s\ge3$,

$$
J_{r,s}=(S_r(X),S_s(Y),Z-XY),\qquad n=(r-1)(s-1).         \tag{D1}
$$

Since $(X-1)S_r(X)=X^r-1$, all columns have periodic labels
$((a+c)\bmod r,(b+c)\bmod s)$. The only relations on one complete
periodic rectangle are its row and column sums: their span is an ideal
of the group algebra and has dimension $r+s-1$. The quotient columns
therefore represent $M^*(K_{r,s})$ over every field. In characteristic
$p$, taking $r,s$ to be powers of $p$ recovers the earlier local
rectangles after translation.

**This matroid duality is prior art.** Martin--Reiner,
[Theorem 1, its proof and Remark 6](https://www-users.cse.umn.edu/~reiner/Papers/cyclotomic.pdf),
give the cyclotomic prime-factor construction. Hidaka--Itoh,
[Theorem 1.2 and Section 5.2, including the following remark](https://arxiv.org/html/2403.10817v3#S5.SS2),
identify the two-circuit tensor product and its network-matrix structure.
Fu--Ma--Wang,
[Proposition 3.5, equation (3.4) and proof](https://arxiv.org/html/2609.29227v1#S3.SS2),
explicitly give the all-field complete-bipartite cographic realization.
Those specified passages were inspected as primary text. The full papers
and the novelty of the additional ordered assertions below have not been
audited. The graph-rank dictionary itself cannot support a novelty claim.

For distinct prime periods over $\mathbb Q$, (D1) gives the known
cyclotomic-field representation via the Chinese remainder map. For
arbitrary periods it is the two-circuit tensor product. These are
established arithmetic and matroid identifications; the additional
question concerns its monomial-weight orderings.

For the remainder take $r=s=q\ge3$. A strict positive weight is
$w=(\alpha,\beta,\gamma)$. If $\gamma>\alpha+\beta$, the standard
set is the ordinary $(q-1)$-by-$(q-1)$ $X,Y$ rectangle. Assume
$\gamma<\alpha+\beta$ otherwise. Every standard monomial then belongs to

$$
\mathcal D_q=\{X^aY^bZ^c:ab=0,\quad
0\le a,b\le q-2,\quad0\le c<q\}.                        \tag{D2}
$$

The heads $XY,X^{q-1},Y^{q-1},Z^q$ explain these exclusions.

## Branching is confined to two fixed vertices

**Theorem 1.** Let $B\subset\mathcal D_q$ be any monomial downset
that is a quotient basis, without assuming coherence. Its omitted label
set $T$ is a spanning tree of $K_{q,q}$. The $Z$-axis of $B$ has length
$q-1$ or $q$. Every tree vertex except row $q-1$ and column $q-1$
has degree at most two.

**Proof.** Write $h$ for the $Z$-axis length. If $h\le q-2$, no
selected label can have both coordinates at least $h$: every monomial
in (D2) has one label coordinate equal to its $Z$ exponent. Thus $T$
contains $K_{q-h,q-h}$, including a four-cycle. This contradicts the
cographic basis criterion. Hence $h\ge q-1$.

For $0\le u<h$, let $b_u$ be the largest $b$ for which $Y^bZ^u$
belongs to $B$, including $b=0$. These integers are nonincreasing and
lie between $0$ and $q-2$; put $b_u=-1$ for $u\ge h$.
Define $a_v$ similarly using $X^aZ^v$.

Compare tree-neighbor sets of consecutive rows $u,u+1$, where
$u\le q-2$. The $Y$-arm's selected labels in row $u$ are the cyclic
interval $u,u+1,\ldots,u+b_u$. Moving to the next row can introduce
at most the single new selected column

$$
v\equiv u+b_u+1\pmod q.                                 \tag{D3}
$$

For each fixed column $v$, the $X$-arm index is $(u-v)\bmod q$.
Incrementing $u$ cannot introduce a selected $X$-arm label, except
at the wrap $v=u+1$. If this was a tree edge in row $u$, then
$b_u=0$; otherwise the old $Y$-arm already selected it. That exception
is therefore the same column (D3). Consequently

$$
|N_T(u)\setminus N_T(u+1)|\le1.                          \tag{D4}
$$

Two rows of a tree share at most one neighbor, since two shared neighbors
give a four-cycle. Equation (D4) implies $\deg_T(u)\le2$. Interchanging
$X,Y$ gives the column assertion. This proves the theorem over every
field, using only the established cographic criterion. $\square$

Deleting the two exceptional vertices leaves a disjoint union of paths
and isolated vertices. This is a necessary structure of **every**
downset basis. It is not a claim that every such tree is a coherent basis
or that there are only polynomially many such downsets.

## A cycle certificate leaves at most $19q$ candidate edges

Let $\delta=\alpha+\beta-\gamma>0$, and write $u,v\in\{0,\ldots,q-1\}$.
The least representative cost has the exact formula

$$
\ell(u,v)=
\begin{cases}
\min\{\alpha u+(\gamma-\alpha)v,
       (\gamma-\beta)u+\beta v+q\beta\},&u\ge v,\\
\min\{\alpha u+(\gamma-\alpha)v+q\alpha,
       (\gamma-\beta)u+\beta v\},&u<v.
\end{cases}                                               \tag{D5}
$$

It is affine on the pieces cut out by three lines

$$
u-v=0,\qquad u-v=q\beta/\delta,\qquad
u-v=-q\alpha/\delta.                                    \tag{D6}
$$

**Theorem 2.** Suppose also $\gamma\ne\alpha,\beta$. Every maximum
spanning tree uses only edges $(u,v)$ on the boundary of the index square
or satisfying $|u-v-t|\le2$ for at least one of the three levels $t$
in (D6). This is a set of at most $19q$ edges, enumerable in $O(q)$
arithmetic operations.

**Proof.** At any other index pair, choose the unit step in each
coordinate in the direction of the affine piece's positive gradient.
Both gradient components are nonzero. The resulting $2\times2$
index rectangle stays within the square and within the same piece:
its differences change by at most two. All three other edge costs in
this four-cycle are strictly larger. The original edge is its unique
minimum and therefore belongs to no maximum spanning tree, by the
ordinary cycle-exchange argument.

There are at most $4q$ boundary pairs. Each width-four strip around a
diagonal level contains at most five integer diagonals, each of length
at most $q$. The three strips add at most $15q$ pairs. Enumeration
uses those diagonal offsets directly; it need not inspect $q^2$ pairs.
$\square$

This proof is an instance of a general local cycle-domination principle:
for a piecewise-affine grid cost, an index rectangle contained in a
nonzero-gradient piece excludes its unique minimum from a maximum tree.
The cycle property is classical. The present bound uses the three
specified lines of the periodic cost, not arbitrary parametric edge data.

## Compressed standard sets

Theorem 2 and a classical spanning-tree algorithm compute the omitted
tree from $O(q)$ edges in $O(q\log q)$ arithmetic operations/comparisons
and $O(q)$ storage. Least representatives are evaluated on demand by
the earlier constant-candidate formula. Store this rule and the omitted
tree, instead of constructing the complete $q^2$-entry array.

This is an $O(\sqrt n)$-size representation of an initial standard set
for $n=(q-1)^2$. A monomial membership query computes its periodic
label, checks whether it is its label's least representative, and tests
whether the label belongs to the omitted tree. A sorted tree-label list
gives $O(\log q)$ comparisons per query. Listing the complete standard
set still requires $\Theta(n)$ output and is a separate task. These
are arithmetic-model bounds, not bit-complexity or dense reduced-basis
output bounds. On the excluded parameter planes, take a specified order
refinement; generic weights suffice for studying all open initial cones.

## A remaining obstruction, even with graph compression

For the bilinear outer-cone theorem, the transformed exponent nodes of
each cost prefix are **consecutive integers**. At nonnegative $\rho$,
the horizontal and vertical slots give an interval $-I,\ldots,J$.
At negative $\rho$, the first horizontal slots give an interval ending
at $-t+L-1$; after all horizontal slots appear, vertical slots extend
the same interval through zero. Their first $L$ binomial Taylor rows
have determinant $1$, after ordering the nodes, since

$$
\det[\tbinom{m+j}{i}]_{0\le i,j<L}=1
\quad(m\in\mathbb Z).                                   \tag{D7}
$$

Unlike an arbitrary separated-node Vandermonde argument, this integer
identity survives every characteristic, including when $L>p$. Thus
the bilinear case of the [outer-cone theorem](rectangle_boundary_fan.md)
and its exact regional count hold over **all** fields. This extension
does not apply to its general $d>1$ intervals.

Take $q=2^k\ge4$ in characteristic two, side $r=s=q-1$, and
$a=q/2$. The primitive coherent direction

$$
(a,-(a-1),-a),\qquad
\prod(|h_i|+1)=a(a+1)^2=\Theta(n^{3/2})                  \tag{D8}
$$

therefore occurs in this very cographic family. Its graph/border
coefficients are $\pm1$, so the same cone also occurs for the
characteristic-independent $J_{q,q}$. A small compressed basis does
not imply a linear primitive-direction product. The number of fan cells
remains a separate question.

## Scope and next counting task

Complete finite fans at $q=3,4,5,8$ have respectively $6,14,26,86$
states. Exact positive cone coverage and compressed witnesses are saved.
These motivate the concrete conjecture
$g(J_{q,q})=2(q-1)(q-2)+2$, which would be only linear in the quotient
length. Neither the branching theorem nor the sparse graph proves it.
No larger numerical sweep is justified before attempting that missing
classification. A proof would exclude this family, not complete the
unrestricted research goal.

### A classical corank screen for a higher-dimensional extension

There is a useful general check before moving to more complicated tensor
matroids. Suppose every initial standard set of a length-$n$ quotient
in three variables lies in one integer box with $M=L_1L_2L_3$ monomials.
Let $e=M-n$. Then the established lattice-polytope vertex bound gives

$$
g(I)=O(1+e^{3/2}\sqrt M).                                \tag{D9}
$$

Indeed use the monomial column matroid on that box. A basis complement
is a rank-$e$ basis in its dual; their profiles differ by reflection and
translation by the sum of the whole box. Every dual profile has coordinate
at most $e(L_i-1)$. If necessary add a unit cube to make the profile
polytope full-dimensional; its normal fan refines the old one, so vertex
count cannot decrease. Its volume is $O(e^3M)$ for $e\ge1$. The same
classical three-dimensional lattice-vertex theorem used in the earlier
state bounds now gives (D9). For $e=0$ there is just one standard set.
This is a direct corollary of matroid duality and the classical volume
bound, not a new lattice theorem. In particular, a common box of size
$O(n)$ already excludes a polynomial separation above $n^2$.

This rules out one tempting nonregular extension without a fan search.
For odd $q\ge3$, take the three-factor algebra
$K[U,V,W]/(S_q(U),S_q(V),S_q(W))$ and the three generators
$X=UV,Y=UW,Z=VW$. Their period map on $(\mathbb Z/q)^3$ has determinant
$-2$ and is invertible. The generator algebra is therefore the full
quotient, of dimension $n=(q-1)^3$. An explicit ideal presentation adds
$X^q-1,Y^q-1,Z^q-1$ and the three relations
$S_q(XY/Z),S_q(XZ/Y),S_q(YZ/X)$; inverses mean powers $q-1$.
Under the period isomorphism these are $S_q(U^2),S_q(V^2),S_q(W^2)$,
which equal the corresponding geometric sums modulo the periods.
This proves the realization over every field.

Every standard set is inside the $q^3$-point coordinate box. Here
$e=q^3-(q-1)^3=O(q^2)$, so (D9) gives

$$
g(I)=O(q^{9/2})=O(n^{3/2}).                               \tag{D10}
$$

The three-circuit tensor matroid is nonregular, as covered by the
[triple-factor obstruction](https://arxiv.org/html/2609.29227v1#S4.SS3).
That extra dependency complexity does not overcome the small-box screen.
No larger fan experiment is justified for this equal-period construction
as a route to the requested separation.

Distinct odd prime periods give a different, defined case: use
$X=\zeta_{p_1}\zeta_{p_2}$,
$Y=\zeta_{p_1}\zeta_{p_3}$,
$Z=\zeta_{p_2}\zeta_{p_3}$ over $\mathbb Q$. These generate
$\mathbb Q(\zeta_{p_1p_2p_3})$; modulo the odd product, squaring is
invertible, so each individual root is recoverable from the pair products.
The monomial labels are $(a+b,a+c,b+c)$ reduced in the three respective
prime coordinates. The known tensor/simplicial matroid specifies their
ranks, but their least representatives vary with the three prices.
The elementary period box now has size $(p_1p_2p_3)^2$, rather than
$O(n)$. Thus (D10) does not close this case. This is a next candidate
requiring a uniform coherence count, not a proposed result or a claim
that no better exclusion applies.

The proof mechanisms now restrict the omitted trees to ordered paths
with two possible branching vertices and a sparse set of affine pieces.
The next requirement is a uniform count of the **coherent** possibilities
as the weights vary, or a proved exclusion from polynomial separation.
The candidate set itself varies with the weights. A sparse graph at one
weight does not give a bound on the number of states over all weights.

No exact all-$q$ fan formula, unrestricted $G_3$ improvement, strongest
current bound, or novelty certificate is asserted here. Prior manuscripts
and the original success criteria remain unchanged.

```bash
python3 -B scripts/run_cographic_structure_checks.py
```

The checks use independent polynomial reduction, exhaustive small downset
bases, dense versus sparse graph optimization, local four-cycle witnesses,
small-characteristic binomial-prefix ranks and actual coherent facets.
The bounded phase references also walk complete small fans independently
using original local multiplication. The corank controls verify an actual
three-generator algebra and the complete polynomial relation space.
The report states the finite scope separately for each part.
