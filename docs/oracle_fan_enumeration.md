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
asserted. Separately, one can use

\[
 O((g+1)^2\log(g+2)+(g+1)nL^2)                            \tag{O2}
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

field operations. The separate rational geometric-operation budget
from (O2) is $O(n^{23/4}L^{7/4})$ and can dominate when fast matrix
multiplication is used. This is an arithmetic bound, not a bit-complexity, benchmark or
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

Equations (O6) and (O10) now prove (O1). For (O2), rebuild a
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

The geometry is not free: its bound is $O(n^6L)$ with the classical
count and $O(n^{23/4}L^{7/4})$ with the shadow count. Combining all
unit-cost operation types therefore gives
$O(n^{\omega+23/8}L^{11/8}+n^{23/4}L^{7/4})$, rather than the
field-operation term alone. No coefficient-bit bound follows.

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
