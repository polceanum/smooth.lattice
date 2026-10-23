# Generic matching is exact for one-step pure arm prefixes in characteristic zero

The first-step transverse operator \(W=1+vH(q)\) admits an exact generic
rank rule for every pure cost prefix in characteristic zero. An explicit
single coefficient choice attains all of these ranks simultaneously.
The proof uses a planar path matrix only after checking the boundary order
of its sources and sinks. It does not replace the coupled two-arm head
rank problem by matching.

Let \(A\ge t\ge2\), and work in
\(R_A=K[q,v]/(q^A,v^t)\). The pure columns and prices are
\[
 U_{i,j}=q^iv^j,\quad \theta+i+\rho j
       \quad(0\le i<A,\ 0\le j<t),
\]
\[
 L_{i,j}=Wq^iv^j,\quad i+\rho j
       \quad(t\le i<A,\ 0\le j<t),                    \tag{P1}
\]
where \(\rho>0\) and \(\theta\in\mathbb R\). A pure prefix retains
every column whose price is below a cutoff avoiding all prices. Its
support rank is the maximum matching size of its row-column support graph.

**Theorem.** In characteristic zero, for
\[
                           H(q)=(1+q)^D,\qquad D\ge0,          \tag{P2}
\]
every pure prefix in (P1) has actual rank equal to its support rank,
for all \(\theta,\rho\) and strict cutoffs. Consequently, for fixed
\(A,t,D\), a nonempty coefficient-open subset of degree-\(D\) polynomials
has this property simultaneously for all pure prefixes. Only coefficients
below \(q^{A-t}\) affect the lower columns.

In particular, over this coefficient-open set every strict pure minimum
basis is exactly its support-greedy basis. The entire pure basis fan can
therefore be computed by matching without coefficient-field elimination.
This is a pure arm statement: it does not identify the coordinates of the
heads in those bases or prove a whole hook-fan count.

The same proof works over the reals for
\(H=\prod_{k=1}^D(1+\beta_kq)\), \(\beta_k>0\). The generic statement
over arbitrary characteristic-zero fields follows from the integer point
(P2), rather than from an order on the field.

## Remove the short long-coordinate rows

Every lower vector is divisible by \(q^t\). Selected upper columns with
\(i<t\) are independent coloops and may be counted separately. On the
remaining rectangle put \(B=A-t\) and replace \(i\) by \(t+i\).
The case \(B=0\) is immediate. It suffices to study
\[
 V=\{(i,j):0\le i<B,\ 0\le j<t\},\qquad c(i,j)=i+\rho j,
\]
and the columns \(I_U\) and \((I+T)_S\), where
\[
 U=\{a\in V:c(a)<\ell\},\qquad
 S=\{a\in V:c(a)<h\},\qquad
 T(q^iv^j)=vH(q)q^iv^j.                               \tag{P3}
\]
The cutoffs are \(h=\text{cutoff}-t\) and \(\ell=h-\theta\).

Eliminating the identity columns gives
\[
 \operatorname{rank}[I_U\mid(I+T)_S]
       =|U|+\operatorname{rank}(I+T)_{V\setminus U,S}.          \tag{P4}
\]
The analogous equality holds for support rank. Indeed a maximum matching
can be made to use every identity column: if its row is occupied by a
lower column, exchange that matched column for the identity column.

If \(\ell\ge h\), then \(S\subseteq U\), so the remaining matrix has
only the \(T\) entries. It splits into independent one-step Toeplitz
blocks. Their matching ranks are attained by (P2), as the elementary path
construction below shows. Thus assume \(\ell<h\).

## Schur elimination gives a boundary transfer matrix

Write \(C=S\setminus U\) and \(O=V\setminus S\). Multiplying the
row and column at level \(j\) by \((-1)^j\) replaces \(I+T\) by
\(I-T\) without changing ranks or support. The diagonal block on \(C\)
is \(I-T_{CC}\), which is invertible because every \(T\)-edge increases
\(j\). Its Schur complement gives
\[
 \operatorname{rank}(I-T)_{V\setminus U,S}=|C|+
                                      \operatorname{rank}P,   \tag{P5}
\]
\[
 P=T_{OU}+T_{OC}(I-T_{CC})^{-1}T_{CU}.                         \tag{P6}
\]
Thus \(P_{o,u}\) sums the weights of directed paths from \(u\in U\)
to \(o\in O\), with every intermediate coarse vertex in \(C\).
Original vertices in \(U\) cannot be internal vertices; original vertices
in \(O\) are terminal. This restriction is essential.

Let \(\kappa\) be the maximum number of vertex-disjoint such paths in
the coarse support graph. Then
\[
 \operatorname{termrank}(I-T)_{V\setminus U,S}=|C|+\kappa.      \tag{P7}
\]
One proof puts distinct indeterminates on all coarse off-diagonal edges.
Every nonzero matching minor has a distinct edge monomial for each
permutation, since the remaining diagonal entries are fixed ones. Hence
its generic rank is its term rank. Apply (P5) over that function field.
The Lindström--Gessel--Viennot cancellation leaves vertex-disjoint path
systems. Distinct disjoint systems have distinct edge monomials in this
acyclic graph, so \(\operatorname{rank}P=\kappa\) there. This proves
(P7) without assuming any positivity for the original repeated coefficients.

## Refine each Toeplitz step into a planar strip

For (P2), replace the edge from level \(j\) to \(j+1\) by \(D\)
small steps. Vertices have coordinates
\[
                  (s,i),\qquad s=j+k/D,
\]
and each small step increases \(s\) by \(1/D\) and increases \(i\)
by zero or one. All edge weights are one. The number of paths from
\((j,i)\) to \((j+1,i+r)\) is \(\binom Dr\), as required.
For \(D=0\), use a single horizontal step with no upward option.
The positive factor weights \(\beta_k\) give the more general real
version in the same way.

Refinement preserves the maximum number \(\kappa\). Any fine disjoint
system gives a coarse disjoint system. Conversely, at one coarse step,
list its input and output levels increasingly as
\(a_1<\cdots<a_m\) and \(b_1<\cdots<b_m\). A feasible interval
matching can be uncrossed to the increasing matching, so
\(0\le b_i-a_i\le D\). Explicit disjoint routes are
\[
                 q_i(k)=\max\{a_i,b_i-D+k\},\quad 0\le k\le D.
                                                                  \tag{P8}
\]
Both sequences inside the maximum increase strictly with \(i\), so
these routes never share a vertex. Reconnect them at the next coarse
level. The endpoint sets are unchanged; since coarse time increases,
this preserves the number of full input-output paths. This also proves
the Toeplitz-block assertion needed when \(\ell\ge h\).

## The cost cuts force compatible path orders

Order the input set \(U\), and separately the output set \(O\), by
increasing coarse time \(j\), breaking ties by decreasing \(i\).
We claim that every fine vertex-disjoint path system pairs its selected
inputs with its selected outputs in these orders.

Extend cost to the drawing by \(c(s,i)=i+\rho s\). It increases
strictly along every small edge. Move \(\ell,h\) slightly, if necessary,
without changing their coarse sets, so neither cut passes through a fine
vertex. Every path crosses each of the two parallel lines
\(c=\ell\), \(c=h\) exactly once.

The first crossing lies in the first coarse step after its input: otherwise
the path would pass through another vertex of \(U\), which is forbidden.
For inputs at different coarse times, these first steps have disjoint
time interiors in the same order. For two inputs in the same step, the
higher path stays strictly above the lower one until that step ends;
therefore it crosses \(c=\ell\) first. Thus input order equals crossing
order along the first line, in increasing \(s\).

The last crossing similarly lies in the last coarse step before the
output, because an \(O\)-vertex cannot be internal. The same argument
shows that output order equals crossing order along \(c=h\).
Between the cuts the paths lie in the convex strip inside the drawing.
Disjoint planar paths cannot reverse their boundary order across that
strip. This proves the claim.

Consequently the ordinary path determinant has no negative surviving
terms. Every minor of \(P\), with the stated row and column orders,
is nonnegative, and it is positive exactly when its endpoints admit a
vertex-disjoint system. All weights in (P2) are one. A system of size
\(\kappa\) therefore gives a positive minor of order \(\kappa\):
\[
                         \operatorname{rank}_{\mathbb Q}P=\kappa.
                                                                  \tag{P9}
\]
Equations (P4)--(P9) prove equality with support rank for every prefix.
Reduction modulo a prime can destroy this positive integer minor; no
positive-characteristic conclusion is being imported.

## Generic coefficients and an original-prefix consequence

For fixed \(A,t,D\), only finitely many distinct pure prefixes can
occur. For each prefix, a maximal matching minor is a polynomial in the
coefficients of \(H\), nonzero at the same integer point (P2). Their
product, together with the nonzero coefficient conditions fixing the
degree and support, defines a nonempty open set over every
characteristic-zero field. All those ranks hold simultaneously there.
Equality of every prefix rank gives equality of the corresponding strict
greedy basis decisions, proving the pure-fan consequence.

There is also a direct consequence for the original central-only hook
with arm operators \(1+vH(q)\) and \(1+vG(q)\). Consider a complete
original cost prefix containing no upper or pure lower Z-arm columns.
Keep every head that is actually cheaper than its cutoff. Those head
columns have independent Z projections: they are distinct corner
monomials multiplied by the unit \(1+vG(q)\). Thus their actual rank
contribution is exactly their number, independently of the pure Y span.
Their constant-term Z entries give the same support matching contribution.
The theorem therefore identifies the rank of this complete original
prefix, with its heads retained, for arbitrary \(G\). The transposed
statement holds as well.

When both arms contribute other columns, those projections can be
contracted or coupled. The pure theorem supplies no replacement for
that remaining simultaneous head-rank test. In particular it proves
neither matching for all original prefixes nor a full-fan bound.

## Verification and priority boundary

Run from the repository root:

```bash
python3 -B scripts/run_charzero_pure_prefix_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checks compare actual entire pure prefix ranks with independent
matchings, compare coarse Schur transfers with direct ranks and fine
vertex-capacitated flows, inspect transfer minors, and verify the ordered
routes (P8). Complete small pure fans and original one-sided prefixes
provide separate controls. Finite checks support the proof; they do not
prove its quantified statement by exhaustion at growing size.

The path determinant formula is classical; see
[Abdallah--McDaniel, Fact 4.2](https://alco.centre-mersenne.org/item/10.5802/alco.414.pdf#page=14).
Their discussion after Lemma 4.4 and Proposition 4.5 also illustrates why
planarity without compatible terminal order permits cancellation.
[Postnikov, Sections 4--5](https://arxiv.org/pdf/math/0609764)
develops the classical boundary-measurement and positivity framework.
The cost-strip boundary argument above verifies the needed ordering for
this particular repeated Toeplitz system. Exact priority of this
specialization has not been established.

The [Frobenius complete-prefix obstruction](frobenius_coherent_prefix_gaps.md)
continues to exclude an all-field generic rule. Arbitrary selected minors
also need not satisfy matching in characteristic zero. The unrestricted
sharp \(G_3\) endpoint, polynomial actual fan excess, consequential
general theorem, current novelty and final-paper requirements remain open.
