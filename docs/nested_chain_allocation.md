# Chain allocation in the explicit nested quotients

This continues [the multiplicative construction](nested_jet_matroids.md).
For its curvilinear family, all basis choices reduce exactly to integer
convex allocation with ascending constraints. This gives a certificate for
each positive initial cone and a monotonicity law for changing weights.
It also bounds the entire equal-valuation family uniformly and counts one
of its strata sharply. The uniform bound excludes a polynomial separation
above the generic scale for that family. These statements do not settle
$G_3(n)$ or establish novelty.

The reduction concerns the fixed characteristic-zero ideal $I_{n,A,B}$
already constructed, with valuation

$$
v(a,b,c)=a+Ab+Bc,\qquad A,B\ge1.
$$

Its length is exactly $n$. Its $n$ distinct monomial residues are a basis
precisely when their increasing valuations satisfy $v_i\le i-1$.
The ambient truncation equals the length here; this specialization is
essential. The results below do not assert the same reduction for
$m_1>1$ or arbitrary polynomial quotients.

## An exact linear description of chain bases

Index the $x$-chains by pairs $j=(b_j,c_j)$ with
$q_j=Ab_j+Bc_j<n$, ordered by nondecreasing $q_j$; break ties arbitrarily.
For nonnegative integers $\ell_j$ with $\sum_j\ell_j=n$, put

$$
D(\ell)=\{(a,b_j,c_j):0\le a<\ell_j\}.
$$

**Proposition 1.** The residues of $D(\ell)$ are a basis if and only if

$$
\sum_{q_i<q}\ell_i\ge q
\quad\text{for every distinct chain start }q<n.                 \tag{1}
$$

Equivalently, $\sum_{i<j}\ell_i\ge q_j$ for every ordered chain $j$.
Constraints after the first chain with a repeated start are redundant.

For necessity, the valuation condition for a full basis says that at least
$q$ selected monomials have valuation less than $q$. All of them come
from chains starting below $q$, giving (1). For sufficiency, concatenate
the chosen chains in their displayed order. If $S$ monomials precede
chain $j$, (1) gives $S\ge q_j$. Its $a$th monomial then has valuation
$q_j+a\le S+a$, its index starting at zero. An assignment of distinct
indices $0,\ldots,n-1$ dominating the valuations implies the same
inequalities after sorting. This proves the basis criterion. It also
implies $\ell_j\le n-q_j$, so no selected monomial is a zero residue.

This does not say that every feasible allocation is a monomial downset.
It describes bases closed under division by $x$. That larger feasible
class is useful and sufficient for optimization.

## Exact weight optimization and cone certificates

Let $w=(w_x,w_y,w_z)>0$ and $s_j=b_jw_y+c_jw_z$. Define

$$
F_w(\ell)=\sum_j\left[w_x\binom{\ell_j}{2}+s_j\ell_j\right].    \tag{2}
$$

The minimum-weight residue basis is exactly a minimizer of (2) under
(1), nonnegativity and total $n$. Indeed, in any basis compressing the
selected exponents in each $x$-chain to $0,\ldots,\ell_j-1$ lowers
valuations and preserves independence. It strictly lowers weight unless
the chain was already a prefix, since $w_x>0$. Every minimum basis is
therefore of this form. A unique minimum basis is the standard monomial
set of its positive monomial initial ideal.

Here is a complete certificate using only the chain lengths. A transfer
$\ell-e_i+e_j$, $i\ne j$, is feasible precisely when $\ell_i>0$ and,
if $q_i<q_j$, every quota crossed by the transfer has integer slack at
least one. Moving to an earlier or equal start is always feasible.
Write

$$
r_i=s_i+w_x(\ell_i-1),\qquad a_j=s_j+w_x\ell_j.                \tag{3}
$$

**Proposition 2.** An allocation is optimal if and only if
$r_i\le a_j$ for all feasible transfers. It is uniquely optimal if and
only if all these inequalities are strict. Thus its open positive
initial cone is given exactly by

$$
w\cdot(\ell_j-\ell_i+1,\ b_j-b_i,\ c_j-c_i)>0                 \tag{4}
$$

for those pairs, together with positivity.

Necessity follows by taking one transfer. For sufficiency, cut the
ordered chains at the quotas tight at $\ell$. Relative to any other
feasible allocation $m$, the cumulative number of missing units at each
such cut is at least the cumulative surplus. Consequently surplus units
can be sent to missing units without moving right across a tight cut:
process recipient blocks from left to right, using donors in that or
later blocks. Every resulting donor-recipient pair is a feasible
transfer at $\ell$. Discrete convexity bounds the cost of removing
multiple units from a chain below by minus its number removed times
$r_i$, and the cost of adding units below by their number times $a_j$.
Summing these inequalities along the transport gives
$F_w(m)-F_w(\ell)\ge0$. Strict transfer inequalities give a positive
difference whenever $m\ne\ell$. Conversely, equality on a feasible
transfer gives another optimum. This proves both assertions.

For example, checking (4) avoids computing a separate residue determinant
for every proposed neighboring basis. It is an exact certificate, not a
claim of a new optimization algorithm or a runtime advantage over the
existing scheduling oracle.

## A monotonicity law for releasing rank constraints

The following is a classical submodular uncrossing argument applied to
this quotient. Let the distinct starts be
$0=q_0<q_1<\cdots<q_{s-1}<n$. Fix $w_x>0$ and base offsets $s_j^0$,
and vary the chain costs by

$$
s_j(\tau)=s_j^0+\tau q_j.                                    \tag{5}
$$

For the actual three variable weights this is the line
$w_y(\tau)=w_y^0+A\tau$, $w_z(\tau)=w_z^0+B\tau$, with $w_x$
fixed; restrict to its positive portion. Let $T_d$ be the total length
of chains starting at $q_d$, and $P_k=\sum_{d<k}T_d$.

**Proposition 3.** If $\tau<\tau'$, every optimal cumulative vector at
$\tau$ is componentwise at most every optimal cumulative vector at
$\tau'$. In particular, a quota $P_k\ge q_k$ can cease to be tight
only once along this line. The number of distinct uniquely optimal
initial sets on the line is at most

$$
1+\sum_{k=1}^{s-1}(n-q_k).                                  \tag{6}
$$

To prove this, minimize within each group for a prescribed total $T$.
The resulting cost $f_d(T)$ is discretely convex: its increments are
the merged, sorted marginal costs $s_j^0+w_xa$, $a\ge0$, of that
group's chains. The feasible cumulative vectors form a finite sublattice:
$P_0=0$, $P_s=n$, $P_k\ge q_k$, and $P_k\le P_{k+1}$.
Apart from a constant, the cost is

$$
H(P)-\tau\sum_{k=1}^{s-1}(q_k-q_{k-1})P_k,\qquad
H(P)=\sum_{d=0}^{s-1}f_d(P_{d+1}-P_d).                       \tag{7}
$$

The function $H$ is submodular on this lattice. For each adjacent
coordinate pair, crossing two vectors replaces the two increments by
two more balanced increments with the same sum; convexity gives the
inequality. Summing gives
$H(P)+H(Q)\ge H(P\wedge Q)+H(P\vee Q)$.

If $P,Q$ are optima at $\tau,\tau'$, compare them with
$P\wedge Q,P\vee Q$ at the respective parameters. Optimality makes
the summed difference nonpositive, while submodularity bounds it below
by

$$
(\tau'-\tau)\sum_{k=1}^{s-1}(q_k-q_{k-1})(P_k-Q_k)_+.
$$

All coefficients are positive, so $P\le Q$. Within each fixed group
total, the objective deciding chain lengths is independent of $\tau$.
Two different uniquely optimal initial sets must therefore have different
cumulative vectors. Their sum of cumulative coordinates increases by at
least one each time and has range at most the sum in (6).

For $n=M^2$, $A=B=M$, this specializes to

$$
P_k\ge Mk\quad(1\le k<M),\qquad
\#\text{strict states on the line}\le1+\frac{M^2(M-1)}2.       \tag{8}
$$

The line estimate is not asserted to improve the existing general
matroid shadow bound. The useful additional information is the ordered,
irreversible release of the rank quotas. It does not bound the number
of interactions as the direction of the line changes.

## A uniform bound for the entire equal-valuation family

**Proposition 4.** For every $n,M\ge1$, the fixed ideal $I_{n,M,M}$
has at most

$$
O\bigl(n^2\log^4(n+1)\bigr)                                  \tag{U1}
$$

positive monomial initial ideals, with a constant independent of $M$.
In particular the large-box family $n=M^2$ cannot give a polynomial
separation above the generic $n^2\log^2(n+1)$ scale. The proof has
three steps: an active exchange is an ordinary cut on a consecutive
degree block; each such cut has a small box in one fixed lattice
coordinate system; all possible facet directions then fit a hyperbolic
cross. The last arrangement argument is classical.

### Every active exchange has a consecutive-block cut certificate

Partition the degree groups at the quotas tight at a feasible allocation.
Include the virtual endpoints zero and total $n$. Each resulting block
has positive total allocation. At an optimal weight define, for a block
$E$, $R_E$ as the largest removal cost $r_i$ among its occupied chains,
and $A_E$ as the smallest addition cost $a_j$ among all its chains.
All unit transfers within a block are feasible, so Proposition 2 gives
$R_E\le A_E$. The same-chain inequality follows from $w_x>0$.
For an earlier block $E$ and later block $F$, all transfers from $F$
to $E$ are feasible, so $R_F\le A_E$.

Suppose a feasible transfer is active, with removal and addition cost
both $t$. If its two chains are in the same block, that block satisfies
$R_E=A_E=t$. Otherwise the donor is in a later block $F$ and the
recipient in an earlier block $E$: a rightward transfer crossing a tight
quota is infeasible. Necessarily $R_F=A_E=t$. Every block $G$ between
them satisfies

$$
R_G\le t\le A_G,
$$

including both endpoints. Thus the union of these consecutive blocks
consists of $N\le n$ cheapest points, with ties allowed, in all the
$x$-chains whose degrees belong to an interval $i\le b+c\le j$.
Every selected point has cost at most $t$, and every unselected point
has cost at least $t$. The removed and added points are distinct
boundary points. Infinite continuation of each $x$-chain is harmless:
its first unselected point has cost $a_j$, and all later costs are larger.

### A cut in the shifted wedge has a small relative box

Work first in $w_y>w_z$. Use the fixed unimodular coordinates

$$
(a,b,c)\longmapsto(a,b,d=b+c),\qquad
\alpha=w_x>0,\quad\beta=w_y-w_z>0,\quad\gamma=w_z>0.
$$

For the degree interval just obtained, translate $e=d-i$. Its entire
ground becomes

$$
\mathcal W_{i,j}=\{(a,b,e)\in\mathbb N_0^3:
                   e\le j-i,\ b\le i+e\}.
$$

The closed cut $C$ has score $\alpha a+\beta b+\gamma e\le T$,
where $T=t-\gamma i>0$. It has at most $N-1$ points strictly below
the cutoff. Coordinatewise floor-halving preserves $\mathcal W_{i,j}$:
$\lfloor b/2\rfloor\le i+\lfloor e/2\rfloor$. It sends every
point of $C$ strictly below $T$, with fibers of size at most eight.
Therefore $|C|\le8(N-1)$.

Let $A,B,D$ be the three actual coordinate counts of its bounding box.
The wedge inequality gives $B-1\le i+D-1$. The following grid lies
inside $C$:

$$
\begin{split}
0\le a&\le\lfloor(A-1)/4\rfloor,\\
0\le b&\le\lfloor(B-1)/4\rfloor,\\
\lfloor(D-1)/4\rfloor\le e&\le\lfloor(D-1)/2\rfloor.
\end{split}
$$

Its $b$ bound is at most $i$ plus its minimum $e$ bound. Each
coordinate maximum defining $A,B,D$ has its individual cost at most
$T$, so every displayed grid point has cost at most
$(1/4+1/4+1/2)T=T$. There are at least $ABD/64$ grid points.
Consequently

$$
ABD\le64|C|\le512(N-1)\le512n.                              \tag{U2}
$$

This is a relative box for an exchange block in the fixed coordinates
$(a,b,d)$, translated in $d$. It is not a small box for the whole
standard set. The earlier $n^{3/2}$ footprint remains valid.

### The direction count gives the uniform fan bound

The difference $h$ of the two tied points in these coordinates satisfies
$\prod_{k=1}^3(|h_k|+1)\le ABD\le512n$. Dividing by its integer
gcd only decreases this product. Every interior facet of a full initial
cone comes from an active feasible-transfer inequality in (4), so its
primitive normal in this one coordinate system belongs to

$$
\{h\in\mathbb Z^3:\prod_k(|h_k|+1)\le512n\}.
$$

This set has $O(n\log^2(n+1))$ elements, by summing
$\lfloor512n/(uv)\rfloor$ over positive $u,v$. The central planes
with these normals have $O(n^2\log^4(n+1))$ regions. They refine
the initial fan on $\alpha,\beta,\gamma>0$, proving (U1) on that
half of the positive weight space. Swap $y,z$ for the other half.
A uniquely optimal set on $w_y=w_z$ persists under a small perturbation,
so the dividing plane adds no exceptional initial ideals.

This proof closes the equal-valuation family as a source of the desired
power separation. It does not give the sharp logarithmic factor for that
family or extend to arbitrary valuation pairs $A,B$. For unequal pairs,
consecutive starts cut an arithmetic strip $Ab+Bc$ rather than the
degree wedge above; (U2) in one common lattice system has not been proved.

## Counting the stratum where all quotas are tight

Continue with $n=M^2$, $A=B=M$, $M\ge2$. If every quota is tight,
then each degree group $d=b+c$, $0\le d<M$, has $T_d=M$.
For $w_y>w_z$, put

$$
\delta=(w_y-w_z)/w_x>0,\qquad p=w_z/w_x>0.
$$

The optimal allocation inside group $d$ consists of the $M$ cheapest
points of the strip $\{(a,b):a\ge0,\ 0\le b\le d\}$ with score
$a+\delta b$. The common cost $pdM$ does not affect this choice.
The case $w_z>w_y$ swaps $y,z$.

**Proposition 5.** The number of positive monomial initial sets for
which all the quotas are tight is

$$
\Theta(M\log(M+1))=\Theta(\sqrt n\log(n+1)).                 \tag{9}
$$

This is a consequence of the classical planar prefix count, not an
unrestricted initial-ideal counting theorem.

First, every joint tuple of strict strip prefixes is realizable. At a
generic $\delta$, let $r_d$ and $a_d$ be the $M$th and $(M+1)$st
strip scores. We have $r_d<a_d$ and $r_j\le r_i$ for $i<j$, since
the strips increase. Choose

$$
0<p<\min_{i<j}\frac{a_i-r_j}{j-i}.                           \tag{10}
$$

All transfers within a group increase cost. All transfers from a later
group $j$ to an earlier group $i$ increase cost by at least
$a_i-r_j-p(j-i)>0$; the reverse transfers are blocked by tight quotas.
Proposition 2 proves unique optimality. Thus the tuple occurs for one
fixed ideal, rather than for separately chosen coefficient matrices.

For the upper bound, suppose a strip prefix changes at
$\delta=p_0/q_0$ in lowest terms. Two tied boundary points differ by
$(p_0,-q_0)$, so the cutoff score is at least $p_0$, and the strip
contains all points below the triangle with vertices
$(0,0),(p_0,0),(0,q_0)$. The number strictly below its hypotenuse is

$$
\frac{p_0q_0+p_0+q_0-1}{2}.
$$

Reflection in the center of the containing rectangle pairs points below
and above that line; its only lattice points are the two endpoints
because $p_0,q_0$ are coprime. A rank-$M$ change has at most $M-1$
points strictly below its cutoff. Consequently every wall in every group
satisfies the same bound

$$
(p_0+1)(q_0+1)\le2M.                                       \tag{11}
$$

There are $O(M\log(M+1))$ such positive integer pairs. Their union
partitions the positive $\delta$ line into that many intervals, with
at most one joint tuple per interval. Swapping $y,z$ costs a factor at
most two. A uniquely optimal set at $w_y=w_z$ persists under a small
perturbation and adds no exception.

For the lower bound, the strip $d=M-1$ has exactly the ordinary planar
$M$-point prefix possibilities: an $M$-point positive prefix never
needs $b\ge M$. The established planar count is
$\Theta(M\log(M+1))$. Different last-group prefixes give different
joint tuples, and (10) realizes each. This proves (9).

Large bounding boxes occur even in this small stratum. Choose generic
$0<\delta<1/(M-1)$, then take $p$ as in (10). Every group includes
its entire row $a=0$, since those $d+1$ scores precede the score of
$(1,0)$ and $d+1\le M$. Group zero contains $x^0,\ldots,x^{M-1}$.
Thus all three axis counts are $M$, and the box has $M^3=n^{3/2}$
points. Neither that footprint nor the all-tight stratum supplies a
polynomial separation above the generic fan scale.

## An ordinary-cut sector and what is left

If $w_y\ge Aw_x$ and $w_z\ge Bw_x$, the ordinary $n$ cheapest
orthant monomials already form a residue basis. A selected monomial of
valuation at least $k$ has cost at least $kw_x$, so the prefix contains
$1,x,\ldots,x^{k-1}$. If there is no such monomial, all $n$ valuations
are less than $k$. Either way the full-basis valuation inequalities hold.
Ties can be resolved by a small perturbation. The quotient's positive
fan in this sector is exactly the ordinary corner-cut fan restricted
to the sector.

For equal valuations, Proposition 4 bounds those interactions uniformly
and excludes the whole family as a power-separation route, including
the partial-release strata. For unequal valuations, the arithmetic
strips between consecutive tight starts remain the concrete obstacle.
A uniform direction bound for their active exchanges, or a proved
large family of different interactions, is the next construction task.
The general $m_1>1$ realization is also preserved. Neither (U1) nor the
classical allocation reductions settle the original extremal goal, and
no new paper follows from these reductions alone.

## Literature boundary and verification

The optimization class is established. The
[Akhil--Sundaresan survey, Sections 2 and 4](https://arxiv.org/html/1608.08000v1)
places ascending constraints in polymatroid optimization and discusses
the integer problem. Those sections, including their exchange criterion
and integer greedy formulation, were inspected. The continuous
[Padakandla--Sundaresan paper](https://arxiv.org/abs/0707.2265)
was inspected only at abstract level for this continuation; it is not
used as an integer optimality theorem.

[Vidal--Jaillet--Maculan, equations (1.1)--(1.4), Theorems 3.1--3.2 and the integer proof](https://www.mit.edu/~jaillet/general/RAP-LB-15.pdf)
study convex allocation with nested constraints for integer and continuous
variables. Reversing our chain order converts cumulative lower bounds to
their upper bounds; a redundant bound $\ell_j\le n$ is sufficient.
These optimization statements are prior art. Their inspected statements
do not by themselves count all initial ideals of the fixed quotient.

The monotonicity argument is an application of classical lattice
comparative statics; the [Topkis 1978 publisher abstract](https://doi.org/10.1287/opre.26.2.305)
describes that general subject. Only that abstract was accessed here;
the complete special-case proof is given above. The planar counting input
and its earlier primary comparisons remain recorded in
[the slab note](nested_jet_matroids.md#uniform-slab-bounds-and-a-sharp-fixed-height-count).
No novelty is inferred from a bounded search or from the explicit proofs.

Run from the repository root:

```bash
python3 -B scripts/run_nested_chain_checks.py
./build.sh
python3 tests/test_smoke.py
```

The report records every finite fixture, exact weights and allocations,
source hashes, and the independent checks. Gaussian elimination uses the
original rational residue columns; the deadline scheduler is a separate
oracle. Exact line arrangements check monotonicity, and a separate strip
enumeration checks (9)--(11). Complete small half-fans check every
positive facet against the consecutive-block cut certificate and (U2),
including lengths not divisible by $M$. Finite agreement is not the
uniform proof.
