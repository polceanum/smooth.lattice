# A complete graph-rank model for prime-power-minus-one rectangles

There is a second tractable characteristic-dependent class inside the
original nonlinear rectangle family. Its side lengths are one **below**
a power of the characteristic, rather than powers themselves. The whole
monomial residue matroid is cographic. This gives an exact rank formula
for every subset and reduces weighted initial standard sets to maximum
spanning trees of a specified complete bipartite graph.

The graph representation and matroid duality used below are classical.
This is an explicit application to these local polynomial quotients, not
a new general spanning-tree algorithm, a proved counting separation, or
a verified novelty claim. The unrestricted goal remains open.

## The quotient and all its linear relations

Let $q_1=p^a,q_2=p^b\ge3$ and let $K$ have characteristic $p$.
After swapping the two variables if necessary, assume $q_1\le q_2$.
Thus $q_1\mid q_2$. Consider

$$
I_{q_1,q_2}=(x^{q_1-1},y^{q_2-1},z-x-y-xy),\quad n=(q_1-1)(q_2-1).          \tag{C1}
$$

As in the [bilinear translation lemma](bilinear_rectangle_fans.md),
$X=1+x,Y=1+y,Z=1+z$ preserves every global initial ideal. The quotient
becomes

$$
K[X,Y]/((X-1)^{q_1-1},(Y-1)^{q_2-1}),\qquad Z=XY.             \tag{C2}
$$

Frobenius gives $X^{q_1}=Y^{q_2}=1$ in (C2). Thus every monomial has the
same column as a residue $e_{uv}=X^uY^v$, $0\le u<q_1$, $0\le v<q_2$.
Moreover,

$$
(X-1)^{q_1-1}=1+X+\cdots+X^{q_1-1},                         \tag{C3}
$$

and similarly for $Y$. Start in the $q_1q_2$-dimensional group algebra
$K[X,Y]/(X^{q_1}-1,Y^{q_2}-1)$. The additional relations in (C2) span precisely
the $q_2$ row sums and $q_1$ column sums of $(e_{uv})$. Multiplication
by $X$ or $Y$ permutes these sums, so their linear span is already an
ideal. It has dimension $q_1+q_2-1$: it is the row space of an oriented
vertex-edge incidence matrix of the connected graph $K_{q_1,q_2}$, after
changing signs on one part. This argument works in characteristic two
as well. The quotient dimension is $q_1q_2-(q_1+q_2-1)=n$, as required.

For completeness, the incidence rank statement needs no characteristic
restriction. The sum of all its rows is zero; after one row is deleted,
the columns of any spanning tree form a triangular matrix with diagonal
entries $\pm1$ after leaf elimination. Thus the remaining $q_1+q_2-1$
rows are independent over every field.

## Exact rank of every residue subset

Identify $(u,v)$ with the edge joining the two corresponding vertices
of $G=K_{q_1,q_2}$. For $S\subseteq E(G)$, let $c(E\setminus S)$ count
connected components of the complementary spanning subgraph, including
isolated vertices. Then the columns in (C2) have rank

$$
\operatorname{rank}(S)=|S|+1-c(E\setminus S).              \tag{C4}
$$

Here is a direct proof. A relation supported on $S$ is a cut-space
vector whose coordinates vanish on every complementary edge. Its vertex
potential is therefore constant on each component of $E\setminus S$.
There are $c(E\setminus S)$ independent component potentials, with
one global constant giving zero cut. The relation space supported on
$S$ has dimension $c(E\setminus S)-1$, proving (C4). For monomial
lists, first merge repeated residue labels; duplicates are parallel
columns and change no rank.

Consequently $S$ is independent exactly when the complementary graph
is connected. It is a basis exactly when its complement is a spanning
tree. This is the usual cographic matroid, obtained here from commuting
polynomial multiplication. It is a description of all column subsets,
unlike the characteristic-zero prefix-only theorem.

## Weighted initial sets by maximum spanning trees

Fix a strict positive weight $w=(\alpha,\beta,\gamma)$. For each
residue $(u,v)$ choose its least-cost monomial representative $m_{uv}$.
It suffices to compare the $q_2$ candidates

$$
X^{(u-c)\bmod q_1}Y^{(v-c)\bmod q_2}Z^c,
\quad 0\le c<q_2.                                        \tag{C5}
$$

Any larger exponent can be reduced using the two periods in the group algebra,
which strictly lowers its positive cost. The
[divisible-period representative lemma](bilinear_rectangle_fans.md) improves
(C5) to a constant number of comparisons per residue. Write
$\ell_{uv}=w\cdot\operatorname{exp}(m_{uv})$.

Let $T$ be a maximum-weight spanning tree of $K_{q_1,q_2}$ with edge
weights $\ell_{uv}$. Then the initial standard set of (C1) is exactly

$$
\{m_{uv}:(u,v)\notin T\}.                                \tag{C6}
$$

Indeed a least-cost monomial basis must use the least representative
of each of its labels. By (C4), its labels are complements of spanning
trees. Minimizing the selected total is equivalent to maximizing the
complementary total. This is precisely minimum-weight matroid basis
selection, hence gives the global initial standard set. Translation
preserves its monomial exponents. At tied weights a refining monomial
order is needed; (C6) asserts the unique state only for strict weights.

With the divisible-period representative rule, a standard comparison-based
spanning-tree algorithm computes (C6) in $O(q_1q_2\log(q_1+q_2))=O(n\log n)$
comparisons and $O(n)$ storage, counting arithmetic comparisons on the
given weights. This is not a bit-complexity bound or a bound for writing
every coefficient of a reduced Gröbner basis.

## Precise established boundary and the next reduction

The underlying cographic identity is established, rather than a proposed
new matroid theorem. See Hidaka--Itoh,
[Theorem 1.2 and Section 5.2 with its following remark](https://arxiv.org/html/2403.10817v3#S5.SS2),
and Fu--Ma--Wang,
[Proposition 3.5, equation (3.4)](https://arxiv.org/html/2609.29227v1#S3.SS2).
The corresponding cyclotomic-prime construction is in Martin--Reiner,
[Theorem 1 and Remark 6](https://www-users.cse.umn.edu/~reiner/Papers/cyclotomic.pdf).
These specific primary passages were inspected; the full papers were
not audited. The [ordered continuation](cographic_order_ideals.md)
extends the polynomial realization to all fields and separates the
remaining coherence count from that known dictionary. For equal periods
it proves a two-vertex branching constraint and a linear-size weighted
candidate graph. No all-period fan count or novelty conclusion follows.

## What this does and does not permit next

This family has an exact rank rule over the full weight orthant. Counting
its initial ideals becomes a parametric spanning-tree problem with the
**specific periodic shortest-representative costs (C5)**. General
parametric spanning-tree lower bounds cannot be imported unless their
weight functions are realized by these costs. The number of trees of
$K_{q_1,q_2}$ is not the number realized by three monomial weights.

The next substantive test is to classify these maximum trees across
the representative chambers, or prove that this structure excludes a
polynomial separation. It is unnecessary to run arbitrary coefficient
searches or infer a growth law from a few fans. The classical graph-rank
identity alone does not establish publication significance.

[Eppstein's primary text](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf),
Introduction and Section 2, studies a single parameter with independently
given affine element weights. Its uniform-matroid/k-set equivalence and
its general graphic examples do not provide the required realization
by (C5), or a two-parameter count for this fixed family. Those introductory
statements were inspected as text; the full lower-bound proofs were not
audited here. This is a scope comparison, not a strongest-current-bound
assessment.

```bash
python3 -B scripts/run_bilinear_cographic_checks.py
```

The checker compares (C4) with direct Pascal-column elimination,
exhaustively for $(q_1,q_2)=(3,3)$ and on deterministic subsets for
larger equal and unequal powers.
It independently checks row/column relations, shortest representatives,
maximum-tree optimality, and standard sets from original local polynomial
multiplication. Small graph cases also enumerate all spanning trees.
These are finite controls for the proof, not a count of every fan.
