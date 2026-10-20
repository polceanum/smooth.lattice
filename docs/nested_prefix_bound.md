# A uniform bound for the whole additive nested family

The [explicit quotients](nested_jet_matroids.md) all have

$$
g(I)=O\bigl(n^2\log^4(n+1)\bigr)                              \tag{1}
$$

in three variables, uniformly in their positive integer valuations and
ambient truncation. Here $n$ is the **actual quotient length**. This includes
unequal valuations and the non-binomial quotients of embedding dimension
three. It supersedes the remaining-count statements for this family in
the earlier construction and [chain-allocation note](nested_chain_allocation.md).

The proof is a more general matroid statement. Its key point is that a
tied feasible exchange forces a whole cheap valuation prefix into the
chosen basis. A small box follows in the original exponent coordinates;
no continued-fraction analysis of the arithmetic strips is needed.
Classical matroid geometry and hyperplane-arrangement counting complete
the argument. Priority of the exact statement has not been established.
The unrestricted $G_3(n)$ endpoint remains unresolved here. This closes
a construction route; it does not complete the research goal.

## The abstract hypothesis

Fix $d\ge2$, a positive vector $m\in\mathbb R_{>0}^d$, and $L>0$.
Let

$$
T=\{a\in\mathbb N_0^d:m\cdot a<L\}.
$$

Let $M$ be a rank-$r\ge1$ matroid on this entire finite ground. Require
the following property, including equality of valuations:

$$
\begin{gathered}
B\text{ a basis},\quad e\in B,\quad f\in T\setminus B,
\quad m\cdot f\le m\cdot e\\
\Longrightarrow\quad B-e+f\text{ is a basis}.                 \tag{2}
\end{gathered}
$$

This is the familiar downward basis property of a Schubert matroid in
the valuation order, with symmetry within each tied valuation class.
It is a hypothesis about all bases, not an assertion about arbitrary
representable matroids or arbitrary polynomial quotients. Sorting an
independence criterion consisting of upper bounds on ordered valuations
immediately gives (2).

For $w>0$ put $c_w(a)=w\cdot a$ and minimize the total cost of a
basis. We count the distinct **uniquely optimal** bases on open regions.
At a degenerate weight there can be several optimal bases; they are not
counted separately as new open regions.

Every optimal basis is a monomial downset. If $a\in B$ has a proper
divisor $b\notin B$, (2) makes $B-a+b$ a basis with strictly smaller
positive cost. Thus every selected monomial has at most $r$ divisors,
and every positive optimum lies in the usual hyperbolic cross
$\prod_i(a_i+1)\le r$. This observation does not bound the box of
the whole selected set by $O(r)$.

## The cheap-prefix lemma

**Lemma 1.** Suppose $B$ is optimal at $w>0$, and $B-e+f$ is a
basis with $e\in B$, $f\notin B$ and
$c_w(e)=c_w(f)=t$. Put

$$
Q=\max(m\cdot e,m\cdot f),\qquad
C=\{a\in\mathbb N_0^d:m\cdot a\le Q,\ c_w(a)\le t\}.
                                                                  \tag{3}
$$

Then $t>0$, $Q<L$, both $e,f$ belong to $C$, and every point of
$C$ of cost strictly less than $t$ belongs to $B$. In particular,

$$
\#\{a\in C:c_w(a)<t\}\le r-1.                              \tag{4}
$$

For every $a\notin B$ with $m\cdot a\le Q$, the replacement
$B-e+a$ is a basis. If $Q=m\cdot e$, this follows directly from
(2). Otherwise apply (2) to the basis $B-e+f$, replacing $f$ by
$a$; the case $a=f$ is the given feasible exchange. Optimality of
$B$ therefore forces $c_w(a)\ge c_w(e)=t$. This proves the
containment. The selected $e$ has cost exactly $t$, so at most $r-1$
selected points can lie strictly below it. Positivity and $e\ne f$
exclude $t=0$, and both endpoint valuations are below $L$.

The argument uses the full prefix below the higher endpoint valuation,
including chains earlier than the active allocation block. That extra
containment is what removes the arithmetic-strip obstacle left by the
previous proof.

## Small relative boxes for every active direction

**Lemma 2.** Under Lemma 1, the closed cut satisfies

$$
|C|\le2^d(r-1),\qquad
\prod_i A_i\le d^d|C|\le(2d)^d(r-1),                       \tag{5}
$$

where $A_i=1+\max_{a\in C}a_i$ are its coordinate counts. Consequently
the raw exchange difference and its primitive integer direction obey

$$
\prod_i(|f_i-e_i|+1)\le(2d)^d(r-1).                        \tag{6}
$$

Coordinatewise floor-halving preserves $m\cdot a\le Q$ and sends
$c_w(a)\le t$ strictly below $t$, since $t>0$. It has fibers of
size at most $2^d$, proving the first bound. Both inequalities in (3)
have nonnegative coefficients, so $C$ is a downset. Its coordinate
extrema occur on the axes. The box

$$
0\le a_i\le\lfloor(A_i-1)/d\rfloor\quad(1\le i\le d)
$$

lies in $C$: each of the $d$ coordinate contributions to either
inequality is at most $1/d$ of its threshold. It contains at least
$\prod_i A_i/d^d$ points. The difference of two points of $C$ has
$|f_i-e_i|+1\le A_i$, and division by its integer gcd cannot increase
the product. This proves (5)--(6).

In three variables the constant in (6) is $216$. It depends on the
rank, not the number of ground monomials, $L$, $m$, or the coefficient
sizes of a realization. Large boxes for individual coherent bases remain
compatible with these small boxes for active exchanges.

## The fan bound and one common weight list

**Theorem 3.** Matroids satisfying (2) have at most

$$
O_d\bigl(r^{d-1}\log^{(d-1)^2}(r+1)\bigr)                  \tag{7}
$$

uniquely optimal bases for positive additive monomial costs. All their
positive fans of a fixed rank are refined by a single arrangement
depending only on $d,r$.

The minimum-basis cone at $B$ is defined by the inequalities
$c_w(f)-c_w(e)\ge0$ over feasible single-element replacements.
This is classical matroid polytope geometry: edges of a matroid base
polytope are single-element exchanges, and checking all improving edges
is sufficient for a linear minimum. Thus every facet meeting $w>0$
has an active feasible exchange. By Lemma 2 its primitive normal belongs to

$$
\mathcal H_d(r)=\{h\in\mathbb Z^d\setminus\{0\}:
                    \prod_i(|h_i|+1)\le(2d)^d r\}.          \tag{8}
$$

For fixed $d$, this set has $O_d(r\log^{d-1}(r+1))$ elements.
For example, summing successively over $|h_i|+1$ bounds the count by
a constant times $r$ times $d-1$ harmonic sums. A central arrangement
of $s$ planes in dimension $d$ has $O_d(s^{d-1})$ regions. The
planes $h\cdot w=0$ for (8), restricted to $w>0$, refine every
positive basis cone and give (7). Rank one has a single positive optimum
and needs no active exchange.

There is a useful uniform consequence beyond an individual count.
Choose one positive integer weight from each full-dimensional chamber
of this common arrangement. The resulting list $W_d(r)$ samples all
positive uniquely optimal bases for **every** $m,L,M$ satisfying (2),
without adapting the weights to the chosen matroid. It has the size in
(7), and weights can be chosen with

$$
\|w\|_\infty=O_d(r^{d-1}).                                  \tag{9}
$$

To see the height bound, include the coordinate planes. Each chamber
closure is a full-dimensional pointed cone. Each extreme ray is the
intersection of $d-1$ independent planes whose integer coefficients
are $O_d(r)$, hence has an integer generator of size $O_d(r^{d-1})$
by cofactors. Sum $d$ independent ray generators. Every chamber
inequality is strict on at least one summand, so the sum lies in its
interior and satisfies (9). This gives a finite shared test list and
an oracle enumeration consequence. It is not a claim of a practical
implementation or priority over existing universal-basis algorithms.

## Application to every explicit nested quotient

For the construction in $d$ variables, set

$$
B=\max(L,2),\quad C=L\sum_{i=0}^{d-1}B^i,\quad
\lambda_i=Cm_i+B^{i-1},\quad
X_i=t^{m_i}e^{\lambda_i t}\pmod{t^L}.                       \tag{10}
$$

Here $L,m_i$ are positive integers. The digit encoding and the ordered
Hermite--Birkhoff argument in the original note work without a restriction
to three variables. The nonzero monomials have $m\cdot a<L$, and
$s$ distinct columns with increasing valuations $v_i$ are independent
exactly when

$$
v_i\le L-s+i-1\quad(1\le i\le s).                         \tag{11}
$$

Lowering one valuation preserves these inequalities, including ties.
The true image length is

$$
n=L-\max_{0\le k\le L}\bigl(k-\#\{a:m\cdot a<k\}\bigr).   \tag{12}
$$

The monomial residue matroid therefore satisfies (2) with rank $r=n$.
Theorem 3 gives (1) for every member in three variables and the bound
in (7) in every fixed dimension. It applies over every characteristic-zero
field; it is not an
all-characteristics realization claim. The abstract matroid theorem itself
does not assume representability over any field.

Every global monomial initial ideal of a finite quotient is represented
by a positive weight on its finite basis and border comparisons. Thus the
count covers the full set of global monomial initial ideals, and the
common list $W_d(n)$ covers the whole explicit family. The existing
generic subfamily $m=(1,\ldots,1),L=n$ supplies the generic corner-cut
lower bound. Together these bounds fix the polynomial growth exponent
$d-1$ for this construction class, with a logarithmic gap remaining for
$d\ge3$. The sharp generic lower result is recorded in the existing
[corner-cut manuscript](../paper/corner_cut_growth/main.tex); its priority
assessment is separate.

This closes both the unequal curvilinear valuation route and the
embedding-dimension-three subfamily as sources of a power separation.
It does not cover all curvilinear embeddings, all local ideals, or all
residue matroids whose vector-space representations happen to be nested
after an arbitrary relabeling. The additive exponent labeling and (2)
are essential.

## Generic fixed-order power series have the same bound

The conclusion covers more than the displayed exponential coefficients.
Let

$$
Y_i=t^{m_i}+\sum_{j=m_i+1}^{L-1}c_{i,j}t^j\pmod{t^L},       \tag{13}
$$

with $Y_i=0$ when $m_i\ge L$, and regard the visible higher coefficients
as parameters over a characteristic-zero field. There is a nonempty
Zariski-open set of these coefficients on which **every** monomial rank
is given by (11), the image length is (12), and (7) bounds all global
monomial initial ideals. This includes generic curvilinear embeddings
with a specified tuple of leading orders and generic finite subalgebras
of a one-parameter truncation with arbitrary positive leading orders.
Nonzero leading coefficients can be normalized by scaling the variables,
which preserves all monomial ranks and monomial initial ideals.

For any coefficients, the column of a monomial of valuation $v$ has zero
entries below row $v$. The dimension argument in the original rank
theorem therefore gives the same necessary inequalities for independence.
There are finitely many subsets of $T$ satisfying them. For each such
subset, choose a matrix minor nonzero at the explicit exponential point
(10); its determinant is a polynomial in the coefficients of (13).
The product of all these chosen minors is a nonzero polynomial. Where
it is nonzero, all the necessary rank inequalities are also sufficient.
The chosen point has leading coefficients one, so this argument already
works in the normalized coefficient space (13). It proves the assertion
over every characteristic-zero field without a random search.

The exceptional set matters. For instance, in truncation $L=5$ with
orders $(1,2,3)$, write
$Y_1=t+u t^2+\cdots$ and $Y_2=t^2+v t^3+\cdots$. The minor of
the columns $Y_1^2,Y_2$ on rows $t^2,t^3$ is $v-2u$.
At the specialization $Y_i=t^{m_i}$ those two full columns coincide,
whereas the ordered exponential choice makes them independent.
A power-separated family in this coefficient model would have to leave
the generic rank locus. Such a specialization can also lower the image
length, so its length must be recomputed. The two-cone example below
shows why a failure of generic rank alone is not enough.

## A concrete boundary of the method

Large primitive exchange products do occur in finite quotients. For
coprime $p,q\ge2$, over any field consider

$$
J_{p,q}=(z,xy,x^{p+1},y^{q+1},x^p-y^q).                    \tag{14}
$$

Its length is $n=p+q$. Besides one common socle line, the nonzero
monomial residue lines are $1,x,\ldots,x^{p-1},y,\ldots,y^{q-1}$;
the socle is represented by either $x^p$ or $y^q$. These lines are
independent. The ideal has exactly two monomial initial ideals,

$$
(z,xy,x^p,y^{q+1}),\qquad (z,xy,x^{p+1},y^q),
$$

separated by $p w_x=q w_y$. The primitive normal $(p,-q,0)$ has
product $(p+1)(q+1)$. For $p=M,q=M+1$, this is quadratic in the
length. At $M=512$ it already exceeds $216(n-1)$.

No positive additive valuation can give (2) on a full sublevel ground
containing both socle representatives and all their lower valuations.
At least one of $p m_x,q m_y$ is at least $m_x+m_y$: the two contrary
strict inequalities would imply
$m_y>(p-1)m_x$ and $m_x>(q-1)m_y$, which are incompatible.
A basis containing that nonzero socle monomial cannot replace it by
the lower-valued zero monomial $xy$. This is an explicit failure of
the hypothesis, not evidence of many initial ideals: (14) has only two.

The next lower construction must combine failures of this single-valuation
order across many coherent choices in one fixed algebra. The preserved
nonlinear rectangles provide a defined candidate with overlapping
dependencies and an existing $O(n^{9/4})$ square-case upper bound.
Their support-matching shortcut is already disproved, so a uniform proof
must use the actual ranks or another exact algebraic description. A new
finite coefficient screen or a normal-height example would not meet the
original counting goal.

## Prior-art scope and reproducible checks

[Benedetti-Velásquez--Knauer, Section 2.1 and Remark 3](https://link.springer.com/article/10.1007/s00493-024-00085-4)
describe the established Schubert, shifted and nested matroid terminology
and the ordered Gale-basis definition. Their order convention is reversed
relative to (2). Those definitions and the paper's stated main results
were inspected; the present argument does not claim a new matroid class.

The classical edge-direction input is restated as
[Arkani-Hamed--Lam--Spradlin, Proposition 2.1](https://link.springer.com/article/10.1007/s00220-021-04041-x)
with attribution to Gelfand--Goresky--MacPherson--Serganova. That proposition
and the accompanying Schubert definitions were inspected. It concerns
the base polytope before projection; the cheap-prefix estimate (4)--(6)
is an additional use of the additive lattice labeling. Hyperbolic-cross
counting and the central-arrangement bound are classical.

The shared-weight-list principle is already in
[Babson--Onn--Thomas, Definition 3.1, Proposition 3.2, Theorems 3.5 and 4.2](https://arxiv.org/pdf/math/0207135).
These statements and their proofs were inspected for this continuation.
Their Hilbert zonotope uses primitive directions from differences of two
rank-$n$ hyperbolic-cross points and gives a common list of size
$O_d(n^{2(d-1)}\log^{2(d-1)^2}n)$ for all finite quotients.
The present restricted list uses the smaller direction set (8).
Thus the common-list construction and universal-basis enumeration method
are prior art; the additional statement to assess is the prefix bound
under hypothesis (2) and its generic fixed-order-jet application.
Their Theorem 2.4 and Corollary 2.5 remain the earlier general
all-quotient comparison in the existing extremal audit. A newly checked
[2026 F4 paper abstract](https://arxiv.org/abs/2603.16378) concerns generic
grevlex computation, conditional on Moreno--Socías, rather than all
initial ideals. Only its abstract was inspected; it is not a replacement
for a complete strongest-current-bound audit. The exact priority of this
restricted theorem and the unrestricted $G_3$ endpoint remain unsettled.

From the repository root:

```bash
python3 -B scripts/run_nested_prefix_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker tests the lowering hypothesis independently with rational
matrices on small full grounds, reuses the eight saved complete fans,
and checks every positive facet against the full cheap-prefix certificate.
Additional low-dimensional matroids test the dimension-independent
statement and a small common weight list. The two-cone algebra (14)
is checked by sparse commuting multiplication operators. Generic
fixed-order coefficient controls also compare full subset ranks and an
explicit exceptional minor. The report
contains exact witnesses, relevant source hashes and commands. These
finite controls support the uniform proof; they do not establish novelty
or complete the unrestricted research goal.
