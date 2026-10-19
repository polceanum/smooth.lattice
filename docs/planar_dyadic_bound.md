# A dyadic bound for planar matroid profiles

This note gives a uniform proof of an $O(n\log(n+1))$ bound on the number
of global monomial initial ideals of every colength-$n$ ideal in two
variables, over every field. The known generic planar corner-cut lower
bound has the same order. The argument applies more generally to
parametric matroid optimization with integer weights in a hyperbolic
cross. It uses greedy exchanges and an elementary distinct-direction
estimate; it does not require the earlier lattice-shadow argument or
the unresolved reading of McMullen's 1975 lattice application.

The proof and exact controls are supplied here. The precise novelty
comparison remains under audit. A failed keyword search is not evidence
of priority. The three-variable target remains open: this planar theorem
closes the proposed planar power-amplification route rather than
establishing the original spatial endpoint.

## Statement for matroids

Let $M$ be a finite rank-$r$ matroid on $E$, and let
$a(e)=(x_e,y_e)\in\mathbb Z_{\ge0}^2$ satisfy

$$
(x_e+1)(y_e+1)\le H,\qquad H\in\mathbb Z_{\ge1}.
$$

Repeated weight vectors are allowed with a fixed tie order. We count
the distinct basis profiles exposed by objectives in an open quadrant,
or the canonical greedy bases on its open parameter intervals. We do
not count every tied optimum at a degenerate objective. Put
$L=\lceil\log_2 H\rceil$, with $L=0$ when $H=1$.

**Theorem.** In each open quadrant the number of profiles is at most

$$
1+L\left(r+3(r^2H)^{1/3}\right). \tag{1}
$$

Consequently the full basis-profile polygon has at most four times
this many vertices. In particular, for $r=H=n$ the positive-quadrant
count is at most

$$
1+4n\lceil\log_2 n\rceil. \tag{2}
$$

The number of ground elements does not appear in (1). Representation
over a field is not assumed. For polygon vertices, matroid duality
also permits replacing $r$ by $\min(r,|E|-r)$: complementary basis
profiles are related by central reflection, which takes a quadrant
to its opposite.

## A monotone exchange path

First use objectives $t x+y$, $t>0$. Away from finitely many crossing
parameters, increasing-cost greedy selection gives a minimum basis.
At a crossing $t_0$, elements with the same cost form consecutive tie
blocks. Immediately before the crossing their order within each tie
block is decreasing in $x$; immediately after, increasing in $x$.
Reverse each block by adjacent swaps.

The standard adjacent-swap property of matroid greedy selection says
that swapping adjacent elements $e,f$ either changes nothing or
replaces the selected $e$ by the previously unselected $f$. To see
this directly, let $P$ be the prefix before the pair and contract its
closure. If the pair is independent, both are selected; if both are
loops, neither is selected; if just one is a loop, the other is selected.
The only changing case is two nonloops that are parallel in this
contraction; then the first of them is selected. Their combined span
after the pair is unchanged, so all later greedy decisions are unchanged.

Every successful swap as $t$ increases satisfies

$$
x_e>x_f,\qquad y_e<y_f,\qquad
\frac{y_f-y_e}{x_e-x_f}=t_0. \tag{3}
$$

Put one token on each element of the initial basis and transfer the
token from $e$ to $f$ at an exchange. There are always $r$ tokens.
Along each token path $x$ strictly decreases and $y$ strictly
increases. Simultaneous crossings cause no problem: all elementary
swaps at $t_0$ have the same ratio in (3). A fixed order among
identical weight vectors never changes.

This is the usual greedy exchange/token viewpoint in parametric
matroid optimization; compare
[Eppstein, Section 4](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf).
The dyadic charging below, not the existence of tokens, is the
quantitative step being investigated here.

## Exchanges between dyadic blocks

Give $x=0$ label zero. For $j=0,\ldots,L-1$, give

$$
2^j\le x<2^{j+1}
$$

label $j+1$. These labels cover every element because $x\le H-1$.
The sum of the token labels is initially at most $rL$, never
increases, and drops by at least one at every exchange between
different blocks. Thus there are at most $rL$ such exchanges, and
at most that many crossing parameters containing one.

There are no within-block exchanges for $x=0$, since the relative
costs there do not depend on $t$. For a positive block $j$, each
token visits the block at most once. The total $x$ decrease of that
token during within-block exchanges is at most $2^j-1$. Also
$y<H/2^j$ throughout the block, so its total $y$ increase in those
exchanges is at most $H/2^j$. Summing over tokens gives budgets

$$
A_j\le r2^j,\qquad B_j\le rH/2^j,\qquad A_jB_j\le r^2H.
\tag{4}
$$

The budgets include only within-block exchanges. Entry and exit
exchanges are already charged to the label potential. A token
cannot leave the block and later return, because its $x$ coordinate
is monotone. This is the point where the common matroid exchange
path matters; unrelated pairs of bases would not give (4).

## An elementary integer-direction estimate

**Lemma.** Suppose $K$ distinct positive slopes have positive integer
vectors $(u_i,v_i)$ in those directions, with
$\sum_i u_i\le A$ and $\sum_i v_i\le B$. Then

$$
K^3\le27AB. \tag{5}
$$

For $K=0$ the assertion is immediate. Otherwise $A,B>0$. Replace
each vector by its primitive positive integer vector, without
increasing either sum. Distinct slopes give distinct primitive
pairs. The sum of $u_i/A+v_i/B$ is at most two. At least $K/3$
pairs therefore satisfy $u_i/A+v_i/B\le3/K$. All those pairs lie in
the positive integer rectangle

$$
1\le u\le3A/K,\qquad1\le v\le3B/K,
$$

which has at most $9AB/K^2$ lattice points. This proves (5).
This is an elementary lattice-chain estimate, not a claim of a new
lattice-polygon vertex theorem.

For each block $j$ and each crossing parameter having within-block
exchanges, sum the positive integer $x$ decreases and $y$ increases
of those exchanges. Equation (3) makes the result a positive integer
vector of slope $t_0$. Different crossing parameters give different
slopes. Their coordinate sums obey (4), so (5) gives at most
$3(r^2H)^{1/3}$ such parameters per block. Every actual basis change
contains either a between-block exchange or a within-block exchange.
Adding the bounds and the initial state proves (1) in the positive
quadrant.

For an arbitrary open quadrant use objectives
$\sigma t x+\tau y$, $t>0$, with $\sigma,\tau\in\{-1,1\}$.
The same adjacent swaps make $x$ monotone in the direction
$-\sigma$ and $y$ monotone in direction $\tau$. Each token still
visits each $x$ block at most once, and the absolute coordinate
changes have ratio $t_0$. The label potential is reversed when
$\sigma=-1$; its total range is still $rL$. The proof of (4)--(5)
is unchanged. This proves all four versions of (1). Every vertex
of a planar polygon, segment or singleton is exposed in an open
quadrant by some objective, so four such counts bound the full
polygon. For $H=1$ all weights are zero and there is one profile.

## Application to every finite planar ideal

Let $I\subset K[x,y]$ have colength $n\ge1$. Its standard monomials
for any global monomial order form an $n$-element downset. If
$x^a y^b$ is standard, all its $(a+1)(b+1)$ divisors are standard.
Every initial standard set is therefore contained in

$$
V_n^2=\{(a,b)\in\mathbb Z_{\ge0}^2:(a+1)(b+1)\le n\}.
$$

The residues of these monomials give a rank-$n$ linear matroid over
$K$. Scanning them in monomial order and selecting independent
columns gives exactly the initial standard set: a nonstandard
monomial reduces to earlier standard monomials, and standard
columns are independent.

Every order on this finite ground induced by a global monomial order
has a strict positive weight representation. One usual proof applies
the strict theorem of alternatives to its exponent comparisons and
the coordinate unit vectors. A nonnegative rational dependence
among those vectors could be cleared to integer coefficients; it
would give equal monomials on the two sides of a product of strict
monomial inequalities, a contradiction. A rational positive
solution can be perturbed to avoid all ground-element ties.
These are the established universal-ground and greedy/state
principles of
[Babson--Onn--Thomas, Proposition 2.2 and Theorem 2.4](https://arxiv.org/pdf/math/0207135).

Distinct monomial initial ideals have distinct standard sets.
Apply (2) with $r=H=n$ to obtain the all-fields bound

$$
g(I)\le1+4n\lceil\log_2 n\rceil. \tag{6}
$$

Generic $n$-point ideals realize all plane corner cuts. The classical
number of these is $\Theta(n\log(n+1))$, due to
Corteel--Rémond--Schaeffer--Thomas (1999); a precise primary
restatement is
[Elizalde--Galván, Theorem 2.3](https://arxiv.org/html/2312.16353v1).
For the generic-ideal dictionary, take independent point-coordinate
parameters. Every maximal evaluation minor on $n$ distinct monomials
is a nonzero polynomial: different determinant permutations give
different monomials in the independently labelled point coordinates.
There are only finitely many such minors on $V_n^2$, so they can be
made simultaneously nonzero over an infinite field. Greedy selection
then chooses precisely the first $n$ monomials in every strict positive
weight order.
Combining that known lower bound with (6) gives

$$
F_2(n)=\Theta(n\log(n+1)),\qquad
G_2(n)=\Theta(n\log(n+1)). \tag{7}
$$

This is an asymptotic statement, not equality of the two maxima
for each $n$, nor a statement that generic ideals maximize the
finite count. The published seven-point example with ten initials
versus eight generic cuts remains valid.

The link to the original smooth-number observations is the same
positive weighted exponent lattice. Algebraic dependencies can
increase particular finite counts, but (7) shows they cannot give
a polynomial gain in this two-variable setting. No consequence for
the distribution of prime numbers is asserted.

## Consequences beyond the count

**Universal bases.** A planar monomial ideal of colength $n$ with
$m$ minimal generators has $n\ge m(m-1)/2$. Sort its generators
by decreasing $x$ and increasing $y$ exponents. Each has total
degree at least $m-1$, since those exponents are distinct nonnegative
integers. Hence all monomials of total degree at most $m-2$ are
standard. Thus $m\le\sqrt{2n}+1$. A reduced basis has one polynomial
per minimal initial generator, and taking the union over (6) gives

$$
|\mathcal U(I)|=O(n^{3/2}\log(n+1)). \tag{8}
$$

The generator packing fact and union argument are established; the
new quantitative input under investigation is (6). Equation (8)
counts polynomials, not arithmetic operations or coefficient bits.

**A realization obstruction.** A planar quotient of length $N$
preserving $\Omega(M^{4/3})$ different parametric bases must have

$$
N=\Omega\!\left(\frac{M^{4/3}}{\log(M+1)}\right). \tag{9}
$$

Indeed $N\log(N+1)=\Omega(M^{4/3})$ by (6); if
$N<M^{4/3}$, substitute $\log(N+1)=O(\log(M+1))$, and otherwise
the claim is immediate. This rules out the proposed realization
with colength $O(M^\alpha)$ for any fixed $\alpha<4/3$, even
before imposing its extra coefficient and multiplication requirements.

**The spatial product route is closed at the power scale.** In the
common-coordinate construction from the main note, planar branch
counts multiply while colengths add. Equations (6)--(7) imply that
every fixed number of such planar branches gives at most
$O_d(N^{d-1}\log^{d-1}(N+1))$ states at spatial length $N$.
Consequently it cannot produce a polynomial gain over the generic
scale. The three-variable goal requires another mechanism, or a
uniform upper proof in three variables. It has not been completed
by (7).

## Prior-art boundary and verification status

| Input or comparison | Scope |
| --- | --- |
| Onn--Sturmfels, *Cutting Corners* (1999), Corollary 5.4 and Proposition 5.6 | The point-configuration maximum and its historical $O(n^{4/3})$ planar bound; generic maximality already fails at seven points. |
| Babson--Onn--Thomas (2003), Proposition 2.2, Theorem 2.4 and Corollary 2.5 | The universal monomial ground, greedy basis dictionary and historical all-ideal $O(n^{4/3})$ planar count. |
| Corteel--Rémond--Schaeffer--Thomas (1999), as restated in Elizalde--Galván Theorem 2.3 | The known matching lower order from generic plane corner cuts. It does not itself state the upper bound for arbitrary quotient matroids. |
| Eppstein, Section 4 | Classical greedy token paths and general parametric matroid examples; the dyadic product-budget application requires its own priority comparison. |
| [Hausbrandt--Ruzika (2025), Corollary 3.4](https://arxiv.org/html/2503.08178v1#S3) | A bound $O(m^{2p})$ for unrestricted $p$-parameter matroids on $m$ elements, from the pairwise comparison arrangement. It does not state the hyperbolic integer-weight bound (1); this one comparison is not a complete priority audit. |

The new argument has a complete written proof and exact controls.
Its current-priority audit is incomplete; the inspected primary
statements above are comparisons, not a guarantee that no later
equivalent theorem exists. The original three-variable objective
and its paper requirements remain active. No author contact,
submission or private correspondence is authorized.

Reproduce the exact controls from the repository root:

```bash
python3 -B scripts/run_planar_dyadic_bound_checks.py
```

The saved report contains **108 quadrant fixtures** with 314 states,
including rank zero, full rank, repeated weights, simultaneous crossings,
linear and uniform matroids, a basis-exchange-verified sparse-paving
matroid, one-block convex chains and actual nonreduced and reduced
quotients. Every comparison interval is replayed by direct greedy
selection, independently of the token updates. Small full profile hulls
are also checked by enumerating every basis.

An exhaustive control covers all **813 labelled rank-two matroids on six
elements**, with 3,252 quadrant sweeps and independently enumerated full
profile hulls. A separate quotient control checks all **526 nonempty point
subsets of $\mathbb F_2^2$ and $\mathbb F_3^2$**, totaling 810 positive
initial bases. The report saves exact witnesses and ledger digests, verifies
the per-token dyadic budgets and primitive-direction inequality, and
checks downward closure for quotient standard sets. These finite controls
support the uniform proof; they do not establish novelty or the spatial
endpoint.
