# Five initial ideals can require quadratic integer weights

The number of initial ideals and the arithmetic size of a weight selecting
one of them are different parameters. The following all-field family has
exactly five initial ideals, but one cone requires integer weights of
quadratic height in the quotient length. Its smallest intersecting rational
projection also has the square-root dependence on primitive ray height
used in the [rational-shadow proof](rational_shadow_bound.md).

This is a uniform obstruction to two proposed shortcuts, not an improvement
to the extremal count $G_3(n)$ or a novelty certificate. The active goal's
sharp-growth or polynomial-separation endpoint is unchanged.

## The quotient and its marked basis

For every integer $k\ge2$ and every field $K$, put

$$
I_k=(xy-y^2,\ y^{k+1}-x^k,\ x^{k+1}-x^k,\ z-x^k)
       \subset K[x,y,z]. \tag{1}
$$

Use maximal weight for leading terms. The four displayed polynomials form
a reduced Gröbner basis whenever $w=(s,t,u)$ belongs to

$$
C_k=\{s>t,\ (k+1)t>ks,\ u>ks\}. \tag{2}
$$

These inequalities imply $t>0$, $s>0$, $u>0$: from $s>t$ and
$(k+1)t>ks>kt$ one gets $t>0$. The heads are
$xy,y^{k+1},x^{k+1},z$, and every tail is a standard monomial.
The only noncoprime head pairs give the S-polynomials

$$
x^{k+1}-y^{k+2},\qquad x^ky-x^ky^2.
$$

In each difference both terms reduce to $x^k$. The other head pairs are
coprime. Buchberger's criterion therefore proves the assertion in every
characteristic. In particular,

$$
D_k=\{1,x,\ldots,x^k,y,\ldots,y^k\},\qquad
n=\dim_K K[x,y,z]/I_k=2k+1. \tag{3}
$$

Conversely, selecting the indicated heads of the first, second and fourth
polynomials forces (2). Thus $C_k$ is the full open positive cone for
this initial ideal, not a small subcone created by unnecessary comparisons.

## Exactly five initial ideals

The monomial residue classes have a field-independent description:

* $X_a=\{x^a\}$ for $0\le a<k$;
* $Y_j=\{x^a y^b:b\ge1,\ a+b=j\}$ for $1\le j\le k$;
* $E$ consists of pure powers $x^a$ with $a\ge k$, mixed monomials
  $x^a y^b$ with $b\ge1$, $a+b\ge k+1$, and every monomial containing $z$.

All monomials in $E$ reduce to $x^k$. The relations give these reductions,
and the independent basis (3) proves that no listed classes are identified.
The residue vectors are consequently the different unit vectors of a
$(2k+1)$-dimensional space. A standard basis selects the cheapest monomial
in each class.

For positive weights, the minimum of $Y_j$ is $y^j$ if $s>t$, and
$x^{j-1}y$ if $s<t$. The minimum cost in $E$ is

$$
\min\{ks,(k+1)t,u\}. \tag{4}
$$

For $s>t$ all three choices in (4) are possible. For $s<t$, the middle
choice cannot win because $(k+1)t>ks$. There are therefore exactly five
distinct standard sets, realized by these strict weights:

| $Y_j$ representative | $E$ representative | Weight $(s,t,u)$ |
|---|---|---|
| $y^j$ | $x^k$ | $(k+2,k+1,(k+1)^2)$ |
| $y^j$ | $y^{k+1}$ | $(k+2,k,k(k+1)+1)$ |
| $y^j$ | $z$ | $(2,1,1)$ |
| $x^{j-1}y$ | $x^k$ | $(1,2,k+1)$ |
| $x^{j-1}y$ | $z$ | $(1,2,1)$ |

Every set also contains $1,x,\ldots,x^{k-1}$. The assumption $k\ge2$
makes the two $Y$ choices different. Boundary tie-breaking supplies no
further initial ideals: a global order has a positive weight representing
the finitely many strict comparisons in its reduced basis. Equivalently,
its comparison of $x$ and $y$ determines all $Y_j$ choices simultaneously.

For another description, the minimum basis-profile cost is

$$
\frac{k(k-1)}2s+kt+
\frac{k(k-1)}2\min(s,t)+\min(ks,(k+1)t,u). \tag{5}
$$

This is the lower support function of a translated segment plus a triangle.
Its restriction to positive objectives has the five cells above. The full
positive fan is thus explicitly classified, rather than estimated by a
sample of weight vectors.

This also exhibits what algebraic dependencies add to the original
smooth-prefix geometry. The standard set $D_k$ contains $x^2,y^2$ but
omits $xy$. Since the exponent of $xy$ is the midpoint of those of
$x^2,y^2$, no strict linear cutoff can include the two endpoints and
exclude $xy$. Thus $D_k$ is not a corner cut or a smooth-number prefix
for any prime bases. The quotient may skip $xy$ because its residue is
already represented by $y^2$. This changes admissible shapes, without
itself producing many different shapes or a prime-distribution theorem.

## Exact integer certificate height

The slack matrix of (2) is unimodular:

$$
V_k=\begin{pmatrix}1&-1&0\\-k&k+1&0\\-k&0&1\end{pmatrix},
\qquad \det V_k=1.
$$

The primitive rays of its closed cone are the inverse columns

$$
r_1=(k+1,k,k(k+1)),\quad r_2=(1,1,k),\quad r_3=(0,0,1). \tag{6}
$$

Every strict integer weight has three positive integer slacks. It is
therefore an integer combination $a r_1+b r_2+c r_3$ with $a,b,c\ge1$.
The coordinatewise smallest is their sum. Hence the exact minimum is

$$
\min_{w\in C_k\cap\mathbb Z^3}\|w\|_\infty
 =(k+1)^2=\frac{(n+1)^2}{4}, \tag{7}
$$

attained at $(k+2,k+1,(k+1)^2)$. In particular, even a unimodular
initial cone of an actual finite quotient need not contain a positive
integer weight of height $O(n)$.

For comparison, the classical cone and determinant argument gives, for
every length-$n$ quotient in $d$ variables and every initial ideal,

$$
\min\|w\|_\infty\le d!n^{d-1}. \tag{8}
$$

Indeed every standard monomial has total degree at most $n-1$, since its
divisor chain contains degree plus one different standard monomials. A
minimal generator of the initial ideal has degree at most $n$: all its
proper divisors are standard. The reduced-basis head-minus-tail vectors
thus have integer coordinates of absolute value at most $n$. Include
the positivity inequalities. The closure of the full-dimensional cone
is pointed. Every primitive extreme ray is given by cofactors of $d-1$
independent tight rows, so its coordinates are at most
$(d-1)!n^{d-1}$ in absolute value. The sum of $d$ independent primitive
rays is a strict positive integer interior point, of height at most
$d!n^{d-1}$. No facet can contain all those independent rays. For $d=1$
the assertion is immediate. This proves (8) using established rational-cone
geometry, rather than a new representation theorem.

For $d=3$, (7) and (8) show that the worst certificate-height exponent is
exactly two along infinitely many lengths, with the upper bound valid for
every length. This is not a pointwise lower claim at every even length,
nor a lower bound on coefficient arithmetic or algorithm running time.

## Exact rational projection radius

Let $q$ be a primitive nonzero integer vector. The plane $q^\perp$ contains
an objective strictly selecting $D_k$ exactly when $q^\perp\cap C_k$ is
nonempty. The exact minimum possible kernel height is

$$
\min_{q^\perp\cap C_k\ne\varnothing}\|q\|_\infty
 =\lfloor k/2\rfloor+1. \tag{9}
$$

To prove the lower bound, first suppose $q_3=0$. Then $q_1,q_2$ have
opposite signs and $s/t=|q_2|/|q_1|\in(1,1+1/k)$. Writing
$|q_2|=|q_1|+j$, with integer $j\ge1$, gives $|q_1|>kj$ and
$\|q\|_\infty\ge k+2$.

Otherwise change the sign of $q$ so that $q_3=-c$, $c\ge1$. With
$Q=\|q\|_\infty$ and $w\in C_k\cap q^\perp$,

$$
ck s<c u=q_1s+q_2t\le Q(s+t)<2Qs.
$$

Thus $Q>k/2$. For the matching upper witness take

$$
Q=\lfloor k/2\rfloor+1,\qquad q=(Q,Q,-1),\qquad
w=(2k+1,2k,Q(4k+1)).
$$

These give $q\cdot w=0$ and all three slacks in (2) positive.
This proves (9). The maximum primitive ray height in (6) is
$h=k(k+1)$, so the minimum radius is $\Theta(\sqrt h)$.
This also applies to the actual finite residue-matroid profile polytope
on $V_n^3=\{a:\prod_i(a_i+1)\le n\}$. That ground contains the four
heads and their tails in (1). Exchanging $y^2$ for $xy$, or $x^k$
for $y^{k+1}$ or $z$, gives three independent profile differences.
The polytope is full dimensional, and these comparisons force its
minimum normal cone at the profile of $D_k$ to be precisely $C_k$.
Conversely every objective in $C_k$ uniquely selects $D_k$ by the class
calculation. Its incident primitive facet normals therefore have
maximum height $h$, up to the sign convention for outward normals.
Consequently the square-root exponent in the rational-shadow visibility
lemma cannot be decreased uniformly, even for these finite quotients.
This does not rule out a smaller adaptively chosen *collection* of
projection planes, or a stronger aggregate estimate for large-normal
facets. It rules out a better universal radius for each individual cone.

## The optimal cover uses only two planes

For this same family, exactly two planes suffice to meet all five open
initial cones. The optimum maximum primitive kernel height of such a
cover is still $Q=\lfloor k/2\rfloor+1$.

For the upper bound use $q_0=(Q,Q,-1)$ and $q_1=(1,0,-1)$. In the order
of the five standard sets in the table, the following objectives meet
the indicated plane and strictly select the corresponding set:

| Plane | Objective |
|---|---|
| $q_0^\perp$ | $(2k+1,2k,Q(4k+1))$ |
| $q_0^\perp$ | $(2k+2,1,Q(2k+3))$ |
| $q_1^\perp$ | $(2,1,2)$ |
| $q_0^\perp$ | $(1,2k+2,Q(2k+3))$ |
| $q_1^\perp$ | $(1,2,1)$ |

For the lower bound on cover size, intersect any one plane through the
origin with the positive weight simplex. This is an interval, if
nonempty. Along that interval the affine part of (5) has no breaks;
$\min(s,t)$ has at most one, and the lower envelope of the three
affine functions in (4) has at most two. Each affine function's winning
set is an interval. Their sum consequently has at most three breaks,
and meets at most four open initial cones. A line lying in a tie wall
cannot add strict cones. Thus one plane cannot cover all five. Equation
(9) gives the lower bound on the maximum kernel height of every cover,
and the displayed cover attains it.

The cover cardinality is therefore constant while its necessary radius
grows linearly in $n$. An aggregate projection method can be much more
economical here than summing all integer directions up to that radius.
There is no conclusion yet that an equally economical cover exists for
arbitrary finite quotients.

## Prior-art boundary and the next counting question

[Fukuda--Jensen--Lauritzen--Thomas, *The generic Gröbner walk*,
arXiv:math/0501345](https://arxiv.org/pdf/math/0501345), Sections 2.2--2.3,
already provide finite positive-weight representation and the reduced-basis
cone inequalities in their rational-order exposition. Lemma 2.1,
Corollary 2.2 and the cone construction were
read with their proofs. Their introduction also explains why a symbolic
walk avoids constructing large perturbation weights. Neither that
representation principle nor symbolic avoidance is a new result here.

[Tran, *A fast algorithm for Gröbner basis conversion and its applications*
(2000)](https://www3.risc.jku.at/research/theorema/Groebner-Bases-Bibliography/gbbib_files/publication_627.pdf),
Lemma 3.1 and Theorem 3.1 with proofs, give degree-based constructions of
representing weights for lexicographic and general admissible orders.
Those statements were inspected, not the entire paper's applications.
Equation (8) is retained as a classical cofactor consequence, without a
new-algorithm claim. The precise priority of the explicit five-cone family
and exact minima (7), (9) is unestablished; these scoped comparisons do not
prove novelty.

For the active extremal problem, a large required integer weight or
projection height therefore cannot be used as evidence of a large fan.
The next upper route must control how many difficult cones can coexist in
one quotient, improving the aggregate $O(n^4/h)$ facet incidence or the
number of selected projections. The next lower route must count many
simultaneous dependencies outside the already bounded binomial class.
The family (1), itself binomial, supplies neither a power separation nor a
reason to start a larger unconstrained search.

## Independent checks

From the repository root run

```bash
python3 -B scripts/run_initial_cone_height_checks.py
./build.sh
python3 tests/test_smoke.py
```

The standard-library checker implements ordinary polynomial division and
all S-pairs independently of the monoid residue rule. It also constructs
integer commuting multiplication maps, verifies the defining equations on
every basis vector, compares matrix powers with polynomial normal forms,
and uses Gaussian selection and exact cone clipping to check the five
cells, their disjoint interiors and full simplex area. Separate integer
cube and primitive-kernel enumerations test the arithmetic minima. The
report saves the exact witnesses, scope and relevant source hashes at
`results/research/initial_cone_height_checks.json`. Finite agreement supports
the uniform proofs; it does not prove them or establish novelty.

The saved controls include 570 ordinary S-pairs over the rationals and
four prime fields, 3,665 independently reduced ground monomials, 11
complete five-cone fans over $\mathbb F_{101}$, 55,415 primitive kernel checks and 171,178
subcritical integer weights. Two-plane cover objectives are checked
against the quotient's Gaussian basis oracle as well as the cone
inequalities. The all-parameter assertions depend on the proofs above.

## The hard cone can be entirely interior

Boundary coverage does not remove the worst individual-cone obstruction.
For every $k\ge3$, change only the last relation of (1), putting

$$
\widetilde I_k=(xy-y^2,\ y^{k+1}-x^k,\ x^{k+1}-x^k,\ z-y^k).
\tag{10}
$$

Its length is again $n=2k+1$ over every field: the same three
relations in $x,y$ have the reference basis (3), and $z=y^k$.
Consider the standard set

$$
\widetilde D_k=\{1,x,\ldots,x^k,y,\ldots,y^{k-1},z\}.
$$

Its reduced basis has heads $xy,y^k,x^{k+1},xz,yz,z^2$, with respective
tails $y^2,z,x^k,x^k,x^k,x^k$. These reductions follow from (10);
their heads generate the complement of $\widetilde D_k$, which has
length $n$. Consequently strict selection of these heads gives this
initial ideal. The first, second and fifth comparisons are necessary,
and imply all the others. Its full open cone is exactly

$$
\widetilde C_k=\{s>t,\ kt>u,\ u>ks-t\}. \tag{11}
$$

The assumption $k\ge3$ ensures $y^2$ is standard and the first
comparison is necessary. Positivity follows from (11): it gives
$kt>ks-t$, hence $t>k(s-t)>0$, and $u>(k-1)t>0$.
The slack matrix and its inverse columns are

$$
\widetilde V_k=\begin{pmatrix}1&-1&0\\0&k&-1\\-k&1&1\end{pmatrix},
\qquad\det\widetilde V_k=1,
$$

$$
\widetilde r_1=(k+1,k,k^2),\quad
\widetilde r_2=(1,1,k-1),\quad
\widetilde r_3=(1,1,k). \tag{12}
$$

All three rays are strictly positive. Thus the closed simplex cell
meets no coordinate boundary, even at a vertex. Positive integer
slacks once more give the exact minimum certificate:

$$
\widetilde w_k=(k+3,k+2,k^2+2k-1),\qquad
\min_{w\in\widetilde C_k\cap\mathbb Z^3}\|w\|_\infty
 =k^2+2k-1=\frac{(n+1)^2}{4}-2. \tag{13}
$$

The minimum primitive kernel height of a plane meeting this cone is
exactly $k$. For the upper bound use
$q=(k,k-1,-2)$ and $\widetilde w_k$ from (13).
For the lower bound suppose $Q=\|q\|_\infty\le k-1$.
If $q_3=0$, the interval $1<s/t<1+1/k$ excludes such an integer
relation, as in the proof of (9). Otherwise orient $q$ so that
$q_3=-c$, $c\ge1$. Then

$$
c(k-1)s<cu=q_1s+q_2t<2Qs\le2(k-1)s,
$$

so $c=1$. Set $a=s-t>0$ and $S=q_1+q_2\in\mathbb Z$.
The two remaining inequalities are

$$
(k-1)t+ka<q_1a+St<kt.
$$

If $S\ge k$, the upper inequality forces $q_1<0$ and hence
$q_2\ge k+1$. If $S\le k-1$, the lower inequality forces $q_1>k$.
Both contradict $Q\le k-1$. This proves the lower bound.

For completeness this quotient has exactly nine initial ideals.
For $1\le j<k$, its $Y_j$ class from the earlier calculation chooses
$y^j$ or $x^{j-1}y$ according as $s>t$ or $s<t$.
The $Y_k$ class also has representative $z$, and the final class $E$
has minimum cost

$$
\min\{ks,(k+1)t,s+u,t+u,2u\}.
$$

When $s>t$, the possible $(Y_k,E)$ representatives are
$(y^k,x^k)$, $(y^k,y^{k+1})$, $(z,x^k)$, $(z,yz)$, $(z,z^2)$.
When $s<t$, they are $(x^{k-1}y,x^k)$, $(z,x^k)$, $(z,xz)$,
$(z,z^2)$. Each occurs on an open positive region, by choosing $u$
on the indicated side of the thresholds and $s/t$ sufficiently close
to, or sufficiently far from, one as required. The two lists give
distinct standard sets since $k\ge3$. Tie-breaking gives no additional
initial ideal. Only the cone (11) is entirely interior.

Hence even after boundary cells are covered, the quadratic integer
certificate exponent and linear projection-radius exponent persist
for an actual finite quotient with constant fan size. This does not
refute a smaller adaptively selected cover or a cubic weighted-incidence
budget. It is a further obstruction to replacing those aggregate
questions by a uniform radius bound on the remaining cells.
The exact replay is in `scripts/run_algebraic_incidence_controls.py`.
