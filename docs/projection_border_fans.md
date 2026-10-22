# Coordinate-quotient fans bound coupling on the surface yz = 0

Let $A=K[x,y,z]/I$ have finite length $n>0$, over any field, with
$yz\in I$. Put

$$
 A_Y=A/zA,\quad A_Z=A/yA,\quad n_Y=\dim A_Y,\quad n_Z=\dim A_Z,
 \quad p_Y=g(A_Y),\quad p_Z=g(A_Z).
$$

The quotient $A_Y$ is a planar algebra in $x,y$, and $A_Z$ is a
planar algebra in $x,z$. Write $g(0)=1$ for the empty standard set.
There are absolute constants such that

$$
 \boxed{g(A)\le C n^2 p_Yp_Z.}                            \tag{1}
$$

Consequently all such quotients with bounded $p_Yp_Z$ have only
$O(n^2)$ initial ideals. This order is sharp in characteristic zero,
even when **both coordinate quotients have just one initial ideal**.
For $p_Yp_Z=n^{o(1)}$ the count is $n^{2+o(1)}$, so coupling these
quotients cannot yield a polynomial separation above the generic
corner-cut scale. A family with $g(A)\ge c n^{2+\epsilon}$ requires
$p_Yp_Z\ge c' n^\epsilon$, and hence at least one of the two
actual planar quotient fans must grow polynomially.

The proof applies with arbitrary overlap $yA\cap zA$, arbitrary
coefficients, nonreduced structure and growing coordinate degrees.
It imposes no row-polynomial hypothesis on the branch annihilators.
The quotient ideals $I+(z)$ and $I+(y)$ are essential: merely
restricting $I$ to a coordinate plane computes a different object.

This is a uniform conditional class bound and an exclusion of
bounded-complexity component constructions. The unrestricted $G_3$
bound is unchanged. The general planar bound only gives a possibly
quadratic product $p_Yp_Z$, so (1) alone does not close the nodal
surface problem. Its exact quantitative priority is unassessed;
the proof uses elementary ideal multiplication and the established
lattice-volume theorem. It does not complete the active research goal.

## A border lift, not an assumed uniform rank rule

Fix any positive global monomial order and let $D,B_Y,B_Z$ be the
standard exponent sets of $A,A_Y,A_Z$. Quotient maps give

$$
                     B_Y\cup B_Z\subseteq D.             \tag{2}
$$

There is a second inclusion that controls every extra monomial:

$$
 D\subseteq \{1,x,\ldots,x^{a-1}\}\cup yB_Y\cup zB_Z,
 \qquad a=\dim K[x]1\le n.                              \tag{3}
$$

Indeed $y(I+(z))\subseteq I$ since $yz\in I$. If a monomial $u$
is in the initial ideal of $I+(z)$, choose a polynomial there with
head $u$. Multiplying it by $y$ gives a polynomial of $I$ with
head $yu$. Thus a standard monomial $yu$ of $A$ must have $u\in B_Y$.
The $z$ argument is identical. Mixed $yz$ monomials vanish, and a
monic $x$-annihilator of degree $a$ supplies the axis bound. This
does not use generic coefficients, a support-matching rank rule,
or identification of a vector contraction with an algebra quotient.

Set $C=B_Y\cup B_Z$ and

$$
 E=\bigl(\{1,x,\ldots,x^{a-1}\}\cup yB_Y\cup zB_Z\bigr)\setminus C.
                                                               \tag{4}
$$

All $y$-positive labels of $E$ lie in $(yB_Y)\setminus B_Y$.
For each fixed $x$-exponent, a planar downset has one vertically
maximal point. Multiplication by $y$ moves it to height equal to
that column's length. Summing over the columns gives exactly

$$
       \sum_{u\in E}u_y=n_Y,\qquad
       \sum_{u\in E}u_z=n_Z.                             \tag{5}
$$

The two parts are disjoint, and the axis adds zero to both sums.
This is a budget over the entire candidate set, not just one
selected basis. It is why arbitrary dependencies among these
candidates cannot create unlimited local fan complexity.

## One fixed represented matroid on each pair of planar cones

Normalize $w_x=1$, $u=w_y/w_x>0$, $v=w_z/w_x>0$. Every open
initial cone of $A_Y$ is an interval in $u$; every cone of $A_Z$
is an interval in $v$. Their products give precisely $p_Yp_Z$
nonempty open regions. On one such region, $B_Y,B_Z,C,E$ are fixed.
Let

$$
                  k=n-|C|.
$$

The actual columns $C$ are independent, since (2) puts them in
every parent standard basis there. Contract them in the represented
residue matroid on $C\cup E$. It has residual rank $k$ and captures
the entire parent fan patch. Every parent optimum is available by
(2)--(3), and no lift of another residual basis can beat the original
global minimum. The compulsory span need not be an ideal; that is
not required for this vector-matroid argument.

All labels in (4) have $x$-exponent at most $n-1$: the quotient
downsets have sizes at most $n$, and $a\le n$. The convex hull $Q$
of residual basis profiles therefore satisfies

$$
 Q\subseteq[0,k(n-1)]\times[0,n_Y]\times[0,n_Z].           \tag{6}
$$

Its bases need not be coherent or downsets. Equation (5) still
controls their last two coordinates because each is a subset of $E$.
Distinct parent standard sets with a strict witness in this region
give distinct exposed vertices of $Q$. The full-dimensional lattice
polytope $Q+[0,1]^3$ has at least as many vertices and volume at most
$(1+k(n-1))(n_Y+1)(n_Z+1)$. The classical Andrews bound in dimension
three gives

$$
 \#\{\text{parent patches in the region}\}
    \le C\sqrt{(1+k(n-1))(n_Y+1)(n_Z+1)}.                \tag{7}
$$

For $k=0$ there is exactly one possible parent basis $C$. For all
$k\le n$, (7) is $O(n^2)$. Every full parent cone meets one of the
open product regions; finitely many quotient walls cannot cover an
open cone. Summing proves (1). A parent cone may meet several
regions; this only overcounts.

There is an exact overlap interpretation of $k$. Write

$$
 h=\dim(yA\cap zA),\quad d=\dim A/(yA+zA),\quad
 r=|B_Y\cap B_Z|.
$$

The intersection of these two monomial sets is a pure-axis prefix.
Its length is the minimum of their axis lengths and is at least $d$.
The dimension identities give

$$
                   k=h+r-d.                              \tag{8}
$$

Thus (7) can be summed with the actual values $h+r-d$ instead of
the coarse $k\le n$. Overlap need not be zero, and the residual rank
is not in general just $h$. Dropping the extra axis term $r-d$
would omit actual dependencies even when the two branch ideals are
disjoint.

## Sharp amplification from trivial quotient fans

The existing [Pascal family](univariate_branch_gluing.md#an-explicit-quadratic-fan-and-its-prime-power-collapse)
has $A=K[t]/t^{2m}$ with
$x=t$, $y=t^m(1+t)^m$, $z=t^m(1+t)^{2m}$. It satisfies $yz=0$,
and $yA=zA=t^mA$. Both coordinate quotients are $K[t]/t^m$ with
$y=z=0$; hence $p_Y=p_Z=1$. The already proved binomial-Vandermonde
calculation gives

$$
                       g(A)=\binom{m+2}{2}
$$

in characteristic zero. Since $n=2m$, this attains the quadratic
power in (1). The same saved family has cubic primitive wall cost.
Together with (1), the general bound $\nu(q)\le n$ per actual
wall and the three-dimensional edge bound give sharp cubic order
for the bounded-$p_Yp_Z$ class. A reduced-basis union has at most
$O(n^{8/3})$ polynomials at fixed $p_Yp_Z$, by the usual minimal-head
packing. No new Pascal identity, construction, minimum universal-basis
bound, or actual edge-length estimate is claimed here.

The lower example is characteristic zero. The previous prime-power
collapses in positive characteristic remain: the all-field upper
theorem does not imply that this one coefficient family attains it
in every field.

There is also a whole nonrow, zero-overlap construction excluded by
the theorem. For arbitrary $R,S\ge2$, $1\le k<R$, $1\le l<S$,
take the product of two curve jets with coordinate maps

$$
 \begin{array}{c|ccc}
 &x&y&z\\\hline
 K[t]/t^R&t&1+t^k&0\\
 K[s]/s^S&1+s&0&1+s^l.
 \end{array}                                             \tag{9}
$$

The $x$ minimal polynomials $x^R,(x-1)^S$ are coprime, so $x$
generates the product by Chinese remainders. The displayed branch
coordinates are units on their own components, giving exactly the
two indicated coordinate quotients and zero branch intersection.
Each planar quotient has two initial ideals. To verify this, translate
its coordinates to $x=t$, $Y=t^k$; translations preserve global
initial monomials. Nonzero monomials of exponent $a+kb<R$ represent
one valuation class. Within a class their cost changes monotonically
with $b$ according to the sign of $w_Y-kw_x$, so just two standard
sets occur. The second component is identical after translating $x$.
This proof works in every characteristic. Thus (9) has $g=O((R+S)^2)$
uniformly in both growing curve degrees, without a generic rank rule.
The finite report's counts 23 in characteristic zero and 19 in
characteristic three at $(R,S,k,l)=(5,4,2,3)$ only illustrate why
the product of the two planar fan sizes is not itself a parent count;
no asymptotic lower bound is inferred from them.

## Prior art and the remaining proof problem

The ideal inclusions behind (2)--(3) are elementary and belong to
standard Gröbner-basis and colon-ideal reasoning. Matroid greedy
optimization, fixed independent-column contraction and the state-profile
dictionary are the established
[Babson--Onn--Thomas framework](https://arxiv.org/pdf/math/0207135).
The quantitative geometric input is Andrews' lattice-volume theorem;
[Dillon, Theorem 1](https://travisadillon.com/documents/volume_vs_vertex_in_lattice_polytopes.pdf)
states the exponent used after cube regularization. The
[row-branch proof](row_branch_fans.md) already uses a linear transverse
budget for a different fixed quotient, and the minimum-overlap proof
uses the same telescoping mechanism at canonical row endpoints.
No new general geometric or matroid theorem is asserted.

For zero overlap, $I=(I+(y))\cap(I+(z))$, so the algebra is the
usual fiber product of its planar projections over their common
axis quotient. The general pullback identity is recorded in
[Ananthnarayan--Avramov--Moore, equations (1.0.2)--(1.0.3)](https://pi.math.cornell.edu/~frankmoore/Papers/connected.pdf),
before their later local and Gorenstein hypotheses. This pullback
description is classical. It does not
make the coordinate quotient fans a sufficient description of the
parent fan, as the sharp amplification above demonstrates for the
broader class with overlap. The fixed candidate matroid still has
the actual residue coefficients and the axis term in (8).

Shifted Popov/interpolation algorithms compute one supplied shift;
the existing [all-shift comparison](oracle_fan_enumeration.md) explains
why that does not give the required all-weight count. A full priority
assessment of (1), (7) and their application remains open. Failed
searches, or checking a single supplied-shift theorem, do not establish
novelty.

This advances the feasibility boundary: arbitrarily chosen coupling
cannot produce a polynomial separation while both actual coordinate
quotient fans have subpolynomial complexity. The remaining problem
includes growing planar quotient fans and their simultaneous border
choices. Summing $O(n^2)$ independently over every quotient pair
does not achieve the desired unrestricted growth. The original
$G_3$ objective, strongest saved bound and success criteria remain
unchanged.

## Reproducible independent checks

```bash
python3 -B scripts/run_projection_border_checks.py
./build.sh
python3 tests/test_smoke.py
```

The [report](../results/research/projection_border_checks.json) constructs
the actual coordinate-ideal quotients from multiplication matrices,
checks their column maps, enumerates both complete quotient fans and
the parent fan, and intersects every joint cone with every parent
cone using exact rational clipping. It checks full area coverage,
(2)--(5), the rank identity (8), all residual moment bounds, and
independent residual ranks and greedy replays on every nonempty patch.

Examples cover nonrow branches with zero and nonzero overlap,
the sharp Pascal family with one-initial coordinate quotients,
characteristic-dependent coefficients, direct products of two
monomial-curve jets with two-initial planar projections, and zero
coordinate quotients. Every computed patch has a saved exact standard
set and rational positive witness. Finite verification supports the
implementation; the uniform proof is above. Neither finite checks
nor the classical inputs establish novelty or complete the goal.
