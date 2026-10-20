# Exact ranks for truncations of two-part hypersurfaces

The shared depth bound leaves room for a polynomially larger fan when
relations propagate through a growing number of layers. This note
makes one such class explicit. It gives an exact homogeneous initial
ideal and a finite free-module description of **all** monomial ranks.
For the explicit cubic, the subsequent [uniform classification](cubic_fan_classification.md)
now proves the complete fan count and a matrix-free rank query over every
field. The more general two-part class is not classified. The free-module
and graded-initial mechanisms are classical; novelty, the main extremal
goal and its publication gates remain open.

Let $S=K[x,y,z]$, $\mathfrak m=(x,y,z)$, and

$$
f=g+h,\qquad 1\le d=\deg g<e=\deg h,
\qquad \gcd(g,h)=1,
$$

where both parts are homogeneous and nonzero. For $s\ge0$ set
$I_s=(f)+\mathfrak m^{s+1}$. All ordinary initial forms below retain
the **highest** ordinary degree, matching global orders. They are
different from the lowest-degree tangent cone at the origin.

## A uniform decomposition, before choosing any monomial order

After a field extension if necessary, choose a linear form $L$ for
which $(g,h,L)$ is a homogeneous system of parameters. Choose
homogeneous lifts $t_1,\ldots,t_{de}$ of a vector basis of
$S/(g,h,L)$, with degrees $\delta_j$. Then

$$
S=\bigoplus_{j=1}^{de}t_j K[g,h,L],\qquad
\sum_j t^{\delta_j}=
\frac{(1-t^d)(1-t^e)}{(1-t)^2}.                          \tag{1}
$$

This is the classical homogeneous free-module construction. To
recall its justification, coprimality makes $g,h$ a regular sequence
in the polynomial ring. A general $L$ completes it to a regular
system of parameters. Graded induction gives generation by the
chosen lifts over $K[g,h,L]$. The regular-sequence Hilbert-series
identity shows that the resulting graded surjection has zero kernel.
It also gives the degree polynomial in (1). Field extension preserves
all original monomial-column ranks and Gröbner bases; it introduces
no extra initial ideals for the original ideal.

**Proposition 1.** The following is an exact vector basis of $S/I_s$:

$$
\{t_j g^a L^c:\ e a+c+\delta_j\le s,\quad a,c\ge0\}.
                                                               \tag{2}
$$

In particular,

$$
n_s=\dim_K S/I_s=
\binom{s+3}{3}-\binom{s-d+3}{3},                         \tag{3}
$$

where the second term is zero for $s<d$. For fixed $d$, $n_s$ has
order $s^2$. The cutoff coefficient in (2) is the **higher degree**
$e$, even though (3) depends on the lower degree $d$.

**Proof.** Modulo $g+h$, a free monomial $t_j g^a h^b L^c$
maps to $(-1)^b t_j g^{a+b}L^c$. For fixed $k=a+b$, its largest
possible ordinary degree is $\delta_j+e k+c$, attained at $a=0$.
Thus the image of $\mathfrak m^{s+1}$ is exactly the span of those
free labels with $\delta_j+e k+c>s$. Distinct remaining labels are
independent in the free module modulo $g+h$. This proves (2).
Taking the generating function in $s$ gives

$$
\sum_{s\ge0}n_s t^s=
\frac{\sum_j t^{\delta_j}}{(1-t)^2(1-t^e)}
=\frac{1-t^d}{(1-t)^4},
$$

which proves (3). $\square$

If $g$ is not a monomial, the canonical monomial part of $I_s$ is
exactly $\mathfrak m^{s+1}$. Otherwise a monomial below the cutoff
could equal a truncated multiple of $f$; its lowest homogeneous
part would be divisible by the nonmonomial $g$, impossible in a
polynomial ring. For $s\ge d$, the parent relation module has
coordinate annihilation depths exactly $s-d+1$: that power kills
all terms, while multiplying $f$ by $x_i^{s-d}$ leaves a nonzero
degree-$s$ part. Thus the depth product is $\Theta(s^3)$ for fixed
$d$, and the preceding depth bound leaves an $n_s^{9/4}$ window.
This class is not excluded by the bounded-depth theorem.

This description preserves the rank of every selected set of original
monomials: expand each one in (1), replace $h$ by $-g$, and discard
exactly the labels excluded by (2). Take the rank of the resulting
columns. No claim that the original monomials themselves belong to
individual free summands is involved.

## The global degree initial ideal is explicit

Write $\mathfrak m^q=S$ for $q\le0$. Then

$$
\boxed{\operatorname{in}_{\rm deg}(I_s)=
(h)+\sum_{k=0}^{\lceil(s+1)/e\rceil}
         g^k\mathfrak m^{\max(0,s+1-e k)}.}               \tag{4}
$$

The right side is homogeneous. To prove inclusion, use
$g^k\equiv(-h)^k\pmod f$: every displayed multiple of $g^k$
belongs to $I_s$ because its corresponding $h^k$ multiple has
degree at least $s+1$. Also $h$ is the highest-degree initial
form of $f$. For equality, use (1) after setting $h=0$.
In the resulting free basis $t_j g^a L^c$, the displayed ideal
kills precisely the labels with $e a+c+\delta_j>s$.
The choice $k=a$ proves one direction; for any other $k\le a$,
the lower bound on its multiplier's degree proves the other.
Consequently its colength is (3). Inclusion and equal colength
prove (4).

For instance, for $g=x^2+y^2+z^2$, $h=xyz$ and $s=5$,
(4) includes $g^2$, of degree four. Merely retaining $h$ and
discarding terms above degree five would miss that relation.
The identity $g^2\equiv h^2\pmod f$ explains it uniformly.

Formula (4) describes refinements near the ordinary-degree ray.
The earlier [positive-grading bound](nonnegative_graded_fans.md)
already supplies a linear bound in $n_s$ for that local star.
Neither fact counts cells elsewhere in the positive weight orthant.

## An explicit six-generator model

Consider

$$
f=x^2+y^2+z^2+xyz,
$$

over a field of characteristic different from two. Set
$U=x^2+y^2+z^2$, $L=x+y+z$. The polynomial ring is free over
$K[U,H,L]$, where $H=xyz$, with basis

$$
1,y,y^2,x,xy,xy^2.                                      \tag{5}
$$

Here is a direct verification, rather than an assumed normal-form
rule. Eliminate $z=L-x-y$. The two monic relations are

$$
\begin{split}
x^2&=-xy-y^2+Lx+Ly-\tfrac12L^2+\tfrac12U,\\
y^3&=Ly^2-\tfrac12(L^2-U)y+H.
\end{split}                                             \tag{6}
$$

Their leading monomials for the auxiliary lex order $x>y>U,H,L$
are $x^2,y^3$, which are relatively prime. They form a monic
Gröbner basis and prove (5). Modulo $f$, replace $H$ by $-U$.
Proposition 1 becomes the concrete basis

$$
\mathcal B_s=\{x^i y^j U^a L^c:
 i<2,\ j<3,\quad i+j+3a+c\le s\},\qquad
|\mathcal B_s|=(s+1)^2.                                 \tag{7}
$$

Multiplication by $x,y,z=L-x-y$ is now given by the bounded
recurrences (6), with $H=-U$. Terms beyond (7) vanish. Each
recurrence raises its filtration level by at least one, so the
discarded part is stable under all three variables. The original
equation and original monomial columns can be checked directly.
This uses $O(s^2)$ module coordinates, compared with the
$\binom{s+3}{3}$ coordinates of the full monomial parent.
Each of the three multiplication matrices has a bounded number of
nonzero entries per column, independently of $s$: (6) reduces only
$x^2$ and $y^3,y^4$ in a single multiplication step. Thus the
operators themselves require $O(n_s)$ stored field coefficients.
It is an exact research rank model, not a new general FGLM or
Gröbner-basis algorithm or a claimed running-time improvement.

The variables $U,L$ are polynomial functions of the original
coordinates. The six generators in (5) are therefore **not six
common price rays** for the original monomial optimization.
For example $z=L-x-y$ already has three components. A rank basis
change does not turn the cost $a w_x+b w_y+c w_z$ into a cost
on one free-module label. The earlier common-ray theorem cannot
be applied to infer a fan count from (7).

## A smaller integral model

For this particular cubic, classical monic division is simpler than
the parameter construction (5)--(7). It also removes the restriction on
the characteristic. More generally, suppose

$$
F=x^d+\sum_{a<d} c_{abc}x^a y^b z^c,
\qquad a+b+c\ge d
$$

for every nonzero displayed coefficient, over any coefficient ring.
The usual monic normal form has $x$-exponent below $d$. Each
replacement of $x^d$ strictly decreases $x$-exponent and never
decreases total degree. Therefore the image of $\mathfrak m^{s+1}$
in that free module is precisely the span of its standard monomials
of degree greater than $s$. This proves, over the coefficient ring,

$$
S/((F)+\mathfrak m^{s+1})
\quad\text{has basis}\quad
\{x^i y^j z^k:i<d,\ i+j+k\le s\}.                    \tag{9}
$$

No division by a field element is required. In particular, for

$$
f_c=x^2+y^2+c z^2+xyz,\qquad c\in\mathbb Z,
$$

(9) gives $(s+1)^2$ original monomials over $\mathbb Z$.
The only multiplication rule beyond shifts and the cutoff is

$$
x(xy^jz^k)=-y^{j+2}z^k-cy^jz^{k+2}-xy^{j+1}z^{k+1}.    \tag{10}
$$

All three multiplication matrices have at most three nonzero entries
per column, so they use $O(n_s)$ integer coefficients. The module has
two channels over $\mathbb Z[y,z]$ before truncation. This supersedes
the six-channel construction as a rank implementation for this cubic;
the previous identities and checks remain valid. It is a classical
monic-division simplification, not a new mathematical complexity bound.
Original monomials with $x$-exponent at least two still have mixed
columns. Nothing here supplies independent original-cost rays.

## Where the unknown count is located

For any $f$ as above, suppose its strict leading monomial $u$
has ordinary degree $d$. Then

$$
\operatorname{in}_w(I_s)=(u)+\mathfrak m^{s+1}.           \tag{8}
$$

Indeed the right side is contained in the initial ideal and has
exactly the colength (3). This proves equality. It counts these
regions completely, without a generic-coefficient assumption.

For the explicit cubic and $s\ge3$, (8) gives exactly the three outer cones
$w_x>w_y+w_z$ and its permutations. All further cells lie in
the open triangle defined by the three strict triangle inequalities,
where the leading monomial of $f$ is $xyz$. The saved examples
have extra walls inside that triangle. Their existence is not
a polynomial lower rate.

The classical profile-volume bound is still only $O(n_s^{9/4})$
for this fixed-degree family: all standard monomials have total
degree at most $s$, so profiles lie in a cube of side $n_s s$,
and $n_s=\Theta(s^2)$. Thus a real quantitative window remains.
The [subsequent proof](cubic_fan_classification.md) closes the explicit
cubic's whole-triangle count at linear growth, using its coupled original
costs. Both (7) and the smaller integral kernel (9)--(10) remain valid
rank models. The classical window is not attained by this cubic.

## A precise boundary-cover question

The reduction below records the previously unresolved question. Its
endpoint inequality is now proved for every cutoff and field in the
[classification note](cubic_fan_classification.md), by matching two exact
price-prefix rank formulas. The conditional argument is retained here
to explain the role of the finite certificates.

Set

$$
J_s=(x^2+y^2+xyz)+\mathfrak m^{s+1},\qquad v=(1,1,0).
$$

For the original cubic $I_s$, the highest-$v$ initial ideal is
exactly $J_s$. The generators give inclusion, and (9) proves equality
of colengths. This equality and the standard refinement lemma imply
that each positive initial cone of $J_s$ lifts to a cone of $I_s$
incident to $v$. They do **not** imply that those incident cones
cover the positive triangle.

By symmetry it is enough to consider $w_x>w_y>w_z>0$ in the triangle.
After rescaling, write

$$
w=(u+\lambda,u,1),\qquad u>1,\quad 0<\lambda<1.         \tag{11}
$$

$J_s$ is homogeneous for $v$. Its degree-$k$ quotient summand has
dimension $2(s-k+1)$ for $1\le k\le s$, and the degree-zero summand
has dimension $s+1$. All of its monomial bases have the same total
$v$-degree. Thus its greedy standard basis in (11) depends only on
$\lambda$, not $u$.

For a direct independent description of each summand, put $t=x/y$
as a **vector-space label**, without inverting the nilpotent variable
in the original quotient. For $L=s-k+1$ its monomial columns are

$$
x^a y^{k-a}z^c\longmapsto t^a z^c
\quad\text{in}\quad
K[t,z]/(t^2+zt+1,z^L),\qquad 0\le a\le k.               \tag{12}
$$

The basis is $1,t$ over $K[z]/z^L$. Formula (12) follows from the
monic relation $x^2=-y^2-xyz$ and preserves the entire degree-$k$
column space. It is a finite block-rank identity, not localization of
$S/J_s$.

Let $B$ be a greedy basis of $J_s$ on an open interval
$(\ell,h)\subset(0,1)$, and write the reduced original border
relations as $x^a-\sum_{b\in B}c_{ab}x^b$. Here $a$ runs over the
minimal excluded border. For each nonzero coefficient put
$z=a-b=(z_x,z_y,z_z)$ and $A=z_x+z_y$. The refinement at $v$
already proves $A\ge0$. The following explicit condition is
sufficient for this **entire** sector, including arbitrarily large $u$:

$$
\boxed{ A+\ell z_x+z_z\ge0,\qquad
        A+h z_x+z_z\ge0\quad\text{for every border tail}.} \tag{13}
$$

Indeed its cone inequality is

$$
z\cdot w=A(u-1)+(A+\lambda z_x+z_z).                     \tag{14}
$$

The second term is affine in $\lambda$, so its endpoint values
control the whole interval. If $A>0$ the first term is strictly
positive for $u>1$. If $A=0$, the corresponding tail is a border tail
of the homogeneous initial $J_s$, and its inequality is strict in
the open interval. Thus all these relations retain their designated
heads. Inclusion of the border ideal and equality of colengths give
the original initial ideal with standard set $B$.

**Conditional consequence.** If (13) holds for every $s$ and every
sector basis over every field, the whole cubic family has only
$O(n_s)$ positive monomial initial ideals. To see this directly,
all sector endpoints come from the boundary border normals,
$z_x+z_y=0$, and have the form $-z_z/z_x$. Standard monomials have
degree at most $s$, and minimal border monomials have degree at most
$s+1$. There are only $O(s^2)$ such distinct fractions in $(0,1)$.
Permuting coordinates covers the six ordered central regions, and
(8) adds just three outer cones. Since $n_s=(s+1)^2$, the conclusion
follows. The nonnegative-grading theorem gives another bound of the
same order for each boundary star. This is a conditional reduction,
not a proof that the cubic family satisfies the hypothesis uniformly.

There is also a uniform parity restriction on the boundary walls:
$z_x+z_z$ is even. Both the $v$-grading and the parity of the sum of
the $x$- and $z$-exponents are preserved by the relation in (12).
Consequently a primitive fraction $p/q$ with opposite parities in
$p,q$ must be represented by an even multiple of $(q,-q,-p)$.
This explains a restriction on candidate Farey fractions; it does
not say that every permitted fraction is an actual wall.

## Complete finite certificates, over every field

The checker now verifies (13) on every boundary sector at the
following cutoffs. This is a continuous cone certificate, not a mesh
of sampled weights.

| Cutoff $s$ | Length $n_s$ | Ordered central sectors | Entire positive fan |
|---|---|---|---|
| 3 | 16 | 1 | 6 |
| 4 | 25 | 2 | 12 |
| 5 | 36 | 2 | 12 |
| 6 | 49 | 3 | 18 |
| 8 | 81 | 5 | 30 |
| 10 | 121 | 7 | 42 |
| 12 | 169 | 10 | 60 |
| 14 | 225 | 12 | 72 |

These finite counts hold over **every field**, including characteristic
two. Every minimal border relation was obtained with integer tail
coefficients and leading coefficient one. An explicit integer
quotient $Q$ verifies that the relation minus $f_cQ$ has only terms
of degree greater than $s$. Reducing those identities into any field
therefore retains all designated heads. Each standard set has the
known length $(s+1)^2$, so inclusion proves equality of the initial
ideals in every characteristic. The same monic argument certifies
the boundary sectors. The saved border generators, rational interval
endpoints, critical normal witnesses and reproduction command
determine the exact finite certificates.

The intervals exhaust $(0,1)$: a next test point is taken just past
the exact preceding endpoint, with step below the minimum possible
gap between fractions of denominators at most $s+1$. Formula (13)
covers every $u>1$. Coordinate permutations and the outer-cone
classification finish the positive fan; lower-dimensional weight
boundaries cannot hide an additional full-dimensional initial cone.
All previous complete fan counts at cutoffs three and four agree.

The arbitrary-cutoff statement is proved in the separate classification
note. That proof excludes the entire full-support diagonal cubic family
as a source of a polynomial separation and supplies an exact rank oracle.
The table is retained as independently reproducible support; it did not
prove the uniform result. Unrestricted $G_3$ and the active goal remain open.

## Prior-art boundary and reproduction

[Lauritzen, Sections 2--3, Lemma 2.1, Proposition 2.3 and
Theorem 3.1 with proofs](https://arxiv.org/pdf/math/0509247)
uses a homogeneous ideal and a truncating degree set. Its theorem
identifies the fan with that of the corresponding homogeneous
ideal plus a monomial truncation. The inhomogeneous ideal $(g+h)$
does not meet that hypothesis; (4) records the additional relations
in this setting. No general truncation-commutation claim is inferred
outside the homogeneous setting.

[Bahloul--Takayama, Section 1 definitions and Section 3 comparison](https://www.numdam.org/item/10.1016/j.crma.2006.11.030.pdf),
and [their full preprint, Theorem 5.1.1 and proof](https://arxiv.org/pdf/math/0412044v2),
use nonpositive weights for the power-series fan. The latter
comparison assumes a positive grading and compares fans on that
nonpositive domain. Neither statement by itself supplies a bound
for all positive global-weight cells of $I_s$. These exact scope
comparisons are not a novelty certificate for (2), (4) or (7).
Their classical free-module ingredients and complete subsequent
priority coverage still need assessment before any paper claim.

The cone inequalities and refinement mechanism are established
Gröbner-fan theory. In [Fukuda--Jensen--Thomas, Proposition 2.6,
Lemma 2.10, Corollary 2.11, Lemma 2.13 and Corollary 2.14, with proofs](https://arxiv.org/pdf/math/0509544),
the reduced-basis inequalities describe a cone and its closure, and
initials of its basis give refinements. Those passages were read;
the general fan-convexity argument in (13)--(14) is not claimed as
a new theorem. Its application still needs the unproved uniform
border inequality.

After $z\mapsto-z$, the cubic itself is the classical Markoff surface
$x^2+y^2+z^2-xyz=0$. [Goldman, introduction, Section 1 and
Section 5.2](https://math.umd.edu/~wmg/RemarkableFamily.pdf) identifies
this member of the relative character-variety family and gives its
Vieta involutions. These specific passages were inspected. In the
present sign convention they are $x\mapsto-x-yz$ and its permutations.
They preserve $f$ and $\mathfrak m$, hence act on every truncated
quotient. This established symmetry is a possible proof tool, not
evidence that (13) is known or true. The cited passages concern
surface geometry and automorphisms, not the count of initials of
these finite truncations; no whole-literature novelty assessment is
inferred from this comparison.

The checker compares full homogeneous initial row spaces with (4)
over the rationals and several prime fields, including two.
The explicit model (7) is checked only where $2$ is invertible.
Commuting cyclic operators, all free-basis vectors, original
equations, strict filtration and independent Macaulay kernels
verify the complete monomial-rank transfer. Existing complete
fan witnesses are also replayed; the additional complete finite
certificates are described above.

The additional integral model is compared against forty full Macaulay
kernels, covering both $f_0$ and $f_1$ over the rationals and four prime
fields, and against every original monomial column of the earlier
six-channel model at the stated cutoffs. Separately, 42 entire ordered
sectors and 1,069 original monic integer border relations certify the
eight all-field fans above. The earlier checks are retained unchanged;
neither selected finite fields nor the integral rank basis alone is
used as proof of characteristic-independent fan counts.

~~~bash
python3 -B scripts/run_hypersurface_truncation_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The unrestricted $G_3$ endpoint, further consequence, strongest-current
audit and honest-paper requirements are unchanged. This is a structural
reduction for the next count, not completion of that goal.
