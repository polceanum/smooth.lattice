# Exact ranks for truncations of two-part hypersurfaces

The shared depth bound leaves room for a polynomially larger fan when
relations propagate through a growing number of layers. This note
makes one such class explicit. It gives an exact homogeneous initial
ideal and a finite free-module description of **all** monomial ranks.
It does not count the entire positive fan or establish novelty. The
free-module and graded-initial mechanisms are classical; the main
extremal goal and its publication gates remain open.

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
Closing it requires a whole-triangle count of original monomial
rank changes. The exact kernel (7), including its coupled
coefficients, is available for that proof. Increasing the small
enumeration or treating the six channels as independent prices
does not close the window.

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

The checker compares full homogeneous initial row spaces with (4)
over the rationals and several prime fields, including two.
The explicit model (7) is checked only where $2$ is invertible.
Commuting cyclic operators, all free-basis vectors, original
equations, strict filtration and independent Macaulay kernels
verify the complete monomial-rank transfer. Existing complete
fan witnesses are replayed; no larger fan search is performed.

~~~bash
python3 -B scripts/run_hypersurface_truncation_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The unrestricted $G_3$ endpoint, further consequence, strongest-current
audit and honest-paper requirements are unchanged. This is a structural
reduction for the next count, not completion of that goal.
