# Rank structure and high-order perturbations in rectangle embeddings

Consider the preserved family

$$
I_F=(x^r,y^s,z-F(x,y)),\qquad
A=K[x,y]/(x^r,y^s),\qquad F=\alpha x+\beta y+H,
\quad \alpha\beta\ne0,
$$

where $H\in(x,y)^2$. Its length is $n=rs$. The previous nested-family
bound does not apply to these monomial residue matroids in general.
This note identifies which rank data are independent of $H$, then gives
an exact decomposition of the ranks when $H$ has high order. The latter
gives a uniform bound on the amount of fan complexity such a perturbation
can add. The final lattice estimate uses classical geometry; its exact
priority is unassessed. The unrestricted extremal endpoint remains
unresolved, and no new paper follows from these supporting results alone.

## All mixed multiplication-image ranks are already fixed

Write $f(y)=F(0,y)$ and define an endomorphism of $A$ by

$$
\phi(x)=F(x,y)-f(y)=xU(x,y),\qquad
\phi(y)=f(y)=yV(y).                                         \tag{1}
$$

The constant terms of $U,V$ are $\alpha,\beta$. Thus the two images
satisfy the rectangular nilpotence relations and induce an invertible
diagonal map on the maximal ideal modulo its square. Equivalently, on
the rectangular reference monomials, $\phi$ is triangular by degree
with diagonal $\alpha^i\beta^j$. It is an automorphism over every
field, preserves every coordinate ideal $(x^a y^b)A$, and sends $x+y$
to $F$.

For any finite set $E$ of triples of nonnegative integers, it follows that

$$
\dim_K\sum_{(a,b,c)\in E}x^a y^b F^c A
=\dim_K\sum_{(a,b,c)\in E}x^a y^b(x+y)^c A.                 \tag{2}
$$

Indeed each generator is changed by the unit $U^aV^b$, which does
not change the ideal it generates. In particular every rank of a mixed
multiplication operator $m_x^a m_y^b m_F^c$ is independent of $H$.
The same holds for the image rank of the map given by any list of these
operators. This is elementary coordinate algebra, not a new Lefschetz
or Jordan-type theorem.

In characteristic zero there is a familiar closed formula. Put
$R=r-a$, $T=s-b$; the rank is zero if either is nonpositive. Otherwise,
using the classical Clebsch--Gordan Jordan blocks on the smaller rectangle,

$$
\operatorname{rank}(m_{x^a y^bF^c})
=\sum_{j=0}^{\min(R,T)-1}\max(R+T-1-2j-c,0).                \tag{3}
$$

The Jordan type is $(R+T-1,R+T-3,\ldots,|R-T|+1)$.
The automorphism argument is valid in all characteristics; (3) is used
here only in characteristic zero. Small-characteristic Jordan blocks
require the established modular classification.

These image ranks are different from the ranks of the individual
monomial columns needed for a global initial ideal. Already for $r=s=2$,
the four columns $1,x,y,z$ have rank three for $F=x+y$ and rank four
for $F=x+y+xy$. In the latter case they are the unique minimum-weight
basis at $w=(1,1,1)$, since all other nonzero monomials have degree at
least two. The two ideals have five and six initial ideals over $\mathbb Q$,
respectively, despite agreeing in every image dimension in (2).
Consequently even all mixed multiplication-image ranks together cannot
determine this fan. A rank description for the desired lower construction
must retain the actual selected column spans.

## An exact decomposition for a high-order coordinate perturbation

The next statement works for any finite standard-graded plane algebra,
over any field. Let

$$
A=K[x,y]/J=\bigoplus_{j=0}^{S}A_j,
$$

where $J$ is homogeneous and $S$ is the last nonzero degree. Take
$F=\ell+H$, with $\ell=\alpha x+\beta y$, $\alpha\beta\ne0$,
and $H\in A_{\ge D}$.
Assume $2\le D\le S$ and $2D>S+1$. Set

$$
t=S-D+1<D,\qquad
U=\left(\bigoplus_{j=0}^{t}A_j\right)
   \oplus\left(\bigoplus_{j=D}^{S}A_j\right),\qquad
R=\dim_K U.                                                \tag{4}
$$

For the graph quotient $I_F=(J,z-F)$, a monomial $x^a y^b z^c$
has ordinary degree $k=a+b+c$. Expansion gives

$$
x^a y^bF^c=x^a y^b\ell^c+\text{terms of degrees at least }k+D-1.
                                                                  \tag{5}
$$

For $c=0$ the correction is zero. If $k>t$, the correction vanishes
in $A$. If $k\le t$, its degree is at least $D$, so the full column
belongs to $U$. Therefore every monomial column belongs to exactly one
of the following designated summands:

- the unchanged homogeneous column in $A_k$, when $t<k<D$;
- the combined subspace $U$, for all other degrees.

Zero columns cause no exception. Each summand is spanned by its pure
$x,y$ monomials. For every finite selected monomial set $E$, its rank
is exactly the sum of the ranks of its designated components:

$$
\operatorname{rank}_F(E)
=\operatorname{rank}_U(E_{\rm low}\cup E_{\rm high})
 +\sum_{t<j<D}\operatorname{rank}_{A_j}(E_j).                \tag{6}
$$

This is an exact matroid direct sum, without a generic-coefficient
assumption. All dependence on $H$ is confined to an $R$-dimensional
problem. The unchanged middle components are collectively controlled
by the homogeneous graph $(J,z-\ell)$.

## The state polytopes give a uniform comparison

Let $P_U$ be the convex hull of the exponent sums of bases for the
combined component, and let $P_j$ be the corresponding polytope for
$A_j$. Equation (6) gives the exact Minkowski decomposition

$$
P_F=P_U+\sum_{t<j<D}P_j.                                   \tag{7}
$$

The middle sum has a positive normal fan coarser than the whole
homogeneous graph fan: its omitted graded components are other
Minkowski summands. The graph of a linear form over a homogeneous plane
ideal has at most six initial ideals. To see this, the largest of $x,y,z$
is eliminated by the linear graph equation. In each of the three cases
the remaining ideal is homogeneous in two variables, whose initial
ideal depends only on the order of those two variables. For a rectangular
base, eliminating $z$ gives one fixed monomial ideal, improving six to
five. This is the elementary homogeneous control already used in the
earlier note.

Every nonzero monomial column has degree at most $S$, so every coordinate
of a basis profile in $P_U$ is between zero and $RS$. The classical
three-dimensional lattice-polytope vertex bound, after adding the unit
cube if needed, gives

$$
f_0(P_U)=O\bigl((1+RS)^{3/2}\bigr).                         \tag{8}
$$

The unit-cube sum preserves at least the original number of vertices
and has volume $O((1+RS)^3)$. Its normal fan refines that of $P_U$.
Intersecting each positive cone of $P_U$ with the at most six middle
cones gives at most six nonempty convex intersections. Hence

$$
g(I_F)=O\bigl((1+RS)^{3/2}\bigr),                          \tag{9}
$$

with an absolute constant, uniformly in the field, coefficients and
homogeneous plane algebra. More precisely, the factor from the unchanged
part is at most $g(J,z-\ell)$, so the argument also supplies a comparison
when that homogeneous graph is known more exactly. The normal-fan and
lattice-vertex steps are established geometry; the decomposition (6)
is the specific structural fact being used here.

## A quantitative restriction on square-rectangle constructions

Return to $A=K[x,y]/(x^q,y^q)$, with $S=2q-2$. If
$H\in(x,y)^D$ and $q\le D\le2q-2$, put

$$
h=2q-D,\qquad 2\le h\le q.
$$

The low degrees $0,\ldots,h-1$ contribute $h(h+1)/2$ dimensions;
the top degrees $2q-h,\ldots,2q-2$ contribute $h(h-1)/2$. Thus

$$
R=h^2,\qquad
g(I_F)=O\bigl(1+q^{3/2}h^3\bigr).                         \tag{10}
$$

The unaffected middle has dimension $q^2-h^2$ and at most five positive
fan cones. This gives a precise compression in actual rank, not an
assumption that the monomial ranks are determined by a Jordan partition.
At $D=q$, (10) recovers the existing $O(n^{9/4})$ square-family
window. For fixed $h$ it gives $O(n^{3/4})$ for the entire fan, with
$n=q^2$. If the first nonlinear degree $D<q$, the old all-coefficient
bound gives the unified statement

$$
g(I_F)=O\left(1+q^{3/2}\min(q,2q-D)^3\right).              \tag{11}
$$

In particular, a family with $g(I_F)\ge c n^{2+\epsilon}$ needs

$$
2q-D=\Omega\bigl(n^{5/12+\epsilon/3}\bigr).                 \tag{12}
$$

For the only possible range $0<\epsilon\le1/4$, this is automatic
when $D<q$; for high-order perturbations it follows from (10).
Perturbations confined to a thinner terminal band cannot provide the
desired power separation. General low-order nonlinear terms remain
outside this compression, and (11) does not improve the unrestricted
$G_3$ exponent or close the whole rectangle route.

## What this changes in the next proof task

Lefschetz and Jordan data in (2)--(3) can be computed exactly, but are
constant throughout the family with nonzero linear coefficients. They
cannot explain a larger fan. The high-order reduction gives exact ranks
and a smaller coefficient-dependent subproblem only in its stated
degree range. It rules out shallow terminal perturbations as an amplifier
and gives a necessary scale for any proposed construction of that type.

A new lower theorem must control the selected monomial spans when
nonlinear terms reach farther into the grading. The earlier generic
coefficient probe and its rational-surface specialization remain distinct;
the known support-matching counterexample still applies. Neither a
larger coefficient screen, an invariant mixed-operator table, nor a
fit to the old four fan counts supplies that theorem. The general
all-fields shadow route remains available in parallel.

## Literature boundary and checks

The characteristic-zero Jordan formula is classical. The inspected
[Iarrobino--Macias Marques--McDaniel preprint, Section 3.2, Lemma 3.7](https://arxiv.org/html/1802.07383v3)
states the Clebsch--Gordan blocks and connects them to the Hilbert function.
We use only its characteristic-zero part here. Its local/graded distinction
and the Jordan-type definitions in Sections 2.1 and 2.4 were also inspected.
The elementary conjugacy (1) supplies the additional all-characteristic
comparison for these particular rectangle embeddings.

[Nicklasson, Sections 1--3, especially Theorems 2.2 and 3.2--3.4](https://link.springer.com/article/10.1007/s13348-017-0209-3)
records the established characteristic-zero input and classifies the
two-variable monomial complete intersections with strong Lefschetz in
positive characteristic using base-$p$ digits. Those statements were
inspected. They concern multiplication ranks, not the full selected-column
matroid or the number of global initial ideals. No new prime-digit or
Lefschetz classification is claimed.

The state-polytope and lattice-vertex input in (8)--(9) is the same
established argument used in
[Babson--Onn--Thomas, Theorem 2.4 and Corollary 2.5](https://arxiv.org/pdf/math/0207135)
and in the existing local Hilbert-function bound. The new parameter is
the exact dimension of the combined low and high degree pieces, rather
than the entire quotient length. Priority of this comparison theorem
remains unestablished; an unsuccessful search does not establish novelty.

Run from the repository root:

```bash
python3 -B scripts/run_rectangle_rank_checks.py
./build.sh
python3 tests/test_smoke.py
```

The report includes exact automorphism matrices, mixed image-rank checks,
a complete small fan comparison with identical operator-rank data,
and checks of the component decomposition on every column of bounded
full monomial grounds. The controls include the nonmonomial homogeneous
quotient $(x^3-y^3,x^2y,xy^2)$ in characteristics zero, two and three;
its two pure cubic residues coincide. Full-matrix and projected greedy selections give
independent controls of the reduced problem. No larger fan-count fit is
used. The core kernels, certification range, prior manuscripts and paused
ML work remain unchanged.
