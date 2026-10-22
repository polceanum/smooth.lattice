# Marked corner cuts and the size of generic reduced-basis output

This note proves a weighted consequence of the exponent-prefix counting
construction. For each fixed dimension $d\ge2$, the total number of
minimal excluded monomials, summed over the distinct $n$-point corner
cuts, has order $n^{d-1/d}$. The logarithms in the unweighted cut count
disappear. For a generic radical point ideal this is exactly the number
of literal monic polynomials in the union of its reduced Gröbner bases.

This is a generic output theorem, not a bound on the unrestricted
extremal function $G_3(n)$ or a family separating it from the generic
scale. Its exact priority is unassessed. The connection with point
ideals and their determinant relations is classical; the weighted
estimate is proved below rather than inferred from those formulas.
The active extremal goal remains incomplete.

## Parameters and statement

Let $\mathcal C_d(n)$ consist of the subsets of $\mathbb N_0^d$ of size
$n$ strictly separated from their complements by a positive linear
score. For $D\in\mathcal C_d(n)$ let $E(D)$ be the coordinatewise
minimal elements of $\mathbb N_0^d\setminus D$. Put

$$
 R_d(n)=\sum_{D\in\mathcal C_d(n)}|E(D)|.                 \tag{1}
$$

Coordinates are labelled. An incidence $(D,a)$ with $a\in E(D)$ is
counted once, regardless of how many weights or internal orders realize
$D$. This differs from counting all successor points of $D$.

**Theorem 1.** For every fixed $d\ge2$,

$$
 R_d(n)=\Theta_d(n^{d-1/d}).                             \tag{2}
$$

Let $F$ be any field and $K=F(t_{ij}:1\le i\le n,1\le j\le d)$, with
all $t_{ij}$ algebraically independent. Let $I$ be the ideal of the $n$
points $(t_{i1},\ldots,t_{id})$ in $K^d$, and let $\mathcal U(I)$ be
the union of the **monic reduced** bases for all global monomial orders.
Then

$$
 |\mathcal U(I)|=R_d(n),\qquad
 \sum_{f\in\mathcal U(I)}|\operatorname{supp}f|
       =(n+1)R_d(n)=\Theta_d(n^{d+1-1/d}).                \tag{3}
$$

The same term count holds if every reduced basis is listed separately.
Thus in three variables there are $\Theta(n^{8/3})$ literal monic
polynomials and $\Theta(n^{11/3})$ nonzero term records. These claims
hold in every characteristic over the displayed field extension. They
do not assert that any fixed finite base field has sufficiently many
generic points without extending the field.

The output convention matters. Different monic polynomials in (3) can
be scalar multiples of the same circuit. Let $\overline{\mathcal U}(I)$
identify nonzero scalar multiples. Then

$$
 \frac{c_d n^{d-1/d}}{(1+\log(n+1))^{d-1}}
 \le |\overline{\mathcal U}(I)|\le C_d n^{d-1/d}.          \tag{4}
$$

Equation (4) preserves the polynomial exponent without relying on a
choice of normalization. No claim about the minimum possible size of
an arbitrary, possibly nonreduced, universal Gröbner basis follows.

## Upper bound: retain the aspect ratio in the border count

The proof refines the active-face dyadic argument in the
[ordered-prefix theorem](../paper/corner_cut_growth/main.tex), Section
“Counting complete ordered prefixes.” Choose a strictly separating
weight, perturbed to give distinct scores on $D$, and put
$t=\max_{a\in D}w\cdot a$. Assume $n>1$. The active set is
$J=\{i:w_i\le t\}$, with $r=|J|\ge1$. All other coordinates vanish
on $D$. In the active coordinates write

$$
 b_i=t/w_i\ge1,\qquad L_i=2^{\lfloor\log_2b_i\rfloor},
 \quad Q=\prod_{i\in J}L_i=2^R.
$$

The lattice cubes covering the simplex of score at most $t$, and its
enclosing coordinate box, give the shell bounds

$$
 n<4^rQ,\qquad Q\le r!n.                                \tag{5}
$$

For fixed $r$, only a bounded number of integers $R$ occur. In a fixed
chart $(J,L)$, all included points lie in
$\prod_{i\in J}\{0,\ldots,2L_i-1\}$. Comparisons on this box involve
at most $4^rQ$ distinct central hyperplanes. Their arrangement in
$\mathbb R^r$ has $O_r(Q^{r-1})$ chambers. Each chamber fixes the
complete order on the box and hence at most one $n$-point cut with a
witness in this chart. Overcounting cuts between charts is harmless.

The excluded heads supported in $J$ form an antichain in
$\prod_{i\in J}\{0,\ldots,2L_i\}$. Deleting any one coordinate
projects that antichain injectively: points with the same projection
would be comparable. Consequently

$$
 |E_J(D)|\le\min_{i\in J}\prod_{j\in J\setminus\{i\}}(2L_j+1)
       \le 3^{r-1}Q/L_{\max}.                            \tag{6}
$$

The remaining heads are precisely the $d-r$ inactive unit vectors.
Indeed, any other excluded point using an inactive coordinate dominates
one of those units and so is not minimal.

For every fixed $r$ the needed dyadic identity is

$$
 \sum_{\substack{j_1+\cdots+j_r=R\\j_i\ge0}}
             2^{-\max_i j_i}=O_r(2^{-R/r}).              \tag{7}
$$

Here is a uniform proof. Write $M=\lceil R/r\rceil+s$. The number of
tuples with maximum $M$ is at most
$\binom{rM-R+r-1}{r-1}$, by putting $u_i=M-j_i\ge0$ and discarding
the restrictions $u_i\le M$. Since $0\le rM-R<r(s+1)$, this is at
most $r^{r-1}(s+2)^{r-1}/(r-1)!$. The sum over $s\ge0$ against
$2^{-s}$ converges. For example, using the rising factorial gives
the valid constant $r^{r-1}(2^{r+1}-2)$ times
$2^{-\lceil R/r\rceil}$, which proves (7).

Multiplying the per-chart cut count by (6), then summing (7), bounds
the active heads at one shell by

$$
 O_r\left(Q^r\sum_{j_1+\cdots+j_r=R}2^{-\max j_i}\right)
       =O_r(n^{r-1/r}).                                  \tag{8}
$$

Inactive units contribute at most
$O_d(n^{r-1}(1+\log(n+1))^{r-1})$. This term is absent when $r=d$;
when $r<d$ it is $O_d(n^{d-1/d})$. Summing the finitely many active
sets and shell values proves the upper bound in (2). For $n=1$,
$R_d(1)=d$ directly.

## Lower bound at every sufficiently large prescribed length

We use one balanced aspect range of the primitive-frame construction,
rather than the logarithmically many ranges used to count all cuts.
The required uniform input is explicit in the
[existing proof](../paper/corner_cut_growth/main.tex), Lemmas
“Primitive frames,” “Fixed-size injectivity,” and the extent estimate:
if all axis lengths are $\ell$, $Q=\ell^d$, $\ell$ is sufficiently
large, and $n/Q\ge T_d$, there are at least $c_d Q^{d-1}$ distinct
$n$-point cuts. Their integer normals satisfy

$$
 (9/10)H/\ell<a_i<(11/10)H/\ell.                        \tag{9}
$$

Here $T_d$ and $c_d>0$ are the fixed constants in that construction;
the hypotheses do not require a subsequence of $n$. To see the
consequence precisely, set $\ell=\lfloor(n/T_d)^{1/d}\rfloor$.
For sufficiently large $n$, $Q\asymp_d n$, $n/Q\ge T_d$, and all
primitive-frame size and aspect hypotheses hold. The fixed-size
injectivity argument distinguishes $\Omega_d(n^{d-1})$ cuts.

The construction orders score ties lexicographically. Its explicit
strict integer witness is
$w=n^d a+(n^{d-1},n^{d-2},\ldots,1)$. By (9), for sufficiently
large $n$ its coordinate ratio is less than two. A further small
generic perturbation, if needed, preserves the cut and this ratio.
Thus at least $c'_d n^{d-1}$ distinct cuts have a strict witness which
can be normalized to

$$
 1\le w_i\le2,\qquad\min_iw_i=1.                       \tag{10}
$$

Every cut with such a witness has many minimal excluded heads. Put
$t=\max_D w\cdot a$. Since $D\subset[0,t]^d$, one has
$t\ge n^{1/d}-1$. Relabel a minimum-weight coordinate as $d$.
Choose each of $a_1,\ldots,a_{d-1}$ independently in
$\{1,\ldots,\lfloor t/(4d)\rfloor\}$. Their score $s$ is less
than $t/2$. Set

$$
 a_d=\lfloor t-s\rfloor+1.
$$

The resulting point has score in $(t,t+1]$ and positive coordinates.
Subtracting any coordinate unit decreases score by at least one, so
every immediate predecessor has score at most $t$ and belongs to
$D$. The point is therefore in $E(D)$. Distinct choices give distinct
heads. For large $n$ this constructs
$\Omega_d(t^{d-1})=\Omega_d(n^{1-1/d})$ heads for each balanced
cut. Multiplication proves the lower bound in (2). Finitely many
small positive lengths can be absorbed into the constants.

This is a weighted consequence of the existing uniform frame proof,
not an independent proof or priority claim for that construction.

## Generic point relations and normalization

For an exponent $a$ let $v_a=(\prod_jt_{ij}^{a_j})_{i=1}^n$ be its
evaluation column. Any $n$ distinct such columns are linearly
independent over $K$. In their determinant, different permutations
assign different exponent vectors to some row of independent variables.
The determinant monomials are therefore distinct with coefficients
$1$ or $-1$ and cannot cancel, including in characteristic two.
This is the determinant argument of Onn and Sturmfels, Lemma 5.1.

Greedy interpolation for any global monomial order consequently
accepts its first $n$ monomials. The finite-weight representation
identifies their exponents with a member $D$ of $\mathcal C_d(n)$.
Conversely each such cut is realized. For each $a\in E(D)$ its
unique reduced relation is

$$
 f_{D,a}=x^a-\sum_{b\in D}\lambda_bx^b.                 \tag{11}
$$

Every $\lambda_b$ is nonzero: otherwise at most $n$ distinct columns
would be dependent. Hence its support is exactly $S=D\cup\{a\}$.
These are the classical reduced-basis formulas in Onn and Sturmfels,
Theorem 5.2 and Corollary 5.3; they are not new here.

It remains to check that counting incidences does not accidentally
count identical monic polynomials. Different supports distinguish
polynomials immediately. On a fixed ordered support $S$ of size
$n+1$, let $c_s$ be the signed determinant of the columns other than
$s$. A common circuit is $f_S=\sum_{s\in S}c_sx^s$, and the monic
polynomial with head $a$ is $f_S/c_a$. Give $t_{ij}$ multidegree
$e_j$, simultaneously for all rows $i$. Then

$$
 \deg c_s=\sum_{u\in S}u-s.
$$

The nonzero polynomials $c_s$ have different multidegrees for different
$s$. They cannot be equal. Thus $f_S/c_a\ne f_S/c_b$ if $a\ne b$,
in every characteristic. This proves the equality and term count in
(3), not just an upper bound obtained by adding basis sizes.

## Identifying nonzero scalar multiples

For a fixed support $S$, every possible deleted head $a$ must be a
vertex of $\operatorname{conv}S$: the strict separator makes it
uniquely largest among $S$. Fix one realization $D=S\setminus\{a\}$.
Then

$$
 \#\operatorname{vert}(\operatorname{conv}S)
 \le1+\#\operatorname{vert}(\operatorname{conv}D).        \tag{12}
$$

The linear integer-score certificate in the
[prefix proof](../paper/corner_cut_growth/main.tex), Lemma “Linear
integer-score certificate,” writes $D$ as exactly the integer points
in a simplex

$$
 P=\{x\ge0:A\cdot x\le S_0\},\quad
 A_i\in\mathbb Z_{>0},\quad S_0=O_d(n),\quad A_i=O_d(n).
$$

This is $d+1$ inequalities in fixed dimension with encoding size
$O_d(1+\log(n+1))$. The classical integer-hull bound of
[Cook, Hartmann, Kannan and McDiarmid, Theorem 2.1, p. 29](https://www.math.uwaterloo.ca/~bico/papers/ippoints.pdf)
therefore gives
$O_d((1+\log(n+1))^{d-1})$ vertices of $\operatorname{conv}D$.
No full-dimensionality hypothesis is needed. Equations (12) and (2)
give (4). Distinct supports are also distinct scalar classes since
every coefficient of (11) is nonzero.

## Relation to prior work and the extremal objective

[Onn and Sturmfels, *Cutting Corners* (1999)](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Onn/Onn4.pdf),
Section 5, already identify generic point ideals with corner cuts,
write their reduced bases by determinants, and obtain a polynomial
universal-basis bound in Corollary 5.4. That corollary multiplies a
general state-vertex bound by the largest staircase generator count;
it does not give the sharp weighted estimate (2). Their genericity
definition suffices for the state dictionary; the purely transcendental
configuration here also guarantees dense tails and distinct monic
normalizations. No new interpolation or state-polyhedron dictionary
is asserted.

[Wagner's thesis (2003)](https://ti.inf.ethz.ch/ew/research/phdtheses/wagner03.pdf),
Section 3.1, reviews that dictionary and universal-basis application.
Its cut-count upper bound, Theorem 3.15, controls unweighted states.
The balanced primitive-frame lower bound and ordered-prefix chart
bound used here are from the saved corner-cut manuscript; their own
priority assessment is not settled by this note. The integer-hull
input in (4) is an established theorem. These exact comparisons do
not rule out the weighted law appearing elsewhere. A complete
priority assessment of that law remains open.

[Elizalde and Galván, inspected preprint v1](https://arxiv.org/html/2312.16353v1),
Definition 2.5, Lemma 2.6 and Section 6, use “addable” for a cell
whose addition **preserves triangularity**. Such cells number one or
two, whereas $E(D)$ here contains every addition preserving the
downset property. Its size can grow. Their generating functions for
triangular-preserving additions therefore count a different marked
parameter. The checker explicitly distinguishes these sets. As in the
[standing version note](interleaving_research.md), relevant results
are assumed unchanged in the 2025 journal version, pending the
deferred version comparison; this is not a novelty argument.

One further comparison candidate is Maruri-Aguilar, Sáenz-de-Cabezón
and Wynn, [*Betti numbers of polynomial hierarchical models for
experimental designs* (2012)](https://doi.org/10.1007/s10472-012-9295-9).
The author's publication page and institutional search record identify
it as a study of extremal Betti numbers and corner-cut models. Its
full text was not inspected in this pass: the indexed institutional
PDF URL returned HTTP 404. An abstract is insufficient to exclude an
aggregate generator-count result there. This remains a specific
priority check, not evidence either for novelty or for equivalence.

The same cuts describe the finite smooth-number exponent prefixes.
Equations (2)--(3) quantify all minimal interpolation relations across
their generic algebraic models. They do not give a theorem about
prime distribution. Nongeneric algebraic dependencies can admit
standard sets which are not ordinary corner cuts, so this argument
does not bound their total border incidence by (2). In particular it
does not close the gap between the generic three-variable fan scale
$\Theta(n^2\log^2(n+1))$ and the current unrestricted upper bound.

Combining (2) with the saved sharp cut count gives a mean reduced-basis
size of order $n^{1-1/d}/\log^{d-1}n$ over distinct generic initial
sets. This is not an average over weights and is not a bound for
every individual shape. The explicit term count is a lower bound for
writing the specified output with separate coefficient/monomial
records. It does not measure coefficient bit lengths, compressed
representations, finding one basis, or all possible universal-basis
algorithms. The polynomial count up to scalars is (4), not the exact
monic count (3).

## Reproducible independent controls

Run from the repository root:

```bash
python3 -B scripts/run_marked_border_checks.py
./build.sh
python3 tests/test_smoke.py
```

The [exact report](../results/research/marked_border_checks.json) is
generated by the first command. It enumerates every downset through
length 24 in dimension two and length 10 in dimension three. An
independent full-box antichain extraction is checked against the
existing immediate-predecessor border routine. Exact rational slope
intervals or polygon clipping test strict separation; their complete
weight coverage is checked at each length. Every retained cut checks
the dyadic shell and all antichain projections. The balanced head
construction is replayed whenever its witness has coordinate ratio
at most two. Every planar cut through length 23 also checks its one
or two triangular-preserving additions against the next complete
level, and the report records a strict difference from the ordinary
minimal border. Inclusion-exclusion shell sums are checked against
independent tuple enumeration on small shells, and against the
proved uniform constant for $1\le r\le8$, $0\le R\le120$.
Six larger balanced cuts in dimensions two through four additionally
test the lower-head construction with at least two choices in each
free coordinate. Their borders are checked by independent score
inequalities, so these controls do not rely on an empty construction
in the small-length enumeration.

The dense-output controls use actual radical point quotients over
$\mathbb F_p(T)$ with points
$p_{i,j}=(iT)^{(n+1)^j}$, for $1\le i\le n$ and $0\le j<d$.
Every required cofactor is explicitly tested nonzero; genericity is
not assumed from a large prime. Integer Bareiss determinants and a
separate modular interpolation solver independently check all
coefficients. Exact powers of $T$ are retained, so counting the monic
polynomials does not conflate numeric specializations with rational
functions. The cases cover $n\le8$ in dimension two, $n\le6$ in
dimension three, and both reported characteristics. They verify the
incidence, scalar-class and term-record identities exactly. Finite
controls are implementation evidence for the proof mechanisms, not
an asymptotic proof, novelty assessment or goal completion.
