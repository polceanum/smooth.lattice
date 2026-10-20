# Periodic covers, coherent arithmetic and a failed amplification route

The remaining mixed-period candidate cannot be treated as a regular,
characteristic-independent ordering problem. This note separates its
periodic lattice cover from its circuit quotient, records an exact
positively coherent basis with index two, and explains why amplifying
that index does not amplify the number of initial ideals. The underlying
duality, good-prime principle and flat substitution are established tools.
None of these results settles the unrestricted extremal goal or establishes
priority of a new counting theorem.

## The full periodic cover has a quadratic bound

Let $r,s,t\ge2$ be pairwise coprime. Put $M=rst$, and define

$$
B=K[U,V,W]/(U^r-1,V^s-1,W^t-1),\qquad
X=UV,\ Y=UW,\ Z=VW.
$$

These generators produce all of $B$ by the same CRT recovery used for
the [circuit quotient](cyclotomic_price_regions.md). A monomial has label
$(a+b\bmod r,a+c\bmod s,b+c\bmod t)$. The label lattice $L$ has
index $M$ in $\mathbb Z^3$. Over every field, including characteristics
dividing a period, a monomial basis contains exactly one representative
of every label. Every difference of basis profiles belongs to $L$.

The exact orders of $X,Y,Z$ are $rs,rt,st$, so all initial standard
sets lie in the box with these side lengths. Their profile polytope lies
in a box of side lengths $M rs,M rt,M st$, of volume $M^5$.
In its affine lattice coset the normalized volume is at most $M^4$.
The classical three-dimensional lattice-vertex theorem therefore gives

$$
g(B)=O(M^2).                                             \tag{Q1}
$$

For a lower-dimensional profile polytope, add a fundamental
parallelepiped of $L$. A triangular reduced lattice basis has coordinate
entries bounded by the corresponding pure-axis orders, since their
axis multiples lie in $L$. The added widths are at most three times
those orders, so the same normalized-volume order holds. Its normal
fan refines the original one. This supplies the regularization without
assuming that a unit cube is a fundamental cell of $L$.

More generally the same argument gives $O(m\sqrt{abc})$ for a
three-generator finite group algebra of order $m$ whose generator orders
are $a,b,c$. Formula (Q1) is a relative-volume corollary, not a new
lattice-vertex result. The lattice-ideal and shortest-coset dictionary
is established in [Babson--Onn--Thomas, Section 5](https://arxiv.org/pdf/math/0207135).
The exact priority of this numerical specialization has not been audited.

The loss of the profile index in $A$ is exact, not merely possible.
For coprime $r,s\ge2$, the minimal polynomial of $X=UV$ in its
coefficient algebra is

$$
F_{r,s}(T)=\frac{(1-T^{rs})(1-T)}{(1-T^r)(1-T^s)},
\qquad \deg F_{r,s}=(r-1)(s-1).                          \tag{Q1a}
$$

Over characteristic zero the $rs$ distinct labels for $UV$ are
the $rs$th roots, and deleting the $U=1$ or $V=1$ labels gives this
polynomial. Its monic integral division has no remainder: equivalently,
use the disjoint factorization of $T^{rs}-1$ into cyclotomic factors
and cancel those dividing $r$ or $s$. The identity survives in the
free integral tensor quotient and hence in every characteristic.
CRT recovery makes $K[X]=K[U,V]$ of the displayed dimension in
every field, so the polynomial remains minimal, including repeated-root
cases. Its constant and linear coefficients are $1,-1$.

If a group grading makes $A$ homogeneous, its pure-$X$ elimination
ideal is homogeneous too. The monic minimal polynomial must have a
single group degree: a different component would be an annihilating
polynomial of smaller degree. Its nonzero constant and linear terms
force $\deg X=0$. Applying the same argument to $Y,Z$ shows that
every variable-induced group grading is trivial. Thus the
[intrinsic profile-difference lattice](profile_lattice_bound.md) is
exactly $\mathbb Z^3$ over every field.
The cover's index-$M$ volume saving cannot be carried into $A$ by
assuming a residual finite grading. Formula (Q1a) is a classical
inclusion-exclusion polynomial; this is an explicit obstruction in
the present construction, not a new cyclotomic-polynomial claim.

## The dual quotient and the cost obstruction

Let $S_d(T)=1+\cdots+T^{d-1}$, and recover the original algebra as

$$
A=B/J,\quad J=(S_r(U),S_s(V),S_t(W)),\quad
n=(r-1)(s-1)(t-1).
$$

The coefficient-of-$1$ pairing $\langle f,g\rangle=[1](fg)$ on $B$
is perfect: group monomials pair with their inverse labels. This uses
no division by $M$. In one factor the annihilator of $S_r(U)$ is
$(U-1)$; multiplication by $S_r$ sends $h$ to $h(1)S_r$. Tensoring
the factors shows

$$
\operatorname{ann}_B(J)=((U-1)(V-1)(W-1)).
$$

Consequently the length-$k=M-n$ algebra

$$
C=B/((U-1)(V-1)(W-1))                                  \tag{Q2}
$$

represents the dual column matroid when its columns are indexed by
**inverse monomials**. This is an explicit algebra form of the known
cyclotomic/simplicial duality, not a new matroid identification; compare
[Martin--Reiner, Theorem 1 and its proof](https://www-users.cse.umn.edu/~reiner/Papers/cyclotomic.pdf).

A basis of $C$ consists of the periodic monomials with at least one
zero coordinate. For positive $i,j,k$, reduction uses

$$
U^iV^jW^k=U^iV^j+U^iW^k+V^jW^k-U^i-V^j-W^k+1.
$$

This also gives a direct field-independent dual representation. In
the $X,Y,Z$ coordinates the fixed polynomial

$$
(Z-XY)(Y-XZ)(X-YZ)                                      \tag{Q3}
$$

vanishes in $C$: its image is a unit times
$(U^2-1)(V^2-1)(W^2-1)$, a multiple of the defining product in
(Q2). No claim that the factors $U+1,V+1,W+1$ are units is needed
in characteristic two.

For a positive generic weight, let $E$ be the standard set of $B$.
The standard set of $A$ is contained in $E$. By ordinary matroid
duality its complement in $E$ is a **maximum**-cost basis of the
inverse-monomial representation of $C$. This is not the positive
minimum-basis problem for all monomials of $C$.

Reflection does not fix the distinction. Set $\rho_i=\max_{a\in E}a_i$.
Multiplication by the common unit of label $\rho$ turns inverse columns
into columns labelled $\rho-a$, and changes maximum costs into minimum
costs on the finite reflected ground. Unless $E$ is a coordinate box,
$\rho\notin E$, so the reflected ground has no constant monomial.
It cannot contain a global initial standard set of $C$. Thus neither
(Q1), the low-degree carrier (Q3), nor a count of global initial ideals
of $C$ can be imported as a count for the omissions without an additional
argument about their actual finite-ground costs.

## A characteristic-sensitive coherent basis

The effect is visible in the original three-prime family. Take
$(r,s,t)=(5,7,11)$, $n=240$, $L=481$, and

$$
w=(L^3+1,L^3+L,L^3+L^2).
$$

Exact ordinary quotient arithmetic gives two strict initial standard sets:
$B_0$ over $\mathbb Q$ and $B_2$ over $\mathbb F_2$. Their only
differences are

$$
B_0\setminus B_2=\{XY^{11}\},\qquad
B_2\setminus B_0=\{Y^{12}\}.                            \tag{Q4}
$$

Relative to the integral tensor basis, $|\det B_0|=2$ and
$|\det B_2|=1$. The rational basis has rank $239$ modulo two.
Its normal form for $Y^{12}$ has denominator two. Clearing denominators
gives an integral relation with leading term $2Y^{12}$; modulo two
its leading monomial is $XY^{11}$. The complete relation and both exact
standard sets are in the [report](../results/research/periodic_quotient_checks.json).
Every minimal border inequality is checked strictly. Their ordinary-degree
histograms agree, so the change is not a change in quotient length or
in that histogram. Two divides none of $5,7,11$: each geometric-sum
polynomial is squarefree in characteristic two. Nilpotents are not
needed for this ordering change.

This is a coherent arithmetic obstruction, not merely a nonunimodular
abstract basis. It refutes the shortcut that all positively ordered bases
in this family are integral or characteristic independent. It says nothing
by itself about which characteristic has more cones, or about a growth law.

For reference, the exact good-prime test in a finite free integral algebra
is elementary. If a rational initial standard set $D$ has index
$\Delta$ in the algebra's integral lattice, it remains that same standard
set modulo $p$ **if and only if** $p\nmid\Delta$. If the determinant
is a unit modulo $p$, its border coordinates reduce and their strict
weight inequalities persist; if not, the proposed columns are dependent.
The good/unlucky-prime framework is established: [Arnold, Theorems
5.11--5.13](https://www3.risc.jku.at/research/theorema/Groebner-Bases-Bibliography/gbbib_files/publication_355.pdf)
was inspected in its homogeneous setting. The finite-free argument above
does not assume the nongraded problem is covered verbatim by that theorem.

## Arbitrarily large index does not amplify the fan

For any ideal $I\subset K[X,Y,Z]$, form
$I^{(d)}=(f(X^d,Y,Z):f\in I)$, $d\ge1$.
The ambient polynomial ring is free over $K[X^d,Y,Z]$, with
basis $1,X,\ldots,X^{d-1}$. Substituting $X^d$ into a Gröbner
basis preserves its leading monomials, reductions and $S$-polynomial
identities under the induced weight $(d w_X,w_Y,w_Z)$. Hence

$$
\operatorname{Std}_w(I^{(d)})
 =\{(d a+i,b,c):(a,b,c)\in
 \operatorname{Std}_{(d w_X,w_Y,w_Z)}(I),\ 0\le i<d\}.
$$

The positive weight map is bijective. In particular

$$
\dim K[X,Y,Z]/I^{(d)}=d\dim K[X,Y,Z]/I,\qquad
g(I^{(d)})=g(I).                                        \tag{Q5}
$$

For a finite free integral algebra, the lifted monomial-basis matrix
is a block diagonal repetition after ordering by exponent residue.
Its index is $\Delta^d$. Applying this to (Q4) gives coherent rational
bases of length $240d$, index $2^d$ and rank defect $d$ modulo two.
Yet the number of initial ideals in each characteristic is unchanged
under this substitution. This is a free-extension consequence, not a
new general flatness theorem or an infinite counting separation. These
extensions are not asserted to remain the original three-prime family.

This excludes every finite monomial self-map of affine three-space,
not just a chosen power substitution. If each variable maps to a
nonconstant monomial and the polynomial-ring map is finite, each
coordinate axis needs a pure power among the three images. Otherwise
the fiber over the origin maps onto a polynomial ring in that axis
coordinate and is not finite. Three axes and three images force a
permutation of pure powers. Conversely such maps are finite and free.
Composing (Q5) in the three coordinates therefore shows that every
such finite monomial pullback preserves the full count. This elementary
toric observation is a mechanism exclusion, not a new classification
claim. A monomial map finite only on a particular unit-valued quotient,
or a general nonmonomial finite map, is outside this argument.

Thus torsion magnitude alone and this natural repetition construction
cannot supply the requested growth advantage. The remaining task is
to control or construct more coherent **cost regions** across multiple
periods, retaining characteristic-sensitive ranks. The quadratic
one-period exclusion, all earlier research and the full unrestricted
goal remain in force.

A concrete remaining target is the number of circuit-quotient
refinements **inside one fixed lattice-cover cone**, with its actual
inverse-label maximum costs. A uniform polylogarithmic bound there,
together with (Q1), would exclude a polynomial separation for this
family. A lower bound would have to count distinct quotient standard
sets; the number of overlay cells of two fans cannot simply be
multiplied. Characteristic differences and torsion orders alone do not
provide that count.

## What can be proved inside a fixed cover cone

Fix a full-period cover standard set $E$ of size $M$, and let
$A$ have length $n=M-k$. The omitted set $O=E\setminus\operatorname{Std}(A)$
is a rank-$k$ maximum basis in the inverse-label dual representation.
It is also an upset in the monomial downset $E$.

This gives a coefficient-independent pruning rule. Define

$$
U_k(E)=\{a\in E: |\{b\in E:b\ge a\}|\le k\}.
$$

Every actual omitted set lies in $U_k(E)$, since it contains all
elements of $E$ above each of its members. Restricting the dual
optimization to this shell preserves its optimum throughout the cover
cone. If $s$ is a maximal monomial of $E$ above $a\in U_k(E)$,
then the whole box $[a,s]$ lies in that upper interval, so

$$
\prod_i(s_i-a_i+1)\le k.
$$

Writing $\sigma$ for the number of maximal monomials of $E$ and
$V_k=\{u\ge0:\prod_i(u_i+1)\le k\}$, we obtain

$$
U_k(E)\subseteq\bigcup_{s\text{ maximal}}(s-V_k),\qquad
|U_k(E)|=O(\sigma k\log^2(k+1)).                         \tag{Q6}
$$

The intersection with nonnegative exponents and $E$ is implicit on
the right. This is an elementary shifted divisor-box argument. It
does not assert that $\sigma$ is small, that every shell subset is
realizable, or that shell size counts coherent bases.

There is also a uniform trace bound. On the line
$w(\alpha)=(\alpha,\beta,\gamma)$, with $\beta,\gamma$ fixed,
the omitted profile maximizes
$\alpha\sum_{a\in O}a_X+\beta\sum a_Y+\gamma\sum a_Z$.
Its integral $X$-slope strictly increases whenever the unique generic
optimum changes. Consequently any interval of this line inside the
fixed cover cone has at most

$$
1+k\bigl(\max_{a\in E}a_X-\min_{a\in E}a_X\bigr)         \tag{Q7}
$$

distinct generic states. Equal-slope affine pieces cannot exchange
optimality on open intervals. Ties at isolated comparison events do
not create additional strict states of the line; other directions
of tie refinement are not counted here.

The existing [rational-shadow theorem](rational_shadow_bound.md) can
likewise be applied with rank $k$, not $n$, because $E\subset V_M$.
For any rank-two integer projection $T$, put
$\Delta(T)=\sum_{i<j}|\det(T_i,T_j)|/\gcd_{i<j}|\det(T_i,T_j)|$.
The projected dual basis polytope has at most

$$
O\!\left(\log(M+1)
 [1+k+(k^2M\Delta(T))^{1/3}]\right)                      \tag{Q8}
$$

vertices. Restrictions to the cover cone cannot increase the number
of uniquely optimal projected vertices. The projection height remains
in the formula: it cannot be dropped for an arbitrary rational price
line. This is an application of the earlier proved bound, not a new
unrestricted exponent, and line traces do not count a two-dimensional
fan.

There is a further pruning rule that uses the actual dual matroid,
not its regularity. Suppose $r,s,t>2$. The three integer directions

$$
v_1=(1,1,-1),\quad v_2=(1,-1,1),\quad v_3=(-1,1,1)
$$

change exactly one period label by two. For a weight avoiding
$w\cdot v_i=0$, put $d_i=\operatorname{sign}(w\cdot v_i)v_i$.
If all eight monomials

$$
a+\epsilon_1d_1+\epsilon_2d_2+\epsilon_3d_3,
\qquad \epsilon_i\in\{0,1\},                            \tag{Q9}
$$

belong to $E$, they form a circuit in the inverse-label complement
representation. Indeed the alternating sum is a product of three
coordinate differences, a multiple of $(U-1)(V-1)(W-1)$ in $B$.
In the tensor product of the three two-point coordinate spaces, its
intersection with the augmentation product is one-dimensional.
Thus the eight columns have rank seven and every seven are independent,
over every field, including characteristic two.

The base corner $a$ is uniquely cheapest among these columns, so the
classical matroid cycle property excludes it from every maximum-cost
dual basis. Each sign sector has a coefficient-independent candidate
set obtained by these eight-corner membership tests. Once $E$ is
given, its construction takes $O(M)$ set-membership operations. A uniform $O(k)$ bound for this cube filter alone is false without
additional geometry, even when all three period thresholds are reached;
an infinite family below proves this. The joint shell/cube size and
sublinear construction remain separate questions. The three sign planes give a constant number of positive
price sectors, but they do not count the quotient's cones.
This extends the earlier four-cycle cost-domination mechanism to the
known simplicial dual representation; it does not assume that all its
circuits are cubes or that the matroid is regular. The useful remaining target is the joint shell/circuit geometry,
with genuinely three-variable states distinguished from planar ones.

### Complete finite traces and their limits

The exact checker fixes the cover at
$w_0=(L^3+1,L^3+L,L^3+L^2)$, with $L=2n+1$, and varies only
$w_X$. The cover's strict binomial borders determine its entire open
line interval. Every pair-comparison event of its finite ground is
enumerated; one sample in every complementary interval determines
the whole generic trace. Each new quotient state is independently
verified by ordinary polynomial arithmetic and strict borders, and
the fixed cover is checked by another shortest-path computation.
Both the upper-shell and cube exclusions are checked against every
generic interval's optimum. Sample cube circuits are independently
verified by their alternating relations and all eight rank-seven tests.

For $(3,5,7)$ there are three quotient states in each checked field
$\mathbb Q,\mathbb F_2$. They stay within an already excluded
one-period class. For $(5,7,11)$ there are thirteen states in each
field, all reaching all three period thresholds. These thirteen
states belong to one fixed cover cone; a cover diagram therefore
does not determine the quotient diagram. The saved trace covers all
310 generic comparison intervals there. Equality of these two finite
counts is not equality of their state sets or full fans. These prices
are bounded in ratio, so the data do not supply a separating family.

Two primary comparisons guard against stronger inferences.
[Hrubeš--Yehudayoff, Theorem 4 and Lemma 5](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.CCC.2021.9)
relate line breakpoints and planar shadows, and bound greedy-like
optimization by pair comparisons. Those statements and proofs were
inspected; the trace enumeration uses that established principle.
[Sabariego--Santos, Section 4 and Theorems 5.1/5.3](https://arxiv.org/pdf/0705.3631)
give unbounded coherent minimum-distance diagrams for triple-loop
networks. Their diagrams first minimize total route length and then
resolve ties, whereas a full positive-weight fan allows varying route
costs. The inspected result does not settle refinements of the present
circuit quotient. Neither that paper nor these finite traces proves
or refutes the desired polylogarithmic refinement bound.

The next step remains a structural count for the full two-dimensional
refinement fan, preferably using the shell and cost constraints
together. Further one-dimensional sweeps or a larger numerical growth
fit cannot establish the unrestricted goal.


## The exact prior-art boundary for box residuals

The box residual identity was already proved as (RC1)--(RC3) in the
[fourteenth feasibility pass](interleaving_research.md#fourteenth-feasibility-pass-the-exact-boundary-of-residual-duality).
Its exact primary precedent is now identified:
Huneke--Ulrich, *Liaison of monomial ideals*,
[Lemma 3.1 in the November 25, 2006 manuscript, Section 3, p. 7](https://www.math.purdue.edu/~bulrich/papers/monomiallinkagenewnew06.pdf).
The definitions and full proof were inspected. The
[published paper](https://doi.org/10.1112/blms/bdl008) is from 2007;
numbering here follows the inspected manuscript.
The lemma states that linkage commutes with taking initials when the
linking polynomials have regular initial monomials.

In the present notation, if $Q\subset I$ and
$\operatorname{in}_wQ=Q_0=(X^{D_1},Y^{D_2},Z^{D_3})$, put
$E=\prod_i\{0,\ldots,D_i-1\}$, $\rho=(D_i-1)_i$, $M=|E|$,
$n=\dim R/I$ and $k=M-n$. The established identity gives

$$
\operatorname{in}_w(Q:I)=Q_0:\operatorname{in}_wI,
\qquad
\operatorname{Std}_w(Q:I)=\rho-(E\setminus\operatorname{Std}_wI).
                                                               \tag{Q10}
$$

For separate monic univariate generators of $Q$, this identifies the
entire positive fans, including nonreduced cases. On a box-initial cone
of a coupled parent, it identifies the restricted fans. This is established
linkage, not a new general theorem. Its original proof uses initial-product
containment and equality of complementary finite lengths.

The [earlier corank screen (D9)--(D10)](cographic_order_ideals.md#a-classical-corank-screen-for-a-higher-dimensional-extension)
already proves $O(n^{3/2})$ for the equal odd-period three-circuit family
over every field. That estimate is retained as a prior classical-volume
exclusion; it is not an advance of this continuation. In an equal-period
$q$ box, the cube count is exactly $q^3-(q-3)^3\le3k$, since each
coordinate span of its eight offsets is three. This small-box fact does
not extend to arbitrary covers.

### Nonbox degeneration has an exact socle defect

The box restriction is essential even before costs are considered. Let
$B$ be a full finite group algebra in its monomial generators, fix a
positive cover cone with standard downset $E$, and write
$B_0=K[X,Y,Z]/\operatorname{in}Q$. Let $\sigma$ be the number of
maximal monomials of $E$. The augmentation ideal $H\subset B$ has
codimension one, and its annihilator is the line spanned by the group
norm $\Theta=\sum_{g\in G}g$. This holds in every characteristic:
annihilation by each $X_i-1$ says that coefficients are constant on the
group generated by the $X_i$.

In $E$ coordinates, $\Theta=\sum_{a\in E}X^a$. At a generic weight
its initial annihilator is thus one-dimensional, spanned by the unique
most expensive monomial of $E$. In contrast, the initial augmentation
ideal is the maximal monomial ideal of $B_0$, whose annihilator is its
socle, of dimension $\sigma$. Consequently

$$
\dim\operatorname{ann}_{B_0}(\operatorname{in}H)
 -\dim\operatorname{in}(\operatorname{ann}_B H)
 =\sigma-1.                                            \tag{Q11}
$$

A finite monomial downset has $\sigma=1$ exactly when it is a box:
with unique maximal exponent $\rho$, every element lies below $\rho$
and down-closure supplies the entire box. Thus every nonbox group cover
has a quotient witnessing failure of the annihilator/initial interchange.
This is an elementary obstruction to extending the established box
linkage argument. It does not say that every particular residual has a
defect, that the defect counts fan regions, or that the monomial socle
is a previously unknown invariant.


## A uniform surface bound for translated-circuit filters

The useful uniform control is a boundary bound, rather than a corank-only
size assertion. Let $E\subset\mathbb N_0^d$ be any finite downset, and
let $d_1,\ldots,d_h\in\mathbb Z^d$ be fixed shifts. Write

$$
p_i=\sum_j\max(d_{j,i},0),\quad
m_i=\sum_j\max(-d_{j,i},0),\quad
F_i=|\pi_i(E)|,
$$

where $\pi_i$ deletes coordinate $i$. If a base $a\in E$ does not
admit all corners $a+\sum_j\epsilon_jd_j$, it lies in the explicit union

$$
\{a\in E:\exists i,\ a_i<m_i\}
 \ \cup\ \{a\in E:a+p\notin E\}.                         \tag{Q12a}
$$

Indeed outside that union all corners are nonnegative and coordinatewise
at most $a+p\in E$, so down-closure includes them. The first set has
at most $\sum_i m_iF_i$ points, by counting its first $m_i$ coordinate
slices. For the second, move from $a$ to $a+p$ one positive coordinate
step at a time. At its first exit from $E$, the preceding point lies on
the top of a coordinate fiber. For each specified step this is an
injective charge to that coordinate boundary, which has exactly $F_i$
points. Thus the second set has at most $\sum_i p_iF_i$ points.
Consequently the full surviving set satisfies

$$
|\text{candidates}|\le
\min\left(|E|,\sum_i\Bigl(\sum_j|d_{j,i}|\Bigr)F_i\right).
                                                               \tag{Q12b}
$$

For the three oriented cube shifts in (Q9), every coordinate's total
absolute variation is three, giving $|\text{candidates}|\le3(F_1+F_2+F_3)$
in every sign sector. The proof is the elementary coordinate-boundary
charging argument. It does not assume convexity, a box cover, a bounded
socle count, or characteristic-independent ranks of other circuits.
The conclusion supplies a union of boundedly many bottom and top fiber
layers; it does not bound their projection sizes by $k$, or count the
coherent bases supported there.

## A cube-only candidate bound fails in the simultaneous-period range

There is also a uniform counterexample to using the small-box candidate
count without its geometry. For every integer $q\ge5$, set

$$
(r,s,t)=(2q+1,q,q+1),\quad
M=q(q+1)(2q+1),\quad n=2q^2(q-1),\quad k=q(5q+1).
$$

These periods are pairwise coprime. With $L=2n+1$, take
$\alpha=tL^3+1$, $\beta=sL^3+L$. Let

$$
T_q=\{(u,0):0\le u\le q\}\cup\{(0,v):1\le v\le q\},
\quad
E=\{(i+su,j+tv,0):0\le i<s,\ 0\le j<t,\ (u,v)\in T_q\}.
                                                               \tag{Q13}
$$

Choose $\gamma>\max_{(a,b,0)\in E}(\alpha a+\beta b)$.
Then $E$ is the full-period standard set at $(\alpha,\beta,\gamma)$.
To prove it, first use only $X,Y$. Equal labels force exponent
differences $(si,tj)$, with $i\equiv j\pmod r$, since
$s\equiv-t\pmod r$. The reduced generators $X^s,Y^t$ are inverse
cyclic generators of order $r$. A representative of label $p$ therefore
costs the smaller of $p\alpha s$ and $(r-p)\beta t$; using both
generators is more expensive because their product is $1$.
The strict threshold

$$
\left\lfloor\frac{r\beta t}{\alpha s+\beta t}\right\rfloor=q
$$

After clearing the positive denominator, the lower and upper strict
inequalities have differences $q(q+1)L^3+(q+1)^2L-q^2>0$ and
$q(q+1)(L^3+1-L)>0$, respectively. This gives the two arms in (Q13). The $s t$ fixed exponent residues give
the complete set of $M$ labels. The chosen $\gamma$ makes any monomial
using $Z$ more expensive than its label's representative in $E$.
All cover inequalities are strict, so this describes an open cover cone.

Every oriented cube in (Q9) has positive span three in the $Z$ coordinate;
none fits in $E$. All $M$ monomials survive, and

$$
\frac{|\text{cube candidates}|}{k}
 =\frac{(q+1)(2q+1)}{5q+1}\longrightarrow\infty.          \tag{Q14}
$$

Nevertheless **every** quotient standard set of size $n$ inside $E$
reaches all three period thresholds. The counts of cover monomials with
$a<s$ or $b<t$ are each $q(q+1)^2<n$ for $q\ge5$.
The number with $a+b<r$ is at most $r(r+1)/2<n$.
Since $c=0$, these are exactly the three relevant pair-degree conditions.
Thus simultaneous threshold crossing alone cannot imply an $O(k)$
cube-candidate bound.

This is not a lower-bound construction for $G_3$: these states all omit
$Z$, so they are planar. The cover has only two maximal monomials,
and the joint upper-shell bound still gives $O(k\log^2(k+1))$.
There is also a general bound for these planar states. For any fixed
length-$n$ ideal $I\subset K[X,Y,Z]$, all its initial standard sets that
omit at least one variable number only

$$
O(n\log(n+1)).                                         \tag{Q15}
$$

For example, a basis omitting $Z$ makes the map $K[X,Y]\to R/I$
surjective. Its fixed kernel $I\cap K[X,Y]$ has colength $n$, and
every such basis is a minimum monomial basis of that planar quotient
at the restricted price $(w_X,w_Y)$. Apply the proved all-field planar
bound and sum over the three omitted-variable choices. Thus a
polynomial separation for the unrestricted target must come from
standard sets containing all three variables, not merely all three
pair-degree thresholds. This is an elimination consequence of the
existing planar theorem, not a new planar count.

Neither the joint filter nor a bound restricted to genuinely
three-variable states is contradicted. The next counting argument must
use that distinction and more than the local cubes; the cube-only
shortcut is closed by a proof, not a finite growth fit.

### Verification scope of this continuation

Norm-quotient frontiers check (Q11) at all nine existing cover prices
over $\mathbb Q$ and $\mathbb F_2$, including box controls with zero
defect. Direct corner membership and separate bottom/top counts check
(Q12) in all four positive sign sectors on nineteen downsets: nine
existing periodic covers, five boxes and five thin covers. Twenty box
checks also verify the exact cube count. Independent shortest paths
verify (Q13) at $q=5,6,7,11,17$; ordinary quotient frontiers and strict
borders check the $q=5$ example in three fields. These checks support
the implementations; the uniform assertions use the proofs above.

The candidate geometry is now bounded by projections, and the stronger
corank-only cube assertion is disproved. The joint shell/circuit count
for genuinely three-variable states and the unrestricted extremal goal
remain open. The exact box linkage input is published prior art; there
is no new sharp exponent, separating family or completed-paper claim.


## A socle-free upper-shell bound

The number of maximal monomials in (Q6) can be removed. For every fixed
$d\ge2$, every nonempty finite downset $E\subset\mathbb N_0^d$ of size
$M$, and every integer $1\le k\le M$, put

$$
u_E(a)=|\{b\in E:b\ge a\}|,\qquad
U_k(E)=\{a\in E:u_E(a)\le k\}.
$$

Then

$$
|U_k(E)|\le C_d M^{(d-1)/d} k^{1/d}
                    (1+\log(M+1))^{(d-1)/d}.             \tag{Q16}
$$

This is an internally proved combinatorial bound and a supporting tool
for the refinement problem. It does not count fan regions or establish
novelty of a published theorem. The proof uses elementary dyadic covering,
coordinate-fiber boundaries and Hölder's inequality. The logarithmic
exponent in this particular $M,k$ form is necessary uniformly; an explicit
monomial-algebra family below shows why. This is not a claim of sharpness
at every fixed pair $(M,k)$; for instance $k=1$ also admits the classical
log-free socle packing bound.

### A small outer cover with bounded total volume

Use the integer intervals
$I_j=\{2^j-1,\ldots,2^{j+1}-2\}$, of size $2^j$, and their
product cells. Let $\mathcal J$ be the indices of cells meeting $E$,
and $\mathcal A$ its coordinatewise maximal indices. If $j\in\mathcal J$,
its lower corner $(2^{j_i}-1)_i$ belongs to $E$. Its anchored divisor
box has $2^{\sum j_i}$ points, so $\sum j_i\le L=\lfloor\log_2M\rfloor$.
Maximal indices with the same first $d-1$ coordinates cannot differ in
the last. Hence

$$
m:=|\mathcal A|\le\binom{L+d-1}{d-1}.
$$

For each $j\in\mathcal A$, take the **outer** anchored box
$B_j=\prod_i\{0,\ldots,2^{j_i+1}-2\}$, with volume $V_j$.
These boxes cover $E$, since any occupied cell lies below a maximal
occupied cell. They need not be contained in $E$.

Their total volume nevertheless satisfies

$$
\sum_{j\in\mathcal A} V_j\le8^dM.                        \tag{Q17}
$$

To see this, let $h_i=\max(j_i-1,0)$. The whole predecessor cell
$\prod_i I_{h_i}$ lies below the occupied cell's lower corner and is
therefore inside $E$. Each outer side is at most four times the
corresponding predecessor-cell side, so $V_j\le4^d\prod_i|I_{h_i}|$.
Different predecessor cells are disjoint, and clamping the indices zero
and one can repeat any cell at most $2^d$ times. Summing proves (Q17).
This total-volume control matters: covering by polylogarithmically many
boxes of volume $O(M)$ separately would lose an extra logarithmic factor.

### An adaptive integer exit step

Let a downset $F$ lie in a box with integer sides $b_1,\ldots,b_d$
and volume $V=\prod_i b_i$. Set

$$
t=(k/V)^{1/d},\qquad h_i=\lfloor tb_i\rfloor.
$$

If $a+h\in F$, the entire integer box $[a,a+h]$ lies in $F$, with
$\prod_i(h_i+1)>t^dV=k$ elements. Thus $U_k(F)$ lies in the exit
boundary $\{a\in F:a+h\notin F\}$. Move by $h_i$ positive unit
steps in each coordinate. The first-exit charge used in (Q12) gives

$$
|U_k(F)|\le\sum_i h_i\prod_{j\ne i}b_j
 \le d k^{1/d}V^{(d-1)/d}.                               \tag{Q18}
$$

Zero steps in thin coordinates are allowed and are essential: rounding
every side upward would add a boundary term as large as $V$. The argument
also works if $k\ge V$, when the bound may simply be trivial.

Apply (Q18) to $F_j=E\cap B_j$. Each original shell point in $B_j$
has at most $k$ upper points in $F_j$, so

$$
|U_k(E)|\le d k^{1/d}\sum_j V_j^{(d-1)/d}
 \le d k^{1/d}m^{1/d}\Bigl(\sum_jV_j\Bigr)^{(d-1)/d}.
$$

Equations (Q17) and the bound on $m$ prove (Q16), including its
uniformity over the shape and socle count of $E$.

### The logarithm cannot be removed in this form

For $L\ge0$, define the dyadic downset

$$
D_L=\{a\ge0:\sum_i\lfloor\log_2(a_i+1)\rfloor\le L\},
\qquad H=2^L.
$$

Its exact cardinality and top-cell mass are

$$
M_L=\sum_{s=0}^L2^s\binom{s+d-1}{d-1},\qquad
T_L=2^L\binom{L+d-1}{d-1},\qquad T_L\le M_L<2T_L.       \tag{Q19}
$$

If $a$ belongs to a cell indexed by $j$ with $\sum j_i=L$, every
$b\in D_L$ above it stays in that same cell: its dyadic indices are
at least $j$ coordinatewise and cannot have larger total sum. There are
only $2^L=H$ points in the cell, so $u_{D_L}(a)\le H$.
Consequently

$$
|U_H(D_L)|\ge T_L>M_L/2.
$$

For fixed $d$, $M_L=\Theta_d(H L^{d-1})$ and
$\log M_L=\Theta_d(L)$. Any uniform estimate with the powers
$M^{(d-1)/d}k^{1/d}$ and a factor $(1+\log M)^\beta$ therefore
requires $\beta\ge(d-1)/d$. The witness is a realizable monomial
quotient over every field, but it has only one initial ideal. Its shell
size is not an extremal fan lower bound. No claim that $D_L$ is a full
periodic group cover is used.

### Algebraic and fixed-cover consequences

For any finite quotient algebra $A=K[x_1,\ldots,x_d]/I$ with a
monomial downset basis $E$, multiplication by $x^a$, $a\in E$, has
rank at least $u_E(a)$. Indeed the distinct standard products
$x^{a+b}$, for $b\in E$ with $a+b\in E$, are independent, and
$b\mapsto a+b$ bijects this set with the upper interval above $a$.
Thus (Q16) bounds how many monomials in **any** such basis can have
multiplication rank at most $k$, uniformly over fields and coefficients.
For a monomial quotient rank equals $u_E(a)$, and (Q19) supplies the
logarithmic obstruction in this algebraic formulation too. The assertion
is about basis monomials, not the span of all arbitrary elements of
small multiplication rank.

For a fixed three-dimensional periodic cover of size $M$ and corank $k$,
all actual omissions lie in the same coefficient-independent set

$$
|U_k(E)|=O\!\left(M^{2/3}k^{1/3}\log^{2/3}(M+1)\right).  \tag{Q20}
$$

For comparable periods tending to infinity, $k=O(M^{2/3})$, so this
is $O(M^{8/9}\log^{2/3}(M+1))=o(M)$, without a restriction on prices
inside the cover cone. More generally it is $o(M)$ whenever
$(k/M)\log^2(M+1)\to0$. Restricting the known inverse-label maximum
basis problem to this shell preserves its optimum at every price in
the cover cone. Given the shell, greedy selection uses at most its
cardinality in independence tests. This is a reduction in coefficient
rank tests, not a sublinear construction of the cover or an unproved
overall running-time improvement. The present checker uses direct
counts; standard orthant-counting methods are a separate implementation
choice. The joint cube/shell ground can only be smaller.

The full two-dimensional refinement count remains missing. A sublinear
candidate list alone gives neither a sharp state count nor a separating
family, and the cover itself varies across prices. Every original
extremal and publication criterion remains in force.


### General downset geometry alone cannot improve the power

There is an obstruction even after combining the shell with the actual
cube test and retaining all three variables. For integers $R\ge2,L\ge1$,
put

$$
E_{R,L}=\{(a,b,c)\ge0:a+\lfloor b/L\rfloor+\lfloor c/L\rfloor<R\}.
$$

This is a three-dimensional monomial downset containing $X,Y,Z$, with

$$
M=L^2\binom{R+2}{3},\qquad k=L^2.
$$

Its top $X$ layer consists of $a=R-i-j-1$, where
$i=\lfloor b/L\rfloor$, $j=\lfloor c/L\rfloor$, $i+j<R$.
It has $T=L^2\binom{R+1}{2}$ points. An upper point of $E_{R,L}$
above such a point must keep $a,i,j$ fixed; its remaining $b,c$ rectangle
has at most $L^2=k$ points. Thus every top-layer point is in $U_k(E)$.
In the all-positive cube sign sector, the eight corners include
$a+2e_X$, which is outside the downset at every top-layer point.
All $T$ points therefore survive the **joint** shell/cube filter.

As $R\to\infty$,

$$
T=\Theta(M^{2/3}k^{1/3}).                               \tag{Q21}
$$

Taking $L=R^3$ gives $k=\Theta(M^{2/3})$ and
$T=\Theta(M^{8/9})$, the same power as the upper-shell bound at the
periodic candidate's natural corank scale. Hence general downset
geometry, even with all variables present and both filters imposed,
cannot give an $O(k\,\operatorname{polylog}M)$ candidate bound or a
smaller uniform power of $M$ at that scale.

These examples are ordinary monomial quotients and have one initial
ideal. They are **not** asserted to be standard domains of a periodic
group cover. Thus this is not a counterexample to a bound using the
actual period labels, tiling and representative-cost constraints.
It pinpoints the information the next proof must use beyond abstract
shell and boundary geometry.

### Prior comparison and exact checks

Related polylogarithmic aspect families are established in
[van der Hoeven--Lecerf, inspected PDF, Lemma 4.1, Lemma 4.2 and Corollary 4.3](https://www.texmacs.org/joris/genamp/genamp.pdf),
Sections 4.1--4.2, pp. 8--10. Those statements and proofs were reopened:
they balance shifted exponents for heterogeneous quasi-reduction.
They are relevant prior art for the dyadic method, not a checked
statement of this upper-interval distribution estimate. A complete
comparison with discrete boundary and multiplication-rank bounds is
still required; unsuccessful searches do not establish novelty.
[Hausbrandt--Ruzika v2, Observation 3.1, Theorems 3.2/3.6 and Corollary 3.4](https://arxiv.org/html/2503.08178v2#S3)
were also rechecked. Their general pair-comparison arrangement gives
$O(s^4)$ regions on a ground of $s$ elements with two price parameters.
Putting $s=O(M^{8/9}\log^{2/3}M)$ in that bound does not improve the
current unrestricted count, even before summing cover cones. The shell
reduction still needs a quotient-sensitive two-dimensional counting
argument.

The checker independently constructs the outer boxes and their fully
contained predecessor cells, computes every upper interval by dominance
counts, and compares it to the rank of the monomial multiplication
matrix through its distinct nonzero columns. Adaptive steps are computed
by exact integer roots and checked against actual exit boundaries.
All thresholds are tested for all 377 downsets of lengths at most nine
in dimension two, seven in dimension three, and five in dimension four,
giving 2,211 threshold checks and 4,161 adaptive box-boundary checks. Separate periodic covers and
dyadic crosses verify the same mechanism; the infinite family has exact
integer cardinality formulas. Five staircase-block fixtures directly check
the joint shell/cube obstruction with all three variables, and five
additional exact formulas check its natural-corank scale. These data support the uniform proof,
not a completed priority audit or the main goal.

## Period gaps give regional bounds even when every cube fails

The preceding obstruction used arbitrary downsets. Actual periodic labels
give more information. The following bounds use that information, and
apply to an unbounded-price part of the original three-generator family.
They also show that ineffective local pruning can occur inside a region
whose entire initial fan is only quadratic.

### A gap-dependent common box

Suppose $r>t>s\ge3$ are odd and pairwise coprime. Define

$$
h=(r-s)/2,\quad e=(r-t)/2,\quad g=\gcd(h,e),\quad
h_0=h/g,\quad e_0=e/g.
$$

In $B$, and hence in $A$, the following identities hold over every field:

$$
X^{e_0s}=Y^{h_0t},\qquad
Z^h=X^{s+h}Y^h,\qquad
Z^e=X^eY^{t+e}.                                        \tag{Q22}
$$

For the first identity, $e_0s-h_0t=(e-h)r/g$; its other two
label coordinates are multiples of their periods. The other identities
use $s+2h=r$ and $t+2e=r$. These are group-label identities,
requiring no squarefreeness or division by the characteristic.

Let $\Omega_Y$ be the region of positive prices satisfying

$$
h_0t w_Y>e_0s w_X,\qquad
h w_Z>(s+h)w_X+h w_Y.
$$

Every standard set of $A$ in this region lies in the **same** box

$$
0\le a<rs,\quad 0\le b<h_0t,\quad 0\le c<h,
\qquad \text{volume }h h_0 M.                           \tag{Q23}
$$

Indeed the first two costly heads in (Q22) are $Y^{h_0t}$ and
$Z^h$, while $X^{rs}=1$ excludes $a\ge rs$. Similarly, in
$\Omega_X$ defined by

$$
e_0s w_X>h_0t w_Y,\qquad
e w_Z>e w_X+(t+e)w_Y,
$$

one common box has sides $e_0s,rt,e$, and volume $e e_0M$.
Both regions are nonempty: choose the first two prices satisfying
their first inequality and then choose $w_Z$ sufficiently large.

The number of initial ideals realized anywhere in each whole region is

$$
g_{\Omega_Y}(A)=O\!\left(n^2\sqrt{h h_0}\right),\qquad
g_{\Omega_X}(A)=O\!\left(n^2\sqrt{e e_0}\right).          \tag{Q24}
$$

To prove this, take the represented monomial matroid on the relevant
common box, of rank $n$. Its basis-profile polytope lies in a box of
volume at most $n^3$ times the ground-box volume. Each positive
initial basis in the specified region uniquely minimizes its profile
there, and hence determines a distinct vertex. Adding a unit cube
handles a lower-dimensional polytope, keeps the enclosing-box volume
of order $n^3$ times the ground-box volume, and refines its normal fan.
The classical lattice-vertex bound
then gives $O(n^{3/2}\sqrt{h h_0 M})$ in the first region.
Since $M\le8n$, this proves (Q24); the other region is identical.

The basis-profile dictionary and volume-to-vertex step are established
in [Babson--Onn--Thomas, Theorem 2.4 and the proof of Corollary 2.5](https://arxiv.org/pdf/math/0207135).
This is a classical geometric corollary using explicit period relations,
not a new lattice-polytope theorem. It also applies to any quotient of
$B$ of dimension comparable to $M$. Bounded period gaps make (Q24)
quadratic; polylogarithmic gaps give $n^2$ times a polylogarithm.
Neither statement covers the complement of the specified price regions.
The result must not be multiplied by a count of cover cones: the common
box already handles every cover and every quotient refinement throughout
each entire region.

### An all-field obstruction with all three generators present

For every integer $q\ge1$, set

$$
(r,s,t)=(6q+7,6q+1,6q+5),\quad D=2rst+1,\quad H=2q+3,
$$

and choose the positive integer prices

$$
(w_X,w_Y,w_Z)=(D^3+1,D^3+D,H D^3+D^2).                  \tag{Q25}
$$

These periods are odd and pairwise coprime: their differences are
$6,2,4$, and $s\equiv1\pmod3$. Here $h=3,e=1,h_0=3,e_0=1$.
In particular

$$
Z^3=X^{s+3}Y^3,\quad X^s=Y^{3t},\quad Z=XY^{t+1}.
$$

The first two give positive strict price margins

$$
3w_Z-(s+3)w_X-3w_Y=2D^3+3D^2-3D-(s+3)>0,
\quad 3t w_Y-s w_X>0.
$$

Thus every cover staircase $E$ and every quotient standard set lie in
the box $a<rs,b<3t,c<3$, of volume $9M$. All relevant costs
are distinct: the three lower base-$D$ digits separate exponents in
the pure-axis box, and their differences have magnitude less than
$D^3$. The prices are inside $\Omega_Y$, with
$w_Z/w_X\sim2q\to\infty$. Formula (Q24) bounds the entire regional
fan by $O(n^2)$, even though the price ratio is unbounded.

Every quotient standard set at (Q25) contains

$$
X,\quad Y,\quad Z,\quad X^r,\quad Y^t                 \tag{Q26}
$$

over **every field**. Here is a direct independence proof, including
characteristics dividing a period. Use the monic tensor basis
$U^iV^jW^k$ with $i<r-1,j<s-1,k<t-1$.

* Before $Z$, the monomials are exactly $X^aY^b$ with
  $a+b\le H$. Because $H<s-1$, their residues
  $U^{a+b}V^aW^b$ are distinct tensor-basis vectors. The residue
  $VW$ of $Z$ is different from all of them. This also proves that
  $X$ and $Y$ are selected.
* Before $Y^t$, every monomial has $a+b\le t$. Those at
  $U$-degree $t$ are precisely $X^aY^{t-a}$, $1\le a\le t$:
  $X$ is cheaper than $Y$, and any positive $Z$ exponent would
  cost more than $Y^t$. There is no $U$ reduction because
  $t=r-2$. It suffices to show
  $1\notin\operatorname{span}\{V^aW^{t-a}:1\le a\le t\}$.
  In a proposed relation, only $a=1$ has $W$ exponent $t-1$,
  which reduces to $-\sum_{j=0}^{t-2}W^j$. At $W$-degree
  $t-2$ the coefficient is $c_2V^2-c_1V$; independence gives
  $c_1=c_2=0$. The remaining nonzero $W$ degrees each occur
  once, and $V$ is a unit, so their coefficients vanish. The
  constant equation is $c_t V^t=1$, impossible since
  $V^t=V^4$ is a nonconstant tensor-basis vector.
* Before $X^r$, every monomial has $a+b\le r-1$ and
  $c\le2$. The only $U$ reduction is
  $U^{r-1}=-\sum_{i=0}^{r-2}U^i$. Let $f\in K[V,W]$
  be the combined coefficient before this reduction at degree
  $r-1$. To cancel $U$-degrees $1$ and $2$, $f$ must belong
  both to the coefficient spans from $a+b=1$ and from $a+b=2$,
  with $c\le2$. The former uses $V,W$ monomials of odd total
  degree $1,3,5$; the latter uses even degrees $2,4,6$.
  All their individual exponents are at most $4$, below both
  monic-reduction thresholds. Their supports are disjoint, so
  $f=0$. At $U$-degree zero the remaining earlier span is
  $\operatorname{span}\{1,VW,V^2W^2\}$. But $X^r$ has residue
  $V^r=V^6$, outside this span. For $q=1$, $V^6$ reduces to
  $-\sum_{i=0}^5V^i$, still nonconstant in every field; for
  $q\ge2$ it is an ordinary nonconstant basis vector.

Each argument considers **all cheaper monomials**, rather than just
one chosen frontier. Independence from that span forces greedy
selection and hence membership in the initial standard set. The
strict inequalities persist on open neighborhoods of (Q25).
The members $X^r,Y^t$ reach the three thresholds
$a+b\ge r,a+c\ge s,b+c\ge t$, and $X,Y,Z$ ensure that no
variable is omitted. Thus this obstruction survives both earlier
restrictions on a possible separating family.

Nevertheless every oriented local cube fails to fit in $E$. Each of
the three cube directions has $Z$ increment of absolute value one;
the eight corners span three $Z$ levels of difference. The cover has
only levels $0,1,2$. Hence the cube-only filter leaves all $M$ points.
Since

$$
M=\Theta(q^3),\quad n=\Theta(q^3),\quad
k=M-n=\Theta(q^2),\quad M/k\longrightarrow\infty,        \tag{Q27}
$$

there is no uniform $O(k)$ cube-only bound even when all generators
are active and every period threshold is reached. The joint upper-shell
bound remains valid; this does not supply a lower bound on omissions or
on the number of initial ideals.

The cases $q=1,2,6$ give genuine prime sets $\{7,11,13\}$,
$\{13,17,19\}$ and $\{37,41,43\}$. Infinitely many prime
instances are **not** asserted. The infinite family uses pairwise
coprime periods, including composite ones. Its identities explain an
arithmetic restriction on these exponent lattices, not a consequence
for prime distribution.

### Verification and the remaining problem

Eight odd coprime period triples check (Q22) by both group labels and
independent raw monic division over $\mathbb Q,\mathbb F_2,
\mathbb F_3,\mathbb F_5$. Sixteen regional price controls verify
strict inequalities and common-box volumes; four also construct the
entire shortest-representative cover. Seven family parameters check
the exact integer formulas. For $q=1,2,3$, ninety prefix tests over
$\mathbb Q$ and characteristics $2,3,5,7,13$ verify that every
target in (Q26) increases the rank of the full cheaper-monomial span.
These include characteristics dividing periods.

For $q=1,2$, the complete binary quotient bases of lengths $720$ and
$3456$ are computed directly, and independently recovered as complements
of maximum-cost inverse-label dual bases. Ordinary monic division then
verifies every standard column and all $35$ minimal-border relations,
with strictly cheaper tails. Exact omissions, integer prices and border
certificates are saved in the report. Direct corner tests confirm that
the full covers of sizes $1001$ and $4199$ survive cube-only pruning.
No finite fan count is inferred from these two bases.

This closes another arithmetic region as a source of polynomial fan
separation, while showing why adding an active-generator assumption
cannot repair the cube-only shortcut. The complement of
$\Omega_X\cup\Omega_Y$, arbitrary period gaps and the unrestricted
$G_3(n)$ target remain open. A next proof must either count the
joint-shell refinements across those remaining regions or construct a
different realizable source of many simultaneous dependencies. Exact
priority, a consequential unrestricted theorem and the paper endpoint
remain unestablished; no new manuscript follows from this corollary.

## A global gap bound and why small cover boxes cannot make it quadratic

The preceding price regions do not exhaust the fan. There is, however,
a global bound from two period identities. It still leaves a power gap
above the generic scale, and its classical geometric proof does not
establish novelty or complete the main goal.

### Four sectors cover every positive weight

Continue with odd coprime periods $r>t>s\ge3$. In addition to the
previous $h,e,g,h_0,e_0$, put $\delta=h-e=(t-s)/2$ and
$\delta_0=\delta/g$. The exact identities

$$
X^{e_0s}=Y^{h_0t},\qquad
X^{\delta_0r}=Z^{h_0t},\qquad
Y^{\delta_0r}Z^{e_0s}=1                                \tag{Q28}
$$

hold in $B$ and every quotient. The first two were obtained by
cancelling shared period coordinates. For the third,
$\delta_0r+e_0s=h_0t$. These are equalities of integral labels,
valid in ramified characteristics as well.

Write $D_X,D_Y,D_Z$ for the exact pure-axis degrees (P3), and
$\alpha=e_0s,\beta=h_0t,\gamma=\delta_0r$. The first two
relations partition positive weight space into four strict sectors.
Each has one common box, obtained by intersecting the pure-axis box
with the following head bounds:

| Expensive head in the first relation | Expensive head in the second | Additional bounds |
|---|---|---|
| $Y^\beta$ | $Z^\beta$ | $b<\beta,c<\beta$ |
| $X^\alpha$ | $Z^\beta$ | $a<\alpha,c<\beta$ |
| $Y^\beta$ | $X^\gamma$ | $b<\beta,a<\gamma$ |
| $X^\alpha$ | $X^\gamma$ | $a<\min(\alpha,\gamma)$ |

Let $V_1,\ldots,V_4$ be these four box volumes, after clipping each
side by its corresponding $D_i$. The common-box profile argument gives
the **global**, all-field bound

$$
g(A)=O\!\left(n^{3/2}\sum_{j=1}^4\sqrt{V_j}\right).    \tag{Q29}
$$

Only four fixed sectors are summed, rather than an unknown number of
periodic-cover cones. A standard set on a comparison wall is represented
by a generic perturbation and adds no extra case. Each counted basis
is a unique minimum on the common ground for its sector. The argument
uses the established state/profile dictionary and lattice-vertex bound
of [Babson--Onn--Thomas, Theorem 2.4 and Corollary 2.5](https://arxiv.org/pdf/math/0207135).
It does not assert a new vertex theorem or a field-independent equality
of fans.

Suppose the periods are comparable, with $r/s$ bounded. Put
$H=h_0$ and use $R=s$ as a size parameter, so $n=\Theta(R^3)$
and $H=O(R)$. In the first three sectors the box volume is
$O(H^2R^4)$; in the last it is $O(HR^5)$. Therefore

$$
g(A)=O\!\left(\sqrt H\,n^{7/3}\right),                \tag{Q30}
$$

with constants depending only on the period comparability bound.
For bounded normalized gaps this is $O(n^{7/3})$ across **all**
positive weights, including every earlier intermediate region.
The first three sectors separately have
$O(H n^{13/6})$ states. The older pure-axis box alone gives
$O(n^{5/2})$; (Q30) saves a power when $H=o(n^{1/3})$.
Fixed ratios between the period gaps are enough: the actual gaps need
not stay bounded. For example three-term arithmetic progressions have
$(h_0,e_0,\delta_0)=(2,1,1)$ whenever the odd periods are coprime.
No statement about the infinitude or distribution of such prime triples
is used.

More generally (Q29) holds for every quotient of $B$ if the pure-axis
sides are replaced by $rs,rt,st$ and its actual dimension is used
as $n$. If that dimension is bounded below by a positive constant
times $M$, (Q30) holds with the corresponding constant. Thus the
global bound is not limited to the particular three-circuit quotient.
Its exact priority among dense periodic-quotient counts remains open;
the restricted exponent does not improve the unrestricted $G_3$ bound.

### A single shortest cover really can need the larger box

The fourth-sector box has volume of order $R^5$ at fixed normalized
gaps. That scale cannot be reduced to $O(M\operatorname{polylog}M)$
for every actual shortest cover by adapting just its enclosing box.
Here is a uniform proof, without fitting numerical extents.

Take $(r,s,t)=(6q+7,6q+1,6q+5)$ and $D=2M+1$, as before,
but now use

$$
w=(tD^3+1,\ 2D^3+D,\ 2D^3+D^2).                       \tag{Q31}
$$

The shortest-representative cover $E$ contains

$$
X^q,\qquad Y^{q^2},\qquad Z^{q^2}.                    \tag{Q32}
$$

To prove this, first compare the integral score
$\ell(a,b,c)=ta+2b+2c$. For a competitor to $X^q$ of score
at most $tq$, its labels force

$$
b+c=tk,\quad b=q-a+ri,\quad c=q-a+sj,
\quad a+2k\le q,
$$

with $i,j,k\ge0$. The nonnegativity of $i,j$ follows from
$0\le q-a\le q<s$. Also $i,j\le k\le q/2$. Substitution gives

$$
s(k-i-j)=2(q-a)+6i-4k.
$$

The right side lies in $[-2q,5q]$, strictly between $-s$ and
$s$, so it is zero. Hence $k=i+j$ and $q-a=2j-i$.
The score bound also gives $q-a\ge2(i+j)$, forcing $i=0$.
All low-score competitors are consequently exactly

$$
(a,b,c)=(q-2j,2j,(s+2)j),\quad0\le j\le\lfloor q/2\rfloor.
$$

They tie in integral score, but the lower price digits in (Q31)
increase by $j(-2+2D+(s+2)D^2)>0$ for $j>0$.

For a competitor to $Y^m$, $m=q^2$, the score bound forces
$b+c=m-tl$ and $0\le a\le2l$, with $0\le l\le m/t<q/6$.
Its other labels give $a+b=m-rp$ and $a+c=sj$, with
$0\le p\le m/r<q/6$ and $j\ge0$. Thus

$$
2a-4l=s(j-p+l)-6p\in[-4l,0].
$$

An integer value $j-p+l\ge1$ makes the right side positive;
a value at most $-1$ makes it less than $-4l$. Therefore
$p=j+l$ and $a=2l-3p\le-l$, forcing $l=p=j=a=0$.
The unique low-score competitor is $Y^m$ itself.

For $Z^m$, the analogous labels are $a+b=rj$ and
$a+c=m-sp$, with $j<q/6$. They give

$$
2a-4l=s(j-p+l)+6j\in[-4l,0].
$$

The same small-multiple argument forces $j-p+l=0$ and $j=0$,
then $a+b=0$ and $a=2l$ force $l=a=b=0$. Only $Z^m$ remains.
The signs and bounds use only integer labels, not coefficient ranks.
Finally the lower-digit cost of each target is less than $D^3$,
whereas one extra unit of integral score costs $D^3$. No higher-score
monomial can defeat a target. This proves (Q32) for the actual prices.

Every coordinate box enclosing this **one** cover must therefore have
volume at least

$$
(q+1)(q^2+1)^2=\Omega(q^5)=\Omega(M^{5/3}).             \tag{Q33}
$$

The fourth-sector cover box is also $O(q^5)$, so this order is
attained up to constants. The assertion is about the full group cover,
not a proof that the quotient basis contains the same axis monomials.
It does not refute a sharper quotient-rank/profile bound or the
sublinear omission-shell estimate.

### The difficult box does not imply a difficult local fan

Prices (Q31) refine the positive integer ray $w_0=(t,2,2)$.
Indeed the lower digits cannot reverse any nonzero integral-score
comparison in the pure-axis ground. All infinitesimal refinements at
this ray are initial ideals of $\operatorname{in}_{w_0}(I)$.
That ideal is homogeneous for the strictly positive grading $w_0$.
The [previous all-field grading theorem](nonnegative_graded_fans.md#all-three-variable-degrees-positive)
therefore bounds their total number by $O(n)$, uniformly in $q$.
This uses the established transitivity of successive initials and
the earlier proved grading bound. It is not a global count of cones
near the ray at a finite distance.

Thus a cover can require a large price-adapted box while the quotient
has only linearly many infinitesimal refinements there. The next count
must use ranks or control how cones meet these arithmetic rays; the
box volume alone loses that information. Proving that every remaining
cone meets a controlled collection of rays would be substantive.
It is **not** established by the examples, and the earlier unrestricted
counterexample to a polylogarithmic projection-plane cover remains valid.

Ten period models check (Q28) by raw monic division over four fields.
Forty strict sector controls verify the clipped boxes, and eight
construct complete covers. Eight parameter choices check (Q33)'s
integer formulas; four complete covers and twelve exhaustive low-score
representation tests independently verify (Q32), including every tie
at $X^q$. All previous evidence remains preserved. The new result is
a global restricted upper bound and a proof-budget obstruction, not
a polynomial fan separation, a novelty certificate or the requested
paper endpoint.

## Exact reproduction

```bash
python3 -B scripts/run_periodic_quotient_checks.py
```

The checker validates shortest-coset binomial borders, the inverse dual
representation in several fields, actual quotient omissions and the
nonbox reflection obstruction. It independently computes the two integer
determinants, modular rank, cleared relation, strict borders and direct
degree-two/three substituted frontiers. Its report records the exact
models and relevant source hashes. The continuation adds complete
generic fixed-cover line traces, independent polynomial checks and
the upper-shell and integral-slope invariants. This is neither a full fan sweep nor
a numerical estimate of asymptotic growth; the current-priority audit
and main endpoint remain incomplete. The residual continuation adds direct
norm-quotient checks of nonbox socle defects, sector-wide projection
budgets and an infinite thin-cover obstruction to cube-only pruning.
The exact box linkage precedent and earlier symmetric-family exclusion
are cross-referenced rather than claimed again. The shell continuation
adds exact dyadic covers, integer exit steps, principal multiplication
counts and logarithmic-obstruction formulas.
The period-gap continuation adds raw monic checks of the arithmetic
identities, regional common boxes, all-cheaper prefix independence and
complete primal/dual bases with strict borders for active slab covers.
The global continuation adds four-sector coverage, exact period identities,
complete low-score representation tests and the large-box obstruction at
one positive arithmetic weight ray.
