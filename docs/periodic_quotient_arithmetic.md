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
given, its construction takes $O(M)$ set-membership operations. No
bound of $O(k)$, polylogarithmic size or sublinear cover construction
is proved. The three sign planes give a constant number of positive
price sectors, but they do not count the quotient's cones.
This extends the earlier four-cycle cost-domination mechanism to the
known simplicial dual representation; it does not assume that all its
circuits are cubes or that the matroid is regular. Bounding the surviving
candidate geometry uniformly is an unresolved part of the full-cone
problem.

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
and main endpoint remain incomplete.
