# Radical lifts of prime digit quotient fans

The prime-digit dependencies do not require a nonreduced algebra or
positive characteristic. For every prime and truncation length, an
explicit selection of roots of unity gives the same complete monomial
rank rule in characteristic zero. Every initial standard set and every
positive initial cone of the unit-power examples transfers to a radical
point ideal of the same length.

In the binary case there is a stronger integral lift. One explicit
integer polynomial preserves the monomial matroid over every field.
Its power-column matrices are totally unimodular, and every reduced
Gröbner basis of its unit-power graph ideal has coefficients in
$\{0,1,-1\}$ and specializes correctly to every characteristic.
These are structural transfer results. The unrestricted count and
polynomial lower separation remain open; exact priority is unsettled.
The simultaneous rank minimum and its complete frequency-set
classification below make the remaining prior-art comparison precise.
The balanced-basis criterion and general cross-field representability
already belong to the sampling and Bhargava-greedoid literature.

## The prime digit rank rule

Let $p$ be prime, $M=p^s$, $0\le n\le M$. The columns
$$
c_a=(\tbinom aj\bmod p)_{0\le j<n},\qquad0\le a<M,
\tag{DL1}
$$
represent powers of $X=1+t$ in $K[t]/t^n$. For $h=M/p$ and
$n=qh+r$, $0\le r<h$, put
$$
l_u=|S\cap(u+h\{0,\ldots,p-1\})|,\qquad U=\{u:l_u\ge q+1\}.
$$
The previously proved [Pascal recursion](oracle_fan_enumeration.md) is
$$
R_{n,M}(S)=\sum_{u=0}^{h-1}\min(q,l_u)+R_{r,h}(U), \tag{DL2}
$$
with $R_{0,M}=0$ and $R_{M,M}(S)=|S|$. Its interpretation as
iterated pointed uniform-matroid composition is
[already recorded](adjacent_rank_fans.md). The digit rule and the
binary series/parallel regularity are established inputs here.

## A Fourier realization with the same ranks

Define frequency sets recursively by
$$
\begin{split}
E_{0,M}&=\varnothing,\qquad E_{M,M}=\{0,\ldots,M-1\},\\
E_{n,M}&=\{a+pb:0\le a<q,\ 0\le b<h\}
       \ \sqcup\ \{q+pb:b\in E_{r,h}\}. \tag{DL3}
\end{split}
$$
Equivalently, $E_{n,M}$ is the base-$p$ digit reversal of
$\{0,\ldots,n-1\}$. Let $L$ have characteristic different from $p$
and contain a primitive $M$th root $\zeta$. Put
$$
v_a=(\zeta^{ia})_{i\in E_{n,M}},\qquad0\le a<M. \tag{DL4}
$$

**Theorem 1.** Every subset of (DL4) has rank $R_{n,M}(S)$.
Its labelled matroid equals (DL1).

**Proof.** Write each column label $a=u+hv$ and each full-block row
frequency as $i=b_0+pb$. Its entry is
$\zeta^{b_0u}(\zeta^p)^{bu}(\zeta^h)^{b_0v}$.
For fixed $b_0<q$, invert the full $h$-point Fourier matrix in $b$.
These rows then separate the fibers $u$, supplying their first $q$
Vandermonde rows on the distinct nodes $(\zeta^h)^v$, up to nonzero
scalars. Their rank on a selected fiber is $\min(q,l_u)$.

The remaining rows have $b_0=q$, $b\in E_{r,h}$. On the kernel of
the first $q$ rows in a fiber, their common extra Vandermonde
functional is nonzero exactly when $l_u\ge q+1$. This follows from
the ordinary $(q+1)$-row Vandermonde determinant. Its image is then
the line with vector
$$
\zeta^{qu}\bigl((\zeta^p)^{bu}\bigr)_{b\in E_{r,h}}.
$$
The residual rank is the smaller Fourier rank on $U$. Induction gives
(DL2). The zero and full cases start it; the full Fourier matrix is
invertible on its distinct nodes. This uses the specified row sets,
not a claim about arbitrary Fourier minors. $\square$

## Actual radical polynomial quotients

Compare $K[X]/((X-1)^n)$ in characteristic $p$ with
$$
A^{\mathrm{red}}=L[X]/G_{n,M}(X),\qquad
G_{n,M}(X)=\prod_{i\in E_{n,M}}(X-\zeta^i). \tag{DL5}
$$
The latter is the algebra of $n$ distinct roots of unity. Powers
repeat modulo $M$ in both algebras, so Theorem 1 identifies ranks of
every finite Laurent monomial list, including repeated residues and
negative exponents.

Map polynomial generators to specified powers of $X$, with one
generator equal to $X$ so the image is the full algebra. The two
relation ideals have the same residue-column matroid with the same
exponent labels. Actual-column greedy gives the same standard set for
every global monomial order. The whole positive fans, their walls and
adjacency agree, rather than only their numbers of cells.

For the earlier unit-power family an explicit radical lift is
$$
(G_{n,M}(X),XY-1,Z-Y^k), \tag{DL6}
$$
with points $(\zeta^i,\zeta^{-i},\zeta^{-ki})$, $i\in E_{n,M}$.
The construction works in any fixed number of polynomial generators.
Finite extensions supply roots over finite fields when necessary;
no realization over a prescribed small prime field is asserted.

Thus a supergeneric lower family among these prime-digit quotients
would already yield a radical characteristic-zero lower family with
no length overhead. The transfer theorem supplies that reduction;
it does not establish such growth. Repeated roots are unnecessary
for the observed ordering patterns throughout this class.

## An integral family with no bad residue primes

There is a simultaneous arithmetic version for every prime $p$.
Take $\mathcal O=\mathbb Z[\zeta]$, with $\zeta$ a primitive $M$th
root, and form the monic quotient $\mathcal O[X]/G_{n,M}$.
It is free of rank $n$. At a residue prime of characteristic
$\ell\ne p$, the image of $\zeta$ still has order $M$: otherwise
its order would divide $M/p$ and its cyclotomic equation would give
$p=0$ in that residue field. At a residue prime over $p$, it becomes
one and $G_{n,M}$ becomes $(X-1)^n$. Theorem 1 and (DL1) therefore
give the same power-column matroid at every residue prime, as well
as in the fraction field.

Every maximal power-basis minor is consequently zero or a unit of
$\mathcal O$: a nonzero minor vanishes at no maximal ideal. All
maximal ideals here have positive residue characteristic, since
$\mathcal O$ is finite over $\mathbb Z$. Cramer's rule now gives
monic reduced bases over $\mathcal O$ for the unit-power graph
ideals. Every nonzero coefficient is a unit, so all residue-prime
specializations preserve both their supports and their leading ideals.
There are no bad residue primes for these bases and fans.

This is a finite flat family with a radical characteristic-zero fiber
and a nilpotent characteristic-$p$ fiber, without losing any of the
coherent monomial orderings. The base ring may require a cyclotomic
extension. The following binary construction removes that extension
and strengthens its unit coefficients to ordinary integer signs.

## A totally unimodular binary lift over the integers

For $n\ge1$ define
$$
F_n(X)=\prod_{j:\,\operatorname{bit}_j(n)=1}(1+X^{2^j})
      =\sum_{a\,\&\,n=a}X^a. \tag{DL7}
$$
It is monic of degree $n$, has constant term one, and reduces to
$(X-1)^n$ in characteristic two. Let $C_{n,M}$ be the integer matrix
of $1,X,\ldots,X^{M-1}$ modulo $F_n$, in the power basis of length
$n$, for a power of two $M\ge n$.

**Theorem 2.** The matrix $C_{n,M}$ is totally unimodular. Its column
matroid over every field is the binary Pascal matroid (DL1).

First prove field independence by induction on $M$. If $n=M$, the
matrix is the identity. If $n\le M/2$, columns repeat in parallel
blocks of period at most $M/2$: for nonpower $n$, $F_n$ divides
$X^{M/2}-1$; at $n=M/2$ it is $1+X^n$ and the repeated block has
a minus sign. This reduces to a smaller period. For even $n=2m$,
$F_n(X)=F_m(X^2)$; even and odd powers give two independent summands
with smaller matrix $C_{m,M/2}$.

For odd $M/2<n<M$, put $r=M-n<M/2$. Binary complementation gives
$$
F_n(X)F_r(-X)=1-X^M. \tag{DL8}
$$
The residue-column matroids modulo these two factors are dual on the
same $M$ labels. To justify the cyclic-code step, in
$L[X]/(X^M-1)$ the kernel of reduction modulo a monic divisor $f$
is $(f)$. For complementary reciprocal divisors $f,g$, the ideal
$(f)$ is the dot-product orthogonal complement of $(g)$: the pairing
is the constant coefficient of $a(X)b(X^{-1})$; inclusion follows
from $fg=0$ and complementary dimensions give equality. Thus the
row spaces of the two reduction matrices are orthogonal complements.
This works even when the characteristic divides $M$.

Both factors in (DL8) are reciprocal up to a sign. Substitution of
$-X$ only rescales rows and columns of the $r$-matrix by signs.
The smaller-period hypothesis and duality finish field independence.
In characteristic two, $X=1+t$ identifies the columns with (DL1).

Every nonzero maximal integer minor now stays nonzero over every
prime field, so has no prime divisor and equals $1$ or $-1$.
The matrix begins with an identity block. Complete any nonzero
smaller minor to a maximal one by adjoining the identity columns of
its complementary rows; none can already occur in a nonzero minor.
This proves total unimodularity. $\square$

## Integral reduced bases and the odd prime limit

The quotient
$$
\mathbb Z[X,Y,Z]/(F_n(X),XY-1,Z-Y^k) \tag{DL9}
$$
is free of rank $n$ over $\mathbb Z$. The constant term makes $X$
a unit. Every Laurent power is a column of $C_{n,M}$ up to sign:
use period $M$ if $n<M$ and signed period $M=n$ if $n$ is a power
of two. Repetitions and signs preserve total unimodularity.

For any standard monomial basis, Cramer's rule expresses each border
monomial with coefficients $0,1,-1$: the denominator basis minor
is a unit and each numerator is another maximal minor. Its normalized
reduced Gröbner basis is therefore monic over $\mathbb Z$, has only
these coefficients, and remains the normalized reduced basis over
every field with unchanged leading ideal. In characteristic zero
(DL7) is a product of distinct dyadic cyclotomic factors, hence
(DL9) is reduced. This flat integral family preserves the whole fan
while its characteristic-two fiber is nonreduced. The conclusion
removes coefficient growth here, not the cost of finding or listing
all monomials of all bases.

The stronger all-field assertion has a definite limit. For any odd
prime $p$, the Pascal matroid of rank $p+1$ on $p^2$ labels has a
$U_{2,p+1}$ minor. After changing the first $p$ rows invertibly,
columns labelled $a+pb$ are $(e_a,b)$. Contract $(e_a,0)$ for
$a\ne0$. The $p$ columns $(e_0,b)$ and one $(e_1,1)$ yield the
$p+1$ distinct lines of a two-dimensional quotient. For $p\ge3$
this contains $U_{2,4}$, impossible over $\mathbb F_2$, whose plane
has only three nonzero vectors. A literal all-field integral lift
preserving these ranks cannot exist. The Fourier construction instead
permits a suitable extension field.

## Universal sampling, simultaneous minima and all equality cases

Here work in characteristic zero, with a primitive $M=p^s$ root.
Call a set $B\subseteq[0,M)$ **balanced** if its counts in the
congruence classes modulo $p^j$ differ by at most one, for every
$0\le j\le s$. The Pascal bases of size $n$ are exactly these sets.
Indeed the ordinary binomial evaluation determinant is the integer
$$
 \det[\tbinom aj]_{0\le j<n,\ a\in B}
 =\frac{\prod_{a<b\in B}(b-a)}{\prod_{j=0}^{n-1}j!}. \tag{DL10}
$$
Its $p$-valuation is a sum, over $j\ge1$, of the excess numbers
of pairs in equal classes modulo $p^j$, relative to the interval
$[0,n)$. At each level the sum of $\binom{c}{2}$ over class counts
is minimized exactly by balanced counts. Every excess is nonnegative;
the determinant is nonzero modulo $p$ exactly when all vanish.

This criterion is established prior art. Osgood--Siripuram--Wu,
Theorem 4 and its Lemmas 1--3, identify the same sets with universal
complex Fourier sampling. Their Schur quotient (13) and valuation
argument (14) give precisely the determinant test above. Thus the
Pascal basis sets are known objects, even though the frequency set
$E_{n,M}$ itself is usually not a universal sampling set.

For a frequency set $J\subseteq[0,M)$ of size $n$, write
$r_J(S)=\operatorname{rank}(\zeta^{ia})_{i\in J,a\in S}$.

**Theorem 3 (simultaneous minimum and equality classification).**
For every $S\subseteq[0,M)$,
$$
 \min_{|J|=n}r_J(S)=R_{n,M}(S). \tag{DL11}
$$
One choice $J=E_{n,M}$ attains every minimum simultaneously.
All frequency sets attaining every minimum have the following exact
recursive description. Write $n=qh+r$, $h=M/p$. If $r=0$, choose
$q$ whole congruence classes modulo $p$. If $r>0$, choose $q$
whole classes, one further class containing $r$ entries, and leave
the others empty. After subtracting that further class's residue
and dividing by $p$, its $r$ entries must satisfy the same rule
at period $h$. Empty and full sets start the recursion.

**Proof of the minimum.** Every balanced $n$-set $B$ is a basis
for every $J$, by the cited universal-sampling theorem. Therefore
$r_J(S)\ge\max_B|B\cap S|=R_{n,M}(S)$. Theorem 1 supplies equality
for $E_{n,M}$. Equivalently, specialize the monic quotient
$\mathbb Z[\zeta][X]/\prod_{i\in J}(X-\zeta^i)$ at the prime
over $p$. It becomes $K[X]/(X-1)^n$, and specialization cannot
increase any power-column rank. This gives the same inequality
without introducing a new uncertainty principle.

By rank-nullity, the largest dimension of a space of signals supported
on a prescribed $S$ and vanishing on $n$ Fourier coordinates is
exactly $|S|-R_{n,M}(S)$. The same $E_{n,M}$ realizes this maximum
for every $S$. This is the corresponding subspace-loss consequence,
not a claim that the scalar uncertainty tradeoff is new.

**Proof of the equality classification.** Suppose the labelled
column matroid is the Pascal matroid. On the $p$ columns
$H=\{0,h,\ldots,(p-1)h\}$ the Fourier rank is the number of
occupied residue classes of $J$ modulo $p$. The Pascal rank there
is $\lceil n/h\rceil$.

The dual Fourier matroid is represented by the rows $-J^c$.
The Pascal matroid is dual to the rank-$(M-n)$ Pascal matroid:
complements of balanced sets are balanced. Thus the number of
occupied classes of $J^c$ is $\lceil(M-n)/h\rceil$.
Let $f,t,e$ be the numbers of full, partial and empty classes.
Then
$$
 f+t=\lceil n/h\rceil,\qquad
 e+t=\lceil(M-n)/h\rceil,\qquad f+t+e=p. \tag{DL12}
$$
If $r=0$ these equations force $t=0,f=q$. Otherwise they force
$t=1,f=q$. This proves the top-level necessity.

For these full classes, perform the Fourier block separation in
Theorem 1. The $q$ full row residues need not be consecutive.
The prime-order Fourier theorem makes their $q$-row matrix, and
that matrix with the partial class's row appended, full spark over
characteristic zero. Hence the rank on any $S$ is
$$
 \sum_u\min(q,l_u)+r_{J'}(\{u:l_u\ge q+1\}), \tag{DL13}
$$
where $J'$ is the scaled partial class. To see necessity of the
smaller equality, take $q$ columns in every fiber and one extra
column exactly in a prescribed set $U$ of fibers. The rank is
$qh+r_{J'}(U)$, while the Pascal rank is $qh+R_{r,h}(U)$.
Thus the smaller matroids must coincide. Conversely (DL13) and
the inductive equality give (DL2) on every subset. This proves
both directions, including $q=0$. $\square$

The number $C_p(n,M)$ of simultaneous minimizing frequency sets is
therefore exactly
$$
 C_p(n,M)=
 \begin{cases}
 1,&n=0\text{ or }n=M,\\
 \binom pq,&r=0,\\
 \binom pq(p-q)C_p(r,h),&0<r<h.
 \end{cases} \tag{DL14}
$$
This counts frequency choices with one common labelled matroid,
not distinct coherent fans. It is not a fan lower bound.

There is an arithmetic interpretation: these are exactly the root
choices whose full power-column matroid is preserved when all roots
coalesce at the prime over $p$. They consequently give the same
unit-power graph fans as the characteristic-$p$ truncation. A
nonminimizing row set gains some bases in characteristic zero;
that alone need not change a particular three-generator fan.

The characteristic-zero scope of Theorem 3 is essential. At period
seven, rows and columns $\{0,1,3\}$ have rank three over $\mathbb C$
but rank two over $\mathbb F_8$, with its primitive root satisfying
$T^3+T+1=0$. The corresponding minor in the power-basis remainder
matrix has absolute norm eight. Thus a clustered set can have a bad
prime other than $p$.
Already at period four the universal row set $\{0,1\}$ has the
power-basis minor $1+i$, of norm two. Universal sampling does not
assert unit minors or the all-residue-prime result for the canonical
construction above.

For another exact priority boundary, maximizing
$-\sum_{a<b\in B}v_p(a-b)$ produces the same balanced sets.
They are a constant-weight Bhargava greedoid. Grinberg's
Theorems 4.7 and 5.1 give their simultaneous Gaussian-elimination
representation over exactly the fields of size at least $p$;
the largest equidistant congruence-tree clique has size $p$.
Thus general representability of this entire flag is established.
The explicit algebraic root choices, their equality classification,
and their specialization details still require their own precise
priority comparison. None resolves the coherent-count problem.

In fact rank minimization does not maximize the coherent count, even
in this class. At $M=8,n=3$ compare rows $\{0,2,4\}$ with
$\{0,1,2\}$ in the graph ideal (DL6). With $k=1$ both have five
initial cones; with $k=2$ they have five and six respectively.
There are only six downsets of size three: the three pure chains
and the origin with two coordinate variables. For $k=1$ the pair
$Y,Z$ is dependent in both models. For $k=2$ the only dependent
candidate in the minimizing model is $\{1,Z,Z^2\}$, whereas all
six are bases in the other model. Every independent candidate can
be made the three cheapest monomials, so these counts are exact.
The checker separately traverses the entire four fans and verifies
their complex matroids by exact Fourier minors. This is a finite
scope guard, not an asymptotic comparison or a fan monotonicity law.

## Literature boundary and verification

Fourier realizations of cyclic representations are established.
[Giansiracusa--Manaker, arXiv:1810.08674v1](https://arxiv.org/pdf/1810.08674),
Section 4.1 before Theorem 4.1.7, describes the selected-row Fourier
matrices and their matroids. Theorem 4.1.7 uses Chebotarev for prime
order; Theorem 4.1.9 gives a necessary minor-vanishing condition at
prime-power order. These passages and proofs were inspected. The
present argument concerns the specified digit-reversed row sets.

[Martin--Reiner, arXiv:math/0402206v1](https://arxiv.org/pdf/math/0402206),
definition of $\mu_n$ and Theorem 1, concern powers in a single
cyclotomic field and simplicial duals. That irreducible quotient differs
from the reducible product (DL7). Fourier transforms, Vandermonde
blocks, cyclic-code duality, regular matroids and monic specialization
are classical.

[Osgood--Siripuram--Wu, arXiv:1204.0992v1](https://arxiv.org/pdf/1204.0992v1),
Theorem 4, its proof and Lemmas 2--3, and Section IV.C's elementary
digit-reversal identity were inspected. Lemma 1's external proof was
not re-audited. The balanced-basis identification above replaces the
earlier incomplete comparison of that theorem with these columns.

[Grinberg, arXiv:2001.05535v5](https://arxiv.org/html/2001.05535v5),
Definitions 1.1, 2.2, 2.5, 2.7 and 4.4, Example 2.4,
Proposition 1.4, and Theorems 4.7 and 5.1
were inspected. The proof of the last two through their valadic
construction was not audited here. The fixed-cardinality matroid
theorem is also in [Grinberg--Petrov, Theorem 7.1](https://arxiv.org/pdf/1909.01965v4).
Its statement and the final deduction from the strong-greedoid
theorem were inspected; the full earlier proof was not audited.

[Delvaux--Van Barel, 2008](https://doi.org/10.1016/j.laa.2008.04.043)
studies prime-power Fourier rank defects and hierarchical blocks.
Only its publisher abstract and metadata were retrieved. The
author-preprint [TW470](https://www.cs.kuleuven.be/publicaties/rapporten/tw/TW470.pdf)
was identified but its full text was unavailable in this audit.
No equivalence with the September 2006 report is assumed. Its
theorems must be compared with (DL11)--(DL14) before any novelty
or publication claim. The complete strongest-current audit remains
open; an unsuccessful access attempt is no evidence of a gap.

In the binary case the root ordering itself is also established:
[Chkifa, arXiv:1503.01731v2](https://arxiv.org/pdf/1503.01731v2),
Section 2, Theorem 2.1 and equations (2.2)--(2.3), describes the
bit-reversed roots as an explicit Leja sequence on the unit disk.
Those definitions and formulas were inspected; their delegated proofs
and the later interpolation estimates were not audited here. This
identifies the roots in (DL5) at $p=2$ with known Leja sections.
It supplies no novelty claim for the root selection, and does not by
itself settle the complete rank or integral-lift comparisons above.

Run from the repository root:

```bash
python3 -B scripts/run_prime_digit_lift_checks.py
./build.sh
python3 tests/test_smoke.py
```

The [report](../results/research/prime_digit_lift_checks.json) records
exact cyclotomic determinant tests, all subset ranks at specified
small periods, cross-characteristic finite Fourier matrices, exhaustive
integer minors, complete positive fans, integral basis specializations,
and odd-prime minor witnesses. It also exhausts every frequency set
and column subset at periods $2,3,4,5,8,9$, checking the simultaneous
rank minima, every equality case, and the recurrence (DL14). The
two bad-prime guards use exact cyclotomic norms and direct
$\mathbb F_8$ arithmetic. The uniform statements rest on the proofs
above. Sharp unrestricted $G_3$, an
infinite polynomial separation, and the full priority audit remain open.
