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

[Osgood--Siripuram--Wu, arXiv:1204.0992v1](https://arxiv.org/pdf/1204.0992),
Theorem 4 characterizes universal Fourier sampling sets by balanced
counts in every prime-power congruence class; Section IV.C uses digit
reversal and identifies the corresponding interval counts. Those
statements, definitions, and the proof of the elementary digit-reversal
identity were inspected. Universal sampling requires every maximal
minor to be nonzero, whereas these frequency sets can deliberately
fail that condition: $E_{5,8}=\{0,1,2,4,6\}$ has parity counts $4,1$.
The paper's full rank and coding comparisons have not been audited.
Precise coverage of (DL3)--(DL9) remains unassessed; these comparisons
establish no novelty claim or publication decision.

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
and odd-prime minor witnesses. The uniform transfer and unimodularity
statements rest on the proofs above. Sharp unrestricted $G_3$, an
infinite polynomial separation, and the full priority audit remain open.
