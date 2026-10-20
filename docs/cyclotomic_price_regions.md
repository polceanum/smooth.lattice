# Price regions for mixed circuit products

The distinct-prime construction from the [cographic note](cographic_order_ideals.md)
survives its global small-box screen. In fact its smallest common coordinate
box has exactly $n^2$ monomials. This note proves two regional exclusions
before any full-fan experiment: bounded price ratios give only $O(n^2)$
states for comparable periods, and each sufficiently cheap-coordinate
region has a coefficient-independent list of fewer than $2n$ weights.
The latter holds for general finite algebra towers, including nonreduced
algebras over small fields. Neither result counts the intervening regions.
The sharp unrestricted $G_3$ problem and the paper criterion remain open.

## The algebra, without a characteristic assumption

Let $r,s,t\ge2$ be pairwise coprime integers, put
$S_m(T)=1+T+\cdots+T^{m-1}$, and use

$$
A=K[U,V,W]/(S_r(U),S_s(V),S_t(W)),\qquad
X=UV,\quad Y=UW,\quad Z=VW.                              \tag{P1}
$$

All statements in this note hold over every field $K$. In particular
they do not use separability or infer a small-characteristic rank from a
characteristic-zero computation. The original cyclotomic field is the
special case $K=\mathbb Q$ with three distinct prime periods.

The variables $U,V,W$ are units of periods dividing $r,s,t$. Choose
$j\equiv1\pmod r$, $j\equiv0\pmod s$. Then $X^j=U$, and a
second Chinese remainder exponent recovers $V$. Finally $W=Y/U$.
Thus $X,Y,Z$ generate
$A$, and their relation ideal has colength

$$
n=(r-1)(s-1)(t-1).                                      \tag{P2}
$$

The periodic labels of a monomial are
$(a+b\bmod r,a+c\bmod s,b+c\bmod t)$. The circuit tensor product,
and its cyclotomic/simplicial interpretation over $\mathbb Q$, are
established background, not a new matroid identification; see
[Martin--Reiner, Theorem 1, its proof and Remark 6](https://www-users.cse.umn.edu/~reiner/Papers/cyclotomic.pdf).
Those passages were read as primary text; a full subsequent-literature
audit of the regional assertions below has not been completed.

## An exact obstruction to the common-box shortcut

The subalgebra $F=K[X]=K[U,V]$ has dimension
$D=(r-1)(s-1)$. Indeed choose $j\equiv1\pmod r$,
$j\equiv0\pmod s$. Then $X^j=U$, and a second Chinese remainder
exponent recovers $V$. The monic minimal polynomial of $X$ therefore
has degree $D$, even if $F$ has zero divisors. Also
$A=F[W]/(S_t(W))$ is free of rank $R=t-1$ over $F$.
Permuting the variables gives the three exact pure-axis degrees

$$
D_X=(r-1)(s-1),\quad D_Y=(r-1)(t-1),\quad
D_Z=(s-1)(t-1).                                         \tag{P3}
$$

Every standard set is contained in the box with these three side
lengths: the corresponding monic univariate polynomials have pure-power
heads in every global monomial order. Conversely, all powers
$1,X,\ldots,X^{D_X-1}$ are standard when $X$ is sufficiently cheap,
as the tower theorem below proves. Each axis attains its corresponding
extent in some positive-weight initial ideal. Consequently the smallest
anchored coordinate box containing **all** initial standard sets has
exactly

$$
D_XD_YD_Z=n^2                                           \tag{P4}
$$

monomials. Allowing a translated coordinate box cannot help: every
standard set contains the origin. This establishes that a common
$O(n)$-box exclusion is unavailable for this family, rather than merely
not supplied by its elementary period presentation. A large box alone
is not a lower bound on the number of initial ideals.

## A general cheap-coordinate tower theorem

Let $F=K[x]/(f(x))$, with $f$ monic of degree $D\ge1$; let
$u,v\in F$ be units. Let
$A=F[T]/(g(T))$, where $g\in F[T]$ is monic of degree $R\ge1$.
Map the polynomial variables $X,Y,Z$ to $x,uT,vT$. They generate
$A$: the polynomial representative of $u^{-1}$ recovers $T$ from $Y$.
The quotient length is $n=DR$. No field-extension or radicality
assumption on $F$ or $A$ is made.

**Theorem.** For positive weights $w=(\alpha,\beta,\gamma)$ satisfying

$$
\min(\beta,\gamma)>(D-1)\alpha,                          \tag{P5}
$$

the standard set is determined by the one parameter
$\rho=(\gamma-\beta)/\alpha$. For every $0\le k<R$, take the
greedy $K$-basis of $F$ among

$$
x^a(v/u)^c,\qquad 0\le a<D,\quad0\le c\le k,
\quad\hbox{ordered by }a+\rho c.                         \tag{P6}
$$

Replace each selected slot $(a,c)$ by $X^aY^{k-c}Z^c$ and take the
union over $k$. For $D,R\ge2$ the total number of different standard
sets in (P5) is at most

$$
2(D-1)(R-1)+2<2n.                                      \tag{P7}
$$

**Proof.** Powers with $a\ge D$ reduce using the monic $f$ to
strictly cheaper powers with the same $Y,Z$ factors. For $k<R$, the
degree-$k$ monomials lie in the independent $K$-space $FT^k$ of
dimension $D$. The $D$ monomials $X^a$ times the $k$th power of the
cheaper of $Y,Z$ form a basis for this space. They all cost at most
$k m+(D-1)\alpha<(k+1)m$, where $m=\min(\beta,\gamma)$.
Every higher-$k$ monomial costs at least $(k+1)m$. Thus each block
finishes before any later block can start; expensive remaining terms
of an already completed block are dependent. All $R$ blocks finish
before a monomial of $Y,Z$ degree $R$ can be selected.

Within block $k$, multiplication by the unit $u^k$ preserves every
linear dependence, and the common factor $T^k$ identifies it with $F$.
Removing the common cost $k\beta$ gives exactly (P6).
The only possible changes of slot order occur at
$\rho=j/h$, with $|j|\le D-1$ and $1\le h\le R-1$.
There are at most $2(D-1)(R-1)+1$ such rational numbers, including
zero, hence at most (P7) complementary intervals. Their ordering is
independent of $k$ and of all algebra coefficients. A monomial initial
ideal at a boundary is obtained by a generic perturbation. The block
separation is strict, and a perturbation in each block changes the
same parameter $\rho$, so it adds no further states. $\square$

If $R=1$ there is just the pure $X$ basis. If $D=1<R$ there are at
most two states, using only the cheaper of $Y,Z$. These cases are
included separately when generating witnesses.

There is also an explicit **common** positive-integer weight list.
Sort the distinct rational breakpoints and use the mediant of each
adjacent pair, together with the two extreme samples $\rho=\pm D$.
For a reduced sample $\rho=A/B$, choose

$$
w=(B,\; DB+\max(-A,0),\; DB+\max(A,0)).                 \tag{P8}
$$

Every weight satisfies (P5). When $D,R\ge2$, the denominators of
finite-interval mediants are at most $2(R-1)$, their numerators have
absolute value at most $2(D-1)$, and all three weight coordinates are
at most $2n$. The extreme samples obey the same bound. Thus the
regional list is independent of $f,g,u,v$, and has $O(n)$ members of
$O(n)$ height. This is a counting and coverage result, not a claimed
fast implementation of all the rank computations.

The linear order cannot be improved uniformly. For $D\ge2$ over
$\mathbb Q$, take

$$
f(x)=\prod_{i=0}^{D-1}(x-2^i),\quad g(T)=T^2,\quad
u=1,\quad v=x^D,\qquad n=2D.                            \tag{P8a}
$$

Here $v$ is a unit. In block $k=1$, the $2D$ columns are the powers
$x^e$, $0\le e<2D$. Evaluating at the $D$ roots $2^i$ turns any
$D$ selected columns into an ordinary Vandermonde matrix on distinct
nodes $2^e$. Thus every $D$ columns are independent. A minimum basis
with $j$ slots from row $c=1$ takes its first $j$ slots and the first
$D-j$ slots from row $c=0$. Its cost, up to a constant, is
$j^2-Dj+\rho j$. Successive changes occur exactly at

$$
\rho=D-2j-1,\qquad 0\le j<D.                            \tag{P8b}
$$

All $D+1=n/2+1$ choices occur strictly inside (P5). Hence the
regional $O(n)$ order is sharp along this infinite family. In fact
these are its whole positive fan: $T^2=0$ kills every later block,
and the two remaining blocks are independent regardless of their
interleaving costs. This is a supporting sharpness construction using
classical interpolation, not a separation above the generic scale or
an independent priority certificate.

For (P1), use $F=K[UV]$, $T=W$, $u=U$, $v=V$. The union of the
three cheap-coordinate regions therefore has $O(n)$ initials over
every field. The proof uses monic algebra and finite slot orders;
it is not an application of a general parametric-tree lower bound.

## Comparable periods and bounded price ratios

Put $S=r+s+t$. Every periodic label has a nonnegative representative
whose three coordinates are at most $S$. To prove this, the kernel
of the label map contains the independent lattice vectors

$$
r(1,1,-1),\qquad s(1,-1,1),\qquad t(-1,1,1).             \tag{P9}
$$

Their determinant is $-4rst$. Start from any integral representative,
subtract the real vector $(S/2,S/2,S/2)$, and express the difference
in this lattice basis. Round its three coefficients to nearest
integers and subtract the resulting lattice vector. Each coordinate
of the residual differs from $S/2$ by at most $S/2$, as required.
Existence of an initial integral representative follows from the
Chinese remainder surjectivity already proved, not from assuming (P9) is
the whole kernel; it has index four in the kernel.

Fix a price-ratio bound $C_w\ge1$. If
$\max w_i/\min w_i\le C_w$, a least representative of any label
has each exponent at most $3C_w S$: its cost is at most
$S\sum_i w_i$. Every standard monomial is a least representative
of its own label, since otherwise an equal column occurs earlier.
Hence all standard sets for these weights fit one box of side
$\lfloor3C_w S\rfloor+1$.

If also $\max(r,s,t)/\min(r,s,t)\le C_p$ with fixed $C_p$, this
box has $M=O_{C_p,C_w}(n)$ elements. Restrict the monomial column
matroid to that box. Each relevant standard set is still its unique
minimum basis. Applying the established three-dimensional
lattice-polytope vertex bound to its rank-$(M-n)$ dual profiles gives

$$
g_{\{\max w_i/\min w_i\le C_w\}}(A)
 =O_{C_p,C_w}(n^2).                                     \tag{P10}
$$

This is a regional application of the
[classical common-box/corank bound](cographic_order_ideals.md#scope-and-next-counting-task),
not a new lattice-vertex theorem. A standard set fitting this box in
one region does not imply that all standard sets fit it globally.

## What remains, and what these arguments do not establish

For comparable periods, a polynomial separation above the generic
$n^2\log^2 n$ scale cannot come from bounded price ratios or the
three regions (P5). It would have to use **unbounded price ratios
outside the cheap-coordinate tower regions**. That remaining set is
nonempty: for balanced periods of size $q$, weights proportional
to $(1,q,q)$ have unbounded ratios and lie outside every tower region
for large $q$, because $D_X$ is of order $q^2$. This example specifies
a regime; it supplies no count or lower bound. Unbalanced periods
require their own analysis as well.

The next mathematical question is whether the cyclic/cyclotomic slot
structure controls this intermediate regime uniformly. A useful
extension would account for overlapping, not yet completed tower
blocks. A search for full fans or isolated large examples should wait
for a bound or a mechanism for many **coherent** exchanges there.
Abstract counts of simplicial bases do not supply that mechanism.

The standard-set/rank dictionary and state-polyhedron framework are
established; see
[Babson--Onn--Thomas, Proposition 2.2 and Theorem 2.4](https://arxiv.org/pdf/math/0207135).
The classical FGLM work supplies the general change-of-order linear
algebra setting ([1993 publisher abstract](https://doi.org/10.1006/jsco.1993.1051));
its full text was not inspected in this continuation. These references
do not by themselves establish novelty or noncoverage of (P5)--(P8).
No strongest-current-bound or complete priority audit follows from the
regional reductions. The main goal is active and unachieved.

## Reproducible controls

```bash
python3 -B scripts/run_cyclotomic_price_checks.py
```

The checker compares the tower rule with direct monic polynomial
multiplication and a separate frontier greedy computation. Complete
rational interval coverage is verified for bounded fixtures, with
strict border-relation certificates for the resulting standard sets.
It also tests nonfield, nonreduced towers with a coefficient-dependent
monic extension, pure-axis degrees, periodic-label representatives and
central-region inclusions. The explicit sharpness family is checked by
ordinary monic reduction as well as its independent Vandermonde formula.
Additional price controls outside (P5) test its whole-fan assertion.
It records finite scope, exact witnesses and
source hashes. Computation is supporting evidence for the proofs above,
not a substitute for them or a full intermediate-region fan count.
