# A fixed equation with exceptional primes at arbitrarily deep cutoffs

Consider the single polynomial
$$
F=x^2+y^2+z^2+xyz^3,\qquad
A_s=\mathbb Z[x,y,z]/(F,\mathfrak m^{s+1}). \tag{P1}
$$
The algebra is free over $\mathbb Z$, of rank $(s+1)^2$. Nevertheless,
every prime changes a coherent standard set at a suitable cutoff. The
change occurs in the remaining expensive-$z$ region of the
[higher-interaction analysis](higher_interaction_fans.md), even though
the equation, its degree and all its coefficients remain fixed.

This is an arithmetic obstruction to transferring rational fan computations
to one fixed finite collection of characteristics uniformly in $s$. It
does **not** prove an increase in the number of initial ideals, change the
unrestricted $G_3$ bound, establish novelty, or complete the active goal.

## An integral reduction preserving equal binary prices

Put $U=x^2+y^2$ and $V=xy$. There is a graded integral decomposition
$$
\mathbb Z[x,y]=
\mathbb Z[U,V]\oplus x\mathbb Z[U,V]\oplus
y\mathbb Z[U,V]\oplus y^2\mathbb Z[U,V]. \tag{P2}
$$
One proof uses the monic relations
$x^2+y^2-U$, $xy-V$, $y^3-Uy+Vx$. With degree order on $x,y$
and $x>y$, over the coefficient ring $\mathbb Z[U,V]$, their leading
monomials are $x^2,xy,y^3$; all S-polynomials reduce to zero.
Thus the four displayed remainders form a free basis over the coefficient
ring. Their degrees are $\delta=0,1,1,2$.

For $F_r=x^2+y^2+z^2+xyz^r$, $r\ge3$, eliminate
$U=-z^2-Vz^r$ in each summand. The image of the ordinary-degree cutoff
in its degree-$\delta$ summand is precisely
$$
(z^aV^b:a+2b>s-\delta).
$$
Indeed replacing $U$ never decreases ordinary degree, and the terms
with no $U$ already generate the displayed image. Therefore the quotient
is the direct sum, as an abelian group, of the four truncated modules
in $z,V$. This also reproves freeness and rank $(s+1)^2$.

At weights $(t,t,1)$, replacing monomials of any fixed binary degree
by their basis in (P2) is an integral unimodular change within one
price tier. It preserves the span of every complete tier. It is not
a price-preserving operation when the prices of $x,y$ differ.

A parent column $z^aU^iV^j$ in the summand of degree $\delta$ has label
$$
q=a+2i-(r-2)j.
$$
Its possible output rows are $z^{q+(r-2)b}V^b$, with
$$
b_0=\max(0,\lceil-q/(r-2)\rceil),\quad
b_1=\lfloor(s-\delta-q)/r\rfloor .
$$
Different pairs $(\delta,q)$ are independent blocks; the two degree-one
summands keep their separate labels. Write $L=b_1-b_0+1$ and omit empty
blocks. After identifying those rows with $1,\xi,\ldots,\xi^{L-1}$,
the columns, up to their signs, are
$$
\xi^J(1+\xi)^i\pmod{\xi^L},\qquad
0\le J<L,\quad
0\le i\le\left\lfloor\frac{q'+(r-2)J}{2}\right\rfloor,\quad
q'=q+(r-2)b_0. \tag{P3}
$$
Their original prices are
$$
q'+2tb_0+\delta t-2(1-t)i+(r-2+2t)J. \tag{P4}
$$
The notation $\xi$ labels rows; this construction does not invert the
nilpotent variable $z$. Increasing upper polynomial degrees in (P3),
unlike the previously solved region, permit nonunit prefix determinants.

## One explicit cyclic quotient

**Theorem 1.** Let $K\ge1$, and set
$$
s=6K+5,\quad H=6K+1,\quad w=(1,1,H),\quad
a_K=\frac1{K+1}\binom{3K+1}{K}. \tag{P5}
$$
In the quotient over any field, $z^2$ increases the rank of the monomials
strictly cheaper than itself if and only if the characteristic does not
divide $a_K$. In characteristic zero the rank increases. The conclusion
holds for every monomial-order refinement of $w$.

More precisely, the block $(\delta,q)=(0,2)$ has dimension $L=2K+2$.
Its cheap columns span a primitive sublattice of rank $L-1$. After adjoining
$z^2$, the block quotient is $\mathbb Z/a_K\mathbb Z$.

**Proof.** Normalize $w$ to $(t,t,1)$ with $t=1/H$. In this block $b_0=0$,
$b_1=2K+1$ and $i\le1+\lfloor j/2\rfloor$. At fixed $j$, the cheapest
column uses the largest $i$. For $j=2v$ its price is $(6v+2)/H$;
for $j=2v+1$ it is $1+(6v+4)/H$.
The columns below price two are exactly these largest-$i$ columns
for $j=0,\ldots,2K$. Decreasing $i$ puts every further column above
price two, except $j=i=0$, which is $z^2$ itself at price two.
In particular no other parent column in this block has that price.

Divide all these columns by the common unit $1+\xi$. With
$T=\xi^2(1+\xi)$, they become
$$
P_{2v}=T^v,\qquad P_{2v+1}=\xi T^v,\qquad 0\le j\le2K. \tag{P6}
$$
Each $P_j$ starts with $\xi^j$ with coefficient one. Hence they give a
primitive rank-$L-1$ lattice. The image of $z^2$ after division is
$W=(1+\xi)^{-1}$.

Let $u\in T\mathbb Z[[T]]$ solve $T=u(1-u)^2$ and set
$C=(1-u)^{-1}$. Then $C=1+TC^3$ and
$$
\frac1{1+\xi}=C(T)-\xi C(T)^2. \tag{P7}
$$
For completeness, $\rho=u-1$ is a root of $\rho^3+\rho^2-T$.
Factoring that polynomial and substituting $T=\xi^2(1+\xi)$ gives
$\xi^2+u\xi+u(u-1)=0$, since $\xi-u+1$ is a unit.
Multiplying the right side of (P7) by $1+\xi$ now gives one.

The expansion in the triangular basis (P6) has coefficient
$-[T^K]C(T)^2=-a_K$ on $P_{2K+1}$. The coefficient formula is classical
Lagrange inversion. It follows over the integers from the recurrence
$C=1+TC^3$; reduction in characteristic dividing $K+1$ does not entail
division by $K+1$ in that field. Eliminating $P_0,\ldots,P_{2K}$ from
$W$ leaves $-a_K\xi^{2K+1}$. The square matrix with these columns
has Smith form
$$
\operatorname{diag}(1,\ldots,1,a_K). \tag{P8}
$$
There is no coupling to another block. Other original monomials tied
with $z^2$ have no further component in this block, because the integral
binary change of basis preserves complete tiers and $z^2$ is its sole
parent column at that price. Thus rational independence persists even
after all other tied columns, whereas divisibility gives dependence
already on strictly cheaper columns. This proves every refinement claim.
$\square$

The whole algebra in (P1) has no integer torsion. The torsion in (P8)
belongs to a quotient by a specified monomial-prefix span, and is not
a statement about local cohomology or the whole prefix's exact Smith form.

There is also an explicit integral membership certificate. If
$c_j$ is the coefficient of $P_j$ in (P7), put
$i_j=1+\lfloor j/2\rfloor$ and $e_j=j\bmod2$. Then
$$
z^2-\sum_{j=0}^{2K}c_j(-1)^{i_j}z^{e_j}U^{i_j}V^j
+a_K z^{2K+3}V^{2K+1}\ \in\ (F,\mathfrak m^{s+1}). \tag{P9}
$$
The middle terms are combinations of strictly cheaper original monomials.
The last monomial has degree exactly $s$. The checker expands (P9) in
$x,y,z$ and replays an integer monic-division witness independently.

## Every prime, with a quadratic length budget

**Corollary 2.** For every prime $p$, a truncation of the single equation
(P1) has a positive generic weight where the characteristic-$p$ initial
standard set differs from the rational one. A cutoff of size $O(p)$
and length $O(p^2)$ suffice.

**Proof.** For $p\ge7$ choose $K=\lceil(p-1)/3\rceil$. Then
$$
2K+2\le p\le3K+1<2p,\qquad K+1<p.
$$
The expression
$$
a_K=\frac{(2K+2)(2K+3)\cdots(3K+1)}{(K+1)K!}
$$
has exactly one numerator factor divisible by $p$ and none in the
denominator; thus $v_p(a_K)=1$. For $p=2$ take $K=1$, $a_1=2$;
for $p=3,5$ take $K=3$, $a_3=30$.
For $p\ge7$, $s\le2p+7$ and $n\le(2p+8)^2$.

To remove ties explicitly, put $B=s+1$, $M=B^3$ and take
$$
\widetilde w=(M+1,\ M+B,\ MH+B^2). \tag{P10}
$$
On the ground $|a|\le s$, the perturbation $a_x+B a_y+B^2a_z$
has distinct integer values between zero and $M-1$. It preserves every
strict original integer-price comparison and distinguishes every pair.
An original tie with $z^2$, other than $z^2$, must have $a_z=1$ and
$a_x+a_y=H$: the alternative $a_z=0$ needs degree $2H>s$.
Its perturbation is at most $BH+B^2<2B^2$, since $H=B-5$.
All such tied monomials therefore precede $z^2$. Theorem 1 applies.
Zero columns beyond the cutoff do not affect the resulting monomial
initial ideal. $\square$

For each fixed cutoff, the set of exceptional primes across all orders
is finite, as in classical universal-denominator theory. The corollary
concerns the union over **different** cutoffs of the same equation.
It does not assert one finite quotient has infinitely many bad primes,
or that a prime remains bad at every larger cutoff.

## Simultaneous failure of any fixed finite set of characteristics

**Corollary 3.** Given any nonempty finite set of primes $S$, there is
one $K$ and one generic weight (P10) for which the monomial prefix through
$z^2$ has strictly smaller rank in **every** characteristic in $S$ than
over $\mathbb Q$.

**Proof.** For $p\ge7$ let $k_p$ be the witness in Corollary 2.
If $K\equiv k_p\pmod p$, then $p\nmid K+1$,
and adding $K$ to $2K+1$ has a carry in the lowest base-$p$ digit:
$2k_p+1<p\le3k_p+1<2p$. Thus $p\mid a_K$.

Here the elementary factorial-valuation identity suffices:
$$
v_p(a_K)=
\sum_{e\ge1}\left(
\left\lfloor\frac{3K+1}{p^e}\right\rfloor
-\left\lfloor\frac K{p^e}\right\rfloor
-\left\lfloor\frac{2K+1}{p^e}\right\rfloor
\right)-v_p(K+1). \tag{P11}
$$
Every summand in parentheses is nonnegative. The $e=1$ summand
is one in the stated residue class.
For $p=5$, use $K\equiv3\pmod5$: the two low digits being added
are three and two, so the same conclusion holds.
For $p=2$, impose $K\equiv1\pmod4$. The $e=1,2$ summands are
both one and $v_2(K+1)=1$, leaving a positive valuation.
For $p=3$, impose $K\equiv3\pmod9$; the $e=2$ summand is one
and $3\nmid K+1$. All these implications hold independently of
the higher digits of $K$.

The Chinese remainder theorem supplies one positive $K$ satisfying
these pairwise coprime congruences for all $p\in S$. It may be chosen
smaller than
$$
M_S=\prod_{p\in S}p^{e_p},\qquad
e_2=e_3=2,\quad e_p=1\ (p\ge5). \tag{P12}
$$
The explicit block loses one rank in every chosen characteristic.
Every other integral block has field rank at most its rational rank,
so no block can compensate for that loss. The generic prefix contains
all original tiers through the price of $z^2$, by the construction
following (P10), and the same conclusion holds there. $\square$

For example $K=93$ works simultaneously for $S=\{2,3,5,7\}$.
The length is less than $(6M_S+6)^2$; this simultaneous bound is
multiplicative in the sampled primes, rather than the single-prime
quadratic budget. Consequently even taking the maximum of ranks over
one fixed finite list of characteristics can undercount the rational
rank on this family. This does not rule out adaptive or randomized
modular methods with verification, or symbolic integer calculations.
The residue and Chinese-remainder arguments are classical; they are
not new number-theoretic principles.

## Primary-source boundary

The coefficient sequence and power-series law are established mathematics.
[Gessel--Xin, arXiv:math/0505217](https://arxiv.org/pdf/math/0505217),
Lemma 3.2, equations (3.3)--(3.5) and the proof on pp. 8--9, and the
coefficient definitions/table on pp. 10--11 were inspected. Their
$b_K$ is exactly $a_K$ above. No new tree enumeration, Lagrange formula,
or general determinant principle is asserted.

[Abbott--Bigatti--Robbiano, arXiv:1801.06112v3](https://arxiv.org/pdf/1801.06112v3),
Section 2.2, Remark 2.3, Definition 2.3 and Proposition 2.1 with its
proof were re-read. A fixed ideal has a finite universal denominator.
Theorems 2.1--2.3 and the order-sensitive bad-prime framework were
previously inspected; Theorem 2.1 delegates its proof to a further source.
Those results do not give a bound uniform over the different $s$ in (P1).
Their examples already exhibit very many bad primes; the phenomenon
itself is not new.

Unbounded prime torsion in a fixed hypersurface also has much older
precedent. [Singh, arXiv:math/0406354](https://arxiv.org/pdf/math/0406354),
Section 4 and its full displayed argument on pp. 8--9 were inspected.
It proves every-prime torsion in a local cohomology module of a
six-variable hypersurface, using a grading and a direct limit. Here the
module is a finite price-prefix quotient of a free three-variable
truncation. No identification with that local cohomology construction
or new local cohomology consequence has been established.

These exact differences identify comparisons, not a novelty certificate.
The specific fixed-equation realization and its quantitative cutoff still
require comparison with the strongest current work before a publication
claim. A keyword search cannot settle that question.

## Reproduction and consequence for the counting goal

~~~bash
python3 -B scripts/run_fixed_equation_prime_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The saved [report](../results/research/fixed_equation_prime_checks.json)
checks two integral models against full Macaulay relations, the signed
block columns, the cyclic matrix, an independent series recurrence,
original monomial ranks over several fields, integral membership
certificates and the generic perturbation. Finite checks support the
uniform proof; they do not establish priority.

For the main problem this eliminates a shortcut: no fixed finite list of
characteristics can represent all these rational prefix ranks uniformly
in the cutoff. It does not exclude symbolic integer methods, randomized
modular methods with verification, or the existing field-independent
upper bounds. The next count must control simultaneous coherent exchanges
in the expensive-$z$ region, retaining their arithmetic. Merely counting
the primes in this corollary would not count cones of any one field's fan.
