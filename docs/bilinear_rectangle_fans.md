# Bilinear rectangles: an exact torsion fan and the characteristic-zero boundary

The preserved two-parameter family includes

$$
I_{r,s}=(x^r,y^s,z-x-y-xy),\qquad n=rs.                    \tag{B1}
$$

This note classifies its full global initial fan when $r,s$ are powers
of the characteristic, and gives an exact interpolation model in other
characteristics. The prime-power class has only $O(n)$ initial ideals;
it cannot supply the active goal's polynomial separation. The proof uses
classical translation, Frobenius, binomial coefficients and lattice-ideal
shortest representatives. Exact priority of the restricted count is
unestablished. Neither this classification nor its interpolation dictionary
settles $G_3(n)$ or completes the goal.

## Translation preserves initials, not arbitrary column ranks

A coordinate translation preserves the leading monomial of every
polynomial for every global monomial order: the extra terms from a
translated monomial are proper divisors of it. Applying the inverse
translation gives equality of initial ideals. Independent nonzero
coordinate scalings have the same property.

Put $X=1+x,Y=1+y,Z=1+z$. The translated quotient is

$$
((X-1)^r,(Y-1)^s,Z-XY).                                   \tag{B2}
$$

Its ordinary monomial column at $(a,b,c)$, in the basis
$(X-1)^i(Y-1)^j$, is exactly

$$
P_{i,j}(a,b,c)=\binom{a+c}{i}\binom{b+c}{j},
\qquad 0\le i<r,\quad0\le j<s.                            \tag{B3}
$$

These are integer binomial coefficients reduced to the field. Formula
(B3) gives every selected-column rank of the **translated** quotient.
It does not identify arbitrary selected-column ranks before and after
translation: $x^r$ is zero in (B1), whereas $X^r$ is a unit in (B2).

There is nevertheless an exact rank comparison on every monomial
downset. Expanding $X^aY^bZ^c=(1+x)^a(1+y)^b(1+z)^c$ gives a
unit triangular change of columns within that downset. Thus the two
column families span spaces of the same dimension on every downset,
including every positive-cost prefix of the universal hyperbolic ground.
Their greedy standard sets and entire global initial fans agree.

In characteristic zero, or when $\operatorname{char}K>\max(r-1,s-1)$,
the falling-factorial change of basis is invertible. Hence (B3) is the
evaluation matrix for the polynomial space

$$
\langle U^iV^j:0\le i<r,\ 0\le j<s\rangle
$$

at the integer nodes $(U,V)=(a+c,b+c)$, followed by an invertible row
change. This is a rectangular tensor Vandermonde interpolation problem,
with repeated nodes when the carrier identifies two monomials. It is not
generic interpolation at unrelated nodes. This dictionary by itself gives
no sharp fan count in characteristic zero.

## The all-field torsion quotient

For integers $1\le r\mid s$, over any field define

$$
T_{r,s}=(X^r-1,Y^s-1,Z-XY).                               \tag{B4}
$$

It has length $rs$, with basis $X^\alpha Y^\beta$,
$0\le\alpha<r,0\le\beta<s$. Every monomial has the unit-vector
residue label

$$
(a,b,c)\longmapsto(a+c\bmod r,b+c\bmod s).                \tag{B5}
$$

Distinct labels are independent in every characteristic, including when
the quotient is nonreduced. A standard set chooses the unique cheapest
representative of each label at a strict weight. Periodicity bounds a
candidate by $a<r,b<s,c<s$; for each $c$ the other exponents are forced.

**Theorem.** The number of global monomial initial ideals of (B4) is

$$
g(T_{r,s})=
\begin{cases}
1,&r=s=1,\\
2,&r=1<s,\\
\binom{r+1}{2}+1,&2\le r=s,\\
\binom{r+2}{2},&2\le r<s.
\end{cases}                                               \tag{B6}
$$

This includes composite $r,s$; primality is not needed for this theorem.
The divisibility hypothesis is used. No assertion about nondivisible
rectangles follows from (B6).

### Equal periods

Write a positive weight as $(u,v,w)$ and put $\delta=u+v-w$.
If $\delta<0$, replacing every $Z$ by $XY$ decreases cost. The one
standard set is the $X,Y$ rectangle, with no $Z$.

Suppose $r=s\ge2$ and $\delta>0$. For a residue $(\alpha,\beta)$
put $d=(\beta-\alpha)\bmod r$. Parameterize its representatives by
$0\le b<r$:

$$
a=(b-d)\bmod r,\qquad c=(\beta-b)\bmod r.
$$

Between wrap points their cost has slope $\delta>0$ in $b$.
The wrap of $c$ increases cost, and the only downward jump is the
wrap of $a$ at $b=d$. Therefore only $b=0$ and, if $d>0$, $b=d$
can minimize. The latter cost minus the former is

$$
-ru+d\delta+rw\,\mathbf1_{d>\beta}.
                                                                  \tag{B7}
$$

Set

$$
A=ru/\delta,\qquad C=r(u-w)/\delta.
$$

For $d\le\beta$, choose $b=d$ exactly when $d<A$; for
$d>\beta$, choose it exactly when $d<C$. Count the integers
$d\in\{1,\ldots,r-1\}$ below these thresholds, calling the counts
$k,\ell$. Since $C<A$ and $C=r-rv/\delta<r$,

$$
0\le\ell\le k\le r-1.                                   \tag{B8}
$$

Every pair is realized by $A=k+3/4,C=\ell+1/4$. Explicit strict
integer weights are

$$
(u,v,w)=(4k+3,\ 4r-4\ell-1,\ 4k-4\ell+2).                \tag{B9}
$$

Different pairs give different standard sets. The residue $(0,d)$
selects $Y^d$ exactly when $d<A$, detecting $k$. At residue
$(r-d,0)$, the pure $X^{r-d}$ is replaced by $Y^dZ^{r-d}$
exactly when $d<C$, detecting $\ell$. There are $r(r+1)/2$ pairs.
None is the rectangle: for $\delta>0$ the residue $(1,1)$ chooses
$Z$ instead of $XY$. Adding the rectangle proves the equal-period case.
Boundary weights introduce ties, not additional open initial cones.

### Unequal divisible periods

Let $2\le r<s$. Since $r\mid s$, the relation $Z^r=Y^r$ holds.
If $w<v$, every cheapest representative has $b<r$; otherwise replacing
$Y^r$ by $Z^r$ improves it. Also $\delta>u>0$. The preceding
two-candidate argument still applies, with cost difference

$$
-ru+d\delta+sw\mathbf1_{d>\beta}.                          \tag{B10}
$$

Now use $A=ru/\delta<r$ and $C=(ru-sw)/\delta<A$. Again the
standard sets are indexed by all pairs (B8). All are realized by

$$
\begin{aligned}
u&=s(4k+3),\\
w&=r(4k-4\ell+2),\\
v&=4rs-u+w.
\end{aligned}                                             \tag{B11}
$$

These weights are positive and have $v>w$. The same two residue tests
distinguish the $r(r+1)/2$ sets.

If $w>v$, every cheapest representative has $c<r$. For $\delta<0$
we obtain the rectangle. For $\delta>0$ one has $\delta<u$. With
$c$ increasing, cost decreases between upward wrap jumps. Once
$c>\alpha$, its excess over the cost at $c=\alpha$ is at least
$ru-(c-\alpha)\delta>0$. Thus if $\beta\ge\alpha$, choose
$c=\alpha$. If $\beta<\alpha$, only $c=\beta$ and $c=\alpha$
remain, and their latter-minus-former difference is

$$
sv-(\alpha-\beta)\delta.                                 \tag{B12}
$$

The threshold $D=sv/\delta$ gives exactly $r$ choices, according to
how many of $1,\ldots,r-1$ are below $D$. Each is realized with

$$
(u,v,w)=(4s,\ 2j+1,\ 2s+2j+1),\qquad 0\le j<r.            \tag{B13}
$$

At residue $(d,0)$, the choice between $X^d$ and
$Y^{s-d}Z^d$ distinguishes these $r$ cases. They differ from the
rectangle because $Z$ is selected at residue $(1,1)$.

No set with $w<v$ coincides with one from $w>v$: every latter set
contains $Y^{s-1}$, while every former set has $b<r$. Here
$s\ge2r$ is used. The total is $r(r+1)/2+r+1=\binom{r+2}{2}$.
Finally $r=1$ reduces to choosing the cheaper of $Y,Z$, proving
the remaining cases of (B6).

## Consequences for actual local rectangles

In characteristic $p$, if $r=p^a,s=p^b$, Frobenius makes (B2)
exactly (B4). Interchange $x,y$ if needed to have $r\le s$.
Equations (B6) therefore classify the whole fan of (B1). In particular
for $q=p^a\ge2$,

$$
g(I_{q,q})=\binom{q+1}{2}+1=\Theta(q^2)=\Theta(n).         \tag{B14}
$$

For unequal sides at least two, the count is
$\binom{\min(r,s)+2}{2}$, still $O(n)$. Scaling reduces every
$F=\alpha x+\beta y+\gamma xy$ with $\alpha\beta\gamma\ne0$
to (B1), so this is an entire coefficient class, without genericity.
The Frobenius step is the actual arithmetic input; it gives no assertion
about primes in intervals, prime distribution or fixed-prime unranking.

The same arithmetic is visible directly in (B3): for $i<p^a$,
$\binom{U}{i}$ modulo $p$ depends only on $U\bmod p^a$.
The Pascal matrix on one full period is unit triangular. Tensoring the
two matrices gives (B5) after an invertible row change. This is the
classical Lucas/Frobenius mechanism, not a new prime-digit identity.

There are explicit consequences beyond the count. The weights (B9),
(B11), (B13) and $(1,1,3)$ form a complete coefficient-independent
list. Their heights are $O(rs)$, improving to $O(r)$ for equal sides.
Given a strict weight and one residue label, the proof computes its
standard monomial from at most two candidates after identifying the
region. Listing a standard set takes $O(rs)$ integer arithmetic
operations. Reduced binomial bases follow from its minimal excluded
border; inverse translation gives the local polynomial bases with the
same heads. These are specialized explicit formulas, not a new general
Gröbner algorithm or a bit-complexity claim.

The characteristic-zero count cannot be read from (B14). Independently
verified rational fans have 6, 16 and 32 states for square sides
2, 3 and 4; the corresponding characteristic-prime-power counts are
4, 7 and 11. No asymptotic law is inferred from these three controls.
The remaining task is a uniform count for the structured prefix ranks
in (B3), or a different defined nonuniform family. Mixed-operator image
ranks and support matching still do not give those ranks.

## Primary comparison and reproducibility

[Sabariego--Santos, arXiv:0705.3631v1](https://arxiv.org/html/0705.3631v1),
Definition 2.1, Remark 3.4 and Theorem 3.5, explicitly connect coherent
minimum distance diagrams with graded initial ideals. Their MDD problem
minimizes total path length and resolves ties; (B6) counts all positive
weight regions. The lattice-ideal and shortest-representative dictionary
is established work. Their stated networks are cyclic; (B4) uses a
product group. The general finite-abelian-group framework is also the
subject of [Gómez--Gutiérrez--Ibeas, 2007](https://doi.org/10.1137/050646056),
whose primary abstract was inspected; its full theorem comparison remains
open. The dictionary's existence is not a novelty claim for this note,
and the exact family count has no completed priority certificate.

The tensor interpolation dictionary is also distinct from the generic
tensor-code rank problem. [Brakensiek--Chen--Dhar--Zhang, ECCC TR25-146](https://eccc.weizmann.ac.il/report/2025/146/)
states a conditional deterministic result for maximally recoverable tensor
codes, assuming a Mason question. Its generic/folded-code statements do
not automatically replace the fixed integer-node ranks in (B3). No
coding-theory algorithm, hardness theorem or general matroid-rank
classification is inferred here from that comparison.

```bash
python3 -B scripts/run_bilinear_rectangle_checks.py
```

The checker uses shortest paths and complete positive-cone coverage
independently of the two-candidate proof, ordinary binomial S-pairs,
direct truncated multiplication and translation expansions, Pascal
residues, and rational border bases for the characteristic-zero controls.
The active unrestricted endpoint, earlier results and manuscripts remain
unchanged. The prime-power bilinear class is excluded as a lower-bound
route; the characteristic-zero interpolation count and general mixed
nonlinear rectangles remain open.
