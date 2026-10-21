# Equal-degree unit ranks and a sharp quadratic-interaction bound

For
$$
I_{s,r}=(x^2+y^2+z^2+xyz^r)+\mathfrak m^{s+1},\qquad n=(s+1)^2,
$$
the full $r=2$ family has $\Theta(n)$ monomial initial ideals, uniformly
over all fields. Every coherent monomial-prefix quotient of its integral
algebra is torsion-free. Thus $r=2$ has no exceptional characteristic at
any cutoff or monomial order, whereas the previously established fixed
$r=3$ family realizes every prime as an exception at some cutoff.

The main ingredient is a complete-tier rank theorem for arbitrary nonzero
binary forms of the same degree. It is not a classification for $r\ge3$
and does not improve unrestricted $G_3(n)$. Priority of these particular
statements remains unverified.

## A rational interval flag

Put $q=1+t$ and work in $\mathbb Z[t]/t^L$, $L\ge1$. For $0\le i\le N$
and $0\le b<L$ consider
$$
t^b q^{-i},\qquad \text{price }Ai+Bb,\quad A,B>0.
$$
The inverse is its finite integral expansion modulo $t^L$.
All columns with $i,b\ge1$ reduce to cheaper columns, using
$$
t^bq^{-i}=t^{b-1}q^{-(i-1)}-t^{b-1}q^{-i}. \tag{U1}
$$
Only the axes $t^b$ and $q^{-i}$, $i\ge1$, are needed.

**Lemma 1.** If a complete price prefix contains axis intervals
$0\le b\le b_*$ and $1\le i\le i_*$, its rank over every field is
$$
\min(L,b_*+i_*+1). \tag{U2}
$$
The first $L$ independent axis columns give an integral unit determinant.

**Proof.** Multiply by the unit $q^{i_*}$. The negative-power columns
give $1,q,\ldots,q^{i_*-1}$; the other columns are
$q^{i_*},tq^{i_*},\ldots,t^{b_*}q^{i_*}$. Their successive degrees are
$0,1,\ldots,b_*+i_*$, with leading coefficients one. Thus their integral
span is exactly the polynomial interval through that degree. Reduction
modulo $t^L$ gives (U2), in every characteristic. $\square$

## A complete-tier theorem for noncoprime forms

Let $d\ge1$, let $Q,H\in K[x,y]$ be nonzero homogeneous forms of
degree $d$, and set
$$
A_s=K[x,y,z]/\big((1+H)z^d+Q,\mathfrak m^{s+1}\big).
$$
Its length is
$$
n_{s,d}=\sum_{e=0}^{\min(d-1,s)}\binom{s-e+2}{2}. \tag{U3}
$$
Indeed $1+H$ is a unit at the finite cutoff, and division by
$z^d=-Q/(1+H)$ never decreases ordinary degree. The monomials $z<d$
and total degree at most $s$ form a reference basis. All their column
entries are polynomials in the coefficients of $Q,H$.

For $0<u\le Z$ give $x,y,z$ prices $(u,u,Z)$.
For $A\ge0$, $0\le e<d$, $A+e\le s$, put
$$
D=\lfloor(s-A-e)/d\rfloor,\quad m_A=\min(d,A+1),\quad C=T-uA-Ze.
$$
Define $f_D(C)=0$ when $C<0$, and otherwise
$$
f_D(C)=\min\left(D+1,\,
1+\min(D,\lfloor C/(du)\rfloor)
 +\min(D,\lfloor C/(dZ)\rfloor)\right). \tag{U4}
$$

**Theorem 2.** Over every field, without a coprimality hypothesis, the
rank of all original monomials of price at most $T$ is
$$
\boxed{R_{s,d}(u,Z,T)=
 \sum_{e=0}^{\min(d-1,s)}\sum_{A=0}^{s-e}m_A f_D(T-uA-Ze).} \tag{U5}
$$
This takes $O(d(s+1))$ exact arithmetic/floor operations and constant
auxiliary storage. Bit complexity and reading the coefficient input are
separate. This is a theorem about complete tiers, not the individual
monomial bases inside the equally priced binary spaces.

### The coprime upper bound

First suppose $Q,H$ are coprime. The homogeneous complete-intersection
decomposition is
$$
K[x,y]=\bigoplus_\nu t_\nu K[Q,H],\qquad
\sum_\nu t^{\deg t_\nu}=P_d(t)^2,\quad P_d(t)=1+\cdots+t^{d-1}.
$$
Generation by lifts of $K[x,y]/(Q,H)$ and equality of Hilbert series
prove freeness. Changing to this homogeneous basis preserves every
complete price-prefix space because $x,y$ have the same price.

Write $\delta=\deg t_\nu$ and $z$-residue $0\le e<d$.
In $t_\nu Q^aH^b z^{di+e}$ put $q=a+i$. The finite block is
$$
H^b(1+H)^{-i}\pmod {H^{B+1}},\qquad
0\le i\le q,\quad 0\le b\le B,
$$
where $B=\lfloor(s-\delta-e)/d\rfloor-q$.
Its price is
$$
u(\delta+dq)+Ze+du b+d(Z-u)i. \tag{U6}
$$
Lemma 1 applies when $Z>u$. Take the cheapest $B+1$ from the two
axis streams $b=0,\ldots,B$ and $i=1,\ldots,q$.
This proves a coefficient- and characteristic-independent generic
rank law. Polynomial coefficient specialization makes this an upper
bound for any pair of forms, since the fixed $z<d$ reference matrices
specialize without changing the algebra's length.

### A canonical lower bound

Extend the field if necessary. A linear change of the equally priced
binary variables can make both pure-$y$ coefficients nonzero. Scale the
other variable by a parameter and specialize it to zero. Rank cannot
increase. The resulting equation is
$$
z^d+a y^d+b y^d z^d,\qquad ab\ne0.
$$
Diagonal rescaling over an algebraic extension makes $a=b=1$;
the needed $d$th roots exist even in inseparable characteristic. Thus
it suffices to compare the generic law with the canonical curve
$$
y^d+z^d+y^d z^d=0 \tag{U7}
$$
with $x$ free. These are field extensions and coefficient
specializations of finite matrices, not localizations at nilpotents.

Fix $x$-exponent $c$ and residues $j,e<d$ of $y,z$. Put $V=y^d$,
$W=z^d$. The relation is $V+W+VW=0$. In the block of dimension
$D+1$, $D=\lfloor(s-c-j-e)/d\rfloor$, interior products reduce to
cheaper products by $V^aW^i=-V^aW^{i-1}-V^{a-1}W^i$.
Retain $1,V,\ldots,V^D$ and $W,\ldots,W^D$.
Since $W=-V/(1+V)$, consecutive powers of $W$ span the same interval
as consecutive negative powers of $1+V$, together with the constant.
The triangular change has diagonal entries $\pm1$.
Lemma 1 therefore gives the cap of the two axis counts.
Summing over $c,j$ with $c+j=A$ gives multiplicity $m_A$, proving
that the canonical law is (U5).

### Matching the two laws uniformly

Here is the necessary identity; finite agreement is not used as proof.
Take irrational $\theta=(Z-u)/u>0$ and track an item's price and its
first ordinary cutoff as $s$ increases. In each generic block an axis
item is born when its position is at most $B+1$.

For an $i$-axis item, $i\ge1$, its position is
$i+\lfloor\theta i\rfloor+1$. Write $A=\delta+d(q-i)$.
Its price and birth cutoff are
$$
uA+Z(e+di),\qquad
A+e+d\big(i+\lfloor(1+\theta)i\rfloor\big). \tag{U8}
$$
Its multiplicity in binary index $A$ has series
$$
\frac{P_d(t)^2}{1-t^d}=\frac{P_d(t)}{1-t}.
$$
These are exactly the canonical $W^i$ items, including multiplicity
from the free exponent $c$ and residue $j$.

For a generic $b$-axis item let $A=\delta+d(q+b)$.
Its price is $uA+Ze$ and its first-cutoff surplus over $A+e$ is
$$
d\min\big(q,\lfloor b/\theta\rfloor\big).
$$
For height at least $h\ge0$, put $B_h=\lceil\theta h\rceil$.
The multiplicity series in $A$ is
$$
\frac{P_d(t)^2t^{d(h+B_h)}}{(1-t^d)^2}
=\frac{t^{d(h+B_h)}}{(1-t)^2}. \tag{U9}
$$
The canonical $V^a$ item has price $u(c+j+da)+Ze$ and surplus
$d\lfloor a/(1+\theta)\rfloor$. Height at least $h$ requires
$a\ge h+B_h$, and its series is
$$
\frac{P_d(t)}{1-t}\frac{t^{d(h+B_h)}}{1-t^d},
$$
also (U9). Equality of every superlevel multiplicity gives equality
of every exact height. Both axis inventories thus match in price and
first cutoff. No column outside the finite parent bound can precede
a selected axis column: on either stream, its own stream already
contains enough cheaper items to reach the block's cap.

Ranks of finite complete tiers are determined by continuous rank-birth
prices (minima of maxima of finitely many linear prices). Passing to
limits extends the identity to rational ratios and $Z=u$.
The specialization upper and lower bounds coincide, proving Theorem 2.
$\square$

Both forms must be nonzero. The price restriction $Z\ge u$ is part of
the theorem. No extension to the other side, to unequal binary prices,
or to $z^2+Q+Hz^r$ with $r\ne2$ is being inferred.

## The complete quadratic-interaction family

Return to $r=2$ and write
$$
J_s=(x^2+z^2+xyz^2)+\mathfrak m^{s+1},\qquad v=(1,-1,1).
$$
The [signed-face reduction](signed_face_reduction.md) proves
$\operatorname{in}_v I_{s,2}=J_s$ integrally.

**Theorem 3.** For every generic positive weight with $y$ cheapest,
the original standard set equals its standard set for $J_s$.
Consequently $g(I_{s,2})=\Theta((s+1)^2)=\Theta(n)$ over every field.

**Proof.** The middle-$z$ region is already
[proved](middle_interaction_fans.md), so suppose $Z>X>Y>0$.
Put
$$
t=(X-Y)/2,\quad u=(X+Y)/2,\quad Z_0=Z-t.
$$
Then
$$
(X,Y,Z)=(u,u,Z_0)+tv,\qquad Z_0-u=Z-X>0. \tag{U10}
$$
Let $B$ be the minimum basis of $J_s$ at this generic weight.
Grading makes it a minimum basis at the seed $(u,u,Z_0)$ as well.
It lifts to a filtered original basis. Every nonzero original border
coefficient has difference $\eta$ satisfying $\eta\cdot v\ge0$.

Theorem 2 applies both to $Q=x^2+y^2,H=xy$ and to
$Q=x^2,H=xy$. The pairs have the same seed price-rank inventory.
Hence $B$ is also a minimum original basis at the seed, and every
original border difference has $\eta\cdot(u,u,Z_0)\ge0$.
Equation (U10) proves its nonnegativity at the target weight.
It is strict when $\eta\cdot v>0$; otherwise the coefficient is the
graded face coefficient and genericity gives strictness. Thus $B$ is
the original standard set throughout this remaining region.

The earlier signed-grading theorem bounds each cheapest-$y$ or
cheapest-$x$ region by $12n+7$ cones. The cheapest-$z$ region has the
previous exact count $2M_{s,2}+4=O(n)$, for $s\ge4$.
Every monomial initial ideal has a generic positive representative,
which can avoid all price-equality planes. These three regions cover
the count. In particular
$$
g(I_{s,2})\le24n+2M_{s,2}+18\quad(s\ge4). \tag{U11}
$$

For the lower bound, the existing cheapest-$z$ wall formula includes
every reduced $p/q\in(0,1)$ with $q\le(s+2)/6$. There are
$\Omega(s^2)$ such fractions: among $1\le p,q\le Q$, pairs with a
common divisor are at most
$Q^2\sum_{m\ge2}m^{-2}\le3Q^2/4$; halve the coprime off-diagonal
pairs. Thus $M_{s,2}=\Omega(s^2)$ uniformly in the field.
The finitely many smaller cutoffs do not affect the asymptotic claim.
$\square$

This excludes the whole $r=2$ family as a polynomial-amplification
route; it does not achieve the unrestricted extremal goal.

## Integral consequences and a complete rank query

**Corollary 4.** For the coefficient-one $I_{s,2}$, every monomial
initial ideal is independent of the field. Every positive-weight
complete-prefix sublattice has a torsion-free quotient in the integral
algebra. Every reduced rational Gröbner basis has integer coefficients,
and their finite union contains a monic Gröbner basis for every order
over $\mathbb Z$, valid after any base change. In particular its
universal denominator is one.

**Proof.** At $Z>u$, Theorem 2 makes the whole-tier ranks of the integral
face $J_s$ equal over $\mathbb Q$ and every $\mathbb F_p$.
Its independent signed-degree blocks can only lose rank on reduction.
Equality of their sum forces equality in each block at every prefix.
Within a block the seed ratio can be chosen irrational; then all
different columns have distinct prices. Thus its generic minimum
monomial basis is the same in every field. Its reference determinant
is an integer nonzero modulo every prime, so is $\pm1$.
The positive-ratio part through ratio one has unit bases by the earlier
Laurent interval lemma; the outer monic $x<2$ case has a unit basis too.
These ranges cover the entire positive
fan of $J_s$, including limits of complete tiers.

Filtered lifting preserves the unit determinant. Theorem 3 and the
previous characteristic-free cheapest-$z$ theorem therefore give an
integral unit basis for every original cone. Normal forms have integer
coefficients and use only no-more-expensive basis columns. Thus every
prefix spans exactly a subset of this unit basis; complete tiers
follow by refinement. Its quotient is free abelian. Border normal forms have integer
coefficients. The integral algebra already has a monic reference
basis, so their integer membership is actual membership in the given
integral ideal, not just in its rational saturation. The border basis
gives a monic Gröbner basis; retaining minimal heads gives the reduced
basis. Finiteness of the fan gives a finite union. The base-change
and denominator statements are classical consequences. $\square$

There is also an exact $O(s+1)$ arithmetic/floor rank query at every
positive original weight, with constant auxiliary storage. The earlier
cheapest-$z$ and middle-$z$ queries cover those regions. For $Z>X\ge Y$
use the signed-face blocks, with $r=2$, from (S3)--(S5). Write
$$
\alpha=X-Z<0,\quad \beta=Y+Z,\quad A=-\alpha>0,
\quad C=\alpha+\beta=X+Y>0.
$$
After $t\mapsto t^{-1}$, reflect $a'=m+j-a$ and split $a'=2h+e$,
$e=0,1$. Put $w=1+t^{-2}$. Up to a common integral unit, the columns
are $w^j(1-w)^{-h}$. For each parity set
$$
j_0=\max(0,e-m),\quad L_e=L-j_0,\quad M=m+j_0-e.
$$
Omit $L_e\le0$. Factoring $w^{j_0}$ leaves
$$
0\le j'<L_e,\qquad 0\le h\le\lfloor(M+j')/2\rfloor. \tag{U12}
$$
Their relative prices are $2Ah+Cj'$, with constant offset
$\alpha m+Cj_0+Ae$ in the original block.

Retain the axis $h=0$ and, for each $h\ge1$, only
$j_h=\max(0,2h-M)<L_e$. Every other column reduces to cheaper columns
by
$$
w^j(1-w)^{-h}
=w^{j-1}(1-w)^{-h}-w^{j-1}(1-w)^{-(h-1)}.
$$
Whenever a boundary column arrives, the axis already reaches its
$j_h$. The boundary columns arrive in consecutive $h$ order. Their
pole at $w=1$ has leading coefficient one. Inductively the span is
exactly all polynomials through degree $J+H$, divided by $(1-w)^H$,
where $J$ is the last cheap axis index and $H$ the last cheap boundary
index. Indeed, on adding a pole, the old numerators are divisible by
$1-w$, and the new numerator $w^{j_h}$ has value one at $w=1$.
On adding an axis index, its leading coefficient is one. Reduction
modulo $w^{L_e}$ therefore gives
$$
\min(L_e,\#\text{ cheap axis items}+\#\text{ cheap boundary items}).
\tag{U13}
$$
The two boundary pieces have prices $2Ah$ when $2h\le M$, and
$2\beta h-CM$ otherwise. Each count is a clipped floor calculation.
There are only $O(s+1)$ signed blocks. The parity with $m=0,e=1$
starts at $j_0=1$ and has cap $L-1$, exactly the missing residue in
the original Laurent model. No localization of the Artinian algebra
or dense multiplication matrix is used.

## Scope, prior art and reproducibility

The homogeneous free decomposition, specialization of matrices,
rational interval flags, signed cone lifting and monic base change
are classical ingredients. The new claims here are their particular
uniform rank identity and application, with priority still open.

[Abbott--Bigatti--Robbiano, arXiv:1801.06112v3](https://arxiv.org/pdf/1801.06112v3),
Theorem 2.2 and its proof, shows why an integral generator set and an
integral reduced basis alone need not define the same modular ideal:
membership coefficients matter. The monic reference quotient above
supplies that membership. Their Section 2.2 defines the universal
denominator and proves finite exceptional sets for each fixed ideal;
it does not assert denominator one for this family. This is a scoped
comparison, not a novelty certificate.

The unrestricted $G_3$ objective, the full $r\ge3$ family count and its
strongest-current priority audit remain open. After the prior uniform
bound and the present exclusion, this family's remaining power-growth
candidate is $3\le r<\sqrt{s}$, in the $z$-largest region only.
The $r=3$ bad-prime witnesses remain valid and cannot be counted as
extra cones.

~~~bash
python3 -B scripts/run_equal_degree_unit_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker independently uses polynomial quotient columns and field
Gaussian elimination, compares the two birth inventories, checks
integral flag determinants and Macaulay kernels, verifies full original
rank queries, and certifies entire finite original regions by integral
border relations. Computation supports the uniform proofs above; it
does not replace them or establish novelty.
