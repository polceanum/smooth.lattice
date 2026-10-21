# Signed face models and an exact extreme cone

Write
$$
I_{s,r}=(x^2+y^2+z^2+xyz^r)+\mathfrak m^{s+1},\qquad
n=(s+1)^2.
$$
The cheapest-$z$ region is already [uniformly classified](higher_interaction_fans.md).
Here are two exact reductions for the remaining analysis. One supplies all
monomial ranks of a signed face; the other classifies the entire $z^2$ cone
for the fixed degree-five equation. Neither proves a global fan bound for
the original higher-interaction family.

## A signed face, with an exact Laurent model

For $r\ge2$ put
$$
v_y=(1,1-r,1),\qquad
J_{s,r}=(x^2+z^2+xyz^r)+\mathfrak m^{s+1}. \tag{S1}
$$
Highest-weight initial forms give
$$
\operatorname{in}_{v_y}I_{s,r}=J_{s,r}. \tag{S2}
$$
Indeed the three retained terms have $v_y$-degree two, while $y^2$
has smaller degree. The cutoff is homogeneous for this grading.
This gives inclusion of the right side. Both quotients have the
integral monic basis $x<2$, total degree at most $s$, of size $n$;
associated grading of a finite filtered vector space preserves length,
so inclusion is equality. This argument permits the negative entry of
$v_y$ and does not turn a signed weight into a global term order.

**Proposition 1.** Every original monomial-column rank of $J_{s,r}$,
at arbitrary positive prices and over every field, is the rank of
independent finite Laurent blocks as follows.

For signed degree $k=a+c-(r-1)b$, set
$$
b_0=\max(0,\lceil-k/(r-1)\rceil),\quad
b_1=\lfloor(s-k)/r\rfloor,\quad L=b_1-b_0+1,\quad
m=k+(r-1)b_0 .
$$
Omit empty blocks. A parent monomial has the unique form
$$
x^a y^{b_0+j}z^{m+(r-1)j-a},\qquad
0\le j<L,\quad 0\le a\le m+(r-1)j. \tag{S3}
$$
Its column is
$$
(-1)^j t^{a-j}(1+t^2)^j
\quad\hbox{in }\ \mathbb Z[t,t^{-1}]/(1+t^2)^L. \tag{S4}
$$
The block has rank $2L$ if $m\ge1$, and $2L-1$ if $m=0$.
In the latter case its image is the submodule
$\mathbb Z+(1+t^2)\mathbb Z[t,t^{-1}]/(1+t^2)^L$.
Its original price is
$$
k w_z+b_0\beta+\alpha a+\beta j,\qquad
\alpha=w_x-w_z,\quad \beta=w_y+(r-1)w_z>0. \tag{S5}
$$

**Proof.** These are row labels, not a localization of the finite algebra
at $z$. To verify them directly, map all parent columns by (S4).
The relation $x^2+z^2+xyz^r$ maps to zero: its three terms contribute
$t^2+1-t(t+t^{-1})=0$, with a common factor.
If its high-degree term crosses the ordinary cutoff, the two remaining
terms have a factor $(1+t^2)^L$ and also map to zero.

The monic reference columns $a=0,1$ at each $j$ give a basis of the
specified target module. Put $q=1+t^2$. The element $t$ is a unit,
with inverse $-t(1+q+\cdots+q^{L-1})$. At each successive $q$-level,
the two columns have residues $(-1)^j t^{-j}$ and
$(-1)^j t^{1-j}$, a unit basis modulo $q$. If $m=0$, only the
constant column is present at $j=0$; every later level has both.
Thus the indicated reference matrix is integrally unimodular.
Its size is exactly the number of monic reference monomials in this
degree, proving the kernel and all rank assertions. Expanding (S3)
gives (S5). $\square$

Consequently the entire positive fan of $J_{s,r}$ depends on the one
ratio
$$
\rho=\frac{w_x-w_z}{w_y+(r-1)w_z}>
-\frac1{r-1}. \tag{S6}
$$
The [signed-grading theorem](nonnegative_graded_fans.md), Theorem 2a,
already gives $g(J_{s,r})\le12n+7$, uniformly in the grading and field.
This is an application of that previous theorem, not a new unconditional
count for $I_{s,r}$.

## The exact remaining transfer question

For a generic positive weight $w$ with $y$ strictly cheapest, let $B$
be its original standard set. The following criterion would close this
whole region:
$$
v_y\in\overline{C_B(I_{s,r})}. \tag{S7}
$$
Equivalently, for every nonzero coefficient of a reduced border relation
headed by $x^a$, with tail $x^b$,
$$
(a-b)\cdot v_y\ge0. \tag{S8}
$$
If (S7) holds, classical successive-initial theory gives
$\operatorname{in}_w J_{s,r}=\operatorname{in}_w I_{s,r}$.
Conversely, if the original standard set is any standard set of $J_{s,r}$,
its residue classes are a graded basis of the associated quotient.
Its lift is a filtered basis, so every original normal form respects
$v_y$-degree; this proves (S8). Thus this is a specific signed-face
transfer problem, rather than an assumed cover by positive equal-price
planes.

**Conditional consequence.** If (S7) holds at every cutoff, field and
weight with $y$ cheapest, the number of original initials in that
region is at most $12n+7$. The exchanged assertion handles $x$ cheapest.
Together with the proved cheapest-$z$ count, this would give
$g(I_{s,r})=O(n)$ uniformly in $s,r$ and the field, closing this entire
family as a polynomial-amplification route.
Every global monomial initial has a generic positive representative on
the finite ground, so counting the three strict cheapest-variable
regions suffices.

**The premise (S7) is not proved.** The Laurent model describes $J$;
it cannot be substituted for the original ranks without this step.
Finite coverage certificates below support the proposed transfer.
Neither these certificates, the prior failure of a single positive
equality-plane cover, nor a keyword search proves the uniform statement.

## A second integral model for the fixed equation

Now fix $r=3$, $U=x^2+y^2$ and $V=xy$. Let $C(T)\in\mathbb Z[[T]]$
be the classical series
$$
C=1+TC^3,\quad
c_i=[T^i]C=\frac1{2i+1}\binom{3i}{i},\quad
a_i=[T^i]C^2=\frac1{i+1}\binom{3i+1}{i}.
$$
In the binary-adic completion, put $T=-UV^2$ and
$$
G=z^2-UV C(T)^2z+UC(T). \tag{S9}
$$
There is the exact factorization
$$
x^2+y^2+z^2+xyz^3=(Vz+C(T)^{-1})G. \tag{S10}
$$
The prefactor has constant term one, so it is a unit at every finite
ordinary cutoff, over the integers and in every field.

To verify (S10), set $u=-UV^2 C^2$, so $1-u=C^{-1}$ and
$u(1-u)^2=-UV^2$. Then
$G=z^2+(u/V)z+u(u-1)/V^2$; all displayed divisions cancel in
$\mathbb Z[[x,y]]$. Multiplication by $Vz+1-u$ gives $Vz^3+z^2+U$.
In particular this is an explicit special instance of preparation,
with no division in a prime field.

The finite truncation $G_s$ is monic in $z$ of degree two; its other
terms have $z$-degree zero or one and ordinary degree at least two.
Division by it decreases $z$-degree and never decreases ordinary
degree. Therefore
$$
\{x^a y^b z^\epsilon:\epsilon<2,\ a+b+\epsilon\le s\} \tag{S11}
$$
is an integral basis of the original quotient. This model preserves
all original monomial ranks, including unequal $x,y$ prices. It is
a finite exact representation; no faster general rank-query bound
follows from it.

## The whole $z^2$ cone in every characteristic

For a field $K$ and $s\ge2$, define
$$
c=\max\{i:6i+2\le s,\ c_i\ne0\text{ in }K\}.
$$
If $s\ge5$, also define
$$
a=\max\{i:6i+5\le s,\ a_i\ne0\text{ in }K\},\qquad
D=\max(c,2a+1);
$$
for $2\le s<5$ put $D=c=0$.
Both maxima, when required, exist because $c_0=a_0=1$.

**Theorem 2.** The full open positive cone of the original monomial
initial ideal $(z^2)+\mathfrak m^{s+1}$ is exactly
$$
w_z>(2D+1)\max(w_x,w_y)+D\min(w_x,w_y). \tag{S12}
$$
For $s<2$ the ideal itself is the cutoff and there is one initial ideal.

**Proof.** The unique normal form of $z^2$ in (S11) is read from $G_s$.
A term from $UC$ has degree $6i+2$, and its largest binary price is
$$
(4i+2)\max(w_x,w_y)+2i\min(w_x,w_y).
$$
A term from $-UV C^2z$ has degree $6i+5$, and its largest price is
$$
(4i+3)\max(w_x,w_y)+(2i+1)\min(w_x,w_y)+w_z.
$$
Within each binary power $U^{i+1}$ the endpoint coefficient is one;
hence this maximum survives exactly when its scalar $c_i$ or $a_i$
survives. This remains true at equal binary prices, as a whole tier.
The two binary degrees are distinct, so no cross-cancellation occurs.

The first maximum is dominated by its largest surviving index $c$;
the second by $a$. Dividing the first inequality by two writes both
requirements in the increasing family
$(2j+1)\max+j\min$, with indices $c$ and $2a+1$.
This proves (S12) is necessary. It is sufficient because $G_s$ then
has leading term $z^2$. Together with the cutoff monomials its
leading ideal contains $(z^2)+\mathfrak m^{s+1}$, which already has
the correct colength $n$. Thus equality holds. $\square$

This identifies exactly how the [every-prime witnesses](fixed_equation_bad_primes.md)
behave. At $s=6K+5$ in characteristic zero, $D=2K+1$.
At equal binary prices one, (S12) starts at $w_z>6K+4$.
If the characteristic divides $a_K$, then $a\le K-1$ and $c\le K$,
so $D\le2K-1$ and the threshold is at most $6K-2$.
The previous weight $(1,1,6K+1)$ is therefore inside the elementary
$z^2$ cone in every bad characteristic, whereas the rational initial
is different. The same holds for its specified generic perturbation.

In particular, every fixed finite set of primes can simultaneously
give this **same simple initial ideal** at the Chinese-remainder
cutoff. The earlier rank obstruction is not evidence of a fan-count
increase. Characteristic changes here can simplify a whole region.

## Scope, prior art and reproducibility

The cone and successive-initial mechanisms are classical.
[Fukuda--Jensen--Thomas](https://arxiv.org/pdf/math/0509544),
Proposition 2.6, Corollary 2.11, Lemmas 2.12--2.13 and Corollary 2.14,
including their displayed proofs on pp. 7--9, were re-read.
They explicitly permit signed vectors in the cone criterion. They
do not supply (S7) for this family.

The series and coefficient formulas are classical
[Gessel--Xin, Lemma 3.2 and pp. 10--11](https://arxiv.org/pdf/math/0505217),
as already audited in the fixed-equation note. The factorization above
is explicit preparation, not a new general preparation theorem.
The signed-grading bound was proved in the earlier repository note;
its exact priority remains open. The precise specializations here
also have no novelty certificate or independent referee assessment.

~~~bash
python3 -B scripts/run_signed_face_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker verifies the Laurent map against a separate monic model
and full Macaulay relations, and the prepared model against the original
integer kernel. It checks the cone threshold by exact original field
ranks and border coefficients. Small signed-fan interval walks are
lifted to original cones and their entire cheapest-variable patches
are checked by rational vertex inequalities, rather than a weight grid.
These are finite certificates for the stated cutoffs and fields, not
a proof of (S7). The unrestricted $G_3$ objective remains unchanged.
