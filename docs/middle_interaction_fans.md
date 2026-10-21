# Noncoprime binary ranks and the entire middle-price region

The remaining transfer question in the [signed-face note](signed_face_reduction.md)
is now proved whenever the interaction variable $z$ is between the other
two prices. The proof uses a whole-tier rank law for arbitrary nonzero
binary forms, including forms with a common factor. The region where $z$
is strictly largest is still unresolved. The unrestricted $G_3$ bound
and the publication assessment are unchanged.

## A Laurent interval lemma

Let $q=1+t^2$, $d,L\ge1$, $m\ge0$, and consider the columns
$$
u_{a,j}=t^{a-j}q^j\pmod {q^L},\qquad
0\le j<L,\quad 0\le a\le m+dj,
$$
of price $\rho a+j$, where $0<\rho<1$.
The target is $\mathbb Z[t,t^{-1}]/q^L$ if $m>0$, and its submodule
$\mathbb Z+q\mathbb Z[t,t^{-1}]/q^L$ if $m=0$.

**Lemma 1.** Retain only these parent columns:
$$
\begin{array}{ll}
j=0,&0\le a\le m,\\
j\ge1,&a=0,\\
j\ge1,&\max(1,m+d(j-1))\le a\le m+dj .
\end{array} \tag{M1}
$$
Every price-prefix rank, over every field, is
$$
\min\bigl(2L-\mathbf1_{m=0},
          \#\{\text{retained columns below the price cutoff}\}\bigr). \tag{M2}
$$
The same holds with complete tiers at $0\le\rho\le1$.
Every generic minimum basis in this model is integral and unimodular.

**Proof.** An omitted column has $j\ge1$ and
$1\le a\le m+d(j-1)-1$. The identity
$$
u_{a,j}=u_{a+1,j-1}+u_{a-1,j-1} \tag{M3}
$$
expresses it in strictly cheaper parent columns. Recursing in $j$
proves equality of the filtered spans.

First suppose $m>0$. The first line of (M1) supplies consecutive
Laurent exponents $0,\ldots,m$. The third line, in its price order,
extends the largest exponent successively through every integer above
$m$; its new endpoint is $a+j$. The second line extends the smallest
exponent through $-1,-2,\ldots$.

A right-hand column at level $j$ has price greater than $j$, so the
left-hand column at that level is already present. Conversely, before
the left-hand column at level $j$ arrives, every positive endpoint
through $j$ is present. An endpoint in the spine costs $\rho j<j$.
For a right-hand endpoint $n=j$, its level $j'<j$ and price
$\rho(j-j')+j'<j$. Thus every new column has coefficient one at its
one new extreme exponent, and all other exponents are in the previous
interval. Its integral span is precisely a consecutive Laurent
polynomial space. Reduction modulo the monic $q^L$, whose constant
coefficient is one, gives rank the smaller of its dimension and $2L$.

If $m=0$, exponent one is missing from the positive endpoint list.
The first left column is $t+t^{-1}$, preceded only by the constant.
After it arrives, the same induction gives the intersection of a
consecutive Laurent space with the condition that its residue modulo
$q$ be constant. This is one primitive integer equation. A space of
width $W$ has dimension $W-1$; its kernel modulo $q^L$ has dimension
$\max(W-2L,0)$. This proves the cap $2L-1$. All pivots in the interval
argument have coefficient one, so no characteristic is excluded.

These arguments apply at generic $\rho$. Rank-birth prices on a finite
ground are continuous minima of maxima of finitely many linear prices.
Passing to limits proves the complete-tier statements, including the
endpoints. $\square$

## The positive whole-tier theorem without a coprimality assumption

Let $Q(Y,Z)$ and $H(Y,Z)$ be nonzero homogeneous binary forms of degrees
$2$ and $r+1$, with $r\ge1$, over any field. Put
$$
A_s=K[X,Y,Z]/(X^2+Q+XH,\mathfrak m^{s+1}).
$$

**Theorem 2.** At prices $(1+\lambda,1,1)$, $\lambda\ge0$, every
complete price-tier rank is given by the earlier universal formula
(H18)--(H19) in [the higher-interaction note](higher_interaction_fans.md).
The forms need not be coprime.

Explicitly, for $0\le\lambda<r$, use channels $\delta=0,\ldots,r+1$
with multiplicities one at the endpoints and two otherwise. In each
channel and for $0\le q\le s-\delta$, take the
$$
L=\min(q,\lfloor(s-\delta-q)/r\rfloor)+1
$$
smallest values from
$$
\kappa_{q,i}=
\max\bigl((1+\lambda)q-2\lambda i,\ (1-r)q+2ri\bigr),
\qquad 0\le i\le q. \tag{M4}
$$
Their prices are $\delta+\kappa_{q,i}$. Count how many are below the
requested cutoff, with multiplicity. For $\lambda\ge r$ the monic
$X$-exponent-below-two basis supplies the complete-tier formula.
This retains the earlier $O((r+1)(s+1))$ exact floor-operation bound
and constant auxiliary storage; coefficient input and bit costs remain
separate.

This is a law for complete tiers. It does not identify the individual
monomial bases at equal binary prices. The extension to arbitrary
nonzero forms is not asserted when $-1<\lambda<0$.

### A degeneration sandwich

A generic pair of forms is coprime. The previously proved theorem gives
its rank by (M4), over every characteristic. In the fixed monic $X$
reference basis, all column entries are polynomials in the form
coefficients. Specialization cannot increase a rank, so (M4) is an
upper bound for every pair.

For the reverse bound, extend the field to an infinite algebraic
extension if necessary; ranks do not change. A linear change of $Y,Z$
can make both pure-$Y$ coefficients nonzero. Such a change preserves
ordinary cutoffs and every whole price-prefix space, since the two
variables have the same price. Replace $Z$ by $\tau Z$ and set
$\tau=0$. The resulting equation is
$$
X^2+aY^2+bXY^{r+1},\qquad ab\ne0.
$$
Its monic reference algebra still has length $(s+1)^2$, and its matrix
is a specialization of the changed original matrix. Rescale $X,Y$
over a further algebraic extension to make $a=b=1$.
Thus the original rank is at least the corresponding rank for
$$
X^2+Y^2+XY^{r+1}, \tag{M5}
$$
with $Z$ a free variable subject to the total cutoff. The required
roots exist even when their degrees are divisible by the characteristic;
only rank preservation under field extension is used.

It remains to prove that this canonical lower bound is (M4).

### The canonical curve and its rank births

For a fixed exponent $c$ of the free variable, split the total
$X,Y$ degree into residues $j$ modulo $r$:
$$
a+b=j+r\ell,\quad 0\le j<r,\quad
L=1+\lfloor(s-c-j)/r\rfloor.
$$
The column of $X^aY^bZ^c$ is
$(-1)^\ell t^{a-\ell}(1+t^2)^\ell$ in a Laurent block.
The reference columns have $a=0,1$, except that $a=1$ is absent at
$j=\ell=0$. Their integral unit matrix verifies this map directly,
without localizing the nilpotent variable $Y$.

Its dimension is $2L-\mathbf1_{j=0}$, its parent bound is
$a\le j+r\ell$, and its price is $c+j+r\ell+\lambda a$.
Lemma 1 applies with $m=j,d=r,\rho=\lambda/r$.
The retained monomials are precisely
$$
Y^bZ^c,\qquad X^aY^bZ^c\quad(a>0,\ 0\le b\le r). \tag{M6}
$$
Within each residue block take the cheapest number equal to its
dimension.
Infinite streams may be used for the birth calculation: the reference
columns $a=0,1$ already fill the block below normalized price $L$,
while every new level $\ell\ge L$ costs at least $L$.

Take irrational $0<\mu=\lambda/r<1$ for the birth comparison.
For an item $X^aY^bZ^c$, $a>0$, put $D=a+b+c$ and write
$a+b=j+r\ell$. Its positive Laurent endpoint is $a+\ell$.
Counting spine/right endpoints and left endpoints before it gives
position
$$
a+2\ell+\lfloor\mu a\rfloor+1-\mathbf1_{j=0}.
$$
It first enters at
$$
s=D+r\psi(a),\qquad
\psi(a)=\left\lfloor\frac{a+\lfloor\mu a\rfloor}{2}\right\rfloor .
\tag{M7}
$$
At fixed $D,a$ its multiplicity is $\min(r+1,D-a+1)$.

The remaining items $Y^bZ^c$ have integer price $D=b+c$.
For surplus at least $rh$, $h\ge1$, put
$$
B_h=\left\lceil\frac{2\mu h}{1-\mu}\right\rceil.
$$
Their first-cutoff surplus is at least $rh$ exactly when
$$
b\ge2h+(r+1)B_h. \tag{M8}
$$
Here is the position calculation. Write $b=j+r\ell$.
The required right endpoint has index $\ell+2h$.
Its level is
$\ell'=\max(0,\lceil(\ell+2h-j)/(r+1)\rceil)$.
It is cheaper than the left item precisely when
$\ell'\le\ell-B_h$. This is equivalent to (M8); (M8) also
forces $\ell\ge B_h$.
Consequently the surplus superlevel has multiplicity
$$
\max(0,D-2h-(r+1)B_h+1), \tag{M9}
$$
and total multiplicity at surplus zero is $D+1$.

These are exactly the birth data of (M4). To recall that comparison,
its nonconstant-price items have (M7), with the same multiplicity.
Its integer-price items have labels
$D=\delta+2u+(r+1)v$ and surplus $r\min(u,\lfloor(1-\mu)v/(2\mu)\rfloor)$.
The channel polynomial $(1+t)(1+\cdots+t^r)$ gives the superlevel series
$$
\frac{(1+t)(1+\cdots+t^r)t^{2h+(r+1)B_h}}
     {(1-t^2)(1-t^{r+1})}
=\frac{t^{2h+(r+1)B_h}}{(1-t)^2},
$$
which is (M9). Thus all price and first-cutoff multisets agree.
Continuity of rank-birth prices gives the complete tiers for every
$0\le\lambda\le r$. The outer monic case handles $\lambda>r$.
The upper and lower bounds coincide, proving Theorem 2. $\square$

The sign restriction matters. With $r=3,s=10$, $Q=Z^2,H=Z^4$,
prices $(2/3,1,1)$ and cutoff $4$, the actual rank is $41$ but
(M4) gives $42$. Dropping coprimality at negative $\lambda$ would be false.

## Uniform transfer when the interaction variable has the middle price

Return to
$$
I_{s,r}=(x^2+y^2+z^2+xyz^r)+\mathfrak m^{s+1},\qquad r\ge2.
$$

**Theorem 3.** Every generic positive weight with
$w_x>w_z>w_y$ selects the same standard set as
$$
J_{s,r}=(x^2+z^2+xyz^r)+\mathfrak m^{s+1}.
$$
These standard sets and all original prefix ranks in this open region
are independent of characteristic. The number of full original cones
meeting this region is at most $12n+7$, with $n=(s+1)^2$.

**Proof.** Write $w=(X,Y,Z)$, set $v=(1,1-r,1)$ and
$$
t=(Z-Y)/r,\qquad q=((r-1)Z+Y)/r,\qquad
\lambda=(X-Z)/q.
$$
Then
$$
w=q(1+\lambda,1,1)+tv. \tag{M10}
$$
If $X>Y+rZ$, the monic $x<2$ case is elementary.
Otherwise $0<\lambda<r$.

Let $B$ be the minimum basis of $J$ at $w$. Since $J$ is $v$-graded,
it is also minimum at $(1+\lambda,1,1)$.
The exact Laurent model has $m+d j$ parent bounds with $d=r-1$
and price ratio
$$
\rho=\frac{X-Z}{Y+(r-1)Z}=\lambda/r\in(0,1).
$$
Lemma 1 proves that its chosen reference determinant is a unit.
As $\operatorname{in}_v I=J$, lifting this graded basis gives an
integral basis of the original quotient. Every nonzero original
border coefficient has difference $\eta$ with $\eta\cdot v\ge0$.

At the seed $(1+\lambda,1,1)$, Theorem 2 applies both to the original
$Q=y^2+z^2,H=yz^r$ and to the face $Q=z^2,H=yz^r$.
The second pair is not coprime, so the earlier theorem alone was
insufficient here. Their entire price-rank inventories now agree.
The lifted $B$ is therefore a minimum-price original basis at the seed.
Any nonzero original border coefficient permits a one-column exchange,
so
$$
\eta\cdot(1+\lambda,1,1)\ge0.
$$
Equation (M10) makes every original border inequality nonnegative at
$w$. It is strict if $\eta\cdot v>0$.
If $\eta\cdot v=0$, that coefficient belongs to the graded face
normal form; genericity of the face basis makes its inequality strict.
Thus $B$ is the original standard set throughout the stated open region.

Finally the previous signed-grading theorem gives
$g(J)\le12n+7$. Transfer injects the original standard sets in this
region into that list. $\square$

The exchanged statement handles $w_y>w_z>w_x$.
Together with the proved cheapest-$z$ region, only $O(n)$ original
cones meet the open region $w_z<\max(w_x,w_y)$. At generic
weights, characteristic-dependent changes are confined to the
$z$-largest region.

### An original rank query in the proved region

The Laurent blocks in (S3)--(S5) now compute original ranks here.
For each block, count the cheap items from (M1) and cap at
$2L-\mathbf1_{m=0}$. For each right-hand offset $b=0,\ldots,r-1$,
its prices form an arithmetic progression in $j=1,\ldots,L-1$:
$$
kZ+b_0\beta+\alpha(m-(r-1)+b)
       +(\beta+\alpha(r-1))j .
$$
Remove the duplicate $a=0,j=1,b=0$ when $m=0$.
Add the spine and left-stream counts. This gives
$O(r^2(s+1))$ exact arithmetic/floor operations with constant
auxiliary storage. The statement is about arithmetic operations,
not a bit-complexity or core-CLI speed claim. Full tiers on the
region's boundary follow by continuity; arbitrary tie refinements
may leave the region.

## What remains for the extremal goal

After exchanging $x,y$ and normalizing $w_z=1$, the unclassified domain
is now the triangle
$$
0<u<1,\qquad 0<\lambda<1-u,\qquad
w=(u+\lambda,u,1). \tag{M11}
$$
The earlier all-field bound already excludes a power amplification
when $r\ge\sqrt{s}$, so only $2\le r<\sqrt{s}$ remains for this route.
The original signed transfer criterion is proved in the middle-price
region, not in (M11).

The unresolved seed is $w_x=w_y\le w_z$. A degeneration to
$z^2+Y^2+Y^2z^r$ gives a lower comparison, but no uniform equality
with the original ranks has been proved. It cannot be imported from
Theorem 2, whose high term is linear in its distinguished variable.
The prepared $z^2$ cone and the earlier bad-prime examples remain
valid inside this triangle. They do not supply a growing fan.

## Prior-art boundary and checks

The generic coprime formula is the earlier repository theorem, with its
proof in the higher-interaction note. The new argument removes that
hypothesis only on the stated nonnegative price-difference side.
Specialization of finite matrices, field extension, diagonal scaling,
Laurent division, cone inequalities and successive initials are classical.
The exact scope of this rank extension and its regional application
still require current-priority comparison; no novelty is inferred
from computation or an unsuccessful search.

A recent superficially related primary source was checked:
[Karadeniz Şen, arXiv:2503.08118v1](https://arxiv.org/pdf/2503.08118).
The introduction's jet-scheme definition, Theorem 1.1 statement,
Section 3 definitions and the $A_n$ setup in Section 4.1 were inspected.
Its fan uses **minimum** weights of a defining hypersurface for embedded
toric resolutions; its jet schemes parameterize truncated arcs.
These are different objects from the global maximum-weight initials
of $(F)+\mathfrak m^{s+1}$. No finite-quotient count or transfer of its
resolution theorem is asserted. The full paper and later version
comparison have not been completed.

~~~bash
python3 -B scripts/run_middle_interaction_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker compares the Laurent flags with independent finite-field
Gaussian elimination, checks the birth multisets, tests the noncoprime
extension against original polynomial columns and a negative-side
counterexample, and certifies entire original middle-price patches by
integer/rational border relations. The uniform statements rest on the
proofs above. The unrestricted $G_3$ objective remains active.
