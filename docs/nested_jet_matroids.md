# Explicit polynomial quotients with a nested rank rule

The earlier controls separated arbitrary representable matroids from residue
matroids of finite polynomial quotients. This note supplies an infinite,
explicit, non-binomial family for which that gap can be closed.
All monomial ranks and the quotient length are known without coefficient
search. The family includes embedding dimension three, as well as a
curvilinear subfamily. It works in characteristic zero and uses classical
interpolation. It does not yet give a new asymptotic bound for $G_3(n)$
or establish novelty.

## A fixed quotient realizes the whole rank rule

Fix an ambient truncation $L\ge1$ and any positive integers
$m=(m_1,m_2,m_3)$. Put

$$
B=\max(L,2),\qquad C=L(1+B+B^2),\qquad
\lambda_i=Cm_i+B^{i-1}.
$$

In $\mathbb Q[t]/(t^L)$ define, using truncated formal exponentials,

$$
X_i=t^{m_i}e^{\lambda_it}.                                 \tag{1}
$$

Let $I_{L,m}$ be the kernel of the map
$\mathbb Q[x,y,z]\longrightarrow\mathbb Q[t]/(t^L)$ sending the variables
to $(X_1,X_2,X_3)$, and let $A$ be its image subalgebra. The polynomial
quotient is exactly $A$, which need not be the whole ambient truncation.
Its three multiplication operators commute on one common cyclic algebra.
It is supported at the origin. No parameter in (1) depends on an objective
weight or a candidate basis. Theorem 1 computes its actual length $n$;
confusing $L$ with that length would invalidate a counting comparison.

For an exponent $u\in\mathbb N_0^3$, write

$$
v(u)=m\cdot u,\qquad \xi(u)=\lambda\cdot u.
$$

Its residue is exactly $t^{v(u)}e^{\xi(u)t}\bmod t^L$. It is zero if
$v(u)\ge L$. On the finite set $T=\{u:v(u)<L\}$, all nodes $\xi(u)$
are distinct and their order refines the order of $v(u)$. Indeed

$$
\xi(u)=Cv(u)+u_1+Bu_2+B^2u_3,
\qquad 0\le u_1+Bu_2+B^2u_3<C.
$$

Every $u_i<L\le B$, so the last expression is an injective base-$B$
encoding within each valuation. This is the reason for the explicit
parameters, rather than a generic choice asserted separately for each basis.

**Theorem 1.** A set of $r$ distinct nonzero monomials, with sorted
valuations $v_1\le\cdots\le v_r$, is independent in $A$ if and only if

$$
v_i\le L-r+i-1\quad(1\le i\le r).                         \tag{2}
$$

For any finite monomial set $S$, including zero residues, its rank is

$$
\operatorname{rk}(S)=
\min\left(|S|,\min_{0\le k\le L}
       \left(L-k+\#\{u\in S:v(u)<k\}\right)\right).        \tag{3}
$$

In particular, writing $\nu_m(k)=\#\{u:m\cdot u<k\}$, the polynomial
quotient's length is exactly

$$
n=L-\max_{0\le k\le L}\bigl(k-\nu_m(k)\bigr).             \tag{3a}
$$

This is a nested, or Schubert, matroid rank rule on the monomial ground.
Those matroid classes and their ordinary linear representations are
established. The point here is the explicit multiplicative realization
(1), not the invention of a new matroid class or a new rank formula.

### The interpolation input, with an elementary proof

We use the following classical ordered Hermite--Birkhoff special case.
Suppose $\xi_1<\cdots<\xi_N$ and
$0\le d_1\le\cdots\le d_N$, with $d_i\le i-1$. The functionals

$$
p\longmapsto p^{(d_i)}(\xi_i),\qquad \deg p<N,
$$

are independent. Here is a short Rolle proof to make the exact special
case explicit. Suppose all functionals vanish on $p$, and let $q\ge1$
be the number of zero derivative orders. The first $q$ nodes are roots
of $p$. Between them Rolle gives $q-1$ distinct roots of $p'$, all
before the remaining original nodes. At those remaining nodes, the
orders for $p'$ are $d_i-1$. The new $N-1$ orders are still monotone
and meet the same inequalities: the order-zero part has length $q-1$,
and at new index $i-1$ a remaining order is at most $i-2$.
Induction implies $p'=0$; its original zero then gives $p=0$.
The case $N=1$ starts the induction.

The broader poisedness theorem is
[Atkinson--Sharma (1969), Theorem 2](https://homepage.math.uiowa.edu/~atkinson/papers/HB_interp.pdf).
Their Pólya condition and conservative incidence matrices include this
monotone special case. That theorem's definitions, statement and proof
were inspected; the elementary interpolation assertion above is classical.

### Proof of the rank theorem

The column of $t^v e^{\xi t}$ in the basis $1,t,\ldots,t^{L-1}$ has
entries $\xi^{j-v}/(j-v)!$ for $0\le j<L$ and $j\ge v$, and zero
otherwise. Multiplying row $j$ by $j!$ turns the transposed column into the functional
$p\mapsto p^{(v)}(\xi)$ on polynomials of degree less than $L$.
For $r=L$, (2) is precisely the interpolation special case just proved.
The determinant is a nonzero rational number, so the same independence
holds over any characteristic-zero field.

For $r<L$ satisfying (2), add $L-r$ auxiliary conditions of derivative
order zero at fresh rational nodes preceding every original node. The
combined $L$ conditions satisfy the monotone inequalities, so their full
matrix is invertible. The original $r$ columns are consequently independent.
These auxiliary columns are a proof device; they need not be monomials of
the quotient.

Conversely, if (2) fails at $i$, the last $r-i+1$ columns all belong to
$t^{v_i}\mathbb Q[t]/(t^L)$, of dimension $L-v_i<r-i+1$.
They are dependent. The case $r>L$ is already dependent by dimension.
This proves (2). For (3), each threshold $k$ gives the displayed rank
upper bound using the ambient subspace $t^k\mathbb Q[t]/(t^L)$.
Let $R$ be the minimum in (3).
The $R$ smallest valuations satisfy (2): otherwise, at a violating $i$,
using $k=v_i$ (or $k=L$ for a zero residue) would make that minimum
strictly less than $R$. Thus $R$ independent columns exist, proving (3).
Applying it to the spanning set $T$ gives (3a); its $k=L$ term already
accounts for $|T|$.

The length defect is bounded for fixed $m$. If $M=\max m_i$, the cube
$0\le u_i\le\lfloor(k-1)/(3M)\rfloor$ supplies at least $(k/(3M))^3$
monomials of valuation less than $k$. Thus $\nu_m(k)\ge k$ for
$k\ge(3M)^{3/2}$. The maximum in (3a) stabilizes once $L$ exceeds
that bound. This is a direct lattice count, not a new asymptotic fan bound.

These quotients need not have embedding dimension one. For $m=(3,4,5)$
and every $L\ge10$, (3a) gives $n=L-2$. The numerical semigroup
generated by $3,4,5$ contains every integer at least three. One monomial
of each such valuation, together with one, gives $\nu_m(k)\ge k-2$;
equality is attained at $k=3$. The generator valuations $3,4,5$ precede
every product of two nonconstant generators. Their images in the maximal
ideal modulo its square are independent and generate it, so the embedding
dimension is exactly three.

This entire family is non-binomial for $L\ge10$. The monomials just
chosen give $n$ different nonzero residue lines. At valuation eight,
$xz$ and $y^2$ have distinct nodes and remain nonproportional modulo
$t^L$, since their next coefficient is visible. Choosing both gives
at least $n+1$ lines of rank $n$. As explained below, this is impossible
for a binomial quotient. Thus the constructive route reaches embedding
dimension three at arbitrarily large lengths.

## Ordinary polynomial generators

Now specialize to $L=n$ and $m=(1,a,b)$, and abbreviate the kernel by
$I_{n,a,b}$. The powers $1,X_1,\ldots,X_1^{n-1}$ have a triangular
coefficient matrix with diagonal one, so $A=\mathbb Q[t]/(t^n)$.
This subfamily is curvilinear and has length exactly $n$. In this
subsection write $(X,Y,Z)=(X_1,X_2,X_3)$.

The construction is a finite polynomial quotient, not a transcendental
coefficient model. It can also be written explicitly as

$$
I_{n,a,b}=(x^n,\ y-P_{a,\lambda_2}(x),\ z-P_{b,\lambda_3}(x)), \tag{4}
$$

where $P_{c,\mu}=0$ if $c\ge n$, and otherwise

$$
P_{c,\mu}(x)=x^c+
 \sum_{j=c+1}^{n-1}
 \frac{(\mu-c\lambda_1)(\mu-j\lambda_1)^{j-c-1}}{(j-c)!}\,x^j.
                                                               \tag{5}
$$

To derive (5), let $t=\tau(x)$ invert $x=t e^{\lambda_1t}$.
Classical Lagrange inversion gives

$$
[x^j]\,\tau(x)^c e^{\mu\tau(x)}
 =\frac1j[t^{j-1}](c t^{c-1}+\mu t^c)e^{(\mu-j\lambda_1)t}.
$$

At $j=c$ this is one; at $j>c$ it simplifies to (5).
See [Gessel, Theorem 2.1.1, equation (2.1.1)](https://arxiv.org/html/1609.05988v1#S2.SS1).
The theorem and its formal-series assumptions were checked. Lagrange
inversion and this coefficient calculation are not new results.
Since (4) has length $n$ and maps to (1), it is exactly the kernel.

For $a\le n-2$, this ideal is not binomial. The $n$ independent powers
of $X$ lie on $n$ different residue lines. $Y$ lies on none of those
lines: different valuations exclude all powers except $X^a$, and the
next coefficient differs by $\lambda_2-a\lambda_1=B-a\ne0$.
Thus there are at least $n+1$ distinct nonzero residue lines of rank $n$.
A binomial quotient's nonzero monomial residues form independent
one-dimensional classes: its generating relations only identify scalar
multiples of monomials, or kill a class. Therefore it cannot have this
rank pattern. This distinguishes the family from the already bounded
binomial subclass.

## Two uniform fan consequences

If $m=(1,1,1)$, every $n$-element monomial downset is a basis. In a
downset, for each $1\le k\le n$, either a monomial of degree at least
$k-1$ supplies a divisor chain of $k$ distinct monomials of degree less
than $k$, or every monomial already has smaller degree. Thus its sorted
degrees meet $v_i\le i-1$, and (2) applies.

At every strict positive objective, the cheapest $n$ monomials form such
a downset. They already constitute a basis, so greedy selection never
leaves that ordinary prefix. Consequently

$$
g(I_{n,1,1})=C_3(n),                                      \tag{6}
$$

with exactly the generic corner-cut positive fan. Large additional linear
dependencies in a quotient need not change that fan at all. The argument
extends to every fixed number of variables using the same base encoding.
This is a comparison consequence, not a spatial power separation.

There are also uniform non-corner-cut initial sets. For every $q\ge2$,
use $n=2q$ and $m=(1,1,q)$. If $E$ is any $q$-element planar corner cut,
then

$$
D=E\times\{0,1\}                                         \tag{7}
$$

is an initial standard set for the same fixed quotient $I_{2q,1,q}$.
The first $q$ valuations are the degrees of $E$, all at most $q-1$;
the other $q$ are those degrees plus $q$. Both blocks satisfy (2).
Every monomial divisible by $z^2$ is zero. Choose positive planar weights
strictly separating $E$, and choose a positive $z$ weight smaller than
their cutoff gap. Then the cheapest $2q$ nonzero monomials are exactly
(7), already an independent set. Hence they are the greedy standard set.

Set (7) is never a corner cut. Choose a largest pure-axis power $u\in E$
of positive degree. Then $uz\in D$, while $u^2,z^2\notin D$, and
$2\exp(uz)=\exp(u^2)+\exp(z^2)$. A strict linear cutoff cannot include
this midpoint and exclude both endpoints. This proves an infinite
realizable departure from smooth-number prefixes, rather than inferring
one from an isolated finite example. It gives $g(I_{2q,1,q})\ge C_2(q)$,
which is not a spatial extremal lower bound of the required strength.

## Balanced valuations give exactly bounded-height lattice prefixes

Let

$$
\mathcal S_h=\mathbb N_0^2\times\{0,\ldots,h-1\},
$$

and let $S_h(n)$ count its distinct $n$-element prefixes under strict
positive linear costs. These are cuts relative to the slab; they need
not be corner cuts of the whole orthant.

**Theorem 2.** For every $h,q\ge1$, the fixed quotient $I_{hq,1,q}$
has exactly the positive prefix fan of $\mathcal S_h$ at cardinality
$n=hq$. Thus $g(I_{hq,1,q})=S_h(hq)$.

To prove it, consider any $n$-element downset $D\subset\mathcal S_h$.
Write its planar layers as $D_c$, with sizes
$a_0\ge\cdots\ge a_{h-1}\ge0$ summing to $hq$. The first $c$
layers have at least $qc$ elements, because their average is at least
the average of all $h$ layers. In each $D_c$, its $j$th monomial in
increasing total degree has degree at most $j-1$, by the divisor-chain
argument above. Assign it the distinct index
$a_0+\cdots+a_{c-1}+j-1$. This index is at least its valuation
$qc+\deg u$. The indices run from zero to $n-1$, so the sorted
valuations satisfy (2). Hence **every** such downset is a basis.
Moreover $z^h=0$ in this quotient. The cheapest $n$ slab monomials
are a downset, already a basis, so greedy selection gives exactly that
prefix at every strict objective. No coefficients are chosen after the
objective is known.

### The same fan has radical point realizations

The slab description is not confined to local jets. For any $n,h$, choose
$h$ distinct elements $\zeta_i$ of an infinite field, and place
$q_i=\lfloor n/h\rfloor+1$ points on each of the first $n\bmod h$
planes $z=\zeta_i$, and $q_i=\lfloor n/h\rfloor$ on the others.
There is a nonempty Zariski-open choice of the planar point coordinates
for which every $n$-element downset in $\mathcal S_h$ is an evaluation
basis. In particular rational, distinct point realizations exist over
$\mathbb Q$.

Here is the determinant proof of that assertion. For one downset, group
its monomials by their $(x,y)$ exponent. Each group has $k\le h$
columns with $z$ powers $0,\ldots,k-1$. Assign the group's $k$ columns
to distinct planes by taking the next $k$ slots in a cyclic list of the
$h$ planes. After all $n$ columns, plane loads are exactly the prescribed
$q_i$. Assign the columns to different point rows in those planes.
In the determinant polynomial in the independent planar row coordinates,
the coefficient of the resulting row-by-row $(x,y)$ monomial is, up to
sign, a product of ordinary Vandermonde determinants in the assigned
distinct $\zeta_i$. It is nonzero. Other assignments of $(x,y)$ groups
give different monomials in the separate row variables, so cannot cancel
this coefficient. Thus each evaluation determinant is a nonzero polynomial.
There are finitely many downsets of size $n$, so the product of these
determinants is nonzero. Over an infinite field its nonvanishing locus,
also requiring distinct points, is nonempty.

The point ideal contains $\prod_i(z-\zeta_i)$, whose head is $z^h$
for every global order. Monomials beyond the slab reduce to cheaper
ones within it. Consequently its positive fan is also exactly the slab
prefix fan. These outside monomials are generally nonzero in the point
quotient; the jet's nilpotence is not asserted for the point ideal.

## Uniform slab bounds and a sharp fixed-height count

For all $n\ge2$ and $h\ge2$ there is an absolute upper bound

$$
S_h(n)=O\!\left(n^2\log(n+1)
             (1+\log(\min(h,n)))\right).                  \tag{8}
$$

For each fixed $h\ge2$, the sharp order is

$$
S_h(n)=\Theta_h(n^2\log(n+1)),                            \tag{9}
$$

whereas $S_1(n)=C_2(n)=\Theta(n\log(n+1))$. In particular no
choice of height in the balanced family gives a polynomial separation
above the generic three-variable scale. This conclusion is uniform in
height; it is not inferred from the finite fan counts.

### The upper bound is a classical geometric corollary

More generally consider a positive linear cut of a box
$\prod_i\{0,\ldots,H_i-1\}$, allowing infinite side bounds. Let $\ell_i$ be the
number of included points on the $i$th coordinate axis, including zero.
Then $D$ lies in their product box, and contains the smaller box with
side lengths $\lceil\ell_i/d\rceil$. Indeed, if the cutting level is
$T$, every point of the smaller box has cost at most
$d^{-1}\sum_i w_i(\ell_i-1)<T$. The ambient box contains these points.
Therefore

$$
n\le\prod_i\ell_i\le d^d n.                              \tag{10}
$$

Put $j_i=\lceil\log_2\ell_i\rceil$. The dyadic product
$2^{\sum j_i}$ is between $n$ and $2^d d^d n$, so its exponent lies
in a fixed-width interval about $\log_2n$. For a slab in dimension
three, $j_3\le\lceil\log_2h\rceil$. There are at most

$$
O\!\left(\log(n+1)(1+\log(\min(h,n)))\right)
$$

possible triples. In one triple all relevant prefix monomials lie in
the fixed box with side counts $2^{j_i}$, intersected with the original
slab. Its volume in lattice-point counts is $O(n)$. Differences of its
points have at most $\prod_i(2^{j_i+1}-1)=O(n)$ possible integer
values. Their central comparison planes have $O(n^2)$ chambers in
three-dimensional weight space. On a chamber the entire box ordering,
and hence the prefix, is fixed. Every strict prefix persists under a
small perturbation avoiding all these comparison planes. Multiplying
the two counts proves (8). The same argument gives the familiar
$O_d(n^{d-1}\log^{d-1}(n+1))$ bound for relative cuts of arbitrary
axis boxes. Intersecting the dyadic box with the original caps is
essential: adding cheaper points beyond a cap would change the prefix.

This uses the established product/dyadic/arrangement mechanism. It is
**not** a new lattice-polytope theorem or an unrestricted ideal bound.
Compare [Wagner's thesis, Section 3.4, especially Observations 3.13--3.14
and its subsequent product encoding](https://ti.inf.ethz.ch/ew/research/phdtheses/wagner03.pdf),
and the repository's existing ordered-prefix proof. The specified
primary argument was inspected. Its unbounded-orthant statement and
the relative-box version are distinguished; the short proof above
supplies the latter as a direct geometric corollary.

### Two layers give an exact prefix-pair formula

At fixed generic positive planar weights let
$\beta_1=0<\beta_2<\cdots$ be the ordered planar monomial costs.
Every two-layer prefix has $a$ lower-layer terms and $b$ upper-layer
terms, with $a+b=n$ and $a\ge b$. Conversely every such pair occurs
for these planar weights: when $b\ge1$, choose the $z$ cost $u$ in

$$
\max(0,\beta_a-\beta_{b+1})
       <u<\beta_{a+1}-\beta_b.                            \tag{11}
$$

This interval is nonempty and makes the first $a$ planar monomials at
height zero and the first $b$ at height one strictly precede every
excluded slab monomial. For $b=0$ take $u>\beta_n$.

Let $\mathcal B_r$ be the finite set of positive planar weight ratios
at which the $r$-element corner cut changes, and put
$\mathcal B_0=\varnothing$. Then

$$
S_2(n)=\sum_{b=0}^{\lfloor n/2\rfloor}
       \left(1+|\mathcal B_{n-b}\cup\mathcal B_b|\right). \tag{12}
$$

For fixed $b$, the pair of planar prefixes changes exactly at that union.
Each pair has a connected parameter interval, so cannot recur in a
disjoint interval. Since $1+|\mathcal B_r|=C_2(r)$ for $r\ge1$,
(12) and the known planar estimate already prove (9) for $h=2$.
The published planar count is due to
[Corteel--Rémond--Schaeffer--Thomas (1999)](https://doi.org/10.1006/aama.1999.0646).
The interleaving and this application of that count do not settle the
priority of (12).

For completeness the lower bound extends to every $h\ge2$, uniformly.
Take $1\le b\le\lfloor(n-2)/5\rfloor$ and $a=n-b$. The planar counting
function $N(T)=\#\{v:w\cdot v<T\}$ obeys $N(2T)\le4N(T)$,
by mapping $v$ to $\lfloor v/2\rfloor$. Hence
$\beta_a\ge2\beta_{b+1}>2\beta_b$. There is a choice in (11)
also satisfying $u>\beta_a/2$ and $u>\beta_b$, because its upper
endpoint exceeds each of these lower requirements. Then $2u$ exceeds
the cost of every selected monomial, so all layers at height at least
two are excluded. For every planar cut of size $a$, choose generic
weights in its parameter interval and its corresponding prefix of size
$b$. This gives at least $C_2(a)$ slab cuts with that upper-layer count.
Different $b$ give different sets. Summing over a linear number of $b$
proves $S_h(n)=\Omega(n^2\log(n+1))$ and completes (9).

This resembles Wagner's established dimension-lifting lower argument
(thesis Lemma 3.5 and Theorem 3.6); the doubling estimate keeps the lift
within two layers. Its restricted sharpness is not the sharp $G_3(n)$
theorem required by the active goal.

### A further algebraic output bound

A slab downset with planar layer sizes $a_c$ has at most
$O(h+\sum_c\sqrt{a_c})=O(h+\sqrt{hn})$ minimal outside monomials.
At height $c<h$ every such monomial projects to a minimal outside
monomial of the planar layer; the latter number is $O(1+\sqrt{a_c})$
by the usual antichain-chain packing bound. At height $h$ only the
possible generator $z^h$ needs to be added. Each reduced basis has
one polynomial per minimal initial generator. Thus the slab point and
balanced-jet ideals satisfy, for fixed $h\ge2$,

$$
|\mathcal U(I)|=O_h(n^{5/2}\log(n+1)).                   \tag{13}
$$

This is an output-cardinality consequence of (9) and a classical union
bound. It is not a new Gröbner-basis algorithm or a field-operation bound.

## General valuations are not covered by the slab argument

The reduction just proved uses $L=hq$ and $m=(1,1,q)$. A general
nested quotient need not have convex initial standard sets. At $L=8$,
$m=(1,3,4)$ and strict weight $(3,2,1)$, its standard set is

$$
\{1,x,x^2,y,y^2,z,xz,yz\}.                               \tag{14}
$$

It includes $x^2,y^2$ and excludes their lattice midpoint $xy$.
The checker independently verifies (14) by rational quotient columns,
border polynomials and all S-pairs; every border inequality has positive
margin. Thus lattice convexity, which would force a bounding-box product
of $O(n)$, cannot be assumed for the whole nested family. This is a
diagnostic counterexample to that proof shortcut, not an asymptotic
lower construction.

### An infinite family also rules out a uniform small bounding box

The midpoint example does not by itself disprove a bounding-box product
bound established by some other method. The following family does.
For $M\ge2$ set $L=n=M^2$, $m=(1,M,M)$ and take the strict weight

$$
w=((M+1)^3,M+1,M+2).                                    \tag{15}
$$

Its initial standard set has exactly $M$ axis points in each coordinate
direction, counting zero. Its bounding box therefore contains
$M^3=n^{3/2}$ lattice points. A uniform $O(n)$ bounding-box product
for all nested initial sets is false.

Here is the proof. All nonzero monomials in $y,z$ have total degree
at most $M-1$. The whole triangle of $\binom{M+1}{2}$ such monomials
is independent: its degree-$d$ columns have valuation $Md$ and meet
(2), since

$$
Md-\binom{d+1}{2}\le M(M-1)/2\quad(0\le d<M).
$$

They all precede $x$ at (15), so all are selected. Every basis also
contains $1,x,\ldots,x^{M-1}$: these are the only monomial columns
with nonzero images needed to span the first $M$ coefficient positions
modulo $t^M$. The $M^2$ monomials
$x^a y^d$, $0\le a,d<M$, have all valuations $0,\ldots,M^2-1$
once each and hence are independent. They all precede $x^M$ at (15).
Thus the entire quotient is spanned before $x^M$ is reached, and no
monomial with $x$ exponent at least $M$ is selected. Since $y^M,z^M$
are zero, all three claimed axis lengths are exactly $M$. Different
$x$ layers are completely separated by (15); within a layer the costs
are distinct because $M+1,M+2$ are coprime and the planar degrees
are less than $M$. The weight is consequently strict on all nonzero
ground monomials, not an unresolved tie.

The layer ranks are also explicit. Let $R_r$ be the rank of all nonzero
monomials with $x$ exponent at most $r$, for $0\le r<M$. Then

$$
R_r=M^2-\delta_r,\qquad
\delta_r=\max_{0\le d\le M}
       \left(Md-(r+1)\binom{d+1}{2}\right).               \tag{16}
$$

At a threshold $k=Md+s$, $0\le s\le M$, the number of those
monomials below $k$ is

$$
(r+1)\binom{d+1}{2}+(d+1)\min(s,r+1).
$$

The deficit $k$ minus this number is nonincreasing up to $s=r+1$ and
then increases. Its maximum is therefore at a multiple of $M$;
applying (3) proves (16). In the lexicographic layer order at (15),
the number of selected monomials in layer zero is $R_0$ and in
layer $r\ge1$ is $R_r-R_{r-1}$. This is an exact formula, not a fit
to the finite data. It identifies a square-root valuation regime with
large coherent footprints and simultaneous rank constraints, without
asserting that its number of initial ideals has a larger power.

The balanced slab subfamily cannot provide the required polynomial
separation and is parked as a source of that lower bound. The family
(15)--(16) remains outside that reduction. The next task is to count
how its rank-defect layers can interact as the objective varies, or to
prove a uniform structural bound for general valuations. Neither the
large footprint nor the fixed-height count is a lower bound of the
required strength. Priority is unassessed and the unrestricted goal
remains open. No paper is proposed on these results alone.

## A combinatorial oracle, and the remaining counting problem

The original smooth-number order still scans exponents by $w\cdot u$,
with $w_i=\log p_i$ for a prime basis. Here the separate additive quantity
$m\cdot u$ imposes the nested rank capacities. The same positive exponent
scan can therefore skip monomials for a precisely known algebraic reason.
Equations (6)--(7) show both possibilities: dependencies can preserve the
generic fan or produce shapes forbidden to ordinary smooth prefixes.
No new prime-distribution consequence follows from this construction.

Independence (2) is equivalent to scheduling unit jobs with deadlines
$L-v(u)$. In minimum-cost greedy order, assign each accepted monomial
the latest available slot at or before its deadline. A failed assignment
is precisely dependence. This is the classical nested scheduling matroid
algorithm. It gives an exact combinatorial optimization oracle on the
established ground $V_n^3$, without floating rank tests or experimental
coefficient choices. It does not enumerate the whole fan or prove a
bound on its number of cells.

Schubert and lattice-path matroids, including their ordinary generic
matrix representations, are explicitly described in
[Benedetti-Velásquez--Knauer (2024), Section 2.1 and Remark 3](https://link.springer.com/article/10.1007/s00493-024-00085-4).
Those portions were inspected. A matroid quotient in that paper is a
map between matroids, not a polynomial quotient with commuting variables.
Neither that terminology nor the classical interpolation input establishes
novelty of (1), nor do they give the spatial fan count required here.

The next proof task is to bound the number of exposed profiles of (2)
with the actual monomial weights, uniformly in $L,m$, charging the true
length (3a). Either a power lower construction or a generic-scale upper
theorem would decide whether
this family can advance the extremal target. Larger finite fan counts or
arbitrary claims about all Schubert matroids would not decide it. The
surviving-multiples structure is present: relations affect whole tails of
the same cyclic algebra, rather than independent socle pairs. The fixed
carrier equations in (4) can have degree $n-1$, so the earlier
fixed-degree curve restriction does not automatically settle the family.

## Exact supporting checks

Run from the repository root:

```bash
python3 -B scripts/run_nested_jet_checks.py
python3 -B scripts/run_nested_slab_checks.py
./build.sh
python3 tests/test_smoke.py
```

The dependency-free checker compares (3) against rational Gaussian
elimination on every subset of several small full valuation grounds,
and separately evaluates permutation determinants. It compares direct
truncated multiplication with exponential columns and the explicit
generators (5) in the curvilinear subfamily, and checks all multiplication
commutators. It independently checks the maximal-ideal-square rank and
non-binomial residue lines in embedding-dimension-three examples. Complete
positive fans are computed by rational polygon clipping. Their cones
are checked both by matrix-derived border relations and by the nested
basis-exchange rule, with ordinary S-pair checks on the smaller fans.
Two full comparison arrangements independently check the fan traversal.
The exact report includes all standard sets, strict objectives and cells.
These are supporting controls, not an asymptotic fit or a priority proof.

The [slab report](../results/research/nested_slab_checks.json) adds exact
rational rank checks on every slab downset in five balanced jet fixtures,
and modular evaluation determinants for every slab downset of eight
radical point fixtures. Nonzero modular determinants also certify the
integer point lifts over $\mathbb Q$. Scalar prefix fans and matrix-derived
cones agree, and four two-layer fans are independently reconstructed from
the planar prefix-pair formula. It checks the small midpoint counterexample
with all S-pairs, the large-footprint and layer-defect formulas, and the
matrix layer ranks in the smallest members. The earlier checker is
replayed with its previous mathematical results unchanged.

These checks accompany the uniform proofs. The upper bound (8) is a
classical geometric corollary; exact priority of the slab realizations,
prefix-pair formula and fixed-height consequences is unassessed. A failed
search is not evidence of novelty. No additional manuscript is created.

All previous manuscripts, smooth-number kernels and the paused ML work
are preserved. The sharp $G_3(n)$ question and the active goal remain open.
