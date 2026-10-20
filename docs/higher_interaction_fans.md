# Higher interaction order: arithmetic changes and a uniform size bound

For any field, consider
$$
F_r=x^2+y^2+z^2+xyz^r,\qquad
I_{s,r}=(F_r)+\mathfrak m^{s+1},\qquad r\ge1. \tag{H1}
$$
The cubic case has the [complete linear classification](cubic_fan_classification.md).
Growing interaction order introduces coherent monomial prefixes with arbitrary
cyclic torsion, hence exceptional characteristics of unbounded size. It does
not by itself give more initial ideals. A uniform bound below excludes a
polynomial gain over the quadratic scale for $r\ge\sqrt{s}$. The unrestricted
extremal goal, novelty and publication significance remain open.

## Integral model and exact grading lattice

Monic division $x^2\mapsto-y^2-z^2-xyz^r$ decreases the $x$ exponent
and never decreases total degree. Consequently
$$
\mathcal B_s=\{x^ay^bz^c:a<2,\ a+b+c\le s\},\qquad n=(s+1)^2 \tag{H2}
$$
is a basis over $\mathbb Z$ and after every field reduction. Each multiplication
column has at most three entries. The canonical monomial part is exactly
$\mathfrak m^{s+1}$: a smaller monomial relation would have lowest homogeneous
part divisible by the nonmonomial $g=x^2+y^2+z^2$. This is impossible even in
characteristic two.

For $s\ge r+2$, the intrinsic exponent-difference lattice is
$$
L_r=\mathbb Z(2,-2,0)+\mathbb Z(2,0,-2)+\mathbb Z(1,-1,-r),
\qquad [\mathbb Z^3:L_r]=4r. \tag{H3}
$$
The equation and cutoff are homogeneous for $\mathbb Z^3/L_r$. Conversely
the lexicographic standard set is (H2), and its reduced basis contains $F_r$,
with all four terms below the cutoff. Those differences generate $L_r$.
The determinant is an integer determinant, also when the field characteristic
divides $4r$. Every vector basis takes a fixed number of columns in each
group degree, so all its exponent-sum profiles lie in one coset of $L_r$.
This uses the established [grading/profile mechanism](profile_lattice_bound.md).

The substitution $X=xz^{r-1},Y=yz^{r-1},Z=z^r$ satisfies
$X^2+Y^2+Z^2+XYZ=z^{2r-2}F_r$. For $r>1$ it collapses the whole plane
$z=0$. It is not finite at the origin, so the earlier finite-pullback
theorem cannot transfer the cubic count. Localizing this nilpotent
quotient at $z$ would give the zero ring.

## A bound uniform in both parameters and the field

**Theorem 1.** If $1\le r\le s-2$, then, with an absolute constant,
$$
g(I_{s,r})=O(n^{9/4}/r). \tag{H4}
$$
If $r>s-2$, the interaction vanishes beyond the cutoff: there are three
initial ideals for $s\ge2$, and one for $s<2$.

**Proof.** Set $e=r+2$, $K=\lceil(s+1)/e\rceil$, and
$D=\min(s,2K-1)$. Since $g^K\equiv(-xyz^r)^K\pmod {F_r}$,
we have $g^K\in I_{s,r}$. In a weight sector where variable $i$ has
strictly largest price, its leading monomial is $x_i^{2K}$, in every
characteristic. Every initial standard set there is contained in
$$
E_i=\{a\ge0:|a|\le s,\ a_i\le D\}. \tag{H5}
$$
The restricted original-column matroid on $E_i$ has rank $n$, since that
sector has an initial standard set. Let $P_i$ be its basis-profile polytope.
Every initial set in the sector is uniquely exposed in $P_i$. A generic
positive objective gives distinct monomial costs; matroid greedy selection
then gives a unique optimal basis. Any global order can be represented by
such an objective for its finite initial ideal, so the three sectors suffice.

The profile box has sides at most $nD,ns,ns$. For a full-dimensional $P_i$,
(H3) bounds its normalized lattice volume by $n^3Ds^2/(4r)$.
The classical three-dimensional lattice-vertex theorem gives
$$
f_0(P_i)=O(n^{3/2}s\sqrt{D/r}). \tag{H6}
$$
If $P_i$ is two-dimensional, its Euclidean area is $O((ns)^2)$ and its
affine integer lattice has covolume at least one. The planar lattice-vertex
bound gives $O((ns)^{2/3})=O(n)$ vertices. A point or segment has at most two.
Thus no zero ambient volume is used for a smaller-dimensional polytope.
Finally $D\le C s/r$ for $r\le s-2$. Substitution in (H6) and summation
over the three sectors prove (H4). $\square$

More generally this argument gives
$$
g(I)=O\big((ns)^{2/3}
+n^{3/2}s\sqrt{\min(s,2K-1)/h}\big) \tag{H7}
$$
for a length-$n$ ideal containing $\mathfrak m^{s+1}$ and
$(x^2+y^2+z^2)^K$ and homogeneous for a finite exponent grading of
index $h$, with $s,K\ge1$. This is a specialization of classical state
geometry and relative lattice volume, not a new lattice-polytope theorem.

In particular $r\ge\sqrt{s}$ implies $g(I_{s,r})=O(n^2)$.
A family with $g(I_{s,r})\ge c n^{2+\epsilon}$ must have
$r=O(n^{1/4-\epsilon})$. No matching lower count or sharp classification
for the remaining parameters follows.

## The boundary remains exactly solvable

The initial ideal at $(1,1,0)$ is
$$
J_{s,r}=(x^2+y^2+xyz^r)+\mathfrak m^{s+1}. \tag{H8}
$$
Inclusion and the two integral monic bases of the same length prove equality.
For binary degree $k\ge1$, put $L=s-k+1$. Its columns are $t^az^c$,
$0\le a\le k$, $0\le c<L$, in
$K[t,z]/(t^2+z^rt+1,z^L)$. Split $c$ modulo $r$: residue $j$ has
$L_j=1+\lfloor(L-1-j)/r\rfloor$ and is the cubic note's binary block
with $Z=z^r$.

For $w=(u+\lambda,u,1)$, $0<\lambda<r$, put $\mu=\lambda/r$.
In each residue select the cheapest $2L_j$ from
$$
\{x^ay^{k-a}z^j:0\le a\le k\},\quad
\{y^kz^{j+rd}:1\le d<L_j\},\quad
\{x^kz^{j+rd}:1\le d<L_j\}. \tag{H9}
$$
The relative prices are $\mu a,d,\mu k+d$. The preceding three-stream
independence proof applies over every field. Include $1,z,\ldots,z^s$
in binary degree zero. This formula describes the boundary, not all
original ranks in (H1).

For $s\ge r+2$, its complete positive fan has exactly
$$
2M_{s,r}+4,\qquad
M_{s,r}=\#\{p/q\in(0,1):\gcd(p,q)=1,\
\epsilon(p,q)((r+2)q+rp)\le2(s+r)\}, \tag{H10}
$$
where $\epsilon=1$ for two odd integers and $2$ otherwise.
In particular the boundary count is $O(1+(s/r)^2)$, uniformly.

A spine/left-arm tie $\mu a=c$ changes the selected set exactly when
$a+c=2L_j$. Write $a=mq,c=mp$. Then $m(p+q)$ is even, and its smallest
possible cutoff is $s=m((r+2)q+rp)/2-r$, attained at $k=a,j=0$.
Increasing $k$ realizes the same tie at every larger cutoff. Minimal $m$
is $\epsilon$. A left/right-arm tie has $k+c+d=2L_j$, $c-d=mp,k=mq$,
with $d\ge1$; its minimal cutoff is at least $r$ larger and introduces
no additional fractions. There is no spine/right-arm tie in $(0,1)$.
Every wall changes one independent block. Different intervals cannot
give the same standard set, by convexity of its initial cone along
this line of weights. The two signs of $w_x-w_y$ already give
distinct sets in degree $k=s$. The two outer cones have leading terms
$x^2,y^2$ and contain $xyz^r$, unlike central sets. This proves the count.

The simplest extension of the cubic boundary cover fails. At $s=10,r=3$,
an entire all-field initial cone contains $(419,68,221)$. Two border
comparisons are
$$
3w_x-3w_y-3w_z\ge0,\qquad -2w_y+2w_z\ge0. \tag{H11}
$$
On $w_x=w_y\ge0$ these force $w_z=w_y=0$. Thus this cone misses
every nonzero nonnegative ray of that plane, including $(1,1,0)$.
The first comparison comes from the nonzero $y^3z^5$ coefficient in
$$
x^3z^2=yz^7+y^3z^5-xz^4-xy^2z^2.
$$
The checker saves the full standard set and integral border relations,
with strict slacks and monic membership certificates. Those certificates
and the length verify the cone over every field. A different cover
would require a new proof.

## Arbitrary cyclic torsion in one coherent price prefix

**Theorem 2.** For any integer $m\ge2$, set
$$
r=2m-1,\ s=m(2m+1),\ D=2m+2,\quad
w=(1,1+1/(8m),1+1/(4m)),\quad
\beta=xy^3z^{2m-2},\ T=w\cdot\beta. \tag{H12}
$$
In $A_{\mathbb Z}=\mathbb Z[x,y,z]/I_{s,r}$ let $V_D$ be the image of
polynomials of degree at most $D$, and $L_<$ the integer span of
monomials of price strictly below $T$. Then
$$
V_D/L_<\cong\mathbb Z^3\oplus\mathbb Z/m\mathbb Z, \tag{H13}
$$
and $[\beta]$ generates its torsion. With $N_D=\binom{2m+5}{3}$,
the strict-prefix rank over any field is
$$
\operatorname{rank}_K\{a:w\cdot a<T\}=
\begin{cases}
N_D-16,&m\ne0\text{ in }K,\\
N_D-17,&m=0\text{ in }K .
\end{cases} \tag{H14}
$$
The rank including the whole tier at $T$ is always $N_D-15$.
For every refinement of this positive weight, $\beta$ is standard
if and only if the characteristic divides $m$.

**Proof.** Write $h=xyz^{2m-1}$ and $g=x^2+y^2+z^2$.
The [two-part truncation formula](hypersurface_truncations.md), equation (4), is
$$
\operatorname{in}_{\rm deg}I_{s,r}
=(h)+\sum_{k\ge0}g^k\mathfrak m^{\max(0,s+1-(2m+1)k)}. \tag{H15}
$$
For $k\le m-1$ the minimum displayed degree is at least $4m>D$.
For $k\ge m$, up through $D$ the contribution is exactly $(g^m)\mathfrak m$.
The kernel on polynomials of degree at most $D$ therefore has thirteen
generators:
$$
F_r,\ xg^m,yg^m,zg^m;\qquad
xF_r,yF_r,zF_r,\quad \{ag^m:|a|=2\}. \tag{H16}
$$
The first four have degree $D-1$, the last nine degree $D$.
Their respective highest forms are independent over every field:
$\gcd(g,h)=1$, and $(h)\cap(g^m)$ first occurs in degree $4m+1>D$.
Formula (H15), by induction on degree, proves spanning. Over $\mathbb Z$
these thirteen columns form a primitive submodule, since their rank is
thirteen modulo every prime. Their rational span is the whole kernel,
so they are also the entire integral kernel. The fieldwise dimension
calculation proves that $V_D$ and its inclusion in $A_{\mathbb Z}$
commute with field reduction.

Now $T=D+(4m-1)/(8m)<D+1$. All monomials below degree $D$ are cheaper,
and all above $D$ are more expensive. At degree $D$, the strict price
condition for $x^ay^bz^{D-a-b}$ is $2a+b\ge6$. Exactly twelve
monomials are missing: $2a+b\le5$.

Project the last nine generators of (H16) onto these twelve coordinates.
Their highest forms are $hx,hy,hz$ and $g^m$ times the six quadratics.
There are eight unit pivots, at $(a,b)$ positions
$$
(0,0),(0,1),(1,0),(0,2),(2,0),(1,2),(2,1),(1,1),
$$
using respectively $z^2g^m,yzg^m,xzg^m,y^2g^m,x^2g^m,hy,hx,hz$.
The corresponding eight-by-eight determinant is one.
After integer elimination the ninth column has the sole entry $m$ at
$(1,3)$, the $\beta$ coordinate: subtracting $hz$ from $xyg^m$ cancels
its $(1,1)$ term, and the only remaining missing coordinate is $m\beta$.
The Smith form is $\operatorname{diag}(1^8,m)$ with three free rows.
All lower-degree terms of (H16) were selected in the prefix.
This proves (H13)--(H14).

Including the tier at $T$ removes the missing rows $(0,5),(1,3),(2,1)$.
The matrix then has seven unit pivots and two zero columns over every
field, proving the stated inclusive rank.

Finally
$$
R=xyg^m-zF_r\in I_{s,r}. \tag{H17}
$$
Indeed $xyg^m\equiv(-1)^mxyh^m$ has degree $s+2$ and vanishes.
The $xyz^{2m}$ terms in (H17) cancel. Its unique largest-price
term is $m\beta$, and all other terms are strictly cheaper.
Thus $\beta$ is excluded when $m\ne0$ in the field. When $m=0$,
the kernel on the inclusive prefix has top parts spanned by $hx$
and $xyg^m-hz$; neither has a nonzero $\beta$ coefficient.
The monomial is independent even after all other columns of price
at most $T$, so every refinement includes it. $\square$

The algebra itself remains integrally free: torsion belongs to the quotient
by a selected prefix sublattice. Taking $m=p$ realizes every prime as an
exceptional characteristic in this four-term coefficient-one family.
For general $m$ the theorem identifies exactly its divisors for this prefix,
not every exceptional prime of the entire ideal.

Here $r\asymp\sqrt{s}$, so Theorem 1 bounds the complete fan by $O(n^2)$
in every field. Arbitrary divisibility in coherent ranks and polynomial fan
amplification are distinct mechanisms.

## A coefficient-independent rank law on an equal-price plane

The preceding boundary calculation now extends to an actual two-dimensional
region of the full fan. First there is a more general rank statement.

**Theorem 3.** Over a field $K$, let $Q,H\in K[y,z]$ be coprime
homogeneous forms of degrees $2,r+1$, with $r\ge1$, and put
$$
A=K[x,y,z]/(x^2+Q+xH,\mathfrak m^{s+1}).
$$
For prices $(1+\lambda,1,1)$, $-1<\lambda<r$, define
$$
c_\delta=\begin{cases}1&\delta=0,r+1,\\2&1\le\delta\le r,\end{cases}
\quad
L_{\delta,q}=1+\min\{q,\lfloor(s-\delta-q)/r\rfloor\},
$$
and, for $0\le i\le q$,
$$
\kappa_{q,i}=\max\{(1+\lambda)q-2\lambda i,\ (1-r)q+2ri\}.
\tag{H18}
$$
The rank of **all original monomials of price at most $T$** is
$$
\boxed{\ \sum_{\delta=0}^{r+1}c_\delta
 \sum_{q=0}^{s-\delta}
 \min\big\{L_{\delta,q},\#\{0\le i\le q:\kappa_{q,i}\le T-\delta\}\big\}\ } .
\tag{H19}
$$
Empty sums are zero. In particular this rank depends only on $r,s,\lambda,T$,
not on the coefficients or field, subject to the stated coprimality.
This is a rank of entire price tiers, not a claim about individual
monomial bases within a tier or the full fan.

**Proof.** The regular sequence $Q,H$ gives the classical homogeneous
free-module decomposition
$$
K[y,z]=\bigoplus_j t_jK[Q,H],\qquad
\sum_j t^{\deg t_j}=(1+t)(1+t+\cdots+t^r). \tag{H20}
$$
Generation by lifts of $K[y,z]/(Q,H)$ and equality of Hilbert series prove
freeness. Thus their degrees have multiplicities $c_\delta$ over any field.
Eliminate $Q=-x^2-xH$. This never decreases ordinary degree when $H$
has degree $r+1$. The exact quotient is consequently
$$
\bigoplus_j t_j\,K[x,H]/
 \langle x^aH^b:a+(r+1)b>s-\deg t_j\rangle. \tag{H21}
$$
At equal $y,z$ prices every homogeneous binary space is present in full
or absent from a prefix. The change of basis (H20) therefore preserves
that prefix span. Its columns are $t_jx^aQ^iH^b$, at price
$\delta+(1+\lambda)a+2i+(r+1)b$.

In auxiliary degree $q=a+2i+b$, elimination gives the finite columns
$\xi^b(1+\xi)^i$ modulo $\xi^{L_{\delta,q}}$, with price
$$
\delta+(1+\lambda)q-2\lambda i+(r-\lambda)b,\qquad
0\le i\le\lfloor(q-b)/2\rfloor. \tag{H22}
$$
This labels a finite vector space; no nilpotent element is inverted.
Columns with $b\ge L_{\delta,q}$ may harmlessly be added, since they
are zero in this space.

For $0<\lambda<r$, the allowed $i$ at successive $b$ form intervals
with nondecreasing lower endpoints and nonincreasing upper endpoints.
After factoring the unit $(1+\xi)^{i_0}$, where
$i_0=\max(0,\lceil((1+\lambda)q-(T-\delta))/(2\lambda)\rceil)$,
the prefix spans all polynomials up to its largest degree. A new degree
increases by at most one and is supplied by a monic column. Degree
$d>q/2$ is first attained at $b=2d-q$, at price
$(1-r)q+2rd$; the initial $b=0$ interval gives the other bound in (H18).
Thus its dimension before the cap is exactly the number of indices in
(H18) below the price threshold.

For $-1<\lambda\le0$, every nonempty $i$ interval starts at zero.
The allowed $b$ are consecutive starting at zero, and their upper $i$
endpoints are nonincreasing. The same monic-degree induction gives the
full polynomial space through its maximum degree. For a prescribed degree
$d=i+b$, its cheapest column has $i=\min(d,q-d)$, since $r+\lambda>0$.
Its price is exactly (H18). These prices are nondecreasing in $d$.
This proves the same rank formula, including $\lambda=0$. Reducing the
full polynomial space, or a unit times it, modulo $\xi^L$ caps the rank
at $L$, in every characteristic. Summing proves (H19). $\square$

For $\lambda>0$, the index set in (H19) is the integer interval
$$
\max\!\left(0,\left\lceil\frac{(1+\lambda)q-(T-\delta)}{2\lambda}\right\rceil\right)
\le i\le
\min\!\left(q,\left\lfloor\frac{T-\delta-(1-r)q}{2r}\right\rfloor\right).
\tag{H23}
$$
For $\lambda<0$ both bounds in (H18) instead give upper bounds on $i$;
for $\lambda=0$ the first is the constant $q$.
Hence (H19) uses $O((r+1)(s+1))$ exact arithmetic/floor operations
and constant auxiliary storage. Hypothesis verification, coefficient
encoding and bit complexity are separate. For $\lambda\ge r$,
the elementary $x$-exponent-below-two basis is minimum-price and supplies
the outer formula; the leading monomial is strictly $x^2$ when $\lambda>r$.

## The full fan when the interaction variable is cheapest

**Theorem 4.** For the original family (H1), every generic positive weight
with $w_z<\min(w_x,w_y)$ selects exactly the boundary standard set in (H9),
after exchanging $x,y$ if necessary. These sets are the same over every
field. For $s\ge r+2$ there are exactly $2M_{s,r}+4$ full initial cones
meeting this open region. If $2\le s<r+2$ there are two, and if $s<2$
there is one.

This counts a specified open region of the original fan, not the whole fan.
At its equal-price boundary the rank formula remains valid for complete
tiers; arbitrary tie refinements need not stay inside the region.

**Proof of the rank identity.** Apply Theorem 3 with
$Q=y^2+z^2,H=yz^r$, and first take irrational $0<\lambda<r$.
Write $\mu=\lambda/r$. The list (H18) is the merge of
$$
q+\lambda a\quad(0\le a\le q,\ a\equiv q\bmod2),\qquad
q+rb\quad(1\le b\le q,\ b\equiv q\bmod2).
$$
An item at position $h$ is retained exactly when
$s\ge\delta+q+r(h-1)$. We compare first cutoffs, not merely finite ranks.

For $a>0$ an item of price $D+\lambda a$, $D=\delta+q$, first occurs at
$$
s=D+r\psi(a),\qquad
\psi(a)=\left\lfloor\frac{a+\lfloor\mu a\rfloor}{2}\right\rfloor.
\tag{H24}
$$
Its channel multiplicity is $\min(r+1,D-a+1)$, since
$$
\frac{(1+t)(1+t+\cdots+t^r)}{1-t^2}
=\frac{1+t+\cdots+t^r}{1-t}.
$$
The boundary's $x$-bearing terms at that price are precisely
$x^ay^bz^c$ with $b+c=D-a$ and either $c<r$ or $b=0$.
They have the same multiplicity. Counting cheaper items in their three
streams gives the same first cutoff (H24), both for the spine
and the vertical $x$ arm.

For the remaining items put
$\alpha=(1-\mu)/2$, $\beta=(1-\mu)/(2\mu)$.
The original items at integer price $D$ have labels
$$
D=\delta+2\ell+(r+1)b,\qquad
s=D+r\min(\ell,\lfloor\beta b\rfloor),\qquad \ell,b\ge0.
\tag{H25}
$$
The case $b=0$ includes the zero-$a$ item of the first list.
Boundary terms $y^kz^c$ with $k+c=D$ first occur at
$$
s=D+r\left\lfloor\min\{\alpha k,\beta\lfloor c/r\rfloor\}\right\rfloor.
\tag{H26}
$$
These first-cutoff multisets agree. Indeed for surplus at least $rj$,
$j\ge1$, put $B_j=\lceil2\mu j/(1-\mu)\rceil$. On the original side
$\ell\ge j,b\ge B_j$. Its generating series is
$$
\frac{(1+t)(1+t+\cdots+t^r)\,
 t^{2j+(r+1)B_j}}{(1-t^2)(1-t^{r+1})}
=\frac{t^{2j+(r+1)B_j}}{(1-t)^2}.
$$
The coefficient at $D$ is $\max(0,D-2j-(r+1)B_j+1)$.
On the boundary side $k\ge2j+B_j,c\ge rB_j$, giving exactly this
number. Total multiplicity at $j=0$ is $D+1$ on both sides.
This proves equality of all price-prefix ranks at $(1+\lambda,1,1)$.
The $k$th rank-birth price is the minimum, over independent $k$-sets,
of their maximum price. It is a continuous finite min-max function.
Equality of these sorted birth prices on the dense irrational parameters
therefore gives the complete-tier equality at rational parameters too.

**Proof of the sector statement.** A boundary basis lifts to an original
basis: the binary-degree filtration of the integral monic model has associated
graded quotient (H8). The change-of-basis matrix is block triangular,
with the same diagonal blocks as the boundary basis. These blocks have
unit determinant over $\mathbb Z$, since they are bases over every field.
Thus the original basis is also integral and field independent.

For any nonzero border coefficient let its exponent difference be
$\eta=(a,b,c)$ and put $A=a+b$. The filtration gives $A\ge0$.
Equality of minimum price ranks proves
$A+\lambda a+c\ge0$ at $(1+\lambda,1,1)$: a nonzero coefficient permits
a one-column basis exchange, so a negative value would improve the basis.
For $w=(u+\lambda,u,1)$, $u>1$,
$$
\eta\cdot w=A(u-1)+(A+\lambda a+c). \tag{H27}
$$
This is strict when $A>0$. If $A=0$, it is the strict boundary comparison
at a parameter between successive walls. Every border inequality is
therefore strict throughout that sector. The two outer $x^2,y^2$ cones
are elementary; exchanging $x,y$ completes coverage and (H10) gives the
exact count. No permutation exchanging $z$ with another variable is used.
$\square$

The same three-stream prefix count now computes an **original** rank for
every $w_z\le\min(w_x,w_y)$, including complete tiers at ties. For each
$(k,j)$ residue block it counts how many items in each arithmetic stream
are cheaper than the price, caps their sum at $2L_j$, then sums the blocks.
This costs $O((r+1)(s+1))$ exact arithmetic/floor operations and constant
extra storage. Three lazy sorted streams also emit one minimum-price basis
using $O(n)$ comparisons and output space. At generic weights this is
the unique initial standard set; at ties it chooses one compatible
refinement from inside the region.

Consequently all standard sets and their prefix lattices arising at strict
weights with $z$ cheapest are unimodular and independent of characteristic.
Any initial set outside the listed family, and any exceptional
characteristic at a generic weight, must occur where
$w_z>\min(w_x,w_y)$. The counterexample (H11) lies there, as does the
cyclic family (H12). The main unresolved count for $2\le r<\sqrt{s}$
is now confined to that remaining region. A universal cover by the
equal-price planes has not been proved.

## Literature boundary and verification

Bad primes and monomial-order dependence are established phenomena.
[Abbott--Bigatti--Robbiano, arXiv:1801.06112v3](https://arxiv.org/pdf/1801.06112v3),
Definitions 2.1--2.3, Theorems 2.1--2.3, Proposition 2.1 and Theorems
3.1--3.2 with their displayed proofs were inspected. Theorem 2.1 delegates
its proof to another paper, not re-audited here. Their universal denominator
already controls exceptional primes across orders. This note instead computes
a specified cyclic prefix quotient; their statements alone do not compute
that matrix or this family's fan count.

[Arnold (2003)](https://www3.risc.jku.at/research/theorema/Groebner-Bases-Bibliography/gbbib_files/publication_355.pdf),
Theorems 5.11--5.13 with proofs, were inspected. They relate unlucky primes,
strong-basis leading coefficients and rational denominators, under the stated
homogeneous setup. They are background, not new principles claimed here.

The state and volume mechanisms are classical.
[Babson--Onn--Thomas, Theorem 2.4, its following extension and Corollary 2.5
with proof](https://arxiv.org/pdf/math/0207135) were re-read. They already
permit collections containing the relevant initial sets and use the classical
lattice-vertex bound. Equations (H4)--(H7) insert a narrower slab and the
actual grading index. Exact priority of these specializations and (H13)
remains unassessed; unsuccessful searches are not evidence of novelty.

For the continuation, [Fukuda--Jensen--Thomas, Proposition 2.6,
Lemma 2.10, Corollary 2.11, Lemma 2.13 and Corollary 2.14](https://arxiv.org/pdf/math/0509544)
with proofs were re-read. They supply the classical cone-inequality and
refinement framework, not the rank identity (H24)--(H26).
[Conca--De Negri--Rossi, Definition 2.4, Theorem 4.11 and
Remarks 4.10/4.17](https://arxiv.org/pdf/0705.3767v2) were inspected at
definition/statement level. Their count concerns Cohen--Macaulay initials
of the homogeneous rational-normal-curve ideal in its standard coordinates.
No transfer of that classification to these finite inhomogeneous quotients
has been established. These scoped comparisons leave priority open.

The checker compares the integral model with full Macaulay kernels, boundary
streams with independent Gaussian ranks and exact interval walks, and the
cone obstruction with monic integer membership. It independently computes
the small Smith matrices and both prefix ranks in original monomial columns.
The extension checks four coprime binary-form pairs in several fields,
including negative and zero distinguished-price difference, both rank
implementations, the generalized cutoff identities, and complete original
sectors by monic border membership and both endpoint inequalities.
These finite controls support the uniform proofs, rather than replacing them.

~~~bash
python3 -B scripts/run_higher_interaction_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The remaining counting issue is the simultaneous set of original-cost
constraints for $2\le r<\sqrt{s}$ where $z$ is not the cheapest variable.
The boundary now counts the proved region in Theorem 4, and cannot be
imported as the whole fan. Changed characteristics likewise cannot be
counted as extra cones without proof. A uniform bound or a coherent
large-fan construction in that remaining regime would change the research
decision. Another fixed-cutoff table would not settle it.
