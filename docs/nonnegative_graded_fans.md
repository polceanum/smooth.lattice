# Linear fans from a nonnegative grading

Interleaving independent tower blocks is not itself a source of additional
initial ideals. The larger obstruction is loss of the grading when a
monomial reaches the tower relation. This note proves a global linear
bound for a much broader class than the completed-block regions: every
finite three-variable quotient homogeneous for a nonnegative integer
grading with at least two positive variable degrees. It also bounds the
local refinement fan at every non-axial nonnegative weight of an arbitrary
finite quotient. A precise prior-art comparison remains incomplete.
These statements do not resolve unrestricted $G_3(n)$ or produce a
polynomial separation above the generic scale.
The signed extension below removes positivity of the grading when all
three variable degrees are nonzero. A zero degree between two opposite
signs requires a different argument and remains outside that extension.

Throughout, $I\subset K[X,Y,Z]$ is proper of colength $n$, $K$ is
arbitrary, and $g(I)$ counts monomial initial ideals for global orders.
Standard sets minimize monomial costs; leading terms have maximum cost.

## The Hilbert budget of a divisor rectangle

Let $(h_\ell)_{\ell\ge0}$ be a finite nonnegative integer sequence of
sum $n$. For coprime positive integers $p,q$, call $(a,b)\in\mathbb N_0^2$
admissible when

$$
\#\{(i,j):0\le i\le a,\ 0\le j\le b,\ pi+qj=\ell\}
 \le h_\ell\quad\hbox{for every }\ell.                    \tag{G1}
$$

**Lemma 1.** There are at most $2n$ admissible rectangles, uniformly in
$p,q$. If one variable has degree zero and the other degree $p>0$,
there are at most $\sum_{j\ge0}h_{pj}\le n$ admissible rectangles.

**Proof.** First take $pa\ge qb$ and write $b=r+pk$, $0\le r<p$.
Then $a\ge qk$. At degree $\ell=pa+qr$, the rectangle contains

$$
(a-qj,r+pj),\qquad 0\le j\le k.
$$

Hence $k<h_\ell$. The map $(a,b)\mapsto(\ell,k)$ is injective:
$\ell\bmod p$ determines $r$ since $q$ is invertible modulo $p$,
then determines $a$, and $k$ determines $b$. This half has at most
$\sum h_\ell=n$ elements. For $pa<qb$, interchange the variables
and use $\ell=qb+p(a\bmod q)$ and $k=\lfloor a/q\rfloor$.
This gives a second injection into the same $n$ slots.

For degrees $(0,p)$, put
$d_b=\min_{0\le j\le b}h_{pj}$. The allowed $a$ satisfy
$a+1\le d_b$, so their number is $\sum_b d_b\le\sum_jh_{pj}$.
$\square$

For noncoprime positive degrees, divide by their gcd and use the
subsequence of $h$ on multiples of that gcd. The bound remains $2n$.
The lemma concerns necessary divisor budgets, not a claim that every
admissible rectangle extends to a quotient basis or coherent staircase.

## Why profile edges have such a rectangle

Use the established monomial residue matroid on a finite universal ground
containing all standard sets. A grading splits it as a direct sum by
degree. Standard sets have exactly $h_\ell=\dim_K(K[X,Y,Z]/I)_\ell$
monomials in degree $\ell$, independently of the global monomial order.

At adjacent exposed standard bases, choose a strong symmetric basis
exchange between their differences at the common minimizing objective.
Both exchanged sets are bases. Optimality in the two directions forces
the exchanged monomials to have equal cost; the direct sum forces equal
degree. The two exchanged monomials belong to the actual endpoint
standard sets. The exchanged profile stays on the same profile edge,
so its nonzero exponent difference is parallel to that edge.

If the difference has positive coordinates $a,b$ in two variables,
one endpoint standard set contains every divisor in that $a$-by-$b$
rectangle. Dividing the difference by its gcd only shrinks the rectangle.
Lemma 1 therefore applies to the positive part of the primitive edge
direction. This step uses genuine standard sets, rather than supposing
that intermediate matroid bases are monomial downsets.

The matroid, greedy and state-profile dictionary is established in
[Babson--Onn--Thomas, Proposition 2.2 and Theorem 2.4](https://arxiv.org/pdf/math/0207135).
Symmetric exchange and the planar convex-polygon facts used below are
classical. The Hilbert-budget count and its exact consequences are the
statements whose priority remains to assess.

## All three variable degrees positive

**Theorem 2.** If $I$ is homogeneous for any positive integer grading
$(d_X,d_Y,d_Z)$, then

$$
g(I)\le12n+6.                                          \tag{G2}
$$

**Proof.** Every basis profile lies in the plane
$d_X\mu_X+d_Y\mu_Y+d_Z\mu_Z=\sum_\ell\ell h_\ell$.
Every vertex can be exposed by a positive objective: add a sufficiently
large multiple of $(d_X,d_Y,d_Z)$ to any exposing vector. Thus the
profile polygon's vertices correspond to standard sets.

An edge direction with three nonzero coordinates has, up to sign, two
positive coordinates and one negative coordinate. Choose the negative
coordinate in three ways. The two positive coordinates form an
admissible rectangle, of which there are at most $2n$ by Lemma 1.
The grading determines the remaining coordinate uniquely. There are
therefore at most $6n$ unoriented full-support edge directions.
With two nonzero coordinates, the grading fixes at most one primitive
direction per coordinate pair, adding at most three.

A convex polygon has at most two edges parallel to any unoriented line.
This gives (G2). A segment or point satisfies it separately. $\square$

The constant is not asserted optimal, and no matching lower bound is
proved here for the strictly positive subclass alone. The bound is
uniform in the numerical grading weights, characteristic, radicality,
generators and coefficient sizes. It does not rely on a lattice-volume
estimate or replace $n$ by a presentation degree.

## Signed degrees: the extension and its exact boundary

Let $h:\mathbb Z\to\mathbb N_0$ have finite support and total $n$.
For opposite pair degrees $p,-q$, normalize $p,q>0$ to be coprime.
An admissible rectangle has at least

$$
1+\min(\lfloor a/q\rfloor,\lfloor b/p\rfloor)
$$

divisors in its corner degree $pa-qb$. The map

$$
(a,b)\longmapsto
\left(pa-qb,\min(\lfloor a/q\rfloor,\lfloor b/p\rfloor)\right)       \tag{G2a}
$$

is injective. Two corners of the same degree differ by $(qt,pt)$;
their two floor quotients, and hence their minimum, both change by
$t$. Each image is a slot of the Hilbert function, so there are at
most $n$ admissible rectangles. Common pair-degree factors and an
overall sign are handled by restricting or reversing the degree index.
Thus any pair of nonzero degrees has budget at most $2n$, even when
the signs differ.

**Theorem 2a.** If $I$ is homogeneous for an integer grading
$d=(d_X,d_Y,d_Z)$ with **all three coordinates nonzero**, then

$$
g(I)\le12n+7.                                          \tag{G2b}
$$

For degrees of one sign, Theorem 2 applies after an overall sign change.
For mixed signs, project the positive cost orthant to the two-dimensional
objective space modulo $\mathbb R d$. Its image is a pointed convex
cone. Indeed a line in the image would lift to two nonzero nonnegative
vectors whose sum belongs to $\mathbb R d$, impossible for mixed $d$.
The positive exposed part of the profile polygon is therefore a connected
chain. Its vertices are actual initial standard sets; other polygon
vertices need not be, and must not be included in the divisor argument.

Each edge of that chain has a symmetric optimal exchange between its
actual endpoint standard sets. As before, its primitive positive pair
is an admissible rectangle. The remaining coordinate is determined
because its grading degree is nonzero. The signed rectangle bound
gives at most $6n+3$ unoriented directions. A convex polygon has at
most two edges per direction, so the chain has at most $12n+7$
vertices. Degenerate profile dimension is separate and satisfies the
same bound. This proves the claim over every field.

The hypothesis about the third degree is necessary **for this encoding
argument**. It is not a proved boundary of the linear conclusion.
For grading $(0,1,-1)$, positive pairs $(c,c)$ in the last two
coordinates do not determine the first coordinate of an exchange.
There is an actual all-field family demonstrating that ambiguity:

$$
I_{a,c}=(XY,XZ,X^{a+1},Y^{c+1},Z^{c+1},X^a-Y^cZ^c),
\quad a,c\ge1.                                         \tag{G2c}
$$

Before the last relation the length is $a+(c+1)^2$; the relation
identifies the two nonzero socle monomials and reduces it by one.
Thus $n=a+(c+1)^2-1$. There are exactly two initial ideals, with
wall direction $(a,-c,-c)$, reduced by its gcd, and wall
$a w_X=c(w_Y+w_Z)$. Both choices are strict positive-weight cones.
For fixed $c$ and coprime varying $a$, the same positive pair $(c,c)$
therefore occurs with different primitive directions. The example
invalidates an automatic extension of the encoding; it does not
disprove a linear bound for $(0,p,-q)$.

There is still useful structure in this excluded case. In any finite
algebra graded by $(0,p,-q)$, $Y,Z$ are nilpotent. Put $T=Y^qZ^p$.
Every grade's nonnegative $Y,Z$ solutions form one chain
$(b_\ell+qj,c_\ell+pj)$, $j\ge0$. Its monomial columns are therefore
a cyclic finite module over $K[X,T]$. Different grades are independent,
but their costs use the same positive parameter
$(q w_Y+p w_Z)/w_X$. Applying a planar bound separately gives the
earlier $O(n\log n)$ bound, not the new linear one. Bounding all
these modules' changes with one shared budget remains a possible
extension; no such budget is asserted here.

## One zero degree and two positive degrees

Normalize the grading to $(0,p,q)$ with $p,q$ coprime and positive.
Put $D=h_0$ and

$$
d^Y_b=\min_{0\le j\le b}h_{pj},\qquad
d^Z_c=\min_{0\le j\le c}h_{qj},\qquad
M=\left\lfloor\frac{n-1}{\max(p,q)}\right\rfloor.
$$

**Theorem 3.** If $M=0$, then $g(I)=1$. Otherwise every active
finite parameter belongs to the coefficient-independent set

$$
\mathcal T(h)=\{0\}
 \cup\{a/c:1\le c\le M,\ 1\le a<d^Y_{qc}\}
 \cup\{-a/c:1\le c\le M,\ 1\le a<d^Z_{pc}\}.             \tag{G3}
$$

In particular

$$
g(I)\le 2+\sum_{c=1}^{M}\bigl((d^Y_{qc}-1)_+
                         +(d^Z_{pc}-1)_+\bigr)
 \le2+2(n-D).                                          \tag{G4}
$$

The finite list of positive integer weights covering its intervals can
be chosen with height at most $4n$.

**Proof.** Write a positive cost as $(\alpha,\beta,\gamma)$.
The fixed grading moment removes one cost coordinate, leaving only
$\rho=(p\gamma-q\beta)/\alpha$. An exchanged pair of equal degree
has exponent difference $(a,qc,-pc)$, after orientation when its
$X$ coordinate and $Y$ coordinate are positive. Its tie requires
$\rho=a/c$. The endpoint staircase contains $X^iY^j$ for
$i\le a,j\le qc$, so $a<d^Y_{qc}$. The opposite sign of $\rho$
uses an $X,Z$ rectangle and $a<d^Z_{pc}$. A zero $X$ difference
gives $\rho=0$. Both endpoints contain respectively at least the
pure powers $Y^{qc}$ and $Z^{pc}$; hence $\max(p,q)c\le n-1$.
This also proves that there is no active parameter if $M=0$.

There is one interval state between successive elements of (G3),
including the two extremes; a boundary monomial refinement is a limit
from one of its sides. This proves the first bound. Each running-minimum
sum is at most the corresponding sum of positive-degree Hilbert values,
which is at most $n-D$, proving the second.

For the height claim, use mediants between the distinct reduced
parameters, and extremes $\rho=\pm D$. A sample $A/B$ has
$|A|\le2D$ and $B\le2(n-1)$. If $M\ge1$, then $p,q\le n-1$.
Choose $\beta_0\in\{1,\ldots,p\}$ satisfying
$q\beta_0\equiv-A\pmod p$, and put
$\gamma_0=(A+q\beta_0)/p$. If $\gamma_0>0$, use
$(\alpha,\beta,\gamma)=(B,\beta_0,\gamma_0)$.
Otherwise add the smallest nonnegative multiple of $(p,q)$ to
$(\beta_0,\gamma_0)$ making $\gamma>0$. Then $\gamma\le q$ and
$\beta<2p+|A|/q\le4n$. In the first case
$\gamma_0\le q+|A|/p\le3n$. Thus all weights are positive and
have height at most $4n$, and give exactly $\rho=A/B$.
If $M=0$, the single weight $(1,1,1)$ suffices. $\square$

This does not assert an $O(n)$ common list simultaneously for every
Hilbert function of length $n$. The list depends on the grading and
the fixed Hilbert function, but not on the particular ideal within
that class.

The linear growth and linear certificate-height orders are both sharp
for this subclass. For $D\ge1$ use the
[explicit length-$2D$ ribbon](cyclotomic_price_regions.md)
$F=\mathbb Q[x]/\prod_{i=0}^{D-1}(x-2^i)$, $T^2=0$,
$X=x,Y=T,Z=x^DT$. It has $D+1$ initial ideals by the ordinary
Vandermonde calculation. Selecting its all-$Y$ degree-one basis
requires $w_Z-w_Y>(D-1)w_X$, hence positive integer height at least
$D+1$. This occurs at arbitrarily large lengths.

Odd lengths can be included without a grading-breaking relation:
give the degree-zero component the polynomial
$x\prod_{i=0}^{D-1}(x-2^i)$ and keep the degree-one component
annihilated by $\prod_{i=0}^{D-1}(x-2^i)$, with all products of
positive-degree elements zero. The same $D+1$ states occur at length
$2D+1$. Thus the maximum over finite quotients admitting a nonnegative
grading with at least two positive coordinates has order $\Theta(n)$
pointwise, in the all-fields maximum, and the support-two subclass
has optimal worst-case integer weight height $\Theta(n)$.

Two zero variable degrees are a real boundary. An arbitrary planar
finite ideal $J(X,Y)$ embeds as $(J,Z)$ with grading $(0,0,1)$;
its known generic count has order $n\log n$. The earlier
[rank-two profile theorem](profile_lattice_bound.md) gives the matching
upper order. No linear theorem for all integer gradings, including
mixed signs or two zeros, is claimed here.

## A local consequence for unrestricted quotients

Let $J$ be any colength-$n$ ideal and let $w\ge0$ be a nonzero
rational weight with at least two positive entries. The associated
graded quotient has the same length, and $\operatorname{in}_w(J)$
is homogeneous for an integer multiple of $w$. Standard initial-form
refinement gives

$$
\operatorname{in}_{w+\varepsilon v}(J)
 =\operatorname{in}_v(\operatorname{in}_w(J))
\quad\hbox{for sufficiently small }\varepsilon>0.         \tag{G5}
$$

For an interior $w$, adding a large multiple of $w$ to $v$ makes
it positive without changing its initial form on the graded ideal.
For a boundary $w$ of support two, retain $v>0$ in the zero coordinate
and apply the same adjustment in the other coordinates. Theorems 2--3
give **$O(n)$ distinct local positive refinements**. Equivalently,
every non-axial nonnegative ray of the positive Gröbner fan has linear
local complexity, uniformly over fields and ray heights. The three
coordinate-axis rays retain the broader $O(n\log n)$ bound.

This is a local statement for arbitrary ideals, not a bound on how
many rays occur or on the full fan. It strengthens the earlier general
rank-two $O(n\log n)$ face estimate for precisely these positive
grading faces. Classical planar support-oracle enumeration uses $O(n)$
greedy calls for the local fan; polynomial arithmetic and coefficient
bit complexity are separate. No unrestricted $G_3$ exponent changes
merely by multiplying this local bound by an unproved number of rays.

## Tower deformations: where extra states must enter

Let $A=F[T]/g$, $F=K[x]/f$, with monic degrees $D,R$ and units
$u,v\in F$, using $X=x,Y=uT,Z=vT$. Put $A_0=F[T]/T^R$.
For **all generic** positive weights the standard set $B_0(\rho)$ of $A_0$
is the union of the greedy slot bases in (P6) of the
[tower note](cyclotomic_price_regions.md), without condition (P5).
The blocks may interleave, but their vector spaces are independent,
and every higher block is zero. This also follows from Theorem 3;
the explicit tower formula gives its smaller, degree-specific list.

The same $B_0(\rho)$ is a vector basis of $A$, since it uses only
$T$ degrees below $R$. Its minimal excluded heads with $Y,Z$ degree
below $R$ have exactly the same normal forms and strict inequalities
in $A$ and $A_0$. The only further inequalities are for minimal
heads of degree $R$. Such a head must have $X$ exponent zero;
otherwise its $X$ predecessor would already be excluded. There are
at most $R+1$ of them, among $Y^bZ^{R-b}$.

Consequently, at a strict slot parameter, $B_0(\rho)$ remains an
initial standard set of $A$ if and only if these heads have strictly cheaper
normal forms in $B_0(\rho)$. This is the standard finite border
criterion, applied to this specific independent-block presentation.
Within a fixed slot phase it gives an intersection of linear
inequalities. It creates at most one surviving cone per $B_0$;
the number of initial sets confined to degrees below $R$ stays $O(n)$.

Conversely, **any** initial standard set of $A$ using only $T$ degrees
below $R$ must be $B_0(\rho)$. Its $DR$ independent vectors fill
all $R$ independent $D$-dimensional blocks; minimizing within them
gives precisely that greedy set. Every additional initial set must
contain a monomial of $Y,Z$ degree at least $R$.

For the original mixed-period construction, the distinction is between
$W^{t-1}=0$ and $S_t(W)=0$. Its still-open growth must involve the
columns reaching the latter relation, not just an interleaving of
incomplete but independent lower blocks. This is a necessary mechanism,
not a count of the resulting states or a lower bound. The next task is
to bound or construct coherent exchanges involving these columns in
the unbounded-price intermediate region and for unbalanced periods.

The original smooth-number observations order exponent vectors by positive
costs such as logarithms of fixed primes. Standard sets impose algebraic
dependencies on that same exponent order. The linear grading theorem
identifies a broad kind of dependency that restricts the ordering family;
the mixed-period relation breaks this grading. It does not yield a prime
distribution theorem or a new certified unranking range.

## An exact prefix formula for the original mixed-period algebra

There is a direct way to study the columns that reach the relation.
Keep the mixed-period notation
$F=K[U,V]/(S_r(U),S_s(V))$, $D=(r-1)(s-1)$ and
$A=F[W]/S_t(W)$, with $X=UV,Y=UW,Z=VW$.
For any finite monomial list $E$, group its elements by
$j=(b+c)\bmod t$ and set

$$
V_j(E)=\operatorname{span}_K
 \{U^{a+b}V^{a+c}:(a,b,c)\in E,\ b+c\equiv j\pmod t\}
 \subseteq F.
$$

Then, over every field,

$$
\operatorname{rank}_A(E)=\sum_{j=0}^{t-1}\dim V_j(E)
                         -\dim\bigcap_{j=0}^{t-1}V_j(E).  \tag{G7}
$$

This is an elementary tensor/circuit rank identity, not a new general
matroid theorem. The kernel of $F^t\to A$, $(v_j)\mapsto\sum v_jW^j$,
is exactly the diagonal copy of $F$: $W^{t-1}=-\sum_{j<t-1}W^j$.
Intersecting that kernel with $\bigoplus V_j$ proves (G7), including
repeated labels and characteristics dividing a period.

The known cographic model makes this more concrete. Let $G=K_{r,s}$,
and let $S_j\subseteq E(G)$ be the distinct coefficient labels
$((a+b)\bmod r,(a+c)\bmod s)$ occurring in group $j$.
Write $\mathcal Z(H)$ for the cycle space of a spanning subgraph
$H$, over $K$, and $c(H)$ for its number of components including
isolated vertices. The established cographic representation identifies
$V_j^\perp$ with $\mathcal Z(G\setminus S_j)$. Hence

$$
\begin{split}
\dim V_j&=|S_j|+1-c(G\setminus S_j),\\
\dim\bigcap_j V_j
 &=D-\dim\left(\sum_j\mathcal Z(G\setminus S_j)\right).    \tag{G8}
\end{split}
$$

The sum in (G8) is a sum of vector spaces, **not** the cycle space
of the union of the graphs. Even two spanning trees can have a union
with cycles while both individual cycle spaces are zero. In $K_{3,3}$,
the tree made of row zero and column zero and the tree made of row one
and column one have a union of eight edges and cycle rank three;
their cycle-space sum has rank zero. Substituting connectivity of the
union would therefore give the wrong rank. This identifies an exact
reason the two-factor connectivity reduction does not extend directly.

For positive prices $w$ use the actual finite cost prefix
$E_C=\{(a,b,c)\ge0:w\cdot(a,b,c)\le C\}$. Formulas (G7)--(G8)
then give its exact rank and, with a tie refinement, its greedy rank
increments. Before $C<(t-1)\min(w_Y,w_Z)$, the last group is empty,
so the intersection is zero and the lower blocks are independent.
Once all groups occur, the growth of their intersection becomes part
of the problem. A new column may still fill an uncompleted lower
span before an intersection grows; no equivalence between those
events is assumed.

The next uniform question is how these specific cost-defined subspaces
and cycle-space sums vary with the three prices. Formula (G7) alone
does not bound the number of their changes, give a new rank algorithm
complexity, or establish a coherent lower family. It turns the remaining
mixed-period question into a defined prefix-intersection problem without
replacing the original extremal endpoint.

### A wider region with at most one complete period

Put $m=\min(w_Y,w_Z)$ and $\alpha=w_X$. If

$$
t m>(D-1)\alpha,                                       \tag{G9}
$$

every standard monomial has $a<D$ and $b+c<t$. The monic degree-$D$
relation of $X$ excludes $a\ge D$ for every positive weight. For the
second assertion, fix a period label $j$. The $D$ monomials $X^a$ times the $j$th power of
the cheaper of $Y,Z$ span its entire coefficient space $F$. They all
cost at most $jm+(D-1)\alpha<(j+t)m$. Every monomial of the same
period label with $b+c\ge j+t$ costs at least $(j+t)m$, so its
coefficient is already in the span of earlier columns in that group.
It cannot be selected by the global greedy process, even when some
earlier columns were rejected through the diagonal relation.

More precisely, let $E_0(\rho)$ contain the $D$ coefficient-greedy
slots in each of the $t$ groups $j=0,\ldots,t-1$, using $a<D$,
$c\le j$ and $\rho=(w_Z-w_Y)/w_X$. Under (G9), every initial
standard set is a subset of this $tD$-element set, and omits exactly
$D$ elements. There are only $O(tD)=O(n)$ possible auxiliary sets
$E_0$; this does not bound how many rank-$D$ omissions can be coherent
within each one.

Condition (G9) is weaker by a factor $t$ than the earlier completed-
block condition $m>(D-1)\alpha$. It permits degree $t-1$, so it does
not exclude the new columns in (G7). It excludes additional complete
periods, rather than assuming they are harmless.

The omissions also give a stronger count than the earlier rectangular
volume estimate. In fact

$$
g_{\text{(G9)}}(A)=O(n^2),                              \tag{G10}
$$

uniformly over the periods, their relative sizes and the field. This
supersedes the valid but weaker $O(n^2\sqrt t)$ estimate. It excludes
this whole price region as a source of a polynomial separation above
the generic scale. An exact classification of the omissions is not
needed for that exclusion.

### The profile-deficit argument

**Lemma 4.** Let a rank-$n$ represented matroid have ground vectors
labelled by monomials $X^aY^bZ^c$ with $0\le a<D$ and $b+c<t$.
Suppose the columns at each exact level $j=b+c$ have rank at most
$D$, and $n=DR$ with $1\le R<t$ and $D\ge2$. Put $s=t-R$.
The number of profiles exposed by generic linear monomial costs is

$$
O\!\left(1+n(t-1)\sqrt{D(D-1)s}\right).                \tag{G11}
$$

For every independent basis $B$, write $m_j=|\{e\in B:b+c=j\}|$.
Its deficits $\delta_j=D-m_j$ are nonnegative and sum to $Ds$.
Consequently

$$
\sum_{e\in B}(b+c)
 =S-\sum_{j=0}^{t-1}j\delta_j,
\qquad S=\frac{Dt(t-1)}2.
$$

The last profile coordinate therefore varies only in
$[S-Ds(t-1),S]$. The integral coordinate change
$(a,b,c)\mapsto(a,c,b+c)$ has determinant $-1$. In these
coordinates the basis-profile polytope lies in an integer translate
of a box of side lengths

$$
n(D-1),\qquad n(t-1),\qquad Ds(t-1).
$$

Its volume is at most $n^2D(D-1)s(t-1)^2$. The established
three-dimensional lattice-vertex estimate gives (G11). Adding a unit
cube handles lower-dimensional polytopes, changes these positive
side bounds by at most a constant factor, and refines the normal fan.
This is an elementary rank-deficit budget followed by a classical
lattice-volume consequence; no new lattice-vertex or matroid-quotient
theorem is claimed.

Apply the lemma to **all** rank-$n$ bases in the finite ground
$a<D,b+c<t$ for $A$. The level-$j$ vectors lie in $FW^j$, of
dimension at most $D$. Every initial standard set in (G9) is a unique
minimum basis in this ground for a suitable generic perturbation, so
its profile is a vertex of that same polytope. Here $R=t-1$, $s=1$,
and the box volume is

$$
n^2D(D-1)(t-1)^2=n^4(1-1/D)\le n^4.
$$

This proves (G10) when $D\ge2$. If $D=1$, $X$ is scalar and
$Y,Z$ are scalar multiples of $W$, giving at most two standard sets.
The thin-degree bound holds for arbitrary independent bases in this
ground, not merely for the previously identified auxiliary slots.
Both ends of its degree interval are attained by some monomial bases:
take all $D$ columns $X^aY^j$ in every level except $j=t-1$ or
except $j=0$. Such bases need not be monomial downsets; this endpoint
observation is not a lower bound for positive-weight initial ideals.

### A general periodic-factor consequence

**Theorem 5.** Let $F=K[x]/f$ have dimension $D$, let $u,v\in F$
be units, and let $g\in F[W]$ be a monic divisor of $W^t-1$ of
degree $R=t-s$, with $1\le s<t$. Set
$A=F[W]/g$, $X=x,Y=uW,Z=vW$, and $n=DR$. In (G9),

$$
g_{\text{(G9)}}(A)
 =O\!\left(1+n(t-1)\sqrt{D(D-1)s}\right)
 =O\!\left(1+n^2\frac{t-1}{R}\sqrt s\right)             \tag{G12}
$$

for $D\ge2$. In particular the bound is $O(n^2\sqrt s)$ if
$R\ge t/2$, and $O_s(n^2)$ for every fixed defect $s$. No
separability, radicality, constant-coefficient assumption on $g$, or
restriction on the characteristic is needed.

Indeed $W^t=1$, so the same cheaper-generator spanning argument
used for (G9) excludes all levels at least $t$. The monic $f$ excludes
$a\ge D$. Every level has rank at most $D$, and the monic $g$
gives length $DR$. Lemma 4 applies. The case $D=1$ again has at
most two standard sets. The original $g=S_t(W)$ has $s=1$ and
therefore satisfies the uniform quadratic bound.

The profile argument actually counts a larger family: all initial
standard sets contained in $a<D,b+c<t$, at **any** positive price.
Each is still a unique minimum basis in that same finite ground, so
the proof does not need (G9) once containment is known. The price
condition is only a sufficient test for containment. The same
observation gives (G12) for every one-period-contained initial set
of the general periodic-factor algebra.

For the original three-period construction, apply this argument with
each of $X,Y,Z$ as the coefficient coordinate. The union of all
initial sets satisfying at least one of

$$
\max_{(a,b,c)\in B}(a+b)<r,\qquad
\max_{(a,b,c)\in B}(a+c)<s,\qquad
\max_{(a,b,c)\in B}(b+c)<t
$$

has only $O(n^2)$ members, uniformly even for unbalanced periods.
Any family giving a polynomial separation must therefore have
asymptotically more than quadratic many bases reaching **all three**
period thresholds. Different monomials may witness the three maxima;
no single monomial is required to reach them simultaneously.

These exclusions do not count the remaining bases. Outside (G9),
actual bases can reach another period; the $(3,5,7)$ controls with
sufficiently cheap $Y$ include $Y^{11}$ and cannot use the $X$-coordinate
one-period ground. They might still fit a different ground, so these
controls alone do not establish simultaneous crossing of all three
thresholds. The remaining candidate mechanism is the variation of
cost-prefix intersections across multiple periods in every coordinate
presentation. No global separation, unrestricted exponent or new
rank-algorithm complexity is asserted.

## Prior-art boundary and verification

The graded/state dictionary, polygon edge facts, symmetric exchange,
flat initial degeneration, initial-of-initial refinement and border
criterion are established inputs. In particular:

* [Altmann--Sturmfels, Proposition 3 and Theorem 8](https://arxiv.org/pdf/math/0209152)
  give the edge-ideal/direction framework. Those primary passages and
  the beginning of Section 3 were inspected; the full paper was not
  audited for an equivalent numerical bound.
* [Kreuzer--Robbiano, Proposition 2.2 and Theorem 2.4](https://arxiv.org/pdf/0710.2641)
  concern Gröbner and border-basis deformations. Those passages and
  the homogeneous-border argument were inspected. They do not make
  arbitrary border bases Gröbner bases, and no monotonicity of fan
  counts under flat deformation is assumed here.
* [Haiman--Sturmfels, Corollary 3.4, Theorem 3.6 and the supportive-set
  definitions](https://arxiv.org/pdf/math/0201271) were inspected. Their
  finite supportive degree sets give the multigraded-Hilbert-scheme
  construction. This is a relevant framework comparison; exact coverage
  of the quantitative rectangle and fan bounds by the full paper or its
  subsequent literature remains to audit. An unsuccessful search is not
  evidence of novelty.

The $2n$ weighted-rectangle injection, the resulting global linear
grading bounds, the support-two certificate-height claim and the local
fan consequence have complete arguments above. Their exact priority
and the strongest-current unrestricted bound remain unresolved.
No new paper, completion claim or publication assessment follows yet.
All earlier manuscripts and the full main objective are retained.

### What the closest inspected statements actually compare

| Primary statement | Established content | Exact remaining comparison |
| --- | --- | --- |
| Babson--Onn--Thomas, Proposition 2.2, Theorem 2.4, Corollary 2.5 | The minimum-basis/state dictionary and a polynomial bound for arbitrary finite quotients. | The dictionary is an input. The stated general bound does not by itself give the Hilbert-sensitive linear bound. |
| [Sturmfels, *Gröbner Bases and Convex Polytopes*, Theorem 2.5 and Lemma 2.6](https://math.berkeley.edu/~bernd/GBCP.pdf) | A state polytope realizes the fan of a positively graded ideal; faces correspond to initial ideals. | The author-hosted indexed primary excerpt was inspected, not the full scanned book. This supplies the classical fan framework, not a checked quantitative-priority conclusion. |
| Altmann--Sturmfels, Proposition 3 and Theorem 8 | Edge ideals and uniqueness of an edge direction between prescribed monomial ideals. | Uniqueness for two endpoints does not count directions over a Hilbert class. The full subsequent quantitative literature remains unassessed. |
| Haiman--Sturmfels, Proposition 3.2, Corollary 3.4, Theorem 3.6 | Finite supportive degree sets and embeddings/constructions of multigraded Hilbert schemes. | A finite supportive set is not the same assertion as a linear rectangle or edge budget. Exact derivability of that numerical budget from other results still needs checking. |
| Kreuzer--Robbiano, Proposition 2.2, Lemma 2.3, Theorem 2.4 | Homogeneous borders, positive weight realizations and conditional flat deformations. | These justify established deformation tools. They do not justify assuming that a flat deformation preserves or bounds the number of fan cones. |
| [Frels--Schmitz, Notation 2.10 and Theorems 2.12--2.13](https://arxiv.org/pdf/1303.3461) | A characteristic-zero lower bound linear in generating degree for generic initial ideals, with a generic-family extension. | Their parameter is generating degree and their count uses generic coordinates. It is not a uniform upper bound in finite quotient length. |

For the normalization check, their explicit Notation 2.10 family is
$J_d=(X^aY^{d-a},Z^d:0\le a\le d-1)$. Every $X^k$ survives,
so that particular family has infinite colength. This does **not** say
that their entire generic-family theorem excludes Artinian ideals.
One can instead truncate at $(X,Y,Z)^{d+1}$. The degree-$d$ component,
where their lower-bound mechanism is proved, remains unchanged, while
the finite length becomes

$$
N_d=\binom{d+2}{3}+\binom{d+2}{2}-(d+1)
   =\frac{d(d+1)(d+5)}6.                                \tag{G6}
$$

Thus a lower bound linear in $d$ remains consistent with a linear
upper bound in $N_d$; it does not supply a matching lower bound for
our strictly positive grading class. Definitions of generic coordinates/fans
and the main degree-$d$ mechanism were inspected; the determinant appendix
was not independently verified.

The formerly unidentified [NSF PDF](https://par.nsf.gov/servlets/purl/10393536)
is Silversmith, *The matroid stratification of the Hilbert scheme of
points on* $\mathbb P^1$, Manuscripta Mathematica 167 (2022),
173--195, [DOI](https://doi.org/10.1007/s00229-021-01280-z).
In [arXiv:1911.03569v3](https://arxiv.org/pdf/1911.03569v3),
Theorem 3.7 describes strata using Schur polynomials; Proposition 5.5
expresses adding a monomial ideal as looped contraction; Theorem 5.11
describes specified planar $T$-graph edges. Those statements and their
relevant definitions were inspected. They concern matroid strata and
edge schemes, not the numerical three-price count proved here.
This closes the bibliographic lead, not the full quantitative-priority
audit; no equivalence with, or independence from, all subsequent work
has been established.

```bash
python3 -B scripts/run_graded_fan_checks.py
```

The controls test rectangle injections independently, construct actual
finite homogeneous quotients from multiplication-closed relation spaces,
enumerate complete small comparison fans independently of the Hilbert
candidate set, verify strict borders and exchanged endpoint rectangles,
test even/odd sharpness and integer weights, replay local refinements
of nonhomogeneous quotients, and compare nilpotent versus geometric-sum
tower relations. They are finite supporting checks, not a uniform proof,
priority certificate or full mixed-period fan search. Additional controls
use actual shared facets for signed full-support gradings and the
zero/opposite-sign boundary family, compare ordinary polynomial ranks
with both subspace-intersection and graph-cycle-sum calculations on full
finite cost balls, and verify the one-period range and exact omitted
slots at specified prices. Characteristics dividing a period are included.
The profile continuation additionally checks the degree-deficit identity,
general monic factors of $W^t-1$ over nonreduced coefficient algebras,
a factor with coefficients in $F$ rather than just $K$, every rank-sized
subset of three small one-period grounds, and actual additional-period
bases outside the price condition. These verify the mechanism and its
boundary, not an asymptotic constant or a global fan count.
