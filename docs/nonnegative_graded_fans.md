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
priority certificate or full mixed-period fan search.
