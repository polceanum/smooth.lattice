# Rational-function curves: generic fans and a uniform quadratic bound

Generic finite samples on the curves $XY=1$, $P(X)Z=Q(X)$ have
$O(n^2\log(n+1))$ monomial initial ideals, uniformly in the degree,
coefficients, root multiplicities and characteristic of the reciprocal
polynomial $P$, provided $P(0)\ne0$, $\gcd(P,Q)=1$ and
$\deg Q\le\deg P$. An exact strip-selection rule and a small moment budget
on the sensitive parts of carrier walls give the bound. Thus this whole
family cannot amplify the generic corner-cut scale by a polynomial factor.
A further exact integer-cost oracle supplies rank and random access
without quotient matrices. Arbitrary special samples are not covered.

This is a whole-family feasibility result, not a resolution of $G_3(n)$,
a polynomial lower bound, or a novelty certificate. The strongest current
unrestricted comparison and the original full research goal remain open.

## The actual quotient and its infinite carrier

Let $F$ be any field, $k,n\ge1$, and take algebraically independent
$t_1,\ldots,t_n$. Over $K=F(t_1,\ldots,t_n)$ put

\[
 f(T)=\prod_{j=1}^n(T-t_j),\qquad
 I_{n,k}=(XY-1,(X-1)^kZ-1,f(X)).                         \tag{H1}
\]

Its length is $n$, with $X=T,Y=T^{-1},Z=(T-1)^{-k}$ in $K[T]/f$.
All denominators are units. It is the ideal of $n$ distinct $K$-rational
points, not an arbitrarily chosen matroid. The carrier is

\[
 J_k=(XY-1,(X-1)^kZ-1),\qquad
 K[X,Y,Z]/J_k=K[T,T^{-1},(T-1)^{-1}].                    \tag{H2}
\]

For the last equality, $(T-1)^{-1}=(T-1)^{k-1}Z$.
The coordinate ring is independent of $k$; its chosen three generators
and their positive cost order are not.

Write $[T^j]=X^j$ for $j\ge0$ and $[T^j]=Y^{-j}$ otherwise.
For $0\le i\le k$ define the actual carrier relation

\[
 G_i=Z\sum_{j=0}^k(-1)^{k-j}\binom kj[T^{j-i}]-Y^i.
                                                               \tag{H3}
\]

It is obtained by multiplying $(X-1)^kZ-1$ by $Y^i$ and reducing
$XY=1$. Consequently it belongs to $J_k$ integrally, in every
characteristic. These are elementary Laurent identities, not new
binomial-coefficient or localization theorems.

## Exactly $k$ carrier cones, over every field

Let $(\alpha,\beta,\gamma)>0$, away from
$(k-i)\alpha=i\beta$, $1\le i<k$. Set

\[
 m=\left\lfloor\frac{k\alpha}{\alpha+\beta}\right\rfloor,
 \qquad a=k-m,\quad b=m+1,\quad a+b=k+1.                 \tag{H4}
\]

**Proposition 1.** The initial ideal of $J_k$ is

\[
 \operatorname{in}_w J_k=(XY,X^aZ,Y^bZ).                \tag{H5}
\]

Its standard monomials are all the pure $X,Y$ powers, counting $1$
once, together with

\[
 Z^c[T^j],\qquad c\ge1,\quad -(b-1)\le j\le a-1.       \tag{H6}
\]

There are $k$ entries in every positive $Z$ layer. In particular the
carrier has exactly $k$ positive initial ideals, independent of $\gamma$
and of the characteristic.

**Proof.** In $G_i$ the two possible largest terms are
$X^{k-i}Z$ and $Y^iZ$. Their coefficients are respectively $1$ and
$(-1)^k$, so they survive in every field. Every intermediate Laurent
power has strictly smaller cost than the larger endpoint. The remaining
term $Y^i$ is smaller than $Y^iZ$. Equations (H3)--(H4) therefore
supply the heads $X^aZ,Y^bZ$, in addition to $XY$.

It remains to prove independence of the proposed standard monomials;
head containment alone would be insufficient. A finite relation among
them in the Laurent ring has the form

\[
 P_0(T)+\sum_{c=1}^C (T-1)^{-kc}P_c(T)=0,
 \qquad T^{b-1}P_c\in K[T]_{\le k-1}\ (c>0),             \tag{H7}
\]

with $P_0$ Laurent. At $T=1$, multiplication by $(T-1)^{kC}$
forces $P_C$ to vanish to order at least $k$. The polynomial
$T^{b-1}P_C$ has degree at most $k-1$, hence is zero. Induction
removes every positive layer; independence of distinct Laurent powers
then removes $P_0$. Polynomial division by the displayed relations
terminates in any refining positive monomial order and leaves only
these independent monomials. This proves (H5), with no generic-rank or
characteristic assumption. All $k$ intervals in (H4) are nonempty and
their head ideals are distinct. $\square$

## Generic evaluation truncates rank, not merely degree

For a finite list of rational functions regular at the generic $t_j$,

\[
 \operatorname{rank}_K(h(t_j))
 =\min\{n,\dim_F\langle h\rangle\}.                     \tag{H8}
\]

Indeed independent functions admit a nonzero evaluation determinant:
choose successive rows over an algebraic closure avoiding the zero
locus of a nonzero rational function. Its infinitude suffices even
when $F$ is finite. Thus that determinant is nonzero in independent
indeterminates. Dependence gives the reverse inequality. This is the
classical generic-evaluation argument; it is not asserted for fixed
numerical evaluation points.

For any strict positive weight, reduction by a carrier Gröbner basis
replaces a monomial by smaller standard monomials. Hence every carrier
cost prefix has rank equal to the number of its carrier-standard
monomials. By (H8), the standard set of (H1) is precisely the $n$
cheapest members of (H6) and of the two pure axis chains.
This gives a complete selection rule, including all $k,n$.

**Corollary 2.** The entire positive initial fan of $I_{n,k}$ is
independent of the characteristic. It is the same selection rule in
every field. The coefficients of reduced bases and ranks of arbitrary
nonprefix subsets can still depend on the characteristic. Equal fan
counts alone would not prove this stronger statement.

The equality on strict weights off the carrier walls proves equality
of the whole finite fans: every open initial cone contains such weights.
The established finite positive-weight representation accounts for
all global monomial orders. Vanishing intermediate binomial coefficients
therefore gives no generic-sample arithmetic amplification here.

## A large-order stabilization, with its full scope

Let $C_n$ be the ideal of $n$ generic points on the cylinder $XY=1$,
with independent parameters $X=s_j,Z=v_j,Y=s_j^{-1}$. Its residue
matroid on canonical monomials $[T^j]Z^c$ is the rank-$n$ truncation
of distinct Laurent monomials. Its standard set consists of the $n$
cheapest such canonical monomials.

**Proposition 3.** If $n\le2k+1$, $I_{n,k}$ and $C_n$ have exactly
the same initial standard set at every strict positive weight, and
therefore the same whole positive fan.

**Proof.** Use (H4), and let
$\tau=\min(a\alpha+\gamma,b\beta+\gamma)$, the cost of the first
carrier head beyond $XY$. The sector inequalities give

\[
 (a-1)\alpha<b\beta,\qquad (b-1)\beta<a\alpha.           \tag{H9}
\]

All $k$ positive-layer monomials with $c=1$ in (H6) cost strictly
less than $\tau$. The same $k$ Laurent powers at $c=0$ do too.
Include in addition $X^a$ if the first head is $X^aZ$, or $Y^b$
otherwise. That axis term costs strictly less than its head.
These $2k+1$ distinct monomials are all carrier-standard and independent.

Every canonical monomial below $\tau$ is carrier-standard, since
all nonstandard canonical monomials are multiples of one of the two
minimal heads. Thus the cylinder's $n$ cheapest canonical monomials
lie strictly below $\tau$ when $n\le2k+1$, and they are independent
in the curve. They are exactly the set selected by (H8). The dense
strict-weight argument extends the identity across both full fans.
$\square$

The threshold is a sufficient uniform bound; no optimal threshold is
claimed. The conclusion does **not** apply to every degree-$n$ polynomial
$f$ avoiding $0,1$, nor to every numerical point sample. Special
samples can lose additional evaluation minors. It also does not say
that the whole residue matroids coincide; only the coherent choices do.

The generic-sample restriction is essential even for reduced quotients.
For $k=n\ge2$ over $\mathbb Q$, take instead
$f(T)=(T-1)^k-2$. It avoids $0,1$, and $Z=1/2$ in its quotient.
Every initial standard set is a two-axis allocation of $n$ Laurent
powers in $X,Y$, so it has exactly $n$ initials. These allocations are
independent: multiplication by a suitable power of the unit $T$ turns
them into $1,T,\ldots,T^{n-1}$; every allocation is coherent.
The generic cylinder has these $n$ initials when $Z$ is expensive and
also the pure $Z$ chain when $Z$ is cheap, hence at least $n+1$.
Thus the stabilization would be false if stated for arbitrary $f$.

## A bound and a remaining construction window

The [two-plane threshold charge](extremal_initial_ideals.md#a-reduction-for-carriers-with-two-coordinate-planes)
applies here with $XY-1$ in place of $(X-1)Z-1$. The two ingredients
used in that proof remain valid: every carrier-standard monomial lies
in the union of the $XZ$ and $YZ$ exponent planes, and each of
$K[X,Z]$, $K[Y,Z]$ generates the finite quotient because $X,Y$ are
inverse units. Its prefix hypothesis is exactly (H8). Applying the
proved rank-truncation charge and rational-shadow estimate to the
$k-1$ carrier walls with primitive normals

\[
 q_i=\frac{(k-i,-i,0)}{\gcd(k,i)},\qquad1\le i<k,
\]

therefore gives

\[
 g(I_{n,k})=O\bigl(n\log(n+1)(n+k^{4/3})\bigr).          \tag{H10}
\]

For the cylinder carrier alone the wall set is empty, so the same
argument gives $g(C_n)=O(n^2\log(n+1))$. By Proposition 3, this
also bounds every $k\ge(n-1)/2$. For smaller $k$, (H10) in particular
excludes polynomial amplification when $k=O(n^{3/4})$. It gives the
uniform all-parameter upper bound $O(n^{7/3}\log(n+1))$.

Bound (H10) alone left the range $n^{3/4}\ll k<(n-1)/2$.
The following joint wall count closes that range and extends to arbitrary
reciprocal denominators. The earlier bound is retained as an intermediate
comparison; it is superseded for this entire generic family by (H11).

## Rational functions regular at zero and infinity have the same generic fan

Let $P,Q\in F[T]$ be coprime, $\deg P=k\ge1$, $P(0)\ne0$ and
$\deg Q\le k$. Replace $(X-1)^kZ-1$ in (H1) by $P(X)Z-Q(X)$
and denote the generic finite-point ideal by $I_{n,P,Q}$. Generic
parameters avoid both $0$ and the zeros of $P$. These rational functions
$Z=Q(T)/P(T)$ have no pole at $0$ or infinity;
no restriction is imposed on the zeros or the other poles. Rescaling
$Z$ permits normalizing $P$ to be monic without changing initial ideals.

**Proposition 4.** Every conclusion about the generic fan in
Propositions 1--3 and Corollary 2 depends only on $k$, not on $P,Q$'s
coefficients, roots or multiplicities. In particular $I_{n,P,Q}$ and
$I_{n,k}$ have exactly the same positive initial fan.

For heads, multiply $ZP(X)-Q(X)$ by $Y^i$ and reduce $XY=1$.
The two endpoint terms from $P$ survive in every field. Every Laurent
power from $Q$ lies between these endpoints, because $\deg Q\le k$;
its cost is smaller than the larger endpoint's cost after the positive
$Z$ cost is added. Thus the heads remain exactly (H5).

For independence, clear the largest power of $P$ in a relation of the
form (H7), replacing $(T-1)^{-kc}$ by $Q(T)^c/P(T)^c$.
Its top Laurent numerator multiplied by $Q^C$ must be divisible by $P$.
Coprimality implies divisibility of that numerator by $P$. A Bézout
identity also shows that adjoining $Q/P$ adjoins $1/P$, so the stated
carrier ideal is the actual rational-function kernel. Multiplication
by $T^{b-1}$ gives a polynomial of degree at most $k-1$; since
$P(0)\ne0$, Laurent divisibility implies ordinary divisibility.
The numerator is zero. Induct on the positive layer, as before.
Generic evaluation then gives the same selection rule. This includes
many distinct poles, repeated poles and inseparable denominators in
positive characteristic. The hypotheses on the degree, coprimality
and constant term are essential to this proof; rational functions
with poles at $0$ or infinity are outside this theorem.

## A degree-independent quadratic bound closes the remaining range

**Theorem 5.** With an absolute constant, over every field and uniformly
in $n$, $k$ and the coefficients of $P,Q$,

\[
 \boxed{g(I_{n,P,Q})=O(n^2\log(n+1)).}                    \tag{H11}
\]

This excludes polynomial amplification throughout the whole generic
rational-function family regular at $0$ and infinity, including the former intermediate range.
The key additional input is a moment budget on the part of a carrier
wall that can affect the finite basis. The general shadow estimate in
(H10) did not retain this restriction.

### The sensitive part of a carrier wall has small $Z$ moment

Fix a carrier wall

\[
 A\alpha=B\beta=t_0,\qquad A,B\ge1,\quad A+B=k.          \tag{H12}
\]

Normalize $\alpha,\beta$ and vary $\gamma>0$. The two carrier strips
exchange their endpoints $X^A$ and $Y^B$ in every positive $Z$ layer;
these endpoints have the same cost. Thus their ordered cost multisets,
and their $n$-th cost $T$, agree on the wall. Call the wall **sensitive**
at this weight if

\[
 T\ge t_0+\gamma.                                     \tag{H13}
\]

If (H13) fails, the first possible nonstandard canonical monomial is
above the finite selection threshold. The finite fan agrees with the
generic cylinder fan on a neighborhood of this weight. This is a
strict-gap assertion, not merely a comparison at the single weight.

Suppose (H13) holds. In the common part of the two strips take

\[
 H=\{1,X,\ldots,X^{\lfloor A/2\rfloor},
             Y,\ldots,Y^{\lfloor B/2\rfloor}\}.
\]

Every member has cost at most $t_0/2$, and
$|H|=1+\lfloor A/2\rfloor+\lfloor B/2\rfloor\ge k/2$.
Let $D$ be any coherent standard set incident to the weight and put
$C=\max_{u\in D}u_Z$. Its monomials have cost at most $T$, hence
$C\gamma\le T$. All the distinct monomials

\[
 hZ^c,\qquad h\in H,\quad0\le c\le\lfloor C/2\rfloor
\]

are in the relevant carrier standard strip. They have cost at most
$t_0/2+C\gamma/2<T$, since $t_0<T$ by (H13). They are independent
and must be selected before the threshold, even when other monomials
tie at $T$. Therefore

\[
 n\ge\frac{k}{2}(\lfloor C/2\rfloor+1)
      \ge\frac{k(C+1)}4,\qquad
 \mu_Z(D)=\sum_{u\in D}u_Z<\frac{4n^2}{k}.              \tag{H14}
\]

An incident basis comes from one of the two carrier strips by a refining
weight. The bound applies to both sides and to every incident coherent
choice, not just one arbitrarily chosen basis.

### Count slope changes once across the whole wall

On the fixed wall, the minimum basis cost is the concave piecewise affine
function

\[
 \phi(\gamma)=\min_D\{
       \alpha\mu_X(D)+\beta\mu_Y(D)+\gamma\mu_Z(D)\}.
\]

Its successive slopes are nonnegative integers $\mu_Z$, in decreasing
order. An interior vertex of the finite fan on this wall exposes a
two-dimensional profile face. Its $Z$ moment cannot be constant:
otherwise that face lies in the intersection of the two independent
planes $\mu_Z=\mathrm{constant}$ and
$\alpha\mu_X+\beta\mu_Y=\mathrm{constant}$, and has dimension at
most one. Thus every such vertex has a positive integer slope drop.

For sensitive vertices both slopes are bounded by (H14). Their drops
occupy disjoint subintervals of $[0,4n^2/k]$, because the slope is
monotone along the **entire** wall. The sensitive set need not be an
interval; no reset of the budget is allowed between its components.
There are at most $4n^2/k$ sensitive vertices per carrier wall. The
$k-1$ walls together contribute at most $4n^2$ vertices.

Every nonsensitive carrier-wall vertex is already a vertex of the
single generic cylinder fan, by the neighborhood identity. The earlier
two-plane threshold charge with the wall-free carrier $(XY-1)$ bounds
all such cylinder vertices by $O(n^2\log(n+1))$. The same charge bounds
the off-carrier vertices of $I_{n,P,Q}$ by $O(n^2\log(n+1))$, without a
factor for the number of sectors. Its boundary contributes
$O(n\log(n+1))$ vertices. Euler's formula for the positive weight
simplex subdivision now proves (H11).

The supporting-face dictionary, concavity of minimum costs, integer slope
counting and the earlier planar charge are established tools. The changed
structural step is (H14), applied only where the carrier can influence
selection. This is a whole-family exclusion, not a sharp count for the
family or an unrestricted bound for $G_3$.

## Exact integer-cost rank and random access without quotient matrices

The same formula gives a concrete computational consequence. Suppose
$w=(\alpha,\beta,\gamma)$ consists of positive integers, and use increasing
$(w\cdot u,u_X,u_Y,u_Z)$ to break ties. In (H4), use the floor even at
a carrier wall; this selects the $Y$ endpoint, consistently with this
lexicographic refinement. For an integer $t\ge0$, the rank of all
monomials of cost at most $t$ in the carrier is exactly

\[
 \begin{split}
 R_k(w,t)={}&1+\lfloor t/\alpha\rfloor+\lfloor t/\beta\rfloor
              +\lfloor t/\gamma\rfloor\\
 &+\sum_{i=1}^{a-1}\max(0,\lfloor(t-i\alpha)/\gamma\rfloor)
  +\sum_{j=1}^{b-1}\max(0,\lfloor(t-j\beta)/\gamma\rfloor).
 \end{split}                                                     \tag{H15}
\]

For the generic length-$n$ quotient its prefix rank is
$\min(n,R_k(w,t))$. This follows by counting the axes and two strips,
with no evaluation or matrix inversion. It is independent of $P,Q$'s
coefficients.

Let $L=\min(a-1,\max(0,\lfloor(t-\gamma)/\alpha\rfloor))$.
The first sum in (H15) is the classical integer floor sum
$\sum_{i=0}^{L-1}\lfloor(t-\alpha-\alpha i)/\gamma\rfloor$;
the second is symmetric. The Euclidean algorithm evaluates each in
$O(\log W)$ integer operations, $W=\max(2,\alpha,\beta,\gamma)$,
after signed quotient extraction. The iteration count is independent of
the strip width and of $t$. The repository's existing exact floor-sum
routine is reused; this arithmetic primitive is not a new algorithm.

**Proposition 6.** For $1\le r\le n$, the $r$-th standard monomial
of any generic $I_{n,P,Q}$ can be found using
$O(\log^2(rW+1))$ integer arithmetic operations, after reading $k$,
without enumerating earlier monomials or constructing quotient matrices.
This is an integer-cost bound, not a field-operation or floating-log audit.

First replace $k$ by $\min(k,r)$: if this changes it, Proposition 3
identifies both first-$r$ problems with the cylinder. Refine weights by

\[
 \widetilde w=r^3w+(r^2,r,1).                            \tag{H16}
\]

Any selected prefix of size $r$ is a downset, so its coordinates are
less than $r$. On those exponents (H16) preserves integer primary
cost differences and implements the prescribed lexicographic tie order
by base-$r$ digits. All refined costs there are distinct. The desired
cost is at most $(r-1)\min_i\widetilde w_i$, using a pure axis chain.
Binary search with (H15) finds the first $t$ with rank at least $r$.
Every exponent of cost at most that bound is less than $r$, so this
cost level contains exactly one carrier-standard monomial.

Recover it from the three pure-axis candidates and the two bounded
Diophantine equations
$\widetilde\alpha i+\widetilde\gamma c=t$ and
$\widetilde\beta j+\widetilde\gamma c=t$, with positive mixed
exponents and the strip bounds. Extended Euclid gives the admissible
residue class, and endpoint rounding finds its bounded solution.
Uniqueness on the relevant box removes enumeration. All intermediate
arithmetic, after clamping $k$, uses $O(\log(rW+1))$-bit integers;
the displayed operation bound is not a unit-cost bit-complexity claim.

This is an implemented research oracle for **generic rational-function
quotients regular at $0$ and infinity, with exact integer costs**. It does
not compute their coefficients,
apply to special samples, replace the public smooth-number CLI, or extend
its certification range for logarithmic prime weights.

## Primary comparison and verification

[Onn--Sturmfels, *Cutting Corners* (1999), Corollary 5.4 and the paragraph
before Proposition 5.6](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Onn/Onn4.pdf)
explicitly formulate the point-configuration maximum and its historical
volume bound; special configurations can have more initials than generic
ambient points. [Babson--Onn--Thomas, Proposition 2.2, Theorem 2.4,
Corollary 2.5 and its proof](https://arxiv.org/pdf/math/0207135) extend
the framework to all length-$n$ quotients and give $G_3(n)=O(n^3)$.
These primary statements were rechecked for scope. The repository's
stronger unrestricted shadow bound has a complete internal argument,
but its strongest-current-priority comparison remains unresolved.
Failed searches and an unavailable original paper do not establish novelty.

The Laurent reductions, evaluation determinant and state-polytope
framework are classical. The old McMullen comparison remains recorded
for the previous continuous-shadow argument; it is not used as an
unproved input in the carrier, counting or random-access proofs. No priority
claim is made for the specialized carrier, stabilization, sensitive-wall
or random-access combination. The classical
Euclidean floor-sum primitive can also be inspected in the maintained
[AtCoder Library implementation](https://github.com/atcoder/ac-library/blob/master/atcoder/internal_math.hpp),
function `floor_sum_unsigned`. Only that function was compared; its input
limits are not imported into the Python integer implementation.

[Jeannerod--Neiger--Schost--Villard, Section 1.1, Problem 1,
Definitions 1.1--1.2 and Theorem 1.3](https://arxiv.org/html/1602.00651v1#S1.SS1)
starts with vector data and a multiplication matrix in standard Jordan
representation, and computes a shifted Popov interpolation basis in
$\widetilde O(m^{\omega-1}\sigma)$ field operations. The definitions
and theorem statement were inspected; the full algorithm proof was not
audited in this comparison. Our oracle returns exponents under a generic
rank rule, and neither computes those coefficient relations nor improves
that general-data bound. For special samples an approximant reduction
must preserve the cost order and count **all** relevant shifts; a
single-shift algorithm does not by itself bound this fan.

Reproduce from the repository root:

```bash
python3 -B scripts/run_pole_order_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker compares ordinary polynomial division and all S-pairs with
independently cleared Laurent numerators, checks coefficient ranks in
several characteristics, and enumerates complete small generic finite
fans. A separate cylinder oracle verifies the stabilization. General numerators and denominators,
whole-wall slope budgets, coefficient independence and exact random access
are independently checked as well. These are
supporting controls; the proofs above cover unbounded $n,k$ and all fields.
No kernel, certification range, previous manuscript or paused ML result
is changed. The original extremal goal remains active.
