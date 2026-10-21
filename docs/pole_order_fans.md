# Growing pole order: exact carrier fans and generic stabilization

The three-pole example has quadratically many initial ideals. A possible
amplification is to replace its third simple pole by a pole of order $k$.
This note gives an exact description of the resulting **generic-sample**
problem. Its whole positive fan is independent of the characteristic,
and for $n\le2k+1$ it stabilizes to the generic cylinder fan. Thus neither
vanishing binomial coefficients nor increasing $k$ indefinitely amplifies
this generic family. Arbitrary special samples are not covered.

This is a reduction for an actual infinite family, not a resolution of
$G_3(n)$, a polynomial lower bound, or a novelty certificate. The strongest
current unrestricted comparison remains open.

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

A polynomial gain, if one exists in this family, must consequently
occur along $k/n^{3/4}\to\infty$ with $2k+1<n$. There is no lower
bound proving such a gain, and no numerical growth fit is made.
The task is now a count of genuine coherent cuts in that remaining
range, not a search for larger binomial coefficients or more pole orders.
The selected rule supplies a concrete route to that count; if a uniform
quadratic bound holds, this whole generic family must be excluded too.

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
unproved input in (H3)--(H10). No priority claim is made for the specialized
carrier and stabilization combination.

Reproduce from the repository root:

```bash
python3 -B scripts/run_pole_order_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker compares ordinary polynomial division and all S-pairs with
independently cleared Laurent numerators, checks coefficient ranks in
several characteristics, and enumerates complete small generic finite
fans. A separate cylinder oracle verifies the stabilization. These are
supporting controls; the proofs above cover unbounded $n,k$ and all fields.
No kernel, certification range, previous manuscript or paused ML result
is changed. The original extremal goal remains active.
