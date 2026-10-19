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

All previous manuscripts, smooth-number kernels and the paused ML work
are preserved. The sharp $G_3(n)$ question and the active goal remain open.
