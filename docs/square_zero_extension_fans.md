# Square-zero extensions: the minors that tangent weights omit

The remaining bounded-coordinate wall problem already includes arbitrary
quotients with $x^2=0$. Their extension data give an exact algebraic
description, but counting planar tangent weights alone cannot bound their
spatial wall directions. An all-field family below has a coherent wall
outside the entire tangent-weight set of its planar ideal. Its fundamental
circuit has unbounded size although its deformation uses only two weights.

These are a reduction and a proof-route obstruction. The extension
description, tangent-space framework and determinant criterion are
classical. No improved unrestricted $G_3$ bound, collective cubic budget,
polynomial lower separation or novelty claim follows here.

## The complete extension data

Put $B=K[y,z]$ and $R=B[x]/(x^2)$. Ideals $I\subset R$ with finite
quotient correspond exactly to triples
$$
 J\subseteq K_0\subset B,\qquad
 \phi\in\operatorname{Hom}_B(J,B/K_0),                    \tag{E1}
$$
where $J,K_0$ have finite colength and
$$
 I=\{f-xg:f\in J,\quad g\bmod K_0=\phi(f)\}.
                                                               \tag{E2}
$$
The zero quotient $B/K_0$ is allowed. This is a classification of the
given coordinate presentation, not of algebras up to isomorphism.

Indeed, let $J$ be the image of $I$ modulo $x$ and let
$K_0=\{g:xg\in I\}$. Multiplying a lift $f-xg\in I$ by $x$
shows $J\subseteq K_0$. Two lifts have the same $g\bmod K_0$,
and multiplication by $B$ shows that this defines $\phi$.
Conversely (E2) is stable under multiplication by $B$ and by $x$,
using $J\subseteq K_0$. These constructions are inverse.

There is an exact sequence of $B$-modules
$$
 0\longrightarrow B/K_0\xrightarrow{g\mapsto xg}R/I
   \longrightarrow B/J\longrightarrow0.                 \tag{E3}
$$
Thus $n=m+s$, where $m=\dim B/J$, $s=\dim B/K_0\le m$, and
the multiplication operator $M_x$ has rank $s$. The algebra is free
over $K[x]/(x^2)$ precisely when $J=K_0$: nilpotent Jordan blocks
have sizes one or two, and freeness is equivalent to $m=s$.
This prevents treating arbitrary square-zero quotients as flat ones.

## The exact minor model in a monomial flat fiber

Suppose now $J=K_0=M$ is monomial, with standard diagram $D$ of size
$m$. Let $D_2$ be the standard diagram of $M^2$. Every standard set
of the extension lies in
$$
 (\{0\}\times D_2)\ \sqcup\ (\{1\}\times D),             \tag{E4}
$$
because $x^2,xM,M^2\subset I$. The columns $(0,d)$ for $d\in D$
are independent coloops: all other nonzero columns have residues in
$x(B/M)$. In particular these plain $D$ monomials are compulsory
for every order, not just for one cone.

After contracting those columns, the residue matrix is exactly
$$
 [\,\operatorname{Id}_{D}\mid\Phi\,],\qquad
 \Phi_{v,u}=[v]\phi(u),\quad v\in D,\quad u\in D_2\setminus D.
                                                               \tag{E5}
$$
The identity column indexed by $v$ represents $xv$; a column indexed
by $u$ represents the plain monomial $u$. Choose $U\subseteq D_2\setminus D$
and $T\subseteq D$ with $|U|+|T|=m$, and put $S=D\setminus T$.
The resulting monomial set is a basis exactly when
$$
                  \det\Phi[S,U]\ne0.                    \tag{E6}
$$
This follows by eliminating the identity columns. Its profile is
$$
 \mu(D\sqcup xD)+\sum_{u\in U}(0,u)-\sum_{v\in S}(1,v).
                                                               \tag{E7}
$$
The positive lower normal fan of the convex hull of these profiles
is the entire initial-ideal fan. This is the established residue-matroid
dictionary applied to the finite complete ground (E4). Merely counting
the nonzero minors does not count the positively exposed profiles.

The same argument covers monomial nested endpoints $J\subseteq K_0$:
replace $D_2$ by the diagram of $JK_0$, keep the compulsory diagram
of $J$, and use the diagram of $K_0$ to index the identity rows.
The residual rank is then $s$, not $m$. This follows from
$xK_0,JK_0\subset I$ and (E3).

Every matrix entry is linear in the extension parameters. The minors
in (E6) have unbounded order. Their cancellations and simultaneous
exposure are the coupling information absent from a list of tangent
weights. Formula (E6) makes that missing information explicit; it does
not supply the required collective bound.

## Two tangent weights and an arbitrarily long coherent circuit

For integers $a,t\ge2$, over every field, consider
$$
 I_{a,t}=(x^2,\ y^a,\ z^t-xy^{a-1}(1+z)).              \tag{E8}
$$
Here $M=(y^a,z^t)$ and
$\phi(y^a)=0$, $\phi(z^t)=y^{a-1}(1+z)$.
The syzygy is respected, so (E1)--(E3) give length $n=2at$.
Set $u=xy^{a-1}$ and, for positive weights, put
$$
        \lambda=(w_x+(a-1)w_y)/w_z.
$$
Let $F$ consist of all plain $y^iz^j$ with $i<a,j<t$, and all
$xy^iz^j$ with $i<a-1,j<t$. There are exactly three initial
standard sets:
$$
\begin{array}{c|c}
\text{range}&\text{standard set}\\\hline
0<\lambda<t-1&F\sqcup\{uz^j:0\le j<t\}\\
t-1<\lambda<2t-1&F\sqcup\{u,z^t,z^{t+1},\ldots,z^{2t-2}\}\\
2t-1<\lambda&F\sqcup\{z^t,z^{t+1},\ldots,z^{2t-1}\}.
\end{array}                                                  \tag{E9}
$$

Here is a length-and-leading-term proof of the whole fan. The first
range has the pure-power heads $x^2,y^a,z^t$ in (E8), whose
complement has size $2at$. Always $xz^t,yz^t,z^{2t}\in I$.
Since $z^{2t}=0$ in the quotient, $1+z$ is a unit there, and
$$
       u=\sum_{j=0}^{t-1}(-1)^jz^{t+j}.                  \tag{E10}
$$
In the middle range use the heads
$x^2,y^a,xz^t,yz^t,uz,z^{2t-1}$. The head $uz$ comes from
(E8), and the last head from (E10). In the last range replace
$uz,z^{2t-1}$ by $u,z^{2t}$. These head ideals have exactly the
displayed complements of size $2at$. Since they are contained in
the corresponding initial ideals and the quotient length agrees,
they are equal. All inequalities are strict inside each range,
and the two boundary planes are disjoint in the positive weight
cone. Refinement at a boundary selects a neighboring standard set.
Thus there are no additional initial ideals.
Equivalently, the established finite positive-weight realization of
global monomial orders lets the open-cell classification count every
monomial initial ideal.

The two actual wall directions are
$$
 (1,a-1,-(t-1)),\qquad (1,a-1,-(2t-1)).                  \tag{E11}
$$
For the rectangular diagram of $M$, every planar tangent weight has
$z$ entry at least $-t$. The second wall is therefore outside the
entire tangent-weight set, not only the two weights present in $\phi$.
In the convention tail exponent minus head exponent, the lifted
weights present in $\phi$ are
$q_0=(1,a-1,-t)$ and $q_1=(1,a-1,1-t)$; the extra direction is
$2q_0-q_1$.

More explicitly, the only variable block has upper columns
$A_j=uz^j=e_j$ and lower columns
$B_j=z^{t+j}=e_j+e_{j+1}$, with $e_t=0$.
It has the fundamental circuit
$$
                 A_0=\sum_{j=0}^{t-1}(-1)^j B_j.        \tag{E12}
$$
Every coefficient is a unit in every field. The $t$ columns $B_j$
are a unit triangular basis, so this circuit has exactly $t+1$
members. At the second wall, the common $B_0,\ldots,B_{t-2}$
are cheaper than the tied $A_0,B_{t-1}$; $F$ supplies the fixed
coloops. This is actual positive exposure of that long circuit.
It proves unbounded propagation length, with only two tangent weights
at the starting planar ideal. It is not a large-fan construction:
the entire fan has three cells.

## A classical box gate for nilpotent lifts

The long circuits above do not make rectangular extensions candidates for
a polynomial separation. A classical moment-volume argument excludes a
larger class, without any restriction on the deformation coefficients.

Let \(A=K[x,y,z]/I\) have length \(n>0\), suppose \(x^h=0\) in
\(A\), and write \(m=\dim_K A/xA\). Suppose monic univariate
polynomials \(P(y),Q(z)\), of degrees \(r,s\ge1\), annihilate
\(A/xA\). No flatness or monomial-kernel assumption is required. Then
\[
 x^h,\quad P(y)^h,\quad Q(z)^h\in I.                 \tag{E13}
\]
Indeed \(P(y),Q(z)\in xA\), so their \(h\)-th powers vanish.
Consequently every global initial standard set lies in the common box
with sides \(h,hr,hs\). Coordinate chains in a standard set \(D\)
give \(\mu_i(D)\le n(d_i-1)/2\) for these three sides \(d_i\).
The convex hull of the actual standard-set moments therefore has
volume \(O(n^3h^3rs)\), after adding the integral unit cube if needed.
Its vertices include all the distinct initial standard sets. The
classical three-dimensional lattice-vertex estimate yields
\[
              g(I)=O(n^{3/2}h^{3/2}\sqrt{rs}).          \tag{E14}
\]
The constant is absolute, over every field. This is the anisotropic
form of [Babson--Onn--Thomas, proof of Corollary 2.5](https://arxiv.org/pdf/math/0207135),
already recorded as (MB1)--(MB2) in the research note; it is not a new
lattice-polytope theorem. Actual undirected fan adjacency is planar,
and each primitive wall cost is at most \(n\), so also
\(V(I)=O(n^{5/2}h^{3/2}\sqrt{rs})\).

For a convenient feasibility parameter put \(\eta=rs/m\ge1\).
The multiplication filtration gives
\[
 m\le n\le hm,\qquad
 g(I)=O(h^{3/2}n^2\sqrt\eta).                         \tag{E15}
\]
For the upper length bound, multiplication by \(x^j\) maps
\(A/xA\) onto \(x^jA/x^{j+1}A\). If \(A\) is free over
\(K[x]/(x^h)\), then \(n=hm\), and the sharper substitution in
(E14) is
\[
                 g(I)=O(hn^2\sqrt\eta).                \tag{E16}
\]
Thus fixed \(h\) and \(rs=O(m)\) exclude a polynomial gain over
the generic quadratic scale, including nonflat extensions. At fixed
\(h\), a family with \(g(I)\ge c n^{2+\epsilon}\) requires
\(\eta=\Omega_{c,h}(n^{2\epsilon})\). Large \(\eta\) is only
a necessary condition, not a construction or a sufficient test.

In particular, take arbitrary polynomials \(f,g\in K[x,y,z]\):
\[
 C=K[x,y,z]/(x^h,\ P(y)-xf,\ Q(z)-xg).               \tag{E17}
\]
This has \(n=hrs\) and is free over \(K[x]/(x^h)\), with basis
\(y^az^b\), \(a<r,b<s\). Here is an elementary flatness proof.
Use monic univariate division on \(P,Q\); every error term gains an
\(x\) factor, so division terminates modulo \(x^h\), proving spanning.
For independence, work modulo \(x^h\) and suppose a polynomial supported
on this proposed basis is \(a(P-xf)+b(Q-xg)\). If both multipliers
are divisible by \(x^k\), the \(x^k\)-coefficient has zero rectangular
remainder and is itself rectangular, so it vanishes and
\(a_kP+b_kQ=0\).
Since the coprime polynomials \(P,Q\) have only their Koszul syzygy,
\(a_k=cQ,b_k=-cP\). Subtract \(x^kc(Q-xg)\) from \(a\) and
add \(x^kc(P-xf)\) to \(b\), preserving the represented polynomial
while raising the coefficient order. Repeating through \(h\) orders
proves independence. This is the elementary regular-sequence lifting
argument; the general flatness input is classical, for example
[Stacks Project, Lemma 10.128.6](https://stacks.math.columbia.edu/tag/00R3).

For (E17), (E16) gives \(g(C)=O(hn^2)\) and
\(V(C)=O(hn^3)\). Fixed nilpotence order therefore closes every such
two-relation deformation, including arbitrary multipliers, composite
orders, nonlocal planar fibers and exceptional characteristics. Adding
further relations that preserve the rectangular special fiber can lose
flatness, but (E15) still gives \(O_h(n^2)\). The earlier unit grid
is the special case \(h=2,P=y^r,Q=z^s,g=0\); its multiplier need
not be a unit for this exclusion.

These are classical feasibility consequences, not a new unrestricted
counting theorem or evidence of publication significance. They do not
apply with small \(\eta\) to arbitrary planar endpoints, and do not
identify a compulsory vector-space contraction with an algebra quotient.
The checker compares two terminating division orders with independent
ideal-span elimination in the larger carrier (E13), then traverses the
complete small fans. Its examples support the algebra and moment
inequalities; the uniform bounds rest on the proof above.

## Literature comparison and the next proof step

[Chaput--Evain, Section 2.1 and Proposition 2.1](https://aif.centre-mersenne.org/item/10.5802/aif.2955.pdf),
printed pages 1210--1211, give the established significant-cleft and
arm/leg descriptions of the tangent weights. These passages and the
proposition's proof were inspected. For a box $e$ of $D$, the weights
are $(\operatorname{arm}(e),-\operatorname{leg}(e)-1)$ and
$(-\operatorname{arm}(e)-1,\operatorname{leg}(e))$, with multiplicity.
This describes first-order deformations of the planar ideal, not all
walls of the resulting three-variable quotient.

[Hering--Maclagan, Remark 4.4, Theorem 4.5 and Definition 4.2](https://arxiv.org/pdf/1110.1861),
PDF pages 14--15, describe significant arrows and the recursive
polynomials for the planar cells. The definition and worked Example
4.3, as well as those statements, were inspected; the delegated Evain
proof of Theorem 4.5 was not re-audited. Those recursions already
contain products of arrow parameters. Nothing here contradicts their
theorems or establishes a new deformation principle.

The saved complete characteristic-two thick-hook fan of length 102
has 153 distinct wall directions with nonzero $x$ entry. Its planar
ideal has 102 tangent parameters and 84 distinct tangent weights.
The checker independently recomputes these counts and the Hom-space
nullity, retaining this as finite evidence against a parameter-count
cap, without inferring an asymptotic growth law.

Rectangular nilpotent lifts are excluded by (E13)--(E17), so larger
fan searches there are unnecessary. The useful next target is a
collective bound on positively exposed
minors of the compatible matrices (E5), or a realizable family in
which those minors yield a polynomial coherent separation. A bound
on tangent weights alone cannot replace that step. Arbitrary
square-zero quotients also allow $J\ne K_0$ and nonmonomial endpoints;
the flat monomial model must not silently be substituted for them.
The unrestricted $G_3$ objective and all publication gates remain open.

Run from the repository root:

```bash
python3 -B scripts/run_square_zero_extension_checks.py
./build.sh
python3 tests/test_smoke.py
```
