# Coordinate projection tags and their exact limit

Let $I\subset K[x,y,z]$ have finite positive colength $n$, over any
field. A coordinate tag can refine its fan at length at most $2n$.
The refinement is exactly computable: it adds one planar elimination
fan. All three planar elimination fans are invariant under tagging,
so repeated tags cannot provide an indefinitely iterated amplification.
After the three coordinate directions have been used, the fan stabilizes.

This settles the proposed iterative mechanism, without settling the
unrestricted extremal question. A single tag can still have a large
gain. The prime-digit computations below give exact finite gains of
quadratic size; they prove no asymptotic upper or lower growth rate.

## One tag over every field

Write $R=K[x,y,z]$, $J=I\cap K[x,y]$, and
$m=\dim_K K[x,y]/J\le n$. For $k\ge1$ put
$$
T_{z,k}(I)=JR+z^k I. \tag{PT1}
$$
For every global monomial order, with its induced order on $K[x,y]$,
$$
\begin{split}
\dim_K R/T_{z,k}(I)&=n+km,\\
\operatorname{in}T_{z,k}(I)&=(\operatorname{in}J)R
                         +z^k\operatorname{in}I. \tag{PT2}
\end{split}
$$
In particular the standard set is
$$
D_{T_{z,k}(I)}=
 \bigsqcup_{j=0}^{k-1} z^jD_J\ \sqcup\ z^kD_I. \tag{PT3}
$$
Here the two displayed portions are disjoint by their $z$ exponents.

**Proof.** The variable $z$ is regular on
$R/JR=(K[x,y]/J)[z]$. Since $JR\subset I$,
$(T_{z,k}(I):z^k)=I$. Multiplication by $z^k$ gives
$$
0\longrightarrow R/I\longrightarrow R/T_{z,k}(I)
 \longrightarrow R/(JR,z^k)\longrightarrow0,
$$
proving the length formula. Leading terms give inclusion of the
right side of (PT2) in the left. Since
$(\operatorname{in}J)R\subset\operatorname{in}I$, the right side
has precisely the standard set (PT3), of size $n+km$. Equal finite
colengths prove equality. $\square$

Both $D_I$ and $D_J$ can be recovered from (PT3). The positive fan is
therefore exactly the common refinement
$$
\mathcal F(T_{z,k}(I))=\mathcal F(I)\vee\mathcal F(J). \tag{PT4}
$$
The planar fan is pulled back to the three weight coordinates. Only
pairs of open cones that intersect occur; its cell count is not a
product of the two counts. The fan is independent of $k$.

## Exact gain from a pencil of lines

Let $Q_J$ be the distinct primitive positive planar wall directions,
embedded as $q=(a,-b,0)$ with $a,b>0$. For an open parent cone $C$
let $C$ count toward $T_I(q)$ exactly when its interior meets
$q^\perp$. In the normalized weight simplex, these wall lines all
meet at the boundary vertex $e_z$, and have no interior intersection
with one another. Thus
$$
g(T_{z,k}(I))=g(I)+\sum_{q\in Q_J}T_I(q). \tag{PT5}
$$
Indeed each line that crosses the interior of a convex parent cell
adds one cell. Different pencil lines add their increments separately.
A line coincident with a parent wall adds nothing there. If $P_C$ is
the closed parent polygon, the exact test is
$$
\min_{v\in P_C}q\cdot v<0<\max_{v\in P_C}q\cdot v. \tag{PT6}
$$

The existing [rational-shadow bound](rational_shadow_bound.md) bounds
the number of parent cells cut by a rational plane by
$O(nL_n(1+\|q\|_1^{1/3}))$, where $L_s=\log(s+1)$.
The [planar wall budget](wall_cost_products.md) gives
$$
|Q_J|=O(mL_m),\qquad \sum_{q\in Q_J}\|q\|_1\le m(m-1).
$$
For $m\ge2$, Holder's inequality therefore yields
$$
g(T_{z,k}(I))-g(I)
 =O\left(nm^{4/3}L_nL_m^{2/3}\right). \tag{PT7}
$$
For $m=1$ the gain is zero. This is an upper bound on the added
cells, not a new unrestricted upper bound for $G_3$ or a separating
lower construction. It leaves open a one-step gain above quadratic.

## The projection fans never change

Take $k=1$ and write $T_z=T_{z,1}$. If
$J_{xz}=I\cap K[x,z]$ and $J_x=I\cap K[x]$, then
$$
\begin{split}
T_z(I)\cap K[x,y]&=J,\\
T_z(I)\cap K[x,z]&=J_xK[x,z]+zJ_{xz}. \tag{PT8}
\end{split}
$$
For the second equality, decompose a polynomial as $f_0(x)+zf_1(x,z)$.
Its membership in $T_z(I)$ implies $f_0\in J_x$ by setting $z=0$.
Subtract $f_0$ and use $(T_z(I):z)=I$ to obtain $f_1\in J_{xz}$.
The reverse inclusion is immediate. The $yz$ formula is identical.

The second ideal in (PT8) is the two-variable version of (PT1) with
a univariate base. A univariate ideal has one initial ideal for
every order. Formula (PT3) therefore recovers exactly the original
fan of $J_{xz}$. Consequently **all three planar projection fans
are invariant** under each coordinate tag.

It follows that any finite sequence of $x,y,z$ tags has fan
$$
\mathcal F(I)\vee\bigvee_{i\text{ used}}\mathcal F(J_{[3]\setminus\{i\}}).
\tag{PT9}
$$
Repeating a direction adds no cones. Once each direction has been
used, further tags change length but not the fan. This argument is
uniform in the field, coefficients, nilpotents, and original fan.

## Completing all coordinate projections

For $S\subseteq\{x,y,z\}$ let $R_S=K[S]$,
$J_S=I\cap R_S$, and $x_S=\prod_{i\in S}i$. Put
$J_\varnothing=0$ and $\dim R_\varnothing/J_\varnothing=1$.
Tagging once in each direction, in any order, gives
$$
I^\sharp=\sum_{\varnothing\ne S\subseteq\{x,y,z\}}x_SJ_SR. \tag{PT10}
$$
For every monomial order its standard set and length are
$$
\begin{split}
D_{I^\sharp}&=\{1\}\ \sqcup\!
  \bigsqcup_{\varnothing\ne S\subseteq\{x,y,z\}}x_SD_{J_S},\\
\dim R/I^\sharp&=1+\sum_{\varnothing\ne S\subseteq\{x,y,z\}}
                     \dim R_S/J_S\le8n. \tag{PT11}
\end{split}
$$
The pieces have different exact coordinate supports, so they are
disjoint. To verify their dimension independently, let
$A_S=R_S/J_S\subset R/I$. On $V=\bigoplus_S A_S$ define the action
of coordinate $i$ on a summand $A_S$ as multiplication by its
residue if $i\in S$, and as the natural inclusion into $A_{S\cup\{i\}}$
if $i\notin S$. These operators commute. The vector $1$ in
$A_\varnothing$ is cyclic, and $x_Sf$, $f\in J_S$, kills it.
Thus (PT10) annihilates $V$. The proposed standard monomials span
the quotient, and their images form a basis of $V$, proving (PT11).
Equivalently the formulas follow by successively applying (PT2)
and (PT8). The same argument in fixed dimension $d$ costs at most
$2^dn$ length.

For a radical ideal of points with every coordinate nonzero,
$I^\sharp$ is the ideal of all their coordinate-face projections,
including the origin. Different exact supports give disjoint face
sets, and the count is (PT11). Nonzero coordinates are necessary:
$(x,y,z)^\sharp=(x^2,y^2,z^2)$ has length eight, whereas the
set-theoretic face projections of the origin have length one.

## Prime digit tests

In characteristic two take
$$
I_{n,k}=((X-1)^n,XY-1,Z-Y^k).
$$
With $M$ the least power of two at least $n$, a monomial $X^aY^bZ^c$
has residue column $(\binom{a-b-kc\bmod M}{j}\bmod2)_{0\le j<n}$.
The $XY$ elimination ideal has its $n$ hooks, with distinct wall
lines $(j,-(n-j),0)^\perp$, $1\le j<n$. Formula (PT5) therefore
gives the exact count after a $Z$ tag, at length $2n$.

| $n$ | $k$ | $g(I_{n,k})$ | Tagged count | Gain |
|---|---|---|---|---|
| 5 | 3 | 13 | 27 | 14 |
| 10 | 3 | 39 | 102 | 63 |
| 10 | 5 | 43 | 112 | 69 |
| 21 | 5 | 129 | 408 | 279 |
| 21 | 13 | 173 | 485 | 312 |
| 42 | 13 | 501 | 1710 | 1209 |
| 42 | 21 | 497 | 1822 | 1325 |

These are fixed finite cases. The ratios gain/$n^2$ range from
$14/25$ to $1325/1764$. Their bounded range establishes no uniform
estimate or polynomial separation. No larger search is justified
by these data alone. The remaining question for this operation is
a structural lower bound for the pencil-crossing sum in (PT5).

## Prior art and verification

The algebraic operation $a+fI$ and the inclusion plus equal Hilbert
function method are classical.
[Migliore--Nagel, Applications of Liaison](https://academicweb.nd.edu/~jmiglior/MN14-webpg.pdf),
Section 8, Definition 8.1, Theorem 8.2, Remark 8.3 and Lemma 8.4,
give the homogeneous basic double link, its Hilbert function, its
module isomorphism, and an initial-ideal criterion. Their stated
Cohen--Macaulay and codimension hypotheses apply to the homogeneous
version here with $a=JR$ and $f=z^k$; $R/a$ is a polynomial ring
over a finite algebra, hence is Cohen--Macaulay of dimension one.
The affine proof of (PT2) avoids importing homogeneous assumptions.
The theorem proofs cited there are delegated to earlier papers;
those earlier proofs were not re-audited for this comparison.
Neither this operation nor its length argument is claimed new.
Exact priority of the combined quantitative statements is unsettled.

Run from the repository root:

```bash
python3 -B scripts/run_projection_tag_checks.py
./build.sh
python3 tests/test_smoke.py
```

The [saved report](../results/research/projection_tag_checks.json)
includes complete parent cone witnesses, exact hit vectors, commands
and source hashes. The checker independently uses direct point
evaluation, nonreduced CRT columns, nilpotent residue columns,
and iterated tag modules. It compares whole cones, not sampled
weights alone. All six coordinate tag orders and repeated-tag
stabilization are checked on specified finite examples. The uniform
statements rest on the proofs above. No unrestricted $G_3$ bound,
polynomial lower separation, current novelty claim, or completed
goal follows from this construction.
