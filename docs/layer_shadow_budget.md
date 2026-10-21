# Layered shadows and a collective sign-class wall budget

The finite-height constraint on one coordinate improves every relevant
rational shadow, even when the other two entries of its kernel are large.
For a three-variable quotient of length $n$ containing a monic degree-$h$
polynomial in $x$, every actual wall direction has a shadow with
$O(hn\log(n+1))$ vertices. This gives a cubic polylogarithmic collective
cost bound for all walls whose $y,z$ entries have the same sign or one
is zero, at each fixed $h$. Walls with opposite nonzero $y,z$ signs
remain uncontrolled by that budget.

The proof extends the saved rational-shadow argument to translated
layers. Its token paths, lattice-direction estimate, state-polytope
dictionary and perimeter charge are established inputs. It supplies no
new unrestricted $G_3$ bound, polynomial lower separation, sharpness
claim for the new sign-class budget, or novelty certificate.

## A translated-layer shadow theorem

Let a rank-$r\ge1$ matroid have labels $a_e=(x_e,y_e,z_e)$ in
$\mathbb N_0^3$ with
$$
 0\le x_e<h,\qquad (x_e+1)(y_e+1)(z_e+1)\le H.
 \tag{LS1}
$$
Replace $h$ by $\min(h,H)$, put $L=1+\lfloor\log_2 H\rfloor$,
and let $P$ be its basis-profile polytope. Write $s_q(P)$ for
the number of vertices in its projection with primitive kernel $q$.

**Theorem 1.** If $q_x\ne0$, put $s=|q_x|$. Uniformly in the other
two coordinates of $q$,
$$
 s_q(P)=O\!\left(L\left[rh+(r^2Hs)^{1/3}h^{2/3}\right]\right).
 \tag{LS2}
$$
If $q_x=0$, then
$$
 s_q(P)\le2\bigl(r(h-1)+1\bigr). \tag{LS3}
$$
The bounds allow repeated labels, arbitrary matroids and all objective
quadrants. They concern exposed profiles, not every tied basis.

### Covering projected layers

Orient $q$ with $q_x=s>0$ and use the integer projection
$$
 A a=(sy-q_yx,\ sz-q_zx).
$$
The lattice $A\mathbb Z^3$ has covolume $s$: up to signs, its three
maximal minors are $s^2,sq_z,sq_y$, with gcd $s$ because $q$ is primitive.
The $x=j$ layer is a translate of $s\mathbb Z^2$ and satisfies
$(y+1)(z+1)\le H_j=\lfloor H/(j+1)\rfloor$.

For $0\le k< L_j=1+\lfloor\log_2 H_j\rfloor$, put
$B_y=2^{k+1}$ and $B_z=\lfloor H_j/2^k\rfloor$.
The rectangles
$$
 (-jq_y,-jq_z)+
 (-\tfrac12,sB_y-\tfrac12)\times
 (-\tfrac12,sB_z-\tfrac12) \tag{LS4}
$$
cover the projected layer, have area at most $2s^2H_j$, and have
no projected ground point on their boundary. There are at most $hL$
rectangles and $4hL$ boundary lines.

Assign every projected point to the first rectangle containing it,
without retaining its original layer as part of the label. This is
essential when translated layers overlap. All points assigned to a
rectangle belong to the global lattice of covolume $s$, not necessarily
to the layer sublattice of covolume $s^2$.

### Token and direction budgets

In a quadrant sweep of planar matroid optimization, the classical
adjacent greedy exchanges give $r$ convex token paths. Along every
linear functional each path has at most two monotone phases.
Consequently each boundary line in (LS4) is crossed at most twice
per token. At most $8rhL$ exchanges change rectangle labels.

Within a rectangle of widths $W_1,W_2$, the absolute coordinate
variation of each token is at most $2W_i$, even if the token leaves
and returns. Sum the internal exchange vectors at each parameter.
They have distinct ray directions, with coordinate budgets at most
$2rW_1,2rW_2$. The saved lattice-direction estimate bounds their
number by
$$
 6+(1728r^2W_1W_2/s)^{1/3}
 \le6+(3456r^2sH_j)^{1/3}. \tag{LS5}
$$
Here the covolume is that of the global projected lattice. Charging
every genuine basis-profile change to either kind of exchange gives
the following explicit bound, after four quadrant sweeps:
$$
 4\left[1+8rhL+
 \sum_{j=0}^{h-1}L_j
       \left(6+\left\lceil(3456r^2sH_j)^{1/3}\right\rceil\right)\right].
 \tag{LS6}
$$
Finally
$\sum_{j=0}^{h-1}H_j^{1/3}
 \le H^{1/3}\sum_{j=1}^h j^{-1/3}
 \le\tfrac32H^{1/3}h^{2/3}$.
This proves (LS2). The geometry and token argument are exactly those
of the [rational-shadow proof](rational_shadow_bound.md); the smaller
normalized rectangle area is the additional quantitative input.

For (LS3), project to $(x,-q_zy+q_yz)$. Every basis profile has
integer first coordinate between zero and $r(h-1)$. At each such
level a convex polygon has at most two vertices. Summing levels
proves the assertion, including segments and points.

## Consequences for actual capped quotients

Let $I\subset K[x,y,z]$ have positive colength $n$ and contain a
monic univariate degree-$h$ polynomial in $x$, over any field.
Use the actual residue matroid on
$$
 E=\{a\in\mathbb N_0^3:x<h,\ \prod_i(a_i+1)\le n\}.
 \tag{LS7}
$$
It has rank $n$ and retains every initial standard set. Every actual
positive profile edge is parallel to an exchange of two ground labels.
Thus a primitive wall direction with $q_x\ne0$ has
$1\le|q_x|\le h-1$. Substituting $r=H=n$ in (LS2) gives
$$
 s_q(P)=O(hn\log(n+1)) \tag{LS8}
$$
for every primitive direction with $0<|q_x|<h$, whether or not it is
an actual wall direction; (LS3) covers $q_x=0$ too. No cyclicity,
generic rank, reducedness or separability assumption is used.

Both the number $p_I(q)$ of parallel undirected wall pieces and
the number $T_I(q)$ of open parent cones strictly cut by $q^\perp$
are $O(s_q(P))$. For the latter, a uniquely exposed vertex in the
projected polygon can expose at most two endpoints of its fiber.
For the former a parallel exposed edge lies in such a vertex fiber.
This is the existing fiber argument, including collapsed edges.
In particular
$$
 p_I(q),T_I(q)=O(hn\log(n+1)) \tag{LS9}
$$
for every actual wall direction. This uniform bound removes the
dependence on the other two entries of $q$; it is not the stronger
unproved residual-rank assertion $p_I(q)=O(n-\nu(q)+1)$.

## A cubic polylogarithmic budget for one full-support sign class

Let $V_{\rm same}(I)$ sum the primitive costs $\nu(q)$ over actual
undirected walls with $q_yq_z\ge0$, retaining multiplicities.
For a mixed direction in this class, orient $q_x>0$, so
$q=(c,-a,-b)$ with $a,b\ge0$, not both zero, and $1\le c<h$.
The exact threshold is
$$
 \nu(q)=c+(a+1)(b+1)-1\le n. \tag{LS10}
$$
For each $c$, there are at most
$D(n)=\sum_{u=1}^n\lfloor n/u\rfloor\le n(1+\log n)$
possible pairs; imposing primitivity can only reduce the count.
Their total primitive cost is at most $nD(n)$.
Combining with (LS9) proves
$$
 V_{\rm same}(I)=O(h^2n^3\log^2(n+1)). \tag{LS11}
$$
There is also a cubic coordinate-level budget for $q_x=0$.
Every positive standard profile has integer $x$-moment between zero
and $n(h-1)$, and total degree moment at most $\binom n2$.
At one fixed $x$-moment, horizontal edges remain on the positive
lower chain of the hull of profiles at that level. Their primitive
cost $|q_y|+|q_z|$ is no larger than their coordinate displacement.
Perimeter telescoping gives at most $n(n-1)$ per level, hence
$$
 V_{x=0}(I)\le\bigl(n(h-1)+1\bigr)n(n-1)=O(hn^3).
 \tag{LS12}
$$
This is the coordinate-level extension of the established
[planar weighted perimeter budget](wall_cost_products.md).

At fixed $h$, (LS11)--(LS12) cover every intermediate cost in these
two specified families. The missing walls have $q_x\ne0$ and
opposite nonzero $y,z$ entries, for example $(c,a,-b)$.
Their threshold is $(c+1)(a+1)+b-1$, with a linear rather than
hyperbolic restriction on $a,b$. Applying (LS9) separately to all
of them does not give a cubic budget.

The loss already occurs at $h=2$. All directions $(1,a,-b)$ with
$a,b\ge1$ are primitive, and their cost is $2a+b+1$.
Consequently
$$
 \#\{(a,b):2a+b+1\le n\}=\Theta(n^2),\qquad
 \sum_{2a+b+1\le n}(2a+b+1)(n-2a-b)=\Theta(n^4).
 \tag{LS14}
$$
For the lower estimate, take
$\lceil n/16\rceil\le a\le\lfloor n/8\rfloor$ and
$\lceil n/4\rceil\le b\le\lfloor n/3\rfloor$.
For sufficiently large $n$ these give $\Omega(n^2)$ directions
with both cost and residual rank proportional to $n$; the upper
estimate follows from $O(n^2)$ pairs and summands $O(n^2)$.
Thus even an individual $O_h(r\log(r+1))$ multiplicity estimate,
where $r=n-\nu(q)+1$, would only give a quartic polylogarithmic
envelope on this class. These are admissible directions, not walls
proved to coexist in one quotient. The missing step must constrain
coexistence or jointly charge their multiplicities.

## Coordinate tagging cannot supply a new power gain here

The [coordinate-tag formula](projection_tag_fans.md) adds a pencil
of planar elimination walls, with gain $\sum_q T_I(q)$.
If the tag is in the bounded $x$ direction, its added normals
have $q_x=0$. The planar wall count is $O(n\log(n+1))$ and
(LS3) gives gain $O(hn^2\log(n+1))$.

For a tag in $y$ or $z$, its planar elimination ideal still contains
the same univariate relation in $x$. Its profiles lie in a rectangle
of side lengths $O(hn),O(n^2)$. The classical lattice-polygon vertex
bound gives $O(h^{1/3}n)$ planar walls. Their nonzero $x$ entries
are at most $h-1$, so (LS2) gives gain
$$
 O(h^{4/3}n^2\log(n+1)). \tag{LS13}
$$
At fixed $h$, a tag therefore adds only a quadratic polylogarithmic
number of cells to any parent in this class. The parent fan itself
is not bounded quadratically by this argument. No claim is made
about growing $h$, including the cyclic unit-power families.

## Prior-art boundary and independent checks

[Onn--Rozenblit, arXiv:1208.5639v1](https://arxiv.org/pdf/1208.5639v1),
Theorem 2.3 and its proof, use a coordinate-height bound and a
zonotope refinement; the projected-edge step is explicit on PDF
pages 5--6. Those passages were inspected. Theorem 4.1 and its proof
on PDF pages 10--11 give the matroid specialization through unit
coordinate exchanges and greedy optimization; those passages were
also inspected. The new estimate instead uses translated
layer rectangles and their normalized areas. This comparison is
not a strongest-current or novelty certificate.

Small shadows alone do not justify a quadratic total fan count.
[Lagarias--Luo--Padrol](https://ems.press/journals/lem/articles/16381)
gives a logarithmic-over-logarithmic worst-case largest-shadow
guarantee for general three-dimensional polytopes. Its published
abstract was inspected here; the proof is not an input to (LS2).
No additional relation between all shadows and the number of
vertices is assumed for these quotient polytopes.

Run from the repository root:

```bash
python3 -B scripts/run_layer_shadow_checks.py
./build.sh
python3 tests/test_smoke.py
```

The report records complete quadrant sweeps, independent full
basis-profile hulls, rectangle-label and covolume checks, token
boundary crossings and per-rectangle direction budgets, and actual
finite quotient fans. Large translations and overlapping layer
rectangles test the essential scope restrictions. The uniform
statements rest on the proof, not a growth fit. The unrestricted
endpoint, further proof requirements and publication gates remain
unchanged.
