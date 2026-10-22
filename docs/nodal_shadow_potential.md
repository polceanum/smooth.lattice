# Global coordinate shadows and a bounded coupling potential on yz = 0

Let $A=K[x,y,z]/I$ have length $n>0$, over any field, with $yz\in I$.
Write $A_Y=A/zA$, $A_Z=A/yA$, their lengths as $n_Y,n_Z$, and their
numbers of initial ideals as $p_Y,p_Z$, with $g(0)=1$.
The [coordinate-border reduction](projection_border_fans.md) bounds each
pair of quotient cones separately. Here are controls that hold across
the changing encodings:

1. Every generic coordinate price line crosses at most
   $1+4n\lceil\log_2 n\rceil$ parent initial cones.
2. After subtracting the two quotient costs, both transverse slopes lie
   in $[0,n_Y]\times[0,n_Z]$. The resulting continuous function is
   nonnegative, coordinatewise increasing and supermodular.
3. Mixed parent walls have total lattice-weighted incidence at most
   $2n_Y(p_Y-1)+2n_Z(p_Z-1)$ with the moving quotient grid. Thus this
   part of the interaction is $O(n^2\log(n+1))$ collectively.

The last bound counts walls meeting the grid, including their positive
interior endpoints. It does **not** count all junctions inside the grid
regions. The unrestricted $G_3(n)$ bound is unchanged. An explicit
all-field, growing-length family below also shows that the residual
cost need not be concave. It cannot generally be replaced by the
minimum support function of one fixed residual polytope.

These are supporting structural controls and a failed-route certificate,
not a solution of the extremal problem or a publication assessment.
Exact quantitative priority remains unassessed.

## A matroid lemma with unrestricted constant-price elements

Let $M$ be any finite rank-$r$ matroid. Give each element a cost
$a_e+t c_e$, $t>0$. Elements with $c_e>0$ have integer labels
$a_e\ge0$, $c_e\ge1$ satisfying

$$
                    (a_e+1)(c_e+1)\le H.                \tag{1}
$$

The other elements have $c_e=0$ and **arbitrary real** $a_e$;
their number and magnitudes are unrestricted. Fix a tie order. With
$L=\lceil\log_2 H\rceil$, the number of canonical greedy bases on
open parameter intervals, and hence of distinct minimum profiles, is at most

$$
                 1+L\{r+3(r^2H)^{1/3}\}.               \tag{2}
$$

For $H=1$ put $L=0$: all prices are constant. Representation over a
field is not assumed. Degenerate objectives do not contribute every
tied optimal basis to this count.

**Proof.** At a cost crossing, reverse tied groups by adjacent swaps.
The greedy adjacent-swap property either leaves the basis unchanged or
exchanges the earlier selected element for the later unselected one.
Transfer one of $r$ tokens at each successful swap. Increasing $t$
strictly decreases $c$ and strictly increases $a$ along each token path.
An element with $c=0$ absorbs its token permanently. In particular,
the set of selected constant-price elements is nested by inclusion.

Give $c=0$ label zero, and $2^j\le c<2^{j+1}$ label $j+1$,
$0\le j<L$. The token-label sum decreases, so at most $rL$
successful swaps cross blocks, including all absorptions. A token
visits each positive block at most once. Within block $j$ its total
$c$ decrease is at most $2^j-1$, and its total $a$ increase is at
most $H/2^j$. Sum these two budgets over the tokens:

$$
                 C_j\le r2^j,\quad A_j\le rH/2^j,
                 \quad C_jA_j\le r^2H.                 \tag{3}
$$

For each distinct within-block crossing choose one successful exchange
there. Its positive integer direction $(\Delta c,\Delta a)$ has
slope $t=\Delta a/\Delta c$. The elementary distinct-direction
lemma in the [planar dyadic proof](planar_dyadic_bound.md) bounds the
number of these slopes by $3(C_jA_j)^{1/3}$. Sum over $L$ blocks
and the cross-block events. Arbitrary constant-price intercepts enter
only at absorption; (3) never uses them. This proves (2), including
simultaneous crossings by adjacent resolution. $\square$

The token construction is classical: [Eppstein, Section 4, Theorem 4
and proof, printed pp. 8--9](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf)
already uses it for general parametric matroids. Neither token passing
nor the unchanged-price monotonicity is claimed as a new principle.
The displayed lemma extends this repository's hyperbolic-label charging
to an unrestricted absorbing axis; the cited theorem is not itself the
quantitative assertion (2). No exhaustive priority claim follows.

## Apply the lemma to the whole nodal quotient

The standard monomials of any length-$n$ quotient lie in the universal
hyperbolic cross $\prod_i(a_i+1)\le n$. Since $yz=0$, only the two
coordinate planes are needed. Set $w_x=1$, fix any real $u=w_y>0$,
and vary $v=w_z>0$. Every $z$-positive surviving monomial is $x^a z^c$
with cost $a+vc$ and $(a+1)(c+1)\le n$. All other ground elements
have cost $a+ub$, constant in $v$. Apply (2) with $r=H=n$:

$$
               \#\text{canonical sweep bases}
                    \le1+4n\lceil\log_2 n\rceil.        \tag{4}
$$

There is no rational-height factor for $u$. A generic line meets this
many open full-dimensional cones at most. For a line lying on fan
walls, (4) counts canonical choices using one fixed monomial tie order;
it does not count every cone closure touching that line. The symmetric
bound holds when $v$ is fixed and $u$ varies. Small shadows alone do
not bound the full spatial fan by their squared size.

## One global, bounded, supermodular potential

For any of these quotients let $F_C(w)$ be the minimum total monomial
weight of a residue basis. Put $B=A/(yA+zA)$, $d=\dim B$; its
standard basis is $1,x,\ldots,x^{d-1}$, or empty if $d=0$. Define

$$
 H(w)=F_A(w)-F_{A_Y}(w)-F_{A_Z}(w)+\binom d2 w_x.        \tag{5}
$$

The function is continuous, homogeneous and piecewise linear. At any
generic weight let $D,B_Y,B_Z$ be the corresponding standard sets.
The border lift proved in the preceding note gives

$$
 B_Y\cup B_Z\subseteq D\subseteq
 \{1,x,\ldots,x^{a-1}\}\cup yB_Y\cup zB_Z,
 \qquad a=\dim K[x]1.                                  \tag{6}
$$

The entire extra $y$-label budget is $n_Y$, and the extra $z$-label
budget is $n_Z$, by telescoping the column heights of the planar
downsets. Consequently, throughout every full-dimensional patch,

$$
       0\le\partial_{w_y}H\le n_Y,\qquad
       0\le\partial_{w_z}H\le n_Z.                     \tag{7}
$$

The slopes are integers. These bounds are uniform across every change
of the two quotient bases. Continuity gives coordinatewise monotonicity
globally on $w_x=1$, including exceptional lines by limits.

As $v$ increases at fixed generic $u$, only the prices of $z$-positive
elements increase. No selected unchanged-price element can leave:
raise the changed prices one at a time and use the adjacent-swap
property. Hence the selected $y$-positive elements form a nested set,
and $\mu_y(D)$ is nondecreasing in $v$. The subtracted $F_{A_Y}$
has no $v$ dependence, while $F_{A_Z}$ has zero $u$ derivative.
Thus $\partial_uH$ is nondecreasing in $v$. Integrate in $u$ to get

$$
 H(u_2,v_2)+H(u_1,v_1)\ge H(u_1,v_2)+H(u_2,v_1),
 \quad u_1<u_2,\ v_1<v_2.                              \tag{8}
$$

Here $H(u,v)$ abbreviates $H(1,u,v)$. Limits remove all genericity
restrictions in (8). This elementary matroid substitutes argument
explains supermodularity; it is not an additional new optimization
theory claim.

For a primary-source comparison, [Gul--Pesendorfer--Zhang, February 2020
version](https://www.princeton.edu/~fgul/walras2019.pdf), footnote 1 on
printed p. 1, states unchanged-price demand retention as the substitutes
property. Appendix A, printed pp. 26--27, gives the basis-exchange
description and discusses weighted-matroid substitutes. Those inspected
passages establish that this monotonicity is prior art; they do not
state the hyperbolic bound (2), the algebraic residual slope budget (7),
or the grid-incidence count (12). This is a limited comparison of
definitions and displayed arguments, not a full audit of that paper or
of the earlier sources it cites.

## A filtered-rank interpretation, including the overlap

Put $W=yA\cap zA$, $h=\dim W$, and filter $A$ by the spans $V_t$
of monomials of cost at most $t$. Let $r_A,r_Y,r_Z,r_B$ be the
dimensions of their images in the four quotients, and put
$k_t=\dim(V_t\cap W)$. The two images of $V_t$ map onto the same
image in $B$. Their compatible pairs form a space of dimension
$r_Y+r_Z-r_B$. The image of $V_t$ in this space has dimension
$r_A-k_t$. Its codimension is therefore

$$
                  c_t=r_Y+r_Z-r_B-r_A+k_t\ge0.          \tag{9}
$$

The vector-space quotient identity is
$n-n_Y-n_Z+d=h$. Integrating the usual greedy rank-deficiency
formula $F_C=\int_0^\infty(\dim C-r_C(t))\,dt$ gives

$$
                 H(w)=\int_0^\infty(h-k_t+c_t)\,dt\ge0.\tag{10}
$$

These are elementary rank and integration identities. In particular,
when $W=0$, $H$ measures the integrated failure of the two filtered
quotient images to glue. No persistence-theory novelty is asserted.
The universal hyperbolic ground suffices at every threshold: positive
greedy selection has a downset of at most $n$ accepted monomials,
so every rank increase lies in that ground.

## A collective budget on moving-encoding walls

Call a parent wall mixed if its primitive profile-edge direction $q$
satisfies $q_yq_z<0$. Write the actual profile difference as $m q$,
$m\ge1$. Its trace in the $(u,v)$ quadrant has positive slope.
Along a generic fixed-$u$ line, every mixed crossing increases
$\mu_y(D)$ by $m|q_y|$. Every other crossing also has a nonnegative
$y$ jump, by unchanged-price retention.
The quotient $B_Y$ is fixed on that line, so (7) gives

$$
 \sum_{\text{mixed crossings}}m|q_y|\le n_Y.            \tag{11}
$$

Symmetrically a generic fixed-$v$ line has sum at most $n_Z$ with
weight $m|q_z|$. This counts actual lattice jumps, including their
multiplicities, rather than just primitive directions.

There are $p_Y-1$ vertical quotient walls and $p_Z-1$ horizontal
ones. On a vertical quotient wall, take two sufficiently nearby
generic parallel lines. Every mixed edge meeting the quotient wall
at a finite strictly positive point, even just an endpoint, meets
at least one of these lines. Thus its total weighted incidence is
at most $2n_Y$. The horizontal argument gives $2n_Z$. Summing yields

$$
 \boxed{\sum_{\text{quotient-grid incidences}}\!!\text{lattice weight}
       \le2n_Y(p_Y-1)+2n_Z(p_Z-1).}                    \tag{12}
$$

Count a wall separately on each grid line it meets. Endpoints on
the boundary of the positive projective quadrant, including infinity,
are excluded. The planar bound $p_Y,p_Z=O(n\log(n+1))$ makes
(12) an $O(n^2\log(n+1))$ bound. Internal vertices and walls confined
to one grid region remain uncontrolled by (12); they are the next
counting obligation, not a negligible remainder.

## Why a single global residual polytope is invalid

For every $N\ge3$ and every field, take

$$
 A_N=K[x,y,z]/(x^N,xy,y^2,z-x-y),\qquad n=N+1.           \tag{13}
$$

Its basis is $1,x,\ldots,x^{N-1},y$, and $yz=0$.
The ideals $yA_N$ and $zA_N$ have zero intersection. The two actual
coordinate quotients have lengths $2,N$ and two initial ideals each;
the common axis quotient has length one. Write $X,Y,Z$ for the
three prices and $C=\binom N2$. In degree one, any two of $x,y,z$
form a basis. For every higher surviving degree $k$, all nonzero
monomials are parallel to $x^k=z^k$ and have cheapest cost
$k\min(X,Z)$. Thus the parent has exactly four initial sets and

$$
 F_{A_N}=(C-1)\min(X,Z)+\min(X+Y,X+Z,Y+Z),
 \quad F_{A_Y}=\min(X,Y),\quad F_{A_Z}=C\min(X,Z).
$$

Subtract to obtain the uniform formula

$$
                      H_N(X,Y,Z)=\max(0,\min(Y,Z)-X).   \tag{14}
$$

At $(1,1/4,3)$ and $(1,5/4,3)$ this has values $0,1/4$;
at their midpoint $(1,3/4,3)$ it has value zero. This violates
concavity, even with zero overlap and bounded quotient fan sizes.
Adding a linear term cannot repair that violation. Also a quotient
wall can pass through the interior of a parent cone. The parent fan
need not refine either quotient fan.

The correct global object is the bounded supermodular potential
(5), with the nonnegative filtered defect (10). A local residual
matroid remains valid on each pair of quotient cones; assuming those
local models form one global fixed-rank matroid loses (14).

## Reproducible controls and remaining target

Run `python3 -B scripts/run_nodal_shadow_checks.py`. The exact report is
`results/research/nodal_shadow_checks.json`. It checks arbitrary-axis
matroids by adjacent swaps and exhaustive independent basis minima,
including verified nonrepresentable basis-exchange data, loops and
repeated labels. Actual quotient controls use independently commuting
coordinate matrices, complete parent and coordinate-quotient fans,
every nonempty overlay patch, direct greedy coordinate sweeps,
lattice-weighted grid incidences and independent dense filtered ranks.
The counterfamily is proved uniformly above; finite checks over
rationals and small prime fields are supporting controls.

The original extremal target still requires a collective bound on
internal junctions, or an actual algebraic family exceeding the generic
scale by a polynomial factor. Neither a coordinate-shadow estimate,
a small boundary-incidence estimate, nor the existence of the potential
settles that requirement. No new unrestricted count, lower separation,
novelty conclusion or completed paper is claimed here.

## Internal junctions need a second count: a verified failed route

The next test rules out bounding all internal junctions by (12). A
junction can have only vertical and horizontal incident walls, while
both coordinate quotients have one fixed monomial initial ideal.
The issue is not just that the mixed-wall bound has a loose constant.

First, the classical single-interval exchange test extends to two
separable prices. Suppose one group has costs $a_e+b_eu$ and the
other has costs $c_f+d_fv$, with $d_f>0$. Fix a first-group tie
at $u=u_0$, with common cost $c$, and a fixed prefix $C$ of the
other first-group elements. Set

$$
                 D(v)=\{f:c_f+d_fv<c\}.
$$

As $v$ increases these sets shrink in a fixed nested order. The tied
pair exchanges under adjacent greedy reversal exactly when it is
nonloop parallel in $M/(C\cup D(v))$. In the opposite direction,
contraction only decreases ranks: once the pair is nonloop parallel,
it stays parallel until both become loops. If just one element has
become a loop, it cannot return. Thus eligibility is one interval,
possibly empty. No equality of the $d_f$ is needed. This is the same
elementary rank argument used in the
[two-colour threshold proof](two_colour_threshold_braids.md), not a
new exchange theorem. It applies to one specified pair and prefix;
different pairs on one primitive wall direction are not identified.

The interval property does not prevent a vertical exchange segment
from crossing a horizontal exchange segment at different cost levels.
Here is an actual algebra showing this, over **every** field.

Take these eight maximal monomials and their prescribed residue columns:

| Monomial | Column |
| --- | --- |
| $x^4y^3$ | $e_1$ |
| $x^3y^4$ | $e_2$ |
| $y^6$ | $e_3$ |
| $x^9y$ | $e_4$ |
| $z^6$ | $e_1+e_2$ |
| $x^9z$ | $e_1$ |
| $x^3z^5$ | $e_4$ |
| $x^4z^4$ | $e_3+e_4$ |

They are an antichain under divisibility. Let $D$ be the union of
their divisor boxes and let $H$ be the monomial algebra with standard
set $D$. Its length is $36+40-10=66$. All eight tip monomials are
killed by $x,y,z$. Any linear subspace of their span is therefore an
algebra ideal. Quotient by the kernel of the displayed rank-four
column map. This gives a cyclic commutative algebra $A$ of length
$62$, with $yz=0$. All $58$ proper-divisor columns are compulsory
in every monomial basis; the remaining choice is exactly the
rank-four tip matroid.

Both sets of four tip columns, the $Y$ set and the $Z$ set, are
bases of the same four-dimensional space. Killing either coordinate
therefore kills that entire tip space and its own branch's
positive-coordinate core. Both actual coordinate quotients are
monomial algebras. Their lengths are $32,36$, and $p_Y=p_Z=1$.

Set $w=(1,u,v)$. On $99/100<u,v<101/100$, the minimum tip cost is

$$
 6v+\min(4+3u,3+4u)+6u+\min(3+5v,4+4v).                \tag{15}
$$

Indeed the first two columns and the first two $Z$ columns form a
rank-two direct summand. Its cheap column is $z^6$, and either
$Y$ column completes it; $x^9z$ is too expensive. The other summand
has cheap column $y^6$, and either of its $Z$ columns completes it;
$x^9y$ is too expensive. These assertions use determinants $\pm1$,
so they hold in every characteristic.

Equation (15) has four strict regions meeting at $(u,v)=(1,1)$.
The four incident walls are all axis walls. There is no mixed wall
in this neighbourhood. Since the quotient grid is empty, (12) is
zero, yet the parent has an internal four-cone junction. The complete
tip matroid has $25$ bases: each summand has the five pairs other
than its one parallel pair. Exact all-basis halfplane intersection
recovers $15$ positive parent cones and agrees with independent
full-algebra border traversal. This is a proof obstruction, not a
large-fan construction.

### An unbounded version of the obstruction

For every $s\ge1$, the same mechanism gives both $p_Y=p_Z=1$ and
at least $s^2$ pure-axis internal junctions in one open price box.
The following deliberately generous integer scales give an explicit
construction. Put $m=2s$ and

$$
 Q=1000(m+1)^6,\quad G=100(m+1)^4,\quad
 K=10(Q+G+m),\quad D=K(m+2).
$$

For block $1\le i\le m$, put $a_i=iK$, $b_i=D-iK$. Define two
pairs of planar labels:

$$
 T_i=((a_i,b_i),(a_i+Q+i,b_i-Q)),\qquad
 R_i=((a_i,b_i-G),(a_i+Q,b_i-Q+G)).                       \tag{16}
$$

For $i\le s$, use $T_i$ in the $Y$ plane and $R_i$ in the $Z$
plane. Give the two $Y$ tips columns $e_{2i-1},e_{2i}$ and the
two $Z$ tips columns $e_{2i-1}+e_{2i},e_{2i-1}$. For $i>s$,
use $R_i$ in $Y$ and $T_i$ in $Z$, with $Y$ columns
$e_{2i-1},e_{2i}$ and $Z$ columns
$e_{2i},e_{2i-1}+e_{2i}$. Different blocks are direct summands.

All these tips are positive in their branch coordinate. Within a
plane, increasing $x$ strictly decreases the other exponent: the
block spacing exceeds every local displacement. They are therefore
one antichain. Form the divisor algebra and the annihilated tip-space
quotient as above. Its common tip rank is $2m$, both actual coordinate
quotients are monomial, and the construction works over every field.

On the box $1<u,v<1+(m+1)/Q$, the two $T_i$ costs lie between
$D$ and $D+G$. The first $R_i$ cost is below $D$ and the second
is at least $D+G$. The inequalities follow from
$Q>2G$ and $D(m+1)/Q<G/2$, with $m<G/2$.
Thus every vertical block selects its cheap $Z$ tip and the cheaper
of its two $Y$ tips. Its only wall in the box is
$u=1+i/Q$. Every horizontal block selects its cheap $Y$ tip and
the cheaper $Z$ tip, with only wall $v=1+i/Q$.
The box consequently contains exactly $(s+1)^2$ parent regions and
$s^2$ pure-axis junctions, with no incident mixed edges there.

The report supplies compressed antichain, downset-length and direct
summand certificates for this family. It does not enumerate the huge
compulsory cores of the growing examples. The small length-$62$
example is separately checked as a full commuting algebra. These
inefficient exponents are only an obstruction to a proposed proof;
they give no polynomial separation above the generic scale.

The next obligation is therefore a joint count of axis-axis and
mixed junctions. A bound solely on mixed curvature, mixed-wall
boundary incidence, or single exchange eligibility intervals cannot
replace that count. The original extremal target remains unchanged.
In particular this does not refute a possible bound with an additional
$O(n^2)$ term; it identifies a contribution that such a proof must cover.
