# Bounded-coordinate wall shadows and a sharp residual logarithm

Let \(I\subset K[x,y,z]\) have colength \(n\), over any field, and
suppose a monic univariate polynomial of degree \(h\) in \(x\)
belongs to \(I\). Replace \(h\) by \(\min(h,n)\) if necessary.
Every initial standard set then has \(x\)-exponents less than \(h\).
For a primitive actual wall direction \(q\), retain the definitions
\[
 \nu(q)=\prod_i(q_i^++1)+\prod_i(q_i^-+1)-2,\qquad
 r=n-\nu(q)+1
\]
and count all undirected parallel pieces in \(p_I(q)\).

**Theorem 1.** If \(q_x\ne0\), then \(a=|q_x|\le h-1\), and
\[
 p_I(q)=O\!\left(h(1+a^{2/3})\,r(1+\log(r+1))\right).
                                                               \tag{B1}
\]
The constant is absolute. For fixed \(h\), this is
\(O_h(r\log(r+1))\), improving the general residual
\(O(r^{4/3}\log^2(r+1))\) estimate in this class.
No reducedness, cyclicity or generic-rank hypothesis is used.

Let \(V_x^{\rm opp}(I)\) sum \(\nu(q)\) over the actual walls whose
orientation has \(q_x>0,\ q_y,q_z\le0\).
Zero entries are allowed, but \(q\) has both signs.

**Theorem 2.**
\[
 V_x^{\rm opp}(I)
       =O\!\left(h^{8/3}n^3(1+\log(n+1))^2\right).
                                                               \tag{B2}
\]
Together with the [coordinate-slice theorem](coordinate_slice_wall_budget.md),
this controls all walls with \(q_x=0\), and all walls where the
bounded coordinate is alone on one side. In particular it covers
every two-coordinate direction.
At fixed \(h\), the remaining directions have the forms
\((a,b,-c)\) and \((a,-b,c)\), with \(1\le a<h\), \(b,c>0\).
Their costs are linear in \(b,c\) at fixed \(h\); they are not
controlled by (B2).

**Theorem 3.** The logarithm in (B1) cannot generally be removed.
There is a sequence of actual characteristic-zero quotients with
\(x^2=0\), and also radical rational two-plane realizations of the
same fans, for which the single fixed direction
\[
                    q=(1,-1,-1),\qquad \nu(q)=4             \tag{B3}
\]
has \(p_I(q)=\Omega(r\log(r+1))\).
More precisely, for the two-layer prefix families \(I_n\),
\[
 \max_{N/2\le n\le2N+1}p_{I_n}(1,-1,-1)
                         =\Theta(N\log(N+1)).              \tag{B4}
\]
The lower assertion is for the maximum in this interval, not every
individual length. This contrasts with the linear residual bound
when \(h=2\) and the bounded coordinate of \(q\) is zero.

These are restricted wall theorems. They do not prove a full spatial
cost bound, improve unrestricted \(G_3\), give a polynomial fan-count
separation, or establish novelty. The original goal remains open.

## A translated planar hyperbolic-cross lemma

Let a rank-\(r\) matroid have integer planar labels in a union of \(t\)
translates of
\[
 m V_H^2,\qquad V_H^2=\{(u,v)\ge0:(u+1)(v+1)\le H\},
 \qquad m,H\in\mathbb Z_{\ge1}.
\]
Anchors can be arbitrarily far apart and labels may repeat. Set
\(\ell=1+\lfloor\log_2 H\rfloor\).
In any open objective quadrant the number of exposed profiles is at
most
\[
 1+t\ell\left(4r+3(r^2m^2H)^{1/3}\right).                 \tag{B5}
\]
Four such expressions bound the full profile polygon.

Here are the details needed for the translation step; merely applying
a bound separately to matroid restrictions would not suffice.
Cover each unscaled cross by the \(\ell\) integer rectangles
\[
 2^j-1\le u\le2^{j+1}-2,\qquad
       0\le v\le\lfloor H/2^j\rfloor-1.
\]
Their side-length products are at most \(H\), and after scaling
at most \(m^2H\). Surround each translated rectangle by vertical
and horizontal lines a half-unit outside its integer endpoints.
There are at most \(4t\ell\) such lines, none containing a label.

Sweep a positive ratio in any one quadrant. The usual adjacent-swap
matroid greedy argument gives \(r\) token paths, with each planar
coordinate monotone in its prescribed direction. Each token crosses
each covering line at most once. Thus at most \(4rt\ell\) successful
exchanges cross a line.

The endpoints of any other exchange are in the same cell of this
line arrangement. A rectangle containing one contains both, since
otherwise one of its surrounding lines separates them. Assign the
exchange to one such rectangle. For each token and rectangle, the
sum of the assigned absolute coordinate changes is bounded by the
two rectangle side lengths. This follows from global coordinate
monotonicity even if exchanges are assigned to other rectangles
in between.

For a rectangle of sides \(W,Z\), sum exchanges at each parameter.
They have one positive slope at that parameter. The elementary
distinct-integer-direction estimate gives at most
\[
            3(r^2WZ)^{1/3}\le3(r^2m^2H)^{1/3}
\]
different such parameters. Sum over all rectangles and add the
line-crossing charge. This proves (B5), including segments, point
polygons, overlapping covers and repeated labels. It is a translated
cover version of the [existing planar dyadic argument](planar_dyadic_bound.md).
It does not split the common rank into ranks of the covered pieces.

## One common residual matroid, then one planar projection

Orient \(q\) so that \(q_x=a>0\). By the
[compulsory-core theorem](wall_coexistence.md), all endpoint standard
sets of all these walls contain the same literal independent set
\[
 C_q=(\operatorname{Div}(q^+)\cup\operatorname{Div}(q^-))
                                      \setminus\{q^+,q^-\}.
\]
Its rank is \(\nu(q)-1\). Contract it once. The residual rank is \(r\),
and its labels lie in at most six translates of \(V_r^3\).
Restrict that ground to \(0\le\alpha_x<h\).
Every counted endpoint still supplies a residual basis, so the
restricted contraction retains rank \(r\).

Project by the integer map
\[
 L_q(\alpha)=(a\alpha_y-q_y\alpha_x,\
                         a\alpha_z-q_z\alpha_x).           \tag{B6}
\]
Its real kernel is \(\mathbb Rq\). For a translate
\(\beta+V_r^3\), and a fixed \(\alpha_x=k<h\), a projected label is
in
\[
 (a\beta_y-q_yk,\ a\beta_z-q_zk)+aV_r^2.                   \tag{B7}
\]
Indeed the two remaining nonnegative increments have product
budget at most \(r\). Thus there are at most \(6h\) translated
scaled planar crosses. Their anchors need not be small.

Every original parallel wall remains the same exposed edge in
this single residual polytope, translated by the fixed core
profile. The support-function agreement is local at each wall;
no agreement at unrelated weights is asserted.
Each such edge projects to a different vertex under (B6):
the fiber over one exposed shadow vertex is a single segment.
Apply (B5) with \(H=r,m=a,t\le6h\), in four quadrants, to prove
(B1). The bound \(a<h\) follows from two different selected
monomials in one wall-grading fiber, both with \(x\)-exponents
less than \(h\); their difference is a nonzero multiple of \(q\).

## Summing an entire sign family

For the directions counted by \(V_x^{\rm opp}\), write
\(q=(a,-b,-c)\), with \(1\le a<h\), \(b,c\ge0\). Their threshold is
\[
             \nu(q)=a+(b+1)(c+1)-1.
\]
For each fixed \(a\), the number with \(\nu(q)\le n\) is at most
\(O(n(1+\log(n+1)))\), by the ordinary two-factor divisor count.
Dropping primitivity only enlarges this number. Each actual piece
has cost at most \(n\), and (B1) is at most
\[
 O\!\left(h(1+a^{2/3})n(1+\log(n+1))\right).
\]
Sum over \(a<h\), using
\(\sum_{a<h}(1+a^{2/3})=O(h^{5/3})\). This proves (B2).
Unlike a bound for distinct directions, it retains every parallel
piece. It includes the full intermediate-cost range of this sign
family.

## A fixed direction already needs the residual logarithm

Use the actual two-layer prefix fan on
\(\{0,1\}\times\mathbb N_0^2\).
For every \(n\), it has a local characteristic-zero realization:
in \(\mathbb Q[u]/u^n\), put \(s=\lceil n/2\rceil\),
\(C=n(1+n+n^2)\), and
\[
 x=u^s e^{(Cs+1)u},\qquad
 y=u e^{(C+n)u},\qquad z=u e^{(C+n^2)u}.                  \tag{B8}
\]
Then \(x^2=0\), \(y\) generates the algebra, and every size-\(n\)
slab downset is a basis. Its lower layer has at least \(s\)
monomials; assigning the lower layer its initial row indices and
the upper layer the remaining ones verifies the same monotone
Hermite conditions as in the [balanced slab proof](nested_jet_matroids.md).
This also covers odd \(n\).
The established generic two-plane construction gives a radical
realization over \(\mathbb Q\) at every \(n\).

On the positive plane \(w_x=w_y+w_z\), the projected labels are
\[
 (i,j)\quad\hbox{from }y^iz^j,\qquad
 (i+1,j+1)\quad\hbox{from }xy^iz^j.
\]
Thus an axis label has one representative and an interior label
has two. Take any ordinary planar free prefix \(E\) of size \(m\),
with a unique next label \(v=(i,j)\), \(i,j\ge1\).
Let \(A(E)\) count the labels of \(E\) on the union of the axes,
including the origin. Fully selecting \(E\) selects \(2m-A(E)\)
slab monomials. At
\[
                         n(E)=2m-A(E)+1                  \tag{B9}
\]
the next label is only half selected. A small perturbation of
\(w_x\) swaps precisely
\[
                    y^iz^j,\qquad xy^{i-1}z^{j-1}.
\]
It gives an actual \(q=(1,-1,-1)\) wall of \(I_{n(E)}\).
Different pairs \((E,v)\) give different pieces: their two endpoint
sets recover both the fully selected labels \(E\) and the last
partial label \(v\).

For each \(m\), every ordinary planar corner cut has an open
interval with a unique next label, so there are at least \(C_2(m)\)
such prefix/next-label pairs. Pairs whose next label is an axis
contribute only \(O(N^2)\) through \(m\le N\).
Indeed, for next label \(y^a\), the cheaper set
\(\{(i,j):i+\rho j<a\}\) increases by inclusion as \(\rho\)
decreases, starts with at least \(a\) members, and has at most
\(N-a+1\) different values of size at most \(N\).
Sum over \(1\le a\le N\) and both axes to obtain at most \(N(N+1)\).

The published \(C_2(m)=\Theta(m\log(m+1))\), summed over
\(N/2\le m\le N\), therefore leaves
\(\Omega(N^2\log(N+1))\) interior pairs. By (B9) their lengths
lie between \(N/2\) and \(2N+1\). Averaging over those lengths
proves the lower half of (B4); (B1), with \(h=2,a=1,r=n-3\),
proves its upper half.
This is an infinite actual-ideal argument, not a fitted finite law.

## Verification and priority boundary

Run

~~~bash
python3 -B scripts/run_bounded_coordinate_wall_checks.py
~~~

The report is results/research/bounded_coordinate_wall_checks.json.
It replays translated-cover token paths in four quadrants, checks
the residual covers on complete actual fans, and independently
sorts the planar prefix/next-label pairs and their actual slab
wall endpoints. Small characteristic-zero determinants verify
the odd-length extension as well.

The classical inputs are greedy matroid exchanges and their
convex token paths, the integer-direction estimate, the existing
compulsory-core reduction, and the ordered Hermite slab realization.
The planar corner-cut lower count is due to
Corteel--Rémond--Schaeffer--Thomas, as restated in
[Elizalde--Galván, Theorem 2.3](https://arxiv.org/html/2312.16353v1).
The translated cover is an extension of the already recorded dyadic
charge, not a new general parametric-matroid exponent.
Exact priority of (B1)--(B4) is unassessed.

Here (B2) is a restricted corollary of the residual shadow estimate,
not an independent unrestricted improvement. Optimal dependence on
a growing \(h\) is not asserted. The essential conclusions are the
residual \(r\) in (B1) and the sharp fixed-direction lower statement
(B3)--(B4).

The remaining fixed-degree problem concerns the two full-support
sign patterns stated after (B2). The unrestricted problem also has
unbounded coordinate degrees. Neither can be inferred from the
present restricted sums, from a linear residual cap, or from the
number of all tied bases at a wall. The original extremal endpoint
and publication gates remain unchanged.
