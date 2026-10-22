# Sharp cone size, including balanced weights

An individual initial cone can be much smaller than the number of cones
suggests. This note rules out a pointwise projection-sampling shortcut to
the unrestricted counting problem. It strengthens the earlier
[relative-margin examples](initial_cone_stability.md) by keeping all
coordinate weights comparable and shrinking every direction of the cone.
The construction and the cofactor bounds below give sharp orders in
every fixed dimension. Exact priority remains unassessed.

Let \(d\ge2\) be fixed. Normalize positive weights by
\(\sum_i w_i=1\). Measure volume in the coordinates
\((w_1,\ldots,w_{d-1})\), relative to the full weight simplex, and
diameter in the Euclidean metric of this simplex. Over all positive
initial cones of colength-\(n\) ideals, over all fields, the worst orders
are
\[
 \inf_C\operatorname{vol}_{\rm rel}(C)
       =\Theta_d(n^{-d(d-1)}),\qquad
 \inf_C\operatorname{diam}(C)=\Theta_d(n^{-d}).             \tag{B1}
\]
The same orders hold for logarithmic volume and diameter in the chart
\(\log(w_1/w_d),\ldots,\log(w_{d-1}/w_d)\).
For every sufficiently large length, the upper examples have exactly
\(2^d-1\) initial ideals. Their small cone remains inside a fixed
compact subset of the positive simplex. Spherical area and angular
diameter have the corresponding orders as well.

In particular, seven initial ideals in three variables can include a
balanced cone of relative area \(\Theta(n^{-6})\) and diameter
\(\Theta(n^{-3})\). Small cone size is not evidence of a large fan.
No unrestricted improvement in \(G_3(n)\), power separation or completed
research goal follows from (B1).

## The universal lower bounds are cofactor consequences

In a length-\(n\) initial staircase, every standard monomial has degree
at most \(n-1\), and every minimal excluded monomial has degree at
most \(n\). Thus the exponent differences defining its cone have
coordinates of absolute value at most \(n\). Include the coordinate
positivity inequalities. A primitive extreme ray of the closed cone
is obtained from \((d-1)\)-by-\((d-1)\) cofactors, so its nonnegative
integer coordinates are at most \((d-1)!n^{d-1}\).

Choose \(d\) independent such rays \(r_j\), and write
\(s_j=\sum_i r_{j,i}\le d!n^{d-1}\). Their normalized simplex lies
in the cone section, and its relative volume is
\[
 \frac{|\det(r_1,\ldots,r_d)|}{\prod_j s_j}
       \ge(d!n^{d-1})^{-d}.                               \tag{B2}
\]
Its determinant is a nonzero integer. Volume bounded by a constant
times diameter to the power \(d-1\) proves the diameter lower bound.
This argument concerns each full-dimensional positive cone, including
ones whose closure touches the boundary of the weight simplex.

The inverse logarithmic chart is the softmax map and has a globally
bounded derivative in fixed dimension. Its Jacobian in the displayed
simplex coordinates is \(\prod_i w_i\); equivalently the forward
Jacobian is \(1/\prod_i w_i\ge d^d\). Consequently (B2) also gives
the logarithmic lower orders. The sphere/simplex radial maps have
bounded Jacobians and derivatives in both directions on the positive
simplex. These are elementary geometric consequences, not new general
integer-cone or integral-geometric theorems.

## A unimodular family with a balanced kernel

For \(d\ge3\), let \(A\) have rows
\(e_i-e_{i+1}\), with cyclic indices. Its adjugate is the all-ones
matrix. The following fixed integer matrices specify \(B\) explicitly.
Let \(U\) have columns \(e_i-e_d\) for \(i<d\), followed by \(e_d\).
Let \(V_0\) have rows \(e_i-e_{i+1}\) for \(i<d\), followed by \(e_d\).
Both are unimodular. Let \(J\) have ones on its superdiagonal and
zeroes elsewhere. Finally let \(S\) have rows \(e_{i+1}\) for \(i<d\)
and last row \((-1)^{d-1}e_1\), so \(\det S=1\).
Put
\[
 B=US^{-1}V_0,\qquad
 M(k)=B+kA=U(I+kJ)S^{-1}V_0.                              \tag{B3}
\]
The equality follows from \(JS^{-1}=\operatorname{diag}(1,\ldots,1,0)\).
Since \(J\) is nilpotent, \(\det M(k)=1\) for every \(k\).
Its inverse \(R(k)\) is an integral polynomial matrix of degree
\(d-1\), with leading coefficient
\(\operatorname{adj}(A)\), the all-ones matrix. Hence all entries
are positive and \(\Theta_d(k^{d-1})\) for sufficiently large \(k\).

The cone \(M(k)w>0\) lies in the positive orthant, and its rays are
the columns of \(R(k)\). Put \(s_j=\sum_i R_{ij}(k)\).
Their normalized simplex has exact relative volume
\[
                         \frac1{\prod_j s_j}
                           =\Theta_d(k^{-d(d-1)}).        \tag{B4}
\]
All its vertices tend to \((1/d,\ldots,1/d)\). Jacobi's complementary
minor identity bounds every two-by-two minor of \(R(k)\) by
\(O_d(k^{d-2})\), since \(\det M=1\). Dividing the wedge product of
two inverse columns by their lengths therefore bounds their angle by
\(O_d(k^{-d})\). Normalization has the same order on the positive
simplex. Thus the diameter is \(O_d(k^{-d})\), and (B4) supplies
the matching lower order. The logarithmic and spherical maps have
bounded distortion on this balanced region.

## Realization by actual polynomial quotients

Each row \(q_i(k)\) of \(M(k)\) has one long positive coordinate
at \(i\), one long negative coordinate at \(i+1\), and bounded other
coordinates. Write \(q_i=q_i^+-q_i^-\). Choose a fixed integer
\(T>2\max_{ij}|B_{ij}|\), and put
\[
 t_i=Te_{i+2},\qquad u_i=q_i^++t_i,\quad v_i=q_i^-+t_i.  \tag{B5}
\]
For large \(k\), these \(2d\) monomials form an antichain. Endpoints
with different long coordinates cannot dominate one another. At a
fixed long coordinate \(j\), the two endpoints have distinct large
tags in coordinates \(j+2\) and \(j+1\), both different from \(j\).
Their differences dominate all bounded entries of \(B\), proving
incomparability also within that pair.

Let \(D_0\) be the union of their divisor boxes, and let \(J_0\) be
the monomial ideal with standard set \(D_0\). Its length is
\(\Theta_d(k)\): each box has one long side and bounded remaining
sides. All endpoints are in the socle of \(K[x]/J_0\). Therefore
\[
                I_k=J_0+(x^{u_i}-x^{v_i}:1\le i\le d)     \tag{B6}
\]
is an actual quotient of length \(|D_0|-d=\Theta_d(k)\), over
every field. The disjoint endpoint relations have independent rows,
even in characteristic two. All other monomial residue classes are
singletons. An initial standard set chooses the cheaper endpoint
of each pair and keeps every singleton. Thus the all-positive
orientation has exactly the cone \(M(k)w>0\).

Every sign pattern except all negative occurs at positive weights.
Indeed choose slacks with the prescribed signs, making any positive
one sufficiently large, and multiply them by the positive inverse
\(R(k)\). All-negative slacks cannot give a positive weight.
The resulting standard sets are different, proving the exact
\(2^d-1\) count. Products of the relations by every variable vanish
modulo \(J_0\); no unsupported commuting-multiplication assumption
is needed. This is the established socle realization used in the
earlier stability note, applied to (B3).

The downsets \(D_0(k)\) increase with \(k\), with positive gaps in
their lengths bounded by a constant depending only on \(d\). A pure
coordinate-power chain can be extended past the last such power in
\(D_0(k)\) to fill each gap. Every old endpoint has a nonzero tag in
another coordinate, so none divides an added pure power. All old
relations stay in the socle; the added monomials are forced independent
columns and leave the fan unchanged. This gives every sufficiently
large length. Dimension two follows from the existing four-endpoint
construction in the stability note and its padding argument. This
completes (B1) for every fixed dimension.

## A smaller explicit three-variable witness

For \(k\ge4\), the following matrix and its positive inverse suffice:
\[
 M_k=\begin{pmatrix}k-1&-k&0\\1&k-1&-k\\-k&0&k+1\end{pmatrix},
\quad
 R_k=\begin{pmatrix}
 k^2-1&k^2+k&k^2\\
 k^2-k-1&k^2-1&k^2-k\\
 k^2-k&k^2&k^2-k+1
 \end{pmatrix}.                                          \tag{B7}
\]
Use the pairs
\[
\begin{array}{c|c}
u_i&v_i\\\hline
(k-1,1,1)&(0,k+1,1)\\
(2,k-1,1)&(1,0,k+1)\\
(0,1,k+1)&(k,1,0).
\end{array}                                               \tag{B8}
\]
The divisor union has \(13k-6\) elements, so (B6) has length
\(n=13k-9\). Its slices at first coordinates \(0,1,2\) have sizes
\(4k+4,3k,2k\); those at \(3,\ldots,k-1\) have size four, and the
last slice at \(k\) has size two. Their sum proves the count.
The column sums in (B7) are
\[
 s_1=3k^2-2k-2,\quad s_2=3k^2+k-1,\quad s_3=3k^2-2k+1.  \tag{B9}
\]
The small cone has exact relative simplex area \(1/(s_1s_2s_3)\).
Adding \(0\le r\le12\) pure powers \(z^{k+2},\ldots,z^{k+r+1}\)
gives every length at least 43 with the same seven-cone fan.

## Sampling consequence and literature boundary

For independent uniform samples in the positive weight simplex, the
waiting time to see the small cone in (B7) has mean \(s_1s_2s_3\).
Constant success probability requires \(\Theta(n^6)\) samples, even
though there are only seven cones. In dimension \(d\), the analogous
order is \(n^{d(d-1)}\).

An independent isotropically random hyperplane through the origin
hits the interior of the small cone with probability
\(\Theta_d(n^{-d})\). For any two ray directions making angle
\(\theta\), their scalar products with a uniform normal have different
signs with probability \(\theta/\pi\), by rotation symmetry. A plane
hits the cone precisely when its products with the rays have both signs.
One pair gives a lower bound; the union over the fixed number of pairs
gives the upper bound. Thus independent isotropic projection planes
require \(\Theta_d(n^d)\) trials for constant detection probability.
For three rays the exact probability is
\((\theta_{12}+\theta_{23}+\theta_{31})/(2\pi)\), because a mixed
triple of signs has exactly two unequal pairs.

These statements concern the specified sampling distributions. They
do not bound adaptive, algebra-dependent sampling, exact fan traversal,
or the number of planes in a carefully selected cover.
[Dimitrova's 2010 publisher abstract](https://doi.org/10.1007/s11786-010-0039-0)
already describes stochastic estimation of cone volumes; its full
subscription text has not been inspected here. The comparison does
not assign that method an unverified worst-case claim.
[Hassin--Tamir's primary introduction, p. 362](https://www.tau.ac.il/~hassin/matroids_89.pdf)
optimizes convex functions of two basis sums, equivalently scanning
one scalar slope. Its title is not a theorem about two affine cost
parameters. The newer multi-parametric comparison remains the one
recorded in the [two-colour note](two_colour_threshold_braids.md).
No exhaustive priority conclusion follows from these scoped readings.

Run `python3 -B scripts/run_balanced_cone_checks.py` for polynomial
matrix identities, actual multiplication, original reduced-basis
S-pairs, entire seven-cone comparisons and exact geometry. Finite checks
support the proof; they neither establish its asymptotic scope nor
turn the sampling obstruction into a large-fan lower bound. The main
remaining task is collective control of one entire fan, not a uniform
pointwise visibility estimate.
