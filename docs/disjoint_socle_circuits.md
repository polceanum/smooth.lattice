# Independent socle circuits cannot produce a supergeneric fan

The same two upper monomials can participate in arbitrarily many
persistent ordinary crossings. An explicit family below has
\(\Theta(n^{1/3})\) such crossings, all in a remaining mixed-sign
wall family. Even convexity of the two threshold curves does not
prevent recrossing.

There is also a uniform reason to stop this construction route.
For every field, the entire class of independent relations with
disjoint supports among monomial socle elements has
\[
                 g(I)=O(n^{4/3})                        \tag{S1}
\]
in three variables. This order is attained. If the monomial carrier
has a fixed number of \(x\)-layers, the bound becomes linear in
\(n\), and linear order is attained already with two layers.
Both are below the generic \(n^2\log^2 n\) target scale.

The upper bounds are elementary consequences of established
Minkowski-sum, normal-fan and antichain tools. We use them as a
feasibility closure, without treating a classical corollary or an
unassessed construction as a new publication result. The stronger
[existing socle theorem](initial_cone_stability.md#a-whole-class-bound-socle-only-relations-have-small-fans)
already bounds arbitrary overlapping socle relation spaces by
\(O(n^{3/2})\) over every field. Overlap alone is therefore not an
open separation route. Arbitrary quadratic-coordinate quotients and
relations with surviving nonconstant multiples are outside these socle
closures.

## Exact class and its state polytope

Let \(D\subset\mathbb N_0^3\) be a finite nonempty downset, of size
\(N\), and \(J_D\) its complementary monomial ideal. Choose pairwise
disjoint sets \(E_1,\ldots,E_r\) of maximal elements of \(D\), each
with at least two elements, and nonzero coefficients over \(K\). Set
\[
 I=J_D+(f_1,\ldots,f_r),\qquad
 f_j=\sum_{a\in E_j}c_{j,a}x^a,\qquad c_{j,a}\ne0,
 \quad M=\sum_j|E_j|.                                   \tag{S2}
\]
The classes in each support lie in the socle of \(K[x,y,z]/J_D\).
Thus the span of the \(f_j\) is already an ideal; disjoint supports
make the relations independent. The quotient length is
\[
                       n=N-r,\qquad N\le2n,             \tag{S3}
\]
since \(r\le M/2\le N/2\). Its column matroid is the direct sum of
coloops outside the supports and one circuit on each \(E_j\).

For every global order, a standard set is exactly \(D\) with the
largest monomial of each \(E_j\) deleted. Indeed the proposed initial
ideal contains \(J_D\) and these heads and has colength \(n\).
Its inclusion in the actual initial ideal therefore proves equality.
All proposed sets are downsets, as only maximal elements were removed.

Write \(P_j=\operatorname{conv}(E_j)\). The basis-profile polytope is
exactly
\[
                     \mu(D)-\sum_jP_j.                  \tag{S4}
\]
This is a translate of the Minkowski sum of the reflected \(P_j\),
not subtraction of support functions with an unspecified sign.
Its minimum-cost function is
\[
        w\cdot\mu(D)-\sum_j\max_{a\in E_j}w\cdot a.     \tag{S5}
\]
There is no characteristic-dependent rank test in this formula.

## Counting the fan and a small universal basis

A three-dimensional polytope on \(m\) vertices has at most \(3m\)
edges. Polygons and segments satisfy the same coarse bound. Hence the
normal fans of all \(P_j\), restricted to the positive weight
triangle, have altogether at most \(3M\) straight wall portions.
Their common refinement has \(O(M^2)\) cells: two noncollinear
portions meet at most once; collinear overlaps can be merged, and
incremental insertion gives the usual quadratic arrangement bound.
By (S5), this refinement is precisely the required minimum-basis
fan, up to harmless coincident boundaries. Thus
\[
                              g(I)=O(1+M^2).            \tag{S6}
\]

Here is a self-contained antichain budget with the required dependence
on quotient length. For an antichain \(S\subset D\), fix an integer
\(R\ge1\). At most \(R^2\) points have all coordinates less than
\(R\): projecting away one coordinate is injective on an antichain.
Assign each other point a coordinate at least \(R\), and move it
back zero through \(R-1\) steps along that coordinate. The resulting
points lie in \(D\). Rays assigned the same coordinate are disjoint,
since any collision would make their endpoints comparable. A point
can be counted at most three times, so
\[
                   |S|\le R^2+3N/R=O(N^{2/3})            \tag{S7}
\]
on taking \(R=\lceil N^{1/3}\rceil\). Applying this to the support
union gives \(M=O(n^{2/3})\), proving (S1).

The same estimate bounds the minimal generators of \(J_D\): they
are an antichain of minimal excluded points. For a large coordinate,
move back one through \(R\) steps instead; all these proper divisors
lie in \(D\), and the same ray argument applies.
The set
\[
                         G(J_D)\cup\{f_j:1\le j\le r\} \tag{S8}
\]
is a universal Gröbner basis, of size \(O(n^{2/3})\).
Every reduced Gröbner basis consists of some monomial generators of
\(J_D\) and the \(f_j\) normalized at their chosen heads. Distinct
supports do not reduce one another, and all support monomials are
maximal. Thus the union of reduced bases has at most
\(|G(J_D)|+M=O(n^{2/3})\) polynomials, and at most
\(|G(J_D)|+\sum_j|E_j|^2=O(n^{4/3})\) terms.
These are sparse-output bounds, not bit-complexity claims.

If \(D\) has at most \(h\) layers in \(x\), let \(M_i\) count
support points in layer \(i\). A planar downset containing an antichain
of size \(M_i\) contains at least \(M_i(M_i+1)/2\) points. Hence
\[
                    M^2\le2hN\le4hn.                   \tag{S9}
\]
It follows that \(g(I)=O(hn)\). The layer presentation of \(J_D\)
also bounds its minimal generators by \(O(\sqrt{hN}+h)\), using
the planar antichain estimate and Cauchy--Schwarz. The reduced-basis
union then has \(O(\sqrt{hn}+h)\) polynomials and \(O(hn+h)\)
terms. At each fixed \(h\ge2\), the linear fan power is sharp below.

## Sharp three-variable order

For \(L\ge1\) and \(L<a,b\le2L\), put
\[
 v_{ab}=(2a,2b,9L-2a-2b),\qquad
 z_{ab}=(a,b,9L-a-b+1).                                  \tag{S10}
\]
All \(2L^2\) vectors form an antichain. The \(v\)'s have total
degree \(9L\), the \(z\)'s degree \(9L+1\). No \(v\) divides
any \(z\), since its first coordinate is greater than \(2L\);
the opposite divisibility is impossible by total degree.
Let \(D\) be their divisor downset and impose the independent
binomials \(x^{v_{ab}}-x^{z_{ab}}\).
Then \(n=|D|-L^2=\Theta(L^3)\): a single divisor box gives a cubic
lower bound, and a box of side lengths \(O(L)\) contains \(D\).

All comparison directions are primitive and equal to
\((a,b,-a-b-1)\). Normalize \(w_z=1\). Their \(L^2\) lines are
\[
                  a(w_x-1)+b(w_y-1)=1.                 \tag{S11}
\]
They have \(\Theta(L^4)\) distinct ordinary intersections at positive
weights. To prove the lower bound, choose coprime integers
\(A,B\) in \((L/2,3L/4]\). For every allowed coefficient pair
\((a,b),(a+A,b-B)\) in the grid, their intersection has
\[
 w_x=1+\frac{B}{aB+bA},\qquad
 w_y=1+\frac{A}{aB+bA},\qquad w_z=1.                    \tag{S12}
\]
It is strictly positive. There are \((L-A)(L-B)=\Theta(L^2)\)
such pairs for each difference. The coefficient line contains no
third grid point: coprimality makes its integer points successive
steps of \((A,-B)\), and \(2A,2B>L\) exclude another step within
the grid. Equivalently, the intersection in (S11) has multiplicity
exactly two.

There are \(\Omega(L^2)\) coprime differences in that square.
For an interval of length \(H=\Theta(L)\), the number of pairs
sharing a divisor is at most
\(H^2\sum_{k\ge2}k^{-2}+O(L\log L)\).
The sum is at most \(1/4+\int_2^\infty t^{-2}dt=3/4<1\).
The complement therefore has quadratic order for all sufficiently
large \(L\). This proves the claimed intersection count and
\[
                       g(I_L)=\Theta(L^4)
                                    =\Theta(n^{4/3}).   \tag{S13}
\]
All binomials occur in every reduced basis, so its union has
\(\Omega(L^2)=\Omega(n^{2/3})\) elements too.

The lower order can be extended to all sufficiently large lengths:
choose \(L\) proportional to \(n^{1/3}\) small enough that the
constructed length is at most \(n\), then add coloop monomials on
a coordinate axis. Every paired endpoint has all coordinates
positive, so those added axis points dominate none of the supports.
The socle circuits and full fan are unchanged.

## Sharp linear order with two layers

Fix \(k\ge1\) and put \(d=4k\). Start with both layers equal to
the planar triangle of total degree at most \(d-1\). For
\(0\le i<2k\), put \(\delta_i=0\) if \(i<k\), and one otherwise,
and add the divisor sets of the lower tips
\[
        a_i=(0,2i,4k+\delta_i-2i).
\]
Pair them with the upper socle tips
\(b_i=(1,3k-1-i,k+i)\). All supports are disjoint maximal
monomials. The lower triangle gains \(4k\) points, so the length
after the \(2k\) identifications is
\[
                              n_k=16k^2+6k.             \tag{S14}
\]
The lines have equations
\[
 t=1+\delta_i+(3k-1-3i)(1-u).                            \tag{S15}
\]
Each group of \(k\) lines is a pencil, centered at \((u,t)=(1,1)\)
or \((1,2)\). The first group's slopes \(3k-1-3i\) are positive,
the second's negative. Every cross-group pair therefore meets at a
different positive point with \(u<1\) and \(1<t<2\). There are no
triple intersections away from the two pencil centers. The complete
region count is
\[
                          g(I_k)=k^2+4k-1=\Theta(n_k).   \tag{S16}
\]
The reduced-basis union has at least \(2k=\Theta(\sqrt{n_k})\)
elements, attaining that polynomial-count order as well.
Padding the lower \(y\)-axis preserves the two-layer carrier and all
paired tips, giving the linear lower order at arbitrary large lengths.

## One upper pair recrosses at a polynomial rate

Let \(s\ge1\), \(B=4s^2\), and for \(1\le i\le2s\) put
\[
 p_i=(0,2i+\mathbf1_{i\text{ even}},
                       B-i^2+\mathbf1_{i\text{ odd}}),
 \quad e=(1,1,B),\quad f=(1,0,B+1).                     \tag{S17}
\]
The \(p_i,e,f\) are an antichain. Let \(D\) be their divisor
downset, and in its monomial quotient impose
\[
            x^e=\sum_{i\text{ even}}x^{p_i},\qquad
            x^f=\sum_{i\text{ odd}}x^{p_i}.             \tag{S18}
\]
These are two disjoint socle circuits, over every field, with
\(x^2=0\). Direct row summation gives
\[
                       n_s=(32s^3+30s^2+10s+12)/3.      \tag{S19}
\]
Their two upper-membership threshold curves are the convex functions
\[
 T_e(u)=\max_{i\text{ even}}(2iu-i^2),\qquad
 T_f(u)=\max_{i\text{ odd}}(2iu-i^2).                    \tag{S20}
\]
The upper monomial of a circuit is selected precisely when
\(t<T_e(u)\) or \(t<T_f(u)\), respectively. At each
\[
        u_i=i+\tfrac12,\qquad t_i=i^2+i,
                                  \quad1\le i<2s,       \tag{S21}
\]
exactly the neighboring indices \(i,i+1\) maximize the combined
integer parabola. One belongs to each parity group; each is the
unique maximizer of its own group in a neighborhood. Thus the same
literal upper pair \(\{e,f\}\) has \(2s-1=\Theta(n_s^{1/3})\)
distinct persistent ordinary crossings. The matching partners are
\(p_i,p_{i+1}\), and the two wall directions are
\[
                      (1,-2i,i^2),\quad
                      (1,-2(i+1),(i+1)^2).              \tag{S22}
\]
Both have opposite nonzero \(y,z\) signs, the remaining family from
the earlier coordinate budget. The literal upper labels depend on
\(s\); within each individual quotient the very same pair is reused.
This refutes a uniform constant or polylogarithmic bound on repetitions
per upper pair. It does not assert a growing recurrence count for one
fixed pair of exponent labels independent of \(s\).

The whole fan still has exactly \(4s\) cones. There is one basis
with both upper tips excluded. There are \(s\) bases with only \(e\)
included and \(s\) with only \(f\) included. When both are included,
the active even and odd maxima change at the successive integers
\(2,3,\ldots,2s-1\), giving \(2s-1\) bases, all present at
positive weights. Their sum is \(4s\). Thus polynomial recrossing
alone is not a large-fan mechanism, even in the remaining sign family.

## Decision and prior-art boundary

Park independent disjoint-support socle circuits as a separation route.
The whole class is uniformly subgeneric; two-layer examples are only
linear. Recrossing is real but does not change this decision. A counting
argument by upper pairs must control weighted reuse, while a separating
construction needs non-socle propagation or another mechanism outside
the stronger arbitrary-socle closure. Coupling supports only inside the
socle cannot suffice.

The established state-polyhedron dictionary is in
[Babson--Onn--Thomas, Theorem 2.4](https://arxiv.org/pdf/math/0207135).
Independent socle moves are already in
[Altmann--Sturmfels, Propositions 1 and 3](https://www.math.fu-berlin.de/altmann/PAPER/graph.pdf),
as compared in the [primitive-wall audit](primitive_wall_spectrum.md).
The Minkowski formula for disjoint circuits, planar arrangement bound,
three-dimensional polytope edge bound, and elementary antichain/ray
estimate are the complete mechanisms of (S1) and (S9). We do not
present them as a new general polytope or matroid theorem. Exact priority
of the explicit sharp families and the same-pair recurrence is unassessed;
an unsuccessful search does not establish novelty. No new manuscript or
publication claim is based on this closure.

The full \(G_3\) goal, its quantitative success criterion and the
strongest-current literature assessment remain open and unchanged.

## Exact verification

~~~bash
python3 -B scripts/run_socle_circuit_checks.py
~~~

The report is `results/research/socle_circuit_checks.json`. The checker
constructs actual residue columns and commuting multiplication maps,
contracts only literal coloops, and traverses whole cones from full
border normal forms. It compares those cones with the independent
support-maximum overlay and records exact witnesses for every cone.
It checks the actual reduced-basis union against (S8), verifies the
ray budget for supports and minimal excluded monomials, and independently
counts positive grid intersections and their multiplicities. Finite
checks support the displayed uniform proofs; they do not establish
current novelty or the unrestricted endpoint.
