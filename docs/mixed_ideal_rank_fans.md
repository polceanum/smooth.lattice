# Initial fans controlled by the rank of an actual mixed ideal

Write \(g(I)\) for the number of distinct monomial initial ideals under
global monomial orders, and \(\mathcal U(I)\) for the union of their
reduced Gröbner bases. The quantity \(V(I)\) sums the primitive wall
cost \(\nu(q)=\prod_i(q_i^++1)+\prod_i(q_i^-+1)-2\) over undirected
adjacent maximal cones in the positive weight fan; here
\(q_i^\pm=\max(\pm q_i,0)\) and \(q\) is a primitive wall normal.

Let \(A=K[x_1,x_2,x_3]/I\) have positive finite length \(n\), over
any field. Choose arbitrary monic univariate polynomials \(f_i(x_i)\)
of positive degrees \(d_i\), and put
\[
 T=(f_1f_2,f_1f_3,f_2f_3)A,\quad t=\dim_K T,\quad
 S=d_1d_2d_3,\quad R=d_1d_2+d_1d_3+d_2d_3,\quad
 L_t=1+\log(t+1).
\]
No vanishing assumption is imposed. In particular \(T\) may be nonzero,
nonreduced, or all of \(A\). There is a uniform bound
\[
 g(I)=O\!\left((1+nR+S)^2(1+S+tL_t^2)^4\right).           \tag{M1}
\]
The constant is absolute. At fixed degrees and bounded \(t\), this
gives quadratic fan and cubic total wall-cost upper orders, extending
the [zero mixed-ideal theorem](univariate_cross_fans.md).
The quadratic and cubic powers are attained within the class with
bounded \(t\), since its \(t=0\) subfamily already attains them.
No sharp dependence on a positive value of \(t\) is asserted.

For linear factors there are no internal carrier-order changes. The
stronger statement, independent of \(n\), is
\[
 \boxed{g(I)=O\!\left((1+t)^4L_t^8\right)}
 \qquad\left(f_i=x_i-a_i,\ a_i\in K\right).               \tag{M2}
\]
This applies to every finite quotient and every choice of center.
Thus a family with \(g(I_n)\ge c n^{2+\epsilon}\) necessarily has
\[
 \dim((x-a)(y-b),(x-a)(z-c),(y-b)(z-c))A_n
   =\Omega_{c,\epsilon}\!\left(
             \frac{n^{1/2+\epsilon/4}}{\log^2(n+1)}\right) \tag{M3}
\]
at every center \((a,b,c)\in K^3\).
One may equivalently use the smallest of these integer dimensions
over all centers. This is an unrestricted necessary condition, not a
new worst-case upper exponent for \(G_3(n)\).

Every initial ideal also has
\[
 O\!\left(S+R+(1+t)^{2/3}\right)                           \tag{M4}
\]
minimal generators. In the linear case the union of reduced bases
therefore has \(O((1+t)^{14/3}L_t^8)\) polynomials.
These are polynomial-count bounds. Degrees, term counts and coefficient
bits may still depend on \(n\). The total primitive cost satisfies
\(V(I)\le3n\,g(I)\); no \(n\)-independent cost bound is inferred.

The unrestricted sharp \(G_3\) endpoint, an actual polynomial lower
separation, and exact numerical priority remain open.

## A unitriangular polynomial basis

Every exponent label \(\alpha\in\mathbb N_0^3\) has a unique expression
\(\alpha_i=b_i+d_i u_i\), with \(0\le b_i<d_i\). Define
\[
                 Q_\alpha=x^b\prod_i f_i(x_i)^{u_i}.      \tag{M5}
\]
This polynomial is \(x^\alpha\) plus a combination of its proper
divisors. Conversely, successive monic division expresses \(x^\alpha\)
as \(Q_\alpha\) plus such lower labelled columns.
Both changes are unitriangular on every finite monomial downset.
They preserve the span of every positive weight prefix, over every
field, including inseparable factors.

Consequently the greedy scan of the \(Q_\alpha\)'s, priced by
\(w\cdot\alpha\), has exactly the original monomial acceptance labels.
If \(B\) is an initial standard set, its \(Q_B\) columns are a basis
and, at every strictly exposing positive objective, uniquely minimize
this labelled cost among polynomial-column bases. Any finite ground
containing \(B\) retains that minimum.
This is a filtration-preserving change of vector-space basis,
not a substitution of a different quotient algebra.
It does not preserve the ranks of arbitrary labelled subsets of
monomial columns: an individual nonzero monomial may become a zero
\(Q\)-column. Only the prefix and downset assertions above are used.

## At most \(t\) mixed labels and a small common ground

Let
\[
 J=(x_1^{d_1}x_2^{d_2},x_1^{d_1}x_3^{d_3},
                                      x_2^{d_2}x_3^{d_3}).
\]
For every size-\(n\) order-ideal monomial vector basis \(B\), whether
coherent or not,
\[
                              |B\cap J|\le t.             \tag{M6}
\]
Indeed \(Q_B\) is a basis by unitriangularity, and every \(Q_\alpha\)
with \(\alpha\in J\) belongs to the actual ideal \(T\).
These selected columns are independent in a space of dimension \(t\).
This does not assert that the ordinary monomial columns in \(B\cap J\)
belong to \(T\).

Let \(h_{ij}=d_i e_i+d_j e_j\), and use
\[
 H_t=\{v\ge0:\prod_i(v_i+1)\le t\},\qquad
 E_T=\bigcup_{i<j}(h_{ij}+H_t),                           \tag{M7}
\]
with \(E_T=\varnothing\) at \(t=0\). Every mixed label of every such
basis lies in \(E_T\). To prove this, choose \(h_{ij}\le\alpha\).
The whole interval \([h_{ij},\alpha]\) consists of mixed divisors of
\(\alpha\) and is contained in \(B\). Its cardinality is at most \(t\)
by (M6), proving the product bound in (M7).
In particular
\[
                        q:=|E_T|\le3|H_t|=O(tL_t^2).      \tag{M8}
\]
There is no factor \(S\) in this count. Ordinary divisor closure of
the standard labels supplies this stronger ground.

Use four carrier groups from the ordinary universal length-\(n\)
ground: the box \(\alpha_i<d_i\), and the three tails with exactly
one coordinate at least its \(d_i\). Add \(E_T\) as a fifth group,
using \(Q\)-columns throughout. This finite ground contains every
initial basis. Some body labels may lie outside the ordinary universal
ground; including them cannot improve the global labelled minimum
established above. Their span has rank at most \(t\), and need not
equal the whole of \(T\).

## Bounded corank modulo an actual algebra ideal

Pass to \(\bar A=A/T\), of length \(n-t\), allowing zero length.
In this algebra the \(f_i\)'s annihilate each other pairwise.
The branch-generation and rank-nullity proof from the cross theorem
gives, for \(\bar Y_i=f_i\bar A\),
\[
 \sum_i\dim\bar Y_i-\dim\sum_i\bar Y_i\le2S.
\]
The core has rank at most \(S\), and each tail \(Q\)-column is
directly in the corresponding branch. If \(m_i\) are the four
carrier-group ranks in \(A\), lifting ranks through \(A\to\bar A\)
therefore gives
\[
                           \sum_{i=0}^3m_i\le n+3S+3t.    \tag{M9}
\]
All four quotient images span \(\bar A\): body columns vanish there,
and the chosen ground contains an entire basis of \(A\).

Fix the four carrier internal orders and greedily compress those groups.
Keep every body column; no body-order arrangement is needed.
Their combined rank is \(n\), with \(M=\sum m_i+q\) columns and corank
\[
                         \delta=M-n\le3S+3t+q.           \tag{M10}
\]
Use the rank-\(\delta\) dual. Descending greedy compression within
each of the four carrier lists retains at most \(\delta\) candidates
per list. Keep the \(q\) body elements in the dual candidate set too.
All maximum dual bases, hence all minimum original bases, are captured
on at most
\[
                           4\delta+q=O(S+tL_t^2)          \tag{M11}
\]
labels. Greedy rejection and duality are classical; the actual ideal
rank is the additional quantitative input.

As before, the four carrier internal orders change on
\(O(nR+S)\) planes. Their positive arrangement has quadratically many
chambers. Inside each, (M11) gives at most
\(O((1+S+tL_t^2)^4)\) cells by the pair-comparison arrangement.
This proves (M1).

## A global finite matroid encoding for linear factors

If \(d_1=d_2=d_3=1\), the carrier groups are \(1\) and the three
positive pure-power rays. Their positive cost orders never change.
The two compressions above are then fixed for the entire positive
weight cone. They give one dual matroid on \(O((1+t)L_t^2)\) labelled
elements. Its maximum-cost basis determines the original standard set
by complementation in one fixed compressed union.
This is an exact finite matroid encoding of the whole fan.
It is not an algebra quotient of dimension \(O(t)\).
An arrangement of its pair-comparison planes proves (M2).
Using \(t\le n\) to bound the logarithm and rearranging proves (M3).

From commuting multiplication matrices and a cyclic vector, at fixed
\(t\) this linear-factor encoding and all original reduced bases can
be constructed in \(O_t(n^3)\) field-and-geometric-arithmetic operations
by dense elimination, a dual nullspace, and the constant-size candidate
arrangement. Preparing \(T\), the ground and the original border
normal forms is included at that fixed parameter. Bit complexity is
excluded. No faster-than-known implementation or optimal-time claim
is made. At fixed degrees and \(t\), the corresponding chamber method
has the earlier \(O_{\boldsymbol d,t}(n^5)\) arithmetic envelope.
The verification script recomputes ranks; it does not implement the
optimized fixed-parameter arithmetic construction.

## Minimal generators with a small mixed part

A minimal excluded monomial outside \(J\) lies in the finite core or
at one of the \(R\) carrier ray frontiers. There are at most \(S+R\)
such monomials, even though \(B\) may also contain body labels.

For a minimal excluded monomial \(\alpha\in J\), choose one
\(h=h_{ij}\le\alpha\) and form the downset
\[
                         D_h=\{v\ge0:h+v\in B\}.
\]
It has size at most \(t\), and \(\alpha-h\) is a minimal excluded
monomial of \(D_h\). An ordinary three-dimensional downset of size at
most \(t\) has \(O(1+t^{2/3})\) minimal excluded monomials. For
\(t>0\), set \(h=\lceil t^{1/3}\rceil\). For each coordinate, the heads
whose value in that coordinate is below \(h\) lie on \(h\) planar
slices. The elementary planar boundary bound and Cauchy--Schwarz give
\(O(h+\sqrt{ht})\) heads across these slices. For the remaining heads,
the \(h\) predecessors obtained by decreasing the first coordinate
are inside the downset. Those chains are disjoint, since two minimal
heads with the same other coordinates would be comparable. There
are at most \(t/h\) such heads. The total is \(O(t^{2/3})\).
The empty downset is handled separately. Summing over the three
anchors proves (M4). The usual union over initials proves the stated
reduced-basis cardinality bound.

## Verification and boundary

Run

~~~bash
python3 -B scripts/run_mixed_ideal_rank_checks.py
~~~

The report is results/research/mixed_ideal_rank_checks.json.
It separately checks ordinary monomial columns, polynomial expansion
columns and evaluated-operator columns; full positive prefix ranks;
unitriangular determinants on downsets; the actual mixed ideal's
multiplication invariance; quotient branch and corank bounds; and
the finite body-ground cover. Full dual cone patches are compared with
original border normal forms, and merged complete fans with an
independent direct traversal. Positive-error examples include translated
centers, inseparable factors, nonlocal points and a coherent staircase
which actually selects a body label.
Explicit guards also exhibit changed singleton ranks and ordinary
monomial body columns outside \(T\), while the corresponding
polynomial body columns remain in \(T\).

The monic free-module basis, unitriangular filtration argument,
rank-nullity, translated divisor cover, matroid duality and arrangements
are established tools. The
[Babson--Onn--Thomas state framework](https://arxiv.org/pdf/math/0207135)
and the greedy setup in
[Eppstein, Sections 1.2 and 4](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf)
are prior inputs. The two-parameter \(O(m^4)\) arrangement bound is
also stated in
[Hausbrandt--Ruzika, Corollary 3.4](https://arxiv.org/html/2503.08178v1);
that counting step is classical. The issue here is the reduction to
\(O((1+t)L_t^2)\) labels. Exact numerical priority of the mixed-rank bound
and its consequences remains unassessed.

The rank of an actual multiplication-stable ideal is essential.
A compulsory-core span cannot replace \(T\), and the number of pair
generators is not its dimension. Large mixed rank is necessary under
(M3), not sufficient for a large fan. When \(t\) is proportional to
\(n\), (M2) does not improve the existing unrestricted exponent.
The original sharp endpoint, actual polynomial-separation alternative,
current-priority audit and honest-paper requirements are unchanged.
