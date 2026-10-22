# Two-colour thresholds improve the overlap-rank bound

Let a rank-\(r\) matroid have a ground partition \(E_0\sqcup E_1\),
with \(m_i=|E_i|\), and costs
\[
 c_e(u,t)=a_eu+b_e+\mathbf1_{e\in E_1}t.
\]
Write
\[
 h=r(E_0)+r(E_1)-r,\qquad
 P=\binom{m_0}{2}+\binom{m_1}{2}.
\]
The number of full-dimensional regions with different uniquely optimal
bases, or different exposed profiles, is
\[
                         O\bigl((h+1)(P+1)+h^2\bigr).       \tag{1}
\]
The bound holds on the whole real \((u,t)\)-plane, hence on any convex
subdomain. No representability or integer labels are required.
In particular two labelled bases of a rank-\(h\) matroid have
\(O((1+h)^3)\) optimal regions in this cost family.

Applied to the existing [four-quotient reduction](quadratic_intersection_fans.md),
this improves its overlap-rank factor from \(h^{10/3}\) to \(h^3\):
for a length-\(n\) quotient containing a scalar monic quadratic in one
coordinate,
\[
       g(I)=O\bigl(n(1+\log(n+1))(1+h)^3\bigr).             \tag{2}
\]
The established \(O(n^{7/3})\) class bound and
\(O((h+1)(1+n^{4/3}))\) slice bound should still be used when smaller.
Thus (2) does not improve unrestricted \(G_3\), the sharp fixed-\(h\)
order, or the strongest saved necessary overlap for a polynomial
separation. Exact priority of (1)--(2) remains unassessed.

## One within-colour tie has one interval of eligibility

First perturb the coefficients generically, within the displayed
two-colour cost family. Same-colour equalities occur at distinct
\(u\)-values, no three relevant threshold curves meet away from them,
and coincidences not forced by the matroid are absent. A finite list
of previously strict optimal witnesses survives a sufficiently small
perturbation. Therefore an upper bound in this case also bounds
exposed profiles and strict regions in degenerate data.

At one such event \(u_0\), let \(e,f\in E_0\) have equal cost
\(c\). Let \(C\) contain all strictly cheaper lower-colour elements.
As \(t\) varies, the upper elements cheaper than \(c\) are successive
prefixes \(D_j\) of one fixed ordered list. Put
\(N_j=M/(C\cup D_j)\). Greedy choice exchanges \(e,f\) when their
order is reversed exactly if they are a nonloop parallel pair in
\(N_j\):
\[
 r(CD_je)-r(CD_j)=r(CD_jf)-r(CD_j)
               =r(CD_jef)-r(CD_j)=1.                       \tag{3}
\]
Here juxtaposition denotes union, not a product of ideals.

The indices satisfying (3) form one consecutive interval, possibly
empty or unbounded at one end. Contraction only decreases ranks on
\(\{e,f\}\). Once the pair is nonloop parallel, a later contraction
either keeps it so or makes both elements loops. It cannot become
independent again, and one member cannot become a loop while the
other remains nonloop. A member that was already a loop cannot
later become nonloop. These facts prove the interval assertion for
arbitrary matroids.

Thus the vertical wall where this same-colour exchange is active has
one connected support interval in \(t\). It can be subdivided by
other walls; the statement concerns its support, not its adjacency
multiplicity. The same argument applies to a tie in \(E_1\), with
the direction of the prefix sweep reversed.

## There are at most two threshold bends per event

Inside an open interval of fixed internal orders, greedily compress
each colour to an ordered independent basis. The established two-list
rank identity gives \(h\) matched pairs and a fixed unmatched set.
The optimal cost is
\[
  \ell(u,t)+\sum_{(e,f)\text{ matched}}
                         \min(a_eu+b_e,a_fu+b_f+t).        \tag{4}
\]
Its nonvertical walls are the \(h\) threshold lines
\(t=(a_e-a_f)u+b_e-b_f\).
This use of two-list greedy matching is classical; it is not a new
matroid representation theorem.

At \(u_0\), away from the two endpoints of the interval (3), a small
neighbourhood either has no same-colour exchange or has that exchange
throughout a neighbourhood in \(t\). In the latter case, reversing
the tied pair changes greedy cost by \(c_e-c_f\), independently of
the local value of \(t\). Consequently every threshold curve there
continues with the same affine formula. Other threshold curves may
cross the vertical wall, but they do not bend there.

For completeness, subtract the two support functions on either side
of \(u_0\). In a small interval of \(t\) away from an endpoint,
their greedy bases differ by either nothing or precisely
\(e\leftrightarrow f\). Their difference is therefore zero or
\(c_e-c_f\). The coefficient of \(t\) and its breakpoints are
unchanged. A breakpoint can change affine slope only at an endpoint
of the exchange interval. This proves the stated continuation rule.

The \(h\) lines can therefore be joined across successive order
events into \(h\) continuous piecewise-affine graphs, continuing
straight through ordinary intersections. They have at most
\[
                              B\le2P                       \tag{5}
\]
bends in total. Along the entire \(t\)-axis there are exactly \(h\)
thresholds counting multiplicity: the feasible upper-colour counts
are the consecutive integers \(r-r(E_0),\ldots,r(E_1)\).
Equation (4) also proves this without an integrality assertion about
arbitrary projected sections.

## Count crossings after counting bends

Two piecewise-affine graphs with \(B_i,B_j\) bends have at most
\(1+B_i+B_j\) proper intersections: partition the real axis at all
their bends and solve one affine equality on each remaining interval.
Genericity removes overlaps and intersection/bend coincidences.
Summing gives
\[
             X\le\binom h2+(h-1)B.                         \tag{6}
\]
The arrangement of the \(h\) complete graphs has at most
\(1+h+X\) regions. Add the at most \(P\) vertical exchange support
intervals. Each meets every threshold graph at most once, so it
increases the number of regions by at most \(h+1\). Hence, for
\(h\ge1\), the following deliberately loose explicit bound suffices:
\[
  1+h+\binom h2+(3h+1)P.                                  \tag{7}
\]
For \(h=0\), the two colour spans are independent, and at most
\(1+P\) vertical regions remain. This proves (1). Restricting a
region of strict optimality to a convex subdomain cannot split it,
since such a region is an intersection of affine halfplanes.

In the four-quotient application, each open planar strip carries
one fixed rank-\(h\) matroid on two \(h\)-element bases, with original
binary first-coordinate labels. There are \(O(n\log(n+1))\) such
strips. Apply (1) once on each, obtaining (2). The common weight list
for **all** compatible matrices still uses the earlier comparison
arrangement; (1) counts the fan of one fixed matrix, not one sampler
that works for every matrix.

## Verification and boundary of the result

~~~bash
python3 -B scripts/run_two_colour_braid_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker uses exact union ranks to recover every matching, checks
the interval (3) independently at each event, follows all threshold
pieces and compares their regions with a separate enumeration of
bases and strict halfplane feasibility, and every original quotient
S-pair. It includes represented
matroids in several characteristics, a nonrepresentable matroid,
loops, parallel columns and actual quotient-column restrictions.
The abstract eight-element example independently satisfies basis
exchange and violates Ingleton's inequality by \(15<16\), so it
requires no assumed linear representation. The classical inequality
is displayed in [Dougherty--Freiling--Zeger, introduction, page 2](https://arxiv.org/pdf/0910.0284v3);
that statement was inspected, without auditing their later
five-variable classification. Finite controls test the proof; they do not establish its uniform
scope, a sharp cubic lower bound or exact novelty.

The [two-list proof](rational_shadow_bound.md#two-lists-give-a-partial-matching)
records the comparison with Gusfield's two-ordered-set matroid problem.
[Hausbrandt--Ruzika, Corollary 3.4](https://arxiv.org/html/2503.08178v2#S3)
gives \(O(m^{2p})\) regions for unrestricted \(p\)-parameter matroids;
its statement and pairwise-arrangement proof were inspected. At two
parameters it gives \(O(m^4)\), rather than the two-colour bound above.
That comparison is not a complete priority audit.

The unrestricted collective obligation is still to control actual
crossings across the moving encodings. A cubic overlap factor alone
does not provide a cubic global wall budget, a polynomially
supergeneric actual family, or completion of the extremal goal.
