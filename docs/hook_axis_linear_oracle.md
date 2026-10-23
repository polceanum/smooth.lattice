# A shift inequality makes axis-hook profile queries linear in thickness

For the actual zero-outer axis hooks of the
[birth-curve theorem](hook_axis_birth_curves.md), the compiled profile
oracle needs only
\[
                         O(t+\log A)                         \tag{AL1}
\]
exact rational and integer operations per generic query. A nested-basis
compilation also reduces preprocessing to
\[
                    O(A^3+A^2t^2)\text{ field operations},
                    \qquad O(A^2)\text{ field storage}.
\]
No coefficient field operations occur after compilation.
It also returns a convex open refinement patch using at most \(4t+4\)
integer inequalities. Polygon construction costs additional rational
operations. These are arithmetic bounds, without a bit-complexity or
optimal basis-output assertion.

Here \(A=a-t\ge t\ge2\), the length is \(n=4at-2t^2\), and
\[
 I=(x^2,y^a,z^a,y^tz^t-x(y^tF(y)+z^tG(z)))
\]
over any field. Both polynomials may be nonunits or zero. The linear
query result uses their shift structure; the larger arbitrary replicated
matrix family (AB3) still has the quadratic-head oracle from the earlier
note. The improved compilation applies to that larger family as well.
The actual fan-count bounds and the unrestricted extremal question are
unchanged.

## Multiplication orders the head gains

The earlier birth checker already sampled the normalized Krylov
interlacing below. Its complete shift proof and its linear-merge
consequence are made explicit here.

In one factor put \(V=K[q]/q^A\). The pure vectors are \(q^r\)
of price \(\theta+r\), \(0\le r<A\), and \(q^rW\) of price
\(r\), \(t\le r<A\). Write \(E_\theta(\rho)\) for their span
at prices at most \(\rho\), and \(h_i=q^iW\).
Multiplication by \(q\) gives
\[
 qE_\theta(\rho)\subseteq E_\theta(\rho+1),\qquad
 q\langle h_0,\ldots,h_{i-1}\rangle
                \subseteq\langle h_0,\ldots,h_i\rangle.       \tag{AL2}
\]
The last vector of either pure family maps to zero. Thus no wraparound
or missing boundary column is needed in (AL2).

Use the whole-prefix births \(\alpha_i(\theta)\) from (AB4).
Multiplying a representation at the birth of \(h_i\) proves
\[
                  \alpha_{i+1}(\theta)\le\alpha_i(\theta)+1.  \tag{AL3}
\]
If \(h_i\) is in the earlier head span, multiplication shows the same
for \(h_{i+1}\). Absent births therefore form a final segment. With
absent gains interpreted as \(-\infty\),
\[
 \alpha_0(\theta)\ge\alpha_1(\theta)-1\ge\cdots
                         \ge\alpha_{t-1}(\theta)-(t-1).       \tag{AL4}
\]
This proof is uniform in the coefficients, field and price parameter.
It uses full spans, including every cheaper pure column and the earlier
heads. It does not replace those spans by support matching. Equality
between gains in one arm is allowed.

An arbitrary replicated factor can have zero
lower vectors, \(h_0=U_0\), \(h_1=U_3\). Its two gains are
\(\theta\) and \(\theta+2\); (AL4) fails. Such data have no
asserted axis-hook realization.

## Compile by nested pure bases

This part needs only the two pure price classes of (AB3). The upper
vectors form a basis and have prices \(q+\theta\); lower vectors have
prices \(q\), with distinct integral labels within each class.
At a nonintegral parameter every pure price is distinct. A greedy basis
contains an upper vector exactly when it is independent of all cheaper
vectors. Its cheaper upper vectors stay fixed as \(\theta\) increases,
and its cheaper lower span increases. Thus an upper vector can leave
the greedy basis only once. The cheaper-upper span for a fixed lower
vector decreases, so a lower vector can enter only once. There are at
most \(A-t\) exchanges altogether.

Initialize before all price crossings with the full upper basis. Keep
the coordinates of every pure vector and every head in the current
basis. At an integral crossing \(k\), exchange a retained upper
\(U_r\) for an omitted lower \(L_q\) whenever
\(q-r=k\) and its fundamental coefficient on \(U_r\) is nonzero.
Repeat until no such pair remains. Every exchange has zero cost at
\(k\), so it preserves minimum cost there and decreases the upper count.

For completeness, the usual fundamental-circuit criterion says a basis
is minimum precisely when an omitted vector is at least as expensive as
every retained vector with a nonzero fundamental coefficient. For a
basis minimum at \(k\), a violation at \(k+1/2\) can only have an
omitted lower and retained upper. Their integral label difference must
equal \(k\). Consequently the stopping test gives the unique greedy
basis throughout \((k,k+1)\). This proves correctness even when many
exchanges occur at one crossing, without assuming generic coefficients.

Each basis replacement updates all pure and head coordinates by the
ordinary rank-one pivot formula in \(O(A^2)\) field operations.
There are \(O(A)\) replacements and \(O(A^3)\) initial field work.
Candidate searches inspect \(O(A^2)\) label/coefficient pairs per
replacement or interval. They require \(O(A^3)\) scalar comparisons
and field-zero tests over the complete compilation.

At each of the \(O(A)\) price intervals, form the coefficient rows on
the heads and recompute their decreasing-price, minimum-row flag pivots.
Eliminating \(O(A)\) rows of dimension at most \(t\) costs
\(O(At^2)\) field operations. This adds \(O(A^2t^2)\) in total.
The resulting pure bases, head birth labels and merged interval records
are exactly those of the earlier compiler. Field storage is \(O(A^2)\),
and the retained integer records use \(O(A^2)\) space. For
\(t^2\le A\), the field bound becomes \(O(A^3)\).

## Merge two ordered lists instead of inspecting every pair

At \(w=(x,y,z)>0\), set
\[
 \theta=(x-tz)/y,\quad\phi=(x-ty)/z,\quad
 a_i=y(\alpha_i(\theta)-i),\quad
 b_j=z(\beta_j(\phi)-j).
\]
Both lists are nonincreasing by (AL4). The entire correction in (AB7)
is \(\sum_{i,j}\max\{a_i,b_j,0\}\).
Within the two compiled records their integer gradients are
\[
 A_i=(\delta_i,q_i-i,-t\delta_i),\qquad
 B_j=(\epsilon_j,-t\epsilon_j,r_j-j).                         \tag{AL5}
\]
Discard nonpositive entries, and merge the two remaining ordered lists.
Keep the original head-index order when entries from the same arm are
equal. Reject a tie between positive entries from different arms.

When \(a_i\) is reached after \(s\) positive Z entries, it wins
exactly \(t-s\) pairs. When \(b_j\) is reached after \(r\)
positive Y entries, it wins exactly \(t-r\) pairs. Therefore the
gradient of the correction is obtained by adding
\[
                 (t-s)A_i\quad\hbox{or}\quad(t-r)B_j          \tag{AL6}
\]
at each merge step. Nonpositive and absent entries of the opposite arm
are already included in these counts. Pairs with neither entry positive
contribute zero. This counts every pair exactly once. Subtract (AL6)
from the two pure profiles and add the fixed profile (AC6).

The two record searches cost \(O(\log A)\). Reading the gains,
checking their order, merging and adding the three-coordinate integer
vectors cost \(O(t)\). No sorting or field elimination is performed
at query time. Formula (AB7), rather than a separate rank approximation,
proves that the result is the original minimum profile.

## A linear number of inequalities determines a whole patch

Keep at most four record-endpoint inequalities. For every finite gain,
keep its sign inequality; this contributes at most \(2t\) normals.
Finally keep the comparison between each consecutive pair in the merged
positive list that comes from different arms. There are at most
\(2t-1\) such comparisons.

On the intersection with the two factor records, the within-arm ordering
holds everywhere by (AL4). The consecutive cross-arm comparisons
therefore preserve the entire merged order by transitivity. Together
with the sign inequalities they preserve every pair winner in (AB7).
The returned profile is consequently valid on the whole polygon, not
only at its query. Conversely this sign pattern and merged word have
exactly these inequalities, so their signature identifies one convex
open patch. Queries on a sign boundary may be rejected even when the
full minimum is differentiable there; this is an additional refinement.
Every finite gain-zero is an integral factor-price line, as proved in
the birth-curve note. Thus the earlier complete arrangement still refines
these patches, and
\[
                         K_{\rm lin}=O(A^2+At^4).             \tag{AL7}
\]

All comparison normals have height \(O(a)\). The existing exact
adjacent-weight construction and signature walk apply. A full traversal
uses \(O(K_{\rm lin})\) patch calls, with \(O(t+\log A)\) scalar
work in each decision. Rational polygon clipping and profile merging
are additional costs. Merge the refinement patches by their original
profiles as in the compilation note. The resulting full normal cones
and initial ideals are unchanged. This does not bound the refined count
by the actual-count exponent \(17/8\).

## Verification and limits

~~~bash
python3 -B scripts/run_hook_axis_linear_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker verifies (AL4) throughout complete compiled price intervals,
including both tails, over rational and finite fields with nonunit and
zero controls. Direct full-prefix ranks independently check selected
birth thresholds. Every complete incremental factor table is compared
with fresh elimination, including arbitrary replicated factors with
nonidentity upper bases and dependent lower vectors. Every new full fan
is compared by its exact profile
polygons with the earlier compiled report; selected entire original
quotient fans are also replayed. Each linear patch is checked against
all inequalities from the quadratic oracle. An explicit arbitrary
replicated factor guards the specialization in (AL2).

Ordered merging and the greedy/dual framework are classical algorithmic
inputs. Parametric arrangement traversal and basis exchange updates
already appear in [Hausbrandt--Ruzika, Section 3 and Theorem 3.6](https://arxiv.org/html/2503.08178v2).
The state-polytope dictionary and a broader universal-basis algorithm
are given by [Babson--Onn--Thomas, Theorems 2.4 and 4.2](https://arxiv.org/pdf/math/0207135).
The elementary shift argument above identifies why ordered merging
applies to the actual axis hooks, and the two price classes give the
displayed incremental compilation bound. Strongest-current priority of
this combined specialization is unestablished. It improves the compiler
and profile oracle, without a new asymptotic fan count, a polynomial
lower excess, a smooth-number unranking speed claim or completion of
the original extremal objective.
