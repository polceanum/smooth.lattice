# One cross-arm exchange can disappear and return

The remaining cross-arm walls in the [full hook reduction](hook_bridge_reduction.md)
cannot be counted by giving each literal pair of monomials one connected
interval of eligibility. An actual length-102 hook quotient has a single
pure cross-arm pair active on three disjoint intervals. These intervals
contain four different wall pieces. The direction itself has seven pieces.

This is a finite obstruction to a proposed counting argument, together
with an exact direction-query algorithm. It is not a growing lower family,
an improved full-hook or unrestricted exponent, or a novelty claim. The
original extremal objective and the remaining collective estimate are
unchanged.

## Query a direction using actual field ranks

Let a finite represented matroid have distinct integer labels
\(E\subset\{0,1\}\times\mathbb Z^2\), and let
\(q=(1,q_y,q_z)\). Parametrize the positive objectives on its plane by
\[
             w(u)=(-q_yu-q_z,u,1),\qquad w(u)>0.                 \tag{R1}
\]
Their domain is an open interval, possibly unbounded. Put
\(d=(-q_y,1,0)\). At the left endpoint \(u_0\), run actual matroid
greedy twice with lexicographic element costs
\[
       (w(u_0)\cdot a,\ d\cdot a,\ +a_x),\qquad
       (w(u_0)\cdot a,\ d\cdot a,\ -a_x).                     \tag{R2}
\]
The derivative selects the state immediately to the right; the last
entry selects the two extrema of the first moment among permanent ties.
A deterministic element order resolves any remaining ties.

For either selected basis, its fundamental circuits give all inequalities
\(w\cdot(a-b)\ge0\) for minimum weight. Substitution of (R1) gives
linear rational inequalities in \(u\). Intersect them for the two bases
and the positive domain. The right endpoint is the next exact event.
Repeat there using (R2). On the resulting open interval, let the basis
profiles be \(\mu_-\) and \(\mu_+\), ordered by first coordinate. Then
\[
                  \mu_+-\mu_-=h q,\qquad h\ge0.               \tag{R3}
\]
Indeed both bases have the same affine cost throughout the interval, so
their difference annihilates both vectors spanning \(q^\perp\).
It has the displayed direction, and the first-moment choices give its
sign. If \(h>0\), the minimizing profile face is a nondegenerate segment
parallel to \(q\). If \(h=0\), it is a point. Thus the positive intervals
with \(h>0\) count exactly the parallel edge pieces, retaining subdivisions
at intervening fan vertices.

The procedure terminates: all events belong to the finite arrangement
of affine element-cost comparisons on (R1). It uses two greedy and
fundamental-circuit computations per shadow interval. It also visits
intervals where \(h=0\); no complexity bound depending only on the number
of parallel edges is asserted. This is classical parametric greedy
optimization with an explicit degeneracy convention, not a new
parametric-matroid complexity theorem. No support matching replaces field
rank.

At an interior parameter avoiding extra projected-label ties, a literal
pair \(e,f=e+q\) is exchanged exactly when it is nonloop parallel after
contracting the columns strictly cheaper than its common cost. If \(C\)
is that set, the condition is
\[
 r(Ce)=r(Cf)=r(Cef)=r(C)+1.                                  \tag{R4}
\]
Every grading fiber has at most two labels because \(a_x\in\{0,1\}\).
This also identifies the literal exchanges from the two extreme bases.
Along (R1), the cheaper set need not be nested. Consequently (R4) can
hold, fail, and hold again. The connected-interval proof for a
[same-colour tie under a scalar colour shift](two_colour_threshold_braids.md)
uses nested cheaper prefixes; that hypothesis is absent here.

## The exact recurrent pair

Use the compatible hook family
\[
 M=(y^{10},y^3z^3,z^{10}),\quad C=K[y,z]/M,\quad
 A_\phi=K[x,y,z]/(x^2,\ g-x\phi(g):g\in M).
\]
It has length \(102\). The notation means the usual compatible flat
extension; its three generator images and their allowed independent
coefficients are those of the [generic-hook construction](generic_hook_fans.md).
The full coefficient list of a common integer specialization is recorded
in the accompanying report. Take
\[
 q=(1,-6,4),\quad w(u)=(6u-4,u,1),\quad u>2/3,
 \quad e=y^7z^3,\quad f=xy z^7.                              \tag{R5}
\]
The pair has difference \(q\), and this direction is outside the exact
same-arm/bridge inventory of the hook reduction. Its codimension-one wall support has the following three open intervals,
with subdivision vertices restored inside each interval:
\[
 (5/6,6/7),\qquad(8/7,6/5),\qquad(9/7,4/3).                  \tag{R6}
\]
The middle component is split at \(7/6\), giving four pieces for this
one pair. The complete \(q\)-wall inventory has seven pieces:

| Open interval in \(u\) | Lower monomial | Upper monomial |
|---|---|---|
| \((5/6,6/7)\) | \(y^7z^3\) | \(xy z^7\) |
| \((8/9,1)\) | \(y^6z^4\) | \(xz^8\) |
| \((8/7,7/6)\) | \(y^7z^3\) | \(xy z^7\) |
| \((7/6,6/5)\) | \(y^7z^3\) | \(xy z^7\) |
| \((9/7,4/3)\) | \(y^7z^3\) | \(xy z^7\) |
| \((5/3,7/4)\) | \(y^6z^4\) | \(xz^8\) |
| \((7/4,9/5)\) | \(y^6z^4\) | \(xz^8\) |

The query traverses 33 shadow intervals in total. These statements hold
for the recorded specialization over \(\mathbb F_{1000003}\) and over
\(\mathbb Q\). They also hold at generic coefficients in every
characteristic, over an infinite field or a sufficiently large finite
extension. No assertion for every specialization is intended.

## Certificate and implications

```bash
python3 -B scripts/run_hook_wall_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker reconstructs the original quotient columns and commuting
multiplication, performs (R1)--(R3) over the actual prime field, and checks
(R4) separately for every active pair. A separate lower-envelope
implementation uses all profiles in the previously certified complete
89-cone and 948-cone generic fans; it agrees on every open interval,
including those without a parallel edge. Every selected support cone is
also checked on its entire claimed interval by its fundamental inequalities.
Small arbitrary represented
matroids over \(\mathbb Q,\mathbb F_2,\mathbb F_3,\mathbb F_5\) are
also compared against independently enumerated lists of all bases.

For the rational specialization, every selected actual basis has nonzero
integer minor modulo \(1000003\), hence nonzero integer determinant.
It is also a minimum support-matching basis. Every actual independent
basis is support-independent, so the support lower bound is attained at
all query intervals by that same rational quotient. This proves the
complete rational query without assuming that every support-matching
basis survives.

For the all-characteristic generic statement, the checker independently
replays the saved integer determinant and content-prime exclusion
certificates for every coherent basis used in the query. Each relevant
minor polynomial is nonzero modulo every prime. Their product is nonzero
over the corresponding polynomial ring, giving one common nonempty open
set where the whole query is attained. This uses the existing finite
certificates and supplies no general convolution-rank theorem.

The example rules out connectedness for a fixed literal cross-arm pair,
not an \(O(t^2\operatorname{polylog}a)\) bound on its direction's total
multiplicity. It also does not show that multiplicity grows at fixed
\(t\). The next bound must control repeated eligibility and subdivisions
collectively, or use a different invariant. The retained full-hook box
bound, and the open unrestricted extremal target, are unchanged.
