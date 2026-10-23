# Rational leading coordinates bound fans without expanding long inverses

Long forward and inverse propagation do not force a large hook fan.
A compact rational representation supplies a price-dependent elimination
and an actual coordinate separating every possible middle-arm exchange.
For two arms whose numerator and denominator have a pure long leading
term on their heavy-coordinate sector, the entire compatible fan is
quadratically bounded when their degrees are \(O(t)\). Both outer
images are arbitrary, over every field.
The numerator may be a nonunit. A filtered relation between outer
sources removes any additional axis-order or image-membership assumption.

The geometric-series arm \(1+v/(1-q)\) is included with rational degree
one, although both its expanded operator and inverse have degree
\(A-1\). This closes a long-propagation family left outside the previous
unit-window closure. The unrestricted sharp count, an actual polynomial
fan excess, exact current novelty and the original paper criteria remain open.

## Compact rational hypotheses and the full bound

Use the [original compatible hook](hook_bridge_reduction.md)
\[
 a\ge2t,\quad t\ge2,\quad A=a-t,\quad n=4at-2t^2.
\]
For each arm use local long and short variables \(q,v\), and suppose
its central operator is
\[
 W=P/Q\quad\text{in }K[q,v]/(q^A,v^t).                     \tag{RL1}
\]
The division is truncated formal division. Require \(Q(0,0)\ne0\),
allow \(P(0,0)=0\), and take integers \(p,r<A\) such that
\[
 \begin{array}{c}
 \operatorname{supp}P\subseteq\{(i,j):i+j\le p\},\quad
                         [q^p]P\ne0,\\
 \operatorname{supp}Q\subseteq\{(i,j):i+j\le r\},\quad
                         [q^r]Q\ne0.
 \end{array}                                               \tag{RL2}
\]
The two arms may use different polynomials and degrees. Put
\(d=\max(p_Y,r_Y,p_Z,r_Z)\). Every compatible coefficient of both
outer images is arbitrary. Then, over every field,
\[
                         g(I_\phi)=O(a^2(t+d)^2)
                                   =O(n^2(1+d/t)^2).       \tag{RL3}
\]
For fixed \(a,t,d_0\), one coefficient-independent complete positive
integer weight list for the entire class with \(d\le d_0\) has
\[
                O(a^2(t+d_0)^2)\text{ weights of height }O(a^2).
                                                               \tag{RL4}
\]
No coprimality, convergence or characteristic-zero assumption is made.
The leading and constant coefficients must remain nonzero in the field.
Earlier bounds can be retained when smaller.

## One arm on its heavy-coordinate sector

Choose Y and normalize the objective to \((\lambda,u,1)\), with
\(u>1\). Under (RL2), \(q^p\) and \(q^r\) are the unique
highest-weight monomials of \(P,Q\): any other monomial has
\(iu+j<pu\), or \(iu+j<ru\), respectively. Their pure long
leading terms survive every short-coordinate truncation.

Let \(U_{q,j}\) be the original upper Y column of label
\((1,t+q,j)\), and \(L_{q,j}\) the central source of label
\((0,t+q,t+j)\). The latter is pure when \(q\ge t\).
The actual relation \(QW=P\), shifted by \(q^kv^j\), is
\[
       \sum_{i,h}Q_{ih}L_{k+i,j+h}
                  =\sum_{i,h}P_{ih}U_{k+i,j+h}.           \tag{RL5}
\]
Terms outside the arm rectangle vanish; every nonzero lower term is
pure if \(k\ge t\). All original column labels and costs are retained.
Define the price plane
\[
                         \lambda_0=t+(r-p)u.             \tag{RL6}
\]

If \(\lambda>\lambda_0\), an upper input \(q\ge t+p\)
is a loop modulo strictly cheaper columns. Solve (RL5) for the
\(P\)-leading upper term, taking \(k=q-p\). Every other upper
term is strictly cheaper, and the lower terms cost at most
\(t-\lambda+(r-p)u<0\) more than that upper term.
If \(\lambda<\lambda_0\), a lower central input
\(q\ge t+r\) is likewise a loop, by solving for the leading
\(Q\)-term. The other lower terms are cheaper and all upper
terms have price difference at most
\(\lambda-t+(p-r)u<0\).

Repeated elimination decreases the same-colour long input, and every
new column is strictly cheaper. Thus the deletions preserve the span
of every cost prefix as well as every minimum basis throughout each
open price region. The price plane itself is added to the wall inventory.

## A single actual coordinate excludes a middle exchange

Apply multiplication by the unit \(Q\) to the chosen Y row space,
leaving Z unchanged. Upper columns become shifts of \(Q\), and
central columns become shifts of \(P\). This is an invertible row
operation, not an algebra coordinate change. It preserves the far-Y
subspace whose original long labels are at least \(a-t\).

Suppose a coherent pure cross-arm exchange uses an upper Y input
\((q,j)\) with \(q\ge t+p\). Away from (RL6), the loop test
forces \(\lambda<\lambda_0\). Consider the transformed row
\((q+r,j)\). Its coefficient in this upper column is \([q^r]Q\ne0\).
Every other upper column containing that coordinate is more expensive:
its input weight is strictly larger than \(qu+j\). A central
column containing it has cost at least the upper column's cost plus
\(\lambda_0-\lambda>0\), since the largest weight in \(P\) is
\(pu\). These statements include central bridge columns.

If the row survives below the far band, and the common exchanged cost
\(\tau<au\), it appears in no other column of the entire
\(\tau\)-prefix. Indeed \(y^a\)-sources are too expensive;
\(z^a\)-bridge projections are in the far band by compatibility;
and opposite pure columns have zero Y component. The upper column is
therefore a coloop of that prefix and cannot exchange with an opposite
pure column.

For a lower central input \(q\ge t+r\), the surviving region is
\(\lambda>\lambda_0\). The same argument uses row
\((q+p,j)\), with nonzero coefficient \([q^p]P\). Every other
central column containing it is more expensive, and every upper column
containing it has cost at least \(\lambda-\lambda_0>0\) more.
It too is a prefix coloop under the same two geometric conditions.

A prefix coloop must belong to the basis of the prefix in every minimum
basis. Additional cost ties do not permit its removal by a coherent
exchange. Nothing is asserted about columns more expensive than \(\tau\).

## Endpoint strips on the heavy sector

Let \(M\) be the original long label of the chosen Y column.
The small inputs excluded from the loop tests have \(M<2t+d\).
For every other central or upper coherent exchange, the coordinate
argument shows that either
\[
                          M\ge a-t-d
                  \quad\text{or}\quad \tau\ge au.        \tag{RL7}
\]
An eligible long upper column has
\(\lambda\le\lambda_0\le t+du\), hence
\(\tau\le u(M+d)+2t-1\). A lower central column has
\(\tau\le uM+2t-1\). Since \(u>1\), (RL7) places both types
within distance at most \(d+2t\) of \(a\).

A pure outer Y source has label \((0,a+i,j)\), \(i\ge t\).
The filtered outer-source identity proved below makes it a loop modulo
strictly cheaper original columns whenever
\(i\ge\max(2t,t+p)\). This argument permits every compatible outer
image, including ones outside the multiplication image of \(W\).
Consequently every actual pure cross-arm exchange in this sector,
away from (RL6), satisfies
\[
                  M\le3t+d\quad\text{or}\quad
                                  |M-a|\le3t+d.          \tag{RL8}
\]
No regular-prefix hypothesis or small expanded inverse has been used.

## Glue the two heavy sectors

For \(w_y>w_z\), use Y. For \(w_z>w_y\), use Z with long-to-short
ratio \(w_z/w_y>1\). Add the separating plane \(w_y=w_z\).
For a pure pair write \(b=M-r_0\), \(c=K-s_0\), where the
original short labels \(r_0,s_0<2t\). Equation (RL8) puts every
actual cross-arm wall, except the two price planes, in one of four strips:
\[
 \begin{array}{c}
 |b|\le6t+2d,\quad |b-a|\le6t+2d,\\
 |c|\le6t+2d,\quad |c-a|\le6t+2d.
 \end{array}                                               \tag{RL9}
\]
These have \(O(a(t+d))\) original difference vectors. Comparisons
within an arm or involving a bridge add \(O(at)\) planes.
For a degree cap include every possible price plane with
\(r-p\in[-d_0,d_0]\), on both arms. There are only \(O(d_0)\)
of these. The common arrangement thus has \(O(a(t+d_0))\) planes
and \(O(a^2(t+d_0)^2)\) positive chambers, proving (RL3)--(RL4).

The classical minimum-basis face and exchange-edge facts ensure that
every positive wall is detected by a coherent exchange; see
[Feichtner--Sturmfels, Section 2](https://arxiv.org/pdf/math/0411260v1).
Original normals have height \(O(a)\). Primitive cross products and
sums of three independent extreme rays give the stated integer weights.
No optimal total construction or rank-evaluation cost is claimed.

## Two long-propagation consequences

On both arms let
\[
                W=1+\frac{v}{1-q}
                  =\frac{1-q+v}{1-q}.                     \tag{RL10}
\]
The numerator and denominator meet (RL2) with \(p=r=1\), over
every field. Both expanded \(W\) and \(W^{-1}\) have long degree
\(A-1\) when \(t\ge2\), but the full compatible fan is \(O(n^2)\)
with arbitrary outer maps. A common quadratic integer list works
uniformly over their coefficients. With zero outer maps, the previous
[interface analysis](hook_tail_interface.md) gives
\(b=t^2\), \(\chi=t-1\), \(\rho=(A-1)(t-1)\) and maximal geometric
prefix gaps. The compact rational argument improves that family's
named \(O(n^2t^3)\) ceiling to \(O(n^2)\). No sharpness or
strongest-current claim is inferred from that comparison.

There is also a one-arm consequence with an arbitrary opposite arm.
If \(P,Q\) depend only on \(q\), their pure leading monomials
are uniquely highest for every \(u>0\). The same proof then gives
endpoint or high-cost windows on the whole positive sector. Small inputs
already have long label below \(2t+d\); the remaining endpoint or
high-cost alternatives imply (OW5) with \(d\) in place of \(D\).
The eligible upper price is at most \(t+du\). Reusing
the [translated-hyperbola count](hook_outer_unit_windows.md), with
\(d\) in place of the expanded propagation bound, yields
\[
             g(I_\phi)=O\bigl(n^2(1+d/t+\log(a/t))^2\bigr)
                                                               \tag{RL11}
\]
with all compatible outer images and opposite coefficients arbitrary.
Both expanded degrees can be \(A-1\), for example
\(W=(1+q)/(1-q)\) in characteristic different from two.
The two-arm theorem (RL3) permits mixed short coefficients in both
rational representations and removes this logarithmic allowance.

The scalar one-arm corollary stated here retains its unit hypothesis:
its proof uses the all-ratio price bound on an individual outer column.
The new outer relation below uses \(u>1\) and does not by itself
remove that hypothesis from (RL11).

## A filtered outer relation removes the nonunit loss

Write \(F=f_{0,Y}/y^t\in R\) for the arbitrary matching outer
component. Let \(O_{i,j}\) denote its original pure outer column,
with label \((0,a+i,j)\), when \(i\ge t\).
Since \(QW=P\), multiplication in the truncated rectangle gives
\[
 \sum_{h,l}P_{hl}O_{k+h,j+l}
       =\sum_{r,s}(QF)_{rs}L_{k+r,j+s},\qquad k\ge t.    \tag{RL12}
\]
Terms outside the rectangle are zero. Every nonzero column in this
identity is an original pure column: both long inputs are at least
\(k\ge t\). In particular arbitrary cross components of the outer
images vanish in these shifts. The Z identity is its transpose.

Take an outer input
\[
                         i\ge\max(2t,t+p),              \tag{RL13}
\]
and set \(k=i-p\) in (RL12). The coefficient of \(O_{i,j}\)
is the nonzero coefficient \([q^p]P\). Every other outer term is
strictly cheaper on \(u>1\), because \(q^p\) is the unique
highest-weight monomial of \(P\). Its long input is also strictly
smaller than \(i\). Every nonzero central term costs at most
\((a-1)u+2t-1\), whereas the target costs \((a+i)u+j\).
The first bound is strictly smaller when \(i\ge2t\).
Solving (RL12) therefore expresses the target in strictly cheaper
original columns, for every \(\lambda>0\). Recursing only on the
outer terms terminates by decreasing long input. All cost-prefix
spans and all minimum bases are preserved throughout the heavy sector.

No division by \(P\), inverse of \(W\), or individual membership
\(O_{i,j}\in WR\) is used. In fact early outer columns can lie
outside the central span. The relation nonetheless proves the endpoint
strip (RL8), and hence the full degree-only bounds (RL3)--(RL4),
for arbitrary nonunit numerators as stated at the start.

For example, both central arms may be
\[
                       W=q^\nu+v,\qquad1\le\nu<A.        \tag{RL14}
\]
Their arbitrary-outer full fan satisfies
\(g=O(n^2(1+\nu/t)^2)\). Thus even a growing axis order
\(\nu=O(t)\) gives a quadratic fan. The earlier axis-power argument
gave only \(O(n^2(1+\nu)^2)\) for these arbitrary outer images.
Its membership certificate remains valid, but is unnecessary for the
stronger filtered bound. This closes that growing-axis-order candidate;
it supplies neither an unrestricted bound nor a polynomial lower excess.

## Verification and scope

~~~bash
python3 -B scripts/run_hook_rational_leading_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checks reconstruct actual quotients and commuting multiplication,
verify denominator clearing against every original central column,
certify original-label loop identities and whole open-region cone
patches, and inspect complete cost prefixes for the separating row.
Complete original fans test both heavy sectors and their price-plane
exceptions. A shared integer list is checked across different fields
and rational coefficients. Long geometric propagation is measured
separately from compact rational degree.
Nonunit controls verify matching outer-image factorizations against
every original shifted pure outer column, including nonzero outer images.
Additional controls use arbitrary nonunit outer images outside that
matching-image class. They retain the earlier axis-power certificate
and central-span reconstructions as controls, and independently verify
(RL12) with the actual outer coefficients. The filtered relations and
strict original prices are checked even before the old axis-power
membership threshold. Rank reduction against all cheaper columns and
complete original fans provide separate checks. Growing-axis examples
use \(\nu\) proportional to \(t\). Early outer columns are also
checked to lie outside the central span.

Rational formal division, triangular leading-term elimination, matroid
faces and plane arrangements are classical. The quantitative coherent
application is proved here; exact strongest-current novelty and
sharpness remain unestablished. General operators need not have small
rational representations satisfying (RL2). Large degree relative to
thickness, failure of a pure long leading representation, and the
unrestricted sharp extremal problem remain open. In particular the
mixed-leading units \(1+q^\nu v\) do not satisfy these hypotheses
through their displayed numerator.
