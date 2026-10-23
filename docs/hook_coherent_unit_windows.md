# Coherent unit windows give quadratic fans despite mixed support components

Connected support is larger than coherent interaction. A chosen unit arm
can be normalized in two ways, retaining the original monomial labels.
Contracting strictly cheaper identity columns in the appropriate chart
confines every actual cross-arm exchange to a short strip. This closes
the mixed-direction family from the monomial-component analysis, even
with arbitrary compatible coefficients on the opposite arm.

The bound is uniform over all fields. It gives a common quadratic list
of integer weights for every fixed bound on forward and inverse
propagation. The unrestricted sharp extremal count, an actual polynomial
lower separation, strongest-current novelty and the original paper
requirements remain open.

## The unit propagation theorem

Use the [original flat hook model](hook_bridge_reduction.md)
\[
 a\ge2t,\quad t\ge2,\quad A=a-t,\quad n=4at-2t^2,
\]
with residual row spaces \(T_Y\oplus T_Z\) and original groups
\(E_Y,E_Z,H\). Assume the Y components of both outer images vanish.
The central Y image is \(y^tW(q,v)\), where
\[
 R=K[q,v]/(q^A,v^t),\qquad W\in R^*.
\]
Here \(q\) is the long coordinate measured from \(y^t\), and
\(v\) is the short coordinate. Every Z component of all three images,
including the outer images, is an arbitrary compatible coefficient.
Let
\[
 d_+=\deg_q W,\quad d_-=\deg_q W^{-1},\quad D=\max(d_+,d_-).
                                                               \tag{CW1}
\]
The inverse and degrees are taken in the displayed truncated ring.
For every such actual quotient,
\[
             g(I_\phi)=O\bigl(a^2(t+D)^2\bigr)
                      =O\bigl(n^2(1+D/t)^2\bigr).       \tag{CW2}
\]
The constants are absolute. Transposition permits either arm to be
chosen, when that arm meets the stated hypotheses.

For fixed \(a,t,D_0\), one coefficient-independent complete positive
integer weight list works for every field and all compatible maps in
this class with \(D\le D_0\):
\[
 |\mathcal W|=O\bigl(a^2(t+D_0)^2\bigr),\qquad
                      \|w\|_\infty=O(a^2).             \tag{CW3}
\]
No optimal construction time or total cost of the quotient-rank
computations is asserted.

## Two charts preserve the actual exchange question

Write a positive objective as \(w=(\lambda,u,1)\).
An upper Y column with input \((q,j)\) has original label
\((1,t+q,j)\), and a pure lower central column has label
\((0,t+q,t+j)\). The latter is pure exactly when \(q\ge t\).
Upper inputs with \(q<t\) already have long label less than \(2t\).
All other nonzero pure lower Y columns are absent, by the chosen-arm
outer-image hypothesis.

Call \(e\leftrightarrow f\) a coherent exchange if two minimum
bases at a positive objective differ by precisely these columns.
The argument includes objectives with additional ties. Put
\[
 H_\pm(u)=(t-1+d_\pm)u+t-1.                              \tag{CW4}
\]
For a pure cross-arm coherent exchange whose Y input is \((q,j)\)
with \(q\ge t\),
\[
 \begin{array}{ll}
 qu+j\le H_+(u),&e\text{ is upper},\\
 qu+j\le H_-(u),&e\text{ is lower}.
 \end{array}                                                   \tag{CW5}
\]
These are necessary eligibility bounds, not bounds for every vector.

For an upper column, use the original identity-row chart. Every Y
entry of a bridge comes from \(W\) times an input with long coordinate
less than \(t\). Its input-row cost \(q'u+j'\) is at most
\(H_+(u)\). If \(qu+j>H_+(u)\), all upper identities on these
rows are strictly cheaper than \(e\). Contract their span. Every
bridge now has zero Y projection, so the whole represented matroid
splits into its Y part and its Z part. It cannot exchange \(e\)
with an opposite pure column.

For a lower column, apply \(W^{-1}\) to the Y row space and leave
Z unchanged. This is an invertible row operation, with every original
column label and cost unchanged. Pure lower inputs with \(q'\ge t\)
are now identity columns. Central bridges have Y support only at input
rows with \(q'<t\); the two outer Y bridge projections are zero.
Upper columns with inputs \(q'\ge t\) stay in the tail \(q'\ge t\).
Upper columns with inputs \(q'<t\) have every output cost at most
\(H_-(u)\).

If \(qu+j>H_-(u)\), contract all tail identity columns whose
input-row cost is less than \(qu+j\). Their original lower-column
costs are strictly less than the exchanged lower column's cost.
After this contraction, every head-input upper column has only head
Y support: all its tail outputs have been killed. The entire represented
matroid splits into the surviving Y tail and the head-Y-plus-Z part.
Again the proposed cross-arm exchange is impossible.

For completeness, contraction by a strictly cheaper independent set is
valid even if that particular set is not contained in a minimum basis.
Every minimum basis contains a basis of the full strictly cheaper span.
Contract that span, which contains the chosen set. The direct-sum split
survives this further contraction, and any remaining basis must have the
fixed rank of each summand. A single exchange between summands is
impossible. No contraction here is asserted to be an algebra quotient.

## A price bound turns eligibility into short strips

An upper Y column with \(q\ge t\) is a linear combination of pure
lower central columns, using \(W^{-1}\) at that input. Each such
source has input \((q+r,j+s)\), with
\(0\le r\le d_-\), \(0\le s<t\). Its original cost minus
the upper column's cost is
\[
                         t-\lambda+ru+s.
\]
Consequently the upper column is a loop modulo strictly cheaper
columns if \(\lambda>2t-1+d_-u\). Every eligible upper Y column
therefore satisfies
\[
                         \lambda\le2t-1+d_-u.           \tag{CW6}
\]

Let \(M\) be the Y long label, \(K\) the opposite Z long label,
and let their short labels be \(s,r<2t\), respectively. Their equal
costs give
\[
 K=(e_x-f_x)\lambda+u(M-r)+s.                            \tag{CW7}
\]
From (CW5), for either type of Y column,
\[
                        uM\le(2t-1+D)u+t-1.             \tag{CW8}
\]
If the Y column is lower, \(e_x-f_x\le0\), so (CW7) and
(CW8) bound \(K\) directly. If it is upper, use (CW6) when
\(f_x=0\); for \(f_x=1\) the price cancels. In all cases,
\[
                         K\le(2t+2D)u+4t.              \tag{CW9}
\]
Thus every coherent pure cross-arm exchange satisfies
\[
                  \boxed{\ M\le3t+D\quad\text{or}\quad
                                  K\le6t+2D.\ }         \tag{CW10}
\]
Indeed for \(u\ge1\), (CW8) bounds \(M\); for \(u\le1\),
(CW9) bounds \(K\). The previously omitted upper inputs \(q<t\)
already satisfy the first alternative. Same-layer cross-arm exchanges
are included; no assumption of a unique tied pair was used.

## Full fan and common integer compilation

All comparisons inside one arm or involving a bridge have \(O(at)\)
distinct planes, by the original arm and three translated bridge boxes.
For a cross-arm exchange, (CW10) confines one long label to
\(O(t+D)\). The other transverse difference has \(O(a)\) values,
and the first difference belongs to \(\{-1,0,1\}\). Thus a common
inventory of \(O(a(t+D))\) original comparison planes contains every
actual positive wall.

To see why it suffices to consider exchanges, a minimum-basis face of
a represented matroid is itself the basis face of its greedy minors.
Any two minimum bases are connected by single minimum-basis exchanges.
At an interior point of a codimension-one wall avoiding other distinct
comparison planes, a nontrivial such exchange lies on that wall's
plane. If it is not an arm or bridge comparison, it obeys (CW10).
This is the usual greedy and basis-exchange argument; no determinant
is replaced by support matching.

An arrangement of \(L\) planes through the origin has \(O(L^2)\)
positive three-dimensional chambers. It refines the entire positive
initial fan for every coefficient choice in the stated class. This
proves (CW2), and gives the common list in (CW3). Original label
differences have entries less than \(2a\). Add coordinate boundary
planes, take primitive integer extreme rays as cross products, and
sum three independent rays of each full chamber. The sum lies strictly
inside that chamber and has height \(O(a^2)\).

## Arbitrary mixed support in a bounded cone

Suppose the chosen unit arm has the form
\[
 W=\gamma+\sum_{j\ge1}\sum_{0\le i\le Cj} c_{ij}q^iv^j,
                 \qquad\gamma\ne0,\quad C\ge0,          \tag{CW11}
\]
in the truncated ring. All displayed coefficients are arbitrary.
The support cone is closed under multiplication. The finite geometric
series for \(W^{-1}\) therefore stays in the same cone, so
\[
                 D\le C(t-1),\qquad
                 g(I_\phi)=O\bigl(n^2(1+C)^2\bigr).      \tag{CW12}
\]
For every fixed \(C\), the whole class has a quadratic fan and a
common quadratic integer weight list, with arbitrary compatible opposite
arm and outer coefficients. There may be several independent monomial
steps and a marked support component of length comparable to \(A\).

In particular
\[
                  W=(1+qv)^{t-1}+\lambda q^2v            \tag{CW13}
\]
has \(C=2\), for every field and every parameter, and hence a full
\(O(n^2)\) fan bound. In characteristic zero or \(p>2t\), when
both arms have this operator, \(\lambda\ne0\), and \(A\ge2t+2\),
the [component analysis](hook_monomial_orbits.md) proves simultaneously
\[
 b=t^2,\quad \chi_Y,\chi_Z\ge\lfloor t^2/4\rfloor,
 \quad\rho=(A-1)(t-1),\quad\sigma_Y,\sigma_Z\ge A-t-1.
\]
The present proof closes that large-rank, long-component family by
coherent eligibility, rather than inferring its fan count from those
parameters. It also permits an arbitrary opposite arm.

Conversely, in this chosen-arm setting, an actual family with
\(g\ge c n^{2+\epsilon}\) must have
\(D/t=\Omega(n^{\epsilon/2})\). Large forward or inverse long
propagation, or a failure of the chosen-arm outer-image hypothesis, is
therefore necessary. This is not a lower-bound construction.

## Verification and priority boundary

~~~bash
python3 -B scripts/run_hook_coherent_unit_window_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker reconstructs the original compatible quotient and commuting
multiplication. It separately verifies the polynomial inverse against
actual central columns, original-label source identities, strictly cheaper
contractions and the global direct-sum splits. Complete original fans
check the necessary bounds on actual minimum-basis exchanges, including
degeneracies. A common arrangement is reused across different coefficient
choices and fields. These finite checks support the uniform proof; they
do not replace it or establish the unrestricted endpoint.

Greedy minors, direct sums, unit inversion and plane arrangements are
classical. The exchange-edge and face facts used above are reviewed in
[Feichtner--Sturmfels, Section 2](https://arxiv.org/pdf/math/0411260v1).
The coherent-window application and its quantitative bound
are proved here. Exact strongest-current novelty and sharpness remain
unestablished. The unrestricted sharp count, an actual polynomial excess
and the original final-paper requirements are still open.
