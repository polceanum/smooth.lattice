# A filtered rational relation controls arbitrary terminal hook coefficients

The rational relation in the [leading-window theorem](hook_rational_leading_windows.md)
need only hold below a terminal band. The omitted coefficients can be
arbitrary, with mixed short-variable dependence, over every field.
This closes the separated-support axis calibration with a compact low
part and arbitrary coefficients in its last four long degrees.

Use the compatible hook with \(A=a-t\ge t\ge2\), length
\(n=4at-2t^2\), and local rectangular rings
\(R=K[q,v]/(q^A,v^t)\). Both outer images are arbitrary compatible
maps. On each arm assume there are \(P,Q\in R\), an integer
\(0\le h\le A\), and degrees \(p,r<A\) such that
\[
 Q(0,0)\ne0,\qquad QW-P\in q^{A-h}R,                         \tag{FR1}
\]
\[
 \operatorname{supp}P\subseteq\{(i,j):i+j\le p\},\quad
 [q^p]P\ne0,\qquad
 \operatorname{supp}Q\subseteq\{(i,j):i+j\le r\},\quad
 [q^r]Q\ne0.                                                 \tag{FR2}
\]
Different arms may have different data. Put \(d=\max(p_Y,r_Y,p_Z,r_Z)\)
and \(H=\max(h_Y,h_Z)\). Then
\[
 \boxed{\quad g(I)=O\bigl(a^2(t+d+H)^2\bigr)
                   =O\bigl(n^2(1+(d+H)/t)^2\bigr).\quad}       \tag{FR3}
\]
For fixed \(a,t,d_0,H_0\), the whole coefficient class has one
complete positive integer weight list of size
\(O(a^2(t+d_0+H_0)^2)\), with weight height \(O(a^2)\).
No optimal construction or rank-evaluation complexity is asserted.
The numerator may be nonunit. Formal truncated division requires no
analytic convergence, coprimality or characteristic-zero assumption.

In particular \(d,H=O(t)\) gives a quadratic full-fan bound with
arbitrary outer images. This is a family theorem; sharp unrestricted
\(G_3\), an infinite polynomial lower excess, strongest-current priority
and the original final-paper goal remain unresolved.

## Exact shifted relations survive the terminal error

Fix the Y arm and write \(E=QW-P\). Put \(m=\max(t,h)\).
Every shift \(q^kE\) vanishes for \(k\ge h\). Consequently, for
the original pure columns of labels
\[
 U_{q,j}=(1,t+q,j),\qquad L_{q,j}=(0,t+q,t+j),
\]
the exact original-column identity is
\[
 \sum Q_{ih}L_{k+i,j+h}=\sum P_{ih}U_{k+i,j+h},\qquad k\ge m.
                                                                    \tag{FR4}
\]
Terms truncated out of the rectangle vanish. Since \(k\ge t\),
every nonzero central source is pure and its opposite-arm projection
vanishes. Unshifted head inputs can retain a nonzero denominator error.

Normalize a heavy-sector weight as \((\lambda,u,1)\), \(u>1\),
and let \(\lambda_0=t+(r-p)u\). The pure long leading terms in
(FR2) uniquely maximize the weights of \(P,Q\). Solving (FR4) for
that leading term shows:

- if \(\lambda>\lambda_0\), every upper input \(q\ge m+p\)
  is spanned by strictly cheaper original columns;
- if \(\lambda<\lambda_0\), every central lower input
  \(q\ge m+r\) is spanned by strictly cheaper original columns.

The same-colour substitutions strictly lower long input and terminate.
Thus the span of every complete cheaper prefix is preserved. This does
not assume any generic Toeplitz minors or separate marginal-rank test.

## A middle separating row sees neither terminal nor outer errors

Apply multiplication by the unit \(Q\) to the chosen Y row space.
Upper columns become shifts of \(Q\), and central columns become
shifts of \(P+E\). The transformed error is supported in long rows
at least \(A-h\). Compatible opposite outer-bridge projections lie
in rows at least \(A-t\). Below row \(A-m\), neither occurs.

Suppose a coherent pure cross-arm exchange uses an upper input
\((q,j)\) with \(q\ge m+p\). The loop test forces
\(\lambda<\lambda_0\), away from the price plane. If
\(q+r<A-m\), the transformed row \((q+r,j)\) has a nonzero
coefficient in that upper column. Every other upper containing it is
more expensive, by the leading property of \(Q\). Every central
column containing it is more expensive by at least
\(\lambda_0-\lambda>0\), by the leading property of \(P\).
This includes bridge columns; the terminal part cannot reach the row.

If the common cost is \(\tau<au\), matching outer sources are
also too expensive, opposite pure columns have no Y projection, and
opposite outer bridges miss the row. The target is a coloop of the
entire \(\tau\)-prefix. It cannot take part in the proposed
cross-arm exchange, including an exchange on a face with other ties.

For a central lower input \(q\ge m+r\), use row \((q+p,j)\)
on the surviving side \(\lambda>\lambda_0\). The same argument
applies if \(q+p<A-m\) and \(\tau<au\). Thus every remaining
middle exchange has a row in the terminal band, or has \(\tau\ge au\).

Writing the original long label as \(M=t+q\), a terminal row gives
\(M\ge a-m-d\). For a high-cost upper exchange, the loop test
also gives \(\lambda\le t+du\), whence
\(\tau\le u(M+d)+2t-1\); therefore \(M>a-d-2t\).
For a high-cost lower exchange, \(\tau<Mu+2t\) similarly gives
\(M>a-2t\). Small inputs already have \(M<t+m+d\).

## Arbitrary matching outer images obey the same filtered identity

Let \(F\in R\) be the matching outer component, and let
\(O_{i,j}\) have original label \((0,a+i,j)\). Since
\(q^kFE=0\) for \(k\ge h\), the actual relation is
\[
 \sum P_{ih}O_{k+i,j+h}
       =\sum (QF)_{rs}L_{k+r,j+s},\qquad k\ge m.                \tag{FR5}
\]
It allows arbitrary \(F\), including components outside the image
of multiplication by \(W\). Other outer components disappear in
these pure shifts by compatibility.

For \(i\ge\max(2t,m+p)\), put \(k=i-p\) and solve (FR5)
for \(O_{i,j}\). Every other outer term is strictly cheaper and
has smaller long input. Every nonzero central term costs at most
\((a-1)u+2t-1\), less than the target \((a+i)u+j\) when \(u>1\).
Thus late outer inputs are loops modulo the entire cheaper span.
Eligible outer columns satisfy \(|M-a|\le2t+h+d\).

Combining the three arguments, every pure cross-arm exchange on this
heavy sector, except its price plane, satisfies
\[
                 M\le3t+h+d\quad\hbox{or}\quad
                                |M-a|\le3t+h+d.              \tag{FR6}
\]
No exact identity at the early head inputs was assumed.

## A common arrangement counts the whole fan

Use Y when \(w_y>w_z\), and Z when \(w_z>w_y\), adding the
plane \(w_y=w_z\). Original short labels are less than \(2t\).
It follows from (FR6) that every pure cross-arm direction lies in
one of the four strips
\[
 |b|\le6t+2d+2H,\quad |b-a|\le6t+2d+2H,\quad
 |c|\le6t+2d+2H,\quad |c-a|\le6t+2d+2H,                       \tag{FR7}
\]
up to the original orientation and signs of its two long differences.
These are strips of original labelled differences, before primitive
normalization; dividing a zero-first-coordinate normal by its gcd need
not preserve an endpoint strip. There are \(O(a(t+d+H))\) such
integer directions. Comparisons
within an arm or involving a bridge add \(O(at)\), as in the
original hook inventory. Include the two price-plane pencils with
\(r-p\in[-d_0,d_0]\) and the sector plane.

Every actual wall is represented by a coherent minimum-face exchange,
so these planes refine the entire positive fan. Their arrangement has
\(O(a^2(t+d_0+H_0)^2)\) chambers, including all parallel pieces
and subdivisions. Original normals have height \(O(a)\); primitive
cross products and three-ray sums give the stated integer weight list.
This uses the established state-polytope and matroid-face dictionary,
not an assumption that one supporting direction is one wall.

## The previous growth probes are calibration families

The nearby and truncated-dense axis probes have polynomial degree
\(O(t)\). They were already covered by the exact rational-leading
theorem, taking \(Q=1\). Their large expanded inverse degrees and
prefix gaps do not reopen that count question.

The separated probes have, on each arm, a low polynomial of degree
\(O(t)\) and additional terms among \(q^{A-4},\ldots,q^{A-1}\).
They meet (FR1) with \(h\le4\le t\), over every field, uniformly
in the low and terminal coefficients while the selected leading term
survives. Therefore (FR3) bounds this entire growing family by
\(O(n^2)\). The finite fans serve as algorithm calibrations. Their
counts do not establish an asymptotic growth rate or a polynomial
lower excess.

## A degree gate identifies a remaining axis stratum

There is a simple obstruction to repeatedly using the compact-filtered
class. Let a scalar polynomial \(W(q)\) have degree \(0<D<A\),
and suppose it admits (FR1)--(FR2), even with mixed \(P,Q\).
Then
\[
                           d+h\ge\min(D,A-D).                 \tag{FR8}
\]
Indeed set \(v=0\). The pure leading terms ensure that the resulting
polynomials \(P_0,Q_0\) have degrees \(p,r\), and \(Q_0\) is
nonzero. If both \(d<D\) and \(d+h<A-D\), the ordinary polynomial
\(Q_0W-P_0\) has degree at most \(D+d<A-h\), while all its
coefficients below \(A-h\) vanish. It is zero. But then
\(p=D+r\ge D>d\), a contradiction.

The bound is sharp for \(W=1+q^D\): either take \(P=W,Q=1,h=0\),
or take \(P=Q=1,h=A-D\). A large filtered parameter is not evidence
of a large fan. For example, take \(A=2D\), \(D\ge t\), zero
outer images and \(W_Y=W_Z=1+q^D\). Partition each arm into
two-row chains with long step \(D\). Each head joins exactly one
Y chain to one Z chain; distinct heads use disjoint chains. After the
common cost \(iu+j\) is subtracted, every joined component has the
same seven original columns and the same price formulas. The unjoined
Y and Z chains each have four columns and one fixed price pattern.
There are only three constant-size component types. Their common
comparison arrangement therefore has \(O(1)\) cones, even as the
minimum filtered degree \(D\) grows. This direct argument concerns
these pure-axis, length-two chains; it does not silently invoke the
earlier theorem whose monomial steps have positive short degree.

For an axis polynomial of degree \(D\) with \(W(0)\ne0\),
every prefix minor with \(t+D<r<A\) vanishes: its Toeplitz matrix
has entries \(W_{r-t+i-j}\), with zero diagonal and zero entries
below it. Hence its longest prefix gap is at least
\(\max(0,A-t-D-1)\). Also a nonconstant polynomial unit satisfies
\(\deg W^{-1}\ge A-D\), since an ordinary product of degree less
than \(A\) cannot equal one unless both factors are constant.

Consequently dense degree-\(\lfloor A/2\rfloor\) axis strata in
the intermediate aspect have growing forward/inverse degrees, forced
large prefix gaps and growing filtered degree. They provide a specific
remaining class for the collective crossing problem. These are
necessary scope checks, not a fan lower bound or an assertion that
every coefficient specialization escapes every earlier theorem.

## Verification and prior-art boundary

~~~bash
python3 -B scripts/run_hook_filtered_rational_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker reconstructs original commuting quotients, verifies that
early denominator residuals are nonzero, checks the surviving exact
pure and outer identities, and reduces selected targets against all
strictly cheaper columns. Entire transformed prefixes independently
verify the unique separating row. Complete original fans and complete
compiled axis fans supply separate finite controls. Cases include
mixed terminal coefficients, arbitrary outer maps, nonunit numerators,
rational coefficients and small characteristic.

Formal rational approximation, filtered multiplication and triangular
elimination are classical. The related Toeplitz/Padé determinant and
zero-block framework is discussed by
[Jedynak--Gilewicz, Sections 1--2](https://doi.org/10.1155/2013/185648).
The coherent minimum-face input is classical
[Feichtner--Sturmfels, Section 2](https://arxiv.org/pdf/math/0411260v1).
This note gives the explicit terminal-error allowance in the existing
window proof; exact strongest-current priority is unestablished. It
does not bound arbitrary large filtered degree, remove the long-leading
hypothesis or settle the unrestricted extremal objective.
