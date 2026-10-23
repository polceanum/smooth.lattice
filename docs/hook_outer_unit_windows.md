# Outer images add a translated hyperbola to coherent unit windows

The restriction on the chosen arm's outer images can be removed. An
actual unit inverse still confines coherent exchanges to the two ends
of that arm, except for a translated hyperbolic inventory. The whole
compatible hook fan has a quadratic polylogarithmic bound whenever
one central unit has bounded forward and inverse propagation relative
to thickness. All opposite and outer coefficients are allowed.

An infinite actual family with the identity as its chosen central
operator shows why the outer case needs an additional inventory.
One quotient has an increasing number of coherent outer-source
partners, and their translated products have order \(at\) at fixed
thickness. This is not a polynomial excess over the generic fan scale.
The unrestricted sharp count, exact strongest-current novelty and the
original final-paper requirements remain open.

## Full compatible-hook bound

Use the [original compatible flat hook](hook_bridge_reduction.md), with
\[
 t\ge2,\quad a\ge2t,\quad A=a-t,\quad n=4at-2t^2.
\]
Assume only that the central Y operator \(W\) is a unit in
\(R=K[q,v]/(q^A,v^t)\). Every compatible coefficient of both outer
images and the entire opposite central image is arbitrary. Put
\[
 D=\max\{\deg_qW,\deg_qW^{-1}\},\qquad
                     H=1+\log(a/t).                       \tag{OW1}
\]
Over every field,
\[
 \begin{split}
 L&=O\bigl(a(D+tH)\bigr),\\
 g(I_\phi)&=O\bigl(a^2(D+tH)^2\bigr)
           =O\bigl(n^2(D/t+H)^2\bigr).                    \tag{OW2}
 \end{split}
\]
Here \(L\) is the size of one common original comparison-plane
inventory containing every positive wall. Either unit central arm
can be chosen, taking the smaller resulting bound.

For fixed \(a,t,D_0\), the same inventory works for all fields and
all compatible coefficients with \(D\le D_0\). It gives a common
complete positive integer weight list of size
\(O(a^2(D_0+tH)^2)\) and height \(O(a^2)\).
No optimal construction time or total quotient-rank cost is asserted.

In particular the support-cone class
\[
 W=\gamma+\sum_{j\ge1,\ 0\le i\le Cj}c_{ij}q^iv^j,
                         \qquad\gamma\ne0,               \tag{OW3}
\]
has \(D\le C(t-1)\), because the cone is closed under the finite
unit-inverse series. Its full compatible hook fan is
\(O(n^2(1+C+\log(a/t))^2)\), with arbitrary outer images.
The number of independent monomial steps is unrestricted.

## Only a prefix must split

Normalize the positive objective to \((\lambda,u,1)\).
A coherent exchange means that two minimum bases differ by one pair;
additional cost ties are allowed. When its common cost is \(\tau\),
only the restriction to columns of cost at most \(\tau\) is needed.
After contracting strictly cheaper independent columns, a direct-sum
split of that entire prefix forbids a coherent exchange across its
summands. Indeed every minimum basis contains a basis of the full
\(\tau\)-prefix and of its strictly cheaper subspace. The fixed ranks
in the direct sum cannot change by one exchange between summands.
The full represented ground need not split.

Write \(d_+=\deg_q W\), \(d_-=\deg_qW^{-1}\), and \(S=t+D\).
The [unit-window price identity](hook_coherent_unit_windows.md)
remains valid without restrictions on outer images: an upper Y input
\(q\ge t\) is a combination of available pure central columns.
It is a loop modulo strictly cheaper columns whenever
\(\lambda>2t-1+d_-u\). Thus eligibility requires
\[
                         \lambda\le2t-1+d_-u.             \tag{OW4}
\]
No row operation changes any original label or price.

A nonzero pure Y column is either an upper identity, a pure central
source, or an outer \(y^a\)-source with input shift \(i\ge t\).
The \(z^a\)-source has no pure Y columns. This is the original
three-box partition, not a specialization of its groups.

## Central exchanges stay near the ends

Suppose the Y column is upper or central, with original long label
\(M<a\), and long input \(q=M-t\ge t\). Then every coherent
pure cross-arm exchange satisfies
\[
 uM\le(4t+2D)u+4t
 \quad\text{or}\quad
 u(a-M)\le(4t+2D)u+4t.                                  \tag{OW5}
\]
Here is a prefix proof, retaining every outer image. Assume both
inequalities fail, and write \(B=4t+2D\).

For an upper column use the original row-identity chart. Its price
satisfies (OW4), so
\(\tau\le u(M+D)+3t<au\).
Consequently no \(y^a\)-source, including its bridge box, belongs
to the \(\tau\)-prefix. Contract all strictly cheaper Y upper
identities. Central bridges have Y support only at long labels at
most \(2t-1+d_+\), and these projections are killed by the first
failed inequality. Every central source in the prefix has input
long label at most \(M+D+2t/u\), hence output long label at most
\(M+2D+2t/u<a-t\). Every upper Y identity in the prefix also
has long label less than \(a-t\).

By compatibility, the Y projection of any \(z^a\)-bridge is
supported at long labels at least \(a-t\). Therefore the whole
remaining \(\tau\)-prefix splits into Y rows below \(a-t\) and
the far-Y-plus-Z rows. The exchanged Y identity lies in the former
part; an opposite pure column lies in the latter. This is impossible.

For a lower central column apply \(W^{-1}\) on the Y rows only.
Pure central sources with \(q\ge t\) become identities on the
Y tail. Again \(\tau\le uM+2t<au\), so no \(y^a\)-source is
in the prefix. Contract the strictly cheaper tail identities.
The tail outputs of every head-input upper column are killed, since
the first failed inequality exceeds
\((2t-1+d_-)u+t-1\). Central bridges have only head Y inputs.
Every upper column in the prefix has input long label at most
\(M+2t/u\), and inverse output long label less than \(a-t\).
The normalized \(z^a\)-bridge projections stay at long labels
at least \(a-t\), because a polynomial unit inverse preserves that
long-coordinate ideal. The prefix now splits into the middle Y rows
\(2t\le y<a-t\) and the head-Y-plus-far-Y-plus-Z rows. The proposed
exchange again crosses summands. This proves (OW5).

If an upper input has \(q<t\), its label already has \(M<2t\).
Thus it belongs to the short inventory without a window argument.

## Outer sources and the translated direction inventory

A pure outer Y source has label \((0,a+i,j)\), with
\(t\le i<A\), \(0\le j<t\), and lies in the Y tail. The pure
central columns span that entire tail, since \(W\) is a unit.
Their largest cost is at most \((a-1)u+2t-1\). Hence eligibility
of the outer source requires
\[
                (i+1)u+j\le2t-1,\qquad iu\le2t.          \tag{OW6}
\]
Otherwise those central columns are all strictly cheaper and make
it a loop.

For a pure cross-arm pair write its chosen Y label as
\((e_x,M,s)\), and its opposite Z label as \((f_x,r,K)\), with
\(r,s<2t\). Equal costs give
\[
 K=(e_x-f_x)\lambda+u(M-r)+s.                            \tag{OW7}
\]
For a central or upper Y column with \(q\ge t\), (OW4) implies
\(K\le u(M+D)+4t\). For an outer Y column, (OW6) gives
\(K\le au+3t\).

Set \(b=M-r\), \(d=K-s\). Every actual pure cross-arm wall lies
in the following coefficient-independent inventory:
\[
 \boxed{\quad
 b\le16S\quad\text{or}\quad d\le16S\quad\text{or}\quad
 |b-a|\le16S\quad\text{or}\quad |b-a|d\le16at.
 \quad}                                                   \tag{OW8}
\]
All coordinates stay in the original bounded label box; negative
\(b,d\) are included by the first two alternatives.

To verify the constants, the first alternative of (OW5), when
\(u\ge1\), bounds \(M\le8S\), and when \(u<1\) bounds
\(K\le12S\). The second alternative bounds
\(a-M\le8S\) if \(u\ge1\). If \(u<1\), put
\(\delta=a-M\). When \(\delta>B\), it gives
\[
                    (\delta-B)(K-4S)\le4at
\]
whenever the second factor is positive. Outside the three strips in
(OW8), \(\delta-B\ge(a-b)/2\) and \(K-4S\ge d/2\), giving
\((a-b)d\le16at\). The outer-source case follows from
\(i(K-3t)\le2at\): outside the same strips,
\((b-a)d\le4at\). These arguments include same-layer exchanges.

The three strips have \(O(aS)\) possible difference vectors. For
the translated hyperbola, \(|b-a|,d=O(a)\) and their product is
\(O(at)\). Split the first coordinate at a constant multiple of
\(t\); the remaining harmonic sum gives
\(O(at(1+\log(a/t)))\) pairs. Add the \(O(at)\) arm and bridge
comparison planes. An arrangement of \(L\) central planes has
\(O(L^2)\) positive chambers, proving (OW2).

The common list follows from that same inventory: primitive extreme
rays are cross products of original normals of size \(O(a)\), and
summing three independent rays gives an interior integer weight of
height \(O(a^2)\). The usual greedy exchange-edge argument ensures
that every actual positive wall lies in the inventory, even with
coincident comparisons. These classical matroid face facts are
reviewed in [Feichtner--Sturmfels, Section 2](https://arxiv.org/pdf/math/0411260v1).

## An infinite actual outer-source family reaches the product scale

Fix \(m\ge3\), and let
\[
 t=2,\quad a=8m^2,\quad A=a-2,\quad K=m,\quad d=A-1-K.
\]
In characteristic zero or \(p>2A\), take
\[
 \begin{split}
 f_0&=y^2(1+y)^{A-1}+yz^{a-1},\\
 f_1&=y^2+z^2(1+z^d),\qquad f_2=0.
 \end{split}                                               \tag{OW9}
\]
All expressions mean their residues in the hook quotient. The map is
compatible by the original allowed supports. Its chosen central
operator is \(W_Y=1\), so \(D_Y=0\).

In this one fixed quotient, for every integer \(m\le s\le2m\),
put
\[
 \begin{split}
 E_s&=(0,a+s,0),\quad F=(1,1,K+2),\\
 u_s&=\frac{K+2}{a-s-3/2},\qquad
 \lambda_s=(2s+1/2)u_s.
 \end{split}                                               \tag{OW10}
\]
The positive costs of \(E_s,F\) agree. Equivalently an integer
witness is
\[
              ((4s+1)(m+2),\ 2(m+2),\ 2a-2s-3).
\]
For all these parameters,
\(0<\lambda_s<1\), \(0<\lambda_s-2u_s<1\), and \(su_s<1\).
The cheap upper Y identities in short degree zero leave exactly the
last \(s\) long rows. Let \(G=(1+q)^{A-1}\) and
\(g_i=q^iG\) restricted to these rows.
The other cheaper outer sources give \(g_1,\ldots,g_{s-1}\).
The head source is \(g_0+z_{A-1}\), where \(z_{A-1}\) is the
upper Z row with short degree one. All central Y identities that
are cheaper are already eliminated by cheaper upper identities.

The cheap Z sources and upper identities give
\(F+z_{A-1}\), and kill all their other relevant entries.
Indeed the Z sources have inputs through \(K\), while the next
input is too expensive. The upper row at \(K\), short degree one,
is exactly \(F\). Central head helpers have cheaper Y entries;
all other cheap columns are on disjoint coordinates or have the
already listed projections.

The two square matrices on the last \(s\) rows with columns
\((g_0,\ldots,g_{s-1})\) and \((g_1,\ldots,g_s)\)
have nonzero determinants. They are the classical binomial Toeplitz
products
\[
 \det\left[\binom{A-1}{r+i-j}\right]_{0\le i,j<s}
 =\prod_{i=0}^{s-1}
       \frac{(A-1+i)!\,i!}{(r+i)!\,(A-1-r+i)!},            \tag{OW11}
\]
with \(r=A-s\) and \(r=A-s-1\), respectively.
See [Krattenthaler, Section 2.3, (2.17)](https://www.mat.univie.ac.at/~slc/s/s42kratt.pdf#page=11).
All factorial arguments are less than \(2A\). Modulo
\(g_1,\ldots,g_{s-1}\), both \(g_0,g_s\) are nonzero in a
one-dimensional space. The head source and Z source therefore make
\(E_s,F\) nonzero and parallel modulo the full strictly cheaper span.
An explicit functional on the last Y rows and the two relevant Z
rows annihilates every cheaper original column, verifying nonloops.
The strict inequalities persist on an open segment of the tie line,
which supplies an actual coherent wall after avoiding other tie planes.

The \(m+1\) distinct primitive directions are
\[
                     (1,1-a-s,K+2).
\]
Their translated products are
\((s-1)(m+2)=\Theta(a)=\Theta(at)\) at fixed thickness.
Thus a uniform product allowance \(O(t^2)\) is false with arbitrary
outer images, even when the chosen central operator is the identity.
This does not prove that the logarithm in (OW2) is necessary, and the
\(\Omega(\sqrt n)\) direction count is below the original extremal target.

## Verification and limits

~~~bash
python3 -B scripts/run_hook_outer_unit_window_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

Checks use original compatible quotient columns and commuting
multiplication, both invertible row charts, whole prefix direct-sum
identities, complete original fans, and a common integer list reused
across fields. Separate binomial minors, explicit source relations and
a full streamed cheaper-prefix functional verify the outer family at
larger lengths. Finite checks support the uniform proofs; they do not
establish a global extremal exponent or exact novelty.

A power-exceeding family in the compatible-hook class with a chosen
unit central arm must now have
\(D/t=\Omega(n^{\epsilon/2})\) if
\(g\ge c n^{2+\epsilon}\). Nonzero outer images alone no longer
provide an escape from the bounded-propagation closure. Nonunit central
arms, long propagation on every available unit arm, the unrestricted
sharp count and strongest-current priority remain unresolved.
