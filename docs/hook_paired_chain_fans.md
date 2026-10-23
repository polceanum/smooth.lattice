# Paired binomial chains have a sharp quadratic thickness count

Long inverse propagation alone does not make the mixed-leading unit
candidate large. With zero outer images, binomial central arms whose
long steps are at least the thickness have only a quadratic number of
initial ideals in their maximal chain lengths. The bound is independent
of the number of chains and of their positions along the long arms.
It is sharp along an infinite actual family.

Nonzero outer images require a separate argument. Exact original
quotients below have coherent walls outside every support plane of the
paired-chain inventory. Their full fans are also missed by the common
weight list for zero outer images. Thus the simplified chain model is
not a replacement for the remaining arbitrary-outer problem. On the
balanced sector, those couplings are retained by an exact bounded-column
extension, giving a separate count linear in arm length at fixed thickness.
Neither an unrestricted extremal bound nor exact current novelty is established.

## The class and its quantitative consequence

Use the [compatible hook](hook_bridge_reduction.md), with
\[
 a\ge2t,\quad t\ge2,\quad A=a-t,\quad n=4at-2t^2.
\]
Both outer images \(f_0,f_2\) are zero. In the two central arms take
\[
 W_i=\alpha_i+\beta_i q^{\nu_i}v,
 \qquad \alpha_i\beta_i\ne0,\quad t\le\nu_i<A,
 \qquad i\in\{Y,Z\}.                                      \tag{PC1}
\]
Work over any field; the two long steps and coefficients may differ.
Put
\[
 h_i=\min\{t,\lceil A/\nu_i\rceil\}.
\]
Then the whole positive initial fan satisfies
\[
                       g(I)=O((h_Y+h_Z)^2)=O(t^2).         \tag{PC2}
\]
For fixed \(a,t,\nu_Y,\nu_Z\), its fan is the same over all fields
and all nonzero choices of the four displayed coefficients. One common
complete positive integer weight list has \(O((h_Y+h_Z)^2)\) members
of height \(O(a^2)\). These are consequences for the full fan and its
weight compilation, not claims about optimal bit complexity.

If \(\nu_Y=\nu_Z=\nu\ge2t\) and \(a\ge t\nu+2t\), then
\[
                            g(I)=\Theta(t^2).             \tag{PC3}
\]
For example \(\nu=2t,\ a=2t^2+2t\) gives \(n=8t^3+6t^2\)
and \(g=\Theta(n^{2/3})\). The expanded inverse degree can have
order \(a\). It does not supply a lower bound on fan growth.

## The original matrix splits into matched chain pairs

Partition the Y rectangle into chains with step \((\nu_Y,1)\),
and the Z rectangle with step \((\nu_Z,1)\). An input \((q,j)\)
has chain root
\[
 (q-k\nu_i,j-k),\qquad k=\min(\lfloor q/\nu_i\rfloor,j).
\]
Because \(\nu_i\ge t\), each input \((q_0,j_0)\) with
\(q_0<t\) is a root, and distinct such inputs lie in distinct
chains. The central bridge with original label
\((0,t+q_0,t+j_0)\) joins precisely the Y chain rooted at
\((q_0,j_0)\) to the Z chain rooted at \((j_0,q_0)\).
There is one bridge in each of these \(t^2\) paired components.
Every other chain is an independent component. Zero outer columns
are loops. This is an actual row-support direct sum of the original
residue matroid, not just a decomposition of its central image.

On a chain with rows \(e_0,\ldots,e_{m-1}\), the upper columns
are \(U_r=e_r\). The lower columns are
\(L_r=\alpha e_r+\beta e_{r+1}\), with
\(L_{m-1}=\alpha e_{m-1}\). Rescale rows to obtain the usual
path incidence columns and edges to a common ground vertex.
Thus their minimum bases are classical minimum spanning trees of
a path with a star. Every rescaling is invertible in every field.

An unpaired chain has all its \(L_r\). A paired chain has only
\(L_1,\ldots,L_{m-1}\); its missing \(L_0\) is the projection
of the one shared bridge. In particular \(U_0\) is a coloop of
the pure chain, and the remaining rows form another path with a star.

## Each chain has only three minimum-basis forms

Normalize the positive weight to \((\lambda,u,1)\), and set
\[
 c_Y=\lambda-t,\quad s_Y=\nu_Yu+1,
 \qquad c_Z=\lambda-tu,\quad s_Z=\nu_Z+u.                  \tag{PC4}
\]
After subtracting the same affine cost from every column of one
component, an unpaired chain has costs
\(L_r:rs\), \(U_r:c+rs\). Subtracting a common cost preserves
all minimum bases because their cardinalities agree.

Away from the displayed equalities, its unique minimum basis is
\[
 \begin{array}{c|c}
 c<-s&U_0,\ldots,U_{m-1}\\
 -s<c<(m-1)s&U_0,L_0,\ldots,L_{m-2}\\
 c>(m-1)s&L_0,\ldots,L_{m-1}.
 \end{array}                                               \tag{PC5}
\]
The first two descriptions coincide when \(m=1\). To check (PC5),
if \(c<-s\), each path edge costs more than both adjacent star
edges, so is discarded by greedy choice. If \(c>-s\), each path
edge is the unique cheapest edge crossing the cut between its suffix
of path vertices and the ground plus its prefix: all later star edges
and the final ground edge cost more. Thus it belongs to the minimum
tree. Contracting the path leaves a choice between
the cheapest star edge \(U_0\) and the final edge
\(L_{m-1}\). Their costs give the other threshold.

For a paired pure chain, retain its compulsory \(U_0\) and apply
(PC5) to the tail starting at row one. The bridge projection is
\(\alpha e_0+\beta e_1\) when \(m\ge2\), and
\(\alpha e_0\) when \(m=1\). Its fundamental circuit in that
pure minimum basis has the following most expensive member:
\[
 \begin{array}{c|c}
 m=1&U_0\quad\text{of cost }c\\
 m\ge2,\ c<(m-2)s&U_1\quad\text{of cost }c+s\\
 m\ge2,\ (m-2)s<c<(m-1)s&L_{m-1}\quad\text{of cost }(m-1)s\\
 m\ge2,\ c>(m-1)s&U_0\quad\text{of cost }c.
 \end{array}                                               \tag{PC6}
\]
All relevant coefficients are nonzero: the representation on the
tail path uses consecutive nonzero powers of \(-\beta/\alpha\).
There is no characteristic-dependent cancellation.

The paired component's bridge has normalized cost zero. Start with
the two pure minimum bases. If the most expensive member of the
bridge's fundamental circuit costs more than zero, insert the bridge
and remove that member. Otherwise retain the pure bases. This is
classical greedy single-element insertion. One direct proof is that
an old column rejected by greedy remains in the span of cheaper old
columns after an element is added. Thus the new optimum is contained
in the old optimum plus the new element; the circuit gives precisely
the allowable removals. This also proves (PC6) determines the entire
paired optimum, rather than only a possible exchange.

The two pure cost offsets agree with the bridge offset: its original
cost is \((t+q_0)u+t+j_0\), equal to the Y root offset plus
\(t\) and to the Z root offset plus \(tu\). Root locations
therefore cancel from (PC4)--(PC6). Only the two chain lengths remain.

## Count restricted walls rather than their complete supporting lines

All changes within (PC5)--(PC6), and all comparisons with zero, lie on
\[
 \lambda=t+ks_Y\quad(-1\le k<h_Y),\qquad
 \lambda=tu+ks_Z\quad(-1\le k<h_Z).                       \tag{PC7}
\]
Comparing an upper critical member on one side with a lower critical
member on the other adds the four families
\[
 \lambda=t+js_Z-\epsilon s_Y,
 \quad 1\le j<h_Z,\ \epsilon\in\{0,1\},
\]
\[
 \lambda=tu+is_Y-\eta s_Z,
 \quad 1\le i<h_Y,\ \eta\in\{0,1\}.                      \tag{PC8}
\]
Comparing two upper critical members gives at most four vertical
lines, from \(c_Y+\epsilon s_Y=c_Z+\eta s_Z\).
There are \(L=O(h_Y+h_Z)\) full lines altogether.

The remaining comparison is \(is_Y=js_Z\), with
\(1\le i<h_Y,\ 1\le j<h_Z\). It can matter only when both
lower critical members are active. Its support is therefore restricted to
\[
 \begin{split}
 u&=\frac{j\nu_Z-i}{i\nu_Y-j},\\
 \lambda&\in[t+(i-1)s_Y,t+is_Y]
             \cap[tu+(j-1)s_Z,tu+js_Z].
 \end{split}                                               \tag{PC9}
\]
Empty or point intersections add no wall segment. Positive denominators
and numerators follow from \(\nu_Y,\nu_Z\ge t\).

There are \(O(h_Yh_Z)\) such vertical segments. Each has length
at most \(\min(s_Y,s_Z)\) in the \(\lambda\) coordinate.
At a fixed \(u\), each of the six nonvertical line families in
(PC7)--(PC8) has spacing either \(s_Y\) or \(s_Z\).
Consequently a segment meets at most twelve of these lines in its
relative interior. Vertical segments do not cross each other
transversely; collinear overlaps only reduce the count. The full
vertical lines likewise add no transverse cuts on them.
Thus the complete arrangement has at most
\[
                1+\frac{L(L+1)}2+13(h_Y-1)(h_Z-1)
                         =O((h_Y+h_Z)^2)                 \tag{PC10}
\]
open regions. The minimum-basis formulas are constant in each region,
proving (PC2). Extending every segment to its full supporting line
would lose this argument and give an unnecessary higher power.

Every segment endpoint lies on a full line from (PC7). Within a cell
of the full-line arrangement, it therefore cuts a complete vertical
chord or misses the interior. The refined cells remain convex.
All homogeneous comparison normals have integer height \(O(a)\):
\((h_i-1)\nu_i<A\). In the positive projective triangle, extreme
rays are cross products of two such normals. Summing three independent
primitive rays gives an interior integer weight of height \(O(a^2)\).
This proves the common-list statement. The formulas show coefficient
and field independence without assuming that arbitrary matroids can
be transferred between fields.

## Sharpness uses genuine independent components

Assume the hypotheses of (PC3). For every \(1\le m\le t\),
the root \((A-m\nu,0)\) lies in an unpaired chain of length \(m\):
its long coordinate is at least \(t\). Its actual wall from (PC5) is
\(\lambda=t+(m-1)s_Y\). The transposed chain gives the actual
wall \(\lambda=tu+(m-1)s_Z\).

Take \(k,l\in\{1,\ldots,t-1\}\). The Y and Z walls with these
indices meet at
\[
 u=\frac{l\nu-t-k}{k\nu-t-l}>0,
 \qquad \lambda=t+k(\nu u+1)>0.                          \tag{PC11}
\]
Both numerator and denominator are at least one when \(\nu\ge2t\).
Different Y lines do not meet for positive \(u\), nor do different
Z lines, so these are \((t-1)^2\) distinct proper intersections.
The corresponding chains are actual direct summands. Their independent
minimum-basis changes give four local choices at every intersection.
The common refinement therefore has \(\Omega(t^2)\) regions.
Together with (PC2), this proves (PC3).

## Balanced weights retain arbitrary outer maps in a fixed extension

Keep the central hypotheses (PC1), but now let both compatible outer
images be arbitrary. Write \(g_{\rm bal}\) for the number of open
fan cells meeting \(1/2<w_y/w_z<2\). There is an exact description
on \(O((h_Y+h_Z)^2)\) regions by fixed represented matroids of
dual rank at most \(8t^2\). In particular
\[
 g_{\rm bal}
 =O\bigl(a t^3(h_Y+h_Z)^2\bigr)
 =O\bigl(n t^2(h_Y+h_Z)^2\bigr)
 =O(nt^4).                                                 \tag{PC12}
\]
This is a sector count. The established
[box ceiling](hook_bridge_reduction.md) and other saved bounds should
still be retained when smaller; no stronger whole-fan or critical-aspect
endpoint follows merely from (PC12).

First delete outer sources whose original input along their own arm
is at least \(4t\). Every such source is pure. Since its central
operator is a unit, its vector belongs to the span of the original
pure central sources. With arm ratio \(v>1/2\), all those central
columns cost at most \((a-1)v+2t-1\), strictly less than the outer
cost \((a+i)v+j\) for \(i\ge4t\). Thus these deletions preserve
all minimum bases throughout the balanced sector, for every value of
\(w_x\). Retain the fixed outer label set \(H\), of size
\[
                     k=|H|\le2(4t)t=8t^2.                \tag{PC13}
\]
Its actual vectors keep every compatible outer coefficient, including
their arbitrary projections into the chains. They are not replaced
by ranks or support graphs.

Choose a region of the zero-outer fan and let \(B_0\) be its unique
minimum residual basis. This is the explicit basis computed above.
Every old column outside \(B_0\) is a linear combination of strictly
cheaper members of \(B_0\) throughout this entire region: otherwise
its fundamental circuit would permit a cheaper basis. The new outer
columns cannot invalidate that dependence. Hence every minimum of
the full quotient is a basis in the one fixed matrix on
\(B_0\cup H\). All original labels and costs are retained. The
matrix has rank \(R=2At\) and \(R+k\) columns; its dual has
rank \(k\). Complementing a minimum basis gives a maximum dual basis.
This is an exact matroid encoding of the whole cone patch, not a
claim that the dual is an algebra quotient or has only \(O(t^2)\)
labels. Its ground still has \(R+k\) labels.

For completeness, the classical lattice-vertex argument gives the
count in (PC12). A dual basis has \(k\) original labels in
\([0,1]\times[0,2a]\times[0,2a]\), so every basis profile lies in
\([0,k]\times[0,2ak]^2\). Add the unit cube to its integer convex
hull. The sum has full dimension, retains at least as many vertices
by normal-fan refinement, and has volume \(O(a^2(1+k)^3)\).
The classical three-dimensional lattice-vertex bound therefore gives
\(O(a(1+k)^{3/2})=O(at^3)\) possible unique optima per encoding.
A convex original region cannot split one dual normal cone into
several pieces. Multiplying by the \(O((h_Y+h_Z)^2)\) backbone
regions proves (PC12), including lower-dimensional profile hulls.
The vertex exponent is classical: see
[Andrews](https://doi.org/10.1090/S0002-9947-1963-0143105-7) and
the primary proof retelling in
[Dillon, Theorem 1](https://arxiv.org/pdf/2606.30856).

## Verification and the remaining outer coupling

~~~bash
python3 -B scripts/run_hook_paired_chain_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker reconstructs original quotients and commuting multiplication,
verifies the exact support partition, compares the explicit component
formulas with independent field elimination, and enumerates complete
original fans. It constructs the partial-wall arrangement and its common
integer list, checking every segment's constant crossing bound separately.
Coefficients, fields and the two long steps vary.
Original unpaired summands also give four independently verified strict
minimum bases at each tested lower-bound junction. The arbitrary-outer
controls test every late-source loop, the complete original cone patch
against the fixed extension, and maximum dual bases against complementary
minimum bases, including characteristic zero.

With nonzero outer maps, two saved original fixtures have 192 cones
instead of 18 at \((t,a,\nu_Y,\nu_Z)=(2,7,2,2)\), and 541
instead of 28 at \((3,11,3,3)\). The report gives actual coherent
exchange witnesses outside every supporting plane in the zero-outer
inventory, as well as missing cones under its common list. These finite
guards demonstrate the failure of that encoding. They establish no
asymptotic outer lower bound, let alone a polynomial excess over the
generic corner-cut scale.

Graphic chain minima, greedy insertion, matroid faces and elementary
arrangement counting are classical. The face and exchange dictionary is
in [Feichtner--Sturmfels, Section 2](https://arxiv.org/pdf/math/0411260v1).
The exact current priority of this combined class result remains
unestablished. The remaining candidate must retain nonzero outer maps
that couple the chains; steps below \(t\) can also join multiple
central bridges within one chain. The unrestricted sharp \(G_3\),
an actual polynomial fan excess and the original final-paper requirements
remain open.
