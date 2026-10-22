# Generic hook couplings can have linear positive colour width

Put \(a=t^2+1\), \(M=(y^a,y^tz^t,z^a)\),
\(C=K[y,z]/M\), \(m=2ta-t^2\), and \(n=2m\).
For a flat square-zero extension \(A_\phi\) of this planar quotient,
let \(k(B)\) be the number of upper-layer monomials in an actual initial
standard set \(B\). This is its \(x\)-moment, since \(x^2=0\).
Its **positive colour width** is the maximum of \(k(B)\) minus its minimum,
over all global monomial initial ideals.

**Theorem.** For every \(t\ge32\), a nonempty Zariski open locus in the
complete independent extension parameters has
\[
 \min_B k(B)=t^2,\qquad \max_B k(B)>m/8,\qquad
                    \operatorname{width}_x\ge n/32.       \tag{1}
\]
The locus contains an explicit integral coefficient choice over every
field. Thus a sublinear positive colour width, in particular \(O(t^2)\),
cannot be imposed on generic growing hooks to improve a profile-volume
estimate. The width is \(\Theta(n)\), since every standard set has at
most \(m\) upper monomials.

The proof needs only one unit minor and does not assume that the full
support-matching model survives at growing thickness. It does not bound
the number of generic cones. A sparse witness below has exactly twelve
cones despite its linear width. Neither statement settles unrestricted
\(G_3\), gives a polynomial lower separation, or establishes novelty.

## The open conditions and the large-price endpoint

Use the [complete independent-parameter model](generic_hook_fans.md).
Write \(f_0,f_1,f_2\) for the images of \(y^a,y^tz^t,z^a\).
Let \(c_0=[y^t]f_0\), \(c_2=[z^t]f_2\), and impose
\(c_0c_2\ne0\). Their images generate the ideal
\[
                  T=(y^t,z^t)C,\qquad \dim T=m-t^2.       \tag{2}
\]
Indeed \(f_0,f_2\) are \(y^t,z^t\) multiplied by a two-by-two
matrix which is invertible modulo \((y,z)\). Its determinant is a
unit in the local Artinian algebra \(C\), so that matrix is invertible.
All three images belong to \(T\), proving (2).

The plain monomials therefore span a space of dimension \(2m-t^2\),
and the whole quotient has dimension \(2m\). A vector basis needs at
least \(t^2\) upper columns, and some vector basis uses exactly that
many. The finite ground consists of the diagram of \(M^2\) in the
plain layer and the diagram \(D\) of \(M\) in the upper layer.
Every variable label there has planar degree less than \(4a\).
Consequently prices \((w_x,w_y,w_z)=(4am+1,1,1)\) force the minimum
upper count \(t^2\). A sufficiently small generic perturbation gives
an open positive initial cone with that count.

## A unit-minor competitor at the balanced price

Set
\[
 L=a-3t-1=t^2-3t,\quad h=2t+1+\lfloor L/2\rfloor,\quad
 r=a-h-t=\lceil L/2\rceil,\quad \ell=2tr.                \tag{3}
\]
For \(t\ge4\), one has \(2t<h<a-t\).
At the integral specialization
\[
             f_0=y^t,\qquad f_1=y^h+z^h,\qquad f_2=z^t,  \tag{4}
\]
all adjacent syzygies hold, in every characteristic.

The plain \(D\) columns are compulsory coloops: every other nonzero
column has only an upper component. Start with the vector basis
\(D\sqcup xD\). For
\(t\le i<a-h\), \(0\le j<t\), replace
\[
 xy^{h+i}z^j\ \longmapsto\ y^{t+i}z^{t+j},\qquad
 xy^jz^{h+i}\ \longmapsto\ y^{t+j}z^{t+i}.               \tag{5}
\]
The new plain column in the first replacement is exactly the old
upper column at (4): its other image term \(xy^iz^{h+j}\) vanishes
modulo \(xM\). The second replacement is symmetric. The \(\ell\)
rows are distinct, so the corresponding residual minor \(\Delta_t\)
has determinant \(\pm1\) at (4). Impose also \(\Delta_t\ne0\).
These replacements then give an available vector basis throughout
that open locus. It need not be a monomial downset or an initial set;
only its availability as a competitor is used.

At price \((1,1,1)\), subtract the compulsory plain \(D\) cost.
The all-upper reference cost is exactly
\[
 C_0=\sum_{v\in D}(1+|v|)=ta(a+t)-t^3
                         =t^5+t^4+t^3+t^2+t.            \tag{6}
\]
Every replacement in (5) saves \(h-2t+1\). Thus the minimum residual
cost is at most
\[
 C_0-G,\qquad G=2t\lceil L/2\rceil(\lfloor L/2\rfloor+2)
                           \ge\tfrac t2(L^2-1).         \tag{7}
\]
The conditions \(c_0c_2\Delta_t\ne0\) define a nonempty open set:
(4) makes this product a unit over every field. No assertion about
all coherent support minors is needed.

## Too few upper columns costs more

Let \(U\) be the selected plain columns outside \(D\), of size
\(q=m-k(B)\). For an integer \(s<a\), at most \(2t(s+1)\)
possibly nonzero such monomials have planar degree at most \(s\).
The pure-power generators have degree \(a\). All other candidates
are divisible by \(y^tz^t\), but cannot have both exponents at least
\(2t\), since \(y^{2t}z^{2t}\in M^2\). They lie in two strips
with a coordinate in \([t,2t-1]\), which proves the count.

If their degrees are ordered as \(d_1\le\cdots\le d_q\), then
\(d_j\ge j/(2t)-1\). This follows from the count when \(d_j<a\);
otherwise use \(j\le q\le m\le2ta\). Consequently
\[
          \operatorname{cost}(U)\ge q(q+1)/(4t)-q.       \tag{8}
\]
For \(t\ge32\), (6)--(7) give
\[
 \frac{C_0-G}{t^5}
 \le\frac{32}{31}-\frac{841}{2048}+\frac1{2\cdot32^4}
 <\frac23.                                               \tag{9}
\]
If \(k(B)\le m/8\), then \(q\ge7m/8\ge441t^3/256\), while
\(q\le2t^3\). Equation (8) would instead give
\[
 \frac{\operatorname{cost}(U)}{t^5}
 \ge\frac{193969}{262144}>\frac7{10}.                    \tag{10}
\]
Upper-column costs are nonnegative. Thus (9)--(10) contradict
optimality. Every optimum at the balanced price has \(k(B)>m/8\).
A small generic positive perturbation preserves this strict exclusion
and supplies an actual initial cone. Combining with the large-price
endpoint proves (1), since \(m/8-t^2\ge m/16=n/32\).

## The sparse witness has exactly twelve cones

More generally take integers \(t\ge1\), \(2t<h<a-t\), and the same
sparse specialization (4). It has length \(4ta-2t^2\) and exactly
\[
                              g=12                       \tag{11}
\]
over every field. This is a count for the sparse family, not the open
locus in (1).

Here is a complete small encoding. Put \(\alpha=a-t\), and write
prices as \(X,Y,Z>0\). Apart from compulsory columns and additive
costs, its residue matrix splits into five types of blocks:

* \(t(h-t)\) copies each of the rank-one lists
  \((X,\alpha Y)\), \((X,\alpha Z)\);
* \(t(a-h-t)\) copies each of the rank-one lists
  \((X,\alpha Y,tZ-(h-t)Y)\) and its symmetric counterpart;
* \(t^2\) rank-two triangle blocks. Their two parallel ground lists
  cost \((X+hY,(\alpha+h)Y)\) and
  \((X+hZ,(\alpha+h)Z)\); the connecting edge costs \(t(Y+Z)\).

The negative coefficients in this encoding result from subtracting
one common cost per block. Its original monomial costs are positive.
Every nonzero maximal minor of a block is a unit, including in
characteristic two. Positive multiplicities preserve the common
refinement of these block fans.

Normalize \(Y=1\), \(Z=u\), \(X=v\). For \(u>1\), the triangle
chooses its connecting edge and the cheaper first ground list; the
other ground list is more expensive than both. The high \(z\)-arm
also has a fixed choice. The remaining low-arm walls are \(v=\alpha\)
and \(v=\alpha u\). In the bottom band \(0<v<\alpha\), the line
\(v=t(u+1)-h\) cuts it into two regions. In each of the middle and
top bands, the ray \(u=(\alpha+h)/t-1\) gives two regions.
All six are nonempty because \(\alpha>h>2t\). Symmetry gives six
more for \(u<1\). The triangle changes choice at \(u=1\), so the
two sides have distinct standard sets. There are no additional open
cones, proving (11).

The four upper counts on these cones are exactly
\[
                t^2,\quad th,\quad2th,\quad t(a+h-t).    \tag{12}
\]
Thus a large positive width does not itself supply a large fan.

## Verification and proof boundary

```bash
python3 -B scripts/run_hook_colour_width_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker compares the entire small block fan with independent
original-column border traversal in characteristic zero, two and 101,
and verifies original commuting multiplication. It separately checks
the replacement columns, the full-support candidate bound, exact degree
sums and rational margins. Large compressed cases verify the unit
column map and block encoding, not a full generic fan.

The residue-matroid and minimum-basis dictionary is the established
[Babson--Onn--Thomas framework](https://arxiv.org/pdf/math/0207135).
The unit generation argument, matching/triangle representation and
finite cost comparisons are elementary existing tools. The inspected
[Reissig, Theorem 3.2](https://arxiv.org/pdf/1305.4277) concerns whole
structured Toeplitz block truncations; it does not supply the selected
convolution-minor theorem needed for growing generic fan counts and
is not used to assume those ranks. Exact priority of this combined
width statement is unassessed.

The unrestricted count, uniform generic-hook exposure count and
nongeneric coefficient strata remain unresolved. Correlated profile
budgets remain possible; (1) only rules out reducing the colour-width
factor to sublinear order on this generic locus. The original extremal
objective and its completion criteria remain unchanged.
