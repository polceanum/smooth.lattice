# Near-quadratic fans with arbitrary scalar matching outer images

Let \(K\) be any field, \(A=a-t\ge t\ge2\), and let
\(F,G,H,J\in K[q]/q^A\) be arbitrary. Consider
\[
\begin{split}
 I=(\,x^2,\;&y^a-xy^tH(y),\ z^a-xz^tJ(z),\\
          &y^tz^t-x(y^tF(y)+z^tG(z))\,).
\end{split}                                                     \tag{SO1}
\]
These are compatible flat hooks of length \(n=4at-2t^2\).
The central arms and the matching outer images depend only on their
long variables; the opposite-arm outer images are zero. Nonunits,
zero operators, arbitrary degrees and every characteristic are allowed.
Then
\[
 \boxed{\quad g(I)=O\!\left(A^2t^2\log(2+A/t)\right)
                 =O\!\left(n^2\log(2+A/t)\right).\quad}           \tag{SO2}
\]
The same bound holds for a complete refinement of the positive fan.
There is a coefficient-dependent complete positive integer weight list
of this cardinality and height \(O(a^2)\).

This extends the [quadratic scalar-axis theorem](hook_axis_lattice_windows.md)
to arbitrary matching scalar outer images. It does not assume a small
polynomial or rational degree, regular minors, unit axes or a fixed
thickness. Transverse coefficients, nonzero opposite-arm outer images,
the sharp unrestricted \(G_3(n)\) and an infinite polynomial excess
remain unresolved. Exact strongest-current novelty is not established
by the bounded primary comparison below.

## Three pure streams retain the full two-factor rank formula

The [original hook reduction](hook_bridge_reduction.md) applies to
(SO1). Its compatibility conditions hold because every cross-multiple
of the displayed scalar matching images vanishes in the planar hook
diagram. The two copies of that diagram give the stated length.
After the fixed diagram and upper corner are removed, the Y row space is
\(K[q]/q^A\otimes K^t\), with pure sources
\[
\begin{array}{c|c|c}
 \text{family}&\text{long vector}&\text{original label}\\ \hline
 \text{upper}&q^r&(1,t+r,j),\quad 0\le r<A\\
 \text{central}&q^rF&(0,t+r,t+j),\quad t\le r<A\\
 \text{outer}&q^rH&(0,a+r,j),\quad 0\le r<A .
\end{array}                                                     \tag{SO3}
\]
Every vector is tensored with the short coordinate \(e_j\),
\(0\le j<t\). Transpose the last two label coordinates for the
Z arm, replacing \((F,H)\) by \((G,J)\). The only coupled columns are
the central heads
\[
       (q^iF\otimes e_j,\ e_i\otimes q^jG),
             \qquad \text{label }(0,t+i,t+j),\quad i,j<t.
\]
All outer columns are pure in this class.

For \(w=(x,y,z)>0\), define
\[
\begin{array}{ll}
 \theta=(x-tz)/y,& \eta=A-tz/y,\\
 \phi=(x-ty)/z,& \zeta=A-ty/z .
\end{array}                                                     \tag{SO4}
\]
In the Y factor, prices are \(\theta+r,r,\eta+r\) for the
three families. Define the actual birth \(\alpha_i(\theta,\eta)\)
of \(q^iF\) modulo earlier heads using the whole pure price prefix.
Absent heads have birth \(-\infty\). Use \(\beta_j(\phi,\zeta)\)
for Z. The two-flag tensor-prefix proof permits an arbitrary pure
source family in each factor, so its unchanged exact formula is
\[
 P_Y(w)+P_Z(w)-
 \sum_{i,j<t}\max\{y(\alpha_i-i),z(\beta_j-j),0\},                 \tag{SO5}
\]
plus the fixed linear translation
\((t^2,at(a+t-2)/2,at(a+t-2)/2)\).
Here the two \(P\)'s are the full pure minimum-basis costs from (SO3).
See the [prefix-rank proof](hook_axis_birth_curves.md); the new
outer family is included among its pure sources, not discarded.

Price orders change only when one of
\[
       \theta,\eta,\theta-\eta,\phi,\zeta,\phi-\zeta
                                                        \tag{SO6}
\]
crosses an integer in a range of size \(O(A)\).
These are integer planes in the original weights. Zero-gain boundaries
also lie on the integer \(\theta,\eta,\phi,\zeta\) pencils.
A central birth has label \(r\ge t>i\), so its gain cannot vanish.

## A bounded-shift lemma for polynomial lattices

Let \(\mathcal L\subset K[q]^m\) be a full-rank polynomial lattice.
Fix shifts \(\sigma^0\) with pairwise distinct fractional parts and
a weak shifted Popov basis \(g_1,\ldots,g_m\). Write its weighted
row degrees as \(d_1,\ldots,d_m\). Distinct leading positions give
the classical predictable-degree identity
\[
  \deg_{\sigma^0}\left(\sum f_jg_j\right)
                   =\max_j(d_j+\deg f_j).                       \tag{SO7}
\]
The largest summand has a leading position the others cannot cancel.
Distinct fractional parts avoid ties; ordinary shifted reduced bases
give the same identity more generally.

Impose an affine constraint on the vector modulo \(q^t\), with a
nonempty feasible set whose prescribed residues exclude zero.
Let \(\gamma_\sigma\) be the minimum \(\sigma\)-degree of a
lattice vector satisfying it. Suppose
\(\|\sigma-\sigma^0\|_\infty\le E\). Then
\[
       \gamma_\sigma\in
          \bigcup_{j=1}^m[d_j-E,\ d_j+t-1+E].                    \tag{SO8}
\]
Indeed choose a minimizing \(v=\sum f_jg_j\), and replace each
\(f_j\) by \(f_j\bmod q^t\), obtaining \(\bar v\). Every component
modulo \(q^t\) is unchanged. By (SO7),
\(\deg_{\sigma^0}\bar v\le\deg_{\sigma^0}v\), and the former belongs
to \(\{d_j+k:0\le k<t\}\). Shift perturbations change the degree of
every fixed vector by at most \(E\). Minimality therefore gives
\[
 \deg_{\sigma^0}\bar v-E
   \le\deg_{\sigma^0}v-E
   \le\gamma_\sigma
   \le\deg_\sigma\bar v
   \le\deg_{\sigma^0}\bar v+E.
\]
This proves (SO8). The argument preserves all dependencies and needs
no monotonicity of a multi-parameter pivot table.

For the three-stream factor take
\[
 \mathcal L_{F,H}=\{(u,b,c):u-Fb-Hc\in q^AK[q]\},
 \qquad \sigma=(\theta,0,\eta).                                 \tag{SO9}
\]
Its basis \((q^A,0,0),(F,1,0),(H,0,1)\) has determinant \(q^A\).
Impose \(b\bmod q^t=q^i+\sum_{j<i}c_jq^j\). The exact minimum is
\[
                   \gamma_i=\max(i,\alpha_i).                  \tag{SO10}
\]
The low coefficients of \(b\) describe the head and earlier heads,
its higher coefficients describe central sources, and \(c\) describes
outer sources. The weighted bound on \(u\) describes upper sources.
Reducing all three polynomials modulo \(q^A\) preserves the constraint
and congruence and never increases a degree. This proves both
directions of (SO10), including absent heads.

Consequently, on any parameter box in which both \(\theta,\eta\)
range over intervals of width at most \(t\), (SO8) with \(m=3,E=t\)
places every positive birth in three intervals of width \(3t\).
Subtract the corresponding source shift and the head index \(i<t\).
There are \(O(t)\) possible integer labels \((f,r-i)\) for all
positive births on the entire box, where \(f\) specifies the family.
One fixed base point with fractional shifts \(1/3,0,2/3\) suffices.
The same statement applies to \((\phi,\zeta)\).

The positive Y gain gradients in original weight coordinates are
\[
 (1,e,-t),\qquad (0,e,0),\qquad (0,A+e,-t),\quad e=r-i,           \tag{SO11}
\]
and the Z gradients transpose their last two coordinates.
Thus \(O(t)\) possible labels on both arms give \(O(t^2)\) comparison
planes for (SO5), including comparisons with zero.

## When all uppers come last, only a two-stream kernel is needed

Suppose a factor has shifts \(T,S\), with \(T>2A\) and \(S\le A\).
Every upper is then more expensive than every central and outer
source. Put
\(\nu_F=\operatorname{ord}_qF\), \(\nu_H=\operatorname{ord}_qH\),
with the order of zero interpreted as \(A\). The full lower-source
span is
\[
        q^mK[q]/q^A,\qquad m=\min(A,t+\nu_F,\nu_H).              \tag{SO12}
\]
Modulo this span the remaining nonzero heads are supported on the
upper interval
\(\nu_F\le r<\min(A,\nu_F+t)\). Therefore every upper birth
uses one of \(O(t)\) labels, independently of \(T\).

Every other finite birth occurs before all uppers, and is a
lower-only birth. It is governed by
\[
        \mathcal K_{F,H}=\{(b,c):Fb+Hc\in q^AK[q]\},
                         \qquad \sigma=(0,S).                  \tag{SO13}
\]
If \(\nu=\min(\nu_F,\nu_H)<A\), factor out \(q^\nu\).
One of the remaining factors is a unit modulo \(q^{A-\nu}\);
truncated inversion gives a free rank-two relation lattice of
determinant degree \(A-\nu\). Apply (SO8) to the same low-\(b\)
constraint. If \(A-\nu<t\), use that smaller cutoff; the remaining
head suffix is zero. Thus when \(S\) ranges over width \(t\), all
lower-only positive births also have \(O(t)\) labels.
If both operators are zero there are no nonzero heads.
These statements include the pure-basis changes: before the uppers,
orders depend only on integer \(S\), while the upper complement to
(SO12) is a fixed coordinate prefix.

## Partition the positive weights without paying a third large grid

By transposition take \(y\ge z\), normalize \(w=(s,1,\rho)\), and put
\[
 0<\rho\le1,\quad
 \theta=s-t\rho,\quad\eta=A-t\rho,\quad
 \phi=(s-t)/\rho,\quad\zeta=A-t/\rho,\quad
                         \rho_0=t/(4A).                        \tag{SO14}
\]
In particular \(\theta>-t\) and \(\eta\in[A-t,A]\).
An interval for any displayed ratio pulls back to two linear
inequalities in \(w\); all regions below are convex projective
polygons. All finite interval grids have width at most \(t\).

**Very unbalanced weights, \(\rho<\rho_0\).**
If \(\phi<-A\), every Z upper has negative factor price and all
Z gains are negative. Formula (SO5) has only Y-versus-zero maxima.
The \(O(A)\) global price and sign pencils (SO6) refine this entire
region, giving \(O(A^2)\) cells.

Otherwise \(\phi\ge-A\), so \(s\ge t-A\rho\). Every Z outer
precedes every upper and central source: its largest price is
\(\zeta+A-1<-2A-1\), and
\[
 \phi-(\zeta+A-1)=s/\rho-2A+1>A+1.
\]
Contract their full span \(q^{\nu_J}K[q]/q^A\).
Their minimum basis is a fixed outer prefix; the remaining factor is
exactly the scalar upper/central factor for
\(G\bmod q^{\nu_J}\) at cutoff \(\nu_J\), with absent heads omitted.
This contraction is used after proving those columns strictly cheaper.
It preserves all positive births. If \(J=0\), the cutoff is \(A\).

Grid \(\phi\in[-A,2A]\), adding its upper tail. Grid
\(\theta\in[-t,2A]\), adding its upper tail.
Each has \(O(A/t)\) intervals. The Y shift \(\eta\) already varies
over width \(t\). Use the rank-three lemma for a finite Y interval,
(SO13) for its tail, and the rank-two scalar lemma for Z. Each
coarse interval pair has \(O(t)\) positive labels on each arm.
There are \(O((A/t)^2)\) such pairs.

**Uppers last in Z, \(\rho\ge\rho_0,\ \phi>2A\).**
Here \(-3A\le\zeta\le A\). Grid \(\zeta\) and grid \(\theta\)
up to \(2A\), adding the Y upper tail. Again there are only
\(O((A/t)^2)\) pairs. Formula (SO13) supplies Z labels, and
(SO9) or (SO13) supplies Y labels. Both arms have \(O(t)\)
positive labels per coarse pair.

**The remaining region, \(\rho\ge\rho_0,\ \phi\le2A\).**
Now \(-4A<\phi\le2A\). Split \(\rho\) dyadically into
\([r,2r]\), truncated at one, starting at \(\rho_0\).
There are \(O(\log(2+A/t))\) bands. On one such band,
\[
 \#\{\zeta\text{-bins}\}=O(1+1/r),\qquad
 \#\{\phi\text{-bins}\}=O(A/t).                                 \tag{SO15}
\]
Also \(\theta=t+(\phi-t)\rho\), so
\[
 -t<\theta\le t+4Ar,\qquad
              \#\{\theta\text{-bins}\}=O(1+Ar/t).
\]
The number of possible triples of bins in this band is at most
\[
 O\!\left((1+1/r)(A/t)(1+Ar/t)\right)
                          =O((A/t)^2),                         \tag{SO16}
\]
because \(t/(4A)\le r\le1\) and \(A\ge t\).
Each triple pulls back to one convex polygon, possibly empty.
The Y pair \((\theta,\eta)\) and Z pair \((\phi,\zeta)\)
each vary by at most \(t\), so (SO8) again gives \(O(t)\)
positive labels on each arm. Counting all three grids without the
dyadic restriction would incorrectly introduce an extra large factor.

## Count full arrangements and recover integer test weights

On every coarse polygon outside the first inactive region, include
the comparisons in (SO11), and the integer pencils (SO6) crossing
its factor intervals. There are \(O(t^2)\) comparison planes and
only \(O(t)\) additional pencils. For Y the difference
\(\theta-\eta\) ranges over width at most \(2t\); the same holds for
Z when both intervals are finite. In an upper-last tail, only the
lower-order pencil remains. In the contracted extreme region, Z
has only the upper/central pencil.

These planes refine every term of (SO5). Pure bases and head labels
are fixed between price pencils; zero crossings are on the indicated
integer pencils; every remaining change compares two positive gains.
Thus their complete line arrangement has \(O(t^4)\) cells in a coarse
polygon. All contacts, intersections, parallel pieces and repeated
wall portions are retained. Summing (SO16) over the dyadic bands
and adding the other regions gives
\[
 O(A^2)+O((A/t)^2)\,O(t^4)\,O(\log(2+A/t)),
\]
which proves (SO2), including the complete refinement bound.
The plane \(y=z\) and transposition only change the constant.

Every normal has integer height \(O(a)\). Integer bin endpoints
have size \(O(A)\); dyadic boundaries
\(\rho=2^kt/(4A)\) also have integral normals of height \(O(A)\).
Primitive cross products give extreme rays of height \(O(a^2)\).
The sum of three independent primitive nonnegative extreme rays is
a positive integer interior weight of that same order.
One per full-dimensional refinement cell gives the complete list
claimed above. It depends on the actual lattice bases and hence on
the coefficients. No optimal basis-output or bit-complexity claim
is made.

## Verification and primary-source boundary

~~~bash
python3 -B scripts/run_hook_scalar_outer_lattice_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker constructs rank-three relation lattices independently of
finite source-prefix elimination, checks exact determinants and
congruences, compares short-lattice and original birth costs, and
perturbs both shifts independently. It separately checks the
two-stream lower kernel, including nonunits and infeasible constraints.
Complete rational fan subdivisions cover projective area one.
Every interior facet belongs to its complete local inventory,
small integer weights reproduce the profiles, and convexly merged
profile cones agree with independently enumerated original quotient
fans. Rational and finite-field fixtures include nonzero outers,
zero outers, absent heads and nonunit axes.
Larger exact factor queries at \((A,t)=(96,2),(64,3)\) check
nonvacuous label exclusions; they are not complete original-fan
enumerations. Separate dyadic grid products check (SO16) across
growing aspect ratios.

The classical inputs are polynomial relation modules and shifted
predictable degree:
[Labahn--Neiger--Vu--Zhou, Sections 2.2 and 2.5--2.6](https://arxiv.org/pdf/2202.09329),
and leading-position lattice reduction:
[Mulders--Storjohann, Section 2](https://cs.uwaterloo.ca/~astorjoh/jscpopov.pdf).
The bounded-shift argument is explicitly derived from those properties.
The inspected statements do not themselves give (SO2), its dyadic
weight partition or the stated complete weight list. This is a bounded
comparison, not a strongest-current novelty certificate.
The previously proved restricted family results and the original
unrestricted extremal and publication requirements remain in force.
