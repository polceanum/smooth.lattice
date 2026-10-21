# Hurwitz rank tables for every quotient of a fixed quartic curve

Let \(K\) be any field and \(F\in K[t]\) be monic of degree \(n>0\),
with \(F(0)F(1)\ne0\). Put
\[
 I_F=(xy-1,(x-1)^2z-1,F(x))\subset K[x,y,z].                 \tag{1}
\]
Uniformly in \(F\), its multiplicities and the characteristic,
\[
 \begin{split}
 \Omega(n^2)&\le g(I_F)\le O(n^2\log(n+1)),\\
 \Omega(n^3)&\le V(I_F)\le O(n^3\log(n+1)),\\
 \Omega(n^2)&\le|\mathcal U(I_F)|\le O(n^2\log(n+1)).
 \end{split}                                                \tag{2}
\]
Here \(V\) sums primitive costs over all actual undirected wall pieces,
with multiplicity, and \(\mathcal U\) is the union of reduced Gröbner
bases. The constants are absolute; lower bounds apply for sufficiently
large \(n\). The quotient length is \(n\).

These are all finite quotients of the one fixed curve: its ring is
\(K[t,t^{-1},(t-1)^{-1}]\), a localization of a principal ideal domain.
Every nonzero finite-colength ideal has a unique monic generator avoiding
\(0,1\). No generic evaluation, reducedness, splitting or characteristic-zero
hypothesis is used.
The projective closure has degree four: the birational parametrization
\((1,t,t^{-1},(t-1)^{-2})\) has common pole divisor
\((\infty)+(0)+2(1)\), of degree four, with no residual base point.

This improves the general fixed-curve \(O_C(n^{14/5})\) envelope on this
class and extends the [generic-sample comparison](pole_order_fans.md) to
every quotient of this quartic carrier. It excludes this entire class as a
polynomial-separation route above the generic corner-cut scale. It does
not resolve unrestricted \(G_3\), remove the logarithm, or establish a
publication result. Hurwitz matrices are classical; exact priority of this
quantitative application remains unassessed.

## The three-ray case was already elementary

For a fixed monomial carrier decomposed into \(t\ge1\) disjoint
divisibility rays and a finite set \(E\), every standard downset is
determined by its subset of \(E\) and its ray lengths. Their sum is
\(n-|D\cap E|\). The number of possibilities is at most
\[
                   2^{|E|}\binom{n+t-1}{t-1}.              \tag{3}
\]
Sum over the finitely many curve carriers. With at most three rays this
is already \(O_C(n^2)\); with no rays only finitely many lengths occur.
A proposed three-list face argument is unnecessary. The curve (1) has
four rays in both carriers, so this count alone gives cubic order.

## Four lengths and an exact all-field rank formula

For positive costs \((\alpha,\beta,\gamma)\), its curve initials are
\[
 \begin{cases}
 (xy,x^2z,yz),&\alpha<\beta,\\
 (xy,xz,y^2z),&\beta<\alpha.
 \end{cases}                                                \tag{4}
\]
Use \(xy-1,\ x^2z-2xz+z-1,\ yz+xz-2z-y\) in the first chart.
In the second use \(xy-1,\ yz+xz-2z-y,\ y^2z-2yz+z-y^2\).
Their S-pairs reduce to zero integrally. Principal parts at
\(0,1,\infty\) also verify independence of the proposed carriers.
Thus (4) remains valid when \(2=0\) in \(K\).

Every initial standard set in one chart has the form
\[
 D^v_{a,b,c,d}=\{1,x,\ldots,x^a,y,\ldots,y^b,z,\ldots,z^c,
                                      vz,\ldots,vz^d\},      \tag{5}
\]
where \(v=x\) in the first chart, \(v=y\) in the second. All lengths
are nonnegative, \(c\ge d\), and \(d>0\) requires \(a\ge1\) for
\(v=x\), or \(b\ge1\) for \(v=y\). Put
\[
 m=1+a+b+c+d,\quad h=c-d,\quad \ell=\max(0,m+h-n).
\]
Expand \(G_b(s)=F(1+s)/(1+s)^b=\sum_{j\ge0}g_j(b)s^j\), setting
negative-index coefficients to zero, and define
\[
 H_{h,\ell}(b)=\bigl(g_{2i+1-j}(b)\bigr)_{
                        0\le i<h,\ 0\le j<\ell}.           \tag{6}
\]
**Rank theorem.** In either chart, the residue columns of (5) have rank
\[
                 m-\ell+\operatorname{rank}_K H_{h,\ell}(b).
                                                               \tag{7}
\]
In particular, a size-\(n\) allocation is a basis precisely when
\[
 \Delta_h(b):=\det H_{h,h}(b)\ne0,\qquad \Delta_0(b)=1.       \tag{8}
\]
Its feasibility depends only on \(b,h\), not separately on \(a,c,d\).

**Proof.** Put \(s=t-1\). The rational span of (5) equals
\[
 \langle1,t,\ldots,t^a,t^{-1},\ldots,t^{-b},
 s^{-1},\ldots,s^{-2d},s^{-2d-2},s^{-2d-4},\ldots,s^{-2c}\rangle.
                                                               \tag{9}
\]
For the \(x\) tag, use \(tz^j=s^{-2j}+s^{-2j+1}\). For the \(y\)
tag, partial fractions of \(t^{-1}s^{-2j}\) have top terms
\(s^{-2j}-s^{-2j+1}\), lower principal parts, and a multiple of
\(t^{-1}\). Induction for \(j\le d\) proves (9). The downset condition
supplies the needed core \(x\) or \(y\); the top coefficients are units
in every characteristic.

Multiplication by the unit \(t^bs^{2c}\) modulo \(F\) identifies (9)
with polynomials \(p(t)\) of degree at most \(m+h-1\) satisfying
\[
 [s^{2i+1}]\frac{p(1+s)}{(1+s)^b}=0,\qquad 0\le i<h.        \tag{10}
\]
These are exactly the missing odd principal parts at \(1\). Partial
fraction uniqueness proves the converse characterization and dimension
\(m\). The kernel modulo \(F\) is \(p=FQ\), with \(\deg Q<\ell\).
Expanding \(Q\) in powers of \(s\), (10) has matrix (6). Its kernel
has dimension \(\ell-\operatorname{rank}H\), proving (7). No root
separation or division by a derivative is involved. \(\square\)

## Feasibility transport gives a hyperbolic candidate budget

Take an actual size-\(n\) standard set (5), inside a strict chart.
Normalize \(\gamma=1\), set \(t=\alpha/\gamma,\ u=\beta/\gamma\),
\(r=\min(t,u)\), and \(M=n-1-b\). The tag costs \(r\). Up to a term
independent of \(a,h\), its basis cost is
\[
 \Phi(a,h)=\frac{2t+1}{4}a^2+\frac{t-r-M-1}{2}a+
                                      \frac{h^2}{4}-\frac r2h.
                                                               \tag{11}
\]
These divisions are real arithmetic on exponents, not arithmetic in \(K\).

Suppose \(a\ge3\), excluding \(v=y,b=0\); that exception has \(d=0\)
and only \(n\) allocations in total. The shift
\[
                   (a,b,c,d)\longmapsto(a-2,b,c+1,d+1)     \tag{12}
\]
preserves the downset conditions and \(b,h,n\), hence preserves being a
basis by (8). Comparing costs gives
\[
 (2t+1)a\le M+2+t+r,\qquad
                         2t(a-1)\le M+2-a\le n+1.         \tag{13}
\]
Every size-\(n\) allocation with \(h=0\) is a basis. If \(h\) is even,
compare with \((a,0)\) in the coordinates of (11). If \(h\) is odd,
compare with \((a+1,0)\); it is permitted since \(a+h\le M\).
The increase of the \(a\)-part in the latter comparison is at most
\(3t/2+3/4\), by (13). In both cases,
\[
 h^2-2rh\le6t+3,\qquad h<2t+3,\qquad (a-1)h<4n.           \tag{14}
\]
For the middle implication use \(r\le t\) and
\((2t+3)^2-2t(2t+3)>6t+3\). The last uses (13) and \(a\le n-1\).

For each \(b<n\), the cases \(a\le2\) contribute \(O(n)\) allocations.
All others satisfy
\[
 3\le a\le n-1-b,\quad
 0\le h\le\min(n-1-b-a,\lfloor4n/(a-1)\rfloor),\quad
                    h\equiv n-1-b-a\pmod2.                \tag{15}
\]
The harmonic sum bounds these pairs by \(O(n\log(n+1))\). Summing over
\(b\) and both charts proves the upper fan bound in (2). Every open
cone meets a strict chart; the line \(\alpha=\beta\) adds no open
initial ideals. This is a count of actual candidate standard sets,
with all joint dependencies retained through (8), rather than a count
of wall carriers that ignores their intersections.

There is a sharper diagnosis: the entire chart \(\alpha<\beta\)
has \(\Theta(n^2)\) initials, uniformly in \(F\). To prove its upper
bound fix \(b\) and the parity of \(h\). Then \(a\) has fixed parity,
feasible \(h\)'s form one fixed set by (8), and (11) has \(r=t\).
Its only parameter is \(t\); the \(u\)-term is constant. As \(t\)
increases through unique minima, \(a\) cannot increase and \(h\)
cannot decrease. Here is an elementary uncrossing proof.

For two minima at \(t_1<t_2\), if both coordinates move oppositely
to the asserted directions, the slope \(a^2/2-h/2\) increases,
contradicting optimality. If they move in the same coordinate direction,
swap their \(a\) coordinates. Both cross-pairs are feasible because
\(a+h\le M\), the parities agree, and membership in the feasible
\(h\)-set is unchanged. The only permitted \(a=0\) point is \(h=M\);
it causes no exception to those swaps. Separability of (11) now
contradicts optimality in the coordinate moving in the wrong direction.
Thus at most \(2M+1\) distinct pairs appear per parity. Sum over \(b\).
The lower witnesses below are all in this chart. The unresolved logarithm
in (2) is confined to the other chart; an arbitrary hole set at each
fixed \(b\) is insufficient to remove it.

## Quadratic output and cubic wall cost are forced over every field

Choose integer pairs
\[
 2n/5\le a\le n/2,\quad n/10\le b\le n/5,\quad
 a+b\equiv n-1\pmod2,\qquad c=(n-1-a-b)/2.                 \tag{16}
\]
For sufficiently large \(n\), there are \(\Theta(n^2)\) such pairs,
\(0<c<a-1\) and \(b<a-1\). At
\[
 \gamma=1,\qquad\alpha=\frac{c}{a-1/2},\qquad
                         \beta=\frac{\alpha+c}{b+1},     \tag{17}
\]
we have \(\alpha<\beta\). Exactly \(n-1\) carrier monomials cost less
than \(\alpha+c\): they form \(D^x_{a,b,c,c}\setminus\{xz^c\}\).
Exactly \(xz^c,y^{b+1}\) tie at that threshold. All other carrier
monomials cost more, as follows from
\[
 a\alpha<\alpha+c<(a+1)\alpha,\quad
 c<\alpha+c<c+1,\quad b\alpha<c.
\]
The \(h=0\) allocation is always independent. Increasing \(\beta\)
slightly exposes it strictly, proving the first lower bound in (2).

The other allocation replaces \(xz^c\) by \(y^{b+1}\) and has
\(b'=b+1,h'=1\). It is independent exactly when
\[
                       F'(1)-(b+1)F(1)\ne0.               \tag{18}
\]
When (18) holds, (17) lies inside an actual wall with primitive direction
\[
 q=(1,-(b+1),c),\qquad \nu(q)=2c+b+2=n-a+1=\Theta(n).      \tag{19}
\]
Different pairs give different directions. In characteristic other than
two, at most one of two consecutive admissible \(b\)'s, differing by
two, fails (18). In characteristic two, pairing consecutive \(a\)'s
switches admissible \(b\)-parity; at least one row in each pair has no
failures. If \(F'(1)/F(1)\) is outside the prime field there are none.
Thus \(\Theta(n^2)\) walls of cost \(\Theta(n)\) remain uniformly over
all fields, proving the cubic lower cost.

On the \(h=0\) side the reduced polynomial headed by \(y^{b+1}\)
has nonzero \(xz^c\) coefficient, since that exact exchange is possible.
For fixed \(b\), its largest tagged \(z\) exponent is therefore \(c\);
larger tagged powers are absent from that basis. Different \(c\)'s
give different reduced polynomials, and different \(b\)'s give different
heads. This proves the quadratic lower reduced-basis union.

Each initial ideal is generated by the three fixed carrier heads and
at most four ray frontiers, so each reduced basis has at most seven
polynomials. The positive three-variable adjacency graph is planar and
simple. Every wall has \(\nu(q)\le n\), by the
[primitive threshold theorem](primitive_wall_spectrum.md). Consequently
\(V\le3ng\) and \(|\mathcal U|\le7g\), proving the upper bounds in (2).
These are output-cardinality and primitive-cost bounds, excluding bits.

## Classical inputs, enumeration and the remaining question

The transpose of (6) at \(\ell=h\) is the leading Hurwitz block of
the ascending coefficient sequence \(g_j(b)\); reversing a finite
polynomial gives the descending coefficient convention.
[Holtz--Tyaglov, Sections 1.1--1.3](https://arxiv.org/pdf/0912.4703)
give the established structured-matrix, resultant and continued-fraction
framework. Its complex root-location or positivity conclusions are not
assumed over arbitrary fields. The coefficient pattern is also the
\(M=2\) case of [Holtz--Khrushchev--Kushel, Section 3](https://arxiv.org/pdf/1506.07379).
The matrix type, Euclidean algorithm and minor identities are prior
tools; no exact fan-count priority conclusion follows from this comparison.
More precisely, formula (1.37) and Definition 1.16 in Holtz--Tyaglov
give exactly the square coefficient pattern in (8). Their Sections
1.1 and 1.3 supply the surrounding rational-function and Euclidean
framework, rather than an independently verified fan-count theorem.
The introduction of
[Conca--De Negri--Rossi](https://msp.org/ant/2007/1-3/ant-v1-n3-p01-p.pdf)
counts Cohen--Macaulay initials of the homogeneous rational normal curve,
with its embedding degree varying. That stated parameter and object
differ from finite quotients of the fixed affine curve (1). Only its
introduction and definitions were compared here; a complete theorem-level
priority audit, including possible reductions between the problems, remains
required before any novelty claim.

For an exact enumeration reduction, generate (15) and the small and
exceptional cases, retain precisely those passing (8), and form their
profile hull in each chart. Restrict each minimum fan to its chart.
Every retained allocation is a genuine quotient basis and every coherent
basis is retained, so this recovers the full fan. There are
\(O(n^2\log(n+1))\) candidates before any field-dependent filtering.

All tests (8) come from leading principal ranks of \(n\) matrices of
size at most \(n\). The classical rank-profile algorithms of
[Dumas--Pernet--Sultan](https://arxiv.org/pdf/1601.01798)
give \(O(n^{\omega+1})\) field operations, for any feasible dense
matrix exponent \(\omega>2\), plus exact hull arithmetic.
Consecutive coefficient series satisfy \(G_{b+1}=G_b/(1+s)\), so
preparation through degree \(2n\) costs \(O(n^2)\) field operations.
No new matrix algorithm, optimized implementation or bit bound is claimed.

The next question is whether the joint Hurwitz table removes the logarithm,
or admits coherent examples of that order. Higher pole orders require
their own joint minor analysis; generic rank truncation cannot be assumed.
The unrestricted collective budget and original \(G_3\) goal remain open.

## Exact controls

~~~bash
python3 -B scripts/run_quartic_hurwitz_checks.py
~~~

The report is results/research/quartic_hurwitz_checks.json. Independent
polynomial-remainder columns check (7), including rectangular blocks
and allocation sizes above and below \(n\). Whole border-cone fans are
compared on complete chart patches with minima of all feasible allocation
profiles. The shift, candidate inventory and lower-wall witnesses are
checked separately, in characteristic zero and finite fields, with
repeated factors and exhaustive small binary moduli. These controls
support the proofs, not current novelty or the unrestricted endpoint.
