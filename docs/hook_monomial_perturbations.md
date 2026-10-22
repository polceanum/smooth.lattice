# Monomial perturbations localize the valuation-kernel exceptions

A large valuation loss in an all-vector uncertainty estimate need not
survive in the fan. For central-only hooks whose two off-axis perturbations
are monomials times units, primal kernel exceptions are ineligible and dual
kernel exceptions force an exchanged exponent into a short translated
strip. This removes the large-valuation factor from the fan bound, uniformly
in both valuations and in the remaining unit coefficients.

This is a restricted-family closure. The sharp unrestricted \(G_3\) question,
an actual polynomial lower separation, exact priority and the original
final-paper criteria remain open.

## Family and uniform bound

Let \(t\ge2\), \(a\ge2t\), \(A=a-t\), \(L=a-2t\), and
\[
 I=(x^2,y^a,z^a,y^tz^t-xf),\qquad n=4at-2t^2.             \tag{MP1}
\]
The two outer generator images are zero. Write \(f=y^tU+z^tV\) in
the rectangular arm modules
\[
 R_Y=K[y,z]/(y^A,z^t),\qquad R_Z=K[y,z]/(y^t,z^A).
\]
Assume each nonscalar arm has the form
\[
 U=\alpha+y^{\nu_Y}z^{d_Y}B_Y,\qquad
 V=\beta+z^{\nu_Z}y^{d_Z}B_Z,                             \tag{MP2}
\]
where \(\alpha,\beta\ne0\), \(B_Y,B_Z\) are units, \(0\le\nu_i<A\),
and \(1\le d_i<t\). All coefficients of these units are allowed.
This includes arbitrary higher transverse coefficients divisible by the
same long-variable power. It does not assume that a first transverse
coefficient has nonzero constant term.

Put \(G_Y(y)=B_Y(y,0)\), \(G_Z(z)=B_Z(0,z)\). Let \(H_i\) be the
maximum unit Toeplitz prefix gap of \(G_i\) at the positive lengths
\(A-\nu_i,L-\nu_i\), omitting nonpositive lengths. The gap convention
and regular determinants are those of (TP3), with endpoints \(0,s\).
Define
\[
                         \Theta=3t+H_Y+H_Z.              \tag{MP3}
\]

**Theorem.** Over every field and uniformly in these coefficients,
\[
 D_{\rm all}(I)=O(a\Theta(1+\log(a+1))),\qquad
 g(I)=O(na\Theta(1+\log(n+1))^2).                         \tag{MP4}
\]
If the residual unit prefixes are regular, then
\[
                         g(I)=O(n^2\log^2(n+1))          \tag{MP5}
\]
for all valuations, including those arbitrarily close to \(A\).
If either arm is scalar, the earlier direct-sum bound applies without
conditions on the other arm. Retain the older box ceiling when smaller.

For fixed \(t,a,\nu_Y,\nu_Z\) and an upper bound \(\Theta_0\), a common
coefficient-independent complete list has
\[
 O(a^2\Theta_0^2(1+\log(a+1))^2)
\]
positive integer weights of height \(O(a^2)\). The valuation parameters
are fixed for this common list; no common list of this size for all
valuation choices at once is asserted.

The regular locus is nonempty in every characteristic over an infinite
field or a sufficiently large finite extension. Explicit axis units are
\[
 G_i(q)=(1+q)^{A-\nu_i-1}
\]
in characteristic zero or \(p>2A\), by the classical binomial Toeplitz
product. Higher transverse unit coefficients remain arbitrary.

## Strong uncertainty away from a monomial annihilator

Use long and transverse variables \(q,v\), and write
\[
 R_m=K[q,v]/(q^m,v^t),\qquad
 N=q^\nu v^d B,\qquad W=\gamma+N,\quad\gamma\ne0,\ B\text{ a unit}.
\]
When \(\nu<m\),
\[
 J_m=\operatorname{ann}_{R_m}(N)
                  =(q^{m-\nu},v^{t-d})R_m.               \tag{MP6}
\]
When \(\nu\ge m\), \(N=0\) and \(J_m=R_m\). These are actual
monomial subspaces because multiplication by \(B\) is invertible.

Suppose \(2\nu\ge m\) and \(H\notin J_m\). Put \(s=m-\nu>0\) and
let \(h\) be the image of \(H\) in
\[
                R_m/J_m=K[q,v]/(q^s,v^{t-d}).
\]
On this quotient \(W\) acts as the scalar \(\gamma\), since \(\nu\ge s\).
Thus \(\deg_qh\le\min(\deg_qH,\deg_qWH)\), while
\[
 NH=q^\nu v^d(Bh),\qquad
 \nu+\deg_q(Bh)\le\max(\deg_qH,\deg_qWH).
\]
The unit \(B\), whose transverse axis is \(G(q)\), satisfies
\[
                   \deg_qh+\deg_q(Bh)\ge s-1-h_s(G).     \tag{MP7}
\]
Indeed each regular univariate prefix minor becomes a block triangular
matrix on the transverse module, with identical invertible diagonal
blocks. Moving to the next regular index gives the stated gap loss,
by the same kernel contradiction as (TP7).
Consequently
\[
       \deg_qH+\deg_q(WH)\ge m-1-h_s(G),\qquad H\notin J_m.
                                                               \tag{MP8}
\]
The maximum-weight version follows from nonnegative exponents.
For \(H\in J_m\), \(WH=\gamma H\); no degree uncertainty is claimed.

Apply this at \(m=A,L\) when \(\nu\ge A/2\).
When \(\nu<A/2\), the previous valuation theorem already has
\[
 \max(e_A(q^\nu G),e_L(q^\nu G))\le H_i+t.
                                                               \tag{MP9}
\]
Thus small valuations incur only an \(O(t)\) allowance, and it remains
to treat the exact annihilator exceptions for the large valuations.

## Primal annihilator exceptions cannot exchange across arms

Use the original central-source and upper-identity columns, with no
outer images. At an eligible pure cross-arm pair, the source relation
and cost maxima are the ones in (TP12)--(TP13), or their same-layer
versions. If its \(Y\)-source restriction lies in \(J_A\), then
\(UH_Y=\alpha H_Y\). Because \(J_A\) is monomial, the leading
monomial of that restriction is also in \(J_A\).

If the exchanged \(Y\)-label is a lower source, the strict cost relation
forces \(\lambda<t\); that source is parallel to its cheaper upper
identity. If it is an upper identity, the relation forces \(\lambda>t\);
it is parallel to its cheaper pure central source. Either makes the
column a loop modulo cheaper columns. The same argument on the \(Z\)
arm uses the threshold \(tu\). All long exchanged labels have the
matching pure source available.

Hence actual long exchanges exclude \(J_A\) on both large-valuation
arms. The primal lower estimates (TP14) hold with losses at most
\(H_i+t\), for all three exchange types. This excludes the sharp
multiplier-kernel degree witnesses from the primal wall relations.
It does not assume that every vector outside the scalar transverse
tail has a large degree sum.

## Dual exceptions lie in a strip of width \(t\)

Outside the known \(O(at)\) directions and the bounded-product inventory
\(|bc|\le64a\Theta\), the cheap-overlap argument from (TP10)--(TP11)
applies. Reverse a separating functional on each rectangular arm as
in (TP15). Its polynomial \(K\) has long degree less than \(L\).
If \(K\notin J_L\), (MP8) gives the earlier reverse bound.

Suppose instead \(K\in J_L\) on an arm of valuation \(\nu\ge A/2\).
Split it by its monomial support:
\[
 K=K_0+K_1,\qquad K_0\in J_A,\qquad
 \operatorname{supp}K_1\subset
 \{(i,j):\max(0,L-\nu)\le i<A-\nu,\quad j<t-d\}.          \tag{MP10}
\]
The long interval has width at most \(t=A-L\). Here \(K_0,K_1\)
are still supported in long degrees below \(L\).
The full-length perturbation kills \(K_0\), while it kills \(K_1\)
after truncation to length \(L\).

The functional component corresponding to \(K_0\) annihilates every
cheaper original column. For upper columns this follows by taking a
subset of the original functional's zero row coefficients. On pure
central columns the length-\(L\) transpose product is the scalar
\(\gamma(K_0+K_1)\); the two supports are disjoint, so their zero
constraints separate coordinate by coordinate. On bridge central
columns the scalar part has no long degrees at least \(L\), and the
perturbation kills \(K_0\) even at full length \(A\). The opposite
arm contributes zero.

This component must also vanish on the exchanged label on its own arm:
otherwise it would separate the two nonzero parallel vectors modulo the
cheaper-column span. The full functional is nonzero there, so that
value comes from \(K_1\). Reflection of its long interval gives
\[
                \nu+t\le M<\min(a,\nu+2t),               \tag{MP11}
\]
where \(M=Y\) on the \(Y\)-arm and \(M=S\) on the \(Z\)-arm.
This reasoning covers both lower and upper exchanged labels: the
relevant coefficient is in the scalar transpose product for a source,
and in the functional itself for an upper identity.
If \(\nu\ge L\), all length-\(L\) vectors are annihilator exceptions,
but (MP10)--(MP11) remain valid with their possibly shorter interval.

## Count the extra strips, then use the existing fan argument

For a mixed pair write \(b=Y-R\), \(c=S-Z\) as in the transverse
proof. The partner short coordinate lies in \([0,2t)\). Thus (MP11)
forces
\[
                    |b-\nu_Y|<2t
       \quad\hbox{or}\quad |c-\nu_Z|<2t.                 \tag{MP12}
\]
There are \(O(at)\) raw label differences in these additional strips.
For same-layer directions normalize to a primitive vector only after
listing the raw differences; division need not preserve their translated
centers, but cannot increase their number.

Outside these strips no large-valuation dual exception remains. The four
old inequalities therefore hold with losses at most \(H_Y+t,H_Z+t\).
Their cross-layer diagonal strip has width \(20\Theta\). The same-layer
strip has width \(6\Theta\); its small-upper-threshold case uses the
primal estimates and original source boxes, without a dual assumption.
Together with the bounded-product inventory, these give (MP4).
All parallel subdivisions are retained through the original direction
theorem. The common-weight construction adds (MP12) to its arrangement;
the number and height estimates follow as before.

## Verification and boundary

~~~bash
python3 -B scripts/run_hook_monomial_perturbation_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checks use original commuting quotient matrices, full-subspace kernel
ranks for (MP8), and actual coherent exchanges. In dual-exception cases
they reconstruct the full-annihilator functional component and verify its
vanishing on every cheaper column before checking (MP11).
Validated null-row additions to some separating functionals exercise
nonzero removable components without changing either exchanged value.
Finite checks support the uniform proof; they are not a fan growth law.

Monomial annihilators, unit Toeplitz ranks, separating functionals and the
linked fan/arrangement bounds are the inputs. The explicit regular points
use the classical
[binomial determinant evaluation](https://www.mat.univie.ac.at/~slc/s/s42kratt.pdf#page=11).
Exact priority of this combined closure remains unassessed. Perturbations
with several incomparable lowest monomials, long residual-unit prefix gaps,
arbitrary outer images and the unrestricted endpoint remain open.
