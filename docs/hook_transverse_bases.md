# Transverse units preserve the hook direction bound

The constant base of the monomial-perturbation theorem can be replaced by
an arbitrary unit depending only on the transverse variable. This covers
some perturbations with several incomparable lowest monomials and
nonmonomial annihilators. The useful decomposition separates the transverse
unit from the mixed perturbation; a projection by whole long-coordinate
fibers preserves the original cost constraints.

An explicit all-field family below has a coherent wall whose primal source
lies in the nonmonomial annihilator of the full off-axis perturbation.
Thus the earlier primal exclusion does not extend to that annihilator
without a new argument.

The restricted closure and this obstruction do not settle unrestricted
\(G_3\), produce polynomial excess over the generic scale, or establish
exact novelty.

## Family and count

Use the central-only hook
\[
 I=(x^2,y^a,z^a,y^tz^t-xf),\quad
 n=4at-2t^2,\quad t\ge2,\quad a\ge2t,\quad A=a-t,\quad L=a-2t.
                                                               \tag{TB1}
\]
In its two rectangular arm modules assume
\[
 \begin{split}
 U&=C_Y(z)+y^{\nu_Y}z^{d_Y}B_Y(y,z),\\
 V&=C_Z(y)+z^{\nu_Z}y^{d_Z}B_Z(z,y),
 \end{split}                                                   \tag{TB2}
\]
where \(C_Y,C_Z\) have nonzero constant coefficients, \(B_Y,B_Z\)
are units, \(0\le\nu_i<A\), and \(1\le d_i<t\).
The higher coefficients of the \(C_i\)'s are arbitrary.
The mixed terms retain their displayed common long-variable factors.
Both outer generator images are zero.

Let \(H_i\) be the maximum unit Toeplitz prefix gap of
\(G_i(q)=B_i(q,0)\) at positive lengths \(A-\nu_i,L-\nu_i\).
Set \(\Theta=3t+H_Y+H_Z\), as in (MP3).

**Theorem.** Uniformly in all these coefficients, valuations and fields,
\[
 D_{\rm all}(I)=O(a\Theta(1+\log(a+1))),\qquad
                 g(I)=O(na\Theta(1+\log(n+1))^2).         \tag{TB3}
\]
In particular short gaps \(H_Y+H_Z=O(t)\), including regular residual
prefixes, give \(g(I)=O(n^2\log^2(n+1))\).
Keep the earlier box bound when smaller.
The common complete weight list from the monomial-perturbation proof
extends with the same cardinality and height, fixing
\(t,a,\nu_Y,\nu_Z,\Theta_0\):
\[
 |\mathcal W|=O(a^2\Theta_0^2(1+\log(a+1))^2),\qquad
                         \|w\|_\infty=O(a^2).            \tag{TB4}
\]
No common list of this size across all valuation choices at once is asserted.

If one arm depends only on its transverse variable, the near-quadratic
bound holds with no condition on the other constant-axis arm. Its pure
long-coordinate fibers are direct summands; all cross-arm labels left after
removing them have short long-coordinate displacement.

## Long-degree uncertainty survives the base change

Write one arm as \(W=C(v)+N\) in
\(R_m=K[q,v]/(q^m,v^t)\), with \(N=q^\nu v^d B\).
Multiplication by \(C(v)\) and by its inverse preserves long degree
exactly: the coefficient of the largest \(q\)-power is multiplied by a
unit of \(K[v]/(v^t)\). Hence
\[
 \deg_q(WH)=\deg_q(C^{-1}WH),\qquad
 C^{-1}W=1+q^\nu v^d(C^{-1}B).                            \tag{TB5}
\]
The axis unit of \(C^{-1}B\) is \(G(q)/C(0)\), whose prefix
regularity and gaps are unchanged.
All long-degree uncertainty bounds from (MP8)--(MP9) therefore apply.
Their maximum-weight lower bounds follow from nonnegative exponents.
This normalization is used for degree estimates; no equality of the two
monomial fans under the row transformation is asserted.

For large valuation \(2\nu\ge A\), its exceptions are
\[
 J_m=\operatorname{ann}(N)=(q^{m-\nu},v^{t-d})R_m,         \tag{TB6}
\]
with \(J_m=R_m\) if \(\nu\ge m\). On that subspace
\(WH=C(v)H\). For smaller valuations the previous defect allowance is
at most \(H_i+t\).

## Primal exceptions are excluded fiber by fiber

Suppose the source restriction \(H\) of an eligible pure cross-arm
exchange belongs to \(J_A\). Write
\(H=\sum_i q^i H_i(v)\).
Each long-coordinate fiber separately belongs to \(J_A\), and \(C(v)\)
does not mix these fibers.

Take the fiber of the exchanged long coordinate. All its source labels
are pure on that arm, since the exchange is outside the known short-label
inventory. Its vector identity is multiplication by \(C(v)\).
If the exchanged column is a lower source, the other source labels in
this fiber and every upper label in its image are strictly cheaper,
by the original exchange relation. This makes the lower column a loop
modulo cheaper columns. If it is an upper identity, all source labels
in the fiber and all the other upper labels in its image are cheaper,
giving the same contradiction.

Thus the primal exception is still ineligible. This argument uses a
fiber relation and does not claim that an individual source is parallel
to an upper monomial when \(C(v)\) is nonconstant.

## Dual exceptions use a whole-fiber projection

At a long exchange with cheap overlap rows, let \(K\), of long degree
less than \(L\), be the reversed separating functional from (TP15).
If \(NK\ne0\) modulo \(q^L\), the reverse uncertainty estimate applies.
If \(K\in J_L\), split it using only the long-coordinate interval:
\[
 K=K_0+K_1,\qquad
 \operatorname{supp}K_1\subset
 \{(i,j):\max(0,L-\nu)\le i<A-\nu,\ 0\le j<t\}.           \tag{TB7}
\]
Keep every short-coordinate coefficient in each retained fiber.
The remaining fibers form \(K_0\). A fiber below \(L-\nu\) in
\(J_L\) is divisible by \(v^{t-d}\); a fiber at or above \(A-\nu\)
is killed by the long cutoff. Consequently \(NK_0=0\) at full length
\(A\), and \(NK_1=0\) at length \(L\).

The original input constraints on \(K\) and output constraints on
\(WK=C(v)K\) both survive this projection. Their fibers are disjoint,
and \(C(v)\) preserves the long coordinate. Thus \(K_0\) defines a
functional annihilating every cheaper upper and pure central column.
It also annihilates bridge sources: its scalar-base product has no long
degrees at least \(L\), and its mixed perturbation vanishes at full length.
The opposite arm contributes zero.

It cannot separate the exchanged pair. The exchanged value must come from
the \(K_1\) interval, so the same translated strip follows:
\[
                      \nu+t\le M<\min(a,\nu+2t).         \tag{TB8}
\]
The projection need not preserve individual short-coordinate coefficients;
it preserves both weight filtrations through the whole-fiber split.

The extra raw direction strips still have \(O(at)\) elements. All other
primal and dual inequalities have losses at most \(H_i+t\).
The previous bounded-product, diagonal-strip, same-layer and parallel-piece
arguments now prove (TB3)--(TB4).

For a wholly transverse arm \(W=C(v)\), every pure long-coordinate
fiber has only its own lower and upper columns and is a direct summand.
All comparisons within such fibers have zero long-coordinate displacement.
After removing them the remaining label inventory has \(O(at)\)
directions. This proves the stated degenerate case.

## A coherent primal relation in a nonmonomial annihilator

Choose \(t=d+2e\) with \(d,e\ge1\), \(A\ge2t\), \(a=A+t\),
and \(\nu=A-t\). In the two local arm coordinates set
\[
 U=1+z^{d+e}+y^\nu z^d,\qquad
 V=1+y z^{A-1-e}.                                       \tag{TB9}
\]
These are instances of (TB2). Their residual prefix gaps are at most
\(t-1\) and \(e\), so (TB3) is near-quadratic even when \(t\) grows.

In the \(Y\)-arm,
\[
 (U-1)(y^\nu-z^e)=y^{2\nu}z^d-z^{d+2e}=0.              \tag{TB10}
\]
Both monomials separately have nonzero image under \(U-1\).
Thus this annihilator is nonmonomial. In fact it is
\((z^e-y^\nu)\): that ideal annihilates \(U-1\), its quotient
has dimension \(eA\), and the monic high-transverse coefficients show
that multiplication by \(U-1\) has rank \(eA\).

Consider the original labels
\[
 E=(0,t+\nu,t)=(0,A,t),\qquad F=(1,1,a-1),\qquad
                  F-E=(1,-(A-1),A-1).                   \tag{TB11}
\]
For every \(e<\lambda<t\), put
\[
                  w=(\lambda,u,1),\qquad
                  u=1+\lambda/(A-1)>1.                 \tag{TB12}
\]
Their costs are equal. The global source polynomial \(y^\nu-z^e\)
has \(Y\)-restriction (TB10) and \(Z\)-restriction \(-z^e\).
Its residue identity expresses \(E+F\) using the strictly cheaper
central bridge \((0,t,t+e)\) and upper labels
\[
                    (1,t+\nu,0),\quad(1,t,e),\quad(1,0,t+e).
\]
Hence the two columns are parallel modulo cheaper columns.

They are nonzero there. Take a functional equal to \(1\) on the upper
\(Y\)-row \((y^\nu,z^{d+e})\), equal to \(-1\) on the upper
\(Z\)-row \((y,z^{A-1})\), and zero on other reference rows.
Its only shared cheaper source contribution is the bridge above, where
the two values cancel. The other columns on which it is nonzero
are either \(E,F\) or more expensive. Its \(Y\)-upper row costs more
than the pair by \(\lambda-e>0\); its matching pure \(Z\)-source
costs more by \(tu-\lambda>0\). Thus it annihilates every cheaper
column and is nonzero on both \(E,F\).

At a generic point of this open interval this is an actual coherent wall,
over every field. The assertion is an interval of eligibility, not a full
fan count or a claim about all weights outside that interval. It refutes
the inference that a primal restriction in
\(\operatorname{ann}(U-U(0,0))\) must always be ineligible.
The useful monomial annihilator in the new proof belongs to the mixed
operator \(U-C_Y(z)\), which is a different operator.

## Verification and priority

~~~bash
python3 -B scripts/run_hook_transverse_base_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

Checks reconstruct original commuting quotient columns, test long-degree
and weighted subspace kernels with nonconstant transverse bases, verify
whole-fiber functional splits on coherent walls, and independently certify
(TB9)--(TB12), including its nonmonomial annihilator and two adjacent
strict minima. Validated null-row additions exercise nonzero removable
fiber components. No finite example is extrapolated into a growth law.

The unit-prefix input retains the classical
[binomial determinant evaluation](https://www.mat.univie.ac.at/~slc/s/s42kratt.pdf#page=11).
The degree-preserving unit multiplication, monomial mixed kernels,
original rank separation and the linked direction/arrangement bounds are
the inputs. Exact priority of this combined extension remains unassessed.
Several incomparable mixed monomials after removal of the transverse
unit, long residual prefix gaps, arbitrary outer images and unrestricted
sharp growth remain unresolved.
