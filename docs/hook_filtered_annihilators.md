# Filtered annihilator exceptions are greedy basis differences

The scalar-annihilator exceptions in a central hook admit an exact
description for arbitrary mixed operators. No monomial-times-unit
factorization is required. Each exception uses one of at most \(t^2\)
labels in a difference of two planar greedy bases. The dual reduction
preserves original cheaper-column constraints by a filtered linear
projection, rather than by deleting selected coefficients.

This reduces the unresolved mixed case. A small set at each weight does
not count cross-arm walls across weights. No unrestricted exponent
improvement, polynomial lower separation or exact novelty is asserted.

## Four bases for one arm

Put \(A=a-t\), \(L=A-t>0\), and
\[
 R_A=K[q,v]/(q^A,v^t),\qquad W=\alpha+N,\quad
 \alpha\ne0,\quad N\in(v)R_A.                         \tag{FA1}
\]
The mixed polynomial \(N\) is arbitrary; write \(T_A\) for its
multiplication map. Order monomials increasingly by
\(\rho(q^iv^j)=iu+j\), with \(u>0\) separating the rectangle's
monomials. Use this same order on subsets. Set
\[
 E_A=[0,A)\times[0,t),\quad P=[t,A)\times[0,t),\quad
 E_L=[0,L)\times[0,t).
\]
Define four increasing-cost greedy column bases:

| Basis | Input labels | Output map |
|---|---|---|
| \(B\) | \(E_A\) | \(T_A\) |
| \(B_P\) | \(P\) | \(T_A\) |
| \(B_+\) | \(E_L\) | \(T_A\), retaining all output rows |
| \(B_-\) | \(E_L\) | \(T_L\), discarding output rows with \(q\ge L\) |

Here \(T_L\) is multiplication in \(R_L\). Commutation with \(q^t\)
identifies \(B_P=q^tB_-\), with the same translated label order.
Define
\[
 S_P=B_P\setminus B,\qquad S_D=B_+\setminus B_-.        \tag{FA2}
\]
Then
\[
 B\cap P\subseteq B_P,\quad B_-\subseteq B_+,\qquad
 |S_P|\le t^2,\quad |S_D|\le t^2.                     \tag{FA3}
\]
For the first inclusion, dependence on earlier inputs in \(P\)
implies dependence on earlier full inputs. Adding the \(t^2\) inputs
outside \(P\) removes at most that many selected old inputs; precisely
\[
 |S_P|\le t^2-(\operatorname{rank}T_A-\operatorname{rank}T_L).
\]
For the second, independence after output projection implies
independence before it. The rank difference is at most the number
\(t^2\) of discarded output rows. These are all-field matroid statements.

## Primal membership uses actual pure sources

Consider the central-only ideal
\[
 I=(x^2,y^a,z^a,y^tz^t-xf),\quad a>2t,\qquad
 n=4at-2t^2,                                          \tag{FA4}
\]
with nonzero constant axis coefficients in both arm units and zero
outer generator images. Use the original quotient columns of the
[bridge reduction](hook_bridge_reduction.md).

Take an eligible pure cross-arm exchange at a positive wall weight,
with no internal arm ties. Eligibility means both exchanged columns
are nonzero and parallel modulo strictly cheaper original columns.
On one arm suppose its primal restriction \(H\) satisfies \(NH=0\).
Let \(p\) be the local label of that arm's exchanged column, whether
a lower source or an upper identity. Assume its local long coordinate
is at least \(t\), so it is pure. Then
\[
                              p\in S_P.               \tag{FA5}
\]
Explicitly, a Y source \((0,Y,Z)\) has local label
\((Y-t,Z-t)\), and a Y upper identity \((1,Y,Z)\) has
\((Y-t,Z)\). The Z labels interchange the long and short coordinates:
\((S-t,R-t)\) for \((0,R,S)\), and \((S-t,R)\) for
\((1,R,S)\). On Z use parameter \(1/u\) and multiply local
weights by \(u\). Its reflected label in (FA8) is
\((A-1-p_q,t-1-p_v)\).

Indeed \(WH=\alpha H\). The original source/output inequalities
make \(p\) the unique maximum-weight monomial of \(H\). Hence
\(T_Ap\) is a combination of earlier full inputs, so \(p\notin B\).

If also \(p\notin B_P\), there is a relation
\(G=p+\text{earlier pure inputs}\) in \(\ker T_A\).
Thus \(WG=\alpha G\) uses only that pure arm. For a lower exchanged
source, the primal inequalities imply that its corresponding upper
identity is cheaper: \(\lambda<t\) on Y and \(\lambda<tu\) on Z.
Every other source and image of \(G\) is cheaper too. For an upper
exchange these two conditions reverse, making every source and every
other upper image cheaper. Either case makes the exchanged column
a loop modulo cheaper original columns, a contradiction.

The argument uses an actual pure-source relation. The kernel may be
nonmonomial; it does not exclude all primal annihilator relations.

## The dual projection preserves every cheaper constraint

Assume the overlap rows are cheaper than the exchanged pair, as in
(TP11) of the [transverse-prefix proof](hook_transverse_prefix.md).
Reverse an actual separating functional's arm rectangle to obtain
\(K\in R_L\), as in (TP15), and suppose \(T_LK=0\).
There is a unique \(K_1\) supported on \(B_+\) satisfying
\[
 T_AK_1=T_AK,\qquad K_0=K-K_1\in\ker T_A.              \tag{FA6}
\]
Greedy reduction has the filtration property
\[
                  \rho(K_0),\rho(K_1)\le\rho(K).        \tag{FA7}
\]
Reduce each input monomial separately: its representation in the
greedy basis uses only inputs no later than itself. Linearity proves
(FA7) for every vector. Also \(T_LK_1=T_LK_0=0\), and both inputs
have long degree less than \(L\).

Let \(p\) be the reflected local label of this arm's exchanged
column. Then
\[
                (K_1)_p=K_p\ne0,\qquad p\in S_D.        \tag{FA8}
\]
The original-column justification is necessary. Since
\(WK=\alpha K\bmod q^L\), the input and product inequalities in
(TP16) are two bounds on the same maximum weight. They hold for
\(K_0\) by (FA7). At full length \(WK_0=\alpha K_0\).
After reversal its functional therefore annihilates every cheaper
upper identity and pure central source. It annihilates every bridge
too: its scalar product has no long degrees at least \(L\), and its
full nonscalar product is zero. Its opposite-arm part is zero, and
the outer-source columns are zero by hypothesis.

This functional annihilates all cheaper original columns and the
opposite exchanged column. Parallelism forces it to vanish on its own
exchanged column too. Thus \((K_0)_p=0\). The input/product equality
makes \(p\) the leading monomial of \(K\), hence of \(K_1\).
A leading monomial in a nonzero kernel relation for \(T_L\) cannot
belong to its increasing greedy basis. Since \(K_1\) is supported on
\(B_+\), this proves (FA8).

Both arguments apply to cross-layer and same-layer exchanges under
the stated no-internal-ties hypotheses; the dual also requires cheap
overlap. Directions in the internal-arm inventory are already retained.
The same-layer small-upper-threshold branch is not covered here.

## Aggregate information and the remaining question

All four bases are constant between the rational numbers
\[
       u=j/i,\qquad 1\le i<A,\quad1\le j<t.             \tag{FA9}
\]
There are at most \((A-1)(t-1)\) such numbers. Include their
reciprocals for the other arm. This gives a coefficient-independent
partition into \(O(at)\) open intervals. On each interval each arm has
at most \(2t^2\) exceptional labels, so there are \(O(at^3)\)
interval/label incidences. Endpoints belong to the retained internal
direction inventory.

This compiles the annihilator exceptions covered by (FA5) and (FA8)
from exact planar ranks.
It does not give that many fan cells: one exceptional label may
interact with many opposite-arm labels, and \((\lambda,u)\) still
varies. Recurrent pieces must be retained. The restricted-input
matroid defining \(B_+\) is not declared to be an algebra quotient.

The next question is how often these moving basis differences
participate in simultaneous nonzero parallel exchanges, using the
actual coupling and nonscalar uncertainty together. Replacing it by
the pointwise \(t^2\) estimate leaves the counting problem unproved.

## A genuinely mixed nonmonomial control

Choose \(t=d+2e\), \(d,e\ge1\), \(A\ge2t\), \(1\le\mu<t\),
\(\nu=A-t+\mu\), and
\[
 U=1+q^\mu v^{d+e}+q^\nu v^d,\qquad
 V=1+q^{A-1-e}v.                                      \tag{FA10}
\]
Both incomparable corners of \(U-1\) have positive long valuation.
Its kernel contains \(H=q^{A-t}-v^e\), while both monomials separately
have nonzero image. In fact
\[
 \operatorname{ann}(U-1)=(q^{A-\mu},v^e-q^{A-t})\subset R_A,
 \qquad \operatorname{rank}T_A=e(A-\mu).                \tag{FA11}
\]
The displayed ideal annihilates \(U-1\); its quotient has dimension
\(e(A-\mu)\). This is also the multiplication rank. The images of
\(q^iv^j\), \(i<A-\mu\), \(j<e\), have distinct lowest-long-degree
terms \(q^{i+\mu}v^{j+d+e}\), giving a triangular independent set.

The original pair \(E=(0,A,t)\), \(F=(1,1,a-1)\) remains eligible
throughout \(e<\lambda<t\), \(u=1+\lambda/(A-1)\).
The cheaper circuit in (TB11)--(TB12) uses source \(q^{A-t}-v^e\).
A separator has coefficients \(+1\) on the upper Y row
\((A+\mu,d+e)\), and \(-1\) on the upper Z row \((1,a-1)\).
The first upper-row cost exceeds the pair's cost by
\(\mu u+\lambda-e>0\); the matching pure Z source exceeds it by
\(tu-\lambda>0\). All other nonzero evaluations are on the pair or
more expensive columns, except the shared cheaper bridge, which
cancels. Thus eligibility holds over every field. At generic points
of this interval it gives an actual wall. No full fan count or
assertion outside the interval is inferred.

## A same-layer wall with a nonmonomial dual exception

A complementary family exercises the truncated-kernel branch. Choose
\(t=d+2e\), \(1\le\mu<t\), and an integer \(A\) such that
\[
 3t-\mu<A<4t-\mu,\qquad h=A-2t,\qquad hA/(2t)>e.
\]
Set \(a=A+t\), \(L=A-t\), \(r=5t-A-\mu-1\), and
\[
 U=1+q^\mu v^{d+e}+q^{h+\mu}v^d,\qquad
 V=1+q^{A-2e}v^{2t-1-r}.                              \tag{FA12}
\]
These inequalities give \(t\le r<2t-1\), so both axis coefficients
are constant. With \(K_Y=q^h-v^e\),
\[
 N_YK_Y=q^{\mu+2h}v^d\ne0\ \text{in }R_A,\qquad
 N_YK_Y=0\ \text{in }R_L.                              \tag{FA13}
\]
Each monomial of \(K_Y\) separately has nonzero image in \(R_L\).
Thus this is a nonmonomial truncated annihilator exception rather than
a vector supported on a monomial kernel.

Take
\[
 E=(1,3t-1,t-1),\quad F=(1,t-1,a-1),\quad
 w=(\lambda,A/(2t),1),\quad0<\lambda<t.                 \tag{FA14}
\]
Their costs agree. The central bridge
\(B=(0,r,t+2e-1)\) has residue equal to the sum of \(E,F\)
and the two upper identities
\((1,r,2e-1)\), \((1,r-t,t+2e-1)\). Those three other columns
are strictly cheaper, proving parallelism modulo cheaper columns.

For nonzero quotient rank, reverse \(K_Y\) on Y and the constant
\(K_Z=-1\) on Z to obtain a functional. Its upper evaluations are
\(1,-1\) on \(E,F\); the other nonzero Y upper evaluation lies above
them by \(hu-e>0\). Its central evaluations cancel on \(B\).
The remaining three central evaluations exceed the pair's cost by
\(t-\lambda\), \(hu-e+t-\lambda\), and \(tu-\lambda\).
Thus the functional kills all cheaper original columns. This proves
eligibility throughout the interval over every field. The common
overlap upper rows are cheaper as well. Whenever the wall has no
internal arm ties, both dual memberships (FA8) apply. Such choices
include the exact fixtures in the checker; no count is extrapolated
from them.

## Reproduction and priority

~~~bash
python3 -B scripts/run_hook_filtered_annihilator_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

Checks cover complete filtered kernel spaces, with their dimensions
independently checked by dense elimination, original commuting quotient
columns, and both memberships on actual walls. Nonzero removable dual
components are included only after checking their original-column
constraints. The two-positive-valuation family is checked in
characteristic zero and small and larger finite characteristics.

Greedy deletion, output projection and prefix-preserving basis reduction
are classical linear/matroid operations; compare
[Eppstein, Sections 1.2 and 4](https://ics.uci.edu/~eppstein/pubs/Epp-DCG-98.pdf).
The linked hook notes retain the quotient and state-polytope inputs.
This application has no exact priority certificate. Arbitrary outer
images, aggregate mixed-wall counts, unrestricted sharp growth and
the final-paper requirements remain open.
