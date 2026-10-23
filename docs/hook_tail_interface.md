# The far-arm part of the bridges controls a full hook fan

The bridge rank can be replaced by the rank of just its far-arm
projection. This gives a collective bound even when departure from
every transverse comparison is large. A bounded rational recurrence
on one arm makes this smaller interface explicit; the other arm can
remain arbitrary.

Use the original compatible flat hook model, over every field, with
\[
 M=(y^a,y^tz^t,z^a),\quad a\ge2t,\quad t\ge2,\qquad
 n=4at-2t^2,\quad A=a-t.
\]
In its residual row space \(T_Y\oplus T_Z\), let
\[
 F_Y=\langle Y_{q,v}:t\le q<A,\ 0\le v<t\rangle,
 \qquad F_Z=\langle Z_{q,v}:t\le q<A,\ 0\le v<t\rangle.
\]
These are original upper rows whose long exponent is at least \(2t\).
For the full three-block bridge set \(H\), define
\[
 \chi_Y=\operatorname{rank}(\pi_{F_Y}H),\quad
 \chi_Z=\operatorname{rank}(\pi_{F_Z}H),\quad
 \chi=\min(\chi_Y,\chi_Z).                              \tag{TI1}
\]
All compatible outer images are allowed. Then
\[
 \boxed{\quad g(I_\phi)=O(n^2(1+\chi)^3).\quad}          \tag{TI2}
\]
Since projection cannot increase rank, \(\chi\le b=\operatorname{rank}H\).
Combining the [departure-rank theorem](hook_rank_departure.md) gives
\(g(I_\phi)=O(n^2(1+\min(\chi,\rho))^3)\).
This is a uniform parameter theorem, not a new unrestricted worst-case
exponent or an actual polynomial fan separation. Exact priority and
the original final-paper requirements remain unresolved.

## The comparison matrix retains the complicated pure arm

Fix Y. In the original residue matrix \(V\), erase only the entries
of bridge columns in \(F_Y\), obtaining \(V_0\). Upper identities
and all pure-arm columns are unchanged. In particular both matrices
have their original full row rank \(R=2At\), and
\[
                    \operatorname{rank}(V-V_0)=\chi_Y.  \tag{TI3}
\]
This is an auxiliary vector matrix. It need not be the monomial
matrix of an algebra or a compatible hook map.

Every pure lower Y column is supported in \(F_Y\). Indeed its long
source shift is at least \(t\); all allowed Y image exponents are
nonnegative after the tag \(y^t\) is removed. Upper Y identities are
either in \(F_Y\) or in the short core \(q<t\). Thus \(V_0\)
splits as the direct sum of:

- a far-Y component containing only labels of the original group \(E_Y\);
- a component on the remaining rows, containing the core identities,
  the Z group, and the three modified bridge blocks.

The first component's greedy basis is fixed whenever the internal
\(E_Y\) order is fixed, regardless of its actual coefficients or rank
structure. No parallel-pair assumption or triangular inversion is used.

Every remaining nonzero label lies in
\[
 \{0,1\}\times[0,2t)\times[0,2a)
 \quad\cup\quad
 \{0\}\times[a,a+t)\times[0,t).                         \tag{TI4}
\]
The second box is the outer bridge \(H_0\); \(H_1,H_2\) lie in
the first box. Differences within and between these boxes have
\(O(at)\) values: their Y variation has width \(O(t)\), including
the translated second box, and their Z variation has width \(O(a)\).
Add the \(O(at)\) original three-group internal-order planes. A
coefficient-independent arrangement of \(O(n)\) planes therefore
fixes the entire comparison optimum and those orders on \(O(n^2)\)
convex chambers.

Apply the [common-quotient compression lemma](hook_rank_departure.md#a-common-quotient-compression-lemma)
to \(V,V_0\): quotient their output by \(\operatorname{im}(V-V_0)\).
On each chamber a fixed compulsory set and a minor of rank at most
\(3\chi_Y\) on at most \(9\chi_Y\) original labels recover the
whole original cone patches. The saved two-colour theorem counts
that minor by \(O((1+\chi_Y)^3)\). Transposition proves (TI2).
The auxiliary projection is not asserted to bound every cheaper
prefix's full arm coupling; the short core remains in the comparison
model and is counted through its label geometry.

## Rational recurrences bound the central interface

For a central-only Y image \(y^tW_Y\), assume that its two outer
Y images vanish; all Z components, including both outer images, may
be arbitrary. Suppose in \(R_A=K[q,v]/(q^A,v^t)\) that
\[
 W_Y=C(v)+N,\qquad N=v^\delta P(q,v)/Q(q),\quad
 Q(0)\ne0,\quad \deg_qP,\deg Q\le d,\quad0\le\delta<t.  \tag{TI5}
\]
Division means truncated formal-series division, so there is no
analytic convergence hypothesis. No unit assumption on \(C\) or
\(W_Y\) is required. Then
\[
 \chi_Y\le\min\{t^2,\ t(A-t),\ d(t-\delta)\},\qquad
        g(I_\phi)=O(n^2(1+d(t-\delta))^3).               \tag{TI6}
\]
This also includes \(d=0\), when the central Y arm is transverse
and its far bridge projection vanishes.

Here is a direct finite-matrix proof of the recurrence bound. Split
long degrees at \(t\), writing L for \([0,t)\) and H for \([t,A)\).
Multiplication by a polynomial in \(q,v\) is lower block triangular
in increasing long degree. From \(QN=v^\delta P\),
\[
 Q_{HH}N_{HL}=(v^\delta P)_{HL}-Q_{HL}N_{LL}.            \tag{TI7}
\]
The matrix \(Q_{HH}\) is invertible because \(Q(0)\ne0\).
Every vector on the right is supported in the first \(d\) long
rows above the cut and in short degrees at least \(\delta\).
Their span has dimension at most \(d(t-\delta)\), also when
the rectangle ends earlier. Thus \(N_{HL}\) has at most that rank.
The transverse part \(C(v)\) has no HL block. The original central
bridges give exactly the \(t^2\) low inputs, so their far-Y projection
is precisely \(N_{HL}\); the outer Y projections are zero by hypothesis.
The other two bounds in (TI6) are its column and row dimensions.

The relation between rational recurrences and finite Hankel rank is
classical; see [Bacher, Remarks 9.1--9.2](https://msp.org/ant/2010/4-7/ant-v4-n7-p05-p.pdf#page=22).
For one scalar series, reversing the low input order turns the
off-diagonal Toeplitz block into a shifted Hankel block. The precise
cut and short-degree bound used here are proved in (TI7); no new
rational-series rank theorem is claimed.

## An all-field family with large previous parameters

Take zero outer maps, \(A>t\), and on both arms put
\[
                  W=1+v/(1-q)\quad\bmod(q^A,v^t).       \tag{TI8}
\]
Each off-diagonal long block is an all-ones matrix tensored with
multiplication by \(v\). Therefore, over every field,
\[
                      b=t^2,\qquad \chi=t-1.            \tag{TI9}
\]
The bridge rank follows from the injective Y projection: \(W\) is
a unit. The new full bound is \(O(n^2t^3)\), rather than the coarse
\(O(n^2t^6)\) ceiling. At \(a=t^{23}\) these orders are
\(O(n^{17/8})\) and \(O(n^{9/4})\), respectively. Smaller saved
coefficient-specific bounds should still be retained.

This improvement is not explained by small departure rank. For every
transverse \(C(v)\), if its constant coefficient is not one, \(W-C\)
is a unit. If its constant coefficient is one but its coefficient
of \(v\) is not one, \(W-C\) is \(v\) times a unit, of rank
\(A(t-1)\). Otherwise its first \(v\)-coefficient is
\(q/(1-q)\). The \((A-1)(t-1)\) images of
\(q^iv^j\), \(i<A-1,j<t-1\), are independent by taking their
lowest nonzero short-degree coefficients successively. Choosing
\(C=1+v\) attains this lower bound. Consequently
\[
                         \rho=(A-1)(t-1).               \tag{TI10}
\]
This is of order \(n\) when \(A/t\) grows. Both earlier parameters
\(b\) and \(\rho\) can therefore grow while the smaller interface
still yields a polynomial improvement in the displayed regime.

Nor are the scalar Toeplitz prefixes regular. For
\(F(q)=(1-q)^{-1}\bmod q^m\), \(m\ge2\), every coefficient is one.
The prefix determinant at \(1\le r\le m-2\) has two equal final
rows and vanishes; the indices \(0,m-1,m\) are regular. Its longest
nonregular run is exactly \(m-2\). Thus bounded prefix gaps do not
explain (TI9), and the active transverse width is \(t-1\).

For the general compatible hook class, an actual polynomial excess
\(g\ge c n^{2+\epsilon}\) now requires
\(\chi_Y,\chi_Z=\Omega(n^{\epsilon/3})\), as well as the previous
departure-rank condition. Large interface rank is necessary, not
sufficient. The simultaneous large-interface, large-departure range
still contains the unresolved collective problem.

## Verification and priority limits

~~~bash
python3 -B scripts/run_hook_tail_interface_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker reconstructs original quotients and commuting multiplication,
tests the interface projections and auxiliary direct sums, and compares
entire original cone patches with the small minors. Independent polynomial
division and multiplication check (TI7). The geometric-series family
includes finite-field enumeration of all transverse comparison polynomials
at small sizes, exact bridge and interface ranks, and defective-prefix
checks. Complete original fans give separate coverage controls.

The auxiliary entry-erased matrices are explicitly distinguished from
algebra quotients. Common-quotient compression and the two-colour count
are retained with their existing attribution; rational/Hankel rank is
classical. The comparison does not establish strongest-current exact
novelty of this hook parameter theorem, sharpness, an unrestricted
exponent, an actual polynomial fan separation or a final paper.
