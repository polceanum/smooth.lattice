# Uniform fan bounds from two pole transports

Let \(K\) have characteristic \(p>0\), let \(q=p^e\ge2\), and let
\(F\in K[t]\) be monic of degree \(n>0\), with \(F(0)F(1)\ne0\).
For \(I_{F,q}=(xy-1,(x-1)^qz-1,F(x))\),
\[
 \begin{split}
 g(I_{F,q})&=O(q^9n^2(1+\log(n+1))),\\
 V(I_{F,q})&=O(q^9n^3(1+\log(n+1))),\\
 |\mathcal U(I_{F,q})|&=O(q^{10}n^2(1+\log(n+1))).
 \end{split}                                                \tag{1}
\]
The constants are absolute: no \(q,F,K\) or multiplicity dependence is
hidden. No splitting, perfectness or generic evaluation is required.
For every integer \(q\ge1\), over every field, the same fan estimate holds
for open initial cones meeting \(w_x<w_y\). The full estimate uses the prime-power
identity in the field characteristic.

Consequently a Frobenius-pole family with
\(g(I_{F_n,q_n})\ge c n^{2+\epsilon}\) requires
\[
 q_n=\Omega_{c,\epsilon}
       ((n^\epsilon/\log(n+1))^{1/9}).                       \tag{2}
\]
Fixed or subpolynomial pole orders cannot produce a polynomial separation
above the generic corner-cut scale through special quotient coefficients.
Nine is a loose upper exponent, not a sharpness assertion. This is a bound
with a growing parameter, rather than an \(O_C\) estimate with uncontrolled
carrier dependence.

For \(q=2\) in characteristic two, there is the sharper uniform result
\[
 g(I_{F,2})=\Theta(n^2),\qquad V(I_{F,2})=\Theta(n^3),\qquad
                         |\mathcal U(I_{F,2})|=\Theta(n^2). \tag{3}
\]
The unrestricted \(G_3\) question, a polynomial lower separation and exact
priority remain open. Frobenius, principal parts, greedy optimization and
the lattice-vertex theorem are existing tools. No novelty or paper-level
significance follows from this combination without its priority audit.

## Carriers and a universal complete-pole basis

The subsequent [weighted completion theorem](ordinary_pole_completion.md)
removes the characteristic restriction with a slightly weaker explicit
degree dependence, and covers all \(P/Q\) regular at zero and infinity,
including quotients where \(P\) is a nonunit. It uses a costed basis
completion; the exact-transport rank failure below remains valid.

Use \(x=t,y=t^{-1},z=(t-1)^{-q}\), normalize \(w_z=1\), and write
the other prices as \(\alpha,\beta\). In the strict chart
\[
 m<q\alpha/(\alpha+\beta)<m+1,\quad 0\le m<q,\quad a=q-m,
\]
the carrier ideal is \((xy,x^az,y^{m+1}z)\), by the integral
[pole-carrier proof](pole_order_fans.md). An initial downset has form
\[
 D(A,B,c)=\{1,x,\ldots,x^A,y,\ldots,y^B\}
   \cup\{[t^j]z^r:1\le r\le c_j,\ -m\le j<a\},\quad
                         n=1+A+B+\sum_jc_j.                \tag{4}
\]
Here \([t^j]=x^j\) for \(j\ge0\), and \(y^{-j}\) otherwise.
Ray lengths decrease away from \(j=0\); a positive \(c_j\) requires its
axis tag in the core. Set \(T=\min(\alpha,\beta)\) and
\(\tau_j=j\alpha\) for \(j\ge0\), \(\tau_j=-j\beta\) otherwise.
The chart inequalities imply
\[
                         0\le\tau_j\le qT.                \tag{5}
\]
Indeed \(m\beta\le a\alpha\) and
\((a-1)\alpha<(m+1)\beta\), so each extreme tag is bounded by
both \(q\alpha\) and \(q\beta\).

If \(A\ge a-1,B\ge m\) and all \(c_j=C\), the rational span is the
full space of functions with pole divisor
\[
                         A(\infty)+B(0)+qC(1).             \tag{6}
\]
At each positive layer the \(q\) consecutive Laurent tags supply all
principal parts of orders one through \(q\). Multiplication by \(t^m\)
turns their numerators into \(1,t,\ldots,t^{q-1}\), a unit triangular
basis in powers of \(t-1\). Induct on the layer, removing allowed \(0\)
and \(\infty\) parts. No factorial division is used. At size \(n\),
clearing the denominator identifies the span with polynomials of degree
less than \(n\), so it is a basis modulo every \(F\).

For \(n\ge q\), start with \(A=a-1,B=m\), add the remainder of \(n-q\)
modulo \(q\) to the cheaper axis, and use equal ray lengths. This basis
has all element costs at most \(L=2qT+n/q\). Matroid greedy minimizes
the largest chosen cost as well as the total. Thus every actual basis has
\[
                         A,B\le2q+n/(qT).                  \tag{7}
\]
This conclusion holds for every integer pole order in every characteristic.

## Two exact rank-preserving transports

Write \(V(A,B,c)\) for the rational span in (4). If
\(A\ge a-1+q,B\ge m\), then
\[
 zV(A,B,c)=V(A-q,B,c+\mathbf1).                             \tag{8}
\]
The axis space maps to the full space
\(L((A-q)\infty+B0+q1)\), supplying the whole new base layer.
Old positive layers advance one step. This proves equality of spans.

When \(q=p^e\) in characteristic \(p\), \((t/(t-1))^q=z+1\).
If \(A\ge a-1,B\ge m+q\), then
\[
 (z+1)V(A,B,c)=V(A,B-q,c+\mathbf1).                         \tag{9}
\]
The axis space maps to \(L(A\infty+(B-q)0+q1)\). The remaining columns
are \(t^j(z^{r+1}+z^r)\); induction using the new base layer gives
exactly the extended rays. Both multipliers are units modulo every \(F\).
Hence both identities preserve all residue ranks, including deficient ranks.

The characteristic hypothesis is necessary for this second transport.
For \(q=2\) in characteristic zero, take \(F=1+(t-1)+(t-1)^5\).
In the \(y\)-tag chart, \(\{1,y,y^2,y^3,z\}\) is a basis but its
proposed shift \(\{1,y,z,z^2,yz\}\) has rank four. The
[quartic Hurwitz tests](quartic_hurwitz_fans.md) are respectively \(-2,0\).

## Missing base tags have only a small collective profile box

Assume \(n\ge4q\). If \(A<a-1\) or \(B<m\), a required tag of cost
at most \(qT\) was rejected. All carrier monomials of cost at most \(qT\)
have poles bounded by
\[
                         q(\infty)+q(0)+q\lfloor qT\rfloor(1).
\]
This space injects modulo \(F\) whenever \(2q+q^2T<n\): cleared
numerators have degree less than \(n\). All its carrier monomials would
then be independent and accepted. Thus a missing tag implies
\(T\ge n/(2q^2)\), and (7) gives both \(A,B\le4q\).
Its actual basis profiles therefore lie in a box with sides
\(4qn,4qn,n^2\). Each remains a vertex of the convex hull of this subset
of profiles. Add a unit cube for a flat hull. The classical
three-dimensional lattice-vertex bound gives \(O(qn^2)\) such states,
collectively across all charts. If \(n<4q\), the established
unrestricted \(O(n^3)\) bound is already \(O(qn^2)\).

## A variance budget without enumerating all ray shapes

Now \(A\ge a-1,B\ge m\). Let \(R\) be the cheaper axis length; keep
the other fixed. Put \(S=\sum c_j=M-R,\ v_j=c_j-S/q\),
\(P=\sum\tau_j\). Up to a term independent of \(R,v\), the cost is
\[
 g(R)+E(v),\quad
 g(R)=\frac{qT+1}{2q}R^2+
       \left(\frac T2-\frac Mq-\frac12-\frac Pq\right)R,
 \quad E(v)=\tfrac12\sum v_j^2+\sum\tau_jv_j.                \tag{10}
\]
These divisions are real exponent costs, not arithmetic in \(K\).

If \(R\) exceeds its required base length by at least \(q\), apply
(8) or (9). Optimality gives \(g(R)\le g(R-q)\), hence
\(g'(R)\le(qT+1)/2\). Otherwise \(R<2q\), giving \(g'(R)\le2qT\).
Raise \(R\) by \(\delta=S\bmod q\) and balance all ray lengths at
\(\lfloor S/q\rfloor\). This is a valid universal basis (6).
Its \(g\)-cost increases by at most \(4q^2T+q\). Therefore
\[
 E(v)\le4q^2T+q,\qquad
                    \|v\|_2\le8q^{3/2}\max(1,T).           \tag{11}
\]
For the norm bound use Cauchy--Schwarz and
\(\|\tau\|_2\le q^{3/2}T\). For arbitrary \(q,K\), this argument
still works on \(\alpha\le\beta\), using (8) alone.

At fixed \(A,B\), consider the integer tag moments
\[
                         X=\sum_{j>0}jc_j,\quad
                         Y=\sum_{j<0}(-j)c_j.
\]
Each differs from its equal-ray value by at most
\(8q^3\max(1,T)\), since the coefficient vectors have norm at most
\(q^{3/2}\). For fixed \(A,B,X,Y\), only the positive-\(z\) moment
\(\sum c_j(c_j+1)/2\) remains in the cost. Its smallest feasible value
is independent of the weights. A tie gives the same whole profile and
cannot expose two different bases at a generic positive objective.
Thus there is at most one actual state per \(A,B,X,Y\).
This avoids enumerating \(q-1\) independent-looking ray differences.

## Count scales with an explicit pole-order power

For \(0<T<1\), there are \(O(q^6)\) moment pairs per \(A,B\), at most
\(n^2\) axis pairs and \(q\) charts. For \(L\le T<2L\), \(1\le L<n\),
the count per chart is at most
\[
 \left(2q+1+\frac n{qL}\right)^2(33q^3L)^2
 \le17424q^8n^2,                                          \tag{12}
\]
since \((2q+1)L+n/q\le4qn\). At most
\(\lceil\log_2 n\rceil\) dyadic intervals are needed. For \(T\ge n\),
(7) gives \(O(q^2)\) axis pairs and the trivial \(0\le X,Y\le qn\)
gives \(O(q^2n^2)\) pairs per axis pair. Together with the missing-tag
bound, these estimates prove (1).

Each reduced basis has at most \(q+5\) polynomials: three carrier heads
and at most \(q+2\) frontiers. Also \(V\le3ng\) by planar adjacency
and the primitive wall threshold. This proves the other estimates in
(1) for polynomial counts rather than coefficient bits.

## The quartic logarithm disappears in characteristic two

In the quartic notation,
\(G_{b+2}(s)=G_b(s)/(1+s^2)\) in characteristic two. This even unit
acts on the odd-coefficient Hurwitz block by unit lower triangular row
operations. Thus \(\Delta_h(b+2)=\Delta_h(b)\) for every \(h,b,F\).
The \(x\)-cheaper chart is already quadratic over every field.
In the other chart fix the pure \(x\)-length and \(h\)-parity. The
\(y\)-length has fixed parity, so its feasible \(h\)-set is independent
of that length. Its cost has the same separated quadratic form, with
the axis prices exchanged. The quartic uncrossing argument gives
\(O(n)\) minima per fixed \(x\)-length and parity, proving the quadratic
upper bound. The all-field quadratic lower fan/basis-union and cubic
lower wall cost in the quartic note prove (3). No Frobenius-power
condition is imposed on \(F\).

## Inputs, controls and remaining scope

[Babson--Onn--Thomas, Proposition 2.2 and Theorem 2.4](https://arxiv.org/pdf/math/0207135)
supply the established profile/initial dictionary; Corollary 2.5 gives
the small-\(n\) fallback. The box step uses the classical
[Andrews theorem](https://doi.org/10.1090/S0002-9947-1963-0143105-7);
[Dillon, Theorem 1 and proof](https://arxiv.org/pdf/2606.30856) supplies
an inspected primary proof, including its relative-lattice induction.
No new lattice theorem, Frobenius identity or greedy algorithm is claimed.
Exact priority of the transported spans and uniform pole-order bound is
not established by these comparisons.
The proof of Proposition 1.2 in
[Caviglia--Sbarra](https://people.dm.unipi.it/sbarra/pdfs/Research/rev-regbounds-0.3.pdf)
also records classical commutation of Frobenius powers with ordinary
initial ideals, before its generic-initial conclusion for homogeneous ideals.
Our modulus \(F\) is arbitrary and the whole ideal is not assumed to be a
Frobenius power. That commutation statement is prior background; it does
not alone give the variance or growing-pole count proved here. Only the
proposition and proof were compared, not the regularity arguments.

~~~bash
python3 -B scripts/run_frobenius_pole_checks.py
~~~

The report is results/research/frobenius_pole_checks.json. Cleared
numerator matrices check rational-function spans before reduction modulo
\(F\). Independent remainder columns check ranks, complete fans,
both transports, universal bases, bottleneck and variance budgets and
conditional profile injections. The non-Frobenius failure is retained.
Exact characteristic-two controls check Hurwitz periodicity and both
one-parameter envelopes. Finite checks support proofs, not novelty.

The estimates in this note do not cover both price regions at arbitrary
non-Frobenius orders or general rational numerators. The
[weighted completion theorem](ordinary_pole_completion.md) now supplies
that extension without assuming a second exact transport. All fixed
curves and unrestricted \(G_3\) remain outside these results. A next
broad step must control joint-rank defects across a growing carrier or
across unrestricted walls. The original goal and its publication criteria
remain active.
