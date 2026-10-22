# Central prefix minors control the outer hook arms

The [hook bridge reduction](hook_bridge_reduction.md) leaves two kinds of
pure lower columns on each arm. The following uniform result puts the
outer kind within the existing quadratic wall budget on an explicit
nonempty open parameter locus, over every characteristic. It leaves the
central-generator wall count unresolved.

Let \(t\ge2\), \(a\ge2t\), \(B=K[y,z]\), and
\[
 M=(y^a,y^tz^t,z^a),\quad C=B/M,\quad A_\phi\text{ the compatible flat
 square-zero extension},\quad n=4at-2t^2.
\]
Write \(f_i=\phi(g_i)\) for the three generator images. The coefficients
of \(f_0,f_2\) remain arbitrary throughout. No full support-matching rank
hypothesis is imposed.

## A uniform open condition involving only the central image

Put \(A=a-t\). For \(t\le r\le A\), define two labelled rectangles
\[
\begin{split}
 P_Y(r)&=\{(1,i,j):t\le i<t+r,\ 0\le j<t\},\\
 Q_Y(r)&=\{(0,i,j):2t\le i<a+t-r,\ t\le j<2t\}.             \tag{CP1}
\end{split}
\]
They contain \(tr\) and \(t(A-r)\) columns. All their images lie in
\(T_Y=y^tC\), of dimension \(tA\). The plain columns are shifted
images of \(f_1\); their opposite-arm component vanishes. After the
identity columns \(P_Y(r)\) are eliminated, let \(\Delta_Y(r)\) be the
square minor on rows
\[
                  t+r\le i<a,\qquad0\le j<t.                 \tag{CP2}
\]
The empty minor at \(r=A\) is one. Define \(P_Z,Q_Z,\Delta_Z\) by
interchanging \(y,z\). Call \(f_1\) prefix regular when
\[
          \prod_{r=t}^{A}\Delta_Y(r)\Delta_Z(r)\ne0.          \tag{CP3}
\]
Then \(P_Y(r)\cup Q_Y(r)\) is an actual basis of \(T_Y\) for every
\(r\), and similarly on the other arm.

Every determinant in (CP3) is a primitive polynomial over the integers.
For \(\Delta_Y(r)\), specialize \(f_1=y^r\). Each column in \(Q_Y(r)\)
then maps to a different row in (CP2), in the same order, so the matrix
is the identity. The analogous specialization for \(Z\) is \(f_1=z^r\).
Thus each determinant has an integer value one and remains a nonzero
polynomial in every characteristic. Their product defines a nonempty
open set over every infinite field and over sufficiently large finite
extensions. Its total degree is
\[
                     t(a-2t)(a-2t+1).                       \tag{CP4}
\]
An extension of a finite field with size greater than (CP4) suffices by
the usual polynomial grid bound. When \(a=2t\), (CP3) is empty and no
pure outer or central lower columns remain.

For \(a>2t\), the first determinants are
\[
 \Delta_Y(t)=([y^t]f_1)^{t(a-2t)},\qquad
 \Delta_Z(t)=([z^t]f_1)^{t(a-2t)},                            \tag{CP5}
\]
so these two coefficients are nonzero
on (CP3). Nonvanishing of those coefficients alone does not imply all
of (CP3).

There is a common explicit point in characteristic zero and in every
characteristic \(p>2(a-2t)\). Put \(L=a-2t\) and choose
\[
             f_1=y^t(1+y)^L+z^t(1+z)^L.                     \tag{CP5a}
\]
For \(d=r-t\) and \(s=L-d=A-r\), the \(Y\) minor is \(t\) copies
of the univariate Toeplitz matrix
\(\bigl[\binom L{d+i-j}\bigr]_{0\le i,j<s}\). Its classical determinant, recorded in
[Krattenthaler, Section 2.3, (2.16)--(2.18)](https://www.mat.univie.ac.at/~slc/s/s42kratt.pdf#page=11),
is
\[
 D_s(d)=\prod_{j=0}^{s-1}
        \frac{(L+j)!\,j!}{(d+j)!\,(L-d+j)!}.                 \tag{CP5b}
\]
Thus \(\Delta_Y(r)=\Delta_Z(r)=D_s(d)^t\ne0\) in those
characteristics. One elementary verification of (CP5b) uses the
Desnanot--Jacobi recurrence
\(D_s(d)D_{s-2}(d)=D_{s-1}(d)^2-D_{s-1}(d-1)D_{s-1}(d+1)\),
with \(D_0=1\), \(D_1(d)=\binom Ld\), and the triangular boundary
values at \(d=0,L\). The displayed product satisfies this recurrence.
All factorial arguments here are at most \(2L-1\); reduction modulo
\(p>2L\) therefore preserves nonvanishing. The open-set argument
above still covers smaller characteristics by allowing finite extensions.
The determinant evaluation and recurrence are classical inputs.

## Outer-source exchanges have a short transverse displacement

An outer pure \(Y\)-column and an upper pure \(Z\)-column have labels
\[
 e=(0,Y,Z),\quad a+t\le Y<2a-t,\ 0\le Z<t,
 \qquad f=(1,R,S),\quad0\le R<t,\ t\le S<a.
\]
Their difference is \(q=f-e=(1,-b,c)\), with
\(b=Y-R>0\) and \(c=S-Z>0\).

**Theorem.** On the locus (CP3), if this pair is exchanged on an actual
positive wall, then
\[
                              c\le6t.                       \tag{CP6}
\]
The transposed assertion holds for an outer \(Z\)-column and an upper
\(Y\)-column. The constant six is not claimed to be optimal.

Normalize the common wall cost by \(w_z=1\), and write
\[
 w=(\lambda,u,1),\quad \lambda=bu-c>0,\quad
                   \tau=Yu+Z=\lambda+Ru+S.                  \tag{CP7}
\]
Suppose \(c>6t\). Then \(a>6t\) and \(S>2t\). The central pure
\(Z\)-columns with labels
\[
                  t\le y<2t,\qquad2t\le z<a
\]
span the entire upper \(Z\)-tail \(z^{2t}C\). Indeed their image is
multiplication by the unit \(f_{1,Z}/z^t\) on that rectangular module,
using (CP5). Their maximum cost is \((2t-1)u+a-1\).
If \(\tau\) exceeded this value, \(f\) would already be in the span
of strictly cheaper columns and could not be exchanged. Consequently
\[
 u\le\frac{a-1-Z}{Y-2t+1}
       \le\frac{a-1}{a-t+1}<2.                              \tag{CP8}
\]
On the other hand \(b<2a\), so positivity in (CP7) gives
\[
             u>c/b>3t/a>2t/(a-t).
\]
Set \(r=\lfloor2t/u\rfloor\). These inequalities imply
\(t\le r<A\) and \((r+1)u>2t\).

Every column of the basis \(P_Y(r)\cup Q_Y(r)\) is now strictly
cheaper than \(e\). For the plain rectangle its maximum cost is
\((a+t-r-1)u+2t-1\), and
\[
 \tau-\max_{Q_Y(r)}w\cdot v
      \ge(r+1)u+Z-2t+1>0.                                  \tag{CP9}
\]
For the upper rectangle its maximum cost is
\(\lambda+(t+r-1)u+t-1\), and
\[
 \begin{split}
 \tau-\max_{P_Y(r)}w\cdot v
  &= S-(t+r-1-R)u-t+1\\
  &\ge S-tu-ru-t+1>S-5t+1>0.                              \tag{CP10}
 \end{split}
\]
Here \(tu<2t\) and \(ru\le2t\). By (CP3) these cheap columns span
\(T_Y\), so they span the actual image of \(e\), independently of
\(f_0\). This contradicts its exchangeability and proves (CP6).
The argument uses strict cheaper-column spans, not a support relaxation.

## Only central source columns remain outside a quadratic budget

The pure lower \(Y\)-columns are exactly the two rectangles
\[
\begin{array}{c|c}
\text{outer source }y^a &[a+t,2a-t)\times[0,t),\\
\text{central source }y^tz^t &[2t,a)\times[t,2t),
\end{array}                                                  \tag{CP11}
\]
and their \(Z\) counterparts are transposes. Add to the bridge note's
known inventory \(\mathcal N\) all directions
\[
 (1,-b,c),\ (1,b,-c),\quad1\le b,c<2a,
 \qquad\min(b,c)\le8t.               \tag{CP12}
\]
The enlarged inventory \(\mathcal N_*\) still has \(O(at)=O(n)\)
elements. Every outer pure cross-arm exchange has its direction there
by (CP6). The actual residual rank and labels still satisfy the
[translated-layer shadow theorem](layer_shadow_budget.md), so each
direction has \(O(n\log(n+1))\) parallel pieces. The
[finite polytope sweep](cross_layer_direction_budget.md) therefore gives
\[
 g(I_\phi)\le O(n^2\log(n+1))+E_{\rm central}(I_\phi),        \tag{CP13}
\]
where the remaining actual cross-layer pieces, outside \(\mathcal N_*\),
can only exchange a central pure lower column with an upper column on
the opposite arm. Both nonzero transverse entries of a remaining direction
have magnitude greater than \(8t\). All parallel subdivisions are retained. The controlled
pieces have total primitive one-norm length
\(O(a n^2\log(n+1))\le O(n^3\log(n+1))\).

This is a uniform localization on the explicit central-image open set,
including arbitrary choices of the other two generator images. It proves
neither an upper bound for \(E_{\rm central}\) nor a whole-family exponent
improvement. The complement of (CP3) remains part of the original
all-quotient problem. The previous monic-box ceiling is retained.

## Long central-pair eligibility uses only two generator images

There is a further pointwise reduction, with a stronger short-coordinate
cutoff. For any pure lower \(Y\)-column \(e=(0,Y,Z)\) and upper
\(Z\)-column \(f=(1,R,S)\), put \(b=Y-R\), \(c=S-Z\) as above.
On (CP3), if \(c>8t\) and this pair is exchanged at a positive weight,
then its common cost satisfies
\[
                              \tau<a u.                    \tag{CP14}
\]
Suppose instead \(\tau\ge a u\). As in (CP8), the central \(Z\)-tail
would make \(f\) a loop unless
\(\tau\le(2t-1)u+a-1\). This and \(\tau\ge a u\) imply \(u<2\).
Also \(b<2a\), so \(u>c/b>4t/a\), and \(a>8t\). Set
\(r=t+\lfloor2t/u\rfloor\). Then
\[
 t\le r<t+a/2<a-t,\qquad (r+1-t)u>2t.
\]
The basis (CP1) is again strictly cheaper than \(e\), since
\[
 \begin{split}
 \tau-\max_{Q_Y(r)}w\cdot v
   &\ge(r+1-t)u-2t+1>0,\\
 \tau-\max_{P_Y(r)}w\cdot v
   &\ge S-tu-ru-t+1>S-7t+1>0.
 \end{split}
\]
Here \(ru\le tu+2t\) and \(tu<2t\). This proves (CP14).

For a central lower \(Y\)-column, its vector depends only on \(f_1\),
while the upper \(Z\)-column is an identity vector. Every column from
the outer generator \(y^a\) has cost at least \(a u\). Thus, at weights
avoiding extra projected-label ties, the rank criterion for exchanging
this long central pair is independent of \(f_0\): below that cost no
\(f_0\)-column is present, and at or above it the pair is inactive for
every \(f_0\). Its eligibility depends only on \(f_1,f_2\). The transposed
orientation depends only on \(f_1,f_0\).

This is an assertion about literal pair eligibility. Other, more expensive
columns can still change the rest of an optimal basis and subdivide a
persistent wall. It is not an assertion that the parallel piece count or
the whole fan is independent of either outer image.

## Verification and priority boundary

```bash
python3 -B scripts/run_hook_central_prefix_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker tests the identity specializations independently against
literal monomial multiplication, verifies the actual prefix minors at
fixed field points, reconstructs original commuting quotient columns,
and checks both cheaper-span certificates in (CP8)--(CP10) on nonvacuous
long-arm cases. Unit-only central images supply excluded-stratum controls.
It also compares exact cheaper-column quotient ranks after changing only
the outer image in (CP14), exercising the below-head case and both
inactivity certificates. All those sampled long pairs are inactive; this
is a control for the stated reduction, not evidence that long central
walls survive or are absent uniformly.
Finite checks support the uniform argument; they are not a growth law.

Determinant specialization, units in finite local modules, matroid greedy
closure and the linked shadow/sweep bounds are established inputs. The
polynomial grid bound retains the source in the
[generic-hook certificate](generic_hook_fans.md). Exact priority of this
combined localization is unassessed. It supplies no general coherent
convolution-rank theorem, unrestricted improvement, or completion of the
active extremal goal.
