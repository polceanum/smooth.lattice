# Superregular polynomial axes have quadratic fans at every degree

The dense middle-degree axis stratum has a quadratic full-fan bound
on an explicit nonempty coefficient-open locus. Its long inverse,
large terminal prefix gap and large filtered rational degree do not
prevent this bound. The essential extra condition is nonsingularity of
all structurally possible minors of the banded Toeplitz matrices.
Full coefficient support alone is insufficient.

Put \(A=a-t\ge t\ge2\), \(n=4at-2t^2\), and take zero outer
images with central image
\[
                  y^tW_Y(y)+z^tW_Z(z).
\]
Let \(D_Y,D_Z\in[0,A-1]\) be the polynomial degrees. For either
arm write \(W=\sum_{d=0}^D c_dq^d\), and form
\[
             T_W=(c_{r-s})_{0\le r,s<A},\qquad
             c_d=0\text{ outside }[0,D].                       \tag{SR1}
\]
Assume that every minor with sorted row and column indices
\(r_1<\cdots<r_m\), \(s_1<\cdots<s_m\), and
\(0\le r_j-s_j\le D\) is nonzero. Call this banded
superregularity. In particular \(c_0c_D\ne0\).
Then, uniformly in the two degrees, coefficients and field,
\[
                             g(I)=O(n^2).                     \tag{SR2}
\]
For fixed \(a,t,D_Y,D_Z\), one coefficient-independent complete
positive integer weight list has \(O(n^2)\) members of height
\(O(a^2)\). The entire minimum-profile function and its normal fan
depend only on these integers. They can be constructed without field
elimination once the hypotheses are known. No efficient procedure for
certifying all the superregular minors, optimal basis-output cost or
bit-complexity bound is asserted.

The locus is nonempty over every infinite field and over sufficiently
large finite extensions in every characteristic. It does not contain
every coefficient specialization. The final reduction below extends the
bound to generic fixed initial orders and zero arms. Arbitrary outer
images and general Toeplitz-degenerate axes remain outside this theorem.
The sharp unrestricted extremal target and exact strongest-current
novelty remain unresolved.

## Why the coefficient locus exists in every characteristic

For this equal-width interval pattern a perfect matching exists exactly
when the sorted diagonal is permitted. A crossing pair in a matching
can be uncrossed while both differences remain in \([0,D]\); repeated
uncrossing gives the sorted diagonal. Thus other minors are structurally
zero.

In a permitted minor specialize \(c_d=X^{d^2}\). Its sorted-diagonal
term has the unique smallest exponent, since
\[
 \sum_j(r_j-s_{\pi(j)})^2
 =\sum_jr_j^2+\sum_js_j^2-2\sum_jr_js_{\pi(j)}
\]
is strictly minimized at the identity permutation by the strict
rearrangement inequality. This term has coefficient one. Every permitted
minor is consequently a nonzero polynomial in every characteristic.
Their finite product defines a nonempty open locus. Over a finite prime
field it has a point over some finite extension; there is no claim that
the prime field itself is large enough.

Over \(\mathbb Q\), the explicit specialization \(c_d=2^{d^2}\)
works simultaneously: factor the least power of two from a minor; the
remaining integer is congruent to one modulo two. This supplies exact
integer fixtures, rather than a floating-point normality test.

Superregular Toeplitz matrices and their finite-extension existence
are classical; see [Gluesing-Luerssen--Rosenthal--Smarandache,
Definition 3.3, Lemma 3.8 and Theorem 3.9](https://math.nd.edu/assets/209984/01580796.pdf).
That source treats the lower-triangular pattern. The argument above
states explicitly the additional upper bandwidth used here.

## The entire one-factor prefix has a closed rank formula

Use upper columns \(q^r\) of price \(\theta+r\), lower columns
\(q^rW\), \(t\le r<A\), of price \(r\), and heads
\(h_i=q^iW\), \(0\le i<t\). At a threshold \(\rho\), let
\(u\) be the number of cheaper-or-equal upper columns and \(s\)
the number of cheaper-or-equal lower columns. The pure span includes
the first \(u\) identity rows. After these rows are contracted, append
the first \(i\) heads. The remaining Toeplitz columns have input indices
\([0,i)\cup[t,t+s)\) and output rows \([u,A)\).

Superregularity makes the rank equal to the largest interval matching.
Greedy matching in increasing input order first places
\[
 h=\max(0,\min(i,A-u,D+i-u))
\]
heads in consecutive rows starting at \(u\). The tail starts at
\(b=\max(t,u+h)\). Hence the rank of the full original prefix is
\[
 u+R_i(u,s),\qquad
 R_i=h+\max(0,\min(s,A-b,D+t+s-b)).                            \tag{SR3}
\]
The interval matching is justified here by actual nonzero minors. It
is not a replacement for rank at a degenerate coefficient point.

Adding head \(i\) is dependent precisely when
\[
                  u>D+i\quad\hbox{or}\quad
                         \min(u,D)+s+i\ge A.                 \tag{SR4}
\]
One can see this directly from (SR3). If the new head has an available
row, its first new matched row increases \(\min(u,D)+i\) by one.
It displaces a tail match precisely when that endpoint reaches
\(A-s\). If the head has no available row, either \(u>D+i\)
or the available output dimension has already reached the second
condition. These alternatives include the truncated high-degree case.

Let \(T_k(\theta)\) be the \(k\)-th smallest price in
\[
 \{\theta,\theta+1,\ldots,\theta+D-1\}
                   \cup\{t,t+1,\ldots,A-1\};                 \tag{SR5}
\]
take \(T_k=+\infty\) when the list is too short. The whole-prefix
births from (AB4) are exactly
\[
 \boxed{\quad\alpha_i(\theta)=
 \min\left\{\begin{array}{l}
       \theta+D+i\quad\text{if }D+i<A,\\
       T_{A-i}(\theta)
       \end{array}\right\}.                                 \tag{SR6}
\]
An unavailable first option is omitted. The first event is the first
condition of (SR4); the order statistic is the second. This proves the
formula for nonintegral \(\theta\), and continuity gives every real
parameter. In particular the first option must be omitted when its
upper index has been truncated; replacing it by a smaller upper index
would give a wrong birth.

## Pure bases and head births need no coefficient elimination

Every greedy pure basis is an upper prefix and a lower prefix. If a
column fails greedy acceptance, multiply its cheaper dependence by
\(q\); the next column of the same type also fails. The shifted
cheaper span is contained in the next cheaper span. Each accepted set
therefore has the form
\[
 U_0,\ldots,U_{r-1},\quad L_t,\ldots,L_{t+A-r-1}.
\]
The available values of \(r\) are exactly
\[
                 \{t,\ldots,\min(A,t+D)\}\cup\{A\}.          \tag{SR7}
\]
Indeed the remaining square Toeplitz minor has constant diagonal
displacement \(r-t\); it is permitted when \(r-t\le D\),
and strictly triangular with zero diagonal otherwise. The empty minor
at \(r=A\) is one.

Choose the value of \(r\) minimizing
\[
 r\theta+\frac{r(r-1)}2+(A-r)t+rac{(A-r)(A-r-1)}2.            \tag{SR8}
\]
Together with (SR6) this constructs every compiled price interval and
the full minimum-profile oracle from degrees alone. The
[two-flag tensor-prefix identity](hook_axis_birth_curves.md) then
retains every coupled head dependency. No matching assertion is made
for the full two-arm residue matrix. The
[ordered-gain merge](hook_axis_linear_oracle.md) still gives
\(O(t+\log A)\) scalar profile queries after compilation.

## Each birth is an endpoint case or an interior balance

Work off integral factor-price lines, and write a selected birth as
\(\alpha_i=\delta\theta+q\), with \(\delta\in\{0,1\}\).
If the first option of (SR6) wins, \(q-i=D\). Otherwise inspect
the order statistic's two prefix counts at its winning price.

For an upper winning price, a zero lower count gives
\(q-i=A-2i-1\), a full lower count gives \(q-i=t-2i-1\),
and a full upper count gives \(q-i=D-i-1\). For a lower winning
price, a full upper count gives \(q-i=A+t-D-2i-1\), and a full
lower count gives \(q-i=A-i-1\). A zero upper count cannot reach
rank \(A-i>A-t\). Thus every endpoint case satisfies
\[
          \operatorname{dist}(q-i,\{0,D,A-D,A\})\le2t.         \tag{SR9}
\]

In every other case both counts are positive and unsaturated. Put
\(k=\lfloor\theta\rfloor\). Their sum is \(A-i\). At an
upper winner this gives \(2q=A+t-i-k-2\); at a lower winner it
gives \(2q=A+t-i+k-1\). Consequently the normalized gain obeys
\[
       2(\alpha_i(\theta)-i)-\theta=A+t+\epsilon_i,
                         \qquad |\epsilon_i|\le4t.           \tag{SR10}
\]
Such an interior event also has \(|\theta|\le A\), directly
from its nonempty unsaturated upper and lower price ranges.

## Count the full interaction inventory

Normalize weights to \((\lambda,u,1)>0\), with
\(\theta=(\lambda-t)/u\), \(\phi=\lambda-tu\).
Besides the two factor-price pencils, every additional candidate wall
has a positive equality
\[
             u(\alpha_i(\theta)-i)=\beta_j(\phi)-j.           \tag{SR11}
\]
If either birth is an endpoint case, one of its two long comparison
coordinates is within \(4t\) of one of the four centers in (SR9)
for that arm, allowing the normal's orientation and its short shift
by \(t\). Original comparison heights are at most \(a\).
These strips supply \(O(At)\) integer planes.

If both births are interior, substitute (SR10) in (SR11) and use
\(u\theta-\phi=t(u-1)\). Then
\[
          (A+2t+\epsilon_Y)u=A+2t+\epsilon_Z.                 \tag{SR12}
\]
For \(A>32t\), this implies \(|u-1|<10t/A\). For mixed
upper/lower labels write the raw comparison as \((1,b,-c)\) or
\((1,-b,c)\). Using \(\phi-\theta=(u-1)(\theta-t)\),
(SR10) gives
\[
                              |b+c-a|\le64t.                  \tag{SR13}
\]
For equal label types, the raw normal is \((0,b,-c)\), with
\(0<b,c\le a\), and (SR12) instead gives \(|b-c|\le64t\).
Each diagonal band again has \(O(At)\) integer normals. If
\(A\le32t\), the entire comparison box already has
\(O(A^2)=O(At)\) planes, so no division estimate is needed.

Count raw original differences before primitive normalization. Add all
factor and gain-zero pencils, only \(O(A)\) further planes. We have
obtained \(O(At)\) planes containing every actual wall. Their complete
positive arrangement has \(O(A^2t^2)=O(n^2)\) chambers, proving
(SR2). This retains intersections, parallel pieces and recurrent
supporting directions. It uses no noncrossing assumption.

All normals have height \(O(a)\). Primitive cross products and
three independent extreme-ray sums give the common integer list of
height \(O(a^2)\), by the established arrangement construction.

## Fixed initial orders and zero axes

The same quadratic theorem holds for a polynomial supported on an
interval \([\nu,D]\), \(0\le\nu\le D<A\), when every minor
whose sorted diagonal differences lie in \([\nu,D]\) is nonzero.
This is again a nonempty open locus in the fixed interval stratum:
apply the same least-exponent argument with
\(c_d=X^{(d-\nu)^2}\). The common integer list now fixes both
initial orders as well as both degrees.

Write \(W=q^\nu P\), \(B=A-\nu\). Upper rows below \(\nu\)
are coloops of the pure factor and are unused by any head. Remove them,
and translate the other row indices by \(\nu\). The effective
upper parameter is \(\theta'=\theta+\nu\), the polynomial degree
is \(D-\nu\), and the row dimension is \(B\).

If \(B\ge t\), apply the one-factor formulas with \((B,t)\).
If \(B<t\), every pure lower column is zero; heads with index
at least \(B\) are also zero. Apply the formulas with \((B,B)\)
to the remaining heads, and retain the absent suffix. Map every upper
birth label back by adding \(\nu\). Only the one-factor labels are
transferred: the original replicated profile still uses the original
short size \(t\).

For \(B\ge t\), an interior birth satisfies
\[
 2(\alpha_i-i)-\theta=B+t+\nu+\epsilon_i=A+t+\epsilon_i.
\]
Thus the interior interaction proof is unchanged. The endpoint labels
are within \(O(t)\) of
\(\{0,\nu,D,A-D,A-\nu,A\}\), still a fixed number of strips.
For \(B<t\) all finite upper birth labels lie in a band of width
less than \(t\); every birth is already an endpoint case.
Both the count and degree-and-order-only oracle follow.

For a zero arm, all its lower columns and head projections are zero;
its upper columns contribute a fixed affine profile. There are no
two-arm gain equalities involving that arm. The remaining pencils
also obey the claimed quadratic count. These reductions do not assert
that nongeneric nonunit coefficients satisfy superregularity.

## Scope and verification

~~~bash
python3 -B scripts/run_hook_superregular_axis_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The checker exhausts all minors at stated small sizes, including
structurally zero ones, and independently tests full pure-plus-head
prefix ranks. Whole factor tables agree with fresh field elimination.
Entire positive interaction segments are checked against the linear
inventory. Original quotient replay and exact whole-fan comparisons
check the degree-only construction separately from its scalar formulas.

An explicit characteristic-two polynomial with all coefficients nonzero
has a vanished admissible minor and a birth different from (SR6).
This fixture is excluded by the hypothesis. Generic nonunit and zero
controls separately verify the initial-order reduction with original
labels. Arbitrary coefficient degenerations remain open, as do
nonzero outer maps for this degree-independent theorem.

Superregularity, structural rank and interval matching are classical
inputs. The explicit whole-prefix specialization, degree-only oracle
and bounded-width interaction inventory are derived here. Their exact
strongest-current priority is unestablished. This generic polynomial
stratum does not settle sharp unrestricted \(G_3\), produce an infinite
actual polynomial excess or complete the original goal.
