# Thin radical point cones with nearly full compulsory cores

The [earlier sharp stability theorem](initial_cone_stability.md)
already disproves pointwise thickness as a route to the unrestricted
count: all-field local quotients have margin \(\Theta(n^{-3})\)
and logarithmic area \(\Theta(n^{-4})\).
Here the obstruction is extended to an explicit radical
characteristic-zero point family and simultaneous intermediate-cost cores.

Over every field of characteristic two, the family below has a full
positive initial cone whose logarithmic area is \(\Theta(n^{-3})\).
Its best relative inequality margin is \(\Theta(n^{-2})\).
All three of its walls have costs proportional to \(n\), and the
union of their compulsory cores has size \(n-2\).
A radical characteristic-zero point ideal of the same length has
exactly the same whole fan, so the obstruction also holds for
distinct points.

This gives a distinct-point version of the known obstruction, not
stronger worst-case stability exponents, a polynomial lower separation
or an improved unrestricted \(G_3(n)\) bound. No whole-fan asymptotic
count or novelty is established.

## An actual cyclic quotient

Let \(j\ge0\), and put

\[
 h=16j+4,\quad L=8h,\quad k=5h+1,\quad n=L+10.
\]

Over a field \(K\) of characteristic two take

\[
 I_h=\big((x^L-1)(x-1)^{10},\ xy-1,\ z-y^k\big).
                                                               \tag{T1}
\]

The quotient is \(K[x]/((x^L-1)(x-1)^{10})\), with \(y=x^{-1}\)
and \(z=x^{-k}\). The inverse exists because the univariate
polynomial has nonzero constant term. Its actual length is \(n\).
For \(h=4\) this is the local quotient with \((x-1)^{42}=0\);
larger \(h\) need not give a local algebra.

The marked standard set is

\[
 \begin{split}
 B_h={}&\{1,x,\ldots,x^{5h+3}\}\\
 &{}\cup\{y,\ldots,y^{3h}\}\\
 &{}\cup\{z,xz,x^2z,yz,y^2z,y^3z\}.                       \tag{T2}
 \end{split}
\]

Its cardinality is \(8h+10=n\). The minimal complement has heads

\[
 xy,\quad x^{5h+4},\quad y^{3h+1},\quad
 x^3z,\quad y^4z,\quad z^2.                              \tag{T3}
\]

## Six explicit relations

In characteristic two the following polynomials belong to \(I_h\):

\[
\begin{aligned}
g_0={}&xy+1,\\
g_1={}&x^{5h+4}+y^3z+x^{3h-4}+y^{3h-4},\\
g_2={}&y^{3h+1}+x^{5h-1}+z+x^{3h-1}+x^2z
                +x^{3h+1}+y^{3h-1}+x^{5h+1},\\
g_3={}&x^3z+x^{3h+2}+y^{3h}+x^{5h}
                +y^{3h-2}+x^{5h+2}+xz+x^{3h},\\
g_4={}&y^4z+x^{3h-5}+x^{5h+3}+y^{3h-3},\\
g_5={}&z^2+y^{2h+2}+x^{5h+2}+y^{3h-2}.                   \tag{T4}
\end{aligned}
\]

At \(h=4\), the two displayed \(y\)-terms in \(g_5\) cancel.
Every remaining tail in (T4) belongs to (T2).

For a direct membership proof put \(\varepsilon=x^L-1\).
Since \(L\) is divisible by \(32\), \((x-1)^{10}\mid x^L-1\).
In the quotient,

\[
 \varepsilon^2=0,\qquad
 x^{r+tL}=x^r(1+t\varepsilon),\qquad
 \varepsilon a(x)=0\ \Longleftrightarrow\
                      (x-1)^{10}\mid a(x).              \tag{T5}
\]

Here \(t\) is any integer, reduced modulo two in the coefficient.
Multiplication by \(x^{16}\) acts as the identity modulo
\((x-1)^{10}\), and \(h\equiv4\pmod {16}\).
Reducing each term of (T4) by (T5), the ordinary residue terms
cancel in pairs. The two nontrivial residual identities reduce to

\[
 \begin{split}
 x^3+x^5+x^{11}+x^{13}
       &=x^3(x+1)^{10},\\
 x^4+x^6+x^{12}+x^{14}
       &=x^4(x+1)^{10}.                                 \tag{T6}
 \end{split}
\]

The other residual pairs have equal exponents modulo \(16\).
This proves membership uniformly in \(j\), without extrapolating
from the case \(h=4\).

There is also a useful independent rank description. Modulo
\(\varepsilon\), the monomials in (T2) represent every residue
class of \(x\) modulo \(L\). Exactly ten classes occur twice:

\[
 \{5h,\ldots,5h+3\}\ \cup\
                    \{3h-4,\ldots,3h+1\}.              \tag{T7}
\]

The pair differences are \(\varepsilon x^a\) at these ten labels.
Their exponents modulo \(16\) are \(4,\ldots,13\); these ten
consecutive powers form a basis modulo \((x-1)^{10}\).
Thus (T2) is an actual vector basis as well.

## The full positive cone

Write a positive weight as \((s,t,u)\), and set

\[
 b=3h+1,\quad c=5h-1,\quad d=5h+4.
\]

The full open cone for (T2) is exactly

\[
 C_h=\{bt>2s+u,\quad u>cs,\quad ds>3t+u\}.                \tag{T8}
\]

These inequalities imply

\[
 \frac35<\frac st<\frac{3h+1}{5h+1}<1,\qquad u>cs.       \tag{T9}
\]

They make the first term of every row of (T4) its head. The potentially
largest tails can be checked as follows:

* For \(g_1\), (T8) handles \(y^3z\); the pure \(x\)-tail is
  smaller by its exponent, and \(ds>(3h-4)t\) follows from \(s/t>3/5\).
* For \(g_2\), the first two inequalities give
  \(bt>(5h+1)s\); the remaining pure \(x\)-tails are no larger,
  while \(y^{3h-1}\) is smaller and \(x^2z\) is the defining comparison.
* For \(g_3\), \(3s+u>(5h+2)s\) handles every pure \(x\)-tail,
  and this is greater than \(3ht\) by (T9). The tail \(xz\)
  is smaller by \(2s\).
* For \(g_4\), \(4t+u-(5h+3)s>4(t-s)>0\);
  its largest pure \(y\)-tail is smaller since
  \(u/t>3h-3/5\).
* For \(g_5\), \(2u>(10h-2)s>(5h+2)s\), and
  \(2u/t>6h-6/5>3h-2\). At \(h\ge4\) the largest displayed
  pure \(y\)-exponent is \(3h-2\).

The heads (T3) leave exactly \(n\) standard monomials. Membership,
this head inclusion, and equality of finite colengths prove that
(T4), after cancellation, is the reduced Gröbner basis.
Conversely its comparisons with \(y^3z,x^2z,x^{5h+2}\) force
all three inequalities in (T8). Hence this is the whole cone,
not a narrow subcone imposed by redundant comparisons.

Its primitive inward facet normals are

\[
 q_1=(-2,b,-1),\quad q_2=(-c,0,1),\quad q_3=(d,-3,-1),
 \qquad \det(q_1,q_2,q_3)=2.                             \tag{T10}
\]

All its closed extreme rays are strictly positive:

\[
 \begin{split}
 r_1&=(b,5h+1,bc),\\
 r_2&=(3,5,3c),\\
 r_3&=\tfrac12(b+3,5h+6,db-6).
                                                               \tag{T11}
 \end{split}
\]

These are primitive. In particular the cone does not reach the
boundary of the positive weight simplex.

## Exact height and relative margin

Every strict integer weight has integer slacks
\(q_i\cdot w\ge1\). The inverse of the row matrix in (T10) is
entrywise positive. Setting all three slacks equal to one gives
the coordinatewise smallest strict integer weight

\[
 w_*=(3h+4,\ 5h+6,\ 15h^2+17h-3).
\]

Therefore

\[
 \min_{w\in C_h\cap\mathbb Z^3}\|w\|_\infty
                       =15h^2+17h-3=\Theta(n^2).         \tag{T12}
\]

For a cone with these irredundant normals define its relative margin by

\[
 \delta(C_h)=\sup_{w\in C_h}
          \min_i\frac{q_i\cdot w}{|q_i|\cdot w}.
                                                               \tag{T13}
\]

Absolute value is coordinatewise. This margin controls independent
multiplicative perturbations of all three prices. The weight \(w_*\)
gives

\[
 \delta(C_h)\ge\frac1{30h^2+64h+31}.
\]

For the other direction (T8)--(T9) give

\[
 0<u-cs<5s-3t<2s/b.
\]

Consequently the \(q_2\) margin is smaller than \(1/(bc)\), and

\[
 \frac1{30h^2+64h+31}
       \le\delta(C_h)<\frac1{(3h+1)(5h-1)}
       =O(n^{-2}).                                     \tag{T14}
\]

Thus no uniform lower bound \(a/n\), with fixed \(a>0\), holds
for every actual initial cone.

## Logarithmic area

Take ratios \(p=s/t,\ q=u/t\). The cone section is the open triangle
with vertices

\[
 \left(\frac35,\frac{3c}5\right),\quad
 \left(\frac b{5h+1},\frac{bc}{5h+1}\right),\quad
 \left(\frac{b+3}{5h+6},\frac{db-6}{5h+6}\right).
\]

Its ordinary area is exactly

\[
 A_h=\frac2{5(5h+1)(5h+6)}.                              \tag{T15}
\]

Under \((p,q)\mapsto(\log p,\log q)\) the area density is
\(1/(pq)\). Throughout this triangle \(h<pq<3h\). Its logarithmic
area \(\mathcal A_h\) therefore obeys

\[
 \frac{A_h}{3h}<\mathcal A_h<\frac{A_h}{h},
 \qquad \mathcal A_h=\Theta(n^{-3}).                    \tag{T16}
\]

This rules out a uniform \(\Omega(n^{-2})\) logarithmic area per
cone, even when every coordinate ratio in the cone has polynomial
size. Counting cones requires more than applying that proposed
pointwise lower bound.

## The cores nearly fill one basis

The [wall threshold](primitive_wall_spectrum.md) gives the three costs

\[
 \nu(q_1)=3h+6,\qquad \nu(q_2)=5h,\qquad
 \nu(q_3)=5h+11.                                       \tag{T17}
\]

All are fixed positive fractions of \(n\) asymptotically. In particular
they lie in the remaining intermediate-cost range for large \(h\).
Let \(C_{q_i}\) be the literal compulsory monomial cores. Directly,

\[
 C_{q_1}\cup C_{q_2}\cup C_{q_3}
                   =B_h\setminus\{x^2z,y^3z\}.           \tag{T18}
\]

Thus their union has rank \(n-2\) in this same actual quotient.
Small residual rank for several incident cores does not imply
relative margin of order \(1/n\), or logarithmic area of order
\(1/n^2\). This does not invalidate counting with shared cores;
it invalidates that pointwise thickness conclusion from them.
The union is asserted at the marked cone. It is not asserted to
survive in all its neighbors or to define an algebra quotient.

The three costs in (T17) do not contradict the conjectured collective
cubic polylogarithmic wall budget. Neither a single thin cone nor its
three walls yields a large fan.
The inverse-unit relation \(xy=1\) also places this family in the
unresolved single-bilinear-defect class. Its anchored \(x,y\) box
with side counts \(2h,2h\) has rank \(4h-1<n\) and cardinality
\(4h^2\): its Laurent powers run from \(-(2h-1)\) to \(2h-1\)
and are independent in the length-\(L\) group quotient. The
linear nonsaturating-box regime does not apply.

## Transfer to radical characteristic-zero points

The transfer here needs only a fixed ten-row Fourier calculation.
Let \(\omega\) be a primitive \(2L\)-th root of unity over
\(\mathbb Q\), and let \(\eta=\omega^h\), of order \(16\).
Put

\[
 E=\{0,2,4,6,8,10,12,14,1,9\}.
\]

Take the \(n=L+10\) distinct parameters

\[
 T=\{\omega^{2a}:0\le a<L\}
       \ \sqcup\ \{\omega\eta^e:e\in E\},                \tag{T19}
\]

and their points \((t,t^{-1},t^{-k})\) over
\(\mathbb Q(\omega)\). The portions are disjoint: the first has
\(t^L=1\), and the second has \(t^L=-1\).

For any subset of power labels modulo \(2L\), let \(U\) be the
set of occupied pairs \(\{a,a+L\}\), \(0\le a<L\), and \(D\)
the labels of the pairs in which both powers are present.
The characteristic-two quotient (T1) has column rank

\[
 |U|+\operatorname{rank}
       \left(\binom a\ell\bmod2\right)_
                         {0\le\ell<10,\ a\in D}.         \tag{T20}
\]

Indeed the first \(L\) residue coordinates distinguish the pairs,
and the differences within double pairs are
\(\varepsilon x^a\). For the radical evaluation matrix, invert
the full \(L\)-point Fourier block on the first portion of (T19).
It separates the same pairs. On the last ten rows, the difference
within a double pair is a nonzero scalar times
\((\eta^{ea})_{e\in E}\). Thus its rank is

\[
 |U|+\operatorname{rank}(\eta^{ea})_{e\in E,\ a\in D}.
                                                               \tag{T21}
\]

The two residual ranks agree for every list \(D\), including repeated
labels modulo \(16\). To see this, group the sixteen labels as
\(\{u,u+8\}\), \(0\le u<8\). The first eight Pascal rows, or the
eight even Fourier rows, separate these pairs. In a double pair,
the two remaining Pascal rows give \((1,u\bmod2)\); the two odd
Fourier rows give, up to a nonzero scalar,
\((1,(-1)^u)\). In either case the residual rank is the number
of parities represented among the double-pair labels. Hence both
matrices have the same rank formula

\[
 \#\{\text{occupied pairs modulo }8\}
       +\#\{\text{parities of double-pair labels}\}.      \tag{T22}
\]

Powers have period \(2L\) in both full algebras. Equations
(T20)--(T22) identify the labelled matroid of every Laurent
monomial list. Greedy selection therefore identifies every
standard set, and the entire positive fans agree.
The cone (T8), its thinness, exact integer height, costs and core
union all transfer without a length overhead.
This uses roots in a cyclotomic field; it is not a claim that
all the points are rational over \(\mathbb Q\).

## Verification and mathematical boundary

Run:

~~~bash
python3 -B scripts/run_thin_initial_cone_checks.py
~~~

The report is results/research/thin_initial_cone_checks.json.
It includes symbolic polynomial certificates valid for every
\(j\ge0\), as well as independent finite companion-division,
commuting-matrix, rank, S-pair, cone and core checks.
All \(2^{16}\) subsets of the residual sixteen-column matrices are
also compared over the binary and Fourier fields. Two direct radical
point fixtures independently check the marked basis and its full cone.
It does not infer the uniform statements from a finite numerical fit.

The cyclic quotient, square-zero extension, consecutive-power basis,
reduced-basis criterion, determinant rays and change-of-variables
area formula are classical. The
[earlier height obstruction](initial_cone_height.md) already proves
the sharp quadratic certificate-height exponent. The new distinction
here is an entirely interior positive cone with the stated relative
margin, logarithmic area and simultaneous intermediate-cost cores.
The exact quantitative priority of this combination is unassessed;
unsuccessful searches are not novelty evidence.
The Pascal/Fourier transfer mechanism is also recorded in the
[prime-digit lift theorem](prime_digit_lifts.md). Here the explicit
\(L+10\) row set handles periods \(2L\) that need not be powers of two.
The [Babson--Onn--Thomas framework](https://arxiv.org/pdf/math/0207135)
supplies the established connection between standard bases and
positive profile cones. This example does not change the unrestricted
extremal endpoint, its priority requirements, or the publication gate.
