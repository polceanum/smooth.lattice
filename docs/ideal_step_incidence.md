# Large incidence and wall-length jumps under deletion of one point

An ideal composition series does not support a uniform quadratic
increment bound for square-root facet incidence. This fails already
for the fixed rational curve
\[
 xy=1,\qquad (x-1)z=1.
\]
In the nested radical point ideals
\[
 I_n=(xy-1,(x-1)z-1,F_n(x)),\qquad
 F_n(T)=\prod_{j=2}^{n+1}(T-j)\quad\text{over }\mathbb Q,
                                                        \tag{1}
\]
one has \(I_n\subset I_{n-1}\), and the kernel of the quotient map
has dimension one. Nevertheless there are positive constants and
infinite arithmetic progressions of lengths for which
\[
 U(I_n)-U(I_{n-1})\ge c n^3,
 \qquad\text{or respectively}\qquad
 U(I_n)-U(I_{n-1})\le-c n^3.                              \tag{2}
\]
Here the established incidence quantity is
\[
 U(I)=\sum_{w>0\text{ primitive interior fan ray}}
                  \lceil\sqrt{\|w\|_\infty}\rceil d(w),
\]
where \(d(w)\) counts incident full-dimensional positive cones.
In particular even the one-sided positive increment is not uniformly
\(O(n^2\log^C(n+1))\) for any fixed \(C\).

By contrast, the primitive wall-length quantity
\(W(I)=\sum_e\|q_e\|_1\), counting each actual undirected wall
piece once, satisfies in this same family
\[
 \boxed{\quad W(I_n)-W(I_{n-1})
       =3\sum_{\substack{d\mid n\\d\ge2}}d\varphi(d)
       \le3n(n-1).\quad}                                \tag{3}
\]
Equality holds when \(n\) is prime. The failure in (2) alone therefore
does not refute a quadratic increment conjecture for \(W\). The product
construction below does refute it: there are actual radical ideal steps
with wall-length jumps of both signs of order \(n^3\), or
\(n^3\log n\). Suitable chains have quartic total variation.
These facts do not refute a cubic-polylogarithmic global budget for either
quantity. No such unrestricted budget is proved here.
The unrestricted extremal initial-ideal problem remains open.

## Exact fan and actual ideal step

The [three-pole fan theorem](weighted_facet_incidence.md#quartic-incidence-is-sharp-for-actual-finite-quotients)
shows that (1) has the standard sets
\[
 D_{a,b,c}=\{1,x,\ldots,x^a,y,\ldots,y^b,z,\ldots,z^c\},
       \qquad a+b+c=n-1.
\]
Every such allocation occurs, so \(g(I_n)=\binom{n+1}{2}\) and
\(g(I_n)-g(I_{n-1})=n\). Its interior rays are exactly
\[
 w_{A,B,C}=\frac{(BC,AC,AB)}{\gcd(BC,AC,AB)},\qquad
 A,B,C\ge1,\quad A+B+C\in\{n,n+1\}.                     \tag{4}
\]
Each has precisely three incident cells.

The ideal step removes the point with parameter \(n+1\).
The residue of \(F_{n-1}(x)\), normalized by its nonzero value there,
is its separator: it vanishes at all other points and equals one at
that point. Its span is an actual common coordinate eigenline and an
algebra ideal. Thus no arbitrary vector-space contraction or generic
matroid realization is used in (2).

The common rays of the length-\(n\) and length-\(n-1\) fans are
precisely the level \(A+B+C=n\) in (4). Proportional integer triples
on two levels have primitive coordinate sum dividing both levels.
Different levels among \(n-1,n,n+1\) have gcd at most two, whereas
a positive primitive triple has sum at least three. Hence the other
levels cannot overlap. The longer fan has \(\binom n2\) new rays,
and the shorter has \(\binom{n-2}2\) rays absent from the longer.
Counting newly positioned rays therefore discards substantial cancellation.

## The arithmetic mean of primitive normal heights

Put
\[
 H(A,B,C)=\frac{\max(AB,AC,BC)}{\gcd(AB,AC,BC)},\qquad
 S(N)=\sum_{\substack{A,B,C\ge1\\A+B+C=N}}
                      \lceil\sqrt{H(A,B,C)}\rceil.
\]
Then, exactly,
\[
 U(I_n)=3(S(n)+S(n+1)),\quad
 U(I_n)-U(I_{n-1})=3(S(n+1)-S(n-1)).                    \tag{5}
\]
The dependence of \(S(N)\) on the prime factors of \(N\) cannot be
replaced by a constant leading coefficient.

Let
\[
 J=\int_{\substack{a,b>0\\a+b<1}}
       \sqrt{\max(ab,a(1-a-b),b(1-a-b))}\,da\,db>0.
\]
For a prime \(p\), define
\[
 \begin{split}
 L_p(0)&=1-\frac{3(1-p^{-1/2})}{p^2(1-p^{-5/2})},\\
 L_p(e)&=1-p^{-2}+p^{-3}L_p(e-1),\quad e\ge1,\\
 L_p(\infty)&=\frac{1-p^{-2}}{1-p^{-3}}
                  =1-\frac1{p^2+p+1},\\
 \mathcal C(N)&=\prod_p L_p(v_p(N)).
 \end{split}                                             \tag{6}
\]
All factors are positive, increase with \(e\), and lie between
\(L_p(0)\) and \(L_p(\infty)\). The products converge uniformly
in \(N\), and their common positive lower bound is
\(\mathcal C_0=\prod_p L_p(0)\).

The following asymptotic holds along the entire sequence of integers:
\[
                   S(N)=J\mathcal C(N)N^3+o(N^3).         \tag{7}
\]
The error divided by \(N^3\) tends to zero without restricting the
factorization of \(N\). In particular \(U(I_n)=\Theta(n^3)\)
for all sufficiently large lengths, including the oscillating subsequences.

Here is a direct proof of (6)--(7). If the three \(p\)-adic valuations
are sorted, the valuation of \(\gcd(AB,AC,BC)\) is the sum of the
two smallest. When \(p\nmid N\), exactly two entries are divisible
by \(p^k\) with minimum valuation exactly \(k\ge1\) with probability
\(3(1-p^{-2})p^{-2k}\), under uniform residues of \(A,B\) and
\(C=N-A-B\). Otherwise the gcd valuation is zero, with probability
\(1-3p^{-2}\). Summing their weights \(p^{-k/2}\) gives \(L_p(0)\).
When \(p\mid N\), either not all three entries are divisible by
\(p\), giving gcd valuation zero, or all three are, with probability
\(p^{-2}\). Dividing the entries by \(p\) divides the reciprocal
square-root gcd weight by \(p\) and replaces \(N\) by \(N/p\).
This proves the recurrence. It also proves convergence to the displayed
limit. The strict inequality \(L_p(0)<L_p(\infty)\) follows directly
from the first formula; for \(p\ge3\), use
\(3(1-p^{-1/2})>1\), and check \(p=2\) separately.

For completeness the passage from local averages to (7) needs no
unproved independence of integer divisibility events. Discard the
ceiling, with error at most \(\binom{N-1}2=O(N^2)\).
Keep only primes \(p\le B\), and cap the individual valuations at
\(K\). This is a periodic weight in \(A,B,N\), of fixed modulus
\(\prod_{p\le B}p^K\). The Chinese remainder theorem factors its
residue average into the finite versions of (6). On every residue
class, the ordinary Riemann sum of the continuous bounded integrand
in \(J\) has the same limit. There are finitely many classes, so
this limit is uniform in the residue of \(N\).

The omitted weights have a uniform small tail. For any integer \(d\),
the number of triples with \(d\) dividing a specified pair is exactly
\[
                   \binom{\lfloor(N-1)/d\rfloor}{2}
                         \le\frac{N^2}{2d^2}.           \tag{8}
\]
A prime can affect the gcd only if it divides at least two entries.
Since \(\sqrt{\max(AB,AC,BC)}\le N/2\), the discarded large-prime
contribution is at most
\[
                \frac{3N^3}{4}\sum_{p>B}p^{-2}.
\]
Capping valuations adds at most
\(3N^3\sum_{p\le B}p^{-2K}/4\). First let \(N\) grow for fixed
\(B,K\), then let \(K\) and \(B\) grow. The identical tail bounds
apply to local averages. This proves (7), including its uniformity.

## Infinite cubic jumps of both signs

Let \(Q\) be the product of the odd primes at most \(101\), and
take \(N=Q(2k+1)\). Both \(N\) and \(N\pm2\) are odd.
Every odd prime at most \(101\) divides \(N\), and none divides
\(N\pm2\). In particular
\[
 \mathcal C(N)\ge\mathcal C_0\frac{L_3(1)}{L_3(0)}
                     >\frac{29}{27}\mathcal C_0.         \tag{9}
\]
For the last elementary inequality, the formula for \(L_3(0)\)
and \(1/\sqrt3<3/5\) give \(L_3(0)<6/7\); the recurrence
then gives \(L_3(1)/L_3(0)>29/27\).

On the other hand, for \(p\ge3\),
\[
 \log\frac{L_p(\infty)}{L_p(0)}
 \le-\log(1-3/p^2)\le\frac{9}{2p^2}.
\]
It follows that
\[
 \mathcal C(N\pm2)\le\mathcal C_0\exp(9/200)
                 \le\frac{200}{191}\mathcal C_0
                 <\frac{29}{27}\mathcal C_0.             \tag{10}
\]
Here \(\sum_{p>101}p^{-2}\le1/100\), and
\(e^x\le(1-x)^{-1}\) for \(0<x<1\).

Equations (7), (9), (10) show
\(S(N)-S(N-2)=\Omega(N^3)\) and
\(S(N+2)-S(N)=-\Omega(N^3)\). Apply (5) at lengths \(N-1\)
and \(N+1\), respectively, to prove (2). The universal quadratic
height bound and \(\Theta(n^2)\) rays also give an \(O(n^3)\)
upper bound on the absolute jump, so these jump orders are sharp in
this family. A small finite instance is
\(U(I_{104})-U(I_{103})=71364\), whereas
\(U(I_{106})-U(I_{105})=-28866\). These values illustrate the
formula, rather than prove its asymptotic consequence.

## Primitive wall length has a different recurrence

For each unordered pair of axes, walls exchange \(x_i^A\) with
\(x_j^B\), where \(A,B\ge1\), \(A+B\le n\). The third-axis
allocation is \(n-A-B\). Each pair gives one wall piece, with
primitive one-norm \((A+B)/\gcd(A,B)\). Thus
\[
 W(I_n)=3\sum_{\substack{A,B\ge1\\A+B\le n}}
                   \frac{A+B}{\gcd(A,B)}
       =3\sum_{s=2}^{n}\sum_{\substack{d\mid s\\d\ge2}}
                   d\varphi(d).                         \tag{11}
\]
The second equality groups by \(s=A+B\) and
\(d=s/\gcd(A,s)\). Taking a difference proves (3), since
\(\sum_{d\mid n,d\ge2}\varphi(d)=n-1\) and \(d\le n\).
The prime case proves sharpness of its increment constant.

One may also sum (11):
\[
 W(I_n)=\frac{\zeta(3)}{\zeta(2)}n^3+O(n^2\log(n+1)).    \tag{12}
\]
Indeed Möbius inversion gives
\(\sum_{d\le X}d\varphi(d)=X^3/(3\zeta(2))+O(X^2\log(X+1))\).
Sum this at \(X=n/k\), and multiply by three. Removing the omitted
divisor \(d=1\) changes the result only by \(3n\).
All these are elementary classical totient manipulations.

The distinction is operational for the unrestricted proof obligation.
An induction paying a quadratic polylogarithmic amount for the positive
change of \(U\) at every ideal step is false, even on a fixed smooth
curve. A successful global cubic \(U\) proof needs a more global
charge or cancellation. The sufficient primitive-length budget
\(U\le2W\) remains available; (3) checks one sharp family but proves
no bound for general adjacent ideals. The original unrestricted
\(G_3\) endpoint and actual polynomial-separation alternative are unchanged.

## Primitive wall length also has large one-point jumps

The quadratic recurrence (3) is specific to the three-pole family. It
cannot be used as a uniform induction hypothesis. The existing
[common-coordinate product](wall_cost_products.md#the-exact-wall-cost-of-a-common-coordinate-fan-product)
transfers a planar wall-length oscillation to a single actual ideal step.

Work over \(\mathbb Q\). Write
\[
 f_m(x)=\prod_{j=1}^m(x-j),\qquad
 T_m=(f_m(x),xy-1),\qquad
 a_m=W(T_m)=\sum_{d\mid m}d\varphi(d)-1.                 \tag{13}
\]
In particular \(a_1=0\), \(g(T_m)=m\) and
\(0\le a_m\le m(m-1)\). Fix \(r\ge2\) and a planar radical
point ideal \(J_r\subset\mathbb Q[x,z]\) of length \(r\), containing
\(f_r\), with all its point coordinates \(z\ne0\). Set
\(b_r=g(J_r)\), \(c_r=W(J_r)\), and, for \(1\le m\le r\),
\[
 H_{m,r}=(f_r,yz,yT_m,zJ_r),\qquad n_{m,r}=2r+m.
                                                               \tag{14}
\]
The all-order product formula gives
\[
 g(H_{m,r})=m b_r,\qquad
 W(H_{m,r})=V(H_{m,r})=b_r a_m+m c_r.                    \tag{15}
\]
Every wall direction has a zero coordinate, so primitive one-norm and
the divisor threshold coincide. Each factor wall has a distinct piece
for every cone of the other factor; these multiplicities must be retained.

The points of (14) are the \(r\) axis points \((j,0,0)\), the
\(m\) points \((j,j^{-1},0)\), and the \(r\) lifted points of
\(J_r\) with \(y=0\). They are distinct. The product length equals
their cardinality, proving that (14) is their radical ideal.
Moreover \(H_{m,r}\subset H_{m-1,r}\) for \(m\ge2\): only
\((m,m^{-1},0)\) is deleted. A normalized separator is
\[
 s_m=\frac{m\,y f_{m-1}(x)}{f_{m-1}(m)}.
\]
It is one at that point and zero at every other point. Its span is
an actual ideal eigenline, not a vector-space contraction. Thus
\[
 W(H_{m,r})-W(H_{m-1,r})
                  =b_r(a_m-a_{m-1})+c_r.                \tag{16}
\]
More generally, for any nested planar branch with fixed common annihilator,
the product increment is \(b_r\,\Delta W_{\rm branch}
+c_r\,\Delta g_{\rm branch}\). The hook branch has \(\Delta g=1\)
and makes the arithmetic in (16) explicit.

### Both signs have cubic size, already on one fixed reduced curve

Put \(F(m)=a_m+1\). This is the classical multiplicative divisor sum
\(\sum_{d\mid m}d\varphi(d)\), with
\[
 F(\ell^e)=\frac{\ell^{2e+1}+1}{\ell+1}
                         \quad(\ell\text{ prime}).
\]
Since \(F(q)/q^2\le1\) and
\(F(2^e)/2^{2e}=2/3+1/(3\cdot2^{2e})\le3/4\) for \(e\ge1\),
every even \(q\) satisfies \(a_q\le3q^2/4-1\).
For each odd prime \(p\), therefore,
\[
 a_p-a_{p-1}\ge\frac{(p+1)^2}{4},\qquad
 a_{p+1}-a_p\le-\frac{p^2}{4}+\frac{5p}{2}-\frac14
                         \le-\frac{p^2}{8}\quad(p\ge23).
                                                               \tag{17}
\]

First take \(J_r=(f_r,xz-1)\). Then \(b_r=r,c_r=a_r\), and every
point of every ideal (14) lies on the single fixed reduced curve
\[
 \{y=z=0\}\ \cup\ \{z=0,xy=1\}\ \cup\ \{y=0,xz=1\}.
                                                               \tag{18}
\]
At \(r=m=p\), (16)--(17) give a positive jump at least
\(p(p+1)^2/4\). At \(r=m=p+1\), they give a negative jump at most
\(-(p+1)p^2/16\) for \(p\ge23\), using
\(c_r\le r(r-1)\). Both are of cubic order in \(n_{m,r}=3r\).
All individual wall costs in this fixed-curve family satisfy
\(W(H_{m,r})\le2r^3\). Large signed jumps are therefore compatible
with its cubic global budget and quadratic fan count.

For a logarithmically larger example take
\[
                J_r=(f_r,z-x^r).
\]
Its positive nodes give the generic planar corner-cut fan: on \(V_r^2\)
the exponents \(i+rj\) are distinct, and the corresponding generalized
Vandermonde minors do not vanish. Hence
\(b_r=C_2(r)=\Theta(r\log(r+1))\), while
\(0\le c_r\le r(r-1)\). The same two prime choices in (17)
give positive and negative jumps of order \(n^3\log(n+1)\).
The upper order follows directly from (16) and \(a_m\le m(m-1)\).
These varying-degree embeddings are not asserted to lie on a fixed curve.

Consequently **no** uniform one-step positive or absolute bound
\(O(n^2\log^C(n+1))\), for fixed \(C\), holds for \(W\), even for
radical rational point ideals. Primes are used only to select easy
large oscillations in a classical divisor sum; this gives no new result
about the distribution of primes.

### Quartic variation along actual deletion chains

Even replacing a stepwise estimate by a uniform total-positive-increment
estimate along arbitrary ideal composition chains is insufficient.
The same family gives an order-one-power loss.

The elementary divisor-sum averages are
\[
 \begin{split}
 A(X):=\sum_{m\le X}a_m
       &=\frac{\zeta(3)}{3\zeta(2)}X^3+O(X^2\log(X+1)),\\
 A_{\mathrm{odd}}(X):=\sum_{\substack{m\le X\\m\text{ odd}}}a_m
       &=\frac{7\zeta(3)}{36\zeta(2)}X^3+O(X^2\log(X+1)).
 \end{split}                                             \tag{19}
\]
Here is a direct proof of the parity constant, without a prime-density
theorem or an inference from numerical data. Since
\(d\varphi(d)=\sum_{e\mid d}\mu(e)d^2/e\), summing first over odd
\(d\le Y\) gives
\[
 \begin{split}
 B_{\mathrm{odd}}(Y):=\sum_{\substack{d\le Y\\d\text{ odd}}}d\varphi(d)
 &=\sum_{\substack{e\le Y\\e\text{ odd}}}\mu(e)e
                 \sum_{\substack{k\le Y/e\\k\text{ odd}}}k^2\\
 &=\frac{Y^3}{6}\sum_{e\text{ odd}}\frac{\mu(e)}{e^2}
                      +O(Y^2\log(Y+1))\\
 &=\frac{2Y^3}{9\zeta(2)}+O(Y^2\log(Y+1)).
 \end{split}
\]
The sum over odd squares is \(Y^3/6+O(Y^2)\). The absolutely convergent
Euler product is \(4/(3\zeta(2))\); extending its finite sum contributes
only \(O(Y^2)\). Now sum \(B_{\mathrm{odd}}(X/k)\) over odd \(k\).
The leading multiplier is \(\sum_{k\text{ odd}}k^{-3}=7\zeta(3)/8\).
Subtracting the omitted divisor-one terms changes the result by \(O(X)\).
This proves the second formula in (19); the unrestricted version proves
the first and agrees with (12).

For even \(r\), put
\[
 L_r=\sum_{j=1}^{r/2}(a_{2j-1}-a_{2j})
    =2A_{\mathrm{odd}}(r)-A(r)
    =\kappa r^3+O(r^2\log(r+1)),\quad
                \kappa=\frac{\zeta(3)}{18\zeta(2)}>0.     \tag{20}
\]
Write \(\Delta_m=W(H_{m,r})-W(H_{m-1,r})\), \(2\le m\le r\),
and let \(P_r=\sum\max(\Delta_m,0)\),
\(N_r=\sum\max(-\Delta_m,0)\). Summing just the even steps gives
\[
 N_r\ge-\sum_{m\text{ even}}\Delta_m
        =b_r L_r-\frac r2c_r.                            \tag{21}
\]
For either choice of \(J_r\) above this is \(\Omega(b_r r^3)\).
Also \(P_r-N_r=W(H_{r,r})-W(H_{1,r})\ge0\), so \(P_r\ge N_r\).
Finally \(\sum_{m=2}^r|a_m-a_{m-1}|\le2A(r)=O(r^3)\), and (16)
gives the matching upper bound. Thus
\[
 \boxed{\quad P_r,N_r,\ \sum_{m=2}^r|\Delta_m|
                           =\Theta(b_r r^3).\quad}       \tag{22}
\]
For even \(r\to\infty\), these are actual consecutive point-deletion steps of lengths between
\(2r+1\) and \(3r\). With maximum length \(N=3r\), their total
positive, total negative and absolute variations are \(\Theta(N^4)\)
on the fixed curve (18), and \(\Theta(N^4\log(N+1))\) for the generic
branch. The corresponding largest individual \(W\) has orders
\(N^3\) and \(N^3\log(N+1)\): (19) supplies an \(m\le r\) with
\(a_m=\Omega(r^2)\). The variation loses a sharp factor of order \(N\)
on these families. Extending a displayed chain by further point deletions
gives a full ideal composition series and preserves its variation lower bound.

This rules out a uniform cubic-polylogarithmic bound on total positive
increments over **every** ideal chain. It does not rule out a carefully
chosen chain with less variation for each ideal, a cancellation argument,
or the global cubic-polylogarithmic \(W\) budget. The product class already
has that global bound and stays below the generic fan-count scale.
The unrestricted \(G_3\) question and exact priority remain open.

## Actual edge lengths: planar cancellation and a spatial obstruction

Keeping integer edge multiples removes the arithmetic oscillation in
(13), but does not restore a uniform quadratic-step induction. Define
\[
 \mathcal L(I)=\sum_{e=\{B,C\}}
                         \|\mu(B)-\mu(C)\|_1,           \tag{23}
\]
over actual undirected positive walls. Thus \(W\le\mathcal L\).
The [visible skeleton](weighted_facet_incidence.md#square-root-incidence-and-the-visible-edge-skeleton)
\(\ell_+\) counts only edges incident to strictly positive two-faces;
\(\ell_+\le\mathcal L\), and either cubic-polylogarithmic bound would
suffice for the existing incidence reduction. Neither is assumed here.

### A sharp planar one-step bound

For a planar length-\(m\) quotient \(J\), let \(B^0,B^\infty\)
be its two limiting positive standard sets as \(w_x/w_y\) tends to
zero and infinity. Monotonicity of the two profile coordinates gives
\[
 \mathcal L(J)=\mu_x(B^0)-\mu_x(B^\infty)
                         +\mu_y(B^\infty)-\mu_y(B^0).   \tag{24}
\]
This is exactly the classical perimeter telescoping already used in the
[planar primitive-cost bound](wall_cost_products.md#a-sharp-planar-cost-bound).

If \(J\subset J'\) have lengths \(m,m-1\), use the same two limiting
global lexicographic orders on both ideals. The actual one-dimensional
quotient rule deletes one monomial from each parent standard set; write
the two deleted exponents as \(u,v\). They are staircase corners, with
\(|u|,|v|\le m-1\). Subtracting (24) gives the exact identity
\[
 \mathcal L(J)-\mathcal L(J')=u_x-u_y-v_x+v_y,
 \qquad |\mathcal L(J)-\mathcal L(J')|\le2(m-1).          \tag{25}
\]
The bound is sharp over every field: \(J_m=(x^m,y-x)\) has
\(\mathcal L(J_m)=m(m-1)\), so its increment over \(J_{m-1}\)
is \(2(m-1)\). It has two initial sets when \(m\ge2\).
No radicality or genericity is used in (25).

In the three-pole family (1), each profile edge exchanges one axis tip
\(x_i^A\) for \(x_j^B\). Keeping its integer multiple gives
\[
 \mathcal L(I_n)=3\sum_{A,B\ge1,\ A+B\le n}(A+B)
               =n^3-n,\qquad
 \mathcal L(I_n)-\mathcal L(I_{n-1})=3n(n-1).            \tag{26}
\]
For (14), the product identity also holds with actual lengths:
\[
 \mathcal L(H_{m,r})=b_r\mathcal L(T_m)+m\mathcal L(J_r).
                                                               \tag{27}
\]
Both branches used above have cyclic coordinates, so their limiting
standard sets are the two pure-power chains and their actual lengths
are respectively \(m(m-1)\) and \(r(r-1)\). Hence their product's
increment is \(2b_r(m-1)+r(r-1)\), always positive. For example the
length-36 and length-35 specialization with primitive increment \(-332\)
has actual lengths \(3168,2772\), with increment \(396\).
This cancellation is real, but specific to these changing hook branches.

### Growing cone multiplicity defeats the spatial step estimate

Here are all-field families, with no number-theoretic specialization,
whose actual lengths have jumps of **both** signs of order \(n^{7/3}\).
Only one dimension is removed, and their fans have only \(O(n^{1/3})\)
cones. Thus even \(\mathcal L\) has no uniform positive or absolute
\(O(n^2\log^C(n+1))\) increment for fixed \(C\).

For an integer \(s\ge1\), first use the convex antichain
\[
 S_s^+=\{(s-i,i^2):0\le i\le s\},\qquad
 D_s^+=\operatorname{Div}(S_s^+),\qquad
 |D_s^+|=\frac{s(s+1)(2s+1)}6+s+1.
\]
In \(K[x,y]\), let \(A_s^+\) be the ideal complementary to this
downset, together with relations identifying every tip with the first
tip. Let \(B_s^+\) additionally kill that common tip. Then
\[
 m_s^+=\operatorname{length}(A_s^+)
       =\frac{s(s+1)(2s+1)}6+1,
 \qquad \operatorname{length}(B_s^+)=m_s^+-1.
\]
All tips are maximal, so their span is socle and these are actual ideals
over every field. The parent tips form one parallel residue class.
Every parent standard set is \((D_s^+\setminus S_s^+)\cup\{t_i\}\),
where \(t_i=(s-i,i^2)\). It is coherent precisely when that tip is
cheapest. All \(s+1\) choices occur: the consecutive thresholds are
\(w_x/w_y=2i+1\), \(0\le i<s\). The child has one standard set.
The consecutive profile differences are primitive, and
\[
 g(A_s^+)=s+1,\quad g(B_s^+)=1,\quad
 \mathcal L(A_s^+)=W(A_s^+)=s(s+1),\quad
                         \mathcal L(B_s^+)=0.            \tag{28}
\]

For the opposite sign use the concave antichain
\[
 S_s^-=\{(i,s^2-i^2):0\le i\le s\},\quad
 D_s^-=\operatorname{Div}(S_s^-),\quad
 m_s^-=(s+1)(s^2+1)-\frac{s(s+1)(2s+1)}6.
\]
Let \(A_s^-\) be the complementary monomial ideal, and \(B_s^-\)
its quotient by the single relation \(\sum_{t\in S_s^-}x^{t_x}y^{t_y}\).
Their lengths are \(m_s^-,m_s^--1\). The parent has one standard set;
the child omits the most expensive tip. The same \(s+1\) thresholds
give \(g(B_s^-)=s+1\) and \(\mathcal L(B_s^-)=s(s+1)\).
The relation is nonzero in every characteristic, and its span is socle.

For either sign put \(m=m_s^\pm\), \(f=x^m\), and take the second
branch \(Q_m=(x^m,z-x)\). It has two planar initials and
\(\mathcal L(Q_m)=m(m-1)\), while \(W(Q_m)=2\). Form
\[
 P_s^\pm=(f,yz,yA_s^\pm,zQ_m),\qquad
 C_s^\pm=(f,yz,yB_s^\pm,zQ_m).                           \tag{29}
\]
The common annihilator belongs to both first branches, since their
largest possible \(x\)-exponent is \(s<m\). The product proof gives
actual lengths \(n=3m,n-1\) and an actual one-dimensional kernel.
The two full fan counts are \(2(s+1),2\) for the positive sign and
\(2,2(s+1)\) for the negative sign. From the same product identity,
\[
 \boxed{\quad
 \mathcal L(P_s^\pm)-\mathcal L(C_s^\pm)
              =\pm\{s\,m(m-1)+2s(s+1)\}
              =\pm\Theta(n^{7/3}).\quad}                \tag{30}
\]
The upper and lower orders follow from \(m_s^\pm=\Theta(s^3)\).
The large contribution is the second branch's long profile edge copied
across the \(s+1\) first-branch cones. Primitive normalization of that
edge is constant. In fact the variable fan has
\(W=2(s+1)^2\) and the fixed fan has \(W=2\), both small here.
This mechanism is therefore independent of the gcd oscillation in (17).

Every wall in the two-factor variable fan touches an interior grid
intersection, so its \(\ell_+\) equals \(\mathcal L\). The fixed
two-cone fan has no positive two-face and \(\ell_+=0\). Consequently
the same order and both signs occur for \(\ell_+\), though its exact
increments differ from (30). Neither a cubic-polylogarithmic global
\(\mathcal L\) or \(\ell_+\) bound, nor a cubic total-variation bound
for these actual lengths, is disproved by (30).

### Actual lengths can exceed primitive cost quadratically

There is a separate cost to keeping the multiples. For every length-\(n\)
quotient, the staircase moment cap \(\sum_i\mu_i\le\binom n2\)
gives \(\|\mu(B)-\mu(C)\|_1\le n(n-1)\) at each edge. Since a
primitive mixed direction has one-norm at least two,
\[
 W(I)\le\mathcal L(I)\le\frac{n(n-1)}2W(I).              \tag{31}
\]
The quadratic order of this loss is attained even with a fixed fan.
For the existing characteristic-zero five-cone family
\(I_k=(x^2,y^k,z-x-y)\), \(k\ge2\), \(n=2k\), the five
actual edge lengths are
\[
 2,\quad2(k-1),\quad2k,\quad2(k-1)^2,\quad2k(k-1).
\]
They follow by summing the five explicit standard sets in
[the classified fan](adjacent_rank_fans.md#a-sharp-linear-loss-from-retaining-local-core-length).
Thus \(\mathcal L(I_k)=n^2-n+2\), while \(W(I_k)=10\).
An actual-length budget is a stronger sufficient target, with a sharp
quadratic possible loss, rather than an automatically easier replacement
for primitive length. The main global count remains unresolved.

### A collective section problem, without an ideal-chain recurrence

There is a precise way to ask for a global actual-length budget. Let
\(P=\operatorname{conv}\{\mu(B):B\text{ a positive initial standard set}\}\).
For coordinate \(i\) and integer \(k\), let \(c_i(k)\) count actual
positive profile edges whose endpoints lie on opposite sides of
\(x_i=k+1/2\). Integer endpoints give the elementary identity
\[
              \mathcal L(I)=\sum_{i=1}^3\sum_k c_i(k).   \tag{32}
\]
Each crossing point is a vertex of the genuine section
\(P\cap\{x_i=k+1/2\}\): a plane transverse to a polytope edge exposes
that intersection as a section vertex. Distinct edges give distinct
crossings, since the plane contains no integer vertex. Thus \(c_i(k)\)
is at most the number of section vertices.

Every coordinate moment is between zero and \(\binom n2\), so a uniform
\(O(n\log^C(n+1))\) bound on these section vertices would imply the
desired cubic-polylogarithmic \(\mathcal L\) bound. Even an averaged
bound on the \(c_i(k)\) would suffice. This is an **unproved** possible
next task, not a theorem about arbitrary quotient sections.
It retains the final fan and makes no charge along ideal deletions.
These sections usually have rational vertices and are not planar quotient
state polytopes. The planar \(G_2\) theorem cannot simply be applied to
them. They also differ from the hull of profiles lying at one *integer*
coordinate level in the [coordinate-slice theorem](coordinate_slice_wall_budget.md).

### A sharp cubic global budget with two quadratic coordinates

One substantial restriction does permit a uniform actual-length bound.
If a length-\(n\) quotient contains monic scalar quadratics in \(x\)
and \(y\), then, over every field,
\[
       g(I)\le(n+1)^2,\qquad
                  \mathcal L(I)\le3(n+1)^3.              \tag{33}
\]
Both powers are sharp. This is a specialization of the existing
[one-box rank estimate](weighted_facet_incidence.md#a-rank-sensitive-estimate-inside-one-box),
with a direct positive-fan proof and explicit constants below; it is
not claimed as a new general polytope estimate.

Every standard label has first two exponents in \(\{0,1\}\) and last
exponent at most \(n-1\). Thus the first two profile coordinates are
integers in \([0,n]\). Every positively exposed vertex or edge point
is the unique *lowest* point of \(P\) over its \((x,y)\)-projection:
a lower point in the same fiber would have smaller cost because \(w_z>0\).
Projection is therefore injective on the union of the positive vertices
and edge interiors. There are at most \((n+1)^2\) positive vertices.

Write a positive edge displacement as \(k_eq_e\), with primitive
\(q_e\). The projected matroid-edge principle gives
\((q_e)_x,(q_e)_y\in\{-1,0,1\}\), at least one nonzero, and
\(|(q_e)_z|\le n-1\). Indeed a nonzero ground-label difference parallel
to the edge has these bounds, and a pure \(z\) direction cannot be
orthogonal to a positive cost. It follows that \(\|q_e\|_1\le n+1\).
The \(k_e-1\) interior lattice points on each edge project to distinct
integer pairs, distinct also from projected vertices. Hence, with \(V\)
positive vertices and \(E\) actual positive edges,
\[
  V+\sum_e(k_e-1)\le(n+1)^2,\qquad E\le3V,
  \qquad\sum_e k_e\le3(n+1)^2.
\]
The graph inequality is ordinary planarity, with segments and polygons
handled separately. Multiplying by the primitive norm bound proves (33).

The [four-line parity-point family](parity_point_fan_updates.md) makes
the cubic power sharp and gives a stronger characteristic-zero update
example than (30). In that note's notation, every profile edge of
\(J_s=I(S_s)\) and \(I_s=I(P_s)\) has lattice multiple one. For
\(J_s\) an edge flips one opposite pair, with displacement
\((1,\pm1,k)\). For \(I_s\), the support function adds half the
minimum absolute pair cost. At a new wall the profile displacement is
half a signed sum of the two competing pair differences. Its first two
coordinates are in \(\{-1,0,1\}\) with at least one unit; if both
vanished the displacement would be a pure \(z\) direction and would
not give a positive wall. Thus it too is primitive. Consequently
\[
 \mathcal L(J_s)=W(J_s),\quad \mathcal L(I_s)=W(I_s),\quad
 \mathcal L(I_s)-\mathcal L(J_s)
                  =\frac{2s^3+3s^2-8s-6}{3}.           \tag{34}
\]
Together with (33), cubic order is the sharp worst-case positive or
absolute one-point actual-length update order in the two-quadratic class.

Both signs occur on those same four fixed lines. Put
\(G_s=\{0,1\}^2\times\{0,\ldots,2s-1\}\),
\(E_s=G_s\setminus S_s\) and \(F_s=G_s\setminus P_s\).
Their radical ideals have lengths \(4s,4s-1\), with one point deleted.
The established grid-complement correspondence sends a standard set
\(B\) to \(G_s\setminus(v-B)\), \(v=(1,1,2s-1)\). Its profile is
\(\mu(G_s)-|B|v+\mu(B)\), a fixed translation, so it preserves actual
edge lengths as well as the fan. Therefore the length-\(4s\) to
length-\(4s-1\) increment is the negative of (34), and its cone-count
increment is \(-(s^2-s-1)\). This is an application of existing duality,
not a new complementary-ideal theorem; see
[Dimitrova--He--Robbiano--Stigler, Theorem 4.8 and Corollary 4.11](https://arxiv.org/pdf/1704.03680v2)
and the [existing box-colon proof](interleaving_research.md#an-exact-colon-identity-on-every-box-initial-cone).

The all-field construction (30), the sharp restricted global estimate
(33) and the characteristic-zero point updates (34) have different
scopes. None controls unrestricted three-variable quotients or supplies
a supergeneric count. In particular (33) cannot be used after an
arbitrary rational section or vector-space contraction.

## Reproduction and prior-art boundary

~~~bash
python3 -B scripts/run_ideal_step_incidence_checks.py
./build.sh
python3 tests/test_smoke.py
~~~

The report is `results/research/ideal_step_incidence_checks.json`.
The checker compares complete actual point-ideal fans from a separate
border traversal with (4), verifies separator eigenlines, whole-cone
S-pairs, ray intersections and the exact wall recurrence. Exhaustive
finite residue averages check the local factors without floating
arithmetic. Larger exact integer sums illustrate both jump signs.
Finite checks do not establish (7) or replace its proof.
The product extension independently compares all branch-product cones
with complete spatial point-ideal border fans, checks ordinary S-pairs
and the separator eigenline, and reconstructs whole wall inventories.
It includes \(W(H_{11,12})=2156\), \(W(H_{12,12})=1824\) over
\(\mathbb F_{101}\), a verified specialization of the fixed-curve
example. Exact gcd, multiplicative, parity and Möbius identities are
checked separately through larger parameter ranges. Formula (22) is an
analytic theorem, not a fit to those data.
The actual-length extension verifies (26)--(27) on those complete fans,
then constructs both socle-product signs over characteristics zero, two
and three. It independently checks commuting multiplication, every
monomial successor on a finite spanning ground, the lifted ideal
generators, the one-dimensional socle kernel, all product cones and
ordinary S-pairs. For example the positive length-45 to length-44 step
has actual lengths \(864,210\), while the negative length-78 to
length-77 step has lengths \(650,2624\), in each tested field.
The closed formulas, rather than these finite values, prove (30).

The curve fan, eigenline interpretation and incidence definitions are
preserved established inputs from the cited notes. Arithmetic factors
depending on the prime divisors of a composition total are classical:
[Bubboloni--Luca--Spiga, Theorem 4](https://arxiv.org/pdf/1202.1670),
[Thomas, Theorem 1.1](https://arxiv.org/pdf/2001.12001), and
[Tóth, Theorem 2.2](https://arxiv.org/pdf/2006.04500) give
unweighted coprimality-composition asymptotics. Their displayed theorem
statements and the residue/divisibility setup were inspected; their
full error-term proofs were not independently audited. The reciprocal
square-root gcd weight is derived explicitly above; no claim that it
escapes the classical multiplicative-function framework is made.
Exact priority of the weighted formulation and its fan consequence
remains unassessed. The sharp unrestricted \(G_3\) question, an actual
polynomial lower separation and a global cubic budget remain unresolved.

For the extension, the arithmetic function itself is established:
[Tóth, equations (13) and (51)](https://cs.uwaterloo.ca/journals/JIS/VOL13/Toth/toth10.pdf)
identify \(F(m)/m\) as the mean element order of a cyclic group and give
the equivalent prime-power formula. Those displayed formulas and their
setup were inspected; no new number-theoretic function is claimed.
The uniform \(C_2(r)=\Theta(r\log r)\) input is the classical
Corteel--Rémond--Schaeffer--Thomas result restated in
[Elizalde--Galván, Theorem 2.3](https://arxiv.org/pdf/2312.16353v1),
whose statement was inspected. The standing article/preprint access
assumption is recorded in [the research note](interleaving_research.md).
The all-order product formula predates this continuation in this repository;
[Conca--Römer, Theorem 1.2](https://arxiv.org/pdf/0707.0678) concerns
generic reverse-lexicographic initials of standard graded fibre products,
rather than the fixed-coordinate, all-order wall variation here. Its
statement and introductory definitions were inspected. These comparisons
do not establish priority of the elementary product consequence (16)--(22),
or justify a standalone paper. They do justify abandoning a uniform
quadratic-increment or arbitrary-chain-positive-variation proof of the
global budget. No runtime lower bound for fan-update algorithms is asserted.

For (23)--(32), the perimeter telescoping, standard-set/profile dictionary,
socle relation-space construction and all-order fan product are established
inputs already cited in this repository. The planar deletion identity is
their direct consequence, and the spatial families expose the separate
fan-multiplicity term in (27). The section identity is elementary segment
counting; it supplies an open sufficient condition, not a section bound.
Exact priority of these combined observations has not been established.
No claim of a new general polytope theorem, algorithmic lower bound,
unrestricted exponent improvement or standalone paper follows.
For (33)--(34), the box rank bound and the parity-family classification
are established inputs explicitly linked above. The radical grid
duality in Theorem 4.8 and Corollary 4.11 of the cited primary article
was re-read at its statements and radical-grid proof. Their role is
fan-preserving duality, not a new bound on unrestricted complexity.
