# Cubic incidence jumps under deletion of one point

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
Equality holds when \(n\) is prime. Thus the failure in (2) does not
refute a quadratic increment conjecture for \(W\), or a cubic global
budget for either quantity. No such unrestricted budget is proved here.
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
