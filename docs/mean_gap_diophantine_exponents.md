# Mean relative-gap exponents and Diophantine approximation

## Status and the smooth-number question

This note proves that average gap exponents in a positive exponent lattice
have limiting extremes equal to its uniform and ordinary linear-form
exponents. For every pair of distinct primes the averaged normalized gap
exponent has **liminf exactly 1**. Its limit exists if and only if the ratio
of prime logarithms has irrationality exponent 2; any limit is therefore 1.
This does not evaluate the irrationality exponent of that ratio or prove
convergence for a particular prime pair.

[Languasco--Luca--Moree--Togbe, 2025, Section 4.3](https://doi.org/10.1007/s12188-025-00293-9)
asks whether the mean of

\[
 \rho_j=\frac{-\log(n_{j+1}/n_j-1)}{\log\log n_j}
\]

converges for the increasing two-prime smooth sequence. The definition,
Lemma 4.1, Question 4.1 and Remark 4.1 were inspected in the published
version. The theorem below identifies both limiting extremes of this mean.
The final comparison separates this result from classical inputs and
records the scope of the literature assessment. The proof has been
internally audited; it has not undergone independent specialist review.

## Definitions and main statement

Let \(k\geq2\), and let \(w=(w_1,\ldots,w_k)\) have positive coordinates
linearly independent over \(\mathbb Q\). Here
\(\mathbb N=\{0,1,2,\ldots\}\). Arrange

\[
 \{w\cdot e:e\in\mathbb N^k\}=\{0=t_1<t_2<\cdots\}.
\]

For \(t_j\geq2\), put

\[
 \rho_j=\frac{-\log(\exp(t_{j+1}-t_j)-1)}{\log t_j},
 \qquad
 A_w(T)=\frac{1}{N_*(T)}\sum_{2\leq t_j\leq T}\rho_j,
\]

where \(N_*(T)=\#\{j:2\leq t_j\leq T\}\). All logarithms are natural.
Changing the fixed lower cutoff discards finitely many terms only.

For \(H\geq1\), define

\[
 \delta_w(H)=\min_{0<\|z\|_\infty\leq H, z\in\mathbb Z^k}|w\cdot z|,
 \qquad L_w(H)=-\log\delta_w(H),
\]

and define the exponents, allowing value \(+\infty\), by

\[
 \widehat\omega(w)=\liminf_{H\to\infty}\frac{L_w(H)}{\log H},
 \qquad
 \omega(w)=\limsup_{H\to\infty}\frac{L_w(H)}{\log H}.
\]

**Theorem.** Under these hypotheses,

\[
 \boxed{\liminf_{T\to\infty}A_w(T)=\widehat\omega(w),
 \qquad\limsup_{T\to\infty}A_w(T)=\omega(w).}              \tag{1}
\]

Thus averaging over the whole ordered lattice prefix preserves both
Diophantine exponents. The averages converge precisely when the two
exponents agree, with the same convention for an infinite limit.

## Counting and reciprocal logarithms

Let \(C_w(T)=\#\{e\geq0:w\cdot e\leq T\}\). Comparing unit cubes with
the positive simplex gives

\[
 C_w(T)=\frac{T^k}{k!\prod_iw_i}+O_w((T+1)^{k-1}).          \tag{2}
\]

The same estimate gives

\[
 \frac{1}{N_*(T)}\sum_{2\leq t_j\leq T}\frac1{\log t_j}
                  =\frac{1+o_w(1)}{\log T}.               \tag{3}
\]

Indeed, the terms below \(T/(\log T)^2\) form an
\(O_w((\log T)^{-2k})\) fraction, and each reciprocal logarithm is at
most \(1/\log2\). On the remaining terms the denominators lie between
\(\log T-2\log\log T\) and \(\log T\). The lower-tail contribution is
\(o(1/\log T)\), proving (3).

Pigeonholing the \((\lfloor H\rfloor+1)^k\) distinct numbers
\(w\cdot e\), \(0\leq e_i\leq\lfloor H\rfloor\), inside an interval of
length \(O_w(H)\) proves

\[
 \delta_w(H)\ll_w H^{1-k}.                                \tag{4}
\]

In particular \(\delta_w(H)\to0\), and
\(k-1\leq\widehat\omega(w)\leq\omega(w)\).

## Upper bound on the average

Write \(w_{\min}=\min_i w_i\). Multiplication by the corresponding
generator gives \(t_{j+1}\leq t_j+w_{\min}\). If \(t_j\leq T\), both
exponent vectors therefore have coordinates at most
\((T+w_{\min})/w_{\min}\). For a fixed constant \(c_w>0\) and all
sufficiently large \(T\), their difference has sup norm at most \(c_wT\).
Consequently

\[
 t_{j+1}-t_j\geq\delta_w(c_wT),\qquad
 \exp(t_{j+1}-t_j)-1\geq\delta_w(c_wT).
\]

For every term included in the average,
\(\rho_j\leq L_w(c_wT)/\log t_j\). Equation (3) gives

\[
 A_w(T)\leq(1+o_w(1))\frac{L_w(c_wT)}{\log T}.             \tag{5}
\]

Fixed rescaling of \(H\) leaves the liminf and limsup defining the
exponents unchanged. Taking a subsequence attaining the finite liminf of
\(L_w(H)/\log H\) proves \(\liminf A_w\leq\widehat\omega(w)\).
If \(\omega(w)<\infty\), the eventual upper bound in its definition
proves \(\limsup A_w\leq\omega(w)\). Infinite upper bounds need no
argument. The multiplicative error in (5) is harmless on each selected
subsequence with a finite limiting exponent.

## Translation turns a small form into a density-one gap bound

Choose \(z\) attaining \(\delta_w(H)\), and replace it by its negative
if needed so that \(w\cdot z=\delta_w(H)>0\). Write
\(z=z^+-z^-\), and set \(a_z=w\cdot z^-\leq H\sum_iw_i\).

For every \(e\geq z^-\), the vector \(e+z\) is nonnegative. It supplies
a later lattice point at distance \(\delta_w(H)\), whether or not that
point is the actual successor. Thus

\[
 0<t_{j+1}-t_j\leq\delta_w(H)
       \quad\text{whenever }e_j\geq z^-.                   \tag{6}
\]

The number of such points up to height \(T\) is **exactly**
\(C_w(T-a_z)\), by translation. If \(H/T\to0\), equation (2) shows
that their fraction among the points in the average is

\[
 f(H,T)=1-O_w(H/T+1/T).                                   \tag{7}
\]

For large \(H\), (4) and \(\exp x-1\leq2x\) for small positive \(x\)
give, on these points,

\[
 \rho_j\geq\frac{L_w(H)-\log2}{\log t_j}
          \geq\frac{L_w(H)-\log2}{\log T}.                 \tag{8}
\]

All other included terms have a fixed lower bound \(-K_w\): their
relative gaps are at most \(\exp(w_{\min})-1\), and
\(\log t_j\geq\log2\). Combining (7)--(8) yields

\[
 A_w(T)\geq
 f(H,T)\frac{L_w(H)-\log2}{\log T}-K_w(1-f(H,T)).           \tag{9}
\]

Set \(H=T/(\log T)^2\). The eventual lower bound defining
\(\widehat\omega(w)\), together with \(\log H/\log T\to1\), proves
\(\liminf A_w\geq\widehat\omega(w)\), including the infinite case.

For the limsup, choose heights \(H_\nu\) along which
\(L_w(H_\nu)/\log H_\nu\) tends to \(\omega(w)\), and take
\(T_\nu=H_\nu(\log H_\nu)^2\). Equations (7)--(9) prove
\(\limsup A_w\geq\omega(w)\). This remains valid when that limsup is
infinite. Together with (5), these prove (1). \(\square\)

## Finite comparison bounds

The two estimates behind the limiting theorem have finite versions. Let
\(I(T)=\{j:2\leq t_j\leq T\}\) and
\(Q(T)=|I(T)|^{-1}\sum_{j\in I(T)}1/\log t_j\). The upper estimate is
exactly

\[
 A_w(T)\leq L_w(H_T)Q(T),\qquad
 H_T=\left\lfloor\frac{T+w_{\min}}{w_{\min}}\right\rfloor. \tag{12}
\]

For any positive form \(\varepsilon=w\cdot z<\log2\), let \(f\) be
the actual fraction of \(e_j\), \(j\in I(T)\), satisfying \(e_j\geq z^-\).
Then

\[
 A_w(T)\geq
 f\frac{-\log(\exp\varepsilon-1)}{\log T}-K_w(1-f),        \tag{13}
\]

where \(K_w=\max(0,\log(\exp w_{\min}-1))/\log2\).
The fraction is computed by translation, with the fixed lower-cutoff
points removed. These bounds require neither a rate of Diophantine
approximation nor an assumption about future best approximants.

For an explicit counting bound, put \(W=\sum_iw_i\) and
\(V=1/(k!\prod_iw_i)\). The same unit-cube argument gives

\[
 VT^k\leq C_w(T)\leq V(T+W)^k\quad(T\geq0).
\]

Consequently the number of points excluded by a move of coefficient
height at most \(H\) is uniformly an \(O_w(H/T+1/T)\) fraction.
This uniformity allows the minimizing form to change arbitrarily with
\(H\); no single move is assumed to persist as \(T\) grows.

## Two generators and the published mean-gap question

For \(k=2\), put \(\theta=w_1/w_2\). Let

\[
 \mu_{\rm irr}(\theta)=\sup\{\tau:
       |\theta-p/q|<q^{-\tau}\text{ for infinitely many }p,q\}
\]

be its ordinary irrationality exponent. The classical one-dimensional
approximation uses \(p\in\mathbb Z\), \(q\in\mathbb Z_{\geq1}\).
The classical one-dimensional
approximation identities are

\[
 \widehat\omega(w)=1,\qquad
 \omega(w)=\mu_{\rm irr}(\theta)-1.                        \tag{10}
\]

For clarity, the sup-norm bound on both coefficients in \(\delta_w\)
does not change these exponents. If
\(\psi(Q)=\min_{1\leq q\leq Q}\|q\theta\|\), then for a fixed
\(c_\theta>1\) and all large \(H\),

\[
 w_2\psi(H)\leq\delta_w(H)
                 \leq w_2\psi(H/c_\theta).
\]

The ordinary exponent identity follows directly by dividing
\(|q\theta-p|\) by \(q\). For the uniform exponent, Dirichlet
approximation gives \(\psi(Q)\ll1/Q\). If \(q_n,q_{n+1}\) are
consecutive continued-fraction denominators, best approximation gives

\[
 \psi(q_{n+1}-1)=\|q_n\theta\|>
                     \frac1{q_n+q_{n+1}}>\frac1{2q_{n+1}}.
\]

These prove the uniform identity in (10). Thus the mean has

\[
 \boxed{\liminf A_w=1,\qquad
        \limsup A_w=\mu_{\rm irr}(\theta)-1.}              \tag{11}
\]

For distinct primes \(p,q\), take \(w=(\log p,\log q)\).
The heights are the logarithms of the smooth integers. Passing from a
height cutoff to the first \(N\) terms changes no limiting extreme,
because the heights are distinct. The finitely many discarded terms also
do not affect the conclusion.

The published paper conditions its mean on \(\rho_j\geq0\). There are
only finitely many excluded terms: one-sided approximation supplies,
on each pure-prime axis, a fixed improving ratio less than 2; its
multiples cover the complement of a finite exponent box. Hence the
conditioned mean has the same limiting extremes. Equation (11) gives

\[
 \lim_N\mu(p,q;N)\text{ exists}
 \quad\Longleftrightarrow\quad
 \mu_{\rm irr}(\log p/\log q)=2,
\]

and any such limit equals 1. This is a reduction of the convergence
question; it supplies no estimate proving the right-hand assertion for a
particular prime pair. A limit strictly larger than 1 is incompatible
with the theorem, even though finite numerical means can exceed 1.

## A separate proof of the sharp two-generator liminf

The unconditional conclusion \(\liminf A_w=1\) can also be verified
without the general transfer theorem. Relabel the weights so that
\(0<w_1<w_2\), and put \(\theta=w_1/w_2\).

First, the successor ratios tend to 1. For any fixed \(\eta>0\),
one-sided approximation gives a move using a sufficiently large power on
either pure axis with ratio less than \(1+\eta\). All multiples of those
axis powers admit the move, and their complement is a finite box. In
particular the relative gaps \(r_j=\exp(t_{j+1}-t_j)-1\) are eventually
less than 1.

After discarding those finitely many exceptions, their logarithms have the
sign needed for \(\rho_j\geq-\log r_j/\log T\). Also
\(r_j\leq c_w(t_{j+1}-t_j)\), since the log gaps are at most \(w_1\).
Their sum up to height \(T\) is at most \(c_w(T+w_1)\), by telescoping.
Jensen's inequality and (2) therefore give

\[
 A_w(T)\geq
 \frac{\log N_*(T)-\log(c_w(T+w_1))}{\log T}+o(1)
                         =1+o(1).
\]

For the reverse bound, let \(q_n,q_{n+1}\) be consecutive convergent
denominators of \(\theta\), and take
\(T_n=w_1(q_{n+1}-2)/2\). Both exponent vectors in any successor pair
starting below this height have first coordinate at most
\((T_n+w_1)/w_1=q_{n+1}/2\). Their first-coordinate difference has
absolute value less than \(q_{n+1}\). Best approximation gives

\[
 t_{j+1}-t_j\geq w_2\|q_n\theta\|
                    >\frac{w_2}{2q_{n+1}}\gg_w\frac1{T_n}.
\]

If that coordinate difference is zero, the gap is at least \(w_2\), so
the same bound holds. Hence \(-\log r_j\leq\log T_n+O_w(1)\)
for every included pair. Equation (3) now gives
\(A_w(T_n)\leq1+o(1)\). This proves the sharp liminf independently
of the translation-based lower estimate.

## Concentration, typical weights and failure of universal convergence

**Corollary.** If
\(\widehat\omega(w)=\omega(w)=\kappa<\infty\), then for every
\(\eta>0\),

\[
 \frac{\#\{j:2\leq t_j\leq T,\ |\rho_j-\kappa|>\eta\}}
      {N_*(T)}\longrightarrow0.                           \tag{14}
\]

In fact the mean of \(|\rho_j-\kappa|^r\) tends to zero for every fixed
\(r>0\).

**Proof.** The ordinary exponent bounds all sufficiently late individual
\(\rho_j\) above by \(\kappa+\eta/2\), using the coefficient-height
bound on each successor pair. With \(H=T/(\log T)^2\), the uniform
exponent and (6)--(8) bound \(\rho_j\) below by
\(\kappa-\eta/2\) on a fraction tending to 1. This proves (14).
The individual exponents are eventually bounded above, and are bounded
below by \(-K_w\), so the same argument gives convergence of every
fixed positive absolute moment. \(\square\)

For Lebesgue almost every positive weight vector, the common exponent is
\(k-1\). Here is the classical metric argument, included to make this
consequence self-contained. On a compact positive box, the slab
\(|w\cdot z|<h^{-(k-1+\eta)}\), \(\|z\|_\infty=h\), has measure
\(O(h^{-k-\eta})\), by integrating in a coordinate with coefficient
\(h\). There are \(O(h^{k-1})\) such integer vectors. The sum of these
measures converges, and Borel--Cantelli gives \(\omega(w)\leq k-1\)
almost everywhere. Equation (4) gives the reverse and uniform bounds.
Countably many compact boxes cover the positive quadrant. Thus the mean
and the concentrating individual exponent are \(k-1\) almost everywhere.
This full-measure statement gives no conclusion for the particular
countable collection of prime-logarithm vectors.

Convergence is nevertheless false for general real generators. For any
\(1<\gamma<\infty\), choose an irrational continued fraction whose
successive coefficients satisfy
\(a_{n+1}=\max(1,\lfloor q_n^{\gamma-1}\rfloor)\) after an initial
positive prefix. Then \(q_{n+1}\sim q_n^\gamma\), and the classical
best-approximation formulas give

\[
 \omega(1,\theta)=\limsup_n\frac{\log q_{n+1}}{\log q_n}
                      =\gamma,\qquad
 \widehat\omega(1,\theta)=1.
\]

The mean gap exponent for the semigroup generated by \(e\) and
\(e^\theta\) therefore has liminf 1 and limsup \(\gamma\).
Bounded partial quotients give the convergent case \(\gamma=1\).
Choosing \(a_{n+1}=q_n^n\) instead gives limsup \(+\infty\).
These are examples for real generators, not counterexamples for primes.

## Reproducible finite checks

Run from the repository root:

```bash
python3 -B scripts/run_mean_gap_checks.py
./build.sh
python3 tests/test_smoke.py
```

The checker uses genuinely irrational models in dimensions 2, 3 and 4:
square-root weights, prime-logarithm ratios, and explicitly specified
infinite continued fractions with bounded, polynomially growing or
Liouville coefficients. Signs and floors are separated by rational
enclosures. Logarithms of integers use the positive atanh series with an
explicit tail bound; continued fractions use cylinder endpoints.

Minimum forms are found by nearest-integer elimination and independently
checked by exhaustive signed boxes at the smaller heights. The checker
tests all successor lower bounds, the exact translated-simplex counts,
and every feasible translated gap bound. Prime-log orderings are also
compared against exact integer products. A rationally dependent model
must be rejected. The source-bound report is
`results/research/mean_gap_checks.json`.

The limiting identities are consequences of the analytic proof, not of
these finite checks. Decimal evaluations of the means and bounds are
diagnostic; their transcendental values are not interval-certified.

## Primary comparison and contribution boundary

The retained contribution is the passage from minimum linear forms to
**averaged** relative-gap exponents, uniformly in the dimension, and the
sharp two-prime mean-liminf and convergence equivalence. The following
are prior art or comparison results, not new claims here.

| Primary source and inspected scope | Established result | Boundary of the comparison |
|---|---|---|
| [Languasco--Luca--Moree--Togbe, 2025](https://doi.org/10.1007/s12188-025-00293-9), Sections 4.2--4.3, including the complete mean discussion | Defines this same normalized mean; Lemma 4.1 gives eventual constant bounds; Question 4.1 asks whether the limit exists; Remark 4.1 identifies the conditioned/unconditioned comparison | The sharp liminf and irrationality-exponent equivalence strengthen the statement for that same statistic. The full convergence question for prime pairs is reduced, not decided. The finite table uses neighbor windows; the theorem here concerns full prefixes. |
| [Berczes--Dujella--Hajdu, 2014](https://math.unideb.hu/sites/default/files/inline-files/sunit_submit.pdf), introduction and Theorems 2.1--2.3 | Log-hyperplane geometry, effective neighbors, continued-fraction structure and occurrence information for two primes | The lattice interpretation and neighbor theory are established inputs. They are not presented as discoveries of this theorem. |
| [Fan--Queffelec--Queffelec, 2024](https://ems.press/content/serial-article-files/47663?nt=1), Section 2.1 statements and the proof of Theorem 2.4 and its subsequent remark, pp. 70--75 | Irrationality-dependent counting and pointwise relative-gap bounds for the Furstenberg set; infinitely many convergent gaps of inverse-logarithmic size | These are pointwise/counting statements; they do not identify this prefix mean's two limiting extremes. |
| [Weiss, 2023](https://doi.org/10.1007/s00010-023-00958-x), definition of the gap-ratio function and the Kronecker-sequence discussion around Theorem 2.2 | Gap quantiles on the circle, three-gap lengths and multiplicities, and subsequence bounds | These statistics and the finite circle model differ from the present positive-cone prefix mean. Three-gap and best-approximation facts remain classical. |
| [Chudnovsky, arXiv:2610.28079v1](https://arxiv.org/html/2610.28079v1#S11), Theorem 11.1 and Corollary 11.5 statements | A recent preprint claims an upper bound of 3 for the irrationality exponent of rational-logarithm ratios | That claim is not an exponent-two result or a mean-gap theorem. Its full proof and reported formalization have not been independently checked here, and no bound from it is imported. |

The positive comparison with the published mean-gap question supports
the contribution assessment; unsuccessful searches are not evidence of
priority. No exhaustive all-literature priority guarantee or independent
peer-review claim is made. Classical counting, Borel--Cantelli,
continued-fraction identities and the interpretation of logarithms as
exponent weights are explicitly separated from the transfer theorem.

## Internal proof audit

The upper estimate uses the smallest linear form available at a fixed
height. The lower estimate propagates one small form to almost every
interior exponent vector by translation. The audit checks the following
points explicitly.

* Positivity bounds the coefficient height of **both** members of every
  successor pair, including successors beyond the cutoff.
* The varying denominator is handled by the reciprocal-log estimate (3),
  not by replacing every \(\log t_j\) with \(\log T\) without a tail bound.
* The translated feasible count is exact before removing the fixed lower
  cutoff; the latter removes at most a fixed number of points, uniformly
  in the chosen move.
* The infimum and supremum cases use different rescalings of height. The
  multiplicative errors are applied to finite-exponent subsequences or
  to positive quantities tending to infinity; no bounded-exponent
  assumption is hidden in the infinite cases.
* The two-generator specialization respects the sup-norm bound on both
  coefficients, the paper's nonnegative-term conditioning and all finite
  cutoff changes. A separate Jensen/continued-fraction argument recovers
  its unconditional sharp liminf.
* Rational dependence is excluded, prime-logarithm exponents are not
  evaluated, and full-measure or real-generator corollaries are not
  transferred to a fixed prime pair.
