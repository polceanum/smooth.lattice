# Finite certificates for relative gaps between smooth numbers

## Question and status

For a finite prime set \(P\), let \(S_P\) be its positive smooth integers,
including 1, and let \(s_P(m)\) be the smallest member of \(S_P\) strictly
larger than \(m\in S_P\). Given a rational \(\alpha>1\), consider the **entire**
exceptional set

\[
 D_\alpha=\{e\in\mathbb N^k:s_P(P^e)>\alpha P^e\}.
\]

The question is to describe this finite set, count its elements, and find its
largest smooth integer without generating the whole ordered prefix.
The following reduction gives an integer-only algorithm for three primes and
a small certificate for a global numerical example. Exact priority of the
combined reduction and algorithm is **unestablished**. This is a research
result for assessment, not a claim to have settled a recognized open problem.

The earlier random-access kernels are unchanged. This calculation concerns
all large relative gaps, rather than an individual prescribed rank.

## Divisibility and coordinate faces

Write \(R(e)=s_P(P^e)/P^e\). For the vector \(e^{(i)}\) obtained by setting
coordinate \(i\) to zero,

\[
 \boxed{R(e)=\min_{1\leq i\leq k}R(e^{(i)})}.                 \tag{1}
\]

Also \(R(e+f)\leq R(e)\) for every nonnegative \(f\).

**Proof.** Multiplying a successor candidate by \(P^f\) preserves its
ratio, proving monotonicity. Put \(m=P^e\), \(n=s_P(m)\), and divide both
by their gcd. Write the reduced pair as \(A<B\), so that \(A\mid m\),
\(\gcd(A,B)=1\), and \(B/A=R(e)\). Some prime \(p_i\) divides \(B\), and
it does not divide \(A\). Thus \(A\mid P^{e^{(i)}}\), and
\(P^{e^{(i)}}B/A\) is a smooth successor candidate. This proves
\(R(e^{(i)})\leq R(e)\). Monotonicity gives the reverse inequality for
every face. This also covers a reduced denominator equal to 1. \(\square\)

In particular \(D_\alpha\) is closed under taking divisors, and membership
is determined by its coordinate-face sections. Its complement is a monomial
ideal whose minimal generators each have at least one zero coordinate.
Indeed, any good exponent vector has a good face projection by (1).

For \(k\geq2\), every axis eventually becomes good. For distinct primes
\(p,q\), the irrationality of \(\log p/\log q\) supplies positive integers
\(a,b\) with \(1<q^b/p^a\leq\alpha\), by one-sided continued-fraction
approximation. Such a pair bounds the bad exponents on the \(p\)-axis.
Monotonicity then bounds \(D_\alpha\) inside a finite box. This is an
elementary finiteness argument, not a new quantitative Diophantine bound.

## An exact three-prime construction

Let \(P=\{p_0,p_1,p_2\}\) and assume
\(1<\alpha=u/v<\min P\). Find the first exponent \(U_i\geq1\) for which
there is a product of the other two primes satisfying

\[
 p_i^{U_i}<p_j^b p_\ell^c\leq (u/v)p_i^{U_i}.                \tag{2}
\]

This is exactly the first good exponent on that axis. To see this, cancel
common powers of \(p_i\) in any improving candidate. An earlier good axis
point would therefore give a witness of the form (2) at an exponent no
larger than that point. The restriction \(\alpha<\min P\) rules out an
improvement by a nontrivial integer multiplier.

Every bad point now lies in \(0\leq e_i<U_i\). Within this box a minimal
good point on the \((i,j)\)-face has both coordinates positive and has a
coprime improving numerator which is a pure power of the remaining prime.
Consequently it suffices to test

\[
 p_i^a p_j^b<p_\ell^c\leq (u/v)p_i^a p_j^b,
 \qquad 0\leq a<U_i,\quad 0\leq b<U_j.                    \tag{3}
\]

Only the smallest pure power strictly greater than the denominator needs
testing. Coordinatewise minima of these witnesses give all mixed minimal
generators; the three axis generators complete the list.

For each face retain its nonincreasing frontier

\[
 f_{ij}(a)=\min\bigl(\{b:\text{some mixed witness }(a',b')
                         \leq(a,b)\}\cup\{U_j\}\bigr).
\]

Then the full bad set is described exactly by

\[
\begin{split}
 &0\leq a<U_0,\quad 0\leq b<U_1,\quad 0\leq c<U_2,\\
 &b<f_{01}(a),\quad c<f_{02}(a),\quad c<f_{12}(b).          \tag{4}
\end{split}
\]

For each permitted \((a,b)\), its largest allowed third coordinate is
\(\min(f_{02}(a),f_{12}(b))-1\). Summing these column lengths counts
\(D_\alpha\). Comparing the corresponding integer products gives its
largest member. One exact successor query then determines the last
exceptional consecutive pair.

For fixed primes, put \(B=U_0+U_1+U_2\). Constructing the axis witnesses,
scanning the three faces, finding the maximum, and making the final
successor query take \(O_P(B^2)\) integer arithmetic operations. The integers
involved have \(O_P(B)\) bits. This is an arithmetic-operation bound, not a
unit-cost bit-complexity claim. The face frontiers use \(O(B)\) integers;
listing all maximal bad points adds its output storage.

For comparison, the ordered prefix through the last exceptional start has
\(\Theta_P(B^3)\) terms as \(B\) grows. Its logarithmic height is at most
\(\sum_i(U_i-1)\log p_i=O_P(B)\), and at least
\(\max_i(U_i-1)\log p_i=\Omega_P(B)\), since each preceding axis point is
bad. The usual simplex count gives the assertion. Thus the face algorithm
removes one exponent dimension relative to generating that full prefix.
This comparison does not establish superiority over the best available
specialized algorithm.

## A five-generator certificate

Take \(P=\{2,3,5\}\) and \(\alpha=1000001/1000000\). The entire set of
starts whose successor ratio exceeds \(\alpha\) is

\[
\begin{split}
 &0\leq a<573,\qquad 0\leq b<1412,\qquad 0\leq c<852,\\
 &\neg(a\geq174\ \text{and}\ b\geq1175),\\
 &\neg(b\geq1085\ \text{and}\ c\geq60).                    \tag{5}
\end{split}
\]

The five positive witnesses are the exact integer inequalities

\[
\begin{array}{rcl}
2^{573}&<&3^{237}5^{85}\leq\alpha\,2^{573},\\
3^{1412}&<&2^{399}5^{792}\leq\alpha\,3^{1412},\\
5^{852}&<&2^{1460}3^{327}\leq\alpha\,5^{852},\\
2^{174}3^{1175}&<&5^{877}\leq\alpha\,2^{174}3^{1175},\\
3^{1085}5^{60}&<&2^{1859}\leq\alpha\,3^{1085}5^{60}.
\end{array}                                                   \tag{6}
\]

Their multiples are all good starts. The independent negative checks
exclude every other primitive axis or face witness in the remaining box;
this is essential to the equality in (5).

Inclusion-exclusion gives the exact number of exceptional gaps:

\[
573\cdot1412\cdot852
-399\cdot237\cdot852
-573\cdot327\cdot792
+399\cdot237\cdot792
=\boxed{535261140}.
\]

The maximal exponent vectors of (5) are

\[
(173,1411,59),\quad(572,1084,851),\quad(572,1174,59).
\]

The largest exceptional consecutive pair is

\[
 m=2^{572}3^{1084}5^{851},\qquad
 n=2^{1459}3^{1648}5^{84}.                                  \tag{7}
\]

Both adjacency and the strict inequality \(n>\alpha m\) are checked with
integers. Thus for real interval starts the sharp coverage infimum for
\([x,\alpha x)\) is \(n/\alpha\): every \(x>n/\alpha\) is covered,
while \(x=n/\alpha\) is not. For integer starts the exact threshold here is

\[
 M=\left\lfloor\frac{1000000\,n}{1000001}\right\rfloor+1.
\]

Every integer \(x\geq M\) has a 5-smooth integer in \([x,\alpha x)\);
the interval at \(x=M-1\) is empty. The additional condition
\(\lfloor n/\alpha\rfloor>m\) is verified. In other cases this condition
can fail: a real exceptional gap need not contain an integer interval
start. The checker then reports an upper bound and leaves the sharp
integer threshold unspecified.

## Verification and reproduction

From the repository root, run:

```bash
python3 -B scripts/run_smooth_gap_certificate.py --suite --den 1000000 \
  --out results/research/smooth_gap_certificate_checks.json
python3 -B scripts/run_smooth_gap_certificate.py \
  --verify-file results/research/smooth_gap_certificate_checks.json
./build.sh
python3 tests/test_smoke.py
```

The certificate verifier does not trust the constructed frontier arrays.
It checks the positive near-product witnesses, tests every uncovered axis
and face point independently, counts the complement by inclusion-exclusion,
and finds its maximal points from a separate corner characterization:
every coordinate of a maximal standard monomial is one less than a positive
coordinate of some minimal generator. An independent two-outer-loop interval
test proves that (7) is consecutive and counts its exact starting rank.

The suite also compares every gap membership against exact sequential DP
for 30 varied-prime controls and a larger \(\alpha=1001/1000\) control.
Deleting a necessary generator and corrupting a positive witness must each
cause rejection. Assertions must remain enabled; do not run with Python's
optimization flags. These checks concern the stated finite certificates;
the uniform face identity and complexity bound have the proofs above.

## Primary comparisons and unresolved questions

* [Tijdeman--Meijer, 1974](https://www.numdam.org/item/CM_1974__29_3_273_0/),
  pp. 273--279: division of adjacent terms by a common factor is already
  explicit, and Theorem 1 characterizes two-generator ratios by one-sided
  convergents. The introduction and surrounding statements were directly
  inspected. Identity (1) and the certificate construction use that elementary
  division principle. They do not settle the paper's prescribed-support
  infinitude conjecture for larger prime sets or supply the requested
  higher-dimensional continued-fraction characterization.
* [Berczes--Dujella--Hajdu, 2014](https://math.unideb.hu/sites/default/files/inline-files/sunit_submit.pdf),
  Theorems 2.1--2.3: the log-hyperplane formulation, effective general
  neighbour determination, and efficient two-prime neighbour algorithms are
  established prior art. The introductory and theorem statements were
  directly inspected; no full-proof or all-literature audit is claimed.
* [Languasco--Luca--Moree--Togbe](https://arxiv.org/abs/2309.12806),
  Sections 3--4 of the inspected arXiv version: multiplicative interval
  coverage and computation of the last exceptional integer start are
  established two-prime problems. The three-prime certificate is compared
  with these tasks, without claiming that the task itself is new.

The useful next questions concern how many minimal gap witnesses can occur,
whether their changes admit stronger uniform bounds, and what those changes
say about the frequency and size of relative gaps. No asymptotic estimate
for \(U_i\) as \(\alpha\downarrow1\), new Diophantine approximation bound,
zeta-zero theorem, or new gap-distribution limit is established here.
