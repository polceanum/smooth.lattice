# Certified generic hook fans in every characteristic

For the flat square-zero extensions of
\[
 M_t=(y^{t^2+1},y^tz^t,z^{t^2+1}),
 \qquad n=4t^3-2t^2+4t,
\]
the complete **generic** fan counts at thickness two and three are
\[
             g_{\rm gen}(M_2)=89,\qquad g_{\rm gen}(M_3)=948.       \tag{1}
\]
Both counts hold in every characteristic. Here generic means on a nonempty
Zariski open subset of the independent extension parameters, over an
infinite field or a sufficiently large finite extension. It does not mean every coefficient choice over a given finite
field. The earlier [zero-one fixtures](square_zero_extension_fans.md#exact-arithmetic-calibration-of-the-compatible-hook-model)
retain their different counts, including 709 in characteristic zero
and 546 in characteristic two at length 102.

These are exact finite calibrations of the unresolved compatible-minor
route. They do not prove any law as \(t\) grows, give a polynomial
separation above the generic corner-cut scale, improve unrestricted
\(G_3\), or establish novelty. In particular a support-matching fan is
not an upper bound on the number of cones after arbitrary coefficient
specialization.

## The independent parameters

Put \(a=t^2+1\), and let \(D\) be the standard diagram of \(M_t\).
The images \(f_0,f_1,f_2\) of its three minimal generators can be chosen
independently with supports
\[
\begin{array}{c|c}
f_0&\{(i,j)\in D:i\ge t\text{ or }j\ge a-t\}\\
f_1&\{(i,j)\in D:i\ge t\text{ or }j\ge t\}\\
f_2&\{(i,j)\in D:i\ge a-t\text{ or }j\ge t\}.
\end{array}                                                       \tag{2}
\]
These are precisely the conditions from the adjacent syzygies. For
example, the two sides of \(z^tf_0=y^{a-t}f_1\) have disjoint possible
nonzero monomial supports, since \(a\ge2t\); hence each vanishes.
The other adjacent syzygy is symmetric. This proves (2) over every
field and over the integral polynomial ring in the parameters.
There are \(2|D|=n\) independent parameters.

The [existing determinant model](square_zero_extension_fans.md#the-exact-minor-model-in-a-monomial-flat-fiber)
contracts the compulsory plain diagram \(D\), leaving
\([\operatorname{Id}_D\mid\Phi]\). Every nonzero entry of \(\Phi\)
is one of the parameters in (2). Its repeated appearances encode
multiplication; they are not replaced by independent matrix entries.

## A support relaxation and a survival criterion

Keep the nonzero entry positions and take the corresponding transversal
matroid: a set of columns is independent when it can be matched to
distinct rows. Every actual matrix basis is a basis of this relaxation,
by determinant expansion. The relaxation itself is always a matroid;
it is not asserted to arise from an algebra.

Let \(\mathcal B\) be its complete set of positively coherent bases,
with the original exponent labels. For \(B\in\mathcal B\), eliminate
its selected identity columns and write \(\Delta_B\) for the remaining
square minor of \(\Phi\), as a polynomial over \(\mathbb Z\).
If every \(\Delta_B\) is nonzero over a coefficient field, the actual
and relaxed **entire** positive fans agree on the open parameter locus
\[
                         \prod_{B\in\mathcal B}\Delta_B\ne0.       \tag{3}
\]
Indeed the relaxation's unique minimum at any generic positive weight
is available in the actual matrix. Every actual basis was already a
relaxed basis, so neither a different minimum nor an extra cone can
occur. This argument certifies whole cones and their labels, rather
than selected weights. It requires the complete relaxed fan.

For the two fixtures the checker traverses that whole fan using exact
matching exchanges and rational polygons. Every relaxed optimum is a
monomial downset, and all bases in the supplied certificate are
recovered. The support cones cover the positive projective triangle
with total area one.

## Excluding every content prime exactly

The following elementary certificate avoids inferring characteristic
zero or generic behavior from a large-prime computation. If
\(\Delta\in\mathbb Z[c_1,\ldots,c_s]\), its coefficient content divides
every integer value \(\Delta(v)\). Let a finite collection of exact
integer evaluations have nonzero gcd \(d\). Only primes dividing
\(d\) can make \(\Delta\) identically zero after reduction.
For each such prime \(p\), one nonzero evaluation in any finite
extension of \(\mathbb F_p\) proves that this reduction is not the
zero polynomial. Excluding all these primes proves that \(\Delta\)
is primitive. The selected evaluations are explicit witnesses; no
probability of correctness is used.

The certificate applies this argument to every coherent minor:

| Length | Coherent minors | Sum of degrees | Largest degree | Content-prime exclusions |
|---:|---:|---:|---:|---:|
| 32 | 89 | 664 | 12 | 63 |
| 102 | 948 | 26084 | 42 | 730 |

The possible primes left by the integer gcds are \(2,5,7\) at length
32 and \(2,3,5,7,11,31\) at length 102. Exact evaluations in the
reported finite fields exclude every one for every relevant minor.
The field constructions are checked by verifying that every nonzero
element of the stated polynomial quotient is invertible.
All 1037 coherent determinant polynomials are therefore primitive.
Equation (3) is nonempty over an algebraic closure in every
characteristic, proving (1). This conclusion concerns the generic
parameter locus; it preserves all exceptional specializations.

For a finite starting field \(\mathbb F_q\), an extension with
\(q^e>664\), respectively \(q^e>26084\), suffices for simultaneous
nonvanishing. The product in (3) has at most that total degree, so the
classical polynomial grid bound gives a nonzero value on the larger
field. This is an existence bound for a common specialization, not a
claim that the separate minor witnesses already form one.

## Reproduction and remaining problem

```bash
python3 -B scripts/run_generic_hook_checks.py
./build.sh
python3 tests/test_smoke.py
```

The [certificate and verification report](../results/research/generic_hook_checks.json)
contains every coherent basis, all integer evaluation data, finite-field
moduli and nonzero determinant witnesses. The checker reconstructs both
whole support fans, replays the arithmetic certificates, and separately
checks original commuting multiplication and all ground columns.
Integer determinants use established fraction-free elimination; selected
ones are checked by independent exact rational Gaussian elimination.
The polynomial grid estimate is the classical input in
[Schwartz, Section 1, Lemma 1 and Corollary 1](https://doi.org/10.1145/322217.322225).
None of matching, generic specialization, determinant content or finite
field arithmetic is claimed as new.

A growing-thickness proof must still establish which coherent support
bases survive symbolically and then count their exposure uniformly.
Alternatively it must retain the further rank constraints when they
fail. The tests above cannot supply that theorem. Nongeneric coefficient
strata remain part of the original all-quotient extremal problem even
if a generic rank theorem is obtained. The unrestricted target and
precise current-priority assessment remain open.
