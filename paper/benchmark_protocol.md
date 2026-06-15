# Benchmark Protocol Notes

## Priority Target: X+Y vs Layer-Compressed k=5

The next defensible comparison is the first-five-primes fixed-rank task:

```text
P = (2,3,5,7,11)
N = 10^12
```

The target claim to test is deliberately narrow:

```text
For fixed P=(2,3,5,7,11) and N=10^12, the layer-compressed
rank/unrank solver is faster and/or lower-memory than a practical
adaptive Cartesian-sum X+Y selector, while returning an exponent vector
certified by the independent interval auditor.
```

This should not be generalized to all smooth-number algorithms or all
fixed-prime rank/unrank tasks without more comparisons.

## Harness

Use:

```bash
python3 scripts/run_xplusy_vs_layer5.py
```

Default outputs are written to:

```text
results/local/xplusy_vs_layer5_<timestamp>/
```

Expected files:

- `report.json`: machine-readable full provenance and raw command results.
- `summary.csv`: compact table suitable for later aggregation.
- `report.md`: human-readable benchmark summary.

The harness records:

- git commit and dirty status;
- compiler path, version, and flags;
- platform, CPU, and memory when available;
- exact command lines;
- stdout/stderr and return codes;
- wall time and peak RSS when `/usr/bin/time` supports it;
- parsed method metrics such as reported seconds, band size, exponent vector,
  and certification fields.

## Compared Methods

### `xplusy_adaptive`

Command shape:

```bash
bin/smooth_xplusy_baseline bench 2,3,5,7,11 1000000000000 0
```

This is a practical adaptive Cartesian-sum value-selection baseline over a
materialized meet-in-the-middle split.

Important caveat: this row currently measures selection of the target log value,
not full exponent-vector reconstruction. Therefore, if the layer-compressed full
unrank solver beats this row, the result is conservative. If it loses, the loss
must be reported plainly.

### `layer_compressed`

Command shape:

```bash
bin/smooth_layer_compressed_general nth 2,3,5,7,11 1000000000000 50000
```

This solver returns the exponent vector and exact-value digit count. It uses
floating-log search/counting internally, so the raw solver result is experimental
until independently audited.

### `interval_audit`

Command shape:

```bash
bin/smooth_interval_audit_exps_k6 2,3,5,7,11 <layer_exps> 1000000000000
```

This certifies the returned exponent vector by independently checking:

```text
count_le(P, value(exps)) == N
```

Only rows with `rank_certified=true` should be used as certified correctness
evidence.

## Paper Reporting Rules

For any table derived from this harness, report:

- whether `xplusy_adaptive` is selection-only;
- whether `layer_compressed` was independently audited;
- wall time and peak RSS, not only internal solver timers;
- compiler and flags;
- machine metadata;
- raw artifact path and git commit;
- negative results, including cases where `X+Y` wins or certification fails.

Do not write "best known" from this comparison alone. A safe phrasing, if the
data supports it, is:

```text
On this fixed five-prime random-access task, the layer-compressed solver
outperforms our practical adaptive X+Y materialized-MITM baseline under the
recorded hardware/compiler conditions, and the returned rank is independently
certified.
```

## Certified Five-Prime Suite

The broader paper claim should be tested on all five-prime subsets of
`{2,3,5,7,11,13}` at `N=10^12`:

```text
(2,3,5,7,11)
(2,3,5,7,13)
(2,3,5,11,13)
(2,3,7,11,13)
(2,5,7,11,13)
(3,5,7,11,13)
```

Use:

```bash
python3 scripts/run_five_prime_suite.py
```

Clean suite artifact:

```text
results/benchmarks/five_prime_suite_1e12/
```

Observed result at commit `87da37090b939d217a5ee2a51e2c08101d5d13ac` on the
recorded macOS/x86_64 Apple-clang machine:

- 6/6 cases completed.
- 6/6 returned exponent vectors were independently certified.
- Layer-compressed full unranking won 6/6 wall-time comparisons against
  adaptive materialized `X+Y` value selection.
- X+Y/layer wall-time ratio ranged from 1.1868 to 1.9125, with mean 1.4963.

The resulting paper claim is:

```text
Across all six five-prime subsets of {2,3,5,7,11,13} at N=10^12,
layer-compressed full unranking outperformed our practical adaptive
materialized X+Y value-selection baseline under the recorded hardware/compiler
conditions, and every returned exponent vector was independently certified.
```

This is still not a "best known" claim. It is a certified benchmark-suite win
against one serious practical comparator family.

## Full Materialized X+Y Unrank Suite

The next stricter comparator is `smooth_xplusy_full_unrank`, a materialized
meet-in-the-middle baseline that returns an exponent vector rather than only
selecting a target log value. It:

- materializes both split sides with log sums and packed exponent vectors;
- uses adaptive rank counting to bracket the target rank;
- enumerates a final candidate band;
- exact-sorts that band by multiprecision integer value;
- reports the selected exponent vector.

Use:

```bash
python3 scripts/run_full_xplusy_suite.py
```

Default outputs are written to:

```text
results/local/full_xplusy_suite_<timestamp>/
```

The harness compares `xplusy_full_unrank` with `layer_compressed` on the same
six five-prime subsets at `N=10^12`. It audits both returned exponent vectors
with `smooth_interval_audit_exps_k6` and records whether both methods returned
the same vector.

This supports a cleaner paper comparison than the value-selection baseline
because both rows perform the same top-level task:

```text
unrank(P, N) -> exponent vector.
```

A safe claim, if the clean suite data supports it, is:

```text
Across all six five-prime subsets of {2,3,5,7,11,13} at N=10^12,
layer-compressed full unranking outperformed our practical full materialized
X+Y unrank baseline under the recorded hardware/compiler conditions. Both
methods returned the same exponent vector in every case, and both outputs were
independently rank-certified.
```

Clean suite artifact:

```text
results/benchmarks/full_xplusy_suite_1e12/
```

Observed result at commit `5314fd8a56fb7ca760046076e1f7f168fd48386a` on the
recorded macOS/x86_64 Apple-clang machine:

- 6/6 cases completed.
- 6/6 cases returned matching exponent vectors.
- 6/6 layer outputs and 6/6 full-X+Y outputs were independently certified.
- Layer-compressed full unranking won 6/6 wall-time comparisons.
- Full-X+Y/layer wall-time ratio ranged from 1.9672 to 2.7248, with mean
  2.3273.
- The layer solver used the analytic asymptotic bracket in all six cases.

This supports the safe claim above for this benchmark suite and recorded
machine/compiler. It is still not, by itself, a broad "best known" claim.

## Certified k=8 Sums-Only Instance

The first paper-grade higher-k target is:

```text
P = (2,3,5,7,11,13,17,19)
N = 10^12
```

Use:

```bash
python3 scripts/run_k8_certificate.py
```

Default outputs are written to:

```text
results/local/k8_certificate_<timestamp>/
```

The harness runs the exploratory high-k sums-only MITM unranker, then audits the
returned exponent vector with the independent interval-log rank auditor. The
auditor uses integer log intervals and exact big-integer resolution of boundary
ambiguities; it does not reuse the sums-only solver's floating reconstruction
path.

A safe claim, if the clean artifact supports it, is:

```text
For P=(2,3,5,7,11,13,17,19) and N=10^12, the high-k sums-only
solver returned an exponent vector whose rank was independently certified by
the interval-log auditor.
```

This is a correctness/certification claim for a fixed higher-dimensional
instance. It is not a broad performance or best-known claim.

## Analytic Count Residual Probe

The next analytic-oracle diagnostic evaluates the asymptotic lattice-count
approximation at certified target exponent vectors:

```bash
python3 scripts/run_analytic_count_probe.py
```

Default outputs are written to:

```text
results/local/analytic_count_probe_<timestamp>/
```

The harness records, for each target vector:

- the analytic count estimate from the truncated generating-function expansion;
- the analytic derivative at the same log height;
- the current layer counter's count at that log height;
- the residual against the independently certified rank `N`.

This is deliberately not yet an algorithmic claim. It is a measurement stage for
the proposed analytic-lattice count oracle. A future oracle would need an error
envelope and exact boundary correction before it can replace exact counting in a
published solver.

Clean probe artifact:

```text
results/benchmarks/analytic_count_probe_1e12/
```

Observed result at commit `f244aeafa0235adbb5f9441d4ff8c01b84e938e4`:

- 8/8 probes completed.
- The six k=5 residuals were within 167 ranks at `N=10^12`.
- The k=6 first-six-primes residual was about 453 ranks.
- The k=8 first-eight-primes residual was about 1.98e6 ranks, or 2e-6
  relative.
- The floating-log layer counter matched the certified rank in 6/8 probes and
  was off by one on two boundary k=5 probes.

This supports studying analytic boundary-correction bands, but it also
reinforces that exact-boundary correctness claims need the interval auditor.

## Analytic Boundary-Correction Band Probe

The next prototype uses the analytic estimate as an actual band center:

```bash
python3 scripts/run_analytic_band_probe.py
```

For each target rank, the harness:

1. solves the analytic count equation for `T`;
2. computes `half_width = rank_radius / analytic_derivative(T)`;
3. checks the endpoint counts with the current layer counter;
4. enumerates and exact-sorts the band only if its endpoint count is below the
   configured candidate cap.

This tests the intended oracle shape:

```text
analytic count chooses a small correction band;
exact enumeration recovers the target only inside that band.
```

The current form is still experimental. Endpoint counts use the floating-log
layer counter, and recovered vectors still need independent interval auditing
before supporting correctness claims.

Clean band-probe artifact:

```text
results/benchmarks/analytic_band_probe_1e12/
```

Observed result at commit `4bd49bec50aab80c87e599438f325f55814e6222`:

- 8/8 targets were inside the analytic band.
- 7/8 cases were enumerated under the configured cap and recovered the expected
  exponent vector.
- The six k=5 cases used radius 1000 and produced bands of 1980-2032
  candidates.
- The k=6 first-six-primes case used radius 2000 and produced a band of 4005
  candidates.
- The k=8 first-eight-primes case used radius 2.5e6 and contained the target,
  but the band had 5,000,120 candidates and was not enumerated under the
  200,000-candidate cap.

This is the first positive evidence for the analytic-oracle architecture at
k=5/k=6 and the first concrete bottleneck for k=8: the residual model is good
enough to contain the target, but not yet sharp enough to make high-k boundary
enumeration cheap.

## Residual-Corrected Analytic Band Probe

The next probe asks whether one exact count at the analytic center is enough to
remove most of the residual:

```bash
python3 scripts/run_analytic_band_corrected_probe.py
```

For each target rank, the harness:

1. solves the analytic count equation for an initial center `T0`;
2. computes the exact layer count `C(T0)`;
3. shifts the center by `(N - C(T0)) / analytic_derivative(T0)`;
4. checks a smaller rank-radius band around the corrected center;
5. enumerates and exact-sorts that band if it is under the configured cap.

This tests the more ambitious oracle shape:

```text
analytic count predicts the center;
one exact count estimates the residual;
exact enumeration handles only the corrected boundary band.
```

Clean corrected-band artifact:

```text
results/benchmarks/analytic_band_corrected_probe_1e12/
```

Observed result at commit `f965d9cb677d26aa83b7b35c8412a90f856c8fc2`:

- 8/8 targets were inside the corrected analytic band.
- 8/8 cases were enumerated under the 200,000-candidate cap and recovered the
  expected exponent vector.
- The six k=5 cases used radius 100 and produced bands of 193-204 candidates.
- The k=6 first-six-primes case used radius 500 and produced a band of 997
  candidates.
- The k=8 first-eight-primes case used radius 10,000 and produced a band of
  19,853 candidates, compared with the previous unenumerated 5,000,120-candidate
  uncorrected band.

This is stronger evidence for the analytic-oracle architecture, especially at
k=8. It is still not a proof: endpoint counts use the floating-log layer
counter, and the recovered vectors must be independently interval-audited before
supporting correctness claims.

The audit-aware version of the same harness runs
`bin/smooth_interval_audit_exps` on every recovered vector:

```bash
python3 scripts/run_analytic_band_corrected_probe.py \
  --out-dir results/benchmarks/analytic_band_corrected_certified_probe_1e12
```

Clean certified corrected-band artifact:

```text
results/benchmarks/analytic_band_corrected_certified_probe_1e12/
```

Observed result at commit `11687a9ef1b183ff88ce3fea53bca4c281315228`:

- 8/8 targets were inside the corrected analytic band.
- 8/8 cases were enumerated under the 200,000-candidate cap.
- 8/8 recovered the expected exponent vector.
- 8/8 recovered vectors were independently rank-certified by the interval-log
  auditor with `certified_count_le=N`.
- The k=8 first-eight-primes case used radius 10,000, had a 19,853-candidate
  corrected band, and certified in 18.2 wall seconds on the benchmark machine.

This is now certification-grade evidence for this fixed target suite: the
analytic and floating-log machinery propose and enumerate the candidate, while
the independent interval auditor supplies the rank certificate. It still does
not establish a general analytic error bound or a broad best-known claim.

## Corrected-Oracle Random Suite

The next broader suite samples deterministic pseudo-random prime subsets from
the current audited prime universe `{2,3,5,7,11,13,17,19}` and compares three
full unrank paths:

1. residual-corrected analytic band plus interval audit;
2. sums-only MITM;
3. full materialized X+Y unrank.

The default command is:

```bash
python3 scripts/run_corrected_oracle_random_suite.py \
  --out-dir results/benchmarks/corrected_oracle_random_suite_1e9_1e12
```

The default seed is `20260609`. The default suite uses k=5, k=6, and k=8 cases
at `N=10^9` and `N=10^12`, taking two pseudo-random subsets per `(k,N)` when
that many subsets exist. The k=8 row is necessarily the full audited universe.

Pass criteria for a row are deliberately strict:

- the corrected analytic band must recover an exponent vector;
- the interval auditor must certify that vector with `certified_count_le=N`;
- sums-only MITM and full materialized X+Y must return the same exponent vector.

Clean random-suite artifact:

```text
results/benchmarks/corrected_oracle_random_suite_1e9_1e12/
```

Observed result at commit `e53d69b364d568a49e42784892289a15aeb58742`:

- 10/10 rows completed.
- 10/10 corrected analytic oracle outputs were independently rank-certified.
- 10/10 sums-only and full materialized X+Y outputs matched the certified
  corrected-oracle vector.
- Corrected oracle used less peak RSS than full materialized X+Y in 10/10 rows,
  with mean full-X+Y/corrected RSS ratio about 7.90.
- Wall time was mixed: corrected oracle won 4/10 rows against full materialized
  X+Y and 3/10 against sums-only MITM. The large k=5, `N=10^12` cases favor the
  corrected oracle strongly; k=6/k=8 currently favor the simpler MITM baselines
  in wall time.

This is broader certified evidence that the residual-corrected oracle is a real
algorithmic path, and it suggests a narrow memory-efficiency claim against full
materialized X+Y. It is also a negative speed result for the current k=6/k=8
implementation, so the paper should not present the corrected oracle as a
general wall-time winner.

## Residual-Corrected Full X+Y Speed Suite

The corrected-oracle random suite compares the layer-compressed corrected
oracle against other methods. To isolate the speed contribution of residual
correction itself, the full materialized X+Y benchmark also has a corrected
mode:

```bash
bin/smooth_xplusy_full_unrank analytic-band-corrected \
  primes_csv N rank_radius max_candidates
```

This mode uses the same materialized X+Y groups and exponent packing as
`bin/smooth_xplusy_full_unrank nth`, but replaces the adaptive
interpolation/bisection count loop with:

1. analytic center `T0`;
2. one exact full-X+Y count `C(T0)`;
3. residual shift `(N-C(T0))/analytic_derivative(T0)`;
4. a small corrected boundary band, exact-sorted to recover the exponent vector.

The reproducible suite is:

```bash
python3 scripts/run_xplusy_corrected_speed_suite.py \
  --out-dir results/benchmarks/xplusy_corrected_speed_suite_1e9_1e12
```

Clean speed artifact:

```text
results/benchmarks/xplusy_corrected_speed_suite_1e9_1e12/
```

Observed result at commit `4a581fe2fb495997442c4e233fc01ad172a2d4df`:

- 10/10 rows completed.
- 10/10 corrected-mode outputs were independently rank-certified and matched
  the adaptive full-X+Y exponent vector.
- Corrected mode won 9/10 wall-time comparisons.
- Corrected mode won 10/10 reported in-process time comparisons.
- Mean adaptive/corrected wall-time ratio was about 1.33.
- Mean adaptive/corrected reported-time ratio was about 1.82.
- Corrected mode used less peak RSS in 10/10 rows, with mean adaptive/corrected
  RSS ratio about 1.34.
- The one wall-time loss was a small `k=5, N=10^9` case where process/cache
  overhead dominated; the corrected mode still won that row on reported
  in-process time.

This supports a narrow speed claim: residual correction accelerates this
adaptive full materialized X+Y unrank implementation while preserving
independent certification of the corrected output. It does not imply a speed
win over sums-only MITM, layer-compressed MITM, or every known X+Y selection
algorithm.

## Residual-Corrected Sums-Only MITM Speed Suite

The same residual-correction idea is also implemented inside the sums-only
MITM kernel:

```bash
bin/smooth_sums_only_scalable analytic-band-corrected \
  primes_csv N rank_radius max_candidates
```

This mode uses the same sums-only MITM split and reconstruction machinery as
`bin/smooth_sums_only_scalable nth`, but replaces the adaptive
interpolation/bisection search with:

1. analytic center `T0`;
2. one exact sums-only MITM count `C(T0)`;
3. residual shift `(N-C(T0))/analytic_derivative(T0)`;
4. a corrected boundary band, sorted to recover the exponent vector.

The reproducible suite is:

```bash
python3 scripts/run_sums_corrected_speed_suite.py \
  --out-dir results/benchmarks/sums_corrected_speed_suite_1e9_1e12
```

Clean speed artifact:

```text
results/benchmarks/sums_corrected_speed_suite_1e9_1e12/
```

Observed result at commit `be4e24a8f4adf0947e568dcc57f4797b0b35ca60`:

- 10/10 rows completed.
- 10/10 corrected-mode outputs were independently rank-certified and matched
  the adaptive sums-only MITM exponent vector.
- Corrected mode won 9/10 wall-time comparisons.
- Corrected mode won 10/10 reported in-process time comparisons.
- Corrected mode used less peak RSS in 10/10 rows.
- Mean adaptive/corrected wall-time ratio was about 1.33.
- Mean adaptive/corrected reported-time ratio was about 1.61.
- Mean adaptive/corrected RSS ratio was about 1.14.
- The one wall-time loss was a small `k=5, N=10^9` row where the corrected
  process wall time was anomalously larger despite winning reported in-process
  time.

This supports a narrow speed claim against the prior adaptive sums-only MITM
path. Together with the full-X+Y corrected suite, it shows that the residual
correction is not merely helping one data representation. It still does not
establish a broad best-known theorem or a comparison against every published
selection method.

## Certified Portfolio Suite

Single-method suites are useful for isolating implementation effects, but they
do not answer the practical question "which implemented method should be used
for this case?" The portfolio harness runs all applicable full-unrank methods,
requires successful methods to agree on the exponent vector, audits that vector
when possible, and records the fastest wall-time and reported-time methods:

```bash
python3 scripts/run_certified_portfolio_suite.py \
  --out-dir results/benchmarks/certified_portfolio_suite_1e9_1e12
```

Clean portfolio artifact:

```text
results/benchmarks/certified_portfolio_suite_1e9_1e12/
```

Observed result at commit `b8357ac6e09c27b781892b9a304d72c184d2abae`:

- 14/14 rows had solver agreement among successful implemented methods.
- 14/14 rows were independently interval-certified.
- There were no solver-disagreement rows and no audit-blocked rows.
- Wall-time winners were: `xplusy_corrected` in 5 rows, `layer_corrected` in 4
  rows, `beatty3` in 3 rows, and `xplusy_adaptive` in 2 rows.
- Reported-time winners were: `xplusy_corrected` in 6 rows, `layer_corrected` in
  4 rows, and `beatty3` in 4 rows.

This is the current implemented-method scoreboard. It should be the default
artifact cited when discussing practical fastest-certified behavior in this
repository. It still does not justify broad best-known language: the remaining
gates are a faithful soft-heap/Frederickson-Johnson-style `X+Y` comparison,
clearer Barvinok/fixed-dimensional lattice-count positioning, and arbitrary-prime
certification beyond the current fixed k<=12 audited prime universe.

## Heap Frontier Baseline Suite

The heap/frontier baseline is a canonical sequential generator. It stores a
trace parent for each generated exponent vector and expands children only in
nondecreasing prime-index order; this gives each exponent vector one path while
preserving value-order extraction by a log-keyed priority queue. It returns a
full exponent vector, not just a selected value.

Default small-to-medium suite:

```bash
python3 scripts/run_heap_frontier_baseline_suite.py \
  --out-dir results/benchmarks/heap_frontier_baseline_suite_1e5_1e6
```

Stress suite:

```bash
python3 scripts/run_heap_frontier_baseline_suite.py \
  --out-dir results/benchmarks/heap_frontier_baseline_suite_1e7 \
  --dp-max-n 10000000 \
  --case k3_N1e7:2,3,5:10000000 \
  --case k5_N1e7:2,3,5,7,11:10000000 \
  --case k6_N1e7:2,3,5,7,11,13:10000000 \
  --case k8_N1e7:2,3,5,7,11,13,17,19:10000000
```

The harness compares heap/frontier generation with the DP pointer baseline when
enabled and with the current solver (`beatty3` for k=3, `sums_adaptive`
otherwise). A row is certified only when successful methods agree on the
exponent vector and the independent interval auditor proves
`count_le(candidate)=N`.

Clean default artifact:

```text
results/benchmarks/heap_frontier_baseline_suite_1e5_1e6/
```

Observed result at commit `b26e3cdf865b15754d8b189dafe4e7cc36e19f80`:

- 8/8 rows over k=3,5,6,8 and N=10^5,10^6 had solver agreement.
- 8/8 rows were independently interval-certified.
- Wall-time winners were current_beatty3 in 2 rows, current_sums_adaptive in 5
  rows, and dp_pointer in 1 row.
- Mean heap/current wall-time ratio was 4.5654.
- Mean heap/DP wall-time ratio was 1.9641.

Clean stress artifact:

```text
results/benchmarks/heap_frontier_baseline_suite_1e7/
```

Observed result at commit `6526d68a0adcd45b0041a1dcfe8da1c6a9fe1f56`:

- 4/4 rows over k=3,5,6,8 at N=10^7 had solver agreement.
- 4/4 rows were independently interval-certified.
- Current solvers won all 4 wall-time comparisons.
- Mean heap/current wall-time ratio was 76.5315.
- Mean heap/DP wall-time ratio was 2.9536.

This is a legitimate sequential-frontier comparison. It is not evidence against
soft-heap `X+Y`, Frederickson-Johnson sorted-matrix selection, or
Barvinok-style fixed-dimensional lattice counting.

## Analytic-Bracket Layer Hybrid

The layer-compressed solver now uses a non-MITM analytic seed for large-rank
initial bracketing. The seed comes from the small-`s` expansion of

```text
product_i (1 - exp(-s log p_i))^-1 / s,
```

which is the Laplace transform of the fixed-prime count function. This produces
an asymptotic estimate for the log height `T` satisfying `count_le(P,T) ~= N`.

The approximation is not used as proof. It only proposes a bracket; the solver
still checks and repairs the bracket with its exact layer-count routine, and the
final exponent vector must still be independently audited before being called
certified. Benchmark artifacts record `leading_seed`, `analytic_seed`, and
`analytic_bracket`.

In the updated clean full-X+Y suite artifact, the analytic-bracket solver:

- used analytic bracketing in 6/6 cases;
- kept final rank gaps between 334 and 366;
- reduced the mean final layer band from 11286.8 in the previous full-X+Y
  artifact to 345.8;
- improved reported layer time in 5/6 cases, with mean reported layer time
  changing from 0.8093s to 0.7443s.

These are optimization results, not independent correctness evidence.

## Sorted-Matrix / LOH Workbench

The exploratory sorted-matrix workbench probes two adjacent comparator ideas:

- a range-pruned block counter over the sorted Cartesian-sum matrix;
- a Mirzaian-Arjomandi square sorted-matrix value-selection probe, adapted to
  rectangular MITM splits by padding the shorter side;
- a Kaplan/Frederickson-Johnson-style Mat-Select2 exponential-block row-sorted
  selector using an exact binary heap for the Mat-Select1 primitive;
- an LOH-style output top-k probe at a capped rank.

Use:

```bash
python3 scripts/run_sorted_matrix_workbench.py
```

Default outputs are written to:

```text
results/local/sorted_matrix_workbench_<timestamp>/
```

The default cases are the same six five-prime subsets of
`{2,3,5,7,11,13}` at `N=10^12`. The report compares the internal time of the
range-pruned block counter against the ordinary linear saddleback count used by
the adaptive X+Y selector.

Important caveats:

- The range-pruned block counter is sorted-matrix/FJ-style diagnostic code, not
  a faithful Frederickson-Johnson selection implementation.
- The Mirzaian-Arjomandi row selects a log value only. It does not reconstruct
  an exponent vector and does not certify a rank by itself.
- The Mat-Select2 row implements the exponential-block row-sorted selector, but
  the primitive that is soft-heap based in the asymptotic paper algorithm is
  implemented here with an exact binary heap. It is therefore a serious
  algorithmic bridge, not a full soft-heap time-bound comparison.
- The LOH row is an output-style top-k probe. When `N_probe != N`, it is not a
  random-access unranking comparator for the full target rank.
- These results should be reported as diagnostic or negative evidence unless a
  future patch implements the full published algorithm being invoked.

Clean workbench artifact:

```text
results/benchmarks/sorted_matrix_workbench_1e12/
```

Observed result at commit `17a7e40ae790ae9017ab8274204eaaf107add212` on the
recorded macOS/x86_64 Apple-clang machine:

- 6/6 cases completed.
- The Mirzaian-Arjomandi value selector passed exhaustive small validation:
  5185/5185 cases, 0 failures, maximum delta 0.
- The range-pruned block counter beat the ordinary linear saddleback count in
  1/6 cases.
- That win was narrow: block/linear ratio 0.9115.
- The mean block/linear internal time ratio was 1.3467, so the probe was slower
  on average.
- The Mirzaian-Arjomandi value selector matched the adaptive selected log in
  6/6 cases, but won 0/6 timing comparisons. Its mean MA/linear internal time
  ratio was 9.4042.
- The LOH row used `N_probe=10^6`, not the full target rank `N=10^12`.

This is a validated negative/mixed result. The Mirzaian-Arjomandi row is a
faithful value-selection comparator for the sorted-matrix problem, and it does
not strengthen the practical X+Y baseline at the headline five-prime target
scale. It still does not discharge the full Frederickson-Johnson or soft-heap
comparison obligations for full unranking.

Mat-Select2 heap-primitive first-k artifact:

```bash
python3 scripts/run_sorted_matrix_workbench.py \
  --out-dir results/benchmarks/sorted_matrix_matselect2_heap_firstk_1e12 \
  --N 1000000000000 \
  --prime-set 2,3,5,7,11 \
  --prime-set 2,3,5,7,11,13 \
  --prime-set 2,3,5,7,11,13,17,19
```

Observed result:

- Mat-Select2 heap-primitive exhaustive validation passed 3233/3233 cases.
- The first k=5 target matched the selected log exactly.
- It won 0/1 comparable timing rows; Mat-Select2/current linear ratio was
  12.7823.
- The first k=6 and k=8 rows were skipped by the active-row cap.

This is negative evidence for the exact-heap Mat-Select2 bridge. It also gives
the next implementation target a precise shape: replace the exact heap primitive
with a real soft heap before making any soft-heap claim.

## Mirzaian-Arjomandi Full-Unrank Comparator

The sorted-matrix workbench above tests value selection only. The full-unrank
comparator wraps the Mirzaian-Arjomandi sorted-matrix selector in the same
reconstruction discipline used by the materialized X+Y unranker:

1. build both MITM sides with exponent packs;
2. select the target Cartesian-sum log value with the Mirzaian-Arjomandi
   square sorted-matrix selector, padding the shorter rectangular side;
3. open a narrow analytic-width band around that selected value;
4. enumerate candidate exponent vectors in the band;
5. exact-sort the candidate band by multiprecision integer value;
6. independently audit the returned exponent vector.

Use:

```bash
python3 scripts/run_ma_full_unrank_suite.py \
  --out-dir results/benchmarks/ma_full_unrank_suite_1e12
```

Default outputs are written to:

```text
results/local/ma_full_unrank_suite_<timestamp>/
```

The default benchmark compares `xplusy_ma_full` against
`xplusy_corrected` on the six five-prime subsets of `{2,3,5,7,11,13}` at
`N=10^12`, with `rank_radius=250`, `max_candidates=200000`,
`ma_max_n=30000000`, and `ma_max_middle=200000000`.

Clean comparator artifact:

```text
results/benchmarks/ma_full_unrank_suite_1e12/
```

Observed result on the recorded macOS/x86_64 Apple-clang machine:

- 6/6 cases completed.
- 6/6 returned the same exponent vector as analytic-corrected X+Y.
- 6/6 matching vectors were independently interval-rank certified.
- The Mirzaian-Arjomandi full-unrank path won 0/6 wall-time comparisons.
- The MA/corrected wall-time ratio ranged from 4.4781 to 7.1840, with mean
  5.7000.
- Peak RSS for the MA path ranged from about 1.43GB to 2.02GB on these rows,
  versus about 0.36GB to 0.58GB for analytic-corrected X+Y.

This is a rigorous negative result for the implemented MA full-unrank wrapper:
it is a correctness-complete comparator, but it is not a speed improvement over
the current analytic-corrected X+Y implementation at the tested five-prime
`N=10^12` target. It also does not discharge the full
Frederickson-Johnson or soft-heap comparison obligations, because the selected
published algorithm here is Mirzaian-Arjomandi value selection plus our exact
band reconstruction.

Clean first-k comparator artifact:

```text
results/benchmarks/ma_full_unrank_first_k_1e12/
```

Command:

```bash
python3 scripts/run_ma_full_unrank_suite.py \
  --out-dir results/benchmarks/ma_full_unrank_first_k_1e12 \
  --prime-set 2,3,5,7,11 \
  --prime-set 2,3,5,7,11,13 \
  --prime-set 2,3,5,7,11,13,17,19 \
  --rank-radius 250 \
  --max-candidates 200000 \
  --ma-max-n 30000000 \
  --ma-max-middle 200000000
```

Observed result at commit `7daf269404c1672a69da6f6d2edb2ff53eb648e3`:

- 3/3 rows completed.
- 3/3 returned the same exponent vector as analytic-corrected X+Y.
- 3/3 matching vectors were independently interval-rank certified.
- The Mirzaian-Arjomandi full-unrank path won 0/3 wall-time comparisons.
- MA/current wall-time ratios were 3.5898 for k=5, 4.0593 for k=6, and
  3.7086 for k=8, with mean 3.7859.

This is the current strongest implemented published-selector checkpoint for
"where are we versus a best-known-style comparator?" It supports a narrow claim
against this MA sorted-matrix selector wrapper only.

## Best-Known Comparator Gate

The comparator gate is the current "no hiding" status artifact for the serious
baseline obligations:

```bash
python3 scripts/run_best_known_comparator_gate.py \
  --out-dir results/benchmarks/best_known_comparator_gate \
  --normaliz-scale-digits 1 \
  --normaliz-timeout-seconds 5 \
  --normaliz-max-cases 3
```

Clean artifact:

```text
results/benchmarks/best_known_comparator_gate/
```

Observed result in the gate artifact:

- Mirzaian-Arjomandi full-unrank gate passed.
- Output-sensitive `X+Y`/LOH probe gate executed at N=10^6.
- Full Frederickson-Johnson gate remained open/not implemented.
- Soft-heap `X+Y` gate advanced to
  `semantic_probe_validated_not_selector_integrated`.
- Barvinok-style external count gate found installed tools, exported rational
  simplex inputs, passed PyNormaliz toy simplex counts, and timed out on all
  certified-target rationalized simplexes under the configured cap.

Numerical summary:

- MA first-k comparison: 3/3 completed and certified, 0/3 MA wall-time wins,
  mean MA/current wall ratio 3.8757.
- Output-sensitive probe: 3/3 completed, range-pruning won 2/3 against linear
  saddleback count at N=10^6 with mean block/linear ratio 0.8825; MA value
  selection won 0/3 with mean MA/linear ratio 2.1749. Mat-Select2 heap-primitive
  matched the selected log in 3/3 rows, won 0/3 timing comparisons, and had mean
  Mat-Select2/linear ratio 9.6202.
- Soft-sequence-heap probe: 3/3 validation rows passed at
  `n=2048, epsilon=0.25`. The timing probe at `n=20000` was negative for speed:
  the vector-backed soft-sequence heap was 304.6520x slower than a binary heap.
- Barvinok/Normaliz statuses: `barvinok_count` returned -11 on toy 1D and 2D
  ISL inputs. PyNormaliz passed 2/2 known rational-vertex toy simplex counts,
  but timed out on 3/3 certified-target rationalized simplexes at the 5-second
  per-count cap.

The gate is intentionally not a broad best-known claim. It is the current
auditable checklist for what has been compared, what has failed, and what still
requires implementation.

## Iterative Corrected High-k Sums-Only Suite

The one-step residual correction used in earlier analytic-band experiments can
still leave a rank bias for higher dimensions. The high-k suite therefore uses
multiple exact residual corrections:

1. start from the analytic lattice estimate;
2. count exactly at the current estimate using the sums-only MITM counter;
3. apply a Newton-style rank residual correction using the analytic derivative;
4. repeat for the configured number of refinement steps;
5. open a narrow final band, reconstruct exponent vectors by parent traces, and
   exact-sort the band by multiprecision integer value;
6. independently audit the corrected exponent vector.

Use:

```bash
python3 scripts/run_iterative_corrected_highk_suite.py \
  --out-dir results/benchmarks/iterative_corrected_highk_suite_1e12
```

Default outputs are written to:

```text
results/local/iterative_corrected_highk_<timestamp>/
```

The clean default artifact compares adaptive sums-only MITM against iterative
corrected sums-only MITM at `N=10^12` using one shared policy:
`rank_radius=25`, `max_candidates=200000`, and `refine_steps=4`.

The default cases are the first eight, first ten, and first twelve primes.

Clean artifact:

```text
results/benchmarks/iterative_corrected_highk_suite_1e12/
```

Observed result at commit `e0825a6f7ed1cbadaf823897f2b918037a6c7be9`:

- 3/3 rows completed.
- 3/3 corrected rows matched adaptive sums-only exponent vectors.
- 3/3 corrected rows were independently interval-rank certified.
- Corrected mode won 3/3 wall-time comparisons.
- Adaptive/corrected wall-time ratios were 1.8605 for `k=8`, 1.7140 for
  `k=10`, and 2.3339 for `k=12`, with mean 1.9695.
- Final exact bands had 52 candidates for `k=8`, 50 candidates for `k=10`,
  and 38 candidates for `k=12`.

This is a current-code speed claim against the implemented adaptive sums-only
MITM unranker. It is not a broad best-known claim: the comparator set still
does not include full soft-heap X+Y, full Frederickson-Johnson ranking/selection,
or a Barvinok-style fixed-dimensional lattice-count implementation.
