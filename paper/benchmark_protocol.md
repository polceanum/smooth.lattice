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
- The LOH row is an output-style top-k probe. When `N_probe != N`, it is not a
  random-access unranking comparator for the full target rank.
- These results should be reported as diagnostic or negative evidence unless a
  future patch implements the full published algorithm being invoked.

Clean workbench artifact:

```text
results/benchmarks/sorted_matrix_workbench_1e12/
```

Observed result at commit `44175948b9b9bf3010f5fd388e2d926180fe5f45` on the
recorded macOS/x86_64 Apple-clang machine:

- 6/6 cases completed.
- The range-pruned block counter beat the ordinary linear saddleback count in
  2/6 cases.
- Those two wins were narrow: block/linear ratios 0.9634 and 0.9907.
- The mean block/linear internal time ratio was 1.2765, so the probe was slower
  on average.
- The Mirzaian-Arjomandi value selector matched the adaptive selected log in
  6/6 cases, but won 0/6 timing comparisons. Its mean MA/linear internal time
  ratio was 9.0509.
- The LOH row used `N_probe=10^6`, not the full target rank `N=10^12`.

This is a negative/mixed result. It should be used to document that the current
sorted-matrix probes do not strengthen the X+Y baseline at the headline
five-prime target scale.
