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

Observed result at commit `13cba09636b22e139046442b9faa4ea0f5330377` on the
recorded macOS/x86_64 Apple-clang machine:

- 6/6 cases completed.
- 6/6 cases returned matching exponent vectors.
- 6/6 layer outputs and 6/6 full-X+Y outputs were independently certified.
- Layer-compressed full unranking won 6/6 wall-time comparisons.
- Full-X+Y/layer wall-time ratio ranged from 1.4703 to 3.1276, with mean
  2.3309.

This supports the safe claim above for this benchmark suite and recorded
machine/compiler. It is still not, by itself, a broad "best known" claim.
