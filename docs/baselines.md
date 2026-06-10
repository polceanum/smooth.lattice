# Baselines and comparator families

This project targets fixed-prime random access, not general smooth-number counting.

## Included baselines

- `benchmarks/dp_pointer_baseline.cpp`: standard DP/pointer method used for Hamming/super-ugly-number generation. It is a strong sequential-generation baseline but has O(Nk) time and O(N) memory if it must return rank N.
- `benchmarks/smooth_xplusy_baseline.cpp`: practical adaptive Cartesian-sum value-selection baseline.
- `benchmarks/smooth_xplusy_full_unrank.cpp`: practical materialized Cartesian-sum full-unrank baseline. It stores exponent packs on both MITM sides, narrows to a log-value band, exact-sorts the candidate band by multiprecision integer value, and returns an exponent vector. Its modes include adaptive X+Y, residual-corrected analytic X+Y, and a Mirzaian-Arjomandi sorted-matrix selector wrapped with the same exact reconstruction step.
- `benchmarks/smooth_xplusy_fj_loh_workbench.cpp`: exploratory sorted-matrix/range-pruning, Mirzaian-Arjomandi value-selection, and LOH-style probes. Use `scripts/run_sorted_matrix_workbench.py` for reproducible artifacts; this is not a faithful Frederickson-Johnson or soft-heap implementation.

## Current priority comparison

The next narrow target is the first-five-primes case:

```text
P = (2,3,5,7,11)
N = 10^12
```

Use `python3 scripts/run_xplusy_vs_layer5.py` to compare the practical
adaptive `X+Y` value-selection baseline against the layer-compressed full
unrank solver. The harness records command lines, compiler metadata, wall time,
peak RSS when available, raw stdout/stderr, and interval-audit status.

The `X+Y` row is selection-only, not full exponent reconstruction. A win by the
layer-compressed solver is therefore a conservative result; a loss should be
reported as such.

The stricter follow-up target is the full materialized `X+Y` unrank baseline:

```bash
python3 scripts/run_full_xplusy_suite.py
```

This compares two full exponent-vector unrankers across the same six five-prime
suite. The harness audits both returned vectors independently and requires the
vectors to match before counting a case as certified comparison evidence.

Current clean artifact:

```text
results/benchmarks/xplusy_vs_layer5_1e12/
results/benchmarks/five_prime_suite_1e12/
results/benchmarks/full_xplusy_suite_1e12/
results/benchmarks/sorted_matrix_workbench_1e12/
results/benchmarks/ma_full_unrank_suite_1e12/
```

On the recorded macOS/x86_64 Apple-clang run at commit
`161445507617a9435f9baadf4e70a3679d9e8d9a`, layer-compressed full unrank took
1.194730s wall time versus 2.029418s for adaptive materialized `X+Y` value
selection, with the returned exponent vector independently certified.

On the broader six-case suite over all five-prime subsets of
`{2,3,5,7,11,13}` at `N=10^12`, layer-compressed full unrank won all six
wall-time comparisons against the same adaptive materialized `X+Y` value
selection baseline, and all six returned exponent vectors were independently
certified.

On the stricter six-case suite against full materialized `X+Y` unranking,
recorded at commit `5314fd8a56fb7ca760046076e1f7f168fd48386a`,
layer-compressed full unrank won all six wall-time comparisons. Both methods
returned the same exponent vector in every case, and both outputs were
independently certified. The full-X+Y/layer wall-time ratio ranged from 1.9672
to 2.7248, with mean 2.3273. In this artifact the layer solver used the
analytic asymptotic bracket in all six cases; this is an optimization seed, not
a certificate.

Use `python3 scripts/run_sorted_matrix_workbench.py` to record exploratory
sorted-matrix/range-pruning, Mirzaian-Arjomandi value-selection, and LOH probes.
These rows are useful negative or diagnostic evidence, but they do not discharge
the "full Frederickson-Johnson" or "full soft-heap X+Y" comparison obligations.

In the clean six-case `N=10^12` workbench artifact at commit
`17a7e40ae790ae9017ab8274204eaaf107add212`, the range-pruned block counter
beat the ordinary linear saddleback count in only 1/6 cases, and only narrowly.
The mean block/linear internal time ratio was 1.3467, so this probe does not
improve the practical X+Y comparator. The Mirzaian-Arjomandi value selector
passed exhaustive small validation on 5185/5185 cases, matched the adaptive
selected log in 6/6 large cases, but won 0/6 timing comparisons, with mean
MA/linear internal time ratio 9.4042. The LOH row is a capped
`N_probe=10^6` top-k probe, not a full-rank `N=10^12` random-access comparator.

Use `python3 scripts/run_ma_full_unrank_suite.py` for the full
Mirzaian-Arjomandi comparator. This is no longer merely a value-selection
probe: it selects the X+Y log value, reconstructs exponent vectors, exact-sorts
the candidate band, and audits matching vectors. In the six-case `N=10^12`
artifact, the MA full-unrank path matched the analytic-corrected X+Y output and
was independently certified in 6/6 cases, but won 0/6 wall-time comparisons.
The mean MA/corrected wall-time ratio was 5.7000, so this comparator is useful
negative evidence rather than a speed path.

## Not yet fully implemented

- Full Frederickson-Johnson sorted matrix selection.
- Full soft-heap X+Y selection implementations.
- Barvinok-style fixed-dimensional lattice-point counting.
- Bernstein-style general y-smooth counting/listing algorithms.

Claims should be limited accordingly.
