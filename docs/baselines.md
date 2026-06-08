# Baselines and comparator families

This project targets fixed-prime random access, not general smooth-number counting.

## Included baselines

- `benchmarks/dp_pointer_baseline.cpp`: standard DP/pointer method used for Hamming/super-ugly-number generation. It is a strong sequential-generation baseline but has O(Nk) time and O(N) memory if it must return rank N.
- `benchmarks/smooth_xplusy_baseline.cpp`: practical adaptive Cartesian-sum value-selection baseline.
- `benchmarks/smooth_xplusy_full_unrank.cpp`: practical materialized Cartesian-sum full-unrank baseline. It stores exponent packs on both MITM sides, narrows to a log-value band, exact-sorts the candidate band by multiprecision integer value, and returns an exponent vector.
- `benchmarks/smooth_xplusy_fj_loh_workbench.cpp`: exploratory sorted-matrix/range-pruning and LOH-style probes.

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
recorded at commit `13cba09636b22e139046442b9faa4ea0f5330377`,
layer-compressed full unrank won all six wall-time comparisons. Both methods
returned the same exponent vector in every case, and both outputs were
independently certified. The full-X+Y/layer wall-time ratio ranged from 1.4703
to 3.1276, with mean 2.3309.

## Not yet fully implemented

- Full Frederickson-Johnson sorted matrix selection.
- Full soft-heap X+Y selection implementations.
- Barvinok-style fixed-dimensional lattice-point counting.
- Bernstein-style general y-smooth counting/listing algorithms.

Claims should be limited accordingly.
