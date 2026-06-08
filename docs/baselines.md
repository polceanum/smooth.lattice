# Baselines and comparator families

This project targets fixed-prime random access, not general smooth-number counting.

## Included baselines

- `benchmarks/dp_pointer_baseline.cpp`: standard DP/pointer method used for Hamming/super-ugly-number generation. It is a strong sequential-generation baseline but has O(Nk) time and O(N) memory if it must return rank N.
- `benchmarks/smooth_xplusy_baseline.cpp`: practical adaptive Cartesian-sum value-selection baseline.
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

## Not yet fully implemented

- Full Frederickson-Johnson sorted matrix selection.
- Full soft-heap X+Y selection implementations.
- Barvinok-style fixed-dimensional lattice-point counting.
- Bernstein-style general y-smooth counting/listing algorithms.

Claims should be limited accordingly.
