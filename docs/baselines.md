# Baselines and comparator families

This project targets fixed-prime random access, not general smooth-number counting.

## Included baselines

- `benchmarks/dp_pointer_baseline.cpp`: standard DP/pointer method used for Hamming/super-ugly-number generation. It is a strong sequential-generation baseline but has O(Nk) time and O(N) memory if it must return rank N.
- `benchmarks/smooth_xplusy_baseline.cpp`: practical adaptive Cartesian-sum value-selection baseline.
- `benchmarks/smooth_xplusy_fj_loh_workbench.cpp`: exploratory sorted-matrix/range-pruning and LOH-style probes.

## Not yet fully implemented

- Full Frederickson-Johnson sorted matrix selection.
- Full soft-heap X+Y selection implementations.
- Barvinok-style fixed-dimensional lattice-point counting.
- Bernstein-style general y-smooth counting/listing algorithms.

Claims should be limited accordingly.
