# X+Y vs Layer-Compressed Benchmark

- Timestamp: `2026-06-08T20:28:49.366807+00:00`
- Git commit: `161445507617a9435f9baadf4e70a3679d9e8d9a`
- Git dirty: `False`
- Primes: `2,3,5,7,11`
- N: `1000000000000`
- CXX: `g++`
- CXXFLAGS: `-O3 -std=c++17`

The `xplusy_adaptive` row measures adaptive Cartesian-sum value selection only.
The `layer_compressed` row measures full exponent-vector unranking.

| method | return | wall s | max RSS KB | reported s | exps | certified |
|---|---:|---:|---:|---:|---|---|
| build_core | 0 | 32.458152 | 171740 |  | `` |  |
| build_xplusy_baseline | 0 | 5.046990 | 151560 |  | `` |  |
| xplusy_adaptive | 0 | 2.029418 | 153444 | 0.696378 | `` |  |
| layer_compressed | 0 | 1.194730 | 23264 | 0.703781 | `[1052, 26, 33, 53, 4]` |  |
| interval_audit | 0 | 6.127383 | 997012 |  | `[1052, 26, 33, 53, 4]` | True |

## Comparison

- `xplusy_wall_over_layer_wall`: `1.698641670208074`
- `xplusy_reported_total_over_layer_reported`: `1.7977353750669598`
- `same_rank`: `True`
- `same_primes`: `True`
- `rank_certified`: `True`
- `certified_count_le`: `1000000000000`
