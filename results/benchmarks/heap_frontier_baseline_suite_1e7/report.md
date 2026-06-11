# Heap Frontier Baseline Suite

- Timestamp: `2026-06-11T15:46:09.209119+00:00`
- Git commit: `6526d68a0adcd45b0041a1dcfe8da1c6a9fe1f56`
- Git dirty: `False`
- Cases: `4`
- Certified cases: `4`
- Solver-agreement cases: `4`
- Solver-disagreement cases: `0`
- Wall-time winners: `{'current_beatty3': 1, 'current_sums_adaptive': 3}`
- Mean heap/current wall ratio: `76.53146049779198`
- Mean heap/DP wall ratio: `2.953644796217955`

The heap row is a canonical frontier generator that returns the full
exponent vector. A row is certified only when the heap, DP baseline when
enabled, and current solver agree on the exponent vector and the
independent interval auditor proves `count_le(candidate)=N`.

| label | k | N | fastest | heap s | DP s | current s | heap/current | cert | exps |
|---|---:|---:|---|---:|---:|---:|---:|---|---|
| `k3_N1e7` | 3 | 10000000 | `current_beatty3` | 2.186698 | 0.798331 | 0.557295 | 3.9238 | True | `80,92,162` |
| `k5_N1e7` | 5 | 10000000 | `current_sums_adaptive` | 1.877599 | 0.644027 | 0.471126 | 3.9853 | True | `9,26,15,8,6` |
| `k6_N1e7` | 6 | 10000000 | `current_sums_adaptive` | 2.087263 | 0.626187 | 0.013827 | 150.9563 | True | `35,3,6,7,0,5` |
| `k8_N1e7` | 8 | 10000000 | `current_sums_adaptive` | 2.344829 | 0.829500 | 0.015923 | 147.2605 | True | `8,2,3,3,7,1,0,2` |

## Claim Status

This suite discharges the repository's heap/frontier baseline obligation
for the tested small-to-medium ranks. It is not a broad best-known
claim, and it does not replace the still-open soft-heap, full
Frederickson-Johnson, or Barvinok-style comparator obligations.
