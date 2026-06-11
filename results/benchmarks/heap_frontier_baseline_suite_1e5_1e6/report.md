# Heap Frontier Baseline Suite

- Timestamp: `2026-06-11T15:11:45.035924+00:00`
- Git commit: `b26e3cdf865b15754d8b189dafe4e7cc36e19f80`
- Git dirty: `False`
- Cases: `8`
- Certified cases: `8`
- Solver-agreement cases: `8`
- Solver-disagreement cases: `0`
- Wall-time winners: `{'current_beatty3': 2, 'current_sums_adaptive': 5, 'dp_pointer': 1}`
- Mean heap/current wall ratio: `4.565355705191498`
- Mean heap/DP wall ratio: `1.964112652773901`

The heap row is a canonical frontier generator that returns the full
exponent vector. A row is certified only when the heap, DP baseline when
enabled, and current solver agree on the exponent vector and the
independent interval auditor proves `count_le(candidate)=N`.

| label | k | N | fastest | heap s | DP s | current s | heap/current | cert | exps |
|---|---:|---:|---|---:|---:|---:|---:|---|---|
| `k3_N1e5` | 3 | 100000 | `current_beatty3` | 0.684588 | 0.516886 | 0.448630 | 1.5260 | True | `96,1,13` |
| `k5_N1e5` | 5 | 100000 | `dp_pointer` | 0.024340 | 0.013995 | 0.517746 | 0.0470 | True | `4,10,3,5,2` |
| `k6_N1e5` | 6 | 100000 | `current_sums_adaptive` | 0.024840 | 0.019108 | 0.012077 | 2.0568 | True | `0,1,7,5,2,0` |
| `k8_N1e5` | 8 | 100000 | `current_sums_adaptive` | 0.025308 | 0.018966 | 0.015881 | 1.5936 | True | `15,1,1,0,1,0,1,1` |
| `k3_N1e6` | 3 | 1000000 | `current_beatty3` | 0.125638 | 0.047210 | 0.015937 | 7.8836 | True | `55,47,64` |
| `k5_N1e6` | 5 | 1000000 | `current_sums_adaptive` | 0.155668 | 0.074191 | 0.026521 | 5.8696 | True | `3,12,11,5,5` |
| `k6_N1e6` | 6 | 1000000 | `current_sums_adaptive` | 0.245933 | 0.089303 | 0.033362 | 7.3717 | True | `9,22,1,0,3,1` |
| `k8_N1e6` | 8 | 1000000 | `current_sums_adaptive` | 0.254774 | 0.101844 | 0.025040 | 10.1745 | True | `10,7,1,0,0,0,1,4` |

## Claim Status

This suite discharges the repository's heap/frontier baseline obligation
for the tested small-to-medium ranks. It is not a broad best-known
claim, and it does not replace the still-open soft-heap, full
Frederickson-Johnson, or Barvinok-style comparator obligations.
