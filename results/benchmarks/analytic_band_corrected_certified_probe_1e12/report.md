# Residual-Corrected Analytic Boundary Band Probe

- Timestamp: `2026-06-09T22:05:14.037528+00:00`
- Git commit: `11687a9ef1b183ff88ce3fea53bca4c281315228`
- Git dirty: `False`
- Cases: `8`
- Completed cases: `8`
- Target inside band: `8`
- Enumerated cases: `8`
- Recovered expected vectors: `8`
- Independently rank-certified vectors: `8`
- Max band count: `19853`
- Mean band count: `2756.625`

The band center starts at the analytic solution for rank N. The harness
then performs one exact layer count at that center, shifts the center by
`(N - count(center)) / analytic_derivative(center)`, and checks a smaller
rank-radius band around the corrected center. Endpoint counts and band
enumeration still use the floating-log layer counter. Each recovered vector
is then passed to the independent interval-log rank auditor; only rows with
`rank_certified=True` should be used as certified correctness evidence.

| label | k | radius | center error | band count | inside | recovered | certified | wall s | audit s |
|---|---:|---:|---:|---:|---|---|---|---:|---:|
| `k5_P2_3_5_7_11` | 5 | 100 | -68 | 201 | True | True | True | 1.173644 | 6.863999 |
| `k5_P2_3_5_7_13` | 5 | 100 | 60 | 204 | True | True | True | 0.386204 | 5.766177 |
| `k5_P2_3_5_11_13` | 5 | 100 | -14 | 204 | True | True | True | 0.412296 | 4.766443 |
| `k5_P2_3_7_11_13` | 5 | 100 | 51 | 193 | True | True | True | 0.441604 | 4.507973 |
| `k5_P2_5_7_11_13` | 5 | 100 | -158 | 197 | True | True | True | 0.292441 | 6.114585 |
| `k5_P3_5_7_11_13` | 5 | 100 | -89 | 204 | True | True | True | 0.445531 | 7.094919 |
| `k6_first6` | 6 | 500 | 452 | 997 | True | True | True | 2.081170 | 8.849259 |
| `k8_first8` | 8 | 10000 | 1979534 | 19853 | True | True | True | 3.555122 | 18.198369 |
