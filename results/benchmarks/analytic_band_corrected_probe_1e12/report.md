# Residual-Corrected Analytic Boundary Band Probe

- Timestamp: `2026-06-09T21:28:43.799727+00:00`
- Git commit: `f965d9cb677d26aa83b7b35c8412a90f856c8fc2`
- Git dirty: `False`
- Cases: `8`
- Completed cases: `8`
- Target inside band: `8`
- Enumerated cases: `8`
- Recovered expected vectors: `8`
- Max band count: `19853`
- Mean band count: `2756.625`

The band center starts at the analytic solution for rank N. The harness
then performs one exact layer count at that center, shifts the center by
`(N - count(center)) / analytic_derivative(center)`, and checks a smaller
rank-radius band around the corrected center. Endpoint counts and band
enumeration still use the floating-log layer counter, so this remains an
experimental oracle probe rather than a correctness certificate.

| label | k | radius | center error | band count | inside | enumerated | recovered | wall s |
|---|---:|---:|---:|---:|---|---|---|---:|
| `k5_P2_3_5_7_11` | 5 | 100 | -68 | 201 | True | True | True | 1.262695 |
| `k5_P2_3_5_7_13` | 5 | 100 | 60 | 204 | True | True | True | 0.394480 |
| `k5_P2_3_5_11_13` | 5 | 100 | -14 | 204 | True | True | True | 0.411428 |
| `k5_P2_3_7_11_13` | 5 | 100 | 51 | 193 | True | True | True | 0.420623 |
| `k5_P2_5_7_11_13` | 5 | 100 | -158 | 197 | True | True | True | 0.336689 |
| `k5_P3_5_7_11_13` | 5 | 100 | -89 | 204 | True | True | True | 0.419742 |
| `k6_first6` | 6 | 500 | 452 | 997 | True | True | True | 1.989259 |
| `k8_first8` | 8 | 10000 | 1979534 | 19853 | True | True | True | 3.490245 |
