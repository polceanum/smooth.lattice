# Analytic Boundary Band Probe

- Timestamp: `2026-06-09T19:13:32.981117+00:00`
- Git commit: `4bd49bec50aab80c87e599438f325f55814e6222`
- Git dirty: `False`
- Cases: `8`
- Completed cases: `8`
- Target inside band: `8`
- Enumerated cases: `7`
- Recovered expected vectors: `7`
- Max band count: `5000120`
- Mean band count: `627014.625`

The band center is the analytic solution for rank N. The half-width is
`rank_radius / analytic_derivative`; endpoints are checked with the current
layer counter. Cases are enumerated only when the endpoint band count is at
or below `max_candidates`. This is an experimental correction-band probe,
not a certified analytic oracle.

| label | k | radius | band count | inside | enumerated | recovered | wall s | cands |
|---|---:|---:|---:|---|---|---|---:|---:|
| `k5_P2_3_5_7_11` | 5 | 1000 | 1989 | True | True | True | 0.623578 | 1989 |
| `k5_P2_3_5_7_13` | 5 | 1000 | 1980 | True | True | True | 0.195691 | 1980 |
| `k5_P2_3_5_11_13` | 5 | 1000 | 2001 | True | True | True | 0.222866 | 2001 |
| `k5_P2_3_7_11_13` | 5 | 1000 | 2032 | True | True | True | 0.210268 | 2032 |
| `k5_P2_5_7_11_13` | 5 | 1000 | 1988 | True | True | True | 0.176107 | 1988 |
| `k5_P3_5_7_11_13` | 5 | 1000 | 2002 | True | True | True | 0.197228 | 2002 |
| `k6_first6` | 6 | 2000 | 4005 | True | True | True | 0.966114 | 4005 |
| `k8_first8` | 8 | 2500000 | 5000120 | True | False | False | 1.158039 | 0 |
