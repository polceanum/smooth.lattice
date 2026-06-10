# Residual-Corrected Full X+Y Speed Suite

- Timestamp: `2026-06-10T05:31:05.528455+00:00`
- Git commit: `4a581fe2fb495997442c4e233fc01ad172a2d4df`
- Git dirty: `False`
- Seed: `20260609`
- Cases: `10`
- Completed cases: `10`
- Certified matching cases: `10`
- Corrected wall-time wins: `9`
- Corrected reported-time wins: `10`
- Corrected memory wins: `10`
- Mean wall ratio adaptive/corrected: `1.3323180717903615`
- Mean reported ratio adaptive/corrected: `1.815214100132841`
- Mean RSS ratio adaptive/corrected: `1.343461310292311`
- Max corrected band count: `19970`

This suite compares two modes of the same full materialized X+Y unranker.
The corrected mode uses the analytic residual correction to choose a small
boundary band, then exact-sorts that band. The adaptive mode is the prior
interpolation/bisection full-X+Y unrank path. Corrected outputs are audited
with the independent interval-rank auditor; adaptive outputs must match the
certified corrected exponent vector.

| label | k | N | corr wall | adaptive wall | wall ratio | corr band | adaptive band | calls | cert | exps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| `k5_N1000000000_case0_P3_5_7_17_19` | 5 | 1000000000 | 0.818732 | 0.141819 | 0.173217 | 525 | 31255 | 9 | True | `[55, 36, 39, 13, 29]` |
| `k5_N1000000000_case1_P3_7_11_13_19` | 5 | 1000000000 | 0.062649 | 0.136345 | 2.176335 | 504 | 31134 | 9 | True | `[38, 77, 16, 23, 16]` |
| `k6_N1000000000_case0_P2_3_7_11_13_17` | 6 | 1000000000 | 0.050352 | 0.069904 | 1.388289 | 2967 | 5076 | 9 | True | `[14, 33, 29, 5, 10, 6]` |
| `k6_N1000000000_case1_P2_3_5_11_13_19` | 6 | 1000000000 | 0.036740 | 0.049699 | 1.352732 | 2995 | 4693 | 9 | True | `[106, 6, 4, 0, 3, 20]` |
| `k8_N1000000000_case0_P2_3_5_7_11_13_17_19` | 8 | 1000000000 | 0.101818 | 0.108875 | 1.069307 | 19970 | 17594 | 12 | True | `[5, 2, 28, 5, 5, 1, 2, 1]` |
| `k5_N1000000000000_case0_P2_3_5_11_17` | 5 | 1000000000000 | 1.163798 | 1.682264 | 1.445494 | 513 | 2030 | 12 | True | `[400, 112, 101, 93, 74]` |
| `k5_N1000000000000_case1_P3_5_11_13_17` | 5 | 1000000000000 | 1.502254 | 2.264301 | 1.507269 | 495 | 3428 | 12 | True | `[124, 147, 147, 57, 149]` |
| `k6_N1000000000000_case0_P3_5_7_13_17_19` | 6 | 1000000000000 | 0.909675 | 1.087378 | 1.195348 | 2987 | 5900 | 12 | True | `[83, 125, 58, 14, 23, 34]` |
| `k6_N1000000000000_case1_P2_5_7_11_13_19` | 6 | 1000000000000 | 0.729475 | 1.128937 | 1.547602 | 3059 | 5914 | 12 | True | `[161, 12, 56, 88, 21, 14]` |
| `k8_N1000000000000_case0_P2_3_5_7_11_13_17_19` | 8 | 1000000000000 | 1.626452 | 2.386960 | 1.467588 | 19854 | 1112 | 13 | True | `[75, 28, 9, 16, 3, 22, 5, 1]` |
