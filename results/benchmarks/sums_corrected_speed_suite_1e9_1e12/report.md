# Residual-Corrected Sums-Only MITM Speed Suite

- Timestamp: `2026-06-10T06:44:53.034393+00:00`
- Git commit: `be4e24a8f4adf0947e568dcc57f4797b0b35ca60`
- Git dirty: `False`
- Seed: `20260609`
- Cases: `10`
- Completed cases: `10`
- Certified matching cases: `10`
- Corrected wall-time wins: `9`
- Corrected reported-time wins: `10`
- Corrected memory wins: `10`
- Mean wall ratio adaptive/corrected: `1.333458489201051`
- Mean reported ratio adaptive/corrected: `1.6141397780332163`
- Mean RSS ratio adaptive/corrected: `1.1398223445454323`
- Max corrected band count: `19970`

This suite compares two modes of the same sums-only MITM kernel. The
corrected mode uses the analytic residual correction to center a small
rank band and then sorts only that band. The adaptive mode is the prior
interpolation/bisection sums-only unrank path. Corrected outputs are
checked with the independent interval-rank auditor, and adaptive outputs
must match the certified corrected exponent vector.

| label | k | N | corr wall | adaptive wall | wall ratio | corr band | adaptive band | calls | cert | exps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| `k5_N1000000000_case0_P3_5_7_17_19` | 5 | 1000000000 | 0.609537 | 0.047656 | 0.078183 | 525 | 557 | 10 | True | `[55, 36, 39, 13, 29]` |
| `k5_N1000000000_case1_P3_7_11_13_19` | 5 | 1000000000 | 0.033992 | 0.045286 | 1.332249 | 504 | 466 | 11 | True | `[38, 77, 16, 23, 16]` |
| `k6_N1000000000_case0_P2_3_7_11_13_17` | 6 | 1000000000 | 0.035306 | 0.054362 | 1.539736 | 2967 | 5076 | 9 | True | `[14, 33, 29, 5, 10, 6]` |
| `k6_N1000000000_case1_P2_3_5_11_13_19` | 6 | 1000000000 | 0.051270 | 0.052683 | 1.027558 | 2995 | 4693 | 9 | True | `[106, 6, 4, 0, 3, 20]` |
| `k8_N1000000000_case0_P2_3_5_7_11_13_17_19` | 8 | 1000000000 | 0.062753 | 0.100129 | 1.595604 | 19970 | 17594 | 12 | True | `[5, 2, 28, 5, 5, 1, 2, 1]` |
| `k5_N1000000000000_case0_P2_3_5_11_17` | 5 | 1000000000000 | 0.714788 | 1.137313 | 1.591120 | 513 | 2030 | 12 | True | `[400, 112, 101, 93, 74]` |
| `k5_N1000000000000_case1_P3_5_11_13_17` | 5 | 1000000000000 | 0.886818 | 1.476247 | 1.664656 | 495 | 3428 | 12 | True | `[124, 147, 147, 57, 149]` |
| `k6_N1000000000000_case0_P3_5_7_13_17_19` | 6 | 1000000000000 | 0.729222 | 1.090382 | 1.495268 | 2987 | 5900 | 12 | True | `[83, 125, 58, 14, 23, 34]` |
| `k6_N1000000000000_case1_P2_5_7_11_13_19` | 6 | 1000000000000 | 0.680456 | 1.049856 | 1.542872 | 3059 | 5914 | 12 | True | `[161, 12, 56, 88, 21, 14]` |
| `k8_N1000000000000_case0_P2_3_5_7_11_13_17_19` | 8 | 1000000000000 | 1.354856 | 1.988034 | 1.467339 | 19854 | 1112 | 13 | True | `[75, 28, 9, 16, 3, 22, 5, 1]` |
