# Certified Corrected-Oracle Random Suite

- Timestamp: `2026-06-09T23:14:13.606086+00:00`
- Git commit: `e53d69b364d568a49e42784892289a15aeb58742`
- Git dirty: `False`
- Seed: `20260609`
- Cases: `10`
- Completed cases: `10`
- Certified corrected-oracle cases: `10`
- All-method same-exponent cases: `10`
- Corrected wall-time wins vs sums-only: `3`
- Corrected wall-time wins vs full X+Y: `4`
- Corrected memory wins vs full X+Y: `10`
- Mean wall ratio sums/corrected: `1.2991960920797336`
- Mean wall ratio full X+Y/corrected: `2.000232086650819`
- Mean RSS ratio full X+Y/corrected: `7.895893008836071`
- Max corrected band count: `19970`

Cases are deterministic pseudo-random subsets of the supported audit prime
universe. The corrected analytic oracle is independently interval-audited;
the sums-only and full materialized X+Y baselines must return the same
exponent vector to count as matching the certified result.

| label | k | N | band | cert | sums same | X+Y same | corr wall | sums wall | X+Y wall | X+Y/corr RSS | exps |
|---|---:|---:|---:|---|---|---|---:|---:|---:|---:|---|
| `k5_N1000000000_case0_P3_5_7_17_19` | 5 | 1000000000 | 525 | True | True | True | 0.806470 | 0.705755 | 0.742781 | 8.126 | `[55, 36, 39, 13, 29]` |
| `k5_N1000000000_case1_P3_7_11_13_19` | 5 | 1000000000 | 504 | True | True | True | 0.036360 | 0.051490 | 0.154154 | 10.045 | `[38, 77, 16, 23, 16]` |
| `k6_N1000000000_case0_P2_3_7_11_13_17` | 6 | 1000000000 | 2967 | True | True | True | 0.051786 | 0.044973 | 0.056859 | 2.152 | `[14, 33, 29, 5, 10, 6]` |
| `k6_N1000000000_case1_P2_3_5_11_13_19` | 6 | 1000000000 | 2995 | True | True | True | 0.061327 | 0.052678 | 0.056498 | 2.254 | `[106, 6, 4, 0, 3, 20]` |
| `k8_N1000000000_case0_P2_3_5_7_11_13_17_19` | 8 | 1000000000 | 19970 | True | True | True | 0.121037 | 0.081264 | 0.113011 | 2.075 | `[5, 2, 28, 5, 5, 1, 2, 1]` |
| `k5_N1000000000000_case0_P2_3_5_11_17` | 5 | 1000000000000 | 513 | True | True | True | 0.364100 | 1.187604 | 1.721467 | 21.533 | `[400, 112, 101, 93, 74]` |
| `k5_N1000000000000_case1_P3_5_11_13_17` | 5 | 1000000000000 | 495 | True | True | True | 0.414671 | 1.467992 | 2.251960 | 25.868 | `[124, 147, 147, 57, 149]` |
| `k6_N1000000000000_case0_P3_5_7_13_17_19` | 6 | 1000000000000 | 2987 | True | True | True | 2.191503 | 1.097639 | 1.130048 | 1.778 | `[83, 125, 58, 14, 23, 34]` |
| `k6_N1000000000000_case1_P2_5_7_11_13_19` | 6 | 1000000000000 | 3059 | True | True | True | 2.253830 | 1.086017 | 1.187812 | 2.093 | `[161, 12, 56, 88, 21, 14]` |
| `k8_N1000000000000_case0_P2_3_5_7_11_13_17_19` | 8 | 1000000000000 | 19853 | True | True | True | 3.581938 | 1.852949 | 2.462018 | 3.036 | `[75, 28, 9, 16, 3, 22, 5, 1]` |
