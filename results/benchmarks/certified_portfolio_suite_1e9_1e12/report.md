# Certified Portfolio Suite

- Timestamp: `2026-06-10T20:11:03.859826+00:00`
- Git commit: `b8357ac6e09c27b781892b9a304d72c184d2abae`
- Git dirty: `False`
- Cases: `14`
- Certified cases: `14`
- Solver-agreement cases: `14`
- Audit-blocked agreement cases: `0`
- Solver-disagreement cases: `0`
- Wall-time winners: `{'beatty3': 3, 'layer_corrected': 4, 'xplusy_adaptive': 2, 'xplusy_corrected': 5}`
- Reported-time winners: `{'beatty3': 4, 'layer_corrected': 4, 'xplusy_corrected': 6}`

This suite is a portfolio benchmark over implemented full-unrank methods.
A row is certified only when successful methods agree on the exponent vector
and the independent interval auditor proves `count_le(candidate)=N`.

| label | k | N | fastest wall | wall sec | fastest reported | reported sec | cert | audit blocked | exps |
|---|---:|---:|---|---:|---|---:|---|---|---|
| `k3_N1000000000_case0_P2_3_5` | 3 | 1000000000 | `xplusy_adaptive` | 0.172574 | `beatty3` | 0.044463 | True | False | `1334,335,404` |
| `k3_N1000000000_case1_P2_3_7` | 3 | 1000000000 | `beatty3` | 0.073372 | `beatty3` | 0.046157 | True | False | `902,261,595` |
| `k5_N1000000000_case0_P3_5_7_17_19` | 5 | 1000000000 | `layer_corrected` | 0.033939 | `layer_corrected` | 0.018068 | True | False | `55,36,39,13,29` |
| `k5_N1000000000_case1_P3_7_11_13_19` | 5 | 1000000000 | `layer_corrected` | 0.032905 | `layer_corrected` | 0.015197 | True | False | `38,77,16,23,16` |
| `k6_N1000000000_case0_P2_3_7_11_13_17` | 6 | 1000000000 | `xplusy_adaptive` | 0.052374 | `xplusy_corrected` | 0.027211 | True | False | `14,33,29,5,10,6` |
| `k6_N1000000000_case1_P2_3_5_11_13_19` | 6 | 1000000000 | `xplusy_corrected` | 0.039154 | `xplusy_corrected` | 0.021004 | True | False | `106,6,4,0,3,20` |
| `k8_N1000000000_case0_P2_3_5_7_11_13_17_19` | 8 | 1000000000 | `xplusy_corrected` | 0.079592 | `xplusy_corrected` | 0.057587 | True | False | `5,2,28,5,5,1,2,1` |
| `k3_N1000000000000_case0_P2_3_5` | 3 | 1000000000000 | `beatty3` | 0.628805 | `beatty3` | 0.605709 | True | False | `1126,16930,40` |
| `k3_N1000000000000_case1_P2_3_7` | 3 | 1000000000000 | `beatty3` | 0.494665 | `beatty3` | 0.467605 | True | False | `20814,2928,1578` |
| `k5_N1000000000000_case0_P2_3_5_11_17` | 5 | 1000000000000 | `layer_corrected` | 0.332698 | `layer_corrected` | 0.305411 | True | False | `400,112,101,93,74` |
| `k5_N1000000000000_case1_P3_5_11_13_17` | 5 | 1000000000000 | `layer_corrected` | 0.553059 | `layer_corrected` | 0.503198 | True | False | `124,147,147,57,149` |
| `k6_N1000000000000_case0_P3_5_7_13_17_19` | 6 | 1000000000000 | `xplusy_corrected` | 0.748380 | `xplusy_corrected` | 0.677723 | True | False | `83,125,58,14,23,34` |
| `k6_N1000000000000_case1_P2_5_7_11_13_19` | 6 | 1000000000000 | `xplusy_corrected` | 0.721063 | `xplusy_corrected` | 0.626598 | True | False | `161,12,56,88,21,14` |
| `k8_N1000000000000_case0_P2_3_5_7_11_13_17_19` | 8 | 1000000000000 | `xplusy_corrected` | 1.439791 | `xplusy_corrected` | 1.376875 | True | False | `75,28,9,16,3,22,5,1` |

## Best-Known Claim Gate

This artifact is not, by itself, a broad best-known claim. It is the
current implemented portfolio with per-row certification status.
Rows marked `audit blocked` had solver agreement but exceeded the
current interval auditor's numeric range. Remaining comparator gates
before using best-known language include a faithful soft-heap/selection
baseline, clearer Barvinok/fixed-dimensional lattice-count positioning,
and broader independent certification beyond the current audited prime
universe.
