# Certified Portfolio Suite

- Timestamp: `2026-06-10T18:45:17.746380+00:00`
- Git commit: `251606b886061eac4944af272728207ceba7b366`
- Git dirty: `False`
- Cases: `14`
- Certified cases: `12`
- Solver-agreement cases: `14`
- Audit-blocked agreement cases: `2`
- Solver-disagreement cases: `0`
- Wall-time winners: `{'beatty3': 1, 'layer_corrected': 4, 'xplusy_adaptive': 1, 'xplusy_corrected': 6}`
- Reported-time winners: `{'beatty3': 1, 'layer_corrected': 4, 'sums_corrected': 1, 'xplusy_corrected': 6}`

This suite is a portfolio benchmark over implemented full-unrank methods.
A row is certified only when successful methods agree on the exponent vector
and the independent interval auditor proves `count_le(candidate)=N`.

| label | k | N | fastest wall | wall sec | fastest reported | reported sec | cert | audit blocked | exps |
|---|---:|---:|---|---:|---|---:|---|---|---|
| `k3_N1000000000_case0_P2_3_5` | 3 | 1000000000 | `xplusy_adaptive` | 0.183933 | `sums_corrected` | 0.048871 | True | False | `1334,335,404` |
| `k3_N1000000000_case1_P2_3_7` | 3 | 1000000000 | `beatty3` | 0.063579 | `beatty3` | 0.041485 | True | False | `902,261,595` |
| `k5_N1000000000_case0_P3_5_7_17_19` | 5 | 1000000000 | `layer_corrected` | 0.031123 | `layer_corrected` | 0.016300 | True | False | `55,36,39,13,29` |
| `k5_N1000000000_case1_P3_7_11_13_19` | 5 | 1000000000 | `layer_corrected` | 0.035107 | `layer_corrected` | 0.014288 | True | False | `38,77,16,23,16` |
| `k6_N1000000000_case0_P2_3_7_11_13_17` | 6 | 1000000000 | `xplusy_corrected` | 0.043746 | `xplusy_corrected` | 0.021881 | True | False | `14,33,29,5,10,6` |
| `k6_N1000000000_case1_P2_3_5_11_13_19` | 6 | 1000000000 | `xplusy_corrected` | 0.037561 | `xplusy_corrected` | 0.017733 | True | False | `106,6,4,0,3,20` |
| `k8_N1000000000_case0_P2_3_5_7_11_13_17_19` | 8 | 1000000000 | `xplusy_corrected` | 0.109875 | `xplusy_corrected` | 0.073714 | True | False | `5,2,28,5,5,1,2,1` |
| `k3_N1000000000000_case0_P2_3_5` | 3 | 1000000000000 | `beatty3` | 0.559972 | `beatty3` | 0.533659 | False | True | `1126,16930,40` |
| `k3_N1000000000000_case1_P2_3_7` | 3 | 1000000000000 | `beatty3` | 0.497679 | `beatty3` | 0.474644 | False | True | `20814,2928,1578` |
| `k5_N1000000000000_case0_P2_3_5_11_17` | 5 | 1000000000000 | `layer_corrected` | 0.340595 | `layer_corrected` | 0.325717 | True | False | `400,112,101,93,74` |
| `k5_N1000000000000_case1_P3_5_11_13_17` | 5 | 1000000000000 | `layer_corrected` | 0.368788 | `layer_corrected` | 0.347340 | True | False | `124,147,147,57,149` |
| `k6_N1000000000000_case0_P3_5_7_13_17_19` | 6 | 1000000000000 | `xplusy_corrected` | 0.699677 | `xplusy_corrected` | 0.627076 | True | False | `83,125,58,14,23,34` |
| `k6_N1000000000000_case1_P2_5_7_11_13_19` | 6 | 1000000000000 | `xplusy_corrected` | 0.762898 | `xplusy_corrected` | 0.677287 | True | False | `161,12,56,88,21,14` |
| `k8_N1000000000000_case0_P2_3_5_7_11_13_17_19` | 8 | 1000000000000 | `xplusy_corrected` | 1.546891 | `xplusy_corrected` | 1.463303 | True | False | `75,28,9,16,3,22,5,1` |

## Best-Known Claim Gate

This artifact is not, by itself, a broad best-known claim. It is the
current implemented portfolio with per-row certification status.
Rows marked `audit blocked` had solver agreement but exceeded the
current interval auditor's numeric range. Remaining comparator gates
before using best-known language include a faithful soft-heap/selection
baseline, clearer Barvinok/fixed-dimensional lattice-count positioning,
and broader independent certification beyond the current audited prime
universe.
