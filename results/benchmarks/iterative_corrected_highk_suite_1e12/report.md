# Iterative Corrected High-k Suite

- Timestamp: `2026-06-10T23:44:14.676493+00:00`
- Git commit: `e0825a6f7ed1cbadaf823897f2b918037a6c7be9`
- Git dirty: `False`
- N: `1000000000000`
- Cases: `3`
- Same-exponent certified cases: `3`
- Corrected wall-time wins: `3`
- Mean adaptive/corrected wall ratio: `1.9694800865681135`

The corrected row applies a shared iterative policy: four exact residual
corrections and a rank-radius-25 final exact-sorted band.
Every corrected output is independently interval-rank certified.

| label | k | radius | refine | adaptive wall s | corrected wall s | ratio | band | cands | audit s | certified | exps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| `k8_first8` | 8 | 25 | 4 | 2.034849 | 1.093697 | 1.860523 | 52 | 52 | 4.645253 | True | `[75, 28, 9, 16, 3, 22, 5, 1]` |
| `k10_first10` | 10 | 25 | 4 | 3.436380 | 2.004847 | 1.714036 | 50 | 50 | 7.694640 | True | `[71, 29, 3, 1, 10, 2, 1, 2, 0, 3]` |
| `k12_first12` | 12 | 25 | 4 | 7.970556 | 3.415151 | 2.333881 | 38 | 38 | 11.134256 | True | `[27, 2, 9, 1, 2, 3, 5, 4, 1, 0, 1, 6]` |
