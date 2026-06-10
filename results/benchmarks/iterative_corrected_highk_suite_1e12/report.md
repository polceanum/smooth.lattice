# Iterative Corrected High-k Suite

- Timestamp: `2026-06-10T22:46:18.937459+00:00`
- Git commit: `dee2ec6e5ce7ebfe8a644746eb3afe0560bcbd99`
- Git dirty: `False`
- N: `1000000000000`
- Cases: `2`
- Same-exponent certified cases: `2`
- Corrected wall-time wins: `2`
- Mean adaptive/corrected wall ratio: `2.270174667714107`

The corrected row applies multiple exact residual corrections to the
analytic lattice estimate before opening the final exact-sorted band.
Every corrected output is independently interval-rank certified.

| label | k | radius | refine | adaptive wall s | corrected wall s | ratio | band | cands | audit s | certified | exps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| `k10_first10` | 10 | 50 | 2 | 5.771289 | 3.006165 | 1.919818 | 99 | 99 | 11.560988 | True | `[71, 29, 3, 1, 10, 2, 1, 2, 0, 3]` |
| `k12_first12` | 12 | 25 | 3 | 12.643785 | 4.824893 | 2.620531 | 45 | 45 | 16.282043 | True | `[27, 2, 9, 1, 2, 3, 5, 4, 1, 0, 1, 6]` |
