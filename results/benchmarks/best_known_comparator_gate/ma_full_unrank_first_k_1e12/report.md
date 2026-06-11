# MA Full X+Y Unrank Suite

- Timestamp: `2026-06-11T21:55:44.154906+00:00`
- Git commit: `d3e2f716a92ef6e8270bb5b1bc346205f69f77f1`
- Git dirty: `False`
- N: `1000000000000`
- Cases: `3`
- Completed cases: `3`
- Same-exponent certified cases: `3`
- MA wall-time wins: `0`
- Mean MA/corrected wall ratio: `3.270836895451455`

Both rows return full exponent vectors. `xplusy_ma_full` uses the
Mirzaian-Arjomandi sorted-matrix selector to choose the X+Y log value,
then reconstructs and exact-sorts a narrow exponent band. Matching rows
are checked by the independent interval-rank auditor.

| P | corrected wall s | MA wall s | ratio | MA phase s | corrected RSS KB | MA RSS KB | MA n | MA band | certified | exps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| `2,3,5,7,11` | 2.429897 | 8.749338 | 3.600703 | 7.040366 | 428056 | 1763784 | 18731384 | 494 | True | `[1052, 26, 33, 53, 4]` |
| `2,3,5,7,11,13` | 1.110450 | 3.190004 | 2.872714 | 2.276573 | 208416 | 700160 | 5166940 | 503 | True | `[55, 126, 27, 54, 2, 52]` |
| `2,3,5,7,11,13,17,19` | 1.713758 | 5.722401 | 3.339094 | 3.935002 | 365272 | 1055432 | 8512638 | 480 | True | `[75, 28, 9, 16, 3, 22, 5, 1]` |

## Claim Status

This is a published sorted-matrix/X+Y selector comparison wrapped into
full exponent-vector unranking and independent rank auditing. It is a
narrow best-known-style comparator for the sorted-matrix value-selection
subproblem, not a claim against soft-heap top-k/output algorithms, full
Frederickson-Johnson ranking/selection, or Barvinok-style lattice
counting.
