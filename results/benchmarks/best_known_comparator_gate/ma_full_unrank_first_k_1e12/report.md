# MA Full X+Y Unrank Suite

- Timestamp: `2026-06-15T22:23:34.668529+00:00`
- Git commit: `3c852e9b88e955b6e4e8de6d55b9c90a2476c703`
- Git dirty: `False`
- N: `1000000000000`
- Cases: `3`
- Completed cases: `3`
- Same-exponent certified cases: `3`
- MA wall-time wins: `0`
- Mean MA/corrected wall ratio: `3.8331755528502947`

Both rows return full exponent vectors. `xplusy_ma_full` uses the
Mirzaian-Arjomandi sorted-matrix selector to choose the X+Y log value,
then reconstructs and exact-sorts a narrow exponent band. Matching rows
are checked by the independent interval-rank auditor.

| P | corrected wall s | MA wall s | ratio | MA phase s | corrected RSS KB | MA RSS KB | MA n | MA band | certified | exps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| `2,3,5,7,11` | 1.798728 | 7.412142 | 4.120769 | 5.992245 | 428008 | 1763776 | 18731384 | 494 | True | `[1052, 26, 33, 53, 4]` |
| `2,3,5,7,11,13` | 0.675758 | 2.564762 | 3.795383 | 1.792432 | 208328 | 700092 | 5166940 | 503 | True | `[55, 126, 27, 54, 2, 52]` |
| `2,3,5,7,11,13,17,19` | 1.355144 | 4.855989 | 3.583374 | 3.361469 | 365220 | 1055424 | 8512638 | 480 | True | `[75, 28, 9, 16, 3, 22, 5, 1]` |

## Claim Status

This is a published sorted-matrix/X+Y selector comparison wrapped into
full exponent-vector unranking and independent rank auditing. It is a
narrow best-known-style comparator for the sorted-matrix value-selection
subproblem, not a claim against soft-heap top-k/output algorithms, full
Frederickson-Johnson ranking/selection, or Barvinok-style lattice
counting.
