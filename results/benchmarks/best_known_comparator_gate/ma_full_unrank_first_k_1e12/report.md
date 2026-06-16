# MA Full X+Y Unrank Suite

- Timestamp: `2026-06-16T06:13:11.814112+00:00`
- Git commit: `3721e8168de0872ba378f9ece667e2e52df273e8`
- Git dirty: `False`
- N: `1000000000000`
- Cases: `3`
- Completed cases: `3`
- Same-exponent certified cases: `3`
- MA wall-time wins: `0`
- Mean MA/corrected wall ratio: `3.795348433408204`

Both rows return full exponent vectors. `xplusy_ma_full` uses the
Mirzaian-Arjomandi sorted-matrix selector to choose the X+Y log value,
then reconstructs and exact-sorts a narrow exponent band. Matching rows
are checked by the independent interval-rank auditor.

| P | corrected wall s | MA wall s | ratio | MA phase s | corrected RSS KB | MA RSS KB | MA n | MA band | certified | exps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| `2,3,5,7,11` | 1.847217 | 7.412623 | 4.012860 | 5.969936 | 427984 | 1763760 | 18731384 | 494 | True | `[1052, 26, 33, 53, 4]` |
| `2,3,5,7,11,13` | 0.693548 | 2.645139 | 3.813922 | 1.859122 | 208340 | 700112 | 5166940 | 503 | True | `[55, 126, 27, 54, 2, 52]` |
| `2,3,5,7,11,13,17,19` | 1.422181 | 5.061917 | 3.559263 | 3.460213 | 365212 | 1055444 | 8512638 | 480 | True | `[75, 28, 9, 16, 3, 22, 5, 1]` |

## Claim Status

This is a published sorted-matrix/X+Y selector comparison wrapped into
full exponent-vector unranking and independent rank auditing. It is a
narrow best-known-style comparator for the sorted-matrix value-selection
subproblem, not a claim against soft-heap top-k/output algorithms, full
Frederickson-Johnson ranking/selection, or Barvinok-style lattice
counting.
