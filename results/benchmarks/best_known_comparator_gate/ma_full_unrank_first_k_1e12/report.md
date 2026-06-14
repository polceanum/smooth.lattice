# MA Full X+Y Unrank Suite

- Timestamp: `2026-06-14T22:27:33.492120+00:00`
- Git commit: `ed10e2d9677350cb7dfb919c2156b5065a7c2605`
- Git dirty: `True`
- N: `1000000000000`
- Cases: `3`
- Completed cases: `3`
- Same-exponent certified cases: `3`
- MA wall-time wins: `0`
- Mean MA/corrected wall ratio: `3.5279214370995486`

Both rows return full exponent vectors. `xplusy_ma_full` uses the
Mirzaian-Arjomandi sorted-matrix selector to choose the X+Y log value,
then reconstructs and exact-sorts a narrow exponent band. Matching rows
are checked by the independent interval-rank auditor.

| P | corrected wall s | MA wall s | ratio | MA phase s | corrected RSS KB | MA RSS KB | MA n | MA band | certified | exps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| `2,3,5,7,11` | 1.603780 | 6.200188 | 3.865984 | 4.860684 | 427948 | 1763780 | 18731384 | 494 | True | `[1052, 26, 33, 53, 4]` |
| `2,3,5,7,11,13` | 0.615101 | 2.147090 | 3.490627 | 1.476059 | 208340 | 700096 | 5166940 | 503 | True | `[55, 126, 27, 54, 2, 52]` |
| `2,3,5,7,11,13,17,19` | 1.371509 | 4.426070 | 3.227153 | 2.910884 | 365228 | 1055476 | 8512638 | 480 | True | `[75, 28, 9, 16, 3, 22, 5, 1]` |

## Claim Status

This is a published sorted-matrix/X+Y selector comparison wrapped into
full exponent-vector unranking and independent rank auditing. It is a
narrow best-known-style comparator for the sorted-matrix value-selection
subproblem, not a claim against soft-heap top-k/output algorithms, full
Frederickson-Johnson ranking/selection, or Barvinok-style lattice
counting.
