# Certified Full X+Y Five-Prime Suite

- Timestamp: `2026-06-08T22:38:45.402216+00:00`
- Git commit: `13cba09636b22e139046442b9faa4ea0f5330377`
- Git dirty: `False`
- N: `1000000000000`
- Cases: `6`
- Same-exponent certified cases: `6`
- Layer wall-time wins: `6`
- Mean wall ratio full X+Y/layer: `2.3309099897793213`

Both compared rows return full exponent vectors. Each successful case audits both outputs
with the independent interval-rank auditor and checks that the exponent vectors match.

| P | layer wall s | full X+Y wall s | ratio | layer RSS KB | X+Y RSS KB | same exps | both certified | exps |
|---|---:|---:|---:|---:|---:|---|---|---|
| `2,3,5,7,11` | 1.340342 | 2.949791 | 2.200775 | 23276 | 460136 | True | True | `[1052, 26, 33, 53, 4]` |
| `2,3,5,7,13` | 0.715908 | 1.966326 | 2.746618 | 23728 | 448784 | True | True | `[205, 279, 119, 131, 16]` |
| `2,3,5,11,13` | 1.149269 | 1.689815 | 1.470340 | 31256 | 416108 | True | True | `[254, 220, 258, 4, 52]` |
| `2,3,7,11,13` | 0.764076 | 1.518400 | 1.987237 | 29056 | 387596 | True | True | `[291, 331, 90, 84, 28]` |
| `2,5,7,11,13` | 0.770630 | 1.890280 | 2.452901 | 26244 | 479616 | True | True | `[106, 306, 164, 53, 32]` |
| `3,5,7,11,13` | 0.879796 | 2.751641 | 3.127590 | 30872 | 623608 | True | True | `[273, 219, 5, 219, 5]` |
