# Certified Full X+Y Five-Prime Suite

- Timestamp: `2026-06-09T06:07:43.168058+00:00`
- Git commit: `5314fd8a56fb7ca760046076e1f7f168fd48386a`
- Git dirty: `False`
- N: `1000000000000`
- Cases: `6`
- Same-exponent certified cases: `6`
- Layer wall-time wins: `6`
- Mean wall ratio full X+Y/layer: `2.327315382062943`

Both compared rows return full exponent vectors. Each successful case audits both outputs
with the independent interval-rank auditor and checks that the exponent vectors match.
The `layer_compressed` row records whether the asymptotic analytic bracket was used.

| P | layer wall s | full X+Y wall s | ratio | layer calls | layer gap | analytic | layer RSS KB | X+Y RSS KB | same exps | both certified | exps |
|---|---:|---:|---:|---:|---:|---|---:|---:|---|---|---|
| `2,3,5,7,11` | 1.276416 | 2.628903 | 2.059597 | 11 | 340 | True | 23036 | 460168 | True | True | `[1052, 26, 33, 53, 4]` |
| `2,3,5,7,13` | 0.696232 | 1.782958 | 2.560867 | 11 | 336 | True | 23536 | 448772 | True | True | `[205, 279, 119, 131, 16]` |
| `2,3,5,11,13` | 0.811143 | 1.595716 | 1.967243 | 11 | 334 | True | 24988 | 416116 | True | True | `[254, 220, 258, 4, 52]` |
| `2,3,7,11,13` | 0.750894 | 1.530324 | 2.038003 | 11 | 366 | True | 22824 | 387612 | True | True | `[291, 331, 90, 84, 28]` |
| `2,5,7,11,13` | 0.704921 | 1.920789 | 2.724829 | 11 | 341 | True | 19096 | 479564 | True | True | `[106, 306, 164, 53, 32]` |
| `3,5,7,11,13` | 0.942377 | 2.462765 | 2.613353 | 11 | 358 | True | 22808 | 623608 | True | True | `[273, 219, 5, 219, 5]` |
