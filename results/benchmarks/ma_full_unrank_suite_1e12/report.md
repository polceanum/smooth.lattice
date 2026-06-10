# MA Full X+Y Unrank Suite

- Timestamp: `2026-06-10T20:54:52.095537+00:00`
- Git commit: `472bb6d21f866315de868606ab163db53d49572f`
- Git dirty: `False`
- N: `1000000000000`
- Cases: `6`
- Completed cases: `6`
- Same-exponent certified cases: `6`
- MA wall-time wins: `0`
- Mean MA/corrected wall ratio: `5.69995979563249`

Both rows return full exponent vectors. `xplusy_ma_full` uses the
Mirzaian-Arjomandi sorted-matrix selector to choose the X+Y log value,
then reconstructs and exact-sorts a narrow exponent band. Matching rows
are checked by the independent interval-rank auditor.

| P | corrected wall s | MA wall s | ratio | MA phase s | corrected RSS KB | MA RSS KB | MA n | MA band | certified | exps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| `2,3,5,7,11` | 2.306325 | 10.327856 | 4.478058 | 8.076410 | 428088 | 1427020 | 18731384 | 494 | True | `[1052, 26, 33, 53, 4]` |
| `2,3,5,7,13` | 1.462174 | 8.762559 | 5.992829 | 7.048237 | 417460 | 1717700 | 18234176 | 496 | True | `[205, 279, 119, 131, 16]` |
| `2,3,5,11,13` | 1.289603 | 6.511644 | 5.049339 | 4.882822 | 386356 | 1711596 | 16774575 | 500 | True | `[254, 220, 258, 4, 52]` |
| `2,3,7,11,13` | 1.297363 | 6.583541 | 5.074554 | 5.174579 | 360480 | 1589164 | 15549479 | 504 | True | `[291, 331, 90, 84, 28]` |
| `2,5,7,11,13` | 1.469983 | 9.438766 | 6.421002 | 7.814691 | 445448 | 1839028 | 19543443 | 486 | True | `[106, 306, 164, 53, 32]` |
| `3,5,7,11,13` | 2.019782 | 14.510070 | 7.183977 | 12.229094 | 579736 | 2020080 | 25758466 | 505 | True | `[273, 219, 5, 219, 5]` |
