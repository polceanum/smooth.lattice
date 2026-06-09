# Sorted-Matrix / LOH X+Y Workbench

- Timestamp: `2026-06-09T05:27:41.188294+00:00`
- Git commit: `44175948b9b9bf3010f5fd388e2d926180fe5f45`
- Git dirty: `False`
- N: `1000000000000`
- Cases: `6`
- Completed cases: `6`
- Range-pruning wins over linear saddleback count: `2`
- Mean block/linear internal time ratio: `1.2765487937018578`
- Mirzaian-Arjomandi probe wins over linear saddleback count: `0`
- Mean Mirzaian-Arjomandi/linear internal time ratio: `9.050859617138036`

The `block_rank` rows are sorted-matrix range-pruning probes, not a faithful
Frederickson-Johnson implementation. The `ma_select_probe` row adapts the
Mirzaian-Arjomandi square sorted-matrix selector by padding the shorter MITM
side and selects only a log value. The `loh_topk_probe` row uses a capped
output-style top-k rank (`N_probe`) and is not evidence for random access at
the full rank when `N_probe != N`.

| P | linear total s | best block s | block/linear | MA s | MA/linear | MA delta | wall s | RSS KB | LOH probe N | LOH s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `2,3,5,7,11` | 1.062584 | 1.023656 | 0.963365 | 6.819068 | 6.417439 | 0 | 14.393216 | 1475764 | 1000000 | 0.043760 |
| `2,3,5,7,13` | 0.657743 | 0.810564 | 1.232342 | 5.562616 | 8.457127 | 0 | 9.453442 | 1441480 | 1000000 | 0.033762 |
| `2,3,5,11,13` | 0.712644 | 0.706035 | 0.990726 | 7.858682 | 11.027500 | 0 | 11.620346 | 1854044 | 1000000 | 0.029931 |
| `2,3,7,11,13` | 0.583441 | 0.786766 | 1.348493 | 4.107915 | 7.040840 | 0 | 7.733516 | 1348480 | 1000000 | 0.036617 |
| `2,5,7,11,13` | 0.749615 | 1.038837 | 1.385827 | 6.779213 | 9.043593 | 0 | 11.973836 | 1544148 | 1000000 | 0.038218 |
| `3,5,7,11,13` | 0.888561 | 1.544799 | 1.738540 | 10.945879 | 12.318658 | 0 | 18.697194 | 2031084 | 1000000 | 0.048415 |
