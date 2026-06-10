# Sorted-Matrix / LOH X+Y Workbench

- Timestamp: `2026-06-10T20:16:47.142378+00:00`
- Git commit: `17a7e40ae790ae9017ab8274204eaaf107add212`
- Git dirty: `False`
- N: `1000000000000`
- Cases: `6`
- Completed cases: `6`
- Range-pruning wins over linear saddleback count: `1`
- Mean block/linear internal time ratio: `1.3467233146527213`
- Mirzaian-Arjomandi probe wins over linear saddleback count: `0`
- Mean Mirzaian-Arjomandi/linear internal time ratio: `9.404243248473163`
- Mirzaian-Arjomandi exhaustive validation cases: `5185`
- Mirzaian-Arjomandi exhaustive validation failures: `0`
- Mirzaian-Arjomandi exhaustive validation max delta: `0.0`

The `block_rank` rows are sorted-matrix range-pruning probes, not a faithful
Frederickson-Johnson implementation. The `ma_select_probe` row adapts the
Mirzaian-Arjomandi square sorted-matrix selector by padding the shorter MITM
side and selects only a log value. The `loh_topk_probe` row uses a capped
output-style top-k rank (`N_probe`) and is not evidence for random access at
the full rank when `N_probe != N`.

| P | linear total s | best block s | block/linear | MA s | MA/linear | MA delta | wall s | RSS KB | LOH probe N | LOH s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `2,3,5,7,11` | 0.660537 | 1.004498 | 1.520729 | 6.534540 | 9.892769 | 0 | 11.495860 | 1475752 | 1000000 | 0.028856 |
| `2,3,5,7,13` | 0.690694 | 0.629582 | 0.911521 | 4.566296 | 6.611171 | 0 | 7.982556 | 1441576 | 1000000 | 0.027748 |
| `2,3,5,11,13` | 0.531501 | 0.711579 | 1.338810 | 6.023064 | 11.332178 | 0 | 9.558280 | 1854040 | 1000000 | 0.033946 |
| `2,3,7,11,13` | 0.457047 | 0.684666 | 1.498021 | 4.437861 | 9.709857 | 0 | 7.802844 | 1348528 | 1000000 | 0.038534 |
| `2,5,7,11,13` | 0.673154 | 1.062560 | 1.578480 | 6.700831 | 9.954380 | 0 | 11.690500 | 1544220 | 1000000 | 0.032943 |
| `3,5,7,11,13` | 0.888153 | 1.094896 | 1.232779 | 7.926858 | 8.925104 | 0 | 13.443623 | 2030928 | 1000000 | 0.033831 |
