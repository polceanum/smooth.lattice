# Sorted-Matrix / LOH X+Y Workbench

- Timestamp: `2026-06-08T23:22:36.225372+00:00`
- Git commit: `674746fc15fd5f1cad31773bee81ba5a13ff13ea`
- Git dirty: `False`
- N: `1000000000000`
- Cases: `6`
- Completed cases: `6`
- Range-pruning wins over linear saddleback count: `2`
- Mean block/linear internal time ratio: `1.2542418501101371`

The `block_rank` rows are sorted-matrix range-pruning probes, not a faithful
Frederickson-Johnson implementation. The `loh_topk_probe` row uses a capped
output-style top-k rank (`N_probe`) and is not evidence for random access at
the full rank when `N_probe != N`.

| P | linear total s | best block s | block/linear | best leaf | wall s | RSS KB | LOH probe N | LOH s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `2,3,5,7,11` | 0.621935 | 0.972390 | 1.563491 | 2048 | 4.701049 | 183904 | 1000000 | 0.046326 |
| `2,3,5,7,13` | 0.598310 | 0.750283 | 1.254004 | 8192 | 3.532873 | 180316 | 1000000 | 0.042101 |
| `2,3,5,11,13` | 0.717543 | 0.703218 | 0.980036 | 8192 | 3.612273 | 169040 | 1000000 | 0.050068 |
| `2,3,7,11,13` | 0.598486 | 0.585925 | 0.979012 | 8192 | 3.052502 | 159444 | 1000000 | 0.033767 |
| `2,5,7,11,13` | 0.595853 | 0.936939 | 1.572433 | 8192 | 4.320058 | 190268 | 1000000 | 0.053414 |
| `3,5,7,11,13` | 0.947450 | 1.114651 | 1.176475 | 8192 | 5.354844 | 238248 | 1000000 | 0.045728 |
