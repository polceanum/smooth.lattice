# Sorted-Matrix / LOH X+Y Workbench

- Timestamp: `2026-06-11T21:55:44.543104+00:00`
- Git commit: `d3e2f716a92ef6e8270bb5b1bc346205f69f77f1`
- Git dirty: `True`
- N: `1000000`
- Cases: `3`
- Completed cases: `3`
- Range-pruning wins over linear saddleback count: `3`
- Mean block/linear internal time ratio: `0.8458775068607555`
- Mirzaian-Arjomandi probe wins over linear saddleback count: `0`
- Mean Mirzaian-Arjomandi/linear internal time ratio: `2.2590266549340225`
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
| `2,3,5,7,11` | 0.000471 | 0.000414 | 0.878981 | 0.001244 | 2.641189 | 0 | 0.072726 | 31452 | 1000000 | 0.047909 |
| `2,3,5,7,11,13` | 0.001012 | 0.000770 | 0.760870 | 0.001454 | 1.436759 | 0 | 0.068597 | 31728 | 1000000 | 0.038262 |
| `2,3,5,7,11,13,17,19` | 0.001037 | 0.000931 | 0.897782 | 0.002799 | 2.699132 | 0 | 0.066691 | 31896 | 1000000 | 0.040184 |
