# Sorted-Matrix / LOH X+Y Workbench

- Timestamp: `2026-06-14T22:27:33.893385+00:00`
- Git commit: `ed10e2d9677350cb7dfb919c2156b5065a7c2605`
- Git dirty: `True`
- N: `1000000`
- Cases: `3`
- Completed cases: `3`
- Range-pruning wins over linear saddleback count: `3`
- Mean block/linear internal time ratio: `0.8500797477344003`
- Mirzaian-Arjomandi probe wins over linear saddleback count: `0`
- Mean Mirzaian-Arjomandi/linear internal time ratio: `2.020283917554494`
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
| `2,3,5,7,11` | 0.000689 | 0.000562 | 0.815675 | 0.001353 | 1.963716 | 0 | 0.068465 | 31336 | 1000000 | 0.038560 |
| `2,3,5,7,11,13` | 0.000952 | 0.000826 | 0.867647 | 0.001295 | 1.360294 | 0 | 0.063629 | 31664 | 1000000 | 0.036452 |
| `2,3,5,7,11,13,17,19` | 0.001330 | 0.001153 | 0.866917 | 0.003640 | 2.736842 | 0 | 0.068529 | 31928 | 1000000 | 0.037350 |
