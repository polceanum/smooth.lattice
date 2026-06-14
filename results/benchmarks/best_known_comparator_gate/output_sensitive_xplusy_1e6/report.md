# Sorted-Matrix / LOH X+Y Workbench

- Timestamp: `2026-06-14T22:34:51.532771+00:00`
- Git commit: `26366bf06aa59d88f29d3d3264578eda7992b045`
- Git dirty: `True`
- N: `1000000`
- Cases: `3`
- Completed cases: `3`
- Range-pruning wins over linear saddleback count: `3`
- Mean block/linear internal time ratio: `0.8690779105459225`
- Mirzaian-Arjomandi probe wins over linear saddleback count: `0`
- Mean Mirzaian-Arjomandi/linear internal time ratio: `1.8715664428483259`
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
| `2,3,5,7,11` | 0.000598 | 0.000529 | 0.884615 | 0.001114 | 1.862876 | 0 | 0.067483 | 31408 | 1000000 | 0.042269 |
| `2,3,5,7,11,13` | 0.000933 | 0.000866 | 0.928189 | 0.001335 | 1.430868 | 0 | 0.067304 | 31632 | 1000000 | 0.037976 |
| `2,3,5,7,11,13,17,19` | 0.001508 | 0.001198 | 0.794430 | 0.003500 | 2.320955 | 0 | 0.070687 | 31800 | 1000000 | 0.038652 |
