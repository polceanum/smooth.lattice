# Sorted-Matrix / LOH X+Y Workbench

- Timestamp: `2026-06-15T22:23:35.125465+00:00`
- Git commit: `3c852e9b88e955b6e4e8de6d55b9c90a2476c703`
- Git dirty: `True`
- N: `1000000`
- Cases: `3`
- Completed cases: `3`
- Range-pruning wins over linear saddleback count: `2`
- Mean block/linear internal time ratio: `0.9566811987145698`
- Mirzaian-Arjomandi probe wins over linear saddleback count: `0`
- Mean Mirzaian-Arjomandi/linear internal time ratio: `2.007773448900695`
- Mat-Select2 heap-primitive wins over linear saddleback count: `0`
- Mean Mat-Select2 heap-primitive/linear internal time ratio: `10.167390045375946`
- Mirzaian-Arjomandi exhaustive validation cases: `5185`
- Mirzaian-Arjomandi exhaustive validation failures: `0`
- Mirzaian-Arjomandi exhaustive validation max delta: `0.0`
- Mat-Select2 heap-primitive exhaustive validation cases: `3233`
- Mat-Select2 heap-primitive exhaustive validation failures: `0`
- Mat-Select2 heap-primitive exhaustive validation max delta: `0.0`

The `block_rank` rows are sorted-matrix range-pruning probes, not a faithful
Frederickson-Johnson implementation. The `ma_select_probe` row adapts the
Mirzaian-Arjomandi square sorted-matrix selector by padding the shorter MITM
side and selects only a log value. The `matselect2_heap_probe` row implements
the exponential-block selector from the Kaplan/Frederickson-Johnson line of
algorithms with an exact binary heap for the primitive that is soft-heap based
in the asymptotic paper algorithm. The `loh_topk_probe` row uses a capped
output-style top-k rank (`N_probe`) and is not evidence for random access at
the full rank when `N_probe != N`.

| P | linear total s | best block s | block/linear | MA s | MA/linear | Mat2 s | Mat2/linear | Mat2 delta | wall s | RSS KB | LOH probe N | LOH s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `2,3,5,7,11` | 0.000481 | 0.000512 | 1.064449 | 0.001196 | 2.486486 | 0.003734 | 7.762994 | 0 | 0.071988 | 31492 | 1000000 | 0.045896 |
| `2,3,5,7,11,13` | 0.000894 | 0.000778 | 0.870246 | 0.001160 | 1.297539 | 0.009513 | 10.640940 | 0 | 0.070075 | 31864 | 1000000 | 0.034402 |
| `2,3,5,7,11,13,17,19` | 0.001191 | 0.001114 | 0.935348 | 0.002667 | 2.239295 | 0.014409 | 12.098237 | 0 | 0.081758 | 31808 | 1000000 | 0.035475 |
