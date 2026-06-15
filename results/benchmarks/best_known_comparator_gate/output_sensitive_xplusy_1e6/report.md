# Sorted-Matrix / LOH X+Y Workbench

- Timestamp: `2026-06-15T22:51:44.528191+00:00`
- Git commit: `aaef2d09004798d1719614567c71f6f325f0cfe6`
- Git dirty: `True`
- N: `1000000`
- Cases: `3`
- Completed cases: `3`
- Range-pruning wins over linear saddleback count: `2`
- Mean block/linear internal time ratio: `0.9531453682122057`
- Mirzaian-Arjomandi probe wins over linear saddleback count: `0`
- Mean Mirzaian-Arjomandi/linear internal time ratio: `2.848539503864853`
- Mat-Select2 heap-primitive wins over linear saddleback count: `0`
- Mean Mat-Select2 heap-primitive/linear internal time ratio: `9.701120596922914`
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
| `2,3,5,7,11` | 0.000461 | 0.000480 | 1.041215 | 0.001341 | 2.908894 | 0.004066 | 8.819957 | 0 | 0.070372 | 31504 | 1000000 | 0.043054 |
| `2,3,5,7,11,13` | 0.000826 | 0.000775 | 0.938257 | 0.001454 | 1.760291 | 0.007367 | 8.918886 | 0 | 0.069434 | 31780 | 1000000 | 0.037839 |
| `2,3,5,7,11,13,17,19` | 0.001133 | 0.000997 | 0.879965 | 0.004392 | 3.876434 | 0.012876 | 11.364519 | 0 | 0.081690 | 31960 | 1000000 | 0.038988 |
