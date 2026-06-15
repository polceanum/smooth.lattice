# Sorted-Matrix / LOH X+Y Workbench

- Timestamp: `2026-06-15T22:42:09.470152+00:00`
- Git commit: `2263264ccb755499cf7b7c9b9a74fc3a334975c5`
- Git dirty: `True`
- N: `1000000`
- Cases: `3`
- Completed cases: `3`
- Range-pruning wins over linear saddleback count: `3`
- Mean block/linear internal time ratio: `0.8825121289281279`
- Mirzaian-Arjomandi probe wins over linear saddleback count: `0`
- Mean Mirzaian-Arjomandi/linear internal time ratio: `2.174870800317018`
- Mat-Select2 heap-primitive wins over linear saddleback count: `0`
- Mean Mat-Select2 heap-primitive/linear internal time ratio: `9.620167656669393`
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
| `2,3,5,7,11` | 0.000497 | 0.000449 | 0.903421 | 0.000800 | 1.609658 | 0.003803 | 7.651911 | 0 | 0.071455 | 31504 | 1000000 | 0.042071 |
| `2,3,5,7,11,13` | 0.000796 | 0.000780 | 0.979899 | 0.002333 | 2.930905 | 0.008184 | 10.281407 | 0 | 0.071392 | 31860 | 1000000 | 0.035434 |
| `2,3,5,7,11,13,17,19` | 0.001442 | 0.001102 | 0.764216 | 0.002861 | 1.984050 | 0.015757 | 10.927184 | 0 | 0.084232 | 31908 | 1000000 | 0.037855 |
