# Sorted-Matrix / LOH X+Y Workbench

- Timestamp: `2026-06-15T22:10:49.409709+00:00`
- Git commit: `992256f6ea695dbbe0f495bcbac2f37b279574b6`
- Git dirty: `True`
- N: `1000000000000`
- Cases: `3`
- Completed cases: `3`
- Range-pruning wins over linear saddleback count: `0`
- Mean block/linear internal time ratio: `1.290216118017839`
- Mirzaian-Arjomandi probe wins over linear saddleback count: `0`
- Mean Mirzaian-Arjomandi/linear internal time ratio: `5.887569848479676`
- Mat-Select2 heap-primitive wins over linear saddleback count: `0`
- Mean Mat-Select2 heap-primitive/linear internal time ratio: `13.0590293637715`
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
| `2,3,5,7,11` | 0.619573 | 0.944710 | 1.524776 | 5.934710 | 9.578710 | 8.091022 | 13.059029 | 0 | 18.613790 | 1475720 | 1000000 | 0.018353 |
| `2,3,5,7,11,13` | 0.467507 | 0.541168 | 1.157561 | 1.488711 | 3.184361 | 0.000000 |  | 0 | 4.205567 | 511380 | 1000000 | 0.029940 |
| `2,3,5,7,11,13,17,19` | 0.797101 | 0.947204 | 1.188311 | 3.905507 | 4.899639 | 0.000000 |  | 0 | 8.651832 | 1224820 | 1000000 | 0.032929 |
