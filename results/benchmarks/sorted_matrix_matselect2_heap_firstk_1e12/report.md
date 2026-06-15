# Sorted-Matrix / LOH X+Y Workbench

- Timestamp: `2026-06-15T22:17:12.893215+00:00`
- Git commit: `a58e6f286e9eabc895ce36010fc1c4063dd506bf`
- Git dirty: `False`
- N: `1000000000000`
- Cases: `3`
- Completed cases: `3`
- Range-pruning wins over linear saddleback count: `0`
- Mean block/linear internal time ratio: `1.3127183898439034`
- Mirzaian-Arjomandi probe wins over linear saddleback count: `0`
- Mean Mirzaian-Arjomandi/linear internal time ratio: `5.828476470644191`
- Mat-Select2 heap-primitive wins over linear saddleback count: `0`
- Mean Mat-Select2 heap-primitive/linear internal time ratio: `12.782335609420839`
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
| `2,3,5,7,11` | 0.628564 | 0.941070 | 1.497175 | 5.956272 | 9.475999 | 8.034516 | 12.782336 | 0 | 18.522556 | 1475740 | 1000000 | 0.020230 |
| `2,3,5,7,11,13` | 0.480929 | 0.587816 | 1.222251 | 1.448888 | 3.012686 | 0.000000 |  | 0 | 4.207385 | 511392 | 1000000 | 0.028840 |
| `2,3,5,7,11,13,17,19` | 0.782624 | 0.953807 | 1.218730 | 3.910572 | 4.996744 | 0.000000 |  | 0 | 8.625675 | 1224760 | 1000000 | 0.028789 |
