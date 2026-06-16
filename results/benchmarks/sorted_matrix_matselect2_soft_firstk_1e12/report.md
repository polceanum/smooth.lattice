# Sorted-Matrix / LOH X+Y Workbench

- Timestamp: `2026-06-16T06:07:26.921351+00:00`
- Git commit: `799473167e787b7ac205e27525d7acaef79dcf05`
- Git dirty: `False`
- N: `1000000000000`
- Cases: `3`
- Completed cases: `3`
- Range-pruning wins over linear saddleback count: `0`
- Mean block/linear internal time ratio: `1.2836936332563065`
- Mirzaian-Arjomandi probe wins over linear saddleback count: `0`
- Mean Mirzaian-Arjomandi/linear internal time ratio: `5.835195566983406`
- Mat-Select2 heap-primitive wins over linear saddleback count: `0`
- Mean Mat-Select2 heap-primitive/linear internal time ratio: `13.764194582533719`
- Mat-Select2 soft-selector wins over linear saddleback count: `0`
- Mean Mat-Select2 soft-selector/linear internal time ratio: `66.91679515968262`
- Mirzaian-Arjomandi exhaustive validation cases: `5185`
- Mirzaian-Arjomandi exhaustive validation failures: `0`
- Mirzaian-Arjomandi exhaustive validation max delta: `0.0`
- Mat-Select2 heap-primitive exhaustive validation cases: `3233`
- Mat-Select2 heap-primitive exhaustive validation failures: `0`
- Mat-Select2 heap-primitive exhaustive validation max delta: `0.0`
- Mat-Select2 soft-selector exhaustive validation cases: `3233`
- Mat-Select2 soft-selector exhaustive validation failures: `0`
- Mat-Select2 soft-selector exhaustive validation max delta: `0.0`

The `block_rank` rows are sorted-matrix range-pruning probes, not a faithful
Frederickson-Johnson implementation. The `ma_select_probe` row adapts the
Mirzaian-Arjomandi square sorted-matrix selector by padding the shorter MITM
side and selects only a log value. The `matselect2_heap_probe` row implements
the exponential-block selector from the Kaplan/Frederickson-Johnson line of
algorithms with an exact binary heap for the primitive that is soft-heap based
in the asymptotic paper algorithm. The `matselect2_soft_probe` row uses the
same exponential-block structure but replaces the row-list primitive with the
validated soft-sequence-heap Soft-Select bridge; it remains a practical probe,
not a proof of the published time bound.
The `loh_topk_probe` row uses a capped
output-style top-k rank (`N_probe`) and is not evidence for random access at
the full rank when `N_probe != N`.

| P | linear total s | best block s | block/linear | MA s | MA/linear | Mat2 heap s | Mat2 heap/linear | Mat2 soft s | Mat2 soft/linear | Mat2 soft delta | wall s | RSS KB | LOH probe N | LOH s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `2,3,5,7,11` | 0.631033 | 1.002461 | 1.588603 | 6.061583 | 9.605810 | 8.685661 | 13.764195 | 42.226706 | 66.916795 | 0 | 61.651052 | 1475704 | 1000000 | 0.028383 |
| `2,3,5,7,11,13` | 0.491876 | 0.572366 | 1.163639 | 1.547035 | 3.145173 | 0.000000 |  | 0.000000 |  | 0 | 4.347479 | 511380 | 1000000 | 0.027605 |
| `2,3,5,7,11,13,17,19` | 0.857769 | 0.942550 | 1.098839 | 4.078352 | 4.754604 | 0.000000 |  | 0.000000 |  | 0 | 9.026434 | 1224704 | 1000000 | 0.033390 |
