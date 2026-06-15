# Sorted-Matrix / LOH X+Y Workbench

- Timestamp: `2026-06-15T22:19:35.665821+00:00`
- Git commit: `43e2c5e3a3ab3be2b50cb81826b005241df82f51`
- Git dirty: `True`
- N: `1000000`
- Cases: `3`
- Completed cases: `3`
- Range-pruning wins over linear saddleback count: `2`
- Mean block/linear internal time ratio: `0.8435155867637523`
- Mirzaian-Arjomandi probe wins over linear saddleback count: `0`
- Mean Mirzaian-Arjomandi/linear internal time ratio: `1.7874088092643674`
- Mat-Select2 heap-primitive wins over linear saddleback count: `0`
- Mean Mat-Select2 heap-primitive/linear internal time ratio: `10.231261425814914`
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
| `2,3,5,7,11` | 0.000479 | 0.000464 | 0.968685 | 0.001119 | 2.336117 | 0.004652 | 9.711900 | 0 | 0.074515 | 31380 | 1000000 | 0.045781 |
| `2,3,5,7,11,13` | 0.000742 | 0.000751 | 1.012129 | 0.001228 | 1.654987 | 0.009108 | 12.274933 | 0 | 0.068305 | 31812 | 1000000 | 0.036223 |
| `2,3,5,7,11,13,17,19` | 0.001870 | 0.001028 | 0.549733 | 0.002564 | 1.371123 | 0.016282 | 8.706952 | 0 | 0.083315 | 31916 | 1000000 | 0.035970 |
