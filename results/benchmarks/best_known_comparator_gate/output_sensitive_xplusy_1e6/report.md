# Sorted-Matrix / LOH X+Y Workbench

- Timestamp: `2026-06-16T06:13:12.280976+00:00`
- Git commit: `3721e8168de0872ba378f9ece667e2e52df273e8`
- Git dirty: `True`
- N: `1000000`
- Cases: `3`
- Completed cases: `3`
- Range-pruning wins over linear saddleback count: `3`
- Mean block/linear internal time ratio: `0.7815497376562345`
- Mirzaian-Arjomandi probe wins over linear saddleback count: `0`
- Mean Mirzaian-Arjomandi/linear internal time ratio: `1.4495965444620005`
- Mat-Select2 heap-primitive wins over linear saddleback count: `0`
- Mean Mat-Select2 heap-primitive/linear internal time ratio: `10.362578266023231`
- Mat-Select2 soft-selector wins over linear saddleback count: `0`
- Mean Mat-Select2 soft-selector/linear internal time ratio: `109.8413096886784`
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
| `2,3,5,7,11` | 0.000531 | 0.000461 | 0.868173 | 0.000835 | 1.572505 | 0.005258 | 9.902072 | 0.053879 | 101.467043 | 0 | 0.139881 | 32504 | 1000000 | 0.053059 |
| `2,3,5,7,11,13` | 0.000858 | 0.000791 | 0.921911 | 0.001231 | 1.434732 | 0.010887 | 12.688811 | 0.110422 | 128.696970 | 0 | 0.181112 | 33216 | 1000000 | 0.032455 |
| `2,3,5,7,11,13,17,19` | 0.001906 | 0.001057 | 0.554565 | 0.002557 | 1.341553 | 0.016195 | 8.496852 | 0.189380 | 99.359916 | 0 | 0.283123 | 34192 | 1000000 | 0.047309 |
