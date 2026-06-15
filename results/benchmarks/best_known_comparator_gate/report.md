# Best-Known Comparator Gate

- Timestamp: `2026-06-15T22:50:57.277149+00:00`
- Git commit: `aaef2d09004798d1719614567c71f6f325f0cfe6`
- Git dirty: `False`

## Gate Summary

- Mirzaian-Arjomandi full-unrank gate: `passed`
- Output-sensitive X+Y/LOH probe gate: `executed`
- Full Frederickson-Johnson gate: `open_not_implemented`
- Soft-heap X+Y gate: `semantic_probe_validated_not_selector_integrated`
- Barvinok-style external count gate: `pynormaliz_toy_passed_target_counts_timed_out`

## Published Sorted-Matrix Selector

- Artifact: `results/benchmarks/best_known_comparator_gate/ma_full_unrank_first_k_1e12`
- Cases: `3`
- Certified same-exponent cases: `3`
- MA wall-time wins: `0`
- Mean MA/current wall ratio: `3.6952038158297444`

## Output-Sensitive X+Y Probe

- Artifact: `results/benchmarks/best_known_comparator_gate/output_sensitive_xplusy_1e6`
- Cases: `3`
- Completed cases: `3`
- Mean block/linear ratio: `0.9531453682122057`
- Mean MA/linear ratio: `2.848539503864853`
- Mat-Select2 heap-primitive comparable cases: `3`
- Mat-Select2 heap-primitive wins: `0`
- Mean Mat-Select2 heap-primitive/linear ratio: `9.701120596922914`

This is an output-sized probe at a feasible rank. It does not establish a
full-rank random-access comparison when `N_probe` is much smaller than the
paper target rank. The Mat-Select2 heap-primitive row follows the
Kaplan/Frederickson-Johnson exponential-block selector with an exact binary
heap primitive; it is not a soft-heap time-bound implementation.

## Soft-Heap Probe

- Status: `semantic_probe_validated_not_selector_integrated`
- Validation rows: `3`
- Timing row: `{'n': 20000, 'epsilon': 0.25, 'soft_sec': 0.024353, 'binary_heap_sec': 0.002724, 'corrupt_returns': 19348, 'soft_over_binary_heap': 8.940161527165932}`

The soft-heap row is a data-structure semantics probe. It checks the
corruption-set/witness-set invariants and simultaneous corruption bound for
a soft sequence heap, then records a small timing comparison against a
binary heap. It is not yet the row-sorted or `X+Y` selector from the
Kaplan/Frederickson-Johnson paper.

## Barvinok-Style Lattice Counting

- Available local tools: `{'barvinok_count': '/usr/local/Caskroom/miniforge/base/envs/smooth-lattice-count/bin/barvinok_count', 'normaliz': '/usr/local/Caskroom/miniforge/base/envs/smooth-lattice-count/bin/normaliz'}`
- Export directory: `results/benchmarks/best_known_comparator_gate/barvinok_inputs`
- Exported cases: `3`
- External count smoke statuses: `{'barvinok_count_1d_isl': -11, 'barvinok_count_2d_isl': -11, 'normaliz_polytope_number_lattice_points': 1}`
- PyNormaliz probe artifact: `results/benchmarks/best_known_comparator_gate/normaliz_count_probe`
- PyNormaliz toy cases passed: `2/2`
- PyNormaliz target counts completed: `0/3`
- PyNormaliz target timeouts: `3`

The exported `.ine` files and PyNormaliz rational-vertex inputs are
external Barvinok/LattE/Normaliz-style comparator inputs. The target
PyNormaliz counts use rationalized log simplexes; they are performance
comparisons, not correctness certificates for the original irrational-log
rank problem.

## Bottom Line

The repository currently has a clean, certified comparison against a
published Mirzaian-Arjomandi sorted-matrix selector wrapper. The
Mat-Select2 heap-primitive bridge is implemented and negative on the
current probe. A soft-sequence-heap semantic prototype now exists, but
a fast selector-integrated soft heap remains an open obligation.
The external Normaliz path has executable toy-count validation and bounded
certified-target attempts.
