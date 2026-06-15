# Best-Known Comparator Gate

- Timestamp: `2026-06-15T22:22:47.786221+00:00`
- Git commit: `3c852e9b88e955b6e4e8de6d55b9c90a2476c703`
- Git dirty: `False`

## Gate Summary

- Mirzaian-Arjomandi full-unrank gate: `passed`
- Output-sensitive X+Y/LOH probe gate: `executed`
- Full Frederickson-Johnson gate: `open_not_implemented`
- Soft-heap X+Y gate: `open_not_implemented`
- Barvinok-style external count gate: `pynormaliz_toy_passed_target_counts_timed_out`

## Published Sorted-Matrix Selector

- Artifact: `results/benchmarks/best_known_comparator_gate/ma_full_unrank_first_k_1e12`
- Cases: `3`
- Certified same-exponent cases: `3`
- MA wall-time wins: `0`
- Mean MA/current wall ratio: `3.8331755528502947`

## Output-Sensitive X+Y Probe

- Artifact: `results/benchmarks/best_known_comparator_gate/output_sensitive_xplusy_1e6`
- Cases: `3`
- Completed cases: `3`
- Mean block/linear ratio: `0.9566811987145698`
- Mean MA/linear ratio: `2.007773448900695`
- Mat-Select2 heap-primitive comparable cases: `3`
- Mat-Select2 heap-primitive wins: `0`
- Mean Mat-Select2 heap-primitive/linear ratio: `10.167390045375946`

This is an output-sized probe at a feasible rank. It does not establish a
full-rank random-access comparison when `N_probe` is much smaller than the
paper target rank. The Mat-Select2 heap-primitive row follows the
Kaplan/Frederickson-Johnson exponential-block selector with an exact binary
heap primitive; it is not a soft-heap time-bound implementation.

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
current probe. A true soft-heap implementation remains an open obligation.
The external Normaliz path has executable toy-count validation and bounded
certified-target attempts.
