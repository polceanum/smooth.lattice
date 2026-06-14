# Best-Known Comparator Gate

- Timestamp: `2026-06-14T22:34:05.412709+00:00`
- Git commit: `26366bf06aa59d88f29d3d3264578eda7992b045`
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
- Mean MA/current wall ratio: `3.5771829880074257`

## Output-Sensitive X+Y Probe

- Artifact: `results/benchmarks/best_known_comparator_gate/output_sensitive_xplusy_1e6`
- Cases: `3`
- Completed cases: `3`
- Mean block/linear ratio: `0.8690779105459225`
- Mean MA/linear ratio: `1.8715664428483259`

This is an output-sized probe at a feasible rank. It does not establish a
full-rank random-access comparison when `N_probe` is much smaller than the
paper target rank.

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
published Mirzaian-Arjomandi sorted-matrix selector wrapper. Full
Frederickson-Johnson and an actual soft-heap implementation remain open
obligations. The external Normaliz path now has executable toy-count
validation and bounded certified-target attempts.
