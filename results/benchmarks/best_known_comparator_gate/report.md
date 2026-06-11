# Best-Known Comparator Gate

- Timestamp: `2026-06-11T21:54:45.240021+00:00`
- Git commit: `d3e2f716a92ef6e8270bb5b1bc346205f69f77f1`
- Git dirty: `False`

## Gate Summary

- Mirzaian-Arjomandi full-unrank gate: `passed`
- Output-sensitive X+Y/LOH probe gate: `executed`
- Full Frederickson-Johnson gate: `open_not_implemented`
- Soft-heap X+Y gate: `open_not_implemented`
- Barvinok-style external count gate: `tool_available_smoke_failed_exports_ready`

## Published Sorted-Matrix Selector

- Artifact: `results/benchmarks/best_known_comparator_gate/ma_full_unrank_first_k_1e12`
- Cases: `3`
- Certified same-exponent cases: `3`
- MA wall-time wins: `0`
- Mean MA/current wall ratio: `3.270836895451455`

## Output-Sensitive X+Y Probe

- Artifact: `results/benchmarks/best_known_comparator_gate/output_sensitive_xplusy_1e6`
- Cases: `3`
- Completed cases: `3`
- Mean block/linear ratio: `0.8458775068607555`
- Mean MA/linear ratio: `2.2590266549340225`

This is an output-sized probe at a feasible rank. It does not establish a
full-rank random-access comparison when `N_probe` is much smaller than the
paper target rank.

## Barvinok-Style Lattice Counting

- Available local tools: `{'barvinok_count': '/usr/local/Caskroom/miniforge/base/envs/smooth-lattice-count/bin/barvinok_count', 'normaliz': '/usr/local/Caskroom/miniforge/base/envs/smooth-lattice-count/bin/normaliz'}`
- Export directory: `results/benchmarks/best_known_comparator_gate/barvinok_inputs`
- Exported cases: `3`
- External count smoke statuses: `{'barvinok_count_1d_isl': -11, 'barvinok_count_2d_isl': -11, 'normaliz_polytope_number_lattice_points': 1}`

The exported `.ine` files are rational simplex inputs for external
Barvinok/LattE/Normaliz-style tools. The gate records tool availability and
tiny smoke commands separately; a failed smoke is treated as a blocked
external-tool comparison, not as mathematical evidence.

## Bottom Line

The repository currently has a clean, certified comparison against a
published Mirzaian-Arjomandi sorted-matrix selector wrapper. Full
Frederickson-Johnson, an actual soft-heap implementation, and an external
Barvinok-family count run remain open obligations.
