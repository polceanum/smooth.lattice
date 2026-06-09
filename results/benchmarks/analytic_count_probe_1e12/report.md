# Analytic Lattice Count Probe

- Timestamp: `2026-06-09T14:08:51.981889+00:00`
- Git commit: `f244aeafa0235adbb5f9441d4ff8c01b84e938e4`
- Git dirty: `False`
- Cases: `8`
- Completed cases: `8`
- Layer-count matches expected rank: `6`
- Max |expected - analytic|: `1979764.706268`
- Mean |expected - analytic|: `247583.15563725`
- Max absolute relative residual: `2e-06`

These probes evaluate the asymptotic analytic count approximation at
previously certified exponent vectors. `expected_N` is the certified rank
for the vector; `layer_count` is the current floating-log layer counter at
the same log height. This is diagnostic evidence only, not a proof of an
analytic error bound.

| label | k | residual expected-analytic | relative residual | derivative | layer-expected | wall s | A | Base |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `k5_P2_3_5_7_11` | 5 | -71.083878 | -0.000000000 | 5391076448.621829 | 0 | 1.144917 | 91837 | 241925 |
| `k5_P2_3_5_7_13` | 5 | 61.185536 | 0.000000000 | 5318948810.991179 | 0 | 0.202394 | 88206 | 248507 |
| `k5_P2_3_5_11_13` | 5 | -20.278783 | -0.000000000 | 5101341257.229017 | -1 | 0.253751 | 77825 | 270095 |
| `k5_P2_3_7_11_13` | 5 | 35.163448 | 0.000000000 | 4911279582.834452 | 0 | 0.197780 | 83950 | 241069 |
| `k5_P2_5_7_11_13` | 5 | -166.388496 | -0.000000000 | 4550180340.370217 | 0 | 0.160959 | 97784 | 191789 |
| `k5_P3_5_7_11_13` | 5 | -93.612146 | -0.000000000 | 4149775994.535069 | -1 | 0.172085 | 117565 | 230611 |
| `k6_first6` | 6 | 452.826543 | 0.000000000 | 12809777586.262388 | 0 | 1.430192 | 61042 | 1415628 |
| `k8_first8` | 8 | 1979764.706268 | 0.000002000 | 36858875270.849808 | 0 | 1.599638 | 457792 | 1717623 |
