# Baselines and comparator families

This project targets fixed-prime random access, not general smooth-number counting.

## Included baselines

- `benchmarks/dp_pointer_baseline.cpp`: standard DP/pointer method used for Hamming/super-ugly-number generation. It is a strong sequential-generation baseline but has O(Nk) time and O(N) memory if it must return rank N.
- `benchmarks/heap_frontier_baseline.cpp`: canonical duplicate-free heap/frontier generation. It returns the full exponent vector by expanding children in nondecreasing prime-index order, so each exponent vector has one frontier path. It is a serious sequential-generation comparator, not a random-access method.
- `benchmarks/smooth_xplusy_baseline.cpp`: practical adaptive Cartesian-sum value-selection baseline.
- `benchmarks/smooth_xplusy_full_unrank.cpp`: practical materialized Cartesian-sum full-unrank baseline. It stores exponent packs on both MITM sides, narrows to a log-value band, exact-sorts the candidate band by multiprecision integer value, and returns an exponent vector. Its modes include adaptive X+Y, residual-corrected analytic X+Y, and a Mirzaian-Arjomandi sorted-matrix selector wrapped with the same exact reconstruction step.
- `benchmarks/smooth_xplusy_fj_loh_workbench.cpp`: exploratory sorted-matrix/range-pruning, Mirzaian-Arjomandi value-selection, Kaplan/Frederickson-Johnson-style Mat-Select2 heap-primitive, and LOH-style probes. Use `scripts/run_sorted_matrix_workbench.py` for reproducible artifacts; this is not a full soft-heap implementation.

## Current priority comparison

The next narrow target is the first-five-primes case:

```text
P = (2,3,5,7,11)
N = 10^12
```

Use `python3 scripts/run_xplusy_vs_layer5.py` to compare the practical
adaptive `X+Y` value-selection baseline against the layer-compressed full
unrank solver. The harness records command lines, compiler metadata, wall time,
peak RSS when available, raw stdout/stderr, and interval-audit status.

The `X+Y` row is selection-only, not full exponent reconstruction. A win by the
layer-compressed solver is therefore a conservative result; a loss should be
reported as such.

The stricter follow-up target is the full materialized `X+Y` unrank baseline:

```bash
python3 scripts/run_full_xplusy_suite.py
```

This compares two full exponent-vector unrankers across the same six five-prime
suite. The harness audits both returned vectors independently and requires the
vectors to match before counting a case as certified comparison evidence.

Current clean artifact:

```text
results/benchmarks/xplusy_vs_layer5_1e12/
results/benchmarks/five_prime_suite_1e12/
results/benchmarks/full_xplusy_suite_1e12/
results/benchmarks/sorted_matrix_workbench_1e12/
results/benchmarks/ma_full_unrank_suite_1e12/
results/benchmarks/ma_full_unrank_first_k_1e12/
results/benchmarks/best_known_comparator_gate/
results/benchmarks/heap_frontier_baseline_suite_1e5_1e6/
results/benchmarks/heap_frontier_baseline_suite_1e7/
```

On the recorded macOS/x86_64 Apple-clang run at commit
`161445507617a9435f9baadf4e70a3679d9e8d9a`, layer-compressed full unrank took
1.194730s wall time versus 2.029418s for adaptive materialized `X+Y` value
selection, with the returned exponent vector independently certified.

On the broader six-case suite over all five-prime subsets of
`{2,3,5,7,11,13}` at `N=10^12`, layer-compressed full unrank won all six
wall-time comparisons against the same adaptive materialized `X+Y` value
selection baseline, and all six returned exponent vectors were independently
certified.

On the stricter six-case suite against full materialized `X+Y` unranking,
recorded at commit `5314fd8a56fb7ca760046076e1f7f168fd48386a`,
layer-compressed full unrank won all six wall-time comparisons. Both methods
returned the same exponent vector in every case, and both outputs were
independently certified. The full-X+Y/layer wall-time ratio ranged from 1.9672
to 2.7248, with mean 2.3273. In this artifact the layer solver used the
analytic asymptotic bracket in all six cases; this is an optimization seed, not
a certificate.

Use `python3 scripts/run_sorted_matrix_workbench.py` to record exploratory
sorted-matrix/range-pruning, Mirzaian-Arjomandi value-selection,
Mat-Select2 heap-primitive, and LOH probes. These rows are useful negative or
diagnostic evidence, but they do not discharge the "full soft-heap X+Y"
comparison obligation.

In the clean six-case `N=10^12` workbench artifact at commit
`17a7e40ae790ae9017ab8274204eaaf107add212`, the range-pruned block counter
beat the ordinary linear saddleback count in only 1/6 cases, and only narrowly.
The mean block/linear internal time ratio was 1.3467, so this probe does not
improve the practical X+Y comparator. The Mirzaian-Arjomandi value selector
passed exhaustive small validation on 5185/5185 cases, matched the adaptive
selected log in 6/6 large cases, but won 0/6 timing comparisons, with mean
MA/linear internal time ratio 9.4042. The LOH row is a capped
`N_probe=10^6` top-k probe, not a full-rank `N=10^12` random-access comparator.

Use `python3 scripts/run_ma_full_unrank_suite.py` for the full
Mirzaian-Arjomandi comparator. This is no longer merely a value-selection
probe: it selects the X+Y log value, reconstructs exponent vectors, exact-sorts
the candidate band, and audits matching vectors. In the six-case `N=10^12`
artifact, the MA full-unrank path matched the analytic-corrected X+Y output and
was independently certified in 6/6 cases, but won 0/6 wall-time comparisons.
The mean MA/corrected wall-time ratio was 5.7000, so this comparator is useful
negative evidence rather than a speed path.

The clean first-k MA full-unrank artifact extends this published sorted-matrix
comparator to the first 5, 6, and 8 primes at `N=10^12`:

```text
results/benchmarks/ma_full_unrank_first_k_1e12/
```

Observed at commit `7daf269404c1672a69da6f6d2edb2ff53eb648e3`, all 3 rows
completed, matched analytic-corrected X+Y exponent vectors, and were
independently interval-certified. MA won 0/3 wall-time comparisons; the mean
MA/corrected wall-time ratio was 3.7859. This is the current cleanest
best-known-style comparator checkpoint in the repository, but it is still
narrowly about Mirzaian-Arjomandi sorted-matrix value selection wrapped into
full unranking.

The Mat-Select2 heap-primitive artifact adds a direct
Kaplan/Frederickson-Johnson-style exponential-block row-sorted selector probe:

```text
results/benchmarks/sorted_matrix_matselect2_heap_firstk_1e12/
```

It passed 3233/3233 exhaustive small validation cases. On the first
k=5, N=10^12 target it matched the selected log exactly but was 13.0590x slower
than the current linear saddleback selector. The first k=6 and k=8 rows were
skipped by the active-row cap. This is useful negative evidence for the
exponential-block selector with an exact binary heap primitive, but it is still
not a full soft-heap time-bound implementation.

Use `python3 scripts/run_heap_frontier_baseline_suite.py` for the full
heap/frontier comparator. The heap baseline returns exponent vectors, and the
harness checks agreement with the DP pointer baseline where enabled and with the
current solver, then independently interval-audits the agreed vector.

Clean default artifact:

```text
results/benchmarks/heap_frontier_baseline_suite_1e5_1e6/
```

Observed result at commit `b26e3cdf865b15754d8b189dafe4e7cc36e19f80`:

- 8/8 rows over k=3,5,6,8 and N=10^5,10^6 were solver-agreement rows.
- 8/8 rows were independently interval-certified.
- Wall-time winners were current_beatty3 in 2 rows, current_sums_adaptive in 5
  rows, and dp_pointer in 1 row.
- Mean heap/current wall-time ratio was 4.5654; mean heap/DP wall-time ratio was
  1.9641.

Clean stress artifact:

```text
results/benchmarks/heap_frontier_baseline_suite_1e7/
```

Observed result at commit `6526d68a0adcd45b0041a1dcfe8da1c6a9fe1f56`:

- 4/4 rows over k=3,5,6,8 at N=10^7 were solver-agreement rows.
- 4/4 rows were independently interval-certified.
- Current solvers won all 4 wall-time comparisons.
- Mean heap/current wall-time ratio was 76.5315; mean heap/DP wall-time ratio
  was 2.9536.

This supports a narrow claim against a canonical sequential heap/frontier
generator on the tested fixed-prime ranks. It does not imply a comparison
against soft-heap `X+Y`, Frederickson-Johnson sorted-matrix selection, or
Barvinok-style lattice counting.

Use `python3 scripts/run_best_known_comparator_gate.py` to produce the current
serious-comparator status dashboard. In the current artifact, the MA
full-unrank gate passed, the output-sensitive `X+Y` probe executed at N=10^6,
full FJ and soft-heap gates remained open/not implemented, and the
Barvinok/Normaliz external count path passed PyNormaliz toy simplex validation
but timed out on certified-target rationalized simplexes under the configured
cap.

## Not yet fully implemented

- Full soft-heap time-bound sorted matrix / X+Y selection.
- Full soft-heap X+Y selection implementations.
- Competitive target-scale Barvinok/Normaliz-style fixed-dimensional
  lattice-point counting.
- Bernstein-style general y-smooth counting/listing algorithms.

Claims should be limited accordingly.
