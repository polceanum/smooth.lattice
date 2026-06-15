# Best-Known Comparator Status

This project should use "best known" language only for a named, narrow
algorithmic subproblem with a reproducible artifact. The current strongest
implemented published comparator is the Mirzaian-Arjomandi sorted-matrix
selection algorithm wrapped into full fixed-prime unranking.

## Current Implemented Published Comparator

Comparator:

- A. Mirzaian and E. Arjomandi, "Selection in X+Y and matrices with sorted rows
  and columns", Information Processing Letters 20(1), 1985.
  https://doi.org/10.1016/0020-0190(85)90123-1

Why it is relevant:

- A fixed-prime MITM split produces sorted lists `A` and `B`.
- Selecting the N-th Cartesian sum in `A+B` is the sorted-matrix selection
  subproblem.
- The implemented `ma-full` mode selects the log value with the
  Mirzaian-Arjomandi selector, reconstructs a narrow exponent band, exact-sorts
  the band by multiprecision integer value, and independently audits the final
  exponent vector.

Validation:

- `bin/smooth_xplusy_fj_loh_workbench validate-ma`
- Latest local validation: 5185/5185 exhaustive small Cartesian-sum cases
  passed, maximum delta 0.

## Clean First-k Result

Artifact:

```text
results/benchmarks/ma_full_unrank_first_k_1e12/
```

Observed at commit `7daf269404c1672a69da6f6d2edb2ff53eb648e3`:

- P=(2,3,5,7,11), N=10^12: certified same vector, MA/current wall ratio 3.5898.
- P=(2,3,5,7,11,13), N=10^12: certified same vector, MA/current wall ratio 4.0593.
- P=(2,3,5,7,11,13,17,19), N=10^12: certified same vector, MA/current wall ratio 3.7086.
- Aggregate: 3/3 completed, 3/3 same-exponent certified, 0/3 MA wall-time wins,
  mean MA/current wall ratio 3.7859.

This supports a narrow statement:

> On first-prime fixed-prime random-access targets for k=5,6,8 at N=10^12, the
> current analytic-corrected X+Y unranker returns the same independently
> certified exponent vectors as a Mirzaian-Arjomandi sorted-matrix selection
> wrapper and is faster in all tested rows on the recorded machine.

## Still Open

This does not justify broad "best known" or "state of the art" language. The
open comparator obligations remain:

- full Frederickson-Johnson ranking/selection implementation;
- soft-heap X+Y/top-k implementation where the output-sensitive regime is
  appropriate;
- Barvinok-style fixed-dimensional lattice counting or a theoretical argument
  explaining why it is not the right practical comparator for these targets;
- arbitrary-prime independent certification beyond the current k<=12 fixed
  audited universe.

## Comparator Gate Artifact

Artifact:

```text
results/benchmarks/best_known_comparator_gate/
```

Observed at commit `d3e2f716a92ef6e8270bb5b1bc346205f69f77f1`:

- Mirzaian-Arjomandi full-unrank gate: passed.
- Output-sensitive `X+Y`/LOH probe gate at N=10^6: executed.
- Full Frederickson-Johnson gate: open, not implemented.
- Soft-heap `X+Y` gate: open, not implemented.
- Barvinok-style external count gate: PyNormaliz toy lattice-count validation
  passed, but the certified-target rationalized simplexes timed out under the
  bounded gate; no external-count speed claim is made.
- Kaplan/Frederickson-Johnson-style Mat-Select2 heap-primitive: implemented as
  an exact row-sorted selector probe. It passed exhaustive small validation, but
  on the first k=5, N=10^12 row it was 13.0590x slower than the current linear
  saddleback selector; k=6 and k=8 were skipped by active-row cap.

Gate details:

- MA first-k comparison: 3/3 completed and certified, 0/3 MA wall-time wins,
  mean MA/current wall ratio 3.5772.
- Output-sensitive probe: 3/3 completed, range-pruning beat linear saddleback
  in 3/3 rows at N=10^6 with mean block/linear ratio 0.8691; MA value selection
  won 0/3 with mean MA/linear ratio 1.8716.
- Barvinok/Normaliz availability: `barvinok_count`, `normaliz`, and
  `PyNormaliz` were installed in the `smooth-lattice-count` conda environment.
  `barvinok_count` still segfaulted on toy ISL sets, but the PyNormaliz
  rational-vertex path passed 2/2 known toy simplex counts. The three
  certified-target rationalized simplexes at one decimal digit all timed out at
  the gate's 5-second per-count cap.
- Mat-Select2 heap-primitive artifact:
  `results/benchmarks/sorted_matrix_matselect2_heap_firstk_1e12/`. Validation
  passed 3233/3233 small cases; the single comparable first-k target was exact
  but slower by a factor of 13.0590.

The gate is a status artifact, not a victory lap. It says exactly where the
repository stands today: one published sorted-matrix selector comparison is
cleanly beaten on certified rows; the Mat-Select2 heap-primitive path is now a
documented negative comparator; full soft-heap comparisons still need real
implementation, while the off-the-shelf Normaliz path is presently a documented
timeout rather than a competitive count oracle.
