# Related Work And Comparator Positioning

This note is for paper preparation. It separates comparator families that are
implemented and benchmarked in this repository from literature baselines that
must be discussed but are not yet implemented faithfully.

## Fixed-Prime Random Access Versus General Smooth Counting

The target problem here is fixed-prime random access:

```text
unrank(P, N) -> exponent vector e
```

where `P` is fixed and finite and `N` is a rank in the ordered sequence of
`P`-smooth numbers. This is not the same problem as estimating or computing the
general smooth-counting function `Psi(x, y)`, and it is not the same problem as
finding smooth parts of arbitrary input integers.

Bernstein's smooth-number work belongs in the paper as important adjacent
number-theoretic algorithmics, but the problem statement is different. His
"smooth parts" algorithm takes a finite prime set and a sequence of integers and
finds the largest `P`-smooth divisor of each integer. That is a batch
recognition/factor-extraction problem, not fixed-prime rank/unrank.

Source to cite:

- D. J. Bernstein, "How to find smooth parts of integers", 2004.
  https://cr.yp.to/factorization.html

## Sequential Generation Baselines

The standard super-ugly-number dynamic program is a serious common baseline for
generating the first `N` fixed-prime smooth numbers. It keeps one pointer per
prime and advances the minimum candidate. For random access to rank `N`, this is
sequential generation: `O(Nk)` time and `O(N)` stored values if the final value
must be recovered without overflow.

Safe comparison language:

```text
For fixed-prime random access at large N, the certified rank/unrank solver is
faster than the standard O(Nk) pointer-generation baseline under the recorded
conditions.
```

Do not imply that the method dominates DP generation when the task is to list all
prefix values.

Implemented comparator:

- `benchmarks/dp_pointer_baseline.cpp`

## Cartesian-Sum Selection: X+Y

Meet-in-the-middle rank/unrank for fixed primes naturally produces Cartesian
sums. If `A` and `B` are sorted lists of partial log weights, then rank counting
or value selection asks about the sorted multiset:

```text
A + B = {a + b : a in A, b in B}.
```

Classical selection in `X+Y` and sorted matrices is therefore the most important
algorithmic baseline family after sequential generation.

### Johnson-Mizoguchi

Johnson and Mizoguchi gave an `O(n log n)` algorithm for selecting the `K`-th
element of `X+Y`, with extensions to sums of more than two sets. This is
theoretically relevant because it treats Cartesian-sum selection directly.

Source to cite:

- D. B. Johnson and T. Mizoguchi, "Selecting the Kth Element in X+Y and
  X1+X2+...+Xm", SIAM Journal on Computing 7(2), 1978.
  https://doi.org/10.1137/0207013

Current implementation status:

- Not implemented faithfully.
- The materialized `X+Y` comparators in this repository use adaptive rank counts
  plus final-band selection instead.

### Mirzaian-Arjomandi

Mirzaian and Arjomandi gave an `O(n)` selection algorithm for an `n x n` matrix
with sorted rows and columns, and they note that sorted `X+Y` forms such a
matrix. This is a serious value-selection comparator for square Cartesian sums.

Source to cite:

- A. Mirzaian and E. Arjomandi, "Selection in X+Y and matrices with sorted rows
  and columns", Information Processing Letters 20(1), 1985.
  https://doi.org/10.1016/0020-0190(85)90123-1

Current implementation status:

- Implemented as `ma_select_probe` in
  `benchmarks/smooth_xplusy_fj_loh_workbench.cpp`.
- Rectangular MITM splits are embedded in a square matrix by padding the shorter
  side with negative infinity after transforming kth-smallest log selection into
  kth-largest selection on negated weights.
- This row selects a log value only. It does not reconstruct an exponent vector
  or independently certify a returned smooth number.
- Clean `N=10^12` five-prime suite artifact:
  `results/benchmarks/sorted_matrix_workbench_1e12/`. The probe matched the
  adaptive selected log in 6/6 cases but was slower in 6/6 cases, with mean
  MA/linear internal time ratio 9.0509.

### Frederickson-Johnson

Frederickson and Johnson studied selection and ranking in `X+Y` and matrices
with sorted columns, giving algorithms whose complexity depends on rank and is
sublinear in the full Cartesian product size for relevant regimes. They also
developed generalized selection/ranking for sorted matrices.

Sources to cite:

- G. N. Frederickson and D. B. Johnson, "The Complexity of Selection and Ranking
  in X+Y and Matrices with Sorted Columns", Journal of Computer and System
  Sciences 24(2), 1982.
  https://doi.org/10.1016/0022-0000(82)90048-4
- G. N. Frederickson and D. B. Johnson, "Generalized Selection and Ranking:
  Sorted Matrices", SIAM Journal on Computing 13(1), 1984.
  https://doi.org/10.1137/0213002

Current implementation status:

- Not implemented faithfully.
- `benchmarks/smooth_xplusy_fj_loh_workbench.cpp` contains a range-pruned
  sorted-matrix counting probe. It is FJ-style diagnostic code, not the published
  FJ algorithm.
- The clean workbench result at
  `results/benchmarks/sorted_matrix_workbench_1e12/` is negative/mixed: the
  range-pruned block counter won only 2/6 cases, narrowly, and was slower on
  average.

Acceptance criteria before claiming an FJ comparison:

- Implement the actual selection/ranking algorithm from the paper, including its
  reduction phases, not only a quadtree/block pruning heuristic.
- Validate on small Cartesian sums against exhaustive sorting and against the
  existing adaptive `X+Y` count.
- Measure wall time, memory, and result equality on the six-case five-prime
  suite.
- State clearly whether the implementation returns only a selected log value or
  a full exponent vector.

### Soft-Heap X+Y

Kaplan, Kozma, Zamir, and Zwick give simpler optimal algorithms for selection
from heaps, sorted lists, row-sorted matrices, and `X+Y` using soft heaps. This
is a serious modern theoretical comparator, especially for output-sensitive
selection of the `k` smallest items.

Source to cite:

- H. Kaplan, L. Kozma, O. Zamir, and U. Zwick, "Selection from Heaps,
  Row-Sorted Matrices, and X+Y Using Soft Heaps", SOSA 2019.
  https://doi.org/10.4230/OASIcs.SOSA.2019.5

Current implementation status:

- Not implemented.
- The current LOH probe is an output-style top-k diagnostic capped at
  `N_probe=10^6`; it is not a full-rank `N=10^12` random-access comparator.

Acceptance criteria before claiming a soft-heap comparison:

- Implement or vendor a soft heap with documented error/corruption semantics.
- Reproduce the paper's selection primitive on ordinary arrays/lists.
- Adapt it to the `X+Y` split used by the smooth-number MITM baseline.
- Record whether it returns only a selected value, the set of `k` smallest
  values, or a full exponent vector.

## Fixed-Dimensional Lattice Counting

For fixed prime set size `k`, the count

```text
#{e in N^k : sum_i e_i log(p_i) <= T}
```

is a fixed-dimensional lattice-point counting problem under one linear
inequality with irrational coefficients. Barvinok-style algorithms are therefore
important theoretical context because they give polynomial-time lattice-point
counting in fixed dimension for rational polyhedra.

Source to cite:

- A. I. Barvinok, "A Polynomial Time Algorithm for Counting Integral Points in
  Polyhedra When the Dimension is Fixed", Mathematics of Operations Research
  19(4), 1994.
  https://doi.org/10.1287/moor.19.4.769

Current implementation status:

- Not implemented.
- Direct application requires rational interval treatment of logarithms or a
  proof that rationalized bounds preserve the target rank decision.
- The independent interval auditor is currently the repository's proof-oriented
  substitute for headline examples, limited to `k <= 6` and primes in
  `{2,3,5,7,11,13}`.

## Implemented Baseline Set For Current Claims

Current paper claims may compare against:

- DP pointer sequential generation;
- heap/frontier generation where present;
- adaptive materialized `X+Y` value selection;
- Mirzaian-Arjomandi sorted-matrix value selection;
- full materialized `X+Y` unranking;
- the sorted-matrix/range-pruning workbench as negative diagnostic evidence.

Current paper claims may not say:

- "beats Frederickson-Johnson";
- "beats soft-heap X+Y";
- "beats Barvinok-style fixed-dimensional counting";
- "state of the art";
- "best known";

unless those comparator obligations are implemented, benchmarked, and documented
with the same artifact discipline as the current X+Y suites.
