# Next Implementation Milestones

## Milestone 1: Remove uncertified floating dependencies from core published solvers
- Replace floating-log search brackets with interval-certified brackets.
- Permit floating estimates only to propose candidate intervals, never as proof.

## Milestone 2: Arbitrary-prime log-bound generator
- Generate rational intervals for log p for any input prime p.
- Store/generated bounds with proof of error terms.
- Integrate into C++ auditor.

## Milestone 3: General k<=12 auditor
- Extend split auditor from k<=8 to k<=12 for the supported first-prime universe.
- Use MITM interval pairs with exact ambiguity resolution.
- Test first 8, first 10, and first 12 primes at N=10^12 under one shared correction policy.

## Milestone 4: Reproducible benchmark suite
- Benchmarks: DP pointer, heap frontier, materialized MITM, adaptive X+Y, layer-compressed, sums-only, parent-pointer.
- Report wall time, peak RSS, result vector, certification status.
- Include raw logs and CSV outputs.

## Milestone 5: Comparator obligation closure
- Decide whether the paper will implement a faithful Frederickson-Johnson or
  soft-heap `X+Y` comparator. If yes, satisfy the acceptance checklist in
  `paper/related_work.md`.
- If no, say so explicitly in the paper: current implemented comparators are DP,
  heap/frontier where present, adaptive materialized `X+Y`, full materialized
  `X+Y` unrank, and diagnostic sorted-matrix/LOH workbench probes.
- Keep the negative sorted-matrix/range-pruning artifact in the benchmark
  section rather than silently dropping it.

## Milestone 6: Paper draft
- Submit after theorem statements, certified k<=12 benchmarks, and a clear comparison section are complete.
