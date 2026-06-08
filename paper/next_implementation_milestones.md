# Next Implementation Milestones

## Milestone 1: Remove uncertified floating dependencies from core published solvers
- Replace floating-log search brackets with interval-certified brackets.
- Permit floating estimates only to propose candidate intervals, never as proof.

## Milestone 2: Arbitrary-prime log-bound generator
- Generate rational intervals for log p for any input prime p.
- Store/generated bounds with proof of error terms.
- Integrate into C++ auditor.

## Milestone 3: General k<=8 auditor
- Extend split auditor from k<=6 to k<=8.
- Use MITM interval pairs with exact ambiguity resolution.
- Test first 8 primes at N=10^12.

## Milestone 4: Reproducible benchmark suite
- Benchmarks: DP pointer, heap frontier, materialized MITM, adaptive X+Y, layer-compressed, sums-only, parent-pointer.
- Report wall time, peak RSS, result vector, certification status.
- Include raw logs and CSV outputs.

## Milestone 5: Paper draft
- Submit after theorem statements, certified k<=6/8 benchmarks, and a clear comparison section are complete.
