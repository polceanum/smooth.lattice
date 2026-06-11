# Repository Audit

Date: 2026-06-08

Repository: `polceanum/smooth.lattice`

Audited commit: `22512d6599b776cad90236def34d7b0efb88942d`

Worktree status at audit time: dirty. `tests/test_smoke.py` contains an uncommitted small-N fuzz-test patch.

## Post-Audit Updates

This file records the original Phase 1 audit context. Subsequent hardening has
changed several items that were listed as gaps:

- the interval auditor now supports `k <= 12` for primes in
  `{2,3,5,7,11,13,17,19,23,29,31,37}`;
- `scripts/run_certified_portfolio_suite.py` records a certified portfolio over
  implemented unrankers;
- `benchmarks/heap_frontier_baseline.cpp` and
  `scripts/run_heap_frontier_baseline_suite.py` add a canonical full-vector
  heap/frontier comparator with clean certified artifacts;
- CI and local smoke tests were repaired after the initial Boost dependency
  failure.

For current paper claims, prefer `paper/claims_matrix.md`,
`paper/benchmark_protocol.md`, and `paper/hard_results_table.md` over the
historical limitations below.

## Scope

This audit followed the requested reading order:

1. `AGENTS.md`
2. `README.md`
3. `docs/`
4. `paper/`
5. `results/`
6. `src/`
7. `benchmarks/`
8. `tests/`

The repository is correctly framed around fixed-prime smooth-number random access:

```text
unrank(P, N) -> exponent vector
```

It mostly avoids broader claims about prime generation, factoring, and general `Psi(x,y)` smooth-number counting.

## Local Environment

Observed local toolchain:

```text
Python: Python 3.13.3
C++ compiler path: /usr/bin/g++
C++ compiler: Apple clang version 21.0.0 (clang-2100.0.123.102)
System: Darwin MacBookPro.Home 25.4.0 x86_64
```

CPU and memory metadata commands were attempted with:

```bash
sysctl -n machdep.cpu.brand_string
sysctl -n hw.memsize
```

Both failed in the current sandbox with:

```text
sysctl: sysctl fmt -1 1024 1: Operation not permitted
```

## Build And Test Results

Initial command:

```bash
./build.sh
```

Result: failed.

Failure:

```text
src/smooth_3prime_beatty_ranker.cpp:13:10: fatal error: 'boost/multiprecision/cpp_int.hpp' file not found
```

After the health patch in this working tree, the same command fails earlier with an explicit dependency diagnostic:

```text
ERROR: Boost.Multiprecision headers were not found by the active C++ compiler.
```

Initial command:

```bash
python3 tests/test_smoke.py
```

Result: failed at the same build step, before tests could exercise the algorithms.

Failure:

```text
src/smooth_3prime_beatty_ranker.cpp:13:10: fatal error: 'boost/multiprecision/cpp_int.hpp' file not found
subprocess.CalledProcessError: Command '['/Users/mike/Work/smooth.lattice/build.sh']' returned non-zero exit status 1.
```

After the health patch, this command also stops at the explicit Boost preflight. Interpretation: no local correctness claim can be made from this audit run. The immediate blocker is missing Boost.Multiprecision headers for the active compiler environment.

## Current Code Wiring

User-facing entry point:

- `smooth_toolkit.py`

Core binaries built by `build.sh`:

- `src/smooth_3prime_beatty_ranker.cpp` -> `bin/smooth_3prime_beatty_ranker`
- `src/smooth_layer_compressed_general.cpp` -> `bin/smooth_layer_compressed_general`
- `src/smooth_sums_only_scalable.cpp` -> `bin/smooth_sums_only_scalable`
- `src/smooth_interval_audit_exps_k6.cpp` -> `bin/smooth_interval_audit_exps_k6`

Benchmark binaries built by scripts:

- `benchmarks/dp_pointer_baseline.cpp`
- `benchmarks/smooth_xplusy_baseline.cpp`
- `benchmarks/smooth_xplusy_fj_loh_workbench.cpp`

No source file in `src/` appears completely orphaned; all four are referenced by `build.sh` and `smooth_toolkit.py`. Benchmark files are referenced by scripts and docs.

The manifest mentions `.github/workflows/ci.yml`. It was not present at the start of the audit; this working tree adds a CI workflow that installs `g++` and `libboost-dev`, then runs the normal build and smoke commands.

## Algorithm Notes

### `smooth_3prime_beatty_ranker.cpp`

Implements a special 3-prime ranker using a two-prime Beatty/floor-sum inner count. It uses `long double` logs and small epsilons for search and counting, then sorts final-band candidates by exact `cpp_int` values.

Status: promising experimental fast path. Not proof-grade without an independent audit of the returned candidate.

### `smooth_layer_compressed_general.cpp`

Implements a layer-compressed MITM method that represents a 3-prime side as shifted 2-prime layers. It uses floating log sums and tolerances in count/band construction, with exact `cpp_int` sorting for the final candidate band.

Status: useful high-performance solver, especially for 5-prime cases. Treat as experimental unless its returned candidate is independently audited.

### `smooth_sums_only_scalable.cpp`

Implements higher-k sums-only MITM with delayed exponent reconstruction. Reconstruction uses floating-sum DFS tolerances, including a loosened fallback tolerance.

Status: useful exploratory/high-k solver. It should not be a paper-grade correctness path without independent audit or exact reconstruction validation.

### `smooth_interval_audit_exps_k6.cpp`

Implements an independent interval-rank auditor for `k <= 6` and primes in `{2,3,5,7,11,13}` using hard-coded rational-ish scaled log intervals and exact `cpp_int` ambiguity resolution.

Important limitation: it aborts if sorted `B` intervals overlap, because the robust overlap fallback is not implemented.

Status: strongest correctness component in the repository for supported prime sets. Needs broader log-bound generation and overlap handling before claims are generalized.

## Correctness Status

Supported certified path:

- Fast solver returns candidate exponent vector.
- `smooth_interval_audit_exps_k6` verifies `count_le(P, x) == N`.
- This supports `k <= 6` and primes drawn from `{2,3,5,7,11,13}`.

Experimental paths:

- 3-prime Beatty fast ranker before audit.
- 5-prime layer-compressed fast ranker before audit.
- Higher-k sums-only MITM, especially for `k > 6`.

Current smoke tests compare selected small ranks against a heap oracle and audit one small 6-prime case, but the test command did not run locally because the build failed.

An uncommitted test patch currently expands small-N heap-oracle coverage across varied prime sets, including `k=1`, non-prefix prime sets, and `k=7`.

## Benchmark And Result Provenance

Stored certified examples:

- `results/certified/unrank_5prime_1e12.json`
- `results/certified/audit_5prime_1e12.json`
- `results/certified/unrank_6prime_1e12.json`
- `results/certified/audit_6prime_1e12.json`

These include raw command output and certification fields, but they do not include full machine/compiler metadata.

Stored benchmark summaries:

- `results/benchmarks/beat_dp_results.txt`
- `results/benchmarks/xplusy_summary.csv`
- `results/benchmarks/smooth_xplusy_fj_loh_summary.csv`
- `results/benchmarks/smooth_layer_compressed_general_summary.csv`
- `results/benchmarks/workbench_summary.csv`

Credibility gaps:

- Some CSV rows reference absolute paths under `/mnt/data/smooth_best_known_workbench/...`, not paths in this repository.
- `workbench_summary.csv` includes methods such as `smooth_hybrid_mitm_interpolation_tight` and `smooth_parent_pointer_scalable`, but matching source files are not present.
- `scripts/run_dp_comparison.sh` currently regenerates only a small local comparison for `N=10^6` and `N=10^7`; it does not regenerate the larger curated benchmark table.
- Result files do not consistently include git commit, compiler flags, CPU/RAM, or raw stderr.
- There is no single benchmark harness covering DP, heap, generic MITM, adaptive `X+Y`, current best solver, and certification status.

Conclusion: stored benchmarks are useful research logs, but they are not yet a fully reproducible paper artifact.

## Experimental Or Unsafe Areas

- Floating tolerances in fast counting/search paths are acceptable for experiments but must not be used as proof.
- `smooth_sums_only_scalable.cpp` reconstructs exponents from floating sums, then computes exact value only after reconstruction.
- The interval auditor has hard-coded log bounds only up to prime 13.
- The interval auditor has no overlap fallback if interval widths overlap in `B`.
- `Pack` encodings assume 16-bit local exponents; overflow is checked, but any future larger-`N` claims should mention this implementation limit or lift it.
- Benchmark scripts and curated results are not yet generated from one reproducible command.

## Health Patch Added In This Working Tree

This working tree adds:

- a Boost.Multiprecision preflight in `build.sh`;
- `.github/workflows/ci.yml` for build and smoke-test coverage on Ubuntu;
- a reproducibility note documenting the Boost dependency;
- expanded deterministic small-N heap-oracle tests in `tests/test_smoke.py`.

## Recommended Next Patch

The smallest next patch that improves correctness credibility is:

1. Install/provide Boost locally and rerun `./build.sh` plus `python3 tests/test_smoke.py`.
2. Add property tests for exponent-stream self-similarity and zero-position face recursion.
3. Start a reproducible benchmark harness that records git commit, compiler flags, CPU/RAM, exact command, raw stdout/stderr, and certification status.
