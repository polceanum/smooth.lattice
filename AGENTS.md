# AGENTS.md

Guidance for Codex or other coding agents working in this repository.

## Project goal

This is a research prototype for fixed-prime smooth-number random access. The main task is:

```text
Given a finite prime set P and 1-based rank N, return the exponent vector of the N-th P-smooth number.
```

Do not reframe this as prime generation, factoring, or general smooth-number counting unless explicitly asked.

## Core architecture

- `smooth_toolkit.py` is the user-facing CLI/API wrapper.
- `src/smooth_3prime_beatty_ranker.cpp` is the special 3-prime kernel.
- `src/smooth_layer_compressed_general.cpp` is the 5-prime layer-compressed kernel.
- `src/smooth_sums_only_scalable.cpp` is the general/higher-k adaptive MITM fallback.
- `src/smooth_interval_audit_exps_k6.cpp` is the independent proof/audit kernel for `k <= 6` and primes drawn from `{2,3,5,7,11,13}`.

## Build and test commands

Use these commands before claiming a change works:

```bash
./build.sh
python3 tests/test_smoke.py
```

For certified examples:

```bash
python3 smooth_toolkit.py audit 2,3,5,7,11 1000000000000 --timeout 300
python3 smooth_toolkit.py audit 2,3,5,7,11,13 1000000000000 --timeout 300
```

The second audit may take several seconds and use hundreds of MB of RAM.

## Correctness rules

- Do not claim a fast result is certified unless the independent auditor returns `rank_certified=true`.
- Do not broaden the certified range without adding rigorous log bounds and passing an audit test.
- Avoid increasing floating-point tolerances to hide failures.
- If final-band selection changes, add an exact-integer validation path or a small-N fuzz test against DP generation.
- Preserve 1-based indexing: `N=1` corresponds to exponent vector all zeroes.

## Benchmark rules

- Keep raw benchmark output under `results/benchmarks/`.
- Include compiler flags, CPU/RAM if available, command used, and method name.
- Separate tasks clearly:
  - sequential generation of first `N` terms;
  - fixed-prime random access / unrank;
  - count/rank audit;
  - general `y`-smooth counting.
- Do not claim “state of the art” unless comparing against a serious published implementation or a faithful implementation of the relevant algorithm family.

## Code style

- Current C++ kernels are single-file prototypes. Prefer minimal, reviewable changes.
- Keep all core kernels C++17-compatible.
- Do not add mandatory external dependencies beyond standard C++17, Boost.Multiprecision headers, and Python 3 unless discussed.
- Python wrapper should remain dependency-free.

## Known limitations

- Fast high-k kernels still use floating log sums for search.
- The independent auditor currently supports only `k <= 6`, primes subset `{2,3,5,7,11,13}`.
- Results beyond audited cases should be labelled experimental.

## Suggested next tasks

1. Extend rigorous log-bound generation beyond prime 13.
2. Add a certified audit path for `k=8`.
3. Add fuzz tests comparing small `N` against exact DP generation for varied prime sets.
4. Refactor C++ kernels into shared library components only after the algorithms stabilize.
5. Prepare a reproducible paper artifact with raw logs and exact commands.
