# smooth.lattice

Fixed-prime smooth-number rank/unrank experiments based on recursive exponent lattices.

A `P`-smooth number for a fixed prime set `P = {p1, ..., pk}` is a number of the form

```text
p1^e1 * p2^e2 * ... * pk^ek,  ei >= 0.
```

This repository focuses on **random access**:

```text
Given P and a 1-based rank N, return the exponent vector of the N-th P-smooth number.
```

It is not intended to replace general `y`-smooth counting/listing algorithms, and it is not a prime-generation or factoring algorithm.

## Current contents

Core kernels:

| File | Purpose |
|---|---|
| `src/smooth_3prime_beatty_ranker.cpp` | Fast special ranker for 3-prime fixed sets using a two-prime floor-sum/Beatty-style inner count. |
| `src/smooth_layer_compressed_general.cpp` | Layer-compressed hybrid MITM kernel, strongest for 5-prime cases in our tests. |
| `src/smooth_sums_only_scalable.cpp` | Sums-only/adaptive MITM fallback for other and higher prime counts. |
| `src/smooth_interval_audit_exps_k6.cpp` | Independent interval-rank auditor for `k <= 6` and primes in `{2,3,5,7,11,13}`. |
| `smooth_toolkit.py` | Python CLI/API wrapper choosing kernels automatically. |

Benchmark comparators:

| File | Purpose |
|---|---|
| `benchmarks/dp_pointer_baseline.cpp` | Standard DP/pointer “super ugly number” baseline. |
| `benchmarks/smooth_xplusy_baseline.cpp` | Practical Cartesian-sum / adaptive `X+Y` baseline. |
| `benchmarks/smooth_xplusy_fj_loh_workbench.cpp` | FJ/LOH-style comparator workbench. |

Documentation and paper material:

```text
docs/       claim matrix, certification protocol, hard result table
paper/      LaTeX skeleton and paper-planning notes
results/    selected certified examples and benchmark summaries
```

## Build

Requirements:

```text
C++17 compiler, tested with g++
Python 3.10+
Boost.Multiprecision headers, usually included with a normal g++/libboost-dev installation
```

Build the core kernels:

```bash
./build.sh
```

Run smoke tests:

```bash
python3 tests/test_smoke.py
```

## Basic usage

Unrank a 5-prime fixed smooth number:

```bash
python3 smooth_toolkit.py unrank 2,3,5,7,11 1000000000000
```

Audit a certified case:

```bash
python3 smooth_toolkit.py audit 2,3,5,7,11,13 1000000000000 --timeout 300
```

Compute an exact integer value from an exponent vector:

```bash
python3 smooth_toolkit.py value 2,3,5,7,11,13 55,126,27,54,2,52
```

## Certified examples

The current independent auditor supports prime subsets of:

```text
{2,3,5,7,11,13}
```

with `k <= 6`.

Certified result:

```text
P = (2,3,5,7,11,13)
N = 10^12
exps = [55,126,27,54,2,52]
rank_certified = true
certified_count_le = 1000000000000
```

The auditor proves that exactly `N` `P`-smooth numbers are `<=` the returned candidate.

## Important caveats

1. The fastest kernels use floating-point logarithms for search plus exact integer comparison in the final local band. Treat them as high-performance experimental kernels unless independently audited.
2. The proof-grade audit path currently stops at `k <= 6` for primes `{2,3,5,7,11,13}`.
3. Benchmark results depend on hardware, compiler, and memory pressure. Use the included scripts to regenerate local numbers.
4. This repository targets fixed-prime rank/unrank, not general `Psi(x,y)` smooth-number counting.

## Current best method selection

The wrapper uses this default policy:

```text
k = 3: 3-prime Beatty/floor-sum ranker
k = 5: layer-compressed MITM
other k >= 2: sums-only/adaptive MITM
```

This is an implementation policy, not a theorem of optimality.

## Repo status

Research prototype. The intended next step is to make all headline results independently certifiable and to compare against stronger published Cartesian-sum and fixed-dimensional lattice-counting implementations.

## License

No open-source license has been selected yet. Add a license before publishing if you want others to reuse the code.
