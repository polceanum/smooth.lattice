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
paper/      quantitative manuscript and separate research notes
results/    selected certified examples and benchmark summaries
```

Current research returns to the original occurrence-count and prime-exponent
interleaving observations. The [face reconstruction study](docs/occurrence_reconstruction.md)
proves that face prefixes involving O(log n) coordinates determine an n-term
full prefix, and that logarithmic dimension is necessary in the worst case,
including for prime bases. Exact checks include two four-prime systems whose
proper faces agree through 1,000 terms while their full words differ at rank 10.
Novelty and publication significance of this new study remain unresolved.
The [goal record](results/research/occurrence_reconstruction_goal.json) states
the next proof and primary-literature gates. The separate ML track is paused;
its results are preserved without a claimed ML advantage.
The [single-stream continuation](docs/occurrence_reconstruction.md#one-exponent-stream-determines-the-normalized-ordering)
also proves normalized-weight recovery from one infinite exponent stream or
its binary zero mask, an explicit finite height-error bound with calibration,
and absolute prime-basis uniqueness for three or more primes using the
classical six exponentials theorem. Arbitrarily long finite stream prefixes
can nevertheless hide a rank-6 ordering change. Priority and an efficient
recovery algorithm remain unestablished.
The [finite-recovery continuation](docs/occurrence_reconstruction.md#finite-recovery-under-a-quantitative-separation-promise)
supplies an explicit sufficient mask length under stated weight bounds and
a Diophantine separation promise, including generator identification and
all observation costs. Its bound is polynomial in numerical rank for fixed
class parameters but can be enormous. This is a finite-information result;
priority, practical recovery and an unranking advantage remain unresolved.
The [sharp dimension continuation](docs/occurrence_reconstruction.md#a-sharp-dimension-threshold-for-finite-interleaving)
determines the worst-case required face dimension up to a bounded additive
error: `(log n + 0.5 log log n) / log(3 sqrt(3)/2) + O(1)`, including for
prime bases. A uniform composition-layer moment bound explains the sharp
interaction-rank scale. Priority and publication significance remain open.

The [mathematical research note](docs/interleaving_research.md) targets general
structure and certification for fixed-prime rank/unrank. It separates known
signature-sequence theory, general occurrence-transform laws, a restricted
cycle classification, and an obstruction to extending finite pairwise patterns.
It also derives a theoretical generic unranking baseline from published log
separation and lattice-counting results. New work must add something beyond
that baseline. Exact research checks and witnesses are saved; novelty and prime
distribution applications remain unestablished.

The [completed quantitative goal](results/research/quantitative_goal.json)
required a new quantitative theorem, a verified literature gap and a useful
mathematical consequence. Its original target was Wagner's published
logarithmic gap in the number of three-dimensional corner cuts, the shapes
underlying smooth exponent prefixes.
The existing proofs and closest primary results have been compared at the
theorem level, including the full Vasiliev--Pavlov finite-ordering paper and
Wagner's stronger face-flag bound. The final internal assessment supports a
focused submission under the documented literature comparison and the user's
explicit final-version assumption.

The quantitative manuscript proves the matching pointwise count
C_d(n)=Theta_d(n^(d-1)(log n)^(d-1)) in every fixed dimension d>=2,
including the original three-dimensional target. It counts primitive normal
directions and uses two integer lattice patches to prevent collisions at a
prescribed size. Exact witnesses test the delicate step for 2,280 ordered pairs of
normals at five boundary heights each. The fresh internal proof audit found
no unresolved logical step. The earlier cumulative construction is preserved
in the note. The count also gives order log n for the mean number of
completions of a coordinate slice, addressing Wagner's suggested route. The
[anonymous mathematical manuscript](paper/corner_cut_growth/main.tex) gives
the complete proof, explicit constants, source comparisons and an
AI-assistance disclosure separately from the earlier structural note. It also
proves the same growth order for complete ordered prefixes and a sharp
worst-case prime-size bound: for fixed dimension and lower prime bound,
O(n)-bit primes realize every admissible n-term pattern, and some patterns
require Omega(n) bits. Known finite-ordering, representation and qualitative
prime-realization results are credited explicitly. This supplies no fixed-prime
unranking speedup or new prime-gap estimate.
The [PDF](output/pdf/corner_cut_growth.pdf) is compiled and visually checked.
The [acceptance audit](results/research/corner_cut_growth_audit.json) records
the mathematical contribution and the reading scope of each comparison.
The novelty position is a documented research assessment; no independent
referee assessment or journal acceptance is claimed.

The [next research plan](docs/interleaving_research.md#higher-impact-research-plan-universal-algebraic-complexity)
aims to explain how algebraic relations change the number of monomial initial
ideals of arbitrary finite polynomial algebras. The
[new goal record](results/research/structural_complexity_goal.json) separates
the full extremal ambition from the first literature and proof-feasibility
gates. Success requires a uniform asymptotic advance and a mathematical
consequence; examples or a roadmap alone do not complete it. This is a new
research objective, not an additional theorem of the existing manuscript.

The existing structural results are preserved as a separate note about the
geometry and containment intervals of pyramidal partitions. Applying Edelman's
region-poset theorem gives a cutting-layer formula for their Mobius function
in every dimension, addressing an observation in an Elizalde--Galvan preprint.
It also gives a subpartition-counting recurrence over removable exposed faces,
with at most linearly many contributing faces at each three-dimensional top.
This is a local structural bound, not a claimed faster overall enumeration method.
An infinite family also refutes the coordinate-slice recognition conjecture
in Galvan's 2023 thesis: triangular sections need not fit one cutting plane.
The [paper draft](paper/pyramidal_mobius.tex), exact checks of 27,100 intervals,
and 1,280 slice certificates are saved. These do not meet the new quantitative
goal. The draft is not submission-ready:
the final published version and
later literature still need comparison. The preceding integer prefix-compression
bound is retained as background; no faster rank/unrank method is claimed.

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
