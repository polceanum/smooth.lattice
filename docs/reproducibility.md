# Reproducibility checklist

Before using a benchmark in a paper or public claim, record:

- git commit hash;
- compiler and flags;
- CPU model and RAM;
- exact command line;
- raw stdout/stderr;
- whether the result was independently audited;
- whether the method used floating log search.

Recommended commands:

```bash
./build.sh
python3 tests/test_smoke.py
bash scripts/run_certified_examples.sh
bash scripts/run_dp_comparison.sh
python3 scripts/run_xplusy_vs_layer5.py
python3 scripts/run_full_xplusy_suite.py
python3 scripts/run_sorted_matrix_workbench.py
python3 scripts/run_heap_frontier_baseline_suite.py
python3 scripts/run_ma_full_unrank_suite.py
python3 scripts/run_best_known_comparator_gate.py
```

Dependency notes:

- The C++ kernels require a C++17 compiler and Boost.Multiprecision headers.
- On Debian/Ubuntu, install `g++` and `libboost-dev`.
- On macOS with Homebrew, install `boost`; if the compiler cannot find it automatically, pass the include path through `CXXFLAGS`.
- `build.sh` runs a small Boost header preflight before compiling the kernels.

For headline claims, keep JSON or text output under `results/local/` and copy final curated results to `results/certified/` or `results/benchmarks/`.

The `run_xplusy_vs_layer5.py` harness writes `report.json`, `summary.csv`, and
`report.md` under `results/local/xplusy_vs_layer5_<timestamp>/`. Its default
case is `P=(2,3,5,7,11), N=10^12`; it compares adaptive Cartesian-sum value
selection against the layer-compressed full unrank solver and, unless disabled,
audits the returned exponent vector with the independent interval auditor.

The full-X+Y and sorted-matrix workbench harnesses follow the same artifact
format. The sorted-matrix workbench is diagnostic: its Mirzaian-Arjomandi row is
value selection only, its LOH row is capped, and neither is a full exponent-vector
unrank certificate.

The heap-frontier harness also writes `report.json`, `summary.csv`, and
`report.md`. It compares a canonical full-vector frontier generator against the
DP pointer baseline and current solver, then interval-audits matching rows.

For Barvinok-style external comparisons, the optional local setup used for the
comparator-gate artifact was:

```bash
conda create -n smooth-lattice-count -c conda-forge barvinok normaliz pynormaliz -y
```

The gate records tool versions, smoke-command outputs, rationalized `.ine`
inputs, and a bounded PyNormaliz count probe. The PyNormaliz probe must pass
known toy simplexes before any target result is reported. Certified-target
Normaliz counts use rationalized log simplexes and are comparator measurements,
not correctness certificates for the original irrational-log rank problem.
