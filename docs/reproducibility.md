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
```

For headline claims, keep JSON or text output under `results/local/` and copy final curated results to `results/certified/` or `results/benchmarks/`.
