# Claims Matrix

| Claim | Current status | What is needed |
|---|---|---|
| Exponent-stream self-similarity | Proven theorem, easy | Write formal proof |
| Zero positions as ranked faces | Proven theorem, easy | Write formal proof |
| Fast 3-prime random access | Strong computational + theory | Certify arbitrary triples, state complexity |
| 5-prime layer-compressed MITM advantage | Strong empirical, certified benchmark for first five primes | More prime sets, full cert path, memory/time bounds |
| Residual correction inside X+Y/MITM unrankers | Clean certified speed evidence for both full materialized X+Y and sums-only MITM on deterministic k=5/k=6/k=8, N=10^9/10^12 suites | Artifacts: `results/benchmarks/xplusy_corrected_speed_suite_1e9_1e12/` and `results/benchmarks/sums_corrected_speed_suite_1e9_1e12/`; keep the claim implementation-specific |
| Implemented portfolio scoreboard | Clean portfolio artifact: 14/14 solver-agreement rows and 14/14 interval-certified rows across k=3,5,6,8 at N=10^9 and N=10^12 | Artifact: `results/benchmarks/certified_portfolio_suite_1e9_1e12/`; use this as the current fastest-implemented scoreboard, not broad best-known language |
| Current solvers beat canonical heap/frontier generation | Clean certified artifacts: default N=10^5/10^6 suite had 8/8 matching and certified rows, mean heap/current wall ratio 4.5654; N=10^7 stress suite had 4/4 matching and certified rows, current solvers won 4/4, mean heap/current wall ratio 76.5315 | Artifacts: `results/benchmarks/heap_frontier_baseline_suite_1e5_1e6/` and `results/benchmarks/heap_frontier_baseline_suite_1e7/`; narrow sequential-frontier baseline claim only |
| Current solver beats a published sorted-matrix/X+Y selector wrapper on first k=5,6,8 targets | Clean certified first-k MA artifact: 3/3 completed, 3/3 independently rank-certified, 0/3 MA wall-time wins; mean MA/current wall ratio 3.7859 | Artifact: `results/benchmarks/ma_full_unrank_first_k_1e12/`; narrow Mirzaian-Arjomandi sorted-matrix comparison only |
| Serious-comparator gate status | Clean gate artifact: MA full-unrank gate passed; output-sensitive X+Y/LOH probe executed at N=10^6; full FJ and soft-heap gates remain open/not implemented; Barvinok/Normaliz tools installed but smoke commands failed and only rational inputs were exported | Artifact: `results/benchmarks/best_known_comparator_gate/`; status artifact, not a broad best-known claim |
| 6-prime certified benchmark | Achieved for first six primes at N=10^12 | Generalize auditor to arbitrary k/primes |
| Beats standard Super Ugly Number DP baseline | Achieved, large margins | Present as random-access vs sequential-generation distinction |
| Best known for fixed-prime random access | Not yet | Compare against full X+Y/FJ/soft-heap/Barvinok implementations or give theoretical argument |
| Useful in factoring/crypto directly | Not established | Avoid claim |
| General smooth-number counting improvement | Not established | Avoid claim unless compared with Bernstein-style algorithms |
