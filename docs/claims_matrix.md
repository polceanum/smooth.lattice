# Claims Matrix

| Claim | Current status | What is needed |
|---|---|---|
| Exponent-stream self-similarity | Proven theorem, easy | Write formal proof |
| Zero positions as ranked faces | Proven theorem, easy | Write formal proof |
| Fast 3-prime random access | Strong computational + theory | Certify arbitrary triples, state complexity |
| 5-prime layer-compressed MITM advantage | Strong empirical, certified benchmark for first five primes | More prime sets, full cert path, memory/time bounds |
| Residual correction inside X+Y/MITM unrankers | Clean certified speed evidence for both full materialized X+Y and sums-only MITM on deterministic k=5/k=6/k=8, N=10^9/10^12 suites | Artifacts: `results/benchmarks/xplusy_corrected_speed_suite_1e9_1e12/` and `results/benchmarks/sums_corrected_speed_suite_1e9_1e12/`; keep the claim implementation-specific |
| Implemented portfolio scoreboard | Clean portfolio artifact: 14/14 solver-agreement rows and 14/14 interval-certified rows across k=3,5,6,8 at N=10^9 and N=10^12 | Artifact: `results/benchmarks/certified_portfolio_suite_1e9_1e12/`; use this as the current fastest-implemented scoreboard, not broad best-known language |
| 6-prime certified benchmark | Achieved for first six primes at N=10^12 | Generalize auditor to arbitrary k/primes |
| Beats standard Super Ugly Number DP baseline | Achieved, large margins | Present as random-access vs sequential-generation distinction |
| Best known for fixed-prime random access | Not yet | Compare against full X+Y/FJ/soft-heap/Barvinok implementations or give theoretical argument |
| Useful in factoring/crypto directly | Not established | Avoid claim |
| General smooth-number counting improvement | Not established | Avoid claim unless compared with Bernstein-style algorithms |
