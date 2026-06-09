# Claims Matrix

| Claim | Current status | What is needed |
|---|---|---|
| Exponent-stream self-similarity | Proven theorem, easy | Write formal proof |
| Zero positions as ranked faces | Proven theorem, easy | Write formal proof |
| Fast 3-prime random access | Strong computational + theory | Certify arbitrary triples, state complexity |
| 5-prime layer-compressed MITM advantage | Strong empirical, certified benchmark for first five primes | More prime sets, full cert path, memory/time bounds |
| Layer-compressed vs adaptive X+Y on P=(2,3,5,7,11), N=10^12 | Clean certified benchmark achieved | Artifact: `results/benchmarks/xplusy_vs_layer5_1e12/`; report selection-only X+Y caveat and hardware/compiler conditions |
| Layer-compressed vs adaptive X+Y on all five-prime subsets of {2,3,5,7,11,13}, N=10^12 | Clean certified suite win: 6/6 certified, 6/6 layer wall-time wins | Artifact: `results/benchmarks/five_prime_suite_1e12/`; report selection-only X+Y caveat and hardware/compiler conditions |
| Layer-compressed vs full materialized X+Y unrank on all five-prime subsets of {2,3,5,7,11,13}, N=10^12 | Clean certified suite win: 6/6 matching vectors, 6/6 both independently certified, 6/6 layer wall-time wins | Artifact: `results/benchmarks/full_xplusy_suite_1e12/`; report hardware/compiler conditions and avoid broader "best known" wording |
| Analytic asymptotic bracketing improves layer-compressed unranking | Clean optimization evidence: 6/6 analytic brackets, mean final band 11286.8 -> 345.8, reported layer time improved 5/6 cases | Artifact: `results/benchmarks/full_xplusy_suite_1e12/`; non-proof seed only, correctness still comes from audits |
| Analytic lattice count residuals near certified targets | Initial diagnostic artifact: 8/8 probes completed; k=5 residuals within 167 ranks at N=10^12, k=8 first-eight residual about 1.98M ranks / 2e-6 relative | Artifact: `results/benchmarks/analytic_count_probe_1e12/`; diagnostic only, no rigorous error bound yet |
| Analytic boundary-correction bands can recover certified targets | Initial prototype: 8/8 targets inside analytic band; 7/8 enumerated and recovered expected vectors; k=8 inside but too wide under 200k candidate cap | Artifact: `results/benchmarks/analytic_band_probe_1e12/`; experimental, endpoint counts are not certificates |
| One exact residual correction makes analytic bands enumerably small on tested k=5/k=6/k=8 targets | Clean experimental artifact: 8/8 targets inside, enumerated, and recovered expected vectors; k=8 first-eight band shrank to 19,853 candidates under the 200k cap | Artifact: `results/benchmarks/analytic_band_corrected_probe_1e12/`; still uses floating-log endpoint counts and needs interval-audited certification |
| Practical sorted-matrix/FJ-style range pruning improves X+Y selection at N=10^12 | Negative/mixed: 2/6 narrow wins, 4/6 losses, mean block/linear ratio 1.2765 | Artifact: `results/benchmarks/sorted_matrix_workbench_1e12/`; diagnostic only, not a faithful Frederickson-Johnson implementation |
| Mirzaian-Arjomandi sorted-matrix value selection improves practical X+Y selection at N=10^12 | Negative: 6/6 exact log matches, 0/6 timing wins, mean MA/linear ratio 9.0509 | Artifact: `results/benchmarks/sorted_matrix_workbench_1e12/`; value-selection only, not full unrank |
| 6-prime certified benchmark | Achieved for first six primes at N=10^12 | Generalize auditor to arbitrary prime sets |
| 8-prime certified sums-only random access | Achieved for first eight primes at N=10^12: high-k sums-only solver output independently rank-certified | Artifact: `results/benchmarks/k8_certificate_1e12/`; still a fixed-instance certificate, not a broad best-known claim |
| Beats standard Super Ugly Number DP baseline | Achieved, large margins | Present as random-access vs sequential-generation distinction |
| Best known for fixed-prime random access | Not yet | Compare against full X+Y/FJ/soft-heap/Barvinok implementations or give theoretical argument |
| Useful in factoring/crypto directly | Not established | Avoid claim |
| General smooth-number counting improvement | Not established | Avoid claim unless compared with Bernstein-style algorithms |
