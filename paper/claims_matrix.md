# Claims Matrix

| Claim | Current status | What is needed |
|---|---|---|
| Exponent-stream self-similarity | Proven theorem, easy | Write formal proof |
| Zero positions as ranked faces | Proven theorem, easy | Write formal proof |
| Fast 3-prime random access | Strong computational + theory | Certify arbitrary triples, state complexity |
| 5-prime layer-compressed MITM advantage | Strong empirical, certified benchmark for first five primes | More prime sets, full cert path, memory/time bounds |
| Layer-compressed vs adaptive X+Y on P=(2,3,5,7,11), N=10^12 | Clean certified benchmark achieved | Artifact: `results/benchmarks/xplusy_vs_layer5_1e12/`; report selection-only X+Y caveat and hardware/compiler conditions |
| Layer-compressed vs adaptive X+Y on all five-prime subsets of {2,3,5,7,11,13}, N=10^12 | Suite protocol added | Run `scripts/run_five_prime_suite.py`; use exact win/loss/certification counts |
| 6-prime certified benchmark | Achieved for first six primes at N=10^12 | Generalize auditor to arbitrary k/primes |
| Beats standard Super Ugly Number DP baseline | Achieved, large margins | Present as random-access vs sequential-generation distinction |
| Best known for fixed-prime random access | Not yet | Compare against full X+Y/FJ/soft-heap/Barvinok implementations or give theoretical argument |
| Useful in factoring/crypto directly | Not established | Avoid claim |
| General smooth-number counting improvement | Not established | Avoid claim unless compared with Bernstein-style algorithms |
