# Hard Result Set for Fixed-Prime Smooth-Number Random Access

## Target problem

Given a fixed ordered prime set P = {p1,...,pk} and a rank N, return the exponent vector e = (e1,...,ek) of the N-th P-smooth number p1^e1 ... pk^ek in increasing numerical order, with 1 counted as rank 1.

## Results that should be theorem/proof level

1. **Log-lattice equivalence**
   Sorting P-smooth numbers is equivalent to sorting lattice points e in N^k by the linear form sum_i e_i log(p_i).

2. **Exponent-stream self-similarity**
   For each prime p_i, if E_i is the exponent stream in increasing P-smooth order, then deleting zeroes from E_i and subtracting 1 gives E_i again.

3. **Zero positions as ranked lower-dimensional faces**
   Zeroes in the p_i-exponent stream are exactly the P\{p_i}-smooth numbers ranked inside the full P-smooth order.

4. **Two-prime counting / floor-sum reduction**
   C_{p,q}(T) = #{(a,b) in N^2 : a log p + b log q <= T} can be written as a sum of floors and evaluated by a Euclidean/continued-fraction style floor-sum routine. This is the rigorous core behind the 3-prime speedup.

5. **Three-prime random-access algorithm**
   C_{p,q,r}(T) = sum_c C_{p,q}(T - c log r). Combined with a certified final selection band, this yields a sublinear-in-N random-access algorithm for fixed triples.

6. **Layer-compressed five-prime decomposition**
   For a split A + (B2 + t*w), a three-prime MITM side can be represented as shifted layers of a two-prime base. The count identity is:
   C_3(X) = sum_{x in B2, x <= X} (floor((X-x)/w)+1).
   Using quotient/residue decomposition avoids materializing the full 3D side.

7. **Certification theorem**
   If the independent interval-rank auditor proves count_le(x)=N and exact comparison verifies the candidate value x, then the returned exponent vector is the N-th P-smooth number.

## Computational benchmark results already obtained

- 3 primes, P=(2,3,5), N=10^12: exps [1126,16930,40], fast ranker about 0.31-0.43s in recent builds.
- 5 primes, P=(2,3,5,7,11), N=10^12: exps [1052,26,33,53,4], layer-compressed solver about 0.37-0.73s depending build; certified by independent auditor.
- 5 primes, P=(2,3,5,7,11), N=10^12: clean `xplusy_vs_layer5` artifact at commit `161445507617a9435f9baadf4e70a3679d9e8d9a` on macOS/x86_64 with Apple clang 21. Layer-compressed full unrank wall time 1.194730s and peak RSS 23,264 KB; adaptive materialized X+Y value selection wall time 2.029418s and peak RSS 153,444 KB. Ratio: X+Y wall / layer wall = 1.6986. Returned exps [1052,26,33,53,4] certified with count_le=10^12.
- 5-prime suite, all six five-prime subsets of {2,3,5,7,11,13}, N=10^12: clean artifact at commit `87da37090b939d217a5ee2a51e2c08101d5d13ac`. All 6 cases completed, all 6 exponent vectors certified, and layer-compressed full unrank beat adaptive materialized X+Y value selection in all 6 wall-time comparisons. X+Y/layer wall-time ratios ranged from 1.1868 to 1.9125, mean 1.4963.
- 6 primes, P=(2,3,5,7,11,13), N=10^12: exps [55,126,27,54,2,52], fast solver about 0.55-0.61s; independent auditor certified count_le=N.
- Analytic lattice-count residual probe, certified k=5/k=6/k=8 target vectors
  at N=10^12: clean artifact at commit
  `f244aeafa0235adbb5f9441d4ff8c01b84e938e4`. All 8 probes completed. The six
  k=5 residuals were within 167 ranks of the certified rank, the k=6 residual
  was about 453 ranks, and the k=8 first-eight-primes residual was about
  1.98e6 ranks, or 2e-6 relative. The layer counter matched the certified rank
  in 6/8 probes and was off by one on two boundary k=5 probes, reinforcing that
  certified claims must use the interval auditor rather than floating-log layer
  counts at exact boundaries. Artifact:
  `results/benchmarks/analytic_count_probe_1e12/`.
- Analytic boundary-correction band probe, same certified k=5/k=6/k=8 target
  ranks at N=10^12: clean artifact at commit
  `4bd49bec50aab80c87e599438f325f55814e6222`. All 8 targets were inside the
  analytic band. The six k=5 cases recovered the expected vector from bands of
  1980-2032 candidates; the k=6 case recovered from a band of 4005 candidates.
  The k=8 first-eight-primes target was inside the band, but the band contained
  5,000,120 candidates and was not enumerated under the 200,000-candidate cap.
  Artifact: `results/benchmarks/analytic_band_probe_1e12/`.
- Residual-corrected analytic boundary band probe, same certified target ranks
  at N=10^12: clean artifact at commit
  `f965d9cb677d26aa83b7b35c8412a90f856c8fc2`. One exact count at the analytic
  center was used to shift the center before enumeration. All 8 targets were
  inside the corrected band, all 8 were enumerated under the 200,000-candidate
  cap, and all 8 recovered the expected vector. The six k=5 bands had 193-204
  candidates, the k=6 band had 997 candidates, and the k=8 first-eight-primes
  band had 19,853 candidates. Artifact:
  `results/benchmarks/analytic_band_corrected_probe_1e12/`.
- Certified residual-corrected analytic boundary band probe, same target suite:
  clean artifact at commit `11687a9ef1b183ff88ce3fea53bca4c281315228`. The
  audit-aware harness recovered all 8 expected vectors and independently
  interval-certified all 8 with `certified_count_le=N`. The k=8 first-eight
  case had a 19,853-candidate corrected band and certified in 18.2 wall seconds
  on the benchmark machine. Artifact:
  `results/benchmarks/analytic_band_corrected_certified_probe_1e12/`.
- Certified corrected-oracle deterministic random suite: clean artifact at
  commit `e53d69b364d568a49e42784892289a15aeb58742`. The suite sampled k=5,
  k=6, and k=8 cases at N=10^9 and N=10^12. All 10 corrected-oracle outputs
  were independently rank-certified and matched both sums-only MITM and full
  materialized X+Y unrank. Corrected oracle used less peak RSS than full X+Y in
  10/10 rows, with mean full-X+Y/corrected RSS ratio about 7.90. Wall-time wins
  were mixed: 4/10 vs full X+Y and 3/10 vs sums-only. Artifact:
  `results/benchmarks/corrected_oracle_random_suite_1e9_1e12/`.
- Residual-corrected full materialized X+Y speed suite: clean artifact at
  commit `4a581fe2fb495997442c4e233fc01ad172a2d4df`. This compares two modes
  of the same full-X+Y implementation: adaptive unrank versus residual-corrected
  analytic band. All 10 corrected outputs were independently rank-certified and
  matched the adaptive full-X+Y vector. Corrected mode won 9/10 wall-time
  comparisons and 10/10 reported in-process timing comparisons, with mean
  adaptive/corrected wall-time ratio about 1.33 and mean reported-time ratio
  about 1.82. Artifact:
  `results/benchmarks/xplusy_corrected_speed_suite_1e9_1e12/`.
- Residual-corrected sums-only MITM speed suite: clean artifact at commit
  `be4e24a8f4adf0947e568dcc57f4797b0b35ca60`. This compares two modes of the
  same sums-only MITM implementation: adaptive unrank versus residual-corrected
  analytic band. All 10 corrected outputs were independently rank-certified and
  matched the adaptive sums-only vector. Corrected mode won 9/10 wall-time
  comparisons, 10/10 reported in-process timing comparisons, and 10/10 peak-RSS
  comparisons, with mean adaptive/corrected ratios about 1.33 for wall time,
  1.61 for reported time, and 1.14 for RSS. This artifact predates exact
  final-band sorting in the sums-only path, so rerun before using it as a
  current-code headline. Artifact:
  `results/benchmarks/sums_corrected_speed_suite_1e9_1e12/`.
- Certified implemented portfolio suite: clean artifact at commit
  `b8357ac6e09c27b781892b9a304d72c184d2abae`. Across k=3,5,6,8 cases at
  N=10^9 and N=10^12, all 14 rows had solver agreement and all 14 rows were
  independently interval-certified. Wall-time winners were xplusy_corrected in
  5 rows, layer_corrected in 4 rows, beatty3 in 3 rows, and xplusy_adaptive in
  2 rows. Artifact:
  `results/benchmarks/certified_portfolio_suite_1e9_1e12/`.
- Sorted-matrix/range-pruning diagnostic, same six five-prime subsets,
  N=10^12: the range-pruned block counter beat the linear saddleback X+Y count
  in only 1/6 cases, with mean block/linear internal time ratio 1.3467. The
  Mirzaian-Arjomandi value-selection probe passed exhaustive small validation
  on 5185/5185 cases, matched the adaptive selected log in 6/6 large cases, but
  won 0/6 timing comparisons, with mean MA/linear internal time ratio 9.4042.
  This is a validated negative/mixed result, not a headline improvement.
  Artifact: `results/benchmarks/sorted_matrix_workbench_1e12/`.
- 8 primes, P=(2,3,5,7,11,13,17,19), N=10^12: exps
  [75,28,9,16,3,22,5,1], high-k sums-only solver output independently
  certified with count_le=N by the k<=8 interval auditor. Artifact:
  `results/benchmarks/k8_certificate_1e12/`.
- Higher-k exploratory: sums-only MITM handled k=10,12 at N=10^12 within
  seconds in local checks, but proof-grade certification beyond k=8 is not yet
  implemented.

## Claims to avoid for now

- Do not claim general best-known algorithm for smooth numbers.
- Do not claim replacement for Bernstein-style general y-smooth counting.
- Do not claim impact on prime generation, factoring, or twin primes.
- Do not claim proof-grade correctness for k>8 high-k variants until
  arbitrary-prime interval auditing is complete.
