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
