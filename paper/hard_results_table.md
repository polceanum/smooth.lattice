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
- Five-prime suite, all six subsets of {2,3,5,7,11,13}, N=10^12:
  layer-compressed full unranking beat full materialized X+Y unranking in 6/6
  wall-time comparisons on the recorded macOS/x86_64 Apple-clang machine. Both
  methods returned the same exponent vector in every case, and both outputs were
  independently certified. Artifact: `results/benchmarks/full_xplusy_suite_1e12/`.
- Sorted-matrix/range-pruning diagnostic, same six five-prime subsets,
  N=10^12: the range-pruned block counter beat the linear saddleback X+Y count
  in only 2/6 cases, with mean block/linear internal time ratio 1.2765. The
  Mirzaian-Arjomandi value-selection probe matched the adaptive selected log in
  6/6 cases but won 0/6 timing comparisons, with mean MA/linear internal time
  ratio 9.0509. This is a negative/mixed result, not a headline improvement.
  Artifact: `results/benchmarks/sorted_matrix_workbench_1e12/`.
- 6 primes, P=(2,3,5,7,11,13), N=10^12: exps [55,126,27,54,2,52], fast solver about 0.55-0.61s; independent auditor certified count_le=N.
- Higher-k exploratory: sums-only MITM handled k=8,10,12 at N=10^12 within seconds in the container, but proof-grade certification beyond k=6 is not yet implemented.

## Claims to avoid for now

- Do not claim general best-known algorithm for smooth numbers.
- Do not claim replacement for Bernstein-style general y-smooth counting.
- Do not claim impact on prime generation, factoring, or twin primes.
- Do not claim proof-grade correctness for fastest high-k variants until arbitrary-prime interval auditing is complete.
