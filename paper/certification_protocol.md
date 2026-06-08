# Certification Protocol

## Required for paper-quality results

1. **Exact problem statement**
   Fix P and N. Rank starts at 1, with 1 as the first P-smooth number.

2. **Independent solver/auditor separation**
   - Fast solver returns exponent vector e.
   - Auditor independently verifies rank of x = prod p_i^e_i.
   - The auditor should not reuse the solver's candidate-generation logic except for reading the candidate vector.

3. **Rank certificate**
   For candidate x, prove:
   count_le(P,x) = N.
   Since P-smooth numbers are unique by exponent vector, this proves x is the N-th element.

4. **Interval-log counting**
   Use rational intervals for log(p_i). For each split partial sum, compute lower and upper rational bounds. Count pairs that are certainly <= log x, reject pairs certainly > log x, and resolve only ambiguous boundary cases using exact big-integer comparison.

5. **Ambiguity log**
   Record:
   - number of partial sums on each side,
   - number of certainly counted pairs,
   - number of ambiguous pairs,
   - number of exact comparisons used,
   - final certified count.

6. **Fuzz tests**
   For small N and prime sets, compare against exact DP/pointer generation:
   - all N <= 10,000 for selected prime sets,
   - random N up to 1,000,000,
   - random prime subsets up to k=6.

7. **No floating tolerance in certified path**
   Any floating tolerance may remain in experimental high-performance solvers, but not in the certified auditor.

8. **Reproducibility**
   Provide:
   - compiler version,
   - command lines,
   - CPU/RAM info,
   - raw stdout logs,
   - generated JSON certificates.
