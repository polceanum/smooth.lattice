# Claim Language For The Paper

The paper should make claims about mathematical objects, algorithms,
certificates, and measured comparisons. The repository is reproducibility
evidence, not the subject of the scientific claim.

## Claim Types

### Theorem-level claims

Use for statements intended to be proven in the paper.

- Fixed-prime smooth-number ordering is equivalent to ordering lattice points
  in `N^k` by the linear form `sum_i e_i log p_i`.
- For any coordinate exponent stream, deleting zeroes and subtracting one from
  the remaining entries recovers the same stream.
- Zero positions in one coordinate correspond to lower-dimensional faces ranked
  inside the full lattice order.
- If an independent rank-count certificate proves `count_le(P, x) = N`, then
  `x` is the `N`-th `P`-smooth number.

### Algorithmic claims

Use for method descriptions and complexity or structure arguments.

- The three-prime algorithm reduces the count subproblem to two-prime
  floor-sum counts.
- The five-prime layer-compressed algorithm represents one three-dimensional
  meet-in-the-middle side as shifted two-dimensional layers.
- The sums-only high-`k` variant delays exponent reconstruction and reduces
  stored state, but remains experimental unless independently audited.

### Certified computational claims

Use when a returned exponent vector is independently certified.

```text
For P=(2,3,5,7,11) and N=10^12, the layer-compressed algorithm returned
e=(1052,26,33,53,4). The independent interval-rank auditor verified
count_le(P, prod_i p_i^e_i) = 10^12.
```

### Experimental comparison claims

Use when timings and memory are measured under recorded conditions.

```text
On the recorded macOS/x86_64 machine with Apple clang 21, for
P=(2,3,5,7,11) and N=10^12, layer-compressed full unranking took 1.194730s
wall time, while the adaptive materialized X+Y value selector took 2.029418s.
```

Always report:

- hardware and compiler conditions;
- exact prime set and rank;
- whether the baseline returns a value, a rank count, or a full exponent vector;
- whether the algorithm output was independently certified;
- negative results and failures.

## Avoid

Do not write:

- "the repository proves";
- "best known" without a specific benchmark class and comparator set;
- "state of the art" without a literature-quality comparison;
- "smooth-number algorithms" when the claim is only about fixed-prime
  random-access unranking;
- "certified" unless `rank_certified=true` or an equivalent proof is recorded.

## Safe Current Wording

The current one-case benchmark supports:

```text
For P=(2,3,5,7,11) and N=10^12, layer-compressed full unranking
outperformed our practical adaptive materialized X+Y value-selection baseline
on the recorded machine/compiler, and the returned exponent vector was
independently certified.
```

The certified five-prime suite supports:

```text
Across all six five-prime subsets of {2,3,5,7,11,13} at N=10^12,
layer-compressed full unranking outperformed the practical adaptive
materialized X+Y value-selection baseline on the recorded machine/compiler.
All six returned exponent vectors were independently certified.
```

Attach the caveat that the `X+Y` baseline is value selection only, not full
exponent reconstruction, and report the exact hardware/compiler conditions.
