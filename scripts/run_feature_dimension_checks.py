#!/usr/bin/env python3
"""Exact algebraic certificates, not additional model experiments.

Uniform statements are proved/cited in docs/interleaving_research.md.
Finite controls do not establish novelty or a trained-model advantage.
"""
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from math import comb
from pathlib import Path
import random

ROOT = Path(__file__).resolve().parents[1]


def lower_bound(size, dimension, power):
    return comb(power + dimension - 1, dimension) * size - (power - 1) * comb(power + dimension - 1, dimension - 1)


def long_simplex(size, dimension):
    zero = (0,) * dimension
    return {zero, *(tuple(int(i == j) for i in range(dimension)) for j in range(1, dimension)),
            *((i,) + (0,) * (dimension - 1) for i in range(1, size - dimension + 1))}


def sumset(a, b):
    return {tuple(x + y for x, y in zip(u, v)) for u in a for v in b}


def rank(matrix):
    a = [[F(x) for x in row] for row in matrix]
    r = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = a[r][j]
        a[r] = [x / scale for x in a[r]]
        for i in range(r + 1, len(a)):
            scale = a[i][j]
            a[i] = [x - scale * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def main():
    fixtures = []
    for d in range(1, 6):
        for size in (d + 1, d + 3, d + 5):
            a = long_simplex(size, d)
            assert len(a) == size and rank([list(x) for x in a]) == d
            current = {(0,) * d}
            for power in range(1, 7):
                current = sumset(current, a)
                assert len(current) == lower_bound(size, d, power)
                if power in (2, 4, 6):
                    fixtures.append({"dimension": d, "base_width": size, "power": power,
                                     "sumset_cardinality": len(current), "published_sharp_bound": lower_bound(size, d, power)})
    widths = []
    for d in (2, 3, 8, 32, 63):
        a = long_simplex(64, d)
        width = len(sumset(a, a))
        assert width == lower_bound(64, d, 2)
        widths.append({"base_width": 64, "independent_directions": d, "minimum_quadratic_width": width})
    box = set(product(range(8), repeat=2))
    assert len(box) == 64 and len(sumset(box, box)) == 225
    assert lower_bound(64, 3, 2) == 250 > 225
    full_eight = long_simplex(29, 8)
    assert len(sumset(full_eight, full_eight)) == 225

    # A=(0,1), positive unit coefficients, normalized degree-two kernel.
    epsilon = F(1, 10)
    points = [-epsilon, epsilon]
    features = [[1 / (1 + x*x), x / (1 + x*x), x*x / (1 + x*x)] for x in points]
    weights = [1, 2, 1]
    kernel = [[sum(c * x * y for c, x, y in zip(weights, u, v)) for v in features] for u in features]
    mu = F(9801, 10201)
    assert kernel == [[1, mu], [mu, 1]] and rank(kernel) == 2
    state = [F(0)] * 3
    for key, value in zip(features, (F(1), F(-1))):
        prior = sum(c * x * h for c, x, h in zip(weights, key, state))
        state = [h + x * (value - prior) for h, x in zip(state, key)]
    reads = [sum(c * x * h for c, x, h in zip(weights, key, state)) for key in features]
    assert reads == [1 - mu - mu*mu, F(-1)] and reads[0] < 0

    corners = list(product((-1, 1), repeat=2))
    additive = [[1, x, y] for x, y in corners]
    with_binding = [[1, x, y, x*y] for x, y in corners]
    assert rank(additive) == 3 and rank(with_binding) == 4

    # Independent projectors cannot be pushed through a scalar sum sketch.
    before = [(1, 0), (0, 1)]
    assert [sum(x) for x in before] == [1, 1]
    after = [x[0] for x in before]
    assert after == [1, 0]

    rng = random.Random(20261018)
    retrieval_checks = 0
    for n in range(2, 17):
        for k in range(1, min(n, 4) + 1):
            left = [rng.randrange(-5, 6) for _ in range(n)]
            right = [rng.randrange(-5, 6) for _ in range(n)]
            top_left = sorted(range(n), key=lambda i: (-left[i], i))[:k]
            top_right = sorted(range(n), key=lambda j: (-right[j], j))[:k]
            ordering = lambda pair: (-left[pair[0]] - right[pair[1]], pair[0], pair[1])
            exact = sorted(product(range(n), repeat=2), key=ordering)[:k]
            factored = sorted(product(top_left, top_right), key=ordering)[:k]
            assert exact == factored
            retrieval_checks += 1

    sources = ["scripts/run_feature_dimension_checks.py", "docs/interleaving_research.md"]
    report = {"schema_version": 1, "command": "python3 scripts/run_feature_dimension_checks.py",
              "scope": "Exact finite certificates for established rank/sumset principles and elementary failure controls; no new model training or novelty claim.",
              "published_higher_sumset_source": "https://arxiv.org/html/2409.05638v1#S7",
              "kernel_rank_source": "https://arxiv.org/html/2106.02755v1#S3.SS1",
              "product_key_source": "https://arxiv.org/html/1907.05242v2#S3.SS1",
              "sharp_sumset_fixtures": fixtures, "quadratic_dimension_frontier": widths,
              "original_box": {"base_width": 64, "directions": 2, "collected_width": 225},
              "fixed_64_to_225_maximum_affine_dimension": 2,
              "full_eight_direction_example": {"base_width": 29, "collected_width": 225,
                                               "condition": "Long simplex; all eight coordinate monomials occur, but higher-degree features are strongly anisotropic."},
              "exact_delta_interference": {"keys": [str(x) for x in points], "desired_values": [1, -1],
                                           "kernel_correlation": str(mu), "kernel_rank": 2,
                                           "gram_determinant": str(1 - mu*mu), "reads_after_two_writes": [str(x) for x in reads],
                                           "first_read_decimal": float(reads[0]), "arithmetic": "fractions.Fraction"},
              "block_additive_binding_counter": {"additive_feature_rank": 3, "parity_augmented_rank": 4,
                                                  "points": corners, "desired_parity": [x*y for x, y in corners]},
              "arbitrary_channel_gate_counter": {"compressed_before": [1, 1], "compressed_after_same_gate": after},
              "exact_product_key_topk_checks_including_ties": retrieval_checks,
              "new_architecture_advantage_established": False,
              "source_sha256": {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sources}}
    target = ROOT / "results/research/full_dimensional_factorization_checks.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"sharp_sumset_checks": len(fixtures), "dimension_frontier": widths,
                      "first_read_despite_full_rank": str(reads[0]), "product_key_topk_checks": retrieval_checks}, indent=2))


if __name__ == "__main__":
    main()
