#!/usr/bin/env python3
"""Exact certificates for the first memory/access/retention proof iteration.

These check elementary baselines and obstructions, not novelty or ML quality.
All arithmetic relevant to the certificates is rational. No training is run.
"""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import random

ROOT = Path(__file__).resolve().parents[1]


def identity(n):
    return [[Q(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def multiply(a, b):
    return [[sum((x*y for x, y in zip(row, col)), Q(0))
             for col in transpose(b)] for row in a]


def apply(a, x):
    return [sum((u*v for u, v in zip(row, x)), Q(0)) for row in a]


def plus(a, b):
    return [[x+y for x, y in zip(u, v)] for u, v in zip(a, b)]


def scaled(a, factor):
    return [[factor*x for x in row] for row in a]


def inverse(a):
    n = len(a)
    out = [[Q(x) for x in row] + unit for row, unit in zip(a, identity(n))]
    for j in range(n):
        pivot = next(i for i in range(j, n) if out[i][j])
        out[j], out[pivot] = out[pivot], out[j]
        c = out[j][j]
        out[j] = [x/c for x in out[j]]
        for i in range(n):
            if i != j:
                c = out[i][j]
                out[i] = [x-c*y for x, y in zip(out[i], out[j])]
    return [row[n:] for row in out]


def rank(a):
    a = [[Q(x) for x in row] for row in a]
    r = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        c = a[r][j]
        a[r] = [x/c for x in a[r]]
        for i in range(r+1, len(a)):
            c = a[i][j]
            a[i] = [x-c*y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def norm_inf(a):
    return max(sum(abs(x) for x in row) for row in a)


def full_spark_controls():
    fixtures = []
    for n in range(2, 7):
        m = n+2
        f = [[Q((i+1)**j) for j in range(m)] for i in range(n)]
        checked = 0
        for target in range(n):
            for cols in combinations(range(m), n-1):
                deleted = [[f[i][j] for j in cols] for i in range(n) if i != target]
                assert rank(deleted) == n-1
                checked += 1
            square = [row[:n] for row in f]
            unit = [Q(i == target) for i in range(n)]
            update = apply(inverse(square), unit)
            assert all(update) and apply(square, update) == unit
        fixtures.append({"keys": n, "cells": m, "row_deleted_minors_checked": checked,
                         "minimum_exact_overwrite_support": n})
    return fixtures


def leverage_controls():
    fixtures = []
    for n in range(2, 6):
        f = [[Q((i+1)**j) for j in range(n+2)] for i in range(n)]
        for t in range(1, n+1):
            c = [row[:t] for row in f]
            ct = transpose(c)
            normal_inverse = inverse(multiply(ct, c))
            for target in range(n):
                u0 = apply(normal_inverse, c[target])
                ell = sum(x*y for x, y in zip(c[target], u0))
                assert 0 < ell <= 1
                u = [x/ell for x in u0]
                reads = apply(c, u)
                assert reads[target] == 1
                energy = sum(reads[i]**2 for i in range(n) if i != target)
                assert energy == (1-ell)/ell
                # Independent constrained-quadratic stationarity check.
                other = [row for i, row in enumerate(c) if i != target]
                gradient = apply(multiply(transpose(other), other), u)
                assert gradient == [energy*x for x in c[target]]
                fixtures.append({"keys": n, "write_support": t, "target": target,
                                 "leverage": str(ell), "optimal_interference_energy": str(energy)})
    return fixtures


def hadamard_controls():
    fixtures, checked = [], 0
    h = [[1]]
    while len(h) < 8:
        h = [row+row for row in h] + [row+[-x for x in row] for row in h]
        n = len(h)
        assert multiply(transpose(h), h) == scaled(identity(n), n)
        for t in range(1, n+1):
            for cols in combinations(range(n), t):
                for target in range(n):
                    u = [Q(h[target][j], t) for j in cols]
                    reads = apply([[row[j] for j in cols] for row in h], u)
                    assert reads[target] == 1
                    assert sum(reads[i]**2 for i in range(n) if i != target) == Q(n, t)-1
                    checked += 1
            fixtures.append({"keys_and_cells": n, "write_support": t,
                             "sharp_interference_energy": str(Q(n, t)-1)})
    return fixtures, checked


def local_correction_controls():
    rng = random.Random(20261019)
    fixtures = []
    for n in (5, 9, 17):
        # Fixed sparse invertible feature matrix, no private coordinates.
        # F=I+(1/16)*adjacency(cycle); all coordinates occur in three keys.
        f = identity(n)
        for i in range(n):
            f[i][(i-1) % n] += Q(1, 16)
            f[i][(i+1) % n] += Q(1, 16)
        g = multiply(f, transpose(f))
        alpha = 1/max(g[i][i] for i in range(n))
        a = plus(identity(n), scaled(g, -alpha))
        rho = norm_inf(a)
        assert rho < 1
        exact_encoding = multiply(transpose(f), inverse(g))
        assert multiply(f, exact_encoding) == identity(n)
        assert all(all(row[j] for row in exact_encoding) for j in range(n))
        power, polynomial = identity(n), [[Q(0) for _ in range(n)] for _ in range(n)]
        for depth in range(4):
            polynomial = plus(polynomial, power)
            p = scaled(polynomial, alpha)
            b = multiply(transpose(f), p)
            next_power = multiply(power, a)
            assert multiply(f, b) == plus(identity(n), scaled(next_power, -1))
            # All binary assignments for n<=9; fixed-seed samples for n=17.
            corners = (list(product((-1, 1), repeat=n)) if n <= 9 else
                       [tuple(rng.choice((-1, 1)) for _ in range(n)) for _ in range(96)])
            worst = Q(0)
            for values in corners:
                reads = apply(f, apply(b, values))
                error = max(abs(x-y) for x, y in zip(reads, values))
                worst = max(worst, error)
                assert error <= rho**(depth+1)
            values, state = [Q(0)]*n, [Q(0)]*n
            max_stream_error = Q(0)
            for _ in range(160):
                target, new_value = rng.randrange(n), Q(rng.randrange(-8, 9), 8)
                delta = new_value-values[target]
                values[target] = new_value
                state = [x+row[target]*delta for x, row in zip(state, b)]
                assert state == apply(b, values)
                error = max(abs(x-y) for x, y in zip(apply(f, state), values))
                max_stream_error = max(max_stream_error, error)
                assert error <= rho**(depth+1)
            # Fixed-point coefficients and bounded integer accumulators avoid
            # a hidden, horizon-growing floating-point summation error.
            bits = 12
            b_integer = [[(x*(1 << bits)).numerator // (x*(1 << bits)).denominator
                          for x in row] for row in b]
            b_quantized = scaled(b_integer, Q(1, 1 << bits))
            assert all(abs(x-y) < Q(1, 1 << bits)
                       for row, rounded in zip(b, b_quantized) for x, y in zip(row, rounded))
            maximum_coefficient = max(abs(x) for row in b_integer for x in row)
            numerator_bound = n*maximum_coefficient*8
            signed_accumulator_bits = numerator_bound.bit_length()+1
            quantized_bound = rho**(depth+1) + norm_inf(f)*n*Q(1, 1 << bits)
            values_integer, accumulators = [0]*n, [0]*n
            for _ in range(160):
                target, value = rng.randrange(n), rng.randrange(-8, 9)
                delta = value-values_integer[target]
                values_integer[target] = value
                accumulators = [x+row[target]*delta for x, row in zip(accumulators, b_integer)]
                assert accumulators == apply(b_integer, values_integer)
                assert all(abs(x) <= numerator_bound for x in accumulators)
                reads = apply(f, [Q(x, 8*(1 << bits)) for x in accumulators])
                assert max(abs(x-Q(y, 8)) for x, y in zip(reads, values_integer)) <= quantized_bound
            fixtures.append({"keys": n, "feature_width": n, "feature_row_support": 3,
                             "correction_depth": depth, "contraction_rho": str(rho),
                             "uniform_in_horizon_error_bound": str(rho**(depth+1)),
                             "binary_value_assignments_checked": len(corners),
                             "binary_assignments_exhaustive": n <= 9,
                             "maximum_checked_binary_error": str(worst),
                             "maximum_checked_stream_error": str(max_stream_error),
                             "stream_overwrites": 160,
                             "fixed_point_stream_overwrites": 160,
                             "fixed_point_coefficient_fractional_bits": bits,
                             "fixed_point_accumulator_signed_bits": signed_accumulator_bits,
                             "fixed_point_uniform_error_bound": str(quantized_bound),
                             "stored_encoding_nonzeros": sum(bool(x) for row in b for x in row),
                             "maximum_updated_state_cells": max(sum(bool(row[j]) for row in b) for j in range(n)),
                             "exact_interpolating_encoding_nonzeros": n*n})
            power = next_power
    return fixtures


def dynamic_key_controls():
    """Local coefficient updates versus an independent full recomputation.

    Dense Gram formation is an oracle here; this is not a runtime benchmark.
    The selected-column construction itself traverses sparse supports.
    """
    rng = random.Random(20261020)
    fixtures = []
    n, alpha, rho, bits = 11, Q(4, 5), Q(1, 2), 12

    def gram(f, active):
        g = multiply(f, transpose(f))
        for i in range(n):
            if not active[i]:
                g[i][i] += 1/alpha
        return g

    def encoding(f, g, depth):
        a = plus(identity(n), scaled(g, -alpha))
        assert norm_inf(a) <= rho
        power, total = identity(n), identity(n)
        for _ in range(depth):
            power = multiply(power, a)
            total = plus(total, power)
        return multiply(transpose(f), scaled(total, alpha))

    def local_column(f, g, depth, target):
        a = plus(identity(n), scaled(g, -alpha))
        columns = [{i: a[i][j] for i in range(n) if a[i][j]} for j in range(n)]
        power, total = {target: Q(1)}, {target: alpha}
        for _ in range(depth):
            following = {}
            for j, x in power.items():
                for i, coefficient in columns[j].items():
                    following[i] = following.get(i, Q(0)) + coefficient*x
            power = {i: x for i, x in following.items() if x}
            for i, x in power.items():
                total[i] = total.get(i, Q(0)) + alpha*x
        out = [Q(0)]*n
        for i, x in total.items():
            for j, coefficient in enumerate(f[i]):
                if coefficient:
                    out[j] += coefficient*x
        return out

    def rounded(b):
        return [[(x*(1 << bits)).numerator // (x*(1 << bits)).denominator
                 for x in row] for row in b]

    for depth in range(3):
        active = [True]*n
        f = identity(n)
        for i in range(n):
            f[i][(i-1) % n] = f[i][(i+1) % n] = Q(1, 16)
        values = [rng.randrange(-4, 5) for _ in range(n)]
        g = gram(f, active)
        b = encoding(f, g, depth)
        bi = rounded(b)
        state, integer_state = apply(b, [Q(x, 4) for x in values]), apply(bi, values)
        changed_max, visited_max = 0, 0
        counts = {"insertions": 0, "deletions": 0, "replacements": 0}
        for step in range(36):
            target = rng.randrange(n)
            new_active = rng.randrange(4) != 0
            # Force both occupancy transitions into every certificate.
            if step in (0, 1):
                target, new_active = 0, step == 1
            updated_f = [row[:] for row in f]
            updated_f[target] = [Q(0)]*n
            if new_active:
                updated_f[target][target] = Q(1)
                for j in ((target-1) % n, (target+1) % n):
                    updated_f[target][j] = Q(rng.randrange(3), 32)
            updated_active = active[:]
            updated_active[target] = new_active
            updated_g = gram(updated_f, updated_active)
            expected_b = encoding(updated_f, updated_g, depth)
            # Union of old and new overlap graphs, excluding self loops.
            neighborhood = {target}
            for _ in range(depth):
                neighborhood |= {j for i in tuple(neighborhood) for j in range(n)
                                 if i != j and (g[i][j] or updated_g[i][j])}
            assert all(expected_b[i][j] == b[i][j]
                       for j in range(n) if j not in neighborhood for i in range(n))
            updated_b = [row[:] for row in b]
            for j in neighborhood:
                column = local_column(updated_f, updated_g, depth, j)
                assert column == [row[j] for row in expected_b]
                for i, x in enumerate(column):
                    updated_b[i][j] = x
            assert updated_b == expected_b
            updated_bi = rounded(updated_b)
            for j in neighborhood:
                for i in range(n):
                    state[i] += (updated_b[i][j]-b[i][j])*Q(values[j], 4)
                    integer_state[i] += (updated_bi[i][j]-bi[i][j])*values[j]
            new_value = rng.randrange(-4, 5) if new_active else 0
            delta = new_value-values[target]
            for i in range(n):
                state[i] += updated_b[i][target]*Q(delta, 4)
                integer_state[i] += updated_bi[i][target]*delta
            values[target] = new_value
            assert state == apply(updated_b, [Q(x, 4) for x in values])
            assert integer_state == apply(updated_bi, values)
            assert all(abs(x) <= n*(1 << (bits+2))*4 for x in integer_state)
            error = max(abs(x-Q(y, 4)) for x, y in zip(apply(updated_f, state), values))
            assert error <= rho**(depth+1)
            quantized_state = [Q(x, 4*(1 << bits)) for x in integer_state]
            quantized_error = max(abs(x-Q(y, 4)) for x, y in zip(apply(updated_f, quantized_state), values))
            assert quantized_error <= rho**(depth+1)+Q(9, 8)*n*Q(1, 1 << bits)
            counts['replacements' if active[target] == new_active else
                   'insertions' if new_active else 'deletions'] += 1
            changed_max = max(changed_max, sum(any(updated_b[i][j] != b[i][j]
                                                 for i in range(n)) for j in range(n)))
            visited_max = max(visited_max, len(neighborhood))
            f, g, b, bi, active = updated_f, updated_g, updated_b, updated_bi, updated_active
        fixtures.append({"slot_capacity": n, "correction_depth": depth,
                         "fixed_alpha": str(alpha), "uniform_rho_bound": str(rho),
                         "feature_row_support_cap": 3, "gram_off_diagonal_degree_cap": 4,
                         "operations": 36, "operation_counts": counts,
                         "maximum_changed_encoding_columns": changed_max,
                         "maximum_visited_union_graph_ball": visited_max,
                         "quantized_coefficients_fractional_bits": bits,
                         "certificate_only_not_runtime_benchmark": True})
    return fixtures


def main():
    spark = full_spark_controls()
    leverage = leverage_controls()
    hadamard, hadamard_count = hadamard_controls()
    local = local_correction_controls()
    dynamic = dynamic_key_controls()
    sources = ["scripts/run_retention_tradeoff_checks.py", "docs/interleaving_research.md"]
    report = {"schema_version": 1,
              "command": "python3 scripts/run_retention_tradeoff_checks.py",
              "scope": "Exact algebraic certificates for elementary access/retention bounds and a classical fixed-feature polynomial correction baseline. No new neural training, publication claim or superiority over the strongest prior art.",
              "arithmetic": "fractions.Fraction",
              "full_spark_overwrite_fixtures": spark,
              "leverage_interference_fixtures": leverage,
              "hadamard_sharp_frontier": hadamard,
              "hadamard_target_support_cases": hadamard_count,
              "local_polynomial_correction_fixtures": local,
              "dynamic_key_local_correction_fixtures": dynamic,
              "uniform_proof_location": "docs/interleaving_research.md, Associative-memory goal: first proof iteration",
              "new_publishable_advantage_established": False,
              "summary": {"row_deleted_minors": sum(x["row_deleted_minors_checked"] for x in spark),
                          "leverage_optimizer_cases": len(leverage),
                          "hadamard_target_support_cases": hadamard_count,
                          "fixed_feature_correction_cases": len(local),
                          "total_stream_overwrites": sum(x["stream_overwrites"] + x["fixed_point_stream_overwrites"] for x in local),
                          "dynamic_key_edit_cases": len(dynamic),
                          "dynamic_key_edits": sum(x["operations"] for x in dynamic)},
              "source_sha256": {f: sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    target = ROOT/"results/research/retention_tradeoff_checks.json"
    target.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
