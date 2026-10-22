#!/usr/bin/env python3
"""Actual prefix gaps, all-vector rank bounds and full longest-gap cone controls."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_generic_hook_checks import bareiss
from run_hook_bridge_checks import COUNTS as REPLAY_COUNTS, actual, model, normals, polygon
from run_hook_central_prefix_checks import rectangles
from run_hook_prefix_uncertainty_checks import (COUNTS as UNIT_COUNTS, scalar_interval,
                                               tail_rank, unit_columns)
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import area2, dot

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def prefix_data(data, columns, co, t, a, arm, p):
    index = {v: i for i, v in enumerate(data[0])}
    axis = {y - t: c for (g, y, z), c in zip(data[2], co) if g == 1 and not z} if arm == 'Y' else {
            z - t: c for (g, y, z), c in zip(data[2], co) if g == 1 and not y}
    regular, records = [], []
    for r in range(t, a - t + 1):
        upper, plain, rows = rectangles(t, a, r, arm)
        d, s = r - t, a - t - r
        determinant = bareiss([[axis.get(d + i - j, 0) for j in range(s)] for i in range(s)]) % p
        rank = len(field_basis((*upper, *plain), columns.__getitem__, p)[0])
        assert (rank == t * (a - t)) == bool(determinant)
        if len(plain) <= 32:
            matrix = [[columns[v].get(index[row], 0) for v in plain] for row in rows]
            assert bareiss(matrix) % p == pow(determinant, t, p)
            COUNTS['independent_full_prefix_vs_axis_determinants'] += 1
        if determinant:
            regular.append(r)
        records.append({'r': r, 'axis_determinant_mod_p': determinant, 'actual_prefix_basis_rank': rank})
        COUNTS['actual_prefix_regularity_rank_tests'] += 1
    assert regular[0] == t and regular[-1] == a - t
    gap = max((v - u - 1 for u, v in zip(regular, regular[1:])), default=0)
    for r in range(t, a - t + 1):
        following = next(v for v in regular if v >= r)
        assert r <= following <= r + gap
        COUNTS['actual_next_regular_index_distance_controls'] += 1
    return regular, gap, records


def field_controls():
    rng = Random(2026102241)
    records = []
    cases = ((2, 12, 3, 5, 101), (3, 20, 10, 9, 101), (2, 20, 13, 14, 101),
             (3, 13, 0, 3, 101), (2, 8, 2, 4, 11),
             (2, 12, None, None, 2), (3, 13, None, None, 3))
    for t, a, kY, kZ, p in cases:
        data, L = model(t, a), a - 2 * t
        co = [rng.randrange(p) for _ in data[2]]
        if kY is not None:
            assert p > 2 * L
            for i, (g, y, z) in enumerate(data[2]):
                if g == 1 and not z:
                    co[i] = comb(kY, y - t) % p if t <= y <= t + kY else 0
                elif g == 1 and not y:
                    co[i] = comb(kZ, z - t) % p if t <= z <= t + kZ else 0
        for parameter in ((1, t, 0), (1, 0, t)):
            co[data[2].index(parameter)] = 1
        columns = actual(data, co, p, replay=True)
        arms = []
        for arm, exponent in (('Y', kY), ('Z', kZ)):
            regular, gap, prefixes = prefix_data(data, columns, co, t, a, arm, p)
            if exponent is not None:
                expected = sorted(set(range(t, t + exponent + 1)) | {a - t})
                assert regular == expected and gap == max(0, L - exponent - 1)
                COUNTS['independent_binomial_axis_gap_predictions'] += 1
                for prefix in prefixes:
                    d, s = prefix['r'] - t, a - t - prefix['r']
                    value = Q(1)
                    if d > exponent and s:
                        value = Q(0)
                    elif s:
                        for j in range(s):
                            value *= Q(factorial(exponent + j) * factorial(j),
                                       factorial(d + j) * factorial(exponent - d + j))
                    assert value.denominator == 1 and int(value) % p == prefix['axis_determinant_mod_p']
                    COUNTS['independent_general_binomial_axis_products'] += 1
            rows, values = unit_columns(data, columns, co, t, a, arm, p)
            coordinate, lower = (0 if arm == 'Y' else 1), L - 1 - gap
            for D in range(max(0, lower)):
                source = tuple(h for h in rows if h[coordinate] <= D)
                for E in range(lower - D):
                    assert tail_rank(rows, values, source, lambda v: v[coordinate] > E, p) == len(source)
                    COUNTS['all_vector_gap_degree_uncertainty_tail_ranks'] += 1
            if exponent is not None and exponent <= L - 1:
                h = (0, t - 1) if arm == 'Y' else (t - 1, 0)
                output_degree = max(rows[i][coordinate] for i in values[h])
                assert output_degree == exponent == lower
                COUNTS['sharp_gap_loss_original_vector_witnesses'] += 1
            arms.append({'arm': arm, 'regular_indices': regular, 'maximum_gap': gap,
                         'degree_sum_lower_bound': lower, 'prefixes': prefixes})
        records.append({'t': t, 'a': a, 'length': 2 * len(data[0]), 'field_characteristic': p,
                        'axis_exponents': [kY, kZ], 'coefficients': co, 'arms': arms,
                        'effective_width': t + sum(arm['maximum_gap'] for arm in arms)})
    return records


def scalar_controls():
    rng = Random(2026102243)
    records = []
    for t, hY, hZ in ((2, 1, 2), (3, 5, 4), (4, 10, 2)):
        T = t + hY + hZ
        a, A = 512 * T + 1, 512 * T + 1 - t
        feasible, total = 0, 0
        for trial in range(2500):
            b = rng.randrange(a // 4, 3 * a // 4)
            c = a - b + rng.randrange(-20 * T, 20 * T + 1) if trial % 2 else rng.randrange(a // 4, 3 * a // 4)
            R, Z = rng.randrange(t), rng.randrange(t, 2 * t)
            Y, S = b + R, c + Z
            if not (0 < c and Y < a and S < a and b * c > 64 * a * T):
                continue
            total += 1
            K, J = Y + R, S + Z
            conditions = ((b, -c, True), (-Y, a - Z, True), (a - Y, -Z, True),
                          (K - a + 1 + hY, J - t, True), (K - t, J - a + 1 + hZ, True),
                          (a + 2 * t - 1 + hY - K, 3 * t - 2 - J, False),
                          (3 * t - 2 - K, a + 2 * t - 1 + hZ - J, False))
            interval = scalar_interval(conditions)
            COUNTS['independent_gap_scalar_interval_eliminations'] += 1
            if interval is None:
                continue
            lo, hi = interval
            assert Q(a - 3 * t + 1 - hZ, a + t - 1 + hY) <= lo < hi <= Q(a + t - 1 + hZ, a - 3 * t + 1 - hY)
            u = (lo + hi) / 2
            assert Q(a - 3 * T, a + T) < u < Q(a + T, a - 3 * T) < Q(65, 61)
            F = K * u + J
            assert a * abs(u - 1) < 5 * T and a - T - 1 < F < a + 6 * T
            assert abs(F - b - c) < 12 * T and abs(b + c - a) < 20 * T
            beta = 2 * t - 1 - Z
            rY = a + t - Y + int(Q(beta) // u)
            rZ = int((a + t - 1 - Z - (Y - 2 * t + 1) * u) // 1) + 1
            assert t <= rY < A and t <= rZ < A
            # Moving to any regular index within the gap is bounded by
            # these worst upper-prefix costs; the plain rectangles shrink.
            nextY, nextZ = min(A, rY + hY), min(A, rZ + hZ)
            tau, lam = Y * u + Z, b * u - c
            gapY = tau - (lam + (t + nextY - 1) * u + t - 1)
            gapZ = tau - (lam + (t - 1) * u + t + nextZ - 1)
            assert gapY >= F - a * u - (2 * t - 1 + hY) * u - 3 * t + 2
            assert gapZ >= F - a - (3 * t - 2) * u - 2 * t + 1 - hZ
            feasible += 1
            COUNTS['nonvacuous_gap_scalar_diagonal_strip_controls'] += 1
        assert total and feasible
        records.append({'t': t, 'a': a, 'hY': hY, 'hZ': hZ, 'effective_width': T,
                        'sampled_high_product_pairs': total, 'feasible_necessary_intervals': feasible,
                        'scope': 'Exact rational necessary-inequality controls, not actual coherent wall assertions.'})
    return records


def five_cone_controls():
    records = []
    for t, a in ((2, 5), (2, 12), (3, 10), (4, 17)):
        data = model(t, a)
        co = [int(parameter in ((1, t, 0), (1, 0, t))) for parameter in data[2]]
        for p in (0, 2, 101):
            columns = actual(data, co, p, replay=True)
            weights = ((Q(t, 2), Q(1), Q(1)), (Q(3 * t, 4), Q(1, 2), Q(1)),
                       (Q(2 * t), Q(1, 2), Q(1)), (Q(3 * t, 2), Q(2), Q(1)),
                       (Q(3 * t), Q(2), Q(1)))
            bases, polygons = [], []
            for w in weights:
                B, _ = field_basis(sorted(columns, key=lambda v: (dot(w, v), v)), columns.__getitem__, p)
                N = normals(B, columns, p)
                assert all(dot(w, q) > 0 for q in N)
                poly = polygon(N)
                bases.append(frozenset(B));polygons.append(poly)
                COUNTS['independent_original_full_normal_cones'] += 1
            assert len(set(bases)) == 5 and sum(map(area2, polygons)) == 1
            COUNTS['independent_five_cone_whole_fan_coverages'] += 1
            records.append({'t': t, 'a': a, 'length': 2 * len(data[0]), 'field_characteristic': p,
                            'central_image': 'y^t+z^t', 'outer_images': [0, 0],
                            'maximum_gaps': [a - 2 * t - 1] * 2,
                            'actual_cones': 5, 'sum_normalized_cone_areas': '1',
                            'bases': [sorted(B) for B in bases]})
    return records


def head_cost_controls():
    rng = Random(2026102247)
    records = []
    for t, hY, hZ in ((2, 1, 2), (3, 5, 4), (4, 10, 2)):
        T = t + hY + hZ
        a, A = 512 * T + 1, 512 * T + 1 - t
        local = Counter()
        for trial in range(150):
            Y, Z, R = rng.randrange(a + t, 2 * a - t), rng.randrange(t), rng.randrange(t)
            S = rng.randrange(8 * T + Z + 1, a)
            b, c = Y - R, S - Z
            for u in (Q(c, b) + Q(1, 10007), Q(3) + Q(1, 10007)):
                tau = Y * u + Z
                if tau > (2 * t - 1) * u + a - 1:
                    local['outer_opposite_tail_cost_controls'] += 1
                else:
                    assert u < 2
                    r0 = int(Q(2 * t) // u)
                    r = min(A, r0 + hY)
                    assert t <= r0 <= r <= A
                    assert tau - ((a + t - r - 1) * u + 2 * t - 1) > 0
                    assert S - (t + r - 1 - R) * u - t + 1 > 0
                    local['outer_shifted_prefix_cost_controls'] += 1
            Y, Z, R = rng.randrange(a - 2 * t - 2, a), rng.randrange(t, 2 * t), rng.randrange(t)
            S = rng.randrange(a // 2, a)
            b, c = Y - R, S - Z
            assert b * c > 64 * a * T
            for u in (Q(c, b) + Q(1, 10007), Q(a, Y) + Q(1, 10007), Q(3) + Q(1, 10007)):
                tau = Y * u + Z
                if tau < min(a * u, a):
                    local['central_below_both_heads_cost_controls'] += 1
                elif tau > (2 * t - 1) * u + a - 1:
                    local['central_opposite_tail_cost_controls'] += 1
                elif tau >= a * u:
                    assert u < 2
                    r = t + int(Q(2 * t) // u) + hY
                    assert t <= r < A
                    assert tau - ((a + t - r - 1) * u + 2 * t - 1) > 0
                    assert S - (t + r - 1 - R) * u - t + 1 > 0
                    local['central_own_head_shifted_prefix_cost_controls'] += 1
                else:
                    assert tau >= a
                    r = t + int((2 * t * u) // 1) + hZ
                    assert t <= r < A
                    assert tau - ((2 * t - 1) * u + a + t - r - 1) > 0
                    assert S - (t - 1 - R) * u - t - r + 1 > 0
                    local['central_opposite_head_shifted_prefix_cost_controls'] += 1
        COUNTS.update(local)
        records.append({'t': t, 'a': a, 'hY': hY, 'hZ': hZ, 'counts': dict(sorted(local.items())),
                        'scope': 'Exact rational cost bounds at the worst allowed next-prefix shift; regular-basis ranks are checked separately in original-field controls.'})
    for name in ('outer_shifted_prefix_cost_controls', 'central_below_both_heads_cost_controls',
                 'central_own_head_shifted_prefix_cost_controls',
                 'central_opposite_head_shifted_prefix_cost_controls'):
        assert COUNTS[name], name
    return records


def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    fields = field_controls()
    scalars = scalar_controls()
    head_costs = head_cost_controls()
    fans = five_cone_controls()
    files = {'docs/hook_prefix_gaps.md', 'scripts/run_hook_prefix_gap_checks.py',
             'docs/hook_prefix_uncertainty.md', 'docs/hook_central_prefix.md',
             'docs/hook_central_skeleton.md', 'docs/hook_bridge_reduction.md',
             'docs/cross_layer_direction_budget.md', 'docs/layer_shadow_budget.md'}
    pending = ['scripts/run_hook_prefix_gap_checks.py']
    while pending:
        filename = pending.pop()
        for node in ast.walk(ast.parse((ROOT / filename).read_text())):
            names = ([node.module] if isinstance(node, ast.ImportFrom) else
                     [item.name for item in node.names] if isinstance(node, ast.Import) else [])
            for name in names:
                target = 'scripts/' + str(name) + '.py'
                if str(name).startswith('run_') and (ROOT / target).exists() and target not in files:
                    files.add(target);pending.append(target)
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_hook_prefix_gap_checks.py',
              'scope': 'Uniform gap-parameter direction and fan bound for every compatible hook with two leading central units; arbitrary singular prefix minors and outer images are retained. Long gaps and vanished leading coefficients remain unresolved at the original target scale.',
              'counts': dict(sorted(COUNTS.items())), 'imported_unit_column_counts': dict(sorted(UNIT_COUNTS.items())),
              'original_quotient_replay_counts': dict(sorted(REPLAY_COUNTS.items())),
              'actual_prefix_gap_controls': fields, 'scalar_controls': scalars,
              'head_and_outer_cost_controls': head_costs,
              'longest_gap_five_cone_controls': fans,
              'limits': ['The parameter bound does not prove a near-quadratic count when gaps are large.',
                         'Maximal gaps can coexist with exactly five cones; they do not establish a lower separation.',
                         'The scalar controls are not actual-field wall constructions.',
                         'Vanished leading coefficients, unrestricted sharp growth, exact priority and the original final-paper goal remain open.'],
              'source_sha256': {filename: hashlib.sha256((ROOT / filename).read_bytes()).hexdigest() for filename in sorted(files)}}
    (ROOT / 'results/research/hook_prefix_gap_checks.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'counts': report['counts'], 'source_hashes': len(files)}))


if __name__ == '__main__':
    main()
