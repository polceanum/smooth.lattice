#!/usr/bin/env python3
"""Original-field checks for the large-product central-only wall reduction."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_hook_bridge_checks import COUNTS as REPLAY_COUNTS, actual, model
from run_hook_central_prefix_checks import (COUNTS as PREFIX_COUNTS, eligibility,
                                           explicit_coefficients, prefix_ranks,
                                           rectangles, transpose)
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import dot

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def inventory_controls():
    records = []
    for t, a in ((2, 4), (2, 5), (2, 13), (3, 17), (2, 48),
                 (3, 101), (4, 257), (32, 1025), (64, 4097)):
        bound, budget = 2 * a - 1, 16 * a * t
        positive = sum(min(bound, budget // b) for b in range(1, bound + 1))
        count = 4 * bound + 1 + 4 * positive
        if a <= 48:
            brute = sum(abs(d * e) <= budget for d, e in
                        product(range(-bound, bound + 1), repeat=2))
            assert count == brute
            COUNTS['independent_signed_inventory_enumerations'] += 1
        dyadic_ceiling = 4 * bound + 1 + 4 * budget * (1 + bound.bit_length())
        assert count <= dyadic_ceiling
        assert 2 * a * t <= 4 * a * t - 2 * t * t <= 4 * a * t
        COUNTS['hyperbolic_inventory_uniform_bound_controls'] += 1
        records.append({'t': t, 'a': a, 'length': 4 * a * t - 2 * t * t,
                        'signed_hyperbolic_directions': count,
                        'integer_dyadic_ceiling': dyadic_ceiling,
                        'previous_near_inventory_additional_bound': 256 * a * t})
    return records


def multiply(first, second, width, height, p):
    result = Counter()
    for (i, j), c in first.items():
        for (r, s), d in second.items():
            if i + r < width and j + s < height:
                result[i + r, j + s] += c * d
    return {v: c % p for v, c in result.items() if c % p}


def inverse_unit(unit, width, height, p):
    constant = pow(unit[0, 0], -1, p)
    inverse = {(0, 0): constant}
    for i, j in sorted(product(range(width), range(height)), key=lambda v: (sum(v), v)):
        if (i, j) == (0, 0):
            continue
        value = -constant * sum(c * inverse.get((i - r, j - s), 0)
                                for (r, s), c in unit.items() if (r, s) != (0, 0)
                                and r <= i and s <= j) % p
        if value:
            inverse[i, j] = value
    assert multiply(unit, inverse, width, height, p) == {(0, 0): 1}
    return inverse


def module_controls():
    rng = Random(2026102229)
    records = []
    for t, a, p in ((2, 5, 2), (3, 10, 3), (2, 13, 101), (3, 17, 101)):
        data, A = model(t, a), a - t
        co = [rng.randrange(p) for _ in data[2]]
        for parameter in ((1, t, 0), (1, 0, t)):
            co[data[2].index(parameter)] = 1
        columns = actual(data, co, p, replay=True)
        f1 = {(i, j): c for (g, i, j), c in zip(data[2], co) if g == 1 and c}
        U = {(i - t, j): c for (i, j), c in f1.items() if i >= t}
        V = {(i, j - t): c for (i, j), c in f1.items() if j >= t}
        invU, invV = inverse_unit(U, A, t, p), inverse_unit(V, t, A, p)

        def normalize(column):
            Y, Z = {}, {}
            for row, c in column.items():
                i, j = data[0][row]
                if i >= t:
                    Y[i - t, j] = c
                elif j >= t:
                    Z[i, j - t] = c
                else:
                    assert False, 'The residual vectors have no common-box coordinates.'
            return multiply(Y, invU, A, t, p), multiply(Z, invV, t, A, p)

        sources = tuple((i, j) for i, j in product(range(A), repeat=2) if i < t or j < t)
        central = tuple((0, t + i, t + j) for i, j in sources)
        rank = len(field_basis(central, columns.__getitem__, p)[0])
        assert rank == len(sources) == 2 * t * A - t * t
        COUNTS['actual_entire_central_column_independence'] += 1
        for h, label in zip(sources, central):
            Y, Z = normalize(columns[label])
            assert Y == ({h: 1} if h[0] < A and h[1] < t else {})
            assert Z == ({h: 1} if h[0] < t and h[1] < A else {})
            COUNTS['independent_original_inverse_unit_column_identities'] += 1
        for shift in ((A, 0), (t, t), (0, A)):
            literal = {(i + shift[0], j + shift[1]): c for (i, j), c in f1.items()
                       if (i + shift[0], j + shift[1]) in set(data[0])}
            assert not literal
            COUNTS['independent_monomial_annihilator_generator_checks'] += 1
        box = tuple(product(range(t), repeat=2))
        box_index = {v: i for i, v in enumerate(box)}
        quotient = {}
        for label in data[4][0] + data[4][1]:
            if not label[0]:
                continue
            Y, Z = normalize(columns[label])
            quotient[label] = {box_index[v]: (Y.get(v, 0) - Z.get(v, 0)) % p
                               for v in box if (Y.get(v, 0) - Z.get(v, 0)) % p}
        assert len(field_basis(tuple(quotient), quotient.__getitem__, p)[0]) == t * t
        assert 2 * t * A - rank == t * t
        COUNTS['actual_common_box_cokernel_dimensions'] += 1
        records.append({'t': t, 'a': a, 'field_characteristic': p,
                        'central_coefficients': [[i, j, c] for (i, j), c in sorted(f1.items())],
                        'annihilator_generators': [[A, 0], [t, t], [0, A]],
                        'central_column_rank': rank, 'cokernel_dimension': t * t,
                        'scope': 'Only the two leading units are assumed; other prefix minors may vanish.'})
    return records


def certificate(columns, t, a, e, f, u, arm, p):
    Y, Z, R, S = e[1], e[2], f[1], f[2]
    b, c = Y - R, S - Z
    assert 2 * t <= Y < a and t <= Z < 2 * t and 0 <= R < t and t <= S < a
    assert b * c > 16 * a * t and b > 16 * t and c > 16 * t
    w = (b * u - c, u, Q(1))
    tau = Y * u + Z
    assert min(w) > 0 and tau == dot(w, f)
    if tau < min(a * u, a):
        return None, w, e, f, 'below_both_outer_heads'
    if tau > (2 * t - 1) * u + a - 1:
        chosen = tuple((0, i, j) for i, j in product(range(t, 2 * t), range(2 * t, a)))
        target, rank = f, t * (a - 2 * t)
        branch = 'opposite_tail_strict_closure'
    elif tau >= a:
        assert u <= Q(a - 1 - Z, Y - 2 * t + 1) < Q(a, b - 2 * t) < Q(4 * a, 3 * b)
        r = t + int((2 * t * u) // 1)
        assert t <= r < a - t and Q(3 * t) * u < Q(c, 4)
        upper, plain, _ = rectangles(t, a, r, 'Z')
        chosen, target, rank = (*upper, *plain), f, t * (a - t)
        assert r + 1 - t - (2 * t - 1) * u > 0
        assert S - 3 * t * u - 2 * t + 1 > 0
        branch = 'opposite_prefix_strict_closure'
    else:
        # CP14 supplies the already established own-head exclusion.
        assert tau >= a * u and u < 2
        r = t + int(Q(2 * t) // u)
        assert t <= r < a - t
        upper, plain, _ = rectangles(t, a, r, 'Y')
        chosen, target, rank = (*upper, *plain), e, t * (a - t)
        branch = 'own_prefix_strict_closure'
    if arm == 'Z':
        chosen = tuple(map(transpose, chosen))
        target = transpose(target)
    test_w = transpose(w) if arm == 'Z' else w
    assert all(dot(test_w, v) < tau for v in chosen)
    basis, pivots = field_basis(chosen, columns.__getitem__, p)
    assert len(basis) == rank
    remainder, _ = field_reduce(columns[target], pivots, p)
    assert not remainder
    COUNTS[branch] += 1
    COUNTS['independent_original_strict_closure_certificates'] += 1
    return chosen, w, e, f, branch


def field_controls():
    rng = Random(2026102227)
    records = []
    configurations = [
        (2, 48, 101, ((42, 2, 0, 42), (43, 2, 0, 42), (42, 2, 1, 43),
                      (44, 3, 0, 45), (47, 2, 0, 47))),
        (3, 64, 257, ((59, 3, 0, 61), (60, 3, 0, 60), (60, 3, 1, 61),
                      (61, 4, 0, 63), (63, 3, 0, 63)))]
    for t, a, p, pairs in configurations:
        data = model(t, a)
        coefficients = explicit_coefficients(data, t, a, rng, p)
        changed = [rng.randrange(p) if parameter[0] != 1 else value
                   for parameter, value in zip(data[2], coefficients)]
        central = [value if parameter[0] == 1 else 0
                   for parameter, value in zip(data[2], coefficients)]
        first = prefix_ranks(data, coefficients, t, a, p)
        second = actual(data, changed, p, replay=True)
        only_central = actual(data, central, p, replay=True)
        maps = first, second, only_central
        checks = []
        for Y, Z, R, S in pairs:
            b, c = Y - R, S - Z
            e, f = (0, Y, Z), (1, R, S)
            for u in (Q(c, b) + Q(1, 10007), Q(a, Y) + Q(1, 10007), Q(3) + Q(1, 10007)):
                for arm in ('Y', 'Z'):
                    chosen, w, E, F, branch = certificate(first, t, a, e, f, u, arm, p)
                    if arm == 'Z':
                        w, E, F = transpose(w), transpose(E), transpose(F)
                    results = [eligibility(columns, E, F, w, p) for columns in maps]
                    assert len({value for value, _ in results}) == 1
                    assert results[0][1] == results[1][1] == results[2][1]
                    if chosen is None:
                        cheap = results[0][1]
                        assert all(v[0] or (v[1] < a and v[2] < a) for v in cheap)
                        assert all(first[v] == second[v] == only_central[v] for v in (*cheap, E, F))
                        COUNTS['both_outer_images_absent_from_cheap_rank_test'] += 1
                    else:
                        assert not results[0][0]
                        assert all(first[v] == second[v] == only_central[v] for v in chosen)
                    COUNTS['actual_two_outer_image_eligibility_invariance'] += 1
                    direction = tuple(x - y for x, y in zip(F, E))
                    outside = direction not in data[7] and abs(direction[1] * direction[2]) > 16 * a * t
                    COUNTS['actual_pair_tests_outside_known_inventory'] += outside
                    checks.append({'arm': arm, 'e': E, 'f': F, 'u_before_transposition': str(u),
                                   'direction': direction, 'outside_known_inventory': outside,
                                   'certificate': branch, 'eligible_in_all_three_maps': results[0][0]})
        records.append({'t': t, 'a': a, 'length': 2 * len(data[0]), 'field_characteristic': p,
                        'common_central_image': 'y^t*(1+y)^(a-2t)+z^t*(1+z)^(a-2t)',
                        'first_coefficients': coefficients, 'changed_coefficients': changed,
                        'central_only_coefficients': central, 'cases': checks,
                        'eligible_sampled_pairs': sum(case['eligible_in_all_three_maps'] for case in checks)})
    return records


def controlled_inventory_guard():
    t, a, p = 3, 10, 1000003
    data = model(t, a)
    saved = json.loads((ROOT / 'results/research/generic_hook_checks.json').read_text())
    co = saved['families'][1]['independent_reference_coefficients']
    central = [value if parameter[0] == 1 else 0 for parameter, value in zip(data[2], co)]
    first, second = actual(data, co, p, replay=True), actual(data, central, p, replay=True)
    index = {v: i for i, v in enumerate(data[0])}
    for arm in ('Y', 'Z'):
        for r in range(t, a - t + 1):
            _, plain, rows = rectangles(t, a, r, arm)
            idx = {index[v] for v in rows}
            assert len(field_basis(plain, lambda v: {i: c for i, c in first[v].items()
                                                   if i in idx}, p)[0]) == len(plain)
            COUNTS['retained_guard_actual_prefix_regularity_minors'] += 1
    e, f, w = (0, 7, 3), (1, 1, 7), (Q(15, 14), Q(71, 84), Q(1))
    assert eligibility(first, e, f, w, p)[0]
    assert not eligibility(second, e, f, w, p)[0]
    direction = tuple(x - y for x, y in zip(f, e))
    assert abs(direction[1] * direction[2]) <= 16 * a * t
    COUNTS['positive_literal_pair_rank_control'] += 1
    COUNTS['controlled_inventory_needed_guard'] += 1
    return {'t': t, 'a': a, 'length': 102, 'field_characteristic': p,
            'coefficient_source': 'results/research/generic_hook_checks.json:families[1].independent_reference_coefficients',
            'prefix_regular_by_actual_ranks': True, 'lower': e, 'upper': f,
            'weight': list(map(str, w)), 'direction': direction,
            'eligible_full_extension': True, 'eligible_central_only': False,
            'scope': 'The controlled direction inventory cannot be omitted; no growing-family claim.'}


def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    inventories = inventory_controls()
    modules = module_controls()
    fields = field_controls()
    guard = controlled_inventory_guard()
    for name in ('opposite_tail_strict_closure', 'opposite_prefix_strict_closure',
                 'own_prefix_strict_closure', 'both_outer_images_absent_from_cheap_rank_test',
                 'actual_pair_tests_outside_known_inventory'):
        assert COUNTS[name], name
    files = {'docs/hook_central_skeleton.md', 'scripts/run_hook_central_skeleton_checks.py',
             'docs/hook_central_prefix.md', 'docs/hook_bridge_reduction.md',
             'docs/layer_shadow_budget.md', 'docs/cross_layer_direction_budget.md',
             'docs/hook_wall_recurrence.md', 'results/research/generic_hook_checks.json'}
    pending = ['scripts/run_hook_central_skeleton_checks.py']
    while pending:
        filename = pending.pop()
        for node in ast.walk(ast.parse((ROOT / filename).read_text())):
            names = ([node.module] if isinstance(node, ast.ImportFrom) else
                     [item.name for item in node.names] if isinstance(node, ast.Import) else [])
            for name in names:
                target = 'scripts/' + str(name) + '.py'
                if str(name).startswith('run_') and (ROOT / target).exists() and target not in files:
                    files.add(target)
                    pending.append(target)
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_hook_central_skeleton_checks.py',
              'scope': 'Uniform large-product literal exchange and geometric wall-locus reduction to the central-only extension on the prefix-regular locus. Distinct direction inventories differ only within an O(n log(n+1)) set.',
              'counts': dict(sorted(COUNTS.items())),
              'prefix_check_counts': dict(sorted(PREFIX_COUNTS.items())),
              'original_quotient_replay_counts': dict(sorted(REPLAY_COUNTS.items())),
              'inventory_controls': inventories, 'actual_field_controls': fields,
              'central_module_controls': modules, 'controlled_inventory_guard': guard,
              'limits': ['Geometric wall-locus equality does not identify parallel piece subdivisions or full fan counts.',
                         'Finite inactive-pair controls do not prove uniform central-wall absence.',
                         'The central-only direction count and prefix-degenerate coefficient strata remain open.',
                         'No new unrestricted exponent, actual infinite polynomial separation, current novelty, or goal completion is asserted.'],
              'source_sha256': {filename: hashlib.sha256((ROOT / filename).read_bytes()).hexdigest()
                                for filename in sorted(files)}}
    path = ROOT / 'results/research/hook_central_skeleton_checks.json'
    path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'counts': report['counts'],
                      'source_hashes': len(files)}))


if __name__ == '__main__':
    main()
