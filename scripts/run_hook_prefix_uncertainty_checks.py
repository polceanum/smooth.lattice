#!/usr/bin/env python3
"""All-vector tail-rank and exact scalar checks for prefix-regular hook counts."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_hook_bridge_checks import COUNTS as REPLAY_COUNTS, actual, model
from run_hook_central_prefix_checks import (COUNTS as PREFIX_COUNTS,
                                           explicit_coefficients, prefix_ranks,
                                           rectangles)
from run_signed_face_checks import field_basis

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def unit_columns(data, columns, coefficients, t, a, arm, p):
    A = a - t
    rows = tuple(product(range(A), range(t))) if arm == 'Y' else tuple(product(range(t), range(A)))
    index = {v: i for i, v in enumerate(rows)}
    unit = {(y - t, z): c for (g, y, z), c in zip(data[2], coefficients)
            if g == 1 and y >= t and c} if arm == 'Y' else {
            (y, z - t): c for (g, y, z), c in zip(data[2], coefficients)
            if g == 1 and z >= t and c}
    values = {}
    for h in rows:
        literal = {index[v]: c % p for shift, c in unit.items()
                   if (v := (h[0] + shift[0], h[1] + shift[1])) in index and c % p}
        label = (0, t + h[0], t + h[1])
        original = {}
        for r, c in columns[label].items():
            y, z = data[0][r]
            if arm == 'Y' and y >= t:
                original[index[y - t, z]] = c
            elif arm == 'Z' and z >= t:
                original[index[y, z - t]] = c
        assert literal == original
        values[h] = literal
        COUNTS['independent_original_rectangular_unit_columns'] += 1
    return rows, values


def tail_rank(rows, values, source, tail, p):
    allowed = {i for i, row in enumerate(rows) if tail(row)}
    restricted = {h: {i: c for i, c in values[h].items() if i in allowed} for h in source}
    return len(field_basis(source, restricted.__getitem__, p)[0])


def degree_controls():
    rng = Random(2026102237)
    records = []
    for t, a, p in ((2, 8, 11), (2, 12, 101), (3, 13, 101), (2, 20, 101), (3, 20, 101)):
        data, L = model(t, a), a - 2 * t
        coefficients = explicit_coefficients(data, t, a, rng, p)
        # Axis polynomials stay binomial; every off-axis central coefficient
        # varies independently. The original full prefix matrices are replayed.
        for i, (g, y, z) in enumerate(data[2]):
            if g == 1 and y and z:
                coefficients[i] = rng.randrange(p)
            elif g != 1:
                coefficients[i] = 0
        columns = prefix_ranks(data, coefficients, t, a, p)
        tested = []
        for arm in ('Y', 'Z'):
            rows, values = unit_columns(data, columns, coefficients, t, a, arm, p)
            coordinate = 0 if arm == 'Y' else 1
            for D in range(L - 1):
                source = tuple(h for h in rows if h[coordinate] <= D)
                for E in range(L - 1 - D):
                    assert tail_rank(rows, values, source, lambda v: v[coordinate] > E, p) == len(source)
                    COUNTS['all_vector_degree_uncertainty_tail_ranks'] += 1
                    tested.append((D, E))
            for u in (Q(1, 3), Q(1), Q(7, 3)):
                weight = lambda v: u * v[0] + v[1]
                threshold = u * (L - 1) if arm == 'Y' else Q(L - 1)
                for d in sorted({weight(h) for h in rows if weight(h) < threshold}):
                    source = tuple(h for h in rows if weight(h) <= d)
                    e = threshold - d - Q(1, 10007)
                    assert tail_rank(rows, values, source, lambda v: weight(v) > e, p) == len(source)
                    COUNTS['all_vector_weighted_uncertainty_tail_ranks'] += 1
        records.append({'t': t, 'a': a, 'field_characteristic': p, 'length': 2 * len(data[0]),
                        'central_coefficients': [[y, z, c] for (g, y, z), c in zip(data[2], coefficients)
                                                 if g == 1 and c],
                        'degree_sum_lower_bound': L - 1,
                        'degree_subspace_tests_both_arms': len(tested),
                        'scope': 'Each tail rank excludes every nonzero vector in the full tested source subspace.'})
    return records


def excluded_stratum():
    t, a, p = 2, 12, 101
    data = model(t, a)
    co = [int(parameter in ((1, t, 0), (1, 0, t))) for parameter in data[2]]
    columns = actual(data, co, p, replay=True)
    index = {v: i for i, v in enumerate(data[0])}
    for arm in ('Y', 'Z'):
        rows, values = unit_columns(data, columns, co, t, a, arm, p)
        coordinate = 0 if arm == 'Y' else 1
        source = tuple(h for h in rows if h[coordinate] == 0)
        assert tail_rank(rows, values, source, lambda v: v[coordinate] > 0, p) == 0
        _, plain, target = rectangles(t, a, t + 1, arm)
        idx = {index[v] for v in target}
        rank = len(field_basis(plain, lambda v: {i: c for i, c in columns[v].items() if i in idx}, p)[0])
        assert rank < len(plain)
        COUNTS['retained_unit_only_degree_uncertainty_failures'] += 1
    return {'t': t, 'a': a, 'field_characteristic': p,
            'central_image': 'y^t+z^t', 'witness': 'H=1 has both degrees zero',
            'required_degree_sum': a - 2 * t - 1,
            'scope': 'Leading units alone do not justify the new family bound.'}


def scalar_interval(inequalities):
    lo, hi = Q(0), None
    for alpha, beta, strict in inequalities:
        # Only nonempty open intervals are relevant to generic wall points.
        # For nonconstant inequalities, strictness changes their endpoints
        # but not their interior. Constant zero needs separate treatment.
        if not alpha:
            if beta < 0 or (strict and not beta):
                return None
            continue
        bound = Q(-beta, alpha)
        if alpha > 0:
            lo = max(lo, bound)
        else:
            hi = bound if hi is None else min(hi, bound)
    return None if hi is None or lo >= hi else (lo, hi)


def scalar_controls():
    records = []
    for t, a in ((2, 80), (2, 160), (3, 128), (3, 192)):
        candidates, feasible, maximum_error = 0, 0, 0
        sample = None
        rho = Q(a + t - 1, a - 3 * t + 1)
        for b, c in product(range(1, a), repeat=2):
            if b * c <= 16 * a * t:
                continue
            for R, Z in product(range(t), range(t, 2 * t)):
                Y, S = b + R, c + Z
                if Y >= a or S >= a:
                    continue
                candidates += 1
                K, J = Y + R, S + Z
                inequalities = ((b, -c, True), (-Y, a - Z, True), (a - Y, -Z, True),
                                (K - a + 1, J - t, True), (K - t, J - a + 1, True),
                                (a + 2 * t - 1 - K, 3 * t - 2 - J, False),
                                (3 * t - 2 - K, a + 2 * t - 1 - J, False))
                interval = scalar_interval(inequalities)
                COUNTS['independent_exact_scalar_interval_eliminations'] += 1
                if interval is None:
                    continue
                lo, hi = interval
                assert Q(1, 1) / rho <= lo < hi <= rho
                u = (lo + hi) / 2
                F = K * u + J
                assert max((a - 1) * u + t, a + t * u - 1) < F
                assert F <= min((a + 2 * t - 1) * u + 3 * t - 2,
                                a + (3 * t - 2) * u + 2 * t - 1)
                beta = 2 * t - 1 - Z
                rY = a + t - Y + int(Q(beta) // u)
                rZ = int((a + t - 1 - Z - (Y - 2 * t + 1) * u) // 1) + 1
                assert t <= rY < a - t and t <= rZ < a - t
                assert (Y - a - t + rY + 1) * u - beta > 0
                assert Y * u + Z - ((2 * t - 1) * u + a + t - rZ - 1) > 0
                assert a * abs(u - 1) < 5 * t and a - 1 < F < a + 6 * t
                assert abs(F - b - c) < 12 * t and abs(b + c - a) < 20 * t
                maximum_error = max(maximum_error, abs(b + c - a))
                feasible += 1
                COUNTS['nonvacuous_scalar_diagonal_strip_controls'] += 1
                sample = sample or {'b': b, 'c': c, 'R': R, 'Z': Z, 'u': str(u)}
        assert candidates
        records.append({'t': t, 'a': a, 'candidate_literal_pairs': candidates,
                        'feasible_necessary_intervals': feasible, 'maximum_abs_sum_minus_a': maximum_error,
                        'sample': sample,
                        'scope': 'Complete elimination of the necessary scalar inequalities at these sizes; feasible intervals are not actual-field wall claims.'})
    return records


def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    fields = degree_controls()
    excluded = excluded_stratum()
    scalars = scalar_controls()
    assert COUNTS['nonvacuous_scalar_diagonal_strip_controls']
    files = {'docs/hook_prefix_uncertainty.md', 'scripts/run_hook_prefix_uncertainty_checks.py',
             'docs/hook_central_prefix.md', 'docs/hook_central_skeleton.md',
             'docs/hook_bridge_reduction.md', 'docs/cross_layer_direction_budget.md',
             'docs/layer_shadow_budget.md', 'docs/bounded_coordinate_wall_shadows.md',
             'docs/coordinate_slice_wall_budget.md'}
    pending = ['scripts/run_hook_prefix_uncertainty_checks.py']
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
              'command': 'python3 -B scripts/run_hook_prefix_uncertainty_checks.py',
              'scope': 'Uniform O(n^2 log^2(n+1)) fan bound and cubic polylogarithmic collective wall budgets on the entire central prefix-regular hook locus, with arbitrary outer images. The prefix-degenerate strata and unrestricted endpoint remain open.',
              'counts': dict(sorted(COUNTS.items())),
              'prefix_check_counts': dict(sorted(PREFIX_COUNTS.items())),
              'original_quotient_replay_counts': dict(sorted(REPLAY_COUNTS.items())),
              'all_vector_degree_controls': fields, 'excluded_stratum_control': excluded,
              'complete_scalar_controls': scalars,
              'limits': ['No uniform bound for prefix-degenerate central images is proved.',
                         'Feasible scalar intervals and finite inactive queries are not actual coherent wall constructions.',
                         'No unrestricted sharp growth, infinite polynomial separation or exact current novelty is established.',
                         'The original goal and its final-paper criteria remain incomplete.'],
              'source_sha256': {filename: hashlib.sha256((ROOT / filename).read_bytes()).hexdigest()
                                for filename in sorted(files)}}
    (ROOT / 'results/research/hook_prefix_uncertainty_checks.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'counts': report['counts'], 'source_hashes': len(files)}))


if __name__ == '__main__':
    main()
