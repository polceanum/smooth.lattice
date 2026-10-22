#!/usr/bin/env python3
"""Exact parallel-wall queries and a literal-pair recurrence certificate."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations
from math import gcd
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_generic_hook_checks import (
    Field, bareiss, evaluate_minor, exchange_normals, matching_basis, minor_labels, model as generic_model,
    prime_factors,
)
from run_hook_bridge_checks import COUNTS as BRIDGE_COUNTS, actual, model, normals
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import dot

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def domain(q):
    assert q[0] == 1
    if q[1] > 0:
        return (Q(0), Q(-q[2], q[1])) if q[2] < 0 else None
    if q[1] < 0:
        return max(Q(0), Q(-q[2], q[1])), None
    return (Q(0), None) if q[2] < 0 else None


def moment(B):
    return tuple(sum(a[i] for a in B) for i in range(3))


def parallel_query(columns, q, p):
    """All exact shadow intervals on q-perp in the positive orthant.

    Columns are actual vectors over Q (p=0) or the prime field F_p.
    The two returned bases have minimum and maximum first moment.
    """
    bounds = domain(q)
    if bounds is None:
        return []
    lo, hi = bounds
    derivative = (-q[1], 1, 0)
    rows = []
    while hi is None or lo < hi:
        w = (-q[1] * lo - q[2], lo, Q(1))
        end, bases = hi, []
        for sign in (-1, 1):
            order = sorted(columns, key=lambda a: (dot(w, a), dot(derivative, a), sign * a[0], a))
            B, _ = field_basis(order, columns.__getitem__, p)
            for v in normals(B, columns, p):
                alpha, beta = dot(derivative, v), v[2] - q[2] * v[0]
                value = alpha * lo + beta
                assert value >= 0 and (value or alpha >= 0)
                if alpha < 0:
                    bound = Q(-beta, alpha)
                    assert bound > lo
                    end = bound if end is None else min(end, bound)
            bases.append(frozenset(B))
        mu = list(map(moment, bases))
        delta = tuple(x - y for x, y in zip(*mu))
        assert delta[0] >= 0 and delta == tuple(delta[0] * x for x in q)
        rows.append({'lo': lo, 'hi': end, 'upper': bases[0], 'lower': bases[1], 'multiplicity': delta[0]})
        COUNTS['actual_field_shadow_intervals'] += 1
        if end is None:
            break
        assert end > lo
        lo = end
    return rows


def profile_envelope(bases, q):
    """Independent lower envelope of an explicit list of basis profiles."""
    bounds = domain(q)
    if bounds is None:
        return []
    left, right = bounds
    groups = defaultdict(set)
    for B in bases:
        mu = moment(B)
        groups[(mu[1] - q[1] * mu[0], mu[2] - q[2] * mu[0])].add(mu)
    # Within a slope, only the least intercept can occur on the lower hull.
    best = {}
    for a, b in groups:
        best[a] = min(best.get(a, b), b)
    hull = []
    for a, b in sorted(best.items(), reverse=True):
        start = None
        while hull:
            a0, b0, start0 = hull[-1]
            start = Q(b - b0, a0 - a)
            if start0 is None or start > start0:
                break
            hull.pop()
        hull.append((a, b, start if hull else None))
    out = []
    for i, (a, b, start) in enumerate(hull):
        end = hull[i + 1][2] if i + 1 < len(hull) else None
        lo = left if start is None else max(left, start)
        hi = right if end is None else end if right is None else min(right, end)
        if hi is not None and lo >= hi:
            continue
        mus = groups[(a, b)]
        low, high = min(mus, key=lambda v: v[0]), max(mus, key=lambda v: v[0])
        delta = tuple(x - y for x, y in zip(high, low))
        assert delta == tuple(delta[0] * x for x in q)
        out.append((lo, hi, low, high))
    return out


def canonical(rows):
    return [(r['lo'], r['hi'], moment(r['lower']), moment(r['upper'])) for r in rows]


def brute_controls():
    rng = Random(2026102209)
    for p in (0, 2, 3, 5):
        for trial in range(5):
            labels = rng.sample([(x, y, z) for x in range(2) for y in range(5) for z in range(5)], 8)
            columns = {a: {i: c for i in range(3) if (c := rng.randrange(3))} for a in labels}
            if p:
                columns = {a: {i: c % p for i, c in v.items() if c % p} for a, v in columns.items()}
            rank = len(field_basis(labels, columns.__getitem__, p)[0])
            bases = [frozenset(B) for B in combinations(labels, rank)
                     if len(field_basis(B, columns.__getitem__, p)[0]) == rank]
            for q in ((1, -2, 1), (1, 2, -3), (1, 0, -1)):
                rows = parallel_query(columns, q, p)
                assert canonical(rows) == profile_envelope(bases, q)
                COUNTS['independent_all_basis_envelope_equalities'] += 1
            COUNTS['abstract_actual_rank_models'] += 1


def selected_content_checks(family, indices):
    """Recheck only the coherent minors used by the direction certificate."""
    D, index, _, _, parameters, columns, _ = generic_model(family['thickness'])
    Bs = [frozenset(map(tuple, B)) for B in family['support_fan_bases']]
    minors = {i: minor_labels(Bs[i], index, columns) for i in indices}
    gcds = {i: 0 for i in indices}
    for evaluation in family['integer_evaluations']:
        for i in indices:
            claimed = evaluation['determinants'][i]
            if claimed is None:
                assert gcds[i] == 1
                continue
            value = bareiss(evaluate_minor(minors[i], evaluation['coefficients']))
            assert value == claimed
            gcds[i] = gcd(gcds[i], value)
            COUNTS['independent_selected_integer_minor_evaluations'] += 1
    obligations = set()
    for i, value in gcds.items():
        assert value == family['integer_gcds'][i] and value
        factors = prime_factors(value)
        assert factors == family['prime_divisors'][i]
        obligations.update((i, p) for p in factors)
    for certificate in family['finite_field_certificates']:
        if not any(p == certificate['characteristic'] for _, p in obligations):
            continue
        F = Field(certificate['characteristic'], certificate['degree'], certificate['modulus'])
        for evaluation in certificate['evaluations']:
            for i, claimed in evaluation['nonzero_determinants']:
                if (i, F.p) not in obligations:
                    continue
                value = F.determinant(evaluate_minor(minors[i], evaluation['coefficients']))
                assert value == claimed and value
                obligations.remove((i, F.p))
                COUNTS['independent_selected_content_prime_exclusions'] += 1
    assert not obligations
    COUNTS['all_characteristic_selected_coherent_minors'] += len(indices)


def hook_control(family, q, expected):
    t, a, p = family['thickness'], family['thickness'] ** 2 + 1, 1000003
    data = model(t, a)
    columns = actual(data, family['independent_reference_coefficients'], p, replay=True)
    rows = parallel_query(columns, q, p)
    fixed_upper = data[5]
    Bs = [frozenset(map(tuple, B)) for B in family['support_fan_bases']]
    assert canonical(rows) == profile_envelope([B - fixed_upper for B in Bs], q)
    COUNTS['independent_complete_saved_fan_envelope_equalities'] += 1
    lookup = {B: i for i, B in enumerate(Bs)}
    used = {lookup[r[key] | fixed_upper] for r in rows for key in ('upper', 'lower')}
    selected_content_checks(family, used)
    # Every selected support-greedy basis is independent at the common integer
    # specialization modulo p. Its integer minor is nonzero over Q as well.
    # Thus the same complete envelope is attained by that rational quotient.
    _, _, _, _, _, support, _ = generic_model(t)
    pairs = defaultdict(list)
    record = []
    for r in rows:
        u = (r['lo'] + r['hi']) / 2 if r['hi'] is not None else r['lo'] + 1
        w = (-q[1] * u - q[2], u, Q(1))
        for key, sign in (('upper', -1), ('lower', 1)):
            B, _ = matching_basis(sorted(support, key=lambda a: (dot(w, a), sign * a[0], a)), support, len(data[0]))
            assert B == r[key] | fixed_upper
            for v in exchange_normals(B, support, len(data[0])):
                alpha, beta = -q[1] * v[0] + v[1], v[2] - q[2] * v[0]
                assert alpha * r['lo'] + beta >= 0
                assert alpha >= 0 if r['hi'] is None else alpha * r['hi'] + beta >= 0
                COUNTS['independent_entire_support_interval_inequalities'] += 1
            COUNTS['independent_support_optima_with_actual_minor_survival'] += 1
        low, high = r['lower'] - r['upper'], r['upper'] - r['lower']
        assert len(low) == len(high) == r['multiplicity']
        for e in low:
            f = tuple(x + y for x, y in zip(e, q))
            assert f in high
            ey, ez, _ = data[4]
            assert (e in ey and f in ez) or (e in ez and f in ey)
            COUNTS['literal_pure_cross_arm_pair_checks'] += 1
            cheap = tuple(v for v in columns if dot(w, v) < dot(w, e))
            rank = lambda E: len(field_basis(E, columns.__getitem__, p)[0])
            base = rank(cheap)
            assert rank((*cheap, e)) == rank((*cheap, f)) == rank((*cheap, e, f)) == base + 1
            pairs[(e, f)].append((r['lo'], r['hi']))
            COUNTS['independent_literal_pair_parallel_rank_tests'] += 1
        record.append({'lo': str(r['lo']), 'hi': str(r['hi']) if r['hi'] is not None else None,
                       'multiplicity': r['multiplicity'],
                       'lower_profile': moment(r['lower']), 'upper_profile': moment(r['upper']),
                       'lower_certificate_basis_index': lookup[r['lower'] | fixed_upper],
                       'upper_certificate_basis_index': lookup[r['upper'] | fixed_upper],
                       'exchanged_pairs': [(e, tuple(x + y for x, y in zip(e, q))) for e in sorted(low)]})
    assert sum(bool(r['multiplicity']) for r in rows) == expected
    recurrences = []
    for (e, f), intervals in sorted(pairs.items()):
        components = []
        for lo, hi in intervals:
            if components and components[-1][1] == lo:
                components[-1] = (components[-1][0], hi)
            else:
                components.append((lo, hi))
        recurrences.append({'lower': e, 'upper': f, 'pieces': len(intervals),
                            'components': [[str(lo), str(hi) if hi is not None else None] for lo, hi in components]})
    assert q not in data[7]
    return {'t': t, 'a': a, 'length': 2 * len(data[0]), 'q': q,
            'prime_field_characteristic': p, 'common_integer_coefficients': family['independent_reference_coefficients'],
            'same_complete_query_over_Q': True,
            'generic_query_in_every_characteristic': True,
            'shadow_intervals': len(rows), 'parallel_pieces': expected,
            'all_characteristic_coherent_minor_indices': sorted(used),
            'pair_recurrences': recurrences, 'intervals': record}


def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    brute_controls()
    zero = model(2, 4)
    for p in (0, 2, 3):
        columns = actual(zero, [0] * len(zero[2]), p)
        rows = parallel_query(columns, (1, -3, 2), p)
        assert len(rows) == 1 and not rows[0]['multiplicity']
        COUNTS['actual_zero_map_query_controls'] += 1
    saved = json.loads((ROOT / 'results/research/generic_hook_checks.json').read_text())
    records = [hook_control(saved['families'][0], (1, -7, 2), 1),
               hook_control(saved['families'][1], (1, -6, 4), 7)]
    target = next(r for r in records[1]['pair_recurrences'] if r['lower'] == (0, 7, 3))
    assert target['upper'] == (1, 1, 7)
    assert target['components'] == [['5/6', '6/7'], ['8/7', '6/5'], ['9/7', '4/3']]
    assert target['pieces'] == 4
    files = {'docs/hook_wall_recurrence.md', 'scripts/run_hook_wall_checks.py',
             'docs/hook_bridge_reduction.md', 'docs/generic_hook_fans.md',
             'docs/two_colour_threshold_braids.md', 'results/research/generic_hook_checks.json'}
    pending = ['scripts/run_hook_wall_checks.py']
    while pending:
        f = pending.pop()
        for node in ast.walk(ast.parse((ROOT / f).read_text())):
            names = ([node.module] if isinstance(node, ast.ImportFrom) else
                     [v.name for v in node.names] if isinstance(node, ast.Import) else [])
            for name in names:
                target = 'scripts/' + str(name) + '.py'
                if str(name).startswith('run_') and (ROOT / target).exists() and target not in files:
                    files.add(target)
                    pending.append(target)
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_hook_wall_checks.py',
              'scope': 'Exact positive parallel-wall query and a recurring literal pure cross-arm exchange in actual hook quotients. No growing multiplicity theorem, improved asymptotic count, novelty claim or goal completion.',
              'counts': dict(sorted(COUNTS.items())),
              'original_quotient_replay_counts': dict(sorted(BRIDGE_COUNTS.items())), 'hook_queries': records,
              'limits': ['A literal pair can have disconnected eligibility; counting pairs does not count wall pieces.',
                         'The query traverses all shadow intervals, including those without a q-wall; its running time is not bounded just by the number of parallel pieces.',
                         'All-characteristic generic assertions are restricted to these finite selected coherent minors; arbitrary support-matching ranks are not certified.',
                         'The growing-thickness total of remaining cross-arm pieces is still unbounded at the target scale.'],
              'source_sha256': {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT / 'results/research/hook_wall_checks.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'counts': report['counts'], 'source_hashes': len(files)}))


if __name__ == '__main__':
    main()
