#!/usr/bin/env python3
"""Actual-field certificates for central prefix regularity and outer-arm exclusion."""
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
from run_hook_bridge_checks import COUNTS as REPLAY_COUNTS, actual, model
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import dot

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def transpose(a):
    return a[0], a[2], a[1]


def rectangles(t, a, r, arm):
    assert t <= r <= a - t
    upper = tuple((1, i, j) for i, j in product(range(t, t + r), range(t)))
    plain = tuple((0, i, j) for i, j in product(range(2 * t, a + t - r), range(t, 2 * t)))
    rows = tuple((i, j) for i, j in product(range(t + r, a), range(t)))
    if arm == 'Z':
        upper, plain = tuple(map(transpose, upper)), tuple(map(transpose, plain))
        rows = tuple((j, i) for i, j in rows)
    assert len(upper) + len(plain) == t * (a - t) and len(plain) == len(rows)
    return upper, plain, rows


def specialize_controls():
    cases = []
    for t, a in ((2, 4), (2, 5), (2, 20), (2, 32), (3, 10), (3, 25), (4, 17), (4, 33)):
        data = model(t, a)
        index = {v: i for i, v in enumerate(data[0])}
        for r, arm in product(range(t, a - t + 1), ('Y', 'Z')):
            upper, plain, tail = rectangles(t, a, r, arm)
            parameter = data[2].index((1, r, 0) if arm == 'Y' else (1, 0, r))
            for label in upper:
                assert data[3][label] == {index[label[1:]]: None}
            for label, row in zip(plain, tail):
                # Parameter extraction and literal multiplication are independent
                # ways to evaluate this monomial specialization.
                image = {j: 1 for j, s in data[3][label].items() if s == parameter}
                shift = label[1] - t, label[2] - t
                literal = (r + shift[0], shift[1]) if arm == 'Y' else (shift[0], r + shift[1])
                assert literal == row and image == {index[row]: 1}
                COUNTS['independent_literal_monomial_specializations'] += 1
            assert len(set(v[1:] for v in upper) | set(tail)) == t * (a - t)
            COUNTS['primitive_prefix_minor_identity_specializations'] += 1
        L = a - 2 * t
        assert sum(2 * t * (a - t - r) for r in range(t, a - t + 1)) == t * L * (L + 1)
        cases.append({'t': t, 'a': a, 'length': 2 * len(data[0]), 'central_product_degree': t * L * (L + 1)})
    return cases


def binomial_value(L, d):
    s, value = L - d, Q(1)
    for j in range(s):
        value *= Q(factorial(L + j) * factorial(j), factorial(d + j) * factorial(L - d + j))
    assert value.denominator == 1
    return int(value)


def explicit_coefficients(data, t, a, rng, p):
    coefficients = [rng.randrange(p) for _ in data[2]]
    L = a - 2 * t
    for i, (g, y, z) in enumerate(data[2]):
        if g == 1:
            coefficients[i] = (comb(L, y - t) if z == 0 and t <= y <= a - t else
                               comb(L, z - t) if y == 0 and t <= z <= a - t else 0) % p
    return coefficients


def prefix_ranks(data, coefficients, t, a, p):
    columns = actual(data, coefficients, p, replay=True)
    index = {v: i for i, v in enumerate(data[0])}
    L = a - 2 * t
    for r in range(t, a - t + 1):
        d, size = r - t, a - t - r
        T = [[comb(L, k) if 0 <= (k := d + i - j) <= L else 0 for j in range(size)]
             for i in range(size)]
        exact = bareiss(T)
        assert exact == binomial_value(L, d) and exact % p
        COUNTS['independent_binomial_determinants_vs_product'] += 1
        for arm in ('Y', 'Z'):
            upper, plain, rows = rectangles(t, a, r, arm)
            row_indices = {index[row] for row in rows}
            restricted = {v: {i: c for i, c in columns[v].items() if i in row_indices} for v in plain}
            B, _ = field_basis(plain, restricted.__getitem__, p)
            assert len(B) == len(plain)
            full, _ = field_basis((*upper, *plain), columns.__getitem__, p)
            assert len(full) == t * (a - t)
            if len(plain) <= 32:
                matrix = [[columns[v].get(index[row], 0) for v in plain] for row in rows]
                assert bareiss(matrix) % p == pow(exact, t, p)
                COUNTS['independent_full_prefix_determinants_vs_tensor_formula'] += 1
            COUNTS['actual_prefix_minor_full_ranks'] += 1
            COUNTS['actual_entire_arm_basis_ranks'] += 1
    return columns


def exclude_case(data, columns, t, a, Y, Z, R, S, u, arm, p):
    b, c = Y - R, S - Z
    assert a + t <= Y < 2 * a - t and 0 <= Z < t and 0 <= R < t and t <= S < a
    assert c > 6 * t and b > 0 and u > Q(c, b)
    e, f = (0, Y, Z), (1, R, S)
    w = (b * u - c, u, Q(1))
    if arm == 'Z':
        e, f, w = transpose(e), transpose(f), transpose(w)
    tau = dot(w, e)
    assert tau == dot(w, f) and min(w) > 0
    if Y * u + Z > (2 * t - 1) * u + a - 1:
        chosen = tuple((0, i, j) for i, j in product(range(t, 2 * t), range(2 * t, a)))
        chosen = tuple(map(transpose, chosen)) if arm == 'Z' else chosen
        target = f
        branch = 'opposite_tail_cheap_span'
    else:
        assert Q(3 * t, a) < u < 2 and u > Q(2 * t, a - t)
        r = int(Q(2 * t) // u)
        assert t <= r < a - t
        U, P, _ = rectangles(t, a, r, arm)
        chosen, target = (*U, *P), e
        assert (r + 1) * u + Z - 2 * t + 1 > 0
        assert S - t * u - r * u - t + 1 > 0
        branch = 'own_arm_prefix_cheap_span'
    assert all(dot(w, v) < tau for v in chosen)
    B, pivots = field_basis(chosen, columns.__getitem__, p)
    remainder, _ = field_reduce(columns[target], pivots, p)
    assert not remainder
    assert len(B) == (t * (a - 2 * t) if branch == 'opposite_tail_cheap_span' else t * (a - t))
    COUNTS[branch] += 1
    COUNTS['independent_original_strict_cheap_span_exclusions'] += 1
    return {'arm': arm, 'b': b, 'c': c, 'u': str(u), 'branch': branch,
            'chosen_rank': len(B), 'strict_cost_gap': str(tau - max(dot(w, v) for v in chosen))}


def field_controls():
    rng = Random(2026102219)
    records = []
    for t, a in ((2, 20), (2, 32), (3, 25), (4, 33)):
        p = 101
        assert p > 2 * (a - 2 * t)
        data = model(t, a)
        coefficients = explicit_coefficients(data, t, a, rng, p)
        columns = prefix_ranks(data, coefficients, t, a, p)
        exclusions = []
        for trial in range(12):
            Y, Z, R = rng.randrange(a + t, 2 * a - t), rng.randrange(t), rng.randrange(t)
            S = rng.randrange(6 * t + Z + 1, a) if 6 * t + Z + 1 < a else a - 1
            if S - Z <= 6 * t:
                Z, S = 0, a - 1
            b, c = Y - R, S - Z
            # These two samples exercise the direct opposite-tail branch and
            # the nontrivial prefix branch close to the positivity boundary.
            for u in (Q(c, b) + Q(1, 1000), Q(3)):
                for arm in ('Y', 'Z'):
                    exclusions.append(exclude_case(data, columns, t, a, Y, Z, R, S, u, arm, p))
        records.append({'t': t, 'a': a, 'length': 2 * len(data[0]), 'field_characteristic': p,
                        'central_image': 'y^t*(1+y)^(a-2t)+z^t*(1+z)^(a-2t)',
                        'common_coefficients': coefficients, 'prefix_minor_pairs': a - 2 * t + 1,
                        'exclusions': exclusions})
    return records



def eligibility(columns, e, f, w, p):
    tau = dot(w, e)
    assert tau == dot(w, f)
    assert {v for v in columns if dot(w, v) == tau} == {e, f}
    cheap = tuple(v for v in columns if dot(w, v) < tau)
    _, pivots = field_basis(cheap, columns.__getitem__, p)
    E, _ = field_reduce(columns[e], pivots, p)
    F, _ = field_reduce(columns[f], pivots, p)
    if not E or not F:
        return False, cheap
    i = next(iter(E))
    if i not in F:
        return False, cheap
    ratio = E[i] * pow(F[i], -1, p) % p
    return all((E.get(j, 0) - ratio * F.get(j, 0)) % p == 0 for j in set(E) | set(F)), cheap


def own_outer_invisibility_controls():
    t, a, p = 2, 32, 101
    data = model(t, a)
    rng = Random(2026102223)
    co = explicit_coefficients(data, t, a, rng, p)
    changed = [rng.randrange(p) if v[0] == 0 else c for v, c in zip(data[2], co)]
    first, second = actual(data, co, p), actual(data, changed, p)
    records = []
    for trial in range(12):
        Y, Z, R = rng.randrange(2 * t, a), rng.randrange(t, 2 * t), rng.randrange(t)
        S = rng.randrange(Z + 8 * t + 1, a)
        e, f = (0, Y, Z), (1, R, S)
        b, c = Y - R, S - Z
        for u in (Q(c, b) + Q(1, 10007), Q(2) + Q(1, 10007)):
            w = (b * u - c, u, Q(1))
            tau = dot(w, e)
            before, cheap = eligibility(first, e, f, w, p)
            after, cheap2 = eligibility(second, e, f, w, p)
            assert cheap == cheap2 and before == after
            if tau >= a * u:
                assert not before
                r = t + int(Q(2 * t) // u)
                # The opposite-tail branch may already prove inactivity.
                if tau > (2 * t - 1) * u + a - 1:
                    COUNTS['long_central_inactive_by_opposite_tail'] += 1
                else:
                    assert u < 2 and t <= r < a - t
                    U, P, _ = rectangles(t, a, r, 'Y')
                    assert all(dot(w, v) < tau for v in (*U, *P))
                    COUNTS['long_central_inactive_by_prefix_basis'] += 1
            else:
                assert all(v[1] < a or v[0] for v in cheap)
                assert all(first[v] == second[v] for v in cheap)
                COUNTS['long_central_cheap_columns_ignore_own_outer_image'] += 1
            COUNTS['actual_central_pair_eligibility_outer_image_invariance'] += 1
            records.append({'e': e, 'f': f, 'u': str(u), 'eligible': before,
                            'common_cost_below_own_outer_head': tau < a * u})
    return {'t': t, 'a': a, 'field_characteristic': p,
            'fixed_central_image': 'y^t*(1+y)^(a-2t)+z^t*(1+z)^(a-2t)',
            'first_coefficients': co, 'changed_coefficients': changed, 'cases': records,
            'eligible_sampled_pairs': sum(record['eligible'] for record in records),
            'scope': 'Finite invariance controls; no uniform central-wall survival or absence claim.'}


def excluded_strata_controls():
    records = []
    for t, a, p in ((2, 16, 2), (3, 20, 3), (4, 33, 101)):
        data = model(t, a)
        co = [int(v in ((1, t, 0), (1, 0, t))) for v in data[2]]
        C = actual(data, co, p)
        index = {v: i for i, v in enumerate(data[0])}
        ranks = []
        for r, arm in product((t, t + 1), ('Y', 'Z')):
            _, plain, rows = rectangles(t, a, r, arm)
            idx = {index[v] for v in rows}
            ev = lambda v: {i: c for i, c in C[v].items() if i in idx}
            rank = len(field_basis(plain, ev, p)[0])
            assert (rank == len(plain)) == (r == t)
            ranks.append({'r': r, 'arm': arm, 'rank': rank, 'minor_order': len(plain)})
            COUNTS['unit_coefficients_do_not_imply_prefix_regularity'] += 1
        records.append({'t': t, 'a': a, 'field_characteristic': p, 'ranks': ranks})
    return records


def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    specializations = specialize_controls()
    fields = field_controls()
    excluded = excluded_strata_controls()
    invisibility = own_outer_invisibility_controls()
    assert COUNTS['own_arm_prefix_cheap_span'] and COUNTS['opposite_tail_cheap_span']
    files = {'docs/hook_central_prefix.md', 'scripts/run_hook_central_prefix_checks.py',
             'docs/hook_bridge_reduction.md', 'docs/generic_hook_fans.md',
             'docs/layer_shadow_budget.md', 'docs/cross_layer_direction_budget.md',
             'docs/hook_wall_recurrence.md'}
    pending = ['scripts/run_hook_central_prefix_checks.py']
    while pending:
        f = pending.pop()
        for node in ast.walk(ast.parse((ROOT / f).read_text())):
            names = ([node.module] if isinstance(node, ast.ImportFrom) else
                     [x.name for x in node.names] if isinstance(node, ast.Import) else [])
            for name in names:
                target = 'scripts/' + str(name) + '.py'
                if str(name).startswith('run_') and (ROOT / target).exists() and target not in files:
                    files.add(target)
                    pending.append(target)
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_hook_central_prefix_checks.py',
              'scope': 'Uniform outer-source localization under explicit central prefix regularity, on a nonempty open set in every characteristic. The remaining central pieces and the excluded strata are not bounded at the goal scale.',
              'counts': dict(sorted(COUNTS.items())),
              'original_quotient_replay_counts': dict(sorted(REPLAY_COUNTS.items())),
              'identity_specializations': specializations, 'actual_field_controls': fields,
              'excluded_stratum_controls': excluded, 'own_outer_invisibility_controls': invisibility,
              'limits': ['No arbitrary support-matching rank theorem is assumed.',
                         'The central prefix condition is not implied just by its two leading coefficients.',
                         'The transverse bound6t is not asserted optimal.',
                         'The central-source cross-arm count, arbitrary coefficient strata, unrestricted sharp endpoint and current priority remain open.'],
              'source_sha256': {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT / 'results/research/hook_central_prefix_checks.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'counts': report['counts'], 'source_hashes': len(files)}))


if __name__ == '__main__':
    main()
