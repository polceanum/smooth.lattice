#!/usr/bin/env python3
"""Original-matrix sector recursion and coherent axis-data contrast checks."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import comb
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_generic_hook_checks import bareiss
from run_hook_bridge_checks import normals, polygon
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_signed_face_checks import field_basis
from run_square_zero_extension_checks import extension, matrix_controls
from run_structural_complexity_checks import area2, dot

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def clean(v, p):
    return {i: c % p if p else Q(c) for i, c in v.items() if (c % p if p else c)}


def hook(t, a, f, p, replay=True):
    assert t >= 1 and a >= 2 * t
    D = tuple(v for v in product(range(a), repeat=2) if v[0] < t or v[1] < t)
    assert all(v in D and (v[0] >= t or v[1] >= t) for v in f)
    n, ev, ground = extension(D, D, ((a, 0), (t, t), (0, a)), ({}, f, {}), p)
    assert n == 4 * a * t - 2 * t * t
    reference = tuple((x, *v) for x in range(2) for v in D)
    if replay:
        matrix_controls(n, ev, reference, ground, p, COUNTS)
    return {'t': t, 'a': a, 'n': n, 'D': D, 'ev': ev, 'ground': ground,
            'reference': reference, 'columns': {v: ev(v) for v in ground}, 'p': p}


def lift(v):
    return v[0], v[1] + 1, v[2] + 1


def parent_models(t, a, f, p, alpha=1, beta=1):
    child = hook(t, a, f, p)
    shifted = {(i + 1, j + 1): c for (i, j), c in f.items()}
    full = dict(shifted)
    full[t + 1, 0], full[0, t + 1] = alpha, beta
    parent, dropped = hook(t + 1, a + 1, full, p), hook(t + 1, a + 1, shifted, p)
    assert parent['n'] == dropped['n'] == child['n'] + 4 * a + 2
    halo = frozenset((x, i, j) for x in range(2) for i, j in product(range(a + 1), repeat=2)
                     if not i or not j)
    assert len(halo) == 4 * a + 2
    pindex = {v: i for i, v in enumerate(parent['reference'])}
    rowmap = {i: pindex[lift(v)] for i, v in enumerate(child['reference'])}
    s = t + 1
    for v in parent['ground']:
        difference = Counter(parent['ev'](v))
        difference.subtract(dropped['ev'](v))
        expected = Counter()
        helpers = []
        if not v[0] and v[1] >= s and v[2] >= s:
            helpers = [(alpha, (1, v[1], v[2] - s)), (beta, (1, v[1] - s, v[2]))]
            for coefficient, helper in helpers:
                for i, c in parent['ev'](helper).items():
                    expected[i] += coefficient * c
        assert clean(difference, p) == clean(expected, p)
        COUNTS['independent_triangular_axis_column_identities'] += 1
        if v[1] and v[2]:
            before = (v[0], v[1] - 1, v[2] - 1)
            assert dropped['ev'](v) == {rowmap[i]: c for i, c in child['ev'](before).items()}
        elif v in halo:
            assert dropped['ev'](v) == {pindex[v]: 1}
        else:
            assert not dropped['ev'](v)
        COUNTS['independent_padding_row_map_identities'] += 1
    return child, parent, dropped, halo


def greedy(data, w):
    B, _ = field_basis(sorted(data['ground'], key=lambda v: (dot(w, v), v)), data['columns'].__getitem__, data['p'])
    assert len(B) == data['n']
    return frozenset(B)


def patch_controls(t, a, f, p, label, alpha=1, beta=1):
    child, parent, dropped, halo = parent_models(t, a, f, p, alpha, beta)
    s = t + 1
    sector = {(-1, s, 0), (-1, 0, s)}
    target_area = area2(polygon(sector))
    base = parent['n'] + 1
    queue = [(1, base, base * base)]
    visited, mapped = {}, set()
    while queue:
        w = queue.pop()
        assert all(dot(q, w) > 0 for q in sector)
        B = greedy(child, w)
        if B in visited:
            continue
        NC = normals(B, child['columns'], p)
        assert all(dot(q, w) > 0 for q in NC)
        BP = halo | frozenset(map(lift, B))
        assert len(BP) == parent['n'] and greedy(dropped, w) == BP
        NP = normals(BP, parent['columns'], p)
        assert all(dot(q, w) > 0 for q in NP)
        pc, pp = polygon(NC | sector), polygon(NP | sector)
        assert area2(pc) == area2(pp)
        assert all(dot(q, v) >= 0 for v in pc for q in NP | sector)
        assert all(dot(q, v) >= 0 for v in pp for q in NC | sector)
        visited[B] = pc; mapped.add(BP)
        COUNTS['independent_entire_original_sector_cone_equalities'] += 1
        for left, right in zip(pc, pc[1:] + pc[:1]):
            if any(left[i] == right[i] == 0 for i in range(3)):
                continue
            if any(dot(q, left) == dot(q, right) == 0 for q in sector):
                continue
            q = next(q for q in NC if dot(q, left) == dot(q, right) == 0)
            new, _ = exact_adjacent_weight(parent['n'], left, right, q)
            assert all(dot(v, new) > 0 for v in sector)
            queue.append(new)
    assert len(mapped) == len(visited) and sum(map(area2, visited.values())) == target_area
    COUNTS['independent_complete_sector_fan_coverages'] += 1
    return {'label': label, 't': t, 'a': a, 'field_characteristic': p,
            'child_coefficients': [[i, j, c] for (i, j), c in sorted(f.items())],
            'axis_coefficients': [alpha, beta], 'child_length': child['n'], 'parent_length': parent['n'],
            'matched_sector_cones': len(visited), 'normalized_sector_area': str(target_area),
            'parent_maximum_gaps': [max(0, a - 2 * t - 2)] * 2,
            'scope': 'Complete finite fan-patch equality inside the open low-price sector, not equality of the whole fans.'}


def lower_family(a, p):
    A = a - 1
    f = {(1 + i, 0): comb(A - 1, i) for i in range(A)}
    f[0, 1], f[0, A] = 1, 1
    if p:
        assert p > 2 * A
        f = {v: c % p for v, c in f.items() if c % p}
    child, parent, dropped, halo = parent_models(1, a, f, p)
    # All possible source/identity prefixes occurring in the pure Y merge.
    cols = child['columns']
    for r in range(A + 1):
        for k in range(A):
            if k and r < k + 1:
                continue
            upper = tuple((1, 1 + j, 0) for j in range(r))
            source = tuple((0, 1 + i, 1) for i in range(1, k + 1))
            assert len(field_basis((*upper, *source), cols.__getitem__, p)[0]) == min(A, r + k)
            COUNTS['independent_all_price_prefix_binomial_ranks'] += 1
    states, records = set(), []
    hs = [h for h in range(2, A - 1) if (h - A) % 2 == 0]
    for h in hs:
        u, lam = Q(2, 2 * h + 1), Q(1, 2 * h + 1)
        assert 0 < lam < u < Q(1, 2)
        Bs = []
        for change in (-Q(1, 4), Q(1, 4)):
            w = (lam + change * u, u, Q(1))
            # The formula for the arm optimum is independently certified by
            # the full original-field basis inequalities below.
            choices = [(1, j + 1, 0) for j in range(A)] + [(0, i + 1, 1) for i in range(1, A)]
            BY = sorted(choices, key=lambda v: (dot(w, v), v))[:A]
            fixed = {(0, *v) for v in child['D']} | {(1, 0, 0)}
            B = frozenset(fixed | set(BY) | {(1, 0, j + 1) for j in range(A - 1)} | {(0, 1, 1)})
            assert len(B) == child['n']
            N = normals(B, child['columns'], p)
            assert all(dot(q, w) > 0 for q in N)
            BP = halo | frozenset(map(lift, B))
            NP = normals(BP, parent['columns'], p)
            assert all(dot(q, w) > 0 for q in NP)
            states.add(BP);Bs.append(BP)
            COUNTS['independent_original_adjacent_family_minima'] += 1
        i, j = (A - h) // 2, (A + h) // 2
        e, v = lift((0, i + 1, 1)), lift((1, j + 1, 0))
        assert Bs[0] - Bs[1] == {v} and Bs[1] - Bs[0] == {e}
        assert tuple(x - y for x, y in zip(v, e)) == (1, h, -1)
        COUNTS['independent_coherent_lifted_family_walls'] += 1
        records.append({'h': h, 'wall_weight': list(map(str, (lam, u, Q(1)))),
                        'parent_lower': e, 'parent_upper': v, 'direction': [1, h, -1]})
    assert len(hs) == A // 2 - 1 and len(states) == A // 2
    # The axis determinants are independent of the off-axis coefficients.
    D = parent['D'];index = {v: i for i, v in enumerate(D)}
    for r in range(2, a):
        plain = tuple((0, y, z) for y, z in product(range(4, a + 3 - r), range(2, 4)))
        rows = tuple((y, z) for y, z in product(range(2 + r, a + 1), range(2)))
        assert len(plain) == len(rows)
        if len(plain) <= 40:
            matrix = [[parent['ev'](v).get(len(D) + index[row], 0) for v in plain] for row in rows]
            value = bareiss(matrix)
            assert (value % p if p else value) == int(r in (2, a - 1))
            COUNTS['independent_parent_constant_axis_prefix_determinants'] += 1
    return {'a': a, 'field_characteristic': p, 'child_length': child['n'], 'parent_length': parent['n'],
            'parent_maximum_gaps': [a - 4, a - 4], 'certified_distinct_parent_states': len(states),
            'coherent_walls': records, 'comparison_family_full_cones': 5,
            'axis_only_approximation_factor_squared_lower_bound': str(Q(len(states), 5)),
            'scope': 'Actual finite coherent witnesses for the uniform linear lower bound; no full parent count or superquadratic separation is claimed.'}


def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    rng = Random(2026102257)
    random_f = {(i, j): rng.randrange(101) for i, j in product(range(5), repeat=2)
                if (i < 2 or j < 2) and (i >= 2 or j >= 2)}
    zero_leading = {(i, j): rng.randrange(2) for i, j in product(range(7), repeat=2)
                    if i and j and (i < 2 or j < 2) and (i >= 2 or j >= 2)}
    explicit = {(1 + i, 0): comb(6, i) for i in range(7)}
    explicit[0, 1], explicit[0, 7] = 1, 1
    patches = [patch_controls(2, 5, random_f, 101, 'arbitrary-central-image', 2, 3),
               patch_controls(2, 7, zero_leading, 2, 'zero-leading-central-image'),
               patch_controls(3, 10, {}, 101, 'zero-child-image'),
               patch_controls(1, 8, explicit, 101, 'explicit-binomial-child')]
    families = [lower_family(a, p) for a, p in ((8, 0), (12, 0), (20, 0),
                                               (8, 101), (12, 101), (20, 101), (30, 101))]
    files = {'docs/hook_axis_sector_reduction.md', 'scripts/run_hook_axis_sector_checks.py',
             'docs/canonical_fan_deflation.md', 'docs/hook_prefix_gaps.md',
             'docs/hook_central_skeleton.md', 'docs/hook_bridge_reduction.md'}
    pending = ['scripts/run_hook_axis_sector_checks.py']
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
              'command': 'python3 -B scripts/run_hook_axis_sector_checks.py',
              'scope': 'Uniform exact low-price sector recursion into a maximal-axis-gap parent, and an actual infinite linear lower bound contrasting with a five-cone family with identical axis determinants. The superquadratic target and full outside-sector counts remain open.',
              'counts': dict(sorted(COUNTS.items())), 'complete_sector_controls': patches,
              'coherent_lower_family_controls': families,
              'limits': ['The recursion identifies only the stated open weight sector, not the entire parent fan.',
                         'Equal axis determinants do not assert equal full prefix ranks.',
                         'The proved linear lower bound does not establish a superquadratic separation.',
                         'Exact priority, unrestricted sharp growth and the original final-paper goal remain unresolved.'],
              'source_sha256': {filename: hashlib.sha256((ROOT / filename).read_bytes()).hexdigest() for filename in sorted(files)}}
    (ROOT / 'results/research/hook_axis_sector_checks.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'counts': report['counts'], 'source_hashes': len(files)}))


if __name__ == '__main__':
    main()
