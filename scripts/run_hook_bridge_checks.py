#!/usr/bin/env python3
"""Original-field rank and full cone-patch checks for the hook bridge reduction."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import gcd, prod
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_generic_hook_checks import model as saved_model, exchange_normals
from run_signed_face_checks import field_basis, field_reduce
from run_square_zero_extension_checks import extension, matrix_controls
from run_structural_complexity_checks import clip, area2, dot

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()
TRIANGLE = [tuple(Q(i == j) for i in range(3)) for j in range(3)]


def model(t, a):
    assert t >= 2 and a >= 2 * t
    D = tuple((i, j) for i, j in product(range(a), repeat=2) if i < t or j < t)
    index = {v: i for i, v in enumerate(D)}
    G = ((a, 0), (t, t), (0, a))
    allowed = (tuple(v for v in D if v[0] >= t or v[1] >= a - t),
               tuple(v for v in D if v[0] >= t or v[1] >= t),
               tuple(v for v in D if v[0] >= a - t or v[1] >= t))
    parameters = tuple((g, *v) for g, group in enumerate(allowed) for v in group)
    pi = {s: i for i, s in enumerate(parameters)}
    shifts = (tuple(product(range(a - t), range(t))),
              tuple((i, j) for i, j in product(range(a - t), repeat=2) if i < t or j < t),
              tuple(product(range(t), range(a - t))))
    columns = {(1, *v): {index[v]: None} for v in D}
    for g in range(3):
        for shift in shifts[g]:
            label = (0, G[g][0] + shift[0], G[g][1] + shift[1])
            assert label not in columns
            columns[label] = {index[w]: pi[(g, *v)] for v in allowed[g]
                             if (w := tuple(x + y for x, y in zip(v, shift))) in index}
            assert columns[label]
    Y = {index[v] for v in D if v[0] >= t}
    Z = {index[v] for v in D if v[1] >= t}
    ey = tuple(v for v, c in columns.items() if set(c) <= Y)
    ez = tuple(v for v, c in columns.items() if set(c) <= Z)
    H = tuple(v for v, c in columns.items() if set(c) & Y and set(c) & Z)
    box = frozenset(v for v in columns if v[0] and v[1] < t and v[2] < t)
    fixed = frozenset((0, *v) for v in D) | box
    l = t * (a - t)
    expected = {(0, a + i, j) for i, j in product(range(t), repeat=2)}
    expected |= {(0, t + i, t + j) for i, j in product(range(t), repeat=2)}
    expected |= {(0, i, a + j) for i, j in product(range(t), repeat=2)}
    assert set(H) == expected and len(H) == 3 * t * t
    assert len(ey) == len(ez) == 3 * l - 2 * t * t
    assert set(ey) | set(ez) | set(H) | box == set(columns)
    assert not Y & Z and len(Y) == len(Z) == l
    near = set()
    for arm in (ey, ez):
        for u in arm:
            if not u[0]:
                continue
            for v in (*arm, *H):
                if not v[0]:
                    near.add(tuple(x - y for x, y in zip(u, v)))
    assert len(near) <= 256 * a * t
    assert all(prod(x + 1 for x in v) <= 8 * a * t for v in columns)
    if a == t * t + 1:
        old = saved_model(t)
        assert D == old[0] and parameters == old[4] and columns == old[5]
        COUNTS['independent_original_carrier_model_equalities'] += 1
    COUNTS['exact_potential_bridge_partitions'] += 1
    return D, G, parameters, columns, (ey, ez, H), box, fixed, near


def actual(data, coefficients, p, replay=False):
    D, G, parameters, columns, groups, box, fixed, near = data
    def clean(v):
        return {i: c % p if p else Q(c) for i, c in v.items() if (c % p if p else c)}
    # Keep zero columns; specializing a possible bridge to a loop is allowed.
    values = {v: clean({i: 1 if s is None else coefficients[s] for i, s in c.items()})
              for v, c in columns.items()}
    residual = {v: values[v] for group in groups for v in group}
    if replay:
        tails = [{(i, j): coefficients[s] for s, (g0, i, j) in enumerate(parameters)
                  if g0 == g and coefficients[s]} for g in range(3)]
        n, ev, ground = extension(D, D, G, tails, p)
        assert n == 2 * len(D)
        for v, c in values.items():
            assert ev(v) == {len(D) + i: x for i, x in c.items()}
            COUNTS['independent_original_quotient_column_equalities'] += 1
        reference = tuple((0, *v) for v in D) + tuple((1, *v) for v in D)
        matrix_controls(n, ev, reference, ground, p, COUNTS)
    return residual


def normals(B, columns, p, maximize=False):
    _, pivots = field_basis(B, columns.__getitem__, p)
    assert len(pivots) == len(B)
    out = set()
    for a, v in columns.items():
        if a in B:
            continue
        rem, co = field_reduce(v, pivots, p)
        assert not rem
        for b, coefficient in co.items():
            if not coefficient:
                continue
            q = tuple(x - y for x, y in zip(b if maximize else a, a if maximize else b))
            d = gcd(*q)
            assert d
            out.add(tuple(x // d for x in q))
    return out


def polygon(N):
    poly = list(TRIANGLE)
    for q in N:
        poly = clip(poly, q, 1)
    assert poly and area2(poly) > 0
    return poly


def compress(data, columns, w, p):
    groups = [tuple(sorted(group, key=lambda v: (dot(w, v), v))) for group in data[4]]
    assert len({dot(w, v) for v in columns}) == len(columns)
    by, py = field_basis(groups[0], columns.__getitem__, p)
    bz, pz = field_basis(groups[1], columns.__getitem__, p)
    k = len(groups[2])
    U = (*by, *bz, *groups[2])
    reference, pivots = field_basis((*by, *bz), columns.__getitem__, p)
    assert len(reference) == len(by) + len(bz)
    co = []
    for v in groups[2]:
        rem, coeff = field_reduce(columns[v], pivots, p)
        assert not rem
        co.append(coeff)
    dual = {v: {j: (-co[j].get(v, 0)) % p if p else -co[j].get(v, 0)
                for j in range(k) if co[j].get(v, 0)} for v in reference}
    dual.update({v: {j: 1} for j, v in enumerate(groups[2])})
    for j in range(k):
        row = Counter()
        for v in U:
            for i, c in columns[v].items():
                row[i] += dual[v].get(j, 0) * c
        assert all(not (x % p if p else x) for x in row.values())
        COUNTS['independent_dual_kernel_row_identities'] += 1
    small = []
    for group in (by, bz, groups[2]):
        keep, _ = field_basis(tuple(reversed(group)), dual.__getitem__, p)
        assert len(keep) <= k
        small.extend(keep)
    assert len(small) <= 3 * k and len(set(small)) == len(small)
    B, _ = field_basis(sorted(columns, key=lambda v: dot(w, v)), columns.__getitem__, p)
    restricted, _ = field_basis(sorted(U, key=lambda v: dot(w, v)), columns.__getitem__, p)
    best, _ = field_basis(sorted(U, key=lambda v: -dot(w, v)), dual.__getitem__, p)
    tiny, _ = field_basis(sorted(small, key=lambda v: -dot(w, v)), dual.__getitem__, p)
    assert len(best) == len(tiny) == k
    assert set(B) == set(restricted) == set(U) - set(best)
    assert set(best) == set(tiny)
    original_N = normals(B, columns, p)
    tiny_N = normals(tiny, {v: dual[v] for v in small}, p, True)
    order_N = {tuple(y - x for x, y in zip(v, u)) for group in groups
               for v, u in zip(group, group[1:])}
    assert all(dot(w, q) > 0 for q in original_N | tiny_N | order_N)
    A = polygon(original_N | order_N)
    T = polygon(tiny_N | order_N)
    assert all(dot(v, q) >= 0 for v in A for q in tiny_N)
    assert all(dot(v, q) >= 0 for v in T for q in original_N)
    assert area2(A) == area2(T)
    COUNTS['independent_entire_cone_patch_equalities'] += 1
    COUNTS['actual_dual_compressions'] += 1
    return frozenset(B), {'field_characteristic': p, 'residual_rank': len(reference),
                           'decision_rank': k, 'retained_labels': len(small),
                           'patch_doubled_area': str(area2(A))}


def strict_weight(poly, labels):
    for seed in range(1, 1000):
        w = tuple(sum(Q(seed + 1) ** i * v[j] for i, v in enumerate(poly)) for j in range(3))
        if min(w) > 0 and len({dot(w, a) for a in labels}) == len(labels):
            return w
    raise AssertionError('No strict weight in a full-dimensional polygon.')


def generic_control():
    saved = json.loads((ROOT / 'results/research/generic_hook_checks.json').read_text())['families'][0]
    data = model(2, 5)
    p = 1000003
    columns = actual(data, saved['independent_reference_coefficients'], p, True)
    old = saved_model(2)
    assert all(len(field_basis(tuple(v for v in map(tuple, B) if v not in data[5]), columns.__getitem__, p)[0]) == 12
               for B in saved['support_fan_bases'])
    COUNTS['common_specialization_whole_generic_basis_survivals'] += 89
    pieces = set()
    records = []
    total_area = Q(0)
    for raw in saved['support_fan_bases']:
        B = frozenset(map(tuple, raw))
        N = exchange_normals(B, old[5], len(old[0]))
        poly = polygon(N)
        total_area += area2(poly)
        w = strict_weight(poly, columns)
        selected, record = compress(data, columns, w, p)
        assert selected == B - data[5]
        records.append(record)
        for u, v in zip(poly, poly[1:] + poly[:1]):
            if any(u[i] == v[i] == 0 for i in range(3)):
                continue
            q = next(q for q in N if dot(q, u) == dot(q, v) == 0)
            if q[0]:
                q = q if q[0] > 0 else tuple(-x for x in q)
                pieces.add((q, tuple(sorted((u, v)))))
    assert total_area == 1
    COUNTS['complete_generic_original_fan_area_equalities'] += 1
    remaining = {(q, edge) for q, edge in pieces if q not in data[7]}
    assert len(pieces) == 124 and len(remaining) == 2
    assert {q for q, _ in remaining} == {(1, -7, 2), (1, 2, -7)}
    COUNTS['complete_finite_generic_fan_replays'] += 1
    return {'thickness': 2, 'a': 5, 'length': 32, 'field_characteristic': p,
            'common_coefficients': saved['independent_reference_coefficients'],
            'initials': 89, 'known_direction_inventory': len(data[7]),
            'cross_layer_pieces': len(pieces), 'remaining_cross_arm_pieces': len(remaining),
            'remaining_directions': sorted({q for q, _ in remaining}), 'patches': records}


def special_controls():
    rng = Random(20261022)
    records = []
    for t, a, p, kind in [(2, 4, 0, 'zero'), (2, 8, 2, 'random'),
                           (3, 10, 3, 'random'), (4, 17, 2, 'sparse'),
                           (4, 17, 3, 'sparse')]:
        data = model(t, a)
        coefficients = [0] * len(data[2])
        if kind == 'random':
            coefficients = [rng.randrange(p) for _ in coefficients]
        if kind == 'sparse':
            nonzero = {(0, 0, 15), (0, 1, 13), (1, 1, 12), (1, 2, 13), (1, 14, 1)}
            coefficients = [int(label in nonzero) for label in data[2]]
        columns = actual(data, coefficients, p, True)
        if kind == 'sparse':
            gap = json.loads((ROOT / 'results/research/coherent_convolution_checks.json').read_text())
            bad = frozenset(map(tuple, gap['support_cone']['basis'])) - data[6]
            assert len(bad) == 104
            assert len(field_basis(tuple(bad), columns.__getitem__, p)[0]) == 103
            COUNTS['independent_retained_sparse_rank_gap_controls'] += 1
        for trial in range(3):
            w = ((39, 43, 109) if kind == 'sparse' and trial == 0 else
                 tuple(rng.randrange(1, 1000000) for _ in range(3)))
            while len({dot(w, v) for v in columns}) != len(columns):
                w = tuple(rng.randrange(1, 1000000) for _ in range(3))
            B, record = compress(data, columns, w, p)
            if kind == 'zero':
                assert all(v[0] == 1 for v in B)
            records.append({'t': t, 'a': a, 'kind': kind, **record})
    return records


def main():
    inventory = []
    for t, a in [(2, 4), (2, 5), (2, 12), (3, 10), (4, 17), (6, 37)]:
        data = model(t, a)
        inventory.append({'t': t, 'a': a, 'length': 2 * len(data[0]),
                          'pure_arm_columns': list(map(len, data[4][:2])),
                          'bridges': len(data[4][2]), 'known_directions': len(data[7])})
    generic = generic_control()
    special = special_controls()
    files = {'docs/hook_bridge_reduction.md', 'scripts/run_hook_bridge_checks.py',
             'docs/generic_hook_fans.md', 'docs/square_zero_extension_fans.md',
             'docs/layer_shadow_budget.md', 'docs/cross_layer_direction_budget.md',
             'docs/univariate_cross_fans.md', 'docs/planar_overlap_endomorphisms.md',
             'results/research/generic_hook_checks.json', 'results/research/coherent_convolution_checks.json'}
    pending = ['scripts/run_hook_bridge_checks.py']
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
              'command': 'python3 -B scripts/run_hook_bridge_checks.py',
              'scope': 'Exact conditional small-dual and global direction localization for every compatible flat hook specialization. No new full-family exponent, support-matching rank theorem, current novelty claim or goal completion.',
              'counts': dict(sorted(COUNTS.items())), 'potential_inventories': inventory,
              'complete_generic_fixture': generic, 'specialization_controls': special,
              'proof_limits': ['The rank-3t^2 matrix is conditional on internal orders and changes between regions.',
                               'The central-generator image and actual ranks are retained; sparse symbolic rank failure is not relaxed to matching.',
                               'Known directions contribute O(n^2 log n) cells; the remaining actual cross-arm piece count is unproved.',
                               'Small single-convolution minor experiments are private exploratory data, not a generic rank lemma.'],
              'source_sha256': {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT / 'results/research/hook_bridge_checks.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'counts': report['counts'], 'source_hashes': len(files)}))


if __name__ == '__main__':
    main()
