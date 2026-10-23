#!/usr/bin/env python3
"""Exact common-quotient compression for original hook residue matrices."""
from collections import Counter
from fractions import Fraction as Q
from itertools import permutations
from math import gcd
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_hook_bridge_checks import model, actual, normals, polygon, strict_weight, COUNTS as REPLAY
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import dot, area2
from run_oracle_enumeration_checks import exact_adjacent_weight

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def clean(vector, p):
    return {i: (c % p if p else Q(c)) for i, c in vector.items() if (c % p if p else c)}


def rank(labels, columns, p):
    return len(field_basis(tuple(labels), columns.__getitem__, p)[0])


def dual(labels, columns, p):
    J, piv = field_basis(tuple(labels), columns.__getitem__, p)
    L = tuple(v for v in labels if v not in J)
    co = [field_reduce(columns[v], piv, p)[1] for v in L]
    out = {v: clean({i: -row.get(v, 0) for i, row in enumerate(co)}, p) for v in J}
    out.update({v: {i: 1} for i, v in enumerate(L)})
    assert rank(labels, out, p) == len(L)
    for j in range(len(L)):
        summed = Counter()
        for v in labels:
            for i, c in columns[v].items():
                summed[i] += out[v].get(j, 0)*c
        assert not clean(summed, p)
        COUNTS['independent_common_quotient_dual_row_identities'] += 1
    return out


def transfer(columns, baseline, partition, w, p, cache):
    groups = [tuple(sorted(g, key=lambda v: dot(w, v))) for g in partition]
    assert len({dot(w, v) for v in columns}) == len(columns)
    B0, _ = field_basis(sorted(baseline, key=lambda v: dot(w, v)), baseline.__getitem__, p)
    R = len(B0); assert rank(columns, columns, p) == R
    difference = {}
    for v in columns:
        out = Counter(columns[v]); out.subtract(baseline[v]); difference[v] = clean(out, p)
    db, dp = field_basis(tuple(difference), difference.__getitem__, p); r = len(db)
    projected = {v: field_reduce(columns[v], dp, p)[0] for v in columns}
    for v in columns:
        assert projected[v] == field_reduce(baseline[v], dp, p)[0]
        COUNTS['independent_common_projection_column_equalities'] += 1
    assert rank(columns, projected, p) == R-r
    dq = dual(tuple(B0), projected, p); D = []
    for group in groups:
        keep, _ = field_basis(tuple(v for v in reversed(group) if v in B0), dq.__getitem__, p)
        assert len(keep) <= r; D.extend(keep)
    assert len(D) <= 3*r
    C = tuple(v for v in B0 if v not in D)
    BQ, _ = field_basis(sorted(projected, key=lambda v: dot(w, v)), projected.__getitem__, p)
    assert set(BQ) <= set(B0) and set(C) <= set(BQ)
    cp = field_basis(C, columns.__getitem__, p)[1]
    assert len(cp) == len(C)
    residual = {v: field_reduce(columns[v], cp, p)[0] for v in columns if v not in C}
    small = []
    for group in groups:
        keep, _ = field_basis(tuple(v for v in group if v not in C), residual.__getitem__, p)
        assert len(keep) <= len(D); small.extend(keep)
    tiny = {v: residual[v] for v in small}
    assert len(tiny) <= 9*r and rank(tiny, tiny, p) == len(D) <= 3*r
    BM, _ = field_basis(sorted(columns, key=lambda v: dot(w, v)), columns.__getitem__, p)
    BS, _ = field_basis(sorted(tiny, key=lambda v: dot(w, v)), tiny.__getitem__, p)
    assert set(BQ) <= set(BM) and set(BM) == set(C)|set(BS)
    orders = {tuple(b-a for a, b in zip(v, u)) for g in groups for v, u in zip(g, g[1:])}
    O = normals(B0, baseline, p)|orders
    N = normals(BM, columns, p); T = normals(BS, tiny, p)
    original = polygon(N|O); compressed = polygon(T|O)
    assert all(dot(v, q) >= 0 for v in original for q in T)
    assert all(dot(v, q) >= 0 for v in compressed for q in N)
    assert area2(original) == area2(compressed) > 0
    key = (frozenset(B0), tuple(groups))
    signature = (frozenset(C), tuple(sorted((v, tuple(sorted(c.items()))) for v, c in tiny.items())))
    if key in cache:
        assert cache[key] == signature
        COUNTS['independent_whole_region_encoding_stability_checks'] += 1
    else:
        cache[key] = signature
    COUNTS['independent_entire_original_vs_small_minor_patches'] += 1
    return {'departure_rank': r, 'original_rank': R, 'forced_labels': len(C),
            'small_rank': len(D), 'small_labels': len(tiny),
            'patch_area2': str(area2(original))}, polygon(O)


def plane(q):
    d = gcd(*q); assert d
    q = tuple(v//d for v in q)
    return q if next(v for v in q if v) > 0 else tuple(-v for v in q)


def comparison(data, coefficients, arm, C):
    t = data[1][1][0]
    out = []
    for c, (g, y, z) in zip(coefficients, data[2]):
        long, short = (y, z) if arm == 0 else (z, y)
        if long < t:
            out.append(c)
        else:
            out.append(C[short] if g == 1 and long == t else 0)
    return out


def fiber_controls(data, baseline, arm, t, a, w, p):
    removed = set(); chosen = set()
    for q in range(t, a-t):
        rows = {i for i, v in enumerate(data[0])
                if v[0] == t+q and v[1] < t} if arm == 0 else {
                    i for i, v in enumerate(data[0]) if v[1] == t+q and v[0] < t}
        labels = {(1, t+q, j) if arm == 0 else (1, j, t+q) for j in range(t)}
        labels |= {(0, t+q, t+j) if arm == 0 else (0, t+j, t+q) for j in range(t)}
        assert {v for v, c in baseline.items() if rows & set(c)} <= labels
        assert all(set(baseline[v]) <= rows for v in labels)
        B, _ = field_basis(sorted(labels, key=lambda v: dot(w, v)), baseline.__getitem__, p)
        assert len(B) == t
        chosen.update(B); removed.update(labels)
        COUNTS['independent_isolated_transverse_fibers'] += 1
    rest = tuple(v for v, c in baseline.items() if v not in removed and c)
    for v in rest:
        x, y, z = v if arm == 0 else (v[0], v[2], v[1])
        assert (x in (0, 1) and y < 2*t and z < 2*a) or (x == 0 and a <= y < a+t and z < t)
    RB, _ = field_basis(sorted(rest, key=lambda v: dot(w, v)), baseline.__getitem__, p)
    BB, _ = field_basis(sorted(baseline, key=lambda v: dot(w, v)), baseline.__getitem__, p)
    assert chosen|set(RB) == set(BB)
    inventory = {plane(tuple(x-y for x, y in zip(v, u)))
                 for j, v in enumerate(rest) for u in rest[:j]}
    for g in data[4]:
        inventory |= {plane(tuple(x-y for x, y in zip(v, u))) for j, v in enumerate(g) for u in g[:j]}
    inventory |= {plane((1, 0, -j) if arm == 0 else (1, -j, 0)) for j in range(1, 2*t)}
    assert len(inventory) <= 1024*a*t
    COUNTS['independent_comparison_box_and_plane_inventories'] += 1
    return len(inventory)


def weight_controls(labels, rng):
    b = 4*max(max(v) for v in labels)+5
    out = [tuple(map(Q, v)) for v in permutations((1, b, b*b))]
    while len(out) < 14:
        w = tuple(Q(rng.randrange(1, 10000001)) for _ in range(3))
        if len({dot(w, v) for v in labels}) == len(labels):
            out.append(w)
    return out


def complete_fan(columns, n, p):
    queue = [((1, 2*n+1, (2*n+1)**2), None)]; cells = {}
    while queue:
        w, old = queue.pop()
        B = frozenset(field_basis(sorted(columns, key=lambda v: (dot(w, v), v)), columns.__getitem__, p)[0])
        if old is not None:
            assert B != old
        if B in cells:
            continue
        N = normals(B, columns, p); poly = polygon(N)
        assert all(dot(w, q) > 0 for q in N)
        cells[B] = poly
        for left, right in zip(poly, poly[1:]+poly[:1]):
            if any(left[i] == right[i] == 0 for i in range(3)):
                continue
            q = next(q for q in N if dot(q, left) == dot(q, right) == 0)
            nw, _ = exact_adjacent_weight(n, left, right, q); queue.append((nw, B))
    assert sum(map(area2, cells.values())) == 1
    COUNTS['independent_complete_original_fans'] += 1
    return cells


def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    previous = json.loads((ROOT/'results/research/hook_high_price_partner_checks.json').read_text())
    assert previous['status'] == 'passed'
    for f, sha in previous['source_sha256'].items():
        assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest() == sha
    rng = Random(2026102320); cases = []
    specs = [(2, 5, 101, 0, [1, 1], 1, True), (2, 5, 0, 0, [0, 1], -1, True),
             (3, 8, 2, 0, [0, 1, 0], 1, False), (3, 10, 3, 0, [0, 0, 0], 2, False),
             (3, 9, 3, 1, [1, 2, 1], 1, False), (4, 12, 101, 0, [1, 1, 1, 1], 0, False),
             (6, 17, 3, 0, [1, 0, 0, 0, 0, 0], 1, False),
             (3, 6, 2, 1, [0, 1, 0], 1, False)]
    for t, a, p, arm, C, short_start, full in specs:
        data = model(t, a); co = []
        for g, y, z in data[2]:
            long, short = (y, z) if arm == 0 else (z, y)
            if long < t:
                c = rng.randrange(p or 7)
            else:
                c = C[short] if g == 1 and long == t else 0
                if short_start >= 0 and long == a-1 and short >= short_start:
                    c += 1
            co.append(c)
        cb = comparison(data, co, arm, C)
        for outer in (0, 2):
            assert any((c % p if p else c) for c, parameter in zip(co, data[2]) if parameter[0] == outer)
        columns = actual(data, co, p, True); baseline = actual(data, cb, p, True)
        cache = {}; patches = []; weights = weight_controls(columns, rng)
        inventories = []
        for w in weights:
            row, region = transfer(columns, baseline, data[4], w, p, cache); patches.append(row)
            other = strict_weight(region, columns)
            row, _ = transfer(columns, baseline, data[4], other, p, cache); patches.append(row)
        for w in weights[:6]:
            inventories.append(fiber_controls(data, baseline, arm, t, a, w, p))
        fan = complete_fan(columns, 4*a*t-2*t*t, p) if full else None
        if fan:
            for B, poly in fan.items():
                w = strict_weight(poly, columns)
                row, _ = transfer(columns, baseline, data[4], w, p, cache); patches.append(row)
                assert frozenset(field_basis(sorted(columns, key=lambda v: dot(w, v)), columns.__getitem__, p)[0]) == B
            COUNTS['complete_fan_cones_with_exact_transfer_patches'] += len(fan)
        cases.append({'t': t, 'a': a, 'characteristic': p, 'comparison_arm': 'Y' if arm == 0 else 'Z',
                      'transverse_base': C, 'nonzero_outer_images_allowed': True,
                      'actual_bridge_rank': rank(data[4][2], columns, p),
                      'patches': patches, 'plane_inventory_max': max(inventories),
                      'complete_original_fan_cones': len(fan) if fan else None})
        assert patches[0]['departure_rank'] == (t-short_start if short_start >= 0 else 0)
        print('Departure hook case passed', t, a, p, arm, 'rank', patches[0]['departure_rank'],
              'full cones', len(fan) if fan else None, flush=True)
    abstract = []
    labels = tuple((j % 2, 3*j % 17, (j*j+5) % 23) for j in range(12))
    groups = tuple(tuple(v for j, v in enumerate(labels) if j % 3 == i) for i in range(3))
    for p, target in ((0, 1), (2, 1), (3, 2), (101, 2)):
        while True:
            baseline = {v: clean({i: int(i == j) if j < 5 else rng.randrange(p or 7)
                                  for i in range(5)}, p) for j, v in enumerate(labels)}
            U = [[rng.randrange(p or 7) for _ in range(target)] for i in range(5)]
            V = [[rng.randrange(p or 7) for v in labels] for j in range(target)]
            columns = {v: clean({i: baseline[v].get(i, 0)+sum(U[i][k]*V[k][j] for k in range(target))
                                 for i in range(5)}, p) for j, v in enumerate(labels)}
            if rank(labels, columns, p) == 5:
                break
        cache = {}; patches = []
        for w in weight_controls(labels, rng):
            row, region = transfer(columns, baseline, groups, w, p, cache); patches.append(row)
            transfer(columns, baseline, groups, strict_weight(region, columns), p, cache)
            assert row['departure_rank'] <= target
        abstract.append({'field_characteristic': p, 'rank': 5, 'target_update_rank': target, 'patches': patches})
    files = {'docs/hook_rank_departure.md', 'scripts/run_hook_rank_departure_checks.py',
             'docs/hook_bridge_reduction.md', 'docs/hook_prefix_bridge_rank.md',
             'docs/two_colour_threshold_braids.md', 'docs/hook_high_price_partners.md',
             'results/research/hook_high_price_partner_checks.json'}
    pending = ['scripts/run_hook_rank_departure_checks.py']
    while pending:
        f = pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names = ([node.module] if isinstance(node, ast.ImportFrom) else
                     [v.name for v in node.names] if isinstance(node, ast.Import) else [])
            for name in names:
                target = 'scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:
                    files.add(target); pending.append(target)
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_hook_rank_departure_checks.py',
              'scope': 'Uniform whole-hook fan parameter bound from actual low-rank departure from one transverse arm. All compatible outer images and all fields allowed. No unrestricted exponent improvement or exact novelty claimed.',
              'uniform_result': {'fan_bound': 'O(n^2*(1+rho)^3)', 'combined_bridge_bound': 'O(n^2*(1+min(b,rho))^3)',
                                 'common_partition_chambers': 'O(n^2)', 'small_rank_per_chamber': '<=3*r',
                                 'small_original_labels_per_chamber': '<=9*r',
                                 'necessary_departure_for_polynomial_excess': 'rho=Omega(n^(epsilon/3))'},
              'counts': dict(sorted(COUNTS.items())), 'original_quotient_replay_counts': dict(REPLAY),
              'hook_cases': cases, 'abstract_matrix_controls': abstract,
              'limits': ['The low-rank update and common-quotient framework is classical.',
                         'The exact encoding depends on the actual coefficients and chamber.',
                         'No common weight sampler of the fan-count size is inferred.',
                         'Large departure and bridge ranks, unrestricted sharp growth, polynomial separation and exact priority remain open.'],
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_rank_departure_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'counts': report['counts'], 'source_hashes': len(files)}))


if __name__ == '__main__':
    main()
