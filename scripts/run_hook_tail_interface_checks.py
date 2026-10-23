#!/usr/bin/env python3
"""Far-arm bridge projections and rational recurrence interface checks."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_hook_bridge_checks import model, actual, strict_weight, COUNTS as REPLAY
from run_hook_rank_departure_checks import transfer, weight_controls, rank, clean, plane, complete_fan, comparison, COUNTS as TRANSFER
from run_hook_transverse_prefix_checks import mul
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import dot

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def auxiliary(data, columns, t, arm):
    far = {i for i, v in enumerate(data[0]) if v[arm] >= 2*t}
    H = set(data[4][2])
    base = {v: {i: c for i, c in vector.items() if not (v in H and i in far)}
            for v, vector in columns.items()}
    projected = {v: {i: c for i, c in columns[v].items() if i in far} for v in H}
    return base, far, projected


def split_controls(data, columns, base, far, t, a, arm, w, p):
    far_labels = tuple(v for v, c in base.items() if c and set(c) <= far)
    core = tuple(v for v, c in base.items() if c and v not in far_labels)
    assert set(far_labels) <= set(data[4][arm])
    assert all(not (set(base[v]) & far) for v in core)
    left = field_basis(sorted(far_labels, key=lambda v: dot(w, v)), base.__getitem__, p)[0]
    right = field_basis(sorted(core, key=lambda v: dot(w, v)), base.__getitem__, p)[0]
    whole = field_basis(sorted(base, key=lambda v: dot(w, v)), base.__getitem__, p)[0]
    assert set(left)|set(right) == set(whole)
    COUNTS['independent_auxiliary_direct_sum_minimum_equalities'] += 1
    # Use ALL possible supports, not just columns nonzero in this coefficient case.
    schema = {v: {i: 1 for i in data[3][v]} for v in columns}
    possible, _, _ = auxiliary(data, schema, t, arm)
    near = tuple(v for v, c in possible.items() if c and not set(c) <= far)
    for v in near:
        x, y, z = v if arm == 0 else (v[0], v[2], v[1])
        assert (x in (0, 1) and y < 2*t and z < 2*a) or (x == 0 and a <= y < a+t and z < t)
    inventory = {plane(tuple(x-y for x, y in zip(v, u)))
                 for j, v in enumerate(near) for u in near[:j]}
    for g in data[4]:
        inventory |= {plane(tuple(x-y for x, y in zip(v, u))) for j, v in enumerate(g) for u in g[:j]}
    assert len(inventory) <= 1024*a*t
    COUNTS['coefficient_independent_possible_support_plane_inventories'] += 1
    # An erased bridge matrix is allowed to fail the original module relation.
    # Record such a witness explicitly instead of treating it as an algebra.
    index = {v: i for i, v in enumerate(data[0])}
    for label in data[4][2]:
        neighbor = tuple(v+int(j == arm+1) for j, v in enumerate(label))
        if neighbor not in base:
            continue
        shifted = Counter()
        for i, c in base[label].items():
            v = tuple(x+int(j == arm) for j, x in enumerate(data[0][i]))
            if v in index:
                shifted[index[v]] += c
        if clean(shifted, p) != base[neighbor]:
            COUNTS['explicit_auxiliary_noncommuting_source_controls'] += 1
            break
    return len(inventory)


def series_operator(A, t, delta, denominator, numerator, C, p):
    q0 = denominator[0] % p if p else Q(denominator[0])
    assert q0
    inv = pow(int(q0), -1, p) if p else 1/q0
    values = {}
    for i in range(A):
        for j in range(t-delta):
            c = numerator.get((i, j), 0)-sum(denominator[k]*values.get((i-k, j), 0)
                                             for k in range(1, min(i+1, len(denominator))))
            c = c*inv
            if p:
                c %= p
            if c:
                values[i, j] = c
    N = {(i, j+delta): c for (i, j), c in values.items()}
    W = Counter(N)
    W.update({(0, j): c for j, c in enumerate(C)})
    return clean(W, p), clean(N, p)


def rational_block(data, columns, t, a, arm, delta, d, denominator, numerator, N, p):
    A = a-t; Qpoly = {(i, 0): c for i, c in enumerate(denominator) if c}
    Ppoly = clean({(i, j+delta): c for (i, j), c in numerator.items() if i < A and j+delta < t}, p)
    assert mul(Qpoly, N, A, t, p) == Ppoly
    base, far, projection = auxiliary(data, columns, t, arm)
    blocks = {}; rhs_columns = {}
    for i, j in product(range(t), repeat=2):
        full = mul(N, {(i, j): 1}, A, t, p)
        low = {v: c for v, c in full.items() if v[0] < t}
        high = {v: c for v, c in full.items() if v[0] >= t}
        lhs = mul(Qpoly, high, A, t, p)
        rhs = Counter(mul(Ppoly, {(i, j): 1}, A, t, p))
        rhs.subtract(mul(Qpoly, low, A, t, p))
        rhs = clean({v: c for v, c in rhs.items() if v[0] >= t}, p)
        assert lhs == rhs
        assert all(t <= q < t+d and v >= delta for q, v in rhs)
        blocks[i, j] = high; rhs_columns[i, j] = rhs
        label = (0, t+i, t+j) if arm == 0 else (0, t+j, t+i)
        original = {(data[0][row][arm]-t, data[0][row][1-arm]): c
                    for row, c in projection[label].items()}
        assert original == high
        COUNTS['original_bridge_vs_rational_cross_cut_columns'] += 1
        COUNTS['independent_full_rational_block_equations'] += 1
    chi = rank(tuple(projection), projection, p)
    assert rank(tuple(blocks), blocks, p) == rank(tuple(rhs_columns), rhs_columns, p) == chi
    assert chi <= min(t*t, t*(A-t), d*(t-delta))
    COUNTS['independent_rational_cut_rank_bounds'] += 1
    return chi


def geometric_controls(data, co, columns, t, a, p):
    A = a-t; expected = (A-1)*(t-1); candidates = []
    if p and p**t <= 32:
        candidates = list(product(range(p), repeat=t))
        exhaustive = True
    else:
        candidates = [tuple([1, 1]+[0]*(t-2)), tuple([0]*t), tuple([1]+[0]*(t-1)),
                      tuple([1, 1]+[1]*(t-2))]
        exhaustive = False
    for arm in (0, 1):
        ranks = []
        for C in candidates:
            base = actual(data, comparison(data, co, arm, C), p)
            diff = {}
            for v in columns:
                c = Counter(columns[v]); c.subtract(base[v]); diff[v] = clean(c, p)
            value = rank(tuple(diff), diff, p)
            assert value >= expected
            ranks.append(value)
            COUNTS['geometric_transverse_departure_rank_controls'] += 1
        assert min(ranks) == expected
    for m in (A, A-t):
        if m < 2:
            continue
        regular = []
        for r in range(m+1):
            h = m-r
            matrix = {j: {i: 1 for i in range(h) if 0 <= r+i-j < m} for j in range(h)}
            if rank(tuple(matrix), matrix, p) == h:
                regular.append(r)
        assert set(regular) == {0, m-1, m}
        COUNTS['independent_maximal_geometric_prefix_gap_controls'] += 1
    assert rank(data[4][2], columns, p) == t*t
    return {'least_departure_proved': expected, 'finite_comparisons_exhaustive': exhaustive,
            'comparison_polynomials_per_arm_checked': len(candidates), 'bridge_rank': t*t}


def run_case(t, a, p, kind, arm, rng, spec=None, full=False):
    data = model(t, a); A = a-t; rational = spec is not None
    if rational:
        delta, d, denominator, numerator, C = spec
        W, N = series_operator(A, t, delta, denominator, numerator, C, p)
    co = []
    for g, y, z in data[2]:
        long, short = (y, z) if arm == 0 else (z, y)
        if rational and long >= t:
            c = W.get((long-t, short), 0) if g == 1 else 0
        elif kind == 'geometric':
            c = int(g == 1 and ((y, z) in ((t, 0), (0, t)) or
                               (y >= t and z == 1) or (z >= t and y == 1)))
        else:
            c = rng.randrange(p or 7)
        co.append(c)
    columns = actual(data, co, p, True)
    base, far, projected = auxiliary(data, columns, t, arm)
    chi = rank(tuple(projected), projected, p)
    assert chi <= rank(data[4][2], columns, p)
    if rational:
        assert chi == rational_block(data, columns, t, a, arm, delta, d, denominator, numerator, N, p)
    if kind == 'geometric':
        assert chi == t-1
        geom = geometric_controls(data, co, columns, t, a, p)
    else:
        geom = None
    cache = {}; patches = []; weights = weight_controls(columns, rng)
    for w in weights:
        row, region = transfer(columns, base, data[4], w, p, cache)
        assert row['departure_rank'] == chi; patches.append(row)
        row, _ = transfer(columns, base, data[4], strict_weight(region, columns), p, cache)
        patches.append(row)
    inventory = split_controls(data, columns, base, far, t, a, arm, weights[-1], p)
    fan = complete_fan(columns, 4*a*t-2*t*t, p) if full else None
    if fan:
        for B, poly in fan.items():
            w = strict_weight(poly, columns)
            row, _ = transfer(columns, base, data[4], w, p, cache); patches.append(row)
            assert frozenset(field_basis(sorted(columns, key=lambda v: dot(w, v)), columns.__getitem__, p)[0]) == B
        COUNTS['complete_original_fan_cones_with_interface_patches'] += len(fan)
    print('Interface case passed', kind, t, a, p, arm, 'chi', chi,
          'full cones', len(fan) if fan else None, flush=True)
    return {'kind': kind, 't': t, 'a': a, 'characteristic': p, 'arm': 'Y' if arm == 0 else 'Z',
            'tail_interface_rank': chi, 'original_bridge_rank': rank(data[4][2], columns, p),
            'uniform_plane_inventory': inventory, 'patches': patches, 'geometric_controls': geom,
            'complete_original_fan_cones': len(fan) if fan else None}


def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    previous = json.loads((ROOT/'results/research/hook_rank_departure_checks.json').read_text())
    assert previous['status'] == 'passed'
    for f, sha in previous['source_sha256'].items():
        assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest() == sha
    rng = Random(2026102322); cases = []
    for t, a, p, arm, full in ((2, 5, 101, 0, True), (3, 8, 2, 1, False),
                               (3, 9, 0, 0, False), (3, 6, 3, 0, False)):
        cases.append(run_case(t, a, p, 'arbitrary_outer', arm, rng, full=full))
    for t, a, p, arm, full in ((2, 7, 2, 0, True), (3, 10, 2, 1, False),
                               (4, 14, 2, 0, False), (3, 8, 0, 1, False)):
        cases.append(run_case(t, a, p, 'geometric', arm, rng, full=full))
    specs = [(3, 10, 3, 0, (1, 1, [1, -1], {(0, 0): 1, (1, 1): 1}, [0, 1, 1])),
             (4, 13, 101, 1, (1, 2, [1, 1, 1], {(0, 0): 1, (1, 1): 2, (2, 2): 1}, [1, 0, 1, 0])),
             (5, 17, 2, 0, (2, 2, [1, 1, 1], {(0, 0): 1, (1, 1): 1, (2, 2): 1}, [0, 1, 1, 0, 0])),
             (3, 9, 0, 1, (0, 2, [2, 1, -1], {(0, 0): 1, (1, 1): 2, (2, 2): 1}, [0, 1, 0])),
             (3, 6, 2, 0, (1, 1, [1, -1], {(0, 0): 1}, [1, 0, 0])),
             (3, 10, 101, 1, (1, 0, [2], {(0, 0): 1, (0, 1): 1}, [1, 0, 1]))]
    for t, a, p, arm, spec in specs:
        cases.append(run_case(t, a, p, 'rational', arm, rng, spec=spec))
    assert COUNTS['explicit_auxiliary_noncommuting_source_controls']
    files = {'docs/hook_tail_interface.md', 'scripts/run_hook_tail_interface_checks.py',
             'docs/hook_rank_departure.md', 'docs/hook_bridge_reduction.md',
             'docs/two_colour_threshold_braids.md', 'results/research/hook_rank_departure_checks.json'}
    pending = ['scripts/run_hook_tail_interface_checks.py']
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
              'command': 'python3 -B scripts/run_hook_tail_interface_checks.py',
              'scope': 'Uniform whole-hook far-arm interface-rank bound, with an all-field rational-recurrence specialization. Auxiliary erased matrices are not asserted to be algebras. No unrestricted exponent or exact novelty claimed.',
              'uniform_result': {'fan_bound': 'O(n^2*(1+min(chi_Y,chi_Z))^3)',
                                 'combined_departure_bound': 'O(n^2*(1+min(chi_Y,chi_Z,rho))^3)',
                                 'rational_interface_bound': 'chi_Y<=min(t^2,t*(A-t),d*(t-delta))',
                                 'geometric_family': 'b=t^2, chi=t-1, rho=(A-1)*(t-1)',
                                 'necessary_interface_for_polynomial_excess': 'chi_Y,chi_Z=Omega(n^(epsilon/3))'},
              'counts': dict(sorted(COUNTS.items())), 'transfer_counts': dict(sorted(TRANSFER.items())),
              'original_quotient_replay_counts': dict(REPLAY), 'cases': cases,
              'limits': ['Common-quotient compression and rational/Hankel rank are classical inputs.',
                         'The exact small minor and compulsory set depend on coefficients and chamber.',
                         'The rational specialization requires vanishing outer components on the chosen arm.',
                         'Large interface and departure ranks, unrestricted sharp growth, polynomial separation and exact priority remain open.'],
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_tail_interface_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'counts': report['counts'],
                      'transfer_counts': report['transfer_counts'], 'source_hashes': len(files)}))


if __name__ == '__main__':
    main()
