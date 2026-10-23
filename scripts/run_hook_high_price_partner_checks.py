#!/usr/bin/env python3
"""Original-field high-price exchanges and independent Toeplitz controls."""
from collections import Counter
from fractions import Fraction as Q
from math import comb, factorial
from pathlib import Path
import ast
import hashlib
import json

from run_hook_filtered_annihilator_checks import assemble
from run_hook_bridge_checks import normals, polygon, COUNTS as REPLAY
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import dot, area2
from run_oracle_enumeration_checks import exact_adjacent_weight

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def scalar(c, p):
    return c % p if p else Q(c)


def bareiss(matrix):
    """Integer determinant, independently of the represented-rank routines."""
    a = [list(row) for row in matrix]
    if not a:
        return 1
    old = 1; sign = 1
    for k in range(len(a)-1):
        pick = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pick is None:
            return 0
        if pick != k:
            a[k], a[pick] = a[pick], a[k]; sign = -sign
        pivot = a[k][k]
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                value = a[i][j]*pivot-a[i][k]*a[k][j]
                assert value % old == 0
                a[i][j] = value//old
            a[i][k] = 0
        old = pivot
    return sign*a[-1][-1]


def minor_product(D, h, r):
    out = Q(1)
    for j in range(h):
        out *= Q(factorial(D+j)*factorial(j),
                 factorial(r+j)*factorial(D-r+j))
    assert out.denominator == 1
    return int(out)


def strict_wall_weight(A, K, cols):
    h = A-1-K
    effective = {(A*v[0]+v[1], (1-h)*v[0]+v[2]) for v, c in cols.items() if c}
    lo = Q(K-1, A-3); hi = Q(K, A-3)
    for seed in range(1, 5000):
        u = ((seed+1)*lo+hi)/(seed+2)
        if len({u*b+c for b, c in effective}) == len(effective):
            return (A*u+1-h, u, Q(1))
    raise AssertionError('No generic wall witness found.')


def projected_quotient(data, cols, A, K, w, p):
    """A separately assembled row map for the complete cheaper quotient."""
    h = A-1-K; lam, u, _ = w; tau = A*u+3
    R = max(q for q in range(A) if 2*u+2+q < tau)
    assert h < K <= R < A-1
    G = [{q-K-1: scalar(comb(A-1, q-j), p)
          for q in range(K+1, A) if 0 <= q-j < A and scalar(comb(A-1, q-j), p)}
         for j in range(A)]
    dimension = h+A-R-1
    rows = {}
    for r, (y, z) in enumerate(data[0]):
        v = {}
        if y >= 2:
            if (y, z) == (A, 1):
                v = {i: scalar(-c, p) for i, c in G[0].items()}
        elif z >= 2:
            q = z-2
            if y == 1 and q > K:
                v = {q-K-1: scalar(1, p)}
            elif y == 0 and q >= 2:
                v = ({i: scalar(-c, p) for i, c in G[q].items()}
                     if q <= R else {h+q-R-1: scalar(1, p)})
        rows[r] = v

    images = {}
    for label, vector in cols.items():
        out = Counter()
        for r, c in vector.items():
            for i, v in rows[r].items():
                out[i] += c*v
        images[label] = {i: scalar(c, p) for i, c in out.items() if scalar(c, p)}
    prefix = list(range(1, h))
    _, piv = field_basis(prefix, G.__getitem__, p)
    cheap = [v for v in cols if dot(w, v) < tau]
    for v in cheap:
        assert not field_reduce(images[v], piv, p)[0]
        COUNTS['complete_cheaper_column_projection_identities'] += 1
    original_rank = len(field_basis(cheap, cols.__getitem__, p)[0])
    image_rank = len(field_basis(tuple(images), images.__getitem__, p)[0])
    assert image_rank == dimension and len(piv) == h-1
    assert 4*A-original_rank == image_rank-len(piv) == A-R
    COUNTS['independent_entire_cheaper_quotient_kernel_equalities'] += 1
    return cheap, G, original_rank, R


def check_wall(data, cols, A, K, p, order_cache):
    h = A-1-K; E = (0, A, 3); F = (1, 0, h+2)
    w = strict_wall_weight(A, K, cols); lam, u, _ = w; tau = dot(w, E)
    assert Q(1, 2) < u < 1 and lam > Q(A, 5)+2 > 2
    assert dot(w, F) == tau == A*u+3
    cheap, G, r, R = projected_quotient(data, cols, A, K, w, p)
    # These are the four separate layer/arm lists, not combined lists.
    order = tuple(tuple(sorted((v for v in data[4][arm] if v[0] == layer and cols[v]),
                               key=lambda v: dot((0, u, 1), v)))
                  for arm in (0, 1) for layer in (0, 1))
    if order_cache:
        assert order == order_cache[0]
    else:
        order_cache.append(order)
    COUNTS['fixed_four_planar_list_order_controls'] += 1
    for start, offset in ((1, K), (0, K+1)):
        matrix = [[comb(A-1, q-j) if 0 <= q-j < A else 0
                   for j in range(start, start+h)] for q in range(K+1, A)]
        determinant = bareiss(matrix)
        assert determinant == minor_product(A-1, h, offset) > 0
        assert not p or determinant % p
        COUNTS['independent_integer_determinants_vs_classical_product'] += 1
    assert len(field_basis(range(1, h), G.__getitem__, p)[0]) == h-1
    assert len(field_basis(range(h), G.__getitem__, p)[0]) == h
    assert len(field_basis(range(1, h+1), G.__getitem__, p)[0]) == h
    _, piv = field_basis(cheap, cols.__getitem__, p)
    e = field_reduce(cols[E], piv, p)[0]
    f = field_reduce(cols[F], piv, p)[0]
    assert e and f and len(field_basis((E, F), lambda v: e if v == E else f, p)[0]) == 1
    COUNTS['original_full_cheaper_nonloop_parallel_certificates'] += 1
    # Extract a separating functional from the ORIGINAL cheap-row quotient.
    # Its reversed Y leading label must be the same dual exception q.
    free_row = next(iter(e)); reflected = {}
    for row, (y, z) in enumerate(data[0]):
        if y < 2:
            continue
        coefficient = field_reduce({row: scalar(1, p)}, piv, p)[0].get(free_row, 0)
        if coefficient:
            reflected[A-1-(y-2), 1-z] = coefficient
    assert reflected == {(1, 0): e[free_row]}
    assert A-1 >= A-2 and scalar(e[free_row], p)
    COUNTS['independent_original_functionals_with_one_fixed_dual_exception'] += 1
    costs = sorted({dot(w, v) for v, c in cols.items() if c})
    gap = min(b-a for a, b in zip(costs, costs[1:]))
    delta = min(gap/4, (lam-2)/4)
    bases = []
    for sign in (-1, 1):
        ww = (lam+sign*delta, u, Q(1))
        B = frozenset(field_basis(sorted(cols, key=lambda v: (dot(ww, v), v)),
                                 cols.__getitem__, p)[0])
        assert len(B) == 4*A
        bases.append(B)
    assert bases[0]-bases[1] == {F} and bases[1]-bases[0] == {E}
    q = (1, -A, h-1)
    assert tuple(-v for v in q) in normals(bases[0], cols, p)
    assert q in normals(bases[1], cols, p)
    COUNTS['independent_adjacent_original_basis_and_wall_normal_checks'] += 1
    return {'K': K, 'partner_h': h, 'u': str(u), 'lambda': str(lam),
            'lambda_over_scalar_high_boundary': str(lam/2),
            'direction': q, 'cheaper_rank': r, 'last_cheap_Z_lower_head': R}


def complete_fan(cols, A, p):
    n = 8*A+8; queue = [((1, 2*n+1, (2*n+1)**2), None)]; cells = {}
    while queue:
        w, old = queue.pop()
        B = frozenset(field_basis(sorted(cols, key=lambda v: (dot(w, v), v)),
                                 cols.__getitem__, p)[0])
        if old is not None:
            assert B != old
        if B in cells:
            continue
        N = normals(B, cols, p); poly = polygon(N)
        assert area2(poly) > 0 and all(dot(w, q) > 0 for q in N)
        cells[B] = poly
        for left, right in zip(poly, poly[1:]+poly[:1]):
            if any(left[i] == right[i] == 0 for i in range(3)):
                continue
            q = next(q for q in N if dot(q, left) == dot(q, right) == 0)
            nw, _ = exact_adjacent_weight(n, left, right, q)
            queue.append((nw, B))
    assert sum(map(area2, cells.values())) == 1
    COUNTS['independent_complete_original_fans'] += 1
    COUNTS['independent_complete_original_fan_cones'] += len(cells)
    return {'initial_cones': len(cells), 'area2': '1'}


def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    previous = json.loads((ROOT/'results/research/hook_passive_band_checks.json').read_text())
    assert previous['status'] == 'passed'
    for f, sha in previous['source_sha256'].items():
        assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest() == sha
    cases = []
    for A, p, full in ((20, 101, True), (24, 0, True), (36, 1000003, False),
                       (48, 101, False), (64, 1000003, False)):
        assert not p or p > 2*A
        units = [{(0, 0): 1, (A-2, 1): 1},
                 {(0, 0): 1, **{(q, 1): comb(A-1, q) for q in range(A)}}]
        data, co, cols = assemble(2, A+2, units, p)
        operator = {(q, v): ({(q+A-2, v+1): scalar(1, p)}
                             if q+A-2 < A and v+1 < 2 else {})
                    for q in range(A) for v in range(2)}
        plus = field_basis(tuple(operator), operator.__getitem__, p)[0]
        truncated = {v: {r: c for r, c in vector.items() if r[0] < A-2}
                     for v, vector in operator.items()}
        minus = field_basis(tuple(truncated), truncated.__getitem__, p)[0]
        assert set(plus) == {(0, 0), (1, 0)} and not minus
        COUNTS['independent_two_label_dual_exception_sets'] += 1
        order_cache = []
        walls = [check_wall(data, cols, A, K, p, order_cache)
                 for K in range((3*A+4)//5, 3*A//4+1)]
        assert len({tuple(v['direction']) for v in walls}) == len(walls)
        assert len(walls) >= Q(3*A, 20)-1
        fan = complete_fan(cols, A, p) if full else None
        cases.append({'A': A, 'length': 8*A+8, 'characteristic': p,
                      'distinct_partners_certified': len(walls), 'walls': walls,
                      'complete_original_fan': fan})
        print('High-price case passed', A, p, 'partners', len(walls),
              'full cones', fan['initial_cones'] if fan else None, flush=True)
    # The stated characteristic condition is material to the minor proof.
    bad = minor_product(23, 8, 15)
    assert bad % 2 == 0
    COUNTS['small_characteristic_minor_degeneration_controls'] += 1
    files = {'docs/hook_high_price_partners.md', 'scripts/run_hook_high_price_partner_checks.py',
             'docs/hook_bridge_reduction.md', 'docs/hook_exception_partners.md',
             'docs/hook_clipped_threshold_budget.md', 'docs/hook_passive_band_fans.md',
             'results/research/hook_passive_band_checks.json'}
    pending = ['scripts/run_hook_high_price_partner_checks.py']
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
              'command': 'python3 -B scripts/run_hook_high_price_partner_checks.py',
              'scope': 'Sharp linear one-column partner growth in an actual thickness-two constant-axis central hook at prices of order n. No unrestricted count improvement, polynomial fan separation or exact novelty claimed.',
              'uniform_result': {'length': '8*A+8', 'field': 'characteristic zero or p>2*A',
                                 'partner_lower_bound': 'ceil(3*A/5)<=K<=floor(3*A/4), at least 3*A/20-1 distinct partners',
                                 'original_label': [0, 'A', 3], 'partner': [1, 0, 'A+1-K'],
                                 'interval': '(K-1)/(A-3)<u<K/(A-3)',
                                 'wall': 'lambda=A*u+1-h, h=A-1-K',
                                 'price': 'lambda>A/5+2>2*max(u,1)',
                                 'direction': '(1,-A,h-1)', 'four_planar_orders': 'fixed on 1/2<u<1'},
              'counts': dict(sorted(COUNTS.items())), 'original_quotient_replay_counts': dict(REPLAY),
              'cases': cases,
              'limits': ['The determinant evaluation and greedy exchange criterion are classical.',
                         'The full fan is already O(n^2) by the saved passive-band theorem.',
                         'This refutes extending a local low-price partner cap, not a collective full-fan bound.',
                         'Exact priority and the original unrestricted objective remain unresolved.'],
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_high_price_partner_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'counts': report['counts'], 'source_hashes': len(files)}))


if __name__ == '__main__':
    main()
