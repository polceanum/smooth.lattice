#!/usr/bin/env python3
"""Exact joint-branch ranks, complete fans and sharp quadric-cross examples."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import factorial, gcd
from pathlib import Path
import hashlib
import json

from run_coordinate_slice_budget_checks import downset, edges_from_polygons, primitive, profile
from run_interleaving_research import minimal_excluded_border
from run_primitive_wall_checks import threshold
from run_signed_face_checks import field_basis, normals
from run_structural_complexity_checks import area2, clip, dot, finite_staircases

ROOT = Path(__file__).resolve().parents[1]


def normalized(v, p):
    return {i: c % p if p else Q(c) for i, c in v.items() if (c % p if p else c)}


def apply(M, v, p):
    out = defaultdict(int)
    for j, c in v.items():
        for i, d in M[j].items():
            out[i] += c*d
    return normalized(out, p)


def rank(columns, p):
    return len(field_basis(tuple(range(len(columns))), lambda i: columns[i], p)[0])


def point_model(points, p):
    points = tuple(tuple(a % p if p else a for a in v) for v in points)
    assert len(set(points)) == len(points)
    assert all(y*z == 0 for _, y, z in points)
    n = len(points)

    def ev(a):
        return normalized({i: (pow(x, a[0])*pow(y, a[1])*pow(z, a[2]))
                           for i, (x, y, z) in enumerate(points)}, p)

    M = [tuple(normalized({i: v[j]}, p) for i, v in enumerate(points)) for j in range(3)]
    return n, ev, M, p, (0, -1, 1), {'kind': 'radical_axis_points', 'points': points, 'field': p}


def gaussian_point_model():
    # Q(i)^3, with x=i and each component on one of the two axis branches.
    nodes = (((1, 1), (0, 0)), ((2, 1), (0, 0)), ((0, 0), (3, 2)))
    n = 2*len(nodes)

    def mul(a, b):
        return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]

    def power(a, k):
        out = (1, 0)
        for _ in range(k):
            out = mul(out, a)
        return out

    def ev(a):
        out = {}
        for j, (y, z) in enumerate(nodes):
            v = mul(mul(power((0, 1), a[0]), power(y, a[1])), power(z, a[2]))
            out.update({2*j+i: c for i, c in enumerate(v) if c})
        return normalized(out, 0)

    M = []
    for coord in range(3):
        columns = []
        for y, z in nodes:
            a, b = ((0, 1), y, z)[coord]
            j = len(columns)
            columns.extend((normalized({j: a, j+1: b}, 0),
                            normalized({j: -b, j+1: a}, 0)))
        M.append(tuple(columns))
    return n, ev, M, 0, (1, 0, 1), {'kind': 'irreducible_quadratic_Q_i', 'nodes': nodes, 'field': 0}


def gluing_model(lam, p):
    # Both pure tips are maximal; the two mixed relations use their full span.
    D0 = {(0, 0, 0)} | {(0, j, 0) for j in range(1, 4)} | {(0, 0, j) for j in range(1, 4)}
    D1 = {(1, 0, 0), (1, 1, 0), (1, 0, 1)}
    D = D0 | D1
    U = ((1, 1, 0), (1, 0, 1))
    V = ((0, 3, 0), (0, 0, 3))
    reference = sorted(D-set(U))
    index = {a: i for i, a in enumerate(reference)}

    def ev(a):
        if a == U[0]:
            return normalized({index[V[0]]: 1, index[V[1]]: 1}, p)
        if a == U[1]:
            return normalized({index[V[0]]: 1, index[V[1]]: lam}, p)
        return {index[a]: 1} if a in index else {}

    M = [tuple(ev(tuple(e+(i == j) for j, e in enumerate(a))) for a in reference)
         for i in range(3)]
    return len(reference), ev, M, p, (0, 0, 1), {
        'kind': 'nonreduced_two_tip_gluing', 'lambda': lam, 'field': p,
        'relations': ['xy=y^3+z^3', 'xz=y^3+lambda*z^3']}


def sharp_model(L):
    r = 2*L
    n = 2+2*r
    C = r*(1+r)
    lx, lt = C*L+1, C+r

    @lru_cache(None)
    def branch(a, b):
        v, node = L*a+b, lx*a+lt*b
        return normalized({i: Q(node**(i-v), factorial(i-v))
                           for i in range(v, r)}, 0) if v < r else {}

    def insert(v, start):
        return {start+i: c for i, c in v.items()}

    def ev(a):
        x, y, z = a
        if x >= 2 or y and z:
            return {}
        if y:
            return insert(branch(x, y-1), 2)
        if z:
            return insert(branch(x, z-1), 2+r)
        return {x: Q(1)}

    X, T = branch(1, 0), branch(0, 1)

    def multiply(v, shift):
        return {i+shift: c for i, c in v.items() if i+shift < r}

    M = []
    for coord in range(3):
        columns = [{} for _ in range(n)]
        if coord == 0:
            columns[0] = {1: Q(1)}
        elif coord == 1:
            columns[0], columns[1] = insert({0: Q(1)}, 2), insert(X, 2)
        else:
            columns[0], columns[1] = insert({0: Q(1)}, 2+r), insert(X, 2+r)
        for start in (2, 2+r):
            value = X if coord == 0 else T if (coord == 1 and start == 2) or (coord == 2 and start == 2+r) else {}
            for j in range(r):
                columns[start+j] = insert(multiply(value, j), start)
        M.append(tuple(columns))
    return n, ev, M, 0, (0, 0, 1), {
        'kind': 'sharp_two_branch_product', 'L': L, 'field': 0,
        'branch_length': r, 'C': C, 'lambda_x': lx, 'lambda_t': lt}


def carrier(n):
    return tuple(sorted({(0, 0, 0), (1, 0, 0)} |
                        {(i, j, 0) for i in range(2) for j in range(1, n)} |
                        {(i, 0, j) for i in range(2) for j in range(1, n)}))


def algebra_controls(model, counts):
    n, ev, M, p, f, _ = model
    for i, j in combinations(range(3), 2):
        for c in range(n):
            assert apply(M[i], M[j][c], p) == apply(M[j], M[i][c], p)
            counts['independent_commuting_matrix_columns'] += 1
    for c in range(n):
        x2 = apply(M[0], M[0][c], p)
        value = defaultdict(int, x2)
        for i, d in M[0][c].items():
            value[i] += f[1]*d
        value[c] += f[0]
        assert not normalized(value, p)
        assert not apply(M[1], M[2][c], p)
        counts['independent_quadratic_and_cross_annihilations'] += 1
    E = carrier(n)
    assert rank([ev(a) for a in E], p) == n
    for a in E:
        for i in range(3):
            successor = tuple(e+(i == j) for j, e in enumerate(a))
            assert apply(M[i], ev(a), p) == ev(successor)
            counts['independent_raw_column_successor_replays'] += 1
    Y, Z = list(M[1]), list(M[2])
    C = [ev((0, 0, 0)), ev((1, 0, 0))]
    mY, mZ = rank(Y, p), rank(Z, p)
    s = mY+mZ-rank(Y+Z, p)
    ker_y = mY-rank([apply(M[1], v, p) for v in Y], p)
    ker_z = mZ-rank([apply(M[2], v, p) for v in Z], p)
    dY = mY+rank(Z+C, p)-n
    assert 0 <= s <= min(ker_y, ker_z)
    assert max(ker_y, ker_z) <= 2
    assert 0 <= dY <= 4
    assert rank(Y+Z+C, p) == n
    counts['independent_joint_branch_intersection_bounds'] += 1
    return {'y_ideal_rank': mY, 'z_ideal_rank': mZ, 'branch_intersection_rank': s,
            'kernel_y_on_y_ideal': ker_y, 'kernel_z_on_z_ideal': ker_z,
            'y_intersection_z_plus_core': dY}


def shape_candidates(n, mY):
    for tY in range(max(0, mY-4), min(n-1, mY)+1):
        yield frozenset({(0, 0, 0)} |
                        {(0, j, 0) for j in range(1, tY+1)} |
                        {(0, 0, j) for j in range(1, n-tY)})
    for tY in range(max(0, mY-4), min(n-2, mY)+1):
        tZ = n-2-tY
        for b in range(tY//2+1):
            for d in range(tZ//2+1):
                a, c = tY-b, tZ-d
                yield frozenset({(0, 0, 0), (1, 0, 0)} |
                                {(0, j, 0) for j in range(1, a+1)} |
                                {(1, j, 0) for j in range(1, b+1)} |
                                {(0, 0, j) for j in range(1, c+1)} |
                                {(1, 0, j) for j in range(1, d+1)})


def full_fan(model, counts):
    n, ev, M, p, f, metadata = model
    controls = algebra_controls(model, counts)
    candidates = set(shape_candidates(n, controls['y_ideal_rank']))
    assert len(candidates) <= Q(5)+Q(5, 16)*(n+2)**2
    cells, independent, ledger = {}, set(), []
    for B in sorted(candidates, key=lambda B: sorted(B)):
        assert len(B) == n and downset(B)
        counts['joint_rank_shape_candidates'] += 1
        if len(field_basis(tuple(sorted(B)), ev, p)[0]) != n:
            continue
        independent.add(B)
        tY = sum(bool(a[1]) for a in B)
        assert controls['y_ideal_rank']-4 <= tY <= controls['y_ideal_rank']
        N, _ = normals(tuple(sorted(B)), ev, p)
        poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for q in N:
            poly = clip(poly, q, 1)
        if not area2(poly):
            continue
        assert len(minimal_excluded_border(tuple(B))) <= 6
        w = tuple(sum(v[i] for v in poly)/len(poly) for i in range(3))
        E = carrier(n)
        greedy = field_basis(sorted(E, key=lambda a: (dot(w, a), a)), ev, p)[0]
        assert frozenset(greedy) == B
        cells[B] = poly
        ledger.append({'standard_set': sorted(B), 'profile': profile(B),
                       'positive_cone': [[str(a) for a in v] for v in poly]})
        counts['independent_complete_cone_and_greedy_equalities'] += 1
    assert sum(map(area2, cells.values())) == 1
    edges = edges_from_polygons(cells, counts)
    total = 0
    for pair in edges:
        B, C = tuple(pair)
        q = primitive(tuple(a-b for a, b in zip(profile(B), profile(C))))
        assert threshold(q) <= n
        total += threshold(q)
    assert total <= 3*n*len(cells)
    assert len(cells) <= len(independent) <= Q(5)+Q(5, 16)*(n+2)**2
    if n <= 8:
        exhaustive = set()
        for B in finite_staircases(n, 3):
            B = frozenset(B)
            if any(a[0] >= 2 or a[1] and a[2] for a in B):
                continue
            if len(field_basis(tuple(sorted(B)), ev, p)[0]) == n:
                exhaustive.add(B)
            counts['independent_exhaustive_carrier_downsets'] += 1
        assert exhaustive == independent
        counts['independent_exhaustive_rank_window_equalities'] += 1
    if metadata['kind'] == 'nonreduced_two_tip_gluing':
        assert controls['branch_intersection_rank'] == 2
    if metadata['kind'] == 'sharp_two_branch_product':
        L = metadata['L']
        assert len(cells) == len(independent) == (L+1)**2
        assert len(edges) == 2*L*(L+1)
        assert total == 2*L*(L+1)**2
        for B in ({(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)},
                  {(0, 0, 0), (1, 0, 0), (0, 0, 1), (1, 0, 1)}):
            assert len(field_basis(tuple(sorted(B)), ev, p)[0]) == 4
            counts['independent_full_rank_other_pair_boxes'] += 1
        box = tuple(product(range(2), range(L), range(2*L+1)))
        assert rank([ev(a) for a in box], p) == n-2
        controls['quadratic_nonsaturating_box'] = {'side_counts': [2, L, 2*L+1],
                                                 'cardinality': len(box), 'rank': n-2}
        counts['independent_superlinear_nonsaturating_boxes'] += 1
    counts['complete_actual_quadric_cross_fans'] += 1
    return {'length': n, **metadata, 'joint_rank_controls': controls,
            'shape_candidates': len(candidates), 'order_ideal_vector_bases': len(independent),
            'initial_ideals': len(cells), 'undirected_walls': len(edges),
            'total_primitive_wall_cost': total, 'complete_area': '1', 'cells': ledger}


def main():
    models = [
        point_model(((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 2, 0),
                     (0, 0, 3), (1, 0, 4)), 0),
        point_model(((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0),
                     (0, 0, 1), (1, 0, 1)), 2),
        point_model(((0, 1, 0), (1, 2, 0), (0, 3, 0), (1, 0, 1),
                     (0, 0, 2), (1, 0, 3), (0, 0, 0)), 101),
        gaussian_point_model(),
    ]
    models += [gluing_model(lam, p) for p in (0, 2, 3, 101) for lam in (1, 2)]
    models += [sharp_model(L) for L in range(1, 7)]
    counts, rows = Counter(), []
    for model in models:
        row = full_fan(model, counts)
        rows.append(row)
        print('complete quadric cross:', row['kind'], row['length'],
              row['initial_ideals'], row['total_primitive_wall_cost'], flush=True)
    sources = ['docs/quadric_cross_fans.md', 'scripts/run_quadric_cross_checks.py',
               'docs/nested_jet_matroids.md', 'docs/wall_cost_products.md',
               'docs/primitive_wall_spectrum.md', 'docs/bilinear_rank_defects.md',
               'scripts/run_coordinate_slice_budget_checks.py',
               'scripts/run_signed_face_checks.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py']
    out = {'schema_version': 1, 'all_checks_passed': True,
           'scope': 'Uniform sharp class bounds for arbitrary finite quotients containing yz '
                    'and a monic quadratic in x, over every field. Complete finite fans test '
                    'joint ideal intersections, all-basis rank windows and sharp characteristic-zero '
                    'products. No unrestricted G3 result, polynomial lower separation or novelty established.',
           'uniform_theorem': {'initial_ideals': 'g(I)<=5+5*(n+2)^2/16',
                               'all_order_ideal_bases': 'The same bound, with no coherence assumption.',
                               'wall_cost': 'V(I)<=3*n*g(I)=O(n^3)',
                               'universal_basis': '|U(I)|<=6*g(I)=O(n^2)',
                               'sharp_family': 'N=4L+2, g=(L+1)^2, V=2L*(L+1)^2',
                               'joint_invariant': 'Y=yA, Z=zA; dim(Y intersection Z)<=2 '
                                                  'and dim(Y intersection (Z+C))<=4, dim C<=2.'},
           'counters': dict(counts), 'complete_actual_fans': rows,
           'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}}
    path = ROOT/'results/research/quadric_cross_checks.json'
    path.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'all_checks_passed': True, 'counters': dict(counts),
                      'report': str(path.relative_to(ROOT))}), flush=True)


if __name__ == '__main__':
    main()
