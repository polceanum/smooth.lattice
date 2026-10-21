#!/usr/bin/env python3
"""Exact Laurent rank translation, parity paths and complete quadratic fans."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import comb, factorial, gcd
from pathlib import Path
import hashlib
import json

from run_coordinate_slice_budget_checks import downset, edges_from_polygons, primitive, profile
from run_quadric_cross_checks import apply, carrier, normalized, rank
from run_primitive_wall_checks import threshold
from run_signed_face_checks import field_basis, normals
from run_structural_complexity_checks import area2, clip, dot, finite_staircases

ROOT = Path(__file__).resolve().parents[1]


def series_model(n, p=0, exponential=True):
    q = (n+1)//2
    if exponential:
        invfact = [pow(factorial(i), -1, p) if p else Q(1, factorial(i)) for i in range(n)]

        @lru_cache(None)
        def ev(a):
            v, node = q*a[0], 3*n*a[0]+a[1]-a[2]
            return normalized({i: node**(i-v)*invfact[i-v] for i in range(v, n)}, p) if v < n else {}
    else:
        X = {q: 1} if q < n else {}
        Y = {i: 1 for i in (0, 1) if i < n}
        Z = {i: (-1)**i for i in range(n)}

        def multiply(A, B):
            out = defaultdict(int)
            for i, c in A.items():
                for j, d in B.items():
                    if i+j < n:
                        out[i+j] += c*d
            return normalized(out, p)

        @lru_cache(None)
        def ev(a):
            out = {0: 1}
            for value, exponent in zip((X, Y, Z), a):
                for _ in range(exponent):
                    out = multiply(out, value)
            return normalized(out, p)

    values = [ev(tuple(int(i == j) for i in range(3))) for j in range(3)]
    M = [tuple({i+j: c for i, c in value.items() if i+j < n} for j in range(n))
         for value in values]
    return n, ev, M, p, (0, 0, 1), {
        'kind': 'sharp_exponential_Laurent' if exponential else 'special_polynomial_Laurent',
        'length': n, 'field': p, 'x_valuation': q, 'upper_node_shift': 3*n if exponential else 0}


def unit_points(nodes, p):
    points = [(x, y, pow(y, -1, p) if p else Q(1, y)) for x, y in nodes]
    n = len(points)

    def ev(a):
        return normalized({i: x**a[0]*y**a[1]*z**a[2]
                           for i, (x, y, z) in enumerate(points)}, p)

    M = [tuple(normalized({i: value[j]}, p) for i, value in enumerate(points)) for j in range(3)]
    return n, ev, M, p, (0, -1, 1), {
        'kind': 'radical_unit_pair_points', 'nodes_x_y': nodes, 'field': p}


def linear_points():
    points = ((0, 2, 1), (1, 3, 2), (0, 4, 3), (1, 5, 4))
    assert all(y-z == 1 and x*x == x for x, y, z in points)

    def ev(a):
        return normalized({i: x**a[0]*y**a[1]*z**a[2]
                           for i, (x, y, z) in enumerate(points)}, 0)

    M = [tuple(normalized({i: value[j]}, 0) for i, value in enumerate(points))
         for j in range(3)]
    return len(points), ev, M, 0, (0, -1, 1), {
        'kind': 'radical_linear_pair_points', 'points': points, 'field': 0}


def gaussian_units():
    nodes = ((1, 1), (2, 1), (3, 2))

    def mul(A, B):
        return A[0]*B[0]-A[1]*B[1], A[0]*B[1]+A[1]*B[0]

    def power(A, k):
        out = (Q(1), Q(0))
        for _ in range(k):
            out = mul(out, A)
        return out

    triples = [((0, 1), A, (Q(A[0], A[0]**2+A[1]**2), Q(-A[1], A[0]**2+A[1]**2)))
               for A in nodes]

    def ev(a):
        out = {}
        for j, triple in enumerate(triples):
            v = (Q(1), Q(0))
            for A, k in zip(triple, a):
                v = mul(v, power(A, k))
            out.update({2*j+i: c for i, c in enumerate(v) if c})
        return out

    M = []
    for coordinate in range(3):
        columns = []
        for triple in triples:
            a, b = triple[coordinate]
            j = len(columns)
            columns.extend((normalized({j: a, j+1: b}, 0),
                            normalized({j: -b, j+1: a}, 0)))
        M.append(tuple(columns))
    return 6, ev, M, 0, (1, 0, 1), {'kind': 'irreducible_quadratic_units_Q_i',
                                    'y_nodes': nodes, 'field': 0}


def algebra_controls(model, counts):
    n, ev, M, p, f, _ = model
    for i, j in combinations(range(3), 2):
        for c in range(n):
            assert apply(M[i], M[j][c], p) == apply(M[j], M[i][c], p)
            counts['independent_commuting_matrix_columns'] += 1
    for c in range(n):
        out = defaultdict(int, apply(M[0], M[0][c], p))
        for i, d in M[0][c].items():
            out[i] += f[1]*d
        out[c] += f[0]
        assert not normalized(out, p)
        assert apply(M[1], M[2][c], p) == {c: 1}
        counts['independent_quadratic_and_inverse_columns'] += 1
    E = carrier(n)
    assert rank([ev(a) for a in E], p) == n
    for a in E:
        for i in range(3):
            b = tuple(e+(i == j) for j, e in enumerate(a))
            assert apply(M[i], ev(a), p) == ev(b)
            counts['independent_raw_monomial_successor_replays'] += 1


def basis_shape(n, u, delta, d):
    ell = n-u
    b, a, c = delta+d, ell-1-delta-d, u-1-d
    B = {(0, j, 0) for j in range(a+1)} | {(0, 0, j) for j in range(1, b+1)}
    if u:
        B |= {(1, j, 0) for j in range(c+1)} | {(1, 0, j) for j in range(1, d+1)}
    B = frozenset(B)
    assert len(B) == n and downset(B)
    return B


def parity_path(n, u, D, parity, counts):
    M = n-2*u
    points = sorted((delta+2*d, delta) for delta in D if delta % 2 == parity for d in range(u))
    if not points:
        return [], set()
    P = set(points)
    for A, B in combinations(points, 2):
        assert tuple(map(min, zip(A, B))) in P and tuple(map(max, zip(A, B))) in P
        counts['independent_parity_sublattice_pairs'] += 1
    coefficients = {A: ((A[0]+1)**2+A[1]**2, -2*(n*(A[0]+1)+M*A[1])) for A in points}
    events = {Q(0), Q(1)}
    for A, B in combinations(points, 2):
        a, b = coefficients[A]
        c, d = coefficients[B]
        if b != d:
            t = Q(c-a, b-d)
            if 0 < t < 1:
                events.add(t)
    ordered = sorted(events)
    path = []
    for left, right in zip(ordered, ordered[1:]):
        t = (left+right)/2
        costs = {A: a+b*t for A, (a, b) in coefficients.items()}
        low = min(costs.values())
        winners = [A for A, c in costs.items() if c == low]
        A = (min(P[0] for P in winners), min(P[1] for P in winners))
        assert A in winners
        if not path or path[-1] != A:
            path.append(A)
    assert all(A[0] <= B[0] and A[1] <= B[1] for A, B in zip(path, path[1:]))
    assert len(path) <= n-u
    counts['independent_monotone_parity_paths'] += 1
    return path, set(path)


def interval_bases(model, counts):
    n, ev, _, p, _, _ = model
    candidates, controls, support = {}, [], {}
    r0 = rank([ev((0, j, 0)) for j in range(n)], p)
    for b in range(n):
        B = basis_shape(n, 0, b, 0)
        assert rank([ev(a) for a in B], p) == r0
        if r0 == n:
            candidates[B] = (0, b, 0)
    for u in range(1, n//2+1):
        ell, M = n-u, n-2*u
        D = []
        for delta in range(M+1):
            normal = [(0, j, 0) for j in range(ell)] + [(1, j, 0) for j in range(delta, delta+u)]
            rk = rank([ev(a) for a in normal], p)
            for d in range(u):
                B = basis_shape(n, u, delta, d)
                assert rank([ev(a) for a in B], p) == rk
                counts['independent_Laurent_rank_translations'] += 1
                if rk == n:
                    candidates[B] = (u, delta, d)
            if rk == n:
                D.append(delta)
        paths = []
        for parity in (0, 1):
            path, exposed = parity_path(n, u, D, parity, counts)
            support[u, parity] = exposed
            paths.append(path)
        controls.append({'upper_count': u, 'feasible_offsets': D, 'parity_paths': paths})
    return candidates, controls, support


def extract_fan(model, candidates, counts, parity_support=None):
    n, ev, _, p, _, metadata = model
    cells, ledger = {}, []
    for B, parameters in sorted(candidates.items(), key=lambda item: sorted(item[0])):
        assert len(field_basis(tuple(sorted(B)), ev, p)[0]) == n
        N, _ = normals(tuple(sorted(B)), ev, p)
        poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for q in N:
            poly = clip(poly, q, 1)
        if not area2(poly):
            continue
        w = tuple(sum(v[i] for v in poly)/len(poly) for i in range(3))
        greedy = field_basis(sorted(carrier(n), key=lambda a: (dot(w, a), a)), ev, p)[0]
        assert frozenset(greedy) == B
        if parity_support is not None:
            u, delta, d = parameters
            if u:
                assert (delta+2*d, delta) in parity_support[u, delta % 2]
                t = w[1]/(w[1]+w[2])
                ell, M, s = n-u, n-2*u, delta+2*d
                b, a, c = delta+d, ell-1-delta-d, u-1-d
                raw = t*(a*(a+1)+c*(c+1))/2+(1-t)*(b*(b+1)+d*(d+1))/2
                shifted = ((s+1-n*t)**2+(delta-M*t)**2)/4
                # The omitted term depends only on n,u,t; compare another point.
                for P, param in candidates.items():
                    v, e, h = param
                    if v != u:
                        continue
                    bb, aa, cc, ss = e+h, ell-1-e-h, u-1-h, e+2*h
                    raw2 = t*(aa*(aa+1)+cc*(cc+1))/2+(1-t)*(bb*(bb+1)+h*(h+1))/2
                    shifted2 = ((ss+1-n*t)**2+(e-M*t)**2)/4
                    assert raw-raw2 == shifted-shifted2
                    assert P == B or raw < raw2
                counts['independent_unique_parity_cost_optima'] += 1
        cells[B] = poly
        ledger.append({'standard_set': sorted(B), 'parameters': parameters, 'profile': profile(B),
                       'positive_cone': [[str(a) for a in v] for v in poly]})
        counts['independent_complete_cone_greedy_replays'] += 1
    assert sum(map(area2, cells.values())) == 1
    edges = edges_from_polygons(cells, counts)
    value = 0
    for pair in edges:
        B, C = tuple(pair)
        q = primitive(tuple(a-b for a, b in zip(profile(B), profile(C))))
        assert threshold(q) <= n
        value += threshold(q)
    assert len(cells) <= n*n and value <= 3*n*len(cells)
    counts['complete_actual_quadratic_Laurent_fans'] += 1
    return {'length': n, **metadata, 'order_ideal_vector_bases': len(candidates),
            'initial_ideals': len(cells), 'undirected_walls': len(edges),
            'total_primitive_cost': value, 'complete_area': '1', 'cells': ledger}


def full_unit(model, counts):
    algebra_controls(model, counts)
    candidates, controls, support = interval_bases(model, counts)
    out = extract_fan(model, candidates, counts, support)
    n, ev, _, p, _, metadata = model
    if n <= 7:
        exhaustive = set()
        for B in finite_staircases(n, 3):
            B = frozenset(B)
            if any(a[0] >= 2 or a[1] and a[2] for a in B):
                continue
            if rank([ev(a) for a in B], p) == n:
                exhaustive.add(B)
        assert exhaustive == set(candidates)
        counts['independent_exhaustive_interval_rank_equalities'] += 1
    if metadata['kind'] == 'sharp_exponential_Laurent' and p == 0:
        expected = n+sum(n-gcd(n, u) for u in range(1, n//2+1))
        bases = n+sum(u*(n-2*u+1) for u in range(1, n//2+1))
        assert out['initial_ideals'] == expected and len(candidates) == bases
        counts['independent_exact_sharp_fan_counts'] += 1
    out['normalized_rank_and_parity_controls'] = controls
    return out


def translated_ev(model, a, b):
    n, ev, M, p, f, metadata = model

    @lru_cache(None)
    def shifted(alpha):
        out = defaultdict(int)
        x, y, z = alpha
        for i in range(y+1):
            for j in range(z+1):
                coeff = comb(y, i)*a**(y-i)*comb(z, j)*b**(z-j)
                for k, c in ev((x, i, j)).items():
                    out[k] += coeff*c
        return normalized(out, p)

    matrices = [M[0]]
    for index, c in ((1, a), (2, b)):
        columns = []
        for j, col in enumerate(M[index]):
            out = defaultdict(int, col)
            out[j] += c
            columns.append(normalized(out, p))
        matrices.append(tuple(columns))
    return n, shifted, matrices, p, f, {**metadata, 'translation': [a, b]}


def affine_controls(model, original, counts):
    transformed = translated_ev(model, 2, -3)
    n, ev, _, p, _, _ = transformed
    for cell in original['cells']:
        B = frozenset(map(tuple, cell['standard_set']))
        assert rank([ev(a) for a in B], p) == n
        N, _ = normals(tuple(sorted(B)), ev, p)
        poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for q in N:
            poly = clip(poly, q, 1)
        assert set(poly) == {tuple(map(Q, v)) for v in cell['positive_cone']}
        counts['independent_affine_bilinear_full_cone_equalities'] += 1
    return {'length': n, 'translation': [2, -3], 'entire_fan_cones': len(original['cells'])}


def main():
    counts, rows, affine = Counter(), [], []
    models = [series_model(n) for n in range(1, 13)]
    models += [series_model(n, p) for n, p in ((4, 7), (5, 11), (6, 101), (7, 101))]
    models += [series_model(n, p, False) for n, p in ((4, 0), (6, 0), (6, 2), (7, 3))]
    models += [unit_points(((0, 1), (1, 2), (0, 3), (1, 5), (0, 7), (1, 11)), 0),
               unit_points(((0, 1), (1, 2), (0, 3), (1, 4), (0, 5)), 7),
               gaussian_units()]
    for model in models:
        row = full_unit(model, counts)
        rows.append(row)
        if ((row['length'] == 5 and row['field'] == 0)
                or row['kind'] == 'irreducible_quadratic_units_Q_i'):
            affine.append(affine_controls(model, row, counts))
        print('complete quadratic Laurent:', row['kind'], row['length'], row['field'],
              row['initial_ideals'], flush=True)
    # A larger rank-only fixture separates the cubic basis count from coherent count.
    large = series_model(20)
    algebra_controls(large, counts)
    candidates, _, _ = interval_bases(large, counts)
    assert len(candidates) == 20+sum(u*(21-2*u) for u in range(1, 11)) == 405 > 20**2
    rank_only = {'length': 20, 'order_ideal_vector_bases': 405,
                 'entire_spatial_fan_enumerated': False}
    # The linear-pair case is checked without assuming inverse translation.
    linear = linear_points()
    n, ev, _, p, _, _ = linear
    all_shapes = {basis_shape(n, 0, b, 0): (0, b, 0) for b in range(n)}
    for u in range(1, n//2+1):
        for delta in range(n-2*u+1):
            for d in range(u):
                all_shapes[basis_shape(n, u, delta, d)] = (u, delta, d)
    feasible = {B: param for B, param in all_shapes.items() if rank([ev(a) for a in B], p) == n}
    linear_row = extract_fan(linear, feasible, counts)
    assert linear_row['initial_ideals'] <= n+2
    sources = ['docs/quadratic_laurent_fans.md', 'scripts/run_quadratic_laurent_checks.py',
               'docs/quadric_cross_fans.md', 'scripts/run_quadric_cross_checks.py',
               'docs/nested_jet_matroids.md', 'docs/bilinear_rank_defects.md',
               'scripts/run_coordinate_slice_budget_checks.py',
               'scripts/run_signed_face_checks.py', 'scripts/run_structural_complexity_checks.py']
    out = {'schema_version': 1, 'all_checks_passed': True,
           'scope': 'Uniform all-field quadratic-anchor inverse-pair theorem and bilinear-defect '
                    'corollary; sharp characteristic-zero prefix family. Complete fans, parity paths '
                    'and rank translation are checked separately from the larger rank-only fixture. '
                    'No unrestricted G3 endpoint, polynomial separation or novelty established.',
           'uniform_theorem': {'unit_pair': 'g<=n^2, V<=3*n^3, |U|<=6*n^2',
                               'bilinear_defect_corollary': 'Monic quadratic in x plus rank(1,y,z,yz)<=3 '
                                                           'gives sharp quadratic fan and cubic cost order.',
                               'sharp_count': 'n+sum_(u=1)^floor(n/2)(n-gcd(n,u))',
                               'coherence_essential': 'The number of order-ideal vector bases can be cubic.'},
           'priority_limits': 'Monotone lattice optimization is classical; Topkis (1978), '
                              'doi:10.1287/opre.26.2.305, publisher abstract inspected. '
                              'The proof here is self-contained. Exact numerical class-bound '
                              'and sharp-count priority remains unassessed.',
           'counters': dict(counts), 'complete_actual_unit_fans': rows,
           'affine_translation_controls': affine, 'linear_pair_control': linear_row,
           'larger_rank_only_fixture': rank_only,
           'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}}
    path = ROOT/'results/research/quadratic_laurent_checks.json'
    path.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'all_checks_passed': True, 'counters': dict(counts),
                      'report': str(path.relative_to(ROOT))}), flush=True)


if __name__ == '__main__':
    main()
