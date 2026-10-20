#!/usr/bin/env python3
"""Exact controls for profile lattices and the refined shadow budgets."""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import gcd, prod
from pathlib import Path
from random import Random
import hashlib
import json

from run_nested_jet_checks import rank_columns
from run_planar_dyadic_bound_checks import convex_hull
from run_rational_shadow_bound_checks import cover, det2, direction_control, project, sweep
from run_rectangle_rank_checks import rank, rectangle
from run_structural_complexity_checks import cross, determinant, dot, primitive

ROOT = Path(__file__).resolve().parents[1]


def difference(a, b):
    return tuple(x-y for x, y in zip(a, b))


def profile(ground, basis):
    return tuple(sum(ground[j][i] for j in basis) for i in range(len(ground[0])))


def pair_minors(rows):
    return [det2(a, b) for a, b in combinations(zip(*rows), 2)]


def oracle_polygon(ground, independent):
    """Classical planar support-oracle enumeration; no full basis list used."""
    queries = []
    def query(weight):
        secondary = (-weight[1], weight[0])
        chosen = set()
        for j in sorted(range(len(ground)), key=lambda j: (dot(weight, ground[j]),
                                                           dot(secondary, ground[j]), j)):
            if independent(frozenset(chosen | {j})):
                chosen.add(j)
        point = profile(ground, chosen)
        queries.append({'primary_weight': weight, 'secondary_weight': secondary,
                        'basis': sorted(chosen), 'profile': point})
        return point
    left, right = query((1, 0)), query((-1, 0))
    assert left != right
    vertices = {left, right}
    pending = [(left, right), (right, left)]
    edges = []
    while pending:
        a, b = pending.pop()
        delta = difference(b, a)
        c = query((delta[1], -delta[0]))
        height = det2(delta, difference(c, a))
        assert height >= 0
        if height:
            assert c not in vertices
            vertices.add(c)
            pending.extend(((a, c), (c, b)))
        else:
            edges.append((a, b))
    successor = dict(edges)
    assert len(successor) == len(vertices)
    polygon = [left]
    while successor[polygon[-1]] != left:
        polygon.append(successor[polygon[-1]])
        assert len(polygon) <= len(vertices)
    assert set(polygon) == vertices and len(queries) == 2*len(vertices)
    return {'vertices': polygon, 'support_queries': queries}


def normalization_controls():
    rng = Random(82109)
    records = []
    for d in range(2, 6):
        for _ in range(25):
            while True:
                u, v = [tuple(rng.randrange(-5, 6) for _ in range(d)) for j in range(2)]
                if rank_columns([u, v]) == 2:
                    break
            plucker = pair_minors((u, v))
            gram = dot(u, u)*dot(v, v)-dot(u, v)**2
            assert gram == sum(x*x for x in plucker)
            # These three generators span precisely Z*u+Z*v: their coefficient
            # minors are 6, 2 and -3, whose gcd is one.
            generators = [tuple(2*x for x in u), tuple(3*x for x in v),
                          tuple(x+y for x, y in zip(u, v))]
            base_eta = Fraction(sum(map(abs, plucker)), gram)
            controls = []
            for transform in (((1, 0), (0, 1)), ((2, 0), (0, 3)), ((1, 2), (3, 7))):
                rows = [tuple(a*x+b*y for x, y in zip(u, v)) for a, b in transform]
                images = [project(rows, g) for g in generators]
                c = gcd(*(det2(a, b) for a, b in combinations(images, 2)))
                assert c == abs(det2(*transform))*gram
                eta = Fraction(sum(map(abs, pair_minors(rows))), c)
                assert eta == base_eta <= 1
                controls.append({'row_transform': transform, 'exchange_covolume': c, 'eta': str(eta)})
            records.append({'dimension': d, 'basis_rows': [u, v], 'plucker_coordinates': plucker,
                            'gram_determinant': gram, 'row_basis_controls': controls})
    return records


def rank_two_control(label, ground, columns, prime=0, bases=None, face=False):
    ground = tuple(ground)
    r = rank(columns, prime)
    @lru_cache(None)
    def independent(ids):
        return (len(ids) <= r and rank([columns[j] for j in sorted(ids)], prime) == len(ids)
                and (not face or any(ids <= b for b in bases)))
    if bases is None:
        assert len(ground) <= 12
        bases = [frozenset(b) for b in combinations(range(len(ground)), r) if independent(frozenset(b))]
    assert bases and all(len(b) == r and independent(b) for b in bases)
    basis_set = set(bases)
    profiles = sorted({profile(ground, b) for b in bases})
    generators = sorted({difference(p, profiles[0]) for p in profiles} - {(0,)*len(ground[0])})
    rows = next((u, v) for u, v in combinations(generators, 2) if rank_columns([u, v]) == 2)
    assert all(rank_columns([*rows, u]) == 2 for u in generators)
    projected_generators = [project(rows, u) for u in generators]
    c = gcd(*(det2(u, v) for u, v in combinations(projected_generators, 2)))
    assert c > 0
    minors = pair_minors(rows)
    eta = Fraction(sum(map(abs, minors)), c)
    assert 0 < eta <= 1
    exchanges = set()
    for b in bases:
        for e in b:
            for f in set(range(len(ground)))-b:
                if b-{e}|{f} in basis_set:
                    exchanges.add(difference(ground[f], ground[e]))
    assert all(rank_columns([*rows, u]) == 2 for u in exchanges)
    exchange_c = gcd(*(det2(project(rows, u), project(rows, v)) for u, v in combinations(exchanges, 2)))
    assert exchange_c == c
    H = max(prod(x+1 for x in u) for u in ground)
    boxes = cover(rows, H)
    records = []
    direction_checks = 0
    for signs in product((-1, 1), repeat=2):
        record = sweep(ground, rows, r, independent, H, signs)
        local = defaultdict(lambda: defaultdict(lambda: [0, 0]))
        pivots_checked = 0
        for wall in record['walls']:
            for pivot in wall['pivots']:
                e, f = pivot['from'], pivot['to']
                delta = difference(ground[f], ground[e])
                assert delta in exchanges
                assert rank_columns([*rows, delta]) == 2
                k = pivot['label_from']
                if k == pivot['label_to']:
                    vector = project(boxes[k]['T'], project(rows, delta))
                    for i in range(2):
                        local[k][Fraction(wall['parameter'])][i] += vector[i]
                pivots_checked += 1
        finer = []
        for k, events in sorted(local.items()):
            vectors = [v for _, v in sorted(events.items())]
            control = direction_control(vectors, c*boxes[k]['minor'])
            assert all(v <= 2*r*w for v, w in zip(control['absolute_budgets'], boxes[k]['widths']))
            finer.append({'label': k, 'scaled_determinant_over_exchange_covolume': str(Fraction(boxes[k]['D'], c)),
                          **control})
            direction_checks += 1
        record['exchange_lattice_covolume'] = c
        record['exchange_normalized_height'] = str(eta)
        record['every_successful_swap_checked'] = pivots_checked
        record['finer_internal_controls'] = finer
        records.append(record)
    exposed = {tuple(s['profile']) for rec in records for s in rec['states']}
    hull = convex_hull({project(rows, p) for p in profiles})
    assert exposed == set(hull)
    oracle = oracle_polygon([project(rows, u) for u in ground], independent)
    assert set(oracle['vertices']) == set(hull)
    # The projected vertex list injects back to the original profile plane.
    assert len({project(rows, p) for p in profiles}) == len(profiles)
    return {'model': label, 'characteristic': prime, 'rank': r, 'H': H,
            'ground': ground, 'full_bases': [sorted(b) for b in bases], 'full_profiles': profiles,
            'projection_rows': rows, 'projection_minors': minors,
            'exchange_covolume': c, 'exchange_height': str(eta),
            'full_exchange_directions': sorted(exchanges), 'full_profile_hull': hull,
            'support_oracle_enumeration': oracle,
            'boxes_with_fractional_D_over_c': sum(box['D'] < c for box in boxes),
            'complete_quadrant_sweeps': records, 'finer_direction_checks': direction_checks}


def matroid_controls():
    fixtures = []
    for M in (2, 7, 17):
        u, v = (1, M, 0), (0, 1, M)
        ground = [(0, 0, 0), u, v, tuple(a+b for a, b in zip(u, v))]
        fixtures.append(rank_two_control(f'rank_one_coarse_plane_M{M}', ground, [(1,)]*4))
    for d, u, v in ((3, (1, -1, 0), (0, 1, -1)),
                    (3, (2, 0, -2), (0, 3, -3)),
                    (4, (1, 2, 0, 1), (0, 1, 3, 1))):
        offsets = [(6,)*d, tuple(7+i for i in range(d))]
        ground = []
        for offset, directions in zip(offsets, ((u, v), (tuple(2*x for x in u), tuple(2*x for x in v)))):
            ground.extend([offset, tuple(a+b for a, b in zip(offset, directions[0])),
                           tuple(a+b for a, b in zip(offset, directions[1]))])
        assert all(min(u) >= 0 for u in ground)
        assert rank_columns([difference(a, ground[0]) for a in ground]) > 2
        columns = [(1, 0)]*3+[(0, 1)]*3
        fixtures.append(rank_two_control(f'partition_offset_planes_d{d}_u{u}', ground, columns))
    return fixtures


def algebra_controls():
    fixtures = []
    for p in (0, 2, 3):
        data = rectangle(2, 2, {(1, 0): 1, (0, 1): 1}, p)
        ground = [u for u in data['ground'] if any(data['columns'][u])]
        fixtures.append(rank_two_control(f'linear_rectangle_char{p}', ground,
                                        [data['columns'][u] for u in ground], p))
    # A homogeneous algebra of embedding dimension three: all cubics vanish,
    # all squares vanish, and xy=xz+yz.
    reference = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1))
    ground = [*reference, (1, 1, 0)]
    for p in (0, 2, 3):
        columns = [tuple(int(i == j) for i in range(6)) for j in range(6)]
        columns.append((0, 0, 0, 0, 1, 1))
        column_map = dict(zip(ground, columns))
        zero = (0,)*6
        operators = []
        for k in range(3):
            operators.append([column_map.get(tuple(a+int(i == k) for i, a in enumerate(u)), zero)
                              for u in reference])
        def apply(matrix, column):
            result = tuple(sum(v[i]*c for v, c in zip(matrix, column)) for i in range(6))
            return tuple(x % p for x in result) if p else result
        for left, right in combinations(operators, 2):
            for j in range(6):
                assert apply(left, right[j]) == apply(right, left[j])
        for u, expected in zip(ground, columns):
            replay = columns[0]
            for op, k in zip(operators, u):
                for _ in range(k):
                    replay = apply(op, replay)
            assert replay == expected
        for op in operators:
            assert all(apply(op, column) == zero for column in op)
        for j in range(6):
            xy = apply(operators[0], operators[1][j])
            xz = apply(operators[0], operators[2][j])
            yz = apply(operators[1], operators[2][j])
            assert all((a-b-c) % p == 0 if p else a == b+c for a, b, c in zip(xy, xz, yz))
        record = rank_two_control(f'graded_embedding_dimension_three_char{p}', ground, columns, p)
        record['commuting_multiplication_columns'] = operators
        record['commutator_column_checks'] = 18
        fixtures.append(record)
    for weights, cutoff in (((3, 4, 5), 11), ((3, 5, 7), 16)):
        ground = sorted(u for u in product(*(range((cutoff-1)//w+1) for w in weights)) if dot(u, weights) < cutoff)
        values = sorted({dot(u, weights) for u in ground})
        columns = [tuple(int(value == dot(u, weights)) for value in values) for u in ground]
        groups = [[j for j, u in enumerate(ground) if dot(u, weights) == value] for value in values]
        bases = [frozenset(b) for b in product(*groups)]
        assert len(bases) < 20000
        fixtures.append(rank_two_control(f'truncated_monomial_curve_{weights}_L{cutoff}', ground, columns, bases=bases))
    return fixtures


def facet_controls():
    records = []
    for prime in (0, 2, 3):
        data = rectangle(2, 2, {(1, 0): 1, (0, 1): 1, (1, 1): 1}, prime)
        ground = [u for u in data['ground'] if any(data['columns'][u])]
        columns = [data['columns'][u] for u in ground]
        bases = [frozenset(b) for b in combinations(range(len(ground)), 4)
                 if rank([columns[j] for j in b], prime) == 4]
        profiles = sorted({profile(ground, b) for b in bases})
        assert rank_columns([difference(p, profiles[0]) for p in profiles]) == 3
        planes = set()
        for a, b, c in combinations(profiles, 3):
            normal = cross(difference(b, a), difference(c, a))
            divisor = gcd(*normal)
            if not divisor:
                continue
            normal = tuple(x//divisor for x in normal)
            level = dot(normal, a)
            values = [dot(normal, p) for p in profiles]
            if min(values) < level < max(values):
                continue
            if level == max(values):
                normal, level = tuple(-x for x in normal), -level
            planes.add((normal, level))
        for normal, level in sorted(planes):
            face_bases = [b for b in bases if dot(normal, profile(ground, b)) == level]
            for b in face_bases:
                for c in face_bases:
                    for e in b-c:
                        assert any(b-{e}|{f} in face_bases for f in c-b)
            control = rank_two_control(f'full_rank_rectangle_facet_char{prime}_normal{normal}',
                                       ground, columns, prime, face_bases, face=True)
            control['supporting_plane'] = {'inward_normal': normal, 'level': level}
            records.append(control)
    return records


def enlarged_ground_controls():
    answer = []
    for n, p in ((2, 0), (2, 2), (3, 3), (5, 0)):
        # K[x,y,z]/(x^n-1,y,z). The extra head x^n is essential.
        reference = [tuple(int(i == j) for i in range(n)) for j in range(n)]
        columns = reference+[reference[0]]
        bases = [b for b in combinations(range(n+1), n) if rank([columns[j] for j in b], p) == n]
        profiles = sorted({sum(b) for b in bases})
        assert len(bases) == 2 and profiles[1]-profiles[0] == n
        assert prod((n+1, 1, 1)) > n and prod((n+1, 1, 1)) <= 2*n
        answer.append({'length': n, 'characteristic': p, 'head': [n, 0, 0],
                       'old_ground_profile_rank': 0, 'enlarged_ground_profile_rank': 1,
                       'profile_difference_generator': [n, 0, 0], 'all_profiles': profiles})
    return answer


def rectangle_index_controls():
    records = []
    fixtures = [(2, (2,), 0), (3, (3,), 0), (3, (4,), 0), (4, (4,), 0),
                (5, (5,), 0), (4, (3, 5), 0), (4, (2, 5), 0),
                (3, (3,), 2), (3, (4,), 3), (5, (5,), 2), (5, (5,), 3)]
    for q, degrees, p in fixtures:
        F = {(a, b): 1 for a, b in product(range(q), repeat=2) if a+b in degrees}
        assert F
        F.update({(1, 0): 1, (0, 1): 1})
        data = rectangle(q, q, F, p)
        nonlinear = [(a, b) for a, b in data['F'] if a+b >= 2]
        h = gcd(*(a+b-1 for a, b in nonlinear))
        support = [(1, 0), (0, 1), *nonlinear]
        exchange_generators = [(-a, -b, 1) for a, b in support]
        determinants = [determinant([exchange_generators[0], exchange_generators[1], w])
                        for w in exchange_generators[2:]]
        assert all(v.denominator == 1 for v in determinants)
        index_witnesses = [abs(int(v)) for v in determinants]
        assert gcd(*index_witnesses) == h
        assert rank_columns(exchange_generators) == 3
        fundamental = 0
        for u in data['ground']:
            for cell, value in zip(data['cells'], data['columns'][u]):
                if value:
                    assert (sum(u)-sum(cell)) % h == 0
                    fundamental += 1
        reference = [tuple(int(i == j) for i in range(q*q)) for j in range(q*q)]
        explicit_swaps = []
        for cell in support:
            j = data['index'][cell]
            chosen = reference[:]
            chosen[j] = data['columns'][(0, 0, 1)]
            assert rank(chosen, p) == q*q
            explicit_swaps.append({'tail': cell, 'profile_difference': (-cell[0], -cell[1], 1)})
        D, S = min(degrees), 2*q-2
        component = None
        if D >= q:
            rows = [i for i, cell in enumerate(data['cells']) if sum(cell) <= S-D+1 or sum(cell) >= D]
            assert len(rows) == (2*q-D)**2
            assert all(data['index'][cell] in rows for cell in support)
            for cell in support:
                j = rows.index(data['index'][cell])
                chosen = [tuple(int(i == k) for i in range(len(rows))) for k in range(len(rows))]
                chosen[j] = tuple(data['columns'][(0, 0, 1)][i] for i in rows)
                assert rank(chosen, p) == len(rows)
            component = {'rank': len(rows), 'same_exact_index': h,
                         'fan_bound_scale_squared': str(Fraction((len(rows)*S)**3, h))}
        n = q*q
        volume = Fraction(n**3*(q-1)**2*S, 8)
        records.append({'side': q, 'length': n, 'characteristic': p, 'nonlinear_degrees': degrees,
                        'intrinsic_lattice_index': h, 'nonzero_fundamental_exchanges': fundamental,
                        'exact_index_determinants': index_witnesses, 'verified_swaps': explicit_swaps,
                        'moment_box_volume_over_index': str(volume/h), 'high_order_component': component})
    return records


def main():
    normalization = normalization_controls()
    abstract = matroid_controls()
    algebras = algebra_controls()
    facets = facet_controls()
    enlarged = enlarged_ground_controls()
    rectangles = rectangle_index_controls()
    all_rows = abstract+algebras+facets
    summary = {'complete_rank_two_profile_models': len(all_rows),
               'independent_lattice_normalizations': sum(len(r['row_basis_controls']) for r in normalization),
               'full_bases_enumerated': sum(len(f['full_bases']) for f in all_rows),
               'full_basis_exchange_directions': sum(len(f['full_exchange_directions']) for f in all_rows),
               'complete_quadrant_sweeps': sum(len(f['complete_quadrant_sweeps']) for f in all_rows),
               'replayed_successful_swaps': sum(r['every_successful_swap_checked'] for f in all_rows for r in f['complete_quadrant_sweeps']),
               'finer_internal_direction_budgets': sum(f['finer_direction_checks'] for f in all_rows),
               'fractional_box_covolume_controls': sum(f['boxes_with_fractional_D_over_c'] for f in all_rows),
               'enlarged_ground_fixtures': len(enlarged), 'rectangle_index_fixtures': len(rectangles),
               'full_rank_quotient_facets': len(facets),
               'rectangle_fundamental_exchanges': sum(r['nonzero_fundamental_exchanges'] for r in rectangles),
               'high_order_exact_component_indices': sum(r['high_order_component'] is not None for r in rectangles)}
    summary['support_oracle_calls'] = sum(len(r['support_oracle_enumeration']['support_queries']) for r in all_rows)
    summary['independently_enumerated_profile_vertices'] = sum(len(r['support_oracle_enumeration']['vertices']) for r in all_rows)
    assert summary['fractional_box_covolume_controls'] > 0
    sources = ['scripts/run_profile_lattice_checks.py', 'docs/profile_lattice_bound.md',
               'scripts/run_rational_shadow_bound_checks.py', 'scripts/run_rectangle_rank_checks.py',
               'scripts/run_nested_jet_checks.py', 'scripts/run_planar_dyadic_bound_checks.py',
               'scripts/run_structural_complexity_checks.py']
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_profile_lattice_checks.py',
              'scope': 'Uniform exchange-lattice proof gives O_d(n log(n+1)) for rank-at-most-two intrinsic profile lattices. Finite grading index sharpens rectangle volume bounds. No unrestricted G_3 exponent improvement or completed novelty assessment is claimed.',
              'summary': summary, 'lattice_normalization_controls': normalization,
              'abstract_profile_models': abstract, 'graded_algebra_models': algebras,
              'full_rank_quotient_facet_models': facets,
              'enlarged_ground_controls': enlarged, 'rectangle_index_controls': rectangles,
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              'unrestricted_endpoint_achieved': False, 'novelty_established': False, 'goal_complete': False}
    (ROOT/'results/research/profile_lattice_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', **summary}, indent=2))


if __name__ == '__main__':
    main()
