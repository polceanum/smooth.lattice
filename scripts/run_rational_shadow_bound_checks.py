#!/usr/bin/env python3
"""Exact controls for rational shadows and short-normal vertex visibility.

The uniform proofs are in docs/rational_shadow_bound.md. These finite
controls do not establish asymptotic sharpness, novelty, or goal completion.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import gcd, isqrt, prod
from pathlib import Path
from random import Random
import hashlib
import json

from run_extremal_initial_ideals_checks import rectangle_model
from run_planar_dyadic_bound_checks import convex_hull, matrix_oracle
from run_structural_complexity_checks import cross, determinant, dot, hyperbolic_cross, primitive

ROOT = Path(__file__).resolve().parents[1]


def det2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def project(A, a):
    return tuple(dot(row, a) for row in A)


def projection_data(A):
    columns = tuple(zip(*A))
    minors = [det2(a, b) for a, b in combinations(columns, 2)]
    delta = gcd(*minors)
    assert delta > 0
    return delta, sum(abs(x) for x in minors) // delta


def kernel_projection(q):
    q = primitive(q)
    pivot = next(i for i, x in enumerate(q) if x)
    rows = []
    for j in range(3):
        if j != pivot:
            row = [0] * 3
            row[pivot], row[j] = q[j], -q[pivot]
            rows.append(tuple(row))
    assert all(dot(row, q) == 0 for row in rows)
    assert primitive(cross(*rows)) == q
    assert projection_data(rows)[1] == sum(abs(x) for x in q)
    return tuple(rows)


def cover(A, H):
    d = len(A[0])
    L = H.bit_length() - 1
    delta, height = projection_data(A)
    result = []
    for indices in product(range(L + 1), repeat=d - 1):
        if sum(indices) > L:
            continue
        sides = tuple(2 ** (j + 1) - 1 for j in indices) + (max(1, H // 2 ** sum(indices) - 1),)
        assert min(sides) >= 1 and prod(b + 1 for b in sides) <= 2 ** d * H
        vectors = [tuple(b * x for x in col) for b, col in zip(sides, zip(*A))]
        i, j = max(combinations(range(d), 2), key=lambda pair: abs(det2(*(vectors[k] for k in pair))))
        if det2(vectors[i], vectors[j]) < 0:
            i, j = j, i
        a, b = tuple(zip(*A))[i], tuple(zip(*A))[j]
        minor = det2(a, b)
        D = sides[i] * sides[j] * minor
        assert D > 0 and all(abs(det2(v, b)) <= sides[i] * minor
                             and abs(det2(a, v)) <= sides[j] * minor for v in vectors)
        T = ((b[1], -b[0]), (-a[1], a[0]))
        widths = (2 * d * sides[i] * minor + 1, 2 * d * sides[j] * minor + 1)
        assert det2(*T) == minor
        assert prod(widths) <= (2 * d + 1) ** 2 * D * minor
        assert D <= 2 ** d * H * delta * height
        result.append({'sides': sides, 'T': T, 'D': D, 'minor': minor, 'widths': widths})
    return result


def shared_lines(boxes):
    result = set()
    for box in boxes:
        for row, width in zip(box['T'], box['widths']):
            g = gcd(*row)
            sign = 1 if next(x for x in row if x) > 0 else -1
            normal = tuple(sign * x // g for x in row)
            for s in (-1, 1):
                result.add((normal, Fraction(s * width, 2 * g)))
    return sorted(result)


def direction_control(vectors, covolume):
    """Vectors are ordered rays in a pointed cone, including oblique lattices."""
    k = len(vectors)
    U, V = (sum(abs(v[i]) for v in vectors) for i in range(2))
    if k <= 1:
        return {'directions': k, 'absolute_budgets': [U, V], 'short_vectors': vectors}
    assert U > 0 and V > 0
    signs = [det2(a, b) for a, b in zip(vectors, vectors[1:])]
    assert all(x > 0 for x in signs) or all(x < 0 for x in signs)
    small = [v for v in vectors if Fraction(abs(v[0]), U) + Fraction(abs(v[1]), V) <= Fraction(3, k)]
    assert 3 * len(small) >= k
    twice_area = sum(abs(det2(a, b)) for a, b in zip(small, small[1:]))
    assert all(abs(det2(a, b)) >= covolume and det2(a, b) % covolume == 0
               for a, b in zip(small, small[1:]))
    assert twice_area >= (len(small) - 1) * covolume
    assert twice_area * k * k <= 72 * U * V, 'triangles fit inside the centered short-vector rectangle'
    if k >= 6:
        assert k ** 3 * covolume <= 432 * U * V
    return {'directions': k, 'absolute_budgets': [U, V], 'short_vectors': small,
            'short_triangle_twice_area': twice_area, 'lattice_covolume': covolume}


def sweep(original, A, rank, independent, H, signs):
    assert all(prod(x + 1 for x in a) <= H and min(a) >= 0 for a in original)
    ground = tuple(project(A, a) for a in original)
    n = len(ground)
    boxes = cover(A, H)
    lines = shared_lines(boxes)
    d = len(A[0])
    assert len(lines) <= 2 * d * (d - 1) * H.bit_length()
    delta, height = projection_data(A)
    labels = [next(i for i, b in enumerate(boxes)
                   if all(2 * abs(x) < width for x, width in zip(project(b['T'], v), b['widths']))) for v in ground]
    sigma, tau = signs

    def greedy(t):
        selected = set()
        for j in sorted(range(n), key=lambda j: (sigma * ground[j][0] * t + tau * ground[j][1], j)):
            if independent(frozenset(selected | {j})):
                selected.add(j)
        assert len(selected) == rank
        return selected

    events = defaultdict(set)
    for i, j in combinations(range(n), 2):
        dx, dy = (a - b for a, b in zip(ground[i], ground[j]))
        if sigma * tau * dx * dy < 0:
            events[Fraction(-tau * dy, sigma * dx)].update((i, j))
    events = sorted(events.items())
    first = events[0][0] / 2 if events else Fraction(1)
    basis = greedy(first)
    order = sorted(range(n), key=lambda j: (tau * ground[j][1], sigma * ground[j][0], j))
    positions = {j: i for i, j in enumerate(order)}
    paths = [[ground[j]] for j in sorted(basis)]
    token = {j: i for i, j in enumerate(sorted(basis))}
    local = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    variation = defaultdict(lambda: [0, 0])
    cross_moves = 0
    states, walls = [], []

    def record(t):
        assert basis == greedy(t), 'independent direct greedy replay on every interval'
        if not states or states[-1]['basis'] != sorted(basis):
            states.append({'basis': sorted(basis), 'profile': [sum(ground[j][i] for j in basis) for i in range(2)],
                           'parameter': str(t)})

    record(first)
    for event_index, (t, involved) in enumerate(events):
        groups = defaultdict(list)
        for j in involved:
            groups[sigma * ground[j][0] * t + tau * ground[j][1]].append(j)
        pivots = []
        for group in groups.values():
            locations = sorted(positions[j] for j in group)
            assert locations == list(range(locations[0], locations[-1] + 1))
            for target, j in enumerate(sorted(group, key=lambda j: (sigma * ground[j][0], j)), locations[0]):
                while positions[j] > target:
                    pos = positions[j]
                    e, f = order[pos - 1], j
                    if e in basis and f not in basis and independent(frozenset(basis - {e} | {f})):
                        step = tuple(b - a for a, b in zip(ground[e], ground[f]))
                        assert -sigma * step[0] > 0 and tau * step[1] > 0
                        assert Fraction(tau * step[1], -sigma * step[0]) == t
                        k = token.pop(e)
                        token[f] = k
                        paths[k].append(ground[f])
                        before, after = labels[e], labels[f]
                        if before != after:
                            cross_moves += 1
                            assert any((dot(row, ground[e]) < level) != (dot(row, ground[f]) < level)
                                       for row, level in lines)
                        else:
                            image = project(boxes[before]['T'], step)
                            for i in range(2):
                                local[before][t][i] += image[i]
                                variation[k, before][i] += abs(image[i])
                        basis.remove(e)
                        basis.add(f)
                        pivots.append({'token': k, 'from': e, 'to': f, 'label_from': before, 'label_to': after})
                    order[pos - 1], order[pos] = f, e
                    positions[f], positions[e] = pos - 1, pos
        if pivots:
            walls.append({'parameter': str(t), 'pivots': pivots})
        following = events[event_index + 1][0] if event_index + 1 < len(events) else t + 2
        record((t + following) / 2)
    assert len(states) == len(walls) + 1
    boundary_crossings = 0
    for path in paths:
        for row, level in lines:
            values = [dot(row, v) for v in path]
            assert all(x != level for x in values)
            changes = [b - a for a, b in zip(values, values[1:]) if b != a]
            assert sum(a * b < 0 for a, b in zip(changes, changes[1:])) <= 1, 'two monotone phases'
            side = [v < level for v in values]
            count = sum(a != b for a, b in zip(side, side[1:]))
            assert count <= 2, 'a convex token path crosses a ground-avoiding line at most twice'
            boundary_crossings += count
    assert cross_moves <= boundary_crossings <= 2 * rank * len(lines)
    controls = []
    for label, rows in sorted(local.items()):
        box = boxes[label]
        for k in range(rank):
            assert all(x <= 2 * w for x, w in zip(variation[k, label], box['widths'])), 'returning to a box does not restart its budget'
        vectors = [v for _, v in sorted(rows.items())]
        control = direction_control(vectors, delta * box['minor'])
        assert all(x <= 2 * rank * w for x, w in zip(control['absolute_budgets'], box['widths']))
        control.update({'label': label, 'box': box})
        controls.append(control)
    assert len(walls) <= cross_moves + sum(len(rows) for rows in local.values())
    return {'signs': signs, 'states': states, 'walls': walls, 'token_paths': paths,
            'projection': A, 'lattice_covolume': delta, 'normalized_height': height,
            'cover_boxes': len(boxes), 'distinct_boundary_lines': len(lines), 'cross_label_moves': cross_moves,
            'boundary_crossings': boundary_crossings, 'internal_controls': controls,
            'all_intervals_replayed': True}


def matroid_controls():
    q_values = ((0, 0, 1), (1, 1, 1), (2, 3, 5), (1, -2, 3), (7, 11, 13))
    models = []
    for H, rank in ((1, 0), (1, 1), (3, 2), (4, 4), (6, 6), (8, 5)):
        ground = tuple(sorted(hyperbolic_cross(H, 3)))
        models.append((f'uniform_H{H}_r{rank}', H, ground, rank, lambda s, r=rank: len(s) <= r))
        random = Random(100 * H + rank)
        columns = [tuple(int(i == j) for i in range(rank)) if j < rank
                   else tuple(random.randrange(101) for _ in range(rank)) for j in range(len(ground))]
        models.append((f'linear_H{H}_r{rank}', H, ground, rank, matrix_oracle(columns, rank, 101)))
    for r, s, p in ((2, 2, 2), (2, 2, 101), (3, 2, 101)):
        coefficients = {(1, 0): 1, (0, 1): 1, (1, 1): 3}
        ground, columns, _ = rectangle_model(r, s, coefficients, p)
        models.append((f'local_rectangle_{r}_{s}_char{p}', r * s, ground, r * s,
                       matrix_oracle([columns[a] for a in ground], r * s, p)))
    ground = tuple((8 - i, i * i, 0) for i in range(9))
    models.append(('rank_one_convex_chain', 128, ground, 1, lambda s: len(s) <= 1))
    ground = ((0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0),
              (0, 1, 0), (1, 1, 0), (0, 0, 1), (0, 1, 1))
    pairs = [frozenset((2 * i, 2 * i + 1)) for i in range(4)]
    forbidden = {a | b for a, b in combinations(pairs, 2) if a | b != pairs[2] | pairs[3]}
    bases = {frozenset(s) for s in combinations(range(8), 4)} - forbidden
    assert all(any(B - {e} | {f} in bases for f in C - B) for B in bases for C in bases for e in B - C)
    models.append(('sparse_paving_basis_exchange_verified', 4, ground, 4,
                   lambda s: len(s) < 4 or len(s) == 4 and s not in forbidden))
    fixtures, ledger = [], []
    all_sweeps = interval_states = hull_controls = 0
    for label, H, ground, rank, oracle in models:
        for q in q_values:
            A = kernel_projection(q)
            rows = [sweep(ground, A, rank, oracle, H, signs) for signs in product((-1, 1), repeat=2)]
            exposed = {tuple(state['profile']) for row in rows for state in row['states']}
            if len(ground) <= 12:
                profiles = {tuple(sum(project(A, ground[j])[i] for j in b) for i in range(2))
                            for b in combinations(range(len(ground)), rank) if oracle(frozenset(b))}
                assert exposed == set(convex_hull(profiles)), 'independent enumeration of every basis-profile hull'
                hull_controls += 1
            all_sweeps += 4
            interval_states += sum(len(row['states']) for row in rows)
            ledger.append([label, q, [len(row['states']) for row in rows],
                           [row['cross_label_moves'] for row in rows],
                           [[c['directions'] for c in row['internal_controls']] for row in rows]])
            # Save representative complete witnesses and all nontrivial internal budgets.
            if q == (2, 3, 5) or label == 'rank_one_convex_chain':
                fixtures.append({'model': label, 'H': H, 'rank': rank, 'ground': ground, 'quadrants': rows})
    return {'models': len(models), 'projection_kernel_vectors': q_values,
            'quadrant_sweeps': all_sweeps, 'interval_states': interval_states,
            'independent_full_basis_hulls': hull_controls,
            'ledger_sha256': hashlib.sha256(json.dumps(ledger, sort_keys=True).encode()).hexdigest(),
            'exact_fixtures': fixtures}


def short_relation(w, Q):
    if w[2] == 0:
        return (0, 0, 1)
    for a, b in product(range(-Q, Q + 1), repeat=2):
        numerator = -w[0] * a - w[1] * b
        if numerator % w[2] == 0:
            c = numerator // w[2]
            if abs(c) <= Q and (a or b or c):
                return primitive((a, b, c))
    raise AssertionError('pigeonhole short orthogonal relation must exist')


def ceiling_sqrt(n):
    return isqrt(n) + int(isqrt(n) ** 2 != n)


def short_relation_controls():
    ledger = []
    for w in product(range(-5, 6), repeat=3):
        if not any(w):
            continue
        W = sum(abs(x) for x in w)
        Q = ceiling_sqrt(W)
        assert (Q + 1) ** 3 > Q * W + 1
        q = short_relation(w, Q)
        assert dot(w, q) == 0 and 0 < max(abs(x) for x in q) <= Q and gcd(*q) == 1
        assert projection_data(kernel_projection(q))[1] == sum(abs(x) for x in q)
        ledger.append([w, q, Q])
    return {'integer_vectors_checked': len(ledger),
            'ledger_sha256': hashlib.sha256(json.dumps(ledger).encode()).hexdigest(),
            'examples': ledger[::max(1, len(ledger) // 12)]}


def geometry_controls():
    cube = tuple(product((0, 1), repeat=3))
    fixtures = [('cube', cube),
                ('skew_cube', [project(((1, 3, 0), (0, 1, 2), (0, 0, 1)), a) for a in cube]),
                ('octahedron', [(s * int(i == j) for i in range(3)) for j in range(3) for s in (-1, 1)]),
                ('non_simple_pyramid', [(0, 0, 0), (3, 0, 0), (3, 3, 0), (0, 3, 0), (1, 1, 4)]),
                ('integer_cloud', [(0, 0, 0), (4, 0, 0), (0, 4, 0), (0, 0, 4), (3, 2, 1), (2, 3, 3), (1, 2, 4)])]
    result, visible = [], 0
    for name, supplied in fixtures:
        points = sorted({tuple(a) for a in supplied})
        planes = {}
        for a, b, c in combinations(points, 3):
            normal = cross(tuple(x - y for x, y in zip(b, a)), tuple(x - y for x, y in zip(c, a)))
            if not any(normal):
                continue
            normal = primitive(normal)
            level = dot(normal, a)
            sides = [dot(normal, p) - level for p in points]
            if min(sides) < 0 < max(sides):
                continue
            if max(sides) > 0:
                normal, level = tuple(-x for x in normal), -level
            planes[normal + (level,)] = [p for p in points if dot(normal, p) == level]
        facets, vertices = [], set()
        for plane, face in planes.items():
            normal, level = plane[:3], plane[3]
            drop = max(range(3), key=lambda i: abs(normal[i]))
            axes = [i for i in range(3) if i != drop]
            image = {tuple(p[i] for i in axes): p for p in face}
            polygon = [image[p] for p in convex_hull(image)]
            vertices.update(polygon)
            multiples = []
            for i in range(1, len(polygon) - 1):
                area_vector = cross(tuple(a - b for a, b in zip(polygon[i], polygon[0])),
                                    tuple(a - b for a, b in zip(polygon[i + 1], polygon[0])))
                ratio = Fraction(area_vector[drop], normal[drop])
                assert ratio.denominator == 1 and ratio != 0
                assert area_vector == tuple(ratio * x for x in normal)
                multiples.append(abs(int(ratio)))
            assert len(polygon) <= 3 * sum(multiples), 'facet vertices paid by normalized lattice area'
            facets.append({'normal': normal, 'level': level, 'vertices': polygon, 'twice_normalized_area': sum(multiples)})
        certificates = []
        for v in sorted(vertices):
            incident = [f['normal'] for f in facets if v in f['vertices']]
            chosen = next(t for t in combinations(incident, 3) if determinant(t))
            h = max(abs(x) for normal in incident for x in normal)
            w = tuple(sum(z[i] for z in chosen) for i in range(3))
            assert any(w) and sum(abs(x) for x in w) <= 9 * h
            assert all(dot(w, v) > dot(w, p) for p in points if p != v), 'sum of independent normals strictly exposes'
            Q = ceiling_sqrt(9 * h)
            q = short_relation(w, Q)
            A = kernel_projection(q)
            assert project(A, v) in convex_hull(project(A, p) for p in points), 'vertex visible in independently computed short rational shadow'
            visible += 1
            certificates.append({'vertex': v, 'incident_normals': incident, 'chosen_normals': chosen,
                                 'h': h, 'strict_objective': w, 'short_kernel': q, 'kernel_bound': Q})
        thresholds = []
        for h in (1, 2, 4, 8, 16):
            long = [f for f in facets if max(abs(x) for x in f['normal']) > h]
            bad = {v for f in long for v in f['vertices']}
            assert len(bad) <= sum(len(f['vertices']) for f in long)
            thresholds.append({'h': h, 'bad_vertices': len(bad), 'good_vertices': len(vertices) - len(bad)})
        result.append({'name': name, 'points': points, 'facets': facets, 'vertex_visibility': certificates, 'thresholds': thresholds})
    return {'polytope_fixtures': len(result), 'strict_short_shadow_vertex_certificates': visible, 'fixtures': result}


def lattice_direction_controls():
    rows = []
    for k in (1, 2, 6, 12, 24):
        for T in (((1, 0), (0, 1)), ((13, 8), (8, 5)), ((2, 3), (0, 7))):
            vectors = [project(T, (-i, i * i)) for i in range(1, k + 1)]
            control = direction_control(vectors, abs(det2(*T)))
            control['matrix'] = T
            rows.append(control)
    path = [(8 - i, i * i) for i in range(9)]
    values = [6 * x + y for x, y in path]
    differences = [b - a for a, b in zip(values, values[1:])]
    sides = [Fraction(v) < Fraction(87, 2) for v in values]
    assert sum(a * b < 0 for a, b in zip(differences, differences[1:])) == 1
    assert sum(a != b for a, b in zip(sides, sides[1:])) == 2
    return {'cases': len(rows), 'fixtures': rows,
            'twice_crossed_line_control': {'convex_path': path, 'functional': [6, 1],
                                          'level': '87/2', 'values': values, 'crossings': 2}}


def floor_cuberoot(n):
    lower, upper = 0, 1 << ((n.bit_length() + 2) // 3)
    while lower < upper:
        middle = (lower + upper + 1) // 2
        if middle ** 3 <= n:
            lower = middle
        else:
            upper = middle - 1
    assert lower ** 3 <= n < (lower + 1) ** 3
    return lower


def dyadic_cover_controls():
    fixtures, pair_checks, ground_checks = [], 0, 0
    scale = 1 << 12
    for d, H in product((2, 3, 4), (1, 2, 4, 7, 16, 31)):
        A = (tuple(i + 1 for i in range(d)), tuple((i + 1) ** 2 + 1 for i in range(d)))
        boxes = cover(A, H)
        lines = shared_lines(boxes)
        assert len(lines) <= 2 * d * (d - 1) * H.bit_length()
        for original in hyperbolic_cross(H, d):
            v = project(A, original)
            assert any(all(2 * abs(x) < w for x, w in zip(project(b['T'], v), b['widths'])) for b in boxes)
            ground_checks += 1
        for i, j in combinations(range(d), 2):
            upper = sum(floor_cuberoot(b['sides'][i] * b['sides'][j] * scale ** 3) + 1 for b in boxes)
            lower_rhs = 2 * 5 ** (d - 2) * H.bit_length() * floor_cuberoot(H * scale ** 3)
            assert upper <= lower_rhs, 'exact rational enclosures for the uniform weighted box sum'
            pair_checks += 1
        U = ((2, 3), (0, 7))
        transformed = tuple(tuple(sum(U[i][k] * A[k][j] for k in range(2)) for j in range(d)) for i in range(2))
        delta, height = projection_data(A)
        assert projection_data(transformed) == (14 * delta, height)
        fixtures.append({'dimension': d, 'H': H, 'boxes': len(boxes), 'distinct_boundary_lines': len(lines),
                         'line_bound': 2 * d * (d - 1) * H.bit_length(), 'projection': A,
                         'normalized_height': height})
    return {'fixtures': fixtures, 'weighted_coordinate_pair_sums': pair_checks,
            'projected_ground_cover_checks': ground_checks,
            'cube_root_enclosures': 'exact rational upper and lower bounds with denominator 4096; no floating-point comparison'}


def main():
    if not __debug__:
        raise RuntimeError('Run without -O: exact assertions are required.')
    matroids = matroid_controls()
    short = short_relation_controls()
    geometry = geometry_controls()
    directions = lattice_direction_controls()
    dyadic = dyadic_cover_controls()
    power, log_power = Fraction(9, 8), Fraction(-3, 8)
    assert 1 + Fraction(5, 3) * power == 4 - power == Fraction(23, 8)
    assert 1 + Fraction(5, 3) * log_power == -log_power == Fraction(3, 8)
    assert Fraction(23, 8) + Fraction(2, 3) == Fraction(85, 24)
    sources = ('scripts/run_rational_shadow_bound_checks.py', 'scripts/run_planar_dyadic_bound_checks.py',
               'scripts/run_extremal_initial_ideals_checks.py', 'scripts/run_structural_complexity_checks.py')
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_rational_shadow_bound_checks.py',
              'uniform_proof': 'docs/rational_shadow_bound.md',
              'spatial_count_bound': 'G_3(n)=O(n^(23/8)*log(n+1)^(3/8))',
              'universal_basis_cardinality_bound': 'O(n^(85/24)*log(n+1)^(3/8))',
              'height_dependent_matroid_shadows': matroids,
              'short_integer_relations': short, 'lattice_polytope_visibility': geometry,
              'oblique_lattice_direction_controls': directions,
              'shared_lines_and_weighted_dyadic_sums': dyadic,
              'sharp_asymptotic_growth_established': False, 'novelty_established': False, 'goal_complete': False,
              'scope': 'Exact supporting checks for a uniform rational-shadow theorem and a three-variable upper bound. These do not settle the sharp extremal target, establish current novelty, or imply a running-time bound.',
              'source_sha256': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in sources}}
    out = ROOT / 'results/research/rational_shadow_bound_checks.json'
    out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'quadrant_sweeps': matroids['quadrant_sweeps'],
                      'independent_full_hulls': matroids['independent_full_basis_hulls'],
                      'short_integer_vectors': short['integer_vectors_checked'],
                      'strict_vertex_visibility_certificates': geometry['strict_short_shadow_vertex_certificates'],
                      'oblique_direction_cases': directions['cases'], 'goal_complete': False,
                      'output': str(out.relative_to(ROOT))}))


if __name__ == '__main__':
    main()
