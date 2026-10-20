#!/usr/bin/env python3
"""Exact controls for rational shadows and short-normal vertex visibility.

The uniform proofs are in docs/rational_shadow_bound.md. These finite
controls do not establish asymptotic sharpness, novelty, or goal completion.
"""
from collections import defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations, product
from math import comb, gcd, isqrt, prod
from pathlib import Path
from random import Random
import hashlib
import json

from run_extremal_initial_ideals_checks import rectangle_model
from run_planar_dyadic_bound_checks import convex_hull, matrix_oracle
from run_structural_complexity_checks import cross, determinant, dot, hyperbolic_cross, primitive
from run_structural_complexity_checks import area2, clip
from run_structural_complexity_checks import cyclic_univariate_model

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
                ('octahedron', [tuple(s * int(i == j) for i in range(3)) for j in range(3) for s in (-1, 1)]),
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


def product_directions(T):
    directions = set()
    for a in range(T):
        for b in range(T//(a+1)):
            for c in range(T//((a+1)*(b+1))):
                for signs in product((-1,1),repeat=3):
                    q = (signs[0]*a,signs[1]*b,signs[2]*c)
                    if any(q) and gcd(*q) == 1:
                        directions.add(primitive(q))
    return sorted(directions)


def wall_product_controls(geometry):
    moments = []
    scale = 4096
    for T in (1,2,4,8,16,32,64,128):
        directions = product_directions(T)
        assert all(prod(abs(x)+1 for x in q) <= T and gcd(*q) == 1 for q in directions)
        upper = sum(scale+floor_cuberoot(sum(abs(x) for x in q)*scale**3)+1 for q in directions)
        lower_rhs = 1536*T*floor_cuberoot(T*scale**3)
        assert upper <= lower_rhs, 'exact cube-root enclosures verify the weighted product-direction moment'
        moments.append({'T':T,'primitive_directions_up_to_sign':len(directions),
                        'weighted_moment_upper_numerator':upper,'denominator':scale,
                        'uniform_rhs_lower_numerator':lower_rhs})

    fixtures,ledger,examples = [],[],[]
    for item in geometry['fixtures']:
        vertices = sorted({tuple(p) for f in item['facets'] for p in f['vertices']})
        edges = {tuple(sorted((tuple(a),tuple(b)))) for f in item['facets']
                 for a,b in zip(f['vertices'],f['vertices'][1:]+f['vertices'][:1])}
        fixtures.append((item['name'],vertices,edges))
    fixtures += [('segment',[(0,0,0),(4,8,0)],{((0,0,0),(4,8,0))}),
                 ('point',[(3,2,1)],set()),
                 ('planar_square',[(0,0,0),(2,0,0),(2,3,0),(0,3,0)],
                  {((0,0,0),(2,0,0)),((2,0,0),(2,3,0)),((0,3,0),(2,3,0)),((0,0,0),(0,3,0))})]
    counts = {'polytope_kernel_checks':0,'shadow_vertex_fibers':0,'collapsed_edge_fibers':0,
              'two_endpoint_oracle_checks':0,'positive_endpoint_objectives':0}
    summaries = []
    for name,vertices,edges in fixtures:
        kernels = {primitive(q) for q in product(range(-2,3),repeat=3) if any(q)}
        kernels |= {primitive(tuple(x-y for x,y in zip(a,b))) for a,b in edges}
        widths = [max(p[i] for p in vertices)-min(p[i] for p in vertices) for i in range(3)]
        R = max(widths)
        for q in sorted(kernels):
            A = kernel_projection(q)
            hull = convex_hull(project(A,p) for p in vertices)
            recovered = set()
            for j,y in enumerate(hull):
                if len(hull) >= 3:
                    before,after = hull[j-1],hull[(j+1)%len(hull)]
                    u,v = tuple(a-b for a,b in zip(y,before)),tuple(a-b for a,b in zip(after,y))
                    objective = (u[1]+v[1],-u[0]-v[0])
                elif len(hull) == 2:
                    objective = tuple(a-b for a,b in zip(y,hull[1-j]))
                else:
                    objective = (1,0)
                assert all(dot(objective,y)>dot(objective,z) for z in hull if z!=y)
                w = tuple(sum(objective[k]*A[k][i] for k in range(2)) for i in range(3))
                fiber = {p for p in vertices if project(A,p)==y}
                assert 1 <= len(fiber) <= 2, 'a vertex fiber has at most two original vertices, not necessarily one'
                epsilon = Fraction(1,2*(R+1)*sum(abs(x) for x in q)+1)
                endpoints = set()
                for sign in (-1,1):
                    perturbed = tuple(w[i]+sign*epsilon*q[i] for i in range(3))
                    endpoint = max(vertices,key=lambda p:dot(perturbed,p))
                    assert all(dot(perturbed,endpoint)>dot(perturbed,p) for p in vertices if p!=endpoint)
                    endpoints.add(endpoint)
                    counts['two_endpoint_oracle_checks'] += 1
                    if min(w)>0:
                        assert min(perturbed)>0, 'the explicit perturbation preserves positive objectives'
                        counts['positive_endpoint_objectives'] += 1
                assert endpoints==fiber, 'two ordinary support calls recover the complete exposed fiber'
                recovered |= endpoints
                counts['shadow_vertex_fibers'] += 1
                if len(fiber)==2:
                    counts['collapsed_edge_fibers'] += 1
                    if len(examples)<8:
                        examples.append({'fixture':name,'kernel':q,'shadow_vertex':y,
                                         'fiber':sorted(fiber),'objective':w,'epsilon':str(epsilon)})
                ledger.append((name,q,y,sorted(fiber)))
            incident = {p for edge in edges if primitive(tuple(x-y for x,y in zip(*edge)))==q for p in edge}
            assert incident <= recovered and len(incident)<=2*len(hull), 'every parallel-edge endpoint is covered by a shadow fiber'
            counts['polytope_kernel_checks'] += 1
        summaries.append({'name':name,'vertices':len(vertices),'edges':len(edges),'kernels':len(kernels)})

    box_checks = 0
    for H in (3,7,16,31):
        ground = tuple(hyperbolic_cross(H,3))
        for box in cover(((1,0,0),(0,1,1)),H):
            points = [a for a in ground if all(x<=b for x,b in zip(a,box['sides']))]
            for a,b in combinations(points,2):
                q = primitive(tuple(x-y for x,y in zip(a,b)))
                assert prod(abs(x)+1 for x in q)<=prod(x+1 for x in box['sides'])<=8*H
                box_checks += 1

    axis_controls = []
    for n in (6,9,12,16):
        r = n-1
        for a in range(1,r-1):
            for b in range(1,r-a):
                allocation = (a,b,r-a-b)
                normals = [tuple((allocation[i]+1 if k==i else 0)-(allocation[j] if k==j else 0)
                                 for k in range(3)) for i,j in product(range(3),repeat=2) if i!=j]
                polygon = [tuple(Fraction(i==j) for i in range(3)) for j in range(3)]
                for normal in normals:
                    polygon = clip(polygon,normal,1)
                assert polygon and area2(polygon)>0 and all(min(p)>0 for p in polygon)
                actual = {primitive(normal) for normal in normals if any(
                    not dot(normal,p) and not dot(normal,q) for p,q in zip(polygon,polygon[1:]+polygon[:1]))}
                assert actual=={primitive(v) for v in normals}, 'all six transfer comparisons are actual interior facets'
                products = [prod(abs(x)+1 for x in v) for v in actual]
                axis_controls.append({'n':n,'allocation':allocation,'interior_facets':len(actual),
                                      'minimum_primitive_wall_product':min(products)})
    sieves = []
    G = 4096
    for n in (16,32,64,128,256):
        lo,hi = (n+3)//4,3*n//8
        allocations = [(a,b,n-1-a-b) for a in range(lo,hi+1) for b in range(lo,hi+1)]
        assert len(allocations)*256>=n*n
        good = 0
        for allocation in allocations:
            assert all(n<=8*a and 2*a<=n for a in allocation)
            pairs = [(allocation[i]+1,allocation[j]) for i,j in product(range(3),repeat=2) if i!=j]
            if all(gcd(a,b)<=G for a,b in pairs):
                good += 1
                assert all(64*G*G*(a//gcd(a,b)+1)*(b//gcd(a,b)+1)>=n*n for a,b in pairs)
        assert (len(allocations)-good)*G<=6*n*n and good*2048>=5*n*n
        # Check the integer tail inequality used in the uniform sieve at nontrivial cutoffs too.
        for cutoff in (1,2,4,8):
            tail = sum((n//d)**2 for d in range(cutoff+1,n+1))
            assert cutoff*tail<=n*n
        sieves.append({'n':n,'central_rectangle_allocations':len(allocations),
                       'good_allocations_at_G_4096':good,'gcd_cutoff':G})
    return {'weighted_direction_moments':moments,'projection_fiber_summary':counts,
            'projection_fiber_fixtures':summaries,'collapsed_fiber_examples':examples,
            'fiber_ledger_sha256':hashlib.sha256(json.dumps(ledger).encode()).hexdigest(),
            'same_box_direction_checks':box_checks,'complete_central_axis_cells':axis_controls,
            'uniform_axis_gcd_sieve_controls':sieves,
            'scope':'Uniform partial count O(n*T^(4/3)*log(n+1)) for initial cells incident to an interior wall with primitive direction product at most T. Exact fiber recovery includes endpoints of collapsed edges. The known three-pole quotient has quadratically many interior cells whose every primitive wall product is Omega(n^2). No bound for all remaining cells, unrestricted G_3 improvement, new enumeration implementation or novelty certificate.'}


def ordered_flag_controls():
    """Exact rank-oracle identities on fixed internal orders; no fan sweep."""
    def compile_lists(groups,raw_rank):
        @cache
        def rank(mask):
            return raw_rank(mask)
        original = groups
        compressed = []
        for group in groups:
            mask,kept = 0,[]
            for e in group:
                if rank(mask | (1<<e)) > len(kept):
                    mask |= 1<<e
                    kept.append(e)
            compressed.append(tuple(kept))
        groups = tuple(compressed)+((),)*(3-len(compressed))
        prefixes = []
        for group in groups:
            row = [0]
            for e in group:
                row.append(row[-1] | (1<<e))
            prefixes.append(row)
        r = rank(sum(1<<e for group in original for e in group))
        assert max(map(len,groups),default=0)<=r
        terms = {tuple([e]):1 for group in groups for e in group}
        pairs = []
        for a,b in combinations(range(3),2):
            support = []
            for i,e in enumerate(groups[a],1):
                for j,f in enumerate(groups[b],1):
                    delta = (rank(prefixes[a][i]|prefixes[b][j-1])+
                             rank(prefixes[a][i-1]|prefixes[b][j])-
                             rank(prefixes[a][i-1]|prefixes[b][j-1])-
                             rank(prefixes[a][i]|prefixes[b][j]))
                    assert delta in (0,1)
                    if delta:
                        support.append((i,j))
                        terms[(e,f)] = -1
            assert len({i for i,j in support})==len({j for i,j in support})==len(support)
            assert len(support)==len(groups[a])+len(groups[b])-rank(prefixes[a][-1]|prefixes[b][-1])
            pairs.append({'groups':[a,b],'matching':support})
        intervals = []
        length = len(groups[2])
        for i,e in enumerate(groups[0],1):
            for j,f in enumerate(groups[1],1):
                base = prefixes[0][i-1]|prefixes[1][j-1]

                def gain(k,mask):
                    context = base|prefixes[2][k]
                    return rank(context|mask)-rank(context)

                def first(mask,threshold):
                    low,high = 0,length+1
                    while low<high:
                        middle = (low+high)//2
                        if middle==length+1 or gain(middle,mask)<=threshold:
                            high = middle
                        else:
                            low = middle+1
                    return low

                birth = first((1<<e)|(1<<f),1)
                death = min(first(1<<e,0),first(1<<f,0))
                if birth<death:
                    intervals.append({'pair':[i,j],'birth':birth,'death':death})
                    if birth:
                        terms[(e,f,groups[2][birth-1])] = -1
                    if death<=length:
                        terms[(e,f,groups[2][death-1])] = 1
        calls = rank.cache_info().misses
        assert sum(len(t)==3 for t in terms)<=2*r*r
        assert sum(terms.values())==r
        return {'groups':groups,'terms':terms,'pair_matchings':pairs,'intervals':intervals,
                'preprocessing_rank_calls':calls,'rank':rank,'total_rank':r}

    def shuffles(groups):
        if not any(groups):
            yield ()
        for i,group in enumerate(groups):
            if group:
                rest = groups[:i]+(group[1:],)+groups[i+1:]
                for suffix in shuffles(rest):
                    yield (group[0],)+suffix

    def verify(data,groups,orders):
        prefixes_checked = queries = 0
        compressed,terms,rank = data['groups'],data['terms'],data['rank']
        for lengths in product(*(range(len(group)+1) for group in compressed)):
            selected = {e for group,k in zip(compressed,lengths) for e in group[:k]}
            mask = sum(1<<e for e in selected)
            assert rank(mask)==sum(c for term,c in terms.items() if set(term)<=selected)
            prefixes_checked += 1
        for order in orders:
            costs = {e:i+1 for i,e in enumerate(order)}
            predicted = defaultdict(int)
            for term,c in terms.items():
                predicted[max(term,key=lambda e:(costs[e],e))] += c
            assert all(c in (0,1) for c in predicted.values()), 'signed derivatives give an actual basis indicator'
            mask,basis = 0,set()
            for e in order:
                if rank(mask|(1<<e))>len(basis):
                    mask |= 1<<e
                    basis.add(e)
            assert basis=={e for e,c in predicted.items() if c}
            assert sum(costs[e] for e in basis)==sum(c*max(costs[e] for e in term) for term,c in terms.items())
            queries += 1
        return prefixes_checked,queries

    def labels(prefix,maximum):
        if len(prefix)==6:
            if maximum>=2:
                yield tuple(prefix)
        else:
            for c in range(maximum+2):
                yield from labels(prefix+[c],max(maximum,c))

    two_cases = three_cases = prefixes_checked = queries = 0
    ledger = []
    for colors in labels([],0):
        def rank(mask):
            return min(2,len({colors[i] for i in range(6) if mask>>i&1 and colors[i]}))
        for groups in (((0,1,2),(3,4,5)),((0,1),(2,3),(4,5))):
            data = compile_lists(groups,rank)
            p,q = verify(data,groups,shuffles(groups))
            prefixes_checked += p; queries += q
            two_cases += len(groups)==2; three_cases += len(groups)==3
            ledger.append((colors,groups,sorted(data['terms'].items())))
    assert two_cases==three_cases==813
    linear = []
    for r,p in product((3,4),(2,3,101)):
        rng = Random(7400+100*r+p)
        columns = [tuple(int(i==j) for i in range(r)) if j<r else
                   tuple(rng.randrange(p) for i in range(r)) for j in range(9)]
        independent = matrix_oracle(columns,r,p)

        def rank(mask):
            selected = frozenset()
            for e in range(9):
                if mask>>e&1 and independent(selected|{e}):
                    selected = selected|{e}
            return len(selected)

        groups = tuple(tuple(range(a,9,3)) for a in range(3))
        data = compile_lists(groups,rank)
        pcount,qcount = verify(data,groups,shuffles(groups))
        prefixes_checked += pcount; queries += qcount
        linear.append({'rank':r,'field':p,'original_lists':groups,'compressed_lists':data['groups'],
                       'preprocessing_rank_calls':data['preprocessing_rank_calls'],
                       'triple_terms':sum(len(t)==3 for t in data['terms']),'queries':qcount})

    triangle_groups = ((0,),(1,),(2,3))
    triangle_colors = (1,2,3,1)
    triangle = compile_lists(triangle_groups,lambda mask:min(2,len({triangle_colors[i] for i in range(4) if mask>>i&1})))
    pcount,qcount = verify(triangle,triangle_groups,shuffles(triangle_groups))
    prefixes_checked += pcount; queries += qcount
    assert triangle['terms'][(0,1,2)]==-1 and triangle['terms'][(0,1,3)]==1
    assert triangle['intervals']==[{'pair':[1,1],'birth':1,'death':2}]
    uniform = []
    for r in range(1,9):
        groups = tuple(tuple(range(a*r,(a+1)*r)) for a in range(3))
        data = compile_lists(groups,lambda mask:min(r,mask.bit_count()))
        triple = {t:c for t,c in data['terms'].items() if len(t)==3}
        assert len(triple)==r*r
        for term,c in triple.items():
            level = term[0]+1+term[1]-r+1+term[2]-2*r+1
            assert (level,c) in ((r+1,-1),(r+2,1)), 'the full quadratic correction is present in the uniform prefix model'
        pcount,_ = verify(data,groups,())
        prefixes_checked += pcount
        allocations = []
        for a in range(r+1):
            for b in range(r-a+1):
                allocation = a,b,r-a-b
                costs = {e:Fraction(2*(e%r+1),2*allocation[e//r]+1) for e in range(3*r)}
                order = tuple(sorted(costs,key=lambda e:(costs[e],e)))
                _,qcount = verify_query_only(data,order,costs)
                queries += qcount
                actual = sorted(order[:r])
                assert tuple(sum(e//r==axis for e in actual) for axis in range(3))==allocation
                allocations.append(allocation)
        assert len(allocations)==comb(r+2,2)
        uniform.append({'rank':r,'triple_terms':len(triple),'negative_triple_terms':sum(c<0 for c in triple.values()),
                        'positive_triple_terms':sum(c>0 for c in triple.values()),
                        'preprocessing_rank_calls':data['preprocessing_rank_calls'],
                        'axis_allocation_queries':len(allocations),
                        'maximum_states_per_two_list_region':r+1,
                        'regions_required_for_whole_axis_fan_at_least':(r+3)//2})
    direct_prefixes = 0
    for n,p in product(range(3,7),(0,2,3,5)):
        modulus = [1,-1]+[0]*(n-2)+[1]
        images = [(0,1)+(0,)*(n-2),(1,)+(0,)*(n-2)+(-1,),(0,)+(-1,)*(n-1)]
        ground = {(0,0,0)} | {tuple(k if axis==i else 0 for i in range(3))
                               for axis in range(3) for k in range(1,n)}
        columns = cyclic_univariate_model(modulus,images,ground)
        for lengths in product(range(n),repeat=3):
            selected = {(0,0,0)} | {tuple(k if axis==i else 0 for i in range(3))
                                      for axis in range(3) for k in range(1,lengths[axis]+1)}
            pivots = {}
            for a in sorted(selected):
                vector = [x%p if p else Fraction(x) for x in columns[a]]
                for i,row in sorted(pivots.items()):
                    c = vector[i]
                    if c:
                        vector = [(x-c*y)%p if p else x-c*y for x,y in zip(vector,row)]
                i = next((i for i,x in enumerate(vector) if x),None)
                if i is not None:
                    inverse = pow(vector[i],-1,p) if p else 1/vector[i]
                    pivots[i] = [(x*inverse)%p if p else x*inverse for x in vector]
            assert len(pivots)==min(n,1+sum(lengths)), 'ordinary quotient columns give the asserted three-axis prefix ranks'
            direct_prefixes += 1
    return {'all_rank_two_matroids_on_six_elements':two_cases,'three_list_rank_two_models':three_cases,
            'prefix_rank_reconstructions':prefixes_checked,'direct_greedy_queries':queries,
            'linear_flag_models':linear,'uniform_quadratic_correction_models':uniform,
            'parallel_interval_example':{'lists':triangle_groups,'parallel_classes':triangle_colors,
                                         'intervals':triangle['intervals'],
                                         'signed_terms':[[term,c] for term,c in sorted(triangle['terms'].items())]},
            'direct_three_pole_prefix_field_checks':direct_prefixes,
            'rank_two_ledger_sha256':hashlib.sha256(json.dumps(ledger).encode()).hexdigest(),
            'scope':'Two fixed internal orders have a direct matroid-interleaving prior-art boundary. Three-list prefix ranks admit signed pair/triple max formulas, with at most 2*r^2 triple terms and O(r^2*log(r+1)) rank-oracle preprocessing after list compression. Uniform prefix rank already has exactly r^2 triple terms; ordinary three-pole quotient arithmetic verifies that prefix model after contracting the constant. This does not identify its full residue matroid with a uniform matroid. Any two-list region covers at most n states of its positive fan, so the quadratic large-wall subset requires Omega(n) such regions. No polylogarithmic cover, unrestricted G_3 bound, strongest-current or novelty certificate.'}


def verify_query_only(data,order,costs):
    predicted = defaultdict(int)
    for term,c in data['terms'].items():
        predicted[max(term,key=lambda e:(costs[e],e))] += c
    assert all(c in (0,1) for c in predicted.values())
    mask,basis = 0,set()
    for e in order:
        if data['rank'](mask|(1<<e))>len(basis):
            mask |= 1<<e
            basis.add(e)
    assert basis=={e for e,c in predicted.items() if c}
    assert sum(costs[e] for e in basis)==sum(c*max(costs[e] for e in term) for term,c in data['terms'].items())
    return basis,1


def main():
    if not __debug__:
        raise RuntimeError('Run without -O: exact assertions are required.')
    matroids = matroid_controls()
    short = short_relation_controls()
    geometry = geometry_controls()
    directions = lattice_direction_controls()
    dyadic = dyadic_cover_controls()
    walls = wall_product_controls(geometry)
    flags = ordered_flag_controls()
    power, log_power = Fraction(9, 8), Fraction(-3, 8)
    assert 1 + Fraction(5, 3) * power == 4 - power == Fraction(23, 8)
    assert 1 + Fraction(5, 3) * log_power == -log_power == Fraction(3, 8)
    assert Fraction(23, 8) + Fraction(2, 3) == Fraction(85, 24)
    sources = ('docs/rational_shadow_bound.md','scripts/run_rational_shadow_bound_checks.py', 'scripts/run_planar_dyadic_bound_checks.py',
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
              'partial_wall_product_count':walls,
              'ordered_flag_interaction_controls':flags,
              'sharp_asymptotic_growth_established': False, 'novelty_established': False, 'goal_complete': False,
              'scope': 'Exact supporting checks for a uniform rational-shadow theorem and a three-variable upper bound. A partial wall-product count, explicit projection-fiber recovery and a uniform large-wall exception family are also checked. Fixed-list matching, signed three-list ranks and the obstruction to a small region cover are independently verified. These do not settle the sharp extremal target, establish current novelty, or imply an implemented running-time bound.',
              'source_sha256': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in sources}}
    out = ROOT / 'results/research/rational_shadow_bound_checks.json'
    out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'quadrant_sweeps': matroids['quadrant_sweeps'],
                      'independent_full_hulls': matroids['independent_full_basis_hulls'],
                      'short_integer_vectors': short['integer_vectors_checked'],
                      'strict_vertex_visibility_certificates': geometry['strict_short_shadow_vertex_certificates'],
                      'oblique_direction_cases': directions['cases'], 'goal_complete': False,
                      'wall_product_fiber_checks':walls['projection_fiber_summary']['shadow_vertex_fibers'],
                      'ordered_flag_greedy_checks':flags['direct_greedy_queries'],
                      'output': str(out.relative_to(ROOT))}))


if __name__ == '__main__':
    main()
