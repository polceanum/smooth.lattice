#!/usr/bin/env python3
"""Exact controls for weighted facet incidence and an area-budget obstruction."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import gcd, isqrt
from pathlib import Path
import hashlib
import json

from run_extremal_initial_ideals_checks import cone_for_basis, rectangle_model
from run_initial_cone_height_checks import (
    add, clean, integer_ray, leading, normal_form, subtract_shift,
)
from run_planar_dyadic_bound_checks import convex_hull
from run_adjacent_rank_checks import rank
from run_oracle_enumeration_checks import hull_facets
from run_structural_complexity_checks import (
    area2, cross, determinant, dot, greedy, hyperbolic_cross,
)

ROOT = Path(__file__).resolve().parents[1]


def fan_incidence(n, ground, columns, residue, witnesses, prime=101):
    Q = 12 * n * n + 1
    planes = [tuple(Q if i == j else -1 for i in range(3)) for j in range(3)]
    incidence, states, polygons = defaultdict(set), [], []
    for index, witness in enumerate(witnesses):
        basis = set(map(tuple, witness['basis']))
        assert greedy(ground, columns, witness['strict_weight'], prime) == basis
        polygon, normals = cone_for_basis(basis, residue, prime)
        assert all(abs(x) <= n for z in normals for x in z)
        rays = [integer_ray(p) for p in polygon]
        triple = next(t for t in combinations(rays, 3) if determinant(t))
        certificate = tuple(sum(r[i] for r in triple) for i in range(3))
        assert 0 < min(certificate) and max(certificate) <= 6 * n * n
        assert all(dot(z, certificate) > 0 for z in normals)
        for ray in rays:
            if min(ray) > 0:
                incidence[ray].add(index)
        boundary = next(((i, r) for r in rays for i in range(3) if r[i] == 0), None)
        state = {'basis': sorted(basis), 'strict_certificate': certificate}
        if boundary:
            i, ray = boundary
            q = planes[i]
            a, b = -dot(q, ray), dot(q, certificate)
            assert a > 0 and b > 0
            w = tuple(a * x + b * y for x, y in zip(certificate, ray))
            factor = gcd(*w)
            w = tuple(x // factor for x in w)
            assert dot(q, w) == 0 and min(w) > 0
            assert all(dot(z, w) > 0 for z in normals)
            assert greedy(ground, columns, w, prime) == basis
            state['boundary_cover'] = {'plane': i, 'objective': w}
        states.append(state)
        polygons.append(polygon)
    assert sum(area2(p) for p in polygons) == 1
    assert len({tuple(map(tuple, s['basis'])) for s in states}) == len(states)
    T = sum(max(w) * len(cells) for w, cells in incidence.items())
    return {'length': n, 'initial_ideals': len(states), 'weighted_interior_incidence': T,
            'interior_rays': [{'weight': w, 'height': max(w), 'cells': sorted(cells)}
                              for w, cells in sorted(incidence.items())],
            'boundary_cover_planes': planes,
            'boundary_cells_checked': sum('boundary_cover' in s for s in states), 'states': states}


def quartic_model(k):
    n = 4 * k
    reference = tuple((a, 0, 0) for a in range(1, 2 * k + 1)) + tuple((0, b, 0) for b in range(2 * k))
    indices = {a: i for i, a in enumerate(reference)}

    def reduced(a):
        x, y, z = a
        if x:
            if z or y >= 2:
                return None
            exponent = x + k * y
            return (exponent, 0, 0) if exponent <= 2 * k else None
        exponent = y + k * z
        return (0, exponent, 0) if exponent < 2 * k else None

    def residue(a):
        target = reduced(a)
        return tuple(int(target == b) for b in reference)

    maps = [tuple(None if reduced(add(a, tuple(int(i == j) for i in range(3)))) is None
                  else indices[reduced(add(a, tuple(int(i == j) for i in range(3))))]
                  for a in reference) for j in range(3)]

    def image(a, start):
        for i, power in enumerate(a):
            for _ in range(power):
                start = None if start is None else maps[i][start]
        return start

    for i, j in combinations(range(3), 2):
        assert all((None if maps[j][v] is None else maps[i][maps[j][v]]) ==
                   (None if maps[i][v] is None else maps[j][maps[i][v]]) for v in range(n))
    relations = [{(1, 1, 0): 1, (k + 1, 0, 0): -1}, {(2 * k + 1, 0, 0): 1},
                 {(0, 2 * k, 0): 1}, {(0, 0, 1): 1, (0, k, 0): -1}]
    for g in relations:
        for start in range(n):
            result = [0] * n
            for a, c in g.items():
                target = image(a, start)
                if target is not None:
                    result[target] += c
            assert not any(result)
    assert all(image(a, indices[(0, 0, 0)]) == i for i, a in enumerate(reference))
    ground = tuple(hyperbolic_cross(n, 3))
    columns = {a: residue(a) for a in ground}
    weight = (1, k + 1, k * (k + 1) + 1)
    s_pairs = 0
    for characteristic in (0, 2, 3, 101):
        basis = [clean(g, characteristic) for g in relations]
        heads = [leading(g, weight) for g in basis]
        for i, j in combinations(range(4), 2):
            common = tuple(max(a, b) for a, b in zip(heads[i], heads[j]))
            first = {add(a, tuple(x - y for x, y in zip(common, heads[i]))): c
                     for a, c in basis[i].items()}
            pair = subtract_shift(first, basis[j], tuple(x - y for x, y in zip(common, heads[j])),
                                  1, characteristic)
            assert not normal_form(pair, basis, weight, characteristic)
            s_pairs += 1
    rational_basis = [clean(g, 0) for g in relations]
    for a in ground:
        expected = {} if reduced(a) is None else {reduced(a): 1}
        assert normal_form({a: 1}, rational_basis, weight, 0) == expected
        target = image(a, indices[(0, 0, 0)])
        assert columns[a] == tuple(int(i == target) for i in range(n))

    # Enumerate all profile sums by actual residue classes, without cost selection.
    classes = defaultdict(list)
    for a in ground:
        if any(columns[a]):
            classes[columns[a]].append(a)
    assert len(classes) == n
    profiles = {(0, 0, 0)}
    for choices in classes.values():
        profiles = {add(a, b) for a in profiles for b in choices}
    base = (k * (2 * k + 1), k * (2 * k - 1), 0)
    v1, v2 = (-k, 1, 0), (0, -k, 1)
    expected = {tuple(base[i] + a * v1[i] + b * v2[i] for i in range(3))
                for a, b in product(range(k + 1), repeat=2)}
    assert profiles == expected
    vertices = [tuple(base[i] + a * k * v1[i] + b * k * v2[i] for i in range(3))
                for a, b in product(range(2), repeat=2)]
    assert set(convex_hull([p[:2] for p in profiles])) == {p[:2] for p in vertices}
    w = (1, k, k * k)
    assert cross(v1, v2) == w
    assert len(profiles) == (k + 1)**2
    witnesses = []
    for cost in [weight, (2, 1, k + 1), (1, k + 1, 1), (2, 1, 1)]:
        selected = greedy(ground, columns, cost, 101)
        witnesses.append({'basis': sorted(selected), 'strict_weight': cost})
    fan = fan_incidence(n, ground, columns, residue, witnesses)
    assert fan['weighted_interior_incidence'] == 4 * k * k
    return {'k': k, 'length': n, 'ground_normal_forms': len(ground), 's_pairs': s_pairs,
            'basis_count': 2**(2 * k), 'profile_count': len(profiles),
            'profile_vertices': vertices, 'primitive_normal': w,
            'facet_area_squared': k**4 * dot(w, w), 'facet_lattice_area': k * k,
            'regularized_geometry': regularized_geometry(vertices, k), 'fan': fan}


def regularized_geometry(vertices, k):
    points = sorted({add(a, b) for a in vertices for b in product(range(2), repeat=3)})
    facets = {}
    for a, b, c in combinations(points, 3):
        normal = cross(tuple(x - y for x, y in zip(b, a)), tuple(x - y for x, y in zip(c, a)))
        factor = gcd(*normal)
        if not factor:
            continue
        normal = tuple(x // factor for x in normal)
        values = [dot(normal, p) - dot(normal, a) for p in points]
        if min(values) < 0 < max(values):
            continue
        if min(values) >= 0:
            normal = tuple(-x for x in normal)
        offset = dot(normal, a)
        facets[(normal, offset)] = [p for p in points if dot(normal, p) == offset]
    large, bad = [], set()
    for (normal, offset), face in facets.items():
        pivot = max(range(3), key=lambda i: abs(normal[i]))
        axes = [i for i in range(3) if i != pivot]
        projected = {(p[axes[0]], p[axes[1]]): p for p in face}
        hull = convex_hull(list(projected))
        if 2 * max(map(abs, normal)) > k * k:
            expected = (1, k, k * k)
            assert normal in (expected, tuple(-x for x in expected))
            assert len(hull) == 4
            doubled_area = abs(sum(p[0] * q[1] - p[1] * q[0]
                                   for p, q in zip(hull, hull[1:] + hull[:1])))
            area_squared = Fraction(doubled_area**2 * dot(normal, normal), 4 * normal[pivot]**2)
            assert area_squared == k**4 * dot(normal, normal)
            large.append({'normal': normal, 'offset': offset, 'vertices': [projected[p] for p in hull],
                          'area_squared': str(area_squared)})
            bad.update(projected[p] for p in hull)
        else:
            assert max(map(abs, normal)) <= k
    assert len(large) == 2 and len(bad) == 8
    return {'facets': len(facets), 'threshold': str(Fraction(k * k, 2)),
            'large_normal_facets': large, 'bad_vertices': len(bad)}


def ceil_sqrt(n):
    r = isqrt(n)
    return r+int(r*r != n)


def polygon_perimeter(poly):
    return sum(sum(abs(a-b) for a, b in zip(u, v))
               for u, v in zip(poly, poly[1:]+poly[:1]))


def corner_incidence_check(poly, normal, counts):
    h = max(map(abs, normal))
    perimeter = polygon_perimeter(poly)
    for j, v in enumerate(poly):
        u = tuple(a-b for a, b in zip(poly[j-1], v))
        z = tuple(a-b for a, b in zip(poly[(j+1) % len(poly)], v))
        q = cross(u, z)
        assert any(q) and cross(q, normal) == (0, 0, 0)
        assert max(map(abs, q)) >= h
        assert h <= sum(map(abs, u))*sum(map(abs, z))
        counts['independent_integer_corner_normal_bounds'] += 1
    assert len(poly)*ceil_sqrt(h) <= 2*perimeter
    counts['independent_facet_half_incidence_perimeter_bounds'] += 1
    return perimeter


def boxed_matroid_skeleton_controls():
    counts = defaultdict(int)
    labels = ((0, 0, 0), (3, 0, 0), (0, 4, 0),
              (0, 0, 2), (2, 1, 1), (1, 3, 1))
    graphic = tuple(tuple(int(i == a)-int(i == b) for i in range(3))
                    for a, b in combinations(range(4), 2))
    cases = [('uniform_three_six', tuple((1, a, a*a) for a in range(6)), 0),
             ('graphic_K4', graphic, 2),
             ('three_parallel_pairs', ((1, 0, 0), (1, 0, 0), (0, 1, 0),
                                       (0, 1, 0), (0, 0, 1), (0, 0, 1)), 3),
             ('uniform_two_six', tuple((1, a) for a in range(6)), 0)]
    fixtures = []
    for name, columns, p in cases:
        r = rank(columns, p)
        profiles = {tuple(sum(labels[j][i] for j in B) for i in range(3))
                    for B in combinations(range(6), r)
                    if rank([columns[j] for j in B], p) == r}
        facets = hull_facets(profiles)
        edges, vertices, U = set(), set(), 0
        for (normal, _), face in facets.items():
            pivot = max(range(3), key=lambda i: abs(normal[i]))
            axes = [i for i in range(3) if i != pivot]
            mapping = {tuple(v[i] for i in axes): v for v in face}
            poly = [mapping[v] for v in convex_hull(mapping)]
            corner_incidence_check(poly, normal, counts)
            U += len(poly)*ceil_sqrt(max(map(abs, normal)))
            vertices.update(poly)
            edges.update(tuple(sorted((a, b))) for a, b in zip(poly, poly[1:]+poly[:1]))
        sides = tuple(max(a[i] for a in labels) for i in range(3))
        vertical = max(range(3), key=lambda i: sides[i])
        axes = [i for i in range(3) if i != vertical]
        project = lambda v: tuple(v[i] for i in axes)
        interior = defaultdict(set)
        vertical_edges = []
        primitive_lengths, nonvertical_k = [], 0
        for a, b in edges:
            difference = tuple(x-y for x, y in zip(b, a))
            k = gcd(*difference)
            q = tuple(x//k for x in difference)
            assert all(abs(x) <= bound for x, bound in zip(q, sides))
            assert any(cross(q, tuple(x-y for x, y in zip(c, d))) == (0, 0, 0)
                       for c, d in combinations(labels, 2) if c != d)
            if all(not q[i] for i in axes):
                vertical_edges.append((a, b))
            else:
                nonvertical_k += k
                primitive_lengths.append(sum(map(abs, q)))
                for j in range(1, k):
                    v = tuple(x+j*z for x, z in zip(a, q))
                    interior[project(v)].add(v)
        assert all(len(v) <= 2 for v in interior.values())
        assert nonvertical_k == sum(map(len, interior.values()))+len(edges)-len(vertical_edges)
        grid = (r*sides[axes[0]]+1)*(r*sides[axes[1]]+1)
        assert len(vertices) <= 2*grid and len(edges) <= 3*len(vertices)
        assert sum(map(len, interior.values())) <= 2*grid
        shadow = set(convex_hull([project(v) for v in vertices]))
        assert all(project(a) == project(b) and project(a) in shadow for a, b in vertical_edges)
        assert len({project(a) for a, _ in vertical_edges}) == len(vertical_edges)
        length = sum(sum(abs(x-y) for x, y in zip(a, b)) for a, b in edges)
        assert U <= 4*length
        counts['independent_boxed_matroid_skeletons'] += 1
        counts['nonvertical_integer_edge_steps'] += nonvertical_k
        fixtures.append({'model': name, 'characteristic': p, 'rank': r,
                         'ground_labels': labels, 'columns': columns,
                         'box_sides': sides, 'profile_vertices': sorted(vertices),
                         'edges': sorted(edges), 'total_l1_edge_length': length,
                         'all_facet_half_incidence': U,
                         'rank_squared_box_volume': r*r*sides[0]*sides[1]*sides[2],
                         'vertical_edges': len(vertical_edges)})
    return {'counts': dict(counts), 'fixtures': fixtures}


def grid_zonotope_controls():
    counts, fixtures = defaultdict(int), []
    for K in range(1, 5):
        grid = tuple(product(range(K, 4*K+1), range(6*K, 9*K+1)))
        G = tuple((a, -b, 1) for a, b in grid)
        planes = defaultdict(set)
        incident = [set() for _ in G]
        for i, j in combinations(range(len(G)), 2):
            normal = cross(G[i], G[j])
            g = gcd(*normal)
            normal = tuple(x//g for x in normal)
            if next(x for x in normal if x) < 0:
                normal = tuple(-x for x in normal)
            planes[normal].update((i, j))
            incident[i].add(normal); incident[j].add(normal)
        positive = {w: S for w, S in planes.items() if min(w) > 0}
        bucketed = {}
        for u, v in product(range(1, 3*K+1), repeat=2):
            if gcd(u, v) != 1:
                continue
            buckets = defaultdict(set)
            for i, (a, b) in enumerate(grid):
                value = b*u-a*v
                if value > 0:
                    buckets[value].add(i)
            for value, S in buckets.items():
                if len(S) >= 2:
                    bucketed[v, u, value] = S
        assert positive == bucketed
        selected = defaultdict(int)
        for a, b, u, v in product(range(K, 2*K+1), range(6*K, 7*K+1),
                                   range(K, 2*K+1), range(K, 2*K+1)):
            if gcd(u, v) == 1:
                w = (v, u, b*u-a*v)
                assert 2*K*K <= w[2] <= 13*K*K and w in positive
                assert len(positive[w]) <= 4
                selected[w] += 1
        assert max(selected.values()) <= 6
        positive_edges, witnesses = {}, []
        U, all_plane_perimeters = 0, 0
        translation = (0, sum(b for a, b in grid), 0)
        for w, S in sorted(planes.items()):
            assert S == {i for i, g in enumerate(G) if dot(w, g) == 0}
            pivot = max(range(3), key=lambda i: abs(w[i]))
            axes = [i for i in range(3) if i != pivot]
            hull = [(0, 0)]
            for i in S:
                q = tuple(G[i][j] for j in axes)
                hull = convex_hull(hull+[tuple(a+b for a, b in zip(v, q)) for v in hull])
            poly = []
            for v in hull:
                point = [0, 0, 0]
                for j, value in zip(axes, v):point[j] = value
                numerator = -sum(w[j]*point[j] for j in axes)
                assert numerator % w[pivot] == 0
                point[pivot] = numerator//w[pivot]
                poly.append(tuple(point))
            assert len(poly) == 2*len(S)
            if len(S) <= 6:
                exhaustive = {(0, 0)}
                for i in S:
                    q = tuple(G[i][j] for j in axes)
                    exhaustive |= {tuple(a+b for a, b in zip(v, q)) for v in exhaustive}
                assert set(convex_hull(exhaustive)) == set(hull)
                counts['independent_subset_sum_face_polygons'] += 1
            perimeter = corner_incidence_check(poly, w, counts)
            all_plane_perimeters += perimeter
            if min(w) <= 0:continue
            U += len(poly)*ceil_sqrt(max(w))
            base = tuple(translation[j]+sum(g[j] for g in G if dot(w, g) < 0)
                         for j in range(3))
            actual = [add(base, v) for v in poly]
            assert all(dot(w, x) == dot(w, base) for x in actual)
            for a, b in zip(actual, actual[1:]+actual[:1]):
                positive_edges[tuple(sorted((a, b)))] = sum(abs(x-y) for x, y in zip(a, b))
            if w in selected and len(witnesses) < 4:
                witnesses.append({'primitive_weight': w, 'zero_generators': [G[i] for i in sorted(S)],
                                  'face_vertices': actual, 'incidences': len(poly)})
        skeleton = sum(2*len(C)*sum(map(abs, g)) for C, g in zip(incident, G))
        assert skeleton == all_plane_perimeters
        visible_length = sum(positive_edges.values())
        assert U <= 4*visible_length <= 4*skeleton
        R = translation[1]
        n = ceil_sqrt(R)
        distinct_labels = {(a, 0, 1) for a, b in grid} | {(0, b, 0) for a, b in grid}
        assert len(distinct_labels) == 2*(3*K+1) < len(G)
        assert max(abs(x) for g in G for x in g) <= n
        assert max(abs(x) for w in planes for x in w) <= n*n
        assert len(G) > n
        counts['complete_primitive_grid_zonotopes'] += 1
        counts['independent_positive_facet_line_buckets'] += len(positive)
        counts['selected_positive_high_normal_facets'] += len(selected)
        fixtures.append({'grid_scale': K, 'matroid_rank': len(G),
                         'ground_elements': 2*len(G), 'distinct_exponent_labels': len(distinct_labels),
                         'hyperbolic_label_parameter': 9*K+2,
                         'profile_box_width': R, 'geometric_n_not_colength': n,
                         'unoriented_facet_planes': len(planes), 'positive_facets': len(positive),
                         'selected_high_normal_facets': len(selected),
                         'square_root_positive_incidence': U,
                         'total_l1_edge_length': skeleton, 'visible_l1_edge_length': visible_length,
                         'rank_squared_times_H': len(G)**2*(9*K+2),
                         'exact_positive_facet_witnesses': witnesses,
                         'all_positive_normals_sha256': hashlib.sha256(
                             json.dumps(sorted(positive)).encode()).hexdigest()})
    return {'counts': dict(counts), 'fixtures': fixtures,
            'scope': 'Classical partition-matroid projections of rank Theta(K^2), not length-n polynomial ideals. The geometric n=Theta(K^(3/2)) is smaller than that rank. Repeated exponent labels on different ground elements also prevent a monomial-column realization at the larger rank. No counterexample to the actual cubic budget is claimed.'}


def main():
    if not __debug__:
        raise RuntimeError('Exact assertions must be enabled.')
    models = [quartic_model(k) for k in range(2, 9)]
    old = json.loads((ROOT / 'results/research/extremal_initial_ideals_checks.json').read_text())
    rectangles = []
    for row in old['rectangular_complete_intersections']:
        r, s = row['rectangle']
        coefficients = {(a, b): c for a, b, c in row['polynomial_terms']}
        ground, columns, residue = rectangle_model(r, s, coefficients, row['characteristic'])
        rectangles.append(fan_incidence(r * s, ground, columns, residue, row['witnesses'], row['characteristic']))
    paths = ['scripts/run_weighted_facet_incidence_checks.py', 'docs/weighted_facet_incidence.md',
             'scripts/run_initial_cone_height_checks.py', 'scripts/run_extremal_initial_ideals_checks.py',
             'scripts/run_structural_complexity_checks.py', 'scripts/run_planar_dyadic_bound_checks.py',
             'scripts/run_interleaving_research.py', 'results/research/extremal_initial_ideals_checks.json',
             'scripts/run_adjacent_rank_checks.py', 'scripts/run_oracle_enumeration_checks.py']
    boxes = boxed_matroid_skeleton_controls()
    zonotopes = grid_zonotope_controls()
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_weighted_facet_incidence_checks.py',
              'scope': 'Uniform area obstruction and exact boundary-cover controls. No improved weighted incidence bound, sharp extremal count or novelty certificate is inferred.',
              'summary': {'quartic_fixtures': len(models), 'ordinary_s_pairs': sum(m['s_pairs'] for m in models),
                          'ground_normal_forms': sum(m['ground_normal_forms'] for m in models),
                          'complete_fans': len(models) + len(rectangles),
                          'boundary_cells_checked': sum(m['fan']['boundary_cells_checked'] for m in models) +
                                                    sum(m['boundary_cells_checked'] for m in rectangles)},
              'quartic_area_family': models, 'existing_rectangle_controls': rectangles,
              'boxed_matroid_skeletons': boxes, 'rank_sensitive_grid_zonotopes': zonotopes,
              'square_root_incidence_scope': 'U<=4 times visible l1 skeleton length is proved. The O(r^2*a*b*c) skeleton bound applies to labels in one box. A classical zonotope control violates a cubic estimate in the geometric profile-box parameter but has larger matroid rank, and is not an actual ideal counterexample. No global cubic budget, new unrestricted exponent or novelty is established.',
              'source_sha256': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}}
    (ROOT / 'results/research/weighted_facet_incidence_checks.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'old_summary': report['summary'], 'boxed_skeletons': boxes['counts'],
                      'grid_zonotopes': zonotopes['counts']}, indent=2))


if __name__ == '__main__':
    main()
