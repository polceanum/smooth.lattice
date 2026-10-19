#!/usr/bin/env python3
"""Exact controls for weighted facet incidence and an area-budget obstruction."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import gcd
from pathlib import Path
import hashlib
import json

from run_extremal_initial_ideals_checks import cone_for_basis, rectangle_model
from run_initial_cone_height_checks import (
    add, clean, integer_ray, leading, normal_form, subtract_shift,
)
from run_planar_dyadic_bound_checks import convex_hull
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


def main():
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
             'scripts/run_interleaving_research.py', 'results/research/extremal_initial_ideals_checks.json']
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_weighted_facet_incidence_checks.py',
              'scope': 'Uniform area obstruction and exact boundary-cover controls. No improved weighted incidence bound, sharp extremal count or novelty certificate is inferred.',
              'summary': {'quartic_fixtures': len(models), 'ordinary_s_pairs': sum(m['s_pairs'] for m in models),
                          'ground_normal_forms': sum(m['ground_normal_forms'] for m in models),
                          'complete_fans': len(models) + len(rectangles),
                          'boundary_cells_checked': sum(m['fan']['boundary_cells_checked'] for m in models) +
                                                    sum(m['boundary_cells_checked'] for m in rectangles)},
              'quartic_area_family': models, 'existing_rectangle_controls': rectangles,
              'source_sha256': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}}
    (ROOT / 'results/research/weighted_facet_incidence_checks.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report['summary'], indent=2))


if __name__ == '__main__':
    main()
