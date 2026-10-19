#!/usr/bin/env python3
"""Exact controls for a five-cone quotient with quadratic weight height.

The all-field proofs are in docs/initial_cone_height.md. These finite checks
do not prove novelty or improve the extremal initial-ideal count.
"""
from fractions import Fraction
from itertools import combinations, product
from math import gcd, lcm
from pathlib import Path
import hashlib
import json

from run_extremal_initial_ideals_checks import cone_for_basis
from run_structural_complexity_checks import (
    area2, clip, determinant, dot, greedy, hyperbolic_cross,
)

ROOT = Path(__file__).resolve().parents[1]


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def generators(k):
    return [
        {(1, 1, 0): 1, (0, 2, 0): -1},
        {(0, k + 1, 0): 1, (k, 0, 0): -1},
        {(k + 1, 0, 0): 1, (k, 0, 0): -1},
        {(0, 0, 1): 1, (k, 0, 0): -1},
    ]


def reference_basis(k):
    return tuple((a, 0, 0) for a in range(k + 1)) + tuple((0, b, 0) for b in range(1, k + 1))


def class_representative(exponent, k):
    a, b, c = exponent
    if c or (b == 0 and a >= k) or (b and a + b >= k + 1):
        return (k, 0, 0)
    return (0, a + b, 0) if b else (a, 0, 0)


def coefficient(value, characteristic):
    return value % characteristic if characteristic else Fraction(value)


def clean(poly, characteristic):
    return {a: coefficient(c, characteristic) for a, c in poly.items()
            if coefficient(c, characteristic)}


def leading(poly, weight):
    return max(poly, key=lambda a: (dot(weight, a), a))


def subtract_shift(poly, other, shift, factor, characteristic):
    result = dict(poly)
    for a, c in other.items():
        u = add(a, shift)
        result[u] = result.get(u, 0) - factor * c
    return clean(result, characteristic)


def normal_form(poly, basis, weight, characteristic):
    pending, remainder = clean(poly, characteristic), {}
    marked = [(leading(g, weight), g) for g in basis]
    steps = 0
    while pending:
        steps += 1
        assert steps < 100000, 'strict positive weighted division must terminate'
        u = leading(pending, weight)
        for v, g in marked:
            if all(a >= b for a, b in zip(u, v)):
                inverse = pow(g[v], -1, characteristic) if characteristic else 1 / g[v]
                pending = subtract_shift(pending, g, tuple(a - b for a, b in zip(u, v)),
                                         pending[u] * inverse, characteristic)
                break
        else:
            remainder[u] = pending.pop(u)
    return remainder


def polynomial_controls(k):
    weight = (k + 2, k + 1, (k + 1)**2)
    expected_heads = [(1, 1, 0), (0, k + 1, 0), (k + 1, 0, 0), (0, 0, 1)]
    cases = []
    for characteristic in (0, 2, 3, 5, 101):
        basis = [clean(g, characteristic) for g in generators(k)]
        heads = [leading(g, weight) for g in basis]
        assert heads == expected_heads
        for i, g in enumerate(basis):
            assert g[heads[i]] == 1
            assert all(not all(a >= b for a, b in zip(u, v))
                       for u in g for j, v in enumerate(heads) if j != i or u != heads[i])
        pairs = 0
        for i, j in combinations(range(len(basis)), 2):
            a, b = heads[i], heads[j]
            common = tuple(max(x, y) for x, y in zip(a, b))
            first = {add(u, tuple(x - y for x, y in zip(common, a))): c
                     for u, c in basis[i].items()}
            s_pair = subtract_shift(first, basis[j], tuple(x - y for x, y in zip(common, b)),
                                    1, characteristic)
            assert not normal_form(s_pair, basis, weight, characteristic)
            pairs += 1
        assert all(normal_form({a: 1}, basis, weight, characteristic) == {a: 1}
                   for a in reference_basis(k))
        cases.append({'characteristic': characteristic, 's_pairs': pairs})

    # Ordinary rational polynomial division is independent of the residue formula.
    basis = [clean(g, 0) for g in generators(k)]
    ground = tuple(hyperbolic_cross(2 * k + 1, 3))
    for a in ground:
        assert normal_form({a: 1}, basis, weight, 0) == {class_representative(a, k): 1}
    return {'k': k, 'length': 2 * k + 1, 'fields': cases, 'ground_normal_forms': len(ground)}


def quotient_model(k):
    reference = reference_basis(k)
    indices = {a: i for i, a in enumerate(reference)}
    maps = [tuple(indices[class_representative(add(a, tuple(int(i == j) for i in range(3))), k)]
                  for a in reference) for j in range(3)]
    for i, j in combinations(range(3), 2):
        assert all(maps[i][maps[j][v]] == maps[j][maps[i][v]] for v in range(len(reference)))

    def image(a, start=0):
        for i, power in enumerate(a):
            for _ in range(power):
                start = maps[i][start]
        return start

    assert all(image(a) == i for i, a in enumerate(reference)), 'the module is cyclic'
    for g in generators(k):
        for start in range(len(reference)):
            result = [0] * len(reference)
            for a, c in g.items():
                result[image(a, start)] += c
            assert not any(result), 'relations annihilate the entire integer module'

    def residue(a):
        index = image(a)
        assert reference[index] == class_representative(a, k)
        return tuple(int(i == index) for i in range(len(reference)))

    return reference, maps, residue


def integer_ray(point):
    scale = lcm(*(x.denominator for x in point))
    vector = tuple(int(x * scale) for x in point)
    factor = gcd(*vector)
    return tuple(x // factor for x in vector)


def fan_controls(k):
    reference, maps, residue = quotient_model(k)
    fixed = {(a, 0, 0) for a in range(k)}
    ys = {(0, j, 0) for j in range(1, k + 1)}
    mixed = {(j - 1, 1, 0) for j in range(1, k + 1)}
    specifications = [
        (ys, (k, 0, 0), (k + 2, k + 1, (k + 1)**2)),
        (ys, (0, k + 1, 0), (k + 2, k, k * (k + 1) + 1)),
        (ys, (0, 0, 1), (2, 1, 1)),
        (mixed, (k, 0, 0), (1, 2, k + 1)),
        (mixed, (0, 0, 1), (1, 2, 1)),
    ]
    ground = tuple(hyperbolic_cross(2 * k + 1, 3))
    columns = {a: residue(a) for a in ground}
    states, polygons = [], []
    for middle, last, weight in specifications:
        expected = fixed | middle | {last}
        selected = greedy(ground, columns, weight, 101)
        assert selected == expected and len(selected) == 2 * k + 1
        polygon, normals = cone_for_basis(selected, residue, 101)
        assert all(dot(z, weight) > 0 for z in normals)
        assert all(tuple(x - int(i == j) for i, x in enumerate(a)) in selected
                   for a in selected for j in range(3) if a[j])
        rays = [integer_ray(p) for p in polygon]
        assert all(max(ray) <= 2 * (2 * k + 1)**2 for ray in rays)
        independent = next(triple for triple in combinations(rays, 3) if determinant(triple))
        certificate = tuple(sum(a[i] for a in independent) for i in range(3))
        assert min(certificate) > 0 and all(dot(z, certificate) > 0 for z in normals)
        assert max(certificate) <= 6 * (2 * k + 1)**2
        states.append({'standard_set': sorted(selected), 'strict_weight': weight,
                       'normals': normals, 'primitive_rays': rays,
                       'integer_certificate': certificate,
                       'simplex_area_fraction': str(area2(polygon)),
                       'polygon': [[str(x) for x in p] for p in polygon]})
        polygons.append((polygon, normals))
    assert len({tuple(map(tuple, s['standard_set'])) for s in states}) == 5
    assert sum(area2(p) for p, _ in polygons) == 1
    for (polygon, _), (_, normals) in combinations(polygons, 2):
        overlap = polygon
        for z in normals:
            overlap = clip(overlap, z, 1)
        assert not overlap, 'different cone interiors are disjoint'
    first_rays = {(k + 1, k, k * (k + 1)), (1, 1, k), (0, 0, 1)}
    assert set(map(tuple, states[0]['primitive_rays'])) == first_rays
    profile_exchanges = []
    original = set(reference)
    for removed, inserted in [((0, 2, 0), (1, 1, 0)),
                              ((k, 0, 0), (0, k + 1, 0)),
                              ((k, 0, 0), (0, 0, 1))]:
        exchanged = original - {removed} | {inserted}
        assert inserted in columns and len(exchanged) == len(reference)
        assert len({columns[a] for a in exchanged}) == len(reference)
        difference = tuple(b - a for a, b in zip(removed, inserted))
        profile_exchanges.append({'removed': removed, 'inserted': inserted,
                                  'profile_difference': difference})
    assert determinant([row['profile_difference'] for row in profile_exchanges]) == 1
    Q = k // 2 + 1
    covering_planes = [(Q, Q, -1), (1, 0, -1)]
    covering_weights = [
        (0, (2 * k + 1, 2 * k, Q * (4 * k + 1))),
        (0, (2 * k + 2, 1, Q * (2 * k + 3))),
        (1, (2, 1, 2)),
        (0, (1, 2 * k + 2, Q * (2 * k + 3))),
        (1, (1, 2, 1)),
    ]
    for state, (plane, weight) in zip(states, covering_weights):
        assert dot(covering_planes[plane], weight) == 0
        assert min(weight) > 0 and all(dot(z, weight) > 0 for z in state['normals'])
        assert greedy(ground, columns, weight, 101) == set(map(tuple, state['standard_set']))
        state['cover_plane'] = plane
        state['cover_objective'] = weight
    line_hit_histogram = {}
    for q in product(range(-3, 4), repeat=3):
        if not any(q) or gcd(*q) != 1 or next(x for x in q if x) < 0:
            continue
        hits = 0
        for state in states:
            values = [dot(q, r) for r in state['primitive_rays']]
            hits += int(min(values) < 0 < max(values))
        assert hits <= 4
        line_hit_histogram[hits] = line_hit_histogram.get(hits, 0) + 1
    return {'k': k, 'length': len(reference), 'initial_ideals': 5,
            'multiplication_maps': maps, 'full_simplex_area': '1',
            'independent_profile_exchanges': profile_exchanges,
            'two_plane_cover': covering_planes,
            'small_plane_strict_cone_hit_histogram': line_hit_histogram, 'states': states}


def arithmetic_controls(k):
    V = ((1, -1, 0), (-k, k + 1, 0), (-k, 0, 1))
    rays = ((k + 1, k, k * (k + 1)), (1, 1, k), (0, 0, 1))
    assert determinant(V) == determinant(rays) == 1
    assert [[dot(v, r) for r in rays] for v in V] == [[int(i == j) for j in range(3)] for i in range(3)]
    weight = tuple(sum(r[i] for r in rays) for i in range(3))
    assert weight == (k + 2, k + 1, (k + 1)**2)
    assert tuple(dot(v, weight) for v in V) == (1, 1, 1)
    Q = k // 2 + 1
    q = (Q, Q, -1)
    visible = (2 * k + 1, 2 * k, Q * (4 * k + 1))
    assert dot(q, visible) == 0 and all(dot(v, visible) > 0 for v in V)
    tested = 0
    for candidate in product(range(-Q, Q + 1), repeat=3):
        if not any(candidate) or gcd(*candidate) != 1 or next(x for x in candidate if x) < 0:
            continue
        values = [dot(candidate, r) for r in rays]
        meets_interior = min(values) < 0 < max(values)
        if max(map(abs, candidate)) < Q:
            assert not meets_interior, 'independent exhaustive subcritical kernel check'
        tested += 1
    integer_tests = 0
    if k <= 6:
        bound = (k + 1)**2 - 1
        for candidate in product(range(1, bound + 1), repeat=3):
            assert not all(dot(v, candidate) > 0 for v in V)
            integer_tests += 1
    return {'k': k, 'length': 2 * k + 1, 'slack_matrix': V, 'primitive_rays': rays,
            'minimal_integer_weight': weight, 'minimal_integer_height': max(weight),
            'minimal_kernel_radius': Q, 'kernel_witness': q, 'visible_objective': visible,
            'primitive_kernels_examined': tested, 'subcritical_integer_weights_examined': integer_tests}


def main():
    polynomial = [polynomial_controls(k) for k in range(2, 21)]
    fans = [fan_controls(k) for k in range(2, 13)]
    arithmetic = [arithmetic_controls(k) for k in range(2, 25)]
    paths = ['scripts/run_initial_cone_height_checks.py', 'docs/initial_cone_height.md',
             'scripts/run_extremal_initial_ideals_checks.py', 'scripts/run_structural_complexity_checks.py',
             'scripts/run_interleaving_research.py']
    report = {
        'status': 'passed',
        'command': 'python3 -B scripts/run_initial_cone_height_checks.py',
        'scope': 'Finite exact supporting controls for the all-field uniform proof. No extremal fan-count improvement, novelty certificate or goal completion is inferred.',
        'summary': {'polynomial_fixtures': len(polynomial),
                    'ordinary_s_pairs': sum(t['s_pairs'] for row in polynomial for t in row['fields']),
                    'ground_normal_forms': sum(row['ground_normal_forms'] for row in polynomial),
                    'complete_fans': len(fans), 'initial_ideals_checked': sum(row['initial_ideals'] for row in fans),
                    'arithmetic_fixtures': len(arithmetic),
                    'primitive_kernels_examined': sum(row['primitive_kernels_examined'] for row in arithmetic),
                    'subcritical_integer_weights_examined': sum(row['subcritical_integer_weights_examined'] for row in arithmetic)},
        'polynomial_controls': polynomial, 'complete_fans': fans, 'arithmetic_controls': arithmetic,
        'source_sha256': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
    }
    target = ROOT / 'results/research/initial_cone_height_checks.json'
    target.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report['summary'], indent=2))
    print('Saved results/research/initial_cone_height_checks.json')


if __name__ == '__main__':
    main()
