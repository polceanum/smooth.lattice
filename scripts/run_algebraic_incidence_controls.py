#!/usr/bin/env python3
"""Exact negative controls for geometric approaches to quotient incidence.

The uniform-axis matroid is not a quotient residue matroid. Different
three-pole quotient matroids realize exactly its positive fan, disproving
the proposed cubic incidence budget. Nine-cone quotients check a second
shortcut. No stronger extremal count or novelty certificate is inferred. See the two linked notes.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import gcd
from pathlib import Path
import hashlib
import json

from run_extremal_initial_ideals_checks import cone_for_basis, integer_weight
from run_initial_cone_height_checks import (
    add, class_representative, clean, generators, integer_ray, leading,
    normal_form, reference_basis, subtract_shift,
)
from run_structural_complexity_checks import (
    area2, clip, cyclic_univariate_model, determinant, dot, greedy,
    hyperbolic_cross, univariate_power_mod, univariate_product_mod,
)
from run_interleaving_research import minimal_excluded_border

ROOT = Path(__file__).resolve().parents[1]


def axis_controls(n):
    r = n - 1
    axes = tuple(tuple(j if i == axis else 0 for i in range(3))
                 for axis in range(3) for j in range(1, r + 1))
    ground = ((0, 0, 0),) + axes
    assert all((a[0] + 1) * (a[1] + 1) * (a[2] + 1) <= n for a in ground)
    incidence, total_area, selections, samples = defaultdict(set), Fraction(0), 0, []
    # Distinct Vandermonde parameters represent U_(r,3r), plus a coloop.
    columns = {(0, 0, 0): (1,) + (0,) * r}
    columns.update({a: (0,) + tuple(pow(t, j, 101) for j in range(r))
                    for t, a in enumerate(axes, 1)})
    assert 3 * r < 101
    for a in range(r + 1):
        for b in range(r - a + 1):
            allocation = (a, b, r - a - b)
            basis = {(0, 0, 0)} | {point for point in axes
                      if next(x for x in point if x) <= allocation[next(i for i, x in enumerate(point) if x)]}
            w = integer_weight(tuple(Fraction(2, 2 * count + 1) for count in allocation))
            selected = {(0, 0, 0)} | set(sorted(axes, key=lambda p: (dot(w, p), p))[:r])
            assert basis == selected and len(basis) == n
            assert all(tuple(x - int(i == j) for i, x in enumerate(p)) in basis
                       for p in basis for j in range(3) if p[j])
            if n <= 7:
                assert greedy(ground, columns, w, 101) == basis
                selections += 1
            polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
            for i, j in product(range(3), repeat=2):
                if allocation[i] < r and allocation[j] and i != j:
                    normal = tuple((allocation[i] + 1 if k == i else 0) -
                                   (allocation[j] if k == j else 0) for k in range(3))
                    assert dot(normal, w) > 0
                    polygon = clip(polygon, normal, 1)
            assert polygon and area2(polygon) > 0
            total_area += area2(polygon)
            index = len(samples)
            samples.append((allocation, w))
            for p in polygon:
                ray = integer_ray(p)
                if min(ray) > 0:
                    incidence[ray].add(index)
    assert total_area == 1
    expected = {}
    for level in (n, n + 1):
        for A in range(1, level - 1):
            for B in range(1, level - A):
                C = level - A - B
                vector = (B * C, A * C, A * B)
                factor = gcd(*vector)
                ray = tuple(x // factor for x in vector)
                assert ray not in expected
                expected[ray] = (A, B, C)
    assert set(incidence) == set(expected)
    assert all(len(cells) == 3 for cells in incidence.values())
    T = sum(max(w) * len(cells) for w, cells in incidence.items())
    assert len(samples) == n * (n + 1) // 2 and len(incidence) == (n - 1)**2
    maximum = max(incidence, key=max)
    return {'rank': n, 'nonloop_ground_size': len(ground), 'positive_cells': len(samples),
            'all_positive_bases_downsets': True, 'complete_simplex_double_area': str(total_area),
            'interior_rays': len(incidence), 'all_interior_degrees': 3,
            'weighted_interior_incidence': T, 'maximum_ray_height': max(maximum),
            'maximum_ray': maximum, 'threshold_triple': expected[maximum],
            'independent_matrix_selections': selections,
            'quotient_realizable': False,
            'obstruction': 'xy is a loop, y is in the span of x through x^(n-1), but y^2 is a nonloop.',
            'allocation_witnesses': [{'counts': a, 'strict_weight': w} for a, w in samples[:3]]}


def arithmetic_controls(n):
    assert n % 6 == 3
    triples, coprime, contribution = 0, 0, 0
    for a in range(1, n, 6):
        if not n <= 4 * a <= 2 * n:
            continue
        for b in range(1, n - a, 6):
            c = n - a - b
            if c % 6 != 1 or not n <= 4 * b <= 2 * n or not n <= 4 * c <= 2 * n:
                continue
            triples += 1
            if gcd(a, b) == gcd(a, c) == gcd(b, c) == 1:
                coprime += 1
                normal = (b * c, a * c, a * b)
                assert gcd(*normal) == 1 and 16 * max(normal) >= n * n
                contribution += 3 * max(normal)
    assert coprime > 0
    return {'rank': n, 'central_congruence_triples': triples,
            'pairwise_coprime_triples': coprime,
            'exact_partial_weighted_incidence': contribution,
            'scope': 'Finite arithmetic control of the uniform sieve proof, not an asymptotic fit.'}


def pair_checks(polynomials, weight, characteristic):
    basis = [clean(g, characteristic) for g in polynomials]
    heads = [leading(g, weight) for g in basis]
    for i, j in combinations(range(len(basis)), 2):
        common = tuple(max(a, b) for a, b in zip(heads[i], heads[j]))
        first = {add(a, tuple(x - y for x, y in zip(common, heads[i]))): c
                 for a, c in basis[i].items()}
        pair = subtract_shift(first, basis[j], tuple(x - y for x, y in zip(common, heads[j])),
                              1, characteristic)
        assert not normal_form(pair, basis, weight, characteristic)
    return len(basis) * (len(basis) - 1) // 2


def interior_controls(k):
    reference = reference_basis(k)
    index = {a: i for i, a in enumerate(reference)}

    def reduced(a):
        return class_representative((a[0], a[1] + k * a[2], 0), k)

    maps = [tuple(index[reduced(add(a, tuple(int(i == j) for i in range(3))))]
                  for a in reference) for j in range(3)]
    for i, j in combinations(range(3), 2):
        assert all(maps[i][maps[j][v]] == maps[j][maps[i][v]] for v in range(len(reference)))

    def image(a, start=0):
        for i, power in enumerate(a):
            for _ in range(power):
                start = maps[i][start]
        return start

    assert all(image(a) == i for i, a in enumerate(reference))
    polynomials = generators(k)[:3] + [{(0, 0, 1): 1, (0, k, 0): -1}]
    for g in polynomials:
        for start in range(len(reference)):
            vector = [0] * len(reference)
            for a, c in g.items():
                vector[image(a, start)] += c
            assert not any(vector)

    def residue(a):
        target = image(a)
        assert reference[target] == reduced(a)
        return tuple(int(i == target) for i in range(len(reference)))

    hard_weight = (k + 3, k + 2, k * k + 2 * k - 1)
    hard = [{(1, 1, 0): 1, (0, 2, 0): -1},
            {(0, k, 0): 1, (0, 0, 1): -1},
            {(k + 1, 0, 0): 1, (k, 0, 0): -1},
            {(1, 0, 1): 1, (k, 0, 0): -1},
            {(0, 1, 1): 1, (k, 0, 0): -1},
            {(0, 0, 2): 1, (k, 0, 0): -1}]
    reference_weight = (k + 2, k + 1, k * (k + 1) + 1)
    ground = tuple(hyperbolic_cross(2 * k + 1, 3))
    columns = {a: residue(a) for a in ground}
    s_pairs, reductions = 0, 0
    for characteristic in (0, 2, 3, 5, 101):
        s_pairs += pair_checks(polynomials, reference_weight, characteristic)
        s_pairs += pair_checks(hard, hard_weight, characteristic)
        for a in ground:
            expected = reduced(a)
            assert normal_form({a: 1}, [clean(g, characteristic) for g in polynomials],
                               reference_weight, characteristic) == {expected: 1}
            alternative = (0, 0, 1) if expected == (0, k, 0) else expected
            assert normal_form({a: 1}, [clean(g, characteristic) for g in hard],
                               hard_weight, characteristic) == {alternative: 1}
            reductions += 2
    weights = [(k + 2, k + 1, k * (k + 1) + 1),
               (k + 2, k, k * k + k + 1), hard_weight, (k + 2, 1, 2), (3, 2, 1),
               (1, 2, 3 * k), (2, 3, 2 * k - 1), (2, 3, 3), (2, 3, 1)]
    states, polygons = [], []
    for w in weights:
        basis = greedy(ground, columns, w, 101)
        polygon, normals = cone_for_basis(basis, residue, 101)
        assert all(dot(z, w) > 0 for z in normals)
        states.append({'basis': sorted(basis), 'strict_weight': w,
                       'primitive_rays': [integer_ray(p) for p in polygon],
                       'simplex_double_area': str(area2(polygon))})
        polygons.append((polygon, normals))
    assert len({tuple(map(tuple, s['basis'])) for s in states}) == 9
    assert sum(area2(p) for p, _ in polygons) == 1
    for (polygon, _), (_, normals) in combinations(polygons, 2):
        overlap = polygon
        for z in normals:
            overlap = clip(overlap, z, 1)
        assert not overlap
    hard_rays = [(k + 1, k, k * k), (1, 1, k - 1), (1, 1, k)]
    assert set(map(tuple, states[2]['primitive_rays'])) == set(hard_rays)
    assert [i for i, s in enumerate(states) if all(min(p) > 0 for p in s['primitive_rays'])] == [2]
    slack = [(1, -1, 0), (0, k, -1), (-k, 1, 1)]
    assert determinant(slack) == 1
    assert [[dot(z, ray) for ray in hard_rays] for z in slack] == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    assert tuple(sum(ray[i] for ray in hard_rays) for i in range(3)) == hard_weight
    q = (k, k - 1, -2)
    assert gcd(*q) == 1 and dot(q, hard_weight) == 0
    plane_checks = 0
    for candidate in product(range(-k, k + 1), repeat=3):
        if not any(candidate) or gcd(*candidate) != 1 or next(x for x in candidate if x) < 0:
            continue
        values = [dot(candidate, ray) for ray in hard_rays]
        if max(map(abs, candidate)) < k:
            assert not min(values) < 0 < max(values)
        plane_checks += 1
    height, weight_checks = max(hard_weight), 0
    for s in range(1, height + 1):
        for t in range(1, s):
            u = k * s - t + 1
            if u < k * t:
                assert max(s, t, u) >= height
                weight_checks += 1
    return {'k': k, 'length': len(reference), 'initial_ideals': 9,
            'unique_interior_cell': 2, 'minimum_integer_weight': hard_weight,
            'minimum_integer_height': height, 'minimum_kernel_height': k,
            'kernel_witness': q, 'multiplication_maps': maps, 'states': states,
            'ordinary_s_pairs': s_pairs, 'ordinary_normal_forms': reductions,
            'primitive_plane_controls': plane_checks, 'integer_weight_controls': weight_checks}


def field_inverse(basis, columns, characteristic):
    n = len(basis)
    scalar = (lambda x: x % characteristic) if characteristic else Fraction
    rows = [[scalar(columns[a][i]) for a in basis] + [scalar(i == j) for j in range(n)]
            for i in range(n)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scale = pow(rows[j][j], -1, characteristic) if characteristic else 1 / rows[j][j]
        rows[j] = [scalar(x * scale) for x in rows[j]]
        for i in range(n):
            if i != j:
                scale = rows[i][j]
                rows[i] = [scalar(x - scale * y) for x, y in zip(rows[i], rows[j])]
    return [row[n:] for row in rows]


def three_pole_controls(n, characteristic, points=False):
    ground = tuple(hyperbolic_cross(n, 3))
    unit = (1,) + (0,) * (n - 1)
    if points:
        nodes = tuple(range(2, n + 2))
        assert not characteristic or characteristic > n + 1

        def residue(a):
            if characteristic:
                return tuple(pow(t, a[0], characteristic) * pow(pow(t, -1, characteristic), a[1], characteristic) *
                             pow(pow(t - 1, -1, characteristic), a[2], characteristic) % characteristic for t in nodes)
            return tuple(Fraction(t)**(a[0] - a[1]) * Fraction(t - 1)**(-a[2]) for t in nodes)

        columns = {a: residue(a) for a in ground}
        assert determinant([[t**j for j in range(n)] for t in nodes]) != 0
    else:
        modulus = [1, -1] + [0] * (n - 2) + [1]
        images = [(0, 1) + (0,) * (n - 2), (1,) + (0,) * (n - 2) + (-1,),
                  (0,) + (-1,) * (n - 1)]
        columns = cyclic_univariate_model(modulus, images, ground)
        assert univariate_product_mod(images[0], images[1], modulus) == unit
        difference = tuple(a - b for a, b in zip(images[0], unit))
        assert univariate_product_mod(difference, images[2], modulus) == unit

        def residue(a):
            vector = unit
            for image, exponent in zip(images, a):
                vector = univariate_product_mod(vector, univariate_power_mod(image, exponent, modulus), modulus)
            return vector

        assert all(columns[a] == residue(a) for a in ground)
    scalar = (lambda x: x % characteristic) if characteristic else Fraction
    incidence, total_area, pairs, reductions, states, witness = defaultdict(set), Fraction(0), 0, 0, 0, None
    for a in range(n):
        for b in range(n - a):
            counts = (a, b, n - 1 - a - b)
            selected = {(0, 0, 0)} | {tuple(j if i == axis else 0 for i in range(3))
                        for axis, count in enumerate(counts) for j in range(1, count + 1)}
            w = integer_weight(tuple(Fraction(2, 2 * count + 1) for count in counts))
            assert greedy(ground, columns, w, characteristic) == selected
            ordered = sorted(selected)
            inverse = field_inverse(ordered, columns, characteristic)

            def coefficients(u):
                vector = residue(u)
                return [scalar(sum(x * y for x, y in zip(row, vector))) for row in inverse]

            basis, normals = [], set()
            for u in minimal_excluded_border(ordered):
                polynomial = {u: scalar(1)}
                for v, c in zip(ordered, coefficients(u)):
                    if c:
                        polynomial[v] = -c
                        normal = tuple(x - y for x, y in zip(u, v))
                        factor = gcd(*normal)
                        normals.add(tuple(x // factor for x in normal))
                assert leading(polynomial, w) == u
                assert all(dot(w, u) > dot(w, v) for v in polynomial if v != u)
                basis.append(polynomial)
            polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
            for z in sorted(normals):
                polygon = clip(polygon, z, 1)
            assert polygon and area2(polygon) > 0
            total_area += area2(polygon)
            pairs += pair_checks(basis, w, characteristic)
            for u in ground:
                expected = {v: c for v, c in zip(ordered, coefficients(u)) if c}
                assert normal_form({u: 1}, basis, w, characteristic) == expected
                reductions += 1
            # Independent axis-threshold cone, without polynomial reductions.
            threshold_polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
            for i, j in product(range(3), repeat=2):
                if counts[i] < n - 1 and counts[j] and i != j:
                    z = tuple((counts[i] + 1 if k == i else 0) - (counts[j] if k == j else 0) for k in range(3))
                    threshold_polygon = clip(threshold_polygon, z, 1)
            rays = {integer_ray(p) for p in polygon}
            assert rays == {integer_ray(p) for p in threshold_polygon}
            for ray in rays:
                if min(ray) > 0:
                    incidence[ray].add(states)
            states += 1
            if min(counts) > 0 and witness is None:
                witness = {'axis_counts': counts, 'strict_weight': w, 'basis': ordered,
                           'primitive_rays': sorted(rays),
                           'reduced_basis': [[[u, str(c)] for u, c in sorted(g.items())] for g in basis]}
    assert total_area == 1 and states == n * (n + 1) // 2
    assert len(incidence) == (n - 1)**2 and all(len(v) == 3 for v in incidence.values())
    T = sum(max(ray) * len(cells) for ray, cells in incidence.items())
    expected_T = 0
    for level in (n, n + 1):
        for A in range(1, level - 1):
            for B in range(1, level - A):
                C = level - A - B
                vector = (B * C, A * C, A * B)
                expected_T += 3 * max(vector) // gcd(*vector)
    assert T == expected_T
    return {'length': n, 'characteristic': characteristic,
            'model': 'distinct_rational_points' if points else 't^n-t+1_quotient',
            'initial_ideals': states, 'interior_rays': len(incidence),
            'weighted_interior_incidence': T, 'complete_simplex_double_area': str(total_area),
            'all_cones_match_independent_axis_thresholds': True,
            'ordinary_s_pairs': pairs, 'ordinary_normal_forms': reductions,
            'interior_cell_witness': witness}


def axis_plane_controls(n, start, end):
    r = n - 1
    axes = tuple(tuple(j if i == axis else 0 for i in range(3))
                 for axis in range(3) for j in range(1, r + 1))
    assert all(a + b > 0 for a, b in zip(start, end))
    assert min(start) == min(end) == 0
    rates = [Fraction(v, u + v) for u, v in zip(start, end)]
    assert len(set(rates)) == 3
    rank = {i: sorted(range(3), key=lambda i: rates[i]).index(i) for i in range(3)}
    parameters = {Fraction(0)}
    for a, b in combinations(axes, 2):
        difference = tuple(x - y for x, y in zip(a, b))
        initial, speed = dot(difference, start), dot(difference, end)
        if speed and -Fraction(initial, speed) > 0:
            parameters.add(-Fraction(initial, speed))
    parameters = sorted(parameters)
    samples = [(a + b) / 2 for a, b in zip(parameters, parameters[1:])]
    samples.append(parameters[-1] + 1)
    states, potentials = [], []
    for t in samples:
        weight = tuple(u + t * v for u, v in zip(start, end))
        selected = sorted(axes, key=lambda a: dot(a, weight))[:r]
        counts = tuple(sum(a[i] > 0 for a in selected) for i in range(3))
        potential = sum(rank[i] * count for i, count in enumerate(counts))
        if not states or states[-1] != counts:
            assert not potentials or potential < potentials[-1]
            states.append(counts)
            potentials.append(potential)
    assert len(states) <= 2 * n - 1
    kernel = (start[1] * end[2] - start[2] * end[1],
              start[2] * end[0] - start[0] * end[2],
              start[0] * end[1] - start[1] * end[0])
    factor = gcd(*kernel)
    kernel = tuple(x // factor for x in kernel)
    assert dot(kernel, start) == dot(kernel, end) == 0
    return {'length': n, 'kernel': kernel, 'boundary_rays': [start, end],
            'complete_comparison_intervals': len(samples), 'strict_cells_met': len(states),
            'uniform_upper_bound': 2 * n - 1, 'allocation_path': states,
            'strictly_decreasing_type_potential': potentials,
            'selected_plane_cover_lower_bound': (n * (n + 1) // 2 + 2 * n - 2) // (2 * n - 1)}


def main():
    axes = [axis_controls(n) for n in (4, 5, 6, 7, 9, 12, 15, 21)]
    arithmetic = [arithmetic_controls(n) for n in (99, 303, 603)]
    interior = [interior_controls(k) for k in range(3, 11)]
    curves = [three_pole_controls(n, characteristic) for n in (4, 5, 6, 7)
              for characteristic in (0, 2, 3, 5, 101)]
    curves += [three_pole_controls(n, characteristic, points=True) for n in (5, 7)
               for characteristic in (0, 101)]
    planes = [axis_plane_controls(n, start, end) for n in (4, 5, 6, 7, 9, 12, 15, 21)
              for start, end in [((1, 1, 0), (1, 0, 1)), ((2, 1, 0), (0, 1, 2)),
                                 ((1, 0, 2), (0, 1, 3)), ((5, 3, 0), (7, 0, 3))]]
    paths = ['scripts/run_algebraic_incidence_controls.py', 'docs/weighted_facet_incidence.md',
             'docs/initial_cone_height.md', 'scripts/run_extremal_initial_ideals_checks.py',
             'scripts/run_structural_complexity_checks.py', 'scripts/run_initial_cone_height_checks.py',
             'scripts/run_interleaving_research.py']
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_algebraic_incidence_controls.py',
              'scope': 'Different three-pole quotient matroids realize the uniform-axis positive fan and give quartic incidence, disproving the cubic budget. Nine-cone quotients check interior radius. No stronger G_3 count or novelty certificate is inferred.',
              'summary': {'complete_axis_fans': len(axes), 'axis_cells': sum(a['positive_cells'] for a in axes),
                          'matrix_selection_controls': sum(a['independent_matrix_selections'] for a in axes),
                          'arithmetic_fixtures': len(arithmetic), 'complete_quotient_fans': len(interior),
                          'ordinary_s_pairs': sum(a['ordinary_s_pairs'] for a in interior),
                          'ordinary_normal_forms': sum(a['ordinary_normal_forms'] for a in interior),
                          'primitive_plane_controls': sum(a['primitive_plane_controls'] for a in interior),
                          'complete_three_pole_fans': len(curves),
                          'three_pole_initial_ideals': sum(a['initial_ideals'] for a in curves),
                          'three_pole_s_pairs': sum(a['ordinary_s_pairs'] for a in curves),
                          'three_pole_normal_forms': sum(a['ordinary_normal_forms'] for a in curves),
                          'complete_positive_plane_traces': len(planes),
                          'plane_comparison_intervals': sum(a['complete_comparison_intervals'] for a in planes)},
              'uniform_axis_matroids': axes, 'arithmetic_lower_bound_controls': arithmetic,
              'interior_quotients': interior, 'three_pole_quotients': curves,
              'positive_plane_traces': planes,
              'source_sha256': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}}
    (ROOT / 'results/research/algebraic_incidence_controls.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report['summary'], indent=2))


if __name__ == '__main__':
    main()
