#!/usr/bin/env python3
"""Exact negative controls for geometric approaches to quotient incidence.

The uniform-axis matroid is not a quotient residue matroid. Different
three-pole quotient matroids realize exactly its positive fan, disproving
the proposed cubic incidence budget. Nine-cone quotients check a second
shortcut. Scaled-component controls verify uniform carrier bases and
generic sampling ranks. No unrestricted extremal count or novelty
certificate is inferred. See the linked research notes.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import comb, gcd, isqrt
from heapq import nlargest
from pathlib import Path
from random import Random
import hashlib
import json

from run_extremal_initial_ideals_checks import (
    adjacent_weights, cone_for_basis, integer_weight, primitive,
)
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


def two_curve_symbolic_controls():
    """Universal certificates over Z[s,p], with parameter exponents last."""
    zero = (0, 0, 0, 0, 0)
    g = [
        {(1, 0, 1, 0, 0): 1, (0, 0, 1, 0, 0): -1, zero: -1},
        {(2, 2, 0, 0, 0): 1, (1, 1, 0, 1, 0): -1, (0, 0, 0, 0, 1): 1},
        {(1, 2, 0, 0, 0): 1, (0, 2, 1, 0, 0): 1, (0, 2, 0, 0, 0): 1,
         (0, 1, 1, 1, 0): -1, (0, 1, 0, 1, 0): -1, (0, 0, 1, 0, 1): 1},
        {(0, 2, 2, 0, 0): 1, (0, 2, 1, 0, 0): 2, (0, 2, 0, 0, 0): 1,
         (0, 1, 2, 1, 0): -1, (0, 1, 1, 1, 0): -1, (0, 0, 2, 0, 1): 1},
    ]

    def relation(terms):
        value = {}
        for sign, index, shift in terms:
            value = subtract_shift(value, g[index], shift, -sign, 0)
        assert not value

    relation([(1, 2, zero), (-1, 1, (0, 0, 1, 0, 0)),
              (1, 0, (1, 2, 0, 0, 0)), (1, 0, (0, 2, 0, 0, 0)),
              (-1, 0, (0, 1, 0, 1, 0))])
    relation([(1, 3, zero), (-1, 2, (0, 0, 1, 0, 0)),
              (1, 0, (0, 2, 0, 0, 0))])
    relation([(1, 1, zero), (-1, 2, (1, 0, 0, 0, 0)), (1, 2, zero),
              (1, 0, (0, 2, 0, 0, 0)), (-1, 0, (0, 1, 0, 1, 0)),
              (1, 0, (0, 0, 0, 0, 1))])
    cases = []
    for side, indices in [(1, (0, 2, 3)), (-1, (0, 1, 2))]:
        weight = (5, 3, 2, 0, 0) if side == 1 else (2, 3, 5, 0, 0)
        basis = [g[i] for i in indices]

        def universally_lower(a, b):
            dx, dy, dz = (x - y for x, y in zip(b[:3], a[:3]))
            coefficients = (dx, dy, dx + dz) if side == 1 else (dz, dy, dx + dz)
            return min(coefficients) >= 0 and max(coefficients) > 0

        heads = [leading(f, weight) for f in basis]
        assert all(h[3:] == (0, 0) and f[h] == 1 for h, f in zip(heads, basis))
        assert all(universally_lower(a, h) for h, f in zip(heads, basis) for a in f if a != h)
        certificates = []
        for i, j in combinations(range(3), 2):
            common = tuple(max(a, b) for a, b in zip(heads[i], heads[j]))
            first = {add(a, tuple(x - y for x, y in zip(common, heads[i]))): c
                     for a, c in basis[i].items()}
            pending = subtract_shift(first, basis[j], tuple(x - y for x, y in zip(common, heads[j])), 1, 0)
            steps = []
            while pending:
                a = leading(pending, weight)
                assert universally_lower(a, common)
                for index, (h, f) in enumerate(zip(heads, basis)):
                    if all(x >= y for x, y in zip(a, h)):
                        shift, coefficient = tuple(x - y for x, y in zip(a, h)), pending[a]
                        steps.append({'basis_index': indices[index], 'shift': shift,
                                      'coefficient': str(coefficient)})
                        pending = subtract_shift(pending, f, shift, coefficient, 0)
                        break
                else:
                    raise AssertionError('nonzero formal S-pair remainder')
            certificates.append({'pair': [indices[i], indices[j]], 'reductions': steps})
        cases.append({'domain': 'w_x>w_z>0' if side == 1 else 'w_z>w_x>0',
                      'basis_indices': indices, 'heads': [h[:3] for h in heads],
                      'formal_s_pair_certificates': certificates})
    return {'coefficient_ring': 'Z[s,p]', 'parameter_order': ['s', 'p'],
            'generators': [[[a, c] for a, c in sorted(f.items())] for f in g],
            'ideal_identity_checks': 3, 'domains': cases,
            'all_characteristics_and_all_positive_weights_in_each_domain': True}


def two_curve_cut_fan(n, side, quota):
    """Enumerate only this proved six-ray model, without quotient ranks."""
    rays = [((1, 0, 0), (1, 0, 0)), ((0, 1, 0), (0, 1, 0)),
            ((0, 0, 1), (0, 0, 1)), ((1, 1, 0), (1, 0, 0)),
            ((0, 2, 1), (0, 1, 0)), ((0, 1, 1), (0, 0, 1))]
    if side < 0:
        rays = [(tuple(reversed(a)), tuple(reversed(b))) for a, b in rays]
    ground = ((0, 0, 0),) + tuple(add(a, tuple(j * x for x in b)) for a, b in rays for j in range(n))
    comparisons = sorted({primitive(tuple(x - y for x, y in zip(a, b)))
                          for a, b in combinations(ground, 2)})
    todo, seen = [(5, 3, 2) if side > 0 else (2, 3, 5)], {}
    def select(weight):
        chosen, secondary = [], 0
        for a in sorted(ground, key=lambda a: (dot(a, weight), a)):
            mixed = sum(x > 0 for x in a) > 1
            if not mixed or secondary < quota:
                chosen.append(a)
                secondary += int(mixed)
            if len(chosen) == n:
                break
        return tuple(sorted(chosen)), secondary

    while todo:
        weight = todo.pop()
        if side * (weight[0] - weight[2]) <= 0:
            continue
        selected, secondary = select(weight)
        if selected in seen:
            continue
        counts = []
        for start, direction in rays:
            count = 0
            while add(start, tuple(count * x for x in direction)) in selected:
                count += 1
            counts.append(count)
        assert sum(counts) == n - 1
        last = [(i, add(a, tuple((k - 1) * x for x in b)))
                for i, ((a, b), k) in enumerate(zip(rays, counts)) if k]
        first = [(i, add(a, tuple(k * x for x in b))) for i, ((a, b), k) in enumerate(zip(rays, counts))]
        normals = {tuple(y - x for x, y in zip(a, b)) for i, a in last for j, b in first
                   if not (secondary == quota and i < 3 <= j)}
        normals.add((side, 0, -side))
        polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
        for normal in sorted(normals):
            polygon = clip(polygon, normal, 1)
        assert polygon and area2(polygon)
        strict = integer_weight(tuple(sum(p[i] for p in polygon) / len(polygon) for i in range(3)))
        assert min(strict) > 0 and all(dot(z, strict) > 0 for z in normals)
        assert select(strict)[0] == selected
        primary = sum(sum(x > 0 for x in a) <= 1 for a in selected)
        assert primary > n - primary
        if secondary == quota:
            prices = strict if side > 0 else tuple(reversed(strict))
            r = sum(counts[:3])
            assert r == n - quota - 1 and r >= quota
            axes = [tuple(j if i == axis else 0 for i in range(3))
                    for axis in range(3) for j in range(1, quota + 1)]
            ordinary = sorted(axes, key=lambda a: (dot(a, prices), a))[:quota]
            smaller = tuple(sum(a[i] > 0 for a in ordinary) for i in range(3))
            assert all(abs(r * d - quota * a) < 6 * r for a, d in zip(counts[:3], smaller))
            delta_x, delta_z = counts[3] - smaller[0], counts[5] - smaller[2]
            assert 0 <= delta_x <= 1 and 0 <= delta_z <= 1
            assert counts[4] == smaller[1] - delta_x - delta_z
        seen[selected] = {'weight': strict, 'ray_counts': counts, 'polygon': polygon,
                          'quota_saturated': secondary == quota}
        todo.extend(adjacent_weights(polygon, sorted(normals), comparisons))
    assert sum(area2(v['polygon']) for v in seen.values()) == Fraction(1, 2)
    return seen


def column_prefix_ranks(ordered, columns, characteristic):
    pivots, ranks = {}, []
    for a in ordered:
        vector = [x % characteristic for x in columns[a]] if characteristic else list(map(Fraction, columns[a]))
        for i, row in sorted(pivots.items()):
            coefficient = vector[i]
            vector = [x - coefficient * y for x, y in zip(vector, row)]
            if characteristic:
                vector = [x % characteristic for x in vector]
        index = next((i for i, x in enumerate(vector) if x), None)
        if index is not None:
            inverse = pow(vector[index], -1, characteristic) if characteristic else 1 / vector[index]
            pivots[index] = tuple(x * inverse % characteristic if characteristic else x * inverse for x in vector)
        ranks.append(len(pivots))
    return ranks


def two_curve_point_controls(n, component_sizes, pieces, characteristic):
    m, alpha = component_sizes[0], tuple(range(1, len(component_sizes) + 1))
    assert sum(component_sizes) == n
    blocks, offset = [], 0
    for size in component_sizes:
        blocks.append((offset, offset + size))
        offset += size
    ground = tuple(hyperbolic_cross(n, 3))
    assert all(a in ground for selected in pieces for a in selected)
    rng = Random(20261020 + n + characteristic)
    for attempt in range(1, 31):
        nodes = rng.sample(range(2, min(characteristic, 97) if characteristic else 97), n)

        def residue(a):
            if characteristic:
                return tuple(pow(t, a[0] - a[1], characteristic) * pow(c, a[1], characteristic) *
                             pow(t - 1, -a[2], characteristic) % characteristic
                             for c, (lo, hi) in zip(alpha, blocks) for t in nodes[lo:hi])
            return tuple(Fraction(t)**(a[0] - a[1]) * c**a[1] * Fraction(t - 1)**(-a[2])
                         for c, (lo, hi) in zip(alpha, blocks) for t in nodes[lo:hi])

        columns = {a: residue(a) for a in ground}
        if all(column_prefix_ranks(selected, columns, characteristic)[-1] == n for selected in pieces):
            break
    else:
        raise AssertionError('bounded search failed to find generic sample witnesses')
    # Independent rational-function columns after a common denominator.
    A, B, C = (max(a[i] for a in ground) for i in range(3))
    functional = [{} for _ in alpha]
    for a in ground:
        polynomial = [0] * (A + B + C + 1)
        for j in range(C - a[2] + 1):
            polynomial[a[0] - a[1] + B + j] = comb(C - a[2], j) * (-1)**(C - a[2] - j)
        for component, c in enumerate(alpha):
            functional[component][a] = tuple(c**a[1] * x for x in polynomial)
    union = {a: sum((model[a] for model in functional), ()) for a in ground}
    total_area, pairs, reductions, prefix_checks = Fraction(0), 0, 0, 0
    scalar = (lambda x: x % characteristic) if characteristic else Fraction
    for selected, data in pieces.items():
        weight = data['weight']
        assert greedy(ground, columns, weight, characteristic) == set(selected)
        ordered = sorted(ground, key=lambda a: (dot(a, weight), a))
        models = (functional[0], functional[1], union, columns) if len(alpha) == 2 else (union, columns)
        ranks = [column_prefix_ranks(ordered, model, characteristic) for model in models]
        for i in range(len(ordered)):
            if i + 1 == len(ordered) or dot(ordered[i], weight) != dot(ordered[i + 1], weight):
                if len(alpha) == 2:
                    d1, d2, d, actual = (r[i] for r in ranks)
                    assert actual == min(n, component_sizes[1] + d1, m + d2, d)
                else:
                    assert len(set(component_sizes)) == 1
                    d, actual = (r[i] for r in ranks)
                    assert actual == min(n, d)
                prefix_checks += 1
        inverse = field_inverse(selected, columns, characteristic)

        def coordinates(a):
            return [scalar(sum(x * y for x, y in zip(row, residue(a)))) for row in inverse]

        polynomials, normals = [], set()
        for a in minimal_excluded_border(selected):
            polynomial = {a: scalar(1)}
            for b, c in zip(selected, coordinates(a)):
                if c:
                    polynomial[b] = -c
                    normals.add(tuple(x - y for x, y in zip(a, b)))
            assert all(dot(a, weight) > dot(b, weight) for b in polynomial if b != a)
            polynomials.append(polynomial)
        polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
        for normal in sorted(normals):
            polygon = clip(polygon, normal, 1)
        assert polygon and area2(polygon)
        total_area += area2(polygon)
        pairs += pair_checks(polynomials, weight, characteristic)
        for a in ground:
            expected = {b: c for b, c in zip(selected, coordinates(a)) if c}
            assert normal_form({a: 1}, polynomials, weight, characteristic) == expected
            reductions += 1
    assert total_area == 1
    return {'length': n, 'component_sizes': component_sizes, 'characteristic': characteristic,
            'alpha': alpha, 'nodes_by_component': [nodes[lo:hi] for lo, hi in blocks],
            'generic_sample_attempts': attempt, 'complete_initial_ideals': len(pieces),
            'complete_simplex_double_area': str(total_area), 'ordinary_s_pairs': pairs,
            'ordinary_normal_forms': reductions, 'generic_rank_formula_prefix_checks': prefix_checks}


def two_curve_controls():
    symbolic = two_curve_symbolic_controls()
    fans, points = [], []
    for component_sizes in [(2, 2), (3, 3), (4, 4), (3, 1), (4, 2), (5, 1), (5, 3)]:
        n = sum(component_sizes)
        domains = [two_curve_cut_fan(n, side, min(component_sizes)) for side in (1, -1)]
        pieces = dict(domains[0])
        pieces.update(domains[1])
        fans.append({'length': n, 'component_sizes': component_sizes,
                     'positive_domain_piece_counts': [len(d) for d in domains],
                     'distinct_cut_sets': len(pieces),
                     'quota_saturated_sets': sum(d['quota_saturated'] for d in pieces.values()),
                     'cut_witnesses': [{'basis': selected, 'weight': data['weight']}
                                       for selected, data in sorted(pieces.items())]})
        points.extend(two_curve_point_controls(n, component_sizes, pieces, p) for p in (0, 101, 1009))
    special = tuple((0, j, 0) for j in range(4))
    special_columns = {a: tuple(Fraction(c, t)**a[1] for c, nn in [(1, (2, 3)), (2, (4, 5))] for t in nn)
                       for a in special}
    special_rank = column_prefix_ranks(special, special_columns, 0)[-1]
    assert special_rank == 3
    return {'scope': 'Generic samples of arbitrary positive component sizes on two scaled three-pole curves have Theta(n^2) initial ideals. Special samples are not covered; no unrestricted G_3 improvement or novelty claim.',
            'symbolic_curve_bases': symbolic, 'complete_six_ray_cut_fans': fans,
            'generic_point_quotients': points,
            'nongeneric_negative_control': {'alpha': [1, 2], 'nodes_by_component': [[2, 3], [4, 5]],
                                            'proposed_basis': special, 'rank': special_rank,
                                            'obstruction': 'The first node on each component has y=1/2.'}}


def multicurve_allowed(a, r, weight):
    if not r or a[0] * a[2]:
        return False
    ell = r * weight[0] // (weight[0] + weight[2])
    return a[1] < r or (a[0] < r - ell and a[2] < ell + 1)


def multicurve_symbolic_controls():
    cases = []
    for r in range(1, 8):
        zero = (0,) * (r + 3)
        relation = {(1, 0, 1) + (0,) * r: 1, (0, 0, 1) + (0,) * r: -1, zero: -1}
        f = {(r, r, 0) + (0,) * r: 1}
        for j in range(r):
            a = [j, j, 0] + [0] * r
            a[3 + j] = 1
            f[tuple(a)] = 1
        h = [f]
        for _ in range(r):
            h.append(normal_form({add(a, (0, 0, 1) + (0,) * r): c for a, c in h[-1].items()},
                                 [relation], (1, 1, 1) + (0,) * r, 0))
        for ell in range(r):
            weight = (2 * ell + 1, 2 * r + 1, 2 * r - 2 * ell - 1) + (0,) * r
            basis = [relation, h[ell], h[ell + 1]]
            heads = [leading(g, weight) for g in basis]
            assert [a[:3] for a in heads] == [(1, 0, 1), (r - ell, r, 0), (0, r, ell + 1)]

            def lower(a, b):
                dx, dy, dz = (x - y for x, y in zip(b[:3], a[:3]))
                values = (dy, ell * dx + (r - ell) * dz,
                          (ell + 1) * dx + (r - ell - 1) * dz)
                return min(values) >= 0 and max(values) > 0

            assert all(lower(a, head) for head, g in zip(heads, basis) for a in g if a != head)
            steps = 0
            for i, j in combinations(range(3), 2):
                common = tuple(max(a, b) for a, b in zip(heads[i], heads[j]))
                first = {add(a, tuple(x - y for x, y in zip(common, heads[i]))): c
                         for a, c in basis[i].items()}
                pending = subtract_shift(first, basis[j], tuple(x - y for x, y in zip(common, heads[j])), 1, 0)
                while pending:
                    a = leading(pending, weight)
                    assert lower(a, common)
                    for head, g in zip(heads, basis):
                        if all(x >= y for x, y in zip(a, head)):
                            pending = subtract_shift(pending, g, tuple(x - y for x, y in zip(a, head)), pending[a], 0)
                            steps += 1
                            break
                    else:
                        raise AssertionError('nonzero universal multicurve S-pair')
            cases.append({'components': r, 'sector': ell, 'heads': [a[:3] for a in heads],
                          'formal_s_pairs': 3, 'universally_lower_reduction_steps': steps})
    return cases


def floor_sum_integer(n, modulus, slope, offset):
    """Classical Euclidean floor sum; Python integers, including signed inputs."""
    qa, slope = divmod(slope, modulus)
    qb, offset = divmod(offset, modulus)
    answer = qa * n * (n - 1) // 2 + qb * n
    while True:
        qa, slope = divmod(slope, modulus)
        qb, offset = divmod(offset, modulus)
        answer += qa * n * (n - 1) // 2 + qb * n
        height = slope * n + offset
        if height < modulus:
            return answer
        n, offset = divmod(height, modulus)
        modulus, slope = slope, modulus


def surface_prefix_count(weight, threshold):
    if threshold < 0:
        return 0
    x, y, z = weight
    rows = threshold // y + 1
    return (rows + floor_sum_integer(rows, x, -y, threshold) +
            rows + floor_sum_integer(rows, z, -y, threshold) - rows)


def subcurve_prefix_counts(r, weight, threshold):
    x, y, z = weight
    axis_cost, next_x, next_z, total, ranks = 0, x, z, 0, []
    for s in range(1, r + 1):
        level = threshold - (s - 1) * y
        axis_count = 0 if level < 0 else 1 + level // x + level // z
        total += max(0, axis_count - s + 1)
        total += max(0, (threshold - axis_cost) // y - s + 1)
        ranks.append(total)
        if next_x < next_z:
            axis_cost, next_x = next_x, next_x + x
        else:
            axis_cost, next_z = next_z, next_z + z
    return ranks


def component_capacity_profile(sizes):
    assert sizes and min(sizes) >= 1
    n = sum(sizes)
    capacities, prefix = [], 0
    for size in nlargest(min(len(sizes), isqrt(n - 1)), sizes):
        prefix += size
        capacities.append(prefix)
    return n, capacities


def generic_multicurve_prefix_rank(profile, weight, threshold):
    n, capacities = profile
    ranks = subcurve_prefix_counts(len(capacities), weight, threshold)
    return min([n, surface_prefix_count(weight, threshold)] +
               [n - capacity + rank for capacity, rank in zip(capacities, ranks)])


def multicurve_rank_controls():
    def partitions(n, maximum):
        if n == 0:
            yield ()
        for first in range(min(n, maximum), 0, -1):
            for rest in partitions(n - first, first):
                yield (first,) + rest

    floor_stress = 0
    for n, modulus, slope, offset in product(range(21), (1, 2, 7, 101, 10**20 + 3),
                                            (-10**20, -3, 0, 5, 10**20 + 7),
                                            (-10**20 - 1, -1, 0, 11, 10**20 + 9)):
        assert floor_sum_integer(n, modulus, slope, offset) == sum((slope * i + offset) // modulus for i in range(n))
        floor_stress += 1
    profiles = [p for n in range(1, 9) for p in partitions(n, n)]
    rank_checks, hook_checks, floor_checks = 0, 0, 0
    for weight in product(range(1, 5), repeat=3):
        for threshold in range(13):
            ground = tuple((a, b, c) for a, b, c in product(range(threshold + 1), repeat=3)
                           if not a * c and dot((a, b, c), weight) <= threshold)
            assert surface_prefix_count(weight, threshold) == len(ground)
            floor_checks += 1
            direct = [sum(multicurve_allowed(a, r, weight) for a in ground) for r in range(1, 9)]
            assert direct == subcurve_prefix_counts(8, weight, threshold)
            increments = [direct[0]] + [b - a for a, b in zip(direct, direct[1:])]
            assert all(a >= b for a, b in zip(increments, increments[1:]))
            hook_checks += 8
            for sizes in profiles:
                n, prefix, expected = sum(sizes), 0, min(sum(sizes), direct[len(sizes) - 1])
                for size, d in zip(sizes, direct):
                    prefix += size
                    expected = min(expected, n - prefix + d)
                assert generic_multicurve_prefix_rank(component_capacity_profile(sizes), weight, threshold) == expected
                rank_checks += 1
    # Independent rational-function numerators, all component subsets.
    matrix_checks = 0
    for r, prime in product((3, 5), (7, 101)):
        for weight in ((1, 1, 1), (1, 2, 3), (3, 2, 1)):
            for threshold in range(7):
                ground = tuple((a, b, c) for a, b, c in product(range(threshold + 1), repeat=3)
                               if not a * c and dot((a, b, c), weight) <= threshold)
                B, C = (max(a[i] for a in ground) for i in (1, 2))
                degree = max(a[0] for a in ground) + B + C
                functional = [{} for _ in range(r)]
                for a in ground:
                    polynomial = [0] * (degree + 1)
                    for j in range(C - a[2] + 1):
                        polynomial[a[0] - a[1] + B + j] = comb(C - a[2], j) * (-1)**(C - a[2] - j)
                    for i in range(r):
                        functional[i][a] = tuple((i + 1)**a[1] * c % prime for c in polynomial)
                for size in range(1, r + 1):
                    expected = sum(multicurve_allowed(a, size, weight) for a in ground)
                    for subset in combinations(range(r), size):
                        columns = {a: sum((functional[i][a] for i in subset), ()) for a in ground}
                        assert column_prefix_ranks(ground, columns, prime)[-1] == expected
                        matrix_checks += 1
    return {'capacity_partitions': len(profiles), 'compressed_rank_queries': rank_checks,
            'signed_large_integer_floor_sums': floor_stress,
            'direct_hook_counts': hook_checks, 'surface_floor_sum_counts': floor_checks,
            'independent_functional_subset_ranks': matrix_checks}


def multicurve_point_prefix_controls():
    prime, fixtures = 1009, []
    queries = [(weight, threshold) for weight in ((1, 1, 1), (1, 2, 3), (3, 2, 1)) for threshold in range(7)]
    for sizes in ((15, 15, 15), (30, 10, 5), (5, 30, 10), (20, 1, 1, 1, 1, 1, 1, 1)):
        n, profile = sum(sizes), component_capacity_profile(sizes)
        rng = Random(20261020 + n + len(sizes))
        query_grounds = {q: tuple((a, b, c) for a, b, c in product(range(q[1] + 1), repeat=3)
                                  if not a * c and dot((a, b, c), q[0]) <= q[1]) for q in queries}
        ground = set().union(*map(set, query_grounds.values()))
        for attempt in range(1, 31):
            nodes = rng.sample(range(2, prime), n)
            blocks, offset = [], 0
            for size in sizes:
                blocks.append(nodes[offset:offset + size])
                offset += size
            columns = {a: tuple(pow(t, a[0] - a[1], prime) * pow(j, a[1], prime) *
                                pow(t - 1, -a[2], prime) % prime
                                for j, nn in enumerate(blocks, 1) for t in nn) for a in ground}
            results = []
            for weight, threshold in queries:
                actual = column_prefix_ranks(query_grounds[(weight, threshold)], columns, prime)[-1]
                expected = generic_multicurve_prefix_rank(profile, weight, threshold)
                if actual != expected:
                    break
                results.append({'weight': weight, 'threshold': threshold, 'evaluation_rank': actual,
                                'surface_prefix_count': surface_prefix_count(weight, threshold)})
            else:
                break
        else:
            raise AssertionError('no generic finite-prefix evaluation witness found')
        assert len(results) == len(queries)
        defects = [q for q in results if q['evaluation_rank'] < min(n, q['surface_prefix_count'])]
        assert defects, sizes
        fixtures.append({'component_sizes': sizes, 'length': n, 'characteristic': prime,
                         'nodes_by_component': blocks, 'generic_sample_attempts': attempt,
                         'retained_capacity_constraints': len(profile[1]),
                         'direct_prefix_rank_queries': len(results), 'queries': results,
                         'rank_defect_witness': defects[0]})
    return fixtures


def multicurve_balanced_fan(r, m):
    n = r * m
    ground = tuple(a for a in hyperbolic_cross(n, 3) if not a[0] * a[2])
    comparisons = {primitive(tuple(x - y for x, y in zip(a, b))) for a, b in combinations(ground, 2)}
    comparisons.update(primitive((r - h, 0, -h)) for h in range(1, r))
    pieces, domain_counts = {}, []
    for ell in range(r):
        todo, seen = [(2 * ell + 1, 3 * r + 1, 2 * r - 2 * ell - 1)], {}
        allowed = tuple(a for a in ground if a[1] < r or (a[0] < r - ell and a[2] < ell + 1))
        lower, upper = (r - ell, 0, -ell), (-(r - ell - 1), 0, ell + 1)
        while todo:
            weight = todo.pop()
            if dot(lower, weight) <= 0 or dot(upper, weight) <= 0:
                continue
            selected = tuple(sorted(sorted(allowed, key=lambda a: (dot(a, weight), a))[:n]))
            if selected in seen:
                continue
            normals = {tuple(y - x for x, y in zip(a, b)) for a in selected for b in allowed if b not in selected}
            normals.update((lower, upper))
            polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
            for normal in sorted(normals):
                polygon = clip(polygon, normal, 1)
            assert polygon and area2(polygon)
            strict = integer_weight(tuple(sum(p[i] for p in polygon) / len(polygon) for i in range(3)))
            assert all(dot(z, strict) > 0 for z in normals)
            if n <= r * r:
                assert selected == tuple(sorted(sorted(ground, key=lambda a: (dot(a, strict), a))[:n]))
            seen[selected] = {'weight': strict, 'polygon': polygon}
            todo.extend(adjacent_weights(polygon, sorted(normals), sorted(comparisons)))
        assert sum(area2(v['polygon']) for v in seen.values()) == Fraction(1, r)
        domain_counts.append(len(seen))
        pieces.update(seen)
    return pieces, {'components': r, 'points_per_component': m, 'length': n,
                    'carrier_domain_piece_counts': domain_counts, 'distinct_initial_ideals': len(pieces),
                    'surface_fan_equality': n <= r * r,
                    'cut_witnesses': [{'basis': selected, 'weight': data['weight']}
                                      for selected, data in sorted(pieces.items())]}


def multicurve_controls():
    formal, rank_controls = multicurve_symbolic_controls(), multicurve_rank_controls()
    fans, points = [], []
    for r, m in ((3, 2), (3, 3), (4, 2), (5, 1)):
        pieces, fan = multicurve_balanced_fan(r, m)
        fans.append(fan)
        points.extend(two_curve_point_controls(r * m, (m,) * r, pieces, p) for p in (0, 1009))
    prefixes = multicurve_point_prefix_controls()
    return {'scope': 'Uniform carrier and generic sampling-rank formulas for any number of scaled components. Balanced samples with r^2>=n have exactly the generic surface fan. Intermediate balanced counts and unrestricted G_3 remain unresolved; no novelty certificate.',
            'formal_parameter_sectors': formal, 'rank_controls': rank_controls,
            'complete_balanced_model_fans': fans, 'generic_point_quotients': points,
            'finite_prefix_rank_fixtures': prefixes}


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
    unions = two_curve_controls()
    multicurves = multicurve_controls()
    paths = ['scripts/run_algebraic_incidence_controls.py', 'docs/weighted_facet_incidence.md',
             'docs/initial_cone_height.md', 'scripts/run_extremal_initial_ideals_checks.py',
             'scripts/run_structural_complexity_checks.py', 'scripts/run_initial_cone_height_checks.py',
             'scripts/run_interleaving_research.py', 'docs/extremal_initial_ideals.md',
             'docs/generic_multiaffine_bound.md']
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_algebraic_incidence_controls.py',
              'scope': 'Three-pole quotients give quartic incidence and generic two-curve unions have quadratic fans. Arbitrarily many scaled components have an exact generic prefix-rank formula with only square-root many capacity constraints. The full extremal count and priority remain open.',
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
                          'plane_comparison_intervals': sum(a['complete_comparison_intervals'] for a in planes),
                          'two_curve_formal_s_pairs': 6,
                          'two_curve_complete_cut_fans': len(unions['complete_six_ray_cut_fans']),
                          'two_curve_complete_point_fans': len(unions['generic_point_quotients']),
                          'two_curve_quota_allocation_checks': sum(a['quota_saturated_sets'] for a in unions['complete_six_ray_cut_fans']),
                          'two_curve_rank_formula_prefix_checks': sum(a['generic_rank_formula_prefix_checks'] for a in unions['generic_point_quotients']),
                          'two_curve_s_pairs': sum(a['ordinary_s_pairs'] for a in unions['generic_point_quotients']),
                          'two_curve_normal_forms': sum(a['ordinary_normal_forms'] for a in unions['generic_point_quotients']),
                          'multicurve_formal_sectors': len(multicurves['formal_parameter_sectors']),
                          'multicurve_formal_s_pairs': sum(a['formal_s_pairs'] for a in multicurves['formal_parameter_sectors']),
                          'multicurve_complete_point_fans': len(multicurves['generic_point_quotients']),
                          'multicurve_prefix_rank_checks': sum(a['generic_rank_formula_prefix_checks'] for a in multicurves['generic_point_quotients']),
                          'multicurve_normal_forms': sum(a['ordinary_normal_forms'] for a in multicurves['generic_point_quotients']),
                          'multicurve_direct_point_rank_queries': sum(a['direct_prefix_rank_queries'] for a in multicurves['finite_prefix_rank_fixtures']),
                          **multicurves['rank_controls']},
              'uniform_axis_matroids': axes, 'arithmetic_lower_bound_controls': arithmetic,
              'interior_quotients': interior, 'three_pole_quotients': curves,
              'positive_plane_traces': planes,
              'generic_two_curve_unions': unions,
              'generic_multicurve_unions': multicurves,
              'source_sha256': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}}
    (ROOT / 'results/research/algebraic_incidence_controls.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report['summary'], indent=2))


if __name__ == '__main__':
    main()
