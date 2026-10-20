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
from run_planar_dyadic_bound_checks import convex_hull

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


def hook_threshold_group(r, n, A, B, level):
    planar = [(a, b, 0) for a in range(level // A + 1) for b in range(level // B + 1)
              if a < r or b < r]
    cheap = [a for a in planar if A * a[0] + B * a[1] < level]
    tied = [a for a in planar if A * a[0] + B * a[1] == level]
    k, gain = len(cheap), len(tied)
    if k >= n or gain < 2:
        return None
    capped_gain = min(gain, n - k)
    changing_ranks = []
    for take in range(1, min(gain - 1, n - k) + 1):
        j = k + take
        before = set(cheap + sorted(tied, key=lambda a: a[0])[:take])
        after = set(cheap + sorted(tied, key=lambda a: a[0], reverse=True)[:take])
        assert before != after and len(before) == len(after) == j
        ground = [(a, b, 0) for a, b in hyperbolic_cross(j, 2) if a < r or b < r]
        scale = 2 * j + 1
        for basis, sign in ((before, 1), (after, -1)):
            weight = (scale * A + sign, scale * B, 0)
            assert set(sorted(ground, key=lambda a: (dot(a, weight), a))[:j]) == basis
        changing_ranks.append(j)
    assert capped_gain <= 2 * len(changing_ranks)

    low_z, high_z = Fraction(level, n + 2), Fraction(level + 1)
    carriers = {Fraction(A * (r - h), h) for h in range(1, r)}
    events = {low_z, high_z} | carriers
    events.update(Fraction(level - B * b, c) for b in range(level // B + 1)
                  for c in range(1, n + 3) if level - B * b > 0)
    events = sorted(z for z in events if low_z <= z <= high_z)
    samples = sorted(set(events) | {(a + b) / 2 for a, b in zip(events, events[1:])})
    strict_ranks, closed_ranks, critical, vertices = [], [], [], []
    for z in samples:
        denominator = z.denominator
        weight = (A * denominator, B * denominator, z.numerator)
        strict = subcurve_prefix_counts(r, weight, level * denominator - 1)[-1]
        closed = subcurve_prefix_counts(r, weight, level * denominator)[-1]
        strict_ranks.append(strict)
        closed_ranks.append(closed)
        critical.append(strict < n <= closed)
        if z not in events or z in carriers:
            continue
        xy = [a for a in tied if multicurve_allowed(a, r, (A, B, z))]
        yz = []
        for b in range(level // B + 1):
            c = Fraction(level - B * b) / z
            if c > 0 and c.denominator == 1 and multicurve_allowed((0, b, int(c)), r, (A, B, z)):
                yz.append((0, b, int(c)))
        assert closed - strict == len(xy) + len(yz)
        if strict < n < closed and len(xy) >= 2 and yz:
            vertices.append({'z_price': str(z), 'strict_rank': strict, 'closed_rank': closed,
                             'xy_threshold_points': xy, 'positive_z_threshold_points': yz})
    assert all(a >= b for a, b in zip(strict_ranks, strict_ranks[1:]))
    assert all(a >= b for a, b in zip(closed_ranks, closed_ranks[1:]))
    indices = [i for i, value in enumerate(critical) if value]
    assert not indices or indices == list(range(indices[0], indices[-1] + 1))
    assert len(vertices) <= capped_gain + 2 <= 3 * capped_gain
    return {'components': r, 'length': n, 'xy_price': [A, B], 'threshold': level,
            'planar_strict_rank': k, 'planar_tie_rank': gain, 'capped_tie_rank': capped_gain,
            'changing_planar_truncation_ranks': changing_ranks,
            'complete_transverse_samples': len(samples), 'off_carrier_vertices': vertices}


def hook_carrier_trace(r, n, A, B):
    assert gcd(A, B) == 1 and (r * A) % (A + B) == 0
    ground = tuple(a for a in hyperbolic_cross(n, 3) if multicurve_allowed(a, r, (A, 1, B)))
    projected = {a: (A * a[0] + B * a[2], a[1]) for a in ground}
    H = max(A, B) * n
    assert all((p + 1) * (b + 1) <= H for p, b in projected.values())
    events = {Fraction(0)}
    for a, b in combinations(ground, 2):
        dp, dy = (x - y for x, y in zip(projected[a], projected[b]))
        if dp * dy < 0:
            events.add(-Fraction(dp, dy))
    events = sorted(events)
    samples = [(a + b) / 2 for a, b in zip(events, events[1:])] + [events[-1] + 1]
    path, bases = [], {}
    for t in samples:
        weight = (A * t.denominator, t.numerator, B * t.denominator)
        selected = tuple(sorted(sorted(ground, key=lambda a: (dot(a, weight), a))[:n]))
        profile = tuple(sum(projected[a][i] for a in selected) for i in range(2))
        if not path or path[-1]['profile'] != profile:
            if path:
                assert path[-1]['profile'][0] < profile[0] and path[-1]['profile'][1] > profile[1]
            path.append({'y_price': str(t), 'profile': profile})
        bases[selected] = weight
    cube = 0
    while cube**3 < n * n * H:
        cube += 1
    bound = 1 + (H - 1).bit_length() * (n + 3 * cube)
    assert len(path) <= bound
    return {'components': r, 'length': n, 'primitive_axis_prices': [A, B],
            'projected_hyperbolic_bound': H, 'complete_comparison_intervals': len(samples),
            'distinct_projected_profiles': len(path), 'explicit_planar_upper_bound': bound,
            'profile_path': path}, bases


def hook_active_local_control(trace_bases):
    r, n, level, prime = 3, 36, 6, 1000003
    pool = [a for a in product(range(level + 1), repeat=3)
            if not a[0] * a[2] and sum(a) <= level and multicurve_allowed(a, r, (1, 1, 1))]
    cheap, tied = [a for a in pool if sum(a) < level], [a for a in pool if sum(a) == level]
    assert len(cheap) == 34 and len(tied) == 9
    take = n - len(cheap)
    profiles = defaultdict(list)
    for choice in combinations(tied, take):
        profiles[tuple(sum(a[i] for a in choice) for i in range(2))].append(choice)
    hull = convex_hull(profiles)
    ground = tuple(hyperbolic_cross(n, 3))
    local = {}
    for index, p in enumerate(hull):
        previous, following = hull[index - 1], hull[(index + 1) % len(hull)]
        dx, dy = -(p[1] - previous[1]) - (following[1] - p[1]), (p[0] - previous[0]) + (following[0] - p[0])
        assert all(dx * (q[0] - p[0]) + dy * (q[1] - p[1]) > 0 for q in hull if q != p)
        assert len(profiles[p]) == 1
        selected = tuple(sorted(cheap + list(profiles[p][0])))
        scale = 100 * n * max(1, abs(dx), abs(dy))
        weight = (scale + dx, scale + dy, scale)
        allowed = [a for a in ground if multicurve_allowed(a, r, weight)]
        assert tuple(sorted(sorted(allowed, key=lambda a: (dot(a, weight), a))[:n])) == selected
        local[selected] = weight
    all_bases = dict(trace_bases)
    all_bases.update(local)
    rng = Random(20261020 + n)
    for attempt in range(1, 31):
        nodes = rng.sample(range(2, prime), n)

        def residue(a):
            return tuple(pow(t, a[0] - a[1], prime) * pow(j, a[1], prime) * pow(t - 1, -a[2], prime) % prime
                         for j in range(1, r + 1) for t in nodes[12 * (j - 1):12 * j])

        columns = {a: residue(a) for a in ground}
        if all(column_prefix_ranks(selected, columns, prime)[-1] == n for selected in all_bases):
            break
    else:
        raise AssertionError('no active local/trace point witness found')
    for selected, weight in local.items():
        assert greedy(ground, columns, weight, prime) == set(selected)
    for selected, weight in list(trace_bases.items())[::max(1, len(trace_bases) // 4)]:
        assert greedy(ground, columns, weight, prime) == set(selected)
    selected, weight = next(iter(local.items()))
    inverse = field_inverse(selected, columns, prime)

    def coordinates(a):
        return [sum(x * y for x, y in zip(row, residue(a))) % prime for row in inverse]

    polynomials = []
    for a in minimal_excluded_border(selected):
        polynomial = {a: 1}
        for b, c in zip(selected, coordinates(a)):
            if c:
                polynomial[b] = -c
        assert all(dot(a, weight) > dot(b, weight) for b in polynomial if b != a)
        polynomials.append(polynomial)
    pairs = pair_checks(polynomials, weight, prime)
    for a in ground:
        expected = {b: c for b, c in zip(selected, coordinates(a)) if c}
        assert normal_form({a: 1}, polynomials, weight, prime) == expected
    return {'components': r, 'length': n, 'characteristic': prime, 'points_per_component': 12,
            'nodes_by_component': [nodes[i:i + 12] for i in range(0, n, 12)],
            'generic_sample_attempts': attempt, 'threshold_weight': [1, 1, 1], 'threshold': level,
            'strict_carrier_rank': len(cheap), 'closed_carrier_rank': len(cheap) + len(tied),
            'local_exposed_bases': len(local), 'all_local_and_trace_basis_rank_checks': len(all_bases),
            'ordinary_s_pairs': pairs, 'ordinary_normal_forms': len(ground),
            'local_basis_witnesses': [{'basis': basis, 'strict_weight': weight} for basis, weight in local.items()]}


def hook_charging_controls():
    groups = []
    for r, n in ((1, 8), (2, 16), (3, 36), (4, 64)):
        for A, B, level in product(range(1, 4), range(1, 4), range(1, 13)):
            if gcd(A, B) == 1:
                group = hook_threshold_group(r, n, A, B, level)
                if group is not None:
                    groups.append(group)
    assigned = defaultdict(set)
    for group in groups:
        key = (group['components'], group['length'], *group['xy_price'])
        ranks = set(group['changing_planar_truncation_ranks'])
        assert not (assigned[key] & ranks)
        assigned[key].update(ranks)
    traces, bases = [], {}
    for A, B in ((1, 2), (2, 1)):
        trace, selected = hook_carrier_trace(3, 36, A, B)
        traces.append(trace)
        bases.update(selected)
    active = hook_active_local_control(bases)
    return {'scope': 'Uniform O(n^2*log(n+1)) balanced multicomponent fan bound. Threshold charges share a monotone rank interval across all carrier sectors. The general unequal-size and unrestricted G_3 counts remain open; priority is not certified.',
            'planar_tie_groups': groups, 'complete_carrier_traces': traces,
            'active_local_point_quotient': active}


def unequal_phase_controls():
    """Compare nested quotas, active-prefix indices and the floor-sum oracle."""
    def partitions(n, maximum):
        if not n:
            yield ()
        for first in range(min(n, maximum), 0, -1):
            for rest in partitions(n - first, first):
                yield (first,) + rest

    counts = defaultdict(int)
    examples = []
    for n in range(1, 9):
        ground = tuple(hyperbolic_cross(n, 3))
        for sizes in partitions(n, n):
            r, capacities = len(sizes), [0]
            for size in sizes:
                capacities.append(capacities[-1] + size)
            for base_weight in product(range(1, 4), repeat=3):
                # Preserve integer cost gaps and resolve all ties on V_n^3.
                weight = tuple((n**3 + 1) * b + d for b, d in
                               zip(base_weight, (n*n, n, 1)))
                allowed = {a: tuple(multicurve_allowed(a, s, weight)
                                    for s in range(r + 1)) for a in ground}
                ordered = sorted((a for a in ground if allowed[a][r]),
                                 key=lambda a: dot(a, weight))
                assert len({dot(a, weight) for a in ordered}) == len(ordered)
                ranks, outside, selected = [0] * (r + 1), [0] * (r + 1), []
                smallest, largest, phases = r, r, [r]
                for a in ordered:
                    feasible = all(outside[s] + (not allowed[a][s]) <= n - capacities[s]
                                   for s in range(r + 1))
                    predicted = allowed[a][smallest]
                    assert feasible == predicted
                    old_rank = n - max(c - d for c, d in zip(capacities, ranks))
                    ranks = [d + member for d, member in zip(ranks, allowed[a])]
                    values = [c - d for c, d in zip(capacities, ranks)]
                    active = [s for s in range(r + 1) if values[s] == max(values)]
                    assert min(active) <= smallest and max(active) <= largest
                    assert n - max(values) - old_rank == int(predicted)
                    assert generic_multicurve_prefix_rank(component_capacity_profile(sizes), weight,
                                                          dot(a, weight)) == n - max(values)
                    if feasible:
                        selected.append(a)
                        outside = [count + (not member) for count, member in zip(outside, allowed[a])]
                    else:
                        assert min(active) == smallest
                    if min(active) != smallest:
                        assert outside[min(active)] == n - capacities[min(active)]
                        phases.append(min(active))
                    smallest, largest = min(active), max(active)
                    counts['prefix_decisions'] += 1
                    if len(selected) == n:
                        break
                assert len(selected) == n and phases[-1] == 0
                selected_set = set(selected)
                assert all(tuple(v - int(i == j) for i, v in enumerate(a)) in selected_set
                           for a in selected for j in range(3) if a[j])
                for s, t in zip(phases, phases[1:]):
                    block = [a for a in ordered if allowed[a][s] and not allowed[a][t]]
                    take = capacities[s] - capacities[t]
                    assert selected_set.intersection(block) == set(block[:take])
                    cutoff = dot(block[take - 1], weight)
                    # Search the entire two-plane cost prefix, including points
                    # outside the universal hyperbolic ground used above.
                    full = [(a, b, c) for b in range(cutoff//weight[1] + 1)
                            for a, c in ([(a, 0) for a in range((cutoff - b*weight[1])//weight[0] + 1)] +
                                         [(0, c) for c in range(1, (cutoff - b*weight[1])//weight[2] + 1)])]
                    assert {a for a in full if multicurve_allowed(a, s, weight) and
                            not multicurve_allowed(a, t, weight)} == set(block[:take])
                    counts['frozen_blocks'] += 1
                # The same extremal monotonicity must hold without resolving ties.
                previous = (r, r)
                for threshold in range(13):
                    direct = [0] + subcurve_prefix_counts(r, base_weight, threshold)
                    values = [c - d for c, d in zip(capacities, direct)]
                    active = [s for s in range(r + 1) if values[s] == max(values)]
                    assert min(active) <= previous[0] and max(active) <= previous[1]
                    previous = min(active), max(active)
                    counts['closed_tie_thresholds'] += 1
                counts['price_profiles'] += 1
                if sizes == (4, 2, 1) and base_weight in ((1, 2, 3), (3, 2, 1)):
                    examples.append({'sizes': sizes, 'weight': weight,
                                     'decreasing_active_indices': phases, 'selected': sorted(selected)})
    return {'scope': 'A fixed carrier-flag domain has the same positive greedy outputs as a classical nested matroid. This is not an identification of the full residue matroid on arbitrary subsets or a fan-count bound.',
            'counts': dict(counts), 'phase_witnesses': examples}


def hook_capacity_jump_controls():
    """Direct Ferrers hooks versus every capacity profile through length 12."""
    def partitions(n, maximum):
        if not n:
            yield ()
        for first in range(min(n, maximum), 0, -1):
            for rest in partitions(n - first, first):
                yield (first,) + rest

    profiles = [sizes for n in range(1, 13) for sizes in partitions(n, n)]
    queries, slow, raw = 0, 0, 0
    for weight in product(range(1, 5), repeat=3):
        x, y, z = weight
        for threshold in range(16):
            axes = sorted([0] + [a*x for a in range(1, threshold//x + 1)] +
                          [c*z for c in range(1, threshold//z + 1)])
            hooks = [0] * 12
            for i, cost in enumerate(axes):
                for b in range((threshold - cost)//y + 1):
                    if min(i, b) < 12:
                        hooks[min(i, b)] += 1
            assert all(a >= b + 2 for a, b in zip(hooks, hooks[1:]) if b)
            counts = [0]
            for h in hooks:
                counts.append(counts[-1] + h)
            assert counts[1:] == subcurve_prefix_counts(12, weight, threshold)
            raw += 1
            for sizes in profiles:
                r, capacities = len(sizes), [0]
                for size in sizes:
                    capacities.append(capacities[-1] + size)
                jumps = [s for s in range(1, r) if sizes[s - 1] - sizes[s] > 2]
                values = [capacity - count for capacity, count in zip(capacities, counts)]
                assert max(values) == max(values[s] for s in [0, *jumps, r])
                n = sum(sizes)
                retained = [s for s in [*jumps, r] if s*s < n]
                reduced = min([n, surface_prefix_count(weight, threshold)] +
                              [n - capacities[s] + counts[s] for s in retained])
                assert reduced == n - max(values)
                if not jumps:
                    assert max(values) == max(0, values[-1])
                    slow += 1
                queries += 1
    return {'scope': 'Only size drops greater than two can supply additional capacity constraints; combine this exact reduction with the earlier square-root cutoff. With no such drops, uniform carrier truncation and the existing quadratic-power fan bound apply. Diagonal-hook facts are classical; no unrestricted bound or exact priority claim.',
            'counts': {'hook_grid_prefixes': raw, 'capacity_jump_rank_queries': queries,
                       'slow_drop_rank_queries': slow},
            'lengths': [1, 12], 'integer_prices': [1, 4], 'integer_costs': [0, 15]}


def unequal_threshold_family(t):
    """Explicit two-cost parallelogram, uniformly for integer t >= 1."""
    assert t >= 1
    sizes = (31*t + 14, 31*t + 3, 31*t - 9)
    n = sum(sizes)
    first, middle, outer = 30*t + 14, 30*t + 8, 30*t + 1
    primary = ({(0, 0, 0)} | {(a, 0, 0) for a in range(1, first//2 + 1)} |
               {(0, b, 0) for b in range(1, first//3 + 1)} |
               {(0, 0, c) for c in range(1, first//5 + 1)})
    second = ({(a, 1, 0) for a in range(1, 15*t + 3)} |
              {(1, b, 0) for b in range(2, 10*t + 2)} |
              {(0, 1, c) for c in range(1, 6*t + 1)})
    third = ({(a, 2, 0) for a in range(2, 15*t - 2)} |
             {(2, b, 0) for b in range(3, 10*t - 1)} |
             {(0, 2, c) for c in range(1, 6*t - 1)})
    flip_pairs = [((1, 10*t + 2, 0), (0, 1, 6*t + 1)),
                  ((2, 10*t - 1, 0), (0, 2, 6*t - 1))]
    common = primary | second | third
    assert [len(primary), len(second) + 1, len(third) + 1] == list(sizes)
    assert len(common) == n - 2
    directions = [tuple(a - b for a, b in zip(*pair)) for pair in flip_pairs]
    assert all(gcd(*q) == 1 and all(q) for q in directions)
    cross = tuple(directions[0][(i+1) % 3] * directions[1][(i+2) % 3] -
                  directions[0][(i+2) % 3] * directions[1][(i+1) % 3] for i in range(3))
    assert cross == tuple(-(2*t + 1)*v for v in (2, 3, 5))
    cases = []
    for s2, s3 in product((-1, 1), repeat=2):
        dx, dy = (10*t + 1)*s3 - (10*t - 3)*s2, 2*s2 - s3
        scale = 10*n*max(1, abs(dx), abs(dy)) + 1
        weight = (2*scale + dx, 3*scale + dy, 5*scale)
        selected = common | {flip_pairs[0][int(s2 > 0)], flip_pairs[1][int(s3 > 0)]}
        assert [dot(q, weight) for q in directions] == [(10*t + 5)*s2, (10*t + 5)*s3]
        assert all(tuple(v - int(i == j) for i, v in enumerate(a)) in selected
                   for a in selected for j in range(3) if a[j])
        cases.append((weight, tuple(sorted(selected))))
    checks = [(outer - 1, n - 22), (outer, n - 20), (middle - 1, n - 8),
              (middle, n - 6), (first - 1, n - 1), (first, n)]
    assert all(generic_multicurve_prefix_rank(component_capacity_profile(sizes), (2, 3, 5), C) == R
               for C, R in checks)
    return sizes, common, flip_pairs, directions, cases, checks


def unequal_point_controls(t):
    sizes, common, flip_pairs, directions, cases, checks = unequal_threshold_family(t)
    n, prime, rng = sum(sizes), 1000003, Random(20261020 + t)
    basis_ground = common | {a for pair in flip_pairs for a in pair}
    for attempt in range(1, 31):
        nodes = rng.sample(range(2, prime), n)
        components, start = [], 0
        for j, size in enumerate(sizes, 1):
            components.extend((j, node) for node in nodes[start:start + size])
            start += size

        def residue(a):
            return tuple(pow(node, a[0] - a[1], prime) * pow(j, a[1], prime) *
                         pow(node - 1, -a[2], prime) % prime for j, node in components)

        columns = {a: residue(a) for a in basis_ground}
        try:
            inverses = [field_inverse(selected, columns, prime) for _, selected in cases]
        except StopIteration:
            continue
        break
    else:
        raise AssertionError('no simultaneous unequal-component point witness')
    certificates, total_border, strict_coefficients = [], 0, 0
    for (weight, selected), inverse in zip(cases, inverses):
        border = sorted(minimal_excluded_border(selected))
        for a in border:
            column = columns.setdefault(a, residue(a))
            coordinates = [sum(x*y for x, y in zip(row, column)) % prime for row in inverse]
            assert all(sum(coordinates[j] * columns[b][i] for j, b in enumerate(selected)) % prime == column[i]
                       for i in range(n))
            assert all(dot(a, weight) > dot(b, weight) for b, c in zip(selected, coordinates) if c)
            strict_coefficients += sum(bool(c) for c in coordinates)
        total_border += len(border)
        certificates.append({'weight': weight, 'flip_choices': [next(i for i, a in enumerate(pair) if a in selected)
                                                               for pair in flip_pairs],
                             'minimal_initial_generators': border})
    # Exact raw evaluation ranks, independently of the nested/floor-sum oracle.
    ranks = []
    for threshold, expected in checks:
        pool = tuple((a, b, c) for b in range(threshold//3 + 1)
                     for a, c in ([(a, 0) for a in range((threshold - 3*b)//2 + 1)] +
                                  [(0, c) for c in range(1, (threshold - 3*b)//5 + 1)]))
        direct_columns = {a: columns.get(a) or residue(a) for a in pool}
        actual = column_prefix_ranks(pool, direct_columns, prime)[-1]
        assert actual == expected
        ranks.append({'cost': threshold, 'rank': actual})
    start, nodes_by_component = 0, []
    for size in sizes:
        nodes_by_component.append(nodes[start:start + size])
        start += size
    return {'parameter': t, 'length': n, 'component_sizes': sizes, 'characteristic': prime,
            'component_scalars': [1, 2, 3], 'nodes_by_component': nodes_by_component,
            'generic_sample_attempts': attempt, 'common_selected_monomials': sorted(common),
            'two_boundary_pairs': flip_pairs, 'primitive_exchange_directions': directions,
            'four_local_certificates': certificates, 'direct_prefix_ranks': ranks,
            'matrix_bases': len(cases), 'minimal_border_relations': total_border,
            'strictly_lower_nonzero_coefficients': strict_coefficients,
            'certificate': 'Four invertible evaluation bases and vanishing minimal-border relations with strictly lower tails prove the four actual initial ideals by inclusion and equal colength. No whole-fan enumeration or full residue-matroid equivalence is claimed.'}


def hook_area(weight, s, cost):
    """Continuous triangle intersected with the width-s diagonal hook strip."""
    x, y, z = weight
    alpha = Fraction(1, x) + Fraction(1, z)
    kappa = y + 1/alpha
    return alpha/(2*y) * (max(0, cost)**2 - max(0, cost - kappa*s)**2)


def linear_profile_limit(weight, slope, tau):
    x, y, z = weight
    alpha = Fraction(1, x) + Fraction(1, z)
    kappa, gamma = y + 1/alpha, alpha + Fraction(1, y)
    curvature = gamma*kappa
    if slope <= curvature:
        return min(Fraction(slope, 2), alpha/(2*y)*tau*tau)
    if tau <= kappa:
        return alpha/(2*y)*tau*tau
    if tau >= slope/gamma:
        return Fraction(slope, 2)
    return Fraction(slope, 2) - (slope - gamma*tau)**2/(2*(slope - curvature))


def hook_step_bound(weight):
    x, y, z = weight
    ceil_ratio = lambda a, b: (a + b - 1)//b
    return 2 + ceil_ratio(y, x) + ceil_ratio(y, z) + ceil_ratio(min(x, z), y)


def steep_rank_controls():
    prices = [(1, 1, 1), (2, 1, 2), (3, 2, 6), (5, 3, 7), (1, 4, 1), (7, 2, 3)]
    counts, witnesses = defaultdict(int), []
    for weight in prices:
        x, y, z = weight
        alpha = Fraction(1, x) + Fraction(1, z)
        kappa, gamma, q = y + 1/alpha, alpha + Fraction(1, y), alpha*y
        curvature, D = gamma*kappa, hook_step_bound(weight)
        assert curvature >= 4
        for r in (1, 2, 4, 8, 16, 32, 64, 128):
            axes = sorted([0] + [j*x for j in range(1, r + 1)] + [j*z for j in range(1, r + 1)])
            assert all(Fraction(i, 1)/alpha <= axes[i] <= Fraction(i + 1, 1)/alpha for i in range(r + 1))
            for slope, last in product((1, 3, 4, 5, 9, 12), (1, 4)):
                sizes = [last + slope*(r - s) for s in range(1, r + 1)]
                n, capacities = sum(sizes), [0]
                for size in sizes:
                    capacities.append(capacities[-1] + size)
                for tau in (0, kappa/2, kappa, (kappa + slope/gamma)/2, slope/gamma, 2*slope/gamma):
                    cost = int(r*tau)
                    ranks = [0] + subcurve_prefix_counts(r, weight, cost)
                    hooks = [b - a for a, b in zip(ranks, ranks[1:])]
                    assert all(0 <= a - b <= D for a, b in zip(hooks, hooks[1:]))
                    for s, rank in enumerate(ranks):
                        assert hook_area(weight, s, cost - 1/alpha) <= rank <= hook_area(weight, s, cost + kappa)
                        counts['area_sandwiches'] += 1
                    values = [a - b for a, b in zip(capacities, ranks)]
                    rank = n - max(values)
                    predicted = r*r*linear_profile_limit(weight, slope, Fraction(cost, r))
                    error = (abs(last - Fraction(slope, 2))*r +
                             alpha/(2*y)*(2*cost*kappa + kappa*kappa) +
                             max(Fraction(0), slope - curvature)/8)
                    assert abs(rank - predicted) <= error
                    assert generic_multicurve_prefix_rank(component_capacity_profile(sizes), weight, cost) == rank
                    counts['finite_error_rank_queries'] += 1
                    if slope > curvature and kappa*r <= cost <= slope/gamma*r:
                        center = (slope*r - gamma*cost)/(slope - curvature)
                        radius = 1 + (abs(last - curvature) + 2 + 1/q)/(slope - curvature)
                        active = [s for s in range(r + 1) if values[s] == max(values)]
                        assert all(abs(s - center) <= radius for s in active)
                        counts['active_index_localizations'] += 1
    for r in (4, 8, 12, 16, 32, 64, 128, 256):
        n, cost, parity = (9*r*r - 7*r)//2, 13*r//4, (r//4) % 2
        sizes = [1 + 9*(r - s) for s in range(1, r + 1)]
        ranks, capacity, values = subcurve_prefix_counts(r, (2, 1, 2), cost), 0, [0]
        for size, rank in zip(sizes, ranks):
            capacity += size
            values.append(capacity - rank)
        active = [s for s, v in enumerate(values) if v == max(values)]
        deficit = Fraction(5*r*r, 8) - Fraction(11*r, 4) + 2 + parity
        assert active == [r//2 - 1] and max(values) == deficit and ranks[-1] >= n
        assert generic_multicurve_prefix_rank(component_capacity_profile(sizes), (2, 1, 2), cost) == n - deficit
        witnesses.append({'components': r, 'length': n, 'cost': cost, 'rank': int(n - deficit),
                          'carrier_truncation_rank': n, 'exact_deficit': int(deficit), 'unique_active_index': active[0]})
    return {'counts': dict(counts), 'exact_linear_deficit_family': witnesses,
            'scope': 'An explicit O(r) finite error proves the linear-profile rank limit and its curvature transition at slope four. Rank loss and many cost phases do not count coherent initial ideals or establish a polynomial fan separation.'}


def steep_model_basis(r, last, base_weight):
    sizes = [last + 9*(r - s) for s in range(1, r + 1)]
    digit_base = 4*(2*r + max(sizes)) + 1
    scale = digit_base**3 + 1
    weight = tuple(scale*w + d for w, d in zip(base_weight, (digit_base**2, digit_base, 1)))
    assert max(weight) < 2*min(weight) and hook_step_bound(weight) < 9
    cutoff = (2*r + max(sizes))*max(weight)
    pool = [(a, b, c) for b in range(cutoff//weight[1] + 1)
            for a, c in ([(a, 0) for a in range((cutoff - b*weight[1])//weight[0] + 1)] +
                         [(0, c) for c in range(1, (cutoff - b*weight[1])//weight[2] + 1)])]
    assert all(max(a) < digit_base for a in pool)
    ordered = sorted(pool, key=lambda a: dot(a, weight))
    assert len({dot(a, weight) for a in ordered}) == len(ordered)
    flags = {a: [multicurve_allowed(a, s, weight) for s in range(r + 1)] for a in ordered}
    layers = [[a for a in ordered if flags[a][s] and not flags[a][s - 1]] for s in range(1, r + 1)]
    selected = {a for layer, size in zip(layers, sizes) for a in layer[:size]}
    n, capacities, counts, outside, quota_selected = sum(sizes), [0], [0]*(r + 1), [0]*(r + 1), set()
    for size in sizes:
        capacities.append(capacities[-1] + size)
    phases = [r]
    for a in ordered:
        if not flags[a][r]:
            continue
        feasible = all(outside[s] + (not flags[a][s]) <= n - capacities[s] for s in range(r + 1))
        if feasible:
            quota_selected.add(a)
            outside = [v + (not member) for v, member in zip(outside, flags[a])]
        counts = [v + member for v, member in zip(counts, flags[a])]
        values = [v - count for v, count in zip(capacities, counts)]
        active = values.index(max(values))
        if active != phases[-1]:
            assert active == phases[-1] - 1
            phases.append(active)
        if len(quota_selected) == n:
            break
    assert quota_selected == selected and phases == list(range(r, -1, -1))
    assert all(tuple(v - int(i == j) for i, v in enumerate(a)) in selected
               for a in selected for j in range(3) if a[j])
    return sizes, weight, tuple(sorted(selected)), phases


def steep_point_control(r):
    bases = [steep_model_basis(r, 1, w) for w in ((11, 13, 17), (13, 17, 11), (17, 11, 13))]
    sizes, n, prime, rng = bases[0][0], sum(bases[0][0]), 1000003, Random(20261020 + 100*r)
    for attempt in range(1, 31):
        nodes = rng.sample(range(2, prime), n)
        components, start = [], 0
        for j, size in enumerate(sizes, 1):
            components.extend((j, node) for node in nodes[start:start + size])
            start += size

        def residue(a):
            return tuple(pow(node, a[0] - a[1], prime)*pow(j, a[1], prime)*
                         pow(node - 1, -a[2], prime) % prime for j, node in components)

        columns = {a: residue(a) for _, _, selected, _ in bases for a in selected}
        try:
            inverses = [field_inverse(selected, columns, prime) for _, _, selected, _ in bases]
        except StopIteration:
            continue
        break
    else:
        raise AssertionError('no simultaneous steep-component point witness')
    certificates, border_checks = [], 0
    for (_, weight, selected, phases), inverse in zip(bases, inverses):
        border = sorted(minimal_excluded_border(selected))
        for a in border:
            column = residue(a)
            coordinates = [sum(x*y for x, y in zip(row, column)) % prime for row in inverse]
            assert all(sum(coordinates[j]*columns[b][i] for j, b in enumerate(selected)) % prime == column[i]
                       for i in range(n))
            assert all(dot(a, weight) > dot(b, weight) for b, c in zip(selected, coordinates) if c)
        border_checks += len(border)
        certificates.append({'weight': weight, 'selected': selected, 'decreasing_phases': phases,
                             'minimal_initial_generators': border})
    point_rank = None
    if r == 4:
        cost = 13
        pool = [(a, b, c) for b in range(cost + 1)
                for a, c in ([(a, 0) for a in range((cost - b)//2 + 1)] +
                             [(0, c) for c in range(1, (cost - b)//2 + 1)])]
        actual = column_prefix_ranks(pool, {a: residue(a) for a in pool}, prime)[-1]
        assert actual == 56
        point_rank = {'weight': [2, 1, 2], 'cost': cost, 'rank': actual, 'carrier_truncation_rank': n}
    start, nodes_by_component = 0, []
    for size in sizes:
        nodes_by_component.append(nodes[start:start + size])
        start += size
    return {'components': r, 'length': n, 'component_sizes': sizes, 'characteristic': prime,
            'component_scalars': list(range(1, r + 1)), 'nodes_by_component': nodes_by_component,
            'generic_sample_attempts': attempt, 'positive_basis_certificates': certificates,
            'matrix_bases': len(bases), 'minimal_border_relations': border_checks,
            'direct_prefix_rank': point_rank,
            'scope': 'Separate weights in one actual finite point ideal verify the layer decomposition. Invertible bases and strictly lower vanishing border relations prove the initials by equal colength. No complete fan or arbitrary-subset matroid identification.'}


def layer_wall_controls():
    """Exact wall curves, the arithmetic intersection filter and break counts."""
    pair_checks, segment_checks, wall_checks, fixtures = 0, 0, 0, []
    for sizes in ((10, 1), (19, 10, 1), (28, 19, 10, 1), (7, 4, 1)):
        r, L = len(sizes), max(m + s - 1 for s, m in enumerate(sizes, 1))
        curves = [(s, i, m + s - 1 - i, m) for s, m in enumerate(sizes, 1)
                  for i in range(s, m + s - 1)]

        def breaks(index):
            return {Fraction(c, a) for level in (index, index + 1)
                    for a in range(1, level) for c in (level - a,)}

        def axis_line(index, rho):
            # Independent closed two-list selection formula, checked below
            # against an explicitly sorted merge at every wall sample.
            ceil = lambda t: -((-t.numerator)//t.denominator)
            a, c = ceil(index/(1 + rho)), ceil(index*rho/(1 + rho))
            return (a, 0) if a*rho < c else (0, c)

        def line(curve, rho):
            s, i, k, _ = curve
            a, c = axis_line(i, rho)
            b, d = axis_line(s - 1, rho)
            return Fraction(a - b, k), Fraction(c - d, k)

        breakpoints = {curve: breaks(curve[1]) | breaks(curve[0] - 1) for curve in curves}
        assert all(len(breakpoints[curve]) <= 2*(curve[1] + curve[0] - 1) for curve in curves)
        candidates, actual_pairs, meetings = defaultdict(int), 0, 0
        for first, second in combinations(curves, 2):
            s, i, k, m = first
            t, j, ell, nn = second
            if s == t:
                continue
            pair_checks += 1
            necessary = abs(m*ell - nn*k) <= k + ell
            candidates[min(s, t), max(s, t)] += int(necessary)
            cuts = [Fraction(0)] + sorted(breakpoints[first] | breakpoints[second])
            found, complexity = False, 0
            for left, right in zip(cuts, cuts[1:] + [None]):
                probe = (left + right)/2 if right is not None else left + 1
                a, b = line(first, probe)
                c, d = line(second, probe)
                if a != c:
                    root = (d - b)/(a - c)
                    if root > 0 and left <= root and (right is None or root <= right):
                        assert necessary
                        found, complexity = True, complexity + 1
                elif b == d:
                    assert necessary
                    found, complexity = True, complexity + 1
                segment_checks += 1
            if found:
                actual_pairs += 1
                meetings += complexity
            assert complexity <= 1 + len(breakpoints[first]) + len(breakpoints[second])
        for (s, t), count in candidates.items():
            assert count <= 3*sizes[s - 1] + 2*sizes[t - 1]
        # Check the exact quota at wall points independently of the area filter.
        for curve in curves:
            cuts = [Fraction(0)] + sorted(breakpoints[curve])
            for left, right in zip(cuts, cuts[1:] + [None]):
                rho = (left + right)/2 if right is not None else left + 1
                a, b = line(curve, rho)
                eta = a*rho + b
                if eta <= 0:
                    continue
                s, i, k, m = curve
                assert i - s + k == m - 1
                axis = sorted([Fraction(0)] + [j*rho for j in range(1, L + 1)] +
                              [Fraction(j) for j in range(1, L + 1)])
                threshold = axis[i] - axis[s - 1]
                assert threshold == k*eta
                strict = (sum(v - axis[s - 1] < threshold for v in axis[s - 1:]) +
                          sum(j*eta < threshold for j in range(1, m + 1)))
                assert strict == m - 1
                wall_checks += 1
        fixtures.append({'sizes': sizes, 'curve_count': len(curves), 'maximum_axis_index': L,
                         'arithmetically_possible_curve_pairs': sum(candidates.values()),
                         'meeting_curve_pairs_including_overlaps': actual_pairs,
                         'piecewise_meeting_intervals_including_shared_endpoints': meetings})
    return {'curve_pair_checks': pair_checks, 'affine_piece_pair_checks': segment_checks,
            'exact_wall_quota_checks': wall_checks, 'fixtures': fixtures,
            'scope': 'All positive wall curves of formal independent layer selection are checked. The uniform arithmetic filter and piece count prove an upper bound; these fixtures are not full point-ideal fan counts.'}


def steep_component_controls():
    models = [steep_model_basis(r, last, weight) for r in range(2, 9) for last in (1, 3)
              for weight in ((11, 13, 17), (11, 17, 13), (13, 11, 17),
                             (13, 17, 11), (17, 11, 13), (17, 13, 11))]
    points = [steep_point_control(r) for r in (4, 5)]
    return {'rank_limit': steep_rank_controls(), 'layer_wall_counting': layer_wall_controls(),
            'fixed_price_model_bases': len(models),
            'singleton_phase_blocks': sum(len(model[3]) - 1 for model in models),
            'point_ideals': points,
            'scope': 'Sufficiently steep gaps give every singleton phase and an exact positive layer sum. An arithmetic wall-curve argument bounds all all-tight selections by O(n^2) for fixed linear profiles, including all axis orders and aspects. Pooled phases and unrestricted G_3 remain open.'}


def pooled_block_controls():
    """Complete one-parameter trace and short-block exchange certificates."""
    r, slope, H, R = 6, 5, 2, 3
    sizes = [1 + slope*(r - s) for s in range(1, r + 1)]
    n, capacities = sum(sizes), [0]
    for size in sizes:
        capacities.append(capacities[-1] + size)
    axes = sorted([(0, (0, 0, 0))] + [(101*j, (j, 0, 0)) for j in range(1, n)] +
                  [(103*j, (0, 0, j)) for j in range(1, n)])[:n]
    assert len({v for v, _ in axes}) == n
    grid = [(i, b) for i in range(n) for b in range(n)]
    differences = {v - u for u, _ in axes for v, _ in axes if v > u}
    events = sorted({Fraction(48), Fraction(56)} |
                    {Fraction(d, k) for d in differences for k in range(1, n) if 48*k < d < 56*k})
    blocks = [(t, s, capacities[s] - capacities[t]) for t in range(r)
              for s in range(t + 1, min(r, t + 2*H) + 1)]
    counters, signatures, candidates, phase_patterns = defaultdict(int), {}, {}, set()
    observed_curves, observed_lines = set(), set()
    max_width, previous, previous_phase = 0, None, None

    def ordered_grid(y):
        keys = {p: axes[p[0]][0]*y.denominator + p[1]*y.numerator for p in grid}
        return sorted(grid, key=keys.__getitem__), keys

    def block_head(ordered, keys, t, s, take):
        head = []
        cutoff = None
        for p in ordered:
            if t <= min(p) < s:
                if cutoff is not None and keys[p] > cutoff:
                    break
                head.append(p)
                if len(head) == take:
                    cutoff = keys[p]
        assert len(head) >= take
        return head

    for left, right in zip(events, events[1:]):
        y = (left + right)/2
        ordered, keys = ordered_grid(y)
        assert len(set(keys.values())) == len(grid)
        counts, outside, chosen, phases = [0]*(r + 1), [0]*(r + 1), [], [r]
        for p in ordered:
            layer = min(p)
            if layer >= r:
                continue
            expected = layer < phases[-1]
            feasible = all(outside[s] + int(layer >= s) <= n - capacities[s] for s in range(r + 1))
            assert expected == feasible
            if feasible:
                chosen.append(p)
                outside = [v + int(layer >= s) for s, v in enumerate(outside)]
            counts = [v + int(layer < s) for s, v in enumerate(counts)]
            values = [a - d for a, d in zip(capacities, counts)]
            active = values.index(max(values))
            if active != phases[-1]:
                phases.append(active)
            if len(chosen) == n:
                break
        assert len(chosen) == n and phases[-1] == 0
        selected = tuple(sorted(chosen))
        width = max(s - t for s, t in zip(phases, phases[1:]))
        max_width = max(max_width, width)
        cutoffs = []
        for s, t in zip(phases, phases[1:]):
            take = capacities[s] - capacities[t]
            prefix = block_head(ordered, keys, t, s, take)[:take]
            assert set(prefix) == {p for p in chosen if t <= min(p) < s}
            cutoffs.append(max(keys[p] for p in prefix))
            counters['canonical_block_prefix_checks'] += 1
        assert all(a < b for a, b in zip(cutoffs, cutoffs[1:]))
        signature = []
        for t, s, take in blocks:
            prefix = block_head(ordered, keys, t, s, take)[:take]
            signature.extend((tuple(sorted(prefix)), tuple(sorted(prefix[:-1]))))
        if width <= H:
            signature = tuple(signature)
            assert signatures.setdefault(signature, selected) == selected
            counters['bounded_width_intervals'] += 1
            if width > 1:
                counters['genuinely_pooled_intervals'] += 1
                candidates.setdefault(selected, (y, phases))
        if previous is not None and previous != selected:
            counters['actual_basis_transitions'] += 1
        if previous_phase is not None and previous_phase != phases:
            counters['phase_partition_transitions'] += 1
        previous, previous_phase = selected, phases
        phase_patterns.add(tuple(phases))
        counters['complete_cost_order_intervals'] += 1
    assert max_width <= H and counters['genuinely_pooled_intervals'] > 0
    assert counters['phase_partition_transitions'] > 0
    assert counters['actual_basis_transitions'] > 0

    for y in events[1:-1]:
        ordered, keys = ordered_grid(y)
        at_event = []
        for t, s, take in blocks:
            head = block_head(ordered, keys, t, s, take)
            h = s - t
            for target in (take, take - 1):
                if not target:
                    continue
                cost = keys[head[target - 1]]
                strict = sum(keys[p] < cost for p in head)
                tied = [p for p in head if keys[p] == cost]
                if not strict < target < strict + len(tied):
                    continue
                assert len(tied) <= 3*h
                for first, second in combinations(tied, 2):
                    if (first[0] < s) != (second[0] < s):
                        col, row = (first, second) if first[0] < s else (second, first)
                        i, b1 = col
                        j, b0 = row
                        k = b1 - b0
                        assert k >= 1 and j >= s
                        N = j + k
                        assert abs(Fraction(N - s) - Fraction(target, h)) <= h + 3*h*R + 4
                        observed_curves.add((i, j, k))
                        at_event.append((N - i, k))
                        counters['cross_row_column_walls'] += 1
                    else:
                        d = first[1] - second[1]
                        assert d and abs(d) <= h*R
                        a = axes[first[0]][1][0] - axes[second[0]][1][0]
                        c = axes[second[0]][1][2] - axes[first[0]][1][2]
                        if a < 0 or (not a and c < 0):
                            a, c, d = -a, -c, -d
                        assert a >= 0 and c >= 0
                        assert (a*101 - c*103)*y.denominator + d*y.numerator == 0
                        observed_lines.add((a, c, d))
                        counters['small_y_coefficient_walls'] += 1
        for (mu, k), (nu, ell) in combinations(at_event, 2):
            assert abs(mu*ell - nu*k) <= k + ell
            counters['generalized_curve_intersection_filters'] += 1

    def axis_line(index, rho):
        ceil = lambda t: -((-t.numerator)//t.denominator)
        a, c = ceil(index/(1 + rho)), ceil(index*rho/(1 + rho))
        return (a, 0) if a*rho < c else (0, c)

    Q = max(j for _, j, _ in observed_curves)
    # The narrow trace need not encounter a line steep enough for the
    # monotonicity branch; include explicit admissible small-d lines too.
    test_lines = observed_lines | {(a, c, d) for a in range(max(1, Q - 2), Q + 1)
                                   for c in range(max(1, Q - 2), Q + 1) for d in (-1, 1)}
    for (i, j, k), (a, c, d) in product(observed_curves, test_lines):
        if a*k <= 2*abs(d)*Q:
            continue
        breaks = {Fraction(1, R), Fraction(R)}
        for index in (i, j):
            breaks.update(Fraction(v, u) for level in (index, index + 1)
                          for u in range(1, level) for v in (level - u,)
                          if Fraction(1, R) < Fraction(v, u) < R)
        roots = set()
        cuts = sorted(breaks)
        for left, right in zip(cuts, cuts[1:]):
            rho = (left + right)/2
            x1, z1 = axis_line(i, rho)
            x2, z2 = axis_line(j, rho)
            slope_difference = Fraction(x2 - x1, k) + Fraction(a, d)
            intercept_difference = Fraction(z2 - z1, k) - Fraction(c, d)
            assert slope_difference
            root = -intercept_difference/slope_difference
            if left <= root <= right:
                roots.add(root)
        assert len(roots) <= 1
        counters['steep_line_curve_intersection_checks'] += 1
    assert counters['steep_line_curve_intersection_checks'] > 0

    # Global logarithmic-coordinate branch: the product of the two line
    # coefficients replaces a compact aspect bound. No logarithms are
    # evaluated; roots and the normalized derivatives are rational.
    for (i, j, k), (a, c, d) in product(observed_curves, test_lines):
        if a*c*k*k <= d*d*(Q + 1)**2:
            continue
        cuts = [Fraction(0)] + sorted({Fraction(v, u) for index in (i, j)
                                       for level in (index, index + 1)
                                       for u in range(1, level) for v in (level - u,)})
        roots = set()
        for left, right in zip(cuts, cuts[1:] + [None]):
            rho = (left + right)/2 if right is not None else left + 1
            x1, z1 = axis_line(i, rho)
            x2, z2 = axis_line(j, rho)
            first_derivative = x1*rho if x1 else -Fraction(z1, 1)/rho
            second_derivative = x2*rho if x2 else -Fraction(z2, 1)/rho
            assert abs(first_derivative) <= i + 1 and abs(second_derivative) <= j + 1
            assert abs((second_derivative - first_derivative)/k) <= Fraction(2*(Q + 1), k)
            slope_difference = Fraction(x2 - x1, k) + Fraction(a, d)
            intercept_difference = Fraction(z2 - z1, k) - Fraction(c, d)
            if not slope_difference:
                assert intercept_difference
                continue
            root = -intercept_difference/slope_difference
            if root > 0 and left <= root and (right is None or root <= right):
                roots.add(root)
        assert len(roots) <= 1
        counters['global_logarithmic_line_curve_checks'] += 1
    assert counters['global_logarithmic_line_curve_checks'] > 0
    for limit in (4, 8, 16, 32):
        harmonic = sum(Fraction(1, a) for a in range(1, limit + 1))
        for k, d in product(range(1, limit + 1), range(1, 5)):
            threshold = Fraction(d*d*(limit + 1)**2, k*k)
            direct = sum(a*c <= threshold for a in range(1, limit + 1) for c in range(1, limit + 1))
            assert direct <= threshold*harmonic
            counters['integer_product_pair_counts'] += 1

    chosen_cases = list(candidates.items())
    assert len(chosen_cases) >= 3
    cases = [chosen_cases[j] for j in sorted({0, len(chosen_cases)//2, len(chosen_cases) - 1})]
    prime, rng = 1000003, Random(20261020 + 4300)

    def exponent(p):
        a = axes[p[0]][1]
        return (a[0], p[1], a[2])

    for attempt in range(1, 31):
        nodes = rng.sample(range(2, prime), n)
        components, start = [], 0
        for j, size in enumerate(sizes, 1):
            components.extend((j, node) for node in nodes[start:start + size])
            start += size

        def residue(a):
            return tuple(pow(node, a[0] - a[1], prime)*pow(j, a[1], prime)*
                         pow(node - 1, -a[2], prime) % prime for j, node in components)

        bases = [tuple(sorted(exponent(p) for p in selected)) for selected, _ in cases]
        columns = {a: residue(a) for selected in bases for a in selected}
        try:
            inverses = [field_inverse(selected, columns, prime) for selected in bases]
        except StopIteration:
            continue
        break
    else:
        raise AssertionError('no simultaneous pooled point witness')
    certificates = []
    for selected, inverse, (_, (y, phases)) in zip(bases, inverses, cases):
        weight = (101*y.denominator, y.numerator, 103*y.denominator)
        border = sorted(minimal_excluded_border(selected))
        for a in border:
            column = residue(a)
            coordinates = [sum(x*z for x, z in zip(row, column)) % prime for row in inverse]
            assert all(sum(coordinates[j]*columns[b][i] for j, b in enumerate(selected)) % prime == column[i]
                       for i in range(n))
            assert all(dot(a, weight) > dot(b, weight) for b, value in zip(selected, coordinates) if value)
        counters['point_matrix_bases'] += 1
        counters['point_border_relations'] += len(border)
        certificates.append({'weight': weight, 'selected': selected, 'canonical_phases': phases,
                             'minimal_initial_generators': border})
    start, nodes_by_component = 0, []
    for size in sizes:
        nodes_by_component.append(nodes[start:start + size])
        start += size
    return {'counts': dict(counters), 'component_sizes': sizes, 'length': n,
            'fixed_axis_prices': [101, 103], 'y_price_interval': [48, 56],
            'complete_breakpoints': [str(y) for y in events],
            'maximum_canonical_block_width': max_width,
            'certificate_block_width': 2*H, 'aspect_bound': R,
            'distinct_short_block_signatures': len(signatures),
            'distinct_actual_bases': len(set(signatures.values())),
            'distinct_phase_patterns': sorted(phase_patterns),
            'point_ideal': {'characteristic': prime, 'component_scalars': list(range(1, r + 1)),
                            'nodes_by_component': nodes_by_component, 'generic_sample_attempts': attempt,
                            'positive_basis_certificates': certificates},
            'scope': 'Complete one-dimensional trace, not a whole fan. Short-block rank-q and rank-(q-1) signatures, critical walls, compact and global logarithmic line-curve branches and integer-product counts are checked. Three actual pooled bases have exact point certificates. Uniform bounded-width theorems are proved separately.'}


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
    charging = hook_charging_controls()
    unequal = {'phase_controls': unequal_phase_controls(),
               'hook_capacity_compression': hook_capacity_jump_controls(),
               'simultaneous_two_cost_point_faces': [unequal_point_controls(t) for t in (1, 2)]}
    steep = steep_component_controls()
    pooled = pooled_block_controls()
    paths = ['scripts/run_algebraic_incidence_controls.py', 'docs/weighted_facet_incidence.md',
             'docs/initial_cone_height.md', 'scripts/run_extremal_initial_ideals_checks.py',
             'scripts/run_structural_complexity_checks.py', 'scripts/run_initial_cone_height_checks.py',
             'scripts/run_interleaving_research.py', 'docs/extremal_initial_ideals.md',
             'docs/generic_multiaffine_bound.md', 'docs/planar_dyadic_bound.md',
             'docs/rational_shadow_bound.md', 'scripts/run_planar_dyadic_bound_checks.py']
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_algebraic_incidence_controls.py',
              'scope': 'Short-block certificates give O_(H,R)(L^4) states on compact aspects and O_(H,u,v)(L^4 log(L+1)) on a bounded harmonic-ratio window, without aspect bounds. Fixed supercritical margins of linear profiles therefore have quadratic-power counts globally. Growing widths, critical/subcritical regimes, unrestricted G_3 and exact priority remain unresolved.',
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
                          **multicurves['rank_controls'],
                          'hook_charged_planar_groups': len(charging['planar_tie_groups']),
                          'hook_transverse_rank_samples': sum(a['complete_transverse_samples'] for a in charging['planar_tie_groups']),
                          'hook_charged_spatial_vertices': sum(len(a['off_carrier_vertices']) for a in charging['planar_tie_groups']),
                          'hook_carrier_trace_intervals': sum(a['complete_comparison_intervals'] for a in charging['complete_carrier_traces']),
                          'hook_active_local_bases': charging['active_local_point_quotient']['local_exposed_bases'],
                          'hook_active_matrix_bases': charging['active_local_point_quotient']['all_local_and_trace_basis_rank_checks'],
                          'hook_active_s_pairs': charging['active_local_point_quotient']['ordinary_s_pairs'],
                          'hook_active_normal_forms': charging['active_local_point_quotient']['ordinary_normal_forms'],
                          **{'unequal_' + k: v for k, v in unequal['phase_controls']['counts'].items()},
                          **{'unequal_' + k: v for k, v in unequal['hook_capacity_compression']['counts'].items()},
                          'unequal_point_matrix_bases': sum(a['matrix_bases'] for a in unequal['simultaneous_two_cost_point_faces']),
                          'unequal_point_border_relations': sum(a['minimal_border_relations'] for a in unequal['simultaneous_two_cost_point_faces']),
                          'unequal_direct_prefix_rank_queries': sum(len(a['direct_prefix_ranks']) for a in unequal['simultaneous_two_cost_point_faces']),
                          **{'steep_' + k: v for k, v in steep['rank_limit']['counts'].items()},
                          'steep_model_bases': steep['fixed_price_model_bases'],
                          'steep_singleton_phase_blocks': steep['singleton_phase_blocks'],
                          'steep_point_matrix_bases': sum(a['matrix_bases'] for a in steep['point_ideals']),
                          'steep_point_border_relations': sum(a['minimal_border_relations'] for a in steep['point_ideals']),
                          **{'steep_' + k: steep['layer_wall_counting'][k] for k in
                             ('curve_pair_checks', 'affine_piece_pair_checks', 'exact_wall_quota_checks')},
                          **{'pooled_' + k: v for k, v in pooled['counts'].items()}},
              'uniform_axis_matroids': axes, 'arithmetic_lower_bound_controls': arithmetic,
              'interior_quotients': interior, 'three_pole_quotients': curves,
              'positive_plane_traces': planes,
              'generic_two_curve_unions': unions,
              'generic_multicurve_unions': multicurves,
              'balanced_hook_threshold_charging': charging,
              'unequal_component_phases': unequal,
              'steep_component_profiles': steep,
              'pooled_component_blocks': pooled,
              'source_sha256': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}}
    (ROOT / 'results/research/algebraic_incidence_controls.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report['summary'], indent=2))


if __name__ == '__main__':
    main()
