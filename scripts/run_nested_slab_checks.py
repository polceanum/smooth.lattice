#!/usr/bin/env python3
"""Exact slab classification, prefix-pair, and large-footprint controls.

All uniform proofs and the classical geometric boundary are in
docs/nested_jet_matroids.md. No asymptotic growth is fitted to these cases.
"""
from fractions import Fraction
from itertools import combinations, product
from math import comb, gcd, prod
from pathlib import Path
import hashlib
import json

from run_extremal_initial_ideals_checks import (
    adjacent_weights, cone_for_basis, integer_weight, inverse_mod,
)
from run_initial_cone_height_checks import add, clean, leading, normal_form, subtract_shift
from run_interleaving_research import minimal_excluded_border
from run_nested_jet_checks import (
    border_basis, downsets, independent, model, rank_columns,
    schedule_greedy, valuation_rank,
)
from run_structural_complexity_checks import area2, clip, dot, greedy, hyperbolic_cross, primitive

ROOT = Path(__file__).resolve().parents[1]


def s_pair_checks(polynomials, weight, characteristic):
    heads = [leading(g, weight) for g in polynomials]
    for i, j in combinations(range(len(polynomials)), 2):
        common = tuple(max(x, y) for x, y in zip(heads[i], heads[j]))
        first = {add(a, tuple(x - y for x, y in zip(common, heads[i]))): c
                 for a, c in polynomials[i].items()}
        pair = subtract_shift(first, polynomials[j],
                              tuple(x - y for x, y in zip(common, heads[j])), 1, characteristic)
        assert not normal_form(pair, polynomials, weight, characteristic)
    return comb(len(polynomials), 2)


def planar_flags(n):
    ground = tuple(hyperbolic_cross(n + 1, 2))
    boundaries = sorted({Fraction(b[1] - a[1], a[0] - b[0])
                         for a, b in combinations(ground, 2)
                         if (a[0] - b[0]) * (a[1] - b[1]) < 0})
    intervals = list(zip([Fraction(0)] + boundaries, boundaries + [None]))
    words = []
    for left, right in intervals:
        t = left + 1 if right is None else (left + right) / 2
        ordered = tuple(sorted(ground, key=lambda a: dot(a, (t, 1)))[:n])
        assert len({dot(a, (t, 1)) for a in ground}) == len(ground)
        words.append((t, ordered))
    walls = {r: set() for r in range(n + 1)}
    for i, boundary in enumerate(boundaries):
        for r in range(1, n + 1):
            if set(words[i][1][:r]) != set(words[i + 1][1][:r]):
                walls[r].add(boundary)
    expected, bases, samples = 0, set(), []
    for b in range(n // 2 + 1):
        a = n - b
        pairs = {(frozenset(word[:a]), frozenset(word[:b])) for _, word in words}
        assert len(pairs) == 1 + len(walls[a] | walls[b])
        expected += len(pairs)
        bases |= {frozenset({u + (0,) for u in first} | {u + (1,) for u in second})
                  for first, second in pairs}
        samples.append({'upper_layer_count': b, 'common_prefix_pairs': len(pairs),
                        'lower_rank_walls': len(walls[a]), 'upper_rank_walls': len(walls[b]),
                        'shared_walls': len(walls[a] & walls[b])})
    assert len(bases) == expected
    return expected, bases, {'length': n, 'full_planar_comparison_intervals': len(words),
                            'two_layer_formula': expected, 'prefix_pair_counts': samples}


def point_fixture(n, h, characteristic=1000003):
    q, r = divmod(n, h)
    levels = [i for i in range(h) for _ in range(q + int(i < r))]
    points = [(pow(2, j + 1, characteristic), pow(3, (j + 1)**2 + j + 2, characteristic), z)
              for j, z in enumerate(levels)]
    assert len(points) == len(set(points)) == n

    def residue(a):
        return tuple(prod(pow(x, b, characteristic) for x, b in zip(point, a)) % characteristic
                     for point in points)

    states = [D for D in downsets(n, 3) if max(a[2] for a in D) < h]
    for D in states:
        assert inverse_mod([[residue(a)[i] for a in sorted(D)] for i in range(n)],
                           characteristic) is not None
    # These are integer coordinate lifts. A nonzero modular determinant also
    # certifies the corresponding integer determinant is nonzero over Q.
    return points, residue, len(states)


def modular_border(basis, residue, characteristic):
    basis = sorted(basis)
    rows = inverse_mod([[residue(a)[i] for a in basis] for i in range(len(basis))], characteristic)
    assert rows is not None
    result = []
    for u in minimal_excluded_border(basis):
        column = residue(u)
        coefficients = [sum(x * y for x, y in zip(row, column)) % characteristic for row in rows]
        result.append(clean({u: 1, **{a: -c for a, c in zip(basis, coefficients) if c}}, characteristic))
    return result


def slab_fan(n, h):
    points, residue, downset_count = point_fixture(n, h)
    prime = 1000003
    ground = tuple(a for a in hyperbolic_cross(n + 1, 3) if a[2] < h)
    columns = {a: residue(a) for a in ground}
    comparisons = sorted({primitive(tuple(x - y for x, y in zip(a, b)))
                          for a, b in combinations(ground, 2)
                          if min(x - y for x, y in zip(a, b)) < 0 < max(x - y for x, y in zip(a, b))})
    pending, visited = [(1, n + 2, (n + 2)**2)], {}
    s_pairs, matrix_cones, box_certificates = 0, 0, []
    while pending:
        w = pending.pop()
        ordered = sorted(ground, key=lambda a: (dot(a, w), a))
        D = frozenset(ordered[:n])
        if D in visited:
            continue
        assert greedy(ground, columns, w, prime) == set(D)
        border = [u for u in minimal_excluded_border(sorted(D)) if u[2] < h]
        normals = sorted({tuple(x // gcd(*difference) for x in difference)
                          for u in border for v in D
                          for difference in [tuple(a - b for a, b in zip(u, v))]})
        assert all(dot(w, normal) > 0 for normal in normals)
        polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
        for normal in normals:
            polygon = clip(polygon, normal, 1)
        assert polygon
        matrix_polygon, _ = cone_for_basis(D, residue, prime)
        assert set(matrix_polygon) == set(polygon)
        matrix_cones += 1
        if n <= 6:
            polynomials = modular_border(D, residue, prime)
            assert set(leading(g, w) for g in polynomials) == set(minimal_excluded_border(sorted(D)))
            s_pairs += s_pair_checks(polynomials, w, prime)
        widths = tuple(1 + max(a[i] for a in D) for i in range(3))
        assert n <= prod(widths) <= 27 * n
        inner = tuple((width + 2) // 3 for width in widths)
        assert all(a in D for a in product(*(range(side) for side in inner)))
        dyadic = tuple((side - 1).bit_length() for side in widths)
        assert n <= 2**sum(dyadic) <= 216 * n
        box_certificates.append({'axis_counts': widths, 'dyadic_exponents': dyadic})
        visited[D] = (polygon, w)
        pending.extend(adjacent_weights(polygon, normals, comparisons))
    assert sum(area2(polygon) for polygon, _ in visited.values()) == 1
    formula = None
    if h == 2:
        expected, pair_bases, formula = planar_flags(n)
        assert expected == len(visited) and pair_bases == set(visited)
    return {'length': n, 'slab_height': h, 'complete_positive_fan': len(visited),
            'complete_simplex_double_area': '1', 'point_field_characteristic': prime,
            'integer_point_lifts_also_valid_over_Q': True, 'points': points,
            'every_slab_downset_matrix_bases': downset_count,
            'independent_matrix_cones': matrix_cones, 'ordinary_s_pairs': s_pairs,
            'two_layer_prefix_formula': formula,
            'box_certificates': box_certificates,
            'witnesses': [{'standard_set': sorted(D), 'strict_weight': w,
                           'closed_cell': [[str(x) for x in p] for p in polygon]}
                          for D, (polygon, w) in sorted(visited.items(), key=lambda item: sorted(item[0]))]}


def balanced_jet_controls():
    fixtures = []
    for h, q in ((1, 4), (2, 2), (2, 3), (3, 2), (4, 2)):
        n = h * q
        data = model(n, (1, 1, q))
        states = [D for D in downsets(n, 3) if max(a[2] for a in D) < h]
        for D in states:
            sizes = [sum(a[2] == c for a in D) for c in range(h)]
            assert all(sum(sizes[:c]) >= q * c for c in range(h))
            index = 0
            for c in range(h):
                layer = sorted((a for a in D if a[2] == c), key=lambda a: (a[0] + a[1], a))
                for a in layer:
                    assert dot(a, (1, 1, q)) <= index
                    index += 1
            assert index == n and independent([dot(a, (1, 1, q)) for a in D], n)
            assert rank_columns([data['column'](a) for a in D]) == n
        fixtures.append({'length': n, 'height': h, 'q': q,
                         'every_slab_downset_rational_basis': len(states),
                         'explicit_valuation_index_certificates': len(states)})
    return fixtures


def general_controls():
    data = model(8, (1, 3, 4))
    w = (3, 2, 1)
    D = schedule_greedy(data['ground'], w, (1, 3, 4), 8)
    expected = {(0, 0, 0), (1, 0, 0), (2, 0, 0), (0, 1, 0),
                (0, 2, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)}
    assert D == expected and greedy(data['ground'], {a: data['column'](a) for a in data['ground']}, w, 0) == expected
    polynomials, normals = border_basis(D, data['residue'])
    assert min(dot(w, normal) for normal in normals) == 1
    pair_count = s_pair_checks(polynomials, w, 0)
    assert (2, 0, 0) in D and (0, 2, 0) in D and (1, 1, 0) not in D
    footprint = []
    for M in (2, 3, 4, 6, 8, 12):
        n, m = M * M, (1, M, M)
        ground = tuple((a, b, c) for b in range(M) for c in range(M - b)
                       for a in range(n - M * (b + c)))
        weight = ((M + 1)**3, M + 1, M + 2)
        D = schedule_greedy(ground, weight, m, n)
        assert len(D) == n
        widths = tuple(1 + max(a[i] for a in D) for i in range(3))
        assert widths == (M, M, M)
        assert all((a, 0, 0) in D for a in range(M))
        assert all((0, b, c) in D for b in range(M) for c in range(M - b))
        defects = [max(M * d - (r + 1) * comb(d + 1, 2) for d in range(M + 1))
                   for r in range(M)]
        ranks = [n - delta for delta in defects]
        assert ranks[0] == comb(M + 1, 2) and ranks[-1] == n
        layer_sizes = [sum(a[0] == r for a in D) for r in range(M)]
        assert layer_sizes == [ranks[0]] + [ranks[i] - ranks[i - 1] for i in range(1, M)]
        for r in range(M):
            assert valuation_rank([dot(a, m) for a in ground if a[0] <= r], n) == ranks[r]
        matrix_layers = 0
        if M <= 4:
            actual = model(n, m)
            columns = {a: actual['column'](a) for a in ground}
            assert greedy(ground, columns, weight, 0) == set(D)
            for r in range(M):
                assert rank_columns([column for a, column in columns.items() if a[0] <= r]) == ranks[r]
                matrix_layers += 1
            _, inequalities = border_basis(D, actual['residue'])
            assert all(dot(weight, normal) > 0 for normal in inequalities)
        footprint.append({'M': M, 'length': n, 'valuations': m, 'strict_weight': weight,
                          'axis_counts': widths, 'box_lattice_points': M**3,
                          'box_to_length_ratio': M, 'exact_layer_defects': defects,
                          'selected_layer_sizes': layer_sizes, 'independent_matrix_layer_ranks': matrix_layers,
                          'standard_set': sorted(D)})
    return {'nonconvex_fixture': {'length': 8, 'valuations': [1, 3, 4], 'strict_weight': w,
                                 'standard_set': sorted(expected), 'included_endpoints': [(2, 0, 0), (0, 2, 0)],
                                 'excluded_midpoint': [1, 1, 0], 'minimum_border_margin': 1,
                                 'ordinary_s_pairs': pair_count},
            'large_footprint_family': footprint}


def main():
    balanced = balanced_jet_controls()
    fans = [slab_fan(n, h) for n, h in ((4, 2), (5, 2), (6, 2), (8, 2),
                                       (6, 3), (8, 4), (7, 7), (6, 1))]
    general = general_controls()
    paths = ['scripts/run_nested_slab_checks.py', 'docs/nested_jet_matroids.md',
             'scripts/run_nested_jet_checks.py', 'scripts/run_extremal_initial_ideals_checks.py',
             'scripts/run_initial_cone_height_checks.py', 'scripts/run_interleaving_research.py',
             'scripts/run_structural_complexity_checks.py']
    summary = {'balanced_jet_all_downset_bases': sum(a['every_slab_downset_rational_basis'] for a in balanced),
               'radical_point_all_downset_bases': sum(a['every_slab_downset_matrix_bases'] for a in fans),
               'complete_slab_fans': len(fans), 'slab_initial_ideals': sum(a['complete_positive_fan'] for a in fans),
               'independent_two_layer_formulas': sum(a['slab_height'] == 2 for a in fans),
               'ordinary_s_pairs': sum(a['ordinary_s_pairs'] for a in fans) + general['nonconvex_fixture']['ordinary_s_pairs'],
               'large_footprint_fixtures': len(general['large_footprint_family']),
               'independent_large_footprint_layer_ranks': sum(a['independent_matrix_layer_ranks'] for a in general['large_footprint_family'])}
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_nested_slab_checks.py',
              'scope': 'Uniform slab classification and fixed-height asymptotic count have explicit proofs, with their classical geometric ingredients identified. General nested valuations admit nonconvex sets and coherent boxes of n^(3/2) points. Neither fixture counts nor footprint sizes give a spatial power lower bound or a novelty certificate.',
              'summary': summary, 'balanced_jet_controls': balanced, 'slab_fans': fans,
              'general_valuation_controls': general,
              'source_sha256': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
              'new_unrestricted_counting_bound': False, 'novelty_established': False, 'goal_complete': False}
    destination = ROOT / 'results/research/nested_slab_checks.json'
    destination.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': 'passed', **summary,
                      'slab_counts': [[a['length'], a['slab_height'], a['complete_positive_fan']] for a in fans]}, indent=2))


if __name__ == '__main__':
    main()
