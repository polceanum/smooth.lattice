#!/usr/bin/env python3
"""Exact controls for polynomial quotients with a nested rank rule.

The construction and proof are in docs/nested_jet_matroids.md. These
checks do not establish an asymptotic fan count or a novelty claim.
"""
from fractions import Fraction
from itertools import combinations, permutations, product
from math import factorial, gcd, prod
from pathlib import Path
import hashlib
import json

from run_extremal_initial_ideals_checks import adjacent_weights, integer_weight
from run_initial_cone_height_checks import add, clean, leading, normal_form, subtract_shift
from run_interleaving_research import minimal_excluded_border
from run_structural_complexity_checks import (
    area2, clip, dot, greedy, hyperbolic_cross, positive_chambers, primitive,
)

ROOT = Path(__file__).resolve().parents[1]


def rank_columns(columns):
    """Rational Gaussian elimination, independent of the valuation rule."""
    pivots = {}
    for column in columns:
        v = list(map(Fraction, column))
        for index, row in sorted(pivots.items()):
            c = v[index]
            v = [x - c * y for x, y in zip(v, row)]
        index = next((i for i, x in enumerate(v) if x), None)
        if index is not None:
            c = v[index]
            pivots[index] = [x / c for x in v]
    return len(pivots)


def valuation_rank(values, n):
    return min([len(values)] + [n - k + sum(v < k for v in values)
                               for k in range(n + 1)])


def independent(values, n):
    return all(v <= n - len(values) + i for i, v in enumerate(sorted(values)))


def product_truncated(a, b):
    n = len(a)
    return tuple(sum(a[j] * b[i - j] for j in range(i + 1)) for i in range(n))


def model(n, multiplicities):
    """n is the ambient truncation; the image length may be smaller."""
    assert n >= 1 and min(multiplicities) >= 1
    B = max(n, 2)
    digits = tuple(B**i for i in range(3))
    C = n * sum(digits)
    lambdas = tuple(C * m + d for m, d in zip(multiplicities, digits))
    zero = (Fraction(0),) * n
    one = (Fraction(1),) + zero[1:]

    def exponential_column(a):
        v, node = dot(a, multiplicities), dot(a, lambdas)
        return tuple(Fraction(node**(j - v), factorial(j - v)) if j >= v else Fraction(0)
                     for j in range(n))

    generators = [exponential_column(tuple(int(i == j) for i in range(3)))
                  for j in range(3)]

    def residue(a):
        answer = one
        for generator, power in zip(generators, a):
            for _ in range(power):
                answer = product_truncated(answer, generator)
        return answer

    full_ground = tuple(a for a in product(range(n), repeat=3)
                        if dot(a, multiplicities) < n)
    ordered = sorted(full_ground, key=lambda a: dot(a, lambdas))
    assert len({dot(a, lambdas) for a in full_ground}) == len(full_ground)
    assert [dot(a, multiplicities) for a in ordered] == sorted(dot(a, multiplicities)
                                                                            for a in ordered)
    assert all(exponential_column(a) == residue(a) for a in full_ground)
    length = valuation_rank([dot(a, multiplicities) for a in full_ground], n)
    assert rank_columns([exponential_column(a) for a in full_ground]) == length
    ideal_square = [dot(a, multiplicities) for a in full_ground if sum(a) >= 2]
    embedding_dimension = length - 1 - valuation_rank(ideal_square, n)
    assert embedding_dimension == length - 1 - rank_columns(
        [exponential_column(a) for a in full_ground if sum(a) >= 2])
    relation_coefficients = []
    if multiplicities[0] == 1:
        assert length == n
        assert rank_columns([residue((j, 0, 0)) for j in range(n)]) == n
        # Check the explicit curvilinear ideal by independent composition.
        powers = [one]
        for _ in range(n):
            powers.append(product_truncated(powers[-1], generators[0]))
        for m, lam, target in zip(multiplicities[1:], lambdas[1:], generators[1:]):
            coefficients = [Fraction(0)] * n
            if m < n:
                coefficients[m] = Fraction(1)
            for j in range(m + 1, n):
                coefficients[j] = Fraction((lam - m * lambdas[0]) *
                                           (lam - j * lambdas[0])**(j - m - 1), factorial(j - m))
            replay = tuple(sum(c * column[i] for c, column in zip(coefficients, powers))
                           for i in range(n))
            assert replay == target
            relation_coefficients.append([str(c) for c in coefficients])
        assert powers[n] == zero
    # Direct multiplication operators commute on the entire reference basis.
    commutators = 0
    for f, g in combinations(generators, 2):
        for j in range(n):
            e = tuple(Fraction(i == j) for i in range(n))
            assert product_truncated(f, product_truncated(g, e)) == \
                product_truncated(g, product_truncated(f, e))
            commutators += 1
    return {'n': n, 'length': length, 'embedding_dimension': embedding_dimension,
            'multiplicities': multiplicities, 'lambdas': lambdas,
            'ground': full_ground, 'residue': residue, 'column': exponential_column,
            'generator_coefficients': relation_coefficients, 'commutators': commutators}


def schedule_greedy(ground, weight, multiplicities, n):
    """Classical latest-free-slot scheduling, independent of Gaussian elimination."""
    free, selected = list(range(1, n + 1)), set()
    for a in sorted(ground, key=lambda a: (dot(a, weight), a)):
        deadline = n - dot(a, multiplicities)
        available = [slot for slot in free if slot <= deadline]
        if available:
            free.remove(max(available))
            selected.add(a)
    return frozenset(selected)


def inverse(matrix):
    n = len(matrix)
    rows = [[Fraction(x) for x in row] + [Fraction(i == j) for j in range(n)]
            for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scale = rows[j][j]
        rows[j] = [x / scale for x in rows[j]]
        for i in range(n):
            if i != j:
                scale = rows[i][j]
                rows[i] = [x - scale * y for x, y in zip(rows[i], rows[j])]
    return [row[n:] for row in rows]


def border_basis(basis, residue):
    basis = sorted(basis)
    n = len(basis)
    dimension = len(residue(basis[0]))
    matrix = [[residue(a)[i] for a in basis] for i in range(dimension)]
    selected_rows, selected = [], []
    for i, row in enumerate(matrix):
        if rank_columns(selected + [row]) > len(selected):
            selected.append(row)
            selected_rows.append(i)
    assert len(selected) == n
    rows = inverse(selected)
    polynomials, normals = [], set()
    for u in minimal_excluded_border(basis):
        column = residue(u)
        coefficients = [sum(x * column[i] for x, i in zip(row, selected_rows)) for row in rows]
        assert tuple(sum(c * residue(a)[i] for c, a in zip(coefficients, basis))
                     for i in range(dimension)) == column
        polynomial = clean({u: Fraction(1), **{a: -c for a, c in zip(basis, coefficients) if c}}, 0)
        polynomials.append(polynomial)
        for a, c in zip(basis, coefficients):
            if c:
                normal = tuple(x - y for x, y in zip(u, a))
                factor = gcd(*normal)
                normals.add(tuple(x // factor for x in normal))
    return polynomials, normals


def fan_control(truncation, multiplicities, exhaustive=False):
    data = model(truncation, multiplicities)
    n = data['length']
    residue = data['residue']
    ground = tuple(a for a in hyperbolic_cross(n, 3) if dot(a, multiplicities) < truncation)
    columns = {a: residue(a) for a in ground}
    comparisons = sorted({primitive(tuple(x - y for x, y in zip(a, b)))
                          for a, b in combinations(ground, 2)
                          if min(x - y for x, y in zip(a, b)) < 0 < max(x - y for x, y in zip(a, b))})
    pending, visited = [(1, n + 1, (n + 1)**2)], {}
    spairs, matrix_decisions, exchange_controls = 0, 0, 0
    while pending:
        w = pending.pop()
        basis = schedule_greedy(ground, w, multiplicities, truncation)
        if basis in visited:
            continue
        assert len(basis) == n
        assert greedy(ground, columns, w, 0) == set(basis)
        assert schedule_greedy(data['ground'], w, multiplicities, truncation) == basis
        assert all(tuple(x - int(i == j) for i, x in enumerate(a)) in basis
                   for a in basis for j in range(3) if a[j])
        matrix_decisions += 1
        polynomials, normals = border_basis(basis, residue)
        assert all(dot(w, normal) > 0 for normal in normals)
        assert set(leading(g, w) for g in polynomials) == set(minimal_excluded_border(sorted(basis)))
        polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
        for normal in sorted(normals):
            polygon = clip(polygon, normal, 1)
        assert polygon and area2(polygon) > 0
        # A second cone construction uses the exact nested-matroid exchange rule.
        exchange_polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
        for u in sorted(set(ground) - basis):
            for v in sorted(basis):
                if independent([dot(a, multiplicities) for a in (basis - {v}) | {u}], truncation):
                    normal = tuple(x - y for x, y in zip(u, v))
                    assert dot(w, normal) > 0
                    exchange_polygon = clip(exchange_polygon, normal, 1)
                    exchange_controls += 1
        assert set(exchange_polygon) == set(polygon)
        if n <= 6 or multiplicities[0] > 1:
            for i, j in combinations(range(len(polynomials)), 2):
                g, h = polynomials[i], polynomials[j]
                a, b = leading(g, w), leading(h, w)
                common = tuple(max(x, y) for x, y in zip(a, b))
                first = {add(u, tuple(x - y for x, y in zip(common, a))): c for u, c in g.items()}
                pair = subtract_shift(first, h, tuple(x - y for x, y in zip(common, b)), 1, 0)
                assert not normal_form(pair, polynomials, w, 0)
                spairs += 1
        visited[basis] = (polygon, w)
        pending.extend(adjacent_weights(polygon, sorted(normals), comparisons))
    assert sum(area2(p) for p, _ in visited.values()) == 1
    chamber_count = 0
    if exhaustive:
        weights = positive_chambers(comparisons)
        assert {schedule_greedy(ground, w, multiplicities, truncation) for w in weights} == set(visited)
        chamber_count = len(weights)
    return {'ambient_truncation': truncation, 'length': n,
            'embedding_dimension': data['embedding_dimension'],
            'valuations': multiplicities, 'exponential_parameters': data['lambdas'],
            'explicit_y_z_coefficients': data['generator_coefficients'],
            'full_valuation_ground_size': len(data['ground']), 'active_universal_ground_size': len(ground),
            'initial_ideal_count': len(visited), 'complete_simplex_double_area': '1',
            'commuting_multiplication_checks': data['commutators'],
            'independent_matrix_selections': matrix_decisions, 'ordinary_s_pairs': spairs,
            'independent_exchange_inequalities': exchange_controls,
            'full_comparison_chambers': chamber_count,
            'witnesses': [{'standard_set': sorted(b), 'strict_weight': w,
                           'closed_cell_vertices': [[str(c) for c in p] for p in poly]}
                          for b, (poly, w) in sorted(visited.items(), key=lambda item: sorted(item[0]))]}


def rank_controls():
    fixtures = []
    for n, valuations in ((4, (1, 2, 3)), (5, (1, 2, 3)), (5, (1, 3, 4)),
                          (10, (3, 4, 5)), (12, (4, 5, 6))):
        data = model(n, valuations)
        ground = data['ground']
        columns = {a: data['column'](a) for a in ground}
        ranks = determinants = 0
        for mask in range(1 << len(ground)):
            selected = [a for i, a in enumerate(ground) if mask >> i & 1]
            values = [dot(a, valuations) for a in selected]
            actual = rank_columns([columns[a] for a in selected])
            assert actual == valuation_rank(values, n)
            assert (actual == len(selected)) == independent(values, n)
            ranks += 1
            if n == 4 and len(selected) == n:
                determinant = sum((-1)**sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n)) *
                                  prod(columns[a][p[j]] for j, a in enumerate(selected))
                                  for p in permutations(range(n)))
                assert bool(determinant) == (actual == n)
                determinants += 1
        # Loops outside the valuation simplex also obey the rank formula.
        assert rank_columns([data['column']((n, 0, 0))]) == valuation_rank([n], n) == 0
        lines = {tuple(c / next(x for x in column if x) for c in column)
                 for column in columns.values() if any(column)}
        nonbinomial = len(lines) > data['length']
        assert nonbinomial
        if valuations[0] == 1 and n >= valuations[1] + 2:
            extra = (0, 1, 0)
            assert all(rank_columns([columns[extra], data['column']((j, 0, 0))]) == 2
                       for j in range(n))
        fixtures.append({'ambient_truncation': n, 'length': data['length'],
                         'embedding_dimension': data['embedding_dimension'],
                         'valuations': valuations, 'all_ground_subsets': ranks,
                         'permutation_determinants': determinants,
                         'nonbinomial_residue_class_control': nonbinomial,
                         'distinct_nonzero_residue_lines': len(lines),
                         'commuting_multiplication_checks': data['commutators']})
    return fixtures


def downsets(n, dimension):
    states = {frozenset({(0,) * dimension})}
    for _ in range(1, n):
        states = {state | {a} for state in states
                  for a in minimal_excluded_border(sorted(state))}
    return sorted(states, key=lambda state: sorted(state))


def consequence_controls():
    generic = []
    for n in range(1, 8):
        data = model(n, (1, 1, 1))
        states = downsets(n, 3)
        for state in states:
            assert independent([sum(a) for a in state], n)
            assert rank_columns([data['column'](a) for a in state]) == n
        generic.append({'length': n, 'all_spatial_downsets': len(states),
                        'every_downset_a_matrix_basis': True})
    products = []
    for q in range(2, 7):
        data = model(2 * q, (1, 1, q))
        fixtures = []
        for state in downsets(q, 2):
            low, high = Fraction(0), None
            feasible = True
            for u in minimal_excluded_border(sorted(state)):
                for v in state:
                    dx, dy = u[0] - v[0], u[1] - v[1]
                    if dx > 0:
                        low = max(low, Fraction(-dy, dx))
                    elif dx < 0:
                        bound = Fraction(-dy, dx)
                        high = bound if high is None else min(high, bound)
                    elif dy <= 0:
                        feasible = False
            if not feasible or (high is not None and low >= high):
                continue
            ratio = low + 1 if high is None else (low + high) / 2
            planar = (ratio, Fraction(1))
            gap = min(dot(u, planar) for u in minimal_excluded_border(sorted(state))) - \
                max(dot(v, planar) for v in state)
            assert gap > 0
            w = integer_weight(planar + (gap / 2,))
            basis = frozenset(a + (c,) for a in state for c in (0, 1))
            assert schedule_greedy(data['ground'], w, (1, 1, q), 2 * q) == basis
            assert rank_columns([data['column'](a) for a in basis]) == 2 * q
            u = max((a for a in state if bool(a[0]) != bool(a[1])), key=sum)
            assert tuple(2 * x for x in u) not in state
            midpoint = u + (1,)
            outside = [tuple(2 * x for x in u) + (0,), (0, 0, 2)]
            assert midpoint in basis and all(a not in basis for a in outside)
            assert tuple(2 * x for x in midpoint) == add(*outside)
            fixtures.append({'planar_standard_set': sorted(state), 'strict_weight': w,
                             'included_midpoint': midpoint, 'excluded_endpoints': outside})
        products.append({'length': 2 * q, 'valuations': [1, 1, q],
                         'planar_corner_cuts_lifted': len(fixtures), 'fixtures': fixtures})
    return {'generic_curvilinear': generic, 'non_corner_cut_products': products}


def main():
    ranks = rank_controls()
    consequences = consequence_controls()
    fans = [fan_control(n, v, exhaustive=n <= 5)
            for n, v in ((4, (1, 2, 3)), (5, (1, 2, 3)), (6, (1, 2, 3)),
                         (8, (1, 2, 3)), (8, (1, 3, 5)), (7, (1, 1, 1)),
                         (10, (3, 4, 5)), (12, (4, 5, 6)))]
    paths = ['scripts/run_nested_jet_checks.py', 'docs/nested_jet_matroids.md',
             'scripts/run_extremal_initial_ideals_checks.py', 'scripts/run_initial_cone_height_checks.py',
             'scripts/run_interleaving_research.py', 'scripts/run_structural_complexity_checks.py']
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_nested_jet_checks.py',
              'scope': 'Explicit characteristic-zero polynomial-quotient realizations of valuation-nested matroids, including embedding dimension three. Classical interpolation and scheduling inputs are not novelty claims. No spatial growth law is inferred from finite fan counts.',
              'summary': {'exhaustive_subset_ranks': sum(r['all_ground_subsets'] for r in ranks),
                          'permutation_determinants': sum(r['permutation_determinants'] for r in ranks),
                          'complete_positive_fans': len(fans),
                          'initial_ideals_checked': sum(r['initial_ideal_count'] for r in fans),
                          'ordinary_s_pairs': sum(r['ordinary_s_pairs'] for r in fans),
                          'exchange_inequalities': sum(r['independent_exchange_inequalities'] for r in fans),
                          'commuting_multiplication_checks': sum(r['commuting_multiplication_checks'] for r in ranks + fans),
                          'all_downset_matrix_bases': sum(r['all_spatial_downsets'] for r in consequences['generic_curvilinear']),
                          'non_corner_cut_product_bases': sum(r['planar_corner_cuts_lifted'] for r in consequences['non_corner_cut_products'])},
              'rank_controls': ranks, 'fan_controls': fans,
              'uniform_consequence_controls': consequences,
              'source_sha256': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
              'goal_complete': False, 'asymptotic_count_proved': False, 'novelty_established': False}
    destination = ROOT / 'results/research/nested_jet_checks.json'
    destination.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': 'passed', **report['summary'],
                      'fan_counts': [[f['length'], f['valuations'], f['initial_ideal_count']] for f in fans]}, indent=2))


if __name__ == '__main__':
    main()
