#!/usr/bin/env python3
"""Exact controls for mixed-image invariance and high-order rank splitting."""
from fractions import Fraction
from itertools import combinations, permutations, product
from math import prod
from pathlib import Path
import hashlib
import json

from run_nested_jet_checks import rank_columns
from run_structural_complexity_checks import dot, greedy, positive_chambers, primitive

ROOT = Path(__file__).resolve().parents[1]


def rank(columns, prime):
    if not prime:
        return rank_columns(columns)
    pivots = {}
    for column in columns:
        vector = [x % prime for x in column]
        for i, row in sorted(pivots.items()):
            scale = vector[i]
            vector = [(x-scale*y) % prime for x, y in zip(vector, row)]
        i = next((j for j, x in enumerate(vector) if x), None)
        if i is not None:
            scale = pow(vector[i], -1, prime)
            pivots[i] = [x*scale % prime for x in vector]
    return len(pivots)


def algebra(cells, coefficients, prime=0):
    cells = tuple(sorted(cells))
    index = {u: i for i, u in enumerate(cells)}
    n, S = len(cells), max(map(sum, cells))
    assert all(tuple(x-int(i == j) for i, x in enumerate(u)) in index
               for u in cells for j in range(2) if u[j])
    def reduce(poly):
        return {u: (c % prime if prime else Fraction(c)) for u, c in poly.items()
                if u in index and (c % prime if prime else c)}
    def multiply(left, right):
        answer = {}
        for (a, b), c in left.items():
            for (i, j), d in right.items():
                u = (a+i, b+j)
                if u in index:
                    answer[u] = answer.get(u, 0)+c*d
        return reduce(answer)
    def power(poly, k):
        answer = {(0, 0): 1}
        for _ in range(k):
            answer = multiply(answer, poly)
        return answer
    F = reduce(coefficients)
    assert (0, 0) not in F
    powers = [power(F, k) for k in range(S+2)]
    assert not powers[-1]
    def vector(poly):
        return tuple(poly.get(u, 0) for u in cells)
    def residue(exponent):
        a, b, c = exponent
        if c > S:
            return (0,)*n
        return vector(reduce({(i+a, j+b): value for (i, j), value in powers[c].items()}))
    def operator(poly):
        return tuple(vector(multiply(poly, {u: 1})) for u in cells)
    def apply(matrix, column):
        result = tuple(sum(v[i]*c for v, c in zip(matrix, column)) for i in range(n))
        return tuple(x % prime for x in result) if prime else result
    ground = tuple((a, b, c) for a, b in cells for c in range(S-a-b+1))
    columns = {u: residue(u) for u in ground}
    # Independent replay using multiplication matrices.
    operators = [operator({(1, 0): 1}), operator({(0, 1): 1}), operator(F)]
    one = vector({(0, 0): 1})
    for u in ground:
        replay = one
        for matrix, k in zip(operators, u):
            for _ in range(k):
                replay = apply(matrix, replay)
        assert replay == columns[u]
    return {'cells': cells, 'index': index, 'length': n, 'socle_degree': S,
            'F': F, 'prime': prime, 'ground': ground, 'columns': columns,
            'multiply': multiply, 'power': power, 'vector': vector,
            'operator': operator, 'apply': apply, 'residue': residue}


def rectangle(r, s, coefficients, prime=0):
    return algebra(product(range(r), range(s)), coefficients, prime)


def automorphism_controls():
    fixtures = []
    for r, s, prime in ((2, 2, 0), (3, 3, 0), (3, 4, 0), (3, 3, 2), (4, 3, 3), (4, 4, 101)):
        F = {(i, j): (i+2)*(j+3)+1 for i in range(r) for j in range(s) if i or j}
        F[(1, 0)] = F[(0, 1)] = 1
        data = rectangle(r, s, F, prime)
        linear = rectangle(r, s, {(1, 0): 1, (0, 1): 1}, prime)
        cells = data['cells']
        pure_y = {u: c for u, c in data['F'].items() if u[0] == 0}
        x_image = {u: c for u, c in data['F'].items() if u[0]}
        phi = [data['vector'](data['multiply'](data['power'](x_image, i),
                                               data['power'](pure_y, j))) for i, j in cells]
        assert rank(phi, prime) == r*s
        for cell, column in zip(cells, phi):
            for row, value in zip(cells, column):
                if sum(row) <= sum(cell):
                    assert value == int(row == cell)
        for column, target in zip(phi, linear['operator']({(1, 0): 1, (0, 1): 1})):
            assert data['apply'](data['operator'](data['F']), column) == data['apply'](phi, target)
        ideal_checks = 0
        for a in range(r):
            for b in range(s):
                chosen = [column for u, column in zip(cells, phi) if u[0] >= a and u[1] >= b]
                assert all(not c or (v[0] >= a and v[1] >= b)
                           for column in chosen for v, c in zip(cells, column))
                assert rank(chosen, prime) == (r-a)*(s-b)
                ideal_checks += 1
        images = []
        for a in range(r):
            for b in range(s):
                for c in range(r+s):
                    left = [data['residue']((a+i, b+j, c)) for i, j in cells]
                    right = [linear['residue']((a+i, b+j, c)) for i, j in cells]
                    actual = rank(left, prime)
                    assert actual == rank(right, prime)
                    formula = None
                    if not prime:
                        R, T = r-a, s-b
                        formula = sum(max(R+T-1-2*j-c, 0) for j in range(min(R, T)))
                        assert actual == formula
                    images.append({'monomial': [a, b, c], 'rank': actual,
                                   'characteristic_zero_formula': formula})
        fixtures.append({'rectangle': [r, s], 'length': r*s, 'characteristic': prime,
                         'F_terms': [[*u, int(c)] for u, c in sorted(data['F'].items())],
                         'automorphism_columns': [[str(c) for c in column] for column in phi],
                         'coordinate_ideal_checks': ideal_checks, 'mixed_image_ranks': images})
    return fixtures


def determinant(columns):
    n = len(columns)
    answer = Fraction(0)
    for order in permutations(range(n)):
        sign = (-1)**sum(order[i] > order[j] for i in range(n) for j in range(i+1, n))
        answer += sign*prod(columns[j][order[j]] for j in range(n))
    return answer


def small_fan_and_image_controls():
    data = [rectangle(2, 2, {(1, 0): 1, (0, 1): 1}),
            rectangle(2, 2, {(1, 0): 1, (0, 1): 1, (1, 1): 1})]
    ground = tuple(u for u in data[0]['ground'] if any(data[0]['columns'][u]) or any(data[1]['columns'][u]))
    normals = sorted({primitive(tuple(a-b for a, b in zip(u, v))) for u, v in combinations(ground, 2)
                      if min(a-b for a, b in zip(u, v)) < 0 < max(a-b for a, b in zip(u, v))})
    weights = positive_chambers(normals)
    fan_records, determinant_count = [], 0
    for label, model in zip(('linear', 'nonlinear'), data):
        bases = []
        for selected in combinations(ground, 4):
            value = determinant([model['columns'][u] for u in selected])
            assert bool(value) == (rank_columns([model['columns'][u] for u in selected]) == 4)
            if value:
                bases.append(frozenset(selected))
            determinant_count += 1
        witnesses = {}
        for weight in weights:
            best = min(sum(dot(weight, u) for u in basis) for basis in bases)
            selected = [basis for basis in bases if sum(dot(weight, u) for u in basis) == best]
            assert len(selected) == 1
            basis = selected[0]
            assert greedy(ground, model['columns'], weight, 0) == set(basis)
            witnesses.setdefault(basis, weight)
        fan_records.append({'label': label, 'initial_ideal_count': len(witnesses),
                            'complete_comparison_intervals': len(weights),
                            'witnesses': [{'basis': sorted(b), 'strict_weight': w}
                                          for b, w in sorted(witnesses.items(), key=lambda pair: sorted(pair[0]))]})
    assert [f['initial_ideal_count'] for f in fan_records] == [5, 6]
    image_controls = []
    for prime in (0, 2, 3):
        linear = rectangle(2, 2, {(1, 0): 1, (0, 1): 1}, prime)
        nonlinear = rectangle(2, 2, {(1, 0): 1, (0, 1): 1, (1, 1): 1}, prime)
        ideals = [{u: [model['residue']((u[0]+i, u[1]+j, u[2])) for i, j in model['cells']]
                   for u in ground} for model in (linear, nonlinear)]
        differences = 0
        for mask in range(1 << len(ground)):
            selected = [u for i, u in enumerate(ground) if mask >> i & 1]
            ranks = [rank([column for u in selected for column in image[u]], prime) for image in ideals]
            assert ranks[0] == ranks[1]
            column_ranks = [rank([model['columns'][u] for u in selected], prime)
                            for model in (linear, nonlinear)]
            differences += column_ranks[0] != column_ranks[1]
        generators = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
        assert [rank([model['columns'][u] for u in generators], prime)
                for model in (linear, nonlinear)] == [3, 4]
        image_controls.append({'characteristic': prime, 'all_operator_list_subsets': 1 << len(ground),
                               'subsets_with_different_column_rank': differences,
                               'generator_column_ranks': [3, 4]})
    return {'ground': ground, 'permutation_determinants': determinant_count,
            'complete_rational_fans': fan_records, 'all_joint_image_subset_controls': image_controls}


def projected_greedy(points, columns, weight, prime):
    selected = []
    for point in sorted(points, key=lambda u: (dot(weight, u), u)):
        if rank([columns[u] for u in selected+[point]], prime) > len(selected):
            selected.append(point)
    return set(selected)


def high_order_controls():
    fixtures = []
    specifications = [(q, h, prime) for q, h, prime in
                      ((3, 2, 0), (4, 2, 0), (5, 3, 0), (6, 3, 0), (4, 2, 2), (5, 3, 3))]
    for q, h, prime in specifications:
        D = 2*q-h
        cells = tuple(product(range(q), repeat=2))
        F = {(1, 0): 1, (0, 1): 1}
        F.update({(i, j): (i+2)*(j+3)+1 for i, j in cells if i+j >= D})
        fixtures.append(decomposition_control(cells, F, D, prime,
                                              {'rectangle': [q, q], 'terminal_band_parameter': h}))
    cells = [(i, j) for i in range(4) for j in range(4) if i < 2 or j < 2]
    fixtures.append(decomposition_control(cells, {(1, 0): 1, (0, 1): 1, (3, 1): 2, (1, 3): 3},
                                          4, 0, {'base_monomial_ideal': ['x^4', 'y^4', 'x^2*y^2']}))
    for prime in (0, 2, 3):
        F = {(1, 0): 1, (0, 1): 1, (3, 0): 1}
        data = nonmonomial_homogeneous_model(F, prime)
        linear = nonmonomial_homogeneous_model({(1, 0): 1, (0, 1): 1}, prime)
        fixtures.append(decomposition_control(data['cells'], F, 3, prime,
                                              {'base_homogeneous_ideal': ['x^3-y^3', 'x^2*y', 'x*y^2']},
                                              data=data, linear=linear))
    return fixtures


def nonmonomial_homogeneous_model(F, prime):
    """A graded quotient with y^3=x^3 and all mixed cubics zero."""
    cells = tuple(sorted([(a, b) for a in range(3) for b in range(3-a)]+[(3, 0)]))
    index = {u: i for i, u in enumerate(cells)}
    def reduce(poly):
        answer = {}
        for u, c in poly.items():
            v = (3, 0) if u == (0, 3) else u
            if v in index:
                answer[v] = answer.get(v, 0)+c
        return {u: (c % prime if prime else Fraction(c)) for u, c in answer.items()
                if (c % prime if prime else c)}
    def raw_multiply(left, right):
        answer = {}
        for (a, b), c in left.items():
            for (i, j), d in right.items():
                answer[(a+i, b+j)] = answer.get((a+i, b+j), 0)+c*d
        return answer
    def vector(poly):
        return tuple(poly.get(u, 0) for u in cells)
    def operator(poly):
        return tuple(vector(reduce(raw_multiply(poly, {u: 1}))) for u in cells)
    def apply(matrix, column):
        values = tuple(sum(v[i]*c for v, c in zip(matrix, column)) for i in range(len(cells)))
        return tuple(x % prime for x in values) if prime else values
    operators = [operator({(1, 0): 1}), operator({(0, 1): 1}), operator(F)]
    for left, right in combinations(operators, 2):
        for j in range(len(cells)):
            assert apply(left, right[j]) == apply(right, left[j])
    # Direct polynomial expansion followed by the explicit homogeneous relation,
    # checked against repeated application of the commuting operators.
    ground = tuple((a, b, c) for a in range(4) for b in range(4-a) for c in range(4-a-b))
    powers = [{(0, 0): 1}]
    for _ in range(3):
        powers.append(raw_multiply(powers[-1], F))
    columns = {u: vector(reduce({(i+u[0], j+u[1]): value for (i, j), value in powers[u[2]].items()}))
               for u in ground}
    for u in ground:
        replay = vector({(0, 0): 1})
        for matrix, k in zip(operators, u):
            for _ in range(k):
                replay = apply(matrix, replay)
        assert replay == columns[u]
    assert columns[(3, 0, 0)] == columns[(0, 3, 0)]
    assert not any(columns[(2, 1, 0)]) and not any(columns[(1, 2, 0)])
    return {'cells': cells, 'length': len(cells), 'socle_degree': 3,
            'F': reduce(F), 'prime': prime, 'ground': ground, 'columns': columns}


def decomposition_control(cells, F, D, prime, description, data=None, linear=None):
    if data is None:
        data = algebra(cells, F, prime)
        linear = algebra(cells, {(1, 0): 1, (0, 1): 1}, prime)
    S, n = data['socle_degree'], data['length']
    t = S-D+1
    assert 2*D > S+1
    components = {'combined': [i for i, u in enumerate(data['cells']) if sum(u) <= t or sum(u) >= D]}
    components.update({str(j): [i for i, u in enumerate(data['cells']) if sum(u) == j]
                       for j in range(t+1, D)})
    R = len(components['combined'])
    if 'terminal_band_parameter' in description:
        assert R == description['terminal_band_parameter']**2
    assigned = {u: str(sum(u)) if t < sum(u) < D else 'combined' for u in data['ground']}
    projected = {}
    for u in data['ground']:
        target = assigned[u]
        rows = components[target]
        column = data['columns'][u]
        assert all(not c or i in rows for i, c in enumerate(column))
        if target != 'combined':
            assert column == linear['columns'][u]
        projected[u] = tuple(column[i] for i in rows)
    for label, rows in components.items():
        points = [u for u in data['ground'] if assigned[u] == label]
        assert rank([projected[u] for u in points], prime) == len(rows)
    scale = S+1
    weights = sorted(set(permutations((1, scale, scale**2))))
    weights += [(1, 1, 1), (3, 4, 5), (4*scale**2, 2*scale-1, 1)]
    witnesses = []
    for weight in weights:
        full = greedy(data['ground'], data['columns'], weight, prime)
        split = set()
        groups = {}
        for label in components:
            points = [u for u in data['ground'] if assigned[u] == label]
            chosen = projected_greedy(points, projected, weight, prime)
            split.update(chosen)
            assert len(chosen) == len(components[label])
            groups[label] = sorted(chosen)
        assert split == full and len(full) == n
        assert rank([data['columns'][u] for u in full], prime) == n
        witnesses.append({'weight': weight, 'basis': sorted(full), 'component_bases': groups})
    # Every actual row coordinate count is used; no guessed rank parameter.
    assert sum(map(len, components.values())) == n
    return {**description, 'characteristic': prime, 'length': n, 'socle_degree': S,
            'first_allowed_nonlinear_degree': D, 'low_degree_cutoff': t, 'combined_rank': R,
            'unchanged_middle_rank': n-R, 'all_ground_columns_checked': len(data['ground']),
            'component_row_indices': components,
            'F_terms': [[*u, int(c)] for u, c in sorted(data['F'].items())],
            'fan_bound_scale_squared': (1+R*S)**3,
            'selection_controls': witnesses}


def main():
    automorphisms = automorphism_controls()
    small = small_fan_and_image_controls()
    splits = high_order_controls()
    summary = {'automorphism_fixtures': len(automorphisms),
               'coordinate_ideal_checks': sum(f['coordinate_ideal_checks'] for f in automorphisms),
               'mixed_multiplication_rank_comparisons': sum(len(f['mixed_image_ranks']) for f in automorphisms),
               'complete_small_rational_fans': len(small['complete_rational_fans']),
               'independent_permutation_determinants': small['permutation_determinants'],
               'all_joint_image_subset_checks': sum(f['all_operator_list_subsets'] for f in small['all_joint_image_subset_controls']),
               'high_order_decomposition_fixtures': len(splits),
               'all_decomposition_columns': sum(f['all_ground_columns_checked'] for f in splits),
               'full_vs_projected_selection_controls': sum(len(f['selection_controls']) for f in splits)}
    sources = ['scripts/run_rectangle_rank_checks.py', 'docs/rectangle_rank_structure.md',
               'scripts/run_nested_jet_checks.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_extremal_initial_ideals_checks.py', 'scripts/run_initial_cone_height_checks.py',
               'scripts/run_interleaving_research.py']
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_rectangle_rank_checks.py',
              'scope': 'All-characteristic mixed-image invariance and exact rank decomposition for high-order coordinate perturbations; a classical lattice-vertex corollary bounds the restricted fans. Jordan data do not determine selected-column ranks. No unrestricted G_3 improvement, power-separated family or novelty certificate is claimed.',
              'summary': summary, 'automorphism_controls': automorphisms,
              'identical_image_data_different_fans': small, 'high_order_rank_decompositions': splits,
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              'entire_rectangle_family_closed': False, 'unrestricted_endpoint_achieved': False,
              'novelty_established': False, 'goal_complete': False}
    (ROOT/'results/research/rectangle_rank_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', **summary}, indent=2))


if __name__ == '__main__':
    main()
