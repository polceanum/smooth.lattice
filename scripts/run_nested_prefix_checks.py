#!/usr/bin/env python3
"""Exact controls for the uniform additive nested-matroid prefix bound."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import factorial, gcd, prod
from pathlib import Path
import hashlib
import json

from run_nested_jet_checks import independent, model, product_truncated, rank_columns, valuation_rank
from run_structural_complexity_checks import dot, primitive

ROOT = Path(__file__).resolve().parents[1]


def certificate(ground, basis, weight, m, removed, added, is_basis):
    dimension, rank = len(m), len(basis)
    assert removed in basis and added not in basis
    assert is_basis((basis-{removed}) | {added})
    threshold = max(dot(m, removed), dot(m, added))
    level = dot(weight, removed)
    assert level == dot(weight, added) > 0
    prefix = [u for u in ground if dot(m, u) <= threshold]
    # Test the full replacement prefix, including candidates earlier than
    # either endpoint's allocation block.
    for u in prefix:
        if u not in basis:
            assert is_basis((basis-{removed}) | {u})
    closed = {u for u in prefix if dot(weight, u) <= level}
    strict = {u for u in closed if dot(weight, u) < level}
    assert strict <= basis and len(strict) <= rank-1
    assert removed in closed and added in closed
    halves = Counter(tuple(x//2 for x in u) for u in closed)
    assert set(halves) <= strict and max(halves.values()) <= 2**dimension
    assert len(closed) <= 2**dimension*(rank-1)
    box = tuple(max(u[i] for u in closed)+1 for i in range(dimension))
    inner = set(product(*(range((side-1)//dimension+1) for side in box)))
    assert inner <= closed and dimension**dimension*len(inner) >= prod(box)
    assert prod(box) <= dimension**dimension*len(closed) <= (2*dimension)**dimension*(rank-1)
    difference = tuple(a-b for a, b in zip(added, removed))
    raw_product = prod(abs(x)+1 for x in difference)
    direction = primitive(difference)
    assert raw_product <= prod(box)
    assert prod(abs(x)+1 for x in direction) <= raw_product
    return {'valuation_threshold': threshold, 'cost_threshold': str(level),
            'closed_cut': sorted(closed), 'strict_cut_size': len(strict),
            'closed_cut_box': box, 'box_volume': prod(box),
            'removed': removed, 'added': added, 'primitive_direction': direction,
            'raw_difference_product': raw_product}


def dimensional_model(L, m):
    """The same explicit realization, independently extended to d=2,4."""
    d, B = len(m), max(L, 2)
    C = L*sum(B**i for i in range(d))
    lambdas = tuple(C*mi+B**i for i, mi in enumerate(m))
    ground = tuple(u for u in product(range(L), repeat=d) if dot(m, u) < L)
    def column(u):
        v, node = dot(m, u), dot(lambdas, u)
        return tuple(Fraction(node**(i-v), factorial(i-v)) if i >= v else Fraction(0)
                     for i in range(L))
    columns = {u: column(u) for u in ground}
    generators = [column(tuple(int(i == j) for i in range(d))) for j in range(d)]
    one = (Fraction(1),)+(Fraction(0),)*(L-1)
    for u in ground:
        replay = one
        for generator, power in zip(generators, u):
            for _ in range(power):
                replay = product_truncated(replay, generator)
        assert replay == columns[u]
    length = rank_columns(list(columns.values()))
    formula = L-max(k-sum(dot(m, u) < k for u in ground) for k in range(L+1))
    assert length == formula
    return {'ground': ground, 'columns': columns, 'length': length, 'lambdas': lambdas}


def table_control(ground, m, bounds, columns=None):
    rank = len(bounds)
    bases, checked = set(), 0
    for selected in combinations(ground, rank):
        criterion = all(v <= bound for v, bound in zip(sorted(dot(m, u) for u in selected), bounds))
        if columns is not None:
            assert (rank_columns([columns[u] for u in selected]) == rank) == criterion
        if criterion:
            bases.add(frozenset(selected))
        checked += 1
    assert bases
    exchange_axioms = lowering = prefix_replacements = 0
    for basis in bases:
        for other in bases:
            for removed in basis-other:
                assert any((basis-{removed}) | {added} in bases for added in other-basis)
                exchange_axioms += 1
        for removed in basis:
            for added in set(ground)-basis:
                replacement = (basis-{removed}) | {added}
                if dot(m, added) <= dot(m, removed):
                    assert replacement in bases
                    lowering += 1
                if replacement in bases:
                    threshold = max(dot(m, removed), dot(m, added))
                    for candidate in set(ground)-basis:
                        if dot(m, candidate) <= threshold:
                            assert (basis-{removed}) | {candidate} in bases
                            prefix_replacements += 1
    return bases, {'ground_size': len(ground), 'rank': rank, 'sorted_valuation_bounds': bounds,
                   'candidate_bases': checked, 'matrix_basis_tests': checked if columns is not None else 0,
                   'bases': len(bases), 'basis_exchange_axioms': exchange_axioms,
                   'lowering_replacements': lowering, 'full_prefix_replacements': prefix_replacements}


def hypothesis_controls():
    fixtures = []
    for L, m in ((4, (1, 2, 3)), (5, (1, 2, 3)), (10, (3, 4, 5)), (12, (4, 5, 6)),
                 (5, (1, 2)), (3, (1, 1, 2, 2))):
        data = model(L, m) if len(m) == 3 else dimensional_model(L, m)
        ground, rank = data['ground'], data['length']
        columns = {u: data['column'](u) for u in ground} if len(m) == 3 else data['columns']
        bounds = tuple(L-rank+i for i in range(rank))
        _, control = table_control(ground, m, bounds, columns)
        fixtures.append({'ambient_truncation': L, 'valuations': m, 'dimension': len(m),
                         'actual_quotient_length': rank, **control})
    return fixtures


def complete_fan_controls():
    old = json.loads((ROOT/'results/research/nested_jet_checks.json').read_text())
    fixtures = []
    for fan in old['fan_controls']:
        L, m, rank = fan['ambient_truncation'], tuple(fan['valuations']), fan['length']
        data = model(L, m)
        ground = data['ground']
        value = {u: dot(m, u) for u in ground}
        def is_basis(selected):
            return independent([value[u] for u in selected], L)
        facets, active, matrices, witnesses, directions = 0, 0, 0, [], set()
        for witness in fan['witnesses']:
            basis = frozenset(map(tuple, witness['standard_set']))
            assert rank_columns([data['column'](u) for u in basis]) == rank
            polygon = [tuple(map(Fraction, v)) for v in witness['closed_cell_vertices']]
            for left, right in zip(polygon, polygon[1:]+polygon[:1]):
                weight = tuple((a+b)/2 for a, b in zip(left, right))
                if min(weight) <= 0:
                    continue
                records = []
                for removed in sorted(basis):
                    for added in ground:
                        if added in basis or dot(weight, removed) != dot(weight, added):
                            continue
                        replacement = (basis-{removed}) | {added}
                        if not is_basis(replacement):
                            continue
                        record = certificate(ground, basis, weight, m, removed, added, is_basis)
                        assert rank_columns([data['column'](u) for u in replacement]) == rank
                        matrices += 1
                        active += 1
                        directions.add(tuple(record['primitive_direction']))
                        records.append(record)
                assert records
                facets += 1
                witnesses.append({'basis': sorted(basis), 'positive_facet_weight': list(map(str, weight)),
                                  'active_exchange_certificates': records})
        fixtures.append({'ambient_truncation': L, 'valuations': m, 'actual_quotient_length': rank,
                         'embedding_dimension': data['embedding_dimension'],
                         'complete_initial_ideals': len(fan['witnesses']),
                         'positive_cell_facets': facets, 'active_exchanges': active,
                         'independent_exchange_matrix_tests': matrices,
                         'primitive_directions': sorted(directions), 'witnesses': witnesses})
    return fixtures


def common_planar_weights(rank):
    H = 16*rank
    rays = {(q, p) for p in range(1, H) for q in range(1, H//(p+1))
            if gcd(p, q) == 1 and (p+1)*(q+1) <= H}
    ordered = [(0, 1)] + sorted(rays, key=lambda r: Fraction(r[0], r[1])) + [(1, 0)]
    weights = [tuple(a+b for a, b in zip(left, right)) for left, right in zip(ordered, ordered[1:])]
    assert all(min(w) > 0 and max(w) <= 2*H for w in weights)
    return weights


def abstract_controls():
    fixtures = []
    for L, m, bounds in ((4, (1, 1), (0, 1, 3)), (7, (2, 3), (0, 4, 6)),
                          (3, (1, 1, 2, 2), (1, 2, 2))):
        dimension, rank = len(m), len(bounds)
        ground = tuple(u for u in product(range(L), repeat=dimension) if dot(m, u) < L)
        bases, control = table_control(ground, m, bounds)
        def optima(weight):
            totals = {basis: sum(dot(weight, u) for u in basis) for basis in bases}
            best = min(totals.values())
            return {basis for basis, value in totals.items() if value == best}
        witnesses = []
        for weight in (tuple([1]*dimension), tuple(range(1, dimension+1)),
                       tuple(range(dimension, 0, -1))):
            for basis in sorted(optima(weight), key=lambda b: sorted(b)):
                for removed in sorted(basis):
                    for added in ground:
                        if added in basis or dot(weight, removed) != dot(weight, added):
                            continue
                        if (basis-{removed}) | {added} in bases:
                            record = certificate(ground, basis, weight, m, removed, added,
                                                 lambda selected: selected in bases)
                            witnesses.append({'basis': sorted(basis), 'weight': weight, **record})
        common = None
        if dimension == 2:
            walls = sorted({Fraction(v[1]-u[1], u[0]-v[0]) for u, v in combinations(ground, 2)
                            if u[0] != v[0] and Fraction(v[1]-u[1], u[0]-v[0]) > 0})
            probes = [walls[0]/2]+[(a+b)/2 for a, b in zip(walls, walls[1:])]+[walls[-1]+1] \
                if walls else [Fraction(1)]
            exact_states, sampled_states = set(), set()
            for ratio in probes:
                states = optima((ratio, Fraction(1)))
                assert len(states) == 1
                exact_states.update(states)
            weights = common_planar_weights(rank)
            for weight in weights:
                states = optima(weight)
                assert len(states) == 1
                sampled_states.update(states)
            assert exact_states == sampled_states
            common = {'rank': rank, 'universal_arrangement_weights': weights,
                      'complete_ground_comparison_intervals': len(probes),
                      'distinct_positive_bases': len(exact_states),
                      'state_list': [sorted(b) for b in sorted(exact_states, key=lambda b: sorted(b))]}
        fixtures.append({'dimension': dimension, 'ground_valuation_cap': L, 'valuations': m,
                         'polynomial_realization_asserted': False, **control,
                         'active_exchange_certificates': witnesses, 'common_planar_list_control': common})
    return fixtures


def two_cone_control():
    M = 512
    p, q, n = M, M+1, 2*M+1
    assert gcd(p, q) == 1
    # Reference basis: 1,x,...,x^p,y,...,y^(q-1); x^p=y^q.
    basis = [(0, 0, 0)]+[(i, 0, 0) for i in range(1, p+1)]+[(0, j, 0) for j in range(1, q)]
    index = {u: i for i, u in enumerate(basis)}
    assert len(index) == n
    def residue(u):
        a, b, c = u
        if c or (a and b) or a > p or b > q:
            return None
        if b == q:
            return index[(p, 0, 0)]
        return index[u]
    operators = [[residue(tuple(x+int(i == j) for i, x in enumerate(u))) for u in basis]
                 for j in range(3)]
    checks = 0
    for X, Y in combinations(operators, 2):
        for i in range(n):
            xy = None if Y[i] is None else X[Y[i]]
            yx = None if X[i] is None else Y[X[i]]
            assert xy == yx
            checks += 1
    for a in range(p+2):
        for b in range(q+2):
            if a and b:
                assert residue((a, b, 0)) is None
    assert residue((p, 0, 0)) == residue((0, q, 0))
    wall_weight = (q, p, 1)
    normal = primitive((p, -q, 0))
    assert dot(normal, wall_weight) == 0
    direction_product = prod(abs(x)+1 for x in normal)
    assert direction_product > 216*(n-1)
    strict = ((q+1, p, 1), (q, p+1, 1))
    assert dot(strict[0], (p, 0, 0)) > dot(strict[0], (0, q, 0))
    assert dot(strict[1], (p, 0, 0)) < dot(strict[1], (0, q, 0))
    return {'p': p, 'q': q, 'length': n,
            'ideal': ['z', 'xy', f'x^{p+1}', f'y^{q+1}', f'x^{p}-y^{q}'],
            'initial_ideal_count': 2, 'primitive_wall_normal': normal,
            'normal_product': direction_product, 'nested_theorem_product_bound': 216*(n-1),
            'wall_weight': wall_weight, 'strict_weights': strict,
            'sparse_commuting_operator_checks': checks,
            'claim_limit': 'This binomial control violates the lowering hypothesis. Its quadratic primitive direction product does not imply a large initial fan.'}


def generic_jet_controls():
    L, m = 5, (1, 2, 3)
    original = model(L, m)
    ground = original['ground']
    good = ((0, 1, 2, 3, 5), (0, 0, 1, 7, 11), (0, 0, 0, 1, 13))
    special = ((0, 1, 0, 0, 0), (0, 0, 1, 0, 0), (0, 0, 0, 1, 0))
    fixtures = []
    for label, generators in (('nonexponential_integer_point', good), ('monomial_specialization', special)):
        columns = {}
        for u in ground:
            value = (Fraction(1),)+(Fraction(0),)*(L-1)
            for generator, power in zip(generators, u):
                for _ in range(power):
                    value = product_truncated(value, generator)
            columns[u] = value
        checked = differences = 0
        for mask in range(1 << len(ground)):
            selected = [u for i, u in enumerate(ground) if mask >> i & 1]
            actual = rank_columns([columns[u] for u in selected])
            expected = valuation_rank([dot(m, u) for u in selected], L)
            assert actual <= expected
            differences += actual != expected
            checked += 1
        x_squared, y = columns[(2, 0, 0)], columns[(0, 1, 0)]
        minor = x_squared[2]*y[3]-y[2]*x_squared[3]
        assert minor == generators[1][3]-2*generators[0][2]
        pair_rank = rank_columns([x_squared, y])
        if label == 'nonexponential_integer_point':
            assert differences == 0 and pair_rank == 2 and minor == 3
        else:
            assert differences > 0 and pair_rank == 1 and minor == 0
        fixtures.append({'label': label, 'ambient_truncation': L, 'valuations': m,
                         'generator_coefficients': generators, 'all_subset_ranks_checked': checked,
                         'rank_pattern_differences': differences,
                         'actual_image_length': rank_columns(list(columns.values())),
                         'minor_Y1_squared_Y2_rows_2_3': str(minor), 'pair_rank': pair_rank})
    # The image length itself can drop on the exceptional coefficient locus.
    L, m = 16, (4, 6, 13)
    data = model(L, m)
    monomial_columns = [tuple(Fraction(i == dot(m, u)) for i in range(L)) for u in data['ground']]
    special_length = rank_columns(monomial_columns)
    assert data['length'] == 9 and special_length == 8
    return {'full_subset_controls': fixtures,
            'image_length_drop_control': {'ambient_truncation': L, 'valuations': m,
                                         'generic_image_length': data['length'],
                                         'monomial_specialization_length': special_length,
                                         'full_nonzero_ground_size': len(data['ground'])},
            'scope': 'One explicit normalized nonexponential Taylor point has every predicted subset rank; a special monomial point has extra dependencies. The general nonempty-open assertion uses the finite nonzero-minor product proof, not these two points.'}


def main():
    hypotheses = hypothesis_controls()
    fans = complete_fan_controls()
    abstract = abstract_controls()
    boundary = two_cone_control()
    jets = generic_jet_controls()
    summary = {'matrix_basis_tests': sum(f['matrix_basis_tests'] for f in hypotheses),
               'matrix_realization_fixtures': len(hypotheses),
               'lowering_replacements': sum(f['lowering_replacements'] for f in hypotheses+abstract),
               'full_prefix_replacements': sum(f['full_prefix_replacements'] for f in hypotheses+abstract),
               'basis_exchange_axioms': sum(f['basis_exchange_axioms'] for f in hypotheses+abstract),
               'saved_complete_fans': len(fans), 'saved_initial_ideals': sum(f['complete_initial_ideals'] for f in fans),
               'positive_cell_facets': sum(f['positive_cell_facets'] for f in fans),
               'active_exchange_certificates': sum(f['active_exchanges'] for f in fans),
               'independent_exchange_matrix_tests': sum(f['independent_exchange_matrix_tests'] for f in fans),
               'abstract_flag_fixtures': len(abstract),
               'abstract_active_exchange_certificates': sum(len(f['active_exchange_certificates']) for f in abstract),
               'complete_common_weight_list_controls': sum(f['common_planar_list_control'] is not None for f in abstract),
               'boundary_sparse_commutators': boundary['sparse_commuting_operator_checks'],
               'generic_jet_full_subset_rank_checks': sum(f['all_subset_ranks_checked'] for f in jets['full_subset_controls']),
               'exceptional_length_drop_controls': 1}
    sources = ['scripts/run_nested_prefix_checks.py', 'docs/nested_prefix_bound.md',
               'scripts/run_nested_jet_checks.py', 'results/research/nested_jet_checks.json',
               'scripts/run_structural_complexity_checks.py', 'scripts/run_extremal_initial_ideals_checks.py',
               'scripts/run_initial_cone_height_checks.py', 'scripts/run_interleaving_research.py']
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_nested_prefix_checks.py',
              'scope': 'Uniform prefix saturation and direction bounds for matroids on full additive-valuation sublevels whose bases permit lowering valuation, with the entire explicit nested polynomial family and generic fixed-order Taylor jets as applications. All valuations and true quotient lengths are covered by the proof; finite controls do not certify novelty or the unrestricted G_3 endpoint.',
              'summary': summary, 'matrix_hypothesis_controls': hypotheses,
              'saved_complete_fan_controls': fans, 'abstract_flag_controls': abstract,
              'two_cone_boundary_control': boundary,
              'generic_fixed_order_jet_controls': jets,
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              'entire_explicit_nested_family_bound': 'O_d(n^(d-1)*log(n+1)^((d-1)^2)), with actual length n',
              'entire_nested_family_power_separation_excluded': True,
              'generic_fixed_order_jets_power_separation_excluded': True,
              'unrestricted_endpoint_achieved': False, 'novelty_established': False, 'goal_complete': False}
    (ROOT/'results/research/nested_prefix_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', **summary}, indent=2))


if __name__ == '__main__':
    main()
