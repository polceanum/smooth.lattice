#!/usr/bin/env python3
"""Exact controls for the classical greedy/logarithmic ray reduction."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import comb, prod
from pathlib import Path
from random import Random
import hashlib
import json

from run_structural_complexity_checks import area2, clip, cross, dot, greedy, primitive
from run_extremal_initial_ideals_checks import adjacent_weights, cone_for_basis, rectangle_model

ROOT = Path(__file__).resolve().parents[1]


def determinant_mod(columns, p):
    n = len(columns)
    if not n:
        return 1
    rows = [[columns[j][i] % p for j in range(n)] for i in range(n)]
    value = 1
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            rows[pivot], rows[j] = rows[j], rows[pivot]
            value = -value
        x = rows[j][j]
        value = value*x % p
        inv = pow(x, -1, p)
        for i in range(j + 1, n):
            scale = rows[i][j]*inv % p
            for k in range(j, n):
                rows[i][k] = (rows[i][k] - scale*rows[j][k]) % p
    return value % p


def complete_fan(colors, factors, bases, rank):
    """Enumerate every basis, then its full positive cone by all competitors."""
    assert len(colors) == len(factors) and min(factors) > 0
    ground = [tuple(a if i == c else 0 for i in range(3))
              for c, a in zip(colors, factors)]
    profiles = [tuple(sum(ground[e][i] for e in b) for i in range(3)) for b in bases]
    profile_by_basis = dict(zip(bases, profiles))
    least_product, product_minimizers = {}, defaultdict(set)
    for b in bases:
        c = tuple(sum(colors[e] == i for e in b) for i in range(3))
        q = prod(factors[e] for e in b)
        if c not in least_product or q < least_product[c]:
            least_product[c], product_minimizers[c] = q, {b}
        elif q == least_product[c]:
            product_minimizers[c].add(b)
    states, owners, total_area, exhaustive_areas, seen_profiles = [], {}, Fraction(0), {}, set()
    for b, p in zip(bases, profiles):
        if p in seen_profiles:
            continue
        seen_profiles.add(p)
        polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
        for q in profiles:
            if q != p:
                polygon = clip(polygon, tuple(x-y for x, y in zip(q, p)), 1)
                if not polygon:
                    break
        if not polygon or not area2(polygon):
            continue
        w = tuple(sum(v[i] for v in polygon)/len(polygon) for i in range(3))
        assert min(w) > 0
        values = [dot(w, q) for q in profiles]
        assert dot(w, p) == min(values) and all(v > min(values) for v, q in zip(values, profiles) if q != p)
        # Independently replay greedy using extendibility to the exhaustive bases.
        chosen = frozenset()
        for e in sorted(range(len(ground)), key=lambda e: (dot(w, ground[e]), e)):
            candidate = chosen | {e}
            if any(candidate <= other for other in bases):
                chosen = candidate
        assert profile_by_basis[chosen] == p
        b = chosen
        c = tuple(sum(colors[e] == i for e in b) for i in range(3))
        assert c not in owners, 'one positive optimum per ray-count vector'
        assert b in product_minimizers[c]
        assert all(profile_by_basis[other] == p for other in product_minimizers[c]), 'exact product comparison within its count fiber'
        owners[c] = b
        exhaustive_areas[b] = area2(polygon)
        total_area += area2(polygon)
        states.append({'basis': sorted(b), 'ray_counts': c, 'profile': p,
                       'strict_weight': [str(x) for x in w],
                       'product_coefficient': least_product[c]})
    assert total_area == 1, 'complete positive simplex covered by exhaustive basis cones'
    assert len(states) <= comb(rank + 2, 2)
    # A separate traversal uses single exchanges and at most six ray-pair
    # inequalities, rather than every competing basis-profile inequality.
    base_set = set(bases)
    comparison = sorted({primitive(tuple(x-y for x, y in zip(a, b)))
                         for a, b in combinations(ground, 2) if a != b})
    scale = max(factors) + 1
    pending, visited, calls, tests = [(1, scale, scale*scale)], {}, 0, 0
    maximum_facets = 0
    while pending:
        w = pending.pop()
        chosen = frozenset()
        for e in sorted(range(len(ground)), key=lambda e: (dot(w, ground[e]), e)):
            candidate = chosen | {e}
            if any(candidate <= other for other in bases):
                chosen = candidate
        calls += 1
        if chosen in visited:
            continue
        tight = {}
        for e in chosen:
            for f in set(range(len(ground))) - chosen:
                tests += 1
                if chosen - {e} | {f} not in base_set:
                    continue
                i, j = colors[e], colors[f]
                if i == j:
                    assert factors[e] <= factors[f]
                    continue
                ratio = Fraction(factors[f], factors[e])
                if (i, j) not in tight or ratio < tight[i, j][0]:
                    tight[i, j] = (ratio, tuple(x-y for x, y in zip(ground[f], ground[e])))
        normals = [normal for _, normal in tight.values()]
        assert len(normals) <= 6 and all(dot(w, z) > 0 for z in normals)
        polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
        for z in normals:
            polygon = clip(polygon, z, 1)
        assert polygon and area2(polygon) == exhaustive_areas[chosen]
        assert all(dot(v, q) >= dot(v, profile_by_basis[chosen])
                   for v in polygon for q in profiles), 'every traversal-cone vertex satisfies all competing-basis inequalities'
        visited[chosen] = polygon
        neighbors = adjacent_weights(polygon, normals, comparison)
        maximum_facets = max(maximum_facets, len(neighbors))
        assert len(neighbors) <= 6
        pending.extend(neighbors)
    assert set(visited) == set(exhaustive_areas)
    assert calls <= 1 + 6*len(visited)
    return {'rank': rank, 'ground_size': len(colors), 'independent_bases': len(bases),
            'positive_states': len(states), 'composition_bound': comb(rank + 2, 2),
            'simplex_double_area': str(total_area), 'states': states,
            'independent_traversal_support_calls': calls, 'single_exchange_queries': tests,
            'maximum_interior_facets': maximum_facets}


def abstract_controls():
    counts, ledger, fixtures = defaultdict(int), [], []

    def add(label, colors, factors, bases, rank, save=False):
        row = complete_fan(colors, factors, bases, rank)
        counts['complete_abstract_fans'] += 1
        counts['exhaustive_independent_bases'] += len(bases)
        counts['positive_profile_and_product_checks'] += row['positive_states']
        counts['independent_traversal_support_calls'] += row['independent_traversal_support_calls']
        counts['single_exchange_queries'] += row['single_exchange_queries']
        ledger.append([label, row])
        if save:
            fixtures.append({'model': label, **row})

    # All 813 labelled rank-two matroids on six elements, including loops.
    def labels(prefix, maximum):
        if len(prefix) == 6:
            if maximum >= 2:
                yield tuple(prefix)
        else:
            for c in range(maximum + 2):
                yield from labels(prefix + [c], max(c, maximum))
    number = 0
    for parallel in labels([], 0):
        bases = [frozenset((i, j)) for i, j in combinations(range(6), 2)
                 if parallel[i] and parallel[j] and parallel[i] != parallel[j]]
        add('rank_two_' + ''.join(map(str, parallel)), (0, 0, 1, 1, 2, 2),
            (1, 3, 2, 5, 1, 7), bases, 2)
        number += 1
    assert number == 813
    counts['all_labelled_rank_two_matroids'] = number
    add('uniform_rank_three_with_repeated_ray_weights', (0, 0, 0, 1, 1, 2, 2),
        (1, 1, 3, 1, 2, 1, 4), [frozenset(b) for b in combinations(range(7), 3)], 3, save=True)
    colors, factors = (0, 0, 0, 1, 1, 1, 2, 2, 2), (1, 2, 7, 1, 4, 9, 2, 3, 11)
    for seed in range(5):
        random = Random(700 + seed)
        columns = [tuple(int(i == j) for i in range(3)) if j < 3
                   else tuple(random.randrange(101) for i in range(3)) for j in range(9)]
        bases = [frozenset(b) for b in combinations(range(9), 3)
                 if determinant_mod([columns[i] for i in b], 101)]
        add('rank_three_linear_' + str(seed), colors, factors, bases, 3, save=seed == 0)
    pairs = [frozenset((2*i, 2*i+1)) for i in range(4)]
    excluded = {pairs[i] | pairs[j] for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3))}
    bases = [frozenset(b) for b in combinations(range(8), 4) if frozenset(b) not in excluded]
    base_set = set(bases)
    for b, other in product(bases, repeat=2):
        for e in b - other:
            assert any(b - {e} | {f} in base_set for f in other - b)
    add('rank_four_sparse_paving', (0, 0, 0, 1, 1, 2, 2, 2),
        (1, 2, 5, 1, 4, 1, 3, 7), bases, 4, save=True)
    direct = []
    rows = ((2, 3, 5), (7, 11, 13), (17, 19, 23), (29, 31, 37))
    for r in (2, 3, 4):
        colors = tuple(i for _ in range(r) for i in range(3))
        factors = tuple(a for row in rows[:r] for a in row)
        bases = [frozenset(3*j + c for j, c in enumerate(choice))
                 for choice in product(range(3), repeat=r)]
        row = complete_fan(colors, factors, bases, r)
        assert row['positive_states'] < 3**r
        direct.append({'component_ranks': [1]*r, 'product_of_separate_counts': 3**r, **row})
        counts['direct_sum_common_fans'] += 1
    return {'counts': dict(counts), 'ledger_sha256': hashlib.sha256(json.dumps(ledger, sort_keys=True).encode()).hexdigest(),
            'exact_abstract_fixtures': fixtures, 'direct_sum_witnesses': direct}


def compositions(n, d):
    if d == 1:
        yield (n,)
    else:
        for a in range(n + 1):
            for rest in compositions(n-a, d-1):
                yield (a,) + rest


def pole_controls():
    counts, fixtures = defaultdict(int), []
    p = 1009
    for d, n in product((2, 3, 4), (1, 2, 4, 7, 10)):
        poles, nodes = tuple(range(d-1)), tuple(range(d + 1, d + n + 1))
        points = [(t,) + tuple(pow(t-a, -1, p) for a in poles) for t in nodes]
        axis = [tuple(j if k == i else 0 for k in range(d))
                for i in range(d) for j in range(1, n)]
        cols = {a: tuple(prod(pow(x, e, p) for x, e in zip(point, a)) % p for point in points)
                for a in [(0,)*d] + axis}
        seen = set()
        for c in compositions(n-1, d):
            basis = [(0,)*d] + [a for a in axis if next(x for x in a if x) <= c[next(i for i, x in enumerate(a) if x)]]
            assert len(basis) == n and determinant_mod([cols[a] for a in basis], p)
            w = tuple(Fraction(2, 2*a + 1) for a in c)
            assert all(dot(w, a) < 1 for a in basis)
            assert all(dot(w, a) > 1 for a in axis if a not in basis)
            seen.add(c)
            counts['all_pole_allocation_determinants'] += 1
        assert len(seen) == comb(n+d-2, d-1)
        for point in points:
            for i, a in enumerate(poles, 1):
                assert (point[0]*point[i] - 1 - a*point[i]) % p == 0
            for i, j in combinations(range(1, d), 2):
                assert ((poles[i-1]-poles[j-1])*point[i]*point[j] - point[i]+point[j]) % p == 0
            counts['partial_fraction_point_checks'] += 1
        counts['sharp_pole_models'] += 1
        fixtures.append({'dimension': d, 'length': n, 'characteristic': p, 'nodes': nodes, 'finite_poles': poles,
                         'all_initial_allocations': len(seen), 'sharp_bound': 'binomial(n+d-2,d-1)'})
    return {'counts': dict(counts), 'exact_models': fixtures}


def mixed_control():
    p = 101
    coefficients = {(i, j): ((i+1)**2+3*(j+1)+5*i*j) % p
                    for i, j in product(range(2), range(3)) if i+j}
    ground, columns, residue = rectangle_model(2, 3, coefficients, p)
    common = {(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 0), (1, 0, 0)}
    records = []
    for mixed, w in (((0, 1, 1), (20, 15, 13)), ((1, 0, 1), (31, 51, 38))):
        basis = common | {mixed}
        assert greedy(ground, columns, w, p) == basis
        assert determinant_mod([columns[a] for a in sorted(basis)], p)
        polygon, normals = cone_for_basis(basis, residue, p)
        assert area2(polygon) > 0 and all(dot(w, z) > 0 for z in normals)
        axis_counts = tuple(sum(a[i] > 0 and sum(x > 0 for x in a) == 1 for a in basis) for i in range(3))
        records.append({'basis': sorted(basis), 'strict_weight': w, 'pure_axis_counts': axis_counts,
                        'profile': [sum(a[i] for a in basis) for i in range(3)],
                        'strict_border_comparisons': len(normals)})
    assert records[0]['pure_axis_counts'] == records[1]['pure_axis_counts'] == (1, 1, 2)
    assert records[0]['profile'] != records[1]['profile']
    return {'length': 6, 'characteristic': p, 'ideal': '(x^2,y^3,z-F)',
            'F_terms': [{'exponents': a, 'coefficient': b} for a, b in sorted(coefficients.items())],
            'strict_initial_witnesses': records,
            'scope': 'An existing actual quotient has different initial ideals with identical pure-axis counts. Mixed monomials invalidate a general extension of the ray-count certificate; this is not an asymptotic lower bound.'}


def margin_controls():
    """Historical diagnostics for a packing premise now refuted by actual quotients."""
    old = json.loads((ROOT/'results/research/extremal_initial_ideals_checks.json').read_text())
    checked, witnesses = 0, []
    for model in old['rectangular_complete_intersections']:
        n = model['quotient_length']
        scale = Fraction(1, 16*n)
        for state in model['witnesses']:
            polygon = [tuple(map(Fraction, v)) for v in state['cone_polygon']]
            center = tuple(sum(v[i] for v in polygon)/len(polygon) for i in range(3))
            normals = []
            for a, b in zip(polygon, polygon[1:] + polygon[:1]):
                z = cross(a, b)
                if dot(z, center) < 0:
                    z = tuple(-x for x in z)
                assert dot(z, center) > 0
                normals.append(z)
            eroded = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
            for z in normals:
                eroded = clip(eroded, tuple(x-scale*abs(x) for x in z), 1)
            assert eroded and area2(eroded) > 0
            w = tuple(sum(v[i] for v in eroded)/len(eroded) for i in range(3))
            assert min(w) > 0 and all(dot(z, w) > scale*sum(abs(x)*y for x, y in zip(z, w)) for z in normals)
            checked += 1
        witnesses.append({'length': n, 'complete_initial_cells': len(model['witnesses']),
                          'relative_margin_checked': str(scale)})
    hard = []
    for k in (3, 4, 10, 100, 10000, 10**12):
        n, w = 2*k+1, (4*k+1, 4*k, 4*k*k-2*k)
        normals = ((1, -1, 0), (0, k, -1), (-k, 1, 1))
        slacks = [Fraction(dot(z, w), sum(abs(x)*y for x, y in zip(z, w))) for z in normals]
        assert slacks == [Fraction(1, 8*k+1), Fraction(1, 4*k-1), Fraction(1, 8*k+3)]
        assert min(slacks) >= Fraction(1, 4*n)
        hard.append({'k': k, 'length': n, 'strict_weight': w,
                     'relative_slacks': list(map(str, slacks)), 'uniform_lower_margin': str(Fraction(1, 4*n))})
    needles = []
    for M in (2, 3, 5, 10, 100, 10000):
        n = 10*(M+5)
        lower, upper = Fraction(M, M+1), Fraction(M+1, M+2)
        lower_z, upper_z = Fraction(M+3, M+4), Fraction(M+4, M+5)
        assert lower < upper and lower_z < upper_z
        s = Fraction(1, 16*n)
        # Eroding both x/y inequalities scales their ratio bounds in opposite directions.
        ratio = (1+s)/(1-s)
        fails = lower*ratio >= upper/ratio
        if M >= 100:
            assert fails
        # Four binary parallel classes plus n-4 distinct axis coloops.
        pairs = (((M+1, 0, 0), (0, M, 0)),
                 ((M+2, 0, 0), (0, M+1, 0)),
                 ((M+4, 0, 0), (0, 0, M+3)),
                 ((M+5, 0, 0), (0, 0, M+4)))
        occupied = {a for pair in pairs for a in pair}
        assert len(occupied) == 8
        available = [tuple(j if i == c else 0 for i in range(3))
                     for c in range(3) for j in range(1, n) if tuple(j if i == c else 0 for i in range(3)) not in occupied]
        coloops = available[:n-4]
        assert len(coloops) == n-4
        assert all(prod(x+1 for x in a) <= n for a in occupied | set(coloops))
        selected = set(coloops) | {pairs[0][1], pairs[1][0], pairs[2][1], pairs[3][0]}
        assert len(selected) == n and (M+2, 0, 0) in selected and (M+1, 0, 0) not in selected
        w = (Fraction(1), 2/(lower+upper), 2/(lower_z+upper_z))
        assert all(dot(w, chosen) < dot(w, other) for chosen, other in
                   ((pairs[0][1], pairs[0][0]), (pairs[1][0], pairs[1][1]),
                    (pairs[2][1], pairs[2][0]), (pairs[3][0], pairs[3][1])))
        needles.append({'M': M, 'rank_and_hyperbolic_bound': n, 'parallel_pairs': pairs,
                        'coloop_count': n-4, 'selected_basis_is_downset': False,
                        'margin_one_over_16n_fails': fails,
                        'log_cell_area_upper': str(Fraction(1, M*(M+2)*(M+3)*(M+5)))})
    return {'complete_actual_quotient_cells_checked': checked, 'quotient_models': witnesses,
            'infinite_hard_cone_margin_witnesses': hard, 'abstract_matroid_needles': needles,
            'scope': 'These historical feasibility diagnostics remain valid, but the uniform premise has since been refuted by the actual all-field socle families in docs/initial_cone_stability.md. Saved rectangle fans pass margin 1/(16n), and the older interior hard-cone family has a 1/(4n) witness. The abstract matroid needles are retained as earlier rank/lattice-only obstructions; the stronger actual-quotient refutation is checked separately.'}


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    abstract, poles, mixed = abstract_controls(), pole_controls(), mixed_control()
    margins = margin_controls()
    sources = ('scripts/run_ray_scaling_checks.py', 'docs/ray_scaling_bound.md',
               'scripts/run_structural_complexity_checks.py', 'scripts/run_extremal_initial_ideals_checks.py',
               'results/research/extremal_initial_ideals_checks.json', 'docs/initial_cone_stability.md')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_ray_scaling_checks.py',
              'counts': {**abstract['counts'], **poles['counts'],
                         'relative_margin_actual_initial_cones': margins['complete_actual_quotient_cells_checked'],
                         'hard_relative_margin_identities': len(margins['infinite_hard_cone_margin_witnesses']),
                         'abstract_relative_margin_needles': len(margins['abstract_matroid_needles']),
                         'abstract_margin_one_over_16n_failures': sum(r['margin_one_over_16n_fails'] for r in margins['abstract_matroid_needles'])},
              'abstract_ray_matroids': abstract, 'sharp_pole_quotients': poles, 'mixed_quotient_obstruction': mixed,
              'conditional_relative_margin_route': margins,
              'scope': 'Classical greedy invariance gives a sharp positive common-ray ceiling and bounded-facet traversal; the partition case matches published tropical coarse types. The relative-margin packing implication remains conditional, and its universal Omega(1/n) premise is now refuted by actual all-field quotient families. Earlier finite diagnostics are preserved as historical controls. No unrestricted exponent, novelty certificate or completed goal is claimed.',
              'uniform_relative_margin_hypothesis_proved': False,
              'uniform_relative_margin_hypothesis_status': 'refuted_by_actual_all_field_quotients',
              'novelty_established': False, 'goal_complete': False,
              'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}}
    (ROOT/'results/research/ray_scaling_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
