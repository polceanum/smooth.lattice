#!/usr/bin/env python3
"""Exact controls for Hilbert budgets, full graded fans and tower relations."""
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from hashlib import sha256
from itertools import combinations, product
from math import ceil, gcd
from pathlib import Path
import json

from run_cyclotomic_price_checks import (border_certificate, cyclotomic_model,
    frontier_basis, insert, scalar, sharp_model, solver, subtract, tensor, tower_basis,
    univariate_powers)
from run_planar_dyadic_bound_checks import convex_hull
from run_rectangle_rank_checks import rectangle
from run_structural_complexity_checks import dot, greedy, hyperbolic_cross, positive_chambers, primitive

ROOT = Path(__file__).resolve().parents[1]


def require(test, message):
    if not test:
        raise AssertionError(message)


@cache
def rectangle_counts(p, q, a, b):
    return Counter(p*i+q*j for i in range(a+1) for j in range(b+1))


def rectangles(h, p, q):
    n = sum(h.values())
    return [(a, b) for a in range(n) for b in range(n//(a+1))
            if all(value <= h.get(level, 0) for level, value in rectangle_counts(p, q, a, b).items())]


def rectangle_injection(h, p, q, cells):
    common = gcd(p, q)
    p, q = p//common, q//common
    images = [set(), set()]
    for a, b in cells:
        if p*a >= q*b:
            k, r = divmod(b, p)
            level, side = p*a+q*r, 0
            points = [(a-q*j, r+p*j) for j in range(k+1)]
        else:
            k, r = divmod(a, q)
            level, side = q*b+p*r, 1
            points = [(r+q*j, b-p*j) for j in range(k+1)]
        require(all(0 <= i <= a and 0 <= j <= b and p*i+q*j == level for i, j in points),
                'all injection witnesses are actual divisors in one degree')
        require(k < h.get(common*level, 0), 'injection slot is inside its Hilbert budget')
        require((level, k) not in images[side], 'each half-map is injective')
        images[side].add((level, k))
    require(len(cells) <= 2*sum(h.values()), 'weighted two-positive rectangle bound')
    return len(cells)


def compositions(n, width):
    if width == 1:
        yield (n,)
    else:
        for first in range(n+1):
            for rest in compositions(n-first, width-1):
                yield (first,)+rest


def injection_controls():
    cases = rectangles_checked = zero_cases = 0
    worst = {}
    for n in range(1, 8):
        for values in compositions(n, 6):
            h = dict(enumerate(values))
            for p, q in product(range(1, 5), repeat=2):
                cells = rectangles(h, p, q)
                rectangle_injection(h, p, q, cells)
                cases += 1
                rectangles_checked += len(cells)
                key = p, q
                if key not in worst or Q(len(cells), n) > Q(worst[key]['rectangles'], worst[key]['n']):
                    worst[key] = {'weights': key, 'n': n, 'hilbert': values,
                                  'rectangles': len(cells), 'witnesses': cells}
            for p in range(1, 5):
                cells = rectangles(h, 0, p)
                require(len(cells) <= sum(h.get(p*j, 0) for j in range(6)), 'zero-positive rectangle bound')
                zero_cases += 1
    return {'arbitrary_hilbert_sequences': 'All weak compositions of n<=7 into six degrees; zeros and nonrealizable sequences included.',
            'positive_weight_cases': cases, 'rectangle_injection_witnesses': rectangles_checked,
            'zero_positive_cases': zero_cases, 'largest_observed_ratios': list(worst.values())}


def ambient_quotient(domain, relations, grading, p, label):
    domain = tuple(sorted(domain))
    index = {a: i for i, a in enumerate(domain)}
    require((0, 0, 0) in index, 'proper monomial ambient algebra')
    require(all(tuple(x-(i == j) for j, x in enumerate(a)) in index
                for a in domain for i in range(3) if a[i]), 'ambient monomial downset')
    require(all(len({dot(grading, a) for a in f}) == 1 for f in relations), 'homogeneous defining relations')
    pivots = {}
    for f in relations:
        for shift in domain:
            vector = [0]*len(domain)
            for a, coefficient in f.items():
                target = tuple(x+y for x, y in zip(a, shift))
                if target in index:
                    vector[index[target]] += coefficient
            insert(pivots, vector, p)
    free = [i for i in range(len(domain)) if i not in pivots]

    def reduce(vector):
        value = {i: scalar(c, p) for i, c in enumerate(vector) if scalar(c, p)}
        for i, row in sorted(pivots.items()):
            if i in value:
                subtract(value, row, value[i], p)
        return tuple(value.get(i, 0) for i in free)

    columns = {a: reduce([int(i == index[a]) for i in range(len(domain))]) for a in domain}
    zero = (0,)*len(free)

    def residue(a):
        return columns.get(tuple(a), zero)

    # Independently check that the relation space is closed under every
    # variable, rather than treating an arbitrary vector matroid as an algebra.
    closure = 0
    for row in pivots.values():
        for axis in range(3):
            shifted = [0]*len(domain)
            for i, coefficient in row.items():
                a = tuple(x+(axis == j) for j, x in enumerate(domain[i]))
                if a in index:
                    shifted[index[a]] += coefficient
            require(not any(reduce(shifted)), 'relation space is multiplication closed')
            closure += 1
    h = Counter(dot(grading, domain[i]) for i in free)
    require(all(not c or dot(grading, a) == dot(grading, domain[i])
                for a, column in columns.items() for i, c in zip(free, column)), 'quotient coordinates preserve degrees')
    return {'label': label, 'n': len(free), 'p': p, 'grading': grading, 'h': dict(h),
            'ground': domain, 'residue': residue, 'relations': relations,
            'closure_checks': closure}


def positive_weight(rho, p, q):
    A, B = rho.numerator, rho.denominator
    beta = next(b for b in range(1, p+1) if (A+q*b) % p == 0)
    gamma = (A+q*beta)//p
    if gamma <= 0:
        k = (-gamma)//q+1
        beta, gamma = beta+p*k, gamma+q*k
    result = (B, beta, gamma)
    require(min(result) > 0 and Q(p*gamma-q*beta, B) == rho, 'positive integer weight realizes its exact parameter')
    return result


def samples(points, extreme):
    points = sorted(set(points))
    if not points:
        return [Q(0)]
    return [-Q(extreme)]+[Q(a.numerator+b.numerator, a.denominator+b.denominator)
                            for a, b in zip(points, points[1:])]+[Q(extreme)]


def hilbert_parameters(h, p, q):
    n, D = sum(h.values()), h.get(0, 0)
    M = (n-1)//max(p, q)
    if not M:
        return set(), [(1, 1, 1)], 0
    points = {Q(0)} if M else set()
    count = 0
    for c in range(1, M+1):
        dy = min(h.get(p*j, 0) for j in range(q*c+1))
        dz = min(h.get(q*j, 0) for j in range(p*c+1))
        count += max(0, dy-1)+max(0, dz-1)
        points.update(Q(a, c) for a in range(1, dy))
        points.update(-Q(a, c) for a in range(1, dz))
    require(count <= 2*(n-D), 'linear Hilbert parameter budget')
    weights = [positive_weight(r, p, q) for r in samples(points, max(1, D))]
    require(len(weights) <= 2+count and max(map(max, weights)) <= 4*n,
            'linear common list and linear certificate height')
    return points, weights, count


def directions_from_hilbert(h, grading):
    n, directions, pair_records = sum(h.values()), set(), []
    for i, j in combinations(range(3), 2):
        k = next(k for k in range(3) if k not in (i, j))
        cells = rectangles(h, grading[i], grading[j])
        rectangle_injection(h, grading[i], grading[j], cells)
        pair_records.append({'coordinates': (i, j), 'admissible_rectangles': cells})
        for a, b in cells:
            if not a or not b:
                continue
            total = grading[i]*a+grading[j]*b
            if total % grading[k]:
                continue
            v = [0, 0, 0]
            v[i], v[j], v[k] = a, b, -total//grading[k]
            directions.add(primitive(v))
        v = [0, 0, 0]
        factor = gcd(grading[i], grading[j])
        v[i], v[j] = grading[j]//factor, -grading[i]//factor
        directions.add(primitive(v))
    require(len(directions) <= 6*n+3, 'linear unoriented edge-direction budget')
    return directions, pair_records


def rank_basis(model, basis):
    pivots = {}
    for a in basis:
        insert(pivots, model['residue'](a), model['p'])
    return len(pivots)


def exchanged_pair(model, left, right, weight):
    for a in sorted(left-right):
        for b in sorted(right-left):
            if dot(weight, a) != dot(weight, b):
                continue
            if (rank_basis(model, left-{a}|{b}) == model['n'] and
                    rank_basis(model, right-{b}|{a}) == model['n']):
                require(dot(model['grading'], a) == dot(model['grading'], b), 'exchange has the same degree')
                return a, b
    raise AssertionError('a symmetric optimal exchange must exist')


def full_fan(model):
    grading, p, n, residue = model['grading'], model['p'], model['n'], model['residue']
    active = [a for a in model['ground'] if any(residue(a))]
    columns = {a: residue(a) for a in active}
    normals = {primitive(tuple(x-y for x, y in zip(a, b)))
               for a, b in combinations(active, 2) if dot(grading, a) == dot(grading, b)
               and min(x-y for x, y in zip(a, b)) < 0 < max(x-y for x, y in zip(a, b))}
    if grading[0] == 0:
        gp, gq = grading[1:]
        require(gcd(gp, gq) == 1, 'normalized support-two grading')
        points = sorted({Q(-gp*v[0], v[2]) for v in normals if v[2]})
        rho = samples(points, max([abs(t) for t in points]+[0])+1)
        weights = [positive_weight(t, gp, gq) for t in rho]
        allowed, common, budget = hilbert_parameters(model['h'], gp, gq)
    else:
        points, rho = None, None
        weights = positive_chambers(sorted(normals))
        allowed, pair_records = directions_from_hilbert(model['h'], grading)
    unique, interval_bases = {}, []
    for weight in weights:
        dense = frozenset(greedy(active, columns, weight, p))
        direct = frontier_basis(residue, n, weight, p)
        require(len(dense) == n and dense == direct, 'dense and frontier computations agree')
        counts = Counter(dot(grading, a) for a in dense)
        require(dict(counts) == model['h'], 'each initial staircase has the quotient Hilbert function')
        interval_bases.append(dense)
        if dense not in unique:
            unique[dense] = {'weight': weight, 'standard_set': sorted(dense),
                             'strict_border_relations': border_certificate(dense, residue, weight, p)}
    exchanges = []
    if grading[0] == 0:
        for t, left, right in zip(points, interval_bases, interval_bases[1:]):
            if left == right:
                continue
            require(t in allowed, 'every actual wall is in the Hilbert list')
            wall_weight = positive_weight(t, gp, gq)
            a, b = exchanged_pair(model, left, right, wall_weight)
            exchanges.append({'rho': str(t), 'from': a, 'to': b})
        common_bases = {frontier_basis(residue, n, w, p) for w in common}
        require(common_bases == set(unique), 'Hilbert-dependent short weight list covers the independent full fan')
        require(len(unique) <= 2+budget, 'global support-two linear count')
        extra = {'hilbert_parameter_candidates': [str(t) for t in sorted(allowed)],
                 'common_short_weights': common, 'parameter_budget': budget}
    else:
        by_profile = {(sum(a[0] for a in basis), sum(a[1] for a in basis)): basis for basis in unique}
        require(len(by_profile) == len(unique), 'exposed profiles identify their standard sets')
        polygon = convex_hull(list(by_profile))
        edges = list(zip(polygon, polygon[1:]+polygon[:1])) if len(polygon) > 1 else []
        for u, v in edges:
            dx, dy = v[0]-u[0], v[1]-u[1]
            objective = (-dy, dx, 0) if len(polygon) > 2 else (0, 0, 0)
            C = max([ceil(Q(-x, d)) for x, d in zip(objective, grading)]+[0])+1
            wall_weight = tuple(x+C*d for x, d in zip(objective, grading))
            left, right = by_profile[u], by_profile[v]
            require(all(sum(dot(wall_weight, a) for a in basis) >= sum(dot(wall_weight, a) for a in left)
                        for basis in unique), 'polygon edge objective is globally minimizing')
            a, b = exchanged_pair(model, left, right, wall_weight)
            direction = primitive(tuple(x-y for x, y in zip(a, b)))
            require(direction in allowed, 'actual edge exchange is Hilbert admissible')
            exchanges.append({'edge_profiles': (u, v), 'from': a, 'to': b, 'primitive_direction': direction})
        require(len(unique) <= 12*n+6, 'global positive-grading linear bound')
        extra = {'profile_polygon': polygon, 'candidate_edge_directions': sorted(allowed),
                 'rectangle_budget_controls': pair_records}
    report = {'label': model['label'], 'field': p, 'grading': grading, 'n': n,
              'hilbert_function': model['h'], 'ambient_domain': model['ground'],
              'ideal_relations': [[[*a, str(c)] for a, c in sorted(f.items())] for f in model.get('relations', [])],
              'relation_closure_checks': model.get('closure_checks', 0),
              'complete_comparison_chambers': len(weights), 'initial_ideal_count': len(unique),
              'full_positive_comparison_coverage': True, 'witnesses': list(unique.values()),
              'endpoint_exchange_controls': exchanges, **extra}
    return report


def graded_fixtures():
    fixtures = []
    for p in (0, 2, 3):
        data = [
            ('nonfree_equal', (0, 1, 1),
             [(a, b, c) for a in range(4) for b in range(4) for c in range(4)
              if b+c <= 3 and a < 4-b-c],
             [{(0, 2, 0): 1, (1, 0, 2): 1}, {(2, 1, 0): 1, (1, 0, 1): -1}]),
            ('unequal_1_2', (0, 1, 2),
             [(a, b, c) for a in range(3) for b in range(6) for c in range(3) if b+2*c <= 5],
             [{(0, 0, 1): 1, (1, 2, 0): -1}]),
            ('unequal_2_3', (0, 2, 3),
             [(a, b, c) for a in range(2) for b in range(4) for c in range(3) if 2*b+3*c <= 7],
             [{(0, 3, 0): 1, (1, 0, 2): -1}]),
            ('large_grading_constant_fan', (0, 101, 103),
             [(a, b, c) for a in range(2) for b in range(2) for c in range(2) if b+c <= 1], []),
            ('positive_equal', (1, 1, 1),
             [a for a in product(range(4), repeat=3) if sum(a) <= 3],
             [{(2, 0, 0): 1, (0, 1, 1): 1}, {(0, 2, 0): 1, (1, 0, 1): 1},
              {(0, 0, 2): 1, (1, 1, 0): 1}]),
            ('positive_1_2_3', (1, 2, 3),
             [a for a in product(range(8), repeat=3) if dot(a, (1, 2, 3)) <= 7],
             [{(0, 2, 0): 1, (1, 0, 1): 1}, {(0, 0, 2): 1, (0, 3, 0): 1},
              {(3, 0, 0): 1, (0, 0, 1): 1}]),
            ('positive_2_3_5', (2, 3, 5),
             [a for a in product(range(7), repeat=3) if dot(a, (2, 3, 5)) <= 12],
             [{(0, 0, 1): 1, (1, 1, 0): -1}, {(3, 0, 0): 1, (0, 2, 0): 1}])]
        for label, grading, domain, relations in data:
            model = ambient_quotient(domain, relations, grading, p, label)
            record = full_fan(model)
            fixtures.append(record)
            print(json.dumps({'graded_fan_passed': label, 'field': p, 'states': record['initial_ideal_count']}), flush=True)
        for gp, gq in ((5, 7), (13, 17)):
            domain = ([(a, b, 0) for a in range(2) for b in range(gq+1)]
                      +[(a, 0, c) for a in range(2) for c in range(1, gp+1)])
            model = ambient_quotient(domain, [{(0, gq, 0): 1, (1, 0, gp): -1}],
                                      (0, gp, gq), p, 'large_active_grading')
            fixtures.append(full_fan(model))
    return fixtures


def padded_sharp(D, p, odd):
    original = sharp_model(D, p)
    if not odd:
        return {**original, 'p': p}
    coefficients = [0]+original['f_coefficients']
    powers = [(1,)+(0,)*D]
    def degree_zero(k):
        while len(powers) <= k:
            old = powers[-1]
            powers.append(tuple(scalar((old[i-1] if i else 0)-old[-1]*coefficients[i], p)
                                for i in range(D+1)))
        return powers[k]
    def residue(a):
        x, y, z = a
        if y+z == 0:
            return degree_zero(x)+(0,)*D
        if y+z == 1:
            vector = original['residue'](a)
            return (0,)*(D+1)+tuple(vector[2*i+1] for i in range(D))
        return (0,)*(2*D+1)
    return {'n': 2*D+1, 'p': p, 'residue': residue}


def sharpness_controls():
    records = []
    for D, odd, p in product((1, 2, 3, 5), (False, True), (0, 101)):
        model = padded_sharp(D, p, odd)
        model.update({'ground': tuple(hyperbolic_cross(model['n'], 3)), 'grading': (0, 1, 1),
                      'h': {0: D+int(odd), 1: D}, 'label': 'odd_padded_ribbon' if odd else 'even_ribbon'})
        rec = full_fan(model)
        require(rec['initial_ideal_count'] == D+1, 'pointwise linear lower family')
        extreme = next(w for w in rec['witnesses'] if all(a[2] == 0 for a in w['standard_set']))
        require(extreme['weight'][2]-extreme['weight'][1] > (D-1)*extreme['weight'][0],
                'extreme certificate has the required linear height lower bound')
        small_weights = 0
        for weight in product(range(1, D+1), repeat=3):
            strictly_selects = all(dot(weight, relation['head']) > dot(weight, term['exponent'])
                                   for relation in extreme['strict_border_relations'] for term in relation['tail'])
            require(not strictly_selects, 'no smaller integer weight strictly selects the extreme cone')
            small_weights += 1
        witness = (1, 1, D+1)
        require(frontier_basis(model['residue'], model['n'], witness, p) == frozenset(map(tuple, extreme['standard_set'])),
                'the linear integer height lower bound is attained')
        border_certificate(frozenset(map(tuple, extreme['standard_set'])), model['residue'], witness, p)
        rec['D'], rec['odd_padding'], rec['minimum_extreme_weight_height'] = D, odd, D+1
        rec['subcritical_integer_weights_checked'] = small_weights
        records.append(rec)
    return records


def local_refinement_controls():
    records = []
    for side, p, grading in product((2, 3), (2, 3), ((1, 1, 1), (1, 1, 2))):
        original = rectangle(side, side, {(1, 0): 1, (0, 1): 1, (1, 1): 1}, p)
        n = side*side
        L = (2*n+1)**3
        refinement = (1, 2*n+1, (2*n+1)**2)
        weight = tuple(L*d+v for d, v in zip(grading, refinement))
        basis = frontier_basis(original['residue'], n, weight, p)
        borders = border_certificate(basis, original['residue'], weight, p)
        initial = []
        for record in borders:
            poly = {tuple(record['head']): 1}
            for term in record['tail']:
                poly[tuple(term['exponent'])] = -int(term['coefficient'])
            level = max(dot(grading, a) for a in poly)
            require(dot(grading, record['head']) == level, 'chosen order refines the requested weight')
            initial.append({a: c for a, c in poly.items() if dot(grading, a) == level})
        domain = [a for a in product(range(2*side-1), repeat=3) if sum(a) <= 2*side-2]
        model = ambient_quotient(domain, initial, grading, p, 'initial_form_of_bilinear_rectangle')
        require(model['n'] == n, 'initial form preserves the actual quotient length')
        rec = full_fan(model)
        for witness in rec['witnesses']:
            v = witness['weight']
            C = 4*n*max(v)+1
            near = tuple(C*d+x for d, x in zip(grading, v))
            require(frontier_basis(original['residue'], n, near, p) == frozenset(map(tuple, witness['standard_set'])),
                    'ordinary original quotient agrees with the local initial-of-initial refinement')
        records.append({'side': side, 'field': p, 'ray': grading,
                        'local_refinements_checked': rec['initial_ideal_count'], 'initial_form_fan': rec})
    return records


def tower_relation_controls():
    records = []
    for periods, p in product(((3, 5, 7), (3, 5, 11)), (0, 2, 3)):
        model = cyclotomic_model(periods, p)
        D, R, n = (periods[0]-1)*(periods[1]-1), periods[2]-1, model['n']
        powers = [univariate_powers(m, p) for m in periods[:2]]
        def nilpotent_residue(a):
            x, y, z = a
            if y+z >= R:
                return (0,)*n
            return tensor(powers[0](x+y), powers[1](x+z),
                          tuple(int(i == y+z) for i in range(R)))
        for template in ((1, 1, 1), (1, 3, 3), (1, 5, 7), (1, 100, 101)):
            L = 2*n+1
            weight = tuple(x*L**3+L**i for i, x in enumerate(template))
            rho = Q(weight[2]-weight[1], weight[0])
            graded = tower_basis(model, 0, D, R, rho)
            require(frontier_basis(nilpotent_residue, n, weight, p) == graded,
                    'ordinary nilpotent arithmetic verifies the global tower formula without block separation')
            ordinary = frontier_basis(model['residue'], n, weight, p)
            certificate = border_certificate(graded, lambda a: model['residue'](a), weight, p) if ordinary == graded else None
            ordered = sorted(graded)
            solve = solver([model['residue'](a) for a in ordered], p)
            outside = {tuple(x+(axis == j) for j, x in enumerate(a)) for a in graded for axis in range(3)}-graded
            heads = [a for a in outside if all(tuple(x-(axis == j) for j, x in enumerate(a)) in graded
                                             for axis in range(3) if a[axis])]
            failures, reaching = [], []
            for head in sorted(heads):
                coefficients = solve(model['residue'](head))
                bad = [ordered[j] for j in coefficients if dot(weight, head) <= dot(weight, ordered[j])]
                if head[1]+head[2] < R:
                    require(not bad, 'independent lower blocks introduce no deformation failure')
                else:
                    require(head[1]+head[2] == R and head[0] == 0, 'only pure degree-R heads can reach the relation')
                    reaching.append(head)
                    if bad:
                        failures.append({'head': head, 'noncheaper_tails': bad})
            require((ordinary == graded) == (not failures), 'relation-reaching border criterion is exact')
            if ordinary != graded:
                require(any(a[1]+a[2] >= R for a in ordinary), 'every additional state contains a relation-reaching monomial')
            records.append({'periods': periods, 'field': p, 'n': n, 'R': R, 'weight': weight,
                            'graded_standard_set': sorted(graded), 'deformed_standard_set': sorted(ordinary),
                            'relation_reaching_heads': reaching, 'failed_inequalities': failures,
                            'strict_unchanged_border_relations': certificate})
    require(any(r['failed_inequalities'] for r in records), 'control includes an actual failed extrapolation')
    return records


def main():
    injections = injection_controls()
    fans = graded_fixtures()
    sharpness = sharpness_controls()
    local = local_refinement_controls()
    towers = tower_relation_controls()
    all_fans = fans+sharpness+[m['initial_form_fan'] for m in local]
    sources = ('docs/nonnegative_graded_fans.md', 'scripts/run_graded_fan_checks.py',
               'scripts/run_cyclotomic_price_checks.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_planar_dyadic_bound_checks.py', 'scripts/run_rectangle_rank_checks.py')
    summary = {'arbitrary_positive_hilbert_budget_cases': injections['positive_weight_cases'],
               'explicit_rectangle_injections': injections['rectangle_injection_witnesses'],
               'zero_positive_hilbert_budget_cases': injections['zero_positive_cases'],
               'complete_graded_model_fans': len(fans), 'sharp_even_odd_fans': len(sharpness),
               'subcritical_integer_weights_checked': sum(m['subcritical_integer_weights_checked'] for m in sharpness),
               'complete_local_initial_form_fans': len(local),
               'independent_comparison_chambers': sum(m['complete_comparison_chambers'] for m in all_fans),
               'standard_sets_certified': sum(m['initial_ideal_count'] for m in all_fans),
               'strict_border_relations': sum(len(w['strict_border_relations']) for m in all_fans for w in m['witnesses']),
               'endpoint_symmetric_exchanges': sum(len(m['endpoint_exchange_controls']) for m in all_fans),
               'relation_space_closure_checks': sum(m['relation_closure_checks'] for m in all_fans),
               'ordinary_local_refinements': sum(m['local_refinements_checked'] for m in local),
               'tower_deformation_controls': len(towers),
               'actual_relation_reaching_changes': sum(bool(m['failed_inequalities']) for m in towers)}
    report = {'schema_version': 1, 'command': 'python3 -B scripts/run_graded_fan_checks.py',
              'scope': 'Finite controls of uniform global nonnegative-grading and local-refinement proofs. Complete fans are small homogeneous fixtures and are independently enumerated from all monomial comparisons. No full mixed-period fan search, unrestricted G_3 resolution or novelty certificate.',
              'source_sha256': {s: sha256((ROOT/s).read_bytes()).hexdigest() for s in sources},
              'hilbert_budget_controls': injections, 'actual_graded_quotients': fans,
              'pointwise_sharpness_controls': sharpness, 'local_refinement_controls': local,
              'tower_relation_controls': towers, 'summary': summary}
    (ROOT/'results/research/graded_fan_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(summary, sort_keys=True))


if __name__ == '__main__':
    main()
