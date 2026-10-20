#!/usr/bin/env python3
"""Exact controls for Hilbert budgets, full graded fans and tower relations."""
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from hashlib import sha256
from itertools import combinations, product
from math import ceil, gcd, lcm
from pathlib import Path
import json

from run_cyclotomic_price_checks import (border_certificate, cyclotomic_model,
    frontier_basis, insert, scalar, sharp_model, solver, subtract, tensor, tower_basis,
    univariate_powers)
from run_planar_dyadic_bound_checks import convex_hull
from run_rectangle_rank_checks import rectangle
from run_structural_complexity_checks import clip, dot, greedy, hyperbolic_cross, positive_chambers, primitive

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


def signed_rectangle_injection(h, p, q, cells):
    if p < 0:
        h, p, q = {-k: v for k, v in h.items()}, -p, -q
    if q > 0:
        return rectangle_injection(h, p, q, cells)
    require(p > 0 > q, 'opposite pair signs')
    common = gcd(p, -q)
    p, q = p//common, -q//common
    images = set()
    for a, b in cells:
        level, k = p*a-q*b, min(a//q, b//p)
        divisors = [(a-q*j, b-p*j) for j in range(k+1)]
        require(all(0 <= i <= a and 0 <= j <= b and p*i-q*j == level for i, j in divisors),
                'opposite-sign corner fiber is an actual divisor chain')
        require(k < h.get(common*level, 0), 'corner fiber fits the Hilbert slot')
        require((level, k) not in images, 'opposite-sign corner map is injective')
        images.add((level, k))
    require(len(cells) <= sum(h.values()), 'opposite-sign rectangle budget is n')
    return len(cells)


def signed_injection_controls():
    cases = corners = 0
    largest = {}
    for n in range(1, 7):
        for values in compositions(n, 7):
            h = dict(zip(range(-3, 4), values))
            for p, q in product(range(1, 5), range(-4, 0)):
                cells = rectangles(h, p, q)
                signed_rectangle_injection(h, p, q, cells)
                cases += 1
                corners += len(cells)
                if (p, q) not in largest or Q(len(cells), n) > Q(largest[p, q]['rectangles'], largest[p, q]['n']):
                    largest[p, q] = {'pair_degrees': (p, q), 'n': n, 'hilbert': h,
                                     'rectangles': len(cells), 'corners': cells}
    return {'scope': 'All weak compositions of total n<=6 on degrees -3 through 3, and pair degrees p=1..4,q=-4..-1.',
            'opposite_sign_budget_cases': cases, 'explicit_corner_injections': corners,
            'largest_observed_ratios': list(largest.values())}


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


def signed_directions(h, grading):
    require(all(grading), 'the omitted coordinate must have nonzero degree')
    allowed, records, n = set(), [], sum(h.values())
    for i, j in combinations(range(3), 2):
        k = next(k for k in range(3) if k not in (i, j))
        cells = rectangles(h, grading[i], grading[j])
        signed_rectangle_injection(h, grading[i], grading[j], cells)
        records.append({'coordinates': (i, j), 'admissible_rectangles': cells})
        for a, b in cells:
            total = grading[i]*a+grading[j]*b
            if not a or not b or total % grading[k] or total//grading[k] <= 0:
                continue
            vector = [0, 0, 0]
            vector[i], vector[j], vector[k] = a, b, -total//grading[k]
            allowed.add(primitive(vector))
        common = gcd(grading[i], grading[j])
        vector = [0, 0, 0]
        vector[i], vector[j] = grading[j]//common, -grading[i]//common
        if min(vector) < 0 < max(vector):
            allowed.add(primitive(vector))
    require(len(allowed) <= 6*n+3, 'signed full-support grading direction budget')
    return allowed, records


def signed_full_fan(model, use_linear_theorem=True):
    """Positive cones and their shared facets; no invisible polygon edge is used."""
    g, p, n, residue = model['grading'], model['p'], model['n'], model['residue']
    active = [a for a in model['ground'] if any(residue(a))]
    columns = {a: residue(a) for a in active}
    normals = {primitive(tuple(x-y for x, y in zip(a, b)))
               for a, b in combinations(active, 2) if dot(g, a) == dot(g, b)
               and min(x-y for x, y in zip(a, b)) < 0 < max(x-y for x, y in zip(a, b))}
    weights, unique = positive_chambers(sorted(normals)), {}
    for weight in weights:
        basis = frozenset(greedy(active, columns, weight, p))
        require(basis == frontier_basis(residue, n, weight, p), 'independent signed positive-fan computations')
        require(dict(Counter(dot(g, a) for a in basis)) == model['h'], 'signed Hilbert function is preserved')
        if basis not in unique:
            border = border_certificate(basis, residue, weight, p)
            inequalities = {tuple(x-y for x, y in zip(f['head'], t['exponent']))
                            for f in border for t in f['tail']}
            polygon = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
            for normal in sorted(inequalities):
                polygon = clip(polygon, normal, 1)
                require(polygon, 'strict signed positive cone survives border clipping')
            unique[basis] = {'weight': weight, 'standard_set': sorted(basis),
                             'strict_border_relations': border, 'polygon': polygon,
                             'inequalities': inequalities}
    wall_map = {}
    for basis, record in unique.items():
        polygon = record['polygon']
        for a, b in zip(polygon, polygon[1:]+polygon[:1]):
            middle = tuple((x+y)/2 for x, y in zip(a, b))
            if min(middle) <= 0:
                continue
            normal = next(v for v in record['inequalities'] if dot(v, a) == dot(v, b) == 0)
            key = (primitive(normal), tuple(sorted((a, b))))
            wall_map.setdefault(key, []).append(basis)
    allowed, rectangle_records = signed_directions(model['h'], g) if use_linear_theorem else (None, [])
    adjacency, exchanges = Counter(), []
    for (normal, endpoints), bases in wall_map.items():
        require(len(bases) == 2, 'each interior facet has its two actual positive endpoint bases')
        middle = tuple((x+y)/2 for x, y in zip(*endpoints))
        scale = lcm(*(x.denominator for x in middle))
        weight = tuple(int(x*scale) for x in middle)
        left, right = bases
        a, b = exchanged_pair(model, left, right, weight)
        direction = primitive(tuple(x-y for x, y in zip(a, b)))
        require(direction == normal, 'endpoint exchange is parallel to the actual positive facet')
        if use_linear_theorem:
            require(direction in allowed, 'signed actual direction has a Hilbert-budget encoding')
        adjacency[left] += 1
        adjacency[right] += 1
        exchanges.append({'strict_positive_wall_weight': weight, 'from': a, 'to': b,
                          'primitive_direction': direction})
    require(len(exchanges) == len(unique)-1 and max(adjacency.values(), default=0) <= 2,
            'the positive visible part is a chain')
    if use_linear_theorem:
        require(len(unique) <= 12*n+7, 'global signed full-support linear bound')
    witnesses = []
    for record in unique.values():
        witnesses.append({'weight': record['weight'], 'standard_set': record['standard_set'],
                          'strict_border_relations': record['strict_border_relations'],
                          'positive_cone_polygon': [[str(x) for x in a] for a in record['polygon']]})
    return {'label': model['label'], 'field': p, 'grading': g, 'n': n, 'hilbert_function': model['h'],
            'ambient_domain': model['ground'],
            'relations': [[[*a, str(c)] for a, c in sorted(f.items())] for f in model['relations']],
            'relation_closure_checks': model['closure_checks'],
            'complete_comparison_chambers': len(weights), 'initial_ideal_count': len(unique),
            'full_positive_comparison_coverage': True,
            'linear_theorem_applied': use_linear_theorem, 'witnesses': witnesses,
            'actual_positive_wall_exchanges': exchanges, 'rectangle_budget_controls': rectangle_records,
            'candidate_directions': sorted(allowed) if allowed is not None else None}


def signed_fan_controls():
    records = []
    for p in (0, 2, 3):
        data = [
            ('signed_nonlinear_1_1_minus1', (1, 1, -1), tuple(product(range(3), range(3), range(2))),
             [{(1, 0, 0): 1, (0, 1, 0): -1, (0, 2, 1): -1}]),
            ('signed_nonlinear_1_minus1_2', (1, -1, 2), tuple(product(range(4), range(4), range(2))),
             [{(0, 0, 1): 1, (2, 0, 0): -1, (3, 1, 0): -1}]),
            ('signed_nonlinear_2_minus3_5', (2, -3, 5), tuple(product(range(6), range(3), range(2))),
             [{(0, 0, 1): 1, (4, 1, 0): -1}])]
        for label, grading, domain, relations in data:
            model = ambient_quotient(domain, relations, grading, p, label)
            records.append(signed_full_fan(model))
            print(json.dumps({'signed_full_fan_passed': label, 'field': p}), flush=True)
    return records


def zero_opposite_boundary_controls():
    records = []
    for a, c, p in product((3, 5, 7), (2,), (0, 2, 3)):
        domain = [(i, 0, 0) for i in range(1, a+1)]+[(0, j, k) for j in range(c+1) for k in range(c+1)]
        model = ambient_quotient(domain, [{(a, 0, 0): 1, (0, c, c): -1}],
                                  (0, 1, -1), p, 'zero_opposite_socle_boundary')
        rec = signed_full_fan(model, False)
        require(model['n'] == a+(c+1)**2-1, 'socle identification has the exact stated length')
        expected = {frozenset(domain)-{(a, 0, 0)}, frozenset(domain)-{(0, c, c)}}
        require({frozenset(map(tuple, w['standard_set'])) for w in rec['witnesses']} == expected,
                'the actual all-field boundary family has precisely its two proved states')
        require(rec['actual_positive_wall_exchanges'][0]['primitive_direction'] == primitive((a, -c, -c)),
                'fixed positive pair does not determine the remaining coordinate')
        rec['parameters'], rec['positive_pair'] = (a, c), (c, c)
        records.append(rec)
    return records


def primary_parameter_controls():
    records = []
    for d in range(3, 13):
        generators = [(a, d-a, 0) for a in range(d)]+[(0, 0, d)]
        require(all(a[1] or a[2] for a in generators), 'Notation 2.10 leaves the whole pure-X polynomial subalgebra')
        require(all(not any(all(u <= v for u, v in zip(g, (k, 0, 0))) for g in generators)
                    for k in range(3*d+1)), 'explicit pure-X survival controls')
        finite = [a for a in product(range(d+1), repeat=3) if sum(a) <= d
                  and not any(all(u <= v for u, v in zip(g, a)) for g in generators)]
        require(len(finite) == d*(d+1)*(d+5)//6, 'exact Artinian truncation length, distinct from generator degree')
        records.append({'generator_degree': d, 'generators': generators,
                        'untruncated_colength': 'infinite', 'truncation': '(X,Y,Z)^(d+1)',
                        'truncated_length': len(finite), 'degree_d_dimension': sum(sum(a) == d for a in finite)})
    return {'primary_source': 'https://arxiv.org/pdf/1303.3461',
            'scope': 'Only the exact parameter/finite-length normalization is independently checked. The published generic-initial lower bound and determinant appendix are not independently re-proved.',
            'fixtures': records}


def annihilator(vectors, width, p):
    pivots = {}
    for v in vectors:
        insert(pivots, v, p)
    basis = []
    for free in range(width):
        if free in pivots:
            continue
        vector = [0]*width
        vector[free] = 1
        for i, row in sorted(pivots.items(), reverse=True):
            vector[i] = scalar(-sum(value*vector[j] for j, value in row.items() if i != j), p)
        require(all(not scalar(dot(v, vector), p) for v in vectors), 'annihilator equations hold exactly')
        basis.append(tuple(vector))
    return basis


def graph_cycle_basis(r, s, edges, p):
    all_edges = tuple(product(range(r), range(s)))
    edge_index = {e: i for i, e in enumerate(all_edges)}
    parent = list(range(r+s))
    forest = [[] for _ in parent]
    chords = []
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for edge in sorted(edges):
        u, v = edge[0], r+edge[1]
        a, b = find(u), find(v)
        if a == b:
            chords.append(edge)
        else:
            parent[a] = b
            forest[u].append((v, edge))
            forest[v].append((u, edge))
    components = len({find(i) for i in range(r+s)})
    result = []
    for edge in chords:
        start, end = r+edge[1], edge[0]
        reached, queue = {start: None}, [start]
        for u in queue:
            for v, e in forest[u]:
                if v not in reached:
                    reached[v] = (u, e)
                    queue.append(v)
        require(end in reached, 'fundamental cycle has a forest path')
        vector = [0]*len(all_edges)
        vector[edge_index[edge]] = 1
        v = end
        while v != start:
            u, e = reached[v]
            vector[edge_index[e]] += 1 if u < r else -1
            v = u
        require(all(sum(vector[edge_index[e]] for e in all_edges if e[0] == u) == 0 for u in range(r))
                and all(sum(vector[edge_index[e]] for e in all_edges if e[1] == v) == 0 for v in range(s)),
                'fundamental vectors satisfy oriented incidence equations')
        result.append(tuple(scalar(x, p) for x in vector))
    require(len(result) == len(edges)-(r+s)+components, 'graph cycle dimension includes isolated vertices')
    return result, components


def tensor_prefix_controls():
    records = []
    for periods, p in product(((3, 5, 7), (3, 5, 11)), (0, 2, 3, 7)):
        r, s, t = periods
        model, D = cyclotomic_model(periods, p), (r-1)*(s-1)
        powers = [univariate_powers(m, p) for m in periods[:2]]
        all_edges = set(product(range(r), range(s)))
        for weight, cutoff in (((1, 1, 1), t+1), ((1, 3, 3), 3*t-1)):
            monomials = [(a, b, c) for a in range(cutoff//weight[0]+1)
                         for b in range((cutoff-weight[0]*a)//weight[1]+1)
                         for c in range((cutoff-weight[0]*a-weight[1]*b)//weight[2]+1)]
            monomials.sort(key=lambda a: (dot(weight, a), a))
            groups, labels, direct, position = [{} for _ in range(t)], [set() for _ in range(t)], {}, 0
            checks, different_union = [], None
            for cost in range(cutoff+1):
                while position < len(monomials) and dot(weight, monomials[position]) <= cost:
                    a, b, c = monomials[position]
                    insert(direct, model['residue']((a, b, c)), p)
                    j = (b+c) % t
                    coefficient = tensor(powers[0](a+b), powers[1](a+c))
                    insert(groups[j], coefficient, p)
                    labels[j].add(((a+b) % r, (a+c) % s))
                    position += 1
                component_bases = [list(basis.values()) for basis in groups]
                # Convert sparse row bases to ordinary D-vectors.
                component_bases = [[tuple(row.get(i, 0) for i in range(D)) for row in rows]
                                   for rows in component_bases]
                normals = [annihilator(rows, D, p) for rows in component_bases]
                joined = [v for group in normals for v in group]
                intersection = annihilator(joined, D, p)
                for v in intersection:
                    require(all(not insert(dict(group), v, p) for group in groups), 'intersection vector belongs to every coefficient span')
                dims = [len(group) for group in groups]
                predicted = sum(dims)-len(intersection)
                require(predicted == len(direct), 'tensor prefix formula equals direct ordinary polynomial rank')
                cycles, union = [], set()
                for label_set, dimension in zip(labels, dims):
                    complement = all_edges-label_set
                    cycle_vectors, components = graph_cycle_basis(r, s, complement, p)
                    require(dimension == len(label_set)+1-components, 'coefficient prefix has the established cographic rank')
                    cycles += cycle_vectors
                    union |= complement
                cycle_pivots = {}
                for v in cycles:
                    insert(cycle_pivots, v, p)
                require(D-len(cycle_pivots) == len(intersection), 'sum of complement cycle spaces equals the annihilator sum')
                union_cycles, _ = graph_cycle_basis(r, s, union, p)
                require(len(cycle_pivots) <= len(union_cycles), 'cycle-space sum is a subspace of union cycles')
                if len(cycle_pivots) != len(union_cycles) and different_union is None:
                    different_union = {'cost': cost, 'coefficient_labels': [sorted(g) for g in labels],
                                       'cycle_sum_rank': len(cycle_pivots), 'union_cycle_rank': len(union_cycles)}
                if cost < (t-1)*min(weight[1:]):
                    require(not groups[t-1] and not intersection, 'before the last period level arrives the blocks are independent')
                checks.append({'cost_cutoff': cost, 'monomials': position,
                               'coefficient_dimensions': dims, 'intersection_dimension': len(intersection),
                               'intersection_basis': [[str(x) for x in v] for v in intersection],
                               'cycle_space_sum_rank': len(cycle_pivots), 'direct_polynomial_rank': len(direct)})
            records.append({'periods': periods, 'field': p, 'weight': weight, 'maximum_cost': cutoff,
                            'scope': 'Full finite cost balls, not a sampled monomial ground or a full fan count.',
                            'prefixes': checks, 'first_cycle_union_mismatch': different_union})
    counterexamples = []
    for p in (0, 2, 3):
        trees = [{(u, v) for u, v in product(range(3), repeat=2) if u == a or v == a} for a in (0, 1)]
        cycle_sets = [graph_cycle_basis(3, 3, tree, p)[0] for tree in trees]
        union_cycles = graph_cycle_basis(3, 3, set.union(*trees), p)[0]
        require(not any(cycle_sets) and len(union_cycles) == 3, 'two trees do not generate the cycle space of their union')
        counterexamples.append({'field': p, 'trees': [sorted(tree) for tree in trees],
                                'cycle_space_sum_rank': 0, 'union_cycle_rank': 3})
    require(any(m['first_cycle_union_mismatch'] for m in records), 'actual cost prefixes also distinguish the two cycle operations')
    return {'prefix_models': records, 'cycle_union_counterexamples': counterexamples}


def one_period_region_controls():
    records = []
    for periods, p in product(((3, 5, 7), (3, 5, 11)), (0, 2, 3, 7)):
        r, s, t = periods
        model, D = cyclotomic_model(periods, p), (r-1)*(s-1)
        n, scale = model['n'], 2*model['n']+1
        for template in ((2, 3, 3), (9, 10, 10), (1, 10, 11)):
            # Small distinct perturbations resolve comparisons of degree <= n.
            weight = tuple(v*scale**3+scale**i for i, v in enumerate(template))
            alpha, m = weight[0], min(weight[1:])
            require(t*m > (D-1)*alpha, 'strict one-period price condition')
            basis = frontier_basis(model['residue'], n, weight, p)
            auxiliary = tower_basis(model, 0, D, t, Q(weight[2]-weight[1], alpha))
            require(len(auxiliary) == t*D, 'each of the t auxiliary groups has rank D')
            require(basis <= auxiliary, 'ordinary quotient basis is contained in independent coefficient slots')
            require(len(auxiliary-basis) == D, 'the diagonal relation removes exactly D slots')
            require(all(a < D and b+c < t for a, b, c in basis), 'all actual standard monomials satisfy the one-period range')
            pivots = {}
            for a in auxiliary:
                insert(pivots, model['residue'](a), p)
            require(len(pivots) == n, 'auxiliary columns have full ordinary quotient rank')
            border = border_certificate(basis, model['residue'], weight, p)
            records.append({'periods': periods, 'field': p, 'n': n, 'D': D,
                            'price_template': template, 'strict_weight': weight,
                            'older_complete_block_condition': m > (D-1)*alpha,
                            'one_period_condition': True, 'standard_set': sorted(basis),
                            'auxiliary_slots': sorted(auxiliary), 'omitted_slots': sorted(auxiliary-basis),
                            'strict_border_relations': border})
    require(any(not m['older_complete_block_condition'] for m in records),
            'controls include the wider region beyond the completed-block hypothesis')
    require(any(max(b+c for a, b, c in m['standard_set']) == m['periods'][2]-1 for m in records),
            'controls actually reach the coupled last level')
    return {'scope': 'Exact ordinary polynomial arithmetic and strict border controls at specified prices; not a complete regional fan enumeration.',
            'models': records}


def periodic_nilpotent_tower(D, t, defect, p, translated=False):
    """F=K[x]/x^D, with an independently divided monic factor of W^t-1."""
    one = (1,)+(0,)*(D-1)

    def mul(a, b):
        return tuple(scalar(sum(a[j]*b[i-j] for j in range(i+1)), p) for i in range(D))

    def power(a, k):
        result = one
        for _ in range(k):
            result = mul(result, a)
        return result

    u = tuple(int(i in (0, 1)) for i in range(D))
    v = tuple(int(i in (0, 2)) for i in range(D))
    if translated:
        require(defect == 1 and power(u, t) == one, 'translated factor has period t in F')
        g = [power(u, t-1-i) for i in range(t)]
        factor = [tuple(scalar(-a, p) for a in u), one]
    else:
        require(t % defect == 0, 'geometric factor uses a divisor of t')
        g = [one if i % defect == 0 else (0,)*D for i in range(t-defect+1)]
        factor = [tuple(scalar(-a, p) for a in one)]+[(0,)*D]*(defect-1)+[one]
    R = len(g)-1
    for k in range(t+1):
        coefficient = tuple(scalar(sum(mul(g[i], factor[k-i])[j]
                                       for i in range(len(g)) if 0 <= k-i < len(factor)), p)
                            for j in range(D))
        expected = tuple(scalar(-a, p) for a in one) if k == 0 else one if k == t else (0,)*D
        require(coefficient == expected, 'ordinary F-polynomial product is exactly W^t-1')

    def multiply_f(vector, f):
        result = [0]*(D*R)
        for k in range(R):
            coefficients = mul(tuple(vector[i*R+k] for i in range(D)), f)
            for i, a in enumerate(coefficients):
                result[i*R+k] = a
        return tuple(result)

    def multiply_w(vector):
        result = [0]*(D*R)
        last = tuple(vector[i*R+R-1] for i in range(D))
        for k in range(R):
            tail = mul(last, g[k])
            for i in range(D):
                result[i*R+k] = scalar((vector[i*R+k-1] if k else 0)-tail[i], p)
        return tuple(result)

    unit = (1,)+(0,)*(D*R-1)
    w_period = unit
    for _ in range(t):
        w_period = multiply_w(w_period)
    require(w_period == unit, 'monic-division arithmetic independently verifies W^t=1')

    @cache
    def residue(a):
        if not any(a):
            return unit
        i = next(i for i, x in enumerate(a) if x)
        vector = residue(tuple(x-(i == j) for j, x in enumerate(a)))
        if i == 0:
            return tuple(vector[(j-1)*R+k] if j else 0 for j in range(D) for k in range(R))
        return multiply_w(multiply_f(vector, u if i == 1 else v))

    inverse_u = tuple(scalar((-1)**i, p) for i in range(D))
    require(mul(u, inverse_u) == one, 'coefficient unit inverse is verified')
    h = mul(v, inverse_u)

    @cache
    def slot(a, c, cheap):
        require(cheap == 0, 'specified nilpotent coefficient coordinate')
        h_power = power(h, c)
        return tuple(h_power[i-a] if i >= a else 0 for i in range(D))

    return {'field': p, 'n': D*R, 'D': D, 'R': R, 'period': t, 'defect': defect,
            'translated_factor': translated, 'g_coefficients': g, 'residue': residue, 'slot': slot}


def profile_deficit_witness(basis, D, t, n):
    require(len(basis) == n and all(0 <= a < D and b >= 0 and c >= 0 and b+c < t for a, b, c in basis),
            'rank-sized monomial set belongs to the full one-period ground')
    R, remainder = divmod(n, D)
    require(not remainder and 0 < R < t, 'rank is D*(t-defect)')
    counts = Counter(b+c for a, b, c in basis)
    require(all(0 <= j < t and count <= D for j, count in counts.items()),
            'actual independent columns obey every level rank cap')
    deficits = [D-counts[j] for j in range(t)]
    missing, center = D*t-n, D*t*(t-1)//2
    degree = sum(b+c for a, b, c in basis)
    require(sum(deficits) == missing, 'the total rank deficit is fixed')
    require(degree == center-sum(j*d for j, d in enumerate(deficits)),
            'degree profile equals the balanced center minus the weighted deficit')
    require(center-missing*(t-1) <= degree <= center, 'profile lies in the thin integer slab')
    volume = n*n*(D-1)*(t-1)*missing*(t-1)
    if R == t-1:
        require(volume <= n**4, 'circuit-quotient profile-box volume is at most n^4')
    return {'level_counts': [counts[j] for j in range(t)], 'level_deficits': deficits,
            'profile_in_unimodular_coordinates': [sum(a for a, b, c in basis), sum(c for a, b, c in basis), degree],
            'degree_profile_interval': [center-missing*(t-1), center],
            'integer_profile_box_volume_bound': volume}


def profile_balance_controls(one_period):
    circuit_records = []
    for record in one_period['models']:
        D, t, n = record['D'], record['periods'][2], record['n']
        witness = profile_deficit_witness(record['standard_set'], D, t, n)
        require(Counter(b+c for a, b, c in record['auxiliary_slots']) == {j: D for j in range(t)},
                'auxiliary slot levels have exactly their full rank')
        circuit_records.append({'periods': record['periods'], 'field': record['field'],
                                'strict_weight': record['strict_weight'], **witness})
    general, exhaustive = [], []
    fixtures = [(D, t, s, p, False) for D, t, s, fields in
                ((2, 3, 1, (0, 2, 3)), (3, 6, 2, (0, 2, 3)),
                 (4, 6, 3, (0, 2, 3)), (5, 10, 2, (0, 2, 5))) for p in fields]
    fixtures.append((3, 3, 1, 3, True))
    for D, t, defect, p, translated in fixtures:
        model = periodic_nilpotent_tower(D, t, defect, p, translated)
        n, scale = model['n'], 2*model['n']+1
        for template in ((t, D, D), (1, D, D+1), (t, D, 2*D)):
            weight = tuple(x*scale**3+scale**i for i, x in enumerate(template))
            require(t*min(weight[1:]) > (D-1)*weight[0], 'general periodic tower is in the one-period region')
            basis = frontier_basis(model['residue'], n, weight, p)
            auxiliary = tower_basis(model, 0, D, t, Q(weight[2]-weight[1], weight[0]))
            require(basis <= auxiliary and len(auxiliary-basis) == D*defect,
                    'general monic-divisor quotient has exactly D*defect omissions')
            witness = profile_deficit_witness(basis, D, t, n)
            general.append({'D': D, 'period': t, 'defect': defect, 'field': p, 'n': n,
                            'translated_factor': translated,
                            'monic_factor_coefficients': [[str(c) for c in row] for row in model['g_coefficients']],
                            'strict_weight': weight, 'standard_set': sorted(basis),
                            'strict_border_relations': border_certificate(basis, model['residue'], weight, p),
                            **witness})
        if D == 2 and t == 3:
            ground = [(a, j-c, c) for a in range(D) for j in range(t) for c in range(j+1)]
            columns = {a: model['residue'](a) for a in ground}
            tested = bases = 0
            degree_min, degree_max = n*t, 0
            for subset in combinations(ground, n):
                pivots = {}
                for a in subset:
                    insert(pivots, columns[a], p)
                tested += 1
                if len(pivots) == n:
                    witness = profile_deficit_witness(subset, D, t, n)
                    degree = witness['profile_in_unimodular_coordinates'][2]
                    degree_min, degree_max = min(degree_min, degree), max(degree_max, degree)
                    bases += 1
            center = D*t*(t-1)//2
            require((degree_min, degree_max) == (center-D*(t-1), center),
                    'both degree-slab endpoints are attained by arbitrary monomial bases')
            exhaustive.append({'D': D, 'period': t, 'defect': defect, 'field': p,
                               'all_subsets_checked': tested, 'independent_bases': bases,
                               'degree_profile_range': [degree_min, degree_max],
                               'scope': 'Every rank-sized subset of the full finite one-period ground; bases need not be monomial downsets.'})
    boundary = []
    for p in (0, 2, 3, 7):
        model, periods = cyclotomic_model((3, 5, 7), p), (3, 5, 7)
        D, t, scale = 8, 7, 2*model['n']+1
        weight = tuple(x*scale**3+scale**i for i, x in enumerate((1000, 1, 1000)))
        basis = frontier_basis(model['residue'], model['n'], weight, p)
        require(t*min(weight[1:]) <= (D-1)*weight[0] and max(b+c for a, b, c in basis) >= t,
                'actual bases outside the region can enter an additional period')
        boundary.append({'periods': periods, 'field': p, 'strict_weight': weight,
                         'standard_set': sorted(basis), 'maximum_level': max(b+c for a, b, c in basis),
                         'strict_border_relations': border_certificate(basis, model['residue'], weight, p)})
    periods, p = (3, 5, 7), 101
    model = cyclotomic_model(periods, p)
    scale = 2*model['n']+1
    weight = tuple(x*scale**3+scale**i for i, x in enumerate((2, 1, 1)))
    basis = frontier_basis(model['residue'], model['n'], weight, p)
    pair_maxima = [max(a+b for a, b, c in basis), max(a+c for a, b, c in basis),
                   max(b+c for a, b, c in basis)]
    require(all(a >= b for a, b in zip(pair_maxima, periods)),
            'an actual positive-price basis can reach all three period thresholds')
    simultaneous = {'periods': periods, 'field': p, 'strict_weight': weight,
                    'standard_set': sorted(basis), 'pair_degree_maxima': pair_maxima,
                    'strict_border_relations': border_certificate(basis, model['residue'], weight, p),
                    'scope': 'One exact nonempty-boundary witness, at bounded price ratios; not an infinite lower family or a separation.'}
    return {'scope': 'Exact rank-cap and profile-slab identities, general nonreduced periodic factors, exhaustive finite-ground bases, and actual outside-region controls. These do not constitute a global fan enumeration or a numerical test of an asymptotic constant.',
            'circuit_profile_witnesses': circuit_records, 'general_periodic_factor_controls': general,
            'exhaustive_ground_controls': exhaustive, 'additional_period_boundary_controls': boundary,
            'simultaneous_period_boundary_control': simultaneous}


def main():
    injections = injection_controls()
    fans = graded_fixtures()
    sharpness = sharpness_controls()
    local = local_refinement_controls()
    towers = tower_relation_controls()
    signed_injections = signed_injection_controls()
    signed_fans = signed_fan_controls()
    zero_opposite = zero_opposite_boundary_controls()
    primary = primary_parameter_controls()
    prefixes = tensor_prefix_controls()
    one_period = one_period_region_controls()
    balance = profile_balance_controls(one_period)
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
    summary.update({'opposite_sign_hilbert_budget_cases': signed_injections['opposite_sign_budget_cases'],
                    'opposite_sign_corner_injections': signed_injections['explicit_corner_injections'],
                    'complete_signed_full_support_fans': len(signed_fans),
                    'zero_opposite_boundary_fans': len(zero_opposite),
                    'signed_and_boundary_comparison_chambers': sum(m['complete_comparison_chambers'] for m in signed_fans+zero_opposite),
                    'signed_and_boundary_states': sum(m['initial_ideal_count'] for m in signed_fans+zero_opposite),
                    'signed_and_boundary_border_relations': sum(len(w['strict_border_relations']) for m in signed_fans+zero_opposite for w in m['witnesses']),
                    'actual_positive_signed_exchanges': sum(len(m['actual_positive_wall_exchanges']) for m in signed_fans+zero_opposite),
                    'primary_parameter_normalizations': len(primary['fixtures'])})
    summary.update({'mixed_period_cost_prefix_models': len(prefixes['prefix_models']),
                    'ordinary_and_graph_prefix_ranks': sum(len(m['prefixes']) for m in prefixes['prefix_models']),
                    'coefficient_subspace_ranks': sum(len(w['coefficient_dimensions']) for m in prefixes['prefix_models'] for w in m['prefixes']),
                    'cycle_union_counterexamples': len(prefixes['cycle_union_counterexamples'])})
    summary.update({'one_period_region_models': len(one_period['models']),
                    'models_outside_older_block_condition': sum(not m['older_complete_block_condition'] for m in one_period['models']),
                    'one_period_strict_border_relations': sum(len(m['strict_border_relations']) for m in one_period['models'])})
    summary.update({'balanced_circuit_profile_controls': len(balance['circuit_profile_witnesses']),
                    'general_periodic_factor_price_controls': len(balance['general_periodic_factor_controls']),
                    'exhaustive_profile_subsets': sum(m['all_subsets_checked'] for m in balance['exhaustive_ground_controls']),
                    'exhaustive_independent_profile_bases': sum(m['independent_bases'] for m in balance['exhaustive_ground_controls']),
                    'additional_period_boundary_controls': len(balance['additional_period_boundary_controls'])})
    report = {'schema_version': 1, 'command': 'python3 -B scripts/run_graded_fan_checks.py',
              'scope': 'Finite controls of uniform grading bounds, local refinements, exact mixed-period prefix ranks, the one-period price restriction and the rank-deficit profile mechanism for its quadratic count. General nonreduced periodic factors and actual additional-period boundaries are included. Complete small graded fans use independent comparisons and actual facets. No full mixed-period fan search, unrestricted G_3 resolution or novelty certificate.',
              'source_sha256': {s: sha256((ROOT/s).read_bytes()).hexdigest() for s in sources},
              'hilbert_budget_controls': injections, 'actual_graded_quotients': fans,
              'pointwise_sharpness_controls': sharpness, 'local_refinement_controls': local,
              'tower_relation_controls': towers,
              'signed_hilbert_budget_controls': signed_injections, 'signed_full_support_controls': signed_fans,
              'zero_opposite_proof_boundary': zero_opposite, 'primary_parameter_controls': primary,
              'mixed_period_prefix_rank_controls': prefixes,
              'one_period_price_controls': one_period,
              'profile_balance_controls': balance,
              'summary': summary}
    (ROOT/'results/research/graded_fan_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(summary, sort_keys=True))


if __name__ == '__main__':
    main()
