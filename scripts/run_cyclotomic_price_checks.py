#!/usr/bin/env python3
"""Exact regional controls; no intermediate-region or unrestricted fan claim."""
from fractions import Fraction as Q
from functools import cache
from hashlib import sha256
from heapq import heappop, heappush
from itertools import product
from math import floor, prod
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


def require(test, message):
    if not test:
        raise AssertionError(message)


def scalar(value, p):
    return value % p if p else Q(value)


def inverse(value, p):
    return pow(value, -1, p) if p else Q(1)/value


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def sparse(vector, p):
    return {i: scalar(v, p) for i, v in enumerate(vector) if scalar(v, p)}


def subtract(target, source, factor, p):
    for i, value in source.items():
        z = scalar(target.get(i, 0)-factor*value, p)
        if z:
            target[i] = z
        else:
            target.pop(i, None)


def insert(pivots, vector, p):
    vector = sparse(vector, p)
    while vector:
        i = min(vector)
        if i in pivots:
            subtract(vector, pivots[i], vector[i], p)
        else:
            f = inverse(vector[i], p)
            pivots[i] = {j: scalar(f*x, p) for j, x in vector.items()}
            return True
    return False


def solver(columns, p):
    pivots = {}
    for j, column in enumerate(columns):
        vector, coordinates = sparse(column, p), {j: scalar(1, p)}
        while vector:
            i = min(vector)
            if i in pivots:
                pv, pc = pivots[i]
                f = vector[i]
                subtract(vector, pv, f, p)
                subtract(coordinates, pc, f, p)
            else:
                f = inverse(vector[i], p)
                pivots[i] = ({k: scalar(f*x, p) for k, x in vector.items()},
                             {k: scalar(f*x, p) for k, x in coordinates.items()})
                break
        require(vector, 'the proposed standard columns are independent')

    def solve(column):
        vector, coordinates = sparse(column, p), {}
        while vector:
            i = min(vector)
            require(i in pivots, 'standard columns span the target residue')
            pv, pc = pivots[i]
            f = vector[i]
            subtract(vector, pv, f, p)
            subtract(coordinates, pc, -f, p)
        return coordinates
    return solve


def frontier_basis(residue, dimension, weight, p):
    """Ordinary quotient arithmetic, without a tower/label order assumption."""
    heap, queued, pivots, basis = [(0, (0, 0, 0))], {(0, 0, 0)}, {}, []
    while heap and len(basis) < dimension:
        _, a = heappop(heap)
        if insert(pivots, residue(a), p):
            basis.append(a)
            for i in range(3):
                b = tuple(x+(i == j) for j, x in enumerate(a))
                if b not in queued:
                    queued.add(b)
                    heappush(heap, (dot(weight, b), b))
    require(len(basis) == dimension, 'frontier reaches the full quotient rank')
    return frozenset(basis)


def border_certificate(basis, residue, weight, p):
    ordered = sorted(basis)
    solve = solver([residue(a) for a in ordered], p)
    outside = {tuple(x+(i == j) for j, x in enumerate(a))
               for a in basis for i in range(3)}-basis
    heads = sorted(a for a in outside if all(
        tuple(x-(i == j) for j, x in enumerate(a)) in basis
        for i in range(3) if a[i]))
    require(all(tuple(x-(i == j) for j, x in enumerate(a)) in basis
                for a in basis for i in range(3) if a[i]), 'monomial downset')
    certificates = []
    for a in heads:
        coordinates = solve(residue(a))
        require(all(dot(weight, a) > dot(weight, ordered[j]) for j in coordinates),
                'every minimal border head is strictly more expensive than its tail')
        certificates.append({'head': a, 'tail': [
            {'exponent': ordered[j], 'coefficient': str(c)}
            for j, c in sorted(coordinates.items())]})
    return certificates


def univariate_powers(m, p):
    """Repeated monic division by S_m, without a period shortcut."""
    powers = [tuple([1]+[0]*(m-2))]

    def power(k):
        while len(powers) <= k:
            a = powers[-1]
            powers.append(tuple(scalar((a[i-1] if i else 0)-a[-1], p)
                                for i in range(m-1)))
        return powers[k]
    return power


def circuit_vector(m, exponent, p):
    exponent %= m
    return tuple(scalar(-1 if exponent == m-1 else int(i == exponent), p)
                 for i in range(m-1))


def tensor(*vectors):
    return tuple(prod(values) for values in product(*vectors))


def cyclotomic_model(periods, p):
    r, s, t = periods
    powers = [univariate_powers(m, p) for m in periods]

    @cache
    def residue(a):
        i, j, k = a
        return tensor(powers[0](i+j), powers[1](i+k), powers[2](j+k))

    def slot(a, c, cheap):
        pair = ((r, s), (r, t), (s, t))[cheap]
        return tensor(circuit_vector(pair[0], a-c, p),
                      circuit_vector(pair[1], a+c, p))
    return {'periods': periods, 'field': p, 'n': prod(m-1 for m in periods),
            'residue': residue, 'slot': slot}


def generic_tower(p, coupled):
    """F=K[x]/x^3; u=1+x, v=1+x^2; T^3=xT+1 or T^3=0."""
    D = R = 3

    def multiply_x(vector):
        return tuple(scalar(vector[(i-1)*R+k] if i else 0, p)
                     for i in range(D) for k in range(R))

    def multiply_t(vector):
        out = [0]*(D*R)
        for i in range(D):
            for k in range(R):
                value = vector[i*R+k]
                if k < R-1:
                    out[i*R+k+1] += value
                elif coupled:
                    out[i*R] += value
                    if i < D-1:
                        out[(i+1)*R+1] += value
        return tuple(scalar(x, p) for x in out)

    def multiply_generator(vector, index):
        if not index:
            return multiply_x(vector)
        shifted = multiply_x(vector)
        if index == 2:
            shifted = multiply_x(shifted)
        return multiply_t(tuple(scalar(a+b, p) for a, b in zip(vector, shifted)))

    @cache
    def residue(a):
        if not any(a):
            return tuple([1]+[0]*(D*R-1))
        i = next(i for i, x in enumerate(a) if x)
        predecessor = tuple(x-(i == j) for j, x in enumerate(a))
        return multiply_generator(residue(predecessor), i)

    def mul_f(a, b):
        return tuple(scalar(sum(a[j]*b[i-j] for j in range(i+1)), p)
                     for i in range(D))

    # (1+x^2)/(1+x) = 1-x+2x^2 mod x^3.
    h = tuple(scalar(x, p) for x in (1, -1, 2))
    h_powers = [(1, 0, 0)]
    for _ in range(R-1):
        h_powers.append(mul_f(h_powers[-1], h))

    def slot(a, c, cheap):
        require(cheap == 0, 'generic model uses its specified cheap coordinate')
        return tuple(h_powers[c][i-a] if i >= a else 0 for i in range(D))
    return {'field': p, 'n': D*R, 'D': D, 'R': R, 'residue': residue,
            'slot': slot, 'relation': 'T^3=x*T+1' if coupled else 'T^3=0'}


def sharp_model(D, p):
    coefficients = [1]
    for i in range(D):
        root = 2**i
        new = [0]*(len(coefficients)+1)
        for j, value in enumerate(coefficients):
            new[j] -= root*value
            new[j+1] += value
        coefficients = new
    powers = [tuple([1]+[0]*(D-1))]

    def power(k):
        while len(powers) <= k:
            previous = powers[-1]
            powers.append(tuple(scalar((previous[i-1] if i else 0)
                                       -previous[-1]*coefficients[i], p) for i in range(D)))
        return powers[k]

    @cache
    def residue(a):
        i, j, k = a
        if j+k >= 2:
            return (0,)*(2*D)
        vector = power(i+D*k)
        return tuple(vector[h] if b == j+k else 0 for h in range(D) for b in range(2))

    def slot(a, c, cheap):
        require(cheap == 0, 'specified sharpness tower coordinate')
        return power(a+D*c)
    require(len({scalar(2**e, p) for e in range(2*D)}) == 2*D,
            'all Vandermonde nodes remain distinct in the selected field')
    return {'field': p, 'n': 2*D, 'D': D, 'R': 2, 'residue': residue,
            'slot': slot, 'f_coefficients': coefficients}


def breakpoint_samples(D, R):
    if R == 1:
        return [], [Q(0)]
    points = sorted({Q(j, h) for j in range(-(D-1), D) for h in range(1, R)})
    samples = [-Q(D)]+[Q(a.numerator+b.numerator, a.denominator+b.denominator)
                            for a, b in zip(points, points[1:])]+[Q(D)]
    require(len(samples) == len(points)+1, 'one witness per complementary interval')
    require(all(a < rho < b for a, rho, b in zip(points, samples[1:], points[1:])),
            'every finite mediant is strictly inside its interval')
    require(len(samples) <= 2*(D-1)*(R-1)+2, 'uniform linear interval bound')
    return points, samples


def tower_basis(model, cheap, D, R, rho):
    p, result = model['field'], set()
    other = [i for i in range(3) if i != cheap]
    for k in range(R):
        pivots = {}
        slots = sorted(product(range(D), range(k+1)), key=lambda a: (a[0]+rho*a[1], a))
        require(len({a+rho*c for a, c in slots}) == len(slots), 'strict slot ordering')
        for a, c in slots:
            if insert(pivots, model['slot'](a, c, cheap), p):
                exponent = [0, 0, 0]
                exponent[cheap], exponent[other[0]], exponent[other[1]] = a, k-c, c
                result.add(tuple(exponent))
                if len(pivots) == D:
                    break
        require(len(pivots) == D, 'each tower block spans exactly F')
    return frozenset(result)


def check_region(model, cheap, D, R, complete):
    points, samples = breakpoint_samples(D, R)
    if not complete:
        samples = sorted({samples[i] for i in (0, len(samples)//4, len(samples)//2,
                                              3*len(samples)//4, len(samples)-1)})
    unique, coverage = {}, []
    for rho in samples:
        A, B = rho.numerator, rho.denominator
        local = (B, D*B+max(-A, 0), D*B+max(A, 0))
        weight = [0, 0, 0]
        order = [cheap]+[i for i in range(3) if i != cheap]
        for i, j in enumerate(order):
            weight[j] = local[i]
        require(min(local[1:]) > (D-1)*local[0], 'strict tower-region inequality')
        require(max(weight) <= 2*model['n'], 'linear-height integer witnesses')
        basis = tower_basis(model, cheap, D, R, rho)
        direct = frontier_basis(model['residue'], model['n'], weight, model['field'])
        require(basis == direct, 'tower slots equal independent ordinary polynomial greedy')
        key = tuple(sorted(basis))
        if key not in unique:
            unique[key] = {'weight': weight, 'standard_set': key,
                           'strict_border_relations': border_certificate(
                               basis, model['residue'], weight, model['field'])}
        coverage.append({'rho': str(rho), 'weight': weight,
                         'standard_set_sha256': sha256(json.dumps(key).encode()).hexdigest()})
    require(all(tuple(i if j == cheap else 0 for j in range(3)) in next(iter(unique))
                for i in range(D)), 'the full pure-coordinate extent occurs')
    return {'cheap_coordinate': cheap, 'D': D, 'R': R, 'complete_regional_coverage': complete,
            'breakpoints': [str(x) for x in points] if complete else None,
            'coverage': coverage, 'distinct_standard_sets_checked': len(unique),
            'witnesses': list(unique.values())}


def shortest_representatives(periods, weight):
    steps = ((1, 1, 0), (1, 0, 1), (0, 1, 1))
    best, heap = {(0, 0, 0): (0, (0, 0, 0))}, [(0, (0, 0, 0), (0, 0, 0))]
    while heap:
        cost, exponent, node = heappop(heap)
        if best[node] != (cost, exponent):
            continue
        for i, step in enumerate(steps):
            target = tuple((x+y) % m for x, y, m in zip(node, step, periods))
            a = tuple(x+(i == j) for j, x in enumerate(exponent))
            candidate = (cost+weight[i], a)
            if target not in best or candidate < best[target]:
                best[target] = candidate
                heappush(heap, (*candidate, target))
    require(len(best) == prod(periods), 'all coprime-period labels are generated')
    return best


def label(a, periods):
    x, y, z = a
    return ((x+y) % periods[0], (x+z) % periods[1], (y+z) % periods[2])


def bounded_representative(a, periods):
    S = sum(periods)
    v = [Q(x)-Q(S, 2) for x in a]
    coefficients = ((v[0]+v[1])/(2*periods[0]),
                    (v[0]+v[2])/(2*periods[1]),
                    (v[1]+v[2])/(2*periods[2]))
    k = [floor(x+Q(1, 2)) for x in coefficients]
    vectors = ((1, 1, -1), (1, -1, 1), (-1, 1, 1))
    result = tuple(a[j]-sum(periods[i]*k[i]*vectors[i][j] for i in range(3))
                   for j in range(3))
    require(min(result) >= 0 and max(result) <= S, 'constructive common S-box representative')
    require(label(result, periods) == label(a, periods), 'rounding preserves periodic labels')
    return result


def period_controls(periods):
    raw = shortest_representatives(periods, (1, 1, 1))
    bounded = {node: bounded_representative(a, periods) for node, (_, a) in raw.items()}
    degrees = [(periods[i]-1)*(periods[j]-1) for i, j in ((0, 1), (0, 2), (1, 2))]
    n = prod(m-1 for m in periods)
    require(prod(degrees) == n*n, 'exact minimum common-box volume')
    ranks = []
    for p in (0, 2, 3):
        model = cyclotomic_model(periods, p)
        for axis, D in enumerate(degrees):
            pivots = {}
            for k in range(D+1):
                a = tuple(k if axis == i else 0 for i in range(3))
                independent = insert(pivots, model['residue'](a), p)
                require(independent == (k < D), 'ordinary pure-axis minimal polynomial degree')
            ranks.append({'field': p, 'axis': axis, 'degree': D})
        for node, a in bounded.items():
            expected = tensor(*(circuit_vector(m, x, p) for m, x in zip(periods, node)))
            require(model['residue'](a) == expected, 'ordinary monic reduction equals label column')
    central = []
    model = cyclotomic_model(periods, 2)
    L = 2*n+1
    for template in ((1, 1, 1), (5, 3, 4), (4, 5, 3)):
        weight = tuple(x*L**3+L**i for i, x in enumerate(template))
        best = shortest_representatives(periods, weight)
        bound = Q(sum(periods)*sum(weight), min(weight))
        require(all(max(a) <= bound for _, a in best.values()), 'least-representative cost bound')
        basis = frontier_basis(model['residue'], n, weight, 2)
        require(all(best[label(a, periods)][1] == a for a in basis),
                'ordinary standard monomials are least label representatives')
        require(max(max(a) for a in basis) <= bound, 'central common-box inclusion')
        central.append({'weight': weight, 'coordinate_bound': str(bound),
                        'standard_set': sorted(basis)})
    return {'periods': periods, 'n': n, 'axis_degrees': degrees,
            'minimum_common_box_size': n*n, 'raw_labels_checked': len(raw),
            'bounded_representatives': [{'label': node, 'exponent': a} for node, a in bounded.items()],
            'ordinary_axis_ranks': ranks, 'central_inclusion_controls': central}


def main():
    models = []
    for periods, fields, complete in (((2, 3, 5), (0, 2, 3), True),
                                     ((3, 5, 7), (0, 2, 3, 5, 7), True),
                                     ((3, 5, 11), (0, 3), False),
                                     ((5, 7, 9), (2, 3), False)):
        for p in fields:
            model = cyclotomic_model(periods, p)
            r, s, t = periods
            dimensions = (((r-1)*(s-1), t-1), ((r-1)*(t-1), s-1), ((s-1)*(t-1), r-1))
            regions = [check_region(model, axis, D, R, complete)
                       for axis, (D, R) in enumerate(dimensions)]
            models.append({'periods': periods, 'field': p, 'n': model['n'], 'regions': regions})
            print(json.dumps({'cyclotomic_model_passed': periods, 'field': p}), flush=True)
    generic = []
    for p, coupled in product((0, 2, 3), (False, True)):
        model = generic_tower(p, coupled)
        generic.append({'field': p, 'F': 'K[x]/x^3', 'u': '1+x', 'v': '1+x^2',
                        'g': model['relation'], 'n': model['n'],
                        'region': check_region(model, 0, model['D'], model['R'], True)})
    sharpness = []
    for D, p in product((2, 3, 5, 8), (0, 101)):
        model = sharp_model(D, p)
        region = check_region(model, 0, D, 2, True)
        require(region['distinct_standard_sets_checked'] == D+1,
                'sharpness family has exactly D+1 regional standard sets')
        # Recover the predicted j from its weight and compare every monomial,
        # independently of the rank computations used to choose the slots.
        for witness in region['witnesses']:
            w = witness['weight']
            rho = Q(w[2]-w[1], w[0])
            costs = [j*j-D*j+rho*j for j in range(D+1)]
            j = costs.index(min(costs))
            expected = {(a, 0, 0) for a in range(D)}
            expected |= {(a, 1, 0) for a in range(D-j)}
            expected |= {(a, 0, 1) for a in range(j)}
            require(set(map(tuple, witness['standard_set'])) == expected,
                    'closed Vandermonde/quadratic-cost formula matches the entire basis')
        outside = []
        for weight in ((2*D+1, 1, 2), (2*D+1, 2, 1),
                       (3, 1, 3*D+1), (3, 3*D+1, 1)):
            require(min(weight[1:]) <= (D-1)*weight[0], 'control lies outside the tower region')
            rho = Q(weight[2]-weight[1], weight[0])
            costs = [j*j-D*j+rho*j for j in range(D+1)]
            require(costs.count(min(costs)) == 1, 'strict closed-form sharpness state')
            j = costs.index(min(costs))
            expected = frozenset({(a, 0, 0) for a in range(D)}
                                 |{(a, 1, 0) for a in range(D-j)}
                                 |{(a, 0, 1) for a in range(j)})
            require(frontier_basis(model['residue'], 2*D, weight, p) == expected,
                    'whole-fan formula survives costs outside the tower condition')
            outside.append({'weight': weight, 'standard_set': sorted(expected),
                            'strict_border_relations': border_certificate(expected, model['residue'], weight, p)})
        sharpness.append({'D': D, 'n': 2*D, 'field': p,
                          'f_coefficients': model['f_coefficients'],
                          'exact_breakpoints': list(range(1-D, D, 2)), 'region': region,
                          'outside_region_formula_controls': outside})
    periods = [period_controls(p) for p in ((2, 3, 5), (3, 5, 7), (3, 5, 11), (5, 7, 9))]
    regions = ([region for model in models for region in model['regions']]
               +[m['region'] for m in generic]+[m['region'] for m in sharpness])
    summary = {'cyclotomic_models': len(models), 'nonreduced_general_towers': len(generic),
               'explicit_linear_sharpness_models': len(sharpness),
               'sharpness_checks_outside_tower_region': sum(len(m['outside_region_formula_controls']) for m in sharpness),
               'complete_bounded_regions': sum(r['complete_regional_coverage'] for r in regions),
               'direct_polynomial_greedy_comparisons': sum(len(r['coverage']) for r in regions)+sum(len(m['outside_region_formula_controls']) for m in sharpness),
               'strict_border_certificates': (sum(len(w['strict_border_relations']) for r in regions for w in r['witnesses'])
                                              +sum(len(w['strict_border_relations']) for m in sharpness for w in m['outside_region_formula_controls'])),
               'distinct_regional_standard_sets_checked': sum(r['distinct_standard_sets_checked'] for r in regions),
               'period_box_models': len(periods), 'bounded_label_representatives': sum(m['raw_labels_checked'] for m in periods),
               'ordinary_pure_axis_rank_controls': sum(len(m['ordinary_axis_ranks']) for m in periods),
               'central_inclusion_controls': sum(len(m['central_inclusion_controls']) for m in periods)}
    sources = ('scripts/run_cyclotomic_price_checks.py', 'docs/cyclotomic_price_regions.md')
    report = {'schema_version': 1, 'command': 'python3 -B scripts/run_cyclotomic_price_checks.py',
              'scope': 'Exact finite controls for the uniform tower, minimum-box and central-region proofs. Complete coverage is regional and limited to the specified fixtures. No intermediate-region fan count, unrestricted bound, separation or priority certificate.',
              'source_sha256': {s: sha256((ROOT/s).read_bytes()).hexdigest() for s in sources},
              'cyclotomic_controls': models, 'nonreduced_tower_controls': generic,
              'linear_sharpness_controls': sharpness,
              'period_and_box_controls': periods, 'summary': summary}
    (ROOT/'results/research/cyclotomic_price_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(summary, sort_keys=True))


if __name__ == '__main__':
    main()
