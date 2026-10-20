#!/usr/bin/env python3
"""Exact controls for periodic lattice covers, dual quotients and bad primes."""
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from hashlib import sha256
from itertools import product
from math import lcm, prod
from pathlib import Path
import json

from run_cyclotomic_price_checks import (border_certificate, circuit_vector,
    cyclotomic_model, frontier_basis, insert, label, scalar, shortest_representatives,
    tensor)

ROOT = Path(__file__).resolve().parents[1]


def require(test, message):
    if not test:
        raise AssertionError(message)


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def bareiss(columns):
    a = [[int(x) for x in row] for row in columns]
    require(all(Q(x).denominator == 1 for row in columns for x in row), 'integer determinant input')
    previous = sign = 1
    for k in range(len(a)-1):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        current = a[k][k]
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                value = a[i][j]*current-a[i][k]*a[k][j]
                require(value % previous == 0, 'fraction-free determinant divisions are exact')
                a[i][j] = value//previous
            a[i][k] = 0
        previous = current
    return sign*a[-1][-1]


def complement_model(periods, p):
    boundary = [g for g in product(*(range(m) for m in periods)) if 0 in g]
    index = {g: i for i, g in enumerate(boundary)}

    @cache
    def column(g):
        g = tuple(x % m for x, m in zip(g, periods))
        if g in index:
            return tuple(int(i == index[g]) for i in range(len(boundary)))
        i, j, k = g
        result = [0]*len(boundary)
        for node, coefficient in (((i,j,0),1), ((i,0,k),1), ((0,j,k),1),
                                  ((i,0,0),-1), ((0,j,0),-1), ((0,0,k),-1), ((0,0,0),1)):
            result[index[node]] += coefficient
        return tuple(scalar(x, p) for x in result)

    def residue(a):
        return column(label(a, periods))
    return {'n': len(boundary), 'field': p, 'column': column, 'residue': residue,
            'boundary_basis': boundary}


def maximum_basis(ground, column, rank, weight, p):
    pivots, selected = {}, set()
    for a in sorted(ground, key=lambda e: dot(weight, e), reverse=True):
        if insert(pivots, column(a), p):
            selected.add(a)
            if len(pivots) == rank:
                break
    require(len(selected) == rank, 'dual optimization reaches its full rank')
    return frozenset(selected)


def pair_minimal_polynomial(r, s):
    numerator = [0]*(r*s+2)
    numerator[0], numerator[1], numerator[r*s], numerator[r*s+1] = 1,-1,-1,1
    divisor = [0]*(r+s+1)
    divisor[0], divisor[r], divisor[s], divisor[r+s] = 1,-1,-1,1
    quotient = [0]*((r-1)*(s-1)+1)
    for i in range(len(quotient)-1,-1,-1):
        coefficient = numerator[i+r+s]
        quotient[i] = coefficient
        for j,value in enumerate(divisor):
            numerator[i+j] -= coefficient*value
    require(not any(numerator) and quotient[0] == quotient[-1] == 1 and quotient[1] == -1,
            'coprime pair polynomial has exact integral division and constant/linear coefficients 1,-1')
    return quotient


def intrinsic_index_controls():
    records = []
    for periods in ((2,3,5),(3,4,5),(3,5,7),(5,7,11)):
        for p in (0,2,3):
            model = cyclotomic_model(periods,p)
            for axis,(i,j) in enumerate(((0,1),(0,2),(1,2))):
                polynomial = pair_minimal_polynomial(periods[i],periods[j])
                pivots, total = {}, [0]*model['n']
                for k,coefficient in enumerate(polynomial):
                    exponent = tuple(k if axis == h else 0 for h in range(3))
                    column = model['residue'](exponent)
                    independent = insert(pivots,column,p)
                    require(independent == (k < len(polynomial)-1), 'the exact pair polynomial is minimal in every checked field')
                    for h,x in enumerate(column):
                        total[h] += coefficient*x
                require(all(not scalar(x,p) for x in total), 'pair polynomial vanishes by ordinary monic arithmetic')
                records.append({'periods':periods,'field':p,'axis':axis,'monic_minimal_coefficients':polynomial,
                                'degree':len(polynomial)-1,'constant_coefficient':1,'linear_coefficient':-1})
    return records


def periodic_and_dual_controls():
    lattice, dual = [], []
    carrier = {(0,0,0): 1}
    for factor in ({(0,0,1):1,(1,1,0):-1}, {(0,1,0):1,(1,0,1):-1},
                   {(1,0,0):1,(0,1,1):-1}):
        result = Counter()
        for a, x in carrier.items():
            for b, y in factor.items():
                result[tuple(i+j for i, j in zip(a,b))] += x*y
        carrier = {a: x for a, x in result.items() if x}
    for periods in ((2,3,5), (3,4,5), (3,5,7)):
        M, n = prod(periods), prod(m-1 for m in periods)
        scale, labels = 2*M+1, tuple(product(*(range(m) for m in periods)))
        orders = [lcm(periods[i], periods[j]) for i,j in ((0,1),(0,2),(1,2))]
        require(prod(orders) == M*M, 'pair-product pure orders have product M^2')
        reference = None
        for template in ((1,1,1), (7,1,3), (1,9,2)):
            weight = tuple(x*scale**3+scale**i for i,x in enumerate(template))
            best = shortest_representatives(periods, weight)
            ground = frozenset(a for cost,a in best.values())
            require(len(ground) == M and {label(a,periods) for a in ground} == set(labels),
                    'full-period standard set has one representative of every label')
            require(all(0 <= a[i] < orders[i] for a in ground for i in range(3)), 'pure period bounds')
            require(all(tuple(x-(i==j) for j,x in enumerate(a)) in ground
                        for a in ground for i in range(3) if a[i]), 'full-period representatives form a downset')
            outside = {tuple(x+(i==j) for j,x in enumerate(a)) for a in ground for i in range(3)}-ground
            heads = sorted(a for a in outside if all(tuple(x-(i==j) for j,x in enumerate(a)) in ground
                                                   for i in range(3) if a[i]))
            relations = []
            for a in heads:
                tail = best[label(a,periods)][1]
                require(dot(weight,a) > dot(weight,tail), 'strict binomial border certificate')
                relations.append({'head':a,'tail':tail})
            profile = tuple(sum(a[i] for a in ground) for i in range(3))
            if reference is None:
                reference = profile
            require(label(tuple(x-y for x,y in zip(profile,reference)),periods) == (0,0,0),
                    'profile differences lie in the index-M label lattice')
            lattice.append({'periods':periods,'length':M,'pure_axis_orders':orders,
                            'label_lattice_index':M,'strict_weight':weight,'standard_set':sorted(ground),
                            'strict_binomial_border':relations,'basis_profile':profile,
                            'normalized_profile_box_volume_bound':M**4})
            for p in (0,2,3):
                A, C = cyclotomic_model(periods,p), complement_model(periods,p)
                inverse_columns = [C['column'](tuple((-x)%m for x,m in zip(g,periods))) for g in labels]
                a_columns = [tensor(*(circuit_vector(m,x,p) for m,x in zip(periods,g))) for g in labels]
                require(C['n'] == M-n, 'complement quotient has exact complementary length')
                for i in range(n):
                    for j in range(C['n']):
                        require(not scalar(sum(a[i]*c[j] for a,c in zip(a_columns,inverse_columns)),p),
                                'inverse-monomial complement rows are orthogonal to the circuit rows')
                for g in labels:
                    result = [0]*C['n']
                    for epsilon in product((0,1), repeat=3):
                        node = tuple((x+e)%m for x,e,m in zip(g,epsilon,periods))
                        coefficient = (-1)**(3-sum(epsilon))
                        for i,x in enumerate(C['column'](node)):
                            result[i] += coefficient*x
                    require(all(not scalar(x,p) for x in result), 'every translate of (U-1)(V-1)(W-1) vanishes')
                require(all(not scalar(sum(coefficient*C['residue'](a)[i] for a,coefficient in carrier.items()),p)
                            for i in range(C['n'])), 'the fixed degree-six carrier vanishes in C')
                basis = frontier_basis(A['residue'],n,weight,p)
                require(basis <= ground, 'quotient standard set is contained in its full-period cover standard set')
                inverse = lambda a: C['column'](tuple((-x)%m for x,m in zip(label(a,periods),periods)))
                omitted = ground-basis
                require(maximum_basis(ground,inverse,C['n'],weight,p) == omitted,
                        'actual omissions are the maximum basis in the known dual matroid')
                rho = tuple(max(a[i] for a in ground) for i in range(3))
                reflected = {tuple(x-y for x,y in zip(rho,a)) for a in ground}
                reflection_basis = {tuple(x-y for x,y in zip(rho,a)) for a in omitted}
                pivots = {}
                for a in reflection_basis:
                    insert(pivots,C['residue'](a),p)
                require(len(pivots) == C['n'], 'reflected omissions form a complement-quotient basis')
                selected = maximum_basis(reflected,C['residue'],C['n'],tuple(-x for x in weight),p)
                require(selected == reflection_basis, 'reflection turns maximum cost into finite-ground minimum cost')
                nonbox = rho not in ground
                if nonbox:
                    require((0,0,0) not in reflected, 'a nonbox reflected ground cannot contain any global standard set')
                dual.append({'periods':periods,'field':p,'length_A':n,'length_C':C['n'],
                             'strict_weight':weight,'standard_set_A':sorted(basis),'omitted_cover_monomials':sorted(omitted),
                             'reflection_corner':rho,'reflected_dual_basis':sorted(reflection_basis),
                             'nonbox_reflection_has_no_constant':nonbox,
                             'strict_A_border_count':len(border_certificate(basis,A['residue'],weight,p))})
    require(any(x['nonbox_reflection_has_no_constant'] for x in dual), 'the reflection obstruction occurs in actual covers')
    return lattice, dual


def arithmetic_controls():
    periods, n, L = (5,7,11), 240, 481
    weight = tuple(L**3+L**i for i in range(3))
    models = {p:cyclotomic_model(periods,p) for p in (0,2,3)}
    bases = {p:frontier_basis(m['residue'],n,weight,p) for p,m in models.items()}
    dets = {p:bareiss([models[0]['residue'](a) for a in sorted(b)]) for p,b in bases.items()}
    require(abs(dets[0]) == 2 and abs(dets[2]) == 1, 'actual coherent bases have indices two and one')
    require(bases[3] == bases[0], 'the determinant-two basis survives characteristic three')
    require(bases[0]-bases[2] == {(1,11,0)} and bases[2]-bases[0] == {(0,12,0)},
            'characteristic two makes the exact stated one-column change')
    pivots = {}
    for a in bases[0]:
        insert(pivots,models[2]['residue'](a),2)
    require(len(pivots) == n-1, 'the rational basis loses exactly one rank modulo two')
    borders = {p:border_certificate(bases[p],models[p]['residue'],weight,p) for p in models}
    relation = next(f for f in borders[0] if tuple(f['head']) == (0,12,0))
    denominator = lcm(*(Q(t['coefficient']).denominator for t in relation['tail']))
    require(denominator == 2, 'the relevant rational border has denominator exactly two')
    polynomial = {(0,12,0):2}
    for term in relation['tail']:
        polynomial[tuple(term['exponent'])] = -int(denominator*Q(term['coefficient']))
    require(all(not sum(c*models[0]['residue'](a)[i] for a,c in polynomial.items()) for i in range(n)),
            'cleared border relation vanishes by ordinary integral quotient arithmetic')
    leading_two = max((a for a,c in polynomial.items() if c%2), key=lambda a:dot(weight,a))
    require(leading_two == (1,11,0), 'modulo two the integral relation changes its leading monomial')
    require(Counter(map(sum,bases[0])) == Counter(map(sum,bases[2])), 'ordinary degree distributions agree in this witness')
    lifts = []
    for degree in (2,3):
        lifted_weight = (weight[0],degree*weight[1],degree*weight[2])
        for p in (0,2):
            original = models[p]
            @cache
            def residue(a):
                q,remainder = divmod(a[0],degree)
                vector = original['residue']((q,a[1],a[2]))
                return tuple(x if i == remainder else 0 for i in range(degree) for x in vector)
            predicted = frozenset((degree*a+i,b,c) for a,b,c in bases[p] for i in range(degree))
            actual = frontier_basis(residue,degree*n,lifted_weight,p)
            require(actual == predicted, 'independent monomial-pullback frontier equals the lifted original basis')
            if p == 2:
                deficient = {}
                for a,b,c in bases[0]:
                    for i in range(degree):
                        insert(deficient,residue((degree*a+i,b,c)),2)
                require(len(deficient) == (n-1)*degree, 'direct lifted rational columns have the proved modular rank defect')
            for f in borders[p]:
                head = (degree*f['head'][0],*f['head'][1:])
                tails = [((degree*t['exponent'][0],*t['exponent'][1:]),Q(t['coefficient'])) for t in f['tail']]
                require(all(dot(lifted_weight,head)>dot(lifted_weight,a) for a,c in tails), 'every pulled-back border remains strict')
            lifts.append({'degree_in_X':degree,'field':p,'length':degree*n,'strict_weight':lifted_weight,
                          'standard_set_sha256':sha256(json.dumps(sorted(actual)).encode()).hexdigest(),
                          'predicted_integral_index_of_Q_basis':2**degree,
                          'Q_basis_rank_mod_two':(n-1)*degree,
                          'scope':'Direct lifted frontier and pulled-back strict borders; index growth and unchanged full fan follow from the uniform free-module proof.'})
    return {'periods':periods,'length':n,'strict_weight':weight,
            'Q_standard_set':sorted(bases[0]),'F2_standard_set':sorted(bases[2]),
            'integer_determinants':{str(p):d for p,d in dets.items()},'Q_basis_rank_mod_two':239,
            'Q_only':sorted(bases[0]-bases[2]),'F2_only':sorted(bases[2]-bases[0]),
            'cleared_integer_relation':[{'exponent':a,'coefficient':c} for a,c in sorted(polynomial.items())],
            'leading_monomial_Q':(0,12,0),'leading_monomial_F2':leading_two,
            'strict_border_counts':{str(p):len(v) for p,v in borders.items()},
            'strict_border_sha256':{str(p):sha256(json.dumps(v,sort_keys=True).encode()).hexdigest() for p,v in borders.items()},
            'ordinary_degree_histogram':dict(sorted(Counter(map(sum,bases[0])).items())),
            'pullback_controls':lifts}


def main():
    lattice, dual = periodic_and_dual_controls()
    index_controls = intrinsic_index_controls()
    arithmetic = arithmetic_controls()
    sources = ('docs/periodic_quotient_arithmetic.md','scripts/run_periodic_quotient_checks.py',
               'scripts/run_cyclotomic_price_checks.py')
    report = {'schema_version':1,'command':'python3 -B scripts/run_periodic_quotient_checks.py',
              'scope':'Exact controls for the periodic lattice cover, known dual quotient, characteristic-sensitive coherent bases and monomial-pullback obstruction. No full fan enumeration, asymptotic separation or novelty certificate.',
              'source_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              'full_period_lattice_controls':lattice,'inverse_monomial_dual_controls':dual,
              'trivial_intrinsic_grading_controls':index_controls,
              'coherent_arithmetic_witness':arithmetic,
              'summary':{'full_period_cover_prices':len(lattice),'dual_quotient_price_field_controls':len(dual),
                         'nonbox_reflection_obstructions':sum(x['nonbox_reflection_has_no_constant'] for x in dual),
                         'all_field_pair_minimal_polynomial_controls':len(index_controls),
                         'coherent_bad_prime_witnesses':1,'strict_arithmetic_border_checks':sum(arithmetic['strict_border_counts'].values()),
                         'direct_monomial_pullback_controls':len(arithmetic['pullback_controls'])}}
    (ROOT/'results/research/periodic_quotient_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['summary'],sort_keys=True))


if __name__ == '__main__':
    main()
