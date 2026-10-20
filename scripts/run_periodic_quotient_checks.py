#!/usr/bin/env python3
"""Exact periodic-cover, residual-fan and coherent-arithmetic controls."""
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from hashlib import sha256
from heapq import heappop, heappush
from itertools import combinations, product
from math import comb, gcd, lcm, prod
from pathlib import Path
import json

from run_cyclotomic_price_checks import (border_certificate, circuit_vector,
    cyclotomic_model, frontier_basis, insert, label, scalar, shortest_representatives,
    tensor, univariate_powers)
from run_interleaving_research import minimal_excluded_border
from run_nested_jet_checks import downsets

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


def fixed_cover_trace(periods, fields):
    M, n = prod(periods), prod(m-1 for m in periods)
    scale = 2*n+1
    center = tuple(scale**3+scale**i for i in range(3))
    best = shortest_representatives(periods,center)
    ground = frozenset(a for cost,a in best.values())
    k = M-n
    upper_sizes = {a:sum(all(x<=y for x,y in zip(a,b)) for b in ground) for a in ground}
    shell = frozenset(a for a,count in upper_sizes.items() if count <= k)
    socles = frozenset(a for a in ground if all(tuple(x+(i==j) for j,x in enumerate(a)) not in ground for i in range(3)))
    directions = ((1,1,-1),(1,-1,1),(-1,1,1))
    cube_cache = {}
    def cube_candidates(weight):
        values = [dot(weight,v) for v in directions]
        if any(not x for x in values):
            return ground,[]
        signs = tuple(1 if x>0 else -1 for x in values)
        if signs not in cube_cache:
            steps = [tuple(sign*x for x in v) for sign,v in zip(signs,directions)]
            cubes = []
            for a in sorted(ground):
                corners = [tuple(a[i]+sum(e[j]*steps[j][i] for j in range(3)) for i in range(3))
                           for e in product((0,1),repeat=3)]
                if all(b in ground for b in corners):
                    require(len({label(b,periods) for b in corners}) == 8, 'local cube has eight distinct period labels')
                    cubes.append(corners)
            cube_cache[signs] = (ground-{corners[0] for corners in cubes},cubes)
        return cube_cache[signs]
    beta,gamma = center[1:]
    lower,upper,cover_relations = Q(0),None,[]
    for head in minimal_excluded_border(sorted(ground)):
        tail = best[label(head,periods)][1]
        dx,dy,dz = (a-b for a,b in zip(head,tail))
        constant = beta*dy+gamma*dz
        require(dx*center[0]+constant > 0, 'center lies strictly inside its cover cone')
        if dx > 0:
            lower = max(lower,Q(-constant,dx))
        elif dx < 0:
            bound = Q(-constant,dx)
            upper = bound if upper is None else min(upper,bound)
        else:
            require(constant > 0, 'constant cover inequalities stay strict on the line')
        cover_relations.append({'head':head,'tail':tail})
    require(lower < center[0] and (upper is None or center[0] < upper), 'nonempty full open cover-line interval')
    points = set()
    for a,b in combinations(ground,2):
        dx,dy,dz = (x-y for x,y in zip(a,b))
        if not dx:
            require(beta*dy+gamma*dz, 'parallel costs are distinct throughout the trace')
            continue
        crossing = Q(-beta*dy-gamma*dz,dx)
        if lower < crossing and (upper is None or crossing < upper):
            points.add(crossing)
    endpoints = [lower]+sorted(points)
    if upper is not None:
        endpoints.append(upper)
    intervals = list(zip(endpoints,endpoints[1:]))
    if upper is None:
        intervals.append((endpoints[-1],None))
    results = []
    for p in fields:
        A,C = cyclotomic_model(periods,p),complement_model(periods,p)
        columns = {a:C['column'](tuple((-x)%m for x,m in zip(label(a,periods),periods))) for a in ground}
        previous,states,coverage,last_slope = None,[],[],None
        for left,right in intervals:
            alpha = left+1 if right is None else (left+right)/2
            weight = (alpha.numerator,beta*alpha.denominator,gamma*alpha.denominator)
            omitted = maximum_basis(ground,columns.__getitem__,M-n,weight,p)
            cube_ground,cubes = cube_candidates(weight)
            require(omitted <= cube_ground, 'a strict-minimum cube-circuit column never enters the actual maximum basis')
            require(omitted <= shell, 'actual omissions lie in the coefficient-independent upper k-shell')
            require(maximum_basis(shell & cube_ground,columns.__getitem__,M-n,weight,p) == omitted,
                    'upper-shell and cube pruning preserve the actual dual optimum')
            basis = ground-omitted
            if basis != previous:
                require(frozenset(a for cost,a in shortest_representatives(periods,weight).values()) == ground,
                        'independent shortest-path computation confirms the fixed cover at each new quotient state')
                actual = frontier_basis(A['residue'],n,weight,p)
                require(actual == basis, 'ordinary polynomial frontier independently confirms each trace state')
                border = border_certificate(basis,A['residue'],weight,p)
                slope = sum(a[0] for a in omitted)
                for cube in cubes:
                    require(all(dot(weight,cube[0])<dot(weight,a) for a in cube[1:]), 'the excluded cube corner is uniquely cheapest')
                for cube in cubes[:3]:
                    circuit = [columns[a] for a in cube]
                    require(all(not scalar(sum((-1)**sum(e)*v[i] for e,v in zip(product((0,1),repeat=3),circuit)),p)
                                for i in range(C['n'])), 'alternating cube relation vanishes in ordinary complement coordinates')
                    for deleted in range(8):
                        pivots = {}
                        for i,v in enumerate(circuit):
                            if i != deleted:
                                insert(pivots,v,p)
                        require(len(pivots) == 7, 'every checked cube has every seven columns independent')
                for a in omitted:
                    require(all(b in omitted for b in ground if all(x<=y for x,y in zip(a,b))),
                            'omissions are an upset of the fixed standard downset')
                    require(all(prod(y-x+1 for x,y in zip(a,s)) <= k
                                for s in socles if all(x<=y for x,y in zip(a,s))),
                            'each containing socle supplies a shifted hyperbolic-cross bound')
                require(last_slope is None or slope > last_slope, 'distinct optimum pieces have strictly increasing integral omission slope')
                pair_maxima = [max(a+b for a,b,c in basis),max(a+c for a,b,c in basis),max(b+c for a,b,c in basis)]
                states.append({'first_comparison_interval_left':str(left),'alpha_witness':str(alpha),
                               'strict_weight':weight,'standard_set':sorted(basis),'omission_X_profile':slope,
                               'pair_degree_maxima':pair_maxima,
                               'all_three_thresholds_reached':all(a>=b for a,b in zip(pair_maxima,periods)),
                               'cube_candidates':len(cube_ground),'joint_shell_cube_candidates':len(shell & cube_ground),
                               'strict_cube_exclusions':len(cubes),
                               'strict_border_count':len(border),
                               'strict_border_sha256':sha256(json.dumps(border,sort_keys=True).encode()).hexdigest()})
                previous,last_slope = basis,slope
            coverage.append({'left':str(left),'right':str(right) if right is not None else 'infinity',
                             'state_index':len(states)-1})
        spread = max(a[0] for a in ground)-min(a[0] for a in ground)
        require(len(states) <= 1+k*spread, 'uniform corank integer-slope trace bound')
        results.append({'field':p,'complete_trace':True,'comparison_intervals':coverage,
                        'distinct_states':states,'integer_slope_state_bound':1+k*spread})
    return {'periods':periods,'length':n,'cover_length':M,'cover_standard_set':sorted(ground),
            'fixed_beta_gamma':[beta,gamma],
            'corank':k,'upper_k_shell_size':len(shell),'cover_socle_count':len(socles),
            'open_alpha_interval':[str(lower),str(upper) if upper is not None else 'infinity'],
            'strict_cover_binomial_border':cover_relations,'all_pair_crossings_in_interval':len(points),
            'fields':results,
            'scope':'Complete exact generic comparison-interval traces inside one fixed cover cone, with upper-shell pruning. Other directions of tie refinement, the full two-dimensional refinement fan and asymptotic growth are not counted.'}



def thin_cover_controls():
    records = []
    for q in (5,6,7,11,17):
        r,s,t = 2*q+1,q,q+1
        periods = (r,s,t)
        M,n = prod(periods),prod(m-1 for m in periods)
        scale = 2*n+1
        alpha,beta = t*scale**3+1,s*scale**3+scale
        h = (r*beta*t)//(alpha*s+beta*t)
        require(h==q, 'strict balanced arm threshold has the uniform stated value')
        ground = frozenset((i+s*u,j+t*v,0) for i in range(s) for j in range(t)
                           for u,v in [(u,0) for u in range(h+1)]+[(0,v) for v in range(1,r-h)])
        gamma = 1+max(alpha*a+beta*b for a,b,c in ground)
        weight = (alpha,beta,gamma)
        require(len(ground)==M and {a for cost,a in shortest_representatives(periods,weight).values()}==ground,
                'independent shortest paths confirm the exact infinite-family thin cover formula')
        socles = [a for a in ground if all(tuple(x+(i==j) for j,x in enumerate(a)) not in ground for i in range(3))]
        require(len(socles)==2, 'the two-arm cover has exactly two maximal monomials')
        small_counts = [sum(a+b<r for a,b,c in ground),sum(a+c<s for a,b,c in ground),sum(b+c<t for a,b,c in ground)]
        require(all(count<n for count in small_counts), 'every length-n subset of this cover reaches each period threshold')
        require(M-n==q*(5*q+1), 'thin-family corank is quadratic while every cube candidate survives')
        direct = []
        if q==5:
            for p in (0,2,3):
                A = cyclotomic_model(periods,p)
                basis = frontier_basis(A['residue'],n,weight,p)
                border = border_certificate(basis,A['residue'],weight,p)
                require(basis<=ground and all(max((a+b,a+c,b+c)[i] for a,b,c in basis)>=periods[i] for i in range(3)),
                        'direct quotient arithmetic confirms a simultaneous-threshold state in the thin cover')
                direct.append({'field':p,'standard_set':sorted(basis),'strict_border_count':len(border),
                               'strict_border_sha256':sha256(json.dumps(border,sort_keys=True).encode()).hexdigest()})
        records.append({'parameter':q,'periods':periods,'length':n,'cover_length':M,'corank':M-n,
                        'strict_weight':weight,'cover_socles':sorted(socles),'below_threshold_ground_counts':small_counts,
                        'cube_candidates':M,'candidate_to_corank_ratio':str(Q(M,M-n)),
                        'cover_standard_set_sha256':sha256(json.dumps(sorted(ground)).encode()).hexdigest(),
                        'direct_quotient_controls':direct,
                        'scope':'An infinite-family obstruction to a uniform O(k) cube-only candidate bound, even when all three pair-degree thresholds are reached. Z is absent. The two-socle shell bound survives; the quotient states are planar, not a separating family.'})
    return records



def nonbox_linkage_controls(covers):
    records = []
    for cover in covers:
        periods,weight = tuple(cover['periods']),tuple(cover['strict_weight'])
        ground = frozenset(tuple(a) for a in cover['standard_set'])
        M = len(ground)
        socles = sorted(a for a in ground if all(tuple(x+(i==j) for j,x in enumerate(a)) not in ground for i in range(3)))
        require(socles, 'every finite cover has at least one maximal monomial')
        head = max(ground,key=lambda a:dot(weight,a))
        require(head in socles, 'the norm element has a unique socle head at strict positive prices')
        labels = list(product(*(range(m) for m in periods)))
        index = {g:i for i,g in enumerate(labels[:-1])}
        for p in (0,2):
            @cache
            def residue(a):
                g = label(a,periods)
                return tuple(scalar(int(i==index[g]),p) for i in range(M-1)) if g in index else (scalar(-1,p),)*(M-1)
            # This is the ordinary quotient by the group norm, a one-dimensional ideal.
            basis = frontier_basis(residue,M-1,weight,p)
            require(basis==ground-{head}, 'direct norm-quotient arithmetic confirms the one-dimensional initial annihilator')
            require(all(not scalar(sum(residue(a)[i] for a in ground),p) for i in range(M-1)),
                    'the sum over distinct group labels vanishes, including in ramified characteristic')
            border = border_certificate(basis,residue,weight,p)
            records.append({'periods':periods,'field':p,'strict_weight':weight,'cover_length':M,
                            'cover_socles':socles,'box_cover':len(socles)==1,'norm_initial_monomial':head,
                            'initial_annihilator_dimension':1,'annihilator_of_initial_point_ideal_dimension':len(socles),
                            'dimension_defect':len(socles)-1,'strict_border_count':len(border)})
    return records



def cube_surface_controls(covers, thin):
    """Check the dimension-independent coordinate-boundary charging argument."""
    grounds = [('periodic_cover',tuple(c['periods']),frozenset(tuple(a) for a in c['standard_set'])) for c in covers]
    grounds += [('common_period_box',(q,q,q),frozenset(product(range(q),repeat=3))) for q in (3,5,7,11,17)]
    for c in thin:
        q = c['parameter']
        r,s,t = c['periods']
        ground = frozenset((i+s*u,j+t*v,0) for i in range(s) for j in range(t)
                           for u,v in [(u,0) for u in range(q+1)]+[(0,v) for v in range(1,q+1)])
        grounds.append(('thin_cover',tuple(c['periods']),ground))
    records = []
    for model,periods,ground in grounds:
        require(all(tuple(x-(i==j) for j,x in enumerate(a)) in ground for a in ground for i,x in enumerate(a) if x),
                'the geometric surface control is applied to a downset')
        projections = [len({tuple(x for j,x in enumerate(a) if j!=i) for a in ground}) for i in range(3)]
        for weight in ((3,3,3),(7,1,3),(1,7,3),(1,3,7)):
            directions = ((1,1,-1),(1,-1,1),(-1,1,1))
            steps = [tuple((1 if dot(weight,v)>0 else -1)*x for x in v) for v in directions]
            negative = tuple(sum(max(-v[i],0) for v in steps) for i in range(3))
            positive = tuple(sum(max(v[i],0) for v in steps) for i in range(3))
            low = {a for a in ground if any(a[i]<negative[i] for i in range(3))}
            high = {a for a in ground if tuple(x+y for x,y in zip(a,positive)) not in ground}
            surface = low|high
            candidates = {a for a in ground if not all(tuple(a[i]+sum(e[j]*steps[j][i] for j in range(3)) for i in range(3)) in ground
                                                      for e in product((0,1),repeat=3))}
            require(candidates<=surface, 'cube survivors lie in the explicit bottom/top surface union')
            require(len(low)<=sum(x*y for x,y in zip(negative,projections)), 'bottom slices obey their projection budgets')
            require(len(high)<=sum(x*y for x,y in zip(positive,projections)), 'diagonal top layers obey the telescoping projection budgets')
            require(len(surface)<=3*sum(projections), 'the full sector candidate set has a uniform projection-area bound')
            if model=='common_period_box':
                q = periods[0]
                require(len(candidates)==q**3-(q-3)**3 and len(candidates)<=3*(q**3-(q-1)**3),
                        'box cube count equals the exact formula and is at most three times corank')
            records.append({'model':model,'periods':periods,'ground_size':len(ground),'weight_sector_witness':weight,
                            'projection_sizes':projections,'negative_spans':negative,'positive_spans':positive,
                            'bottom_union_size':len(low),'top_union_size':len(high),'surface_union_size':len(surface),
                            'cube_candidates':len(candidates),'projection_bound':3*sum(projections)})
    return records



def integer_root(value, degree):
    low, high = 0, 1
    while high**degree <= value:
        high *= 2
    while high-low > 1:
        middle = (low+high)//2
        if middle**degree <= value:
            low = middle
        else:
            high = middle
    return low


def dyadic_outer_cover(ground):
    d,M = len(next(iter(ground))),len(ground)
    indices = {tuple((x+1).bit_length()-1 for x in a) for a in ground}
    maximal = sorted(j for j in indices if not any(j!=h and all(x<=y for x,y in zip(j,h)) for h in indices))
    require(len(maximal)<=comb(M.bit_length()-1+d-1,d-1), 'maximal dyadic cells obey the logarithmic antichain bound')
    predecessors,boxes = Counter(),[]
    for j in maximal:
        lower = tuple(2**x-1 for x in j)
        require(lower in ground and 2**sum(j)<=M, 'every occupied dyadic cell has its full anchored predecessor box inside the downset')
        h = tuple(max(x-1,0) for x in j)
        predecessor = tuple(product(*(range(2**x-1,2**(x+1)-1) for x in h)))
        require(set(predecessor)<=ground, 'each half-scale predecessor cell is completely inside the ground')
        predecessors[h] += 1
        sides = tuple(2**(x+1)-1 for x in j)
        require(prod(sides)<=4**d*len(predecessor), 'outer box volume is bounded by its half-scale interior cell')
        boxes.append(sides)
    require(max(predecessors.values())<=2**d, 'clamped predecessor indices have bounded multiplicity')
    require(sum(prod(b) for b in boxes)<=8**d*M, 'the logarithmic outer cover has total volume O_d(M)')
    require(all(any(all(x<y for x,y in zip(a,b)) for b in boxes) for a in ground), 'all original monomials are covered')
    return maximal,boxes


def shell_volume_controls(covers):
    records,aggregates,sharpness = [],[],[]
    def check(ground, ks, keep):
        d,M = len(next(iter(ground))),len(ground)
        maximal,boxes = dyadic_outer_cover(ground)
        # Direct dominance counts and monomial multiplication images are independent constructions.
        counts = {a:sum(all(x<=y for x,y in zip(a,b)) for b in ground) for a in ground}
        require(all(counts[a]==sum(tuple(x+y for x,y in zip(a,b)) in ground for b in ground) for a in ground),
                'principal monomial multiplication rank equals the upper-interval count')
        entries = []
        total_volume = sum(prod(b) for b in boxes)
        for k in ks:
            shell = {a for a,u in counts.items() if u<=k}
            budget,box_rows = 0,[]
            for sides in boxes:
                V = prod(sides)
                steps = tuple(integer_root((k*x**d)//V,d) for x in sides)
                require(prod(x+1 for x in steps)>k, 'adaptive integer box step has more than k divisors')
                require(all(x**d*V<=k*y**d for x,y in zip(steps,sides)), 'all real-root step budgets are verified by integer powers')
                inside = {a for a in ground if all(x<y for x,y in zip(a,sides))}
                boundary = {a for a in inside if tuple(x+y for x,y in zip(a,steps)) not in inside}
                require(shell & inside <= boundary, 'the actual shell in each box lies in its computed exit boundary')
                charge = sum(h*(V//b) for h,b in zip(steps,sides))
                require(len(boundary)<=charge, 'coordinate-fiber exit charging obeys its exact integer budget')
                budget += len(boundary)
                if keep:
                    box_rows.append({'sides':sides,'integer_steps':steps,'exit_boundary_size':len(boundary),'fiber_charge':charge})
            require(len(shell)<=budget, 'the box union covers every low-rank monomial')
            require(len(shell)**d<=d**d*k*len(boxes)*total_volume**(d-1), 'the summed Holder bound is checked without floating roots')
            entries.append({'corank_threshold':k,'shell_size':len(shell),'exit_boundary_sum':budget,'box_checks':box_rows})
        if keep:
            return {'dimension':d,'length':M,'ground':sorted(ground),'maximal_dyadic_indices':maximal,
                    'outer_box_sides':boxes,'total_outer_volume':total_volume,
                    'upper_interval_histogram':dict(sorted(Counter(counts.values()).items())), 'threshold_checks':entries}
        return len(entries),sum(len(boxes) for k in ks)
    for d,maximum in ((2,9),(3,7),(4,5)):
        models,thresholds,box_checks = 0,0,0
        digest = sha256()
        for n in range(1,maximum+1):
            for ground in downsets(n,d):
                t,b = check(ground,range(1,n+1),False)
                models += 1; thresholds += t; box_checks += b
                digest.update(json.dumps(sorted(ground)).encode())
        aggregates.append({'dimension':d,'maximum_length':maximum,'all_downsets':models,
                           'all_corank_threshold_checks':thresholds,'adaptive_box_boundary_checks':box_checks,
                           'grounds_sha256':digest.hexdigest()})
    for cover in covers:
        ground = frozenset(tuple(a) for a in cover['standard_set'])
        M = len(ground)
        k = M-prod(m-1 for m in cover['periods'])
        result = check(ground,sorted({1,k,max(1,M//2),M}),True)
        result.update({'model':'periodic_cover','periods':cover['periods'],'quotient_corank':k})
        records.append(result)
    for d in (2,3,4):
        for L in range(5):
            ground = frozenset(tuple(a) for j in product(range(L+1),repeat=d) if sum(j)<=L
                               for a in product(*(range(2**h-1,2**(h+1)-1) for h in j)))
            M = sum(2**h*comb(h+d-1,d-1) for h in range(L+1))
            H,top = 2**L,2**L*comb(L+d-1,d-1)
            require(len(ground)==M and top<=M<2*top, 'dyadic-cross length and top-cell formulas agree with the enumerated ground')
            result = check(ground,sorted({1,H,M}),True)
            actual = next(c['shell_size'] for c in result['threshold_checks'] if c['corank_threshold']==H)
            require(actual>=top, 'at least half the dyadic-cross algebra consists of monomials of multiplication rank at most 2^L')
            result.update({'model':'dyadic_cross','dyadic_level':L,'top_cell_mass':top})
            records.append(result)
        for L in (8,16,32,64):
            H,top = 2**L,2**L*comb(L+d-1,d-1)
            M = sum(2**h*comb(h+d-1,d-1) for h in range(L+1))
            require(top<=M<2*top, 'the infinite-family exact formulas retain a positive top-cell proportion')
            sharpness.append({'dimension':d,'dyadic_level':L,'length':M,'rank_threshold':H,'guaranteed_shell_size':top,
                              'powered_no_log_ratio':str(Q(top**d,M**(d-1)*H))})
    joint,formulas = [],[]
    directions = ((1,1,-1),(1,-1,1),(-1,1,1))
    for R,L in ((2,1),(2,2),(3,2),(4,3),(5,2)):
        ground = frozenset((a,b,c) for b in range(R*L) for c in range(R*L)
                           for a in range(max(0,R-b//L-c//L)))
        top = frozenset((R-b//L-c//L-1,b,c) for b in range(R*L) for c in range(R*L) if b//L+c//L<R)
        M,k,T = L*L*comb(R+2,3),L*L,L*L*comb(R+1,2)
        require(len(ground)==M and len(top)==T, 'staircase-block length and top-layer formulas')
        require(all(tuple(int(i==j) for j in range(3)) in ground for i in range(3)), 'every obstruction uses all three coordinate generators')
        for a in top:
            upper = sum(all(x<=y for x,y in zip(a,b)) for b in ground)
            images = sum(tuple(x+y for x,y in zip(a,b)) in ground for b in ground)
            require(upper==images and upper<=k, 'every top-layer monomial has multiplication rank at most k')
            require(not all(tuple(a[i]+sum(e[j]*directions[j][i] for j in range(3)) for i in range(3)) in ground
                            for e in product((0,1),repeat=3)), 'every top-layer monomial survives the actual oriented cube test')
        joint.append({'height':R,'block_side':L,'length':M,'rank_threshold':k,'guaranteed_joint_candidates':T,
                      'ground':sorted(ground),'top_layer':sorted(top),
                      'powered_scale_ratio':str(Q(T**3,M*M*k)),
                      'scope':'An arbitrary-downset obstruction using all three variables; it is a monomial quotient, not an asserted periodic group cover or a fan-count lower bound.'})
    for R in (2,4,8,16,32):
        L = R**3
        M,k,T = L*L*comb(R+2,3),L*L,L*L*comb(R+1,2)
        formulas.append({'height':R,'block_side':L,'length':M,'rank_threshold':k,'guaranteed_joint_candidates':T,
                         'powered_scale_ratio':str(Q(T**3,M*M*k))})
    return {'exhaustive_downsets':aggregates,'explicit_controls':records,'infinite_family_integer_formulas':sharpness,
            'joint_cube_shell_obstructions':joint,'joint_cube_shell_integer_formulas':formulas}


def period_gap_controls():
    """Period identities, regional boxes, and an active three-generator obstruction."""
    def sparse_residue(periods, p):
        powers = [univariate_powers(m, p) for m in periods]

        @cache
        def residue(a):
            # Raw exponents and repeated monic division, not a label shortcut.
            u, v, w = a[0]+a[1], a[0]+a[2], a[1]+a[2]
            factors = [[(i, x) for i, x in enumerate(power(e)) if x]
                       for power, e in zip(powers, (u, v, w))]
            return {((i*(periods[1]-1)+j)*(periods[2]-1)+k): scalar(x*y*z, p)
                    for (i,x), (j,y), (k,z) in product(*factors)}
        return residue

    def add_sparse(pivots, original, p):
        v = dict(original)
        while v:
            h = min(v)
            if h not in pivots:
                scale = pow(v[h], -1, p) if p else Q(1)/v[h]
                pivots[h] = {i: scalar(x*scale, p) for i,x in v.items()}
                return True
            factor = v[h]
            for i,x in pivots[h].items():
                value = scalar(v.get(i,0)-factor*x, p)
                if value:
                    v[i] = value
                else:
                    v.pop(i,None)
        return False

    relations = []
    for periods in ((13,7,11),(19,13,17),(25,19,23),(43,37,41),
                    (37,17,29),(35,23,29),(41,17,23),(21,5,11)):
        r,s,t = periods
        require(r>t>s>=3 and all(m%2 for m in periods) and
                all(gcd(a,b)==1 for a,b in combinations(periods,2)), 'odd coprime gap hypotheses')
        h,e = (r-s)//2,(r-t)//2
        g = gcd(h,e)
        h0,e0 = h//g,e//g
        pairs = (((e0*s,0,0),(0,h0*t,0)),
                 ((0,0,h),(s+h,h,0)),((0,0,e),(e,t+e,0)))
        fields = []
        for p in (0,2,3,5):
            powers = [univariate_powers(m,p) for m in periods]
            for a,b in pairs:
                require(label(a,periods)==label(b,periods), 'gap binomial has equal group labels')
                for power,i,j in zip(powers,(a[0]+a[1],a[0]+a[2],a[1]+a[2]),
                                              (b[0]+b[1],b[0]+b[2],b[1]+b[2])):
                    require(power(i)==power(j), 'independent raw monic division verifies every gap identity')
            fields.append(p)
        M,n = prod(periods),prod(m-1 for m in periods)
        scale = 2*M+1
        region_controls = []
        for region in ('Y','X'):
            x,y = (h0*t,e0*s+1) if region=='Y' else (h0*t+1,e0*s)
            z = ((s+h)*x+h*y)//h+1 if region=='Y' else (e*x+(t+e)*y)//e+1
            weight = (x*scale**3+1,y*scale**3+scale,z*scale**3+scale**2)
            if region=='Y':
                margins = (h0*t*weight[1]-e0*s*weight[0],
                           h*weight[2]-(s+h)*weight[0]-h*weight[1])
                sides = (r*s,h0*t,h)
                factor = h*h0
            else:
                margins = (e0*s*weight[0]-h0*t*weight[1],
                           e*weight[2]-e*weight[0]-(t+e)*weight[1])
                sides = (e0*s,r*t,e)
                factor = e*e0
            require(min(margins)>0 and prod(sides)==factor*M, 'strict regional prices and exact common-box volume')
            record = {'region':region,'weight':weight,'strict_margins':margins,
                      'common_box_sides':sides,'box_to_cover_volume_ratio':factor}
            if periods in ((13,7,11),(19,13,17)):
                cover = shortest_representatives(periods,weight)
                ground = [a for cost,a in cover.values()]
                require(all(all(a[i]<sides[i] for i in range(3)) for a in ground), 'actual shortest cover is inside the regional box')
                record['direct_cover_length'] = len(ground)
            region_controls.append(record)
        relations.append({'periods':periods,'length':n,'cover_length':M,'half_gaps':[h,e],
                          'primitive_half_gaps':[h0,e0],'binomials':pairs,
                          'monic_division_fields':fields,'regions':region_controls})

    prefixes, formulas, complete = [],[],[]
    directions = ((1,1,-1),(1,-1,1),(-1,1,1))
    for q in (1,2,3,6,16,32,64):
        periods = r,s,t = 6*q+7,6*q+1,6*q+5
        M,n = prod(periods),prod(m-1 for m in periods)
        scale,H = 2*M+1,2*q+3
        weight = (scale**3+1,scale**3+scale,H*scale**3+scale**2)
        margins = (3*weight[2]-(s+3)*weight[0]-3*weight[1],
                   3*t*weight[1]-s*weight[0])
        require(min(margins)>0 and prod((r*s,3*t,3))==9*M, 'uniform slab family lies strictly in the quadratic regional bound')
        formulas.append({'parameter':q,'periods':periods,'length':n,'cover_length':M,
                         'corank':M-n,'cover_to_corank_ratio':str(Q(M,M-n)),
                         'weight':weight,'strict_relation_margins':margins,
                         'common_box_sides':[r*s,3*t,3]})
        if q not in (1,2,3):
            continue
        for p in (0,2,3,5,7,13):
            residue = sparse_residue(periods,p)
            for target in ((1,0,0),(0,1,0),(0,0,1),(r,0,0),(0,t,0)):
                cost = dot(weight,target)
                earlier = [a for a in product(range(cost//weight[0]+1),
                                             range(cost//weight[1]+1),range(cost//weight[2]+1))
                           if dot(weight,a)<cost]
                pivots = {}
                for a in sorted(earlier,key=lambda a:dot(weight,a)):
                    add_sparse(pivots,residue(a),p)
                rank = len(pivots)
                require(add_sparse(pivots,residue(target),p), 'target is independent of every cheaper monomial by direct monic arithmetic')
                prefixes.append({'parameter':q,'periods':periods,'field':p,'target':target,
                                 'all_cheaper_monomials':len(earlier),'cheaper_span_rank':rank,
                                 'rank_with_target':len(pivots)})
        if q==3:
            continue

        # A complete binary calculation uses an independent tensor-label model.
        def primal_bits(a):
            u,v,w = label(a,periods)
            first = range(r-1) if u==r-1 else (u,)
            second = range(s-1) if v==s-1 else (v,)
            last = (1<<(t-1))-1 if w==t-1 else 1<<w
            result = 0
            for i in first:
                for j in second:
                    result ^= last<<((i*(s-1)+j)*(t-1))
            return result

        def add_bits(pivots, v):
            while v:
                h = v.bit_length()-1
                if h not in pivots:
                    pivots[h] = v
                    return True
                v ^= pivots[h]
            return False

        heap,queued,pivots,basis = [(0,(0,0,0))],{(0,0,0)},{},set()
        while heap and len(basis)<n:
            _,a = heappop(heap)
            if add_bits(pivots,primal_bits(a)):
                basis.add(a)
                for i in range(3):
                    b = tuple(x+(i==j) for j,x in enumerate(a))
                    if b not in queued:
                        queued.add(b)
                        heappush(heap,(dot(weight,b),b))
        require(len(basis)==n, 'complete binary frontier reaches the true quotient dimension')
        best = shortest_representatives(periods,weight)
        ground = frozenset(a for cost,a in best.values())
        require(len({cost for cost,a in best.values()})==M, 'cover costs are distinct')
        require(basis<=ground and max(a[2] for a in ground)<3, 'complete quotient basis is in the three-layer cover')
        oriented = [tuple((1 if dot(weight,v)>0 else -1)*x for x in v) for v in directions]
        corners = [tuple(sum(e[j]*oriented[j][i] for j in range(3)) for i in range(3))
                   for e in product((0,1),repeat=3)]
        require(all(not all(tuple(a[i]+v[i] for i in range(3)) in ground for v in corners)
                    for a in ground), 'every actual oriented cube test fails on the slab cover')
        boundary = [g for g in product(*(range(m) for m in periods)) if 0 in g]
        index = {g:i for i,g in enumerate(boundary)}

        def inverse_dual_bits(a):
            i,j,k = tuple(-x%m for x,m in zip(label(a,periods),periods))
            if 0 in (i,j,k):
                return 1<<index[i,j,k]
            result = 0
            for node in ((i,j,0),(i,0,k),(0,j,k),(i,0,0),(0,j,0),(0,0,k),(0,0,0)):
                result ^= 1<<index[node]
            return result

        dual_pivots,omitted = {},set()
        for a in sorted(ground,key=lambda a:dot(weight,a),reverse=True):
            if add_bits(dual_pivots,inverse_dual_bits(a)):
                omitted.add(a)
                if len(omitted)==M-n:
                    break
        require(len(omitted)==M-n and basis==ground-omitted, 'independent inverse-label maximum basis exactly complements the direct primal frontier')
        require(all(a in basis for a in ((1,0,0),(0,1,0),(0,0,1),(r,0,0),(0,t,0))), 'all three generators and every pair-degree threshold occur in the full basis')
        require(all(tuple(x-(i==j) for j,x in enumerate(a)) in basis
                    for a in basis for i in range(3) if a[i]), 'the complete computed basis is a downset')

        # Independently divide every standard and minimal-border monomial, then
        # solve its normal form. No inferred initial ideal from sampled prices.
        direct = sparse_residue(periods,2)
        ordered = sorted(basis)
        solve_pivots = {}
        for j,a in enumerate(ordered):
            raw = direct(a)
            v = sum(1<<i for i in raw)
            require(v==primal_bits(a), 'ordinary monic division agrees with every standard label column')
            coordinates = 1<<j
            while v:
                h = v.bit_length()-1
                if h not in solve_pivots:
                    solve_pivots[h] = (v,coordinates)
                    break
                pv,pc = solve_pivots[h]
                v ^= pv
                coordinates ^= pc
            require(v!=0, 'normal-form solver independently verifies standard-column independence')
        certificates = []
        for a in sorted(minimal_excluded_border(ordered)):
            v,coordinates = sum(1<<i for i in direct(a)),0
            while v:
                h = v.bit_length()-1
                require(h in solve_pivots, 'standard columns span every border residue')
                pv,pc = solve_pivots[h]
                v ^= pv
                coordinates ^= pc
            tail = [ordered[j] for j in range(n) if (coordinates>>j)&1]
            require(all(dot(weight,b)<dot(weight,a) for b in tail), 'every actual minimal border has a strictly cheaper exact normal form')
            certificates.append({'head':a,'tail':tail})
        complete.append({'parameter':q,'periods':periods,'field':2,'length':n,'cover_length':M,
                         'weight':weight,'cover_sha256':sha256(repr(sorted(ground)).encode()).hexdigest(),
                         'standard_set_sha256':sha256(repr(ordered).encode()).hexdigest(),
                         'maximum_exponents':[max(a[i] for a in basis) for i in range(3)],
                         'maximum_pair_degrees':[max(a[0]+a[1] for a in basis),max(a[0]+a[2] for a in basis),max(a[1]+a[2] for a in basis)],
                         'cube_surviving_candidates':M,'inverse_dual_omissions':sorted(omitted),
                         'strict_border_certificates':certificates})
    return {'gap_relation_controls':relations,'uniform_family_formulas':formulas,
            'all_cheaper_prefix_controls':prefixes,'complete_binary_controls':complete,
            'scope':'The regional count is a classical lattice-volume corollary using explicit period identities. The all-field active-generator obstruction defeats cube-only corank pruning, not the joint shell bound. No infinite prime constellation, full fan count, unrestricted improvement or novelty certificate is asserted.'}


def global_gap_controls():
    """Four global price sectors and a uniform obstruction to a small cover box."""
    sectors = []
    for periods in ((13,7,11),(19,13,17),(25,19,23),(43,37,41),
                    (37,17,29),(35,23,29),(41,17,23),(21,5,11),
                    (13,5,9),(29,17,23)):
        r,s,t = periods
        require(r>t>s>=3 and all(m%2 for m in periods) and
                all(gcd(a,b)==1 for a,b in combinations(periods,2)), 'global gap hypotheses')
        h,e = (r-s)//2,(r-t)//2
        g = gcd(h,e)
        h0,e0,d0 = h//g,e//g,(h-e)//g
        pairs = (((e0*s,0,0),(0,h0*t,0)),
                 ((d0*r,0,0),(0,0,h0*t)),
                 ((0,d0*r,e0*s),(0,0,0)))
        for p in (0,2,3,5):
            powers = [univariate_powers(m,p) for m in periods]
            for a,b in pairs:
                require(label(a,periods)==label(b,periods), 'global pair/product identities have equal labels')
                for power,i,j in zip(powers,(a[0]+a[1],a[0]+a[2],a[1]+a[2]),
                                              (b[0]+b[1],b[0]+b[2],b[1]+b[2])):
                    require(power(i)==power(j), 'raw monic arithmetic checks global identities in every control field')
        M,n = prod(periods),prod(m-1 for m in periods)
        pure = ((r-1)*(s-1),(r-1)*(t-1),(s-1)*(t-1))
        scale = 2*M+1
        for expensive_x_y,expensive_x_z in product((False,True),repeat=2):
            # Independent Y and Z prices realize each strict orientation.
            x = 2*h0*t
            y = e0*s if expensive_x_y else 4*e0*s
            z = d0*r if expensive_x_z else 4*d0*r
            weight = (x*scale**3+1,y*scale**3+scale,z*scale**3+scale**2)
            require((e0*s*weight[0]>h0*t*weight[1])==expensive_x_y and
                    (d0*r*weight[0]>h0*t*weight[2])==expensive_x_z, 'each entire global price sector is nonempty')
            sides = list(pure)
            sides[0 if expensive_x_y else 1] = min(sides[0 if expensive_x_y else 1],
                                                  e0*s if expensive_x_y else h0*t)
            sides[0 if expensive_x_z else 2] = min(sides[0 if expensive_x_z else 2],
                                                  d0*r if expensive_x_z else h0*t)
            cover_sides = [r*s,r*t,s*t]
            cover_sides[0 if expensive_x_y else 1] = min(cover_sides[0 if expensive_x_y else 1],
                                                        e0*s if expensive_x_y else h0*t)
            cover_sides[0 if expensive_x_z else 2] = min(cover_sides[0 if expensive_x_z else 2],
                                                        d0*r if expensive_x_z else h0*t)
            row = {'periods':periods,'normalized_half_gaps':[h0,e0,d0],
                   'x_more_expensive_in_xy_relation':expensive_x_y,
                   'x_more_expensive_in_xz_relation':expensive_x_z,
                   'weight':weight,'quotient_length':n,'pure_axis_degrees':pure,
                   'global_sector_box_sides':sides,'global_sector_box_volume':prod(sides),
                   'dense_cover_sector_box_sides':cover_sides,'dense_cover_sector_box_volume':prod(cover_sides)}
            if periods in ((13,7,11),(19,13,17)):
                E = {a for cost,a in shortest_representatives(periods,weight).values()}
                require(all(all(a[i]<cover_sides[i] for i in range(3)) for a in E), 'complete cover fits its whole-sector box')
                for p in (0,2):
                    model = cyclotomic_model(periods,p)
                    # Independent raw arithmetic on every excluded pure head.
                    for head,tail in pairs[:2]:
                        require(model['residue'](head)==model['residue'](tail), 'global-sector heads have exact quotient identities')
                row['direct_cover_length'] = len(E)
            sectors.append(row)

    obstruction = []
    for q in (1,2,3,6,10,16,32,64):
        periods = r,s,t = 6*q+7,6*q+1,6*q+5
        M = prod(periods)
        scale = 2*M+1
        weight = (t*scale**3+1,2*scale**3+scale,2*scale**3+scale**2)
        targets = ((q,0,0),(0,q*q,0),(0,0,q*q))
        lower_volume = (q+1)*(q*q+1)**2
        row = {'parameter':q,'periods':periods,'cover_length':M,'weight':weight,
               'unperturbed_ray':[t,2,2],'guaranteed_axis_monomials':targets,
               'guaranteed_box_volume':lower_volume,
               'box_to_cover_lower_ratio':str(Q(lower_volume,M))}
        if q in (1,2,3,6):
            best = shortest_representatives(periods,weight)
            E = {a for cost,a in best.values()}
            require(all(best[label(a,periods)][1]==a for a in targets), 'every uniform box-obstruction target is the actual unique shortest representative')
            sides = [1+max(a[i] for a in E) for i in range(3)]
            require(prod(sides)>=lower_volume, 'exact full cover obeys the uniform box-volume lower bound')
            row.update({'actual_cover_box_sides':sides,'actual_cover_box_volume':prod(sides),
                        'actual_box_to_cover_ratio':str(Q(prod(sides),M))})
            # Enumerate all nonnegative representations of cost at most each
            # target's unperturbed score. This checks the uniform proof's
            # low-score alternatives, including every exact tie at X^q.
            score = (t,2,2)
            alternatives = []
            for target in targets:
                budget = dot(score,target)
                candidates = [a for a in product(range(budget//t+1),range(budget//2+1),range(budget//2+1))
                              if dot(score,a)<=budget and label(a,periods)==label(target,periods)]
                expected = ([target] if target[0]==0 else
                            [(q-2*j,2*j,(s+2)*j) for j in range(q//2+1)])
                require(set(candidates)==set(expected), 'complete low-score representation enumeration agrees with the all-q proof')
                require(all(dot(weight,a)>=dot(weight,target) for a in candidates), 'digit perturbation selects the target among every low-score tie')
                alternatives.append({'target':target,'all_low_score_representations':sorted(candidates)})
            row['independent_low_score_controls'] = alternatives
        obstruction.append(row)
    return {'global_sector_controls':sectors,'critical_ray_box_obstruction':obstruction,
            'scope':'A classical common-box vertex bound now covers all positive weights, with O(n^(7/3)) for comparable periods with bounded normalized gaps, not a quadratic count or unrestricted improvement. The uniform box obstruction concerns the full group cover, not an asserted quotient basis or fan-count lower bound. Its infinitesimal perturbations lie at a single positive-ray degeneration; the already proved local linear refinement bound remains valid.'}


def main():
    lattice, dual = periodic_and_dual_controls()
    index_controls = intrinsic_index_controls()
    arithmetic = arithmetic_controls()
    traces = [fixed_cover_trace((3,5,7),(0,2)),fixed_cover_trace((5,7,11),(0,2))]
    require(any(len(f['distinct_states']) > 1 for t in traces for f in t['fields']),
            'a fixed original three-period cover cone has actual quotient refinements')
    nonbox_controls = nonbox_linkage_controls(lattice)
    thin_controls = thin_cover_controls()
    surface_controls = cube_surface_controls(lattice,thin_controls)
    shell_controls = shell_volume_controls(lattice)
    gap_controls = period_gap_controls()
    global_controls = global_gap_controls()
    sources = ('docs/periodic_quotient_arithmetic.md','scripts/run_periodic_quotient_checks.py',
               'scripts/run_cyclotomic_price_checks.py','scripts/run_interleaving_research.py',
               'scripts/run_nested_jet_checks.py')
    report = {'schema_version':1,'command':'python3 -B scripts/run_periodic_quotient_checks.py',
              'scope':'Exact periodic-cover, arithmetic, fixed-cover trace and candidate-filter controls. Direct norm quotients verify the nonbox socle defect. Sector-wide surface budgets, thin-cover obstructions, socle-free dyadic shell bounds and their logarithmic obstruction are checked. Period-gap identities, common regional boxes and active-generator slab obstructions are independently checked by monic division and primal/dual bases. Four global price sectors and a uniform large-box obstruction at a positive arithmetic ray are checked. No unrestricted asymptotic separation, sharp nonbox refinement-fan count or novelty certificate.',
              'source_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              'full_period_lattice_controls':lattice,'inverse_monomial_dual_controls':dual,
              'trivial_intrinsic_grading_controls':index_controls,
              'coherent_arithmetic_witness':arithmetic,
              'fixed_cover_trace_controls':traces,
              'nonbox_linkage_defects':nonbox_controls,
              'cube_surface_controls':surface_controls,
              'thin_cover_obstructions':thin_controls,
              'socle_free_shell_controls':shell_controls,
              'period_gap_controls':gap_controls,
              'global_gap_controls':global_controls,
              'summary':{'full_period_cover_prices':len(lattice),'dual_quotient_price_field_controls':len(dual),
                         'nonbox_reflection_obstructions':sum(x['nonbox_reflection_has_no_constant'] for x in dual),
                         'all_field_pair_minimal_polynomial_controls':len(index_controls),
                         'coherent_bad_prime_witnesses':1,'strict_arithmetic_border_checks':sum(arithmetic['strict_border_counts'].values()),
                         'direct_monomial_pullback_controls':len(arithmetic['pullback_controls'])}}
    report['summary'].update({'complete_fixed_cover_traces':sum(len(t['fields']) for t in traces),
                              'fixed_cover_comparison_intervals':sum(len(f['comparison_intervals']) for t in traces for f in t['fields']),
                              'fixed_cover_distinct_state_checks':sum(len(f['distinct_states']) for t in traces for f in t['fields']),
                              'trace_states_reaching_all_periods':sum(s['all_three_thresholds_reached'] for t in traces for f in t['fields'] for s in f['distinct_states']),
                              'strict_trace_border_checks':sum(s['strict_border_count'] for t in traces for f in t['fields'] for s in f['distinct_states']),
                              'upper_shell_pruned_cover_monomials':sum(t['cover_length']-t['upper_k_shell_size'] for t in traces)})
    report['summary'].update({'linkage_defect_field_controls':len(nonbox_controls),
                              'positive_nonbox_defects':sum(c['dimension_defect']>0 for c in nonbox_controls),
                              'cube_surface_sector_controls':len(surface_controls),
                              'box_cube_formula_controls':sum(c['model']=='common_period_box' for c in surface_controls),
                              'thin_cover_formula_controls':len(thin_controls),
                              'thin_simultaneous_threshold_quotient_controls':sum(len(c['direct_quotient_controls']) for c in thin_controls)})
    report['summary'].update({'socle_free_exhaustive_downsets':sum(c['all_downsets'] for c in shell_controls['exhaustive_downsets']),
                              'socle_free_exhaustive_thresholds':sum(c['all_corank_threshold_checks'] for c in shell_controls['exhaustive_downsets']),
                              'socle_free_adaptive_box_checks':sum(c['adaptive_box_boundary_checks'] for c in shell_controls['exhaustive_downsets']),
                              'socle_free_explicit_fixtures':len(shell_controls['explicit_controls']),
                              'dyadic_sharpness_formula_controls':len(shell_controls['infinite_family_integer_formulas']),
                              'joint_shell_cube_downset_controls':len(shell_controls['joint_cube_shell_obstructions']),
                              'joint_shell_cube_formula_controls':len(shell_controls['joint_cube_shell_integer_formulas'])})
    report['summary'].update({'period_gap_identity_models':len(gap_controls['gap_relation_controls']),
                              'period_gap_regional_boxes':sum(len(c['regions']) for c in gap_controls['gap_relation_controls']),
                              'active_slab_family_formulas':len(gap_controls['uniform_family_formulas']),
                              'all_cheaper_active_prefix_checks':len(gap_controls['all_cheaper_prefix_controls']),
                              'complete_active_slab_bases':len(gap_controls['complete_binary_controls']),
                              'strict_active_slab_borders':sum(len(c['strict_border_certificates']) for c in gap_controls['complete_binary_controls'])})
    report['summary'].update({'global_gap_price_sectors':len(global_controls['global_sector_controls']),
                              'global_gap_direct_covers':sum('direct_cover_length' in c for c in global_controls['global_sector_controls']),
                              'critical_ray_box_formulas':len(global_controls['critical_ray_box_obstruction']),
                              'critical_ray_direct_covers':sum('actual_cover_box_sides' in c for c in global_controls['critical_ray_box_obstruction']),
                              'complete_low_score_target_checks':sum(len(c.get('independent_low_score_controls',[])) for c in global_controls['critical_ray_box_obstruction'])})
    (ROOT/'results/research/periodic_quotient_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['summary'],sort_keys=True))


if __name__ == '__main__':
    main()
