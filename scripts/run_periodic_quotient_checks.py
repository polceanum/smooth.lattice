#!/usr/bin/env python3
"""Exact periodic-cover, residual-fan and coherent-arithmetic controls."""
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from hashlib import sha256
from itertools import combinations, product
from math import lcm, prod
from pathlib import Path
import json

from run_cyclotomic_price_checks import (border_certificate, circuit_vector,
    cyclotomic_model, frontier_basis, insert, label, scalar, shortest_representatives,
    tensor)
from run_interleaving_research import minimal_excluded_border

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
    sources = ('docs/periodic_quotient_arithmetic.md','scripts/run_periodic_quotient_checks.py',
               'scripts/run_cyclotomic_price_checks.py','scripts/run_interleaving_research.py')
    report = {'schema_version':1,'command':'python3 -B scripts/run_periodic_quotient_checks.py',
              'scope':'Exact periodic-cover, arithmetic, fixed-cover trace and candidate-filter controls. Direct norm quotients verify the nonbox socle defect. Sector-wide surface budgets, common-period box cube counts and infinite-family thin-cover formulas are checked. No unrestricted asymptotic separation, full nonbox refinement-fan count or novelty certificate.',
              'source_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              'full_period_lattice_controls':lattice,'inverse_monomial_dual_controls':dual,
              'trivial_intrinsic_grading_controls':index_controls,
              'coherent_arithmetic_witness':arithmetic,
              'fixed_cover_trace_controls':traces,
              'nonbox_linkage_defects':nonbox_controls,
              'cube_surface_controls':surface_controls,
              'thin_cover_obstructions':thin_controls,
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
    (ROOT/'results/research/periodic_quotient_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['summary'],sort_keys=True))


if __name__ == '__main__':
    main()
