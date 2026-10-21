#!/usr/bin/env python3
"""Exact controls for growing-pole carrier fans and generic stabilization.

No growth fit, novelty claim, or unrestricted extremal result is inferred.
Generic finite-quotient cones are computed by carrier dependence and
Cramer's rule, rather than by pretending numerical samples are generic.
"""
from functools import cache
from fractions import Fraction
from itertools import combinations
from math import comb, gcd
from pathlib import Path
import hashlib
import json

from run_interleaving_research import minimal_excluded_border
from run_initial_cone_height_checks import clean, leading, normal_form, subtract_shift
from run_extremal_initial_ideals_checks import adjacent_weights
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross, primitive

ROOT = Path(__file__).resolve().parents[1]


def coefficient(x, p):
    return x % p if p else Fraction(x)


def relations(k, p):
    result = [clean({(1, 1, 0): 1, (0, 0, 0): -1}, p)]
    for i in range(k + 1):
        g = {(0, i, 0): -1}
        for j in range(k + 1):
            exponent = j - i
            a = (max(exponent, 0), max(-exponent, 0), 1)
            g[a] = (-1) ** (k-j) * comb(k, j)
        result.append(clean(g, p))
    return result


def sector(k, w):
    alpha, beta, _ = w
    assert k * alpha % (alpha + beta), 'carrier wall excluded in this control'
    m = k * alpha // (alpha + beta)
    return k - m, m + 1


def is_standard(u, a, b):
    x, y, z = u
    return not x * y and (not z or (x < a and y < b))


def cleared_columns(points, k, p):
    """Ordinary polynomial coefficients after one common denominator.

    None denotes the cylinder: Laurent monomials in two independent
    indeterminates, represented by distinct unit columns.
    """
    points = tuple(sorted(set(points)))
    if k is None:
        labels = sorted({(a-b, c) for a, b, c in points})
        return {u: tuple(int((u[0]-u[1], u[2]) == v) for v in labels) for u in points}
    shift = max((b-a for a, b, c in points), default=0)
    shift = max(0, shift)
    level = max((c for a, b, c in points), default=0)
    degree = max((a-b+shift+k*(level-c) for a, b, c in points), default=0)
    columns = {}
    for u in points:
        a, b, c = u
        offset, power = a-b+shift, k*(level-c)
        vector = [0] * (degree+1)
        for j in range(power+1):
            vector[offset+j] = coefficient((-1)**(power-j)*comb(power,j), p)
        columns[u] = tuple(vector)
    return columns


def echelon(vectors, p):
    pivots = {}
    for vector in vectors:
        v = [coefficient(x,p) for x in vector]
        for j, row in sorted(pivots.items()):
            factor = v[j]
            if factor:
                v = [coefficient(x-factor*y,p) for x,y in zip(v,row)]
        j = next((i for i,x in enumerate(v) if x), None)
        if j is not None:
            inverse = pow(v[j],-1,p) if p else 1/v[j]
            pivots[j] = tuple(coefficient(x*inverse,p) for x in v)
    return pivots


def rank(vectors,p):
    return len(echelon(vectors,p))


def solver(columns,p):
    n = len(columns)
    pivots = {}
    for index, column in enumerate(columns):
        vector = [coefficient(x,p) for x in column]
        coords = [coefficient(index==j,p) for j in range(n)]
        for j,(row,combination) in sorted(pivots.items()):
            factor = vector[j]
            if factor:
                vector = [coefficient(x-factor*y,p) for x,y in zip(vector,row)]
                coords = [coefficient(x-factor*y,p) for x,y in zip(coords,combination)]
        j = next((i for i,x in enumerate(vector) if x), None)
        assert j is not None
        inv = pow(vector[j],-1,p) if p else 1/vector[j]
        pivots[j] = (tuple(coefficient(x*inv,p) for x in vector),
                     tuple(coefficient(x*inv,p) for x in coords))
    def solve(column):
        vector = [coefficient(x,p) for x in column]
        coords = [coefficient(0,p)] * n
        for j,(row,combination) in sorted(pivots.items()):
            factor = vector[j]
            if factor:
                vector = [coefficient(x-factor*y,p) for x,y in zip(vector,row)]
                coords = [coefficient(x+factor*y,p) for x,y in zip(coords,combination)]
        return coords,not any(vector)
    return solve


def truncated_greedy(ground,columns,w,n,p):
    selected, pivots = [],{}
    for u in sorted(ground,key=lambda a:(dot(a,w),a)):
        v = [coefficient(x,p) for x in columns[u]]
        for j,row in sorted(pivots.items()):
            factor=v[j]
            if factor:
                v=[coefficient(x-factor*y,p) for x,y in zip(v,row)]
        j=next((i for i,x in enumerate(v) if x),None)
        if j is not None:
            inv=pow(v[j],-1,p) if p else 1/v[j]
            pivots[j]=tuple(coefficient(x*inv,p) for x in v)
            selected.append(u)
            if len(selected)==n:
                return frozenset(selected)
    raise AssertionError('insufficient rank')


@cache
def generic_cone(B,k,p):
    border=tuple(minimal_excluded_border(B))
    columns=cleared_columns(B+border,k,p)
    solve=solver([columns[u] for u in B],p)
    normals=set()
    for u in border:
        coords,in_span=solve(columns[u])
        # Independent carrier functions have every n-subset independent in
        # the generic evaluation truncation. Dependent ones retain exactly
        # their ordinary carrier normal form.
        support=[v for v,c in zip(B,coords) if c] if in_span else B
        for v in support:
            z=tuple(a-b for a,b in zip(u,v))
            factor=gcd(*z)
            normals.add(tuple(x//factor for x in z))
    polygon=[tuple(Fraction(i==j) for i in range(3)) for j in range(3)]
    for z in sorted(normals):
        polygon=clip(polygon,z,1)
        if not polygon:
            break
    return polygon,sorted(normals)


def formula_basis(n,k,w,ground):
    if k is None:
        candidates=ground
    else:
        a,b=sector(k,w)
        candidates=[u for u in ground if is_standard(u,a,b)]
    return frozenset(sorted(candidates,key=lambda u:(dot(u,w),u))[:n])


def full_fan(n,k,p=0):
    ground=tuple(u for u in hyperbolic_cross(n,3) if not u[0]*u[1])
    comparison=sorted({primitive(tuple(x-y for x,y in zip(a,b)))
                       for a,b in combinations(ground,2)
                       if min(x-y for x,y in zip(a,b))<0<max(x-y for x,y in zip(a,b))})
    # Mixed radix avoids monomial ties; avoid a carrier wall as well.
    base=n+2
    while k is not None and k % (base+1)==0:
        base+=1
    pending=[(1,base,base**2)]
    visited={}
    while pending:
        w=pending.pop()
        if k is not None and k*w[0] % (w[0]+w[1])==0:
            # A generic perturbation within the strict comparison chamber.
            scale=10*(n+1)*sum(w)+1
            w=tuple(scale*x+(i==0) for i,x in enumerate(w))
        B=formula_basis(n,k,w,ground)
        if B in visited:
            continue
        poly,normals=generic_cone(tuple(sorted(B)),k,p)
        assert poly and area2(poly)>0 and all(dot(z,w)>0 for z in normals), (n,k,p,w,sorted(B),[z for z in normals if dot(z,w)<=0],poly)
        visited[B]=(poly,w)
        pending.extend(adjacent_weights(poly,normals,comparison))
    assert sum(area2(poly) for poly,w in visited.values())==1
    return visited


def carrier_controls():
    rows=[]; pairs=0; replay=0; ranks=0
    for k in (1,2,3,4,5,7):
        for m in range(k):
            # k*alpha/(alpha+beta)=m+1/2, independently of gamma.
            alpha,beta=2*m+1,2*k-2*m-1
            for gamma in (1,2*k+1):
                w=(alpha,beta,gamma);a,b=sector(k,w)
                assert (a,b)==(k-m,m+1)
                for p in (0,2,3,5):
                    G=relations(k,p)
                    heads=[leading(g,w) for g in G]
                    assert (1,1,0) in heads and (a,0,1) in heads and (0,b,1) in heads
                    for g in G:
                        columns=cleared_columns(g,k,p)
                        assert all(coefficient(sum(c*columns[u][i] for u,c in g.items()),p)==0
                                   for i in range(len(next(iter(columns.values())))))
                        replay+=1
                    for i,j in combinations(range(len(G)),2):
                        u,v=heads[i],heads[j]
                        common=tuple(max(x,y) for x,y in zip(u,v))
                        factor=pow(G[i][u],-1,p) if p else 1/G[i][u]
                        first={tuple(x+y-z for x,y,z in zip(t,common,u)):coefficient(c*factor,p)
                               for t,c in G[i].items()}
                        factor=pow(G[j][v],-1,p) if p else 1/G[j][v]
                        pair=subtract_shift(first,G[j],tuple(x-y for x,y in zip(common,v)),factor,p)
                        assert not normal_form(pair,G,w,p)
                        pairs+=1
                    finite=tuple(u for u in hyperbolic_cross(2*k+3,3) if is_standard(u,a,b))
                    columns=cleared_columns(finite,k,p)
                    assert rank(columns.values(),p)==len(finite)
                    ranks+=1
                    rows.append({'k':k,'sector':m,'weights':w,'characteristic':p,
                                 'heads':[(1,1,0),(a,0,1),(0,b,1)],
                                 'independent_standard_columns':len(finite)})
    return {'cases':len(rows),'ordinary_s_pairs':pairs,'cleared_relation_replays':replay,
            'independent_carrier_rank_checks':ranks,'witnesses':rows}


def finite_controls():
    rows=[]; gaussian=0; evaluation=0; stable=0; same_field=0; exhaustive=0
    for n,k in ((4,1),(4,2),(5,1),(5,2),(6,2),(6,3),(8,2),(8,3),(8,4)):
        visited=full_fan(n,k)
        ground=tuple(u for u in hyperbolic_cross(n,3) if not u[0]*u[1])
        fields={p:cleared_columns(ground,k,p) for p in (0,2,3,101)}
        # Degree-D polynomials inject into D+1 ordinary evaluations, over
        # Q here. This independent representation checks coefficient ranks.
        rational=fields[0]
        size=len(next(iter(rational.values())))
        evaluated={u:tuple(sum(c*Fraction(t)**i for i,c in enumerate(v))
                          for t in range(2,size+2)) for u,v in rational.items()}
        witnesses=[]
        for B,(poly,w) in sorted(visited.items(),key=lambda item:sorted(item[0])):
            assert B==formula_basis(n,k,w,ground)
            for p,columns in fields.items():
                assert truncated_greedy(ground,columns,w,n,p)==B
                gaussian+=1
            assert truncated_greedy(ground,evaluated,w,n,0)==B
            evaluation+=1
            witnesses.append({'basis':sorted(B),'weight':w,
                              'cone_polygon':[[str(x) for x in vertex] for vertex in poly]})
        if n<=5:
            from run_structural_complexity_checks import finite_staircases
            feasible=set()
            for B in finite_staircases(n,3):
                if any(u[0]*u[1] for u in B):
                    continue
                columns=cleared_columns(B,k,0)
                if rank(columns.values(),0)!=n:
                    continue
                poly,_=generic_cone(tuple(sorted(B)),k,0)
                exhaustive+=1
                if poly and area2(poly)>0:
                    feasible.add(B)
            assert feasible==set(visited)
        if n<=6:
            mod2=full_fan(n,k,2)
            assert set(mod2)==set(visited)
            assert all(area2(mod2[B][0])==area2(visited[B][0]) and
                       all(dot(z,v)>=0 for v in mod2[B][0]
                           for z in generic_cone(tuple(sorted(B)),k,0)[1])
                       for B in visited)
            same_field+=len(visited)
        cylinder=None
        if n<=2*k+1:
            cylinder=full_fan(n,None)
            assert set(cylinder)==set(visited)
            assert all(area2(cylinder[B][0])==area2(visited[B][0]) and
                       all(dot(z,v)>=0 for v in cylinder[B][0]
                           for z in generic_cone(tuple(sorted(B)),k,0)[1])
                       for B in visited)
            stable+=len(visited)
        rows.append({'n':n,'k':k,'base_field_for_full_fan':'Q',
                     'generic_initial_ideals':len(visited),
                     'generic_cylinder_equal':cylinder is not None,'witnesses':witnesses})
    return {'models':len(rows),'independent_coefficient_greedy_checks':gaussian,
            'unisolvent_evaluation_greedy_checks':evaluation,
            'complete_characteristic_two_cone_equalities':same_field,
            'complete_cylinder_cone_equalities':stable,
            'exhaustive_staircase_cone_tests':exhaustive,'fixtures':rows}



def special_sample_controls():
    from run_structural_complexity_checks import univariate_product_mod, univariate_power_mod
    rows=[]
    for k in (2,3,5):
        modulus=[Fraction((-1)**(k-j)*comb(k,j)) for j in range(k+1)]
        modulus[0]-=2
        x=(Fraction(0),Fraction(1))+(Fraction(0),)*(k-2)
        y=tuple(-modulus[i+1]/modulus[0] for i in range(k))
        z=(Fraction(1,2),)+(Fraction(0),)*(k-1)
        unit=(Fraction(1),)+(Fraction(0),)*(k-1)
        assert univariate_product_mod(x,y,modulus)==unit
        assert univariate_power_mod((-1,1),k,modulus)==tuple(2*c for c in unit)
        ground=tuple(hyperbolic_cross(k,3))
        columns={}
        for u in ground:
            value=unit
            for image,power in zip((x,y,z),u):
                value=univariate_product_mod(value,univariate_power_mod(image,power,modulus),modulus)
            columns[u]=value
        allocations=[]
        for a in range(k):
            b=k-1-a
            w=(2*b+1,2*a+1,1)
            B=truncated_greedy(ground,columns,w,k,0)
            expected={(j,0,0) for j in range(a+1)}|{(0,j,0) for j in range(1,b+1)}
            assert B==expected
            allocations.append({'basis':sorted(B),'weight':w})
        cylinder=formula_basis(k,None,(k+1,(k+1)**2,1),ground)
        assert cylinder==frozenset((0,0,j) for j in range(k))
        assert truncated_greedy(ground,columns,(k+1,(k+1)**2,1),k,0)!=cylinder
        rows.append({'n':k,'k':k,'field':'Q','special_modulus_low_to_high':[str(c) for c in modulus],
                     'constant_Z':'1/2','actual_initial_ideals_by_uniform_proof':k,
                     'generic_cylinder_lower_bound':k+1,'allocations':allocations})
    return {'counterexamples':len(rows),'scope':'Explicit reduced special samples refute extending generic stabilization to every f.','fixtures':rows}


def main():
    if not __debug__:
        raise RuntimeError('Run without -O: assertions are required.')
    carrier=carrier_controls()
    finite=finite_controls()
    special=special_sample_controls()
    sources=('docs/pole_order_fans.md','scripts/run_pole_order_checks.py',
             'docs/extremal_initial_ideals.md','docs/rational_shadow_bound.md',
             'scripts/run_initial_cone_height_checks.py','scripts/run_extremal_initial_ideals_checks.py',
             'scripts/run_structural_complexity_checks.py','scripts/run_interleaving_research.py')
    report={'schema_version':1,'status':'passed',
            'command':'python3 -B scripts/run_pole_order_checks.py',
            'uniform_proof':'docs/pole_order_fans.md',
            'scope':'Actual generic n-point quotients of XY=1,(X-1)^k*Z=1 over every characteristic. Numerical special samples are not assumed generic. Laurent identities, generic evaluation and the earlier threshold/shadow tools are classical inputs.',
            'carrier':carrier,'generic_finite_quotients':finite,'special_sample_counterexamples':special,
            'unrestricted_G3_improved':False,'polynomial_lower_family_established':False,
            'novelty_established':False,'goal_complete':False,
            'source_sha256':{s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}}
    out=ROOT/'results/research/pole_order_checks.json'
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','carrier_cases':carrier['cases'],
                      'ordinary_s_pairs':carrier['ordinary_s_pairs'],
                      'generic_models':finite['models'],
                      'generic_counts':[(r['n'],r['k'],r['generic_initial_ideals']) for r in finite['fixtures']],
                      'gaussian_checks':finite['independent_coefficient_greedy_checks'],
                      'unisolvent_checks':finite['unisolvent_evaluation_greedy_checks'],
                      'characteristic_two_cone_equalities':finite['complete_characteristic_two_cone_equalities'],
                      'cylinder_cone_equalities':finite['complete_cylinder_cone_equalities'],
                      'special_sample_counterexamples':special['counterexamples'],
                      'output':str(out.relative_to(ROOT))}))


if __name__=='__main__':
    main()
