#!/usr/bin/env python3
"""Exact controls for rational-function fans and uniform counting.

No growth fit, novelty claim, or unrestricted extremal result is inferred.
Generic finite-quotient cones are computed by carrier dependence and
Cramer's rule, rather than by pretending numerical samples are generic.
"""
from functools import cache
from fractions import Fraction
from itertools import combinations
from math import comb, gcd, isqrt
from pathlib import Path
import hashlib
import json

from run_interleaving_research import minimal_excluded_border
from run_initial_cone_height_checks import clean, leading, normal_form, subtract_shift, integer_ray
from run_extremal_initial_ideals_checks import adjacent_weights
from run_algebraic_incidence_controls import floor_sum_integer
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross, primitive

ROOT = Path(__file__).resolve().parents[1]


def coefficient(x, p):
    return x % p if p else Fraction(x)


def relations(k, p, polynomial=None, numerator=None):
    coefficients=polynomial or tuple((-1)**(k-j)*comb(k,j) for j in range(k+1))
    assert len(coefficients)==k+1 and coefficient(coefficients[0],p) and coefficient(coefficients[-1],p)
    result = [clean({(1, 1, 0): 1, (0, 0, 0): -1}, p)]
    for i in range(k + 1):
        g = {}
        for j,value in enumerate(numerator or (1,)):
            exponent=j-i
            g[(max(exponent,0),max(-exponent,0),0)]=-value
        for j in range(k + 1):
            exponent = j - i
            a = (max(exponent, 0), max(-exponent, 0), 1)
            g[a] = coefficients[j]
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


@cache
def polynomial_power(polynomial, exponent, p):
    def multiply(left,right):
        result=[0]*(len(left)+len(right)-1)
        for i,x in enumerate(left):
            for j,y in enumerate(right):
                result[i+j]+=x*y
        return tuple(x%p for x in result) if p else tuple(result)
    result=(1,)
    base=polynomial
    while exponent:
        if exponent & 1:
            result=multiply(result,base)
        base=multiply(base,base)
        exponent//=2
    return result


def cleared_columns(points, k, p, polynomial=None, numerator=None):
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
    degree = max((a-b+shift+k*(level-c)+(len(numerator)-1 if numerator else 0)*c for a, b, c in points), default=0)
    columns = {}
    for u in points:
        a, b, c = u
        offset, power = a-b+shift, k*(level-c)
        vector = [0] * (degree+1)
        coefficients=(polynomial_power(tuple(polynomial),level-c,p) if polynomial is not None
                      else tuple((-1)**(power-j)*comb(power,j) for j in range(power+1)))
        if numerator is not None:
            factor=polynomial_power(tuple(numerator),c,p)
            expanded=[0]*(len(coefficients)+len(factor)-1)
            for i,x in enumerate(coefficients):
                for j,y in enumerate(factor):
                    expanded[i+j]+=x*y
            coefficients=tuple(coefficient(x,p) for x in expanded)
        for j,value in enumerate(coefficients):
            vector[offset+j] = coefficient(value,p)
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
def generic_cone(B,k,p,polynomial=None,numerator=None):
    border=tuple(minimal_excluded_border(B))
    columns=cleared_columns(B+border,k,p,polynomial,numerator)
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


@cache
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


def carrier_prefix_rank(k,w,t):
    """Exact generic carrier prefix rank, also on lex-refined carrier walls."""
    assert k>=1 and min(w)>0 and all(isinstance(x,int) for x in w)
    if t<0:
        return 0
    alpha,beta,gamma=w
    m=k*alpha//(alpha+beta)
    a,b=k-m,m+1
    total=1+t//alpha+t//beta+t//gamma
    for step,width in ((alpha,a-1),(beta,b-1)):
        length=min(width,max(0,(t-gamma)//step))
        total+=floor_sum_integer(length,gamma,-step,t-step)
    return total


def bounded_two_axis_solution(total,step,vertical,width):
    if width<1 or total<step+vertical:
        return None
    divisor=gcd(step,vertical)
    if total%divisor:
        return None
    modulus=vertical//divisor
    residue=(total//divisor*pow(step//divisor,-1,modulus))%modulus if modulus>1 else 0
    first=residue if residue>=1 else residue+modulus
    last=min(width,(total-vertical)//step)
    if first>last:
        return None
    assert first+modulus>last,'refined cost is unique in the relevant exponent box'
    return first,(total-step*first)//vertical


def integer_unrank(r,k,w):
    """Research oracle, with increasing (integer cost, exponent tuple) order.

    This does not change the public smooth-number CLI or certify logarithms.
    """
    assert r>=1 and k>=1 and min(w)>0 and all(isinstance(x,int) for x in w)
    k=min(k,r)
    refined=tuple(r**3*x+y for x,y in zip(w,(r*r,r,1)))
    low,high=0,(r-1)*min(refined)
    queries=0
    while low<high:
        middle=(low+high)//2
        queries+=1
        if carrier_prefix_rank(k,refined,middle)>=r:
            high=middle
        else:
            low=middle+1
    t=low
    alpha,beta,gamma=refined
    m=k*alpha//(alpha+beta)
    a,b=k-m,m+1
    candidates=set()
    for axis,step in enumerate(refined):
        if t%step==0:
            u=tuple(t//step if i==axis else 0 for i in range(3))
            candidates.add(u)
    for axis,step,width in ((0,alpha,a-1),(1,beta,b-1)):
        solution=bounded_two_axis_solution(t,step,gamma,width)
        if solution is not None:
            horizontal,vertical=solution
            candidates.add((horizontal,0,vertical) if axis==0 else (0,horizontal,vertical))
    assert len(candidates)==1
    u=candidates.pop()
    assert max(u)<r and dot(refined,u)==t and is_standard(u,a,b)
    assert carrier_prefix_rank(k,refined,t-1)==r-1
    assert carrier_prefix_rank(k,refined,t)==r
    return u,{'refined_weights':refined,'cost':t,'binary_rank_queries':queries}


def explicit_standard_points(k,w,max_cost):
    alpha,beta,gamma=w
    m=k*alpha//(alpha+beta)
    a,b=k-m,m+1
    result={(i,0,0) for i in range(max_cost//alpha+1)}
    result|={(0,j,0) for j in range(1,max_cost//beta+1)}
    for c in range(1,max_cost//gamma+1):
        result.add((0,0,c))
        result.update((i,0,c) for i in range(1,a) if i*alpha+c*gamma<=max_cost)
        result.update((0,j,c) for j in range(1,b) if j*beta+c*gamma<=max_cost)
    return result


def integer_access_controls():
    from itertools import product
    prefix=0;access=0;large=[];max_queries=0
    for k in range(1,8):
        for w in product(range(1,5),repeat=3):
            points=explicit_standard_points(k,w,40)
            counts=[sum(dot(w,u)<=t for u in points) for t in range(41)]
            for t,expected in enumerate(counts):
                assert carrier_prefix_rank(k,w,t)==expected
                prefix+=1
    weights=((1,1,1),(1,2,3),(5,3,2),(2,7,1),(13,11,7),(1,13,2),
             (6,10,15),(100,1,7),(1,100,7))
    for k in (1,2,3,4,7,23,1000):
        for w in weights:
            # This independently sorts the infinite-carrier standard points
            # below the pure-axis upper bound. It does not use rank queries.
            points=sorted(explicit_standard_points(k,w,82*min(w)),key=lambda u:(dot(w,u),u))
            for r in (1,2,5,17,50,83):
                u,certificate=integer_unrank(r,k,w)
                assert u==points[r-1]
                access+=1
                max_queries=max(max_queries,certificate['binary_rank_queries'])
    def equal_degree_count(k,d):
        if d<0:
            return 0
        a=k-k//2;b=k//2+1
        A=min(a-1,max(0,d-1));B=min(b-1,max(0,d-1))
        return 1+3*d+(A+B)*d-A*(A+1)//2-B*(B+1)//2
    def equal_degree_unrank(k,r):
        a=k-k//2;b=k//2+1
        low,high=0,r
        while low<high:
            mid=(low+high)//2
            if equal_degree_count(k,mid)>=r:
                high=mid
            else:
                low=mid+1
        d=low;position=r-equal_degree_count(k,d-1)
        if d==0:
            return (0,0,0)
        y_mixed=min(b-1,d-1)+1
        if position<=y_mixed:
            j=position-1
            return (0,j,d-j)
        position-=y_mixed
        if position==1:
            return (0,d,0)
        position-=1
        x_mixed=min(a-1,d-1)
        if position<=x_mixed:
            return (position,0,d-position)
        assert position==x_mixed+1
        return (d,0,0)
    for r,k in ((10**12,1),(10**12,7),(10**12,10**6),(10**12,10**10),
                (10**18,10**16),(10**18,10**24),(10**18+123,37)):
        u,certificate=integer_unrank(r,k,(1,1,1))
        assert u==equal_degree_unrank(k,r)
        large.append({'rank':r,'k':k,'weights':[1,1,1],'monomial':u,
                      'independent_check':'closed ordinary-degree Hilbert counts and direct within-degree order',
                      **certificate})
    return {'exhaustive_small_prefix_counts':prefix,'independent_sort_unrank_checks':access,
            'maximum_small_binary_queries':max_queries,'large_closed_formula_checks':len(large),
            'large_witnesses':large,
            'scope':'Exact integer costs and prescribed lex ties, generic rational-function quotients regular at zero and infinity only. Classical Euclidean floor sums; no bit-complexity, finite-special-sample or prime-log certification claim.'}


def whole_wall_controls():
    rows=[];sensitive=0;insensitive=0;monotone=0;incident=0;off=0
    for n,k in ((5,1),(5,2),(6,2),(8,2),(8,3),(8,4)):
        visited=full_fan(n,k)
        cylinder=full_fan(n,None)
        ground=tuple(u for u in hyperbolic_cross(n,3) if not u[0]*u[1])
        rays={integer_ray(v) for poly,w in visited.values() for v in poly if min(v)>0}
        per_wall={i:[] for i in range(1,k)}
        insensitive_rows=[]
        for ray in sorted(rays):
            matches=[i for i in range(1,k) if (k-i)*ray[0]==i*ray[1]]
            if not matches:
                off+=1
                continue
            assert len(matches)==1
            i=matches[0];A=k-i;B=i
            strip_a,strip_b=A,B+1
            canonical=[u for u in ground if is_standard(u,strip_a,strip_b)]
            canonical.sort(key=lambda u:(dot(u,ray),u))
            T=dot(canonical[n-1],ray)
            t0=A*ray[0]
            face=[basis for basis in visited
                  if all(dot(z,ray)>=0 for z in generic_cone(tuple(sorted(basis)),k,0)[1])]
            moments=[sum(u[2] for u in basis) for basis in face]
            assert max(moments)>min(moments)
            event={'ray':ray,'nth_cost':T,'carrier_head_cost':t0+ray[2],
                   'Z_moment_before':max(moments),'Z_moment_after':min(moments),
                   'drop':max(moments)-min(moments)}
            if T>=t0+ray[2]:
                sensitive+=1
                H={(j,0,0) for j in range(A//2+1)}|{(0,j,0) for j in range(1,B//2+1)}
                assert 2*len(H)>=k
                for basis in face:
                    C=max(u[2] for u in basis)
                    cheap={(h[0],h[1],c) for h in H for c in range(C//2+1)}
                    assert cheap<=basis and all(dot(ray,u)<T for u in cheap)
                    assert k*(C+1)<=4*n and k*sum(u[2] for u in basis)<4*n*n
                    incident+=1
                event['incident_bases_checked']=len(face)
            else:
                insensitive+=1
                cylinder_face=[basis for basis in cylinder
                               if all(dot(z,ray)>=0 for z in generic_cone(tuple(sorted(basis)),None,0)[1])]
                assert set(face)==set(cylinder_face)
                insensitive_rows.append(event)
            event['sensitive']=T>=t0+ray[2]
            per_wall[i].append(event)
        saved=[]
        for i,events in per_wall.items():
            events.sort(key=lambda e:Fraction(e['ray'][2],e['ray'][0]))
            for previous,current in zip(events,events[1:]):
                assert previous['Z_moment_after']>=current['Z_moment_before']
                monotone+=1
            total=sum(e['drop'] for e in events if e['sensitive'])
            assert k*total<=4*n*n
            saved.append({'wall_i':i,'sensitive_integer_drop_sum':total,'events':events})
        rows.append({'n':n,'k':k,'walls':saved,'insensitive_cylinder_vertices':insensitive_rows})
    return {'models':len(rows),'sensitive_vertices':sensitive,'insensitive_cylinder_vertex_equalities':insensitive,
            'incident_completed_rectangles':incident,'whole_wall_monotonicity_pairs':monotone,
            'off_carrier_vertices':off,'fixtures':rows}


def denominator_controls():
    rows=[];s_pairs=0;replays=0;prefix=0;cone_equalities=0
    polynomials=[]
    for k in (2,3,5):
        polynomials.extend(((1,)+(0,)*(k-1)+(1,),tuple(1 for _ in range(k+1)),
                            (1,1)+(0,)*(k-2)+(1,)))
    for P in polynomials:
        k=len(P)-1
        for m in range(k):
            w=(2*m+1,2*k-2*m-1,k+1)
            a,b=sector(k,w)
            for p in (0,2,3):
                G=relations(k,p,P)
                heads=[leading(g,w) for g in G]
                assert (a,0,1) in heads and (0,b,1) in heads
                for g in G:
                    columns=cleared_columns(g,k,p,P)
                    assert all(coefficient(sum(c*columns[u][j] for u,c in g.items()),p)==0
                               for j in range(len(next(iter(columns.values())))))
                    replays+=1
                for i,j in combinations(range(len(G)),2):
                    u,v=heads[i],heads[j];common=tuple(max(x,y) for x,y in zip(u,v))
                    factor=pow(G[i][u],-1,p) if p else 1/G[i][u]
                    first={tuple(x+y-z for x,y,z in zip(t,common,u)):coefficient(c*factor,p)
                           for t,c in G[i].items()}
                    factor=pow(G[j][v],-1,p) if p else 1/G[j][v]
                    pair=subtract_shift(first,G[j],tuple(x-y for x,y in zip(common,v)),factor,p)
                    assert not normal_form(pair,G,w,p)
                    s_pairs+=1
                n=6
                ground=tuple(u for u in hyperbolic_cross(n,3) if not u[0]*u[1])
                columns=cleared_columns(ground,k,p,P)
                assert truncated_greedy(ground,columns,w,n,p)==formula_basis(n,k,w,ground)
                prefix+=1
                rows.append({'polynomial_low_to_high':P,'sector':m,'weights':w,'characteristic':p})
    for n,k in ((4,2),(5,2),(6,3)):
        for P in polynomials:
            if len(P)!=k+1:
                continue
            for p in (0,2,3):
                for basis,(poly,w) in full_fan(n,k).items():
                    alternate,normals=generic_cone(tuple(sorted(basis)),k,p,P)
                    assert alternate and area2(alternate)==area2(poly)
                    assert all(dot(z,v)>=0 for v in poly for z in normals)
                    cone_equalities+=1
    numerator_rows=[];numerator_prefix=0;numerator_cones=0;numerator_pairs=0
    for P in polynomials:
        k=len(P)-1
        for Q in ((0,1),(0,)*k+(1,)):
            # P(0)=1 makes every monomial numerator coprime over every field.
            for p in (0,2,3):
                w=(1,2*k+1,k+3)
                G=relations(k,p,P,Q);heads=[leading(g,w) for g in G]
                for i,j in combinations(range(len(G)),2):
                    u,v=heads[i],heads[j];common=tuple(max(x,y) for x,y in zip(u,v))
                    factor=pow(G[i][u],-1,p) if p else 1/G[i][u]
                    first={tuple(x+y-z for x,y,z in zip(t,common,u)):coefficient(c*factor,p)
                           for t,c in G[i].items()}
                    factor=pow(G[j][v],-1,p) if p else 1/G[j][v]
                    pair=subtract_shift(first,G[j],tuple(x-y for x,y in zip(common,v)),factor,p)
                    assert not normal_form(pair,G,w,p)
                    numerator_pairs+=1
                for n in (4,6):
                    ground=tuple(u for u in hyperbolic_cross(n,3) if not u[0]*u[1])
                    columns=cleared_columns(ground,k,p,P,Q)
                    for w in ((1,2*k+1,k+3),(2*k+1,1,2*k+3),(3,4,11)):
                        if k*w[0]%(w[0]+w[1])==0:
                            continue
                        assert truncated_greedy(ground,columns,w,n,p)==formula_basis(n,k,w,ground)
                        numerator_prefix+=1
                if k in (2,3):
                    for basis,(poly,w) in full_fan(4,k).items():
                        alternate,normals=generic_cone(tuple(sorted(basis)),k,p,P,Q)
                        assert alternate and area2(alternate)==area2(poly)
                        assert all(dot(z,v)>=0 for v in poly for z in normals)
                        numerator_cones+=1
                numerator_rows.append({'P':P,'Q':Q,'characteristic':p})
    return {'denominator_polynomials':len(polynomials),'carrier_field_sector_cases':len(rows),
            'ordinary_s_pairs':s_pairs,'cleared_relation_replays':replays,
            'independent_generic_prefix_checks':prefix,'complete_generic_cone_equalities':cone_equalities,
            'nonconstant_numerator_cases':len(numerator_rows),'nonconstant_numerator_s_pairs':numerator_pairs,
            'nonconstant_numerator_prefix_checks':numerator_prefix,'nonconstant_numerator_cone_equalities':numerator_cones,
            'numerator_fixtures':numerator_rows,
            'fixtures':rows,'scope':'Monic degree-k integer denominators with constant one, sparse and geometric-sum cases, and coprime nonconstant monomial numerators of degree one or k, in several characteristics. The uniform proof covers all coprime Q/P regular at zero and infinity.'}



def main():
    if not __debug__:
        raise RuntimeError('Run without -O: assertions are required.')
    carrier=carrier_controls()
    finite=finite_controls()
    special=special_sample_controls()
    walls=whole_wall_controls()
    denominators=denominator_controls()
    access=integer_access_controls()
    sources=('docs/pole_order_fans.md','scripts/run_pole_order_checks.py',
             'docs/extremal_initial_ideals.md','docs/rational_shadow_bound.md',
             'scripts/run_initial_cone_height_checks.py','scripts/run_extremal_initial_ideals_checks.py',
             'scripts/run_structural_complexity_checks.py','scripts/run_interleaving_research.py',
             'scripts/run_algebraic_incidence_controls.py')
    report={'schema_version':1,'status':'passed',
            'command':'python3 -B scripts/run_pole_order_checks.py',
            'uniform_proof':'docs/pole_order_fans.md',
            'scope':'Generic n-point quotients of XY=1,P(X)*Z=Q(X), coprime Q/P of degree Q<=degree P with P(0) nonzero, in every characteristic. A uniform O(n^2*log(n+1)) bound closes the whole generic family. Numerical special samples are not assumed generic. Laurent identities, generic evaluation, monotone support slopes and Euclidean floor sums are classical inputs.',
            'carrier':carrier,'generic_finite_quotients':finite,'special_sample_counterexamples':special,
            'sensitive_wall_integer_budgets':walls,'general_denominator_controls':denominators,
            'exact_integer_rank_unrank':access,
            'uniform_family_bound':'O(n^2*log(n+1))',
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
                      'sensitive_vertices':walls['sensitive_vertices'],
                      'general_polynomial_cone_equalities':denominators['complete_generic_cone_equalities'],
                      'integer_prefix_checks':access['exhaustive_small_prefix_counts'],
                      'integer_unrank_checks':access['independent_sort_unrank_checks'],
                      'large_closed_formula_checks':access['large_closed_formula_checks'],
                      'output':str(out.relative_to(ROOT))}))


if __name__=='__main__':
    main()
