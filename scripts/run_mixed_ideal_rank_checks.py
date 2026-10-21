#!/usr/bin/env python3
"""Exact polynomial-prefix, actual-ideal and full-fan mixed-rank controls."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import prod
from pathlib import Path
import hashlib
import json

from run_coordinate_slice_budget_checks import downset, edges_from_polygons, primitive, profile
from run_interleaving_research import minimal_excluded_border
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_planar_dyadic_bound_checks import convex_hull
from run_primitive_wall_checks import threshold
from run_projection_tag_checks import compare_fans, fan, polygon_key
from run_quadric_cross_checks import apply, gluing_model, normalized, rank
from run_signed_face_checks import field_basis, normals as border_normals
from run_structural_complexity_checks import (area2, clip, determinant, dot,
    finite_staircases, hyperbolic_cross, univariate_product_mod, univariate_power_mod)
from run_univariate_cross_checks import (TRIANGLE, center, chambers, extended_sharp,
                                        gaussian_model, point_model, polynomial_matrix)

ROOT = Path(__file__).resolve().parents[1]


def tensor_model(r, s, p, extra):
    reference = tuple(product(range(r), range(s)))
    index = {a:j for j,a in enumerate(reference)}
    M = []
    for coord in range(3):
        columns = []
        for a,b in reference:
            out = defaultdict(int)
            if coord in (0,2) and a+1 < r: out[index[a+1,b]] += 1
            if coord in (1,2) and b+1 < s: out[index[a,b+1]] += 1
            if coord == 2 and extra and a+1 < r and b+1 < s: out[index[a+1,b+1]] += 1
            columns.append(normalized(out,p))
        M.append(tuple(columns))
    @lru_cache(None)
    def ev(a):
        v={index[0,0]:1}
        for A,k in zip(M,a):
            for _ in range(k):v=apply(A,v,p)
        return v
    return r*s,ev,M,p,((0,1),)*3,{
        'kind':'coupled_tensor','r':r,'s':s,'extra_uv':extra,'field':p}


def long_gluing(s,p):
    old,ev0,M0,_,_,metadata=gluing_model(2,p)
    n=old+s-2
    unit_index=next(iter(ev0((0,0,0))))
    x_index=next(iter(ev0((1,0,0))))
    other=[j for j in range(old) if j not in (unit_index,x_index)]
    mapping={unit_index:0,x_index:1,**{j:s+k for k,j in enumerate(other)}}
    def shift(v):return {mapping[j]:c for j,c in v.items()}
    @lru_cache(None)
    def ev(a):
        if not a[1] and not a[2]:return {a[0]:1} if a[0]<s else {}
        return shift(ev0(a))
    M=[]
    for i in range(3):
        columns=[{} for _ in range(n)]
        if not i:
            for j in range(s-1):columns[j]={j+1:1}
        else:
            columns[0]=shift(M0[i][unit_index])
            columns[1]=shift(M0[i][x_index])
        for j in other:columns[mapping[j]]=shift(M0[i][j])
        M.append(tuple(columns))
    return n,ev,M,p,((0,1),)*3,{
        **metadata,'kind':'long_axis_fixed_mixed_rank','core_length':s,'field':p}


def shifted_rectangle(p):
    n=9
    modulus=(0,)*n+(1,)
    X=(1,1)+(0,)*(n-2)
    Y=(1,0,0,1)+(0,)*(n-4)
    Z=(0,)*n
    images=(X,Y,Z)
    M=[tuple(normalized(dict(enumerate(univariate_product_mod(
        image,tuple(int(i==j) for i in range(n)),modulus))),p) for j in range(n))
       for image in images]
    @lru_cache(None)
    def ev(a):
        v=(1,)+(0,)*(n-1)
        for image,k in zip(images,a):
            v=univariate_product_mod(v,univariate_power_mod(image,k,modulus),modulus)
        return normalized(dict(enumerate(v)),p)
    return n,ev,M,p,((1,-2,1),(1,-2,1),(0,1)),{
        'kind':'coherent_selected_body_rectangle','field':p,'expected_mixed_rank':1}


def polynomial_columns(model):
    n,mono,M,p,fs,_=model
    degrees=tuple(len(f)-1 for f in fs)
    F=tuple(polynomial_matrix(A,f,p) for A,f in zip(M,fs))
    @lru_cache(None)
    def univariate(e,i):
        power,b=divmod(e,degrees[i])
        values=[1]
        for _ in range(power):
            out=[0]*(len(values)+degrees[i])
            for j,c in enumerate(values):
                for k,d in enumerate(fs[i]):out[j+k]+=c*d
            values=out
        return tuple((j+b,c) for j,c in enumerate(values) if c)
    @lru_cache(None)
    def expansion(a):
        out={tuple(item[0] for item in triple):prod(item[1] for item in triple)
             for triple in product(*(univariate(e,i) for i,e in enumerate(a)))}
        assert out[a]==1
        assert all(b==a or (b!=a and all(i<=j for i,j in zip(b,a))) for b in out)
        return out
    @lru_cache(None)
    def evaluate(a):
        out=defaultdict(int)
        for b,c in expansion(a).items():
            for j,d in mono(b).items():out[j]+=c*d
        return normalized(out,p)
    def operator_value(a):
        quotient,remainder=zip(*(divmod(e,d) for e,d in zip(a,degrees)))
        v=mono(tuple(remainder))
        for A,k in zip(F,quotient):
            for _ in range(k):v=apply(A,v,p)
        return v
    return degrees,F,expansion,evaluate,operator_value


def compile_fan(n,mono,evaluate,p,carrier,body,counts,actual,S,t):
    # Zero polynomial columns are global loops and can be discarded first.
    carrier=[tuple(a for a in group if evaluate(a)) for group in carrier]
    body=tuple(a for a in body if evaluate(a))
    E=tuple(sorted(set(body).union(*map(set,carrier))))
    assert rank([evaluate(a) for a in E],p)==n
    inventory={primitive(tuple(x-y for x,y in zip(a,b)))
               for group in carrier for a,b in combinations(group,2)}
    inventory={q for q in inventory if min(q)<0<max(q)}
    order_cells=chambers(inventory)
    m_carrier=sum(rank([evaluate(a) for a in group],p) for group in carrier)
    m=m_carrier+len(body)
    delta=m-n
    assert 0<=delta<=3*S+3*t+len(body)
    bound=1+max(max(a) for a in E)
    Hn=tuple(hyperbolic_cross(n,3))
    patches=defaultdict(list);ledger=[]
    for chamber in order_cells:
        w=center(chamber,inventory)
        lists=[field_basis(sorted(group,key=lambda a:(dot(w,a),a)),evaluate,p)[0]
               for group in carrier]
        M=tuple(a for L in lists for a in L)+body
        assert len(M)==m and len(set(M))==m
        @lru_cache(None)
        def independent(T):
            return len(T)<=delta and rank([evaluate(a) for a in M if a not in T],p)==n
        retained=set(body)
        for L in lists:
            D=set()
            for a in reversed(L):
                if independent(frozenset(D|{a})):D.add(a)
            assert len(D)<=delta
            retained.update(D)
        K=tuple(sorted(retained))
        assert len(K)<=4*delta+len(body)
        comparisons={primitive(tuple(x-y for x,y in zip(a,b))) for a,b in combinations(K,2)}
        queue=[center(chamber,inventory|comparisons)]
        visited=set();local=[]
        while queue:
            w=queue.pop()
            D=set()
            for a in sorted(K,key=lambda a:(dot(w,a),a),reverse=True):
                if independent(frozenset(D|{a})):D.add(a)
            assert len(D)==delta
            D=frozenset(D)
            if D in visited:continue
            visited.add(D)
            B=frozenset(set(M)-D)
            assert len(B)==n and downset(B)
            assert B==frozenset(field_basis(sorted(Hn,key=lambda a:(dot(w,a),a)),mono,p)[0])
            normals={tuple(x-y for x,y in zip(a,b)) for a in D for b in K if b not in D
                     and independent(frozenset(D-{a}|{b}))}
            poly=list(chamber)
            for q in normals:poly=clip(poly,q,1)
            assert area2(poly) and all(dot(q,w)>0 for q in normals)
            original,_=border_normals(tuple(sorted(B)),mono,p)
            replay=list(chamber)
            for q in original:replay=clip(replay,q,1)
            assert polygon_key(poly)==polygon_key(replay)
            counts['independent_entire_polynomial_dual_original_border_patches']+=1
            local.append(poly);patches[B].append(poly)
            for a,b in zip(poly,poly[1:]+poly[:1]):
                if any(dot(q,a)==dot(q,b)==0 for q in inventory):continue
                if any(a[i]==b[i]==0 for i in range(3)):continue
                q=next(q for q in normals if dot(q,a)==dot(q,b)==0)
                neighbor,_=exact_adjacent_weight(bound,a,b,q)
                queue.append(tuple(map(Q,neighbor)))
        assert sum(map(area2,local))==area2(chamber)
        assert len(local)<=1+len(comparisons)*(len(comparisons)+1)//2
        counts['complete_carrier_order_chambers']+=1
        counts['independent_primal_dual_carrier_compressions']+=4
        ledger.append({'positive_chamber':[[str(c) for c in v] for v in chamber],
                       'dual_candidate_labels':sorted(K),'full_cone_patches':len(local)})
    merged={}
    for B,parts in patches.items():
        vertices=convex_hull([(v[0],v[1]) for part in parts for v in part])
        poly=[(a,b,1-a-b) for a,b in vertices]
        assert area2(poly)==sum(map(area2,parts))
        merged[B]=poly
    assert sum(map(area2,merged.values()))==1
    compare_fans(actual,merged,counts)
    return {'carrier_order_planes':len(inventory),'carrier_order_chambers':len(order_cells),
            'nonzero_body_columns_kept':len(body),'carrier_compressed_columns':m_carrier,
            'combined_columns':m,'dual_rank':delta,
            'largest_dual_candidate_count':max(len(v['dual_candidate_labels']) for v in ledger),
            'global_encoding':len(order_cells)==1,'chambers':ledger}


def mixed(a,degrees):
    return sum(i>=d for i,d in zip(a,degrees))>=2


def border(D):
    return minimal_excluded_border(tuple(sorted(D))) if D else [(0,0,0)]


def basis_controls(B,degrees,t,body,expansion,evaluate,p,counts):
    ordered=sorted(B,key=lambda a:(sum(a),a))
    assert rank([evaluate(a) for a in ordered],p)==len(B)
    for a in B:assert set(expansion(a))<=B
    assert determinant([[expansion(a).get(b,0) for b in ordered] for a in ordered])==1
    counts['independent_unitriangular_basis_determinants']+=1
    chosen={a for a in B if mixed(a,degrees)}
    assert len(chosen)<=t and chosen<=body
    anchors=[tuple(degrees[k] if k in pair else 0 for k in range(3))
             for pair in combinations(range(3),2)]
    for a in chosen:
        for h in anchors:
            if all(i>=j for i,j in zip(a,h)):
                assert prod(i-j+1 for i,j in zip(a,h))<=len(chosen)
                counts['independent_relative_body_divisor_boxes']+=1
    heads=set(map(tuple,border(B)))
    relative=set()
    for h in anchors:
        D={tuple(a-b for a,b in zip(alpha,h)) for alpha in B
           if all(a>=b for a,b in zip(alpha,h))}
        assert len(D)<=t and downset(D)
        relative|={tuple(a+b for a,b in zip(h,beta)) for beta in border(D)}
    assert {a for a in heads if mixed(a,degrees)}<=relative
    S=prod(degrees);R=sum(degrees[i]*degrees[j] for i,j in combinations(range(3),2))
    assert sum(not mixed(a,degrees) for a in heads)<=S+R
    counts['independent_mixed_relative_border_covers']+=1
    return len(chosen)


def run_model(model,counts):
    n,mono,M,p,fs,metadata=model
    degrees,F,expansion,evaluate,op=polynomial_columns(model)
    S=prod(degrees)
    Tcols=[apply(F[i],v,p) for i,j in combinations(range(3),2) for v in F[j]]
    t=rank(Tcols,p)
    for A in M:
        assert rank(Tcols+[apply(A,v,p) for v in Tcols],p)==t
        counts['independent_actual_mixed_ideal_invariance']+=1
    for i,j in combinations(range(3),2):
        for k in range(n):
            assert apply(M[i],M[j][k],p)==apply(M[j],M[i][k],p)
            counts['independent_commuting_matrix_columns']+=1
    anchors=[tuple(degrees[k] if k in pair else 0 for k in range(3))
             for pair in combinations(range(3),2)]
    body={tuple(a+b for a,b in zip(h,v)) for h in anchors for v in hyperbolic_cross(t,3)} if t else set()
    assert len(body)<=3*len(hyperbolic_cross(t,3)) if t else not body
    Hn=tuple(sorted(hyperbolic_cross(n,3)))
    carrier=[[] for _ in range(4)]
    for a in Hn:
        if mixed(a,degrees):continue
        group=next((i+1 for i,(v,d) in enumerate(zip(a,degrees)) if v>=d),0)
        carrier[group].append(a)
    ground=set(body).union(*map(set,carrier))
    D=set(Hn)
    for a in ground:D.update(product(*(range(v+1) for v in a)))
    for a in D:
        assert evaluate(a)==op(a)
        counts['independent_expansion_operator_polynomial_columns']+=1
        for i in range(3):
            b=tuple(v+(i==j) for j,v in enumerate(a))
            assert apply(M[i],mono(a),p)==mono(b)
            counts['independent_raw_monomial_successor_replays']+=1
    changed_singletons=sum(bool(mono(a))!=bool(evaluate(a)) for a in D)
    ordinary_body_excess=rank(Tcols+[mono(a) for a in body],p)-t
    counts['guard_changed_singleton_ranks']+=changed_singletons
    counts['guard_ordinary_monomial_body_rank_above_actual_ideal']+=ordinary_body_excess
    assert rank(Tcols+[evaluate(a) for a in body],p)==t
    def quotient_rank(columns):return rank(Tcols+list(columns),p)-t
    branch=[quotient_rank(A) for A in F]
    kernels=[branch[i]-quotient_rank([apply(F[i],v,p) for v in F[i]]) for i in range(3)]
    assert max(kernels)<=S
    assert sum(branch)-quotient_rank([v for A in F for v in A])<=2*S
    quotient_groups=[quotient_rank([evaluate(a) for a in group]) for group in carrier]
    assert sum(quotient_groups)<=n-t+3*S
    assert rank(Tcols+[evaluate(a) for group in carrier for a in group],p)==n
    assert sum(rank([evaluate(a) for a in group],p) for group in carrier)<=n+3*S+3*t
    counts['independent_quotient_branch_and_lifted_corank_bounds']+=1
    actual=fan(n,mono,p,counts)
    max_body=0
    for B,poly in actual.items():
        max_body=max(max_body,basis_controls(B,degrees,t,body,expansion,evaluate,p,counts))
        w=center(poly)
        order=sorted(D,key=lambda a:(dot(w,a),a))
        positions={a:i for i,a in enumerate(order)}
        for a in order:
            assert all(positions[b]<=positions[a] for b in expansion(a))
            counts['independent_prefix_triangular_span_inclusions']+=1
        ordinary=field_basis(order,mono,p)[0]
        transformed=field_basis(order,evaluate,p)[0]
        assert ordinary==transformed
        counts['independent_full_prefix_rank_equalities']+=len(order)
    if n<=6:
        for B in finite_staircases(n,3):
            B=frozenset(B)
            if rank([mono(a) for a in B],p)==n:
                basis_controls(B,degrees,t,body,expansion,evaluate,p,counts)
                counts['independent_exhaustive_order_ideal_basis_controls']+=1
    compilation=compile_fan(n,mono,evaluate,p,carrier,body,counts,actual,S,t)
    if degrees==(1,1,1):
        assert compilation['carrier_order_planes']==0 and compilation['global_encoding']
        counts['independent_global_linear_factor_encodings']+=1
    if metadata['kind']=='coherent_selected_body_rectangle':
        assert t==1 and max_body==1
        rectangle=frozenset((a,b,0) for a,b in product(range(3),repeat=2))
        assert rectangle in actual
        counts['independent_positive_error_coherent_body_rectangles']+=1
    if metadata['kind']=='long_axis_fixed_mixed_rank':
        assert t==2
    edges=edges_from_polygons(actual,counts)
    V=sum(threshold(primitive(tuple(a-b for a,b in zip(profile(B),profile(C)))))
          for B,C in map(tuple,edges))
    return {**metadata,'length':n,'univariate_factors':fs,'degrees':degrees,'mixed_ideal_rank':t,
            'mixed_ground_columns':len(body),'mixed_ground_outside_universal':len(body-set(Hn)),
            'changed_monomial_polynomial_singleton_ranks':changed_singletons,
            'ordinary_body_rank_above_actual_ideal':ordinary_body_excess,
            'largest_selected_body':max_body,'quotient_branch_ranks':branch,
            'quotient_branch_kernels':kernels,'initial_ideals':len(actual),
            'undirected_walls':len(edges),'total_primitive_cost':V,'compilation':compilation,
            'full_cones':[{'standard_set':sorted(B),
                           'positive_cone':[[str(c) for c in v] for v in poly]}
                          for B,poly in actual.items()]}


def main():
    counts=Counter();models=[]
    models.extend([long_gluing(3,0),long_gluing(8,0),long_gluing(4,2)])
    models.extend(tensor_model(2,2,p,True) for p in (0,2,101))
    models.extend(tensor_model(3,2,p,False) for p in (0,3))
    points=((1,2,3),(2,3,5),(3,5,7))
    for p in (0,7):
        n,ev,M,_,_,meta=point_model(points,p)
        for fs in (((-1,1),(-2,1),(-3,1)),
                   ((1,-2,1),(-2,1),(-3,1))):
            models.append((n,ev,M,p,fs,{**meta,'kind':'translated_nonlocal_positive_error'}))
    models.extend(shifted_rectangle(p) for p in (0,2,3))
    models.append(gaussian_model())
    models.append(extended_sharp(1,4))
    models.append(point_model(((2,0,0),(3,0,0),(0,4,0),(0,0,5)),0,True))
    rows=[]
    for model in models:
        row=run_model(model,counts);rows.append(row)
        print('complete mixed-rank:',row['kind'],row['length'],row['field'],
              't',row['mixed_ideal_rank'],'g',row['initial_ideals'],
              'K',row['compilation']['largest_dual_candidate_count'],flush=True)
    assert counts['guard_changed_singleton_ranks']>0
    assert counts['guard_ordinary_monomial_body_rank_above_actual_ideal']>0
    sources=['docs/mixed_ideal_rank_fans.md','scripts/run_mixed_ideal_rank_checks.py',
             'docs/univariate_cross_fans.md','scripts/run_univariate_cross_checks.py',
             'scripts/run_coordinate_slice_budget_checks.py','scripts/run_interleaving_research.py',
             'scripts/run_oracle_enumeration_checks.py','scripts/run_planar_dyadic_bound_checks.py',
             'scripts/run_primitive_wall_checks.py','scripts/run_projection_tag_checks.py',
             'scripts/run_quadric_cross_checks.py','scripts/run_signed_face_checks.py',
             'scripts/run_structural_complexity_checks.py']
    out={'schema_version':1,'all_checks_passed':True,
         'scope':'Actual nonzero mixed-ideal rank theorem, unitriangular polynomial-column '
                 'ground and complete dual fan encoding. Linear-factor fan counts depend '
                 'only on mixed rank, not ambient quotient length. No sharp t exponent, '
                 'unrestricted G3 improvement, polynomial separation or novelty established.',
         'counters':dict(counts),'complete_actual_fans':rows,
         'source_sha256':{s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}}
    target=ROOT/'results/research/mixed_ideal_rank_checks.json'
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_checks_passed':True,'counters':dict(counts),
                      'report':str(target.relative_to(ROOT))}),flush=True)


if __name__=='__main__':
    main()
