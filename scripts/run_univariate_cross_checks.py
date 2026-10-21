#!/usr/bin/env python3
"""Complete two-stage rank/dual compilation for univariate cross quotients."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import prod
from pathlib import Path
from random import Random
import hashlib
import json

from run_coordinate_slice_budget_checks import edges_from_polygons, primitive, profile
from run_curve_carrier_shadow_checks import curve_model
from run_interleaving_research import minimal_excluded_border
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_planar_dyadic_bound_checks import convex_hull
from run_primitive_wall_checks import threshold
from run_projection_tag_checks import compare_fans, fan, polygon_key
from run_quadric_cross_checks import apply, gluing_model, normalized, rank, sharp_model
from run_signed_face_checks import field_basis, normals as border_normals
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]
TRIANGLE = tuple(tuple(Q(i == j) for i in range(3)) for j in range(3))


def center(poly, avoid=()):
    for k in range(1, 2+max(1, len(avoid))*max(1, len(poly))):
        coeffs = [(k+1)**j for j in range(len(poly))]
        w = tuple(sum(c*v[i] for c, v in zip(coeffs, poly))/sum(coeffs) for i in range(3))
        if all(dot(q, w) for q in avoid):
            return w
    raise AssertionError('nonzero polynomials have only finitely many roots')


def chambers(normals):
    cells = [list(TRIANGLE)]
    for q in sorted(normals):
        out = []
        for poly in cells:
            sides = min(dot(q, v) for v in poly), max(dot(q, v) for v in poly)
            if sides[0] < 0 < sides[1]:
                for sign in (-1, 1):
                    part = clip(poly, q, sign)
                    assert area2(part)
                    out.append(part)
            else:
                out.append(poly)
        cells = out
    assert sum(map(area2, cells)) == 1
    return cells


def compress_compile(n, ev, p, groups, counts, actual=None):
    E = tuple(sorted(set().union(*map(set, groups))))
    assert sum(map(len, groups)) == len(E) and rank([ev(a) for a in E], p) == n
    bound = 1+max(max(a) for a in E)
    inventory = {primitive(tuple(x-y for x, y in zip(a, b)))
                 for group in groups for a, b in combinations(group, 2)}
    inventory = {q for q in inventory if min(q) < 0 < max(q)}
    cuts = chambers(inventory)
    m = sum(rank([ev(a) for a in group], p) for group in groups)
    delta = m-n
    merged, ledger = defaultdict(list), []
    for chamber in cuts:
        w = center(chamber, inventory)
        lists = [field_basis(sorted(group, key=lambda a: (dot(w, a), a)), ev, p)[0]
                 for group in groups]
        M = tuple(a for L in lists for a in L)
        assert len(M) == m and len(set(M)) == m

        @lru_cache(None)
        def dual_independent(T):
            if len(T) > delta:
                return False
            return rank([ev(a) for a in M if a not in T], p) == n

        retained = []
        for L in lists:
            B = set()
            for a in reversed(L):
                if dual_independent(frozenset(B|{a})):
                    B.add(a)
            assert len(B) <= delta
            retained.extend(B)
        K = tuple(sorted(retained))
        assert len(K) <= len(groups)*delta
        counts['independent_primal_and_dual_group_compressions'] += len(groups)
        comparisons = {primitive(tuple(x-y for x, y in zip(a, b))) for a, b in combinations(K, 2)}
        first = center(chamber, comparisons|inventory)
        queue, visited, local_patches = [first], set(), []
        chamber_edges = [(a, b) for a, b in zip(chamber, chamber[1:]+chamber[:1])]
        while queue:
            w = queue.pop()
            D = set()
            for a in sorted(K, key=lambda a: (dot(w, a), a), reverse=True):
                if dual_independent(frozenset(D|{a})):
                    D.add(a)
            assert len(D) == delta
            D = frozenset(D)
            if D in visited:
                continue
            visited.add(D)
            B = frozenset(set(M)-D)
            assert len(B) == n
            assert B == frozenset(field_basis(sorted(E, key=lambda a: (dot(w, a), a)), ev, p)[0])
            local_normals = {tuple(x-y for x, y in zip(a, b))
                             for a in D for b in K if b not in D
                             and dual_independent(frozenset(D-{a}|{b}))}
            poly = list(chamber)
            for q in local_normals:
                poly = clip(poly, q, 1)
            assert area2(poly)
            assert all(dot(q, w) > 0 for q in local_normals)
            if actual is not None:
                N, _ = border_normals(tuple(sorted(B)), ev, p)
                independent = list(chamber)
                for q in N:
                    independent = clip(independent, q, 1)
                assert polygon_key(independent) == polygon_key(poly)
                counts['independent_entire_dual_border_patch_equalities'] += 1
            merged[B].append(poly)
            local_patches.append(poly)
            for a, b in zip(poly, poly[1:]+poly[:1]):
                # Boundary portions of this internal-order chamber are not crossed.
                if any(dot(q, a) == dot(q, b) == 0 for q in inventory):
                    continue
                if any(a[i] == b[i] == 0 for i in range(3)):
                    continue
                q = next(q for q in local_normals if dot(q, a) == dot(q, b) == 0)
                neighbor, _ = exact_adjacent_weight(bound, a, b, q)
                assert all(dot(primitive(tuple(x-y for x, y in zip(c, d))), neighbor) != 0
                           for L in lists for c, d in combinations(L, 2))
                queue.append(tuple(map(Q, neighbor)))
        assert sum(map(area2, local_patches)) == area2(chamber)
        assert len(local_patches) <= 1+len(comparisons)*(len(comparisons)+1)//2
        counts['complete_internal_order_chambers'] += 1
        ledger.append({'chamber': [[str(c) for c in v] for v in chamber],
                       'compressed_columns': m, 'corank': delta,
                       'dual_candidates': sorted(K), 'dual_patches': len(local_patches)})
    combined = {}
    for B, patches in merged.items():
        vertices = convex_hull([(v[0], v[1]) for patch in patches for v in patch])
        combined[B] = [(a, b, 1-a-b) for a, b in vertices]
        assert area2(combined[B]) == sum(map(area2, patches))
    assert sum(map(area2, combined.values())) == 1
    if actual is not None:
        compare_fans(actual, combined, counts)
    return combined, {'internal_comparison_planes': len(inventory),
                      'internal_chambers': len(cuts), 'compressed_union_rank': n,
                      'compressed_union_columns': m, 'corank': delta, 'chambers': ledger}


def point_model(points, p, linear=False):
    n = len(points)
    def ev(a):
        return normalized({i: prod(pow(v, e) for v, e in zip(point, a))
                           for i, point in enumerate(points)}, p)
    M = [tuple(normalized({i: point[j]}, p) for i, point in enumerate(points)) for j in range(3)]
    fs=((0,1),)*3 if linear else ((0,-1,1),)*3
    return n, ev, M, p, fs, {
        'kind':'linear_coordinate_cross_points' if linear else 'coordinate_grid_line_points',
        'points':points,'field':p}


def gaussian_model():
    points = (((2,1),(0,1),(0,1)), ((0,1),(3,1),(0,1)))
    def mul(a,b):
        return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]
    def power(a,k):
        out=(Q(1),Q(0))
        for _ in range(k):out=mul(out,a)
        return out
    def ev(a):
        out={}
        for j, point in enumerate(points):
            value=(Q(1),Q(0))
            for v,k in zip(point,a):value=mul(value,power(v,k))
            out.update({2*j+i:c for i,c in enumerate(value) if c})
        return out
    M=[]
    for coord in range(3):
        columns=[]
        for point in points:
            a,b=point[coord];j=len(columns)
            columns.extend((normalized({j:a,j+1:b},0),normalized({j:-b,j+1:a},0)))
        M.append(tuple(columns))
    return 4,ev,M,0,((1,0,1),)*3,{'kind':'irreducible_quadratic_grid_Q_i','field':0}


def extended_sharp(L,s):
    old,ev0,M0,p,_,meta=sharp_model(L)
    n=s+4*L
    def index(j):return j+s-2 if j>=2 else j
    def shift(v):return {index(j):c for j,c in v.items()}
    def ev(a):
        if not a[1] and not a[2]:
            return {a[0]:Q(1)} if a[0]<s else {}
        return shift(ev0(a))
    M=[]
    for i in range(3):
        cols=[{} for _ in range(n)]
        if i==0:
            for j in range(s-1):cols[j]={j+1:Q(1)}
        else:
            for j in range(2):cols[j]=shift(M0[i][j])
        for j in range(2,old):cols[index(j)]=shift(M0[i][j])
        M.append(tuple(cols))
    return n,ev,M,p,((0,0,1),(0,1),(0,1)),{
        **meta,'kind':'sharp_product_extended_core','core_length':s,'length':n}


def lifted_gluing(p):
    old,ev0,M0,_,_,meta=gluing_model(2,p)
    n=2*old
    def ev(a):
        q,r=divmod(a[0],2)
        return {r*old+j:c for j,c in ev0((q,a[1],a[2])).items()}
    M=[]
    for i in range(3):
        cols=[]
        for block in range(2):
            for j in range(old):
                if i==0:
                    cols.append({old+j:1} if block==0 else dict(M0[0][j]))
                else:
                    cols.append({block*old+k:c for k,c in M0[i][j].items()})
        M.append(tuple(cols))
    return n,ev,M,p,((0,0,0,0,1),(0,1),(0,1)),{
        **meta,'kind':'four_layer_coupled_gluing','field':p}


def polynomial_matrix(M,f,p):
    n=len(M);out=[]
    for j in range(n):
        vector={j:1};value=defaultdict(int)
        for c in f:
            for i,d in vector.items():value[i]+=c*d
            vector=apply(M,vector,p)
        out.append(normalized(value,p))
    return tuple(out)


def algebra_controls(model, groups, counts):
    n,ev,M,p,fs,_=model
    degrees=tuple(len(f)-1 for f in fs);S=prod(degrees)
    F=[polynomial_matrix(A,f,p) for A,f in zip(M,fs)]
    for i,j in combinations(range(3),2):
        for k in range(n):
            assert apply(M[i],M[j][k],p)==apply(M[j],M[i][k],p)
            assert not apply(F[i],F[j][k],p)
            counts['independent_commuting_and_cross_annihilation_columns']+=1
    C=[ev(a) for a in product(*(range(d) for d in degrees))]
    mY=[rank(list(A),p) for A in F]
    kernel=[mY[i]-rank([apply(F[i],v,p) for v in F[i]],p) for i in range(3)]
    assert max(kernel)<=S
    allY=[v for A in F for v in A]
    assert rank(C+allY,p)==n
    assert sum(mY)-rank(allY,p)<=2*S
    for i in range(3):
        generated=[]
        others_indices=[j for j in range(3) if j!=i]
        for a,b in product(range(degrees[others_indices[0]]),range(degrees[others_indices[1]])):
            exponent=[0,0,0]
            exponent[others_indices[0]]=a;exponent[others_indices[1]]=b
            v=apply(F[i],ev(tuple(exponent)),p)
            for _ in range(n):
                generated.append(v)
                v=apply(M[i],v,p)
        assert rank(generated,p)==mY[i]
        counts['independent_univariate_branch_module_generations']+=1
        others=[v for j,A in enumerate(F) if j!=i for v in A]
        intersection=mY[i]+rank(others,p)-rank(list(F[i])+others,p)
        assert intersection<=kernel[i]
        assert rank(C+list(F[i])+[ev(a) for a in groups[i+1]],p)==rank(C+list(F[i]),p)
        counts['independent_branch_kernel_intersection_and_tail_containments']+=1
    delta=sum(rank([ev(a) for a in group],p) for group in groups)-n
    assert 0<=delta<=6*S
    if all(f[:-1]==tuple(0 for _ in f[:-1]) for f in fs):
        assert delta<=3*S
    for a in set().union(*map(set,groups)):
        for i in range(3):
            b=tuple(v+(i==j) for j,v in enumerate(a))
            assert apply(M[i],ev(a),p)==ev(b)
            counts['independent_raw_monomial_successor_replays']+=1
    return {'degrees':degrees,'box_size':S,'branch_ranks':mY,
            'branch_kernel_dimensions':kernel,'compressed_corank':delta}


def abstract_controls(counts):
    rows=[]
    E=((0,0,0),(2,0,0),(0,1,0),(0,3,0),(1,0,1),(1,0,3))
    groups=(E[:2],E[2:4],E[4:])
    for n in (3,4,5):
        rng=Random(83+n)
        cols=[tuple(int(i==j) for i in range(n)) if j<n else
              tuple(rng.randrange(101) for _ in range(n)) for j in range(len(E))]
        ev=lambda a:{i:c for i,c in enumerate(cols[E.index(a)]) if c}
        original={}
        for B in combinations(E,n):
            if rank([ev(a) for a in B],101)!=n:continue
            N={tuple(x-y for x,y in zip(b,a)) for a in B for b in E if b not in B
               and rank([ev(v) for v in set(B)-{a}|{b}],101)==n}
            poly=list(TRIANGLE)
            for q in N:poly=clip(poly,q,1)
            if area2(poly):original[frozenset(B)]=poly
        assert sum(map(area2,original.values()))==1
        compiled,ledger=compress_compile(n,ev,101,groups,counts)
        compare_fans(original,compiled,counts)
        rows.append({'rank':n,'ground':E,'initial_profiles':len(compiled),'compilation':ledger})
    return rows


def inverse_overlap_guards(counts):
    rows=[]
    for n,p in ((4,0),(6,0),(5,2)):
        _,ev0,_,_,_=curve_model(n,p,True)
        ev=lambda a: ev0(a) if not a[2] else {}
        zero=(0,0,0)
        groups=([zero],[(j,0,0) for j in range(1,n)],
                [(0,j,0) for j in range(1,n)],[])
        assert [rank([ev(a) for a in group],p) for group in groups]==[1,n-1,n-1,0]
        actual=fan(n,ev,p,counts)
        compiled,ledger=compress_compile(n,ev,p,groups,counts,actual)
        assert ledger['corank']==n-1 and len(compiled)==n
        assert ev((1,1,0))==ev(zero)
        counts['independent_fixed_curve_unbounded_corank_guards']+=1
        rows.append({'length':n,'field':p,'curve':'xy=1,z=0',
                     'initial_ideals':len(compiled),'compressed_corank':n-1})
    return rows


def main():
    counts=Counter();abstract=abstract_controls(counts);rows=[]
    models=[extended_sharp(L,L+3) for L in (1,2,3)]
    for p in (0,2,3,101):
        base=gluing_model(2,p)
        models.append((*base[:4],((0,0,1),(0,1),(0,1)),base[-1]))
    models.append(lifted_gluing(101))
    points=((2,0,0),(3,1,0),(0,2,1),(1,0,4))
    models.extend([point_model(points,0),point_model(points,7),gaussian_model()])
    axes=((2,0,0),(3,0,0),(0,4,0),(0,0,5))
    models.extend([point_model(axes,0,True),point_model(axes,7,True)])
    for model in models:
        n,ev,M,p,fs,metadata=model
        degrees=tuple(len(f)-1 for f in fs)
        E=tuple(sorted(a for a in hyperbolic_cross(n,3)
                       if sum(v>=d for v,d in zip(a,degrees))<=1))
        groups=[[] for _ in range(4)]
        for a in E:
            group=next((i+1 for i,(v,d) in enumerate(zip(a,degrees)) if v>=d),0)
            groups[group].append(a)
        invariants=algebra_controls(model,groups,counts)
        actual=fan(n,ev,p,counts)
        compiled,ledger=compress_compile(n,ev,p,groups,counts,actual)
        S=prod(degrees)
        R=sum(degrees[i]*degrees[j] for i,j in combinations(range(3),2))
        for B in actual:
            assert len(minimal_excluded_border(tuple(B)))<=3+S+R
            counts['independent_constant_minimal_generator_bounds']+=1
        if metadata['kind']=='linear_coordinate_cross_points':
            assert ledger['internal_chambers']==1 and len(actual)==3
        edges=edges_from_polygons(actual,counts)
        value=sum(threshold(primitive(tuple(a-b for a,b in zip(profile(B),profile(C)))))
                  for B,C in map(tuple,edges))
        if metadata['kind']=='sharp_product_extended_core':
            L,s=metadata['L'],metadata['core_length']
            assert len(actual)==(L+1)**2 and value==2*L*(L+1)**2
            unit=ev((0,0,0));minimal=[]
            for A in M:
                v=unit
                for j in range(1,n+1):
                    v=apply(A,v,p)
                    if not v:
                        minimal.append(j);break
            assert minimal==[s,2*L+1,2*L+1] and min(minimal)>=3
            for pair in ((0,1),(0,2)):
                box=[tuple(int(i==pair[0])*a+int(i==pair[1])*b for i in range(3))
                     for a,b in product(range(2),repeat=2)]
                assert rank([ev(a) for a in box],p)==4
            counts['independent_growing_coordinate_degree_sharp_families']+=1
        rows.append({**metadata,'length':n,'initial_ideals':len(actual),
                     'undirected_walls':len(edges),'total_primitive_cost':value,
                     'rank_invariants':invariants,'compilation':ledger,
                     'full_cones':[{'standard_set':sorted(B),
                                    'positive_cone':[[str(c) for c in v] for v in poly]}
                                   for B,poly in actual.items()]})
        print('complete univariate cross:',metadata['kind'],n,p,len(actual),
              'corank',invariants['compressed_corank'],flush=True)
    guards=inverse_overlap_guards(counts)
    sources=['docs/univariate_cross_fans.md','scripts/run_univariate_cross_checks.py',
             'docs/curve_carrier_shadows.md','docs/quadric_cross_fans.md',
             'scripts/run_curve_carrier_shadow_checks.py',
             'scripts/run_coordinate_slice_budget_checks.py','scripts/run_oracle_enumeration_checks.py',
             'scripts/run_planar_dyadic_bound_checks.py','scripts/run_primitive_wall_checks.py',
             'scripts/run_projection_tag_checks.py','scripts/run_quadric_cross_checks.py',
             'scripts/run_signed_face_checks.py','scripts/run_structural_complexity_checks.py']
    out={'schema_version':1,'all_checks_passed':True,
         'scope':'Uniform univariate pair-product class: sharp quadratic fan and cubic full-cost '
                 'orders at fixed degrees, with arbitrary joint quotient ranks. Exact whole-cone '
                 'two-stage primal/dual compilation, not finite growth extrapolation. '
                 'No unrestricted G3 endpoint, polynomial lower separation or novelty established.',
         'counters':dict(counts),'abstract_matroid_controls':abstract,
         'complete_actual_quotient_fans':rows,
         'fixed_curve_unbounded_corank_guards':guards,
         'source_sha256':{s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}}
    target=ROOT/'results/research/univariate_cross_checks.json'
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_checks_passed':True,'counters':dict(counts),
                      'report':str(target.relative_to(ROOT))}),flush=True)


if __name__=='__main__':
    main()
