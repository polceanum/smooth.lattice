#!/usr/bin/env python3
"""Exact four-quotient compression and full quadratic-coordinate fans."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd, lcm
from pathlib import Path
import json
import hashlib

ROOT=Path(__file__).resolve().parents[1]
from run_signed_face_checks import field_basis, field_reduce
from run_projection_tag_checks import compare_fans, fan
from run_quadratic_coordinate_fan_checks import (TRIANGLE, joint_matching,
    quadratic_tensor, square_zero_fixtures)
from run_quadric_cross_checks import normalized
from run_planar_layer_obstruction_checks import (actual_fan as paired_fan,
                                                model as paired_model)
from run_canonical_fan_deflation_checks import marked_basis
from run_algebraic_incidence_controls import pair_checks
from run_planar_dyadic_bound_checks import convex_hull
from run_structural_complexity_checks import area2, clip, cross, dot, hyperbolic_cross


@lru_cache(None)
def universal(n,d):
    if not d:return ((),)
    return tuple((a,)+b for a in range(n) for b in universal(n//(a+1),d-1))


def common_samples(anchors,lo,hi,n,counts):
    poly=clip(list(TRIANGLE),(0,lo.denominator,-lo.numerator),1)
    if hi is not None:poly=clip(poly,(0,-hi.denominator,hi.numerator),1)
    normals=set()
    for a,b in combinations(anchors,2):
        v=tuple(x-y for x,y in zip(a,b));d=gcd(*v)
        v=tuple(x//d for x in v)
        if next(x for x in v if x)<0:v=tuple(-x for x in v)
        normals.add(v)
    cells=[poly]
    for normal in sorted(normals):
        cells=[v for cell in cells for sign in (-1,1)
               if (v:=clip(cell,normal,sign)) and area2(v)>0]
    samples=[]
    for cell in cells:
        triple=next(v for v in combinations(cell,3) if dot(v[0],cross(v[1],v[2])))
        rays=[]
        for v in triple:
            d=lcm(*(x.denominator for x in v))
            a=tuple(int(d*x) for x in v);g=gcd(*a)
            rays.append(tuple(x//g for x in a))
        w=tuple(sum(v[j] for v in rays) for j in range(3))
        assert min(w)>0 and max(w)<=6*n*n
        assert all(dot(w,q) for q in normals)
        assert lo<w[1]/Q(w[2]) and (hi is None or w[1]/Q(w[2])<hi)
        samples.append(w)
        counts['independent_conditional_sampler_integer_height_certificates']+=1
    return samples


def compressed(n,ev,p,counts,complete_ground=None):
    ground=universal(n,2)
    complete_ground=universal(n,3) if complete_ground is None else tuple(complete_ground)
    raw={j:tuple((j,*a) for a in ground) for j in (0,1)}
    bases={j:field_basis(raw[j],ev,p)[0] for j in (0,1)}
    pivots={j:field_basis(bases[j],ev,p)[1] for j in (0,1)}
    h=len(bases[0])+len(bases[1])-n
    assert h>=0
    def quotient(j,a):return field_reduce(ev(a),pivots[1-j],p)[0]
    qev={j:lambda a,j=j:quotient(j,a) for j in (0,1)}
    unique=all(sum(bool(evaluate(a)) for a in raw[j])==rank
               for j,evaluate,rank in ((0,ev,len(bases[0])),(1,ev,len(bases[1])),
                 (0,qev[0],len(bases[0])-h),(1,qev[1],len(bases[1])-h)))
    cuts=([Q(0),None] if unique else [Q(0)]+
          sorted({Q(b,a) for a,b in product(range(1,n),repeat=2)})+[None])
    strips=[]
    histories=[set() for _ in range(4)]
    for lo,hi in zip(cuts,cuts[1:]):
        u=lo+1 if hi is None else (lo+hi)/2
        state=[]
        for j in (0,1):
            ordered=sorted(raw[j],key=lambda a:(u*a[1]+a[2],a))
            B=tuple(field_basis(ordered,ev,p)[0])
            C=tuple(field_basis(ordered,qev[j],p)[0])
            assert set(C)<=set(B) and len(B)-len(C)==h
            state.extend((frozenset(B),frozenset(C)))
        key=tuple(state)
        for history,S in zip(histories,key):history.add(S)
        if strips and strips[-1]['key']==key:strips[-1]['hi']=hi
        else:strips.append({'lo':lo,'hi':hi,'key':key})
    assert len(strips)<=1+sum(len(history)-1 for history in histories)
    assert len(strips)<=1+8*n*(n-1).bit_length()
    patches=defaultdict(list);sample_bases=set();all_samples=set()
    for strip in strips:
        B0,C0,B1,C1=strip['key']; C=C0|C1
        anchors=tuple(sorted((B0-C0)|(B1-C1)))
        assert len(C)==n-h and len(anchors)==2*h
        if 0<h<n:
            assert (0,0,0) in C
            counts['proper_forced_span_contains_unit_and_is_not_ideal']+=1
        _,cp=field_basis(tuple(sorted(C)),ev,p)
        assert len(cp)==n-h
        def rem(a):return field_reduce(ev(a),cp,p)[0]
        left=tuple(sorted(B0-C0));right=tuple(sorted(B1-C1))
        _,hp=field_basis(left,rem,p);assert len(hp)==h
        def small_ev(a):
            r,co=field_reduce(rem(a),hp,p);assert not r
            return {left.index(b):v for b,v in co.items()}
        assert len(field_basis(right,small_ev,p)[0])==h
        lo,hi=strip['lo'],strip['hi']
        for w in common_samples(anchors,lo,hi,n,counts):
            selected=frozenset(C|set(field_basis(sorted(anchors,key=lambda a:(dot(w,a),a)),small_ev,p)[0]))
            direct=frozenset(field_basis(sorted(complete_ground,key=lambda a:(dot(w,a),a)),ev,p)[0])
            assert selected==direct
            sample_bases.add(selected);all_samples.add(w)
            counts['independent_conditional_sampler_vs_original_greedy']+=1
        inner=sorted({Q(b[2]-a[2],a[1]-b[1])
                      for group in (left,right) for a,b in combinations(group,2)
                      if a[1]!=b[1] and lo<Q(b[2]-a[2],a[1]-b[1]) and
                      (hi is None or Q(b[2]-a[2],a[1]-b[1])<hi)})
        local=[lo]+inner+[hi]
        for l,r in zip(local,local[1:]):
            u=l+1 if r is None else (l+r)/2
            L=tuple(sorted(left,key=lambda a:(u*a[1]+a[2],a)))
            U=tuple(sorted(right,key=lambda a:(u*a[1]+a[2],a)))
            pairs=joint_matching(L,U,h,small_ev,p,counts) if h else set()
            assert len(pairs)==h
            poly=clip(list(TRIANGLE),(0,l.denominator,-l.numerator),1)
            if r is not None:poly=clip(poly,(0,-r.denominator,r.numerator),1)
            cells=[poly]
            for i,j in sorted(pairs):
                q=tuple(x-y for x,y in zip(L[i],U[j]))
                cells=[v for cell in cells for sign in (-1,1) if (v:=clip(cell,q,sign))]
            for cell in cells:
                w=tuple(sum(v[j] for v in cell)/len(cell) for j in range(3))
                selected=frozenset(C|{L[i] if dot(w,L[i])<dot(w,U[j]) else U[j] for i,j in pairs})
                direct=frozenset(field_basis(sorted(complete_ground,key=lambda a:(dot(w,a),a)),ev,p)[0])
                assert selected==direct
                assert C<=direct<=B0|B1
                patches[selected].extend(cell)
                counts['independent_rank_h_cell_vs_original_greedy']+=1
        counts['fixed_rank_h_encodings']+=1
    result={}
    for B,vertices in patches.items():
        hull=convex_hull([(v[0],v[1]) for v in vertices])
        result[B]=[(a,b,1-a-b) for a,b in hull]
    assert sum(area2(poly) for poly in result.values())==1
    assert sample_bases==set(result)
    return result,{'n':n,'h':h,'four_planar_counts':list(map(len,histories)),
                   'fixed_encodings':len(strips),'initial_ideals':len(result),
                   'four_planar_models_have_unique_bases':unique,
                   'conditional_sample_count':len(all_samples),
                   'largest_sample_coordinate':max(max(w) for w in all_samples)},histories,all_samples


def algebra_controls(n,ev,p,constant,linear,counts):
    ground=universal(n,2)
    labels={j:tuple((j,*a) for a in ground) for j in (0,1)}
    B={j:field_basis(labels[j],ev,p)[0] for j in (0,1)}
    P={j:field_basis(B[j],ev,p)[1] for j in (0,1)}
    h=len(B[0])+len(B[1])-n
    def quo(a):return field_reduce(ev(a),P[0],p)[0]
    C,CP=field_basis(B[1],quo,p)
    vectors=[]
    for a in B[1]:
        if a in C:continue
        rem,co=field_reduce(quo(a),CP,p);assert not rem
        w=defaultdict(int,ev(a))
        for b,c in co.items():
            for j,v in ev(b).items():w[j]-=c*v
        vectors.append(normalized(w,p))
    _,WP=field_basis(tuple(range(len(vectors))),lambda j:vectors[j],p)
    assert len(WP)==h
    for w in vectors:
        assert not field_reduce(w,P[0],p)[0] and not field_reduce(w,P[1],p)[0]
    reference,RP=field_basis(B[0]+B[1],ev,p);assert len(reference)==n
    def move(v,i):
        rem,co=field_reduce(v,RP,p);assert not rem
        out=defaultdict(int)
        for a,c in co.items():
            successor=tuple(x+int(i==j) for j,x in enumerate(a))
            for j,d in ev(successor).items():out[j]+=c*d
        return normalized(out,p)
    for a in reference:
        v=ev(a)
        for i,j in combinations(range(3),2):
            assert move(move(v,i),j)==move(move(v,j),i)
            counts['independent_original_coordinate_commutations']+=1
        rel=defaultdict(int,move(move(v,0),0))
        for j,c in move(v,0).items():rel[j]+=linear*c
        for j,c in v.items():rel[j]+=constant*c
        assert not normalized(rel,p)
        counts['independent_actual_quadratic_relation_replays']+=1
    for w in vectors:
        for i in range(3):
            assert not field_reduce(move(w,i),WP,p)[0]
            counts['independent_intersection_is_actual_algebra_ideal']+=1
    counts['independent_four_cyclic_quotient_length_identities']+=1
    assert len(B[0])+len(B[1])+(len(B[0])-h)+(len(B[1])-h)==2*n
    return h


def confluent(m,h):
    n=m+2*h;M=n
    G=[1]
    for j in range(2,m+2):
        out=[0]*(len(G)+1)
        for k,v in enumerate(G):out[k]+=(1-j)*v;out[k+1]+=v
        G=out
    @lru_cache(None)
    def ev(a):
        x,y,z=a;E=y+M*z
        if x>=2:return {}
        if x==0:
            out={k:comb(E,k) for k in range(min(E,2*h-1)+1)}
            out.update({2*h+j-2:j**E for j in range(2,m+2)})
            return normalized(out,0)
        return normalized({k:sum(G[j]*comb(E,k-h-j) for j in range(len(G))
                                  if 0<=k-h-j<=E) for k in range(h,2*h)},0)
    return n,ev


def main():
    counts=Counter();rows=[]
    for n,d in product(range(1,11),(2,3)):
        assert set(universal(n,d))==set(hyperbolic_cross(n,d))
        counts['independent_exact_hyperbolic_ground_replays']+=1
    cases=list(square_zero_fixtures())
    for p,r,a,b in ((0,3,1,0),(0,3,0,0),(2,3,1,1)):
        n,ev=quadratic_tensor(p,r,a,b)
        cases.append(({'model':'quadratic','char':p,'constant':a,'linear':b},n,ev))
    for m,h in ((1,1),(3,1),(2,2),(3,2),(4,1),(2,3),(5,2)):
        n,ev=confluent(m,h);cases.append(({'model':'confluent','char':0,'m':m,'h':h},n,ev))
    for p in (0,7):
        from run_projection_tag_checks import point_model
        n,ev=point_model(((0,1,2),(1,2,3),(1,3,5)),p)
        cases.append(({'model':'split_radical_points','char':p,'constant':0,'linear':-1},n,ev))
    guard_indices=defaultdict(list);guard_carriers={}
    for k in (1,2):
        for swap in (False,True):
            perm=[j^int(swap) for j in range(2*k)]
            n,D,_,unit_ev,_=paired_model(k,perm,counts)
            def ev(a,unit_ev=unit_ev):
                value=unit_ev(a)
                return {} if value is None else {value:1}
            cases.append(({'model':'same_four_quotients_guard','char':2,'k':k,'swap':swap},n,ev))
            guard_carriers[k,swap]=D
    guard_histories={};guard_samples={}
    for label,n,ev in cases:
        p=label.get('char',label.get('characteristic',0))
        h=algebra_controls(n,ev,p,label.get('constant',0),label.get('linear',0),counts)
        carrier=guard_carriers.get((label.get('k'),label.get('swap')))
        compiled,row,histories,samples=compressed(n,ev,p,counts,carrier)
        assert h==row['h']
        if label['model']=='same_four_quotients_guard':
            k=label['k'];swap=label['swap'];perm=[j^int(swap) for j in range(2*k)]
            data=paired_fan(k,perm,'swapped' if swap else 'identity',counts)
            _,D,pairs,_,_=paired_model(k,perm,counts)
            fixed=D-{a for pair in pairs for a in pair}
            actual={frozenset(fixed|{tuple(a) for a in cell['chosen_relation_tips']}):
                    [tuple(Q(x) for x in v) for v in cell['positive_cone']]
                    for cell in data['cells']}
        else:actual=fan(n,ev,p,counts)
        compare_fans(actual,compiled,counts)
        for B,poly in actual.items():
            G,_=marked_basis(B,ev,p)
            w=tuple(sum(v[j] for v in poly)/len(poly) for j in range(3))
            counts['independent_original_reduced_basis_S_pairs']+=pair_checks(G,w,p)
        if label['model']=='confluent':
            assert h==label['h']
            lower={B for B in histories[0]}
            planar=universal(n,2)
            cuts=[Q(0)]+sorted({Q(b,a) for a,b in product(range(1,n),repeat=2)})+[None]
            free=set()
            for lo,hi in zip(cuts,cuts[1:]):
                u=lo+1 if hi is None else (lo+hi)/2
                B=frozenset((0,*a) for a in sorted(planar,key=lambda a:(u*a[0]+a[1],a))[:n])
                assert len(field_basis(tuple(B),ev,p)[0])==n
                free.add(B)
                counts['independent_confluent_free_prefix_rank_certificates']+=1
            assert free==lower
            assert lower<={frozenset(a for a in B if a[0]==0) for B in actual if all(a[0]==0 for a in B)}
            counts['independent_sharp_family_planar_fan_lifts']+=len(lower)
        if label['model']=='same_four_quotients_guard':
            k=label['k'];guard_histories[k,label['swap']]=histories
            guard_samples[k,label['swap']]=samples
            guard_indices[k].append(len(actual))
        rows.append(dict(row,label=label))
        print(rows[-1],flush=True)
    for k,values in guard_indices.items():
        assert guard_histories[k,False]==guard_histories[k,True]
        assert guard_samples[k,False]==guard_samples[k,True]
        assert values==[2,(k+1)**2]
        counts['same_four_planar_quotients_still_require_intersection_matrix']+=1
    sources=['docs/quadratic_intersection_fans.md','scripts/run_quadratic_intersection_checks.py',
             'docs/quadratic_coordinate_fan_compilation.md','docs/planar_dyadic_bound.md',
             'docs/mixed_rank_shadow_fans.md','scripts/run_signed_face_checks.py',
             'scripts/run_projection_tag_checks.py','scripts/run_quadratic_coordinate_fan_checks.py',
             'scripts/run_quadric_cross_checks.py','scripts/run_planar_layer_obstruction_checks.py',
             'scripts/run_canonical_fan_deflation_checks.py','scripts/run_algebraic_incidence_controls.py',
             'scripts/run_planar_dyadic_bound_checks.py','scripts/run_structural_complexity_checks.py',
             'scripts/run_oracle_enumeration_checks.py','scripts/run_interleaving_research.py',
             'scripts/run_square_zero_extension_checks.py']
    report={'all_checks_passed':True,'scope':'Four cyclic planar quotients fix n-h compulsory columns on O(n log n) vertical strips. The remaining exact encoding has rank h on 2h labels. Conditional weight lists are common to every compatible coupling and have certified integer height at most 6*n^2. This gives a parameterized class count, sharp planar order at fixed h, and no unrestricted G3 improvement or novelty certificate.',
            'rows':rows,'counts':dict(counts),
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    target=ROOT/'results/research/quadratic_intersection_checks.json'
    target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_checks_passed':True,'models':len(rows),'counts':dict(counts),
                      'report':str(target.relative_to(ROOT))}),flush=True)


if __name__=='__main__':main()
