#!/usr/bin/env python3
"""Exact fan padding, arbitrary relation propagation and weighted-reuse guards."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import prod
from pathlib import Path
import hashlib
import json

from run_coordinate_slice_budget_checks import downset, edges_from_polygons, primitive, profile
from run_algebraic_incidence_controls import pair_checks
from run_interleaving_research import minimal_excluded_border
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_primitive_wall_checks import divisors, threshold
from run_projection_tag_checks import compare_fans
from run_quadric_cross_checks import normalized, rank
from run_signed_face_checks import field_basis, field_reduce
from run_socle_circuit_checks import (TRIANGLE, factor_fan, multiplication_controls,
                                      recross_model, scalar)
from run_structural_complexity_checks import (area2, clip, dot, univariate_power_mod,
                                              univariate_product_mod)
from run_univariate_cross_checks import center, gaussian_model, point_model

ROOT=Path(__file__).resolve().parents[1]


def shift(a,v):
    return tuple(x+y for x,y in zip(a,v))


def maximal(D):
    return {a for a in D if all(tuple(x+int(i==j) for j,x in enumerate(a)) not in D
                               for i in range(len(a)))}


def generic_model(D,E,rows,char):
    D,E=set(D),tuple(sorted(E))
    assert set(E)<=D and downset(D)
    vectors=[normalized({j:row.get(a,0) for j,a in enumerate(E)},char) for row in rows]
    _,pivots=field_basis(tuple(range(len(vectors))),lambda j:vectors[j],char)
    rho=len(pivots)
    positions={a:j for j,a in enumerate(E)}
    for row in rows:
        for i in range(3):
            moved={}
            for a,c in row.items():
                b=tuple(x+int(i==j) for j,x in enumerate(a))
                if b in D and scalar(c,char):
                    assert b in positions
                    moved[positions[b]]=scalar(c,char)
            assert not field_reduce(moved,pivots,char)[0]
    free={j for j in range(len(E)) if j not in pivots}
    reference=tuple(sorted((D-set(E))|{E[j] for j in free}))
    index={a:j for j,a in enumerate(reference)}
    block={a:field_reduce({j:1},pivots,char)[0] for j,a in enumerate(E)}
    columns={a:{index[a]:scalar(1,char)} for a in D-set(E)}
    for a in E:
        assert set(block[a])<=free
        columns[a]={index[E[j]]:c for j,c in block[a].items()}
    for row in rows:
        replay=defaultdict(int)
        for a,c in row.items():
            for j,d in columns[a].items():replay[j]+=c*d
        assert not normalized(replay,char)
    return {'D':D,'S':set(E),'fixed':D-set(E),'basis':reference,'polynomials':rows,
            'rho':rho,'n':len(reference),'char':char,'ev':lambda a:columns.get(a,{}),
            'block_ev':lambda a:block.get(a,{}),'variable':E}


def full_fan(A,counts):
    n,p=A['n'],A['char'];ev=A['block_ev'];ground=A['variable']
    queue=[(1,n+1,(n+1)**2)];out={}
    while queue:
        w=queue.pop()
        chosen,pivots=field_basis(sorted(ground,key=lambda a:(dot(w,a),a)),ev,p)
        assert len(chosen)==len(A['S'])-A['rho']
        B=frozenset(A['fixed']|set(chosen))
        if B in out:continue
        assert len(B)==n and downset(B)
        normals=set()
        for head in minimal_excluded_border(tuple(B)):
            rem,co=field_reduce(ev(head),pivots,p);assert not rem
            replay=defaultdict(int)
            for j,c in A['ev'](head).items():replay[j]+=c
            for a,c in co.items():
                for j,d in A['ev'](a).items():replay[j]-=c*d
                q=tuple(x-y for x,y in zip(head,a))
                assert dot(w,q)>0;normals.add(q)
            assert not normalized(replay,p)
            counts['independent_full_original_border_normal_forms']+=1
        poly=list(TRIANGLE)
        for q in normals:poly=clip(poly,q,1)
        assert area2(poly)>0;out[B]=poly
        for a,b in zip(poly,poly[1:]+poly[:1]):
            if any(a[i]==b[i]==0 for i in range(3)):continue
            q=next(q for q in normals if dot(q,a)==dot(q,b)==0)
            neighbor,_=exact_adjacent_weight(n,a,b,q)
            queue.append(neighbor)
    assert sum(area2(poly) for poly in out.values())==1
    counts['complete_independent_quotient_fans']+=1
    return out


def projected_count(D,S):
    return len({tuple(x for i,x in enumerate(a) if i not in S) for a in D})


def padded(A,v,counts):
    Dv=set().union(*(set(divisors(shift(a,v))) for a in maximal(A['D'])))
    box=tuple(product(*(range(x+1) for x in v)))
    assert Dv=={shift(a,b) for a in A['D'] for b in box}
    assert maximal(Dv)=={shift(a,v) for a in maximal(A['D'])}
    count=0
    for k in range(4):
        for S in combinations(range(3),k):
            count+=projected_count(A['D'],S)*prod(v[i] for i in S)
    assert count==len(Dv)
    counts['independent_full_extrusion_projection_formulas']+=1
    rows=[{shift(a,v):c for a,c in row.items()} for row in A['polynomials']]
    B=generic_model(Dv,{shift(a,v) for a in A['S']},rows,A['char'])
    assert B['rho']==A['rho'] and B['n']==len(Dv)-A['rho']
    return B


def basis_image(B,A,Av,v):
    return frozenset(Av['fixed']|{shift(a,v) for a in B&A['S']})


def compare_padding(A,v,counts,exhaustive=True):
    Av=padded(A,v,counts)
    multiplication_controls(A,counts);multiplication_controls(Av,counts)
    before=full_fan(A,counts);after=full_fan(Av,counts)
    delta=tuple(x-y-A['rho']*z for x,y,z in zip(profile(Av['D']),profile(A['D']),v))
    expected={basis_image(B,A,Av,v):poly for B,poly in before.items()}
    compare_fans(after,expected,counts)
    for B in before:
        assert profile(basis_image(B,A,Av,v))==shift(profile(B),delta)
        counts['independent_initial_profile_translation_equalities']+=1
    basis_count=0
    if exhaustive:
        k=len(A['S'])-A['rho']
        for selected in combinations(A['variable'],k):
            if rank([A['block_ev'](a) for a in selected],A['char'])!=k:continue
            B=frozenset(A['fixed']|set(selected));Bv=basis_image(B,A,Av,v)
            assert len(Bv)==Av['n'] and downset(Bv)==downset(B)
            assert rank([Av['block_ev'](shift(a,v)) for a in selected],A['char'])==k
            assert profile(Bv)==shift(profile(B),delta)
            basis_count+=1
            counts['independent_exhaustive_all_basis_profile_translations']+=1
    C=set.intersection(*(set(B) for B in before))
    Cv=set.intersection(*(set(B) for B in after))
    assert Cv==Av['fixed']|{shift(a,v) for a in C&A['S']}
    assert A['n']-len(C)==Av['n']-len(Cv)
    edges=edges_from_polygons(before,counts);edgesv=edges_from_polygons(after,counts)
    assert {frozenset(basis_image(B,A,Av,v) for B in e) for e in edges}==edgesv
    costs=Counter()
    for B,D in map(tuple,edges):
        diff=tuple(x-y for x,y in zip(profile(B),profile(D)))
        bv,dv=basis_image(B,A,Av,v),basis_image(D,A,Av,v)
        assert tuple(x-y for x,y in zip(profile(bv),profile(dv)))==diff
        costs[threshold(primitive(diff))]+=1
        counts['independent_primitive_cost_and_edge_multiple_translations']+=1
    return Av,after,{'field':A['char'],'shift':v,'original_length':A['n'],
                     'padded_length':Av['n'],'relation_rank':A['rho'],
                     'initial_ideals':len(before),'undirected_walls':len(edges),
                     'profile_translation':delta,'global_common_core_size':len(Cv),
                     'global_common_core_deficit':Av['n']-len(Cv),
                     'exhaustively_checked_vector_bases':basis_count,
                     'wall_cost_multiplicities':dict(sorted(costs.items())),
                     'total_primitive_cost':sum(c*m for c,m in costs.items())}


def arbitrary_relations(counts):
    D={(x,y,z) for x in range(2) for y in range(4-x) for z in range(4-x-y)}
    D.add((0,0,4))  # Overlapping rows now mix different total degrees.
    E=tuple(sorted(maximal(D)))
    assert len(E)==7
    coefficients=((1,1,0,1,0,1,0),(0,1,1,0,1,0,1),(1,0,1,1,1,1,1))
    rows=[{a:c for a,c in zip(E,row) if c} for row in coefficients]
    out=[]
    for p in (0,2,3,101):
        A=generic_model(D,E,rows,p)
        for v in ((0,2,1),(1,1,2)):
            _,_,record=compare_padding(A,v,counts)
            record['model']='overlapping_socle_relation_space';out.append(record)
    # A dependent row and a vanishing characteristic-two row test actual rank.
    rows=[{E[0]:1,E[1]:1},{E[1]:2,E[2]:2},{E[0]:1,E[1]:1}]
    for p in (0,2,3):
        A=generic_model(D,E,rows,p)
        _,_,record=compare_padding(A,(0,1,1),counts)
        record['model']='field_dependent_relation_rank';out.append(record)
    return out


def original_length(s):
    return (32*s**3+30*s*s+10*s+12)//3


def padded_length(s,A,B=0):
    return original_length(s)+(8*s*s+4)*A+(4*s+4)*B+2*A*B


def recross_controls(counts):
    out=[]
    for s,p in product(range(1,5),(0,2,3)):
        old,tips,upper=recross_model(s,p)
        A=generic_model(old['D'],old['S'],old['polynomials'],p)
        assert A['n']==original_length(s)
        Av,actual,record=compare_padding(A,(0,s,0),counts)
        assert Av['n']==padded_length(s,s)==(56*s**3+30*s*s+22*s+12)//3
        compare_fans(actual,factor_fan(Av),counts)
        assert len(actual)==4*s and record['global_common_core_deficit']==2*s
        shifted_upper=[shift(a,(0,s,0)) for a in upper]
        profiles={B:profile(B) for B in actual}
        witnesses=[]
        for i in range(1,2*s):
            u,t=Q(2*i+1,2),Q(i*i+i);w=(t,u,Q(1))
            best=min(dot(w,v) for v in profiles.values())
            incident=[B for B,value in profiles.items() if dot(w,value)==best]
            assert len(incident)==4
            assert {tuple(a in B for a in shifted_upper) for B in incident}==set(product((False,True),repeat=2))
            common=set.intersection(*(set(B) for B in incident))
            assert len(common)==Av['n']-2
            ev=shift(upper[0],(0,s,0));fv=shift(upper[1],(0,s,0))
            volumes=[prod(c+1 for c in a) for a in (ev,fv)]
            assert volumes==[2*(s+2)*(4*s*s+1),2*(s+1)*(4*s*s+2)]
            for j in (i,i+1):
                q=(1,-2*j,j*j)
                assert threshold(q)==2*j*j+2*j+1
            witnesses.append({'index':i,'u':str(u),'t':str(t),'common_face_core_deficit':2,
                              'primitive_costs':[2*j*j+2*j+1 for j in (i,i+1)]})
            counts['independent_same_pair_persistent_crossing_witnesses']+=1
        assert len(witnesses)==2*s-1
        record.update({'model':'padded_recrossing_pair','s':s,'upper_pair':shifted_upper,
                       'upper_divisor_volumes':volumes,'persistent_crossings':witnesses,
                       'sharp_two_envelope_upper_bound':2*s-1})
        out.append(record)
    return out


def symbolic_large_controls(counts):
    out=[]
    for mode,values in (('linear_padding',(256,512,1024,2048)),
                        ('near_full_boxes',(32,64,128,256))):
        for k in values:
            s=k if mode=='linear_padding' else k**4
            A=s if mode=='linear_padding' else k**5
            N=padded_length(s,A)
            e=(1,A+1,4*s*s);f=(1,A,4*s*s+1)
            ve,vf=(prod(c+1 for c in a) for a in (e,f))
            assert max(ve,vf)<=N
            samples=[]
            for i in (s,min(s+1,2*s-1),2*s-1):
                # Polynomial identity, valid for arbitrary integer j:
                # t_i - (2*j*u_i-j*j) = (j-i)*(j-i-1).
                left=[i*i+i,-(2*i+1),1]
                right=[i*(i+1),-i-(i+1),1]
                assert left==right
                nu=2*i*i+2*i+1
                assert nu**5>N**3
                assert (N-nu)**7>=N**3
                samples.append({'index':i,'primitive_cost':nu,
                                'above_N_to_three_fifths':True,
                                'below_N_minus_N_to_three_sevenths':True})
                counts['exact_large_intermediate_cost_witnesses']+=1
            if mode=='near_full_boxes':
                assert N-ve==original_length(s)+2*A-16*s*s-4
                assert N-vf==original_length(s)-8*s*s-4
            out.append({'mode':mode,'parameter':k,'s':s,'y_padding':A,'length':N,
                        'upper_divisor_volumes':[ve,vf],
                        'global_common_core_deficit':2*s,'face_common_core_deficit':2,
                        'persistent_crossings_formula':2*s-1,'whole_fan_count_formula':4*s,
                        'selected_exact_witnesses':samples,
                        'scope':'Closed formulas and factorized crossing witnesses only; no full fan or quotient matrices enumerated at this length.'})
    return out


def truncated_relations(D,generators):
    rows=[]
    for g in generators:
        for a in D:
            row={shift(a,b):c for b,c in g.items() if shift(a,b) in D}
            if row:rows.append(row)
    return rows


def propagating_controls(counts):
    out=[]
    for t,p in product((2,3),(0,2,3)):
        D=set(product(range(2),range(2),range(2*t)))
        g={(0,0,t):1,(1,1,0):-1,(1,1,1):-1}
        assert any(shift(a,(0,1,0)) in D for a in g)
        counts['guard_nonzero_relation_multiples_before_quotient']+=1
        A=generic_model(D,D,truncated_relations(D,[g]),p)
        assert A['n']==4*t
        Av,actual,record=compare_padding(A,(1,1,0),counts,False)
        assert len(actual)==3
        record['model']='nonzero_multiples_of_square_zero_relation'
        record['long_circuit_parameter']=t
        out.append(record)
    return out


def poly_product(f,g):
    out=[0]*(len(f)+len(g)-1)
    for i,c in enumerate(f):
        for j,d in enumerate(g):out[i+j]+=c*d
    return tuple(out)


def poly_power(f,k):
    out=(1,)
    for _ in range(k):out=poly_product(out,f)
    return out


def multivariate_product(f,g):
    out=defaultdict(int)
    for a,c in f.items():
        for b,d in g.items():out[shift(a,b)]+=c*d
    return dict(out)


def fadic_control(source,fs,v,counts):
    n,source_ev,_,p,_,meta=source
    degrees=tuple(len(f)-1 for f in fs)
    D=tuple(product(*(range(d) for d in degrees)))
    base,pivots=field_basis(sorted(D,key=lambda a:(sum(a),a)),source_ev,p)
    assert len(base)==n
    index={a:j for j,a in enumerate(base)}
    @lru_cache(None)
    def old_ev(a):
        rem,co=field_reduce(source_ev(a),pivots,p);assert not rem
        return {index[b]:c for b,c in co.items()}
    original={'D':set(D),'S':set(D),'fixed':set(),'basis':base,
              'n':n,'rho':len(D)-n,'char':p,'ev':old_ev,'block_ev':old_ev,
              'variable':D,'polynomials':[]}
    for i,f in enumerate(fs):
        for a in base:
            replay=defaultdict(int)
            for j,c in enumerate(f):
                for k,d in old_ev(tuple(x+j*int(i==h) for h,x in enumerate(a))).items():
                    replay[k]+=c*d
            assert not normalized(replay,p)
            counts['independent_supplied_univariate_annihilator_columns']+=1
    relations=[]
    for a in D:
        if a in index:continue
        row={a:scalar(1,p)}
        row.update({base[j]:-c for j,c in old_ev(a).items()})
        relations.append(row)
    moduli=tuple(poly_power(f,k+1) for f,k in zip(fs,v))
    factors=tuple(poly_power(f,k) for f,k in zip(fs,v))
    lengths=tuple(len(f)-1 for f in moduli)
    Dv=tuple(product(*(range(d) for d in lengths)))
    positions={a:j for j,a in enumerate(Dv)}
    def tensor(polys):
        out={}
        for terms in product(*(tuple((j,c) for j,c in enumerate(f) if c) for f in polys)):
            a=tuple(j for j,c in terms)
            out[positions[a]]=prod(c for j,c in terms)
        return normalized(out,p)
    @lru_cache(None)
    def image(a):
        return tensor([univariate_product_mod(f,(0,)*b+(1,),m)
                       for f,b,m in zip(factors,a,moduli)])
    new_relations=[]
    for row in relations:
        result=defaultdict(int)
        for a,c in row.items():
            for j,d in image(a).items():result[j]+=c*d
        new_relations.append(normalized(result,p))
    _,new_pivots=field_basis(tuple(range(len(new_relations))),lambda j:new_relations[j],p)
    assert len(new_pivots)==len(D)-n
    reference=tuple(a for j,a in enumerate(Dv) if j not in new_pivots)
    new_index={positions[a]:j for j,a in enumerate(reference)}
    def reduce_vector(raw):
        rem,_=field_reduce(raw,new_pivots,p)
        assert set(rem)<=set(new_index)
        return {new_index[j]:c for j,c in rem.items()}
    @lru_cache(None)
    def ev(a):
        return reduce_vector(tensor([univariate_power_mod((0,1),b,m)
                                     for b,m in zip(a,moduli)]))
    @lru_cache(None)
    def Q_ev(a):
        u,b=zip(*(divmod(x,d) for x,d in zip(a,degrees)))
        return reduce_vector(tensor([univariate_product_mod(poly_power(f,k),(0,)*j+(1,),m)
                                     for f,k,j,m in zip(fs,u,b,moduli)]))
    N=len(reference)
    padded_model={'D':set(Dv),'S':set(Dv),'fixed':set(),'basis':reference,
                  'n':N,'rho':len(new_relations),'char':p,'ev':ev,'block_ev':ev,
                  'variable':Dv,'polynomials':[]}
    assert N==n+len(D)*(prod(k+1 for k in v)-1)
    multiplication_controls(original,counts);multiplication_controls(padded_model,counts)
    before=full_fan(original,counts);after=full_fan(padded_model,counts)
    s=tuple(d*k for d,k in zip(degrees,v))
    halo=set(Dv)-{shift(a,s) for a in D}
    predicted={frozenset(halo|{shift(a,s) for a in B}):poly for B,poly in before.items()}
    compare_fans(after,predicted,counts)
    delta=tuple(a-b-(len(D)-n)*c for a,b,c in zip(profile(Dv),profile(D),s))
    Fpoly={tuple(a for a,c in terms):prod(c for a,c in terms)
           for terms in product(*(tuple((j,c) for j,c in enumerate(f) if c) for f in factors))}
    pure=[{tuple(j if h==i else 0 for h in range(3)):c for j,c in enumerate(f) if c}
          for i,f in enumerate(moduli)]
    old_union=set();lifted_union=set()
    for B,poly in before.items():
        Bv=frozenset(halo|{shift(a,s) for a in B})
        assert profile(Bv)==shift(profile(B),delta)
        order=sorted(Dv,key=lambda a:(dot(center(poly),a),a))
        assert field_basis(order,ev,p)[0]==field_basis(order,Q_ev,p)[0]
        counts['independent_entire_univariate_padding_prefix_scans']+=len(order)
        weight=center(poly)
        _,bp=field_basis(sorted(B),old_ev,p)
        G=[]
        for head in minimal_excluded_border(tuple(B)):
            rem,co=field_reduce(old_ev(head),bp,p);assert not rem
            g={head:scalar(1,p)};g.update({a:-c for a,c in co.items()})
            G.append(normalized(g,p))
        lifted=[normalized(g,p) for g in pure]+[
            normalized(multivariate_product(Fpoly,g),p) for g in G]
        counts['independent_lifted_Groebner_basis_S_pairs']+=pair_checks(lifted,weight,p)
        heads=[max(g,key=lambda a:(dot(weight,a),a)) for g in lifted]
        assert all(any(all(a>=b for a,b in zip(head,h)) for h in heads)
                   for head in minimal_excluded_border(tuple(Bv)))
        old_union.update(tuple(sorted(g.items())) for g in G)
        lifted_union.update(tuple(sorted(g.items())) for g in lifted)
    assert len(lifted_union)<=len(old_union)+3
    C=set.intersection(*(set(B) for B in before))
    Cv=set.intersection(*(set(B) for B in after))
    assert Cv==halo|{shift(a,s) for a in C}
    assert N-len(Cv)==n-len(C)
    edges=edges_from_polygons(before,counts)
    assert {frozenset(halo|{shift(a,s) for a in B}) for B in before}==set(after)
    for B,D0 in map(tuple,edges):
        bv=halo|{shift(a,s) for a in B};dv=halo|{shift(a,s) for a in D0}
        assert tuple(x-y for x,y in zip(profile(bv),profile(dv)))==tuple(
            x-y for x,y in zip(profile(B),profile(D0)))
        counts['independent_univariate_padding_edge_vector_equalities']+=1
    rank_changes=0
    if len(Dv)<=12:
        for pair in combinations(Dv,2):
            rank_changes+=rank([ev(a) for a in pair],p)!=rank([Q_ev(a) for a in pair],p)
    counts['guard_univariate_padding_arbitrary_subset_rank_changes']+=rank_changes
    return {**meta,'field':p,'annihilators':fs,'padding_powers':v,'label_shift':s,
            'original_length':n,'padded_length':N,'monic_box_size':len(D),
            'halo_columns':len(halo),'initial_ideals':len(before),
            'profile_translation':delta,'global_common_core_deficit':n-len(C),
            'original_reduced_basis_union_polynomials':len(old_union),
            'lifted_universal_basis_polynomials':len(lifted_union),
            'arbitrary_subset_rank_changes':rank_changes,
            'scope':'Complete original and padded monomial border fans, evaluated polynomial prefixes, module multiplication and exact profile translation. No equality of arbitrary-subset monomial and polynomial matroids is assumed.'}


def univariate_controls(counts):
    out=[]
    points=((1,2,3),(2,3,5),(3,5,7))
    fs=tuple(poly_product(poly_product((-a[i],1),(-b[i],1)),(-c[i],1))
             for i in range(3) for a,b,c in [(points[0],points[1],points[2])])
    for p in (0,7):
        out.append(fadic_control(point_model(points,p),fs,(1,0,0),counts))
    gaussian_ann=(poly_product((5,-4,1),(1,0,1)),
                  poly_product((10,-6,1),(1,0,1)),(1,0,1))
    out.append(fadic_control(gaussian_model(),gaussian_ann,(0,0,1),counts))
    for p in (0,2):
        one=point_model(((1,2,3),),p)
        out.append(fadic_control(one,((1,-2,1),(-2,1),(-3,1)),(1,0,0),counts))
    assert counts['guard_univariate_padding_arbitrary_subset_rank_changes']>0
    return out


def main():
    counts=Counter()
    generic=arbitrary_relations(counts)
    print('overlapping socle padding:',len(generic),'controls passed',flush=True)
    recross=recross_controls(counts)
    print('padded complete recrossing fans:',len(recross),'passed',flush=True)
    large=symbolic_large_controls(counts)
    propagating=propagating_controls(counts)
    print('propagating monomial padding:',len(propagating),'passed',flush=True)
    univariate=univariate_controls(counts)
    print('all-field univariate padding:',len(univariate),'passed',flush=True)
    sources=['docs/initial_fan_padding.md','scripts/run_initial_fan_padding_checks.py',
             'docs/disjoint_socle_circuits.md','scripts/run_socle_circuit_checks.py',
             'scripts/run_coordinate_slice_budget_checks.py','scripts/run_interleaving_research.py',
             'scripts/run_oracle_enumeration_checks.py','scripts/run_primitive_wall_checks.py',
             'scripts/run_projection_tag_checks.py','scripts/run_quadric_cross_checks.py',
             'scripts/run_signed_face_checks.py','scripts/run_structural_complexity_checks.py',
             'scripts/run_univariate_cross_checks.py','scripts/run_algebraic_incidence_controls.py',
             'scripts/run_initial_cone_height_checks.py','scripts/run_extremal_initial_ideals_checks.py']
    out={'schema_version':1,'all_checks_passed':True,
         'scope':'Exact monomial padding including propagating relations, all-field '
                 'univariate padding for arbitrary finite quotients, and an '
                 'infinite-family obstruction to absolute-divisor-volume reuse charges. '
                 'Intermediate primitive costs and near-full common face cores are retained. '
                 'Independent socle fans remain subgeneric; no unrestricted G3 improvement, '
                 'polynomial lower separation, or novelty established.',
         'counters':dict(counts),'overlapping_relation_space_controls':generic,
         'complete_padded_recrossing_fans':recross,
         'large_formula_and_witness_controls':large,
         'propagating_relation_controls':propagating,
         'univariate_padding_controls':univariate,
         'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    target=ROOT/'results/research/initial_fan_padding_checks.json'
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_checks_passed':True,'counts':dict(counts),
                      'report':str(target.relative_to(ROOT))}),flush=True)


if __name__=='__main__':
    main()
