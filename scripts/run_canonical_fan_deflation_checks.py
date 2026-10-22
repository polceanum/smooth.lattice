#!/usr/bin/env python3
"""Exact order-independent factor extraction, actual colon modules and full fans."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import comb
from pathlib import Path
import hashlib
import json

from run_algebraic_incidence_controls import pair_checks
from run_coordinate_slice_budget_checks import edges_from_polygons, primitive, profile
from run_initial_fan_padding_checks import (full_fan, generic_model, padded,
                                           truncated_relations)
from run_interleaving_research import minimal_excluded_border
from run_primitive_wall_checks import threshold
from run_projection_tag_checks import compare_fans, fan, point_model
from run_quadric_cross_checks import apply, normalized, rank
from run_signed_face_checks import field_basis, field_reduce
from run_socle_circuit_checks import factor_fan, multiplication_controls, recross_model
from run_structural_complexity_checks import (cross, dot, hyperbolic_cross,
                                              univariate_power_mod, univariate_product_mod)
from run_univariate_cross_checks import center

ROOT=Path(__file__).resolve().parents[1]
ZERO=(0,0,0)


def add(a,b):
    return tuple(x+y for x,y in zip(a,b))


def divides(a,b):
    return all(x<=y for x,y in zip(a,b))


def marked_basis(B,ev,p):
    _,pivots=field_basis(tuple(sorted(B)),ev,p)
    G=[]
    for head in minimal_excluded_border(tuple(B)):
        rem,co=field_reduce(ev(head),pivots,p);assert not rem
        g={head:1};g.update({a:-c for a,c in co.items()})
        G.append(normalized(g,p))
    return G,pivots


def common_factor(G):
    support=[a for g in G if len(g)>1 for a in g]
    return tuple(min(a[i] for a in support) for i in range(3)) if support else None


def kernel_factor(n,ev,p,counts):
    E=tuple(a for a in sorted(hyperbolic_cross(n+1,3),key=lambda a:(sum(a),a)) if ev(a))
    B,pivots=field_basis(E,ev,p)
    assert len(B)==n
    support=set();rows=0
    for a in E:
        if a in B:continue
        rem,co=field_reduce(ev(a),pivots,p);assert not rem and co
        support.add(a);support.update(co);rows+=1
        counts['independent_nonzero_column_kernel_relations']+=1
    value=tuple(min(a[i] for a in support) for i in range(3)) if support else None
    counts['complete_universal_finite_kernel_factor_extractions']+=1
    return value,len(E),rows


def divided_basis(G,v):
    H=[]
    for g in G:
        if len(g)==1:
            a,c=next(iter(g.items()))
            H.append({tuple(max(x-y,0) for x,y in zip(a,v)):c})
        else:
            assert all(divides(v,a) for a in g)
            H.append({tuple(x-y for x,y in zip(a,v)):c for a,c in g.items()})
    return H


def evaluate_polynomial(g,ev,p):
    result=defaultdict(int)
    for a,c in g.items():
        for j,d in ev(a).items():result[j]+=c*d
    return normalized(result,p)


def actual_fan(A,counts):
    return full_fan(A,counts) if 'D' in A else fan(A['n'],A['ev'],A['char'],counts)


def control(A,counts):
    n,ev,p=A['n'],A['ev'],A['char']
    actual=actual_fan(A,counts)
    first=next(iter(actual))
    G0,pivots=marked_basis(first,ev,p)
    v=common_factor(G0)
    from_kernel,ground_size,kernel_rows=kernel_factor(n,ev,p,counts)
    assert v==from_kernel
    if v is None:
        assert len(actual)==1 and all(len(g)==1 for g in G0)
        counts['independent_trivial_monomial_cases']+=1
        return {'kind':A['kind'],'field':p,'length':n,'monomial':True,
                'initial_ideals':1,'finite_nonzero_ground':ground_size}
    image_labels=tuple(sorted(a for a in first if divides(v,a)))
    halo=frozenset(set(first)-set(image_labels))
    assert v in image_labels
    m=len(image_labels)
    assert 1<=m<=n and len(halo)==n-m
    assert len(halo)>=product_volume(v)-1>=sum(v)
    _,image_pivots=field_basis(image_labels,ev,p)
    image_index={a:j for j,a in enumerate(image_labels)}
    @lru_cache(None)
    def core_ev(a):
        rem,co=field_reduce(ev(add(v,a)),image_pivots,p)
        assert not rem
        return {image_index[b]:c for b,c in co.items()}
    core_basis=tuple(tuple(x-y for x,y in zip(a,v)) for a in image_labels)
    assert all(core_ev(a)=={j:Q(1)} for j,a in enumerate(core_basis))
    assert rank([ev(add(v,b)) for b in first],p)==m
    counts['independent_actual_multiplication_image_ranks']+=1
    matrices=[]
    for i in range(3):
        columns=[]
        for a,b in zip(core_basis,image_labels):
            step=tuple(int(i==j) for j in range(3))
            rem,co=field_reduce(ev(add(b,step)),pivots,p)
            assert not rem and set(co)<=set(image_labels)
            restricted={image_index[c]:d for c,d in co.items()}
            assert restricted==core_ev(add(a,step))
            columns.append(restricted)
            counts['independent_principal_invariant_multiplication_blocks']+=1
        matrices.append(columns)
    for i,j in combinations(range(3),2):
        for k in range(m):
            assert apply(matrices[i],matrices[j][k],p)==apply(matrices[j],matrices[i][k],p)
            counts['independent_colon_coordinate_commutations']+=1
    if any(v):
        nonideal=False
        for h in halo:
            for i in range(3):
                step=tuple(int(i==j) for j in range(3))
                rem,co=field_reduce(ev(add(h,step)),pivots,p);assert not rem
                nonideal|=bool(set(co)-set(halo))
        assert nonideal
        counts['guard_halo_span_is_not_the_quotient_kernel']+=1
    core=fan(m,core_ev,p,counts)
    predicted={frozenset(set(halo)|{add(v,a) for a in B}):poly for B,poly in core.items()}
    compare_fans(actual,predicted,counts)
    delta=add(profile(halo),tuple(m*x for x in v))
    marked=[]
    for B,poly in actual.items():
        G,_=marked_basis(B,ev,p)
        assert common_factor(G)==v
        w=center(poly)
        counts['independent_original_reduced_basis_S_pairs']+=pair_checks(G,w,p)
        H=divided_basis(G,v)
        assert len(H)<=len(G)
        for h in H:assert not evaluate_polynomial(h,core_ev,p)
        counts['independent_divided_colon_basis_S_pairs']+=pair_checks(H,w,p)
        Bcore=frozenset(tuple(x-y for x,y in zip(a,v)) for a in B if divides(v,a))
        assert Bcore in core and set(B)-{add(v,a) for a in Bcore}==set(halo)
        assert profile(B)==add(profile(Bcore),delta)
        heads=[max(h,key=lambda a:(dot(w,a),a)) for h in H]
        assert all(any(divides(a,b) for a in heads) for b in minimal_excluded_border(tuple(Bcore)))
        Gcore,_=marked_basis(Bcore,core_ev,p)
        assert common_factor(Gcore)==ZERO
        counts['independent_every_order_factor_and_idempotence_checks']+=1
        marked.append({'original_basis_polynomials':len(G),'divided_basis_polynomials':len(H)})
    again,_,_=kernel_factor(m,core_ev,p,counts)
    assert again==ZERO
    edges=edges_from_polygons(actual,counts);core_edges=edges_from_polygons(core,counts)
    assert {frozenset(set(halo)|{add(v,a) for a in B}) for B in core}==set(actual)
    assert len(edges)==len(core_edges)
    costs=[]
    for B,C in map(tuple,edges):
        b=frozenset(tuple(x-y for x,y in zip(a,v)) for a in B if divides(v,a))
        c=frozenset(tuple(x-y for x,y in zip(a,v)) for a in C if divides(v,a))
        assert frozenset((b,c)) in core_edges
        difference=tuple(x-y for x,y in zip(profile(B),profile(C)))
        assert difference==tuple(x-y for x,y in zip(profile(b),profile(c)))
        nu=threshold(primitive(difference));assert nu<=m
        costs.append(nu)
        counts['independent_edge_vector_and_effective_cost_checks']+=1
    if A.get('expected_factor') is not None:assert v==A['expected_factor']
    if A.get('expected_core_length') is not None:assert m==A['expected_core_length']
    return {'kind':A['kind'],'field':p,'length':n,'factor_exponents':v,'core_length':m,
            'removed_halo_labels':sorted(halo),'profile_translation':delta,
            'initial_ideals':len(actual),'undirected_walls':len(edges),
            'largest_primitive_cost':max(costs,default=0),'total_primitive_cost':sum(costs),
            'finite_nonzero_ground':ground_size,'finite_kernel_relations':kernel_rows,
            'marked_basis_checks':marked,
            'scope':'The actual colon is evaluated through the invariant multiplication image. Its whole fan is traversed separately and compared cone by cone.'}


def product_volume(a):
    out=1
    for x in a:out*=x+1
    return out


def univariate_source(modulus,images,p,kind):
    n=len(modulus)-1
    @lru_cache(None)
    def ev(a):
        result=(1,)+(0,)*(n-1)
        for image,k in zip(images,a):
            result=univariate_product_mod(result,univariate_power_mod(image,k,modulus),modulus)
        return normalized(dict(enumerate(result)),p)
    return {'n':n,'ev':ev,'char':p,'kind':kind}


def nonlocal_halo(a,p):
    source=univariate_source((0,0,0,0,1),((1,1),(0,0,1),(0,)),p,'core')
    D=tuple(product(range(4),range(2),(0,)))
    base,pivots=field_basis(D,source['ev'],p);assert len(base)==4
    index={b:j for j,b in enumerate(base)}
    halo=tuple(product(range(a),range(2),(0,)))
    halo_index={b:j for j,b in enumerate(halo)}
    @lru_cache(None)
    def ev(b):
        if b[2] or b[1]>=2:return {}
        if b[0]<a:return {halo_index[b]:1}
        c=(b[0]-a,b[1],0)
        rem,co=field_reduce(source['ev'](c),pivots,p);assert not rem
        return {len(halo)+index[d]:x for d,x in co.items()}
    return {'n':2*a+4,'ev':ev,'char':p,'kind':'nonlocal_propagating_halo',
            'expected_factor':(a,0,0),'expected_core_length':4}


def fixtures(counts):
    out=[]
    D={(0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,2,0)}
    for p in (0,2,3):
        A=generic_model(D,{(1,1,0),(0,2,0)},[{(1,1,0):1,(0,2,0):-1}],p)
        A.update(kind='nonideal_halo_example',expected_factor=(0,1,0),expected_core_length=2)
        out.append(A)
    for a,p in product((1,3,7),(0,2,3)):out.append(nonlocal_halo(a,p))
    for a,p in product((1,4),(0,2,3)):
        # x^a*(x-1)*(x-2), with field-dependent zero root in char 2.
        modulus=(0,)*a+(2,-3,1)
        A=univariate_source(modulus,((0,1),(0,),(0,)),p,'nonlocal_univariate_factor')
        A.update(expected_factor=(a+int(p==2),0,0),expected_core_length=1 if p==2 else 2)
        out.append(A)
    for t,p in product((2,3),(0,2,3)):
        D=set(product(range(2),range(2),range(2*t)))
        g={(0,0,t):1,(1,1,0):-1,(1,1,1):-1}
        A=generic_model(D,D,truncated_relations(D,[g]),p)
        Av=padded(A,(1,2,0),counts)
        Av.update(kind='padded_propagating_relation',expected_factor=(1,2,0),
                  expected_core_length=A['n'])
        out.append(Av)
    for s,p in product((1,2),(0,2)):
        old,_,_=recross_model(s,p)
        A=generic_model(old['D'],old['S'],old['polynomials'],p)
        Av=padded(A,(0,s,0),counts)
        Av.update(kind='padded_repeated_upper_pair',expected_factor=(0,s,0),
                  expected_core_length=A['n'])
        out.append(Av)
    for p in (0,7):
        n,ev=point_model(((1,2,3),(2,3,5),(3,5,7)),p)
        out.append({'n':n,'ev':ev,'char':p,'kind':'primitive_radical_points',
                    'expected_factor':ZERO,'expected_core_length':n})
    for p in (0,2):
        out.append(dict(univariate_source((-1,3,-3,1),((0,1),(2,),(3,)),p,
                                         'univariate_padding_not_removed'),
                        expected_factor=ZERO,expected_core_length=3))
    mono=generic_model(set(product(range(2),range(2),(0,))),set(),[],0)
    mono['kind']='trivial_monomial_case';out.append(mono)
    return out


def failure_guards(counts):
    D={(0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,2,0)}
    A=generic_model(D,{(1,1,0),(0,2,0)},[{(1,1,0):1,(0,2,0):-1}],0)
    cells=full_fan(A,counts);B=next(iter(cells));G,_=marked_basis(B,A['ev'],0)
    assert common_factor(G)==(0,1,0)
    redundant={(2,0,0):1,(0,1,1):1}
    assert not evaluate_polynomial(redundant,A['ev'],0)
    assert common_factor(G+[redundant])==ZERO
    counts['guard_redundant_nonreduced_input_changes_factor']+=1
    # Arbitrary colon by x changes the fan of x=y^2, y^3=0.
    source=univariate_source((0,0,0,1),((0,0,1),(0,1),(0,)),0,'bad_arbitrary_colon')
    cells=fan(3,source['ev'],0,counts)
    G,_=marked_basis(next(iter(cells)),source['ev'],0)
    assert common_factor(G)==ZERO and len(cells)==2
    assert rank([source['ev'](add((1,0,0),a)) for a in next(iter(cells))],0)==1
    counts['guard_arbitrary_colon_changes_the_fan']+=1
    return {'redundant_input_factor':[0,0,0],'correct_reduced_factor':[0,1,0],
            'arbitrary_colon_original_cones':2,'arbitrary_colon_core_length':1,
            'scope':'No arbitrary common-subspace or arbitrary colon replacement is inferred.'}


def primitive_reuse_guard(s,A,p,counts):
    """One extra disjoint circuit defeats the global core-size repair."""
    old,_,_=recross_model(s,p)
    original=generic_model(old['D'],old['S'],old['polynomials'],p)
    base=padded(original,(0,A,0),counts)
    B=4*s*s
    a,b=(2,0,0),(0,0,B+2)
    assert a not in base['D'] and b not in base['D']
    D=base['D']|{a,b}
    rows=base['polynomials']+[{a:1,b:-1}]
    model=generic_model(D,base['S']|{a,b},rows,p)
    n=model['n']
    expected=(32*s**3+30*s*s+10*s+12)//3+(8*s*s+4)*A+1
    assert n==expected==base['n']+1
    multiplication_controls(model,counts)
    cells=full_fan(model,counts)
    compare_fans(cells,factor_fan(model),counts)
    assert 4*s<=len(cells)<=8*s
    assert kernel_factor(n,model['ev'],p,counts)[0]==ZERO
    for basis,poly in cells.items():
        G,_=marked_basis(basis,model['ev'],p)
        assert common_factor(G)==ZERO
        counts['primitive_guard_full_reduced_basis_S_pairs']+=pair_checks(G,center(poly),p)
        counts['primitive_guard_every_marked_factor_zero']+=1
    edges=edges_from_polygons(cells,counts)
    junctions=[]
    for i in range(1,2*s):
        t=i*(i+1);u=Q(2*i+1,2);w=(t,u,1)
        assert t!=2*s*s+1
        incident=[basis for basis,poly in cells.items()
                  if all(dot(cross(a,b),w)>=0 for a,b in zip(poly,poly[1:]+poly[:1]))]
        assert len(incident)==4
        common=set.intersection(*map(set,incident))
        assert len(common)==n-2
        upper={(1,A+1,B),(1,A,B+1)}
        lower={(0,A+2*j+int(j%2==0),B-j*j+int(j%2==1)) for j in (i,i+1)}
        assert set.union(*map(set,incident))-common==upper|lower
        local_edges=[edge for edge in edges if set(edge)<=set(incident)]
        assert len(local_edges)==4
        directions=set()
        for C,E in map(tuple,local_edges):
            q=primitive(tuple(x-y for x,y in zip(profile(C),profile(E))))
            if q[0]<0:q=tuple(-x for x in q)
            directions.add(q)
        expected_directions={(1,-2*j,j*j) for j in (i,i+1)}
        assert directions==expected_directions
        costs=sorted(threshold(q) for q in directions)
        assert costs==[2*i*i+2*i+1,2*(i+1)**2+2*(i+1)+1]
        junctions.append({'i':i,'weight':[t,str(u),1],'incident_cells':4,
                          'actual_common_core_size':len(common),
                          'primitive_directions':sorted(directions),'primitive_costs':costs})
        counts['primitive_guard_preserved_persistent_four_cell_junctions']+=1
    # At any nonzero center one translated coordinate is a unit. An
    # essential annihilator containing 1 certifies monomial gcd one.
    nilpotence=[max(a[j] for a in D)+1 for j in range(3)]
    centers=(list(product(range(p),repeat=3)) if p else
             [ZERO,(1,0,0),(0,1,0),(0,0,1),(1,-1,2)])
    for c in centers:
        if c==ZERO:continue
        j=next(j for j in range(3) if c[j])
        N=nilpotence[j];value=c[j]
        polynomial=normalized({k:comb(N,k)*value**(N-k) for k in range(N+1)},p)
        assert polynomial[0]
        def ev_power(k):
            out=defaultdict(int)
            for l in range(k+1):
                e=tuple(l*int(j==h) for h in range(3))
                for z,v in model['ev'](e).items():
                    out[z]+=comb(k,l)*(-value)**(k-l)*v
            return normalized(out,p)
        assert all(ev_power(k) for k in polynomial)
        total=defaultdict(int)
        for k,v in polynomial.items():
            for z,vv in ev_power(k).items():total[z]+=v*vv
        assert not normalized(total,p)
        counts['primitive_guard_nonzero_center_essential_constant_certificates']+=1
    return {'s':s,'padding':A,'field':p,'length':n,
            'core_length_at_every_center':n,'initial_ideals':len(cells),
            'persistent_crossings':len(junctions),
            'upper_divisor_volumes':[2*(A+2)*(B+1),2*(A+1)*(B+2)],
            'junctions':junctions,'all_marked_factors_one':True,
            'translated_centers_checked':len(centers)-1,
            'scope':'Every-center core length follows uniformly from the local unit argument. The finite translated annihilators and independent full-fan comparisons support that proof; they do not enumerate centers over infinite fields.'}


def primitive_reuse_arithmetic(counts):
    rows=[]
    for s,A in [(s,s) for s in (1024,4096,16384)]+[(k**4,k**5) for k in (64,128,256)]:
        n=(32*s**3+30*s*s+10*s+12)//3+(8*s*s+4)*A+1
        volumes=[2*(A+2)*(4*s*s+1),2*(A+1)*(4*s*s+2)]
        assert max(volumes)<=n
        costs=[2*i*i+2*i+1 for i in (s,(3*s)//2,2*s-1)]
        assert all(cost**5>n**3 and (n-cost)**7>=n**3 for cost in costs)
        counts['primitive_guard_exact_middle_cost_scale_controls']+=len(costs)
        rows.append({'s':s,'padding':A,'length':n,'persistent_crossings':2*s-1,
                     'upper_divisor_volumes':volumes,'sampled_bulk_costs':costs,
                     'scope':'Exact integer formula checks; no huge quotient fan is enumerated.'})
    return rows


def main():
    counts=Counter();rows=[]
    for A in fixtures(counts):
        row=control(A,counts);rows.append(row)
        print('canonical factor:',row['kind'],row['field'],row['length'],
              row.get('factor_exponents'),'core',row.get('core_length'),
              'g',row['initial_ideals'],flush=True)
    guards=failure_guards(counts)
    primitive_guards=[primitive_reuse_guard(s,A,p,counts)
                      for s,A,p in ((1,1,0),(1,1,2),(2,2,0),(2,2,2),(2,5,3))]
    primitive_arithmetic=primitive_reuse_arithmetic(counts)
    sources=['docs/canonical_fan_deflation.md','scripts/run_canonical_fan_deflation_checks.py',
             'docs/initial_fan_padding.md','scripts/run_initial_fan_padding_checks.py',
             'docs/rational_shadow_bound.md','docs/wall_coexistence.md',
             'scripts/run_algebraic_incidence_controls.py','scripts/run_coordinate_slice_budget_checks.py',
             'scripts/run_interleaving_research.py','scripts/run_primitive_wall_checks.py',
             'scripts/run_projection_tag_checks.py','scripts/run_quadric_cross_checks.py',
             'scripts/run_signed_face_checks.py','scripts/run_socle_circuit_checks.py',
             'scripts/run_structural_complexity_checks.py','scripts/run_univariate_cross_checks.py']
    out={'schema_version':1,'all_checks_passed':True,
         'scope':'Canonical monomial factor independent of global order, matrix-only extraction, '
                 'actual colon algebra and invariant multiplication block, complete fan preservation, '
                 'and effective length in existing wall/count bounds. This does not remove all '
                 'univariate padding or improve unrestricted G3. A one-relation extension retains '
                 'padded same-pair recrossing while making the core full length at every translation. '
                 'Exact novelty is unassessed.',
         'counters':dict(counts),'complete_quotient_controls':rows,'failure_guards':guards,
         'primitive_core_reuse_guards':primitive_guards,
         'primitive_core_middle_cost_arithmetic':primitive_arithmetic,
         'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    target=ROOT/'results/research/canonical_fan_deflation_checks.json'
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_checks_passed':True,'counters':dict(counts),
                      'report':str(target.relative_to(ROOT))}),flush=True)


if __name__=='__main__':
    main()
