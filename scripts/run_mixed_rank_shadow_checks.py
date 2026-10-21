#!/usr/bin/env python3
"""Exact compulsory-carrier contraction, mixed covers, shadows and full fans."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd, prod
from pathlib import Path
from random import Random
import hashlib
import json

from run_coordinate_slice_budget_checks import downset, edges_from_polygons, primitive, profile
from run_mixed_ideal_rank_checks import (long_gluing, mixed, polynomial_columns,
                                       shifted_rectangle, tensor_model)
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_planar_dyadic_bound_checks import convex_hull, matrix_oracle
from run_primitive_wall_checks import threshold
from run_projection_tag_checks import compare_fans, fan
from run_quadric_cross_checks import apply, rank
from run_rational_shadow_bound_checks import (cover, direction_control, kernel_projection,
                                              project, projection_data, shared_lines)
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross
from run_univariate_cross_checks import TRIANGLE, center, extended_sharp, gaussian_model, point_model

ROOT=Path(__file__).resolve().parents[1]
AXES=[((0,0,0),tuple(int(i==j) for i in range(3))) for j in range(3)]


def meeting_axes(L,M,p):
    n=L+M
    @lru_cache(None)
    def ev(a):
        x,y,z=a
        if z or (x and y):return {}
        if not x and not y:return {0:1}
        if x:return {x:1} if x<=L else {}
        if y<M:return {L+y:1}
        return {L:1} if y==M else {}
    reference=((0,0,0),)+tuple((j,0,0) for j in range(1,L+1))+tuple(
        (0,j,0) for j in range(1,M))
    matrices=[tuple(ev(tuple(v+int(i==j) for j,v in enumerate(a))) for a in reference)
              for i in range(3)]
    assert tuple(ev(a) for a in reference)==tuple({j:1} for j in range(n))
    return n,ev,matrices,p,((0,1),)*3,{
        'kind':'long_axes_meeting_in_socle','field':p,'x_length':L,'y_length':M}


def hybrid_sweep(labels,r,oracle,H,A,signs,counts):
    points=tuple(project(A,a) for a in labels)
    # Choose strip widths against ALL ground points, including body points.
    strips=[]
    for start,direction in AXES:
        c,v=project(A,start),project(A,direction)
        if any(v):
            normal=(-v[1],v[0]);level=dot(normal,c)
            distances=[abs(dot(normal,x)-level) for x in points if dot(normal,x)!=level]
            eps=min(distances)/Q(2) if distances else Q(1)
            lines=[(normal,level-eps),(normal,level+eps)]
        else:
            distances=[abs(x[i]-c[i]) for x in points for i in range(2) if x[i]!=c[i]]
            eps=min(distances)/Q(2) if distances else Q(1)
            lines=[(tuple(int(i==j) for i in range(2)),c[j]+s*eps)
                   for j in range(2) for s in (-1,1)]
        strips.append({'point':c,'direction':v,'lines':lines})
    def contains_strip(c,v):
        return all(dot(row,v)>low and dot(row,v)<high
                   for (row,low),(_,high) in zip(c['lines'][::2],c['lines'][1::2]))
    boxes=cover(A,max(1,H))
    boundaries=list(set(shared_lines(boxes))|{line for c in strips for line in c['lines']})
    assert all(dot(row,v)!=level for row,level in boundaries for v in points)
    assert len(boundaries)<=12*max(1,H).bit_length()+12
    assigned=[]
    for a,v in zip(labels,points):
        j=next((j for j,c in enumerate(strips) if contains_strip(c,v)),None)
        if j is not None:
            c=strips[j];step=tuple(x-y for x,y in zip(v,c['point']))
            assert not any(step) if not any(c['direction']) else (
                step[0]*c['direction'][1]==step[1]*c['direction'][0])
            assigned.append(('axis',j))
        else:
            assert prod(x+1 for x in a)<=H
            j=next(j for j,b in enumerate(boxes) if all(
                2*abs(x)<width for x,width in zip(project(b['T'],v),b['widths'])))
            assigned.append(('body',j))
    sigma,tau=signs
    def greedy(t):
        B=set()
        for j in sorted(range(len(labels)),key=lambda j:(sigma*t*points[j][0]+tau*points[j][1],j)):
            if oracle(frozenset(B|{j})):B.add(j)
        assert len(B)==r
        return B
    events=defaultdict(set)
    for i,j in combinations(range(len(labels)),2):
        dx,dy=(a-b for a,b in zip(points[i],points[j]))
        if sigma*tau*dx*dy<0:events[-Q(tau*dy,sigma*dx)].update((i,j))
    events=sorted(events.items())
    first=events[0][0]/2 if events else Q(1)
    B=greedy(first)
    order=sorted(range(len(labels)),key=lambda j:(tau*points[j][1],sigma*points[j][0],j))
    positions={j:i for i,j in enumerate(order)}
    paths=[[points[j]] for j in sorted(B)]
    tokens={j:i for i,j in enumerate(sorted(B))}
    local=defaultdict(lambda:defaultdict(lambda:[0,0]))
    variation=defaultdict(lambda:[0,0]);line_events=defaultdict(set)
    changing=set();states=set();cross_moves=0
    def replay(t):
        assert B==greedy(t)
        states.add(tuple(sum(points[j][i] for j in B) for i in range(2)))
        counts['independent_every_interval_greedy_replays']+=1
    replay(first)
    for k,(t,involved) in enumerate(events):
        groups=defaultdict(list)
        for j in involved:groups[sigma*t*points[j][0]+tau*points[j][1]].append(j)
        for group in groups.values():
            slots=sorted(positions[j] for j in group)
            assert slots==list(range(slots[0],slots[-1]+1))
            for target,f in enumerate(sorted(group,key=lambda j:(sigma*points[j][0],j)),slots[0]):
                while positions[f]>target:
                    pos=positions[f];e=order[pos-1]
                    if e in B and f not in B and oracle(frozenset(B-{e}|{f})):
                        step=tuple(b-a for a,b in zip(points[e],points[f]))
                        assert -sigma*step[0]>0 and tau*step[1]>0
                        assert Q(tau*step[1],-sigma*step[0])==t
                        token=tokens.pop(e);tokens[f]=token;paths[token].append(points[f])
                        left,right=assigned[e],assigned[f]
                        if left!=right:
                            assert any((dot(row,points[e])<level)!=(dot(row,points[f])<level)
                                       for row,level in boundaries)
                            cross_moves+=1
                        elif left[0]=='axis':
                            v=strips[left[1]]['direction']
                            assert any(v) and step[0]*v[1]==step[1]*v[0]
                            line_events[left[1]].add(t)
                        else:
                            image=project(boxes[left[1]]['T'],step)
                            for i in range(2):
                                local[left[1]][t][i]+=image[i]
                                variation[token,left[1]][i]+=abs(image[i])
                        changing.add(t);B.remove(e);B.add(f)
                    order[pos-1],order[pos]=f,e;positions[f],positions[e]=pos-1,pos
        following=events[k+1][0] if k+1<len(events) else t+2
        replay((t+following)/2)
    crossings=0
    for path in paths:
        for row,level in boundaries:
            vals=[dot(row,x)-level for x in path]
            assert all(vals)
            moves=[b-a for a,b in zip(vals,vals[1:]) if b!=a]
            assert sum(a*b<0 for a,b in zip(moves,moves[1:]))<=1
            crossed=sum((a<0)!=(b<0) for a,b in zip(vals,vals[1:]))
            assert crossed<=2;crossings+=crossed
            counts['independent_combined_cover_token_boundary_checks']+=1
    assert cross_moves<=crossings<=2*r*len(boundaries)
    assert all(len(v)<=1 for v in line_events.values())
    delta,height=projection_data(A);internal=[]
    for j,rows in sorted(local.items()):
        b=boxes[j]
        for token in range(r):
            assert all(x<=2*w for x,w in zip(variation[token,j],b['widths']))
        vectors=[v for _,v in sorted(rows.items())]
        ctl=direction_control(vectors,delta*b['minor'])
        assert all(x<=2*r*w for x,w in zip(ctl['absolute_budgets'],b['widths']))
        internal.append({'box':j,'directions':len(vectors),
                         'absolute_budgets':ctl['absolute_budgets'],
                         'covolume':delta*b['minor']})
        counts['independent_body_lattice_direction_and_return_budgets']+=1
    assert len(changing)<=cross_moves+sum(map(len,line_events.values()))+sum(map(len,local.values()))
    counts['complete_mixed_ground_quadrant_sweeps']+=1
    return states,{'signs':signs,'changing_parameters':len(changing),
                   'shadow_profiles':len(states),'cover_boundary_lines':len(boundaries),
                   'cross_label_exchanges':cross_moves,'boundary_crossings':crossings,
                   'internal_axis_events':sum(map(len,line_events.values())),
                   'internal_body_controls':internal,'normalized_projection_height':height}


def full_shadow(labels,r,oracle,H,q,counts,exhaustive=False):
    A=kernel_projection(q);vertices=set();traces=[]
    for signs in product((-1,1),repeat=2):
        values,trace=hybrid_sweep(labels,r,oracle,H,A,signs,counts)
        vertices.update(values);traces.append(trace)
    if exhaustive:
        values={tuple(sum(project(A,labels[j])[i] for j in B) for i in range(2))
                for B in combinations(range(len(labels)),r) if oracle(frozenset(B))}
        assert vertices==set(convex_hull(values))
        counts['independent_exhaustive_entire_shadow_hulls']+=1
    counts['complete_full_mixed_shadows']+=1
    return vertices,traces


def matroid_fan(labels,r,oracle,P,counts):
    comparisons={primitive(tuple(x-y for x,y in zip(a,b))) for a,b in combinations(labels,2)}
    comparisons={q for q in comparisons if min(q)<0<max(q)}
    queue=[center(TRIANGLE,comparisons)];out={}
    bound=1+max((max(a) for a in labels),default=1)
    while queue:
        w=queue.pop();B=set()
        for j in sorted(range(len(labels)),key=lambda j:(dot(w,labels[j]),labels[j])):
            if oracle(frozenset(B|{j})):B.add(j)
        assert len(B)==r
        original=frozenset(set(P)|{labels[j] for j in B})
        if original in out:continue
        normals={tuple(x-y for x,y in zip(labels[j],labels[i]))
                 for i in B for j in range(len(labels)) if j not in B
                 and oracle(frozenset(B-{i}|{j}))}
        assert all(dot(q,w)>0 for q in normals)
        poly=list(TRIANGLE)
        for q in normals:poly=clip(poly,q,1)
        assert area2(poly)>0
        out[original]=poly
        for a,b in zip(poly,poly[1:]+poly[:1]):
            if any(a[i]==b[i]==0 for i in range(3)):continue
            q=next(q for q in normals if dot(q,a)==dot(q,b)==0)
            w,_=exact_adjacent_weight(bound,a,b,q);queue.append(tuple(map(Q,w)))
    assert sum(area2(p) for p in out.values())==1
    counts['complete_contracted_matroid_fans']+=1
    return out


def actual_control(model,counts):
    n,mono,M,p,fs,meta=model
    degrees,F,_,ev,operator_ev=polynomial_columns(model)
    assert degrees==(1,1,1)
    Tcols=[apply(F[i],v,p) for i,j in combinations(range(3),2) for v in F[j]]
    t=rank(Tcols,p)
    for A in M:
        assert rank(Tcols+[apply(A,v,p) for v in Tcols],p)==t
        counts['independent_actual_mixed_ideal_invariance']+=1
    for i,j in combinations(range(3),2):
        for k in range(n):
            assert apply(M[i],M[j][k],p)==apply(M[j],M[i][k],p)
            counts['independent_commuting_matrix_columns']+=1
    body={tuple(a+b for a,b in zip(h,v)) for h in ((1,1,0),(1,0,1),(0,1,1))
          for v in hyperbolic_cross(t,3)} if t else set()
    assert all(prod(x+1 for x in a)<=4*t for a in body)
    body=tuple(sorted(a for a in body if ev(a)))
    groups=[[] for _ in range(4)]
    for a in sorted(hyperbolic_cross(n,3),key=lambda a:(sum(a),a)):
        if mixed(a,degrees):continue
        j=next((i+1 for i,x in enumerate(a) if x),0)
        groups[j].append(a)
    C_lists=[field_basis(group,ev,p)[0] for group in groups]
    C=tuple(a for group in C_lists for a in group)
    ground=C+body;m_C=len(C);b0=rank([ev(a) for a in body],p)
    for a in ground:
        assert ev(a)==operator_ev(a)
        counts['independent_polynomial_expansion_operator_replays']+=1
        for i,A in enumerate(M):
            successor=tuple(v+int(i==j) for j,v in enumerate(a))
            assert apply(A,mono(a),p)==mono(successor)
            counts['independent_raw_monomial_successor_replays']+=1
    assert rank([ev(a) for a in ground],p)==n and b0<=t
    assert n-t<=m_C<=n+3+3*t
    h=m_C-n+b0
    assert 0<=h<=3+4*t
    @lru_cache(None)
    def dual_independent(S):
        return rank([ev(a) for a in ground if a not in S],p)==n
    dual_C=set()
    for a in C:
        if dual_independent(frozenset(dual_C|{a})):dual_C.add(a)
    assert len(dual_C)==h
    counts['independent_exact_dual_carrier_restriction_ranks']+=1
    Fkeep=set()
    for group in C_lists:
        B=set()
        for a in reversed(group):
            if dual_independent(frozenset(B|{a})):B.add(a)
        assert len(B)<=h;Fkeep.update(B)
    P=frozenset(set(C)-Fkeep);labels=tuple(sorted(Fkeep|set(body)))
    r=n-len(P)
    assert len(Fkeep)<=4*h and r==b0-h+len(Fkeep)<=9+13*t
    assert rank([ev(a) for a in P],p)==len(P)
    selected,pivots=field_basis(tuple(sorted(P)),ev,p)
    assert len(selected)==len(P)
    residual=[field_reduce(ev(a),pivots,p)[0] for a in labels]
    assert rank(residual,p)==r
    @lru_cache(None)
    def oracle(S):
        return len(S)<=r and rank([residual[j] for j in S],p)==len(S)
    rng=Random(1701+n+p)
    for _ in range(32):
        S=frozenset(j for j in range(len(labels)) if rng.randrange(3)==0)
        assert oracle(S)==(rank([ev(a) for a in P]+[ev(labels[j]) for j in S],p)==len(P)+len(S))
        counts['independent_contraction_rank_replays']+=1
    actual=fan(n,mono,p,counts)
    assert all(P<=B<=set(ground) for B in actual)
    predicted=matroid_fan(labels,r,oracle,P,counts)
    compare_fans(actual,predicted,counts)
    counts['independent_fixed_compulsory_carrier_sets']+=1
    nonideal=any(rank([ev(a) for a in P]+[apply(A,ev(a),p) for a in P],p)>len(P) for A in M)
    counts['guard_nonideal_forced_carrier_spans']+=int(nonideal)
    outside=0
    for w in product((-1,1),repeat=3):
        B=set()
        objective=tuple(v*(i+2) for i,v in enumerate(w))
        for j in sorted(range(len(labels)),key=lambda j:(dot(objective,labels[j]),j)):
            if oracle(frozenset(B|{j})):B.add(j)
        assert len(B)==r
        if not downset(set(P)|{labels[j] for j in B}):outside+=1
    counts['guard_other_objective_nondownset_bases']+=outside
    directions=defaultdict(set)
    for a,b in combinations(labels,2):
        q=primitive(tuple(x-y for x,y in zip(a,b)))
        if not min(q)<0<max(q):continue
        kind='body_body' if a in body and b in body else (
            'carrier_carrier' if a in Fkeep and b in Fkeep else 'carrier_body')
        directions[kind].add(q)
        if kind=='body_body':
            assert threshold(q)<=8*t-2
            counts['independent_body_direction_threshold_checks']+=1
    all_directions=set().union(*directions.values()) if directions else set()
    assert len(directions['carrier_carrier'])<=comb(len(Fkeep),2)
    assert len(directions['carrier_body'])<=len(Fkeep)*len(body)
    edges=edges_from_polygons(actual,counts)
    by_direction=defaultdict(list)
    for B,D in map(tuple,edges):
        q=primitive(tuple(x-y for x,y in zip(profile(B),profile(D))))
        assert q in all_directions
        by_direction[q].append((B,D))
    assert len(actual)<=1+len(edges)
    qs={(1,0,0),(0,0,1),(1,-2,3),(1000001,-1000000,1)}
    qs.update(by_direction)
    shadow_rows=[]
    for q in sorted(qs):
        vertices,traces=full_shadow(labels,r,oracle,4*t,q,counts,
                                   len(labels)<=12 and comb(len(labels),r)<=3000)
        A=kernel_projection(q)
        fibers={project(A,tuple(x-y for x,y in zip(profile(B),profile(P))))
                for B,D in by_direction[q]}
        assert len(fibers)==len(by_direction[q]) and fibers<=vertices
        strict=[B for B,poly in actual.items() if min(dot(q,v) for v in poly)<0<max(dot(q,v) for v in poly)]
        strict_images=[project(A,tuple(x-y for x,y in zip(profile(B),profile(P)))) for B in strict]
        assert len(set(strict_images))==len(strict_images) and set(strict_images)<=vertices
        counts['independent_parallel_fiber_and_strict_cut_injections']+=1
        shadow_rows.append({'kernel':q,'full_shadow_vertices':len(vertices),
                            'parallel_walls':len(fibers),'strictly_cut_cones':len(strict),
                            'quadrants':traces})
    counts['guard_axis_labels_above_mixed_height']+=int(any(
        prod(x+1 for x in a)>max(1,4*t) for a in Fkeep))
    if meta['kind']=='long_axes_meeting_in_socle':
        assert t==0 and r==1 and len(actual)==2 and len(edges)==1
        assert set(Fkeep)=={(meta['x_length'],0,0),(0,meta['y_length'],0)}
        divisor=gcd(meta['x_length'],meta['y_length'])
        assert sum(threshold(q)*len(v) for q,v in by_direction.items())==n//divisor
        counts['guard_primitive_cost_gcd_reduction']+=int(divisor>1)
        counts['independent_fixed_rank_unbounded_label_and_cost_guards']+=1
    return {**meta,'field':p,'length':n,'linear_factors':fs,'mixed_ideal_rank':t,
            'carrier_columns':m_C,'body_columns':len(body),'body_rank':b0,
            'dual_carrier_restriction_rank':h,'retained_carrier_columns':len(Fkeep),
            'fixed_compulsory_labels':sorted(P),'contracted_ground_labels':labels,
            'contracted_rank':r,'initial_ideals':len(actual),'undirected_walls':len(edges),
            'nonideal_compulsory_span':nonideal,'other_objective_nondownset_bases':outside,
            'direction_inventory':{k:len(v) for k,v in directions.items()},
            'shadows':shadow_rows}


def synthetic_controls(counts):
    rows=[]
    for H,r,length in ((1,0,17),(1,2,17),(2,3,1000001),(4,3,101)):
        body={tuple(a+b for a,b in zip(h,v)) for h in ((1,1,0),(1,0,1),(0,1,1))
              for v in hyperbolic_cross(H,3)}
        labels=tuple(sorted(body|{tuple(k*length if i==j else 0 for i in range(3))
                                 for j in range(3) for k in (1,2)}))
        rng=Random(H*11)
        columns=[tuple(int(i==j) for i in range(r)) if j<r else
                 tuple(rng.randrange(101) for _ in range(r)) for j in range(len(labels))]
        for kind,oracle in (('uniform',lambda S,r=r:len(S)<=r),
                            ('linear_F101',matrix_oracle(columns,r,101))):
            for q in ((1,0,0),(1,-2,3),(1000001,-1000000,1)):
                vertices,traces=full_shadow(labels,r,oracle,4*H,q,counts,
                                           len(labels)<=12)
                rows.append({'kind':kind,'body_parameter':H,'rank':r,
                             'axis_spacing':length,'kernel':q,'ground_size':len(labels),
                             'shadow_vertices':len(vertices),'quadrants':traces})
    return rows


def direction_controls(counts):
    rows=[]
    for H in (1,2,4,8,16,32,64,128):
        total=sum(max(0,H+1-u*v) for u in range(1,H+2) for v in range(1,(H+1)//u+1))
        divisor_sum=sum((H+1)//u for u in range(1,H+2))
        assert total<=(H+1)*divisor_sum
        counts['independent_mixed_direction_inventory_sums']+=1
        rows.append({'threshold':H,'oriented_overcount':total,'divisor_sum':divisor_sum})
    return rows


def main():
    counts=Counter();models=[]
    models.extend(meeting_axes(L,M,p) for L,M,p in ((3,4,0),(64,9,0),(48,9,2)))
    models.extend(long_gluing(s,p) for s,p in ((3,0),(8,0),(64,0),(4,2)))
    models.extend(tensor_model(2,2,p,True) for p in (0,2,101))
    models.extend(tensor_model(3,2,p,False) for p in (0,3))
    for p in (0,7):
        n,ev,M,_,_,meta=point_model(((1,2,3),(2,3,5),(3,5,7)),p)
        models.append((n,ev,M,p,((-1,1),(-2,1),(-3,1)),meta))
    for p in (0,2,3):
        n,ev,M,_,_,meta=shifted_rectangle(p)
        models.append((n,ev,M,p,((-1,1),(-1,1),(0,1)),meta))
    for model in (gaussian_model(),extended_sharp(1,4)):
        n,ev,M,p,_,meta=model
        models.append((n,ev,M,p,((0,1),)*3,meta))
    models.append(point_model(((2,0,0),(3,0,0),(0,4,0),(0,0,5)),0,True))
    rows=[]
    for model in models:
        row=actual_control(model,counts);rows.append(row)
        print('complete mixed-rank shadow:',row['kind'],row['length'],row['field'],
              't',row['mixed_ideal_rank'],'r',row['contracted_rank'],
              'fixed',len(row['fixed_compulsory_labels']),'g',row['initial_ideals'],flush=True)
    synthetic=synthetic_controls(counts)
    inventory=direction_controls(counts)
    assert counts['guard_nonideal_forced_carrier_spans']>0
    assert counts['guard_other_objective_nondownset_bases']>0
    assert counts['guard_axis_labels_above_mixed_height']>0
    sources=['docs/mixed_rank_shadow_fans.md','scripts/run_mixed_rank_shadow_checks.py',
             'docs/mixed_ideal_rank_fans.md','docs/rational_shadow_bound.md',
             'docs/curve_carrier_shadows.md','scripts/run_mixed_ideal_rank_checks.py',
             'scripts/run_curve_carrier_shadow_checks.py','scripts/run_rational_shadow_bound_checks.py',
             'scripts/run_coordinate_slice_budget_checks.py','scripts/run_oracle_enumeration_checks.py',
             'scripts/run_planar_dyadic_bound_checks.py','scripts/run_primitive_wall_checks.py',
             'scripts/run_projection_tag_checks.py','scripts/run_quadric_cross_checks.py',
             'scripts/run_signed_face_checks.py','scripts/run_structural_complexity_checks.py',
             'scripts/run_univariate_cross_checks.py']
    out={'schema_version':1,'all_checks_passed':True,
         'scope':'Fixed compulsory carrier contraction of rank O(1+t), mixed hyperbolic/axis '
                 'shadow budgets, complete original fan equality, and direction inventory. '
                 'The exponent 1/3 is the published Dey parametric-matroid bound. '
                 'No unrestricted G3 improvement, sharp t exponent, polynomial separation, '
                 'novelty or optimized implementation is established.',
         'counters':dict(counts),'complete_actual_fans':rows,
         'synthetic_arbitrary_matroid_shadows':synthetic,'direction_inventory_controls':inventory,
         'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    target=ROOT/'results/research/mixed_rank_shadow_checks.json'
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_checks_passed':True,'counters':dict(counts),
                      'report':str(target.relative_to(ROOT))}),flush=True)


if __name__=='__main__':
    main()
