#!/usr/bin/env python3
"""Exact shift-order and linear-merge checks for actual axis-hook fans."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter,defaultdict
from random import Random
import hashlib,json

from run_hook_axis_compilation_checks import (axis_factor,compile_factor,compile_tables,lookup,
    compiled_decision,hull,centroid,reference_controls)
from run_hook_axis_birth_checks import direct_births,factor,assignment
from run_hook_bridge_checks import polygon
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_structural_complexity_checks import dot,area2
from run_signed_face_checks import field_basis,field_reduce

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def incremental_factor(pure,head,A,t,p):
    """Nested pure bases, with one rank-one coordinate update per exchange."""
    B,piv=field_basis(tuple((1,q) for q in range(A)),pure.__getitem__,p);assert len(B)==A
    B=set(B);D={};H=[]
    for e,v in pure.items():
        rem,rep=field_reduce(v,piv,p);assert not rem;D[e]=rep
    for v in head:
        rem,rep=field_reduce(v,piv,p);assert not rem;H.append(rep)
    states=[];changes=[];lowers=sorted(e for e in pure if e[0]==0)
    for k in range(-A-1,A+1):
        while k!=-A-1:
            uppers=sorted(b for b in B if b[0]==1)
            pair=next(((e,b) for e in lowers if e not in B for b in uppers
                       if e[1]-b[1]==k and D[e].get(b,0)),None)
            if pair is None:break
            e,b=pair;v=dict(D[e]);inverse=pow(v[b],-1,p) if p else Q(1)/v[b]
            def update(rep):
                c=rep.get(b,0)
                if not c:return dict(rep)
                lam=(c*inverse)%p if p else c*inverse
                out={r:c for r,c in rep.items() if r!=b}
                for r,c in v.items():
                    if r==b:continue
                    value=out.get(r,0)-lam*c
                    if p:value%=p
                    if value:out[r]=value
                    else:out.pop(r,None)
                out[e]=lam
                return out
            D={label:update(rep) for label,rep in D.items()};H=[update(rep) for rep in H]
            B.remove(b);B.add(e);changes.append((k,e,b))
            assert all(D[label]=={label:1} for label in B)
        theta=Q(2*k+1,2)
        C={b:{i:c for i,rep in enumerate(H) if(c:=rep.get(b,0))} for b in B}
        birth=assignment(C,{e:theta*e[0]+e[1] for e in B},p)
        signature=(tuple(sorted(B)),tuple(sorted(birth.items())))
        if states and states[-1]['signature']==signature:
            states[-1]['hi']=None if k==A else k+1
            continue
        sx=sum(x for x,q in B);sq=sum(q for x,q in B)
        states.append({'lo':None if k==-A-1 else k,'hi':None if k==A else k+1,
                       'signature':signature,'birth':birth,
                       'pure_profile':(t*sx,t*(A*t+sq),t*t*(A-sx)+A*t*(t-1)//2)})
    assert len(changes)<=A-t
    assert len({b for k,e,b in changes})==len(changes)==len({e for k,e,b in changes})
    COUNTS['incremental_factor_compilations']+=1
    COUNTS['rank_one_pure_basis_exchanges']+=len(changes)
    return {'states':states,'breaks':[S['hi'] for S in states[:-1]]},changes


def linear_decision(compiled,w):
    """Requires the actual axis-shift factors, rather than arbitrary AB3 data."""
    t=compiled['t'];x,y,z=map(Q,w);assert min(x,y,z)>0
    indices=[];states=[];N=set()
    for arm,(table,theta) in enumerate(zip(compiled['factors'],[(x-t*z)/y,(x-t*y)/z])):
        k,S=lookup(table,theta);indices.append(k);states.append(S)
        if S['lo'] is not None:N.add((1,-S['lo'],-t) if not arm else(1,-t,-S['lo']))
        if S['hi'] is not None:N.add((-1,S['hi'],t) if not arm else(-1,t,S['hi']))
    py,pz=states[0]['pure_profile'],states[1]['pure_profile']
    profile=[compiled['translation'][k]+py[k]+pz[(0,2,1)[k]] for k in range(3)]
    lists=[];signatures=[];steps=0
    for arm,S in enumerate(states):
        positive=[];signs=[];previous=None;absent=False
        for i in range(t):
            steps+=1
            if i not in S['birth']:
                absent=True;signs.append(None);continue
            assert not absent,'An absent birth must be followed only by absent births.'
            delta,q=S['birth'][i]
            v=(delta,q-i,-t*delta) if not arm else(delta,-t*delta,q-i)
            value=dot(w,v)
            if previous is not None:assert value<=previous,'The axis shift order is required.'
            previous=value
            if value==0:raise ValueError('Query lies on a gain-sign boundary.')
            sign=1 if value>0 else -1;signs.append(sign);N.add(tuple(sign*c for c in v))
            if value>0:positive.append((value,i,v))
        signatures.append(tuple(signs));lists.append(positive)
    i=j=0;word=[];last=None
    while i<len(lists[0]) or j<len(lists[1]):
        steps+=1
        if i<len(lists[0]) and j<len(lists[1]):
            av,bv=lists[0][i][0],lists[1][j][0]
            if av==bv:raise ValueError('Query lies on a positive cross-arm comparison.')
            arm=0 if av>bv else 1
        else:arm=0 if i<len(lists[0]) else 1
        _,head,v=lists[arm][i if not arm else j]
        multiplicity=t-(j if not arm else i)
        for k in range(3):profile[k]-=multiplicity*v[k]
        if last is not None and last[0]!=arm:
            q=tuple(a-b for a,b in zip(last[1],v));assert any(q);N.add(q)
        word.append((arm,head));last=(arm,v)
        if not arm:i+=1
        else:j+=1
    assert steps<=4*t and len(N)<=4*t+4
    assert all(dot(w,q)>0 for q in N)
    COUNTS['linear_profile_queries']+=1
    COUNTS['maximum_linear_loop_steps']=max(COUNTS['maximum_linear_loop_steps'],steps)
    COUNTS['maximum_patch_normal_count']=max(COUNTS['maximum_patch_normal_count'],len(N))
    return (tuple(indices),tuple(signatures),tuple(word)),tuple(profile),N


def walk_linear(compiled):
    n,t=compiled['n'],compiled['t'];L=2*n+1
    queue=[((1,L,L*L),None)];cells={};areas=defaultdict(Q);vertices=defaultdict(list)
    while queue:
        w,previous=queue.pop();key,profile,N=linear_decision(compiled,w)
        assert previous is None or key!=previous
        if key in cells:continue
        P=polygon(N);assert area2(P)>0
        # Independently inspect every quadratic winner inequality on the entire new patch.
        _,old_profile,old_N=compiled_decision(compiled,centroid(P))
        assert profile==old_profile and all(dot(q,v)>=0 for q in old_N for v in P)
        COUNTS['whole_linear_patches_inside_quadratic_patches']+=1
        cells[key]=(profile,P);areas[profile]+=area2(P);vertices[profile].extend(P)
        for x,y in zip(P,P[1:]+P[:1]):
            if any(x[i]==y[i]==0 for i in range(3)):continue
            q=next(q for q in N if dot(q,x)==dot(q,y)==0)
            nw,_=exact_adjacent_weight(n,x,y,q);queue.append((nw,key))
    assert sum(areas.values())==1
    fan={v:hull(points) for v,points in vertices.items()}
    assert all(area2(fan[v])==a for v,a in areas.items())
    assert sum(map(area2,fan.values()))==1
    COUNTS['complete_linear_refinement_covers']+=1
    return fan,cells


def interval_order_controls(rng):
    cases=[]
    for p in [0,2,3,101]:
        for A,t in [(4,2),(7,3),(10,4),(15,5)]:
            for mode in ['zero','unit','nonunit','dense','sparse']:
                W={} if mode=='zero' else {0:1} if mode=='unit' else {
                    i:rng.randrange(p or 7) for i in range(t//2 if mode=='nonunit' else 0,A)
                    if mode!='sparse' or rng.randrange(3)==0}
                W={q:(c%p if p else Q(c)) for q,c in W.items() if c}
                pure,head=axis_factor(W,A,t);table=compile_factor(pure,head,A,t,p);states=0
                incremental,changes=incremental_factor(pure,head,A,t,p)
                assert incremental==table
                COUNTS['entire_incremental_factor_tables_equal_fresh_elimination']+=1
                for S in table['states']:
                    seen_absent=False
                    for i in range(t):
                        if i not in S['birth']:seen_absent=True;continue
                        assert not seen_absent
                        if not i:continue
                        delta,q=S['birth'][i];old_delta,old_q=S['birth'][i-1]
                        slope=delta-old_delta;intercept=q-old_q-1
                        if S['lo'] is None and S['hi'] is None:assert slope==0 and intercept<=0
                        if S['lo'] is None:assert slope>=0
                        else:assert slope*S['lo']+intercept<=0
                        if S['hi'] is None:assert slope<=0
                        else:assert slope*S['hi']+intercept<=0
                        COUNTS['whole_interval_shift_gain_inequalities']+=1
                    states+=1
                # These are complete original pure-prefix rank tests, not only assigned labels.
                theta=Q(-3,2);B,C,assigned=factor(pure,head,theta,p)
                direct=direct_births(pure,head,theta,p)
                assert direct==[theta*assigned[i][0]+assigned[i][1] if i in assigned else None for i in range(t)]
                assert all(direct[i] is None or direct[i-1] is not None and direct[i]-i<=direct[i-1]-(i-1) for i in range(1,t))
                COUNTS['independent_whole_prefix_birth_values']+=t
                cases.append({'A':A,'t':t,'characteristic':p,'mode':mode,'complete_intervals':states,
                              'rank_one_exchanges':len(changes),'upper_and_lower_labels_never_reenter':True})
    COUNTS['complete_shift_factor_fixtures']+=len(cases)
    return cases


def scope_counterexample():
    A,t,p=4,2,101;pure={(1,q):{q:1} for q in range(A)}
    pure.update({(0,q):{} for q in range(t,A)})
    head=[{0:1},{3:1}];births=direct_births(pure,head,Q(1,2),p)
    assert births==[Q(1,2),Q(7,2)] and births[1]-1>births[0]
    return {'A':A,'t':t,'characteristic':p,'theta':'1/2','births':list(map(str,births)),
            'head_upper_indices':[0,3],'pure_lowers':'all zero',
            'scope':'An arbitrary replicated factor fails the axis shift order. No algebra realization claimed.'}


def arbitrary_incremental_controls(rng):
    cases=[]
    for p in [0,2,3,101]:
        for A,t in [(5,2),(8,3),(11,4)]:
            for mode in ['dense','dependent']:
                pure={(1,q):{**{i:c for i in range(q) if(c:=rng.randrange(p or 5))},q:1} for q in range(A)}
                v={i:c for i in range(A) if(c:=rng.randrange(p or 5))}
                pure.update({(0,q):dict(v) if mode=='dependent' else {i:c for i in range(A) if(c:=rng.randrange(p or 5))} for q in range(t,A)})
                head=[{i:c for i in range(A) if(c:=rng.randrange(p or 5))} for j in range(t)]
                head[0]=dict(pure[0,t])
                fresh=compile_factor(pure,head,A,t,p);incremental,changes=incremental_factor(pure,head,A,t,p)
                assert fresh==incremental
                COUNTS['arbitrary_replicated_incremental_tables_equal_fresh_elimination']+=1
                cases.append({'A':A,'t':t,'characteristic':p,'mode':mode,'upper_basis':'nonidentity triangular basis',
                              'head_in_lower_span':True,'rank_one_exchanges':len(changes),
                              'scope':'Arbitrary replicated factor. No algebra realization claimed.'})
    return cases


def serialization(fan):
    return [{'profile':v,'polygon':[[str(x) for x in w] for w in fan[v]]} for v in sorted(fan)]


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    prior_path='results/research/hook_axis_compilation_checks.json'
    prior=json.loads((ROOT/prior_path).read_text());assert prior['status']=='passed'
    for f,h in prior['source_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
    rng=Random(20261023123);controls=interval_order_controls(rng)
    arbitrary=arbitrary_incremental_controls(rng);cases=[]
    for ci,old in enumerate(prior['axis_cases']):
        t,a,p=old['t'],old['a'],old['characteristic'];A=a-t
        Ws=[{int(q):int(c) for q,c in W.items()} for W in old['operators']]
        compiled=compile_tables(t,a,[axis_factor(W,A,t) for W in Ws],p)
        for arm,W in enumerate(Ws):
            table,changes=incremental_factor(*axis_factor(W,A,t),A,t,p)
            assert table==compiled['factors'][arm]
            compiled['factors'][arm]=table
            COUNTS['entire_incremental_factor_tables_equal_fresh_elimination']+=1
        fan,cells=walk_linear(compiled)
        digest=hashlib.sha256(json.dumps(serialization(fan),sort_keys=True).encode()).hexdigest()
        assert digest==old['normal_cones_sha256'] and len(fan)==old['actual_initial_ideal_cones']
        COUNTS['exact_entire_profile_fan_equalities_with_prior_report']+=1
        reference=reference_controls(t,a,Ws,p,compiled,fan,cells,True) if ci<2 else None
        for _ in range(30):
            w=(rng.randrange(1,10000),rng.randrange(1,1000),rng.randrange(1,1000))
            try:_,profile,N=linear_decision(compiled,w)
            except ValueError:continue
            _,old_profile,old_N=compiled_decision(compiled,w);assert profile==old_profile
            COUNTS['fresh_linear_vs_quadratic_profile_queries']+=1
        cases.append({'t':t,'a':a,'A':A,'length':compiled['n'],'characteristic':p,
                      'linear_refinement_patches':len(cells),'prior_refinement_patches':old['refinement_patches'],
                      'actual_initial_ideal_cones':len(fan),'exact_total_area':'1',
                      'normal_cones_sha256':digest,'independent_original_fan_replay':reference})
        print('Linear axis fan passed',t,A,'patches',len(cells),'actual',len(fan),flush=True)
    files=set(prior['source_sha256'])|{prior_path,'docs/hook_axis_linear_oracle.md','scripts/run_hook_axis_linear_checks.py'}
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_axis_linear_checks.py',
            'uniform_result':{'whole_prefix_shift_inequality':'alpha_(i+1)(theta) <= alpha_i(theta)+1; absent births form a suffix',
                              'compiled_profile_query':'O(t+log A) rational/integer operations with no field elimination',
                              'patch_normals':'at most 4*t+4 integer inequalities, with additional rational polygon-construction cost',
                              'compilation':'O(A^3+A^2*t^2) field operations and O(A^2) field storage, using at most A-t rank-one pure basis exchanges per factor.',
                              'scope':'Actual zero-outer axis hooks over all fields. The arbitrary replicated matrix specialization is not asserted.'},
            'counts':dict(COUNTS),'whole_interval_shift_controls':controls,'axis_fans':cases,
            'arbitrary_replicated_incremental_controls':arbitrary,
            'arbitrary_replicated_scope_counterexample':scope_counterexample(),
            'limits':['No new asymptotic fan-count exponent or infinite polynomial lower excess.',
                      'Finite exact full fans are compared with the earlier report; two original quotient fans are freshly replayed.',
                      'Sign inequalities can introduce additional refinement patches; the actual full profile fans agree exactly.',
                      'No optimal bit-complexity, basis-output or smooth-number unranking speed claim.',
                      'The unrestricted extremal goal, strongest-current novelty and final paper remain incomplete.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_axis_linear_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}),flush=True)


if __name__=='__main__':main()
