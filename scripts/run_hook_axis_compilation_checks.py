#!/usr/bin/env python3
"""Compile whole axis-hook fans from exact one-factor birth states."""
from pathlib import Path
from fractions import Fraction as Q
from bisect import bisect_right
from collections import Counter,defaultdict
from itertools import product
from random import Random
import ast,hashlib,json

from run_hook_axis_birth_checks import factor,encoded
from run_hook_bridge_checks import model,actual,normals,polygon,COUNTS as ORIGINAL_REPLAY
from run_hook_rank_departure_checks import complete_fan,plane,rank
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_signed_face_checks import field_basis,field_reduce
from run_structural_complexity_checks import dot,area2

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def fixed_profile(t,a):
    return (t*t,a*t*(a+t-2)//2,a*t*(a+t-2)//2)


def axis_factor(W,A,t):
    pure={(1,q):{q:1} for q in range(A)}
    pure.update({(0,q):{q+j:c for j,c in W.items() if q+j<A and c} for q in range(t,A)})
    head=[{i+j:c for j,c in W.items() if i+j<A and c} for i in range(t)]
    return pure,head


def compile_factor(pure,head,A,t,p):
    states=[]
    for k in range(-A-1,A+1):
        theta=Q(2*k+1,2);B,C,birth=factor(pure,head,theta,p)
        assert len(B)==A
        signature=(tuple(sorted(B)),tuple(sorted(birth.items())))
        if states and states[-1]['signature']==signature:
            states[-1]['hi']=None if k==A else k+1
            continue
        sx=sum(x for x,q in B);sq=sum(q for x,q in B)
        pure_profile=(t*sx,t*(A*t+sq),t*t*(A-sx)+A*t*(t-1)//2)
        states.append({'lo':None if k==-A-1 else k,'hi':None if k==A else k+1,
                       'signature':signature,'birth':dict(birth),'pure_profile':pure_profile})
        COUNTS['retained_factor_states']+=1
    assert states[0]['lo'] is None and states[-1]['hi'] is None
    assert all(S['hi']==T['lo'] for S,T in zip(states,states[1:]))
    COUNTS['all_integer_interval_factor_queries']+=2*A+2
    return {'states':states,'breaks':[S['hi'] for S in states[:-1]]}


def compile_tables(t,a,sources,p,translation=None):
    A=a-t;assert A>=t>=2
    return {'t':t,'a':a,'A':A,'n':4*a*t-2*t*t,
            'translation':fixed_profile(t,a) if translation is None else tuple(translation),
            'factors':[compile_factor(*source,A,t,p) for source in sources]}


def lookup(compiled,theta):
    k=bisect_right(compiled['breaks'],theta);S=compiled['states'][k]
    if not ((S['lo'] is None or S['lo']<theta) and(S['hi'] is None or theta<S['hi'])):
        raise ValueError('Query lies on a compiled factor boundary.')
    return k,S


def compiled_decision(compiled,w):
    t=compiled['t'];x,y,z=map(Q,w);assert min(x,y,z)>0
    indices=[];states=[]
    for table,theta in zip(compiled['factors'],[(x-t*z)/y,(x-t*y)/z]):
        k,S=lookup(table,theta);indices.append(k);states.append(S)
    py,pz=states[0]['pure_profile'],states[1]['pure_profile']
    profile=[compiled['translation'][k]+py[k]+pz[(0,2,1)[k]] for k in range(3)]
    N=set();choices=[]
    for arm,S in enumerate(states):
        if S['lo'] is not None:N.add((1,-S['lo'],-t) if not arm else(1,-t,-S['lo']))
        if S['hi'] is not None:N.add((-1,S['hi'],t) if not arm else(-1,t,S['hi']))
    for i,j in product(range(t),repeat=2):
        options={(0,0,0)}
        if i in states[0]['birth']:
            xx,q=states[0]['birth'][i];options.add((xx,q-i,-t*xx))
        if j in states[1]['birth']:
            xx,q=states[1]['birth'][j];options.add((xx,-t*xx,q-j))
        best=max(sorted(options),key=lambda e:dot(w,e));choices.append(best)
        for k in range(3):profile[k]-=best[k]
        for option in options:
            q=tuple(a-b for a,b in zip(best,option))
            if any(q):N.add(q)
    if any(dot(w,q)==0 for q in N):raise ValueError('Query lies on a nontrivial maximum comparison.')
    assert all(dot(w,q)>0 for q in N)
    COUNTS['compiled_profiles_without_field_elimination']+=1
    return (tuple(indices),tuple(choices)),tuple(profile),N


def compiled_patch(compiled,w):
    key,profile,N=compiled_decision(compiled,w)
    return key,profile,N,polygon(N)


def hull(points):
    points=sorted(set(points));assert len(points)>=3
    def turn(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    lower=[]
    for p in points:
        while len(lower)>=2 and turn(lower[-2],lower[-1],p)<=0:lower.pop()
        lower.append(p)
    upper=[]
    for p in reversed(points):
        while len(upper)>=2 and turn(upper[-2],upper[-1],p)<=0:upper.pop()
        upper.append(p)
    out=lower[:-1]+upper[:-1];assert len(out)>=3 and area2(out)>0
    return tuple(out)


def walk_compiled(compiled,verbose=False):
    n,t,A=compiled['n'],compiled['t'],compiled['A'];L=2*n+1
    queue=[((1,L,L*L),None)];cells={};profile_areas=defaultdict(Q);vertices=defaultdict(list)
    field_queries_before=COUNTS['all_integer_interval_factor_queries']
    while queue:
        w,previous=queue.pop();key,profile,N,P=compiled_patch(compiled,w)
        assert previous is None or previous!=key
        if key in cells:continue
        cells[key]=(profile,P);profile_areas[profile]+=area2(P);vertices[profile].extend(P)
        for x,y in zip(P,P[1:]+P[:1]):
            if any(x[i]==y[i]==0 for i in range(3)):continue
            q=next(q for q in N if dot(q,x)==dot(q,y)==0)
            nw,_=exact_adjacent_weight(n,x,y,q);queue.append((nw,key))
        if verbose and len(cells)%250==0:print('Live compiled axis fan',t,A,len(cells),'patches',len(vertices),'profiles',flush=True)
    assert COUNTS['all_integer_interval_factor_queries']==field_queries_before
    assert sum(profile_areas.values())==1
    fan={profile:hull(points) for profile,points in vertices.items()}
    assert all(area2(fan[profile])==area for profile,area in profile_areas.items())
    assert sum(map(area2,fan.values()))==1
    COUNTS['independent_complete_compiled_covers']+=1
    COUNTS['merged_profile_normal_polygons']+=len(fan)
    COUNTS['redundant_refinement_patches_removed']+=len(cells)-len(fan)
    return fan,cells


def axis_model(t,a,Ws,p):
    data=model(t,a);co=[]
    for g,y,z in data[2]:
        arm=0 if y>=t else 1;q,j=(y-t,z) if not arm else(z-t,y)
        co.append(Ws[arm].get(q,0) if g==1 and not j else 0)
    fixed=tuple(sum(e[k] for e in data[6]) for k in range(3));assert fixed==fixed_profile(t,a)
    return data,actual(data,co,p,True)


def centroid(P):return tuple(sum(v[k] for v in P)/len(P) for k in range(3))


def reference_controls(t,a,Ws,p,compiled,fan,cells,full=True):
    data,cols=axis_model(t,a,Ws,p);translation=compiled['translation']
    if full:
        original=complete_fan(cols,compiled['n'],p)
        reference={tuple(translation[k]+sum(e[k] for e in B) for k in range(3)):(B,hull(P)) for B,P in original.items()}
        assert len(reference)==len(original) and set(reference)==set(fan)
        for profile,(B,P) in reference.items():
            assert set(P)==set(fan[profile])
            COUNTS['independent_entire_merged_original_normal_cones']+=1
        for profile,P in cells.values():
            B,_=reference[profile];N=normals(B,cols,p)
            assert all(dot(q,v)>=0 for q in N for v in P)
        COUNTS['independent_full_original_axis_fan_equalities']+=1
        return {'complete_original_cones':len(original),'whole_normal_cones_equal':True}
    selected=list(cells.values())[::max(1,len(cells)//5)][:6]
    for profile,P in selected:
        # A refinement-cell centroid is strict for every nonzero comparison used by it.
        w=centroid(P);B=field_basis(sorted(cols,key=lambda e:(dot(w,e),e)),cols.__getitem__,p)[0]
        expected=tuple(translation[k]+sum(e[k] for e in B) for k in range(3));assert profile==expected
        COUNTS['independent_large_original_quotient_profile_spot_checks']+=1
    return {'complete_original_cones':None,'independent_original_profile_spot_checks':len(selected),
            'scope':'Original quotient and commuting multiplication replayed. Full original fan not independently enumerated at this size.'}


def abstract_controls(rng):
    A,t,p=6,2,101;a=A+t;sources=[];cols={}
    for arm in range(2):
        pure={(1,q):{q:1} for q in range(A)}
        pure.update({(0,q):{i:c for i in range(A) if (c:=rng.randrange(p))} for q in range(t,A)})
        head=[{i:c for i in range(A) if (c:=rng.randrange(p))} for _ in range(t)]
        if not arm:head[0]=dict(pure[0,t])
        sources.append((pure,head))
        for (x,q),v in pure.items():
            for j in range(t):
                e=(x,t+q,t*(1-x)+j) if not arm else(x,t*(1-x)+j,t+q)
                cols[e]={i*t+j:c for i,c in v.items()} if not arm else{A*t+j*A+i:c for i,c in v.items()}
    for i,j in product(range(t),repeat=2):
        cols[0,t+i,t+j]={**{r*t+j:c for r,c in sources[0][1][i].items()},**{A*t+i*A+r:c for r,c in sources[1][1][j].items()}}
    compiled=compile_tables(t,a,sources,p,(0,0,0));fan,cells=walk_compiled(compiled)
    original=complete_fan(cols,compiled['n'],p)
    reference={tuple(sum(e[k] for e in B) for k in range(3)):hull(P) for B,P in original.items()}
    assert set(reference)==set(fan) and all(set(reference[v])==set(fan[v]) for v in fan)
    COUNTS['independent_complete_non_Toeplitz_replicated_fans']+=1
    return {'long_dimension':A,'short_dimension':t,'characteristic':p,'refinement_patches':len(cells),'actual_profile_cones':len(fan),
            'scope':'Arbitrary replicated matrix, including a head in the pure lower span. No algebra realization claimed.'}


def interior_remainder_controls(compiled,fan,Ws,p):
    t,a=compiled['t'],compiled['a'];assert (t,a,p)==(4,17,1009)
    data,cols=axis_model(t,a,Ws,p);known={plane(q) for q in data[7]}
    w0=(Q(9,44),Q(17,44),Q(18,44));qs=[(1,-9,8),(1,9,-9)]
    theta=(w0[0]-t*w0[2])/w0[1];phi=(w0[0]-t*w0[1])/w0[2]
    assert theta==Q(-63,17) and phi==Q(-59,18)
    assert theta.denominator>1 and phi.denominator>1 and dot((1,-t,-t),w0)
    incident=sum(w0 in P for P in fan.values());assert incident==4
    facets=[]
    for q in qs:
        assert q not in known and dot(q,w0)==0
        edges=set()
        for P in fan.values():
            for x,y in zip(P,P[1:]+P[:1]):
                if (x==w0 or y==w0) and dot(q,x)==dot(q,y)==0:edges.add(tuple(sorted((x,y))))
        assert len(edges)==2
        for edge in sorted(edges):
            mid=tuple((x+y)/2 for x,y in zip(*edge));bases=[]
            for sign in (-1,1):
                w,_=exact_adjacent_weight(compiled['n'],*edge,tuple(sign*x for x in q))
                B=frozenset(field_basis(sorted(cols,key=lambda e:(dot(w,e),e)),cols.__getitem__,p)[0]);N=normals(B,cols,p)
                assert all(dot(w,r)>0 for r in N) and all(dot(r,v)>=0 for r in N for v in edge)
                profile=tuple(compiled['translation'][i]+sum(e[i] for e in B) for i in range(3))
                assert set(hull(polygon(N)))==set(fan[profile]);bases.append(B)
                COUNTS['independent_interior_remainder_whole_neighbor_cones']+=1
            assert bases[0]!=bases[1]
            e,f=((0,10,5),(1,1,13)) if q==(1,-9,8) else((0,4,9),(1,13,0))
            left=next(B for B in bases if e in B and f not in B)
            piv=field_basis(tuple(left),cols.__getitem__,p)[1]
            rem,representation=field_reduce(cols[f],piv,p)
            assert not rem and representation.get(e,0)
            exchanged=(left-{e})|{f}
            assert len(field_basis(tuple(exchanged),cols.__getitem__,p)[0])==len(left)
            assert sum(dot(mid,v) for v in exchanged)==sum(dot(mid,v) for v in left)
            COUNTS['independent_interior_remainder_single_minimum_face_exchanges']+=1
            assert e[0]!=f[0] and plane(tuple(x-y for x,y in zip(e,f)))==q
            assert (e in data[4][0] and f in data[4][1])or(e in data[4][1] and f in data[4][0])
            tau=dot(mid,e);assert dot(mid,f)==tau
            cheap=tuple(v for v in cols if dot(mid,v)<tau);r=rank(cheap,cols,p)
            assert r==(97 if q==(1,-9,8) else 80)
            increments=[rank((*cheap,*extra),cols,p)-r for extra in [(e,),(f,),(e,f)]]
            assert increments==[1,1,1]
            COUNTS['independent_interior_remainder_entire_cheaper_prefix_ranks']+=1
            facets.append({'primitive_plane':q,'endpoints':[[str(v) for v in w] for w in edge],
                           'exchange':sorted((e,f)),'strictly_cheaper_rank':r,'single_and_joint_increments':increments,
                           'whole_neighbor_basis_symmetric_difference':len(bases[0]^bases[1]),
                           'single_exchange_on_minimum_face_verified':True})
    return {'t':t,'a':a,'characteristic':p,'common_weight':list(map(str,w0)),
            'theta':str(theta),'phi':str(phi),'incident_original_cones':incident,'facets':facets,
            'scope':'Two original remainder directions cross inside a factor-order cell, away from the singular price line. This finite four-cone certificate does not prove an asymptotic lower rate.'}


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    prior=json.loads((ROOT/'results/research/hook_axis_birth_checks.json').read_text());assert prior['status']=='passed'
    for f,h in prior['source_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
    rng=Random(2026102353);cases=[];interior_crossing=None
    configs=[(2,8,[{0:1,1:3,4:2},{0:1,2:4,5:6}],101,True),
             (3,10,[{1:1,4:1},{2:1,5:1}],2,True),
             (2,9,[{}, {0:1,3:1}],101,True),
             (3,12,[{0:1},{0:1}],0,True),
             (4,14,[{0:1,6:1,7:1},{0:1,8:2,9:3}],1009,True),
             (4,17,[{0:1,5:2,6:3,11:5},{0:1,7:7,8:11,12:13}],1009,True),
             (6,26,[{0:1,8:1,9:1},{0:1,10:2,11:3}],1009,False),
             (8,40,[{0:1,10:1,11:1},{0:1,12:2,13:3}],1009,False)]
    for t,a,Ws,p,full in configs:
        A=a-t;Ws=[{q:(c%p if p else Q(c)) for q,c in W.items() if q<A and(c%p if p else c)} for W in Ws]
        compiled=compile_tables(t,a,[axis_factor(W,A,t) for W in Ws],p)
        fan,cells=walk_compiled(compiled,True);reference=reference_controls(t,a,Ws,p,compiled,fan,cells,full)
        if (t,a,p)==(4,17,1009):interior_crossing=interior_remainder_controls(compiled,fan,Ws,p)
        for _ in range(12):
            w=(rng.randrange(1,10000),rng.randrange(1,1000),rng.randrange(1,1000))
            theta=Q(w[0]-t*w[2],w[1]);phi=Q(w[0]-t*w[1],w[2])
            if theta.denominator==1 or phi.denominator==1:continue
            try:key,profile,N=compiled_decision(compiled,w)
            except ValueError as error:
                assert 'boundary' in str(error) or 'comparison' in str(error)
                continue
            expected,original_N,F=encoded(t,a,Ws,p,w)
            assert profile==tuple(expected[k]+compiled['translation'][k] for k in range(3))
            COUNTS['independent_fresh_field_vs_compiled_profile_queries']+=1
        serialization=[{'profile':v,'polygon':[[str(x) for x in w] for w in fan[v]]} for v in sorted(fan)]
        case={'t':t,'a':a,'A':A,'length':compiled['n'],'characteristic':p,'operators':[{str(q):str(c) for q,c in W.items()} for W in Ws],
              'compiled_factor_states':[len(F['states']) for F in compiled['factors']],
              'refinement_patches':len(cells),'actual_initial_ideal_cones':len(fan),'exact_total_area':'1',
              'normal_cones_sha256':hashlib.sha256(json.dumps(serialization,sort_keys=True).encode()).hexdigest(),
              'reference':reference}
        if len(fan)<=100:case['complete_profile_normal_cones']=serialization
        cases.append(case);print('Compiled axis fan passed',t,A,'patches',len(cells),'actual',len(fan),flush=True)
    abstract=abstract_controls(rng)
    assert interior_crossing is not None
    files={'docs/hook_axis_fan_compilation.md','scripts/run_hook_axis_compilation_checks.py','docs/hook_axis_birth_curves.md',
           'docs/oracle_fan_enumeration.md','results/research/hook_axis_birth_checks.json'}
    pending=['scripts/run_hook_axis_compilation_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=([node.module] if isinstance(node,ast.ImportFrom) else[v.name for v in node.names] if isinstance(node,ast.Import) else[])
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and(ROOT/target).exists() and target not in files:files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_axis_compilation_checks.py',
            'scope':'Exact coefficient-dependent compilation and full actual profile-fan traversal for zero-outer axis hooks. Refinement cells are merged by actual profiles; complete original normal cones are checked independently at small sizes.',
            'uniform_result':{'field_precomputation':'O(A) one-factor price intervals, O(A^4) field operations and O(A^2) field storage',
                              'profile_queries_after_compilation':'O(t^2+log A) rational/integer cost operations, no field elimination',
                              'complete_fan':'O(A^4) total field operations; O(K) exact patch-oracle queries with K=O(A^2+A*t^4); no optimal rational-operation or bit-complexity claim',
                              'merging':'The convex hull of all patches with one original profile has exactly their summed area and recovers one full original normal cone.',
                              'scope_guard':'No coefficient-independent common list, new unrestricted exponent, actual asymptotic lower excess or strongest-current novelty claim.'},
            'counts':dict(sorted(COUNTS.items())),'original_quotient_replay_counts':dict(ORIGINAL_REPLAY),
            'axis_cases':cases,'arbitrary_replicated_matrix_control':abstract,
            'interior_original_remainder_crossing':interior_crossing,
            'limits':['Nonzero outer images and arbitrary transverse terms are outside the two-factor formula.',
                      'Larger cases use a complete exact compiled cover and original profile spot checks, not an independent full original fan enumeration.',
                      'Finite fan counts are not a growing-thickness rate or evidence of polynomial excess.',
                      'The unrestricted sharp G_3 target, exact strongest-current priority and final original-goal paper remain incomplete.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_axis_compilation_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}),flush=True)


if __name__=='__main__':main()
