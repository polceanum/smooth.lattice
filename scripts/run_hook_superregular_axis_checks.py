#!/usr/bin/env python3
"""Exact generic banded-Toeplitz birth and whole axis-fan checks."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from itertools import combinations,product
from random import Random
import ast,hashlib,json

from run_hook_axis_compilation_checks import axis_factor,compile_factor,fixed_profile,axis_model,hull,centroid
from run_hook_axis_linear_checks import walk_linear,linear_decision
from run_hook_axis_birth_checks import factor,direct_births,curves
from run_hook_bridge_checks import normals
from run_hook_rank_departure_checks import complete_fan,plane
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import dot,area2

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def degree_factor(A,t,D):
    assert 0<=D<A
    regular=sorted(set(range(t,min(A,t+D)+1))|{A});states=[]
    for k in range(-A-1,A+1):
        theta=Q(2*k+1,2)
        def cost(r):
            s=A-r
            return r*theta+r*(r-1)//2+s*t+s*(s-1)//2
        r=min(regular,key=cost);s=A-r
        B=tuple(sorted([(1,q) for q in range(r)]+[(0,q) for q in range(t,t+s)]))
        merged=sorted([(1,q) for q in range(D)]+[(0,q) for q in range(t,A)],key=lambda e:theta*e[0]+e[1])
        birth={}
        for i in range(t):
            options=([(1,D+i)] if D+i<A else[])+([merged[A-i-1]] if A-i<=len(merged) else[])
            birth[i]=min(options,key=lambda e:theta*e[0]+e[1])
            assert birth[i] in B
        signature=(B,tuple(sorted(birth.items())))
        if states and states[-1]['signature']==signature:
            states[-1]['hi']=None if k==A else k+1;continue
        sx=r;sq=sum(q for x,q in B)
        states.append({'lo':None if k==-A-1 else k,'hi':None if k==A else k+1,'signature':signature,
                       'birth':birth,'pure_profile':(t*sx,t*(A*t+sq),t*t*(A-sx)+A*t*(t-1)//2)})
    return {'states':states,'breaks':[S['hi'] for S in states[:-1]]}


def degree_compile(t,A,Ds):
    return {'t':t,'A':A,'a':A+t,'n':4*(A+t)*t-2*t*t,'translation':fixed_profile(t,A+t),
            'factors':[degree_factor(A,t,D) for D in Ds]}


def interval_factor(A,t,spec):
    if spec is None:
        B=tuple((1,q) for q in range(A));birth={};base=[{'lo':None,'hi':None,'signature':(B,()),'birth':birth}]
    else:
        nu,D=spec;assert 0<=nu<=D<A
        effective=A-nu;small=min(t,effective)
        source=degree_factor(effective,small,D-nu);base=[]
        for S in source['states']:
            B=tuple(sorted([(1,q) for q in range(nu)]+[(x,q+nu if x else q) for x,q in S['signature'][0]]))
            birth={i:(x,q+nu if x else q) for i,(x,q) in S['birth'].items()}
            base.append({'lo':None if S['lo'] is None else S['lo']-nu,
                         'hi':None if S['hi'] is None else S['hi']-nu,
                         'signature':(B,tuple(sorted(birth.items()))),'birth':birth})
    for S in base:
        sx=sum(x for x,q in S['signature'][0]);sq=sum(q for x,q in S['signature'][0])
        S['pure_profile']=(t*sx,t*(A*t+sq),t*t*(A-sx)+A*t*(t-1)//2)
    return {'states':base,'breaks':[S['hi'] for S in base[:-1]]}


def interval_compile(t,A,specs):
    return {'t':t,'A':A,'a':A+t,'n':4*(A+t)*t-2*t*t,'translation':fixed_profile(t,A+t),
            'factors':[interval_factor(A,t,S) for S in specs]}


def matching_rank(A,t,D,u,s,i):
    h=max(0,min(i,A-u,D+i-u));b=max(t,u+h)
    return h+max(0,min(s,A-b,D+t+s-b))


def admitted(raw,A,t,Ds):
    x,y,z=raw;a=A+t
    def centers(S):
        if S is None:return [0,A]
        if isinstance(S,tuple):
            nu,D=S;return [0,nu,D,A-D,A-nu,A]
        return [0,S,A-S,A]
    CY,CZ=map(centers,Ds)
    if any(abs(abs(y)-v)<=4*t for v in CY)or any(abs(abs(z)-v)<=4*t for v in CZ):return True
    if x==0:return y*z<0 and abs(y+z)<=64*t
    return abs(abs(y-z)-a)<=64*t


def all_minors(W,A,D,p,nu=0):
    nontrivial=trivial=0
    for m in range(1,A+1):
        for R,C in product(combinations(range(A),m),repeat=2):
            cols={j:{r:W.get(r-c,0) for r in R if W.get(r-c,0)} for j,c in enumerate(C)}
            rank=len(field_basis(tuple(cols),cols.__getitem__,p)[0])
            allowed=all(nu<=r-c<=D for r,c in zip(R,C))
            assert (rank==m)==allowed,(A,D,p,R,C)
            if allowed:nontrivial+=1
            else:trivial+=1
    COUNTS['independent_nontrivial_Toeplitz_minors']+=nontrivial
    COUNTS['independent_structurally_zero_Toeplitz_minors']+=trivial
    return {'nontrivial_nonzero':nontrivial,'structurally_zero':trivial}


def field_controls(W,A,t,D,p,entire=False):
    pure,head=axis_factor(W,A,t);fresh=compile_factor(pure,head,A,t,p);predicted=degree_factor(A,t,D)
    assert fresh==predicted
    COUNTS['entire_degree_tables_equal_field_elimination']+=1
    for u,s,i in product(range(A+1),range(A-t+1),range(t+1)):
        if not entire and (u+s+i)%7:continue
        sources=[{q:1} for q in range(u)]+head[:i]+[pure[0,q] for q in range(t,t+s)]
        rank=len(field_basis(tuple(range(len(sources))),sources.__getitem__,p)[0])
        assert rank==u+matching_rank(A,t,D,u,s,i)
        COUNTS['independent_entire_pure_plus_head_prefix_ranks']+=1
    for theta in [Q(-2*A-1,2),Q(-3,2),Q(3,2),Q(2*A+1,2)]:
        expected=direct_births(pure,head,theta,p);B,C,birth=factor(pure,head,theta,p)
        assert expected==[theta*birth[i][0]+birth[i][1] for i in range(t)]
        COUNTS['independent_entire_prefix_births']+=t
    return {'A':A,'t':t,'D':D,'characteristic':p,'operators':{str(q):str(c) for q,c in W.items()}}


def birth_pieces(table,t):
    out=[[] for i in range(t)]
    for S in table['states']:
        for i,(x,q) in S['birth'].items():
            P={'lo':None if S['lo'] is None else Q(S['lo']),'hi':None if S['hi'] is None else Q(S['hi']),'x':x,'q':q}
            if out[i] and(out[i][-1]['x'],out[i][-1]['q'])==(x,q):out[i][-1]['hi']=P['hi']
            else:out[i].append(P)
    return out


def arc_controls():
    cases=[]
    for A,t,Ds in [(6,2,(0,5)),(8,3,(2,4)),(12,4,(6,6)),(14,3,(3,12)),(20,2,(10,10)),(144,2,(72,72))]:
        C=degree_compile(t,A,Ds);maps=curves(*[birth_pieces(F,t) for F in C['factors']],A,t);count=0
        for pieces in maps.values():
            for P in pieces:
                i,j=P['indices'];x,q=P['Y']['x'],P['Y']['q'];z,r=P['Z']['x'],P['Z']['q']
                raw=(x-z,q-i+t*z,-t*x-r+j)
                assert admitted(raw,A,t,Ds),(A,t,Ds,P)
                CY=[0,Ds[0],A-Ds[0],A];CZ=[0,Ds[1],A-Ds[1],A]
                if not(any(abs(abs(raw[1])-v)<=4*t for v in CY)or any(abs(abs(raw[2])-v)<=4*t for v in CZ)):
                    COUNTS['positive_interaction_arcs_outside_endpoint_strips']+=1
                COUNTS['whole_positive_interaction_arcs_in_linear_inventory']+=1;count+=1
        cases.append({'A':A,'t':t,'degrees':Ds,'entire_positive_straight_arcs':count})
    return cases


def fan_control(A,t,Ds,p,full=False,Ws=None,intervals=False):
    if Ws is None:Ws=[{q:2**(q*q) for q in range(D+1)} for D in Ds]
    C=interval_compile(t,A,Ds) if intervals else degree_compile(t,A,Ds)
    fan,cells=walk_linear(C)
    data,cols=axis_model(t,A+t,Ws,p)
    if full:
        original=complete_fan(cols,C['n'],p)
        expected={tuple(C['translation'][k]+sum(e[k] for e in B) for k in range(3)):hull(P) for B,P in original.items()}
        assert set(expected)==set(fan) and all(set(expected[v])==set(fan[v]) for v in fan)
        COUNTS['independent_complete_original_degree_fans']+=1
    else:
        for profile,P in list(cells.values())[::max(1,len(cells)//5)][:6]:
            w=centroid(P);B=field_basis(sorted(cols,key=lambda e:(dot(w,e),e)),cols.__getitem__,p)[0]
            assert profile==tuple(C['translation'][k]+sum(e[k] for e in B) for k in range(3))
            COUNTS['independent_large_original_degree_profiles']+=1
    serial=[{'profile':v,'polygon':[[str(x) for x in w] for w in fan[v]]} for v in sorted(fan)]
    for profile,P in cells.values():
        _,_,N=linear_decision(C,centroid(P))
        for left,right in zip(P,P[1:]+P[:1]):
            if any(left[k]==right[k]==0 for k in range(3)):continue
            raw=next(q for q in N if dot(q,left)==dot(q,right)==0)
            assert admitted(raw,A,t,Ds)
            COUNTS['complete_refinement_facets_in_linear_inventory']+=1
    return {'A':A,'t':t,'axis_order_and_degree_specs':[S if isinstance(S,tuple) or S is None else(0,S) for S in Ds],
            'characteristic':p,'actual_initial_ideal_cones':len(fan),
            'refinement_patches':len(cells),'exact_total_area':'1','entire_original_fan_equal':full,
            'canonical_profile_fan_sha256':hashlib.sha256(json.dumps(serial,sort_keys=True).encode()).hexdigest()}


def interval_controls():
    cases=[]
    for A,t,spec,p in [(5,2,(1,3),0),(6,3,(4,5),0),(8,3,(2,5),0),(5,3,(4,4),2),(6,2,(2,3),2),(5,2,None,2)]:
        W={} if spec is None else {q:2**((q-spec[0])**2) if not p else 1 for q in range(spec[0],spec[1]+1)}
        pure,head=axis_factor(W,A,t);fresh=compile_factor(pure,head,A,t,p);predicted=interval_factor(A,t,spec)
        assert fresh==predicted,(A,t,spec,p)
        COUNTS['entire_degree_and_order_tables_equal_field_elimination']+=1
        for theta in [Q(-2*A-1,2),Q(-3,2),Q(3,2),Q(2*A+1,2)]:
            direct=direct_births(pure,head,theta,p)
            B,V,birth=factor(pure,head,theta,p)
            assert direct==[theta*birth[i][0]+birth[i][1] if i in birth else None for i in range(t)]
            COUNTS['independent_nonunit_entire_prefix_births']+=t
        item={'A':A,'t':t,'degree_and_order':spec,'characteristic':p,'operators':{str(q):str(c) for q,c in W.items()}}
        if spec is not None and A<=6:item['all_minors']=all_minors(W,A,spec[1],p,spec[0])
        cases.append(item)
    return cases


def singular_scope_control():
    A,t,D,p=6,2,3,2;W={q:1 for q in range(D+1)};theta=Q(-1,2)
    pure,head=axis_factor(W,A,t);_,_,birth=factor(pure,head,theta,p)
    actual=theta*birth[1][0]+birth[1][1];predicted=Q(3);assert actual==Q(7,2)
    prefix=[v for e,v in pure.items() if theta*e[0]+e[1]<=predicted]+head[:1]
    r=len(field_basis(tuple(range(len(prefix))),prefix.__getitem__,p)[0])
    after=prefix+[head[1]];rr=len(field_basis(tuple(range(len(after))),after.__getitem__,p)[0])
    assert (r,rr)==(5,6)
    minor=[[W.get(r-c,0) for c in [0,1]] for r in [1,2]];assert minor==[[1,1],[1,1]]
    COUNTS['independent_full_support_singular_scope_guards']+=1
    return {'A':A,'t':t,'D':D,'characteristic':p,'coefficients':[1,1,1,1],
            'admissible_minor_rows':[1,2],'admissible_minor_columns':[0,1],'minor':minor,
            'theta':str(theta),'head_index':1,'generic_predicted_birth':str(predicted),'actual_birth':str(actual),
            'entire_original_prefix_rank_before_and_after_head':[r,rr],
            'scope':'Every polynomial coefficient is nonzero, but an admissible minor vanishes. The superregular theorem and degree-only oracle do not apply.'}


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    prior=json.loads((ROOT/'results/research/hook_axis_linear_checks.json').read_text());assert prior['status']=='passed'
    for f,h in prior['source_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
    rng=Random(20261023133);controls=[]
    for A,t,D in [(4,2,0),(4,2,2),(5,2,4),(6,3,3),(8,3,4),(10,4,9)]:
        W={q:2**(q*q) for q in range(D+1)}
        item=field_controls(W,A,t,D,0,A<=6)
        if A<=6:item['all_minors']=all_minors(W,A,D,0)
        controls.append(item)
    W={0:1,1:1};item=field_controls(W,5,2,1,2,True);item['all_minors']=all_minors(W,5,1,2);controls.append(item)
    for A,t,D in [(5,2,2),(6,2,4)]:
        for _ in range(100):
            W={q:rng.randrange(1,65537) for q in range(D+1)}
            try:certificate=all_minors(W,A,D,65537)
            except AssertionError:continue
            break
        else:raise AssertionError('No explicitly certified finite-field fixture found.')
        item=field_controls(W,A,t,D,65537,True);item['all_minors']=certificate;controls.append(item)
    finite=next(item for item in controls if(item['A'],item['D'],item['characteristic'])==(5,2,65537))
    WF={int(q):int(c) for q,c in finite['operators'].items()}
    arcs=arc_controls();print('Degree flags and arc inventories passed',flush=True)
    fans=[fan_control(4,2,(2,3),0,True),fan_control(5,2,(2,2),65537,True,[WF,WF]),
          fan_control(5,2,(2,2),0,True),
          fan_control(8,3,(4,4),0),fan_control(12,4,(6,6),0)]
    assert fans[1]['canonical_profile_fan_sha256']==fans[2]['canonical_profile_fan_sha256']
    COUNTS['entire_same_degree_profile_fan_equalities_across_fields']+=1
    print('Generic degree fans passed',fans,flush=True)
    interval=interval_controls()
    more=[fan_control(5,2,((1,3),(1,3)),0,True,[{q:2**((q-1)**2) for q in range(1,4)}]*2,True),
          fan_control(5,3,((4,4),(4,4)),2,True,[{4:1}]*2,True),
          fan_control(5,2,(None,(0,2)),0,True,[{},{q:2**(q*q) for q in range(3)}],True)]
    fans.extend(more);guard=singular_scope_control()
    files={'docs/hook_superregular_axis_fans.md','scripts/run_hook_superregular_axis_checks.py',
           'docs/hook_axis_birth_curves.md','docs/hook_axis_linear_oracle.md','results/research/hook_axis_linear_checks.json'}
    pending=['scripts/run_hook_superregular_axis_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=([node.module] if isinstance(node,ast.ImportFrom) else[v.name for v in node.names] if isinstance(node,ast.Import) else[])
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and(ROOT/target).exists() and target not in files:files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_superregular_axis_checks.py',
            'uniform_result':'Zero-outer polynomial axis hooks with banded-superregular Toeplitz factors have g=O(n^2), uniformly in all degrees, initial orders, coefficients and fields satisfying the hypotheses. Zero arms are allowed.',
            'coefficient_open_scope':'The banded-superregular locus is geometrically nonempty in every characteristic, with points over every infinite field or sufficiently large finite extensions. Full support or a fixed small prime field is not asserted sufficient.',
            'consequences':{'common_integer_list':'O(n^2) weights of height O(a^2), for fixed a,t and both initial orders/degrees',
                            'degree_and_order_only_profiles':'Pure prefix bases and entire head births depend only on integer degrees and initial orders; no field elimination is required after the hypotheses are known.',
                            'scope_guard':'No efficient certification of all minors, optimal basis-output or bit-complexity claim.'},
            'counts':dict(COUNTS),'factor_controls':controls,'nonunit_and_zero_controls':interval,'positive_interaction_arc_controls':arcs,'complete_degree_fans':fans,
            'full_support_singular_counterexample':guard,
            'limits':['Toeplitz-degenerate coefficients and nonzero outer images are outside this theorem.',
                      'Larger original quotient checks are profile spot checks, not independent entire original fans.',
                      'Generic polynomial strata at every degree are covered; no new unrestricted G3 exponent, infinite actual polynomial lower excess or strongest-current novelty is established.',
                      'The original extremal goal and final-paper gates remain incomplete.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_superregular_axis_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}),flush=True)


if __name__=='__main__':main()
