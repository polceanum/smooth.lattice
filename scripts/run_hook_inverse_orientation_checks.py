#!/usr/bin/env python3
"""Lower-column windows controlled only by the chosen central inverse."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from random import Random
import ast,hashlib,json

from run_hook_bridge_checks import model,actual,normals,polygon,COUNTS as REPLAY
from run_hook_rank_departure_checks import complete_fan,plane,rank,weight_controls,COUNTS as FAN
from run_hook_coherent_unit_window_checks import inverse,multiply,transform,lower
from run_hook_outer_unit_window_checks import admitted
from run_signed_face_checks import field_basis,field_reduce
from run_structural_complexity_checks import dot,area2
from run_oracle_enumeration_checks import exact_adjacent_weight

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def build(t,a,spec,p,rng,zero=False):
    data=model(t,a);A=a-t;invs=[];operators=[]
    for S in spec:
        if S is None:invs.append(None);operators.append(None);continue
        S={v:(c%p if p else Q(c)) for v,c in S.items() if v[0]<A and v[1]<t and (c%p if p else c)}
        W=inverse(S,A,t,p);assert multiply(W,S,A,t,p)=={(0,0):1}
        invs.append(S);operators.append(W)
    co=[]
    for g,y,z in data[2]:
        arm=0 if y>=t else 1;q,j=(y-t,z) if arm==0 else(z-t,y)
        if g==1 and operators[arm] is not None:c=operators[arm].get((q,j),0)
        elif g==1:c=0 if (q,j)==(0,0) else rng.randrange(p or 7)
        else:c=0 if zero else rng.randrange(p or 7)
        co.append(c)
    cols=actual(data,co,p,True);charts=[]
    for arm,S in enumerate(invs):
        if S is None:charts.append(None);continue
        mat=transform(data,cols,t,a,arm,S,p)
        for q,j in product(range(A),range(t)):
            label=lower(t,arm,q,j)
            pair=(t+q,j) if arm==0 else(j,t+q);row=data[0].index(pair)
            own={i:c for i,c in mat[label].items() if data[0][i][arm]>=t}
            assert own=={row:1}
            if q>=t:assert own==mat[label]
            COUNTS['independent_original_lower_identity_chart_checks']+=1
        charts.append(mat)
    return data,co,cols,invs,operators,charts


def weight(lam,u,arm):return (lam,u,Q(1)) if arm==0 else(lam,Q(1),u)


def prefix_controls(data,cols,chart,t,a,arm,d,zero,p):
    A=a-t;states=[]
    for u in (Q(1,2),Q(1),Q(3)):
        M=a-1 if zero else a//2;j=0;q=M-t;e=lower(t,arm,q,j)
        w=weight(Q(7,3),u,arm);tau=dot(w,e)
        if zero:
            if q*u+j<=(t-1+d)*u+t-1:continue
        else:
            B=4*t+2*d
            if not (u*M>B*u+4*t and u*(a-M)>B*u+4*t):continue
            assert tau<a*u
        removed=set();labels=[]
        for i,xy in enumerate(data[0]):
            if xy[arm]<2*t:continue
            v=lower(t,arm,xy[arm]-t,xy[1-arm])
            if dot(w,v)<tau:
                assert chart[v]=={i:1};removed.add(i);labels.append(v)
        projected={v:{i:c for i,c in col.items() if i not in removed} for v,col in chart.items()}
        inside={i for i,xy in enumerate(data[0]) if xy[arm]>=2*t and (zero or xy[arm]<a-t)}
        relevant=projected if zero else {v:c for v,c in projected.items() if dot(w,v)<=tau}
        assert projected[e] and set(projected[e])<=inside
        assert all(not (set(c)&inside) or set(c)<=inside for c in relevant.values())
        if not states:
            prefix=tuple(v for v in chart if dot(w,v)<=tau)
            assert rank(prefix,chart,p)==len(removed)+rank(prefix,projected,p)
            COUNTS['independent_original_prefix_contraction_rank_equalities']+=1
        higher_mix=sum(bool(set(c)&inside) and not set(c)<=inside for v,c in projected.items() if dot(w,v)>tau)
        if not zero:COUNTS['higher_price_mixed_columns_not_claimed_to_split']+=higher_mix
        COUNTS['independent_entire_inverse_only_lower_splits']+=1
        states.append({'arm':arm,'zero_outer':zero,'weight':list(map(str,w)),'target':e,'contracted_rows':len(removed),'checked_entire_matrix':zero,'higher_price_mixed_columns':higher_mix})
    return states


def directed_inventory(data,t,a,d,arm,zero):
    out=set()
    for e in data[4][arm]:
        if e[0] or (zero and e[arm+1]>=a):continue
        for f in data[4][1-arm]:
            if not f[0]:continue
            b=e[arm+1]-f[arm+1];c=f[2-arm]-e[2-arm]
            ok=(b<=6*t+2*d or c<=6*t+2*d) if zero else admitted(t,a,d,b,c)
            if ok:out.add(plane(tuple(x-y for x,y in zip(f,e))))
    H=1+(a//t).bit_length()
    assert len(out)<=1024*a*(d+t*(1 if zero else H))
    COUNTS['coefficient_independent_oriented_plane_inventories']+=1
    return out


def fan_controls(data,cols,invs,t,a,p,zero):
    fan=complete_fan(cols,4*a*t-2*t*t,p)
    ds=[max(i for i,j in S) if S else None for S in invs]
    inventories=[directed_inventory(data,t,a,d,arm,zero) if d is not None else None for arm,d in enumerate(ds)]
    known={plane(v) for v in data[7]};checks=[0,0];uncontrolled=0;rempieces=0
    for B,poly in fan.items():
        mids=[];N=normals(B,cols,p)
        for x,y in zip(poly,poly[1:]+poly[:1]):
            if any(x[i]==y[i]==0 for i in range(3)):continue
            q=plane(next(q for q in N if dot(q,x)==dot(q,y)==0))
            if q[0] and all(S is not None for S in invs):assert q in known|inventories[0]|inventories[1]
            if q[0] and q not in known:rempieces+=1
            mids.append(tuple((s+r)/2 for s,r in zip(x,y)))
        piv=field_basis(tuple(B),cols.__getitem__,p)[1]
        for f in cols:
            if f in B:continue
            rem,rep=field_reduce(cols[f],piv,p);assert not rem
            for e,c0 in rep.items():
                if not c0 or e[0]==f[0]:continue
                if not ((e in data[4][0] and f in data[4][1])or(e in data[4][1] and f in data[4][0])):continue
                lo,up=(e,f) if not e[0] else(f,e);arm=0 if lo in data[4][0] else 1
                for w in mids:
                    if dot(w,e)!=dot(w,f):continue
                    if invs[arm] is None:uncontrolled+=1;continue
                    d=ds[arm];u=w[arm+1]/w[2-arm];M=lo[arm+1];s=lo[2-arm];K=up[2-arm]
                    if M>=a:
                        assert not zero and (M-a+1)*u+s<=2*t-1
                    elif zero:
                        assert (M-t)*u+s-t<=(t-1+d)*u+t-1
                        assert M<=3*t+d or K<=4*t+d
                    else:
                        B0=4*t+2*d
                        assert u*M<=B0*u+4*t or u*(a-M)<=B0*u+4*t
                    b=M-up[arm+1];c=K-s
                    assert (b<=6*t+2*d or c<=6*t+2*d) if zero else admitted(t,a,d,b,c)
                    assert plane(tuple(x-y for x,y in zip(up,lo))) in inventories[arm]
                    checks[arm]+=1;COUNTS['independent_actual_oriented_coherent_exchange_checks']+=1
    COUNTS['independent_complete_nonzero_wall_inventory_checks']+=len(fan)
    return {'complete_original_cones':len(fan),'oriented_exchange_checks':checks,'opposite_nonunit_orientation_not_claimed':uncontrolled,'directed_planes':[len(P) if P is not None else None for P in inventories],'known_planes':len(known),'remainder_wall_occurrences_with_both_incident_cones':rempieces}


def remainder_crossing():
    """Replay four entire facets; other walls may meet at their common vertex."""
    t,a,p=3,16,1009;n=4*a*t-2*t*t;data=model(t,a)
    Ws=[{(0,0):1,(3,1):771,(4,1):191,(6,2):985,(7,2):544,(8,2):815},
        {(0,0):1,(3,1):730,(4,1):221,(6,2):430,(7,2):986,(8,2):570}]
    co=[]
    for g,y,z in data[2]:
        arm=0 if y>=t else 1;q,j=(y-t,z) if arm==0 else(z-t,y)
        co.append(Ws[arm].get((q,j),0) if g==1 else 0)
    cols=actual(data,co,p,True);v=(Q(3,5),Q(1,5),Q(1,5))
    known={plane(q) for q in data[7]};segments=[]
    for W in Ws:
        assert all(3*j<=i<=4*j for i,j in W)
        assert all(3*j<=i<=4*j for i,j in inverse(W,a-t,t,p))
    star_basis=field_basis(sorted(cols,key=lambda x:(dot(v,x),x)),cols.__getitem__,p)[0]
    optimum=sum(dot(v,x) for x in star_basis)
    fixtures=[((1,-6,3),(12,5,6)),((1,-6,3),(24,7,6)),
              ((1,3,-6),(12,6,5)),((1,3,-6),(24,6,7))]
    for q,endpoint in fixtures:
        assert q not in known
        end=tuple(Q(x,sum(endpoint)) for x in endpoint)
        assert min(end)>0 and dot(q,v)==dot(q,end)==0
        midpoint=tuple((x+y)/2 for x,y in zip(v,end));bases=[]
        for sign in (-1,1):
            w,_=exact_adjacent_weight(n,v,end,tuple(sign*x for x in q))
            B=frozenset(field_basis(sorted(cols,key=lambda x:(dot(w,x),x)),cols.__getitem__,p)[0])
            N=normals(B,cols,p);P=polygon(N)
            assert all(dot(w,x)>0 for x in N) and v in P and end in P
            assert sum(dot(v,x) for x in B)==optimum
            assert all(dot(x,midpoint)>=0 for x in N)
            bases.append(B);COUNTS['independent_crossing_entire_adjacent_normal_polygons']+=1
        assert len(bases[0]^bases[1])==2
        e=next(iter(bases[0]-bases[1]));f=next(iter(bases[1]-bases[0]))
        assert e[0]!=f[0] and plane(tuple(x-y for x,y in zip(e,f)))==q
        assert (e in data[4][0] and f in data[4][1])or(e in data[4][1] and f in data[4][0])
        tau=dot(midpoint,e);assert dot(midpoint,f)==tau
        cheap=tuple(x for x in cols if dot(midpoint,x)<tau);r=rank(cheap,cols,p)
        increments=[rank((*cheap,*extra),cols,p)-r for extra in [(e,),(f,),(e,f)]]
        assert r==52 and increments==[1,1,1]
        COUNTS['independent_actual_remainder_crossing_segment_prefix_ranks']+=1
        segments.append({'primitive_plane':q,'other_endpoint':list(map(str,end)),
                         'midpoint':list(map(str,midpoint)),'exchange':sorted((e,f)),
                         'strictly_cheaper_rank':r,'single_and_joint_increments':increments,
                         'adjacent_bases_differ_by_one_pair':True})
    for i in (0,2):
        left=tuple(Q(x) for x in segments[i]['other_endpoint'])
        right=tuple(Q(x) for x in segments[i+1]['other_endpoint'])
        assert any((x-z)*(y-z)<0 for x,y,z in zip(left,right,v))
    q,r=fixtures[0][0],fixtures[2][0]
    assert any(q[i]*r[j]!=q[j]*r[i] for i in range(3) for j in range(i))
    return {'t':t,'a':a,'characteristic':p,'zero_outer_images':True,
            'operators':[[[i,j,c] for (i,j),c in sorted(W.items())] for W in Ws],
            'common_interior_vertex':list(map(str,v)),'segments':segments,
            'scope':'Both remainder directions continue through a common interior vertex. Other walls meet there; neither an isolated four-cone crossing nor an asymptotic lower count is asserted.'}


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    old=json.loads((ROOT/'results/research/hook_two_orbit_checks.json').read_text());assert old['status']=='passed'
    for f,h in old['source_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
    rng=Random(2026102342);cases=[]
    def mixed(t):return {(0,0):1,(t,1):-1,(t+1,1):-1}
    def axis(t):return {(0,0):1,(1,0):-1,(t,1):1}
    configs=[(2,8,(mixed(2),mixed(2)),101,False,True),
             (3,11,(mixed(3),axis(3)),2,False,True),
             (2,8,(axis(2),mixed(2)),0,False,True),
             (3,9,(mixed(3),None),3,False,True),
             (3,12,(axis(3),axis(3)),101,True,True),
             (4,25,(mixed(4),mixed(4)),1009,True,True),
             (2,80,(axis(2),None),101,False,False),
             (3,100,(axis(3),mixed(3)),3,False,False),
             (3,50,(axis(3),axis(3)),2,True,False)]
    for t,a,spec,p,zero,full in configs:
        data,co,cols,invs,Ws,charts=build(t,a,spec,p,rng,zero);arms=[]
        for arm,S in enumerate(invs):
            if S is None:arms.append(None);continue
            d=max(i for i,j in S);forward=max(i for i,j in Ws[arm]);patches=prefix_controls(data,cols,charts[arm],t,a,arm,d,zero,p)
            if (1,0) in S:
                assert forward==a-t-1;COUNTS['independent_long_forward_small_inverse_examples']+=1
            arms.append({'inverse_long_degree':d,'expanded_forward_long_degree':forward,'prefix_split_certificates':patches})
        fan=fan_controls(data,cols,invs,t,a,p,zero) if full else None
        cases.append({'t':t,'a':a,'characteristic':p,'zero_outer_images':zero,'specified_inverse_arms':arms,'full_fan':fan})
        print('Inverse orientation passed',t,a,p,'zero_outer',zero,'cones',fan['complete_original_cones'] if fan else None,flush=True)
    crossing=remainder_crossing()
    files={'docs/hook_inverse_orientation.md','scripts/run_hook_inverse_orientation_checks.py','docs/hook_coherent_unit_windows.md','docs/hook_outer_unit_windows.md','docs/hook_bridge_reduction.md','docs/layer_shadow_budget.md','docs/cross_layer_direction_budget.md','results/research/hook_two_orbit_checks.json'}
    pending=['scripts/run_hook_inverse_orientation_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=([node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else [])
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_inverse_orientation_checks.py',
            'scope':'All-field inverse-only oriented lower-source windows with arbitrary compatible outer maps. Both central arms unit give O(n^2*(1+Dinv/t+log(a/t))*log(n+1)) full fans, and O(n^2*(1+Dinv/t)*log(n+1)) with zero outers. No forward-degree hypothesis, unrestricted endpoint or exact novelty claimed.',
            'uniform_result':{'one_arm':'Only the chosen central arm must be unit; its inverse long degree controls lower-on-that-arm versus opposite-upper faces. Opposite central and all outer images are arbitrary.',
                              'arbitrary_outer_lower_window':'u*M<=(4t+2d)*u+4t or u*(a-M)<=(4t+2d)*u+4t, d=deg_q W^-1; outer sources satisfy (i+1)*u+j<=2t-1',
                              'zero_outer_lower_window':'(M-t)*u+(s-t)<=(t-1+d)*u+t-1; hence M<=3t+d or opposite upper K<=4t+d',
                              'one_oriented_direction_inventory':'O(a*(d+t*(1+log(a/t)))), height O(a); O(a*(t+d)) if the chosen outer components are zero',
                              'two_unit_full_fan':'O(n^2*(1+Dinv/t+log(a/t))*log(n+1)), all compatible outers; Dinv=max inverse degrees',
                              'two_unit_zero_outer_full_fan':'O(n^2*(1+Dinv/t)*log(n+1))',
                              'geometric_family':'Both W=(1-M)^-1, M=v*q^nu*(1+q), have inverse long degree nu+1; nu=O(t) is quadratic polylogarithmic even with arbitrary compatible outer maps.',
                              'necessary_condition':'A unit-arm family with g>=c*n^(2+epsilon) needs Dinv/t=Omega(n^epsilon/log(n+1)); long forward propagation alone cannot produce a polynomial excess.',
                              'counting_input':'The established translated-layer shadow theorem bounds all parallel pieces in each primitive first-coordinate-one direction by O(n*log(n+1)); the finite polytope sweep counts the full fan.'},
            'counts':dict(sorted(COUNTS.items())),'fan_counts':dict(FAN),'original_quotient_replay_counts':dict(REPLAY),'cases':cases,'actual_remainder_crossing':crossing,
            'limits':['Both central units are required for the whole-fan theorem; one unit controls only its specified lower-column orientation.',
                      'The common plane inventory covers nonzero-first-coordinate directions; no complete coefficient-independent weight list of the displayed size is claimed.',
                      'The all-field geometric/geometric family is closed, but plus/plus and power/geometric models can have large inverse degree and are not closed by this corollary.',
                      'Existing stronger bounds remain available when smaller; logarithmic sharpness is not asserted.',
                      'Exact strongest-current novelty, sharp unrestricted G_3, actual polynomial excess and the original final-paper requirements remain unresolved.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_inverse_orientation_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}),flush=True)


if __name__=='__main__':main()
