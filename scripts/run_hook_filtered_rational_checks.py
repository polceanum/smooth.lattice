#!/usr/bin/env python3
"""Original-prefix certificates for rational relations with terminal errors."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from random import Random
import ast,hashlib,json

from run_hook_bridge_checks import model,actual,normals,polygon,COUNTS as REPLAY
from run_hook_rank_departure_checks import clean,rank,complete_fan,plane
from run_hook_coherent_unit_window_checks import scalar,multiply,inverse,transform,lower,upper
from run_hook_axis_compilation_checks import axis_model,axis_factor,compile_tables,hull
from run_hook_axis_linear_checks import incremental_factor,walk_linear
from run_signed_face_checks import field_basis,field_reduce
from run_structural_complexity_checks import dot,area2

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def add(P,E,p):return clean({v:P.get(v,0)+E.get(v,0) for v in set(P)|set(E)},p)


def specification(P,Qpoly,A,t,h,p,rng):
    P={v:scalar(c,p) for v,c in P.items() if scalar(c,p)}
    Qpoly={v:scalar(c,p) for v,c in Qpoly.items() if scalar(c,p)}
    dp=max(i+j for i,j in P);dq=max(i+j for i,j in Qpoly)
    assert dp<A and dq<A and P.get((dp,0)) and Qpoly.get((dq,0)) and Qpoly.get((0,0))
    E={(i,j):scalar(rng.randrange(p or 7),p) for i,j in product(range(A-h,A),range(t))}
    E={v:c for v,c in E.items() if c}
    if h:E[A-h,0]=scalar(1,p)
    total=clean({v:P.get(v,0)+E.get(v,0) for v in set(P)|set(E)},p)
    W=multiply(total,inverse(Qpoly,A,t,p),A,t,p)
    product_QW=multiply(Qpoly,W,A,t,p)
    residual=clean({v:product_QW.get(v,0)-P.get(v,0) for v in set(P)|set(product_QW)},p)
    assert residual==E and all(i>=A-h for i,j in residual)
    if h:assert residual
    return {'P':P,'Q':Qpoly,'W':W,'E':E,'p':dp,'r':dq,'h':h}


def build(t,A,p,specs,rng):
    a=A+t;data=model(t,a);co=[]
    for g,y,z in data[2]:
        arm=0 if y>=t else 1;q,j=(y-t,z) if not arm else(z-t,y)
        if g==1:c=specs[arm]['W'].get((q,j),0)
        elif ((g==0 and not arm)or(g==2 and arm)) and(q,j)==(0,0):c=1
        else:c=rng.randrange(p or 7)
        co.append(scalar(c,p))
    cols=actual(data,co,p,True)
    charts=[transform(data,cols,t,a,arm,S['Q'],p) for arm,S in enumerate(specs)]
    for arm,S in enumerate(specs):
        for q,j in product(range(A),range(t)):
            e=lower(t,arm,q,j)
            got={(data[0][r][arm]-t,data[0][r][1-arm]):c for r,c in charts[arm].get(e,{}).items() if data[0][r][arm]>=t}
            expected=multiply(add(S['P'],S['E'],p),{(q,j):1},A,t,p)
            assert got==expected
            COUNTS['original_denominator_chart_columns']+=1
        m=max(t,S['h'])
        for q in range(min(m,A)):
            error=multiply(S['E'],{(q,0):1},A,t,p)
            if error:
                COUNTS['nonzero_early_denominator_errors']+=1
                if q>=t:COUNTS['nonzero_early_pure_denominator_errors']+=1
        for w in [(31,7,1),(19,1,3)]:
            B=field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)[0]
            BB=field_basis(sorted(charts[arm],key=lambda v:(dot(w,v),v)),charts[arm].__getitem__,p)[0]
            assert B==BB and normals(B,cols,p)==normals(BB,charts[arm],p)
            COUNTS['whole_original_normal_cones_preserved_by_chart']+=1
    return data,cols,charts


def combination(terms,cols,p):
    out=Counter()
    for e,c in terms.items():
        for r,v in cols.get(e,{}).items():out[r]+=c*v
    return clean(out,p)


def pure_relation(t,A,arm,S,k,j,p):
    terms=Counter()
    for (i,h),c in S['Q'].items():terms[lower(t,arm,k+i,j+h)]+=c
    for (i,h),c in S['P'].items():terms[upper(t,arm,k+i,j+h)]-=c
    return clean(terms,p)


def wgt(lam,u,arm):return (lam,u,Q(1)) if not arm else(lam,Q(1),u)


def identities_and_loops(data,cols,t,A,arm,S,p):
    m=max(t,S['h']);a=A+t;identities=0;prefix_targets=[]
    for k,j in product(range(m,A),range(t)):
        terms=pure_relation(t,A,arm,S,k,j,p)
        assert not combination(terms,cols,p);identities+=1
    COUNTS['exact_filtered_pure_source_identities']+=identities
    for kind in ['upper','lower']:
        degree=S['p'] if kind=='upper' else S['r']
        for q in sorted({m+degree,A-1}):
            if q>=A:continue
            e=upper(t,arm,q,0) if kind=='upper' else lower(t,arm,q,0)
            terms=pure_relation(t,A,arm,S,q-degree,0,p);assert terms.get(e)
            for u in [Q(1001,1000),Q(3,2),Q(3)]:
                plane_value=t+(S['r']-S['p'])*u
                if kind=='upper':lam=max(Q(1),plane_value+1)
                elif plane_value>0:lam=plane_value/2
                else:continue
                w=wgt(lam,u,arm)
                assert all(dot(w,f)<dot(w,e) for f,c in terms.items() if f!=e and cols.get(f))
                COUNTS['strict_original_pure_loop_certificates']+=1
                prefix_targets.append((e,w))
    head=(0,a,0) if not arm else(0,0,a)
    F={(data[0][r][arm]-t,data[0][r][1-arm]):c for r,c in cols[head].items() if data[0][r][arm]>=t}
    QF=multiply(S['Q'],F,A,t,p)
    start=max(2*t,m+S['p']);outer_ids=0
    for i,j in product(range(start,A),range(t)):
        def outer(q,v):return (0,a+q,v) if not arm else(0,v,a+q)
        terms=Counter();k=i-S['p'];e=outer(i,j)
        for (r,s),c in S['P'].items():terms[outer(k+r,j+s)]+=c
        for (r,s),c in QF.items():terms[lower(t,arm,k+r,j+s)]-=c
        terms=clean(terms,p);assert terms.get(e)==S['P'][S['p'],0]
        assert not combination(terms,cols,p);outer_ids+=1
        for u in [Q(1001,1000),Q(2)]:
            w=wgt(Q(31),u,arm)
            assert all(dot(w,f)<dot(w,e) for f,c in terms.items() if f!=e and cols.get(f))
        if j==0 and i in {start,A-1}:prefix_targets.append((e,wgt(Q(31),Q(1001,1000),arm)))
    COUNTS['exact_filtered_arbitrary_outer_identities']+=outer_ids
    for e,w in prefix_targets[:3]+prefix_targets[-3:]:
        cheap=tuple(v for v in cols if dot(w,v)<dot(w,e))
        piv=field_basis(cheap,cols.__getitem__,p)[1]
        assert not field_reduce(cols[e],piv,p)[0]
        COUNTS['independent_entire_cheaper_prefix_loop_ranks']+=1
    return {'m':m,'pure_identities':identities,'outer_threshold':start,'outer_identities':outer_ids}


def separating_rows(data,cols,chart,t,A,arm,S,p):
    a=A+t;m=max(t,S['h']);checks=[]
    for kind in ['upper','lower']:
        shift=S['r'] if kind=='upper' else S['p']
        cutoff=S['p'] if kind=='upper' else S['r']
        inputs=[q for q in range(m+cutoff,A) if q+shift<A-m]
        for q in inputs[::max(1,len(inputs)//2)][:2]:
            for u in [Q(1001,1000),Q(3,2),Q(3)]:
                pv=t+(S['r']-S['p'])*u
                if kind=='upper':
                    if pv<=0:continue
                    lam=pv/2
                else:lam=max(Q(1),pv+1)
                w=wgt(lam,u,arm)
                e=upper(t,arm,q,0) if kind=='upper' else lower(t,arm,q,0)
                tau=dot(w,e)
                if tau>=a*u:continue
                row=data[0].index((t+q+shift,0) if not arm else(0,t+q+shift))
                assert chart[e].get(row)
                prefix=[v for v in cols if dot(w,v)<=tau]
                assert [v for v in prefix if chart[v].get(row)]==[e]
                cheap=tuple(v for v in cols if dot(w,v)<tau)
                assert rank((*cheap,e),cols,p)==rank(cheap,cols,p)+1
                COUNTS['independent_entire_prefix_unique_separating_rows']+=1
                checks.append({'type':kind,'input':q,'weight':list(map(str,w)),
                               'separating_local_row':q+shift,'prefix_size':len(prefix)})
    return checks


def full_fan_control(data,cols,t,A,specs,p):
    n=4*(A+t)*t-2*t*t;fan=complete_fan(cols,n,p)
    assert sum(area2(P) for P in fan.values())==1
    a=A+t
    d=max(max(S['p'],S['r']) for S in specs);h=max(S['h'] for S in specs);width=6*t+2*d+2*h
    exceptional={plane((1,-(S['r']-S['p']),-t) if arm==0 else(1,-t,-(S['r']-S['p']))) for arm,S in enumerate(specs)}
    exceptional.add((0,1,-1));facets=0
    # Test the original labelled differences before taking primitive normals.
    # Dividing an x=0 normal by its gcd need not preserve an endpoint strip.
    Y,Z=map(set,data[4][:2]);labels=tuple(cols);inventory=set(exceptional)
    for pos,e in enumerate(labels):
        for f in labels[:pos]:
            if (e in Y and f in Z)or(e in Z and f in Y):
                y,z=(e,f) if e in Y else(f,e)
                b=y[1]-z[1];c=z[2]-y[2]
                if not any(min(abs(v),abs(v-a))<=width for v in(b,c)):continue
            inventory.add(plane(tuple(x-y for x,y in zip(e,f))))
    for B,P in fan.items():
        N=normals(B,cols,p)
        for x,y in zip(P,P[1:]+P[:1]):
            if any(x[i]==y[i]==0 for i in range(3)):continue
            q=plane(next(q for q in N if dot(q,x)==dot(q,y)==0));facets+=1
            assert q in inventory
    COUNTS['independent_complete_original_quotient_fans']+=1
    COUNTS['original_normal_facets_in_filtered_inventory']+=facets
    serial=[{'basis':sorted(B),'polygon':[[str(c) for c in v] for v in hull(P)]} for B,P in sorted(fan.items(),key=lambda z:sorted(z[0]))]
    return {'complete_original_cones':len(fan),'exact_total_area':'1','facets_checked':facets,
            'canonical_fan_sha256':hashlib.sha256(json.dumps(serial,sort_keys=True).encode()).hexdigest()}


def axis_calibrations():
    cases=[]
    for t,A in [(2,40),(4,13),(6,32)]:
        p=1009;nu=t+2
        Ws=[{0:1,nu:2,A-3:3,A-1:5},{0:1,nu+1:7,A-4:11,A-2:13}]
        h=4;low=[{q:c for q,c in W.items() if q<A-h} for W in Ws]
        d=max(max(W) for W in low)
        assert d<=3*t
        compiled=compile_tables(t,A+t,[axis_factor(W,A,t) for W in Ws],p)
        for arm,W in enumerate(Ws):compiled['factors'][arm]=incremental_factor(*axis_factor(W,A,t),A,t,p)[0]
        fan,cells=walk_linear(compiled);data,cols=axis_model(t,A+t,Ws,p)
        selected=list(cells.values())[::max(1,len(cells)//5)][:6]
        for profile,P in selected:
            w=tuple(sum(v[k] for v in P)/len(P) for k in range(3))
            B=field_basis(sorted(cols,key=lambda e:(dot(w,e),e)),cols.__getitem__,p)[0]
            assert profile==tuple(compiled['translation'][k]+sum(e[k] for e in B) for k in range(3))
            COUNTS['independent_calibration_original_profiles']+=1
        COUNTS['complete_separated_axis_calibration_fans']+=1
        cases.append({'t':t,'A':A,'h':h,'compact_low_degree':d,'operators':Ws,
                      'refinement_patches':len(cells),'actual_initial_ideal_cones':len(fan),'exact_total_area':'1',
                      'scope':'Complete compiled cover and independent original profile spot checks. Not an independently enumerated entire original fan or an asymptotic rate.'})
        print('Filtered axis calibration passed',t,A,len(fan),flush=True)
    return cases


def filtered_degree_gate(rng):
    sharp=[]
    for p in [0,2,3,101]:
        for A in [6,8,12]:
            for D in range(1,A):
                W={i:scalar(rng.randrange(p or 5),p) for i in range(D+1)}
                W={i:c for i,c in W.items() if c};W[0]=W[D]=scalar(1,p)
                bound=min(D,A-D)
                for d in range(bound):
                    for h in range(bound-d):
                        # Coefficients of QW above degree d and below A-h must vanish.
                        # Fix Q(0)=1. These equations would be necessary even allowing
                        # mixed P,Q; setting v=0 preserves their pure leading terms.
                        rows=range(d+1,A-h)
                        cols={j:{i:W.get(i-j,0) for i in rows if W.get(i-j,0)} for j in range(1,d+1)}
                        rhs={i:scalar(-W.get(i,0),p) for i in rows if W.get(i,0)}
                        piv=field_basis(tuple(cols),cols.__getitem__,p)[1]
                        assert field_reduce(rhs,piv,p)[0]
                        COUNTS['independent_forbidden_filtered_degree_systems']+=1
                # Two exact binomial representations attain the two sides of the bound.
                assert min(D+(0),0+(A-D))==bound
                binomial={(0,0):scalar(1,p),(D,0):scalar(1,p)}
                assert multiply({(0,0):1},binomial,A,2,p)==binomial
                error=clean({v:c-(1 if v==(0,0) else 0) for v,c in binomial.items()},p)
                assert error=={(D,0):scalar(1,p)} and D==A-(A-D)
                COUNTS['exact_sharp_binomial_filtered_degree_certificates']+=2
                sharp.append({'A':A,'degree':D,'characteristic':p,'minimum_d_plus_h':bound})
    return {'uniform_gate':'d+h >= min(D,A-D), including mixed long-leading P,Q by setting v=0',
            'sharp_binomial_controls':sharp,
            'scope':'A parameter lower bound, not a fan lower bound. Length-two pure-axis binomial chains have only three constant-size component types.'}


def binomial_component_controls():
    out=[]
    for t,D,p in [(2,6,101),(3,9,2),(4,12,101)]:
        A=2*D;Ws=[{0:1,D:1}]*2
        compiled=compile_tables(t,A+t,[axis_factor(W,A,t) for W in Ws],p)
        for arm,W in enumerate(Ws):compiled['factors'][arm]=incremental_factor(*axis_factor(W,A,t),A,t,p)[0]
        fan,cells=walk_linear(compiled);data,cols=axis_model(t,A+t,Ws,p)
        # Independently form actual row-support components and normalize every
        # original column price by the component's short-root translation.
        remaining={r for v in cols.values() for r in v};types=Counter()
        while remaining:
            rows={min(remaining)};labels=set();changed=True
            while changed:
                changed=False
                for e,v in cols.items():
                    if e not in labels and rows.intersection(v):
                        labels.add(e);before=len(rows);rows.update(v);changed|=len(rows)!=before
            remaining.difference_update(rows)
            assert len(labels) in [4,7]
            # Pure Y roots have original offset (q0,j), pure Z (i,q0).
            # A joined component uses (i,j), common to its Y and Z roots.
            yrows=[data[0][r] for r in rows if data[0][r][0]>=t]
            zrows=[data[0][r] for r in rows if data[0][r][1]>=t]
            if yrows and zrows:offset=(min(y for y,z in yrows)-t,min(z for y,z in yrows))
            elif yrows:offset=(min(y for y,z in yrows)-t,min(z for y,z in yrows))
            else:offset=(min(y for y,z in zrows),min(z for y,z in zrows)-t)
            normalized=tuple(sorted(((e[0],e[1]-offset[0],e[2]-offset[1]),
                                      tuple(sorted(((data[0][r][0]-offset[0],data[0][r][1]-offset[1]),str(c)) for r,c in cols[e].items())))
                                     for e in labels))
            types[normalized]+=1
            COUNTS['independent_actual_length_two_row_components']+=1
        assert len(types)==3
        COUNTS['constant_component_type_binomial_controls']+=1
        out.append({'t':t,'A':A,'D':D,'characteristic':p,'filtered_degree_lower_bound':D,
                    'component_types':len(types),'type_sizes':sorted(len(T) for T in types),
                    'actual_profile_cones':len(fan),'complete_compiled_area':'1',
                    'scope':'Actual row components independently checked. The constant-type proof gives O(1), independently of the growing filtered-degree parameter.'})
    return out


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    prior=json.loads((ROOT/'results/research/hook_axis_linear_checks.json').read_text());assert prior['status']=='passed'
    for f,h in prior['source_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
    rng=Random(20261023128);cases=[]
    configs=[(2,6,101,1,False,True),(2,7,2,2,True,True),
             (2,28,101,3,False,False),(3,28,3,2,True,False),(2,12,0,1,True,False)]
    for t,A,p,h,nonunit,full in configs:
        P={(2,0):1,(0,1):1};Qpoly={(0,0):1,(1,0):1,(0,1):1}
        if not nonunit:P[0,0]=1
        specs=[specification(P,Qpoly,A,t,h,p,rng) for _ in range(2)]
        data,cols,charts=build(t,A,p,specs,rng)
        identities=[identities_and_loops(data,cols,t,A,arm,S,p) for arm,S in enumerate(specs)]
        rows=[separating_rows(data,cols,charts[arm],t,A,arm,S,p) for arm,S in enumerate(specs)]
        full_report=full_fan_control(data,cols,t,A,specs,p) if full else None
        cases.append({'t':t,'A':A,'characteristic':p,'terminal_width':h,'nonunit_numerator':nonunit,
                      'specifications':[{key:[{'q':i,'v':j,'coefficient':str(c)} for(i,j),c in sorted(S[key].items())]
                                         for key in ['P','Q','E']} for S in specs],
                      'arbitrary_outer_images':True,'filtered_identity_checks':identities,
                      'whole_prefix_separating_rows':rows,'complete_original_fan':full_report})
        print('Filtered rational quotient passed',t,A,p,'separating rows',sum(map(len,rows)),flush=True)
    assert COUNTS['nonzero_early_denominator_errors'] and COUNTS['nonzero_early_pure_denominator_errors']
    assert COUNTS['independent_entire_prefix_unique_separating_rows']
    axis=axis_calibrations();degree_gate=filtered_degree_gate(rng);binomial=binomial_component_controls()
    files={'docs/hook_filtered_rational_windows.md','scripts/run_hook_filtered_rational_checks.py',
           'docs/hook_rational_leading_windows.md','docs/hook_axis_linear_oracle.md','results/research/hook_axis_linear_checks.json'}
    pending=['scripts/run_hook_filtered_rational_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=([node.module] if isinstance(node,ast.ImportFrom) else[v.name for v in node.names] if isinstance(node,ast.Import) else[])
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and(ROOT/target).exists() and target not in files:files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_filtered_rational_checks.py',
            'uniform_result':'If each arm has a long-leading rational relation QW-P supported in its last h long rows, with compact degree d, then g=O(a^2*(t+d+H)^2), all fields and arbitrary compatible outer images. A common complete integer list has the same size and height O(a^2).',
            'mechanism':'Shifts k>=max(t,h) annihilate the residual exactly. Below long row A-max(t,h), transformed terminal/head errors and opposite outer bridges cannot hit the separating row.',
            'counts':dict(COUNTS),'original_quotient_replay_counts':dict(REPLAY),'quotients':cases,'separated_axis_calibrations':axis,
            'filtered_degree_scope_gate':degree_gate,
            'constant_component_type_scope_controls':binomial,
            'limits':['The full unrestricted G_3 endpoint, actual infinite polynomial lower excess, exact strongest-current novelty and original-goal final paper remain open.',
                      'The old nearby and truncated-dense probes were already covered by the compact exact rational theorem. The separated probes are now quadratically bounded by the terminal-error allowance.',
                      'Finite complete covers and profile spot checks are reported separately from independent entire original fans.',
                      'Large filtered degree relative to thickness and lack of a pure long leading representation remain outside the quadratic class.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_filtered_rational_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}),flush=True)


if __name__=='__main__':main()
