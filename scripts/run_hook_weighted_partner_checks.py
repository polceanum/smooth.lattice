#!/usr/bin/env python3
"""Whole joint-prefix frontiers and original mixed-hook partner checks."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import comb
from pathlib import Path
import ast
import hashlib
import json

from run_hook_filtered_annihilator_checks import assemble, bases, dense_nullspace, combination
from run_hook_transverse_prefix_checks import clean, mul, wall_weight
from run_hook_monomial_perturbation_checks import reflected, row_functional
from run_hook_bridge_checks import normals, COUNTS as REPLAY
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import dot

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def val(c,p):
    return c%p if p else Q(c)


def reduce_plain(vector,piv,p):
    out=clean(vector,p)
    for j,row in sorted(piv.items()):
        c=out.get(j,0)
        if c:
            for k,a in row.items():
                out[k]=val(out.get(k,0)-c*a,p)
                if not out[k]:del out[k]
    return out


def add_plain(vector,piv,p):
    row=reduce_plain(vector,piv,p)
    if not row:return False
    j=min(row);inv=pow(int(row[j]),-1,p) if p else 1/Q(row[j])
    piv[j]={k:val(c*inv,p) for k,c in row.items()}
    return True


def own_sectors(t):
    cuts=sorted({Q(b,c) for b in range(1,t) for c in range(1,b+1)})
    return list(zip(cuts,cuts[1:]+[None]))


def projected(data,cols,N,t,a,E,branch,J,p):
    A=a-t;L=a-2*t;Y,Z=E[1:]
    left,right=J;u=(left+right)/2 if right is not None else left+1
    local=(Y-t,Z-t);rho=lambda v:v[0]*u+v[1]
    if branch=='primal':
        if local not in bases(N,A,t,u,p)[-2]:return None
        fixed=tuple(v for v in cols if
          (not v[0] and 2*t<=v[1]<a and t<=v[2]<2*t and dot((0,u,1),v)<dot((0,u,1),E))
          or (v[0] and v[1]>=t and rho((v[1]-t,v[2]))<=rho(local))
          or (not v[0] and t<=v[1]<2*t and t<=v[2]<2*t))
        _,piv=field_basis(fixed,cols.__getitem__,p)
        values={v:field_reduce(c,piv,p)[0] for v,c in cols.items()}
        meta={'fixed':fixed}
    else:
        target=(A-1-local[0],t-1-local[1])
        ground=tuple(v for v in product(range(L),range(t)) if rho(v)<=rho(target))
        maps={v:mul(N,{v:1},L,t,p) for v in ground}
        rank,kernel=dense_nullspace(maps,ground,p)
        assert len(kernel)==len(ground)-rank
        rows=[row_functional(h,0,data,A,t) for h in kernel]
        zrows=[r for r,(y,z) in enumerate(data[0]) if y<t]
        zi={r:i for i,r in enumerate(zrows)}
        values={}
        for v,c in cols.items():
            out={i:val(sum(ll.get(r,0)*x for r,x in c.items()),p)
                 for i,ll in enumerate(rows)}
            out.update({len(rows)+zi[r]:x for r,x in c.items() if r in zi})
            values[v]=clean(out,p)
        meta={'functional_dimension':len(kernel),'target':target}
    if not values[E]:return None
    D=tuple(v for v in cols if v[0] and v[1]<t)
    if branch=='primal':
        pp={}
        for v in D:add_plain(values[v],pp,p)
        assert not reduce_plain(values[E],pp,p)
        COUNTS['primal_original_Z_span_projections']+=1
    C=tuple(v for v in cols if not v[0] and t<=v[1]<2*t and values[v])
    return values,C,D,meta


def frontier(data,cols,N,t,a,E,branch,J,p):
    result=projected(data,cols,N,t,a,E,branch,J,p)
    if result is None:return None
    values,C,D,meta=result;Y,Z=E[1:];A=a-t
    lo=J[0];cap=Q(a,Y-(t if branch=='primal' else 2*t)+1)
    hi=min(J[1],cap) if J[1] is not None else cap
    if lo>=hi:return {'records':[],'partners':set(),'projection':result,'J':J}
    entries=sorted((Q(v[2]-Z,Y-v[1]),v) for v in C)
    cuts={lo,hi}
    cuts.update(x for x,v in entries if lo<x<hi)
    cuts.update(Q(i,j) for i in range(1,A) for j in range(1,t) if lo<Q(i,j)<hi)
    cuts=sorted(cuts);active=[];samples=[];cp={};at=0;last=None;born=False;dead=False
    for left,right in zip(cuts,cuts[1:]):
        u=(left+right)/2
        while at<len(entries) and entries[at][0]<u:
            add_plain(values[entries[at][1]],cp,p);at+=1
        if not reduce_plain(values[E],cp,p):
            dead=True;continue
        assert not dead
        piv=dict(cp);F=None
        for v in sorted(D,key=lambda v:(v[1]*u+v[2],v)):
            if add_plain(values[v],piv,p) and not reduce_plain(values[E],piv,p):
                F=v;break
        if F is None:
            assert not born
            continue
        born=True;b=F[1]*u+F[2]
        if last:
            assert b<=last[1]+(t-1)*(u-last[0])
            COUNTS['independent_threshold_Lipschitz_inequalities']+=1
        last=(u,b);samples.append((u,b))
        R,S=F[1:]
        if branch=='primal':
            aa=max(left,Q(S-Z,Y-R));bb=min(right,Q(S-Z+t,Y-R))
        else:
            aa=max(left,Q(S-Z+t,Y-R));bb=min(right,Q(S-Z,Y-R-t))
        if aa<bb:
            assert F[2]>=2*t
            active.append({'lo':aa,'hi':bb,'partner':F})
        COUNTS['whole_joint_prefix_frontier_cells']+=1
    partners={r['partner'] for r in active}
    bound=4*t*t if branch=='primal' else Q(16*a*t*t,Y)
    assert len(partners)<=bound
    if active:
        l=min(r['lo'] for r in active);h=max(r['hi'] for r in active);width=h-l
        geometric=Q(t,Y-t+1) if branch=='primal' else Q(t*a,(Y-2*t+1)**2)
        assert width<=geometric
        lower_events=len({x for x,v in entries if l<x<h})
        upper_events=len({Q(i,j) for i in range(1,A) for j in range(1,t) if l<Q(i,j)<h})
        assert lower_events<=t*(Y*width+1)
        assert upper_events<=Q(t*(t-1),2)*width+t-1
        COUNTS['entire_active_window_event_budgets']+=1
    return {'records':active,'partners':partners,'projection':result,'J':J}


def actual_certificate(data,cols,N,t,a,E,F,w,branch,p,front,H=None):
    lam,u,_=w;tau=dot(w,E)
    assert E[1]>=3*t and u>1 and tau==dot(w,F)
    assert (0<lam<t) if branch=='primal' else (t<lam<t*u)
    cheap=tuple(v for v in cols if dot(w,v)<tau)
    B,piv=field_basis(cheap,cols.__getitem__,p)
    assert field_reduce(cols[E],piv,p)[0] and field_reduce(cols[F],piv,p)[0]
    _,pF=field_basis((*B,F),cols.__getitem__,p)
    assert not field_reduce(cols[E],pF,p)[0]
    assert front and any(r['lo']<u<r['hi'] and r['partner']==F for r in front['records'])
    values,C,D,meta=front['projection']
    lower=tuple(v for v in C if dot(w,v)<tau)
    upper=tuple(v for v in D if dot(w,v)<tau)
    PB,pp=field_basis((*lower,*upper),values.__getitem__,p)
    assert field_reduce(values[E],pp,p)[0] and field_reduce(values[F],pp,p)[0]
    _,pf=field_basis((*PB,F),values.__getitem__,p)
    assert not field_reduce(values[E],pf,p)[0]
    if branch=='primal':
        assert H and not mul(N,H,a-t,t,p)
        local=(E[1]-t,E[2]-t)
        assert max(H,key=lambda v:v[0]*u+v[1])==local
        assert all(dot(w,v)<tau for v in meta['fixed'])
        COUNTS['actual_nonmonomial_primal_projection_certificates']+=1
    else:
        assert tau-lam>max((2*t-1)*u+t-1,(t-1)*u+2*t-1)
        assert all(not values[v] for v in cheap if
                   (v[0] and v[1]>=t) or (not v[0] and v[1]>=2*t))
        allupper=tuple(v for v in cols if v[0])
        _,pe=field_basis((*PB,E,*allupper),values.__getitem__,p)
        ell={data[0].index(v[1:]):field_reduce(values[v],pe,p)[1].get(E,0) for v in allupper}
        evaluate=lambda v:val(sum(ell.get(i,0)*c for i,c in cols[v].items()),p)
        assert all(not evaluate(v) for v in cheap)
        assert evaluate(E) and evaluate(F)
        KY=reflected(ell,data,a-t,t)[0]
        assert KY and not mul(N,KY,a-2*t,t,p)
        target=meta['target'];rho=lambda v:v[0]*u+v[1]
        assert max(map(rho,KY))==rho(target)
        COUNTS['actual_dual_restricted_space_certificates']+=1
    costs=sorted({dot(w,v) for v in cols})
    eps=min(lam/2,min(y-x for x,y in zip(costs,costs[1:]))/4)
    adj=[]
    for sign in (-1,1):
        wp=(lam+sign*eps,u,Q(1))
        S,_=field_basis(sorted(cols,key=lambda v:(dot(wp,v),v)),cols.__getitem__,p)
        assert all(dot(wp,q)>0 for q in normals(S,cols,p));adj.append(set(S))
        COUNTS['actual_adjacent_strict_minima']+=1
    assert F in adj[0] and E not in adj[0] and E in adj[1] and F not in adj[1]
    COUNTS['actual_coherent_partner_certificates']+=1
    return {'lower':E,'upper':F,'weight':list(map(str,w)),'branch':branch}


def run_case(t,a,units,p,E,branch,witnesses,H=None):
    data,co,cols=assemble(t,a,units,p)
    N={v:c for v,c in units[0].items() if v!=(0,0)}
    fronts=[frontier(data,cols,N,t,a,E,branch,J,p) for J in own_sectors(t)]
    partners=set().union(*(r['partners'] for r in fronts if r))
    bound=4*t**4 if branch=='primal' else Q(16*a*t**4,E[1])
    assert len(partners)<=bound
    cert=[]
    for F,w in witnesses:
        f=next(r for r in fronts if r and r['J'][0]<w[1] and (r['J'][1] is None or w[1]<r['J'][1]))
        cert.append(actual_certificate(data,cols,N,t,a,E,F,w,branch,p,f,H))
    return {'t':t,'a':a,'field_characteristic':p,'lower':E,'branch':branch,
            'coefficients':list(map(int,co)),'candidate_distinct_partners':len(partners),
            'candidate_frontiers':[{'sector':[str(v) if v is not None else None for v in r['J']],
                  'records':[{'lo':str(v['lo']),'hi':str(v['hi']),'partner':v['partner']} for v in r['records']]}
                  for r in fronts if r],'actual_certificates':cert}


def original_controls():
    records=[]
    for t,e,A,mu,p in [(3,1,12,1,0),(3,1,12,1,2),(3,1,12,2,101),(5,2,20,2,101)]:
        d=t-2*e;a=A+t;nu=A-t+mu
        U={(0,0):1,(mu,d+e):1,(nu,d):1};V={(0,0):1,(A-1-e,1):1}
        data,_,cols=assemble(t,a,(U,V),p)
        E=(0,A,t);F=(1,1,a-1);q=(1,-(A-1),A-1)
        w=wall_weight((Q(e),1+Q(e,A-1),Q(1)),(Q(t),1+Q(t,A-1),Q(1)),cols,q)
        bridge=(0,t,t+e);gy=(1,A,0);yl=(1,t,e);zl=(1,0,t+e)
        assert not combination(cols,{E:1,F:1,bridge:-1,gy:-1,yl:1,zl:1},p)
        H={(A-t,0):1,(0,e):val(-1,p)}
        records.append(run_case(t,a,(U,V),p,E,'primal',[(F,w)],H))
    # Three incomparable positive-valuation corners in both arms.
    t=4;a=16;A=a-t;p=101;units=[]
    for corners in [((1,3,1),(4,2,2),(8,1,1)),((2,3,1),(5,2,3),(10,1,1))]:
        N=Counter()
        for nu,d,c in corners:
            for i in range(A-nu):N[nu+i,d]+=c*comb(A-nu-1,i)
        units.append({(0,0):1,**clean(N,p)})
    E=(0,12,5);F=(1,1,14)
    ws=[(Q(1971,463),Q(558,463),Q(1)),(Q(1559,323),Q(406,323),Q(1))]
    records.append(run_case(t,a,units,p,E,'dual',[(F,w) for w in ws]))
    for m,p in [(4,101),(5,1009)]:
        A=2*m*m+m+2;a=A+2;t=2
        U={(0,0):1,(m,1):1};V={(0,0):1,**{(i,1):comb(A-1,i) for i in range(A)}}
        _,_,cols=assemble(t,a,(U,V),p)
        E=(0,m+2,3);ws=[]
        for j in range(2,m+1):
            K=m*m+j;J=A-1-K;u=Q(2*K-1,2*m);lam=(m+1)*u+1-J
            eps=min(Q(1,4*m),(lam-2)/(4*(m+1)),(2*u-lam)/(4*(m-1)))
            w=wall_weight((lam-(m+1)*eps,u-eps,Q(1)),(lam+(m+1)*eps,u+eps,Q(1)),cols,(1,-(m+1),J-1))
            ws.append(((1,1,J+2),w))
        records.append(run_case(t,a,(U,V),p,E,'dual',ws))
    return records


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    cases=original_controls()
    files={'scripts/run_hook_weighted_partner_checks.py','docs/hook_weighted_partners.md',
           'docs/hook_filtered_annihilators.md','docs/hook_exception_partners.md',
           'docs/hook_bridge_reduction.md','docs/cross_layer_direction_budget.md'}
    pending=['scripts/run_hook_weighted_partner_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=[node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else []
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:
                    files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_weighted_partner_checks.py',
            'scope':'General primal low-price and dual middle-price weighted distinct-partner bounds for arbitrary mixed central hooks in the stated orientation. No complete fan or unrestricted exponent bound.',
            'counts':dict(sorted(COUNTS.items())),'original_quotient_replay_counts':dict(REPLAY),'cases':cases,
            'limits':['Both outer images are zero.','The lower long axis has the larger coordinate weight.',
                      'Only the stated primal low-price and dual middle-price classes are covered.',
                      'Actual eligibility may recur; candidate thresholds count distinct partners, not connected pieces.',
                      'Exact priority and the original unrestricted endpoint remain open.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_weighted_partner_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}))


if __name__=='__main__':main()
