#!/usr/bin/env python3
"""Original-rank fiber projections with arbitrary transverse unit bases."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import comb
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_hook_bridge_checks import model,actual,normals,polygon,COUNTS as REPLAY
from run_hook_transverse_prefix_checks import mul,clean,degree,wall_weight
from run_hook_monomial_perturbation_checks import reflected,row_functional
from run_hook_transverse_valuation_checks import effective
from run_signed_face_checks import field_basis,field_reduce
from run_structural_complexity_checks import dot

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def inverse_base(C,t,p):
    inv={0:pow(int(C[0]),-1,p) if p else 1/Q(C[0])}
    for j in range(1,t):
        inv[j]=-inv[0]*sum(C.get(i,0)*inv[j-i] for i in range(1,j+1))
        if p:inv[j]%=p
    return {(0,j):c for j,c in inv.items() if c}


def assemble(t,a,units,p,replay=True):
    data=model(t,a);U,V=units
    f={(t+i,j):c for (i,j),c in U.items()}
    f.update({(j,t+i):c for (i,j),c in V.items()})
    co=[f.get((i,j),0) if g==1 else 0 for g,i,j in data[2]]
    return data,co,actual(data,co,p,replay=replay)


def fixture(t,a,ds,vs,p,rng):
    A=a-t;bases=[];units=[];Ns=[]
    for d,nu in zip(ds,vs):
        s=A-nu
        B={(i,0):comb(s-1,i) for i in range(s)}
        B.update({(i,j):rng.randrange(p or 101) for i,j in product(range(s),range(1,t-d))})
        N=clean({(i+nu,j+d):c for (i,j),c in B.items()},p)
        C={j:rng.randrange(p or 101) for j in range(1,t)}
        C[0]=1;C[t-1]=1
        U=Counter(N)
        for j,c in C.items():U[0,j]+=c
        bases.append(C);Ns.append(N);units.append(clean(U,p))
    data,co,cols=assemble(t,a,units,p)
    return data,co,cols,units,Ns,bases


def vector_controls():
    rng=Random(2026102301);records=[]
    for t,a,ds,vs,p in [(3,12,(1,2),(6,7),0),(4,16,(2,1),(8,10),101),
                        (5,20,(1,3),(10,13),101),(4,15,(1,2),(2,3),101),
                        (4,16,(2,1),(8,9),2),(5,18,(1,3),(8,11),3)]:
        data,co,cols,units,Ns,bases=fixture(t,a,ds,vs,p,rng);A=a-t;mods=[]
        for arm,U,N,C,d,nu in zip(('Y','Z'),units,Ns,bases,ds,vs):
            inv=inverse_base(C,t,p)
            normalized=mul(inv,U,A,t,p)
            assert normalized.get((0,0))==1
            F={i:c for (i,j),c in N.items() if j==d}
            metas=[effective(F,m,p) for m in sorted({A,a-2*t}-{0})]
            H=max((v['unit_gap'] or 0 for v in metas),default=0)
            loss=H+(t if 2*nu<A else 0)
            for h in product(range(A),range(t)):
                hy,hz=h if arm=='Y' else h[::-1]
                original={}
                for row,c in cols[(0,t+hy,t+hz)].items():
                    i,j=data[0][row]
                    if arm=='Y' and i>=t:original[i-t,j]=c
                    elif arm=='Z' and j>=t:original[j-t,i]=c
                assert original==mul(U,{h:1},A,t,p)
                COUNTS['independent_original_transverse_base_columns']+=1
            for meta in metas:
                m=meta['m'];s=max(0,m-nu);ell=max(0,m-1-loss)
                rows=tuple(product(range(m),range(t)));idx={v:i for i,v in enumerate(rows)}
                inside=lambda v:v[0]>=s or v[1]>=t-d
                vals={v:{idx[w]:c for w,c in mul(U,{v:1},m,t,p).items()} for v in rows}
                for D in range(ell):
                    source=tuple(v for v in rows if v[0]<=D)
                    for E in range(ell-D):
                        rank=len(field_basis(source,lambda v:{i:c for i,c in vals[v].items() if rows[i][0]>E},p)[0])
                        expected=sum(inside(v) and v[0]<=min(D,E) for v in source)
                        assert len(source)-rank==expected
                        # Degree-preserving row normalization is checked on the
                        # entire same source subspace, not just basis vectors.
                        normvals={v:{idx[w]:c for w,c in mul(normalized,{v:1},m,t,p).items()} for v in source}
                        nrank=len(field_basis(source,lambda v:{i:c for i,c in normvals[v].items() if rows[i][0]>E},p)[0])
                        assert rank==nrank
                        COUNTS['all_vector_long_degree_normalization_and_exception_kernels']+=1
                for u in (Q(1,3),Q(1),Q(7,3)):
                    threshold=u*ell
                    for cap in sorted({degree({v:1},u) for v in rows if degree({v:1},u)<threshold}):
                        source=tuple(v for v in rows if degree({v:1},u)<=cap)
                        outcap=threshold-cap-Q(1,10007)
                        tail=lambda v:{i:c for i,c in vals[v].items() if degree({rows[i]:1},u)>outcap}
                        rank=len(field_basis(source,tail,p)[0])
                        Jsource=tuple(v for v in source if inside(v))
                        jrank=len(field_basis(Jsource,tail,p)[0])
                        assert len(source)-rank==len(Jsource)-jrank
                        COUNTS['all_vector_weighted_fiber_exception_containments']+=1
                mods.append({'arm':arm,'m':m,'valuation':nu,'mixed_order':d,'unit_gap_bound':H,
                             'degree_loss_allowance':loss,'degree_sum_lower_bound':ell})
        records.append({'t':t,'a':a,'field_characteristic':p,'coefficients':list(map(int,co)),
                        'transverse_bases':[list(C.items()) for C in bases],'modules':mods})
    return records


def certificate(data,cols,units,Ns,t,a,ds,vs,p,w,e,f):
    lam,u,_=w;A=a-t;L=a-2*t;tau=dot(w,e);assert tau==dot(w,f)
    Q0=tau-lam
    if Q0<=max((2*t-1)*u+t-1,(t-1)*u+2*t-1):return None
    cheap=tuple(v for v in cols if dot(w,v)<tau)
    W,piv=field_basis(cheap,cols.__getitem__,p)
    re,_=field_reduce(cols[e],piv,p);rf,_=field_reduce(cols[f],piv,p)
    if not re or not rf:return None
    _,pf=field_basis((*W,f),cols.__getitem__,p)
    rem,rep=field_reduce(cols[e],pf,p)
    if rem:return None
    assert rep.get(f)
    H=Counter()
    if not e[0]:
        H[e[1]-t,e[2]-t]=1
        for v,c in rep.items():
            if not v[0]:H[v[1]-t,v[2]-t]-=c
    else:
        for v,c in rep.items():
            if not v[0]:H[v[1]-t,v[2]-t]+=c
    H=clean(H,p)
    hs=[{h:c for h,c in H.items() if h[0]<A and h[1]<t},
        {(j,i):c for (i,j),c in H.items() if i<t and j<A}]
    assert all(mul(N,h,A,t,p) for N,h in zip(Ns,hs))
    F=tau+Q0
    assert F>(a+t-1)*u+t and F>a+t*u+t-1
    upper=tuple(v for v in cols if v[0])
    full,piv=field_basis((*W,e,*upper),cols.__getitem__,p);assert len(full)==2*t*A
    ell={data[0].index(v[1:]):field_reduce(cols[v],piv,p)[1].get(e,0) for v in upper}
    evaluate=lambda ll,v:sum(ll.get(i,0)*c for i,c in cols[v].items())%p
    assert all(not evaluate(ell,v) for v in cheap)
    assert evaluate(ell,e)==1 and evaluate(ell,f)
    Ks=reflected(ell,data,A,t)
    for arm,(K,N,nu) in enumerate(zip(Ks,Ns,vs)):
        if mul(N,K,L,t,p):continue
        lo,hi=max(0,L-nu),A-nu
        if any(not(lo<=i<hi) for i,j in K):continue
        for q,j in product(range(L),range(t)):
            if lo<=q<hi or mul(N,{(q,j):1},A,t,p):continue
            null=row_functional({(q,j):1},arm,data,A,t)
            if any(evaluate(null,v) for v in cheap) or evaluate(null,e) or evaluate(null,f):continue
            for row,c in null.items():ell[row]=(ell.get(row,0)+c)%p
            assert all(not evaluate(ell,v) for v in cheap)
            assert evaluate(ell,e)==1 and evaluate(ell,f)
            COUNTS['validated_nullrow_functional_augmentations']+=1
            break
    Ks=reflected(ell,data,A,t);branches=[]
    for arm,(K,N,nu) in enumerate(zip(Ks,Ns,vs)):
        assert K and all(i<L for i,j in K)
        if mul(N,K,L,t,p):
            branches.append('nonscalar');continue
        lo,hi=max(0,L-nu),A-nu
        K1={v:c for v,c in K.items() if lo<=v[0]<hi}
        K0={v:c for v,c in K.items() if v not in K1}
        assert not mul(N,K0,A,t,p) and not mul(N,K1,L,t,p)
        ll0=row_functional(K0,arm,data,A,t);ll1=row_functional(K1,arm,data,A,t)
        own=e if not arm else f;other=f if not arm else e
        assert all(not evaluate(ll0,v) for v in cheap)
        assert not evaluate(ll0,own) and not evaluate(ll0,other)
        assert K1 and evaluate(ll1,own)==evaluate(ell,own)
        M=own[1] if not arm else own[2]
        assert nu+t<=M<min(a,nu+2*t)
        branches.append('whole-fiber-strip')
        COUNTS['original_whole_fiber_dual_strip_certificates']+=1
        if K0:COUNTS['nonzero_removable_whole_fiber_components']+=1
    if branches==['nonscalar','nonscalar']:
        assert F<(a+2*t-1)*u+3*t-2 and F<a+(3*t-2)*u+2*t-1
    values=sorted({dot(w,v) for v in cols});gap=min(y-x for x,y in zip(values,values[1:]))
    cross=e[0]!=f[0];eps=min((lam if cross else u)/2,gap/(4 if cross else 8*a))
    bases=[]
    for sign in (-1,1):
        wp=(lam+sign*eps,u,Q(1)) if cross else (lam,u+sign*eps,Q(1))
        B,_=field_basis(sorted(cols,key=lambda v:(dot(wp,v),v)),cols.__getitem__,p)
        assert all(dot(wp,q)>0 for q in normals(B,cols,p));bases.append(set(B))
    visible=(f in bases[0] and e not in bases[0] and e in bases[1] and f not in bases[1]) if cross else (e in bases[0] and f not in bases[0] and f in bases[1] and e not in bases[1])
    if not visible:return None
    COUNTS['actual_coherent_fiber_pair_certificates']+=1
    COUNTS['actual_adjacent_original_minima']+=2
    return {'weight':list(map(str,w)),'Y_label':e,'Z_label':f,'dual_branches':branches,'cross_layer':cross}


def fan_controls():
    rng=Random(2026102303);records=[]
    for t,a,ds,vs in [(3,30,(1,1),(23,26)),(4,32,(1,2),(21,25)),(3,30,(1,1),(18,20))]:
        p=101;data,co,cols,units,Ns,bases=fixture(t,a,ds,vs,p,rng)
        prefixes=[]
        for N,d in zip(Ns,ds):
            F={i:c for (i,j),c in N.items() if j==d}
            for m in (a-t,a-2*t):
                meta=effective(F,m,p);assert not meta['unit_length'] or not meta['unit_gap'];prefixes.append(meta)
        seen=set();done=set();witnesses=[]
        for u,fac in product([Q(4,5),Q(9,10),Q(19,20),Q(21,20),Q(11,10),Q(5,4),Q(vs[1],vs[0])],[Q(1,3),Q(4,3),Q(3),Q(6)]):
            w=(fac*t*min(u,1)+Q(1,100003),u,Q(1))
            B,_=field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)
            if frozenset(B) in seen:continue
            seen.add(frozenset(B));NN=normals(B,cols,p);poly=polygon(NN)
            for left,right in zip(poly,poly[1:]+poly[:1]):
                if any(left[i]==right[i]==0 for i in range(3)):continue
                q=next(q for q in NN if dot(q,left)==dot(q,right)==0)
                if q[0]<0 or (not q[0] and q[1]>0):q=tuple(-x for x in q)
                if not(q[1]<0<q[2]) or q in done:continue
                ww=wall_weight(left,right,cols,q);found=False
                for e in cols:
                    if not(2*t<=e[1]<a and (t<=e[2]<2*t if not e[0] else e[2]<t)):continue
                    if q[0] and e[0]:continue
                    for scale in ((1,) if q[0] else range(1,a)):
                        f=tuple(x+scale*y for x,y in zip(e,q))
                        if f not in cols or not(2*t<=f[2]<a):continue
                        if q[0] and not(f[0]==1 and f[1]<t):continue
                        if not q[0] and not(f[0]==e[0] and (t<=f[1]<2*t if not f[0] else f[1]<t)):continue
                        cert=certificate(data,cols,units,Ns,t,a,ds,vs,p,ww,e,f)
                        if cert:witnesses.append(cert);done.add(q);found=True;break
                    if found:break
                if len(witnesses)>=6:break
            if len(witnesses)>=6:break
        records.append({'t':t,'a':a,'field_characteristic':p,'valuations':vs,'coefficients':co,
                        'residual_prefix_controls':prefixes,'sampled_original_cones':len(seen),'actual_wall_certificates':witnesses})
    assert sum(len(r['actual_wall_certificates']) for r in records)
    assert COUNTS['original_whole_fiber_dual_strip_certificates']
    assert COUNTS['nonzero_removable_whole_fiber_components']
    return records


def counter_controls():
    records=[]
    for t,e,A,p in [(3,1,12,0),(3,1,12,2),(3,1,12,101),(5,2,20,2),(5,2,20,101),(7,3,32,101)]:
        d=t-2*e;nu=A-t;a=A+t
        U={(0,0):1,(0,d+e):1,(nu,d):1};V={(0,0):1,(A-1-e,1):1}
        data,co,cols=assemble(t,a,(U,V),p)
        Nt={v:c for v,c in U.items() if v!=(0,0)}
        H={(nu,0):1,(0,e):-1}
        assert not mul(Nt,H,A,t,p)
        assert all(mul(Nt,{v:1},A,t,p) for v in H)
        rows=tuple(product(range(A),range(t)))
        rank=len(field_basis(rows,lambda v:mul(Nt,{v:1},A,t,p),p)[0])
        assert rank==e*A
        E=(0,A,t);F=(1,1,a-1);q=(1,-(A-1),A-1)
        w=wall_weight((Q(e),1+Q(e,A-1),Q(1)),(Q(t),1+Q(t,A-1),Q(1)),cols,q)
        lam,u,_=w;assert e<lam<t
        tau=dot(w,E);assert tau==dot(w,F)
        bridge=(0,t,t+e);gy=(1,t+nu,0);gylow=(1,t,e);gzlow=(1,0,t+e)
        identity=Counter(cols[E]);identity.update(cols[F])
        for v,c in [(bridge,-1),(gy,-1),(gylow,1),(gzlow,1)]:
            for row,value in cols[v].items():identity[row]+=c*value
        assert not clean(identity,p)
        assert all(dot(w,v)<tau for v in (bridge,gy,gylow,gzlow))
        idx={v:i for i,v in enumerate(data[0])}
        ell={idx[t+nu,d+e]:1,idx[1,t+A-1]:-1}
        ev=lambda v:sum(ell.get(i,0)*c for i,c in cols[v].items())
        mod=lambda c:c%p if p else c
        assert all(not mod(ev(v)) for v in cols if dot(w,v)<tau)
        assert mod(ev(E)) and mod(ev(F))
        vals=sorted({dot(w,v) for v in cols});eps=min(lam/2,min(y-x for x,y in zip(vals,vals[1:]))/4)
        bs=[]
        for sign in (-1,1):
            wp=(lam+sign*eps,u,Q(1))
            B,_=field_basis(sorted(cols,key=lambda v:(dot(wp,v),v)),cols.__getitem__,p)
            assert all(dot(wp,r)>0 for r in normals(B,cols,p));bs.append(set(B))
            COUNTS['counterfamily_original_adjacent_minima']+=1
        assert F in bs[0] and E not in bs[0] and E in bs[1] and F not in bs[1]
        COUNTS['counterfamily_actual_nonmonomial_primal_annihilator_walls']+=1
        records.append({'t':t,'e':e,'d':d,'A':A,'field_characteristic':p,'weight':list(map(str,w)),
                        'lower':E,'upper':F,'direction':q,'annihilator_generator':'v^e-q^nu',
                        'rank_full_off_axis_multiplier':rank})
    return records


def transverse_only_controls():
    rng=Random(2026102307);records=[]
    for t,a,p in [(3,10,0),(3,10,2),(4,14,101)]:
        A=a-t
        U={(0,0):1}
        U.update({(0,j):rng.randrange(p or 101) for j in range(1,t)})
        U[0,t-1]=1
        V={(0,0):1}
        V.update({(i,j):rng.randrange(p or 101) for i,j in product(range(A),range(1,t))})
        data,co,cols=assemble(t,a,(clean(U,p),clean(V,p)),p)
        groups=[];removed=set();index={v:i for i,v in enumerate(data[0])}
        for q in range(t,A):
            rows={index[t+q,j] for j in range(t)}
            group=tuple(v for v,c in cols.items() if set(c)&rows)
            assert all(set(cols[v])<=rows and v[1]==t+q for v in group)
            assert len(field_basis(group,cols.__getitem__,p)[0])==t
            groups.append(group);removed.update(group)
            COUNTS['independent_transverse_only_direct_fibers']+=1
        rest={v:c for v,c in cols.items() if c and v not in removed}
        assert all(v[1]<2*t for v in rest)
        for lam,u in [(Q(t,3)+Q(1,101),Q(701,503)),(Q(3*t)+Q(1,103),Q(307,719))]:
            w=(lam,u,Q(1));parts=[]
            for group in groups:
                part,_=field_basis(sorted(group,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)
                parts.extend(part)
            B,_=field_basis(sorted(rest,key=lambda v:(dot(w,v),v)),rest.__getitem__,p)
            full,_=field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)
            assert set(parts)|set(B)==set(full)
            assert all(dot(w,q)>0 for q in normals(full,cols,p))
            COUNTS['independent_transverse_only_original_minimum_decompositions']+=1
        records.append({'t':t,'a':a,'field_characteristic':p,'direct_fibers':len(groups),'remaining_nonzero_labels':len(rest)})
    return records


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    fields=vector_controls();print('transverse-base all-vector controls passed',flush=True)
    fans=fan_controls();print('original whole-fiber wall controls passed',flush=True)
    counters=counter_controls();print('nonmonomial primal-annihilator family passed',flush=True)
    transverse=transverse_only_controls()
    files={'scripts/run_hook_transverse_base_checks.py','docs/hook_transverse_bases.md','docs/hook_monomial_perturbations.md',
           'docs/hook_transverse_valuation.md','docs/hook_transverse_prefix.md','docs/hook_bridge_reduction.md',
           'docs/cross_layer_direction_budget.md'};pending=['scripts/run_hook_transverse_base_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=[node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else []
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_transverse_base_checks.py',
            'scope':'Uniform all-valuation closure for transverse-unit bases plus monomial-times-unit mixed perturbations, and an actual infinite nonmonomial primal-annihilator wall family. No unrestricted exponent improvement or polynomial lower separation.',
            'counts':dict(sorted(COUNTS.items())),'original_quotient_replay_counts':dict(REPLAY),
            'all_vector_controls':fields,'original_fan_controls':fans,'nonmonomial_counterfamily_controls':counters,
            'transverse_only_controls':transverse,
            'limits':['Outer generator images are zero.','Mixed terms retain one monomial factor after removal of the transverse unit.',
                      'Row normalization preserves long degree; equality of transformed monomial fans is not asserted.',
                      'Directly verified null-row additions to some separators exercise nonzero removable whole-fiber components.',
                      'Counterfamily eligibility is proved on the stated interval, without a full fan count or claim outside it.',
                      'Exact priority and the original unrestricted endpoint remain open.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_transverse_base_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}))


if __name__=='__main__':main()
