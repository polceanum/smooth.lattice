#!/usr/bin/env python3
"""Original-rank annihilator splitting and coherent translated-strip checks."""
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
from run_hook_transverse_valuation_checks import effective
from run_signed_face_checks import field_basis,field_reduce
from run_structural_complexity_checks import dot

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def fixture(t,a,ds,vs,p,rng,ks=None):
    data=model(t,a);A=a-t;units=[];perturbations=[];Bs=[]
    ks=ks or tuple(A-v-1 for v in vs)
    for d,v,k in zip(ds,vs,ks):
        s=A-v
        B={(i,0):comb(k,i) for i in range(min(k,s-1)+1)}
        B.update({(i,j):rng.randrange(p or 101) for i,j in product(range(s),range(1,t-d))})
        B=clean(B,p)
        N={(i+v,j+d):c for (i,j),c in B.items()}
        U={(0,0):1,**N}
        Bs.append(B);perturbations.append(N);units.append(U)
    U,V=units
    f={(t+i,j):c for (i,j),c in U.items()}
    f.update({(j,t+i):c for (i,j),c in V.items()})
    co=[f.get((i,j),0) if g==1 else 0 for g,i,j in data[2]]
    cols=actual(data,co,p,replay=True)
    return data,co,cols,units,perturbations,Bs


def all_vector_controls():
    rng=Random(2026102293);records=[]
    cases=[(3,14,(1,2),(7,8),0,None),(3,18,(1,1),(9,12),101,None),
           (4,20,(2,1),(10,13),101,None),(5,22,(1,3),(12,15),101,None),
           (4,16,(1,2),(8,10),2,(2,1)),(4,18,(2,1),(9,11),3,(2,1))]
    for t,a,ds,vs,p,ks in cases:
        data,co,cols,units,Ns,Bs=fixture(t,a,ds,vs,p,rng,ks);A=a-t;modules=[]
        for arm,U,N,B,d,nu in zip(('Y','Z'),units,Ns,Bs,ds,vs):
            assert 2*nu>=A
            F={i:c for (i,j),c in N.items() if j==d}
            for h in product(range(A),range(t)):
                hy,hz=h if arm=='Y' else h[::-1]
                original={}
                for row,c in cols[(0,t+hy,t+hz)].items():
                    i,j=data[0][row]
                    if arm=='Y' and i>=t:original[i-t,j]=c
                    elif arm=='Z' and j>=t:original[j-t,i]=c
                assert original==mul(U,{h:1},A,t,p)
                COUNTS['independent_original_factorized_unit_columns']+=1
            for m in sorted({A,a-2*t}-{0}):
                meta=effective(F,m,p);s=max(0,m-nu);h=meta['unit_gap'] or 0
                ell=m-1-h;rows=tuple(product(range(m),range(t)));idx={v:i for i,v in enumerate(rows)}
                inside=lambda v:v[0]>=s or v[1]>=t-d
                nv={v:{idx[w]:c for w,c in mul(N,{v:1},m,t,p).items()} for v in rows}
                uv={v:{idx[w]:c for w,c in mul(U,{v:1},m,t,p).items()} for v in rows}
                rank=len(field_basis(rows,nv.__getitem__,p)[0])
                assert rank==s*(t-d)
                assert all(not nv[v] for v in rows if inside(v))
                COUNTS['independent_full_monomial_annihilator_rank_equalities']+=1
                for D in range(ell):
                    source=tuple(v for v in rows if v[0]<=D)
                    for E in range(ell-D):
                        rank=len(field_basis(source,lambda v:{i:c for i,c in uv[v].items() if rows[i][0]>E},p)[0])
                        expected=sum(inside(v) and v[0]<=min(D,E) for v in source)
                        assert len(source)-rank==expected
                        COUNTS['all_vector_annihilator_exception_kernel_dimensions']+=1
                for u in (Q(1,3),Q(1),Q(7,3)):
                    threshold=u*ell
                    for cap in sorted({degree({v:1},u) for v in rows if degree({v:1},u)<threshold}):
                        source=tuple(v for v in rows if degree({v:1},u)<=cap)
                        outcap=threshold-cap-Q(1,10007)
                        rank=len(field_basis(source,lambda v:{i:c for i,c in uv[v].items() if degree({rows[i]:1},u)>outcap},p)[0])
                        expected=sum(inside(v) and degree({v:1},u)<=min(cap,outcap) for v in source)
                        assert len(source)-rank==expected
                        COUNTS['all_vector_weighted_annihilator_exception_kernels']+=1
                if 0<s<m and 2*s<m-1-h:
                    witness={(s,t-d-1):1}
                    assert not mul(N,witness,m,t,p)
                    assert mul(U,witness,m,t,p)==witness
                    COUNTS['previous_sharp_degree_witnesses_in_exact_annihilator']+=1
                modules.append({'arm':arm,'m':m,'valuation':nu,'annihilator_long_threshold':s,
                                'transverse_order':d,'unit_prefix_gap':h,'nonscalar_degree_sum_lower_bound':ell,
                                'old_effective_defect':meta['effective_defect']})
        records.append({'t':t,'a':a,'field_characteristic':p,'coefficients':list(map(int,co)),'modules':modules})
    return records


def reflected(ell,data,A,t):
    out=[{},{}]
    for row,c in ell.items():
        if not c:continue
        i,j=data[0][row]
        if i>=t:out[0][A-1-(i-t),t-1-j]=c
        else:out[1][A-1-(j-t),t-1-i]=c
    return out


def row_functional(K,arm,data,A,t):
    result={}
    index={v:i for i,v in enumerate(data[0])}
    for (q,v),c in K.items():
        long,short=A-1-q,t-1-v
        row=(t+long,short) if not arm else (short,t+long)
        result[index[row]]=c
    return result


def pair_certificate(data,cols,units,Ns,t,a,ds,vs,p,w,e,f):
    lam,u,_=w;A=a-t;L=a-2*t;tau=dot(w,e)
    assert tau==dot(w,f);Q0=tau-lam
    if Q0<=max((2*t-1)*u+t-1,(t-1)*u+2*t-1):return None
    cheap=tuple(v for v in cols if dot(w,v)<tau)
    W,piv=field_basis(cheap,cols.__getitem__,p)
    re,_=field_reduce(cols[e],piv,p);rf,_=field_reduce(cols[f],piv,p)
    if not re or not rf:return None
    _,pf=field_basis((*W,f),cols.__getitem__,p)
    rem,rep=field_reduce(cols[e],pf,p)
    if rem:return None
    assert rep.get(f);H=Counter()
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
    D=tau-t*(u+1);F=tau+Q0
    for arm,h in enumerate(hs):
        weight=u if not arm else 1/u;scale=1 if not arm else u
        source=degree(h,weight)*scale
        image=degree(mul(units[arm],h,A,t,p),weight)*scale
        original=e if not arm else f
        if not original[0]:assert source==D and image<Q0-(t*u if not arm else t)
        else:assert source<D and image==Q0-(t*u if not arm else t)
        assert mul(Ns[arm],h,A,t,p),'Primal annihilator exception must be ineligible'
    assert F>(a+t-1)*u+t and F>a+t*u+t-1
    upper=tuple(v for v in cols if v[0])
    full,piv=field_basis((*W,e,*upper),cols.__getitem__,p);assert len(full)==2*t*A
    ell={data[0].index(v[1:]):field_reduce(cols[v],piv,p)[1].get(e,0) for v in upper}
    evaluate=lambda ell,v:sum(ell.get(i,0)*c for i,c in cols[v].items())%p
    assert all(not evaluate(ell,v) for v in cheap)
    assert evaluate(ell,e)==1 and evaluate(ell,f)
    Ks=reflected(ell,data,A,t)
    # Exercise a genuinely nonzero removable component when an independently
    # validated null row is available; this does not change either pair value.
    for arm,(K,N,nu,d) in enumerate(zip(Ks,Ns,vs,ds)):
        if mul(N,K,L,t,p) or any(i>=A-nu or j>=t-d for i,j in K):continue
        found=False
        for q,j in product(range(L),range(t)):
            if q<A-nu and j<t-d:continue
            null=row_functional({(q,j):1},arm,data,A,t)
            if any(evaluate(null,v) for v in cheap) or evaluate(null,e) or evaluate(null,f):continue
            for row,c in null.items():ell[row]=(ell.get(row,0)+c)%p
            COUNTS['validated_separator_nullrow_augmentations']+=1
            found=True;break
        if found:
            assert all(not evaluate(ell,v) for v in cheap)
            assert evaluate(ell,e)==1 and evaluate(ell,f)
    Ks=reflected(ell,data,A,t);branches=[]
    for arm,(K,N,nu,d) in enumerate(zip(Ks,Ns,vs,ds)):
        assert K and all(i<L for i,j in K)
        own=e if not arm else f
        other=f if not arm else e
        if mul(N,K,L,t,p):
            branches.append('nonscalar')
            continue
        K0={v:c for v,c in K.items() if v[0]>=A-nu or v[1]>=t-d}
        K1={v:c for v,c in K.items() if v not in K0}
        assert not mul(N,K0,A,t,p)
        assert not mul(N,K1,L,t,p)
        assert K1 and all(max(0,L-nu)<=i<A-nu and j<t-d for i,j in K1)
        ell0=row_functional(K0,arm,data,A,t)
        ell1=row_functional(K1,arm,data,A,t)
        assert all(not evaluate(ell0,v) for v in cheap)
        assert not evaluate(ell0,other) and not evaluate(ell0,own)
        assert evaluate(ell1,own)==evaluate(ell,own) and evaluate(ell1,own)
        M=own[1] if not arm else own[2]
        assert nu+t<=M<min(a,nu+2*t)
        branches.append('annihilator-strip')
        COUNTS['validated_full_annihilator_functional_splittings']+=1
        if K0:COUNTS['nonzero_full_annihilator_components']+=1
        if nu<L:COUNTS['nontrivial_truncated_annihilator_exceptions']+=1
        else:COUNTS['zero_perturbation_truncation_exceptions']+=1
    if branches==['nonscalar','nonscalar']:
        assert F<(a+2*t-1)*u+3*t-2 and F<a+(3*t-2)*u+2*t-1
    # Strict adjacent minima independently establish a visible original wall.
    values=sorted({dot(w,v) for v in cols})
    gap=min(y-x for x,y in zip(values,values[1:]))
    cross=e[0]!=f[0]
    eps=min((lam if cross else u)/2,gap/(4 if cross else 8*a))
    bases=[]
    for sign in (-1,1):
        wp=(lam+sign*eps,u,Q(1)) if cross else (lam,u+sign*eps,Q(1))
        B,_=field_basis(sorted(cols,key=lambda v:(dot(wp,v),v)),cols.__getitem__,p)
        assert all(dot(wp,q)>0 for q in normals(B,cols,p))
        bases.append(set(B))
    if cross:
        visible=f in bases[0] and e not in bases[0] and e in bases[1] and f not in bases[1]
    else:
        visible=e in bases[0] and f not in bases[0] and f in bases[1] and e not in bases[1]
    if not visible:return None
    COUNTS['actual_original_primal_annihilator_exclusions']+=2
    COUNTS['actual_original_adjacent_strict_minima']+=2
    COUNTS['actual_coherent_pair_certificates']+=1
    return {'weight':list(map(str,w)),'Y_label':e,'Z_label':f,'cross_layer':cross,
            'dual_branches':branches,'raw_direction':tuple(j-i for i,j in zip(e,f)),
            'primal_source_terms':len(H)}


def fan_controls():
    rng=Random(2026102297);records=[]
    cases=[(3,30,(1,1),(23,26)),(4,32,(1,2),(21,25)),(3,30,(1,1),(18,20))]
    for t,a,ds,vs in cases:
        p=101;data,co,cols,units,Ns,Bs=fixture(t,a,ds,vs,p,rng)
        assert all(2*v>=a-t for v in vs)
        prefix_controls=[]
        for arm,N,d,nu in zip(('Y','Z'),Ns,ds,vs):
            F={i:c for (i,j),c in N.items() if j==d}
            for m in (a-t,a-2*t):
                meta=effective(F,m,p)
                assert not meta['unit_length'] or meta['unit_gap']==0
                if meta['unit_length']:COUNTS['independent_actual_regular_residual_unit_prefixes']+=1
                if meta['effective_defect']>t:COUNTS['large_previous_valuation_defect_controls']+=1
                prefix_controls.append({'arm':arm,**meta})
        seen=set();done=set();witnesses=[]
        us=sorted({Q(4,5),Q(9,10),Q(19,20),Q(21,20),Q(11,10),Q(5,4),Q(vs[1],vs[0]),Q(vs[1]+1,vs[0]+1)})
        for u,fac in product(us,[Q(1,3),Q(4,3),Q(3),Q(6)]):
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
                        witness=pair_certificate(data,cols,units,Ns,t,a,ds,vs,p,ww,e,f)
                        if witness:
                            witnesses.append(witness);done.add(q);found=True;break
                    if found:break
                if len(witnesses)>=8:break
            if len(witnesses)>=8:break
        records.append({'t':t,'a':a,'field_characteristic':p,'first_transverse_orders':ds,'valuations':vs,
                        'coefficients':co,'sampled_original_cones':len(seen),'actual_wall_certificates':witnesses,
                        'residual_unit_prefix_controls':prefix_controls,
                        'scope':'Actual original-field primal exclusions and dual strip certificates, not complete fan counts.'})
    assert sum(len(r['actual_wall_certificates']) for r in records)
    assert COUNTS['validated_full_annihilator_functional_splittings']
    assert COUNTS['nontrivial_truncated_annihilator_exceptions']
    assert COUNTS['nonzero_full_annihilator_components']
    return records


def low_valuation_controls():
    for A in range(3,161):
        for t in range(2,min(A,10)+1):
            for nu in range((A-1)//2+1):
                for m in {A,A-t}-{0}:
                    if nu>=m:
                        assert m-1<t
                    else:
                        s=m-nu
                        for h in {0,s//2,s-1}:
                            e=max(h,2*nu-m-1)
                            assert e<=h+t
                    COUNTS['independent_low_valuation_defect_enclosures']+=1


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    fields=all_vector_controls();print('annihilator and nonscalar uncertainty controls passed',flush=True)
    fans=fan_controls();print('original dual strip and primal exclusion controls passed',flush=True)
    low_valuation_controls()
    files={'scripts/run_hook_monomial_perturbation_checks.py','docs/hook_monomial_perturbations.md',
           'docs/hook_transverse_valuation.md','docs/hook_transverse_prefix.md','docs/hook_bridge_reduction.md',
           'docs/cross_layer_direction_budget.md'};pending=['scripts/run_hook_monomial_perturbation_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=[node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else []
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_monomial_perturbation_checks.py',
            'scope':'Near-quadratic closure for constant-axis central hooks with monomial-times-unit arm perturbations and regular residual unit prefixes, uniformly in all valuations. No unrestricted improvement or lower separation.',
            'counts':dict(sorted(COUNTS.items())),'original_quotient_replay_counts':dict(REPLAY),
            'all_vector_controls':fields,'original_fan_controls':fans,
            'limits':['Both outer generator images are zero.','Higher transverse terms must retain the same long-variable monomial factor.',
                      'The common complete weight list fixes the two valuation parameters.','Finite fan controls exercise cheap-overlap original walls; no complete fan count is inferred.',
                      'Some separating functionals are augmented by directly verified null rows to exercise nonzero removable components.',
                      'Exact priority and the original unrestricted target remain unresolved.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_monomial_perturbation_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}))


if __name__=='__main__':main()
