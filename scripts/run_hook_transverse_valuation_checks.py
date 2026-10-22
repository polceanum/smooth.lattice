#!/usr/bin/env python3
"""Original-field kernel and fan checks for nonunit transverse coefficients."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import comb
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_generic_hook_checks import bareiss
from run_hook_bridge_checks import actual,model,normals,polygon,COUNTS as REPLAY
from run_hook_transverse_prefix_checks import (
    clean,mul,degree,wall_weight,relation_control,same_layer_relation,
    COUNTS as WALL_COUNTS)
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import dot

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def fixture(t,a,ds,vs,p,rng,ks=None):
    data=model(t,a);A=a-t;units=[]
    ks=ks or tuple(A-v-1 for v in vs)
    for d,v,k in zip(ds,vs,ks):
        W={(0,0):1}
        for i in range(v,A):
            if i-v<=k:W[i,d]=comb(k,i-v)
        for i,j in product(range(A),range(d+1,t)):
            W[i,j]=rng.randrange(p or 101)
        units.append(clean(W,p))
    U,V=units
    f={(t+i,j):c for (i,j),c in U.items()}
    f.update({(j,t+i):c for (i,j),c in V.items()})
    co=[f.get((i,j),0) if g==1 else 0 for g,i,j in data[2]]
    cols=actual(data,co,p,replay=True)
    return data,co,cols,U,V,ks


def effective(F,m,p):
    F={i:c for i,c in F.items() if i<m and c}
    if not F:return {'m':m,'valuation':m,'unit_length':0,'unit_gap':None,'effective_defect':m-1,'lower_bound':0,'regular_indices':[]}
    v=min(F);s=m-v;G={i-v:c for i,c in F.items()}
    regular=[]
    for r in range(s+1):
        matrix=[[G.get(r+i-j,0) for j in range(s-r)] for i in range(s-r)]
        val=bareiss(matrix)
        rank=len(field_basis(range(s-r),lambda j:{i:matrix[i][j] for i in range(s-r) if matrix[i][j]},p)[0])
        assert bool(val%p if p else val)==(rank==s-r)
        if rank==s-r:regular.append(r)
        COUNTS['independent_residual_unit_prefix_ranks']+=1
    assert regular[0]==0 and regular[-1]==s
    h=max((j-i-1 for i,j in zip(regular,regular[1:])),default=0)
    e=max(h,2*v-m-1)
    assert 0<=e<=m-1
    return {'m':m,'valuation':v,'unit_length':s,'unit_gap':h,'effective_defect':e,'lower_bound':m-1-e,'regular_indices':regular}


def degree_controls():
    rng=Random(2026102281);records=[]
    cases=[
        (2,10,(1,1),(1,2),0,None),
        (3,14,(1,2),(4,4),101,None),
        (3,16,(2,1),(6,5),101,None),
        (4,18,(1,3),(8,7),101,None),
        (3,16,(1,2),(10,12),101,None),
        (3,15,(2,1),(2,3),101,(2,3)),
        (3,12,(1,2),(2,4),2,(3,2)),
        (4,14,(2,3),(4,6),3,(3,2)),
    ]
    for t,a,ds,vs,p,ks0 in cases:
        data,co,cols,U,V,ks=fixture(t,a,ds,vs,p,rng,ks0);A=a-t;arms=[]
        for arm,W,d,v,k in zip(('Y','Z'),(U,V),ds,vs,ks):
            F={i:c for (i,j),c in W.items() if j==d}
            assert not F.get(0) and min(F)==v
            COUNTS['genuinely_nonunit_first_transverse_polynomials']+=1
            for h in product(range(A),range(t)):
                hy,hz=h if arm=='Y' else h[::-1]
                original={}
                for row,c in cols[(0,t+hy,t+hz)].items():
                    i,j=data[0][row]
                    if arm=='Y' and i>=t:original[i-t,j]=c
                    elif arm=='Z' and j>=t:original[j-t,i]=c
                assert original==mul(W,{h:1},A,t,p)
                COUNTS['independent_original_nonunit_rectangular_columns']+=1
            for m in sorted({A,a-2*t}-{0}):
                meta=effective(F,m,p);ell=meta['lower_bound']
                rows=tuple(product(range(m),range(t)));idx={h:i for i,h in enumerate(rows)}
                values={h:{idx[z]:c for z,c in mul(W,{h:1},m,t,p).items()} for h in rows}
                for D in range(ell):
                    source=tuple(h for h in rows if h[0]<=D)
                    for E in range(ell-D):
                        rank=len(field_basis(source,lambda h:{i:c for i,c in values[h].items() if rows[i][0]>E},p)[0])
                        assert len(source)-rank==d*(min(D,E)+1)
                        COUNTS['all_vector_nonunit_conditional_kernel_dimensions']+=1
                for u in (Q(1,3),Q(1),Q(7,3)):
                    threshold=u*ell
                    for cap in sorted({degree({h:1},u) for h in rows if degree({h:1},u)<threshold}):
                        source=tuple(h for h in rows if degree({h:1},u)<=cap)
                        outcap=threshold-cap-Q(1,10007)
                        rank=len(field_basis(source,lambda h:{i:c for i,c in values[h].items() if degree({rows[i]:1},u)>outcap},p)[0])
                        expected=sum(h[1]>=t-d and degree({h:1},u)<=min(cap,outcap) for h in source)
                        assert len(source)-rank==expected
                        COUNTS['all_vector_weighted_nonunit_conditional_kernels']+=1
                if not p or p>2*A:
                    if v>=m:
                        predicted=0;H={(0,t-d-1):1};mechanism='truncated-first-coefficient'
                    else:
                        s=m-v;B=v+min(k,s-1);predicted=min(B,2*s)
                        if B<=2*s:H={(0,t-d-1):1};mechanism='nonzero-product'
                        else:H={(s,t-d-1):1};mechanism='multiplier-kernel'
                    assert ell==predicted
                    WH=mul(W,H,m,t,p)
                    assert min(j for i,j in H)<t-d
                    assert max(i for i,j in H)+max(i for i,j in WH)==ell
                    COUNTS['sharp_'+mechanism.replace('-','_')+'_witnesses']+=1
                    meta['sharp_witness']={'input':[[i,j,c] for (i,j),c in H.items()],'mechanism':mechanism}
                if meta['unit_length'] and not meta['unit_gap'] and 0<meta['valuation']<=(m+1)/2:
                    assert meta['effective_defect']==0
                    COUNTS['nonunit_half_length_no_loss_controls']+=1
                if meta['unit_length'] and 2*meta['valuation']>m+1:
                    assert meta['effective_defect']>=2*meta['valuation']-m-1
                    COUNTS['over_half_length_defect_controls']+=1
                meta['arm']=arm;meta['first_transverse_order']=d;arms.append(meta)
        records.append({'t':t,'a':a,'field_characteristic':p,'valuations':vs,'first_transverse_orders':ds,
                        'coefficients':list(map(int,co)),'modules':arms})
    assert COUNTS['sharp_multiplier_kernel_witnesses'] and COUNTS['sharp_truncated_first_coefficient_witnesses']
    return records


def fan_controls():
    rng=Random(2026102283);records=[]
    for t,a,ds,vs in [(3,30,(1,1),(4,5)),(4,32,(1,2),(8,9))]:
        p=101;data,co,cols,U,V,ks=fixture(t,a,ds,vs,p,rng)
        for W,d in zip((U,V),ds):
            F={i:c for (i,j),c in W.items() if j==d}
            assert min(F)>0
            assert all(effective(F,m,p)['effective_defect']==0 for m in (a-t,a-2*t))
        seen=set();done=set();cross=[];same=[]
        for u,fac in product([Q(4,5),Q(9,10),Q(19,20),Q(21,20),Q(11,10),Q(5,4)],[Q(1,3),Q(4,3),Q(3),Q(6)]):
            w=(fac*t*min(u,1)+Q(1,100003),u,Q(1))
            B,_=field_basis(sorted(cols,key=lambda h:(dot(w,h),h)),cols.__getitem__,p)
            if frozenset(B) in seen:continue
            seen.add(frozenset(B));NN=normals(B,cols,p);poly=polygon(NN)
            for left,right in zip(poly,poly[1:]+poly[:1]):
                if any(left[i]==right[i]==0 for i in range(3)):continue
                q=next(q for q in NN if dot(q,left)==dot(q,right)==0)
                if q[0]<0 or (not q[0] and q[1]>0):q=tuple(-x for x in q)
                if not(q[1]<0<q[2]) or q in done:continue
                if q[0] and len(cross)>=4:continue
                if not q[0] and len(same)>=4:continue
                ww=wall_weight(left,right,cols,q);lam,uu,_=ww
                found=False
                for e in cols:
                    if not(2*t<=e[1]<a and (t<=e[2]<2*t if not e[0] else e[2]<t)):continue
                    if q[0] and e[0]:continue
                    for scale in ((1,) if q[0] else range(1,a)):
                        f=tuple(x+scale*y for x,y in zip(e,q))
                        if f not in cols or not(2*t<=f[2]<a):continue
                        if q[0]:
                            if not(f[0]==1 and f[1]<t):continue
                            Q0=f[1]*uu+f[2]
                            if Q0<=max((2*t-1)*uu+t-1,(t-1)*uu+2*t-1):continue
                            witness=relation_control(data,cols,U,V,t,a,*ds,p,ww,e,f,(0,0))
                        else:
                            if not(f[0]==e[0] and (t<=f[1]<2*t if not f[0] else f[1]<t)):continue
                            witness=same_layer_relation(data,cols,U,V,t,a,*ds,p,ww,e,f)
                        if witness:
                            (cross if q[0] else same).append(witness);done.add(q);found=True;break
                    if found:break
            if len(cross)>=4 and len(same)>=4:break
        records.append({'t':t,'a':a,'field_characteristic':p,'valuations':vs,'first_transverse_orders':ds,
                        'coefficients':co,'sampled_original_cones':len(seen),'actual_cross_layer_walls':cross,
                        'actual_same_layer_walls':same,'scope':'Actual original-field controls with nonunit first coefficients and zero effective defect, not complete fan counts.'})
    assert sum(len(x['actual_cross_layer_walls']) for x in records)
    assert sum(len(x['actual_same_layer_walls']) for x in records)
    return records


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    fields=degree_controls();print('nonunit all-vector and sharpness controls passed',flush=True)
    fans=fan_controls();print('nonunit original coherent wall controls passed',flush=True)
    files={'scripts/run_hook_transverse_valuation_checks.py','docs/hook_transverse_valuation.md',
           'docs/hook_transverse_prefix.md','docs/hook_bridge_reduction.md','docs/cross_layer_direction_budget.md'}
    pending=['scripts/run_hook_transverse_valuation_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=[node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else []
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_transverse_valuation_checks.py',
            'scope':'All constant-axis central-only hooks admit the valuation-sensitive effective-defect bound. Nonunit valuations up to half the relevant lengths incur no loss when residual unit prefixes are regular. No unrestricted exponent improvement or polynomial lower separation.',
            'counts':dict(sorted(COUNTS.items())),'imported_actual_wall_counts':dict(WALL_COUNTS),'original_quotient_replay_counts':dict(REPLAY),
            'field_controls':fields,'original_fan_controls':fans,
            'limits':['Outer generator images are zero.','Sharpness concerns the conditional degree estimate, not the fan bound.',
                      'Original-field fan controls test the no-loss nonunit locus; positive-defect cases have full-subspace kernel controls and the uniform proof.',
                      'No complete fan count, unrestricted growth law or exact novelty is inferred.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_transverse_valuation_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}))


if __name__=='__main__':main()
