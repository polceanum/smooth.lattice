#!/usr/bin/env python3
"""Exact rank-adaptive hook patches and cheaper-prefix coupling controls."""
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_hook_bridge_checks import model,actual,normals,polygon,strict_weight,COUNTS as REPLAY
from run_signed_face_checks import field_basis,field_reduce
from run_structural_complexity_checks import dot,area2

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def rank(labels,cols,p):
    return len(field_basis(tuple(labels),cols.__getitem__,p)[0])


def compress(data,cols,w,p):
    groups=[tuple(sorted(g,key=lambda v:dot(w,v))) for g in data[4]]
    assert len({dot(w,v) for v in cols})==len(cols)
    by,bz,bh=[field_basis(g,cols.__getitem__,p)[0] for g in groups]
    b=len(bh);U=(*by,*bz,*bh)
    ref,piv=field_basis((*by,*bz),cols.__getitem__,p)
    assert len(ref)==len(by)+len(bz)
    co=[field_reduce(cols[v],piv,p)[1] for v in bh]
    dual={v:{j:(-row.get(v,0))%p if p else -row.get(v,0)
             for j,row in enumerate(co) if row.get(v,0)} for v in ref}
    dual.update({v:{j:1} for j,v in enumerate(bh)})
    for j in range(b):
        summed=Counter()
        for v in U:
            for i,c in cols[v].items():summed[i]+=dual[v].get(j,0)*c
        assert all(not (c%p if p else c) for c in summed.values())
        COUNTS['independent_rank_adaptive_dual_rows']+=1
    small=[]
    for group in (by,bz,bh):
        keep,_=field_basis(tuple(reversed(group)),dual.__getitem__,p)
        assert len(keep)<=b;small.extend(keep)
    assert len(small)<=3*b and rank(small,dual,p)==b
    chosen,_=field_basis(sorted(cols,key=lambda v:dot(w,v)),cols.__getitem__,p)
    deleted,_=field_basis(sorted(small,key=lambda v:-dot(w,v)),dual.__getitem__,p)
    assert set(chosen)==set(U)-set(deleted)
    fullN=normals(chosen,cols,p)
    tinyN=normals(deleted,{v:dual[v] for v in small},p,True)
    orders={tuple(y-x for x,y in zip(v,u)) for group in groups
            for v,u in zip(group,group[1:])}
    P=polygon(fullN|orders);T=polygon(tinyN|orders)
    assert all(dot(v,q)>=0 for v in P for q in tinyN)
    assert all(dot(v,q)>=0 for v in T for q in fullN)
    assert area2(P)==area2(T)
    COUNTS['independent_entire_rank_adaptive_patch_equalities']+=1
    return {'bridge_rank':b,'possible_bridge_columns':len(groups[2]),
            'small_dual_labels':len(small),'patch_area2':str(area2(P))}


def arm_data(data,cols):
    Y={i for i,(y,z) in enumerate(data[0]) if y>=data[1][1][0]}
    Z={i for i,(y,z) in enumerate(data[0]) if z>=data[1][1][1]}
    hy={v:{i:c for i,c in cols[v].items() if i in Y} for v in data[4][2]}
    hz={v:{i:c for i,c in cols[v].items() if i in Z} for v in data[4][2]}
    return Y,Z,hy,hz


def interface(S,cols,hy,hz,p):
    k=rank(S,hy,p)+rank(S,hz,p)-rank(S,cols,p)
    assert k>=0
    return k


def unknown_walls(data,cols,p,source):
    cells={}
    for raw in source['support_fan_bases']:
        B=frozenset(tuple(v) for v in raw if tuple(v) in cols)
        N=normals(B,cols,p);poly=polygon(N)
        assert len(B)==rank(cols,cols,p)
        cells[B]=(poly,N)
    assert len(cells)==89 and sum(area2(v[0]) for v in cells.values())==1
    Y,Z,hy,hz=arm_data(data,cols);found={};patches=[]
    for B,(poly,NN) in cells.items():
        patches.append(compress(data,cols,strict_weight(poly,cols),p))
        for left,right in zip(poly,poly[1:]+poly[:1]):
            if any(left[i]==right[i]==0 for i in range(3)):continue
            q=next(q for q in NN if dot(q,left)==dot(q,right)==0)
            if not q[0]:continue
            if q[0]<0:q=tuple(-v for v in q)
            if q in data[7]:continue
            edge=tuple(sorted((left,right)))
            if edge in found:continue
            w=tuple((x+y)/2 for x,y in zip(left,right))
            costs=sorted({dot(w,v) for v in cols})
            eps=min(w[0]/2,min(y-x for x,y in zip(costs,costs[1:]))/4)
            bases=[]
            for sign in (-1,1):
                wp=(w[0]+sign*eps,w[1],w[2])
                S,_=field_basis(sorted(cols,key=lambda v:(dot(wp,v),v)),cols.__getitem__,p)
                assert all(dot(wp,n)>0 for n in normals(S,cols,p));bases.append(set(S))
            e,f=tuple(bases[0]-bases[1]),tuple(bases[1]-bases[0])
            assert len(e)==len(f)==1;e,f=e[0],f[0]
            assert (e in data[4][0] and f in data[4][1]) or (e in data[4][1] and f in data[4][0])
            tau=dot(w,e);assert tau==dot(w,f)
            cheap=tuple(v for v in cols if dot(w,v)<tau)
            cb,cp=field_basis(cheap,cols.__getitem__,p)
            assert field_reduce(cols[e],cp,p)[0] and field_reduce(cols[f],cp,p)[0]
            _,cf=field_basis((*cb,f),cols.__getitem__,p)
            assert not field_reduce(cols[e],cf,p)[0]
            H=tuple(v for v in data[4][2] if dot(w,v)<tau)
            k=interface(H,cols,hy,hz,p);assert k>0
            iy=tuple(v for v in cols if v[0] and v[1]>=2)
            iz=tuple(v for v in cols if v[0] and v[2]>=2)
            d=rank((*cheap,*iy),cols,p)+rank((*cheap,*iz),cols,p)-rank(cheap,cols,p)-rank(cols,cols,p)
            assert 1<=d<=k
            found[edge]={'direction':q,'exchanged_labels':[e,f],'weight':list(map(str,w)),
                         'cheaper_bridge_labels':H,'prefix_coupling_rank':k,'remaining_arm_intersection_rank':d}
            COUNTS['independent_original_unknown_wall_certificates']+=1
    assert {tuple(v['direction']) for v in found.values()}=={(1,-7,2),(1,2,-7)}
    return list(found.values()),patches


def power_control(base,coeff,s,p):
    data=model(2*s,5*s)
    lookup={(g,s*i,s*j):c for (g,i,j),c in zip(base[2],coeff)}
    co=[lookup.get(v,0) for v in data[2]]
    cols=actual(data,co,p,True);Y,Z,hy,hz=arm_data(data,cols);H=data[4][2]
    rb=rank(H,cols,p);ry=rank(H,hy,p);rz=rank(H,hz,p)
    assert (rb,ry,rz)==(12*s*s,6*s*s,6*s*s) and interface(H,cols,hy,hz,p)==0
    for v in H:
        tag=(v[1]%s,v[2]%s);old=(0,v[1]//s,v[2]//s)
        actual_block={(y//s,z//s):c for i,c in cols[v].items() for y,z in [data[0][i]]}
        assert all((data[0][i][0]%s,data[0][i][1]%s)==tag for i in cols[v])
        original=actual(base,coeff,p)
        assert actual_block=={base[0][i]:c for i,c in original[old].items()}
        COUNTS['independent_full_bridge_residue_block_equalities']+=1
    cut=Q(10*s-1,2)
    prefix=tuple(v for v in H if v[1]+v[2]<cut)
    expected={(0,2*s+i,2*s+j) for i in range(s) for j in range(s-i)}
    assert set(prefix)==expected
    k=interface(prefix,cols,hy,hz,p)
    assert k==s*(s+1)//2 and k>=Q(32*s*s,64)
    for q in ((1,-7*s,2*s),(1,2*s,-7*s)):
        assert q not in data[7]
    COUNTS['exact_infinite_family_prefix_specializations']+=1
    return {'scale':s,'length':32*s*s,'bridge_rank':rb,'projected_bridge_ranks':[ry,rz],
            'full_span_coupling_rank':0,'prefix_cost_upper_bound':str(cut),'prefix_coupling_rank':k,
            'coefficients':co,'fan_transport':'Direct residue blocks at transformed weights (w_x,s*w_y,s*w_z); fan count 89 follows analytically.'}


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    saved=json.loads((ROOT/'results/research/hook_bridge_checks.json').read_text())
    for f,sha in saved['source_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==sha
    source=json.loads((ROOT/'results/research/generic_hook_checks.json').read_text())['families'][0]
    base=model(2,5);coeff=saved['complete_generic_fixture']['common_coefficients'];p=1000003
    cols=actual(base,coeff,p,True)
    walls,patches=unknown_walls(base,cols,p,source)
    print('All 89 original cone patches and both unknown walls passed.',flush=True)
    rng=Random(20261023);controls=[]
    for t,a,char,kind in [(2,5,0,'zero'),(3,10,2,'zero'),(2,6,0,'one'),
                          (3,10,2,'one'),(4,17,3,'one'),(2,6,2,'random'),
                          (3,10,3,'random'),(3,10,101,'random')]:
        data=model(t,a)
        if kind=='zero':co=[0]*len(data[2])
        elif kind=='one':co=[int(g==1 and ((i,j)==(a-1,t-1) or (i,j)==(t-1,a-1))) for g,i,j in data[2]]
        else:co=[rng.randrange(char) for v in data[2]]
        values=actual(data,co,char,True)
        for trial in range(3):
            w=tuple(rng.randrange(1,1000000) for i in range(3))
            while len({dot(w,v) for v in values})<len(values):w=tuple(rng.randrange(1,1000000) for i in range(3))
            row=compress(data,values,w,char)
            if kind in ('zero','one'):assert row['bridge_rank']==int(kind=='one')
            controls.append({'t':t,'a':a,'field_characteristic':char,'kind':kind,**row})
    powers=[power_control(base,coeff,s,p) for s in (1,2,3)]
    files={'docs/hook_prefix_bridge_rank.md','scripts/run_hook_prefix_bridge_rank_checks.py',
           'docs/hook_bridge_reduction.md','docs/hook_barrier_sectors.md','docs/two_colour_threshold_braids.md',
           'results/research/hook_bridge_checks.json','results/research/generic_hook_checks.json'}
    pending=['scripts/run_hook_prefix_bridge_rank_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=[node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else []
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:
                    files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_prefix_bridge_rank_checks.py',
            'scope':'Rank-adaptive exact dual compression and an infinite actual obstruction to replacing cheaper-prefix coupling by full-span coupling. No worst-case thickness or unrestricted exponent improvement, polynomial separation or current novelty claim.',
            'counts':dict(sorted(COUNTS.items())),'original_quotient_replay_counts':dict(REPLAY),
            'rank_adaptive_patch_controls':controls,'complete_generic_original_patches':patches,
            'original_unknown_walls':walls,'free_power_controls':powers,
            'uniform_bounds':{'full_hook':'O(n^2*(1+b)^3), b=rank(H)',
                              'infinite_full_span_coupling':0,'infinite_prefix_coupling':'s*(s+1)/2 >= n_s/64, n_s=32*s^2',
                              'infinite_fan_count':89,'fan_count_transport_proved_analytically':True},
            'limits':['Finite controls do not prove the uniform statements or novelty.',
                      'The full bridge rank bound is a corollary of saved greedy and two-colour inputs.',
                      'Large prefix coupling in the power family does not produce fan growth.',
                      'No aggregate across moving prefixes or unrestricted endpoint is established.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_prefix_bridge_rank_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}))


if __name__=='__main__':main()
