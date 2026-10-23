#!/usr/bin/env python3
"""Exact passive-layer elimination and complete original hook-fan patches."""
from collections import Counter
from fractions import Fraction as Q
from itertools import permutations
from math import gcd,comb
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_hook_filtered_annihilator_checks import assemble
from run_hook_bridge_checks import normals,polygon,COUNTS as REPLAY
from run_hook_prefix_bridge_rank_checks import compress as small_decision
from run_hook_transverse_prefix_checks import mul
from run_signed_face_checks import field_basis,field_reduce
from run_structural_complexity_checks import dot,area2,clip
from run_oracle_enumeration_checks import exact_adjacent_weight

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()
TRIANGLE=[tuple(Q(i==j) for i in range(3)) for j in range(3)]
PHASES=tuple(permutations(range(3)))


def val(c,p):
    return c%p if p else Q(c)


def phase_normals(t,order):
    scales=(t,1,1);out=set()
    for i,j in zip(order,order[1:]):
        q=[0,0,0];q[i]=-scales[j];q[j]=scales[i];out.add(tuple(q))
    return out


def phase_key(t,w):
    v=(w[0]/Q(t),w[1],w[2]);assert len(set(v))==3
    return tuple(sorted(range(3),key=lambda i:v[i]))


def choose_weight(poly,labels,start=1):
    labels=set(labels)
    for seed in range(start,start+5000):
        w=tuple(sum(Q(seed+1)**i*v[j] for i,v in enumerate(poly)) for j in range(3))
        if min(w)>0 and len({dot(w,v) for v in labels})==len(labels):return w
    raise AssertionError('No strict weight found in the full-dimensional patch.')


def witness_labels(t,cols):
    return set(cols)|{(1,0,0),(0,t,0),(0,0,t),(0,1,0),(0,0,1)}


def passive_reduce(data,cols,t,a,s,w,p):
    A=a-t;B=set(range(s))|set(range(t-s,t));P=set(range(t))-B;m=len(B)
    index={v:i for i,v in enumerate(data[0])}
    values={v:dict(c) for v,c in cols.items() if c};fixed=[];removed=set()

    def component(labels,rows,r):
        labels=tuple(labels);rows=set(rows)
        assert all(v not in removed for v in labels)
        assert {v for v,c in values.items() if rows&set(c)}==set(labels)
        assert all(set(values[v])<=rows for v in labels)
        keep=tuple(sorted(labels,key=lambda v:dot(w,v))[:r])
        qb,qp=field_basis(keep,values.__getitem__,p);assert len(qb)==r==len(rows)
        for v in labels:
            values[v]=field_reduce(values[v],qp,p)[0]
            assert not values[v]
        fixed.extend(keep);removed.update(labels)
        COUNTS['independent_isolated_scalar_components']+=1

    def leaf(h,u,row,target):
        assert h not in removed and u not in removed
        assert {v for v,c in values.items() if c.get(row)}=={h,u}
        early=min((h,u),key=lambda v:dot(w,v));late=u if early==h else h
        original=dict(values[early])
        inv=pow(int(original[row]),-1,p) if p else 1/Q(original[row])
        # Eliminate at the passive row; a generic leading-row pivot would
        # mix the arm coordinates and conceal the geometric invariant.
        for v in (h,u):
            c=val(values[v].get(row,0)*inv,p)
            vv=Counter(values[v]);vv.subtract({r:c*x for r,x in original.items()})
            values[v]={r:val(x,p) for r,x in vv.items() if val(x,p)}
        assert not values[early] and values[late]
        assert set(values[late])<=target and dot(w,early)<dot(w,late)
        fixed.append(early);removed.add(early)
        COUNTS['independent_cost_compatible_leaf_eliminations']+=1

    RY={r for r,(y,z) in enumerate(data[0]) if y>=t and z in B}
    RZ={r for r,(y,z) in enumerate(data[0]) if z>=t and y in B}
    for q in range(t,A):
        for v in P:
            component(((1,t+q,v),(0,t+q,t+v)),(index[t+q,v],),1)
            component(((1,v,t+q),(0,t+v,t+q)),(index[v,t+q],),1)
    for i in P:
        for j in P:
            component(((0,t+i,t+j),(1,t+i,j),(1,i,t+j)),
                      (index[t+i,j],index[i,t+j]),2)
    for i in B:
        for j in P:leaf((0,t+i,t+j),(1,t+i,j),index[t+i,j],RZ)
    for i in P:
        for j in B:leaf((0,t+i,t+j),(1,i,t+j),index[i,t+j],RY)
    active={v:c for v,c in values.items() if c and v not in removed}
    assert len(fixed)==2*A*(t-m) and len(set(fixed))==len(fixed)
    assert all(set(c)<=RY|RZ for c in active.values())
    EY=tuple(v for v,c in active.items() if set(c)<=RY)
    EZ=tuple(v for v,c in active.items() if set(c)<=RZ)
    H=tuple(v for v in active if v not in EY and v not in EZ)
    assert set(H)=={(0,t+i,t+j) for i in B for j in B}
    assert len(H)==m*m<=4*s*s
    assert len(field_basis(EY,active.__getitem__,p)[0])==len(field_basis(EZ,active.__getitem__,p)[0])==A*m
    assert len(active)==4*A*m-m*m
    S=B|{t+v for v in B}
    assert all(v[0] in (0,1) and 0<=v[1]<a and v[2] in S for v in EY)
    assert all(v[0] in (0,1) and 0<=v[2]<a and v[1] in S for v in EZ)
    COUNTS['exact_active_rank_and_label_inventories']+=1
    return active,tuple(fixed),(EY,EZ,H),m


def plane(q):
    d=gcd(*q);assert d
    q=tuple(x//d for x in q)
    return q if next(x for x in q if x)>0 else tuple(-x for x in q)


def check_patch(data,cols,t,a,s,w,p,cache):
    order=phase_key(t,w);active,fixed,groups,m=passive_reduce(data,cols,t,a,s,w,p)
    signature=(frozenset(fixed),tuple(sorted((v,tuple(sorted(c.items()))) for v,c in active.items())))
    if order in cache:
        assert cache[order]==signature
        COUNTS['independent_same_phase_matrix_stability_checks']+=1
    else:cache[order]=signature
    raw,_=field_basis(sorted(cols,key=lambda v:dot(w,v)),cols.__getitem__,p)
    small,_=field_basis(sorted(active,key=lambda v:dot(w,v)),active.__getitem__,p)
    assert set(raw)==set(fixed)|set(small)
    assert len(small)==2*(a-t)*m
    NN=normals(raw,cols,p);MM=normals(small,active,p);PP=phase_normals(t,order)
    assert all(dot(w,q)>0 for q in NN|MM|PP)
    original=polygon(NN|PP);reduced=polygon(MM|PP)
    assert all(dot(v,q)>=0 for v in original for q in MM)
    assert all(dot(v,q)>=0 for v in reduced for q in NN)
    assert area2(original)==area2(reduced)
    COUNTS['independent_entire_original_vs_active_cone_patches']+=1
    decision=small_decision([None,None,None,None,groups],active,w,p)
    assert decision['bridge_rank']==m*m and decision['small_dual_labels']<=3*m*m
    inventory=set()
    for group in groups:
        for i,v in enumerate(group):
            for u in group[:i]:inventory.add(plane(tuple(x-y for x,y in zip(v,u))))
    inventory|={plane(q) for q in ((1,0,-t),(1,-t,0),(0,1,-1))}
    assert len(inventory)<=512*a*s
    COUNTS['independent_internal_plane_inventories']+=1
    return {'scalar_phase':order,'original_residual_rank':len(raw),'active_rank':len(small),
            'fixed_selected_labels':len(fixed),'active_columns':len(active),
            'active_bridges':m*m,'internal_planes':len(inventory),
            'patch_area2':str(area2(original)),'small_decision':decision}


def whole_fan(data,cols,t,a,s,p,cache):
    n=4*a*t-2*t*t
    queue=[((1,2*n+1,(2*n+1)**2),None)];cells={};patches=[]
    while queue:
        w,old=queue.pop()
        B=frozenset(field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)[0])
        if old is not None:assert B!=old
        if B in cells:continue
        N=normals(B,cols,p);poly=polygon(N)
        assert all(dot(w,q)>0 for q in N)
        cells[B]=(poly,N)
        for left,right in zip(poly,poly[1:]+poly[:1]):
            if any(left[i]==right[i]==0 for i in range(3)):continue
            q=next(q for q in N if dot(q,left)==dot(q,right)==0)
            nw,_=exact_adjacent_weight(n,left,right,q);queue.append((nw,B))
    assert sum(area2(poly) for poly,N in cells.values())==1
    total=Q(0);phase_area=Counter()
    for B,(poly,N) in cells.items():
        parts=Q(0)
        for order in PHASES:
            piece=poly
            for q in phase_normals(t,order):piece=clip(piece,q,1)
            if not piece or not area2(piece):continue
            w=choose_weight(piece,witness_labels(t,cols))
            assert phase_key(t,w)==order
            row=check_patch(data,cols,t,a,s,w,p,cache)
            assert Q(row['patch_area2'])==area2(piece)
            patches.append(row);parts+=area2(piece);phase_area[str(order)]+=area2(piece)
        assert parts==area2(poly);total+=parts
    assert total==1 and len(phase_area)==6
    COUNTS['complete_original_fans_with_all_phase_patches']+=1
    return {'initial_cones':len(cells),'phase_patches':len(patches),'area2':'1',
            'phase_area2':{k:str(v) for k,v in sorted(phase_area.items())},'patches':patches}


def one_case(t,a,s,units,p,kind,complete=False):
    assert 1<=s<=t
    for W in units:
        N={v:c for v,c in W.items() if v!=(0,0)}
        assert all(v>=t-s for q,v in N)
        if 2*s<=t:
            assert not mul(N,N,a-t,t,p)
            COUNTS['independent_square_zero_off_axis_products']+=1
    data,co,cols=assemble(t,a,units,p);cache={};patches=[]
    assert len(field_basis(data[4][2],cols.__getitem__,p)[0])==t*t
    COUNTS['original_bridge_rank_controls']+=1
    for order in PHASES:
        poly=polygon(phase_normals(t,order))
        for seed in (1,19):
            w=choose_weight(poly,witness_labels(t,cols),seed)
            assert phase_key(t,w)==order
            patches.append(check_patch(data,cols,t,a,s,w,p,cache))
    full=whole_fan(data,cols,t,a,s,p,cache) if complete else None
    print('Passive-band case passed',t,a,s,p,kind,'full cones',full['initial_cones'] if full else None,flush=True)
    return {'t':t,'a':a,'s':s,'field_characteristic':p,'kind':kind,'length':4*a*t-2*t*t,
            'original_bridge_rank':t*t,
            'coefficients':list(map(int,co)),'six_phase_patch_controls':patches,'complete_original_fan':full}


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    previous=json.loads((ROOT/'results/research/hook_clipped_threshold_checks.json').read_text())
    assert previous['status']=='passed'
    for f,sha in previous['source_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==sha
    cases=[];rng=Random(20261023)
    for t,a,s,p,kind in [(3,8,1,101,'dense'),(4,12,1,2,'dense'),
                         (6,14,2,101,'incomparable'),(4,10,1,0,'binomial'),
                         (8,17,1,101,'scalar'),(8,17,1,3,'dense'),
                         (5,12,2,2,'dense'),(4,10,3,101,'broad_band')]:
        units=[]
        for arm in (0,1):
            axis=(2+arm) if p not in (2,3) else 1
            W={(0,0):axis}
            if kind in ('dense','broad_band'):
                for v in range(t-s,t):
                    for q in range(a-t):W[q,v]=rng.randrange(p)
            elif kind=='incomparable':
                W[a-t-2,t-2]=1
                for q in range(a-t):W[q,t-1]=comb(a-t-1,q)*(arm+1)
            elif kind=='binomial':
                for q in range(a-t):W[q,t-1]=comb(a-t-1,q)*(arm+1)
            units.append({v:(c%p if p else c) for v,c in W.items() if (c%p if p else c)})
        cases.append(one_case(t,a,s,units,p,kind,complete=(t,a,kind) in ((3,8,'dense'),(8,17,'scalar'))))
    files={'docs/hook_passive_band_fans.md','scripts/run_hook_passive_band_checks.py',
           'docs/hook_bridge_reduction.md','docs/hook_prefix_bridge_rank.md','docs/two_colour_threshold_braids.md',
           'docs/hook_clipped_threshold_budget.md','results/research/hook_clipped_threshold_checks.json'}
    pending=['scripts/run_hook_passive_band_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=[node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else []
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:
                    files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_passive_band_checks.py',
            'scope':'Full positive-fan encoding and transverse-width bound for constant-axis central hooks with zero outer images. No prefix regularity or factorization of the long coefficients is assumed. No new worst-case unrestricted exponent, actual polynomial separation or exact novelty claimed.',
            'counts':dict(sorted(COUNTS.items())),'original_quotient_replay_counts':dict(REPLAY),'cases':cases,
            'uniform_result':{'fan_bound':'O(a^2*s^8)=O(n^2*s^8/t^2)',
                              'active_short_indices':'[0,s) union [t-s,t), m=min(t,2s)',
                              'rank_per_scalar_sector':'2*(a-t)*m <= 4*(a-t)*s',
                              'columns_per_scalar_sector':'4*(a-t)*m-m^2',
                              'possible_mixed_columns':'m^2 <= 4*s^2',
                              'scalar_sectors':6,'internal_order_planes':'O(a*s)',
                              'necessary_width_for_polynomial_excess':'s=Omega(t^(1/4)*n^(epsilon/8))'},
            'limits':['Both outer images are zero and both axis units are constant and nonzero.',
                      'The full parameter bound recovers the previous worst-case order when s is comparable to t.',
                      'The fixed chosen span is a vector contraction, not necessarily an algebra ideal.',
                      'The small decision matrix retains the actual coefficients; no universal sampler of the fan-count size follows.',
                      'Exact priority, sharpness, arbitrary outer images and unrestricted G_3 remain open.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_passive_band_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}))


if __name__=='__main__':main()
