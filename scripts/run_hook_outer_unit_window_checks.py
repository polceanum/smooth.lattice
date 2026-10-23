#!/usr/bin/env python3
"""Full compatible outer maps, exact prefix splits and streamed coherent walls."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import comb, lcm
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_hook_bridge_checks import model, actual, normals, COUNTS as REPLAY
from run_hook_rank_departure_checks import rank, clean, plane, complete_fan, COUNTS as FAN
from run_hook_coherent_unit_window_checks import inverse, transform, lower, upper, scalar, rational_plane, cross
from run_hook_monomial_orbit_checks import common_list
from run_hook_high_price_partner_checks import bareiss, minor_product
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import dot

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def build(t,a,arm,W,p,rng,replay=True):
    data=model(t,a);A=a-t
    W={v:scalar(c,p) for v,c in W.items() if v[0]<A and v[1]<t and scalar(c,p)}
    co=[]
    for g,y,z in data[2]:
        long,short=(y,z) if arm==0 else (z,y)
        co.append(W.get((long-t,short),0) if g==1 and long>=t else rng.randrange(p or 7))
    assert any(c for c,(g,y,z) in zip(co,data[2]) if g==0 and (y,z)[arm]>=t)
    assert any(c for c,(g,y,z) in zip(co,data[2]) if g==2 and (y,z)[arm]>=t)
    cols=actual(data,co,p,replay)
    inv=inverse(W,A,t,p);trans=transform(data,cols,t,a,arm,inv,p)
    for q,j in product(range(A),range(t)):
        v=lower(t,arm,q,j);pair=(t+q,j) if arm==0 else (j,t+q)
        i=data[0].index(pair)
        yy={r:c for r,c in trans[v].items() if data[0][r][arm]>=t}
        assert yy=={i:scalar(1,p)}
        if q>=t:assert trans[v]==yy
        COUNTS['independent_actual_unit_central_source_identities']+=1
    for w in ((1,2,19),(31,7,1)):
        B=field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)[0]
        BB=field_basis(sorted(trans,key=lambda v:(dot(w,v),v)),trans.__getitem__,p)[0]
        assert B==BB and normals(B,cols,p)==normals(BB,trans,p)
        COUNTS['independent_whole_cone_outer_row_chart_equalities']+=1
    return data,co,cols,W,inv,trans


def prefix_controls(data,cols,trans,t,a,arm,D,p):
    B=4*t+2*D;states=[]
    for is_upper,matrix in ((True,cols),(False,trans)):
        for u in (Q(1,3),Q(1),Q(3)):
            M=a//2
            if not (u*M>B*u+4*t and u*(a-M)>B*u+4*t):continue
            prices=((2*t-1+D*u)/2,) if is_upper else (Q(t,3),Q(5*a))
            for lam in prices:
                w=(lam,u,Q(1)) if arm==0 else (lam,Q(1),u)
                e=upper(t,arm,M-t,0) if is_upper else lower(t,arm,M-t,0)
                tau=dot(w,e);assert tau<a*u
                J=set();chosen=[]
                for i,pair in enumerate(data[0]):
                    y,j=pair[arm],pair[1-arm]
                    if y<t or (not is_upper and y<2*t):continue
                    v=upper(t,arm,y-t,j) if is_upper else lower(t,arm,y-t,j)
                    if dot(w,v)<tau:
                        assert matrix[v]=={i:scalar(1,p)}
                        J.add(i);chosen.append(v)
                inside={i for i,pair in enumerate(data[0]) if (t if is_upper else 2*t)<=pair[arm]<a-t}
                projected={v:{i:c for i,c in col.items() if i not in J} for v,col in matrix.items()}
                prefix={v:c for v,c in projected.items() if dot(w,v)<=tau}
                assert projected[e] and set(projected[e])<=inside
                for c in prefix.values():
                    assert not (set(c)&inside) or set(c)<=inside
                outside_prefix_mix=sum(bool(set(c)&inside) and not set(c)<=inside
                                       for v,c in projected.items() if dot(w,v)>tau)
                COUNTS['expensive_mixed_columns_not_claimed_to_split']+=outside_prefix_mix
                assert len(field_basis(chosen,matrix.__getitem__,p)[0])==len(J)
                if len(states)<2:
                    keep=[v for v in matrix if dot(w,v)<=tau]
                    assert rank(tuple(keep),matrix,p)==rank(tuple(keep),projected,p)+len(J)
                    COUNTS['independent_whole_prefix_contraction_rank_equalities']+=1
                COUNTS['independent_entire_cost_prefix_direct_sum_checks']+=1
                states.append({'upper':is_upper,'u':str(u),'lambda':str(lam),'contracted_rank':len(J),
                               'prefix_columns':len(prefix),'expensive_mixed_columns':outside_prefix_mix})
    return states


def admitted(t,a,D,b,d):
    S=t+D
    return b<=16*S or d<=16*S or abs(b-a)<=16*S or abs(b-a)*d<=16*a*t


def inventory(data,t,a,arm,D):
    Y=set(data[4][arm]);Z=set(data[4][1-arm]);labels=tuple(v for g in data[4] for v in g);out=set()
    for i,v in enumerate(labels):
        for f in labels[:i]:
            e,z=(v,f) if v in Y and f in Z else (f,v)
            if e in Y and z in Z:
                b=e[arm+1]-z[arm+1];d=z[2-arm]-e[2-arm]
                if not admitted(t,a,D,b,d):continue
            out.add(plane(tuple(x-y for x,y in zip(v,f))))
    COUNTS['coefficient_independent_outer_plane_inventories']+=1
    return out


def fan_controls(data,cols,t,a,arm,D,p):
    fan=complete_fan(cols,4*a*t-2*t*t,p);planes=inventory(data,t,a,arm,D)
    Y,Z=set(data[4][arm]),set(data[4][1-arm]);central=outer=0
    for basis,poly in fan.items():
        mids=[]
        for left,right in zip(poly,poly[1:]+poly[:1]):
            if any(left[i]==right[i]==0 for i in range(3)):continue
            assert rational_plane(cross(left,right)) in planes
            mids.append(tuple((x+y)/2 for x,y in zip(left,right)))
            COUNTS['complete_original_outer_wall_plane_checks']+=1
        _,piv=field_basis(sorted(basis),cols.__getitem__,p)
        for f in cols:
            if f in basis:continue
            rem,co=field_reduce(cols[f],piv,p);assert not rem
            for e,c in co.items():
                if not c or not ((e in Y and f in Z) or (f in Y and e in Z)):continue
                for w in mids:
                    if dot(w,e)!=dot(w,f):continue
                    y,z=(e,f) if e in Y else (f,e)
                    M,s=y[arm+1],y[2-arm];K,r=z[2-arm],z[arm+1]
                    u=w[arm+1]/w[2-arm];lam=w[0]/w[2-arm]
                    b,d=M-r,K-s;assert admitted(t,a,D,b,d)
                    if not y[0] and M>=a:
                        i=M-a;assert i>=t and (i+1)*u+s<=2*t-1
                        assert K<=a*u+3*t;outer+=1
                    else:
                        if M>=2*t:
                            assert u*M<=(4*t+2*D)*u+4*t or u*(a-M)<=(4*t+2*D)*u+4*t
                            if y[0]:assert lam<=2*t-1+D*u
                            assert K<=u*(M+D)+4*t
                        central+=1
    COUNTS['actual_coherent_central_or_upper_exchange_checks']+=central
    COUNTS['actual_coherent_outer_source_exchange_checks']+=outer
    return {'complete_original_cones':len(fan),'central_or_upper_exchanges':central,
            'outer_source_exchanges':outer,'common_planes':len(planes)}


@lru_cache(None)
def binomial_coefficients(D,p):
    if not p:return tuple(Q(comb(D,i)) for i in range(D+1))
    out=[1]
    for i in range(D):out.append(out[-1]*(D-i)*pow(i+1,p-2,p)%p)
    return tuple(out)


def product_mod(D,s,r,p):
    if not p:return minor_product(D,s,r)
    facts=[1]
    for i in range(1,D+s):facts.append(facts[-1]*i%p)
    out=1
    for i in range(s):
        out=out*facts[D+i]*facts[i]*pow(facts[r+i]*facts[D-r+i]%p,p-2,p)%p
    return out


def witness(m,s,p,original=None):
    t=2;a=8*m*m;A=a-t;K=m;d=A-1-K;D=A-1;r=A-s
    assert m>=3 and m<=s<=2*m and (not p or p>2*A)
    G=binomial_coefficients(D,p)
    def coefficient(i):return G[i] if 0<=i<A else 0
    matrices=[[[int(coefficient(rr+i-j)) for j in range(s)] for i in range(s)] for rr in (r,r-1)]
    determinants=[]
    for matrix,rr in zip(matrices,(r,r-1)):
        det=bareiss(matrix);expected=product_mod(D,s,rr,p)
        assert (det%p==expected!=0 if p else det==expected>0)
        determinants.append(str(det));COUNTS['independent_outer_binomial_integer_minors']+=1
    columns={j:{i:scalar(coefficient(r+i-j),p) for i in range(s) if coefficient(r+i-j)} for j in range(s)}
    _,piv=field_basis(range(s),columns.__getitem__,p);assert len(piv)==s
    ell=[field_reduce({i:scalar(1,p)},piv,p)[1].get(0,0) for i in range(s)]
    alpha=scalar(sum(ell[i]*coefficient(r+i-s) for i in range(s)),p);assert alpha
    source={i:scalar(coefficient(r+i-s),p) for i in range(s)}
    rem,representation=field_reduce(source,piv,p);assert not rem and representation[0]==alpha
    for j in range(s):assert scalar(sum(ell[i]*coefficient(r+i-j) for i in range(s)),p)==int(j==0)
    wx,wy,wz=(4*s+1)*(m+2),2*(m+2),2*a-2*s-3
    w=(wx,wy,wz);u=Q(wy,wz);lam=Q(wx,wz);E=(0,a+s,0);F=(1,1,K+2);tau=dot(w,E)
    assert tau==dot(w,F) and 0<lam<1 and 0<lam-2*u<1 and s*u<1
    yy={(t+r+i,0):ell[i] for i in range(s)}
    zz={(1,K+2):scalar(1,p),(1,a-1):scalar(-1,p)}
    def functional(y,z):return yy.get((y,z),0)+zz.get((y,z),0)
    values={i:scalar(sum(ell[k]*coefficient(r+k-i) for k in range(s)),p) for i in range(s+1)}
    checked=0;potential=0
    def check(v,value):
        nonlocal checked,potential
        potential+=1
        if dot(w,v)<tau:
            assert not scalar(value,p),(m,s,v,value)
            checked+=1
    # Stream all potential residual columns, plus the fixed upper box,
    # without materializing a quadratic-size convolution representation.
    rows=((i,j) for i,j in product(range(t),repeat=2))
    for y,z in rows:check((1,y,z),functional(y,z))
    for q,j in product(range(A),range(t)):
        check((1,t+q,j),functional(t+q,j));check((1,j,t+q),functional(j,t+q))
    for i,j in product(range(A),range(t)):
        v=(0,a+i,j)
        if dot(w,v)<tau:
            assert j==0 and i<s
            value=values[i]+functional(1+i,a-1+j)
        else:value=0
        check(v,value)
    shifts=((i,j) for i,j in product(range(A),range(t)))
    for i,j in shifts:
        check((0,t+i,t+j),functional(t+i,j)+functional(i,t+j)+functional(i,t+d+j))
    for i,j in product(range(t),range(t,A)):
        check((0,t+i,t+j),functional(t+i,j)+functional(i,t+j)+functional(i,t+d+j))
    for i,j in product(range(t),range(A)):check((0,i,a+j),0)
    assert potential==6*A*t
    assert functional(1,K+2)==1 and values[s]==alpha
    # Explicit original-source relation: E-alpha*F is a combination of
    # cheaper outer sources, one central Z source and cheaper upper rows.
    assert all(dot(w,(0,a+i,0))<tau for i in representation)
    assert dot(w,(0,3,K+2))<tau
    correction=0
    for q in range(A):
        value=scalar(coefficient(q-s)-sum(c*coefficient(q-i) for i,c in representation.items()),p)
        if q>=r:assert not value
        elif value:
            assert dot(w,(1,t+q,0))<tau;correction+=1
    if original is not None:
        data,cols=original
        cheap=[v for v in cols if dot(w,v)<tau];B,bp=field_basis(cheap,cols.__getitem__,p)
        assert field_reduce(cols[E],bp,p)[0] and field_reduce(cols[F],bp,p)[0]
        _,fp=field_basis((*B,F),cols.__getitem__,p)
        rem,co=field_reduce(cols[E],fp,p);assert not rem and co.get(F)==alpha
        full_ell={i:scalar(functional(y,z),p) for i,(y,z) in enumerate(data[0]) if functional(y,z)}
        for v,col in cols.items():
            value=scalar(sum(c*full_ell.get(i,0) for i,c in col.items()),p)
            if dot(w,v)<tau:assert not value
        COUNTS['independent_original_field_outer_parallelism_controls']+=1
    b=a+s-1;dd=K+2;assert (s-1)*(m+2)==abs(b-a)*dd
    if abs(b-a)>32 and dd>32:
        assert b>32 and admitted(t,a,0,b,dd)
        COUNTS['actual_outer_walls_outside_all_three_fixed_strips']+=1
    COUNTS['full_streamed_strict_cheaper_column_checks']+=checked
    COUNTS['full_streamed_potential_original_columns']+=potential
    COUNTS['explicit_original_source_parallelism_relations']+=1
    return {'m':m,'s':s,'a':a,'length':4*a*t-2*t*t,'characteristic':p,'integer_weight':w,
            'primitive_direction':[1,1-a-s,K+2],'translated_product':(s-1)*(m+2),
            'strictly_cheaper_columns_checked':checked,'potential_columns_streamed':potential,
            'nonzero_parallelism_factor':str(alpha),'integer_minors':determinants,
            'upper_prefix_correction_columns':correction,'full_original_matrix_replayed':original is not None}


def witness_original(m,p):
    t=2;a=8*m*m;A=a-t;K=m;d=A-1-K;data=model(t,a);G=binomial_coefficients(A-1,p);co=[]
    for g,y,z in data[2]:
        c=0
        if g==0:
            if y>=t and z==0:c=G[y-t]
            elif (y,z)==(1,a-1):c=1
        elif g==1:
            if (y,z)==(t,0):c=1
            elif y==0 and z in (t,t+d):c=1
        co.append(c)
    return data,actual(data,co,p,True)


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    old=json.loads((ROOT/'results/research/hook_coherent_unit_window_checks.json').read_text());assert old['status']=='passed'
    for f,h in old['source_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
    rng=Random(2026102330);cases=[]
    specs=[(2,80,0,{(0,0):1,(1,1):1,(2,1):1},101,False),
           (3,100,1,{(0,0):1,(1,1):1,(2,1):1,(2,2):1},3,False),
           (3,12,0,{(0,0):2,(1,1):1,(2,1):1},0,False),
           (3,6,1,{(0,0):1,(2,1):1},2,False),
           (2,7,0,{(0,0):1,(1,1):1,(2,1):1},101,True),
           (3,8,1,{(0,0):1,(1,1):1,(2,1):1},2,True)]
    for t,a,arm,W,p,full in specs:
        data,co,cols,W,inv,trans=build(t,a,arm,W,p,rng)
        D=max(max(i for i,j in W),max(i for i,j in inv))
        states=prefix_controls(data,cols,trans,t,a,arm,D,p)
        fan=fan_controls(data,cols,t,a,arm,D,p) if full else None
        cases.append({'t':t,'a':a,'arm':arm,'characteristic':p,'D':D,'prefix_checks':states,'full_fan':fan})
        print('Outer unit case passed',t,a,arm,p,'D',D,flush=True)
    witnesses=[]
    for m,p in ((3,1000003),(4,1000003),(3,0)):
        original=witness_original(m,p)
        for s in range(m,2*m+1):witnesses.append(witness(m,s,p,original))
        print('Original outer wall family passed',m,p,flush=True)
    for m in (12,40):
        for s in (m,m+1,2*m):witnesses.append(witness(m,s,1000003))
        print('Streamed outer wall family passed',m,flush=True)
    t=2;a=5;arm=0;D0=2;data=model(t,a);planes=inventory(data,t,a,arm,D0);weights=common_list(planes,a);samplers=[]
    for p,W in ((2,{(0,0):1,(1,1):1,(2,1):1}),(0,{(0,0):2,(0,1):1,(2,1):1})):
        data,co,cols,W,inv,trans=build(t,a,arm,W,p,rng)
        assert max(i for i,j in inv)<=D0
        fan=complete_fan(cols,4*a*t-2*t*t,p);seen=set()
        for w in weights:
            basis=frozenset(field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)[0])
            assert all(dot(w,q)>0 for q in normals(basis,cols,p));seen.add(basis)
        assert seen==set(fan);COUNTS['independent_full_outer_common_sampler_equalities']+=1
        samplers.append({'characteristic':p,'complete_original_cones':len(fan),'weights_reused':len(weights)})
    files={'docs/hook_outer_unit_windows.md','scripts/run_hook_outer_unit_window_checks.py',
           'docs/hook_bridge_reduction.md','docs/hook_coherent_unit_windows.md',
           'results/research/hook_coherent_unit_window_checks.json'};pending=['scripts/run_hook_outer_unit_window_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=([node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else [])
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_outer_unit_window_checks.py',
            'scope':'All compatible outer images, one unit central arm, every field. Prefix separation gives a translated-hyperbola inventory and a full quadratic-polylog closure at bounded propagation. No unrestricted endpoint or exact novelty claimed.',
            'uniform_result':{'D':'max(deg_q W,deg_q W^(-1))','H':'1+log(a/t)',
                              'full_fan':'O(a^2*(D+t*H)^2)=O(n^2*(D/t+H)^2)',
                              'common_integer_weights':'O(a^2*(D0+t*H)^2), height O(a^2)',
                              'central_window':'u*M<=(4t+2D)*u+4t or u*(a-M)<=(4t+2D)*u+4t',
                              'outer_source':'(i+1)*u+j<=2t-1',
                              'direction_inventory':'b<=16(t+D), d<=16(t+D), |b-a|<=16(t+D), or |b-a|*d<=16at',
                              'outer_family':'t=2,a=8m^2, W_Y=1; m+1 coherent E_s partners of one F, translated product Theta(at), char 0 or p>2A'},
            'counts':dict(sorted(COUNTS.items())),'fan_counts':dict(sorted(FAN.items())),
            'original_quotient_replay_counts':dict(REPLAY),'cases':cases,'outer_family_witnesses':witnesses,
            'common_arrangement':{'planes':len(planes),'integer_weights':len(weights),'degree_cap':D0},'common_sampler':samplers,
            'limits':['The chosen central arm is a unit; arbitrary compatible coefficients in both outer images are allowed.',
                      'Only the relevant cost prefix is asserted to split; expensive columns can mix its parts.',
                      'The streamed large witnesses verify an explicit full residual-prefix functional and original source relation, not a complete fan enumeration.',
                      'The translated product allowance is necessary in order; no logarithmic lower bound or polynomial fan excess is claimed.',
                      'Exact strongest-current novelty, unrestricted sharp growth and the original final paper remain unresolved.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_outer_unit_window_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}),flush=True)


if __name__=='__main__':main()
