#!/usr/bin/env python3
"""Original-rank primal/dual checks for transverse-prefix hook bounds."""
from collections import Counter,defaultdict
from fractions import Fraction as Q
from itertools import product
from math import comb,factorial,gcd,lcm
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_generic_hook_checks import bareiss
from run_hook_bridge_checks import actual,model,normals,polygon,COUNTS as REPLAY
from run_hook_prefix_uncertainty_checks import scalar_interval
from run_signed_face_checks import field_basis,field_reduce
from run_structural_complexity_checks import dot

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def clean(v,p):
    return {i:c%p if p else Q(c) for i,c in v.items() if (c%p if p else c)}


def mul(U,H,m,t,p):
    out=Counter()
    for (i,j),c in U.items():
        for (k,l),v in H.items():
            if i+k<m and j+l<t:out[i+k,j+l]+=c*v
    return clean(out,p)


def degree(H,u):
    assert H
    return max(i*u+j for i,j in H)


def fixture(t,a,dY,dZ,p,rng,exponents=None):
    data=model(t,a);A=a-t
    ks=exponents or (A-1,A-1)
    U={(0,0):1};V={(0,0):1}
    for out,d,k in ((U,dY,ks[0]),(V,dZ,ks[1])):
        for i in range(A):out[i,d]=comb(k,i) if i<=k else 0
        for i,j in product(range(A),range(d+1,t)):out[i,j]=rng.randrange(p or 101)
    U,V=clean(U,p),clean(V,p)
    f={(t+i,j):c for (i,j),c in U.items()}
    f.update({(j,t+i):c for (i,j),c in V.items()})
    co=[f.get((y,z),0) if g==1 else 0 for g,y,z in data[2]]
    cols=actual(data,co,p,replay=True)
    return data,co,cols,U,V


def gap(F,m,p):
    regular=[]
    for r in range(m+1):
        value=bareiss([[F.get(r+i-j,0) for j in range(m-r)] for i in range(m-r)])
        rank=len(field_basis(range(m-r),lambda j:{i:F.get(r+i-j,0) for i in range(m-r) if F.get(r+i-j,0)},p)[0])
        assert bool(value%p if p else value)==(rank==m-r)
        if rank==m-r:regular.append(r)
        COUNTS['independent_transverse_toeplitz_determinants']+=1
    assert regular[0]==0 and regular[-1]==m
    h=max((v-u-1 for u,v in zip(regular,regular[1:])),default=0)
    return h,regular


def degree_controls():
    rng=Random(2026102267);records=[]
    cases=[(2,8,1,1,0,None),(2,12,1,1,101,None),
           (3,13,1,2,101,None),(4,14,2,3,101,None),
           (3,14,1,1,101,(5,7)),(3,10,1,2,2,(3,4)),
           (4,12,3,2,3,(4,6))]
    for t,a,dY,dZ,p,ks in cases:
        data,co,cols,U,V=fixture(t,a,dY,dZ,p,rng,ks);A=a-t;L=a-2*t
        arms=[]
        for arm,W,d in [('Y',U,dY),('Z',V,dZ)]:
            # Replay each unit on the original commuting quotient columns.
            for long,short in product(range(A),range(t)):
                h=(long,short) if arm=='Y' else (short,long)
                original={}
                for row,c in cols[(0,t+h[0],t+h[1])].items():
                    y,z=data[0][row]
                    if arm=='Y' and y>=t:original[y-t,z]=c
                    elif arm=='Z' and z>=t:original[z-t,y]=c
                assert original==mul(W,{(long,short):1},A,t,p)
                COUNTS['independent_original_rectangular_unit_columns']+=1
            F={i:c for (i,j),c in W.items() if j==d};hs=[]
            for m in sorted({A,L}-{0}):
                h,regular=gap(F,m,p);hs.append(h);lower=m-1-h
                if ks is None:
                    assert h==0
                    for r in range(m):
                        val=Q(1)
                        for j in range(m-r):
                            val*=Q(factorial(A-1+j)*factorial(j),factorial(r+j)*factorial(A-1-r+j))
                        assert val.denominator==1 and (int(val)%p if p else int(val))
                        COUNTS['independent_binomial_transverse_products']+=1
                rows=tuple(product(range(m),range(t)));idx={v:i for i,v in enumerate(rows)}
                vals={v:{idx[w]:c for w,c in mul(W,{v:1},m,t,p).items()} for v in rows}
                for D in range(max(0,lower)):
                    source=tuple(v for v in rows if v[0]<=D)
                    for E in range(lower-D):
                        rank=len(field_basis(source,lambda v:{i:c for i,c in vals[v].items() if rows[i][0]>E},p)[0])
                        assert len(source)-rank==d*(min(D,E)+1)
                        COUNTS['all_vector_exceptional_tail_kernel_dimensions']+=1
                for u in [Q(1,3),Q(1),Q(7,3)]:
                    threshold=u*lower
                    for cap in sorted({degree({v:1},u) for v in rows if degree({v:1},u)<threshold}):
                        source=tuple(v for v in rows if degree({v:1},u)<=cap)
                        outcap=threshold-cap-Q(1,10007)
                        rank=len(field_basis(source,lambda v:{i:c for i,c in vals[v].items() if degree({rows[i]:1},u)>outcap},p)[0])
                        expected=sum(v[1]>=t-d and degree({v:1},u)<=min(cap,outcap) for v in source)
                        assert len(source)-rank==expected
                        COUNTS['all_vector_weighted_exceptional_tail_kernels']+=1
                # Matrix transpose versus a separate reflected multiplication.
                for i,j in product(range(m),range(t)):
                    K={(m-1-i,t-1-j):1};WK=mul(W,K,m,t,p)
                    for h0,h1 in product(range(m),range(t)):
                        expected=vals[h0,h1].get(idx[i,j],0)
                        assert WK.get((m-1-h0,t-1-h1),0)==expected
                        COUNTS['independent_reflected_transpose_matrix_entries']+=1
                arms.append({'arm':arm,'module_length':m,'first_transverse_order':d,'regular_indices':regular,'maximum_gap':h,'conditional_degree_sum_lower_bound':lower})
        records.append({'t':t,'a':a,'field_characteristic':p,'coefficients':list(map(int,co)),'arms':arms})
    return records


def wall_weight(left,right,labels,q):
    for k in range(1,2000):
        w=tuple(x+k*y for x,y in zip(left,right))
        if min(w)<=0:continue
        tied=defaultdict(list)
        for v in labels:tied[dot(w,v)].append(v)
        coordinate=next(i for i,c in enumerate(q) if c)
        if all(all(all((v[i]-g[0][i])*q[coordinate]==q[i]*(v[coordinate]-g[0][coordinate]) for i in range(3)) for v in g[1:]) for g in tied.values()):
            return tuple(x/w[2] for x in w)
    raise AssertionError('No generic facet point found')


def relation_control(data,cols,U,V,t,a,dY,dZ,p,w,e,f,hs):
    lam,u,_=w;Y,Z=e[1:];R,S=f[1:];A=a-t;L=a-2*t
    tau=dot(w,e);assert tau==dot(w,f)
    cheap=tuple(v for v in cols if dot(w,v)<tau)
    W,piv=field_basis(cheap,cols.__getitem__,p)
    re,_=field_reduce(cols[e],piv,p);rf,_=field_reduce(cols[f],piv,p)
    if not re or not rf:return None
    _,pf=field_basis((*W,f),cols.__getitem__,p)
    rem,rep=field_reduce(cols[e],pf,p)
    if rem:return None
    assert rep.get(f)
    H=Counter({(Y-t,Z-t):1})
    for v,c in rep.items():
        if not v[0]:H[v[1]-t,v[2]-t]-=c
    H=clean(H,p)
    HY={h:c for h,c in H.items() if h[0]<A and h[1]<t}
    HZ={(j,i):c for (i,j),c in H.items() if i<t and j<A}
    UHY,VHZ=mul(U,HY,A,t,p),mul(V,HZ,A,t,p)
    D=(Y-t)*u+Z-t;Q0=R*u+S;F=(Y+R)*u+S+Z
    assert degree(HY,u)==D and degree(UHY,u)<Q0-t*u
    assert degree(HZ,1/u)*u<D and degree(VHZ,1/u)*u==Q0-t
    assert min(j for i,j in HY)<t-dY and min(j for i,j in HZ)<t-dZ
    COUNTS['actual_primal_nonscalar_tail_certificates']+=2
    # A coordinate in an extended original-column basis is the separating dual functional.
    upper=tuple(v for v in cols if v[0])
    B,piv=field_basis((*W,e,*upper),cols.__getitem__,p)
    assert len(B)==2*t*A
    ell={data[0].index(v[1:]):field_reduce(cols[v],piv,p)[1].get(e,0) for v in upper}
    ev=lambda v:sum(c*ell.get(i,0) for i,c in cols[v].items())%p
    assert all(ev(v)==0 for v in W) and ev(e)==1 and ev(f)!=0
    KY={};KZ={}
    for row,c in ell.items():
        if not c:continue
        i,j=data[0][row]
        if i>=t:KY[A-1-(i-t),t-1-j]=c
        else:KZ[A-1-(j-t),t-1-i]=c
    assert all(i<L for i,j in KY) and all(i<L for i,j in KZ)
    UKY,VKZ=mul(U,KY,L,t,p),mul(V,KZ,L,t,p)
    assert degree(KY,u)<(a-1)*u+t-1-Q0
    assert degree(UKY,u)==(A-1)*u+t-1-D
    assert degree(KZ,1/u)*u==(t-1)*u+a-1-Q0
    assert degree(VKZ,1/u)*u<(t-1)*u+A-1-D
    assert min(j for i,j in KY)<t-dY and min(j for i,j in KZ)<t-dZ
    COUNTS['actual_dual_nonscalar_tail_certificates']+=2
    hY,hZ=hs
    assert F>(a+t-1-hY)*u+t and F>a+t*u+t-1-hZ
    assert F<(a+2*t-1+hY)*u+3*t-2 and F<a+(3*t-2)*u+2*t-1+hZ
    COUNTS['actual_original_primal_dual_four_inequalities']+=1
    # Two strict adjacent original minima independently witness this wall.
    costs=sorted({dot(w,v) for v in cols});eps=min(lam/2,min(y-x for x,y in zip(costs,costs[1:]))/4)
    bases=[]
    for sign in (-1,1):
        wp=(lam+sign*eps,u,Q(1));assert min(wp)>0
        bb,_=field_basis(sorted(cols,key=lambda v:(dot(wp,v),v)),cols.__getitem__,p)
        nn=normals(bb,cols,p);assert all(dot(wp,q)>0 for q in nn)
        bases.append(set(bb));COUNTS['actual_original_adjacent_strict_minima']+=1
    assert f in bases[0] and e not in bases[0] and e in bases[1] and f not in bases[1]
    return {'weight':list(map(str,w)),'lower':e,'upper':f,'direction':tuple(x-y for x,y in zip(f,e)),
            'primal_transverse_valuations':[min(j for i,j in HY),min(j for i,j in HZ)],
            'dual_transverse_valuations':[min(j for i,j in KY),min(j for i,j in KZ)],
            'primal_relation_terms':len(H),'dual_nonzero_rows':sum(bool(c) for c in ell.values())}


def wall_controls():
    rng=Random(2026102269);records=[]
    for t,a,dY,dZ in [(2,12,1,1),(3,30,1,1),(4,32,2,1)]:
        p=101;data,co,cols,U,V=fixture(t,a,dY,dZ,p,rng)
        seen=set();done=set();witnesses=[]
        for u,fac in product([Q(4,5),Q(9,10),Q(19,20),Q(21,20),Q(11,10),Q(5,4)], [Q(1,3),Q(4,3),Q(3),Q(6)]):
            w=(fac*t*min(u,1)+Q(1,100003),u,Q(1))
            B,_=field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)
            if frozenset(B) in seen:continue
            seen.add(frozenset(B));poly=polygon(normals(B,cols,p));NN=normals(B,cols,p)
            for left,right in zip(poly,poly[1:]+poly[:1]):
                if any(left[i]==right[i]==0 for i in range(3)):continue
                q=next(q for q in NN if dot(q,left)==dot(q,right)==0)
                if q[0]<0:q=tuple(-x for x in q)
                if not(q[0]==1 and q[1]<0<q[2]) or q in done:continue
                ww=wall_weight(left,right,cols,q);lam,uu,_=ww
                for e in cols:
                    if e[0] or not(2*t<=e[1]<a and t<=e[2]<2*t):continue
                    f=tuple(x+y for x,y in zip(e,q))
                    if f not in cols or not(0<=f[1]<t and 2*t<=f[2]<a):continue
                    Q0=f[1]*uu+f[2]
                    if Q0<=max((2*t-1)*uu+t-1,(t-1)*uu+2*t-1):continue
                    witness=relation_control(data,cols,U,V,t,a,dY,dZ,p,ww,e,f,(0,0))
                    if witness:
                        witnesses.append(witness);done.add(q);break
                if len(witnesses)>=8:break
            if len(witnesses)>=8:break
        records.append({'t':t,'a':a,'first_transverse_orders':[dY,dZ],'field_characteristic':p,
                        'coefficients':co,'sampled_original_cones':len(seen),'actual_wall_certificates':witnesses,
                        'scope':'Actual original-field wall and dual-functional controls; no complete fan count is inferred.'})
    assert sum(len(v['actual_wall_certificates']) for v in records)
    return records


def scalar_controls():
    rng=Random(2026102271);records=[]
    for t,hY,hZ in [(2,0,0),(3,4,2),(4,2,7)]:
        T=t+hY+hZ;a=512*T+1;total=feasible=0
        for _ in range(3000):
            b=rng.randrange(a//4,3*a//4);c=a-b+rng.randrange(-20*T,20*T+1)
            R,Z=rng.randrange(t),rng.randrange(t,2*t);Y,S=b+R,c+Z
            if not(0<c and Y<a and S<a and b*c>64*a*T):continue
            total+=1;K,J=Y+R,S+Z
            inequalities=[(b,-c,True),(2*t-1-Y,a-1-Z,False),
                          (K-a-t+1+hY,J-t,True),(K-t,J-a-t+1+hZ,True),
                          (a+2*t-1+hY-K,3*t-2-J,True),(3*t-2-K,a+2*t-1+hZ-J,True)]
            interval=scalar_interval(inequalities);COUNTS['independent_scalar_necessary_interval_eliminations']+=1
            if interval is None:continue
            lo,hi=interval;u=(lo+hi)/2;F=K*u+J
            assert Q(a-2*t+1-hZ,a+t-1+hY)<=lo<hi<=Q(a+t-1+hZ,a-2*t+1-hY)
            assert a*abs(u-1)<4*T and a-T<F<a+6*T
            assert abs(F-b-c)<11*T and abs(b+c-a)<20*T
            assert R*u+S>max((2*t-1)*u+t-1,(t-1)*u+2*t-1)
            feasible+=1;COUNTS['nonvacuous_scalar_strip_and_overlap_controls']+=1
        assert total and feasible
        records.append({'t':t,'a':a,'hY':hY,'hZ':hZ,'candidate_intervals':total,'feasible_necessary_intervals':feasible,
                        'scope':'Exact scalar controls only, not coherent wall constructions.'})
    for t,hY,hZ in [(2,0,0),(3,4,2),(4,2,7)]:
        T=t+hY+hZ;a=512*T+1;feasible=0
        for _ in range(3000):
            b=rng.randrange(a//4,3*a//4);c=b+rng.randrange(-12*T,12*T+1)
            if not(0<c<a and b*c>64*a*T):continue
            u=Q(c,b)
            lower=max((a+t-1-hY)*u+t,a+t*u+t-1-hZ)
            upper=min((a+2*t)*u+4*t,a+4*t*u+2*t)
            COUNTS['independent_same_layer_small_threshold_scalar_candidates']+=1
            if lower>=upper:continue
            assert Q(a-3*t-1-hZ,a+t)<u<Q(a+t,a-3*t-1-hY)
            assert a*abs(u-1)<6*T and abs(b-c)<6*T
            feasible+=1;COUNTS['nonvacuous_same_layer_small_threshold_scalar_controls']+=1
        assert feasible
        records.append({'branch':'same-layer-small-upper-threshold','t':t,'a':a,'hY':hY,'hZ':hZ,
                        'feasible_necessary_intervals':feasible,'scope':'Necessary scalar compatibility, not actual coherent wall claims.'})
    return records


def same_layer_relation(data,cols,U,V,t,a,dY,dZ,p,w,e,f):
    lam,u,_=w;A=a-t;L=a-2*t;tau=dot(w,e)
    assert tau==dot(w,f) and e[0]==f[0]
    Q0=tau-lam;D=tau-t*(u+1);F=tau+Q0
    cheap=tuple(v for v in cols if dot(w,v)<tau)
    W,piv=field_basis(cheap,cols.__getitem__,p)
    re,_=field_reduce(cols[e],piv,p);rf,_=field_reduce(cols[f],piv,p)
    if not re or not rf:return None
    _,pivf=field_basis((*W,f),cols.__getitem__,p)
    rem,rep=field_reduce(cols[e],pivf,p)
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
    HY={h:c for h,c in H.items() if h[0]<A and h[1]<t}
    HZ={(j,i):c for (i,j),c in H.items() if i<t and j<A}
    P=[degree(HY,u),degree(mul(U,HY,A,t,p),u),
       degree(HZ,1/u)*u,degree(mul(V,HZ,A,t,p),1/u)*u]
    if not e[0]:assert P[0]==D and P[2]==D and P[1]<Q0-t*u and P[3]<Q0-t
    else:assert P[0]<D and P[2]<D and P[1]==Q0-t*u and P[3]==Q0-t
    assert min(j for i,j in HY)<t-dY and min(j for i,j in HZ)<t-dZ
    assert F>(a+t-1)*u+t and F>a+t*u+t-1
    COUNTS['actual_same_layer_primal_uncertainty_certificates']+=1
    branch='cheap-overlap'
    if Q0<=max((2*t-1)*u+t-1,(t-1)*u+2*t-1):
        assert not e[0]
        assert F<(a+2*t)*u+4*t and F<a+4*t*u+2*t
        branch='small-upper-threshold'
        COUNTS['actual_same_layer_small_threshold_certificates']+=1
    else:
        upper=tuple(v for v in cols if v[0])
        B,fullpiv=field_basis((*W,e,*upper),cols.__getitem__,p);assert len(B)==2*t*A
        ell={data[0].index(v[1:]):field_reduce(cols[v],fullpiv,p)[1].get(e,0) for v in upper}
        KY={};KZ={}
        for row,c in ell.items():
            if not c:continue
            i,j=data[0][row]
            if i>=t:KY[A-1-(i-t),t-1-j]=c
            else:KZ[A-1-(j-t),t-1-i]=c
        assert all(i<L for i,j in KY) and all(i<L for i,j in KZ)
        R=[degree(KY,u),degree(mul(U,KY,L,t,p),u),
           degree(KZ,1/u)*u,degree(mul(V,KZ,L,t,p),1/u)*u]
        bounds=[(a-1)*u+t-1-Q0,(A-1)*u+t-1-D,
                (t-1)*u+a-1-Q0,(t-1)*u+A-1-D]
        if not e[0]:assert R[0]<bounds[0] and R[1]==bounds[1] and R[2]<bounds[2] and R[3]==bounds[3]
        else:assert R[0]==bounds[0] and R[1]<bounds[1] and R[2]==bounds[2] and R[3]<bounds[3]
        assert min(j for i,j in KY)<t-dY and min(j for i,j in KZ)<t-dZ
        assert F<(a+2*t-1)*u+3*t-2 and F<a+(3*t-2)*u+2*t-1
        COUNTS['actual_same_layer_dual_uncertainty_certificates']+=1
    values=sorted({dot(w,v) for v in cols})
    eps=min(u/2,min(y-x for x,y in zip(values,values[1:]))/(8*a))
    bases=[]
    for sign in (-1,1):
        wp=(lam,u+sign*eps,Q(1))
        B,_=field_basis(sorted(cols,key=lambda v:(dot(wp,v),v)),cols.__getitem__,p)
        assert all(dot(wp,q)>0 for q in normals(B,cols,p))
        bases.append(set(B))
        COUNTS['actual_same_layer_adjacent_strict_minima']+=1
    if not(e in bases[0] and f not in bases[0] and f in bases[1] and e not in bases[1]):return None
    return {'weight':list(map(str,w)),'lower_cost_side_label':e,'other_label':f,'layer':e[0],'branch':branch}


def same_layer_controls(records):
    out=[]
    for record in records[1:]:
        t,a=record['t'],record['a'];dY,dZ=record['first_transverse_orders'];p=record['field_characteristic']
        data=model(t,a);co=record['coefficients'];cols=actual(data,co,p)
        U={(i-t,j):c for (g,i,j),c in zip(data[2],co) if g==1 and i>=t and c}
        V={(j-t,i):c for (g,i,j),c in zip(data[2],co) if g==1 and j>=t and c}
        witnesses=[];seen=set();done=set()
        for u,fac in product([Q(4,5),Q(9,10),Q(19,20),Q(21,20),Q(11,10),Q(5,4)], [Q(1,3),Q(4,3),Q(3),Q(6)]):
            w=(fac*t*min(u,1)+Q(1,100003),u,Q(1))
            B,_=field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)
            if frozenset(B) in seen:continue
            seen.add(frozenset(B));NN=normals(B,cols,p);poly=polygon(NN)
            # Primitive integer representatives for a full-dimensional original cone.
            rays=[]
            for v in poly:
                den=lcm(*(q.denominator for q in v));r=tuple(int(q*den) for q in v);g=gcd(*r);r=tuple(q//g for q in r)
                assert max(r)<=8*a*a and r[1]<=4*a and r[2]<=4*a
                rays.append(r)
            for triple in __import__('itertools').combinations(rays,3):
                if bareiss([list(v) for v in triple]):
                    ww=tuple(sum(v[i] for v in triple) for i in range(3))
                    assert min(ww)>0 and all(dot(ww,q)>0 for q in NN)
                    COUNTS['independent_integer_cone_weight_height_certificates']+=1
                    break
            else:raise AssertionError('Full cone needs independent extreme rays')
            for left,right in zip(poly,poly[1:]+poly[:1]):
                if any(left[i]==right[i]==0 for i in range(3)):continue
                q=next(q for q in NN if dot(q,left)==dot(q,right)==0)
                if q[0] or not q[1] or not q[2]:continue
                if q[1]>0:q=tuple(-x for x in q)
                if not q[2]>0 or q in done:continue
                ww=wall_weight(left,right,cols,q)
                for e in cols:
                    if not(2*t<=e[1]<a and (t<=e[2]<2*t if not e[0] else e[2]<t)):continue
                    # A same-layer primitive direction may scale an exchange.
                    for scale in range(1,a):
                        f=tuple(x+scale*y for x,y in zip(e,q))
                        if f not in cols or f[0]!=e[0] or not(2*t<=f[2]<a and (t<=f[1]<2*t if not f[0] else f[1]<t)):continue
                        witness=same_layer_relation(data,cols,U,V,t,a,dY,dZ,p,ww,e,f)
                        if witness:witnesses.append(witness);done.add(q);break
                    if q in done:break
                if len(witnesses)>=6:break
            if len(witnesses)>=6:break
        out.append({'t':t,'a':a,'actual_same_layer_walls':witnesses})
    assert sum(len(v['actual_same_layer_walls']) for v in out)
    return out


def scalar_arm_controls():
    rng=Random(2026102273);records=[]
    for t,a,p,scalar_arm in [(2,9,0,'Y'),(3,12,2,'Y'),(3,12,101,'Z')]:
        data=model(t,a);A=a-t
        f={(t,0):1,(0,t):1}
        for i,j in product(range(1,a),repeat=2):
            if (i<t or j<t) and (i>=t or j>=t):
                arm='Y' if i>=t else 'Z'
                if arm!=scalar_arm:f[i,j]=rng.randrange(p or 101)
        co=[f.get((i,j),0) if g==1 else 0 for g,i,j in data[2]]
        cols=actual(data,co,p,replay=True);pairs=[];removed=set()
        for i,j in product(range(t,A),range(t)):
            row=(t+i,j) if scalar_arm=='Y' else (j,t+i)
            h=(i,j) if scalar_arm=='Y' else (j,i)
            upper=(1,*row);plain=(0,t+h[0],t+h[1]);r=data[0].index(row)
            assert cols[upper]==cols[plain]=={r:1}
            assert {v for v,c in cols.items() if c.get(r)}=={upper,plain}
            pairs.append((upper,plain));removed.update((upper,plain))
            COUNTS['independent_scalar_arm_direct_parallel_summands']+=1
        coordinate=1 if scalar_arm=='Y' else 2
        remaining={v:c for v,c in cols.items() if c and v not in removed}
        assert all(v[coordinate]<2*t for v in remaining)
        for lam,u in [(Q(t,2),Q(1)),(Q(3*t),Q(1)),(Q(3*t,2),Q(2))]:
            w=(lam,u,Q(1));fixed={min(pair,key=lambda v:(dot(w,v),v)) for pair in pairs}
            B,_=field_basis(sorted(remaining,key=lambda v:(dot(w,v),v)),remaining.__getitem__,p)
            full,_=field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)
            assert set(B)|fixed==set(full)
            assert all(dot(w,q)>=0 for q in normals(full,cols,p))
            COUNTS['independent_scalar_arm_original_minimum_decompositions']+=1
        records.append({'t':t,'a':a,'field_characteristic':p,'scalar_arm':scalar_arm,'direct_parallel_summands':len(pairs),'remaining_ground_labels':len(remaining)})
    return records


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    fields=degree_controls();print('degree and transpose controls passed',flush=True)
    walls=wall_controls();print('actual primal/dual wall controls passed',flush=True)
    same_layers=same_layer_controls(walls);print('same-layer controls passed',flush=True)
    scalars=scalar_controls();scalar_arms=scalar_arm_controls()
    files={'scripts/run_hook_transverse_prefix_checks.py','docs/hook_transverse_prefix.md',
           'docs/hook_bridge_reduction.md','docs/hook_central_skeleton.md','docs/hook_prefix_gaps.md',
           'docs/cross_layer_direction_budget.md'};pending=['scripts/run_hook_transverse_prefix_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=[node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else []
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_transverse_prefix_checks.py',
            'scope':'Uniform all-direction and fan bounds for central-only constant-axis hooks from first transverse prefix gaps, including all positive weight sectors. A common complete bounded-height weight list follows on the regular locus. No unrestricted improvement or polynomial lower separation.',
            'counts':dict(sorted(COUNTS.items())),'original_quotient_replay_counts':dict(REPLAY),
            'field_controls':fields,'actual_wall_controls':walls,'same_layer_controls':same_layers,'scalar_controls':scalars,'scalar_arm_controls':scalar_arms,
            'limits':['Outer generator images are zero.','Nonunit first transverse coefficients are not covered.','All-vector uncertainty excludes explicitly checked scalar transverse tails.','The same-layer small-upper-threshold branch has scalar checks and a uniform proof; no original-field witness for that branch is claimed.','Exact priority and the original unrestricted goal remain unresolved.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_transverse_prefix_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}))


if __name__=='__main__':main()
