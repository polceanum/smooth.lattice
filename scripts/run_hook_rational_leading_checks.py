#!/usr/bin/env python3
"""Actual rational source identities, filtered charts and coherent leading rows."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_hook_bridge_checks import model, actual, normals, polygon, COUNTS as REPLAY
from run_hook_rank_departure_checks import clean, plane, rank, complete_fan, COUNTS as FAN
from run_hook_coherent_unit_window_checks import scalar, multiply, inverse, transform, lower, upper, rational_plane, cross
from run_hook_outer_unit_window_checks import admitted as translated_admitted
from run_hook_monomial_orbit_checks import common_list
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import dot, area2

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def pair(P,Qpoly,A,t,p):
    def trim(f):return {v:scalar(c,p) for v,c in f.items() if v[0]<A and v[1]<t and scalar(c,p)}
    P,Qpoly=trim(P),trim(Qpoly)
    assert P and Qpoly.get((0,0))
    dp=max(i+j for i,j in P);dq=max(i+j for i,j in Qpoly)
    assert dp<A and dq<A and P.get((dp,0)) and Qpoly.get((dq,0))
    W=multiply(P,inverse(Qpoly,A,t,p),A,t,p)
    assert multiply(Qpoly,W,A,t,p)==P
    COUNTS['independent_truncated_rational_operator_equalities']+=1
    return P,Qpoly,W,dp,dq


def axis_order(P):return min(i for (i,j),c in P.items() if j==0 and c)


def build(t,a,pairs,p,rng,zero_outer=False,matched_nonunit=True):
    data=model(t,a);spec=[]
    for P,Qpoly in pairs:
        spec.append(pair(P,Qpoly,a-t,t,p) if P is not None else None)
    factors=[None,None];outer_images=[None,None]
    for arm,s in enumerate(spec):
        if s is not None and not s[2].get((0,0)) and matched_nonunit:
            B={} if zero_outer else {(i,j):scalar(rng.randrange(p or 7),p)
                                     for i,j in product(range(a-t),range(t))}
            if not zero_outer:B[0,0]=scalar(1,p)
            factors[arm]=B;outer_images[arm]=multiply(s[2],B,a-t,t,p)
    co=[]
    for g,y,z in data[2]:
        arm=0 if y>=t else 1
        q,j=(y-t,z) if arm==0 else (z-t,y)
        if g==1 and spec[arm] is not None:c=spec[arm][2].get((q,j),0)
        elif g!=1 and zero_outer:c=0
        elif (g==0 and arm==0) or (g==2 and arm==1):
            c=outer_images[arm].get((q,j),0) if outer_images[arm] is not None else rng.randrange(p or 7)
            if spec[arm] is not None and not spec[arm][2].get((0,0)) and not matched_nonunit and (q,j)==(0,0):c=1
        else:c=rng.randrange(p or 7)
        co.append(c)
    cols=actual(data,co,p,True);charts=[]
    for arm,s in enumerate(spec):
        if s is None:charts.append(None);continue
        P,Qpoly,W,dp,dq=s;mat=transform(data,cols,t,a,arm,Qpoly,p)
        if factors[arm] is not None:
            for v in data[4][arm]:
                if v[0] or v[arm+1]<a:continue
                i,j=v[arm+1]-a,v[2-arm]
                expected=multiply(W,multiply(factors[arm],{(i,j):1},a-t,t,p),a-t,t,p)
                got={(data[0][r][arm]-t,data[0][r][1-arm]):c for r,c in cols[v].items()}
                assert got==expected
                COUNTS['independent_nonunit_original_outer_image_factorizations']+=1
            COUNTS['nonunit_arms_with_matching_outer_image_condition']+=1
        for q,j in product(range(a-t),range(t)):
            v=lower(t,arm,q,j)
            expected=multiply(P,{(q,j):1},a-t,t,p)
            got={(data[0][i][arm]-t,data[0][i][1-arm]):c for i,c in mat[v].items() if data[0][i][arm]>=t}
            assert got==expected
            COUNTS['independent_original_central_denominator_clearing_identities']+=1
        for w in ((1,3,17),(31,7,1)):
            B=field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)[0]
            BB=field_basis(sorted(mat,key=lambda v:(dot(w,v),v)),mat.__getitem__,p)[0]
            assert B==BB and normals(B,cols,p)==normals(BB,mat,p)
            COUNTS['independent_whole_normal_cone_denominator_chart_equalities']+=1
        charts.append(mat)
    return data,co,cols,spec,charts


def weights(lam,u,arm):return (lam,u,Q(1)) if arm==0 else (lam,Q(1),u)


def nonunit_outer_controls(data,cols,t,a,arm,spec,p):
    """Axis-power certificate and original late-source reconstruction."""
    P,Qpoly,W,dp,dq=spec;nu=axis_order(P);assert nu>0
    L=a-2*t
    if L<=0:return {'axis_order':nu,'late_sources':0,'early_source_outside':False}
    def trim(f):return {v:c for v,c in f.items() if v[0]<L and v[1]<t}
    P0=trim(P);Q0=trim(Qpoly);W0=trim(W)
    axis={v:c for v,c in P0.items() if v[1]==0}
    R=trim({(i-nu,0):c for (i,j),c in P.items() if j==0});assert R.get((0,0))
    B={v:scalar(-c,p) for v,c in P0.items() if v[1]}
    one={(0,0):scalar(1,p)};ap=[one];bp=[one]
    for k in range(1,t):
        ap.append(multiply(ap[-1],axis,L,t,p));bp.append(multiply(bp[-1],B,L,t,p))
    total=Counter()
    for k in range(t):
        for v,c in multiply(ap[t-1-k],bp[k],L,t,p).items():total[v]+=c
    ir=inverse(R,L,t,p);irpower=one
    for k in range(t):irpower=multiply(irpower,ir,L,t,p)
    certificate=multiply(Q0,multiply(clean(total,p),irpower,L,t,p),L,t,p)
    target={(nu*t,0):scalar(1,p)} if nu*t<L else {}
    assert multiply(W0,certificate,L,t,p)==target
    COUNTS['independent_nonunit_axis_power_image_certificates']+=1
    pure=[v for v in data[4][arm] if not v[0] and v[arm+1]<a and v[2-arm]>=t]
    e0=(0,a+t,0) if arm==0 else (0,0,a+t)
    early=False
    if e0 in cols:
        early=rank(tuple(pure+[e0]),cols,p)==rank(tuple(pure),cols,p)+1
        assert early
        COUNTS['arbitrary_nonunit_early_outer_columns_outside_central_span']+=1
    count=0
    for e in data[4][arm]:
        if e[0] or e[arm+1]<a+t+nu*t:continue
        poly={(data[0][i][arm]-2*t,data[0][i][1-arm]):c for i,c in cols[e].items()}
        assert all(q>=nu*t for q,j in poly)
        rest={(q-nu*t,j):c for (q,j),c in poly.items()}
        H=multiply(certificate,rest,L,t,p)
        assert multiply(W0,H,L,t,p)==poly
        out=Counter();w=weights(Q(1),Q(2),arm)
        for (q,j),c in H.items():
            v=lower(t,arm,t+q,j)
            assert v in pure and dot(w,v)<dot(w,e)
            for i,z in cols[v].items():out[i]+=c*z
        assert clean(out,p)==cols[e]
        count+=1;COUNTS['independent_arbitrary_nonunit_late_outer_source_relations']+=1
    return {'axis_order':nu,'late_sources':count,'early_source_outside':early}


def source_identity(data,cols,t,a,arm,P,Qpoly,dp,dq,is_upper,q,j,w,p):
    k=q-(dp if is_upper else dq)
    assert k>=t
    e=upper(t,arm,q,j) if is_upper else lower(t,arm,q,j)
    terms=Counter()
    for (i,h),c in P.items():terms[upper(t,arm,k+i,j+h)]+=c
    for (i,h),c in Qpoly.items():terms[lower(t,arm,k+i,j+h)]-=c
    terms=clean(terms,p);assert terms.get(e)
    total=Counter();used=0
    for v,c in terms.items():
        col=cols.get(v,{})
        if col and v!=e:
            assert dot(w,v)<dot(w,e),(e,v,w)
            used+=1
        for i,x in col.items():total[i]+=c*x
    assert not clean(total,p)
    COUNTS['independent_original_label_filtered_loop_identities']+=1
    return used


def kept_columns(data,cols,t,a,arm,dp,dq,above):
    drop=set()
    for v in data[4][arm]:
        q=v[arm+1]-t
        if above and v[0] and q>=t+dp:drop.add(v)
        if not above and not v[0] and v[arm+1]<a and v[2-arm]>=t and q>=t+dq:drop.add(v)
    return {v:c for v,c in cols.items() if v not in drop},drop


def open_patch(data,cols,t,a,arm,spec,above,p):
    P,Qpoly,W,dp,dq=spec;u=Q(3000010,2000003);z=t+(dq-dp)*u
    if not above and z<=0:return None
    eps=Q(1,10000019);lam=max(z,0)+1+eps if above else z/2+eps
    w=weights(lam,u,arm)
    if not above:assert lam<z
    assert len({dot(w,v) for v in cols})==len(cols)
    core,dropped=kept_columns(data,cols,t,a,arm,dp,dq,above)
    B=field_basis(sorted(cols,key=lambda v:dot(w,v)),cols.__getitem__,p)[0]
    BB=field_basis(sorted(core,key=lambda v:dot(w,v)),core.__getitem__,p)[0]
    assert B==BB and not set(B)&dropped
    heavy=(0,1,-1) if arm==0 else (0,-1,1)
    price=(1,dp-dq,-t) if arm==0 else (1,-t,dp-dq)
    region={heavy,price if above else tuple(-x for x in price)}
    N=normals(B,cols,p);T=normals(BB,core,p)
    old=polygon(N|region);new=polygon(T|region)
    assert all(dot(v,q)>=0 for v in old for q in T)
    assert all(dot(v,q)>=0 for v in new for q in N)
    assert area2(old)==area2(new)>0
    COUNTS['independent_entire_open_price_region_fan_patches']+=1
    return {'above':above,'deleted_columns':len(dropped),'patch_area2':str(area2(old))}


def chart_controls(data,cols,mat,t,a,arm,spec,p,scalar_case=False):
    P,Qpoly,W,dp,dq=spec;A=a-t;controls=[]
    ratios=(Q(1,7),Q(2,3),Q(3,2),Q(3)) if scalar_case else (Q(3,2),Q(2),Q(5))
    for u in ratios:
        z=t+(dq-dp)*u
        for is_upper in (True,False):
            threshold=dp if is_upper else dq
            q=max(t+threshold,A//2)
            if q>=A:continue
            # In the opposite price region the leading identity deletes this colour.
            lam=max(z,0)+Q(1,3) if is_upper else z/2
            if lam>0:
                w=weights(lam,u,arm)
                for j in (0,t-1):source_identity(data,cols,t,a,arm,P,Qpoly,dp,dq,is_upper,q,j,w,p)
            # In the surviving region, one actual row can make the entire prefix coloop.
            lam=z/2 if is_upper else max(z,0)+Q(1,3)
            if lam<=0:continue
            shift=dq if is_upper else dp
            q=A//2
            if q<t+threshold or q+shift>=A-t:continue
            for j in (0,t-1):
                e=upper(t,arm,q,j) if is_upper else lower(t,arm,q,j);w=weights(lam,u,arm);tau=dot(w,e)
                if tau>=a*u:continue
                target=(t+q+shift,j) if arm==0 else (j,t+q+shift)
                row=data[0].index(target);assert mat[e].get(row)
                prefix=[v for v in cols if dot(w,v)<=tau]
                assert all(not mat[v].get(row) for v in prefix if v!=e)
                if len(controls)<2:
                    rest=[v for v in prefix if v!=e]
                    assert rank(tuple(prefix),cols,p)==rank(tuple(rest),cols,p)+1
                    COUNTS['independent_original_prefix_coloop_rank_equalities']+=1
                COUNTS['independent_entire_prefix_leading_coordinate_certificates']+=1
                controls.append({'upper':is_upper,'u':str(u),'lambda':str(lam),'target_row':target,'prefix_columns':len(prefix)})
    patches=[x for above in (False,True) if (x:=open_patch(data,cols,t,a,arm,spec,above,p)) is not None]
    return controls,patches


def strip_pair(t,a,d,b,c,nu=0):
    R=6*t+2*d+2*t*nu
    return abs(b)<=R or abs(b-a)<=R or abs(c)<=R or abs(c-a)<=R


def inventory(data,t,a,d,scalar_arm=None,nu=0):
    Y,Z=set(data[4][0]),set(data[4][1]);labels=tuple(v for g in data[4] for v in g);out=set()
    for i,v in enumerate(labels):
        for f in labels[:i]:
            e,z=(v,f) if v in Y and f in Z else (f,v)
            if e in Y and z in Z:
                b=e[1]-z[1];c=z[2]-e[2]
                ok=strip_pair(t,a,d,b,c,nu) if scalar_arm is None else translated_admitted(t,a,d,b,c) if scalar_arm==0 else translated_admitted(t,a,d,c,b)
                if not ok:continue
            out.add(plane(tuple(x-y for x,y in zip(v,f))))
    out.add((0,1,-1))
    for delta in range(-d,d+1):
        for arm in (0,1):out.add(plane((1,delta,-t) if arm==0 else (1,-t,delta)))
    COUNTS['coefficient_independent_compact_rational_plane_inventories']+=1
    return out


def fan_controls(data,cols,t,a,spec,p,scalar_arm=None,matched_nonunit=True):
    fan=complete_fan(cols,4*a*t-2*t*t,p)
    d=max(max(s[3:]) for s in spec if s is not None)
    nus=[axis_order(s[0]) if s is not None and not matched_nonunit else 0 for s in spec];nu=max(nus)
    planes=inventory(data,t,a,d,scalar_arm,nu)
    Y,Z=set(data[4][0]),set(data[4][1]);checks=0;price=0;outer=0
    for B,poly in fan.items():
        mids=[]
        for left,right in zip(poly,poly[1:]+poly[:1]):
            if any(left[i]==right[i]==0 for i in range(3)):continue
            assert rational_plane(cross(left,right)) in planes
            mids.append(tuple((x+y)/2 for x,y in zip(left,right)))
            COUNTS['complete_original_rational_wall_plane_checks']+=1
        _,piv=field_basis(sorted(B),cols.__getitem__,p)
        for f in cols:
            if f in B:continue
            rem,co=field_reduce(cols[f],piv,p);assert not rem
            for e,c in co.items():
                if not c or not ((e in Y and f in Z) or (f in Y and e in Z)):continue
                for w in mids:
                    if dot(w,e)!=dot(w,f):continue
                    arm=scalar_arm if scalar_arm is not None else (0 if w[1]>w[2] else 1)
                    if scalar_arm is None and w[1]==w[2]:continue
                    chosen,other=(e,f) if e in data[4][arm] else (f,e)
                    P,Qpoly,W,dp,dq=spec[arm];u=w[arm+1]/w[2-arm];lam=w[0]/w[2-arm];z=t+(dq-dp)*u
                    if lam==z:price+=1;continue
                    M,j=chosen[arm+1],chosen[2-arm];K=other[2-arm]
                    if not chosen[0] and M>=a:
                        i=M-a
                        if nus[arm]:assert u>1 and i<t*max(2,nus[arm]+1)
                        else:assert (i+1)*u+j<=2*t-1
                        if u>1:assert abs(M-a)<=3*t+d+t*nus[arm]
                        outer+=1
                    else:
                        q=M-t
                        if chosen[0] and q>=t+dp:assert lam<z
                        if not chosen[0] and q>=t+dq:assert lam>z
                        shift=dq if chosen[0] else dp
                        if q>=t+(dp if chosen[0] else dq):
                            assert q+shift>=a-2*t or dot(w,chosen)>=a*w[arm+1]
                        if u>1:assert M<=3*t+d or abs(M-a)<=3*t+d
                        else:
                            assert scalar_arm is not None
                            assert u*M<=(4*t+2*d)*u+4*t or u*(a-M)<=(4*t+2*d)*u+4*t
                    b=chosen[arm+1]-other[arm+1];cc=other[2-arm]-chosen[2-arm]
                    assert translated_admitted(t,a,d,b,cc) if scalar_arm is not None else strip_pair(t,a,d,e[1]-f[1] if e in Y else f[1]-e[1],f[2]-e[2] if f in Z else e[2]-f[2],nu)
                    checks+=1
    COUNTS['independent_actual_compact_rational_coherent_exchanges']+=checks
    COUNTS['independent_rational_price_plane_exchange_exceptions']+=price
    COUNTS['independent_actual_rational_outer_exchange_checks']+=outer
    return {'complete_original_cones':len(fan),'coherent_exchange_checks':checks,'price_plane_exceptions':price,'outer_exchange_checks':outer,'common_planes':len(planes)}


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    old=json.loads((ROOT/'results/research/hook_outer_unit_window_checks.json').read_text());assert old['status']=='passed'
    for f,h in old['source_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
    geo=({(0,0):1,(1,0):-1,(0,1):1},{(0,0):1,(1,0):-1})
    two=({(0,0):1,(2,0):1,(1,1):1,(0,2):1},{(0,0):1,(1,0):-1,(0,1):1})
    reverse=({(0,0):2,(1,0):1,(0,1):1},{(0,0):1,(2,0):1,(0,2):1})
    cube=({(0,0):1,(3,0):1,(2,1):1,(1,2):1},{(0,0):1,(2,0):1,(1,1):1,(0,2):1})
    scalar_spec=({(0,0):1,(1,0):1},{(0,0):1,(1,0):-1})
    nonunit=({(1,0):1,(0,1):1},{(0,0):1,(1,0):-1})
    nonunit_two=({(2,0):1,(1,1):1,(0,2):1},{(0,0):1,(1,0):-1})
    rng=Random(2026102333);cases=[]
    configs=[(3,40,(geo,geo),101,False,False,None),
             (3,20,(geo,geo),2,True,False,None),
             (4,30,(two,reverse),101,False,False,None),
             (3,18,(cube,geo),3,False,False,None),
             (3,12,(reverse,geo),0,False,False,None),
             (3,6,(two,geo),2,False,False,None),
             (2,8,(geo,geo),101,False,True,None),
             (3,8,(two,geo),2,False,True,None),
             (2,40,(scalar_spec,(None,None)),101,False,False,0),
             (2,8,(scalar_spec,(None,None)),5,False,True,0),
             (3,18,(nonunit,nonunit),2,True,False,None),
             (3,18,(nonunit,nonunit),101,False,False,None),
             (3,18,(geo,nonunit_two),3,False,False,None),
             (2,7,(nonunit,nonunit),101,False,True,None),
             (3,18,(nonunit,nonunit),2,False,False,None,False),
             (3,18,(nonunit_two,geo),3,False,False,None,False),
             (3,14,(nonunit,nonunit),0,False,False,None,False),
             (2,8,(nonunit,nonunit),101,False,True,None,False)]
    for config in configs:
        t,a,pairs,p,zero_outer,full,scalar_arm=config[:7]
        matched=config[7] if len(config)>7 else True
        data,co,cols,spec,charts=build(t,a,pairs,p,rng,zero_outer,matched);arms=[]
        for arm,s in enumerate(spec):
            if s is None:arms.append(None);continue
            is_scalar=all(j==0 for f in s[:2] for i,j in f)
            certificates,patches=chart_controls(data,cols,charts[arm],t,a,arm,s,p,is_scalar)
            unit=bool(s[2].get((0,0)))
            outer_control=nonunit_outer_controls(data,cols,t,a,arm,s,p) if not unit and not matched else None
            expanded=max(i for i,j in s[2])
            inverse_degree=max(i for i,j in inverse(s[2],a-t,t,p)) if unit else None
            if pairs[arm]==geo:
                assert expanded==inverse_degree==a-t-1
                COUNTS['independent_long_geometric_forward_inverse_controls']+=1
            arms.append({'numerator_degree':s[3],'denominator_degree':s[4],
                         'expanded_forward_degree':expanded,'expanded_inverse_degree':inverse_degree,
                         'scalar_coefficients':is_scalar,'unit':unit,'axis_order':axis_order(s[0]),
                         'matching_outer_restricted':not unit and matched,'arbitrary_outer_certificate':outer_control,
                         'leading_prefix_certificates':certificates,'whole_region_patches':patches})
        fan=fan_controls(data,cols,t,a,spec,p,scalar_arm,matched) if full else None
        cases.append({'t':t,'a':a,'characteristic':p,'zero_outer_images':zero_outer,'scalar_arm':scalar_arm,
                      'matching_nonunit_outer_images':matched,'arms':arms,'full_fan':fan})
        print('Rational leading case passed',t,a,p,'scalar_arm',scalar_arm,flush=True)
    t=2;a=5;d0=1;nu0=1;data=model(t,a);planes=inventory(data,t,a,d0,nu=nu0);weights_list=common_list(planes,a);samplers=[]
    for p,pairs,matched in ((2,(geo,geo),True),(0,(({(0,0):2,(1,0):1,(0,1):1},{(0,0):1,(1,0):-1}),geo),True),
                           (3,(nonunit,nonunit),True),(101,(nonunit,nonunit),False)):
        data,co,cols,spec,charts=build(t,a,pairs,p,rng,matched_nonunit=matched)
        assert max(max(s[3:]) for s in spec)<=d0
        fan=complete_fan(cols,4*a*t-2*t*t,p);seen=set()
        for w in weights_list:
            B=frozenset(field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)[0])
            assert all(dot(w,q)>0 for q in normals(B,cols,p));seen.add(B)
        assert seen==set(fan);COUNTS['independent_compact_rational_common_list_equalities']+=1
        samplers.append({'characteristic':p,'matching_nonunit_outer_images':matched,
                         'complete_original_cones':len(fan),'weights_reused':len(weights_list)})
    files={'docs/hook_rational_leading_windows.md','scripts/run_hook_rational_leading_checks.py',
           'docs/hook_bridge_reduction.md','docs/hook_coherent_unit_windows.md','docs/hook_outer_unit_windows.md',
           'docs/hook_tail_interface.md','results/research/hook_outer_unit_window_checks.json'}
    pending=['scripts/run_hook_rational_leading_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=([node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else [])
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_rational_leading_checks.py',
            'scope':'All-field compact rational leading-coordinate theorem, including arbitrary nonunit outer images with an axis-order allowance. Unit or matching-image arms have the sharper degree-only bound. Scalar rational unit one-arm corollary permits an arbitrary opposite arm. No unrestricted endpoint or exact novelty claimed.',
            'uniform_result':{'hypotheses':'Q has nonzero constant; P,Q have total degrees p,r<A and nonzero pure q^p,q^r coefficients; numerator axis order is finite',
                              'two_arm_full_fan':'O(a^2*(t+d)^2)=O(n^2*(1+d/t)^2) for unit or matching-image arms',
                              'heavy_endpoint_strip':'M<=3t+d or |M-a|<=3t+d for unit/matching arms; add t*nu to the outer band otherwise; price planes are retained',
                              'price_plane':'lambda=t+(r-p)*u',
                              'common_integer_weights':'O(a^2*(t+d0+t*nu0)^2), height O(a^2); nu0=0 for unit or matching-image arms',
                              'scalar_one_arm':'O(n^2*(1+d/t+log(a/t))^2), arbitrary compatible opposite and outer images',
                              'nonunit_extension':'The same quadratic bound holds with Q unit and matching outer components f0Y=W_Y*B_Y, f2Z=W_Z*B_Z on each nonunit arm; cross outer components remain arbitrary.',
                              'arbitrary_nonunit_outer_bound':'O(a^2*(t+d+t*nu)^2)=O(n^2*(1+d/t+nu)^2), nu=max ord_q P(q,0)',
                              'axis_power_certificate':'q^(nu*t) belongs to W*R_L; outer shifts i>=t+nu*t lie in the pure central source span',
                              'geometric_family':'W=1+v/(1-q) has compact degree 1 and expanded forward/inverse degrees A-1; both such arms give a quadratic full fan'},
            'counts':dict(sorted(COUNTS.items())),'fan_counts':dict(sorted(FAN.items())),
            'original_quotient_replay_counts':dict(REPLAY),'cases':cases,
            'common_arrangement':{'planes':len(planes),'integer_weights':len(weights_list),'degree_cap':d0,'axis_order_cap':nu0},'common_sampler':samplers,
            'limits':['Both arms require the compact rational triangle hypotheses for the quadratic theorem; one such arm controls its heavy sector.',
                      'Nonunit arbitrary outer images use the axis-order allowance; matching outer images retain the sharper degree-only count. Unit arms allow every compatible outer coefficient.',
                      'The scalar rational one-arm corollary permits arbitrary opposite coefficients and has the displayed logarithmic allowance.',
                      'Loop identities and prefix coloops retain original labels; denominator multiplication is an invertible row operation, not an algebra coordinate change.',
                      'Small expanded or inverse degree is not assumed; general units need not have compact representations of this form.',
                      'Strongest-current novelty, sharp unrestricted growth, actual polynomial fan separation and the original final paper remain unresolved.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_rational_leading_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}),flush=True)


if __name__=='__main__':main()
