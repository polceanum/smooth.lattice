#!/usr/bin/env python3
"""Exact bridge-intersection flags and unrestricted low-wedge partners."""
from collections import Counter
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import ast
import hashlib
import json

from run_hook_weighted_partner_checks import own_sectors,reduce_plain,add_plain
from run_hook_filtered_annihilator_checks import assemble,dense_nullspace,combination
from run_hook_transverse_prefix_checks import clean,mul,wall_weight
from run_hook_monomial_perturbation_checks import reflected
from run_hook_bridge_checks import normals,polygon,COUNTS as REPLAY
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_signed_face_checks import field_basis,field_reduce
from run_structural_complexity_checks import dot,clip,area2

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def flag_states(data,cols,t,a,E,J,p):
    eps,Y,Z=E;left,right=J;u=(left+right)/2 if right is not None else left+1
    planar=lambda v:v[1]*u+v[2]
    own=lambda v:(v[1]>=t and v[2]<t) if v[0] else (2*t<=v[1]<a and t<=v[2]<2*t)
    LY=tuple(v for v in cols if v[0]==eps and own(v) and planar(v)<planar(E))
    flag=tuple(sorted((v for v in cols if v[0]!=eps and own(v)),key=lambda v:(planar(v),v)))
    bridges=tuple(v for v in cols if not v[0] and t<=v[1]<2*t and t<=v[2]<2*t)
    _,lyp=field_basis(LY,cols.__getitem__,p)
    yrows={i for i,(y,z) in enumerate(data[0]) if y>=t}
    bv={i:field_reduce({r:c for r,c in cols[v].items() if r in yrows},lyp,p)[0]
        for i,v in enumerate((*bridges,E))}
    d=len(field_basis(tuple(bv),bv.__getitem__,p)[0])
    assert d<=t*t+1
    fp={};bp={}
    for v in bv.values():add_plain(v,bp,p)
    critical=[];dimension=0;retained={i:dict(v) for i,v in bv.items()}
    for j,v in enumerate(flag,1):
        vv=field_reduce(cols[v],lyp,p)[0]
        add_plain(vv,fp,p);add_plain(vv,bp,p)
        now=d+len(fp)-len(bp)
        assert dimension<=now<=dimension+1
        reduced={i:reduce_plain(vv0,fp,p) for i,vv0 in bv.items()}
        rank,kernel=dense_nullspace(reduced,tuple(bv),p)
        if now>dimension:
            critical.append(j);retained=reduced;dimension=now
        else:
            assert rank==len(field_basis(tuple(retained),retained.__getitem__,p)[0])
            for h in kernel:assert not combination(retained,h,p)
            COUNTS['independent_omitted_prefix_kernel_equalities']+=1
        COUNTS['original_intersection_flag_prefixes']+=1
    assert len(critical)<=d
    states=[]
    for i,pos in enumerate([0,*critical]):
        fixed=(*LY,*flag[:pos],*bridges)
        _,qp=field_basis(fixed,cols.__getitem__,p)
        values={v:field_reduce(c,qp,p)[0] for v,c in cols.items()}
        if not values[E]:continue
        zp={}
        for v in cols:
            if v[0] and v[1]<t:add_plain(values[v],zp,p)
        if reduce_plain(values[E],zp,p):continue
        C=tuple(v for v in cols if v[0]==eps and v[1]<2*t and values[v]
                and ((not eps and t<=v[1] and v[2]>=2*t) or (eps and v[1]<t)))
        D=tuple(v for v in cols if v[0]!=eps and
                ((not v[0] and t<=v[1]<2*t and v[2]>=2*t) or (v[0] and v[1]<t)))
        low=flag[pos-1] if pos else None
        nxt=critical[i] if i<len(critical) else None
        high=flag[nxt-1] if nxt is not None else None
        states.append({'fixed':fixed,'values':values,'C':C,'D':D,'low':low,'high':high,'prefix':pos})
    assert len(states)<=t*t+2
    COUNTS['critical_bridge_flag_states']+=len(states)
    return states


def state_frontier(state,t,a,E,J,kind,p):
    eps,Y,Z=E;values=state['values'];C,D=state['C'],state['D']
    beta=(1+eps)*t-1
    lo=J[0];cap=Q(a,Y-(t if kind=='lower_low' else 2*t)+1)
    hi=min(J[1],cap) if J[1] is not None else cap
    if lo>=hi:return []
    entries=sorted((Q(v[2]-Z,Y-v[1]),v) for v in C)
    cuts={lo,hi};cuts.update(q for q,v in entries if lo<q<hi)
    cuts.update(Q(i,j) for i in range(1,a) for j in range(1,t) if lo<Q(i,j)<hi)
    cuts=sorted(cuts);cp={};at=0;last=None;born=False;dead=False;active=[]
    for left,right in zip(cuts,cuts[1:]):
        u=(left+right)/2
        while at<len(entries) and entries[at][0]<u:
            add_plain(values[entries[at][1]],cp,p);at+=1
        if not reduce_plain(values[E],cp,p):dead=True;continue
        assert not dead
        pp=dict(cp);F=None
        for v in sorted(D,key=lambda v:(v[1]*u+v[2],v)):
            if add_plain(values[v],pp,p) and not reduce_plain(values[E],pp,p):
                F=v;break
        if F is None:assert not born;continue
        born=True;b=F[1]*u+F[2]
        if last:
            assert b<=last[1]+beta*(u-last[0])
            COUNTS['independent_barrier_threshold_inequalities']+=1
        last=(u,b);R,S=F[1:]
        if kind=='lower_low':
            aa=max(left,Q(S-Z,Y-R));bb=min(right,Q(S-Z+t,Y-R))
        elif kind=='lower_middle':
            aa=max(left,Q(S-Z+t,Y-R));bb=min(right,Q(S-Z,Y-R-t))
        else:
            aa=max(left,Q(S-Z,Y-R+t));bb=min(right,Q(S-Z,Y-R))
        if state['low']:
            q=state['low'];assert q[1]>R
            bb=min(bb,Q(S-q[2],q[1]-R))
        if state['high']:
            q=state['high'];assert q[1]>R
            aa=max(aa,Q(S-q[2],q[1]-R))
        if aa<bb:active.append({'lo':aa,'hi':bb,'partner':F,'state':state})
        COUNTS['whole_barrier_frontier_cells']+=1
    bound=4*t*t if kind=='lower_low' else Q(16*a*t*t,Y)
    assert len(active)<=bound
    if active:
        aa=min(v['lo'] for v in active);bb=max(v['hi'] for v in active)
        maxwidth=Q(t,Y-t+1) if kind=='lower_low' else Q(t*a,(Y-2*t+1)**2)
        assert bb-aa<=maxwidth
        COUNTS['whole_barrier_active_window_bounds']+=1
    return active


def verify_original(data,cols,t,a,E,record,J,kind,p,strict=False):
    F=record['partner'];u=(record['lo']+record['hi'])/2
    lam=(E[1]-F[1])*u+E[2]-F[2]
    if E[0]:lam=-lam
    q=tuple(f-e for e,f in zip(E,F))
    if q[0]<0:q=tuple(-v for v in q)
    w=(lam,u,Q(1));tau=dot(w,E);assert tau==dot(w,F)
    assert 0<lam<t*u and u>1
    cheap=tuple(v for v in cols if dot(w,v)<tau)
    B,pp=field_basis(cheap,cols.__getitem__,p)
    assert field_reduce(cols[E],pp,p)[0] and field_reduce(cols[F],pp,p)[0]
    _,pf=field_basis((*B,F),cols.__getitem__,p)
    assert not field_reduce(cols[E],pf,p)[0]
    state=record['state'];values=state['values']
    assert all(dot(w,v)<tau for v in state['fixed'])
    cheap_proj=tuple(v for v in (*state['C'],*state['D']) if dot(w,v)<tau)
    PB,pv=field_basis(cheap_proj,values.__getitem__,p)
    assert field_reduce(values[E],pv,p)[0] and field_reduce(values[F],pv,p)[0]
    _,pF=field_basis((*PB,F),values.__getitem__,p)
    assert not field_reduce(values[E],pF,p)[0]
    COUNTS['entire_state_original_pair_rank_equalities']+=1
    if strict:
        if not E[0]:
            _,representation=field_reduce(cols[E],pf,p)
            H=Counter({(E[1]-t,E[2]-t):1})
            for v,c in representation.items():
                if not v[0]:H[v[1]-t,v[2]-t]-=c
            HY={v:c for v,c in clean(H,p).items() if v[0]<a-t and v[1]<t}
            NY={(data[0][r][0]-t,data[0][r][1]):c for r,c in cols[0,t,t].items()
                if data[0][r][0]>=t and data[0][r]!=(t,0)}
            if mul(NY,HY,a-t,t,p):COUNTS['nonscalar_original_primal_relations']+=1
            allupper=tuple(v for v in cols if v[0])
            _,pe=field_basis((*B,E,*allupper),cols.__getitem__,p)
            ell={data[0].index(v[1:]):field_reduce(cols[v],pe,p)[1].get(E,0) for v in allupper}
            KY=reflected(ell,data,a-t,t)[0]
            if mul(NY,KY,a-2*t,t,p):COUNTS['nonscalar_original_dual_separators']+=1
        l,r=record['lo'],record['hi'];delta=(r-l)/10
        u0=(l+r)/2
        def point(uu):
            la=(E[1]-F[1])*uu+E[2]-F[2]
            return (-la if E[0] else la,uu,Q(1))
        w=wall_weight(point(u0-delta),point(u0+delta),cols,q);lam,u,_=w
        costs=sorted({dot(w,v) for v in cols})
        small=min(lam/2,min(y-x for x,y in zip(costs,costs[1:]))/4)
        adj=[]
        for sign in (-1,1):
            wp=(lam+sign*small,u,Q(1))
            S,_=field_basis(sorted(cols,key=lambda v:(dot(wp,v),v)),cols.__getitem__,p)
            assert all(dot(wp,v)>0 for v in normals(S,cols,p));adj.append(set(S))
            COUNTS['original_adjacent_strict_minima']+=1
        if E[0]:
            assert E in adj[0] and F not in adj[0] and F in adj[1] and E not in adj[1]
        else:
            assert F in adj[0] and E not in adj[0] and E in adj[1] and F not in adj[1]
        COUNTS['unrestricted_original_wall_certificates']+=1
    return {'heavy_label':E,'partner':F,'kind':kind,'state_prefix':state['prefix'],
            'weight':list(map(str,w))}


def run_case(t,a,units,p,E,kinds,expected):
    assert E[1]>=(3+E[0])*t
    data,co,cols=assemble(t,a,units,p);records=[];total=Counter();found=set()
    for J in own_sectors(t):
        states=flag_states(data,cols,t,a,E,J,p)
        for kind in kinds:
            for state in states:
                active=state_frontier(state,t,a,E,J,kind,p)
                for i,record in enumerate(active):
                    strict=(kind,record['partner']) not in found
                    cert=verify_original(data,cols,t,a,E,record,J,kind,p,strict)
                    records.append(cert);found.add((kind,record['partner']))
                    total[kind]+=1
    for F,w,kind in expected:
        assert (kind,F) in found
    bound={k:8*t**6 if k=='lower_low' else Q(32*a*t**6,E[1]) for k in kinds}
    for k,c in total.items():assert c<=bound[k]
    print('barrier case passed',t,a,E,p,'states/pieces',sum(total.values()),flush=True)
    return {'t':t,'a':a,'field_characteristic':p,'heavy_label':E,
            'coefficients':list(map(int,co)),'candidate_piece_counts':dict(total),
            'original_rank_certificates':records}


def cases():
    out=[]
    source=json.loads((ROOT/'results/research/hook_transverse_prefix_checks.json').read_text())
    for t,a,E,F,w in [(3,30,(0,16,3),(1,1,18),(Q(15,64),Q(65,64),Q(1))),
                      (4,32,(0,17,5),(1,3,18),(Q(107,65),Q(68,65),Q(1)))]:
        old=next(c for c in source['actual_wall_controls'] if c['t']==t and c['a']==a)
        # Recover the exact arbitrary higher coefficients from the saved source.
        from run_hook_bridge_checks import model
        data=model(t,a);f={(i,j):c for (g,i,j),c in zip(data[2],old['coefficients']) if g==1 and c}
        U={(i-t,j):c for (i,j),c in f.items() if i>=t}
        V={(j-t,i):c for (i,j),c in f.items() if j>=t}
        out.append(run_case(t,a,(U,V),101,E,['lower_low','lower_middle'],[(F,w,'lower_low')]))
    for p in (0,2):
        t=3;A=12;a=15
        U={(0,0):1,(1,2):1,(10,1):1};V={(0,0):1,(10,1):1}
        out.append(run_case(t,a,(U,V),p,(0,12,3),['lower_low','lower_middle'],[]))
    t=4;a=16;A=12;units=[]
    for corners in [((1,3,1),(4,2,2),(8,1,1)),((2,3,1),(5,2,3),(10,1,1))]:
        N=Counter()
        for nu,d,c in corners:
            for i in range(A-nu):N[nu+i,d]+=c*comb(A-nu-1,i)
        units.append({(0,0):1,**clean(N,101)})
    out.append(run_case(t,a,units,101,(0,12,5),['lower_low','lower_middle'],[]))
    m=6;A=80;a=82;t=2
    U={(0,0):1,(m,1):1};V={(0,0):1,**{(i,1):comb(A-1,i) for i in range(A)}}
    out.append(run_case(t,a,(U,V),1009,(0,8,3),['lower_low','lower_middle'],[]))
    out.append(run_case(t,a,(U,V),1009,(1,8,1),['upper_wedge'],
                       [((0,2,39),(Q(4,5),Q(31,5),Q(1)),'upper_wedge')]))
    assert COUNTS['unrestricted_original_wall_certificates'] and COUNTS['independent_omitted_prefix_kernel_equalities']
    assert COUNTS['nonscalar_original_primal_relations'] and COUNTS['nonscalar_original_dual_separators']
    return out




def lower_hull_changes(lines):
    best={}
    for a,b in lines:best[a]=min(best.get(a,b),b)
    hull=[]
    for a,b in sorted(best.items(),reverse=True):
        start=None
        while hull:
            aa,bb,old=hull[-1];start=Q(b-bb,aa-a)
            if old is None or start>old:break
            hull.pop()
        hull.append((a,b,start if hull else None))
    return {s for a,b,s in hull if s is not None and s>1}


def complete_sector_control():
    t=2;a=12;A=a-t;p=101
    U={(0,0):1,**{(i,1):comb(A-1,i) for i in range(A)}}
    V={(0,0):1,**{(i,1):2*comb(A-1,i) for i in range(A)}}
    data,co,cols=assemble(t,a,(U,V),p)
    n=2*len(data[0]);queue=[((1,2*n+1,(2*n+1)**2),None)];cells={}
    triangle=[tuple(Q(i==j) for i in range(3)) for j in range(3)]
    wedge=((0,1,-1),(-1,t,0))
    def restricted(poly):
        for q in wedge:poly=clip(poly,q,1)
        return poly
    while queue:
        w,old=queue.pop()
        B=frozenset(field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)[0])
        if old is not None:assert B!=old
        if B in cells:continue
        NN=normals(B,cols,p);assert all(dot(w,q)>0 for q in NN)
        poly=polygon(NN);cells[B]=(poly,NN)
        for left,right in zip(poly,poly[1:]+poly[:1]):
            if any(left[i]==right[i]==0 for i in range(3)):continue
            q=next(q for q in NN if dot(q,left)==dot(q,right)==0)
            nw,_=exact_adjacent_weight(n,left,right,q);queue.append((nw,B))
    assert sum(area2(v[0]) for v in cells.values())==1
    sector_cells=[];edges={}
    boundaries=((1,0,0),(0,1,-1),(-1,t,0),(0,0,1))
    for B,(poly,NN) in cells.items():
        piece=restricted(poly)
        if not piece or not area2(piece):continue
        sector_cells.append(piece)
        for left,right in zip(piece,piece[1:]+piece[:1]):
            if any(dot(q,left)==dot(q,right)==0 for q in boundaries):continue
            q=next(q for q in NN if dot(q,left)==dot(q,right)==0)
            edges[tuple(sorted((left,right)))]=q
    assert sum(map(area2,sector_cells))==area2(restricted(triangle))
    nonvertical=sum(bool(q[0]) for q in edges.values())
    boundary_points=set()
    for (left,right),q in edges.items():
        if q[0]:continue
        for v in (left,right):
            if v[0]==0:boundary_points.add(('zero',v))
            if -v[0]+t*v[1]==0:boundary_points.add(('top',v))
    assert len(sector_cells)<=1+nonvertical+len(boundary_points)
    profiles=[tuple(sum(v[i] for v in B) for i in range(3)) for B in cells]
    changes={'zero':lower_hull_changes([(y,z) for x,y,z in profiles]),
             'top':lower_hull_changes([(t*x+y,z) for x,y,z in profiles])}
    for name,v in boundary_points:
        assert v[2] and v[1]/v[2] in changes[name]
        COUNTS['independent_vertical_boundary_shadow_matches']+=1
    COUNTS['complete_original_positive_cones']+=len(cells)
    COUNTS['complete_original_low_wedge_cells']+=len(sector_cells)
    COUNTS['complete_sector_nonvertical_edge_pieces']+=nonvertical
    COUNTS['complete_sector_sweep_equalities']+=1
    return {'t':t,'a':a,'field_characteristic':p,'coefficients':list(map(int,co)),
            'complete_positive_cones':len(cells),'low_wedge_cells':len(sector_cells),
            'nonvertical_edge_pieces':nonvertical,'vertical_boundary_points':len(boundary_points),
            'boundary_changes':{k:sorted(map(str,v)) for k,v in changes.items()},
            'area_coverage_verified':True}



def vertical_boundary_control():
    # An independent rank-one represented-matroid fan exercises the
    # boundary term with no nonvertical edges. It is not a hook quotient.
    t=2;labels=((0,4,0),(0,0,8));columns={v:{0:1} for v in labels}
    assert len(field_basis(labels,columns.__getitem__,101)[0])==1
    triangle=[tuple(Q(i==j) for i in range(3)) for j in range(3)]
    wedge=((0,1,-1),(-1,t,0));pieces=[]
    for q in ((0,1,-2),(0,-1,2)):
        poly=clip(triangle,q,1)
        for side in wedge:poly=clip(poly,side,1)
        pieces.append(poly)
    region=triangle
    for side in wedge:region=clip(region,side,1)
    assert sum(map(area2,pieces))==area2(region) and all(area2(v)>0 for v in pieces)
    points=((Q(0),Q(2,3),Q(1,3)),(Q(4,7),Q(2,7),Q(1,7)))
    for point in points:
        assert point[1]/point[2] in lower_hull_changes([(4,0),(0,8)])
        COUNTS['independent_vertical_boundary_shadow_matches']+=1
    assert len(pieces)<=1+len(points)
    COUNTS['nontrivial_vertical_only_sweep_controls']+=1
    return {'kind':'abstract rank-one represented matroid','labels':labels,
            'sector_cells':2,'nonvertical_pieces':0,'boundary_points':2}


def full_family_comparison():
    # Validate the original reports rather than silently refreshing their
    # source hashes. The uniform bound is the analytical corollary BS0.
    inputs={}
    for name in ('hook_bridge_checks','two_colour_braid_checks'):
        path='results/research/'+name+'.json'
        report=json.loads((ROOT/path).read_text())
        assert report.get('status')=='passed' or report.get('all_checks_passed') is True
        hashes=report['source_sha256']
        assert hashes and all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==sha
                              for f,sha in hashes.items())
        inputs[path]={'source_hashes_verified':len(hashes),
                      'sha256':hashlib.sha256((ROOT/path).read_bytes()).hexdigest()}
        COUNTS['verified_full_family_comparison_inputs']+=1
    return {'full_family_bound':'O(n^2*t^6)',
            'scope':'Every compatible flat hook map, including arbitrary outer images, over every field.',
            'derivation':'O(n^2) convex bridge-order chambers; fixed rank k=3t^2 dual on at most 3k original binary-layer labels; two-colour bound O(k^3) per chamber.',
            'sector_count_is_asymptotic_improvement':False,
            'other_full_family_ceilings':['O(n^(7/3))','O(n^(3/2)*a)'],
            'inputs':inputs}

def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    comparison=full_family_comparison()
    models=cases()
    sector=complete_sector_control()
    vertical=vertical_boundary_control()
    files={'scripts/run_hook_barrier_sector_checks.py','docs/hook_barrier_sectors.md',
           'docs/hook_weighted_partners.md','docs/hook_bridge_reduction.md',
           'docs/two_colour_threshold_braids.md',
           'docs/cross_layer_direction_budget.md','docs/planar_dyadic_bound.md',
           'results/research/hook_bridge_checks.json','results/research/two_colour_braid_checks.json',
           'results/research/hook_transverse_prefix_checks.json'}
    pending=['scripts/run_hook_barrier_sector_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=[node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else []
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:
                    files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_barrier_sector_checks.py',
            'scope':'Exact critical bridge flags and weighted cross-layer partner bounds in w_x<t*max(w_y,w_z), for arbitrary mixed constant-axis units and zero outer images. The sector count is dominated by the existing full-hook O(n^2*t^6) corollary; no unrestricted G_3 improvement.',
            'counts':dict(sorted(COUNTS.items())),'original_quotient_replay_counts':dict(REPLAY),'cases':models,
            'full_family_comparison':comparison,
            'complete_sector_control':sector,
            'independent_vertical_boundary_control':vertical,
            'limits':['The exact weighted barrier model assumes both outer images zero and the stated price wedge.',
                      'A weighted high-price or arbitrary-outer-image barrier theorem is not proved; those cases already have the coarse full-hook O(n^2*t^6) count.',
                      'The restricted fan count does not improve the saved full-fan ceilings.',
                      'Thickness powers are not asserted optimal.','Exact priority and the original unrestricted endpoint remain open.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_barrier_sector_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}))


if __name__=='__main__':main()
