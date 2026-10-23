#!/usr/bin/env python3
"""Exact paired-chain minima, partial-wall arrangements, and outer-map guards."""
from collections import Counter,defaultdict
from fractions import Fraction as Q
from itertools import product,combinations
from pathlib import Path
from random import Random
import ast,hashlib,json

from run_hook_bridge_checks import model,actual,normals,polygon,strict_weight,COUNTS as REPLAY
from run_hook_rank_departure_checks import complete_fan,plane,rank,dual,COUNTS as FAN
from run_hook_coherent_unit_window_checks import lower,upper
from run_hook_monomial_orbit_checks import ray,det
from run_signed_face_checks import field_basis,field_reduce
from run_structural_complexity_checks import dot,clip,area2

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()
TRIANGLE=[tuple(Q(i==j) for i in range(3)) for j in range(3)]


def build(t,a,nus,p,coefficients=((1,1),(1,1)),outer=False):
    assert min(nus)>=t and max(nus)<a-t
    data=model(t,a);rng=Random(2026102336+t+a+p);co=[]
    for g,y,z in data[2]:
        arm=0 if y>=t else 1;q,j=(y-t,z) if arm==0 else(z-t,y)
        alpha,beta=coefficients[arm]
        c=(alpha if (q,j)==(0,0) else beta if (q,j)==(nus[arm],1) else 0) if g==1 else rng.randrange(p or 5) if outer else 0
        co.append(c)
    assert all(c%p if p else c for pair in coefficients for c in pair)
    return data,co,actual(data,co,p,True)


def root(q,j,nu):
    k=min(q//nu,j)
    return q-k*nu,j-k


def chains(t,a,nus):
    out=[]
    for arm,nu in enumerate(nus):
        rows=defaultdict(list)
        for q,j in product(range(a-t),range(t)):rows[root(q,j,nu)].append((q,j))
        out.append({r:tuple(sorted(v)) for r,v in rows.items()})
    return out


def partition_controls(data,cols,t,a,nus,cc):
    row_to_part={}
    for arm,C in enumerate(cc):
        for r,nodes in C.items():
            key=('pair',r[0],r[1]) if r[0]<t and arm==0 else ('pair',r[1],r[0]) if r[0]<t else ('single',arm,*r)
            for q,j in nodes:
                xy=(t+q,j) if arm==0 else(j,t+q)
                row_to_part[data[0].index(xy)]=key
    blocks=defaultdict(list)
    for v,col in cols.items():
        if not col:continue
        parts={row_to_part[r] for r in col};assert len(parts)==1
        blocks[next(iter(parts))].append(v)
        COUNTS['independent_original_support_direct_summand_checks']+=1
    assert sum(k[0]=='pair' for k in blocks)==t*t
    for q,j in product(range(t),repeat=2):
        k=('pair',q,j);H=(0,t+q,t+j)
        assert H in blocks[k] and cc[0][q,j][0]==(q,j) and cc[1][j,q][0]==(j,q)
        assert len(blocks[k])==2*len(cc[0][q,j])+2*len(cc[1][j,q])-1
    return blocks


def segment_basis(nodes,arm,c,s,t):
    m=len(nodes)
    if not m:return []
    if c < -s:return [upper(t,arm,*v) for v in nodes]
    if c < (m-1)*s:return [upper(t,arm,*nodes[0])]+[lower(t,arm,*v) for v in nodes[:-1]]
    return [lower(t,arm,*v) for v in nodes]


def critical(nodes,arm,c,s,t):
    m=len(nodes)
    if m==1:return upper(t,arm,*nodes[0])
    if c < (m-2)*s:return upper(t,arm,*nodes[1])
    return upper(t,arm,*nodes[0]) if c>(m-1)*s else lower(t,arm,*nodes[-1])


def formula_basis(t,cc,nus,w):
    out=[];cs=(w[0]-t*w[2],w[0]-t*w[1]);ss=(nus[0]*w[1]+w[2],nus[1]*w[2]+w[1])
    for arm,C in enumerate(cc):
        for r,nodes in C.items():
            if r[0]>=t:out.extend(segment_basis(nodes,arm,cs[arm],ss[arm],t))
    for q,j in product(range(t),repeat=2):
        nodes=(cc[0][q,j],cc[1][j,q]);B=[]
        for arm in (0,1):
            B.append(upper(t,arm,*nodes[arm][0]))
            B.extend(segment_basis(nodes[arm][1:],arm,cs[arm],ss[arm],t))
        e=max((critical(nodes[arm],arm,cs[arm],ss[arm],t) for arm in (0,1)),key=lambda v:dot(w,v))
        H=(0,t+q,t+j)
        if dot(w,H)<dot(w,e):B.remove(e);B.append(H)
        out.extend(B)
    assert len(out)==len(set(out))
    return frozenset(out)


def Y(t,nu,k):return (1,-k*nu,-t-k)
def Z(t,nu,k):return (1,-t-k,-k*nu)


def inventory(t,a,nus):
    ny,nz=nus;hy,hz=(min(t,1+(a-t-1)//v) for v in nus)
    full={plane(Y(t,ny,k)) for k in range(-1,hy)}|{plane(Z(t,nz,k)) for k in range(-1,hz)}
    for eps,eta in product((0,1),repeat=2):full.add(plane((0,t+eps*ny-eta,eps-t-eta*nz)))
    for j,eps in product(range(1,hz),(0,1)):full.add(plane((1,eps*ny-j,eps-t-j*nz)))
    for i,eta in product(range(1,hy),(0,1)):full.add(plane((1,eta-t-i*ny,eta*nz-i)))
    seg=[]
    for i,j in product(range(1,hy),range(1,hz)):
        u=Q(j*nz-i,i*ny-j);sy=ny*u+1;sz=nz+u
        lo=max(t+(i-1)*sy,t*u+(j-1)*sz);hi=min(t+i*sy,t*u+j*sz)
        if lo>=hi:continue
        restrictions=(Y(t,ny,i-1),tuple(-v for v in Y(t,ny,i)),Z(t,nz,j-1),tuple(-v for v in Z(t,nz,j)))
        q=plane((0,i*ny-j,i-j*nz))
        cuts={-(b*u+c)/v for v,b,c in full if v and lo<-(b*u+c)/v<hi}
        assert len(cuts)<=12
        COUNTS['independent_partial_wall_constant_crossing_bounds']+=1
        seg.append((q,restrictions,u,lo,hi,len(cuts)))
    return full,seg,(hy,hz)


def arrangement(t,a,nus):
    full,segments,hs=inventory(t,a,nus);cells=[TRIANGLE]
    for q in sorted(full):
        nxt=[]
        for cell in cells:
            for sign in (-1,1):
                part=clip(cell,q,sign)
                if part and area2(part):nxt.append(part)
        cells=nxt
    base_count=len(cells)
    for q,restrictions,u,lo,hi,cuts in segments:
        nxt=[];splits=0
        for cell in cells:
            pos=clip(cell,q,1);neg=clip(cell,q,-1)
            if pos and neg and area2(pos) and area2(neg):
                chord=tuple(v for v in pos if dot(q,v)==0)
                assert len(set(chord))==2
                mid=tuple(sum(v[j] for v in chord)/len(chord) for j in range(3))
                if all(dot(r,mid)>=0 for r in restrictions):
                    assert all(dot(r,v)>=0 for r in restrictions for v in chord)
                    nxt.extend((pos,neg));splits+=1;continue
            nxt.append(cell)
        # Previously added vertical segments are parallel or collinear.
        assert splits<=cuts+1
        cells=nxt
    assert sum(map(area2,cells))==1
    L=len(full);assert len(cells)<=1+L*(L+1)//2+13*(hs[0]-1)*(hs[1]-1)
    weights=[]
    for cell in cells:
        rays=tuple(map(ray,cell));triple=next(v for v in combinations(rays,3) if det(*v))
        w=tuple(sum(v[j] for v in triple) for j in range(3))
        assert min(w)>0 and max(w)<=100*a*a
        assert all(dot(w,q)!=0 for q in full)
        weights.append(w)
    COUNTS['independent_complete_partial_wall_arrangements']+=1
    return cells,weights,{'h':hs,'full_planes':L,'restricted_vertical_segments':len(segments),'base_cells':base_count,'integer_weights':len(weights),'max_height':max(map(max,weights))}


def local_formula_controls():
    rng=Random(2026102337)
    for p in (0,2,3,101):
        for m in range(1,10):
            for z in range(30):
                c=Q(rng.randrange(-40,41),7);s=Q(rng.randrange(1,40),11)
                if c in (-s,(m-1)*s):continue
                labels=[(kind,i) for kind in (0,1) for i in range(m)]
                cols={(1,i):{i:1} for i in range(m)}
                cols.update({(0,i):{i:1,**({i+1:1} if i+1<m else {})} for i in range(m)})
                cost=lambda v:v[1]*s+(c if v[0] else 0)
                B=set(field_basis(sorted(labels,key=cost),cols.__getitem__,p)[0])
                expected={(1,i) for i in range(m)} if c<-s else {(1,0)}|{(0,i) for i in range(m-1)} if c<(m-1)*s else {(0,i) for i in range(m)}
                assert B==expected
                COUNTS['independent_graphic_chain_minimum_bases']+=1


def outer_guard(t,a,nus,p,weights,baseline_cones):
    data,co,cols=build(t,a,nus,p,outer=True);fan=complete_fan(cols,4*a*t-2*t*t,p)
    seen={frozenset(field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)[0]) for w in weights}
    assert seen<set(fan)
    full,segments,hs=inventory(t,a,nus);all_planes=full|{z[0] for z in segments};witness=None
    for B,poly in fan.items():
        N=normals(B,cols,p)
        for x,y in zip(poly,poly[1:]+poly[:1]):
            if any(x[i]==y[i]==0 for i in range(3)):continue
            q=next(q for q in N if dot(q,x)==dot(q,y)==0)
            if not q[0] or plane(q) in all_planes:continue
            w=tuple((r+s)/2 for r,s in zip(x,y));piv=field_basis(tuple(B),cols.__getitem__,p)[1]
            for f in cols:
                if f in B:continue
                rem,rep=field_reduce(cols[f],piv,p);assert not rem
                for e,c in rep.items():
                    if not c or plane(tuple(r-s for r,s in zip(f,e)))!=plane(q) or dot(w,e)!=dot(w,f):continue
                    BB=(set(B)-{e})|{f};assert rank(tuple(BB),cols,p)==len(B)
                    assert all(dot(w,z)>=0 for z in N)
                    witness={'exchange':[e,f],'weight':list(map(str,w)),'wall_plane':plane(q),'outside_every_baseline_support_plane':True}
                    break
                if witness:break
            if witness:break
        if witness:break
    assert witness
    COUNTS['independent_actual_outer_coupling_exchange_guards']+=1
    return {'t':t,'a':a,'nus':nus,'characteristic':p,'baseline_cones':baseline_cones,'actual_outer_cones':len(fan),'baseline_weights':len(weights),'outer_bases_seen_by_baseline_weights':len(seen),'actual_exchange_witness':witness}


def lower_junction_controls(t):
    nu=2*t;a=t*nu+2*t;p=101
    data,co,cols=build(t,a,(nu,nu),p);cc=chains(t,a,(nu,nu));blocks=partition_controls(data,cols,t,a,(nu,nu),cc)
    vertices=set();certificates=0
    for k,l in product(range(1,t),repeat=2):
        u=Q(l*nu-t-k,k*nu-t-l);lam=t+k*(nu*u+1)
        assert min(u,lam)>0 and (lam,u) not in vertices;vertices.add((lam,u))
        keys=(('single',0,a-t-(k+1)*nu,0),('single',1,a-t-(l+1)*nu,0))
        labels=tuple(v for key in keys for v in blocks[key]);sub={v:cols[v] for v in labels}
        seen=set();eps=Q(1,100*a**4)
        for sy,sz in product((-1,1),repeat=2):
            du=(sz-sy)*eps/(k*nu-t-l);dl=sy*eps+k*nu*du
            w=(lam+dl,u+du,Q(1));assert min(w)>0
            B=frozenset(field_basis(sorted(sub,key=lambda v:(dot(w,v),v)),sub.__getitem__,p)[0])
            assert all(dot(w,q)>0 for q in normals(B,sub,p));seen.add(B)
        assert len(seen)==4
        COUNTS['independent_actual_direct_summand_four_region_junctions']+=1;certificates+=1
    return {'t':t,'a':a,'nu':nu,'positive_distinct_junctions':certificates,'four_strict_original_component_bases_per_junction':True}


def balanced_outer_controls(t,a,nus,p):
    data,co,cols=build(t,a,nus,p,outer=True);_,_,base=build(t,a,nus,p);cc=chains(t,a,nus)
    outer=tuple(v for v in cols if not v[0] and (v[1]>=a or v[2]>=a))
    H=tuple(v for v in outer if (v[1]-a if v[1]>=a else v[2]-a)<4*t)
    late=tuple(v for v in outer if v not in H);assert len(H)<=8*t*t
    for arm in (0,1):
        pure=tuple(v for v in data[4][arm] if not v[0] and v[arm+1]<a)
        piv=field_basis(pure,cols.__getitem__,p)[1]
        for e in late:
            if e not in data[4][arm]:continue
            assert not field_reduce(cols[e],piv,p)[0]
            for u in (Q(1,2),Q(1),Q(2)):
                w=(Q(1),u,Q(1))
                assert max(dot(w,v) for v in pure)<dot(w,e)
            COUNTS['independent_balanced_late_outer_prefix_loop_checks']+=1
    patches=[];balance={(0,2,-1),(0,-1,2)}
    for u,lam in ((Q(5001,10000),Q(1,31)),(Q(10001,10000),Q(t)+Q(1,37)),(Q(19999,10000),Q(2*a)+Q(1,41))):
        w=(lam,u,Q(1));B0=formula_basis(t,cc,nus,w)
        assert B0==frozenset(field_basis(sorted(base,key=lambda v:(dot(w,v),v)),base.__getitem__,p)[0])
        core=tuple(sorted(B0))+H;tiny={v:cols[v] for v in core}
        B=frozenset(field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)[0])
        BT=frozenset(field_basis(sorted(tiny,key=lambda v:(dot(w,v),v)),tiny.__getitem__,p)[0]);assert B==BT
        N0=normals(B0,base,p)|balance;N=normals(B,cols,p);T=normals(BT,tiny,p)
        original=polygon(N0|N);restricted=polygon(N0|T)
        assert all(dot(v,q)>=0 for v in original for q in T)
        assert all(dot(v,q)>=0 for v in restricted for q in N)
        assert area2(original)==area2(restricted)>0
        dd=dual(core,tiny,p);D=frozenset(field_basis(sorted(dd,key=lambda v:(-dot(w,v),v)),dd.__getitem__,p)[0])
        assert len(D)==len(H) and set(core)-set(D)==set(B)
        assert normals(D,dd,p,maximize=True)==T
        patches.append({'weight':list(map(str,w)),'original_rank':len(B0),'retained_outer_columns':len(H),
                        'dual_rank':len(D),'dual_labels':len(core),'whole_patch_area2':str(area2(original))})
        COUNTS['independent_entire_balanced_outer_extension_and_dual_patches']+=1
    return {'t':t,'a':a,'nus':nus,'characteristic':p,'late_original_outer_columns':len(late),'patches':patches}


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    previous=json.loads((ROOT/'results/research/hook_rational_leading_checks.json').read_text());assert previous['status']=='passed'
    for f,h in previous['source_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
    local_formula_controls();cases=[];guards=[];rng=Random(2026102338)
    configs=[(2,7,(2,2),101,((1,1),(1,1)),True),(3,11,(3,3),2,((1,1),(1,1)),True),
             (2,12,(4,4),3,((1,2),(2,1)),True),(3,18,(6,6),101,((2,3),(5,7)),True),
             (2,8,(2,3),0,((2,-3),(-5,7)),True),(3,12,(3,4),3,((1,2),(2,1)),True),
             (5,32,(5,7),2,((1,1),(1,1)),False)]
    for t,a,nus,p,coefficients,full in configs:
        data,co,cols=build(t,a,nus,p,coefficients);cc=chains(t,a,nus);blocks=partition_controls(data,cols,t,a,nus,cc)
        cells,weights,stats=arrangement(t,a,nus);seen=set()
        for w in weights:
            B=frozenset(field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)[0])
            assert B==formula_basis(t,cc,nus,w)
            assert all(dot(w,q)>0 for q in normals(B,cols,p));seen.add(B)
            COUNTS['independent_original_formula_and_strict_cone_checks']+=1
        fan=complete_fan(cols,4*a*t-2*t*t,p) if full else None
        if full:assert seen==set(fan)
        cases.append({'t':t,'a':a,'nus':nus,'characteristic':p,'coefficients':coefficients,'direct_summands':len(blocks),'arrangement':stats,'original_bases_seen':len(seen),'complete_original_cones':len(fan) if full else None})
        if (t,a,nus,p) in ((2,7,(2,2),101),(3,11,(3,3),2)):guards.append(outer_guard(t,a,nus,p,weights,len(fan)))
        print('Paired chains passed',t,a,nus,p,'cones',len(fan) if full else None,'weights',len(weights),flush=True)
    for t in range(2,20):
        for nus in ((t,t),(2*t,3*t+1)):
            inventory(t,(t+1)*max(nus)+t,nus)
    lower_checks=[lower_junction_controls(t) for t in (2,3,5)]
    balanced=[balanced_outer_controls(t,a,nus,p) for t,a,nus,p in ((2,30,(2,4),101),(3,16,(3,5),2),(2,10,(2,3),0))]
    files={'docs/hook_paired_chain_fans.md','scripts/run_hook_paired_chain_checks.py','docs/hook_bridge_reduction.md','docs/hook_monomial_orbits.md','results/research/hook_rational_leading_checks.json'}
    pending=['scripts/run_hook_paired_chain_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=([node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else [])
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_paired_chain_checks.py',
            'scope':'Zero outer maps give O((hY+hZ)^2) full fans and common weights, sharply quadratic in t. Arbitrary compatible outer maps on the balanced sector have exact O(t^2)-column extensions of this backbone and fixed dual-rank encodings, yielding O(a*t^3*(hY+hZ)^2) regions by the classical lattice-vertex bound. Original outer couplings are retained and separately checked. No unrestricted endpoint or exact novelty claimed.',
            'uniform_result':{'central_arms':'W_i=alpha_i+beta_i*q^nu_i*v, alpha_i*beta_i!=0, t<=nu_i<A','outer_maps':'both zero',
                              'chain_lengths':'h_i=min(t,ceil(A/nu_i))','full_fan':'O((hY+hZ)^2)','common_integer_weights':'O((hY+hZ)^2), height O(a^2), all fields and nonzero binomial coefficients',
                              'sharp_family':'nuY=nuZ=nu>=2t, a>=t*nu+2t: (t-1)^2 positive intersections of independent unmarked-chain walls, hence g=Theta(t^2)',
                              'balanced_arbitrary_outer_extension':'On 1/2<wy/wz<2, retain at most 8t^2 original outer columns. On each zero-outer cone all minima are contained in its fixed old basis plus those columns; the dual has rank at most 8t^2.',
                              'balanced_arbitrary_outer_count':'O(a*t^3*(hY+hZ)^2)=O(n*t^2*(hY+hZ)^2), at most O(n*t^4); retain earlier bounds when smaller. The lattice exponent is classical.',
                              'partial_wall_inventory':'O(hY+hZ) full lines and O(hY*hZ) restricted vertical segments; each segment meets at most twelve nonvertical full lines'},
            'counts':dict(sorted(COUNTS.items())),'fan_counts':dict(FAN),'original_quotient_replay_counts':dict(REPLAY),'cases':cases,'sharp_family_junction_checks':lower_checks,'balanced_arbitrary_outer_checks':balanced,'actual_outer_coupling_guards':guards,
            'limits':['The zero outer hypothesis is essential to the proved support partition; the finite guards do not prove an asymptotic outer lower bound.',
                      'Steps below t can connect several head bridges in one chain and are not covered.',
                      'Graphic chain minima, greedy single-column insertion and the matroid face dictionary are classical; exact priority of the combined class bound is unestablished.',
                      'The original unrestricted sharp G_3 or actual polynomial-excess objective and final-paper requirements remain unresolved.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_paired_chain_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}),flush=True)


if __name__=='__main__':main()
