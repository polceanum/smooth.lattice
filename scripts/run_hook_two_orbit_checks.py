#!/usr/bin/env python3
"""Original-field coset components and common planes for two monomial arms."""
from collections import Counter,defaultdict
from fractions import Fraction as Q
from itertools import product,combinations
from pathlib import Path
from random import Random
import ast,hashlib,json

from run_hook_bridge_checks import model,actual,normals,polygon,COUNTS as REPLAY
from run_hook_rank_departure_checks import complete_fan,plane,weight_controls,COUNTS as FAN
from run_hook_monomial_orbit_checks import orbit_key,common_list,potential_core,inventory as one_arm_inventory,COUNTS as ORBIT
from run_hook_coherent_unit_window_checks import multiply
from run_signed_face_checks import field_basis,field_reduce
from run_structural_complexity_checks import dot,area2

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()
OFFSETS=((0,1,1),(1,1,0),(1,0,1))


def build(t,a,steps,gamma,p,replay=True):
    data=model(t,a);co=[]
    for g,y,z in data[2]:
        if g!=1:co.append(0);continue
        arm=0 if y>=t else 1;q,j=(y-t,z) if arm==0 else(z-t,y)
        nu,de=steps[arm];k=j//de
        co.append(gamma[arm][k] if j%de==0 and q==k*nu and k<len(gamma[arm]) else 0)
    return data,co,actual(data,co,p,replay)


def classify(t,steps):
    (ny,dy),(nz,dz)=steps
    if dy==dz==1 and ny>=1 and nz>=1:
        return 'diagonal' if ny==nz==1 else 'unit_short'
    if ny>=2*dy and nz>=2*dz:return 'steep'
    return 'bounded_one_arm'


def partition(data,t,a,steps):
    (ny,dy),(nz,dz)=steps;D=ny*nz-dy*dz
    diagonal=steps==((1,1),(1,1))
    def coset(y,z):
        if diagonal:return (y-z,)
        # Parallel nonprimitive steps are covered by the bounded-slope
        # one-arm theorem; merging their marked core is harmless.
        if not D:return (0,)
        return ((nz*y-dz*z)%abs(D),(-dy*y+ny*z)%abs(D))
    marked=[{orbit_key(q,j,*s) for q,j in product(range(t),repeat=2)} for s in steps]
    def group(arm,q,j):
        orbit=orbit_key(q,j,*steps[arm]);point=(q,j) if arm==0 else(j,q)
        return ('core',*coset(*point)) if orbit in marked[arm] else('chain',arm,*orbit)
    rowparts={}
    for i,(y,z) in enumerate(data[0]):
        if y<t and z<t:continue
        arm=0 if y>=t else 1;q,j=(y-t,z) if arm==0 else(z-t,y)
        rowparts[i]=group(arm,q,j)
    parts=defaultdict(list);inputs={};kinds={}
    for v in data[3]:
        if v in data[5] or (not v[0] and (v[1]>=a or v[2]>=a)):continue
        if v[0]:
            arm=0 if v[1]>=t else 1;q,j=(v[1]-t,v[2]) if arm==0 else(v[2]-t,v[1])
            point=(q,j) if arm==0 else(j,q);key=group(arm,q,j);kind=arm+1
        else:
            point=(v[1]-t,v[2]-t);kind=0
            if point[0]<t and point[1]<t:key=('core',*coset(*point))
            else:
                arm=0 if point[0]>=t else 1;q,j=point if arm==0 else point[::-1]
                key=group(arm,q,j)
        offset=OFFSETS[kind]
        assert v==(offset[0],point[0]+t*offset[1],point[1]+t*offset[2])
        parts[key].append(v);inputs[v]=point;kinds[v]=kind
    COUNTS['coefficient_independent_two_orbit_partitions']+=1
    return parts,rowparts,inputs,kinds,marked,D


def comparison_inventory(data,t,a,steps,parts,inputs,marked,D):
    mode=classify(t,steps);(ny,dy),(nz,dz)=steps
    if mode=='bounded_one_arm':
        arm=0 if ny<=2*dy else 1
        labels=potential_core(data,t,arm,*steps[arm],marked[arm])
        planes=one_arm_inventory(data,labels,arm)
        assert len(planes)<=2000*a*t
        return planes,{'method':mode,'chosen_arm':arm,'planes':len(planes)}
    caps=(2*t,0) if mode=='diagonal' else (7*t,7*t) if mode=='unit_short' else (3*((t+dy-1)//dy),3*((t+dz-1)//dz))
    planes=set();largest=[0,0]
    for key,labels in parts.items():
        for v,f in combinations(labels,2):
            y,z=(x-b for x,b in zip(inputs[v],inputs[f]))
            if mode=='diagonal':
                assert y==z and abs(y)<=caps[0]
                m,n=y,0
            else:
                mm=nz*y-dz*z;nn=-dy*y+ny*z
                assert mm%D==nn%D==0
                m,n=mm//D,nn//D
                assert abs(m)<=caps[0] and abs(n)<=caps[1],(mode,steps,t,key,v,f,m,n,caps)
            largest=[max(largest[0],abs(m)),max(largest[1],abs(n))]
            planes.add(plane(tuple(x-b for x,b in zip(v,f))))
            COUNTS['independent_integer_coset_comparison_coordinates']+=1
    bound=9*(2*caps[0]+1)*(2*caps[1]+1)
    assert len(planes)<=bound
    assert all(max(map(abs,q))<2*a for q in planes)
    return planes,{'method':mode,'lattice_index':abs(D),'basis_coordinate_caps':caps,'largest_observed_coordinates':largest,'planes':len(planes),'coarse_plane_bound':bound}


def support_controls(data,cols,parts,rowparts,t):
    assignment={v:key for key,g in parts.items() for v in g};overlaps=Counter()
    for v,col in cols.items():
        if not col:continue
        assert {rowparts[i] for i in col}=={assignment[v]}
        COUNTS['independent_original_row_coset_direct_sum_checks']+=1
        if not v[0] and v[1]<2*t and v[2]<2*t:
            for i in col:overlaps[i]+=1
    return max(overlaps.values(),default=0)


def basis_controls(data,cols,parts,planes,w,p,decompose=False):
    B=frozenset(field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)[0]);N=normals(B,cols,p)
    assert all(plane(q) in planes for q in N)
    if decompose:
        BB=set();NN=set()
        for labels in parts.values():
            sub={v:cols[v] for v in labels}
            local=field_basis(sorted(sub,key=lambda v:(dot(w,v),v)),sub.__getitem__,p)[0]
            BB.update(local);NN.update(normals(local,sub,p))
        assert BB==set(B) and NN==N
        assert area2(polygon(NN))==area2(polygon(N))>0
        COUNTS['independent_entire_component_vs_original_normal_cones']+=1
    COUNTS['independent_original_basis_plane_inclusions']+=1
    return B


def lattice_controls():
    rng=Random(2026102339)
    for t in range(2,35):
        for z in range(12):
            dy,dz=(rng.randrange(1,t) for _ in range(2))
            ny=2*dy+rng.randrange(20*t);nz=2*dz+rng.randrange(20*t)
            BY=t-1+ny*((t-1)//dy);BZ=t-1+nz*((t-1)//dz);D=ny*nz-dy*dz
            assert D>0
            assert Q(nz*BY+dz*BZ,D)<=Q(3*t,dy)
            assert Q(dy*BY+ny*BZ,D)<=Q(3*t,dz)
            COUNTS['independent_steep_lattice_rectangle_bounds']+=1
    for t in range(2,40):
        for ny,nz in ((1,2),(2,1),(t-1,t+1),(3*t,2*t)):
            BY=t-1+ny*(t-1);BZ=t-1+nz*(t-1);D=ny*nz-1
            if not D:continue
            assert Q(nz*BY+BZ,D)<=7*t and Q(BY+ny*BZ,D)<=7*t
            COUNTS['independent_unit_short_lattice_rectangle_bounds']+=1


def nearby_step_guard(t,a,nu,p):
    steps=((nu,1),(nu,1));data,co,base=build(t,a,steps,([1,1],[1,1]),p)
    parts,rp,inputs,kinds,marked,D=partition(data,t,a,steps)
    planes,stats=comparison_inventory(data,t,a,steps,parts,inputs,marked,D)
    co=[c+int(g==1 and ((y==t+nu+1 and z==1) or(y==1 and z==t+nu+1))) for c,(g,y,z) in zip(co,data[2])]
    cols=actual(data,co,p,True);assignment={v:key for key,group in parts.items() for v in group}
    mixed=sum(bool(col) and len({rp[i] for i in col})>1 for col in cols.values());assert mixed
    fan=complete_fan(cols,4*a*t-2*t*t,p);witness=None;extra=set()
    for B,poly in fan.items():
        N=normals(B,cols,p)
        for x,y in zip(poly,poly[1:]+poly[:1]):
            if any(x[i]==y[i]==0 for i in range(3)):continue
            q=next(q for q in N if dot(q,x)==dot(q,y)==0)
            if not q[0] or plane(q) in planes:continue
            extra.add(plane(q))
            if witness:continue
            w=tuple((s+r)/2 for s,r in zip(x,y));piv=field_basis(tuple(B),cols.__getitem__,p)[1]
            for f in cols:
                if f in B:continue
                rem,rep=field_reduce(cols[f],piv,p);assert not rem
                for e,c in rep.items():
                    if not c or plane(tuple(s-r for s,r in zip(f,e)))!=plane(q) or dot(w,e)!=dot(w,f):continue
                    BB=(set(B)-{e})|{f}
                    assert len(field_basis(tuple(BB),cols.__getitem__,p)[0])==len(B)
                    assert all(dot(w,z)>=0 for z in N)
                    assert assignment[e]!=assignment[f]
                    witness={'exchange':[e,f],'weight':list(map(str,w)),'wall_plane':plane(q),
                             'old_component_keys':[assignment[e],assignment[f]],'exchanges_between_distinct_old_components':True}
                    break
                if witness:break
    assert witness and extra
    COUNTS['independent_actual_nearby_step_cross_component_faces']+=1
    return {'t':t,'a':a,'nu':nu,'characteristic':p,'central_operator':'1+q^nu*v+q^(nu+1)*v',
            'original_columns_crossing_old_parts':mixed,'complete_original_cones':len(fan),
            'distinct_crosslayer_planes_outside_two_orbit_inventory':len(extra),'actual_exchange_witness':witness}


def remainder_face_controls():
    out=[];p=1009
    for t,a,modes,e,f,w in (
            (4,25,('plus','plus'),(1,11,1),(0,7,9),(Q(35,48),Q(23,192),Q(29,192))),
            (6,49,('power','geometric'),(1,13,5),(0,7,17),(Q(9,16),Q(25,96),Q(17,96)))):
        nu=t;A=a-t;data=model(t,a);M={(nu,1):1,(nu+1,1):1};one={(0,0):1};operators=[]
        for mode in modes:
            if mode=='plus':W={**one,**M}
            elif mode=='power':
                W=one
                for j in range(t-1):W=multiply(W,{**one,**M},A,t,p)
            else:
                W={};term=one
                for j in range(t):
                    for v,c in term.items():W[v]=(W.get(v,0)+c)%p
                    term=multiply(term,M,A,t,p)
            operators.append(W)
        co=[]
        for g,y,z in data[2]:
            arm=0 if y>=t else 1;q,j=(y-t,z) if arm==0 else(z-t,y)
            co.append(operators[arm].get((q,j),0) if g==1 else 0)
        cols=actual(data,co,p,True);assert e in data[4][0] and f in data[4][1]
        tau=dot(w,e);assert tau==dot(w,f) and min(w)>0
        cheap=tuple(v for v in cols if dot(w,v)<tau);piv=field_basis(cheap,cols.__getitem__,p)[1]
        er=field_reduce(cols[e],piv,p)[0];fr=field_reduce(cols[f],piv,p)[0]
        assert er and fr and len(field_basis((0,1),{0:er,1:fr}.__getitem__,p)[0])==1
        q=plane(tuple(x-y for x,y in zip(e,f)));assert q not in {plane(v) for v in data[7]}
        COUNTS['independent_entire_prefix_pure_cross_arm_remainder_faces']+=1
        out.append({'t':t,'a':a,'nu':nu,'characteristic':p,'M':'q^nu*(1+q)*v','operator_modes':modes,
                    'exchange':[e,f],'weight':list(map(str,w)),'primitive_plane':q,
                    'strictly_cheaper_rank':len(piv),'each_single_and_pair_increment':1,
                    'outside_original_known_direction_inventory':True})
    return out


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    previous=json.loads((ROOT/'results/research/hook_paired_chain_checks.json').read_text());assert previous['status']=='passed'
    for f,h in previous['source_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
    lattice_controls();rng=Random(2026102340);cases=[];saved=None
    configs=[(3,9,((2,1),(2,1)),([1,1,1],[1,2,3]),101,True,True),
             (4,13,((2,1),(3,1)),([1,1,1,1],[0,1,1,1]),2,True,True),
             (4,13,((4,2),(5,2)),([0,1],[1,2]),3,True,True),
             (3,10,((1,1),(1,1)),([0,1,1],[1,-1,2]),0,True,True),
             (3,10,((1,1),(2,1)),([1,1,1],[1,1,1]),3,True,True),
             (4,12,((2,2),(2,2)),([1,1],[0,1]),3,True,True),
             (3,9,((0,1),(4,2)),([1,2,1],[0,1]),101,True,True),
             (5,14,((3,2),(8,3)),([1,2,3],[0,1]),101,True,True),
             (6,24,((3,1),(4,1)),([1,1,0,1,1,1],[0,1,1,1,0,1]),2,False,True),
             (8,28,((8,3),(9,4)),([1,2,1],[0,1]),3,False,True),
             (16,40,((4,1),(5,1)),([1,1]*8,[0,1]*8),2,False,False)]
    for t,a,steps,gamma,p,full,replay in configs:
        data,co,cols=build(t,a,steps,gamma,p,replay)
        parts,rp,inputs,kinds,marked,D=partition(data,t,a,steps)
        planes,stats=comparison_inventory(data,t,a,steps,parts,inputs,marked,D)
        overlap=support_controls(data,cols,parts,rp,t)
        weights=weight_controls(cols,rng)
        for i,w in enumerate(weights[:6]):basis_controls(data,cols,parts,planes,w,p,i<2)
        fan=complete_fan(cols,4*a*t-2*t*t,p) if full else None
        if fan:
            for B in fan:
                assert all(plane(q) in planes for q in normals(B,cols,p))
                COUNTS['complete_original_fan_cones_in_common_inventory']+=1
        cases.append({'t':t,'a':a,'steps':steps,'polynomial_coefficients':gamma,'characteristic':p,'full_original_quotient_replay':replay,'partition_groups':len(parts),'max_original_head_columns_sharing_one_row':overlap,'inventory':stats,'complete_original_cones':len(fan) if fan else None})
        if (t,a,steps)==(3,9,((2,1),(2,1))):saved=(data,parts,rp,planes)
        print('Two orbit case passed',t,a,steps,p,'method',stats['method'],'planes',len(planes),'cones',len(fan) if fan else None,flush=True)
    data,parts,rp,planes=saved;t=3;a=9;steps=((2,1),(2,1));weights=common_list(planes,a);samplers=[]
    for p,gamma in ((2,([1,1,1],[1,1,0])),(3,([0,1,2],[1,2,1])),(0,([2,-1,3],[0,1,-2])),(101,([0],[0]))):
        data,co,cols=build(t,a,steps,gamma,p);support_controls(data,cols,parts,rp,t)
        fan=complete_fan(cols,4*a*t-2*t*t,p);seen=set()
        for w in weights:
            B=basis_controls(data,cols,parts,planes,w,p)
            assert all(dot(w,q)>0 for q in normals(B,cols,p));seen.add(B)
        assert seen==set(fan)
        COUNTS['independent_all_coefficient_common_list_equalities']+=1
        samplers.append({'characteristic':p,'polynomial_coefficients':gamma,'common_weights':len(weights),'complete_original_cones':len(fan)})
    guards=[nearby_step_guard(t,a,nu,p) for t,a,nu,p in ((3,9,2,101),(4,13,3,3))]
    remainders=remainder_face_controls()
    files={'docs/hook_two_orbit_fans.md','scripts/run_hook_two_orbit_checks.py','docs/hook_bridge_reduction.md','docs/hook_monomial_orbits.md','docs/profile_lattice_bound.md','results/research/hook_paired_chain_checks.json'}
    pending=['scripts/run_hook_two_orbit_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=([node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else [])
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_two_orbit_checks.py',
            'scope':'Both central arms arbitrary polynomials in one monomial each, with positive short steps, and zero outer images: uniformly O(n^2) full fans and a common quadratic integer list, all fields and nonunit or zero polynomials. A two-dimensional row-coset partition handles multiple bridges on one chain. No unrestricted endpoint or exact novelty claimed.',
            'uniform_result':{'central_arms':'W_i=F_i(q^nu_i*v^delta_i), nu_i>=0,1<=delta_i<t, F_i arbitrary','outer_images':'both zero',
                              'all_step_full_fan':'O(n^2), absolute constant','all_step_common_integer_list':'O(n^2), height O(a^2), fixed steps and arbitrary coefficients/fields',
                              'both_steep':'nu_i>=2*delta_i: O(t^4/(delta_Y^2*delta_Z^2)); comparison planes O(ceil(t/delta_Y)*ceil(t/delta_Z))',
                              'unit_short':'delta_Y=delta_Z=1,nu_Y,nu_Z>=1: O(t^4); if both nu_i=1, O(t^2)',
                              'bounded_slope_branch':'If one nu_i/delta_i<2, the earlier one-arm monomial-orbit theorem has beta<=3 and supplies the quadratic count/list',
                              'finite_offset_lattice_lemma':'A represented row direct sum with label offsets in s fixed classes and input differences in a lattice coordinate rectangle [-M,M]x[-N,N] has at most s^2*(2M+1)*(2N+1) comparison planes; its complete arrangement fixes every minimum basis. This is not an assertion of a variable-induced algebra grading.'},
            'counts':dict(sorted(COUNTS.items())),'fan_counts':dict(FAN),'original_quotient_replay_counts':dict(REPLAY),'one_arm_reused_counts':dict(ORBIT),
            'cases':cases,'common_sampler':samplers,'nearby_step_actual_face_guards':guards,'pure_cross_arm_remainder_faces':remainders,'common_arrangement':{'planes':len(planes),'integer_weights':len(weights)},
            'limits':['Both outer images vanish; arbitrary compatible outer maps can join the cosets and are not covered.',
                      'Central supports must lie on one monomial ray per arm. Two nearby long steps can destroy the lattice index and the unmarked-chain split.',
                      'The common inventory works over all fields; it does not claim the actual fan is independent of coefficients or characteristic.',
                      'Nearby-step guard fixtures have actual coherent cross-component faces outside every comparison plane of the theorem, but establish no asymptotic lower growth.',
                      'Two independently reconstructed full-prefix certificates also give actual pure cross-arm faces outside the original known-direction inventory; they establish no rate of growth.',
                      'The largest fixture checks the original residue matrix without replaying the full commuting-multiplication audit; the other fixtures replay actual quotients.',
                      'Lattice cosets, matroid direct sums and arrangements are classical; exact strongest-current novelty, sharp unrestricted growth, actual polynomial excess and the original final-paper requirements remain unresolved.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_two_orbit_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}),flush=True)


if __name__=='__main__':main()
