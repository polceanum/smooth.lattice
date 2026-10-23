#!/usr/bin/env python3
"""Exact transverse head filters and the common polynomial-module rank obstruction."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter,defaultdict
from itertools import product
from math import comb
import hashlib,json

from run_signed_face_checks import field_basis
from run_hook_axis_birth_checks import assignment
from run_hook_scalar_outer_lattice_checks import reduce_rows
from run_hook_two_orbit_checks import build
from run_hook_rank_departure_checks import complete_fan
from run_hook_bridge_checks import polygon
from run_hook_axis_compilation_checks import fixed_profile,hull
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_structural_complexity_checks import dot,area2

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def rank(vectors,p):
    return len(field_basis(tuple(range(len(vectors))),vectors.__getitem__,p)[0])


def chain(t,m):
    coefficients=[comb(t-1,r) for r in range(t)]
    pure={('U',r):{r:1} for r in range(m+1)}
    pure.update({('L',r):{r+k:coefficients[k] for k in range(m-r+1)}
                 for r in range(1,m+1)})
    return pure,{r:coefficients[r] for r in range(m+1)}


def born(pure,head,prices,p):
    for rho in sorted(set(prices.values())):
        vectors=[v for e,v in pure.items() if prices[e]<=rho]
        if rank(vectors+[head],p)==rank(vectors,p):
            return rho
    raise AssertionError('The uppers span.')


def formula(t,j,D,rho,theta=Q(1,3)):
    m=t-1-j;V=D+rho
    return j*rho+(theta+(m//2)*V if m%2==0 else ((m+1)//2)*V)


def prefix_controls():
    fixtures=[]
    for p,t in product((0,101),range(2,9)):
        D=t*t;rho=Q(12*t,12*t+1);theta=Q(1,3);V=D+rho
        for m in range(t):
            pure,head=chain(t,m)
            for u,s in product(range(m+2),range(m+1)):
                vectors=[pure['U',r] for r in range(u)]+[pure['L',r] for r in range(1,s+1)]
                assert rank(vectors,p)==min(m+1,u+s)
                assert rank(vectors+[head],p)==min(m+1,u+s+1)
                COUNTS['independent_complete_chain_prefix_ranks']+=1
            prices={e:j*rho+(theta if e[0]=='U' else 0)+e[1]*V
                    for e in pure for j in [t-1-m]}
            assert born(pure,head,prices,p)==formula(t,t-1-m,D,rho)
            COUNTS['independent_chain_head_births']+=1
        selected=[formula(t,t-1-2*r,D,rho) for r in range((t+1)//2)]
        assert all(b-a==D-rho>t-1 for a,b in zip(selected,selected[1:]))
        fixtures.append({'t':t,'D':D,'characteristic':p,'rho':str(rho),
                         'separated_costs':list(map(str,selected)),
                         'necessary_common_module_rank':len(selected)})
    return fixtures


def arm_sources(t,D):
    A=(t-1)*D+t;c=[comb(t-1,r) for r in range(t)]
    def image(q,j):
        return {(q+r*D)*t+j+r:c[r] for r in range(t-j) if q+r*D<A}
    pure={('U',q,j):{q*t+j:1} for q,j in product(range(A),range(t))}
    pure.update({('L',q,j):image(q,j) for q,j in product(range(t,A),range(t))})
    return A,pure,{(i,j):image(i,j) for i,j in product(range(t),repeat=2)}


def original_arm_controls():
    out=[]
    for p,t in [(0,2),(101,2),(0,3),(101,3),(0,4),(101,4),(101,5)]:
        D=t*t;A,pure,heads=arm_sources(t,D);rho=Q(12*t,12*t+1);theta=Q(1,3)
        prices={e:(theta if e[0]=='U' else 0)+e[1]+e[2]*rho for e in pure}
        for j in range(t):
            H=heads[0,j];previous=[v for (i,k),v in heads.items() if i+k*rho<j*rho]
            alpha=formula(t,j,D,rho)
            for closed in (False,True):
                vectors=previous+[v for e,v in pure.items()
                                   if prices[e]<alpha or closed and prices[e]==alpha]
                r=rank(vectors,p)
                assert rank(vectors+[H],p)==r+(not closed)
                COUNTS['independent_original_arm_birth_sides']+=1
        out.append({'t':t,'A':A,'D':D,'characteristic':p,
                    'all_original_pure_sources_retained':True})
    return out


def module_controls():
    out=[]
    for p,t in [(0,2),(101,3),(0,4),(101,5)]:
        D=t*t;A=(t-1)*D+t;rho=Q(12*t,12*t+1);theta=Q(1,3)
        c=[comb(t-1,r) for r in range(t)];rows=[]
        for j in range(t):
            row=[{} for _ in range(2*t)];row[j]={A:1};rows.append(tuple(row))
        for j in range(t):
            row=[{} for _ in range(2*t)];row[t+j]={0:1}
            for r in range(t-j):row[j+r]={r*D:c[r]}
            rows.append(tuple(row))
        shifts=tuple(theta+j*rho for j in range(t))+tuple(j*rho for j in range(t))
        rows,costs=reduce_rows(rows,shifts,p,A*t)
        columns={};prices={}
        head_order=sorted(range(t*t),key=lambda i:i//t+(i%t)*rho)
        place={old:new for new,old in enumerate(head_order)}
        for j,row in enumerate(rows):
            for k in range(t):
                columns[j,k]={place[q*t+v+k*t]:a for v in range(t)
                              for q,a in row[t+v].items() if q+k<t}
                prices[j,k]=costs[j]+k
        birth=assignment(columns,prices,p)
        for j in range(t):
            assert prices[birth[place[j]]]==formula(t,j,D,rho)
            COUNTS['independent_full_vector_lattice_head_minima']+=1
        separated=[formula(t,t-1-2*r,D,rho) for r in range((t+1)//2)]
        assert all(any(d<=a<=d+t-1 for d in costs) for a in separated)
        assert all(sum(d<=a<=d+t-1 for a in separated)<=1 for d in costs)
        out.append({'t':t,'A':A,'D':D,'characteristic':p,'actual_module_rank':2*t,
                    'necessary_cost_preserving_rank':len(separated),
                    'shifted_row_costs':list(map(str,costs))})
    return out


def profile_decision(t,D,w):
    A=(t-1)*D+t;a=A+t;N=set();orders=[];profile=list(fixed_profile(t,a))
    for arm in (0,1):
        for q,j in product(range(A),range(t)):
            if q>=D and j:continue
            m=min(t-1-j,(A-1-q)//D);labels=[]
            for r in range(m+1):
                labels.append((1,t+q+r*D,j+r))
                if q+r*D>=t:labels.append((0,t+q+r*D,t+j+r))
            if arm:labels=[(e[0],e[2],e[1]) for e in labels]
            order=tuple(sorted(labels,key=lambda e:(dot(w,e),e)));orders.append(order)
            for e,f in zip(order,order[1:]):
                v=tuple(b-a for a,b in zip(e,f))
                if not dot(w,v):raise ValueError('Pure price tie.')
                N.add(v)
            for e in order[:m+1]:
                profile=[a+b for a,b in zip(profile,e)]
    def gain(arm,m):
        options=[(1,r*D,r-t) for r in range(m+1)]+[(0,r*D,r) for r in range(1,m+1)]
        if arm:options=[(v[0],v[2],v[1]) for v in options]
        return sorted(options,key=lambda v:dot(w,v))[m]
    Y=[gain(0,t-1-j) for j in range(t)];Z=[gain(1,t-1-i) for i in range(t)]
    winners=[]
    for i,j in product(range(t),repeat=2):
        options={(0,0,0),Y[j],Z[i]};winner=max(options,key=lambda v:dot(w,v));winners.append(winner)
        for v in options:
            n=tuple(a-b for a,b in zip(winner,v))
            if any(n):
                if not dot(w,n):raise ValueError('Gain tie.')
                N.add(n)
        profile=[a-b for a,b in zip(profile,winner)]
    assert all(dot(w,q)>0 for q in N)
    return (tuple(orders),tuple(winners)),tuple(profile),N


def whole_fan(t,D,p):
    A=(t-1)*D+t;a=A+t;n=4*A*t+2*t*t;L=2*n+1
    queue=[((1,L,L*L),None)];cells={};areas=defaultdict(Q);vertices=defaultdict(list)
    while queue:
        w,previous=queue.pop();key,profile,N=profile_decision(t,D,w);assert key!=previous
        if key in cells:continue
        P=polygon(N);cells[key]=(profile,P);areas[profile]+=area2(P);vertices[profile].extend(P)
        for v,u in zip(P,P[1:]+P[:1]):
            if any(v[k]==u[k]==0 for k in range(3)):continue
            normal=next(q for q in N if dot(q,v)==dot(q,u)==0)
            nw,_=exact_adjacent_weight(n,v,u,normal);queue.append((nw,key))
    assert sum(areas.values())==1
    fan={v:hull(P) for v,P in vertices.items()}
    assert all(area2(fan[v])==area for v,area in areas.items())
    c=[comb(t-1,r) for r in range(t)]
    data,co,cols=build(t,a,((D,1),(D,1)),(c,c),p,True)
    reference=complete_fan(cols,n,p);translation=fixed_profile(t,a)
    ref={tuple(translation[k]+sum(e[k] for e in B) for k in range(3)):hull(P)
         for B,P in reference.items()}
    assert set(ref)==set(fan) and all(set(ref[v])==set(fan[v]) for v in fan)
    COUNTS['independent_complete_original_quotient_fans']+=1
    return {'t':t,'A':A,'D':D,'length':n,'characteristic':p,'actual_cones':len(fan),
            'refinement_patches':len(cells),'exact_total_area':'1','whole_original_fan_equal':True}


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    prior_paths=['results/research/hook_scalar_outer_lattice_checks.json',
                 'results/research/hook_two_orbit_checks.json']
    files=set()
    for f in prior_paths:
        r=json.loads((ROOT/f).read_text());assert r['status']=='passed'
        for name,h in r['source_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
        files|=set(r['source_sha256'])|{f}
    prefixes=prefix_controls();arms=original_arm_controls();modules=module_controls()
    fans=[whole_fan(2,4,0),whole_fan(3,9,101)]
    t,m=3,2;pure,head=chain(t,m);rho=Q(1,19);V=9+rho;theta=Q(1,3)
    prices={e:(theta if e[0]=='U' else 0)+e[1]*V for e in pure}
    actual=born(pure,head,prices,2)
    assert actual==2*V and actual!=theta+V
    guard={'t':3,'characteristic':2,'polynomial':'(1+T)^2=1+T^2',
           'actual_birth':str(actual),'characteristic_zero_formula':str(theta+V),
           'scope':'The binomial nonvanishing hypothesis cannot be omitted.'}
    files|={'docs/transverse_lattice_rank_obstruction.md','scripts/run_transverse_lattice_rank_checks.py'}
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_transverse_lattice_rank_checks.py',
        'uniform_result':{'family':'D=t^2, A=(t-1)D+t, W_Y=W_Z=(1+q^D v)^(t-1), zero outers, characteristic zero or p>2t.',
            'positive_weight':'rho=12t/(12t+1), theta=1/3; both arm upper shifts and all head gains are positive.',
            'exact_head_birth':'For m=t-1-j, alpha_j=j*rho+(theta+(m/2)*(D+rho) if m even else ((m+1)/2)*(D+rho)); 0<theta<D+rho.',
            'common_module_rank_lower_bound':'ceil(t/2) for one cost-preserving polynomial module with a common shifted degree and constraints modulo q^t.',
            'natural_module_rank_upper_bound':'2t for raw birth costs; hence order-optimal for that representation model.',
            'length_relation':'n=4t^4-4t^3+6t^2; necessary rank Omega(n^(1/4)).',
            'global_fan_upper_bound':'O(t^4)=O(n), by the saved steep two-monomial theorem.',
            'scope':'Actual one-arm head filters. No assertion that every birth gives a coherent original wall or that all algorithms require this module rank.'},
        'counts':dict(COUNTS),'complete_chain_prefix_fixtures':prefixes,'original_arm_fixtures':arms,
        'full_vector_lattice_fixtures':modules,'complete_original_fans':fans,'small_characteristic_guard':guard,
        'limits':['Classical binomial Toeplitz determinants, predictable degree and monomial-chain splitting are attributed.',
                  'This rules out a constant-rank extension of one representation model; it is not a new extremal fan bound or polynomial lower excess.',
                  'Exact strongest-current novelty, unrestricted G_3 and final publication assessment remain open.'],
        'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/transverse_lattice_rank_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':dict(COUNTS),'source_hashes':len(files),'fans':fans}),flush=True)


if __name__=='__main__':main()
