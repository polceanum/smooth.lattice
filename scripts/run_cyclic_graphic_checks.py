#!/usr/bin/env python3
"""Independent incidence, arithmetic-lift, hierarchy and complete-fan checks."""
import ast
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import gcd, lcm
from pathlib import Path
import hashlib
import json

from run_univariate_branch_gluing_checks import (
    Frame, apply, clean, dense_rank, jsonable, polynomial,
)
from run_projection_tag_checks import compare_fans, fan
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_planar_dyadic_bound_checks import convex_hull
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross

ROOT=Path(__file__).resolve().parents[1]
TRIANGLE=[tuple(Q(i==j) for i in range(3)) for j in range(3)]


class DSU:
    def __init__(self,n):self.parent=list(range(n));self.components=n
    def find(self,u):
        while self.parent[u]!=u:
            self.parent[u]=self.parent[self.parent[u]];u=self.parent[u]
        return u
    def join(self,u,v):
        u,v=self.find(u),self.find(v)
        if u==v:return False
        self.parent[v]=u;self.components-=1;return True


def canonical(q):
    q=tuple(Q(c) for c in q);scale=lcm(*(c.denominator for c in q))
    values=tuple(int(c*scale) for c in q)
    common=gcd(*map(abs,values))
    if not common:return None
    values=tuple(c//common for c in values)
    return values if next(c for c in values if c)>0 else tuple(-c for c in values)


def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])


def model(p,q,a,b):
    n=p+q;axis=p+q-gcd(p,q)
    @lru_cache(None)
    def ev(e):
        i,j,k=e
        if j and k:return {}
        if j:return {(i+a*j)%p:1}
        if k:return {p+(i+b*k)%q:1}
        return {i%p:1,p+i%q:1}
    def edge(e):
        i,j,k=e
        assert not j*k
        if j:return ((i+a*j)%p,n)
        if k:return (n,p+(i+b*k)%q)
        return (i%p,p+i%q)
    E=tuple(e for e in hyperbolic_cross(n,3) if not e[1]*e[2] and
            (e[0]<axis if not e[1] and not e[2] else
             e[0]<p and e[1]<=p if e[1] else e[0]<q and e[2]<=q))
    return dict(p=p,q=q,a=a,b=b,n=n,axis=axis,ev=ev,edge=edge,E=E)


def mst(H,w):
    D=DSU(H['n']+1);out=[]
    for e in sorted(H['E'],key=lambda e:(dot(e,w),e)):
        if D.join(*H['edge'](e)):out.append(e)
    assert len(out)==H['n'] and D.components==1
    return frozenset(out)


def circuit_normals(H,B):
    adj=defaultdict(list)
    for e in B:
        u,v=H['edge'](e);adj[u].append((v,e));adj[v].append((u,e))
    normals=set()
    for e in H['E']:
        if e in B:continue
        u,v=H['edge'](e);stack=[(u,None,())]
        while stack:
            x,parent,path=stack.pop()
            if x==v:break
            stack.extend((y,x,path+(f,)) for y,f in adj[x] if y!=parent)
        else:raise AssertionError('No tree path.')
        normals.update(tuple(e[i]-f[i] for i in range(3)) for f in path)
    return normals


def graph_fan(H,counts):
    queue=[(1,H['n']+1,(H['n']+1)**2)];cells={}
    while queue:
        w=queue.pop();B=mst(H,w)
        if B in cells:continue
        normals=circuit_normals(H,B)
        assert all(dot(q,w)>0 for q in normals)
        poly=list(TRIANGLE)
        for q in sorted(normals):poly=clip(poly,q,1)
        assert area2(poly)>0
        cells[B]=poly
        for u,v in zip(poly,poly[1:]+poly[:1]):
            if any(u[i]==v[i]==0 for i in range(3)):continue
            q=next(q for q in normals if dot(q,u)==dot(q,v)==0)
            neighbor,_=exact_adjacent_weight(H['n'],u,v,q);queue.append(neighbor)
    assert sum(area2(poly) for poly in cells.values())==1
    counts['complete_graph_circuit_fans']+=1
    return cells


def hierarchy(H):
    n=H['n'];D=DSU(n)
    nodes=[dict(vertices=(i,),birth=None,death=None,children=()) for i in range(n)]
    owner={i:i for i in range(n)}
    for i in range(H['axis']):
        u,v=H['edge']((i,0,0));ru,rv=D.find(u),D.find(v)
        assert ru!=rv
        left,right=owner[ru],owner[rv]
        nodes[left]['death']=nodes[right]['death']=i
        vertices=tuple(sorted(nodes[left]['vertices']+nodes[right]['vertices']))
        nodes.append(dict(vertices=vertices,birth=i,death=None,children=(left,right)))
        D.join(u,v);owner[D.find(u)]=len(nodes)-1
    assert D.components==gcd(H['p'],H['q'])
    assert len(nodes)==H['n']+H['axis']
    return nodes


def describe(H,C):
    p,q=H['p'],H['q'];assert p<=q
    Y={i for i in C['vertices'] if i<p};Z={i-p for i in C['vertices'] if i>=p}
    ell=0 if C['birth'] is None else C['birth']+1
    if not Y:
        assert len(Z)==1
        return dict(kind='Z_singleton',r=next(iter(Z)),L=1)
    if ell<q:
        assert len(Y)==1;r=next(iter(Y))
        expected={s for s in range(ell) if s%p==r}
        assert Z==expected
        return dict(kind='pre_wrap',r=r,L=len(Z))
    stop=ell-q
    r=next(r for r in sorted(Y) if (r-q)%p>=stop)
    expected={(r+t*q)%p for t in range(len(Y))}
    assert expected==Y and Z=={s for s in range(q) if s%p in Y}
    return dict(kind='post_wrap',r=r,L=len(Y))


def lifted_membership(H,desc,branch,i,j):
    """Return an integer witness for (G10)/(G11) or its smaller case."""
    p,q,a,b=(H[k] for k in ('p','q','a','b'))
    r,L,kind=desc['r'],desc['L'],desc['kind']
    if kind=='Z_singleton':
        if branch==1 or (i+b*j-r)%q:return None
        return (i,j,(i+b*j-r)//q)
    if kind=='pre_wrap':
        if branch==1:
            return None if (i+a*j-r)%p else (i,j,(i+a*j-r)//p)
        for t in range(L):
            value=i+b*j-r-t*p
            if not value%q:
                return (i,j,t,value//q)
        return None
    if branch==1:
        for t in range(L):
            value=i+a*j-r-t*q
            if not value%p:
                ell=value//p
                assert i+a*j==r+t*q+ell*p
                return (i,j,t,ell)
        return None
    for t in range(L):
        s=(i+b*j)%q
        if not (s-r-t*q)%p:
            ell=(s-r-t*q)//p;m=(i+b*j-s)//q
            assert 0<=r+t*q+ell*p<q and i+b*j==r+t*q+ell*p+m*q
            return (i,j,t,ell,m)
    return None


def reduced_lift_points(H,desc,branch):
    """Enumerate the three-variable polytope, without residue membership tests."""
    p,q,a,b=(H[k] for k in ('p','q','a','b'))
    r,L,kind=desc['r'],desc['L'],desc['kind']
    def interval(base,step,cap):
        return range(-(base//step),(cap-1-base)//step+1)
    points={}
    def save(i,j,witness):
        assert len(witness)<=3
        points[(i,j)]=witness
    if kind=='Z_singleton':
        if branch==1:return points
        for j in range(1,q+1):
            for m in interval(r-b*j,q,q):save(r+m*q-b*j,j,(j,m))
    elif kind=='pre_wrap':
        if branch==1:
            for j in range(1,p+1):
                for ell in interval(r-a*j,p,p):save(r+ell*p-a*j,j,(j,ell))
        else:
            for j in range(1,q+1):
                for t in range(L):
                    base=r+t*p-b*j
                    for m in interval(base,q,q):save(base+m*q,j,(j,t,m))
    elif branch==1:
        for j in range(1,p+1):
            for t in range(L):
                base=r+t*q-a*j
                for ell in interval(base,p,p):save(base+ell*p,j,(j,t,ell))
    else:
        # G12: r+ell*p is in [-(L-1)*q, q), and h absorbs t+m.
        lower=-(L-1)*q
        for j in range(1,q+1):
            for ell in interval(r-lower,p,L*q):
                s=r+ell*p
                assert lower<=s<q
                for h in interval(s-b*j,q,q):
                    i=s+h*q-b*j
                    t=-(s//q);m=h-t
                    assert 0<=t<L and 0<=s+t*q<q
                    assert i+b*j==r+t*q+ell*p+m*q
                    save(i,j,(j,ell,h))
    return points


def component_candidates(H,C,desc,counts):
    candidates=[];witnesses=[]
    for branch,period,slope in ((1,H['p'],H['a']),(2,H['q'],H['b'])):
        points=[];reduced=reduced_lift_points(H,desc,branch)
        saved=0
        for i in range(period):
            for j in range(1,period+1):
                vertex=(i+slope*j)%period+(H['p'] if branch==2 else 0)
                actual=vertex in C['vertices']
                witness=lifted_membership(H,desc,branch,i,j)
                assert actual==(witness is not None)
                counts['independent_integer_lift_memberships']+=1
                if actual:
                    points.append((i,j))
                    if saved<2:
                        witnesses.append(dict(branch=branch,coordinates=witness,
                                              reduced_coordinates=reduced[(i,j)]))
                        saved+=1
        assert set(points)==set(reduced)
        counts['independent_three_variable_lift_points']+=len(reduced)
        counts['independent_entire_reduced_lift_equalities']+=1
        hull=convex_hull(points) if points else []
        candidates.extend((i,j if branch==1 else 0,j if branch==2 else 0) for i,j in hull)
    assert candidates
    return tuple(sorted(candidates)),witnesses


def inventory(H,nodes,counts):
    planes=set();records=[]
    for C in nodes:
        desc=describe(H,C)
        candidates,witnesses=component_candidates(H,C,desc,counts)
        local=set()
        for a in candidates:
            poly=list(TRIANGLE)
            for b in candidates:
                if a==b:continue
                poly=clip(poly,tuple(b[i]-a[i] for i in range(3)),1)
                if not poly:break
            if not poly:continue
            for u,v in zip(poly,poly[1:]+poly[:1]):
                if any(u[i]==v[i]==0 for i in range(3)):continue
                local.add(canonical(cross(u,v)))
        assert len(local)<=3*len(candidates)+3
        planes.update(local)
        for level in (C['birth'],C['death']):
            if level is None:continue
            for a in candidates:
                normal=canonical((a[0]-level,a[1],a[2]))
                if normal and min(normal)<0<max(normal):planes.add(normal)
        records.append(dict(vertices=C['vertices'],birth=C['birth'],death=C['death'],
                            children=C['children'],description=desc,
                            projected_integer_hull_candidates=candidates,
                            integer_lift_witnesses=witnesses))
    counts['arithmetic_hierarchy_components']+=len(nodes)
    return tuple(sorted(planes)),records


def reconstruct(H,nodes,records,w):
    choices=[min(R['projected_integer_hull_candidates'],key=lambda e:(dot(e,w),e))
             for R in records]
    prices=[dot(e,w) for e in choices];out=set()
    for j,C in enumerate(nodes):
        if ((C['birth'] is None or prices[j]>C['birth']*w[0]) and
                (C['death'] is None or prices[j]<C['death']*w[0])):
            out.add(choices[j])
        if C['children']:
            left,right=C['children'];level=C['birth']
            if prices[left]>level*w[0] or prices[right]>level*w[0]:out.add((level,0,0))
    assert len(out)==H['n']
    return frozenset(out)


def generic_weight(poly,planes):
    for base in range(2,2+len(planes)*len(poly)+2):
        w=tuple(sum(base**j*v[i] for j,v in enumerate(poly)) for i in range(3))
        if all(dot(q,w) for q in planes):return w
    raise AssertionError('No generic positive point in a full-dimensional cell.')


def compile_inventory(H,nodes,records,planes,cells,counts):
    patches=[list(TRIANGLE)]
    for q in planes:
        patches=[part for poly in patches for sign in (-1,1)
                 if (part:=clip(poly,q,sign))]
    assert sum(area2(poly) for poly in patches)==1
    grouped=defaultdict(list)
    for poly in patches:
        w=tuple(sum(v[i] for v in poly) for i in range(3))
        B=reconstruct(H,nodes,records,w)
        assert B==mst(H,w)
        assert all(dot(q,v)>=0 for q in circuit_normals(H,B) for v in poly)
        grouped[B].extend(poly)
    recovered={}
    for B,vertices in grouped.items():
        hull=convex_hull([(v[0],v[1]) for v in vertices])
        recovered[B]=[(x,y,1-x-y) for x,y in hull]
    compare_fans(cells,recovered,counts)
    counts['independent_inventory_arrangement_cells']+=len(patches)
    counts['entire_inventory_compilations']+=1
    return len(patches)


def axis_polynomial(p,q):
    g=gcd(p,q)
    value=[0]*(p+q+1);value[0]=1;value[p]-=1;value[q]-=1;value[p+q]+=1
    divisor=[-1]+[0]*(g-1)+[1];answer=[0]*(p+q-g+1)
    for i in range(len(value)-1,g-1,-1):
        c=value[i];answer[i-g]=c
        for j,d in enumerate(divisor):value[i-g+j]-=c*d
    assert not any(value)
    return tuple(answer)


def algebra_controls(H,char,counts):
    p,q,a,b,n=(H[k] for k in ('p','q','a','b','n'))
    X=tuple({(i+1)%p:1} for i in range(p))+tuple({p+(i+1)%q:1} for i in range(q))
    Y=tuple({(i+a)%p:1} for i in range(p))+({},)*q
    Z=({},)*p+tuple({p+(i+b)%q:1} for i in range(q))
    M=(X,Y,Z);unit=H['ev']((0,0,0))
    for U,V in ((X,Y),(X,Z),(Y,Z)):
        assert all(apply(U,V[i],char)==apply(V,U[i],char) for i in range(n))
    assert dense_rank(tuple(H['ev']((i,p,0)) for i in range(p))+
                      tuple(H['ev']((i,0,q)) for i in range(q)),n,char)==n
    f=axis_polynomial(p,q)
    assert len(f)-1==H['axis'] and not polynomial(X,unit,f,char)
    assert dense_rank(tuple(H['ev']((i,0,0)) for i in range(H['axis'])),n,char)==H['axis']
    for e in H['E']:
        u,v=H['edge'](e)
        incidence=clean({i:c for i,c in ((u,1),(v,-1)) if i<n},char)
        transformed=clean({i:c if i<p else -c for i,c in H['ev'](e).items()},char)
        assert incidence==transformed
        value=unit
        for op,power in zip(M,e):
            for _ in range(power):value=apply(op,value,char)
        assert value==clean(H['ev'](e),char)
        counts['independent_incidence_and_multiplication_columns']+=1
    if n<=5:
        for mask in range(1<<len(H['E'])):
            subset=[e for j,e in enumerate(H['E']) if mask>>j&1]
            D=DSU(n+1)
            for e in subset:D.join(*H['edge'](e))
            assert dense_rank(tuple(H['ev'](e) for e in subset),n,char)==n+1-D.components
            counts['exhaustive_small_subset_graph_ranks']+=1
    return M,unit,f


def hashes():
    paths=['docs/cyclic_graphic_fans.md','scripts/run_cyclic_graphic_checks.py',
           'docs/cographic_order_ideals.md']
    pending=[paths[1]]
    while pending:
        path=pending.pop()
        for node in ast.walk(ast.parse((ROOT/path).read_text())):
            names=([node.module] if isinstance(node,ast.ImportFrom) and node.module else
                   [a.name for a in node.names] if isinstance(node,ast.Import) else [])
            for name in names:
                candidate='scripts/'+name.replace('.','/')+'.py'
                if (ROOT/candidate).is_file() and candidate not in paths:
                    paths.append(candidate);pending.append(candidate)
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def main():
    if not __debug__:raise RuntimeError('Assertions must be enabled.')
    before=hashes();counts=Counter();rows=[]
    # Validate every merge component, including unequal and noncoprime periods.
    for p in range(1,17):
        for q in range(p,17):
            H=model(p,q,0,0)
            for C in hierarchy(H):describe(H,C)
            counts['exhaustive_period_pair_hierarchies']+=1
    fixtures=((1,2,0,1),(2,2,1,1),(2,3,1,1),(3,5,2,3),(4,5,2,3),
              (5,7,2,3),(5,7,4,6),(5,5,2,3),(6,9,2,4),(7,9,3,4),
              (7,11,3,5),(8,9,3,4),(9,11,4,5),(11,13,5,6))
    for p,q,a,b in fixtures:
        H=model(p,q,a,b);nodes=hierarchy(H);planes,records=inventory(H,nodes,counts)
        cells=graph_fan(H,counts)
        for B,poly in cells.items():
            axes=sorted(e[0] for e in B if not e[1] and not e[2])
            assert axes==list(range(len(axes))) and len(axes)<=H['axis']
            w=generic_weight(poly,planes)
            assert reconstruct(H,nodes,records,w)==B==mst(H,w)
            for u,v in zip(poly,poly[1:]+poly[:1]):
                if any(u[i]==v[i]==0 for i in range(3)):continue
                assert canonical(cross(u,v)) in planes
            counts['actual_wall_inventory_and_tree_reconstructions']+=1
        arrangement_cells=None
        if H['n']<=8:arrangement_cells=compile_inventory(H,nodes,records,planes,cells,counts)
        for char in (0,2,3,101):
            M,unit,f=algebra_controls(H,char,counts)
            actual=fan(H['n'],H['ev'],char,counts)
            compare_fans(actual,cells,counts)
        rows.append(dict(p=p,q=q,a=a,b=b,length=H['n'],axis_rank=H['axis'],
                         minimal_axis_polynomial=f,characteristics=(0,2,3,101),
                         hierarchy=records,comparison_inventory=planes,cones=len(cells),
                         full_inventory_arrangement_cells=arrangement_cells,
                         all_standard_sets=[sorted(B) for B in sorted(cells,key=lambda B:sorted(B))]))
        print(json.dumps(dict(p=p,q=q,a=a,b=b,cones=len(cells),planes=len(planes))),flush=True)
    assert before==hashes(),'A source changed during verification.'
    report=dict(command='python3 -B scripts/run_cyclic_graphic_checks.py',status='passed',
                source_sha256=before,counts=dict(counts),uniform_proof='docs/cyclic_graphic_fans.md',
                exact_quotient_witnesses=rows,
                uniform_theorems={'forest_ground_model':'At most 2N merge components; O((N+sum_C vertices(Q_C))^2) unique positive-weight MSTs.',
                                  'integer_lift_specialization':'At most three integer variables after exact elimination and a fixed number of O(log n)-bit inequalities per component imply O(n^2 log^4(n+1)) complete initial ideals for every p,q,a,b and field.'},
                primary_sources=[
                    dict(url='https://www.math.uwaterloo.ca/~bico/papers/ippoints.pdf',scope='Cook--Hartmann--Kannan--McDiarmid Theorem 2.1, p.29; full pages 27--29 inspected. General rational polyhedra, at most m inequalities of row size phi, at most 2*m^D*(6*D^2*phi)^(D-1) integer-hull vertices. This is an established input.'),
                    dict(url='https://link.springer.com/article/10.1007/s00453-022-01024-9',scope='Eppstein abstract and Section 6: single-parameter lower bounds, including bounded-treewidth graphs, do not provide the required two-parameter periodic-weight count.'),
                    dict(url='https://doi.org/10.1006/jagm.1994.1019',scope='Fernandez-Baca--Slutzki publisher abstract inspected: one-parameter polynomial optimization on bounded-treewidth graphs; no two-parameter theorem is imported.')],
                limits=['The uniform logarithmic bound relies on the published integer-hull theorem, not a growth fit or finite maximum.',
                        'The integer-lift checker enumerates bounded small sets and is not an optimized integer-hull implementation.',
                        'The graph representation, Kruskal tests, integer-hull estimate and arrangement bound are classical; exact priority of the combined application is unassessed.',
                        'Every prime characteristic is covered by the incidence proof, including characteristics dividing either period. Finite field comparisons are supporting controls.',
                        'This closes this whole periodic direct-product family to polynomial separation, not arbitrary nongeneric yz=0 quotients or the unrestricted G3 endpoint.',
                        'No sharp logarithmic order, current novelty, publication significance, unrestricted complexity theorem or completed goal is claimed.'])
    (ROOT/'results/research/cyclic_graphic_checks.json').write_text(json.dumps(jsonable(report),indent=2)+'\n')
    print(json.dumps(dict(status='passed',counts=dict(counts))),flush=True)


if __name__=='__main__':main()
