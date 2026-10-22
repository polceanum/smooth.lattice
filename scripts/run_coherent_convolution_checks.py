#!/usr/bin/env python3
"""Exact coherent support failure, symbolic ranks and free power transport."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib
import json
import sys

from run_generic_hook_checks import model, matching_basis, exchange_normals, minor_labels
from run_square_zero_extension_checks import extension, matrix_controls
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import clip, area2, dot

ROOT = Path(__file__).resolve().parents[1]
PROOF = 'docs/coherent_convolution_rank_gap.md'
PARAMETERS = ((0,0,15),(0,1,13),(1,1,12),(1,2,13),(1,14,1))
NAMES = ('A','B','D','E','H')
WEIGHT = (39,43,109)


def fixture():
    D,index,G,allowed,P,columns,fixed = model(4)
    ids = {P.index(p):j for j,p in enumerate(PARAMETERS)}
    columns = {a:{r:None if p is None else ids[p] for r,p in c.items()
                  if p is None or p in ids} for a,c in columns.items()}
    columns = {a:c for a,c in columns.items() if c}
    R = {(0,15),(0,16),(16,3)} | set(product(range(1,4),range(12,17)))
    C = {(0,17+i,j) for i,j in product(range(4),range(2))}
    C |= {(0,4,4+j) for j in range(4)}
    C |= {(0,4+i,4+j) for i,j in product(range(1,3),range(3))}
    B = frozenset({(1,*v) for v in D if v not in R}|C)
    assert len(D)==120 and len(R)==len(C)==18 and len(B)==120
    return D,index,G,columns,fixed,R,C,B


def terms(L):
    """List signed nonzero determinant terms, independently of elimination."""
    out,k = [],len(L)
    def visit(i,used,mon,sign):
        if i==k:
            out.append((sign,tuple(sorted(mon))))
            return
        for j in range(k):
            if j not in used and L[i][j] is not None:
                visit(i+1,used|{j},mon+(L[i][j],),
                      sign*(-1)**sum(c>j for c in used))
    visit(0,set(),(),1)
    return out


def laplace(L):
    """Polynomial-valued subset recurrence, with early zero removal."""
    @lru_cache(None)
    def rec(mask):
        i = mask.bit_count()
        if i==len(L):
            return {():1}
        out = Counter()
        for j,s in enumerate(L[i]):
            if s is None or mask & (1<<j):
                continue
            sign = (-1)**((mask>>(j+1)).bit_count())
            for mon,c in rec(mask|(1<<j)).items():
                out[tuple(sorted(mon+(s,)))] += sign*c
        return {m:c for m,c in out.items() if c}
    return rec(0)


def assignment(rows,columns,weight):
    """Independent integer Hungarian optimization plus checked LP dual."""
    ground = sorted(columns)
    n,m = len(rows),len(ground)
    costs = [dot(weight,a) for a in ground]
    inf = 10**12
    A = [[costs[j] if i in columns[a] else inf for j,a in enumerate(ground)]
         for i in range(n)]
    u,v,p,way = [0]*(n+1),[0]*(m+1),[0]*(m+1),[0]*(m+1)
    for i in range(1,n+1):
        p[0],j0 = i,0
        minimum,used = [inf]*(m+1),[False]*(m+1)
        while True:
            used[j0] = True
            i0,delta,j1 = p[j0],inf,0
            for j in range(1,m+1):
                if used[j]:
                    continue
                cur = A[i0-1][j-1]-u[i0]-v[j]
                if cur<minimum[j]:
                    minimum[j],way[j] = cur,j0
                if minimum[j]<delta:
                    delta,j1 = minimum[j],j
            assert delta<inf//2
            for j in range(m+1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minimum[j] -= delta
            j0 = j1
            if not p[j0]:
                break
        while j0:
            j1 = way[j0]
            p[j0],j0 = p[j1],j1
    pairs = [(p[j]-1,j-1) for j in range(1,m+1) if p[j]]
    assert len(pairs)==n and len({i for i,j in pairs})==n
    assert all(i in columns[ground[j]] for i,j in pairs)
    assert all(a<=0 for a in v[1:])
    assert all(u[i+1]+v[j+1]<=costs[j]
               for j,a in enumerate(ground) for i in columns[a])
    primal = sum(costs[j] for i,j in pairs)
    assert primal==sum(u[1:])+sum(v[1:])
    return frozenset(ground[j] for i,j in pairs), {
        'ground':ground,'matching':pairs,'row_dual':u[1:],
        'column_dual':v[1:],'primal_and_dual_cost':primal}


def support_controls(counts):
    D,index,G,columns,fixed,R,C,B = fixture()
    order = sorted(columns,key=lambda a:(dot(WEIGHT,a),a))
    actual,mate = matching_basis(order,columns,120)
    assert actual==B
    independent,dual = assignment(D,columns,WEIGHT)
    assert independent==B
    counts['independent_assignment_primal_dual_certificates'] += 1
    N = exchange_normals(B,columns,120)
    assert len(N)==46 and min(dot(WEIGHT,q) for q in N)==3
    assert all(dot(WEIGHT,q)>0 for q in N)
    poly = [tuple(Q(i==j) for i in range(3)) for j in range(3)]
    for q in N:
        poly = clip(poly,q,1)
    expected = [(Q(0),Q(2,7),Q(5,7)),(Q(0),Q(5,18),Q(13,18)),
                (Q(5,14),Q(5,28),Q(13,28)),(Q(3,10),Q(1,5),Q(1,2))]
    assert set(poly)==set(expected) and area2(poly)==Q(11,2520)
    facets = [q for q in N if any(dot(q,a)==dot(q,b)==0
              for a,b in zip(poly,poly[1:]+poly[:1]))]
    assert set(facets)=={(-1,-11,5),(0,-5,2),(0,13,-5)}
    counts['whole_cone_exchange_normals_verified'] += len(N)
    _,_,_,_,_,full_columns,_ = model(4)
    full_min,_ = matching_basis(sorted(full_columns,key=lambda a:(dot(WEIGHT,a),a)),
                                full_columns,120)
    full_gain = sum(dot(WEIGHT,a) for a in B)-sum(dot(WEIGHT,a) for a in full_min)
    assert full_min!=B and full_gain==4438
    counts['full_parameter_scope_controls'] += 1
    full = fixed|B
    assert all(tuple(v-(j==i) for j,v in enumerate(a)) in full
               for a in full for i in range(3) if a[i])
    for a,c in columns.items():
        for axis in (1,2):
            b = tuple(x+(j==axis) for j,x in enumerate(a))
            target = {index[w]:s for r,s in c.items()
                      if (w:=tuple(x+(j==axis-1) for j,x in enumerate(D[r]))) in index}
            assert columns.get(b,{})==target
            counts['exact_support_successor_columns'] += 1
    return {'basis':sorted(full),'missing_upper_rows':sorted(R),
            'selected_plain_columns':sorted(C),'weight':WEIGHT,
            'exchange_normals':N,'minimum_primitive_margin':3,
            'cone_vertices':[[str(x) for x in v] for v in poly],
            'doubled_area':str(area2(poly)),'assignment_dual':dual,
            'full_parameter_scope_control':{'same_basis':False,
                                            'strict_cost_improvement':full_gain}}


def symbolic_controls(counts):
    D,index,G,columns,fixed,R,C,B = fixture()
    L = minor_labels(B,index,columns)
    T = terms(L)
    expected = (0,)*7+(1,)+(2,)*8+(3,4)
    assert len(T)==2 and {s for s,m in T}=={-1,1}
    assert all(m==expected for s,m in T) and not laplace(L)
    M = [[v for j,v in enumerate(row) if j!=3]
         for i,row in enumerate(L) if i!=2]
    witness = (0,)*8+(2,)*6+(3,)*2+(4,)
    assert terms(M)==[(1,witness)] and laplace(M)=={witness:1}
    counts['independent_exact_symbolic_determinants'] += 4
    rows,cols,forced = list(range(18)),list(range(18)),[]
    while True:
        pair = next(((i,J[0]) for i in rows
                     if len(J:=[j for j in cols if L[i][j] is not None])==1),None)
        if pair is None:
            pair = next(((I[0],j) for j in cols
                         if len(I:=[i for i in rows if L[i][j] is not None])==1),None)
        if pair is None:
            break
        i,j = pair
        forced.append((i,j,L[i][j]))
        rows.remove(i)
        cols.remove(j)
    cycle = [[L[i][j] for j in cols] for i in rows]
    assert cycle==[[2,None,0,None],[None,2,1,None],
                   [3,None,None,0],[None,3,None,1]]
    assert len(forced)==14
    assert Counter(s for i,j,s in forced)=={0:6,2:7,4:1}
    counts['forced_determinant_leaf_expansions'] += len(forced)
    return {'parameter_names':NAMES,'parameter_labels':PARAMETERS,
            'minor_labels':L,'signed_permutation_terms':T,
            'rank17_minor':{'deleted_row':2,'deleted_column':3,
                            'signed_product':(1,witness)},
            'forced_expansions':forced,'remaining_cycle':cycle}


def algebra_controls(counts):
    D,index,G,columns,fixed,R,C,B = fixture()
    rows = []
    for characteristic,values in ((0,(1,1,1,1,1)),(0,(2,3,5,7,11)),
                                   (2,(1,1,1,1,1)),(3,(1,2,1,2,1)),
                                   (101,(2,3,5,7,11))):
        tails = [{}, {}, {}]
        for (g,i,j),v in zip(PARAMETERS,values):
            tails[g][i,j] = v
        n,ev,ground = extension(D,D,G,tails,characteristic)
        reference = tuple((x,*a) for x in range(2) for a in D)
        matrix_controls(n,ev,reference,ground,characteristic,counts)
        for a,c in columns.items():
            target = {120+r:(1 if s is None else values[s]) for r,s in c.items()}
            if characteristic:
                target = {r:v%characteristic for r,v in target.items() if v%characteristic}
            assert ev(a)==target
            counts['independent_original_extension_columns'] += 1
        basis,_ = field_basis(sorted(fixed|B),ev,characteristic)
        assert len(basis)==239
        optimum,_ = field_basis(sorted(ground,key=lambda a:(dot(WEIGHT,a),a)),ev,characteristic)
        assert len(optimum)==240 and frozenset(optimum)!=fixed|B
        cost = sum(dot(WEIGHT,a) for a in optimum)
        assert cost>sum(dot(WEIGHT,a) for a in fixed|B)
        rows.append({'characteristic':characteristic,'coefficients':values,
                     'support_basis_actual_rank':len(basis),
                     'actual_minimum_cost':cost,
                     'support_minimum_cost':sum(dot(WEIGHT,a) for a in fixed|B)})
        counts['independent_field_rank_and_minimum_controls'] += 1
    return rows


def power_controls(counts):
    D,index,G,columns,fixed,R,C,B = fixture()
    tails = [{(i,j):1 for g,i,j in PARAMETERS if g==k} for k in range(3)]
    _,base,_ = extension(D,D,G,tails,0)
    data = []
    for s in (1,2,3):
        ds = tuple(sorted((s*i+r,s*j+q) for i,j in D for r,q in product(range(s),repeat=2)))
        gs = tuple((s*i,s*j) for i,j in G)
        ts = [{(s*i,s*j):v for (i,j),v in f.items()} for f in tails]
        n,ev,ground = extension(ds,ds,gs,ts,0)
        reference = tuple((x,*a) for x in range(2) for a in ds)
        row_labels = [(j%s,k%s,120*x+index[j//s,k//s]) for x,j,k in reference]
        assert n==240*s*s and len(set(row_labels))==n
        for a in ground:
            x,j,k = a
            r,q = j%s,k%s
            independent = {(r,q,i):v for i,v in base((x,j//s,k//s)).items()}
            assert {row_labels[i]:v for i,v in ev(a).items()}==independent
            counts['independent_power_substitution_column_equalities'] += 1
        expanded = tuple(sorted((x,s*j+r,s*k+q)
                         for x,j,k in fixed|B for r,q in product(range(s),repeat=2)))
        rank,_ = field_basis(expanded,ev,0)
        assert len(rank)==239*s*s
        for a in expanded:
            x,j,k = a
            assert dot((39*s,43,109),a)==s*dot(WEIGHT,(x,j//s,k//s))+43*(j%s)+109*(k%s)
        data.append({'s':s,'length':n,'actual_rank':len(rank),'rank_deficit':s*s})
        counts['expanded_power_substitution_rank_controls'] += 1
    return data


def main():
    if not __debug__:
        raise RuntimeError('Exact assertions must remain enabled.')
    before = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT/'scripts').glob('*.py')}
    before[PROOF] = hashlib.sha256((ROOT/PROOF).read_bytes()).hexdigest()
    counts = Counter()
    support = support_controls(counts)
    symbolic = symbolic_controls(counts)
    fields = algebra_controls(counts)
    powers = power_controls(counts)
    loaded = {Path(m.__file__).resolve() for m in tuple(sys.modules.values())
              if getattr(m,'__file__',None) and Path(m.__file__).resolve().parent==ROOT/'scripts'}
    loaded.add(Path(__file__).resolve())
    sources = [PROOF]+sorted(str(p.relative_to(ROOT)) for p in loaded)
    assert all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==before[f] for f in sources)
    report = {'all_checks_passed':True,'support_cone':support,'symbolic_rank':symbolic,
              'field_controls':fields,'expanded_power_controls':powers,
              'counts':dict(counts),'source_sha256':{f:before[f] for f in sources},
              'scope':'All-field coherent rank obstruction on one five-parameter sparse stratum; free substitution gives a uniform linear deficit and preserves the full fan. No refutation of full-parameter generic hook survival, growing fan count, unrestricted exponent, novelty or goal completion.'}
    (ROOT/'results/research/coherent_convolution_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_checks_passed':True,'counts':dict(counts),'powers':powers}),flush=True)


if __name__=='__main__':
    main()
