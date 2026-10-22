#!/usr/bin/env python3
"""Fixed integer hook fibers: exact fans and classical specialization calibration.

Run from the repository root:
    python3 -B scripts/run_hook_specialization_checks.py

This is a bounded arithmetic control, not an asymptotic lower bound,
a new specialization theorem or completion of the extremal goal.
"""
from pathlib import Path
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import gcd, lcm, isqrt
import hashlib
import json
import sys

from run_coordinate_slice_budget_checks import hook_model, complete_hook, edges_from_polygons
from run_square_zero_extension_checks import extension, matrix_controls
from run_interleaving_research import minimal_excluded_border
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_structural_complexity_checks import area2, clip, dot, determinant
from run_projection_tag_checks import polygon_key
from run_signed_face_checks import field_basis, normals as whole_border_normals

ROOT = Path(__file__).resolve().parents[1]


def model(t, p, counts):
    n, _, _, integer_tails = hook_model(t, 271828, counts)
    R = t*t+1
    G = ((R,0),(t,t),(0,R))
    D = tuple((i,j) for i,j in product(range(R),repeat=2) if i<t or j<t)
    tails = tuple({a:1 for a in f} for f in integer_tails)
    nn, ev, ground = extension(D,D,G,tails,p)
    assert n==nn
    matrix_controls(n,ev,tuple((x,*a) for x in range(2) for a in D),ground,p,counts)
    fixed = frozenset((0,*a) for a in D)
    variable = tuple(a for a in ground if a not in fixed)
    m = len(D)
    def clean(v):
        return {j:(c%p if p else Q(c)) for j,c in v.items() if (c%p if p else c)}
    def col(a):
        v = ev(a)
        assert all(j>=m for j in v)
        return clean({j-m:c for j,c in v.items()})
    def eliminate(v,piv,co=None):
        v=dict(v); co={} if co is None else dict(co)
        while v:
            j=min(v)
            if j not in piv: break
            row, rep = piv[j]
            a=v[j]
            v=clean({k:v.get(k,0)-a*row.get(k,0) for k in v.keys()|row.keys()})
            co=clean({k:co.get(k,0)+a*rep.get(k,0) for k in co.keys()|rep.keys()})
        return v,co
    def greedy(w):
        piv={}; B=[]
        for a in sorted(variable,key=lambda a:(dot(w,a),a)):
            v,_=eliminate(col(a),piv)
            if not v: continue
            j=min(v); inv=pow(v[j],-1,p) if p else 1/v[j]
            piv[j]=(clean({k:c*inv for k,c in v.items()}),{})
            B.append(a)
            if len(B)==m: return fixed|frozenset(B)
        raise AssertionError('rank')
    def normals(B):
        piv={}
        for a in sorted(B-fixed):
            v,co=eliminate(col(a),piv)
            assert v
            j=min(v); inv=pow(v[j],-1,p) if p else 1/v[j]
            rep={k:-c*inv for k,c in co.items()};rep[a]=inv
            piv[j]=(clean({k:c*inv for k,c in v.items()}),clean(rep))
        out=set()
        for a in minimal_excluded_border(tuple(B)):
            v,co=eliminate(col(a),piv)
            assert not v
            for b in co:
                q=tuple(x-y for x,y in zip(a,b));d=gcd(*q)
                out.add(tuple(x//d for x in q))
        return sorted(out)
    return n,greedy,normals,tails


def full(t,p,counts):
    n,greedy,normals,tails=model(t,p,counts)
    queue=[((1,n+1,(n+1)**2),None)];cells={};edges=set();queries=0
    while queue:
        w,old=queue.pop();B=greedy(w);queries+=1
        if old is not None:
            assert old!=B;edges.add(frozenset((old,B)))
        if B in cells:continue
        N=normals(B)
        assert all(dot(q,w)>0 for q in N)
        poly=[tuple(Q(i==j) for i in range(3)) for j in range(3)]
        for q in N: poly=clip(poly,q,1)
        assert area2(poly)>0;cells[B]=poly
        for a,b in zip(poly,poly[1:]+poly[:1]):
            if any(a[i]==b[i]==0 for i in range(3)):continue
            q=next(q for q in N if dot(q,a)==dot(q,b)==0)
            nw,_=exact_adjacent_weight(n,a,b,q)
            queue.append((nw,B))
    assert sum(map(area2,cells.values()))==1
    assert queries==1+2*len(edges)
    assert edges_from_polygons(cells,counts)==edges
    return n,cells,edges


def bareiss(A):
    A=[list(row) for row in A]
    if not A:return 1
    sign=old=1;n=len(A)
    for k in range(n-1):
        i=next((i for i in range(k,n) if A[i][k]),None)
        if i is None:return 0
        if i!=k:A[k],A[i]=A[i],A[k];sign=-sign
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value=A[i][j]*pivot-A[i][k]*A[k][j]
                assert value%old==0
                A[i][j]=value//old
            A[i][k]=0
        old=pivot
    return sign*A[-1][-1]


def primes(v):
    v=abs(v);out=[];p=2
    while p*p<=v:
        if v%p==0:
            out.append(p)
            while v%p==0:v//=p
        p+=1
    if v>1:out.append(v)
    return out


def integer_model(t, tails):
    R=t*t+1
    D=tuple((i,j) for i,j in product(range(R),repeat=2) if i<t or j<t)
    n,ev,ground=extension(D,D,((R,0),(t,t),(0,R)),tails,0)
    return D,n,ev,ground


def basis_minor(B,D,ev):
    m=len(D);index={a:i for i,a in enumerate(D)}
    U=sorted(a for a in B if a[0]==0 and a[1:] not in index)
    T={(a[1],a[2]) for a in B if a[0]}
    S=sorted(set(D)-T)
    A=[[int(ev(u).get(m+index[v],0)) for u in U] for v in S]
    return S,U,A


def independent_cone(B,ev,p):
    N,relations=whole_border_normals(tuple(sorted(B)),ev,p)
    poly=[tuple(Q(i==j) for i in range(3)) for j in range(3)]
    for q in N:poly=clip(poly,q,1)
    return poly,len(relations)


def serialize_tails(tails):
    return [[{'exponent':a,'coefficient':c} for a,c in sorted(f.items())]
            for f in tails]


def report_json(out):
    """Keep each complete cone on one line without changing the JSON data."""
    key='whole_rational_fan_and_integer_determinants'
    payload=dict(out)
    payload['length_102_specialization']=dict(out['length_102_specialization'])
    cones=payload['length_102_specialization'][key]
    marker='__COMPLETE_RATIONAL_CONES__'
    payload['length_102_specialization'][key]=marker
    encoded=json.dumps(payload,indent=2)
    token=json.dumps(marker)
    assert encoded.count(token)==1
    block='[\n'+',\n'.join('      '+json.dumps(c,separators=(',',':')) for c in cones)+'\n    ]'
    encoded=encoded.replace(token,block)
    assert json.loads(encoded)==json.loads(json.dumps(out))
    return encoded+'\n'


def main():
    counts=Counter();rows=[];saved={};tails_by_t={}
    for t in (2,3):
        fans={}
        for p in ((2,0,3,101) if t==2 else (2,0,101)):
            n,cells,edges=full(t,p,counts)
            fans[p]=cells
            rows.append({'thickness':t,'length':n,'characteristic':p,
                         'initial_ideals':len(cells),'undirected_walls':len(edges),
                         'complete_area':str(sum(map(area2,cells.values())))})
            print('complete integer hook:',t,p,n,len(cells),len(edges),flush=True)
        # A second implementation uses bitsets instead of sparse field arithmetic.
        nn,binary,edges,_=complete_hook(t,counts)
        assert nn==n and set(binary)==set(fans[2])
        for B,poly in binary.items():
            assert polygon_key(poly)==polygon_key(fans[2][B])
            counts['independent_binary_entire_cone_equalities']+=1
        _,_,_,tails=model(t,0,counts)
        tails_by_t[t]=serialize_tails(tails)
        saved[t]=fans

    t=3;fans=saved[t];_,_,_,tails=model(t,0,counts)
    D,n,ev,ground=integer_model(t,tails)
    certificates=[];exceptional=set();missing={p:0 for p in (2,101,77719)}
    witness_candidates={};coverage=Q(0)
    for number,(B,poly) in enumerate(fans[0].items()):
        S,U,A=basis_minor(B,D,ev)
        det=bareiss(A)
        assert det and determinant(A)==det
        counts['independent_rational_Gaussian_minor_equalities']+=1
        factors=primes(det);remaining=abs(det)
        for p in factors:
            assert all(p%d for d in range(2,isqrt(p)+1))
            while remaining%p==0:remaining//=p
        assert remaining==1
        exceptional.update(factors)
        independent,num=independent_cone(B,ev,0)
        assert polygon_key(independent)==polygon_key(poly)
        coverage+=area2(independent)
        counts['independent_whole_column_cone_equalities']+=1
        counts['independent_original_border_relations']+=num
        for p in missing:
            lost=det%p==0
            missing[p]+=lost
            if p in fans:
                assert (B in fans[p]) == (not lost)
                counts['modular_basis_survival_equalities']+=1
            if lost and p not in witness_candidates:
                witness_candidates[p]=(B,poly,S,U,A,det)
        certificates.append({'standard_set':sorted(B),
                             'positive_cone':[[str(c) for c in v] for v in poly],
                             'integer_basis_determinant':det,'prime_divisors':factors})
    assert coverage==1
    assert len(certificates)==709 and len(exceptional)==80 and max(exceptional)==77719
    assert missing=={2:545,101:5,77719:1}

    # Test the entire stable fan at a prime below several exceptional primes.
    stable=79
    assert stable not in exceptional
    for B,poly in fans[0].items():
        independent,num=independent_cone(B,ev,stable)
        assert polygon_key(independent)==polygon_key(poly)
        counts['independent_stable_prime_entire_cone_equalities']+=1

    witnesses=[]
    for p,(B,poly,S,U,A,det) in sorted(witness_candidates.items()):
        w=tuple(sum(v[i] for v in poly) for i in range(3))
        den=lcm(*(a.denominator for a in w));w=tuple(int(a*den) for a in w)
        d=gcd(*w);w=tuple(a//d for a in w)
        _,qgreedy,qnormals,_=model(t,0,counts)
        _,pgreedy,_,_=model(t,p,counts)
        pB=pgreedy(w)
        assert qgreedy(w)==B and pB!=B and min(w)>0
        assert all(dot(q,w)>0 for q in qnormals(B))
        assert len(field_basis(tuple(sorted(B)),ev,0)[0])==n
        assert len(field_basis(tuple(sorted(B)),ev,p)[0])<n
        assert len(field_basis(tuple(sorted(pB)),ev,p)[0])==n
        assert determinant(A)==det and det%p==0
        witnesses.append({'characteristic':p,'strict_positive_weight':w,
                          'rational_standard_set':sorted(B),'modular_standard_set':sorted(pB),
                          'minor_rows':S,'minor_columns':U,'integer_minor':A,
                          'integer_determinant':det,'determinant_prime_divisors':primes(det)})
        counts['independent_strict_coherent_rank_loss_witnesses']+=1

    expected={(2,2):(47,82),(2,0):(47,85),(2,3):(47,85),(2,101):(47,85),
              (3,2):(546,1041),(3,0):(709,1505),(3,101):(707,1495)}
    for row in rows:
        assert (row['initial_ideals'],row['undirected_walls'])==expected[row['thickness'],row['characteristic']]
    source_paths={Path(__file__).resolve(),ROOT/'docs/square_zero_extension_fans.md'}
    for module in tuple(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name:
            path=Path(name).resolve()
            if path.suffix=='.py' and path.parent==ROOT/'scripts':source_paths.add(path)
    out={'command':'python3 -B scripts/run_hook_specialization_checks.py',
         'all_checks_passed':True,
         'scope':'Bounded arithmetic calibration in one fixed integral flat-hook family. '
                 'Rank specialization and the greedy-minimum criterion are classical. '
                 'No asymptotic growth inference or new specialization theorem is claimed.',
         'integer_family':{'base_ideal':'(y^(t^2+1),y^t*z^t,z^(t^2+1))',
                           'extension':'(x^2, g_i-x*f_i)',
                           'reference_basis':'D disjoint-union xD; length 4t^3-2t^2+4t',
                           'seed':271828,'tails_by_thickness':tails_by_t},
         'complete_fans':rows,'counters':dict(counts),
         'length_102_specialization':{
             'exact_exceptional_characteristics_for_labelled_fan':sorted(exceptional),
             'rational_bases_lost_modulo':missing,
             'stable_prime_entire_fan_replay':{'characteristic':stable,'initial_ideals':709,
                                              'complete_area':str(coverage)},
             'rational_bases_shared_with_modular':{
                 str(p):len(set(fans[0])&set(fans[p])) for p in (2,101)},
             'whole_rational_fan_and_integer_determinants':certificates,
             'strict_rank_loss_witnesses':witnesses,
             'criterion':'In a finite free integral algebra with an integral monomial '
                         'reference basis, the labelled positive fan changes precisely '
                         'at prime divisors of determinants of rationally coherent bases. '
                         'This does not imply that the number of cones changes at every such prime.'},
         'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in sorted(source_paths)},
         'unrestricted_G3_improved':False,'polynomial_lower_separation_proved':False,
         'novelty_established':False,'goal_complete':False}
    destination=ROOT/'results/research/hook_specialization_checks.json'
    destination.write_text(report_json(out))
    print('verified',len(certificates),'rational cones;',len(exceptional),
          'exceptional primes; stable characteristic',stable,flush=True)


if __name__=='__main__':
    main()
