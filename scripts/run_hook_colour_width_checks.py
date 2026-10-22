#!/usr/bin/env python3
"""Uniform hook colour-width bounds and exact twelve-cone sparse controls."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product, combinations
from pathlib import Path
import hashlib
import json
import sys

from run_structural_complexity_checks import clip, area2, dot, determinant
from run_square_zero_extension_checks import extension, matrix_controls
from run_projection_tag_checks import fan as direct_fan, compare_fans
from run_generic_hook_checks import model as full_hook_model

ROOT = Path(__file__).resolve().parents[1]
PROOF = 'docs/hook_colour_width.md'
TRIANGLE = [tuple(Q(i==j) for i in range(3)) for j in range(3)]

def blocks(A,t,h):
 a=A-t
 X=(1,0,0);Y=(0,a,0);Z=(0,0,a)
 return ((X,Y),(X,Z),(X,Y,(0,t-h,t)),(X,Z,(0,t,t-h)),((1,h,0),(0,a+h,0),(1,0,h),(0,0,a+h),(0,t,t)))
def encoded_fan(A,t,h):
 G=blocks(A,t,h)
 choices=[[(i,) for i in range(len(g))] for g in G[:4]]+[[ij for ij in combinations(range(5),2) if ij not in ((0,1),(2,3))]]
 cells=[]
 for C in product(*choices):
  poly=list(TRIANGLE)
  for g,options,c in zip(G,choices,C):
   score=tuple(sum(g[i][j] for i in c) for j in range(3))
   for other in options:
    q=tuple(sum(g[i][j] for i in other)-score[j] for j in range(3))
    poly=clip(poly,q,1)
  if poly and area2(poly)>0:cells.append((C,poly))
 assert sum(area2(p) for c,p in cells)==1
 return cells


def basis_from_choices(A,t,h,C):
 D={(i,j) for i,j in product(range(A),repeat=2) if i<t or j<t}
 B={(0,*v) for v in D}|{(1,i,j) for i,j in product(range(t),repeat=2)}
 for i in range(t,h):
  for j in range(t):
   B.add(((1,i,j),(0,A+i-t,j))[C[0][0]])
   B.add(((1,j,i),(0,j,A+i-t))[C[1][0]])
 for i in range(h+t,A):
  for j in range(t):
   B.add(((1,i,j),(0,A+i-t,j),(0,t+i-h,t+j))[C[2][0]])
   B.add(((1,j,i),(0,j,A+i-t),(0,t+j,t+i-h))[C[3][0]])
 for a,b in product(range(t),repeat=2):
  options=((1,h+a,b),(0,A+h+a-t,b),(1,a,h+b),(0,a,A+h+b-t),(0,t+a,t+b))
  B.update(options[i] for i in C[4])
 assert len(B)==2*len(D)
 return frozenset(B)


def parameters(t):
    A = t*t+1
    m = 2*t*A-t*t
    L = A-3*t-1
    h = 2*t+1+L//2
    r = A-h-t
    return A,m,h,r


def candidate_controls(counts):
    rows = []
    for t in (4,5,6):
        A,m,h,r = parameters(t)
        D,index,G,allowed,params,columns,fixed = full_hook_model(t)
        assert len(D)==m
        assert sum(1+i+j for i,j in D)==t*A*(A+t)-t**3
        counts['independent_reference_degree_sums'] += 1
        U = sorted((a for a in columns if not a[0]),key=lambda a:(sum(a),a))
        for s in range(A):
            actual = sum(a[1]+a[2]<=s for a in U)
            assert actual<=2*t*(s+1)
            counts['full_parameter_candidate_strip_counts'] += 1
        for j,a in enumerate(U[:m],1):
            assert Q(j,2*t)-1<=a[1]+a[2]
            counts['independent_sorted_candidate_degree_bounds'] += 1
        rows.append({'thickness':t,'length':2*m,'full_support_plain_candidates':len(U)})
    return rows


def unit_minor_controls(counts):
    rows = []
    for t in (4,5,8,16,32):
        A,m,h,r = parameters(t)
        assert 2*t<h<A-t and r>0
        images = set()
        total_gain = 0
        for i,j in product(range(t,A-h),range(t)):
            for a,b in ((i,j),(j,i)):
                u = (t+a,t+b)
                assert max(u)<A
                # Independent sparse multiplication of f1, followed by reduction
                # in the original planar monomial ideal.
                image = {(h+a,b):1,(a,h+b):1}
                image = {v:c for v,c in image.items()
                         if v[0]<A and v[1]<A and (v[0]<t or v[1]<t)}
                assert len(image)==1
                v = next(iter(image))
                assert v not in images
                images.add(v)
                total_gain += 1+sum(v)-sum(u)
                counts['exact_unit_competitor_column_images'] += 1
        expected = ({(i,j) for i,j in product(range(h+t,A),range(t))} |
                    {(j,i) for i,j in product(range(h+t,A),range(t))})
        assert images==expected and len(images)==2*t*r
        L = t*t-3*t
        G = 2*t*r*(h-2*t+1)
        assert total_gain==G and G>=2*t*(L//2)*((L+1)//2)
        rows.append({'thickness':t,'length':2*m,'h':h,'replacement_columns':len(images),
                     'unit_minor_determinant_absolute_value':1,'balanced_cost_saving':G,
                     'scope':'Explicit sparse column identity; no large generic fan traversal.'})
    return rows


def uniform_margin_controls(counts):
    upper = Q(32,31)-Q(841,2048)+Q(1,2*32**4)
    lower = Q(193969,262144)
    assert upper<Q(2,3)<Q(7,10)<lower
    for t in (32,33,64,127,1000):
        A,m,h,r = parameters(t)
        L = t*t-3*t
        G = 2*t*(L//2)*((L+1)//2)
        C0 = t*A*(A+t)-t**3
        assert C0-G<Q(2,3)*t**5
        q = Q(7*m,8)
        assert q*q/(4*t)-m>Q(7,10)*t**5
        assert Q(m,8)-t*t>=Q(m,16)
        assert 4*A*m+1>4*A*m
        counts['exact_large_parameter_margin_controls'] += 1
    return {'uniform_upper_margin':str(upper),'uniform_lower_margin':str(lower),
            'scope':'The rational constants prove the uniform margins; selected integer checks do not replace that proof.'}


def main():
    if not __debug__:
        raise RuntimeError('Exact assertions must remain enabled.')
    before = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT/'scripts').glob('*.py')}
    before[PROOF] = hashlib.sha256((ROOT/PROOF).read_bytes()).hexdigest()
    counts = Counter()
    candidate_rows = candidate_controls(counts)
    unit_rows = unit_minor_controls(counts)
    margins = uniform_margin_controls(counts)
    # The only nontrivial atomic block is a rank-two regular triangle
    # with two pairs of parallel ground columns.
    triangle = ((1,1,0,0,1),(0,0,1,1,1))
    minors = [determinant([[triangle[i][j] for j in J] for i in range(2)])
              for J in combinations(range(5),2)]
    assert all(d in (-1,0,1) for d in minors)
    assert sum(d!=0 for d in minors)==8
    counts['integral_atomic_maximal_minors'] += len(minors)
    fixtures = []
    for A,t,h in ((5,1,3),(15,2,6),(18,3,8),(17,4,11)):
        assert 2*t<h<A-t
        encoded = encoded_fan(A,t,h)
        assert len(encoded)==12
        predicted = {basis_from_choices(A,t,h,C):poly for C,poly in encoded}
        m = 2*t*A-t*t
        ks = [sum(a[0] for a in B) for B in predicted]
        assert set(ks)=={t*t,t*h,2*t*h,t*(A+h-t)}
        D = tuple((i,j) for i,j in product(range(A),repeat=2) if i<t or j<t)
        G = ((A,0),(t,t),(0,A))
        tails = ({(t,0):1},{(h,0):1,(0,h):1},{(0,t):1})
        for p in (0,2,101):
            n,ev,E = extension(D,D,G,tails,p)
            matrix_controls(n,ev,tuple((x,*v) for x in range(2) for v in D),E,p,counts)
            got = direct_fan(n,ev,p,counts)
            compare_fans(got,predicted,counts)
            fixtures.append({'a':A,'thickness':t,'h':h,'characteristic':p,
                             'length':n,'entire_positive_fan':len(got),
                             'upper_count_values':sorted(set(ks))})
            print(json.dumps(fixtures[-1]),flush=True)
    compressed = []
    for t in (8,16,32,64,100):
        A,m,h,r = parameters(t)
        encoded = encoded_fan(A,t,h)
        multiplicities = (t*(h-t),t*(h-t),t*r,t*r,t*t)
        ks = [t*t+sum(M*sum(blocks(A,t,h)[g][i][0] for i in C[g])
                       for g,M in enumerate(multiplicities)) for C,poly in encoded]
        assert len(encoded)==12 and set(ks)=={t*t,t*h,2*t*h,t*(A+h-t)}
        compressed.append({'thickness':t,'length':2*m,'encoded_cones':12,
                           'minimum_upper_count':min(ks),'maximum_upper_count':max(ks),
                           'scope':'Compressed regular-block model; not an expanded large-algebra traversal.'})
        counts['large_compressed_block_fans'] += 1
    loaded = {Path(m.__file__).resolve() for m in tuple(sys.modules.values())
              if getattr(m,'__file__',None) and Path(m.__file__).resolve().parent==ROOT/'scripts'}
    loaded.add(Path(__file__).resolve())
    sources = [PROOF]+sorted(str(p.relative_to(ROOT)) for p in loaded)
    assert all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==before[f] for f in sources)
    report = {'all_checks_passed':True,'candidate_controls':candidate_rows,'unit_minor_controls':unit_rows,
              'uniform_margins':margins,'complete_original_fixtures':fixtures,'compressed_models':compressed,
              'counts':dict(counts),'source_sha256':{f:before[f] for f in sources},
              'loaded_source_stability_checked':True,
              'scope':'Uniform positive colour-width obstruction and exact twelve-cone sparse family. No uniform generic fan count, polynomial supergeneric separation, unrestricted G3 improvement or novelty claim.'}
    (ROOT/'results/research/hook_colour_width_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_checks_passed':True,'counts':dict(counts)}),flush=True)


if __name__=='__main__':
    main()
