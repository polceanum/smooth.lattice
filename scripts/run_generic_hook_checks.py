#!/usr/bin/env python3
"""Replay complete support fans and all-characteristic generic hook certificates.

The data are finite certificates, not a growing-thickness count theorem.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import gcd, isqrt
from pathlib import Path
import hashlib
import json
import sys

from run_structural_complexity_checks import dot, clip, area2, determinant
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_square_zero_extension_checks import extension, matrix_controls

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT/'results/research/generic_hook_checks.json'
PROOF = 'docs/generic_hook_fans.md'
TRIANGLE = [tuple(Q(i == j) for i in range(3)) for j in range(3)]


def model(t):
    assert t >= 2
    a = t*t+1
    D = tuple((i,j) for i,j in product(range(a),repeat=2) if i<t or j<t)
    index = {v:i for i,v in enumerate(D)}
    G = ((a,0),(t,t),(0,a))
    allowed = (tuple(v for v in D if v[0]>=t or v[1]>=a-t),
               tuple(v for v in D if v[0]>=t or v[1]>=t),
               tuple(v for v in D if v[0]>=a-t or v[1]>=t))
    parameters = tuple((g,*v) for g,A in enumerate(allowed) for v in A)
    parameter_index = {s:i for i,s in enumerate(parameters)}
    squares = tuple(tuple(x+y for x,y in zip(g,h)) for g,h in product(G,repeat=2))
    outer = tuple(v for v in product(range(2*a),repeat=2)
                  if not any(all(x<=y for x,y in zip(g,v)) for g in squares))
    columns = {(1,*v):{index[v]:None} for v in D}
    for v in outer:
        if v in index:
            continue
        choices = []
        for g,head in enumerate(G):
            if all(x<=y for x,y in zip(head,v)):
                shift = tuple(x-y for x,y in zip(v,head))
                choices.append({index[r]:parameter_index[(g,*tail)]
                                for tail in allowed[g]
                                if (r:=tuple(x+y for x,y in zip(tail,shift))) in index})
        assert choices and all(c==choices[0] for c in choices)
        if choices[0]:
            columns[(0,*v)] = choices[0]
    fixed = frozenset((0,*v) for v in D)
    return D,index,G,allowed,parameters,columns,fixed


def matching_basis(order, columns, rank):
    mate, B = {}, []
    def augment(a,seen):
        for i in columns[a]:
            if i in seen:
                continue
            seen.add(i)
            if i not in mate or augment(mate[i],seen):
                mate[i] = a
                return True
        return False
    for a in order:
        if augment(a,set()):
            B.append(a)
        if len(B)==rank:
            return frozenset(B),mate
    raise AssertionError('Support ground has deficient rank.')


def exchange_normals(B, columns, rank):
    _,mate = matching_basis(B,columns,rank)
    normals = set()
    for a in columns:
        if a in B:
            continue
        todo,seen,exchange = list(columns[a]),set(),set()
        while todo:
            i = todo.pop()
            if i in seen:
                continue
            seen.add(i)
            b = mate[i]
            exchange.add(b)
            todo.extend(set(columns[b])-seen)
        for b in exchange:
            q = tuple(x-y for x,y in zip(a,b))
            d = gcd(*q)
            normals.add(tuple(x//d for x in q))
    return sorted(normals)


def complete_fan(columns, fixed, rank, counts):
    queue = [((1,2*rank+1,(2*rank+1)**2),None)]
    cells = {}
    while queue:
        w,old = queue.pop()
        B,_ = matching_basis(sorted(columns,key=lambda a:(dot(w,a),a)),columns,rank)
        if old is not None:
            assert B!=old
        if B in cells:
            continue
        N = exchange_normals(B,columns,rank)
        assert all(dot(w,q)>0 for q in N)
        poly = list(TRIANGLE)
        for q in N:
            poly = clip(poly,q,1)
        assert area2(poly)>0
        full = fixed|B
        assert all(tuple(v-(j==i) for j,v in enumerate(a)) in full
                   for a in full for i in range(3) if a[i])
        cells[B] = poly
        counts['whole_support_matroid_cones'] += 1
        for a,b in zip(poly,poly[1:]+poly[:1]):
            if any(a[i]==b[i]==0 for i in range(3)):
                continue
            q = next(q for q in N if dot(q,a)==dot(q,b)==0)
            nw,_ = exact_adjacent_weight(2*rank,a,b,q)
            queue.append((nw,B))
    assert sum(map(area2,cells.values()))==1
    counts['complete_support_fans'] += 1
    return cells


def bareiss(A):
    A = [list(row) for row in A]
    if not A:
        return 1
    sign,old,n = 1,1,len(A)
    for k in range(n-1):
        i = next((i for i in range(k,n) if A[i][k]),None)
        if i is None:
            return 0
        if i!=k:
            A[k],A[i] = A[i],A[k]
            sign = -sign
        pivot = A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                v = A[i][j]*pivot-A[i][k]*A[k][j]
                assert v%old==0
                A[i][j] = v//old
            A[i][k] = 0
        old = pivot
    return sign*A[-1][-1]


def prime_factors(a):
    assert a>0
    factors,p = [],2
    while p*p<=a:
        if a%p==0:
            factors.append(p)
            while a%p==0:
                a//=p
        p += 1
    if a>1:
        factors.append(a)
    return factors


def poly_rem(a,f,p):
    a = list(a)
    while a and not a[-1]:
        a.pop()
    while len(a)>=len(f):
        c,k = a[-1]*pow(f[-1],-1,p)%p,len(a)-len(f)
        for i,x in enumerate(f):
            a[k+i] = (a[k+i]-c*x)%p
        while a and not a[-1]:
            a.pop()
    return a


def poly_multiply(a,b,f,p):
    out = [0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j] = (out[i+j]+x*y)%p
    return poly_rem(out,f,p)


class Field:
    """Exact polynomial quotient; all nonzero elements are checked invertible."""
    def __init__(self,p,d,f):
        assert d>=1 and p>1 and all(p%i for i in range(2,isqrt(p)+1))
        assert len(f)==d+1 and f[-1]==1 and all(0<=x<p for x in f)
        self.p,self.d,self.q,self.modulus = p,d,p**d,tuple(f)
        digits = [tuple((a//(p**i))%p for i in range(d)) for a in range(self.q)]
        encode = lambda a:sum(x*p**i for i,x in enumerate(a))
        self.add = [[encode(tuple((x+y)%p for x,y in zip(a,b))) for b in digits] for a in digits]
        self.neg = [encode(tuple(-x%p for x in a)) for a in digits]
        self.mul = [[encode(poly_multiply(a,b,f,p)) for b in digits] for a in digits]
        self.inv = [0]+[self.mul[a].index(1) for a in range(1,self.q)]
        assert all(self.mul[a][self.inv[a]]==1 for a in range(1,self.q))

    def determinant(self,A):
        A,n,out = [list(r) for r in A],len(A),1
        for j in range(n):
            i = next((i for i in range(j,n) if A[i][j]),None)
            if i is None:
                return 0
            if i!=j:
                A[j],A[i] = A[i],A[j]
                out = self.neg[out]
            a = A[j][j]
            out = self.mul[out][a]
            inv = self.inv[a]
            A[j] = [self.mul[x][inv] for x in A[j]]
            for i in range(j+1,n):
                a = A[i][j]
                if a:
                    a = self.neg[a]
                    A[i] = [self.add[x][self.mul[a][y]] for x,y in zip(A[i],A[j])]
        return out


def minor_labels(B,index,columns):
    rows = sorted(set(range(len(index)))-{index[a[1:]] for a in B if a[0]})
    cols = sorted(a for a in B if not a[0])
    assert len(rows)==len(cols)
    return [[columns[a].get(i) for a in cols] for i in rows]


def evaluate_minor(labels,coefficients):
    return [[0 if s is None else coefficients[s] for s in row] for row in labels]


def check_family(family,counts):
    t = family['thickness']
    D,index,G,allowed,parameters,columns,fixed = model(t)
    m = len(D)
    assert family['length']==2*m and family['parameters']==[list(a) for a in parameters]
    cells = complete_fan(columns,fixed,m,counts)
    Bs = [frozenset(map(tuple,B)) for B in family['support_fan_bases']]
    assert len(Bs)==len(set(Bs)) and set(Bs)==set(cells)
    minors = [minor_labels(B,index,columns) for B in Bs]
    gcds = [0]*len(Bs)
    for trial,e in enumerate(family['integer_evaluations']):
        assert len(e['coefficients'])==len(parameters) and len(e['determinants'])==len(Bs)
        for i,(L,claimed) in enumerate(zip(minors,e['determinants'])):
            if claimed is None:
                assert gcds[i]==1
                continue
            A = evaluate_minor(L,e['coefficients'])
            det = bareiss(A)
            assert det==claimed
            if (t==2 and trial==0) or (trial==0 and i%97==0):
                assert determinant(A)==det
                counts['independent_rational_Gaussian_determinants'] += 1
            gcds[i] = gcd(gcds[i],det)
            counts['exact_integer_determinant_evaluations'] += 1
    assert gcds==family['integer_gcds'] and all(gcds)
    factors = [prime_factors(a) for a in gcds]
    assert factors==family['prime_divisors']
    obligations = {(i,p) for i,P in enumerate(factors) for p in P}
    used_fields = []
    for e in family['finite_field_certificates']:
        F = Field(e['characteristic'],e['degree'],e['modulus'])
        used_fields.append([F.p,F.q])
        counts['independently_validated_finite_fields'] += 1
        for v in e['evaluations']:
            coefficients = v['coefficients']
            assert len(coefficients)==len(parameters) and all(0<=a<F.q for a in coefficients)
            for i,claimed in v['nonzero_determinants']:
                assert (i,F.p) in obligations
                det = F.determinant(evaluate_minor(minors[i],coefficients))
                assert det==claimed and det!=0
                obligations.remove((i,F.p))
                counts['exact_content_prime_exclusions'] += 1
    assert not obligations
    counts['primitive_coherent_minor_polynomials'] += len(Bs)
    values = dict(zip(parameters,family['independent_reference_coefficients']))
    tails = [{a:values[(g,*a)] for a in A} for g,A in enumerate(allowed)]
    n,ev,ground = extension(D,D,G,tails,0)
    matrix_controls(n,ev,tuple((x,*a) for x in range(2) for a in D),ground,0,counts)
    for a,c in columns.items():
        expected = ({m+i:1 for i in c} if a[0] else {m+i:family['independent_reference_coefficients'][s] for i,s in c.items()})
        assert ev(a)==expected
        counts['independent_original_extension_column_equalities'] += 1
    return {'thickness':t,'length':2*m,'generic_initial_ideals_in_every_characteristic':len(Bs),
            'sum_of_coherent_minor_degrees':sum(map(len,minors)),
            'largest_minor_degree':max(map(len,minors)),
            'integer_content_primes_excluded':sorted({p for P in factors for p in P}),
            'finite_fields_used':used_fields,'scope':'Generic coefficients over an infinite field or a sufficiently large finite extension; not every specialization and not a growing-thickness theorem.'}


def compact_report(report):
    report = dict(report)
    report['families'] = [dict(f) for f in report['families']]
    replacements = {}
    for i,f in enumerate(report['families']):
        for key in ('support_fan_bases','integer_evaluations'):
            marker = '__CERTIFICATE_'+str(i)+'_'+key+'__'
            replacements[json.dumps(marker)] = '[\n'+',\n'.join('      '+json.dumps(v,separators=(',',':')) for v in f[key])+'\n    ]'
            f[key] = marker
    out = json.dumps(report,indent=2)
    for key,value in replacements.items():
        assert out.count(key)==1
        out = out.replace(key,value)
    return out+'\n'


def main():
    if not __debug__:
        raise RuntimeError('Exact assertions must remain enabled.')
    before = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT/'scripts').glob('*.py')}
    before[PROOF] = hashlib.sha256((ROOT/PROOF).read_bytes()).hexdigest()
    report = json.loads(TARGET.read_text())
    counts = Counter()
    rows = []
    for family in report['families']:
        row = check_family(family,counts)
        rows.append(row)
        print(json.dumps(row),flush=True)
    loaded = {Path(m.__file__).resolve() for m in tuple(sys.modules.values())
              if getattr(m,'__file__',None) and Path(m.__file__).resolve().parent==ROOT/'scripts'}
    loaded.add(Path(__file__).resolve())
    sources = [PROOF]+sorted(str(p.relative_to(ROOT)) for p in loaded)
    assert all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==before[f] for f in sources)
    report['verification'] = {'all_checks_passed':True,'models':rows,'counts':dict(counts),
                              'source_sha256':{f:before[f] for f in sources},
                              'loaded_source_stability_checked':True,
                              'scope':'Exact finite all-characteristic generic fan certificates. No unrestricted exponent improvement, growing-thickness lower law, exact novelty or goal completion.'}
    encoded = compact_report(report)
    assert json.loads(encoded)==report
    TARGET.write_text(encoded)
    print(json.dumps({'all_checks_passed':True,'counts':dict(counts)}),flush=True)


if __name__=='__main__':
    main()
