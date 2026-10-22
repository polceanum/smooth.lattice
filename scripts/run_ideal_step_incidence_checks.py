#!/usr/bin/env python3
"""Exact adjacent point-ideal fans and arithmetic incidence controls."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd, isqrt, lcm
from pathlib import Path
import hashlib
import json
import sys

from run_projection_tag_checks import fan, point_model, polygon_key
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import cross, dot, hyperbolic_cross
from run_canonical_fan_deflation_checks import marked_basis
from run_algebraic_incidence_controls import pair_checks
from run_structural_complexity_checks import clip
from run_wall_coexistence_checks import complete_planar_fan
from run_coordinate_slice_budget_checks import edges_from_polygons, profile

ROOT = Path(__file__).resolve().parents[1]


def clean(poly):
    poly = list(dict.fromkeys(poly))
    while len(poly) > 2:
        keep = [a for i, a in enumerate(poly)
                if any(cross(tuple(x-y for x, y in zip(a, poly[i-1])),
                             tuple(x-y for x, y in zip(poly[(i+1) % len(poly)], a))))]
        if len(keep) == len(poly):
            break
        poly = keep
    return poly


def primitive(v, unoriented=False):
    den = lcm(*(a.denominator if isinstance(a, Q) else 1 for a in v))
    v = tuple(int(a*den) for a in v)
    d = gcd(*v)
    assert d
    v = tuple(a//d for a in v)
    if unoriented and next(a for a in v if a) < 0:
        v = tuple(-a for a in v)
    return v


def ceil_sqrt(h):
    r = isqrt(h)
    return r + (r*r < h)


def fan_budgets(cells):
    rays, walls, wall_counts = Counter(), {}, Counter()
    for poly in cells.values():
        P = clean(poly)
        for a in P:
            if min(a) > 0:
                rays[primitive(a)] += 1
        for a, b in zip(P, P[1:]+P[:1]):
            if any(a[i] == b[i] == 0 for i in range(3)):
                continue
            edge = tuple(sorted((a, b)))
            q = primitive(cross(a, b), True)
            assert min(q) < 0 < max(q)
            if edge in walls:
                assert walls[edge] == q
            walls[edge] = q
            wall_counts[edge] += 1
    assert set(wall_counts.values()) <= {2}
    U = sum(ceil_sqrt(max(w))*d for w, d in rays.items())
    W = sum(sum(map(abs, q)) for q in walls.values())
    assert U <= 2*W
    return rays, walls, U, W


def actual_edge_length(cells):
    """Recover actual profile jumps independently of primitive fan normals."""
    edges=edges_from_polygons(cells,Counter())
    return sum(sum(abs(a-b) for a,b in zip(profile(B),profile(C)))
               for B,C in map(tuple,edges))


def section_controls(cells,counts):
    """Segment counting and independent sections from every profile pair."""
    from run_planar_dyadic_bound_checks import convex_hull
    edges=edges_from_polygons(cells,Counter())
    moments={B:profile(B) for B in cells};points=tuple(moments.values())
    integral=0;checked=[]
    for axis in range(3):
        events=Counter()
        for B,C in map(tuple,edges):
            low,high=sorted((moments[B][axis],moments[C][axis]))
            if low<high:events[low]+=1;events[high]-=1
        levels=sorted(events);active=0;value=0
        for a,b in zip(levels,levels[1:]):
            active+=events[a];value+=active*(b-a)
        assert value==sum(abs(moments[B][axis]-moments[C][axis]) for B,C in map(tuple,edges))
        integral+=value;counts['independent_half_integer_segment_sums']+=1
        lo=min(a[axis] for a in points);hi=max(a[axis] for a in points)
        if lo==hi:continue
        other=tuple(i for i in range(3) if i!=axis)
        for k in sorted({lo,(lo+hi)//2,hi-1}):
            t=Q(2*k+1,2)
            def intersection(a,b):
                if not min(a[axis],b[axis])<t<max(a[axis],b[axis]):return None
                f=(t-a[axis])/(b[axis]-a[axis])
                return tuple(a[j]+f*(b[j]-a[j]) for j in other)
            full={q for a,b in combinations(points,2) if (q:=intersection(a,b)) is not None}
            hull=set(convex_hull(full))
            positive={q for B,C in map(tuple,edges)
                      if (q:=intersection(moments[B],moments[C])) is not None}
            assert positive<=hull
            counts['independent_full_rational_section_hulls']+=1
            checked.append({'axis':axis,'half_integer_level':str(t),
                            'positive_edge_crossings':len(positive),'full_section_vertices':len(hull)})
    assert integral==actual_edge_length(cells)
    return checked


def triples(N):
    for A in range(1, N-1):
        for B in range(1, N-A):
            yield A, B, N-A-B


def triple_ray(A, B, C):
    return primitive((B*C, A*C, A*B))


@lru_cache(None)
def S(N):
    return sum(ceil_sqrt(max(triple_ray(*t))) for t in triples(N))


def phi(n):
    return sum(gcd(a, n) == 1 for a in range(1, n+1))


def divisor_increment(n):
    return 3*sum(d*phi(d) for d in range(2, n+1) if n % d == 0)


def direct_W(n):
    return 3*sum((a+b)//gcd(a, b) for a in range(1, n)
                 for b in range(1, n+1-a))


def local_histogram(p, K, e):
    """Exact distribution of the two smallest capped valuations' sum."""
    if not K:
        return {0: Q(1)}
    if not e:
        result = {0: 1-Q(3, p*p)}
        for j in range(1, K):
            result[j] = 3*(1-Q(1, p*p))/p**(2*j)
        result[K] = Q(3, p**(2*K))
        return {j: v for j, v in result.items() if v}
    result = {0: 1-Q(1, p*p)}
    result.update({j+2: v/p**2 for j, v in local_histogram(p, K-1, e-1).items()})
    return result


def valuation(a, p, K):
    if not a:
        return K
    e = 0
    while e < K and a % p == 0:
        a //= p
        e += 1
    return e


def exact_local_weight(hist, p):
    # a+b/sqrt(p), represented by two rational numbers.
    out = [Q(0), Q(0)]
    for j, probability in hist.items():
        out[j % 2] += probability/p**(j//2)
    return tuple(out)


def local_controls(counts):
    rows = []
    for p, K in ((2, 5), (3, 3), (5, 2), (7, 2)):
        modulus = p**K
        values = [valuation(a, p, K) for a in range(modulus)]
        for N in range(modulus):
            hist = Counter()
            for a, b in product(range(modulus), repeat=2):
                vals = sorted((values[a], values[b], values[(N-a-b) % modulus]))
                hist[vals[0]+vals[1]] += 1
            hist = {j: Q(v, modulus**2) for j, v in hist.items()}
            expected = local_histogram(p, K, values[N])
            assert hist == expected and sum(hist.values()) == 1
            assert exact_local_weight(hist, p) == exact_local_weight(expected, p)
            counts['independent_complete_local_residue_distributions'] += 1
            counts['independent_residue_pairs'] += modulus**2
        rows.append({'prime': p, 'valuation_cap': K, 'modulus': modulus,
                     'all_total_residues_checked': modulus,
                     'weight_arithmetic': 'exact rational coefficients in Q(sqrt(p))'})
    # Rational certificates in the infinite-progression proof.
    assert Q(1, 3) < Q(3, 5)**2  # 1/sqrt(3)<3/5.
    assert Q(8, 9)/Q(6, 7)+Q(1, 27) == Q(29, 27)
    assert Q(1)/(1-Q(9, 200)) == Q(200, 191) < Q(29, 27)
    counts['exact_progression_gap_certificates'] += 3
    return rows


def arithmetic_controls(counts):
    rows = []
    previous = 0
    for n in range(1, 151):
        W = direct_W(n)
        inc = divisor_increment(n)
        assert W-previous == inc <= 3*n*(n-1)
        assert W == 3*sum(d*phi(d)*(n//d) for d in range(2, n+1))
        if n >= 2 and phi(n) == n-1:
            assert inc == 3*n*(n-1)
            counts['independent_prime_increment_equalities'] += 1
        counts['independent_totient_wall_recurrences'] += 1
        previous = W
    for N in (15, 45, 105, 315, 525, 735, 945):
        lower, center, upper = S(N-2), S(N), S(N+2)
        positive, negative = 3*(center-lower), 3*(upper-center)
        if N >= 105:
            assert positive > 0 > negative
        rows.append({'total': N, 'S_minus_two': lower, 'S_center': center,
                     'S_plus_two': upper,
                     'U_increment_at_length_N_minus_one': positive,
                     'U_increment_at_length_N_plus_one': negative})
        counts['independent_exact_integer_incidence_sums'] += 3
        print(rows[-1], flush=True)
    assert rows[2]['U_increment_at_length_N_minus_one'] == 71364
    assert rows[2]['U_increment_at_length_N_plus_one'] == -28866
    for N in range(3, 41):
        for d in range(2, N+1):
            actual = sum(a % d == b % d == 0 for a, b, c in triples(N))
            m = (N-1)//d
            assert actual == m*(m-1)//2
            counts['independent_divisible_pair_tail_identities'] += 1
    return rows


def actual_controls(counts):
    rows = []
    for p, maximum in ((0, 7), (101, 11)):
        previous = None
        for n in range(2, maximum+1):
            def inverse(a):
                return pow(a, -1, p) if p else Q(1, a)
            points = tuple((j, inverse(j), inverse(j-1)) for j in range(2, n+2))
            _, ev = point_model(points, p)
            # Separator F_{n-1}(x) divided by its last-point value.
            values = []
            for point in points:
                v = 1
                for j in range(2, n+1):
                    v *= point[0]-j
                values.append(v % p if p else v)
            assert all(v == 0 for v in values[:-1]) and values[-1]
            normalizer = inverse(values[-1])
            v = tuple(x*normalizer % p if p else x*normalizer for x in values)
            assert v == (0,)*(n-1)+(1,)
            for axis in range(3):
                assert all((points[j][axis]-points[-1][axis])*v[j] == 0 for j in range(n))
                counts['independent_actual_separator_eigenline_actions'] += 1
            for point in points:
                x, y, z = point
                for value in (x*y-1, (x-1)*z-1, y*z-z+y):
                    assert (value % p if p else value) == 0
                    counts['independent_actual_curve_relations'] += 1
            cells = fan(n, ev, p, counts)
            rays, walls, U, W = fan_budgets(cells)
            assert actual_edge_length(cells)==n**3-n
            predicted = {triple_ray(*t) for N in (n, n+1) for t in triples(N)}
            assert set(rays) == predicted and set(rays.values()) == {3}
            assert len(cells) == comb(n+1, 2) and len(rays) == (n-1)**2
            assert U == 3*(S(n)+S(n+1)) and W == direct_W(n)
            assert len(walls) == 3*comb(n, 2)
            counts['independent_entire_actual_ray_and_wall_inventories'] += 1
            for B, poly in cells.items():
                G, _ = marked_basis(B, ev, p)
                w = tuple(sum(a[i] for a in poly)/len(poly) for i in range(3))
                counts['independent_original_reduced_basis_S_pairs'] += pair_checks(G, w, p)
            if previous is not None:
                oldcells, oldrays, oldev = previous
                shared = {triple_ray(*t) for t in triples(n)}
                assert set(rays) & set(oldrays) == shared
                assert len(set(rays)-set(oldrays)) == comb(n, 2)
                assert len(set(oldrays)-set(rays)) == comb(n-2, 2)
                # Same global order on both actual point quotients.
                ground = tuple(hyperbolic_cross(n, 3))
                for C, poly in oldcells.items():
                    w = tuple(sum(a[i] for a in poly)/len(poly) for i in range(3))
                    B = frozenset(field_basis(sorted(ground, key=lambda a: (dot(w, a), a)), ev, p)[0])
                    assert C < B
                    counts['independent_actual_ideal_step_standard_set_inclusions'] += 1
                counts['independent_exact_adjacent_ray_intersections'] += 1
            previous = cells, rays, ev
            rows.append({'characteristic': p, 'length': n, 'initial_cones': len(cells),
                         'interior_rays': len(rays), 'undirected_wall_pieces': len(walls),
                         'U': U, 'W': W})
            print(rows[-1], flush=True)
    return rows


def product_wall_controls():
    """Whole actual point fans, independently compared with branch products."""
    counts=Counter();rows=[]
    cases=[(kind,r,p,tuple(range(1,r+1))) for kind in ('fixed_curve','generic_branch')
           for r in (2,3,4,5) for p in ((0,101) if kind=='fixed_curve' else (0,))]
    cases.append(('fixed_curve',12,101,(11,12)))
    for kind,r,p,ms in cases:
        inv=lambda j:pow(j,-1,p) if p else Q(1,j)
        last=tuple((j,inv(j) if kind=='fixed_curve' else j**r) for j in range(1,r+1))
        _,ev2=point_model(last,p)
        branch2=complete_planar_fan(r,ev2,p,counts)
        b=len(branch2)
        def planar_W(cells):
            moments=[tuple(sum(a[i] for a in C) for i in range(2)) for C,*_ in cells]
            return sum(sum(map(abs,primitive(tuple(x-y for x,y in zip(B,A)))))
                       for A,B in zip(moments,moments[1:]))
        c=planar_W(branch2)
        assert c<=r*(r-1)
        if kind=='fixed_curve':assert b==r and c==divisor_increment(r)//3
        previous=None
        for m in ms:
            points=(tuple((j,0,0) for j in range(1,r+1))+
                    tuple((j,inv(j),0) for j in range(1,m+1))+
                    tuple((x,0,z) for x,z in last))
            n,ev=point_model(points,p);assert n==2*r+m and len(set(points))==n
            for x,y,z in points:
                f=1;fm=1
                for j in range(1,r+1):f*=x-j
                for j in range(1,m+1):fm*=x-j
                values=[f,y*z,y*fm,x*y*y-y]
                values.append(x*z*z-z if kind=='fixed_curve' else z*z-x**r*z)
                assert all((v%p if p else v)==0 for v in values)
                counts['independent_product_point_generator_annihilations']+=len(values)
            _,ev1=point_model(tuple((j,inv(j)) for j in range(1,m+1)),p)
            branch1=complete_planar_fan(m,ev1,p,counts)
            assert len(branch1)==m and planar_W(branch1)==divisor_increment(m)//3
            cells=fan(n,ev,p,counts)
            rays,walls,U,W=fan_budgets(cells)
            actual_L=actual_edge_length(cells)
            assert actual_L==b*m*(m-1)+m*r*(r-1)
            assert len(cells)==m*b
            predicted={}
            for (D1,normals1,*_),(D2,normals2,*_) in product(branch1,branch2):
                B=frozenset([(j,0,0) for j in range(r)]+
                            [(a,t+1,0) for a,t in D1]+[(a,0,t+1) for a,t in D2])
                poly=[tuple(Q(i==j) for i in range(3)) for j in range(3)]
                for q in (tuple((a,t,0) for a,t in normals1)+
                          tuple((a,0,t) for a,t in normals2)):
                    poly=clip(poly,q,1)
                assert polygon_key(cells[B])==polygon_key(poly)
                predicted[B]=poly
                counts['independent_entire_product_cone_equalities']+=1
            assert set(cells)==set(predicted)
            assert all(0 in q for q in walls.values())
            assert len(walls)==(m-1)*b+m*(b-1)
            assert W==b*(divisor_increment(m)//3)+m*c
            counts['independent_product_wall_inventories']+=1
            for B,poly in cells.items():
                G,_=marked_basis(B,ev,p)
                w=tuple(sum(a[i] for a in poly)/len(poly) for i in range(3))
                counts['independent_product_original_S_pairs']+=pair_checks(G,w,p)
            delta=None
            if previous and previous[0]==m-1:
                oldm,oldW,oldpoints,oldcells=previous
                assert set(points)-set(oldpoints)=={(m,inv(m),0)}
                denominator=1
                for j in range(1,m):denominator*=m-j
                norm=pow(denominator,-1,p) if p else Q(1,denominator)
                selector=[]
                for x,y,z in points:
                    value=m*y*norm
                    for j in range(1,m):value*=x-j
                    selector.append(value%p if p else value)
                target=points.index((m,inv(m),0))
                assert selector==[int(i==target) for i in range(n)]
                for axis in range(3):
                    assert all(((point[axis]-points[target][axis])*s)%p==0 if p else
                               (point[axis]-points[target][axis])*s==0
                               for point,s in zip(points,selector))
                    counts['independent_product_separator_eigenline_actions']+=1
                E=tuple(hyperbolic_cross(n,3))
                for D,poly in oldcells.items():
                    w=tuple(sum(a[i] for a in poly)/len(poly) for i in range(3))
                    B=frozenset(field_basis(sorted(E,key=lambda a:(dot(w,a),a)),ev,p)[0])
                    assert D<B
                    counts['independent_product_adjacent_basis_inclusions']+=1
                delta=W-oldW
                assert delta==b*(divisor_increment(m)-divisor_increment(m-1))//3+c
                counts['independent_product_one_point_W_increments']+=1
            rows.append({'kind':kind,'characteristic':p,'fixed_branch_length':r,
                         'varying_branch_length':m,'length':n,'initial_cones':len(cells),
                         'wall_pieces':len(walls),'W':W,'W_increment':delta,
                         'other_branch_cones':b,'other_branch_W':c,
                         'actual_edge_length':actual_L,
                         'actual_length_increment':None if delta is None else 2*b*(m-1)+r*(r-1),
                         'complete_standard_sets':[sorted(B) for B in sorted(cells,key=lambda B:sorted(B))]})
            previous=m,W,points,cells
            print({k:v for k,v in rows[-1].items() if k!='complete_standard_sets'},flush=True)
    negative=next(row for row in rows if row['fixed_branch_length']==12 and
                  row['varying_branch_length']==12)
    assert negative['W']==1824 and negative['W_increment']==-332
    assert negative['actual_edge_length']==3168 and negative['actual_length_increment']==396
    assert next(row['W'] for row in rows if row['fixed_branch_length']==12 and
                row['varying_branch_length']==11)==2156
    return {'counts':dict(counts),'actual_complete_product_fans':rows,
            'scope':'Actual radical point quotients, complete independent border traversal, every product cone and ordinary S-pairs. The length-36 to length-35 example independently checks a negative W increment; no asymptotic is inferred from these finite cases.'}


def wall_variation_controls():
    """Independent divisor, gcd, multiplicative and odd-square identities."""
    limit=4096;tot=list(range(limit+1));mu=[1]*(limit+1);mu[0]=0
    for ell in range(2,limit+1):
        if tot[ell]==ell:
            for n in range(ell,limit+1,ell):tot[n]-=tot[n]//ell;mu[n]*=-1
            for n in range(ell*ell,limit+1,ell*ell):mu[n]=0
    a=[-1]*(limit+1);a[0]=0
    for d in range(1,limit+1):
        for n in range(d,limit+1,d):a[n]+=d*tot[d]
    counts=Counter();primes=[]
    for n in range(1,limit+1):
        assert 0<=a[n]<=n*(n-1)
        rest=n;F=1
        for ell in range(2,isqrt(n)+1):
            if rest%ell:continue
            e=0
            while rest%ell==0:rest//=ell;e+=1
            F*=(ell**(2*e+1)+1)//(ell+1)
        if rest>1:F*=rest*rest-rest+1
        assert a[n]+1==F
        if n%2==0:assert 4*(a[n]+1)<=3*n*n
        if n<=128:
            assert a[n]==sum(n//gcd(j,n) for j in range(1,n))
            counts['independent_hook_gcd_costs']+=1
        if n>=3 and tot[n]==n-1:
            assert 4*(a[n]-a[n-1])>=(n+1)**2
            if n>=23 and n<limit:assert 8*(a[n+1]-a[n])<=-n*n
            primes.append(n);counts['independent_prime_neighbor_jump_bounds']+=1
        counts['independent_multiplicative_wall_costs']+=1
    for X in range(1,257):
        oddB=sum(d*tot[d] for d in range(1,X+1,2))
        def odd_square_sum(M):
            k=(M+1)//2
            return k*(4*k*k-1)//3
        assert oddB==sum(mu[e]*e*odd_square_sum(X//e) for e in range(1,X+1,2))
        assert sum(a[n] for n in range(1,X+1,2))==sum(
            sum(d*tot[d] for d in range(1,X//k+1,2)) for k in range(1,X+1,2))-(X+1)//2
        counts['independent_exact_odd_Mobius_and_divisor_replays']+=1
    chains=[]
    for r in (16,32,64,128,256,512,1024,2048,4096):
        Ws=[r*a[m]+m*a[r] for m in range(1,r+1)]
        deltas=[v-u for u,v in zip(Ws,Ws[1:])]
        positive=sum(max(v,0) for v in deltas);negative=sum(max(-v,0) for v in deltas)
        L=sum(a[m-1]-a[m] for m in range(2,r+1,2))
        odd=sum(a[m] for m in range(1,r+1,2));total=sum(a[1:r+1])
        assert L==2*odd-total
        assert negative>=r*L-(r//2)*a[r]
        assert positive-negative==Ws[-1]-Ws[0]>=0
        assert positive+negative==sum(map(abs,deltas))
        assert positive+negative<=2*r*total+(r-1)*a[r]
        chains.append({'fixed_branch_length':r,'maximum_length':3*r,
                       'parity_difference':L,'total_planar_W':total,'odd_planar_W':odd,
                       'positive_variation':positive,'negative_variation':negative,
                       'absolute_variation':positive+negative,'maximum_W':max(Ws),
                       'even_step_lower_bound':r*L-(r//2)*a[r]})
        counts['independent_exact_fixed_curve_chain_variations']+=1
    return {'counts':dict(counts),'chain_witnesses':chains,
            'prime_jump_parameters_checked':primes,
            'scope':'Exact finite arithmetic and chain identities only. The parity asymptotic constants and quartic variation orders are proved analytically in the note; no numerical fitting or prime-density assumption is used.'}


def weighted_socle_product_controls():
    """Explicit commuting algebra quotients, not arbitrary projected matroids."""
    from run_interleaving_research import minimal_excluded_border
    from run_socle_circuit_checks import multiplication_controls
    counts=Counter();rows=[]
    for s,kind,p in product((1,2,3),('positive','negative'),(0,2,3)):
        S=tuple((s-i,i*i) if kind=='positive' else (i,s*s-i*i) for i in range(s+1))
        D={(a,b) for x,y in S for a in range(x+1) for b in range(y+1)}
        if kind=='positive':
            m=len(D)-s
            assert m==s*(s+1)*(2*s+1)//6+1
        else:
            m=len(D)
            assert m==(s+1)*(s*s+1)-s*(s+1)*(2*s+1)//6
        L=m
        def model(child):
            if kind=='positive':
                B1=sorted(D-set(S)) if child else sorted((D-set(S))|{S[0]})
                planar_rows=[{a:1,S[0]:-1} for a in S[1:]]
                if child:planar_rows.append({S[0]:1})
            else:
                B1=sorted(D-{max(S)}) if child else sorted(D)
                planar_rows=[{a:1 for a in S}] if child else []
            ref=tuple([(a,0,0) for a in range(L)]+
                      [(a,b+1,0) for a,b in B1]+[(a,0,1) for a in range(L)])
            index={a:j for j,a in enumerate(ref)}
            def ev1(a):
                if a not in D:return {}
                if kind=='positive' and a in S:
                    return {} if child else {S[0]:1}
                if kind=='negative' and child and a==max(S):
                    return {b:(-1)%p if p else Q(-1) for b in S if b!=a}
                return {a:1}
            @lru_cache(None)
            def ev(a):
                x,y,z=a
                if y and z:return {}
                if y:return {index[(b,c+1,0)]:v for (b,c),v in ev1((x,y-1)).items()}
                if z:return {index[(x+z-1,0,1)]:1} if x+z-1<L else {}
                return {index[(x,0,0)]:1} if x<L else {}
            ground=({(a,0,0) for a in range(L)}|
                    {(a,b+1,0) for a,b in D}|
                    {(a,0,z) for z in range(1,L+1) for a in range(L-z+1)})
            A={'n':len(ref),'char':p,'basis':ref,'D':ground,'ev':ev}
            multiplication_controls(A,counts)
            lifted=[{(L,0,0):1},{(0,1,1):1},{(0,0,2):1,(1,0,1):-1}]
            lifted += [{(a,b+1,0):1} for a,b in minimal_excluded_border(tuple(D))]
            lifted += [{(a,b+1,0):v for (a,b),v in row.items()} for row in planar_rows]
            for f in lifted:
                replay=Counter()
                for a,c in f.items():
                    for j,v in ev(a).items():replay[j]+=c*v
                assert all((v%p if p else v)==0 for v in replay.values())
                counts['independent_weighted_product_generator_replays']+=1
            return A
        parent,child=model(False),model(True)
        assert (parent['n'],child['n'])==(3*m,3*m-1)
        kernel={(a,b+1,0):1 for a,b in (S if kind=='negative' else (S[0],))}
        def evaluate(f,A):
            replay=Counter()
            for a,c in f.items():
                for j,v in A['ev'](a).items():replay[j]+=c*v
            return {j:(v%p if p else v) for j,v in replay.items() if (v%p if p else v)}
        assert evaluate(kernel,parent) and not evaluate(kernel,child)
        for axis in range(3):
            shifted={tuple(a[i]+(i==axis) for i in range(3)):c for a,c in kernel.items()}
            assert not evaluate(shifted,parent)
            counts['independent_weighted_product_socle_kernel_actions']+=1
        fans=[]
        for A,is_child in ((parent,False),(child,True)):
            cells=fan(A['n'],A['ev'],p,counts)
            variable=(kind=='positive')!=is_child
            assert len(cells)==(2*(s+1) if variable else 2)
            predicted={}
            for i,side in product(range(s+1) if variable else (None,),range(2)):
                if variable:
                    B1=(D-set(S))|{S[i]} if kind=='positive' else D-{S[i]}
                else:B1=D-set(S) if kind=='positive' else D
                B2={(j,0) if side==0 else (0,j) for j in range(L)}
                B=frozenset([(j,0,0) for j in range(L)]+
                            [(a,b+1,0) for a,b in B1]+[(a,0,b+1) for a,b in B2])
                poly=[tuple(Q(j==k) for j in range(3)) for k in range(3)]
                if variable:
                    if i:poly=clip(poly,(1,-(2*i-1),0),1)
                    if i<s:poly=clip(poly,(-1,2*i+1,0),1)
                poly=clip(poly,(-1,0,1) if side==0 else (1,0,-1),1)
                assert polygon_key(cells[B])==polygon_key(poly)
                predicted[B]=poly
                G,_=marked_basis(B,A['ev'],p)
                w=tuple(sum(a[j] for a in poly)/len(poly) for j in range(3))
                counts['independent_weighted_product_original_S_pairs']+=pair_checks(G,w,p)
                counts['independent_weighted_product_entire_cone_equalities']+=1
            assert set(predicted)==set(cells)
            rays,walls,U,W=fan_budgets(cells)
            actual_L=actual_edge_length(cells)
            expected_L=2*s*(s+1)+(s+1)*L*(L-1) if variable else L*(L-1)
            assert actual_L==expected_L
            assert W==(2*(s+1)**2 if variable else 2)
            assert len(rays)==(s if variable else 0)
            visible=0
            for B,C in map(tuple,edges_from_polygons(cells,Counter())):
                if any(min(v)>0 for v in set(clean(cells[B]))&set(clean(cells[C]))):
                    visible+=sum(abs(a-b) for a,b in zip(profile(B),profile(C)))
            assert visible==(actual_L if variable else 0)
            counts['independent_visible_actual_skeletons']+=1
            sections=section_controls(cells,counts)
            fans.append((cells,actual_L))
        signed=fans[0][1]-fans[1][1]
        expected=s*m*(m-1)+2*s*(s+1)
        assert signed==(expected if kind=='positive' else -expected)
        E=tuple(hyperbolic_cross(3*m,3))
        for B,poly in fans[0][0].items():
            w=tuple(sum(a[j] for a in poly)/len(poly) for j in range(3))
            C=frozenset(field_basis(sorted(E,key=lambda a:(dot(w,a),a)),child['ev'],p)[0])
            assert C<B
            counts['independent_weighted_product_ordered_basis_deletions']+=1
        rows.append({'parameter':s,'kind':kind,'characteristic':p,'branch_length':m,
                     'lengths':[3*m,3*m-1],'tips':S,'downset_size':len(D),
                     'fan_counts':[len(F) for F,_ in fans],
                     'actual_edge_lengths':[v for _,v in fans],'signed_increment':signed,
                     'last_quotient_section_witnesses':sections,
                     'kernel_polynomial':sorted(kernel.items())})
        print(rows[-1],flush=True)
    guards=[]
    for k in (2,3,4,7,12):
        @lru_cache(None)
        def ev(a):
            x,y,z=a
            if x>=2:return {}
            return {(x+j)*k+y+z-j:comb(z,j) for j in range(min(z,1-x)+1) if y+z-j<k}
        cells=fan(2*k,ev,0,counts)
        actual_L=actual_edge_length(cells)
        assert len(cells)==5 and actual_L==(2*k)**2-2*k+2
        guards.append({'parameter':k,'length':2*k,'fan_count':5,'actual_L':actual_L,'W':10})
    return {'counts':dict(counts),'all_field_product_steps':rows,
            'quadratic_primitive_to_actual_loss':guards,
            'scope':'All-field formulas proved in the note; small characteristic-zero, two and three commuting-algebra fans independently checked. A growing-fan multiplier, rather than primitive gcd normalization, gives positive and negative order n^(7/3) actual-length jumps. No global cubic budget or count theorem follows.'}


def two_quadratic_actual_budget_controls():
    """Sharp bounded-box cost and existing grid duality on original algebras."""
    from run_parity_point_fan_checks import predicted_polygons,paired_basis,inventory
    counts=Counter();rows=[]
    def check_bound(cells):
        n=len(next(iter(cells)));edges=edges_from_polygons(cells,Counter())
        assert all(a[0]<=1 and a[1]<=1 for B in cells for a in B)
        projected={profile(B)[:2] for B in cells}
        assert len(projected)==len(cells)<=(n+1)**2
        lengths=[];multiples=[]
        for B,C in map(tuple,edges):
            a,b=profile(B),profile(C);delta=tuple(y-x for x,y in zip(a,b))
            k=gcd(*delta);q=tuple(x//k for x in delta)
            assert 0<max(abs(q[0]),abs(q[1]))<=1 and abs(q[2])<=n-1
            assert sum(map(abs,q))<=n+1
            for j in range(1,k):
                xy=tuple(a[i]+j*q[i] for i in range(2))
                assert xy not in projected;projected.add(xy)
            lengths.append(sum(map(abs,delta)));multiples.append(k)
        assert len(projected)<=(n+1)**2 and sum(multiples)<=3*(n+1)**2
        assert sum(lengths)==actual_edge_length(cells)<=3*(n+1)**3
        counts['independent_two_quadratic_actual_budget_checks']+=1
        counts['independent_injective_lower_fiber_lattice_points']+=len(projected)
        return max(multiples,default=0)
    for s,added in product(range(1,9),(False,True)):
        D=2*s-1;box=tuple(product((0,1),(0,1),range(D+1)))
        original={paired_basis(box,tuple(sum(a[i] for a in P)/len(P) for i in range(3)),added):P
                  for P in predicted_polygons(D,added)}
        assert check_bound(original)==1
        W=inventory(original)[2];L=actual_edge_length(original)
        assert W==L
        counts['independent_complete_parity_profile_unit_multiples']+=len(edges_from_polygons(original,Counter()))
        if s<=3:
            v=(1,1,D);G=set(box)
            dual={frozenset(G-{tuple(v[i]-a[i] for i in range(3)) for a in B}):P for B,P in original.items()}
            original_profiles={B:profile(B) for B in original}
            box_moment=profile(box)
            for B,C in zip(original,dual):
                expected=tuple(box_moment[i]-len(B)*v[i]+original_profiles[B][i] for i in range(3))
                assert profile(C)==expected
            for p in (0,101):
                nodes=tuple(a for a in box if sum(a)%2==0 and (not added or any(a)))
                n,ev=point_model(nodes,p);assert n==4*s-int(added)
                actual=fan(n,ev,p,counts)
                assert set(actual)==set(dual)
                for B,P in actual.items():
                    assert polygon_key(P)==polygon_key(dual[B])
                    Gbasis,_=marked_basis(B,ev,p)
                    w=tuple(sum(a[i] for a in P)/len(P) for i in range(3))
                    counts['independent_complement_point_original_S_pairs']+=pair_checks(Gbasis,w,p)
                    counts['independent_complement_entire_fan_cones']+=1
                assert actual_edge_length(actual)==L and inventory(actual)[2]==W
                check_bound(actual)
                rows.append({'parameter':s,'characteristic':p,'complement_of_point_added':added,
                             'length':n,'fan_count':len(actual),'actual_L':L,'W':W,
                             'original_profile_translation':tuple(box_moment[i]-(4*s+int(added))*v[i] for i in range(3))})
    nilpotent=[]
    for m,p in product((1,2,4),(0,2,3)):
        @lru_cache(None)
        def ev(a):
            x,y,z=a
            return {(x+y)*m+z:1} if x+y<2 and z<m else {}
        cells=fan(2*m,ev,p,counts)
        assert len(cells)==2 and actual_edge_length(cells)==2*m
        assert check_bound(cells)==m
        nilpotent.append({'relation_model':'x-y,x^2,z^m','parameter':m,
                           'characteristic':p,'length':2*m,'actual_L':2*m,'W':2,
                           'edge_lattice_multiple':m})
    return {'counts':dict(counts),'complete_complement_point_fans':rows,
            'nonprimitive_nilpotent_edge_controls':nilpotent,
            'scope':'The global O(n^3) actual-length bound is a box-bound corollary proved in the note, not inferred from checks. The existing parity classification gives primitive actual edges and cubic sharpness; grid complement duality gives the opposite signed one-point update. Original point algebras, rather than only paired models, independently verify the complement fans. No unrestricted count or novelty claim.'}


def main():
    if not __debug__:
        raise RuntimeError('Exact assertions must be enabled.')
    proof_sources = ['docs/ideal_step_incidence.md', 'docs/weighted_facet_incidence.md',
                     'docs/wall_cost_products.md','docs/adjacent_rank_fans.md',
                     'docs/coordinate_slice_wall_budget.md','docs/parity_point_fan_updates.md',
                     'docs/interleaving_research.md']
    before = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT/'scripts').glob('*.py')}
    before.update({f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in proof_sources})
    counts = Counter()
    local = local_controls(counts)
    actual = actual_controls(counts)
    arithmetic = arithmetic_controls(counts)
    product_steps = product_wall_controls()
    variations = wall_variation_controls()
    weighted_steps = weighted_socle_product_controls()
    bounded_pair = two_quadratic_actual_budget_controls()
    loaded = {Path(m.__file__).resolve() for m in tuple(sys.modules.values())
              if getattr(m, '__file__', None) and
              Path(m.__file__).resolve().parent == ROOT/'scripts'}
    loaded.add(Path(__file__).resolve())
    sources = proof_sources + sorted(str(p.relative_to(ROOT)) for p in loaded)
    assert all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest() == before[f] for f in sources)
    report = {'all_checks_passed': True,
              'scope': 'Actual nested radical point ideals and complete exact border fans; finite residue distributions, totient recurrences and larger integer sums. These support, and do not computationally prove, the uniform weighted asymptotic or the infinite cubic-jump theorem. No unrestricted budget, extremal solution or priority certificate is claimed.',
              'actual_point_ideals': actual, 'local_residue_controls': local,
              'integer_incidence_sums': arithmetic, 'counts': dict(counts),
              'product_wall_length_ideal_steps':product_steps,
              'wall_length_chain_variation':variations,
              'actual_edge_length_steps':weighted_steps,
              'two_quadratic_actual_global_budget':bounded_pair,
              'uniform_wall_length_conclusions':{
                  'one_point_jumps':'Both signs have order n^3 on a single fixed reduced curve, or n^3*log(n+1) with a generic planar product branch. The formula retains every actual wall piece.',
                  'deletion_chain_variation':'For even r tending to infinity, actual consecutive deletions of lengths between 2r+1 and 3r have positive, negative and absolute W variation Theta(r^4), or Theta(r^4*log(r+1)) for the generic branch. Every individual W stays respectively O(r^3) or O(r^3*log(r+1)).',
                  'parity_average_constant':'sum_odd a_m - sum_even a_m = zeta(3)/(18*zeta(2))*X^3 + O(X^2*log(X+1)) for even X; derived by an explicit odd-square Mobius sum, without a prime-density hypothesis.',
                  'uniform_quadratic_polylog_increment_refuted':True,
                  'uniform_cubic_polylog_positive_variation_over_every_chain_refuted':True,
                  'existence_of_low_variation_chain_refuted':False,
                  'global_cubic_polylog_W_refuted':False,
                  'unrestricted_G3_improved':False,'novelty_established':False,'goal_complete':False},
              'loaded_source_stability_checked': True,
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    target = ROOT/'results/research/ideal_step_incidence_checks.json'
    target.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'all_checks_passed': True, 'counts': dict(counts),
                      'report': str(target.relative_to(ROOT))}), flush=True)


if __name__ == '__main__':
    main()
