#!/usr/bin/env python3
"""Exact adjacent point-ideal fans and arithmetic incidence controls."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
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
                         'complete_standard_sets':[sorted(B) for B in sorted(cells,key=lambda B:sorted(B))]})
            previous=m,W,points,cells
            print({k:v for k,v in rows[-1].items() if k!='complete_standard_sets'},flush=True)
    negative=next(row for row in rows if row['fixed_branch_length']==12 and
                  row['varying_branch_length']==12)
    assert negative['W']==1824 and negative['W_increment']==-332
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


def main():
    if not __debug__:
        raise RuntimeError('Exact assertions must be enabled.')
    proof_sources = ['docs/ideal_step_incidence.md', 'docs/weighted_facet_incidence.md',
                     'docs/wall_cost_products.md']
    before = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT/'scripts').glob('*.py')}
    before.update({f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in proof_sources})
    counts = Counter()
    local = local_controls(counts)
    actual = actual_controls(counts)
    arithmetic = arithmetic_controls(counts)
    product_steps = product_wall_controls()
    variations = wall_variation_controls()
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
