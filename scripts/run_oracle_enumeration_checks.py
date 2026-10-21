#!/usr/bin/env python3
"""Exact bounded verification of truncated-state oracle enumeration.

The hull routine below intentionally enumerates triples. The complexity
theorem uses the cited three-dimensional hull algorithm, not this checker.
"""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd, lcm
from pathlib import Path
import hashlib
import json

from run_higher_interaction_checks import rank_trace
from run_interleaving_research import minimal_excluded_border
from run_liaison_fan_checks import rectangle_graph
from run_planar_dyadic_bound_checks import convex_hull
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import area2, clip, cross, determinant, dot, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]


def apply(M, v, p):
    out = tuple(sum(c*column[i] for column, c in zip(M, v)) for i in range(len(v)))
    return tuple(c % p for c in out) if p else out


def image_from_matrices(matrices, one, a, p):
    v = one
    for M, e in zip(matrices, a):
        for _ in range(e):
            v = apply(M, v, p)
    return v


def curve_model(n):
    polynomials = ((0, 1), (0, 1, 1), (0, 1, 2))

    @lru_cache(None)
    def ev(a):
        v = [1]+[0]*(n-1)
        for f, e in zip(polynomials, a):
            for _ in range(e):
                out = [0]*n
                for i, c in enumerate(v):
                    for j, d in enumerate(f):
                        if i+j < n:
                            out[i+j] += c*d
                v = out
        return {i: c for i, c in enumerate(v) if c}
    reference = tuple((i, 0, 0) for i in range(n))
    matrices = tuple(tuple(tuple(ev(tuple(a+(i == j) for j, a in enumerate(v))).get(k, 0)
                                       for k in range(n)) for v in reference) for i in range(3))
    return n, ev, matrices, (1,)+(0,)*(n-1)


def rectangle_fixture(F):
    B, ev, _ = rectangle_graph(2, 2, F)
    n = len(B)
    matrices = tuple(tuple(tuple(ev(tuple(a+(i == j) for j, a in enumerate(v))).get(k, 0)
                                       for k in range(n)) for v in B) for i in range(3))
    return n, ev, matrices, (1,)+(0,)*(n-1)


def point_model():
    points = tuple((t, t*t+1, t*t+t+1) for t in (1, 2, 3))
    def ev(a):
        return {i: x**a[0]*y**a[1]*z**a[2] for i, (x, y, z) in enumerate(points)}
    matrices = tuple(tuple(tuple(points[j][i] if j == k else 0 for k in range(3))
                           for j in range(3)) for i in range(3))
    return 3, ev, matrices, (1, 1, 1)


def hull_facets(points):
    """Independent small exact H-representation, merging coplanar triangles."""
    P = tuple(sorted(points))
    out = {}
    for a, b, c in combinations(P, 3):
        z = cross(tuple(x-y for x, y in zip(b, a)), tuple(x-y for x, y in zip(c, a)))
        if not any(z):
            continue
        h = dot(z, a)
        values = [dot(z, v)-h for v in P]
        if min(values) < 0 < max(values) or not any(values):
            continue
        if max(values) > 0:
            z, h = tuple(-x for x in z), -h
        d = gcd(*z)
        z, h = tuple(x//d for x in z), h//d
        out[z, h] = frozenset(v for v in P if dot(z, v) == h)
    return out


def compressed_krylov(matrices, seeds, key, p, counts):
    """Order-faithful small implementation of the published doubling method.

    Field arithmetic here is elementary; the stated complexity uses the
    published fast rank-profile and matrix-multiplication algorithms.
    Labels are (exponent tuple, seed index), avoiding a full grid index.
    """
    r=len(matrices);n=len(seeds[0]);zero=(0,)*r
    rows={(zero,j):tuple(c%p if p else Q(c) for c in v) for j,v in enumerate(seeds)}
    def compress(data):
        labels=sorted(data,key=key)
        selected,_=field_basis(labels,lambda u:{i:c for i,c in enumerate(data[u]) if c},p)
        return {u:data[u] for u in selected}
    rows=compress(rows)
    assert len(rows)<=n
    for axis,M in enumerate(matrices):
        power=M;step=1;rounds=0
        while True:
            old=tuple(rows);candidate=dict(rows)
            for (a,j),vector in rows.items():
                shifted=tuple(x+(step if i==axis else 0) for i,x in enumerate(a))
                u=(shifted,j);value=apply(power,vector,p)
                if u in candidate:
                    assert candidate[u]==value
                candidate[u]=value
            assert len(candidate)<=2*n
            updated=compress(candidate)
            counts['compressed_rank_profile_stages']+=1
            counts['compressed_candidate_rows']+=len(candidate)
            if len(updated)==len(rows) and tuple(updated)!=old:
                counts['equal_rank_label_changes']+=1
            rows=updated;rounds+=1
            if tuple(rows)==old:
                break
            power=tuple(apply(power,column,p) for column in power)
            step*=2
            assert rounds<=n.bit_length()+1
        assert rounds<=n.bit_length()+1
    counts['compressed_basis_calls']+=1
    return rows


def fixed_module(model,p):
    n,ev,matrices,one=model
    chains=tuple((b,c) for b in range(n) for c in range(n//(b+1)))
    seeds=tuple(image_from_matrices(matrices,one,(0,b,c),p) for b,c in chains)
    assert len(chains)==sum(n//j for j in range(1,n+1))
    return chains,seeds


def module_standard_basis(model,p,key,counts,module=None):
    n,ev,matrices,one=model
    chains,seeds=module or fixed_module(model,p)
    def order(label):
        (a,),j=label;b,c=chains[j]
        return key((a,b,c))
    rows=compressed_krylov((matrices[0],),seeds,order,p,counts)
    result={(a[0],*chains[j]) for a,j in rows}
    assert len(result)==n
    for (a,j),vector in rows.items():
        exponent=(a[0],*chains[j])
        assert vector==image_from_matrices(matrices,one,exponent,p)
        counts['compressed_direct_column_replays']+=1
    for j in range(len(chains)):
        exponents={a[0] for a,i in rows if i==j}
        assert exponents==set(range(len(exponents)))
    return frozenset(result)


def prefix_module_controls(cases):
    counts=defaultdict(int);fixtures=[]
    for label,model,p in cases:
        n,ev,matrices,one=model;chains,seeds=fixed_module(model,p)
        rows=[]
        for w in ((1,1,1),(2,3,5),(7,1,2),(1,7,2)):
            for threshold in range(-1,2*n+2):
                lengths=[]
                for b,c in chains:
                    limit=n//((b+1)*(c+1))
                    cost=threshold-w[1]*b-w[2]*c
                    lengths.append(max(0,min(limit,cost//w[0]+1)))
                def order(u):
                    (a,),j=u
                    return (a-lengths[j],j)
                chosen=compressed_krylov((matrices[0],),seeds,order,p,counts)
                predicted=sum(a[0]<lengths[j] for a,j in chosen)
                ground=[a for a in hyperbolic_cross(n,3) if dot(w,a)<=threshold]
                direct,_=rank_trace(ground,ev,p,set())
                assert predicted==len(direct)
                counts['independent_nongeneric_prefix_rank_queries']+=1
                rows.append({'weight':w,'threshold':threshold,'rank':predicted})
        fixtures.append({'model':label,'length':n,'characteristic':p,
                         'fixed_chain_generators':len(chains),'prefix_queries':rows})
    return {'counts':dict(counts),'fixtures':fixtures,
            'scope':'One fixed univariate module per actual commuting quotient, with query shifts from exact chain lengths. This is faithful for special and nonreduced samples, and is not the generic min(n,carrier-rank) formula.'}


def digit_rank(n, modulus, p, residues):
    """Rank of selected columns binom(residue,j), 0 <= j < n, in F_p.

    The modulus is a power of the supplied prime p. This computes ranks,
    rather than inferring them from the nonzero-entry pattern.
    """
    assert 0 <= n <= modulus
    residues = set(residues)
    assert all(0 <= u < modulus for u in residues)
    if n == 0 or not residues:
        return 0
    if n == modulus:
        return len(residues)
    assert modulus % p == 0
    h = modulus//p
    q, r = divmod(n, h)
    fibers = defaultdict(int)
    for u in residues:
        fibers[u % h] += 1
    return sum(min(q, size) for size in fibers.values()) + digit_rank(
        r, h, p, {u for u, size in fibers.items() if size > q})


def digit_basis_count(n, modulus, p):
    if n == 0 or n == modulus:
        return 1
    h = modulus//p
    q, r = divmod(n, h)
    return (digit_basis_count(r, h, p)*comb(p, q)**(h-r)
            *comb(p, q+1)**r)


def digit_rank_enumerator(n, modulus, p, absent, present, rank_weight):
    """Evaluate sum_S absent^(M-|S|) present^|S| rank_weight^rank(S)."""
    if n == 0:
        return (absent+present)**modulus
    if n == modulus:
        return (absent+present*rank_weight)**modulus
    h = modulus//p
    q, r = divmod(n, h)
    low = sum(comb(p, l)*absent**(p-l)*present**l*rank_weight**l
              for l in range(q+1))
    high = sum(comb(p, l)*absent**(p-l)*present**l*rank_weight**q
               for l in range(q+1, p+1))
    return digit_rank_enumerator(r, h, p, low, high, rank_weight)


def digit_code_enumerator(n, modulus, p, field_size, absent, present):
    """Division-free specialization to the repeated-root Hamming enumerator."""
    if n == 0:
        return (absent+field_size*present)**modulus
    if n == modulus:
        return (absent+present)**modulus
    h = modulus//p
    q, r = divmod(n, h)
    low = sum(comb(p, l)*absent**(p-l)*present**l for l in range(q+1))
    high = sum(comb(p, l)*absent**(p-l)*field_size**(l-q-1)*present**l
               for l in range(q+1, p+1))
    return digit_code_enumerator(r, h, p, field_size, low, high)


def cheapest_digit_basis(n, modulus, p, costs):
    """Minimum-cost basis when all residue classes are available."""
    assert len(costs) == modulus and len(set(costs)) == modulus
    if n == 0:
        return frozenset()
    if n == modulus:
        return frozenset(range(modulus))
    h = modulus//p
    q, r = divmod(n, h)
    choices = [sorted((u+v*h for v in range(p)), key=costs.__getitem__)
               for u in range(h)]
    selected = {a for fiber in choices for a in fiber[:q]}
    extra_costs = [costs[fiber[q]] for fiber in choices]
    extras = cheapest_digit_basis(r, h, p, extra_costs)
    selected.update(choices[u][q] for u in extras)
    assert len(selected) == n
    return frozenset(selected)


def integer_binomial(a, j):
    value = 1
    for i in range(j):
        value = value*(a-i)//(i+1)
    return value


def unit_power_model(n, k):
    """Actual quotient ((X-1)^n,XY-1,Z-Y^k), in its t=X-1 basis."""
    images = ((1, 1)+(0,)*max(0, n-2),
              tuple((-1)**j for j in range(n)),
              tuple((-1)**j*comb(k+j-1, j) for j in range(n)))
    matrices = tuple(tuple(tuple(f[i-j] if i >= j else 0 for i in range(n))
                           for j in range(n)) for f in images)

    @lru_cache(None)
    def ev(a):
        exponent = a[0]-a[1]-k*a[2]
        return {j: c for j in range(n) if (c := integer_binomial(exponent, j))}
    return n, ev, matrices, (1,)+(0,)*(n-1)


def prime_digit_controls():
    counts = defaultdict(int)
    base_counts = []
    for p, modulus in ((2, 2), (2, 4), (2, 8), (3, 3), (3, 9), (5, 5)):
        for n in range(modulus+1):
            bases = 0
            rank_sizes = defaultdict(int)
            for mask in range(1 << modulus):
                S = {a for a in range(modulus) if mask >> a & 1}
                def column(a):
                    return {j: comb(a, j) % p for j in range(n) if comb(a, j) % p}
                direct, _ = rank_trace(sorted(S), column, p, set())
                rank = digit_rank(n, modulus, p, S)
                assert rank == len(direct)
                counts['exhaustive_selected_pascal_ranks'] += 1
                bases += rank == n and len(S) == n
                rank_sizes[len(S), rank] += 1
            assert bases == digit_basis_count(n, modulus, p)
            for absent, present, rank_weight in ((1, 1, 1), (2, 3, 5), (3, 2, 0)):
                direct_value = sum(count*absent**(modulus-size)*present**size
                                   *rank_weight**rank
                                   for (size, rank), count in rank_sizes.items())
                assert digit_rank_enumerator(n, modulus, p, absent, present,
                                             rank_weight) == direct_value
                counts['exhaustive_subset_rank_enumerator_equalities'] += 1
            base_counts.append({'prime': p, 'modulus': modulus, 'rank': n,
                                'all_residue_bases': bases})
    for p, modulus in ((2, 8), (3, 9), (2, 32), (3, 27), (5, 25), (7, 49)):
        for n in sorted({1, modulus//3, modulus//2, modulus-1, modulus}):
            for slope in (1, 2, 5):
                costs = [((slope*a*a+3*a+7) % modulus)*modulus+a
                         for a in range(modulus)]
                def column(a):
                    return {j: comb(a, j) % p for j in range(n) if comb(a, j) % p}
                order = sorted(range(modulus), key=costs.__getitem__)
                direct, _ = rank_trace(order, column, p, set())
                predicted = cheapest_digit_basis(n, modulus, p, costs)
                assert predicted == frozenset(direct)
                counts['independent_minimum_cost_digit_bases'] += 1
                for S in (set(order[::2]), set(order[:modulus//3]),
                          {a for a in range(modulus) if a % p == 0}):
                    selected, _ = rank_trace(sorted(S), column, p, set())
                    assert digit_rank(n, modulus, p, S) == len(selected)
                    counts['larger_sparse_pascal_ranks'] += 1
    code_witnesses = []
    for p, modulus, n in ((2, 8, 3), (2, 8, 5), (3, 9, 4), (3, 9, 5), (5, 5, 2)):
        weights = [0]*(modulus+1)
        generator = [(-1)**(n-j)*comb(n, j) for j in range(n+1)]
        for coefficients in product(range(p), repeat=modulus-n):
            word = [0]*modulus
            for i, c in enumerate(coefficients):
                for j, d in enumerate(generator):
                    word[i+j] = (word[i+j]+c*d) % p
            assert all(sum(c*comb(a, j) for a, c in enumerate(word)) % p == 0
                       for j in range(n))
            weights[sum(c != 0 for c in word)] += 1
            counts['independent_polynomial_codeword_annihilations'] += 1
        assert sum(weights) == p**(modulus-n) and weights[0] == 1
        for z in (0, 1, 2, 3):
            predicted = digit_rank_enumerator(n, modulus, p, 1-z, p*z, Q(1, p))
            assert predicted == sum(count*z**weight for weight, count in enumerate(weights))
            assert digit_code_enumerator(n, modulus, p, p, 1-z, z) == predicted
            counts['independent_hamming_enumerator_equalities'] += 1
        code_witnesses.append({'prime': p, 'modulus': modulus,
                               'generator_multiplicity': n, 'weight_enumerator': weights})
    witnesses = []
    for p in (2, 3, 5):
        for n in (3, 5, 7):
            modulus = 1
            while modulus < n:
                modulus *= p
            for k in (2, 4):
                model = unit_power_model(n, k)
                _, ev, matrices, one = model
                E = hyperbolic_cross(n, 3)
                module = fixed_module(model, p)
                for a in E:
                    actual = image_from_matrices(matrices, one, a, p)
                    assert actual == tuple(ev(a).get(j, 0) % p for j in range(n))
                    counts['independent_unit_power_column_replays'] += 1
                for w in ((1, 1, 1), (2, 3, 5), (7, 1, 2), (1, 7, 2)):
                    ordered = sorted(E, key=lambda a: (dot(w, a), a))
                    direct, _ = rank_trace(ordered, ev, p, set())
                    compressed = module_standard_basis(
                        model, p, lambda a: (dot(w, a), a), counts, module)
                    assert compressed == frozenset(direct)
                    counts['unit_power_compressed_basis_equalities'] += 1
                    for threshold in range(-1, 2*n+2):
                        ground = [a for a in E if dot(w, a) <= threshold]
                        selected, _ = rank_trace(ground, ev, p, set())
                        S = {(a-b-k*c) % modulus for a, b, c in ground}
                        assert digit_rank(n, modulus, p, S) == len(selected)
                        counts['actual_unit_power_prefix_ranks'] += 1
                if n == 3 and k == 2:
                    S = {0, 2}
                    witnesses.append({'prime': p, 'length': n, 'modulus': modulus,
                                      'monomials': [(0, 0, 0), (2, 0, 0)],
                                      'residues': sorted(S),
                                      'rank': digit_rank(n, modulus, p, S)})
    # In characteristic two, distinct residue columns can also be
    # dependent. Distinctness alone is not the generic independence rule.
    assert digit_rank(3, 4, 2, {0, 2}) == 2
    assert digit_rank(3, 4, 2, {0}) == 1
    witnesses.append({'prime': 2, 'length': 3, 'modulus': 4,
                      'monomials': [(0, 0, 0), (4, 0, 0)],
                      'integer_powers': [0, 4], 'distinct_residues': [0], 'rank': 1})
    assert digit_rank(5, 8, 2, {0, 2, 4, 6}) == 3
    witnesses.append({'prime': 2, 'length': 5, 'modulus': 8,
                      'monomials': [(a, 0, 0) for a in (0, 2, 4, 6)],
                      'distinct_residues': [0, 2, 4, 6], 'rank': 3})
    return {'counts': dict(counts), 'exhaustive_basis_counts': base_counts,
            'exact_witnesses': witnesses,
            'independent_code_weight_enumerators': code_witnesses,
            'uniform_results': {
                'actual_quotient': 'K[X,Y,Z]/((X-1)^n,XY-1,Z-Y^k), k>=1, characteristic p prime',
                'rank': 'R(n,M,S)=sum_u min(q,l_u)+R(r,M/p,{u:l_u>q}), n=q*(M/p)+r, M a power of p >= n',
                'basis_count': 'B(n,M)=B(r,M/p)*binom(p,q)^(M/p-r)*binom(p,q+1)^r',
                'homogeneous_subset_enumerator': 'H(n,M,A,B,Y)=H(r,M/p,A_prime,B_prime,Y), with the two explicit binomial sums in (O16)',
                'code_weight_enumerator': 'W_C(z)=H(n,M,1-z,Q*z,1/Q) for C=<(V-1)^n> over a finite field of size Q and characteristic p',
                'division_free_circuit': 'O(p*s+log(M+1)) integer-constant arithmetic gates for the complete enumerator, M=p^s; coefficient expansion and integer bit costs are additional',
                'coherent_order_count_established': False},
            'scope': 'Exact rank, minimum-cost basis, subset-rank and code-weight recursions for the truncated Pascal matrix over prime characteristic. Independent controls use arbitrary column subsets, actual unit-power quotients and polynomial-generated codewords. This is repeated-root cyclic-code structure; strongest priority and the comparison with existing weight-distribution algorithms remain open. No novelty, coherent-fan lower bound or improvement of unrestricted G3 is asserted.'}



def monomial_key(w, a):
    M = max(w)
    if M <= 0:
        return (-dot(w, a),)+a
    j = max(i for i in range(3) if w[i] == M)
    d = sum(a)
    return (M*d-dot(w, a),)+tuple(a[i]-d*(i == j) for i in range(3))+(a[j],)


def exhaustive_bases(E, n, ev, p, counts):
    profiles = defaultdict(list)
    for B in combinations(E, n):
        d = determinant([[ev(a).get(i, 0) for a in B] for i in range(n)])
        counts['independent_basis_determinants'] += 1
        if (d % p if p else d):
            mu = tuple(sum(a[i] for a in B) for i in range(3))
            profiles[mu].append(frozenset(B))
    assert profiles
    positive = {}
    for b in profiles:
        polygon = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for c in profiles:
            z = tuple(x-y for x, y in zip(c, b))
            if any(z):
                polygon = clip(polygon, z, 1)
            if not polygon:
                break
        if polygon:
            assert len(profiles[b]) == 1
            weight = tuple(sum(v[i] for v in polygon)/len(polygon) for i in range(3))
            positive[b] = profiles[b][0], polygon, weight
    assert sum(area2(poly) for _, poly, _ in positive.values()) == 1
    counts['independent_positive_state_polytopes'] += 1
    return profiles, positive


def verify_reduced_basis(B, ev, matrices, weight, p, counts):
    chosen, pivots = field_basis(sorted(B), ev, p)
    assert len(chosen) == len(B)
    assert all(tuple(x-(i == j) for j, x in enumerate(a)) in B
               for a in B for i in range(3) if a[i])
    records = []
    for a in sorted(minimal_excluded_border(tuple(B))):
        rem, co = field_reduce(ev(a), pivots, p)
        assert not rem
        assert all(dot(a, weight) > dot(b, weight) for b in co)
        polynomial = {a: Q(1)}
        polynomial.update({b: (-c % p if p else -c) for b, c in co.items()})
        for i in range(len(B)):
            vector = tuple(int(k == i) for k in range(len(B)))
            replay = [0]*len(B)
            for v, c in polynomial.items():
                for k, d in enumerate(image_from_matrices(matrices, vector, v, p)):
                    replay[k] += c*d
            assert all(not (c % p if p else c) for c in replay)
            counts['independent_polynomial_matrix_annihilations'] += 1
        counts['recovered_reduced_basis_polynomials'] += 1
        records.append({'head': a, 'tails': [{'exponent': b, 'coefficient': str(c)}
                                            for b, c in sorted(co.items())]})
    return records


def polar_cone_polygon(normals, strict_weight):
    """Exact planar intersection by polarity and a sorted monotone-chain hull."""
    total = sum(strict_weight)
    center = tuple(Q(x, total) for x in strict_weight)
    constraints = tuple(normals)+((1, 0, 0), (0, 1, 0), (0, 0, 1))
    dual = set()
    for z in constraints:
        slack = dot(z, center)
        assert slack > 0
        dual.add((Q(z[2]-z[0], slack), Q(z[2]-z[1], slack)))
    hull = convex_hull(dual)
    assert len(hull) >= 3
    polygon = []
    for a, b in zip(hull, hull[1:]+hull[:1]):
        determinant2 = a[0]*b[1]-a[1]*b[0]
        assert determinant2
        x = center[0]+(b[1]-a[1])/determinant2
        y = center[1]+(a[0]-b[0])/determinant2
        v = (x, y, 1-x-y)
        assert all(dot(z, v) >= 0 for z in constraints)
        polygon.append(v)
    assert area2(polygon) > 0
    return polygon


def exact_adjacent_weight(n, a, b, normal):
    midpoint = tuple((x+y)/2 for x, y in zip(a, b))
    denominator = lcm(*(x.denominator for x in midpoint))
    step = tuple(sum(normal)-3*x for x in normal)
    assert dot(normal, step) < 0 and sum(step) == 0
    scale = 2*(1+18*n*n)
    raw = tuple(scale*denominator*x+y for x, y in zip(midpoint, step))
    assert all(x.denominator == 1 for x in raw) and min(raw) > 0
    w = tuple(int(x) for x in raw)
    factor = gcd(*w)
    return tuple(x//factor for x in w), midpoint


def adjacency_walk_controls(cases, references):
    counts = defaultdict(int)
    fixtures = []
    for (label, model, p), reference in zip(cases, references):
        n, ev, matrices, one = model
        E = tuple(sorted(hyperbolic_cross(n, 3)))
        known = {frozenset(map(tuple, c['standard_set'])): c
                 for c in reference['all_initial_ideals_and_reduced_bases']}
        queue = [((1, n+1, (n+1)**2), None, None)]
        visited, graph, trace = {}, set(), []
        queries = 0
        while queue:
            weight, previous, edge = queue.pop()
            data = compressed_krylov(matrices, (one,),
                                     lambda u: (dot(weight, u[0]),)+u[0], p, counts)
            B = frozenset(a for a, j in data)
            direct, _ = rank_trace(sorted(E, key=lambda a: (dot(weight, a), a)), ev, p, set())
            assert B == frozenset(direct)
            queries += 1
            counts['direct_walk_basis_queries'] += 1
            assert B in known
            if previous is not None:
                assert B != previous
                graph.add(frozenset((B, previous)))
                records = known[B]['reduced_basis']
                normals = [tuple(x-y for x, y in zip(row['head'], t['exponent']))
                           for row in records for t in row['tails']]
                assert all(dot(z, weight) > 0 for z in normals)
                assert all(dot(z, v) >= 0 for z in normals for v in edge)
                counts['independent_whole_shared_edge_neighbors'] += 1
                trace.append({'from_profile': tuple(sum(a[i] for a in previous) for i in range(3)),
                              'to_profile': tuple(sum(a[i] for a in B) for i in range(3)),
                              'neighbor_weight': weight,
                              'shared_edge': [[str(x) for x in v] for v in edge]})
            if B in visited:
                continue
            records = verify_reduced_basis(B, ev, matrices, weight, p, counts)
            normals = tuple(sorted({tuple(x-y for x, y in zip(row['head'], t['exponent']))
                                    for row in records for t in row['tails']}))
            assert all(max(map(abs, z)) <= n for z in normals)
            assert len(records)**3 <= 343*n*n
            polygon = polar_cone_polygon(normals, weight)
            clipped = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
            for z in normals:
                clipped = clip(clipped, z, 1)
            expected = [tuple(Q(x) for x in v) for v in known[B]['weight_simplex_polygon']]
            assert area2(polygon) == area2(clipped) == area2(expected)
            assert all(dot(z, v) >= 0 for z in normals for v in expected)
            profile = tuple(sum(a[i] for a in B) for i in range(3))
            visited[B] = {'profile': profile, 'standard_set': sorted(B),
                          'strict_weight': weight, 'inequality_count': len(normals),
                          'polygon': [[str(x) for x in v] for v in polygon]}
            counts['polar_hull_clipping_and_state_polytope_equalities'] += 1
            for a, b in zip(polygon, polygon[1:]+polygon[:1]):
                if any(a[i] == b[i] == 0 for i in range(3)):
                    continue
                z = next(z for z in normals if dot(z, a) == dot(z, b) == 0)
                adjacent, midpoint = exact_adjacent_weight(n, a, b, z)
                assert dot(z, adjacent) < 0 and max(adjacent) <= 3000*n**6
                has_other_tie = False
                for q in product(range(-n, n+1), repeat=3):
                    before, after = dot(q, midpoint), dot(q, adjacent)
                    if before:
                        assert before*after > 0
                        counts['all_bounded_integer_sign_preservations'] += 1
                    elif any(q) and cross(q, z) != (0, 0, 0):
                        has_other_tie = True
                counts['edge_midpoints_with_additional_bounded_integer_ties'] += int(has_other_tie)
                queue.append((adjacent, B, (a, b)))
        assert set(visited) == set(known)
        assert queries == 1+2*len(graph) <= 6*len(visited)+1
        assert sum(area2([tuple(Q(x) for x in v) for v in cell['polygon']])
                   for cell in visited.values()) == 1
        supporting = {}
        for cell in visited.values():
            for point in cell['polygon']:
                point = tuple(Q(x) for x in point)
                scale = lcm(*(x.denominator for x in point))
                ray = tuple(int(x*scale) for x in point)
                divisor = gcd(*ray)
                ray = tuple(x//divisor for x in ray)
                assert max(ray) <= 2*n*n
                value = dot(ray, cell['profile'])
                assert ray not in supporting or supporting[ray] == value
                supporting[ray] = value
        old_facets = {tuple(-x for x in row['normal']): -row['offset']
                      for row in reference['cut_polytope_facets']
                      if tuple(row['normal']) != (1, 1, 1)}
        assert supporting == old_facets
        counts['independent_complete_walks_and_state_halfspaces'] += 1
        counts['walk_initial_ideals'] += len(visited)
        fixtures.append({'model': label, 'length': n, 'characteristic': p,
                         'initial_ideals': len(visited), 'undirected_adjacencies': len(graph),
                         'basis_queries': queries, 'complete_cells': list(visited.values()),
                         'directed_neighbor_trace': trace,
                         'state_halfspaces': [{'normal': w, 'minimum': h}
                                             for w, h in sorted(supporting.items())]})
    assert counts['edge_midpoints_with_additional_bounded_integer_ties'] > 0
    return {'counts': dict(counts), 'fixtures': fixtures,
            'scope': 'Complete positive-cone traversal independently compared with exhaustive monomial-basis polytopes, the older full reconstruction and direct ground elimination. The checker uses elementary field arithmetic; fast-field bounds rely on the cited Neiger--Schost backend.'}


def run_case(label, model, p, counts):
    n, ev, matrices, one = model
    E = tuple(sorted(hyperbolic_cross(n, 3)))
    U = set(E)
    U.update(tuple(x+(i == j) for j, x in enumerate(a)) for a in E for i in range(3))
    for i, j in combinations(range(3), 2):
        for k in range(n):
            v = tuple(int(l == k) for l in range(n))
            assert apply(matrices[i], apply(matrices[j], v, p), p) == apply(matrices[j], apply(matrices[i], v, p), p)
            counts['commuting_multiplication_columns'] += 1
    for a in sorted(U):
        v = image_from_matrices(matrices, one, a, p)
        expected = tuple(ev(a).get(i, 0) % p if p else ev(a).get(i, 0) for i in range(n))
        assert v == expected
        counts['independent_cyclic_monomial_columns'] += 1
    profiles, positive = exhaustive_bases(E, n, ev, p, counts)
    module=fixed_module(model,p)
    R = n*n+1
    candidates = set(profiles)
    for b in profiles:
        assert sum(b) < R
        candidates.update(tuple(b[i]+(R-sum(b))*(i == j) for i in range(3)) for j in range(3))
    audit_points, trace = set(), []

    def oracle(w, record=False):
        reference,_=rank_trace(sorted(E,key=lambda a:monomial_key(w,a)),ev,p,set())
        B=module_standard_basis(model,p,lambda a:monomial_key(w,a),counts,module)
        assert B==frozenset(reference)
        if record:
            direct=compressed_krylov(matrices,(one,),lambda u:monomial_key(w,u[0]),p,counts)
            assert frozenset(a for a,j in direct)==B
            counts['multivariate_and_fixed_module_equalities']+=1
        b = tuple(sum(a[i] for a in B) for i in range(3))
        assert b in positive and B == positive[b][0]
        q = b
        if max(w) > 0:
            j = max(i for i in range(3) if w[i] == max(w))
            q = tuple(b[i]+(R-sum(b))*(i == j) for i in range(3))
        independent = max(candidates, key=lambda v: (dot(w, v),)+tuple(-x for x in v))
        assert q == independent
        audit_points.add(q)
        counts['explicit_all_basis_support_comparisons'] += 1
        if record:
            trace.append({'direction': list(map(str, w)), 'vertex': q,
                          'initial_profile_label': b, 'on_cut': sum(q) == R})
        return q, frozenset(B)

    for w in product(range(-2, 3), repeat=3):
        oracle(w)
    for w in ((Q(1, 3), Q(1, 3), Q(-2, 7)), (Q(0), Q(0), Q(-1, 5)),
              (Q(-1, 11), Q(2, 13), Q(2, 13))):
        oracle(w)
    b, B = oracle((-1, -1, -1), True)
    points = {b}
    points.update(tuple(b[i]+(R-sum(b))*(i == j) for i in range(3)) for j in range(3))
    labels, verified = {b: B}, set()
    while True:
        H = hull_facets(points)
        pending = sorted(f for f in H if f not in verified)
        if not pending:
            break
        z, h = pending[0]
        q, B = oracle(z, True)
        if dot(z, q) == h:
            verified.add((z, h))
        else:
            assert dot(z, q) > h and q not in points
            points.add(q)
            if sum(q) < R:
                labels[q] = B
    assert all(dot(z, v) <= h for z, h in H for v in candidates)
    vertices = {a for a in points if any(determinant(N) for N in combinations(
        [z for z, h in H if dot(z, a) == h], 3))}
    assert audit_points <= vertices
    assert set(labels) == set(positive) == {v for v in vertices if sum(v) < R}
    assert all(labels[b] == positive[b][0] for b in labels)
    assert len(vertices) <= 4*len(labels) and len(H) <= 2*len(vertices)-4
    assert len(trace) <= 1+len(vertices)+len(H) <= 12*len(labels)+1
    assert all((z, h) == ((1, 1, 1), R) or max(z) <= 0 for z, h in H)
    counts['reconstruction_support_queries'] += len(trace)
    counts['complete_truncated_state_reconstructions'] += 1
    counts['recovered_initial_ideals'] += len(labels)
    recovered = []
    for b, B in sorted(labels.items()):
        _, polygon, weight = positive[b]
        recovered_cone = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for v in sorted(vertices):
            z = tuple(x-y for x, y in zip(v, b))
            if any(z):
                recovered_cone = clip(recovered_cone, z, 1)
        assert area2(recovered_cone) == area2(polygon) > 0
        assert all(dot(tuple(x-y for x, y in zip(c, b)), w) >= 0
                   for c in profiles for w in recovered_cone)
        assert all(dot(tuple(x-y for x, y in zip(v, b)), w) >= 0
                   for v in vertices for w in polygon)
        counts['independent_recovered_normal_cones'] += 1
        recovered.append({'profile': b, 'standard_set': sorted(B),
                          'positive_weight': list(map(str, weight)),
                          'weight_simplex_polygon': [[str(x) for x in v] for v in polygon],
                          'reduced_basis': verify_reduced_basis(B, ev, matrices, weight, p, counts)})
    return {'model': label, 'length': n, 'characteristic': p, 'cutoff': R,
            'ground_size': len(E), 'column_and_border_ground_size': len(U),
            'independent_monomial_bases': sum(map(len, profiles.values())),
            'independent_profile_count': len(profiles), 'initial_ideals': len(labels),
            'cut_polytope_vertices': sorted(vertices),
            'cut_polytope_facets': [{'normal': z, 'offset': h} for z, h in sorted(H)],
            'reconstruction_queries': len(trace), 'query_trace': trace,
            'all_initial_ideals_and_reduced_bases': recovered}


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    counts = defaultdict(int)
    cases = [('one_point', curve_model(1), 0)]
    cases += [('double_point', curve_model(2), p) for p in (0, 2)]
    cases += [('curvilinear_three', curve_model(3), p) for p in (0, 2, 3)]
    cases += [('nonlinear_rectangle', rectangle_fixture({(1, 0): 1, (0, 1): 1, (1, 1): 1}), p) for p in (0, 2)]
    cases += [('monomial_box', rectangle_fixture({}), 3), ('three_distinct_points', point_model(), 5)]
    fixtures = [run_case(label, model, p, counts) for label, model, p in cases]
    prefix=prefix_module_controls(cases)
    digits = prime_digit_controls()
    walk = adjacency_walk_controls(cases, fixtures)
    assert [f['initial_ideals'] for f in fixtures[:9]] == [1, 3, 3, 6, 5, 6, 6, 4, 1]
    sources = ('docs/oracle_fan_enumeration.md', 'scripts/run_oracle_enumeration_checks.py',
               'docs/rational_shadow_bound.md', 'docs/power_carrier_bounds.md',
               'scripts/run_higher_interaction_checks.py', 'scripts/run_hypersurface_truncation_checks.py',
               'scripts/run_signed_face_checks.py', 'scripts/run_liaison_fan_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_extremal_initial_ideals_checks.py', 'scripts/run_socle_stability_checks.py',
               'scripts/run_higher_seed_checks.py', 'scripts/run_planar_dyadic_bound_checks.py')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_oracle_enumeration_checks.py',
              'counts': dict(counts), 'fixtures': fixtures,
              'faithful_fixed_module_prefix_ranks':prefix,
              'prime_digit_rank_law': digits,
              'direct_positive_cone_traversal': walk,
              'uniform_results': {'truncated_state_vertices': 'at most 4*g',
                                  'vertex_queries': 'at most 12*g+1, including initialization',
                                  'field_operations': 'O(n^3*L^2+(g+1)*n^omega*L), L=1+log(n+1)',
                                  'dense_reference_field_operations':'O((g+1)*n^3*L^2)',
                                  'oracle_reconstruction_geometric_reference': 'O((g+1)^2*log(g+2)+(g+1)*n*L^2)',
                                  'geometric_operations': 'O((g+1)*n^(5/3)*L) exact integer/rational operations and cost comparisons via direct positive-cone traversal',
                                  'direct_traversal_basis_queries': 'at most 6*g+1, including initialization',
                                  'direct_neighbor_weight_height': 'O(n^6), hence O(log(n+1))-bit integer exponents/weights',
                                  'unrestricted_shadow_consequence': 'O(n^(omega+23/8)*L^(11/8)) field operations',
                                  'unrestricted_shadow_geometric_consequence': 'O(n^(109/24)*L^(11/8)) rational geometric operations and comparisons',
                                  'combined_unit_cost_consequence': 'O(n^(omega+23/8)*L^(11/8)) field and geometric arithmetic operations combined, plus one-time preprocessing already dominated asymptotically; not a coefficient-bit bound',
                                  'dense_reference_shadow_consequence':'O(n^(47/8)*L^(19/8)) field operations',
                                  'classical_count_composition': 'O(n^(omega+3)*L) field operations, plus one-time O(n^3*L^2) preprocessing',
                                  'dense_reference_classical_composition':'O(n^6*L^2)',
                                  'given_input': 'One marked reduced Groebner basis, or commuting coordinate matrices with a cyclic vector.',
                                  'union_output': 'All reduced bases; their union may retain repeated polynomials.'},
              'new_unrestricted_G3_count': False, 'novelty_established': False,
              'bit_complexity_bound_established': False, 'practical_speedup_established': False,
              'goal_complete': False,
              'scope': 'Neiger--Schost Theorem 1.7 and Proposition 3.7 apply to arbitrary commuting matrices and monomial orders. Direct positive-cone traversal, using classical Groebner cone descriptions and adjacency, is independently compared with exhaustive state polytopes and the retained older oracle reconstruction. Exact polar hulls remove the previous quadratic output geometry budget. Field arithmetic in this checker is elementary; the theorem uses published fast linear algebra. No new general traversal method, unrestricted count, current-best priority, coefficient-bit cost or practical speedup is established.',
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/'results/research/oracle_enumeration_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'previous_counts': report['counts'], 'direct_walk': walk['counts']}, indent=2))


if __name__ == '__main__':
    main()
