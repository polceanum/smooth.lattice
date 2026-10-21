#!/usr/bin/env python3
"""Exact collective coordinate-slice budgets and actual cubic lower walls."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from heapq import heappop, heappush
from itertools import combinations, product
from math import comb, factorial, gcd
from pathlib import Path
from random import Random
import hashlib
import json

from run_interleaving_research import minimal_excluded_border
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_primitive_wall_checks import threshold
from run_structural_complexity_checks import area2, clip, dot

ROOT = Path(__file__).resolve().parents[1]


def profile(B):
    return tuple(sum(a[i] for a in B) for i in range(3))


def primitive(v):
    d = gcd(*v)
    assert d
    q = tuple(a//d for a in v)
    return q if next(a for a in q if a) > 0 else tuple(-a for a in q)


def downset(B):
    return all(tuple(e-(i == j) for j, e in enumerate(a)) in B
               for a in B for i in range(3) if a[i])


def moment_cap(n, h):
    h = min(n, h)
    return (n//h)*comb(h, 2)+comb(n % h, 2)


def hook_model(t, seed, counts):
    """Actual flat thick-hook deformation over GF(2), with matrix replays."""
    a = b = t
    R = S = t*t+1
    gens = ((R, 0), (a, b), (0, S))
    divides = lambda u, v: all(x <= y for x, y in zip(u, v))
    D = tuple(p for p in product(range(R), range(S))
              if not any(divides(g, p) for g in gens))
    index = {p: i for i, p in enumerate(D)}
    m = len(D)
    allowed = (
        [p for p in D if p[0] >= a or p[1] >= S-b],
        [p for p in D if p[0] >= a or p[1] >= b],
        [p for p in D if p[0] >= R-a or p[1] >= b])
    rng = Random(seed)
    tails = [tuple(p for p in A if rng.randrange(2)) for A in allowed]
    cache = {}

    def ev(exponent):
        if exponent in cache:
            return cache[exponent]
        e, *v = exponent
        v = tuple(v)
        out = 0
        if e < 2 and v in index:
            out = 1 << (index[v]+m*e)
        elif not e:
            images = []
            for g, F in zip(gens, tails):
                if divides(g, v):
                    shift = tuple(x-y for x, y in zip(v, g))
                    image = 0
                    for p in F:
                        q = tuple(x+y for x, y in zip(p, shift))
                        if q in index:
                            image ^= 1 << (m+index[q])
                    images.append(image)
            assert images and len(set(images)) == 1
            out = images[0]
        cache[exponent] = out
        return out

    squares = tuple(tuple(x+y for x, y in zip(g, h))
                    for g in gens for h in gens)
    outer = tuple(p for p in product(range(2*R), range(2*S))
                  if not any(divides(g, p) for g in squares))
    fixed = frozenset((0, *p) for p in D)
    ground = tuple(sorted({(0, *p) for p in outer} | {(1, *p) for p in D}))
    variable = tuple(p for p in ground if p not in fixed)
    basis = [(0, *v) for v in D]+[(1, *v) for v in D]
    matrices = [tuple(ev(tuple(e+(i == j) for j, e in enumerate(p)))
                      for p in basis) for i in range(3)]

    def apply(M, v):
        out = 0
        while v:
            j = (v & -v).bit_length()-1
            out ^= M[j]
            v &= v-1
        return out

    for i, j in combinations(range(3), 2):
        for k in range(2*m):
            assert apply(matrices[i], matrices[j][k]) == apply(matrices[j], matrices[i][k])
            counts['independent_commuting_matrix_columns'] += 1
    for p in ground:
        for i in range(3):
            successor = tuple(e+(i == j) for j, e in enumerate(p))
            assert apply(matrices[i], ev(p)) == ev(successor)
            counts['independent_ground_successor_replays'] += 1
    for k in range(2*m):
        assert not apply(matrices[0], matrices[0][k])
        counts['square_zero_matrix_columns'] += 1

    def greedy(w):
        piv = {}
        chosen = []
        for p in sorted(variable, key=lambda p: (dot(w, p), p)):
            v = ev(p) >> m
            while v:
                j = (v & -v).bit_length()-1
                if j not in piv:
                    piv[j] = v
                    chosen.append(p)
                    break
                v ^= piv[j]
            if len(chosen) == m:
                B = fixed | frozenset(chosen)
                assert downset(B)
                return B
        raise AssertionError('rank-deficient carrier')

    def normals(B):
        labels = tuple(sorted(B-fixed))
        piv = {}
        for i, p in enumerate(labels):
            v, co = ev(p) >> m, 1 << i
            while v:
                j = (v & -v).bit_length()-1
                if j not in piv:
                    piv[j] = (v, co)
                    break
                s, d = piv[j]
                v ^= s
                co ^= d
            assert v
        out = set()
        for head in minimal_excluded_border(tuple(B)):
            assert not (ev(head) & ((1 << m)-1))
            v, co = ev(head) >> m, 0
            while v:
                j = (v & -v).bit_length()-1
                s, d = piv[j]
                v ^= s
                co ^= d
            while co:
                i = (co & -co).bit_length()-1
                delta = tuple(a-b for a, b in zip(head, labels[i]))
                d = gcd(*delta)
                out.add(tuple(a//d for a in delta))
                co &= co-1
        return tuple(sorted(out))

    return 2*m, greedy, normals, tails


def complete_hook(t, counts):
    n, greedy, normals, tails = hook_model(t, 271828, counts)
    queue = [((1, n+1, (n+1)**2), None)]
    cells, edges, queries = {}, set(), 0
    while queue:
        w, previous = queue.pop()
        B = greedy(w)
        queries += 1
        if previous is not None:
            assert previous != B
            edges.add(frozenset((previous, B)))
        if B in cells:
            continue
        N = normals(B)
        assert all(dot(q, w) > 0 for q in N)
        poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for q in N:
            poly = clip(poly, q, 1)
        assert area2(poly) > 0
        cells[B] = poly
        for a, b in zip(poly, poly[1:]+poly[:1]):
            if any(a[i] == b[i] == 0 for i in range(3)):
                continue
            q = next(q for q in N if dot(q, a) == dot(q, b) == 0)
            nextw, _ = exact_adjacent_weight(n, a, b, q)
            queue.append((nextw, B))
    assert sum(map(area2, cells.values())) == 1
    assert queries == 1+2*len(edges)
    counts['complete_coupled_hook_fans'] += 1
    counts['complete_hook_oracle_queries'] += queries
    return n, cells, edges, {'thickness': t, 'characteristic': 2,
                            'seed': 271828, 'compatible_tails': tails}


def edges_from_polygons(cells, counts):
    owners = defaultdict(list)
    assert sum(map(area2, cells.values())) == 1
    for B, poly in cells.items():
        for a, b in zip(poly, poly[1:]+poly[:1]):
            if not any(a[i] == b[i] == 0 for i in range(3)):
                owners[tuple(sorted((a, b)))].append(B)
    edges = set()
    for segment, pair in owners.items():
        assert len(pair) == 2
        delta = tuple(a-b for a, b in zip(profile(pair[0]), profile(pair[1])))
        assert all(dot(w, delta) == 0 for w in segment)
        edges.add(frozenset(pair))
        counts['independent_polygon_adjacency_replays'] += 1
    assert len(edges) == len(owners)
    return edges


def budgets(n, cells, edges, caps, counts):
    """Replay planar exposure, per-level charges and every parallel bound."""
    profiles = {B: profile(B) for B in cells}
    T = comb(n, 2)
    assert len(set(profiles.values())) == len(profiles)
    assert all(downset(B) and len(B) == n and sum(v) <= T
               for B, v in profiles.items())
    records = []
    for pair in edges:
        B, C = sorted(pair, key=profile)
        delta = tuple(a-b for a, b in zip(profiles[B], profiles[C]))
        q = primitive(delta)
        cost = threshold(q)
        assert cost <= n and min(q) < 0 < max(q)
        records.append((profiles[B], profiles[C], q, cost, gcd(*delta)))
    rows = []
    for i, h in enumerate(caps):
        j, k = [a for a in range(3) if a != i]
        low = min(v[i] for v in profiles.values())
        high = max(v[i] for v in profiles.values())
        K = moment_cap(n, h)
        assert high <= K
        by_level = defaultdict(list)
        for A, B, q, cost, multiple in records:
            if q[i]:
                continue
            assert A[i] == B[i] and cost == abs(q[j])+abs(q[k])
            level = A[i]
            planar_w = (abs(q[k]), abs(q[j]))
            values = {v for v in profiles.values() if v[i] == level}
            minimum = min(planar_w[0]*v[j]+planar_w[1]*v[k] for v in values)
            exposed = {v for v in values
                       if planar_w[0]*v[j]+planar_w[1]*v[k] == minimum}
            assert exposed == {A, B}
            by_level[level].append((q, cost, multiple))
            counts['independent_positive_slice_exposures'] += 1
        directions = Counter()
        level_rows = []
        for level, walls in sorted(by_level.items()):
            qs = [q for q, _, _ in walls]
            assert len(set(qs)) == len(qs)
            values = [v for v in profiles.values() if v[i] == level]
            variation_bound = sum(max(v[a] for v in values)-min(v[a] for v in values)
                                  for a in (j, k))
            cost = sum(c for _, c, _ in walls)
            weighted = sum(c*d for _, c, d in walls)
            assert cost <= weighted <= variation_bound <= 2*(T-level)
            directions.update(qs)
            level_rows.append({'moment': level, 'walls': len(walls),
                               'primitive_cost': cost, 'edge_variation': weighted,
                               'coordinate_range_bound': variation_bound,
                               'simplex_bound': 2*(T-level)})
            counts['independent_integer_slice_budgets'] += 1
        parallel = []
        for q, p in sorted(directions.items()):
            assert p <= high-low+1
            row = {'direction': q, 'pieces': p, 'moment_width_bound': high-low+1}
            if h == 2:
                residual = n-threshold(q)+1
                bound = min(residual, n//2)+1
                assert p <= bound
                row.update(residual_rank=residual, residual_bound=bound)
                counts['independent_degree_two_residual_bounds'] += 1
            parallel.append(row)
        total = sum(row['primitive_cost'] for row in level_rows)
        bound = (high-low+1)*(n*(n-1)-low-high)
        cap_bound = (K+1)*(n*(n-1)-K)
        assert total <= bound <= cap_bound
        counts['independent_collective_coordinate_budgets'] += 1
        rows.append({'coordinate': i, 'monic_relation_degree': h,
                     'moment_minimum': low, 'moment_maximum': high,
                     'moment_cap': K, 'wall_cost': total, 'width_bound': bound,
                     'degree_bound': cap_bound, 'integer_slices': level_rows,
                     'walls_outside_this_coordinate_family':
                         sum(bool(q[i]) for _, _, q, _, _ in records),
                     'cost_outside_this_coordinate_family':
                         sum(c for _, _, q, c, _ in records if q[i]),
                     'parallel_directions': parallel})
    ordered = sorted(profiles.values())
    return {'length': n, 'initial_ideals': len(cells), 'undirected_walls': len(edges),
            'total_wall_cost': sum(r[3] for r in records), 'complete_area': '1',
            'coordinate_budgets': rows, 'profiles': ordered,
            'walls': [{'endpoints': [A, B], 'direction': q, 'cost': c,
                       'edge_multiple': d}
                      for A, B, q, c, d in sorted(records)]}


def slab_prefix(n, w):
    zero = (0, 0, 0)
    heap, seen, out = [(0, zero)], {zero}, set()
    while len(out) < n:
        _, a = heappop(heap)
        out.add(a)
        for i in range(3):
            b = tuple(e+(i == j) for j, e in enumerate(a))
            if b[0] < 2 and b not in seen:
                seen.add(b)
                heappush(heap, (dot(w, b), b))
    assert downset(out)
    return frozenset(out)


def bareiss_det(A):
    """Integer determinant; every division is independently checked exact."""
    A = [row[:] for row in A]
    previous, sign = 1, 1
    for k in range(len(A)-1):
        pivot = next((i for i in range(k, len(A)) if A[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            A[k], A[pivot] = A[pivot], A[k]
            sign = -sign
        p = A[k][k]
        for i in range(k+1, len(A)):
            for j in range(k+1, len(A)):
                numerator = A[i][j]*p-A[i][k]*A[k][j]
                quotient, remainder = divmod(numerator, previous)
                assert not remainder
                A[i][j] = quotient
            A[i][k] = 0
        previous = p
    return sign*A[-1][-1]


def actual_exponential_basis(n, B, determinant, counts):
    q = n//2
    C = n*(1+n+n*n)
    lambdas = (C*q+1, C+n, C+n*n)
    columns = sorted((q*a+b+c, dot(lambdas, (a, b, c))) for a, b, c in B)
    assert all(v <= i for i, (v, _) in enumerate(columns))
    assert all(x < y for (_, x), (_, y) in zip(columns, columns[1:]))
    counts['ordered_Hermite_actual_basis_conditions'] += 1
    if not determinant:
        return None
    A = [[0 if r < v else factorial(r)//factorial(r-v)*xi**(r-v)
          for v, xi in columns] for r in range(n)]
    det = bareiss_det(A)
    assert det
    counts['independent_characteristic_zero_integer_determinants'] += 1
    return {'nonzero': True, 'absolute_bit_length': abs(det).bit_length(),
            'sha256_of_signed_decimal': hashlib.sha256(str(det).encode()).hexdigest()}


def sharp_family(m, counts):
    n = 16*m
    pieces, directions, total = [], Counter(), 0
    for k in range(1, 2*m+1):
        for c in range(10*m, 12*m+1):
            if (c-k) % 2 != 1:
                continue
            T = (n-k+c-1)//2
            w = (n+c-3*k, 2, 2*c)
            assert w[0] == 2*(T-k)+1 and c < T < 2*c
            scale = 4*n*n
            weights = [tuple(scale*a+(s if i == 2 else 0)
                             for i, a in enumerate(w)) for s in (-1, 1)]
            B, D = [slab_prefix(n, v) for v in weights]
            pair = {(0, T, 0), (0, T-c, 1)}
            assert B ^ D == pair and len(B-D) == len(D-B) == 1
            assert profile(B)[0] == profile(D)[0] == k
            delta = tuple(a-b for a, b in zip(profile(B), profile(D)))
            direction = primitive(delta)
            assert direction == (0, c, -1) and threshold(direction) == c+1
            assert all(dot(w, a) <= 2*T for a in B | D)
            dets = [actual_exponential_basis(n, A, m <= 2, counts) for A in (B, D)]
            directions[c] += 1
            total += c+1
            counts['independent_actual_lower_wall_swaps'] += 1
            pieces.append({'upper_layer_count': k, 'direction': direction,
                           'cost': c+1, 'exact_wall_weight': w,
                           'endpoint_profiles': [profile(B), profile(D)],
                           'integer_determinants': dets})
    assert len(pieces) == m*(2*m+1)
    assert all(p == m for p in directions.values())
    assert total == m*(2*m+1)*(11*m+1)
    return {'m': m, 'length': n, 'actual_wall_pieces': len(pieces),
            'proved_subset_cost': total, 'parallel_multiplicity_per_direction': m,
            'whole_coordinate_cost_upper_bound':
                (n//2+1)*(n*(n-1)-n//2), 'pieces': pieces}


def main():
    counts = Counter()
    hooks = []
    for t in (2, 3):
        n, cells, edges, metadata = complete_hook(t, counts)
        assert edges_from_polygons(cells, counts) == edges
        row = budgets(n, cells, edges, (2, 2*(t*t+1), 2*(t*t+1)), counts)
        row.update(metadata)
        hooks.append(row)
        print('complete coupled hook:', t, n, len(cells), len(edges), flush=True)
    capped = json.loads((ROOT/'results/research/capped_jet_checks.json').read_text())
    assert capped['all_checks_passed']
    saved = []
    for fixture in capped['complete_quotient_fans']:
        n, r = fixture['length'], fixture['jet_order']
        cells = {frozenset(map(tuple, c['standard_set'])):
                 [tuple(map(Q, v)) for v in c['positive_cone']]
                 for c in fixture['all_cells']}
        edges = edges_from_polygons(cells, counts)
        assert len(edges) == fixture['undirected_walls']
        H = r*len(fixture['nodes'])
        row = budgets(n, cells, edges, (r, H, H), counts)
        row.update(jet_order=r, characteristic=fixture['characteristic'],
                   nodes=fixture['nodes'])
        saved.append(row)
        counts['saved_complete_capped_jet_fans'] += 1
    lower = [sharp_family(m, counts) for m in range(1, 5)]
    sources = ['docs/coordinate_slice_wall_budget.md',
               'scripts/run_coordinate_slice_budget_checks.py',
               'docs/nested_jet_matroids.md', 'docs/wall_cost_products.md',
               'docs/wall_coexistence.md', 'docs/primitive_wall_spectrum.md',
               'docs/interleaving_research.md',
               'results/research/capped_jet_checks.json',
               'scripts/run_interleaving_research.py',
               'scripts/run_oracle_enumeration_checks.py',
               'scripts/run_primitive_wall_checks.py',
               'scripts/run_structural_complexity_checks.py']
    out = {
        'schema_version': 1, 'all_checks_passed': True,
        'scope': 'Exact finite controls for a uniform coordinate-slice theorem. '
                 'Sharp cubic cost concerns only walls with zero bounded-coordinate '
                 'component. No full spatial budget, unrestricted improvement, '
                 'polynomially separating family or novelty established.',
        'uniform_theorem': {
            'parallel': 'p_I(0,a,-b)<=k_max-k_min+1',
            'collective': 'V_x^0<=(k_max-k_min+1)*(n*(n-1)-k_min-k_max)',
            'bounded_coordinate': 'V_x^0<=(K_h(n)+1)*(n*(n-1)-K_h(n))=O(h*n^3)',
            'degree_two_residual': 'p_I(q)<=min(n-nu(q)+1,floor(n/2))+1 for q_x=0',
            'sharpness': 'n=16m; m*(2m+1) actual walls with cost sum '
                         'm*(2m+1)*(11m+1); x^2=0 over Q, also radical two-plane realizations.',
            'proof_not_finite_fit': 'docs/coordinate_slice_wall_budget.md'},
        'counters': dict(counts), 'complete_coupled_hook_fans': hooks,
        'saved_capped_jet_budget_replays': saved, 'sharp_slab_families': lower,
        'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest()
                          for s in sources}}
    path = ROOT/'results/research/coordinate_slice_budget_checks.json'
    path.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'all_checks_passed': True, 'counters': dict(counts),
                      'report': str(path.relative_to(ROOT))}), flush=True)


if __name__ == '__main__':
    main()
