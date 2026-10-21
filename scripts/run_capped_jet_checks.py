#!/usr/bin/env python3
"""Exact capped-prefix and noncurvilinear planar-jet fan controls."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd
from pathlib import Path
from random import Random
import hashlib
import json

from run_oracle_enumeration_checks import exact_adjacent_weight
from run_primitive_wall_checks import threshold
from run_signed_face_checks import field_basis, field_reduce, normals
from run_structural_complexity_checks import (
    area2, clip, dot, finite_staircases, hyperbolic_cross,
)

ROOT = Path(__file__).resolve().parents[1]


def jet_model(r, nodes, p, seed):
    """Actual product of K[u,v]/(u,v)^r; x has arbitrary radical jets."""
    local = tuple((a, b) for a in range(r) for b in range(r-a))
    index = {a: j for j, a in enumerate(local)}
    L = len(local)
    n = L*len(nodes)
    rng = Random(seed)
    radical = [{a: rng.randrange(1, 20) for a in local if any(a)}
               for _ in nodes]
    one = {(0, 0): 1}

    def normalize(f):
        return {a: c % p if p else c for a, c in f.items()
                if (c % p if p else c)}

    def multiply(f, g):
        out = defaultdict(int)
        for a, c in f.items():
            for b, d in g.items():
                key = tuple(x+y for x, y in zip(a, b))
                if sum(key) < r:
                    out[key] += c*d
        return normalize(out)

    powers = []
    for f in radical:
        row = [one]
        for k in range(r):
            row.append(multiply(row[-1], f))
        assert not row[r]
        powers.append(row)

    @lru_cache(None)
    def ev(alpha):
        e, a, b = alpha
        if e >= r:
            return {}
        out = {}
        for j, (A, C) in enumerate(nodes):
            # Closed Taylor coefficients; no factorial division in small char.
            plane = {(i, k): comb(a, i)*A**(a-i)*comb(b, k)*C**(b-k)
                     for i in range(min(a, r-1)+1)
                     for k in range(min(b, r-1-i)+1)}
            value = multiply(powers[j][e], normalize(plane))
            out.update({j*L+index[t]: c for t, c in value.items()})
        return out

    def apply(i, column):
        out = {}
        for j, (A, C) in enumerate(nodes):
            f = {a: column.get(j*L+k, 0) for k, a in enumerate(local)
                 if column.get(j*L+k, 0)}
            coordinate = (radical[j], {(0, 0): A, (1, 0): 1},
                          {(0, 0): C, (0, 1): 1})[i]
            value = multiply(f, coordinate)
            out.update({j*L+index[t]: c for t, c in value.items()})
        return out

    return n, ev, apply, radical


def deterministic_model(r, m, p, seed):
    rng = Random(seed)
    nodes = [(j+1, b) for j, b in enumerate(rng.sample(range(2, 100), m))]
    if p:
        nodes = [(a % p, b % p) for a, b in nodes]
    assert len(set(nodes)) == m
    return nodes, jet_model(r, nodes, p, seed+997)


def algebra_controls(r, nodes, p, model, counts):
    n, ev, apply, _ = model
    for i, j in combinations(range(3), 2):
        for k in range(n):
            assert apply(i, apply(j, {k: 1})) == apply(j, apply(i, {k: 1}))
            counts['independent_commuting_matrix_columns'] += 1
    E = tuple(hyperbolic_cross(n, 3))
    assert len(field_basis(E, ev, p)[0]) == n
    for a in E:
        for i in range(3):
            successor = tuple(e+(i == j) for j, e in enumerate(a))
            assert apply(i, ev(a)) == ev(successor)
            counts['closed_Taylor_vs_matrix_successor_replays'] += 1
    # Replay each asserted carrier polynomial, not just its leading power.
    for i in range(3):
        coefficients = [1]
        roots = [0]*r if not i else [t[i-1] for t in nodes for _ in range(r)]
        for a in roots:
            row = [0]*(len(coefficients)+1)
            for k, c in enumerate(coefficients):
                row[k] -= a*c
                row[k+1] += c
            coefficients = row
        image = defaultdict(int)
        for k, c in enumerate(coefficients):
            alpha = tuple(k if j == i else 0 for j in range(3))
            for t, d in ev(alpha).items():
                image[t] += c*d
        assert all(not (c % p if p else c) for c in image.values())
        counts['independent_monic_carrier_replays'] += 1
    return E


def rectangle_controls(r, m, p, ev, counts):
    """Check minimal required rectangles; all larger ones follow by inclusion."""
    H = r*m
    required = [(a, b) for a in range(r, H+1) for b in range(r, H+1)
                if (a//r)*(b//r) >= m]
    minimal = [(a, b) for a, b in required
               if not any(c <= a and d <= b and (c, d) != (a, b)
                          for c, d in required)]
    n = m*r*(r+1)//2
    rows = []
    for a, b in minimal:
        rectangle = tuple((0, i, j) for i in range(a) for j in range(b))
        rk = len(field_basis(rectangle, ev, p)[0])
        assert rk == n, (r, m, p, a, b, rk, n)
        rows.append({'side_counts': [a, b], 'rank': rk})
        counts['independent_minimal_saturating_rectangles'] += 1
    assert all(any(c <= a and d <= b for c, d in minimal) for a, b in required)
    counts['required_rectangles_covered_by_rank_monotonicity'] += len(required)
    return rows


def prefix_controls(r, m, p, ev, weights, counts):
    H = r*m
    D = tuple(product(range(r), range(H), range(H)))
    n = m*r*(r+1)//2
    limit = min(len(D), 72*m+1 if r == 2 else 16*r**3*m+1)
    E = tuple(hyperbolic_cross(n, 3))
    rows = []
    for w in weights:
        ordered = sorted(D, key=lambda a: (dot(w, a), a))
        selected, _ = field_basis(ordered[:limit], ev, p)
        assert len(selected) == n
        direct, _ = field_basis(sorted(E, key=lambda a: (dot(w, a), a)), ev, p)
        assert set(selected) == set(direct)
        positions = {a: j for j, a in enumerate(ordered)}
        stop = max(positions[a]+1 for a in selected)
        tau = dot(w, ordered[stop-1])
        caps = tuple(min(h-1, tau//(4*cost))
                     for h, cost in zip((r, H, H), w))
        cheap_box = tuple(product(*(range(a+1) for a in caps)))
        assert all(dot(w, a) < tau for a in cheap_box)
        assert len(field_basis(cheap_box, ev, p)[0]) < n
        assert stop <= 64*len(cheap_box)
        rows.append({'weight': list(map(str, w)), 'terms_scanned': stop,
                     'prefix_limit': limit, 'cheap_box_caps': list(map(int, caps)),
                     'cheap_box_cardinality': len(cheap_box)})
        counts['actual_capped_prefix_vs_universal_greedy_replays'] += 1
    return rows


def complete_fan(r, nodes, p, model, counts):
    """Generic exact walker; no curve pole or uniform-rank assumptions."""
    n, ev, _, _ = model
    E = tuple(hyperbolic_cross(n, 3))
    queue = [((1, n+1, (n+1)**2), None)]
    cells, edges, queries, certificates = {}, set(), 0, []
    while queue:
        w, previous = queue.pop()
        chosen, _ = field_basis(sorted(E, key=lambda a: (dot(w, a), a)), ev, p)
        B = frozenset(chosen)
        assert len(B) == n and all(a[0] < r for a in B)
        queries += 1
        if previous is not None:
            assert B != previous
            edges.add(frozenset((B, previous)))
        if B in cells:
            continue
        constraints, _ = normals(tuple(sorted(B)), ev, p)
        assert all(dot(w, q) > 0 for q in constraints)
        poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for q in constraints:
            poly = clip(poly, q, 1)
        assert area2(poly) > 0
        cells[B] = poly
        certificates.append(w)
        for a, b in zip(poly, poly[1:]+poly[:1]):
            if any(a[i] == b[i] == 0 for i in range(3)):
                continue
            q = next(q for q in constraints if dot(q, a) == dot(q, b) == 0)
            neighbor, _ = exact_adjacent_weight(n, a, b, q)
            queue.append((neighbor, B))
    assert sum(map(area2, cells.values())) == 1
    assert queries == 1+2*len(edges)
    assert len(edges) <= max(3, 3*len(cells)-6)
    costs = Counter()
    for pair in edges:
        B, C = tuple(pair)
        delta = tuple(sum(a[j] for a in B)-sum(a[j] for a in C) for j in range(3))
        factor = gcd(*delta)
        q = tuple(e//factor for e in delta)
        cost = threshold(q)
        assert min(q) < 0 < max(q) and cost <= n
        costs[cost] += 1
        counts['independent_complete_fan_wall_costs'] += 1
    if n <= 6:
        exhaustive = set()
        for B in finite_staircases(n, 3):
            counts['independent_exhaustive_downset_candidates'] += 1
            if len(field_basis(tuple(sorted(B)), ev, p)[0]) != n:
                continue
            constraints, _ = normals(tuple(sorted(B)), ev, p)
            poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
            for q in constraints:
                poly = clip(poly, q, 1)
            if area2(poly):
                exhaustive.add(B)
        assert exhaustive == set(cells)
        counts['independent_exhaustive_walk_equalities'] += 1
    prefixes = prefix_controls(r, len(nodes), p, ev, certificates, counts)
    counts['independent_complete_actual_quotient_fans'] += 1
    return {'jet_order': r, 'nodes': nodes, 'characteristic': p, 'length': n,
            'initial_ideals': len(cells), 'undirected_walls': len(edges),
            'primitive_cost_histogram': dict(sorted(costs.items())),
            'total_wall_cost': sum(k*v for k, v in costs.items()),
            'oracle_queries': queries, 'complete_area': '1',
            'prefix_replays': prefixes,
            'all_cells': [{'standard_set': sorted(B),
                           'positive_cone': [list(map(str, v)) for v in poly]}
                          for B, poly in cells.items()]}


def tensor_grids(counts):
    rows = []
    for r, q, s, p in ((2, 1, 3, 0), (2, 2, 2, 2), (2, 2, 3, 3),
                        (3, 2, 2, 2), (3, 2, 2, 5), (4, 2, 2, 2)):
        nodes = list(product(range(q), range(s)))
        n, ev, _, _ = jet_model(r, nodes, p, 100+r)
        rect = tuple((0, a, b) for a in range(r*q) for b in range(r*s))
        chosen, _ = field_basis(rect, ev, p)
        assert len(chosen) == n
        rows.append({'jet_order': r, 'grid': [q, s], 'characteristic': p,
                     'rectangle': [r*q, r*s], 'length': n, 'rank': len(chosen)})
        counts['independent_tensor_Hermite_grid_ranks'] += 1
    return rows


def capped_box_controls(r, nodes, p, model, counts):
    n, ev, _, _ = model
    H = r*len(nodes)
    s, witness = 0, None
    for caps in product(range(r), range(H), range(H)):
        box = tuple(product(*(range(a+1) for a in caps)))
        if len(field_basis(box, ev, p)[0]) < n and len(box) > s:
            s, witness = len(box), caps
        counts['independent_capped_subbox_rank_checks'] += 1
    assert witness is not None
    w = tuple(Q(1, a+1) for a in witness)
    D = tuple(product(range(r), range(H), range(H)))
    subbox = tuple(product(*(range(a//3+1) for a in witness)))
    cheap = tuple(a for a in D if dot(w, a) < 1)
    assert set(subbox) <= set(cheap)
    assert len(subbox)*27 >= s
    assert len(field_basis(cheap, ev, p)[0]) < n
    counts['independent_capped_scanning_lower_witnesses'] += 1
    return {'length': n, 'characteristic': p, 's_D': s,
            'maximizing_nonsaturating_caps': witness,
            'lower_witness_weight': list(map(str, w)),
            'strictly_cheaper_nonsaturating_terms': len(cheap),
            'claim_limit': 'Exact s_D; no maximum scanning count inferred from one fan certificate per cell.'}


def negative_controls(counts):
    nodes = [(1, 3), (2, 7)]
    n, ev, _, _ = jet_model(2, nodes, 0, 13)
    quadratics = tuple((0, a, b) for a in range(3) for b in range(3-a))
    rk = len(field_basis(quadratics, ev, 0)[0])
    assert n == len(quadratics) == 6 and rk == 5
    # The line is z=4y-1; square it and replay directly.
    relation = {(0, 0, 2): 1, (0, 1, 1): -8, (0, 0, 1): 2,
                (0, 2, 0): 16, (0, 1, 0): -8, (0, 0, 0): 1}
    image = defaultdict(int)
    for a, c in relation.items():
        for j, d in ev(a).items():
            image[j] += c*d
    assert not any(image.values())
    nodes = [(j, j) for j in range(1, 10)]
    n, ev, _, _ = jet_model(2, nodes, 0, 29)
    rectangle = tuple((0, a, b) for a in range(6) for b in range(6))
    rk_line = len(field_basis(rectangle, ev, 0)[0])
    assert (6//2)**2 >= 9 and n == 27 and rk_line == 20
    counts['independent_negative_rank_controls'] += 2
    return {'nonuniform_two_point_quadratics': {'length': 6, 'columns': 6, 'rank': rk,
                'vanishing_relation': '((z-4*y+1)^2)'},
            'special_collinear_failure': {'points': 9, 'jet_order': 2,
                'rectangle': [6, 6], 'rank': rk_line, 'length': n,
                'numerical_generic_condition_holds': True}}


def large_defect_small_fans(counts):
    rows = []
    for m in (2, 3, 4, 7):
        n = 2*m
        def ev(a):
            k = m*a[0]+a[1]+a[2]
            return {k: 1} if k < n else {}
        D = tuple(product(range(2), range(n), range(n)))
        s = 0
        for A, B, C in product((1, 2), range(1, n+1), range(1, n+1)):
            box = tuple(product(range(A), range(B), range(C)))
            rk = len(field_basis(box, ev, 0)[0])
            expected = min(n, B+C-1) if A == 1 else (
                n if B+C >= m+1 else 2*(B+C-1))
            assert rk == expected
            if rk < n:
                s = max(s, len(box))
            counts['independent_large_defect_box_ranks'] += 1
        assert s == m*m
        cells = []
        for j in (1, 2):
            for cheap in (False, True):
                B = frozenset(tuple((int(cheap)*i if k == 0 else a if k == j else 0)
                                    for k in range(3))
                              for i in range(2 if cheap else 1)
                              for a in range(m if cheap else n))
                w = [1 if cheap else m+1, 3, 3]
                w[j] = 2 if cheap else 1
                chosen, _ = field_basis(sorted(D, key=lambda a: (dot(w, a), a)), ev, 0)
                assert frozenset(chosen) == B and len(B) == n
                constraints, _ = normals(tuple(sorted(B)), ev, 0)
                poly = [tuple(Q(i == k) for i in range(3)) for k in range(3)]
                for q in constraints:
                    poly = clip(poly, q, 1)
                assert area2(poly) > 0 and all(dot(w, q) > 0 for q in constraints)
                cells.append((B, poly))
        assert len({B for B, _ in cells}) == 4
        assert sum(area2(poly) for _, poly in cells) == 1
        if n <= 6:
            exhaustive = set()
            for B in finite_staircases(n, 3):
                if len(field_basis(tuple(sorted(B)), ev, 0)[0]) != n:
                    continue
                constraints, _ = normals(tuple(sorted(B)), ev, 0)
                poly = [tuple(Q(i == k) for i in range(3)) for k in range(3)]
                for q in constraints:
                    poly = clip(poly, q, 1)
                if area2(poly):
                    exhaustive.add(B)
            assert exhaustive == {B for B, _ in cells}
            counts['independent_large_defect_exhaustive_fans'] += 1
        W = n*n+1
        weight = ((n+1)**3*W, W, W+1)
        ordered = sorted(D, key=lambda a: (dot(weight, a), a))
        chosen, _ = field_basis(ordered, ev, 0)
        stop = max(ordered.index(a)+1 for a in chosen)
        assert stop == n*(n-1)//2+1
        rows.append({'parameter': m, 'length': n, 's_D': s, 'initial_ideals': 4,
                     'complete_area': '1', 'weight': weight,
                     'quadratic_capped_scan': stop})
        counts['independent_large_defect_four_cone_fans'] += 1
    return rows


def main():
    counts = defaultdict(int)
    fans, rectangles, boxes = [], [], []
    for r, m, p in ((2, 1, 0), (2, 2, 0), (2, 2, 1009),
                     (2, 3, 1009), (2, 4, 1009), (2, 5, 1009),
                     (3, 1, 2), (3, 2, 1009)):
        nodes, model = deterministic_model(r, m, p, 12345+31*r+m)
        algebra_controls(r, nodes, p, model, counts)
        minimal = rectangle_controls(r, m, p, model[1], counts)
        fan = complete_fan(r, nodes, p, model, counts)
        fan['radical_coordinate_local_coefficients'] = [
            [{'exponents': a, 'coefficient': c} for a, c in sorted(f.items())]
            for f in model[3]]
        fans.append(fan)
        rectangles.append({'jet_order': r, 'points': m, 'characteristic': p,
                           'minimal_required_rectangles': minimal})
        if m <= 3:
            boxes.append(capped_box_controls(r, nodes, p, model, counts))
    # A larger rank-only control makes the oversampling limit smaller than D.
    r, m, p = 2, 20, 1009
    nodes, model = deterministic_model(r, m, p, 91751)
    minimal = rectangle_controls(r, m, p, model[1], counts)
    weights = [(1, 10**6+3, 10**12+7), (10**12+39, 1, 10**6+9),
               (10**6+7, 10**12+39, 1), (127, 1031, 911),
               (99991, 7, 13), (3, 99991, 17)]
    large = {'jet_order': r, 'nodes': nodes, 'characteristic': p,
             'length': model[0], 'carrier_box_cardinality': r*(r*m)**2,
             'minimal_required_rectangles': minimal,
             'prefix_replays': prefix_controls(r, m, p, model[1], weights, counts)}
    grids = tensor_grids(counts)
    guards = negative_controls(counts)
    defect_fans = large_defect_small_fans(counts)
    sources = ['docs/capped_jet_fans.md', 'scripts/run_capped_jet_checks.py',
               'docs/nested_jet_matroids.md', 'docs/wall_rank_filters.md',
               'scripts/run_signed_face_checks.py', 'scripts/run_oracle_enumeration_checks.py',
               'scripts/run_structural_complexity_checks.py', 'scripts/run_primitive_wall_checks.py']
    report = {'schema_version': 1, 'all_checks_passed': True,
              'scope': 'Finite exact controls for a capped-rank theorem and generic planar fat-point class bound. Full fan counts are not fitted to a growth law. No unrestricted improvement, separating ideal family, class sharpness or novelty established.',
              'uniform_theorem': {'generic_fat_points': 'For n=m*binom(r+1,2), generic planar nodes and every radical third coordinate give g=O(r^2*n^2*log(n+1)*(1+log r)); fixed r>=2 gives O_r(n^2*log(n+1)).',
                                  'field_scope': 'Every infinite field has a nonempty open set by tensor CRT and nonzero minors. Finite-field fixtures verify the explicit rectangular-rank condition; existence over every finite field is not asserted.',
                                  'scan_limit': 'M_D<=min(|D|,16*r^3*m+1)<=32*r*n+1; r=2 also M_D<=min(|D|,72*m+1).',
                                  'class_boundary': 'Not all x^r-containing ideals, not arbitrary node configurations, not uniform matroid ranks.'},
              'counters': dict(counts), 'complete_quotient_fans': fans,
              'rectangular_rank_controls': rectangles, 'capped_box_controls': boxes,
              'larger_rank_only_fixture': large, 'tensor_grid_controls': grids,
              'negative_controls': guards,
              'large_capped_defect_four_initial_fans': defect_fans,
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    target = ROOT/'results/research/capped_jet_checks.json'
    target.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'all_checks_passed': True, 'counters': dict(counts),
                      'complete_fans': [{'r': f['jet_order'], 'n': f['length'],
                                         'p': f['characteristic'], 'initials': f['initial_ideals']}
                                        for f in fans],
                      'larger_scan_max': max(t['terms_scanned'] for t in large['prefix_replays']),
                      'report': str(target.relative_to(ROOT))}, indent=2))


if __name__ == '__main__':
    main()
