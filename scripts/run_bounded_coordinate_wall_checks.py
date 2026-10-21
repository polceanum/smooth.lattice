#!/usr/bin/env python3
"""Translated-cross token charges, residual covers and fixed-direction walls."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
from math import factorial, gcd
from pathlib import Path
from random import Random
import hashlib
import json

from run_coordinate_slice_budget_checks import (
    bareiss_det, complete_hook, downset, primitive, profile, slab_prefix,
)
from run_planar_dyadic_bound_checks import convex_hull, matrix_oracle
from run_primitive_wall_checks import threshold
from run_structural_complexity_checks import dot, hyperbolic_cross
from run_wall_coexistence_checks import core_data, outside_divisors

ROOT = Path(__file__).resolve().parents[1]


def translated_rectangles(H, scale, anchors):
    rects = []
    for a, b in anchors:
        for j in range(H.bit_length()):
            rects.append((a+scale*(2**j-1), a+scale*(2**(j+1)-2),
                          b, b+scale*(H//2**j-1)))
    return rects


def sweep(ground, rank, independent, signs):
    """Adjacent swaps with direct greedy replays; labels may be negative."""
    sigma, tau = signs
    N = len(ground)

    def greedy(weight):
        B = set()
        for j in sorted(range(N), key=lambda j: (dot(weight, ground[j]), j)):
            if independent(frozenset(B | {j})):
                B.add(j)
        assert len(B) == rank
        return B

    events = defaultdict(set)
    for i, j in combinations(range(N), 2):
        dx, dy = (a-b for a, b in zip(ground[i], ground[j]))
        if sigma*tau*dx*dy < 0:
            events[Q(-tau*dy, sigma*dx)].update((i, j))
    events = sorted(events.items())
    order = sorted(range(N), key=lambda j: (tau*ground[j][1], sigma*ground[j][0], j))
    positions = {j: i for i, j in enumerate(order)}
    first = events[0][0]/2 if events else Q(1)
    B = greedy((sigma*first.numerator, tau*first.denominator))
    token = {j: i for i, j in enumerate(sorted(B))}
    paths = [[ground[j]] for j in sorted(B)]
    changes, moves = [], []
    states = [tuple(sum(ground[j][i] for j in B) for i in range(2))]
    for k, (t, involved) in enumerate(events):
        groups = defaultdict(list)
        for j in involved:
            groups[sigma*ground[j][0]*t.numerator+tau*ground[j][1]*t.denominator].append(j)
        changed = False
        for group in groups.values():
            places = sorted(positions[j] for j in group)
            assert places == list(range(places[0], places[-1]+1))
            for loc, j in enumerate(sorted(group, key=lambda j: (sigma*ground[j][0], j)), places[0]):
                while positions[j] > loc:
                    pos = positions[j]
                    e, f = order[pos-1], j
                    if e in B and f not in B and independent(frozenset(B-{e} | {f})):
                        A, D = ground[e], ground[f]
                        dx, dy = sigma*(A[0]-D[0]), tau*(D[1]-A[1])
                        assert dx > 0 and dy > 0 and Q(dy, dx) == t
                        i = token.pop(e)
                        token[f] = i
                        paths[i].append(D)
                        moves.append((t, i, A, D))
                        B.remove(e)
                        B.add(f)
                        changed = True
                    order[pos-1], order[pos] = f, e
                    positions[f], positions[e] = pos-1, pos
        nxt = events[k+1][0] if k+1 < len(events) else t+2
        check = (t+nxt)/2
        assert B == greedy((sigma*check.numerator, tau*check.denominator))
        if changed:
            changes.append(t)
            states.append(tuple(sum(ground[j][i] for j in B) for i in range(2)))
    return states, changes, moves, paths


def cover_sweeps(counts):
    fixtures = []
    for H, r, scale, anchors in (
        (1, 1, 1, ((0, 0), (100, -100))),
        (4, 2, 1, ((0, 0), (1, 1), (-7, 3))),
        (5, 3, 3, ((0, 0), (3, 0))),
        (7, 3, 2, ((-100, 70), (2, -3))),
        (9, 4, 1, ((0, 0), (0, 0))),
        (512, 1, 3, ((1000000, -1000000),)),
    ):
        cross = ([(8, 30), (9, 23), (10, 18)] if H == 512 else
                 sorted(hyperbolic_cross(H, 2)))
        ground = tuple((a+scale*x, b+scale*y) for a, b in anchors for x, y in cross)
        rects = translated_rectangles(H, scale, anchors)
        inside = lambda P, R: R[0] <= P[0] <= R[1] and R[2] <= P[1] <= R[3]
        assert all(any(inside(P, R) for R in rects) for P in ground)
        lines = sorted({(i, Q(2*v+s, 2)) for R in rects
                        for i, v, s in ((0, R[0], -1), (0, R[1], 1),
                                         (1, R[2], -1), (1, R[3], 1))})
        rng = Random(10000*H+100*r+scale)
        columns = [tuple(int(i == j) for i in range(r)) if j < r else
                   tuple(rng.randrange(101) for _ in range(r)) for j in range(len(ground))]
        independent = matrix_oracle(columns, r, 101)
        exposed, quadrant_rows = set(), []
        for signs in product((-1, 1), repeat=2):
            states, changes, moves, paths = sweep(ground, r, independent, signs)
            exposed.update(states)
            crossings = Counter()
            local = defaultdict(lambda: defaultdict(lambda: [0, 0]))
            token_local = defaultdict(lambda: [0, 0])
            crossing_moves = 0
            for t, token, A, D in moves:
                hit = [(i, line) for i, line in lines if (A[i]-line)*(D[i]-line) < 0]
                if hit:
                    crossing_moves += 1
                    crossings.update((token, i, line) for i, line in hit)
                else:
                    j = next(j for j, R in enumerate(rects) if inside(A, R))
                    assert inside(D, rects[j])
                    dx, dy = abs(A[0]-D[0]), abs(A[1]-D[1])
                    local[j][t][0] += dx
                    local[j][t][1] += dy
                    token_local[token, j][0] += dx
                    token_local[token, j][1] += dy
                counts['independent_translated_cover_exchanges'] += 1
            assert all(c == 1 for c in crossings.values())
            assert crossing_moves <= len(crossings) <= r*len(lines) <= 4*r*len(rects)
            budgets = []
            for j, slopes in sorted(local.items()):
                R = rects[j]
                W, Z = R[1]-R[0], R[3]-R[2]
                assert W*Z <= scale*scale*H
                for token in range(r):
                    dx, dy = token_local[token, j]
                    assert dx <= W and dy <= Z
                A = sum(v[0] for v in slopes.values())
                D = sum(v[1] for v in slopes.values())
                K = len(slopes)
                pairs = {(x//gcd(x, y), y//gcd(x, y)) for x, y in slopes.values()}
                assert len(pairs) == K and K**3 <= 27*A*D <= 27*r*r*scale*scale*H
                budgets.append({'rectangle': R, 'parameters': K,
                                'total_coordinate_variation': [A, D]})
                counts['independent_translated_rectangle_budgets'] += 1
            assert len(changes) <= crossing_moves+sum(len(v) for v in local.values())
            quadrant_rows.append({'signs': signs, 'states': len(states),
                                  'crossing_exchanges': crossing_moves,
                                  'rectangle_budgets': budgets,
                                  'token_paths': paths})
            counts['independent_translated_quadrant_sweeps'] += 1
        if r <= 3:
            profiles = {tuple(sum(ground[j][i] for j in B) for i in range(2))
                        for B in combinations(range(len(ground)), r)
                        if independent(frozenset(B))}
            assert exposed == set(convex_hull(profiles))
            counts['independent_exhaustive_translated_profile_hulls'] += 1
        fixtures.append({'H': H, 'rank': r, 'scale': scale, 'anchors': anchors,
                         'ground_size': len(ground), 'quadrants': quadrant_rows})
    assert counts['independent_translated_rectangle_budgets'] > 0
    return fixtures


def ceiling_cube_root(n):
    lo, hi = 0, 1
    while hi**3 < n:
        hi *= 2
    while lo+1 < hi:
        mid = (lo+hi)//2
        if mid**3 >= n:
            hi = mid
        else:
            lo = mid
    return hi


def actual_cover_checks(n, cells, edges, h, counts):
    groups = defaultdict(list)
    cost_opp, cost_zero, cost_remaining = 0, 0, 0
    for pair in edges:
        B, D = tuple(pair)
        q = primitive(tuple(a-b for a, b in zip(profile(B), profile(D))))
        if not q[0]:
            cost_zero += threshold(q)
            continue
        if q[0] < 0:
            q = tuple(-a for a in q)
        groups[q].append((B, D))
        if q[1] <= 0 and q[2] <= 0:
            cost_opp += threshold(q)
        else:
            assert q[1]*q[2] < 0
            cost_remaining += threshold(q)
    rows = []
    for q, walls in sorted(groups.items()):
        a = q[0]
        assert 1 <= a < h
        u, v, core, tips = core_data(q)
        r = n-threshold(q)+1
        project = lambda alpha: (a*alpha[1]-q[1]*alpha[0],
                                 a*alpha[2]-q[2]*alpha[0])
        shadow_points = set()
        endpoints = set()
        for B, D in walls:
            assert core <= B and core <= D
            p, t = project(profile(B)), project(profile(D))
            assert p == t and p not in shadow_points
            shadow_points.add(p)
            endpoints.update((B, D))
        for B in endpoints:
            assert all(alpha[0] < h for alpha in B)
            for alpha in B-core:
                assert outside_divisors(alpha, u, v) <= r
                tip = next(t for t in tips if all(b <= c for b, c in zip(t, alpha)))
                beta = tuple(c-b for b, c in zip(tip, alpha))
                assert (beta[1]+1)*(beta[2]+1) <= r
                k = alpha[0]
                anchor = (a*tip[1]-q[1]*k, a*tip[2]-q[2]*k)
                assert project(alpha) == tuple(c+a*b for c, b in zip(anchor, beta[1:]))
                counts['independent_actual_projected_residual_labels'] += 1
        rectangles = 6*h*r.bit_length()
        bound = 4*(1+rectangles*(4*r+3*ceiling_cube_root(r**3*a*a)))
        assert len(walls) <= bound
        rows.append({'direction': q, 'cost': threshold(q), 'residual_rank': r,
                     'actual_parallel_pieces': len(walls), 'explicit_shadow_bound': bound,
                     'projected_edge_vertices': sorted(shadow_points)})
        counts['independent_actual_parallel_shadow_covers'] += 1
    return {'length': n, 'coordinate_degree': h, 'initial_ideals': len(cells),
            'undirected_walls': len(edges), 'zero_component_cost': cost_zero,
            'opposite_sign_cost': cost_opp, 'remaining_full_support_cost': cost_remaining,
            'directions': rows}


def odd_and_even_basis(n, B, counts):
    s = (n+1)//2
    C = n*(1+n+n*n)
    lambdas = (C*s+1, C+n, C+n*n)
    columns = sorted((s*a+b+c, dot(lambdas, (a, b, c))) for a, b, c in B)
    assert all(v <= i for i, (v, _) in enumerate(columns))
    assert all(x < y for (_, x), (_, y) in zip(columns, columns[1:]))
    counts['independent_all_length_ordered_Hermite_conditions'] += 1
    if n > 12:
        return None
    A = [[0 if row < v else factorial(row)//factorial(row-v)*xi**(row-v)
          for v, xi in columns] for row in range(n)]
    det = bareiss_det(A)
    assert det
    counts['independent_actual_Q_determinants'] += 1
    if n % 2:
        counts['independent_odd_length_Q_determinants'] += 1
    return {'nonzero': True, 'absolute_bit_length': abs(det).bit_length(),
            'sha256': hashlib.sha256(str(det).encode()).hexdigest()}


def minimal_actual_fixtures(counts):
    """Independent monomial quotients check the rank-one residual boundary."""
    examples = [
        (3, [frozenset(((0, 0, 0), (0, 1, 0), (1, 0, 0))),
             frozenset(((0, 0, 0), (0, 1, 0), (0, 2, 0)))],
         lambda a: 0 if a[2] or 2*a[0]+a[1] >= 3 else 1 << (2*a[0]+a[1]),
         'Q[y]/(y^3), x=y^2, z=0'),
        (4, [frozenset(((0, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 0))),
             frozenset(((0, 0, 0), (0, 1, 0), (0, 0, 1), (0, 1, 1)))],
         lambda a: 0 if max(a[0]+a[1], a[0]+a[2]) >= 2 else
                   1 << (2*(a[0]+a[1])+a[0]+a[2]),
         'Q[y,z]/(y^2,z^2), x=yz'),
    ]
    out = []
    for n, bases, ev, model in examples:
        for B in bases:
            columns = [ev(a) for a in sorted(B)]
            assert len(set(columns)) == n and all(c and c & (c-1) == 0 for c in columns)
            assert downset(B)
            counts['independent_minimal_actual_monomial_bases'] += 1
        edges = {frozenset(bases)}
        row = actual_cover_checks(n, {B: None for B in bases}, edges, 2, counts)
        assert row['directions'][0]['residual_rank'] == 1
        row['actual_monomial_model'] = model
        out.append(row)
    return out


def planar_pairs(M, counts):
    E = tuple(sorted(hyperbolic_cross(M+1, 2)))
    events = sorted({Q(d-b, a-c) for (a, b), (c, d) in combinations(E, 2)
                     if (d-b)*(a-c) > 0})
    bounds = [Q(0)]+events+[None]
    pairs, prefixes = {}, defaultdict(set)
    for low, high in zip(bounds, bounds[1:]):
        rho = low+1 if high is None else (low+high)/2
        w = (rho.numerator, rho.denominator)
        ordered = sorted(E, key=lambda a: (dot(w, a), a))[:M+1]
        B = set()
        for m, v in enumerate(ordered):
            if m:
                key = (frozenset(B), v)
                pairs.setdefault(key, w)
                prefixes[m].add(frozenset(B))
            B.add(v)
        counts['independent_complete_planar_order_intervals'] += 1
    assert len(pairs) >= sum(map(len, prefixes.values()))
    axes = [(B, v) for B, v in pairs if not all(v)]
    assert len(axes) <= M*(M+1)
    for coordinate in (0, 1):
        for a in range(1, M+1):
            values = {B for B, v in axes if v[coordinate] == a and v[1-coordinate] == 0}
            assert len(values) <= M-a+1
            seq = sorted(values, key=len)
            assert all(B < D for B, D in zip(seq, seq[1:]))
            counts['independent_nested_axis_next_label_bounds'] += 1
    rows, by_length, seen_walls = [], defaultdict(set), set()
    for (B, v), w in sorted(pairs.items(), key=lambda item: (len(item[0][0]), sorted(item[0][0]), item[0][1])):
        if not all(v):
            continue
        m = len(B)
        axis = sum(not all(p) for p in B)
        n = 2*m-axis+1
        assert m+1 <= n <= 2*m+1 and n >= 4
        full = {(0, i, j) for i, j in B}
        full.update((1, i-1, j-1) for i, j in B if i and j)
        i, j = v
        tips = {(0, i, j), (1, i-1, j-1)}
        expected = {frozenset(full | {tip}) for tip in tips}
        center = (w[0]+w[1], *w)
        scale = 4*n*n
        weights = [(scale*center[0]+sign, scale*center[1], scale*center[2])
                   for sign in (-1, 1)]
        endpoints = {slab_prefix(n, weight) for weight in weights}
        assert endpoints == expected
        assert all(downset(A) for A in endpoints)
        A, D = sorted(endpoints, key=lambda A: sorted(A))
        assert primitive(tuple(a-b for a, b in zip(profile(A), profile(D)))) == (1, -1, -1)
        key = frozenset(endpoints)
        assert key not in seen_walls
        seen_walls.add(key)
        by_length[n].add((B, v))
        dets = [odd_and_even_basis(n, A, counts), odd_and_even_basis(n, D, counts)]
        rows.append({'planar_prefix_size': m, 'planar_prefix': sorted(B), 'next_label': v,
                     'length': n, 'exact_wall_weight': center,
                     'endpoint_profiles': [profile(A), profile(D)], 'determinants': dets})
        counts['independent_injective_fixed_direction_wall_lifts'] += 1
    assert len(seen_walls) == len(rows) == len(pairs)-len(axes)
    # A separate weighted projected-prefix sweep captures the whole fixed-q fan
    # at lengths <=M, including intervals whose last label is not duplicated.
    direct = defaultdict(set)
    for low, high in zip(bounds, bounds[1:]):
        rho = low+1 if high is None else (low+high)/2
        ordered = sorted(E, key=lambda a: (rho*a[0]+a[1], a))
        B, used = set(), 0
        for v in ordered:
            copies = 2 if all(v) else 1
            if copies == 2 and 4 <= used+1 <= M:
                direct[used+1].add((frozenset(B), v))
            used += copies
            B.add(v)
            if used > M:
                break
    for n in range(4, M+1):
        assert direct[n] == by_length[n]
        counts['independent_complete_fixed_direction_fans'] += 1
    return {'maximum_planar_rank': M, 'ordinary_corner_cut_counts':
            {m: len(B) for m, B in sorted(prefixes.items())},
            'all_prefix_next_pairs': len(pairs), 'axis_pairs': len(axes),
            'axis_upper_bound': M*(M+1), 'interior_pairs': len(rows),
            'parallel_piece_counts': {n: len(B) for n, B in sorted(by_length.items())},
            'complete_parallel_lengths': list(range(4, M+1)), 'wall_lifts': rows}


def main():
    counts = Counter()
    translated = cover_sweeps(counts)
    fans = minimal_actual_fixtures(counts)
    for t in (2, 3):
        n, cells, edges, metadata = complete_hook(t, counts)
        row = actual_cover_checks(n, cells, edges, 2, counts)
        row.update(metadata)
        fans.append(row)
        print('residual covers of complete coupled hook:', t, n, flush=True)
    saved = json.loads((ROOT/'results/research/capped_jet_checks.json').read_text())
    assert saved['all_checks_passed']
    from run_coordinate_slice_budget_checks import edges_from_polygons
    for fixture in saved['complete_quotient_fans']:
        cells = {frozenset(map(tuple, C['standard_set'])):
                 [tuple(map(Q, v)) for v in C['positive_cone']]
                 for C in fixture['all_cells']}
        edges = edges_from_polygons(cells, counts)
        row = actual_cover_checks(fixture['length'], cells, edges, fixture['jet_order'], counts)
        row.update(characteristic=fixture['characteristic'], jet_order=fixture['jet_order'])
        fans.append(row)
    lower = planar_pairs(24, counts)
    sources = ['docs/bounded_coordinate_wall_shadows.md',
               'scripts/run_bounded_coordinate_wall_checks.py',
               'docs/coordinate_slice_wall_budget.md', 'docs/planar_dyadic_bound.md',
               'docs/wall_coexistence.md', 'docs/nested_jet_matroids.md',
               'scripts/run_coordinate_slice_budget_checks.py',
               'scripts/run_planar_dyadic_bound_checks.py',
               'scripts/run_wall_coexistence_checks.py',
               'scripts/run_primitive_wall_checks.py',
               'scripts/run_structural_complexity_checks.py',
               'results/research/capped_jet_checks.json']
    out = {'schema_version': 1, 'all_checks_passed': True,
           'scope': 'Uniform bounded-coordinate parallel and singleton-sign cost theorems, '
                    'with exact finite translated-cover, actual-fan and fixed-direction controls. '
                    'The fixed-q lower order is proved by an infinite averaging argument, '
                    'not fitted to the finite counts. No full spatial budget, unrestricted '
                    'G3 improvement, polynomially separating family or novelty established.',
           'uniform_results': {
               'parallel': 'p_I(q)=O(h*(1+abs(q_x)^(2/3))*r*log(r+1)) for q_x nonzero',
               'singleton_sign_budget': 'V_x^opp=O(h^(8/3)*n^3*log(n+1)^2)',
               'fixed_direction_sharpness': 'For q=(1,-1,-1), max over N/2<=n<=2N+1 '
                                           'of p_(I_n)(q)=Theta(N*log(N+1)); x^2=0 '
                                           'over Q, also radical two-plane realizations.',
               'remaining': 'q=(a,b,-c) or (a,-b,c), 1<=a<h, b,c>0; unbounded coordinate degrees.'},
           'counters': dict(counts), 'translated_cross_fixtures': translated,
           'actual_complete_fan_residual_covers': fans, 'fixed_direction_lower_controls': lower,
           'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}}
    target = ROOT/'results/research/bounded_coordinate_wall_checks.json'
    target.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'all_checks_passed': True, 'counters': dict(counts),
                      'lower_pieces': lower['interior_pairs'],
                      'report': str(target.relative_to(ROOT))}), flush=True)


if __name__ == '__main__':
    main()
