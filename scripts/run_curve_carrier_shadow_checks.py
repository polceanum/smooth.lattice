#!/usr/bin/env python3
"""Exact finite-line token budgets and complete fixed-curve quotient fans."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import comb
from pathlib import Path
from random import Random
import hashlib
import json

from run_coordinate_slice_budget_checks import edges_from_polygons, primitive, profile
from run_planar_dyadic_bound_checks import convex_hull, matrix_oracle
from run_projection_tag_checks import (compare_fans, fan, hit_formula, predicted_fan,
                                      tag)
from run_quadric_cross_checks import apply, normalized, rank
from run_rational_shadow_bound_checks import kernel_projection, project
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import (area2, clip, dot, hyperbolic_cross,
                                             univariate_product_mod,
                                             univariate_power_mod)

ROOT = Path(__file__).resolve().parents[1]
TRIANGLE = tuple(tuple(Q(i == j) for i in range(3)) for j in range(3))


def line_cover(ground, source_lines, A):
    """Strips contain all and only the finite ground points on their line."""
    covers = []
    for start, direction in source_lines:
        center, vector = project(A, start), project(A, direction)
        if any(vector):
            normal = (-vector[1], vector[0])
            level = dot(normal, center)
            distances = [abs(dot(normal, v)-level) for v in ground
                         if dot(normal, v) != level]
            eps = Q(min(distances), 2) if distances else Q(1)
            covers.append({'point': center, 'direction': vector,
                           'lines': [(normal, level-eps), (normal, level+eps)]})
        else:
            distances = [abs(v[i]-center[i]) for v in ground for i in range(2)
                         if v[i] != center[i]]
            eps = Q(min(distances), 2) if distances else Q(1)
            covers.append({'point': center, 'direction': vector,
                           'lines': [(tuple(int(i == j) for i in range(2)), center[j]+s*eps)
                                     for j in range(2) for s in (-1, 1)]})

    def contains(c, v):
        return all(dot(row, v) > low and dot(row, v) < high
                   for (row, low), (_, high) in zip(c['lines'][::2], c['lines'][1::2]))

    labels = [next(j for j, c in enumerate(covers) if contains(c, v)) for v in ground]
    for v, j in zip(ground, labels):
        c = covers[j]
        if any(c['direction']):
            a, b = (x-y for x, y in zip(v, c['point']))
            assert a*c['direction'][1] == b*c['direction'][0]
        else:
            assert v == c['point']
    lines = [line for c in covers for line in c['lines']]
    assert all(dot(row, v) != level for row, level in lines for v in ground)
    assert len(lines) <= 4*len(source_lines)
    return covers, labels, lines


def sweep(ground, r, oracle, source_lines, A, signs, counts):
    points = tuple(project(A, a) for a in ground)
    covers, labels, boundaries = line_cover(points, source_lines, A)
    sigma, tau = signs

    def greedy(t):
        B = set()
        for j in sorted(range(len(points)),
                        key=lambda j: (sigma*t*points[j][0]+tau*points[j][1], j)):
            if oracle(frozenset(B|{j})):
                B.add(j)
        assert len(B) == r
        return B

    events = defaultdict(set)
    for i, j in combinations(range(len(points)), 2):
        dx, dy = (a-b for a, b in zip(points[i], points[j]))
        if sigma*tau*dx*dy < 0:
            events[-Q(tau*dy, sigma*dx)].update((i, j))
    events = sorted(events.items())
    first = events[0][0]/2 if events else Q(1)
    B = greedy(first)
    order = sorted(range(len(points)), key=lambda j: (tau*points[j][1], sigma*points[j][0], j))
    positions = {j: i for i, j in enumerate(order)}
    paths = [[points[j]] for j in sorted(B)]
    tokens = {j: i for i, j in enumerate(sorted(B))}
    cross_moves, internal, walls, states = 0, defaultdict(set), [], []

    def record(t):
        assert B == greedy(t)
        value = tuple(sum(points[j][i] for j in B) for i in range(2))
        if not states or states[-1]['profile'] != value:
            states.append({'parameter': str(t), 'profile': value, 'basis': sorted(B)})
        counts['independent_direct_greedy_interval_replays'] += 1

    record(first)
    for k, (t, involved) in enumerate(events):
        groups, pivots = defaultdict(list), []
        for j in involved:
            groups[sigma*t*points[j][0]+tau*points[j][1]].append(j)
        for group in groups.values():
            slots = sorted(positions[j] for j in group)
            assert slots == list(range(slots[0], slots[-1]+1))
            for target, f in enumerate(sorted(group, key=lambda j: (sigma*points[j][0], j)), slots[0]):
                while positions[f] > target:
                    pos = positions[f]
                    e = order[pos-1]
                    if e in B and f not in B and oracle(frozenset(B-{e}|{f})):
                        step = tuple(b-a for a, b in zip(points[e], points[f]))
                        assert -sigma*step[0] > 0 and tau*step[1] > 0
                        assert Q(tau*step[1], -sigma*step[0]) == t
                        token = tokens.pop(e)
                        tokens[f] = token
                        paths[token].append(points[f])
                        if labels[e] != labels[f]:
                            assert any((dot(row, points[e]) < level) !=
                                       (dot(row, points[f]) < level)
                                       for row, level in boundaries)
                            cross_moves += 1
                        else:
                            j = labels[e]
                            v = covers[j]['direction']
                            assert any(v) and step[0]*v[1] == step[1]*v[0]
                            internal[j].add(t)
                        B.remove(e); B.add(f)
                        pivots.append((e, f, token))
                    order[pos-1], order[pos] = f, e
                    positions[f], positions[e] = pos-1, pos
        if pivots:
            walls.append({'parameter': str(t), 'pivots': pivots})
        following = events[k+1][0] if k+1 < len(events) else t+2
        record((t+following)/2)
    crossing_count = 0
    for path in paths:
        for row, level in boundaries:
            values = [dot(row, v)-level for v in path]
            moves = sum((a < 0) != (b < 0) for a, b in zip(values, values[1:]))
            assert moves <= 2
            crossing_count += moves
            counts['independent_token_line_crossing_budgets'] += 1
    assert all(len(v) <= 1 for v in internal.values())
    assert cross_moves <= crossing_count <= 2*r*len(boundaries)
    assert len(walls) <= cross_moves+sum(len(v) for v in internal.values())
    assert len(states) == len(walls)+1 <= 1+8*r*len(source_lines)+len(source_lines)
    counts['complete_four_quadrant_component_sweeps'] += 1
    return {'signs': signs, 'states': states, 'walls': walls,
            'token_paths': paths, 'cross_label_exchanges': cross_moves,
            'boundary_crossings': crossing_count,
            'internal_event_counts': {str(j): len(v) for j, v in internal.items()},
            'cover_members': len(covers), 'ground_avoiding_boundary_lines': len(boundaries)}


def full_shadow(ground, r, oracle, lines, q, counts, exhaustive=False):
    A = kernel_projection(q)
    traces = [sweep(ground, r, oracle, lines, A, signs, counts)
              for signs in product((-1, 1), repeat=2)]
    hull = {v['profile'] for trace in traces for v in trace['states']}
    assert len(hull) <= 4*(1+8*r*len(lines)+len(lines))
    if exhaustive:
        profiles = {tuple(sum(project(A, ground[j])[i] for j in B) for i in range(2))
                    for B in combinations(range(len(ground)), r) if oracle(frozenset(B))}
        assert hull == set(convex_hull(profiles))
        counts['independent_exhaustive_basis_profile_hulls'] += 1
    return hull, traces


def arbitrary_matroids(counts):
    zero = (0, 0, 0)
    axes = [(zero, tuple(int(i == j) for i in range(3))) for j in range(3)]
    skew = [((0, 0, 0), (1, 0, 0)), ((0, 2, 1), (1, 0, 0)),
            ((1, 0, 2), (0, 1, 0)), ((2, 1, 0), (0, 0, 1))]
    qs = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1),
          (2, -3, 5), (101, -103, 107), (1000001, -1000000, 1)]
    rows = []
    for name, lines, length, r in [('axes_zero', axes, 3, 0),
                                  ('axes_rank_one', axes, 3, 1),
                                  ('axes_rank_four', axes, 3, 4),
                                  ('skew_rank_five', skew, 3, 5)]:
        ground = tuple(sorted({tuple(a+k*b for a, b in zip(start, direction))
                               for start, direction in lines for k in range(length)}))
        rng = Random(143+r)
        columns = [tuple(int(i == j) for i in range(r)) if j < r else
                   tuple(rng.randrange(101) for _ in range(r)) for j in range(len(ground))]
        for kind, oracle in [('uniform', lambda S, r=r: len(S) <= r),
                             ('linear_F101', matrix_oracle(columns, r, 101))]:
            for q in qs:
                hull, traces = full_shadow(ground, r, oracle, lines, q, counts, True)
                rows.append({'model': name, 'matroid': kind, 'rank': r, 'ground': ground,
                             'line_starts_and_directions': lines, 'kernel': q,
                             'shadow_vertices': sorted(hull), 'quadrants': traces})
    return rows


def reduce_polynomial(poly, generators, heads, w):
    poly, remainder = dict(poly), {}
    while poly:
        a = max(poly, key=lambda a: (dot(w, a), a))
        for g, head in zip(generators, heads):
            if all(i >= j for i, j in zip(a, head)):
                assert g[head] == 1
                coefficient = poly[a]
                shift = tuple(i-j for i, j in zip(a, head))
                for b, c in g.items():
                    target = tuple(i+j for i, j in zip(shift, b))
                    poly[target] = poly.get(target, 0)-coefficient*c
                    if not poly[target]:
                        del poly[target]
                break
        else:
            remainder[a] = poly.pop(a)
    return remainder


def curve_charts(power=False):
    # C=(xy-1,x+y-z); the five coprime-head Groebner bases give every cone.
    x, y, z, one = (1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    linear = {x: 1, y: 1, z: -1}
    data = [
        ('z_xy', (1, 1, 3), ({z: 1, x: -1, y: -1}, {(1, 1, 0): 1, one: -1}),
         (z, (1, 1, 0)), [(one, x), (one, y)]),
        ('x_y2', (4, 2, 1), (linear, {(0, 2, 0): 1, (0, 1, 1): -1, one: 1}),
         (x, (0, 2, 0)), [(one, z), (y, z)]),
        ('x_yz', (4, 1, 2), (linear, {(0, 1, 1): 1, (0, 2, 0): -1, one: -1}),
         (x, (0, 1, 1)), [(one, y), (one, z)]),
        ('y_x2', (2, 4, 1), (linear, {(2, 0, 0): 1, (1, 0, 1): -1, one: 1}),
         (y, (2, 0, 0)), [(one, z), (x, z)]),
        ('y_xz', (1, 4, 2), (linear, {(1, 0, 1): 1, (2, 0, 0): -1, one: -1}),
         (y, (1, 0, 1)), [(one, x), (one, z)])]
    if power:
        data = [
            ('y_xz', (1, 3, 1),
             ({y: 1, (2, 0, 0): -1}, {(1, 0, 1): 1, one: -1}),
             (y, (1, 0, 1)), [(one, x), (one, z)]),
            ('x_yz2', (4, 1, 1),
             ({x: 1, (0, 1, 1): -1}, {(0, 1, 2): 1, one: -1}),
             (x, (0, 1, 2)), [(one, y), (z, y), (one, z)]),
            ('x2_xz_yz', (2, 2, 1),
             ({(2, 0, 0): 1, y: -1}, {(1, 0, 1): 1, one: -1},
              {(0, 1, 1): 1, x: -1}),
             ((2, 0, 0), (1, 0, 1), (0, 1, 1)), [(one, y), (x, y), (one, z)])]
    original = ({y: 1, (2, 0, 0): -1}, {(1, 0, 1): 1, one: -1}) if power else (
        linear, {(1, 1, 0): 1, one: -1})
    out = []
    for name, w, polys, heads, lines in data:
        assert tuple(max(f, key=lambda a: (dot(w, a), a)) for f in polys) == heads
        for f in polys:
            values = defaultdict(int)
            for (a, b, c), coefficient in f.items():
                if power:
                    values[a+2*b-c] += coefficient
                else:
                    for k in range(c+1):
                        values[a-b+2*k-c] += coefficient*comb(c, k)
            assert not any(values.values())
        for i, j in combinations(range(len(polys)), 2):
            common = tuple(max(a, b) for a, b in zip(heads[i], heads[j]))
            S = defaultdict(int)
            for index, sign in ((i, 1), (j, -1)):
                shift = tuple(a-b for a, b in zip(common, heads[index]))
                for a, c in polys[index].items():
                    S[tuple(b+d for b, d in zip(a, shift))] += sign*c
            assert not reduce_polynomial({a:c for a,c in S.items() if c}, polys, heads, w)
        assert all(not reduce_polynomial(g, polys, heads, w) for g in original)
        poly = list(TRIANGLE)
        for f, head in zip(polys, heads):
            for tail in f:
                if tail != head:
                    poly = clip(poly, tuple(a-b for a, b in zip(head, tail)), 1)
        out.append({'name': name, 'heads': heads, 'polynomials': polys, 'lines': lines,
                    'positive_cone': poly})
    assert sum(area2(v['positive_cone']) for v in out) == 1
    return out


def curve_model(n, p, repeated, power=False):
    f = [comb(n, i)*(-1)**(n-i) for i in range(n+1)] if repeated else [-1]+[0]*(n-1)+[1]
    X = tuple(int(i == 1) for i in range(n))
    assert abs(f[0]) == 1
    Y = tuple(-f[i+1]//f[0] for i in range(n))
    Z = tuple(a+b for a, b in zip(X, Y))
    images = (X, univariate_product_mod(X, X, f), Y) if power else (X, Y, Z)

    @lru_cache(None)
    def ev(a):
        v = (1,)+(0,)*(n-1)
        for image, k in zip(images, a):
            v = univariate_product_mod(v, univariate_power_mod(image, k, f), f)
        return normalized(dict(enumerate(v)), p)

    M = [tuple(normalized(dict(enumerate(univariate_product_mod(
         image, tuple(int(i == j) for i in range(n)), f))), p) for j in range(n))
         for image in images]
    return n, ev, M, p, {'length': n, 'field': p, 'modulus': f,
                         'kind': 'repeated_unit' if repeated else 'unit_roots',
                         'curve': 'x=t,y=t^2,z=t^-1' if power else 'x=t,y=t^-1,z=t+t^-1'}


def curve_controls(counts):
    rows = []
    cases = [(n,p,repeated,False) for n,p,repeated in
             [(3,0,True),(4,0,False),(5,0,True),(6,0,False),(4,2,True),(5,3,True),(5,101,False)]]
    cases += [(n,p,repeated,True) for n,p,repeated in
              [(4,0,True),(6,0,False),(4,2,True),(5,3,True),(6,101,False)]]
    for n, p, repeated, power in cases:
        charts = curve_charts(power)
        _, ev, M, _, metadata = curve_model(n, p, repeated, power)
        for j in range(n):
            if power:
                assert apply(M[0], M[2][j], p) == {j: 1}
                assert apply(M[0], M[0][j], p) == M[1][j]
            else:
                assert apply(M[0], M[1][j], p) == {j: 1}
                assert normalized({i: M[0][j].get(i, 0)+M[1][j].get(i, 0)-M[2][j].get(i, 0)
                                   for i in range(n)}, p) == {}
            for a, b in combinations(range(3), 2):
                assert apply(M[a], M[b][j], p) == apply(M[b], M[a][j], p)
            counts['independent_curve_multiplication_columns'] += 1
        actual = fan(n, ev, p, counts)
        edges = edges_from_polygons(actual, counts)
        E = tuple(sorted(hyperbolic_cross(n, 3)))
        grounds, oracles, assignments, patches = [], [], {}, []
        for chart in charts:
            ground = tuple(a for a in E if not any(all(c >= d for c, d in zip(a, head))
                                                  for head in chart['heads']))
            assert rank([ev(a) for a in ground], p) == n
            oracle = lru_cache(None)(lambda S, G=ground: rank([ev(G[j]) for j in S], p) == len(S))
            grounds.append(ground); oracles.append(oracle)
            local = []
            for B, poly in actual.items():
                part = list(poly)
                for f, head in zip(chart['polynomials'], chart['heads']):
                    for tail in f:
                        if tail != head:
                            part = clip(part, tuple(a-b for a, b in zip(head, tail)), 1)
                if area2(part):
                    assert B <= set(ground)
                    local.append((B, part))
                    counts['independent_complete_curve_carrier_patches'] += 1
            patches.append(local)
        for B, poly in actual.items():
            assert sum(area2(part) for local in patches for label, part in local if label == B) == area2(poly)
            assignments[B] = next(j for j, ground in enumerate(grounds) if B <= set(ground))
        qs = {(1, -1, 0), (2, -3, 5), (1000001, -1000000, 1)}
        edge_directions = {}
        for edge in edges:
            B, D = tuple(edge)
            q = primitive(tuple(a-b for a, b in zip(profile(B), profile(D))))
            edge_directions[edge] = q
            qs.add(q)
        ledger = []
        for q in sorted(qs):
            A = kernel_projection(q)
            shadows = [full_shadow(G, n, oracle, chart['lines'], q, counts, n <= 4)[0]
                       for G, oracle, chart in zip(grounds, oracles, charts)]
            cuts, parallel = defaultdict(set), defaultdict(set)
            for B, poly in actual.items():
                if min(dot(q, v) for v in poly) < 0 < max(dot(q, v) for v in poly):
                    j = assignments[B]; value = project(A, profile(B))
                    assert value in shadows[j] and value not in cuts[j]
                    cuts[j].add(value)
                    counts['independent_strict_cone_shadow_injections'] += 1
            for edge, direction in edge_directions.items():
                if direction == primitive(q):
                    B = min(edge, key=lambda B: sorted(B))
                    j = assignments[B]; value = project(A, profile(B))
                    assert value in shadows[j] and value not in parallel[j]
                    parallel[j].add(value)
                    counts['independent_parallel_edge_shadow_injections'] += 1
            ledger.append({'kernel': q, 'cut_cones': sum(map(len, cuts.values())),
                           'parallel_wall_pieces': sum(map(len, parallel.values())),
                           'carrier_shadow_counts': list(map(len, shadows))})
        tags = []
        if n == 3 or (n == 4 and p == 2):
            for axis in range(3):
                N, tagged, cells, embed = tag(n, ev, p, axis, counts)
                expected = predicted_fan(actual, cells, embed, axis, 1)
                child = fan(N, tagged, p, counts)
                compare_fans(child, expected, counts)
                wall_normals, hits = hit_formula(actual, cells, embed)
                assert len(child) == len(actual)+sum(hits)
                assert len(wall_normals) == len(hits)
                assert all(v <= len(charts)*4*(1+8*n*max(len(c['lines']) for c in charts)
                                              +max(len(c['lines']) for c in charts)) for v in hits)
                tags.append({'axis': axis, 'length': N, 'initial_ideals': len(child),
                             'added_cells': len(child)-len(actual)})
                counts['independent_complete_coordinate_tag_fans'] += 1
        rows.append({**metadata, 'initial_ideals': len(actual), 'walls': len(edges),
                     'carrier_patch_counts': list(map(len, patches)), 'shadow_injections': ledger,
                     'coordinate_tags': tags,
                     'full_cones': [{'standard_set': sorted(B),
                                     'positive_cone': [[str(c) for c in v] for v in poly]}
                                    for B, poly in actual.items()]})
        print('complete curve carrier:', n, p, repeated, len(actual), flush=True)
    return rows


def main():
    counts = Counter()
    matroids = arbitrary_matroids(counts)
    actual = curve_controls(counts)
    sources = ['docs/curve_carrier_shadows.md','scripts/run_curve_carrier_shadow_checks.py',
               'docs/rational_shadow_bound.md','docs/projection_tag_fans.md',
               'scripts/run_coordinate_slice_budget_checks.py','scripts/run_planar_dyadic_bound_checks.py',
               'scripts/run_projection_tag_checks.py','scripts/run_quadric_cross_checks.py',
               'scripts/run_rational_shadow_bound_checks.py','scripts/run_signed_face_checks.py',
               'scripts/run_structural_complexity_checks.py']
    report = {'schema_version': 1, 'all_checks_passed': True,
              'scope': 'Finite-line shadow theorem, arbitrary matroids and actual fixed-curve '
                       'positive fans. No sharp full-curve count, unrestricted G3 bound, '
                       'polynomial lower separation or novelty established.',
              'counters': dict(counts), 'arbitrary_matroid_shadows': matroids,
              'complete_curve_quotient_fans': actual,
              'source_sha256': {s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}}
    path = ROOT/'results/research/curve_carrier_shadow_checks.json'
    path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'all_checks_passed': True,'counters':dict(counts),
                      'report': str(path.relative_to(ROOT))}),flush=True)


if __name__ == '__main__':
    main()
