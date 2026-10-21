#!/usr/bin/env python3
"""Exact controls for translated-layer shadows and the sign-class budget."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import gcd, prod
from pathlib import Path
from random import Random
import hashlib
import json

from run_planar_dyadic_bound_checks import convex_hull, matrix_oracle
from run_primitive_wall_checks import threshold
from run_projection_tag_checks import (compare_fans, fan, hit_formula, point_model,
                                      predicted_fan, projection, tag)
from run_rational_shadow_bound_checks import (det2, direction_control, project,
                                             projection_data, sweep)
from run_structural_complexity_checks import (area2, cyclic_univariate_model, dot,
                                              hyperbolic_cross, univariate_power_mod)

ROOT = Path(__file__).resolve().parents[1]


def primitive(v):
    d = gcd(*v)
    q = tuple(a//d for a in v)
    return q if next(a for a in q if a) > 0 else tuple(-a for a in q)


def cube_ceiling(n):
    lo, hi = 0, 1
    while hi**3 < n:
        hi *= 2
    while lo+1 < hi:
        mid = (lo+hi)//2
        if mid**3 >= n:
            hi = mid
        else:
            lo = mid
    return hi if n else 0


def layer_cover(q, H, h):
    s = q[0]
    assert s > 0 and gcd(*q) == 1
    boxes = []
    for j in range(min(h, H)):
        Hj = H//(j+1)
        for k in range(Hj.bit_length()):
            By, Bz = 2**(k+1), Hj//2**k
            assert By*Bz <= 2*Hj
            center = (-j*q[1], -j*q[2])
            low = tuple(2*c-1 for c in center)
            high = tuple(2*(c+s*b)-1 for c, b in zip(center, (By, Bz)))
            widths = (s*By, s*Bz)
            assert prod(widths) <= 2*s*s*Hj
            boxes.append({'layer': j, 'H_layer': Hj, 'lower_twice': low,
                          'upper_twice': high, 'widths': widths})
    return boxes


def contains(box, v):
    return all(a < 2*c < b for a, b, c in
               zip(box['lower_twice'], box['upper_twice'], v))


def replay_layers(original, q, r, H, h, trace, counts):
    s = q[0]
    A = ((-q[1], s, 0), (-q[2], 0, s))
    assert projection_data(A)[0] == s
    ground = tuple(project(A, a) for a in original)
    boxes = layer_cover(q, H, h)
    lines = sorted({(i, b[key][i]) for b in boxes for i in range(2)
                    for key in ('lower_twice', 'upper_twice')})
    assert len(boxes) <= min(h, H)*H.bit_length()
    assert len(lines) <= 4*min(h, H)*H.bit_length()
    assert all(2*v[i] != level for v in ground for i, level in lines)
    labels = [next(i for i, b in enumerate(boxes) if contains(b, v)) for v in ground]
    counts['independent_translated_cover_assignments'] += len(ground)
    counts['overlapping_origin_layer_labels'] += sum(boxes[label]['layer'] != a[0]
                                                    for a, label in zip(original, labels))
    local = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    variation = defaultdict(lambda: [0, 0])
    cross_moves = 0
    for wall in trace['walls']:
        t = Q(wall['parameter'])
        for move in wall['pivots']:
            e, f, token = move['from'], move['to'], move['token']
            step = tuple(b-a for a, b in zip(ground[e], ground[f]))
            left, right = labels[e], labels[f]
            if left != right:
                assert any((2*ground[e][i] < level) != (2*ground[f][i] < level)
                           for i, level in lines)
                cross_moves += 1
            else:
                for i in range(2):
                    local[left][t][i] += step[i]
                    variation[token, left][i] += abs(step[i])
    boundary_crossings = 0
    for path in trace['token_paths']:
        for i, level in lines:
            sides = [2*v[i] < level for v in path]
            moves = sum(a != b for a, b in zip(sides, sides[1:]))
            assert moves <= 2
            boundary_crossings += moves
            counts['independent_token_boundary_budgets'] += 1
    assert cross_moves <= boundary_crossings <= 2*r*len(lines)
    for token, label in variation:
        assert all(v <= 2*w for v, w in zip(variation[token, label], boxes[label]['widths']))
        counts['independent_token_rectangle_budgets'] += 1
    certificates = []
    for label, rows in sorted(local.items()):
        vectors = [v for _, v in sorted(rows.items())]
        cert = direction_control(vectors, s)
        assert all(v <= 2*r*w for v, w in zip(cert['absolute_budgets'], boxes[label]['widths']))
        assert len(vectors) <= 6+cube_ceiling(3456*r*r*s*boxes[label]['H_layer'])
        counts['independent_global_lattice_direction_budgets'] += 1
        certificates.append({'rectangle': label, 'directions': len(vectors),
                             'coordinate_budgets': cert['absolute_budgets']})
    assert len(trace['walls']) <= cross_moves+sum(len(v) for v in local.values())
    explicit_bound = 1+8*r*min(h, H)*H.bit_length()+sum(
        (H//(j+1)).bit_length()*(6+cube_ceiling(3456*r*r*s*(H//(j+1))))
        for j in range(min(h, H)))
    assert len(trace['states']) <= explicit_bound
    counts['complete_quadrant_layer_sweeps'] += 1
    return {'states': len(trace['states']), 'cross_label_exchanges': cross_moves,
            'boundary_crossings': boundary_crossings, 'explicit_quadrant_bound': explicit_bound,
            'internal_direction_certificates': certificates}


def projection_controls(counts):
    models = []
    for H, h, r in ((1, 1, 0), (1, 1, 1), (3, 2, 2), (4, 2, 3), (6, 3, 5)):
        E = tuple(a for a in sorted(hyperbolic_cross(H, 3)) if a[0] < h)
        models.append((f'uniform_{H}_{h}_{r}', H, h, r, E, lambda S, r=r: len(S) <= r))
        rng = Random(17*H+3*r)
        columns = [tuple(int(i == j) for i in range(r)) if j < r
                   else tuple(rng.randrange(101) for _ in range(r)) for j in range(len(E))]
        models.append((f'linear_{H}_{h}_{r}', H, h, r, E, matrix_oracle(columns, r, 101)))
    kernels = ((1, 0, 0), (1, 1, -1), (1, 1009, -1000003),
               (2, 1, -3), (3, -2, 5), (0, 1, -1009))
    ledger, witnesses = [], []
    for name, H, h, r, E, independent in models:
        for q in kernels:
            A = ((-q[1], q[0], 0), (-q[2], 0, q[0])) if q[0] else (
                (1, 0, 0), (0, -q[2], q[1]))
            rows = [sweep(E, A, r, independent, H, signs) for signs in product((-1, 1), repeat=2)]
            exposed = {tuple(state['profile']) for row in rows for state in row['states']}
            if len(E) <= 12:
                profiles = {tuple(sum(project(A, E[i])[j] for i in B) for j in range(2))
                            for B in combinations(range(len(E)), r) if independent(frozenset(B))}
                assert exposed == set(convex_hull(profiles))
                counts['independent_complete_basis_profile_hulls'] += 1
            if q[0]:
                controls = [replay_layers(E, q, r, H, h, row, counts) for row in rows]
            else:
                assert len(exposed) <= 2*(r*(h-1)+1)
                counts['independent_zero_component_shadow_bounds'] += 1
                controls = [{'states': len(row['states'])} for row in rows]
            ledger.append({'model': name, 'kernel': q, 'shadow_vertices': len(exposed),
                           'quadrants': controls})
            if name == 'linear_4_2_3' and q in ((1, 0, 0), (1, 1009, -1000003), (2, 1, -3)):
                witnesses.append({'model': name, 'kernel': q, 'ground': E,
                                  'quadrant_traces': rows})
    return ledger, witnesses


def geometry_guards(counts):
    q, H, h = (1, 1, -1), 128, 2
    E = ((0, 4, 6), (1, 6, 4))
    A = ((-1, 1, 0), (1, 0, 1))
    points = tuple(project(A, a) for a in E)
    boxes = layer_cover(q, H, h)
    lines = {(i, b[key][i]) for b in boxes for i in range(2)
             for key in ('lower_twice', 'upper_twice')}
    assert E[0][0] != E[1][0] and all(prod(c+1 for c in a) <= H for a in E)
    assert all((2*points[0][i] < level) == (2*points[1][i] < level) for i, level in lines)
    labels = [next(j for j, b in enumerate(boxes) if contains(b, v)) for v in points]
    assert labels[0] == labels[1]
    # The determinant two proves that mixed-layer points need the global
    # covolume two, rather than the individual-layer covolume four.
    lattice_vectors = ((2, 0), (-1, 1))
    assert det2(*lattice_vectors) == 2 and det2(*lattice_vectors) % 4
    counts['essential_layer_label_and_covolume_guards'] += 2
    return {'original_layer_label_failure': {'ground': E, 'projection': A,
             'projected_points': points, 'geometric_labels': labels,
             'crossed_cover_boundaries': 0},
            'global_lattice_guard': {'kernel': [2, 1, -1],
             'vectors': lattice_vectors, 'determinant': 2, 'layer_covolume': 4}}


def actual_models(counts):
    models = []
    for n, p in ((4, 2), (5, 3)):
        exponent = (n+1)//2
        images = ((1, 1), tuple(int(j == exponent) for j in range(exponent+1)), (1, 0, 1, 1))
        F = (0,)*n+(1,)
        E = tuple(hyperbolic_cross(n, 3))
        original = cyclic_univariate_model(F, images, E)
        @lru_cache(None)
        def ev(a, n=n, p=p, exponent=exponent):
            if a[0] >= 2:
                return {}
            shift = a[0]*exponent
            left = univariate_power_mod((1, 1), a[1], (0,)*n+(1,))
            right = univariate_power_mod((1, 0, 1, 1), a[2], (0,)*n+(1,))
            out = [0]*n
            for i, c in enumerate(left):
                for j, d in enumerate(right):
                    if i+j+shift < n:
                        out[i+j+shift] += c*d
            return {i: c % p for i, c in enumerate(out) if c % p}
        for a in E:
            assert ev((a[1], a[0], a[2])) == {j: c % p for j, c in enumerate(original[a]) if c % p}
            counts['independent_companion_vs_convolution_columns'] += 1
        models.append((f'nilpotent_{n}_char_{p}', n, 2, p, ev))
    points = ((0, 1, 0), (0, 0, 1), (1, 2, 3), (1, 3, 5))
    models.append(('radical_two_planes', len(points), 2, 101, point_model(points, 101)[1]))
    points = ((0, 1, 0), (0, 0, 1), (1, 2, 3), (1, 3, 5), (2, 4, 7), (2, 5, 2))
    models.append(('radical_three_planes', len(points), 3, 101, point_model(points, 101)[1]))
    return models


def quotient_controls(counts):
    fixtures = []
    for name, n, h, p, ev in actual_models(counts):
        cells = fan(n, ev, p, counts)
        moments = {B: tuple(sum(a[i] for a in B) for i in range(3)) for B in cells}
        E = tuple(a for a in sorted(hyperbolic_cross(n, 3)) if a[0] < h)
        columns = [tuple(ev(a).get(i, 0) for i in range(n)) for a in E]
        independent = matrix_oracle(columns, n, p)
        owners = defaultdict(list)
        for B, poly in cells.items():
            assert len(B) == n and max(a[0] for a in B) < h
            for a, b in zip(poly, poly[1:]+poly[:1]):
                if not any(a[i] == b[i] == 0 for i in range(3)):
                    owners[tuple(sorted((a, b)))].append(B)
        parallel = Counter()
        costs = defaultdict(int)
        for segment, pair in owners.items():
            assert len(pair) == 2
            q = primitive(tuple(a-b for a, b in zip(moments[pair[0]], moments[pair[1]])))
            assert all(dot(q, v) == 0 for v in segment)
            assert abs(q[0]) < h and threshold(q) <= n
            parallel[q] += 1
            kind = 'zero_x' if not q[0] else 'same_yz' if q[1]*q[2] >= 0 else 'remaining_opposite_yz'
            costs[kind] += threshold(q)
            counts['independent_actual_wall_direction_replays'] += 1
        for q, multiplicity in parallel.items():
            if not q[0]:
                continue
            A = ((-q[1], q[0], 0), (-q[2], 0, q[0]))
            rows = [sweep(E, A, n, independent, n, signs) for signs in product((-1, 1), repeat=2)]
            vertices = {tuple(s['profile']) for row in rows for s in row['states']}
            assert multiplicity <= 2*len(vertices)
            for row in rows:
                replay_layers(E, q, n, n, h, row, counts)
        assert costs['zero_x'] <= (n*(h-1)+1)*n*(n-1)
        # Retain multiplicities when combining the individual shadow bound
        # with the independently enumerated singleton-sign cost envelope.
        shadow_envelope = 4*(1+8*n*min(h, n)*n.bit_length()+sum(
            (n//(j+1)).bit_length()*(6+cube_ceiling(3456*n*n*(h-1)*(n//(j+1))))
            for j in range(min(h, n))))
        D = sum(n//u for u in range(1, n+1))
        assert costs['same_yz'] <= shadow_envelope*n*(h-1)*D
        tag_rows = []
        if n == 4:
            for axis in range(3):
                m, planar, embed = projection(n, ev, p, axis, counts)
                qs, hits = hit_formula(cells, planar, embed)
                for q, hit in zip(qs, hits):
                    assert abs(q[0]) < h
                    if not q[0]:
                        assert hit <= n*(h-1)+1
                    else:
                        q = primitive(q)
                        A = ((-q[1], q[0], 0), (-q[2], 0, q[0]))
                        rows = [sweep(E, A, n, independent, n, signs)
                                for signs in product((-1, 1), repeat=2)]
                        vertices = {tuple(s['profile']) for row in rows for s in row['states']}
                        assert hit <= 2*len(vertices)
                        for row in rows:
                            replay_layers(E, q, n, n, h, row, counts)
                        counts['independent_capped_tag_shadow_bounds'] += 1
                    counts['independent_capped_tag_pencil_hits'] += 1
                N, tagged_ev, _, _ = tag(n, ev, p, axis, counts)
                tagged = fan(N, tagged_ev, p, counts)
                compare_fans(tagged, predicted_fan(cells, planar, embed, axis), counts)
                assert len(tagged) == len(cells)+sum(hits)
                tag_rows.append({'axis': axis, 'projection_length': m,
                                 'gain': sum(hits), 'pencil_directions': qs, 'hits': hits})
        fixtures.append({'model': name, 'length': n, 'coordinate_cap': h,
                         'characteristic': p, 'complete_fan_count': len(cells),
                         'undirected_wall_pieces': len(owners), 'costs_by_family': dict(costs),
                         'tag_checks': tag_rows,
                         'complete_cones': [{'basis': sorted(B),
                            'cone': [list(map(str, v)) for v in poly]} for B, poly in cells.items()]})
    return fixtures


def direction_sum_controls(counts):
    ledger = []
    for n in range(2, 65):
        D = sum(n//u for u in range(1, n+1))
        for h in (2, 3, 4):
            total = number = 0
            for c in range(1, min(h, n)):
                for a in range(n-c+1):
                    for b in range((n-c+1)//(a+1)):
                        if a == b == 0 or gcd(c, a, b) != 1:
                            continue
                        q = (c, -a, -b)
                        cost = threshold(q)
                        assert cost == c+(a+1)*(b+1)-1 <= n
                        total += cost
                        number += 1
                        counts['independent_isolated_sign_direction_costs'] += 1
            assert number <= (h-1)*D and total <= n*(h-1)*D
            ledger.append([n, h, number, total, D])
    return ledger


def opposite_sign_envelopes(counts):
    ledger = []
    for n in range(2, 97):
        number = residual_cost_sum = 0
        interior = []
        for a in range(1, n):
            for b in range(1, n-2*a):
                q = (1, a, -b)
                cost = 2*a+b+1
                assert gcd(*q) == 1 and threshold(q) == cost <= n
                residual = n-cost+1
                assert residual >= 1
                number += 1
                residual_cost_sum += cost*residual
                if ((n+15)//16 <= a <= n//8 and
                        (n+3)//4 <= b <= n//3):
                    assert 3*cost >= n and 3*residual >= n
                    interior.append((a, b))
                counts['independent_opposite_sign_residual_envelopes'] += 1
        # Count by a, independently of threshold evaluation and b iteration.
        assert number == sum(max(0, n-2*a-1) for a in range(1, n))
        assert residual_cost_sum <= n*n*(n+1)*(n+1)//4
        if n >= 32:
            assert len(interior)*768 >= n*n
            assert residual_cost_sum*6912 >= n**4
        ledger.append({'length': n, 'primitive_directions': number,
                       'sum_cost_times_residual_rank': residual_cost_sum,
                       'middle_cost_box_directions': len(interior)})
    return ledger


def main():
    counts = defaultdict(int)
    sweeps, witnesses = projection_controls(counts)
    guards = geometry_guards(counts)
    quotients = quotient_controls(counts)
    direction_sums = direction_sum_controls(counts)
    opposite_envelopes = opposite_sign_envelopes(counts)
    sources = ['docs/layer_shadow_budget.md', 'scripts/run_layer_shadow_checks.py',
               'scripts/run_rational_shadow_bound_checks.py', 'scripts/run_planar_dyadic_bound_checks.py',
               'scripts/run_structural_complexity_checks.py', 'scripts/run_projection_tag_checks.py',
               'scripts/run_primitive_wall_checks.py', 'scripts/run_signed_face_checks.py',
               'scripts/run_oracle_enumeration_checks.py']
    report = {'command': 'python3 -B scripts/run_layer_shadow_checks.py',
              'counts': dict(sorted(counts.items())), 'projection_sweeps': sweeps,
              'exact_trace_witnesses': witnesses, 'essential_scope_guards': guards,
              'actual_quotient_fixtures': quotients, 'isolated_sign_direction_sums': direction_sums,
              'opposite_sign_residual_envelopes': opposite_envelopes,
              'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources},
              'unrestricted_G3_improved': False, 'full_collective_budget_proved': False,
              'polynomial_lower_separation_proved': False, 'novelty_established': False,
              'goal_complete': False}
    (ROOT/'results/research/layer_shadow_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'counts': report['counts'], 'quotients': [{k: f[k] for k in
                     ('model', 'length', 'complete_fan_count', 'costs_by_family')} for f in quotients]}, indent=2))


if __name__ == '__main__':
    main()
