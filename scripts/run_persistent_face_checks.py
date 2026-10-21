#!/usr/bin/env python3
"""Exact common-core and sharp valence checks for persistent two-flag faces."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
from math import gcd
from pathlib import Path
import hashlib
import json

from run_planar_dyadic_bound_checks import convex_hull
from run_primitive_wall_checks import canonical, divisors
from run_projection_tag_checks import compare_fans, fan
from run_quadratic_coordinate_fan_checks import (
    TRIANGLE, compile_fan, quadratic_tensor, socle_crossing, square_zero_fixtures,
)
from run_square_zero_extension_checks import matrix_controls
from run_structural_complexity_checks import area2, clip, dot

ROOT = Path(__file__).resolve().parents[1]


def wall_core(q):
    plus, minus = tuple(max(c, 0) for c in q), tuple(max(-c, 0) for c in q)
    return (divisors(plus)|divisors(minus))-{plus, minus}


def divides(a, b):
    return all(x <= y for x, y in zip(a, b))


def persistent_faces(actual, endpoints, n, counts):
    weights = {w for poly in actual.values() for w in poly if min(w) > 0}
    profiles = {B: tuple(sum(a[i] for a in B) for i in range(3)) for B in actual}
    out = []
    for w in sorted(weights-endpoints):
        best = min(dot(w, v) for v in profiles.values())
        incident = {B: v for B, v in profiles.items() if dot(w, v) == best}
        F = set.intersection(*(set(B) for B in incident))
        hull = convex_hull([(v[0], v[1]) for v in incident.values()])
        assert len(hull) == len(incident) and len(hull) % 2 == 0
        half = len(hull)//2
        assert len({tuple(a+b for a, b in zip(hull[i], hull[(i+half) % len(hull)]))
                    for i in range(len(hull))}) == 1
        by_profile = {(v[0], v[1]): B for B, v in incident.items()}
        groups = defaultdict(list)
        for va, vb in zip(hull, hull[1:]+hull[:1]):
            A, B = by_profile[va], by_profile[vb]
            delta = tuple(b-a for a, b in zip(profiles[A], profiles[B]))
            q = canonical(delta)
            groups[q].append((A, B, delta))
        K = set().union(*(wall_core(q) for q in groups))
        assert K <= F
        counts['simultaneous_compulsory_core_containments'] += len(groups)
        nonvertical = {q: edges for q, edges in groups.items() if q[0]}
        vertical = len(groups)-len(nonvertical)
        assert vertical <= 1
        representatives, total_x = [], 0
        for q, edges in nonvertical.items():
            assert q[0] == 1 and len(edges) == 2
            assert abs(edges[0][2][0]) == abs(edges[1][2][0])
            A, B, delta = edges[0]
            if delta[0] < 0:
                A, B = B, A
            upper_added = [a for a in B-A if a[0] == 1]
            assert upper_added
            representatives.append(min(upper_added))
            total_x += abs(delta[0])
        d = len(nonvertical)
        for a, b in combinations(representatives, 2):
            assert not divides(a, b) and not divides(b, a)
            counts['incomparable_upper_direction_representatives'] += 1
        upper_basis = max(incident, key=lambda B: profiles[B][0])
        assert set(representatives) <= upper_basis
        upper_diagram = {a[1:] for a in upper_basis if a[0] == 1}
        assert len(upper_diagram) >= d*(d+1)//2
        assert all((0, *a) in upper_basis for a in upper_diagram)
        assert d*(d+1) <= n
        width = max(v[0] for v in incident.values())-min(v[0] for v in incident.values())
        assert width == total_x <= n-len(K)
        assert len(incident) == 2*(d+vertical)
        counts['persistent_common_core_width_and_valence_controls'] += 1
        out.append({'weight': list(map(str, w)), 'incident_cones': len(incident),
                    'nonvertical_directions': sorted(nonvertical), 'vertical_directions': vertical,
                    'common_fixed_monomials': sorted(F), 'joint_compulsory_core': sorted(K),
                    'residual_rank': n-len(K), 'x_width': width,
                    'upper_antichain_representatives': representatives,
                    'upper_diagram_size': len(upper_diagram), 'profile_polygon': hull})
    return out


def socle_model(D0, D1, pairs):
    ground = {(0, *a) for a in D0}|{(1, *a) for a in D1}
    upper_tips = {b for a, b in pairs}
    basis = sorted(ground-upper_tips)
    columns = {a: {i: 1} for i, a in enumerate(basis)}
    for a, b in pairs:
        assert a in columns and b not in columns
        columns[b] = columns[a]
    return len(basis), lambda a: dict(columns.get(a, {})), ground


def predicted_socle_fan(ground, pairs):
    cells = [list(TRIANGLE)]
    for a, b in pairs:
        q = tuple(x-y for x, y in zip(a, b))
        cells = [piece for cell in cells for sign in (-1, 1)
                 if (piece := clip(cell, q, sign))]
    fixed = ground-{v for pair in pairs for v in pair}
    out = {}
    for poly in cells:
        w = tuple(sum(v[i] for v in poly)/len(poly) for i in range(3))
        B = frozenset(fixed|{a if dot(w, a) < dot(w, b) else b for a, b in pairs})
        assert B not in out
        out[B] = poly
    assert sum(area2(poly) for poly in out.values()) == 1
    return out


def sharp_family(counts):
    out = []
    for d, char in product(range(2, 9), (0, 2, 3)):
        D1 = {(i, j) for i in range(d) for j in range(d-i)}
        D0 = D1|{(d-i, i) for i in range(d)}
        pairs = [((0, d-i, i), (1, i, d-1-i)) for i in range(d)]
        n, ev, ground = socle_model(D0, D1, pairs)
        assert n == d*(d+1)
        basis = tuple(sorted(ground-{b for a, b in pairs}))
        matrix_controls(n, ev, basis, ground, char, counts)
        actual = fan(n, ev, char, counts)
        compare_fans(actual, predicted_socle_fan(ground, pairs), counts)
        assert len(actual) == 2*d
        rows = persistent_faces(actual, set(), n, counts)
        assert len(rows) == 1 and rows[0]['weight'] == ['1/3']*3
        assert len(rows[0]['nonvertical_directions']) == d
        assert rows[0]['vertical_directions'] == 0
        assert rows[0]['upper_diagram_size'] == d*(d+1)//2
        out.append({'d': d, 'characteristic': char, 'length': n,
                    'complete_fan_cones': len(actual), 'socle_pairs': pairs, 'persistent_face': rows[0]})
    return out


def strict_flag_family(counts):
    out = []
    for b, char in product((3, 4, 6, 7, 9), (0, 2)):
        D0 = {(0, j) for j in range(b+1)}|{(1, 0), (2, 0)}
        D1 = {(0, 0), (1, 0), (0, 1)}
        pairs = [((0, 0, b), (1, 1, 0)), ((0, 2, 0), (1, 0, 1))]
        n, ev, ground = socle_model(D0, D1, pairs)
        assert n == b+4 and gcd(b+1, 3) == 1
        basis = tuple(sorted(ground-{v for a, v in pairs}))
        matrix_controls(n, ev, basis, ground, char, counts)
        actual = fan(n, ev, char, counts)
        compare_fans(actual, predicted_socle_fan(ground, pairs), counts)
        rows = persistent_faces(actual, set(), n, counts)
        assert len(rows) == 1 and rows[0]['residual_rank'] == rows[0]['x_width'] == 2
        u, t = Q(b+1, 3), Q(2*b-1, 3)
        for D in (D0, D1):
            costs = [u*i+j for i, j in D]
            assert len(costs) == len(set(costs))
            counts['strict_planar_flag_crossing_controls'] += 1
        out.append({'b': b, 'characteristic': char, 'length': n,
                    'crossing_u': str(u), 'crossing_t': str(t),
                    'complete_fan_cones': len(actual), 'persistent_face': rows[0]})
    return out


def general_controls(counts):
    models = list(square_zero_fixtures())
    for char, r, c, b in ((2, 3, 1, 1), (2, 3, 1, 0), (3, 2, 1, 0), (0, 3, 1, 0), (0, 2, -1, 0)):
        n, ev = quadratic_tensor(char, r, c, b)
        models.append(({'model': 'quadratic_tensor', 'characteristic': char}, n, ev))
    for char in (0, 2):
        n, ev = socle_crossing()
        models.append(({'model': 'socle_crossing', 'characteristic': char}, n, ev))
    out = []
    for meta, n, ev in models:
        actual = fan(n, ev, meta['characteristic'], counts)
        compiled, _, _, endpoints = compile_fan(n, ev, meta['characteristic'], counts)
        compare_fans(actual, compiled, counts)
        rows = persistent_faces(actual, endpoints, n, counts)
        out.append(dict(meta, length=n, complete_fan_cones=len(actual), persistent_faces=rows))
    return out


def main():
    counts = defaultdict(int)
    controls = general_controls(counts)
    sharp = sharp_family(counts)
    strict = strict_flag_family(counts)
    sources = ['docs/persistent_matching_faces.md', 'scripts/run_persistent_face_checks.py',
               'scripts/run_quadratic_coordinate_fan_checks.py', 'scripts/run_projection_tag_checks.py',
               'scripts/run_primitive_wall_checks.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_planar_dyadic_bound_checks.py', 'scripts/run_square_zero_extension_checks.py']
    report = {'command': 'python3 -B scripts/run_persistent_face_checks.py',
              'counts': dict(sorted(counts.items())), 'general_quotient_controls': controls,
              'sharp_valence_family': sharp, 'crossings_inside_fixed_flags': strict,
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              'unrestricted_G3_improved': False, 'quadratic_class_fan_count_improved': False,
              'collective_cubic_budget_proved': False, 'novelty_established': False, 'goal_complete': False}
    (ROOT/'results/research/persistent_face_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'counts': report['counts'], 'sharp_fans': len(sharp),
                      'strict_flag_fans': len(strict), 'general_fans': len(controls)}, indent=2))


if __name__ == '__main__':
    main()
