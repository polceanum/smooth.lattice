#!/usr/bin/env python3
"""Exact complete fans for the logarithmic cubic wall-cost obstruction."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import gcd
from pathlib import Path
import hashlib
import json

from run_free_wall_checks import first_terms, wall_cost
from run_interleaving_research import minimal_excluded_border
from run_planar_dyadic_bound_checks import convex_hull
from run_signed_face_checks import field_basis, normals as border_normals
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross
from run_wall_coexistence_checks import complete_planar_fan, positive_edges

ROOT = Path(__file__).resolve().parents[1]


def profile(B):
    return tuple(sum(a[i] for a in B) for i in range(len(next(iter(B)))))


def primitive(v):
    d = gcd(*v)
    assert d
    q = tuple(x // d for x in v)
    return q if next(x for x in q if x) > 0 else tuple(-x for x in q)


def interval(normals):
    low, high = Q(0), None
    for a, b in normals:
        if not a:
            assert b > 0
        elif a > 0:
            low = max(low, Q(-b, a))
        else:
            bound = Q(-b, a)
            high = bound if high is None else min(high, bound)
    assert high is None or low < high
    return low, high


def uniform_planar(m, counts):
    """Free cuts from every cost-order interval, without a rank oracle."""
    E = tuple(hyperbolic_cross(m + 1, 2))
    events = sorted({Q(d - b, a - c)
                     for (a, b), (c, d) in combinations(E, 2)
                     if (d - b) * (a - c) > 0})
    endpoints = [Q(0)] + events + [None]
    states = []
    for low, high in zip(endpoints, endpoints[1:]):
        w = low + 1 if high is None else (low + high) / 2
        B = frozenset(sorted(E, key=lambda a: (a[0] * w + a[1], a))[:m])
        assert B == first_terms((w.numerator, w.denominator), m)
        counts['independent_unrestricted_heap_cuts'] += 1
        if states and states[-1][0] == B:
            states[-1][2] = high
        else:
            states.append([B, low, high])
    for B, low, high in states:
        inequalities = {tuple(x - y for x, y in zip(b, a))
                        for a in B for b in E if b not in B}
        assert interval(inequalities) == (low, high)
        counts['independent_complete_free_cut_intervals'] += 1
    assert all(a[2] == b[1] for a, b in zip(states, states[1:]))
    return states


def budget(cells, counts):
    moments = [profile(B) for B, *_ in cells]
    m = len(cells[0][0])
    rows = []
    weighted = 0
    for i, (a, b) in enumerate(zip(moments, moments[1:])):
        delta = tuple(y - x for x, y in zip(a, b))
        assert delta[0] < 0 < delta[1]
        q = primitive(delta)
        cost = sum(map(abs, q))
        k = gcd(*delta)
        w = cells[i][3]
        assert w is not None and dot((w, 1), delta) == 0
        weighted += k * cost
        rows.append({'primitive_direction': q, 'cost': cost, 'edge_multiple': k})
        counts['independent_planar_profile_jumps'] += 1
    assert weighted == moments[0][0] - moments[-1][0] + moments[-1][1] - moments[0][1]
    total = sum(row['cost'] for row in rows)
    assert total <= weighted <= m * (m - 1)
    counts['independent_planar_telescoping_budgets'] += 1
    return total, rows


def family(m, p, counts):
    # Evaluation, rather than a tagged vector-space model, independently
    # realizes a radical quotient on 3m distinct points.
    inv = lambda j: pow(j, -1, p) if p else Q(1, j)
    points = tuple(point for j in range(1, m + 1)
                   for point in ((j, 0, 0), (j, inv(j), 0),
                                 (j, 0, pow(j, m, p) if p else j**m)))
    assert len(set(points)) == 3 * m
    for x, y, z in points:
        values = [x * y * y - y, z * z - x**m * z, y * z]
        assert all(not (v % p if p else v) for v in values)
        assert any(x == j for j in range(1, m + 1))
        counts['independent_point_generator_annihilations'] += 4

    def point_columns(nodes):
        @lru_cache(None)
        def ev(a):
            out = {}
            for i, point in enumerate(nodes):
                value = 1
                for x, e in zip(point, a):
                    value *= pow(x, e, p) if p else x**e
                    if p:
                        value %= p
                if value:
                    out[i] = value
            return out
        return ev

    torus_ev = point_columns(tuple((j, inv(j)) for j in range(1, m + 1)))
    generic_ev = point_columns(tuple((j, pow(j, m, p) if p else j**m)
                                    for j in range(1, m + 1)))
    torus = complete_planar_fan(m, torus_ev, p, counts)
    generic = complete_planar_fan(m, generic_ev, p, counts)
    hooks = {frozenset([(j, 0) for j in range(k + 1)] +
                       [(0, j) for j in range(1, m - k)]) for k in range(m)}
    assert {row[0] for row in torus} == hooks and len(torus) == m
    free = uniform_planar(m, counts)
    assert [(B, lo, hi) for B, _, lo, hi, _ in generic] == [
        (B, lo, hi) for B, lo, hi in free]
    W1, jumps1 = budget(torus, counts)
    W2, jumps2 = budget(generic, counts)
    assert W1 == sum(m // gcd(j, m) for j in range(1, m))
    assert {tuple(row['primitive_direction']) for row in jumps1} == {
        primitive((j, -(m - j))) for j in range(1, m)}
    if all(m % j for j in range(2, m)):
        assert W1 == m * (m - 1)
    counts['independent_hook_cost_formulas'] += 1

    N = len(points)
    ev = point_columns(points)
    E = tuple(hyperbolic_cross(N, 3))
    core = frozenset((j, 0, 0) for j in range(m))
    cells, certificates = [], []
    for D1, normals1, lo1, hi1, w1 in torus:
        for D2, normals2, lo2, hi2, w2 in generic:
            B = core | {(a, b + 1, 0) for a, b in D1} | {
                (a, 0, b + 1) for a, b in D2}
            assert len(B) == N
            w = (Q(1), 1 / w1, 1 / w2)
            selected, _ = field_basis(sorted(E, key=lambda a: (dot(w, a), a)), ev, p)
            assert frozenset(selected) == B
            counts['independent_full_ground_spatial_greedy_replays'] += 1
            actual, relations = border_normals(tuple(sorted(B)), ev, p)
            assert {tuple(row['head']) for row in relations} == set(minimal_excluded_border(tuple(B)))
            triangle = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
            polygon, predicted = triangle[:], triangle[:]
            for normal in actual:
                assert dot(w, normal) > 0
                polygon = clip(polygon, normal, 1)
            for normal in tuple((a, b, 0) for a, b in normals1) + tuple(
                    (a, 0, b) for a, b in normals2):
                predicted = clip(predicted, normal, 1)
            hull = lambda poly: convex_hull([(v[0], v[1]) for v in poly])
            assert hull(polygon) == hull(predicted)
            counts['independent_complete_spatial_cone_equalities'] += 1
            counts['independent_spatial_border_relations'] += len(relations)
            cells.append((B, polygon))
            certificates.append({'profile': profile(B), 'strict_weight': list(map(str, w)),
                                 'positive_cone': [list(map(str, v)) for v in polygon]})
    assert len(cells) == m * len(generic)
    assert sum(area2(poly) for _, poly in cells) == 1
    edges = positive_edges(cells)
    costs = Counter(wall_cost(row['direction']) for row in edges.values())
    value = sum(cost * multiplicity for cost, multiplicity in costs.items())
    assert value == len(generic) * W1 + m * W2
    assert len(edges) == (m - 1) * len(generic) + m * (len(generic) - 1)
    first = [row for row in edges.values() if row['direction'][2] == 0]
    second = [row for row in edges.values() if row['direction'][1] == 0]
    assert len(first) + len(second) == len(edges)
    assert sum(wall_cost(row['direction']) for row in first) == len(generic) * W1
    assert all(wall_cost(row['direction']) <= N for row in edges.values())
    counts['independent_spatial_undirected_wall_costs'] += len(edges)
    counts['independent_complete_spatial_point_fans'] += 1
    return {'branch_length': m, 'length': N, 'characteristic': p,
            'torus_initial_ideals': len(torus), 'generic_planar_initial_ideals': len(generic),
            'torus_wall_cost': W1, 'generic_planar_wall_cost': W2,
            'initial_ideals': len(cells), 'undirected_walls': len(edges),
            'total_primitive_wall_cost': value,
            'first_branch_wall_pieces': len(first),
            'first_branch_wall_cost': len(generic) * W1,
            'wall_cost_histogram': dict(sorted(costs.items())),
            'torus_jumps': jumps1, 'generic_planar_jumps': jumps2,
            'complete_spatial_cones': certificates}


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    counts = defaultdict(int)
    fixtures = [family(m, p, counts) for m, p in (
        (2, 0), (3, 0), (4, 0), (5, 0), (7, 1000000007),
        (8, 1000000007), (11, 1000000007), (13, 1000000007))]
    paths = ['docs/wall_cost_products.md', 'scripts/run_wall_cost_product_checks.py',
             'scripts/run_wall_coexistence_checks.py', 'scripts/run_signed_face_checks.py',
             'scripts/run_structural_complexity_checks.py', 'scripts/run_free_wall_checks.py',
             'scripts/run_planar_dyadic_bound_checks.py', 'scripts/run_interleaving_research.py']
    report = {'schema_version': 1, 'status': 'passed',
              'scope': 'Exact complete finite fans support the uniform planar telescoping bound, existing fan-product cost identity and explicit radical rational family with V=Theta(N^3 log N) along prime branch lengths. Finite-field examples are individually verified specializations. No asymptotic law is fitted, unrestricted upper bound improved, polynomial fan separation proved, or novelty established.',
              'command': 'python3 -B scripts/run_wall_cost_product_checks.py',
              'counts': dict(sorted(counts.items())), 'fixtures': fixtures,
              'source_sha256': {s: hashlib.sha256((ROOT / s).read_bytes()).hexdigest()
                                for s in paths}}
    target = ROOT / 'results/research/wall_cost_product_checks.json'
    target.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'counts': report['counts'],
                      'fixtures': [{k: row[k] for k in ('branch_length', 'length',
                                   'characteristic', 'initial_ideals',
                                   'total_primitive_wall_cost')} for row in fixtures]}, indent=2))


if __name__ == '__main__':
    main()
