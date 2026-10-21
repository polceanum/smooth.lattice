#!/usr/bin/env python3
"""Exact bounded verification of truncated-state oracle enumeration.

The hull routine below intentionally enumerates triples. The complexity
theorem uses the cited three-dimensional hull algorithm, not this checker.
"""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import gcd
from pathlib import Path
import hashlib
import json

from run_higher_interaction_checks import rank_trace
from run_interleaving_research import minimal_excluded_border
from run_liaison_fan_checks import rectangle_graph
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
    R = n*n+1
    candidates = set(profiles)
    for b in profiles:
        assert sum(b) < R
        candidates.update(tuple(b[i]+(R-sum(b))*(i == j) for i in range(3)) for j in range(3))
    audit_points, trace = set(), []

    def oracle(w, record=False):
        B, _ = rank_trace(sorted(E, key=lambda a: monomial_key(w, a)), ev, p, set())
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
    assert [f['initial_ideals'] for f in fixtures[:9]] == [1, 3, 3, 6, 5, 6, 6, 4, 1]
    sources = ('docs/oracle_fan_enumeration.md', 'scripts/run_oracle_enumeration_checks.py',
               'docs/rational_shadow_bound.md', 'docs/power_carrier_bounds.md',
               'scripts/run_higher_interaction_checks.py', 'scripts/run_hypersurface_truncation_checks.py',
               'scripts/run_signed_face_checks.py', 'scripts/run_liaison_fan_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_extremal_initial_ideals_checks.py', 'scripts/run_socle_stability_checks.py',
               'scripts/run_higher_seed_checks.py')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_oracle_enumeration_checks.py',
              'counts': dict(counts), 'fixtures': fixtures,
              'uniform_results': {'truncated_state_vertices': 'at most 4*g',
                                  'vertex_queries': 'at most 12*g+1, including initialization',
                                  'field_operations': 'O((g+1)*n^3*(1+log(n+1))^2)',
                                  'geometric_operations': 'O((g+1)^2*log(g+2)+(g+1)*n*(1+log(n+1))^3)',
                                  'unrestricted_shadow_consequence': 'O(n^(47/8)*(1+log(n+1))^(19/8)) field operations',
                                  'classical_count_composition': 'O(n^6*(1+log(n+1))^2)',
                                  'given_input': 'One marked reduced Groebner basis, or commuting coordinate matrices with a cyclic vector.',
                                  'union_output': 'All reduced bases; their union may retain repeated polynomials.'},
              'new_unrestricted_G3_count': False, 'novelty_established': False,
              'bit_complexity_bound_established': False, 'practical_speedup_established': False,
              'goal_complete': False,
              'scope': 'An exact arithmetic consequence of the internal counting theorem using attributed reconstruction and linear algebra. The checker uses a deliberately slow small hull routine. It does not establish current-best priority, asymptotic sharpness, coefficient bit costs or production performance.',
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/'results/research/oracle_enumeration_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
