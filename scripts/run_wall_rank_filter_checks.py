#!/usr/bin/env python3
"""Exact rank and scope controls; no fitted asymptotic claim."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd, prod
from pathlib import Path
from random import Random
import hashlib
import json

from run_extremal_initial_ideals_checks import cone_for_basis, rectangle_model
from run_free_wall_checks import first_terms
from run_oracle_enumeration_checks import curve_model, point_model, rectangle_fixture
from run_primitive_wall_checks import divisors, threshold
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross
from run_wall_coexistence_checks import core_data, positive_edges

ROOT = Path(__file__).resolve().parents[1]


def rank(columns, p=0):
    """Dense row elimination, independent of the sparse quotient oracle."""
    keys = sorted({k for v in columns for k in v})
    rows = [[Q(v.get(k, 0)) if not p else v.get(k, 0) % p for k in keys]
            for v in columns]
    r = 0
    for j in range(len(keys)):
        pivot = next((i for i in range(r, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        inv = pow(int(rows[r][j]), -1, p) if p else 1 / rows[r][j]
        rows[r] = [(x * inv) % p if p else x * inv for x in rows[r]]
        for i in range(r + 1, len(rows)):
            c = rows[i][j]
            if c:
                rows[i] = [(x-c*y) % p if p else x-c*y
                           for x, y in zip(rows[i], rows[r])]
        r += 1
        if r == len(rows):
            break
    return r


def cheap_hull(q):
    u, v = tuple(max(x, 0) for x in q), tuple(max(-x, 0) for x in q)
    pairs = [(i, j) for i in range(3) for j in range(3) if q[i] > 0 > q[j]]
    return frozenset(a for a in product(*(range(abs(x)+1) for x in q))
                     if a not in (u, v) and all(
                         a[i]*(-q[j])+a[j]*q[i] <= q[i]*(-q[j]) for i, j in pairs))


def models():
    saved = json.loads((ROOT/'results/research/oracle_enumeration_checks.json').read_text())
    for case in saved['direct_positive_cone_traversal']['fixtures']:
        label, n, p = case['model'], case['length'], case['characteristic']
        if label in ('one_point', 'double_point', 'curvilinear_three'):
            model = curve_model(n)
        elif label == 'nonlinear_rectangle':
            model = rectangle_fixture({(1, 0): 1, (0, 1): 1, (1, 1): 1})
        elif label == 'monomial_box':
            model = rectangle_fixture({})
        else:
            assert label == 'three_distinct_points'
            model = point_model()
        cells = [(frozenset(map(tuple, row['standard_set'])),
                  [tuple(Q(x) for x in v) for v in row['polygon']])
                 for row in case['complete_cells']]
        yield label, n, p, model[1], cells
    saved = json.loads((ROOT/'results/research/extremal_initial_ideals_checks.json').read_text())
    for case in saved['rectangular_complete_intersections']:
        r, s = case['rectangle']
        n, p = case['quotient_length'], case['characteristic']
        F = {(a, b): c for a, b, c in case['polynomial_terms']}
        _, _, residue = rectangle_model(r, s, F, p)
        ev = lambda a: {i: c for i, c in enumerate(residue(a)) if c}
        cells = []
        for row in case['witnesses']:
            B = frozenset(map(tuple, row['basis']))
            polygon, _ = cone_for_basis(B, residue, p)
            cells.append((B, polygon))
        yield 'rectangle_%s_%s' % (r, s), n, p, ev, cells


def quotient_controls(counts):
    fixtures = []
    for label, n, p, ev, cells in models():
        rows = []
        for edge in positive_edges(cells).values():
            q, w = edge['direction'], edge['weight']
            u, v, C, _ = core_data(q)
            L = cheap_hull(q)
            H = divisors(tuple(abs(x)//2 for x in q))
            assert C <= L and H <= L
            assert all(dot(w, a) < dot(w, u) for a in L)
            selected, pivots = field_basis(tuple(sorted(C)), ev, p)
            assert len(selected) == len(C)
            rL = rank([ev(a) for a in L], p)
            residual = [field_reduce(ev(a), pivots, p)[0] for a in L-C]
            assert rank(residual, p) == rL-len(C) <= n-threshold(q)
            assert rL < n
            counts['independent_actual_wall_cheap_hull_ranks'] += 1
            counts['independent_strictly_cheaper_monomials'] += len(L)
            counts['independent_core_contracted_rank_equalities'] += 1
            counts['cheap_hull_not_literal_endpoint_controls'] += int(
                any(not L <= B for B in edge['endpoints']))
            rows.append({'direction': q, 'wall_weight': list(map(str, w)),
                         'cheap_hull_size': len(L), 'cheap_hull_rank': rL,
                         'core_size': len(C), 'half_box_size': len(H),
                         'cheap_hull_nullity': len(L)-rL})
        box_rows = []
        if n <= 5:
            caps = tuple(product(range(n), repeat=3))
            ranks = {a: rank([ev(b) for b in divisors(a)], p) for a in caps}
            for a in caps:
                if ranks[a] != n or any(a[i] and ranks[tuple(
                        x-(i == j) for j, x in enumerate(a))] == n for i in range(3)):
                    continue
                assert sum(a) <= n-1
                box_rows.append(a)
            assert box_rows and len(box_rows) <= comb(n+1, 2)
            counts['independent_complete_spanning_box_diagrams'] += 1
            counts['independent_minimal_spanning_boxes'] += len(box_rows)
            if all(ranks[a] == min(n, prod(x+1 for x in a)) for a in caps):
                oversampling(n, ev, p, counts)
        fixtures.append({'model': label, 'length': n, 'characteristic': p,
                         'initial_cells': len(cells), 'walls': rows,
                         'minimal_spanning_boxes': box_rows})
    return fixtures


def oversampling(n, ev, p, counts):
    M = max(1, 64*(n-1))
    rng = Random(n+41)
    witnesses = []
    for _ in range(8):
        while True:
            weights = tuple(rng.randrange(M**4+1, 3*M**4+4) for _ in range(3))
            free = sorted(first_terms(weights, M+1), key=lambda a: (dot(weights, a), a))
            if len({dot(weights, a) for a in free}) == M+1:
                break
        direct, _ = field_basis(tuple(free[:M]), ev, p)
        global_basis, _ = field_basis(tuple(sorted(hyperbolic_cross(n, 3),
                                   key=lambda a: (dot(weights, a), a))), ev, p)
        assert len(direct) == n and set(direct) == set(global_basis)
        stop = max(free.index(a)+1 for a in direct)
        assert stop <= M
        if n > 1:
            tau = max(dot(weights, a) for a in direct)
            caps = tuple(tau//(4*w) for w in weights)
            box = divisors(caps)
            assert all(dot(weights, a) < tau for a in box)
            assert rank([ev(a) for a in box], p) < n and len(box) < n
            cheap = sum(dot(weights, a) <= tau for a in free)
            assert cheap == stop <= 64*len(box)
        counts['independent_uniform_free_prefix_filter_replays'] += 1
        witnesses.append({'weight': weights, 'terms_scanned': stop, 'uniform_prefix_length': M})
    return witnesses


def nonuniform_saturated_points(counts):
    n, p = 8, 1000000007
    parameters = ((2202, 9326), (1034, 4180), (1932, 8118), (7365, 7738),
                  (6220, 3440), (1538, 7994), (465, 6387), (7091, 9953))
    points = []
    for s, t in parameters:
        inv = pow(1+s*s+t*t, -1, p)
        points.append((2*s*inv % p, 2*t*inv % p, (1-s*s-t*t)*inv % p))
    assert len(set(points)) == n
    assert all((x*x+y*y+z*z-1) % p == 0 for x, y, z in points)
    @lru_cache(None)
    def ev(a):
        return {j: prod(pow(x, e, p) for x, e in zip(point, a)) % p
                for j, point in enumerate(points)}
    for a in product(range(n), repeat=3):
        assert rank([ev(b) for b in divisors(a)], p) == min(n, prod(x+1 for x in a))
        counts['independent_nonuniform_point_saturation_boxes'] += 1
    relation_support = ((0, 0, 0), (2, 0, 0), (0, 2, 0), (0, 0, 2))
    assert rank([ev(a) for a in relation_support], p) == 3
    witnesses = oversampling(n, ev, p, counts)
    return {'length': n, 'characteristic': p, 'points': points,
            'nonuniform_relation': 'x^2+y^2+z^2-1',
            'exhaustive_saturation_boxes': n**3, 'filter_witnesses': witnesses}


def diagonal_controls(counts):
    rows = []
    for n in range(2, 17):
        base = 3*n+1
        w = base*base, base*base+base, base*base+1
        tau = (n-1)*base*base
        candidates = tuple(product(range(n), repeat=3))
        before = sorted((a for a in candidates if dot(w, a) <= tau),
                        key=lambda a: dot(w, a))
        assert len({dot(w, a) for a in before}) == len(before)
        assert len(before) == comb(n+1, 3)+1
        s = max(prod(x+1 for x in a) for a in candidates if sum(a) <= n-2)
        h, rem = divmod(n-2, 3)
        assert s == (h+2)**rem*(h+1)**(3-rem)
        assert 1+Q(s, 27) <= len(before) <= 1+64*s
        for p in (0, 2, 3, 5):
            ev = lambda a: {sum(a): 1} if sum(a) < n else {}
            selected, _ = field_basis(tuple(before), ev, p)
            assert set(selected) == {(i, 0, 0) for i in range(n)}
            for axis in range(3):
                weights = tuple(1 if j == axis else n+1+j for j in range(3))
                B, _ = field_basis(tuple(sorted(candidates, key=lambda a: (dot(weights, a), a))), ev, p)
                assert set(B) == {tuple(i if j == axis else 0 for j in range(3)) for i in range(n)}
            counts['independent_diagonal_scanning_and_initial_controls'] += 1
        rows.append({'length': n, 'nonsaturating_box_volume': s,
                     'near_equal_weight': w, 'ordinary_terms_through_completion': len(before),
                     'initial_ideals': 3})
    return rows


def contrast_controls(counts):
    fixtures = []
    for q in ((4, -2, -3), (7, -2, -3), (9, -3, -4)):
        u, v, C, _ = core_data(q)
        L, n = cheap_hull(q), threshold(q)
        assert len(L) > len(C)
        D = divisors(u) | divisors(v)
        labels = {a: i for i, a in enumerate(sorted(D-{v}))}
        socle = lambda a: {labels[u if a == v else a]: 1} if a in D else {}
        assert not L <= C | {u} and not L <= C | {v}
        counts['cheap_hull_not_literal_endpoint_controls'] += 2
        for p in (0, 2, 3, 5):
            assert rank([socle(a) for a in C], p) == n-1
            assert rank([socle(a) for a in L], p) == n-1
            counts['independent_sharp_socle_rank_controls'] += 1
        A, B, p = q[0]+1, (q[0]+1)*(abs(q[1])+1), 1000000007
        ev = lambda a: {j: pow(j+1, a[0]+A*a[1]+B*a[2], p) for j in range(n)}
        assert rank([ev(a) for a in C], p) == n-1
        assert rank([ev(a) for a in L], p) == n
        counts['independent_core_only_false_positive_rejections'] += 1
        cardinality = sum((abs(q[1])*(q[0]-i)//q[0]+1)*
                          (abs(q[2])*(q[0]-i)//q[0]+1)
                          for i in range(q[0]+1))-2
        assert len(L) == cardinality
        fixtures.append({'direction': q, 'minimum_length': n,
                         'cheap_hull_size': len(L), 'socle_rank': n-1,
                         'specialized_point_rank': n, 'point_characteristic': p,
                         'point_exponent_codes': (1, A, B)})
    return fixtures


def hypersurface_controls(counts):
    fixtures = []
    for D, p in ((2, 0), (2, 101), (3, 101)):
        @lru_cache(None)
        def remainder(a):
            i, j, k = a
            if k < D:
                return {a: 1}
            out = defaultdict(int)
            for alpha, sign in (((i, j, k-D), 1), ((i+D, j, k-D), -1),
                                ((i, j+D, k-D), -1)):
                for b, c in remainder(alpha).items():
                    out[b] += sign*c
            return {b: (c % p if p else c) for b, c in out.items()
                    if (c % p if p else c)}
        total = 0
        for a in product(range(D+1), repeat=3):
            expected = prod(x+1 for x in a)-prod(max(x-D+1, 0) for x in a)
            actual = rank([remainder(b) for b in divisors(a)], p)
            assert actual == expected
            counts['independent_hypersurface_box_dimensions'] += 1
            total += 1
        fixtures.append({'coordinate_degree': D, 'characteristic': p,
                         'boxes': total, 'point_length_threshold': (D+1)**3-1})
    return fixtures


def matching(S, neighborhoods):
    occupied = {}
    def add(e, seen):
        for j, row in enumerate(neighborhoods):
            if e in row and j not in seen:
                seen.add(j)
                if j not in occupied or add(occupied[j], seen):
                    occupied[j] = e
                    return True
        return False
    return all(add(e, set()) for e in S)


def factor(neighborhoods, weights, counts):
    m, r = len(weights), len(neighborhoods)
    subsets = [S for k in range(r+1) for S in combinations(range(m), k)]
    for seed in range(1, 100):
        rng = Random(seed)
        columns = [{j: rng.randrange(1, 10000) for j, row in enumerate(neighborhoods)
                    if e in row} for e in range(m)]
        if all(rank([columns[e] for e in S]) == len(S)
               for S in subsets if matching(S, neighborhoods)):
            break
    else:
        raise AssertionError('No verified rational transversal specialization')
    for S in subsets:
        assert (rank([columns[e] for e in S]) == len(S)) == matching(S, neighborhoods)
        counts['independent_matching_representation_subsets'] += 1
    bases = [S for S in combinations(range(m), r) if matching(S, neighborhoods)]
    profiles = {tuple(sum(weights[e][j] for e in S) for j in range(2)): S for S in bases}
    assert len(profiles) == len(bases)
    active = []
    for a in profiles:
        lo, hi = Q(0), None
        for b in profiles:
            intercept, slope = b[0]-a[0], b[1]-a[1]
            if slope > 0:
                lo = max(lo, Q(-intercept, slope))
            elif slope < 0:
                bound = Q(-intercept, slope)
                hi = bound if hi is None else min(hi, bound)
            else:
                assert a == b
        if hi is None or lo < hi:
            witness = lo+1 if hi is None else (lo+hi)/2
            assert all(dot((1, witness), a) < dot((1, witness), b)
                       for b in profiles if b != a)
            active.append(a)
    assert len(active) >= 3
    return profiles, active, columns


def product_controls(counts):
    fixtures = []
    for neighborhoods1, neighborhoods2 in (
            (((0, 1, 3), (1, 2, 4), (0, 2, 5)),
             ((0, 1, 4), (1, 2, 5), (0, 3, 5))),
            (((0, 1, 2), (3, 4, 5)),
             ((0, 1, 3), (1, 2, 4), (0, 2, 5)))):
        weights1 = tuple((2**(5-i), 2**i) for i in range(6))
        weights2 = tuple((3**(5-i), 3**i) for i in range(6))
        P, A, columns1 = factor(neighborhoods1, weights1, counts)
        QP, B, columns2 = factor(neighborhoods2, weights2, counts)
        spatial = {(a[0]+b[0], a[1], b[1]) for a in P for b in QP}
        predicted = {(a[0]+b[0], a[1], b[1]) for a in A for b in B}
        total, found = Q(0), set()
        for a in spatial:
            poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
            for b in spatial:
                if b != a:
                    poly = clip(poly, tuple(y-x for x, y in zip(a, b)), 1)
                    if len(poly) < 3:
                        break
            if len(poly) >= 3 and area2(poly):
                found.add(a)
                total += area2(poly)
                counts['independent_complete_direct_sum_spatial_cones'] += 1
        assert found == predicted and len(found) == len(A)*len(B) and total == 1
        fixtures.append({'first_rank': len(neighborhoods1), 'second_rank': len(neighborhoods2),
                         'first_neighborhoods': neighborhoods1, 'second_neighborhoods': neighborhoods2,
                         'first_labels': weights1, 'second_labels': weights2,
                         'first_columns': columns1, 'second_columns': columns2,
                         'factor_positive_profiles': (len(A), len(B)),
                         'direct_sum_positive_profiles': len(found), 'total_simplex_area2': str(total)})
    return fixtures


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    counts = defaultdict(int)
    fans = quotient_controls(counts)
    nonuniform = nonuniform_saturated_points(counts)
    diagonal = diagonal_controls(counts)
    contrasts = contrast_controls(counts)
    hypersurfaces = hypersurface_controls(counts)
    products = product_controls(counts)
    paths = ('docs/wall_rank_filters.md', 'scripts/run_wall_rank_filter_checks.py',
             'paper/corner_cut_growth/main.tex', 'docs/oracle_fan_enumeration.md',
             'scripts/run_free_wall_checks.py',
             'scripts/run_wall_coexistence_checks.py', 'scripts/run_primitive_wall_checks.py',
             'scripts/run_extremal_initial_ideals_checks.py', 'scripts/run_oracle_enumeration_checks.py',
             'scripts/run_signed_face_checks.py', 'scripts/run_structural_complexity_checks.py',
             'results/research/oracle_enumeration_checks.json',
             'results/research/extremal_initial_ideals_checks.json')
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_wall_rank_filter_checks.py',
              'scope': 'Exact saved actual fans verify necessary cheap-hull and contracted ranks. Socle and point columns contrast core independence with rank saturation. Exhaustive box ranks of a nonuniform quadric point ideal and independent free-prefix filtering support the uniform oversampling proof. Diagonal algebras check sharp worst-case scanning despite three initials. Polynomial division checks hypersurface box dimensions. Finite represented transversal direct sums verify the product identity, not the published asymptotic lower construction. No unrestricted improvement, separating ideal family or novelty established.',
              'uniform_results': {
                  'scanning_parameter': 'For fixed d, the largest nonsaturating anchored box volume s and worst free-prefix completion rank M obey 1+s/d^d<=M<=1+(d+1)^d*s, allowing infinity.',
                  'finiteness': 'Uniform free-prefix completion is finite exactly when every coordinate separately is cyclic; the universal finite-case bound n^d is sharp in power.',
                  'whole_class_count': 'Every quotient with s=O_d(n) has g=O_d(n^(d-1)*log(n+1)^(d-1)), by constant oversampling and the existing ordered-prefix theorem. This has the generic sharp order, using the existing generic lower theorem.',
                  'nonuniform_infinite_class': 'Sufficiently generic n distinct points on an absolutely irreducible hypersurface of positive coordinate degrees d_i have rectangular rank saturation when n<=prod(d_i+1)-1. The characteristic-zero Fermat family has a four-term dependence.',
                  'enumeration_consequence': 'For d=3 this class has the existing oracle-enumeration bound O(n^3*L^2+n^(omega+2)*L^3), excluding coefficient bit costs.',
                  'cheap_hull_obstruction': 'Every actual q-wall has rank(L_q)<=n-1, and its image modulo the compulsory core has rank at most n-nu(q).',
                  'general_matroid_scope': 'A direct sum of two published Eppstein transversal lower families gives Omega(r^(8/3)) positive profiles on Theta(r) elements and rank 2r. Therefore the ordinary spatial k-set bound does not extend to arbitrary represented matroids.',
                  'unrestricted_G3_improved': False, 'polynomial_lower_separation_proved': False,
                  'novelty_established': False, 'goal_complete': False},
              'counts': dict(sorted(counts.items())), 'actual_quotient_fans': fans,
              'core_vs_cheap_hull_contrasts': contrasts, 'hypersurface_dimension_checks': hypersurfaces,
              'nonuniform_rectangular_saturated_point_ideal': nonuniform,
              'diagonal_uniform_scanning_obstruction': diagonal,
              'represented_transversal_products': products,
              'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in paths}}
    (ROOT/'results/research/wall_rank_filter_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'counts': report['counts'],
                      'actual_quotient_fans': len(fans)}, indent=2))


if __name__ == '__main__':
    main()
