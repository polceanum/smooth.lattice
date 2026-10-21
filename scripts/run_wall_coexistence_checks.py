#!/usr/bin/env python3
"""Exact controls for a fixed-core contraction of parallel initial walls."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd, lcm, prod
from pathlib import Path
import hashlib
import json

from run_extremal_initial_ideals_checks import cone_for_basis, rectangle_model
from run_higher_interaction_checks import rank_trace
from run_interleaving_research import minimal_excluded_border
from run_oracle_enumeration_checks import curve_model, point_model, rectangle_fixture
from run_planar_dyadic_bound_checks import convex_hull
from run_primitive_wall_checks import canonical, divisors, threshold
from run_signed_face_checks import field_basis, field_reduce, normals as border_normals
from run_structural_complexity_checks import area2, clip, cross, dot, finite_staircases, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]


def core_data(q):
    u = tuple(max(x, 0) for x in q)
    v = tuple(max(-x, 0) for x in q)
    core = (divisors(u) | divisors(v)) - {u, v}
    tips = tuple(sorted(minimal_excluded_border(tuple(sorted(core)))))
    assert len(core) == threshold(q) - 1 and len(tips) <= 6
    return u, v, frozenset(core), tips


def outside_divisors(a, u, v):
    inside = (prod(min(x, y) + 1 for x, y in zip(a, u)) +
              prod(min(x, y) + 1 for x, y in zip(a, v)) - 1)
    inside -= all(x <= y for x, y in zip(u, a))
    inside -= all(x <= y for x, y in zip(v, a))
    return prod(x + 1 for x in a) - inside


def relative_ground(q, r):
    u, v, core, tips = core_data(q)
    V = tuple(hyperbolic_cross(r, 3))
    union = {tuple(a + b for a, b in zip(h, c)) for h in tips for c in V}
    ground = tuple(sorted(a for a in union if outside_divisors(a, u, v) <= r))
    assert all(a not in core for a in ground)
    return ground, tips, len(union)


def grading_histogram(B, q):
    pivot = next(i for i, x in enumerate(q) if x)
    result = defaultdict(int)
    for a in B:
        degree = tuple(q[pivot] * a[i] - q[i] * a[pivot]
                       for i in range(3) if i != pivot)
        result[degree] += 1
    return dict(result)


def downset_controls(counts):
    cached = {}
    fixtures = []
    for n in range(1, 9):
        states = finite_staircases(n, 3)
        moments = {B: tuple(sum(a[i] for a in B) for i in range(3)) for B in states}
        accepted = 0
        for B, D in combinations(states, 2):
            delta = tuple(x - y for x, y in zip(moments[B], moments[D]))
            if not any(delta):
                continue
            q = canonical(delta)
            if not min(q) < 0 < max(q):
                continue
            if grading_histogram(B, q) != grading_histogram(D, q):
                continue
            accepted += 1
            u, v, core, tips = core_data(q)
            r = n - len(core)
            assert r >= 1 and core <= B & D
            key = (q, r)
            if key not in cached:
                ground, _, union_size = relative_ground(q, r)
                cached[key] = set(ground)
                for a in ground:
                    assert outside_divisors(a, u, v) == len(divisors(a) - core)
                    counts['independent_relative_divisor_cardinalities'] += 1
                assert len(ground) <= len(tips) * len(tuple(hyperbolic_cross(r, 3)))
                counts['independent_core_relative_ground_covers'] += 1
            assert (B | D) - core <= cached[key]
            counts['independent_compulsory_core_pairs'] += 1
            if r == 1:
                assert B | D == core | {u, v}
                counts['independent_saturated_endpoint_pairs'] += 1
        fixtures.append({'length': n, 'all_downsets': len(states), 'compatible_pairs': accepted})
    # Separate exhaustive boxes verify that no points satisfying the relative
    # criterion were omitted by generating only the translated cover.
    for q in ((1, -1, -1), (2, -1, -2), (3, -2, 0), (3, -2, -1)):
        u, v, core, tips = core_data(q)
        for r in (1, 2, 3, 5):
            limit = tuple(max(a[i] for a in core) + r for i in range(3))
            direct = {a for a in product(*(range(x + 1) for x in limit))
                      if a not in core and len(divisors(a) - core) <= r}
            generated, _, _ = relative_ground(q, r)
            assert direct == set(generated)
            counts['independent_complete_relative_ground_boxes'] += 1
    return fixtures


def projected_label(a, q):
    pivot = next(i for i, x in enumerate(q) if x)
    return tuple(q[pivot] * a[i] - q[i] * a[pivot]
                 for i in range(3) if i != pivot)


def profile(B):
    return tuple(sum(a[i] for a in B) for i in range(3))


def planar_support_polygon(ground, columns, p, q, r, counts):
    labels = {a: projected_label(a, q) for a in ground}
    cache = {}

    def support(cost):
        cost = tuple(cost)
        if cost not in cache:
            chosen, _ = rank_trace(sorted(ground, key=lambda a: (dot(cost, labels[a]), labels[a], a)),
                                   columns.__getitem__, p, set())
            assert len(chosen) == r
            value = tuple(sum(labels[a][i] for a in chosen) for i in range(2))
            cache[cost] = value
            counts['residual_planar_greedy_support_queries'] += 1
        return cache[cost]

    initial = {support(c) for c in ((1, 0), (-1, 0), (0, 1), (0, -1))}
    points = set(initial)
    if len(points) > 1:
        a, b = min(points), max(points)
        queue = [(a, b), (b, a)]
        processed = set()
        while queue:
            a, b = queue.pop()
            if (a, b) in processed:
                continue
            processed.add((a, b))
            inward = (a[1] - b[1], b[0] - a[0])
            c = support(inward)
            if dot(inward, c) < dot(inward, a):
                points.add(c)
                queue.extend(((a, c), (c, b)))
    polygon = tuple(convex_hull(points))
    if len(polygon) >= 3:
        for a, b in zip(polygon, polygon[1:] + polygon[:1]):
            inward = (a[1] - b[1], b[0] - a[0])
            assert dot(inward, support(inward)) == dot(inward, a)
    active = tuple(a for a in ground if columns[a])
    if len(active) >= r and comb(len(active), r) <= 3000:
        all_profiles = set()
        for B in combinations(active, r):
            selected, _ = rank_trace(B, columns.__getitem__, p, set())
            counts['independent_residual_basis_eliminations'] += 1
            if len(selected) == r:
                all_profiles.add(tuple(sum(labels[a][i] for a in B) for i in range(2)))
        assert set(polygon) == set(convex_hull(all_profiles))
        counts['independent_exhaustive_residual_shadow_equalities'] += 1
    return polygon, support


def positive_edges(cells):
    moments = [profile(B) for B, _ in cells]
    assert len(set(moments)) == len(moments)
    result = {}
    for i, (B, polygon) in enumerate(cells):
        for a, b in zip(polygon, polygon[1:] + polygon[:1]):
            if any(a[j] == b[j] == 0 for j in range(3)):
                continue
            midpoint = tuple((x + y) / 2 for x, y in zip(a, b))
            assert min(midpoint) > 0
            values = [dot(midpoint, moment) for moment in moments]
            best = min(values)
            winners = [j for j, value in enumerate(values) if value == best]
            assert i in winners and len(winners) == 2
            j = next(j for j in winners if j != i)
            edge = tuple(sorted((i, j)))
            normal = cross(a, b)
            scale = lcm(*(x.denominator for x in normal))
            q = canonical(tuple(int(x * scale) for x in normal))
            assert canonical(tuple(x - y for x, y in zip(moments[i], moments[j]))) == q
            if edge in result:
                assert result[edge]['direction'] == q
            else:
                result[edge] = {'direction': q, 'weight': midpoint, 'endpoints': (B, cells[j][0])}
    return result


def check_fan(label, n, ev, p, cells, counts):
    edges = positive_edges(cells)
    groups = defaultdict(list)
    for edge in edges.values():
        groups[edge['direction']].append(edge)
    fixtures = []
    for q, rows in sorted(groups.items()):
        u, v, core, tips = core_data(q)
        r = n - len(core)
        assert r >= 1 and all(core <= B for row in rows for B in row['endpoints'])
        selected, pivots = field_basis(tuple(sorted(core)), ev, p)
        assert len(selected) == len(core)
        ground, _, union_size = relative_ground(q, r)
        columns = {a: field_reduce(ev(a), pivots, p)[0] for a in ground}
        rank_basis, _ = rank_trace(ground, columns.__getitem__, p, set())
        assert len(rank_basis) == r
        polygon, support = planar_support_polygon(ground, columns, p, q, r, counts)
        collapsed = set()
        for row in rows:
            for B in row['endpoints']:
                residual = B - core
                assert residual <= set(ground) and len(residual) == r
                independent, _ = rank_trace(sorted(residual), columns.__getitem__, p, set())
                assert len(independent) == r
                point = projected_label(profile(residual), q)
                assert point in polygon
                pivot = next(i for i, x in enumerate(q) if x)
                cost = tuple(row['weight'][i] / q[pivot] for i in range(3) if i != pivot)
                assert dot(cost, point) == dot(cost, support(cost))
                counts['independent_original_edge_minima_after_contraction'] += 1
            assert len({projected_label(profile(B - core), q) for B in row['endpoints']}) == 1
            collapsed.add(point)
        assert len(collapsed) == len(rows) <= len(polygon)
        if r == 1:
            assert len(rows) == 1
        outside = next(((B, poly) for B, poly in cells if not core <= B), None)
        guard = None
        if outside:
            B, poly = outside
            weight = tuple(sum(v[i] for v in poly)/len(poly) for i in range(3))
            chosen, _ = rank_trace(sorted(ground, key=lambda a: (dot(weight, a), a)),
                                   columns.__getitem__, p, set())
            assert len(chosen) == r and min(weight) > 0
            forced_profile = tuple(x+y for x, y in zip(profile(core), profile(chosen)))
            global_value, forced_value = dot(weight, profile(B)), dot(weight, forced_profile)
            assert global_value < forced_value
            counts['independent_outside_wall_contraction_cost_gaps'] += 1
            guard = {'positive_weight': list(map(str, weight)),
                     'global_minimum': str(global_value),
                     'forced_core_minimum': str(forced_value)}
        counts['independent_same_direction_contracted_shadows'] += 1
        counts['independent_original_parallel_edges'] += len(rows)
        fixtures.append({'direction': q, 'minimum_colength': threshold(q), 'residual_rank': r,
                         'parallel_walls': len(rows), 'compulsory_core': sorted(core),
                         'translated_anchors': tips, 'translated_union_size': union_size,
                         'relative_ground_size': len(ground), 'residual_shadow_vertices': polygon,
                         'collapsed_original_edges': sorted(collapsed),
                         'outside_wall_scope_guard': guard})
    return {'model': label, 'length': n, 'characteristic': p, 'initial_cells': len(cells),
            'all_positive_walls': len(edges), 'directions': fixtures}


def saved_fan_controls(counts):
    saved = json.loads((ROOT / 'results/research/oracle_enumeration_checks.json').read_text())
    assert saved['status'] == 'passed'
    fixtures = []
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
        _, ev, _, _ = model
        cells = [(frozenset(map(tuple, c['standard_set'])),
                  [tuple(Q(x) for x in v) for v in c['polygon']])
                 for c in case['complete_cells']]
        fixtures.append(check_fan(label, n, ev, p, cells, counts))
    saved = json.loads((ROOT / 'results/research/extremal_initial_ideals_checks.json').read_text())
    assert saved['status'] == 'passed'
    for case in saved['rectangular_complete_intersections']:
        r, s = case['rectangle']
        n, p = case['quotient_length'], case['characteristic']
        F = {(a, b): c for a, b, c in case['polynomial_terms']}
        _, _, residue = rectangle_model(r, s, F, p)

        def ev(a):
            return {i: c for i, c in enumerate(residue(a)) if c}

        cells = []
        for witness in case['witnesses']:
            B = frozenset(map(tuple, witness['basis']))
            polygon, _ = cone_for_basis(B, residue, p)
            cells.append((B, polygon))
        fixtures.append(check_fan('nonlinear_rectangle_%s_%s' % (r, s), n, ev, p, cells, counts))
    return fixtures


def repeated_wall_family(counts):
    fixtures = []
    for m in range(1, 9):
        a, b = m**4 + 1, m**4 + 2
        u, v = (a, 0, 0), (0, b, 0)
        corners = tuple((i, 0, 2*m-i+int(i < m)) for i in range(2*m))
        D = divisors(u) | divisors(v)
        for corner in corners:
            D.update(divisors(corner))
        pairs = [(u, v)] + [(corners[i], corners[2*m-1-i]) for i in range(m)]
        canonical_column = {x: min(left, right) for left, right in pairs for x in (left, right)}
        targets = {canonical_column.get(x, x) for x in D}
        n = a + b + 2*m*m + m
        assert len(D) == a + b + 1 + 2*m*m + 2*m
        assert len(targets) == n
        q = (a, -b, 0)
        assert gcd(a, b) == 1
        _, _, core, _ = core_data(q)
        residual = n - len(core)
        assert residual == 2*m*m + m + 1
        ground, tips, union_size = relative_ground(q, residual)
        for x in D:
            for i in range(3):
                y = tuple(c + (i == j) for j, c in enumerate(x))
                if y in D:
                    assert all(x not in pair for pair in pairs)
                counts['independent_socle_annihilation_edges'] += 1
        # Different representatives of a class are maximal, so multiplication
        # by each variable sends both to zero. Check all commuting squares.
        def image(x, i):
            y = tuple(c + (i == j) for j, c in enumerate(x))
            return canonical_column.get(y, y) if y in D else None

        for x in targets:
            for i, j in combinations(range(3), 2):
                left = image(x, i)
                right = image(x, j)
                assert (None if left is None else image(left, j)) == (None if right is None else image(right, i))
                counts['independent_actual_family_commuting_squares'] += 1
        ratios = sorted(Q(d, d + 1) for d in range(1, 2*m, 2))
        samples = [ratios[0] / 2] + [(x+y)/2 for x, y in zip(ratios, ratios[1:])] + [(ratios[-1]+1)/2]
        states, parallel_edges = set(), set()
        rows = []
        for theta in samples:
            tied = (Q(b), Q(a), b * theta)
            scale = lcm(*(x.denominator for x in tied))
            tied = tuple(int(x * scale) for x in tied)
            endpoints = []
            for sign in (-1, 1):
                w = (tied[0], tied[1] + sign, tied[2])
                assert min(w) > 0 and dot(w, q)
                omitted = {max(pair, key=lambda x: dot(w, x)) for pair in pairs}
                B = frozenset(D - omitted)
                assert len(B) == n and core <= B and (B - core) <= set(ground)
                assert all(tuple(c - (i == j) for j, c in enumerate(x)) in B
                           for x in B for i, c in enumerate(x) if c)
                # This replay uses actual sparse quotient columns on the
                # relative ground; it is independent of the cheaper-pair rule.
                columns = {x: ({canonical_column.get(x, x): 1}
                               if x in D and canonical_column.get(x, x) not in core else {})
                           for x in ground}
                selected, _ = rank_trace(sorted(ground, key=lambda x: (dot(w, x), x)),
                                         columns.__getitem__, 101, set())
                assert selected == set(B - core)
                counts['independent_family_contracted_matrix_bases'] += 1
                states.add(B)
                endpoints.append(B)
            assert canonical(tuple(x-y for x, y in zip(profile(endpoints[0]), profile(endpoints[1])))) == q
            parallel_edges.add(frozenset(endpoints))
            rows.append({'theta': str(theta), 'weight_on_q_wall': tied,
                         'endpoint_profiles': [profile(B) for B in endpoints]})
        assert len(states) == 2*(m+1) and len(parallel_edges) == m+1
        fixtures.append({'m': m, 'length': n, 'primitive_direction': q,
                         'minimum_length': threshold(q), 'residual_rank': residual,
                         'all_positive_cells': len(states), 'parallel_wall_pieces': len(parallel_edges),
                         'socle_pairs': pairs, 'relative_ground_size': len(ground),
                         'all_q_wall_intervals': rows})
    return fixtures


def band_controls(counts):
    fixtures = []
    for n in range(2, 65):
        for R in sorted({1, min(3, n), n//2, n}):
            values = set()
            for a in range(1, n):
                for b in range(1, n-a+1):
                    if n-R < a+b <= n and gcd(a, b) == 1:
                        for i, j in combinations(range(3), 2):
                            q = [0]*3
                            q[i], q[j] = a, -b
                            values.add(tuple(q))
                for b in range(1, n):
                    for c in range(1, max(0, (n-a-b)//(b+1))+1):
                        if n-R < a+b*c+b+c <= n and gcd(a, b, c) == 1:
                            for j in range(3):
                                q = [-b, -c]
                                q.insert(j, a)
                                values.add(canonical(q))
            pairs = sum((n+1)//j - 1 for j in range(2, (n+1)//2+1))
            overcount = 3*R*(pairs+n)
            assert len(values) <= overcount
            assert all(1 <= n-threshold(q)+1 <= R for q in values)
            counts['independent_near_saturation_direction_bands'] += 1
            if n in (8, 16, 32, 64):
                fixtures.append({'length': n, 'residual_rank_cap': R,
                                 'primitive_band_directions': len(values),
                                 'divisor_pair_overcount': overcount})
    return fixtures


def polynomial_remainder(poly, monic):
    poly = list(poly)
    while len(poly) >= len(monic):
        value = poly[-1]
        offset = len(poly)-len(monic)
        for j, c in enumerate(monic):
            poly[offset+j] -= value*c
        assert not poly[-1]
        poly.pop()
    return tuple(poly)+((0,)*(len(monic)-1-len(poly)))


def polynomial_product(left, right):
    out = [0]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i+j] += a*b
    return tuple(out)


def planar_embedding_model(kind, size):
    if kind == 'box':
        s, t = size
        reference = tuple(product(range(s), range(t)))
        indices = {a: i for i, a in enumerate(reference)}
        evaluate = lambda a: {indices[a]: 1} if a in indices else {}
        generators = [{(s, 0): 1}, {(0, t): 1}]
        return reference, evaluate, (0,)*s+(1,), generators
    R = size
    f = (0,)*R+(1,)
    if kind == 'points':
        f = (1,)
        for j in range(1, R+1):
            f = polynomial_product(f, (-j, 1))
        h = (0,)*R+(1,)
    else:
        assert kind == 'reciprocal_jet'
        h = (0,)+(1,)*(R-1)
    h = polynomial_remainder(h, f)

    @lru_cache(None)
    def evaluate(a):
        value = (1,)+(0,)*(R-1)
        for _ in range(a[0]):
            value = polynomial_remainder((0,)+value, f)
        for _ in range(a[1]):
            value = polynomial_remainder(polynomial_product(value, h), f)
        return {i: c for i, c in enumerate(value) if c}
    generators = [{(i, 0): c for i, c in enumerate(f) if c},
                  {(0, 1): 1, **{(i, 0): -c for i, c in enumerate(h) if c}}]
    return tuple((i, 0) for i in range(R)), evaluate, f, generators


def complete_planar_fan(R, evaluate, p, counts):
    cells = []
    for B in finite_staircases(R, 2):
        selected, _ = field_basis(tuple(sorted(B)), evaluate, p)
        if len(selected) != R:
            continue
        normals, _ = border_normals(tuple(sorted(B)), evaluate, p)
        low, high, possible = Q(0), None, True
        for x, z in normals:
            if not x:
                possible &= z > 0
            elif x > 0:
                low = max(low, Q(-z, x))
            else:
                bound = Q(-z, x)
                high = bound if high is None else min(high, bound)
        if not possible or high is not None and high <= low:
            continue
        witness = low+1 if high is None else (low+high)/2
        assert witness > 0 and all(x*witness+z > 0 for x, z in normals)
        counts['independent_planar_border_cones'] += 1
        cells.append((B, normals, low, high, witness))
    cells.sort(key=lambda row: row[2])
    assert cells and cells[0][2] == 0 and cells[-1][3] is None
    assert all(left[3] == right[2] for left, right in zip(cells, cells[1:]))
    return cells


def planar_fan_embedding_controls(counts):
    fixtures = []
    cases = [(kind, R, 'balanced') for kind in ('points', 'reciprocal_jet') for R in range(2, 6)]
    cases += [('box', (2, 2), 'balanced'), ('box', (2, 3), 'balanced')]
    cases += [('points', 3, 'one_y'), ('reciprocal_jet', 4, 'one_y'),
              ('box', (2, 2), 'one_y'), ('points', 3, 'long_core'),
              ('points', 6, 'balanced'), ('points', 8, 'generic_counter')]
    for kind, size, tag_mode in cases:
        reference, evaluate, f, generators = planar_embedding_model(kind, size)
        R, a = len(reference), len(reference)+1
        if tag_mode == 'long_core':
            a = R*R+1
        b = 1 if tag_mode == 'one_y' else a+1
        F = (0,)*(a-len(f)+1)+f
        assert len(F) == a+1 and not F[0] and F[-1] == 1
        reference_H = tuple((j, 0, 0) for j in range(a+1))
        reference_H += tuple((0, j, 0) for j in range(1, b))
        reference_H += tuple((i, 0, j+1) for i, j in reference)
        indices = {alpha: i for i, alpha in enumerate(reference_H)}
        N = a+b+R
        assert len(indices) == N

        @lru_cache(None)
        def ev(alpha):
            i, j, k = alpha
            if j:
                if i or k or j > b:
                    return {}
                if j < b:
                    return {indices[alpha]: 1}
                return {indices[(t, 0, 0)]: c for t, c in enumerate(F) if c}
            if k:
                return {a+b+t: c for t, c in evaluate((i, k-1)).items()}
            reduced = polynomial_remainder((0,)*i+(1,), (0,)+F)
            return {indices[(t, 0, 0)]: c for t, c in enumerate(reduced) if c}

        matrices = [tuple(ev(tuple(x+(i == j) for j, x in enumerate(alpha)))
                          for alpha in reference_H) for i in range(3)]
        def apply(matrix, vector):
            out = defaultdict(int)
            for j, c in vector.items():
                for i, d in matrix[j].items():
                    out[i] += c*d
            return {i: c for i, c in out.items() if c}
        for i, j in combinations(range(3), 2):
            for k in range(N):
                assert apply(matrices[i], matrices[j][k]) == apply(matrices[j], matrices[i][k])
                counts['integer_planar_embedding_commuting_columns'] += 1
        @lru_cache(None)
        def from_matrices(alpha):
            if not any(alpha):
                return {indices[(0, 0, 0)]: 1}
            i = next(i for i, x in enumerate(alpha) if x)
            previous = tuple(x-(i == j) for j, x in enumerate(alpha))
            return apply(matrices[i], from_matrices(previous))
        E = tuple(hyperbolic_cross(N, 3))
        for alpha in E:
            assert ev(alpha) == from_matrices(alpha)
            counts['independent_planar_embedding_cyclic_columns'] += 1
        H_generators = [{(1, 1, 0): 1}, {(0, 1, 1): 1}, {(0, b+1, 0): 1},
                        {(i+1, 0, 0): c for i, c in enumerate(F) if c},
                        {(0, b, 0): -1, **{(i, 0, 0): c for i, c in enumerate(F) if c}}]
        H_generators += [{(i, 0, j+1): c for (i, j), c in g.items()} for g in generators]
        for g in H_generators:
            out = defaultdict(int)
            for alpha, c in g.items():
                for i, d in ev(alpha).items():
                    out[i] += c*d
            assert not any(out.values())
            counts['independent_planar_embedding_generator_annihilations'] += 1
        core = frozenset([(j, 0, 0) for j in range(a)]+[(0, j, 0) for j in range(1, b)])
        q = (a, -b, 0)
        for p in ((0,) if tag_mode == 'generic_counter' else (0, 2, 3, 5)):
            planar = complete_planar_fan(R, evaluate, p, counts)
            if kind == 'reciprocal_jet':
                assert len(planar) == R
            if tag_mode == 'generic_counter':
                assert R == 8 and len(planar) == 10 > R+1
            cells, rows = [], []
            for B, normals, low, high, witness in planar:
                for sign in (-1, 1):
                    keep = (a, 0, 0) if sign == -1 else (0, b, 0)
                    full = core | {keep} | {(i, 0, j+1) for i, j in B}
                    assert len(full) == N
                    price = witness, Q(a, b)*witness*(2 if sign == -1 else Q(1, 2)), Q(1)
                    direct, _ = rank_trace(sorted(E, key=lambda alpha: (dot(price, alpha), alpha)),
                                           ev, p, set())
                    assert direct == full
                    actual, _ = border_normals(tuple(sorted(full)), ev, p)
                    triangle = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
                    polygon, expected = triangle[:], triangle[:]
                    for normal in actual:
                        assert dot(price, normal) > 0
                        polygon = clip(polygon, normal, 1)
                    for normal in tuple((x, 0, z) for x, z in normals)+(tuple(sign*x for x in q),):
                        expected = clip(expected, normal, 1)
                    assert convex_hull([(v[0], v[1]) for v in polygon]) == convex_hull(
                        [(v[0], v[1]) for v in expected])
                    cells.append((full, polygon))
                    counts['independent_whole_embedded_cone_equalities'] += 1
                    rows.append({'planar_standard_set': sorted(B), 'base_side': sign,
                                 'strict_weight': list(map(str, price)),
                                 'full_standard_set': sorted(full)})
            assert sum(area2(poly) for B, poly in cells) == 1
            edges = positive_edges(cells)
            parallel = [row for row in edges.values() if row['direction'] == q]
            assert len(cells) == 2*len(planar) and len(parallel) == len(planar)
            counts['independent_complete_planar_fan_embeddings'] += 1
            fixtures.append({'planar_model': kind, 'planar_size': size, 'tag_mode': tag_mode, 'characteristic': p,
                             'planar_length': R, 'length': N, 'a': a, 'b': b,
                             'F_coefficients': F, 'residual_rank': R+1,
                             'planar_initial_ideals': len(planar), 'initial_ideals': len(cells),
                             'q_parallel_walls': len(parallel), 'all_cones': rows})
    return fixtures


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    counts = defaultdict(int)
    downsets = downset_controls(counts)
    fans = saved_fan_controls(counts)
    assert counts['independent_outside_wall_contraction_cost_gaps'] > 0
    family = repeated_wall_family(counts)
    bands = band_controls(counts)
    embeddings = planar_fan_embedding_controls(counts)
    sources = ('docs/wall_coexistence.md', 'scripts/run_wall_coexistence_checks.py',
               'docs/primitive_wall_spectrum.md', 'scripts/run_primitive_wall_checks.py',
               'docs/rational_shadow_bound.md', 'scripts/run_extremal_initial_ideals_checks.py',
               'scripts/run_oracle_enumeration_checks.py', 'scripts/run_signed_face_checks.py',
               'scripts/run_higher_interaction_checks.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_planar_dyadic_bound_checks.py',
               'docs/planar_dyadic_bound.md', 'docs/liaison_fan_composition.md',
               'docs/extremal_initial_ideals.md',
               'results/research/oracle_enumeration_checks.json',
               'results/research/extremal_initial_ideals_checks.json')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_wall_coexistence_checks.py',
              'counts': dict(counts), 'all_small_downset_pairs': downsets,
              'actual_quotient_contractions': fans, 'growing_actual_parallel_wall_family': family,
              'near_saturation_bands': bands,
              'all_planar_fan_embeddings': embeddings,
              'uniform_results': {
                  'compulsory_core': 'Every endpoint of every q-parallel wall contains C_q=Div(q_plus) union Div(q_minus) minus the two tips.',
                  'single_contraction': 'Rank r=n-nu(q)+1 on E(C_q,r), contained in at most six translates of V_r^3.',
                  'local_contraction_scope': 'Restricted and global support functions agree locally at each q-wall, where both adjacent standard sets contain C_q. Forcing C_q can strictly raise the minimum at other positive weights.',
                  'parallel_wall_bound': 'O(min(r^(4/3)*L_r^2, r*(1+||q||_1^(1/3))*L_r)), L_r=1+log(r+1).',
                  'saturated_wall_bound': 'At most one q-wall when nu(q)=n.',
                  'near_saturation_cell_bound': 'O(n*R^(7/3)*log(n+1)*L_R^2) cells incident to a wall with n-R<nu(q)<=n.',
                  'sublinear_example': 'R=floor(n^(3/4)) gives O(n^(11/4)*log(n+1)^3) for that part of the fan.',
                  'near_saturation_wall_cost': 'V_R=O(n^2*R^(7/3)*log(n+1)*L_R^2); R=floor(n^(3/7)) gives O(n^3*log(n+1)^3) for that part only.',
                  'two_sided_cell_localization': 'O(n^2*log(n+1)^3) cells are incident to a wall with nu(q)<=floor(n^(3/7)) or nu(q)>n-floor(n^(3/7)). A polynomial excess must have every incident cost in the middle.',
                  'two_sided_cubic_budget_localization': 'The portions of V with nu(q)<=floor(n^(3/5)) and nu(q)>n-floor(n^(3/7)) both obey O(n^3*log(n+1)^3). No bound on the complementary middle portion is proved.',
                  'actual_lower_family': 'n=a+b+2*m^2+m, (a,b)=(m^4+1,m^4+2), r=2*m^2+m+1, p_I(a,-b,0)=m+1.',
                  'arbitrary_planar_embedding': 'For any length-R planar J, H=(xy,yz,xF,y^(b+1),F-y^b,zJ), with monic F in J, F(0)=0, deg F=a>=R+1 and gcd(a,b)=1, has length a+b+R, g(H)=2*g(J), and p_H(a,-b,0)=g(J). Its initial formula holds for every global order.',
                  'superlinear_residual_obstruction': 'Planar generic ideals give p_H(q)=Theta(r*log(r+1)), r=R+1; a=R^2+1, b=R^2+2 makes n-nu(q)=R while n=Theta(R^2). Thus a uniform O(r) parallel-wall bound is false.'},
              'new_unrestricted_G3_count': False, 'polynomial_lower_separation_proved': False,
              'novelty_established': False, 'goal_complete': False,
              'scope': 'Uniform proofs use the fixed compulsory monomial core, a relative divisor cover, one rank-preserving matroid contraction, the published Dey bound and the translated arithmetic shadow estimate. The tagged construction embeds every planar fan over every field; its superlinear residual lower family uses the established rational generic planar count. Whole embedded cones are checked independently from quotient columns, including local, nonlocal and characteristic-dependent cases. The middle collective cost budget, an improved unrestricted exponent and exact priority remain unproved.',
              'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sources}}
    (ROOT/'results/research/wall_coexistence_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': report['status'], 'counts': report['counts'],
                      'actual_fans': len(fans), 'growing_family_fixtures': len(family)}, indent=2))


if __name__ == '__main__':
    main()
