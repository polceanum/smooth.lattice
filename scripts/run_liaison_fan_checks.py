#!/usr/bin/env python3
"""Exact finite controls for classical liaison addition and fan budgets."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations, product
from math import gcd, prod
from pathlib import Path
import hashlib
import json

from run_extremal_initial_ideals_checks import adjacent_weights, integer_weight
from run_higher_seed_checks import column_rank
from run_higher_interaction_checks import ground, linear_sum, rank_trace
from run_hypersurface_truncation_checks import poly_mul, relation_rows
from run_interleaving_research import minimal_excluded_border
from run_signed_face_checks import normals
from run_structural_complexity_checks import area2, clip, cross, dot, hyperbolic_cross, primitive

ROOT = Path(__file__).resolve().parents[1]


def rectangle_graph(q, free, F):
    coordinates = [i for i in range(3) if i != free]
    cells = tuple(product(range(q), repeat=2))
    index = {v: i for i, v in enumerate(cells)}
    basis = tuple(tuple(0 if i == free else v[coordinates.index(i)] for i in range(3))
                  for v in cells)

    @lru_cache(None)
    def power(e):
        if not e:
            return {(0, 0): 1}
        out = defaultdict(int)
        for (a, b), c in power(e-1).items():
            for (i, j), d in F.items():
                if a+i < q and b+j < q:
                    out[a+i, b+j] += c*d
        return {v: c for v, c in out.items() if c}

    def ev(v):
        i, j = v[coordinates[0]], v[coordinates[1]]
        return {index[a+i, b+j]: c for (a, b), c in power(v[free]).items()
                if a+i < q and b+j < q}

    generators = [{tuple(q if i == c else 0 for i in range(3)): 1}
                  for c in coordinates]
    relation = {tuple(int(i == free) for i in range(3)): 1}
    for v, c in F.items():
        exponent = tuple(0 if i == free else v[coordinates.index(i)] for i in range(3))
        relation[exponent] = -c
    generators.append(relation)
    return basis, ev, generators


def glue_model(a, b, c, first, second):
    base = tuple(product(range(a), range(b), range(c)))
    B1, ev1, gen1 = first
    B2, ev2, gen2 = second
    tag1, tag2 = (a, 0, 0), (0, b, 0)
    shifted = lambda v, w: tuple(x+y for x, y in zip(v, w))
    basis = base+tuple(shifted(v, tag1) for v in B1)+tuple(shifted(v, tag2) for v in B2)
    index = {v: i for i, v in enumerate(base)}
    assert len(set(basis)) == len(basis) == a*b*c+len(B1)+len(B2)

    def ev(v):
        x, y, z = v
        if z >= c or (x >= a and y >= b):
            return {}
        if x >= a:
            return {len(base)+i: d for i, d in ev1((x-a, y, z)).items()}
        if y >= b:
            return {len(base)+len(B1)+i: d for i, d in ev2((x, y-b, z)).items()}
        return {index[v]: 1}

    generators = [{(0, 0, c): 1}, {(a, b, 0): 1}]
    generators += [{shifted(v, tag1): d for v, d in f.items()} for f in gen1]
    generators += [{shifted(v, tag2): d for v, d in f.items()} for f in gen2]
    assert all(ev(v) == {i: 1} for i, v in enumerate(basis))
    return basis, ev, generators


def model_controls():
    counts, fixtures = defaultdict(int), []
    weights = list(permutations((1, 11, 121)))+[(5, 7, 11), (11, 5, 7), (7, 11, 5)]
    for q in (2, 3):
        F = {(1, 0): 1, (0, 1): 2, (1, 1): 1}
        G = {(1, 0): 2, (0, 1): 1, (1, 1): 3}
        first, second = rectangle_graph(q, 0, F), rectangle_graph(q, 1, G)
        basis, ev, generators = glue_model(q, q, q, first, second)
        E = ground(3*q-2)
        for a in E:
            for component in (first, second):
                B, eval_branch, _ = component
                assert all(sum(B[i]) >= sum(a) for i in eval_branch(a))
                counts['branch_filtration_column_terms_checked'] += len(eval_branch(a))
        maps = [[ev(tuple(x+(i == j) for j, x in enumerate(v))) for v in basis]
                for i in range(3)]
        for i, j in combinations(range(3), 2):
            for k in range(len(basis)):
                assert linear_sum((d, maps[i][h]) for h, d in maps[j][k].items()) == \
                       linear_sum((d, maps[j][h]) for h, d in maps[i][k].items())
                counts['cyclic_commuting_multiplication_columns'] += 1
        # Every defining relation annihilates every cyclic reference vector.
        for f, v in product(generators, basis):
            assert not linear_sum((d, ev(tuple(x+y for x, y in zip(a, v)))) for a, d in f.items())
            counts['integer_generator_membership_replays'] += 1
        for p in (0, 2, 3, 5):
            rows = relation_rows(E, generators, p)
            assert len(E)-len(rows) == len(basis)
            for row in rows:
                image = linear_sum((d, ev(a)) for a, d in zip(E, row) if d)
                assert all(not (d % p if p else d) for d in image.values())
                counts['independent_macaulay_kernel_rows'] += 1
            for w in weights:
                ordered = sorted(E, key=lambda a: (dot(a, w), a))
                B, _ = rank_trace(ordered, ev, p, set())
                branches = []
                for ref, eval_branch, _ in (first, second):
                    D, _ = rank_trace(ordered, eval_branch, p, set())
                    assert len(D) == len(ref)
                    assert sum(map(sum, D)) <= q*q*(q-1)
                    branches.append(D)
                proposed = set(product(range(q), repeat=3))
                proposed |= {(a+q, b, c) for a, b, c in branches[0]}
                proposed |= {(a, b+q, c) for a, b, c in branches[1]}
                assert B == proposed
                moment = tuple(sum(v[i] for v in B) for i in range(3))
                fixed = tuple(sum(v[i] for v in product(range(q), repeat=3))+
                              q**3*int(i in (0, 1)) for i in range(3))
                varying = tuple(moment[i]-fixed[i] for i in range(3))
                assert sum(varying) == sum(sum(map(sum, D)) for D in branches)
                assert sum(varying) <= 2*q*q*(q-1)
                counts['independent_initial_set_and_moment_equalities'] += 1
            counts['independent_full_kernel_models'] += 1
        fixtures.append({'side': q, 'branch_lengths': [q*q, q*q], 'box_length': q**3,
                         'glued_length': len(basis), 'sum_branch_moments': 2*q*q*(q-1),
                         'common_refinement_is_product_claimed': False})
    return dict(counts), fixtures


def full_fan(E, ev, p):
    active = tuple(v for v in E if ev(v))
    comparisons = sorted({primitive(tuple(x-y for x, y in zip(v, w)))
                          for v, w in combinations(active, 2)
                          if min(x-y for x, y in zip(v, w)) < 0 < max(x-y for x, y in zip(v, w))})
    pending, visited = [(1, 11, 121)], {}
    while pending:
        w = pending.pop()
        B, _ = rank_trace(sorted(active, key=lambda a: (dot(a, w), a)), ev, p, set())
        key = frozenset(B)
        if key in visited:
            continue
        Z, _ = normals(B, ev, p)
        polygon = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for z in Z:
            polygon = clip(polygon, z, 1)
            assert polygon
        assert area2(polygon) > 0 and all(dot(z, w) > 0 for z in Z)
        visited[key] = polygon
        pending.extend(adjacent_weights(polygon, Z, comparisons))
    assert sum(map(area2, visited.values())) == 1
    return visited


def identical_fan_controls():
    counts, fixtures = defaultdict(int), []
    branch = rectangle_graph(2, 2, {(1, 0): 1, (0, 1): 1, (1, 1): 1})
    basis, ev, _ = glue_model(2, 2, 3, branch, branch)
    for p, expected in ((0, 6), (2, 4), (3, 6)):
        first = full_fan(ground(2), branch[1], p)
        joined = full_fan(ground(4), ev, p)
        assert len(first) == len(joined) == expected
        for D, poly in first.items():
            B = frozenset(set(product(range(2), range(2), range(3))) |
                          {(a+2, b, c) for a, b, c in D} |
                          {(a, b+2, c) for a, b, c in D})
            assert B in joined and joined[B] == poly
            counts['complete_identical_fan_cone_equalities'] += 1
        fixtures.append({'characteristic': p, 'branch_length': 4, 'glued_length': 20,
                         'branch_initials': expected, 'glued_initials': expected,
                         'product_of_counts_is_actual': False,
                         'branch_standard_sets_and_cones': [
                             {'standard_set': sorted(D),
                              'weight_simplex_polygon': [[str(x) for x in v] for v in poly]}
                             for D, poly in sorted(first.items(), key=lambda item: sorted(item[0]))]})
    return dict(counts), fixtures


def nonmonomial_point_controls():
    counts = defaultdict(int)
    grid = tuple(product((-1, 1), repeat=2))
    first = [(value, y, z) for value, (y, z) in zip(range(2, 6), grid)]
    second = [(x, value, z) for value, (x, z) in zip(range(6, 10), grid)]
    base = list(product((-1, 1), repeat=3))
    points = base+first+second
    assert len(set(points)) == 16
    f, g, h = ({tuple(2*int(i == j) for i in range(3)): 1, (0, 0, 0): -1}
               for j in range(3))

    def graph(free, values):
        other = [i for i in range(3) if i != free]
        polynomial = {tuple(int(i == free) for i in range(3)): Q(1)}
        for v in product((0, 1), repeat=2):
            c = sum(Q(value*s**v[0]*t**v[1], 4)
                    for value, (s, t) in zip(values, grid))
            exponent = tuple(0 if i == free else v[other.index(i)] for i in range(3))
            polynomial[exponent] = -c
        return {v: c for v, c in polynomial.items() if c}

    J1, J2 = (g, h, graph(0, range(2, 6))), (f, h, graph(1, range(6, 10)))
    generators = [h, poly_mul(f, g, 0)]
    generators += [poly_mul(f, v, 0) for v in J1]
    generators += [poly_mul(g, v, 0) for v in J2]
    for polynomial in generators:
        for point in points:
            assert sum(c*prod(x**e for x, e in zip(point, v))
                       for v, c in polynomial.items()) == 0
            counts['nonmonomial_tag_point_membership_replays'] += 1
    for p in (0, 101):
        def model(P):
            return lambda v: {i: prod(pow(x, e, p) if p else x**e for x, e in zip(point, v))
                              for i, point in enumerate(P)}
        E = tuple(hyperbolic_cross(16, 3))
        ev, ev1, ev2 = model(points), model(first), model(second)
        assert column_rank([ev(v) for v in E], p) == 16
        for w in permutations((1, 17, 289)):
            order = sorted(E, key=lambda v: (dot(v, w), v))
            B, _ = rank_trace(order, ev, p, set())
            D1, _ = rank_trace(order, ev1, p, set())
            D2, _ = rank_trace(order, ev2, p, set())
            expected = set(product(range(2), repeat=3))
            expected |= {(x+2, y, z) for x, y, z in D1}
            expected |= {(x, y+2, z) for x, y, z in D2}
            assert len(D1) == len(D2) == 4 and B == expected
            counts['nonmonomial_tag_independent_initial_sets'] += 1
    return dict(counts), {'field_characteristics': [0, 101], 'univariate_tags': ['x^2-1', 'y^2-1', 'z^2-1'],
                          'box_length': 8, 'branch_lengths': [4, 4], 'glued_length': 16,
                          'points': points}


def box_wall_controls():
    counts, fixtures = defaultdict(int), []
    directions = sorted({primitive(v) for v in product(range(-3, 4), repeat=3)
                         if min(v) < 0 < max(v)})
    for v in directions:
        positive, negative = sum(max(x, 0) for x in v), sum(max(-x, 0) for x in v)
        w0 = tuple(negative if x > 0 else positive if x < 0 else 1 for x in v)
        w1 = cross(v, w0)
        assert min(w0) > 0 and dot(w0, v) == dot(w1, v) == 0
        low = tuple(max(-x, 0) for x in v)
        T = dot(w0, low), dot(w1, low)
        box = tuple(product(*(range(abs(x)+1) for x in v)))
        values = [(dot(w0, a), dot(w1, a)) for a in box]
        assert values.count(T) == 2
        assert sum(x < T for x in values) == sum(x > T for x in values) == (len(box)-2)//2
        assert len(box) == prod(abs(x)+1 for x in v) == 2*(sum(x < T for x in values)+1)
        counts['sharp_integer_box_reflection_witnesses'] += 1
        if v in ((1, -1, 0), (3, -2, 1)):
            fixtures.append({'primitive_direction': v, 'symbolic_positive_weight': [w0, w1],
                             'box_cardinality': len(box), 'selected_cardinality': len(box)//2})
    assert len(fixtures) == 2
    # Replay actual cut cones of bounded boxes against the product budget.
    for m, caps in ((3, (2, 2, None)), (5, (2, None, 3)), (7, (None, 3, 2))):
        E = tuple(v for v in hyperbolic_cross(m, 3)
                  if all(c is None or x < c for x, c in zip(v, caps)))
        pending, visited = [frozenset(sorted(E, key=lambda a: (dot(a, (1, 11, 121)), a))[:m])], {}
        while pending:
            B = pending.pop()
            if B in visited:
                continue
            maximal = [v for v in B if not any(tuple(x+(i == j) for i, x in enumerate(v)) in B
                                               for j in range(3))]
            border = [v for v in minimal_excluded_border(tuple(B))
                      if all(c is None or x < c for x, c in zip(v, caps))]
            Z = sorted({tuple(x//gcd(*z) for x in z)
                        for v, w in product(border, maximal)
                        if (z := tuple(x-y for x, y in zip(v, w)))})
            poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
            for z in Z:
                poly = clip(poly, z, 1)
                assert poly
            assert area2(poly) > 0
            visited[B] = poly
            for p, q in zip(poly, poly[1:]+poly[:1]):
                if any(p[i] == q[i] == 0 for i in range(3)):
                    continue
                z = next(z for z in Z if dot(z, p) == dot(z, q) == 0)
                assert prod(abs(x)+1 for x in z) <= 2*m
                w = tuple((x+y)/2 for x, y in zip(p, q))
                neighbor = frozenset(sorted(E, key=lambda a: (dot(a, w), -dot(a, z), a))[:m])
                assert neighbor != B
                pending.append(neighbor)
                counts['actual_box_cut_wall_product_checks'] += 1
        assert sum(map(area2, visited.values())) == 1
        counts['complete_box_cut_fans'] += 1
    return dict(counts), fixtures


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    model, local_fixtures = model_controls()
    identical, fan_fixtures = identical_fan_controls()
    point, point_fixture = nonmonomial_point_controls()
    box, box_fixtures = box_wall_controls()
    sources = ('docs/liaison_fan_composition.md', 'scripts/run_liaison_fan_checks.py',
               'docs/extremal_initial_ideals.md', 'docs/generic_multiaffine_bound.md',
               'scripts/run_extremal_initial_ideals_checks.py', 'scripts/run_higher_seed_checks.py',
               'scripts/run_higher_interaction_checks.py', 'scripts/run_hypersurface_truncation_checks.py',
               'scripts/run_signed_face_checks.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_socle_stability_checks.py')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_liaison_fan_checks.py',
              'counts': {**model, **identical, **point, **box}, 'local_composition_fixtures': local_fixtures,
              'complete_identical_fan_fixtures': fan_fixtures, 'sharp_box_wall_fixtures': box_fixtures,
              'nonmonomial_tag_point_fixture': point_fixture,
              'uniform_results': {
                  'operation_and_additive_length': 'Classical liaison addition in a finite affine chart; length abc+m1+m2.',
                  'all_global_order_initial_formula': '(L3), with exact common refinement and disjoint recovery of both branch standard sets.',
                  'unrestricted_product_of_branch_counts': False,
                  'box_cut_refinement_bound': 'O(M^2*log(M+1)^4), independently of the number of input fans.',
                  'local_moment_bound': 'O((1+S1+S2)^(3/2)), with the fixed box/tag profiles removed as translations.',
                  'balanced_arbitrary_nonlinear_rectangle_gluing': 'N=q^3+2*q^2; g(H)=O(N^(3/2)) over all fields.'},
              'unrestricted_G3_improvement': False, 'novelty_established': False, 'goal_complete': False,
              'scope': 'A classical construction, exact fan-composition tool and two uniform feasibility exclusions. The finite checks do not supply a larger asymptotic lower family or a priority certificate.',
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/'results/research/liaison_fan_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
