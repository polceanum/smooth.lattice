#!/usr/bin/env python3
"""Exact quotient fans and coordinate-projection tag controls; no growth fit."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import permutations, product
from math import comb, gcd
from pathlib import Path
import hashlib
import json

from run_interleaving_research import minimal_excluded_border
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_planar_dyadic_bound_checks import convex_hull
from run_primitive_wall_checks import canonical, threshold
from run_signed_face_checks import field_basis, normals as border_normals
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross
from run_wall_coexistence_checks import complete_planar_fan

ROOT = Path(__file__).resolve().parents[1]
TRIANGLE = tuple(tuple(Q(i == j) for i in range(3)) for j in range(3))


def polygon_key(poly):
    return tuple(convex_hull([(v[0], v[1]) for v in poly]))


def fan(n, evaluate, p, counts):
    """Traverse the actual border cones using residue-column elimination."""
    E = tuple(hyperbolic_cross(n, 3))
    queue = [((1, n+1, (n+1)**2), None)]
    visited, edges, queries = {}, set(), 0
    while queue:
        w, previous = queue.pop()
        B, _ = field_basis(sorted(E, key=lambda a: (dot(w, a), a)), evaluate, p)
        B = frozenset(B)
        assert len(B) == n
        queries += 1
        if previous is not None:
            assert B != previous
            edges.add(frozenset((B, previous)))
        if B in visited:
            continue
        normals, relations = border_normals(tuple(sorted(B)), evaluate, p)
        assert all(dot(q, w) > 0 for q in normals)
        poly = list(TRIANGLE)
        for q in normals:
            poly = clip(poly, q, 1)
        assert area2(poly) > 0
        visited[B] = poly
        counts['direct_field_border_relations'] += len(relations)
        for a, b in zip(poly, poly[1:]+poly[:1]):
            if any(a[i] == b[i] == 0 for i in range(3)):
                continue
            q = next(q for q in normals if dot(q, a) == dot(q, b) == 0)
            neighbor, _ = exact_adjacent_weight(n, a, b, q)
            queue.append((neighbor, B))
    assert sum(area2(poly) for poly in visited.values()) == 1
    assert queries == 1+2*len(edges)
    counts['complete_direct_quotient_fans'] += 1
    return visited


def projection(n, evaluate, p, axis, counts):
    indices = tuple(i for i in range(3) if i != axis)
    def embed(a):
        out = [0]*3
        for i, b in zip(indices, a):
            out[i] = b
        return tuple(out)
    ev = lru_cache(None)(lambda a: evaluate(embed(a)))
    basis, _ = field_basis(tuple(hyperbolic_cross(n, 2)), ev, p)
    m = len(basis)
    cells = complete_planar_fan(m, ev, p, counts)
    return m, cells, embed


def tag(n, evaluate, p, axis, counts, k=1):
    """Columns of the module B[z]/z^k direct-sum z^k A.

    This model includes nonunits. Independent CRT/point columns below
    check it in unit-coordinate and radical cases.
    """
    m, cells, embed = projection(n, evaluate, p, axis, counts)
    @lru_cache(None)
    def ev(a):
        b = list(a)
        if a[axis] < k:
            block = a[axis]
            b[axis] = 0
        else:
            block = k
            b[axis] -= k
        return {(block, j): c for j, c in evaluate(tuple(b)).items()}
    return n+k*m, ev, cells, embed


def predicted_fan(parent, cells, embed, axis, k=1):
    out = {}
    for B, poly in parent.items():
        for D, normals, *_ in cells:
            patch = list(poly)
            for q in normals:
                patch = clip(patch, embed(q), 1)
            if not patch or area2(patch) <= 0:
                continue
            low = set()
            for a in D:
                for j in range(k):
                    b = list(embed(a)); b[axis] = j
                    low.add(tuple(b))
            high = set()
            for a in B:
                b = list(a); b[axis] += k
                high.add(tuple(b))
            assert not low & high
            basis = frozenset(low | high)
            assert basis not in out
            out[basis] = patch
    assert sum(area2(poly) for poly in out.values()) == 1
    return out


def compare_fans(actual, predicted, counts):
    assert set(actual) == set(predicted)
    for B in actual:
        assert polygon_key(actual[B]) == polygon_key(predicted[B])
        counts['independent_entire_cone_equalities'] += 1


def hit_formula(parent, cells, embed):
    walls = [embed((lo.denominator, -lo.numerator))
             for _, _, lo, _, _ in cells if lo]
    hits = [sum(min(dot(q, v) for v in poly) < 0 <
                max(dot(q, v) for v in poly) for poly in parent.values())
            for q in walls]
    return walls, hits


def point_model(points, p):
    @lru_cache(None)
    def ev(a):
        out = {}
        for i, point in enumerate(points):
            value = 1
            for x, e in zip(point, a):
                value *= pow(x, e, p) if p else x**e
                if p:
                    value %= p
            if value:
                out[i] = value
        return out
    return len(points), ev


def local_model(n, p):
    # K[t]/t^n, x=t, y=t+t^2, z=t^2, including nonunit tags.
    @lru_cache(None)
    def ev(a):
        r = a[0]+a[1]+2*a[2]
        return {r+j: (comb(a[1], j) % p if p else comb(a[1], j))
                for j in range(a[1]+1) if r+j < n and
                (comb(a[1], j) % p if p else comb(a[1], j))}
    return n, ev


def monic_local_tag_model(n, k, p):
    # Independent polynomial multiplication in
    # K[x,y,z]/(x^n, y-x-x^2, z^(k+1)-x^2*z^k).
    # The basis is x^i*z^j, i<n, j<=k. No occurrence of the
    # low/high tag decomposition is used in this evaluator.
    generators = (((1, 0),), ((1, 0), (2, 0)), ((0, 1),))
    @lru_cache(None)
    def ev(a):
        value = {(0, 0): 1}
        for power, terms in zip(a, generators):
            for _ in range(power):
                out = defaultdict(int)
                for (i, j), c in value.items():
                    for di, dj in terms:
                        u, v = i+di, j+dj
                        if v > k:
                            u += 2*(v-k)
                            v = k
                        if u < n:
                            out[u, v] += c
                value = {b: (c % p if p else c) for b, c in out.items()
                         if (c % p if p else c)}
        return {(j, i): c for (i, j), c in value.items()}
    return ev


def finite_controls(counts):
    rows = []
    for p, points in ((0, ((1, 1, 1), (2, 3, 4))),
                      (5, ((1, 1, 1), (2, 3, 4))),
                      (0, ((1, 1, 1), (2, 3, 4), (3, 2, 5)))):
        n, ev = point_model(points, p)
        parent = fan(n, ev, p, counts)
        for axis in range(3):
            N, tagged, cells, embed = tag(n, ev, p, axis, counts)
            projected = tuple(sorted({tuple(0 if j == axis else x
                                           for j, x in enumerate(point))
                                      for point in points}))
            actual_N, actual_ev = point_model(tuple(points)+projected, p)
            assert actual_N == N
            actual = fan(N, actual_ev, p, counts)
            predicted = predicted_fan(parent, cells, embed, axis)
            compare_fans(actual, predicted, counts)
            compare_fans(fan(N, tagged, p, counts), actual, counts)
            walls, hits = hit_formula(parent, cells, embed)
            assert len(actual) == len(parent)+sum(hits)
            counts['exact_pencil_cell_gain_formulas'] += 1
            for other in range(3):
                _, before, _ = projection(n, ev, p, other, counts)
                _, after, _ = projection(N, tagged, p, other, counts)
                assert [(row[2], row[3]) for row in before] == [
                    (row[2], row[3]) for row in after]
                counts['entire_projection_fan_invariances'] += 1
            rows.append({'kind': 'points', 'characteristic': p,
                         'points': points, 'tag_axis': axis, 'length': N,
                         'parent_g': len(parent), 'tagged_g': len(actual),
                         'projection_walls': walls, 'strict_cell_hits': hits})
    for p, n, k in ((0, 3, 1), (2, 3, 1), (3, 4, 2)):
        _, ev = local_model(n, p)
        parent = fan(n, ev, p, counts)
        N, tagged, cells, embed = tag(n, ev, p, 2, counts, k)
        monic = monic_local_tag_model(n, k, p)
        for a in hyperbolic_cross(N+1, 3):
            assert monic(a) == tagged(a)
            counts['independent_monic_polynomial_tag_normal_forms'] += 1
        actual = fan(N, monic, p, counts)
        compare_fans(actual, predicted_fan(parent, cells, embed, 2, k), counts)
        rows.append({'kind': 'nilpotent_curve', 'characteristic': p,
                     'parent_length': n, 'tag_power': k, 'length': N,
                     'parent_g': len(parent), 'tagged_g': len(actual)})
    return rows


def completion_controls(counts):
    rows = []
    for p in (0, 5):
        points = ((1, 1, 1), (2, 3, 4))
        n, ev = point_model(points, p)
        face_points = tuple(sorted({tuple(point[i] if mask >> i & 1 else 0
                                           for i in range(3))
                                    for point in points for mask in range(8)}))
        N, direct = point_model(face_points, p)
        actual = fan(N, direct, p, counts)
        parent = fan(n, ev, p, counts)
        for order in permutations(range(3)):
            size, model, predicted = n, ev, parent
            for axis in order:
                size, model, cells, embed = tag(size, model, p, axis, counts)
                predicted = predicted_fan(predicted, cells, embed, axis)
            assert size == N <= 8*n
            compare_fans(actual, predicted, counts)
            # Full-ground rank replays use the iterated tag columns,
            # independently of the direct evaluation of face points.
            E = tuple(hyperbolic_cross(N, 3))
            for B, poly in actual.items():
                w = tuple(sum(v[i] for v in poly) for i in range(3))
                selected, _ = field_basis(sorted(E, key=lambda a: (dot(w, a), a)), model, p)
                assert frozenset(selected) == B
                counts['independent_iterated_tag_full_ground_replays'] += 1
            # A repeated tag changes length but not a single cone.
            _, repeated, embed = projection(size, model, p, order[-1], counts)
            refined = predicted_fan(predicted, repeated, embed, order[-1])
            assert sorted(map(polygon_key, refined.values())) == sorted(map(polygon_key, actual.values()))
            counts['complete_repeated_tag_stabilizations'] += 1
        rows.append({'characteristic': p, 'points': points, 'completion_length': N,
                     'initial_ideals': len(actual), 'all_tag_orders': 6})
    # The zero-coordinate case is nonreduced: I=(x,y,z) completes to
    # (x^2,y^2,z^2), not the one-point ideal of its set-theoretic faces.
    n, ev = point_model(((0, 0, 0),), 0)
    for axis in range(3):
        n, ev, _, _ = tag(n, ev, 0, axis, counts)
    assert n == 8
    box = frozenset(product(range(2), repeat=3))
    selected, _ = field_basis(tuple(sorted(box)), ev, 0)
    assert frozenset(selected) == box
    assert all(not ev(tuple(2 if i == j else 0 for i in range(3))) for j in range(3))
    counts['nonreduced_coordinate_face_negative_controls'] += 1
    return rows


def digit_case(n, k, counts):
    M = 1
    while M < n:
        M *= 2
    mask, cols = (1 << n)-1, [1]
    for _ in range(1, M):
        cols.append((cols[-1] ^ (cols[-1] << 1)) & mask)
    for r, col in enumerate(cols):
        assert col == sum((comb(r, j) & 1) << j for j in range(n))
        counts['independent_binomial_vs_recurrence_columns'] += 1
    E = tuple(hyperbolic_cross(n, 3))
    column = lambda a: cols[(a[0]-a[1]-k*a[2]) % M]
    def greedy(w):
        piv, B = {}, []
        for a in sorted(E, key=lambda a: (dot(w, a), a)):
            v = column(a)
            while v:
                j = (v & -v).bit_length()-1
                if j not in piv:
                    piv[j] = v; B.append(a); break
                v ^= piv[j]
        assert len(B) == n
        return frozenset(B)
    def normals(B):
        order, piv = tuple(sorted(B)), {}
        for i, a in enumerate(order):
            v, co = column(a), 1 << i
            while v:
                j = (v & -v).bit_length()-1
                if j not in piv:
                    piv[j] = (v, co); break
                row, rep = piv[j]; v ^= row; co ^= rep
            assert v
        out = set()
        for a in minimal_excluded_border(order):
            v, co = column(a), 0
            while v:
                j = (v & -v).bit_length()-1
                row, rep = piv[j]; v ^= row; co ^= rep
            for i, b in enumerate(order):
                if co >> i & 1:
                    q = tuple(x-y for x, y in zip(a, b)); scale = gcd(*q)
                    out.add(tuple(x//scale for x in q))
        return tuple(sorted(out))
    queue, visited, edges, queries = [((1, n+1, (n+1)**2), None)], {}, set(), 0
    while queue:
        w, previous = queue.pop()
        B = greedy(w); queries += 1
        if previous is not None:
            assert B != previous
            edges.add(frozenset((B, previous)))
        if B in visited:
            continue
        N, poly = normals(B), list(TRIANGLE)
        assert all(dot(q, w) > 0 for q in N)
        for q in N:
            poly = clip(poly, q, 1)
        assert area2(poly) > 0
        visited[B] = poly
        for a, b in zip(poly, poly[1:]+poly[:1]):
            if any(a[i] == b[i] == 0 for i in range(3)):
                continue
            q = next(q for q in N if dot(q, a) == dot(q, b) == 0)
            neighbor, _ = exact_adjacent_weight(n, a, b, q)
            queue.append((neighbor, B))
    assert sum(map(area2, visited.values())) == 1
    assert queries == 1+2*len(edges)
    for left, right in map(tuple, edges):
        q = canonical(tuple(sum(a[i] for a in left)-sum(a[i] for a in right)
                            for i in range(3)))
        assert threshold(q) <= n
    hits = [sum(min(dot(q, v) for v in poly) < 0 < max(dot(q, v) for v in poly)
                for poly in visited.values())
            for j in range(1, n) for q in [(j, -(n-j), 0)]]
    if n <= 10:
        @lru_cache(None)
        def ev(a):
            v = column(a)
            return {j: 1 for j in range(n) if v >> j & 1}
        compare_fans(fan(n, ev, 2, counts), visited, counts)
        # Independent CRT model: z is a unit, so H=I intersect (J,z).
        @lru_cache(None)
        def crt(a):
            out = {(1, j): c for j, c in ev(a).items()}
            if not a[2]:
                out.update({(0, j): c for j, c in ev(a).items()})
            return out
        actual = fan(2*n, crt, 2, counts)
        _, _, cells, embed = tag(n, ev, 2, 2, counts)
        assert len(cells) == n
        compare_fans(actual, predicted_fan(visited, cells, embed, 2), counts)
        assert len(actual) == len(visited)+sum(hits)
        counts['independent_nonreduced_CRT_entire_tag_fans'] += 1
    counts['complete_prime_digit_parent_fans'] += 1
    witnesses = [{'basis': sorted(B), 'cone': [list(map(str, v)) for v in poly]}
                 for B, poly in sorted(visited.items(), key=lambda row: sorted(row[0]))]
    return {'n': n, 'k': k, 'characteristic': 2, 'pascal_period': M,
            'parent_g': len(visited), 'tagged_length': 2*n,
            'tagged_g': len(visited)+sum(hits), 'gain': sum(hits),
            'gain_divided_by_n_squared': str(Q(sum(hits), n*n)),
            'strict_hook_cell_hits': hits, 'parent_cones': witnesses}


def write_report(report, target):
    # Keep each exact cone witness on one line; expanding its exponent
    # triples individually otherwise obscures the small summary in a
    # quarter-million-line diff. The result remains ordinary JSON.
    summary = {k: v for k, v in report.items() if k != 'prime_digit_cases'}
    text = json.dumps(summary, indent=2)[:-2]+',\n  "prime_digit_cases": [\n'
    cases = []
    for row in report['prime_digit_cases']:
        small = {k: v for k, v in row.items() if k != 'parent_cones'}
        block = json.dumps(small, indent=2)[:-2]+',\n  "parent_cones": [\n'
        block += ',\n'.join('    '+json.dumps(w, separators=(',', ':'))
                            for w in row['parent_cones'])
        block += '\n  ]\n}'
        cases.append('\n'.join('    '+line for line in block.splitlines()))
    text += ',\n'.join(cases)+'\n  ]\n}\n'
    assert json.loads(text) == json.loads(json.dumps(report))
    target.write_text(text)


def main():
    counts = defaultdict(int)
    finite = finite_controls(counts)
    completion = completion_controls(counts)
    digit = [digit_case(n, k, counts) for n, k in
             ((5, 3), (10, 3), (10, 5), (21, 5), (21, 13), (42, 13), (42, 21))]
    source_files = ['docs/projection_tag_fans.md', 'scripts/run_projection_tag_checks.py',
                    'scripts/run_signed_face_checks.py', 'scripts/run_wall_coexistence_checks.py',
                    'scripts/run_structural_complexity_checks.py',
                    'scripts/run_oracle_enumeration_checks.py']
    report = {'command': 'python3 -B scripts/run_projection_tag_checks.py',
              'counts': dict(sorted(counts.items())), 'finite_controls': finite,
              'completion_controls': completion, 'prime_digit_cases': digit,
              'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                for name in source_files},
              'uniform_stabilization_proved': True,
              'unrestricted_G3_improved': False, 'polynomial_lower_separation_proved': False,
              'asymptotic_inference_from_finite_cases': False, 'novelty_established': False,
              'goal_complete': False}
    target = ROOT/'results/research/projection_tag_checks.json'
    write_report(report, target)
    print(json.dumps({'counts': report['counts'],
                      'digit_counts': [{key: row[key] for key in
                         ('n', 'k', 'parent_g', 'tagged_length', 'tagged_g', 'gain')}
                                       for row in digit]}, indent=2))


if __name__ == '__main__':
    main()
