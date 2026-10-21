#!/usr/bin/env python3
"""Exact two-flag compilation, checked against independent full border fans."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
from random import Random
import hashlib
import json

from run_planar_dyadic_bound_checks import convex_hull
from run_projection_tag_checks import compare_fans, fan, point_model
from run_signed_face_checks import field_basis, field_reduce
from run_square_zero_extension_checks import extension, hom_space
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]
TRIANGLE = tuple(tuple(Q(i == j) for i in range(3)) for j in range(3))


def reduce_scalar(c, char):
    return c % char if char else Q(c)


def dense_rank(matrix, char):
    """Independent row elimination, retaining every matrix entry."""
    A = [[reduce_scalar(c, char) for c in row] for row in matrix]
    if not A:
        return 0
    k = 0
    for j in range(len(A[0])):
        pivot = next((i for i in range(k, len(A)) if A[i][j]), None)
        if pivot is None:
            continue
        A[k], A[pivot] = A[pivot], A[k]
        inv = pow(A[k][j], -1, char) if char else 1/A[k][j]
        A[k] = [reduce_scalar(c*inv, char) for c in A[k]]
        for i in range(k+1, len(A)):
            c = A[i][j]
            A[i] = [reduce_scalar(a-c*b, char) for a, b in zip(A[i], A[k])]
        k += 1
        if k == len(A):
            break
    return k


def rank_profile(columns, char):
    """Column-ordered elimination; pivots reveal every leading rank."""
    pivots, positions = {}, set()
    for j, column in enumerate(columns):
        v = [reduce_scalar(c, char) for c in column]
        for i, row in sorted(pivots.items()):
            c = v[i]
            if c:
                v = [reduce_scalar(a-c*b, char) for a, b in zip(v, row)]
        i = next((i for i, c in enumerate(v) if c), None)
        if i is not None:
            inv = pow(v[i], -1, char) if char else 1/v[i]
            pivots[i] = [reduce_scalar(c*inv, char) for c in v]
            positions.add((i, j))
    return positions


def profile_controls(counts):
    rng = Random(271828)
    for char, rows, cols, trial in product((0, 2, 3), range(1, 6), range(1, 6), range(3)):
        A = [[rng.randrange(3) for _ in range(cols)] for _ in range(rows)]
        positions = rank_profile(list(zip(*A)), char)
        assert len({i for i, j in positions}) == len(positions)
        assert len({j for i, j in positions}) == len(positions)
        for k, l in product(range(rows+1), range(cols+1)):
            actual = dense_rank([row[:l] for row in A[:k]], char)
            assert actual == sum(i < k and j < l for i, j in positions)
            counts['independent_leading_matrix_ranks'] += 1
        counts['rank_profile_matrices'] += 1


def joint_matching(A, B, n, ev, char, counts):
    C, pivots = field_basis(A+B, ev, char)
    assert len(C) == n and C[:len(A)] == A
    columns = []
    for b in B:
        rem, co = field_reduce(ev(b), pivots, char)
        assert not rem
        columns.append([co.get(c, 0) for c in reversed(C)])
    positions = rank_profile(columns, char)
    pairs = {(n-i-1, j) for i, j in positions if n-i-1 < len(A)}
    assert len(pairs) == len(A)+len(B)-n
    # A separate union-rank computation checks all prefixes, including zeros.
    ranks = {}
    for i, j in product(range(len(A)+1), range(len(B)+1)):
        labels = A[:i]+B[:j]
        direct = dense_rank([[ev(a).get(k, 0) for a in labels] for k in range(n)], char)
        via_profile = i+sum(r < n-i and c < j for r, c in positions)
        assert direct == via_profile
        ranks[i, j] = direct
        counts['independent_joint_prefix_ranks'] += 1
    direct_pairs = {(i-1, j-1) for i in range(1, len(A)+1) for j in range(1, len(B)+1)
                    if ranks[i, j-1]+ranks[i-1, j]-ranks[i-1, j-1]-ranks[i, j] == 1}
    assert direct_pairs == pairs
    counts['matching_factorizations'] += 1
    return pairs


def prefix_changes(old, new):
    assert len(old) == len(new)
    return sum(frozenset(old[:i]) != frozenset(new[:i]) for i in range(1, len(old)+1))


def upper_key(value):
    return (1, Q(0)) if value is None else (0, value)


def positive_carrier(a, b, lo, hi):
    slope, intercept = a[1]-b[1], a[2]-b[2]
    if slope == 0:
        if intercept <= 0:
            return None
    elif slope > 0:
        lo = max(lo, Q(-intercept, slope))
    else:
        root = Q(-intercept, slope)
        hi = root if hi is None else min(hi, root)
    if hi is not None and lo >= hi:
        return None
    return (1, -slope, -intercept), lo, hi


def merge_carriers(portions):
    out = []
    for lo, hi in sorted(portions, key=lambda v: (v[0], upper_key(v[1]))):
        if out and (out[-1][1] is None or lo <= out[-1][1]):
            if upper_key(hi) > upper_key(out[-1][1]):
                out[-1] = out[-1][0], hi
        else:
            out.append((lo, hi))
    return out


def check_wall_cover(actual, covers, counts):
    seen = set()
    for poly in actual.values():
        for a, b in zip(poly, poly[1:]+poly[:1]):
            if any(a[i] == b[i] == 0 for i in range(3)):
                continue
            cross = (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
            if not cross[0]:
                continue
            q = tuple(c/cross[0] for c in cross)
            assert all(c.denominator == 1 for c in q)
            q = tuple(int(c) for c in q)
            lo, hi = sorted((a[1]/a[2] if a[2] else None,
                             b[1]/b[2] if b[2] else None), key=upper_key)
            assert lo is not None
            edge = q, lo, hi
            if edge in seen:
                continue
            seen.add(edge)
            assert any(left <= lo and upper_key(right) >= upper_key(hi) for left, right in covers[q])
            counts['independent_actual_wall_piece_cover_controls'] += 1


def face_controls(actual, run_endpoints, counts):
    """Recover each positive two-face from the independent complete fan."""
    weights = {w for poly in actual.values() for w in poly if min(w) > 0}
    out = []
    for w in sorted(weights):
        # A profile minimizes this objective iff its cone contains w.
        profiles = {tuple(sum(a[i] for a in B) for i in range(3)) for B in actual}
        best = min(dot(w, v) for v in profiles)
        face = {v for v in profiles if dot(w, v) == best}
        hull = convex_hull([(v[0], v[1]) for v in face])
        assert len(hull) >= 3
        half = len(hull)//2
        sums = {tuple(a+b for a, b in zip(hull[i], hull[(i+half) % len(hull)]))
                for i in range(len(hull))}
        central = len(hull) % 2 == 0 and len(sums) == 1
        if not central:
            assert w in run_endpoints
            counts['noncentral_faces_charged_to_literal_carrier_endpoints'] += 1
        else:
            counts['centrally_symmetric_actual_two_faces'] += 1
        out.append({'weight': [str(c) for c in w], 'profiles': sorted(face),
                    'centrally_symmetric': central,
                    'literal_carrier_endpoint': w in run_endpoints})
    profiles = {tuple(sum(a[i] for a in B) for i in range(3)) for B in actual}
    boundary_normals = set()
    boundary_candidates = {w for poly in actual.values() for w in poly if min(w) == 0}
    for w in boundary_candidates:
        zeros = [i for i, c in enumerate(w) if c == 0]
        if len(zeros) == 2:
            boundary_normals.add(w)
        elif len(zeros) == 1:
            axes = [i for i, c in enumerate(w) if c]
            best = min(dot(w, v) for v in profiles)
            projected = {tuple(v[i] for i in axes) for v in profiles if dot(w, v) == best}
            if len(projected) > 1:
                boundary_normals.add(w)
    assert len(actual) <= 2*(len(out)+len(boundary_normals)+1)-4
    counts['independent_truncation_Euler_inequalities'] += 1
    return out, len(boundary_normals)


def compile_fan(n, ev, char, counts):
    planar = tuple(hyperbolic_cross(n, 2))
    ground = tuple(hyperbolic_cross(n, 3))
    cuts = [Q(0)]+sorted({Q(b, a) for a, b in product(range(1, n), repeat=2)})+[None]
    patches, carriers, literal_runs = defaultdict(list), defaultdict(list), defaultdict(list)
    records, previous = [], None
    births = changes = 0
    lower_histories, upper_histories = defaultdict(set), defaultdict(set)
    for lo, hi in zip(cuts, cuts[1:]):
        u = lo+1 if hi is None else (lo+hi)/2
        ordered = sorted(planar, key=lambda a: (u*a[0]+a[1], a))
        A = field_basis(tuple((0, *a) for a in ordered), ev, char)[0]
        B = field_basis(tuple((1, *a) for a in ordered), ev, char)[0]
        for k in range(1, len(A)+1):
            lower_histories[k].add(frozenset(A[:k]))
        for k in range(1, len(B)+1):
            upper_histories[k].add(frozenset(B[:k]))
        pairs = joint_matching(A, B, n, ev, char, counts)
        labels = {(A[i], B[j]) for i, j in pairs}
        for a, b in labels:
            portion = positive_carrier(a, b, lo, hi)
            if portion is not None:
                direction, left, right = portion
                carriers[direction].append((left, right))
                literal_runs[a, b].append((left, right))
        if previous is None:
            births = len(labels)
        else:
            old_A, old_B, old_labels = previous
            changed = prefix_changes(old_A, A)+prefix_changes(old_B, B)
            newly_born = len(labels-old_labels)
            # A changed prefix can affect only its two neighboring label rows.
            assert newly_born <= 4*changed
            births += newly_born
            changes += changed
            counts['matching_birth_transition_controls'] += 1
        previous = A, B, labels
        poly = list(TRIANGLE)
        poly = clip(poly, (0, lo.denominator, -lo.numerator), 1)
        if hi is not None:
            poly = clip(poly, (0, -hi.denominator, hi.numerator), 1)
        assert poly
        cells = [poly]
        for a, b in sorted(labels):
            normal = tuple(x-y for x, y in zip(a, b))
            cells = [piece for cell in cells for sign in (-1, 1)
                     if (piece := clip(cell, normal, sign))]
        assert len(cells) <= 1+len(pairs)*(len(pairs)+1)//2
        matched = {a for pair in labels for a in pair}
        fixed = (set(A)|set(B))-matched
        for cell in cells:
            w = tuple(sum(v[i] for v in cell)/len(cell) for i in range(3))
            selected = frozenset(fixed | {a if dot(w, a) < dot(w, b) else b for a, b in labels})
            direct = frozenset(field_basis(sorted(ground, key=lambda a: (dot(w, a), a)), ev, char)[0])
            assert selected == direct and len(selected) == n
            assert all(a[0] < 2 for a in selected)
            patches[selected].extend(cell)
            counts['independent_compiled_cell_greedy_equalities'] += 1
        # Two increasing endpoint orders are necessary for an interior crossing.
        for i, j in pairs:
            for k, l in pairs:
                if i < k and j > l:
                    assert u*(A[i][1]-B[j][1])+A[i][2]-B[j][2] < (
                        u*(A[k][1]-B[l][1])+A[k][2]-B[l][2])
                    counts['oppositely_ordered_pair_noncrossing_controls'] += 1
        records.append({'u_lower': str(lo), 'u_upper': str(hi) if hi is not None else 'infinity',
                        'lower_basis': A, 'upper_basis': B, 'matched_label_pairs': sorted(labels),
                        'patch_cells': len(cells)})
    out = {}
    for basis, vertices in patches.items():
        hull = convex_hull([(v[0], v[1]) for v in vertices])
        out[basis] = [(a, b, 1-a-b) for a, b in hull]
    assert sum(area2(poly) for poly in out.values()) == 1
    assert changes == sum(len(v)-1 for v in lower_histories.values())+sum(
        len(v)-1 for v in upper_histories.values())
    assert births <= n+4*changes
    covers = {q: merge_carriers(v) for q, v in carriers.items()}
    cover_size = sum(len(v) for v in covers.values())
    assert cover_size <= births
    run_endpoints = set()
    for (a, b), portions in literal_runs.items():
        for lo, hi in merge_carriers(portions):
            for u in (lo, hi):
                if u is not None and u > 0:
                    t = (a[1]-b[1])*u+a[2]-b[2]
                    if t > 0:
                        run_endpoints.add(tuple(c/(t+u+1) for c in (t, u, Q(1))))
    assert len(run_endpoints) <= 2*births
    counts['compiled_complete_fans'] += 1
    return out, {'matching_carrier_occurrences': births, 'changed_prefix_sets': changes,
                 'nonvertical_support_cover_size': cover_size,
                 'nonvertical_support_cover': [{'direction': q, 'u_intervals': [
                     [str(lo), str(hi) if hi is not None else 'infinity'] for lo, hi in spans]}
                      for q, spans in sorted(covers.items())],
                 'positive_literal_carrier_endpoints': [list(map(str, w)) for w in sorted(run_endpoints)],
                 'universal_u_strips': len(records), 'strips': records}, covers, run_endpoints


def square_zero_fixtures():
    cases = [(set(product(range(2), repeat=2)), {(0, 0)}),
             ({(0, 0), (1, 0), (2, 0), (0, 1)}, {(0, 0), (1, 0)}),
             (set(product(range(2), repeat=2)), set(product(range(2), repeat=2))),
             ({(0, 0), (1, 0), (0, 1)}, set())]
    for case, (D, E) in enumerate(cases):
        for char in (2, 3, 101):
            G, variables, hom = hom_space(D, E, char)
            rng = Random(1234+case*char)
            parameters = [rng.randrange(char) for _ in hom]
            tails = [{} for _ in G]
            for j, (i, v) in enumerate(variables):
                c = sum(a*b[j] for a, b in zip(parameters, hom)) % char
                if c:
                    tails[i][v] = c
            n, ev, _ = extension(D, E, G, tails, char)
            yield {'model': 'square_zero_extension', 'case': case, 'characteristic': char,
                   'relation': [0, 0, 1]}, n, ev


def quadratic_tensor(char, r, constant, linear):
    """F[x,z]/(x^2+linear*x+constant,z^r), y=x+z."""
    @lru_cache(None)
    def xpower(k):
        if k < 2:
            return (int(k == 0), int(k == 1))
        a, b = xpower(k-1)
        return reduce_scalar(-constant*b, char), reduce_scalar(a-linear*b, char)
    from math import comb
    @lru_cache(None)
    def ev(a):
        x, y, z = a
        out = defaultdict(int)
        for j in range(y+1):
            if z+j < r:
                for i, c in enumerate(xpower(x+y-j)):
                    out[i*r+z+j] += comb(y, j)*c
        return {i: reduce_scalar(c, char) for i, c in out.items() if reduce_scalar(c, char)}
    return 2*r, ev


def socle_crossing():
    # On the two-layer downset identify xy=z^2 and xz=y^2.
    # The two actual walls cross at t=u=1, a planar order event.
    basis = ((0, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 0), (1, 1, 0), (1, 0, 1))
    columns = {a: {i: 1} for i, a in enumerate(basis)}
    columns[(0, 2, 0)] = {5: 1}
    columns[(0, 0, 2)] = {4: 1}
    return 6, lambda a: dict(columns.get(a, {}))


def main():
    counts, fixtures = defaultdict(int), []
    profile_controls(counts)
    models = list(square_zero_fixtures())
    for char in (0, 2):
        n, ev = socle_crossing()
        models.append(({'model': 'socle_crossing_at_planar_event', 'characteristic': char,
                        'relation': [0, 0, 1]}, n, ev))
    for char, r, c, b in ((2, 3, 1, 1), (2, 3, 1, 0), (3, 2, 1, 0),
                          (0, 3, 1, 0), (0, 2, -1, 0)):
        n, ev = quadratic_tensor(char, r, c, b)
        models.append(({'model': 'quadratic_tensor', 'characteristic': char,
                        'relation': [c, b, 1], 'z_nilpotence': r}, n, ev))
    for char in (0, 101):
        points = [(0, 0, 0), (0, 1, 2), (1, 0, 1), (1, 2, 3), (1, 3, 2)]
        n, ev = point_model(points, char)
        models.append(({'model': 'distinct_two_plane_points', 'characteristic': char,
                        'relation': [0, -1, 1], 'points': points}, n, ev))
    for x in (0, 1):
        points = [(x, 0, 0)]
        n, ev = point_model(points, 2)
        models.append(({'model': 'one_point_boundary_control', 'characteristic': 2,
                        'relation': [0, -1, 1], 'points': points}, n, ev))
    for meta, n, ev in models:
        char = meta['characteristic']
        c, b, _ = meta['relation']
        for y, z in hyperbolic_cross(n, 2):
            relation = defaultdict(int)
            for x, coefficient in ((0, c), (1, b), (2, 1)):
                for i, value in ev((x, y, z)).items():
                    relation[i] += coefficient*value
            assert all(not reduce_scalar(v, char) for v in relation.values())
            counts['quadratic_relation_column_controls'] += 1
        compiled, record, covers, endpoints = compile_fan(n, ev, char, counts)
        actual = fan(n, ev, char, counts)
        compare_fans(actual, compiled, counts)
        check_wall_cover(actual, covers, counts)
        faces, boundary_facets = face_controls(actual, endpoints, counts)
        record['positive_exposed_two_faces'] = faces
        record['noncompact_state_facets'] = boundary_facets
        fixtures.append(dict(meta, length=n, initial_ideals=len(actual), **record))
    sources = ['docs/quadratic_coordinate_fan_compilation.md',
               'scripts/run_quadratic_coordinate_fan_checks.py',
               'scripts/run_projection_tag_checks.py', 'scripts/run_signed_face_checks.py',
               'scripts/run_square_zero_extension_checks.py',
               'scripts/run_structural_complexity_checks.py',
               'scripts/run_planar_dyadic_bound_checks.py']
    report = {'command': 'python3 -B scripts/run_quadratic_coordinate_fan_checks.py',
              'counts': dict(sorted(counts.items())), 'complete_fan_fixtures': fixtures,
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              'unrestricted_G3_improved': False, 'new_quadratic_class_count_bound': False,
              'collective_cubic_budget_proved': False, 'novelty_established': False,
              'optimized_arithmetic_algorithm_implemented': False, 'goal_complete': False}
    (ROOT/'results/research/quadratic_coordinate_fan_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'counts': report['counts'],
                      'fixtures': [{k: row[k] for k in ('model', 'characteristic', 'length', 'initial_ideals',
                                                       'matching_carrier_occurrences', 'changed_prefix_sets')}
                                   for row in fixtures]}, indent=2))


if __name__ == '__main__':
    main()
