#!/usr/bin/env python3
"""Power witnesses, exact finite fans and the bounded-coordinate carrier gate."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json

from run_extremal_initial_ideals_checks import adjacent_weights
from run_higher_interaction_checks import ground, linear_sum, rank_trace
from run_hypersurface_truncation_checks import poly_mul, relation_rows
from run_interleaving_research import minimal_excluded_border
from run_signed_face_checks import normals
from run_socle_stability_checks import depth_rref
from run_structural_complexity_checks import area2, clip, determinant, dot, primitive
from run_two_interaction_checks import model as two_model

ROOT = Path(__file__).resolve().parents[1]


def add_polys(*terms):
    out = defaultdict(int)
    for coefficient, f in terms:
        for v, c in f.items():
            out[v] += coefficient*c
    return {v: c for v, c in out.items() if c}


def power(f, e):
    out = {(0, 0, 0): 1}
    for _ in range(e):
        out = poly_mul(out, f, 0)
    return out


def local_model(a, s, q, h):
    B = tuple(v for v in ground(s) if v[2] < a)
    index = {v: i for i, v in enumerate(B)}
    F = add_polys((1, q), (1, h))
    assert F[0, 0, a] == 1

    @lru_cache(None)
    def ev(v):
        if sum(v) > s:
            return {}
        if v[2] < a:
            return {index[v]: 1}
        return linear_sum([(-c, ev((v[0]+i, v[1]+j, v[2]+k-a)))
                           for (i, j, k), c in F.items() if (i, j, k) != (0, 0, a)])
    return B, ev, F


def power_controls():
    counts, records = defaultdict(int), []
    for a, s, D in ((2, 7, 5), (2, 10, 4), (3, 10, 6), (2, 6, 3)):
        q = {(a, 0, 0): 1, (0, a, 0): 1, (0, 0, a): 1, (a-1, 1, 0): -2}
        h = {(D-2, 1, 1): 1, (2, 0, D-2): -3, (0, D, 0): 1, (0, 0, D+1): 2}
        e = (s+D)//D
        B, ev, F = local_model(a, s, q, h)
        n = len(B)
        hilbert = [(j+2)*(j+1)//2-((j-a+2)*(j-a+1)//2 if j >= a else 0)
                   for j in range(s+1)]
        assert n == sum(hilbert)
        S = sum(j*c for j, c in enumerate(hilbert))
        negative_h = {v: -c for v, c in h.items()}
        multiplier = add_polys(*[(1, poly_mul(power(q, e-1-j), power(negative_h, j), 0))
                                 for j in range(e)])
        qe = power(q, e)
        replay = poly_mul(F, multiplier, 0)
        difference = add_polys((1, qe), (-1, replay))
        assert difference == power(negative_h, e)
        assert all(sum(v) >= s+1 for v in difference)
        assert all(qe[tuple(a*e*int(i == j) for i in range(3))] == 1 for j in range(3))
        counts['integer_power_membership_witnesses'] += 1
        E = ground(s)
        for p in (0, 2, 3, 5):
            rows = relation_rows(E, [F], p)
            assert len(E)-len(rows) == n
            for row in rows:
                image = linear_sum([(c, ev(v)) for v, c in zip(E, row) if c])
                assert all(not (c % p if p else c) for c in image.values())
                counts['independent_hypersurface_kernel_rows'] += 1
            image = linear_sum([(c, ev(v)) for v, c in qe.items()])
            assert all(not (c % p if p else c) for c in image.values())
            for w in permutations((1, s+2, (s+2)**2)):
                selected, _ = rank_trace(sorted(E, key=lambda v: (dot(v, w), v)), ev, p, set())
                largest = max(range(3), key=lambda i: w[i])
                assert len(selected) == n
                assert all(v[largest] < a*e for v in selected)
                assert sum(map(sum, selected)) <= S
                counts['independent_order_and_moment_controls'] += 1
        records.append({'tangent_degree': a, 'cutoff': s, 'first_perturbation_degree': D,
                        'power': e, 'carrier_degree': a*e, 'length': n,
                        'local_hilbert_function': hilbert, 'local_moment': S,
                        'q': [[v, c] for v, c in sorted(q.items())],
                        'h': [[v, c] for v, c in sorted(h.items())],
                        'exact_multiplier': [[v, c] for v, c in sorted(multiplier.items())]})
    return counts, records


def complete_fan(r, p, counts):
    s = 2*r+1
    B0, ev, _, F = two_model(r, {(2, 0): 1, (0, 2): 1}, {(1, 1): 1}, {(2, 0): 1})
    E = ground(s)
    comparisons = sorted({primitive(tuple(x-y for x, y in zip(v, w)))
                          for v, w in combinations(E, 2)
                          if min(x-y for x, y in zip(v, w)) < 0 < max(x-y for x, y in zip(v, w))})
    pending, seen = [(1, s+2, (s+2)**2)], {}
    rows = relation_rows(E, [F], p) if r == 3 else None
    S = s*(s+1)*(4*s+5)//6
    while pending:
        w = pending.pop()
        B, _ = rank_trace(sorted(E, key=lambda a: (dot(a, w), a)), ev, p, set())
        key = frozenset(B)
        if key in seen:
            continue
        Z, _ = normals(B, ev, p)
        poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for z in Z:
            poly = clip(poly, z, 1)
        assert area2(poly) > 0 and all(dot(z, w) > 0 for z in Z)
        largest = max(range(3), key=lambda i: w[i])
        moment = tuple(sum(v[i] for v in B) for i in range(3))
        assert len(B) == len(B0) and max(v[largest] for v in B) < 4
        assert moment[largest] <= 3*len(B0) and sum(moment) <= S
        if rows is not None:
            order = sorted(range(len(E)), key=lambda i: (dot(E[i], w), E[i]), reverse=True)
            _, pivots = depth_rref(rows, p, order)
            independent = {v for i, v in enumerate(E) if i not in pivots}
            assert independent == B
            counts['independent_macaulay_initial_sets'] += 1
        record = {'minimal_initial_generators': sorted(minimal_excluded_border(tuple(B))),
                  'strict_weight': list(map(str, w)), 'largest_variable': largest,
                  'profile': moment, 'polygon': [[str(x) for x in v] for v in poly]}
        seen[key] = record
        pending.extend(adjacent_weights(poly, Z, comparisons))
        counts['complete_fan_cones_and_slab_caps'] += 1
    assert sum(area2([tuple(map(Q, v)) for v in record['polygon']]) for record in seen.values()) == 1
    counts['complete_positive_fans'] += 1
    return {'r': r, 'cutoff': s, 'length': len(B0), 'characteristic': p,
            'initial_ideals': len(seen), 'normalized_total_area': '1',
            'cones': list(seen.values())}, seen


def missing_wall_control(seen, counts):
    r, s = 3, 7
    q = {(2, 0, 0): 1, (0, 2, 0): 1, (0, 0, 2): 1}
    B0, ev, _, F = two_model(r, {(2, 0): 1, (0, 2): 1}, {(1, 1): 1}, {(2, 0): 1})
    multiplier = {(1, 0, 0): 1, (0, 1, 1): 1, (1, 0, 4): -1}
    G = add_polys((1, poly_mul(q, {(1, 0, 0): 1, (0, 1, 1): 1}, 0)),
                  (1, {(2, 1, 3): 1}), (-1, {(1, 0, 6): 1}))
    replay = poly_mul(F, multiplier, 0)
    assert G == {v: c for v, c in replay.items() if sum(v) <= s}
    weights = [(Q(5, 4), Q(7, 4), Q(1)+Q(sign, 1000)) for sign in (-1, 1)]
    bases = []
    for w in weights:
        B, _ = rank_trace(sorted(ground(s), key=lambda v: (dot(v, w), v)), ev, 0, set())
        assert frozenset(B) in seen
        Z, _ = normals(B, ev, 0)
        assert all(dot(z, w) > 0 for z in Z)
        bases.append(B)
    assert bases[0]-bases[1] == {(1, 0, 6)}
    assert bases[1]-bases[0] == {(2, 1, 3)}
    counts['actual_missing_entry_wall_adjacent_cones'] += 2
    entries = set()
    for pivot in range(3):
        basis = tuple(v for v in ground(s) if v[pivot] < 2)
        index = {v: i for i, v in enumerate(basis)}
        head = tuple(2*int(i == pivot) for i in range(3))

        @lru_cache(None)
        def coordinate_ev(v):
            if sum(v) > s:
                return {}
            if v[pivot] < 2:
                return {index[v]: 1}
            return linear_sum([(-c, coordinate_ev(tuple(x+y-z for x, y, z in zip(v, a, head))))
                               for a, c in F.items() if a != head])

        for v in ground(s):
            assert ev(v) == linear_sum([(c, ev(basis[i])) for i, c in coordinate_ev(v).items()])
            counts['independent_coordinate_normal_form_equalities'] += 1
            for i in coordinate_ev(v):
                z = tuple(x-y for x, y in zip(v, basis[i]))
                if min(z) < 0 < max(z):
                    entries.add(primitive(z))
    assert (1, 1, -3) not in entries
    return {'wall': [1, 1, -3], 'prices_on_wall': ['5/4', '7/4', '1'],
            'relation': [[v, c] for v, c in sorted(G.items())],
            'ideal_multiplier': [[v, c] for v, c in sorted(multiplier.items())],
            'adjacent_prices': [list(map(str, w)) for w in weights],
            'lost_standard': [1, 0, 6], 'gained_standard': [2, 1, 3],
            'all_coordinate_entry_wall_directions': sorted(entries)}


def coordinate_controls(counts):
    # For the chosen polynomial xy, its nonvanishing open condition has no F2 points.
    # This does not claim that a particular ideal needs an extension for every carrier.
    good = []
    for columns in product(tuple(product((0, 1), repeat=3)), repeat=3):
        if int(determinant(columns)) % 2 and all(v[0]*v[1] % 2 for v in columns):
            good.append(columns)
        counts['finite_field_coordinate_matrices_checked'] += 1
    assert not good
    columns = ((1, 1, 0), (1, 2, 0), (1, 1, 1))
    assert determinant(columns) == 1
    assert all(v[0]*v[1] for v in columns)
    assert all(v[0]*v[1] % 3 for v in columns)

    def f4mul(a, b):
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a & 4:
                a ^= 7
        return out
    columns4 = ((1, 1, 0), (1, 2, 0), (1, 1, 1))
    det4 = 0
    for perm in permutations(range(3)):
        term = 1
        for i in range(3):
            term = f4mul(term, columns4[i][perm[i]])
        det4 ^= term
    assert det4 and all(f4mul(v[0], v[1]) for v in columns4)
    counts['coordinate_extension_positive_controls'] += 2
    return {'polynomial': 'x*y', 'good_GL3_F2_matrices_for_this_polynomial': 0,
            'rational_and_F3_columns': columns, 'F4_columns': columns4,
            'F4_modulus': 'alpha^2+alpha+1', 'F4_determinant_encoding': det4,
            'scope': 'Chosen-polynomial open condition only; not an ideal-level necessity claim.'}


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    counts, powers = power_controls()
    fans, first = [], None
    for r, p, expected in ((3, 0, 27), (3, 2, 27), (5, 101, 65)):
        record, seen = complete_fan(r, p, counts)
        assert record['initial_ideals'] == expected
        fans.append(record)
        if (r, p) == (3, 0):
            first = seen
    wall = missing_wall_control(first, counts)
    coordinates = coordinate_controls(counts)
    sources = ('docs/power_carrier_bounds.md', 'scripts/run_power_carrier_checks.py',
               'docs/two_interaction_ranks.md', 'scripts/run_two_interaction_checks.py',
               'scripts/run_extremal_initial_ideals_checks.py', 'scripts/run_higher_interaction_checks.py',
               'scripts/run_hypersurface_truncation_checks.py', 'scripts/run_signed_face_checks.py',
               'scripts/run_socle_stability_checks.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_higher_seed_checks.py')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_power_carrier_checks.py',
              'counts': dict(counts), 'power_membership_fixtures': powers,
              'complete_fan_fixtures': fans, 'missing_entry_wall': wall,
              'coordinate_open_condition_control': coordinates,
              'uniform_results': {'pure_carrier_bound': 'min((S+1)*sqrt(b+1),(b+1)*(S+1)^(2/3)); b=min(S,(m-1)*n).',
                                  'truncated_hypersurface_bound': 'O_a(sqrt(e)*n^2), e=ceil((s+1)/D).',
                                  'fixed_positive_fraction_perturbation_degree': 'O_(a,c)(n^2), for every number of higher terms and every field.',
                                  'two_interaction_first_order_family': 'Excluded as a source of a polynomial gain; q^2 lies in the ideal.',
                                  'generic_linear_coordinate_corollary': 'O(n^(8/3)) after a suitable generic linear change, over an extension if necessary.',
                                  'coordinate_entry_wall_refinement': False},
              'unrestricted_G3_improvement': False, 'novelty_established': False,
              'goal_complete': False,
              'scope': 'Classical lattice geometry applied to an ideal power witness; uniform family exclusion and a coordinate-dependent corollary. Finite complete fans verify examples, not asymptotic sharpness or novelty.',
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/'results/research/power_carrier_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
