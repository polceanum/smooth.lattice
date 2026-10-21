#!/usr/bin/env python3
"""Exact supporting controls for the individual product-shadow obstruction.

No fitted asymptotic or novelty inference. Uses only the standard library.
"""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
from math import gcd, prod
from pathlib import Path
from random import Random
import hashlib
import json

from run_structural_complexity_checks import area2, clip, cross, dot

ROOT = Path(__file__).resolve().parents[1]
L = [[2, -1, -1], [-1, 2, -1], [-1, -1, 2]]
B = [[2, -1, 0], [0, 1, -1], [-3, 2, 0]]
C = [[54, -45, 0], [6, 15, -21], [-69, 48, 12]]
E = [[34, -32, 10], [25, -11, -8], [-53, 49, -14]]
F = [[5, -5, 3], [9, -6, 1], [-6, 7, -5]]
IDENTITY = [[int(i == j) for j in range(3)] for i in range(3)]


def multiply(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(3))
             for j in range(3)] for i in range(3)]


def row_times(p, a):
    return tuple(sum(p[i]*a[i][j] for i in range(3)) for j in range(3))


def mixed(p):
    return min(p) < 0 < max(p)


def norm(p):
    return max(map(abs, p))


def rho(p):
    return prod(abs(x)+1 for x in p)


def determinant(a):
    return dot(a[0], cross(a[1], a[2]))


def adjugate(a):
    return list(map(list, zip(cross(a[1], a[2]), cross(a[2], a[0]), cross(a[0], a[1]))))


def add_scaled(*terms):
    return [[sum(c*a[i][j] for c, a in terms) for j in range(3)] for i in range(3)]


def polynomial_multiply(a, b):
    out = [[[0]*3 for _ in range(3)] for _ in range(len(a)+len(b)-1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] = add_scaled((1, out[i+j]), (1, multiply(x, y)))
    return out


def uniform_identities(counts):
    # Coefficientwise identities over Z[t], independent of evaluating at t.
    inverse = [[[2, 0, 1], [3, 0, 2], [3, -1, 2]],
               [[3, 0, 3], [4, 1, 4], [5, 2, 5]], [[3]*3 for _ in range(3)]]
    identity = polynomial_multiply([B, L], inverse)
    assert identity == [IDENTITY]+[[[0]*3 for _ in range(3)] for _ in range(3)]
    cube = polynomial_multiply(polynomial_multiply([B, L], [B, L]), [B, L])
    assert cube == [F, E, C, add_scaled((9, L))]
    assert determinant(B) == 1  # A polynomial inverse plus this fixes det M = 1.
    counts['coefficientwise_integer_matrix_identities'] += 3
    cross_columns = [cross([L[i][j] for i in range(3)], [C[i][j] for i in range(3)])
                     for j in range(3)]
    assert cross_columns == [(75, 84, 66), (111, 93, 75), (30, 12, 21)]
    assert all(min(p) > 0 for p in cross_columns)
    assert [max(sum(abs(a[i][j]) for i in range(3)) for j in range(3))
            for a in (B, C, E, F)] == [5, 129, 112, 20]
    counts['exact_column_separation_controls'] += 3
    return {'M_coefficients': [B, L], 'inverse_coefficients': inverse,
            'cube_coefficients': cube, 'positive_common_orthogonals': cross_columns}


def endpoints(a):
    assert determinant(a) > 0
    assert all(x > 0 for row in adjugate(a) for x in row)
    assert all(a[i][j] < 0 for i in range(3) for j in range(3) if i != j)
    assert all(a[j][j] > -a[i][j] for i in range(3) for j in range(3) if i != j)
    pure = [tuple(a[i][i]*int(i == j) for j in range(3)) for i in range(3)]
    plane = [tuple(-a[i][j] if i != j else 0 for j in range(3)) for i in range(3)]
    assert all(not all(x <= y for x, y in zip(u, v))
               for u, v in product(pure+plane, repeat=2) if u != v)
    n = sum(a[i][i] for i in range(3))+sum(prod(-a[i][j] for j in range(3) if i != j)
                                         for i in range(3))-2
    return pure, plane, n


def full_fan(a, counts):
    pure, plane, n = endpoints(a)
    records = []
    for signs in product((-1, 1), repeat=3):
        polygon = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for row, sign in zip(a, signs):
            polygon = clip(polygon, row, sign)
        if not polygon:
            assert signs == (-1, -1, -1)
            continue
        w = tuple(sum(v[i] for v in polygon)/len(polygon) for i in range(3))
        assert min(w) > 0 and all(sign*dot(row, w) > 0 for sign, row in zip(signs, a))
        heads = [pure[i] if signs[i] == 1 else plane[i] for i in range(3)]
        tails = [plane[i] if signs[i] == 1 else pure[i] for i in range(3)]
        assert all(dot(h, w) > dot(b, w) for h, b in zip(heads, tails))
        records.append({'comparison_signs': signs, 'removed_heads': heads, 'retained_tails': tails,
                        'strict_weight': list(map(str, w)),
                        'polygon': [[str(x) for x in v] for v in polygon]})
    assert len(records) == 7
    assert sum(area2([tuple(map(Q, v)) for v in r['polygon']]) for r in records) == 1
    counts['complete_exact_positive_fans'] += 1
    counts['initial_cones'] += len(records)
    return {'matrix': a, 'length': n, 'axis_endpoints': pure, 'plane_endpoints': plane,
            'initial_ideals': len(records), 'normalized_area': '1', 'cones': records}


def quotient_control(a, counts):
    fixture = full_fan(a, counts)
    pure, plane, n = endpoints(a)
    D = set()
    for endpoint in pure+plane:
        D.update(product(*(range(v+1) for v in endpoint)))
    assert len(D) == n+3
    for v in pure+plane:
        assert all(tuple(x+int(i == j) for i, x in enumerate(v)) not in D for j in range(3))
    counts['explicit_divisor_union_monomials'] += len(D)
    for record in fixture['cones']:
        heads, tails = map(tuple, record['removed_heads']), map(tuple, record['retained_tails'])
        replacement = dict(zip(heads, tails))
        basis = D-set(replacement)
        assert len(basis) == n

        def ev(v):
            return replacement.get(v, v) if v in D else None

        actions = {v: tuple(ev(tuple(x+int(i == j) for i, x in enumerate(v)))
                            for j in range(3)) for v in basis}

        def apply(v, j):
            return actions[v][j] if v is not None else None

        # Each defining socle relation is killed, and the cyclic monomials
        # have the proposed independent integral basis normal forms.
        assert ev((0, 0, 0)) == (0, 0, 0)
        for h, b in replacement.items():
            assert ev(h) == ev(b) and all(apply(ev(h), j) is None for j in range(3))
            counts['integral_socle_relations'] += 1
        for v in D:
            for j in range(3):
                inc = tuple(x+int(i == j) for i, x in enumerate(v))
                assert apply(ev(v), j) == ev(inc)
                counts['independent_cyclic_monomial_transitions'] += 1
        for v in basis:
            for i, j in combinations(range(3), 2):
                assert apply(apply(v, i), j) == apply(apply(v, j), i)
                counts['sparse_commuting_multiplication_checks'] += 1
            assert all(tuple(x-int(i == j) for i, x in enumerate(v)) in basis
                       for j in range(3) if v[j])
        # All actions are unit columns or zero. Integral verification
        # therefore also verifies reduction in every characteristic.
    fixture['explicit_integral_quotient_verification'] = True
    return fixture


def visibility_controls(a, counts, radius):
    inverse = adjugate(a)
    rays = list(zip(*inverse))
    tested, visible, boundary_only = 0, 0, 0
    for q in product(range(-radius, radius+1), repeat=3):
        if not any(q) or gcd(*q) != 1:
            continue
        p = row_times(q, inverse)
        # Independent geometric test against the three exact simplex rays.
        levels = [Q(dot(q, r), sum(r)) for r in rays]
        assert mixed(p) == (min(levels) < 0 < max(levels))
        if min(levels) == 0 or max(levels) == 0:
            boundary_only += int(not mixed(p))
        visible += int(mixed(p))
        tested += 1
    counts['exhaustive_primitive_plane_controls'] += tested
    counts['strict_visibility_controls'] += visible
    counts['boundary_only_planes_excluded'] += boundary_only
    return {'coordinate_radius': radius, 'primitive_vectors': tested,
            'visible_vectors': visible, 'boundary_only_excluded': boundary_only,
            'scope': 'Complete coordinate cube, not complete product minimum.'}


def family_controls(counts):
    rng = Random(172903)
    records = []
    for t in (10, 1000, 1001, 10000, 10**12, 10**30):
        m = add_scaled((t, L), (1, B))
        a = multiply(multiply(m, m), m)
        assert determinant(m) == determinant(a) == 1
        inverse = adjugate(a)
        assert multiply(a, inverse) == IDENTITY
        assert a == add_scaled((9*t**3, L), (t*t, C), (t, E), (1, F))
        assert [sum(row) for row in a] == [9*t*t+12*t+3, 6*t+4, -9*t*t-18*t-4]
        fixture = full_fan(a, counts)
        q = row_times((1, -1, 0), a)
        assert q == (27*t**3+48*t*t+9*t-4, -27*t**3-60*t*t-21*t+1, 21*t*t+18*t+2)
        assert gcd(*q) == 1 and row_times(q, inverse) == (1, -1, 0)
        # A strict interior point on q^perp, rather than a wall-only witness.
        w = tuple(sum(row) for row in inverse)
        assert min(w) > 0 and dot(q, w) == 0
        assert tuple(dot(row, w) for row in a) == (1, 1, 1)
        if t >= 1000:
            assert 192*t**6 <= fixture['length'] <= 301*t**6
            assert rho(q) <= 18000*t**8
            ps = [p for p in product(range(-3, 4), repeat=3) if mixed(p)]
            for scale in (t//1000, t, t*t):
                for _ in range(60):
                    p = tuple(rng.randint(-scale-3, scale+3) for _ in range(3))
                    if mixed(p):
                        ps.append(p)
            for p in ps:
                v = p
                for _ in range(3):
                    u = row_times(v, m)
                    assert mixed(u) and norm(u) >= (t-5)*norm(v)
                    v = u
                assert v == row_times(p, a)
                assert abs(sum(v))*2 < norm(v)
                assert rho(v)*256000000 >= t**8
                if norm(p)*1000 <= t:
                    assert min(map(abs, v))*2 >= t*t
                    counts['small_integral_row_coordinate_separations'] += 1
                else:
                    assert norm(v)*8000 > t**4
                    counts['large_integral_row_growth_controls'] += 1
                counts['mixed_row_expansion_and_product_controls'] += 1
        fixture['primitive_visible_plane'] = {'q': q, 'product': rho(q),
                                              'strict_interior_weight': w}
        if t in (10, 1000):
            fixture['exhaustive_visibility'] = visibility_controls(a, counts, 8)
        records.append(fixture)
    return records


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    counts = defaultdict(int)
    identities = uniform_identities(counts)
    small = quotient_control([[2, -1, -1], [-1, 3, -1], [-1, -1, 4]], counts)
    small['exhaustive_visibility'] = visibility_controls(small['matrix'], counts, 8)
    m5 = add_scaled((5, L), (1, B))
    quotient = quotient_control(multiply(m5, m5), counts)
    family = family_controls(counts)
    sources = ('docs/multiplicative_shadow_obstruction.md',
               'scripts/run_multiplicative_shadow_checks.py',
               'scripts/run_structural_complexity_checks.py',
               'docs/rational_shadow_bound.md', 'docs/initial_cone_stability.md')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_multiplicative_shadow_checks.py',
              'counts': dict(counts), 'polynomial_identities': identities,
              'explicit_quotient_fixtures': [small, quotient], 'cubic_family': family,
              'uniform_result': {'fields': 'all', 'family_length': 'Theta(t^6)',
                                 'initial_ideals': 7, 'minimum_plane_product': 'Theta(t^8)=Theta(n^(4/3))',
                                 'lower_constant': '1/256000000', 'upper_constant': 18000,
                                 'parameter_range': 'integer t>=1000'},
              'unrestricted_G3_improvement': False, 'novelty_established': False, 'goal_complete': False,
              'scope': 'Uniform theorem rests on the proof. Exact coefficient identities, full seven-cone fans, sparse integral commuting quotients, visibility cubes and mixed-row controls support it. No exhaustive minimum-product claim, large quotient allocation, asymptotic fit or publication conclusion.',
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/'results/research/multiplicative_shadow_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(counts), indent=2))


if __name__ == '__main__':
    main()
