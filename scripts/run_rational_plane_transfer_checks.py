#!/usr/bin/env python3
"""Exact supporting checks for rational-plane transfer; no growth fitting."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import gcd, isqrt
from pathlib import Path
from random import Random
import hashlib
import json

from run_planar_dyadic_bound_checks import convex_hull, matrix_oracle
from run_rational_shadow_bound_checks import projection_data, sweep
from run_structural_complexity_checks import dot, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]


def canonical(v):
    g = gcd(*v)
    if not g:
        return None
    s = 1 if next(x for x in v if x) > 0 else -1
    return tuple(s*x//g for x in v)


def wedge(u, v):
    return tuple(u[i]*v[j]-u[j]*v[i] for i, j in combinations(range(len(u)), 2))


def height2(p):
    return sum(x*x for x in p)


def contains(p, w):
    d = len(w)
    z = dict(zip(combinations(range(d), 2), p))
    return all(w[i]*z[j, k]-w[j]*z[i, k]+w[k]*z[i, j] == 0
               for i, j, k in combinations(range(d), 3))


def plane_enumeration(d, T, counts):
    # Gauss reduction bounds ALL possible reduced bases at this height.
    bound = isqrt(2*T)+1
    us = [u for u in product(range(-bound, bound+1), repeat=d)
          if any(u) and canonical(u) == u and 3*height2(u)**2 <= 4*T*T]
    planes, pairs = {}, 0
    for u in us:
        U2 = height2(u)
        V2 = 4*T*T//(3*U2)
        radius = isqrt(V2)
        for v in product(range(-radius, radius+1), repeat=d):
            if not U2 <= height2(v) <= V2 or 2*abs(dot(u, v)) > U2:
                continue
            raw = wedge(u, v)
            if gcd(*raw) != 1 or not 0 < height2(raw) <= T*T:
                continue
            p = canonical(raw)
            assert 4*height2(raw) >= 3*U2*height2(v)
            assert contains(p, u) and contains(p, v)
            planes.setdefault(p, (u, v))
            pairs += 1
    if d == 3:
        # Completely separate enumeration by primitive orthogonal normals.
        normals = {canonical(q) for q in product(range(-T, T+1), repeat=3)
                   if any(q) and gcd(*q) == 1 and height2(q) <= T*T}
        from_planes = {canonical((p[2], -p[1], p[0])) for p in planes}
        assert normals == from_planes
        counts['independent_normal_plane_enumerations'] += 1
    if d == 4:
        for p in planes:
            assert p[0]*p[5]-p[1]*p[4]+p[2]*p[3] == 0
            dual = canonical((p[5], -p[4], p[3], p[2], -p[1], p[0]))
            assert dual in planes
            counts['exact_plucker_and_dual_controls'] += 1
    counts['complete_bounded_height_plane_enumerations'] += 1
    counts['reduced_primitive_basis_pairs'] += pairs
    return {'dimension': d, 'height_cutoff': T, 'plane_count': len(planes),
            'count_over_height_to_dimension': str(Q(len(planes), T**d)),
            'primitive_plucker_coordinates': sorted(planes)}


def short_planes(counts):
    fixtures = []
    for d in range(3, 7):
        candidates = list(product((-1, 0, 1), repeat=d))
        ws = [w for w in product((-2, -1, 0, 1, 2), repeat=d)
              if canonical(w) == w and sum(map(abs, w)) <= 4]
        for w in ws:
            ps = [canonical(wedge(w, u)) for u in candidates]
            p = min((p for p in ps if p), key=height2)
            assert contains(p, w)
            # Square and raise to an integer power: no floating roots.
            assert height2(p)**(d-1) <= (4*(d-1))**(d-1)*height2(w)**(d-2)
            counts['exact_short_containing_planes'] += 1
        fixtures.append({'dimension': d, 'primitive_objectives': len(ws),
                         'scope': 'Complete l1<=4 subset of the coordinate cube [-2,2]^d.'})
    return fixtures


def sharpness(counts):
    rng = Random(20417)
    fixtures = []
    for d in range(3, 9):
        for t in (2, 3, 7, 100, 10**12):
            w = tuple(t**i for i in range(d))
            u = (0,)+tuple(t**i for i in range(d-1))
            p = canonical(wedge(w, u))
            assert contains(p, w) and max(map(abs, p)) == t**(d-2)
            assert height2(p) == sum(t**(2*i) for i in range(d-1))
            for _ in range(30):
                v = tuple(rng.randint(-2*t-3, 2*t+3) for _ in range(d))
                a = canonical(tuple(v[i]-v[0]*w[i] for i in range(d)))
                if a is None:
                    continue
                raw = wedge(w, a)
                assert gcd(*raw) == 1
                assert canonical(raw) == canonical(wedge(w, v))
                failures = [j for j in range(1, d-1) if a[j+1] != t*a[j]]
                if failures:
                    j = max(failures)
                    assert w[j]*a[-1]-w[-1]*a[j] == t**(d-2)*(a[j+1]-t*a[j])
                else:
                    assert abs(a[-1]) >= t**(d-2)
                assert max(map(abs, raw)) >= t**(d-2)
                counts['integer_sharpness_recurrence_controls'] += 1
            fixtures.append({'dimension': d, 't': t, 'primitive_objective': w,
                             'containing_plane_plucker': p,
                             'exact_minimum_plucker_infinity_height': t**(d-2)})
    return fixtures


def det_permutations(a):
    n = len(a)
    answer = 0
    for sigma in permutations(range(n)):
        sign = (-1)**sum(sigma[i] > sigma[j] for i in range(n) for j in range(i+1, n))
        term = sign
        for i in range(n):
            term *= a[i][sigma[i]]
        answer += term
    return answer


def matroid_controls(counts):
    rng = Random(91643)
    records = []
    for d in range(3, 7):
        H, rank = 3, 3
        ground = tuple(sorted(hyperbolic_cross(H, d)))
        for label, prime in (('cyclic_quotient', 2), ('linear_matroid', 101)):
            if label == 'cyclic_quotient':
                columns = [tuple(int(sum((i+1)*a[i] for i in range(d)) == j)
                                 for j in range(rank)) for a in ground]
            else:
                columns = [tuple(int(i == j) for i in range(rank)) if j < rank else
                           tuple(rng.randrange(prime) for _ in range(rank)) for j in range(len(ground))]
            oracle = matrix_oracle(columns, rank, prime)
            bases = []
            for b in combinations(range(len(ground)), rank):
                independent = det_permutations([[columns[j][i] for j in b] for i in range(rank)]) % prime != 0
                assert independent == oracle(frozenset(b))
                counts['independent_basis_determinants'] += 1
                if independent:
                    bases.append(b)
            profiles = {tuple(sum(ground[j][i] for j in b) for i in range(d)) for b in bases}
            t = H+1
            u = (0,)+tuple(t**i for i in range(d-1))
            A = (tuple(int(i == 0) for i in range(d)), u)
            delta, Delta = projection_data(A)
            assert delta == 1 and Delta == sum(t**i for i in range(d-1))
            assert Delta*Delta <= (d*(d-1)//2)*height2(canonical(wedge(*A)))
            hulls = []
            for projection in (A, (tuple(2*x for x in A[0]), A[1])):
                pdelta, pheight = projection_data(projection)
                assert pheight == Delta
                assert pdelta == (1 if projection == A else 2)
                points = {tuple(dot(row, v) for row in projection) for v in profiles}
                hull = convex_hull(points)
                quadrants = [sweep(ground, projection, rank, oracle, H, signs)
                             for signs in product((-1, 1), repeat=2)]
                exposed = {tuple(state['profile']) for q in quadrants for state in q['states']}
                assert exposed == set(hull)
                cube = convex_hull(tuple(dot(row, v) for row in projection)
                                   for v in product((0, 1), repeat=d))
                regularized = convex_hull(tuple(x+y for x, y in zip(a, b)) for a in hull for b in cube)
                assert len(regularized) <= len(hull)+2*d
                costs = [dot((1, t), v) for v in hull] if projection == A else None
                if costs:
                    assert costs.count(min(costs)) == 1
                counts['full_projected_matroid_hulls'] += 1
                counts['generic_dimension_quadrant_sweeps'] += 4
                counts['interval_states'] += sum(len(q['states']) for q in quadrants)
                counts['regularizing_cube_shadow_controls'] += 1
                hulls.append({'rows': projection, 'minor_gcd': pdelta, 'normalized_minor_sum': pheight,
                              'independent_profile_hull': hull, 'regularized_vertices': len(regularized)})
            assert len(hulls[0]['independent_profile_hull']) == len(hulls[1]['independent_profile_hull'])
            records.append({'dimension': d, 'model': label, 'characteristic': prime,
                            'H': H, 'rank': rank, 'ground': ground, 'columns': columns,
                            'independent_bases': len(bases), 'projections': hulls})
    return records


def exponents(counts):
    records = []
    for d in range(3, 129):
        beta = Q(3*(d-1), (d+1)*(3*d-5))
        alpha = Q((d-2)*(3*d+1), 3*(d-1))
        e = 2*(d-1)-(2*d-3)*beta
        classical = Q(2*d*(d-1), d+1)
        old = Q((d-1)*(6*d*(d-2)+4), 3*(d*d-d-1))
        assert beta == 1/(1+alpha)
        assert classical-e == Q(d-1, (d+1)*(3*d-5)) > 0
        assert old-e == -Q((d-1)*(3*d*d-19*d+23), 3*(d+1)*(3*d-5)*(d*d-d-1))
        assert (e < old) == (d in (3, 4))
        counts['exact_dimension_exponent_comparisons'] += 1
        records.append({'dimension': d, 'classical_power': str(classical), 'old_Dey_power': str(old),
                        'rational_plane_power': str(e), 'logarithmic_power': str(beta),
                        'smaller_power': str(min(e, old))})
    assert records[1]['rational_plane_power'] == '33/7'
    assert Q(52, 11)-Q(33, 7) == Q(1, 77)
    assert 4-Q(6, 5) == Q(14, 5)
    assert 4-6/(5+Q(1, 3)) == Q(23, 8)
    return records


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    counts = defaultdict(int)
    planes = [plane_enumeration(d, T, counts) for d in (3, 4, 5) for T in (1, 2, 3)]
    short = short_planes(counts)
    sharp = sharpness(counts)
    models = matroid_controls(counts)
    powers = exponents(counts)
    sources = ('docs/rational_plane_transfer.md', 'scripts/run_rational_plane_transfer_checks.py',
               'docs/rational_shadow_bound.md', 'scripts/run_rational_shadow_bound_checks.py',
               'scripts/run_planar_dyadic_bound_checks.py', 'scripts/run_structural_complexity_checks.py')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_rational_plane_transfer_checks.py',
              'counts': dict(counts), 'bounded_height_planes': planes, 'containing_plane_controls': short,
              'infinite_sharpness_family_controls': sharp, 'matroid_projection_controls': models,
              'exact_exponents': powers,
              'results': {'all_fixed_dimensions': True, 'four_variable_power': '33/7',
                          'three_variable_power_improved': False,
                          'unchanged_transfer_budget_minimum_three_variable_power': '14/5',
                          'higher_dimensional_old_bound_retained': 'd>=5',
                          'novelty_established': False, 'goal_complete': False},
              'scope': 'Exact bounded controls support the uniform geometric and matroid proofs. No count fitting, priority certificate, new three-variable exponent or completed publication goal.',
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/'results/research/rational_plane_transfer_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(counts), indent=2))


if __name__ == '__main__':
    main()
