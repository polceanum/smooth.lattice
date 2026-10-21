#!/usr/bin/env python3
"""Exact controls for two bilinear rank defects and the common-factor case."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product, combinations
from math import comb
from pathlib import Path
import hashlib
import json

from run_signed_face_checks import field_basis, field_reduce, normals
from run_structural_complexity_checks import area2, clip, dot, finite_staircases, hyperbolic_cross
from run_wall_rank_filter_checks import rank

ROOT = Path(__file__).resolve().parents[1]


def polynomial_controls(counts):
    def plus(*terms):
        out = defaultdict(int)
        for p in terms:
            for a, c in p.items():
                out[a] += c
        return {a: c for a, c in out.items() if c}
    def times(p, q):
        out = defaultdict(int)
        for a, c in p.items():
            for b, d in q.items():
                out[tuple(x+y for x, y in zip(a, b))] += c*d
        return {a: c for a, c in out.items() if c}
    neg = lambda p: {a: -c for a, c in p.items()}
    symbols = [{tuple(int(i == j) for j in range(9)): 1} for i in range(9)]
    a, b, c, d, e, f, x, y, z = symbols
    X, Y, Z = plus(x, neg(a)), plus(y, neg(b)), plus(z, neg(e))
    F, G = plus(times(X, Y), neg(c)), plus(times(plus(x, neg(d)), Z), neg(f))
    left = plus(times(Y, G), neg(times(Z, F)))
    right = plus(times(plus(a, neg(d)), times(Y, Z)), times(c, Z), neg(times(f, Y)))
    assert left == right
    counts['independent_integral_bilinear_elimination_identities'] += 1
    for p in (2, 3, 5):
        for a, d, c, f in product(range(p), repeat=4):
            assert (not any(((a-d) % p, c, f))) == (a == d and c == f == 0)
            counts['exhaustive_normalized_coefficient_degeneracies'] += 1


def mul(a, b, f, p):
    n = len(f)-1
    out = [0]*(2*n-1)
    for i, c in enumerate(a):
        for j, d in enumerate(b):
            out[i+j] += c*d
    for i in range(len(out)-1, n-1, -1):
        c = out[i]
        for j in range(n):
            out[i-n+j] -= c*f[j]
    return tuple(x % p if p else x for x in out[:n])


def divide_root(f, a, p):
    n = len(f)-1
    q = [0]*n
    q[-1] = 1
    for j in range(n-2, -1, -1):
        q[j] = (f[j+1]+a*q[j+1]) % p if p else f[j+1]+a*q[j+1]
    remainder = f[0]+a*q[0]
    return tuple(q), remainder % p if p else remainder


def reciprocal(f, a, p):
    q, value = divide_root(f, a, p)
    assert value
    scale = -pow(int(value), -1, p) if p else -Q(1, value)
    return tuple(x*scale % p if p else x*scale for x in q)


def as_column(v):
    return {j: c for j, c in enumerate(v) if c}


def monogenic_model(f, y, z, p):
    n = len(f)-1
    one = (1,)+(0,)*(n-1)
    x = (0, 1)+(0,)*(n-2)
    @lru_cache(None)
    def ev(a):
        value = one
        for h, e in zip((x, y, z), a):
            for _ in range(e):
                value = mul(value, h, f, p)
        return as_column(value)
    return n, ev


def inverse_controls(counts):
    rows = []
    for p, n in ((2, 3), (2, 4), (2, 5), (3, 3)):
        total, low = 0, 0
        for coefficients in product(range(p), repeat=n):
            f = coefficients+(1,)
            one, x = (1,)+(0,)*(n-1), (0, 1)+(0,)*(n-2)
            classes = {tuple([c]+[0]*(n-1)): 'constant' for c in range(p)}
            for alpha, beta, gamma, delta in product(range(p), repeat=4):
                if (alpha*delta-beta*gamma) % p == 0:
                    continue
                numerator = (alpha, beta)+(0,)*(n-2)
                if delta:
                    a = -gamma*pow(delta, -1, p) % p
                    _, value = divide_root(f, a, p)
                    if not value:
                        continue
                    inv = tuple(c*pow(delta, -1, p) % p for c in reciprocal(f, a, p))
                else:
                    if not gamma:
                        continue
                    inv = (pow(gamma, -1, p),)+(0,)*(n-1)
                h = mul(numerator, inv, f, p)
                assert h not in classes or classes[h] == 'mobius'
                classes[h] = 'mobius'
            for a in range(p):
                q, value = divide_root(f, a, p)
                if value:
                    continue
                for scalar, scale in product(range(p), range(1, p)):
                    h = tuple((scale*c+(scalar if j == 0 else 0)) % p for j, c in enumerate(q))
                    assert h not in classes or classes[h] == 'rank_one'
                    classes[h] = 'rank_one'
                    residual = tuple((c-(scalar if j == 0 else 0)) % p for j, c in enumerate(h))
                    images = [as_column(mul(residual, tuple(int(j == i) for j in range(n)), f, p))
                              for i in range(n)]
                    assert rank(images, p) == 1
                    counts['independent_rank_one_multiplication_images'] += 1
            for h in product(range(p), repeat=n):
                dimension = rank([as_column(v) for v in (one, x, h, mul(x, h, f, p))], p)
                assert (dimension <= 3) == (h in classes)
                total += 1
                low += dimension <= 3
                counts['exhaustive_monogenic_four_column_rank_tests'] += 1
        rows.append({'characteristic': p, 'length': n, 'monic_polynomials': p**n,
                     'all_elements_tested': total, 'rank_at_most_three': low})
    return rows


def full_fan(label, n, ev, p, counts, expected=None, axis_only=False, plane_axis=False):
    one, x, y, z, xy, xz = (ev(a) for a in ((0, 0, 0), (1, 0, 0), (0, 1, 0),
                                                    (0, 0, 1), (1, 1, 0), (1, 0, 1)))
    assert rank([one, x, y, xy], p) <= 3 and rank([one, x, z, xz], p) <= 3
    E = tuple(hyperbolic_cross(n, 3))
    reference, pivots = field_basis(tuple(sorted(E)), ev, p)
    assert len(reference) == n
    def coordinates(v):
        rem, co = field_reduce(v, pivots, p)
        assert not rem
        return tuple(co.get(a, 0) for a in reference)
    matrices = [tuple(coordinates(ev(tuple(x+(i == j) for j, x in enumerate(a))))
                      for a in reference) for i in range(3)]
    def apply(M, v):
        out = tuple(sum(c*column[i] for column, c in zip(M, v)) for i in range(n))
        return tuple(c % p for c in out) if p else out
    for i, j in combinations(range(3), 2):
        for k in range(n):
            assert apply(matrices[i], matrices[j][k]) == apply(matrices[j], matrices[i][k])
            counts['independent_commuting_operator_columns'] += 1
    @lru_cache(None)
    def replay(a):
        if not any(a):
            return coordinates(ev(a))
        i = next(i for i, e in enumerate(a) if e)
        return apply(matrices[i], replay(tuple(e-(i == j) for j, e in enumerate(a))))
    for a in E:
        assert replay(a) == coordinates(ev(a))
        counts['independent_cyclic_matrix_monomial_replays'] += 1
    cells = []
    for B in finite_staircases(n, 3):
        counts['independent_exhaustive_standard_set_candidates'] += 1
        chosen, _ = field_basis(tuple(sorted(B)), ev, p)
        if len(chosen) != n:
            continue
        constraints, _ = normals(tuple(sorted(B)), ev, p)
        poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for q in constraints:
            poly = clip(poly, q, 1)
        if len(poly) < 3 or not area2(poly):
            continue
        witness = tuple(sum(v[i] for v in poly)/len(poly) for i in range(3))
        assert min(witness) > 0 and all(dot(witness, q) > 0 for q in constraints)
        replay, _ = field_basis(sorted(E, key=lambda a: (dot(witness, a), a)), ev, p)
        assert set(replay) == set(B)
        if axis_only:
            assert all(sum(c != 0 for c in a) <= 1 for a in B)
        cells.append((frozenset(B), poly, witness))
    assert sum(area2(poly) for _, poly, _ in cells) == 1
    if expected is not None:
        assert len(cells) == expected, (label, p, len(cells), expected)
    counts['independent_complete_actual_quotient_fans'] += 1
    groups = []
    if plane_axis:
        assert not xy and not xz
        r = rank([ev((j, 0, 0)) for j in range(1, n+1)], p)
        values = sorted({sum(a[1] == a[2] == 0 and a[0] > 0 for a in B) for B, _, _ in cells})
        assert set(values) <= {r-1, r}
        for a in values:
            C = tuple((j, 0, 0) for j in range(a+1))
            chosen, pivots = field_basis(C, ev, p)
            assert len(chosen) == a+1
            k = n-a-1
            ground = tuple((0, b, c) for b, c in hyperbolic_cross(k+1, 2) if b or c)
            residual = {v: field_reduce(ev(v), pivots, p)[0] for v in ground}
            assert rank(list(residual.values()), p) == k
            found = 0
            for B, _, w in cells:
                if sum(v[1] == v[2] == 0 and v[0] > 0 for v in B) != a:
                    continue
                assert set(C) <= B
                selected, _ = field_basis(sorted(ground, key=lambda v: (dot(w, v), v)), residual.__getitem__, p)
                assert set(selected) == B-set(C)
                counts['independent_fixed_axis_planar_contraction_replays'] += 1
                found += 1
            groups.append({'axis_length': a, 'axis_ideal_rank': r, 'residual_rank': k,
                           'actual_positive_minima': found})
    return {'model': label, 'length': n, 'characteristic': p, 'initial_ideals': len(cells),
            'fixed_axis_classes': groups,
            'complete_cells': [{'standard_set': sorted(B), 'weight': list(map(str, w)),
                                'positive_cone': [list(map(str, v)) for v in poly]}
                               for B, poly, w in cells]}


def monogenic_fans(counts):
    rows = []
    n = 4
    for p in (0, 2, 3, 5):
        # (t^2-t+1)^2 avoids 0,1 and has repeated roots.
        f = (1, -2, 3, -2, 1)
        f = tuple(c % p for c in f) if p else f
        _, ev = monogenic_model(f, reciprocal(f, 0, p), reciprocal(f, 1, p), p)
        rows.append(full_fan('three_poles_repeated_divisor', n, ev, p, counts,
                             expected=comb(n+1, 2), axis_only=True))
        f = (0,)*n+(1,)
        x, tag = (0, 1, 0, 0), (0, 0, 0, 1)
        for label, y, z, expected in (
                ('one_pole', x, x, 3), ('two_poles_duplicate', x, reciprocal(f, 1, p), 2*n-1),
                ('socle_and_second_pole', tag, reciprocal(f, 1, p), 2*n-1),
                ('socle_and_duplicate', tag, x, 4), ('two_proportional_socles', tag, tag, 3)):
            _, ev = monogenic_model(f, y, z, p)
            rows.append(full_fan(label, n, ev, p, counts, expected=expected, axis_only=True))
        f = (1, 0, 0, 0, 1)
        _, ev = monogenic_model(f, (0,)*n, (0,)*n, p)
        rows.append(full_fan('axis_unit_plane_zero', n, ev, p, counts,
                             expected=1, axis_only=True, plane_axis=True))
    return rows


def plane_axis_fans(counts):
    rows = []
    for R, a, p in ((4, 2, 0), (4, 2, 2), (4, 2, 3), (4, 3, 5), (5, 2, 2)):
        n = R+a-1
        @lru_cache(None)
        def ev(alpha):
            i, j, k = alpha
            if i:
                if j or k or i > a:
                    return {}
                return {R-1: 1} if i == a else {R+i-1: 1}
            # y=t, z=t+t^2 in the local planar quotient K[t]/t^R.
            value = [1]+[0]*(R-1)
            for h, e in (((0, 1)+(0,)*(R-2), j), ((0, 1, 1)+(0,)*(R-3), k)):
                for _ in range(e):
                    value = mul(value, h, (0,)*R+(1,), p)
            return as_column(value)
        base = full_fan('plane_axis_socle_gluing_%s_%s' % (R, a), n, ev, p, counts, plane_axis=True)
        rows.append(base)
        shifts = (1, 2, 3)
        # Independent binomial expansion of coordinate translations.
        @lru_cache(None)
        def translated(alpha):
            out = defaultdict(int)
            for beta in product(*(range(e+1) for e in alpha)):
                c = 1
                for e, b, s in zip(alpha, beta, shifts):
                    c *= comb(e, b)*s**(e-b)
                for j, d in ev(beta).items():
                    out[j] += c*d
            return {j: c % p if p else c for j, c in out.items() if (c % p if p else c)}
        moved = full_fan('translated_'+base['model'], n, translated, p, counts)
        base_cells = {tuple(map(tuple, c['standard_set'])): c['positive_cone'] for c in base['complete_cells']}
        moved_cells = {tuple(map(tuple, c['standard_set'])): c['positive_cone'] for c in moved['complete_cells']}
        # Vertex order may change after redundant halfspaces are added.
        assert set(base_cells) == set(moved_cells)
        for B in base_cells:
            assert set(map(tuple, base_cells[B])) == set(map(tuple, moved_cells[B]))
        counts['independent_whole_fan_affine_shift_equalities'] += 1
        rows.append(moved)
    # A direct product tests the invertible T branch of the shared-dimension
    # proof, independently of the local socle identifications above.
    for p in (0, 2, 3):
        r, R = 2, 3
        n = r+R
        @lru_cache(None)
        def ev(alpha):
            i, j, k = alpha
            if i and (j or k):
                return {}
            if i:
                return {i % r: 1}
            value = [1]+[0]*(R-1)
            for h, exponent in (((0, 1, 0), j), ((1, 1, 0), k)):
                for _ in range(exponent):
                    value = mul(value, h, (0,)*R+(1,), p)
            out = {r+j: c for j, c in enumerate(value) if c}
            if not j and not k:
                out[0] = 1
            return out
        rows.append(full_fan('disjoint_axis_plane_components', n, ev, p, counts, plane_axis=True))
    n, p = 8, 1000000007
    @lru_cache(None)
    def ev(alpha):
        i, j, k = alpha
        return {} if i else {t-1: pow(t, j+n*k, p) for t in range(1, n+1)}
    rows.append(full_fan('scalar_axis_generic_planar_points', n, ev, p, counts,
                         expected=10, plane_axis=True))
    return rows


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    counts = defaultdict(int)
    polynomial_controls(counts)
    inverse = inverse_controls(counts)
    fans = monogenic_fans(counts)+plane_axis_fans(counts)
    paths = ('docs/bilinear_rank_defects.md', 'scripts/run_bilinear_rank_defect_checks.py',
             'docs/planar_dyadic_bound.md', 'scripts/run_wall_rank_filter_checks.py',
             'scripts/run_signed_face_checks.py', 'scripts/run_structural_complexity_checks.py')
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_bilinear_rank_defect_checks.py',
              'scope': 'Integral bilinear elimination, exhaustive finite monogenic four-column inverse tests, complete actual fans and the exceptional plane-plus-axis contractions. These support the uniform all-field proof, not an unrestricted G3 theorem, a power-separating family or novelty.',
              'counts': dict(sorted(counts.items())), 'monogenic_inverse_rank_tests': inverse,
              'complete_actual_fans': fans,
              'uniform_results': {'two_pair_defects': 'rank(1,x,y,xy)<=3 and rank(1,x,z,xz)<=3 imply g=O(n^2), over all fields and arbitrary finite quotients.',
                 'dichotomy': 'Either all three pair defects hold, a coordinate is scalar, or a coordinatewise affine shift puts xy,xz in I.',
                 'common_factor_case': 'For I containing xy,xz, only two x-axis lengths occur and one fixed-core planar matroid per length gives g=O(n*log(n+1)).',
                 'cyclic_inverse': 'With cyclic x and n>=3, the low-rank coordinates are constants, unit-denominator Mobius functions, or scalars plus rank-one annihilator elements.',
                 'sharpness': 'The earlier three-pole family attains binom(n+1,2), over every field, including repeated divisors.',
                 'unrestricted_G3_improved': False, 'polynomial_lower_separation_proved': False,
                 'novelty_established': False, 'goal_complete': False},
              'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in paths}}
    (ROOT/'results/research/bilinear_rank_defect_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'counts': report['counts'], 'complete_actual_fans': len(fans)}, indent=2))


if __name__ == '__main__':
    main()
