#!/usr/bin/env python3
"""Exact controls for numerical-recall capacity, not a novelty or ML test."""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

from run_retention_tradeoff_checks import inverse, multiply, plus, scaled, transpose
from run_robust_reader_checks import unit_circle
from run_robust_frontier_checks import key, matrix, pairs

ROOT = Path(__file__).resolve().parents[1]


def dot(x, y):
    return sum((a*b for a, b in zip(x, y)), Q(0))


def mv(a, x):
    return [dot(row, x) for row in a]


def outer(x):
    return [[a*b for b in x] for a in x]


def total(matrices):
    answer = scaled(matrices[0], Q(0))
    for a in matrices:
        answer = plus(answer, a)
    return answer


def trace(a):
    return sum((row[i] for i, row in enumerate(a)), Q(0))


def evaluate(s, b, c, x):
    return dot(x, mv(s, x))+dot(b, x)+c


def pseudoinverse(g):
    """Rational column-basis construction; verify all Penrose identities."""
    a, pivots, row = [x[:] for x in g], [], 0
    for column in range(len(g)):
        pivot = next((i for i in range(row, len(g)) if a[i][column]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        value = a[row][column]
        a[row] = [x/value for x in a[row]]
        for i in range(len(g)):
            if i != row:
                value = a[i][column]
                a[i] = [x-value*y for x, y in zip(a[i], a[row])]
        pivots.append(column)
        row += 1
    if not pivots:
        return scaled(g, Q(0)), 0
    basis = [[x[j] for j in pivots] for x in g]
    left = multiply(inverse(multiply(transpose(basis), basis)), transpose(basis))
    middle = multiply(multiply(left, g), transpose(left))
    gp = multiply(multiply(transpose(left), inverse(middle)), left)
    gg, pg = multiply(g, gp), multiply(gp, g)
    assert multiply(gg, g) == g and multiply(pg, gp) == gp
    assert transpose(gg) == gg and transpose(pg) == pg
    return gp, len(pivots)


def symmetric(p, rng):
    a = [[Q(rng.randrange(-8, 9), 7) for _ in range(p)] for _ in range(p)]
    return scaled(plus(a, transpose(a)), Q(1, 2))


def frames():
    circle = unit_circle()[:6]
    mixed = [[Q(int(i == j)) for j in range(3)] for i in range(3)]
    mixed += [[Q(3, 5), Q(4, 5), Q(0)], [Q(0), Q(3, 5), Q(4, 5)],
              [Q(4, 5), Q(0), Q(3, 5)]]
    answer = [circle, mixed, [[Q(i, 4), Q(0), Q(0), Q(0)] for i in (-3, -1, 1, 3)],
              [list(x)+[Q(0)]*3 for x in circle], [[Q(1, 5), Q(2, 5), Q(0)]]*3]
    for bits in (4, 20, 80):
        answer.append([[Q(1, 5)+Q(i, 3), Q(j, 2**bits), Q(1, 10), Q(0)]
                       for i, j in product((-1, 1), repeat=2)])
    return answer


def centered_controls():
    rng, assignments, cauchy, differences, ranks = random.Random(70031), 0, 0, 0, []
    for points in frames():
        n, p = len(points), len(points[0])
        assert all(dot(u, u) <= 1 for u in points)
        mean = [sum(u[j] for u in points)/n for j in range(p)]
        xs = [[x-m for x, m in zip(u, mean)] for u in points]
        g = total([outer(x) for x in xs])
        gp, rank = pseudoinverse(g)
        ranks.append(rank)
        leverage = [dot(x, mv(gp, x)) for x in xs]
        assert sum(leverage) == rank and trace(g) <= n
        expected = sum(dot(x, x)*ell for x, ell in zip(xs, leverage))/4+Q(3, 4*n)*trace(g)
        assert expected <= rank+1
        accumulated = Q(0)
        states = []
        for _ in range(3):
            s, b, c = symmetric(p, rng), [Q(rng.randrange(-4, 5), 9) for _ in range(p)], Q(13)
            t = [2*x+y for x, y in zip(mv(s, mean), b)]
            energy = sum(dot([2*v+w for v, w in zip(mv(s, u), b)],
                             [2*v+w for v, w in zip(mv(s, u), b)]) for u in points)
            assert energy == 4*trace(multiply(multiply(s, s), g))+n*dot(t, t)
            for u in points:
                for j in range(p):
                    h = [Q(int(k == j), 17) for k in range(p)]
                    positive, negative = [x+y for x, y in zip(u, h)], [x-y for x, y in zip(u, h)]
                    assert evaluate(s, b, c, positive)-evaluate(s, b, c, negative) == 2*dot(h, [2*v+w for v, w in zip(mv(s, u), b)])
                    differences += 1
            states.append((s, b, c, t, energy))
        correlation = Q(0)
        for labels in product((-1, 1), repeat=n):
            average = Q(sum(labels), n)
            z = [Q(y)-average for y in labels]
            m = total([scaled(outer(x), weight) for x, weight in zip(xs, z)])
            v = [sum(weight*x[j] for x, weight in zip(xs, z)) for j in range(p)]
            dual2 = trace(multiply(gp, multiply(m, m)))/4+dot(v, v)/n
            assert dual2 >= 0
            accumulated += dual2
            assert sum(abs(weight) for weight in z) == n*(1-average*average)
            correlation += sum(weight*y for weight, y in zip(z, labels))
            for s, b, c, t, energy in states:
                actual = sum(weight*evaluate(s, b, c, u) for weight, u in zip(z, points))
                assert actual == trace(multiply(s, m))+dot(t, v)
                assert actual*actual <= energy*dual2
                cauchy += 1
            assignments += 1
        assert accumulated/2**n == expected
        assert correlation/2**n == n-1
    return {'frames': len(ranks), 'centered_ranks': ranks, 'sign_assignments': assignments,
            'weighted_Cauchy_controls': cauchy, 'ambient_central_difference_identities': differences}


def numerical_fixtures():
    states = balls = 0
    a, epsilon, delta = Q(1, 2), Q(1, 8), Q(1, 32)
    for p in range(2, 7):
        points = [[Q(int(i == j)) for j in range(p)] for i in range(p)]
        for labels, c in product(product((-1, 1), repeat=p), (Q(-1, 100), Q(1, 100))):
            b = [Q(j+1, 1000) for j in range(p)]
            s = [[(a*labels[i]-b[i]-c)*(i == j) for j in range(p)] for i in range(p)]
            norm = max(abs(s[i][i]) for i in range(p))
            variation = 2*norm*delta+norm*delta**2+sum(abs(x) for x in b)*delta
            assert variation < epsilon
            for u, label in zip(points, labels):
                assert evaluate(s, b, c, u) == a*label
                grad = [2*x+y for x, y in zip(mv(s, u), b)]
                assert dot(grad, grad) <= (epsilon/delta)**2
                balls += 1
            # The general upper bound, applied to this explicitly certified family.
            assert Q((p-1)**2, p) <= (p+1)*epsilon**2/(delta**2*(a-epsilon)**2)
            states += 1
    # Block lower construction at a fixed 1/4 accuracy for all bounded labels.
    counts = 0
    for p, delta in product((2, 3, 7, 16, 64, 257), (Q(1, 16), Q(1, 32), Q(1, 64))):
        r = min(p, 1+int(1/(256*delta**2)))
        assert len(pairs(p, r)) >= Q(p*(r-1), 4)
        counts += 1
    lower_error = Q(35, 192)+Q(25, 6144)
    assert lower_error == Q(1145, 6144) < Q(1, 4)
    return {'whole_ball_states': states, 'certified_key_balls': balls,
            'fixed_accuracy_capacity_counts': counts, 'block_lower_error': str(lower_error)}


def spherical_controls():
    rng, moments, differences = random.Random(3107), 0, 0
    unit_frames = [frame for frame in frames() if all(dot(u, u) == 1 for u in frame)]
    for points in unit_frames:
        n, p = len(points), len(points[0])
        g = total([outer(u) for u in points])
        gp, rank = pseudoinverse(g)
        average = Q(0)
        for labels in product((-1, 1), repeat=n):
            m = total([scaled(outer(u), y) for u, y in zip(points, labels)])
            dual2 = trace(multiply(gp, multiply(m, m)))
            average += dual2
            s = symmetric(p, rng)
            energy = sum(dot(mv(s, u), mv(s, u)) for u in points)
            assert trace(multiply(s, m))**2 <= energy*dual2
            moments += 1
        assert average/2**n == rank
    rotations = [(Q(15, 17), Q(8, 17)), (Q(99, 101), Q(20, 101)),
                 (Q(255, 257), Q(32, 257))]
    for u, (co, si), radial in product(unit_circle()[:6], rotations, (Q(-2**64), Q(0), Q(2**64))):
        s = symmetric(2, rng)
        shifted = [[s[i][j]-radial*(i == j) for j in range(2)] for i in range(2)]
        w, delta = [-u[1], u[0]], 2*si
        positive, negative = ([co*x+side*si*y for x, y in zip(u, w)] for side in (1, -1))
        assert dot(positive, positive) == dot(negative, negative) == 1
        assert sum((x-y)**2 for x, y in zip(positive, u)) <= delta**2
        fp, fm = (evaluate(shifted, [Q(0), Q(0)], radial, x) for x in (positive, negative))
        center = dot(u, mv(s, u))
        epsilon = max(abs(fp-center), abs(fm-center))
        tangent = dot(w, mv(s, u))
        assert fp-fm == 4*co*si*tangent
        assert tangent*tangent <= Q(4, 3)*epsilon**2/delta**2
        assert dot(mv(s, u), mv(s, u)) <= center**2+Q(4, 3)*epsilon**2/delta**2
        differences += 1
    # A six-address block family exceeds p addresses and obeys both numerical bounds.
    p, addresses, a, epsilon, delta = 4, pairs(4, 4), Q(1, 2), Q(1, 8), Q(1, 24)
    for labels in product((-1, 1), repeat=len(addresses)):
        s = matrix(p, {pair: a*y for pair, y in zip(addresses, labels)})
        assert a*(Q(35, 12)*delta*2+Q(25, 24)*delta**2*3) < epsilon
        assert all(dot(key(p, pair), mv(s, key(p, pair))) == a*y for pair, y in zip(addresses, labels))
        assert len(addresses) <= p*((a+epsilon)**2+Q(4, 3)*epsilon**2/delta**2)/(a-epsilon)**2
    return {'weighted_sign_moment_assignments': moments,
            'spherical_odd_difference_and_radial_shift_controls': differences,
            'six_address_block_labelings': 2**len(addresses)}


def clipping_controls():
    checked = 0
    for a, epsilon, raw in product((Q(1, 4), Q(1, 2)), (Q(1, 16), Q(1, 8)),
                                  (Q(i, 8) for i in range(-24, 25))):
        for sign in (-1, 1):
            clipped = max(Q(-1), min(Q(1), raw))
            if abs(clipped-sign*a) <= epsilon:
                assert raw == clipped and abs(raw-sign*a) <= epsilon
            checked += 1
    # Binary rail recall can be exact after clipping despite nonzero gradient.
    assert all(max(Q(-1), min(Q(1), 8*x)) == sign
               for sign in (-1, 1) for x in (sign-Q(1, 4), sign, sign+Q(1, 4)))
    return {'interior_value_implications': checked, 'binary_rail_gradient_counterexample': True}


def projector(a):
    """Compute a column-space projector independently of a coefficient choice."""
    gp, rank = pseudoinverse(multiply(transpose(a), a))
    p = multiply(multiply(a, gp), transpose(a))
    assert p == transpose(p) and multiply(p, p) == p
    assert multiply(p, a) == a and trace(p) == rank
    return p, rank


def value_and_half_tangent_gradient(points, s, b, c):
    values, derivatives = [], []
    for u in points:
        assert dot(u, u) == 1
        values.append(evaluate(s, b, c, u))
        v = [x+y/2 for x, y in zip(mv(s, u), b)]
        radial = dot(u, v)
        derivatives.extend(x-radial*y for x, y in zip(v, u))
    return values+derivatives


def sphere_jet_matrix(points, affine):
    p, columns = len(points[0]), []
    for i in range(p):
        for j in range(i, p):
            s = [[Q(0)]*p for _ in range(p)]
            s[i][j] = s[j][i] = Q(1)
            columns.append(value_and_half_tangent_gradient(points, s, [Q(0)]*p, Q(0)))
    if affine:
        for i in range(p):
            b = [Q(int(i == j)) for j in range(p)]
            columns.append(value_and_half_tangent_gradient(points, scaled(s, Q(0)), b, Q(0)))
        columns.append(value_and_half_tangent_gradient(points, scaled(s, Q(0)), [Q(0)]*p, Q(1)))
    return transpose(columns)


def affine_sphere_projection_controls():
    rng, cases, assignments, inequalities = random.Random(2047), [], 0, 0
    points_list = [f for f in frames() if all(dot(u, u) == 1 for u in f)]
    points_list.append([key(4, pair) for pair in pairs(4, 4)])
    for bits in (4, 24):
        parameters = [Q(j, 2**bits) for j in (-2, -1, 1, 2)]
        points_list.append([[(1-t*t)/(1+t*t), 2*t/(1+t*t), Q(0)] for t in parameters])
    points_list.append([[Q(1), Q(0), Q(0)]]*3)
    for points in points_list:
        n, p = len(points), len(points[0])
        a0, a = sphere_jet_matrix(points, False), sphere_jet_matrix(points, True)
        p0, rank0 = projector(a0)
        proj, rank = projector(a)
        difference = plus(proj, scaled(p0, Q(-1)))
        assert multiply(proj, p0) == multiply(p0, proj) == p0
        assert difference == transpose(difference) and multiply(difference, difference) == difference
        assert all(0 <= difference[i][i] <= 1 for i in range(len(difference)))
        # The constant jet is already the quadratic jet of the identity matrix.
        assert trace(difference) == rank-rank0 <= p
        _, key_rank = pseudoinverse(total([outer(u) for u in points]))
        value_trace0 = sum(p0[i][i] for i in range(n))
        value_trace = sum(proj[i][i] for i in range(n))
        assert value_trace0 <= key_rank
        assert value_trace <= value_trace0+rank-rank0 <= key_rank+p
        states = []
        for _ in range(2):
            s, b = symmetric(p, rng), [Q(rng.randrange(-9, 10), 11) for _ in range(p)]
            output = value_and_half_tangent_gradient(points, s, b, Q(2**32))
            assert mv(proj, output) == output
            states.append(output)
        moment = Q(0)
        for labels in product((-1, 1), repeat=n):
            e = list(map(Q, labels))+[Q(0)]*(n*p)
            dual2 = dot(e, mv(proj, e))
            moment += dual2
            for output in states:
                assert dot(e, output)**2 <= dot(output, output)*dual2
                inequalities += 1
            assignments += 1
        assert moment/2**n == value_trace
        cases.append({'keys': n, 'dimension': p, 'key_rank': key_rank,
                      'quadratic_jet_rank': rank0, 'affine_jet_rank': rank,
                      'quadratic_value_trace': str(value_trace0),
                      'affine_value_trace': str(value_trace)})
    return {'cases': cases, 'independent_sign_moments': assignments,
            'projected_Cauchy_controls': inequalities}


def affine_sphere_difference_controls():
    rng, count = random.Random(1729), 0
    angles = [(Q(15, 17), Q(8, 17)), (Q(99, 101), Q(20, 101)),
              (Q(255, 257), Q(32, 257))]
    for u, (co, si), radial, bscale in product(unit_circle()[:6], angles,
                                              (Q(-2**64), Q(0), Q(2**64)),
                                              (Q(0), Q(1), Q(2**64))):
        delta, w = 2*si, [-u[1], u[0]]
        assert co*co+si*si == 1 and 0 < delta <= 1
        base = symmetric(2, rng)
        s = [[base[i][j]-radial*(i == j) for j in range(2)] for i in range(2)]
        b = [Q(3, 5)*bscale, Q(4, 5)*bscale]
        center = evaluate(s, b, radial, u)
        raw, odd = [], []
        for c, z in ((co, si), (2*co*co-1, 2*co*si)):
            qs = [[c*x+side*z*y for x, y in zip(u, w)] for side in (1, -1)]
            assert all(dot(q, q) == 1 and sum((x-y)**2 for x, y in zip(q, u)) <= delta**2 for q in qs)
            values = [evaluate(s, b, radial, q) for q in qs]
            raw.extend(values)
            odd.append((values[0]-values[1])/(2*z))
        weight = 2*(1+co)/(1+2*co)
        derivative = dot([2*x+y for x, y in zip(mv(s, u), b)], w)
        assert derivative == weight*odd[0]-(weight-1)*odd[1]
        assert 2*weight+(weight-1)/co == 2+1/co <= 4
        epsilon = max(abs(x-center) for x in raw)
        assert derivative**2 <= (4*epsilon/delta)**2
        count += 1
    states = keys = 0
    a, epsilon, delta = Q(1, 2), Q(1, 8), Q(1, 32)
    for p in range(2, 7):
        for labels, c in product(product((-1, 1), repeat=p), (Q(-2**64), Q(0), Q(2**64))):
            b = [Q(j+1, 1000) for j in range(p)]
            s = [[(a*labels[i]-b[i]-c)*(i == j) for j in range(p)] for i in range(p)]
            # On the sphere the radial coefficient cancels the bias exactly.
            norm = max(abs(a*labels[i]-b[i]) for i in range(p))
            assert 2*norm*delta+norm*delta**2+sum(abs(x) for x in b)*delta < epsilon
            points = [[Q(int(i == j)) for j in range(p)] for i in range(p)]
            output = value_and_half_tangent_gradient(points, s, b, c)
            assert output[:p] == [a*y for y in labels]
            ell2 = (a+epsilon)**2+4*epsilon**2/delta**2
            assert dot(output, output) <= p*ell2
            assert p <= 2*p*ell2/(a-epsilon)**2
            states, keys = states+1, keys+p
    return {'two_angle_tangent_derivative_identities': count,
            'whole_cap_affine_states': states, 'certified_affine_key_caps': keys}


def subspace_sum_counterexample():
    cases = []
    # Scalar second-coordinate measurement, Euclidean norm, two one-dimensional spaces.
    for eta in (Q(1, 2), Q(1, 16), Q(1, 2**32)):
        p0, _ = projector([[Q(1)], [Q(0)]])
        p1, _ = projector([[Q(1)], [eta]])
        full, _ = projector([[Q(1), Q(1)], [Q(0), eta]])
        k0, k1, combined = p0[1][1], p1[1][1], full[1][1]
        assert k0 == 0 and k1 == eta*eta/(1+eta*eta) and combined == 1
        assert combined > 2*(k0+k1)
        difference = plus(full, scaled(p0, Q(-1)))
        assert multiply(difference, difference) == difference and trace(difference) == 1
        cases.append({'eta': str(eta), 'first_subspace_score': str(k0),
                      'second_subspace_score': str(k1), 'sum_score': str(combined),
                      'twice_individual_score_sum': str(2*(k0+k1)),
                      'orthogonal_residual_rank': 1})
    return {'source': 'https://arxiv.org/html/2306.00945v2',
            'statement_checked': 'Lemma E.1(iii), also visually verified on PDF page 42',
            'scope': 'Counterexample to the unqualified subspace-sum inequality; no audit of consequences for the other results.',
            'cases': cases}


def main():
    sources = ['scripts/run_numerical_capacity_checks.py', 'scripts/run_retention_tradeoff_checks.py',
               'scripts/run_robust_reader_checks.py', 'scripts/run_robust_frontier_checks.py',
               'docs/interleaving_research.md']
    report = {'schema_version': 1, 'command': 'python3 -B scripts/run_numerical_capacity_checks.py',
              'scope': 'Uniform numerical recall with two-sided value error, not sign-only shattering or an architecture advantage. Exact identities support the written proof; novelty is unresolved.',
              'centered': centered_controls(), 'numerical': numerical_fixtures(),
              'sphere': spherical_controls(), 'clipping': clipping_controls(),
              'affine_sphere_projection': affine_sphere_projection_controls(),
              'affine_sphere_differences': affine_sphere_difference_controls(),
              'subspace_sum_counterexample': subspace_sum_counterexample(),
              'novelty_established': False, 'strongest_method_advantage_established': False,
              'source_sha256': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    (ROOT/'results/research/numerical_capacity_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k in (
        'centered', 'numerical', 'sphere', 'clipping', 'affine_sphere_projection',
        'affine_sphere_differences')}, sort_keys=True))


if __name__ == '__main__':
    main()
