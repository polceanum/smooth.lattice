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


def main():
    sources = ['scripts/run_numerical_capacity_checks.py', 'scripts/run_retention_tradeoff_checks.py',
               'scripts/run_robust_reader_checks.py', 'scripts/run_robust_frontier_checks.py',
               'docs/interleaving_research.md']
    report = {'schema_version': 1, 'command': 'python3 -B scripts/run_numerical_capacity_checks.py',
              'scope': 'Uniform numerical recall with two-sided value error, not sign-only shattering or an architecture advantage. Exact identities support the written proof; novelty is unresolved.',
              'centered': centered_controls(), 'numerical': numerical_fixtures(),
              'sphere': spherical_controls(), 'clipping': clipping_controls(),
              'novelty_established': False, 'strongest_method_advantage_established': False,
              'source_sha256': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    (ROOT/'results/research/numerical_capacity_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k in ('centered', 'numerical', 'sphere', 'clipping')}, sort_keys=True))


if __name__ == '__main__':
    main()
