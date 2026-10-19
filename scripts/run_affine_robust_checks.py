#!/usr/bin/env python3
"""Exact controls for affine-quadratic robust label-capacity proofs.

The uniform proofs, domain distinctions and priority limits are in the
research note. No search optimizer, training, timing or novelty claim.
"""
from fractions import Fraction as Q
from itertools import product
from hashlib import sha256
import json
from pathlib import Path

from run_robust_reader_checks import unit_circle

ROOT = Path(__file__).resolve().parents[1]
SPHERE_GAP_CONSTANT = 64*660**2


def dot(x, y):
    return sum((a*b for a, b in zip(x, y)), Q(0))


def mv(s, x):
    return [dot(row, x) for row in s]


def value(s, b, c, x):
    return dot(x, mv(s, x))+dot(b, x)+c


def cap_measure(u, beta, sh, ch):
    """An explicit rational measure on a two-dimensional spherical cap."""
    delta, cosine, sine = 2*sh, 1-2*sh*sh, 2*sh*ch
    alpha = delta*delta/4
    assert sh*sh+ch*ch == 1 and 0 < delta < 1
    ub = dot(u, beta)
    if -ub >= cosine:
        measure = [(1-alpha, u), (alpha, [-x for x in beta])]
    else:
        w = [-u[1], u[0]]
        mass = (1+ub)/2
        imbalance = -alpha*dot(beta, w)/sine
        measure = [(1-mass, u), ((mass+imbalance)/2,
                   [cosine*x+sine*y for x, y in zip(u, w)]),
                   ((mass-imbalance)/2,
                   [cosine*x-sine*y for x, y in zip(u, w)])]
    assert all(weight >= 0 and dot(q, q) == 1 and
               sum((x-y)**2 for x, y in zip(q, u)) <= delta*delta
               for weight, q in measure)
    assert sum(weight for weight, _ in measure) == 1
    mean = [sum(weight*q[j] for weight, q in measure) for j in range(2)]
    assert mean == [(1-alpha)*x-alpha*y for x, y in zip(u, beta)]
    assert sum(weight*sum((x-y)**2 for x, y in zip(q, u))
               for weight, q in measure) == 2*alpha*(1+ub)
    return alpha, measure


def moment_controls():
    points, checks, measures, derivatives = unit_circle()[::2], 0, 0, 0
    angles = [(Q(5, 13), Q(12, 13)), (Q(7, 25), Q(24, 25)),
              (Q(9, 41), Q(40, 41)), (Q(80, 1601), Q(1599, 1601))]
    for u, beta, (sh, ch) in product(points, points, angles):
        alpha, measure = cap_measure(u, beta, sh, ch)
        measures += 1
        for width, rotation, radial, magnitude in product(
                (Q(0), Q(1), Q(20)), points[:2], (Q(-2**64), Q(0), Q(2**64)),
                (Q(1), Q(100))):
            a, z = rotation
            # Orthogonally conjugate diag(-width/2,+width/2), then add radial*I.
            t = [[width*(z*z-a*a)/2, -width*a*z],
                 [-width*a*z, width*(a*a-z*z)/2]]
            s = [[t[i][j]+radial*(i == j) for j in range(2)] for i in range(2)]
            b = [magnitude*x for x in beta]
            change = sum(weight*value(s, b, Q(17), q) for weight, q in measure)-value(s, b, Q(17), u)
            assert abs(change+alpha*magnitude*(1+dot(u, beta))) <= 4*alpha*width
            checks += 1
            w = [-u[1], u[0]]
            second = 2*dot(w, mv(s, w))-2*dot(u, mv(s, u))-dot(b, u)
            fourth = 8*dot(u, mv(s, u))-8*dot(w, mv(s, w))+dot(b, u)
            assert abs(second) <= 2*(width+magnitude)
            assert abs(fourth) <= 8*(width+magnitude)
            derivatives += 1
    return {'cap_measures': measures, 'exact_centered_moment_controls': checks,
            'geodesic_derivative_bound_controls': derivatives}


def transfer_controls():
    """Check the sphere proof's residual identities and conservative constants."""
    points, count = unit_circle()[::2], 0
    for lo, hi, beta, u, v in product((Q(-3), Q(0)), (Q(1), Q(4)), points, points, points):
        if u == v:
            continue
        s, b = [[lo, Q(0)], [Q(0), hi]], [Q(3)*x for x in beta]
        fu, fv = value(s, b, Q(0), u), value(s, b, Q(0), v)
        if fu < fv:
            u, v, fu, fv = v, u, fv, fu
        gap, d = fu-fv, 1-dot(u, v)
        lu, lv = dot(u, mv(s, u))+dot(b, u)/2, dot(v, mv(s, v))+dot(b, v)/2
        ru = [x+y/2-lu*z for x, y, z in zip(mv(s, u), b, u)]
        rv = [x+y/2-lv*z for x, y, z in zip(mv(s, v), b, v)]
        wu, wv = [-u[1], u[0]], [-v[1], v[0]]
        au = dot(wu, mv(s, wu))-lu
        av = dot(wv, mv(s, wv))-lv
        eta = max(gap, sum(abs(x) for x in ru), sum(abs(x) for x in rv),
                  -3*au/4, 3*av/4, Q(0))
        assert d*(lu-lv) == gap-dot(u, rv)+dot(v, ru)
        assert abs(lu-lv) <= 3*eta/d
        z = [x-y for x, y in zip(u, v)]
        az = [x-lu*y for x, y in zip(mv(s, z), z)]
        assert az == [x-y+(lu-lv)*w for x, y, w in zip(ru, rv, v)]
        assert max(abs(lo-lu), abs(hi-lu)) <= 161*eta/d**3
        assert hi-lo+3 <= 660*eta/d**3
        count += 1
    algebra = 0
    for delta, width, magnitude in product((Q(1, 10), Q(1, 4), Q(1, 2), Q(3, 4)),
                                          (Q(0), Q(1, 100), Q(1), Q(8), Q(9), Q(10), Q(100)),
                                          (Q(0), Q(1, 100), Q(1), Q(8), Q(9), Q(10), Q(100))):
        if width+magnitude == 0:
            continue
        r, alpha = delta*delta, delta*delta/4
        gap = max(r*width, alpha*(2*magnitude-9*width)/(1-alpha))+Q(1, 10**6)
        assert gap > r*(width+magnitude)/10
        algebra += 1
    assert SPHERE_GAP_CONSTANT == 27878400
    return {'residual_transfer_cases': count, 'full_ball_gap_algebra_cases': algebra,
            'sphere_gap_constant': SPHERE_GAP_CONSTANT}


def affine_fixtures():
    count = keys = 0
    delta, bnorm = Q(1, 10), Q(1, 20)
    for p in range(2, 7):
        b = [Q(3, 100), Q(4, 100)]+[Q(0)]*(p-2)
        assert dot(b, b) == bnorm*bnorm
        for labels in product((-1, 1), repeat=p):
            if len(set(labels)) == 1:
                continue
            for radial in (Q(-1, 10), Q(0), Q(1, 10)):
                s = [[(labels[i]+radial)*(i == j) for j in range(p)] for i in range(p)]
                c = Q(1, 100)-radial
                values = [labels[i]+b[i]+Q(1, 100) for i in range(p)]
                # Whole-ball certificate from ||S||<=11/10, ||b||=1/20.
                center_floor = 1-bnorm-Q(1, 100)
                variation = 2*Q(11, 10)*delta+Q(11, 10)*delta**2+bnorm*delta
                assert center_floor > variation
                lower = min(v for v, label in zip(values, labels) if label > 0)
                upper = max(v for v, label in zip(values, labels) if label < 0)
                tau, width = (lower+upper)/2, Q(2)
                total = width+bnorm
                assert lower-upper > delta**2*total/10
                constant = c+radial-tau
                weight2 = (p+dot(b, b)+constant*constant)/total**2
                assert weight2 <= Q(p, 4)+2
                for i, label in enumerate(labels):
                    assert label*(values[i]-tau)/total > delta**2/20
                count, keys = count+1, keys+p
    return {'mixed_label_whole_ball_instances': count, 'certified_key_balls': keys}


def lift_controls():
    frames = [unit_circle()[:6],
              [[Q(int(i == j)) for j in range(3)] for i in range(3)] +
              [[Q(3, 5), Q(4, 5), Q(0)], [Q(0), Q(3, 5), Q(4, 5)], [Q(4, 5), Q(0), Q(3, 5)]]]
    count = 0
    for points in frames:
        features = [[a*b for a in u for b in u]+list(u)+[Q(1)] for u in points]
        assert all(dot(x, x) == 3 for x in features)
        moment = Q(0)
        for labels in product((-1, 1), repeat=len(points)):
            x = [sum(label*f[j] for label, f in zip(labels, features)) for j in range(len(features[0]))]
            moment += dot(x, x)
            count += 1
        assert moment/Q(2**len(points)) == 3*len(points)
    # Unit-sphere robustness of this same form does not imply full-ball robustness.
    s, b = [[Q(1), Q(0)], [Q(0), Q(1)]], [Q(1, 100), Q(0)]
    assert value(s, b, Q(-1), [Q(1), Q(0)]) > 0
    assert value(s, b, Q(-1), [Q(9, 10), Q(0)]) == Q(-181, 1000)
    return {'affine_lift_sign_assignments': count,
            'sphere_only_versus_full_ball_counterexample': True}


def main():
    sources = ['scripts/run_affine_robust_checks.py', 'scripts/run_robust_reader_checks.py',
               'docs/interleaving_research.md']
    report = {'schema_version': 1, 'command': 'python3 -B scripts/run_affine_robust_checks.py',
              'scope': 'Exact supporting controls for full-ball and spherical affine-quadratic capacity proofs; no priority, runtime or architecture-advantage claim.',
              'moments': moment_controls(), 'transfer': transfer_controls(),
              'fixtures': affine_fixtures(), 'lift': lift_controls(),
              'new_theorem_beyond_prior_art_established': False,
              'strongest_method_advantage_established': False,
              'source_sha256': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    (ROOT/'results/research/affine_robust_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('moments', 'transfer', 'fixtures', 'lift')}, sort_keys=True))


if __name__ == '__main__':
    main()
