#!/usr/bin/env python3
"""Rational controls for robust quadratic recall and a quartic binding table.

Finite controls support the uniform proof in the research note. They are not
a shattering optimizer, a novelty test, or a neural/runtime benchmark.
"""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product, combinations
import json
from pathlib import Path
import random

ROOT = Path(__file__).resolve().parents[1]


def unit_circle():
    points = {(-Q(1), Q(0))}
    for t in (Q(-3), Q(-2), Q(-1), Q(-1, 2), Q(-1, 3), Q(0),
              Q(1, 3), Q(1, 2), Q(1), Q(2), Q(3)):
        points.add(((1-t*t)/(1+t*t), 2*t/(1+t*t)))
    assert all(x*x+y*y == 1 for x, y in points)
    return sorted(points)


def normalization_controls():
    checked = mixed_geometries = 0
    for a_root, b_root in product(range(1, 4), repeat=2):
        a, b = Q(a_root*a_root), Q(b_root*b_root)
        for delta in (Q(1, 10), Q(1, 5), Q(1, 3), Q(2, 5)):
            r = delta*delta
            valid = []
            for x, y in unit_circle():
                value = -a*x*x+b*y*y
                # The two zero lines have directions (sqrt(b), +/-sqrt(a)).
                # Exact distance to their union certifies the entire ball.
                distance2 = min((a_root*x-b_root*y)**2,
                                (a_root*x+b_root*y)**2)/(a+b)
                if value and distance2 > r:
                    valid.append((x, y, value))
            if not any(v > 0 for _, _, v in valid) or not any(v < 0 for _, _, v in valid):
                continue
            mixed_geometries += 1
            assert a > b*r and b > a*r
            gamma, bound = r/2, 1-r/2
            shift = (a-b)*r/(2*(a+b))
            eigenvalues = [-a/(a+b)-shift, b/(a+b)-shift]
            assert max(abs(v) for v in eigenvalues) <= bound
            for x, y, value in valid:
                if value > 0:
                    assert value-2*a*delta*abs(x)-a*r > 0
                    assert value/(a+b) > a*r/(a+b)
                else:
                    assert value+2*b*delta*abs(y)+b*r < 0
                    assert -value/(a+b) > b*r/(a+b)
                for scale in (Q(1), Q(2**80), Q(1, 2**60)):
                    # Scaling the unbounded reader changes no normalized result.
                    normalized = (scale*value)/(scale*(a+b))-shift
                    assert (1 if value > 0 else -1)*normalized > gamma
                    checked += 1
    return {'mixed_two_eigenvalue_geometries': mixed_geometries,
            'whole_ball_normalization_checks_including_rescalings': checked}


def rademacher_controls():
    frames = [unit_circle()[:8],
              [[Q(int(j == i)) for j in range(3)] for i in range(3)] +
              [[Q(3, 5), Q(4, 5), Q(0)], [Q(0), Q(3, 5), Q(4, 5)]],
              [[Q(int(j == i)) for j in range(4)] for i in range(4)] +
              [[Q(1, 5), Q(2, 5), Q(2, 5), Q(4, 5)]]]
    assignments = 0
    for frame in frames:
        n, p = len(frame), len(frame[0])
        assert all(sum(x*x for x in u) == 1 for u in frame)
        expected = Q(0)
        for signs in product((-1, 1), repeat=n):
            x = [[sum(s*u[i]*u[j] for s, u in zip(signs, frame))
                  for j in range(p)] for i in range(p)]
            expected += sum(v*v for row in x for v in row)
            assignments += 1
        assert expected/Q(2**n) == n
    return {'frames': len(frames), 'sign_assignments': assignments,
            'mean_squared_frobenius_norm_equals_key_count': True}


def quartic_read(q, values):
    masses = {(i, j): q[i]**2*q[j]**2 for i, j in values}
    denominator = sum(masses.values(), Q(0))
    assert denominator > 0
    return sum(masses[pair]*v for pair, v in values.items())/denominator


def pair_key(p, pair):
    q = [Q(0)]*p
    q[pair[0]], q[pair[1]] = Q(3, 5), Q(4, 5)
    return q


def top_two(q):
    # Sorting is only a proof oracle. The largest two can be found in O(p).
    return tuple(sorted(sorted(range(len(q)), key=lambda j: abs(q[j]), reverse=True)[:2]))


def binding_controls():
    delta = Q(1, 10)
    target_floor = (Q(3, 5)-delta)**2*(Q(4, 5)-delta)**2
    other_ceiling = (1+delta)**2*delta**2+delta**4/2
    bound = 2*other_ceiling/target_floor
    assert bound == Q(243, 1225) and bound < Q(1, 5)
    assignments = reads = 0
    for p in (2, 3, 4):
        pairs = list(combinations(range(p), 2))
        for labels in product((-Q(1), Q(1)), repeat=len(pairs)):
            values = dict(zip(pairs, labels))
            assignments += 1
            for pair in pairs:
                key = pair_key(p, pair)
                probes = [key]
                for coordinate, sign in product(range(p), (-1, 1)):
                    q = key[:]
                    q[coordinate] += sign*delta
                    probes.append(q)
                for sign in (-1, 1):
                    probes.append([x+sign*delta/p for x in key])
                for q in probes:
                    assert sum((x-y)**2 for x, y in zip(q, key)) <= delta**2
                    assert top_two(q) == pair
                    assert abs(quartic_read(q, values)-values[pair]) <= bound
                    reads += 1
    p, rng = 8, random.Random(9031)
    values = {pair: Q(0) for pair in combinations(range(p), 2)}
    for step in range(96):
        pair = rng.choice(list(values))
        before = values.copy()
        values[pair] = Q(0) if step % 4 == 0 else Q(rng.randrange(-4, 5), 4)
        assert all(values[k] == before[k] for k in values if k != pair)
        assert all(quartic_read(pair_key(p, k), values) == v for k, v in values.items())
        # A full-dimensional fresh query, not reduced to two latent scalars.
        query = [Q(1, j+1) for j in range(p)]
        answer = quartic_read(query, values)
        assert -1 <= answer <= 1
    return {'exhaustive_label_assignments': assignments, 'noisy_query_controls': reads,
            'point_replacements_including_zero_clears': 96,
            'fixed_radius': str(delta), 'uniform_value_error_bound': str(bound),
            'indexed_top_two_baseline_matches_all_tested_key_neighborhoods': True}


def main():
    sources = ['scripts/run_robust_reader_checks.py', 'docs/interleaving_research.md']
    report = {'schema_version': 1,
              'command': 'python3 -B scripts/run_robust_reader_checks.py',
              'scope': 'Exact finite controls for the uniform cap-capacity proof; no novelty, runtime or strongest-method claim.',
              'quadratic_capacity_bound': 'n <= p*(2/delta^2-1)^2 for fixed unit keys robustly realizing every binary assignment on full Euclidean radius-delta balls.',
              'normalization': normalization_controls(), 'rademacher': rademacher_controls(),
              'binding': binding_controls(),
              'new_theorem_beyond_prior_art_established': False,
              'strongest_method_advantage_established': False,
              'source_sha256': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    (ROOT/'results/research/robust_reader_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('normalization', 'rademacher', 'binding')}, sort_keys=True))


if __name__ == '__main__':
    main()
