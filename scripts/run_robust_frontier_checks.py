#!/usr/bin/env python3
"""Exact capacity-frontier and literature-quantifier controls.

These are standard binding tables and explicit counterfamilies, not a new
architecture, a neural experiment, a sharp upper bound or a novelty test.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from hashlib import sha256
import json
from pathlib import Path
import random

ROOT = Path(__file__).resolve().parents[1]
SCALE = Q(25, 24)
ERROR = Q(195, 512)


def dot(x, y):
    return sum((a*b for a, b in zip(x, y)), Q(0))


def quad(s, q):
    support = [(i, x) for i, x in enumerate(q) if x]
    return sum((x*y*s[i][j] for i, x in support for j, y in support), Q(0))


def clip(x):
    return max(Q(-1), min(Q(1), x))


def pairs(p, r):
    return [pair for start in range(0, p, r)
            for pair in combinations(range(start, min(start+r, p)), 2)]


def matrix(p, values):
    s = [[Q(0)]*p for _ in range(p)]
    for (i, j), v in values.items():
        s[i][j] = s[j][i] = SCALE*v
    return s


def key(p, pair):
    q = [Q(0)]*p
    q[pair[0]], q[pair[1]] = Q(3, 5), Q(4, 5)
    return q


def direct(values, q):
    return 2*SCALE*sum((v*q[i]*q[j] for (i, j), v in values.items()), Q(0))


def block_controls():
    rng = random.Random(12031)
    configurations = [(4, 2, Q(1, 8)), (6, 3, Q(1, 16)), (4, 4, Q(1, 24)),
                      (10, 5, Q(1, 16)), (20, 10, Q(1, 24))]
    states = exact_keys = noisy = 0
    for p, r, delta in configurations:
        addresses = pairs(p, r)
        assert delta**2*(r-1) <= Q(1, 64)
        labelings = (list(product((-1, 1), repeat=len(addresses))) if len(addresses) <= 6 else
                     [tuple([1]*len(addresses)), tuple([-1]*len(addresses))] +
                     [tuple(rng.choice((-1, 1)) for _ in addresses) for _ in range(6)])
        for labels in labelings:
            values = dict(zip(addresses, map(Q, labels)))
            s = matrix(p, values)
            assert max(sum(abs(x) for x in row) for row in s) <= SCALE*(r-1)
            for pair in addresses:
                u = key(p, pair)
                assert dot(u, u) == 1 and quad(s, u) == values[pair]
                su = [dot(row, u) for row in s]
                assert dot(su, su) <= Q(35, 24)**2*(r-1)
                exact_keys += 1
            for pair in {addresses[0], addresses[len(addresses)//2], addresses[-1]}:
                u = key(p, pair)
                probes = []
                for coordinate, sign in product(range(p), (-1, 1)):
                    q = u[:]
                    q[coordinate] += sign*delta
                    probes.append(q)
                probes += [[x+sign*delta/p for x in u] for sign in (-1, 1)]
                for q in probes:
                    assert sum((x-y)**2 for x, y in zip(q, u)) <= delta**2
                    value = quad(s, q)
                    assert value == direct(values, q)
                    assert abs(value-values[pair]) <= ERROR
                    assert abs(clip(value)-values[pair]) <= ERROR
                    top = tuple(sorted(sorted(range(p), key=lambda j: abs(q[j]), reverse=True)[:2]))
                    assert top == pair
                    noisy += 1
            states += 1
    assert Q(35, 12)/8+Q(25, 24)/64 == ERROR < Q(1, 2)
    count_controls = 0
    for p, delta in product((2, 3, 4, 7, 16, 64, 257),
                            (Q(1, 8), Q(1, 16), Q(1, 32), Q(1, 64))):
        r = min(p, 1+int(1/(64*delta**2)))
        count = len(pairs(p, r))
        assert count >= Q(p*(r-1), 4)
        assert count <= Q(p*(r-1), 2)
        count_controls += 1
    p, r = 12, 4
    values = {pair: Q(0) for pair in pairs(p, r)}
    for step in range(96):
        chosen = rng.choice(list(values))
        old = matrix(p, values)
        values[chosen] = Q(0) if step % 4 == 0 else Q(rng.randrange(-4, 5), 4)
        new = matrix(p, values)
        assert all(old[i][j] == new[i][j] for i in range(p) for j in range(p)
                   if (i, j) not in (chosen, chosen[::-1]))
        assert all(quad(new, key(p, pair)) == v for pair, v in values.items())
        q = [Q(j+1, p*p) for j in range(p)]
        assert quad(new, q) == direct(values, q)
    return {'configurations': len(configurations), 'value_assignments': states,
            'exact_key_reads_and_gradient_norm_checks': exact_keys,
            'noisy_queries_and_indexed_baseline_checks': noisy,
            'capacity_count_controls': count_controls, 'point_overwrites_and_zero_clears': 96,
            'uniform_error_bound': str(ERROR)}


def opposed_cap_controls():
    cases = []
    for k in range(1, 9):
        t = Q(1, 2**k)
        a, s = (1-t*t)/(1+t*t), 2*t/(1+t*t)
        d, delta = 1-a, t/4
        gap, size = d*d, 1+2*a
        positive_floor = (d-delta**2/2)**2-d*d/4
        negative_ceiling = (delta**2/2+s*delta)**2-d*d/4
        assert positive_floor >= Q(705, 1024)*d*d > 0
        assert negative_ceiling <= -Q(175, 1024)*d*d < 0
        ratio = gap/(size*delta**4)
        assert ratio == 1024/((1+t*t)*(3-t*t))
        assert 256 <= ratio <= Q(1024, 3)
        u, v = [Q(1), Q(0), Q(0)], [a, s, Q(0)]
        small = t/32
        c_rot, s_rot = (1-small*small)/(1+small*small), 2*small/(1+small*small)
        for center, directions, sign in ((u, [[Q(0), Q(1), Q(0)], [Q(0), Q(0), Q(1)]], 1),
                                         (v, [[-s, a, Q(0)], [Q(0), Q(0), Q(1)]], -1)):
            for direction, side in product(directions, (-1, 1)):
                q = [c_rot*x+side*s_rot*y for x, y in zip(center, direction)]
                assert dot(q, q) == 1
                assert sum((x-y)**2 for x, y in zip(q, center)) <= delta**2
                assert sign*((q[0]-a)**2-d*d/4) > 0
        cases.append({'t': str(t), 'radius': str(delta),
                      'gap_over_coefficient_size_and_radius_fourth': str(ratio),
                      'gap_over_coefficient_size_and_radius_squared': str(gap/(size*delta**2))})
    return cases


def constant_controls():
    rows = []
    for bits in (1, 4, 16, 64):
        epsilon = Q(1, 2**bits)
        # f=epsilon*x_1^2+x_2^2 on the ball; minimizer zero, test x=e_1.
        assert epsilon <= Q(1, 2)
        for a, b in ((Q(4, 5), Q(3, 5)), (Q(99, 101), Q(20, 101)),
                     (Q(9999, 10001), Q(200, 10001))):
            assert a*a+b*b == 1
            f = epsilon*b*b
            gradient2 = 4*epsilon**2*a*a*b*b
            ratio4 = f**3/gradient2**2
            assert ratio4 == b*b/(16*epsilon*a**4)
            rows.append({'epsilon': str(epsilon), 'sphere_point': [str(a), str(b), '0'],
                         'necessary_spherical_KL_constant_fourth_power': str(ratio4),
                         'necessary_half_error_bound_constant_squared': str(1/epsilon),
                         'necessary_quarter_error_bound_constant_fourth_power': str(1/epsilon)})
    return rows


def main():
    sources = ['scripts/run_robust_frontier_checks.py', 'docs/interleaving_research.md']
    report = {'schema_version': 1,
              'command': 'python3 -B scripts/run_robust_frontier_checks.py',
              'scope': 'Exact standard block-table lower bounds and counterfamilies for uniform constants/gap powers; no new architecture, sharp capacity or novelty claim.',
              'blocks': block_controls(), 'opposed_spherical_caps': opposed_cap_controls(),
              'problem_dependent_constants': constant_controls(),
              'new_theorem_beyond_prior_art_established': False,
              'strongest_method_advantage_established': False,
              'source_sha256': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    (ROOT/'results/research/robust_frontier_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'blocks': report['blocks'], 'opposed_cap_families': len(report['opposed_spherical_caps']),
                      'constant_controls': len(report['problem_dependent_constants'])}, sort_keys=True))


if __name__ == '__main__':
    main()
