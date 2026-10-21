#!/usr/bin/env python3
"""Exact controls for the minimum-colength primitive wall spectrum.

No finite test enumerates all ideals. The uniform necessity and rigidity
arguments are in docs/primitive_wall_spectrum.md.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd, isqrt, prod
from pathlib import Path
import hashlib
import json

from run_higher_interaction_checks import rank_trace
from run_structural_complexity_checks import cross, dot, finite_staircases, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]


def canonical(q):
    g = gcd(*q)
    assert g
    q = tuple(x // g for x in q)
    return q if next(x for x in q if x) > 0 else tuple(-x for x in q)


def threshold(q):
    return prod(max(x, 0) + 1 for x in q) + prod(max(-x, 0) + 1 for x in q) - 2


def divisors(a):
    return set(product(*(range(x + 1) for x in a)))


def mobius_sieve(n):
    mu = [1] * (n + 1)
    primes = []
    least = [0] * (n + 1)
    mu[0] = 0
    for a in range(2, n + 1):
        if not least[a]:
            least[a] = a
            primes.append(a)
            mu[a] = -1
        for p in primes:
            if a * p > n:
                break
            least[a * p] = p
            if a % p == 0:
                mu[a * p] = 0
                break
            mu[a * p] = -mu[a]
    return mu


def divisible_triples(n, h):
    M = n // h
    value = 0
    for B in range(1, max(0, (M - 2) // (h + 1)) + 1):
        C = max(0, (M - B - 1) // (h * B + 1))
        value += C * (M - B) - (h * B + 1) * C * (C + 1) // 2
    return value


def spectrum_count(n):
    mu = mobius_sieve(isqrt(n))
    full = sum(mu[h] * divisible_triples(n, h) for h in range(1, len(mu)))
    phi = list(range(n + 1))
    for p in range(2, n + 1):
        if phi[p] == p:
            for a in range(p, n + 1, p):
                phi[a] -= phi[a] // p
    pair = sum(phi[2:])
    return {'length': n, 'full_support_per_singleton_sign': full,
            'two_support_per_coordinate_pair': pair, 'unoriented_directions': 3 * (full + pair)}


@lru_cache(None)
def progression_product_histogram(k, h, N):
    data = [0] * (N + 1)
    data[1] = 1
    for _ in range(k):
        following = [0] * (N + 1)
        for a, count in enumerate(data):
            if count:
                for b in range(h + 1, N // a + 1, h):
                    following[a * b] += count
        data = following
    return tuple(data)


def fixed_dimension_count(n, dimension):
    N = n + 2
    mu = mobius_sieve(n)
    oriented = 0
    for support in range(2, dimension + 1):
        for positive in range(1, support):
            negative = support - positive
            primitive = 0
            for h in range(1, n + 1):
                if (h + 1)**positive + (h + 1)**negative > N:
                    break
                if not mu[h]:
                    continue
                left = progression_product_histogram(positive, h, N)
                right = progression_product_histogram(negative, h, N)
                prefix, total = [], 0
                for count in right:
                    total += count
                    prefix.append(total)
                value = sum(count * prefix[N - a] for a, count in enumerate(left))
                primitive += mu[h] * value
            oriented += comb(dimension, support) * comb(support, positive) * primitive
    assert oriented % 2 == 0
    return oriented // 2


def fixed_dimension_controls(counts):
    fixtures = []
    for dimension in range(2, 6):
        limit = 8 if dimension <= 4 else 5
        for n in range(1, limit + 1):
            direct = {canonical(q) for q in product(range(1 - n, n), repeat=dimension)
                      if min(q) < 0 < max(q) and gcd(*q) == 1 and threshold(q) <= n}
            value = fixed_dimension_count(n, dimension)
            assert value == len(direct)
            counts['independent_fixed_dimension_direction_cubes'] += 1
            fixtures.append({'dimension': dimension, 'length': n,
                             'unoriented_directions': value})
    for n in range(1, 129):
        assert fixed_dimension_count(n, 3) == spectrum_count(n)['unoriented_directions']
        counts['independent_progression_and_three_variable_formula_equalities'] += 1
    # Exact rational coefficients, rather than numerical zeta estimates or a fit.
    leading = [{'dimension': d,
                'leading_constant_times_zeta_d': str(Fraction(comb(2*d-2, d-1),
                                                             4*prod(range(1, d-1))))}
               for d in range(2, 7)]
    return {'complete_bounded_cube_counts': fixtures, 'leading_rational_factors': leading}


def counting_controls(counts):
    fixtures = []
    for n in range(1, 65):
        direct = set()
        for a in range(1, n):
            for b in range(1, n - a + 1):
                if a + b <= n and gcd(a, b) == 1:
                    for i, j in combinations(range(3), 2):
                        q = [0, 0, 0]
                        q[i], q[j] = a, -b
                        direct.add(tuple(q))
            for b in range(1, n):
                for c in range(1, max(0, (n - a - b) // (b + 1)) + 1):
                    assert a + b * c + b + c <= n
                    if gcd(a, b, c) == 1:
                        for j in range(3):
                            q = [-b, -c]
                            q.insert(j, a)
                            direct.add(canonical(q))
        result = spectrum_count(n)
        assert len(direct) == result['unoriented_directions']
        assert all(threshold(q) <= n for q in direct)
        for h in range(1, isqrt(n) + 1):
            brute = sum(1 for a in range(h, n + 1, h)
                        for b in range(h, n + 1, h)
                        for c in range(h, max(0, (n - a - b) // (b + 1)) + 1, h))
            assert divisible_triples(n, h) == brute
            counts['independent_divisible_triple_equalities'] += 1
        counts['independent_complete_small_direction_counts'] += 1
        if n <= 12 or n in (16, 32, 64):
            fixtures.append(result)
    larger = [spectrum_count(n) for n in (128, 256, 1024, 4096, 16384, 65536)]
    counts['exact_larger_direction_counts'] += len(larger)
    return {'small_exact_counts': fixtures, 'larger_exact_counts': larger}


def model(q, padding):
    u = tuple(max(x, 0) for x in q)
    v = tuple(max(-x, 0) for x in q)
    D = divisors(u) | divisors(v)
    forbidden = {next(i for i, x in enumerate(a) if x) for a in (u, v)
                 if sum(x > 0 for x in a) == 1}
    axis = next(i for i in range(len(q)) if i not in forbidden)
    old = max(a[axis] for a in D if sum(x > 0 for x in a) <= 1)
    for j in range(1, padding + 1):
        D.add(tuple(old + j if i == axis else 0 for i in range(len(q))))
    representative = min(u, v)
    basis = tuple(sorted(D - {max(u, v)}))
    indices = {a: j for j, a in enumerate(basis)}

    def target(a):
        if a not in D:
            return None
        return representative if a in (u, v) else a

    def ev(a):
        b = target(a)
        return {} if b is None else {indices[b]: 1}

    maps = tuple(tuple(None if (b := target(tuple(x + (i == j) for j, x in enumerate(a)))) is None
                       else indices[b] for a in basis) for i in range(len(q)))
    return u, v, D, basis, ev, maps, axis


def realization_controls(counts):
    directions = {canonical(q) for q in product(range(-3, 4), repeat=3)
                  if min(q) < 0 < max(q) and gcd(*q) == 1}
    # Higher-dimensional cases exercise the same proof and the padding-axis rule.
    directions.update({(1, 1, -1, -1), (2, -1, 0, 3), (1, -2, -1, 1),
                       (2, 1, -1, -2, 0), (1, 0, -1, 0, 0)})
    fixtures = []
    for q in sorted(directions, key=lambda q: (len(q), q)):
        for padding in (0, 3):
            u, v, D, basis, ev, maps, axis = model(q, padding)
            n = len(basis)
            assert n == threshold(q) + padding
            assert all(tuple(x - (i == j) for j, x in enumerate(a)) in D
                       for a in D for i, x in enumerate(a) if x)
            # Multiplication commutes on every quotient basis vector, including
            # the identified socle class. Every basis vector is a monomial in 1.
            for i, j in combinations(range(len(q)), 2):
                for a in range(n):
                    left = None if maps[i][a] is None else maps[j][maps[i][a]]
                    right = None if maps[j][a] is None else maps[i][maps[j][a]]
                    assert left == right
                    counts['commuting_integer_multiplication_replays'] += 1
            one = basis.index((0,) * len(q))
            for a in D:
                index = one
                for i, e in enumerate(a):
                    for _ in range(e):
                        index = None if index is None else maps[i][index]
                assert index is not None and ev(a) == {index: 1}
                counts['cyclic_monomial_column_replays'] += 1
            P, N = sum(u), sum(v)
            tie = tuple(N if x > 0 else P if x < 0 else 1 for x in q)
            assert min(tie) > 0 and dot(tie, q) == 0
            pivot = next(i for i, x in enumerate(q) if x)
            weights = [tuple(2 * x + sign * (i == pivot) for i, x in enumerate(tie))
                       for sign in (-1, 1)]
            records = []
            for w in weights:
                assert min(w) > 0 and dot(w, q)
                omitted = u if dot(w, q) > 0 else v
                expected = D - {omitted}
                assert all(tuple(x - (i == j) for j, x in enumerate(a)) in expected
                           for a in expected for i, x in enumerate(a) if x)
                # The finite ground is independently generated from its divisor
                # product budget, rather than using D as the input ground.
                E = tuple(hyperbolic_cross(n, len(q)))
                for p in (0, 2, 3, 5):
                    selected, _ = rank_trace(sorted(E, key=lambda a: (dot(w, a), a)), ev, p, set())
                    assert set(selected) == expected
                    counts['independent_field_greedy_standard_sets'] += 1
                profile = tuple(sum(a[i] for a in expected) for i in range(len(q)))
                records.append({'strict_weight': w, 'profile': profile, 'omitted': omitted})
            assert canonical(tuple(a - b for a, b in zip(records[0]['profile'], records[1]['profile']))) == q
            # Independently count grading classes of the compulsory divisor set.
            minimal = divisors(u) | divisors(v)
            classes = []
            def same_degree(a, b):
                delta = tuple(x - y for x, y in zip(a, b))
                return all(delta[i] * q[pivot] == delta[pivot] * q[i] for i in range(len(q)))
            for a in minimal:
                if not any(same_degree(a, b) for b in classes):
                    classes.append(a)
            assert len(classes) == threshold(q)
            assert all(sum(same_degree(a, b) for b in minimal) == (2 if a in (u, v) else 1)
                       for a in minimal)
            counts['exact_compulsory_grading_footprints'] += 1
            if q in ((1, -1, -1), (2, -1, -1), (3, -2, -3), (1, 1, -1, -1)):
                fixtures.append({'direction': q, 'minimum_length': threshold(q), 'padding': padding,
                                 'length': n, 'padding_axis': axis,
                                 'compulsory_divisors': sorted(minimal), 'both_cones': records})
    return {'primitive_directions_constructed': len(directions), 'selected_exact_models': fixtures}


def saved_fan_controls(counts):
    report = json.loads((ROOT / 'results/research/oracle_enumeration_checks.json').read_text())
    assert report['status'] == 'passed'
    fixtures = []
    for case in report['direct_positive_cone_traversal']['fixtures']:
        n = case['length']
        directions = set()
        for row in case['directed_neighbor_trace']:
            a, b = [tuple(Fraction(x) for x in p) for p in row['shared_edge']]
            normal = cross(a, b)
            scale = 1
            for x in normal:
                scale = scale * x.denominator // gcd(scale, x.denominator)
            q = canonical(tuple(int(x * scale) for x in normal))
            delta = tuple(x - y for x, y in zip(row['to_profile'], row['from_profile']))
            assert canonical(delta) == q
            assert min(q) < 0 < max(q) and threshold(q) <= n
            assert dot(q, a) == dot(q, b) == 0
            midpoint = tuple((x + y) / 2 for x, y in zip(a, b))
            assert min(midpoint) > 0
            directions.add(q)
            counts['independent_saved_positive_adjacency_budget_checks'] += 1
        fixtures.append({'model': case['model'], 'length': n, 'characteristic': case['characteristic'],
                         'wall_directions': sorted(directions),
                         'minimum_colengths': [threshold(q) for q in sorted(directions)]})
    return fixtures


def independent_staircase_controls(counts):
    """Test the lower-bound and equality lemmas without quotient matrices.

    Equal quotient-grading Hilbert functions are necessary for an edge;
    this does not assume that every such pair is an actual edge.
    """
    fixtures = []
    for n in range(1, 9):
        states = finite_staircases(n, 3)
        moments = {B: tuple(sum(a[i] for a in B) for i in range(3)) for B in states}
        accepted, rigid = 0, 0
        for B, C in combinations(states, 2):
            delta = tuple(x - y for x, y in zip(moments[B], moments[C]))
            if not any(delta):
                continue
            q = canonical(delta)
            if not min(q) < 0 < max(q):
                continue
            pivot = next(i for i, x in enumerate(q) if x)

            def histogram(state):
                out = defaultdict(int)
                for a in state:
                    degree = tuple(q[pivot] * a[i] - q[i] * a[pivot]
                                   for i in range(3) if i != pivot)
                    out[degree] += 1
                return dict(out)

            counts['independent_downset_pair_grading_tests'] += 1
            if histogram(B) != histogram(C):
                continue
            accepted += 1
            assert n >= threshold(q)
            counts['independent_equal_hilbert_function_thresholds'] += 1
            if n == threshold(q):
                u = tuple(max(x, 0) for x in q)
                v = tuple(max(-x, 0) for x in q)
                D = divisors(u) | divisors(v)
                assert {B, C} == {frozenset(D - {u}), frozenset(D - {v})}
                rigid += 1
                counts['independent_minimum_length_rigidity_pairs'] += 1
        fixtures.append({'length': n, 'all_monomial_downsets': len(states),
                         'equal_grading_hilbert_pairs': accepted, 'minimum_rigid_pairs': rigid})
    return fixtures


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    counts = defaultdict(int)
    counting = counting_controls(counts)
    dimensions = fixed_dimension_controls(counts)
    models = realization_controls(counts)
    saved = saved_fan_controls(counts)
    staircases = independent_staircase_controls(counts)
    sources = ('docs/primitive_wall_spectrum.md', 'scripts/run_primitive_wall_checks.py',
               'scripts/run_higher_interaction_checks.py', 'scripts/run_structural_complexity_checks.py',
               'results/research/oracle_enumeration_checks.json')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_primitive_wall_checks.py',
              'counts': dict(counts), 'direction_counting': counting, 'actual_quotients': models,
              'fixed_dimension_counting': dimensions,
              'independent_saved_fans': saved,
              'independent_monomial_downset_hilbert_controls': staircases,
              'uniform_results': {
                  'minimum_colength': 'prod(q_plus+1)+prod(q_minus+1)-2, for primitive mixed q',
                  'minimum_rigidity': 'The wall degeneration is J_D+(x^q_plus-lambda*x^q_minus), lambda nonzero; D is the union of the two divisor boxes.',
                  'padding': 'Every larger exact length in dimension at least three, over every field.',
                  'common_fan': 'The exact common positive refinement of all length-n fans is the arrangement of all primitive mixed q with the minimum colength at most n.',
                  'three_variable_spectrum': 'D3(n) ~ 3*n^2*log(n)/(2*zeta(3)).',
                  'fixed_dimension_spectrum': 'D_d(n) ~ binom(2*d-2,d-1)*n^2*log(n)^(d-2)/(4*(d-2)!*zeta(d)); exact length for d>=3, length at most n for d=2.',
                  'universal_shadow_envelope_sum': 'Theta(n^(10/3)*log(n+1)^2), a bound on the stated envelope budget, not actual shadow complexity.'},
              'new_unrestricted_G3_count': False, 'polynomial_lower_separation_proved': False,
              'novelty_established': False, 'goal_complete': False,
              'scope': 'The necessity and rigidity arguments are uniform proofs in the note. Finite controls replay explicit all-field cyclic commuting quotient models, exact primitive counts, and independently extracted directions from saved fan geometry. They do not enumerate all ideals, prove a collective incidence budget, or establish priority.',
              'source_sha256': {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in sources}}
    (ROOT / 'results/research/primitive_wall_checks.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'counts': report['counts'],
                      'selected_counts': counting['small_exact_counts'][-3:]}, indent=2))


if __name__ == '__main__':
    main()
