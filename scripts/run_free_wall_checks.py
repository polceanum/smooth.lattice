#!/usr/bin/env python3
"""Exact free corner-cut witnesses versus algebraic wall thresholds.

The uniform proof is in docs/free_wall_threshold.md. Finite witnesses do
not establish novelty or a new count of coexisting initial ideals.
"""
from collections import defaultdict
from fractions import Fraction
from heapq import heappop, heappush
from itertools import product
from math import gcd, lcm, prod
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def canonical(q):
    assert gcd(*q) == 1 and min(q) < 0 < max(q)
    return q if next(x for x in q if x) > 0 else tuple(-x for x in q)


def first_terms(w, size):
    """Independent heap enumeration on the entire free nonnegative lattice."""
    zero = (0,)*len(w)
    queue = [(0, zero)]
    seen = {zero}
    selected = []
    for _ in range(size):
        _, a = heappop(queue)
        selected.append(a)
        for i in range(len(w)):
            b = tuple(x+(i == j) for j, x in enumerate(a))
            if b not in seen:
                seen.add(b)
                heappush(queue, (dot(w, b), b))
    return frozenset(selected)


def socle_model(q, counts):
    u = tuple(max(x, 0) for x in q)
    v = tuple(max(-x, 0) for x in q)
    D = set(product(*(range(x+1) for x in u))) | set(product(*(range(x+1) for x in v)))
    n = prod(x+1 for x in u)+prod(x+1 for x in v)-2
    assert len(D) == n+1
    B = tuple(sorted(D-{u}))
    index = {a: i for i, a in enumerate(B)}

    def column(a):
        return index[v] if a == u else index.get(a)

    maps = [tuple(column(tuple(x+(i == j) for j, x in enumerate(a))) for a in B)
            for i in range(len(q))]
    for i in range(len(q)):
        for j in range(len(q)):
            for k in range(n):
                left = None if maps[j][k] is None else maps[i][maps[j][k]]
                right = None if maps[i][k] is None else maps[j][maps[i][k]]
                assert left == right
                counts['socle_multiplication_commutators'] += 1
    for a in D:
        image = index[(0,)*len(q)]
        for i, exponent in enumerate(a):
            for _ in range(exponent):
                image = None if image is None else maps[i][image]
        assert image == column(a)
        counts['socle_cyclic_column_replays'] += 1
    assert column(u) == column(v)
    for a in (u, v):
        assert all(tuple(x+(i == j) for j, x in enumerate(a)) not in D
                   for i in range(len(q)))
    for omitted in (u, v):
        state = D-{omitted}
        assert len(state) == n
        assert all(tuple(x-(i == j) for j, x in enumerate(a)) in state
                   for a in state for i, exponent in enumerate(a) if exponent)
    return n, u, v


def run_case(q, counts):
    assert len(q) == 3
    n, u, v = socle_model(q, counts)
    A, B = prod(x+1 for x in u), prod(x+1 for x in v)
    p, s = sum(x > 0 for x in q), sum(x < 0 for x in q)
    denominators = [p*x if x > 0 else -s*x if x < 0 else 1 for x in q]
    scale = lcm(*denominators)
    original = tuple(scale//den if x else 4*scale for x, den in zip(q, denominators))
    assert dot(original, u) == dot(original, v) == scale
    i, j = [i for i, x in enumerate(q) if x][:2]
    h = tuple(q[j] if k == i else -q[i] if k == j else 0 for k in range(3))
    assert dot(h, q) == 0 and any(h)
    coarse = tuple(3*den if x else 0 for x, den in zip(q, denominators))
    denominator = 1+10*sum(abs(c)*(bound+a+1) for c, bound, a in zip(h, coarse, u))
    w = tuple(denominator*a+b for a, b in zip(original, h))
    assert min(w) > 0 and dot(w, q) == 0
    assert all(abs(Fraction(b, denominator)) < Fraction(a, 4)
               for a, b in zip(original, h))
    threshold = dot(w, u)
    ground = tuple(product(*(range(threshold//a+1) for a in w)))
    cheaper = frozenset(a for a in ground if dot(w, a) < threshold)
    tied = frozenset(a for a in ground if dot(w, a) == threshold)
    assert tied == {u, v}
    size = len(cheaper)+1
    half_box = tuple(product(*(range(abs(x)//2+1) for x in q)))
    assert all(dot(w, a) < threshold for a in half_box)
    assert len(half_box)*8 >= A*B and size > len(half_box)
    assert size <= 6**3*A*B
    perturb = 1+10*sum(abs(x)*(bound+a+1) for x, bound, a in zip(q, coarse, u))
    sides = []
    for sign, keep in ((1, v), (-1, u)):
        price = tuple(perturb*a+sign*x for a, x in zip(w, q))
        assert min(price) > 0
        selected = first_terms(price, size)
        assert selected == cheaper | {keep}
        assert all(tuple(x-(i == j) for j, x in enumerate(a)) in selected
                   for a in selected for i, exponent in enumerate(a) if exponent)
        assert all(a[i] <= 2*3*abs(q[i]) if q[i] else a[i] == 0
                   for a in selected for i in range(3))
        sides.append({'sign': sign, 'strict_weight': price, 'selected_endpoint': keep})
        counts['independent_unrestricted_heap_cut_replays'] += 1
    counts['exact_two_endpoint_wall_witnesses'] += 1
    counts['direct_free_threshold_box_columns'] += len(ground)
    counts['universal_half_box_columns'] += len(half_box)
    return {'primitive_direction': q, 'minimum_algebraic_colength': n,
            'positive_divisor_count': A, 'negative_divisor_count': B,
            'positive_endpoint': u, 'negative_endpoint': v,
            'exact_wall_weight': w, 'common_cost': threshold,
            'free_cut_size_witness': size, 'half_box_lower_count': len(half_box),
            'adjacent_free_cuts': sides}


def primitive(q):
    divisor = gcd(*q)
    assert divisor
    return canonical(tuple(x//divisor for x in q))


def wall_cost(q):
    return prod(max(x, 0)+1 for x in q)+prod(max(-x, 0)+1 for x in q)-2


def saved_wall_budgets(counts):
    artifact = json.loads((ROOT/'results/research/oracle_enumeration_checks.json').read_text())
    fixtures = []
    for case in artifact['direct_positive_cone_traversal']['fixtures']:
        n, edges = case['length'], {}
        for edge in case['directed_neighbor_trace']:
            start, end = tuple(edge['from_profile']), tuple(edge['to_profile'])
            q = primitive(tuple(x-y for x, y in zip(start, end)))
            a, b = [tuple(Fraction(x) for x in row) for row in edge['shared_edge']]
            normal = (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
                      a[0]*b[1]-a[1]*b[0])
            denominator = lcm(*(x.denominator for x in normal))
            assert primitive(tuple(int(denominator*x) for x in normal)) == q
            cost = wall_cost(q)
            assert 2 <= cost <= n and sum(map(abs, q)) <= cost
            pair = tuple(sorted((start, end)))
            if pair in edges:
                assert edges[pair] == (q, cost)
            edges[pair] = q, cost
        value = sum(cost for q, cost in edges.values())
        assert case['initial_ideals'] <= len(edges)+1
        for threshold in range(1, n+2):
            small = sum(cost <= threshold for q, cost in edges.values())
            large = len(edges)-small
            assert large*threshold <= value
            assert case['initial_ideals'] <= 1+small+Fraction(value, threshold)
            counts['independent_wall_threshold_partitions'] += 1
        counts['independent_undirected_saved_wall_costs'] += len(edges)
        fixtures.append({'model': case['model'], 'length': n,
                         'characteristic': case['characteristic'],
                         'initial_ideals': case['initial_ideals'],
                         'undirected_walls': len(edges), 'primitive_wall_budget': value,
                         'wall_costs': [{'from_profile': pair[0], 'to_profile': pair[1],
                                         'primitive_direction': q, 'cost': cost}
                                        for pair, (q, cost) in sorted(edges.items())]})
    return fixtures


def three_pole_budgets(counts):
    fixtures = []
    for n in range(2, 33):
        allocations = [(a, b, n-1-a-b) for a in range(n) for b in range(n-a)]
        moments = {a: tuple(x*(x+1)//2 for x in a) for a in allocations}
        edges = set()
        for a in allocations:
            for i, length in enumerate(a):
                if not length:
                    continue
                for j in range(3):
                    if i == j:
                        continue
                    b = tuple(x-(k == i)+(k == j) for k, x in enumerate(a))
                    edges.add(tuple(sorted((a, b))))
        budget = 0
        for a, b in edges:
            q = primitive(tuple(x-y for x, y in zip(moments[a], moments[b])))
            assert sum(x != 0 for x in q) == 2
            cost = wall_cost(q)
            assert cost <= n
            budget += cost
            counts['independent_three_pole_allocation_edges'] += 1
        formula = 3*sum((A+B)//gcd(A, B) for A in range(1, n)
                        for B in range(1, n-A+1))
        assert budget == formula and len(edges) == 3*n*(n-1)//2
        counts['exact_three_pole_budget_equalities'] += 1
        if n in (2, 3, 4, 8, 16, 32):
            fixtures.append({'length': n, 'initial_ideals': len(allocations),
                             'undirected_walls': len(edges), 'primitive_wall_budget': budget})
    return fixtures


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    counts = defaultdict(int)
    directions = {(1, -1, 0), (2, 0, -3), (0, 3, -4), (1, -2, -3),
                  (3, -2, -2), (2, -3, 1), (1, -1, -1)}
    directions.update((k, -k-1, 0) for k in range(1, 13))
    directions.update((k*k, -k, -k-1) for k in range(1, 9))
    fixtures = [run_case(q, counts) for q in sorted({canonical(q) for q in directions})]
    budgets = saved_wall_budgets(counts)
    poles = three_pole_budgets(counts)
    sources = ('docs/free_wall_threshold.md', 'scripts/run_free_wall_checks.py',
               'docs/primitive_wall_spectrum.md', 'docs/rational_shadow_bound.md',
               'docs/ray_scaling_bound.md', 'results/research/oracle_enumeration_checks.json')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_free_wall_checks.py',
              'counts': dict(counts), 'fixtures': fixtures,
              'saved_fan_wall_budgets': budgets, 'actual_three_pole_budgets': poles,
              'uniform_results': {
                  'free_first_wall': 'For fixed d>=2, 2^(-d)*A(q)*B(q) < mu_d(q) <= (2*d)^d*A(q)*B(q).',
                  'algebraic_comparison': 'The minimum algebraic colength is A(q)+B(q)-2; its uniform proof and priority boundary are in docs/primitive_wall_spectrum.md.',
                  'sharp_delay': 'max_(nu(q)<=n) mu_d(q) = Theta_d(n^2); at-most length in d=2, exact length also in d>=3.',
                  'full_support_family': 'q_k=(k^2,-k,-k-1), nu=2*k^2+3*k+1, mu_3=Theta(k^4).',
                  'collective_parameter_bound': 'g(I)=O(1+n*L^2+n^(3/10)*V(I)^(7/10)*L^(3/5)), where V sums primitive minimum-colength costs over undirected positive adjacencies.',
                  'unproved_sufficient_condition': 'A uniform V=O(n^3*L^C) would imply G3=O(n^(12/5)*L^((7*C+6)/10)). The cubic hypothesis is not proved.',
                  'actual_cubic_control': 'The existing three-pole family has V=3*sum_(A,B>=1,A+B<=n) (A+B)/gcd(A,B)=Theta(n^3).'},
              'new_unrestricted_G3_count': False, 'polynomial_lower_separation_proved': False,
              'novelty_established': False, 'goal_complete': False,
              'scope': 'Independent direct lattice threshold enumeration and unrestricted heap replay verify finite actual free walls. Socle unit-column maps replay commuting cyclic algebras, valid over every field. The free witness size need not be the minimum mu; the uniform inequalities rest on the proof, not fitted ratios. Saved fan wall costs and all three-pole allocation edges through n=32 are checked independently. The parameterized bound retains V; its cubic hypothesis, a new unrestricted exponent, current priority and a universal-prefix algorithm are not established.',
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    (ROOT/'results/research/free_wall_threshold_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'counts': report['counts'], 'fixtures': len(fixtures)}, indent=2))


if __name__ == '__main__':
    main()
