#!/usr/bin/env python3
"""Independent exact controls for the uniform cubic classification proof."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from heapq import merge
from itertools import permutations, product
from math import gcd
from pathlib import Path
import hashlib
import json

from run_hypersurface_truncation_checks import (
    integral_border_relations, integral_cubic_model, simplex, sparse_basis,
)
from run_structural_complexity_checks import dot, greedy

ROOT = Path(__file__).resolve().parents[1]


def cap(value, maximum):
    return max(0, min(value, maximum))


def floor(value):
    return value.numerator//value.denominator if isinstance(value, F) else value//1


def ceil(value):
    return -floor(-value)


def ordered_stream_basis(s, lam):
    """The boundary basis in ordered coordinates; merge work is O((s+1)^2)."""
    chosen = [(0, 0, c) for c in range(s+1)]
    for k in range(1, s+1):
        L = s-k+1
        streams = (
            [(lam*a, (a, k-a, 0)) for a in range(k+1)],
            [(F(c), (0, k, c)) for c in range(1, L)],
            [(lam*k+c, (k, 0, c)) for c in range(1, L)],
        )
        for j, (_, exponent) in enumerate(merge(*streams)):
            if j == 2*L:
                break
            chosen.append(exponent)
    assert len(chosen) == (s+1)**2
    return frozenset(chosen)


def cubic_standard_set(s, weight):
    order = sorted(range(3), key=lambda i: (-weight[i], i))
    wx, wy, wz = (F(weight[i]) for i in order)
    assert wz > 0
    if wx >= wy+wz:
        ordered = frozenset((a, b, c) for a in range(min(1, s)+1)
                            for b in range(s-a+1) for c in range(s-a-b+1))
    else:
        ordered = ordered_stream_basis(s, (wx-wy)/wz)
    return frozenset(tuple(a[order.index(i)] for i in range(3)) for a in ordered)


def cubic_prefix_rank(s, weight, price):
    """Original monomial price-prefix rank, O(s) exact operations, all fields."""
    wx, wy, wz = sorted(map(F, weight), reverse=True)
    price = F(price)
    assert wz > 0
    if wx >= wy+wz:
        return sum(cap(floor((price-wx*a-wy*b)/wz)+1, s-a-b+1)
                   for a in range(min(1, s)+1) for b in range(s-a+1))
    u, lam, T = wy/wz, (wx-wy)/wz, price/wz
    result = cap(floor(T)+1, s+1)
    for k in range(1, s+1):
        L, tau = s-k+1, T-u*k
        spine = (k+1 if tau >= 0 else 0) if not lam else cap(floor(tau/lam)+1, k+1)
        result += min(2*L, spine+cap(floor(tau), L-1)+cap(floor(tau-lam*k), L-1))
    return result


def isosceles_channel_rank(s, lam, price):
    """Independent four-channel formula, before the rank-multiset identity."""
    lam, price = F(lam), F(price)
    assert 0 < lam < 1
    result = 0
    for delta in (0, 1, 1, 2):
        S, T = s-delta, price-delta
        for q in range(S+1):
            L = min(q, S-q)+1
            low = max(0, ceil(((1+lam)*q-T)/(2*lam)))
            high = min(q, floor(T/2))
            result += cap(1+high-low, L)
    return result


def cubic_walls(s):
    return tuple(sorted({F(p, q) for q in range(2, s+2) for p in range(1, q)
                         if gcd(p, q) == 1 and
                         (1 if p % 2 and q % 2 else 2)*(3*q+p) <= 2*(s+1)}))


def block_controls():
    counts = defaultdict(int)
    for k, L in product(range(1, 10), repeat=2):
        ground = tuple(product(range(k+1), range(L)))
        columns = {}
        for a, c in ground:
            vector = {(0, c): 1}
            for _ in range(a):
                following = defaultdict(int)
                for (i, j), value in vector.items():
                    if not i:
                        following[1, j] += value
                    else:
                        following[0, j] -= value
                        if j+1 < L:
                            following[1, j+1] -= value
                vector = {b: value for b, value in following.items() if value}
            columns[a, c] = tuple(vector.get((i, j), 0) for i in range(2) for j in range(L))
        for lam, p in product((F(1, 37), F(7, 37), F(19, 37), F(31, 37)), (0, 2, 3, 5)):
            linear = greedy(ground, columns, (lam, F(1)), p)
            arms = [(a, 0) for a in range(k+1)]+[(0, c) for c in range(1, L)]+[(k, c) for c in range(1, L)]
            predicted = set(sorted(arms, key=lambda a: (dot(a, (lam, F(1))), a))[:2*L])
            assert linear == predicted
            for T in {lam*a+c for a, c in ground}:
                actual = sum(lam*a+c <= T for a, c in linear)
                expected = min(2*L, cap(floor(T/lam)+1, k+1)
                               +cap(floor(T), L-1)+cap(floor(T-lam*k), L-1))
                assert actual == expected
                counts['binary_block_price_prefixes'] += 1
            counts['binary_block_field_bases'] += 1
    return dict(counts)


def original_rank_controls():
    counts, fixtures = defaultdict(int), []
    lams = (F(1, 37), F(7, 37), F(19, 37), F(31, 37))
    for s in (0, 1, 2, 3, 4, 6, 8, 10):
        ground = simplex(s)
        basis, evaluate = integral_cubic_model(s)
        columns = {}
        for a in ground:
            column = evaluate(a)
            columns[a] = tuple(column.get(k, 0) for k in range(len(basis)))
        weights = [(1+lam, F(1), F(1)) for lam in lams]
        weights += [(F(17, 13)+lam, F(17, 13), F(1)) for lam in lams]
        weights += [(F(5), F(1), F(2)), (F(2), F(1), F(1)), (F(1), F(1), F(1))]
        weights += [(w[2], w[0], w[1]) for w in weights[:4]]
        for p, weight in product((0, 2, 3, 5, 101), weights):
            chosen = greedy(ground, columns, weight, p)
            assert len(chosen) == (s+1)**2
            stream_set = cubic_standard_set(s, weight)
            # Comparing tier multiplicities also covers weights with unavoidable ties.
            assert Counter(dot(a, weight) for a in chosen) == Counter(dot(a, weight) for a in stream_set)
            if weight[0] > weight[1] > weight[2] and weight[0] < weight[1]+weight[2]:
                assert chosen == stream_set
                counts['strict_original_standard_set_equalities'] += 1
            for T in sorted({dot(a, weight) for a in ground}):
                actual = sum(dot(a, weight) <= T for a in chosen)
                assert actual == cubic_prefix_rank(s, weight, T)
                if weight[1] == weight[2] == 1 and 1 < weight[0] < 2:
                    assert actual == isosceles_channel_rank(s, weight[0]-1, T)
                    counts['independent_four_channel_prefix_ranks'] += 1
                counts['original_field_price_prefix_ranks'] += 1
            counts['original_field_basis_models'] += 1
        fixtures.append({'cutoff': s, 'length': len(basis), 'characteristics': [0, 2, 3, 5, 101],
                         'weights_per_field': len(weights)})
    return {'counts': dict(counts), 'fixtures': fixtures,
            'scope': 'Exact Gaussian ranks in original monomial columns, across entire price tiers, against both independent formulas and the stream standard sets. The original integral kernel is separately verified against full Macaulay kernels in the retained truncation report.'}


def cutoff_multiset_controls():
    counts = defaultdict(int)
    for s, lam in product((1, 2, 3, 5, 7, 11, 17, 31, 60), (F(1, 127), F(37, 127), F(79, 127), F(113, 127))):
        original = []
        for delta in (0, 1, 1, 2):
            for q in range(s-delta+1):
                L = min(q, s-delta-q)+1
                kappa = sorted(max(F(2*i), (1+lam)*q-2*lam*i) for i in range(q+1))
                original.extend(delta+value for value in kappa[:L])
        boundary = [dot(a, (1+lam, F(1), F(1))) for a in ordered_stream_basis(s, lam)]
        assert Counter(original) == Counter(boundary)
        counts['cutoff_price_multiset_identities'] += 1
        counts['matched_price_items'] += len(original)
    # Directly check the superlevel identity used in the uniform proof.
    for D, lam in product(range(2, 70), (F(1, 127), F(37, 127), F(79, 127), F(113, 127))):
        boundary = Counter(floor(min((1-lam)*k/2, (1-lam)*(D-k)/(2*lam))) for k in range(1, D))
        original = Counter()
        for delta in (0, 1, 1, 2):
            if D-delta < 2 or (D-delta) % 2:
                continue
            m = (D-delta)//2
            for b in range(1, m+1):
                original[min(m-b, floor((1-lam)*b/(2*lam)))] += 1
        assert original == boundary
        counts['mixed_plane_superlevel_multisets'] += 1
    return dict(counts)


def fan_controls():
    saved = json.loads((ROOT/'results/research/hypersurface_truncation_checks.json').read_text())
    counts, fixtures = defaultdict(int), []
    by_cutoff = {f['cutoff']: f for f in saved['all_field_cubic_boundary_certificates']['fixtures']}
    for s in range(3, 21):
        walls = cubic_walls(s)
        endpoints = (F(0),)+walls+(F(1),)
        states = set()
        for lo, hi in zip(endpoints, endpoints[1:]):
            B = ordered_stream_basis(s, (lo+hi)/2)
            for order in permutations(range(3)):
                states.add(frozenset(tuple(a[i] for i in order) for a in B))
        for i in range(3):
            states.add(frozenset(a for a in simplex(s) if a[i] < 2))
        assert len(states) == 6*(len(walls)+1)
        if s in by_cutoff:
            old = by_cutoff[s]
            old_walls = tuple(F(v['lambda_interval'][1]) for v in old['ordered_central_sectors'][:-1])
            assert walls == old_walls
            assert len(states) == old['complete_positive_fan_cells']
            counts['retained_complete_fan_comparisons'] += 1
        fixtures.append({'cutoff': s, 'length': (s+1)**2, 'positive_fan_cells': len(states),
                         'positive_ordered_walls': [str(v) for v in walls]})
        counts['exact_stream_fan_counts'] += 1
    # Fresh monic certificates at odd cutoffs absent from the retained table.
    for s in (7, 9, 11, 13):
        _, evaluate = integral_cubic_model(s)
        endpoints = (F(0),)+cubic_walls(s)+(F(1),)
        for lo, hi in zip(endpoints, endpoints[1:]):
            B = ordered_stream_basis(s, (lo+hi)/2)
            chosen, pivots = sparse_basis(sorted(B), evaluate)
            assert len(chosen) == len(B)
            normals, stats = integral_border_relations(chosen, pivots, evaluate, s, 1)
            assert all(x+y >= 0 and min(x+y+x*lo+z, x+y+x*hi+z) >= 0 for x, y, z in normals)
            counts['fresh_odd_entire_sector_certificates'] += 1
            counts['fresh_odd_monic_integer_relations'] += stats[0]
    return {'counts': dict(counts), 'fixtures': fixtures,
            'scope': 'The exact wall formula and coordinate-permutation count are checked against every retained all-field complete fan, with fresh original monic certificates at four intervening odd cutoffs. No fitted asymptotic is used.'}


def main():
    if not __debug__:
        raise RuntimeError('Run without -O: exact assertions are required.')
    block = block_controls()
    original = original_rank_controls()
    multisets = cutoff_multiset_controls()
    fans = fan_controls()
    sources = ('scripts/run_cubic_fan_checks.py', 'docs/cubic_fan_classification.md',
               'scripts/run_hypersurface_truncation_checks.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_socle_stability_checks.py',
               'results/research/hypersurface_truncation_checks.json')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_cubic_fan_checks.py',
              'counts': {**block, **original['counts'], **multisets, **fans['counts']},
              'binary_block_controls': block, 'original_monomial_rank_controls': original,
              'cutoff_multiset_controls': multisets, 'fan_controls': fans,
              'uniform_family_theorem': 'For s>=3 over every field, the truncated full-support diagonal Markoff cubic has exactly 6*(1+M_s) monomial initial ideals, with epsilon(p,q)*(3*q+p)<=2*(s+1). Its count is 3*n/pi^2+O(sqrt(n)*log(n+1)).',
              'algorithmic_consequence': 'Original monomial price-prefix ranks use O(sqrt(n)) exact arithmetic/floor operations and constant auxiliary storage; a standard set uses O(n) comparisons. These are operation counts for the explicit family, not general bit-complexity or state-of-the-art performance claims.',
              'scope': 'Finite independent exact controls support the complete uniform proof in the note. The theorem closes this deep-relation family as a polynomial-amplification route; unrestricted G_3, novelty and publication gates remain open.',
              'unrestricted_G3_improvement': False, 'novelty_established': False, 'goal_complete': False,
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    (ROOT/'results/research/cubic_fan_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
