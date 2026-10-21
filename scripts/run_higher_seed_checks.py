#!/usr/bin/env python3
"""Exact higher-seed binomial reductions; the uniform core identity is open."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd
from pathlib import Path
import hashlib
import json

from run_higher_interaction_checks import ground, linear_sum, rank_trace, smith_units
from run_hypersurface_truncation_checks import relation_rows
from run_structural_complexity_checks import determinant, dot

ROOT = Path(__file__).resolve().parents[1]


def core_ground(s, a):
    return tuple((i, j, c) for i in range(s//a+1)
                 for j in range(s//a-i+1) for c in range(s-a*(i+j)+1))


def core_model(s, a, r, canonical=False):
    rows = tuple((b, d) for b in range(s//a+1) for d in range(s-a*b+1))
    index = {v: i for i, v in enumerate(rows)}

    @lru_cache(None)
    def evaluate(v):
        i, j, c = v
        out = {}
        if canonical:
            if i == 0:
                return {index[j, c]: 1}
            for h in range((s-a*(i+j)-c)//r+1):
                out[index[j, c+a*i+r*h]] = (-1)**(i+h)*comb(i+h-1, h)
        else:
            for h in range(i+1):
                b, d = j+h, c+a*i+(r-a)*h
                if a*b+d <= s:
                    out[index[b, d]] = (-1)**i*comb(i, h)
        return out
    return rows, evaluate


def column_rank(columns, p):
    pivots = {}
    for col in columns:
        v = {k: (c % p if p else F(c)) for k, c in col.items()
             if (c % p if p else c)}
        while v:
            k = min(v)
            if k not in pivots:
                inv = pow(v[k], -1, p) if p else 1/v[k]
                pivots[k] = {j: c*inv % p if p else c*inv for j, c in v.items()}
                break
            c = v[k]
            for j, d in pivots[k].items():
                value = v.get(j, 0)-c*d
                if p:
                    value %= p
                if value:
                    v[j] = value
                else:
                    v.pop(j, None)
    return len(pivots)


def regular_blocks(s, a, r):
    for k in range(-(r-a)*(s//a), s+1):
        b0 = max(0, -(k//(r-a)))
        b1 = (s-k)//r
        if b0 <= b1:
            yield k, b0, b1-b0+1


def ceil(x):
    return -((-x.numerator)//x.denominator)


def core_rank(s, a, r, u, Z, T, p, canonical=False):
    if s < 0:
        return 0
    u, Z, T = map(F, (u, Z, T))
    if not (r > a >= 1 and 0 < u < Z):
        raise ValueError('This block query requires r>a>=1 and 0<u<Z.')
    A = a*(Z-u)
    out = 0
    if not canonical:
        B = a*u+(r-a)*Z
        for k, b0, L in regular_blocks(s, a, r):
            def columns():
                for j in range(b0, b0+L):
                    lo = max(0, ceil((k*Z+B*j-T)/A))
                    hi = (k+(r-a)*j)//a
                    for i in range(lo, hi+1):
                        yield {j-b0+h: (-1)**i*comb(i, h)
                               for h in range(min(i+1, L-j+b0))}
            out += column_rank(columns(), p)
    else:
        for b in range(s//a+1):
            for k in range(min(r, s-a*b+1)):
                L = (s-a*b-k)//r+1
                def columns():
                    for j in range(L):
                        lo = max(0, ceil((a*u*b+k*Z+r*Z*j-T)/A))
                        hi = (k+r*j)//a
                        for i in range(lo, hi+1):
                            if i == 0:
                                yield {j: 1}
                            else:
                                yield {j+h: (-1)**(i+h)*comb(i+h-1, h)
                                       for h in range(L-j)}
                out += column_rank(columns(), p)
    return out


def form_model(s, a, r, Q, H):
    """Independent finite original z<a recurrence, with increasing-degree tails."""
    basis = tuple(v for v in ground(s) if v[2] < a)
    index = {v: i for i, v in enumerate(basis)}

    @lru_cache(None)
    def evaluate(v):
        x, y, z = v
        if sum(v) > s:
            return {}
        if z < a:
            return {index[v]: 1}
        tails = [(-c, evaluate((x+i, y+j, z-a))) for (i, j), c in Q.items()]
        tails += [(-c, evaluate((x+i, y+j, z-a+r))) for (i, j), c in H.items()]
        return linear_sum(tails)
    return basis, evaluate


def summed_rank(s, a, r, u, Z, T, p, canonical=False):
    return sum(m*core_rank(s-delta, a, r, u, Z, T-delta*u, p, canonical)
               for delta, m in Counter(i+j for i in range(a) for j in range(a)).items())


def block_controls():
    counts = defaultdict(int)
    for s, a, r in ((6, 1, 2), (8, 2, 3), (11, 2, 4), (12, 3, 5), (17, 2, 3)):
        E = core_ground(s, a)
        for canonical in (False, True):
            rows, ev = core_model(s, a, r, canonical)
            if not canonical:
                maps = {k: (b0, L) for k, b0, L in regular_blocks(s, a, r)}
                for i, j, c in E:
                    k = c+a*i+(a-r)*j
                    b0, L = maps[k]
                    actual = {(rows[n][0]-b0): x for n, x in ev((i, j, c)).items()}
                    expected = {j-b0+h: (-1)**i*comb(i, h)
                                for h in range(min(i+1, L-j+b0))}
                    assert actual == expected
                    counts['regular_integral_block_column_equalities'] += 1
            else:
                for i, b, c in E:
                    k, j = (a*i+c) % r, (a*i+c)//r
                    L = (s-a*b-k)//r+1
                    actual = {((rows[n][1]-k)//r): x for n, x in ev((i, b, c)).items()}
                    expected = ({j: 1} if not i else
                                {j+h: (-1)**(i+h)*comb(i+h-1, h) for h in range(L-j)})
                    assert actual == expected
                    counts['canonical_integral_block_column_equalities'] += 1
            for p, u in product((0, 2, 3, 5), (F(2, 7), F(5, 7), F(17, 19))):
                cost = lambda v: a*u*(v[0]+v[1])+v[2]
                B, _ = rank_trace(sorted(E, key=lambda v: (cost(v), v)), ev, p, set())
                for T in (F(0), F(2), F(s, 3), F(2*s, 3), F(s)):
                    assert core_rank(s, a, r, u, 1, T, p, canonical) == sum(cost(v) <= T for v in B)
                    counts['independent_core_block_rank_queries'] += 1
        if s <= 8:
            for canonical, p in product((False, True), (0, 2, 3)):
                rows, ev = core_model(s, a, r, canonical)
                f = {(1, 0, 0): 1, (0, 0, a): 1,
                     (1, 0, r) if canonical else (0, 1, r): 1}
                relations = relation_rows(E, [f], p)
                assert len(E)-len(relations) == len(rows)
                for rel in relations:
                    image = linear_sum((x, ev(v)) for v, x in zip(E, rel) if x)
                    assert all(not (x % p if p else x) for x in image.values())
                    counts['independent_core_macaulay_rows'] += 1
                counts['independent_core_full_kernel_equalities'] += 1
    return dict(counts)


def form_controls():
    counts, fixtures = defaultdict(int), []
    pairs = ((2, 3, {(2, 0): 1, (0, 2): 1}, {(1, 1): 1}, True),
             (2, 3, {(2, 0): 1}, {(1, 1): 1}, False),
             (2, 4, {(0, 2): 1}, {(0, 2): 1}, False),
             (2, 5, {(2, 0): 1}, {(0, 2): 1}, True),
             (3, 5, {(3, 0): 1, (0, 3): 1}, {(2, 1): 1}, True),
             (3, 4, {(3, 0): 1}, {(2, 1): 1}, False),
             (1, 3, {(1, 0): 1}, {(0, 1): 1}, True),
             (1, 2, {(1, 0): 1}, {(1, 0): 1}, False))
    for label, (a, r, Q, H, regular) in enumerate(pairs):
        s = a+r+4
        E, (basis, ev) = ground(s), form_model(s, a, r, Q, H)
        assert len(basis) == sum(comb(s-e+2, 2) for e in range(min(a, s+1)))
        for p, u in product((0, 2, 3, 5), (F(2, 7), F(5, 7), F(17, 19))):
            cost = lambda v: dot(v, (u, u, F(1)))
            B, _ = rank_trace(sorted(E, key=lambda v: (cost(v), v)), ev, p, set())
            for T in (F(0), F(2), F(s, 2), F(s)):
                actual = sum(cost(v) <= T for v in B)
                low = summed_rank(s, a, r, u, 1, T, p, True)
                high = summed_rank(s, a, r, u, 1, T, p)
                assert low <= actual <= high
                if regular:
                    assert actual == high
                counts['independent_original_rank_sandwiches'] += 1
                counts['finite_matching_bounds'] += int(low == high)
            counts['independent_original_field_models'] += 1
        fixtures.append({'pair': label, 'degree': a, 'interaction': r, 'cutoff': s,
                         'regular': regular, 'Q': [[*v, c] for v, c in Q.items()],
                         'H': [[*v, c] for v, c in H.items()]})
    return dict(counts), fixtures


@lru_cache(None)
def binary_module(x, y):
    """Integral basis 1,x,y,x^2 over U=x^2+y^2,V=xy."""
    def shift(v, i, j):
        return {(token, a+i, b+j): c for (token, a, b), c in v.items()}
    if x and y:
        return shift(binary_module(x-1, y-1), 0, 1)
    if x >= 3:
        return linear_sum(((1, shift(binary_module(x-2, 0), 1, 0)),
                           (-1, shift(binary_module(x-3, 1), 0, 1))))
    if y == 2:
        return {(0, 1, 0): 1, (3, 0, 0): -1}
    if y >= 3:
        return linear_sum(((1, shift(binary_module(0, y-2), 1, 0)),
                           (-1, shift(binary_module(1, y-3), 0, 1))))
    token = 1 if x == 1 else 3 if x == 2 else 2 if y == 1 else 0
    return {(token, 0, 0): 1}


def integral_original_controls():
    counts = defaultdict(int)
    for s in (4, 7, 10):
        degrees = (0, 1, 1, 2)
        models = [core_model(s-d, 2, 3) for d in degrees]
        rows = tuple((token, *v) for token, (basis, _) in enumerate(models) for v in basis)
        index = {v: i for i, v in enumerate(rows)}

        @lru_cache(None)
        def mapped(v):
            x, y, z = v
            out = defaultdict(int)
            for (token, i, j), c in binary_module(x, y).items():
                basis, ev = models[token]
                for row, d in ev((i, j, z)).items():
                    out[index[token, *basis[row]]] += c*d
            return {k: c for k, c in out.items() if c}

        Q, H = {(2, 0): 1, (0, 2): 1}, {(1, 1): 1}
        basis, original = form_model(s, 2, 3, Q, H)
        assert len(basis) == len(rows) == (s+1)**2
        for v in ground(s):
            assert mapped(v) == linear_sum((c, mapped(basis[i])) for i, c in original(v).items())
            counts['original_integral_core_column_equalities'] += 1
        M = [[mapped(v).get(i, 0) for v in basis] for i in range(len(rows))]
        units, residual = smith_units(M)
        assert units == len(rows) and all(not x for row in residual for x in row)
        counts['original_integral_reference_unit_pivots'] += units
    return dict(counts)


def smith_diagonal(matrix):
    """Euclidean unimodular row/column operations on a small residual."""
    A = [list(row) for row in matrix]
    if not A or not A[0]:
        return []
    m, n, out = len(A), len(A[0]), []
    for k in range(min(m, n)):
        pos = next(((i, j) for i in range(k, m) for j in range(k, n) if A[i][j]), None)
        if pos is None:
            break
        i, j = pos
        A[k], A[i] = A[i], A[k]
        for row in A:
            row[k], row[j] = row[j], row[k]
        while True:
            done = False
            for i in range(k+1, m):
                if A[i][k]:
                    q = A[i][k]//A[k][k]
                    A[i] = [x-q*y for x, y in zip(A[i], A[k])]
                    if A[i][k]:
                        A[k], A[i] = A[i], A[k]
                    done = True
                    break
            if done:
                continue
            for j in range(k+1, n):
                if A[k][j]:
                    q = A[k][j]//A[k][k]
                    for row in A:
                        row[j] -= q*row[k]
                    if A[k][j]:
                        for row in A:
                            row[k], row[j] = row[j], row[k]
                    done = True
                    break
            if done:
                continue
            pos = next(((i, j) for i in range(k+1, m) for j in range(k+1, n)
                        if A[i][j] % A[k][k]), None)
            if pos is None:
                break
            i, j = pos
            A[k] = [x+y for x, y in zip(A[k], A[i])]
        out.append(abs(A[k][k]))
    assert all(not b % a for a, b in zip(out, out[1:]))
    return out


def smith_signature(rows, ev, selected):
    # Selected reference unit columns are eliminated by deleting their rows.
    units = {next(iter(ev(v))) for v in selected if v[0] == 0}
    remaining = [i for i in range(len(rows)) if i not in units]
    columns = [ev(v) for v in selected if v[0]]
    if not remaining or not columns:
        return len(units), []
    M = [[col.get(i, 0) for col in columns] for i in remaining]
    q, residual = smith_units(M)
    diagonal = smith_diagonal(residual)
    return len(units)+q+sum(x == 1 for x in diagonal), [x for x in diagonal if x > 1]


def smith_controls():
    counts, fixtures = defaultdict(int), []
    small = ([[6, 0], [0, 10]], [[2, 4, 6], [4, 2, 8]],
             [[0, 6, 9], [3, 0, 12], [15, 6, 3]], [[0, 0], [0, 0]],
             [[2, 0, 4], [0, 6, 9], [4, 12, 26]])
    for M in small:
        diag = smith_diagonal(M)
        divisor = 1
        for k in range(1, min(len(M), len(M[0]))+1):
            actual = 0
            for I in combinations(range(len(M)), k):
                for J in combinations(range(len(M[0])), k):
                    value = determinant([[M[i][j] for j in J] for i in I])
                    assert value.denominator == 1
                    actual = gcd(actual, abs(int(value)))
            if k <= len(diag):
                divisor *= diag[k-1]
                assert actual == divisor
            else:
                assert actual == 0
            counts['independent_determinantal_divisor_checks'] += 1
    for s, a, r, H in ((11, 2, 3, F(7)), (17, 2, 3, F(13)),
                        (23, 2, 3, F(19)), (19, 2, 4, F(31, 7)), (21, 3, 5, F(17, 3))):
        E = core_ground(s, a)
        models = [core_model(s, a, r, canonical) for canonical in (False, True)]
        price = lambda v: a*(v[0]+v[1])+H*v[2]
        prices = sorted(set(map(price, E)))
        tests = sorted({2*H, F(s), s*H}|set(prices[::max(1, len(prices)//20)]))
        prefixes = []
        for T in tests:
            selected = [v for v in E if price(v) <= T]
            signatures = [smith_signature(rows, ev, selected) for rows, ev in models]
            assert signatures[0] == signatures[1]
            counts['finite_integral_smith_signature_matches'] += 1
            prefixes.append({'price_cutoff': str(T), 'unit_invariants': signatures[0][0],
                             'nonunit_invariants': signatures[0][1]})
        fixtures.append({'cutoff': s, 'degree': a, 'interaction': r,
                         'weight': [str(a), str(a), str(H)], 'prefixes': prefixes})
    return dict(counts), fixtures


def trace_inventory(E, ev, p, cost, n):
    pivots, out = {}, Counter()
    for v in E:
        if p == 2:
            x = sum(1 << k for k, c in ev(v).items() if c % 2)
            while x:
                k = (x & -x).bit_length()-1
                if k not in pivots:
                    pivots[k] = x
                    out[cost(v)] += 1
                    break
                x ^= pivots[k]
        else:
            x = {k: c % p for k, c in ev(v).items() if c % p}
            while x:
                k = min(x)
                if k not in pivots:
                    inv = pow(x[k], -1, p)
                    pivots[k] = {j: c*inv % p for j, c in x.items()}
                    out[cost(v)] += 1
                    break
                c = x[k]
                for j, d in pivots[k].items():
                    value = (x.get(j, 0)-c*d) % p
                    if value:
                        x[j] = value
                    else:
                        x.pop(j, None)
        if len(pivots) == n:
            break
    return out


def entire_seed_controls():
    counts, fixtures = defaultdict(int), []
    for s, a, r, p in ((41, 2, 3, 2), (35, 2, 3, 3), (35, 2, 3, 5)):
        E = core_ground(s, a)
        models = [core_model(s, a, r, canonical) for canonical in (False, True)]
        n = len(models[0][0])
        walls = sorted({F(1), F(s+1)}|{F(a*k, b) for k in range(1, s//a+1)
                                       for b in range(1, min(s, a*k-1)+1)})
        intervals = []
        for step, (lo, hi) in enumerate(zip(walls, walls[1:])):
            H = (lo+hi)/2
            price = lambda v: a*(v[0]+v[1])+H*v[2]
            ordered = sorted(E, key=lambda v: (price(v), v))
            inventories = [trace_inventory(ordered, ev, p, price, n) for _, ev in models]
            if step == 0:
                for (_, ev), inventory in zip(models, inventories):
                    B, _ = rank_trace(ordered, ev, p, set())
                    assert inventory == Counter(map(price, B))
                    counts['independent_complete_trace_controls'] += 1
            assert inventories[0] == inventories[1]
            assert sum(inventories[0].values()) == n
            serial = [[str(t), m] for t, m in sorted(inventories[0].items())]
            digest = hashlib.sha256(json.dumps(serial).encode()).hexdigest()
            intervals.append({'ratio_interval': [str(lo), 'infinity' if step == len(walls)-2 else str(hi)],
                              'sample_ratio': str(H), 'inventory_sha256': digest})
            counts['finite_entire_seed_price_intervals'] += 1
        fixtures.append({'cutoff': s, 'degree': a, 'interaction': r, 'characteristic': p,
                         'core_dimension': n, 'intervals': intervals,
                         'scope': 'Every open ordering interval for these parameters; whole boundary tiers follow from matching rank-birth prices by continuity. No uniform parameter conclusion.'})
        print('Complete finite core seed', s, a, r, p, len(intervals), flush=True)
    return dict(counts), fixtures


def all_field_witness(rows, evaluate, selected):
    M = [[evaluate(v).get(i, 0) for v in selected] for i in range(len(rows))]
    units, residual = smith_units(M)
    assert all(not x for row in residual for x in row)
    return {'rank_over_every_field': units, 'unit_invariants': units,
            'residual_zero': True}


def scope_controls():
    counts = defaultdict(int)
    s, a, r = 15, 2, 3
    E = core_ground(s, a)
    bounds = [13, 1, 0]+[-1]*5
    selected = [v for v in E if v[2] <= bounds[v[0]+v[1]]]
    assert len(selected) == 21
    signatures = [all_field_witness(*core_model(s, a, r, canonical), selected)
                  for canonical in (False, True)]
    assert [x['rank_over_every_field'] for x in signatures] == [21, 20]
    for j, c in ((2, 0), (0, 13)):
        assert (j, 0, c) in selected
    assert (1, 0, 2) not in selected
    # (1,0,2) is the stated strict convex combination of U^2,z^13,1.
    assert F(1, 2)*2 == 1 and F(2, 13)*13 == 2 and F(1, 2)+F(2, 13) < 1
    staircase = {'cutoff': s, 'degree': a, 'interaction': r,
                 'selected_columns': selected, 'regular': signatures[0], 'canonical': signatures[1],
                 'coherent_positive_prefix': False}
    counts['all_field_noncoherent_rank_counterexamples'] += 1
    s, w, T = 5, (F(9), F(9), F(4)), F(20)
    selected = [v for v in ground(s) if dot(v, w) <= T]
    pairs = (({(2, 0): 1, (0, 2): 1}, {(1, 1): 1}),
             ({(0, 2): 1}, {(0, 2): 1}))
    signatures = [all_field_witness(*form_model(s, 2, 3, Q, H), selected) for Q, H in pairs]
    assert [x['rank_over_every_field'] for x in signatures] == [15, 14]
    price = {'cutoff': s, 'degree': 2, 'interaction': 3,
             'weight': list(map(str, w)), 'price_cutoff': str(T),
             'selected_original_columns': selected,
             'regular': signatures[0], 'canonical': signatures[1], 'within_u_less_than_Z': False}
    counts['all_field_price_scope_counterexamples'] += 1
    return dict(counts), {'noncoherent_staircase': staircase, 'opposite_price_side': price}


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    blocks = block_controls()
    integral = integral_original_controls()
    forms, form_fixtures = form_controls()
    smith, smith_fixtures = smith_controls()
    intervals, interval_fixtures = entire_seed_controls()
    scope, scope_fixtures = scope_controls()
    sources = ('docs/higher_seed_binomial_reduction.md', 'scripts/run_higher_seed_checks.py',
               'docs/equal_degree_unit_fans.md', 'docs/middle_interaction_fans.md',
               'docs/signed_face_reduction.md', 'docs/higher_interaction_fans.md',
               'scripts/run_higher_interaction_checks.py',
               'scripts/run_hypersurface_truncation_checks.py',
               'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_socle_stability_checks.py')
    report = {
        'status': 'passed', 'command': 'python3 -B scripts/run_higher_seed_checks.py',
        'counts': {**blocks, **integral, **forms, **smith, **intervals, **scope},
        'original_form_fixtures': form_fixtures, 'finite_smith_matches': smith_fixtures,
        'entire_finite_seed_intervals': interval_fixtures, 'scope_counterexamples': scope_fixtures,
        'uniform_results': {
            'regular_binary_forms': 'For r>a>=1, coprime degree-a binary Q,H have the exact regular-core rank sum (B3) at complete tiers, all fields; characteristic dependence is retained.',
            'noncoprime_bracket': 'For all nonzero degree-a Q,H their seed ranks lie between the explicit canonical and regular binomial sums.',
            'original_integral_decomposition': 'For the coefficient-one quadratic pair the seed prefix lattices split into G_s, two G_(s-1), and G_(s-2), including ambient free lattices.',
            'finite_block_query': 'At exact rational prices: O_a(s^3*(1+s/r)^2) field operations, O_a(s^2) additional rational arithmetic/floor operations, and O_a((1+s/r)^2) field elements of scratch space; binomial generation and bit costs separate.',
            'sufficient_remaining_identity': '(B8), equality of regular/canonical core ranks at every coherent cutoff with 0<u<Z.',
            'identity_equivalence': 'Uniform original regular/canonical seed equality is equivalent to (B8), by inverting the finite homogeneous-channel convolution.',
            'conditional_consequence': 'For a=2, (B8) would transfer every remaining original largest-z cone and give a uniform full-family O(n) bound.',
            'arbitrary_staircase_equality': False,
            'arbitrary_price_equality': False},
        'uniform_core_identity_proved': False,
        'r_ge_3_full_signed_transfer_proved': False,
        'general_r_full_family_linear_bound_proved': False,
        'unrestricted_G3_improvement': False, 'novelty_established': False, 'goal_complete': False,
        'scope': 'Exact uniform binomial reductions and bounds, with finite rank/Smith matches only. No uniform core rank identity, full r>=3 fan count, unrestricted improvement or priority certificate.',
        'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/'results/research/higher_seed_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
