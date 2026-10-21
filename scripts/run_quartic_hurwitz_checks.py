#!/usr/bin/env python3
"""Exact all-field quartic-curve ranks, coherent budgets and lower wall costs."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import comb
from pathlib import Path
import hashlib
import json

from run_curve_carrier_shadow_checks import reduce_polynomial
from run_interleaving_research import minimal_excluded_border
from run_primitive_wall_checks import canonical, threshold
from run_projection_tag_checks import fan, polygon_key
from run_quadric_cross_checks import apply, normalized, rank
from run_signed_face_checks import field_basis, field_reduce, normals
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]


def scalar(c, p):
    return c % p if p else Q(c)


def inverse(c, p):
    return pow(c % p, -1, p) if p else 1/Q(c)


def dense_rank(rows, width, p):
    """A separate row elimination for Taylor matrices, not residue columns."""
    A = [[scalar(c, p) for c in row] for row in rows]
    pivot = 0
    for j in range(width):
        i = next((i for i in range(pivot, len(A)) if A[i][j]), None)
        if i is None:
            continue
        A[pivot], A[i] = A[i], A[pivot]
        inv = inverse(A[pivot][j], p)
        A[pivot] = [scalar(c*inv, p) for c in A[pivot]]
        for i in range(pivot+1, len(A)):
            c = A[i][j]
            if c:
                A[i] = [scalar(a-c*b, p) for a, b in zip(A[i], A[pivot])]
        pivot += 1
        if pivot == len(A):
            break
    return pivot


def model(f, p):
    """Polynomial remainders build actual commuting quotient multiplication."""
    f = tuple(scalar(c, p) for c in f)
    n, f1 = len(f)-1, scalar(sum(f), p)
    assert n > 0 and f[-1] == 1 and f[0] and f1
    unit = (scalar(1, p),)+(0,)*(n-1)

    def mul(a, b):
        v = [0]*(max(n, len(a)+len(b)-1))
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                v[i+j] = scalar(v[i+j]+x*y, p)
        for i in range(len(v)-1, n-1, -1):
            c = v[i]
            for j in range(n):
                v[i-n+j] = scalar(v[i-n+j]-c*f[j], p)
        return tuple(v[:n])

    def power(a, k):
        out = unit
        while k:
            if k & 1:
                out = mul(out, a)
            a = mul(a, a)
            k //= 2
        return out

    X = mul(unit, (0, 1))
    Y = tuple(scalar(-f[i+1]*inverse(f[0], p), p) for i in range(n))
    H = tuple(scalar(-sum(f[i+1:])*inverse(f1, p), p) for i in range(n))
    Z = mul(H, H)
    images = (X, Y, Z)

    @lru_cache(None)
    def ev(a):
        out = unit
        for image, k in zip(images, a):
            out = mul(out, power(image, k))
        return normalized(dict(enumerate(out)), p)

    matrices = [tuple(normalized(dict(enumerate(mul(
        image, tuple(int(i == j) for i in range(n))))), p) for j in range(n))
        for image in images]
    return n, ev, matrices


def allocation(a, b, c, d, tag):
    B = {(0, 0, 0)}
    B.update((i, 0, 0) for i in range(1, a+1))
    B.update((0, j, 0) for j in range(1, b+1))
    B.update((0, 0, k) for k in range(1, c+1))
    B.update((int(tag == 0), int(tag == 1), k) for k in range(1, d+1))
    return frozenset(B)


def valid(a, b, c, d, tag):
    return c >= d and (not d or (a, b)[tag] >= 1)


def allocations(n, tag, bounded=False):
    for b in range(n):
        M = n-1-b
        for a in range(M+1):
            end = M-a
            if bounded and a >= 3 and not (tag == 1 and b == 0):
                end = min(end, 4*n//(a-1))
            for h in range((M-a) % 2, end+1, 2):
                c, d = (M-a+h)//2, (M-a-h)//2
                if valid(a, b, c, d, tag):
                    yield a, b, c, d


def profile(B):
    return tuple(sum(a[i] for a in B) for i in range(3))


def rank_oracle(f, p, counts):
    n = len(f)-1

    @lru_cache(None)
    def series(b, length):
        shifted = [sum(f[j]*comb(j, k) for j in range(k, n+1)) for k in range(length)]
        denom = [(-1)**k*comb(b+k-1, k) if b else int(k == 0) for k in range(length)]
        return tuple(scalar(sum(shifted[j]*denom[k-j] for j in range(k+1)), p)
                     for k in range(length))

    @lru_cache(None)
    def matrix_rank(b, h, ell):
        g = series(b, 2*h)
        rows = [[g[2*i+1-j] if 0 <= 2*i+1-j < len(g) else 0
                 for j in range(ell)] for i in range(h)]
        out = dense_rank(rows, ell, p)
        counts['independent_Taylor_matrix_eliminations'] += 1
        return out

    def predicted(a, b, c, d):
        m, h = 1+a+b+c+d, c-d
        ell = max(0, m+h-n)
        return m-ell+matrix_rank(b, h, ell)

    return predicted


def carrier_controls(counts):
    one = (0, 0, 0)
    xy = {(1, 1, 0): 1, one: -1}
    x2z = {(2, 0, 1): 1, (1, 0, 1): -2, (0, 0, 1): 1, one: -1}
    xzyz = {(1, 0, 1): 1, (0, 1, 1): 1, (0, 0, 1): -2, (0, 1, 0): -1}
    y2z = {(0, 2, 1): 1, (0, 1, 1): -2, (0, 0, 1): 1, (0, 2, 0): -1}
    for tag, w, polys, heads in [
        (0, (1, 2, 1), (xy, x2z, xzyz), ((1, 1, 0), (2, 0, 1), (0, 1, 1))),
        (1, (2, 1, 1), (xy, xzyz, y2z), ((1, 1, 0), (1, 0, 1), (0, 2, 1)))]:
        assert tuple(max(f, key=lambda a: (dot(w, a), a)) for f in polys) == heads
        for f in (xy, x2z):
            assert not reduce_polynomial(f, polys, heads, w)
        for i, j in combinations(range(3), 2):
            common = tuple(max(a, b) for a, b in zip(heads[i], heads[j]))
            S = defaultdict(int)
            for k, sign in ((i, 1), (j, -1)):
                shift = tuple(a-b for a, b in zip(common, heads[k]))
                for a, c in polys[k].items():
                    S[tuple(x+y for x, y in zip(a, shift))] += sign*c
            assert not reduce_polynomial({a: c for a, c in S.items() if c}, polys, heads, w)
            counts['integral_carrier_S_pairs'] += 1
    return True


def fixture(f, p, label, counts, rectangular=False):
    f = tuple(scalar(c, p) for c in f)
    n, ev, matrices = model(f, p)
    predicted = rank_oracle(f, p, counts)
    for j in range(n):
        assert apply(matrices[0], matrices[1][j], p) == {j: scalar(1, p)}
        # (x-1)^2 z = 1, verified as actual endomorphisms.
        z = matrices[2][j]
        xz = apply(matrices[0], z, p)
        x2z = apply(matrices[0], xz, p)
        assert normalized({i: x2z.get(i, 0)-2*xz.get(i, 0)+z.get(i, 0)
                           for i in range(n)}, p) == {j: scalar(1, p)}
        for a, b in combinations(range(3), 2):
            assert apply(matrices[a], matrices[b][j], p) == apply(matrices[b], matrices[a][j], p)
        counts['independent_quotient_multiplication_columns'] += 1
    if rectangular:
        for a, b, c, d in product(range(4), repeat=4):
            for tag in (0, 1):
                if not valid(a, b, c, d, tag):
                    continue
                B = allocation(a, b, c, d, tag)
                assert rank([ev(v) for v in B], p) == predicted(a, b, c, d)
                counts['independent_rectangular_rank_equalities'] += 1
    candidates, lookup, shifts = [], [], 0
    for tag in (0, 1):
        local, dictionary = [], {}
        for lengths in allocations(n, tag):
            a, b, c, d = lengths
            feasible = predicted(a, b, c, d) == n
            B = allocation(*lengths, tag)
            if n <= 6:
                assert rank([ev(v) for v in B], p) == predicted(a, b, c, d)
                counts['independent_all_allocation_rank_equalities'] += 1
            if a >= 3 and not (tag == 1 and b == 0):
                assert predicted(a-2, b, c+1, d+1) == predicted(a, b, c, d)
                shifts += 1
            if feasible:
                local.append((B, profile(B), lengths))
                dictionary[B] = lengths
        candidates.append(local)
        lookup.append(dictionary)
    for b in range(n):
        for parity in (0, 1):
            points = [(value[0], value[2], lengths)
                      for _, value, lengths in candidates[0]
                      if lengths[1] == b and (lengths[2]-lengths[3]) % 2 == parity]
            path = []
            for X, Z, lengths in points:
                low, high, possible = Q(0), None, True
                for U, V, other in points:
                    if other == lengths:
                        continue
                    slope, limit = X-U, V-Z
                    if not slope:
                        if limit <= 0:
                            possible = False
                            break
                    elif slope > 0:
                        bound = Q(limit, slope)
                        high = bound if high is None else min(high, bound)
                    else:
                        low = max(low, Q(limit, slope))
                if possible and (high is None or low < high):
                    path.append((low, lengths[0], lengths[2]-lengths[3]))
            path.sort()
            assert len(path) <= 2*(n-1-b)+1
            assert all(a0 >= a1 and h0 <= h1 for (_, a0, h0), (_, a1, h1)
                       in zip(path, path[1:]))
            counts['independent_complete_fixed_b_parity_envelopes'] += 1
    actual = fan(n, ev, p, counts)
    differences = []
    for B, poly in actual.items():
        reduced_normals, relations = normals(tuple(sorted(B)), ev, p)
        assert len(relations) <= 7
        center = tuple(sum(v[i] for v in poly)/len(poly) for i in range(3))
        assert all(dot(q, center) > 0 for q in reduced_normals)
        for tag in (0, 1):
            patch = clip(list(poly), (1, -1, 0), -1 if tag == 0 else 1)
            if not area2(patch):
                continue
            assert B in lookup[tag]
            a, b, c, d = lookup[tag][B]
            assert a <= 2 or (tag == 1 and b == 0) or (a-1)*(c-d) < 4*n
            base = profile(B)
            # Full patches, including their edges, are checked against every
            # feasible allocation, not just a finite set of sampled weights.
            for _, value, _ in candidates[tag]:
                delta = tuple(x-y for x, y in zip(value, base))
                assert all(dot(delta, v) >= 0 for v in patch)
                counts['independent_entire_patch_profile_comparisons'] += 1
            counts['independent_entire_quartic_chart_patches'] += 1
            w = tuple(sum(v[i] for v in patch)/len(patch) for i in range(3))
            ground = [a for a in hyperbolic_cross(n, 3) if
                      not (a[0] and a[1]) and
                      (not a[2] or (a[0] < (2 if tag == 0 else 1) and
                                    a[1] < (1 if tag == 0 else 2)))]
            generic = frozenset(sorted(ground, key=lambda a: (dot(w, a), a))[:n])
            if generic != B:
                generic_rank = rank([ev(a) for a in generic], p)
                assert generic_rank < n
                if len(differences) < 2:
                    differences.append({'tag': 'xy'[tag], 'weight': list(map(str, w)),
                                        'actual_lengths': (a, b, c, d),
                                        'generic_basis': sorted(generic), 'generic_rank': generic_rank})
                counts['verified_failures_of_generic_rank_truncation'] += 1
    digest = hashlib.sha256(json.dumps(
        [(sorted(B), polygon_key(poly)) for B, poly in sorted(actual.items(), key=lambda v: sorted(v[0]))],
        default=str).encode()).hexdigest()
    return {'label': label, 'length': n, 'characteristic': p, 'modulus': list(f),
            'initial_ideals': len(actual), 'complete_fan_sha256': digest,
            'feasible_allocation_counts': list(map(len, candidates)),
            'feasibility_preserving_shifts': shifts, 'generic_rule_failure_witnesses': differences}


def lower_walls(n, p, counts):
    f = [1, -1]+[0]*(n-2)+[1]
    n, ev, _ = model(f, p)
    rows = []
    for a in range((2*n+4)//5, n//2+1):
        for b in range((n+9)//10, n//5+1):
            if (n-1-a-b) % 2:
                continue
            c = (n-1-a-b)//2
            if not (0 < c < a-1 and b < a-1):
                continue
            alpha = Q(2*c, 2*a-1)
            beta = (alpha+c)/(b+1)
            w = (alpha, beta, Q(1))
            B = allocation(a, b, c, c, 0)
            tip, head = (1, 0, c), (0, b+1, 0)
            alternate = B-{tip}|{head}
            assert alpha < beta and alpha < 1
            E = [e for e in hyperbolic_cross(n, 3) if
                 not (e[0] and e[1]) and (not e[2] or (e[0] < 2 and e[1] < 1))]
            assert {e for e in E if dot(w, e) < alpha+c} == B-{tip}
            assert {e for e in E if dot(w, e) == alpha+c} == {head, tip}
            selected, pivots = field_basis(tuple(sorted(B)), ev, p)
            assert len(selected) == n
            rem, co = field_reduce(ev(head), pivots, p)
            assert not rem
            condition = scalar((n-1)-(b+1), p) != 0
            assert bool(co.get(tip, 0)) == condition
            assert rank([ev(v) for v in alternate], p) == n-int(not condition)
            nn, _ = normals(tuple(sorted(B)), ev, p)
            active = {canonical(q) for q in nn if dot(q, w) == 0}
            q = (1, -(b+1), c)
            assert all(dot(q0, w) >= 0 for q0 in nn)
            assert active == ({q} if condition else set())
            if condition:
                assert threshold(q) == n-a+1
                assert len(minimal_excluded_border(tuple(B))) <= 7
                assert head in minimal_excluded_border(tuple(B))
                strict = [dot(q0, w)/(4*(1+abs(q0[1]))) for q0 in nn if dot(q0, w)]
                eps = min(strict+[beta/4])
                plus = (alpha, beta+eps, Q(1))
                minus = (alpha, beta-eps, Q(1))
                an, _ = normals(tuple(sorted(alternate)), ev, p)
                assert all(dot(q0, plus) > 0 for q0 in nn)
                assert all(dot(q0, minus) > 0 for q0 in an)
                counts['independent_strict_lower_wall_sides'] += 2
            rows.append({'a': a, 'b': b, 'c': c, 'weight': list(map(str, w)),
                         'actual_wall': condition, 'direction': q,
                         'primitive_cost': threshold(q) if condition else None})
            counts['independent_lower_wall_and_reduced_polynomial_checks'] += 1
    assert rows and any(v['actual_wall'] for v in rows)
    return {'length': n, 'characteristic': p, 'modulus': f, 'witnesses': rows,
            'distinct_forced_directions': sum(v['actual_wall'] for v in rows),
            'forced_total_cost': sum(v['primitive_cost'] or 0 for v in rows)}


def inventory_controls():
    rows = []
    for n in (16, 32, 64, 128, 256):
        raw = sum(sum(1 for _ in allocations(n, tag)) for tag in (0, 1))
        kept = sum(sum(1 for _ in allocations(n, tag, True)) for tag in (0, 1))
        upper = 2*n*(5*n+4*n*(1+(n-1).bit_length()))+n
        assert kept <= raw and kept <= upper
        if n >= 32:
            assert kept < raw
        rows.append({'length': n, 'all_downset_allocations': raw,
                     'hyperbolic_candidate_allocations': kept, 'explicit_harmonic_upper': upper})
    return rows


def large_budget_controls(counts):
    """Test the nontrivial budget using full monomial-column greedy scans."""
    rows = []
    weights = ((1, 2, 1), (2, 1, 1), (10, 20, 1), (20, 10, 1),
               (50, 60, 1), (60, 50, 1), (1000, 1, 10), (1, 1000, 10))
    for n, p in ((32, 2), (32, 101), (64, 2), (64, 101)):
        f = [1, -1]+[0]*(n-2)+[1]
        _, ev, _ = model(f, p)
        E = tuple(hyperbolic_cross(n, 3))
        for w in weights:
            B, _ = field_basis(sorted(E, key=lambda e: (dot(w, e), e)), ev, p)
            B = frozenset(B)
            tag = int(w[1] < w[0])
            a = max((v[0] for v in B if not v[1] and not v[2]), default=0)
            b = max((v[1] for v in B if not v[0] and not v[2]), default=0)
            c = max((v[2] for v in B if not v[0] and not v[1]), default=0)
            d = max((v[2] for v in B if v[:2] == ((1, 0) if tag == 0 else (0, 1))), default=0)
            assert B == allocation(a, b, c, d, tag) and len(B) == n
            h, M = c-d, n-1-b
            t, u = Q(w[0], w[2]), Q(w[1], w[2])
            r = min(t, u)

            def phi(aa, hh):
                return (2*t+1)*aa*aa/4+(t-r-M-1)*aa/2+Q(hh*hh, 4)-r*hh/2

            if a >= 3 and not (tag == 1 and b == 0):
                shifted = allocation(a-2, b, c+1, d+1, tag)
                assert rank([ev(v) for v in shifted], p) == n
                shift_cost = dot(w, profile(shifted))-dot(w, profile(B))
                assert shift_cost == w[2]*(phi(a-2, h)-phi(a, h)) >= 0
                assert (2*t+1)*a <= M+2+t+r
                assert 2*t*(a-1) <= n+1
                aa = a+h % 2
                cc = (M-aa)//2
                balanced = allocation(aa, b, cc, cc, tag)
                assert rank([ev(v) for v in balanced], p) == n
                balanced_cost = dot(w, profile(balanced))-dot(w, profile(B))
                assert balanced_cost == w[2]*(phi(aa, 0)-phi(a, h)) >= 0
                assert h*h-2*r*h <= 6*t+3 and h < 2*t+3 and (a-1)*h < 4*n
                counts['independent_large_length_transport_and_balancing_checks'] += 1
            rows.append({'length': n, 'characteristic': p, 'weight': w,
                         'tag': 'xy'[tag], 'actual_lengths': (a, b, c, d),
                         'budget_product': (a-1)*h, 'budget_limit': 4*n})
            counts['independent_large_length_full_ground_greedy_scans'] += 1
    return rows


def main():
    counts, cases = defaultdict(int), []
    carrier_controls(counts)
    cases.extend([([-2, 1], 0, 'length_one', True),
                  ([-2, 1], 3, 'length_one_finite_field', True),
                  ([1, 1, 1], 2, 'length_two_binary', True)])
    for n in range(3, 11):
        cases.append(([comb(n, i)*(-2)**(n-i) for i in range(n+1)], 0, 'jet_at_2', n <= 6))
        cases.append(([1, -1]+[0]*(n-2)+[1], 0, 'tn_minus_t_plus_1', n <= 6))
    for n in (5, 8):
        for b in (0, 1, 3):
            f = [comb(n, i)*(-1)**(n-i) for i in range(n+1)]
            f[0] += 1-b
            f[1] += b
            if f[0]:
                cases.append((f, 0, 'prescribed_first_Hurwitz_zero', True))
    for n in range(3, 9):
        for middle in product(range(2), repeat=n-1):
            f = (1,)+middle+(1,)
            if sum(f) % 2:
                cases.append((f, 2, 'all_binary_unit_moduli', n <= 4))
    for p in (3, 101):
        for n in (4, 7):
            cases.append(([comb(n, i)*(-2)**(n-i) for i in range(n+1)],
                          p, 'repeated_point_finite_characteristic', n == 4))
    rows = []
    for f, p, label, rectangle in cases:
        rows.append(fixture(f, p, label, counts, rectangle))
    lower = [lower_walls(n, p, counts) for n, p in ((24, 0), (32, 2), (32, 3), (32, 101))]
    inventory = inventory_controls()
    large = large_budget_controls(counts)
    sources = ['docs/quartic_hurwitz_fans.md', 'scripts/run_quartic_hurwitz_checks.py',
               'scripts/run_curve_carrier_shadow_checks.py', 'scripts/run_interleaving_research.py',
               'scripts/run_primitive_wall_checks.py', 'scripts/run_projection_tag_checks.py',
               'scripts/run_quadric_cross_checks.py', 'scripts/run_signed_face_checks.py',
               'scripts/run_structural_complexity_checks.py']
    report = {'command': 'python3 -B scripts/run_quartic_hurwitz_checks.py',
              'scope': 'All-field uniform proof for every quotient of one fixed quartic curve. '
                       'Finite controls support the exact rank, transport, candidate and lower-wall arguments.',
              'counts': dict(sorted(counts.items())), 'complete_fans': rows,
              'forced_wall_witnesses': lower, 'candidate_inventory': inventory,
              'large_length_budget_controls': large,
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              'class_upper_fan_order': 'n^2 log(n+1)', 'class_lower_fan_order': 'n^2',
              'class_upper_wall_cost_order': 'n^3 log(n+1)', 'class_lower_wall_cost_order': 'n^3',
              'logarithm_sharpness_established': False, 'unrestricted_G3_improved': False,
              'polynomial_lower_separation_proved': False, 'novelty_established': False,
              'goal_complete': False}
    destination = ROOT/'results/research/quartic_hurwitz_checks.json'
    destination.write_text(json.dumps(report, indent=2, default=str)+'\n')
    print(json.dumps({'counts': report['counts'], 'complete_fans': len(rows),
                      'whole_fan_cones': sum(v['initial_ideals'] for v in rows),
                      'candidate_inventory': inventory}, indent=2))


if __name__ == '__main__':
    main()
