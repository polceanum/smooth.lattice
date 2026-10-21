#!/usr/bin/env python3
"""Exact pole transports, variance budgets and Frobenius fan controls."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import comb
from pathlib import Path
import hashlib
import json

from run_pole_order_checks import cleared_columns, polynomial_power, rank as polynomial_rank
from run_projection_tag_checks import fan, polygon_key
from run_quartic_hurwitz_checks import dense_rank, model as quartic_model, scalar
from run_quadric_cross_checks import apply, normalized, rank
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]


def frobenius(q, p):
    if q == 1:
        return True
    if not p:
        return False
    while q % p == 0:
        q //= p
    return q == 1


def model(f, p, q):
    """Use the independently checked square-pole columns to recover H=1/(t-1).

    H^(2r+1)=(t-1)*H^(2r+2). No root extraction or field division by two.
    """
    n, old, old_matrices = quartic_model(f, p)

    @lru_cache(None)
    def ev(v):
        a, b, c = v
        e = q*c
        if e % 2 == 0:
            return old((a, b, e//2))
        left, right = old((a+1, b, (e+1)//2)), old((a, b, (e+1)//2))
        return normalized({i: left.get(i, 0)-right.get(i, 0) for i in range(n)}, p)

    matrices = old_matrices[:2]+[tuple(ev((i, 0, 1)) for i in range(n))]
    return n, ev, matrices


def allocation(A, B, heights):
    D = {(0, 0, 0)}
    D.update((i, 0, 0) for i in range(1, A+1))
    D.update((0, i, 0) for i in range(1, B+1))
    D.update((max(j, 0), max(-j, 0), r)
             for j, c in heights.items() for r in range(1, c+1))
    return frozenset(D)


def parse(D, q, m):
    A = max((v[0] for v in D if not v[1] and not v[2]), default=0)
    B = max((v[1] for v in D if not v[0] and not v[2]), default=0)
    heights = {j: max((v[2] for v in D if v[:2] == (max(j, 0), max(-j, 0))), default=0)
               for j in range(-m, q-m)}
    assert D == allocation(A, B, heights)
    assert all(heights.get(j-1, 0) >= c for j, c in heights.items() if j > 0)
    assert all(heights.get(j+1, 0) >= c for j, c in heights.items() if j < 0)
    assert all(not c or (A >= j if j >= 0 else B >= -j) for j, c in heights.items())
    return A, B, heights


def profile(D):
    return tuple(sum(v[i] for v in D) for i in range(3))


def remainder(poly, f, p):
    v = list(poly)+[0]*max(0, len(f)-1-len(poly))
    n = len(f)-1
    for i in range(len(v)-1, n-1, -1):
        c = v[i]
        for j in range(n):
            v[i-n+j] = scalar(v[i-n+j]-c*f[j], p)
    return normalized(dict(enumerate(v[:n])), p)


def span_rank(points, all_columns, p):
    return polynomial_rank([all_columns[v] for v in points], p)


def transport_controls(counts):
    rows = []
    orders = ((1, 0), (2, 0), (2, 2), (3, 0), (3, 2), (3, 3),
              (4, 2), (5, 5), (8, 2), (9, 3))
    for q, p in orders:
        for m in range(q):
            a = q-m
            A, B = a+q, m+q+1
            for kind in range(3):
                heights = {j: (1 if kind == 0 else max(0, (2 if kind == 1 else 3)-abs(j)))
                           for j in range(-m, a)}
                D = allocation(A, B, heights)
                C = max(heights.values())
                three = ((A, 0, 0), (0, B, 0), (0, 0, C))
                source = cleared_columns(three, q, p)
                f = tuple(scalar(sum(row[i] for row in source.values()), p)
                          for i in range(len(next(iter(source.values())))))
                assert f[-1] == 1 and f[0] and scalar(sum(f), p)
                n, ev, _ = model(f, p, q)
                old_rank = rank([ev(v) for v in D], p)
                deficit = q*C-sum(heights.values())
                if deficit:
                    assert n >= len(D) and old_rank < len(D)
                    counts['genuinely_deficient_quotient_rank_fixtures'] += 1
                shifted = {j: c+1 for j, c in heights.items()}
                checks = []
                for axis in (0, 1):
                    if axis == 1 and not frobenius(q, p):
                        continue
                    E = allocation(A-q if axis == 0 else A, B-q if axis == 1 else B, shifted)
                    multiplied = frozenset((v[0]+(q if axis else 0), v[1], v[2]+1) for v in D)
                    columns = cleared_columns(E | multiplied, q, p)
                    rr = (span_rank(E, columns, p), span_rank(multiplied, columns, p),
                          span_rank(E | multiplied, columns, p))
                    assert rr == (len(D),)*3
                    assert rank([ev(v) for v in E], p) == old_rank
                    counts['independent_transported_rational_span_equalities'] += 1
                    counts['independent_transported_residue_rank_equalities'] += 1
                    checks.append({'axis': 'xy'[axis], 'rational_span_rank': rr[0],
                                   'residue_rank': old_rank})
                rows.append({'pole_order': q, 'characteristic': p, 'chart': m,
                             'axes': (A, B), 'heights': heights, 'monomials': len(D),
                             'quotient_length': n, 'modulus': f, 'checks': checks})
    f = (-1, 6, -10, 10, -5, 1)
    _, ev, _ = model(f, 0, 2)
    D = allocation(0, 3, {-1: 0, 0: 1})
    E = allocation(0, 1, {-1: 1, 0: 2})
    assert (rank([ev(v) for v in D], 0), rank([ev(v) for v in E], 0)) == (5, 4)
    counts['verified_non_Frobenius_transport_failures'] += 1
    return rows, {'pole_order': 2, 'characteristic': 0, 'modulus': f,
                  'before': sorted(D), 'after': sorted(E), 'ranks': (5, 4)}


def budget(D, w, q, m, f, p, ev, counts, cache):
    n = len(D)
    A, B, heights = parse(D, q, m)
    alpha, beta = Q(w[0], w[2]), Q(w[1], w[2])
    T, a = min(alpha, beta), q-m
    tau = {j: (j*alpha if j >= 0 else -j*beta) for j in heights}
    assert all(0 <= price <= q*T for price in tau.values())
    if n < q:
        return 'small_length', None
    cheap = int(beta < alpha)
    C, extra = divmod(n-q, q)
    A0, B0 = a-1+extra*(cheap == 0), m+extra*(cheap == 1)
    balanced0 = allocation(A0, B0, {j: C for j in heights})
    key = m, cheap
    if key not in cache:
        assert len(balanced0) == n and rank([ev(v) for v in balanced0], p) == n
        cache[key] = True
        counts['independent_universal_complete_pole_bases'] += 1
    maximum = max(dot((alpha, beta, 1), v) for v in balanced0)
    assert maximum <= 2*q*T+Q(n, q)
    assert max(dot((alpha, beta, 1), v) for v in D) <= maximum
    assert max(A, B) <= 2*q+Q(n, q)/T
    counts['independent_greedy_bottleneck_and_axis_budgets'] += 1
    if A < a-1 or B < m:
        if n >= 4*q:
            assert n <= 2*q+q*q*T
            assert T >= Q(n, 2*q*q) and max(A, B) <= 4*q
            mu = profile(D)
            assert max(mu[:2]) <= 4*q*n and mu[2] <= n*n
            counts['independent_missing_tag_collective_box_controls'] += 1
        return 'missing_tag', None
    if cheap and not frobenius(q, p):
        return 'outside_full_transport_hypotheses', None
    R, other, required = (A, B, a-1) if cheap == 0 else (B, A, m)
    M, S = n-1-other, sum(heights.values())
    assert S == M-R
    v = {j: Q(c)-Q(S, q) for j, c in heights.items()}
    P = sum(tau.values())
    energy = sum(x*x for x in v.values())/2+sum(tau[j]*v[j] for j in heights)

    def g(r):
        return Q(q*T+1, 2*q)*r*r+(T/2-Q(M, q)-Q(1, 2)-P/q)*r

    derivative = Q(q*T+1, q)*R+T/2-Q(M, q)-Q(1, 2)-P/q
    if R >= required+q:
        shifted = allocation(A-q if cheap == 0 else A, B-q if cheap else B,
                             {j: c+1 for j, c in heights.items()})
        assert rank([ev(e) for e in shifted], p) == n
        difference = dot((alpha, beta, 1), profile(shifted))-dot((alpha, beta, 1), profile(D))
        assert difference == g(R-q)-g(R) >= 0
        assert derivative <= (q*T+1)/2
        counts['independent_actual_coherent_transport_costs'] += 1
    else:
        assert R < 2*q and derivative <= 2*q*T
    C, extra = divmod(S, q)
    balanced = allocation(A+extra*(cheap == 0), B+extra*(cheap == 1),
                          {j: C for j in heights})
    assert len(balanced) == n and rank([ev(e) for e in balanced], p) == n
    increase = g(R+extra)-g(R)
    assert increase <= 4*q*q*T+q
    difference = dot((alpha, beta, 1), profile(balanced))-dot((alpha, beta, 1), profile(D))
    assert difference == increase-energy >= 0
    assert energy <= 4*q*q*T+q
    norm2 = sum(x*x for x in v.values())
    U = max(Q(1), T)
    assert norm2 <= 64*q**3*U*U
    X = sum(j*c for j, c in heights.items() if j > 0)
    Y = sum(-j*c for j, c in heights.items() if j < 0)
    assert abs(X-Q(S, q)*sum(j for j in heights if j > 0)) <= 8*q**3*U
    assert abs(Y-Q(S, q)*sum(-j for j in heights if j < 0)) <= 8*q**3*U
    counts['independent_exact_balancing_variance_and_moment_controls'] += 1
    return 'variance', {'axes': (A, B), 'heights': heights, 'minimum_price': str(T),
                        'variance_squared': str(norm2), 'shape_energy': str(energy),
                        'balanced_cost_increase': str(increase)}


def fixture(f, p, q, label, counts):
    n, ev, matrices = model(f, p, q)
    for j in range(n):
        assert apply(matrices[0], matrices[1][j], p) == {j: scalar(1, p)}
        v = matrices[2][j]
        for _ in range(q):
            x = apply(matrices[0], v, p)
            v = normalized({i: x.get(i, 0)-v.get(i, 0) for i in range(n)}, p)
        assert v == {j: scalar(1, p)}
        for i, k in combinations(range(3), 2):
            assert apply(matrices[i], matrices[k][j], p) == apply(matrices[k], matrices[i][j], p)
        counts['independent_quotient_multiplication_columns'] += 1
    actual = fan(n, ev, p, counts)
    grounds = tuple(hyperbolic_cross(n, 3))
    cleared = cleared_columns(grounds, q, p)
    # Multiplying by one common unit preserves the rank of every selected set.
    replay = {e: remainder(poly, f, p) for e, poly in cleared.items()}
    cache, scope, keys = {}, defaultdict(int), {}
    chart_counts = defaultdict(int)
    for D, poly in actual.items():
        assert rank([replay[v] for v in D], p) == n
        counts['independent_cleared_numerator_basis_rank_replays'] += 1
        covered = Q(0)
        for m in range(q):
            patch = clip(list(poly), (q-m, -m, 0), 1)
            patch = clip(patch, (-(q-m-1), m+1, 0), 1)
            if not area2(patch):
                continue
            covered += area2(patch)
            w = tuple(sum(v[i] for v in patch)/len(patch) for i in range(3))
            A, B, heights = parse(D, q, m)
            X = sum(j*c for j, c in heights.items() if j > 0)
            Y = sum(-j*c for j, c in heights.items() if j < 0)
            key = (m, A, B, X, Y)
            assert key not in keys or keys[key] == D
            keys[key] = D
            which, _ = budget(D, w, q, m, f, p, ev, counts, cache)
            scope[which] += 1
            chart_counts[m] += 1
            counts['independent_complete_chart_patches_and_profile_injections'] += 1
        assert covered == area2(poly)
        counts['independent_complete_carrier_chart_area_covers'] += 1
    digest = hashlib.sha256(json.dumps(
        [(sorted(D), polygon_key(P)) for D, P in sorted(actual.items(), key=lambda v: sorted(v[0]))],
        default=str).encode()).hexdigest()
    return {'label': label, 'pole_order': q, 'characteristic': p, 'length': n,
            'modulus': f, 'initial_ideals': len(actual), 'chart_patch_counts': dict(chart_counts),
            'budget_scopes': dict(scope), 'complete_fan_sha256': digest}


def large_controls(counts):
    rows = []
    for q, p in ((3, 3), (4, 2), (5, 5), (8, 2), (9, 3)):
        n = max(32, 6*q)
        f = [1, -1]+[0]*(n-2)+[1]
        _, ev, _ = model(f, p, q)
        E, cache = tuple(hyperbolic_cross(n, 3)), {}
        prices = (Q(1, 16), Q(1), Q(n, 2*q*q), Q(n, 4), Q(2*n))
        for m in range(q):
            ratio = Q(2*m+1, 2)
            other = q-ratio
            low = min(ratio, other)
            for T in prices:
                w = (ratio*T/low, other*T/low, Q(1))
                D, _ = field_basis(sorted(E, key=lambda e: (dot(w, e), e)), ev, p)
                D = frozenset(D)
                assert len(D) == n
                which, record = budget(D, w, q, m, f, p, ev, counts, cache)
                rows.append({'length': n, 'pole_order': q, 'characteristic': p,
                             'chart': m, 'weight': list(map(str, w)), 'scope': which,
                             'variance_witness': record})
                counts['independent_large_length_full_ground_greedy_scans'] += 1
    return rows


def parity_envelope(points):
    path = []
    for Y, Z, b, h in points:
        low, high, possible = Q(0), None, True
        for U, V, bb, hh in points:
            if (bb, hh) == (b, h):
                continue
            slope, limit = Y-U, V-Z
            if slope == 0:
                if limit <= 0:
                    possible = False
                    break
            elif slope > 0:
                bound = Q(limit, slope)
                high = bound if high is None else min(high, bound)
            else:
                low = max(low, Q(limit, slope))
        if possible and (high is None or low < high):
            path.append((low, b, h))
    path.sort()
    assert all(b0 >= b1 and h0 <= h1 for (_, b0, h0), (_, b1, h1) in zip(path, path[1:]))
    return len(path)


def binary_hurwitz_controls(counts):
    rows = []
    for n in range(2, 9):
        for middle in product(range(2), repeat=n-1):
            f = (1,)+middle+(1,)
            if sum(f) % 2 == 0:
                continue
            g = [sum(f[j]*comb(j, k) for j in range(k, n+1)) % 2 for k in range(2*n)]
            matrices, feasible = [], []
            for b in range(n+2):
                H = [[g[2*i+1-j] if 0 <= 2*i+1-j < len(g) else 0
                      for j in range(n)] for i in range(n)]
                matrices.append(H)
                feasible.append([dense_rank([r[:h] for r in H[:h]], h, 2) == h for h in range(n)])
                if b >= 2:
                    assert all(H[i][j] == sum(matrices[b-2][r][j] for r in range(i+1)) % 2
                               for i in range(n) for j in range(n))
                    assert feasible[b] == feasible[b-2]
                    counts['independent_even_unit_Hurwitz_matrix_and_minor_periods'] += 1
                nxt = []
                for k, coefficient in enumerate(g):
                    nxt.append((coefficient-(nxt[-1] if k else 0)) % 2)
                g = nxt
            envelopes = 0
            for a in range(n):
                for parity in (0, 1):
                    points = []
                    for b in range(n-a):
                        total = n-1-a-b
                        for h in range(total % 2, total+1, 2):
                            if h % 2 != parity or not feasible[b][h]:
                                continue
                            c, d = (total+h)//2, (total-h)//2
                            if d and b == 0:
                                continue
                            Y = b*(b+1)//2+d
                            Z = c*(c+1)//2+d*(d+1)//2
                            points.append((Y, Z, b, h))
                    length = parity_envelope(points)
                    assert length <= 2*(n-1-a)+1
                    envelopes += 1
                    counts['independent_entire_second_chart_parity_envelopes'] += 1
            rows.append({'length': n, 'binary_modulus': f, 'checked_second_chart_envelopes': envelopes})
    return rows


def main():
    counts, cases = defaultdict(int), []
    transport, failure = transport_controls(counts)
    for q, p in ((1, 0), (2, 2), (3, 3), (4, 2), (5, 5), (8, 2), (9, 3),
                 (2, 0), (3, 0), (3, 2)):
        for n in (4, 6, 8):
            cases.append(([1, -1]+[0]*(n-2)+[1], p, q, 'unit_modulus'))
        if p and p != 2:
            n = 7
            cases.append(([comb(n, i)*(-2)**(n-i) for i in range(n+1)], p, q, 'jet_at_2'))
        elif p == 2:
            f = polynomial_power((1, 1, 1), 4, p)
            cases.append((f, p, q, 'repeated_irreducible_factor'))
    cases += [([1, -1]+[0]*14+[1], 3, 3, 'larger_complete_fan'),
              ([1, -1]+[0]*18+[1], 2, 4, 'larger_complete_fan')]
    fans = [fixture(f, p, q, label, counts) for f, p, q, label in cases]
    large = large_controls(counts)
    binary = binary_hurwitz_controls(counts)
    sources = ['docs/frobenius_pole_transport.md', 'docs/quartic_hurwitz_fans.md',
               'scripts/run_frobenius_pole_checks.py',
               'scripts/run_pole_order_checks.py', 'scripts/run_quartic_hurwitz_checks.py',
               'scripts/run_projection_tag_checks.py', 'scripts/run_quadric_cross_checks.py',
               'scripts/run_signed_face_checks.py', 'scripts/run_structural_complexity_checks.py']
    report = {'command': 'python3 -B scripts/run_frobenius_pole_checks.py',
              'counts': dict(sorted(counts.items())), 'rational_and_residue_transports': transport,
              'non_Frobenius_failure': failure, 'complete_fans': fans,
              'large_length_budget_controls': large, 'binary_Hurwitz_periodicity': binary,
              'uniform_full_fan_bound': 'C*q^9*n^2*(1+log(n+1)), q=p^e in characteristic p',
              'arbitrary_q_scope': 'Same estimate on the x-cheaper part; other part not established.',
              'binary_quartic_sharp_orders': {'fan': 'n^2', 'wall_cost': 'n^3', 'basis_union': 'n^2'},
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              'pole_order_power_sharp': False, 'unrestricted_G3_improved': False,
              'polynomial_lower_separation_proved': False, 'novelty_established': False,
              'goal_complete': False}
    (ROOT/'results/research/frobenius_pole_checks.json').write_text(json.dumps(report, indent=2, default=str)+'\n')
    print(json.dumps({'counts': report['counts'], 'complete_fans': len(fans),
                      'whole_fan_cones': sum(v['initial_ideals'] for v in fans),
                      'pole_orders': sorted({v['pole_order'] for v in fans}),
                      'binary_moduli': len(binary)}, indent=2))


if __name__ == '__main__':
    main()
