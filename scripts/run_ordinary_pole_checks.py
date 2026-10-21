#!/usr/bin/env python3
"""Exact ordinary-pole basis completions and field-independent fan budgets."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import comb
from pathlib import Path
import hashlib
import json

from run_frobenius_pole_checks import allocation, model, parse, profile, remainder
from run_pole_order_checks import cleared_columns, polynomial_power, rank as polynomial_rank, solver
from run_projection_tag_checks import fan, polygon_key
from run_quadric_cross_checks import rank
from run_quartic_hurwitz_checks import scalar, inverse
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]


def rational_model(f, p, numerator, denominator):
    """Build 1/Q by solving multiplication by Q, without restrictions at t=1."""
    f, P, Qpoly = (tuple(scalar(c, p) for c in v)
                   for v in (f, numerator, denominator))
    n, q = len(f)-1, len(Qpoly)-1
    assert n > 0 and q > 0 and f[-1] == Qpoly[-1] == 1
    assert f[0] and Qpoly[0] and len(P) <= q+1
    pq = [remainder((0,)*j+P, Qpoly, p) for j in range(q)]
    assert rank(pq, p) == q, 'P and Q must be relatively prime'
    one = (scalar(1, p),)+(0,)*(n-1)

    def mul(left, right):
        out = [0]*max(n, len(left)+len(right)-1)
        for i, x in enumerate(left):
            for j, y in enumerate(right):
                out[i+j] = scalar(out[i+j]+x*y, p)
        for i in range(len(out)-1, n-1, -1):
            for j in range(n):
                out[i-n+j] = scalar(out[i-n+j]-out[i]*f[j], p)
        return tuple(out[:n])

    X = mul(one, (0, 1))
    Y = tuple(scalar(-f[i+1]*inverse(f[0], p), p) for i in range(n))
    Qmod, Pmod = mul(one, Qpoly), mul(one, P)
    columns = [mul(Qmod, (0,)*j+(1,)) for j in range(n)]
    H, in_span = solver(columns, p)(one)
    H = tuple(H)
    assert in_span and mul(Qmod, H) == one and mul(X, Y) == one
    Z = mul(Pmod, H)
    assert mul(Qmod, Z) == Pmod
    images, units = (X, Y, Z), (H, mul(mul(one, (0,)*q+(1,)), H))

    @lru_cache(None)
    def power(image, exponent):
        out = one
        while exponent:
            if exponent & 1:
                out = mul(out, image)
            image = mul(image, image)
            exponent //= 2
        return out

    @lru_cache(None)
    def ev(v):
        out = one
        for image, exponent in zip(images, v):
            out = mul(out, power(image, exponent))
        return {i: c for i, c in enumerate(out) if c}

    def unit_apply(vector, axis):
        out = mul(tuple(vector.get(i, 0) for i in range(n)), units[axis])
        return {i: c for i, c in enumerate(out) if c}

    matrices = [tuple({k: c for k, c in enumerate(mul(image, (0,)*j+(1,))) if c}
                      for j in range(n)) for image in images]
    return n, ev, matrices, unit_apply


def downset(D):
    return all(tuple(c-(i == j) for j, c in enumerate(v)) in D
               for v in D for i in range(3) if v[i])


def completion(D, w, q, m, ev, p, counts, axis=1, unit=None):
    A, B, heights = parse(D, q, m)
    n = len(D)
    assert A >= q-m-1+q*(axis == 0) and B >= m+q*(axis == 1)
    ell, top = min(heights.values()), max(heights.values())
    excess = sum(c-ell for c in heights.values())
    newA, newB = A-q*(axis == 0), B-q*(axis == 1)
    W = allocation(newA, newB, {j: ell+1 for j in heights})
    old_core = allocation(A, B, {j: ell for j in heights})
    assert len(W) == len(old_core)
    assert rank([ev(v) for v in W], p) == len(W)
    moved = ([unit(ev(v), axis) for v in D-old_core] if unit else
             [ev((v[0]+axis*q, v[1], v[2]+1)) for v in D-old_core])
    assert rank([ev(v) for v in W]+moved, p) == n
    pool = allocation(newA, newB, {j: top+1 for j in heights})
    order = tuple(sorted(W))+tuple(sorted(pool-W, key=lambda v: (dot(w, v), v)))
    E, _ = field_basis(order, ev, p)
    E = frozenset(E)
    assert len(E) == n and W <= E
    literal = allocation(newA, newB, {j: c+1 for j, c in heights.items()})
    T = min(Q(w[0], w[2]), Q(w[1], w[2]))
    literal_cost, completed_cost = dot(w, profile(literal)), dot(w, profile(E))
    upper = literal_cost+w[2]*q*T*excess
    assert completed_cost <= upper
    for r in range(ell+1, top+1):
        bound = w[2]*(r+1+q*T)
        prefix = W | {v for v in pool if dot(w, v) <= bound}
        need = len(W)+sum(max(0, min(c, r)-ell) for c in heights.values())
        assert rank([ev(v) for v in prefix], p) >= need
        counts['independent_completion_budget_prefix_ranks'] += 1
    literal_rank = rank([ev(v) for v in literal], p)
    if literal_rank < n:
        counts['dependent_literal_shifts_repaired_by_completion'] += 1
        counts['dependent_literal_shifts_axis_'+str(axis)] += 1
    if not downset(E):
        counts['genuine_non_downset_vector_completions'] += 1
    counts['independent_forced_core_basis_completions'] += 1
    counts['independent_forced_completions_axis_'+str(axis)] += 1
    return E, {'axis': axis, 'literal_rank': literal_rank, 'completed_rank': n, 'excess_columns': excess,
               'literal_cost': str(literal_cost), 'completion_cost': str(completed_cost),
               'cost_upper': str(upper), 'completion_is_downset': downset(E)}


def rational_span_controls(counts):
    rows = []
    for q, p in ((1, 0), (2, 0), (3, 0), (3, 2), (4, 3), (5, 2), (6, 5), (7, 3)):
        for m in range(q):
            A, B = 2*q-m, m+q+1
            for kind in (1, 2):
                heights = {j: max(0, kind+1-abs(j)) for j in range(-m, q-m)}
                ell = min(heights.values())
                D0 = allocation(A, B, {j: ell for j in heights})
                W = allocation(A, B-q, {j: ell+1 for j in heights})
                moved = frozenset(tuple(x+y for x, y in zip(v, (q, 0, 1))) for v in D0)
                all_columns = cleared_columns(W | moved, q, p)
                r0 = polynomial_rank([all_columns[v] for v in W], p)
                assert r0 == len(W)
                assert polynomial_rank([all_columns[v] for v in W | moved], p) == r0
                counts['independent_complete_core_rational_span_equalities'] += 1
                for r in range(ell+1, max(heights.values())+1):
                    core = allocation(q-m-1, m, {j: r+1 for j in heights})
                    moved = {(max(j, 0)+q, max(-j, 0), r+1)
                             for j, c in heights.items() if c >= r}
                    cols = cleared_columns(core | moved, q, p)
                    rr = polynomial_rank([cols[v] for v in core], p)
                    assert rr == len(core)
                    assert polynomial_rank([cols[v] for v in core | moved], p) == rr
                    counts['independent_excess_column_expansion_spans'] += 1
                rows.append({'pole_order': q, 'characteristic': p, 'chart': m,
                             'axes': (A, B), 'heights': heights, 'complete_core_rank': r0})
    return rows


def rational_specs():
    return [
        (0, (-1, 1), (1, 1, 1), (-1, 1), 'numerator_zero_at_quotient_support'),
        (2, (1, 1), (1, 1, 1), (1, 1), 'irreducible_poles_nilpotent_numerator'),
        (2, (1, 1), polynomial_power((1, 1, 1), 2, 2), (1, 1),
         'repeated_irreducible_poles_nilpotent_numerator'),
        (101, (1, 1), (6, -5, 1), (1, 1), 'distinct_poles_nilpotent_numerator'),
        (2, (1, 1, 1), (1, 1, 0, 1), None, 'cubic_denominator'),
        (3, (-1, 0, 1), (1, 2, 0, 1), (-1, 1), 'cubic_poles_nonunit_numerator'),
        (2, (1, 1), (1, 0, 1, 0, 0, 1), (1, 1), 'quintic_denominator'),
        (101, (2, -1, 1), (1, 1, 1), None, 'equal_numerator_denominator_degrees'),
    ]


def general_span_controls(counts):
    """Clear rational functions before quotienting, with numerator powers separate."""
    def terms(D):
        return {(a-b, c, c) for a, b, c in D}

    def columns(points, P, Qpoly, p):
        shift = max(0, -min(v[0] for v in points))
        C = max(v[2] for v in points)
        out = {}
        for v in points:
            k, r, denominator = v
            left, right = polynomial_power(tuple(P), r, p), polynomial_power(tuple(Qpoly), C-denominator, p)
            row = [0]*(len(left)+len(right)-1)
            for i, x in enumerate(left):
                for j, y in enumerate(right):
                    row[i+j] = scalar(row[i+j]+x*y, p)
            out[v] = [0]*(k+shift)+row
        width = max(map(len, out.values()))
        return {v: row+[0]*(width-len(row)) for v, row in out.items()}

    rows = []
    for p, P, Qpoly, _, label in rational_specs():
        q = len(Qpoly)-1
        for m, axis in product(range(q), (0, 1)):
            A, B = 2*q-m, m+q+1
            heights = {j: max(0, 3-abs(j)) for j in range(-m, q-m)}
            ell = min(heights.values())
            old = terms(allocation(A, B, {j: ell for j in heights}))
            W = terms(allocation(A-q*(axis == 0), B-q*(axis == 1),
                                 {j: ell+1 for j in heights}))
            moved = {(k+axis*q, r, d+1) for k, r, d in old}
            cols = columns(W | moved, P, Qpoly, p)
            rr = polynomial_rank([cols[v] for v in W], p)
            assert rr == len(W) and polynomial_rank([cols[v] for v in W | moved], p) == rr
            counts['independent_general_rational_core_span_equalities'] += 1
            for r in range(ell+1, max(heights.values())+1):
                core = terms(allocation(q-m-1, m, {j: r+1 for j in heights}))
                excess = {(j+axis*q, r, r+1) for j, c in heights.items() if c >= r}
                cols = columns(core | excess, P, Qpoly, p)
                rr = polynomial_rank([cols[v] for v in core], p)
                assert rr == len(core) and polynomial_rank([cols[v] for v in core | excess], p) == rr
                counts['independent_general_rational_excess_span_inclusions'] += 1
            rows.append({'label': label, 'characteristic': p, 'numerator': P,
                         'denominator': Qpoly, 'chart': m, 'axis': axis, 'core_rank': len(W)})
    return rows


def budget(D, w, q, m, f, p, ev, counts, cache, unit=None):
    n = len(D)
    A, B, heights = parse(D, q, m)
    alpha, beta = Q(w[0], w[2]), Q(w[1], w[2])
    weights = (alpha, beta, Q(1))
    T, a = min(alpha, beta), q-m
    tau = {j: (j*alpha if j >= 0 else -j*beta) for j in heights}
    assert all(0 <= v <= q*T for v in tau.values())
    if n < q:
        return 'small_length', None
    cheap = int(beta < alpha)
    C, extra = divmod(n-q, q)
    standard = allocation(a-1+extra*(cheap == 0), m+extra*(cheap == 1),
                          {j: C for j in heights})
    if (m, cheap) not in cache:
        assert len(standard) == n and rank([ev(v) for v in standard], p) == n
        cache[m, cheap] = True
        counts['independent_unconditional_complete_pole_bases'] += 1
    assert max(dot(weights, v) for v in D) <= max(dot(weights, v) for v in standard)
    assert max(A, B) <= 2*q+Q(n, q)/T
    if A < a-1 or B < m:
        if n >= 4*q:
            assert n <= 2*q+q*q*T and max(A, B) <= 4*q
            assert max(profile(D)[:2]) <= 4*q*n
            counts['independent_collective_boundary_box_controls'] += 1
        return 'missing_tag', None
    R, other, required = (A, B, a-1) if cheap == 0 else (B, A, m)
    S, M = sum(heights.values()), n-1-other
    v = {j: Q(c)-Q(S, q) for j, c in heights.items()}
    norm2 = sum(x*x for x in v.values())
    P = sum(tau.values())
    energy = norm2/2+sum(tau[j]*v[j] for j in heights)
    excess = S-q*min(heights.values())
    assert excess*excess <= q*q*norm2

    def g(r):
        return Q(q*T+1, 2*q)*r*r+(T/2-Q(M, q)-Q(1, 2)-P/q)*r

    derivative = Q(q*T+1, q)*R+T/2-Q(M, q)-Q(1, 2)-P/q
    repaired = None
    if R >= required+q:
        if cheap == 0 and unit is None:
            shifted = allocation(A-q, B, {j: c+1 for j, c in heights.items()})
            assert rank([ev(e) for e in shifted], p) == n
            assert dot(weights, profile(shifted))-dot(weights, profile(D)) == g(R-q)-g(R) >= 0
            assert derivative <= (q*T+1)/2
            counts['independent_exact_forward_transport_costs'] += 1
        else:
            E, repaired = completion(D, weights, q, m, ev, p, counts, cheap, unit)
            assert dot(weights, profile(D)) <= dot(weights, profile(E))
            literal = allocation(A-q*(cheap == 0), B-q*(cheap == 1),
                                 {j: c+1 for j, c in heights.items()})
            assert dot(weights, profile(literal))-dot(weights, profile(D)) == g(R-q)-g(R)
            assert g(R-q)-g(R)+q*T*excess >= 0
            assert derivative <= (q*T+1)/2+T*excess
    else:
        assert R < 2*q and derivative <= 2*q*T
    C, extra = divmod(S, q)
    balanced = allocation(A+extra*(cheap == 0), B+extra*(cheap == 1),
                          {j: C for j in heights})
    assert rank([ev(e) for e in balanced], p) == n
    increase = g(R+extra)-g(R)
    assert increase <= 4*q*q*T+q+q*T*excess
    assert dot(weights, profile(balanced))-dot(weights, profile(D)) == increase-energy >= 0
    assert energy <= 4*q*q*T+q+q*T*excess
    U = max(Q(1), T)
    assert norm2 <= 64*q**4*U*U
    X = sum(j*c for j, c in heights.items() if j > 0)
    Y = sum(-j*c for j, c in heights.items() if j < 0)
    assert (X-Q(S, q)*sum(j for j in heights if j > 0))**2 <= 64*q**7*U*U
    assert (Y-Q(S, q)*sum(-j for j in heights if j < 0))**2 <= 64*q**7*U*U
    counts['independent_ordinary_variance_and_two_moment_budgets'] += 1
    return 'ordinary_variance', repaired


def fixture(f, p, q, label, counts, numerator=None, denominator=None):
    if denominator is None:
        n, ev, matrices = model(f, p, q)
        unit = None
    else:
        n, ev, matrices, unit = rational_model(f, p, numerator, denominator)
        counts['independent_general_rational_quotient_models'] += 1
        counts['nonunit_numerator_quotient_models'] += rank(matrices[2], p) < n
        counts['models_with_F_one_zero'] += not scalar(sum(f), p)
    actual = fan(n, ev, p, counts)
    cols = cleared_columns(tuple(hyperbolic_cross(n, 3)), q, p, denominator, numerator)
    replay = {v: remainder(c, f, p) for v, c in cols.items()}
    cache, keys, scopes, repairs = {}, {}, Counter(), []
    for D, poly in actual.items():
        assert rank([replay[v] for v in D], p) == n
        counts['independent_cleared_numerator_basis_replays'] += 1
        area = Q(0)
        for m in range(q):
            patch = clip(list(poly), (q-m, -m, 0), 1)
            patch = clip(patch, (-(q-m-1), m+1, 0), 1)
            if not area2(patch):
                continue
            area += area2(patch)
            w = tuple(sum(v[i] for v in patch)/len(patch) for i in range(3))
            A, B, heights = parse(D, q, m)
            key = m, A, B, sum(j*c for j, c in heights.items() if j > 0), sum(-j*c for j, c in heights.items() if j < 0)
            assert key not in keys or keys[key] == D
            keys[key] = D
            which, record = budget(D, w, q, m, f, p, ev, counts, cache, unit)
            scopes[which] += 1
            if record and record['literal_rank'] < n and len(repairs) < 3:
                repairs.append(dict(record, weight=list(map(str, w)), axes=(A, B), heights=heights))
            counts['independent_whole_chart_patches_and_tuple_injections'] += 1
        assert area == area2(poly)
        counts['independent_exact_chart_area_covers'] += 1
    digest = hashlib.sha256(json.dumps([(sorted(D), polygon_key(P)) for D, P in
        sorted(actual.items(), key=lambda v: sorted(v[0]))], default=str).encode()).hexdigest()
    return {'pole_order': q, 'characteristic': p, 'length': n, 'label': label,
            'modulus': f, 'numerator': numerator or (1,),
            'denominator': denominator or polynomial_power((-1, 1), q, p),
            'numerator_unit_mod_F': rank(matrices[2], p) == n,
            'initial_ideals': len(actual), 'budget_scopes': scopes,
            'dependent_shift_repair_witnesses': repairs, 'complete_fan_sha256': digest}


def larger_controls(counts):
    rows = []
    for q, p in ((2, 101), (3, 2), (4, 3), (5, 3), (6, 5), (7, 2), (8, 3)):
        n = max(32, 4*q+4)
        f = [1, -1]+[0]*(n-2)+[1]
        _, ev, _ = model(f, p, q)
        ground, cache = tuple(hyperbolic_cross(n, 3)), {}
        for m in range(q):
            ratio, other = Q(2*m+1, 2), q-Q(2*m+1, 2)
            for T in (Q(1, 16), Q(1), Q(n, 4), Q(2*n)):
                w = (ratio*T/min(ratio, other), other*T/min(ratio, other), Q(1))
                D, _ = field_basis(sorted(ground, key=lambda e: (dot(w, e), e)), ev, p)
                which, record = budget(frozenset(D), w, q, m, f, p, ev, counts, cache)
                rows.append({'length': n, 'pole_order': q, 'characteristic': p,
                             'chart': m, 'weight': list(map(str, w)), 'scope': which,
                             'completion': record})
                counts['independent_larger_ordinary_order_greedy_scans'] += 1
    return rows


def general_larger_controls(counts):
    rows = []
    for p, P, Qpoly, root, label in rational_specs():
        if not p:
            continue
        q, n = len(Qpoly)-1, 24
        f = (polynomial_power(root, n, p) if root else
             tuple(Qpoly)+(0,)*(n-q-1)+(1,))
        _, ev, _, unit = rational_model(f, p, P, Qpoly)
        ground, cache = tuple(hyperbolic_cross(n, 3)), {}
        for m in range(q):
            ratio, other = Q(2*m+1, 2), q-Q(2*m+1, 2)
            for T in (Q(1, 16), Q(1), Q(n, 4), Q(2*n)):
                w = (ratio*T/min(ratio, other), other*T/min(ratio, other), Q(1))
                D, _ = field_basis(sorted(ground, key=lambda e: (dot(w, e), e)), ev, p)
                scope, record = budget(frozenset(D), w, q, m, f, p, ev, counts, cache, unit)
                rows.append({'label': label, 'length': n, 'characteristic': p,
                             'numerator': P, 'denominator': Qpoly, 'modulus': f,
                             'chart': m, 'weight': list(map(str, w)), 'scope': scope,
                             'completion': record})
                counts['independent_general_rational_larger_scans'] += 1
    return rows


def lower_controls(counts):
    rows, arithmetic = [], []
    specs = [(q, p, None, None, None) for q, p in ((2, 101), (3, 2), (4, 3), (5, 3))]
    specs += [(len(denominator)-1, p or 101, P, denominator, root)
              for p, P, denominator, root, _ in rational_specs()[:4]]
    for q, p, P, denominator, root in specs:
        n = max(64, 16*q)
        if denominator is None:
            f = [1, -1]+[0]*(n-2)+[1]
            _, ev, _ = model(f, p, q)
        else:
            f = polynomial_power(root, n, p)
            _, ev, _, _ = rational_model(f, p, P, denominator)
        witnesses = []
        for A in range((2*n+4)//5, 9*n//20+1):
            for B in range((2*n+4)//5, 9*n//20+1):
                if (n-1-A-B) % q:
                    continue
                C = (n-1-A-B)//q
                level = Q(C)+Q(3, 4)
                alpha, beta = level/Q(2*A+1, 2), level/Q(2*B+1, 2)
                if (q*alpha/(alpha+beta)).denominator == 1:
                    alpha += min(alpha, beta)/(1000*q*n)
                m = int(q*alpha/(alpha+beta))
                D = allocation(A, B, {j: C for j in range(-m, q-m)})
                assert len(D) == n and q*min(alpha, beta) < Q(3, 4)
                assert A*alpha < level < (A+1)*alpha and B*beta < level < (B+1)*beta
                assert all(C+(j*alpha if j >= 0 else -j*beta) < level < C+1
                           for j in range(-m, q-m))
                assert rank([ev(v) for v in D], p) == n
                witnesses.append({'axes': (A, B), 'ray_length': C, 'chart': m,
                                  'weight': list(map(str, (alpha, beta, 1)))})
                counts['independent_forced_complete_pole_lower_witnesses'] += 1
        assert witnesses
        rows.append({'pole_order': q, 'characteristic': p, 'length': n,
                     'numerator': P, 'denominator': denominator, 'modulus': f,
                     'witnesses': witnesses})
    for q, factor in product(range(1, 17), (100, 200, 1000)):
        n = factor*q
        low, high = (2*n+4)//5, 9*n//20
        number = sum((high-(n-1-A-low) % q-low)//q+1
                     for A in range(low, high+1))
        brute = sum((n-1-A-B) % q == 0 for A in range(low, high+1)
                    for B in range(low, high+1)) if n <= 400 else None
        if brute is not None:
            assert number == brute
        assert number*q*4000 >= n*n
        arithmetic.append({'pole_order': q, 'length': n, 'forced_axis_pairs': number})
    return rows, arithmetic


def main():
    counts, cases = defaultdict(int), []
    spans = rational_span_controls(counts)
    general_spans = general_span_controls(counts)
    for q, p in ((1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (3, 2),
                 (4, 3), (5, 3), (6, 2), (7, 3), (8, 3), (2, 2)):
        for n in (4, 7):
            cases.append(([1, -1]+[0]*(n-2)+[1], p, q, 'unit_modulus'))
        n = 6
        f = [comb(n, i)*(-1)**(n-i) for i in range(n+1)]
        f[1] += 1
        cases.append((f, p, q, 'prescribed_Taylor_relation'))
    cases += [((-1, 6, -10, 10, -5, 1), 0, 2, 'literal_rank_failure_modulus'),
              ([1, -1]+[0]*14+[1], 2, 3, 'larger_non_Frobenius_complete_fan'),
              ([1, -1]+[0]*18+[1], 3, 4, 'larger_non_Frobenius_complete_fan'),
              (polynomial_power((1, 1, 1), 4, 2), 2, 5, 'repeated_irreducible')]
    fans = [fixture(f, p, q, label, counts) for f, p, q, label in cases]
    for p, P, denominator, root, label in rational_specs():
        q = len(denominator)-1
        for n in sorted({max(4, q+1), 7}):
            f = (polynomial_power(root, n, p) if root else
                 tuple(denominator)+(0,)*(n-q-1)+(1,))
            fans.append(fixture(f, p, q, label, counts, P, denominator))
        fans.append(fixture(P, p, q, label+'_zero_z', counts, P, denominator))
    f = (-1, 6, -10, 10, -5, 1)
    _, ev, _ = model(f, 0, 2)
    D = allocation(0, 3, {-1: 0, 0: 1})
    E, counter = completion(D, (Q(2), Q(1), Q(1)), 2, 1, ev, 0, counts)
    assert counter['literal_rank'] == 4 and not counter['completion_is_downset']
    counter.update({'modulus': f, 'weight': (2, 1, 1), 'basis_before': sorted(D),
                    'repaired_basis': sorted(E)})
    p, P, denominator, root, _ = rational_specs()[1]
    f = polynomial_power(root, 7, p)
    _, ev, _, unit = rational_model(f, p, P, denominator)
    D, weight = allocation(4, 1, {0: 1, 1: 0}), (3, 10, 13)
    B, _ = field_basis(sorted(hyperbolic_cross(7, 3), key=lambda e: (dot(weight, e), e)), ev, p)
    assert frozenset(B) == D
    E, forward_counter = completion(D, weight, 2, 0, ev, p, counts, 0, unit)
    assert forward_counter['literal_rank'] == 6 and not forward_counter['completion_is_downset']
    forward_counter.update({'characteristic': p, 'modulus': f, 'numerator': P,
                            'denominator': denominator, 'weight': weight,
                            'basis_before': sorted(D), 'repaired_basis': sorted(E)})
    larger = larger_controls(counts)
    general_larger = general_larger_controls(counts)
    lower, arithmetic = lower_controls(counts)
    sources = ['docs/ordinary_pole_completion.md', 'scripts/run_ordinary_pole_checks.py',
               'docs/frobenius_pole_transport.md', 'docs/pole_order_fans.md',
               'scripts/run_frobenius_pole_checks.py', 'scripts/run_pole_order_checks.py',
               'scripts/run_projection_tag_checks.py', 'scripts/run_quartic_hurwitz_checks.py',
               'scripts/run_quadric_cross_checks.py', 'scripts/run_signed_face_checks.py',
               'scripts/run_structural_complexity_checks.py']
    report = {'command': 'python3 -B scripts/run_ordinary_pole_checks.py',
              'uniform_scope': 'Every field; monic Q of degree q>=1 with Q(0)!=0; deg P<=q, gcd(P,Q)=1; monic F avoiding0 and Q. Repeated/non-rational poles, nonreduced F and nonunit/zero P modulo F are allowed.',
              'fan_upper': 'C*q^10*n^2*(1+log(n+1)), absolute C',
              'fan_lower': 'c*n^2/q for n>=100q, absolute positive c',
              'counts': dict(sorted(counts.items())), 'rational_completion_span_controls': spans,
              'general_rational_completion_span_controls': general_spans,
              'complete_fans': fans, 'literal_rank_failure_repaired': counter,
              'general_forward_rank_failure_repaired': forward_counter,
              'larger_ordinary_order_scans': larger, 'common_lower_basis_witnesses': lower,
              'general_rational_larger_scans': general_larger,
              'uniform_lower_count_controls': arithmetic,
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              'exact_inverse_transport_asserted': False, 'pole_order_power_sharp': False,
              'general_exact_forward_transport_asserted': False,
              'numerator_coordinate_invertibility_assumed': False,
              'completion_required_to_be_downset': False,
              'logarithmic_gap_closed': False, 'unrestricted_G3_improved': False,
              'polynomial_lower_separation_proved': False, 'novelty_established': False,
              'goal_complete': False}
    (ROOT/'results/research/ordinary_pole_checks.json').write_text(json.dumps(report, indent=2, default=str)+'\n')
    print(json.dumps({'counts': report['counts'], 'complete_fans': len(fans),
                      'whole_fan_cones': sum(v['initial_ideals'] for v in fans),
                      'repaired_counterexample': counter}, indent=2))


if __name__ == '__main__':
    main()
