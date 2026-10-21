#!/usr/bin/env python3
"""Independent integral and polynomial controls for coherent seed transport."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import product
from math import gcd
from pathlib import Path
import hashlib
import json

from run_higher_seed_checks import (
    ceil, column_rank, core_ground, core_model, form_model, summed_rank,
)
from run_higher_interaction_checks import ground, monic_model, rank_trace
from run_signed_face_checks import normals, signed_monic
from run_structural_complexity_checks import dot
from run_three_point_rank_checks import cokernel_signature, source_columns

ROOT = Path(__file__).resolve().parents[1]


def mul(v, w, L):
    out = defaultdict(int)
    for i, a in v.items():
        for j, b in w.items():
            if i+j < L:
                out[i+j] += a*b
    return {i: a for i, a in out.items() if a}


def q_power(e, L):
    v = {0: 1}
    step = {0: 1, 1: 1} if e >= 0 else {j: (-1)**j for j in range(L)}
    for _ in range(abs(e)):
        v = mul(v, step, L)
    return v


def fractional_linear(v, L):
    step = {j+1: -c for j, c in q_power(-1, L).items() if j+1 < L}
    power, out = {0: 1}, defaultdict(int)
    for i in range(L):
        for j, c in power.items():
            out[j] += v.get(i, 0)*c
        power = mul(power, step, L)
    return {j: c for j, c in out.items() if c}


def R_columns(L, a, r, u, k, tau):
    A, B, d = a*(1-u), r-a*(1-u), r-a
    return source_columns(L, [(j, max(0, ceil((B*j-tau)/A)), (k+d*j)//a)
                              for j in range(L)])


def axes(L, h):
    return [q_power(-e, L) for e in range(1, h+1)]


def flag_columns(L, a, r, u, rho, theta, h):
    return R_columns(L, a, r, u, rho, theta)+axes(L, h)


def same_integer_span(L, first, second):
    # Equal free ranks and torsion orders for both lattices and their sum
    # force both finite-index inclusions to have index one.
    signature = cokernel_signature(L, first)
    assert signature == cokernel_signature(L, second)
    assert signature == cokernel_signature(L, first+second)
    return signature


def span_controls():
    counts = defaultdict(int)
    for a, d, L, u in product((1, 2, 3), (1, 2, 4), (1, 3, 6),
                              (F(1, 7), F(3, 7), F(6, 7))):
        r, A = a+d, a*(1-u)
        B = r-A
        for rho in range(a):
            for theta in (-A, F(0), A, B, B+A, B+3):
                H = max(0, theta//A)
                for e in (0, 1, 3):
                    cols = R_columns(L, a, r, u, rho+a*e, theta-A*e)
                    transformed = [mul(v, q_power(-e, L), L) for v in cols]
                    target = flag_columns(L, a, r, u, rho, theta, min(e, H))
                    same_integer_span(L, transformed, target)
                    counts['independent_integral_shift_span_equalities'] += 1
            for k, tau in product((rho, rho+a, rho+3*a), (-A, F(0), B, F(2*r))):
                levels = [(j, max(0, ceil((r*j-tau)/A)), (k+r*j)//a)
                          for j in range(L)]
                canonical = source_columns(L, levels, True)
                J = max(0, min(L-1, tau//r))
                transformed = [mul(fractional_linear(v, L), q_power(J, L), L)
                               for v in canonical]
                target = R_columns(L, a, r, u, k+a*J, tau-A*J)
                same_integer_span(L, transformed, target)
                counts['independent_canonical_regular_span_equalities'] += 1
            D, child_rho = divmod(rho+d, a)
            for theta, h in product((B, B+F(1, 3), B+3*A), range(L+2)):
                parent = flag_columns(L, a, r, u, rho, theta, h)
                if L == 1:
                    assert cokernel_signature(L, parent) == {
                        'free_dimension': 0, 'nonunit_invariants': []}
                else:
                    # Remove the selected unit constant, divide rows by t,
                    # then apply the claimed child row automorphism.
                    reduced = [{j-1: c for j, c in v.items() if j} for v in parent]
                    reduced = [mul(v, q_power(-D, L-1), L-1) for v in reduced]
                    child = flag_columns(L-1, a, r, u, child_rho,
                                         theta-B+A*D, D+h)
                    same_integer_span(L-1, reduced, child)
                    assert cokernel_signature(L, parent) == cokernel_signature(L-1, child)
                counts['independent_integral_transport_equalities'] += 1
    return dict(counts)


def state(s, a, r, u, T, m, h):
    rho, E = m % a, m//a
    L, A = (s-m)//r+1, a*(1-u)
    theta = T-u*m-(1-u)*rho
    k = rho+a*h
    return {'m': m, 'h': h, 'rho': rho, 'E': E, 'L': L, 'A': A,
            'B': r-A, 'theta': theta, 'k': k, 'b': E-h,
            'slack': theta-A*h}


def pairing_controls():
    counts, fixtures = defaultdict(int), []
    for s, a, d, u in product((0, 1, 5, 11, 23, 47), (1, 2, 3, 4),
                              (1, 2, 5), (F(1, 7), F(3, 7), F(6, 7))):
        r = a+d
        for T in (F(0), F(s, 3), F(2*s, 3), F(s+1)):
            positive, negative = [], set()
            for m in range(s+1):
                L, rho, E = (s-m)//r+1, m % a, m//a
                M = min(E+1, (r-1-rho)//a+1)
                for h in range(1, L+1):
                    z = state(s, a, r, u, T, m, h)
                    minus = z['b'] >= 0 and z['k'] >= r and z['slack'] >= 0
                    if minus:
                        negative.add((m, h))
                    plus = []
                    for e in range(1, L):
                        if z['b']+e >= 0 and 0 <= z['k']-a*e < r and z['B']*e <= z['slack']:
                            plus.append(e)
                            positive.append((m, h, e))
                    direct = 0
                    if z['A']*h <= z['theta']:
                        direct = sum(ell+max(0, min(L-1, (z['theta']-z['A']*ell)//r)) >= h
                                     for ell in range(M))
                        direct -= max(0, M-1-h)+int(E >= h)
                    assert direct == len(plus)-int(minus)
                    counts['independent_coefficient_identities'] += 1
            visited, killed = set(), 0
            for m, h, e in positive:
                z = state(s, a, r, u, T, m, h)
                target_m, target_k = m+r*e, z['k']+d*e
                target_h = target_k//a
                target = state(s, a, r, u, T, target_m, target_h)
                assert target_m <= s and target_k//r == e
                # Verify every transport hypothesis, including the last
                # step into a zero flag when the dimension runs out.
                at_m, at_h = m, h
                for _ in range(e):
                    at = state(s, a, r, u, T, at_m, at_h)
                    assert at['theta'] >= at['B']
                    at_h += (at['rho']+d)//a
                    at_m += r
                assert (at_m, at_h) == (target_m, target_h)
                if target_h > target['L']:
                    killed += 1
                    counts['zero_boundary_positive_terms'] += 1
                else:
                    key = target_m, target_h
                    assert key in negative and key not in visited
                    visited.add(key)
                    # Replay the unique inverse independently.
                    back_e = target['k']//r
                    back_m = target_m-r*back_e
                    back_k = target['k']-d*back_e
                    assert (back_m, back_k//a, back_e) == (m, h, e)
                    counts['bijectively_matched_defect_terms'] += 1
            assert visited == negative
            counts['complete_pairing_instances'] += 1
            if s == 23 and a == 2 and d == 1 and u == F(3, 7):
                fixtures.append({'cutoff': s, 'degree': a, 'interaction': r,
                                 'binary_price': str(u), 'price_cutoff': str(T),
                                 'positive_terms': len(positive), 'negative_terms': len(negative),
                                 'zero_boundary_terms': killed})
    return dict(counts), fixtures


def whole_core_controls():
    counts, fixtures = defaultdict(int), []
    # A same-degree pairing would be false even over the integers.
    L, a, r, u, T = 6, 2, 3, F(9, 17), F(17, 3)
    A = a*(1-u)
    regular = R_columns(L, a, r, u, 0, T)
    canonical = source_columns(L, [(j, max(0, ceil((r*j-T)/A)), r*j//a)
                                   for j in range(L)], True)
    assert same_integer_span(L, regular, [{j: 1} for j in range(5)]) == {
        'free_dimension': 1, 'nonunit_invariants': []}
    assert cokernel_signature(L, canonical) == {
        'free_dimension': 0, 'nonunit_invariants': []}
    counts['all_field_same_degree_pairing_obstructions'] += 1
    for s, a, r in ((5, 1, 2), (9, 2, 3), (11, 2, 3), (14, 2, 4),
                     (17, 2, 3), (13, 3, 5), (15, 3, 7), (12, 4, 5)):
        E = core_ground(s, a)
        models = [core_model(s, a, r, c) for c in (False, True)]
        assert len(models[0][0]) == len(models[1][0])
        for u, T in product((F(1, 7), F(3, 7), F(6, 7)),
                             (F(2), F(s, 3), F(2*s, 3), F(s))):
            selected = [v for v in E if a*u*(v[0]+v[1])+v[2] <= T]
            signatures = [cokernel_signature(len(rows), [ev(v) for v in selected])
                          for rows, ev in models]
            assert signatures[0] == signatures[1]
            counts['independent_whole_core_smith_equalities'] += 1
            for p in (0, 2, 3, 5, 7):
                ranks = [column_rank([ev(v) for v in selected], p) for rows, ev in models]
                assert ranks[0] == ranks[1]
                counts['independent_whole_core_field_equalities'] += 1
            if signatures[0]['nonunit_invariants']:
                fixtures.append({'cutoff': s, 'degree': a, 'interaction': r,
                                 'binary_price': str(u), 'price_cutoff': str(T),
                                 'cokernel': signatures[0]})
    return dict(counts), fixtures


def original_controls():
    counts = defaultdict(int)
    for a, r, s in ((1, 3, 7), (2, 3, 9), (2, 5, 11), (3, 4, 10)):
        pairs = [({(a, 0): 1, (0, a): 1}, {(a-1, 1): 1}),
                 ({(a, 0): 1}, {(a-1, 1): 1}),
                 ({(0, a): 1}, {(0, a): 1})]
        E = ground(s)
        for Q, H in pairs:
            rows, ev = form_model(s, a, r, Q, H)
            for u, T, p in product((F(1, 7), F(6, 7)),
                                    (F(2), F(s, 2), F(s)), (0, 2, 3, 5)):
                selected = [v for v in E if u*(v[0]+v[1])+v[2] <= T]
                actual = column_rank([ev(v) for v in selected], p)
                assert actual == summed_rank(s, a, r, u, F(1), T, p)
                assert actual == summed_rank(s, a, r, u, F(1), T, p, True)
                counts['independent_noncoprime_original_seed_queries'] += 1
    # Direct original/face monic recurrences, with no core rank formula.
    for s, r, p, u in product((5, 9, 13), (3, 5), (0, 2, 3),
                              (F(1009, 3001), F(2017, 3001))):
        t = F(1, 1000037)
        seed, v = (u, u, F(1)), (1, 1-r, 1)
        w = tuple(seed[i]+t*v[i] for i in range(3))
        assert w[2] > w[0] > w[1] > 0
        E = ground(s)
        costs = [dot(a, w) for a in E]
        assert len(set(costs)) == len(E)
        ordered = sorted(E, key=lambda a: dot(a, w))
        _, original = monic_model(s, r)
        _, face = signed_monic(s, r)
        B, _ = rank_trace(ordered, original, p, set())
        C, _ = rank_trace(ordered, face, p, set())
        assert B == C
        for eta in normals(B, original, p)[0]:
            assert dot(eta, v) >= 0 and dot(eta, seed) >= 0 and dot(eta, w) > 0
            counts['independent_original_border_inequalities'] += 1
        counts['independent_largest_z_face_basis_equalities'] += 1
    return dict(counts)


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    span = span_controls()
    pairing, pairing_fixtures = pairing_controls()
    core, smith_fixtures = whole_core_controls()
    original = original_controls()
    sources = ('docs/higher_seed_transport.md', 'scripts/run_higher_seed_transport_checks.py',
               'docs/higher_seed_binomial_reduction.md', 'docs/three_point_rank_ideals.md',
               'docs/signed_face_reduction.md', 'docs/middle_interaction_fans.md',
               'docs/equal_degree_unit_fans.md',
               'docs/nonnegative_graded_fans.md', 'docs/higher_interaction_fans.md',
               'scripts/run_higher_seed_checks.py', 'scripts/run_three_point_rank_checks.py',
               'scripts/run_higher_interaction_checks.py', 'scripts/run_signed_face_checks.py',
               'scripts/run_structural_complexity_checks.py',
               'scripts/run_hypersurface_truncation_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_socle_stability_checks.py')
    report = {'status': 'passed',
              'command': 'python3 -B scripts/run_higher_seed_transport_checks.py',
              'counts': {**span, **pairing, **core, **original},
              'pairing_fixtures': pairing_fixtures,
              'nonunit_whole_core_fixtures': smith_fixtures,
              'uniform_results': {
                  'regular_canonical_integral_cokernels_isomorphic': True,
                  'all_field_core_rank_identity_B8': True,
                  'all_nonzero_equal_degree_form_seed_rank_law': True,
                  'original_full_signed_transfer_r_ge_2': True,
                  'original_full_family_fan_bound': '25*n+14, n=(s+1)^2, uniformly in r>=2 and field',
                  'arbitrary_positive_prices_or_downward_staircases': False,
                  'integral_universality_for_arbitrary_integer_forms': False},
              'unrestricted_G3_improvement': False, 'novelty_established': False,
              'goal_complete': False,
              'scope': 'Uniform proofs close the whole higher-interaction family as an amplification route. Independent finite checks support the integral transport, exact global cancellation and original polynomial consequences. No unrestricted endpoint, new prime-distribution theorem or priority certificate.',
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/'results/research/higher_seed_transport_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
