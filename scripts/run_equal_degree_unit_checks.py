#!/usr/bin/env python3
"""Exact controls for equal-degree unit ranks and the full r=2 family."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb
from pathlib import Path
import hashlib
import json

from run_higher_interaction_checks import (
    border_certificate, ground, linear_sum, monic_model, rank_trace,
    sector_prefix_rank, wall_formula,
)
from run_hypersurface_truncation_checks import relation_rows
from run_middle_interaction_checks import blocks, cap_floor, middle_rank, raw_laurent
from run_signed_face_checks import normals, signed_monic
from run_structural_complexity_checks import area2, clip, determinant, dot

ROOT = Path(__file__).resolve().parents[1]


def unit_rank(s, d, u, Z, T):
    """Complete original tier rank, Q,H nonzero of degree d, 0<u<=Z."""
    u, Z, T = map(F, (u, Z, T))
    if not 0 < u <= Z:
        raise ValueError('This theorem requires 0 < u <= Z.')
    out = 0
    for e in range(min(d, s+1)):
        for A in range(s-e+1):
            D = (s-A-e)//d
            tau = T-u*A-Z*e
            if tau >= 0:
                amount = cap_floor(tau, d*u, D+1)
                amount += cap_floor(tau, d*Z, D, 1)
                out += min(d, A+1)*min(D+1, amount)
    return out


def generic_inventory(s, d, u, Z):
    out = Counter()
    channels = Counter(a+b for a in range(d) for b in range(d))
    for delta, multiplicity in channels.items():
        for e in range(min(d, s-delta+1)):
            for q in range((s-delta-e)//d+1):
                B = (s-delta-e)//d-q
                base = u*(delta+d*q)+Z*e
                costs = [base+d*u*b for b in range(B+1)]
                costs += [base+d*(Z-u)*i for i in range(1, q+1)]
                for cost in sorted(costs)[:B+1]:
                    out[cost] += multiplicity
    return out


def canonical_inventory(s, d, u, Z):
    out = Counter()
    for e in range(min(d, s+1)):
        for A in range(s-e+1):
            D, base = (s-A-e)//d, u*A+Z*e
            costs = [base+d*u*a for a in range(D+1)]
            costs += [base+d*Z*i for i in range(1, D+1)]
            for cost in sorted(costs)[:D+1]:
                out[cost] += min(d, A+1)
    return out


def unit_model(s, d, Q, H):
    """Direct original polynomial recurrence; no axis/rank formula."""
    basis = tuple(a for a in ground(s) if a[2] < d)
    index = {a: i for i, a in enumerate(basis)}

    @lru_cache(None)
    def evaluate(a):
        x, y, z = a
        if sum(a) > s:
            return {}
        if z < d:
            return {index[a]: 1}
        tails = [(-c, evaluate((x+i, y+j, z-d))) for (i, j), c in Q.items()]
        tails += [(-c, evaluate((x+i, y+j, z))) for (i, j), c in H.items()]
        return linear_sum(tails)
    return basis, evaluate


def inventory_controls():
    counts = defaultdict(int)
    for s, d, u in product((2, 5, 9, 17, 31, 60), range(1, 9),
                            (F(2, 17), F(9, 17), F(16, 17))):
        generic = generic_inventory(s, d, u, F(1))
        canonical = canonical_inventory(s, d, u, F(1))
        assert generic == canonical
        expected = sum(comb(s-e+2, 2) for e in range(min(d, s+1)))
        assert sum(generic.values()) == expected
        for T in (F(0), F(s, 7), F(3*s, 4), F(s)):
            assert unit_rank(s, d, u, F(1), T) == \
                   sum(v for cost, v in canonical.items() if cost <= T)
            counts['independent_floor_inventory_queries'] += 1
        counts['exact_generic_canonical_inventories'] += 1
        counts['matched_price_birth_items'] += expected
    # Literal birth superlevels, separate from comparing finite inventories.
    for d, theta in product(range(1, 6), (F(337, 1009), F(1703, 1009))):
        channels = Counter(a+b for a in range(d) for b in range(d))
        for A in range(41):
            heights = []
            for delta, mult in channels.items():
                if A < delta or (A-delta) % d:
                    continue
                k = (A-delta)//d
                for b in range(k+1):
                    h = min(k-b, (b/theta).numerator//(b/theta).denominator)
                    heights += [h]*mult
            other = []
            for j in range(min(d, A+1)):
                for a in range((A-j)//d+1):
                    # c=A-j-da is then nonnegative and uniquely specified.
                    value = F(a)/(1+theta)
                    other.append(value.numerator//value.denominator)
            assert Counter(heights) == Counter(other)
            counts['literal_birth_height_multisets'] += 1
    return dict(counts)


def form_controls():
    counts, fixtures = defaultdict(int), []
    pairs = []
    for d in range(1, 5):
        pairs.extend((
            (d, {(d, 0): 1, (0, d): 1}, {(d-1, 1): 1}),
            (d, {(d, 0): 1}, {(d-1, 1): 1}),
            (d, {(0, d): 1}, {(0, d): 1}),
            (d, {(d, 0): 1}, {(0, d): 1}),
        ))
    for label, (d, Q, H) in enumerate(pairs):
        s, E = d+6, ground(d+6)
        basis, evaluate = unit_model(s, d, Q, H)
        assert len(basis) == sum(comb(s-e+2, 2) for e in range(min(d, s+1)))
        for u, p in product((F(2, 7), F(5, 7), F(17, 19), F(1)), (0, 2, 3, 5)):
            w = u, u, F(1)
            selected, _ = rank_trace(sorted(E, key=lambda a: (dot(a, w), a)),
                                     evaluate, p, set())
            actual = Counter(dot(a, w) for a in selected)
            assert actual == canonical_inventory(s, d, u, F(1))
            counts['original_polynomial_field_models'] += 1
            for T in (F(0), F(2), F(7, 2), F(s, 2), F(s)):
                assert unit_rank(s, d, u, F(1), T) == \
                       sum(v for cost, v in actual.items() if cost <= T)
                counts['original_polynomial_floor_queries'] += 1
        if label % 4 in (0, 2):
            cutoff, E0 = 5, ground(5)
            B0, ev = unit_model(cutoff, d, Q, H)
            f = {(0, 0, d): 1}
            f.update({(i, j, 0): c for (i, j), c in Q.items()})
            f.update({(i, j, d): c for (i, j), c in H.items()})
            for p in (0, 2, 3):
                rows = relation_rows(E0, [f], p)
                assert len(E0)-len(rows) == len(B0)
                for row in rows:
                    image = linear_sum((c, ev(a)) for a, c in zip(E0, row) if c)
                    assert all(not (c % p if p else c) for c in image.values())
                    counts['independent_macaulay_kernel_rows'] += 1
                counts['independent_full_kernel_equalities'] += 1
        fixtures.append({'pair': label, 'degree': d, 'cutoff': s,
                         'Q': [[*a, c] for a, c in sorted(Q.items())],
                         'H': [[*a, c] for a, c in sorted(H.items())]})
    # Nonzero hypotheses are tested separately rather than silently generalized.
    failures = []
    for Q, H, label in (({}, {(1, 1): 1}, 'Q=0'),
                         ({(2, 0): 1}, {}, 'H=0')):
        s, d, u = 8, 2, F(5, 7)
        _, ev = unit_model(s, d, Q, H)
        w = u, u, F(1)
        B, _ = rank_trace(sorted(ground(s), key=lambda a: (dot(a, w), a)), ev, 0, set())
        costs = sorted({dot(a, w) for a in ground(s)})
        bad = next(T for T in costs if sum(dot(a, w) <= T for a in B) != unit_rank(s, d, u, 1, T))
        failures.append({'case': label, 'cutoff': s, 'degree': d,
                         'weight': list(map(str, w)), 'price_cutoff': str(bad),
                         'actual_rank': sum(dot(a, w) <= bad for a in B),
                         'inapplicable_formula_rank': unit_rank(s, d, u, 1, bad)})
    return dict(counts), fixtures, failures


def negative_block_rank(m, L, alpha, beta, T):
    """Independent parity caps in the signed r=2 block, -beta<alpha<0."""
    alpha, beta, T = map(F, (alpha, beta, T))
    assert -beta < alpha < 0
    A, C = -alpha, alpha+beta
    out = 0
    for e in (0, 1):
        j0 = max(0, e-m)
        Le, M = L-j0, m+j0-e
        if Le <= 0:
            continue
        tau = T-alpha*m-C*j0-A*e
        axis = cap_floor(tau, C, Le)
        maximum = (M+Le-1)//2
        first = min(maximum, M//2)
        boundary = cap_floor(tau, 2*A, first, 1)
        boundary += cap_floor(tau+C*M, 2*beta, maximum-first, first+1)
        out += min(Le, axis+boundary)
    return out


def negative_candidates(m, L, alpha, beta):
    chosen = []
    for e in (0, 1):
        j0 = max(0, e-m)
        Le, M = L-j0, m+j0-e
        if Le <= 0:
            continue
        pairs = [(m+j-e, j) for j in range(j0, L)]
        for h in range(1, (M+Le-1)//2+1):
            j = j0+max(0, 2*h-M)
            pairs.append((m+j-2*h-e, j))
        assert len(pairs) == len(set(pairs))
        chosen += sorted(pairs, key=lambda a: (alpha*a[0]+beta*a[1], a))[:Le]
    assert len(chosen) == len(set(chosen)) == 2*L-int(m == 0)
    return chosen


def negative_flag_controls():
    counts = defaultdict(int)
    for m, L in product((0, 1, 2, 5, 8), (1, 2, 4, 6)):
        cols = raw_laurent(m, 1, L)
        for rho in (-F(1, 101), -F(53, 101), -F(97, 101)):
            price = lambda a: rho*a[0]+a[1]
            chosen = negative_candidates(m, L, rho, F(1))
            rows = [i for i in range(2*L) if m or i != 1]
            matrix = [[cols[a].get(i, 0) for a in chosen] for i in rows]
            assert abs(determinant(matrix)) == 1
            counts['negative_flag_integer_unit_determinants'] += 1
            for p in (0, 2, 3, 5):
                B, _ = rank_trace(sorted(cols, key=lambda a: (price(a), a)),
                                  cols.__getitem__, p, set())
                for T in {price(a) for a in cols}:
                    assert negative_block_rank(m, L, rho, 1, T) == sum(price(a) <= T for a in B)
                    counts['negative_flag_independent_prefix_ranks'] += 1
                assert set(chosen) == B
                counts['negative_flag_field_bases'] += 1
    return dict(counts)


def full_quadratic_rank(s, weight, T):
    X, Y = sorted(map(F, weight[:2]), reverse=True)
    Z, T = F(weight[2]), F(T)
    if min(X, Y, Z) <= 0:
        raise ValueError('Positive prices are required.')
    if Z <= Y:
        return sector_prefix_rank(s, 2, (X, Y, Z), T)
    if Z <= X:
        return middle_rank(s, 2, (X, Y, Z), T)
    alpha, beta = X-Z, Y+Z
    return sum(negative_block_rank(m, L, alpha, beta, T-k*Z-b0*beta)
               for k, b0, L, m in blocks(s, 2))


def full_rank_controls():
    counts, fixtures = defaultdict(int), []
    weights = ((5, 2, 11), (2, 5, 11), (11, 2, 5), (2, 11, 5),
               (11, 5, 2), (5, 11, 2), (1, 1, 7), (3, 2, 3), (1, 1, 1))
    for s in (3, 6, 9, 13):
        E = ground(s)
        _, ev = monic_model(s, 2)
        for p in (0, 2, 3, 5):
            for weight in weights:
                w = tuple(map(F, weight))
                B, _ = rank_trace(sorted(E, key=lambda a: (dot(a, w), a)), ev, p, set())
                costs = sorted({dot(a, w) for a in E})
                tests = {F(-1), F(0), max(costs)}
                for i, T in enumerate(costs):
                    if any(dot(a, w) == T for a in B):
                        tests.add(T)
                        if i:
                            tests.add(costs[i-1])
                for T in tests:
                    assert full_quadratic_rank(s, w, T) == sum(dot(a, w) <= T for a in B)
                    counts['full_original_independent_rank_queries'] += 1
                counts['full_original_field_models'] += 1
        fixtures.append({'cutoff': s, 'length': (s+1)**2,
                         'weights': [list(w) for w in weights],
                         'characteristics': [0, 2, 3, 5]})
        print('Full original r=2 rank controls', s, flush=True)
    return dict(counts), fixtures


def entire_patch_controls():
    """Exact all-field original certificates, not sampled cost inequalities."""
    counts, fixtures = defaultdict(int), []
    for s in (4, 8, 12):
        E = ground(s)
        _, face = signed_monic(s, 2)
        _, original = monic_model(s, 2)
        lower, rho, delta = F(-1), -F(1)+F(1, 100*(s+1)**3), F(1, 100*(s+1)**3)
        patches, seen, total_area = [], set(), F(0)
        while True:
            Y = (1+rho)/(2*(1-rho)) if rho < 0 else F(1, 2)
            w = 1+rho*(Y+1), Y, F(1)
            assert min(w[0], w[2]) > w[1] > 0
            B, _ = rank_trace(sorted(E, key=lambda a: (dot(a, w), a)), face, 0, set())
            key = frozenset(B)
            assert key not in seen
            seen.add(key)
            face_normals, _ = normals(sorted(B), face, 0)
            lo, hi = lower, None
            for x, y, z in face_normals:
                assert x-y+z == 0
                if x > 0:
                    lo = max(lo, -F(y, x))
                elif x < 0:
                    hi = min(hi, -F(y, x)) if hi is not None else -F(y, x)
                else:
                    assert y > 0
            assert lo < rho and (hi is None or rho < hi)
            original_normals, relations, sizes = border_certificate(
                B, original, s, 2, w, save=True)
            assert all(dot(eta, (1, -1, 1)) >= 0 for eta in original_normals)
            polygon = [tuple(F(i == j) for i in range(3)) for j in range(3)]
            for eta in tuple(face_normals)+((1, -1, 0), (0, -1, 1)):
                polygon = clip(polygon, eta, 1)
            assert polygon and area2(polygon) > 0
            assert all(dot(eta, v) >= 0 for eta in original_normals for v in polygon)
            total_area += area2(polygon)
            counts['entire_all_field_original_patches'] += 1
            counts['original_integer_border_relations'] += sizes[0]
            counts['original_integer_membership_terms'] += sizes[1]
            counts['entire_patch_rational_inequalities'] += len(original_normals)*len(polygon)
            patches.append({'rho_interval': [str(lo), str(hi) if hi is not None else 'infinity'],
                            'standard_set': sorted(B),
                            'vertices': [list(map(str, v)) for v in polygon],
                            'original_integer_relations': relations})
            if hi is None:
                break
            lower, rho = hi, hi+delta
        assert total_area == F(1, 3)
        assert len(seen) <= 12*(s+1)**2+7
        counts['entire_cheapest_y_regions'] += 1
        fixtures.append({'cutoff': s, 'order': 2, 'length': (s+1)**2,
                         'normalized_region_area': str(total_area),
                         'all_fields': True, 'patches': patches,
                         'coarse_global_bound': 24*(s+1)**2+2*len(wall_formula(s, 2))+18})
        print('Entire all-field original r=2 region',s,len(seen),flush=True)
    return dict(counts), fixtures


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    inventory = inventory_controls()
    forms, form_fixtures, failures = form_controls()
    flags = negative_flag_controls()
    original, original_fixtures = full_rank_controls()
    patches, patch_fixtures = entire_patch_controls()
    sources = ('docs/equal_degree_unit_fans.md', 'scripts/run_equal_degree_unit_checks.py',
               'docs/middle_interaction_fans.md', 'docs/signed_face_reduction.md',
               'docs/higher_interaction_fans.md', 'docs/nonnegative_graded_fans.md',
               'docs/fixed_equation_bad_primes.md', 'docs/cubic_fan_classification.md',
               'scripts/run_higher_interaction_checks.py',
               'scripts/run_middle_interaction_checks.py', 'scripts/run_signed_face_checks.py',
               'scripts/run_hypersurface_truncation_checks.py',
               'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_socle_stability_checks.py')
    report = {
        'status': 'passed', 'command': 'python3 -B scripts/run_equal_degree_unit_checks.py',
        'counts': {**inventory, **forms, **flags, **original, **patches},
        'original_form_fixtures': form_fixtures,
        'nonzero_hypothesis_counterexamples': failures,
        'original_full_rank_fixtures': original_fixtures,
        'entire_original_region_certificates': patch_fixtures,
        'uniform_results': {
            'equal_degree_nonzero_forms': 'Complete tiers for ((1+H)*z^d+Q)+m^(s+1), deg(Q)=deg(H)=d, prices (u,u,Z), 0<u<=Z, are given by (U5), all fields, without coprimality.',
            'equal_degree_rank_query': 'O(d*(s+1)) exact arithmetic/floor operations, constant auxiliary storage; bit complexity separate.',
            'full_r2_signed_transfer': True,
            'full_r2_fan_count': 'Theta(n), uniformly in s and field, n=(s+1)^2',
            'full_r2_generic_standard_sets': 'independent of characteristic, with integral unit reference determinants',
            'full_r2_coherent_prefix_quotients': 'torsion-free over Z',
            'full_r2_universal_denominator': 1,
            'full_r2_integral_universal_basis': 'The union of the integral reduced bases contains a monic Groebner basis for every order and survives arbitrary base change.',
            'full_r2_rank_query': 'O(s+1) exact arithmetic/floor operations and constant auxiliary storage, at every positive weight.',
            'first_bad_prime_interaction_exponent': 'r=3 within the coefficient-one family: prior r=1 and present r=2 have none; the prior fixed r=3 result realizes every prime at some cutoff.',
            'remaining_power_candidate': '3<=r<sqrt(s), in the largest-z region only.'},
        'r_ge_3_full_signed_transfer_proved': False,
        'general_r_full_family_linear_bound_proved': False,
        'unrestricted_G3_improvement': False, 'novelty_established': False,
        'goal_complete': False,
        'scope': 'Uniform proofs of the equal-degree complete-tier law and full r=2 classification, with independent finite controls. No original r>=3 full fan theorem, unrestricted improvement or priority certificate.',
        'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/'results/research/equal_degree_unit_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
