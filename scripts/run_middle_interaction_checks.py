#!/usr/bin/env python3
"""Exact controls for noncoprime positive ranks and original middle-price cones."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import hashlib
import json

from run_higher_interaction_checks import (
    binary_form_model, border_certificate, channel_price_multiset,
    equal_binary_rank, ground, linear_sum, monic_model, rank_trace,
)
from run_hypersurface_truncation_checks import relation_rows
from run_signed_face_checks import laurent_model, normals, signed_monic
from run_structural_complexity_checks import area2, clip, determinant, dot

ROOT = Path(__file__).resolve().parents[1]


def retained(m, d, L):
    return ([(a, 0) for a in range(m+1)]
            +[(0, j) for j in range(1, L)]
            +[(a, j) for j in range(1, L)
              for a in range(max(1, m+d*(j-1)), m+d*j+1)])


def raw_laurent(m, d, L):
    """Separate polynomial-t arithmetic, not the q-adic binomial formula."""
    N = 2*L
    modulus = [comb(L, i//2) if not i % 2 else 0 for i in range(N+1)]

    def reduce(v):
        v = list(v)
        while len(v) > N:
            c, shift = v[-1], len(v)-N-1
            if c:
                for i, x in enumerate(modulus):
                    v[shift+i] -= c*x
            v.pop()
        return v+[0]*(N-len(v))

    def mul(a, b):
        v = [0]*(len(a)+len(b)-1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                v[i+j] += x*y
        return reduce(v)

    one, t, q = reduce([1]), reduce([0, 1]), reduce([1, 0, 1])
    qpowers = [one]
    for _ in range(L):
        qpowers.append(mul(qpowers[-1], q))
    inverse = [-x for x in mul(t, [sum(a[i] for a in qpowers[:L]) for i in range(N)])]
    assert mul(t, inverse) == one
    powers = {0: one}
    for e in range(1, m+d*(L-1)+1):
        powers[e] = mul(powers[e-1], t)
    for e in range(1, L):
        powers[-e] = mul(powers[1-e], inverse)
    columns = {(a, j): mul(powers[a-j], qpowers[j])
               for j in range(L) for a in range(m+d*j+1)}
    return {a: {i: c for i, c in enumerate(v) if c} for a, v in columns.items()}


def flag_controls():
    counts = defaultdict(int)
    for m, d, L in product((0, 1, 3, 6), (1, 2, 3, 5), (1, 2, 4, 6)):
        columns = raw_laurent(m, d, L)
        candidates = retained(m, d, L)
        cap = 2*L-int(m == 0)
        for rho in (F(0), F(7, 101), F(53, 101), F(97, 101), F(1)):
            cost = lambda a: rho*a[0]+a[1]
            ordered = sorted(columns, key=lambda a: (cost(a), a))
            for p in (0, 2, 3, 5):
                B, _ = rank_trace(ordered, columns.__getitem__, p, set())
                for T in {cost(a) for a in columns}:
                    assert sum(cost(a) <= T for a in B) == \
                           min(cap, sum(cost(a) <= T for a in candidates))
                    counts['independent_laurent_prefix_ranks'] += 1
                counts['independent_laurent_field_models'] += 1
            if 0 < rho < 1:
                chosen = sorted(candidates, key=lambda a: (cost(a), a))[:cap]
                rows = [i for i in range(2*L) if not (m == 0 and i == 1)]
                matrix = [[columns[a].get(i, 0) for a in chosen] for i in rows]
                assert abs(determinant(matrix)) == 1
                counts['independent_integer_flag_determinants'] += 1
    return dict(counts)


def canonical_inventory(s, r, lam):
    out = Counter()
    for c in range(s+1):
        for j in range(min(r-1, s-c)+1):
            L = (s-c-j)//r+1
            costs = sorted(lam*a+r*ell for a, ell in retained(j, r, L))
            out.update(c+j+v for v in costs[:2*L-int(j == 0)])
    assert sum(out.values()) == (s+1)**2
    return out


def canonical_birth_controls():
    counts = defaultdict(int)
    fixtures = []
    for s, r, mu in product((4, 7, 13, 21, 37, 60), range(1, 9),
                            (F(2, 17), F(9, 17), F(16, 17))):
        inventory = canonical_inventory(s, r, r*mu)
        assert inventory == channel_price_multiset(s, r, r*mu)
        counts['canonical_generic_price_multisets'] += 1
        counts['canonical_matched_rank_birth_items'] += sum(inventory.values())
    for r, mu in product(range(1, 6), (F(337, 1009), F(673, 1009))):
        lam = r*mu
        for a in range(1, 13):
            for b in range(r+1):
                j, ell = (a+b) % r, (a+b)//r
                price = a+b+lam*a
                limit = int((price-j)/r)+3
                cheaper = sum(j+r*k+lam*x < price
                              for x, k in retained(j, r, limit))
                position = cheaper+1
                Lfirst = (position+int(j == 0)+1)//2
                actual = j+r*(Lfirst-1)
                psi = (a+(mu*a).numerator//(mu*a).denominator)//2
                assert actual == a+b+r*psi
                counts['literal_nonconstant_birth_positions'] += 1
        heights = []
        for b in range(46):
            j, ell = b % r, b//r
            limit = ell+3
            # Any right item below the left item's price has level < ell.
            cheaper = sum(j+r*k+lam*x < b
                          for x, k in retained(j, r, limit))
            position = cheaper+1
            Lfirst = (position+int(j == 0)+1)//2
            actual = j+r*(Lfirst-1)
            assert (actual-b) % r == 0
            heights.append((actual-b)//r)
        for D in range(46):
            for h in range(1, 9):
                x = 2*mu*h/(1-mu)
                Bh = -((-x.numerator)//x.denominator)
                assert sum(v >= h for v in heights[:D+1]) == \
                       max(0, D-2*h-(r+1)*Bh+1)
                counts['literal_integer_birth_superlevels'] += 1
        fixtures.append({'order': r, 'mu': str(mu),
                         'pure_curve_degree_birth_surpluses': heights})
    return dict(counts), fixtures


def form_controls():
    counts, fixtures = defaultdict(int), []
    pairs = [
        (1, {(0, 2): 1}, {(0, 2): 1}),
        (2, {(1, 1): 1}, {(3, 0): 1}),
        (2, {(1, 1): 1}, {(2, 1): 1, (1, 2): 1}),
        (3, {(0, 2): 1}, {(1, 3): 1}),
        (3, {(0, 2): 1}, {(0, 4): 1}),
        (3, {(2, 0): 1, (0, 2): 1}, {(1, 3): 1}),
        (4, {(2, 0): 1, (1, 1): 2, (0, 2): 1},
             {(i, 5-i): comb(5, i) for i in range(6)}),
    ]
    for label, (r, Q, H) in enumerate(pairs):
        s = r+7
        E = ground(s)
        _, evaluate = binary_form_model(s, Q, H)
        for mu, p in product((F(0), F(2, 17), F(9, 17), F(16, 17), F(1), F(r+1, r)),
                             (0, 2, 3, 5, 7)):
            lam = r*mu
            w = 1+lam, F(1), F(1)
            B, _ = rank_trace(sorted(E, key=lambda a: (dot(a, w), a)), evaluate, p, set())
            actual = Counter(dot(a, w) for a in B)
            assert actual == channel_price_multiset(s, r, lam)
            for T in (F(2), F(s, 2), F(s), max(actual)):
                assert sum(v for cost, v in actual.items() if cost <= T) == \
                       equal_binary_rank(s, r, lam, T)
                counts['noncoprime_floor_rank_queries'] += 1
            counts['noncoprime_original_field_models'] += 1
        fixtures.append({'pair': label, 'order': r, 'cutoff': s,
                         'Q': [[*a, c] for a, c in sorted(Q.items())],
                         'H': [[*a, c] for a, c in sorted(H.items())]})
    for r, Q, H in pairs[:3]:
        for p in (0, 2, 3):
            s = 6
            E = ground(s)
            basis, evaluate = binary_form_model(s, Q, H)
            f = {(2, 0, 0): 1}
            f.update({(0, *a): c for a, c in Q.items()})
            f.update({(1, *a): c for a, c in H.items()})
            relations = relation_rows(E, [f], p)
            assert len(E)-len(relations) == len(basis)
            for row in relations:
                image = linear_sum((c, evaluate(a)) for a, c in zip(E, row) if c)
                assert all(not (c % p if p else c) for c in image.values())
                counts['noncoprime_independent_macaulay_rows'] += 1
            counts['noncoprime_full_kernel_equalities'] += 1
    counterexamples = []
    for Q, H, lam, T, expected, formula in (
        ({(0, 2): 1}, {(0, 4): 1}, -F(1, 3), F(4), 41, 42),
        ({}, {(1, 3): 1}, F(1, 7), F(32, 7), 34, 35),
        ({(0, 2): 1}, {}, F(1, 7), F(16, 7), 9, 10),
    ):
        s, r = 10, 3
        _, evaluate = binary_form_model(s, Q, H)
        w = 1+lam, F(1), F(1)
        prefix = sorted((a for a in ground(s) if dot(a, w) <= T),
                        key=lambda a: (dot(a, w), a))
        B, _ = rank_trace(prefix, evaluate, 0, set())
        assert len(B) == expected
        assert equal_binary_rank(s, r, lam, T) == formula
        counterexamples.append({'Q': [[*a, c] for a, c in sorted(Q.items())],
                                'H': [[*a, c] for a, c in sorted(H.items())],
                                'cutoff': s, 'order': r,
                                'weight': list(map(str, w)), 'price': str(T),
                                'parent_prefix_size': len(prefix), 'actual_rank': expected,
                                'out_of_scope_formula': formula,
                                'selected_monomials': sorted(B)})
        counts['scope_counterexamples'] += 1
    return dict(counts), fixtures, counterexamples


def blocks(s, r):
    for k in range(-(r-1)*s, s+1):
        b0 = max(0, -(k//(r-1)))
        L = (s-k)//r-b0+1
        if L > 0:
            yield k, b0, L, k+(r-1)*b0


def middle_basis(s, r, w):
    X, Y, Z = w
    assert X > Z > Y > 0
    if X > Y+r*Z:
        return frozenset(a for a in ground(s) if a[0] < 2)
    alpha, beta = X-Z, Y+(r-1)*Z
    out = []
    for k, b0, L, m in blocks(s, r):
        candidates = [(a, b0+j, m+(r-1)*j-a) for a, j in retained(m, r-1, L)]
        out.extend(sorted(candidates, key=lambda a: (dot(a, w), a))[:2*L-int(m == 0)])
    assert len(out) == len(set(out)) == (s+1)**2
    return frozenset(out)


def cap_floor(T, step, length, start=0):
    if step == 0:
        return length if T >= 0 else 0
    value = T/step
    return max(0, min(length, value.numerator//value.denominator-start+1))


def middle_rank(s, r, w, T):
    X, Y, Z = map(F, w)
    T = F(T)
    assert X >= Z >= Y > 0
    alpha, beta = X-Z, Y+(r-1)*Z
    if alpha >= beta:
        out = 0
        for a in range(min(1, s)+1):
            for b in range(s-a+1):
                out += cap_floor(T-X*a-Y*b, Z, s-a-b+1)
        return out
    out = 0
    for k, b0, L, m in blocks(s, r):
        tau = T-k*Z-b0*beta
        amount = cap_floor(tau, alpha, m+1)
        amount += cap_floor(tau, beta, L-1, 1)
        for b in range(r):
            offset = alpha*(m-(r-1)+b)
            amount += cap_floor(tau-offset, beta+alpha*(r-1), L-1, 1)
        if m == 0 and L >= 2 and tau >= beta:
            amount -= 1
        out += min(2*L-int(m == 0), amount)
    return out


def original_middle_controls():
    counts, fixtures = defaultdict(int), []
    for s, r in ((4, 2), (6, 2), (8, 3), (11, 3), (15, 3), (12, 5)):
        E = ground(s)
        _, evaluate = monic_model(s, r)
        Bdigit, scale = s+1, (s+1)**3
        weights = [tuple(F(scale*x+d) for x, d in zip(w, (1, Bdigit, Bdigit**2)))
                   for w in ((521, 67, 193), (277, 53, 151), (419, 68, 221))]
        for w in weights:
            predicted = middle_basis(s, r, w)
            for p in (0, 2, 3, 5, 7):
                actual, _ = rank_trace(sorted(E, key=lambda a: (dot(a, w), a)),
                                       evaluate, p, set())
                assert actual == predicted
                prices = sorted({dot(a, w) for a in E})
                if s > 8:
                    prices = prices[::max(1, len(prices)//19)]+prices[-1:]
                for T in prices:
                    assert sum(dot(a, w) <= T for a in actual) == middle_rank(s, r, w, T)
                    counts['independent_original_middle_floor_queries'] += 1
                counts['independent_original_middle_field_bases'] += 1
        # Complete tiers at two boundary planes and the outer boundary.
        for w, p in product(((F(3), F(1), F(3)), (F(3), F(1), F(1)),
                             (F(r+1), F(1), F(1))), (0, 2, 3)):
            B, _ = rank_trace(sorted(E, key=lambda a: (dot(a, w), a)), evaluate, p, set())
            for T in {dot(a, w) for a in E}:
                assert sum(dot(a, w) <= T for a in B) == middle_rank(s, r, w, T)
                counts['closed_middle_whole_tier_queries'] += 1
        fixtures.append({'cutoff': s, 'order': r, 'weights': [list(map(str, w)) for w in weights],
                         'characteristics': [0, 2, 3, 5, 7]})
    return dict(counts), fixtures


def entire_interval_controls():
    counts, fixtures = defaultdict(int), []
    for s, r in ((8, 2), (10, 3), (14, 3), (12, 5)):
        reference, original = monic_model(s, r)
        face_reference, face = signed_monic(s, r)
        assert reference == face_reference
        row_degrees = [dot(a, (1, 1-r, 1)) for a in reference]
        laurent_rows, laurent, _ = laurent_model(s, r)
        delta = F(1, 20*(s+1)**2)
        lower, rho = F(0), delta
        total_area, records, seen = F(0), [], set()
        while True:
            seed = (1+r*rho, F(1), F(1))
            t = F(1, 2*r)
            w = tuple(a+t*b for a, b in zip(seed, (1, 1-r, 1)))
            B = middle_basis(s, r, w)
            assert B not in seen
            seen.add(B)
            face_Z, _ = normals(tuple(sorted(B)), face, 0)
            assert all(dot(eta, (1, 1-r, 1)) == 0 for eta in face_Z)
            lo, hi = lower, None
            for x, y, z in face_Z:
                if x > 0:
                    lo = max(lo, -F(y, x))
                elif x < 0:
                    hi = min(hi, -F(y, x)) if hi is not None else -F(y, x)
                else:
                    assert y > 0
            assert lo < rho and (hi is None or rho < hi)
            # Direct small integer determinants of every graded block.
            for k in {dot(a, (1, 1-r, 1)) for a in B}:
                chosen = sorted(a for a in B if dot(a, (1, 1-r, 1)) == k)
                rows = [i for i, row in enumerate(laurent_rows) if row[0] == k]
                assert len(rows) == len(chosen)
                matrix = [[laurent(a).get(i, 0) for a in chosen] for i in rows]
                assert abs(determinant(matrix)) == 1
                counts['middle_graded_unit_determinants'] += 1
            for a in B:
                degree = dot(a, (1, 1-r, 1))
                image = original(a)
                assert all(row_degrees[i] <= degree for i in image)
                assert {i: c for i, c in image.items() if row_degrees[i] == degree} == face(a)
                counts['original_signed_filtered_lifts'] += 1
            original_Z, relations, stats = border_certificate(B, original, s, r, w, save=True)
            assert all(dot(eta, (1, 1-r, 1)) >= 0 for eta in original_Z)
            polygon = [tuple(F(i == j) for i in range(3)) for j in range(3)]
            for eta in face_Z+((1, 0, -1), (0, -1, 1)):
                polygon = clip(polygon, eta, 1)
            assert polygon and area2(polygon) > 0
            assert all(dot(eta, vertex) >= 0 for eta in original_Z for vertex in polygon)
            total_area += area2(polygon)
            counts['entire_original_middle_interval_certificates'] += 1
            counts['middle_original_integer_membership_witnesses'] += stats[0]
            counts['middle_entire_patch_inequalities'] += len(original_Z)*len(polygon)
            records.append({'rho_interval': [str(lo), str(hi) if hi is not None else 'infinity'],
                            'strict_middle_weight': list(map(str, w)), 'standard_set': sorted(B),
                            'patch_vertices': [list(map(str, vertex)) for vertex in polygon],
                            'original_normals': sorted(original_Z),
                            'original_integer_border_relations': relations,
                            'valid_over_every_field': True})
            if hi is None:
                break
            lower, rho = hi, hi+delta
        assert total_area == F(1, 6)
        fixtures.append({'cutoff': s, 'order': r, 'length': len(reference),
                         'normalized_patch_area': str(total_area), 'intervals': records,
                         'scope': 'Entire x>z>y region, over every field, at this cutoff; uniform proof is in the note.'})
        print('Certified entire all-field middle region',s,r,flush=True)
    return dict(counts), fixtures


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    flags = flag_controls()
    births, birth_fixtures = canonical_birth_controls()
    forms, form_fixtures, counterexamples = form_controls()
    originals, original_fixtures = original_middle_controls()
    intervals, interval_fixtures = entire_interval_controls()
    sources = ('docs/middle_interaction_fans.md', 'scripts/run_middle_interaction_checks.py',
               'docs/higher_interaction_fans.md', 'scripts/run_higher_interaction_checks.py',
               'docs/signed_face_reduction.md', 'scripts/run_signed_face_checks.py',
               'docs/nonnegative_graded_fans.md', 'scripts/run_hypersurface_truncation_checks.py',
               'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_socle_stability_checks.py')
    report = {
        'status': 'passed', 'command': 'python3 -B scripts/run_middle_interaction_checks.py',
        'counts': {**flags, **births, **forms, **originals, **intervals},
        'canonical_birth_fixtures': birth_fixtures,
        'noncoprime_form_fixtures': form_fixtures,
        'scope_counterexamples': counterexamples,
        'independent_original_middle_fixtures': original_fixtures,
        'entire_all_field_middle_interval_certificates': interval_fixtures,
        'uniform_results': {
            'laurent_positive_flag': 'rank=min(2*L-int(m=0), cheap retained candidates), all fields, integer unit generic bases',
            'noncoprime_rank_extension': 'The earlier complete-tier formula holds for every nonzero Q,H of binary degrees 2,r+1 at prices (1+lambda,1,1), lambda>=0, without coprimality.',
            'negative_lambda_extension': False,
            'original_middle_region': 'For w_x>w_z>w_y the original standard set is its signed-face standard set, over every field.',
            'middle_region_count': 'At most 12*n+7 original cones, uniformly in s,r and field, r>=2',
            'all_nonlargest_z_regions': 'O(n) original cones altogether; characteristic changes at generic weights can occur only with z strictly largest.',
            'original_middle_rank_query': 'O(r^2*(s+1)) exact arithmetic/floor operations and constant extra storage',
            'remaining_domain': 'After x,y exchange and w_z=1: 0<u<1, 0<lambda<1-u; only 2<=r<sqrt(s) remains as a power-amplification candidate.'},
        'full_signed_transfer_proved': False, 'original_full_family_linear_bound_proved': False,
        'unrestricted_G3_improvement': False, 'novelty_established': False, 'goal_complete': False,
        'scope': 'Uniform positive rank extension and middle-region transfer have proofs. Exact controls are independent supporting evidence. The z-largest region, full original family count, unrestricted G_3 endpoint and current-priority assessment remain open.',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    (ROOT/'results/research/middle_interaction_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
