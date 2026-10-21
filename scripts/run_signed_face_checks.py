#!/usr/bin/env python3
"""Exact signed-face ranks, prepared normal forms and finite region certificates."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, gcd
from pathlib import Path
import hashlib
import json

from run_higher_interaction_checks import ground, linear_sum, monic_model, rank_trace
from run_hypersurface_truncation_checks import (
    poly_mul, relation_rows, sparse_basis, sparse_basis_reduce,
)
from run_interleaving_research import minimal_excluded_border
from run_structural_complexity_checks import area2, clip, dot

ROOT = Path(__file__).resolve().parents[1]


def sign(n):
    return -1 if n % 2 else 1


def signed_equation(r):
    return {(2, 0, 0): 1, (0, 0, 2): 1, (1, 1, r): 1}


def signed_monic(s, r):
    basis = tuple(a for a in ground(s) if a[0] < 2)
    index = {a: i for i, a in enumerate(basis)}

    @lru_cache(None)
    def evaluate(a):
        x, y, z = a
        if sum(a) > s:
            return {}
        if x < 2:
            return {index[a]: 1}
        return linear_sum((-1, evaluate(b))
                          for b in ((x-2, y, z+2), (x-1, y+1, z+r)))
    return basis, evaluate


def laurent_model(s, r):
    blocks = {}
    for k in range(-(r-1)*s, s+1):
        b0 = max(0, - (k//(r-1)))
        b1 = (s-k)//r
        if b0 <= b1:
            blocks[k] = b0, b1-b0+1, k+(r-1)*b0
    rows = tuple((k, j, e) for k, (_, L, m) in blocks.items()
                 for j in range(L) for e in (0, 1) if m or j or not e)
    index = {row: i for i, row in enumerate(rows)}
    assert len(rows) == (s+1)**2

    @lru_cache(None)
    def evaluate(a):
        x, y, z = a
        if sum(a) > s:
            return {}
        k = x+z-(r-1)*y
        b0, L, m = blocks[k]
        j = y-b0
        assert 0 <= j < L and 0 <= x <= m+(r-1)*j
        h, e = divmod(x-j, 2)
        out, choose = {}, 1
        for d in range(L-j):
            if choose:
                out[index[k, j+d, e]] = sign(j+h+d)*choose
            choose = choose*(h-d)//(d+1)
        return out
    return rows, evaluate, blocks


def model_controls():
    counts = defaultdict(int)
    for s, r in product(range(2, 11), (2, 3, 4, 5)):
        basis, monic = signed_monic(s, r)
        rows, laurent, blocks = laurent_model(s, r)
        assert len(basis) == len(rows)
        for a in ground(s):
            assert laurent(a) == linear_sum(
                (c, laurent(basis[i])) for i, c in monic(a).items())
            counts['signed_integral_column_equalities'] += 1
        # Integer unit elimination on the reference basis proves more than
        # rational independence: no field is lost in the Laurent row labels.
        pivots = {}
        ordered = sorted(basis, key=lambda a: (
            a[0]+a[2]-(r-1)*a[1], a[1],
            (a[0]-a[1]+blocks[a[0]+a[2]-(r-1)*a[1]][0]) % 2))
        for a in ordered:
            v = dict(laurent(a))
            for j, row in sorted(pivots.items()):
                c = v.get(j, 0)
                if c:
                    v = linear_sum(((1, v), (-c, row)))
            j = min(v)
            assert abs(v[j]) == 1
            pivots[j] = {i: c*v[j] for i, c in v.items()}
            counts['laurent_reference_unit_pivots'] += 1
        for a in ground(s):
            x, y, z = a
            k = x+z-(r-1)*y
            b0, _, _ = blocks[k]
            for w in ((F(3), F(1), F(2)), (F(1), F(1), F(7))):
                alpha, beta = w[0]-w[2], w[1]+(r-1)*w[2]
                assert dot(a, w) == k*w[2]+b0*beta+alpha*x+beta*(y-b0)
                counts['signed_original_price_identities'] += 1
    for s, r, p in product((5, 7), (2, 3, 4), (0, 2, 3, 7)):
        E = ground(s)
        _, evaluate, _ = laurent_model(s, r)
        rows = relation_rows(E, [signed_equation(r)], p)
        assert len(E)-len(rows) == (s+1)**2
        for row in rows:
            image = linear_sum((c, evaluate(a)) for a, c in zip(E, row) if c)
            assert all(not (c % p if p else c) for c in image.values())
            counts['signed_independent_macaulay_rows'] += 1
        counts['signed_full_kernel_equalities'] += 1
    return dict(counts)


def series(N):
    c = [1]
    for i in range(1, N+1):
        c.append(sum(c[j]*c[k]*c[i-1-j-k]
                     for j in range(i) for k in range(i-j)))
    a = [sum(c[j]*c[i-j] for j in range(i+1)) for i in range(N+1)]
    assert all(c[i]*(2*i+1) == comb(3*i, i) for i in range(N+1))
    assert all(a[i]*(i+1) == comb(3*i+1, i) for i in range(N+1))
    return c, a


def add_binary_power(out, coefficient, u, v, z):
    for i in range(u+1):
        out[2*i+v, 2*(u-i)+v, z] += coefficient*comb(u, i)


def prepared_polynomial(s, c, a):
    out = defaultdict(int, {(0, 0, 2): 1})
    for i in range(len(c)):
        if 6*i+2 <= s:
            add_binary_power(out, sign(i)*c[i], i+1, 2*i, 0)
        if 6*i+5 <= s:
            add_binary_power(out, -sign(i)*a[i], i+1, 2*i+1, 1)
    return {alpha: value for alpha, value in out.items() if value and sum(alpha) <= s}


def prepared_model(s, c, a):
    basis = tuple(alpha for alpha in ground(s) if alpha[2] < 2)
    index = {alpha: i for i, alpha in enumerate(basis)}
    tails = prepared_polynomial(s, c, a)
    tails.pop((0, 0, 2), None)

    @lru_cache(None)
    def evaluate(alpha):
        x, y, z = alpha
        if sum(alpha) > s:
            return {}
        if z < 2:
            return {index[alpha]: 1}
        return linear_sum((-v, evaluate((x+b[0], y+b[1], z-2+b[2])))
                          for b, v in tails.items())
    return basis, evaluate


def preparation_controls():
    counts = defaultdict(int)
    c, a = series(16)
    for s in range(2, 25):
        G = prepared_polynomial(s, c, a)
        unit = defaultdict(int, {(0, 0, 0): 1, (1, 1, 1): 1})
        for i in range(len(c)):
            if 6*i+6 <= s:
                add_binary_power(unit, sign(i)*a[i], i+1, 2*i+2, 0)
        replay = {alpha: v for alpha, v in poly_mul(dict(unit), G, 0).items()
                  if sum(alpha) <= s}
        expected = {alpha: 1 for alpha in ((2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 3))
                    if sum(alpha) <= s}
        assert replay == expected
        counts['integer_preparation_factorizations'] += 1
    for s in (2, 4, 6, 8, 10, 12, 15):
        basis, original = monic_model(s, 3)
        zbasis, prepared = prepared_model(s, c, a)
        assert len(basis) == len(zbasis)
        for alpha in ground(s):
            assert prepared(alpha) == linear_sum(
                (v, prepared(basis[i])) for i, v in original(alpha).items())
            counts['prepared_original_integer_column_equalities'] += 1
    return dict(counts), c, a


def cone_index(s, p, c, a):
    last_c = max(i for i in range(len(c))
                 if 6*i+2 <= s and (c[i] % p if p else c[i]))
    last_a = max((i for i in range(len(a))
                  if 6*i+5 <= s and (a[i] % p if p else a[i])), default=None)
    return max(last_c, 2*last_a+1 if last_a is not None else -1)


def extreme_cone_controls(c, a):
    counts, fixtures = defaultdict(int), []
    cases = list(product((2, 5, 8, 11), (0, 2, 3, 7)))+[(17, 7), (23, 3), (29, 11)]
    for s, p in cases:
        _, evaluate = monic_model(s, 3)
        D = cone_index(s, p, c, a)
        Bz = {alpha for alpha in ground(s) if alpha[2] < 2}
        tested = []
        for wx, wy in ((F(2), F(1)), (F(1), F(2)), (F(1), F(1))):
            threshold = (2*D+1)*max(wx, wy)+D*min(wx, wy)
            for side in (-1, 1):
                w = wx, wy, threshold+F(side, 101)
                B, _ = rank_trace(sorted(ground(s), key=lambda alpha: (dot(alpha, w), alpha)),
                                  evaluate, p, set())
                assert (B == Bz) == (side == 1)
                counts['independent_original_extreme_cone_tests'] += 1
                tested.append({'weight': list(map(str, w)), 'z_square_initial': side == 1})
        fixtures.append({'cutoff': s, 'characteristic': p, 'D': D, 'weights': tested})
    # Uniform old witnesses now have the entire elementary initial ideal.
    for K in range(1, 10):
        s = 6*K+5
        assert cone_index(s, 0, c, a) == 2*K+1
        for p in (2, 3, 5, 7, 11, 13, 17, 19):
            if a[K] % p == 0:
                assert cone_index(s, p, c, a) <= 2*K-1
                counts['prime_witness_entire_cone_inequalities'] += 1
    return dict(counts), fixtures


def field_basis(E, evaluate, p):
    if not p:
        return sparse_basis(E, evaluate)
    pivots, B = {}, []
    for alpha in E:
        v, co = field_reduce(evaluate(alpha), pivots, p)
        if not v:
            continue
        j = min(v)
        inv = pow(v[j], -1, p)
        rep = {b: -c*inv % p for b, c in co.items() if c % p}
        rep[alpha] = inv
        pivots[j] = ({k: c*inv % p for k, c in v.items()}, rep)
        B.append(alpha)
    return tuple(B), pivots


def field_reduce(vector, pivots, p):
    if not p:
        return sparse_basis_reduce(vector, pivots)
    v = {j: c % p for j, c in vector.items() if c % p}
    co = defaultdict(int)
    for j, (row, rep) in sorted(pivots.items()):
        t = v.get(j, 0)
        if not t:
            continue
        for k, c in row.items():
            v[k] = (v.get(k, 0)-t*c) % p
            if not v[k]:
                del v[k]
        for b, c in rep.items():
            co[b] = (co[b]+t*c) % p
    return v, {b: c for b, c in co.items() if c % p}


def normals(B, evaluate, p):
    selected, pivots = field_basis(B, evaluate, p)
    assert len(selected) == len(B)
    out, relations = set(), []
    for alpha in minimal_excluded_border(tuple(B)):
        rem, co = field_reduce(evaluate(alpha), pivots, p)
        assert not rem
        polynomial = {alpha: 1}
        polynomial.update({b: -c % p if p else -c for b, c in co.items()})
        replay = linear_sum((c, evaluate(b)) for b, c in polynomial.items())
        assert all(not (c % p if p else c) for c in replay.values())
        for b in co:
            eta = tuple(x-y for x, y in zip(alpha, b))
            divisor = gcd(*eta)
            out.add(tuple(x//divisor for x in eta))
        relations.append({'head': alpha,
                          'tails': [{'exponent': b, 'coefficient': str(c)}
                                    for b, c in sorted(co.items())]})
    return tuple(sorted(out)), relations


def positive_weight(rho, r):
    beta = F(r) if rho >= 0 else (F(r-1)-1/rho)/2
    w = 1+rho*beta, beta-(r-1), F(1)
    assert min(w) > 0
    return w


def finite_transfer_controls():
    counts, fixtures = defaultdict(int), []
    # Full interval walks and entire rational patches at these specified
    # cutoffs/fields. No extrapolation to other cutoffs or characteristics.
    for s, r, p in ((4, 2, 0), (6, 2, 0), (6, 2, 2),
                    (8, 3, 0), (8, 3, 3), (11, 3, 2), (17, 3, 7)):
        _, face = signed_monic(s, r)
        _, original = monic_model(s, r)
        E = ground(s)
        delta = F(1, 20*(s+1)**2*(r+1))
        lower, rho = -F(1, r-1), -F(1, r-1)+delta
        patches, visited = [], set()
        total_area = F(0)
        while True:
            w = positive_weight(rho, r)
            B, _ = field_basis(sorted(E, key=lambda alpha: (dot(alpha, w), alpha)), face, p)
            key = frozenset(B)
            assert key not in visited
            visited.add(key)
            Z, _ = normals(B, face, p)
            assert all(dot(eta, (1, 1-r, 1)) == 0 for eta in Z)
            lo, hi = lower, None
            for x, y, z in Z:
                if x > 0:
                    lo = max(lo, -F(y, x))
                elif x < 0:
                    hi = min(hi, -F(y, x)) if hi is not None else -F(y, x)
                else:
                    assert y > 0
            assert lo < rho and (hi is None or rho < hi)
            original_Z, relations = normals(B, original, p)
            assert all(dot(eta, (1, 1-r, 1)) >= 0 for eta in original_Z)
            polygon = [tuple(F(i == j) for i in range(3)) for j in range(3)]
            for eta in Z+((1, -1, 0), (0, -1, 1)):
                polygon = clip(polygon, eta, 1)
            assert polygon and area2(polygon) > 0
            assert all(dot(eta, vertex) >= 0 for eta in original_Z for vertex in polygon)
            total_area += area2(polygon)
            counts['finite_signed_face_intervals'] += 1
            counts['finite_original_border_relations'] += len(relations)
            counts['finite_entire_patch_inequalities'] += len(original_Z)*len(polygon)
            patches.append({'rho_interval': [str(lo), str(hi) if hi is not None else 'infinity'],
                            'face_standard_set': sorted(B),
                            'original_normals': original_Z,
                            'patch_vertices': [list(map(str, vertex)) for vertex in polygon],
                            'original_border_relations': relations})
            if hi is None:
                break
            lower, rho = hi, hi+delta
        assert total_area == F(1, 3)
        assert len(visited) <= 12*(s+1)**2+7
        counts['finite_complete_cheapest_y_regions'] += 1
        fixtures.append({'cutoff': s, 'order': r, 'characteristic': p,
                         'length': (s+1)**2, 'normalized_patch_area': str(total_area),
                         'intervals': patches,
                         'scope': 'Entire open region with y strictly cheapest at generic weights, for this specified cutoff and field only; boundary vertices are closed-patch certificates.'})
        print('Certified finite cheapest-y region',s,r,p,flush=True)
    return dict(counts), fixtures


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    models = model_controls()
    prepared, c, a = preparation_controls()
    extremes, cone_fixtures = extreme_cone_controls(c, a)
    transfer, transfer_fixtures = finite_transfer_controls()
    sources = ('docs/signed_face_reduction.md', 'scripts/run_signed_face_checks.py',
               'docs/nonnegative_graded_fans.md', 'docs/higher_interaction_fans.md',
               'docs/fixed_equation_bad_primes.md', 'scripts/run_higher_interaction_checks.py',
               'scripts/run_hypersurface_truncation_checks.py',
               'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_socle_stability_checks.py')
    report = {
        'status': 'passed', 'command': 'python3 -B scripts/run_signed_face_checks.py',
        'counts': {**models, **prepared, **extremes, **transfer},
        'extreme_cone_fixtures': cone_fixtures,
        'finite_signed_transfer_certificates': transfer_fixtures,
        'uniform_results': {
            'signed_face': 'in_(1,1-r,1)(I_s,r)=(x^2+z^2+x*y*z^r)+m^(s+1), all fields, r>=2',
            'laurent_block_dimensions': '2*L, or 2*L-1 when m=0; integral unit reference matrix',
            'signed_face_fan_bound': '12*n+7 by the previous signed-grading theorem',
            'fixed_cubic_preparation': 'F=(V*z+C^(-1))*(z^2-U*V*C^2*z+U*C), T=-U*V^2, C=1+T*C^3',
            'extreme_z_square_cone': 'w_z>(2*D+1)*max(w_x,w_y)+D*min(w_x,w_y), with characteristic-dependent last surviving indices',
            'prime_witness_consequence': 'All selected bad characteristics have the same elementary initial (z^2)+m^(s+1) at the CRT witness, while the rational initial differs.'},
        'uniform_signed_transfer_proved': False,
        'original_global_linear_fan_bound_proved': False,
        'unrestricted_G3_improvement': False, 'novelty_established': False,
        'goal_complete': False,
        'scope': 'The exact rank models and extreme cone have uniform proofs. Entire signed-transfer patches are certified only for the saved finite cutoffs/fields; neither this computation nor equal-plane samples prove the remaining uniform signed cone criterion.',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    (ROOT/'results/research/signed_face_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
