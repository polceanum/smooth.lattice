#!/usr/bin/env python3
"""Exact truncation and rank controls; no asymptotic fan claim from samples."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import comb
from pathlib import Path
import hashlib
import json

from run_socle_stability_checks import depth_rref
from run_structural_complexity_checks import greedy

ROOT = Path(__file__).resolve().parents[1]


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def poly_mul(f, g, p):
    out = defaultdict(int)
    for a, x in f.items():
        for b, y in g.items():
            out[add(a, b)] += x*y
    return {a: x % p if p else x for a, x in out.items() if (x % p if p else x)}


def simplex(s):
    return tuple(a for a in product(range(s+1), repeat=3) if sum(a) <= s)


def relation_rows(ground, seeds, p):
    index = {a: i for i, a in enumerate(ground)}
    rows = []
    for f in seeds:
        for b in ground:
            row = [0]*len(ground)
            for a, c in f.items():
                if add(a, b) in index:
                    row[index[add(a, b)]] += c
            if any(row):
                rows.append(row)
    return depth_rref(rows, p)[0]


def degree_initial_rows(rows, ground, p):
    # Dense descending degree elimination is independent of the formula.
    order = sorted(range(len(ground)), key=lambda j: (sum(ground[j]), ground[j]), reverse=True)
    echelon, _ = depth_rref(rows, p, order)
    initial = []
    for row in echelon:
        degree = max(sum(a) for a, c in zip(ground, row) if c)
        initial.append([c if sum(a) == degree else 0 for a, c in zip(ground, row)])
    return depth_rref(initial, p)[0]


def homogeneous_formula(g, h, s, p):
    d, e = next(iter(map(sum, g))), next(iter(map(sum, h)))
    assert all(sum(a) == d for a in g) and all(sum(a) == e for a in h) and d < e
    seeds, power = [h], {(0, 0, 0): 1}
    for k in range(1, s//e+2):
        power = poly_mul(power, g, p)
        level = max(0, s+1-e*k)
        for a in simplex(level):
            if sum(a) == level:
                seeds.append({add(a, b): c for b, c in power.items()})
    return seeds


def initial_formula_controls():
    fixtures, counts = [], defaultdict(int)
    pairs = [
        ({(2, 0, 0): 1, (0, 2, 0): 1, (0, 0, 2): 1}, {(1, 1, 1): 1}),
        ({(1, 0, 0): 1, (0, 1, 0): 1}, {(0, 0, 3): 1}),
        ({(2, 0, 0): 1, (0, 2, 0): 1}, {(0, 0, 4): 1}),
    ]
    for label, (g, h) in enumerate(pairs):
        d, e = sum(next(iter(g))), sum(next(iter(h)))
        for p, s in product((0, 2, 3, 5, 101), range(e, 8)):
            ground = simplex(s)
            f = dict(g)
            f.update(h)
            rows = relation_rows(ground, [f], p)
            expected_length = comb(s+3, 3)-comb(s-d+3, 3)
            assert len(ground)-len(rows) == expected_length
            actual = degree_initial_rows(rows, ground, p)
            proposed = relation_rows(ground, homogeneous_formula(g, h, s, p), p)
            assert actual == proposed, 'equality of whole homogeneous initial row spaces'
            hilbert = [0]*(s+1)
            for row in actual:
                degrees = {sum(a) for a, c in zip(ground, row) if c}
                assert len(degrees) == 1
                hilbert[degrees.pop()] += 1
            hilbert = [comb(j+2, 2)-v for j, v in enumerate(hilbert)]
            assert sum(hilbert) == expected_length
            # Direct membership of every monomial tests the canonical cover.
            if label == 0:
                echelon, pivots = depth_rref(rows, p)
                for j in range(len(ground)):
                    row = [int(i == j) for i in range(len(ground))]
                    for v, pivot in zip(echelon, pivots):
                        c = row[pivot]
                        row = [(x-c*y) % p if p else x-c*y for x, y in zip(row, v)]
                    assert any(row)
                    counts['canonical_cover_nonzero_monomials'] += 1
            fixtures.append({'pair': label, 'characteristic': p, 'cutoff': s,
                             'lower_degree': d, 'upper_degree': e,
                             'ambient_size': len(ground), 'length': expected_length,
                             'relation_rank': len(rows), 'global_degree_hilbert_function': hilbert})
            counts['full_initial_row_space_identities'] += 1
            counts['initial_relation_rows'] += len(rows)
    return {'counts': dict(counts), 'fixtures': fixtures,
            'scope': 'Full row-space equality against independent degree elimination, not just equality of lengths. Pairs are coprime over every checked field; finite checks support the uniform regular-pair proof.'}


def compressed_model(s, p):
    """Six free generators over K[U,L], truncated at 3*a+c+i+j<=s."""
    assert p != 2
    basis = tuple((i, j, a, c) for i, j in product(range(2), range(3))
                  for a in range(s//3+1) for c in range(s+1) if i+j+3*a+c <= s)
    index = {v: k for k, v in enumerate(basis)}
    half = pow(2, -1, p) if p else Fraction(1, 2)

    def clean(v):
        return {k: x % p if p else x for k, x in v.items() if (x % p if p else x)}

    def sum_terms(terms):
        out = defaultdict(int)
        for coefficient, vector in terms:
            for k, x in vector.items():
                out[k] += coefficient*x
        return clean(out)

    def reduce_y(i, j, a, c):
        if i+j+3*a+c > s:
            return {}
        if j < 3:
            return {index[(i, j, a, c)]: 1}
        # y^3=L*y^2-L^2*y/2+U*y/2-U.
        return sum_terms([(1, reduce_y(i, j-1, a, c+1)),
                          (-half, reduce_y(i, j-2, a, c+2)),
                          (half, reduce_y(i, j-2, a+1, c)),
                          (-1, reduce_y(i, j-3, a+1, c))])

    def x_column(v):
        i, j, a, c = v
        if i == 0:
            return reduce_y(1, j, a, c)
        # x^2=-x*y-y^2+L*x+L*y-L^2/2+U/2.
        return sum_terms([(-1, reduce_y(1, j+1, a, c)),
                          (-1, reduce_y(0, j+2, a, c)),
                          (1, reduce_y(1, j, a, c+1)),
                          (1, reduce_y(0, j+1, a, c+1)),
                          (-half, reduce_y(0, j, a, c+2)),
                          (half, reduce_y(0, j, a+1, c))])

    X = [x_column(v) for v in basis]
    Y = [reduce_y(i, j+1, a, c) for i, j, a, c in basis]
    L = [reduce_y(i, j, a, c+1) for i, j, a, c in basis]
    U = [reduce_y(i, j, a+1, c) for i, j, a, c in basis]
    Z = [sum_terms([(1, ell), (-1, x), (-1, y)]) for ell, x, y in zip(L, X, Y)]

    def apply(M, vector):
        return sum_terms([(coefficient, M[k]) for k, coefficient in vector.items()])

    def evaluate(a):
        v = {index[(0, 0, 0, 0)]: 1}
        for M, power in zip((X, Y, Z), a):
            for _ in range(power):
                v = apply(M, v)
        return v

    return basis, (X, Y, Z), L, U, apply, evaluate, sum_terms


def compressed_controls():
    counts, fixtures = defaultdict(int), []
    for p, s in product((0, 3, 5, 101), (2, 3, 4, 6, 10)):
        basis, maps, L, U, apply, evaluate, sum_terms = compressed_model(s, p)
        n = len(basis)
        assert n == (s+1)**2
        sparsity = max(len(column) for M in maps for column in M)
        assert sparsity <= 32, 'bounded recurrence support, independent of the cutoff'
        counts['maximum_multiplication_column_support'] = max(counts['maximum_multiplication_column_support'], sparsity)
        for i, j in combinations(range(3), 2):
            for k in range(n):
                assert apply(maps[i], maps[j][k]) == apply(maps[j], maps[i][k])
                counts['commutator_columns'] += 1
        for k, v in enumerate(basis):
            unit = {k: 1}
            g = sum_terms([(1, apply(M, M[k])) for M in maps])
            h = apply(maps[0], apply(maps[1], maps[2][k]))
            assert g == U[k] and sum_terms([(1, g), (1, h)]) == {}
            assert sum_terms([(1, M[k]) for M in maps]) == L[k]
            counts['coordinate_function_identities'] += 3
            for M in maps:
                assert all(sum(basis[j][:2])+3*basis[j][2]+basis[j][3]
                           >= sum(v[:2])+3*v[2]+v[3]+1 for j in M[k])
                counts['strict_filtration_columns'] += 1
            # Independently obtain every module basis vector from original variables.
            i, j, a, c = v
            vector = evaluate((i, j, 0))
            for _ in range(a):
                vector = sum_terms([(1, apply(M, apply(M, vector))) for M in maps])
            for _ in range(c):
                vector = sum_terms([(1, apply(M, vector)) for M in maps])
            assert vector == unit
            counts['cyclic_free_basis_vectors'] += 1
        if s <= 6:
            ground = simplex(s)
            f = {(2, 0, 0): 1, (0, 2, 0): 1, (0, 0, 2): 1, (1, 1, 1): 1}
            rows = relation_rows(ground, [f], p)
            columns = [evaluate(a) for a in ground]
            for row in rows:
                assert not sum_terms([(c, v) for c, v in zip(row, columns) if c])
                counts['independent_macaulay_relations'] += 1
            dense = [[v.get(j, 0) for v in columns] for j in range(n)]
            assert len(depth_rref(dense, p)[0]) == n
            counts['independent_full_kernel_equalities'] += 1
            # Actual original monomials are not each a single free-channel element.
            mixed = next(a for a, v in zip(ground, columns) if len(v) > 1)
        else:
            mixed = (0, 0, 1)
        fixtures.append({'characteristic': p, 'cutoff': s, 'length': n,
                         'free_module_channels': 6, 'ordinary_parent_size': comb(s+3, 3),
                         'maximum_multiplication_column_support': sparsity,
                         'mixed_original_monomial': mixed,
                         'mixed_column': [{'channel_label': basis[k], 'coefficient': str(c)}
                                          for k, c in sorted(evaluate(mixed).items())]})
    return {'counts': dict(counts), 'fixtures': fixtures,
            'scope': 'The six-channel model is checked by commuting cyclic operators, original equations, nilpotent filtration, every free-basis vector and independent full Macaulay kernels. It preserves all ranks of original monomial subsets, not only a list of sampled bases. Characteristic two is covered by the general row-space theorem, not these formulas using 1/2.'}


def saved_fan_controls():
    saved = json.loads((ROOT/'results/research/socle_stability_checks.json').read_text())
    records, count = [], 0
    for fixture in saved['propagating_relation_depth_controls']['actual_multilevel_models']:
        if fixture['characteristic'] not in (3, 5, 101) or fixture['model'] not in (
                'simplex_S3_degree2_seeds1', 'simplex_S4_degree2_seeds1'):
            continue
        if 'complete_positive_fan' not in fixture:
            continue
        s = 3 if 'S3' in fixture['model'] else 4
        p = fixture['characteristic']
        basis, maps, L, U, apply, evaluate, sum_terms = compressed_model(s, p)
        ground = simplex(s)
        columns = {a: tuple(evaluate(a).get(k, 0) for k in range(len(basis))) for a in ground}
        for cell in fixture['complete_positive_fan']:
            chosen = greedy(ground, columns, cell['strict_weight'], p)
            assert chosen == set(map(tuple, cell['standard_basis']))
            count += 1
        records.append({'cutoff': s, 'characteristic': p,
                        'saved_complete_fan_cells_replayed': len(fixture['complete_positive_fan'])})
    return {'original_basis_replays': count, 'fixtures': records,
            'scope': 'Replay existing complete fans in an independently constructed algebra. No larger fan search or extrapolated asymptotic count.'}


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    initial, compressed, fans = initial_formula_controls(), compressed_controls(), saved_fan_controls()
    sources = ('scripts/run_hypersurface_truncation_checks.py', 'docs/hypersurface_truncations.md',
               'scripts/run_socle_stability_checks.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py', 'results/research/socle_stability_checks.json')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_hypersurface_truncation_checks.py',
              'counts': {**initial['counts'], **compressed['counts'],
                         'saved_original_basis_replays': fans['original_basis_replays']},
              'homogeneous_initial_formula': initial, 'six_channel_rank_model': compressed,
              'existing_fan_replays': fans, 'new_asymptotic_fan_bound': False,
              'novelty_established': False, 'goal_complete': False,
              'scope': 'Exact truncation and all-subset rank reductions for a defined deep hypersurface class. The finite-free and graded initial ingredients are classical mechanisms. Neither a full positive-fan count, polynomial separation, unrestricted improvement nor completed priority assessment is asserted.',
              'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}}
    (ROOT/'results/research/hypersurface_truncation_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
