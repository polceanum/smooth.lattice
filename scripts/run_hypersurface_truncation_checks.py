#!/usr/bin/env python3
"""Exact truncation and rank controls; no asymptotic fan claim from samples."""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb
from pathlib import Path
import hashlib
import json

from run_socle_stability_checks import depth_rref
from run_structural_complexity_checks import dot, greedy
from run_interleaving_research import minimal_excluded_border

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


def integral_cubic_model(s, z_square=1):
    """Monic division over Z, with x-exponent <2 and total degree <=s."""
    basis = tuple(a for a in simplex(s) if a[0] < 2)
    index = {a: i for i, a in enumerate(basis)}

    @lru_cache(None)
    def evaluate(a):
        x, y, z = a
        if x+y+z > s:
            return {}
        if x < 2:
            return {index[a]: 1}
        out = defaultdict(int)
        for b, c in (((x-2, y+2, z), -1),
                     ((x-2, y, z+2), -z_square),
                     ((x-1, y+1, z+1), -1)):
            for k, value in evaluate(b).items():
                out[k] += c*value
        return {k: value for k, value in out.items() if value}

    return basis, evaluate


def sparse_basis_reduce(vector, pivots):
    """Exact column elimination, retaining a representation in chosen columns."""
    vector = {k: Fraction(c) for k, c in vector.items() if c}
    coefficients = defaultdict(Fraction)
    for k, (row, representation) in sorted(pivots.items()):
        value = vector.get(k, 0)
        if not value:
            continue
        for j, c in row.items():
            vector[j] = vector.get(j, 0)-value*c
            if not vector[j]:
                del vector[j]
        for a, c in representation.items():
            coefficients[a] += value*c
    return vector, {a: c for a, c in coefficients.items() if c}


def sparse_basis(points, evaluate):
    pivots, chosen = {}, []
    for a in points:
        vector, coefficients = sparse_basis_reduce(evaluate(a), pivots)
        if not vector:
            continue
        k = min(vector)
        value = vector[k]
        representation = {b: -c/value for b, c in coefficients.items()}
        representation[a] = 1/value
        pivots[k] = ({j: c/value for j, c in vector.items()}, representation)
        chosen.append(a)
    return tuple(chosen), pivots


def integral_membership_witness(polynomial, s, z_square):
    """Return Q with polynomial-f*Q in m^(s+1), checking multiplication anew."""
    work, quotient = dict(polynomial), defaultdict(int)
    f = {(2, 0, 0): 1, (0, 2, 0): 1, (1, 1, 1): 1}
    if z_square:
        f[(0, 0, 2)] = z_square
    while True:
        a = max((a for a, c in work.items() if c and a[0] >= 2 and sum(a) <= s), default=None)
        if a is None:
            break
        value = work[a]
        b = (a[0]-2, a[1], a[2])
        quotient[b] += value
        for exponent, c in f.items():
            v = add(b, exponent)
            work[v] = work.get(v, 0)-value*c
            if not work[v]:
                del work[v]
    assert all(sum(a) > s for a in work), 'integer ideal membership, including cutoff terms'
    reconstructed = poly_mul(f, quotient, 0)
    for a in set(polynomial) | set(reconstructed):
        if sum(a) <= s:
            assert polynomial.get(a, 0) == reconstructed.get(a, 0)
    return {a: c for a, c in quotient.items() if c}


def integral_border_relations(chosen, pivots, evaluate, s, z_square):
    normals, relations, coefficients, witness_terms = set(), 0, 0, 0
    chosen_set = set(chosen)
    assert all(tuple(x-(i == j) for i, x in enumerate(a)) in chosen_set
               for a in chosen for j in range(3) if a[j])
    for a in minimal_excluded_border(chosen):
        remainder, representation = sparse_basis_reduce(evaluate(a), pivots)
        assert not remainder
        assert all(c.denominator == 1 for c in representation.values()), 'monic integer certificate'
        polynomial = {a: 1}
        for b, c in representation.items():
            polynomial[b] = -int(c)
            normal = tuple(x-y for x, y in zip(a, b))
            assert any(normal)
            normals.add(normal)
        witness = integral_membership_witness(polynomial, s, z_square)
        relations += 1
        coefficients += len(representation)
        witness_terms += len(witness)
    return normals, (relations, coefficients, witness_terms)


def integral_model_controls():
    counts, fixtures = defaultdict(int), []
    for s, z_square in product((2, 3, 4, 6, 10), (0, 1)):
        basis, evaluate = integral_cubic_model(s, z_square)
        assert len(basis) == (s+1)**2
        maps = [[evaluate(add(a, tuple(int(i == j) for i in range(3))))
                 for a in basis] for j in range(3)]
        assert max(len(v) for M in maps for v in M) <= 3
        counts['integral_operator_columns'] += 3*len(basis)
        for p in (0, 2, 3, 5, 101):
            ground = simplex(s)
            if s <= 6:
                f = {(2, 0, 0): 1, (0, 2, 0): 1, (1, 1, 1): 1}
                if z_square:
                    f[(0, 0, 2)] = 1
                rows = relation_rows(ground, [f], p)
                assert len(ground)-len(rows) == len(basis)
                for row in rows:
                    result = defaultdict(int)
                    for a, c in zip(ground, row):
                        for k, value in evaluate(a).items():
                            result[k] += c*value
                    assert not any(value % p if p else value for value in result.values())
                    counts['integral_independent_macaulay_relations'] += 1
                # The original standard monomials themselves give an identity submatrix.
                assert all(evaluate(a) == {k: 1} for k, a in enumerate(basis))
                counts['integral_independent_full_kernel_equalities'] += 1
            if z_square and p != 2:
                free, _, _, _, _, previous, _ = compressed_model(s, p)
                transfer = [previous(a) for a in basis]
                for a in ground:
                    result = defaultdict(int)
                    for k, c in evaluate(a).items():
                        for j, value in transfer[k].items():
                            result[j] += c*value
                    result = {j: (value % p if p else value)
                              for j, value in result.items() if (value % p if p else value)}
                    assert result == previous(a)
                    counts['old_and_integral_column_equalities'] += 1
        fixtures.append({'cutoff': s, 'z_square_coefficient': z_square,
                         'length': len(basis), 'coefficient_ring': 'Z',
                         'maximum_operator_column_support': max(len(v) for M in maps for v in M)})
    return {'counts': dict(counts), 'fixtures': fixtures,
            'scope': 'A classical monic-division basis gives a smaller integral model. Full independent Macaulay kernels cover both cubics over Q and F_2,F_3,F_5,F_101; explicit column maps preserve the previous six-channel model. No cost-ray or new asymptotic theorem follows from this basis alone.'}


def cubic_boundary_certificates():
    counts, fixtures = defaultdict(int), []
    for s in (3, 4, 5, 6, 8, 10, 12, 14):
        ground = simplex(s)
        basis, boundary = integral_cubic_model(s, 0)
        _, original = integral_cubic_model(s, 1)
        epsilon = Fraction(1, 10*(s+1)**3)
        lam, sectors, all_bases = epsilon, [], set()
        while lam < 1:
            weight = (s+1+lam, s+1, Fraction(1))
            chosen, pivots = sparse_basis(sorted(ground, key=lambda a: (dot(weight, a), a)), boundary)
            assert len(chosen) == len(basis)
            chosen_original, original_pivots = sparse_basis(chosen, original)
            assert chosen_original == chosen
            z0, stats0 = integral_border_relations(chosen, pivots, boundary, s, 0)
            z1, stats1 = integral_border_relations(chosen, original_pivots, original, s, 1)
            lo, hi = Fraction(0), Fraction(1)
            for x, y, z in z0:
                assert x+y == 0, 'boundary grading is preserved exactly'
                assert (x+z) % 2 == 0, 'boundary parity grading is preserved'
                if x > 0:
                    lo = max(lo, Fraction(-z, x))
                elif x < 0:
                    hi = min(hi, Fraction(-z, x))
                else:
                    assert z > 0
            assert lo < lam < hi
            assert not sectors or lo == Fraction(sectors[-1]['lambda_interval'][1])
            endpoint_witnesses = []
            for end in (lo, hi):
                value, normal = min((x+y+x*end+z, (x, y, z)) for x, y, z in z1)
                endpoint_witnesses.append({'lambda': str(end), 'minimum_slack': str(value),
                                           'normal_attaining_minimum': normal})
            for x, y, z in z1:
                assert x+y >= 0, 'every lifted cone contains the boundary ray'
                for end in (lo, hi):
                    value = x+y+x*end+z
                    assert value >= 0, 'entire ordered sector, at u=1 and both lambda endpoints'
                counts['whole_sector_normal_checks'] += 1
            for perm in permutations(range(3)):
                all_bases.add(frozenset(tuple(a[i] for i in perm) for a in chosen))
            for prefix, stats in (('boundary', stats0), ('original', stats1)):
                for name, value in zip(('integer_border_relations', 'integer_tail_coefficients', 'integer_quotient_witness_terms'), stats):
                    counts[prefix+'_'+name] += value
            sectors.append({'lambda_interval': [str(lo), str(hi)],
                            'strict_weight': [str(v) for v in (Fraction(2)+(lo+hi)/2, Fraction(2), Fraction(1))],
                            'standard_basis_border': sorted(minimal_excluded_border(chosen)),
                            'boundary_normal_count': len(z0),
                            'lambda_lower_facet_witnesses': sorted((x, y, z) for x, y, z in z0
                                                                  if x > 0 and Fraction(-z, x) == lo),
                            'lambda_upper_facet_witnesses': sorted((x, y, z) for x, y, z in z0
                                                                  if x < 0 and Fraction(-z, x) == hi),
                            'original_normal_count': len(z1),
                            'original_normals_sha256': hashlib.sha256(json.dumps(sorted(z1)).encode()).hexdigest(),
                            'minimum_boundary_ray_slack': min(x+y for x, y, z in z1),
                            'endpoint_witnesses': endpoint_witnesses,
                            'monic_integer_certificates': True,
                            'original_integer_border_relations': stats1[0]})
            lam = hi+epsilon
        assert sectors[0]['lambda_interval'][0] == '0' and sectors[-1]['lambda_interval'][1] == '1'
        # The three outer cones have exactly the elementary bases (x_i^2,m^(s+1)).
        for i in range(3):
            all_bases.add(frozenset(a for a in ground if a[i] < 2))
        fixtures.append({'cutoff': s, 'length': len(basis), 'coefficient_ring': 'Z',
                         'valid_over_every_field': True,
                         'exact_interval_walk_step': str(epsilon),
                         'complete_positive_fan_cells': len(all_bases),
                         'ordered_central_sectors': sectors})
        counts['all_field_complete_positive_fans'] += 1
        counts['certified_ordered_central_sectors'] += len(sectors)
        counts['complete_positive_fan_cells'] += len(all_bases)
    return {'counts': dict(counts), 'fixtures': fixtures,
            'scope': 'Complete continuous positive fans at the eight stated cutoffs, over every field: monic integer border relations, exact integer membership witnesses, rational cone endpoints and symmetry/outer-cone coverage. No mesh or selected-characteristic inference is used. The uniform boundary-cover conjecture and an asymptotic fan count for the cubic remain unproved.'}


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
    integral, boundary = integral_model_controls(), cubic_boundary_certificates()
    sources = ('scripts/run_hypersurface_truncation_checks.py', 'docs/hypersurface_truncations.md',
               'scripts/run_socle_stability_checks.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py', 'results/research/socle_stability_checks.json')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_hypersurface_truncation_checks.py',
              'counts': {**initial['counts'], **compressed['counts'],
                         'saved_original_basis_replays': fans['original_basis_replays'],
                         **integral['counts'], **boundary['counts']},
              'homogeneous_initial_formula': initial, 'six_channel_rank_model': compressed,
              'existing_fan_replays': fans, 'integral_two_channel_rank_model': integral,
              'all_field_cubic_boundary_certificates': boundary,
              'new_asymptotic_fan_bound': False,
              'novelty_established': False, 'goal_complete': False,
              'scope': 'Exact truncation and all-subset rank reductions, plus complete positive fans at eight stated cubic cutoffs over every field by monic integer certificates. The finite-free and monic-division ingredients are classical mechanisms. No uniform full-fan count, polynomial separation, unrestricted improvement or completed priority assessment is asserted.',
              'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}}
    (ROOT/'results/research/hypersurface_truncation_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
