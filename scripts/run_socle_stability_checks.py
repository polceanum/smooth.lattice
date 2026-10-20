#!/usr/bin/env python3
"""Exact all-field socle constructions and sharp relative-cone stability controls."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from math import gcd, lcm, prod
from pathlib import Path
from random import Random
import hashlib
import json

from run_interleaving_research import minimal_excluded_border
from run_initial_cone_height_checks import clean, normal_form
from run_algebraic_incidence_controls import pair_checks
from run_structural_complexity_checks import area2, clip, cross, determinant, dot, greedy, hyperbolic_cross
from run_extremal_initial_ideals_checks import integer_weight

ROOT = Path(__file__).resolve().parents[1]


def endpoints(d, k):
    if d == 2:
        return ((k+2, 0), (0, k+1)), ((2, k-1), (k+1, 1))
    return ((k+2, 0, 0), (0, k, 0), (0, 0, 4*k+3)), ((2, k-1, 1), (k+1, 0, 1), (1, 1, 2))


def model(d, k, padding=0):
    heads, tails = endpoints(d, k)
    tops = heads + tails
    ambient = set().union(*(set(product(*(range(a+1) for a in v))) for v in tops))
    if padding:
        if d == 2:
            assert k >= 7 and 0 <= padding <= 4
            ambient.update((3, j) for j in range(2, padding+2))
        else:
            assert k >= 12 and 0 <= padding <= 11
            ambient.update((3, j, 0) for j in range(1, padding+1))
    assert all(not any(tuple(a+(i == j) for j, a in enumerate(v)) in ambient
                       for i in range(d)) for v in tops)
    assert all(tuple(a-(i == j) for j, a in enumerate(v)) in ambient
               for v in ambient for i in range(d) if v[i])
    expected = (5*k+1 if d == 2 else 12*k+5) + padding
    assert len(ambient) == expected
    reference = tuple(sorted(ambient - set(heads)))
    index = {a: i for i, a in enumerate(reference)}
    replacements = dict(zip(heads, tails))

    def image(a):
        return index.get(replacements.get(a, a), -1)

    maps = [tuple(image(tuple(a+(i == j) for j, a in enumerate(v))) for v in reference) for i in range(d)]
    for i, j in combinations(range(d), 2):
        assert all((-1 if maps[j][v] < 0 else maps[i][maps[j][v]]) ==
                   (-1 if maps[i][v] < 0 else maps[j][maps[i][v]]) for v in range(len(reference)))

    def evaluate(a, start=None):
        v = index[(0,)*d] if start is None else start
        for i, power in enumerate(a):
            for _ in range(power):
                if v < 0:
                    return -1
                v = maps[i][v]
        return v

    assert all(evaluate(a) == image(a) for a in ambient)
    assert all(evaluate(a) == i for i, a in enumerate(reference))
    return heads, tails, ambient, reference, image, maps, evaluate


def p_mul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def symbolic_matrices():
    matrices = [
        (((0, 1), (1, -1)), ((-1, -1), (0, 1))),
        (((0, 1), (1, -1), (-1,)),
         ((-1, -1), (0, 1), (-1,)),
         ((-1,), (-1,), (1, 4))),
    ]
    inverses = [
        (((0, 1), (-1, 1)), ((1, 1), (0, 1))),
        (((-1, 1, 4), (0, -3, 4), (-1, 2)),
         ((2, 5, 4), (-1, 1, 4), (1, 2)),
         ((1, 2), (-1, 2), (1,))),
    ]
    checks, records = 0, []
    for M, V in zip(matrices, inverses):
        d = len(M)
        for A, B in ((M, V), (V, M)):
            for i, j in product(range(d), repeat=2):
                coefficients = [0]*4
                for t in range(d):
                    for degree, value in enumerate(p_mul(A[i][t], B[t][j])):
                        coefficients[degree] += value
                assert coefficients == [int(i == j), 0, 0, 0]
                checks += 1
        det = [0]*4
        for perm in permutations(range(d)):
            term = [(-1)**sum(perm[i] > perm[j] for i in range(d) for j in range(i+1, d))]
            for i in range(d):
                term = p_mul(term, M[i][perm[i]])
            for degree, value in enumerate(term):
                det[degree] += value
        assert det == [1, 0, 0, 0]
        records.append({'dimension': d, 'matrix_coefficients_ascending': M, 'inverse_coefficients_ascending': V,
                        'determinant_coefficients': det})
    return {'matrix_product_entry_identities': checks, 'polynomial_certificates': records}


def family_controls():
    counts, fixtures = defaultdict(int), []
    for d, k in product((2, 3), (3, 4, 6, 8)):
        heads, tails, ambient, reference, image, maps, evaluate = model(d, k)
        n = len(reference)
        M = [tuple(a-b for a, b in zip(u, v)) for u, v in zip(heads, tails)]
        w = (2*k-1, 2*k+1) if d == 2 else (8*k*k-2, 8*k*k+8*k+2, 4*k+1)
        assert determinant(M) == 1 and [dot(z, w) for z in M] == [1]*d
        border = minimal_excluded_border(reference)
        G = [{u: 1, v: -1} for u, v in zip(heads, tails)] + [{a: 1} for a in border if a not in heads]
        assert all(len(g) <= 2 for g in G)
        ground = tuple(sorted(set(hyperbolic_cross(n, d)) | set(border)))
        fields = []
        for p in (0, 2, 3, 5, 101):
            pairs = pair_checks(G, w, p)
            for a in ground:
                expected = {} if image(a) < 0 else {reference[image(a)]: 1}
                assert normal_form({a: 1}, [clean(g, p) for g in G], w, p) == expected
                counts['ordinary_monomial_normal_forms'] += 1
            counts['ordinary_s_pairs'] += pairs
            fields.append({'characteristic': p, 's_pairs': pairs, 'ground_normal_forms': len(ground)})
        # Relations annihilate each basis vector of the integer module.
        for g in G:
            for start in range(n):
                column = defaultdict(int)
                for a, c in g.items():
                    v = evaluate(a, start)
                    if v >= 0:
                        column[v] += c
                assert not any(column.values())
                counts['integer_relation_columns'] += 1
        counts['commuting_multiplication_models'] += 1
        if d == 3:
            cells = []
            for signs in product((-1, 1), repeat=3):
                polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
                for z, sign in zip(M, signs):
                    polygon = clip(polygon, z, sign)
                if not polygon or not area2(polygon):
                    assert signs == (-1, -1, -1)
                    continue
                weight = integer_weight(tuple(sum(v[i] for v in polygon)/len(polygon) for i in range(3)))
                chosen_heads = tuple(u if sign > 0 else v for u, v, sign in zip(heads, tails, signs))
                chosen_tails = tuple(v if sign > 0 else u for u, v, sign in zip(heads, tails, signs))
                B = tuple(sorted(ambient - set(chosen_heads)))
                H = [{u: 1, v: -1} for u, v in zip(chosen_heads, chosen_tails)]
                H += [{a: 1} for a in minimal_excluded_border(B) if a not in chosen_heads]
                counts['all_cone_rational_s_pairs'] += pair_checks(H, weight, 0)
                if k <= 4:
                    columns = {a: tuple(int(i == image(a)) for i in range(n)) for a in ground}
                    assert greedy(ground, columns, weight, 101) == set(B)
                    counts['independent_gaussian_initial_bases'] += 1
                cells.append({'signs': signs, 'strict_weight': weight, 'simplex_double_area': str(area2(polygon)),
                              'standard_basis': B, 'reduced_basis_heads': sorted(g for g in minimal_excluded_border(B))})
            assert len(cells) == 7 and sum(Fraction(c['simplex_double_area']) for c in cells) == 1
            counts['complete_seven_cone_fans'] += 1
        else:
            low, high = Fraction(k-1, k), Fraction(k, k+1)
            cells = []
            for left, right, ratio in ((0, low, low/2), (low, high, (low+high)/2),
                                       (high, None, high+1)):
                weight = integer_weight((ratio, Fraction(1)))
                signs = tuple(1 if dot(z, weight) > 0 else -1 for z in M)
                assert signs != (-1, -1)
                chosen_heads = tuple(u if sign > 0 else v for u, v, sign in zip(heads, tails, signs))
                chosen_tails = tuple(v if sign > 0 else u for u, v, sign in zip(heads, tails, signs))
                B = tuple(sorted(ambient-set(chosen_heads)))
                H = [{u: 1, v: -1} for u, v in zip(chosen_heads, chosen_tails)]
                H += [{a: 1} for a in minimal_excluded_border(B) if a not in chosen_heads]
                counts['all_cone_rational_s_pairs'] += pair_checks(H, weight, 0)
                if k <= 4:
                    columns = {a: tuple(int(i == image(a)) for i in range(n)) for a in ground}
                    assert greedy(ground, columns, weight, 101) == set(B)
                    counts['independent_gaussian_initial_bases'] += 1
                cells.append({'ratio_interval': [str(left), None if right is None else str(right)],
                              'strict_weight': weight, 'signs': signs, 'standard_basis': B})
            assert len({tuple(row['signs']) for row in cells}) == 3
            assert low < Fraction(w[0], w[1]) < high
            counts['complete_three_cone_fans'] += 1
        fixtures.append({'dimension': d, 'k': k, 'length': n, 'ambient_downset_size': len(ambient),
                         'heads': heads, 'tails': tails, 'matrix': M, 'integer_witness': w,
                         'field_checks': fields, 'positive_initial_ideals': 3 if d == 2 else 7,
                         'exact_cells': cells})
    return {'counts': dict(counts), 'exact_family_fixtures': fixtures}


def asymptotic_controls():
    witnesses, padding = [], []
    for k in (3, 4, 10, 100, 10000, 10**12):
        M = ((k, 1-k, -1), (-k-1, k, -1), (-1, -1, 4*k+1))
        V = ((4*k*k+k-1, 4*k*k-3*k, 2*k-1),
             (4*k*k+5*k+2, 4*k*k+k-1, 2*k+1), (2*k+1, 2*k-1, 1))
        assert all(dot(row, tuple(V[t][j] for t in range(3))) == int(i == j)
                   for i, row in enumerate(M) for j in range(3))
        assert all(x > 0 for row in V for x in row)
        w = tuple(sum(row) for row in V)
        denominators = tuple(sum(abs(a)*b for a, b in zip(row, w)) for row in M)
        assert denominators == (16*k**3-4*k-1, 16*k**3+16*k*k+4*k-1, 32*k*k+16*k+1)
        assert max(denominators) == denominators[1]
        dual = (2*k+1, 2*k-1, 1)
        assert tuple(sum(dual[i]*M[i][j] for i in range(3)) for j in range(3)) == (0, 0, 1)
        dual_abs = tuple(sum(dual[i]*abs(M[i][j]) for i in range(3)) for j in range(3))
        assert dual_abs == (4*k*k+2*k, 4*k*k-2*k, 8*k+1)
        assert dot(dual_abs, (2*k-1, 2*k+1, 1)) == 16*k**3+4*k+1
        rays = [tuple(Fraction(V[i][j], V[2][j]) for i in range(3)) for j in range(3)]
        twice_area = abs(determinant([tuple(v[i] for v in rays) for i in range(3)]))
        assert twice_area == Fraction(1, 4*k*k-1)
        assert all(v[0] >= 2*k-1 and v[1] >= 2*k+1 and v[0]+v[1] <= 4*k+1 for v in rays)
        lo, hi = Fraction(1, max(denominators)), Fraction(1, 16*k**3+4*k+1)
        assert 0 < lo < hi
        witnesses.append({'k': k, 'length': 12*k+2, 'positive_primitive_rays': list(zip(*V)),
                          'minimum_integer_height': max(w), 'relative_margin_lower': str(lo),
                          'relative_margin_upper': str(hi),
                          'log_area_lower': str(Fraction(2, (4*k*k-1)*(4*k+1)**2)),
                          'log_area_upper': str(Fraction(1, 2*(4*k*k-1)**2)),
                          'planar_length': 5*k-1, 'planar_relative_margin_lower': str(Fraction(1, 4*k*k)),
                          'planar_relative_margin_upper': str(Fraction(1, 4*(k*k-1))),
                          'planar_log_width_lower': str(Fraction(1, k*k)),
                          'planar_log_width_upper': str(Fraction(1, k*k-1))})
    for d, k in ((2, 7), (2, 12), (3, 12), (3, 16)):
        for q in range(5 if d == 2 else 12):
            _, _, ambient, reference, _, _, _ = model(d, k, q)
            n = (5*k-1 if d == 2 else 12*k+2) + q
            assert len(reference) == n
            padding.append({'dimension': d, 'k': k, 'padding': q, 'length': n,
                            'ambient_downset_size': len(ambient), 'unchanged_positive_cones': 3 if d == 2 else 7})
    return {'symbolic_scale_witnesses': witnesses, 'all_length_padding_checks': padding}


def socle_budget_controls():
    states = {frozenset({(0, 0, 0)})}
    checked, chains, commutators, examples = 0, 0, 0, []
    for size in range(1, 11):
        for D in states:
            T = {a for a in D if all(tuple(x+(i == j) for j, x in enumerate(a)) not in D for i in range(3))}
            core = D - T
            for i in range(3):
                occupied = set()
                for a in T:
                    chain = {tuple(x-(j == i)*step for j, x in enumerate(a)) for step in range(1, a[i]+1)}
                    assert chain <= core and not occupied & chain
                    occupied.update(chain)
                    chains += len(chain)
                assert sum(a[i] for a in T) <= len(core)
            checked += 1
        following = set()
        for D in states:
            for a in set(minimal_excluded_border(tuple(D))):
                following.add(D | {a})
        states = following
    # Arbitrary linear identifications in an irregular six-dimensional socle.
    _, _, D, _, _, _, _ = model(3, 3)
    T = sorted(a for a in D if all(tuple(x+(i == j) for j, x in enumerate(a)) not in D for i in range(3)))
    core = D - set(T)
    for p, rank in product((2, 3, 101), range(1, 6)):
        random = Random(p*100+rank)
        columns = [tuple(int(i == j) for i in range(rank)) if j < rank
                   else tuple(random.randrange(p) for _ in range(rank)) for j in range(len(T))]
        fixed = tuple(sorted(core))
        reference = fixed + tuple(T[:rank])
        indices, top_indices = {a: i for i, a in enumerate(fixed)}, {a: i for i, a in enumerate(T)}
        q, n = len(fixed), len(reference)

        def residue(a):
            if a in indices:
                return tuple(int(i == indices[a]) for i in range(n))
            if a in top_indices:
                return (0,)*q + columns[top_indices[a]]
            return (0,)*n

        maps = [tuple(residue(tuple(x+(i == j) for j, x in enumerate(a))) for a in reference) for i in range(3)]

        def apply(i, vector):
            return tuple(sum(c*maps[i][j][row] for j, c in enumerate(vector) if c) % p for row in range(n))

        for i, j in combinations(range(3), 2):
            for v in range(n):
                assert apply(i, maps[j][v]) == apply(j, maps[i][v])
                commutators += 1
        for a in D:
            vector = residue((0, 0, 0))
            for i, power in enumerate(a):
                for _ in range(power):
                    vector = apply(i, vector)
            assert vector == residue(a)
        bases = []
        for B in combinations(range(len(T)), rank):
            # Direct exact modular row elimination on each proposed socle basis.
            A = [[columns[j][i] for j in B] for i in range(rank)]
            full = True
            for j in range(rank):
                pivot = next((i for i in range(j, rank) if A[i][j] % p), None)
                if pivot is None:
                    full = False
                    break
                A[pivot], A[j] = A[j], A[pivot]
                inv = pow(A[j][j], -1, p)
                for i in range(j+1, rank):
                    scale = A[i][j]*inv % p
                    A[i] = [(x-scale*y) % p for x, y in zip(A[i], A[j])]
            if full:
                profile = tuple(sum(T[e][i] for e in B) for i in range(3))
                assert max(profile) <= len(core)
                bases.append({'socle_elements': B, 'residual_profile': profile})
        assert bases
        examples.append({'characteristic': p, 'socle_rank': rank, 'quotient_length': len(core)+rank,
                         'fixed_core_size': len(core), 'socle_columns': columns, 'all_independent_socle_bases': bases})
    return {'downsets_checked': checked, 'proper_chain_points_checked': chains,
            'linear_socle_quotient_models': len(examples), 'linear_socle_commutator_columns': commutators,
            'exact_linear_socle_models': examples}


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    symbolic, family, scales, budget = symbolic_matrices(), family_controls(), asymptotic_controls(), socle_budget_controls()
    sources = ('scripts/run_socle_stability_checks.py', 'docs/initial_cone_stability.md',
               'scripts/run_initial_cone_height_checks.py', 'scripts/run_algebraic_incidence_controls.py',
               'scripts/run_extremal_initial_ideals_checks.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py')
    counts = {**family['counts'], 'symbolic_matrix_entry_identities': symbolic['matrix_product_entry_identities'],
              'all_length_padding_models': len(scales['all_length_padding_checks']),
              'socle_budget_downsets': budget['downsets_checked'],
              'socle_budget_chain_points': budget['proper_chain_points_checked'],
              'arbitrary_linear_socle_models': budget['linear_socle_quotient_models'],
              'linear_socle_commutator_columns': budget['linear_socle_commutator_columns']}
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_socle_stability_checks.py', 'counts': counts,
              'symbolic_matrix_certificates': symbolic, 'all_field_family_controls': family,
              'scale_and_padding_controls': scales, 'general_socle_profile_budget': budget,
              'uniform_relative_margin_hypothesis_status': 'refuted_by_actual_all_field_quotients',
              'sharp_worst_relative_margin': {'dimension_2': 'Theta(n^-2)', 'dimension_3': 'Theta(n^-3)'},
              'unrestricted_G3_bound_changed': False, 'novelty_established': False, 'goal_complete': False,
              'scope': 'Uniform proofs show sharp worst-cone relative stability in two and three variables, with constant-fan actual local quotients over every field. The proposed Omega(1/n) premise and an individual log-area floor Omega(1/n^2) are false. Socle-only quotients have a classical lattice-profile count O(n^(3/2)); hard cone conditioning is not a large-fan lower bound. Exact controls do not establish novelty or the main extremal endpoint.',
              'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}}
    (ROOT/'results/research/socle_stability_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(counts, indent=2))


if __name__ == '__main__':
    main()
