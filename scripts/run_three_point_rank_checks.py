#!/usr/bin/env python3
"""Integral binary ideal controls for the open higher-seed rank identity."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import permutations, product
from math import comb, gcd
from pathlib import Path
from random import Random
import hashlib
import json

from run_higher_seed_checks import ceil, column_rank, core_rank, regular_blocks, smith_diagonal
from run_higher_interaction_checks import linear_sum, smith_units
from run_structural_complexity_checks import determinant

ROOT = Path(__file__).resolve().parents[1]


def expand_product(exponents):
    """Repeated polynomial multiplication, independent of the binomial formula."""
    x, y, q = exponents
    out = {(x, y): 1}
    for _ in range(q):
        nxt = defaultdict(int)
        for (a, b), c in out.items():
            nxt[a+1, b] += c
            nxt[a, b+1] += c
        out = dict(nxt)
    return out


def ideal_columns(generators, D):
    out = []
    for exponents in generators:
        degree = sum(exponents)
        if degree > D:
            continue
        polynomial = expand_product(exponents)
        for h in range(D-degree+1):
            out.append({a+h: c for (a, b), c in polynomial.items()})
    return out


def hilbert_value(generators, D, p):
    return D+1-column_rank(ideal_columns(generators, D), p)


def component_ideal(L, levels, negative=False):
    levels = [(j, lo, hi) for j, lo, hi in levels if 0 <= lo <= hi]
    M = max((hi for _, _, hi in levels), default=0) if negative else 0
    transformed = [(j, M-hi, M-lo) for j, lo, hi in levels] if negative else levels
    D = max([L-1]+[j+hi for j, lo, hi in transformed])
    generators = [(L, 0, 0), (0, D+1, 0)]
    generators += [(j, D-j-hi, lo) for j, lo, hi in transformed]
    assert all(min(v) >= 0 for v in generators)
    return D, generators


def source_columns(L, levels, negative=False):
    out = []
    for j, lo, hi in levels:
        for i in range(lo, hi+1):
            if j >= L:
                out.append({})
            elif not negative:
                out.append({j+h: comb(i, h) for h in range(min(i+1, L-j))})
            elif not i:
                out.append({j: 1})
            else:
                out.append({j+h: (-1)**h*comb(i+h-1, h) for h in range(L-j)})
    return out


def cokernel_signature(nrows, columns):
    if not columns:
        return {'free_dimension': nrows, 'nonunit_invariants': []}
    M = [[col.get(i, 0) for col in columns] for i in range(nrows)]
    units, residual = smith_units(M)
    diagonal = smith_diagonal(residual)
    return {'free_dimension': nrows-units-len(diagonal),
            'nonunit_invariants': [x for x in diagonal if x > 1]}


def arbitrary_interval_controls():
    counts, fixtures = defaultdict(int), []
    rng = Random(170310)
    for case in range(48):
        L = 1+case % 7
        levels = []
        for j in range(L+1):
            lo = rng.randrange(7)
            hi = rng.randrange(9)
            if lo <= hi:
                levels.append((j, lo, hi))
        for negative in (False, True):
            raw = source_columns(L, levels, negative)
            D, generators = component_ideal(L, levels, negative)
            homogeneous = ideal_columns(generators, D)
            signature = cokernel_signature(L, raw)
            assert signature == cokernel_signature(D+1, homogeneous)
            counts['independent_integral_component_signatures'] += 1
            for p in (0, 2, 3, 5, 7):
                assert column_rank(raw, p)+hilbert_value(generators, D, p) == L
                counts['independent_component_field_ranks'] += 1
            if case < 8:
                for perm in permutations(range(3)):
                    permuted = [tuple(v[i] for i in perm) for v in generators]
                    assert cokernel_signature(D+1, ideal_columns(permuted, D)) == signature
                    counts['integral_three_form_permutation_signatures'] += 1
            if case < 4:
                fixtures.append({'L': L, 'levels': levels, 'negative': negative,
                                 'degree': D, 'generators_x_y_xplusy': generators,
                                 'cokernel': signature})
    return dict(counts), fixtures


def actual_blocks(s, a, r, u, Z, T, negative):
    A = a*(Z-u)
    if not negative:
        B = a*u+(r-a)*Z
        for k, b0, L in regular_blocks(s, a, r):
            levels = [(j-b0, max(0, ceil((k*Z+B*j-T)/A)), (k+(r-a)*j)//a)
                      for j in range(b0, b0+L)]
            yield {'signed_degree': k, 'starting_V_exponent': b0}, L, levels
    else:
        for b in range(s//a+1):
            for k in range(min(r, s-a*b+1)):
                L = (s-a*b-k)//r+1
                levels = [(j, max(0, ceil((a*u*b+k*Z+r*Z*j-T)/A)), (k+r*j)//a)
                          for j in range(L)]
                yield {'free_V_exponent': b, 'z_residue': k}, L, levels


def coherent_controls():
    counts, fixtures = defaultdict(int), []
    for s, a, r in ((8, 2, 3), (11, 2, 3), (14, 2, 4), (12, 3, 5)):
        for u, T in product((F(2, 7), F(5, 7)), (F(2), F(s, 2), F(s))):
            for negative in (False, True):
                blocks = list(actual_blocks(s, a, r, u, F(1), T, negative))
                for p in (0, 2, 3, 5):
                    total_rank, total_hilbert, ambient = 0, 0, 0
                    for label, L, levels in blocks:
                        D, generators = component_ideal(L, levels, negative)
                        h = hilbert_value(generators, D, p)
                        total_hilbert += h
                        total_rank += L-h
                        ambient += L
                        counts['coherent_block_hilbert_components'] += 1
                    direct = core_rank(s, a, r, u, F(1), T, p, negative)
                    assert total_rank == direct
                    assert total_rank+total_hilbert == ambient
                    counts['independent_core_hilbert_sum_queries'] += 1
                    if p == 0 and T == F(s, 2):
                        fixtures.append({'cutoff': s, 'degree': a, 'interaction': r,
                                         'binary_price': str(u), 'z_price': '1',
                                         'price_cutoff': str(T), 'canonical': negative,
                                         'ambient_dimension': ambient,
                                         'rank': total_rank, 'component_dimension_sum': total_hilbert})
    return dict(counts), fixtures


def shifted(poly, x, y):
    return {(a+x, b+y): c for (a, b), c in poly.items()}


def arithmetic_block_controls():
    counts = defaultdict(int)
    block = next((label, L, levels) for label, L, levels in
                 actual_blocks(11, 2, 3, F(1), F(7), F(14), False)
                 if label['signed_degree'] == 2)
    label, L, levels = block
    raw = source_columns(L, levels)
    assert raw == [{0: 1}, {0: 1, 1: 1}, {1: 1, 2: 1}, {2: 1, 3: 2}]
    D, generators = component_ideal(L, levels)
    assert D == 4
    polynomials = [expand_product(v) for v in generators]
    x4, y5, y3, g1, g2 = polynomials
    x2y2 = linear_sum(((1, g1), (-1, shifted(y3, 1, 0))))
    two_x3y = linear_sum(((1, g2), (-1, x4), (-1, x2y2)))
    assert x2y2 == {(2, 2): 1} and two_x3y == {(3, 1): 2}
    for poly in polynomials:
        for (x, y), c in poly.items():
            assert x >= 4 or y >= 3 or (x >= 2 and y >= 2) or (x >= 3 and y >= 1 and c % 2 == 0)
    signature = cokernel_signature(D+1, ideal_columns(generators, D))
    assert signature == {'free_dimension': 0, 'nonunit_invariants': [2]}
    assert signature == cokernel_signature(L, raw)
    minimal = ((4, 0), (0, 3), (2, 2), (3, 1))
    assert all(not (x <= X and y <= Y) for i, (x, y) in enumerate(minimal)
               for j, (X, Y) in enumerate(minimal) if i != j)
    counts['coherent_integral_ideal_equalities'] += 1
    counts['coherent_cyclic_torsion_components'] += 1
    triple = [(2, 0, 0), (0, 2, 0), (0, 0, 2)]
    fields = []
    for p in (0, 2, 3, 5, 7):
        least = next(d for d in (2, 3)
                     if len(ideal_columns(triple, d)) > column_rank(ideal_columns(triple, d), p))
        gap = 6-2*least
        length = sum(hilbert_value(triple, d, p) for d in range(4))
        assert 4*length == 12+gap**2
        expected = 1 if p == 2 else 0
        assert hilbert_value(generators, D, p) == expected
        assert hilbert_value(triple, 2, p) == expected
        assert cokernel_signature(3, ideal_columns(triple, 2)) == signature
        fields.append({'characteristic': p, 'triple_least_syzygy_degree': least,
                       'triple_syzygy_gap': gap, 'triple_colength': length,
                       'block_component_dimension': expected})
        counts['independent_three_form_syzygy_gap_controls'] += 1
    fixture = {'cutoff': 11, 'degree': 2, 'interaction': 3,
               'binary_price': '1', 'z_price': '7', 'price_cutoff': '14',
               'block': label, 'L': L, 'selected_levels': levels,
               'binary_degree': D, 'generators_x_y_xplusy': generators,
               'integer_monomial_generators': [[4, 0, 1], [0, 3, 1], [2, 2, 1], [3, 1, 2]],
               'component_cokernel': signature, 'field_controls': fields,
               'whole_algebra_equivalence_with_three_form_cell': False}
    return dict(counts), fixture


def support_counterexample():
    label, L, levels = next((label, L, levels) for label, L, levels in
                           actual_blocks(17, 2, 3, F(9, 17), F(1), F(17, 3), True)
                           if label['free_V_exponent'] == 1 and label['z_residue'] == 2)
    cols = source_columns(L, levels, True)
    assert L == len(cols) == 5
    relation = (1, -1, -1, 0, 0)
    assert not linear_sum(zip(relation, cols))
    matched_rows = (0, 4, 1, 3, 2)
    assert len(set(matched_rows)) == 5
    assert all(cols[j].get(i, 0) for j, i in enumerate(matched_rows))
    selected = (0, 2, 3, 4)
    minors = [determinant([[cols[j].get(i, 0) for j in selected] for i in I])
              for I in ((0, 1, 2, 3), (0, 1, 2, 4))]
    assert minors == [F(2), F(-7)] and gcd(*(int(x) for x in minors)) == 1
    for p in (0, 2, 3, 5, 7, 11):
        assert column_rank(cols, p) == 4
    signature = cokernel_signature(L, cols)
    assert signature == {'free_dimension': 1, 'nonunit_invariants': []}
    D, generators = component_ideal(L, levels, True)
    assert signature == cokernel_signature(D+1, ideal_columns(generators, D))
    fixture = {'cutoff': 17, 'degree': 2, 'interaction': 3,
               'binary_price': '9/17', 'z_price': '1', 'price_cutoff': '17/3',
               'block': label, 'L': L, 'selected_levels': levels,
               'columns_up_to_unit_signs': [[[i, c] for i, c in sorted(v.items())] for v in cols],
               'integer_kernel_vector': relation, 'support_matching_rows': matched_rows,
               'four_column_indices': selected, 'four_minors': list(map(int, minors)),
               'rank_over_every_field': 4, 'support_matching_size': 5,
               'component_cokernel': signature}
    return {'coherent_all_field_support_rank_counterexamples': 1}, fixture


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    arbitrary, interval_fixtures = arbitrary_interval_controls()
    coherent, coherent_fixtures = coherent_controls()
    arithmetic, arithmetic_fixture = arithmetic_block_controls()
    support, support_fixture = support_counterexample()
    sources = ('docs/three_point_rank_ideals.md', 'scripts/run_three_point_rank_checks.py',
               'docs/higher_seed_binomial_reduction.md', 'scripts/run_higher_seed_checks.py',
               'scripts/run_higher_interaction_checks.py',
               'scripts/run_hypersurface_truncation_checks.py',
               'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_socle_stability_checks.py')
    report = {
        'status': 'passed', 'command': 'python3 -B scripts/run_three_point_rank_checks.py',
        'counts': {**arbitrary, **coherent, **arithmetic, **support},
        'interval_component_fixtures': interval_fixtures,
        'coherent_hilbert_sum_fixtures': coherent_fixtures,
        'coherent_arithmetic_block': arithmetic_fixture,
        'coherent_support_counterexample': support_fixture,
        'uniform_results': {
            'integral_component_reduction': 'Each selected block cokernel is exactly a homogeneous component of an explicit Artinian ideal generated by products of x,y,x+y, preserving the free rank and nonunit Smith invariant factors; unit pivot counts may differ.',
            'core_identity_equivalence': '(B8) is equivalent to the equality of the two sums of Hilbert-function components in (P4).',
            'three_form_symmetry': 'Permutations of x,y,x+y preserve component cokernels via integral linear changes, up to unit signs.',
            'support_matching_shortcut': False,
            'coherent_arithmetic_example': 'A genuine regular block is a degree-four Z/2 component; a classical three-form cell shares the local torsion but has a different whole algebra.',
            'one_triple_theorem_covers_all_blocks': False,
            'p_fractal_or_Frobenius_application_established': False},
        'uniform_core_identity_proved': False,
        'r_ge_3_full_signed_transfer_proved': False,
        'general_r_full_family_linear_bound_proved': False,
        'unrestricted_G3_improvement': False, 'novelty_established': False, 'goal_complete': False,
        'scope': 'A classical-tool supporting reduction with exact arithmetic and a coherent all-field obstruction to support matching. No uniform Hilbert-function identity, full original count, new p-fractal theorem or priority certificate.',
        'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/'results/research/three_point_rank_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
