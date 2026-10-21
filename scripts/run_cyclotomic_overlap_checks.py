#!/usr/bin/env python3
"""Exact overlap, coherent-integrality and power-substitution controls.

No full count for distinct coprime periods or unrestricted G3 follows.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import floor, gcd, prod
from pathlib import Path
import hashlib
import json

from run_cyclotomic_price_checks import (
    border_certificate, cyclotomic_model, dot, frontier_basis,
    generic_tower, insert, solver, tower_basis,
)
from run_interleaving_research import minimal_excluded_border
from run_structural_complexity_checks import area2, clip, determinant, finite_staircases

ROOT = Path(__file__).resolve().parents[1]


def strict_weight(n, template):
    base = n+1
    return tuple(base**4*a+base**i for i, a in enumerate(template))


def restricted_greedy(model, ground, weight):
    pivots, basis = {}, []
    for a in sorted(ground, key=lambda a: (dot(a, weight), a)):
        if insert(pivots, model['residue'](a), model['field']):
            basis.append(a)
    assert len(basis) == model['n']
    return frozenset(basis)


def overlap_control(model, axis, D, R, template, counts):
    n = D*R
    w = strict_weight(n, template)
    other = [i for i in range(3) if i != axis]
    alpha = w[axis]
    m = min(w[i] for i in other)
    expensive = max(other, key=lambda i: w[i])
    delta = abs(w[other[0]]-w[other[1]])
    rho = Q(w[other[1]]-w[other[0]], alpha)
    C = tower_basis(model, axis, D, R, rho)
    low_ground = []
    for a in range(D):
        for b in range(R):
            for c in range(R-b):
                v = [0]*3
                v[axis], v[other[0]], v[other[1]] = a, b, c
                low_ground.append(tuple(v))
    assert restricted_greedy(model, low_ground, w) == C
    B = frontier_basis(model['residue'], n, w, model['field'])
    H = frozenset(a for a in B if sum(a[i] for i in other) >= R)
    assert B-H <= C and len(C-B) == len(H)
    degrees = sorted(sum(a[i] for i in other) for a in B)
    assert all(d >= i//D for i, d in enumerate(degrees))
    J0 = n*(R-1)//2
    E = sum(degrees)-J0
    mx = sum(a[axis] for a in B)
    mz = sum(a[expensive] for a in B)
    assert E >= 0 and 2*(alpha*mx+m*E+delta*mz) <= alpha*n*(D-1)
    kappa = floor(Q((D-1)*alpha, m))
    assert max(degrees) <= R-1+kappa
    M = len(H)
    q, s = divmod(M, D)
    required = D*q*(q+1)//2+s*(q+1)
    assert E >= required and M*M <= 2*D*E
    assert m*M*M <= D*n*(D-1)*alpha
    cert = border_certificate(B, model['residue'], w, model['field'])
    counts['independent_restricted_greedy_equalities'] += 1
    counts['outside_column_deletion_controls'] += 1
    counts['exact_excess_budget_controls'] += 1
    counts['strict_border_relations'] += len(cert)
    counts['positive_overlap_cases'] += bool(H)
    return {'periods': model.get('periods'), 'general_tower_relation': model.get('relation'),
            'characteristic': model['field'], 'axis': axis, 'D': D, 'R': R, 'length': n,
            'weight': w, 'rho': str(rho), 'baseline_tag_moment': J0,
            'excess_tag_moment': E, 'expensive_tag_moment': mz,
            'pure_axis_moment': mx, 'overlap_layers_bound': kappa,
            'high_degree_columns': sorted(H), 'deleted_low_basis_columns': sorted(C-B),
            'high_column_count': M, 'minimum_excess_for_this_count': required,
            'standard_set': sorted(B), 'low_basis': sorted(C)}


def all_overlap_controls(counts):
    rows = []
    templates = ((1, 2, 2), (1, 5, 5), (1, 8, 9), (2, 1, 3), (5, 1, 1), (3, 2, 5))
    for periods in ((2, 3, 5), (3, 5, 7)):
        r, s, t = periods
        dims = (((r-1)*(s-1), t-1), ((r-1)*(t-1), s-1), ((s-1)*(t-1), r-1))
        for p in (0, 2):
            model = cyclotomic_model(periods, p)
            for axis, (D, R) in enumerate(dims):
                for w0 in templates:
                    rows.append(overlap_control(model, axis, D, R, w0, counts))
    for p, coupled in product((0, 2, 3), (False, True)):
        model = generic_tower(p, coupled)
        for w0 in templates:
            rows.append(overlap_control(model, 0, model['D'], model['R'], w0, counts))
    assert counts['positive_overlap_cases'] > 0
    return rows


def small_complete_fans(counts):
    candidates = finite_staircases(8, 3)
    assert len(candidates) == 160
    integer = cyclotomic_model((3, 3, 3), 0)['residue']
    data, rational = [], None
    for p in (0, 2, 3, 5):
        residue = cyclotomic_model((3, 3, 3), p)['residue']
        tally, cells = Counter(), []
        for state in candidates:
            B = tuple(sorted(state))
            d = int(determinant(list(zip(*(integer(a) for a in B)))))
            if not d or (p and d % p == 0):
                continue
            tally['independent_downsets'] += 1
            normal_form, normals, relations = solver([residue(a) for a in B], p), set(), []
            for head in minimal_excluded_border(B):
                tail = normal_form(residue(head))
                assert all(sum(c*residue(B[j])[i] for j, c in tail.items()) % p == residue(head)[i] % p
                           if p else sum(c*residue(B[j])[i] for j, c in tail.items()) == residue(head)[i]
                           for i in range(8))
                for j, coefficient in tail.items():
                    if coefficient:
                        z = tuple(x-y for x, y in zip(head, B[j]))
                        g = gcd(*z)
                        normals.add(tuple(x//g for x in z))
                relations.append({'head': head, 'tail': [
                    {'exponent': B[j], 'coefficient': str(c)} for j, c in sorted(tail.items())]})
            polygon = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
            for z in sorted(normals):
                if not polygon:
                    break
                polygon = clip(polygon, z, 1)
            if not polygon:
                continue
            w = tuple(sum(a[i] for a in polygon)/len(polygon) for i in range(3))
            assert all(dot(w, z) > 0 for z in normals)
            assert frontier_basis(residue, 8, w, p) == state
            tally['coherent_downsets'] += 1
            tally['coherent_integer_index_'+str(abs(d))] += 1
            cells.append({'standard_set': B, 'integer_determinant': d,
                          'strict_weight': list(map(str, w)),
                          'polygon': [[str(x) for x in a] for a in polygon],
                          'normalized_area': str(area2(polygon)),
                          'reduced_border_relations': relations})
        assert sum(Q(c['normalized_area']) for c in cells) == 1
        assert len(cells) == (6 if p == 2 else 9)
        if not p:
            rational = cells
            assert tally['coherent_integer_index_2'] == 3
        counts['complete_small_tensor_fans'] += 1
        counts['exhausted_downsets'] += len(candidates)
        counts['coherent_small_tensor_states'] += len(cells)
        data.append({'periods': [3, 3, 3], 'characteristic': p, 'length': 8,
                     'scope': 'Equal periods; this is not the distinct-prime cyclotomic field.',
                     'counts': dict(tally), 'cones': cells})
    w = strict_weight(8, (41, 32, 32))
    B = frontier_basis(integer, 8, w, 0)
    assert abs(determinant(list(zip(*(integer(a) for a in sorted(B)))))) == 2
    residue2 = cyclotomic_model((3, 3, 3), 2)['residue']
    B2 = frontier_basis(residue2, 8, w, 2)
    assert B != B2
    witness = {'weight': w, 'Q_standard_set': sorted(B), 'F2_standard_set': sorted(B2),
               'Q_border_relations': border_certificate(B, integer, w, 0),
               'F2_border_relations': border_certificate(B2, residue2, w, 2)}
    return data, rational, witness


def lifted_power(k, exponent):
    """Ordinary monic division by T^(2k)+T^k+1, without period lookup."""
    values = [0]*max(2*k, exponent+1)
    values[exponent] = 1
    for j in range(exponent, 2*k-1, -1):
        c = values[j]
        values[j] = 0
        values[j-k] -= c
        values[j-2*k] -= c
    return {i: x for i, x in enumerate(values[:2*k]) if x}


def lifted_residue(k, exponent):
    a, b, c = exponent
    powers = [lifted_power(k, e) for e in (a+b, a+c, b+c)]
    return {ij: prod(P[i] for P, i in zip(powers, ij))
            for ij in product(*(tuple(P) for P in powers))}


def power_lift_controls(cells, witness, counts):
    rows = []
    for k in (1, 3, 5):
        assert gcd(k, 6) in (1, 3)  # Odd k makes the pair-product grading invertible.
        grades = set()
        for offset in product(range(k), repeat=3):
            a, b, c = offset
            grade = ((a+b) % k, (a+c) % k, (b+c) % k)
            assert grade not in grades
            grades.add(grade)
        assert len(grades) == k**3
        for cell in cells:
            B = cell['standard_set']
            expanded = frozenset(tuple(k*x+y for x, y in zip(a, offset))
                                 for a in B for offset in product(range(k), repeat=3))
            assert len(expanded) == 8*k**3
            assert set(minimal_excluded_border(tuple(expanded))) == {
                tuple(k*x for x in relation['head']) for relation in cell['reduced_border_relations']}
            for relation in cell['reduced_border_relations']:
                terms = [(relation['head'], Q(1))]+[(t['exponent'], -Q(t['coefficient']))
                                                           for t in relation['tail']]
                reduced = Counter()
                for exponent, coefficient in terms:
                    for row, value in lifted_residue(k, tuple(k*x for x in exponent)).items():
                        reduced[row] += coefficient*value
                assert not any(reduced.values())
                counts['lifted_whole_polynomial_annihilations'] += 1
            absolute, rank2, blocks = 1, 0, []
            for offset in product(range(k), repeat=3):
                columns, row_labels = [], set()
                for a in B:
                    v = tuple(k*x+y for x, y in zip(a, offset))
                    sparse = lifted_residue(k, v)
                    columns.append(sparse)
                    row_labels.update(sparse)
                assert len(row_labels) == 8
                matrix = [[column.get(row, 0) for column in columns] for row in sorted(row_labels)]
                d = int(determinant(matrix))
                assert abs(d) == abs(cell['integer_determinant'])
                absolute *= abs(d)
                pivots = {}
                for column in zip(*matrix):
                    insert(pivots, column, 2)
                assert len(pivots) == (7 if abs(d) == 2 else 8)
                rank2 += len(pivots)
                blocks.append({'offset': offset, 'determinant': d})
                counts['independent_lifted_integer_blocks'] += 1
            expected = 2**(k**3) if abs(cell['integer_determinant']) == 2 else 1
            assert absolute == expected
            assert rank2 == (7*k**3 if expected > 1 else 8*k**3)
            rows.append({'power': k, 'length': 8*k**3,
                         'base_standard_set': B, 'base_integer_index': abs(cell['integer_determinant']),
                         'lifted_integer_index': str(absolute), 'rank_in_characteristic_two': rank2,
                         'integer_blocks': blocks})
        for relation in witness['F2_border_relations']:
            terms = [(relation['head'], 1)]+[(t['exponent'], -int(t['coefficient']))
                                                        for t in relation['tail']]
            reduced = Counter()
            for exponent, coefficient in terms:
                for row, value in lifted_residue(k, tuple(k*x for x in exponent)).items():
                    reduced[row] += coefficient*value
            assert all(x % 2 == 0 for x in reduced.values())
            counts['lifted_modular_polynomial_annihilations'] += 1
    return rows


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    counts = Counter()
    overlaps = all_overlap_controls(counts)
    fans, rational, witness = small_complete_fans(counts)
    lifts = power_lift_controls(rational, witness, counts)
    sources = ('docs/cyclotomic_overlap.md', 'scripts/run_cyclotomic_overlap_checks.py',
               'docs/cyclotomic_price_regions.md', 'scripts/run_cyclotomic_price_checks.py',
               'scripts/run_structural_complexity_checks.py', 'scripts/run_interleaving_research.py')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_cyclotomic_overlap_checks.py',
              'counts': dict(counts), 'all_weight_overlap_controls': overlaps,
              'complete_equal_period_fans': fans, 'same_weight_characteristic_change': witness,
              'independent_power_lift_blocks': lifts,
              'uniform_results': {
                  'low_basis_deletion': 'B=(C_rho minus T) union H, |T|=|H|, H has tag degree >=R.',
                  'moment_budget': 'alpha*mu_x + m*(mu_tag-n*(R-1)/2) + abs(beta-gamma)*mu_exp <= alpha*n*(D-1)/2.',
                  'high_degree_count': 'If M=|H|=D*q+s, E>=D*q*(q+1)/2+s*(q+1); M^2<=2*D*E<=D*n*(D-1)*alpha/m.',
                  'degree_cutoff': 'tag degree <=R-1+floor((D-1)*alpha/m).',
                  'nilpotent_tag_corollary': 'For g(T)=T^R, the whole positive fan, at all weights, has at most 2*(D-1)*(R-1)+2 states for D,R>=2; C_rho is its exact basis.',
                  'cheap_cyclotomic_regions': 'Coherent bases are integral with index one; pointwise regional initials are the same in every characteristic, by the established two-circuit unimodularity theorem.',
                  'power_lift': 'For every odd k, I_k=I_1(X^k,Y^k,Z^k), length 8*k^3; fan counts 9 outside characteristic two and 6 in characteristic two. Three rational cones have coherent integer index 2^(k^3).',
                  'distinct_prime_full_fan_classified': False,
                  'unrestricted_G3_improved': False, 'novelty_established': False, 'goal_complete': False},
              'scope': 'Uniform overlap and arithmetic proofs with bounded exact controls. Equal-period coherent torsion and its power lifts do not establish an amplification or characteristic change for distinct coprime prime periods. Classical matroid sensitivity, lucky-prime and monomial-substitution principles are not claimed as new.',
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/'results/research/cyclotomic_overlap_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(counts), indent=2))


if __name__ == '__main__':
    main()
