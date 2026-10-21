#!/usr/bin/env python3
"""Actual rank-indistinguishable square-zero quotients and their whole fans."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
from math import gcd
from pathlib import Path
import hashlib
import json

from run_coordinate_slice_budget_checks import (
    downset, edges_from_polygons, primitive, profile,
)
from run_interleaving_research import minimal_excluded_border
from run_primitive_wall_checks import threshold
from run_structural_complexity_checks import area2, clip, dot

ROOT = Path(__file__).resolve().parents[1]


def diagrams(k):
    m = 2*k
    z = [2*m*(m-i)+(m-i)**2 for i in range(m)]
    U = tuple((2*i, z[i]) for i in range(m))
    V = tuple((2*i+1, z[i]) for i in range(m))
    D0 = frozenset((a, b) for y, zz in V
                   for a in range(y+1) for b in range(zz+1))
    D1 = frozenset((a, b) for y, zz in U
                   for a in range(y+1) for b in range(zz+1))
    D = frozenset((0, *a) for a in D0) | frozenset((1, *a) for a in D1)
    assert D1 <= D0 and downset(D)
    assert len(D1) == (8*m**3+7*m-3)//3
    assert len(D0) == (8*m**3+9*m*m+7*m)//3
    return U, V, D0, D1, D


def model(k, perm, counts):
    U, V, D0, D1, D = diagrams(k)
    m = len(U)
    pairs = tuple(((1, *u), (0, *V[perm[i]])) for i, u in enumerate(U))
    tips = {p for pair in pairs for p in pair}
    assert len(tips) == 2*m
    assert all(tuple(a+(i == j) for j, a in enumerate(p)) not in D
               for p in tips for i in range(3))
    reference = sorted(D-{u for u, _ in pairs})
    index = {p: i for i, p in enumerate(reference)}
    replacements = dict(pairs)

    def ev(alpha):
        return index.get(replacements.get(alpha, alpha))

    n = len(reference)
    assert n == len(D)-m == (128*k**3+36*k*k+22*k-3)//3
    M = [tuple(ev(tuple(a+(i == j) for j, a in enumerate(p))) for p in reference)
         for i in range(3)]

    def apply(i, v):
        return None if v is None else M[i][v]

    for i, j in combinations(range(3), 2):
        for c in range(n):
            assert apply(i, M[j][c]) == apply(j, M[i][c])
            counts['independent_commuting_unit_matrix_columns'] += 1
    for p in D:
        for i in range(3):
            successor = tuple(a+(i == j) for j, a in enumerate(p))
            assert apply(i, ev(p)) == ev(successor)
            counts['independent_monomial_successor_replays'] += 1
    for u, v in pairs:
        assert ev(u) == ev(v) is not None
        assert all(apply(i, ev(u)) is None for i in range(3))
        counts['independent_socle_generator_identifications'] += 1
    for p in minimal_excluded_border(tuple(D)):
        assert ev(p) is None
        counts['independent_monomial_kernel_border_replays'] += 1
    lower = [ev((0, *p)) for p in D0]
    upper = [ev((1, *p)) for p in D1]
    assert None not in lower and len(set(lower)) == len(D0)
    assert None not in upper and len(set(upper)) == len(D1)
    for p in minimal_excluded_border(tuple(D0)):
        assert ev((0, *p)) is None
    for p in minimal_excluded_border(tuple(D1)):
        assert ev((1, *p)) is None
    counts['independent_complete_planar_kernel_and_rank_formulas'] += 2
    nonzero_x = {v for v in M[0] if v is not None}
    assert len(nonzero_x) == len(D1)
    assert all(apply(0, v) is None for v in M[0])
    assert n-2*len(D1) == 3*m*m-m+1
    counts['independent_square_zero_Jordan_types'] += 1
    mod_x_diagram = D0-frozenset(V)
    assert all((ev((0, *p)) not in nonzero_x) == (p in mod_x_diagram) for p in D0)
    assert len(mod_x_diagram) == n-len(D1) > len(D1)
    counts['independent_mod_x_planar_quotient_kernels'] += 1
    return n, D, pairs, ev, {
        'lower_planar_length': len(D0), 'upper_planar_length': len(D1),
        'lower_monomial_kernel_generators': minimal_excluded_border(tuple(D0)),
        'upper_monomial_kernel_generators': minimal_excluded_border(tuple(D1)),
        'mod_x_planar_kernel_generators': minimal_excluded_border(tuple(mod_x_diagram)),
        'mod_x_planar_length': len(mod_x_diagram),
        'x_Jordan_blocks': {'size_two': len(D1), 'size_one': n-2*len(D1)},
        'permutation': perm,
    }


def arrangement(normals):
    cells = [[tuple(Q(i == j) for i in range(3)) for j in range(3)]]
    for q in sorted(set(normals)):
        new = []
        for poly in cells:
            for sign in (-1, 1):
                piece = clip(poly, q, sign)
                if area2(piece) > 0:
                    new.append(piece)
        cells = new
    assert sum(map(area2, cells)) == 1
    return cells


def actual_fan(k, perm, kind, counts):
    n, D, pairs, ev, metadata = model(k, perm, counts)
    qs = [primitive(tuple(a-b for a, b in zip(u, v))) for u, v in pairs]
    geometric = arrangement(qs)
    fixed = D-{p for pair in pairs for p in pair}
    actual, ledger = {}, []
    for poly in geometric:
        w = tuple(sum(v[i] for v in poly)/len(poly) for i in range(3))
        assert min(w) > 0 and all(dot(w, q) for q in qs)
        choices = [u if dot(w, tuple(a-b for a, b in zip(u, v))) < 0 else v
                   for u, v in pairs]
        B = fixed | frozenset(choices)
        assert len(B) == n and downset(B)
        columns = {ev(a): a for a in B}
        assert None not in columns and len(columns) == n
        independent = []
        seen = set()
        for a in sorted(D, key=lambda a: (dot(w, a), a)):
            c = ev(a)
            if c not in seen:
                assert c is not None
                seen.add(c)
                independent.append(a)
        assert frozenset(independent) == B
        counts['independent_full_ground_greedy_standard_sets'] += 1
        normals = set()
        for head in minimal_excluded_border(tuple(B)):
            c = ev(head)
            if c is not None:
                tail = columns[c]
                q = tuple(a-b for a, b in zip(head, tail))
                assert dot(w, q) > 0
                d = gcd(*q)
                normals.add(tuple(a//d for a in q))
            counts['independent_actual_border_normal_forms'] += 1
        check = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for q in sorted(normals):
            check = clip(check, q, 1)
        assert set(check) == set(poly)
        actual[B] = check
        counts['independent_complete_actual_cone_equalities'] += 1
        ledger.append({'chosen_relation_tips': choices, 'profile': profile(B),
                       'positive_cone': [[str(a) for a in v] for v in check],
                       'standard_set_sha256': hashlib.sha256(
                           json.dumps(sorted(B), separators=(',', ':')).encode()).hexdigest()})
    assert len(actual) == len(geometric)
    edges = edges_from_polygons(actual, counts)
    directions = Counter()
    costs = Counter()
    weighted_jumps = []
    for pair in edges:
        B, C = tuple(pair)
        delta = tuple(a-b for a, b in zip(profile(B), profile(C)))
        q = primitive(delta)
        directions[q] += 1
        costs[threshold(q)] += 1
        weighted_jumps.append(gcd(*delta))
    value = sum(c*p for c, p in costs.items())
    if kind == 'identity':
        assert len(actual) == 2 and len(edges) == 1 and value == 2
        assert weighted_jumps == [2*k]
    else:
        assert len(actual) == (k+1)**2 and len(edges) == 2*k*(k+1)
        assert value == 2*k*(k+1)*(9*k+5)
        assert set(weighted_jumps) == {1}
        assert all(q[0] and q[1]*q[2] < 0 for q in directions)
        assert all(p == k+1 for p in directions.values())
        for i, j in product(range(k), repeat=2):
            Li, Lj = 8*k-4*i-1, 8*k-4*j-1
            w = (Q(3*Li-Lj, 4), Q(Li+Lj, 4), Q(1))
            assert min(w) > 0
            assert dot((1, 1, -Li), w) == dot((1, -3, Lj), w) == 0
            counts['independent_positive_remaining_wall_crossings'] += 1
    counts['complete_actual_quotient_fans'] += 1
    out = {'k': k, 'kind': kind, 'length': n, **metadata,
           'initial_ideals': len(actual), 'undirected_walls': len(edges),
           'complete_area': '1', 'primitive_wall_cost': value,
           'directions': [{'direction': q, 'pieces': p, 'cost': threshold(q)}
                          for q, p in sorted(directions.items())],
           'zero_component_cost': sum(threshold(q)*p for q, p in directions.items() if not q[0]),
           'singleton_x_sign_cost': sum(threshold(q)*p for q, p in directions.items()
                                       if q[0] and q[1] <= 0 and q[2] <= 0),
           'cells': ledger}
    return out


def main():
    counts, pairs = Counter(), []
    for k in (1, 2, 3, 4):
        identity = list(range(2*k))
        swapped = [i ^ 1 for i in range(2*k)]
        A = actual_fan(k, identity, 'identity', counts)
        B = actual_fan(k, swapped, 'swapped', counts)
        for key in ('lower_planar_length', 'upper_planar_length',
                    'lower_monomial_kernel_generators', 'upper_monomial_kernel_generators',
                    'mod_x_planar_length', 'mod_x_planar_kernel_generators',
                    'x_Jordan_blocks'):
            assert A[key] == B[key]
        assert B['zero_component_cost'] == B['singleton_x_sign_cost'] == 0
        counts['independent_equal_marginal_rank_and_Jordan_pairs'] += 1
        pairs.append({'identity': A, 'swapped': B})
        print('complete rank-indistinguishable pair:', k, A['length'],
              A['initial_ideals'], B['initial_ideals'], flush=True)
    sources = ['docs/planar_layer_rank_obstruction.md',
               'scripts/run_planar_layer_obstruction_checks.py',
               'docs/primitive_wall_spectrum.md', 'docs/rational_shadow_bound.md',
               'docs/square_zero_extension_fans.md',
               'scripts/run_coordinate_slice_budget_checks.py',
               'scripts/run_interleaving_research.py',
               'scripts/run_primitive_wall_checks.py',
               'scripts/run_structural_complexity_checks.py']
    out = {'schema_version': 1, 'all_checks_passed': True,
           'scope': 'Exact whole fans and actual commuting unit-coordinate quotients '
                    'support a uniform all-field rank-indistinguishability theorem. '
                    'The family is binomial and has a linear-volume monic carrier box; '
                    'no supergeneric lower family, unrestricted bound, full mixed-wall '
                    'budget or novelty is established.',
           'uniform_theorem': {
               'length': 'n_k=(128*k^3+36*k^2+22*k-3)/3',
               'same_data': 'Identical labelled planar column matroids, elimination and '
                            'colon ideals, planar quotient modulo x, and square-zero x-Jordan type.',
               'fan_counts': '2 versus (k+1)^2',
               'wall_costs': '2 versus 2*k*(k+1)*(9*k+5)',
               'missing_sign_scope': 'Every swapped-family wall has nonzero x component '
                                     'and opposite nonzero y,z signs; both controlled '
                                     'coordinate budgets vanish.',
               'information_consequence': 'These marginal data cannot determine the fan '
                                          'or cost. Symmetric approximation requires '
                                          'Omega(n_k^(1/3)) for the count and '
                                          'Omega(n_k^(1/2)) for cost.'},
           'counters': dict(counts), 'actual_quotient_pairs': pairs,
           'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}}
    target = ROOT/'results/research/planar_layer_obstruction_checks.json'
    target.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'all_checks_passed': True, 'counters': dict(counts),
                      'report': str(target.relative_to(ROOT))}), flush=True)


if __name__ == '__main__':
    main()
