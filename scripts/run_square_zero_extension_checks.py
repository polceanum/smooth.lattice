#!/usr/bin/env python3
"""Exact square-zero extension, whole-fan and tangent-weight scope controls."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import gcd
from pathlib import Path
from random import Random
import hashlib
import json

from run_interleaving_research import minimal_excluded_border
from run_projection_tag_checks import compare_fans, fan
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import clip, determinant, dot

ROOT = Path(__file__).resolve().parents[1]
TRIANGLE = tuple(tuple(Q(i == j) for i in range(3)) for j in range(3))


def divides(a, b):
    return all(x <= y for x, y in zip(a, b))


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def clean(v, p):
    return {i: c % p if p else c for i, c in v.items() if (c % p if p else c)}


def nullspace(rows, width, p):
    A = [list(v) for v in rows if any(v)]
    pivots, k = [], 0
    for j in range(width):
        i = next((i for i in range(k, len(A)) if A[i][j] % p), None)
        if i is None:
            continue
        A[k], A[i] = A[i], A[k]
        inv = pow(A[k][j] % p, -1, p)
        A[k] = [c*inv % p for c in A[k]]
        for i in range(len(A)):
            if i != k and A[i][j] % p:
                c = A[i][j]
                A[i] = [(x-c*y) % p for x, y in zip(A[i], A[k])]
        pivots.append(j)
        k += 1
    basis = []
    for free in sorted(set(range(width))-set(pivots)):
        v = [int(j == free) for j in range(width)]
        for row, j in zip(A, pivots):
            v[j] = -row[free] % p
        assert all(sum(a*b for a, b in zip(row, v)) % p == 0 for row in rows)
        basis.append(v)
    return basis


def hom_space(D, E, p):
    """All pair syzygies, reduced in the monomial target B/K."""
    D, E = tuple(sorted(D)), tuple(sorted(E))
    G = tuple(minimal_excluded_border(D))
    variables = tuple(product(range(len(G)), E))
    index = {v: i for i, v in enumerate(variables)}
    equations = []
    for i, j in combinations(range(len(G)), 2):
        lcm = tuple(max(a, b) for a, b in zip(G[i], G[j]))
        rows = defaultdict(lambda: [0]*len(variables))
        for k, sign in ((i, 1), (j, -1)):
            shift = tuple(a-b for a, b in zip(lcm, G[k]))
            for v in E:
                target = add(v, shift)
                if target in E:
                    rows[target][index[k, v]] += sign
        equations.extend(rows.values())
    return G, variables, nullspace(equations, len(variables), p)


def extension(D, E, G, tails, p):
    """Monomial nested endpoints, including nonflat square-zero quotients."""
    D, E = tuple(sorted(D)), tuple(sorted(E))
    assert set(E) <= set(D)
    Kgens = tuple(minimal_excluded_border(E)) if E else ((0, 0),)
    product_gens = tuple(add(g, h) for g, h in product(G, Kgens))
    xmax = max(g[0] for g in product_gens)
    zmax = max(g[1] for g in product_gens)
    outer = tuple(v for v in product(range(xmax), range(zmax))
                  if not any(divides(g, v) for g in product_gens))
    lower = {v: i for i, v in enumerate(D)}
    upper = {v: len(D)+i for i, v in enumerate(E)}
    @lru_cache(None)
    def ev(a):
        x, y, z = a
        v = y, z
        if x >= 2:
            return {}
        if x:
            return {upper[v]: 1} if v in upper else {}
        if v in lower:
            return {lower[v]: 1}
        choices = []
        for g, f in zip(G, tails):
            if divides(g, v):
                shift = tuple(a-b for a, b in zip(v, g))
                image = {upper[add(w, shift)]: c for w, c in f.items()
                         if add(w, shift) in upper}
                choices.append(clean(image, p))
        assert choices and all(v == choices[0] for v in choices)
        return choices[0]
    ground = tuple(sorted({(0, *v) for v in outer} | {(1, *v) for v in E}))
    return len(D)+len(E), ev, ground


def matrix_controls(n, ev, basis, ground, p, counts):
    matrices = [[ev(tuple(a+(i == j) for j, a in enumerate(v)))
                 for v in basis] for i in range(3)]
    def apply(M, v):
        out = defaultdict(int)
        for j, a in v.items():
            for k, b in M[j].items():
                out[k] += a*b
        return clean(out, p)
    for i, j in combinations(range(3), 2):
        for k in range(n):
            assert apply(matrices[i], matrices[j][k]) == apply(matrices[j], matrices[i][k])
            counts['independent_commuting_matrix_columns'] += 1
    for k in range(n):
        assert not apply(matrices[0], matrices[0][k])
    for v in ground:
        for i in range(3):
            assert apply(matrices[i], ev(v)) == ev(tuple(a+(i == j) for j, a in enumerate(v)))
            counts['independent_module_successor_columns'] += 1


def minor_value(A, p):
    if not p:
        return determinant(A)
    A = [list(row) for row in A]
    det = 1
    for j in range(len(A)):
        i = next((i for i in range(j, len(A)) if A[i][j] % p), None)
        if i is None:
            return 0
        if i != j:
            A[i], A[j] = A[j], A[i]
            det = -det
        c = A[j][j] % p
        det = det*c % p
        inv = pow(c, -1, p)
        for i in range(j+1, len(A)):
            factor = A[i][j]*inv % p
            A[i] = [(a-factor*b) % p for a, b in zip(A[i], A[j])]
    return det % p


def arbitrary_compatible_models(counts):
    rectangle = lambda a, b: set(product(range(a), range(b)))
    hook = {(0, 0), (1, 0), (0, 1)}
    cases = [(rectangle(2, 2), rectangle(2, 2)), (hook, hook),
             (rectangle(3, 2), rectangle(2, 1)), (hook, {(0, 0)}),
             (rectangle(2, 2), set())]
    fixtures = []
    for case, (D, E) in enumerate(cases):
        for p in (2, 3, 101):
            G, variables, hom = hom_space(D, E, p)
            rng = Random(1703+case*101+p)
            parameters = [rng.randrange(p) for _ in hom]
            coefficients = [sum(c*v[j] for c, v in zip(parameters, hom)) % p
                            for j in range(len(variables))]
            tails = [dict() for _ in G]
            for (i, v), c in zip(variables, coefficients):
                if c:
                    tails[i][v] = c
            n, ev, ground = extension(D, E, G, tails, p)
            basis = tuple([(0, *v) for v in sorted(D)]+[(1, *v) for v in sorted(E)])
            matrix_controls(n, ev, basis, ground, p, counts)
            fixed = frozenset((0, *v) for v in D)
            variable = tuple(v for v in ground if v not in fixed)
            E = set(E)
            for selected in combinations(variable, len(E)):
                U = [v for v in selected if not v[0]]
                T = {(v[1], v[2]) for v in selected if v[0]}
                S = sorted(E-T)
                A = [[ev(u).get(len(D)+sorted(E).index(v), 0) for u in U] for v in S]
                assert len(S) == len(U)
                determinant_nonzero = bool(minor_value(A, p))
                direct = len(field_basis(tuple(fixed)+selected, ev, p)[0]) == n
                assert direct == determinant_nonzero
                counts['independent_full_basis_vs_minor_equalities'] += 1
            cells = fan(n, ev, p, counts)
            assert all(fixed <= B for B in cells)
            fixtures.append({'case': case, 'characteristic': p, 'length': n,
                             'endpoint_lengths': [len(D), len(E)], 'hom_dimension': len(hom),
                             'x_rank': len(E), 'flat': D == E, 'initial_ideals': len(cells)})
    return fixtures


def tangent_weights(D):
    return {(a, -b-1) for i, j in D
            for a in [max(u-i for u, v in D if v == j)]
            for b in [max(v-j for u, v in D if u == i)]} | {
            (-a-1, b) for i, j in D
            for a in [max(u-i for u, v in D if v == j)]
            for b in [max(v-j for u, v in D if u == i)]}


def three_cell_family(counts):
    fixtures = []
    for a, t in ((2, 2), (2, 3), (3, 2), (3, 4), (2, 7)):
        D = set(product(range(a), range(t)))
        G = ((a, 0), (0, t))
        tails = ({}, {(a-1, 0): 1, (a-1, 1): 1})
        q1, q2 = (1, a-1, -(t-1)), (1, a-1, -(2*t-1))
        assert q2[1:] not in tangent_weights(D)
        for p in (0, 2, 3, 101):
            n, ev, ground = extension(D, D, G, tails, p)
            F = {(0, *v) for v in D} | {(1, i, j) for i, j in D if i < a-1}
            A = {(1, a-1, j) for j in range(t)}
            B = {(0, 0, t+j) for j in range(t)}
            standards = [frozenset(F | A),
                         frozenset(F | {(1, a-1, 0)} | (B-{(0, 0, 2*t-1)})),
                         frozenset(F | B)]
            constraints = [(tuple(-v for v in q1),),
                           (q1, tuple(-v for v in q2)), (q2,)]
            expected = {}
            for S, normals in zip(standards, constraints):
                assert len(S) == n
                poly = list(TRIANGLE)
                for q in normals:
                    poly = clip(poly, q, 1)
                expected[S] = poly
            actual = fan(n, ev, p, counts)
            compare_fans(actual, expected, counts)
            assert len(actual) == 3
            summed = defaultdict(int)
            for j in range(t):
                for i, c in ev((0, 0, t+j)).items():
                    summed[i] += (-1)**j*c
            assert clean(summed, p) == ev((1, a-1, 0))
            row_ids = [next(iter(ev((1, a-1, j)))) for j in range(t)]
            C = [[ev((0, 0, t+j)).get(i, 0) for j in range(t)] for i in row_ids]
            assert minor_value(C, p) == 1
            circuit = [(1, a-1, 0)]+[(0, 0, t+j) for j in range(t)]
            for omitted in range(t+1):
                remaining = circuit[:omitted]+circuit[omitted+1:]
                assert len(field_basis(remaining, ev, p)[0]) == t
                counts['independent_long_circuit_minimality_checks'] += 1
            matrix_controls(n, ev, tuple(sorted({(0, *v) for v in D} |
                                                {(1, *v) for v in D})), ground, p, counts)
            counts['all_field_circuit_formula_replays'] += 1
            fixtures.append({'a': a, 't': t, 'characteristic': p, 'length': n,
                             'whole_fan_cells': 3, 'actual_wall_directions': [q1, q2],
                             'long_circuit_size': t+1,
                             'standard_sets': [sorted(S) for S in standards]})
    return fixtures


def saved_hook_guards(counts):
    path = 'results/research/coordinate_slice_budget_checks.json'
    saved = json.loads((ROOT/path).read_text())
    for source, digest in saved['source_sha256'].items():
        assert hashlib.sha256((ROOT/source).read_bytes()).hexdigest() == digest
    rows = []
    for fixture in saved['complete_coupled_hook_fans']:
        t = fixture['thickness']
        R = t*t+1
        D = {v for v in product(range(R), repeat=2) if v[0] < t or v[1] < t}
        weights = tangent_weights(D)
        _, _, hom = hom_space(D, D, 2)
        assert len(hom) == 2*len(D) == fixture['length']
        directions = set()
        for wall in fixture['walls']:
            q = tuple(wall['direction'])
            difference = tuple(a-b for a, b in zip(*wall['endpoints']))
            assert gcd(*q) == 1
            assert difference == tuple(wall['edge_multiple']*v for v in q) or (
                   difference == tuple(-wall['edge_multiple']*v for v in q))
            if q[0]:
                assert q[0] == 1
                directions.add(q)
            counts['saved_exact_wall_profile_replays'] += 1
        outside = sorted(q for q in directions if q[1:] not in weights)
        assert outside
        if t == 3:
            assert len(directions) == 153 > len(hom) == 102
        rows.append({'length': fixture['length'], 'hom_dimension': len(hom),
                     'distinct_tangent_weights': len(weights),
                     'distinct_nonzero_x_wall_directions': len(directions),
                     'wall_directions_outside_tangent_weights': outside,
                     'scope': 'Replayed saved exact complete fan; no asymptotic inference.'})
    return rows


def main():
    counts = defaultdict(int)
    compatible = arbitrary_compatible_models(counts)
    family = three_cell_family(counts)
    guards = saved_hook_guards(counts)
    sources = ['docs/square_zero_extension_fans.md', 'scripts/run_square_zero_extension_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_projection_tag_checks.py',
               'scripts/run_signed_face_checks.py', 'scripts/run_structural_complexity_checks.py',
               'results/research/coordinate_slice_budget_checks.json']
    report = {'command': 'python3 -B scripts/run_square_zero_extension_checks.py',
              'counts': dict(sorted(counts.items())), 'compatible_nested_extensions': compatible,
              'all_field_three_cell_family': family, 'saved_hook_tangent_scope_guards': guards,
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              'unrestricted_G3_improved': False, 'collective_cubic_budget_proved': False,
              'polynomial_lower_separation_proved': False, 'novelty_established': False,
              'goal_complete': False}
    (ROOT/'results/research/square_zero_extension_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'counts': report['counts'], 'saved_guards': guards}, indent=2))


if __name__ == '__main__':
    main()
