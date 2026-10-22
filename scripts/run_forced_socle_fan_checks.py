#!/usr/bin/env python3
"""Actual local socle quotients, whole fans, and forced primitive-cost jumps."""
import ast
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations
from math import comb, gcd
from pathlib import Path
import hashlib
import json

from run_algebraic_incidence_controls import pair_checks
from run_canonical_fan_deflation_checks import marked_basis
from run_coordinate_slice_budget_checks import edges_from_polygons, profile
from run_ideal_step_incidence_checks import fan_budgets, actual_edge_length
from run_projection_tag_checks import fan, compare_fans
from run_planar_dyadic_bound_checks import convex_hull
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]
TRIANGLE = [tuple(Q(i == j) for i in range(3)) for j in range(3)]


def clean(v, p):
    return {i: c % p if p else c for i, c in v.items() if (c % p if p else c)}


def model(N, kind, p, graded=False):
    """Integral branch reduction followed by the explicit socle quotient."""
    assert N >= 2 and kind in ('H', 'A', 'C')
    ylength = N-1 if kind == 'C' else N
    zlength = N if kind == 'H' else N-1
    reference = tuple([(i, 0, 0) for i in range(N)] +
                      [(i, 1, 0) for i in range(ylength)] +
                      [(i, 0, 1) for i in range(zlength)])
    index = {v: i for i, v in enumerate(reference)}

    @lru_cache(None)
    def ev(a):
        x, y, z = a
        if y and z:
            return {}
        if not y and not z:
            return {index[a]: 1} if x < N else {}
        power = y or z
        start = x+power-1
        if start >= N:
            return {}
        result = defaultdict(int)
        for j in range(N-start):
            coefficient = ((-1)**(power-1)*int(j == 0) if graded else
                           (-1)**(power-1+j)*comb(power+j-2, j)
                           if power > 1 else int(j == 0))
            if not coefficient:
                continue
            exponent = (start+j, int(bool(y)), int(bool(z)))
            if exponent not in index:
                if kind == 'A' and z and start+j == N-1:
                    exponent = (N-1, 1, 0)
                else:
                    continue
            result[index[exponent]] += coefficient
        return clean(result, p)

    assert all(ev(v) == {i: 1} for i, v in enumerate(reference))
    return reference, ev


def apply(M, v, p):
    out = defaultdict(int)
    for j, c in v.items():
        for i, d in M[j].items():
            out[i] += c*d
    return clean(out, p)


def matrices_and_socle(N, kind, p, reference, ev, counts, graded=False):
    n = len(reference)
    M = tuple(tuple(ev(tuple(a+int(i == j) for j, a in enumerate(v)))
                    for v in reference) for i in range(3))
    for i, j in combinations(range(3), 2):
        for k in range(n):
            assert apply(M[i], M[j][k], p) == apply(M[j], M[i][k], p)
            counts['independent_commuting_coordinate_columns'] += 1
    for v in hyperbolic_cross(n, 3):
        for axis in range(3):
            successor = tuple(a+int(j == axis) for j, a in enumerate(v))
            assert apply(M[axis], ev(v), p) == ev(successor)
            counts['actual_quotient_successor_columns'] += 1
    stacked = [{(axis, row): c for axis in range(3) for row, c in M[axis][i].items()}
               for i in range(n)]
    socle = n-len(field_basis(range(n), lambda i: stacked[i], p)[0])
    expected = 3 if kind == 'C' else 2 if kind == 'H' or p == 2 else 1
    assert socle == expected
    counts['stacked_multiplication_socle_dimensions'] += 1
    generators = [
        {(N, 0, 0): 1}, {(0, 1, 1): 1},
        {(0, 2, 0): 1, (1, 1, 0): 1},
        {(0, 0, 2): 1, (1, 0, 1): 1}]
    if not graded:
        generators[2][1, 2, 0] = 1
        generators[3][1, 0, 2] = 1
    if kind == 'A':
        generators.append({(N-1, 1, 0): 1, (N-1, 0, 1): -1})
    elif kind == 'C':
        generators.extend(({(N-1, 1, 0): 1}, {(N-1, 0, 1): 1}))
    for f in generators:
        value = defaultdict(int)
        for a, c in f.items():
            for i, d in ev(a).items():
                value[i] += c*d
        assert not clean(value, p)
        counts['displayed_defining_relations'] += 1
    return socle


def planar_prefix(k, b, branch):
    a = k-b
    selected = {(i, 0, 0) for i in range(a)}
    selected.update(tuple(j if axis == branch else 0 for axis in range(3))
                    for j in range(1, b+1))
    normals = []
    if b:
        q = [a, 0, 0]
        q[branch] = -b
        normals.append(tuple(q))
    if a > 1:
        q = [-(a-1), 0, 0]
        q[branch] = b+1
        normals.append(tuple(q))
    return selected, normals


def predicted(N, kind):
    core = {(i, 0, 0) for i in range(N)}
    choices = ([(N, N, None)] if kind == 'H' else
               [(N-1, N-1, None)] if kind == 'C' else
               [(N, N-1, (0, -1, 1)), (N-1, N, (0, 1, -1))])
    cells = {}
    for yl, zl, side in choices:
        for b in range(yl):
            Y, Ynormals = planar_prefix(yl, b, 1)
            for c in range(zl):
                Z, Znormals = planar_prefix(zl, c, 2)
                poly = list(TRIANGLE)
                for q in Ynormals+Znormals+([side] if side else []):
                    poly = clip(poly, q, 1)
                    if not poly:
                        break
                if not poly or area2(poly) <= 0:
                    continue
                B = frozenset(core | {(a, y+1, z) for a, y, z in Y} |
                              {(a, y, z+1) for a, y, z in Z})
                assert B not in cells
                cells[B] = poly
    assert sum(area2(poly) for poly in cells.values()) == 1
    return cells


def direct_a(n):
    return sum(n//gcd(n, b) for b in range(1, n))


def formulas(N, kind):
    a, b = direct_a(N), direct_a(N-1)
    return ({'g': N*N, 'W': 2*N*a, 'L': 2*N*N*(N-1)} if kind == 'H' else
            {'g': (N-1)*(N+2), 'W': N*a+(N+1)*b+4*N-4,
             'L': 2*(N-1)*(N*N+1)} if kind == 'A' else
            {'g': (N-1)**2, 'W': 2*(N-1)*b,
             'L': 2*(N-1)**2*(N-2)})


def actual_controls(counts):
    rows = []
    for N in (2, 3, 4, 5, 7):
        for p in (0, 2, 3, 101):
            cases = {}
            for kind in ('H', 'A', 'C'):
                reference, ev = model(N, kind, p)
                n = len(reference)
                socle = matrices_and_socle(N, kind, p, reference, ev, counts)
                expected = predicted(N, kind)
                direct = fan(n, ev, p, counts)
                compare_fans(direct, expected, counts)
                _, walls, _, W = fan_budgets(direct)
                L = actual_edge_length(direct)
                assert {'g': len(direct), 'W': W, 'L': L} == formulas(N, kind)
                counts['complete_fans_vs_uniform_count_and_cost_formulas'] += 1
                for B, poly in direct.items():
                    G, _ = marked_basis(B, ev, p)
                    w = tuple(sum(v[j] for v in poly)/len(poly) for j in range(3))
                    counts['ordinary_reduced_basis_S_pairs'] += pair_checks(G, w, p)
                cases[kind] = {'length': n, 'socle_dimension': socle,
                               'g': len(direct), 'W': W, 'L': L,
                               'whole_cones': ([{'standard_monomials': sorted(B),
                                                'rational_cone_vertices': [list(map(str, v)) for v in poly]}
                                               for B, poly in expected.items()] if p == 0 else []),
                               'whole_cones_reference': {'N': N, 'characteristic': 0, 'kind': kind}}
            parent_ref, parent_ev = model(N, 'A', p)
            child_ref, child_ev = model(N, 'C', p)
            # The actual quotient sends precisely the common top socle to
            # zero; identify coordinates by their monomial reference labels.
            child_index = {a: i for i, a in enumerate(child_ref)}
            for a in hyperbolic_cross(3*N-1, 3):
                projected = {child_index[parent_ref[i]]: c for i, c in parent_ev(a).items()
                             if parent_ref[i] in child_index}
                assert projected == child_ev(a)
                counts['independent_actual_socle_quotient_columns'] += 1
            A, C = cases['A'], cases['C']
            assert A['g']-C['g'] == 3*N-3
            assert A['L']-C['L'] == 6*N*N-8*N+2
            rows.append({'N': N, 'characteristic': p,
                         'one_dimensional_ideal_quotient_unique': p != 2,
                         'cases': cases,
                         'parent_minus_child': {key: A[key]-C[key] for key in ('g', 'W', 'L')}})
        print(json.dumps({'N': N, 'lengths': [3*N-1, 3*N-2],
                          'parent': formulas(N, 'A'), 'child': formulas(N, 'C')}), flush=True)
    return rows


def graded_prediction(N):
    patches = defaultdict(list)
    for order in permutations(range(3)):
        position = {a: i for i, a in enumerate(order)}
        B = {(i, 0, 0) for i in range(N)}
        for d in range(1, N):
            B.add((d-1, 1, 0) if position[0] < position[1] else (0, d, 0))
            B.add((d-1, 0, 1) if position[0] < position[2] else (0, 0, d))
        branch = 1 if position[1] < position[2] else 2
        top = [0, 0, 0]
        if position[0] < position[branch]:
            top[0], top[branch] = N-1, 1
        else:
            top[branch] = N
        B.add(tuple(top))
        assert len(B) == 3*N-1
        poly = list(TRIANGLE)
        for a, b in zip(order, order[1:]):
            q = tuple(int(j == b)-int(j == a) for j in range(3))
            poly = clip(poly, q, 1)
        assert area2(poly) > 0
        patches[frozenset(B)].extend(poly)
    return {B: [(*v, 1-sum(v)) for v in convex_hull([(a[0], a[1]) for a in vertices])]
            for B, vertices in patches.items()}


def graded_controls(counts):
    rows = []
    for N in (2, 3, 4, 5, 7):
        for p in (0, 2, 3, 101):
            reference, ev = model(N, 'A', p, graded=True)
            other_ref, other = model(N, 'A', p)
            assert reference == other_ref
            n = len(reference)
            socle = matrices_and_socle(N, 'A', p, reference, ev, counts, graded=True)
            T = []
            for j, a in enumerate(reference):
                power = a[1]+a[2]
                image = defaultdict(int)
                for k in range(power+1):
                    for i, c in other((a[0]+k, a[1], a[2])).items():
                        image[i] += comb(power, k)*c
                image = clean(image, p)
                assert image[j] == 1
                assert all(i == j or reference[i][0] > a[0] for i in image)
                T.append(image)
            for a in hyperbolic_cross(n, 3):
                direct = defaultdict(int)
                power = a[1]+a[2]
                for k in range(power+1):
                    for i, c in other((a[0]+k, a[1], a[2])).items():
                        direct[i] += comb(power, k)*c
                assert apply(T, ev(a), p) == clean(direct, p)
                counts['independent_unit_triangular_isomorphism_columns'] += 1
            expected = graded_prediction(N)
            direct = fan(n, ev, p, counts)
            compare_fans(direct, expected, counts)
            _, _, _, W = fan_budgets(direct)
            L = actual_edge_length(direct)
            assert len(direct) == (4 if N == 2 else 6)
            assert W == (8 if N == 2 else 12)
            assert L == 4*N*N-6*N+6
            for B, poly in direct.items():
                G, _ = marked_basis(B, ev, p)
                w = tuple(sum(v[j] for v in poly)/len(poly) for j in range(3))
                counts['graded_presentation_reduced_basis_S_pairs'] += pair_checks(G, w, p)
            rows.append({'N': N, 'characteristic': p, 'length': n,
                         'socle_dimension': socle, 'graded_initial_cones': len(direct),
                         'unit_rescaled_initial_cones': formulas(N, 'A')['g'],
                         'graded_W': W, 'graded_L': L,
                         'isomorphism': ['x -> x', 'Y -> (1+x)y', 'Z -> (1+x)z'],
                         'isomorphism_matrix_determinant': 1,
                         'complete_graded_cones': ([{'standard_monomials': sorted(B),
                                                    'rational_cone_vertices': [list(map(str, v)) for v in poly]}
                                                   for B, poly in expected.items()] if p == 0 else []),
                         'complete_graded_cones_reference': {'N': N, 'characteristic': 0}})
    return rows


def arithmetic_controls(counts):
    bound = 4000
    phi = list(range(bound+1))
    primes = []
    for p in range(2, bound+1):
        if phi[p] == p:
            primes.append(p)
            for j in range(p, bound+1, p):
                phi[j] -= phi[j]//p
    sums = [0]*(bound+1)
    for d in range(1, bound+1):
        for j in range(d, bound+1, d):
            sums[j] += d*phi[d]
    for n in range(1, 301):
        assert sums[n]-1 == direct_a(n)
        counts['independent_totient_sieve_vs_gcd_wall_sums'] += 1
    rows = []
    for p in primes:
        if p == 2 or p+1 > bound:
            continue
        for N in (p, p+1):
            d = N*(sums[N]-1)-(N-3)*(sums[N-1]-1)+4*N-4
            if N == p:
                assert 4*d >= p**3
            elif p >= 47:
                assert 8*d <= -p**3
            counts['independent_prime_forced_jump_inequalities'] += 1
        if p in (3, 5, 7, 47, 101, 997, 1999):
            for N in (p, p+1):
                a, b = formulas(N, 'A'), formulas(N, 'C')
                rows.append({'prime': p, 'N': N, 'parent_length': 3*N-1,
                             'g_parent': a['g'], 'g_child': b['g'],
                             'W_parent': a['W'], 'W_child': b['W'],
                             'W_difference': a['W']-b['W'],
                             'L_difference': 6*N*N-8*N+2,
                             'evidence': 'Exact uniform formulas, not full fan enumeration at this size.'})
    return rows


def hashes():
    paths = ['scripts/run_forced_socle_fan_checks.py', 'docs/forced_socle_fan_steps.md',
             'docs/adjacent_rank_fans.md', 'docs/ideal_step_incidence.md',
             'docs/wall_cost_products.md']
    pending = [p for p in paths if p.endswith('.py')]
    while pending:
        path = pending.pop()
        for node in ast.walk(ast.parse((ROOT/path).read_text())):
            names = ([node.module] if isinstance(node, ast.ImportFrom) and node.module
                     else [a.name for a in node.names] if isinstance(node, ast.Import) else [])
            for name in names:
                candidate = 'scripts/'+name.replace('.', '/')+'.py'
                if (ROOT/candidate).is_file() and candidate not in paths:
                    paths.append(candidate)
                    pending.append(candidate)
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def main():
    if not __debug__:
        raise RuntimeError('Assertions must be enabled.')
    before = hashes()
    counts = Counter()
    cases = actual_controls(counts)
    graded = graded_controls(counts)
    arithmetic = arithmetic_controls(counts)
    assert before == hashes(), 'Source changed during verification.'
    report = {'command': 'python3 -B scripts/run_forced_socle_fan_checks.py',
              'status': 'passed', 'source_sha256': before, 'counts': dict(counts),
              'actual_local_algebras': cases, 'exact_arithmetic_witnesses': arithmetic,
              'same_algebra_coordinate_separation': graded,
              'primary_sources': [{'url': 'https://arxiv.org/pdf/1802.06491v2',
                                   'scope': 'Theorem 3.1 at dimension zero: one-dimensional socle is the classical Artinian local Gorenstein criterion.'},
                                  {'url': 'https://arxiv.org/pdf/math/0207135',
                                   'scope': 'Theorem 2.4: the classical state-polyhedron dictionary; no priority claim for the explicit family follows.'}],
              'uniform_statement': 'For char(K)!=2 the local Gorenstein A_N has one unique ideal quotient C_N of length one less. At odd primes N=p the primitive cost difference is >=p^3/4; at N=p+1, p>=47 it is <=-p^3/8. All fan and cost formulas hold over every field, with a two-dimensional parent socle in characteristic two.',
              'limits': ['Complete proofs in docs/forced_socle_fan_steps.md establish the uniform statements; finite checks only support their implementation.',
                         'No choice of immediate quotient can repair a quadratic-polylog primitive-cost step bound. A cubic total budget along a specially chosen complete chain is not refuted.',
                         'Actual edge length has only the quadratic difference 6*N^2-8*N+2 in this family.',
                         'Both fans have quadratic size; no supergeneric lower bound, improved unrestricted G3 exponent or required complexity consequence follows.',
                         'Local socle criteria, planar inverse-curve fans and ideal deletion are established tools. Exact priority of the assembled forced-step family is unassessed.']}
    (ROOT/'results/research/forced_socle_fan_checks.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
