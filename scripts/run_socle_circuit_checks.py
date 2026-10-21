#!/usr/bin/env python3
"""Independent quotient fans for disjoint socle circuits and repeated pairs."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
from math import gcd
from pathlib import Path
import hashlib
import json

from run_interleaving_research import minimal_excluded_border
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_primitive_wall_checks import divisors
from run_projection_tag_checks import compare_fans
from run_quadratic_coordinate_fan_checks import TRIANGLE
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import area2, clip, cross, dot

ROOT = Path(__file__).resolve().parents[1]


def scalar(a, char):
    return a % char if char else Q(a)


def model(D, polynomials, char):
    S = set().union(*(set(f) for f in polynomials))
    assert sum(map(len, polynomials)) == len(S)
    assert all(len(f) >= 2 and all(scalar(c, char) for c in f.values()) for f in polynomials)
    assert S <= D
    assert all(tuple(c+(i == j) for j, c in enumerate(a)) not in D for a in S for i in range(3))
    # Pick a different pivot per disjoint circuit. The remaining columns are units.
    heads = [max(f) for f in polynomials]
    basis = tuple(sorted(D-set(heads)))
    index = {a: i for i, a in enumerate(basis)}
    columns = {a: {i: scalar(1, char)} for a, i in index.items()}
    for head, f in zip(heads, polynomials):
        inv = pow(scalar(f[head], char), -1, char) if char else 1/Q(f[head])
        columns[head] = {index[a]: scalar(-c*inv, char) for a, c in f.items() if a != head}
    variable = tuple(sorted(S))
    block_basis = tuple(a for a in basis if a in S)
    block_index = {index[a]: i for i, a in enumerate(block_basis)}
    block = {a: {block_index[i]: c for i, c in columns[a].items()} for a in S}
    return {'D': D, 'S': S, 'fixed': D-S, 'basis': basis, 'polynomials': polynomials,
            'n': len(basis), 'char': char, 'ev': lambda a: columns.get(a, {}),
            'block_ev': lambda a: block.get(a, {}), 'variable': variable}


def multiplication_controls(A, counts):
    n, ev, char = A['n'], A['ev'], A['char']
    matrices = [[ev(tuple(c+(i == j) for j, c in enumerate(a))) for a in A['basis']]
                for i in range(3)]
    def apply(M, vector):
        out = defaultdict(int)
        for j, c in vector.items():
            for k, d in M[j].items():
                out[k] += c*d
        return {k: scalar(c, char) for k, c in out.items() if scalar(c, char)}
    for i, j in combinations(range(3), 2):
        for k in range(n):
            assert apply(matrices[i], matrices[j][k]) == apply(matrices[j], matrices[i][k])
            counts['independent_commuting_matrix_columns'] += 1
    for a in A['D']:
        for i in range(3):
            assert apply(matrices[i], ev(a)) == ev(tuple(c+(i == j) for j, c in enumerate(a)))
            counts['independent_monomial_successor_replays'] += 1


def normalized(f, head, char):
    c = scalar(f[head], char)
    inv = pow(c, -1, char) if char else 1/c
    return tuple(sorted((a, scalar(b*inv, char)) for a, b in f.items()))


def exact_fan(A, counts):
    """Residue-column greedy and actual full borders, with coloops contracted."""
    n, ev, char = A['n'], A['block_ev'], A['char']
    queue, cells, edges, polynomials = [((1, n+1, (n+1)**2), None)], {}, set(), set()
    expected_nonmonomial = {normalized(f, a, char) for f in A['polynomials'] for a in f}
    old_border = set(minimal_excluded_border(tuple(A['D'])))
    while queue:
        w, previous = queue.pop()
        variable, pivots = field_basis(sorted(A['variable'], key=lambda a: (dot(w, a), a)), ev, char)
        assert len(variable) == len(A['S'])-len(A['polynomials'])
        B = frozenset(A['fixed']|set(variable))
        if previous is not None:
            assert B != previous
            edges.add(frozenset((B, previous)))
        if B in cells:
            continue
        assert len(B) == n
        normals = set()
        for head in minimal_excluded_border(tuple(B)):
            rem, co = field_reduce(ev(head), pivots, char)
            assert not rem
            f = {head: scalar(1, char)}
            f.update({a: scalar(-c, char) for a, c in co.items()})
            key = tuple(sorted(f.items()))
            if len(f) == 1:
                assert head in old_border
            else:
                assert key in expected_nonmonomial
            polynomials.add(key)
            replay = defaultdict(int)
            for a, c in f.items():
                for k, d in ev(a).items():
                    replay[k] += c*d
            assert all(not scalar(c, char) for c in replay.values())
            for a in co:
                q = tuple(x-y for x, y in zip(head, a))
                divisor = gcd(*q)
                q = tuple(c//divisor for c in q)
                assert dot(q, w) > 0
                normals.add(q)
            counts['independent_full_border_normal_forms'] += 1
        poly = list(TRIANGLE)
        for q in normals:
            poly = clip(poly, q, 1)
        assert poly and area2(poly) > 0
        cells[B] = poly
        for a, b in zip(poly, poly[1:]+poly[:1]):
            if any(a[i] == b[i] == 0 for i in range(3)):
                continue
            q = next(q for q in normals if dot(q, a) == dot(q, b) == 0)
            neighbor, _ = exact_adjacent_weight(n, a, b, q)
            queue.append((neighbor, B))
    assert sum(area2(poly) for poly in cells.values()) == 1
    M = len(A['S'])
    assert len(polynomials) <= len(old_border)+M
    assert sum(len(f) for f in polynomials) <= len(old_border)+sum(len(f)**2 for f in A['polynomials'])
    counts['complete_actual_quotient_fans'] += 1
    counts['exact_universal_reduced_basis_union_controls'] += 1
    return cells, len(polynomials), sum(map(len, polynomials))


def factor_fan(A):
    cells = [(list(TRIANGLE), [])]
    for f in A['polynomials']:
        updated = []
        for poly, heads in cells:
            for head in f:
                patch = list(poly)
                for a in f:
                    if a != head:
                        patch = clip(patch, tuple(x-y for x, y in zip(head, a)), 1)
                if patch:
                    updated.append((patch, heads+[head]))
        cells = updated
    return {frozenset(A['D']-set(heads)): poly for poly, heads in cells}


def width_controls(A, counts):
    N, M = len(A['D']), len(A['S'])
    r = 1
    while r**3 < N:
        r += 1
    for antichain, outside in ((A['S'], False), (set(minimal_excluded_border(tuple(A['D']))), True)):
        small = [a for a in antichain if max(a) < r]
        assert len({a[1:] for a in small}) == len(small) <= r*r
        multiplicity = Counter()
        for a in antichain:
            if max(a) < r:
                continue
            i = next(i for i, c in enumerate(a) if c >= r)
            for k in range(1, r+1) if outside else range(r):
                b = tuple(c-k*(i == j) for j, c in enumerate(a))
                assert b in A['D']
                multiplicity[b] += 1
        assert max(multiplicity.values(), default=0) <= 3
        assert len(antichain)*r <= r**3+3*N
        counts['antichain_disjoint_ray_budget_controls'] += 1
    assert N <= 2*A['n']
    if max(a[0] for a in A['D']) < 2:
        assert M*M <= 4*N
        counts['two_layer_socle_width_controls'] += 1


def recross_model(s, char):
    B = 4*s*s
    tips = [(0, 2*i+int(i % 2 == 0), B-i*i+int(i % 2 == 1)) for i in range(1, 2*s+1)]
    upper = ((1, 1, B), (1, 0, B+1))
    D = set().union(*(divisors(a) for a in tips+list(upper)))
    f = [{upper[0]: 1}, {upper[1]: 1}]
    for i, a in enumerate(tips, 1):
        f[i % 2][a] = -1
    return model(D, f, char), tips, upper


def two_pencil_model(k, char):
    degree = 4*k
    D1 = {(i, j) for i in range(degree) for j in range(degree-i)}
    tips, fs = [], []
    for i in range(2*k):
        delta = int(i >= k)
        lower = (0, 2*i, degree+delta-2*i)
        upper = (1, 3*k-1-i, k+i)
        tips.extend((lower, upper))
        fs.append({lower: 1, upper: -1})
    D = {(x, *a) for x in (0, 1) for a in D1}
    D |= set().union(*(divisors(a) for a in tips))
    return model(D, fs, char), tips


def grid_model(L, char):
    fs, tips = [], []
    for a, b in product(range(L+1, 2*L+1), repeat=2):
        v, w = (2*a, 2*b, 9*L-2*a-2*b), (a, b, 9*L-a-b+1)
        tips.extend((v, w)); fs.append({v: 1, w: -1})
    D = set().union(*(divisors(a) for a in tips))
    return model(D, fs, char), tips


def positive_grid_arrangement(L):
    directions = [(a, b, -a-b-1) for a, b in product(range(L+1, 2*L+1), repeat=2)]
    vertices = defaultdict(set)
    for i, j in combinations(range(len(directions)), 2):
        v = cross(directions[i], directions[j])
        if not v[0]*v[1]*v[2] or not (min(v) > 0 or max(v) < 0):
            continue
        divisor = gcd(*v)
        v = tuple(c//divisor for c in v)
        if min(v) < 0:
            v = tuple(-c for c in v)
        vertices[v].update((i, j))
    regions = 1+len(directions)+sum(len(v)-1 for v in vertices.values())
    ordinary = sum(len(v) == 2 for v in vertices.values())
    certified = 0
    differences = []
    for a, b in product(range(L//2+1, 3*L//4+1), repeat=2):
        if gcd(a, b) == 1:
            certified += (L-a)*(L-b)
            differences.append((a, b))
    assert certified <= ordinary
    return {'L': L, 'lines': len(directions), 'positive_vertices': len(vertices),
            'ordinary_vertices': ordinary, 'complete_positive_regions': regions,
            'certified_primitive_long_difference_vertices': certified,
            'primitive_long_differences': differences}


def run_fixture(A, counts):
    multiplication_controls(A, counts)
    width_controls(A, counts)
    actual, union_size, terms = exact_fan(A, counts)
    predicted = factor_fan(A)
    compare_fans(actual, predicted, counts)
    M = len(A['S'])
    assert len(actual) <= 1+3*M+(3*M)*(3*M-1)//2
    witnesses = []
    for B, poly in actual.items():
        w = tuple(sum(v[i] for v in poly)/len(poly) for i in range(3))
        witnesses.append({'deleted_heads': sorted(A['D']-B), 'weight': list(map(str, w))})
    return actual, {'length': A['n'], 'monomial_length': len(A['D']),
                    'socle_support_size': M, 'independent_relations': len(A['polynomials']),
                    'initial_ideals': len(actual), 'union_reduced_basis_size': union_size,
                    'union_reduced_basis_terms': terms, 'whole_fan_witnesses': witnesses}


def main():
    counts, recross, pencils, grids = defaultdict(int), [], [], []
    for s, char in product(range(1, 5), (0, 2, 3)):
        A, tips, upper = recross_model(s, char)
        actual, record = run_fixture(A, counts)
        assert A['n'] == (32*s**3+30*s*s+10*s+12)//3
        assert len(actual) == 4*s
        profiles = {B: tuple(sum(a[i] for a in B) for i in range(3)) for B in actual}
        crossings = []
        for i in range(1, 2*s):
            u, t = Q(2*i+1, 2), Q(i*i+i)
            w = (t, u, Q(1))
            best = min(dot(w, v) for v in profiles.values())
            incident = [B for B, v in profiles.items() if dot(w, v) == best]
            assert len(incident) == 4
            assert {tuple(a in B for a in upper) for B in incident} == set(product((False, True), repeat=2))
            assert set.intersection(*(set(B) for B in incident)) == A['D']-set(upper)-set(tips[i-1:i+1])
            crossings.append({'u': str(u), 't': str(t), 'active_lower_tips': tips[i-1:i+1]})
            counts['same_upper_pair_distinct_persistent_crossing_controls'] += 1
        recross.append(dict(record, parameter=s, characteristic=char, upper_pair=upper,
                            lower_tips=tips, repeated_crossings=crossings))
    for k, char in product(range(1, 5), (0, 2)):
        A, tips = two_pencil_model(k, char)
        _, record = run_fixture(A, counts)
        assert A['n'] == 16*k*k+6*k
        assert record['initial_ideals'] == k*k+4*k-1
        pencils.append(dict(record, parameter=k, characteristic=char, paired_tips=tips))
    for L, char in product(range(1, 4), (0, 2)):
        A, tips = grid_model(L, char)
        _, record = run_fixture(A, counts)
        arrangement = positive_grid_arrangement(L)
        assert record['initial_ideals'] == arrangement['complete_positive_regions']
        grids.append(dict(record, parameter=L, characteristic=char, paired_tips=tips))
    sources = ['docs/disjoint_socle_circuits.md', 'scripts/run_socle_circuit_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_oracle_enumeration_checks.py',
               'scripts/run_projection_tag_checks.py', 'scripts/run_signed_face_checks.py',
               'scripts/run_primitive_wall_checks.py', 'scripts/run_structural_complexity_checks.py']
    report = {'command': 'python3 -B scripts/run_socle_circuit_checks.py',
              'counts': dict(sorted(counts.items())), 'same_pair_recrossing_family': recross,
              'sharp_two_layer_families': pencils, 'sharp_spatial_families': grids,
              'grid_arrangement_checks': [positive_grid_arrangement(L) for L in range(1, 17)],
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              'unrestricted_G3_improved': False, 'quadratic_coordinate_endpoint_settled': False,
              'polynomial_lower_separation_proved': False, 'novelty_established': False,
              'goal_complete': False}
    (ROOT/'results/research/socle_circuit_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'counts': report['counts'], 'recross_fans': len(recross),
                      'two_layer_fans': len(pencils), 'spatial_fans': len(grids)}, indent=2))


if __name__ == '__main__':
    main()
