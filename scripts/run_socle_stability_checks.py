#!/usr/bin/env python3
"""Exact stability constructions and shared profile budgets for propagating relations."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, gcd, lcm, prod
from pathlib import Path
from random import Random
import hashlib
import json

from run_interleaving_research import minimal_excluded_border
from run_initial_cone_height_checks import clean, normal_form
from run_algebraic_incidence_controls import pair_checks
from run_structural_complexity_checks import area2, clip, cross, determinant, dot, greedy, hyperbolic_cross, primitive
from run_extremal_initial_ideals_checks import adjacent_weights, cone_for_basis, integer_weight, inverse_mod

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


def depth_budget_controls():
    """Test the profile inequality against all subsets, without matroid ranks."""
    counts, digest = defaultdict(int), hashlib.sha256()
    for d, last in ((2, 9), (3, 10), (4, 6)):
        states = {frozenset({(0,)*d})}
        for size in range(1, last+1):
            for D in sorted(states, key=lambda s: tuple(sorted(s))):
                counts['downsets'] += 1
                for h in product((1, 2, 3), repeat=d):
                    T = {a for a in D if all(tuple(x+(h[i] if j == i else 0)
                                                  for j, x in enumerate(a)) not in D for i in range(d))}
                    Q = D-T
                    assert all(tuple(x-(i == j) for j, x in enumerate(a)) in Q
                               for a in Q for i in range(d) if a[i])
                    for i in range(d):
                        fibers = defaultdict(list)
                        for a in D:
                            fibers[a[:i]+a[i+1:]].append(a)
                        for line in fibers.values():
                            tail = sorted(a[i] for a in line if a in T)
                            core = sorted(a[i] for a in line if a in Q)
                            assert core == list(range(len(core)))
                            assert tail == list(range(len(core), len(line))) and len(tail) <= h[i]
                            counts['coordinate_fibers'] += 1
                        values = sorted((a[i] for a in T), reverse=True)
                        partial = 0
                        for r in range(len(T)+1):
                            if r:
                                partial += values[r-1]
                            bound = h[i]*len(Q)+(h[i]-1)*r
                            assert partial <= bound
                            digest.update(repr((d, tuple(sorted(D)), h, i, r, partial, bound)).encode())
                            counts['all_subset_coordinate_maxima'] += 1
                    counts['anisotropic_supports'] += 1
            if size < last:
                states = {D | {a} for D in states for a in minimal_excluded_border(tuple(D))}
    return {'counts': dict(counts), 'deterministic_certificate_sha256': digest.hexdigest(),
            'dimensions_and_maximum_sizes': [[2, 9], [3, 10], [4, 6]],
            'scope': 'For every support and every residual cardinality, compare the proposed bound with the direct sum of the largest coordinate values. This covers all subsets, including non-downsets and nonbases.'}


def depth_rref(rows, p, order=None):
    rows = [[x % p if p else Fraction(x) for x in row] for row in rows]
    if not rows:
        return [], []
    pivots, k = [], 0
    for j in (range(len(rows[0])) if order is None else order):
        pivot = next((i for i in range(k, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[k], rows[pivot] = rows[pivot], rows[k]
        inv = pow(rows[k][j], -1, p) if p else 1/rows[k][j]
        rows[k] = [(x*inv) % p if p else x*inv for x in rows[k]]
        for i in range(len(rows)):
            if i != k and rows[i][j]:
                c = rows[i][j]
                rows[i] = [(x-c*y) % p if p else x-c*y for x, y in zip(rows[i], rows[k])]
        pivots.append(j)
        k += 1
        if k == len(rows):
            break
    return rows[:k], pivots


def depth_module(D, seeds, p, counts):
    D = tuple(sorted(D))
    index = {a: i for i, a in enumerate(D)}
    m = len(D)

    def shifted(poly, b):
        row = [0]*m
        for a, c in poly.items():
            u = tuple(x+y for x, y in zip(a, b))
            if u in index:
                row[index[u]] += c
        return row

    # Generate the full polynomial ideal in the monomial reference algebra.
    rows, pivots = depth_rref([shifted(f, b) for f in seeds for b in D], p)
    free = [j for j in range(m) if j not in pivots]
    n = len(free)
    assert n > 0 and rows

    def reduce(row):
        row = [x % p if p else Fraction(x) for x in row]
        for v, j in zip(rows, pivots):
            c = row[j]
            if c:
                row = [(x-c*y) % p if p else x-c*y for x, y in zip(row, v)]
        return row

    def residue(a):
        row = [int(index.get(a, -1) == j) for j in range(m)]
        reduced = reduce(row)
        return tuple(reduced[j] for j in free)

    support = {a for row in rows for a, c in zip(D, row) if c}
    h = tuple(max(next(t for t in range(1, m+1)
                       if tuple(x+(t if j == i else 0) for j, x in enumerate(a)) not in index)
                  for a in support) for i in range(3))
    T = {a for a in D if all(tuple(x+(h[i] if j == i else 0) for j, x in enumerate(a)) not in index
                            for i in range(3))}
    Q = set(D)-T
    assert support <= T and len(Q) <= n
    assert all(any(residue(a)) for a in D), 'these covers equal the intrinsic monomial part'
    reference = [D[j] for j in free]
    maps = [[residue(tuple(x+(i == j) for j, x in enumerate(a))) for a in reference] for i in range(3)]

    def apply(i, v):
        return tuple(sum(c*maps[i][j][k] for j, c in enumerate(v)) % p if p
                     else sum(c*maps[i][j][k] for j, c in enumerate(v)) for k in range(n))

    for i, j in combinations(range(3), 2):
        for k in range(n):
            assert apply(i, maps[j][k]) == apply(j, maps[i][k])
            counts['commutator_columns'] += 1
    for a in D:
        v = residue((0, 0, 0))
        for i, power in enumerate(a):
            for _ in range(power):
                v = apply(i, v)
        assert v == residue(a)
        counts['cyclic_monomial_evaluations'] += 1
    for row in rows:
        f = dict(zip(D, row))
        for i in range(3):
            unit = tuple(int(j == i) for j in range(3))
            assert not any(reduce(shifted(f, unit)))
            counts['relation_stability_columns'] += 1
            kill = tuple(h[i]*int(j == i) for j in range(3))
            assert not any(shifted(f, kill))
            counts['annihilation_columns'] += 1
    # Minimality is checked in the reference algebra, before quotienting R out.
    for i in range(3):
        witness = tuple((h[i]-1)*int(j == i) for j in range(3))
        assert any(any(shifted(dict(zip(D, row)), witness)) for row in rows)
    return D, rows, residue, n, h, T, Q


def depth_complete_fan(module, p, counts):
    D, rows, residue, n, h, T, Q = module
    columns = {a: residue(a) for a in D}
    comparisons = sorted({primitive(tuple(x-y for x, y in zip(a, b))) for a, b in combinations(D, 2)
                          if min(x-y for x, y in zip(a, b)) < 0 < max(x-y for x, y in zip(a, b))})
    pending, states = [(1, n+1, (n+1)**2)], {}
    while pending:
        w = pending.pop()
        B = frozenset(greedy(D, columns, w, p))
        assert len(B) == n and Q <= B and B-Q <= T
        assert all(tuple(x-(i == j) for j, x in enumerate(a)) in B for a in B for i in range(3) if a[i])
        if B in states:
            continue
        polygon, normals = cone_for_basis(B, residue, p)
        assert area2(polygon) > 0
        order = sorted(range(len(D)), key=lambda j: (dot(w, D[j]), D[j]), reverse=True)
        _, pivots = depth_rref(rows, p, order)
        assert B == set(D)-{D[j] for j in pivots}, 'independent descending relation-row elimination'
        moment = tuple(sum(a[i] for a in B-Q) for i in range(3))
        budget = tuple(h[i]*len(Q)+(h[i]-1)*(n-len(Q)) for i in range(3))
        assert all(a <= b for a, b in zip(moment, budget))
        basis = tuple(sorted(B))
        inv = inverse_mod([[residue(a)[i] for a in basis] for i in range(n)], p)
        G = []
        for u in minimal_excluded_border(basis):
            v = residue(u)
            coefficients = [sum(x*y for x, y in zip(row, v)) % p for row in inv]
            G.append(clean({u: 1, **{a: -c for a, c in zip(basis, coefficients) if c}}, p))
        counts['ordinary_s_pairs_depth_models'] += pair_checks(G, w, p)
        for f in [*({a: 1} for a in minimal_excluded_border(D)),
                  *(dict(zip(D, row)) for row in rows)]:
            assert not normal_form(clean(f, p), G, w, p)
            counts['ideal_generator_normal_forms'] += 1
        states[B] = {'strict_weight': w, 'standard_basis': basis,
                     'residual_profile': moment, 'coordinate_budget': budget,
                     'simplex_double_area': str(area2(polygon))}
        pending.extend(adjacent_weights(polygon, normals, comparisons))
    assert sum(Fraction(s['simplex_double_area']) for s in states.values()) == 1
    counts['complete_depth_fans'] += 1
    counts['depth_fan_cells'] += len(states)
    return sorted(states.values(), key=lambda s: s['standard_basis'])


def relation_depth_controls():
    counts, fixtures = defaultdict(int), []
    f = {(2, 0, 0): 1, (0, 2, 0): 1, (0, 0, 2): 1, (1, 1, 1): 1}
    g = {(1, 1, 0): 1, (0, 1, 1): 1, (1, 0, 1): 1, (3, 0, 0): 1}
    cases = []
    for S, D, seeds in ((3, 2, [f]), (3, 2, [f, g]),
                         (4, 3, [{(3, 0, 0): 1, (0, 3, 0): 1, (0, 0, 3): 1, (1, 1, 2): 1}]),
                         (4, 2, [f])):
        ambient = {a for a in product(range(S+1), repeat=3) if sum(a) <= S}
        cases.append((f'simplex_S{S}_degree{D}_seeds{len(seeds)}', ambient, seeds))
    tops = ((3, 2, 0), (3, 0, 1), (0, 2, 1))
    ambient = set().union(*(set(product(*(range(scale*(u[i]+1)) for i, scale in enumerate((3, 2, 1)))))
                            for u in ((1, 1, 0), (1, 0, 1), (0, 1, 1))))
    cases.append(('stretched_socle_3_2_1', ambient, [dict.fromkeys(tops, 1)]))
    for label, ambient, seeds in cases:
        for p in (0, 2, 3, 5, 101):
            module = depth_module(ambient, seeds, p, counts)
            D, rows, residue, n, h, T, Q = module
            if label == 'stretched_socle_3_2_1':
                assert h == (3, 2, 1) and len(D) == 42 and n == 36
            elif len(seeds) == 1:
                degree = min(sum(a) for a in seeds[0])
                S = max(sum(a) for a in ambient)
                assert len(rows) == comb(S-degree+3, 3) and h == (S-degree+1,)*3
            record = {'model': label, 'characteristic': p, 'ambient_size': len(D), 'quotient_length': n,
                      'relation_dimension': len(rows), 'intrinsic_coordinate_depths': h,
                      'forced_core_size': len(Q), 'residual_rank': n-len(Q),
                      'coordinate_budget': [h[i]*len(Q)+(h[i]-1)*(n-len(Q)) for i in range(3)],
                      'ideal_seeds': [[{'exponent': a, 'coefficient': c} for a, c in sorted(seed.items())] for seed in seeds]}
            if p and (label.startswith('simplex_S3') or p == 101):
                record['complete_positive_fan'] = depth_complete_fan(module, p, counts)
            fixtures.append(record)
            counts['multilevel_quotient_models'] += 1
    deep = []
    for n in (3, 5, 8, 16):
        D = {a for a in product(range(n), repeat=3) if sum(a) < n}
        # The canonical monomial cover of (y-x,z-x,x^n) is m^n.
        for i in range(3):
            relation = ((1, 0, 0), (0, 1, 0))
            assert all(tuple(a[j]+(n-2)*(j == i) for j in range(3)) in D for a in relation)
            assert all(tuple(a[j]+(n-1)*(j == i) for j in range(3)) not in D for a in relation)
        witnesses = []
        for i in range(3):
            w = tuple(1 if j == i else 2+j for j in range(3))
            B = sorted(tuple(k*(j == i) for j in range(3)) for k in range(n))
            for degree in range(n):
                level = [a for a in D if sum(a) == degree]
                assert min(level, key=lambda a: dot(a, w)) == B[degree]
                counts['deep_diagonal_class_minima'] += 1
            witnesses.append({'strict_weight': w, 'standard_basis': B})
        deep.append({'length': n, 'canonical_cover_size': len(D), 'intrinsic_coordinate_depths': [n-1]*3,
                     'number_initial_ideals': 3, 'all_initial_witnesses': witnesses})
    return {'counts': dict(counts), 'actual_multilevel_models': fixtures, 'large_depth_small_fan_controls': deep,
            'scope': 'Full ideal closures in finite monomial algebras; exact reference annihilator depths, commuting cyclic quotient matrices, row/column greedy comparison and ordinary polynomial division. Complete fans are field-specific bounded controls, not evidence of a uniform asymptotic lower rate.'}


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    symbolic, family, scales, budget = symbolic_matrices(), family_controls(), asymptotic_controls(), socle_budget_controls()
    depth_budget, depth = depth_budget_controls(), relation_depth_controls()
    sources = ('scripts/run_socle_stability_checks.py', 'docs/initial_cone_stability.md',
               'scripts/run_initial_cone_height_checks.py', 'scripts/run_algebraic_incidence_controls.py',
               'scripts/run_extremal_initial_ideals_checks.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py')
    counts = {**family['counts'], 'symbolic_matrix_entry_identities': symbolic['matrix_product_entry_identities'],
              'all_length_padding_models': len(scales['all_length_padding_checks']),
              'socle_budget_downsets': budget['downsets_checked'],
              'socle_budget_chain_points': budget['proper_chain_points_checked'],
              'arbitrary_linear_socle_models': budget['linear_socle_quotient_models'],
              'linear_socle_commutator_columns': budget['linear_socle_commutator_columns'],
              **{f'depth_budget_{k}': v for k, v in depth_budget['counts'].items()},
              **depth['counts']}
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_socle_stability_checks.py', 'counts': counts,
              'symbolic_matrix_certificates': symbolic, 'all_field_family_controls': family,
              'scale_and_padding_controls': scales, 'general_socle_profile_budget': budget,
              'anisotropic_profile_budget': depth_budget, 'propagating_relation_depth_controls': depth,
              'uniform_relative_margin_hypothesis_status': 'refuted_by_actual_all_field_quotients',
              'sharp_worst_relative_margin': {'dimension_2': 'Theta(n^-2)', 'dimension_3': 'Theta(n^-3)'},
              'unrestricted_G3_bound_changed': False, 'novelty_established': False, 'goal_complete': False,
              'scope': 'Uniform proofs show sharp worst-cone stability and refute the individual margin/area premises. A shared residual-profile budget extends the socle count to propagating relations: coordinate parent-annihilator depths h_i give O(n^(3/2)*sqrt(h_1*h_2*h_3)) initials, by the classical lattice-vertex theorem. The canonical monomial cover makes this a necessary depth-product test for local polynomial separation. Large depth alone can still give only three initials. Neither this restricted aggregate theorem nor finite agreement establishes novelty, a new unrestricted exponent or the main endpoint.',
              'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}}
    (ROOT/'results/research/socle_stability_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(counts, indent=2))


if __name__ == '__main__':
    main()
