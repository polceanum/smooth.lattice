#!/usr/bin/env python3
"""Exact balanced-cone realizations and checks; no Monte Carlo growth inference."""
import ast
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
from math import prod
from pathlib import Path
import hashlib
import json

from run_algebraic_incidence_controls import pair_checks
from run_initial_cone_height_checks import clean, normal_form
from run_interleaving_research import minimal_excluded_border
from run_projection_tag_checks import compare_fans, fan, polygon_key
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]
TRIANGLE = [tuple(Q(i == j) for i in range(3)) for j in range(3)]


def matmul(A, B):
    return tuple(tuple(sum(a*b for a, b in zip(row, col)) for col in zip(*B))
                 for row in A)


def eye(d):
    return tuple(tuple(int(i == j) for j in range(d)) for i in range(d))


def inverse(A):
    d = len(A)
    rows = [list(map(Q, a))+list(eye(d)[i]) for i, a in enumerate(A)]
    for j in range(d):
        pivot = next(i for i in range(j, d) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scale = rows[j][j]
        rows[j] = [x/scale for x in rows[j]]
        for i in range(d):
            if i != j:
                scale = rows[i][j]
                rows[i] = [x-scale*y for x, y in zip(rows[i], rows[j])]
    assert all(x.denominator == 1 for row in rows for x in row)
    return tuple(tuple(int(x) for x in row[d:]) for row in rows)


def determinant(A):
    if not A:
        return 1
    rows = [list(map(Q, a)) for a in A]
    value = Q(1)
    for j in range(len(A)):
        pivot = next((i for i in range(j, len(A)) if rows[i][j]), None)
        if pivot is None:
            return Q(0)
        if pivot != j:
            rows[pivot], rows[j] = rows[j], rows[pivot]
            value = -value
        scale = rows[j][j]
        value *= scale
        for i in range(j+1, len(A)):
            factor = rows[i][j]/scale
            rows[i] = [x-factor*y for x, y in zip(rows[i], rows[j])]
    return value


def general_matrix(d, counts):
    U = tuple(tuple(int(i == j)-int(i == d-1 and j < d-1)
                    for j in range(d)) for i in range(d))
    V0 = tuple(tuple(int(i == j)-int(i < d-1 and j == i+1)
                     for j in range(d)) for i in range(d))
    J = tuple(tuple(int(j == i+1) for j in range(d)) for i in range(d))
    S = tuple(tuple(int(j == i+1) if i < d-1 else
                    (-1)**(d-1)*int(j == 0) for j in range(d))
              for i in range(d))
    A = tuple(tuple(int(j == i)-int(j == (i+1) % d)
                    for j in range(d)) for i in range(d))
    B = matmul(matmul(U, inverse(S)), V0)
    assert matmul(matmul(matmul(U, J), inverse(S)), V0) == A
    assert determinant(U) == determinant(V0) == determinant(S) == determinant(B) == 1
    power = eye(d)
    coeff = []
    for j in range(d):
        term = matmul(matmul(matmul(inverse(V0), S), power), inverse(U))
        coeff.append(tuple(tuple((-1)**j*x for x in row) for row in term))
        power = matmul(power, J)
    assert not any(x for row in power for x in row)
    assert coeff[-1] == tuple((1,)*d for _ in range(d))
    zero = tuple((0,)*d for _ in range(d))
    for reverse in (False, True):
        for degree in range(d+1):
            total = [list(row) for row in zero]
            for j, C in enumerate((B, A)):
                if 0 <= degree-j < d:
                    left, right = (coeff[degree-j], C) if reverse else (C, coeff[degree-j])
                    for i, row in enumerate(matmul(left, right)):
                        for s, x in enumerate(row):
                            total[i][s] += x
            assert tuple(map(tuple, total)) == (eye(d) if degree == 0 else zero)
            counts['polynomial_matrix_product_coefficients'] += d*d
    return A, B, coeff


def special_matrix(k):
    M = ((k-1, -k, 0), (1, k-1, -k), (-k, 0, k+1))
    V = ((k*k-1, k*k+k, k*k), (k*k-k-1, k*k-1, k*k-k),
         (k*k-k, k*k, k*k-k+1))
    assert matmul(M, V) == matmul(V, M) == eye(3) and determinant(M) == 1
    return M, V


def quotient(tips, padding, counts):
    d = len(tips[0])
    heads, tails = tips[::2], tips[1::2]
    D = set().union(*(set(product(*(range(x+1) for x in a))) for a in tips))
    maximum = max(a[-1] for a in D if not any(a[:-1]))
    D.update((0,)*(d-1)+(maximum+j,) for j in range(1, padding+1))
    assert len(set(tips)) == len(tips)
    assert all(not all(a <= b for a, b in zip(u, v))
               for u, v in product(tips, repeat=2) if u != v)
    assert all(tuple(x+int(i == j) for i, x in enumerate(a)) not in D
               for a in tips for j in range(d))
    assert all(tuple(x-int(i == j) for i, x in enumerate(a)) in D
               for a in D for j in range(d) if a[j])
    reference = tuple(sorted(D-set(heads)))
    index = {a: j for j, a in enumerate(reference)}
    replacement = dict(zip(heads, tails))
    def image(a):
        return index.get(replacement.get(a, a), -1)
    maps = tuple(tuple(image(tuple(x+int(i == j) for i, x in enumerate(a)))
                       for a in reference) for j in range(d))
    for i, j in combinations(range(d), 2):
        for r in range(len(reference)):
            left = -1 if maps[j][r] < 0 else maps[i][maps[j][r]]
            right = -1 if maps[i][r] < 0 else maps[j][maps[i][r]]
            assert left == right
            counts['independent_commuting_multiplication_entries'] += 1
    def evaluate(a, start=None):
        r = index[(0,)*d] if start is None else start
        for j, exponent in enumerate(a):
            for _ in range(exponent):
                if r < 0:
                    return {}
                r = maps[j][r]
        return {} if r < 0 else {r: 1}
    for a in D:
        assert evaluate(a) == {image(a): 1}
        counts['independent_cyclic_monomial_columns'] += 1
    counts['actual_quotient_models'] += 1
    return D, reference, image, maps, evaluate


def sign_witness(M, V, signs):
    positive = [j for j, s in enumerate(signs) if s > 0]
    assert positive and all(x > 0 for row in V for x in row)
    factor = 1+max(sum(row[j] for j, s in enumerate(signs) if s < 0) //
                   sum(row[j] for j in positive) for row in V)
    slack = tuple(factor if s > 0 else -1 for s in signs)
    w = tuple(sum(a*b for a, b in zip(row, slack)) for row in V)
    assert min(w) > 0 and tuple(dot(row, w) for row in M) == slack
    return w


def length_certificate(counts):
    # Each coordinate is (constant, coefficient of k); all minima stay fixed for k>=4.
    tips = (((-1, 1), (1, 0), (1, 0)), ((0, 0), (1, 1), (1, 0)),
            ((2, 0), (-1, 1), (1, 0)), ((1, 0), (0, 0), (1, 1)),
            ((0, 0), (1, 0), (1, 1)), ((0, 1), (1, 0), (0, 0)))
    total = [0]*4
    records = []
    for mask in range(1, 64):
        selected = [tips[j] for j in range(6) if mask >> j & 1]
        minima = []
        term = [(-1)**(len(selected)+1)]
        for i in range(3):
            choice = min((a[i] for a in selected), key=lambda a: (a[0]+4*a[1], a[1]))
            assert all(choice[1] <= a[i][1] and choice[0]+4*choice[1] <= a[i][0]+4*a[i][1]
                       for a in selected)
            minima.append(choice)
            new = [0]*(len(term)+1)
            for j, value in enumerate(term):
                new[j] += (choice[0]+1)*value
                new[j+1] += choice[1]*value
            term = new
        for j, x in enumerate(term):
            total[j] += x
        records.append(dict(tip_mask=mask, minimum_affine_coordinates=minima,
                            signed_divisor_intersection_polynomial=term))
        counts['symbolic_divisor_union_intersections'] += 1
    assert total == [-6, 13, 0, 0]
    return dict(divisor_union_polynomial_ascending=total, subset_certificates=records)


def geometry(M, V, counts):
    d = len(M)
    sums = [sum(row[j] for row in V) for j in range(d)]
    rays = [tuple(Q(row[j], sums[j]) for row in V) for j in range(d)]
    assert all(Q(1, 2*d) < x < Q(2, d) for a in rays for x in a)
    volume = Q(1, prod(sums))
    assert abs(determinant(tuple(zip(*rays)))) == volume
    diameter = max(sum((x-y)**2 for x, y in zip(a, b)) for a, b in combinations(rays, 2))
    for I in combinations(range(d), 2):
        for J in combinations(range(d), 2):
            minor = tuple(tuple(V[i][j] for j in J) for i in I)
            complement = tuple(tuple(M[i][j] for j in range(d) if j not in I)
                               for i in range(d) if i not in J)
            assert determinant(minor) == (-1)**(sum(I)+sum(J))*determinant(complement)
            counts['independent_complementary_inverse_minor_identities'] += 1
    return dict(column_sums=sums, normalized_rays=rays, relative_simplex_volume=volume,
                exact_simplex_hit_waiting_mean=prod(sums), diameter_squared=diameter)


def special_case(k, padding, full, counts):
    M, V = special_matrix(k)
    tips = ((k-1, 1, 1), (0, k+1, 1), (2, k-1, 1), (1, 0, k+1),
            (0, 1, k+1), (k, 1, 0))
    assert tuple(tuple(a-b for a, b in zip(u, v)) for u, v in zip(tips[::2], tips[1::2])) == M
    D, reference, image, maps, ev = quotient(tips, padding, counts)
    n = len(reference)
    assert n == 13*k-9+padding
    predicted, witnesses = {}, []
    for signs in product((-1, 1), repeat=3):
        poly = list(TRIANGLE)
        for q, sign in zip(M, signs):
            poly = clip(poly, q, sign)
        if signs == (-1,)*3:
            assert not poly
            continue
        w = sign_witness(M, V, signs)
        removed = tuple(tips[2*i] if signs[i] > 0 else tips[2*i+1] for i in range(3))
        kept = tuple(tips[2*i+1] if signs[i] > 0 else tips[2*i] for i in range(3))
        B = frozenset(D-set(removed))
        assert len(B) == n and area2(poly) > 0
        border = minimal_excluded_border(tuple(sorted(B)))
        G = [{u: 1, v: -1} for u, v in zip(removed, kept)]
        G += [{a: 1} for a in border if a not in removed]
        if full:
            for char in (0, 2, 3, 101):
                counts['independent_original_reduced_basis_S_pairs'] += pair_checks(G, w, char)
            # Check original polynomial reductions and the independent action model.
            representative = {image(a): a for a in B}
            assert len(representative) == n
            for a in set(D) | set(border):
                expected = {} if image(a) < 0 else {representative[image(a)]: 1}
                assert normal_form({a: 1}, [clean(g, 0) for g in G], w, 0) == expected
                counts['independent_original_monomial_reductions'] += 1
            for g in G:
                for start in range(n):
                    result = Counter()
                    for a, c in g.items():
                        for r, value in ev(a, start).items():
                            result[r] += c*value
                    assert not any(result.values())
                    counts['independent_original_relation_action_entries'] += 1
        predicted[B] = poly
        witnesses.append(dict(signs=signs, strict_integer_weight=w, standard_set=sorted(B),
                              simplex_vertices=poly, relative_simplex_area=area2(poly)))
    assert len(predicted) == 7 and sum(map(area2, predicted.values())) == 1
    geom = geometry(M, V, counts)
    small = next(poly for B, poly in predicted.items() if not set(tips[::2]) & B)
    assert area2(small) == geom['relative_simplex_volume']
    assert polygon_key(small) == polygon_key(geom['normalized_rays'])
    if full:
        for char in (0, 2, 3, 101):
            compare_fans(fan(n, ev, char, counts), predicted, counts)
    counts['entire_seven_cone_sign_arrangements'] += 1
    return dict(k=k, padding=padding, length=n, tips=tips, matrix=M, inverse=V,
                geometry=geom, all_cones=witnesses,
                independent_direct_field_fans=(0, 2, 3, 101) if full else ())


def hashes():
    paths = ['docs/balanced_initial_cones.md', 'scripts/run_balanced_cone_checks.py',
             'docs/initial_cone_stability.md']
    pending = [paths[1]]
    while pending:
        path = pending.pop()
        for node in ast.walk(ast.parse((ROOT/path).read_text())):
            names = ([node.module] if isinstance(node, ast.ImportFrom) and node.module else
                     [a.name for a in node.names] if isinstance(node, ast.Import) else [])
            for name in names:
                candidate = 'scripts/'+name.replace('.', '/')+'.py'
                if (ROOT/candidate).is_file() and candidate not in paths:
                    paths.append(candidate)
                    pending.append(candidate)
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def jsonable(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [jsonable(v) for v in x]
    return x


def main():
    if not __debug__:
        raise RuntimeError('Assertions must be enabled.')
    before = hashes()
    counts = Counter()
    symbolic = length_certificate(counts)
    general = []
    for d in range(3, 9):
        A, B, coefficients = general_matrix(d, counts)
        k = 64
        M = tuple(tuple(B[i][j]+k*A[i][j] for j in range(d)) for i in range(d))
        V = tuple(tuple(sum(k**a*C[i][j] for a, C in enumerate(coefficients))
                        for j in range(d)) for i in range(d))
        assert min(x for row in V for x in row) > 0
        geom = geometry(M, V, counts)
        signs = []
        for pattern in product((-1, 1), repeat=d):
            if pattern != (-1,)*d:
                w = sign_witness(M, V, pattern)
                signs.append(dict(signs=pattern, strict_integer_weight=w))
                counts['independent_all_dimension_sign_witnesses'] += 1
        model_size = None
        if d <= 5:
            T = 1+2*max(abs(x) for row in B for x in row)
            tips = tuple(tuple(max(sign*M[i][j], 0)+T*int(j == (i+2) % d)
                               for j in range(d)) for i in range(d) for sign in (1, -1))
            D, reference, *_ = quotient(tips, 0, counts)
            model_size = len(reference)
        general.append(dict(dimension=d, A=A, B=B, inverse_polynomial_ascending=coefficients,
                            test_parameter=k, geometry=geom, all_sign_witnesses=signs,
                            actual_multiplication_test_length=model_size))
    cases = [special_case(k, 0, True, counts) for k in (4, 7)]
    cases += [special_case(4, r, False, counts) for r in range(1, 13)]
    cases += [special_case(k, 0, False, counts) for k in (16, 32, 64)]
    assert hashes() == before, 'A checked source changed.'
    report = dict(command='python3 -B scripts/run_balanced_cone_checks.py', status='passed',
                  source_sha256=before, counts=dict(counts), uniform_proof='docs/balanced_initial_cones.md',
                  symbolic_divisor_union=symbolic, general_dimension_certificates=general,
                  exact_three_variable_quotients=cases,
                  theorem_scope='Every fixed dimension, every field and every sufficiently large length. Relative cone volume has sharp worst order n^(-d(d-1)); diameter has sharp worst order n^(-d), attained in balanced weights with 2^d-1 cones.',
                  primary_sources=[
                      dict(url='https://doi.org/10.1007/s11786-010-0039-0',
                           scope='Dimitrova publisher abstract: stochastic cone-volume estimation. Full subscription text not inspected; no worst-case assertion is assigned to it.'),
                      dict(url='https://www.tau.ac.il/~hassin/matroids_89.pdf',
                           scope='Hassin--Tamir primary introduction p.362: two objective sums and one scalar slope. This is distinct from two affine cost parameters.')],
                  limits=['Cofactor, adjugate, complementary-minor and sampling identities are classical inputs.',
                          'Exact priority of the combined sharp algebraic realization remains unassessed.',
                          'The sampling bounds concern independent uniform simplex weights or independent isotropic hyperplanes, not adaptive or algebra-dependent choices.',
                          'Seven cones and rare detection do not give a large-fan family or improve unrestricted G3.',
                          'The uniform theorem uses the proof, not finite geometric ratios or Monte Carlo.'])
    (ROOT/'results/research/balanced_cone_checks.json').write_text(json.dumps(jsonable(report), indent=2)+'\n')
    print(json.dumps(dict(status='passed', counts=dict(counts))), flush=True)


if __name__ == '__main__':
    main()
