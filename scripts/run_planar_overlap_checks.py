#!/usr/bin/env python3
"""Exact planar-duality, actual-overlap and sharp gluing-budget checks."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path
import hashlib
import json

import run_complex_basis_checks as original

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def norm(x, p):
    return x % p if p else Q(x)


def reduce_rows(rows, width, p):
    return original.rref(rows, width, p if p else None)


def rank(columns, p):
    return len(reduce_rows(columns, len(columns[0]), p)[1]) if columns else 0


def kernel(rows, width, p):
    R, pivots = reduce_rows(rows, width, p)
    result = []
    for j in range(width):
        if j in pivots:
            continue
        v = [norm(i == j, p) for i in range(width)]
        for row, i in zip(R, pivots):
            v[i] = norm(-row[j], p)
        assert all(norm(sum(x * y for x, y in zip(row, v)), p) == 0 for row in rows)
        result.append(v)
    return result


def apply(A, v, p):
    return tuple(norm(sum(x * y for x, y in zip(row, v)), p) for row in A)


def coordinates(basis, v, p):
    width = len(basis)
    if not width:
        assert not any(v)
        return ()
    rows = [[b[i] for b in basis] + [v[i]] for i in range(len(v))]
    R, pivots = reduce_rows(rows, width + 1, p)
    assert pivots == list(range(width))
    c = tuple(row[-1] for row in R)
    assert tuple(norm(sum(b[i] * t for b, t in zip(basis, c)), p)
                 for i in range(len(v))) == tuple(v)
    return c


def independent(columns, p):
    out = []
    for v in columns:
        if rank(out + [v], p) > len(out):
            out.append(v)
    return out


def partitions(n, cap=None):
    if not n:
        yield ()
    for a in range(min(n, n if cap is None else cap), 0, -1):
        for tail in partitions(n - a, a):
            yield (a,) + tail


def diagram(shape):
    return tuple((i, j) for j, width in enumerate(shape) for i in range(width))


def shift_actions(D):
    index = {a: i for i, a in enumerate(D)}
    actions = []
    for k in range(2):
        A = [[0] * len(D) for _ in D]
        for j, a in enumerate(D):
            target = tuple(x + (i == k) for i, x in enumerate(a))
            if target in index:
                A[index[target]][j] = 1
        actions.append(A)
    return actions


def intertwiner_equations(domain, target):
    a = len(domain[0])
    b = len(target[0])
    rows = []
    for D, T in zip(domain, target):
        for i, j in product(range(b), range(a)):
            row = [0] * (a * b)
            for k in range(a):
                row[i * a + k] += D[k][j]
            for k in range(b):
                row[k * a + j] -= T[i][k]
            rows.append(row)
    return rows


def hom_dimension(domain, target, p):
    width = len(domain[0]) * len(target[0])
    return width - len(reduce_rows(intertwiner_equations(domain, target), width, p)[1])


def heads(E):
    if not E:
        return ((0, 0),)
    E = set(E)
    border = {(a + 1, b) for a, b in E} | {(a, b + 1) for a, b in E}
    return tuple(sorted((a, b) for a, b in border - E
                        if (not a or (a - 1, b) in E) and
                           (not b or (a, b - 1) in E)))


def ideal_hom_dimension(E, D, p):
    """All pair-lcm syzygies, independently of finite-module Hom matrices."""
    G = heads(E)
    index = {a: i for i, a in enumerate(D)}
    width = len(G) * len(D)
    equations = []
    for i, j in combinations(range(len(G)), 2):
        lcm = tuple(max(a, b) for a, b in zip(G[i], G[j]))
        rows = {}
        for k, sign in ((i, 1), (j, -1)):
            delta = tuple(a - b for a, b in zip(lcm, G[k]))
            for t, a in enumerate(D):
                image = tuple(a[l] + delta[l] for l in range(2))
                if image in index:
                    row = rows.setdefault(image, [0] * width)
                    row[k * len(D) + t] += sign
        equations.extend(rows.values())
    return width - len(reduce_rows(equations, width, p)[1])


def planar_identities():
    diagrams = [diagram(shape) for n in range(6) for shape in partitions(n)]
    records = []
    for p in (0, 2, 3):
        for D in diagrams:
            C = shift_actions(D)
            for E in diagrams:
                if len(E) > 4:
                    continue
                F = shift_actions(E)
                value = ideal_hom_dimension(E, D, p)
                dual = hom_dimension(C, F, p)
                assert value == len(D) + dual
                COUNTS['independent_syzygy_and_dual_Hom_identities'] += 1
                if set(D) <= set(E):
                    assert dual <= len(E) and value <= len(D) + len(E)
                    if D == E:
                        assert value == 2 * len(D)
                    COUNTS['fixed_endpoint_extension_parameter_bounds'] += 1
                if set(E) <= set(D):
                    assert value == len(D) + len(E)
                    W = tuple(a for a in D if a not in E)
                    dim_end = hom_dimension(shift_actions(W), shift_actions(W), p)
                    assert dim_end <= len(D) + len(E)
                    COUNTS['planar_ideal_endomorphism_bounds'] += 1
                if 0 < len(D) + len(E) <= 5:
                    H = intertwiner_equations(F, C)
                    width = len(E) * len(D)
                    # H=[D_y;D_z]. The last Koszul differential is [D_z,-D_y].
                    Dy, Dz = H[:width], H[width:]
                    last = [[*a, *[-x for x in b]] for a, b in zip(Dz, Dy)]
                    h2 = width - len(reduce_rows(last, 2 * width, p)[1])
                    assert h2 == hom_dimension(C, F, p)
                    COUNTS['independent_last_Koszul_cokernel_dualities'] += 1
        records.append({'characteristic': p, 'diagram_count': len(diagrams)})
    return records


def actual_overlap(matrices, unit, planar_ground, p, kind):
    n = len(unit)
    def ev(a):
        v = unit
        for A, power in zip(matrices, a):
            for _ in range(power):
                v = apply(A, v, p)
        return v
    L = independent([ev((0, *a)) for a in planar_ground], p)
    U = independent([ev((1, *a)) for a in planar_ground], p)
    assert rank(L + U, p) == n
    rows = [[*[v[i] for v in L], *[-v[i] for v in U]] for i in range(n)]
    relations = kernel(rows, len(L) + len(U), p)
    W = independent([tuple(norm(sum(c * v[i] for c, v in zip(r, L)), p)
                           for i in range(n)) for r in relations], p)
    h = len(W)
    assert len(L) + len(U) - h == n
    actions = []
    for A in matrices:
        images = [coordinates(W, apply(A, v, p), p) for v in W]
        actions.append([[v[i] for v in images] for i in range(h)])
    dim_end = hom_dimension(actions[1:], actions[1:], p)
    assert dim_end <= n
    if kind == 'square_zero':
        assert len(L) - len(U) >= h and dim_end <= n - h
    for i, j in combinations(range(3), 2):
        for k in range(n):
            v = tuple(int(k == l) for l in range(n))
            assert apply(matrices[i], apply(matrices[j], v, p), p) == apply(matrices[j], apply(matrices[i], v, p), p)
            COUNTS['independent_original_commuting_columns'] += 1
    for v in W:
        for A in matrices:
            coordinates(W, apply(A, v, p), p)
            COUNTS['actual_overlap_ideal_stability_checks'] += 1
    assert all(not any(row) for A in actions[1:] for row in A)
    COUNTS['actual_overlap_endomorphism_bounds'] += 1
    return {'length': n, 'lower_dimension': len(L), 'upper_dimension': len(U),
            'overlap_dimension': h, 'endomorphism_dimension': dim_end}


def square_zero_sharp(h, p, shear):
    C = [[int(i == j) + int(shear and j == i + 1) for j in range(h)] for i in range(h)]
    lower = [(0, h - 1 - j, j + 1) for j in range(h)]
    upper = [(1, h - 1 - i, i) for i in range(h)]
    relations = [{(2, 0, 0): 1}, {(0, h, 0): 1}]
    for j in range(h):
        relation = {lower[j]: 1}
        relation.update({upper[i]: -C[i][j] for i in range(h) if C[i][j]})
        relations.append(relation)
    D, free, ev = original.quotient(h + 1, relations, p if p else None)
    n = len(free)
    reference = [D[i] for i in free]
    matrices = []
    for k in range(3):
        columns = [ev(tuple(x + (i == k) for i, x in enumerate(a))) for a in reference]
        matrices.append([[v[i] for v in columns] for i in range(n)])
    data = actual_overlap(matrices, ev((0, 0, 0)), diagram((h + 1,) * (h + 1)), p, 'square_zero')
    assert data['length'] == h * h + h
    assert data['overlap_dimension'] == h
    assert data['endomorphism_dimension'] == h * h == n - h
    data.update(characteristic=p, h=h, shear=shear, gluing_matrix=C)
    COUNTS['sharp_arbitrary_gluing_realizations'] += 1
    return data


def split_sharp(h, p):
    D = tuple((i, j) for i, j in product(range(h), repeat=2) if i + j < h)
    socle = sorted(a for a in D if sum(a) == h - 1)
    V0, V1 = set(socle[:h // 2]), set(socle[h // 2:])
    factors = [tuple(a for a in D if a not in V1), tuple(a for a in D if a not in V0)]
    reference = [(k, a) for k in (0, 1) for a in factors[k]]
    index = {a: i for i, a in enumerate(reference)}
    n = len(reference)
    matrices = [[[0] * n for _ in range(n)] for _ in range(3)]
    unit = [0] * n
    for j, (k, a) in enumerate(reference):
        matrices[0][j][j] = norm(1 if k == 0 else -1, p)
        for l in (0, 1):
            b = tuple(x + (i == l) for i, x in enumerate(a))
            if (k, b) in index:
                matrices[l + 1][index[k, b]][j] = 1
        if a == (0, 0):
            unit[j] = 1
    data = actual_overlap(matrices, tuple(unit), D, p, 'split_unit')
    assert data['length'] == h * h and data['overlap_dimension'] == h
    assert data['endomorphism_dimension'] == h * h
    data.update(characteristic=p, h=h)
    COUNTS['sharp_general_quadratic_realizations'] += 1
    return data


def module_controls():
    records = []
    for p in (0, 2, 3):
        for r in range(1, 4):
            n = 2 * r + 1
            Y = [[0] * n for _ in range(n)]
            Z = [[0] * n for _ in range(n)]
            for i in range(r):
                Y[r + 1 + i][i] = 1
                Z[r + 1 + i][i + 1] = 1
            e = hom_dimension((Y, Z), (Y, Z), p)
            assert e == r * (r + 1) + 1 and 9 * e >= n * n
            columns = [tuple(A[i][j] for i in range(n)) for A in (Y, Z) for j in range(n)]
            top = n - rank(columns, p)
            socle = n - len(reduce_rows(Y + Z, n, p)[1])
            identity = [[int(i == j) for j in range(n)] for i in range(n)]
            carrier = rank([tuple(x for row in A for x in row) for A in (identity, Y, Z)], p)
            assert (top, socle, carrier) == (r + 1, r, 3)
            assert n <= carrier * min(top, socle) and e >= top * socle
            COUNTS['independent_top_socle_and_carrier_budget_checks'] += 1
            if r >= 2:
                assert 3 * e < n * n
                COUNTS['nonsemisimple_linear_denominator_counterexamples'] += 1
            records.append({'characteristic': p, 'length': n, 'endomorphism_dimension': e,
                            'carrier_length': 3, 'T': 9})
            COUNTS['nonsemisimple_module_lower_budget_checks'] += 1
    for m in range(1, 4):
        n = 2 * m
        Y = [[0] * n for _ in range(n)]
        for k in range(m):
            Y[2 * k + 1][2 * k] = 1
            Y[2 * k][2 * k + 1] = 1
            Y[2 * k + 1][2 * k + 1] = 1
        Z = [[0] * n for _ in range(n)]
        e = hom_dimension((Y, Z), (Y, Z), 2)
        assert e == 2 * m * m and e * 2 == n * n
        COUNTS['nonrational_semisimple_residue_field_equalities'] += 1
    return records


def hook_budget():
    records = []
    for t in (2, 3, 4, 6):
        a = t * t + 1
        D = {(i, j) for i, j in product(range(a), repeat=2) if i < t or j < t}
        left = {(i, j) for i, j in D if i >= t}
        right = {(i, j) for i, j in D if j >= t}
        assert not left & right and len(left | right) == 2 * t * (a - t)
        sides = ((a - t, t), (t, a - t))
        blocks = []
        for source in sides:
            row = []
            for target in sides:
                E = tuple(product(range(target[0]), range(target[1])))
                equations = []
                for k in (0, 1):
                    for output in E:
                        equations.append([int(tuple(v[i] + (source[k] if i == k else 0)
                                                        for i in (0, 1)) == output)
                                          for v in E])
                value = len(E) - len(reduce_rows(equations, len(E), 2)[1])
                assert value == min(source[0], target[0]) * min(source[1], target[1])
                row.append(value)
                COUNTS['independent_rectangular_Hom_annihilator_counts'] += 1
            blocks.append(row)
        h = len(left | right)
        n = 2 * len(D)
        e = sum(map(sum, blocks))
        assert e == n - h == 2 * t * a
        if t == 2:
            T = tuple(sorted(left | right))
            for p in (0, 2, 3):
                actions = shift_actions(T)
                assert hom_dimension(actions, actions, p) == e
                COUNTS['independent_whole_hook_commutant_equalities'] += 1
        max_relative_degree = max(i - t + j for i, j in left)
        assert max_relative_degree + 1 == a - 1
        if t >= 3:
            A0 = a - t
            left_images = tuple(product(range(A0 - t, A0), range(t)))
            right_images = tuple(product(range(t), range(A0 - t, A0)))
            assert len(left_images) == len(right_images) == t * t
            for v, w in product(left_images, right_images):
                exponent = tuple(x + y for x, y in zip(v, w))
                assert exponent[0] >= t and exponent[1] >= t
                COUNTS['exact_off_diagonal_zero_compositions'] += 2
        records.append({'t': t, 'a': a, 'length': n, 'overlap_dimension': h,
                        'endomorphism_dimension': e, 'block_Hom_dimensions': blocks,
                        'support_degree': 1, 'Loewy_order': a - 1,
                        'off_diagonal_ideal_square_zero': t >= 3})
        COUNTS['generic_hook_linear_budget_equalities'] += 1
    return records


def main():
    planar = planar_identities()
    square = [square_zero_sharp(h, p, shear) for p in (0, 2, 3)
              for h in range(1, 5) for shear in (False, True)]
    split = [split_sharp(h, p) for p in (0, 3, 5) for h in range(2, 5)]
    modules = module_controls()
    hooks = hook_budget()
    sources = ['docs/planar_overlap_endomorphisms.md', 'scripts/run_planar_overlap_checks.py',
               'scripts/run_complex_basis_checks.py', 'docs/quadratic_intersection_fans.md',
               'docs/profile_section_budget.md', 'docs/square_zero_extension_fans.md',
               'docs/cross_layer_direction_budget.md', 'docs/disjoint_socle_circuits.md',
               'docs/generic_hook_fans.md', 'docs/hook_colour_width.md']
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_planar_overlap_checks.py',
              'scope': 'Finite exact checks of a uniform planar-duality resource reduction and sharp actual overlap realizations. The proof, not fitting, gives the bounded-support/thickness count consequence. No unrestricted exponent, exact novelty or goal completion.',
              'counts': dict(sorted(COUNTS.items())), 'planar_controls': planar,
              'sharp_square_zero': square, 'sharp_split_quadratic': split,
              'nonsemisimple_module_controls': modules,
              'generic_hook_budgets': hooks,
              'proof_scope': {'upper': 'End_B(W)<=n, or <=n-h for square-zero x; actual ideal W only.',
                              'lower': 'End_S(W)>=h^2/T(S), T=sum d_i*ell_i^2; no linear denominator for general nonsemisimple modules.',
                              'fan_consequence': 'If R^r W=0 for radical planar R of colength s, g=O((1+binom(r+1,2)*sqrt(s*n))*(1+n^(4/3))).',
                              'guards': 'A parameter-space dimension is not a direction count. Free matrix mixing costs h^2+h in the square-zero class; no analogous cost asserted for a single prescribed matrix.'},
              'primary_sources': [{'url': 'https://stacks.math.columbia.edu/tag/062F', 'inspected': 'Regular sequence implies Koszul exactness, with proof'},
                                  {'url': 'https://stacks.math.columbia.edu/tag/0621', 'inspected': 'Koszul differential and two-operator complex'},
                                  {'url': 'https://stacks.math.columbia.edu/tag/00NN', 'inspected': 'Lemma 10.106.1 associated graded regular local ring'}],
              'source_sha256': {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in sources}}
    (ROOT / 'results/research/planar_overlap_checks.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'counts': report['counts'], 'source_hashes': len(sources)}))


if __name__ == '__main__':
    main()
