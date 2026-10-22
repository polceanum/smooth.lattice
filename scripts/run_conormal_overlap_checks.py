#!/usr/bin/env python3
"""Exact conormal traces, actual overlap images and all-field rank witnesses."""
import ast
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import hashlib
import json

from run_signed_face_checks import field_basis
from run_square_zero_extension_checks import (
    extension, hom_space, matrix_controls,
)
from run_structural_complexity_checks import determinant

ROOT = Path(__file__).resolve().parents[1]


def partitions(n, cap=None):
    if not n:
        yield ()
        return
    for a in range(min(n, n if cap is None else cap), 0, -1):
        for tail in partitions(n-a, a):
            yield (a,)+tail


def diagram(shape):
    return frozenset((i, j) for j, a in enumerate(shape) for i in range(a))


def heads(D):
    """Independent staircase boundary, in decreasing y-exponent order."""
    rows = [sum(j == z for _, j in D) for z in range(1+max(j for _, j in D))]
    return tuple([(rows[0], 0)] +
                 [(a, j) for j, a in enumerate(rows[1:], 1) if a < rows[j-1]] +
                 [(0, len(rows))])


def plus(a, b):
    return tuple(x+y for x, y in zip(a, b))


def shifted(f, delta, D):
    return {plus(v, delta): c for v, c in f.items() if plus(v, delta) in D}


def valid_map(D, G, tails):
    """Check every pair syzygy integrally, independently of a nullspace."""
    for i, j in combinations(range(len(G)), 2):
        lcm = tuple(max(a, b) for a, b in zip(G[i], G[j]))
        si = tuple(a-b for a, b in zip(lcm, G[i]))
        sj = tuple(a-b for a, b in zip(lcm, G[j]))
        assert shifted(tails[i], si, D) == shifted(tails[j], sj, D)


def arrow_basis(D, G):
    """Fine-degree adjacent equations: connected equalities and pinned zeros.

    No matrix nullspace, determinant or finite-field arithmetic is used.
    An entry is (generator index, output monomial), with coefficient one.
    """
    degrees = sorted({(v[0]-g[0], v[1]-g[1]) for g in G for v in D})
    out = []
    for delta in degrees:
        nodes = {i: plus(g, delta) for i, g in enumerate(G) if plus(g, delta) in D}
        parent = {i: i for i in nodes}

        def root(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        pinned = set()
        for i in range(len(G)-1):
            common = (G[i][0]+delta[0], G[i+1][1]+delta[1])
            if common not in D:
                continue
            present = [j for j in (i, i+1) if j in nodes]
            if len(present) == 2:
                parent[root(i+1)] = root(i)
            elif present:
                pinned.add(present[0])
        bad = {root(i) for i in pinned}
        components = defaultdict(list)
        for i, v in sorted(nodes.items()):
            if root(i) not in bad:
                components[root(i)].append((i, v))
        out.extend(tuple(v) for _, v in sorted(components.items()))
    return tuple(out)


def poly_determinant(A):
    """Independent Leibniz determinant in Z[y,z], for small sparse matrices."""
    n = len(A)
    out = defaultdict(int)
    for sigma in permutations(range(n)):
        coefficient = (-1)**sum(sigma[i] > sigma[j]
                                 for i in range(n) for j in range(i+1, n))
        exponent = (0, 0)
        for i, j in enumerate(sigma):
            if A[i][j] is None:
                coefficient = 0
                break
            c, v = A[i][j]
            coefficient *= c
            exponent = plus(exponent, v)
        out[exponent] += coefficient
    return {v: c for v, c in out.items() if c}


def minors(G, counts):
    r = len(G)-1
    X = [[None]*r for _ in G]
    for i in range(r):
        X[i][i] = (1, (0, G[i+1][1]-G[i][1]))
        X[i+1][i] = (-1, (G[i][0]-G[i+1][0], 0))
    generators = set()
    for k in range(r):
        columns = [i for i in range(r) if i != k]
        for p, q in combinations(range(r+1), 2):
            rows = [i for i in range(r+1) if i not in (p, q)]
            value = poly_determinant([[X[i][j] for j in columns] for i in rows])
            if p <= k < q:
                expected = (G[p][0]-G[k][0]+G[q][0],
                            G[p][1]+G[q][1]-G[k+1][1])
                assert len(value) == 1 and abs(value.get(expected, 0)) == 1
                generators.add(expected)
            else:
                assert not value
            counts['direct_polynomial_submaximal_minor_equalities'] += 1
    return frozenset(generators)


def cofactor_maps(D, G):
    """Differentiate each maximal minor with respect to each matrix entry."""
    r = len(G)-1
    X = [[None]*r for _ in G]
    for i in range(r):
        X[i][i] = (1, (0, G[i+1][1]-G[i][1]))
        X[i+1][i] = (-1, (G[i][0]-G[i+1][0], 0))
    out = []
    for a, b in product(range(r+1), range(r)):
        tails = []
        for i in range(r+1):
            if i == a:
                tails.append({})
                continue
            rows = [j for j in range(r+1) if j not in (i, a)]
            columns = [j for j in range(r) if j != b]
            value = poly_determinant([[X[j][k] for k in columns] for j in rows])
            sign = (-1)**(r+i+a-int(a > i)+b)
            tails.append({v: sign*c for v, c in value.items() if v in D})
        valid_map(D, G, tails)
        out.append(tails)
    return out


def ideal_support(D, generators):
    return frozenset(v for v in D
                     if any(all(a <= b for a, b in zip(g, v)) for g in generators))


def small_controls(counts):
    rows = []
    for n in range(1, 13):
        for shape in partitions(n):
            D = diagram(shape)
            G = heads(D)
            arrows = arrow_basis(D, G)
            trace = frozenset(v for arrow in arrows for _, v in arrow)
            predicted = ideal_support(D, minors(G, counts))
            assert trace == predicted
            # The graph vectors obey all pair syzygies and are disjoint in
            # their coefficient positions, hence integrally independent.
            positions = []
            for arrow in arrows:
                tails = [dict() for _ in G]
                for i, v in arrow:
                    tails[i][v] = 1
                valid_map(D, G, tails)
                positions.extend(arrow)
                counts['integral_arrow_maps_all_pair_syzygies'] += 1
            assert len(set(positions)) == len(positions)
            # Independent derivative/cofactor construction. Its C-span,
            # rather than just its coordinate ideal, must fill Hom.
            differential = cofactor_maps(D, G)
            products = [tuple(shifted(f, v, D) for f in tails)
                        for tails in differential for v in sorted(D)]
            dimensions = []
            for p in (2, 3, 101):
                other_G, variables, dense = hom_space(D, D, p)
                assert set(other_G) == set(G)
                lookup = {(G.index(other_G[i]), v): j
                          for j, (i, v) in enumerate(variables)}
                columns = [{lookup[i, v]: 1 for i, v in a} for a in arrows]
                rank = len(field_basis(range(len(columns)), lambda j: columns[j], p)[0])
                assert rank == len(arrows) == len(dense) == 2*n
                dense_trace = {v for vector in dense for (_, v), c in zip(variables, vector) if c}
                assert dense_trace == trace
                derivative_columns = [{lookup[i, v]: c for i, f in enumerate(tails)
                                       for v, c in f.items()} for tails in products]
                derivative_rank = len(field_basis(range(len(derivative_columns)),
                                                   lambda j: derivative_columns[j], p)[0])
                assert derivative_rank == len(dense)
                counts['independent_cofactor_module_vs_entire_hom_space'] += 1
                dimensions.append(len(dense))
                counts['independent_dense_all_pair_vs_adjacent_graph_spaces'] += 1
            rows.append({'row_lengths': shape, 'length': n, 'minimal_heads': G,
                         'hom_dimensions_in_characteristics_2_3_101': dimensions,
                         'trace_length': len(trace), 'trace_monomials': sorted(trace)})
        print(json.dumps({'checked_diagrams_through_length': n, 'cases': len(rows)}), flush=True)
    assert len(rows) == 271
    return rows


def bareiss(A):
    """Integer fraction-free elimination, separate from rational determinant."""
    A = [list(row) for row in A]
    if not A:
        return 1
    sign, previous = 1, 1
    for k in range(len(A)-1):
        pivot_row = next((i for i in range(k, len(A)) if A[i][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            A[k], A[pivot_row] = A[pivot_row], A[k]
            sign *= -1
        pivot = A[k][k]
        for i in range(k+1, len(A)):
            for j in range(k+1, len(A)):
                value = A[i][j]*pivot-A[i][k]*A[k][j]
                assert value % previous == 0
                A[i][j] = value//previous
            A[i][k] = 0
        previous = pivot
    return sign*A[-1][-1]


def witness(s):
    shape = (9, 9, 9, 7, 7, 4, 4, 4, 4)
    D = frozenset((s*i+a, s*j+b) for i, j in diagram(shape)
                  for a, b in product(range(s), repeat=2))
    G = tuple((s*a, s*b) for a, b in ((9, 0), (7, 3), (4, 5), (0, 9)))
    supports = (((5, 0), (2, 4)), ((3, 3), (7, 0), (0, 7)),
                ((0, 5), (4, 2), (6, 0)), ((0, 6), (2, 4)))
    tails = [{(s*a, s*b): 1 for a, b in support} for support in supports]
    return D, G, tails


def direct_image(D, G, tails, plain):
    """Independent map evaluation, with output rows labelled by D itself."""
    images = []
    for g, tail in zip(G, tails):
        if all(a <= b for a, b in zip(g, plain)):
            images.append(shifted(tail, tuple(b-a for a, b in zip(g, plain)), D))
    assert images and all(f == images[0] for f in images)
    return images[0]


def gap_controls(counts):
    rows, certificates = [], []
    for s in (1, 2, 3):
        D, G, tails = witness(s)
        valid_map(D, G, tails)
        trace = ideal_support(D, minors(G, counts))
        assert len(D) == 57*s*s and len(trace) == 38*s*s
        qtrace = {plus(v, delta) for v in trace for delta in ((s, 0), (0, s))
                  if plus(v, delta) in D}
        qring = {plus(v, delta) for v in D for delta in ((s, 0), (0, s))
                 if plus(v, delta) in D}
        assert len(trace-qtrace) == 5*s*s and len(D-qring) == s*s
        bound = len(trace)-len(trace-qtrace)+len(G)*len(D-qring)
        assert bound == 37*s*s
        n, ev, ground = extension(D, D, G, tails, 0)
        standard = tuple([(0, *v) for v in sorted(D)]+[(1, *v) for v in sorted(D)])
        matrix_controls(n, ev, standard, ground, 0, counts)
        plain = tuple(v for v in ground if v[0] == 0 and v[1:] not in D)
        index = {v: j for j, v in enumerate(sorted(D))}
        columns = []
        for u in plain:
            image = direct_image(D, G, tails, u[1:])
            column = {index[v]: c for v, c in image.items()}
            assert ev(u) == {len(D)+i: c for i, c in column.items()}
            assert set(image) <= trace
            columns.append(column)
            counts['independent_extension_image_columns'] += 1
        ranks = []
        for p in (2, 3, 101):
            selected = field_basis(range(len(columns)), lambda j: columns[j], p)[0]
            assert len(selected) == bound
            # Directly generating the ideal from all four tails is a
            # different column ground from D(M^2) minus D(M).
            products = [shifted(f, v, D) for f in tails for v in sorted(D)]
            direct_rank = len(field_basis(range(len(products)),
                                         lambda j: {index[v]: c for v, c in products[j].items()}, p)[0])
            assert direct_rank == len(selected)
            counts['image_ideal_rank_vs_extension_overlap_rank'] += 1
            ranks.append({'characteristic': p, 'rank': len(selected)})
            if s == 1 and p in (2, 3):
                row_vectors = [{j: columns[c].get(i, 0) for j, c in enumerate(selected)
                                if columns[c].get(i, 0)} for i in range(len(D))]
                selected_rows = field_basis(range(len(D)), lambda i: row_vectors[i], p)[0]
                A = [[columns[c].get(i, 0) for c in selected] for i in selected_rows]
                det = bareiss(A)
                assert det == determinant(A)
                assert det % p
                certificates.append({'pivot_characteristic': p,
                                     'plain_monomial_columns': [plain[c][1:] for c in selected],
                                     'upper_row_monomials': [sorted(D)[i] for i in selected_rows],
                                     'integer_matrix': A, 'determinant': det})
                counts['independent_exact_integer_determinants'] += 1
        if s == 1:
            assert gcd(*(abs(c['determinant']) for c in certificates)) == 1
            assert [c['determinant'] for c in certificates] == [-3, 4]
        rows.append({'scale': s, 'planar_length': len(D), 'extension_length': n,
                     'minimal_heads': G,
                     'tails': [[[a, b, c] for (a, b), c in sorted(f.items())] for f in tails],
                     'trace_length': len(trace),
                     'trace_mod_parameter_ideal_length': len(trace-qtrace),
                     'ring_mod_parameter_ideal_length': len(D-qring),
                     'uniform_capacity_upper_bound': bound, 'image_ranks': ranks,
                     'rectangle_to_length_ratio': [81, 57]})
        print(json.dumps({'scale': s, 'extension_length': n, 'trace': len(trace),
                          'maximum_overlap': bound}), flush=True)
    return {'scaled_actual_extensions': rows, 'all_field_rank_certificates': certificates,
            'uniform_statement': 'For every field and s>=1, max image dimension=37*s^2, trace length=38*s^2, actual extension length=114*s^2.',
            'proof_of_upper_bound': 'CT7 with q=(y^s,z^s): the trace quotient has length 5*s^2 and each of four generators contributes at most s^2.',
            'proof_of_lower_bound': 'The displayed integral map has 37-minors -3 and 4, hence rank at least 37 over every field; free base change multiplies length by s^2.',
            'scope': 'Bounded rectangle ratio already implies O(n^2) fans for every deformation. No supergeneric fan growth follows.'}


def dense_power_controls(counts):
    rows = []
    for d in range(1, 9):
        D = frozenset((a, b) for a in range(d) for b in range(d-a))
        G = heads(D)
        arrows = arrow_basis(D, G)
        trace = frozenset(v for arrow in arrows for _, v in arrow)
        assert trace == frozenset(v for v in D if sum(v) == d-1)
        tails = [{v: 1} for v in sorted(trace)]+[{}]
        assert len(tails) == len(G)
        valid_map(D, G, tails)
        for i, f in enumerate(tails):
            assert all(sum(G[i]) == 1+sum(v) for v in f)
        rows.append({'power': d, 'extension_length': d*(d+1),
                     'trace_length_and_maximum_overlap': d,
                     'all_maps_ordinary_homogeneous': True,
                     'fan_bound': 'Existing positive-grading theorem gives at most 12*n+6.'})
        counts['dense_power_trace_and_ordinary_grading_controls'] += 1
    return rows


def hook_controls(counts):
    rows = []
    for a, b in ((2, 2), (3, 7), (7, 3), (8, 8), (20, 20)):
        D = frozenset({(i, 0) for i in range(a)} | {(0, j) for j in range(b)})
        G = ((a, 0), (1, 1), (0, b))
        tails = [{(1, 0): 1}, {}, {(0, 1): 1}]
        valid_map(D, G, tails)
        T = ideal_support(D, minors(G, counts))
        assert T == D-{(0, 0)}
        n, ev, ground = extension(D, D, G, tails, 0)
        basis = tuple([(0, *v) for v in sorted(D)]+[(1, *v) for v in sorted(D)])
        matrix_controls(n, ev, basis, ground, 0, counts)
        assert not ev((2, 0, 0)) and not ev((0, 2, 2))
        for p in (2, 3, 101):
            plain = tuple(v for v in ground if v[0] == 0)
            lower_rank = len(field_basis(plain, ev, p)[0])
            assert lower_rank-len(D) == len(T)
            counts['hook_actual_attainable_large_overlap_controls'] += 1
        rows.append({'axis_lengths': [a, b], 'extension_length': n,
                     'trace_length_and_attained_overlap': len(T),
                     'rectangle_ratio': [a*b, len(D)],
                     'common_monomial_carrier': [[2, 0, 0], [0, 2, 2]],
                     'fan_bound': 'Existing fixed one-dimensional monomial carrier theorem gives O(n^2) for every deformation.'})
    return rows


def source_hashes():
    paths = ['scripts/run_conormal_overlap_checks.py', 'docs/conormal_overlap.md',
             'docs/square_zero_extension_fans.md', 'docs/profile_section_budget.md',
             'docs/nonnegative_graded_fans.md', 'docs/univariate_cross_fans.md',
             'docs/quadratic_intersection_fans.md']
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
    before = source_hashes()
    counts = Counter()
    small = small_controls(counts)
    gap = gap_controls(counts)
    powers = dense_power_controls(counts)
    hooks = hook_controls(counts)
    assert before == source_hashes(), 'Source files changed during verification.'
    report = {'command': 'python3 -B scripts/run_conormal_overlap_checks.py',
              'status': 'passed', 'source_sha256': before, 'counts': dict(counts),
              'primary_sources': [
                  {'url': 'https://arxiv.org/pdf/1802.06491v2',
                   'scope': 'Definition 2.1 and Remark 2.3: trace and transpose-kernel coordinates.'},
                  {'url': 'https://arxiv.org/pdf/2406.07517v1',
                   'scope': 'Theorems 3.1--3.2 and Lemma 3.3: canonical trace, kernel of X, bidiagonal presentation. The separate conormal proof uses kernel of transpose X.'},
                  {'url': 'https://arxiv.org/pdf/2112.04237v2',
                   'scope': 'Introduction: Teter number counts canonical-module maps needed to cover the canonical trace. Individual image questions are prior work, with a different domain module.'},
                  {'url': 'https://math.berkeley.edu/~robin/math274root.pdf',
                   'scope': 'Theorem 8.5(3), attributed there to Schaps: first-order codimension-two deformations lift the Hilbert-Burch matrix. Differentiating maximal minors gives the conormal trace identity as a classical corollary.'},
                  {'url': 'https://arxiv.org/pdf/2407.07993v4',
                   'scope': 'Section 5, Theorem 5.4 and proof: explicit arbitrary-field matrix parametrization gives every first-order deformation of a planar monomial ideal.'}],
              'small_planar_diagrams': small, 'sharp_nonattainment_family': gap,
              'dense_power_scope_controls': powers,
              'large_overlap_and_rectangle_ratio_guards': hooks,
              'limits': ['The complete uniform proofs are in docs/conormal_overlap.md; finite checks do not prove all-field or all-scale statements.',
                         'The conormal trace identity is a classical determinantal-deformation corollary; canonical trace and the presentation-kernel framework are established prior work. Exact priority of the capacity application and sharp family is unassessed.',
                         'The trace is a sum of all deformation images; it need not equal any single image, including a generic one.',
                         'No unrestricted G3 improvement, polynomial lower separation, consequential general complexity theorem or completed goal follows.']}
    target = ROOT/'results/research/conormal_overlap_checks.json'
    target.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
