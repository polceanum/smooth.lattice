#!/usr/bin/env python3
"""Exact bounded checks for adjacent rank fans and finite group gap fans."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from heapq import heappop, heappush
from itertools import combinations, product
from math import comb, gcd, isqrt
from pathlib import Path
import hashlib
import json

from run_interleaving_research import minimal_excluded_border
from run_oracle_enumeration_checks import (
    apply, digit_rank, digit_rank_enumerator, image_from_matrices,
    verify_reduced_basis,
)
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import area2, clip, cross, determinant, dot

ROOT = Path(__file__).resolve().parents[1]
TRIANGLE = [tuple(Q(i == j) for i in range(3)) for j in range(3)]


def rank(vectors, p=0):
    pivots = {}
    for vector in vectors:
        v = [int(x) % p if p else Q(x) for x in vector]
        for j, row in sorted(pivots.items()):
            c = v[j]
            if c:
                v = [(a-c*b) % p if p else a-c*b for a, b in zip(v, row)]
        j = next((j for j, x in enumerate(v) if x), None)
        if j is not None:
            inverse = pow(v[j], -1, p) if p else 1/v[j]
            pivots[j] = [x*inverse % p if p else x*inverse for x in v]
    return len(pivots)


def primitive(v):
    d = gcd(*v)
    if not d:
        return None
    return tuple(x//d for x in v)


def clean(poly):
    poly = list(dict.fromkeys(poly))
    while len(poly) > 2:
        keep = [a for i, a in enumerate(poly)
                if any(cross(tuple(x-y for x, y in zip(a, poly[i-1])),
                             tuple(x-y for x, y in zip(poly[(i+1) % len(poly)], a))))]
        if len(keep) == len(poly):
            break
        poly = keep
    return poly


def intersect(poly, normals):
    for normal in normals:
        poly = clip(poly, normal, 1)
        if not poly:
            return []
    return clean(poly) if area2(poly) else []


def point_fan(points):
    points = sorted(set(points))
    out = []
    for a in points:
        normals = [tuple(x-y for x, y in zip(b, a)) for b in points if b != a]
        poly = intersect(TRIANGLE, normals)
        if poly:
            out.append({'profile': a, 'polygon': poly, 'normals': normals})
    assert sum(area2(c['polygon']) for c in out) == 1
    return out


def overlay(left, right, counts):
    vertices = {v for cell in left+right for v in cell['polygon']}
    result = []
    for a, b in product(left, right):
        poly = intersect(a['polygon'], b['normals'])
        if not poly:
            continue
        assert set(poly) <= vertices
        counts['complete_adjacent_overlay_cells'] += 1
        result.append((a, b, poly))
    assert sum(area2(poly) for _, _, poly in result) == 1
    assert len(result) <= 4*(len(left)+len(right))+2
    return result


def indicators(B, size):
    return tuple(int(i in B) for i in range(size))


def differences(points):
    points = list(points)
    return [tuple(a-b for a, b in zip(v, points[0])) for v in points[1:]]


def matroid_controls():
    counts = defaultdict(int)
    labels = ((0, 0, 0), (4, 0, 0), (0, 4, 0), (0, 0, 4),
              (2, 1, 0), (0, 2, 1), (1, 0, 2), (1, 2, 2), (2, 2, 1))
    graphic = tuple(tuple(int(i == a)-int(i == b) for i in range(3))
                    for a, b in combinations(range(4), 2))
    cases = [('uniform_three_six', tuple((1, a, a*a) for a in range(6)), 0),
             ('graphic_K4', graphic, 2),
             ('loops_parallel_coloops', ((1, 0, 0), (1, 0, 0), (0, 1, 0),
                                         (0, 1, 0), (0, 0, 1), (0, 0, 0)), 3),
             ('binary_pascal', tuple(tuple(comb(a, j) % 2 for j in range(3))
                                      for a in range(8)), 2),
             ('ternary_pascal', tuple(tuple(comb(a, j) % 3 for j in range(4))
                                       for a in range(9)), 3)]
    fixtures = []
    for name, columns, p in cases:
        size = len(columns)
        r = rank(columns, p)
        bases = {k: [frozenset(B) for B in combinations(range(size), k)
                     if rank([columns[i] for i in B], p) == k]
                 for k in (r-1, r)}
        profiles = {k: {B: tuple(sum(labels[i][j] for i in B) for j in range(3))
                         for B in sets} for k, sets in bases.items()}
        normals = {primitive(tuple(a-b for a, b in zip(x, y)))
                   for x, y in combinations(labels[:size], 2)} - {None}
        directions = set(product((-1, 0, 1), repeat=3))
        for a, b in combinations(normals, 2):
            v = primitive(cross(a, b))
            if v:
                directions.update((v, tuple(-x for x in v)))
        for w in sorted(directions):
            faces = {}
            for k in (r-1, r):
                minimum = min(dot(w, v) for v in profiles[k].values())
                faces[k] = [B for B, v in profiles[k].items() if dot(w, v) == minimum]
            ordered = sorted(range(size), key=lambda i: (dot(w, labels[i]), i))
            selected = []
            for i in ordered:
                if rank([columns[j] for j in selected+[i]], p) > len(selected):
                    selected.append(i)
            last = max(dot(w, labels[i]) for i in selected)
            below = [columns[i] for i in range(size) if dot(w, labels[i]) < last]
            s = r-rank(below, p)
            spans = {k: differences([indicators(B, size) for B in faces[k]])
                     for k in (r-1, r)}
            small, large = (r-1, r) if s == 1 else (r, r-1)
            assert rank(spans[large]+spans[small]) == rank(spans[large])
            counts['top_layer_rank_one_faces' if s == 1 else 'top_layer_higher_rank_faces'] += 1
            images = {k: differences([profiles[k][B] for B in faces[k]]) for k in (r-1, r)}
            assert rank(images[r]+images[r-1]) == max(rank(images[r]), rank(images[r-1]))
            if rank(images[r]) == 2 and rank(images[r-1]) < 2:
                assert s == 1
                counts['rank_one_new_ray_controls'] += 1
            counts['independent_projected_face_span_equalities'] += 1
        fans = {k: point_fan(profiles[k].values()) for k in (r-1, r)}
        common = overlay(fans[r], fans[r-1], counts)
        previous_vertices = {v for cell in fans[r-1] for v in cell['polygon']}
        births = []
        for w in sorted({v for cell in fans[r] for v in cell['polygon']}
                        - previous_vertices):
            if min(w) <= 0:
                continue
            best = min(bases[r], key=lambda B: sum(dot(w, labels[i]) for i in B))
            last = max(dot(w, labels[i]) for i in best)
            below = [columns[i] for i in range(size) if dot(w, labels[i]) < last]
            assert r-rank(below, p) == 1
            face_dimensions = []
            for k in (r-1, r):
                value = min(dot(w, v) for v in profiles[k].values())
                face_dimensions.append(rank(differences(
                    [v for v in profiles[k].values() if dot(w, v) == value])))
            assert face_dimensions[0] < 2 and face_dimensions[1] == 2
            counts['complete_positive_new_ray_rank_one_equalities'] += 1
            births.append({'weight': list(map(str, w)),
                           'previous_and_full_face_dimensions': face_dimensions,
                           'last_layer_rank_increment': 1})
        for a, b, poly in common:
            w = tuple(sum(v[i] for v in poly)/len(poly) for i in range(3))
            best = min(bases[r], key=lambda B: sum(dot(w, labels[i]) for i in B))
            selected_costs = sorted(dot(w, labels[i]) for i in best)
            assert dot(w, a['profile'])-dot(w, b['profile']) == selected_costs[-1]
            counts['independent_last_greedy_cost_equalities'] += 1
        fixtures.append({'model': name, 'characteristic': p, 'rank': r,
                         'ground_size': size, 'full_rank_profile_cells': len(fans[r]),
                         'previous_rank_profile_cells': len(fans[r-1]),
                         'common_cells': len(common), 'exact_positive_new_rays': births})
    # Nonadjacent uniform ranks can have genuinely crossing projected edges.
    points = ((1, 0, 0), (0, 1, 0), (2, 0, 0), (0, 0, 2))
    w = (1, 1, 1)
    face1 = [a for a in points if dot(w, a) == 1]
    face3 = [tuple(sum(points[i][j] for i in B) for j in range(3))
             for B in combinations(range(4), 3)
             if sum(dot(w, points[i]) for i in B) == 4]
    assert rank(differences(face1)) == rank(differences(face3)) == 1
    assert rank(differences(face1)+differences(face3)) == 2
    return {'counts': dict(counts), 'fixtures': fixtures,
            'nonadjacent_counterexample': {'ground_labels': points, 'ranks': [1, 3],
                                            'weight': w, 'face_dimensions': [1, 1],
                                            'sum_face_dimension': 2,
                                            'scope': 'General projected matroid data; no polynomial-ideal realization is claimed.'}}


def quotient_pivot_components(columns, removed, labels, w, p):
    """Independent Gaussian contraction and exhaustive minimal circuits."""
    def reduce(v, rows):
        v = [x % p if p else Q(x) for x in v]
        for j, row in sorted(rows.items()):
            c = v[j]
            v = [(a-c*b) % p if p else a-c*b for a, b in zip(v, row)]
        return tuple(v)
    rows = {}
    for cost in sorted({dot(w, a) for a in labels}):
        level = [i for i, a in enumerate(labels) if dot(w, a) == cost]
        vectors = [reduce(columns[i], rows) for i in level]
        line = reduce(removed, rows)
        if rank(vectors+[line], p) == rank(vectors, p):
            parent = list(range(len(level)))
            def root(i):
                while parent[i] != i:
                    i = parent[i]
                return i
            for size in range(1, len(level)+1):
                for S in combinations(range(len(level)), size):
                    if rank([vectors[i] for i in S], p) == size-1 and all(
                            rank([vectors[i] for i in S if i != j], p) == size-1 for j in S):
                        for i in S[1:]:
                            parent[root(i)] = root(S[0])
            groups = defaultdict(list)
            for i in range(len(level)):
                groups[root(i)].append(i)
            single = any(rank([vectors[i] for i in C]+[line], p)
                         == rank([vectors[i] for i in C], p) for C in groups.values())
            return single, [[level[i] for i in C] for C in groups.values()]
        for vector in vectors:
            v = reduce(vector, rows)
            j = next((i for i, x in enumerate(v) if x), None)
            if j is not None:
                inverse = pow(v[j], -1, p) if p else 1/v[j]
                rows[j] = tuple(x*inverse % p if p else x*inverse for x in v)
    raise AssertionError('The removed vector must enter the full column span.')


def elementary_quotient_controls():
    """Direct matrix projections, not generic truncations."""
    counts = defaultdict(int)
    fixtures = []
    labels = ((0, 0, 0), (4, 0, 0), (0, 4, 0),
              (0, 0, 4), (2, 1, 0), (0, 2, 1))
    graphic = tuple(tuple(int(i == a)-int(i == b) for i in range(3))
                    for a, b in combinations(range(4), 2))
    models = [('three_parallel_classes', ((1, 0, 0), (1, 0, 0), (0, 1, 0),
                                         (0, 1, 0), (0, 0, 1), (0, 0, 1)), 3),
              ('graphic_K4', graphic, 2),
              ('rational_uniform', tuple((1, i, i*i) for i in range(6)), 0)]
    for name, columns, p in models:
        for v in ((1, 0, 0), (1, 1, 0), (1, 1, 1)):
            pivot = next(i for i, x in enumerate(v) if x)
            scale = pow(v[pivot], -1, p) if p else Q(1, v[pivot])
            child = [tuple((a[i]-a[pivot]*scale*v[i]) % p if p
                           else a[i]-a[pivot]*scale*v[i]
                           for i in range(3) if i != pivot) for a in columns]
            assert rank(columns, p) == 3 and rank(child, p) == 2
            bases = [[frozenset(B) for B in combinations(range(6), k)
                      if rank([matrix[i] for i in B], p) == k]
                     for matrix, k in ((columns, 3), (child, 2))]
            profiles = [{B: tuple(sum(labels[i][j] for i in B) for j in range(3))
                         for B in sets} for sets in bases]
            normals = {primitive(tuple(a-b for a, b in zip(x, y)))
                       for x, y in combinations(labels, 2)} - {None}
            directions = set(product((-1, 0, 1), repeat=3))
            for a, b in combinations(normals, 2):
                w = primitive(cross(a, b))
                if w:
                    directions.update((w, tuple(-x for x in w)))
            for w in sorted(directions):
                spans = []
                indicator_spans = []
                for data in profiles:
                    value = min(dot(w, a) for a in data.values())
                    spans.append(differences([a for a in data.values() if dot(w, a) == value]))
                    indicator_spans.append(differences(
                        [indicators(B, 6) for B, a in data.items() if dot(w, a) == value]))
                single, components = quotient_pivot_components(columns, v, labels, w, p)
                small, large = (1, 0) if single else (0, 1)
                assert rank(indicator_spans[large]+indicator_spans[small]) == rank(indicator_spans[large])
                counts['single_component_pivot_faces' if single else 'multiple_component_pivot_faces'] += 1
                assert rank(spans[0]+spans[1]) == max(map(rank, spans))
                counts['independent_elementary_quotient_face_equalities'] += 1
            fans = [point_fan(data.values()) for data in profiles]
            common = overlay(fans[0], fans[1], counts)
            parent_vertices = {a for cell in fans[0] for a in cell['polygon']}
            child_births = []
            for w in sorted({a for cell in fans[1] for a in cell['polygon']}-parent_vertices):
                if min(w) > 0:
                    single, components = quotient_pivot_components(columns, v, labels, w, p)
                    assert not single
                    counts['complete_positive_child_ray_merger_equalities'] += 1
                    child_births.append({'weight': list(map(str, w)),
                                         'pivot_layer_components': components})
            fixtures.append({'model': name, 'characteristic': p, 'removed_line': v,
                             'parent_rank': 3, 'child_rank': 2,
                             'parent_cells': len(fans[0]), 'child_cells': len(fans[1]),
                             'common_cells': len(common), 'exact_positive_child_ray_births': child_births})
    return {'counts': dict(counts), 'fixtures': fixtures}


@lru_cache(None)
def downsets(n):
    if n == 0:
        return (frozenset(),)
    if n == 1:
        return (frozenset({(0, 0, 0)}),)
    out = set()
    for B in downsets(n-1):
        for a in minimal_excluded_border(tuple(B)):
            out.add(B | {a})
    return tuple(sorted(out, key=lambda B: sorted(B)))


def group_fan(order, residue, norm, counts):
    n = order-int(norm)
    out = []
    for B in downsets(n):
        classes = {residue(a): a for a in B}
        if len(classes) != n:
            continue
        normals = set()
        for a in minimal_excluded_border(tuple(sorted(B))):
            r = residue(a)
            support = (classes[r],) if r in classes else tuple(B)
            for b in support:
                normals.add(tuple(x-y for x, y in zip(a, b)))
            raw = [0]*order
            raw[r] += 1
            if r in classes:
                raw[r] -= 1
            else:
                for b in B:
                    raw[residue(b)] += 1
            assert raw == ([1]*order if r not in classes else [0]*order)
            assert norm or not any(raw)
            counts['integral_border_polynomial_annihilations'] += 1
        poly = intersect(TRIANGLE, sorted(normals))
        if poly:
            out.append({'standard_set': B, 'polygon': poly, 'normals': sorted(normals),
                        'profile': tuple(sum(a[i] for a in B) for i in range(3))})
    assert sum(area2(c['polygon']) for c in out) == 1
    return out


def graph_distances(order, residue, w):
    steps = [residue(tuple(int(i == j) for i in range(3))) for j in range(3)]
    # Fixtures use a supplied group-addition table, attached to residue.
    distances = {residue((0, 0, 0)): Q(0)}
    pending = [(Q(0), residue((0, 0, 0)))]
    while pending:
        value, g = heappop(pending)
        if distances[g] != value:
            continue
        for j, step in enumerate(steps):
            target = residue.add(g, step)
            candidate = value+w[j]
            if candidate < distances.get(target, candidate+1):
                distances[target] = candidate
                heappush(pending, (candidate, target))
    assert len(distances) == order
    return distances


def group_controls():
    counts = defaultdict(int)
    fixtures = []
    specs = [(4, (1, 3, 2)), (5, (1, 2, 4)), (7, (1, 2, 5)),
             (8, (1, 3, 5))]
    models = []
    for modulus, steps in specs:
        def residue(a, m=modulus, v=steps):
            return dot(a, v) % m
        residue.add = lambda a, b, m=modulus: (a+b) % m
        models.append((f'cyclic_{modulus}', modulus, steps, residue))
    def klein(a):
        return (a[0]+a[2]) % 2+2*((a[1]+a[2]) % 2)
    klein.add = lambda a, b: a ^ b
    models.append(('klein_four', 4, ((1, 0), (0, 1), (1, 1)), klein))
    for name, order, steps, residue in models:
        full = group_fan(order, residue, False, counts)
        reduced = group_fan(order, residue, True, counts)
        common = overlay(full, reduced, counts)
        parent_vertices = {v for cell in full for v in cell['polygon']}
        merger_witnesses = []
        for w in sorted({v for cell in reduced for v in cell['polygon']}-parent_vertices):
            if min(w) > 0:
                distances = graph_distances(order, residue, w)
                last = max(distances.values())
                tied = sorted(g for g, value in distances.items() if value == last)
                assert len(tied) >= 2
                counts['actual_group_child_ray_merger_equalities'] += 1
                merger_witnesses.append({'weight': list(map(str, w)),
                                        'pivot_layer_residue_components': tied})
        gap_cells = []
        for cell in full:
            B = cell['standard_set']
            for farthest in sorted(B):
                extra = [tuple(a-b for a, b in zip(farthest, v)) for v in B if v != farthest]
                poly = intersect(cell['polygon'], extra)
                if poly:
                    target = B-{farthest}
                    match = [c for c in reduced if c['standard_set'] == target]
                    assert len(match) == 1
                    paired = intersect(cell['polygon'], match[0]['normals'])
                    assert set(paired) == set(poly)
                    gap_cells.append((cell, target, farthest, poly))
        assert len(gap_cells) == len(common)
        for cell, target, farthest, poly in gap_cells:
            w = tuple(sum(v[i] for v in poly)/len(poly) for i in range(3))
            distances = graph_distances(order, residue, w)
            assert {residue(a): dot(w, a) for a in cell['standard_set']} == distances
            gap = max(distances.values())
            assert gap == dot(w, farthest)
            reduced_profile = tuple(sum(a[i] for a in target) for i in range(3))
            assert dot(w, cell['profile'])-dot(w, reduced_profile) == gap
            counts['independent_graph_gap_and_support_differences'] += 1
        for missing in range(order):
            columns = [tuple(int(i == j) for i in range(order-1))
                       if j < order-1 else (-1,)*(order-1)
                       for j in range(order) if j != missing]
            assert abs(determinant(list(zip(*columns)))) == 1
            counts['integral_norm_quotient_basis_determinants'] += 1
        witnesses = [{'full_standard_set': sorted(c['standard_set']),
                      'norm_standard_set': sorted(B), 'farthest': a,
                      'weight_simplex_polygon': [[str(x) for x in v] for v in poly]}
                     for c, B, a, poly in gap_cells]
        fixtures.append({'group': name, 'order': order, 'generator_residues': steps,
                         'group_initial_ideals': len(full), 'norm_initial_ideals': len(reduced),
                         'whole_gap_fan_cells': len(gap_cells),
                         'exact_child_ray_mergers': merger_witnesses,
                         'gap_cells': witnesses})
    return {'counts': dict(counts), 'fixtures': fixtures}


def tensor_controls():
    counts = defaultdict(int)
    for p, modulus in ((2, 4), (2, 8), (3, 9), (5, 5)):
        for n in range(1, modulus):
            h = modulus//p
            q, r = divmod(n, h)
            smaller = [frozenset(B) for B in combinations(range(h), r)
                       if rank([tuple(comb(a, j) % p for j in range(r)) for a in B], p) == r]
            constructed = set()
            for B in smaller:
                choices = [list(combinations((u+v*h for v in range(p)), q+int(u in B)))
                           for u in range(h)]
                for fibers in product(*choices):
                    constructed.add(frozenset(x for fiber in fibers for x in fiber))
            direct = {frozenset(B) for B in combinations(range(modulus), n)
                      if rank([tuple(comb(a, j) % p for j in range(n)) for a in B], p) == n}
            assert constructed == direct
            counts['complete_uniform_tensor_basis_bijections'] += 1
            counts['independent_pascal_bases'] += len(direct)
            subset_ranks = [(len(S), max(len(S & B) for B in direct))
                            for bits in product((0, 1), repeat=modulus)
                            for S in [{i for i, bit in enumerate(bits) if bit}]]
            for x, y in ((2, 3), (3, 2), (0, 0)):
                tutte = sum((x-1)**(n-rank_value)*(y-1)**(size-rank_value)
                            for size, rank_value in subset_ranks)
                predicted = (x-1)**n*digit_rank_enumerator(
                    n, modulus, p, 1, y-1, Q(1, (x-1)*(y-1)))
                assert predicted == tutte
                counts['independent_tutte_rank_enumerator_conversions'] += 1
            # Independent subset restriction verifies rank via the tensor bases.
            for S in (set(range(modulus)), set(range(0, modulus, 2)), set(range(modulus//2))):
                assert max(len(S & B) for B in direct) == digit_rank(n, modulus, p, S)
                counts['tensor_rank_restriction_equalities'] += 1
    return dict(counts)


def socle_algebra(domain, support, columns, p, counts):
    """Construct actual quotient actions, then check them independently."""
    fixed = tuple(sorted(set(domain)-set(support)))
    size = len(fixed)+len(columns[0])
    raw = {a: tuple(int(i == j) for i in range(size))
           for j, a in enumerate(fixed)}
    raw.update({a: (0,)*len(fixed)+tuple(v) for a, v in zip(support, columns)})
    def evaluate(a):
        return {i: c for i, c in enumerate(raw.get(a, ())) if c}
    reference, pivots = field_basis(sorted(domain), evaluate, p)
    assert len(reference) == size
    def coordinates(a):
        rem, co = field_reduce(evaluate(a), pivots, p)
        assert not rem
        return tuple(co.get(b, 0) for b in reference)
    matrices = tuple(tuple(coordinates(tuple(x+(i == j) for j, x in enumerate(a)))
                           for a in reference) for i in range(3))
    one = coordinates((0, 0, 0))
    for i, j in combinations(range(3), 2):
        for k in range(size):
            v = tuple(int(l == k) for l in range(size))
            assert apply(matrices[i], apply(matrices[j], v, p), p) == apply(
                matrices[j], apply(matrices[i], v, p), p)
            counts['socle_commuting_multiplication_columns'] += 1
    for a in sorted(set(domain) | set(minimal_excluded_border(tuple(domain)))):
        assert image_from_matrices(matrices, one, a, p) == coordinates(a)
        counts['socle_independent_cyclic_columns'] += 1
    assert all(not any(coordinates(tuple(x+(i == j) for j, x in enumerate(a))))
               for a in support for i in range(3))
    return evaluate, matrices, coordinates, fixed


def deletion_control(domain, parent, child, removed, p, w, counts):
    B, pivots = field_basis(sorted(domain, key=lambda a: (dot(w, a), a)), parent, p)
    C, _ = field_basis(sorted(domain, key=lambda a: (dot(w, a), a)), child, p)
    rem, co = field_reduce(removed, pivots, p)
    assert not rem and co
    a = max(co, key=lambda b: (dot(w, b), b))
    assert set(C) == set(B)-{a}
    assert all(tuple(x+(i == j) for j, x in enumerate(a)) not in B for i in range(3))
    counts['independent_kernel_deletion_and_socle_equalities'] += 1
    return {'weight': list(map(str, w)), 'deleted': a,
            'parent_profile': tuple(sum(b[i] for b in B) for i in range(3)),
            'child_profile': tuple(sum(b[i] for b in C) for i in range(3)),
            'kernel_normal_form': [{'exponent': b, 'coefficient': str(c)}
                                   for b, c in sorted(co.items())]}


def ideal_flag_controls():
    counts = defaultdict(int)
    S = ((3, 0, 0), (0, 3, 0), (0, 0, 4), (2, 2, 0))
    D = tuple(sorted({a for v in S for a in product(*(range(x+1) for x in v))}))
    column_sets = (((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, -1, 1)),
                   ((1, 0), (1, 0), (0, 1), (0, 1)), ((1,),)*4)
    relations = (((1, -1, 1, -1),),
                 ((1, -1, 0, 0), (0, 0, 1, -1)),
                 ((1, -1, 0, 0), (1, 0, -1, 0), (1, 0, 0, -1)))
    assert len(D) == 15
    fixtures = []
    for p in (0, 2, 3, 7):
        models = []
        fans = []
        face_spans = []
        witness_cells = []
        for k, (columns, rows) in enumerate(zip(column_sets, relations)):
            assert rank(rows, p) == k+1
            assert all(not ((sum(c*v[i] for c, v in zip(row, columns))) % p if p
                            else sum(c*v[i] for c, v in zip(row, columns)))
                       for row in rows for i in range(3-k))
            if k:
                assert rank(list(relations[k-1])+list(rows), p) == len(rows)
            model = socle_algebra(D, S, columns, p, counts)
            models.append(model)
            ev, matrices, coordinates, fixed = model
            bases = [tuple(sorted(set(fixed) | {S[i] for i in B}))
                     for B in combinations(range(4), 3-k)
                     if rank([columns[i] for i in B], p) == 3-k]
            profile_sets = {tuple(sum(a[i] for a in B) for i in range(3)): B for B in bases}
            fan = point_fan(profile_sets)
            fans.append(fan)
            best = min(dot((1, 1, 1), v) for v in profile_sets)
            face_spans.append(differences(v for v in profile_sets if sum(v) == best))
            cells = []
            for cell in fan:
                w = tuple(sum(v[i] for v in cell['polygon'])/len(cell['polygon'])
                          for i in range(3))
                B, _ = field_basis(sorted(D, key=lambda a: (dot(w, a), a)), ev, p)
                assert tuple(sorted(B)) == profile_sets[cell['profile']]
                def ev_coordinates(a):
                    return {i: c for i, c in enumerate(coordinates(a)) if c}
                recovered = verify_reduced_basis(B, ev_coordinates, matrices, w, p, counts)
                record = {'standard_set': sorted(B), 'weight': list(map(str, w)),
                          'polygon': [[str(c) for c in v] for v in cell['polygon']],
                          'reduced_basis_polynomials': len(recovered)}
                if not p:
                    record['reduced_basis'] = recovered
                else:
                    record['reduced_basis_sha256'] = hashlib.sha256(
                        json.dumps(recovered, sort_keys=True).encode()).hexdigest()
                cells.append(record)
                counts['independent_actual_ideal_profile_and_border_cones'] += 1
            witness_cells.append(cells)
        assert [rank(s) for s in face_spans] == [1, 2, 1]
        assert rank(face_spans[0]+face_spans[2]) == 2
        assert [len(fan) for fan in fans] == [4, 4, 3]
        deletions = []
        for j, row in enumerate(((1, -1, 0, 0), (1, 0, -1, 0))):
            parent, child = models[j][0], models[j+1][0]
            removed = defaultdict(int)
            for c, a in zip(row, S):
                for i, value in parent(a).items():
                    removed[i] += c*value
            common = overlay(fans[j], fans[j+1], counts)
            for _, _, poly in common:
                w = tuple(sum(v[i] for v in poly)/len(poly) for i in range(3))
                deletions.append(deletion_control(D, parent, child, dict(removed), p, w, counts))
        outside = [poly for a, b in product(fans[0], fans[2])
                   if (poly := intersect(a['polygon'], b['normals']))]
        assert len(outside) == 9 and sum(area2(poly) for poly in outside) == 1
        crossing = (Q(1, 3),)*3
        assert crossing in {v for poly in outside for v in poly}
        assert crossing not in {v for f in (fans[0], fans[2]) for c in f for v in c['polygon']}
        counts['actual_nonadjacent_transverse_crossings'] += 1
        fixtures.append({'characteristic': p, 'lengths': [14, 13, 12],
                         'outside_face_dimensions': [1, 1], 'middle_face_dimension': 2,
                         'outside_sum_dimension': 2, 'positive_fan_counts': [4, 4, 3],
                         'outside_common_cells': 9, 'crossing_weight': [1, 1, 1],
                         'complete_fans': witness_cells, 'deletions': deletions})
    arity = []
    for L in range(1, 5):
        D = tuple(a for a in product(range(L+1), repeat=3) if sum(a) <= L)
        S = tuple(a for a in D if sum(a) == L)
        m = len(S)
        columns = tuple(tuple(int(i == j) for i in range(m-1))
                        if j < m-1 else (-1,)*(m-1) for j in range(m))
        assert len(D) == comb(L+3, 3) and m == comb(L+2, 2)
        for p in (0, 2, 3, 7):
            ev, _, _, fixed = socle_algebra(D, S, columns, p, counts)
            assert rank(columns, p) == m-1
            assert all(rank([v for j, v in enumerate(columns) if j != i], p) == m-1
                       for i in range(m))
            fan = point_fan(tuple(-a[i] for i in range(3)) for a in S)
            assert len(fan) == 3
            for cell in fan:
                w = tuple(sum(v[i] for v in cell['polygon'])/len(cell['polygon'])
                          for i in range(3))
                B, _ = field_basis(sorted(D, key=lambda a: (dot(w, a), a)), ev, p)
                deleted = max(S, key=lambda a: (dot(w, a), a))
                assert set(B) == set(D)-{deleted}
                counts['simplex_relation_three_cell_equalities'] += 1
            counts['simplex_full_component_mergers'] += 1
        arity.append({'total_degree': L, 'parent_length': len(D), 'child_length': len(D)-1,
                      'merged_components_at_equal_cost': m, 'child_initial_ideals': 3})
    # Nonzero coordinate eigenvalues: actual separators of three rational points.
    points = ((1, 1, 1), (2, 3, 4), (3, 2, 5))
    E = ((0, 0, 0),)+tuple(tuple(j if i == k else 0 for i in range(3))
                           for k in range(3) for j in (1, 2))
    for p in (0, 2, 3, 7):
        def evaluate(a):
            return {i: value for i, point in enumerate(points)
                    if (value := (point[0]**a[0]*point[1]**a[1]*point[2]**a[2]) % p
                        if p else point[0]**a[0]*point[1]**a[1]*point[2]**a[2])}
        for missing in range(3):
            def child(a):
                return {i: c for i, c in evaluate(a).items() if i != missing}
            for w in product(range(1, 4), repeat=3):
                deletion_control(E, evaluate, child, {missing: 1}, p, w, counts)
                counts['rational_point_separator_deletions_including_ties'] += 1
    return {'counts': dict(counts), 'actual_chain_support':
            ((3, 0, 0), (0, 3, 0), (0, 0, 4), (2, 2, 0)),
            'chain_relation_rows': relations, 'complete_chain_fixtures': fixtures,
            'unbounded_merger_arity': arity,
            'scope': 'Actual polynomial quotient actions and border identities, not arbitrary matroid data. Finite controls support the uniform socle and ideal-flag proofs; no unrestricted count follows.'}


def budget_comparison_controls():
    """Primitive edge budgets and a uniform fixed-fan core-length guard."""
    from run_canonical_fan_deflation_checks import (common_factor, evaluate_polynomial,
                                                   marked_basis)
    from run_coordinate_slice_budget_checks import edges_from_polygons, profile
    from run_initial_cone_height_checks import integer_ray
    from run_initial_fan_padding_checks import full_fan, generic_model
    from run_primitive_wall_checks import divisors, threshold
    from run_projection_tag_checks import fan, point_model, polygon_key
    counts=defaultdict(int)
    directions=[q for q in product(range(-3,4),repeat=3)
                if min(q)<0<max(q) and gcd(*q)==1 and next(x for x in q if x)>0]
    for q,r in combinations(directions,2):
        a,b=sum(map(abs,q)),sum(map(abs,r))
        c=cross(q,r)
        assert max(map(abs,c))<=(a-1)*(b-1)
        counts['independent_mixed_integer_cross_product_inequalities']+=1
        if min(c)>0 or max(c)<0:
            w=primitive(c)
            H=max(map(abs,w))
            assert 2*(isqrt(H-1)+1)<=a+b
            counts['positive_primitive_ray_rounding_inequalities']+=1

    def budgets(cells):
        edges=edges_from_polygons(cells,counts)
        directions={edge:primitive(tuple(x-y for x,y in zip(profile(B),profile(C))))
                    for edge in edges for B,C in [tuple(edge)]}
        W=sum(sum(map(abs,q)) for q in directions.values())
        V=sum(threshold(q) for q in directions.values())
        rays={integer_ray(v) for poly in cells.values() for v in poly if min(v)>0}
        U=0
        for w in rays:
            point=tuple(Q(x,sum(w)) for x in w)
            incident={B for B,poly in cells.items() if point in poly}
            local_edges=[edge for edge in edges if set(edge)<=incident]
            assert len(local_edges)==len(incident)>=3
            local_U=(isqrt(max(w)-1)+1)*len(incident)
            assert local_U<=sum(sum(map(abs,directions[e])) for e in local_edges)
            counts['whole_local_polygon_incidence_length_inequalities']+=1
            U+=local_U
        assert U<=2*W<=2*V
        counts['whole_fan_incidence_length_cost_inequalities']+=1
        return {'U':U,'W':W,'V':V,'interior_rays':len(rays),'walls':len(edges)}

    expected=[]
    for order in ((2,None),(1,(2,0)),(1,(0,2)),(0,(2,1)),(0,(1,2))):
        top,rest=order
        pairs=[(top,i) for i in range(3) if i!=top] if rest is None else [(top,rest[0]),rest]
        poly=list(TRIANGLE)
        for i,j in pairs:
            poly=clip(poly,tuple(int(k==i)-int(k==j) for k in range(3)),1)
        expected.append(polygon_key(poly))
    expected=set(expected)
    family=[]
    for k in range(2,13):
        @lru_cache(None)
        def ev(a):
            e,b,c=a
            if e>=2:return {}
            return {(e+j)*k+b+c-j:comb(c,j) for j in range(min(c,1-e)+1)
                    if b+c-j<k}
        cells=fan(2*k,ev,0,counts)
        assert len(cells)==5 and {polygon_key(P) for P in cells.values()}==expected
        assert all(common_factor(marked_basis(B,ev,0)[0])==(0,0,0) for B in cells)
        row=budgets(cells)
        assert row=={'U':5,'W':10,'V':10,'interior_rays':1,'walls':5}
        # Independent explicit standard sets and reduced head complements.
        standards=[{(a,b,0) for a in range(2) for b in range(k)},
                   {(a,0,b) for a in range(2) for b in range(k)},
                   {(0,0,b) for b in range(k+1)}|{(1,0,b) for b in range(k-1)},
                   {(0,a,b) for a in range(k) for b in range(2)},
                   {(0,0,b) for b in range(k+1)}|{(0,1,b) for b in range(k-1)}]
        assert set(cells)=={frozenset(B) for B in standards}
        counts['uniform_five_cone_explicit_standard_set_equalities']+=5
        @lru_cache(None)
        def translated(a):
            e,b,c=a
            if e>=2:return {}
            return {(e+i)*k+b+j:comb(c,i)*comb(c-i,j)
                    for i in range(min(c,1-e)+1)
                    for j in range(min(c-i,k-b-1)+1)}
        shifted=fan(2*k,translated,0,counts)
        assert set(shifted)==set(cells)
        assert all(polygon_key(shifted[B])==polygon_key(cells[B]) for B in cells)
        assert budgets(shifted)=={'U':5,'W':10,'V':10,'interior_rays':1,'walls':5}
        for B in shifted:
            G,_=marked_basis(B,translated,0)
            top=[{a:c for a,c in g.items() if sum(a)==max(map(sum,g))} for g in G]
            assert common_factor(top)==(0,0,0)
            assert all(not evaluate_polynomial(g,ev,0) for g in top)
            counts['translated_full_dimension_positive_initial_basis_checks']+=1
        ordinary=standards[0]-{(0,0,0)}|{(0,0,1)}
        assert len(field_basis(sorted(ordinary),translated,0)[0])==2*k
        graded_profiles=[profile(B) for B in cells]
        assert len({sum(a) for a in graded_profiles})==1
        assert sum(profile(ordinary))==sum(graded_profiles[0])+1
        assert rank([tuple(x-y for x,y in zip(a,graded_profiles[0]))
                     for a in graded_profiles]+[tuple(x-y for x,y in
                         zip(profile(ordinary),graded_profiles[0]))])==3
        row.update(k=k,length=2*k,local_factor=(0,0,0),local_core_length=2*k,M=10*k,
                   translated_full_profile_dimension=3,
                   translated_full_fan_equal=True,
                   standard_sets=[sorted(B) for B in standards],
                   common_full_fan_polygons=[[[str(c) for c in v] for v in P]
                                            for P in sorted(expected)])
        family.append(row)

    radicals=[]
    for n in range(3,8):
        points=tuple((j,j*j+j,j*j*j+2*j) for j in range(1,n+1))
        _,ev=point_model(points,0)
        cells=fan(n,ev,0,counts)
        radicals.append(dict(budgets(cells),length=n,points=points,initial_ideals=len(cells)))
    relaxation=[]
    for k,p in product((2,3,5),(0,2)):
        u,v=(k,k,0),(0,0,k+1)
        D=set(divisors(u))|set(divisors(v))
        model=generic_model(D,{u,v},[{u:1,v:-1}],p)
        cells=full_fan(model,counts);row=budgets(cells)
        assert len(cells)==2 and row['W']==3*k+1 and row['V']==k*k+3*k+1
        assert model['n']==row['V']
        row.update(k=k,characteristic=p,length=model['n'],primitive_direction=(k,k,-k-1))
        relaxation.append(row)
    return {'counts':dict(counts),'fixed_five_cone_core_inflation':family,
            'radical_complete_fan_controls':radicals,'single_wall_length_cost_separation':relaxation,
            'scope':'Uniform elementary comparison U<=2W<=2V, and a verified infinite proof-route guard against replacing U or W by exact local core length. No cubic aggregate or unrestricted exponent follows.'}


def main():
    if not __debug__:
        raise RuntimeError('Exact assertions must be enabled.')
    matroids = matroid_controls()
    quotients = elementary_quotient_controls()
    groups = group_controls()
    tensor = tensor_controls()
    flags = ideal_flag_controls()
    budgets = budget_comparison_controls()
    sources = ('docs/adjacent_rank_fans.md', 'scripts/run_adjacent_rank_checks.py',
               'docs/rational_shadow_bound.md', 'docs/profile_lattice_bound.md',
               'docs/interleaving_research.md',
               'docs/oracle_fan_enumeration.md', 'scripts/run_oracle_enumeration_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_signed_face_checks.py',
               'docs/canonical_fan_deflation.md','scripts/run_canonical_fan_deflation_checks.py',
               'scripts/run_coordinate_slice_budget_checks.py','scripts/run_initial_cone_height_checks.py',
               'scripts/run_initial_fan_padding_checks.py','scripts/run_primitive_wall_checks.py',
               'scripts/run_projection_tag_checks.py')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_adjacent_rank_checks.py',
              'matroid_adjacent_faces': matroids, 'finite_group_gap_fans': groups,
              'elementary_linear_quotients': quotients,
              'prime_digit_tensor_counts': tensor,
              'actual_ideal_flags': flags,
              'incidence_primitive_length_and_local_core_budgets': budgets,
              'uniform_results': {
                  'adjacent_exposed_direction_spaces': 'Nested for every cost and every linear projection; inclusion reverses according as the last greedy rank increment is one or at least two.',
                  'elementary_quotient_exposed_direction_spaces': 'Nested for every actual one-dimensional linear quotient. Only the cost layer where the removed line first enters the prefix span changes; single-component support refines that component, multiple-component support merges components.',
                  'elementary_quotient_child_ray_condition': 'A new child ray requires the removed line to have support in multiple connected components of the pivot cost minor; it is a necessary condition, not a count.',
                  'positive_three_coordinate_overlay': 'max(a,b)<=c<=4*(a+b)+2, including flat and tied profile configurations',
                  'new_ray_necessary_condition': 'A new rank-r interior ray absent at rank r-1 must have last-layer rank increment one.',
                  'threshold_scope': 'The exact last greedy cost is linear on this O(a+b) refinement. The inherited O(k^(23/8)*log(k+1)^(3/8)) coarse power also follows from the previous shadow theorem on a direct sum; it is not a new unrestricted exponent.',
                  'norm_quotient_gap_fan': 'For every full-rank lattice of index m in Z^3, its faithful gap fan is the common refinement of the group ideal and its all-fields length-(m-1) norm quotient.',
                  'corank_one_boundary': 'The earlier binomial-extension envelope theorem already excludes polynomial amplification in this whole family and directly bounds the faithful gap diagram.',
                  'tensor_priority_boundary': 'The Pascal basis recursion is an iterated pointed-uniform two-sum substitution. The general Tutte substitution is classical; the exact code-specific priority remains unassessed.',
                  'actual_kernel_socle_rule': 'After a rank-preserving field extension, every finite algebra has an ideal composition series. A positive-filter initial kernel line is socle; each ordered child deletes the leading supported monomial of the parent kernel normal form, necessarily a staircase corner.',
                  'unbounded_merger_arity': 'A single socle relation in the degree-L simplex merges binom(L+2,2) singleton components at equal weight. The child still has exactly three initial ideals over every field.',
                  'actual_nonadjacent_crossing': 'The explicit all-fields chain of lengths 14,13,12 has face dimensions 1,2,1 at weight (1,1,1), with nonparallel outside edges. Its outside fan overlay has a new interior vertex and nine cells; adjacent nesting does not iterate.',
                  'local_realization_boundary': 'Arbitrary nested relation spaces on a downset maximal-corner support give actual ideal flags. They are already within the full-socle O(n^(3/2)) count, so a new global argument needs a length or multiplication-layer incidence budget.',
                  'unrestricted_incidence_obligation': 'A cubic polylogarithmic bound for U=sum_w ceil(sqrt(height(w)))*incident_cells(w) would give G3=O(n^(33/13)*log(n+1)^((10*C+3)/13)). This budget is unproved. The general conditional reduction is established; no unrestricted improvement follows yet.'},
              'new_unrestricted_G3_count': False, 'polynomial_lower_separation': False,
              'novelty_established': False, 'goal_complete': False,
              'scope': 'Exact finite controls for the uniform proofs. Exhaustive small projected faces and complete positive fans, actual integral group/norm quotients, shortest-path gap checks, uniform-gadget basis bijections and actual socle ideal flags. The new matrix-action and border checks independently verify nonadjacent crossings, unbounded merger arity and separator deletions. No growth fitting, generic-rank replacement, runtime benchmark or priority conclusion.',
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/'results/research/adjacent_rank_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'matroid': matroids['counts'], 'quotients': quotients['counts'],
                      'groups': groups['counts'], 'tensor': tensor,
                      'actual_ideal_flags': flags['counts'],
                      'incidence_budget_comparisons':budgets['counts']}, indent=2))


if __name__ == '__main__':
    main()
