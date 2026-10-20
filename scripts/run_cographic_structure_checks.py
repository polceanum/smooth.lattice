#!/usr/bin/env python3
"""Exact checks for ordered cographic bases and sparse periodic tree costs.

Uniform proofs and primary-source boundaries: docs/cographic_order_ideals.md.
No all-size fan count, novelty certificate or unrestricted separation follows.
"""
from fractions import Fraction as Q
from bisect import bisect_left
from itertools import combinations, combinations_with_replacement, product
from math import ceil, comb, floor, gcd, prod
from pathlib import Path
import hashlib
import json

from run_bilinear_cographic_checks import components, independent_tree_certificate, maximum_tree
from run_bilinear_rectangle_checks import labels, rank_sequence
from run_extremal_initial_ideals_checks import rectangle_fan
from run_initial_cone_height_checks import add, leading, normal_form, subtract_shift
from run_interleaving_research import minimal_excluded_border
from run_rectangle_boundary_checks import predicted_basis, slots
from run_rectangle_rank_checks import rectangle
from run_structural_complexity_checks import clip, cross, determinant, dot, greedy, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]


def require(value, reason):
    if not value:
        raise AssertionError(reason)


def representative(q, edge, weight):
    u, v = edge
    if weight[2] > weight[0]+weight[1]:
        return (u, v, 0)
    d = (v-u) % q
    options = [((-d) % q, 0, v)]
    if d:
        options.append((0, d, (v-d) % q))
    values = [dot(weight, a) for a in options]
    require(values.count(min(values)) == 1, "strict least representative")
    return options[values.index(min(values))]


def affine_cost(q, edge, weight):
    u, v = edge
    alpha, beta, gamma = weight
    if gamma > alpha+beta:
        return alpha*u+beta*v, (alpha, beta)
    first = alpha*u+(gamma-alpha)*v+(q*alpha if u < v else 0)
    second = (gamma-beta)*u+beta*v+(q*beta if u >= v else 0)
    return (first, (alpha, gamma-alpha)) if first <= second else (second, (gamma-beta, beta))


def candidate_edges(q, weight):
    alpha, beta, gamma = weight
    delta = alpha+beta-gamma
    require(delta and gamma not in (alpha, beta), "generic parameter planes avoided")
    result = {(u, v) for u in range(q) for v in (0, q-1)}
    result |= {(u, v) for u in (0, q-1) for v in range(q)}
    if delta > 0:
        for level in (Q(0), Q(q*beta, delta), Q(-q*alpha, delta)):
            for difference in range(ceil(level-2), floor(level+2)+1):
                for u in range(q):
                    v = u-difference
                    if 0 <= v < q:
                        result.add((u, v))
    require(len(result) <= 19*q, "linear candidate-edge bound")
    return result


def compressed(q, weight):
    candidates = candidate_edges(q, weight)
    costs = {e: dot(weight, representative(q, e, weight)) for e in candidates}
    parent, size = list(range(2*q)), [1]*(2*q)
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    tree, remaining = [], 2*q
    for u, v in sorted(candidates, key=lambda e: costs[e], reverse=True):
        a, b = find(u), find(q+v)
        if a != b:
            if size[a] > size[b]:
                a, b = b, a
            parent[a] = b
            size[b] += size[a]
            tree.append((u, v))
            remaining -= 1
    require(remaining == 1, "sparse candidate graph is connected")
    return frozenset(tree), candidates


def compressed_membership(q, ordered_tree, weight, exponent):
    edge = labels(exponent, q, q)
    index = bisect_left(ordered_tree, edge)
    omitted = index < len(ordered_tree) and ordered_tree[index] == edge
    return not omitted and representative(q, edge, weight) == exponent


def weights(q):
    base = 2*q+1
    scale = (2*q+1)*(1+base+base**2)+1
    return [tuple(scale*x+y for x, y in zip(w, (1, base, base**2)))
            for w in ((1, 1, 3), (3, 4, 6), (13, 10, 5), (9, 7, 3), (10, 13, 5), (3, 10, 5))]


def sparse_controls(q):
    records = []
    edges = set(product(range(q), repeat=2))
    for weight in weights(q):
        tree, candidates = compressed(q, weight)
        # The reference deliberately constructs every edge independently.
        dense = {e: dot(weight, representative(q, e, weight)) for e in edges}
        require(len(set(dense.values())) == q*q, "strict edge costs")
        require(tree == maximum_tree(q, q, dense), "sparse versus full graph optimization")
        cycle_checks = independent_tree_certificate(q, q, tree, dense)
        for u, v in edges-candidates:
            cost, gradient = affine_cost(q, (u, v), weight)
            du, dv = (1 if x > 0 else -1 for x in gradient)
            rectangle_edges = ((u+du, v), (u, v+dv), (u+du, v+dv))
            require(cost == dense[u, v] and all(0 <= a < q and 0 <= b < q and dense[a, b] > cost
                    for a, b in rectangle_edges), "explicit strict four-cycle domination")
        records.append({'strict_weight': weight, 'candidate_edges': len(candidates),
                        'dense_edges': q*q, 'excluded_edge_cycle_witnesses': q*q-len(candidates),
                        'omitted_tree': sorted(tree), 'independent_tree_cycles': cycle_checks})
    return {'period': q, 'length': (q-1)**2, 'witnesses': records}


def polynomial_residue(q, characteristic):
    # Ordinary multiplication by X followed by division by S_q(X).
    # No periodic label reduction or graph rank is used in this oracle.
    powers = [(1,)+(0,)*(q-2)]
    def power(k):
        while len(powers) <= k:
            old = powers[-1]
            new = (-old[-1],)+tuple(old[i-1]-old[-1] for i in range(1, q-1))
            powers.append(tuple(x % characteristic for x in new) if characteristic else new)
        return powers[k]
    def residue(a):
        left, right = power(a[0]+a[2]), power(a[1]+a[2])
        return tuple(x*y % characteristic if characteristic else x*y for x in left for y in right)
    return residue


def branch_check(q, basis):
    selected = {labels(a, q, q) for a in basis}
    require(len(selected) == len(basis), "basis has no parallel duplicate label")
    tree = set(product(range(q), repeat=2))-selected
    require(len(tree) == 2*q-1 and components(q, q, tree) == 1, "complementary spanning tree")
    h = sum((0, 0, c) in basis for c in range(q))
    require(h in (q-1, q), "only two possible axis lengths")
    row_neighbors = [{v for a, v in tree if a == u} for u in range(q)]
    column_neighbors = [{u for u, b in tree if b == v} for v in range(q)]
    for axis, neighbors in enumerate((row_neighbors, column_neighbors)):
        for u in range(q-1):
            b = max((j for j in range(q-1) if ((0, j, u) if axis == 0 else (j, 0, u)) in basis), default=-1)
            lost = neighbors[u]-neighbors[u+1]
            require(lost <= {(u+b+1) % q} and len(neighbors[u] & neighbors[u+1]) <= 1,
                    "independent neighbor-inclusion and four-cycle test")
            require(len(neighbors[u]) <= 2, "branching confined to the two final vertices")
    return h


def exhaustive_downsets(q):
    n = (q-1)**2
    residue = polynomial_residue(q, 0)
    found, digest = [], hashlib.sha256()
    for h in range(1, q+1):
        parts = [tuple(reversed(a)) for a in combinations_with_replacement(range(h+1), q-2)]
        by_sum = {}
        for b in parts:
            by_sum.setdefault(sum(b), []).append(b)
        for a in parts:
            for b in by_sum.get(n-h-sum(a), ()):
                basis = frozenset([(0, 0, c) for c in range(h)]+
                        [(i, 0, c) for i, count in enumerate(a, 1) for c in range(count)]+
                        [(0, j, c) for j, count in enumerate(b, 1) for c in range(count)])
                ranks = rank_sequence([residue(u) for u in sorted(basis)], 0)
                if ranks[-1] != n:
                    continue
                branch_check(q, basis)
                record = {'axis_length': h, 'x_arm_heights': a, 'y_arm_heights': b}
                found.append(record)
                digest.update(json.dumps(record, sort_keys=True).encode())
    return {'period': q, 'all_downset_bases_checked': len(found),
            'axis_lengths': sorted({r['axis_length'] for r in found}),
            'exact_height_witnesses': found, 'enumeration_sha256': digest.hexdigest(),
            'scope': 'All downsets inside D_q with n elements, selected by independent ordinary polynomial rank; coherence is not assumed or counted.'}


def geometric_controls(q):
    n = (q-1)**2
    ground = hyperbolic_cross(n, 3)
    records = []
    for field in (0, 2, 3):
        residue = polynomial_residue(q, field)
        columns = {a: residue(a) for a in ground}
        for weight in weights(q):
            tree, _ = compressed(q, weight)
            basis = frozenset(representative(q, e, weight) for e in set(product(range(q), repeat=2))-tree)
            require(set(basis) == greedy(ground, columns, weight, field),
                    "compressed graph output versus ordinary polynomial quotient")
            ordered_tree = sorted(tree)
            require(all(compressed_membership(q, ordered_tree, weight, a) == (a in basis) for a in ground),
                    "compressed membership versus independent polynomial standard sets")
            if weight[2] < weight[0]+weight[1]:
                branch_check(q, basis)
            records.append({'characteristic': field, 'strict_weight': weight, 'standard_set': sorted(basis)})
    return {'period': q, 'independent_polynomial_standard_sets': records}


def border_cone(q, basis, weight, replay=False):
    basis = frozenset(basis)
    by_label = {labels(a, q, q): a for a in basis}
    tree = set(product(range(q), repeat=2))-set(by_label)
    normals, polynomials = set(), []
    for head in minimal_excluded_border(sorted(basis)):
        edge = labels(head, q, q)
        polynomial = {head: 1}
        if edge in by_label:
            polynomial[by_label[edge]] = -1
        else:
            adjacency = [[] for _ in range(2*q)]
            for u, v in tree-{edge}:
                adjacency[u].append(q+v)
                adjacency[q+v].append(u)
            part, stack = {edge[0]}, [edge[0]]
            while stack:
                for other in adjacency[stack.pop()]:
                    if other not in part:
                        part.add(other)
                        stack.append(other)
            sign = int(edge[0] in part)-int(q+edge[1] in part)
            require(abs(sign) == 1, "fundamental cut separates the deleted edge")
            for (u, v), a in by_label.items():
                coefficient = int(u in part)-int(q+v in part)
                if coefficient:
                    polynomial[a] = coefficient*sign
        for a in polynomial:
            if a != head:
                z = tuple(x-y for x, y in zip(head, a))
                factor = gcd(*z)
                normals.add(tuple(x//factor for x in z))
        polynomials.append({a: Q(c) for a, c in polynomial.items()})
    poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
    for normal in sorted(normals):
        require(dot(weight, normal) > 0, "strict ordinary border heads")
        poly = clip(poly, normal, 1)
    require(poly, "a full positive cone")
    checks = 0
    if replay:
        for field in (0, 2, 3):
            residue = polynomial_residue(q, field)
            for polynomial in polynomials:
                values = [sum(c*residue(a)[i] for a, c in polynomial.items()) for i in range((q-1)**2)]
                require(all((x % field if field else x) == 0 for x in values),
                        "fundamental cut border relation versus ordinary polynomial remainders")
                checks += 1
    pairs = 0
    if q <= 4:
        for g, h in combinations(polynomials, 2):
            a, b = leading(g, weight), leading(h, weight)
            common = tuple(max(x, y) for x, y in zip(a, b))
            first = {add(c, tuple(x-y for x, y in zip(common, a))): Q(value) for c, value in g.items()}
            spair = subtract_shift(first, h, tuple(x-y for x, y in zip(common, b)), Q(1), 0)
            require(not normal_form(spair, polynomials, weight, 0), "ordinary rational S-pair")
            pairs += 1
    return normals, checks, pairs


def coherent_wall(q):
    a, r = q//2, q-1
    require(q >= 4 and q & (q-1) == 0, "power of two")
    rho = Q(a-1, a)
    epsilon = Q(1, 10*a*a)
    normal, wall = (a, -(a-1), -a), (r*a+a-1, a, r*a)
    sides, border_checks, pairs = [], 0, 0
    for value in (rho-epsilon, rho+epsilon):
        weight = (r*value.denominator+value.numerator, value.denominator, r*value.denominator)
        tree, _ = compressed(q, weight)
        basis = frozenset(representative(q, e, weight) for e in set(product(range(q), repeat=2))-tree)
        require(basis == predicted_basis(r, r, 1, value), "all-field bilinear outer selection")
        normals, replay, count = border_cone(q, basis, weight, replay=q <= 8)
        zero = [h for h in normals if not dot(h, wall)]
        require(all(dot(h, wall) >= 0 for h in normals) and zero and
                all(cross(h, normal) == (0, 0, 0) for h in zero), "relative interior of an actual shared facet")
        sides.append({'weight': weight, 'basis': sorted(basis), 'active_cut_normals': sorted(zero)})
        border_checks += replay
        pairs += count
    before, after = (set(map(tuple, s['basis'])) for s in sides)
    require(before-after == {(a, 0, t-a) for t in range(a, r)} and
            after-before == {(0, a-1, t) for t in range(a, r)}, "exact parallel coherent swaps")
    return {'period': q, 'length': r*r, 'primitive_direction': normal,
            'direction_product': prod(abs(x)+1 for x in normal), 'shared_wall_weight': wall,
            'sides': sides, 'ordinary_border_replays': border_checks, 'ordinary_S_pairs': pairs}


def binomial_prefix_controls():
    comparisons = determinants = prefixes = 0
    for r, s in ((3, 3), (4, 2), (7, 7)):
        for field in (2, 3):
            model = rectangle(r, s, {(1, 0): 1, (0, 1): 1, (1, 1): 1}, field)
            for rho in (Q(2, 7), Q(5, 7), Q(13, 7), Q(-2, 7), Q(-9, 7)):
                weight = (s*rho.denominator+max(rho.numerator, 0), rho.denominator,
                          s*rho.denominator-min(rho.numerator, 0))
                basis = predicted_basis(r, s, 1, rho)
                require(basis == greedy(model['ground'], model['columns'], weight, field),
                        "outer bilinear formula in small characteristic")
                comparisons += 1
                for t in range(r):
                    order = sorted(slots(t, s, 1), key=lambda u: (rho*u[0]+u[1], u))
                    columns = [tuple(model['residue']((i, j, t-i))[model['index'][t, k]] for k in range(s))
                               for i, j in order]
                    require(rank_sequence(columns, field) == tuple(min(s, k+1) for k in range(len(order))),
                            "every direct leading prefix retains full rank")
                    prefixes += len(order)
                    for length in range(1, min(s, len(order))+1):
                        nodes = sorted(-i if i else j for i, j in order[:length])
                        require(nodes == list(range(nodes[0], nodes[0]+length)), "consecutive integer exponents")
                        def coefficient(m, k):
                            return comb(m, k) if m >= 0 else (-1)**k*comb(k-m-1, k)
                        require(determinant([[coefficient(m, k) for m in nodes] for k in range(length)]) == 1,
                                "integer binomial determinant, including length above characteristic")
                        determinants += 1
    return {'direct_small_characteristic_standard_sets': comparisons,
            'direct_leading_prefix_ranks': prefixes, 'exact_unit_binomial_determinants': determinants}


def triple_period_box_controls():
    records = []
    for q in (3, 5):
        ground = tuple(product(range(q), repeat=3))
        inverse_two = pow(2, -1, q)
        def label(a):
            return ((a[0]+a[1]) % q, (a[0]+a[2]) % q, (a[1]+a[2]) % q)
        require(len({label(a) for a in ground}) == q**3, "pair-product period map is a bijection")
        coordinate_basis = [tuple(inverse_two*x % q for x in (u+v-w, u+w-v, v+w-u))
                            for u, v, w in product(range(q-1), repeat=3)]
        for field in (0, 2, 3):
            def vector(a):
                factors = [tuple(int(i == x) for i in range(q-1)) if x < q-1 else (-1,)*(q-1)
                           for x in label(a)]
                return tuple(x*y*z % field if field else x*y*z for x, y, z in product(*factors))
            columns = {a: vector(a) for a in ground}
            rank = rank_sequence([columns[a] for a in coordinate_basis], field)[-1]
            require(rank == (q-1)**3, "explicit tensor coordinate basis has the stated quotient rank")
            # All translates of the three geometric-sum ideal generators
            # vanish. This checks actual multiplication compatibility.
            relations = 0
            relation_vectors = {}
            indices = {a: i for i, a in enumerate(ground)}
            for offset in ground:
                for step in ((1, 1, -1), (1, -1, 1), (-1, 1, 1)):
                    support = tuple(sorted(tuple((x+k*y) % q for x, y in zip(offset, step)) for k in range(q)))
                    translated = [columns[a] for a in support]
                    values = [sum(col[i] for col in translated) for i in range(rank)]
                    require(all((x % field if field else x) == 0 for x in values),
                            "every translated polynomial generator vanishes")
                    if support not in relation_vectors:
                        vector_relation = [0]*len(ground)
                        for a in support:
                            vector_relation[indices[a]] = 1
                        relation_vectors[support] = vector_relation
                    relations += 1
            relation_rank = rank_sequence(list(relation_vectors.values()), field)[-1]
            require(relation_rank == q**3-rank, "complete polynomial relation space has the right codimension")
            samples = []
            for weight in weights(q)[:2]:
                basis = greedy(ground, columns, weight, field)
                complement = set(ground)-basis
                defect = q**3-len(basis)
                profile = tuple(sum(a[i] for a in complement) for i in range(3))
                require(len(basis) == (q-1)**3 and all(x <= defect*(q-1) for x in profile),
                        "dual-profile coordinate budget for the classical corank screen")
                samples.append({'weight': weight, 'dual_rank': defect, 'dual_profile': profile,
                                'coordinate_budget': defect*(q-1)})
            records.append({'period': q, 'characteristic': field, 'box_size': q**3,
                            'quotient_length': rank, 'defect': q**3-rank,
                            'independent_relation_space_rank': relation_rank,
                            'translated_generators_checked': relations, 'profile_controls': samples})
    return records


def bounded_phase_references():
    records = []
    for q, p in ((3, 3), (4, 2), (5, 5), (8, 2)):
        # Independent original local multiplication and complete cone walk.
        reference = rectangle_fan(q-1, q-1, {(1, 0): 1, (0, 1): 1, (1, 1): 1}, p)
        witnesses = []
        for record in reference['witnesses']:
            basis = frozenset(map(tuple, record['basis']))
            tree = set(product(range(q), repeat=2))-{labels(a, q, q) for a in basis}
            weight = tuple(record['strict_weight'])
            reconstructed = frozenset(representative(q, e, weight) for e in set(product(range(q), repeat=2))-tree)
            require(reconstructed == basis, "compressed exact complete-phase witness")
            witnesses.append({'strict_weight': weight, 'omitted_tree': sorted(tree),
                              'cone_polygon': record['cone_polygon']})
        records.append({'period': q, 'matching_characteristic_reference': p,
                        'initial_ideal_count': reference['initial_ideal_count'],
                        'positive_fan_exhaustive': reference['positive_fan_exhaustive'],
                        'normalized_double_area': reference['normalized_double_area'],
                        'scope': 'A complete finite fan, with all-field transfer from the proved cographic representation. No asymptotic law is inferred.',
                        'compressed_witnesses': witnesses})
        print(json.dumps({'bounded_complete_fan_passed': q}), flush=True)
    return records


def main():
    sparse = []
    for q in (3, 4, 5, 8, 32, 64, 127):
        sparse.append(sparse_controls(q))
        print(json.dumps({'sparse_graph_fixture_passed': q}), flush=True)
    downsets = [exhaustive_downsets(q) for q in (3, 4, 5)]
    geometric = [geometric_controls(q) for q in (3, 4, 5)]
    walls = [coherent_wall(q) for q in (4, 8, 16)]
    prefixes = binomial_prefix_controls()
    triple = triple_period_box_controls()
    phases = bounded_phase_references()
    sources = ['docs/cographic_order_ideals.md', 'scripts/run_cographic_structure_checks.py',
               'scripts/run_bilinear_cographic_checks.py', 'scripts/run_bilinear_rectangle_checks.py',
               'scripts/run_rectangle_boundary_checks.py', 'scripts/run_rectangle_rank_checks.py',
               'scripts/run_extremal_initial_ideals_checks.py',
               'scripts/run_initial_cone_height_checks.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_nested_jet_checks.py']
    summary = {'sparse_graph_models': len(sparse), 'dense_sparse_comparisons': sum(len(r['witnesses']) for r in sparse),
               'explicit_excluded_edge_cycles': sum(w['excluded_edge_cycle_witnesses'] for r in sparse for w in r['witnesses']),
               'independent_tree_cycle_checks': sum(w['independent_tree_cycles'] for r in sparse for w in r['witnesses']),
               'exhaustive_downset_bases': sum(r['all_downset_bases_checked'] for r in downsets),
               'ordinary_polynomial_standard_sets': sum(len(r['independent_polynomial_standard_sets']) for r in geometric),
               'all_field_coherent_wall_fixtures': len(walls),
               'ordinary_border_replays': sum(r['ordinary_border_replays'] for r in walls),
               'ordinary_S_pairs': sum(r['ordinary_S_pairs'] for r in walls),
               'equal_period_triple_models': len(triple),
               'triple_polynomial_generator_translates': sum(r['translated_generators_checked'] for r in triple),
               'triple_dual_profile_controls': sum(len(r['profile_controls']) for r in triple), **prefixes}
    summary.update({'bounded_complete_fans': len(phases),
                    'bounded_complete_fan_states': sum(r['initial_ideal_count'] for r in phases)})
    report = {'schema_version': 1, 'command': 'python3 -B scripts/run_cographic_structure_checks.py',
              'scope': 'Finite independent controls for uniform ordered-tree, sparse-optimization and corank-screen proofs, with bounded complete phase references. No all-size fan asymptotic or unrestricted separation/novelty assertion.',
              'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources},
              'sparse_graph_controls': sparse, 'all_small_downset_bases': downsets,
              'ordinary_polynomial_controls': geometric, 'coherent_wall_controls': walls,
              'classical_corank_screen_controls': triple,
              'bounded_complete_phase_references': phases,
              'small_characteristic_prefix_controls': prefixes, 'summary': summary}
    (ROOT/'results/research/cographic_structure_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(summary, sort_keys=True))


if __name__ == '__main__':
    main()
