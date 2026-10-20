#!/usr/bin/env python3
"""Exact controls for the complete cographic model of (q-1)-side rectangles.

The uniform proof and comparison limits are in docs/bilinear_cographic_family.md.
This script does not count the complete parametric fan or establish novelty.
"""
from collections import Counter
from itertools import combinations, product
from math import comb
from pathlib import Path
from random import Random
import hashlib
import json

from run_bilinear_rectangle_checks import labels, rank_sequence, shortest_paths, two_candidates
from run_rectangle_rank_checks import rectangle
from run_structural_complexity_checks import dot, greedy

ROOT = Path(__file__).resolve().parents[1]


def require(value, reason):
    if not value:
        raise AssertionError(reason)


def components(q1, q2, edges):
    parent = list(range(q1+q2))
    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for u, v in edges:
        a, b = find(u), find(q1+v)
        parent[a] = b
    return len({find(i) for i in range(q1+q2)})


def maximum_tree(q1, q2, costs):
    parent = list(range(q1+q2))
    size = [1]*(q1+q2)
    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    selected = []
    for edge in sorted(costs, key=lambda e: costs[e], reverse=True):
        u, v = edge
        a, b = find(u), find(q1+v)
        if a != b:
            if size[a] > size[b]:
                a, b = b, a
            parent[a] = b
            size[b] += size[a]
            selected.append(edge)
    require(len(selected) == q1+q2-1 and components(q1, q2, selected) == 1, "a spanning tree")
    return frozenset(selected)


def independent_tree_certificate(q1, q2, tree, costs):
    adjacency = [[] for _ in range(q1+q2)]
    for u, v in tree:
        adjacency[u].append((q1+v, costs[u, v]))
        adjacency[q1+v].append((u, costs[u, v]))
    checked = 0
    for u, v in set(costs)-tree:
        stack, seen = [(u, None)], {u}
        found = False
        while stack:
            vertex, smallest = stack.pop()
            if vertex == q1+v:
                require(costs[u, v] < smallest, "strict independent maximum-tree cycle certificate")
                found = True
                break
            for other, cost in adjacency[vertex]:
                if other not in seen:
                    seen.add(other)
                    stack.append((other, cost if smallest is None else min(smallest, cost)))
        require(found, "every complementary edge has a tree path")
        checked += 1
    return checked


def fixture(q1, q2, p):
    for period in (q1, q2):
        power = period
        while power % p == 0:
            power //= p
        require(power == 1, "both periods are powers of the field characteristic")
    edges = tuple(product(range(q1), range(q2)))
    full = set(edges)
    columns = {(u, v): tuple(comb(u, i)*comb(v, j) % p
                            for i in range(q1-1) for j in range(q2-1)) for u, v in edges}
    n = (q1-1)*(q2-1)
    # Pascal columns independently give every defining row/column relation.
    groups = [[columns[u, fixed] for u in range(q1)] for fixed in range(q2)]
    groups += [[columns[fixed, v] for v in range(q2)] for fixed in range(q1)]
    for group in groups:
        require(all(sum(col[j] for col in group) % p == 0 for j in range(n)),
                "direct Pascal row/column relation")
    if (q1, q2) == (3, 3):
        masks = set(range(1 << len(edges)))
    else:
        masks = {(1 << k)-1 for k in range(len(edges)+1)}
        masks |= {((1 << len(edges))-1)^m for m in list(masks)}
        rng = Random(7100+100*q1+q2)
        masks |= {rng.getrandbits(len(edges)) for _ in range(32 if q2 >= 8 else 128)}
    ranks, replay = Counter(), []
    for mask in sorted(masks):
        selected = {edge for i, edge in enumerate(edges) if mask >> i & 1}
        sequence = rank_sequence([columns[edge] for edge in sorted(selected)], p)
        actual = sequence[-1] if sequence else 0
        predicted = len(selected)+1-components(q1, q2, full-selected)
        require(actual == predicted, "every tested subset: Pascal elimination versus complementary graph")
        ranks[actual] += 1
        replay.append([mask, actual])
    trees = None
    if q2 <= 4:
        trees = [frozenset(T) for T in combinations(edges, q1+q2-1) if components(q1, q2, T) == 1]
        require(len(trees) == q1**(q2-1)*q2**(q1-1), "all small complete-bipartite spanning trees enumerated")
    model = rectangle(q1-1, q2-1, {(1, 0): 1, (0, 1): 1, (1, 1): 1}, p)
    base = q2+1
    scale = (q2+1)*(1+base+base**2)+1
    weights = [tuple(scale*x+y for x, y in zip(w, (1, base, base**2)))
               for w in ((1, 1, 3), (3, 4, 6), (13, 10, 5), (9, 7, 3), (10, 13, 5), (3, 10, 5))]
    witnesses = []
    for weight in weights:
        representatives = two_candidates(q1, q2, weight)
        require(representatives == shortest_paths(q1, q2, weight), "two-candidate rule versus independent shortest paths")
        by_label = {labels(a, q1, q2): a for a in representatives}
        costs = {edge: dot(weight, by_label[edge]) for edge in edges}
        require(len(set(costs.values())) == len(edges), "strict residue costs")
        tree = maximum_tree(q1, q2, costs)
        cycles = independent_tree_certificate(q1, q2, tree, costs)
        if trees is not None:
            best = max(sum(costs[e] for e in T) for T in trees)
            require(sum(costs[e] for e in tree) == best, "maximum over every small spanning tree")
        basis = frozenset(by_label[edge] for edge in full-tree)
        require(basis == greedy(model['ground'], model['columns'], weight, p),
                "graph-selected initial set versus original local polynomial multiplication")
        require(len(basis) == n and all(tuple(x-(i == j) for i, x in enumerate(a)) in basis
                for a in basis for j in range(3) if a[j]), "full standard monomial downset")
        witnesses.append({'strict_weight': weight, 'maximum_tree': sorted(tree),
                          'complementary_standard_set': sorted(basis),
                          'maximum_tree_cost': sum(costs[e] for e in tree),
                          'independent_cycle_inequalities': cycles})
    return {'periods': [q1, q2], 'characteristic': p, 'length': n,
            'subset_rank_controls': len(masks), 'all_subsets_exhausted': (q1, q2) == (3, 3),
            'subset_generation': 'All masks for periods (3,3); prefixes, complements and Random(7100+100*q1+q2) masks otherwise.',
            'subset_rank_histogram': dict(sorted(ranks.items())),
            'subset_rank_replay_sha256': hashlib.sha256(json.dumps(replay).encode()).hexdigest(),
            'enumerated_spanning_trees': len(trees) if trees is not None else None,
            'witnesses': witnesses}


def main():
    records = []
    for q1, q2, p in ((3, 3, 3), (4, 4, 2), (5, 5, 5), (8, 8, 2), (3, 9, 3), (4, 8, 2)):
        records.append(fixture(q1, q2, p))
        print(json.dumps({'cographic_fixture_passed': [q1, q2, p]}), flush=True)
    # The prime-power hypothesis cannot simply be removed.
    q, p = 4, 3
    failed_row_sum = [sum(comb(u, i) for u in range(q)) % p for i in range(q-1)]
    require(any(failed_row_sum), "a non-power boundary does not have the asserted row relations")
    sources = ['docs/bilinear_cographic_family.md', 'scripts/run_bilinear_cographic_checks.py',
               'scripts/run_bilinear_rectangle_checks.py', 'scripts/run_rectangle_rank_checks.py',
               'scripts/run_structural_complexity_checks.py', 'scripts/run_nested_jet_checks.py']
    report = {'schema_version': 1, 'command': 'python3 -B scripts/run_bilinear_cographic_checks.py',
              'scope': 'Finite independent controls for the proved whole-matroid graph reduction, not a full fan count, asymptotic separation or novelty claim.',
              'source_sha256': {s: hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources},
              'fixtures': records,
              'non_power_boundary': {'q': q, 'characteristic': p, 'nonzero_row_sum': failed_row_sum},
              'summary': {'complete_cographic_models': len(records),
                          'direct_subset_rank_controls': sum(r['subset_rank_controls'] for r in records),
                          'spanning_tree_standard_sets': sum(len(r['witnesses']) for r in records),
                          'small_spanning_trees_enumerated': sum(r['enumerated_spanning_trees'] or 0 for r in records),
                          'independent_tree_cycle_certificates': sum(w['independent_cycle_inequalities']
                                                                   for r in records for w in r['witnesses'])}}
    (ROOT/'results/research/bilinear_cographic_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['summary'], sort_keys=True))


if __name__ == '__main__':
    main()
