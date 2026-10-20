#!/usr/bin/env python3
"""Exact controls for the chain-allocation reduction and its limited counts."""
from fractions import Fraction
from itertools import combinations
from math import gcd, prod
from pathlib import Path
import hashlib
import json

from run_nested_jet_checks import border_basis, independent, model, rank_columns, schedule_greedy
from run_structural_complexity_checks import area2, clip, dot

ROOT = Path(__file__).resolve().parents[1]


def chains(n, A, B):
    return tuple(sorted((A*b+B*c, b, c) for b in range(n) for c in range(n)
                        if A*b+B*c < n))


def compositions(n, k):
    if k == 1:
        yield (n,)
    else:
        for a in range(n+1):
            for tail in compositions(n-a, k-1):
                yield (a,) + tail


def feasible(lengths, ground):
    prefix = 0
    for length, (q, _, _) in zip(lengths, ground):
        if prefix < q:
            return False
        prefix += length
    return True


def standard_set(lengths, ground):
    return frozenset((a, b, c) for length, (_, b, c) in zip(lengths, ground)
                     for a in range(length))


def allocation(state, ground):
    return tuple(sum((b, c) == (v[1], v[2]) for v in state) for _, b, c in ground)


def cost(lengths, ground, weight):
    return sum(weight[0]*(length*(length-1)//2) + (b*weight[1]+c*weight[2])*length
               for length, (_, b, c) in zip(lengths, ground))


def transfers(lengths, ground):
    """Use a literal feasibility replay, independently of crossed-slack tests."""
    for i, amount in enumerate(lengths):
        if not amount:
            continue
        for j in range(len(ground)):
            if i == j:
                continue
            candidate = list(lengths)
            candidate[i] -= 1
            candidate[j] += 1
            if feasible(candidate, ground):
                raw = (lengths[j]-lengths[i]+1,
                       ground[j][1]-ground[i][1], ground[j][2]-ground[i][2])
                yield i, j, raw


def normals(lengths, ground):
    return {tuple(x//gcd(*raw) for x in raw) for _, _, raw in transfers(lengths, ground)}


def dynamic_allocation(n, ground, weight):
    """Finite prefix DP, independent of matroid scheduling and local exchanges."""
    states = {0: (Fraction(0), (), 1)}
    for q, b, c in ground:
        next_states = {}
        for used, (value, prefix, count) in states.items():
            if used < q:
                continue
            for length in range(n-used+1):
                total = used+length
                candidate = value + weight[0]*(length*(length-1)//2) + \
                    (b*weight[1]+c*weight[2])*length
                entry = (candidate, prefix+(length,), count)
                if total not in next_states or candidate < next_states[total][0]:
                    next_states[total] = entry
                elif candidate == next_states[total][0]:
                    old = next_states[total]
                    next_states[total] = (candidate, min(old[1], entry[1]), min(2, old[2]+count))
        states = next_states
    return states[n]


def exhaustive_controls():
    fixtures = []
    for n, A, B in ((4, 2, 3), (5, 2, 3), (6, 2, 3), (4, 1, 1), (9, 3, 3)):
        ground = chains(n, A, B)
        data = model(n, (1, A, B))
        columns = {(a, b, c): data['column']((a, b, c))
                   for _, b, c in ground for a in range(n)}
        bases, count, transfer_tests = [], 0, 0
        for lengths in compositions(n, len(ground)):
            state = standard_set(lengths, ground)
            actual = rank_columns([columns[u] for u in state]) == n
            criterion = independent([a+A*b+B*c for a, b, c in state], n)
            assert actual == criterion == feasible(lengths, ground)
            count += 1
            if actual:
                base_normals = normals(lengths, ground)
                bases.append((lengths, base_normals))
                for i, amount in enumerate(lengths):
                    if not amount:
                        continue
                    for j in range(len(ground)):
                        if i == j:
                            continue
                        candidate = list(lengths)
                        candidate[i] -= 1
                        candidate[j] += 1
                        q_i, q_j = ground[i][0], ground[j][0]
                        crossed_slack = all(sum(lengths[k] for k, (q, _, _) in enumerate(ground)
                                               if q < cut) > cut
                                            for cut in {q for q, _, _ in ground}
                                            if q_i < cut <= q_j)
                        assert feasible(candidate, ground) == crossed_slack
                        transfer_tests += 1
        weights = ((1, 1, 1), (1, 2, 3), (2, 3, 5), (11, 3, 18),
                   (1, A, B), (3, 2, 1), (4, 17, 5), (5, 7, 36))
        local_tests, optimum_controls = 0, []
        for raw_weight in weights:
            weight = tuple(map(Fraction, raw_weight))
            values = [cost(lengths, ground, weight) for lengths, _ in bases]
            optimum = min(values)
            minimizers = [lengths for (lengths, _), value in zip(bases, values) if value == optimum]
            dp_value, dp_lengths, dp_count = dynamic_allocation(n, ground, weight)
            assert dp_value == optimum and dp_lengths in minimizers
            assert (dp_count == 1) == (len(minimizers) == 1)
            for (lengths, base_normals), value in zip(bases, values):
                gaps = [dot(weight, normal) for normal in base_normals]
                assert all(x >= 0 for x in gaps) == (value == optimum)
                assert all(x > 0 for x in gaps) == (value == optimum and len(minimizers) == 1)
                local_tests += 1
            selected = schedule_greedy(data['ground'], weight, (1, A, B), n)
            assert allocation(selected, ground) in minimizers
            optimum_controls.append({'weight': raw_weight, 'minimum_cost': str(optimum),
                                     'minimizer_count': len(minimizers), 'dp_allocation': dp_lengths})
        fixtures.append({'length': n, 'valuations': [1, A, B], 'chains': ground,
                         'all_compositions_matrix_checked': count, 'feasible_allocations': len(bases),
                         'transfer_feasibility_checks': transfer_tests,
                         'local_optimality_checks': local_tests, 'weight_controls': optimum_controls})
    return fixtures


def existing_fan_controls():
    old = json.loads((ROOT/'results/research/nested_jet_checks.json').read_text())
    fixtures = []
    for fan in old['fan_controls']:
        n = fan['length']
        if fan['valuations'][0] != 1:
            continue
        _, A, B = fan['valuations']
        ground = chains(n, A, B)
        for witness in fan['witnesses']:
            state = frozenset(map(tuple, witness['standard_set']))
            lengths = allocation(state, ground)
            assert standard_set(lengths, ground) == state and feasible(lengths, ground)
            weight = witness['strict_weight']
            inequalities = normals(lengths, ground)
            assert all(dot(weight, normal) > 0 for normal in inequalities)
            polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
            for normal in sorted(inequalities):
                polygon = clip(polygon, normal, 1)
            expected = {tuple(map(Fraction, vertex)) for vertex in witness['closed_cell_vertices']}
            assert area2(polygon) > 0 and set(polygon) == expected
        fixtures.append({'length': n, 'valuations': [1, A, B],
                         'complete_initial_cones_reconstructed': len(fan['witnesses'])})
    return fixtures


def interval_probes(walls):
    if not walls:
        return [Fraction(1)]
    return [walls[0]/2] + [(a+b)/2 for a, b in zip(walls, walls[1:])] + [walls[-1]+1]


def line_controls():
    fixtures = []
    for M in (2, 3, 4):
        n = M*M
        ground = tuple((a, b, c) for b in range(M) for c in range(M-b)
                       for a in range(n-M*(b+c)))
        chain_ground = chains(n, M, M)
        data = model(n, (1, M, M)) if M <= 3 else None
        for delta in (Fraction(17, 101), Fraction(157, 101), Fraction(601, 101)):
            walls = set()
            for u, v in combinations(ground, 2):
                degree_difference = u[1]+u[2]-v[1]-v[2]
                if degree_difference:
                    crossing = (v[0]-u[0]+delta*(v[1]-u[1]))/degree_difference
                    if crossing > 0:
                        walls.add(crossing)
            witnesses, previous, previous_prefix = [], None, None
            probes = interval_probes(sorted(walls))
            for p in probes:
                weight = (Fraction(1), p+delta, p)
                state = schedule_greedy(ground, weight, (1, M, M), n)
                if state == previous:
                    continue
                prefix = tuple(sum(b+c < k for a, b, c in state) for k in range(1, M))
                assert all(x >= M*k for k, x in enumerate(prefix, 1))
                if previous_prefix is not None:
                    assert prefix != previous_prefix
                    assert all(x <= y for x, y in zip(previous_prefix, prefix))
                lengths = allocation(state, chain_ground)
                value, dp_lengths, count = dynamic_allocation(n, chain_ground, weight)
                assert count == 1 and dp_lengths == lengths
                assert value == cost(lengths, chain_ground, weight)
                assert all(dot(weight, normal) > 0 for normal in normals(lengths, chain_ground))
                if data:
                    assert rank_columns([data['column'](u) for u in state]) == n
                witnesses.append({'p': str(p), 'weight': list(map(str, weight)),
                                  'cumulative_prefix': prefix, 'allocation': lengths,
                                  'standard_set': sorted(state)})
                previous, previous_prefix = state, prefix
            bound = 1+n*(M-1)//2
            assert len(witnesses) <= bound
            fixtures.append({'M': M, 'length': n, 'delta': str(delta),
                             'full_comparison_intervals': len(probes),
                             'distinct_strict_states': len(witnesses), 'proved_bound': bound,
                             'independent_matrix_checks': len(witnesses) if data else 0,
                             'witnesses': witnesses})
    return fixtures


def strip(M, d, delta):
    points = sorted((Fraction(a)+delta*b, a, b)
                    for a in range(M+1) for b in range(d+1))
    selected = points[:M]
    lengths = tuple(sum(b == index for _, a, b in selected) for index in range(d+1))
    assert points[M-1][0] < points[M][0]
    return lengths, points[M-1][0], points[M][0]


def asymmetric_line_controls():
    fixtures = []
    for n, A, B in ((6, 2, 3), (8, 1, 3), (9, 2, 5)):
        ground = tuple((a, b, c) for a in range(n) for b in range(n) for c in range(n)
                       if a+A*b+B*c < n)
        chain_ground = chains(n, A, B)
        starts = sorted({q for q, _, _ in chain_ground})
        base = (Fraction(1), Fraction(17, 101), Fraction(31, 103))
        walls = set()
        for u, v in combinations(ground, 2):
            difference = A*(u[1]-v[1])+B*(u[2]-v[2])
            if difference:
                crossing = (dot(base, v)-dot(base, u))/difference
                if crossing > 0:
                    walls.add(crossing)
        previous = previous_prefix = None
        witnesses = []
        probes = interval_probes(sorted(walls))
        for tau in probes:
            weight = (Fraction(1), base[1]+A*tau, base[2]+B*tau)
            state = schedule_greedy(ground, weight, (1, A, B), n)
            if state == previous:
                continue
            prefix = tuple(sum(A*b+B*c < q for a, b, c in state) for q in starts[1:])
            if previous_prefix is not None:
                assert prefix != previous_prefix
                assert all(x <= y for x, y in zip(previous_prefix, prefix))
            lengths = allocation(state, chain_ground)
            _, dp_lengths, count = dynamic_allocation(n, chain_ground, weight)
            assert count == 1 and dp_lengths == lengths
            assert all(dot(weight, normal) > 0 for normal in normals(lengths, chain_ground))
            witnesses.append({'tau': str(tau), 'weight': list(map(str, weight)),
                              'cumulative_prefix': prefix, 'allocation': lengths,
                              'standard_set': sorted(state)})
            previous, previous_prefix = state, prefix
        bound = 1+sum(n-q for q in starts[1:])
        assert len(witnesses) <= bound
        fixtures.append({'length': n, 'valuations': [1, A, B], 'distinct_starts': starts,
                         'base_weight': list(map(str, base)), 'full_comparison_intervals': len(probes),
                         'distinct_strict_states': len(witnesses), 'proved_bound': bound,
                         'witnesses': witnesses})
    return fixtures


def tight_stratum_controls():
    fixtures = []
    for M in range(2, 13):
        candidates = sorted({Fraction(p, q) for p in range(1, M+1) for q in range(1, M)})
        probes = interval_probes(candidates)
        joint, planar, witnesses, actual_walls = set(), set(), [], set()
        previous = None
        chain_ground = chains(M*M, M, M)
        full_ground = tuple((a, b, c) for b in range(M) for c in range(M-b)
                            for a in range(M*M-M*(b+c))) if M <= 4 else None
        for index, delta in enumerate(probes):
            data = [strip(M, d, delta) for d in range(M)]
            signature = tuple(row[0] for row in data)
            planar.add(signature[-1])
            if signature == previous:
                continue
            if previous is not None:
                wall = candidates[index-1]
                assert (wall.numerator+1)*(wall.denominator+1) <= 2*M
                actual_walls.add(wall)
            previous = signature
            if signature in joint:
                raise AssertionError('A joint prefix tuple returned after leaving its interval')
            joint.add(signature)
            p = min((data[i][2]-data[j][1])/(j-i)
                    for i in range(M) for j in range(i+1, M))/2
            assert p > 0
            weight = (Fraction(1), p+delta, p)
            lengths = tuple(signature[b+c][b] for _, b, c in chain_ground)
            assert sum(lengths) == M*M and feasible(lengths, chain_ground)
            assert all(dot(weight, normal) > 0 for normal in normals(lengths, chain_ground))
            state = standard_set(lengths, chain_ground)
            if full_ground:
                assert schedule_greedy(full_ground, weight, (1, M, M), M*M) == state
                _, dp_lengths, count = dynamic_allocation(M*M, chain_ground, weight)
                assert count == 1 and dp_lengths == lengths
            witnesses.append({'delta': str(delta), 'p': str(p),
                              'weight': list(map(str, weight)), 'degree_group_allocations': signature})
        assert len(joint) == len(actual_walls)+1 >= len(planar)
        first = witnesses[0]
        lengths = tuple(first['degree_group_allocations'][b+c][b] for _, b, c in chain_ground)
        state = standard_set(lengths, chain_ground)
        assert Fraction(first['delta']) < Fraction(1, M-1)
        axes = tuple(sum(tuple(a if i == j else 0 for i in range(3)) in state
                         for a in range(M*M)) for j in range(3))
        assert axes == (M, M, M)
        envelope = sum(gcd(p, q) == 1 and (p+1)*(q+1) <= 2*M
                       for p in range(1, M) for q in range(1, M))
        assert len(actual_walls) <= envelope
        fixtures.append({'M': M, 'length': M*M, 'complete_ratio_intervals': len(probes),
                         'distinct_joint_states_wy_above_wz': len(joint),
                         'last_group_ordinary_planar_states': len(planar),
                         'exact_walls': list(map(str, sorted(actual_walls))),
                         'primitive_direction_envelope': envelope,
                         'strict_witnesses_checked_by_original_scheduler': len(witnesses) if full_ground else 0,
                         'large_box_axis_counts': axes, 'witnesses': witnesses})
    return fixtures


def generic_sector_controls():
    fixtures = []
    for n, A, B in ((4, 1, 1), (6, 2, 3), (9, 3, 3), (12, 2, 5)):
        ground = tuple((a, b, c) for a in range(n) for b in range(n) for c in range(n))
        for weight in ((Fraction(1), Fraction(A)+Fraction(17, 101), Fraction(B)+Fraction(29, 103)),
                       (Fraction(1), Fraction(A), Fraction(B))):
            prefix = frozenset(sorted(ground, key=lambda u: (dot(u, weight), u))[:n])
            assert independent([a+A*b+B*c for a, b, c in prefix], n)
            assert schedule_greedy(ground, weight, (1, A, B), n) == prefix
            fixtures.append({'length': n, 'valuations': [1, A, B],
                             'weight': list(map(str, weight)), 'standard_set': sorted(prefix)})
    return fixtures


def block_cut_controls():
    """Enumerate full small half-fans, then replay every active facet exchange."""
    fixtures = []
    for n, M in ((4, 2), (6, 2), (8, 3), (9, 3)):
        ground = chains(n, M, M)
        data = model(n, (1, M, M))
        groups = (n-1)//M+1
        witnesses, cells, facet_count, candidates = [], 0, 0, 0
        for lengths in compositions(n, len(ground)):
            if not feasible(lengths, ground):
                continue
            candidates += 1
            pairs = list(transfers(lengths, ground))
            directions = {tuple(x//gcd(*raw) for x in raw) for _, _, raw in pairs}
            polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
            for h in sorted(directions):
                polygon = clip(polygon, (h[0], h[1], h[1]+h[2]), 1)
                if not polygon:
                    break
            if not polygon or area2(polygon) == 0:
                continue
            cells += 1
            state = standard_set(lengths, ground)
            assert all(tuple(x-int(i == j) for i, x in enumerate(u)) in state
                       for u in state for j in range(3) if u[j])
            # Independent exact coefficient reconstruction of the same whole cone.
            _, coefficient_normals = border_basis(state, data['residue'])
            coefficient_polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
            for h in sorted(coefficient_normals):
                coefficient_polygon = clip(coefficient_polygon, (h[0], h[1], h[1]+h[2]), 1)
            assert set(coefficient_polygon) == set(polygon)
            cuts = [0] + [k for k in range(1, groups)
                          if sum(amount for amount, (q, _, _) in zip(lengths, ground)
                                 if q < M*k) == M*k] + [groups]
            for left, right in zip(polygon, polygon[1:]+polygon[:1]):
                transformed = tuple((a+b)/2 for a, b in zip(left, right))
                if min(transformed) <= 0:
                    continue
                facet_count += 1
                weight = (transformed[0], transformed[1]+transformed[2], transformed[2])
                active = [(i, j, raw) for i, j, raw in pairs if dot(weight, raw) == 0]
                assert active
                for i, j, raw in active:
                    removed = (lengths[i]-1, ground[i][1], ground[i][2])
                    added = (lengths[j], ground[j][1], ground[j][2])
                    donor_degree, recipient_degree = sum(removed[1:]), sum(added[1:])
                    donor_block = next(k for k in range(len(cuts)-1)
                                       if cuts[k] <= donor_degree < cuts[k+1])
                    recipient_block = next(k for k in range(len(cuts)-1)
                                           if cuts[k] <= recipient_degree < cuts[k+1])
                    lo = cuts[min(donor_block, recipient_block)]
                    hi = cuts[max(donor_block, recipient_block)+1]-1
                    selected = {(a, b, b+c-lo) for a, b, c in state if lo <= b+c <= hi}
                    size = len(selected)
                    assert size == (n if hi == groups-1 else M*(hi+1))-M*lo
                    level = dot(weight, removed)-transformed[2]*lo
                    assert level > 0 and dot(weight, removed) == dot(weight, added)
                    closed = {(a, b, e) for a in range(int(level/transformed[0])+1)
                              for e in range(hi-lo+1) for b in range(lo+e+1)
                              if dot(transformed, (a, b, e)) <= level}
                    strict = {u for u in closed if dot(transformed, u) < level}
                    assert selected <= closed and strict <= selected
                    assert len(strict) <= size-1 and len(closed) >= size+1
                    half_images = {}
                    for u in closed:
                        image = tuple(x//2 for x in u)
                        assert image in strict
                        half_images[image] = half_images.get(image, 0)+1
                    assert max(half_images.values()) <= 8
                    assert len(closed) <= 8*(size-1)
                    box = tuple(1+max(u[k] for u in closed) for k in range(3))
                    A, B, D = box
                    inner = {(a, b, e) for a in range((A-1)//4+1)
                             for b in range((B-1)//4+1)
                             for e in range((D-1)//4, (D-1)//2+1)}
                    assert inner <= closed and 64*len(inner) >= prod(box)
                    direction = (raw[0], raw[1], raw[1]+raw[2])
                    assert prod(abs(x)+1 for x in direction) <= prod(box)
                    assert prod(box) <= 64*len(closed) <= 512*(size-1) <= 512*n
                    witnesses.append({'allocation': lengths, 'transformed_weight': list(map(str, transformed)),
                                      'weight': list(map(str, weight)), 'transfer_indices': [i, j],
                                      'removed': removed, 'added': added, 'degree_interval': [lo, hi],
                                      'block_rank': size, 'translated_cutoff': str(level),
                                      'strict_cut_size': len(strict), 'closed_cut_size': len(closed),
                                      'closed_cut_box': box, 'direction': direction,
                                      'direction_product': prod(abs(x)+1 for x in direction)})
        fixtures.append({'length': n, 'M': M, 'feasible_allocations_examined': candidates,
                         'complete_half_fan_cells': cells, 'positive_facets': facet_count,
                         'independent_rational_border_cones': cells,
                         'active_facet_transfers_checked': len(witnesses), 'witnesses': witnesses})
    return fixtures


def main():
    exhaustive = exhaustive_controls()
    fans = existing_fan_controls()
    lines = line_controls()
    asymmetric = asymmetric_line_controls()
    tight = tight_stratum_controls()
    generic = generic_sector_controls()
    blocks = block_cut_controls()
    paths = ['scripts/run_nested_chain_checks.py', 'docs/nested_chain_allocation.md',
             'scripts/run_nested_jet_checks.py', 'results/research/nested_jet_checks.json',
             'scripts/run_structural_complexity_checks.py', 'scripts/run_extremal_initial_ideals_checks.py',
             'scripts/run_initial_cone_height_checks.py', 'scripts/run_interleaving_research.py']
    summary = {'exhaustive_allocation_matrix_checks': sum(f['all_compositions_matrix_checked'] for f in exhaustive),
               'feasible_allocations': sum(f['feasible_allocations'] for f in exhaustive),
               'transfer_feasibility_checks': sum(f['transfer_feasibility_checks'] for f in exhaustive),
               'local_optimality_checks': sum(f['local_optimality_checks'] for f in exhaustive),
               'previous_complete_cones_reconstructed': sum(f['complete_initial_cones_reconstructed'] for f in fans),
               'complete_weight_lines': len(lines),
               'complete_line_comparison_intervals': sum(f['full_comparison_intervals'] for f in lines),
               'strict_line_states': sum(f['distinct_strict_states'] for f in lines),
               'independent_line_matrix_checks': sum(f['independent_matrix_checks'] for f in lines),
               'asymmetric_complete_weight_lines': len(asymmetric),
               'asymmetric_line_comparison_intervals': sum(f['full_comparison_intervals'] for f in asymmetric),
               'asymmetric_strict_line_states': sum(f['distinct_strict_states'] for f in asymmetric),
               'tight_stratum_fixtures': len(tight),
               'strict_tight_stratum_witnesses': sum(f['distinct_joint_states_wy_above_wz'] for f in tight),
               'tight_stratum_scheduler_controls': sum(f['strict_witnesses_checked_by_original_scheduler'] for f in tight),
               'generic_sector_fixtures': len(generic),
               'block_cut_complete_half_fans': len(blocks),
               'block_cut_half_fan_cells': sum(f['complete_half_fan_cells'] for f in blocks),
               'block_cut_independent_border_cones': sum(f['independent_rational_border_cones'] for f in blocks),
               'block_cut_positive_facets': sum(f['positive_facets'] for f in blocks),
               'active_facet_cut_certificates': sum(f['active_facet_transfers_checked'] for f in blocks)}
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_nested_chain_checks.py',
              'scope': 'Exact chain reduction, local exchange cones, monotone cumulative constraints, all-tight stratum count and the consecutive-block cut certificate for the uniform equal-valuation family bound. Classical allocation, arrangement and planar counting inputs are not claimed as new. Unequal valuations, the full nested family and the unrestricted spatial endpoint remain open.',
              'summary': summary, 'exhaustive_controls': exhaustive, 'previous_fan_controls': fans,
              'monotone_line_controls': lines, 'tight_stratum_controls': tight,
              'asymmetric_monotone_line_controls': asymmetric,
              'generic_sector_controls': generic,
              'block_cut_controls': blocks,
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
              'equal_valuation_uniform_bound': 'O(n^2 log(n+1)^4), independent of M',
              'equal_valuation_power_separation_excluded': True,
              'unrestricted_endpoint_achieved': False, 'novelty_established': False, 'goal_complete': False}
    (ROOT/'results/research/nested_chain_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', **summary}, indent=2))


if __name__ == '__main__':
    main()
