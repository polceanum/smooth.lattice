#!/usr/bin/env python3
"""Exact absorbing-axis matroids and global nodal coupling controls."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import gcd
from pathlib import Path
from random import Random
import hashlib
import json
import sys

from run_forced_socle_fan_checks import model as forced_model
from run_projection_border_checks import operators, quotient
from run_projection_tag_checks import fan, polygon_key, TRIANGLE
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import area2, clip, determinant, dot, hyperbolic_cross
from run_univariate_branch_gluing_checks import (
    Frame, apply, clean, dense_kernel, dense_rank, pascal_ray_model,
)

ROOT = Path(__file__).resolve().parents[1]


def compact(value, depth=0):
    if isinstance(value, dict):
        if 'case_name' in value:
            return json.dumps(value, separators=(',', ':'))
        pad = '  '*(depth+1)
        return '{\n'+',\n'.join(pad+json.dumps(k)+': '+compact(v, depth+1)
                               for k, v in value.items())+'\n'+'  '*depth+'}'
    if isinstance(value, (tuple, list)) and any(isinstance(v, (dict, list, tuple)) for v in value):
        pad = '  '*(depth+1)
        return '[\n'+',\n'.join(pad+compact(v, depth+1) for v in value)+'\n'+'  '*depth+']'
    return json.dumps(value)


def absorbing_sweep(labels, rank, independent, H, name, counts):
    """All crossings, adjacent-swap tokens, and independent all-basis minima."""
    assert all(z == 0 or isinstance(x, int) and x >= 0 and z > 0 and
               (x+1)*(z+1) <= H for x, z in labels)
    m, L = len(labels), (H-1).bit_length()
    bases = [frozenset(b) for b in combinations(range(m), rank) if independent(frozenset(b))]
    assert bases

    def greedy(t):
        B = set()
        for j in sorted(range(m), key=lambda j: (labels[j][0]+t*labels[j][1], j)):
            if independent(frozenset(B | {j})):
                B.add(j)
        assert len(B) == rank
        return B

    events = defaultdict(set)
    for i, j in combinations(range(m), 2):
        dx, dz = (a-b for a, b in zip(labels[i], labels[j]))
        if dz and (t := Q(-dx, dz)) > 0:
            events[t].update((i, j))
    events = sorted(events.items())
    first = events[0][0]/2 if events else Q(1)
    order = sorted(range(m), key=lambda j: (labels[j][0], labels[j][1], j))
    positions = {j: i for i, j in enumerate(order)}
    B = greedy(first)
    token = {j: i for i, j in enumerate(sorted(B))}
    moves, local = defaultdict(lambda: [0, 0]), defaultdict(lambda: defaultdict(lambda: [0, 0]))
    paths = [[j] for j in sorted(B)]
    block = lambda j: labels[j][1].bit_length()
    initial = sum(block(j) for j in B)
    cross, drop, states, walls = 0, 0, [], []

    def record(t):
        assert B == greedy(t)
        cost = lambda C: sum(labels[j][0]+t*labels[j][1] for j in C)
        assert cost(B) == min(map(cost, bases)), 'independent exhaustive basis-cost check'
        counts['exhaustive_abstract_basis_minimum_intervals'] += 1
        if not states or states[-1]['basis'] != sorted(B):
            axes = {j for j in B if not labels[j][1]}
            if states:
                assert set(states[-1]['selected_axis_elements']) <= axes
            states.append({'basis': sorted(B), 'selected_axis_elements': sorted(axes),
                           'strict_parameter': str(t)})

    record(first)
    for k, (t, involved) in enumerate(events):
        groups = defaultdict(list)
        for j in involved:
            groups[labels[j][0]+t*labels[j][1]].append(j)
        pivots = []
        for group in groups.values():
            locations = sorted(positions[j] for j in group)
            assert locations == list(range(locations[0], locations[-1]+1))
            for loc, j in enumerate(sorted(group, key=lambda j: (labels[j][1], j)), locations[0]):
                while positions[j] > loc:
                    pos = positions[j]
                    e, f = order[pos-1], j
                    if e in B and f not in B and independent(frozenset(B-{e} | {f})):
                        dz, dx = labels[e][1]-labels[f][1], labels[f][0]-labels[e][0]
                        assert dz > 0 and dx > 0 and Q(dx, dz) == t
                        old, new = block(e), block(f)
                        assert old >= new
                        i = token.pop(e)
                        token[f] = i
                        paths[i].append(f)
                        if old != new:
                            cross += 1
                            drop += old-new
                        else:
                            assert old
                            moves[i, old][0] += dz
                            moves[i, old][1] += dx
                            local[old][t][0] += dz
                            local[old][t][1] += dx
                        B.remove(e)
                        B.add(f)
                        pivots.append([e, f, i])
                    order[pos-1], order[pos] = f, e
                    positions[f], positions[e] = pos-1, pos
        if pivots:
            walls.append({'parameter': str(t), 'pivots': pivots})
        nxt = events[k+1][0] if k+1 < len(events) else t+2
        record((t+nxt)/2)
    assert len(states) == len(walls)+1
    assert cross <= drop == initial-sum(block(j) for j in B) <= rank*L
    controls = []
    for label, rows in sorted(local.items()):
        j = label-1
        for i in range(rank):
            dz, dx = moves[i, label]
            assert dz <= 2**j-1 and dx <= H//2**j
        A, C = (sum(row[a] for row in rows.values()) for a in range(2))
        K = len(rows)
        assert K**3 <= 27*A*C <= 27*rank**2*H
        controls.append({'block': label, 'distinct_parameters': K, 'z_budget_used': A,
                         'x_budget_used': C})
    assert len(walls) <= cross+sum(c['distinct_parameters'] for c in controls)
    assert max(0, len(states)-1-rank*L)**3 <= 27*L**3*rank**2*H
    counts['abstract_adjacent_swaps_with_actual_basis_change'] += sum(len(w['pivots']) for w in walls)
    return {'case_name': name, 'rank': rank, 'H': H, 'labels': labels,
            'states': states, 'walls': walls, 'token_paths': paths,
            'cross_block_moves': cross, 'dyadic_controls': controls,
            'exhaustive_represented_or_abstract_bases': len(bases)}


def abstract_controls(counts):
    labels = ((-10**6, 0), (7, 0), (10**8, 0), (7, 0),
              (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (0, 3), (1, 3))
    rows = []
    for r in (1, 3, 4):
        rows.append(absorbing_sweep(labels, r, lambda s: len(s) <= r, 12,
                                    'uniform_rank_'+str(r), counts))
        for p in (0, 2, 5):
            random = Random(271828+100*r+p)
            columns = [clean({i: int(i == j) for i in range(r)}, p) if j < r else
                       clean({i: random.randrange(-3, 4) for i in range(r)}, p)
                       for j in range(len(labels))]
            columns[-1] = {}  # a genuine loop, including at cheap parameter values
            @lru_cache(None)
            def independent(s):
                return len(s) <= r and dense_rank([columns[j] for j in sorted(s)], r, p) == len(s)
            rows.append(absorbing_sweep(labels, r, independent, 12,
                        'linear_rank_'+str(r)+'_field_'+str(p), counts))
    rows.append(absorbing_sweep(((-10**9, 0), (Q(1, 7), 0), (10**10, 0)), 1,
                                lambda s: len(s) <= 1, 1, 'axis_only_real_labels', counts))
    # The rank-four Vamos sparse-paving matroid: no field representation is used.
    pairs = [frozenset((2*i, 2*i+1)) for i in range(4)]
    forbidden = {a | b for a, b in combinations(pairs, 2) if a | b != pairs[2] | pairs[3]}
    bases = {frozenset(b) for b in combinations(range(8), 4)}-forbidden
    assert all(any(B-{e} | {f} in bases for f in C-B)
               for B in bases for C in bases for e in B-C)
    independent = lambda s: len(s) < 4 or len(s) == 4 and s not in forbidden
    labels_v = ((10**6, 0), (-10**6, 0), (0, 1), (1, 1), (0, 2), (1, 2), (0, 3), (1, 3))
    rows.append(absorbing_sweep(labels_v, 4, independent, 12, 'vamos_basis_exchange_verified', counts))
    return rows


def counter_model(N, p):
    reference = tuple([(i, 0, 0) for i in range(N)]+[(0, 1, 0)])
    @lru_cache(None)
    def ev(a):
        x, y, z = a
        if y:
            return {N: 1} if y == 1 and not x and not z else {}
        if not x and z == 1:
            return clean({1: 1, N: 1}, p)
        return {x+z: 1} if x+z < N else {}
    return reference, ev


def moment(B):
    return tuple(sum(a[i] for a in B) for i in range(3))


def exchange_intervals(labels, columns, p, counts):
    """Two separable prices: independent prefix ranks and actual adjacent swaps."""
    rank = lambda S: len(Frame([columns[j] for j in sorted(S)], p).basis)
    intervals = []
    for axis in (1, 2):
        own = [j for j, a in enumerate(labels) if a[axis]]
        other = [j for j, a in enumerate(labels) if a[3-axis]]
        for e, f in combinations(own, 2):
            a, b = labels[e], labels[f]
            if a[axis] == b[axis]:
                continue
            u = Q(b[0]-a[0], a[axis]-b[axis])
            if u <= 0:
                continue
            c = a[0]+u*a[axis]
            C = {j for j in own if labels[j][0]+u*labels[j][axis] < c}
            events = sorted({Q(c-labels[j][0], labels[j][3-axis]) for j in other
                             if c > labels[j][0]})
            bounds = [Q(0)]+events
            samples = [(a+b)/2 for a, b in zip(bounds, bounds[1:])]+[bounds[-1]+1]
            active = []
            for k, t in enumerate(samples):
                D = {j for j in other if labels[j][0]+t*labels[j][3-axis] < c}
                S = C | D
                r = rank(S)
                yes = rank(S | {e}) == rank(S | {f}) == rank(S | {e, f}) == r+1
                if yes:
                    active.append(k)
                # With no third element at this tie, replay the two greedy
                # orders directly rather than inferring selection from ranks.
                if all(labels[j][0]+u*labels[j][axis] != c for j in set(own)-{e, f}):
                    w = [Q(1), Q(0), Q(0)]
                    w[axis], w[3-axis] = u, t
                    order = sorted(range(len(labels)), key=lambda j: (dot(labels[j], w), j))
                    i, j = order.index(e), order.index(f)
                    assert abs(i-j) == 1
                    before = set(field_basis(order, lambda j: columns[j], p)[0])
                    order[i], order[j] = order[j], order[i]
                    after = set(field_basis(order, lambda j: columns[j], p)[0])
                    assert (before != after) == yes
                    if yes:
                        assert before ^ after == {e, f}
                    counts['independent_separable_exchange_swap_replays'] += 1
                counts['independent_separable_contracted_pair_prefixes'] += 1
            assert not active or active == list(range(active[0], active[-1]+1))
            intervals.append({'axis': axis, 'pair': [e, f], 'fixed_price': str(u),
                              'other_price_breakpoints': list(map(str, events)),
                              'active_open_intervals': active})
            counts['single_interval_separable_exchange_carriers'] += 1
    return intervals


def socle_model(labels, columns, p):
    """An actual monomial downset quotient by an annihilated tip-space kernel."""
    tips = set(labels)
    assert len(tips) == len(labels)
    assert all(not all(a[i] <= b[i] for i in range(3)) for a, b in permutations_pair(tips))
    D = {b for a in tips for b in product(*(range(t+1) for t in a))}
    core = tuple(sorted(D-tips))
    index = {a: i for i, a in enumerate(core)}
    tip_index = {a: j for j, a in enumerate(labels)}
    r = len(columns)//2
    assert columns[:r] == [{i: 1} for i in range(r)]
    @lru_cache(None)
    def ev(a):
        if a in index:
            return {index[a]: 1}
        if a in tip_index:
            return clean({len(core)+i: c for i, c in columns[tip_index[a]].items()}, p)
        return {}
    return core+tuple(labels[:r]), ev, core


def permutations_pair(items):
    return ((a, b) for a in items for b in items if a != b)


def internal_axis_control(counts):
    labels = ((4, 3, 0), (3, 4, 0), (0, 6, 0), (9, 1, 0),
              (0, 0, 6), (9, 0, 1), (3, 0, 5), (4, 0, 4))
    columns = [{i: 1} for i in range(4)]+[{0: 1, 1: 1}, {0: 1}, {3: 1}, {2: 1, 3: 1}]
    minors = {b: determinant([[columns[j].get(i, 0) for j in b] for i in range(4)])
              for b in combinations(range(8), 4)}
    assert all(d in (-1, 0, 1) for d in minors.values())
    integral_bases = {b for b, d in minors.items() if d}
    assert len(integral_bases) == 25
    counts['independent_integral_internal_tip_maximal_minors'] += len(minors)
    cases = []
    for p in (0, 2, 5):
        reference, ev, core = socle_model(labels, columns, p)
        n = len(reference)
        assert n == 62
        M = operators(reference, ev)
        for U, V in combinations(M, 2):
            for j in range(n):
                assert apply(U, V[j], p) == apply(V, U[j], p)
                counts['independent_internal_obstruction_commuting_columns'] += 1
        assert all(not apply(M[1], v, p) for v in M[2])
        ny, ey, _ = quotient(M, ev((0, 0, 0)), 2, p, counts)
        nz, ez, _ = quotient(M, ev((0, 0, 0)), 1, p, counts)
        Y, Z = fan(ny, ey, p, counts), fan(nz, ez, p, counts)
        assert len(Y) == len(Z) == 1 and (ny, nz) == (32, 36)
        assert next(iter(Y)) == frozenset(a for a in core if not a[2])
        assert next(iter(Z)) == frozenset(a for a in core if not a[1])
        parent = fan(n, ev, p, counts)
        assert len(parent) == 15
        tip_cells = {}
        bases = [b for b in combinations(range(8), 4) if len(Frame([columns[j] for j in b], p).basis) == 4]
        assert set(bases) == integral_bases
        for b in bases:
            mu = moment([labels[j] for j in b])
            poly = list(TRIANGLE)
            for other in bases:
                delta = tuple(a-c for a, c in zip(moment([labels[j] for j in other]), mu))
                if any(delta):
                    poly = clip(poly, delta, 1)
                if not poly:
                    break
            if poly and area2(poly) > 0:
                B = frozenset(core) | {labels[j] for j in b}
                assert B in parent and polygon_key(poly) == polygon_key(parent[B])
                tip_cells[b] = poly
        assert len(tip_cells) == len(parent) and sum(area2(v) for v in tip_cells.values()) == 1
        vertex = (Q(1, 3),)*3
        incident = [e for e in edges(parent) if vertex in e[:2]]
        assert len(incident) == 4
        assert all(not e[2][1] or not e[2][2] for e in incident)
        local = []
        for u, v in product((Q(99, 100), Q(101, 100)), repeat=2):
            w = (1, u, v)
            chosen = min(tip_cells, key=lambda b: dot(moment([labels[j] for j in b]), w))
            assert all(dot(moment([labels[j] for j in b]), w) >
                       dot(moment([labels[j] for j in chosen]), w) for b in bases if b != chosen)
            local.append({'weight': list(map(str, w)), 'tip_basis': list(chosen)})
        assert len({tuple(row['tip_basis']) for row in local}) == 4
        counts['actual_internal_axis_junction_obstructions'] += 1
        cases.append({'case_name': 'length62_pure_axis_internal_field_'+str(p),
                      'length': n, 'characteristic': p, 'compulsory_core': core,
                      'tip_labels': labels, 'tip_columns': columns, 'parent_initials': len(parent),
                      'coordinate_quotient_lengths': [ny, nz], 'coordinate_quotient_initials': [1, 1],
                      'complete_tip_basis_enumeration': len(bases), 'local_four_cone_witnesses': local,
                      'integral_tip_maximal_minors': {str(v): sum(d == v for d in minors.values())
                                                     for v in (-1, 0, 1)},
                      'incident_wall_jumps': [e[2] for e in incident],
                      'mixed_quotient_grid_incidence': 0,
                      'separable_exchange_carriers': exchange_intervals(labels, columns, p, counts)})
    return cases


def growing_axis_controls(counts):
    """Compressed exact certificates; enormous cores are specified, not enumerated."""
    rows = []
    for s in (1, 2, 4, 8):
        m = 2*s
        Q0, G = 1000*(m+1)**6, 100*(m+1)**4
        K = 10*(Q0+G+m)
        D = K*(m+2)
        assert Q0 > 2*G and 2*D*(m+1) < G*Q0
        Y, Z, extra = [], [], []
        for i in range(1, m+1):
            x, b = i*K, D-i*K
            tie = [(x, b), (x+Q0+i, b-Q0)]
            bracket = [(x, b-G), (x+Q0, b-Q0+G)]
            yp, zp = (tie, bracket) if i <= s else (bracket, tie)
            Y.extend((a, b, 0) for a, b in yp)
            Z.extend((a, 0, c) for a, c in zp)
            j = 2*(i-1)
            extra.extend(([{j: 1, j+1: 1}, {j: 1}] if i <= s else
                          [{j+1: 1}, {j: 1, j+1: 1}]))
        labels = tuple(Y+Z)
        assert all(not all(a[i] <= b[i] for i in range(3)) for a, b in permutations_pair(labels))
        def plane_area(tips, axis):
            previous, area = -1, 0
            for a in sorted(tips):
                area += (a[0]-previous)*(a[axis]+1)
                previous = a[0]
            return area
        length = plane_area(Y, 1)+plane_area(Z, 2)-min(max(a[0] for a in Y), max(a[0] for a in Z))-1-2*m
        columns = [{i: 1} for i in range(2*m)]+extra
        cut = lambda indices: [Q(1)]+[1+Q(i, Q0) for i in indices]+[1+Q(m+1, Q0)]
        uc, vc = cut(range(1, s+1)), cut(range(s+1, m+1))
        us, vs = [(a+b)/2 for a, b in zip(uc, uc[1:])], [(a+b)/2 for a, b in zip(vc, vc[1:])]
        certificates = []
        for p in (0, 2, 5):
            states = set()
            for u, v in product(us, vs):
                w, selected = (1, u, v), []
                for i in range(m):
                    ids = (2*i, 2*i+1, 2*m+2*i, 2*m+2*i+1)
                    bases = [b for b in combinations(ids, 2) if len(Frame([columns[j] for j in b], p).basis) == 2]
                    ordered = sorted(bases, key=lambda b: dot(moment([labels[j] for j in b]), w))
                    assert dot(moment([labels[j] for j in ordered[0]]), w) < dot(moment([labels[j] for j in ordered[1]]), w)
                    chosen = set(ordered[0])
                    if i < s:
                        assert 2*m+2*i in chosen and len(chosen & {2*i, 2*i+1}) == 1
                    else:
                        assert 2*i in chosen and len(chosen & {2*m+2*i, 2*m+2*i+1}) == 1
                    selected.extend(ordered[0])
                assert len(Frame([columns[j] for j in selected], p).basis) == 2*m
                states.add(tuple(sorted(selected)))
                counts['compressed_growing_axis_grid_cells'] += 1
            assert len(states) == (s+1)**2
            certificates.append({'characteristic': p, 'different_local_bases': len(states),
                                 'pure_axis_internal_junctions': s*s})
        rows.append({'case_name': 'compressed_axis_obstruction_s'+str(s), 'parameter_s': s,
                     'coordinate_quotient_initials': [1, 1], 'mixed_grid_incidence': 0,
                     'length_from_downset_telescope': length, 'Q': Q0, 'G': G, 'K': K, 'D': D,
                     'tip_labels': labels, 'tip_columns': columns, 'field_certificates': certificates,
                     'scope': 'Compressed socle-matroid and antichain certificate for the proved algebra family; '
                              'the huge compulsory monomial core is not enumerated and no full large algebra fan is claimed.'})
    return rows


def energy(cells, w):
    return min(dot(moment(B), w) for B in cells)


def edges(cells):
    result = []
    for (B, P), (C, R) in combinations(cells.items(), 2):
        common = sorted(set(P) & set(R))
        if len(common) < 2:
            continue
        a, b = common[0], common[-1]
        mid = tuple((u+v)/2 for u, v in zip(a, b))
        if not all(t > 0 for t in mid):
            continue
        jump = tuple(u-v for u, v in zip(moment(B), moment(C)))
        assert any(jump) and dot(jump, a) == dot(jump, b) == 0
        scale = gcd(*jump)
        primitive = tuple(t//scale for t in jump)
        result.append((a, b, jump, primitive, scale))
    return result


def intersection(a, b, axis, value):
    da, db = a[axis]-value*a[0], b[axis]-value*b[0]
    if da == db == 0:
        return None  # coincident walls are excluded for generic sweeps
    if da*db > 0:
        return None
    t = da/(da-db)
    point = tuple(x+t*(y-x) for x, y in zip(a, b))
    return point if all(c > 0 for c in point) else None


def filtered_integral(n, ev, p, M, w, expected, counts):
    ground = [a for a in hyperbolic_cross(n, 3) if not a[1]*a[2]]
    Y, Z = (Frame(op, p).basis for op in M[1:])
    joined = Y+Z
    relation, joint_rank = dense_kernel(joined, n, p)
    h = len(joined)-joint_rank
    W = tuple(apply(Y, {j: c for j, c in v.items() if j < len(Y)}, p) for v in relation)
    assert dense_rank(W, n, p) == h
    thresholds = sorted({Q(0)} | {dot(a, w) for a in ground})
    integral, max_defect = Q(0), 0
    certificates = []
    for t, nxt in zip(thresholds, thresholds[1:]):
        V = Frame([ev(a) for a in ground if dot(a, w) <= t], p).basis
        r = len(V)
        ry = dense_rank(V+Z, n, p)-len(Z)
        rz = dense_rank(V+Y, n, p)-len(Y)
        rb = dense_rank(V+joined, n, p)-joint_rank
        k = r+h-dense_rank(V+W, n, p)
        defect = ry+rz-rb-r+k
        assert 0 <= k <= h and defect >= 0
        integrand = h-k+defect
        integral += (nxt-t)*integrand
        max_defect = max(max_defect, defect)
        certificates.append([str(t), str(nxt), r, ry, rz, rb, k, defect])
        counts['independent_filtered_rank_defect_bands'] += 1
    assert integral == expected >= 0
    return {'weight': [str(t) for t in w], 'integral': str(integral), 'overlap_dimension': h,
            'maximum_compatibility_defect': max_defect, 'bands_t_next_r_ry_rz_rb_k_defect': certificates}


def inspect(name, reference, ev, p, counts):
    n = len(reference)
    M = operators(reference, ev)
    assert all(ev(a) == {i: 1} for i, a in enumerate(reference))
    for U, V in combinations(M, 2):
        for j in range(n):
            assert apply(U, V[j], p) == apply(V, U[j], p)
            counts['independent_commuting_parent_columns'] += 1
    assert all(not apply(M[1], v, p) for v in M[2])
    ny, ey, _ = quotient(M, ev((0, 0, 0)), 2, p, counts)
    nz, ez, _ = quotient(M, ev((0, 0, 0)), 1, p, counts)
    parent = fan(n, ev, p, counts)
    Y = fan(ny, ey, p, counts) if ny else {frozenset(): TRIANGLE}
    Z = fan(nz, ez, p, counts) if nz else {frozenset(): TRIANGLE}
    joint_rank = len(Frame(M[1]+M[2], p).basis)
    d = n-joint_rank
    axis_cost = d*(d-1)//2
    potential = lambda w: energy(parent, w)-energy(Y, w)-energy(Z, w)+axis_cost*w[0]
    # Every full-dimensional overlay patch is checked, not just sample gradients.
    from run_signed_face_checks import normals
    from run_structural_complexity_checks import clip, area2
    patches, overlay_area = 0, Q(0)
    for BY, BZ in product(Y, Z):
        constraints = (() if not ny else normals(tuple(sorted(BY)), ey, p)[0]) + \
                      (() if not nz else normals(tuple(sorted(BZ)), ez, p)[0])
        for B, poly in parent.items():
            patch = list(poly)
            for q in constraints:
                patch = clip(patch, q, 1)
                if not patch:
                    break
            if not patch or area2(patch) <= 0:
                continue
            gradient = tuple(u-v-t for u, v, t in zip(moment(B), moment(BY), moment(BZ)))
            assert 0 <= gradient[1] <= ny and 0 <= gradient[2] <= nz
            patches += 1
            overlay_area += area2(patch)
            counts['complete_overlay_residual_gradient_bounds'] += 1
    assert overlay_area == 1, 'complete parent/quotient overlay coverage'
    parent_edges = edges(parent)
    mixed = [e for e in parent_edges if e[2][1]*e[2][2] < 0]
    sweeps = []
    ground = [a for a in hyperbolic_cross(n, 3) if not a[1]*a[2]]
    denominator = 2*n*n+1
    for axis, cap in ((1, ny), (2, nz)):
        moving = 3-axis
        for fixed in (Q(1, denominator), Q(n+1, denominator), 1+Q(1, denominator), 5+Q(1, denominator)):
            events = sorted({pt[moving]/pt[0] for a, b, *_ in parent_edges
                             if (pt := intersection(a, b, axis, fixed)) is not None})
            bounds = [Q(0)]+events
            samples = [(a+b)/2 for a, b in zip(bounds, bounds[1:])]
            samples.append(bounds[-1]+1)
            previous, states = set(), []
            for t in samples:
                w = [Q(1), Q(0), Q(0)]
                w[axis], w[moving] = fixed, t
                B, _ = field_basis(sorted(ground, key=lambda a: (dot(a, w), a)), ev, p)
                assert len(B) == n and frozenset(B) in parent
                unchanged = {a for a in B if not a[moving]}
                assert previous <= unchanged, 'gross substitutes for the unchanged-price elements'
                previous = unchanged
                if not states or states[-1]['basis'] != sorted(B):
                    states.append({'basis': sorted(B), 'parameter': str(t)})
                counts['independent_actual_shadow_greedy_replays'] += 1
            assert len(states) <= 1+4*n*(n-1).bit_length()
            jump_budget = sum(abs(e[2][axis]) for e in mixed
                              if intersection(e[0], e[1], axis, fixed) is not None)
            assert jump_budget <= cap
            sweeps.append({'fixed_axis': axis, 'fixed_price': str(fixed), 'states': states,
                           'mixed_lattice_jump_budget_used': jump_budget, 'budget': cap})
    grids, total = [], 0
    for axis, cells, cap in ((1, Y, ny), (2, Z, nz)):
        ratios = sorted({Q(-e[2][0], e[2][axis]) for e in edges(cells) if e[2][axis]})
        assert len(ratios) == len(cells)-1
        for fixed in ratios:
            incidence = [i for i, e in enumerate(mixed)
                         if intersection(e[0], e[1], axis, fixed) is not None]
            budget = sum(abs(mixed[i][2][axis]) for i in incidence)
            assert budget <= 2*cap
            total += budget
            grids.append({'axis': axis, 'price': str(fixed), 'mixed_edges': incidence,
                          'weighted_incidence': budget, 'bound': 2*cap})
            counts['exact_moving_quotient_grid_line_budgets'] += 1
    assert total <= 2*ny*(len(Y)-1)+2*nz*(len(Z)-1)
    prices = (Q(1, 4), Q(3, 4), Q(5, 4), Q(2), Q(3))
    for u1, u2 in combinations(prices, 2):
        for v1, v2 in combinations(prices, 2):
            a, b = potential((1, u1, v1)), potential((1, u2, v2))
            c, e = potential((1, u1, v2)), potential((1, u2, v1))
            assert min(a, b, c, e) >= 0 and a+b >= c+e
            assert a <= c and a <= e and c <= b and e <= b
            counts['exact_supermodular_and_monotone_rectangles'] += 1
    w = (Q(1), Q(7, 5), Q(11, 6))
    integral = filtered_integral(n, ev, p, M, w, potential(w), counts)
    return {'case_name': name, 'characteristic': p, 'length': n, 'parent_initials': len(parent),
            'coordinate_quotient_lengths': [ny, nz], 'coordinate_quotient_initials': [len(Y), len(Z)],
            'complete_overlay_patches': patches, 'coordinate_sweeps': sweeps,
            'mixed_edges': [{'ends': [[str(t) for t in v] for v in e[:2]], 'jump': e[2],
                             'primitive': e[3], 'lattice_multiplicity': e[4]} for e in mixed],
            'moving_quotient_grid_lines': grids, 'total_weighted_grid_incidence': total,
            'filtered_rank_certificate': integral}, parent, potential


def main():
    counts, cases, counter = Counter(), [], []
    abstract = abstract_controls(counts)
    internal = internal_axis_control(counts)
    growing_axis = growing_axis_controls(counts)
    for p in (0, 2, 5):
        for N in (3, 5, 8):
            reference, ev = counter_model(N, p)
            row, cells, potential = inspect('nonconcave_zero_overlap_N'+str(N), reference, ev, p, counts)
            predicted = {
                frozenset([(i, 0, 0) for i in range(N)]+[(0, 1, 0)]),
                frozenset([(0, 0, i) for i in range(N)]+[(0, 1, 0)]),
                frozenset([(i, 0, 0) for i in range(N)]+[(0, 0, 1)]),
                frozenset([(0, 0, i) for i in range(N)]+[(1, 0, 0)])}
            assert set(cells) == predicted
            assert row['coordinate_quotient_lengths'] == [2, N]
            assert row['coordinate_quotient_initials'] == [2, 2]
            for w in product((Q(1, 4), Q(3, 4), Q(5, 4), Q(2), Q(3)), repeat=3):
                assert potential(w) == max(0, min(w[1:])-w[0])
                counts['independent_counterfamily_global_potential_values'] += 1
            w1, w2, mid = (1, Q(1, 4), 3), (1, Q(5, 4), 3), (1, Q(3, 4), 3)
            vals = [potential(w) for w in (w1, w2, mid)]
            assert vals == [0, Q(1, 4), 0] and vals[2] < (vals[0]+vals[1])/2
            counter.append({'case_name': 'N'+str(N)+'_field_'+str(p),
                            'nonconcavity_weights': [[str(t) for t in w] for w in (w1, w2, mid)],
                            'potential_values': list(map(str, vals))})
            cases.append(row)
        for parameter, kind in ((3, 'H'), (4, 'A'), (4, 'C')):
            reference, ev = forced_model(parameter, kind, p)
            cases.append(inspect('nonrow_forced_'+kind, reference, ev, p, counts)[0])
        for m in (2, 3, 4):
            n, ev = pascal_ray_model(m, p)
            cases.append(inspect('pascal_overlap_m'+str(m), tuple((i, 0, 0) for i in range(n)),
                                 ev, p, counts)[0])
    paths = {ROOT/'scripts/run_nodal_shadow_checks.py', ROOT/'docs/nodal_shadow_potential.md'}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.suffix == '.py' and path.parent == ROOT/'scripts':
                paths.add(path)
    out = {'command': 'python3 -B scripts/run_nodal_shadow_checks.py', 'all_checks_passed': True,
           'scope': 'Global coordinate shadows, bounded residual potential and mixed-grid incidence for yz=0. '
                    'No unrestricted count improvement or novelty claim. Classical greedy exchange is an input.',
           'abstract_matroids': abstract, 'actual_quotients': cases,
           'pure_axis_internal_junction_obstructions': internal,
           'compressed_growing_axis_obstructions': growing_axis,
           'uniform_nonconcavity_counterfamily': counter, 'counters': dict(counts),
           'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in sorted(paths)},
           'unrestricted_G3_improved': False, 'polynomial_lower_separation_proved': False,
           'novelty_established': False, 'goal_complete': False}
    # Fraction-valued arbitrary axis labels are serialized exactly.
    out = json.loads(json.dumps(out, default=str))
    serialized = compact(out)+'\n'
    assert json.loads(serialized) == out
    (ROOT/'results/research/nodal_shadow_checks.json').write_text(serialized)
    print(json.dumps(dict(counts), sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
