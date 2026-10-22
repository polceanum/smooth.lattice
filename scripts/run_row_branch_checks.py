#!/usr/bin/env python3
"""Arbitrary row-polynomial branches: exact coupling and profile budgets."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import sys

from run_minimal_branch_overlap_checks import mul, power, rem, partitions
from run_fitting_branch_core_checks import quotient_by_relations
from run_univariate_branch_gluing_checks import (
    Frame, apply, clean, cyclic_span, dense_kernel, dense_rank, linear,
    pascal_ray_model, pascal_ray_predicted, polynomial,
)
from run_projection_tag_checks import fan, compare_fans
from run_forced_socle_fan_checks import model as forced_model
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import area2, dot, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]


def envelope(drops, P, p):
    """No ordering or divisibility hypothesis on the row drops."""
    factors = tuple(power(P, a, p) for a in drops)
    gs = [None]*len(drops)+[(1,)]
    for j in range(len(drops)-1, -1, -1):
        gs[j] = mul(factors[j], gs[j+1], p)
    h = len(gs[0])-1
    largest = max(drops)
    D = largest*(len(P)-1)
    f = mul(gs[0], power(P, largest, p), p)
    labels = tuple((i, 0, 0) for i in range(len(f)-1))
    labels += tuple((i, j+1 if b == 1 else 0, j+1 if b == 2 else 0)
                    for b in (1, 2) for j, g in enumerate(gs[:-1]) for i in range(len(g)-1))
    index = {a: i for i, a in enumerate(labels)}
    def row(poly, j=0, b=0):
        if b and j == len(drops):
            return {}
        coeffs = rem(poly, gs[j] if b else f, p)
        return clean({index[i, j+1 if b == 1 else 0, j+1 if b == 2 else 0]: c
                      for i, c in enumerate(coeffs) if c}, p)
    X, Y, Z = [], [], []
    for i, a, b in labels:
        branch, j = (1, a-1) if a else (2, b-1) if b else (0, 0)
        X.append(row((0,)*(i+1)+(1,), j, branch))
        Y.append(row((0,)*i+(1,), j+1, 1) if branch == 1 else
                 row((0,)*i+(1,), 0, 1) if not branch else {})
        Z.append(row((0,)*i+(1,), j+1, 2) if branch == 2 else
                 row((0,)*i+(1,), 0, 2) if not branch else {})
    M = tuple(X), tuple(Y), tuple(Z)
    WY = tuple(row(gs[j+1], j, 1) for j in range(len(drops)))
    WZ = tuple(row(gs[j+1], j, 2) for j in range(len(drops)))
    E = frozenset((i, j+1 if b == 1 else 0, j+1 if b == 2 else 0)
                  for b in (1, 2) for j, g in enumerate(gs[:-1])
                  for i in range(len(gs[j+1])-1, len(g)-1))
    E |= frozenset((i, 0, 0) for i in range(h, h+D))
    S = tuple(row((0,)*i+gs[j+1], j, b)
              for b in (1, 2) for j, g in enumerate(gs[:-1])
              for i in range(len(g)-len(gs[j+1])))
    S += tuple(row((0,)*i+gs[0]) for i in range(D))
    return dict(drops=drops, P=P, gs=gs, h=h, D=D, f=f, labels=labels,
                M=M, WY=WY, WZ=WZ, E=E, C=frozenset(labels)-E, S=S,
                unit={index[0, 0, 0]: 1}, row=row)


def coupling(H, p, core, shear):
    X, _, _ = H['M']
    relations = []
    for j, (a, v) in enumerate(zip(H['drops'], H['WY'])):
        image = H['WZ'][j]
        if shear:
            image = linear([(1, image)]+[(i+j+1, polynomial(X, H['WZ'][i],
                        power(H['P'], max(0, H['drops'][i]-a), p), p))
                        for i in range(j)], p)
        relations.append(linear(((1, v), (-1, image)), p))
    if core:
        j = H['drops'].index(max(H['drops']))
        relations.append(linear(((1, H['row'](H['gs'][0])), (-1, H['WY'][j])), p))
    return tuple(relations)


def inspect(drops, P, p, core, shear, complete, counts):
    H = envelope(drops, P, p)
    N = len(H['labels'])
    X, Y, Z = H['M']
    assert cyclic_span(H['M'], H['unit'], p) == N
    direct, rank = dense_kernel(tuple({**Y[i], **{N+j: c for j, c in Z[i].items()}}
                                     for i in range(N)), 2*N, p)
    assert rank+len(H['S']) == N and len(direct) == len(H['S'])
    assert len(Frame(H['S'], p).basis) == len(H['S'])
    F = Frame(H['S'], p)
    for v in direct:
        F.coords(v)
    for v in H['S']:
        assert not apply(Y, v, p) and not apply(Z, v, p)
    # Independent leading terms, all at divisor tails within their row.
    heads = set()
    for v in H['S']:
        head = max((H['labels'][i] for i in v), key=lambda a: a[0])
        assert all(all(b <= a for b, a in zip(H['labels'][i], head)) for i in v)
        heads.add(head)
    assert heads == H['E']
    m = sum(len(g)-1 for g in H['gs'][:-1])
    assert sum(a[1] for a in H['E']) == m == sum(a[2] for a in H['E'])
    counts['independent_joint_annihilator_and_row_budget_checks'] += 1
    T = coupling(H, p, core, shear)
    n, M, unit = quotient_by_relations(H['M'], H['unit'], T, p)
    assert n == N-H['h']-(H['D'] if core else 0)
    assert N <= 2*n
    for A, B in combinations(M, 2):
        for j in range(n):
            assert apply(A, B[j], p) == apply(B, A[j], p)
            counts['independent_commuting_coordinate_columns'] += 1
    @lru_cache(None)
    def ev(a):
        v = unit
        for op, e in zip(M, a):
            for _ in range(e):
                v = apply(op, v, p)
        return v
    blocks = [[a for a in H['labels'] if not a[1] and not a[2]]
              if b == 0 else [a for a in H['labels'] if a[b]] for b in range(3)]
    for block in blocks:
        assert len(field_basis(block, ev, p)[0]) == len(block)
    # Remove the joint annihilator: direct polynomial division gives this
    # fixed quotient's coordinate columns independently of quotient matrices.
    C = tuple(sorted(H['C']))
    ci = {a: i for i, a in enumerate(C)}
    def fixed_ev(a):
        i, y, z = a
        if y and z:
            return {}
        j, branch = (y-1, 1) if y else (z-1, 2) if z else (0, 0)
        if branch and j >= len(drops)-1:
            return {}
        bound = H['gs'][j+1] if branch else H['gs'][0]
        return clean({ci[k, y, z]: c for k, c in enumerate(rem((0,)*i+(1,), bound, p)) if c}, p)
    selected, pivots = field_basis(C, ev, p)
    assert len(selected) == len(C)
    Fparent = Frame(tuple(ev(a) for a in H['labels']), p)
    selected_labels = [H['labels'][i] for i in Fparent.selected]
    for a in H['labels']:
        co = Fparent.coords(ev(a))
        assert linear(((c, fixed_ev(selected_labels[j])) for j, c in co.items()), p) == fixed_ev(a)
        counts['independent_actual_fixed_quotient_columns'] += 1
    weights = ((1, 2, 7), (9, 2, 1), (2, 11, 3), (5, 3, 17))
    states = fan(n, ev, p, counts) if complete else {}
    samples = list(weights)
    samples += [tuple(sum(v[i] for v in poly) for i in range(3)) for poly in states.values()]
    for w in samples:
        B, _ = field_basis(sorted(H['labels'], key=lambda a: (dot(a, w), a)), ev, p)
        B = frozenset(B)
        assert H['C'] <= B <= frozenset(H['labels'])
        variable = B-H['C']
        assert sum(a[1] for a in variable) <= m and sum(a[2] for a in variable) <= m
        assert sum(a[0] for a in variable) <= N*(N-1)//2
        counts['independent_compulsory_core_and_entire_profile_budgets'] += 1
    phi = (len(P)-1)*sum((j+1)*d for j, d in enumerate(sorted(drops, reverse=True)))
    return {'row_drop_exponents': drops, 'base_polynomial': P, 'characteristic': p,
            'core_coupling': core, 'off_diagonal_shear': shear,
            'envelope_length': N, 'quotient_length': n, 'branch_length': m,
            'full_kernel_branch_excess': 2*(m-phi), 'compulsory_labels': C,
            'residual_labels': sorted(H['E']), 'transverse_budget_each_branch': m,
            'complete_initial_ideals': len(states) if complete else None}


def moving_annihilator_controls(counts):
    rows = []
    for p in (0, 2):
        for parameter in (3, 4, 6):
            reference, ev = forced_model(parameter, 'A', p)
            n = len(reference)
            M = tuple(tuple(ev(tuple(a[j]+(j == i) for j in range(3)))
                            for a in reference) for i in range(3))
            L, rank = dense_kernel(tuple({**M[1][i], **{n+j: c for j, c in M[2][i].items()}}
                                        for i in range(n)), 2*n, p)
            child_n, child_M, child_unit = quotient_by_relations(M, ev((0, 0, 0)), L, p)
            assert child_n == rank
            ry, rz = (len(Frame(op, p).basis) for op in child_M[1:])
            assert len(Frame(child_M[1]+child_M[2], p).basis) == ry+rz
            counts['actual_annihilator_quotients_have_zero_branch_overlap'] += 1
            @lru_cache(None)
            def child_ev(a):
                v = child_unit
                for op, e in zip(child_M, a):
                    for _ in range(e):
                        v = apply(op, v, p)
                return v
            parent_fan = fan(n, ev, p, counts)
            child_fan = fan(child_n, child_ev, p, counts)
            my, mz = (len(Frame(op, p).basis) for op in M[1:])
            ground = tuple(hyperbolic_cross(n, 3))
            witnesses = [tuple(sum(a[i] for a in poly) for i in range(3))
                         for poly in list(parent_fan.values())+list(child_fan.values())]
            for w in witnesses:
                order = sorted(ground, key=lambda a: (dot(a, w), a))
                B = frozenset(field_basis(order, ev, p)[0])
                C = frozenset(field_basis(order, child_ev, p)[0])
                assert len(B) == n and len(C) == child_n and C <= B
                E = B-C
                for a in E:
                    assert all(tuple(a[j]+(j == i) for j in range(3)) not in B for i in (1, 2))
                assert sum(a[1] for a in E) <= my and sum(a[2] for a in E) <= mz
                assert len(E) == len(L) and sum(a[0] for a in E) <= len(L)*(n-1)
                counts['nonrow_moving_actual_annihilator_profile_checks'] += 1
            rows.append({'characteristic': p, 'forced_family_parameter': parameter,
                         'parent_length': n, 'actual_annihilator_dimension': len(L),
                         'quotient_length': child_n, 'parent_initial_ideals': len(parent_fan),
                         'quotient_initial_ideals': len(child_fan), 'branch_dimensions': [my, mz]})
    assert any(row['quotient_initial_ideals'] > 1 for row in rows)
    return rows


def main():
    counts, rows = Counter(), []
    for p in (0, 2, 3, 5):
        for drops in ((1, 3), (1, 2, 3), (2, 1, 2)):
            for core, shear in ((False, True), (True, True)):
                rows.append(inspect(drops, (0, 1), p, core, shear,
                                    p in (0, 2) and core, counts))
        rows.append(inspect((1, 2), (1, 1, 1), p, True, True, p == 2, counts))
    asymmetric = []
    for p in (0, 2, 3):
        labels = tuple([(i, 0, 0) for i in range(4)]
                       +[(i, j+1, 0) for j, width in enumerate((4, 3)) for i in range(width)]
                       +[(i, 0, j+1) for j, width in enumerate((4, 2)) for i in range(width)])
        index = {a: i for i, a in enumerate(labels)}
        M = tuple(tuple({index[b]: 1} if (b := tuple(a[j]+(j == k) for j in range(3)))
                         in index else {} for a in labels) for k in range(3))
        E = frozenset(a for a in labels if all(tuple(a[j]+(j == k) for j in range(3))
                                               not in index for k in (1, 2)))
        C = frozenset(labels)-E
        assert sum(a[1] for a in E) == 7 and sum(a[2] for a in E) == 6
        T = ({index[1, 2, 0]: 1, index[2, 0, 1]: -1},)
        n, ops, unit = quotient_by_relations(M, {index[0, 0, 0]: 1}, T, p)
        assert n == 15
        @lru_cache(None)
        def ev(a):
            v = unit
            for op, e in zip(ops, a):
                for _ in range(e):
                    v = apply(op, v, p)
            return v
        for b in range(3):
            block = [a for a in labels if (not a[1] and not a[2])]
            if b:
                block = [a for a in labels if a[b]]
            assert len(field_basis(block, ev, p)[0]) == len(block)
            counts['asymmetric_actual_endpoint_injections'] += 1
        states = fan(n, ev, p, counts)
        for B in states:
            assert C <= B <= frozenset(labels)
            assert sum(a[1] for a in B-C) <= 7 and sum(a[2] for a in B-C) <= 6
            counts['asymmetric_propagating_branch_entire_profiles'] += 1
        asymmetric.append({'characteristic': p, 'branch_row_lengths': [[4, 3], [4, 2]],
                           'envelope_length': len(labels), 'quotient_length': n,
                           'relation': 'x*y^2-x^2*z and all its x multiples',
                           'complete_initial_ideals': len(states)})
    # Exhaust all small monomial diagrams, including repeated row widths.
    for n in range(1, 16):
        for shape in partitions(n):
            D = {(i, j) for j, width in enumerate(shape) for i in range(width)}
            S = {a for a in D if (a[0], a[1]+1) not in D}
            assert sum(j+1 for i, j in S) == n
            counts['independent_all_staircase_telescoping_checks'] += 1
    for p in (0, 2, 5):
        for m in range(1, 6):
            n, ev = pascal_ray_model(m, p)
            if p == 0 or (m in (1, 2, 4) and p == 2) or (m in (1, 5) and p == 5):
                actual = fan(n, ev, p, counts)
                predicted = pascal_ray_predicted(m, collapse=p != 0)
                compare_fans(actual, predicted, counts)
                counts['existing_sharp_pascal_family_full_fans'] += 1
    moving = moving_annihilator_controls(counts)
    paths = {Path(__file__).resolve(), ROOT/'docs/row_branch_fans.md'}
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name:
            path = Path(name).resolve()
            if path.suffix == '.py' and path.parent == ROOT/'scripts':
                paths.add(path)
    out = {'command': 'python3 -B scripts/run_row_branch_checks.py', 'all_checks_passed': True,
           'scope': 'Uniform row-polynomial branch class; a classical lattice-volume consequence '
                    'of a fixed joint-annihilator quotient. No unrestricted result or novelty claim.',
           'cases': rows, 'asymmetric_branches': asymmetric,
           'general_moving_annihilator_quotients': moving, 'counters': dict(counts),
           'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in sorted(paths)},
           'unrestricted_G3_improved': False, 'polynomial_lower_separation_proved': False,
           'novelty_established': False, 'goal_complete': False}
    (ROOT/'results/research/row_branch_checks.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(dict(counts), sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
