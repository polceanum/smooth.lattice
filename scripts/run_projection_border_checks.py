#!/usr/bin/env python3
"""Exact coordinate quotients, border lifts and all joint fan patches."""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import hashlib
import json
import sys

from run_forced_socle_fan_checks import model as forced_model
from run_projection_tag_checks import fan, point_model
from run_signed_face_checks import field_basis, field_reduce, normals
from run_structural_complexity_checks import area2, clip, dot
from run_univariate_branch_gluing_checks import Frame, apply, clean, pascal_ray_model

ROOT = Path(__file__).resolve().parents[1]


def report_json(value, depth=0):
    """Keep each complete joint-cone certificate on one reviewable line."""
    if isinstance(value, dict):
        if 'x_y_projection_basis' in value:
            return json.dumps(value, separators=(',', ':'))
        pad = '  '*(depth+1)
        return '{\n'+',\n'.join(pad+json.dumps(k)+': '+report_json(v, depth+1)
                               for k, v in value.items())+'\n'+'  '*depth+'}'
    if isinstance(value, (list, tuple)) and any(isinstance(v, (dict, list, tuple)) for v in value):
        pad = '  '*(depth+1)
        return '[\n'+',\n'.join(pad+report_json(v, depth+1) for v in value)+'\n'+'  '*depth+']'
    return json.dumps(value)


def operators(reference, ev):
    return tuple(tuple(ev(tuple(a[j]+(i == j) for j in range(3)))
                       for a in reference) for i in range(3))


def evaluator(M, unit, p):
    @lru_cache(None)
    def ev(a):
        v = unit
        for op, exponent in zip(M, a):
            for _ in range(exponent):
                v = apply(op, v, p)
        return v
    return ev


def quotient(M, unit, axis, p, counts):
    """Actual quotient by the entire coordinate-generated algebra ideal."""
    n = len(M[0])
    K = Frame(M[axis], p)
    for v in K.basis:
        for op in M:
            K.coords(apply(op, v, p))
    F = Frame(K.basis+tuple({i: 1} for i in range(n)), p)
    k = len(K.basis)
    def project(v):
        return clean({i-k: c for i, c in F.coords(v).items() if i >= k}, p)
    Q = tuple(tuple(project(apply(op, v, p)) for v in F.basis[k:]) for op in M)
    u = project(unit)
    assert all(not v for v in Q[axis])
    for A, B in combinations(Q, 2):
        assert all(apply(A, B[j], p) == apply(B, A[j], p) for j in range(n-k))
    counts['independent_actual_coordinate_quotients'] += 1
    return n-k, evaluator(Q, u, p), project


def axis_degree(X, unit, p):
    powers, v = [], unit
    while len(Frame(powers+[v], p).basis) > len(powers):
        powers.append(v)
        v = apply(X, v, p)
    return len(powers)


def inspect(name, reference, ev, p, counts):
    n = len(reference)
    M = operators(reference, ev)
    unit = ev((0, 0, 0))
    assert all(ev(a) == {i: 1} for i, a in enumerate(reference))
    for U, V in combinations(M, 2):
        for j in range(n):
            assert apply(U, V[j], p) == apply(V, U[j], p)
            counts['independent_parent_commuting_columns'] += 1
    assert all(not apply(M[1], v, p) for v in M[2])
    assert len(field_basis(reference, ev, p)[0]) == n
    ny, ey, py = quotient(M, unit, 2, p, counts)  # A/zA in x,y
    nz, ez, pz = quotient(M, unit, 1, p, counts)  # A/yA in x,z
    for a in reference:
        assert py(ev(a)) == ey(a) and pz(ev(a)) == ez(a)
        counts['independent_quotient_column_maps'] += 1
    parent = fan(n, ev, p, counts)
    Y = fan(ny, ey, p, counts) if ny else {frozenset(): next(iter(parent.values()))}
    Z = fan(nz, ez, p, counts) if nz else {frozenset(): next(iter(parent.values()))}
    # All supplied examples have nonzero projections. Zero cases are tested
    # separately by replacing the zero fan with the entire weight simplex.
    from fractions import Fraction as Q
    triangle = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
    if not ny:
        Y = {frozenset(): triangle}
    if not nz:
        Z = {frozenset(): triangle}
    my, mz = (len(Frame(op, p).basis) for op in M[1:])
    total = len(Frame(M[1]+M[2], p).basis)
    h, d = my+mz-total, n-total
    a = axis_degree(M[0], unit, p)
    axis = {(i, 0, 0) for i in range(a)}
    joint_area, patch_area, patches, cells = 0, 0, 0, []
    for BY in Y:
        qy = normals(tuple(sorted(BY)), ey, p)[0] if ny else ()
        assert all(not b[2] for b in BY)
        for BZ in Z:
            qz = normals(tuple(sorted(BZ)), ez, p)[0] if nz else ()
            assert all(not b[1] for b in BZ)
            poly = list(triangle)
            for q in qy+qz:
                poly = clip(poly, q, 1)
            assert poly and area2(poly) > 0, 'every pair of planar intervals gives a joint cone'
            joint_area += area2(poly)
            C = BY | BZ
            EY = {tuple(b[i]+(i == 1) for i in range(3)) for b in BY}-C
            EZ = {tuple(b[i]+(i == 2) for i in range(3)) for b in BZ}-C
            E = (axis | EY | EZ)-C
            k = n-len(C)
            r = len(BY & BZ)
            assert all(not b[1] and not b[2] for b in BY & BZ)
            assert k == h+r-d >= 0
            assert sum(b[1] for b in E) == ny and sum(b[2] for b in E) == nz
            assert len(field_basis(tuple(sorted(C)), ev, p)[0]) == len(C)
            assert len(field_basis(tuple(sorted(C | E)), ev, p)[0]) == n
            for b in C | E:
                for i, op in enumerate(M):
                    assert apply(op, ev(b), p) == ev(tuple(b[j]+(i == j) for j in range(3)))
                    counts['independent_candidate_monomial_successors'] += 1
            local = []
            for B, parent_poly in parent.items():
                qparent = normals(tuple(sorted(B)), ev, p)[0]
                patch = list(poly)
                for q in qparent:
                    patch = clip(patch, q, 1)
                    if not patch:
                        break
                if not patch or area2(patch) <= 0:
                    continue
                patch_area += area2(patch)
                assert C <= B <= C | E
                variable = B-C
                assert len(variable) == k
                assert sum(b[0] for b in variable) <= k*(n-1)
                assert sum(b[1] for b in variable) <= ny
                assert sum(b[2] for b in variable) <= nz
                w = tuple(sum(v[i] for v in patch) for i in range(3))
                selected, _ = field_basis(tuple(sorted(C))+
                    tuple(sorted(E, key=lambda b: (dot(b, w), b))), ev, p)
                assert frozenset(selected) == B
                # Independent rank test of the chosen residual columns modulo
                # the compulsory span, without replaying merged greedy.
                _, pivots = field_basis(tuple(sorted(C)), ev, p)
                residual = {b: field_reduce(ev(b), pivots, p)[0] for b in E}
                assert len(Frame(tuple(residual[b] for b in variable), p).basis) == k
                local.append({'basis': sorted(B), 'strict_weight': [str(t) for t in w],
                              'normalized_double_area': str(area2(patch))})
                patches += 1
                counts['independent_complete_patch_border_rank_and_moment_checks'] += 1
            assert local and (k or len(local) == 1)
            cells.append({'x_y_projection_basis': sorted(BY), 'x_z_projection_basis': sorted(BZ),
                          'residual_rank': k, 'axis_intersection_length': r,
                          'all_residual_candidates': sorted(E),
                          'actual_parent_patches': local})
    assert joint_area == patch_area == 1
    counts['complete_joint_and_parent_overlay_coverages'] += 1
    return {'name': name, 'characteristic': p, 'length': n,
            'branch_intersection_dimension': h, 'axis_quotient_length': d,
            'x_axis_minimal_polynomial_degree': a,
            'coordinate_quotient_lengths': [ny, nz],
            'coordinate_quotient_initial_ideals': [len(Y), len(Z)],
            'parent_initial_ideals': len(parent), 'all_joint_cones': len(cells),
            'complete_parent_joint_patches': patches, 'joint_cones': cells}


def main():
    counts, rows = Counter(), []
    for p in (0, 2, 5):
        for parameter, kind in ((3, 'H'), (4, 'A'), (4, 'C')):
            reference, ev = forced_model(parameter, kind, p)
            rows.append(inspect('nonrow_forced_'+kind, reference, ev, p, counts))
        for m in (2, 3, 4):
            n, ev = pascal_ray_model(m, p)
            reference = tuple((i, 0, 0) for i in range(n))
            row = inspect('pascal_high_overlap', reference, ev, p, counts)
            assert row['coordinate_quotient_initial_ideals'] == [1, 1]
            if not p:
                assert row['parent_initial_ideals'] == (m+1)*(m+2)//2
            rows.append(row)
    # An actual zero-overlap direct product of two cyclic monomial-curve
    # jets. Each actual planar quotient has two initial ideals, independently
    # of exponent and jet length (translations preserve global initials).
    for p in (0, 3):
        R, S, k, ell = 5, 4, 2, 3
        n = R+S
        from run_minimal_branch_overlap_checks import mul, power, rem
        P, Q = (0,)*R+(1,), power((-1, 1), S, p)
        fy = tuple([1]+[0]*(k-1)+[1])
        gz = list(power((-1, 1), ell, p)); gz[0] += 1
        @lru_cache(None)
        def ev(a):
            x, y, z = a
            if y and z:
                return {}
            left = {} if z else {i: c for i, c in enumerate(rem(mul((0,)*x+(1,), power(fy, y, p), p), P, p)) if c}
            right = {} if y else {R+i: c for i, c in enumerate(rem(mul((0,)*x+(1,), power(tuple(gz), z, p), p), Q, p)) if c}
            return clean({**left, **right}, p)
        basis, _ = field_basis(tuple((i, 0, 0) for i in range(n)), ev, p)
        assert len(basis) == n
        # Change to the actual x-power reference basis, independently of the
        # component coordinate implementation used above.
        F = Frame(tuple(ev(b) for b in basis), p)
        original = ev
        ev = lru_cache(None)(lambda a: F.coords(original(a)))
        row = inspect('two_monomial_curve_jets', basis, ev, p, counts)
        assert row['branch_intersection_dimension'] == 0
        assert row['coordinate_quotient_initial_ideals'] == [2, 2]
        rows.append(row)
    # Zero coordinate quotients are permitted by the theorem.
    for point in ((1, 0, 1), (0, 1, 0)):
        n, ev = point_model((point,), 0)
        rows.append(inspect('zero_coordinate_quotient', ((0, 0, 0),), ev, 0, counts))
    for p in (0, 7):
        points = ((1, 2, 0), (2, 3, 0), (3, 0, 4), (4, 0, 2), (5, 0, 0))
        n, raw = point_model(points, p)
        reference = tuple((i, 0, 0) for i in range(n))
        F = Frame(tuple(raw(a) for a in reference), p)
        ev = lru_cache(None)(lambda a: F.coords(raw(a)))
        row = inspect('radical_two_planes_with_axis_point', reference, ev, p, counts)
        assert row['branch_intersection_dimension'] == 0 and row['axis_quotient_length'] == 1
        rows.append(row)
    paths = {Path(__file__).resolve(), ROOT/'docs/projection_border_fans.md'}
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name:
            path = Path(name).resolve()
            if path.suffix == '.py' and path.parent == ROOT/'scripts':
                paths.add(path)
    out = {'command': 'python3 -B scripts/run_projection_border_checks.py', 'all_checks_passed': True,
           'scope': 'All-field uniform conditional count for yz=0, not an unrestricted G3 theorem. '
                    'Border multiplication, forced-column contraction and lattice volume are classical inputs.',
           'cases': rows, 'counters': dict(counts),
           'source_sha256': {str(f.relative_to(ROOT)): hashlib.sha256(f.read_bytes()).hexdigest()
                             for f in sorted(paths)},
           'unrestricted_G3_improved': False, 'polynomial_lower_separation_proved': False,
           'novelty_established': False, 'goal_complete': False}
    serialized = report_json(out)+'\n'
    assert json.loads(serialized) == json.loads(json.dumps(out))
    (ROOT/'results/research/projection_border_checks.json').write_text(serialized)
    print(json.dumps(dict(counts), sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
