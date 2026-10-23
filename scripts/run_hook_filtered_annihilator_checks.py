#!/usr/bin/env python3
"""Filtered full-kernel checks and original coherent hook exchanges."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_hook_bridge_checks import model, actual, normals, COUNTS as REPLAY
from run_hook_transverse_prefix_checks import mul, clean, degree, wall_weight
from run_hook_monomial_perturbation_checks import reflected, row_functional
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import dot

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def scalar(c, p):
    return c % p if p else Q(c)


def dense_nullspace(columns, labels, p):
    """Independent elimination, with input coefficients as free variables."""
    rows = sorted(set().union(*(set(columns[v]) for v in labels)))
    matrix = [[scalar(columns[v].get(r, 0), p) for v in labels] for r in rows]
    pivots = []; top = 0
    for col in range(len(labels)):
        pick = next((i for i in range(top, len(matrix)) if matrix[i][col]), None)
        if pick is None:
            continue
        matrix[top], matrix[pick] = matrix[pick], matrix[top]
        inv = pow(int(matrix[top][col]), -1, p) if p else 1 / matrix[top][col]
        matrix[top] = [scalar(c * inv, p) for c in matrix[top]]
        for i in range(len(matrix)):
            if i == top or not matrix[i][col]:
                continue
            c = matrix[i][col]
            matrix[i] = [scalar(a-c*b, p) for a, b in zip(matrix[i], matrix[top])]
        pivots.append(col); top += 1
        if top == len(matrix):
            break
    free = [i for i in range(len(labels)) if i not in pivots]
    kernel = []
    for f in free:
        h = {labels[f]: scalar(1, p)}
        h.update({labels[c]: scalar(-matrix[i][f], p) for i, c in enumerate(pivots)
                  if matrix[i][f]})
        kernel.append(h)
    return len(pivots), kernel


def combination(columns, h, p):
    out = Counter()
    for v, c in h.items():
        for row, value in columns[v].items():
            out[row] += c * value
    return clean(out, p)


def bases(N, A, t, u, p):
    L = A-t; cost = lambda v: v[0]*u+v[1]
    full = tuple(sorted(product(range(A), range(t)), key=lambda v: (cost(v), v)))
    pure = tuple(v for v in full if v[0] >= t)
    small = tuple(v for v in full if v[0] < L)
    values = {v: mul(N, {v: 1}, A, t, p) for v in full}
    truncated = {v: {r: c for r, c in values[v].items() if r[0] < L} for v in small}
    B, _ = field_basis(full, values.__getitem__, p)
    BP, _ = field_basis(pure, values.__getitem__, p)
    plus, piv = field_basis(small, values.__getitem__, p)
    minus, _ = field_basis(small, truncated.__getitem__, p)
    assert set(B) & set(pure) <= set(BP)
    assert set(minus) <= set(plus)
    assert set(BP) == {(i+t, j) for i, j in minus}
    SP = set(BP)-set(B); SD = set(plus)-set(minus)
    assert len(SP) <= t*t-(len(B)-len(BP))
    assert len(SD) == len(plus)-len(minus) <= t*t
    return full, small, values, truncated, B, BP, plus, minus, piv, SP, SD


def project(h, values, piv, p):
    rem, one = field_reduce(combination(values, h, p), piv, p)
    assert not rem
    zero = Counter(h)
    for v, c in one.items():
        zero[v] -= c
    return clean(zero, p), clean(one, p)


def interval_weights(A, t):
    cuts = sorted({Q(j, i) for i in range(1, A) for j in range(1, t)})
    return [cuts[0]/2] + [(a+b)/2 for a, b in zip(cuts, cuts[1:])] + [cuts[-1]+1]


def operator_controls():
    rng = Random(2026102317); records = []
    cases = [(6, 2, 0), (7, 3, 2), (8, 3, 3), (9, 4, 101)]
    for A, t, p in cases:
        N = clean({(i, j): rng.randrange(p or 11)-2
                   for i, j in product(range(A), range(1, t))}, p)
        operators = [N, {(1, t-1): 1, (A-t+1, 1): 1}, {}]
        if t >= 4:
            operators.append({(1, t-1): 1, (3, t-2): 1, (A-3, 1): 1})
        for number, N in enumerate(operators):
            controls = []
            for u in interval_weights(A, t):
                full, small, vals, trunc, B, BP, plus, minus, piv, SP, SD = bases(N, A, t, u, p)
                rank, kernel = dense_nullspace(trunc, small, p)
                assert rank == len(minus) and len(kernel) == len(small)-rank
                rankA, _ = dense_nullspace(vals, full, p)
                assert rankA == len(B)
                cost = lambda v: v[0]*u+v[1]
                for v in small:
                    zero, one = project({v: 1}, vals, piv, p)
                    assert all(cost(q) <= cost(v) for q in (*zero, *one))
                    assert not combination(vals, zero, p)
                    COUNTS['triangular_input_column_projections'] += 1
                for h in kernel:
                    assert not combination(trunc, h, p)
                    zero, one = project(h, vals, piv, p)
                    cap = max(map(cost, h))
                    assert all(cost(q) <= cap for q in (*zero, *one))
                    assert not combination(vals, zero, p)
                    assert not combination(trunc, zero, p)
                    assert not combination(trunc, one, p)
                    if len(zero) > 1 and any(vals[v] for v in zero):
                        COUNTS['nonmonomial_full_kernel_projection_components'] += 1
                    if one:
                        assert max(one, key=cost) in SD
                    COUNTS['independent_dense_kernel_projection_vectors'] += 1
                # Every prefix is checked, rather than just testing individual
                # monomials or a random sample of kernel vectors.
                for size in range(1, len(small)+1):
                    prefix = small[:size]
                    r, _ = dense_nullspace(trunc, prefix, p)
                    assert size-r == sum(v not in minus for v in prefix)
                    COUNTS['all_vector_filtered_prefix_kernel_dimensions'] += 1
                controls.append({'u': str(u), 'full_rank': len(B), 'pure_rank': len(BP),
                                 'restricted_rank': len(plus), 'truncated_rank': len(minus),
                                 'primal_shell': sorted(SP), 'dual_shell': sorted(SD)})
                COUNTS['all_open_internal_order_intervals'] += 1
            records.append({'A': A, 't': t, 'field_characteristic': p,
                            'operator': [[*v, int(c)] for v, c in sorted(N.items())],
                            'operator_number': number, 'interval_controls': controls})
    assert COUNTS['nonmonomial_full_kernel_projection_components']
    return records


def assemble(t, a, units, p):
    data = model(t, a); U, V = units
    f = {(t+i, j): c for (i, j), c in U.items()}
    f.update({(j, t+i): c for (i, j), c in V.items()})
    co = [f.get((i, j), 0) if g == 1 else 0 for g, i, j in data[2]]
    return data, co, actual(data, co, p, replay=True)


def eval_row(ll, v, cols, p):
    return scalar(sum(ll.get(i, 0)*c for i, c in cols[v].items()), p)


def original_controls():
    records = []
    cases = [(3, 1, 12, 1, 0), (3, 1, 12, 1, 2), (3, 1, 12, 2, 101),
             (5, 2, 20, 1, 2), (5, 2, 20, 2, 101), (7, 3, 24, 3, 3),
             (7, 2, 28, 2, 101)]
    for t, e0, A, mu, p in cases:
        d = t-2*e0; a = A+t; nu = A-t+mu; L = A-t
        U = {(0, 0): 1, (mu, d+e0): 1, (nu, d): 1}
        V = {(0, 0): 1, (A-1-e0, 1): 1}
        units = (U, V); Ns = [{v: c for v, c in W.items() if v != (0, 0)} for W in units]
        data, co, cols = assemble(t, a, units, p)
        H = {(A-t, 0): 1, (0, e0): scalar(-1, p)}
        assert not mul(Ns[0], H, A, t, p)
        assert all(mul(Ns[0], {v: 1}, A, t, p) for v in H)
        labels = tuple(product(range(A), range(t)))
        nv = {v: mul(Ns[0], {v: 1}, A, t, p) for v in labels}
        rank, _ = dense_nullspace(nv, labels, p)
        assert rank == e0*(A-mu)
        assert not mul(Ns[0], {(A-mu, 0): 1}, A, t, p)
        COUNTS['nonmonomial_full_annihilator_rank_checks'] += 1
        E = (0, A, t); F = (1, 1, a-1); q = (1, -(A-1), A-1)
        w = wall_weight((Q(e0), 1+Q(e0, A-1), Q(1)),
                        (Q(t), 1+Q(t, A-1), Q(1)), cols, q)
        lam, u, _ = w; tau = dot(w, E)
        assert e0 < lam < t and tau == dot(w, F)
        assert tau-lam > max((2*t-1)*u+t-1, (t-1)*u+2*t-1)
        assert len({i*u+j for i, j in product(range(A), range(t))}) == A*t
        cheap = tuple(v for v in cols if dot(w, v) < tau)
        bridge = (0, t, t+e0); gy = (1, A, 0); yl = (1, t, e0); zl = (1, 0, t+e0)
        relation = {E: 1, F: 1, bridge: -1, gy: -1, yl: 1, zl: 1}
        assert not combination(cols, relation, p)
        assert all(dot(w, v) < tau for v in (bridge, gy, yl, zl))
        idx = {v: i for i, v in enumerate(data[0])}
        ell = {idx[A+mu, d+e0]: 1, idx[1, a-1]: scalar(-1, p)}
        assert all(not eval_row(ell, v, cols, p) for v in cheap)
        assert eval_row(ell, E, cols, p) and eval_row(ell, F, cols, p)
        # This exact primal source comes from the independently verified
        # original quotient circuit, not from a kernel chosen to fit a shell.
        assert max(H, key=lambda v: v[0]*u+v[1]) == (A-t, 0)
        BB = bases(Ns[0], A, t, u, p)
        assert (A-t, 0) in BB[-2]
        COUNTS['actual_original_primal_annihilator_shell_memberships'] += 1
        # Include a nonzero removable component only after checking its
        # original cheaper constraints and both exchanged evaluations.
        added = False
        if e0+1 < L:
            candidate = row_functional({(e0+1, 0): 1}, 1, data, A, t)
            if (all(not eval_row(candidate, v, cols, p) for v in cheap)
                    and not eval_row(candidate, E, cols, p)
                    and not eval_row(candidate, F, cols, p)):
                for row, c in candidate.items():
                    ell[row] = scalar(ell.get(row, 0)+c, p)
                added = True; COUNTS['validated_original_nullrow_additions'] += 1
        branches = []
        for arm, K in enumerate(reflected(ell, data, A, t)):
            assert K and all(i < L for i, j in K)
            if mul(Ns[arm], K, L, t, p):
                branches.append('nonscalar'); continue
            weight = u if not arm else 1/u
            BB = bases(Ns[arm], A, t, weight, p)
            zero, one = project(K, BB[2], BB[8], p)
            cap = degree(K, weight)
            assert (not zero or degree(zero, weight) <= cap)
            assert one and degree(one, weight) <= cap
            ll0 = row_functional(zero, arm, data, A, t)
            ll1 = row_functional(one, arm, data, A, t)
            own, other = (E, F) if not arm else (F, E)
            assert all(not eval_row(ll0, v, cols, p) for v in cheap)
            assert not eval_row(ll0, own, cols, p) and not eval_row(ll0, other, cols, p)
            assert eval_row(ll1, own, cols, p) == eval_row(ell, own, cols, p)
            lead = max(one, key=lambda v: v[0]*weight+v[1])
            assert lead in BB[-1]
            COUNTS['actual_original_dual_filtered_shell_memberships'] += 1
            if zero:
                COUNTS['nonzero_removable_full_kernel_components'] += 1
            branches.append('filtered-shell')
        values = sorted({dot(w, v) for v in cols})
        eps = min(lam/2, min(y-x for x, y in zip(values, values[1:]))/4)
        adjacent = []
        for sign in (-1, 1):
            wp = (lam+sign*eps, u, Q(1))
            B, _ = field_basis(sorted(cols, key=lambda v: (dot(wp, v), v)), cols.__getitem__, p)
            assert all(dot(wp, r) > 0 for r in normals(B, cols, p))
            adjacent.append(set(B)); COUNTS['actual_original_adjacent_strict_minima'] += 1
        assert F in adjacent[0] and E not in adjacent[0]
        assert E in adjacent[1] and F not in adjacent[1]
        COUNTS['actual_original_coherent_mixed_corner_walls'] += 1
        records.append({'t': t, 'e': e0, 'd': d, 'A': A, 'mu': mu, 'nu': nu,
                        'field_characteristic': p, 'weight': list(map(str, w)),
                        'lower': E, 'upper': F, 'full_off_axis_rank': rank,
                        'coefficients': list(map(int, co)), 'dual_branches': branches,
                        'independently_validated_nullrow_added': added})
    assert COUNTS['nonzero_removable_full_kernel_components']
    return records




def same_layer_controls():
    """Nonmonomial truncated kernels on original upper/upper walls."""
    records = []
    for t, e0, A, mu, p in [(3, 1, 10, 1, 0), (3, 1, 10, 1, 2),
                             (5, 2, 17, 1, 101), (7, 2, 23, 2, 3)]:
        d = t-2*e0; a = A+t; L = A-t; h = A-2*t; nu = h+mu
        r = 5*t-A-mu-1; dz = 2*t-1-r
        assert t <= r < 2*t-1 and 0 < dz < t
        u = Q(A, 2*t); assert h*u > e0
        U = {(0, 0): 1, (mu, d+e0): 1, (nu, d): 1}
        V = {(0, 0): 1, (A-2*e0, dz): 1}
        Ns = [{v: c for v, c in W.items() if v != (0, 0)} for W in (U, V)]
        data, co, cols = assemble(t, a, (U, V), p)
        E = (1, 3*t-1, t-1); F = (1, t-1, a-1)
        q = (0, -2*t, A)
        w = wall_weight((Q(1, 1009), u, Q(1)),
                        (t-Q(1, 1009), u, Q(1)), cols, q)
        lam, u, _ = w; tau = dot(w, E)
        assert 0 < lam < t and tau == dot(w, F)
        assert tau-lam > max((2*t-1)*u+t-1, (t-1)*u+2*t-1)
        assert len({i*u+j for i, j in product(range(A), range(t))}) == A*t
        bridge = (0, r, t+2*e0-1)
        iy = (1, r, 2*e0-1); iz = (1, r-t, t+2*e0-1)
        assert not combination(cols, {E: 1, F: 1, bridge: -1, iy: 1, iz: 1}, p)
        assert all(dot(w, v) < tau for v in (bridge, iy, iz))
        KY = {(h, 0): 1, (0, e0): scalar(-1, p)}
        KZ = {(0, 0): scalar(-1, p)}
        assert not mul(Ns[0], KY, L, t, p)
        assert mul(Ns[0], KY, A, t, p) == {(mu+2*h, d): scalar(1, p)}
        assert all(mul(Ns[0], {v: 1}, L, t, p) for v in KY)
        ell = row_functional(KY, 0, data, A, t)
        ell.update(row_functional(KZ, 1, data, A, t))
        cheap = tuple(v for v in cols if dot(w, v) < tau)
        assert all(not eval_row(ell, v, cols, p) for v in cheap)
        assert eval_row(ell, E, cols, p) and eval_row(ell, F, cols, p)
        branches = []
        for arm, K in enumerate((KY, KZ)):
            weight = u if not arm else 1/u
            BB = bases(Ns[arm], A, t, weight, p)
            zero, one = project(K, BB[2], BB[8], p)
            cap = degree(K, weight)
            assert not zero or degree(zero, weight) <= cap
            assert one and degree(one, weight) <= cap
            own, other = (E, F) if not arm else (F, E)
            ll0 = row_functional(zero, arm, data, A, t)
            ll1 = row_functional(one, arm, data, A, t)
            assert all(not eval_row(ll0, v, cols, p) for v in cheap)
            assert not eval_row(ll0, own, cols, p) and not eval_row(ll0, other, cols, p)
            assert eval_row(ll1, own, cols, p) == eval_row(ell, own, cols, p)
            assert max(one, key=lambda v: v[0]*weight+v[1]) in BB[-1]
            COUNTS['actual_original_dual_filtered_shell_memberships'] += 1
            if not arm:
                assert any(mul(Ns[arm], {v: 1}, L, t, p) for v in one)
                COUNTS['actual_original_nonmonomial_dual_exceptions'] += 1
            branches.append({'arm': 'Y' if not arm else 'Z',
                             'retained_terms': len(one), 'removed_terms': len(zero)})
        values = sorted({dot(w, v) for v in cols})
        eps = min(u/2, min(y-x for x, y in zip(values, values[1:]))/(8*a))
        adjacent = []
        for sign in (-1, 1):
            wp = (lam, u+sign*eps, Q(1))
            B, _ = field_basis(sorted(cols, key=lambda v: (dot(wp, v), v)), cols.__getitem__, p)
            assert all(dot(wp, v) > 0 for v in normals(B, cols, p))
            adjacent.append(set(B)); COUNTS['actual_original_adjacent_strict_minima'] += 1
        assert E in adjacent[0] and F not in adjacent[0]
        assert F in adjacent[1] and E not in adjacent[1]
        COUNTS['actual_original_same_layer_coherent_walls'] += 1
        records.append({'t': t, 'e': e0, 'd': d, 'A': A, 'mu': mu, 'nu': nu,
                        'field_characteristic': p, 'weight': list(map(str, w)),
                        'Y_upper': E, 'Z_upper': F, 'dual_branches': branches,
                        'coefficients': list(map(int, co))})
    return records

def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    operators = operator_controls()
    print('all-interval filtered full-kernel controls passed', flush=True)
    originals = original_controls()
    same_layer = same_layer_controls()
    print('original commuting quotient and mixed-corner walls passed', flush=True)
    files = {'scripts/run_hook_filtered_annihilator_checks.py', 'docs/hook_filtered_annihilators.md',
             'docs/hook_transverse_bases.md', 'docs/hook_transverse_prefix.md',
             'docs/hook_bridge_reduction.md', 'docs/planar_dyadic_bound.md'}
    pending = ['scripts/run_hook_filtered_annihilator_checks.py']
    while pending:
        f = pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names = [node.module] if isinstance(node, ast.ImportFrom) else [v.name for v in node.names] if isinstance(node, ast.Import) else []
            for name in names:
                target = 'scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:
                    files.add(target); pending.append(target)
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_hook_filtered_annihilator_checks.py',
              'scope': 'Exact filtered annihilator reduction for arbitrary mixed constant-axis central hook operators; no aggregate fan bound or unrestricted exponent improvement.',
              'counts': dict(sorted(COUNTS.items())), 'original_quotient_replay_counts': dict(REPLAY),
              'operator_controls': operators, 'original_coherent_controls': originals,
              'same_layer_coherent_controls': same_layer,
              'limits': ['Outer images are zero.', 'Both same-layer and cross-layer original-wall controls satisfy no internal ties and cheap overlap.',
                         'The small-upper-threshold same-layer branch is not covered by this reduction.',
                         'An O(t^2) pointwise difference does not count simultaneous or recurrent wall pieces.',
                         'Exact priority and the original unrestricted endpoint remain open.'],
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_filtered_annihilator_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'counts': report['counts'], 'source_hashes': len(files)}))


if __name__ == '__main__':
    main()
