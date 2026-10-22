#!/usr/bin/env python3
"""Exact module minima, rigid equality models and independent full fan checks."""
import ast
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import hashlib
import json

from run_univariate_branch_gluing_checks import (
    Frame, apply, clean, cyclic_span, dense_kernel, dense_rank, jsonable, linear,
    pascal_ray_model, pascal_ray_predicted, polynomial,
)
from run_projection_tag_checks import compare_fans, fan
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import area2, clip, dot

ROOT = Path(__file__).resolve().parents[1]
TRIANGLE = [tuple(Q(i == j) for i in range(3)) for j in range(3)]


def partitions(n, cap=None):
    if not n:
        yield ()
    else:
        for first in range(min(n, cap or n), 0, -1):
            for rest in partitions(n-first, first):
                yield (first,)+rest


def mul(a, b, p):
    out = [0]*(len(a)+len(b)-1)
    for i, c in enumerate(a):
        for j, d in enumerate(b):
            out[i+j] += c*d
    return tuple(c % p if p else Q(c) for c in out)


def power(a, k, p):
    out = (1,)
    for _ in range(k):
        out = mul(out, a, p)
    return out


def rem(a, b, p):
    a = list(a)
    assert b[-1] == 1
    for i in range(len(a)-1, len(b)-2, -1):
        lead = a[i]
        for j, c in enumerate(b):
            a[i-len(b)+1+j] -= lead*c
        if p:
            a = [c % p for c in a]
    return tuple(a[:len(b)-1])


def module_data(ds, p):
    assert ds and all(len(d) > 1 and d[-1] == 1 for d in ds)
    assert all(not any(rem(a, b, p)) for a, b in zip(ds, ds[1:]))
    gs = [None]*len(ds)+[(1,)]
    for j in range(len(ds)-1, -1, -1):
        gs[j] = mul(ds[j], gs[j+1], p)
    h = sum(len(d)-1 for d in ds)
    phi = sum((j+1)*(len(d)-1) for j, d in enumerate(ds))
    return tuple(gs), h, phi


def envelope(ds, p):
    """Row-polynomial model of H_max; no quotient-fan formula is used."""
    gs, h, phi = module_data(ds, p)
    D, r = len(ds[0])-1, len(ds)
    core_poly = mul(gs[0], ds[0], p)
    labels = tuple((i, 0, 0) for i in range(h+D))
    labels += tuple((i, j+1 if branch == 1 else 0,
                     j+1 if branch == 2 else 0)
                    for branch in (1, 2) for j in range(r)
                    for i in range(len(gs[j])-1))
    index = {a: i for i, a in enumerate(labels)}

    def row(poly, j=0, branch=0):
        if branch and j >= r:
            return {}
        bound = gs[j] if branch else core_poly
        coeffs = rem(poly, bound, p)
        return clean({index[i, j+1 if branch == 1 else 0,
                             j+1 if branch == 2 else 0]: c
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
    wlabels = tuple((j, k) for j, d in enumerate(ds) for k in range(len(d)-1))
    WX = []
    for j, k in wlabels:
        reduced = rem((0,)*(k+1)+(1,), ds[j], p)
        WX.append(clean({wlabels.index((j, i)): c
                         for i, c in enumerate(reduced) if c}, p))
    WY = tuple(row((0,)*k+gs[j+1], j, 1) for j, k in wlabels)
    WZ = tuple(row((0,)*k+gs[j+1], j, 2) for j, k in wlabels)
    return dict(ds=ds, gs=gs, h=h, phi=phi, D=D, r=r, labels=labels,
                M=M, WY=WY, WZ=WZ, WX=tuple(WX), wlabels=wlabels,
                g0=row(gs[0]), unit={index[0, 0, 0]: 1})


def parameters(H, p):
    h, X, ds = H['h'], H['WX'], H['ds']
    identity = tuple({i: 1} for i in range(h))
    choices = [('identity', identity, {})]
    u = {i: 1 for i, (_, k) in enumerate(H['wlabels']) if k == 0}
    choices.append(('nonzero_core', identity, u))
    unit = tuple(linear(((1, {i: 1}), (1, X[i])), p) for i in range(h))
    if dense_rank(unit, h, p) == h:
        choices.append(('unit_1_plus_x', unit, u))
    if len(ds) > 1:
        shear = list(identity)
        for i, (j, k) in enumerate(H['wlabels']):
            if j == 0:
                reduced = rem((0,)*k+(1,), ds[1], p)
                extra = {H['wlabels'].index((1, t)): c
                         for t, c in enumerate(reduced) if c}
                shear[i] = linear(((1, shear[i]), (1, extra)), p)
        choices.append(('off_diagonal_shear', tuple(shear), u))
    return choices


def relations(H, theta, u, p):
    h, X = H['h'], H['WX']
    assert dense_rank(theta, h, p) == h
    assert all(apply(theta, X[i], p) == apply(X, theta[i], p) for i in range(h))
    T = [linear(((1, H['WY'][i]), (-1, apply(H['WZ'], theta[i], p))), p)
         for i in range(h)]
    core = linear(((1, H['g0']), (-1, apply(H['WY'], u, p))), p)
    for _ in range(H['D']):
        T.append(core)
        core = apply(H['M'][0], core, p)
    assert not polynomial(H['M'][0], T[h], H['ds'][0], p)
    return tuple(T)


def quotient(H, theta, u, p, counts):
    T = relations(H, theta, u, p)
    N, k = len(H['labels']), H['h']+H['D']
    assert dense_rank(T, N, p) == k
    FT = Frame(T, p)
    for v in T:
        FT.coords(apply(H['M'][0], v, p))
        assert not apply(H['M'][1], v, p) and not apply(H['M'][2], v, p)
    full = Frame(T+tuple({i: 1} for i in range(N)), p)
    assert len(full.basis) == N

    def project(v):
        return clean({i-k: c for i, c in full.coords(v).items() if i >= k}, p)

    refs = full.basis[k:]
    M = tuple(tuple(project(apply(op, v, p)) for v in refs) for op in H['M'])
    unit = project(H['unit'])
    n = N-k
    assert n == 2*H['phi']
    for U, V in combinations(M, 2):
        for i in range(n):
            assert apply(U, V[i], p) == apply(V, U[i], p)
            counts['commuting_coordinate_columns'] += 1

    @lru_cache(None)
    def ev(a):
        v = unit
        for op, e in zip(M, a):
            for _ in range(e):
                v = apply(op, v, p)
        return v

    assert all(not apply(M[1], v, p) for v in M[2])
    assert dense_rank(tuple(ev(a) for a in H['labels']), n, p) == n
    counts['actual_cyclic_quotients'] += 1
    return n, M, unit, ev


def intersection(M, p):
    n = len(M[0])
    Y, Z = Frame(M[1], p).basis, Frame(M[2], p).basis
    kernel, _ = dense_kernel(Y+tuple(linear(((-1, v),), p) for v in Z), n, p)
    W = tuple(apply(Y, {i: c for i, c in v.items() if i < len(Y)}, p)
              for v in kernel)
    assert len(Frame(W, p).basis) == len(W)
    return Y, Z, W


def recover_partitions(X, W, factors, p, counts):
    """Recover primary exponents from ranks, not the construction's rows."""
    F, answer = Frame(W, p), []
    action = tuple(F.coords(apply(X, v, p)) for v in W)
    for f in factors:
        op = tuple(polynomial(action, {i: 1}, f, p) for i in range(len(W)))
        images = tuple({i: 1} for i in range(len(W)))
        nullities = [0]
        for _ in range(len(W)+1):
            images = tuple(apply(op, v, p) for v in images)
            rank = dense_rank(images, len(W), p)
            assert rank == len(Frame(images, p).basis)
            nullities.append(len(W)-rank)
            counts['independent_primary_power_rank_checks'] += 1
            if nullities[-1] == nullities[-2]:
                break
        else:
            raise AssertionError('Primary ranks did not stabilize.')
        deg = len(f)-1
        assert all((b-a) % deg == 0 for a, b in zip(nullities, nullities[1:]))
        heights = [(b-a)//deg for a, b in zip(nullities, nullities[1:])]
        assert all(a >= b for a, b in zip(heights, heights[1:]))
        shape = tuple(sum(c > j for c in heights) for j in range(heights[0]))
        answer.append((f, shape))
    assert sum((len(f)-1)*sum(shape) for f, shape in answer) == len(W)
    return tuple(answer)


def verify_model(name, H, theta, u, factors, p, counts, full_fan=False):
    n, M, unit, ev = quotient(H, theta, u, p, counts)
    Y, Z, W = intersection(M, p)
    assert len(Y) == len(Z) == H['phi'] and len(W) == H['h']
    assert all(not apply(M[1], v, p) and not apply(M[2], v, p) for v in W)
    primary = recover_partitions(M[0], W, factors, p, counts)
    recovered = [(1,)]*max(len(shape) for _, shape in primary)
    for f, shape in primary:
        for j, exponent in enumerate(shape):
            recovered[j] = mul(recovered[j], power(f, exponent, p), p)
    assert tuple(recovered) == tuple(tuple(c % p if p else Q(c) for c in d)
                                     for d in H['ds'])
    assert sum((len(f)-1)*sum((i+1)*v for i, v in enumerate(shape))
               for f, shape in primary) == H['phi']
    # Actual branch annihilators and kernels: this is independent of row ranks.
    for axis, branch in ((1, Y), (2, Z)):
        for j, g in enumerate(H['gs']):
            a = [0, 0, 0]; a[axis] = j+1
            assert not polynomial(M[0], ev(tuple(a)), g, p)
        assert len(branch)-dense_rank(tuple(apply(M[axis], v, p) for v in branch), n, p) == H['h']
    FW = Frame(W, p)
    FW.coords(polynomial(M[0], unit, H['gs'][0], p))
    assert not polynomial(M[0], unit, mul(H['gs'][0], H['ds'][0], p), p)

    C = {(i, 0, 0) for i in range(H['h'])}
    for j in range(H['r']-1):
        for i in range(len(H['gs'][j+1])-1):
            C.update(((i, j+1, 0), (i, 0, j+1)))
    C = tuple(sorted(C))
    E = tuple(sorted(set(H['labels'])-set(C)))
    assert len(C) == n-H['h'] and len(E) == 2*H['h']+H['D']
    assert dense_rank(tuple(ev(a) for a in C)+W, n, p) == n
    assert sum(a[1] for a in E) == sum(a[2] for a in E) == H['phi']
    assert max(a[0] for a in E) <= H['h']+H['D']-1
    row = dict(name=name, characteristic=p, invariant_factors=H['ds'],
               primary_partitions=primary, h=H['h'], r=H['r'], phi=H['phi'],
               length=n, theta_columns=theta, core_value=u,
               common_standard_set=C, residual_ground=E)
    if full_fan:
        actual = fan(n, ev, p, counts)
        predicted = residual_fan(C, E, H['h'], ev, p, counts)
        compare_fans(actual, predicted, counts)
        for Bset in actual:
            assert set(C) <= Bset <= set(H['labels'])
            residual = Bset-set(C)
            moment = tuple(sum(a[i] for a in residual) for i in range(3))
            assert moment[0] <= H['h']*(H['h']+H['D']-1)
            assert moment[1] <= H['phi'] and moment[2] <= H['phi']
        row['cones'] = len(actual)
        row['complete_standard_sets'] = [sorted(Bset) for Bset in sorted(actual, key=lambda s: sorted(s))]
        counts['independent_full_residual_compilations'] += 1
    return row


def residual_fan(C, E, h, ev, p, counts):
    """Enumerate all abstract matroid bases and all fundamental exchanges.

    This does not traverse algebra borders or assume the contraction is an
    algebra. Compare its complete rational polygons with the separate parent
    quotient border traversal.
    """
    _, pivots = field_basis(C, ev, p)
    columns = {a: field_reduce(ev(a), pivots, p)[0] for a in E}
    cells = {}
    for chosen in combinations(E, h):
        F = Frame(tuple(columns[a] for a in chosen), p)
        if len(F.basis) < h:
            continue
        normals = set()
        for a in E:
            if a in chosen:
                continue
            for j in F.coords(columns[a]):
                normals.add(tuple(a[k]-chosen[j][k] for k in range(3)))
        poly = list(TRIANGLE)
        for q in sorted(normals):
            poly = clip(poly, q, 1)
            if not poly or area2(poly) <= 0:
                break
        if poly and area2(poly) > 0:
            cells[frozenset(C+chosen)] = poly
        counts['all_subset_residual_basis_candidates'] += 1
    assert sum(area2(poly) for poly in cells.values()) == 1
    return cells


def rigidity_probe(lam, counts):
    """Exhaust all specified upper-row tails, not all planar ideals."""
    a, p = tuple(sum(lam[j:]) for j in range(len(lam))), 2
    labels = tuple((i, j) for j, size in enumerate(a) for i in range(size))
    index = {v: i for i, v in enumerate(labels)}
    params = tuple((j, k, i) for j in range(len(a)) for k in range(j+1, len(a))
                   for i in range(a[k]))
    good = equality = 0
    for values in product(range(p), repeat=len(params)):
        tails = [dict() for _ in a]
        for (j, k, i), c in zip(params, values):
            if c:
                tails[j][i, k] = c

        @lru_cache(None)
        def ev(i, j):
            if j >= len(a):
                return {}
            if i < a[j]:
                return {index[i, j]: 1}
            return linear(((c, ev(i-a[j]+u, v)) for (u, v), c in tails[j].items()), p)

        X = tuple(ev(i+1, j) for i, j in labels)
        Y = tuple(ev(i, j+1) for i, j in labels)
        if any(apply(X, Y[k], p) != apply(Y, X[k], p) for k in range(len(labels))):
            continue
        good += 1
        W, _ = dense_kernel(Y, len(labels), p)
        shape = recover_partitions(X, W, ((0, 1),), p, counts)[0][1]
        phi = sum((i+1)*v for i, v in enumerate(shape))
        assert phi <= len(labels)
        if phi == len(labels):
            assert not any(values) and shape == lam
            equality += 1
    assert equality == 1
    counts['exhaustive_commuting_planar_perturbations'] += good
    return dict(primary_partition=lam, length=len(labels), coefficient_slots=params,
                all_coefficient_assignments=2**len(params), commuting_models=good,
                equality_models=equality, only_equality_parameters=[0]*len(params))


def parameter_counts(lam, counts):
    """Enumerate all B-linear automorphisms and all u for small F2 modules."""
    H, p = envelope(tuple((0,)*v+(1,) for v in lam), 2), 2
    h, aut, keys = H['h'], 0, set()
    for values in product(range(p), repeat=h*h):
        theta = tuple({i: values[j*h+i] for i in range(h) if values[j*h+i]}
                      for j in range(h))
        if dense_rank(theta, h, p) < h:
            continue
        if any(apply(theta, H['WX'][i], p) != apply(H['WX'], theta[i], p) for i in range(h)):
            continue
        aut += 1
        for coeffs in product(range(p), repeat=h):
            u = {i: c for i, c in enumerate(coeffs) if c}
            T = relations(H, theta, u, p)
            transpose = tuple({j: v.get(i, 0) for j, v in enumerate(T) if v.get(i, 0)}
                              for i in range(len(H['labels'])))
            orthogonal, rank = dense_kernel(transpose, len(T), p)
            assert rank == h+H['D']
            key = tuple(tuple(sorted(v.items())) for v in orthogonal)
            assert key not in keys
            keys.add(key)
    assert len(keys) == aut*2**h
    counts['distinct_minimal_ideal_kernels'] += len(keys)
    return dict(partition=lam, characteristic=p, automorphisms=aut,
                core_parameters=2**h, distinct_coordinate_ideals=len(keys))


def padding_control(H, theta, u, factors, t, p, counts):
    n, M, unit, _ = quotient(H, theta, u, p, counts)
    padded = tuple(op+tuple(
        {n+i+1: 1} if axis == 0 and i+1 < t else
        {n+i: 1} if axis == 1 else {}
        for i in range(t)) for axis, op in enumerate(M))
    joined_unit = linear(((1, unit), (1, {n: 1})), p)
    assert cyclic_span(padded, joined_unit, p) == n+t
    Y, Z, W = intersection(padded, p)
    primary = recover_partitions(padded[0], W, factors, p, counts)
    assert len(W) == H['h'] and len(Y) == H['phi']+t and len(Z) == H['phi']
    assert sum((len(f)-1)*sum((i+1)*v for i, v in enumerate(shape))
               for f, shape in primary) == H['phi']
    counts['all_field_disjoint_padding_controls'] += 1
    return dict(characteristic=p, invariant_factors=H['ds'], added_length=t,
                total_length=n+t, overlap_partitions=primary,
                extra_component='K[x]/x^t, y=1, z=0')


def hashes():
    paths = ['docs/minimal_branch_overlap.md', 'scripts/run_minimal_branch_overlap_checks.py']
    pending = [paths[-1]]
    while pending:
        path = pending.pop()
        for node in ast.walk(ast.parse((ROOT/path).read_text())):
            names = ([node.module] if isinstance(node, ast.ImportFrom) and node.module else
                     [a.name for a in node.names] if isinstance(node, ast.Import) else [])
            for name in names:
                candidate = 'scripts/'+name.replace('.', '/')+'.py'
                if (ROOT/candidate).is_file() and candidate not in paths:
                    paths.append(candidate); pending.append(candidate)
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def main():
    if not __debug__:
        raise RuntimeError('Assertions must be enabled.')
    before, counts, rows, padded_rows = hashes(), Counter(), [], []
    for p in (0, 2, 3, 101):
        for size in range(1, 6):
            for lam in partitions(size):
                H = envelope(tuple((0,)*v+(1,) for v in lam), p)
                for name, theta, u in parameters(H, p):
                    # Full fan equality uses a separate all-subset compiler.
                    full = size <= 3 or (lam in ((2, 2), (2, 1, 1)) and name == 'off_diagonal_shear')
                    rows.append(verify_model(str(lam)+'/'+name, H, theta, u,
                                             ((0, 1),), p, counts, full))
        ds = (mul(power((0, 1), 2, p), (-1, 1), p),
              mul((0, 1), (-1, 1), p), (-1, 1))
        H = envelope(ds, p)
        for name, theta, u in parameters(H, p):
            rows.append(verify_model('two_primary/'+name, H, theta, u,
                                     ((0, 1), (-1, 1)), p, counts))
        name, theta, u = parameters(H, p)[-1]
        for t in (1, 3):
            padded_rows.append(padding_control(H, theta, u,
                                               ((0, 1), (-1, 1)), t, p, counts))
        print(json.dumps(dict(characteristic=p, completed_minimal_models=len(rows))), flush=True)
    f = (1, 1, 1)
    H = envelope((mul(f, f, 2), f), 2)
    for name, theta, u in parameters(H, 2):
        rows.append(verify_model('irreducible_quadratic/'+name, H, theta, u,
                                 (f,), 2, counts))
    padded_rows.append(padding_control(H, theta, u, (f,), 2, 2, counts))
    rigidity = [rigidity_probe(lam, counts) for lam in
                ((1, 1), (2, 1), (2, 2), (1, 1, 1), (2, 1, 1),
                 (2, 2, 1), (2, 2, 2), (3, 2, 1), (3, 3, 1))]
    finite_counts = [parameter_counts(lam, counts) for size in range(1, 4)
                     for lam in partitions(size)]
    pascal = []
    for p in (0, 2, 3, 101):
        for m in range(1, 6):
            n, ev = pascal_ray_model(m, p)
            actual = fan(n, ev, p, counts)
            q = m
            if p:
                while q > 1 and q % p == 0:
                    q //= p
            collapse = p > 0 and q == 1
            if not p or p > 3*m-1 or collapse:
                compare_fans(actual, pascal_ray_predicted(m, collapse), counts)
            M = tuple(tuple(ev(tuple(a+int(i == j) for j, a in enumerate((k, 0, 0))))
                            for k in range(n)) for i in range(3))
            Y, Z, W = intersection(M, p)
            shape = recover_partitions(M[0], W, ((0, 1),), p, counts)[0][1]
            assert len(Y) == len(Z) == len(W) == m and shape == (m,)
            pascal.append(dict(characteristic=p, m=m, length=n, overlap_partition=shape,
                               cones=len(actual), prime_power_collapse=collapse))
    assert before == hashes(), 'Source changed during verification.'
    report = dict(command='python3 -B scripts/run_minimal_branch_overlap_checks.py',
                  status='passed', source_sha256=before, counts=dict(counts),
                  uniform_proof='docs/minimal_branch_overlap.md',
                  minimal_quotient_witnesses=rows, exhaustive_planar_tail_probes=rigidity,
                  finite_field_parameter_counts=finite_counts, pascal_sharpness_controls=pascal,
                  disjoint_padding_controls=padded_rows,
                  primary_sources=[
                      dict(url='https://doi.org/10.1090/S0002-9947-1963-0143105-7',
                           scope='Andrews lattice-polytope vertex bound; an established ingredient.'),
                      dict(url='https://travisadillon.com/documents/volume_vs_vertex_in_lattice_polytopes.pdf',
                           scope='Theorem 1 states the vertex/volume exponent (d-1)/(d+1); author-hosted modern proof.'),
                      dict(url='https://arxiv.org/abs/2304.08550',
                           scope='Commuting Jordan types concern a related, different optimization. Exact comparison of minimum and rigidity statements remains open.'),
                      dict(url='https://pi.math.cornell.edu/~frankmoore/Papers/connected.pdf',
                           scope='Sections 1--2 already use fiber products and diagonal module quotients; no new gluing operation is asserted.')],
                  limits=[
                      'All-field and all-length claims depend on the proof, not these finite computations.',
                      'Primary factors are explicitly supplied; the checker is not a general polynomial factorizer.',
                      'Upper-row probes exhaust only their recorded finite coefficient slots, not all planar ideals.',
                      'The residual compiler enumerates abstract matroid bases and all fundamental exchanges; it does not treat the contraction as an algebra quotient.',
                      'Quadratic full-fan bound applies to minimum length n=2*Phi(W), not arbitrary yz=0 quotients.',
                      'The Pascal quadratic lower bound uses characteristic zero or the stated sufficiently large characteristics; prime-power collapse is recorded.',
                      'Exact priority, publication significance, unrestricted G3 and the active goal remain unresolved.'])
    (ROOT/'results/research/minimal_branch_overlap_checks.json').write_text(
        json.dumps(jsonable(report), indent=2)+'\n')
    print(json.dumps(dict(status='passed', counts=dict(counts))), flush=True)


if __name__ == '__main__':
    main()
