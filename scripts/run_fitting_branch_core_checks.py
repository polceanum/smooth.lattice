#!/usr/bin/env python3
"""Exact canonical quotient maps, nonminimal gluings and excess fan budgets."""
import ast
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import prod
from pathlib import Path
import hashlib
import json

from run_minimal_branch_overlap_checks import (
    envelope, intersection, module_data, mul, parameters, partitions, power,
    quotient, recover_partitions, rem, residual_fan,
)
from run_univariate_branch_gluing_checks import (
    Frame, apply, axis_padding, clean, cyclic_span, dense_kernel, dense_rank, jsonable,
    linear, pascal_ray_model, polynomial,
)
from run_forced_socle_fan_checks import model as forced_model
from run_interleaving_research import minimal_excluded_border
from run_projection_tag_checks import compare_fans, fan
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import dot, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]


def planar(a, values, p):
    r = len(a)
    labels = tuple((i, j) for j, bound in enumerate(a) for i in range(bound))
    index = {v: i for i, v in enumerate(labels)}
    slots = tuple((j, k, i) for j in range(r) for k in range(j+1, r) for i in range(a[k]))
    assert len(values) == len(slots)
    tails = [dict() for _ in a]
    for (j, k, i), c in zip(slots, values):
        if c:
            tails[j][i, k] = c

    @lru_cache(None)
    def ev(i, j):
        if j >= r:
            return {}
        if i < a[j]:
            return {index[i, j]: 1}
        return linear(((c, ev(i-a[j]+u, v)) for (u, v), c in tails[j].items()), p)

    X, Y = tuple(ev(i+1, j) for i, j in labels), tuple(ev(i, j+1) for i, j in labels)
    return dict(a=a, values=values, slots=slots, tails=tails, labels=labels,
                X=X, Y=Y, unit={index[0, 0]: 1}, ev=ev)


def planar_controls(counts):
    summaries, pool, seen = [], [], set()
    for n in range(1, 12):
        for a in partitions(n):
            if len(a) > 4:
                continue
            slots = sum(a[k] for j in range(len(a)) for k in range(j+1, len(a)))
            if slots > 9:
                continue
            commuting = maximum_excess = 0
            for values in product(range(2), repeat=slots):
                C = planar(a, values, 2)
                counts['planar_tail_assignments'] += 1
                X, Y = C['X'], C['Y']
                if any(apply(X, Y[i], 2) != apply(Y, X[i], 2) for i in range(n)):
                    continue
                commuting += 1
                W, _ = dense_kernel(Y, n, 2)
                lam = recover_partitions(X, W, ((0, 1),), 2, counts)[0][1]
                heights = tuple(sum(lam[j:]) for j in range(len(lam)))
                for j, bound in enumerate(a):
                    if j < len(heights):
                        assert bound >= heights[j]
                    for (i, k), c in C['tails'][j].items():
                        assert k >= len(heights) or i >= heights[k]
                        counts['actual_planar_relation_coefficients'] += 1
                phi = sum((i+1)*v for i, v in enumerate(lam))
                excess = n-phi
                assert excess >= 0
                maximum_excess = max(maximum_excess, excess)
                key = (lam, excess, bool(any(values)))
                if (2 <= len(lam) <= 3 and 0 < excess <= 3 and sum(lam) <= 5
                        and phi <= 10 and key not in seen):
                    seen.add(key)
                    C.update(lam=lam, phi=phi, W=W)
                    pool.append(C)
            summaries.append(dict(length=n, row_widths=a, coefficient_slots=slots,
                                  assignments=2**slots, commuting_models=commuting,
                                  maximum_excess=maximum_excess))
            counts['commuting_planar_coefficient_models'] += commuting
    assert counts['planar_tail_assignments'] == 9735
    assert counts['commuting_planar_coefficient_models'] == 1691
    return summaries, pool


def jordan_frame(X, W, lam, p):
    """Find module chains from actual kernels, independently of row labels."""
    F = Frame(W, p)
    action = tuple(F.coords(apply(X, v, p)) for v in W)
    chosen = []
    for size in lam:
        columns = tuple({i: 1} for i in range(len(W)))
        for _ in range(size):
            columns = tuple(apply(action, v, p) for v in columns)
        killed, _ = dense_kernel(columns, len(W), p)
        for v in killed:
            chain = [v]
            for _ in range(size-1):
                chain.append(apply(action, chain[-1], p))
            if dense_rank(tuple(chosen)+tuple(chain), len(W), p) == len(chosen)+size:
                chosen.extend(chain)
                break
        else:
            raise AssertionError('No required Jordan chain.')
    assert len(chosen) == len(W)
    return tuple(apply(W, v, p) for v in chosen)


def quotient_by_relations(M, unit, T, p):
    n = len(M[0])
    # Close only under x; the supplied relations must be killed by y,z.
    old = -1
    while True:
        F = Frame(T, p)
        if len(F.basis) == old:
            break
        old = len(F.basis)
        T = tuple(F.basis)+tuple(apply(M[0], v, p) for v in F.basis)
    T = F.basis
    assert all(not apply(op, v, p) for op in M[1:] for v in T)
    full = Frame(T+tuple({i: 1} for i in range(n)), p)
    k = len(T)

    def project(v):
        return clean({i-k: c for i, c in full.coords(v).items() if i >= k}, p)

    ops = tuple(tuple(project(apply(op, v, p)) for v in full.basis[k:]) for op in M)
    return n-k, ops, project(unit)


def glue(C, p, coupled):
    """Two identical genuine planar branches, glued along their full kernels.

    The optional relation x^h+v_Y+v_Z has v outside ker(y) in general.
    It tests coupling through the axis, beyond a diagonal branch quotient.
    """
    X, Y, unit = C['X'], C['Y'], C['unit']
    m, h = len(X), sum(C['lam'])
    a, v = 0, unit
    while v:
        a += 1; v = apply(X, v, p)
        assert a <= m
    assert all(not polynomial(X, {i: 1}, (0,)*a+(1,), p) for i in range(m))
    N = a+2*m

    def lift(v, offset):
        return {offset+i: c for i, c in v.items()}

    core = tuple({i+1: 1} if i+1 < a else {} for i in range(a))
    HX = core+tuple(lift(v, a) for v in X)+tuple(lift(v, a+m) for v in X)
    xpowers, v = [], unit
    for _ in range(a):
        xpowers.append(v); v = apply(X, v, p)
    HY = tuple(lift(v, a) for v in xpowers)+tuple(lift(v, a) for v in Y)+({},)*m
    HZ = tuple(lift(v, a+m) for v in xpowers)+({},)*m+tuple(lift(v, a+m) for v in Y)
    W, _ = dense_kernel(Y, m, p)
    W = jordan_frame(X, W, C['lam'], p)
    T = tuple(linear(((1, lift(v, a)), (1, lift(v, a+m))), p) for v in W)
    if coupled and h < a:
        target = linear(((-1, xpowers[h]),), p)
        FY = Frame(Y, p)
        co = FY.coords(target)
        v = {FY.selected[j]: c for j, c in co.items()}
        assert apply(Y, v, p) == target
        T += (linear(((1, {h: 1}), (1, lift(v, a)), (1, lift(v, a+m))), p),)
    out = quotient_by_relations((HX, HY, HZ), {0: 1}, T, p)
    assert out[0] == (2*m if coupled else a+2*m-h)
    assert N >= out[0]
    return out


def base_change(data, f, p):
    """Free monic substitution t=f(x), including irreducible/inseparable f."""
    n, M, unit = data
    d = len(f)-1
    assert f[-1] == 1
    X, Y, Z = [], [], []
    for j in range(n):
        for k in range(d):
            if k+1 < d:
                X.append({j*d+k+1: 1})
            else:
                X.append(linear(((1, {i*d: c for i, c in M[0][j].items()}),
                                 (-1, {j*d+i: c for i, c in enumerate(f[:-1]) if c})), p))
            Y.append({i*d+k: c for i, c in M[1][j].items()})
            Z.append({i*d+k: c for i, c in M[2][j].items()})
    return n*d, (tuple(X), tuple(Y), tuple(Z)), {i*d: c for i, c in unit.items()}


def pad(data, t, p):
    n, M, unit = data
    out = tuple(op+tuple({n+i+1: 1} if axis == 0 and i+1 < t else
                         {n+i: 1} if axis == 1 else {} for i in range(t))
                for axis, op in enumerate(M))
    return n+t, out, linear(((1, unit), (1, {n: 1})), p)


def matrix_annihilator(M, p):
    n = len(M)
    if not n:
        return (1,)
    powers, op = [], tuple({i: 1} for i in range(n))
    for k in range(n+1):
        vector = {j*n+i: c for j, v in enumerate(op) for i, c in v.items()}
        F = Frame(powers+[vector], p)
        if len(F.basis) == len(powers):
            co = Frame(powers, p).coords(vector)
            result = tuple(-co.get(i, 0) for i in range(k))+(1,)
            assert all(not polynomial(M, {i: 1}, result, p) for i in range(n))
            return result
        powers.append(vector)
        op = tuple(apply(M, v, p) for v in op)
    raise AssertionError('Matrix annihilator degree exceeds dimension.')


def inspect(name, data, factors, p, counts, full=False):
    n, M, unit = data
    assert cyclic_span(M, unit, p) == n
    for U, V in combinations(M, 2):
        for j in range(n):
            assert apply(U, V[j], p) == apply(V, U[j], p)
            counts['actual_commuting_columns'] += 1
    assert all(not apply(M[1], v, p) for v in M[2])

    @lru_cache(None)
    def ev(a):
        v = unit
        for op, exponent in zip(M, a):
            for _ in range(exponent):
                v = apply(op, v, p)
        return v

    Y, Z, W = intersection(M, p)
    assert W
    primary = recover_partitions(M[0], W, factors, p, counts)
    r = max(len(shape) for _, shape in primary)
    ds = [(1,)]*r
    for f, shape in primary:
        for j, size in enumerate(shape):
            ds[j] = mul(ds[j], power(f, size, p), p)
    gs, h, phi = module_data(tuple(ds), p)
    eY, eZ, s, D = len(Y)-phi, len(Z)-phi, n-2*phi, len(ds[0])-1
    d = n-len(Y)-len(Z)+h
    assert eY >= 0 and eZ >= 0 and s == eY+eZ+d-h and d >= h
    C = {(i, 0, 0) for i in range(h)}
    C.update((i, j, 0) for j in range(1, r) for i in range(len(gs[j])-1))
    C.update((i, 0, j) for j in range(1, r) for i in range(len(gs[j])-1))
    C = tuple(sorted(C)); index = {a: i for i, a in enumerate(C)}
    assert len(C) == 2*phi-h and dense_rank(tuple(ev(a) for a in C), n, p) == len(C)

    @lru_cache(None)
    def canonical(a):
        x, y, z = a
        if y and z or max(y, z) >= r:
            return {}
        coeffs = rem((0,)*x+(1,), gs[max(y, z)], p)
        return clean({index[i, y, z]: c for i, c in enumerate(coeffs) if c}, p)

    ground = tuple(hyperbolic_cross(n, 3))
    B, _ = field_basis(sorted(ground, key=lambda a: (dot(a, (1, n+1, (n+1)**2)), a)), ev, p)
    F = Frame(tuple(ev(a) for a in B), p)
    assert len(F.basis) == n
    target = tuple(canonical(a) for a in B)
    projection = tuple(apply(target, F.coords({i: 1}), p) for i in range(n))
    canon_ops = tuple(tuple(canonical(tuple(a+int(i == j) for j, a in enumerate(v)))
                            for v in C) for i in range(3))
    assert apply(projection, unit, p) == canonical((0, 0, 0))
    for op, cop in zip(M, canon_ops):
        for i in range(n):
            assert apply(projection, op[i], p) == apply(cop, projection[i], p)
            counts['independent_canonical_quotient_map_columns'] += 1
    assert dense_rank(projection, len(C), p) == len(C)
    assert all(not apply(projection, v, p) for v in W)
    V, _ = dense_kernel(projection, len(C), p)
    assert len(V) == h+s
    FV = Frame(tuple(W)+V, p)
    assert len(FV.basis) == h+s

    def residual(v):
        return clean({i-h: c for i, c in FV.coords(v).items() if i >= h}, p)

    quotient_X = tuple(residual(apply(M[0], v, p)) for v in FV.basis[h:])
    ann = matrix_annihilator(quotient_X, p)
    assert len(ann)-1 <= s
    ann_A = mul(mul(gs[0], ds[0], p), ann, p)
    assert len(ann_A)-1 <= h+D+s
    assert not polynomial(M[0], unit, ann_A, p)
    counts['actual_defect_quotient_annihilators'] += 1

    DW = {(i, j) for j in range(r) for i in range(len(gs[j])-1)}
    border = tuple(minimal_excluded_border(tuple(DW)))

    def shell(e):
        if not e:
            return set()
        return {tuple(a[i]+b[i] for i in range(2))
                for a in border for b in hyperbolic_cross(e, 2)}

    allowed = {(i, 0, 0) for i in range(h+D+s)}
    allowed.update((i, j+1, 0) for i, j in DW | shell(eY))
    allowed.update((i, 0, j+1) for i, j in DW | shell(eZ))
    allowed = {a for a in allowed if prod(i+1 for i in a) <= n}
    E = tuple(sorted(allowed-set(C)))
    _, cpivots = field_basis(C, ev, p)
    residual_columns = tuple(field_reduce(ev(a), cpivots, p)[0] for a in E)
    assert dense_rank(residual_columns, n, p) == h+s
    assert dense_rank(tuple(ev(a) for a in allowed), n, p) == n
    candidate_bound = (2*h+D+s+(r+1)*
                       (len(tuple(hyperbolic_cross(eY,2))) if eY else 0))
    candidate_bound += (r+1)*(len(tuple(hyperbolic_cross(eZ,2))) if eZ else 0)
    assert len(E) <= candidate_bound
    counts['actual_whole_fan_fixed_shell_matroids'] += 1
    xbound = (h+s)*(h+D+s-1)
    ybound, zbound = phi+r*eY+eY*(eY+1)//2, phi+r*eZ+eZ*(eZ+1)//2
    row = dict(name=name, characteristic=p, length=n, invariant_factors=ds,
               primary_partitions=primary, h=h, phi=phi, r=r, e_Y=eY, e_Z=eZ,
               axis_quotient_length=d, excess=s, branch_excess=eY+eZ,
               common_standard_set=C,
               residual_ground=E, residual_ground_bound=candidate_bound,
               defect_x_annihilator=ann,
               full_x_annihilator=ann_A, residual_moment_bounds=(xbound,ybound,zbound))
    if full:
        cells = fan(n, ev, p, counts)
        for S, poly in cells.items():
            assert set(C) <= S <= allowed
            R = S-set(C)
            moments = tuple(sum(a[i] for a in R) for i in range(3))
            assert all(u <= v for u, v in zip(moments, (xbound,ybound,zbound)))
            assert len(minimal_excluded_border(tuple(S))) <= 2*r+eY+eZ+4
            w = tuple(sum(v[i] for v in poly) for i in range(3))
            for axis, m, e in ((1,len(Y),eY),(2,len(Z),eZ)):
                def branch(a):
                    return ev((a[0], a[1]+1 if axis == 1 else 0,
                               a[1]+1 if axis == 2 else 0))
                Bbranch, _ = field_basis(sorted(hyperbolic_cross(m,2),
                                                key=lambda a:(a[0]*w[0]+a[1]*w[axis],a)), branch,p)
                Bbranch = set(Bbranch)
                assert DW <= Bbranch and len(Bbranch-DW) == e
                actual_branch = {(a[0],a[axis]-1) for a in S if a[axis]}
                assert actual_branch <= Bbranch
                assert Bbranch-DW <= shell(e)
                counts['independent_full_branch_greedy_containments'] += 1
        if h+s <= 5 and len(E) <= 16:
            predicted = residual_fan(C,E,h+s,ev,p,counts)
            compare_fans(cells,predicted,counts)
            counts['independent_entire_fixed_shell_compilations'] += 1
        row['cones'] = len(cells)
        row['complete_standard_sets'] = [sorted(S) for S in sorted(cells,key=lambda b:sorted(b))]
    counts['actual_fitting_core_reconstructions'] += 1
    return row


def hashes():
    paths = ['docs/fitting_branch_core.md','scripts/run_fitting_branch_core_checks.py',
             'docs/minimal_branch_overlap.md']
    pending = [paths[1]]
    while pending:
        path = pending.pop()
        for node in ast.walk(ast.parse((ROOT/path).read_text())):
            names = ([node.module] if isinstance(node,ast.ImportFrom) and node.module else
                     [a.name for a in node.names] if isinstance(node,ast.Import) else [])
            for name in names:
                candidate = 'scripts/'+name.replace('.','/')+'.py'
                if (ROOT/candidate).is_file() and candidate not in paths:
                    paths.append(candidate);pending.append(candidate)
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def main():
    if not __debug__:
        raise RuntimeError('Assertions must be enabled.')
    before, counts = hashes(), Counter()
    probes,pool = planar_controls(counts)
    print(json.dumps(dict(planar_models=counts['commuting_planar_coefficient_models'],
                          available_nonminimal_branches=len(pool))),flush=True)
    rows, witnesses = [], []
    chosen = pool[:16]
    for j,C in enumerate(chosen):
        witnesses.append(dict(row_widths=C['a'],coefficient_slots=C['slots'],
                              coefficients=C['values'],kernel_partition=C['lam']))
        for coupled in (False,True):
            data=glue(C,2,coupled)
            rows.append(inspect('planar_gluing_'+str(j)+'_coupled='+str(coupled),data,
                                ((0,1),),2,counts,full=j<8))
    for p in (0,2,3,101):
        for lam in ((2,),(1,1),(2,1),(2,1,1)):
            H=envelope(tuple((0,)*v+(1,) for v in lam),p)
            name,theta,u=parameters(H,p)[-1]
            n,M,unit,_=quotient(H,theta,u,p,counts)
            data=n,M,unit
            rows.append(inspect('minimum_'+str(lam),data,((0,1),),p,counts,full=sum(lam)<=3))
            if lam==(1,1):
                for t in (1,3):
                    rows.append(inspect('unit_y_padding_'+str(t),pad(data,t,p),
                                        ((0,1),),p,counts,full=t==1))
                row=inspect('nonlocal_axis_padding_5',axis_padding(data,5,p),
                            ((0,1),),p,counts,full=True)
                assert row['branch_excess']==0 and row['excess']==5
                rows.append(row)
        a=(3,2);C=planar(a,(0,)*2,p)
        W,_=dense_kernel(C['Y'],sum(a),p)
        C['lam']=recover_partitions(C['X'],W,((0,1),),p,counts)[0][1]
        rows.append(inspect('monomial_planar_nonminimal_core_coupling',glue(C,p,True),
                            ((0,1),),p,counts,full=True))
        for N in (3,4):
            reference,ev=forced_model(N,'A',p)
            M=tuple(tuple(ev(tuple(a+int(i==j) for j,a in enumerate(v)))
                          for v in reference) for i in range(3))
            rows.append(inspect('forced_socle_N='+str(N),(len(reference),M,ev((0,0,0))),
                                ((0,1),),p,counts,full=True))
        print(json.dumps(dict(characteristic=p,completed_actual_models=len(rows))),flush=True)
    for f,p,factors in (((1,1,1),2,((1,1,1),)),((0,0,1),2,((0,1),)),
                         ((0,-1,1),101,((0,1),(-1,1)))):
        C=planar((3,2),(0,0),p)
        C['lam']=(2,1)
        data=base_change(glue(C,p,True),f,p)
        rows.append(inspect('free_monic_base_change_'+str(f),data,factors,p,counts,full=True))
    for m in (2,3,4):
        n,ev=pascal_ray_model(m,0)
        M=tuple(tuple(ev((i+1,0,0)) if axis==0 else
                      ev((i,int(axis==1),int(axis==2))) for i in range(n)) for axis in range(3))
        row=inspect('Pascal_sharp_m='+str(m),(n,M,ev((0,0,0))),((0,1),),0,counts,full=True)
        assert row['excess']==0 and row['cones']==(m+1)*(m+2)//2
        rows.append(row)
    assert before==hashes(),'A source changed during verification.'
    report=dict(command='python3 -B scripts/run_fitting_branch_core_checks.py',status='passed',
                source_sha256=before,counts=dict(counts),uniform_proof='docs/fitting_branch_core.md',
                exhaustive_planar_coefficient_controls=probes,nonminimal_branch_witnesses=witnesses,
                actual_quotient_witnesses=rows,
                primary_sources=[
                    dict(url='https://arxiv.org/html/0708.3576v1',scope='Conca--Valla Section 3, Theorem 3.3, Lemma 3.6 and the proof of Step 2: polynomial Hilbert--Burch representations of arbitrary finite planar ideals.'),
                    dict(url='https://arxiv.org/html/2407.07993v4',scope='Oszer Theorem 2.1 gives the classical resolution theorem. Introduction and Section 6 warn against assuming arbitrary matrix perturbations preserve colength.'),
                    dict(url='https://doi.org/10.1090/S0002-9947-1963-0143105-7',scope='The established lattice vertex-volume estimate; no new geometric theorem.'),
                    dict(url='https://travisadillon.com/documents/volume_vs_vertex_in_lattice_polytopes.pdf',scope='Theorem 1 states and proves the vertex-volume exponent used in dimension three.')],
                limits=['Coefficient containment is a classical Hilbert--Burch/Smith-form corollary, strengthening and reclassifying the previous planar minimum lemma.',
                        'The finite upper-row controls exhaust only their stated coefficient slots, not all planar ideals.',
                        'Every claimed quotient map is checked against all multiplication columns; sampled weights are not used to infer a global quotient.',
                        'Complete fan comparisons use exact rational cones. The independent all-subset compiler is used only where the recorded residual size permits it.',
                        'The theorem gives quadratic order for total planar branch excess O(sqrt(n)), irrespective of the extra axis length; it does not bound all yz=0 ideals quadratically.',
                        'No unrestricted G3 improvement, supergeneric polynomial separation, completed priority audit, significance judgment or completed goal follows.'])
    (ROOT/'results/research/fitting_branch_core_checks.json').write_text(json.dumps(jsonable(report),indent=2)+'\n')
    print(json.dumps(dict(status='passed',actual_models=len(rows),counts=dict(counts))),flush=True)


if __name__=='__main__':
    main()
