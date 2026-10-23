#!/usr/bin/env python3
"""Actual row components, monomial chains and common integer weight lists."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, gcd, lcm, prod
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_hook_bridge_checks import model, actual, normals, polygon, strict_weight, COUNTS as REPLAY
from run_hook_rank_departure_checks import rank, clean, plane, weight_controls, complete_fan, comparison, COUNTS as FAN
from run_hook_high_price_partner_checks import bareiss, minor_product
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import dot, clip, area2

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()
TRIANGLE = [tuple(Q(i == j) for i in range(3)) for j in range(3)]


def components(data, cols, t, arm):
    rows = {i for i, v in enumerate(data[0]) if v[arm] >= t}
    parent = {i: i for i in rows}

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i

    for v in data[4][arm]:
        support = tuple(cols[v]); assert set(support) <= rows
        if support:
            first = find(support[0])
            for i in support[1:]:
                parent[find(i)] = first
    marked = {find(i) for v in data[4][2] for i in cols[v] if i in rows}
    far = {i for i in rows if find(i) not in marked}
    COUNTS['independent_actual_support_graph_partitions'] += 1
    return far


def orbit_key(q, v, nu, delta):
    return delta*q-nu*v, v % delta


def orbit_partition(data, cols, t, a, arm, nu, delta):
    A = a-t; h = (t+delta-1)//delta
    marked = {orbit_key(q, v, nu, delta) for q, v in product(range(t), repeat=2)}
    far = {i for i, v in enumerate(data[0]) if v[arm] >= t and
           orbit_key(v[arm]-t, v[1-arm], nu, delta) not in marked}
    for v in data[4][arm]:
        keys = {orbit_key(data[0][i][arm]-t, data[0][i][1-arm], nu, delta) for i in cols[v]}
        assert len(keys) <= 1
        COUNTS['independent_pure_column_orbit_support_checks'] += 1
    for v in data[4][2]:
        assert all(orbit_key(data[0][i][arm]-t, data[0][i][1-arm], nu, delta) in marked
                   for i in cols[v] if data[0][i][arm] >= t)
    assert far <= components(data, cols, t, arm)
    S = {q for k in range(h) for q in range(k*nu, min(A, k*nu+t))}
    for q, v in product(range(A), range(t)):
        if orbit_key(q, v, nu, delta) in marked:
            assert q in S
    beta = (min(Q(h), 1+Q(nu, delta), Q((A+nu-1)//nu), Q(A, t)) if nu else Q(1))
    spread = {q+j for q in S for j in range(-2*t, 2*t+1)}
    assert len(spread) <= 16*t*beta
    COUNTS['independent_chain_band_and_refined_beta_bounds'] += 1
    return far, marked, beta


def core_matrix(data, cols, far, arm):
    omitted = tuple(v for v, c in cols.items() if c and set(c) <= far)
    assert set(omitted) <= set(data[4][arm])
    core = {v: c for v, c in cols.items() if c and v not in omitted}
    assert all(not (set(c) & far) for c in core.values())
    return omitted, core


def check_patch(data, cols, far, arm, w, p):
    omitted, core = core_matrix(data, cols, far, arm)
    B = field_basis(sorted(cols, key=lambda v: (dot(w, v), v)), cols.__getitem__, p)[0]
    fixed = field_basis(sorted(omitted, key=lambda v: (dot(w, v), v)), cols.__getitem__, p)[0]
    chosen = field_basis(sorted(core, key=lambda v: (dot(w, v), v)), core.__getitem__, p)[0]
    assert set(B) == set(fixed)|set(chosen)
    order = sorted(data[4][arm], key=lambda v: dot(w, v))
    assert len({dot(w, v) for v in order}) == len(order)
    O = {tuple(b-a for a, b in zip(v, u)) for v, u in zip(order, order[1:])}
    N = normals(B, cols, p); T = normals(chosen, core, p)
    old = polygon(N|O); new = polygon(T|O)
    assert all(dot(v, q) >= 0 for v in old for q in T)
    assert all(dot(v, q) >= 0 for v in new for q in N)
    assert area2(old) == area2(new) > 0
    COUNTS['independent_entire_original_vs_component_core_patches'] += 1
    return {'omitted_rank': len(fixed), 'core_rank': len(chosen),
            'core_columns': len(core), 'patch_area2': str(area2(old))}


def inventory(data, core_labels, arm):
    labels = tuple(core_labels)
    out = {plane(tuple(x-y for x, y in zip(v, u))) for j, v in enumerate(labels) for u in labels[:j]}
    group = data[4][arm]
    out |= {plane(tuple(x-y for x, y in zip(v, u))) for j, v in enumerate(group) for u in group[:j]}
    return out


def potential_core(data, t, arm, nu, delta, marked):
    out = set(data[4][1-arm])|set(data[4][2])
    for v in data[4][arm]:
        long, short = v[arm+1], v[2-arm]
        if v[0]:
            q, j = long-t, short
        elif short >= t:
            q, j = long-t, short-t
        else:
            continue  # chosen-arm outer pure columns are uniformly zero
        if orbit_key(q, j, nu, delta) in marked:
            out.add(v)
    return out


def ray(v):
    den = lcm(*(x.denominator for x in v))
    values = tuple(int(x*den) for x in v); d = gcd(*values)
    return tuple(x//d for x in values)


def det(a, b, c):
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])


def common_list(planes, a):
    cells = [TRIANGLE]
    for q in sorted(planes):
        nxt = []
        for poly in cells:
            for sign in (-1, 1):
                part = clip(poly, q, sign)
                if part and area2(part):
                    nxt.append(part)
        cells = nxt
    assert sum(map(area2, cells)) == 1
    weights = []
    for poly in cells:
        rays = tuple(map(ray, poly))
        triple = next(v for v in combinations(rays, 3) if det(*v))
        w = tuple(sum(v[j] for v in triple) for j in range(3))
        assert min(w) > 0 and max(w) <= 100*a*a
        assert all(dot(w, q) != 0 for q in planes)
        weights.append(w)
    assert len(set(weights)) == len(weights)
    COUNTS['independent_complete_common_arrangements'] += 1
    COUNTS['independent_integer_sampler_chambers'] += len(weights)
    return weights


def build(t, a, arm, nu, delta, gamma, p, rng, symmetric=False):
    data = model(t, a); W = {(k*nu, k*delta): c for k, c in enumerate(gamma)
                            if k*nu < a-t and k*delta < t}
    co = []
    for g, y, z in data[2]:
        long, short = (y, z) if arm == 0 else (z, y)
        if long >= t:
            c = W.get((long-t, short), 0) if g == 1 else 0
        elif symmetric:
            c = W.get((short-t, long), 0) if g == 1 else 0
        else:
            c = rng.randrange(p or 7)
        co.append(c)
    return data, co, actual(data, co, p, True)


def large_rank_controls(data, co, cols, t, a, p):
    A = a-t; assert A >= 2*t-1 and (not p or p > 2*t)
    expected = (A-1)*(t-1); chis = []
    for arm in (0, 1):
        far = {i for i, v in enumerate(data[0]) if v[arm] >= 2*t}
        projected = {v: {i: c for i, c in cols[v].items() if i in far} for v in data[4][2]}
        chi = rank(tuple(projected), projected, p); chis.append(chi)
        for d in range(1, t):
            rows = [(v+d, v) for v in range(t-d, t)]
            inputs = [(j+d, j) for j in range(t-d)]
            blocks = {}
            for i, j in inputs:
                label = (0, t+i, t+j) if arm == 0 else (0, t+j, t+i)
                block = {(data[0][r][arm]-t, data[0][r][1-arm]): c
                         for r, c in projected[label].items()}
                assert set(block) <= set(rows)
                blocks[i, j] = block
            assert rank(tuple(blocks), blocks, p) == min(d, t-d)
            m = min(d, t-d); r = t-d
            matrix = [[comb(t-1, r+i-j) for j in range(m)] for i in range(m)]
            value = bareiss(matrix)
            assert value == minor_product(t-1, m, r) > 0 and (not p or value % p)
            COUNTS['independent_binomial_interface_block_minors'] += 1
        assert chi == t*t//4
        choices = list(product(range(p), repeat=t)) if p and p**t <= 32 else [tuple([1]+[0]*(t-1)), tuple([0]*t), tuple([1, 1]+[0]*(t-2))]
        values = []
        for C in choices:
            base = actual(data, comparison(data, co, arm, C), p)
            diff = {}
            for v in cols:
                c = Counter(cols[v]); c.subtract(base[v]); diff[v] = clean(c, p)
            value = rank(tuple(diff), diff, p); assert value >= expected; values.append(value)
            COUNTS['independent_large_departure_comparisons'] += 1
        assert min(values) == expected
    assert rank(data[4][2], cols, p) == t*t
    return {'bridge_rank': t*t, 'tail_interface_ranks': chis, 'least_departure_proved': expected}


def mixed_direction_controls(t, a, lam, p, rng):
    """A second step joins marked components; an integer minor retains rank."""
    A = a-t; assert A >= 2*t+2 and lam and (not p or p > 2*t)
    W = {(j, j): comb(t-1, j) for j in range(t)}
    W[2, 1] = lam
    data = model(t, a); co = []
    for g, y, z in data[2]:
        q, j = (y-t, z) if y >= t else (z-t, y)
        co.append(W.get((q, j), 0) if g == 1 else 0)
    cols = actual(data, co, p, True); index = {v: i for i, v in enumerate(data[0])}
    assert all(j <= i <= 2*j for (i, j), c in W.items() if c)
    ranks = []; spreads = []; determinants = []; patches = []
    expected_minor = prod(minor_product(t-1, min(d, t-d), t-d) for d in range(1, t))
    expected_departure = (A-1)*(t-1)
    for arm in (0, 1):
        def row(q, j):
            return index[(t+q, j) if arm == 0 else (j, t+q)]

        def label(q, j):
            return (0, t+q, t+j) if arm == 0 else (0, t+j, t+q)

        far = components(data, cols, t, arm)
        for q in range(t, A-2):
            first, second = cols[label(q, 0)], cols[label(q+1, 0)]
            assert {row(q, 0), row(q+1, 1), row(q+2, 1)} <= set(first)
            assert row(q+2, 1) in second and row(q+1, 0) in second
            COUNTS['independent_mixed_direction_shared_support_steps'] += 1
        assert row(t+1, 1) in cols[label(t-1, 0)] and label(t-1, 0) in data[4][2]
        assert not ({row(q, 0) for q in range(t, A-1)} & far)
        _, core = core_matrix(data, cols, far, arm)
        sigma = len({v[arm+1]-u[arm+1] for v in core if v in data[4][arm]
                     for u in data[4][1-arm]})
        assert sigma >= A-t-1; spreads.append(sigma)
        selected_rows = []; selected_columns = []; row_blocks = []; column_blocks = []
        for d in range(1, t):
            m = min(d, t-d)
            selected_rows += [row(v+d, v) for v in range(t-d, t-d+m)]
            selected_columns += [label(j+d, j) for j in range(m)]
            row_blocks += [d]*m; column_blocks += [d]*m
        matrix = [[int(cols[v].get(i, 0)) for v in selected_columns] for i in selected_rows]
        assert all(not matrix[i][j] for i, d in enumerate(row_blocks)
                   for j, e in enumerate(column_blocks) if d < e)
        determinant = bareiss(matrix)
        assert (determinant % p == expected_minor % p != 0 if p else determinant == expected_minor > 0)
        COUNTS['independent_mixed_direction_combined_integer_minors'] += 1
        determinants.append(str(determinant))
        tail = {v: {i: c for i, c in cols[v].items() if data[0][i][arm] >= 2*t}
                for v in data[4][2]}
        chi = rank(tuple(tail), tail, p)
        assert t*t//4 <= chi <= t*t; ranks.append(chi)
        values = []
        for C in ([1]+[0]*(t-1), [0]*t, [1, 1]+[0]*(t-2)):
            base = actual(data, comparison(data, co, arm, C), p); diff = {}
            for v in cols:
                c = Counter(cols[v]); c.subtract(base[v]); diff[v] = clean(c, p)
            value = rank(tuple(diff), diff, p)
            assert value >= expected_departure; values.append(value)
            COUNTS['independent_mixed_direction_departure_comparisons'] += 1
        assert values[0] == expected_departure
        patches += [check_patch(data, cols, far, arm, w, p) for w in weight_controls(cols, rng)]
    assert rank(data[4][2], cols, p) == t*t
    return {'t': t, 'a': a, 'parameter': lam, 'characteristic': p,
            'bridge_rank': t*t, 'tail_interface_ranks': ranks, 'long_difference_spreads': spreads,
            'combined_integer_determinants': determinants, 'unchanged_minor_product': str(expected_minor),
            'least_departure_proved': expected_departure, 'patches': patches}


def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    old = json.loads((ROOT/'results/research/hook_tail_interface_checks.json').read_text())
    assert old['status'] == 'passed'
    for f, sha in old['source_sha256'].items():
        assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest() == sha
    rng = Random(2026102324); cases = []
    specs = [(3, 12, 0, 3, 2, [1, 1], 2), (4, 16, 0, 4, 3, [0, 1], 101),
             (4, 16, 1, 1, 1, [1, 2, 1, 1], 3), (3, 10, 0, 0, 1, [0, 1, 2], 0),
             (5, 17, 1, 7, 1, [1, 1, 1, 1, 1], 2),
             (7, 22, 0, 6, 4, [1, 1], 3), (3, 9, 1, 20, 1, [1, 1, 1], 2),
             (3, 6, 0, 2, 1, [0, 1, 1], 3)]
    specs += [(t, a, 0, 1, 1, [comb(t-1, k) for k in range(t)], p)
              for t, a, p in ((2, 7, 5), (3, 11, 0), (4, 16, 101), (6, 22, 101), (8, 28, 101))]
    for t, a, arm, nu, delta, gamma, p in specs:
        symmetric = nu == delta == 1 and gamma == [comb(t-1, k) for k in range(t)] and (not p or p > 2*t)
        data, co, cols = build(t, a, arm, nu, delta, gamma, p, rng, symmetric)
        far, marked, beta = orbit_partition(data, cols, t, a, arm, nu, delta)
        _, core = core_matrix(data, cols, far, arm)
        possible = potential_core(data, t, arm, nu, delta, marked)
        assert set(core) <= possible
        sigma = len({v[arm+1]-u[arm+1] for v in core if v in data[4][arm] for u in data[4][1-arm]})
        assert sigma <= 16*t*beta
        planes = inventory(data, possible, arm)
        assert len(planes) <= 4096*a*t*beta
        patches = [check_patch(data, cols, far, arm, w, p) for w in weight_controls(cols, rng)]
        # The canonical actual support graph gives a second, finer decomposition.
        canonical = components(data, cols, t, arm)
        patches.append(check_patch(data, cols, canonical, arm, weight_controls(cols, rng)[-1], p))
        ranks = large_rank_controls(data, co, cols, t, a, p) if symmetric else None
        cases.append({'t': t, 'a': a, 'nu': nu, 'delta': delta, 'coefficients': gamma,
                      'characteristic': p, 'arm': arm, 'beta': str(beta), 'cross_arm_long_differences': sigma,
                      'uniform_planes': len(planes), 'patches': patches, 'large_parameter_controls': ranks})
        print('Orbit case passed', t, a, nu, delta, p, 'beta', beta, flush=True)
    mixed = []
    for t, a, lam, p in ((2, 10, 1, 5), (3, 14, -1, 0), (4, 20, 2, 101), (6, 26, 1, 101)):
        mixed.append(mixed_direction_controls(t, a, lam, p, rng))
        print('Mixed-direction case passed', t, a, lam, p, flush=True)
    general = []
    for t, a, p, arm in ((2, 5, 101, 0), (3, 8, 2, 1), (3, 7, 0, 0)):
        data = model(t, a); co = [rng.randrange(p or 7) for _ in data[2]]
        cols = actual(data, co, p, True); far = components(data, cols, t, arm)
        _, core = core_matrix(data, cols, far, arm)
        sigma = len({v[arm+1]-u[arm+1] for v in core if v in data[4][arm] for u in data[4][1-arm]})
        planes = inventory(data, core, arm); assert len(planes) <= 1024*a*(t+sigma)
        patches = [check_patch(data, cols, far, arm, w, p) for w in weight_controls(cols, rng)]
        general.append({'t': t, 'a': a, 'characteristic': p, 'sigma': sigma, 'patches': patches})
    # One coefficient-independent arrangement is reused for two actual maps.
    t = 2; a = 5; nu = delta = 1; arm = 0; p = 101
    data, co, cols = build(t, a, arm, nu, delta, [1, 1], p, rng)
    far, marked, beta = orbit_partition(data, cols, t, a, arm, nu, delta)
    planes = inventory(data, potential_core(data, t, arm, nu, delta, marked), arm)
    weights = common_list(planes, a); sampler = []
    for gamma in ([1, 1], [0, 1]):
        data, co, cols = build(t, a, arm, nu, delta, gamma, p, rng)
        fan = complete_fan(cols, 4*a*t-2*t*t, p); observed = set()
        for w in weights:
            B = frozenset(field_basis(sorted(cols, key=lambda v: (dot(w, v), v)), cols.__getitem__, p)[0])
            assert all(dot(w, q) > 0 for q in normals(B, cols, p))
            observed.add(B)
        assert observed == set(fan)
        COUNTS['independent_common_sampler_full_fan_equalities'] += 1
        sampler.append({'coefficients': gamma, 'complete_original_fan_cones': len(fan), 'weights_reused': len(weights)})
    files = {'docs/hook_monomial_orbits.md', 'scripts/run_hook_monomial_orbit_checks.py',
             'docs/hook_bridge_reduction.md', 'docs/hook_tail_interface.md',
             'docs/hook_rank_departure.md', 'results/research/hook_tail_interface_checks.json'}
    pending = ['scripts/run_hook_monomial_orbit_checks.py']
    while pending:
        f = pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names = ([node.module] if isinstance(node, ast.ImportFrom) else
                     [v.name for v in node.names] if isinstance(node, ast.Import) else [])
            for name in names:
                target = 'scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:
                    files.add(target); pending.append(target)
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_hook_monomial_orbit_checks.py',
              'scope': 'Actual component bound for all compatible hooks, and a full-fan/common-list monomial-orbit class theorem. No unrestricted exponent, polynomial lower separation or exact novelty claimed.',
              'uniform_result': {'component_bound': 'O(n^2*(1+min(sigma_Y,sigma_Z)/t)^2)',
                                 'monomial_orbit_bound': 'O(n^2*beta^2)',
                                 'beta': 'min(ceil(t/delta),1+nu/delta,ceil(A/nu),A/t), or 1 if nu=0',
                                 'common_integer_weights': 'O(n^2*beta^2), height O(a^2)',
                                 'binomial_family_ranks': 'b=t^2, chi_Y=chi_Z=floor(t^2/4), rho=(A-1)*(t-1)',
                                 'mixed_direction_obstruction': 'W=(1+qv)^(t-1)+lambda*q^2*v, lambda!=0, A>=2t+2, char 0 or p>2t: chi_Y,chi_Z>=floor(t^2/4), sigma_Y,sigma_Z>=A-t-1, b=t^2, rho=(A-1)*(t-1)'},
              'counts': dict(sorted(COUNTS.items())), 'fan_counts': dict(sorted(FAN.items())),
              'original_quotient_replay_counts': dict(REPLAY), 'orbit_cases': cases,
              'mixed_direction_controls': mixed,
              'general_outer_map_controls': general, 'common_sampler': sampler,
              'limits': ['Chosen-arm outer components vanish in the polynomial-orbit class; the opposite arm is arbitrary.',
                         'The general support-component bound itself allows every compatible outer map.',
                         'Several independent monomial steps can join the chains and are not covered by the orbit bound.',
                         'Exact priority, sharpness and the unrestricted goal remain unresolved.'],
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_monomial_orbit_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'counts': report['counts'], 'source_hashes': len(files)}))


if __name__ == '__main__':
    main()
