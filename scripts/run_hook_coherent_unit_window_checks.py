#!/usr/bin/env python3
"""Original-field unit charts, coherent exchange windows and common weight lists."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import comb, gcd, lcm
from pathlib import Path
from random import Random
import ast
import hashlib
import json

from run_hook_bridge_checks import model, actual, normals, COUNTS as REPLAY
from run_hook_rank_departure_checks import rank, clean, plane, complete_fan, COUNTS as FAN
from run_hook_monomial_orbit_checks import common_list
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import dot

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def scalar(c, p):
    return c % p if p else Q(c)


def multiply(f, g, A, t, p):
    out = Counter()
    for (i, j), c in f.items():
        for (r, s), d in g.items():
            if i+r < A and j+s < t:
                out[i+r, j+s] += c*d
    return {v: scalar(c, p) for v, c in out.items() if scalar(c, p)}


def inverse(W, A, t, p):
    gamma = scalar(W.get((0, 0), 0), p); assert gamma
    gi = pow(int(gamma), p-2, p) if p else 1/gamma
    N = {v: scalar(-gi*c, p) for v, c in W.items() if v != (0, 0)}
    term = {(0, 0): scalar(1, p)}; out = Counter(); steps = 0
    while term:
        for v, c in term.items():
            out[v] += gi*c
        term = multiply(term, N, A, t, p); steps += 1
        assert steps <= A+t-1
    out = {v: scalar(c, p) for v, c in out.items() if scalar(c, p)}
    assert multiply(W, out, A, t, p) == {(0, 0): scalar(1, p)}
    COUNTS['independent_truncated_unit_inverse_products'] += 1
    return out


def build(t, a, arm, W, p, rng, symmetric=False):
    A = a-t; data = model(t, a)
    W = {v: scalar(c, p) for v, c in W.items() if v[0] < A and v[1] < t and scalar(c, p)}
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
    return data, co, actual(data, co, p, True), W


def row(data, t, arm, q, j):
    pair = (t+q, j) if arm == 0 else (j, t+q)
    return data[0].index(pair)


def lower(t, arm, q, j):
    return (0, t+q, t+j) if arm == 0 else (0, t+j, t+q)


def upper(t, arm, q, j):
    return (1, t+q, j) if arm == 0 else (1, j, t+q)


def transform(data, cols, t, a, arm, inv, p):
    A = a-t; out = {}; ri = {v: i for i, v in enumerate(data[0])}
    for v, column in cols.items():
        c = Counter()
        for i, value in column.items():
            pair = data[0][i]
            if pair[arm] < t:
                c[i] += value
                continue
            q, j = pair[arm]-t, pair[1-arm]
            for (r, s), factor in inv.items():
                if q+r < A and j+s < t:
                    target = (t+q+r, j+s) if arm == 0 else (j+s, t+q+r)
                    c[ri[target]] += value*factor
        out[v] = clean(c, p)
    return out


def charts(data, cols, t, a, arm, W, inv, p):
    trans = transform(data, cols, t, a, arm, inv, p)
    for q, j in product(range(a-t), range(t)):
        v = lower(t, arm, q, j); i = row(data, t, arm, q, j)
        projected = {r: c for r, c in trans[v].items() if data[0][r][arm] >= t}
        assert projected == {i: scalar(1, p)}
        if q >= t:
            assert trans[v] == projected and v in data[4][arm]
        else:
            assert v in data[4][2]
        COUNTS['independent_actual_central_column_inverse_identities'] += 1
    for v in data[4][2]:
        assert all(data[0][i][arm] < 2*t for i in trans[v] if data[0][i][arm] >= t)
    for w in ((1, 3, 17), (29, 5, 1), (1, 13, 2)):
        B = field_basis(sorted(cols, key=lambda v: (dot(w, v), v)), cols.__getitem__, p)[0]
        BB = field_basis(sorted(trans, key=lambda v: (dot(w, v), v)), trans.__getitem__, p)[0]
        assert B == BB and normals(B, cols, p) == normals(BB, trans, p)
        COUNTS['independent_whole_normal_cone_row_chart_equalities'] += 1
    return trans


def contraction_controls(data, cols, trans, t, a, arm, dp, dm, p):
    A = a-t; states = []; tail = {i for i, v in enumerate(data[0]) if v[arm] >= 2*t}
    yrows = {i for i, v in enumerate(data[0]) if v[arm] >= t}
    for is_upper, matrix, degree, summand in ((True, cols, dp, yrows), (False, trans, dm, tail)):
        for u, j in product((Q(1, 11), Q(1, 2), Q(1), Q(7)), (0, t-1)):
            q = A-1; e = upper(t, arm, q, j) if is_upper else lower(t, arm, q, j)
            if q < t or q*u+j <= (t-1+degree)*u+t-1:
                continue
            lam = Q(2*t+degree)*u+3*t
            w = (lam, u, Q(1)) if arm == 0 else (lam, Q(1), u)
            J = {i for i, v in enumerate(data[0]) if v[arm] >= t and
                 (is_upper or v[arm] >= 2*t) and (v[arm]-t)*u+v[1-arm] < q*u+j}
            chosen = []
            for i in J:
                pair = data[0][i]; q0, j0 = pair[arm]-t, pair[1-arm]
                v = upper(t, arm, q0, j0) if is_upper else lower(t, arm, q0, j0)
                assert matrix[v] == {i: scalar(1, p)} and dot(w, v) < dot(w, e)
                chosen.append(v)
            projected = {v: {i: c for i, c in col.items() if i not in J} for v, col in matrix.items()}
            for col in projected.values():
                assert not (set(col) & summand) or set(col) <= summand
            assert projected[e] and set(projected[e]) <= summand
            assert len(field_basis(chosen, matrix.__getitem__, p)[0]) == len(J)
            assert rank(tuple(projected), projected, p)+len(J) == rank(tuple(cols), cols, p)
            if not is_upper:
                COUNTS['nonidentity_original_source_columns_in_inverse_contractions'] += sum(len(cols[v]) > 1 for v in chosen)
            COUNTS['independent_strictly_cheaper_global_direct_sum_contractions'] += 1
            states.append({'upper': is_upper, 'u': str(u), 'input_short': j, 'contracted_rank': len(J)})
    return states


def price_controls(data, cols, trans, t, a, arm, dm, p):
    A = a-t; checks = 0
    for q, j, u in product(sorted({t, min(A-1, t+1), A-1}), (0, t-1), (Q(1, 3), Q(1), Q(5))):
        if q >= A or q < t:
            continue
        e = upper(t, arm, q, j); lam = 2*t-1+dm*u+1
        w = (lam, u, Q(1)) if arm == 0 else (lam, Q(1), u)
        out = Counter()
        assert trans[e]
        for i, c in trans[e].items():
            pair = data[0][i]; assert pair[arm] >= 2*t
            v = lower(t, arm, pair[arm]-t, pair[1-arm])
            assert v in data[4][arm] and dot(w, v) < dot(w, e)
            for r, d in cols[v].items():
                out[r] += c*d
        assert clean(out, p) == cols[e]
        COUNTS['independent_original_label_strict_price_loop_identities'] += 1
        checks += 1
    return checks


def inventory(data, t, arm, D):
    Y, Z = set(data[4][arm]), set(data[4][1-arm]); labels = tuple(v for g in data[4] for v in g)
    out = set()
    for i, v in enumerate(labels):
        for f in labels[:i]:
            e, z = (v, f) if v in Y and f in Z else (f, v)
            if e in Y and z in Z and e[arm+1] > 3*t+D and z[2-arm] > 6*t+2*D:
                continue
            out.add(plane(tuple(x-y for x, y in zip(v, f))))
    assert len(out) <= 4096*(max(max(v) for v in labels)+1)*(t+D)
    COUNTS['coefficient_independent_coherent_strip_plane_inventories'] += 1
    return out


def cross(a, b):
    return a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]


def rational_plane(q):
    d = lcm(*(x.denominator for x in q))
    return plane(tuple(int(x*d) for x in q))


def actual_fan_controls(data, cols, t, a, arm, dp, dm, p):
    D = max(dp, dm); planes = inventory(data, t, arm, D)
    labels = tuple(cols)
    all_planes = {plane(tuple(x-y for x, y in zip(v, f)))
                  for i, v in enumerate(labels) for f in labels[:i]}
    excluded = len(all_planes-planes)
    if t == 2 and a == 20:
        assert excluded > 0
        COUNTS['elongated_complete_fans_with_nonempty_excluded_plane_inventory'] += 1
    fan = complete_fan(cols, 4*a*t-2*t*t, p); Y, Z = set(data[4][arm]), set(data[4][1-arm])
    exchanges = 0; long_exchanges = 0; strict_rank_checks = 0
    for B, poly in fan.items():
        mids = []
        for left, right in zip(poly, poly[1:]+poly[:1]):
            if any(left[i] == right[i] == 0 for i in range(3)):
                continue
            assert rational_plane(cross(left, right)) in planes
            mids.append(tuple((x+y)/2 for x, y in zip(left, right)))
            COUNTS['complete_original_wall_plane_coverage_checks'] += 1
        _, piv = field_basis(sorted(B), cols.__getitem__, p)
        for f in cols:
            if f in B:
                continue
            rem, co = field_reduce(cols[f], piv, p); assert not rem
            for e, c in co.items():
                if not c or not ((e in Y and f in Z) or (f in Y and e in Z)):
                    continue
                for mid in mids:
                    if dot(mid, e) != dot(mid, f):
                        continue
                    chosen, opposite = (e, f) if e in Y else (f, e)
                    M, s = chosen[arm+1], chosen[2-arm]; K = opposite[2-arm]
                    assert M <= 3*t+D or K <= 6*t+2*D
                    q, j = M-t, s if chosen[0] else s-t
                    if q >= t:
                        short_weight = mid[2-arm]; u = mid[arm+1]/short_weight; lam = mid[0]/short_weight
                        degree = dp if chosen[0] else dm
                        assert q*u+j <= (t-1+degree)*u+t-1
                        if chosen[0]:
                            assert lam <= 2*t-1+dm*u
                        if u >= 1:
                            assert M <= 3*t+D
                        if u <= 1:
                            assert K <= 6*t+2*D
                        long_exchanges += 1
                    if strict_rank_checks < 16:
                        BB = (set(B)-{e})|{f}; assert rank(tuple(BB), cols, p) == len(B)
                        strict_rank_checks += 1
                    exchanges += 1
    COUNTS['independent_actual_coherent_cross_arm_exchanges'] += exchanges
    COUNTS['independent_actual_long_input_window_checks'] += long_exchanges
    COUNTS['independent_original_exchange_basis_rank_checks'] += strict_rank_checks
    return {'complete_original_cones': len(fan), 'coherent_cross_arm_exchanges': exchanges,
            'long_input_window_checks': long_exchanges, 'uniform_planes': len(planes),
            'potential_comparison_planes_excluded': excluded}


def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    previous = json.loads((ROOT/'results/research/hook_monomial_orbit_checks.json').read_text())
    assert previous['status'] == 'passed'
    for f, h in previous['source_sha256'].items():
        assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest() == h
    rng = Random(2026102327); cases = []
    specs = [(3, 24, 0, {(0,0):1,(1,1):1,(2,1):1}, 101, 2, False),
             (4, 20, 0, {**{(j,j):comb(3,j) for j in range(4)},(2,1):2}, 101, 2, False),
             (5, 18, 1, {(0,0):1,(2,1):1,(1,2):2,(3,2):1}, 3, 2, False),
             (4, 16, 0, {(0,0):2,(1,1):1,(0,1):2,(3,2):1}, 0, 2, False),
             (3, 11, 1, {(0,0):1,(1,0):1}, 5, None, False),
             (3, 16, 0, {(0,0):1,(5,1):1}, 2, 5, False),
             (3, 9, 0, {(0,0):1,(2,0):1}, 3, None, False),
             (3, 8, 1, {(0,0):1,(0,1):1,(0,2):1}, 2, 0, False),
             (3, 6, 0, {(0,0):1,(2,1):1}, 3, 2, False),
             (2, 8, 0, {(0,0):1,(1,1):1,(2,1):1}, 101, 2, True),
             (3, 10, 0, {(j,j):comb(2,j) for j in range(3)}|{(2,1):1}, 1000003, 2, True),
             (2, 20, 0, {(0,0):1,(1,1):1,(2,1):1}, 101, 2, 'symmetric')]
    for t, a, arm, W, p, C, full in specs:
        data, co, cols, W = build(t, a, arm, W, p, rng, symmetric=full == 'symmetric')
        inv = inverse(W, a-t, t, p); dp = max(i for i,j in W); dm = max(i for i,j in inv)
        if C is not None:
            assert all(i <= C*j for i,j in W) and max(dp, dm) <= C*(t-1)
            COUNTS['independent_forward_inverse_support_cone_controls'] += 1
        trans = charts(data, cols, t, a, arm, W, inv, p)
        states = contraction_controls(data, cols, trans, t, a, arm, dp, dm, p)
        prices = price_controls(data, cols, trans, t, a, arm, dm, p)
        fan = actual_fan_controls(data, cols, t, a, arm, dp, dm, p) if full else None
        cases.append({'t':t,'a':a,'arm':arm,'characteristic':p,'operator':[[i,j,str(c)] for (i,j),c in sorted(W.items())],
                      'forward_long_degree':dp,'inverse_long_degree':dm,'cone_slope_bound':C,
                      'opposite_mode':'same operator, zero outer maps' if full == 'symmetric' else 'arbitrary compatible coefficients',
                      'strict_contractions':states,'strict_price_identities':prices,'full_fan':fan})
        print('Unit-window case passed', t, a, arm, p, 'degrees', dp, dm, flush=True)
    # One geometric list for different fields and unrelated opposite coefficients.
    t=2; a=5; arm=0; D0=2; data=model(t,a); planes=inventory(data,t,arm,D0)
    weights=common_list(planes,a); samplers=[]
    for p,W in ((2,{(0,0):1,(1,1):1,(2,1):1}), (0,{(0,0):2,(0,1):1,(2,1):1})):
        data, co, cols, W=build(t,a,arm,W,p,rng)
        inv=inverse(W,a-t,t,p);assert max(i for i,j in inv)<=D0
        fan=complete_fan(cols,4*a*t-2*t*t,p);seen=set()
        for w in weights:
            B=frozenset(field_basis(sorted(cols,key=lambda v:(dot(w,v),v)),cols.__getitem__,p)[0])
            assert all(dot(w,q)>0 for q in normals(B,cols,p));seen.add(B)
        assert seen==set(fan)
        COUNTS['independent_common_integer_list_complete_fan_equalities']+=1
        samplers.append({'characteristic':p,'complete_original_cones':len(fan),'weights_reused':len(weights)})
    files={'docs/hook_coherent_unit_windows.md','scripts/run_hook_coherent_unit_window_checks.py',
           'docs/hook_bridge_reduction.md','docs/hook_monomial_orbits.md',
           'results/research/hook_monomial_orbit_checks.json'}
    pending=['scripts/run_hook_coherent_unit_window_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=([node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else [])
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:
                    files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_coherent_unit_window_checks.py',
            'scope':'Uniform all-field full fan and common weight list from one chosen unit arm, with its outer components zero and every opposite component arbitrary. No unrestricted endpoint or exact novelty claimed.',
            'uniform_result':{'D':'max(deg_q W,deg_q W^(-1)) in the truncated arm ring',
                              'full_fan_bound':'O(a^2*(t+D)^2)=O(n^2*(1+D/t)^2)',
                              'coherent_strips':'M<=3t+D or K<=6t+2D',
                              'support_cone_class':'q-degree<=C*short-degree implies D<=C*(t-1), hence O(n^2*(1+C)^2)',
                              'integer_compilation':'O(a^2*(t+D0)^2) common weights, height O(a^2)'},
            'counts':dict(sorted(COUNTS.items())),'fan_counts':dict(sorted(FAN.items())),
            'original_quotient_replay_counts':dict(REPLAY),'cases':cases,'common_sampler':samplers,
            'common_arrangement':{'planes':len(planes),'integer_weights':len(weights),'degree_cap':D0},
            'limits':['The chosen arm must be a unit and its two outer components must vanish.',
                      'Opposite-arm coefficients and both opposite outer components are arbitrary compatible values.',
                      'Vector contractions and invertible row charts keep original labels; they are not asserted to be algebra quotients or coordinate changes.',
                      'The mixed-direction large-rank family is quadratically bounded; its connected support was not an actual fan-count lower bound.',
                      'Strongest-current novelty, unrestricted sharp growth, polynomial lower separation and the original final paper remain unestablished.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_coherent_unit_window_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}),flush=True)


if __name__=='__main__':
    main()
