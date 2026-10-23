#!/usr/bin/env python3
"""Exact polynomial-lattice windows and complete arbitrary-axis hook fans."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter, defaultdict
from bisect import bisect_right
from itertools import combinations, product
from math import gcd, lcm
from random import Random
import hashlib, json

from run_hook_axis_compilation_checks import (axis_factor, compile_tables, lookup,
    compiled_decision, hull, centroid, reference_controls)
from run_hook_axis_linear_checks import linear_decision, serialization
from run_hook_superregular_axis_checks import degree_factor
from run_hook_axis_birth_checks import direct_births, factor
from run_hook_bridge_checks import polygon
from run_hook_rank_departure_checks import plane
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import dot, area2

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def degree(P):
    return max(P, default=-10**9)


def subtract(P, R, shift, scalar, p):
    out = dict(P)
    for k, c in R.items():
        value = out.get(k+shift, 0)-scalar*c
        if p:
            value %= p
        if value:
            out[k+shift] = value
        else:
            out.pop(k+shift, None)
    return out


def multiply(P, R, p):
    out = {}
    for k, c in P.items():
        for j, d in R.items():
            out[k+j] = out.get(k+j, 0)+c*d
    return {k: (c % p if p else c) for k, c in out.items()
            if (c % p if p else c)}


def weighted_degree(row, theta):
    return max(theta+degree(row[0]), degree(row[1]))


def lattice(W, A, theta, p):
    """Two-row shifted reduction, independently of the finite rank compiler."""
    assert theta.denominator > 1
    rows = [({A: 1}, {}), (dict(W), {0: 1})]

    def position(row):
        return 0 if theta+degree(row[0]) > degree(row[1]) else 1

    steps = 0
    while position(rows[0]) == position(rows[1]):
        pos = position(rows[0])
        i, j = (0, 1) if degree(rows[0][pos]) >= degree(rows[1][pos]) else (1, 0)
        d = degree(rows[i][pos])-degree(rows[j][pos])
        v, c = rows[i][pos][degree(rows[i][pos])], rows[j][pos][degree(rows[j][pos])]
        scale = v*pow(c, -1, p) % p if p else Q(v)/c
        before = weighted_degree(rows[i], theta)
        rows[i] = tuple(subtract(a, b, d, scale, p) for a, b in zip(rows[i], rows[j]))
        assert weighted_degree(rows[i], theta) < before
        steps += 1
        assert steps <= 4*A+4
    rows.sort(key=position)
    r, s = degree(rows[0][0]), degree(rows[1][1])
    assert 0 <= r <= A and r+s == A
    determinant = subtract(multiply(rows[0][0], rows[1][1], p),
                           multiply(rows[0][1], rows[1][0], p), 0, 1, p)
    assert set(determinant) == {A} and determinant[A] in (1, -1, p-1)
    for a, b in rows:
        assert all(k >= A for k in subtract(a, multiply(W, b, p), 0, 1, p))
    COUNTS['independent_exact_lattice_determinants_and_congruences'] += 1
    COUNTS['maximum_two_row_reductions'] = max(COUNTS['maximum_two_row_reductions'], steps)
    return rows, r


def independent_least_upper_pivot(W, A, theta, p):
    """Finite linear equations for the first leading upper monomial."""
    for r in range(A):
        bdegree = min(A-1, int((theta+r).__floor__()))
        vectors = [{j: 1} for j in range(r)]
        vectors += [{q+j: c for j, c in W.items() if q+j < A}
                    for q in range(bdegree+1)]
        rank = len(field_basis(tuple(range(len(vectors))), vectors.__getitem__, p)[0])
        vectors.append({r: 1})
        if len(field_basis(tuple(range(len(vectors))), vectors.__getitem__, p)[0]) == rank:
            return r
    return A


def windows(W, A, t, p):
    pivots = []
    for k in range(-A, A):
        _, r = lattice(W, A, Q(2*k+1, 2), p)
        pivots.append(r)
    assert all(x >= y for x, y in zip(pivots, pivots[1:]))
    cuts = {-A, A} | set(range(-A, A+1, t))
    cuts |= {k for k in range(-A+1, A)
             if pivots[k+A-1]//t != pivots[k+A]//t}
    cuts = sorted(cuts)
    result = [{'lo': None, 'hi': -A, 'labels': set(), 'bucket': None}]
    for lo, hi in zip(cuts, cuts[1:]):
        R = t*(pivots[lo+A]//t)
        assert hi-lo <= t and all(R <= pivots[k+A] < R+t for k in range(lo, hi))
        # e=q-i for each active birth (delta,q), with the matching pivot type.
        labels = {(1, e) for e in range(R-t+1, R+2*t-1)}
        labels |= {(0, e) for e in range(A-R-2*t+2, A-R+t)}
        labels = {(d, e) for d, e in labels if -t < e < A}
        assert len(labels) <= 6*t
        result.append({'lo': lo, 'hi': hi, 'labels': labels, 'bucket': R})
    nu = min(W, default=A)
    high = {(1, e) for e in range(nu-t+1, min(A-1, nu+t-1)+1)} if W else set()
    result.append({'lo': A, 'hi': None, 'labels': high, 'bucket': None})
    assert len(result) <= (2*A+t-1)//t+A//t+3
    return {'states': result, 'breaks': [S['hi'] for S in result[:-1]], 'pivots': pivots}


def factor_control(W, A, t, p, independent=True):
    W = {q: (c % p if p else Q(c)) for q, c in W.items() if (c % p if p else c)}
    pure, heads = axis_factor(W, A, t)
    table = windows(W, A, t, p)
    for k in range(-A-1, A+1):
        theta = Q(2*k+1, 2)
        _, S = lookup(table, theta)
        _, _, birth = factor(pure, heads, theta, p)
        if -A <= k < A:
            rows, r = lattice(W, A, theta, p)
            short = {}
            for j, row in enumerate(rows):
                for shift in range(t):
                    price = (1, r+shift) if not j else (0, A-r+shift)
                    short[price] = {q+shift: c for q, c in row[1].items() if q+shift < t}
            gamma = direct_births(short, [{i: 1} for i in range(t)], theta, p)
            if independent and k in {-A, -1, A-1}:
                assert r == independent_least_upper_pivot(W, A, theta, p)
                COUNTS['independent_finite_equation_least_pivots'] += 1
        else:
            gamma = None
        for i in range(t):
            alpha = theta*birth[i][0]+birth[i][1] if i in birth else None
            if gamma is not None:
                assert gamma[i] == (max(Q(i), alpha) if alpha is not None else Q(i))
                COUNTS['short_lattice_vs_original_whole_prefix_births'] += 1
            if alpha is not None and alpha > i:
                delta, q = birth[i]
                assert (delta, q-i) in S['labels']
                if gamma is not None:
                    assert 0 <= q-(r if delta else A-r) < t
                COUNTS['active_births_in_matching_pivot_band'] += 1
        # Two further rational interior values verify the whole affine piece's band.
        for offset in (Q(1, 3), Q(2, 3)):
            theta2 = k+offset
            for i, (delta, q) in birth.items():
                if theta2*delta+q > i:
                    assert (delta, q-i) in S['labels']
                    COUNTS['whole_fine_interval_positive_label_band_checks'] += 1
    COUNTS['entire_axis_window_fixtures'] += 1
    return {'A': A, 't': t, 'characteristic': p,
            'operator': {str(q): str(c) for q, c in sorted(W.items())},
            'coarse_intervals_including_tails': len(table['states']),
            'upper_pivot_total_drop': table['pivots'][0]-table['pivots'][-1]}


def local_inventory(Y, Z, t):
    out = set()
    vectors = [{(d, e, -t*d) for d, e in Y['labels']},
               {(d, -t*d, e) for d, e in Z['labels']}]
    for v in vectors[0] | vectors[1]:
        if any(v):
            out.add(plane(v))
    for v, u in product(*vectors):
        q = tuple(a-b for a, b in zip(v, u))
        if any(q):
            out.add(plane(q))
    for arm, S in enumerate((Y, Z)):
        if S['lo'] is not None and S['hi'] is not None:
            for k in range(S['lo'], S['hi']+1):
                out.add(plane((1, -k, -t) if not arm else (1, -t, -k)))
        for endpoint in (S['lo'], S['hi']):
            if endpoint is not None:
                out.add(plane((1, -endpoint, -t) if not arm else (1, -t, -endpoint)))
    return out


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def small_integer_weight(P, a):
    rays = []
    for v in P:
        denominator = lcm(*(c.denominator for c in v))
        r = tuple(int(c*denominator) for c in v)
        d = gcd(*r)
        rays.append(tuple(c//d for c in r))
    triple = next(T for T in combinations(rays, 3) if dot(T[0], cross(T[1], T[2])))
    w = tuple(sum(r[k] for r in triple) for k in range(3))
    assert min(w) > 0 and max(w) <= 128*a*a
    return w


def walk_windows(compiled, factors):
    n, t, a = compiled['n'], compiled['t'], compiled['a']
    L = 2*n+1
    queue = [((1, L, L*L), None)]
    cells, areas, vertices, cache = {}, defaultdict(Q), defaultdict(list), {}
    while queue:
        w, previous = queue.pop()
        key, profile, N = linear_decision(compiled, w)
        indices = []
        for arm, (table, theta) in enumerate(zip(factors, [(Q(w[0])-t*w[2])/w[1], (Q(w[0])-t*w[1])/w[2]])):
            index, S = lookup(table, theta)
            indices.append(index)
            if S['lo'] is not None:
                N.add((1, -S['lo'], -t) if not arm else (1, -t, -S['lo']))
            if S['hi'] is not None:
                N.add((-1, S['hi'], t) if not arm else (-1, t, S['hi']))
        key = (tuple(indices), key)
        assert key != previous
        if key in cells:
            continue
        P = polygon(N)
        assert area2(P) > 0
        pair = tuple(indices)
        if pair not in cache:
            cache[pair] = local_inventory(*(factors[j]['states'][indices[j]] for j in range(2)), t)
        inventory = cache[pair]
        _, exact_profile, exact_N = compiled_decision(compiled, centroid(P))
        assert profile == exact_profile and all(dot(q, v) >= 0 for q in exact_N for v in P)
        sample = small_integer_weight(P, a)
        assert all(dot(sample, q) > 0 for q in N)
        assert linear_decision(compiled, sample)[1] == profile
        COUNTS['whole_window_patches_and_small_integer_samples'] += 1
        cells[key] = (profile, P)
        areas[profile] += area2(P)
        vertices[profile].extend(P)
        for v, u in zip(P, P[1:]+P[:1]):
            if any(v[k] == u[k] == 0 for k in range(3)):
                continue
            normal = next(q for q in N if dot(q, v) == dot(q, u) == 0)
            assert plane(normal) in inventory, (pair, normal)
            COUNTS['complete_refinement_facets_in_local_quadratic_inventory'] += 1
            nw, _ = exact_adjacent_weight(n, v, u, normal)
            queue.append((nw, key))
    assert sum(areas.values()) == 1
    fan = {profile: hull(points) for profile, points in vertices.items()}
    assert all(area2(fan[v]) == a for v, a in areas.items())
    assert sum(map(area2, fan.values())) == 1
    COUNTS['exact_full_window_covers_and_convex_profile_merges'] += 1
    return fan, cells, cache


def fan_control(A, t, Ws, p, full):
    Ws = [{q: (c % p if p else Q(c)) for q, c in W.items() if (c % p if p else c)} for W in Ws]
    compiled = compile_tables(t, A+t, [axis_factor(W, A, t) for W in Ws], p)
    factors = [windows(W, A, t, p) for W in Ws]
    fan, cells, inventory = walk_windows(compiled, factors)
    reference = reference_controls(t, A+t, Ws, p, compiled, fan, cells, full)
    COUNTS['fresh_entire_original_fans' if full else 'larger_original_spot_check_fans'] += 1
    digest = hashlib.sha256(json.dumps(serialization(fan), sort_keys=True).encode()).hexdigest()
    out = {'A': A, 't': t, 'length': compiled['n'], 'characteristic': p,
           'operators': [{str(q): str(c) for q, c in sorted(W.items())} for W in Ws],
           'coarse_factor_intervals': [len(F['states']) for F in factors],
           'nonempty_coarse_pairs': len(inventory), 'refinement_patches': len(cells),
           'actual_initial_ideal_cones': len(fan), 'exact_total_area': '1',
           'maximum_local_plane_inventory': max(map(len, inventory.values())),
           'canonical_profile_fan_sha256': digest, 'original_quotient_check': reference}
    print('Lattice window fan passed', A, t, p, len(cells), len(fan), flush=True)
    return out


def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    prior_path = 'results/research/hook_superregular_axis_checks.json'
    prior = json.loads((ROOT/prior_path).read_text())
    assert prior['status'] == 'passed'
    for f, h in prior['source_sha256'].items():
        assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest() == h
    rng = Random(20261023187)
    controls = []
    for p, (A, t), mode in product((0, 2, 3, 101), ((5, 2), (8, 3), (12, 4)), ('zero', 'sparse', 'dense', 'nonunit')):
        W = {} if mode == 'zero' else {q: rng.randrange(p or 7) for q in range(t if mode == 'nonunit' else 0, A)
                                      if mode != 'sparse' or rng.randrange(3) == 0}
        controls.append(factor_control(W, A, t, p))
    for p, A in ((2, 5), (3, 3)):
        for coefficients in product(range(p), repeat=A):
            controls.append(factor_control(dict(enumerate(coefficients)), A, 2, p, False))
    middle = [
        [335,109,93,521,395,930,652,814,88,760,844,787,314,559,385,632,322],
        [975,232,985,733,598,104,843,663,738,945,158,114,722,32,642,217,319]]
    for coefficients in middle:
        controls.append(factor_control(dict(enumerate(coefficients)), 32, 8, 1009))
    # This fully supported operator violates the old generic birth formula.
    W, A, t, p, theta = dict.fromkeys(range(4), 1), 6, 2, 2, Q(-1, 2)
    _, _, actual_birth = factor(*axis_factor(W, A, t), theta, p)
    _, generic_state = lookup(degree_factor(A, t, 3), theta)
    actual_price = sum(a*b for a, b in zip((theta, 1), actual_birth[1]))
    generic_price = sum(a*b for a, b in zip((theta, 1), generic_state['birth'][1]))
    assert actual_price == Q(7, 2) and generic_price == 3
    nongeneric_guard = {'A': A, 't': t, 'characteristic': p, 'head': 1,
        'theta': str(theta), 'operator': '1+q+q^2+q^3',
        'actual_birth': str(actual_price), 'superregular_formula_birth': str(generic_price),
        'scope': 'Actual rank degeneracy retained by the lattice proof, not identified with the generic degree table.'}
    COUNTS['nongeneric_full_support_birth_scope_guards'] += 1
    print('Entire lattice factor controls passed', len(controls), flush=True)
    cases = [fan_control(*spec) for spec in [
        (4, 2, [{0: 1, 1: 2, 2: 3}, {0: 2, 1: 1, 3: 2}], 0, True),
        (6, 2, [dict.fromkeys(range(4), 1), {0: 1, 2: 1, 4: 1}], 2, True),
        (5, 2, [{2: 1, 3: 2, 4: 1}, {}], 3, True),
        (5, 3, [{4: 1}, {1: 1, 3: 1}], 2, True),
        (6, 3, [{0: 1, 2: 1, 4: 1}, {0: 1, 3: 1, 5: 1}], 101, True),
        (8, 3, [{0: 1, 1: 1, 2: 1, 3: 1}, {2: 1, 4: 1, 6: 1}], 2, False),
        (13, 4, [{0: 1, 5: 2, 6: 3, 11: 5}, {0: 1, 7: 7, 8: 11, 12: 13}], 1009, False),
    ]]
    files = set(prior['source_sha256']) | {prior_path,
        'docs/hook_axis_lattice_windows.md', 'scripts/run_hook_axis_lattice_checks.py',
        'docs/hook_axis_linear_oracle.md', 'scripts/run_hook_axis_linear_checks.py'}
    report = {'schema_version': 1, 'status': 'passed',
        'command': 'python3 -B scripts/run_hook_axis_lattice_checks.py',
        'uniform_result': {'scope': 'Every zero-outer scalar-axis hook, A=a-t>=t>=2, over every field, including zero and nonunit operators.',
            'fan_and_complete_refinement_bound': 'O(A^2*t^2)=O(n^2), n=4at-2t^2.',
            'active_upper_label': 'q in r(theta)+{0,...,t-1}',
            'active_lower_label': 'q in A-r(theta)+{0,...,t-1}',
            'monotonicity': 'The least upper pivot degree r(theta) is nonincreasing, between zero and A.',
            'coarse_partition': 'O(A/t) intervals per factor; finite interval length at most t and pivot bucket width t.',
            'local_inventory': 'O(t) positive gain gradients, O(t^2) planes, O(t^4) arrangement cells per coarse pair.',
            'consequence': 'Coefficient-dependent O(n^2) complete positive integer weight list, height O(a^2).'},
        'counts': dict(COUNTS), 'entire_axis_controls': controls, 'complete_axis_fans': cases,
        'nongeneric_full_support_guard': nongeneric_guard,
        'limits': ['The elementary polynomial lattice reduction and predictable-degree property are classical.',
            'Finite checks support the uniform proof; they do not establish an asymptotic lower bound.',
            'Larger cases have original quotient spot checks, not independent complete original fan enumeration.',
            'No common coefficient-independent quadratic weight list or coefficient-independent profile formula is asserted.',
            'No unrestricted G_3(n) endpoint, current novelty certificate, bit-complexity or smooth-number speed claim.'],
        'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_axis_lattice_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'counts': dict(COUNTS), 'source_hashes': len(files)}), flush=True)


if __name__ == '__main__':
    main()
