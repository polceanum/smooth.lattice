#!/usr/bin/env python3
"""Original-wall certificates for unbounded partners of one dual label."""
from collections import Counter
from fractions import Fraction as Q
from math import comb, factorial
from pathlib import Path
import ast
import hashlib
import json

from run_hook_filtered_annihilator_checks import assemble, bases, combination
from run_hook_transverse_prefix_checks import mul, clean, wall_weight
from run_hook_monomial_perturbation_checks import reflected
from run_hook_bridge_checks import normals, COUNTS as REPLAY
from run_signed_face_checks import field_basis, field_reduce
from run_generic_hook_checks import bareiss
from run_structural_complexity_checks import dot

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def value(c, p):
    return c % p if p else Q(c)


def binomial(D, r):
    return comb(D, r) if 0 <= r <= D else 0


def toeplitz(D, r, s, p):
    matrix = [[binomial(D, r+i-j) for j in range(s)] for i in range(s)]
    determinant = bareiss(matrix)
    product = Q(1)
    for i in range(s):
        product *= Q(factorial(D+i)*factorial(i),
                     factorial(r+i)*factorial(D-r+i))
    assert determinant == product and determinant > 0
    assert value(determinant, p)
    COUNTS['independent_integer_binomial_determinants'] += 1
    return str(determinant)


def operator_equivalence(data, original, control, G, A, p):
    """Extract N from original sources, independently of helper bases()."""
    t = 2; L = A-t; idx = data[0]
    inv = {(i, 0): (-1)**i*comb(A+i-2, i) for i in range(A)}
    assert mul(G, inv, A, t, p) == {(0, 0): value(1, p)}
    for i in range(A):
        for j in range(t):
            label = (0, t+j, t+i)
            images = []
            for cols in (original, control):
                h = Counter()
                for row, c in cols[label].items():
                    y, z = idx[row]
                    if z >= t:
                        h[z-t, y] += c
                h[i, j] -= 1
                images.append(clean(h, p))
            assert images[0] == mul(G, images[1], A, t, p)
            assert all(v[1] == 1 for v in images[0])
            COUNTS['independent_original_unit_output_equivalences'] += 1
    # The nonzero columns of N_Z=vG are exactly all v^0 inputs, and
    # their distinct lowest q terms prove independence for every subset.
    for length in (A, L):
        assert mul(G, inv, length, t, p) == {(0, 0): value(1, p)}
    COUNTS['invertible_full_and_truncated_output_maps'] += 2


def one_case(m, p):
    t = 2; A = 2*m*m+m+2; a = A+2; L = A-2
    U = {(0, 0): 1, (m, 1): 1}
    G = {(i, 0): comb(A-1, i) for i in range(A)}
    V = {(0, 0): 1, **{(i, 1): c for (i, _), c in G.items()}}
    data, co, cols = assemble(t, a, (U, V), p)
    _, co0, control = assemble(t, a, (U, {(0, 0): 1, (0, 1): 1}), p)
    operator_equivalence(data, cols, control, G, A, p)
    NY = {(m, 1): 1}; NZ = {(i, 1): c for (i, _), c in G.items()}
    E = (0, m+2, 3); target = (A-m-1, 0)
    records = []
    expectedY = {
        4: {(i, 0) for i in range(A-m)},
        5: {(i, 0) for i in range(2, A-m)},
        6: {(i, 0) for i in range(A-m)},
        7: {(i, 0) for i in range(L-m)}}
    expectedZ = {
        4: {(i, 0) for i in range(A)},
        5: {(i, 0) for i in range(2, A)},
        6: {(i, 0) for i in range(L)},
        7: {(i, 0) for i in range(L)}}
    # The two operators have constant selected sets for every positive
    # ratio. Sample both sides of internal order events as well as walls.
    for u in (Q(1, 2*A+1), Q(2, 3), Q(3, 2), Q(m)+Q(1, 3), Q(A+1)):
        yy = bases(NY, A, t, u, p)
        zz = bases(NZ, A, t, 1/u, p)
        z0 = bases({(0, 1): 1}, A, t, 1/u, p)
        for k in range(4, 8):
            assert set(yy[k]) == expectedY[k]
            assert set(zz[k]) == set(z0[k]) == expectedZ[k]
            COUNTS['identical_complete_greedy_basis_sets'] += 1
        assert not yy[-2] and not zz[-2] and not zz[-1]
        assert yy[-1] == {(A-m-2, 0), target}
    for j in range(2, m+1):
        K = m*m+j; J = A-1-K
        u = Q(2*K-1, 2*m); lam = (m+1)*u+1-J
        assert m < u < m+1 and 2 < lam < 2*u
        eps = min(Q(1, 4*m), (lam-2)/(4*(m+1)), (2*u-lam)/(4*(m-1)))
        q = (1, -(m+1), J-1)
        w = wall_weight((lam-(m+1)*eps, u-eps, Q(1)),
                        (lam+(m+1)*eps, u+eps, Q(1)), cols, q)
        lam, u, _ = w
        assert m < u < m+1 and K < m*u+1 < K+1 and 2 < lam < 2*u
        F = (1, 1, J+2); tau = dot(w, E)
        assert tau == dot(w, F)
        assert tau-lam > max(3*u+1, u+3)
        assert q not in data[7]
        # Independent integer determinants for both completing columns.
        determinants = [toeplitz(A-1, J, s, p) for s in (K, K+1)]
        rows = tuple(range(A))
        zcols = {(1, i): {i: value(1, p)} for i in rows}
        zcols.update({(0, i): {k+i: value(c, p) for (k, _), c in G.items()
                               if k+i < A and value(c, p)} for i in range(K+1)})
        zcheap = tuple((1, i) for i in range(J)) + tuple((0, i) for i in range(1, K+1))
        ZB, zp = field_basis(zcheap, zcols.__getitem__, p)
        assert len(ZB) == A-1
        assert field_reduce(zcols[0, 0], zp, p)[0] and field_reduce(zcols[1, J], zp, p)[0]
        _, zpf = field_basis((*ZB, (1, J)), zcols.__getitem__, p)
        functional = {i: field_reduce(zcols[1, i], zpf, p)[1].get((1, J), 0) for i in rows}
        gval = value(sum(functional.get(i, 0)*c for (i, _), c in G.items()), p)
        assert gval and functional[J] == 1
        cheap = tuple(v for v in cols if dot(w, v) < tau)
        B, piv = field_basis(cheap, cols.__getitem__, p)
        assert field_reduce(cols[E], piv, p)[0] and field_reduce(cols[F], piv, p)[0]
        _, pF = field_basis((*B, F), cols.__getitem__, p)
        # Compare every remaining Z coordinate and E as linear vectors,
        # not merely the ranks of the two exchanged columns.
        for i in rows:
            actual, rep = field_reduce(cols[1, 1, i+2], pF, p)
            assert not actual and rep.get(F, 0) == functional[i]
            COUNTS['entire_original_Z_quotient_coordinate_equalities'] += 1
        rem, rep = field_reduce(cols[E], pF, p)
        assert not rem and rep.get(F, 0) == value(-gval, p)
        index = {v: i for i, v in enumerate(data[0])}
        ell = {index[m+2, 1]: value(-gval, p)}
        ell.update({index[1, i+2]: c for i, c in functional.items() if c})
        ev = lambda v: value(sum(ell.get(i, 0)*c for i, c in cols[v].items()), p)
        assert all(not ev(v) for v in cheap)
        assert ev(E) and ev(F)
        KY = reflected(ell, data, A, t)[0]
        assert KY == {target: value(-gval, p)}
        assert not mul(NY, KY, L, t, p)
        assert mul(NY, KY, A, t, p) == {(A-1, 1): value(-gval, p)}
        assert target in bases(NY, A, t, u, p)[-1]
        COUNTS['same_dual_label_original_wall_certificates'] += 1
        # The control has exactly the same annihilator data, but E is
        # already in the original strictly cheaper span.
        relation = {E: 1, (0, 2, 2): -1, (1, 2, 0): 1,
                    (1, 0, 2): 1, (1, 1, 2): 1}
        assert not combination(control, relation, p)
        assert all(dot(w, v) < tau for v in relation if v != E)
        _, pc = field_basis(tuple(v for v in control if dot(w, v) < tau),
                            control.__getitem__, p)
        assert not field_reduce(control[E], pc, p)[0]
        COUNTS['identical_data_control_original_loop_certificates'] += 1
        values = sorted({dot(w, v) for v in cols})
        delta = min(lam/2, min(y-x for x, y in zip(values, values[1:]))/4)
        adjacent = []
        for sign in (-1, 1):
            ww = (lam+sign*delta, u, Q(1))
            S, _ = field_basis(sorted(cols, key=lambda v: (dot(ww, v), v)), cols.__getitem__, p)
            assert all(dot(ww, z) > 0 for z in normals(S, cols, p))
            adjacent.append(set(S)); COUNTS['original_adjacent_strict_minima'] += 1
        assert F in adjacent[0] and E not in adjacent[0]
        assert E in adjacent[1] and F not in adjacent[1]
        records.append({'j': j, 'K': K, 'J': J, 'weight': list(map(str, w)),
                        'lower': E, 'upper': F, 'direction': q,
                        'same_dual_label': target, 'integer_determinants': determinants})
        print('original partner wall passed', m, j, p, flush=True)
    assert len({tuple(v['upper']) for v in records}) == m-1
    return {'m': m, 'A': A, 'a': a, 'length': 8*A+8, 'field_characteristic': p,
            'same_order_interval': [m, m+1], 'partner_lower_bound': m-1,
            'coefficients': list(map(int, co)), 'control_coefficients': list(map(int, co0)),
            'walls': records}




def saturation_controls():
    for A, p in [(5, 0), (7, 0), (8, 101), (11, 101)]:
        G = [comb(A-1, i) for i in range(A)]
        cols = {(1, j): {j: value(1, p)} for j in range(A)}
        cols.update({(0, i): {i+j: value(c, p) for j, c in enumerate(G[:A-i])}
                     for i in range(A)})
        for J in range(A):
            for K in range(A):
                cheap = tuple((1, j) for j in range(J)) + tuple((0, i) for i in range(1, K+1))
                B, _ = field_basis(cheap, cols.__getitem__, p)
                rank = len(B)
                assert rank == min(A, J+K)
                ranks = []
                for added in [((0, 0),), ((1, J),), ((0, 0), (1, J))]:
                    ranks.append(len(field_basis((*B, *added), cols.__getitem__, p)[0])-rank)
                expected = [0, 0, 0] if J+K >= A else [1, 1, 1] if J+K == A-1 else [1, 1, 2]
                assert ranks == expected
                COUNTS['independent_all_prefix_saturation_patterns'] += 1
    largest = 0
    for A in range(7, 90):
        for nu in range(3, min(A-2, 25)):
            for r in range(1, 1+A//nu):
                eligible = []
                for K in range(nu*r+1, min(A, nu*(r+1)+1)):
                    lo = max(Q(r), Q(K-1, nu), Q(A-K, nu+1))
                    hi = min(Q(r+1), Q(K, nu), Q(A-K-2, nu-1))
                    if lo < hi:
                        eligible.append(K)
                        u = (lo+hi)/2
                        J = A-1-K
                        lam = (nu+1)*u+1-J
                        assert J > 0 and 2 < lam < 2*u and K < nu*u+1 < K+1
                assert len(eligible) <= nu
                assert len(eligible) <= 1+Q(2*A*nu, 4*nu*nu-1)
                largest = max(largest, len(eligible))
                COUNTS['exact_partner_interval_upper_bound_controls'] += 1
    return {'largest_partner_count_in_scalar_controls': largest,
            'scope': 'Exact prefix saturation and interval arithmetic; no new original-quotient wall is inferred from a scalar control.'}

def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    cases = [one_case(m, p) for m, p in [(3, 0), (3, 101), (4, 101), (5, 1009)]]
    saturation = saturation_controls()
    files = {'scripts/run_hook_exception_partner_checks.py', 'docs/hook_exception_partners.md',
             'docs/hook_filtered_annihilators.md', 'docs/hook_bridge_reduction.md',
             'docs/real_marginal_rank_fans.md'}
    pending = ['scripts/run_hook_exception_partner_checks.py']
    while pending:
        f = pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names = [node.module] if isinstance(node, ast.ImportFrom) else [v.name for v in node.names] if isinstance(node, ast.Import) else []
            for name in names:
                target = 'scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:
                    files.add(target); pending.append(target)
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_hook_exception_partner_checks.py',
              'scope': 'Sharp square-root local partner order within a binomial thickness-two class; identical four annihilator matroids and kernel spaces can suppress all selected exchanges. No full fan count, general mixed-operator upper bound or unrestricted exponent improvement.',
              'counts': dict(sorted(COUNTS.items())), 'original_quotient_replay_counts': dict(REPLAY),
              'cases': cases, 'saturation_and_local_upper_bound_controls': saturation,
              'limits': ['Coherent family: characteristic zero or p>2A.',
                         'The complete four annihilator data agree in every characteristic, but the binomial wall minors need not.',
                         'Selected walls are not a complete fan count.',
                         'The result is a local aggregate-budget obstruction, not a polynomial lower separation.',
                         'Exact priority and the original unrestricted endpoint remain open.'],
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_exception_partner_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'counts': report['counts'], 'source_hashes': len(files)}))


if __name__ == '__main__':
    main()
