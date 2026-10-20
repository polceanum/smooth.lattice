#!/usr/bin/env python3
"""Exact controls for higher-order hypersurface fans and cyclic price prefixes."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, gcd
from pathlib import Path
import hashlib
import json

from run_hypersurface_truncation_checks import (
    add, poly_mul, relation_rows, sparse_basis, sparse_basis_reduce,
)
from run_interleaving_research import minimal_excluded_border
from run_structural_complexity_checks import dot

ROOT = Path(__file__).resolve().parents[1]


def ground(s):
    return tuple((a, b, c) for a in range(s+1)
                 for b in range(s-a+1) for c in range(s-a-b+1))


def equation(r, z_square=True):
    out = {(2, 0, 0): 1, (0, 2, 0): 1, (1, 1, r): 1}
    if z_square:
        out[0, 0, 2] = 1
    return out


def monic_model(s, r, z_square=True):
    # Generate the small basis directly; never materialize the full s-simplex.
    basis = tuple((a, b, c) for a in range(min(1, s)+1)
                  for b in range(s-a+1) for c in range(s-a-b+1))
    index = {a: i for i, a in enumerate(basis)}

    @lru_cache(None)
    def evaluate(a):
        x, y, z = a
        if sum(a) > s:
            return {}
        if x < 2:
            return {index[a]: 1}
        out = defaultdict(int)
        tails = [(x-2, y+2, z), (x-1, y+1, z+r)]
        if z_square:
            tails.append((x-2, y, z+2))
        for b in tails:
            for k, c in evaluate(b).items():
                out[k] -= c
        return {k: c for k, c in out.items() if c}

    return basis, evaluate


def linear_sum(terms):
    out = defaultdict(int)
    for c, v in terms:
        for k, x in v.items():
            out[k] += c*x
    return {k: x for k, x in out.items() if x}


def rank_trace(points, evaluate, p, checkpoints):
    """Separate sparse Gaussian implementation; no selected-basis formulas."""
    pivots, selected, ranks = {}, set(), {}
    for step, a in enumerate(points, 1):
        vector = {k: (c % p if p else F(c)) for k, c in evaluate(a).items()
                  if (c % p if p else c)}
        for k, row in sorted(pivots.items()):
            c = vector.get(k, 0)
            if not c:
                continue
            for j, value in row.items():
                v = vector.get(j, 0)-c*value
                if p:
                    v %= p
                if v:
                    vector[j] = v
                else:
                    vector.pop(j, None)
        if vector:
            k = min(vector)
            inv = pow(vector[k], -1, p) if p else 1/vector[k]
            pivots[k] = {j: c*inv % p if p else c*inv for j, c in vector.items()}
            selected.add(a)
        if step in checkpoints:
            ranks[step] = len(pivots)
    return selected, ranks


def membership(polynomial, s, r, z_square=True):
    work, quotient = dict(polynomial), defaultdict(int)
    f = equation(r, z_square)
    while True:
        a = max((a for a, c in work.items() if c and a[0] >= 2 and sum(a) <= s),
                default=None)
        if a is None:
            break
        c, b = work[a], (a[0]-2, a[1], a[2])
        quotient[b] += c
        for exponent, value in f.items():
            v = add(b, exponent)
            work[v] = work.get(v, 0)-c*value
            if not work[v]:
                del work[v]
    assert all(sum(a) > s for a in work)
    # Independent multiplication replay of the returned ideal witness.
    replay = poly_mul(f, quotient, 0)
    assert all(polynomial.get(a, 0) == replay.get(a, 0)
               for a in set(polynomial) | set(replay) if sum(a) <= s)
    return {a: c for a, c in quotient.items() if c}


def border_certificate(B, evaluate, s, r, weight, z_square=True, save=False):
    chosen, pivots = sparse_basis(sorted(B), evaluate)
    assert len(chosen) == len(B)
    assert all(tuple(x-(i == j) for j, x in enumerate(a)) in B
               for a in B for i in range(3) if a[i])
    normals, records = set(), []
    total = [0, 0]
    for a in minimal_excluded_border(tuple(sorted(B))):
        rem, coefficients = sparse_basis_reduce(evaluate(a), pivots)
        assert not rem and all(c.denominator == 1 for c in coefficients.values())
        polynomial = {a: 1, **{b: -int(c) for b, c in coefficients.items()}}
        witness = membership(polynomial, s, r, z_square)
        for b in coefficients:
            z = tuple(x-y for x, y in zip(a, b))
            assert dot(z, weight) > 0
            normals.add(z)
        total[0] += 1
        total[1] += len(witness)
        if save:
            records.append({'head': a, 'tails': [{'exponent': b, 'coefficient': int(c)}
                                               for b, c in sorted(coefficients.items())],
                            'multiplier': [{'exponent': b, 'coefficient': c}
                                           for b, c in sorted(witness.items())]})
    return normals, records, total


def model_controls():
    counts = defaultdict(int)
    for s, r in product((0, 1, 2, 4, 6, 8), (1, 2, 3, 5, 9)):
        basis, evaluate = monic_model(s, r)
        assert len(basis) == (s+1)**2
        maps = [[evaluate(add(a, tuple(int(i == j) for i in range(3))))
                 for a in basis] for j in range(3)]
        assert max(len(c) for M in maps for c in M) <= 3
        for i, j in ((0, 1), (0, 2), (1, 2)):
            for k in range(len(basis)):
                assert linear_sum((c, maps[i][v]) for v, c in maps[j][k].items()) == \
                       linear_sum((c, maps[j][v]) for v, c in maps[i][k].items())
                counts['integral_commutator_columns'] += 1
        for a in ground(s):
            for v in evaluate(a):
                z = tuple(x-y for x, y in zip(a, basis[v]))
                assert sum(z) % r == 0
                gamma = -sum(z)//r
                assert (z[1]+gamma) % 2 == (z[2]+r*gamma) % 2 == 0
                counts['graded_column_terms'] += 1
        K = (s+r+2)//(r+2)
        power = {(2*a, 2*b, 2*(K-a-b)): comb(K, a)*comb(K-a, b)
                 for a in range(K+1) for b in range(K-a+1)}
        assert not linear_sum((c, evaluate(a)) for a, c in power.items())
        counts['annihilating_homogeneous_power_identities'] += 1
        counts['integral_models'] += 1
    fixtures = []
    for (s, r), p in product(((4, 1), (4, 2), (6, 3), (6, 4), (8, 2), (8, 5)),
                             (0, 2, 3, 5, 101)):
        basis, evaluate = monic_model(s, r)
        E = ground(s)
        rows = relation_rows(E, [equation(r)], p)
        assert len(rows) == len(E)-len(basis)
        columns = [evaluate(a) for a in E]
        assert all(any(c % p if p else c for c in column.values()) for column in columns)
        counts['canonical_nonzero_monomial_columns'] += len(columns)
        for row in rows:
            v = linear_sum((c, column) for c, column in zip(row, columns) if c)
            assert all(not (c % p if p else c) for c in v.values())
            counts['independent_macaulay_kernel_rows'] += 1
        # The monic coordinates include every coordinate unit, so equal dimensions
        # and annihilation of the full relation space prove kernel equality.
        assert all(evaluate(a) == {i: 1} for i, a in enumerate(basis))
        counts['independent_full_kernel_equalities'] += 1
        fixtures.append({'cutoff': s, 'order': r, 'characteristic': p,
                         'parent_size': len(E), 'length': len(basis)})
    # Check the two independent parameters, including the growing-order regime.
    for s in range(3, 401):
        for r in range(1, s-1):
            K = (s+r+2)//(r+2)
            D = min(s, 2*K-1)
            assert D*r <= 3*s
            counts['uniform_slab_arithmetic_checks'] += 1
    return dict(counts), fixtures


def boundary_basis(s, r, mu):
    lam = r*mu
    chosen = [(0, 0, c) for c in range(s+1)]
    for k in range(1, s+1):
        L = s-k+1
        for j in range(min(r, L)):
            Lj = 1+(L-1-j)//r
            arms = [(a, k-a, j) for a in range(k+1)]
            arms += [(0, k, j+r*d) for d in range(1, Lj)]
            arms += [(k, 0, j+r*d) for d in range(1, Lj)]
            chosen += sorted(arms, key=lambda a: (lam*a[0]+a[2], a))[:2*Lj]
    assert len(chosen) == len(set(chosen)) == (s+1)**2
    return frozenset(chosen)


def wall_formula(s, r):
    return tuple(sorted({F(p, q) for q in range(2, 2*(s+r)//(r+2)+1)
                         for p in range(1, q) if gcd(p, q) == 1 and
                         (1 if p % 2 and q % 2 else 2)*((r+2)*q+r*p) <= 2*(s+r)}))


def boundary_controls():
    counts, fixtures = defaultdict(int), []
    for s in (4, 6, 8, 10, 12):
        for r in sorted({1, 2, 3, 5, s-2}):
            if r > s-2:
                continue
            _, evaluate = monic_model(s, r, False)
            delta = F(1, 8*(s+1)**2)
            mu, left, walls, states = delta, F(0), [], set()
            while mu < 1:
                B = boundary_basis(s, r, mu)
                weight = (1+r*mu, F(1), F(1))
                normals, _, stats = border_certificate(B, evaluate, s, r, weight, False)
                lo, hi = F(0), F(1)
                for x, y, z in normals:
                    assert x+y == 0
                    if x > 0:
                        lo = max(lo, F(-z, r*x))
                    elif x < 0:
                        hi = min(hi, F(-z, r*x))
                    else:
                        assert z >= 0
                assert lo == left and lo < mu < hi
                assert B not in states
                states.add(B)
                counts['boundary_entire_interval_certificates'] += 1
                counts['boundary_integer_membership_witnesses'] += stats[0]
                # Independent field Gaussian replay of each continuous interval.
                for p in (0, 2, 3, 5):
                    chosen, _ = rank_trace(sorted(ground(s), key=lambda a: (dot(a, weight), a)),
                                           evaluate, p, set())
                    assert chosen == B
                    counts['independent_boundary_field_bases'] += 1
                if hi < 1:
                    walls.append(hi)
                left, mu = hi, hi+delta
            assert left == 1 and tuple(walls) == wall_formula(s, r)
            central = states | {frozenset((b, a, c) for a, b, c in B) for B in states}
            assert len(central) == 2*len(states)
            outer = {frozenset(a for a in ground(s) if a[i] < 2) for i in (0, 1)}
            assert len(central | outer) == 2*len(walls)+4
            fixtures.append({'cutoff': s, 'order': r, 'length': (s+1)**2,
                             'boundary_fan_cells': 2*len(walls)+4,
                             'positive_normalized_walls': list(map(str, walls))})
    return dict(counts), fixtures


def plane_obstruction():
    s, r, weight = 10, 3, (419, 68, 221)
    basis, evaluate = monic_model(s, r)
    B, _ = sparse_basis(sorted(ground(s), key=lambda a: (dot(a, weight), a)), evaluate)
    B = frozenset(B)
    normals, records, stats = border_certificate(B, evaluate, s, r, weight, save=True)
    assert (3, -3, -3) in normals and (0, -2, 2) in normals
    assert min(dot(z, weight) for z in normals) > 0
    counts = {'all_field_original_cone_certificates': 1,
              'original_monic_border_membership_witnesses': stats[0]}
    return counts, {'cutoff': s, 'order': r, 'length': len(basis),
                    'strict_weight': weight, 'standard_set': sorted(B),
                    'border_relations': records, 'normals': sorted(normals),
                    'minimum_strict_slack': min(dot(z, weight) for z in normals),
                    'separating_normals': [(3, -3, -3), (0, -2, 2)],
                    'valid_over_every_field': True}


def g_power_coefficient(m, exponent):
    if min(exponent) < 0 or sum(exponent) != 2*m or any(v % 2 for v in exponent):
        return 0
    a, b, _ = (v//2 for v in exponent)
    return comb(m, a)*comb(m-a, b)


def excluded_matrix(m, inclusive=False):
    D = 2*m+2
    bound = 4 if inclusive else 5
    rows = tuple((a, b, D-a-b) for a in range(3) for b in range(bound+1)
                 if 2*a+b <= bound)
    h = (1, 1, 2*m-1)
    quadratics = ((2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0), (0, 1, 1), (0, 0, 2))
    matrix = [[int(a == add(h, tuple(int(i == j) for j in range(3)))) for i in range(3)]
              +[g_power_coefficient(m, tuple(x-y for x, y in zip(a, b))) for b in quadratics]
              for a in rows]
    return rows, matrix


def smith_units(matrix):
    """Unimodular unit elimination; return all nonunit residual entries."""
    A = [list(row) for row in matrix]
    rank = 0
    while rank < min(len(A), len(A[0])):
        pivot = next(((i, j) for i in range(rank, len(A)) for j in range(rank, len(A[0]))
                      if abs(A[i][j]) == 1), None)
        if pivot is None:
            break
        i, j = pivot
        A[i], A[rank] = A[rank], A[i]
        for row in A:
            row[j], row[rank] = row[rank], row[j]
        if A[rank][rank] < 0:
            A[rank] = [-v for v in A[rank]]
        for i in range(len(A)):
            if i != rank:
                c = A[i][rank]
                A[i] = [x-c*y for x, y in zip(A[i], A[rank])]
        for j in range(len(A[0])):
            if j != rank:
                c = A[rank][j]
                for row in A:
                    row[j] -= c*row[rank]
        rank += 1
    return rank, [row[rank:] for row in A[rank:]]


def cyclic_controls():
    counts, smith, original = defaultdict(int), [], []
    for m in (*range(2, 31), 64, 81, 210, 1009, 10**9+7):
        rows, matrix = excluded_matrix(m)
        units, residual = smith_units(matrix)
        assert units == 8 and len(residual[0]) == 1
        assert gcd(*(v[0] for v in residual)) == m
        # The beta column is the cyclic generator, with relation exactly m*beta.
        difference = [row[4]-row[2] for row in matrix]
        assert difference == [m*int((a, b) == (1, 3)) for a, b, c in rows]
        _, inclusive = excluded_matrix(m, True)
        inc_units, inc_residual = smith_units(inclusive)
        assert inc_units == 7 and not any(v for row in inc_residual for v in row)
        counts['exact_cyclic_smith_matrices'] += 1
        counts['exact_inclusive_unit_matrices'] += 1
        smith.append({'m': m, 'strict_shape': [12, 9], 'unit_invariants': 8,
                      'remaining_invariant': m, 'free_cokernel_rank': 3,
                      'inclusive_unit_invariants': 7})
    for m in (2, 3, 4, 5, 6, 7):
        s, r, D = m*(2*m+1), 2*m-1, 2*m+2
        basis, evaluate = monic_model(s, r)
        beta = (1, 3, 2*m-2)
        weight = (F(1), 1+F(1, 8*m), 1+F(1, 4*m))
        T = dot(beta, weight)
        all_points = ground(D)
        cheap = sorted((a for a in all_points if dot(a, weight) < T),
                       key=lambda a: (dot(a, weight), a))
        tier = sorted(a for a in all_points if dot(a, weight) == T and a != beta)
        rest = sorted(a for a in all_points if dot(a, weight) > T)
        points = cheap+tier+[beta]+rest
        assert len(cheap) == comb(2*m+5, 3)-12 and len(tier) == 2
        N = len(points)
        for p in (0, 2, 3, 5, 7, 11):
            selected, trace = rank_trace(points, evaluate, p,
                                         {len(cheap), len(cheap)+2, len(cheap)+3, N})
            bad = bool(p and m % p == 0)
            assert trace[len(cheap)] == N-16-int(bad)
            assert trace[len(cheap)+3] == N-15
            assert trace[N] == N-13
            assert (beta in selected) == bad
            assert trace[len(cheap)+3]-trace[len(cheap)+2] == int(bad)
            counts['independent_original_cyclic_prefixes'] += 1
            original.append({'m': m, 'cutoff': s, 'order': r, 'length': len(basis),
                             'characteristic': p, 'strict_prefix_size': len(cheap),
                             'strict_prefix_rank': trace[len(cheap)],
                             'inclusive_prefix_rank': trace[len(cheap)+3],
                             'beta_standard': bad})
    return dict(counts), {'smith_fixtures': smith, 'original_column_fixtures': original,
                          'scope': 'Integer unit operations on constant-size excluded coefficient matrices, independently compared with sparse Gaussian ranks of original monomials in the full integral quotient. Composite m and primes not dividing m are included.'}


def main():
    if not __debug__:
        raise RuntimeError('Run without -O: exact assertions are required.')
    model_counts, models = model_controls()
    boundary_counts, boundary = boundary_controls()
    cone_counts, cone = plane_obstruction()
    cyclic_counts, cyclic = cyclic_controls()
    counts = {**model_counts, **boundary_counts, **cone_counts, **cyclic_counts}
    sources = ('docs/higher_interaction_fans.md', 'scripts/run_higher_interaction_checks.py',
               'scripts/run_hypersurface_truncation_checks.py', 'scripts/run_interleaving_research.py',
               'scripts/run_structural_complexity_checks.py', 'scripts/run_socle_stability_checks.py',
               'docs/hypersurface_truncations.md', 'docs/profile_lattice_bound.md',
               'docs/cubic_fan_classification.md')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_higher_interaction_checks.py',
              'counts': counts, 'independent_kernel_models': models,
              'boundary_interval_fixtures': boundary, 'equal_xy_plane_obstruction': cone,
              'cyclic_prefix_controls': cyclic,
              'uniform_results': {
                  'family_fan_bound': 'O(n^(9/4)/r), n=(s+1)^2, uniform in s,r and field',
                  'excluded_amplification_regime': 'r>=sqrt(s) implies O(n^2) full positive fan',
                  'boundary_fan': '2*M_(s,r)+4, with epsilon*((r+2)*q+r*p)<=2*(s+r)',
                  'cyclic_prefix_quotient': 'For every integer m>=2, V_D/L_< is Z^3 plus Z/m; beta changes standard membership exactly in characteristics dividing m.',
                  'integer_algebra_itself_free': True},
              'scope': 'Finite independent exact controls support the uniform proofs in the note. A changed standard set is not a count increase. The small-order full fan, unrestricted G_3, current priority and publication significance remain open.',
              'unrestricted_G3_improvement': False, 'novelty_established': False,
              'goal_complete': False,
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    (ROOT/'results/research/higher_interaction_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(counts, indent=2))


if __name__ == '__main__':
    main()
