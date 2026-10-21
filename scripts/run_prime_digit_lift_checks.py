#!/usr/bin/env python3
"""Exact Fourier lifts, integral dyadic lifts, and unchanged quotient fans."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations
from math import comb
from pathlib import Path
import hashlib
import json

from run_adjacent_rank_checks import rank
from run_oracle_enumeration_checks import digit_rank
from run_projection_tag_checks import compare_fans, fan
from run_signed_face_checks import normals as border_normals

ROOT = Path(__file__).resolve().parents[1]


def frequencies(n, M, p):
    if not n:
        return ()
    if n == M:
        return tuple(range(M))
    h = M//p
    q, r = divmod(n, h)
    out = [a+p*b for a in range(q) for b in range(h)]
    out += [q+p*b for b in frequencies(r, h, p)]
    assert len(set(out)) == n and all(0 <= a < M for a in out)
    return tuple(sorted(out))


@lru_cache(None)
def signed_permutations(n):
    return tuple((a, (-1)**sum(a[i] > a[j] for i in range(n) for j in range(i+1, n)))
                 for a in permutations(range(n)))


def root_determinant(rows, cols, M, p):
    """Exact determinant in Q[T]/Phi_(p^s)(T), no numeric roots.

    Phi_(p^s)=1+T^h+...+T^((p-1)h). Leibniz expansion is used
    only for at most four rows in the fixed exhaustive fixtures.
    """
    assert len(rows) == len(cols) <= 4
    out = [0]*M
    for perm, sign in signed_permutations(len(rows)):
        out[sum(a*cols[j] for a, j in zip(rows, perm)) % M] += sign
    h, d = M//p, M-M//p
    for i in range(M-1, d-1, -1):
        value = out[i]
        for j in range(p):
            out[i-d+j*h] -= value
    assert not any(out[d:])
    return tuple(out[:d])


def cyclotomic_norm(poly, M, p):
    h, d = M//p, M-M//p
    columns = []
    for j in range(d):
        out = [0]*j+list(poly)+[0]*(d-1-j)
        for i in range(len(out)-1, d-1, -1):
            value = out[i]
            for k in range(p):
                out[i-d+k*h] -= value
        assert not any(out[d:])
        columns.append(out[:d])
    return determinant([[columns[j][i] for j in range(d)] for i in range(d)])


def exact_fourier_controls(counts):
    fixtures = []
    for p, M in ((2, 2), (3, 3), (2, 4), (5, 5), (2, 8), (3, 9)):
        for n in range(1, M+1):
            I = frequencies(n, M, p)
            bases = []
            pascal = [tuple(comb(a, j) % p for j in range(n)) for a in range(M)]
            base_cols = tuple(range(n)) if n <= M//2 else tuple(range(n, M))
            base_rows = I if n <= M//2 else tuple((-i) % M for i in range(M) if i not in I)
            base_norm = abs(cyclotomic_norm(root_determinant(base_rows, base_cols, M, p), M, p))
            assert base_norm
            for B in combinations(range(M), n):
                # Jacobi's complementary-minor identity for the full Fourier
                # matrix reduces the expansion size when n>M/2.
                if n <= M//2:
                    rows, cols = I, B
                else:
                    rows = tuple((-i) % M for i in range(M) if i not in I)
                    cols = tuple(i for i in range(M) if i not in B)
                root_minor = root_determinant(rows, cols, M, p)
                nonzero = any(root_minor)
                direct = rank([pascal[a] for a in B], p) == n
                assert nonzero == direct
                counts['exact_cyclotomic_vs_Pascal_maximal_minors'] += 1
                if nonzero:
                    assert abs(cyclotomic_norm(root_minor, M, p)) == base_norm
                    counts['independent_cyclotomic_unit_basis_minor_norms'] += 1
                    bases.append(sum(1 << a for a in B))
            for mask in range(1 << M):
                S = tuple(a for a in range(M) if mask >> a & 1)
                fourier_rank = max((mask & B).bit_count() for B in bases)
                actual = rank([pascal[a] for a in S], p)
                assert fourier_rank == actual == digit_rank(n, M, p, S)
                counts['exhaustive_exact_Fourier_subset_rank_equalities'] += 1
            fixtures.append({'prime': p, 'period': M, 'rank': n,
                             'frequencies': I, 'basis_masks': bases,
                             'reference_Fourier_minor_absolute_norm': base_norm,
                             'all_nonzero_power_basis_minor_absolute_norms': 1})
    return fixtures


def modular_fourier_controls(counts):
    fixtures = []
    for p, M, ell in ((2, 16, 17), (2, 32, 97), (3, 27, 109), (5, 25, 101)):
        assert (ell-1) % M == 0
        zeta = next(pow(a, (ell-1)//M, ell) for a in range(2, ell)
                    if pow(pow(a, (ell-1)//M, ell), M//p, ell) != 1)
        assert pow(zeta, M, ell) == 1 and pow(zeta, M//p, ell) != 1
        for n in sorted({1, M//p, M//p+1, M//2, M-1, M}):
            I = frequencies(n, M, p)
            columns = [tuple(pow(zeta, i*a, ell) for i in I) for a in range(M)]
            subsets = {(), tuple(range(M)), tuple(range(n)), tuple(range(0, M, p))}
            for seed in range(160):
                subsets.add(tuple(a for a in range(M)
                                  if ((a+1)*(seed+3)+a*a+seed*seed) % 37 < seed % 36))
            for S in sorted(subsets):
                assert rank([columns[a] for a in S], ell) == digit_rank(n, M, p, S)
                counts['independent_cross_characteristic_Fourier_rank_equalities'] += 1
            fixtures.append({'prime': p, 'period': M, 'rank': n,
                             'lift_characteristic': ell, 'primitive_root': zeta,
                             'frequencies': I, 'subsets': len(subsets)})
    return fixtures


def dyadic_columns(n, M):
    f = [int(j & n == j) for j in range(n+1)]
    columns, v = [], [1]+[0]*(n-1)
    for _ in range(M):
        columns.append(tuple(v))
        coefficient = v[-1]
        v = [0]+v[:-1]
        for j in range(n):
            v[j] -= coefficient*f[j]
    return f, columns


def determinant(A):
    """Integer Bareiss elimination, independent of finite-field rank."""
    A = [list(row) for row in A]
    sign, previous = 1, 1
    for j in range(len(A)-1):
        k = next((i for i in range(j, len(A)) if A[i][j]), None)
        if k is None:
            return 0
        if k != j:
            A[j], A[k] = A[k], A[j]
            sign = -sign
        pivot = A[j][j]
        for i in range(j+1, len(A)):
            for k in range(j+1, len(A)):
                value = pivot*A[i][k]-A[i][j]*A[j][k]
                assert value % previous == 0
                A[i][k] = value//previous
            A[i][j] = 0
        previous = pivot
    return sign*A[-1][-1] if A else 1


def dyadic_minor_controls(counts):
    fixtures = []
    for M in (2, 4, 8, 16):
        for n in range(1, M+1):
            f, columns = dyadic_columns(n, M)
            histogram = Counter()
            for B in combinations(range(M), n):
                d = determinant([[columns[a][j] for a in B] for j in range(n)])
                assert d in (-1, 0, 1)
                histogram[d] += 1
                # This reduction is a triangular change of basis from
                # powers X^a to (1+t)^a modulo t^n.
                expected = rank([tuple(comb(a, j) & 1 for j in range(n)) for a in B], 2)
                assert bool(d) == (expected == n)
                counts['independent_integral_unimodular_maximal_minors'] += 1
            if M <= 8:
                for size in range(1, n+1):
                    for rows in combinations(range(n), size):
                        for B in combinations(range(M), size):
                            assert determinant([[columns[a][j] for a in B] for j in rows]) in (-1, 0, 1)
                            counts['exhaustive_all_size_integral_minors'] += 1
            fixtures.append({'rank': n, 'period': M, 'polynomial': f,
                             'maximal_minor_histogram': dict(sorted(histogram.items()))})
    return fixtures


def dyadic_fan_controls(counts):
    fixtures = []
    for n, k in ((1, 1), (3, 2), (4, 3), (5, 3), (10, 3), (10, 5)):
        M = 1
        while M < n:
            M *= 2
        f, columns = dyadic_columns(n, M)
        def mul(left, right):
            out = [0]*(2*n-1)
            for i, a in enumerate(left):
                for j, b in enumerate(right):
                    out[i+j] += a*b
            for i in range(len(out)-1, n-1, -1):
                value = out[i]
                for j in range(n+1):
                    out[i-n+j] -= value*f[j]
            return tuple(out[:n])
        inverse = tuple(-f[j+1] for j in range(n))
        direct, backwards = (1,)+(0,)*(n-1), (1,)+(0,)*(n-1)
        for e in range(3*M+1):
            for exponent, value in ((e, direct), (-e, backwards)):
                quotient, residue = divmod(exponent, M)
                sign = -1 if n == M and quotient % 2 else 1
                assert value == tuple(sign*c for c in columns[residue])
                counts['independent_signed_Laurent_period_replays'] += 1
            # X modulo F_n is the second column when n>1; at n=1
            # it is -1. Multiplication itself uses polynomial reduction.
            x = (0, 1)+(0,)*(n-2) if n > 1 else (-1,)
            direct, backwards = mul(direct, x), mul(backwards, inverse)
        @lru_cache(None)
        def ev(a):
            exponent = a[0]-a[1]-k*a[2]
            quotient, residue = divmod(exponent, M)
            sign = -1 if n == M and quotient % 2 else 1
            return {j: sign*c for j, c in enumerate(columns[residue]) if c}
        rational = fan(n, ev, 0, counts)
        for p in (2, 3, 5):
            actual = fan(n, ev, p, counts)
            compare_fans(actual, rational, counts)
            for B in rational:
                _, relations = border_normals(tuple(sorted(B)), ev, 0)
                _, reduced = border_normals(tuple(sorted(B)), ev, p)
                assert [r['head'] for r in relations] == [r['head'] for r in reduced]
                for left, right in zip(relations, reduced):
                    terms = {tuple(t['exponent']): Q(t['coefficient']) for t in left['tails']}
                    assert all(c in (-1, 1) for c in terms.values())
                    right_terms = {tuple(t['exponent']): int(t['coefficient']) for t in right['tails']}
                    assert right_terms == {a: int(c) % p for a, c in terms.items()}
                    counts['independent_integral_reduced_basis_specializations'] += 1
        # Recheck characteristic two against the original Pascal columns,
        # rather than only against a reduction of the integer lift.
        @lru_cache(None)
        def pascal(a):
            residue = (a[0]-a[1]-k*a[2]) % M
            return {j: 1 for j in range(n) if comb(residue, j) & 1}
        compare_fans(fan(n, pascal, 2, counts), rational, counts)
        fixtures.append({'rank': n, 'k': k, 'period': M, 'polynomial': f,
                         'initial_ideals_in_every_tested_field': len(rational),
                         'characteristics': [0, 2, 3, 5],
                         'complete_cones': [{'basis': sorted(B),
                              'cone': [list(map(str, v)) for v in poly]}
                                           for B, poly in rational.items()]})
    return fixtures


def odd_prime_controls(counts):
    fixtures = []
    for p in (3, 5, 7):
        n, M = p+1, p*p
        columns = [tuple(comb(a, j) % p for j in range(n)) for a in range(M)]
        C = tuple(range(1, p))
        S = tuple(p*b for b in range(p))+(p+1,)
        assert rank([columns[a] for a in C], p) == p-1
        for a, b in combinations(S, 2):
            assert rank([columns[j] for j in C+(a, b)], p) == p+1
            counts['independent_odd_prime_uniform_minor_pairs'] += 1
        assert rank([columns[j] for j in C+S], p) == p+1
        fixtures.append({'prime': p, 'rank': n, 'period': M,
                         'contracted_labels': C, 'uniform_rank_two_labels': S,
                         'minor': 'U_(2,p+1)',
                         'binary_prime_field_representation_possible': False})
    return fixtures


def main():
    counts = defaultdict(int)
    exact = exact_fourier_controls(counts)
    modular = modular_fourier_controls(counts)
    integral = dyadic_minor_controls(counts)
    fans = dyadic_fan_controls(counts)
    odd = odd_prime_controls(counts)
    names = ['docs/prime_digit_lifts.md', 'scripts/run_prime_digit_lift_checks.py',
             'scripts/run_adjacent_rank_checks.py', 'scripts/run_oracle_enumeration_checks.py',
             'scripts/run_projection_tag_checks.py', 'scripts/run_signed_face_checks.py']
    report = {'command': 'python3 -B scripts/run_prime_digit_lift_checks.py',
              'counts': dict(sorted(counts.items())), 'exact_cyclotomic_fixtures': exact,
              'cross_characteristic_fixtures': modular, 'dyadic_integral_fixtures': integral,
              'complete_fan_fixtures': fans, 'odd_prime_obstructions': odd,
              'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
              'unrestricted_G3_improved': False, 'polynomial_lower_separation_proved': False,
              'novelty_established': False, 'goal_complete': False}
    (ROOT/'results/research/prime_digit_lift_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'counts': report['counts'],
                      'complete_fans': [{key: row[key] for key in
                       ('rank', 'k', 'initial_ideals_in_every_tested_field')} for row in fans]}, indent=2))


if __name__ == '__main__':
    main()
