#!/usr/bin/env python3
"""Exact independent controls for fixed-equation coherent prime obstructions."""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb
from pathlib import Path
import hashlib
import json

from run_higher_interaction_checks import (
    equation, ground, linear_sum, membership, monic_model, rank_trace,
)
from run_hypersurface_truncation_checks import poly_mul, relation_rows
from run_structural_complexity_checks import dot

ROOT = Path(__file__).resolve().parents[1]
CHANNELS = ((0, 0), (1, 0), (0, 1), (0, 2))
DEGREES = (0, 1, 1, 2)


@lru_cache(None)
def binary_reduce(a, b):
    """Separate homogeneous monic binary reduction over Z[U,V]."""
    if a and b:
        return {(ch, i, j+1): c for (ch, i, j), c in binary_reduce(a-1, b-1).items()}
    if a >= 2:
        terms = [(1, binary_reduce(a-2, b), 1, 0),
                 (-1, binary_reduce(a-2, b+2), 0, 0)]
    elif b >= 3:
        terms = [(1, binary_reduce(a, b-2), 1, 0),
                 (-1, binary_reduce(a+1, b-3), 0, 1)]
    else:
        return {(CHANNELS.index((a, b)), 0, 0): 1}
    out = defaultdict(int)
    for sign, row, di, dj in terms:
        for (ch, i, j), c in row.items():
            out[ch, i+di, j+dj] += sign*c
    return {v: c for v, c in out.items() if c}


def channel_model(s, r):
    basis = tuple((ch, a, b) for ch, d in enumerate(DEGREES)
                  for b in range(max(0, (s-d)//2+1))
                  for a in range(s-d-2*b+1))
    index = {v: k for k, v in enumerate(basis)}

    @lru_cache(None)
    def evaluate(alpha):
        x, y, z = alpha
        if sum(alpha) > s:
            return {}
        out = defaultdict(int)
        for (ch, i, j), c in binary_reduce(x, y).items():
            for ell in range(i+1):
                row = ch, z+2*i+(r-2)*ell, j+ell
                if row in index:
                    out[index[row]] += c*(-1)**i*comb(i, ell)
        return {k: c for k, c in out.items() if c}
    return basis, evaluate


def parent_column(s, r, ch, a, i, j):
    out = {}
    for ell in range(i+1):
        A, B = a+2*i+(r-2)*ell, j+ell
        if DEGREES[ch]+A+2*B <= s:
            out[ch, A, B] = (-1)**i*comb(i, ell)
    return out


def model_controls():
    counts = defaultdict(int)
    for s, r in product(range(3, 10), (2, 3, 4, 5)):
        basis, evaluate = channel_model(s, r)
        original_basis, original = monic_model(s, r)
        assert len(basis) == len(original_basis) == (s+1)**2
        E = ground(s)
        for a in E:
            replay = linear_sum((c, evaluate(original_basis[k]))
                                for k, c in original(a).items())
            assert replay == evaluate(a)
            counts['independent_integral_model_columns'] += 1
        # Test every parent block formula, including negative q, without inverting z.
        for ch, d in enumerate(DEGREES):
            for a in range(s-d+1):
                for i in range((s-d-a)//2+1):
                    for j in range((s-d-a-2*i)//2+1):
                        q = a+2*i-(r-2)*j
                        # Python's // is floor; use the explicit integer ceiling.
                        b0 = 0
                        if r > 2:
                            b0 = max(0, -(q//(r-2)))
                        else:
                            assert q >= 0
                        b1 = (s-d-q)//r
                        assert b0 <= j <= b1
                        J, L = j-b0, b1-b0+1
                        proposed = {(ch, q+(r-2)*(b0+k), b0+k):
                                    (-1)**i*comb(i, k-J)
                                    for k in range(J, min(L, J+i+1))}
                        actual = parent_column(s, r, ch, a, i, j)
                        assert actual == proposed
                        for t in (Fraction(1, 7), Fraction(2, 3)):
                            cost = a+t*(d+2*i+2*j)
                            assert cost == q+(r-2+2*t)*j-2*(1-t)*i+d*t
                            counts['signed_block_price_identities'] += 1
                        counts['signed_block_columns'] += 1
    # Full independent Macaulay relation spaces, rather than just lengths.
    for s, r, p in product((5, 7), (3, 4), (0, 2, 3, 101)):
        basis, evaluate = channel_model(s, r)
        E = ground(s)
        rows = relation_rows(E, [equation(r)], p)
        assert len(E)-len(rows) == len(basis)
        for row in rows:
            v = linear_sum((c, evaluate(a)) for a, c in zip(E, row) if c)
            assert all(not (c % p if p else c) for c in v.values())
            counts['independent_macaulay_relation_rows'] += 1
        # Units z^a V^b times each channel are images of original monomials.
        for k, (ch, a, b) in enumerate(basis):
            x, y = CHANNELS[ch]
            assert evaluate((x+b, y+b, a)) == {k: 1}
        counts['full_macaulay_kernel_equalities'] += 1
    return dict(counts)


def mul(a, b, L):
    return [sum(a[i]*b[k-i] for i in range(k+1)
                if i < len(a) and k-i < len(b)) for k in range(L)]


def series_controls():
    N = 64
    C = [1]
    for n in range(1, N+1):
        C.append(sum(C[i]*C[j]*C[n-1-i-j]
                     for i in range(n) for j in range(n-i)))
    A = mul(C, C, N+1)
    assert all(A[k]*(k+1) == comb(3*k+1, k) for k in range(N+1))
    fixtures, counts = [], defaultdict(int)
    for K in range(1, N+1):
        L = 2*K+2
        # Raw binomial columns and integral triangular elimination do not use C.
        columns = [[comb(j//2, k-j) if j <= k <= j+j//2 else 0
                    for k in range(L)] for j in range(L-1)]
        residual = [(-1)**k for k in range(L)]
        coefficients = []
        for j, column in enumerate(columns):
            assert column[:j] == [0]*j and column[j] == 1
            value = residual[j]
            coefficients.append(value)
            residual = [a-value*b for a, b in zip(residual, column)]
            counts['integer_unit_pivots'] += 1
        assert residual == [0]*(L-1)+[-A[K]]
        assert coefficients == [C[j//2] if j % 2 == 0 else -A[j//2]
                                for j in range(L-1)]
        counts['independent_cyclic_smith_blocks'] += 1
        if K <= 10:
            fixtures.append({'K': K, 'dimension': L, 'last_smith_entry': A[K],
                             'unit_smith_entries': L-1,
                             'integer_elimination_coefficients': coefficients})
    return dict(counts), fixtures, C, A


def prime_controls():
    primes = [p for p in range(2, 10001)
              if all(p % d for d in range(2, int(p**0.5)+1))]
    fixtures = []
    for p in primes:
        K = 1 if p == 2 else 3 if p in (3, 5) else (p+1)//3
        # ceil((p-1)/3) = (p+1)//3.
        if p >= 7:
            assert 2*K+2 <= p <= 3*K+1 < 2*p and K+1 < p
        a = comb(3*K+1, K)//(K+1)
        assert a % p == 0 and a % (p*p) != 0
        s = 6*K+5
        if p >= 7:
            assert s <= 2*p+7
        if p <= 101:
            fixtures.append({'prime': p, 'K': K, 'cutoff': s,
                             'length': (s+1)**2, 'valuation': 1,
                             'coefficient': a})
    return {'prime_divisibility_witnesses_through_10000': len(primes)}, fixtures


def certificate(K, C, A):
    s = 6*K+5
    U = {(2, 0, 0): 1, (0, 2, 0): 1}
    powers, power = [], {(0, 0, 0): 1}
    for i in range(K+2):
        powers.append(power)
        power = poly_mul(power, U, 0)
    out = defaultdict(int, {(0, 0, 2): 1})
    for j in range(2*K+1):
        i, e = 1+j//2, j % 2
        c = C[j//2] if not e else -A[j//2]
        for (x, y, z), v in powers[i].items():
            out[x+j, y+j, z+e] -= c*(-1)**i*v
    beta = (2*K+1, 2*K+1, 2*K+3)
    out[beta] += A[K]
    polynomial = {a: c for a, c in out.items() if c}
    witness = membership(polynomial, s, 3)
    return polynomial, witness


def factorial_valuation(n, p):
    out = 0
    while n:
        n //= p
        out += n
    return out


def finite_prime_controls():
    fixtures, matrix_count = [], 0
    samples = ((2,), (3,), (5,), (2, 3), (2, 3, 5), (2, 3, 5, 7),
               (2, 3, 5, 7, 11), (2, 3, 5, 7, 11, 13, 17, 19),
               (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
                43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101))
    for primes in samples:
        K, modulus = 0, 1
        for p in primes:
            residue = 1 if p == 2 else 3 if p in (3, 5) else (p+1)//3
            q = p*p if p in (2, 3) else p
            K += modulus*((residue-K)*pow(modulus, -1, q) % q)
            modulus *= q
        assert 0 < K < modulus
        valuations = {}
        for p in primes:
            v, denominator = 0, K+1
            while denominator % p == 0:
                v += 1
                denominator //= p
            valuations[str(p)] = (factorial_valuation(3*K+1, p)
                                  -factorial_valuation(K, p)
                                  -factorial_valuation(2*K+1, p)-v)
            assert valuations[str(p)] >= 1
        record = {'primes': list(primes), 'K': K, 'CRT_modulus': modulus,
                  'prime_valuations_of_last_smith_entry': valuations}
        if K <= 1000:
            coefficient = comb(3*K+1, K)//(K+1)
            assert all(coefficient % p == 0 for p in primes)
            record['last_smith_entry'] = coefficient
        if len(primes) >= 2 and K <= 100:
            # Independent integer elimination of the actual block, beyond
            # the one-prime arithmetic test; no giant monomial enumeration.
            L = 2*K+2
            residual = [(-1)**j for j in range(L)]
            for j in range(L-1):
                c = residual[j]
                for k in range(j, min(L, j+j//2+1)):
                    residual[k] -= c*comb(j//2, k-j)
            assert residual == [0]*(L-1)+[-record['last_smith_entry']]
            assert all(residual[-1] % p == 0 for p in primes)
            matrix_count += 1
        fixtures.append(record)
    return {'simultaneous_finite_prime_sets': len(fixtures),
            'simultaneous_integer_smith_blocks': matrix_count}, fixtures


def terms(f):
    return [{'exponent': list(a), 'coefficient': c} for a, c in sorted(f.items())]


def original_controls(C, A):
    fixtures, counts = [], defaultdict(int)
    for K in range(1, 7):
        s, H = 6*K+5, 6*K+1
        B, M = s+1, (s+1)**3
        w, perturbed = (1, 1, H), (M+1, M+B, M*H+B*B)
        _, evaluate = monic_model(s, 3)
        E = ground(s)
        digits = {a: a[0]+B*a[1]+B*B*a[2] for a in E}
        assert len(set(digits.values())) == len(E)
        assert min(digits.values()) == 0 and max(digits.values()) < M
        threshold = 2*H
        strict = [a for a in E if dot(a, w) < threshold]
        inclusive_other = [a for a in E if dot(a, w) <= threshold and a != (0, 0, 2)]
        generic = [a for a in E if dot(a, perturbed) < 2*perturbed[2]]
        assert set(generic) == set(inclusive_other)
        assert all(dot(a, perturbed) == M*dot(a, w)+digits[a] for a in E)
        counts['generic_weight_ground_comparisons'] += len(E)
        assert all(a[2] == 1 and a[0]+a[1] == H
                   for a in set(inclusive_other)-set(strict))
        polynomial, witness = certificate(K, C, A)
        beta = (2*K+1, 2*K+1, 2*K+3)
        assert sum(beta) == s and polynomial[beta] == A[K]
        assert all(a in strict for a in polynomial if a not in (beta, (0, 0, 2)))
        counts['integral_expanded_membership_certificates'] += 1
        fields = [0]+[p for p in (2, 3, 5, 7, 11, 13, 17, 19)
                      if A[K] % p == 0]+[101]
        ranks = []
        for label, prefix, weight in (
            ('strict_equal_binary', strict, w),
            ('all_other_ties_then_generic', generic, perturbed),
        ):
            ordered = sorted(prefix, key=lambda a: (dot(a, weight), a))+[(0, 0, 2)]
            for p in fields:
                _, values = rank_trace(ordered, evaluate, p, {len(ordered)-1, len(ordered)})
                before, after = values[len(ordered)-1], values[len(ordered)]
                assert after-before == int(not p or A[K] % p != 0)
                ranks.append({'prefix': label, 'characteristic': p, 'rank_before': before,
                              'rank_after_z_square': after, 'increment': after-before})
                counts['independent_original_monomial_rank_tests'] += 1
        fixtures.append({'K': K, 'cutoff': s, 'length': B*B, 'coefficient': A[K],
                         'equal_binary_weight': list(w), 'generic_weight': list(perturbed),
                         'strict_parent_columns': len(strict),
                         'generic_parent_columns': len(generic), 'field_rank_tests': ranks,
                         'integer_relation': terms(polynomial),
                         'integer_monic_division_multiplier': terms(witness)})
        print('Verified original monomial fixture K='+str(K), flush=True)
    return dict(counts), fixtures


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    model_counts = model_controls()
    series_counts, smith, C, A = series_controls()
    prime_counts, primes = prime_controls()
    simultaneous_counts, simultaneous = finite_prime_controls()
    original_counts, original = original_controls(C, A)
    sources = ('docs/fixed_equation_bad_primes.md', 'scripts/run_fixed_equation_prime_checks.py',
               'docs/higher_interaction_fans.md', 'scripts/run_higher_interaction_checks.py',
               'scripts/run_hypersurface_truncation_checks.py',
               'scripts/run_socle_stability_checks.py', 'scripts/run_structural_complexity_checks.py',
               'scripts/run_interleaving_research.py')
    report = {
        'status': 'passed',
        'command': 'python3 -B scripts/run_fixed_equation_prime_checks.py',
        'counts': {**model_counts, **series_counts, **prime_counts,
                   **simultaneous_counts, **original_counts},
        'cyclic_matrix_fixtures': smith, 'prime_witnesses': primes,
        'simultaneous_prime_witnesses': simultaneous,
        'original_monomial_fixtures': original,
        'uniform_result': {
            'fixed_equation': 'x^2+y^2+z^2+x*y*z^3',
            'cutoff': 's=6*K+5', 'length': '(s+1)^2',
            'block_smith_form': 'diag(1 repeated 2*K+1 times, binom(3*K+1,K)/(K+1))',
            'distinguished_monomial': [0, 0, 2],
            'rank_increment': '1 in characteristic zero or p not dividing the last Smith entry; 0 otherwise',
            'every_prime_occurs': True,
            'every_fixed_finite_prime_set_can_fail_simultaneously': True,
            'simultaneous_prime_budget': 'K<product(p^e_p), with e_2=e_3=2 and other e_p=1, by explicit CRT congruences',
            'prime_cutoff_budget': 's<=2*p+7 and n<=(2*p+8)^2 for p>=7; explicit small-prime witnesses',
            'ties_removed_by_positive_integer_weights': True,
            'whole_integer_algebra_is_free': True},
        'scope': 'Finite controls support the uniform proof. Torsion is in a coherent prefix block quotient, not the whole algebra or a computed local cohomology module. The equation is fixed while the cutoff and weights vary. No additional fan cones or unrestricted G_3 improvement follow.',
        'novelty_established': False, 'unrestricted_G3_improvement': False,
        'goal_complete': False,
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
    (ROOT/'results/research/fixed_equation_prime_checks.json').write_text(
        json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
