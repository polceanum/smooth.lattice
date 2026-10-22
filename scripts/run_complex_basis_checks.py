#!/usr/bin/env python3
"""Exact guards for complex transforms of actual monomial residue matroids.

Classical constructions and finite checks; no new extremal count is certified.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import log, pi, prod
from pathlib import Path
import cmath
import hashlib
import json
import random

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / 'results/research/complex_basis_checks.json'
COUNTS = Counter()


def rref(rows, width, p=None):
    rows = [[Q(x) if p is None else int(x) % p for x in row] for row in rows]
    pivots = []
    for j in range(width):
        k = next((k for k in range(len(pivots), len(rows)) if rows[k][j]), None)
        if k is None:
            continue
        i = len(pivots)
        rows[i], rows[k] = rows[k], rows[i]
        t = 1 / rows[i][j] if p is None else pow(rows[i][j], -1, p)
        rows[i] = [x * t if p is None else x * t % p for x in rows[i]]
        for k in range(len(rows)):
            if k != i and rows[k][j]:
                t = rows[k][j]
                rows[k] = [x - t * y if p is None else (x - t * y) % p
                           for x, y in zip(rows[k], rows[i])]
        pivots.append(j)
    return rows[:len(pivots)], pivots


def monomials(cap):
    return sorted((a for a in product(range(cap), repeat=3) if sum(a) < cap),
                  key=lambda a: (sum(a), a))


def quotient(cap, relations, p=None):
    """Generate the ideal inside K[x,y,z]/m^cap, then row-reduce it."""
    D = monomials(cap)
    index = {a: i for i, a in enumerate(D)}
    rows = []
    for f in relations:
        for shift in D:
            row = [0] * len(D)
            for a, c in f.items():
                target = tuple(x + y for x, y in zip(a, shift))
                if target in index:
                    row[index[target]] += c
            if any(row):
                rows.append(row)
    reduced, pivots = rref(rows, len(D), p)
    free = [i for i in range(len(D)) if i not in pivots]

    def residue(a):
        if a not in index:
            return tuple(0 for _ in free)
        v = [0] * len(D)
        v[index[a]] = 1
        for row, j in zip(reduced, pivots):
            t = v[j]
            v = [x - t * y if p is None else (x - t * y) % p
                 for x, y in zip(v, row)]
        return tuple(v[j] for j in free)

    mats = [tuple(residue(tuple(a[k] + (k == j) for k in range(3)))
                  for a in (D[i] for i in free)) for j in range(3)]

    def apply(M, v):
        return tuple(sum(M[j][i] * v[j] for j in range(len(free)))
                     if p is None else
                     sum(M[j][i] * v[j] for j in range(len(free))) % p
                     for i in range(len(free)))

    for i, j in combinations(range(3), 2):
        for k in range(len(free)):
            e = tuple(int(k == l) for l in range(len(free)))
            assert apply(mats[i], apply(mats[j], e)) == apply(mats[j], apply(mats[i], e))
            COUNTS['independent_commuting_column_checks'] += 1
    for a in D:
        v = residue((0, 0, 0))
        for j in range(3):
            for _ in range(a[j]):
                v = apply(mats[j], v)
        assert v == residue(a)
        COUNTS['original_monomial_column_replays'] += 1
    for f in relations:
        v = [sum(c * residue(a)[i] for a, c in f.items()) for i in range(len(free))]
        if p is None:
            assert not any(v)
        else:
            assert not any(x % p for x in v)
        COUNTS['original_ideal_relation_checks'] += 1
    return D, free, residue


def rank(columns, p=None):
    if not columns:
        return 0
    return len(rref(columns, len(columns[0]), p)[1])


def det(A):
    n = len(A)
    result = Q(0)
    for perm in permutations(range(n)):
        sign = (-1) ** sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
        result += sign * prod(A[i][perm[i]] for i in range(n))
    return result


def padd(A, B, scale=1):
    C = Counter(A)
    for a, c in B.items():
        C[a] += scale * c
    return {a: c for a, c in C.items() if c}


def pmul(A, B):
    C = Counter()
    for a, x in A.items():
        for b, y in B.items():
            C[tuple(i + j for i, j in zip(a, b))] += x * y
    return {a: c for a, c in C.items() if c}


def gram_polynomial(V, labels):
    n, m = len(V), len(labels)
    entries = [[{} for _ in range(n)] for _ in range(n)]
    for i, j in product(range(n), repeat=2):
        for k in range(m):
            entries[i][j] = padd(entries[i][j], {labels[k]: V[i][k] * V[j][k]})
    result = {}
    for perm in permutations(range(n)):
        term = {(0, 0, 0): 1}
        for i in range(n):
            term = pmul(term, entries[i][perm[i]])
        sign = (-1) ** sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
        result = padd(result, term, sign)
    return result


def minor_polynomial(V, labels):
    n = len(V)
    result = {}
    for B in combinations(range(len(labels)), n):
        d = det([[row[j] for j in B] for row in V])
        if d:
            a = tuple(sum(labels[j][k] for j in B) for k in range(3))
            result = padd(result, {a: d * d})
        COUNTS['independent_minor_expansions'] += 1
    return result


@dataclass(frozen=True)
class G:
    real: Q = Q(0)
    imag: Q = Q(0)

    def __add__(self, other):
        other = other if isinstance(other, G) else G(Q(other))
        return G(self.real + other.real, self.imag + other.imag)

    __radd__ = __add__

    def __neg__(self):
        return G(-self.real, -self.imag)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        other = other if isinstance(other, G) else G(Q(other))
        return G(self.real * other.real - self.imag * other.imag,
                 self.real * other.imag + self.imag * other.real)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, G) else G(Q(other))
        den = other.real * other.real + other.imag * other.imag
        assert den
        return self * G(other.real / den, -other.imag / den)

    def serial(self):
        return {'real': str(self.real), 'imag': str(self.imag)}


def peval(P, values):
    result = G()
    for a, c in P.items():
        term = G(Q(c))
        for z, power in zip(values, a):
            for _ in range(power):
                term = term * z
        result += term
    return result


def finite_fan(bases):
    """Intersect every minimum-profile cone with the rational positive triangle."""
    profiles = [tuple(sum(a[k] for a in B) for k in range(3)) for B in bases]
    result = []
    for B, mu in zip(bases, profiles):
        polygon = [(Q(0), Q(0)), (Q(1), Q(0)), (Q(0), Q(1))]
        for nu in profiles:
            d = tuple(x - y for x, y in zip(nu, mu))
            def value(v):
                return d[0] * v[0] + d[1] * v[1] + d[2] * (1 - v[0] - v[1])
            new = []
            for a, b in zip(polygon, polygon[1:] + polygon[:1]):
                va, vb = value(a), value(b)
                if va >= 0:
                    new.append(a)
                if (va < 0 < vb) or (vb < 0 < va):
                    t = va / (va - vb)
                    new.append(tuple(x + t * (y - x) for x, y in zip(a, b)))
            polygon = new
            if not polygon:
                break
        area = abs(sum(a[0] * b[1] - a[1] * b[0]
                       for a, b in zip(polygon, polygon[1:] + polygon[:1]))) if polygon else 0
        if area:
            assert profiles.count(mu) == 1
            result.append({'basis': [list(a) for a in B], 'twice_area': str(area)})
    assert sum(Q(row['twice_area']) for row in result) == 1
    COUNTS['complete_actual_positive_fans'] += 1
    return result


def length_three():
    records = []
    for c in (Q(1, 4), Q(1, 2), Q(1), Q(2), Q(4)):
        D, free, residue = quotient(3, [{(0, 1, 0): 1, (2, 0, 0): -c}, {(0, 0, 1): 1}])
        assert len(free) == 3
        nonzero = [a for a in D if any(residue(a))]
        bases = [B for B in combinations(nonzero, 3) if rank([residue(a) for a in B]) == 3]
        assert set(map(frozenset, bases)) == {
            frozenset(((0, 0, 0), (1, 0, 0), (2, 0, 0))),
            frozenset(((0, 0, 0), (1, 0, 0), (0, 1, 0)))}
        assert len(finite_fan(bases)) == 2
        V = [[residue(a)[i] for a in nonzero] for i in range(3)]
        P = gram_polynomial(V, nonzero)
        # The quotient's row-reduction basis changes the Gram determinant by
        # one global positive factor; compare normalized coefficient ratios.
        assert set(P) == {(3, 0, 0), (1, 1, 0)}
        factor = P[(3, 0, 0)]
        assert P[(1, 1, 0)] / factor == c * c
        assert P == minor_polynomial(V, nonzero)
        COUNTS['actual_quotient_gram_identities'] += 1
        # Exact algebraic root after the substitution t=exp(-s).
        assert 1 + c * c * (-1 / (c * c)) == 0
        errors = []
        for j in range(-2, 3):
            s = complex(2 * log(float(c)), (2 * j + 1) * pi)
            errors.append(abs(1 + float(c * c) * cmath.exp(-s)))
        assert max(errors) < 1e-12
        records.append({'c': str(c), 'length': 3, 'positive_initials': 2,
                        'zero_real_part': 'log(' + str(c * c) + ')',
                        'numerical_identity_crosscheck_max_error': max(errors)})
    P = {(3, 0, 0): 1, (1, 1, 0): 1}
    for values in [(G(-1, 1), G(0, 2), G(0, 1)),
                   (G(1, 2), G(3, -4), G(1))]:
        assert peval(P, values) == G()
        assert all(z.imag > 0 for z in values) or all(z.real > 0 for z in values)
        COUNTS['exact_compressed_half_plane_zeros'] += 1
    H = [[6, 0, 0], [0, 0, 1], [0, 1, 0]]
    vectors = [(1, 0, 0), (0, 1, 1), (0, 1, -1)]
    gram = [[sum(a[i] * H[i][j] * b[j] for i, j in product(range(3), repeat=2))
             for b in vectors] for a in vectors]
    assert gram == [[6, 0, 0], [0, 2, 0], [0, 0, -2]]
    COUNTS['exact_non_lorentzian_inertia_witnesses'] += 1
    return records


def fano():
    labels = [(3, 0, 0), (0, 3, 0), (0, 0, 3), (2, 1, 0),
              (1, 2, 0), (2, 0, 1), (1, 0, 2)]
    relations = [{a: 1} for a in [(0, 2, 1), (0, 1, 2), (1, 1, 1)]]
    for target, sources in [(3, (0, 1)), (4, (0, 2)), (5, (1, 2)), (6, (0, 1, 2))]:
        relations.append({labels[i]: 1 for i in (target,) + sources})
    D, free, residue = quotient(4, relations, 2)
    assert len(D) == 20 and len(free) == 13
    lower = [a for a in D if sum(a) < 3]
    assert len(lower) == 10 and rank([residue(a) for a in lower], 2) == 10
    bits = [1, 2, 4, 3, 5, 6, 7]
    bases, nonbases = [], []
    full_bases = []
    for S in combinations(range(7), 3):
        available = rank([residue(a) for a in lower + [labels[i] for i in S]], 2) == 13
        assert available == bool(bits[S[0]] ^ bits[S[1]] ^ bits[S[2]])
        (bases if available else nonbases).append(S)
        if available:
            full_bases.append(tuple(lower + [labels[i] for i in S]))
        COUNTS['independent_original_fano_basis_tests'] += 1
    assert len(bases) == 28 and len(nonbases) == 7
    for a in product(range(13), repeat=3):
        if prod(i + 1 for i in a) <= 13:
            assert any(residue(a)) == (a in lower or a in labels)
            COUNTS['whole_universal_ground_loop_checks'] += 1
    fan = finite_fan(full_bases)
    parts = [[], [], [], []]
    for S in bases:
        parts[2 * (0 in S) + (1 in S)].append(tuple(i for i in S if i not in (0, 1)))
    values = {i: z for i, z in zip(range(2, 7), (-3, -3, -3, 1, 2))}
    D0, C, B, A = [sum(prod(values[i] for i in S) for S in group) for group in parts]
    assert (A, B, C, D0) == (-3, -9, 11, 9)
    assert B * C - A * D0 == -72
    COUNTS['negative_rayleigh_witnesses'] += 1
    z = {i: G(Q(values[i]), Q(1, 10)) for i in range(2, 7)}
    D0, C, B, A = [sum((prod((z[i] for i in S), start=G(1)) for S in group), start=G())
                   for group in parts]
    z[0] = G(0, 1)
    z[1] = -(B * z[0] + D0) / (A * z[0] + C)
    assert all(z[i].imag > 0 for i in range(7))
    value = sum((prod((z[i] for i in S), start=G(1)) for S in bases), start=G())
    assert value == G()
    COUNTS['exact_fano_all_upper_half_plane_zeros'] += 1
    return {'characteristic': 2, 'length': 13, 'coloops': lower,
            'fano_labels': labels, 'fano_bases': bases, 'fano_dependent_triples': nonbases,
            'positive_initials': len(fan), 'complete_positive_fan': fan,
            'rayleigh_difference': -72, 'all_upper_half_plane_zero': [z[i].serial() for i in range(7)]}


def main():
    rng = random.Random(20261022)
    for n in (2, 3, 4):
        for _ in range(3):
            labels = [tuple(rng.randrange(3) for _ in range(3)) for _ in range(n + 2)]
            V = [[rng.randrange(-2, 3) for _ in labels] for _ in range(n)]
            assert gram_polynomial(V, labels) == minor_polynomial(V, labels)
            COUNTS['independent_polynomial_cauchy_binet_identities'] += 1
    length_three_records = length_three()
    fano_record = fano()
    sources = ['docs/complex_basis_transforms.md', 'scripts/run_complex_basis_checks.py',
               'docs/cross_layer_direction_budget.md', 'docs/disjoint_socle_circuits.md',
               'docs/initial_cone_stability.md']
    report = {
        'schema_version': 1, 'status': 'passed',
        'command': 'python3 -B scripts/run_complex_basis_checks.py',
        'scope': 'Exact finite actual-quotient and polynomial guards for a classical complex-transform route. No new asymptotic count, strongest-current novelty claim, or goal completion.',
        'counts': dict(sorted(COUNTS.items())), 'length_three': length_three_records,
        'fano_quotient': fano_record,
        'uniform_proof_boundary': {
            'determinant_stability': 'Classical in separate variables over complex representations; exact state support.',
            'monomial_substitution': 'Actual length-three counterexample, including ordinary homogenized Lorentzian failure.',
            'zero_strip': 'The uniform arbitrary-real-part statement follows from the explicit entire function, not numerical sampling.',
            'all_fields': 'Branden Theorem 6.6 applied to the explicitly realized Fano minor rules out every stable polynomial with full basis support.',
            'projected_fan': 'No claim that the fan determines the Fano matroid or that all complex methods are ruled out.'},
        'primary_sources': [
            {'url': 'https://arxiv.org/pdf/math/0202034', 'inspected': 'Theorem 8.1 and Corollary 8.2'},
            {'url': 'https://arxiv.org/pdf/math/0605678', 'inspected': 'Theorem 5.6 and Theorem 6.6'},
            {'url': 'https://arxiv.org/pdf/1902.03719', 'inspected': 'Theorems 2.25 and 3.10'}],
        'source_sha256': {s: hashlib.sha256((ROOT / s).read_bytes()).hexdigest() for s in sources}}
    REPORT.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'counts': report['counts'],
                      'fano_positive_initials': fano_record['positive_initials'],
                      'source_hashes': len(sources)}))


if __name__ == '__main__':
    main()
