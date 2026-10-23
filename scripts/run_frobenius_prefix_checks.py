#!/usr/bin/env python3
"""Complete original-prefix Frobenius ranks and explicit support matchings."""
from collections import Counter
from fractions import Fraction as Rat
from itertools import product
from math import comb, factorial
from pathlib import Path
from random import Random
import hashlib
import json

from run_generic_hook_checks import matching_basis, bareiss
from run_hook_bridge_checks import model, actual
from run_hook_rank_departure_checks import complete_fan
from run_structural_complexity_checks import area2, dot

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


class Field:
    def __init__(self, p, extension=False):
        self.p, self.extension = p, extension

    def add(self, a, b):
        return a ^ b if self.extension else (a+b) % self.p if self.p else a+b

    def neg(self, a):
        return a if self.extension else (-a) % self.p if self.p else -a

    def mul(self, a, b):
        if not self.extension:
            return a*b % self.p if self.p else a*b
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a & 4:
                a ^= 7                 # alpha^2 + alpha + 1
        return out

    def inv(self, a):
        assert a
        return (0, 1, 3, 2)[a] if self.extension else pow(a, -1, self.p) if self.p else Rat(1, a)

    def power(self, a, k):
        out = 1
        while k:
            if k & 1:
                out = self.mul(out, a)
            a = self.mul(a, a)
            k >>= 1
        return out


def accumulate(out, v, c, F):
    if not c:
        return
    for i, b in v.items():
        a = F.add(out.get(i, 0), F.mul(c, b))
        if a:
            out[i] = a
        else:
            out.pop(i, None)


def rank(vectors, F):
    if F.p == 2 and not F.extension:
        piv = {}
        for col in vectors:
            v = sum(1 << i for i, c in col.items() if c)
            while v:
                i = v.bit_length()-1
                if i not in piv:
                    piv[i] = v
                    break
                v ^= piv[i]
        return len(piv)
    piv = {}
    for col in vectors:
        v = dict(col)
        while v:
            i = max(v)
            if i not in piv:
                c = F.inv(v[i])
                piv[i] = {j: F.mul(c, a) for j, a in v.items()}
                break
            accumulate(v, piv[i], F.neg(v[i]), F)
    return len(piv)


def multiply(a, b, A, F):
    out = {}
    for i, c in a.items():
        for j, d in b.items():
            if i+j < A:
                accumulate(out, {i+j: d}, c, F)
    return out


def quotient(Q, H, G, F):
    """Direct monomial reduction, independently of the saved hook model."""
    d, t = Q//2, Q+1
    A = t+2*d
    a = A+t
    D = tuple((i, j) for i, j in product(range(a), repeat=2) if i<t or j<t)
    index = {v: i for i, v in enumerate(D)}
    offset = len(D)
    tails = [(t, 0, 1), (0, t, 1)]
    tails += [(t+i, 1, c) for i, c in H.items()]
    tails += [(1, t+i, c) for i, c in G.items()]
    C = {(0, *v): {index[v]: 1} for v in D}
    C.update({(1, *v): {offset+index[v]: 1} for v in D})
    for i, j in product(range(A), repeat=2):
        v = {}
        for y, z, c in tails:
            if (i+y, j+z) in index:
                accumulate(v, {offset+index[i+y, j+z]: c}, 1, F)
        if v:
            C[0, t+i, t+j] = v
    assert 2*len(D) == 4*A*t+2*t*t
    return t, d, A, a, D, index, C


def fixture(Q, p, extension=False, sparse=False):
    remainder = Q
    while remainder % p == 0:
        remainder //= p
    assert Q >= 2 and remainder == 1
    F = Field(p, extension)
    t, d = Q+1, Q//2
    A = t+2*d
    degree = A//2
    rng = Random(20301+Q+100*p+extension)
    size = 4 if extension else p
    H = {i: rng.randrange(1, size) for i in range(degree+1)}
    G = {i: rng.randrange(1, size) for i in range(degree+1)}
    if sparse:
        H, G, degree = {0: 1, 1: 1}, {}, 1
    t, d, A, a, D, index, C = quotient(Q, H, G, F)
    n, r = 2*len(D), t+d-1
    w0 = ((Q+t)*A, 1, A)
    M = 4*a+1
    K = 300*a*M*M
    w = tuple(K*x+y for x, y in zip(w0, (1, M, M*M)))
    cutoff = K*(Rat(w0[0]+t+r)+Rat(1, 3))
    assert len({dot(w, e) for e in C}) == len(C)
    prefix = {e for e in C if dot(w, e)<cutoff}
    assert prefix == {e for e in C if dot(w0, e)<cutoff/K}
    U = {(1, t+i, 0) for i in range(r+1)}
    L = {(0, t+i, t+j) for i, j in product(range(t, A), range(Q))}
    L |= {(0, t+i, t+Q) for i in range(t, r+1)}
    heads = {(0, t+i, t+j) for i, j in product(range(t), repeat=2)}
    residual = U | L | heads
    assert len(residual) == t*(A+1)
    fixed = {(0, *v) for v in D} | {(1, i, j) for i, j in product(range(t), repeat=2)}
    assert prefix-fixed == residual
    assert not residual & fixed
    COUNTS['complete_original_prefix_membership_checks'] += 1

    # An explicit perfect support matching on all original prefix columns.
    matching = {}
    for e in prefix & fixed:
        matching[e] = next(iter(C[e]))
    for e in U:
        matching[e] = len(D)+index[e[1], e[2]]
    for e in L:
        matching[e] = len(D)+index[e[1], e[2]-t]
    for i in range(t, r+1):
        for j in range(Q):
            q, qq = i+min(j, d), i+min(j+1, d)
            matching[0, t+q, t+j] = len(D)+index[t+qq, j+1]
    for i, j in product(range(t), repeat=2):
        matching[0, t+i, t+j] = len(D)+index[i, t+j]
    assert set(matching) == prefix
    assert len(set(matching.values())) == len(prefix)
    assert all(row in C[e] for e, row in matching.items())
    COUNTS['explicit_perfect_matching_edges'] += len(matching)
    if Q <= 9:
        B, _ = matching_basis(sorted(prefix), C, len(prefix))
        assert set(B) == prefix
        COUNTS['independent_maximum_matching_checks'] += 1

    actual_rank = rank([C[e] for e in sorted(prefix)], F)
    assert actual_rank == len(prefix)-d
    COUNTS['independent_full_original_prefix_ranks'] += 1
    powers = [{0: 1}]
    for j in range(Q):
        powers.append(multiply(powers[-1], H, A, F))
    assert powers[Q] == {i*Q: F.power(c, Q) for i, c in H.items() if i*Q<A}
    for i in range(t, r+1):
        relation = dict(C[1, t+i, 0])
        for j in range(Q+1):
            sign = F.neg(1) if j % 2 == 0 else 1
            for k, c in powers[j].items():
                if i+k >= A:
                    continue
                e = (0, t+i+k, t+j)
                assert e in L
                accumulate(relation, C[e], F.mul(sign, c), F)
        assert not relation
        assert (1, t+i, 0) in prefix
        COUNTS['independent_complete_prefix_kernel_relations'] += 1

    B, _ = matching_basis(sorted(C, key=lambda e: dot(w, e)), C, n)
    assert prefix <= B
    B_rank = rank([C[e] for e in B], F)
    assert B_rank <= n-d
    COUNTS['unique_whole_support_minimum_basis_rank_checks'] += 1
    independent_fan = None
    if not extension and Q <= 3:
        data = model(t, a)
        terms = {(t, 0): 1, (0, t): 1}
        terms.update({(t+i, 1): c for i, c in H.items()})
        terms.update({(1, t+i): c for i, c in G.items()})
        co = [terms.get((y, z), 0) if g == 1 else 0 for g, y, z in data[2]]
        old = actual(data, co, p, True)
        for e, v in old.items():
            assert C.get(e, {}) == {len(D)+i: c for i, c in v.items()}
        COUNTS['independent_original_constructor_column_equalities'] += len(old)
        if Q == 2:
            fan = complete_fan(old, n, p)
            assert sum(area2(poly) for poly in fan.values()) == 1
            independent_fan = {'actual_cones': len(fan), 'exact_total_area': '1'}
            COUNTS['independent_complete_original_fans'] += 1
    return {'Q': Q, 'characteristic': p, 'field_size': 4 if extension else p,
            't': t, 'A': A, 'degree': degree, 'length': n, 'sparse_zero_other_arm': sparse,
            'H': H, 'G': G, 'strict_weight': w, 'cutoff': str(cutoff),
            'residual_prefix_columns': len(residual), 'full_prefix_columns': len(prefix),
            'actual_full_prefix_rank': actual_rank, 'support_full_prefix_rank': len(prefix),
            'exact_prefix_nullity': d, 'whole_support_basis_actual_rank': B_rank,
            'whole_support_basis_deficit': n-B_rank,
            'independent_complete_original_fan': independent_fan}


def characteristic_zero_controls():
    out = []
    for Q in (2, 3, 4, 5, 8, 9, 16):
        d, t = Q//2, Q+1
        A, D = t+2*d, (t+2*d)//2
        L = [[comb(D*Q, d+i-j) if 0<=d+i-j<=D*Q else 0
              for j in range(d)] for i in range(d)]
        det = bareiss(L)
        predicted = Rat(1)
        for j in range(d):
            predicted *= Rat(factorial(D*Q+j)*factorial(j),
                             factorial(d+j)*factorial(D*Q-d+j))
        assert det == predicted and det > 0
        COUNTS['independent_characteristic_zero_solid_determinants'] += 1
        if Q <= 4:
            F = Field(0)
            H = {i: comb(D, i) for i in range(D+1)}
            _, _, _, _, base, _, C = quotient(Q, H, H, F)
            r = t+d-1
            w0 = ((Q+t)*A, 1, A)
            tau = Rat(w0[0]+t+r)+Rat(1, 3)
            P = [C[e] for e in C if dot(w0, e)<tau]
            assert rank(P, F) == len(P)
            COUNTS['independent_characteristic_zero_full_prefix_ranks'] += 1
        out.append({'Q': Q, 'degree': D, 'order': d,
                    'solid_determinant': str(det), 'nonzero': True})
    return out


def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    previous_path = 'results/research/transverse_lattice_rank_checks.json'
    previous = json.loads((ROOT/previous_path).read_text())
    assert previous['status'] == 'passed'
    for p, h in previous['source_sha256'].items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == h
    F = Field(2, True)
    for a, b, c in product(range(4), repeat=3):
        assert F.mul(a, F.add(b, c)) == F.add(F.mul(a, b), F.mul(a, c))
        assert F.mul(F.mul(a, b), c) == F.mul(a, F.mul(b, c))
    assert all(F.mul(a, F.inv(a)) == 1 for a in range(1, 4))
    cases = []
    for Q, p, ext in [(2, 2, False), (4, 2, False), (8, 2, False),
                      (16, 2, False), (3, 3, False), (9, 3, False),
                      (5, 5, False), (7, 7, False), (4, 2, True)]:
        case = fixture(Q, p, ext)
        cases.append(case)
        print('Frobenius prefix passed', Q, p, case['field_size'],
              'length', case['length'], 'nullity', case['exact_prefix_nullity'], flush=True)
    cases.append(fixture(4, 2, sparse=True))
    zero = characteristic_zero_controls()
    files = set(previous['source_sha256']) | {previous_path,
             'docs/frobenius_coherent_prefix_gaps.md', 'scripts/run_frobenius_prefix_checks.py'}
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_frobenius_prefix_checks.py',
              'uniform_result': {
                  'parameters': 'Q=p^e>=2, d=floor(Q/2), t=Q+1, A=t+2d, a=A+t',
                  'class': 'Zero outer maps; W_Y=1+vH(q), W_Z=1+vG(q), H_0 H_1 nonzero; arbitrary remaining coefficients over every field of characteristic p.',
                  'complete_residual_prefix_columns': 'm=t(A+1)',
                  'complete_prefix_support_rank': 'm, plus its compulsory coloops',
                  'complete_prefix_actual_rank': 'm-d, plus the same compulsory coloops',
                  'whole_unique_support_minimum_basis': 'Actual rank at most n-d, with n=4At+2t^2=Theta(Q^2).',
                  'scope': 'Uniform growing coherent rank deficit, not a fan-count lower rate or an unrestricted G_3 theorem.'},
              'counts': dict(COUNTS), 'cases': cases,
              'characteristic_zero_controls': zero,
              'limits': ['Frobenius and the binomial Toeplitz determinant are classical.',
                         'The stated coefficient stratum is smaller than the full generic compatible hook model.',
                         'No actual fan growth rate follows from the support basis rank deficit.',
                         'The existing O(n^(3/2)*a) box ceiling is O(n^2) on this parameter sequence, so this family cannot supply a polynomial fan excess.',
                         'Exact combined priority and the original unrestricted endpoint remain open.'],
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(files)}}
    (ROOT/'results/research/frobenius_prefix_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print('Frobenius complete-prefix checks passed', json.dumps(dict(COUNTS), sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
