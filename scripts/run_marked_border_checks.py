#!/usr/bin/env python3
"""Exact finite controls for marked corner cuts and generic dense output.

Run from the repository root:
    python3 -B scripts/run_marked_border_checks.py

Finite enumeration checks mechanisms in the uniform proof. It cannot prove
an asymptotic estimate or novelty, or settle the unrestricted G_3 problem.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial, gcd, lcm
from pathlib import Path
import hashlib
import json
import sys

from run_interleaving_research import minimal_excluded_border
from run_structural_complexity_checks import area2, clip, dot
from run_signed_face_checks import field_basis, field_reduce

ROOT = Path(__file__).resolve().parents[1]


def border(D):
    """Independent antichain extraction from the entire enclosing box."""
    D = frozenset(D)
    d = len(next(iter(D)))
    candidates = sorted(set(product(*(range(max(a[i] for a in D)+2)
                                      for i in range(d))))-D,
                        key=lambda a: (sum(a), a))
    heads = []
    for a in candidates:
        if not any(all(b[i] <= a[i] for i in range(d)) for b in heads):
            heads.append(a)
    assert sorted(heads) == minimal_excluded_border(tuple(sorted(D)))
    return tuple(sorted(heads))


def integral(w):
    den = lcm(*(x.denominator for x in w))
    v = tuple(int(x*den) for x in w)
    g = gcd(*v)
    return tuple(x//g for x in v)


def separator(D, E):
    d = len(next(iter(D)))
    inequalities = {tuple(a[i]-b[i] for i in range(d)) for a in E for b in D}
    if d == 2:
        lo, hi = Q(0), None
        for x, y in inequalities:
            if not x:
                if y <= 0:
                    return None
            elif x > 0:
                lo = max(lo, Q(-y, x))
            else:
                upper = Q(-y, x)
                hi = upper if hi is None else min(hi, upper)
        if hi is not None and hi <= lo:
            return None
        t = lo+1 if hi is None else (lo+hi)/2
        w = integral((t, Q(1)))
        region = (lo, hi)
    else:
        assert d == 3
        poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for normal in sorted(inequalities):
            if min(normal) >= 0:
                continue
            if max(normal) <= 0:
                return None
            poly = clip(poly, normal, 1)
            if not poly:
                return None
        assert area2(poly) > 0
        w = integral(tuple(sum(a[i] for a in poly) for i in range(3)))
        region = poly
    assert min(w) > 0 and all(dot(q, w) > 0 for q in inequalities)
    return w, region


def chart_checks(D, E, w, counts):
    n, d = len(D), len(w)
    if n == 1:
        assert len(E) == d
        return
    t = max(dot(w, a) for a in D)
    J = tuple(i for i in range(d) if w[i] <= t)
    r = len(J)
    L = []
    for i in J:
        h = 1
        while 2*h*w[i] <= t:
            h *= 2
        assert h*w[i] <= t < 2*h*w[i]
        L.append(h)
    volume = 1
    for h in L:
        volume *= h
    assert n < 4**r*volume and volume <= factorial(r)*n
    assert all(all(not a[i] for i in range(d) if i not in J) for a in D)
    assert all(all(a[i] < 2*h for i, h in zip(J, L)) for a in D)
    active_heads = tuple(a for a in E if all(not a[i] for i in range(d) if i not in J))
    inactive_heads = set(E)-set(active_heads)
    assert inactive_heads == {tuple(int(i == j) for i in range(d))
                              for j in range(d) if j not in J}
    for i, h in zip(J, L):
        assert all(a[i] <= 2*h for a in active_heads)
        projections = {tuple(a[j] for j in J if j != i) for a in active_heads}
        assert len(projections) == len(active_heads)
    assert len(active_heads)*max(L) <= 3**(r-1)*volume
    # Check the explicit balanced-cut head construction whenever applicable.
    if max(w) <= 2*min(w):
        k = t//(4*d*min(w))
        minimum = min(range(d), key=w.__getitem__)
        others = [i for i in range(d) if i != minimum]
        constructed = set()
        for v in product(range(1, k+1), repeat=d-1):
            a = [0]*d
            for i, x in zip(others, v):
                a[i] = x
            s = dot(a, w)
            a[minimum] = (t-s)//w[minimum]+1
            a = tuple(a)
            assert a in E
            constructed.add(a)
        assert len(constructed) == k**(d-1)
        counts['balanced_head_construction_checks'] += 1
    counts['dyadic_shell_and_all_antichain_projections'] += 1


def enumerate_cuts(d, limit, counts):
    downsets = {frozenset({(0,)*d})}
    retained, rows = {}, []
    for n in range(1, limit+1):
        valid, supports = [], Counter()
        next_downsets = set()
        for D in sorted(downsets, key=lambda s: tuple(sorted(s))):
            E = border(D)
            counts['independent_full_box_border_comparisons'] += 1
            if n < limit:
                next_downsets.update(D | {a} for a in E)
            result = separator(D, E)
            if result is None:
                continue
            w, region = result
            chart_checks(D, E, w, counts)
            valid.append((D, E, w, region))
            for a in E:
                supports[tuple(sorted(D | {a}))] += 1
        if d == 3:
            coverage = sum(area2(poly) for _, _, _, poly in valid)
            assert coverage == 1, 'complete normalized weight triangle'
        else:
            intervals = sorted((region for _, _, _, region in valid), key=lambda r: r[0])
            assert intervals[0][0] == 0 and intervals[-1][1] is None
            assert all(a[1] == b[0] for a, b in zip(intervals, intervals[1:]))
            coverage = 'all positive slopes'
        counts['complete_projective_weight_coverages'] += 1
        marked = sum(len(E) for _, E, _, _ in valid)
        assert sum(supports.values()) == marked
        assert max(supports.values()) <= n+1
        best = max(valid, key=lambda v: len(v[1]))
        canonical = [(sorted(D), E, w) for D, E, w, _ in valid]
        rows.append({'dimension': d, 'length': n, 'all_downsets': len(downsets),
                     'corner_cuts': len(valid), 'marked_border_incidence': marked,
                     'distinct_supports': len(supports),
                     'support_multiplicity_histogram': dict(sorted(Counter(supports.values()).items())),
                     'max_border_witness': {'standard_set': sorted(best[0]), 'border': best[1],
                                            'strict_weight': best[2]},
                     'complete_weight_coverage': str(coverage),
                     'all_cuts_sha256': hashlib.sha256(json.dumps(canonical, separators=(',', ':')).encode()).hexdigest()})
        retained[d, n] = valid
        downsets = next_downsets
    return retained, rows


def bareiss(A):
    """Integer determinant, independent of the modular interpolation solver."""
    A = [list(row) for row in A]
    n, previous, sign = len(A), 1, 1
    for k in range(n-1):
        j = next((i for i in range(k, n) if A[i][k]), None)
        if j is None:
            return 0
        if j != k:
            A[k], A[j] = A[j], A[k]
            sign *= -1
        pivot = A[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = A[i][j]*pivot-A[i][k]*A[k][j]
                assert numerator % previous == 0
                A[i][j] = numerator//previous
            A[i][k] = 0
        previous = pivot
    return sign*A[-1][-1]


def dense_fixture(d, n, cuts, p, counts):
    # A genuine radical n-point quotient over F_p(T); numeric coefficients
    # are only the constants multiplying the explicitly retained T powers.
    M = n+1
    encode = lambda a: sum(x*M**i for i, x in enumerate(a))
    ground = set().union(*(D | set(E) for D, E, _, _ in cuts))
    assert all(max(a) < M for a in ground), 'injective Kronecker encoding'
    columns = {a: {i: pow(i+1, encode(a), p) for i in range(n)} for a in ground}
    evaluate = columns.__getitem__
    polynomials, supports, cofactors, first_polynomial = set(), set(), {}, {}
    marked = 0
    scalar_duplicates = []
    for D, E, w, _ in cuts:
        B, pivots = field_basis(tuple(sorted(D)), evaluate, p)
        assert len(B) == n
        for a in E:
            marked += 1
            support = tuple(sorted(D | {a}))
            if support not in cofactors:
                co = {}
                for j, b in enumerate(support):
                    rest = support[:j]+support[j+1:]
                    co[b] = (-1)**j*bareiss([[columns[c][i] for c in rest] for i in range(n)]) % p
                    assert co[b], ('fixture not in the required generic open set', d, n, p, support, b)
                    counts['independent_nonzero_integer_cofactor_determinants'] += 1
                cofactors[support] = co
            co = cofactors[support]
            residual, representation = field_reduce(evaluate(a), pivots, p)
            assert not residual and set(representation) == set(D)
            terms = tuple((b, (1 if b == a else -representation[b]) % p,
                           encode(a)-encode(b)) for b in support)
            assert all(c == co[b]*pow(co[a], -1, p) % p for b, c, _ in terms)
            assert all(sum(c*columns[b][i] for b, c, _ in terms) % p == 0 for i in range(n))
            assert all(dot(w, a) > dot(w, b) for b in D)
            assert all(degree+encode(b) == encode(a) for b, _, degree in terms)
            assert terms not in polynomials, 'literal monic polynomials are distinct'
            polynomials.add(terms)
            if support in supports and not scalar_duplicates:
                old_head, old_terms = first_polynomial[support]
                assert old_terms != terms and old_head != a
                scalar_duplicates.append({'common_support': support,
                                          'first_head': old_head, 'new_head': a,
                                          'first_monic_terms_exponent_coefficient_Tpower': old_terms,
                                          'new_monic_terms_exponent_coefficient_Tpower': terms})
            first_polynomial.setdefault(support, (a, terms))
            supports.add(support)
            counts['full_dense_reduced_relations_and_direct_point_substitutions'] += 1
    assert len(polynomials) == marked and len(supports) == len(cofactors)
    assert sum(len(g) for g in polynomials) == (n+1)*marked
    return {'dimension': d, 'length': n, 'characteristic': p,
            'points': 'p_i,j=(i*T)^((n+1)^j), i=1,...,n; j=0,...,d-1',
            'field': 'F_p(T)', 'monic_polynomials': len(polynomials),
            'polynomials_up_to_nonzero_scalars': len(supports),
            'explicit_term_records': (n+1)*marked,
            'one_multiple_normalization_witness': scalar_duplicates}


def bounded_compositions(r, R, maximum):
    if maximum < 0:
        return 0
    total = 0
    for k in range(r+1):
        remainder = R-k*(maximum+1)
        if remainder >= 0:
            total += (-1)**k*comb(r, k)*comb(remainder+r-1, r-1)
    return total


def dyadic_checks(counts):
    rows = []
    for r in range(1, 9):
        constant = r**(r-1)*(2**(r+1)-2)
        largest = Q(0)
        for R in range(121):
            exact = Q(0)
            for maximum in range((R+r-1)//r, R+1):
                number = (bounded_compositions(r, R, maximum)
                          -bounded_compositions(r, R, maximum-1))
                assert 0 <= number <= comb(r*maximum-R+r-1, r-1)
                exact += Q(number, 2**maximum)
            ratio = exact*2**((R+r-1)//r)
            assert ratio <= constant
            largest = max(largest, ratio)
            counts['exact_dyadic_shell_sums'] += 1
        # Independent direct tuple enumeration for small shells.
        for R in range(7):
            direct = sum((Q(1, 2**max(v)) for v in product(range(R+1), repeat=r)
                          if sum(v) == R), Q(0)) if r <= 5 else None
            if direct is not None:
                formula = sum((Q(bounded_compositions(r, R, m)-bounded_compositions(r, R, m-1), 2**m)
                               for m in range((R+r-1)//r, R+1)), Q(0))
                assert direct == formula
                counts['independent_direct_dyadic_sums'] += 1
        rows.append({'dimension': r, 'R_range': [0, 120],
                     'proof_constant': constant, 'max_observed_normalized_sum': str(largest)})
    return rows


def balanced_controls(counts):
    def score_points(w, t):
        if len(w) == 1:
            return [(i,) for i in range(t//w[0]+1)]
        return [(i, *a) for i in range(t//w[0]+1)
                for a in score_points(w[1:], t-i*w[0])]
    rows = []
    for w, t in (((1, 1), 64), ((1, 2), 64), ((1, 1, 1), 36),
                 ((1, 2, 1), 36), ((1, 1, 1, 1), 32), ((1, 2, 2, 1), 32)):
        D = frozenset(score_points(w, t))
        candidates = set(score_points(w, t+max(w)))-D
        E = tuple(sorted(a for a in candidates
                         if all(dot(w, a)-w[i] <= t for i in range(len(w)) if a[i])))
        assert list(E) == minimal_excluded_border(tuple(sorted(D)))
        k = t//(4*len(w))
        assert k >= 2, 'nonempty multidimensional lower-head controls'
        chart_checks(D, E, w, counts)
        rows.append({'weights': w, 'final_score': t, 'length': len(D),
                     'minimal_excluded_heads': len(E),
                     'explicit_lower_head_construction': k**(len(w)-1)})
        counts['nonvacuous_large_balanced_constructions'] += 1
    return rows


def triangular_extension_control(retained, counts):
    witness = None
    for n in range(1, 24):
        following = {D for D, _, _, _ in retained[2, n+1]}
        for D, E, _, _ in retained[2, n]:
            coherent = [a for a in E if D | {a} in following]
            assert 1 <= len(coherent) <= 2, 'triangular-preserving addable cells'
            if len(E) > len(coherent) and witness is None:
                witness = {'standard_set': sorted(D), 'ordinary_minimal_border': E,
                           'corner_cut_preserving_additions': coherent}
            counts['ordinary_border_versus_triangular_extension_controls'] += 1
    assert witness is not None
    return witness


def main():
    counts, retained, rows = Counter(), {}, []
    for d, limit in ((2, 24), (3, 10)):
        cuts, checks = enumerate_cuts(d, limit, counts)
        retained.update(cuts)
        rows.extend(checks)
        print('complete downset and weight enumeration:', d, 'through length', limit, flush=True)
    fixtures = []
    for d, limit in ((2, 8), (3, 6)):
        for n in range(1, limit+1):
            for p in (1000003, 1000033):
                assert p > n and all(p % k for k in range(2, int(p**0.5)+1))
                fixtures.append(dense_fixture(d, n, retained[d, n], p, counts))
    shells = dyadic_checks(counts)
    balanced = balanced_controls(counts)
    triangular = triangular_extension_control(retained, counts)
    paths = {Path(__file__).resolve(), ROOT/'docs/marked_corner_cut_output.md',
             ROOT/'paper/corner_cut_growth/main.tex'}
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name:
            path = Path(name).resolve()
            if path.suffix == '.py' and path.parent == ROOT/'scripts':
                paths.add(path)
    out = {'command': 'python3 -B scripts/run_marked_border_checks.py', 'all_checks_passed': True,
           'scope': 'Exact bounded controls for the uniform marked-border and generic dense-output proof. '
                    'The fixtures are actual radical point ideals over F_p(T). '
                    'Finite agreement is not the asymptotic proof or a novelty assessment.',
           'enumerations': rows, 'dense_point_fixtures': fixtures, 'dyadic_shell_controls': shells,
           'large_balanced_cut_controls': balanced,
           'ordinary_border_is_not_triangular_addability': triangular,
           'counters': dict(counts),
           'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in sorted(paths)},
           'unrestricted_G3_improved': False, 'polynomial_lower_separation_proved': False,
           'novelty_established': False, 'goal_complete': False}
    (ROOT/'results/research/marked_border_checks.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(dict(counts), sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
