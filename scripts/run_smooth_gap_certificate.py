#!/usr/bin/env python3
"""Construct and independently verify finite three-prime gap certificates.

All decisions use integers. See docs/smooth_gap_certificates.md for the proof.
This research checker does not modify or certify the fast rank/unrank kernels.
"""
from __future__ import annotations

import argparse
import bisect
import copy
import hashlib
import json
import math
import time
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path


def powers_through(p: int, x: int) -> list[int]:
    out = [1]
    while out[-1] <= x:
        out.append(out[-1] * p)
    return out


def next_two(x: int, p: int, q: int) -> tuple[int, tuple[int, int]]:
    """Smallest p^a q^b strictly above integer x, with integer comparisons."""
    ps = powers_through(p, x)
    qs = powers_through(q, x)
    best = ps[-1]
    exps = (len(ps) - 1, 0)
    b = len(qs) - 1
    for a, pa in enumerate(ps):
        while b and pa * qs[b - 1] > x:
            b -= 1
        value = pa * qs[b]
        if value < best:
            best, exps = value, (a, b)
    return best, exps


def axes(primes: tuple[int, int, int], num: int, den: int, cap: int):
    bounds, witnesses = [], []
    for i, p in enumerate(primes):
        others = [q for j, q in enumerate(primes) if j != i]
        x = 1
        for a in range(cap + 1):
            value, exps = next_two(x, *others)
            if den * value <= num * x:
                bounds.append(a)
                witnesses.append({'axis': i, 'exponent': a,
                                  'start': str(x), 'end': str(value),
                                  'other_exps': list(exps)})
                break
            x *= p
        else:
            raise RuntimeError(('axis cap exhausted', i, cap))
    return bounds, witnesses


def plane(primes, bounds, i, j, num, den):
    k = 3 - i - j
    p, q, r = primes[i], primes[j], primes[k]
    ui, uj = bounds[i], bounds[j]
    rp = powers_through(r, num * p ** (ui - 1) * q ** (uj - 1) // den)
    front, changes = [], []
    running = uj
    pa = 1
    cells = 0
    for a in range(ui):
        value = pa
        c = bisect.bisect_right(rp, value)
        for b in range(running):
            while rp[c] <= value:
                c += 1
            cells += 1
            if den * rp[c] <= num * value:
                running = b
                changes.append({'exps': [a, b], 'numerator_exponent': c})
                break
            value *= q
        front.append(running)
        pa *= p
    return {'axes': [i, j], 'frontier': front,
            'changes': changes, 'tested_cells': cells}


def next_three(x: int, primes):
    p, q, r = primes
    rp = powers_through(r, x)
    best = rp[-1]
    exps = (0, 0, len(rp) - 1)
    for c, rc in enumerate(rp):
        val, ab = next_two(x // rc, p, q)
        value = val * rc
        if value < best:
            best, exps = value, (*ab, c)
    return best, exps


def solve(primes=(2, 3, 5), num=1001, den=1000, cap=20000):
    if len(set(primes)) != 3 or not den < num < min(primes) * den:
        raise ValueError('Use three distinct primes and 1 < ratio < min(primes).')
    if not all(p >= 2 and all(p % d for d in range(2, math.isqrt(p) + 1))
               for p in primes):
        raise ValueError('Generators must be primes.')
    started = time.perf_counter()
    bounds, witnesses = axes(primes, num, den, cap)
    axis_seconds = time.perf_counter() - started
    planes = [plane(primes, bounds, i, j, num, den)
              for i, j in ((0, 1), (0, 2), (1, 2))]
    f01, f02, f12 = [d['frontier'] for d in planes]
    best, best_exps, count = 0, None, 0
    ps = [primes[0] ** a for a in range(bounds[0])]
    qs = [primes[1] ** b for b in range(bounds[1])]
    rs = [primes[2] ** c for c in range(bounds[2])]
    corners = []
    for a in range(bounds[0]):
        for b in range(f01[a]):
            top = min(f02[a], f12[b]) - 1
            count += top + 1
            if top < 0:
                continue
            value = ps[a] * qs[b] * rs[top]
            if value > best:
                best, best_exps = value, (a, b, top)
            def good(aa, bb, cc):
                return (aa < bounds[0] and bb < bounds[1] and cc < bounds[2]
                        and bb < f01[aa] and cc < f02[aa] and cc < f12[bb])
            if not good(a + 1, b, top) and not good(a, b + 1, top):
                corners.append([a, b, top])
    end, end_exps = next_three(best, primes)
    assert den * end > num * best
    last_empty_integer_start = den * end // num
    # A real gap need not contain an integer interval start. In this exceptional
    # case floor(end/ratio)+1 remains an upper bound, not a sharp threshold.
    sharp_integer = best < last_empty_integer_start
    return {'primes': list(primes), 'ratio': [num, den],
            'axis_bounds': bounds, 'axis_witnesses': witnesses,
            'planes': planes, 'bad_starts': count,
            'maximal_bad_exponents': corners,
            'last_bad_start': str(best), 'last_bad_exps': list(best_exps),
            'last_bad_end': str(end), 'last_bad_end_exps': list(end_exps),
            'real_coverage_infimum': str(Fraction(den * end, num)),
            'last_empty_integer_start': (str(last_empty_integer_start)
                                         if sharp_integer else None),
            'integer_coverage_threshold': (str(last_empty_integer_start + 1)
                                           if sharp_integer else None),
            'integer_coverage_upper_bound': str(last_empty_integer_start + 1),
            'log_start': math.log(best),
            'volume_prefix_estimate': math.log(best) ** 3 /
            (6 * math.prod(math.log(p) for p in primes)),
            'axis_seconds': axis_seconds,
            'seconds': time.perf_counter() - started}


def generators(report):
    out = []
    for i, u in enumerate(report['axis_bounds']):
        e = [0, 0, 0]
        e[i] = u
        out.append(tuple(e))
    for face in report['planes']:
        i, j = face['axes']
        for row in face['changes']:
            e = [0, 0, 0]
            e[i], e[j] = row['exps']
            out.append(tuple(e))
    return out


def covered(e, gens):
    return any(all(a >= b for a, b in zip(e, g)) for g in gens)


def interval_two_nonempty(low, high, p, q):
    """Independent interval test: divisions and bisection, no two-pointer rank."""
    qp = powers_through(q, high)
    pp = 1
    while pp <= high:
        c = bisect.bisect_right(qp, low // pp)
        if pp * qp[c] <= high:
            return True
        pp *= p
    return False


def interval_three_nonempty(low, high, primes):
    """Test (low, high] and count terms <= low when that interval is empty."""
    p, q, r = primes
    pp = powers_through(p, high)
    count_low = 0
    qb = 1
    while qb <= high:
        rc = 1
        while qb * rc <= high:
            base = qb * rc
            a = bisect.bisect_right(pp, low // base)
            count_low += a
            if base * pp[a] <= high:
                return True, None
            rc *= r
        qb *= q
    return False, count_low


def verify(report):
    """Verify a certificate without rerunning its constructor or trusting fronts.

    Negative checks inspect every uncovered coordinate-face point. Counts use
    inclusion-exclusion, and maxima use a separate finite corner theorem.
    """
    primes = report['primes']
    num, den = report['ratio']
    bounds = report['axis_bounds']
    gens = generators(report)
    def value(e):
        return math.prod(p ** a for p, a in zip(primes, e))
    assert len(set(gens)) == len(gens)
    assert all(not covered(g, gens[:i] + gens[i + 1:])
               for i, g in enumerate(gens)), 'Nonminimal generator.'
    axis_checks = plane_checks = 0
    for witness in report['axis_witnesses']:
        i = witness['axis']
        assert witness['exponent'] == bounds[i]
        others = [p for j, p in enumerate(primes) if i != j]
        x = primes[i] ** bounds[i]
        y = math.prod(p ** a for p, a in zip(others, witness['other_exps']))
        assert int(witness['start']) == x and int(witness['end']) == y
        assert x < y and den * y <= num * x
        for a in range(bounds[i]):
            low = primes[i] ** a
            assert not interval_two_nonempty(low, num * low // den, *others)
            axis_checks += 1
    assert sorted(w['axis'] for w in report['axis_witnesses']) == [0, 1, 2]
    for face in report['planes']:
        i, j = face['axes']
        k = 3 - i - j
        for row in face['changes']:
            a, b = row['exps']
            x = primes[i] ** a * primes[j] ** b
            y = primes[k] ** row['numerator_exponent']
            assert x < y and den * y <= num * x, 'Invalid positive witness.'
        rp = powers_through(primes[k],
                            num * primes[i] ** (bounds[i] - 1) *
                            primes[j] ** (bounds[j] - 1) // den)
        for a in range(bounds[i]):
            for b in range(bounds[j]):
                e = [0, 0, 0]
                e[i], e[j] = a, b
                if covered(e, gens):
                    continue
                x = primes[i] ** a * primes[j] ** b
                c = bisect.bisect_right(rp, x)
                assert den * rp[c] > num * x, 'Missing mixed witness.'
                plane_checks += 1
    assert sorted(tuple(f['axes']) for f in report['planes']) == [(0, 1), (0, 2), (1, 2)]
    mixed = [g for g in gens if sum(a > 0 for a in g) == 2]
    count = 0
    for size in range(len(mixed) + 1):
        for subset in combinations(mixed, size):
            lo = [max((g[i] for g in subset), default=0) for i in range(3)]
            count += (-1) ** size * math.prod(max(0, u - a)
                                             for u, a in zip(bounds, lo))
    assert count == report['bad_starts'], 'Incorrect total number of gaps.'
    choices = [sorted({g[i] - 1 for g in gens if g[i] > 0}) for i in range(3)]
    corners = []
    for e in product(*choices):
        if covered(e, gens):
            continue
        if all(covered(tuple(a + (i == j) for j, a in enumerate(e)), gens)
               for i in range(3)):
            corners.append(e)
    assert sorted(corners) == sorted(map(tuple, report['maximal_bad_exponents']))
    worst = max(corners, key=value)
    low = value(worst)
    assert tuple(report['last_bad_exps']) == worst and int(report['last_bad_start']) == low
    high = value(report['last_bad_end_exps'])
    assert int(report['last_bad_end']) == high and den * high > num * low
    interior, rank_low = interval_three_nonempty(low, high - 1, primes)
    assert not interior, 'Not consecutive.'
    cutoff = den * high // num
    assert report['real_coverage_infimum'] == str(Fraction(den * high, num))
    assert report['integer_coverage_upper_bound'] == str(cutoff + 1)
    if cutoff > low:
        assert report['last_empty_integer_start'] == str(cutoff)
        assert report['integer_coverage_threshold'] == str(cutoff + 1)
    else:
        assert report['last_empty_integer_start'] is None
        assert report['integer_coverage_threshold'] is None
    return {'status': 'passed', 'axis_negative_checks': axis_checks,
            'uncovered_face_negative_checks': plane_checks,
            'minimal_generators': [list(g) for g in gens],
            'maximal_bad_exponents': [list(e) for e in sorted(corners)],
            'bad_starts_by_inclusion_exclusion': count,
            'last_bad_start_exact_rank': rank_low,
            'last_pair_independently_verified': True}


def dp_values(primes, high):
    """Independent exact sequential DP, used only on bounded control cases."""
    values = [1]
    at = [0] * len(primes)
    while True:
        nxt = min(p * values[j] for p, j in zip(primes, at))
        if nxt > high:
            break
        values.append(nxt)
        for i, p in enumerate(primes):
            while p * values[at[i]] <= nxt:
                at[i] += 1
    return values


def dp_control(report):
    primes = report['primes']
    num, den = report['ratio']
    values = dp_values(primes, int(report['last_bad_end']))
    bad = [(a, b) for a, b in zip(values, values[1:]) if den * b > num * a]
    assert len(bad) == report['bad_starts']
    assert bad[-1] == (int(report['last_bad_start']), int(report['last_bad_end']))
    gens = generators(report)
    for a, b in zip(values, values[1:]):
        e, rest = [], a
        for p in primes:
            exponent = 0
            while rest % p == 0:
                rest //= p
                exponent += 1
            e.append(exponent)
        assert rest == 1
        assert covered(e, gens) == (den * b <= num * a)
    return {'status': 'passed', 'exact_dp_values': len(values),
            'gap_memberships_compared': len(values) - 1}


def clean_report(report):
    return {k: v for k, v in report.items() if k not in
            ('seconds', 'axis_seconds', 'log_start', 'volume_prefix_estimate')}


def suite(main_den):
    controls = []
    for primes in ((2, 3, 5), (2, 3, 7), (3, 5, 7), (2, 7, 11), (5, 11, 13)):
        for den in (2, 3, 5, 10, 31, 100):
            row = solve(primes, den + 1, den)
            independent = verify(row)
            control = dp_control(row)
            controls.append({'primes': list(primes), 'ratio': [den + 1, den],
                             'bad_starts': row['bad_starts'],
                             'independent_verification': independent['status'],
                             **control})
    sample = solve((2, 3, 5), 1001, 1000)
    sample['independent_verification'] = verify(sample)
    sample['dp_control'] = dp_control(sample)
    rejected = []
    missing = copy.deepcopy(sample)
    missing['planes'][0]['changes'].pop()
    wrong_witness = copy.deepcopy(sample)
    wrong_witness['planes'][0]['changes'][0]['numerator_exponent'] += 1
    for label, bad in (('missing_mixed_generator', missing),
                       ('invalid_near_product_witness', wrong_witness)):
        try:
            verify(bad)
        except AssertionError:
            rejected.append(label)
        else:
            raise AssertionError(('Invalid certificate accepted', label))
    main = solve((2, 3, 5), main_den + 1, main_den)
    main['independent_verification'] = verify(main)
    root = Path(__file__).resolve().parents[1]
    sources = ['docs/smooth_gap_certificates.md', 'scripts/run_smooth_gap_certificate.py']
    return {'schema': 'smooth.lattice.gap_certificate.v1', 'status': 'passed',
            'scope': 'Exact three-prime gap sets and last-gap certificates; novelty unestablished.',
            'main': clean_report(main), 'sequential_control': clean_report(sample),
            'varied_prime_controls': controls, 'negative_certificates_rejected': rejected,
            'source_sha256': {s: hashlib.sha256((root / s).read_bytes()).hexdigest()
                              for s in sources}}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--den', type=int, default=1000)
    ap.add_argument('--primes', type=int, nargs=3, default=[2, 3, 5])
    ap.add_argument('--out', type=Path)
    ap.add_argument('--verify-file', type=Path)
    ap.add_argument('--dp-control', action='store_true')
    ap.add_argument('--suite', action='store_true')
    args = ap.parse_args()
    if args.suite:
        report = suite(args.den)
    else:
        report = (json.loads(args.verify_file.read_text()) if args.verify_file else
                  solve(tuple(args.primes), args.den + 1, args.den))
        report = report.get('main', report)
        report['independent_verification'] = verify(report)
        if args.dp_control:
            report['dp_control'] = dp_control(report)
        report = clean_report(report)
    if args.out:
        args.out.write_text(json.dumps(report, indent=2) + '\n')
    main = report.get('main', report)
    print(json.dumps({k: main[k] for k in
                      ('primes', 'ratio', 'axis_bounds', 'bad_starts',
                       'last_bad_exps', 'last_bad_end_exps')}, indent=2))
    print('plane_cells', sum(x['tested_cells'] for x in main['planes']))
    print('generators', len(generators(main)))
    print('corners', len(main['maximal_bad_exponents']))
    print('independent_verification', main['independent_verification']['status'])
