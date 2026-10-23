#!/usr/bin/env python3
"""Exact finite controls for the mean-gap transfer theorem.

All ordering, minimum-form and translation decisions use rational enclosures
of genuinely irrational weights. Decimal means are diagnostic only. Finite
controls support the proof mechanisms, not the asymptotic theorem itself.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import cmp_to_key, lru_cache
from hashlib import sha256
from itertools import product
from math import isqrt, prod
from pathlib import Path
import json


def dyadic_bounds(lo: Fraction, hi: Fraction, bits: int):
    scale = 1 << bits
    lower = lo.numerator * scale // lo.denominator
    upper = -((-hi.numerator * scale) // hi.denominator)
    return lower, upper


@lru_cache(None)
def atanh_log(z: Fraction, bits: int):
    """Enclose 2*atanh(z), 0 <= z <= 1/3, by its positive series."""
    assert 0 <= z <= Fraction(1, 3)
    if not z:
        return Fraction(0), Fraction(0)
    total, term, m = Fraction(0), z, 0
    tolerance = Fraction(1, 1 << (bits + 16))
    while True:
        total += 2 * term / (2 * m + 1)
        m += 1
        term *= z * z
        tail = 2 * term / ((2 * m + 1) * (1 - z * z))
        if tail <= tolerance:
            return total, total + tail


@lru_cache(None)
def log_integer(n: int, bits: int):
    assert n >= 2
    k = n.bit_length() - 1
    lo2, hi2 = atanh_log(Fraction(1, 3), bits)
    z = Fraction(n - (1 << k), n + (1 << k))
    lo, hi = atanh_log(z, bits)
    return k * lo2 + lo, k * hi2 + hi


@lru_cache(None)
def continued_fraction_bounds(mode: str, bits: int):
    """A cylinder for [0; 2, a_2, ...], with an infinite specified tail.

    quadratic: a_n=2; mu3: a_{n+1}=q_n; mu4: a_{n+1}=q_n^2;
    liouville: a_{n+1}=q_n^(n+1). Every tail is infinite and irrational.
    """
    pprev, p = 1, 0
    qprev, q = 0, 1
    n = 0
    while True:
        if n == 0 or mode == 'quadratic':
            a = 2
        elif mode == 'mu3':
            a = q
        elif mode == 'mu4':
            a = q * q
        elif mode == 'liouville':
            a = q ** (n + 1)
        else:
            raise ValueError(mode)
        pprev, p = p, a * p + pprev
        qprev, q = q, a * q + qprev
        n += 1
        ends = sorted((Fraction(p, q), Fraction(p + pprev, q + qprev)))
        if ends[1] - ends[0] < Fraction(1, 1 << (bits + 4)):
            return ends[0], ends[1]


class Weights:
    def __init__(self, kind, parameters):
        self.kind, self.parameters = kind, tuple(parameters)
        self.dimension = len(parameters) if kind != 'cf' else 2
        self.max_bits = 0
        self.sign_decisions = 0

    @lru_cache(None)
    def bounds(self, bits):
        scale = 1 << bits
        out = [(scale, scale)]
        if self.kind == 'radicals':
            for d in self.parameters[1:]:
                low = isqrt(d * scale * scale)
                out.append((low, low if low * low == d * scale * scale else low + 1))
        elif self.kind == 'primes':
            low0, high0 = log_integer(self.parameters[0], bits)
            for p in self.parameters[1:]:
                low, high = log_integer(p, bits)
                out.append(dyadic_bounds(low / high0, high / low0, bits))
        elif self.kind == 'cf':
            out.append(dyadic_bounds(*continued_fraction_bounds(self.parameters[0], bits), bits))
        elif self.kind == 'dependent_control':
            out.append((2 * scale, 2 * scale))
        else:
            raise ValueError(self.kind)
        assert len(out) == self.dimension and all(low > 0 for low, _ in out)
        return tuple(out)

    def enclosure(self, z, bits):
        low = high = 0
        for a, (lo, hi) in zip(z, self.bounds(bits)):
            low += a * (lo if a >= 0 else hi)
            high += a * (hi if a >= 0 else lo)
        return low, high

    @lru_cache(None)
    def sign(self, z):
        if not any(z):
            return 0
        bits = 16
        while bits <= 4096:
            low, high = self.enclosure(z, bits)
            if low > 0 or high < 0 or low == high == 0:
                self.max_bits = max(self.max_bits, bits)
                self.sign_decisions += 1
                return 1 if low > 0 else -1 if high < 0 else 0
            bits *= 2
        raise RuntimeError('Nonzero sign was not separated by the finite enclosure limit.')

    @lru_cache(None)
    def floor(self, z):
        if not any(z[1:]):
            return z[0]
        bits = 16
        while bits <= 4096:
            low, high = self.enclosure(z, bits)
            if low // (1 << bits) == high // (1 << bits):
                self.max_bits = max(self.max_bits, bits)
                return low // (1 << bits)
            bits *= 2
        raise RuntimeError('Floor was not separated by the finite enclosure limit.')

    def decimal_weights(self):
        scale = Decimal(1 << 256)
        return tuple(Decimal(lo + hi) / (2 * scale) for lo, hi in self.bounds(256))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def positive(w, z):
    s = w.sign(z)
    if s == 0:
        raise ValueError('A nonzero integer relation violates the theorem hypothesis.')
    return z if s > 0 else tuple(-a for a in z)


def minimum_form(w, height):
    """Eliminate the unit weight by exact nearest-integer optimization."""
    best = None
    for rest in product(range(-height, height + 1), repeat=w.dimension - 1):
        z = (0, *rest)
        nearest = w.floor(tuple(-a for a in z))
        candidates = {max(-height, min(height, nearest)),
                      max(-height, min(height, nearest + 1))}
        if not any(rest):
            candidates.update((-1, 1))
        for first in candidates:
            candidate = (first, *rest)
            if not any(candidate):
                continue
            candidate = positive(w, candidate)
            if best is None or w.sign(subtract(candidate, best)) < 0:
                best = candidate
    assert best is not None
    return best


def minimum_form_brute(w, height):
    best = None
    for z in product(range(-height, height + 1), repeat=w.dimension):
        if not any(z):
            continue
        candidate = positive(w, z)
        if best is None or w.sign(subtract(candidate, best)) < 0:
            best = candidate
    return best


def points_through(w, height):
    scale = 1 << 32
    limits = [height * scale // low for low, _ in w.bounds(32)]
    cutoff = (height, *(0 for _ in range(w.dimension - 1)))
    points = [e for e in product(*(range(h + 1) for h in limits))
              if w.sign(subtract(e, cutoff)) <= 0]
    points.sort(key=cmp_to_key(lambda a, b: w.sign(subtract(a, b))))
    assert all(w.sign(subtract(b, a)) > 0 for a, b in zip(points, points[1:]))
    if w.kind == 'primes':
        exact = sorted(points, key=lambda e: prod(p ** a for p, a in zip(w.parameters, e)))
        assert points == exact, 'Log-interval ordering disagrees with integer products.'
    return points


def check_case(kind, parameters, cutoff):
    w = Weights(kind, parameters)
    points = points_through(w, cutoff + 1)
    top = (cutoff, *(0 for _ in range(w.dimension - 1)))
    initial = []
    gaps = {}
    for a, b in zip(points, points[1:]):
        if w.sign(subtract(a, top)) <= 0:
            initial.append(a)
            gaps[a] = subtract(b, a)
    assert initial and (0,) * w.dimension in initial
    height = max(max(e) for e in points)
    smallest = minimum_form(w, height)
    for gap in gaps.values():
        assert w.sign(subtract(gap, smallest)) >= 0
    comparisons = []
    for bound in (2, 4):
        move = minimum_form(w, bound)
        assert move == minimum_form_brute(w, bound)
        negative = tuple(max(0, -a) for a in move)
        feasible = [e for e in initial if all(a >= b for a, b in zip(e, negative))]
        translated_count = sum(w.sign(subtract(add(e, negative), top)) <= 0 for e in points)
        assert len(feasible) == translated_count
        for e in feasible:
            assert all(a >= 0 for a in add(e, move))
            assert w.sign(subtract(gaps[e], move)) <= 0
        comparisons.append({'height': bound, 'positive_minimum_move': list(move),
                            'exact_translated_points': len(feasible),
                            'translation_count_verified': True,
                            'minimum_checked_by_full_signed_box': True})
    with localcontext() as ctx:
        ctx.prec = 65
        real = w.decimal_weights()
        def value(e):
            return sum(Decimal(a) * b for a, b in zip(e, real))
        eligible = [e for e in initial if w.sign(subtract(e, (2, *(0 for _ in range(w.dimension - 1))))) >= 0]
        reciprocal_logs = [1 / value(e).ln() for e in eligible]
        rho = [-(value(gaps[e]).exp() - 1).ln() / value(e).ln() for e in eligible]
        average = sum(rho) / len(rho)
        upper = -value(smallest).ln() * sum(reciprocal_logs) / len(reciprocal_logs)
        assert average <= upper + Decimal('1e-40')
        constant = (Decimal(1).exp() - 1).ln() / Decimal(2).ln()
        for comparison in comparisons:
            move = tuple(comparison['positive_minimum_move'])
            negative = tuple(max(0, -a) for a in move)
            fraction = Decimal(sum(all(a >= b for a, b in zip(e, negative)) for e in eligible)) / len(eligible)
            amplitude = -(value(move).exp() - 1).ln()
            assert amplitude > 0
            lower = fraction * amplitude / Decimal(cutoff).ln() - constant * (1 - fraction)
            assert average >= lower - Decimal('1e-40')
            comparison['diagnostic_lower_mean_bound'] = str(+lower)[:25]
        diagnostics = {'mean': str(+average)[:25], 'upper_mean_bound': str(+upper)[:25],
                       'scope': '65-digit Decimal diagnostics; transcendental values are not interval-certified.'}
    result = {'kind': kind, 'parameters': list(parameters), 'dimension': w.dimension,
              'cutoff': cutoff, 'ordered_prefix_points': len(initial),
              'mean_points': len(eligible), 'upper_form_height': height,
              'positive_upper_minimum_form': list(smallest),
              'exact_upper_gap_checks': len(gaps), 'comparisons': comparisons,
              'diagnostics': diagnostics,
              'exact_prime_product_order_control': kind == 'primes',
              'sign_decisions': w.sign_decisions, 'largest_sign_precision_bits': w.max_bits,
              'status': 'passed'}
    w.sign.cache_clear()
    w.floor.cache_clear()
    w.bounds.cache_clear()
    return result


def run():
    cases = []
    specifications = [('radicals', (1, 2), 80),
                      ('radicals', (1, 2, 3), 20),
                      ('radicals', (1, 2, 3, 6), 8),
                      ('primes', (2, 3), 80),
                      ('primes', (2, 3, 5), 20)]
    specifications += [('cf', (mode,), 80)
                       for mode in ('quadratic', 'mu3', 'mu4', 'liouville')]
    for args in specifications:
        cases.append(check_case(*args))
    dependent = Weights('dependent_control', (1, 2))
    try:
        minimum_form(dependent, 2)
    except ValueError:
        dependent_rejected = True
    else:
        raise AssertionError('Rational dependence was not rejected.')
    root = Path(__file__).resolve().parents[1]
    files = ['docs/mean_gap_diophantine_exponents.md', 'scripts/run_mean_gap_checks.py']
    return {'schema': 'smooth.lattice.mean_gap_checks.v1', 'status': 'passed',
            'verification_scope': 'Exact finite irrational ordering, minimum-form and translation checks. The limiting theorem is proved analytically, not by these computations.',
            'cases': cases, 'dependent_weights_rejected': dependent_rejected,
            'exact_upper_gap_checks': sum(c['exact_upper_gap_checks'] for c in cases),
            'exact_translation_gap_checks': sum(r['exact_translated_points'] for c in cases for r in c['comparisons']),
            'source_sha256': {s: sha256((root / s).read_bytes()).hexdigest() for s in files}}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, default=Path('results/research/mean_gap_checks.json'))
    args = parser.parse_args()
    result = run()
    args.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases', 'source_sha256')}, indent=2))
