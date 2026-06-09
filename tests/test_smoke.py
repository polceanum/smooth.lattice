#!/usr/bin/env python3
"""Dependency-free smoke tests for smooth.lattice."""
from __future__ import annotations
import heapq
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import smooth_toolkit as st


def heap_unrank(primes, n):
    primes = sorted(primes)
    heap = [(1, tuple([0] * len(primes)))]
    seen = {1}
    for rank in range(1, n + 1):
        value, exps = heapq.heappop(heap)
        if rank == n:
            return list(exps), value
        for i, p in enumerate(primes):
            nv = value * p
            if nv not in seen:
                seen.add(nv)
                ne = list(exps)
                ne[i] += 1
                heapq.heappush(heap, (nv, tuple(ne)))
    raise AssertionError('unreachable')


def assert_unrank_matches_heap(primes, n):
    got = st.unrank(primes, n)
    exp, val = heap_unrank(primes, n)
    assert got['exps'] == exp, (primes, n, got, exp)
    assert st.value(primes, got['exps']) == val
    assert st.count_le(primes, val) == n, (primes, n, val)
    assert st.prev_leq(primes, val)['rank'] == n
    assert st.next_geq(primes, val)['rank'] == n


def main():
    subprocess.run([str(ROOT / 'build.sh')], cwd=ROOT, check=True)

    assert st.value([2, 3], [2, 4]) == 324
    assert st.count_le([2, 3, 5], 12) == 10
    assert st.prev_leq([2, 3, 5], 13)['value'] == 12
    assert st.next_geq([2, 3, 5], 13)['value'] == 15

    cases = [
        ([2, 3], 10),
        ([2, 3, 5], 10),
        ([2, 3, 5], 100),
        ([2, 3, 5, 7], 200),
        ([2, 3, 5, 7, 11], 200),
        ([2, 3, 5, 7, 11, 13], 200),
    ]
    for primes, n in cases:
        assert_unrank_matches_heap(primes, n)

    fuzz_cases = [
        ([2], [1, 2, 10, 40]),
        ([3, 5], [1, 2, 15, 60]),
        ([2, 3, 5], [1, 10, 64, 127]),
        ([3, 5, 7], [7, 23, 80]),
        ([2, 5, 11, 13], [3, 31, 90]),
        ([2, 3, 5, 7, 11], [2, 37, 125]),
        ([3, 5, 7, 11, 13], [4, 42, 120]),
        ([2, 3, 5, 7, 11, 13], [1, 50, 150]),
        ([2, 3, 5, 7, 11, 13, 17], [1, 20, 75]),
    ]
    for primes, ranks in fuzz_cases:
        for n in ranks:
            assert_unrank_matches_heap(primes, n)

    cert = st.audit([2, 3, 5, 7, 11, 13], 1000, timeout=60)
    assert cert['rank_certified'] is True, json.dumps(cert, indent=2)
    assert cert['certified_count_le'] == 1000

    cert8 = st.audit([2, 3, 5, 7, 11, 13, 17, 19], 1000, timeout=60)
    assert cert8['rank_certified'] is True, json.dumps(cert8, indent=2)
    assert cert8['certified_count_le'] == 1000

    probe = subprocess.run(
        [
            str(ROOT / 'bin' / 'smooth_layer_compressed_general'),
            'count-probe-exps',
            '2,3,5',
            '1,1,0',
            '6',
        ],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    assert 'analytic_count_probe' in probe.stdout
    assert 'layer_count=6' in probe.stdout
    assert 'expected_N=6' in probe.stdout

    band = subprocess.run(
        [
            str(ROOT / 'bin' / 'smooth_layer_compressed_general'),
            'analytic-band',
            '2,3,5',
            '6',
            '10',
            '1000',
        ],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    assert 'analytic_band' in band.stdout
    assert 'target_inside=true' in band.stdout
    assert 'recovered=true' in band.stdout
    assert 'exps=[1,1,0]' in band.stdout

    corrected_band = subprocess.run(
        [
            str(ROOT / 'bin' / 'smooth_layer_compressed_general'),
            'analytic-band-corrected',
            '2,3,5',
            '6',
            '10',
            '1000',
        ],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    assert 'analytic_band_corrected' in corrected_band.stdout
    assert 'target_inside=true' in corrected_band.stdout
    assert 'recovered=true' in corrected_band.stdout
    assert 'exps=[1,1,0]' in corrected_band.stdout

    print('smoke tests passed')


if __name__ == '__main__':
    main()
