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
        got = st.unrank(primes, n)
        exp, val = heap_unrank(primes, n)
        assert got['exps'] == exp, (primes, n, got, exp)
        assert st.value(primes, got['exps']) == val

    cert = st.audit([2, 3, 5, 7, 11, 13], 1000, timeout=60)
    assert cert['rank_certified'] is True, json.dumps(cert, indent=2)
    assert cert['certified_count_le'] == 1000

    print('smoke tests passed')


if __name__ == '__main__':
    main()
