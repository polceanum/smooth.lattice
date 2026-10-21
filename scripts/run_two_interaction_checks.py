#!/usr/bin/env python3
"""Exact columns and ranks for two interactions; no asymptotic fan claim."""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import gcd
from pathlib import Path
import hashlib
import json

from run_higher_interaction_checks import ground, linear_sum
from run_higher_seed_checks import column_rank
from run_hypersurface_truncation_checks import relation_rows
from run_structural_complexity_checks import determinant

ROOT = Path(__file__).resolve().parents[1]


def binary_mul(f, g):
    out = defaultdict(int)
    for (a, b), c in f.items():
        for (i, j), d in g.items():
            out[a+i, b+j] += c*d
    return {v: c for v, c in out.items() if c}


def power(f, e):
    out = {(0, 0): 1}
    for _ in range(e):
        out = binary_mul(out, f)
    return out


def binary_degree(d):
    return tuple((a, d-a) for a in range(d+1))


def model(r, Q, H, L):
    s = 2*r+1
    basis = tuple(v for v in ground(s) if v[2] < 2)
    index = {v: i for i, v in enumerate(basis)}

    @lru_cache(None)
    def evaluate(v):
        x, y, z = v
        if sum(v) > s:
            return {}
        if z < 2:
            return {index[v]: 1}
        # This recurrence does not use the first-order column formula.
        return linear_sum([(-c, evaluate((x+a, y+b, z-2)))
                           for (a, b), c in Q.items()]
                          +[(-c, evaluate((x+a, y+b, z+r-2)))
                            for (a, b), c in H.items()]
                          +[(-c, evaluate((x+a, y+b, z+r-1)))
                            for (a, b), c in L.items()])

    @lru_cache(None)
    def qpower(e):
        return power({v: -c for v, c in Q.items()}, e)

    def formula(v):
        x, y, z = v
        j, eps = divmod(z, 2)
        k = (r-1)//2
        terms = [(1, qpower(j), eps)]
        if j:
            terms += [(-j, binary_mul(H, qpower(j+k-1+eps)), 1-eps),
                      (-j, binary_mul(L, qpower(j+k)), eps)]
        out = defaultdict(int)
        for coefficient, f, ztail in terms:
            for (a, b), c in f.items():
                v = (x+a, y+b, ztail)
                if sum(v) <= s:
                    out[index[v]] += coefficient*c
        return {v: c for v, c in out.items() if c}

    equation = {(0, 0, 2): 1}
    for f, z in ((Q, 0), (H, r), (L, r+1)):
        equation.update({(a, b, z): c for (a, b), c in f.items()})
    return basis, evaluate, formula, equation


def divides(f, g, p):
    df, dg = sum(next(iter(f))), sum(next(iter(g)))
    coordinates = binary_degree(dg)
    idx = {v: i for i, v in enumerate(coordinates)}
    columns = [{idx[a+i, b+j]: c for (a, b), c in f.items()}
               for i, j in binary_degree(dg-df)]
    target = {idx[v]: c for v, c in g.items()}
    return column_rank(columns, p) == column_rank(columns+[target], p)


def coupled_columns(r, Q, H, L):
    """Independent sparse matrix of the equations (C9), without Q powers."""
    rows = [(0, d, a) for d in range(1, r) for a in range(d+3)]
    rows += [(1, d, a) for d in range(1, r-1) for a in range(d+3)]
    idx = {v: i for i, v in enumerate(rows)}
    out = []
    QL = binary_mul(Q, L)
    for d in range(1, r):
        for a, b in binary_degree(d):
            col = {idx[0, d, a+i]: c for (i, j), c in H.items()}
            if d <= r-2:
                col.update({idx[1, d, a+i]: c for (i, j), c in L.items()})
            out.append(col)
    for e in range(r-2):
        for a, b in binary_degree(e):
            col = {idx[0, e+2, a+i]: -c for (i, j), c in QL.items()}
            if e >= 1:
                col.update({idx[1, e, a+i]: c for (i, j), c in H.items()})
            out.append(col)
    assert len(out) == r*r-r
    return out


def quotient_binary(f, g, p):
    """Solve f*h=g for a homogeneous quadratic h by separate elimination."""
    df, dg = sum(next(iter(f))), sum(next(iter(g)))
    monomials = binary_degree(dg-df)
    rows = []
    for a, b in binary_degree(dg):
        row = [f.get((a-i, b-j), 0) for i, j in monomials]+[g.get((a, b), 0)]
        rows.append([c % p if p else Fraction(c) for c in row])
    pivot = 0
    for j in range(len(monomials)):
        hit = next(i for i in range(pivot, len(rows)) if rows[i][j])
        rows[pivot], rows[hit] = rows[hit], rows[pivot]
        inv = pow(rows[pivot][j], -1, p) if p else 1/rows[pivot][j]
        rows[pivot] = [c*inv % p if p else c*inv for c in rows[pivot]]
        for i, row in enumerate(rows):
            if i != pivot and row[j]:
                c = row[j]
                rows[i] = [(x-c*y) % p if p else x-c*y for x, y in zip(row, rows[pivot])]
        pivot += 1
    assert all(not any(row[:-1]) and not row[-1] for row in rows[pivot:])
    return {v: rows[i][-1] for i, v in enumerate(monomials) if rows[i][-1]}


def relation_witness(r, Q, H, L, p, basis, ev):
    D = quotient_binary(H, binary_mul(Q, L), p)
    P = {(0, 0): 1}
    if r >= 5:
        C = quotient_binary(power(H, 2), binary_mul(Q, power(L, 2)), p)
        for j in range(1, (r-1)//2):
            for v, c in power(C, j).items():
                P[v] = P.get(v, 0)+(-1)**j*c
    polynomial = {(a, b, 2): c for (a, b), c in binary_mul(D, P).items()}
    polynomial.update({(a, b, 3): c for (a, b), c in P.items()})
    nf = linear_sum([(c, ev(v)) for v, c in polynomial.items()])
    for i, c in nf.items():
        v = basis[i]
        assert v[0]+v[1]+2*v[2] <= r+3 or (p and c % p == 0)
        polynomial[v] = polynomial.get(v, 0)-c
    polynomial = {v: (c % p if p else c) for v, c in polynomial.items()
                  if (c % p if p else c)}
    assert polynomial and all(a+b+2*c <= r+3 for a, b, c in polynomial)
    remainder = linear_sum([(c, ev(v)) for v, c in polynomial.items()])
    assert all(c % p == 0 if p else c == 0 for c in remainder.values())
    return [{'exponent': v, 'coefficient': str(c)} for v, c in sorted(polynomial.items())]


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    fixtures = [
        ('independent', {(2, 0): 1, (0, 2): 1}, {(1, 1): 1}, {(2, 0): 1}),
        ('repeated', {(2, 0): 1, (0, 2): 1}, {(1, 1): 1}, {(1, 1): 1}),
        ('characteristic_two', {(2, 0): 1, (0, 2): 1}, {(2, 0): 1, (0, 2): 1}, {(2, 0): 1, (1, 1): 1}),
        ('first_test_only', {(2, 0): 1}, {(2, 0): 1}, {(0, 2): 1}),
        ('mixed_coefficients', {(2, 0): 1, (1, 1): 1, (0, 2): 1}, {(2, 0): 1, (1, 1): -1}, {(2, 0): 1, (1, 1): 2, (0, 2): 1}),
    ]
    counts, records, witnesses = defaultdict(int), [], []
    for r in range(3, 102, 2):
        V = ((2, -2, 0), (2, 0, -2), (1, 1, r-2), (2, 0, r-1))
        minors = [int(determinant(m)) for m in combinations(V, 3)]
        assert minors == [4*r, 4*(r+1), 4, 2*(r+1)] and gcd(*minors) == 4
        counts['full_rank_constant_index_support_lattices'] += 1
    for r, (label, Q, H, L) in product((3, 5, 7, 9, 11), fixtures):
        s, T = 2*r+1, r+3
        B, ev, formula, equation = model(r, Q, H, L)
        E = tuple(v for v in ground(T) if v[0]+v[1]+2*v[2] <= T)
        selected = set(E)
        full = ground(s)
        assert len(B) == (s+1)**2 and all(ev(v) == {i: 1} for i, v in enumerate(B))
        for v in full:
            assert ev(v) == formula(v)
            counts['integer_recurrence_and_closed_column_equalities'] += 1
        M = coupled_columns(r, Q, H, L)
        for p in (0, 2, 3, 5):
            first = divides(H, binary_mul(Q, L), p)
            second = divides(power(H, 2), binary_mul(Q, power(L, 2)), p)
            delta = int(first if r == 3 else second)
            rank = column_rank(map(ev, E), p)
            coupled_rank = column_rank(M, p)
            assert rank == 2*r*r+6*r+13-delta
            assert coupled_rank == r*r-r-delta
            if label == 'independent':
                assert delta == 0
            if label == 'repeated':
                assert delta == 1
            if label == 'characteristic_two':
                assert delta == int(r == 3 or p == 2)
            if label == 'first_test_only':
                assert first and delta == int(r == 3)
            counts['full_prefix_ranks_and_coupled_system_equalities'] += 1
            record = {'family': label, 'r': r, 'cutoff': s, 'length': len(B),
                      'characteristic': p, 'prices': [1, 1, 2], 'threshold': T,
                      'selected_columns': len(E), 'direct_rank': rank,
                      'coupled_system_rank': coupled_rank, 'delta': delta,
                      'H_divides_QL': first, 'H_squared_divides_QL_squared': second}
            if delta:
                witness = relation_witness(r, Q, H, L, p, B, ev)
                counts['explicit_extra_relations_replayed'] += 1
                if r <= 5 and label in ('repeated', 'characteristic_two'):
                    witnesses.append({**record, 'relation': witness})
            # Independent full Macaulay kernel; no normal-form formula in construction.
            if r == 3 or (r == 5 and label in ('independent', 'repeated')):
                rows = relation_rows(full, [equation], p)
                assert len(rows) == len(full)-len(B)
                for row in rows:
                    rem = linear_sum([(c, ev(v)) for v, c in zip(full, row) if c])
                    assert all(c % p == 0 if p else c == 0 for c in rem.values())
                    counts['independent_macaulay_kernel_rows_replayed'] += 1
                outside = [i for i, v in enumerate(full) if v not in selected]
                projected = [{j: row[i] for j, i in enumerate(outside) if row[i]} for row in rows]
                independent_rank = len(E)-len(rows)+column_rank(projected, p)
                assert rank == independent_rank
                counts['independent_macaulay_prefix_ranks'] += 1
                record['independent_macaulay_rank'] = independent_rank
            records.append(record)
    sources = ('docs/two_interaction_ranks.md', 'scripts/run_two_interaction_checks.py',
               'docs/higher_seed_transport.md', 'scripts/run_higher_interaction_checks.py',
               'scripts/run_higher_seed_checks.py', 'scripts/run_hypersurface_truncation_checks.py',
               'scripts/run_socle_stability_checks.py', 'scripts/run_structural_complexity_checks.py')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_two_interaction_checks.py',
              'counts': dict(counts), 'fixtures': records, 'exact_extra_relations': witnesses,
              'uniform_results': {'length': '4*(r+1)^2',
                                  'complete_prefix_rank': '2*r^2+6*r+13-delta',
                                  'r_3_kernel_condition': 'H divides Q*L',
                                  'odd_r_at_least_5_kernel_condition': 'H^2 divides Q*L^2',
                                  'repeated_form_rank_universality': False,
                                  'two_independent_form_rank_gap': 'One, for every odd r>=3 and every field.'},
              'unrestricted_G3_improvement': False, 'fan_count_established': False,
              'novelty_established': False, 'goal_complete': False,
              'scope': 'Uniform exact column and coupled-kernel proofs; finite independent checks. A rank difference does not give a larger total fan or establish publication significance.',
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/'results/research/two_interaction_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['counts'], indent=2))


if __name__ == '__main__':
    main()
