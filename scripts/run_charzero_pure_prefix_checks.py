#!/usr/bin/env python3
"""Exact pure-prefix ranks, strip transfers, fine flows and boundary orders."""
from collections import Counter, deque
from fractions import Fraction as Q
from itertools import product, combinations, groupby
from math import comb
from pathlib import Path
from random import Random
import hashlib
import json

from run_frobenius_prefix_checks import Field, rank
from run_generic_hook_checks import matching_basis, exchange_normals, bareiss
from run_hook_bridge_checks import model, actual, polygon
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_structural_complexity_checks import dot, area2

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def matching_rank(columns):
    mate = {}
    def augment(j, seen):
        for i in columns[j]:
            if i in seen:
                continue
            seen.add(i)
            if i not in mate or augment(mate[i], seen):
                mate[i] = j
                return True
        return False
    return sum(augment(j, set()) for j in range(len(columns)))


def arm(A, t, D):
    C = {(1, i, j): {i*t+j: 1} for i, j in product(range(A), range(t))}
    for i, j in product(range(t, A), range(t)):
        v = {i*t+j: 1}
        if j+1<t:
            v.update({(i+k)*t+j+1: comb(D, k) for k in range(D+1) if i+k<A})
        C[0, i, j] = v
    return C


def pure_prefix_controls(rng):
    out = []
    p = 1000003
    for A, t, D in [(4, 2, 0), (6, 2, 1), (8, 3, 3), (12, 4, 6), (20, 5, 10)]:
        C = arm(A, t, D)
        cases = [(Q(0), Q(1)), (Q(1, 3), Q(1, 7)), (Q(-A), Q(2))]
        cases += [(Q(rng.randrange(-8*A, 8*A), 7)+Q(1, 101),
                   Q(rng.randrange(1, 8*A), 11)+Q(1, 103)) for _ in range(24)]
        for theta, rho in cases:
            cost = lambda e: theta*e[0]+e[1]+rho*e[2]
            E = sorted(C, key=lambda e: (cost(e), e))
            selected = []
            for _, group in groupby(E, cost):
                selected.extend(group)
                vectors = [C[e] for e in selected]
                mr = matching_rank(vectors)
                ar = rank([{i: c % p for i, c in v.items() if c % p} for v in vectors], Field(p))
                if ar != mr:
                    ar = rank(vectors, Field(0))
                    COUNTS['exact_rational_modular_fallbacks'] += 1
                assert ar == mr
                COUNTS['independent_entire_pure_prefix_ranks'] += 1
            B, _ = matching_basis(E, C, A*t)
            assert rank([C[e] for e in B], Field(0)) == A*t
            COUNTS['exact_rational_pure_minimum_bases'] += 1
        out.append({'A': A, 't': t, 'degree': D, 'weights': len(cases)})
    return out


def transfer(B, t, D, rho, lo, hi):
    V = list(product(range(B), range(t)))
    U = {v for v in V if v[0]+rho*v[1]<lo}
    S = {v for v in V if v[0]+rho*v[1]<hi}
    O, C = set(V)-S, S-U
    inputs = sorted(U, key=lambda v: (v[1], -v[0]))
    outputs = sorted(O, key=lambda v: (v[1], -v[0]))
    P = []
    for u in inputs:
        paths, col = {u: 1}, {}
        for j in range(u[1], t-1):
            for q in range(B):
                value = paths.get((q, j), 0)
                if not value:
                    continue
                for k in range(D+1):
                    z = (q+k, j+1)
                    if z[0]>=B:
                        continue
                    c = value*comb(D, k)
                    if z in C:
                        paths[z] = paths.get(z, 0)+c
                    elif z in O:
                        col[z] = col.get(z, 0)+c
        P.append([col.get(v, 0) for v in outputs])
    R = set(V)-U
    M = []
    for q, j in sorted(S):
        v = {(q, j): 1} if (q, j) in R else {}
        if j+1<t:
            v.update({(q+k, j+1): comb(D, k) for k in range(D+1) if (q+k, j+1) in R})
        M.append(v)
    return U, C, O, inputs, outputs, P, M


def fine_flow(B, t, D, U, C, O, inputs=None, outputs=None):
    """Independent vertex-capacitated flow in the refined grid."""
    inputs = U if inputs is None else inputs
    outputs = O if outputs is None else outputs
    N = B*((t-1)*D+1)
    source, sink = 2*N, 2*N+1
    G = [[] for _ in range(2*N+2)]
    def edge(a, b, capacity=1):
        G[a].append([b, capacity, len(G[b])])
        G[b].append([a, 0, len(G[a])-1])
    def node(k, q):
        return k*B+q
    for k in range((t-1)*D+1):
        for q in range(B):
            i = node(k, q)
            edge(2*i, 2*i+1)
            if k%D == 0:
                v = (q, k//D)
                if v in inputs:
                    edge(source, 2*i)
                if v in outputs:
                    edge(2*i+1, sink)
    for k in range((t-1)*D):
        for q in range(B):
            if k%D == 0 and (q, k//D) in O:
                continue
            for r in (q, q+1):
                if r>=B or ((k+1)%D == 0 and (r, (k+1)//D) in U):
                    continue
                edge(2*node(k, q)+1, 2*node(k+1, r))
    total = 0
    while True:
        previous = {source: None}
        queue = deque([source])
        while queue and sink not in previous:
            i = queue.popleft()
            for j, e in enumerate(G[i]):
                if e[1] and e[0] not in previous:
                    previous[e[0]] = (i, j)
                    queue.append(e[0])
        if sink not in previous:
            return total
        i = sink
        while i != source:
            a, j = previous[i]
            e = G[a][j]
            e[1] -= 1
            G[i][e[2]][1] += 1
            i = a
        total += 1


def network_controls(rng):
    for B, t, D in [(4, 3, 1), (5, 3, 2), (6, 4, 3), (8, 5, 4)]:
        for _ in range(65):
            rho = Q(rng.randrange(1, 50), rng.randrange(1, 20))
            lo = Q(rng.randrange(0, 10*(B+int(rho)*t+1)), 10)+Q(1, 101)
            hi = lo+Q(rng.randrange(1, 10*B), 10)+Q(1, 103)
            U, C, O, I, J, P, M = transfer(B, t, D, rho, lo, hi)
            pr = rank([{j: c for j, c in enumerate(v) if c} for v in P], Field(0))
            mr = matching_rank(M)
            ar = rank(M, Field(0))
            flow = fine_flow(B, t, D, U, C, O)
            assert ar == len(C)+pr == len(C)+flow == mr
            COUNTS['independent_exact_schur_rank_and_fine_flow_controls'] += 1
            for a, b in combinations(range(len(I)), 2):
                for c, d in combinations(range(len(J)), 2):
                    assert P[a][c]*P[b][d]-P[a][d]*P[b][c] >= 0
                    COUNTS['exact_nonnegative_transfer_two_minors'] += 1
            if min(len(I), len(J)) < 3:
                continue
            for trial in range(14):
                k = rng.randrange(3, min(len(I), len(J), 7)+1)
                rows = sorted(rng.sample(range(len(I)), k))
                cols = sorted(rng.sample(range(len(J)), k))
                det = bareiss([[P[i][j] for j in cols] for i in rows])
                assert det >= 0
                COUNTS['exact_nonnegative_transfer_higher_minors'] += 1
                if trial == 0:
                    f = fine_flow(B, t, D, U, C, O, {I[i] for i in rows}, {J[j] for j in cols})
                    assert (det > 0) == (f == k)
                    COUNTS['independent_selected_minor_flow_controls'] += 1


def routing_controls(rng):
    for _ in range(700):
        D = rng.randrange(1, 18)
        size = rng.randrange(1, 12)
        a = sorted(rng.sample(range(45), size))
        b = sorted(x+rng.randrange(D+1) for x in a)
        if len(set(b)) != size or not all(x<=y<=x+D for x, y in zip(a, b)):
            continue
        paths = [[max(x, y-D+k) for k in range(D+1)] for x, y in zip(a, b)]
        assert [p[0] for p in paths] == a and [p[-1] for p in paths] == b
        assert all(p[k+1]-p[k] in (0, 1) for p in paths for k in range(D))
        assert all(len({p[k] for p in paths}) == size for k in range(D+1))
        COUNTS['independent_ordered_local_routes'] += 1


def pure_fan(A, t, D):
    # Original positive labels give theta=(w_x-t*w_z)/w_y, rho=w_z/w_y.
    C = {((x, t+i, j) if x else (0, t+i, t+j)): v for (x, i, j), v in arm(A, t, D).items()}
    r = A*t
    todo, cells = [((1, 2*r+1, (2*r+1)**2), None)], {}
    while todo:
        w, old = todo.pop()
        B, _ = matching_basis(sorted(C, key=lambda e: (dot(w, e), e)), C, r)
        assert B != old
        if B in cells:
            continue
        assert rank([C[e] for e in B], Field(0)) == r
        N = exchange_normals(B, C, r)
        assert all(dot(w, q)>0 for q in N)
        P = polygon(N)
        cells[B] = P
        for v, u in zip(P, P[1:]+P[:1]):
            if any(v[k] == u[k] == 0 for k in range(3)):
                continue
            q = next(q for q in N if dot(q, v) == dot(q, u) == 0)
            nw, _ = exact_adjacent_weight(4*(A+t)*t, v, u, q)
            todo.append((nw, B))
    assert sum(area2(P) for P in cells.values()) == 1
    COUNTS['independent_complete_pure_basis_fans'] += 1
    return {'A': A, 't': t, 'degree': D, 'rank': r, 'pure_cones': len(cells), 'exact_total_area': '1'}


def original_prefix_controls(rng):
    for t, A in [(2, 5), (3, 8), (4, 12)]:
        a, D = A+t, A//2
        data = model(t, a)
        H = {i: comb(D, i) for i in range(D+1)}
        G = {i: rng.randrange(-3, 4) for i in range(D+1)}
        terms = {(t, 0): 1, (0, t): 1}
        terms.update({(t+i, 1): c for i, c in H.items()})
        terms.update({(1, t+i): c for i, c in G.items() if c})
        co = [terms.get((y, z), 0) if g == 1 else 0 for g, y, z in data[2]]
        C = actual(data, co, 0, t<=3)
        Y, Z, heads = data[4]
        for rho in [Q(A), Q(A, 2), Q(2*A)]:
            w = (4*A*rho+Q(1, 101), Q(1), rho)
            bound = min(dot(w, e) for e in Z if C[e])
            E = sorted([e for e in C if C[e] and dot(w, e)<bound], key=lambda e: (dot(w, e), e))
            prefix = []
            for _, group in groupby(E, key=lambda e: dot(w, e)):
                prefix.extend(group)
                assert not set(prefix) & set(Z)
                vectors = [C[e] for e in prefix]
                assert rank(vectors, Field(0)) == matching_rank(vectors)
                COUNTS['independent_complete_original_one_sided_prefixes'] += 1


def main():
    if not __debug__:
        raise RuntimeError('Assertions must remain enabled.')
    previous_path = 'results/research/frobenius_prefix_checks.json'
    previous = json.loads((ROOT/previous_path).read_text())
    assert previous['status'] == 'passed'
    for p, h in previous['source_sha256'].items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == h
    rng = Random(2026102341)
    prefixes = pure_prefix_controls(rng)
    network_controls(rng)
    routing_controls(rng)
    fans = [pure_fan(5, 2, 2), pure_fan(7, 3, 3)]
    original_prefix_controls(rng)
    files = set(previous['source_sha256']) | {previous_path,
             'docs/charzero_pure_prefix_matching.md', 'scripts/run_charzero_pure_prefix_checks.py'}
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_charzero_pure_prefix_checks.py',
              'uniform_result': {
                  'class': 'Pure arm columns U=q^i v^j, L=(1+vH(q))*q^i v^j for i>=t, in K[q,v]/(q^A,v^t), characteristic zero.',
                  'rank_rule': 'H=(1+q)^D attains support matching rank for every strict cost prefix, simultaneously for all theta and positive rho.',
                  'generic_consequence': 'A nonempty coefficient-open degree stratum has the same all-prefix ranks and the same entire pure minimum-basis fan as support matching.',
                  'original_prefix_consequence': 'A complete original prefix with no other-arm upper or pure lower columns has the same rank, retaining every selected head and arbitrary other-arm coefficients.',
                  'missing_step': 'Simultaneous head dependencies when both arms contribute other columns; no whole hook-fan count or unrestricted G_3 bound.'},
              'counts': dict(COUNTS), 'pure_prefix_fixtures': prefixes,
              'complete_pure_fans': fans,
              'limits': ['Path determinants and their boundary-positivity principle are classical.',
                         'Modular rank checks attain the structural upper bound and therefore lift to characteristic zero; they do not assert a universal positive-characteristic statement.',
                         'The full coupled original fan is not replaced by its support fan.',
                         'Exact combined priority and the original unrestricted endpoint remain open.'],
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(files)}}
    (ROOT/'results/research/charzero_pure_prefix_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print('Characteristic-zero pure-prefix checks passed', json.dumps(dict(COUNTS), sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
