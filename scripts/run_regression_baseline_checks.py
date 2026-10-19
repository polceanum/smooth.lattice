#!/usr/bin/env python3
"""Exact baseline controls, not an architecture or runtime benchmark.

The local/shared reader is weighted ridge regression with local decoding.
Check a baseline that reconstructs affected C*U and C*V rows on demand,
instead of keeping those two full arrays. Dense rational oracles are charged
only to this proof checker. No training, networking or third-party packages.
"""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import random

from run_retention_tradeoff_checks import identity, inverse, multiply, plus, scaled, transpose
from run_shared_feature_retention_checks import (
    BSCALE, USCALE, VSCALE, DEN, dense_polynomial, local_column,
    local_gram, positive_definite, quantize, state_from_integers,
)

ROOT = Path(__file__).resolve().parents[1]
CACHE_NAMES = ("hn", "tn", "sn", "yn")


def row_products(ci, ui, vi, i):
    """Temporary current rows; no persistent X or Z arrays."""
    support = [(j, c) for j, c in enumerate(ci[i]) if c]
    x = [sum(c*ui[j][k] for j, c in support) for k in range(len(ui[0]))]
    z = [sum(c*vi[j][k] for j, c in support) for k in range(len(vi[0]))]
    return x, z


def add_contribution(cache, b_row, u_row, x_row, z_row, sign):
    r, v = len(u_row), len(z_row)
    for j in range(r):
        for k in range(r):
            cache['hn'][j][k] += sign*u_row[j]*x_row[k]
        for k in range(v):
            cache['tn'][j][k] += sign*u_row[j]*z_row[k]
    for j, coefficient in enumerate(b_row):
        if coefficient:
            for k in range(r):
                cache['yn'][j][k] += sign*coefficient*x_row[k]
            for k in range(v):
                cache['sn'][j][k] += sign*coefficient*z_row[k]


def check_reader(cache, bi, ui, vi, ci, dense_inverse=False):
    """Check normal equations and the independently materialized decoder."""
    expected = state_from_integers(bi, ui, vi, ci)
    assert set(cache) == set(CACHE_NAMES)
    assert all(cache[name] == expected[name] for name in CACHE_NAMES)
    r = len(ui[0])
    h = scaled(multiply(inverse(cache['hn']), cache['tn']), Q(USCALE, VSCALE))
    b, u, v, c = (scaled(bi, Q(1, BSCALE)), scaled(ui, Q(1, USCALE)),
                  scaled(vi, Q(1, VSCALE)), scaled(ci, Q(1, DEN)))
    residual = plus(v, scaled(multiply(u, h), -1))
    # Stationarity of 1/2*(v-Uh)^T*C*(v-Uh) + 1/2*h^T*h.
    assert multiply(transpose(u), multiply(c, residual)) == h
    assert positive_definite(plus(identity(r), multiply(transpose(u), multiply(c, u))))
    coefficients = multiply(c, residual)
    direct = multiply(transpose([left+right for left, right in zip(b, u)]), coefficients)
    s0 = scaled(cache['sn'], Q(1, BSCALE*DEN*VSCALE))
    y = scaled(cache['yn'], Q(1, BSCALE*DEN*USCALE))
    decoder = plus(s0, scaled(multiply(y, h), -1))+h
    assert decoder == direct
    # A fresh query using every shared coordinate and multiple private ones.
    query = [[Q(1, j+2) for j in range(len(b[0])+r)]]
    assert multiply(query, decoder) == multiply(query, direct)
    if dense_inverse:
        surrogate = plus(inverse(c), multiply(u, transpose(u)))
        assert coefficients == multiply(inverse(surrogate), v)


def run_case(n, r, depth):
    rng = random.Random(5107+100*n+10*r+depth)
    v_width, key_steps, value_steps = 2, 18, 8
    active = [True]*n
    bi = scaled(identity(n), BSCALE)
    for i in range(n):
        bi[i][(i-1) % n], bi[i][(i+1) % n] = 2, 3
    ui = [[USCALE*(j == i % r)+(-1)**(i+j)*8 for j in range(r)] for i in range(n)]
    vi = [[rng.randrange(-4, 5) for _ in range(v_width)] for _ in range(n)]
    g = local_gram(scaled(bi, Q(1, BSCALE)), active)
    a = plus(identity(n), scaled(g, -1))
    ci = [[quantize(x) for x in row] for row in dense_polynomial(a, depth)]
    cache = {k: value for k, value in state_from_integers(bi, ui, vi, ci).items()
             if k in CACHE_NAMES}
    check_reader(cache, bi, ui, vi, ci, dense_inverse=True)
    max_ball, max_row_nnz = 0, 0
    operations = {'insertions': 0, 'deletions': 0, 'key_replacements': 0}
    for step in range(key_steps):
        i = (step//2) % n
        new_on = step % 6 != 0
        bnew, unew, vnew = ([row[:] for row in array] for array in (bi, ui, vi))
        onnew = active[:]
        onnew[i] = new_on
        bnew[i], unew[i], vnew[i] = [0]*n, [0]*r, [0]*v_width
        if new_on:
            bnew[i][i] = BSCALE
            bnew[i][(i-1) % n], bnew[i][(i+1) % n] = rng.randrange(1, 5), rng.randrange(1, 5)
            unew[i] = [USCALE*(j == i % r)+rng.randrange(-8, 9) for j in range(r)]
            vnew[i] = [rng.randrange(-4, 5) for _ in range(v_width)]
        gnew = local_gram(scaled(bnew, Q(1, BSCALE)), onnew)
        anew = plus(identity(n), scaled(gnew, -1))
        ball = {i}
        for _ in range(depth):
            ball |= {j for k in tuple(ball) for j in range(n)
                     if k != j and (g[k][j] or gnew[k][j])}
        cnew = [row[:] for row in ci]
        for j in ball:
            column = local_column(anew, j, depth)
            for k, coefficient in enumerate(column):
                cnew[k][j] = coefficient
        oracle = [[quantize(x) for x in row] for row in dense_polynomial(anew, depth)]
        assert cnew == oracle and cnew == transpose(cnew)
        assert positive_definite(scaled(cnew, Q(1, DEN)))
        assert all(ci[j] == cnew[j] for j in range(n) if j not in ball)
        # Both the old and new products are computed from sparse rows, not
        # read from persistent full products. A live implementation buffers
        # these local contributions before applying the edited inputs.
        for j in ball:
            old_x, old_z = row_products(ci, ui, vi, j)
            new_x, new_z = row_products(cnew, unew, vnew, j)
            add_contribution(cache, bi[j], ui[j], old_x, old_z, -1)
            add_contribution(cache, bnew[j], unew[j], new_x, new_z, 1)
        check_reader(cache, bnew, unew, vnew, cnew)
        name = ('key_replacements' if active[i] == new_on else
                'insertions' if new_on else 'deletions')
        operations[name] += 1
        max_ball = max(max_ball, len(ball))
        max_row_nnz = max(max_row_nnz, max(sum(bool(x) for x in row) for row in cnew))
        bi, ui, vi, ci, active, g = bnew, unew, vnew, cnew, onnew, gnew
    for _ in range(value_steps):
        i = rng.choice([j for j, on in enumerate(active) if on])
        replacement = [rng.randrange(-4, 5) for _ in range(v_width)]
        delta = [new-old for new, old in zip(replacement, vi[i])]
        x_i, _ = row_products(ci, ui, vi, i)
        for j in range(r):
            for k in range(v_width):
                cache['tn'][j][k] += x_i[j]*delta[k]
        for j in range(n):
            if ci[j][i]:
                for coordinate, bval in enumerate(bi[j]):
                    if bval:
                        for k in range(v_width):
                            cache['sn'][coordinate][k] += bval*ci[j][i]*delta[k]
        vi[i] = replacement
        check_reader(cache, bi, ui, vi, ci)
    check_reader(cache, bi, ui, vi, ci, dense_inverse=True)
    return {'n': n, 'shared_dimension': r, 'depth': depth, **operations,
            'value_only_overwrites': value_steps, 'max_repair_ball': max_ball,
            'max_filter_row_support': max_row_nnz,
            'persistent_product_scalars_avoided': n*(r+v_width),
            'all_current_caches_and_query_coefficients_equal': True}


def conditioning_controls():
    controls = []
    for r in (1, 2, 4, 8):
        for copies in (2, 4):
            n = copies*r
            u = [[int(j == i % r) for j in range(r)] for i in range(n)]
            h = plus(identity(r), multiply(transpose(u), u))
            c = Q(1+copies, n)
            assert h == scaled(identity(r), c*n)
            # tr(H)<=r+K*n*M^2, with C=I, K=M^2=1: equality here.
            assert r*(c*n-1) == n
            controls.append({'n': n, 'shared_dimension': r, 'c': str(c),
                             'trace_bound_saturated': True})
    return controls


def main():
    cases = [run_case(n, r, depth) for n, r in ((6, 1), (12, 3)) for depth in (0, 2, 4)]
    sources = ['scripts/run_regression_baseline_checks.py',
               'scripts/run_shared_feature_retention_checks.py',
               'scripts/run_retention_tradeoff_checks.py', 'docs/interleaving_research.md']
    report = {
        'schema_version': 1,
        'command': 'python3 -B scripts/run_regression_baseline_checks.py',
        'scope': 'Exact regression-baseline and cache-reconstruction controls; no runtime or novelty claim.',
        'cases': cases,
        'conditioning_controls': conditioning_controls(),
        'summary': {'cases': len(cases), 'key_edits': sum(sum(c[k] for k in
                    ('insertions', 'deletions', 'key_replacements')) for c in cases),
                    'value_only_overwrites': sum(c['value_only_overwrites'] for c in cases),
                    'independent_dense_surrogate_solves': 2*len(cases),
                    'conditioning_equality_controls': 8},
        'new_architecture_advantage_established': False,
        'source_sha256': {s: sha256((ROOT/s).read_bytes()).hexdigest() for s in sources},
    }
    (ROOT/'results/research/regression_baseline_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['summary'], sort_keys=True))


if __name__ == '__main__':
    main()
