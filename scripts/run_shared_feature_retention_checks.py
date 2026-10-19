#!/usr/bin/env python3
"""Exact controls for local correction with a dense shared feature component.

This is a proof checker, with full oracles; not a sublinear runtime benchmark.
The Woodbury identity and polynomial inverse approximation are known methods.
No neural training, dependencies beyond Python, tolerances or network calls.
"""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import random

from run_retention_tradeoff_checks import (
    identity, inverse, multiply, norm_inf, plus, scaled, transpose,
)

ROOT = Path(__file__).resolve().parents[1]
BITS, BSCALE, USCALE, VSCALE = 20, 128, 64, 4
DEN = 1 << BITS
RHO, K, M2 = Q(1, 4), Q(2), Q(3, 2)


def zeros(n, m):
    return [[0]*m for _ in range(n)]


def positive_definite(a):
    """Exact LDL pivots, used only on small shared-component certificates."""
    a = [[Q(x) for x in row] for row in a]
    assert a == transpose(a)
    for j in range(len(a)):
        pivot = a[j][j]
        if pivot <= 0:
            return False
        for i in range(j+1, len(a)):
            for k in range(j+1, len(a)):
                a[i][k] -= a[i][j]*a[j][k]/pivot
    return True


def local_gram(b, active):
    g = multiply(b, transpose(b))
    for i, on in enumerate(active):
        g[i][i] += not on
    return g


def dense_polynomial(a, depth):
    p, power = identity(len(a)), identity(len(a))
    for _ in range(depth):
        # Dense oracle output, independently using sparse left multiplication.
        power = [[sum((x*power[j][k] for j, x in enumerate(row) if x), Q(0))
                  for k in range(len(a))] for row in a]
        p = plus(p, power)
    return p


def quantize(x):
    x = x*DEN
    return x.numerator//x.denominator


def local_column(a, target, depth):
    columns = [{i: row[j] for i, row in enumerate(a) if row[j]}
               for j in range(len(a))]
    power, total = {target: Q(1)}, {target: Q(1)}
    for _ in range(depth):
        following = {}
        for j, x in power.items():
            for i, coefficient in columns[j].items():
                following[i] = following.get(i, Q(0))+coefficient*x
        power = {i: x for i, x in following.items() if x}
        for i, x in power.items():
            total[i] = total.get(i, Q(0))+x
    return [quantize(total.get(i, Q(0))) for i in range(len(a))]


def state_from_integers(bi, ui, vi, ci):
    n, r = len(ui), len(ui[0])
    z, x = multiply(ci, vi), multiply(ci, ui)
    hn = plus(scaled(identity(r), DEN*USCALE**2), multiply(transpose(ui), x))
    tn = multiply(transpose(ui), z)
    sn, yn = multiply(transpose(bi), z), multiply(transpose(bi), x)
    return {"z": z, "x": x, "hn": hn, "tn": tn, "sn": sn, "yn": yn}


def run_case(n, r, depth):
    rng = random.Random(901+100*n+10*r+depth)
    value_width, steps = 2, 12
    active = [True]*n
    bi = scaled(identity(n), BSCALE)
    for i in range(n):
        bi[i][(i-1) % n] = bi[i][(i+1) % n] = 3
    ui = [[USCALE*(j == i % r)+USCALE//8 for j in range(r)] for i in range(n)]
    vi = [[rng.randrange(-4, 5) for _ in range(value_width)] for _ in range(n)]
    b, u = scaled(bi, Q(1, BSCALE)), scaled(ui, Q(1, USCALE))
    g = local_gram(b, active)
    a = plus(identity(n), scaled(g, -1))
    ci = [[quantize(x) for x in row] for row in dense_polynomial(a, depth)]
    state = state_from_integers(bi, ui, vi, ci)
    c_lower = Q(1, 8*r)
    gamma = 1+M2*K/c_lower
    bound = (RHO**(depth+1)+(1+RHO)*n/Q(DEN))*gamma
    max_error, max_ball, max_changed = Q(0), 0, 0
    counts = {"insertions": 0, "deletions": 0, "replacements": 0}
    for step in range(steps):
        target = step % n
        if step == 1:
            target = 0
        new_active = step != 0
        updated_bi, updated_ui, updated_vi = ([row[:] for row in t] for t in (bi, ui, vi))
        updated_active = active[:]
        updated_active[target] = new_active
        updated_bi[target] = [0]*n
        updated_ui[target] = [0]*r
        updated_vi[target] = [0]*value_width
        if new_active:
            updated_bi[target][target] = BSCALE
            for j in ((target-1) % n, (target+1) % n):
                updated_bi[target][j] = rng.randrange(1, 5)
            updated_ui[target] = [USCALE*(j == target % r)+USCALE//8+rng.randrange(-1, 2)
                                  for j in range(r)]
            updated_vi[target] = [rng.randrange(-4, 5) for _ in range(value_width)]
        updated_b = scaled(updated_bi, Q(1, BSCALE))
        updated_g = local_gram(updated_b, updated_active)
        updated_a = plus(identity(n), scaled(updated_g, -1))
        assert norm_inf(updated_a) <= RHO
        expected_ci = [[quantize(x) for x in row]
                       for row in dense_polynomial(updated_a, depth)]
        ball = {target}
        for _ in range(depth):
            ball |= {j for i in tuple(ball) for j in range(n)
                     if i != j and (g[i][j] or updated_g[i][j])}
        changed = {j for j in range(n) if any(expected_ci[i][j] != ci[i][j] for i in range(n))}
        assert changed <= ball
        assert all(expected_ci[i][j] == ci[i][j]
                   for i in range(n) for j in range(n) if i not in ball or j not in ball)
        updated_ci = [row[:] for row in ci]
        for j in ball:
            column = local_column(updated_a, j, depth)
            for i, coefficient in enumerate(column):
                updated_ci[i][j] = coefficient
        assert updated_ci == expected_ci
        expected = state_from_integers(updated_bi, updated_ui, updated_vi, updated_ci)
        new_z, new_x = [row[:] for row in state['z']], [row[:] for row in state['x']]
        for i in ball:
            nonzero = [(j, v) for j, v in enumerate(updated_ci[i]) if v]
            new_z[i] = [sum(v*updated_vi[j][k] for j, v in nonzero) for k in range(value_width)]
            new_x[i] = [sum(v*updated_ui[j][k] for j, v in nonzero) for k in range(r)]
        assert new_z == expected['z'] and new_x == expected['x']
        for i in ball:
            for j in range(r):
                for k in range(r):
                    state['hn'][j][k] += updated_ui[i][j]*new_x[i][k]-ui[i][j]*state['x'][i][k]
                for k in range(value_width):
                    state['tn'][j][k] += updated_ui[i][j]*new_z[i][k]-ui[i][j]*state['z'][i][k]
            cells = {j for j in range(n) if bi[i][j] or updated_bi[i][j]}
            for j in cells:
                for k in range(r):
                    state['yn'][j][k] += updated_bi[i][j]*new_x[i][k]-bi[i][j]*state['x'][i][k]
                for k in range(value_width):
                    state['sn'][j][k] += updated_bi[i][j]*new_z[i][k]-bi[i][j]*state['z'][i][k]
        state['z'], state['x'] = new_z, new_x
        assert state == expected
        c = scaled(updated_ci, Q(1, DEN))
        b, u, v = updated_b, scaled(updated_ui, Q(1, USCALE)), scaled(updated_vi, Q(1, VSCALE))
        assert norm_inf(c) <= K and positive_definite(c)
        assert all(sum(x*x for x in row) <= M2 for row in u)
        hmat = scaled(state['hn'], Q(1, DEN*USCALE**2))
        assert positive_definite(plus(hmat, scaled(identity(r), -c_lower*n)))
        h = scaled(multiply(inverse(state['hn']), state['tn']), Q(USCALE, VSCALE))
        z = scaled(state['z'], Q(1, DEN*VSCALE))
        x = scaled(state['x'], Q(1, DEN*USCALE))
        coefficients = plus(z, scaled(multiply(x, h), -1))
        features = [left+right for left, right in zip(b, u)]
        full_gram = plus(updated_g, multiply(u, transpose(u)))
        # A row with off-diagonal sum >= diagonal forbids every positive
        # scalar alpha from making ||I-alpha*G||_infinity strictly below one.
        assert any(sum(full_gram[i][j] for j in range(n) if j != i) >= full_gram[i][i]
                   for i in range(n) if updated_active[i])
        direct_state = multiply(transpose(features), coefficients)
        s0 = scaled(state['sn'], Q(1, BSCALE*DEN*VSCALE))
        y = scaled(state['yn'], Q(1, BSCALE*DEN*USCALE))
        implicit_state = plus(s0, scaled(multiply(y, h), -1))+h
        assert implicit_state == direct_state
        reads = multiply(features, direct_state)
        # A coordinate maximum, independently of the matrix row-sum norm.
        error = max(abs(reads[i][j]-v[i][j]) for i in range(n) for j in range(value_width))
        residual = plus(multiply(updated_g, c), scaled(identity(n), -1))
        identity_error = multiply(residual, plus(v, scaled(multiply(u, h), -1)))
        assert plus(reads, scaled(v, -1)) == identity_error
        assert error <= bound
        probe = next(i for i, on in enumerate(updated_active) if on)
        query = features[probe][:]
        query[(probe+2) % n] += Q(1, 1024)
        query[n] += Q(1, 1024)
        noisy = multiply([query], direct_state)[0]
        # h_incidence<=3, |B_ij|<=1, K<=2; use M<=2 for the u term.
        lipschitz = 3*K*gamma/Q(1024)+2*K/c_lower/Q(1024)
        assert all(abs(x-w) <= lipschitz for x, w in zip(noisy, reads[probe]))
        assert all(abs(x-w) <= bound+lipschitz for x, w in zip(noisy, v[probe]))
        max_error = max(max_error, error)
        for i in range(n):
            for row in (state['z'][i], state['x'][i]):
                assert all(Q(x).denominator == 1 for x in row)
            assert all(abs(x) <= DEN*VSCALE*K for x in state['z'][i])
            assert all(abs(x) <= DEN*USCALE*K*2 for x in state['x'][i])
        assert max(abs(x) for row in state['hn'] for x in row) <= (1+n*M2*K)*DEN*USCALE**2
        assert max(abs(x) for row in state['tn'] for x in row) <= n*2*K*DEN*USCALE*VSCALE
        if n == 9 and r == 1 and depth == 2 and step == steps-1:
            # Independent small dense inverse oracle for the Woodbury solve.
            precision = plus(inverse(c), multiply(u, transpose(u)))
            assert coefficients == multiply(inverse(precision), v)
        counts['replacements' if active[target] == new_active else
               'insertions' if new_active else 'deletions'] += 1
        max_ball, max_changed = max(max_ball, len(ball)), max(max_changed, len(changed))
        bi, ui, vi, active, ci, g = updated_bi, updated_ui, updated_vi, updated_active, updated_ci, updated_g
    for _ in range(4):
        target = rng.choice([i for i, on in enumerate(active) if on])
        new_value = [rng.randrange(-4, 5) for _ in range(value_width)]
        delta = [x-y for x, y in zip(new_value, vi[target])]
        # Cheap value-only path: H, X and Y are unchanged.
        for i in range(n):
            if ci[i][target]:
                for k in range(value_width):
                    change = ci[i][target]*delta[k]
                    state['z'][i][k] += change
                    for j, coefficient in enumerate(bi[i]):
                        if coefficient:
                            state['sn'][j][k] += coefficient*change
        for j in range(r):
            for k in range(value_width):
                state['tn'][j][k] += state['x'][target][j]*delta[k]
        vi[target] = new_value
        assert state == state_from_integers(bi, ui, vi, ci)
        h = scaled(multiply(inverse(state['hn']), state['tn']), Q(USCALE, VSCALE))
        coefficients = plus(scaled(state['z'], Q(1, DEN*VSCALE)),
                            scaled(multiply(scaled(state['x'], Q(1, DEN*USCALE)), h), -1))
        direct_state = multiply(transpose(features), coefficients)
        reads = multiply(features, direct_state)
        v = scaled(vi, Q(1, VSCALE))
        error = max(abs(reads[i][j]-v[i][j]) for i in range(n) for j in range(value_width))
        assert error <= bound
        max_error = max(max_error, error)
    if depth == 4:
        assert bound < Q(1, 10)
    return {"slot_capacity": n, "shared_feature_dimension": r, "value_width": value_width,
            "polynomial_depth": depth, "coefficient_fractional_bits": BITS,
            "uniform_local_contraction_bound": str(RHO), "uniform_rounded_inverse_norm_bound": str(K),
            "uniform_shared_row_norm_squared_bound": str(M2),
            "shared_matrix_eigenvalue_lower_bound_divided_by_capacity": str(c_lower),
            "history_independent_maximum_error_bound": str(bound),
            "maximum_observed_exact_error": str(max_error), "operations": steps,
            "value_only_overwrites": 4,
            "operation_counts": counts, "maximum_union_graph_ball": max_ball,
            "maximum_changed_polynomial_columns": max_changed,
            "full_gram_scalar_infinity_contraction_impossible_in_controls": True,
            "current_integer_cache_identities_verified": True,
            "not_a_runtime_or_neural_benchmark": True}


def main():
    fixtures = [run_case(n, r, depth) for n, r in ((9, 1), (9, 3), (27, 3))
                for depth in (0, 2, 4)]
    summary = {"cases": len(fixtures), "local_and_shared_key_edits": sum(x['operations'] for x in fixtures),
               "value_only_overwrites": sum(x['value_only_overwrites'] for x in fixtures),
               "shared_dimensions_checked": sorted({x['shared_feature_dimension'] for x in fixtures}),
               "depth_four_uniform_bounds_below_one_tenth": sum(x['polynomial_depth'] == 4 for x in fixtures)}
    sources = ["scripts/run_shared_feature_retention_checks.py", "scripts/run_retention_tradeoff_checks.py",
               "docs/interleaving_research.md"]
    report = {"schema_version": 1, "command": "python3 -B scripts/run_shared_feature_retention_checks.py",
              "arithmetic": "fractions.Fraction and exact bounded canonical integer caches",
              "scope": "Conditional encoded-feature theorem and classical Woodbury/polynomial controls. No strongest-baseline advantage, end-to-end semantic encoder or novelty established.",
              "fixtures": fixtures, "summary": summary,
              "new_publishable_advantage_established": False,
              "source_sha256": {f: sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/"results/research/shared_feature_retention_checks.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
