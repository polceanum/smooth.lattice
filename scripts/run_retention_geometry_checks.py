#!/usr/bin/env python3
"""Exact controls for address geometry and signed quadratic overwrite.

These support the second proof iteration, not a novelty or ML performance claim.
No training, floating-point tolerances, external packages or network calls.
"""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import random

from run_retention_tradeoff_checks import (
    apply, identity, multiply, norm_inf, plus, rank, scaled, transpose,
)

ROOT = Path(__file__).resolve().parents[1]


def dot(x, y):
    return sum((a*b for a, b in zip(x, y)), Q(0))


def trace_product(a, b):
    return sum((a[i][j]*b[j][i] for i in range(len(a))
                for j in range(len(a))), Q(0))


def projector(x):
    length2 = dot(x, x)
    return [[a*b/length2 for b in x] for a in x]


def gram(rows):
    return multiply(rows, transpose(rows))


def address_breadth_controls():
    positive, signed = [], []
    for n in (3, 5, 9):
        for t in (Q(1, 8), Q(1, 4), Q(1, 2)):
            c, s = (1-t*t)/(1+t*t), 2*t/(1+t*t)
            f = scaled(identity(n), c)
            for i in range(n):
                f[i][(i+1) % n] = s
            g = gram(f)
            assert all(g[i][i] == 1 for i in range(n))
            breadth = [sum(row)**2 for row in f]
            overlap = max(sum(row)-1 for row in g)
            # Cauchy bound with the actual row l1 norms, before squaring.
            assert sum(map(sum, f))**2 <= n*sum(map(sum, g))
            assert all(x == breadth[0] for x in breadth)
            assert breadth[0] == 1+overlap and overlap < 1
            assert norm_inf(plus(identity(n), scaled(g, -1))) == overlap
            positive.append({"keys_and_cells": n, "parameter": str(t),
                             "effective_address_width": str(breadth[0]),
                             "off_diagonal_row_sum": str(overlap),
                             "breadth_bound_attained": True})
    h = [[Q(1)]]
    while len(h) < 16:
        h = [row+row for row in h] + [row+[-x for x in row] for row in h]
        n = len(h)
        if n not in (4, 16):
            continue
        f = scaled(h, Q(1, 2 if n == 4 else 4))
        assert gram(f) == identity(n)
        breadth = sum(abs(x) for x in f[0])**2
        assert breadth == n
        assert n < n*breadth  # Positivity cannot be dropped from the bound.
        signed.append({"keys_and_cells": n, "effective_address_width": str(breadth),
                       "gram_is_identity": True,
                       "dense_signed_change_of_basis_not_advantage": True})
    return positive, signed


def additive_psd_controls():
    fixtures = []
    for p in (2, 3, 4):
        for n in sorted({p+1, 2*p, p*p}):
            points = [[Q((i+1)**j) for j in range(p)] for i in range(n)]
            atoms = [projector(x) for x in points]
            g = [[trace_product(a, b) for b in atoms] for a in atoms]
            assert all(g[i][i] == 1 for i in range(n))
            assert all(x >= 0 for row in g for x in row)
            total = [[sum(a[i][j] for a in atoms) for j in range(p)]
                     for i in range(p)]
            assert sum(map(sum, g)) == trace_product(total, total)
            assert sum(map(sum, g)) >= Q(n*n, p)
            overlaps = [sum(row)-1 for row in g]
            assert max(overlaps) >= Q(n, p)-1
            worst = Q(0)
            for target in range(n):
                labels = [Q(1 if j == target else -1) for j in range(n)]
                answer = dot(g[target], labels)/sum(g[target])
                error = abs(answer-1)
                assert error == 2*overlaps[target]/(1+overlaps[target])
                worst = max(worst, error)
                # Exhaust every label extreme in the small controls.
                if n <= 6:
                    for mask in range(1 << n):
                        v = [Q(1 if (mask >> j) & 1 else -1) for j in range(n)]
                        other_error = abs(dot(g[target], v)/sum(g[target])-v[target])
                        assert other_error <= error
            lower = 2*max(Q(0), 1-Q(p, n))
            assert worst >= lower
            fixtures.append({"key_dimension": p, "keys": n,
                             "worst_label_error": str(worst),
                             "uniform_lower_bound": str(lower),
                             "scope": "canonical additive normalized rank-one PSD memory only"})
    return fixtures


def real_mub_bases(p):
    if p == 2:
        # Projectors avoid the irrational 1/sqrt(2) normalization.
        vectors = [identity(2), [[Q(1), Q(1)], [Q(1), Q(-1)]]]
    elif p == 4:
        h = [[Q(1)]]
        while len(h) < p:
            h = [row+row for row in h] + [row+[-x for x in row] for row in h]
        # A sign flip in one coordinate supplies the third real MUB.
        vectors = [identity(p), h, [[-row[0]]+row[1:] for row in h]]
    else:
        raise ValueError("Exact controls use p=2 or p=4, not an existence search")
    return [[projector(row) for row in basis] for basis in vectors]


def dual_write_controls():
    fixtures = []
    for p in (2, 4):
        bases = real_mub_bases(p)
        for a, first in enumerate(bases):
            for b, second in enumerate(bases):
                for i, x in enumerate(first):
                    for j, y in enumerate(second):
                        assert trace_product(x, y) == (Q(i == j) if a == b else Q(1, p))
        reads, writes = [], []
        for basis in bases:
            for atom in basis[:-1]:
                reads.append(atom)
                writes.append(plus(atom, scaled(basis[-1], -1)))
        n = len(reads)
        assert n == len(bases)*(p-1)
        assert n == p*(p+1)//2-1  # Both controls span real traceless matrices.
        assert [[trace_product(a, b) for b in writes] for a in reads] == identity(n)
        assert all(rank(w) == 2 and sum(w[i][i] for i in range(p)) == 0 for w in writes)
        second_moment, state_moment = Q(0), Q(0)
        for mask in range(1 << n):
            signs = [1 if (mask >> j) & 1 else -1 for j in range(n)]
            x = [[sum(signs[k]*reads[k][i][j] for k in range(n))
                  for j in range(p)] for i in range(p)]
            s = [[sum(signs[k]*writes[k][i][j] for k in range(n))
                  for j in range(p)] for i in range(p)]
            x2, s2 = trace_product(x, x), trace_product(s, s)
            assert trace_product(s, x) == n
            assert n*n <= x2*s2
            second_moment += x2
            state_moment += s2
        assert second_moment/(1 << n) == n
        assert state_moment/(1 << n) == 2*n
        denominator = p
        integer_writes = [scaled(w, denominator) for w in writes]
        assert all(x.denominator == 1 for w in integer_writes for row in w for x in row)
        state, labels = scaled(identity(p), 0), [0]*n
        rng = random.Random(7541+p)
        steps, value_cap = 128, 7
        maximum_entry = 0
        for step in range(steps):
            target = step % n if step < n else rng.randrange(n)
            value = rng.randrange(-value_cap, value_cap+1)
            current = trace_product(reads[target], state)/denominator
            assert current == labels[target]
            state = plus(state, scaled(integer_writes[target], value-current))
            labels[target] = value
            expected = [[sum(labels[k]*integer_writes[k][i][j] for k in range(n))
                         for j in range(p)] for i in range(p)]
            assert state == expected
            # Independently restore the full oracle from the claimed packed state.
            coordinates = [(i, j) for i in range(p) for j in range(i, p)
                           if (i, j) != (p-1, p-1)]
            packed = [state[i][j] for i, j in coordinates]
            assert len(packed) == n
            restored = scaled(identity(p), 0)
            for (i, j), coefficient in zip(coordinates, packed):
                restored[i][j] = restored[j][i] = coefficient
            restored[-1][-1] = -sum(restored[i][i] for i in range(p-1))
            assert restored == state
            assert [trace_product(a, state)/denominator for a in reads] == labels
            assert all(x.denominator == 1 for row in state for x in row)
            maximum_entry = max(maximum_entry, max(abs(x) for row in state for x in row))
            assert maximum_entry <= n*denominator*value_cap
        # A full row of +1 labels exposes the reference-direction noise issue.
        reference = bases[0][-1]
        fragile = plus(identity(p), scaled(reference, -p))
        assert all(trace_product(a, fragile) == 1 for a in bases[0][:-1])
        assert all(trace_product(a, fragile) == 0 for basis in bases[1:] for a in basis)
        # q=(3/5)e_0+(4/5)e_ref is exactly a unit vector.
        q = [Q(0)]*p
        q[0], q[-1] = Q(3, 5), Q(4, 5)
        assert dot(q, q) == 1
        noisy_read = dot(q, apply(fragile, q))
        assert noisy_read == 1-p*Q(16, 25)
        # The ordinary delta direction immediately changes a cross-basis read.
        delta_cross_interference = trace_product(bases[1][0], bases[0][0])
        assert delta_cross_interference == Q(1, p)
        fixtures.append({"real_key_dimension": p, "bases": len(bases), "associations": n,
                         "packed_traceless_state_cells_per_value": n,
                         "write_matrix_rank": 2, "exact_stream_overwrites": steps,
                         "rademacher_assignments_checked": 1 << n,
                         "integer_state_scale": denominator, "value_cap": value_cap,
                         "maximum_integer_entry_observed": str(maximum_entry),
                         "integer_entry_bound": n*denominator*value_cap,
                         "ordinary_delta_cross_read_change_for_unit_correction": str(delta_cross_interference),
                         "reference_mixture_read": str(noisy_read),
                         "fixed_codebook_and_algebraic_queries_only": True,
                         "known_dual_frame_control_not_novel_architecture": True})
    return fixtures


def main():
    positive, signed = address_breadth_controls()
    additive = additive_psd_controls()
    dual = dual_write_controls()
    sources = ["scripts/run_retention_geometry_checks.py",
               "scripts/run_retention_tradeoff_checks.py", "docs/interleaving_research.md"]
    summary = {"positive_breadth_equality_cases": len(positive),
               "signed_change_of_basis_controls": len(signed),
               "additive_psd_capacity_cases": len(additive),
               "dual_frame_cases": len(dual),
               "exact_dual_write_overwrites": sum(x["exact_stream_overwrites"] for x in dual),
               "rademacher_assignments_checked": sum(x["rademacher_assignments_checked"] for x in dual)}
    report = {"schema_version": 1,
              "command": "python3 -B scripts/run_retention_geometry_checks.py",
              "arithmetic": "fractions.Fraction; bounded exact integer states in overwrite controls",
              "scope": "Supporting geometry obstructions and known dual-frame controls, not a new theorem beyond prior art or evidence of ML advantage.",
              "positive_address_breadth": positive,
              "signed_basis_escape": signed,
              "normalized_additive_psd_capacity": additive,
              "signed_quadratic_dual_writes": dual,
              "summary": summary,
              "uniform_proof_location": "docs/interleaving_research.md, Associative-memory goal: address geometry and signed writes",
              "new_publishable_advantage_established": False,
              "source_sha256": {f: sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/"results/research/retention_geometry_checks.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
