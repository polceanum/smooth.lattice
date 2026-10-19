#!/usr/bin/env python3
"""Exact controls for common-grid additive heads and correlated lookups.

Counts unfused head-support entries, not an implementation-independent cost.
Layout, fusion and indexed-lookup controls prevent an overstated separation.
No training, external dependencies or network access.
"""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
import random

ROOT = Path(__file__).resolve().parents[1]


def zeros(n):
    return [[Q(0)]*n for _ in range(n)]


def target(n, perm=None):
    perm = list(range(n)) if perm is None else perm
    out = zeros(n)
    for i, j in enumerate(perm):
        out[i][j] = Q(1, n)
    return out


def head(a, b, k):
    """Algebraic top-k softmax equivalent: positive separable raw weights."""
    n = len(a)
    entries = sorted(((a[i]*b[j], i, j) for i in range(n) for j in range(n)),
                     key=lambda x: (-x[0], x[1], x[2]))[:k]
    total = sum(x[0] for x in entries)
    out = zeros(n)
    for value, i, j in entries:
        out[i][j] = value/total
    return out, total


def support(p):
    return {(i, j) for i, row in enumerate(p) for j, x in enumerate(row) if x}


def combine(terms):
    n = len(terms[0][1])
    return [[sum((c*p[i][j] for c, p in terms), Q(0))
             for j in range(n)] for i in range(n)]


def l1_difference(a, b):
    return sum(abs(x-y) for u, v in zip(a, b) for x, y in zip(u, v))


def check_support_bound(terms):
    n = len(terms[0][1])
    positive = [(c, p) for c, p in terms if c > 0]
    negative_support = set().union(*(support(p) for c, p in terms if c < 0))
    hp = len(positive)
    assert hp
    w = combine(terms)
    error = l1_difference(w, target(n))
    good = [i for i in range(n) if w[i][i] >= Q(1, 2*n)]
    assigned = [[] for _ in positive]
    for i in good:
        best = max(range(hp), key=lambda h: positive[h][0]*positive[h][1][i][i])
        assert positive[best][0]*positive[best][1][i][i] >= Q(1, 2*n*hp)
        assigned[best].append(i)
    forced = set()
    for (c, p), group in zip(positive, assigned):
        for i, j in combinations(group, 2):
            threshold = min(p[i][i], p[j][j])
            candidates = [(i, j), (j, i)]
            edge = next((e for e in candidates if p[e[0]][e[1]] >= threshold), None)
            assert edge is not None
            assert c*p[edge[0]][edge[1]] >= Q(1, 2*n*hp)
            assert edge not in forced
            forced.add(edge)
    missing = forced-negative_support
    assert len(missing)*Q(1, 2*n*hp) <= error
    assert len(forced) >= (Q(len(good)**2, hp)-len(good))/2
    cost = sum(len(support(p)) for c, p in terms if c)
    assert cost >= Q(len(good)**2, hp)-2*n*hp*error
    if error <= Q(1, 2):
        assert len(good) >= n*(1-2*error)
        assert cost >= n*n*(1-2*error)**2/hp-2*n*hp*error
    if not error:
        assert cost*hp >= n*n
    return cost, error


def grouped_signed(n, groups):
    assert n % groups == 0
    width = n//groups
    terms = []
    for block in range(groups):
        a, b = [Q(1, 2**(4*n))]*n, [Q(1, 2**(4*n))]*n
        for j in range(width):
            a[block*width+j], b[block*width+j] = Q(1, 2**j), Q(2**j)
        upper, z_upper = head(a, b, width*(width+1)//2)
        terms.append((z_upper/n, upper))
        if width > 1:
            strict, z_strict = head(a, b, width*(width-1)//2)
            terms.append((-z_strict/n, strict))
    assert combine(terms) == target(n)
    cost, error = check_support_bound(terms)
    assert cost == n*n//groups and not error
    return terms


def sparsemax(scores):
    entries = sorted((x for row in scores for x in row), reverse=True)
    prefix, tau = Q(0), None
    for k, x in enumerate(entries, 1):
        prefix += x
        candidate = (prefix-1)/k
        if x > candidate:
            tau = candidate
    out = [[max(Q(0), x-tau) for x in row] for row in scores]
    assert sum(map(sum, out)) == 1
    return out


def rounded(x):
    y = x+Q(1, 2)
    return y.numerator//y.denominator


def main():
    head_cases = 0
    for n in (2, 3):
        for a in product((Q(1), Q(2), Q(4)), repeat=n):
            for b in product((Q(1), Q(2), Q(4)), repeat=n):
                for k in range(1, n*n+1):
                    p, _ = head(a, b, k)
                    for i, j in combinations(range(n), 2):
                        if p[i][i] and p[j][j]:
                            assert max(p[i][j], p[j][i]) >= min(p[i][i], p[j][j])
                    check_support_bound([(Q(1), p)])
                    head_cases += 1
    tight, approximate = [], []
    for n in (2, 4, 6, 8, 12):
        for groups in range(1, n+1):
            if n % groups:
                continue
            terms = grouped_signed(n, groups)
            cost, _ = check_support_bound(terms)
            tight.append({"matching_size": n, "positive_heads": groups,
                          "total_heads": len(terms), "unfused_support_entries": cost,
                          "head_access_product": cost*groups,
                          "fused_nonzero_coefficients": len(support(combine(terms)))})
            noise, _ = head([Q(1)]*n, [Q(1)]*n, n*n)
            perturbed = [(Q(127, 128)*c, p) for c, p in terms]+[(Q(1, 128), noise)]
            _, error = check_support_bound(perturbed)
            approximate.append({"matching_size": n, "groups": groups,
                                "coefficient_l1_error": str(error)})
    rng = random.Random(91403)
    routing_cases = 0
    for n in range(2, 6):
        for perm in permutations(range(n)):
            for eta in (Q(0), Q(1, 16), Q(1, 8)):
                x = [Q(i+1)+rng.choice((-eta, eta)) for i in range(n)]
                y = [Q(0)]*n
                for i, j in enumerate(perm):
                    y[j] = Q(i+1)+rng.choice((-eta, eta))
                scores = [[-Q(4, n)*(a-b)**2 for b in y] for a in x]
                p = sparsemax(scores)
                expected = target(n, perm)
                assert support(p) == support(expected)
                assert l1_difference(p, expected) <= 8*eta*eta
                index = {rounded(value): j for j, value in enumerate(y)}
                assert [index[rounded(value)] for value in x] == list(perm)
                routing_cases += 1
    # Arbitrary table layout defeats an unrestricted interpretation of the bound.
    layout_controls, overwrite_checks = [], 0
    for n in (2, 4, 8):
        physical, _ = head([Q(1)]*n, [Q(2)]+[Q(1)]*(n-1), n)
        logical = [[physical[i][(j-i) % n] for j in range(n)] for i in range(n)]
        assert logical == target(n)
        terms = grouped_signed(n, 1)
        assert len(support(combine(terms))) == n
        values = [[0]*n for _ in range(n)]
        for _ in range(64):
            i, j, v = rng.randrange(n), rng.randrange(n), rng.choice((-1, 1))
            values[i][j] = v
            perm = list(range(n))
            rng.shuffle(perm)
            tags_y = [0]*n
            for row, col in enumerate(perm):
                tags_y[col] = row+1
            inverse_tags = {tag: col for col, tag in enumerate(tags_y)}
            lookup = sum(Q(values[row][inverse_tags[row+1]], n) for row in range(n))
            # Each uniformly sampled column directly identifies its row because
            # the first tag vector is the promised canonical ordering.
            column_expectation = sum(Q(values[tags_y[col]-1][col], n) for col in range(n))
            dense = sum(expected*values[row][col] for row, line in enumerate(target(n, perm))
                        for col, expected in enumerate(line))
            assert lookup == dense == column_expectation
            overwrite_checks += 1
        layout_controls.append({"matching_size": n, "logical_target": "uniform diagonal",
                                "layout": "(i,j) -> (i,(j-i) mod n)",
                                "physical_additive_heads": 1, "physical_active_values": n,
                                "signed_unfused_entries": n*n, "fused_active_values": n,
                                "unrestricted_probe_lower_bound_refuted": True})
    sources = ["scripts/run_correlated_routing_checks.py", "docs/interleaving_research.md"]
    summary = {"additive_head_cases_including_ties": head_cases,
               "tight_signed_decompositions": len(tight), "approximate_mixture_cases": len(approximate),
               "correlated_raw_query_cases": routing_cases,
               "layout_and_fusion_escape_cases": len(layout_controls),
               "dictionary_point_overwrites_and_matching_reads": overwrite_checks}
    report = {"schema_version": 1, "command": "python3 -B scripts/run_correlated_routing_checks.py",
              "arithmetic": "fractions.Fraction",
              "scope": "Common-grid, unfused head-support tradeoff and explicit stronger-baseline escapes. Not an implementation-independent probe lower bound or a new strongest-method advantage.",
              "tight_signed_decompositions": tight, "approximate_mixtures": approximate,
              "layout_and_fusion_escapes": layout_controls, "summary": summary,
              "novelty_established": False, "new_publishable_advantage_established": False,
              "source_sha256": {f: sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/"results/research/correlated_routing_checks.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
