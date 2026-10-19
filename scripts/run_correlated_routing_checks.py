#!/usr/bin/env python3
"""Exact controls for common-grid additive heads and correlated lookups.

Counts unfused head-support entries, not an implementation-independent cost.
Layout, fusion and indexed-lookup controls prevent an overstated separation.
No training, external dependencies or network access.
"""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, permutations, product
from bisect import bisect_right
from math import comb
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


def prefix_join_setup(a, b, lengths):
    """Integer query weights; no value-dependent preprocessing or real oracle."""
    assert len(a) == len(lengths)
    assert all(isinstance(x, int) and x >= 0 for x in (*a, *b))
    assert all(0 <= k <= len(b) for k in lengths)
    column_cdf, running = [], 0
    for weight in b:
        running += weight
        column_cdf.append(running)
    row_cdf, running = [], 0
    for weight, length in zip(a, lengths):
        running += weight*(column_cdf[length-1] if length else 0)
        row_cdf.append(running)
    return row_cdf, column_cdf


def prefix_join_sample(row_cdf, column_cdf, lengths, row_ticket, column_ticket):
    """Tickets are fresh uniform integers in their respective CDF ranges."""
    assert 0 <= row_ticket < row_cdf[-1]
    i = bisect_right(row_cdf, row_ticket)
    length = lengths[i]
    assert 0 <= column_ticket < column_cdf[length-1]
    j = bisect_right(column_cdf, column_ticket, hi=length)
    return i, j


def explicit_join_oracle(a, b, lengths):
    """Independent enumeration of join completions, for certification only.

    A production indexed-join baseline replaces these sums by prefix sums.
    This deliberately slow oracle is not a runtime baseline or benchmark.
    """
    weights = {(i, j): a[i]*b[j]
               for i, length in enumerate(lengths) for j in range(length)}
    total = sum(weights.values())
    if not total:
        return total, {}
    return total, {edge: Q(weight, total) for edge, weight in weights.items()}


def prefix_join_distribution(a, b, lengths):
    row_cdf, column_cdf = prefix_join_setup(a, b, lengths)
    total = row_cdf[-1]
    assert total
    result = {}
    previous = 0
    for i, row_total in enumerate(row_cdf):
        row_mass = row_total-previous
        previous = row_total
        for j in range(lengths[i]):
            conditional = Q(b[j], column_cdf[lengths[i]-1]) if row_mass else Q(0)
            result[i, j] = Q(row_mass, total)*conditional
    return result


def check_weighted_prefix_join():
    cases = zero_mass = row_tickets = child_tickets = 0
    # Zero weights, empty rows/columns, rectangular grids and every monotone
    # prefix mask of these sizes. Exhaustion concerns finite controls only.
    for m, n in ((1, 1), (2, 2), (2, 3), (3, 3)):
        for lengths in combinations_with_replacement(range(n+1), m):
            for a in product(range(3), repeat=m):
                for b in product(range(3), repeat=n):
                    total, oracle = explicit_join_oracle(a, b, lengths)
                    rows, columns = prefix_join_setup(a, b, lengths)
                    assert rows[-1] == total
                    if not total:
                        zero_mass += 1
                        continue
                    assert prefix_join_distribution(a, b, lengths) == oracle
                    # Exhaust all row tickets and all conditional child tickets
                    # independently; do not confuse the two random ranges.
                    counts = [0]*m
                    for ticket in range(total):
                        i = bisect_right(rows, ticket)
                        counts[i] += 1
                        assert ticket < sum(sum(a[r]*b[j] for j in range(lengths[r]))
                                            for r in range(i+1))
                        row_tickets += 1
                    for i, count in enumerate(counts):
                        assert count == sum(a[i]*b[j] for j in range(lengths[i]))
                        if not count:
                            continue
                        child_counts = [0]*n
                        row_ticket = rows[i]-count
                        for ticket in range(columns[lengths[i]-1]):
                            ii, j = prefix_join_sample(rows, columns, lengths, row_ticket, ticket)
                            assert ii == i and j < lengths[i]
                            child_counts[j] += 1
                            child_tickets += 1
                        for j in range(lengths[i]):
                            assert child_counts[j] == b[j]
                            assert Q(count, total)*Q(child_counts[j], sum(child_counts)) == oracle[i, j]
                    cases += 1

    rng = random.Random(720418)
    overwrites = point_reads = fresh_histories = 0
    for lengths in ((2, 2), (1, 2, 3), (0, 2, 3), (1, 1, 2), (0, 1, 3, 3)):
        m, n = len(lengths), max(lengths)
        edges = [(i, j) for i, length in enumerate(lengths) for j in range(length)]
        offsets = [0]
        for length in lengths:
            offsets.append(offsets[-1]+length)
        values = {edge: (Q(0), Q(0), Q(0)) for edge in edges}
        packed = [values[edge] for edge in edges]
        for _ in range(128):
            edge = rng.choice(edges)
            value = tuple(Q(rng.randrange(-8, 9), 8) for _ in range(3))
            values[edge] = value
            packed[offsets[edge[0]]+edge[1]] = value
            a, b = [rng.randrange(3) for _ in range(m)], [rng.randrange(3) for _ in range(n)]
            # Ensure a valid query, retaining other zero-weight controls above.
            a[edge[0]], b[edge[1]] = 1, 1
            total, oracle = explicit_join_oracle(a, b, lengths)
            distribution = prefix_join_distribution(a, b, lengths)
            for d in range(3):
                dense = sum(oracle[e]*values[e][d] for e in edges)
                conditional_mean = sum(distribution[i, j]*packed[offsets[i]+j][d]
                                       for i, j in edges)
                assert dense == conditional_mean
            ea, eb = [0]*m, [0]*n
            ea[edge[0]], eb[edge[1]] = 1, 1
            exact = prefix_join_distribution(ea, eb, lengths)
            assert exact[edge] == 1 and sum(exact.values()) == 1
            assert packed[offsets[edge[0]]+edge[1]] == values[edge]
            overwrites += 1
            point_reads += 1
        # Current payloads may depend on the preceding sample. Conditional on
        # that history, fresh sampling still has exactly the current mean.
        uniform = prefix_join_distribution([1]*m, [1]*n, lengths)
        for preceding in edges:
            state = {edge: Q(int(edge == preceding)) for edge in edges}
            assert sum(uniform[e]*state[e] for e in edges) == Q(1, len(edges))
            assert state[preceding] == 1  # Reuse would instead return one.
            fresh_histories += 1

    independent = prefix_join_distribution([1, 1], [1, 1], (1, 2))
    row_marginal = [sum(p for (r, _), p in independent.items() if r == i) for i in range(2)]
    col_marginal = [sum(p for (_, c), p in independent.items() if c == j) for j in range(2)]
    assert row_marginal[0]*col_marginal[1] == Q(1, 9)  # Invalid pair.
    rejection_controls = []
    for bits in (2, 8, 32, 128):
        weight = 2**bits
        total, _ = explicit_join_oracle([weight, 1], [1, weight], (1, 2))
        acceptance = Q(total, (weight+1)**2)
        assert acceptance == Q(2*weight+1, (weight+1)**2)
        assert 1/acceptance >= Q(weight, 3)
        rejection_controls.append({"weight_bits": bits+1,
                                   "proposal_acceptance": str(acceptance),
                                   "join_tuple_rejections": 0})

    noise_controls, conditioning_controls = 0, []
    for m, n, lengths in ((2, 2, (2, 2)), (3, 3, (1, 2, 3)), (3, 2, (0, 1, 2))):
        for i, length in enumerate(lengths):
            for j in range(length):
                for ticks in (0, 1, 2, 4):
                    t = Q(ticks, 8)
                    a = [(1-t if r == i else t/Q(m-1)) for r in range(m)]
                    b = [(1-t if c == j else t/Q(n-1)) for c in range(n)]
                    # Rational weights can be scaled independently to integers.
                    aa, bb = [int(x*16) for x in a], [int(x*16) for x in b]
                    p = prefix_join_distribution(aa, bb, lengths)
                    assert p[i, j] >= (1-t)**2
                    assert 2*(1-p[i, j]) <= 2*(2*t-t*t)
                    noise_controls += 1
    for bits in (2, 8, 32, 128):
        weight = 2**bits
        first = prefix_join_distribution([1, 0], [1, weight-1], (1, 2))
        second = prefix_join_distribution([weight-1, 1], [0, 1], (1, 2))
        distance = Q(4, weight)  # Combined simplex l1 distance.
        assert first[0, 0] == second[1, 1] == 1
        assert sum(abs(first[e]-second[e]) for e in first) == 2
        conditioning_controls.append({"t": str(Q(1, weight)),
                                      "query_l1_distance": str(distance),
                                      "output_gap": "2", "normalized_mask_mass": str(Q(1, weight))})

    # A fixed sample's guarantee cannot hold for all later payload assignments.
    reused_samples = [(i, i) for i in range(4)]
    adversarial_values = {(i, j): Q(-1 if (i, j) in reused_samples else 1)
                          for i in range(4) for j in range(4)}
    exact_mean = sum(adversarial_values.values())/16
    reused_mean = sum(adversarial_values[e] for e in reused_samples)/4
    assert exact_mean == Q(1, 2) and reused_mean == -1

    # Exact Bernoulli tail controls. The rational expression is an upper bound
    # on the classical exponential bound, not a floating-point approximation.
    tails = 0
    for k in range(1, 25):
        for epsilon in (Q(1, 4), Q(1, 2), Q(3, 4)):
            probability = Q(sum(comb(k, j) for j in range(k+1)
                                if abs(Q(2*j, k)-1) > epsilon), 2**k)
            upper = 2/(1+Q(k)*epsilon**2/128)**64
            assert probability <= upper
            tails += 1
    return {"positive_mass_query_cases": cases, "zero_mass_queries_rejected": zero_mass,
            "exact_row_tickets": row_tickets, "exact_conditional_tickets": child_tickets,
            "canonical_vector_overwrites": overwrites, "point_code_exact_reads": point_reads,
            "adaptive_fresh_sample_histories": fresh_histories,
            "point_query_simplex_noise_controls": noise_controls,
            "independent_marginals_invalid_pair_mass": "1/9",
            "rejection_controls": rejection_controls,
            "vanishing_mass_conditioning_controls": conditioning_controls,
            "reused_sample_counterexample": {"samples": reused_samples, "true_mean": str(exact_mean),
                                             "reused_mean": str(reused_mean), "error": str(exact_mean-reused_mean)},
            "exact_Bernoulli_tail_controls": tails,
            "baseline_verdict": "Identical conditional distribution, payload table, prefix metadata and update semantics as the weighted indexed-join specialization; no new strongest-baseline advantage.",
            "novelty_established": False, "new_publishable_advantage_established": False}


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
    weighted_prefix_join = check_weighted_prefix_join()
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
              "weighted_prefix_join": weighted_prefix_join,
              "novelty_established": False, "new_publishable_advantage_established": False,
              "source_sha256": {f: sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/"results/research/correlated_routing_checks.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(summary, sort_keys=True))
    print(json.dumps(weighted_prefix_join, sort_keys=True))


if __name__ == "__main__":
    main()
