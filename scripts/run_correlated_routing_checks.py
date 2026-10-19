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
from math import comb, gcd
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


def gram_dot(a, b):
    return sum((x*y for x, y in zip(a, b)), Q(0))


def gram_matrix(rows):
    d = len(rows[0])
    return [[sum((row[i]*row[j] for row in rows), Q(0)) for j in range(d)] for i in range(d)]


def gram_quadratic(matrix, q):
    return sum((q[i]*matrix[i][j]*q[j] for i in range(len(q)) for j in range(len(q))), Q(0))


def gram_integer_scale(values):
    denominator = 1
    for value in values:
        denominator = denominator*Q(value).denominator//gcd(denominator, Q(value).denominator)
    return [int(value*denominator) for value in values]


def gram_psd_components(matrix):
    """Rational LDL rank-one terms, including zero pivots; no square roots."""
    d = len(matrix)
    residual = [list(map(Q, row)) for row in matrix]
    components = []
    for k in range(d):
        pivot = residual[k][k]
        assert pivot >= 0
        if not pivot:
            assert all(residual[i][k] == 0 for i in range(d))
            continue
        direction = [Q(0)]*d
        direction[k] = Q(1)
        for i in range(k+1, d):
            direction[i] = residual[i][k]/pivot
        components.append((pivot, direction))
        for i in range(k, d):
            for j in range(k, d):
                residual[i][j] -= pivot*direction[i]*direction[j]
    assert all(not value for row in residual for value in row)
    assert [[sum((s*u[i]*u[j] for s, u in components), Q(0))
             for j in range(d)] for i in range(d)] == matrix
    return components


class BlockGramIndex:
    """Exact finite certificate for the published blocked partial-Gram method."""

    def __init__(self, rows, block=None):
        self.rows = [list(map(Q, row)) for row in rows]
        self.d = len(rows[0])
        self.block = self.d if block is None else block
        self.nodes, self.leaves = [], []
        chunks = [(s, min(s+self.block, len(rows))) for s in range(0, len(rows), self.block)]

        def build(lo, hi, parent):
            index = len(self.nodes)
            start, stop = chunks[lo][0], chunks[hi-1][1]
            node = {'start': start, 'stop': stop, 'parent': parent}
            self.nodes.append(node)
            if hi-lo == 1:
                self.leaves.append(index)
                node['gram'] = gram_matrix(self.rows[start:stop])
            else:
                mid = (lo+hi)//2
                node['children'] = build(lo, mid, index), build(mid, hi, index)
                left, right = (self.nodes[c]['gram'] for c in node['children'])
                node['gram'] = [[left[i][j]+right[i][j] for j in range(self.d)] for i in range(self.d)]
            return index

        self.root = build(0, len(chunks), -1)
        assert len(self.nodes) == 2*len(chunks)-1

    def update_entry(self, row, column, value):
        old = self.rows[row][:]
        delta = Q(value)-old[column]
        node_id, writes = self.leaves[row//self.block], 0
        while node_id != -1:
            node = self.nodes[node_id]
            matrix = node['gram']
            for k in range(self.d):
                matrix[column][k] += delta*old[k]
                matrix[k][column] += delta*old[k]
            matrix[column][column] += delta*delta
            writes += 2*self.d+1
            node_id = node['parent']
        self.rows[row][column] = Q(value)
        return writes

    def validate(self):
        for node in self.nodes:
            assert node['gram'] == gram_matrix(self.rows[node['start']:node['stop']])

    def distribution(self, q):
        total = gram_quadratic(self.nodes[self.root]['gram'], q)
        if not total:
            return {}
        assert total > 0
        result = {i: Q(0) for i in range(len(self.rows))}

        def descend(index, probability):
            node = self.nodes[index]
            mass = gram_quadratic(node['gram'], q)
            assert mass > 0
            if 'children' in node:
                child_masses = [gram_quadratic(self.nodes[c]['gram'], q) for c in node['children']]
                assert sum(child_masses) == mass and min(child_masses) >= 0
                for child, child_mass in zip(node['children'], child_masses):
                    if child_mass:
                        descend(child, probability*child_mass/mass)
            else:
                weights = [gram_dot(self.rows[i], q)**2 for i in range(node['start'], node['stop'])]
                assert sum(weights) == mass
                for i, weight in zip(range(node['start'], node['stop']), weights):
                    result[i] = probability*weight/mass

        descend(self.root, Q(1))
        assert sum(result.values()) == 1
        return result

    def sample_ticket(self, q, ticket):
        assert all(Q(x).denominator == 1 for x in q)
        assert 0 <= ticket < gram_quadratic(self.nodes[self.root]['gram'], q)
        index = self.root
        while 'children' in self.nodes[index]:
            left, right = self.nodes[index]['children']
            left_mass = gram_quadratic(self.nodes[left]['gram'], q)
            if ticket < left_mass:
                index = left
            else:
                ticket -= left_mass
                index = right
        node = self.nodes[index]
        for i in range(node['start'], node['stop']):
            weight = gram_dot(self.rows[i], q)**2
            if ticket < weight:
                return i
            ticket -= weight
        raise AssertionError('Positive-mass ticket must reach a row')


def gram_pair_distribution(a_index, b_index, x):
    """Enumerate the sampler's law for proof checks, not a runtime benchmark."""
    ga, gb = a_index.nodes[0]['gram'], b_index.nodes[0]['gram']
    components = []
    for scale, direction in gram_psd_components(gb):
        q = [gram_dot(row, direction) for row in x]
        weight = scale*gram_quadratic(ga, q)
        assert weight >= 0
        if weight:
            components.append((weight, q))
    total = sum((w for w, q in components), Q(0))
    if not total:
        return total, {}
    integer_weights = gram_integer_scale([weight for weight, q in components])
    assert [Q(w, sum(integer_weights)) for w in integer_weights] == [w/total for w, q in components]
    row_probability = [Q(0)]*len(a_index.rows)
    for weight, q in components:
        distribution = a_index.distribution(q)
        assert distribution == a_index.distribution(gram_integer_scale(q))
        for i, p in distribution.items():
            row_probability[i] += weight*p/total
    result = {(i, j): Q(0) for i in range(len(a_index.rows)) for j in range(len(b_index.rows))}
    for i, row_mass in enumerate(row_probability):
        if not row_mass:
            continue
        q = [sum((x[r][c]*a_index.rows[i][r] for r in range(a_index.d)), Q(0))
             for c in range(b_index.d)]
        for j, p in b_index.distribution(q).items():
            result[i, j] = row_mass*p
    assert sum(result.values()) == 1
    return total, result


def gram_inverse(matrix):
    n = len(matrix)
    rows = [[Q(x) for x in row]+[Q(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for c in range(n):
        pivot = next((r for r in range(c, n) if rows[r][c]), None)
        if pivot is None:
            return None
        rows[c], rows[pivot] = rows[pivot], rows[c]
        scale = rows[c][c]
        rows[c] = [x/scale for x in rows[c]]
        for r in range(n):
            if r != c:
                scale = rows[r][c]
                rows[r] = [x-scale*y for x, y in zip(rows[r], rows[c])]
    return [row[n:] for row in rows]


def gram_leverages(rows):
    for rank in range(min(len(rows), len(rows[0])), 0, -1):
        for columns in combinations(range(len(rows[0])), rank):
            reduced = [[row[c] for c in columns] for row in rows]
            inverse = gram_inverse(gram_matrix(reduced))
            if inverse is not None:
                ell = [gram_quadratic(inverse, row) for row in reduced]
                assert sum(ell) == rank and all(0 <= x <= 1 for x in ell)
                return columns, inverse, ell
    return (), [], [Q(0)]*len(rows)


def check_dynamic_gram_sampler():
    rng = random.Random(221003)
    snapshots = valid_snapshots = key_edits = payload_edits = ticket_checks = probability_checks = 0
    attainment_checks = capacity_checks = 0
    for m, n, d in ((1, 2, 1), (3, 4, 2), (5, 3, 2), (4, 5, 3), (3, 2, 3)):
        aa = BlockGramIndex([[rng.randrange(-2, 3) for _ in range(d)] for _ in range(m)])
        bb = BlockGramIndex([[rng.randrange(-2, 3) for _ in range(d)] for _ in range(n)])
        values = {(i, j): (Q(0), Q(0)) for i in range(m) for j in range(n)}
        for step in range(64):
            index = aa if step % 2 else bb
            r, c = rng.randrange(len(index.rows)), rng.randrange(d)
            index.update_entry(r, c, rng.choice((-2, -1, 0, 0, 1, 2)))
            aa.validate(); bb.validate()
            key_edits += 1
            edge = rng.randrange(m), rng.randrange(n)
            values[edge] = tuple(Q(rng.randrange(-8, 9), 8) for _ in range(2))
            payload_edits += 1
            x = [[Q(rng.randrange(-2, 3)) for _ in range(d)] for _ in range(d)]
            total, sampled = gram_pair_distribution(aa, bb, x)
            raw = {(i, j): sum((aa.rows[i][r]*x[r][c]*bb.rows[j][c]
                               for r in range(d) for c in range(d)), Q(0))**2
                   for i in range(m) for j in range(n)}
            assert total == sum(raw.values())
            snapshots += 1
            if not total:
                continue
            valid_snapshots += 1
            direct = {e: w/total for e, w in raw.items()}
            assert sampled == direct
            for k in range(2):
                assert sum(sampled[e]*values[e][k] for e in values) == sum(direct[e]*values[e][k] for e in values)
            ac, ai, al = gram_leverages(aa.rows)
            bc, bi, bl = gram_leverages(bb.rows)
            for (i, j), probability in sampled.items():
                assert probability <= al[i]*bl[j]
                probability_checks += 1
            q = [rng.randrange(-1, 2) for _ in range(d)]
            mass = gram_quadratic(aa.nodes[0]['gram'], q)
            if mass:
                counts = [0]*m
                for ticket in range(int(mass)):
                    counts[aa.sample_ticket(q, ticket)] += 1
                    ticket_checks += 1
                assert [Q(x, int(mass)) for x in counts] == list(aa.distribution(q).values())
        # Directly attain each nonzero product leverage, also for singular keys.
        ac, ai, al = gram_leverages(aa.rows)
        bc, bi, bl = gram_leverages(bb.rows)
        for i in range(m):
            for j in range(n):
                if not al[i]*bl[j]:
                    continue
                u = [gram_dot(row, [aa.rows[i][c] for c in ac]) for row in ai]
                v = [gram_dot(row, [bb.rows[j][c] for c in bc]) for row in bi]
                x = [[Q(0)]*d for _ in range(d)]
                for r, cr in enumerate(ac):
                    for c, cc in enumerate(bc):
                        x[cr][cc] = u[r]*v[c]
                _, distribution = gram_pair_distribution(aa, bb, x)
                assert distribution[i, j] == al[i]*bl[j]
                worst_values = {edge: Q(1 if edge == (i, j) else -1) for edge in distribution}
                assert 1-sum(distribution[e]*worst_values[e] for e in distribution) == 2*(1-al[i]*bl[j])
                attainment_checks += 1
        for epsilon in (Q(0), Q(1, 4), Q(1, 2), Q(1), Q(3, 2)):
            qualifying = sum(a*b >= 1-epsilon/2 for a in al for b in bl)
            assert qualifying*(1-epsilon/2) <= len(ac)*len(bc)
            capacity_checks += 1

    frames, confidence_controls = [], 0
    for bases in (1, 2, 3, 5):
        rows = []
        for k in range(bases):
            t = Q(1, k+2)
            c, s = (1-t*t)/(1+t*t), 2*t/(1+t*t)
            rows.extend(((c, s), (-s, c)))
        flat = gram_integer_scale([x for row in rows for x in row])
        rows = [flat[2*i:2*i+2] for i in range(len(rows))]
        columns, inverse, ell = gram_leverages(rows)
        assert len(set(map(tuple, rows))) == 2*bases and len(columns) == 2
        assert ell == [Q(1, bases)]*(2*bases)
        p = Q(1, bases*bases)
        for sample_count in (1, 2, 4, 8):
            for epsilon in (Q(0), Q(1, 4), Q(1, 2), Q(1)):
                success = sum((Q(comb(sample_count, hits))*p**hits*(1-p)**(sample_count-hits)
                               for hits in range(sample_count+1)
                               if 2*(1-Q(hits, sample_count)) <= epsilon), Q(0))
                assert p >= (1-epsilon/2)*success
                confidence_controls += 1
        frames.append({'orthogonal_bases_per_factor': bases, 'independent_pair_payloads': len(rows)**2,
                       'query_matrix_coordinates': 4, 'maximum_individual_read_mass': str(Q(1, bases*bases)),
                       'best_worst_payload_error': str(2*(1-Q(1, bases*bases))),
                       'max_integer_key_bits': max(abs(x).bit_length() for row in rows for x in row)})
    # Full-rank bilinear query gives a correlated joint distribution.
    eye = BlockGramIndex([[1, 0], [0, 1]])
    _, correlated = gram_pair_distribution(eye, eye, [[1, 0], [0, 1]])
    assert correlated == {(0, 0): Q(1, 2), (0, 1): Q(0), (1, 0): Q(0), (1, 1): Q(1, 2)}
    zero = BlockGramIndex([[0, 0], [0, 0]])
    assert gram_pair_distribution(zero, eye, [[1, 0], [0, 1]]) == (0, {})
    assert gram_psd_components([[0, 0], [0, 1]]) == [(Q(1), [Q(0), Q(1)])]
    return {'fresh_query_snapshots': snapshots, 'positive_normalizer_snapshots': valid_snapshots,
            'zero_normalizer_snapshots': snapshots-valid_snapshots, 'factor_coordinate_replacements': key_edits,
            'independent_vector_payload_overwrites': payload_edits, 'integer_tree_tickets': ticket_checks,
            'joint_leverage_probability_bounds': probability_checks,
            'exact_maximum_probability_witnesses': attainment_checks,
            'addressability_count_controls': capacity_checks, 'tight_frame_families': frames,
            'exact_finite_sample_confidence_controls': confidence_controls,
            'zero_pivot_and_zero_normalizer_controls': True, 'full_rank_correlated_query_control': True,
            'classification': 'Constructive baseline from published Gram-tree sampling and PSD mixtures; product-leverage address ceiling is a classical projection corollary. No new strongest-method advantage.',
            'novelty_established': False, 'new_publishable_advantage_established': False}


def simplex_routing_weights(query):
    """Known simplex/Lovasz weights, including zero weights at ties."""
    order = sorted(range(len(query)), key=lambda i: (-query[i], i))
    mask, previous, result = 0, Q(1), {}
    for i in order:
        result[mask] = previous-query[i]
        mask |= 1 << i
        previous = query[i]
    result[mask] = previous
    assert min(result.values()) >= 0 and sum(result.values()) == 1
    return result


def check_simplex_routing():
    grid_queries = threshold_tickets = marginal_controls = corner_controls = 0
    for d in range(1, 6):
        for ticks in product(range(5), repeat=d):
            q = [Q(x, 4) for x in ticks]
            weights = simplex_routing_weights(q)
            # Independent common-threshold enumeration, including all ties.
            counts = {}
            for ticket in range(4):
                mask = sum((1 << j) for j, x in enumerate(ticks) if ticket < x)
                counts[mask] = counts.get(mask, 0)+1
                threshold_tickets += 1
            assert {m: w for m, w in weights.items() if w} == {m: Q(c, 4) for m, c in counts.items()}
            for j in range(d):
                assert sum(w for m, w in weights.items() if m & (1 << j)) == q[j]
                marginal_controls += 1
            assert sum(w > 0 for w in weights.values()) <= d+1
            grid_queries += 1
    for d in range(1, 9):
        for target_mask in range(1 << d):
            for delta in (Q(0), Q(1, 16), Q(1, 8), Q(1, 4)):
                q = [1-delta if target_mask & (1 << j) else delta for j in range(d)]
                weights = simplex_routing_weights(q)
                p = weights.get(target_mask, Q(0))
                assert p >= 1-2*delta
                if target_mask not in (0, (1 << d)-1):
                    assert p == 1-2*delta
                assert 2*(1-p) <= 4*delta
                assert sum((1 << j) for j, x in enumerate(q) if x >= Q(1, 2)) == target_mask
                corner_controls += 1

    rng = random.Random(130066)
    updates = 0
    for d in (2, 4, 6):
        values = {m: (Q(0), Q(0), Q(0)) for m in range(1 << d)}
        for _ in range(64):
            target_mask = rng.randrange(1 << d)
            values[target_mask] = tuple(Q(rng.randrange(-8, 9), 8) for _ in range(3))
            ticks = [rng.randrange(9) for _ in range(d)]
            weights = simplex_routing_weights([Q(x, 8) for x in ticks])
            direct = [sum(values[sum((1 << j) for j, x in enumerate(ticks) if ticket < x)][k]
                          for ticket in range(8))/8 for k in range(3)]
            assert direct == [sum(w*values[m][k] for m, w in weights.items()) for k in range(3)]
            point = simplex_routing_weights([Q(bool(target_mask & (1 << j))) for j in range(d)])
            assert point[target_mask] == 1
            assert tuple(sum(w*values[m][k] for m, w in point.items()) for k in range(3)) == values[target_mask]
            updates += 1

    global_slopes, full_rank_cases = [], 0
    for d in range(1, 9):
        q = [Q(d-j, d+1) for j in range(d)]
        delta = Q(1, 8*(d+1))
        gradient = [Q(-2 if j % 2 == 0 else 2) for j in range(d)]
        perturbed = [x+delta*(1 if g > 0 else -1) for x, g in zip(q, gradient)]
        weights, changed = simplex_routing_weights(q), simplex_routing_weights(perturbed)
        alternating = {mask: Q((-1)**r) for r, mask in enumerate(weights)}
        difference = sum(w*alternating[m] for m, w in changed.items())-sum(w*alternating[m] for m, w in weights.items())
        assert difference == 2*d*delta
        global_slopes.append({'dimension': d, 'input_linf_change': str(delta),
                              'output_change': str(difference), 'sharp_global_linf_constant': 2*d})
        # V_mask[j]=2*bit_j(mask)-1 gives R(q)=2*q-1 and Jacobian 2*I.
        masks = list(weights)
        for position in range(d):
            derivative = [2*(int(bool(masks[position+1] & (1 << j)))-int(bool(masks[position] & (1 << j))))
                          for j in range(d)]
            assert derivative == [2*int(j == position) for j in range(d)]
        full_rank_cases += 1
    rank_controls = 0
    for d in range(1, 6):
        for k in range(1, d+3):
            for _ in range(3):
                values = [[Q(rng.randrange(-8, 9), 8) for _ in range(d)] for _ in range(k)]
                gradients = [[Q(rng.randrange(-4, 5), 8) for _ in range(d)] for _ in range(k-1)]
                gradients.append([-sum(row[j] for row in gradients) for j in range(d)])
                jacobian = [[sum(values[h][i]*gradients[h][j] for h in range(k)) for j in range(d)] for i in range(d)]
                columns, _, _ = gram_leverages(jacobian)
                assert len(columns) <= min(d, k-1)
                rank_controls += 1
    symbolic = []
    for d in (16, 32, 64, 128):
        delta, mask = Q(1, 8), sum(1 << j for j in range(0, d, 2))
        weights = simplex_routing_weights([1-delta if mask & (1 << j) else delta for j in range(d)])
        assert weights[mask] == Q(3, 4)
        symbolic.append({'dimension': d, 'stored_vertices_required': str(1 << d),
                         'simplex_target_mass': '3/4', 'independent_product_target_mass': str((1-delta)**d),
                         'simplex_worst_target_error': '1/2', 'table_was_not_allocated': True})
    return {'dyadic_grid_queries': grid_queries, 'exact_common_threshold_tickets': threshold_tickets,
            'coordinate_marginal_controls': marginal_controls, 'corner_noise_and_hard_rounding_controls': corner_controls,
            'vector_payload_overwrites': updates, 'full_jacobian_examples': full_rank_cases,
            'normalized_common_scalar_mixture_rank_controls': rank_controls,
            'global_sensitivity_witnesses': global_slopes, 'large_dimension_symbolic_controls': symbolic,
            'classification': 'Known simplex/Lovasz routing; local corner retention and derivative/read-span observations are elementary consequences. Exponential table storage and arbitrary raw-key encoder are not free.',
            'novelty_established': False, 'new_publishable_advantage_established': False}


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
    dynamic_gram_sampler = check_dynamic_gram_sampler()
    simplex_routing = check_simplex_routing()
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
              "dynamic_gram_sampler": dynamic_gram_sampler,
              "simplex_routing": simplex_routing,
              "novelty_established": False, "new_publishable_advantage_established": False,
              "source_sha256": {f: sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/"results/research/correlated_routing_checks.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(summary, sort_keys=True))
    print(json.dumps(weighted_prefix_join, sort_keys=True))
    print(json.dumps(dynamic_gram_sampler, sort_keys=True))
    print(json.dumps(simplex_routing, sort_keys=True))


if __name__ == "__main__":
    main()
