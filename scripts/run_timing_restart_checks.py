#!/usr/bin/env python3
"""Exact bounded controls for a candidate timing-adaptive conversion.

This does not implement ANN or establish novelty/performance. It optimizes
over all policies in small timing-only games, including biased surviving
seed posteriors, rather than assuming independent queries after survival.
"""

from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def compositions(total, width):
    if width == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, width - 1):
            yield (first,) + rest


def renewal_controls():
    distributions = truncations = horizons = 0
    for m in range(2, 7):
        for r in range(1, m // 2 + 1):
            a = F(r, m)
            for counts in compositions(m, m):
                if any(sum(counts[:s]) > r * s for s in range(1, m + 1)):
                    continue
                probabilities = [F(c, m) for c in counts]
                distributions += 1
                for h in range(1, m + 1):
                    z = sum(probabilities[:h])
                    mean_truncated = sum(
                        p * min(i + 1, h) for i, p in enumerate(probabilities)
                    )
                    assert z <= a * h
                    assert z <= 2 * a * mean_truncated
                    truncations += 1
                expected = [F(0)]
                for q in range(1, 25):
                    value = sum(
                        probabilities[i - 1] * (1 + expected[q - i])
                        for i in range(1, min(q, m) + 1)
                    )
                    expected.append(value)
                    assert value <= 2 * a * q
                    horizons += 1
    # A near-sharp example; no asymptotic theorem is inferred from a decimal.
    m, limit = 16, 256
    expected = [F(0)]
    for q in range(1, limit + 1):
        expected.append(sum(1 + expected[q - i]
                            for i in range(1, min(m, q) + 1)) / m)
    ratio = expected[-1] / F(limit, m)
    assert F(18, 10) < ratio < F(2)
    return {"cdf_distributions": distributions, "truncated_mean_checks": truncations,
            "iid_renewal_horizons": horizons,
            "near_sharp_uniform_example": {"m": m, "horizon": limit,
                                           "ratio_exceptions_to_Qp": str(ratio)},
            "scope": "Finite IID controls support the general analytical epoch proof; they do not exhaust adaptive epoch laws."}


def adaptive_policy_controls():
    cases = []
    posterior_states = 0
    for m in range(2, 6):
        for k in (1, 2):
            seeds = list(product(range(m), repeat=k))
            full = (1 << len(seeds)) - 1
            for r in range(1, m // 2 + 1):
                bad_masks = []
                for subset in combinations(range(m), r):
                    bad_masks.append(sum(1 << i for i, seed in enumerate(seeds)
                                         if all(s in subset for s in seed)))
                p = F(r, m) ** k
                assert all(F(mask.bit_count(), len(seeds)) == p for mask in bad_masks)

                @lru_cache(None)
                def optimum(q, belief):
                    if q == 0:
                        return F(0)
                    size = belief.bit_count()
                    best = F(0)
                    for bad in bad_masks:
                        failures = (belief & bad).bit_count()
                        good = belief & ~bad
                        value = F(failures, size) * (1 + optimum(q - 1, full))
                        if good:
                            value += F(good.bit_count(), size) * optimum(q - 1, good)
                        best = max(best, value)
                    return best

                maxima = []
                for q in range(1, 25):
                    value = optimum(q, full)
                    assert value <= 2 * p * q
                    maxima.append(value / (p * q))
                posterior_states += optimum.cache_info().currsize
                cases.append({"seed_alphabet": m, "copies": k, "bad_set_size": r,
                              "fresh_timeout_probability": str(p),
                              "maximal_ratio_exceptions_to_Qp": str(max(maxima)),
                              "horizons": 24})
    # Surviving cheap replies can make the next failure certain.
    m = 8
    posterior = set(range(m))
    for query in range(m - 1):
        posterior.remove(query)
    conditional_hazard = F(len(posterior & {m - 1}), len(posterior))
    assert conditional_hazard == 1 > F(1, m)
    return {"optimized_games": cases, "belief_horizon_states": posterior_states,
            "conditional_hazard_counterexample": {"fresh_hazard": "1/8",
                                                   "after_seven_ordinary_replies": str(conditional_hazard)},
            "scope": "All timing-adaptive query policies in each finite game, with fresh independent copies only after exceptions."}


def rebuilding_controls():
    # Nonuniform seeds, observed failure-slot count, and a biased accepted law.
    weights = (F(1, 8), F(3, 8), F(1, 4), F(1, 4))
    accepted = (0, 1)
    acceptance = sum(weights[i] for i in accepted)
    checks = 0
    for j in range(1, 17):
        probability_j = (1 - acceptance) ** (j - 1) * acceptance
        for i in accepted:
            joint = (1 - acceptance) ** (j - 1) * weights[i]
            assert joint / probability_j == weights[i] / acceptance
            checks += 1
    # A build-time observation can distinguish the two installed seeds if
    # successful attempts are not padded, even with no failed attempts.
    leaked_times = {0: 1, 1: 2}
    assert len(set(leaked_times.values())) == len(accepted)
    # Conditioning costs must be bounded relative to the accepted seed law.
    atoms = range(10)
    build_cost = [101 if s == 9 else 1 for s in atoms]
    peak = [101 if s == 8 else 1 for s in atoms]
    R, M = F(sum(build_cost), 10), F(sum(peak), 10)
    survivors = [s for s in atoms if build_cost[s] <= 4 * R and peak[s] <= 4 * M]
    assert F(len(survivors), 10) >= F(1, 2)
    conditioning_checks = 0
    for costs in product((0, 1), repeat=10):
        original = F(sum(costs), 10)
        conditioned = F(sum(costs[s] for s in survivors), len(survivors))
        assert conditioned <= 2 * original
        conditioning_checks += 1
    return {"slot_seed_independence_checks": checks,
            "accepted_seed_conditioning_checks": conditioning_checks,
            "unpadded_build_time_seed_leak": True}


def backlog_controls():
    # A first-moment-only counter can need a logarithmic work budget.
    # Outcomes (j,t) are disjoint; t ranges over the j-th dyadic block.
    records = []
    outcomes = 0
    for levels in range(1, 10):
        length = 2 ** levels
        for divisor, expected_busy in ((3, F(1)), (1, F(1, 3))):
            budget = F(levels, divisor)
            mass = busy = F(0)
            means = [F(0)] * length
            for j in range(levels):
                width = 2 ** j
                for position in range(length - 2 * width + 1, length - width + 1):
                    spike = 3 * budget * width
                    probability = 1 / spike
                    counter = F(0)
                    for t in range(1, length + 1):
                        counter = max(F(0), counter + (spike if t == position else 0) - budget)
                    assert counter > 0
                    means[position - 1] += probability * spike
                    mass += probability
                    busy += probability
                    outcomes += 1
            assert mass <= 1 and busy == expected_busy
            assert all(mean <= 1 for mean in means)
            # Remaining probability is the identically zero-cost outcome.
            records.append({"phase_length": length, "budget": str(budget),
                            "per_operation_mean_bound": "1",
                            "final_busy_probability": str(busy)})
    return {"dyadic_spike_outcomes": outcomes, "witnesses": records,
            "scope": "Necessity of log L order for the generic first-moment counter argument, not a lower bound for all indexes or compilers."}


def canonical_lookup_controls():
    key_sets = [(5, (0, 31)), (6, (0, 63)), (8, (0, 31, 227))]
    choices = overwrites = 0
    for d, keys in key_sets:
        assert all((a ^ b).bit_count() > 4 for a, b in combinations(keys, 2))
        payloads = [(i, -i) for i in range(len(keys))]
        for step in range(8):
            i = step % len(keys)
            payloads[i] = (step - 4, 7 - step)
            overwrites += 1
            for q in range(1 << d):
                inner = [i for i, key in enumerate(keys) if (key ^ q).bit_count() <= 1]
                outer = [i for i, key in enumerate(keys) if (key ^ q).bit_count() <= 2]
                assert len(outer) <= 1
                allowed = inner if inner else [None] + outer
                expected = payloads[inner[0]] if inner else None
                for candidate in allowed:
                    actual = (payloads[candidate] if candidate is not None and
                              (keys[candidate] ^ q).bit_count() <= 1 else None)
                    assert actual == expected
                    choices += 1
    return {"raw_query_candidate_choices": choices, "independent_payload_overwrites": overwrites,
            "scope": "Exact finite reduction controls, with no ANN index, timing benchmark, learned encoder or dynamic performance claim."}


def covering_masks(labels, radius):
    return tuple(sum((((label & u).bit_count() & 1) << j)
                     for j, label in enumerate(labels))
                 for u in range(1, 2 ** (radius + 1)))


def covering_memory_controls():
    families = cover_checks = collision_expectations = 0
    for radius, dimension in ((1, 3), (1, 4), (2, 3), (2, 4)):
        alphabet = 2 ** (radius + 1)
        H = alphabet - 1
        sums = [0] * (1 << dimension)
        count = 0
        for labels in product(range(alphabet), repeat=dimension):
            masks = covering_masks(labels, radius)
            assert len(masks) == H  # Indexed masks retain multiplicity.
            count += 1
            for difference in range(1 << dimension):
                collisions = sum((difference & mask) == 0 for mask in masks)
                sums[difference] += collisions
                if difference.bit_count() <= radius:
                    assert collisions >= 1
                    cover_checks += 1
        for difference, total in enumerate(sums):
            assert F(total, count) == F(H, 2 ** difference.bit_count())
            collision_expectations += 1
        families += count

    # Three copies share the raw keys. Removal must precede shared-key edits.
    # Python buckets certify membership only; the proof costs balanced trees.
    snapshots = reads = mean_reads = empty_reads = payload_writes = key_edits = touches = 0
    witness = None
    for scenario in range(8):
        d, radius, outer, capacity, scale = 8, 1, 2, 4, 16
        keys = [0, 1, 254, 255]
        payloads = [(16, 0, -16), (-16, 16, 0), (0, -16, 16), (8, 8, -8)]
        copies = []
        for copy in range(3):
            labels = ([0] * d if scenario == copy == 0 else
                      [((scenario + 1) * (j + 3) + copy * (j * j + 1)) % 4
                       for j in range(d)])
            masks = covering_masks(labels, radius)
            tables = [{} for _ in masks]
            for i, key in enumerate(keys):
                for mask, table in zip(masks, tables):
                    table.setdefault(key & mask, set()).add(i)
            copies.append((masks, tables))

        def edit_key(i, new):
            nonlocal touches
            old = keys[i]
            if old is not None:
                for masks, tables in copies:
                    for mask, table in zip(masks, tables):
                        signature = old & mask
                        table[signature].remove(i)
                        if not table[signature]:
                            del table[signature]
                        touches += 1
            keys[i] = new
            if new is not None:
                for masks, tables in copies:
                    for mask, table in zip(masks, tables):
                        table.setdefault(new & mask, set()).add(i)
                        touches += 1

        def query(copy, q):
            masks, tables = copy
            candidates, visits = set(), 0
            for mask, table in zip(masks, tables):
                bucket = table.get(q & mask, ())
                visits += len(bucket)
                candidates.update(bucket)
            result = tuple((i, (keys[i] ^ q).bit_count()) for i in sorted(candidates)
                           if (keys[i] ^ q).bit_count() <= radius)
            return result, visits

        def blend(neighbors):
            if not neighbors:
                return None
            weights = [(i, radius + 1 - distance) for i, distance in neighbors]
            denominator = scale * sum(a for _, a in weights)
            return tuple(F(sum(a * payloads[i][j] for i, a in weights), denominator)
                         for j in range(3))

        observed_visits = 0
        for step in range(32):
            # Correctness is tested even with unprotected body-cost feedback.
            # No expected-time guarantee is asserted for this unwrapped oracle.
            i = (step + observed_visits) % capacity
            if step % 4 == 0:
                before = repr(copies)
                payloads[i] = tuple(((step + 7) * (j + 2) % 33) - 16 for j in range(3))
                assert repr(copies) == before
                payload_writes += 1
            elif step % 4 == 2:
                edit_key(i, None)
                key_edits += 1
            else:
                center = 0 if i < 2 else 255
                options = [center] + [center ^ (1 << j) for j in range(d)]
                options = [x for x in options if all(k == i or value != x
                                                   for k, value in enumerate(keys))]
                edit_key(i, options[(step + observed_visits) % len(options)])
                key_edits += 1
            snapshots += 1
            for q in range(1 << d):
                direct = tuple((j, (key ^ q).bit_count()) for j, key in enumerate(keys)
                               if key is not None and (key ^ q).bit_count() <= radius)
                density = sum(key is not None and (key ^ q).bit_count() <= outer for key in keys)
                assert density <= 2
                expected = blend(direct)
                mean_reads += len(direct) > 1
                empty_reads += not direct
                for copy in copies:
                    actual, visits = query(copy, q)
                    assert actual == direct and blend(actual) == expected
                    reads += 1
                    if copy is copies[0] and q == (step * 73) % 256:
                        observed_visits = visits
                if witness is None and len(direct) == 2 and expected not in (
                        tuple(F(z, scale) for z in payloads[j]) for j, _ in direct):
                    witness = {"dimension": d, "query": q,
                               "present_keys": [None if key is None else key for key in keys],
                               "neighbors_and_distances": direct,
                               "payload_numerators": payloads,
                               "payload_scale": scale, "exact_output": [str(x) for x in expected]}
    assert witness is not None
    # A changed projection proves why deleting by the new key is unsafe.
    masks = covering_masks([1, 2, 3], 1)
    assert any((0 & mask) != (1 & mask) for mask in masks)
    return {"exhaustive_label_families": families,
            "all_seed_cover_controls": cover_checks,
            "exact_collision_expectations": collision_expectations,
            "dynamic_snapshots": snapshots, "complete_copy_reads": reads,
            "multiple_neighbor_queries": mean_reads, "empty_queries": empty_reads,
            "payload_only_overwrites": payload_writes, "slot_key_or_occupancy_edits": key_edits,
            "indexed_entry_insertions_and_deletions": touches,
            "full_neighborhood_mixture_witness": witness,
            "shared_key_delete_before_mutate_control": True,
            "scope": "Known covering-family algebra and dynamic exact-operator controls; bucket oracle is not a balanced-tree benchmark or an implementation of the timed wrapper.",
            "novelty_established": False}


def all_prefix_restart_controls():
    mgf_cases = 0
    envelope_records = []
    for m in range(2, 7):
        z = F(m - 1, m)
        beta = z * (1 - z ** m)
        assert beta == sum(z ** a for a in range(1, m + 1)) / m
        assert z ** m >= F(1, 4) and beta <= F(3, 4)
        for counts in compositions(m, m):
            if any(sum(counts[:s]) > s for s in range(1, m + 1)):
                continue
            assert sum(F(c, m) * z ** (a + 1) for a, c in enumerate(counts)) <= beta
            mgf_cases += 1
        envelope_records.append({"inverse_fresh_hazard": m, "z": str(z), "beta": str(beta)})
    confidence_games = []
    for alphabet, copies, bad_size in ((2, 1, 1), (4, 1, 1), (4, 2, 2)):
        seeds = list(product(range(alphabet), repeat=copies))
        full = (1 << len(seeds)) - 1
        bad_masks = [sum(1 << i for i, seed in enumerate(seeds)
                         if all(s in subset for s in seed))
                     for subset in combinations(range(alphabet), bad_size)]
        p = F(bad_size, alphabet) ** copies
        m = 1 / p
        assert m.denominator == 1
        m = int(m)
        z, horizon = 1 - p, 16
        beta = z * (1 - z ** m)
        for delta in (F(1, 4), F(1, 8), F(1, 16)):
            @lru_cache(None)
            def crossing_probability(t, belief, failures):
                if delta * z ** t >= beta ** failures:
                    return F(1)
                if t == horizon:
                    return F(0)
                size, best = belief.bit_count(), F(0)
                for bad in bad_masks:
                    count = (belief & bad).bit_count()
                    good = belief & ~bad
                    value = F(count, size) * crossing_probability(t + 1, full, failures + 1)
                    if good:
                        value += F(good.bit_count(), size) * crossing_probability(t + 1, good, failures)
                    best = max(best, value)
                return best

            probability = crossing_probability(0, full, 0)
            assert probability <= delta
            confidence_games.append({"seed_alphabet": alphabet, "copies": copies,
                                     "bad_set_size": bad_size, "horizon": horizon,
                                     "delta": str(delta),
                                     "optimal_probability_of_any_prefix_crossing": str(probability),
                                     "belief_time_failure_states": crossing_probability.cache_info().currsize})
    return {"extremal_cdf_mgf_controls": mgf_cases, "exact_envelopes": envelope_records,
            "optimal_all_prefix_confidence_games": confidence_games,
            "scope": "Exact bounded checks of the classical exponential-supermartingale argument; finite games do not prove novelty or replace the all-prefix analytical proof."}


def covering_spectral_controls():
    from collections import Counter

    independence_checks = 0
    for radius, dimension in product((1, 2), (1, 2, 3)):
        alphabet = 2 ** (radius + 1)
        families = [covering_masks(labels, radius)
                    for labels in product(range(alphabet), repeat=dimension)]
        for u, v in combinations(range(alphabet - 1), 2):
            counts = Counter((masks[u], masks[v]) for masks in families)
            assert len(counts) == 4 ** dimension
            assert set(counts.values()) == {len(families) // (4 ** dimension)}
            independence_checks += 1
    # Pairwise independent does not mean mutually independent.
    triples = [covering_masks((label,), 1) for label in range(4)]
    assert all(a ^ b == c for a, b, c in triples)
    assert F(triples.count((0, 0, 0)), 4) == F(1, 4) != F(1, 8)

    simplex_keys = tuple(sum(((v & j).bit_count() % 2) << (j - 1)
                             for j in range(1, 8)) for v in range(8))
    cases = [(3, 1, (0, 1, 7), tuple(range(8))),
             (3, 2, (0, 3, 5, 6), tuple(range(8))),
             (4, 2, (0, 3, 12, 15), tuple(range(16))),
             (6, 1, (0, 3, 60, 63), tuple(range(64))),
             (7, 1, simplex_keys, (0, 1, 3, 127))]
    records = []
    moment_checks = label_families = nonzero_tails = 0
    for dimension, radius, keys, queries in cases:
        H = 2 ** (radius + 1) - 1
        # Rounding odd half-distances upward gives a rational upper bound
        # on each row sum of K_ij=2^(-distance(x_i,x_j)/2).
        A = max(sum(F(1, 2 ** ((x ^ y).bit_count() // 2)) for y in keys)
                for x in keys)
        cutoff = (2 * H * A).numerator // (2 * H * A).denominator
        per_query = []
        for q in queries:
            supports = [x ^ q for x in keys]
            one_mask = [sum((s & mask) == 0 for s in supports)
                        for mask in range(1 << dimension)]
            mu = sum(F(1, 2 ** s.bit_count()) for s in supports)
            second = sum(F(1, 2 ** (a | b).bit_count())
                         for a in supports for b in supports)
            assert mu <= A and second <= A * mu
            assert 0 <= second - mu * mu <= A * A / 4
            per_query.append((q, one_mask, mu, second))
        sums = [0] * len(queries)
        squares = [0] * len(queries)
        tails = [0] * len(queries)
        count = 0
        for labels in product(range(H + 1), repeat=dimension):
            masks = covering_masks(labels, radius)
            count += 1
            for j, (_, one_mask, _, _) in enumerate(per_query):
                C = sum(one_mask[mask] for mask in masks)
                sums[j] += C
                squares[j] += C * C
                tails[j] += C > cutoff
        label_families += count
        for j, (q, _, mu, second) in enumerate(per_query):
            mean = F(sums[j], count)
            variance = F(squares[j], count) - mean * mean
            assert mean == H * mu
            assert variance == H * (second - mu * mu)
            assert variance <= H * A * A / 4
            tail = F(tails[j], count)
            assert tail <= F(1, 4 * H)
            moment_checks += 1
            if tail:
                nonzero_tails += 1
                if sum(record['dimension'] == dimension for record in records) < 2:
                    records.append({"dimension": dimension, "radius": radius,
                                    "keys": list(keys), "query": q,
                                    "rational_kernel_norm_upper_bound": str(A),
                                    "mean": str(mean), "variance": str(variance),
                                    "probability_C_exceeds_2HA": str(tail),
                                    "chebyshev_bound": str(F(1, 4 * H))})
    assert nonzero_tails
    # Multiplicity is charged: omitting the spectral factor would be false.
    n, H = 8, 3
    visits = [sum(n for mask in covering_masks((label,), 1) if not mask)
              for label in range(4)]
    mean = F(sum(visits), 4)
    variance = F(sum(c * c for c in visits), 4) - mean * mean
    assert mean == F(H * n, 2)
    assert variance == F(H * n * n, 4) > F(H, 4)
    return {"pairwise_mask_independence_controls": independence_checks,
            "all_label_families": label_families,
            "exact_mean_variance_and_tail_controls": moment_checks,
            "nonzero_tail_control_count": nonzero_tails,
            "nonzero_tail_witnesses": records,
            "mutual_independence_negative_control": {"three_zero_masks": "1/4",
                                                       "incorrect_product": "1/8"},
            "sharp_spectral_variance_witness": {"dimension": 1, "radius": 1,
                                               "eight_identical_keys": 0, "query": 1,
                                               "kernel_norm": 8, "mean": str(mean),
                                               "variance": str(variance),
                                               "ratio_to_H_lambda_squared_over_4": "1"},
            "scope": "Exact finite controls of the analytical spectral second-moment proof; no numerical eigenvalue approximation or novelty inference."}


def dyadic_covering_hashes(x, labels, coefficients, radius, bits, extra_bit=True):
    modulus = 2 ** (bits + int(extra_bit))
    size = 2 ** (radius + 1)
    values = [0] * size
    for j, (label, coefficient) in enumerate(zip(labels, coefficients)):
        if (x >> j) & 1:
            values[label] = (values[label] + coefficient) % modulus
    total = sum(values) % modulus
    width = 1
    while width < size:
        for start in range(0, size, 2 * width):
            for j in range(start, start + width):
                a, b = values[j], values[j + width]
                values[j] = (a + b) % modulus
                values[j + width] = (a - b) % modulus
        width *= 2
    differences = [(total - value) % modulus for value in values[1:]]
    assert all(value % 2 == 0 for value in differences)
    return tuple(value // 2 for value in differences)


def covering_fingerprint_controls():
    hash_values = exact_reads = false_hit_reads = union_controls = 0
    cases = [(1, 1, 1), (2, 1, 1), (2, 1, 3), (2, 1, 4),
             (2, 2, 2), (3, 1, 1), (3, 1, 2), (3, 2, 1)]
    for dimension, radius, bits in cases:
        modulus, H = 2 ** bits, 2 ** (radius + 1) - 1
        keys = ((0, 1, 5, 7) if dimension == 3 else (0, (1 << dimension) - 1))
        payloads = tuple((i - 2, 2 * i + 1) for i in range(len(keys)))
        for labels in product(range(H + 1), repeat=dimension):
            masks = covering_masks(labels, radius)
            noise = [0] * (1 << dimension)
            coefficient_count = 0
            for coefficients in product(range(modulus), repeat=dimension):
                coefficient_count += 1
                fingerprints = [dyadic_covering_hashes(x, labels, coefficients,
                                                       radius, bits)
                                for x in range(1 << dimension)]
                for x, actual in enumerate(fingerprints):
                    direct = tuple(sum(coefficient for j, coefficient in enumerate(coefficients)
                                       if ((x & mask) >> j) & 1) % modulus
                                   for mask in masks)
                    assert actual == direct
                    hash_values += H
                for q in range(1 << dimension):
                    candidates = set()
                    extra_collision = False
                    for u, mask in enumerate(masks):
                        for i, key in enumerate(keys):
                            if fingerprints[q][u] == fingerprints[key][u]:
                                candidates.add(i)
                                extra_collision |= (q & mask) != (key & mask)
                    noise[q] += extra_collision
                    false_hit_reads += extra_collision
                    neighbors = {i for i in candidates if (keys[i] ^ q).bit_count() <= radius}
                    reference = {i for i, key in enumerate(keys) if (key ^ q).bit_count() <= radius}
                    assert neighbors == reference

                    def blend(indices):
                        weights = [(i, radius + 1 - (keys[i] ^ q).bit_count()) for i in indices]
                        return (sum(weight for _, weight in weights),
                                tuple(sum(weight * payloads[i][j] for i, weight in weights)
                                      for j in range(2)))

                    assert blend(neighbors) == blend(reference)
                    exact_reads += 1
            for count in noise:
                assert F(count, coefficient_count) <= F(len(keys) * H, modulus)
                union_controls += 1
    universal_pairs = coefficient_seeds = 0
    for dimension, bits in product(range(1, 5), range(1, 4)):
        modulus = 2 ** bits
        pairs = list(combinations(range(1 << dimension), 2))
        collision_counts = [0] * len(pairs)
        seeds = 0
        for coefficients in product(range(modulus), repeat=dimension):
            values = [sum(coefficient for j, coefficient in enumerate(coefficients)
                          if (x >> j) & 1) % modulus for x in range(1 << dimension)]
            for j, (x, y) in enumerate(pairs):
                collision_counts[j] += values[x] == values[y]
            seeds += 1
        assert set(collision_counts) == {seeds // modulus}
        universal_pairs += len(pairs)
        coefficient_seeds += seeds
    good = dyadic_covering_hashes(1, (1,), (1,), 1, 1)
    wrong = dyadic_covering_hashes(1, (1,), (1,), 1, 1, extra_bit=False)
    assert good == (1, 0, 1) and wrong == (0, 0, 0)
    return {"direct_dot_product_vs_modular_FHT_hash_values": hash_values,
            "complete_exact_weighted_reads": exact_reads,
            "reads_with_extra_fingerprint_collisions": false_hit_reads,
            "shared_coefficient_union_bound_controls": union_controls,
            "universal_binary_projection_pairs": universal_pairs,
            "universal_hash_coefficient_seeds": coefficient_seeds,
            "missing_extra_bit_negative_control": {"correct": list(good),
                                                    "incorrect": list(wrong)},
            "scope": "Finite fingerprints/FHT and complete-function controls, not a production balanced-tree index, benchmark, learned block or claim that compressed tables are independent."}


def constant_family_controls():
    controls = 0
    for H in (3, 7, 15, 63, 255, 1023):
        for cost_ratio in list(range(2, 65)) + [2 ** j for j in range(7, 41)]:
            level = 1
            while (2 * H) ** level < 2 * cost_ratio:
                level += 1
            copies = 4 * level
            reciprocal_failure = (2 * H) ** copies
            assert level <= cost_ratio
            assert reciprocal_failure >= copies * cost_ratio
            controls += 1
    records = []
    for c in (2, 3, 4, 8):
        witnesses = []
        for exponent in range(c, 41):
            n, radius = 2 ** exponent, exponent // c
            H = 2 ** (radius + 1) - 1
            cost_ratio = 4 * n
            level = 1
            while (2 * H) ** level < 2 * cost_ratio:
                level += 1
            copies = 4 * level
            assert (2 * H) ** copies >= copies * cost_ratio
            witnesses.append({"log2_n": exponent, "radius": radius, "copies": copies})
        records.append({"fixed_c": c, "largest_copy_count_in_controls": max(w['copies'] for w in witnesses),
                        "copy_counts_by_log2_n": {str(w['log2_n']): w['copies'] for w in witnesses}})
    spectral_examples = []
    for B in (2, 4, 8, 16, 32):
        n = 2 * B
        t = (n - 1).bit_length()
        dimension = B * t + 1
        keys = tuple((((1 << t) - 1) << (i * t)) | (z << (B * t))
                     for i in range(B) for z in range(2))
        degrees = [sum((x ^ y).bit_count() <= 2 * t for y in keys) - 1 for x in keys]
        row_bounds = [sum(F(1, 2 ** ((x ^ y).bit_count() // 2)) for y in keys) for x in keys]
        assert set(degrees) == {B}
        assert max(row_bounds) < 3
        assert sum(x.bit_count() <= t for x in keys) == B
        radius = t // 2
        assert 1 <= radius < t
        assert min((keys[2 * i + z] ^ keys[2 * j + zz]).bit_count()
                   for i, j in combinations(range(B), 2) for z, zz in product(range(2), repeat=2)) > 2 * radius
        assert all(sum((x ^ q).bit_count() <= radius for x in keys) == 2 for q in keys)
        spectral_examples.append({"capacity": n, "dimension": dimension, "outer_radius": t,
                                  "overlap_degree": B, "zero_query_outer_density": B,
                                  "rational_row_norm_upper_bound": str(max(row_bounds)),
                                  "inner_radius": radius, "inner_neighbor_count_at_each_key": 2})
    return {"integer_copy_choice_controls": controls,
            "polynomial_family_parameter_controls": records,
            "spectral_promise_beyond_bounded_overlap": spectral_examples,
            "scope": "Exact parameter arithmetic at an illustrative reset/base cost ratio 4*n, not measured costs or a proof of a practical speedup."}


def xor_covering_hashes(x, labels, coefficients, radius):
    """Certificate implementation; RAM cost proof scans packed input bits."""
    size = 2 ** (radius + 1)
    grouped = [0] * size
    for j, (label, coefficient) in enumerate(zip(labels, coefficients)):
        if (x >> j) & 1:
            grouped[label] ^= coefficient
    prefix = [0]
    for value in grouped:
        prefix.append(prefix[-1] ^ value)
    basis = []
    for bit in range(radius + 1):
        width = 2 ** bit
        value = 0
        for start in range(width, size, 2 * width):
            value ^= prefix[start + width] ^ prefix[start]
        basis.append(value)
    hashes = [0] * size
    for bit, value in enumerate(basis):
        width = 2 ** bit
        for lower in range(width):
            hashes[width + lower] = hashes[lower] ^ value
    return tuple(hashes[1:])


def xor_covering_controls():
    direct_values = moment_controls = joint_seed_families = tail_controls = 0
    positive_covariance = None
    cases = [(1, 1, 2, (0,), (1,)),
             (2, 1, 1, (0, 1, 3), tuple(range(4))),
             (2, 1, 2, (0, 1, 3), tuple(range(4))),
             (2, 2, 3, (0, 3), tuple(range(4))),
             (3, 1, 2, (0, 1, 3, 7), tuple(range(8))),
             (3, 2, 2, (0, 3, 5, 6), (0, 1, 7))]
    for dimension, radius, bits, keys, queries in cases:
        H, M = 2 ** (radius + 1) - 1, 2 ** bits
        row_bound = max(sum(F(1, 2 ** ((x ^ y).bit_count() // 2)) for y in keys)
                        for x in keys)
        moments = []
        for q in queries:
            supports = [x ^ q for x in keys]
            mu = sum(F(1, 2 ** s.bit_count()) for s in supports)
            second = sum(F(1, 2 ** (s | t).bit_count()) for s in supports for t in supports)
            moments.append((mu, second - mu * mu))
        energy = sum(F(1, 2 ** (x ^ y).bit_count()) for x in keys for y in keys)
        sums, squares, cross, first, second_table, tails = ([0] * len(queries) for _ in range(6))
        cutoff = 2 * H * (row_bound + 1)
        seeds = 0
        for labels in product(range(H + 1), repeat=dimension):
            masks = covering_masks(labels, radius)
            for coefficients in product(range(M), repeat=dimension):
                hashes = [xor_covering_hashes(x, labels, coefficients, radius)
                          for x in range(1 << dimension)]
                for x, actual in enumerate(hashes):
                    direct = []
                    for mask in masks:
                        value = 0
                        for j, coefficient in enumerate(coefficients):
                            if ((x & mask) >> j) & 1:
                                value ^= coefficient
                        direct.append(value)
                    assert actual == tuple(direct)
                    direct_values += H
                for j, q in enumerate(queries):
                    counts = [sum(hashes[x ^ q][u] == 0 for x in keys) for u in range(H)]
                    C = sum(counts)
                    sums[j] += C
                    squares[j] += C * C
                    first[j] += counts[0]
                    second_table[j] += counts[1]
                    cross[j] += counts[0] * counts[1]
                    tails[j] += C > cutoff
                seeds += 1
        joint_seed_families += seeds
        beta, alpha = F(M - 1, M), F(M - 1, M * M)
        for j, (mu, V) in enumerate(moments):
            mean = F(sums[j], seeds)
            variance = F(squares[j], seeds) - mean * mean
            covariance = F(cross[j], seeds) - F(first[j] * second_table[j], seeds * seeds)
            assert mean == H * (beta * mu + F(len(keys), M))
            assert covariance == alpha * V
            predicted = H * beta * (1 + F(H - 3, M)) * V + H * alpha * (energy - mu * mu)
            assert variance == predicted
            moment_controls += 1
            if M >= max(H, len(keys)):
                assert variance <= H * (row_bound * row_bound / 2 + row_bound)
                assert F(tails[j], seeds) <= F(1, 2 * H)
                tail_controls += 1
            if dimension == 1:
                assert covariance == F(3, 64) and variance == F(63, 64)
                assert variance != H * F(15, 64)
                positive_covariance = {"dimension": 1, "radius": 1, "bits": 2,
                                       "keys": [0], "query": 1,
                                       "cross_table_covariance": str(covariance),
                                       "actual_total_variance": str(variance),
                                       "incorrect_independent_table_variance": "45/64"}
    operation_counts = []
    for radius in range(1, 13):
        size = 2 ** (radius + 1)
        ranges = sum(size // (2 ** (bit + 1)) for bit in range(radius + 1))
        outputs = sum(2 ** bit for bit in range(radius + 1))
        assert ranges == outputs == size - 1
        operation_counts.append({"radius": radius, "tables": size - 1,
                                 "prefix_interval_queries": ranges, "output_xors": outputs})
    # A valid low-norm state can still have an exceptional fingerprint seed.
    keys = tuple(sum(((v & j).bit_count() % 2) << (j - 1)
                     for j in range(1, 16)) for v in range(16))
    row_bound = max(sum(F(1, 2 ** ((x ^ y).bit_count() // 2)) for y in keys) for x in keys)
    assert row_bound == F(31, 16) < 2
    zero_fingerprints = [xor_covering_hashes(x, (1,) * 15, (0,) * 15, 1) for x in keys]
    visits = sum(fp[u] == 0 for fp in zero_fingerprints for u in range(3))
    assert visits == 48 > 2 * 3 * (2 + 1)
    return {"direct_xor_fingerprint_values": direct_values,
            "joint_covering_coefficient_seed_families": joint_seed_families,
            "exact_mean_variance_covariance_controls": moment_controls,
            "spectral_tail_controls": tail_controls,
            "positive_shared_fingerprint_covariance": positive_covariance,
            "linear_operation_counts": operation_counts,
            "exceptional_seed_witness": {"dimension": 15, "capacity": 16, "radius": 1,
                                         "kernel_norm_upper_bound": "31/16",
                                         "all_coefficients_zero": True, "candidate_visits": visits,
                                         "A_equals_two_budget": 18},
            "scope": "Independent direct sums and exact joint-seed moments support the XOR/direct-bucket proof. Counts are RAM operation identities, not Python performance measurements or a novelty determination."}


def direct_bucket_memory_controls():
    class Index:
        def __init__(self, n, radius, labels, coefficients):
            self.n, self.radius = n, radius
            self.labels, self.coefficients = labels, coefficients
            self.H, self.M = 2 ** (radius + 1) - 1, 2 ** ((n - 1).bit_length())
            self.heads = [[-1] * self.M for _ in range(self.H)]
            self.prev = [[-1] * n for _ in range(self.H)]
            self.next = [[-1] * n for _ in range(self.H)]
            self.cache = [None] * n
            self.keys, self.payloads, self.marked = [None] * n, [(0, 0)] * n, [False] * n

        def edit(self, slot, key):
            old = self.cache[slot]
            if old is not None:
                for u, value in enumerate(old):
                    previous, following = self.prev[u][slot], self.next[u][slot]
                    if previous < 0:
                        self.heads[u][value] = following
                    else:
                        self.next[u][previous] = following
                    if following >= 0:
                        self.prev[u][following] = previous
                    self.prev[u][slot] = self.next[u][slot] = -1
            self.keys[slot] = key
            self.cache[slot] = None
            if key is not None:
                hashed = xor_covering_hashes(key, self.labels, self.coefficients, self.radius)
                for u, value in enumerate(hashed):
                    following = self.heads[u][value]
                    self.prev[u][slot], self.next[u][slot] = -1, following
                    if following >= 0:
                        self.prev[u][following] = slot
                    self.heads[u][value] = slot
                self.cache[slot] = hashed

        def validate(self):
            active = {i for i, x in enumerate(self.keys) if x is not None}
            for u in range(self.H):
                visited = set()
                for value, head in enumerate(self.heads[u]):
                    previous, slot = -1, head
                    while slot >= 0:
                        assert slot not in visited
                        visited.add(slot)
                        assert self.cache[slot][u] == value and self.prev[u][slot] == previous
                        previous, slot = slot, self.next[u][slot]
                assert visited == active
            assert not any(self.marked)

        def read(self, q, budget):
            hashes = xor_covering_hashes(q, self.labels, self.coefficients, self.radius)
            touched, neighbors, visits = [], [], 0
            finished = True
            for u, value in enumerate(hashes):
                slot = self.heads[u][value]
                while slot >= 0:
                    visits += 1
                    if visits > budget:
                        finished = False
                        break
                    if not self.marked[slot]:
                        self.marked[slot] = True
                        touched.append(slot)
                        distance = (self.keys[slot] ^ q).bit_count()
                        if distance <= self.radius:
                            neighbors.append((slot, distance))
                    slot = self.next[u][slot]
                if not finished:
                    break
            for slot in touched:
                self.marked[slot] = False
            if not finished:
                return False, None, visits
            denominator = sum(self.radius + 1 - distance for _, distance in neighbors)
            numerators = tuple(sum((self.radius + 1 - distance) * self.payloads[slot][j]
                                   for slot, distance in neighbors) for j in range(2))
            return finished, (denominator, numerators), visits

    n, B, t, dimension, radius = 16, 8, 4, 65, 1
    queries = edits = payload_writes = complete = aborted = mixtures = snapshots = 0
    for seed in range(16):
        labels = tuple((j * (seed + 1) + seed) % 4 for j in range(dimension))
        coefficients = tuple((j * seed + seed) % n for j in range(dimension))
        index = Index(n, radius, labels, coefficients)
        for slot in range(n):
            index.edit(slot, (((1 << t) - 1) << ((slot // 2) * t)) | ((slot % 2) << (2 * B * t)))
            index.payloads[slot] = (slot - 8, 7 - slot)
        for step in range(48):
            slot = step % n
            if step % 3 == 0:
                before = (repr(index.heads), repr(index.prev), repr(index.next), repr(index.cache))
                index.payloads[slot] = (step - 24, 23 - step)
                assert before == (repr(index.heads), repr(index.prev), repr(index.next), repr(index.cache))
                payload_writes += 1
            elif step % 3 == 1:
                index.edit(slot, None)
                edits += 1
            else:
                block = slot // 2 + (B if (step // n) % 2 else 0)
                key = (((1 << t) - 1) << (block * t)) | ((slot % 2) << (2 * B * t))
                index.edit(slot, key)
                edits += 1
            index.validate()
            snapshots += 1
            active = [x for x in index.keys if x is not None]
            norm_upper = max(sum(F(1, 2 ** ((x ^ y).bit_count() // 2)) for y in active) for x in active)
            assert norm_upper < 3
            raw_queries = [0, (1 << dimension) - 1] + active[:6] + [x ^ 1 for x in active[:6]]
            for q in raw_queries:
                reference_neighbors = [(i, (x ^ q).bit_count()) for i, x in enumerate(index.keys)
                                       if x is not None and (x ^ q).bit_count() <= radius]
                reference = (sum(radius + 1 - distance for _, distance in reference_neighbors),
                             tuple(sum((radius + 1 - distance) * index.payloads[i][j]
                                       for i, distance in reference_neighbors) for j in range(2)))
                mixtures += len(reference_neighbors) > 1
                for budget in (1, 2 * index.H * 4, n * index.H):
                    finished, answer, _ = index.read(q, budget)
                    assert not any(index.marked)
                    if finished:
                        assert answer == reference
                        complete += 1
                    else:
                        assert answer is None and budget < n * index.H
                        aborted += 1
                    queries += 1
    assert aborted and mixtures
    return {"dynamic_snapshots": snapshots, "direct_bucket_reads": queries,
            "completed_exact_weighted_reads": complete, "aborted_reads_with_clean_scratch": aborted,
            "multiple_neighbor_queries": mixtures, "slot_key_or_occupancy_edits": edits,
            "payload_only_writes": payload_writes,
            "scope": "Bounded intrusive-list oracle checks deletion handles, raw-key edits, independent payloads and abort cleanup. Selected seeds are correctness controls, not a runtime-distribution experiment or production implementation."}


def sparse_shell_controls():
    from math import comb

    overlap_moments = factorial_moments = tails = 0
    for weight in range(1, 65):
        dimension = 16 * weight * weight
        probabilities = [F(comb(weight, i) * comb(dimension - weight, weight - i),
                           comb(dimension, weight)) for i in range(weight + 1)]
        assert sum(probabilities) == 1
        assert sum(p * 2 ** i for i, p in enumerate(probabilities)) <= F(16, 15)
        assert sum(p * 4 ** i for i, p in enumerate(probabilities)) <= F(16, 13)
        overlap_moments += 2
        for level in range(weight + 1):
            moment = sum(p * comb(i, level) for i, p in enumerate(probabilities)
                         if i >= level)
            exact = F(comb(weight, level) ** 2, comb(dimension, level))
            assert moment == exact
            assert exact <= F(comb(weight, level), (16 * weight) ** level)
            factorial_moments += 1
        if weight % 2 == 0:
            probability = sum(probabilities[weight // 2:])
            assert probability <= F(3, 8 * weight) ** (weight // 2)
            tails += 1
    probability_controls = 0
    records = []
    for weight in range(16, 129, 2):
        n, root_n = 2 ** weight, 2 ** (weight // 2)
        bad_pair_bound = n * n * F(3, 8 * weight) ** (weight // 2)
        next_bound = (4 * n) ** 2 * F(3, 8 * (weight + 2)) ** ((weight + 2) // 2)
        assert next_bound / bad_pair_bound < F(6, weight + 2) <= F(1, 3)
        bernstein_exponent = F(3 * n, 48 + 4 * root_n)
        next_exponent = F(3 * 4 * n, 48 + 8 * root_n)
        assert next_exponent >= 2 * bernstein_exponent
        assert bernstein_exponent >= 3 * weight
        assert bad_pair_bound <= F(6561, 2 ** 24) < F(1, 2048)
        probability_controls += 1
        if weight in (16, 24, 32, 64, 128):
            records.append({"log2_n": weight, "dimension": 16 * weight * weight,
                            "kernel_norm_upper_bound": 3,
                            "outer_zero_query_density": str(n),
                            "pair_failure_bound": str(bad_pair_bound),
                            "row_failure_bound": str(F(1, n * n)),
                            "combined_failure_bound": str(bad_pair_bound + F(1, n * n))})
    necessity_controls = sharpness_controls = 0
    configurations = ([(1, d) for d in range(2, 5)] +
                      [(2, d) for d in range(2, 6)] +
                      [(3, d) for d in range(3, 6)] +
                      [(4, d) for d in range(4, 6)])
    for weight, dimension in configurations:
        n = 2 ** weight
        supports = [sum(1 << j for j in coordinates)
                    for coordinates in combinations(range(dimension), weight)]
        for multiplicities in compositions(n, len(supports)):
            keys = [support for support, count in zip(supports, multiplicities)
                    for _ in range(count)]
            frequencies = [sum((x >> j) & 1 for x in keys) for j in range(dimension)]
            overlap_sum = sum((x & y).bit_count() for x in keys for y in keys)
            assert overlap_sum == sum(f * f for f in frequencies)
            assert dimension * overlap_sum >= n * n * weight * weight
            rayleigh = F(sum(2 ** (x & y).bit_count() for x in keys for y in keys), n * n)
            # Exact integer-power comparison of R >= 2^(weight^2/dimension).
            assert rayleigh ** dimension >= 2 ** (weight * weight)
            if dimension == weight:
                assert rayleigh == n and rayleigh ** dimension == 2 ** (weight * weight)
                sharpness_controls += 1
            necessity_controls += 1
    return {"exact_hypergeometric_exponential_moments": overlap_moments,
            "factorial_moment_identities_and_bounds": factorial_moments,
            "overlap_tail_controls": tails,
            "uniform_failure_bound_parameter_controls": probability_controls,
            "probability_witnesses": records,
            "constant_weight_dimension_necessity_controls": necessity_controls,
            "dimension_necessity_equality_controls": sharpness_controls,
            "scope": "Exact distributional and Jensen/Rayleigh controls; no sampled keys, numerical eigenvalues, training, performance benchmark or novelty inference. The written proof supplies the uniform theorem."}


def compact_timing_attack_controls():
    from math import comb

    # The written monotonicity proof extends the rational base comparison
    # to all a>=16; bounded exact checks do not substitute for that proof.
    packing_controls = 0
    for a in range(16, 129):
        universe = comb(8 * a, 4 * a)
        ball = sum(comb(8 * a, j) for j in range(2 * a + 1))
        assert universe * (8 * a + 1) >= 2 ** (8 * a)
        assert ball <= F(65536, 729) ** a
        greedy_bound = F(729, 256) ** a / (8 * a + 1)
        assert greedy_bound > 2 ** a + 1
        assert F(729, 512) * F(8 * a + 1, 8 * a + 9) > 1
        packing_controls += 1

    intersection_controls = factorial_controls = 0
    for copies in range(1, 5):
        for width in range(1, 33):
            dimension = 16 * (2 ** copies - 1) * width * width
            probabilities = [F(comb(width, i) * comb(dimension - width, width - i),
                               comb(dimension, width)) for i in range(width + 1)]
            assert sum(probabilities) == 1
            exponential = sum(p * 2 ** (copies * i) for i, p in enumerate(probabilities))
            assert exponential <= F(16, 15)
            expansion = F(0)
            for level in range(width + 1):
                exact = F(comb(width, level) ** 2, comb(dimension, level))
                moment = sum(p * comb(i, level) for i, p in enumerate(probabilities)
                             if i >= level)
                assert moment == exact
                assert exact <= comb(width, level) * F(width, dimension) ** level
                expansion += exact * (2 ** copies - 1) ** level
                factorial_controls += 1
            assert expansion == exponential
            G = 2 ** (copies * width)
            second_moment = 1 + (1 - F(1, G)) * exponential
            assert F(1, second_moment) >= F(15, 31)
            intersection_controls += 1

    # Enumerate every subset assignment and every mask in small games.
    # These games test the averaging/Cauchy--Schwarz step exactly; their
    # dimensions deliberately differ from the asymptotic construction.
    averaging_games = subset_assignments = mask_games = 0
    for width, dimension in ((1, 4), (2, 5)):
        blocks = [sum(1 << j for j in indices)
                  for indices in combinations(range(dimension), width)]
        G = 2 ** width
        positive = first = second = 0
        best = 0
        for assignment in product(blocks, repeat=G):
            hits = 0
            for mask in range(1 << dimension):
                W = sum(mask & block == 0 for block in assignment)
                first += W
                second += W * W
                positive += W > 0
                hits += W > 0
                mask_games += 1
            best = max(best, hits)
            subset_assignments += 1
        total = len(blocks) ** G * 2 ** dimension
        expected_overlap = sum(2 ** (left & right).bit_count()
                               for left in blocks for right in blocks)
        exponential = F(expected_overlap, len(blocks) ** 2)
        expected_square = 1 + (1 - F(1, G)) * exponential
        assert F(first, total) == 1
        assert F(second, total) == expected_square
        assert F(positive, total) >= 1 / expected_square
        assert F(best, 2 ** dimension) >= F(positive, total)
        averaging_games += 1

    # A separately built, small complete raw-key geometry. It is not a
    # sampled large witness, a novelty claim or an implementation benchmark.
    local_words = (0b00001111, 0b00110011, 0b01010101)
    blocks = (tuple(range(8)), tuple(range(6)) + (8, 9), tuple(range(2, 10)))
    private_dimension, coarse_dimension = 10, 18
    dimension, radius = private_dimension + coarse_dimension, 1
    keys, queries = [], []
    for cluster, block in enumerate(blocks):
        coarse = ((1 << 6) - 1) << (6 * cluster)
        words = [coarse | (sum(((z >> j) & 1) << position for j, position in enumerate(block))
                           << coarse_dimension) for z in local_words]
        queries.append(words[0])
        keys.extend(words[1:])
    assert len(set(keys + queries)) == len(keys) + len(queries)
    assert {x.bit_count() for x in keys + queries} == {10}
    distances = [[(x ^ y).bit_count() for y in keys] for x in keys]
    assert all(distance % 2 == 0 for row in distances for distance in row)
    row_bound = max(sum(F(1, 2 ** (distance // 2)) for distance in row)
                    for row in distances)
    assert row_bound < 3
    assert all((x ^ q).bit_count() > radius for x in keys for q in queries)
    masks_checked = collapses = preserved_collisions = completed_reads = 0
    for mask in range(1 << private_dimension):
        labels = [j % 4 for j in range(coarse_dimension)] + [
            ((mask >> j) & 1) | ((j % 2) << 1) for j in range(private_dimension)]
        family_masks = covering_masks(labels, radius)
        coefficients = tuple((3 * j + 1) % 8 for j in range(dimension))
        raw_times, compressed_times = [], []
        collapsed = []
        for cluster, (q, block) in enumerate(zip(queries, blocks)):
            raw = sum((x ^ q) & projection == 0 for projection in family_masks for x in keys)
            q_hash = xor_covering_hashes(q, labels, coefficients, radius)
            key_hashes = [xor_covering_hashes(x, labels, coefficients, radius) for x in keys]
            compressed = sum(h == q_hash[u] for u in range(3) for values in key_hashes
                             for h in (values[u],))
            assert compressed >= raw
            if all(mask & (1 << j) == 0 for j in block):
                assert all((keys[2 * cluster + i] ^ q) & family_masks[0] == 0 for i in range(2))
                assert all(key_hashes[2 * cluster + i][0] == q_hash[0] for i in range(2))
                assert raw >= 2 and compressed >= 2
                collapsed.append(cluster)
                collapses += 1
                preserved_collisions += 2
            raw_times.append(raw)
            compressed_times.append(compressed)
            completed_reads += 1
        # Observe complete times, repeat their maximum: every later answer
        # is still bottom. Candidate visits provide a deterministic clock.
        if collapsed:
            selected_raw = max(range(len(queries)), key=raw_times.__getitem__)
            selected_compressed = max(range(len(queries)), key=compressed_times.__getitem__)
            assert raw_times[selected_raw] >= 2
            assert compressed_times[selected_compressed] >= 2
        masks_checked += 1
    assert collapses > 0
    # A separate conditional-distribution check: the distinguished mask
    # vanishes on local coordinates but is one at the flipped auxiliary
    # coordinate. Coefficients remain independent after that conditioning.
    conditional_seeds = conditional_means = 0
    for bits in (1, 2, 3):
        M = 2 ** bits
        compressed_sum = raw_sum = seeds = 0
        for local_labels in product((0, 2), repeat=3):
            for auxiliary_label in (1, 3):
                labels = local_labels + (auxiliary_label,)
                projections = covering_masks(labels, 1)
                q, local_keys = 8, (3, 5)
                raw_count = sum((x ^ q) & projection == 0
                                for projection in projections for x in local_keys)
                assert all((x ^ q) & projections[0] != 0 for x in local_keys)
                for coefficients in product(range(M), repeat=4):
                    q_hash = xor_covering_hashes(q, labels, coefficients, 1)
                    point_hashes = [xor_covering_hashes(x, labels, coefficients, 1)
                                    for x in local_keys]
                    compressed_sum += sum(values[u] == q_hash[u]
                                          for u in range(3) for values in point_hashes)
                    raw_sum += raw_count
                    seeds += 1
        raw_mean, compressed_mean = F(raw_sum, seeds), F(compressed_sum, seeds)
        assert raw_mean == F(1, 2)
        assert compressed_mean == F(6, M) + (1 - F(1, M)) * raw_mean
        assert compressed_mean <= 2 * 3
        conditional_seeds += seeds
        conditional_means += 1

    # This exhausts the raw distinguished-table algebra behind the fresh
    # weight-preserving query cube. Classification from actual whole-query
    # RAM time is justified by the separate written Markov/union proof.
    cube_seeds = cube_queries = classification_coordinates = 0
    local_keys, local_query = (0b0101, 0b1001), 0b0011
    auxiliary_width, local_dimension = 6, 4
    auxiliary_base = (1 << 3) - 1
    q_base = local_query | (auxiliary_base << local_dimension)
    cube_keys = [x | (auxiliary_base << local_dimension) for x in local_keys]
    for copies in (1, 2):
        for auxiliary_masks in product(range(1 << auxiliary_width), repeat=copies):
            projections = [mask << local_dimension for mask in auxiliary_masks]
            good = []
            for j in range(auxiliary_width):
                q = q_base ^ (1 << (local_dimension + j))
                all_copies_collide = all(all((x ^ q) & projection == 0 for x in cube_keys)
                                         for projection in projections)
                assert all_copies_collide == all(mask & (1 << j) == 0 for mask in auxiliary_masks)
                if all_copies_collide:
                    good.append(j)
                classification_coordinates += 1
            zero_half = [j for j in good if not auxiliary_base & (1 << j)]
            one_half = [j for j in good if auxiliary_base & (1 << j)]
            pairs = list(zip(zero_half, one_half))
            seen = {q_base} | {q_base ^ (1 << (local_dimension + j))
                               for j in range(auxiliary_width)}
            for subset in range(1, 1 << len(pairs)):
                q = q_base
                for j, (zero, one) in enumerate(pairs):
                    if subset & (1 << j):
                        q ^= (1 << (local_dimension + zero)) | (1 << (local_dimension + one))
                assert q not in seen
                seen.add(q)
                assert q.bit_count() == q_base.bit_count() == cube_keys[0].bit_count()
                assert all((x ^ q).bit_count() > 1 for x in cube_keys)
                assert all(all((x ^ q) & projection == 0 for x in cube_keys)
                           for projection in projections)
                cube_queries += 1
            cube_seeds += 1

    # B_s here is an abstract, supplied reader-cost constant, not a measured
    # or asserted constant of a production implementation.
    fresh_parameter_controls = 0
    fresh_witnesses = []
    for copies in range(1, 9):
        for reader_constant in (1, 16, 256):
            witness = None
            for k in range(16, 513, 2):
                t, L = (8 * copies + 1) * k, 8 * k
                U = D = 16 * (2 ** copies - 1) * L * L
                dimension = 8 * t + D + U
                H, m = 2 ** (k // 2 + 1) - 1, 2 ** k
                mu = F(U, 2 ** (copies + 1))
                epsilon = F(8 * reader_constant * U * dimension * H, m) + F(8, mu)
                if mu / 2 >= 2 * t + 1 and epsilon <= F(1, 8):
                    assert F(15, 31) - epsilon > F(1, 3)
                    n, G = 2 ** t, 2 ** (copies * L)
                    assert n == G * m
                    assert G * (U + 2) <= n * n
                    assert 2 ** (2 * t + 1) - 1 >= n * n
                    witness = {"permanent_families": copies, "abstract_reader_constant": reader_constant,
                               "even_k": k, "log2_n": t, "dimension_with_auxiliary_block": dimension,
                               "classification_error_bound": str(epsilon),
                               "slow_fresh_query_success_probability_lower_bound": str(F(15, 31) - epsilon)}
                    fresh_parameter_controls += 1
                    break
            assert witness is not None
            fresh_witnesses.append(witness)
    return {"uniform_code_packing_parameter_controls": packing_controls,
            "exact_intersection_exponential_moments": intersection_controls,
            "factorial_moment_identities_and_bounds": factorial_controls,
            "exhaustive_averaging_games": averaging_games,
            "exhaustive_subset_assignments": subset_assignments,
            "exhaustive_subset_mask_games": mask_games,
            "constant_weight_empty_answer_witness": {
                "keys": keys, "queries": queries, "radius": radius,
                "dimension": dimension, "kernel_row_sum_bound": str(row_bound),
                "private_masks_checked": masks_checked,
                "completed_reads_with_identical_empty_outputs": completed_reads,
                "cluster_mask_collapse_events": collapses,
                "raw_collisions_preserved_by_xor_compression": preserved_collisions},
            "conditional_mask_coefficient_joint_seeds": conditional_seeds,
            "conditional_compressed_candidate_mean_controls": conditional_means,
            "fresh_query_cube_auxiliary_seed_families": cube_seeds,
            "fresh_query_cube_coordinate_classifications": classification_coordinates,
            "fresh_distinct_weight_preserving_colliding_queries": cube_queries,
            "fresh_query_probability_parameter_controls": fresh_parameter_controls,
            "fresh_query_parameter_witnesses": fresh_witnesses,
            "scope": "Exact rational packing/moment controls, exhaustive finite averaging games, a constant-weight geometry, conditional seed distributions and fresh-query cube algebra support the written proofs. Whole-RAM-time classification follows from the written expectation/Markov argument, not an empirical timing experiment. No full-index lower bound, neural model or novelty determination."}


def full_range_baseline_controls():
    # Correlation defeats multiplying a false-projection probability by an
    # unconditional expected replication count. A deterministic cap works.
    correlated_controls = acceptance_controls = 0
    correlation_witnesses = []
    for seeds in range(2, 33):
        delta = F(1, seeds)
        replication = [seeds * seeds] + [0] * (seeds - 1)
        bad = [True] + [False] * (seeds - 1)
        actual = sum(F(r, seeds) for r, event in zip(replication, bad) if event)
        product_of_means = delta * F(sum(replication), seeds)
        cap_bound = delta * max(replication)
        assert actual == seeds and product_of_means == 1
        assert actual == cap_bound and actual > product_of_means
        correlation_witnesses.append({"equally_likely_seeds": seeds,
                                      "bad_event_probability": str(delta),
                                      "actual_bad_replication_mean": str(actual),
                                      "incorrect_product_of_means": str(product_of_means),
                                      "valid_cap_bound": str(cap_bound)})
        correlated_controls += 1
        # Every accepted subset of this finite seed space, including sets
        # correlated with the bad event, obeys the conditioning bound.
        if seeds <= 8:
            for accepted in range(1, 1 << seeds):
                count = accepted.bit_count()
                alpha = F(count, seeds)
                conditional_bad_work = sum(F(replication[i], count) for i in range(seeds)
                                           if accepted & (1 << i) and bad[i])
                assert conditional_bad_work <= cap_bound / alpha
                acceptance_controls += 1

    # All identities/comparisons use squared rational distances. No floating
    # eigenvalues, Gaussian sampling or numerical root tests are involved.
    directions = ((F(1), F(0)), (F(3, 5), F(4, 5)),
                  (F(0), F(1)), (F(-3, 5), F(4, 5)), (F(-1), F(0)))
    radial_identities = inner_controls = outer_controls = eligible = 0
    for r1 in (F(1), F(2)):
        for delta in (F(1, 16), F(1, 8), F(1, 4)):
            h = delta * r1
            for r2 in (2 * r1 + 2 * h, 3 * r1):
                assert r2 > r1 + 4 * h
                for rp, rq in product((F(j, 4) for j in range(1, 13)), repeat=2):
                    R1 = h * (-((-rp) // h))
                    R2 = h * (-((-rq) // h))
                    assert 0 <= R1 - rp < h and 0 <= R2 - rq < h
                    for up, uq in product(directions, repeat=2):
                        original = sum((rp * x - rq * y) ** 2 for x, y in zip(up, uq))
                        rounded = sum((R1 * x - R2 * y) ** 2 for x, y in zip(up, uq))
                        projected = sum((R1 * x - R1 * y) ** 2 for x, y in zip(up, uq))
                        assert projected == R1 / R2 * (rounded - (R1 - R2) ** 2)
                        radial_identities += 1
                        if abs(R1 - R2) > r1 + 2 * h:
                            assert original > r1 * r1
                            continue
                        eligible += 1
                        inner_threshold = R1 / R2 * ((r1 + 2 * h) ** 2 - (R1 - R2) ** 2)
                        outer_threshold = R1 / R2 * ((r2 - 2 * h) ** 2 - (R1 - R2) ** 2)
                        assert outer_threshold > 0
                        if original <= r1 * r1:
                            assert projected <= inner_threshold
                            inner_controls += 1
                        if projected <= outer_threshold:
                            assert original <= r2 * r2
                            outer_controls += 1

    coarse_ball_controls = coarse_sphere_controls = center_controls = 0
    for r1 in (F(1, 2), F(1), F(2)):
        for R in (F(0), F(1, 2), F(1), F(2)):
            r2 = r1 + 2 * R
            for rq in (F(j, 4) for j in range(17)):
                if rq <= R + r1:
                    assert rq + R <= r2
                    coarse_ball_controls += 1
                for rp in (F(j, 4) for j in range(17)):
                    # At q=o, only annuli meeting the inner distance interval
                    # need enumeration. At R=0 the data location is o.
                    if rq == 0 and rp <= r1:
                        assert rp <= r2
                        center_controls += 1
                    if R == 0 and rq <= r1:
                        assert rq <= r2
                        center_controls += 1
            for up, uq in product(directions, repeat=2):
                distance = sum((R * x - R * y) ** 2 for x, y in zip(up, uq))
                assert distance <= (2 * R) ** 2
                coarse_sphere_controls += 1

    # Rational rotations followed by two orthogonal blocks: test all raw
    # differences of binary vectors in four dimensions, including zero.
    projection_controls = projection_contractions = lifting_controls = 0
    rotation_blocks = [((F(1), F(0)), (F(0), F(1))),
                       ((F(3, 5), F(4, 5)), (F(-4, 5), F(3, 5)))]
    for block1, block2 in product(rotation_blocks, repeat=2):
        for z in product((-1, 0, 1), repeat=4):
            transformed = [sum(row[j] * z[j] for j in range(2)) for row in block1]
            transformed += [sum(row[j] * z[2 + j] for j in range(2)) for row in block2]
            norm = sum(v * v for v in z)
            scaled_squares = (2 * sum(v * v for v in transformed[:2]),
                              2 * sum(v * v for v in transformed[2:]))
            assert sum(scaled_squares) == 2 * norm
            assert min(scaled_squares) <= norm
            projection_contractions += min(scaled_squares) < norm
            projection_controls += 1
    for dimension in range(1, 5):
        points = list(product((-1, 0, 1), repeat=dimension))
        U2 = F(dimension)
        for x, y in product(points, repeat=2):
            dx = sum(F(a - b) ** 2 for a, b in zip(x, y))
            norm_difference = sum(a * a - b * b for a, b in zip(x, y))
            assert norm_difference ** 2 <= 4 * U2 * dx
            for e in (F(1, 2), F(1, 4), F(1, 8)):
                R2 = U2 * (1 + 1 / e)
                # The denominator of the squared added-coordinate difference
                # is at least 4*(R^2-U^2); this is its rational upper bound.
                added_coordinate_bound = F(norm_difference ** 2, 4 * (R2 - U2))
                assert added_coordinate_bound <= e * dx
                lifting_controls += 1

    # An abstract complete-list accounting check with covering projections.
    # It retains all returned occurrences before deduplication; a first-hit
    # count is intentionally not substituted for an enumeration count.
    complete_reads = far_occurrence_controls = nontrivial_mixtures = 0
    for dimension in range(2, 5):
        keys = tuple(range(1 << dimension))
        payloads = [(3 * i - 7, 11 - 2 * i) for i in keys]
        radius, outer_radius = 1, 2
        for labels in product(range(4), repeat=dimension):
            projections = covering_masks(labels, radius)
            for q in range(1 << dimension):
                occurrences = [i for projection in projections for i in keys
                               if (i ^ q) & projection == 0]
                near = [i for i in occurrences if (i ^ q).bit_count() <= outer_radius]
                far = [i for i in occurrences if (i ^ q).bit_count() > outer_radius]
                assert len(occurrences) == len(near) + len(far)
                outer = [i for i in keys if (i ^ q).bit_count() <= outer_radius]
                assert len(near) <= len(projections) * len(outer)
                inner = sorted({i for i in occurrences if (i ^ q).bit_count() <= radius})
                oracle = [i for i in keys if (i ^ q).bit_count() <= radius]
                assert inner == oracle
                denominator = sum(radius + 1 - (i ^ q).bit_count() for i in inner)
                numerator = tuple(sum((radius + 1 - (i ^ q).bit_count()) * payloads[i][j]
                                      for i in inner) for j in range(2))
                assert denominator == dimension + 2
                assert numerator == tuple(2 * payloads[q][j] +
                                          sum(payloads[q ^ (1 << k)][j] for k in range(dimension))
                                          for j in range(2))
                nontrivial_mixtures += len(inner) > 1
                complete_reads += 1
                far_occurrence_controls += len(far)
    return {"correlated_bad_replication_counterexamples": correlated_controls,
            "correlation_witnesses": correlation_witnesses[:4],
            "accepted_seed_conditioning_controls": acceptance_controls,
            "exact_radial_projection_identities": radial_identities,
            "eligible_radial_annulus_controls": eligible,
            "radial_inner_preservation_controls": inner_controls,
            "radial_outer_pullback_controls": outer_controls,
            "coarse_ball_outer_certificates": coarse_ball_controls,
            "coarse_sphere_outer_certificates": coarse_sphere_controls,
            "degenerate_center_controls": center_controls,
            "one_sided_orthogonal_block_controls": projection_controls,
            "strict_contraction_witnesses": projection_contractions,
            "lifting_squared_distance_bound_controls": lifting_controls,
            "complete_exact_weighted_list_reads": complete_reads,
            "verified_far_point_occurrences": far_occurrence_controls,
            "multiple_neighbor_mixtures": nontrivial_mixtures,
            "scope": "Exact algebraic boundary checks support the written full-range arithmetic-model baseline derivation. No optimal LSF implementation, Gaussian experiment, finite-word compilation, adaptive-key theorem, measured performance or novelty claim."}


def sparse_signature_controls():
    from math import comb

    class SignatureIndex:
        """Small exact trie oracle; all payloads remain independent words."""
        def __init__(self, capacity, dimension, weight, level, radius):
            self.n, self.d, self.t, self.level, self.r = capacity, dimension, weight, level, radius
            assert 1 <= level <= weight and 0 <= radius <= weight - level
            self.b = max(1, (dimension - 1).bit_length())
            self.B = level * self.b
            self.limit = 1 + capacity * comb(weight, level) * self.B
            self.child = [[0, 0] for _ in range(self.limit)]
            self.leaf = [None] * self.limit
            self.free = list(range(1, self.limit))
            self.free_count = self.limit - 1
            self.keys = [None] * capacity
            self.values = [[(slot - 9, 3 * slot + 1), (7 - 2 * slot, slot + 13)]
                           for slot in range(capacity)]
            self.marked = [False] * capacity
            self.lookup_steps = 0
            self.false_candidates = self.query_gates = self.queries = 0

        def signatures(self, key, prefix_length=None):
            support = [j for j in range(self.d) if key & (1 << j)]
            if prefix_length is not None:
                assert self.level <= prefix_length <= len(support)
                support = support[:prefix_length]
            for subset in combinations(support, self.level):
                packed = 0
                for coordinate in subset:
                    packed = (packed << self.b) | coordinate
                assert packed.bit_length() <= self.B
                yield packed

        def get(self, packed):
            node = 0
            for bit in range(self.B - 1, -1, -1):
                self.lookup_steps += 1
                node = self.child[node][(packed >> bit) & 1]
                if not node:
                    return None
            return self.leaf[node]

        def put(self, packed, slot):
            node = 0
            for bit in range(self.B - 1, -1, -1):
                direction = (packed >> bit) & 1
                following = self.child[node][direction]
                if not following:
                    assert self.free_count > 0
                    self.free_count -= 1
                    following = self.free[self.free_count]
                    assert self.child[following] == [0, 0] and self.leaf[following] is None
                    self.child[node][direction] = following
                node = following
            assert self.leaf[node] is None
            self.leaf[node] = slot

        def remove(self, packed, slot):
            node, path = 0, []
            for bit in range(self.B - 1, -1, -1):
                direction = (packed >> bit) & 1
                following = self.child[node][direction]
                assert following
                path.append((node, direction, following))
                node = following
            assert self.leaf[node] == slot
            self.leaf[node] = None
            for parent, direction, current in reversed(path):
                if self.leaf[current] is not None or self.child[current] != [0, 0]:
                    break
                self.child[parent][direction] = 0
                self.free[self.free_count] = current
                self.free_count += 1

        def replace(self, slot, key):
            if key is not None:
                if not 0 <= key < 1 << self.d or key.bit_count() != self.t:
                    return False
                for packed in self.signatures(key):
                    owner = self.get(packed)
                    if owner is not None and owner != slot:
                        return False
            old = self.keys[slot]
            if old is not None:
                for packed in self.signatures(old):
                    self.remove(packed, slot)
            if key is not None:
                for packed in self.signatures(key):
                    self.put(packed, slot)
            self.keys[slot] = key
            return True

        def geometry_snapshot(self):
            return repr((self.child, self.leaf, self.free, self.free_count, self.keys))

        def validate(self):
            expected = {}
            for slot, key in enumerate(self.keys):
                if key is None:
                    continue
                assert key.bit_count() == self.t
                for packed in self.signatures(key):
                    assert packed not in expected
                    expected[packed] = slot
            active = [x for x in self.keys if x is not None]
            assert all((x & y).bit_count() < self.level for x, y in combinations(active, 2))
            assert all((x ^ y).bit_count() > 2 * self.r for x, y in combinations(active, 2))
            seen, actual = set(), {}
            def walk(node, prefix, depth):
                assert node not in seen
                seen.add(node)
                if self.leaf[node] is not None:
                    assert depth == self.B
                    actual[prefix] = self.leaf[node]
                for bit, following in enumerate(self.child[node]):
                    if following:
                        assert depth < self.B
                        walk(following, (prefix << 1) | bit, depth + 1)
            walk(0, 0, 0)
            assert actual == expected
            unused = self.free[:self.free_count]
            assert len(unused) == len(set(unused))
            assert set(unused).isdisjoint(seen)
            assert len(unused) + len(seen) == self.limit
            assert all(self.child[i] == [0, 0] and self.leaf[i] is None for i in unused)
            assert not any(self.marked)

        def read(self, full_query):
            query = full_query & ((1 << self.d) - 1)
            suffix = (full_query >> self.d) & 1
            self.queries += 1
            if abs(query.bit_count() - self.t) > self.r:
                self.query_gates += 1
                return (0, (0, 0))
            target, touched, lookups, before = None, [], 0, self.lookup_steps
            prefix_length = self.level + (query.bit_count() - self.t + self.r) // 2
            assert prefix_length <= self.r + self.level
            for packed in self.signatures(query, prefix_length):
                lookups += 1
                slot = self.get(packed)
                if slot is None or self.marked[slot]:
                    continue
                self.marked[slot] = True
                touched.append(slot)
                distance = (self.keys[slot] ^ query).bit_count()
                if distance <= self.r:
                    target = (slot, distance)
                    break
                self.false_candidates += 1
            assert lookups <= comb(self.r + self.level, self.level)
            assert self.lookup_steps - before <= self.B * lookups
            for slot in touched:
                self.marked[slot] = False
            assert not any(self.marked)
            if target is None:
                return (0, (0, 0))
            slot, base_distance = target
            neighbors = [(z, base_distance + (z != suffix)) for z in range(2)
                         if base_distance + (z != suffix) <= self.r]
            denominator = sum(self.r + 1 - distance for z, distance in neighbors)
            numerator = tuple(sum((self.r + 1 - distance) * self.values[slot][z][j]
                                  for z, distance in neighbors) for j in range(2))
            return denominator, numerator

    snapshots = reads = mixtures = accepted = rejected = deletions = value_writes = 0
    false_candidates = weight_gates = maximum_nodes = 0
    for capacity, dimension, weight, level, radius in ((3, 6, 2, 1, 1),
                                                      (4, 7, 3, 2, 1),
                                                      (2, 8, 4, 2, 2)):
        index = SignatureIndex(capacity, dimension, weight, level, radius)
        candidates = [sum(1 << j for j in subset) for subset in combinations(range(dimension), weight)]
        for slot in range(capacity):
            assert any(index.replace(slot, key) for key in candidates)
        for step in range(48):
            slot = step % capacity
            if step % 4 == 0:
                before = index.geometry_snapshot()
                index.values[slot][(step // 4) % 2] = (3 * step - 17, 29 - 2 * step)
                assert index.geometry_snapshot() == before
                value_writes += 1
            elif step % 4 == 1:
                assert index.replace(slot, None)
                deletions += 1
            elif step % 4 == 2:
                for offset in range(len(candidates)):
                    key = candidates[(17 * step + offset) % len(candidates)]
                    legal = all(other == slot or old is None or (key & old).bit_count() < level
                                for other, old in enumerate(index.keys))
                    before, payloads = index.geometry_snapshot(), repr(index.values)
                    result = index.replace(slot, key)
                    assert result == legal
                    if result:
                        accepted += 1
                        break
                    assert index.geometry_snapshot() == before and repr(index.values) == payloads
                    rejected += 1
            else:
                key = next((old for other, old in enumerate(index.keys)
                            if other != slot and old is not None), 0)
                before, payloads = index.geometry_snapshot(), repr(index.values)
                assert not index.replace(slot, key)
                assert index.geometry_snapshot() == before and repr(index.values) == payloads
                rejected += 1
            index.validate()
            snapshots += 1
            maximum_nodes = max(maximum_nodes, index.limit - index.free_count)
            for full_query in range(1 << (dimension + 1)):
                neighbors = [(i, z, ((key | (z << dimension)) ^ full_query).bit_count())
                             for i, key in enumerate(index.keys) if key is not None for z in range(2)
                             if ((key | (z << dimension)) ^ full_query).bit_count() <= radius]
                reference = (sum(radius + 1 - distance for i, z, distance in neighbors),
                             tuple(sum((radius + 1 - distance) * index.values[i][z][j]
                                       for i, z, distance in neighbors) for j in range(2)))
                assert index.read(full_query) == reference
                mixtures += len(neighbors) > 1
                reads += 1
        false_candidates += index.false_candidates
        weight_gates += index.query_gates
    assert mixtures and rejected and accepted and false_candidates and weight_gates

    parameter_controls, induction_controls = 0, 0
    for a in range(8, 257):
        assert a ** 4 <= 16 ** 4 * 2 ** a
        assert F(a + 1, a) ** 4 <= F(9, 8) ** 4 < 2
        induction_controls += 1
    records = []
    for t in (256, 258, 510, 512, 514, 1022, 1024, 2048, 4096, 8192, 16384, 65536):
        a = t.bit_length() - 1
        ell = (4 * t + a - 1) // a
        d = 16 * t * t
        b = (d - 1).bit_length()
        B = ell * b
        assert 1 <= ell <= t // 2
        assert ell * a <= 5 * t and b <= 2 * a + 6 and B <= 14 * t
        assert 2 ** (3 * t) * 3 ** ell <= (16 * ell) ** ell
        F_count = comb(t, ell)
        for c in (2, 3, 4, 8, 16):
            radius = t // c
            G_count = comb(radius + ell, ell)
            assert radius <= t - ell and G_count <= F_count
            # Exact binomial bound and rational entropy envelope, with no
            # logarithmic rounding or claim of a practical speedup.
            assert F_count * ell ** ell <= (3 * t) ** ell
            assert F(3 * t, ell) <= F(3 * a, 4)
            assert F_count <= F(3 * a, 4) ** ell
            parameter_controls += 1
            if c == 2:
                records.append({"log2_n": t, "dimension": d, "signature_level": ell,
                                "signature_bits": B, "radius": radius,
                                "per_key_signature_count_log2_upper": F_count.bit_length(),
                                "query_signature_count_log2_upper": G_count.bit_length(),
                                "overlap_failure_probability_bound": "2^(-t)",
                                "finite_formula_warning": "Subpolynomial asymptotics converge slowly; these bit lengths are exact integer upper bounds, not timings."})
    return {"dynamic_trie_snapshots": snapshots, "complete_raw_query_reads": reads,
            "genuine_two_payload_mixtures": mixtures, "accepted_key_replacements": accepted,
            "atomic_rejected_proposals": rejected, "center_deletions": deletions,
            "payload_only_writes": value_writes, "verified_false_signature_candidates": false_candidates,
            "raw_weight_empty_answer_gates": weight_gates, "maximum_live_trie_nodes_in_controls": maximum_nodes,
            "uniform_overlap_induction_controls": induction_controls,
            "random_shell_parameter_and_integer_bound_controls": parameter_controls,
            "finite_formula_witnesses": records,
            "scope": "Exact preallocated-trie oracle checks dynamic admission, atomic rejection, raw verification, scratch cleanup and independent weighted payloads against exhaustive Hamming reads. Integer controls support the written random-shell/asymptotic proof. This is a known signature-dictionary baseline consequence, not a new architecture, optimal implementation, performance benchmark or novelty claim."}


def product_spectral_controls():
    from math import comb

    class ProductMemory:
        """Exact fixed-codebook oracle with a capacity-bounded array heap."""
        def __init__(self, factors):
            self.p = factors
            self.L = (4 * factors - 1).bit_length()
            self.d, self.n = 2 * factors * self.L, 2 ** factors
            self.owner = list(range(self.n))
            self.address = list(range(self.n))
            self.raw = [self.encode(a) for a in range(self.n)]
            self.values = [(3 * slot - 19, 23 - 2 * slot) for slot in range(self.n)]
            self.heap = [None] * self.n
            self.size = 0

        def encode(self, address):
            key, block = 0, (1 << self.L) - 1
            for j in range(self.p):
                half = 0 if address & (1 << j) else self.L
                key |= block << (2 * j * self.L + half)
            return key

        def decode(self, raw):
            if not 0 <= raw < 1 << self.d:
                return None
            block, address = (1 << self.L) - 1, 0
            for j in range(self.p):
                first = (raw >> (2 * j * self.L)) & block
                second = (raw >> ((2 * j + 1) * self.L)) & block
                if first == block and second == 0:
                    address |= 1 << j
                elif first != 0 or second != block:
                    return None
            return address

        def replace(self, slot, raw):
            address = None if raw is None else self.decode(raw)
            if raw is not None and (address is None or
                                    self.owner[address] not in (None, slot)):
                return False
            old = self.address[slot]
            if old is not None:
                self.owner[old] = None
            if address is not None:
                self.owner[address] = slot
            self.address[slot], self.raw[slot] = address, raw
            return True

        def geometry(self):
            return repr((self.owner, self.address, self.raw))

        def validate(self):
            assert self.size == 0 and all(item is None for item in self.heap)
            expected = [None] * self.n
            for slot, raw in enumerate(self.raw):
                if raw is not None:
                    a = self.decode(raw)
                    assert a is not None and expected[a] is None
                    expected[a] = slot
                    assert self.address[slot] == a
                else:
                    assert self.address[slot] is None
            assert self.owner == expected

        def push(self, item):
            assert self.size < self.n
            index = self.size
            self.size += 1
            while index:
                parent = (index - 1) // 2
                if self.heap[parent] <= item:
                    break
                self.heap[index] = self.heap[parent]
                index = parent
            self.heap[index] = item

        def pop(self):
            item = self.heap[0]
            self.size -= 1
            tail = self.heap[self.size]
            self.heap[self.size] = None
            if self.size:
                index = 0
                while 2 * index + 1 < self.size:
                    child = 2 * index + 1
                    if child + 1 < self.size and self.heap[child + 1] < self.heap[child]:
                        child += 1
                    if tail <= self.heap[child]:
                        break
                    self.heap[index] = self.heap[child]
                    index = child
                self.heap[index] = tail
            return item

        def read(self, query, radius):
            assert self.size == 0
            block, base, favored, penalties = (1 << self.L) - 1, 0, 0, []
            for j in range(self.p):
                u = ((query >> (2 * j * self.L)) & block).bit_count()
                v = ((query >> ((2 * j + 1) * self.L)) & block).bit_count()
                d0, d1 = self.L + u - v, self.L - u + v
                base += min(d0, d1)
                if d1 < d0:
                    favored |= 1 << j
                penalties.append((abs(d0 - d1), j))
            if base > radius:
                return (0, (0, 0)), (), (0, 0)
            penalties.sort()
            self.push((base, favored, -1))
            denominator, numerator, emitted, child_tests, peak = 0, [0, 0], [], 0, 1
            while self.size:
                distance, a, last = self.pop()
                # Oracle-only crosscheck, separate from the heap-cost proof.
                assert distance == (query ^ self.encode(a)).bit_count() <= radius
                emitted.append(a)
                slot = self.owner[a]
                if slot is not None:
                    weight = radius + 1 - distance
                    denominator += weight
                    for j in range(2):
                        numerator[j] += weight * self.values[slot][j]
                for j in range(last + 1, self.p):
                    penalty, coordinate = penalties[j]
                    child_tests += 1
                    if distance + penalty > radius:
                        break
                    self.push((distance + penalty, a ^ (1 << coordinate), j))
                    peak = max(peak, self.size)
            assert len(emitted) == len(set(emitted))
            assert child_tests <= 2 * len(emitted) and peak <= len(emitted)
            return (denominator, tuple(numerator)), tuple(emitted), (child_tests, peak)

    spectral_controls = pair_distances = varying_coordinates = 0
    for p in range(1, 9):
        model = ProductMemory(p)
        theta = F(1, 2 ** model.L)
        spectrum = [theta ** a.bit_count() for a in range(model.n)]
        width = 1
        while width < model.n:
            for start in range(0, model.n, 2 * width):
                for offset in range(width):
                    left, right = spectrum[start + offset], spectrum[start + offset + width]
                    spectrum[start + offset] = left + right
                    spectrum[start + offset + width] = left - right
            width *= 2
        assert theta * p <= F(1, 4)
        for a, value in enumerate(spectrum):
            assert value == (1 + theta) ** (p - a.bit_count()) * (1 - theta) ** a.bit_count()
            assert F(3, 4) <= value <= F(4, 3)
            spectral_controls += 1
        assert max(spectrum) / min(spectrum) <= F(16, 9)
        assert {key.bit_count() for key in model.raw} == {p * model.L}
        for bit in range(model.d):
            assert sum((key >> bit) & 1 for key in model.raw) == model.n // 2
            varying_coordinates += 1
        for a, b in product(range(model.n), repeat=2):
            assert (model.raw[a] ^ model.raw[b]).bit_count() == 2 * model.L * (a ^ b).bit_count()
            pair_distances += 1

    reads = emitted_addresses = mixtures = empty_address_visits = child_tests = 0
    snapshots = admitted = rejected = deletions = value_writes = 0
    format_rejections = duplicate_rejections = moved_keys = inserted_keys = 0
    peak_heap = pure_basins = perturbation_queries = 0

    def verify(model, query, radius):
        nonlocal reads, emitted_addresses, mixtures, empty_address_visits, child_tests, peak_heap
        active = [(slot, (raw ^ query).bit_count()) for slot, raw in enumerate(model.raw)
                  if raw is not None and (raw ^ query).bit_count() <= radius]
        reference = (sum(radius + 1 - distance for slot, distance in active),
                     tuple(sum((radius + 1 - distance) * model.values[slot][j]
                               for slot, distance in active) for j in range(2)))
        answer, emitted, counts = model.read(query, radius)
        assert answer == reference
        ambient = [a for a in range(model.n) if (model.encode(a) ^ query).bit_count() <= radius]
        assert sorted(emitted) == ambient
        bound = min(model.n, 2 ** min(model.p, (2 * radius) // model.L) *
                    sum(comb(model.p, j) for j in range(min(model.p, radius // model.L) + 1)))
        assert len(emitted) <= bound
        reads += 1
        emitted_addresses += len(emitted)
        mixtures += len(active) > 1
        empty_address_visits += sum(model.owner[a] is None for a in emitted)
        child_tests += counts[0]
        peak_heap = max(peak_heap, counts[1])
        return answer

    for p in (1, 2):
        model = ProductMemory(p)
        for step in range(24):
            slot = step % model.n
            if step % 4 == 0:
                slot = (step // 4) % model.n
                before = model.geometry()
                model.values[slot] = (5 * step - 31, 47 - 3 * step)
                assert model.geometry() == before
                value_writes += 1
            elif step % 4 == 1:
                assert model.replace(slot, None)
                deletions += 1
            elif step % 4 == 2:
                slot = (step // 4) % model.n
                old = model.address[slot]
                a = next((a for a, owner in enumerate(model.owner) if owner is None), old)
                if a is not None:
                    assert model.replace(slot, model.encode(a))
                    inserted_keys += old is None
                    moved_keys += old is not None and old != a
                    admitted += 1
            else:
                duplicate = next((a for a, owner in enumerate(model.owner)
                                  if owner is not None and owner != slot), None)
                if (step // 4) % 2 and duplicate is not None:
                    raw = model.encode(duplicate)
                    duplicate_rejections += 1
                else:
                    raw = model.encode(0) ^ 1
                    format_rejections += 1
                before, values = model.geometry(), repr(model.values)
                assert not model.replace(slot, raw)
                assert model.geometry() == before and repr(model.values) == values
                rejected += 1
            model.validate()
            snapshots += 1
            radius = (step // 4) % (model.d + 1)
            for query in range(1 << model.d):
                verify(model, query, radius)

    # Exhaust every block-count profile in a third-factor example. Larger
    # radii deliberately exercise ties/mixtures; the subpolynomial theorem
    # only claims its stated r<=p/c regime.
    model = ProductMemory(3)
    for profile in product(range(model.L + 1), repeat=2 * model.p):
        query = sum(((1 << count) - 1) << (j * model.L) for j, count in enumerate(profile))
        for radius in (0, model.p, 2 * model.L, model.d // 2, model.d):
            verify(model, query, radius)

    # Identity submatrix and every perturbation in the proved pure basin,
    # checked directly rather than inferred from the exponential Gram rank.
    for p in range(1, 5):
        model = ProductMemory(p)
        for radius in range(p + 1):
            h, eta = max(0, radius - model.L + 1), min(radius, model.L - 2)
            fixed_flip = sum(1 << (2 * j * model.L) for j in range(h))
            for a, raw in enumerate(model.raw):
                center = raw ^ fixed_flip
                assert (center ^ raw).bit_count() == h
                for b, other in enumerate(model.raw):
                    if b != a:
                        assert (center ^ other).bit_count() >= h + 2 * (model.L - 1)
                pure_basins += 1
                for weight in range(eta + 1):
                    for positions in combinations(range(model.d), weight):
                        query = center ^ sum(1 << j for j in positions)
                        answer = verify(model, query, radius)
                        assert answer[0] > 0
                        assert answer[1] == tuple(answer[0] * v for v in model.values[a])
                        perturbation_queries += 1

    # Actual mixtures inside the fast r<=p/c regime, beyond the exhaustive
    # small-dimension boundary cases. No large Gram matrix is constructed.
    model = ProductMemory(16)
    model.values[0], model.values[3] = (-71, 29), (113, -47)
    value_writes += 2
    midpoint_one = model.encode(0) ^ ((1 << model.L) - 1)
    first = verify(model, midpoint_one, 8)
    assert first == (6, tuple(3 * (model.values[0][j] + model.values[1][j]) for j in range(2)))
    midpoint_two = midpoint_one ^ (((1 << model.L) - 1) << (2 * model.L))
    second = verify(model, midpoint_two, 12)
    assert second == (4, tuple(sum(model.values[a][j] for a in range(4)) for j in range(2)))
    fast_mixture_controls, fast_four_value_controls = 2, 1
    model.validate()
    snapshots += 1

    parameter_controls = overlap_obstructions = 0
    records = []
    for p in list(range(2, 129)) + [256, 512, 1024, 4096, 16384, 65536]:
        L = (4 * p - 1).bit_length()
        assert 4 * p <= 2 ** L
        for c in (2, 3, 4, 8):
            radius, t = p // c, p * L
            J, weak = radius // L, (2 * radius) // L
            term = total = 1
            for j in range(1, J + 1):
                term = term * (p - j + 1) // j
                total += term
            candidate_bound = min(2 ** p, 2 ** weak * total)
            if J:
                assert total * J ** J <= (3 * p) ** J
            if radius >= L:
                assert not (t - L + 1 <= t - radius)
                overlap_obstructions += 1
            parameter_controls += 1
            if c == 2 and p in (16, 32, 64, 256, 1024, 65536):
                records.append({"log2_capacity": p, "block_half_length": L,
                                "raw_dimension": 2 * p * L, "radius": radius,
                                "kernel_eigenvalue_lower_bound": "3/4",
                                "kernel_eigenvalue_upper_bound": "4/3",
                                "worst_ambient_neighbor_count_log2_upper": candidate_bound.bit_length(),
                                "unique_subset_level_impossible": radius >= L,
                                "pure_recall_basin_radius": min(radius, L - 2)})
    assert mixtures and empty_address_visits and pure_basins and overlap_obstructions
    assert inserted_keys and moved_keys and format_rejections and duplicate_rejections
    return {"exact_tensor_kernel_eigenvalues": spectral_controls,
            "encoded_pair_distance_identities": pair_distances,
            "coordinates_verified_nonconstant": varying_coordinates,
            "complete_raw_or_block_profile_reads": reads,
            "ambient_addresses_enumerated_exactly_once": emitted_addresses,
            "independent_value_mixtures": mixtures, "empty_address_visits_charged": empty_address_visits,
            "fast_regime_mixture_controls": fast_mixture_controls,
            "fast_regime_four_independent_value_controls": fast_four_value_controls,
            "heap_child_tests": child_tests, "maximum_heap_entries_in_controls": peak_heap,
            "dynamic_snapshots": snapshots, "admitted_key_replacements": admitted,
            "atomic_invalid_key_rejections": rejected, "key_deletions": deletions,
            "invalid_format_rejections": format_rejections,
            "duplicate_owner_rejections": duplicate_rejections,
            "existing_key_address_moves": moved_keys, "empty_slot_key_insertions": inserted_keys,
            "payload_only_writes": value_writes, "pure_recall_identity_rows": pure_basins,
            "exhaustive_pure_basin_perturbation_queries": perturbation_queries,
            "integer_volume_parameter_controls": parameter_controls,
            "unique_subset_level_obstruction_controls": overlap_obstructions,
            "finite_parameter_witnesses": records,
            "scope": "Exact Walsh/tensor algebra and a preallocated-array product decoder checked against direct raw Hamming reads, independent current payloads and full pure-basin perturbation controls. Large-radius boundary cases test correctness, not the subpolynomial regime. This is a classical product-indexing baseline consequence, not a new ML architecture, benchmark or novelty determination."}


def compound_rebuild_controls():
    """Exact finite controls for the random-reset supporting lemma."""
    geometric_moments = conditional_copy_products = 0
    for fallback, slot, copies in product((0, 1, 3, 11), (1, 2, 7), (1, 2, 5, 9)):
        scale = fallback + 2 * copies * slot
        gamma = 1 + F(1, 16 * scale)
        v = gamma ** slot
        assert v < 2
        worst = v / (2 - v)
        assert v <= 1 / (1 - F(slot, 16 * scale))
        assert worst <= 1 / (1 - F(slot, 8 * scale))
        assert (1 - F(slot, 8 * scale)) ** copies >= 1 - F(copies * slot, 8 * scale)
        assert gamma ** fallback <= 1 / (1 - F(fallback, 16 * scale))
        assert ((1 - F(fallback, 16 * scale))
                * (1 - F(copies * slot, 8 * scale))) >= F(15, 16)
        assert gamma ** fallback * worst ** copies <= F(16, 15)
        for acceptance in (F(1, 2), F(2, 3), F(3, 4), F(1)):
            # Closed infinite geometric sum; no discarded tail probability.
            moment = acceptance * v / (1 - (1 - acceptance) * v)
            assert moment <= worst
            assert gamma ** fallback * moment ** copies <= F(16, 15)
            assert moment ** copies <= F(16, 15)  # Initial installation.
            geometric_moments += 1
        # Different conditional laws for successive copies are bounded by
        # tower conditioning, without an IID copy-count assumption.
        moment_product = F(1)
        for i in range(copies):
            acceptance = (F(1, 2), F(2, 3), F(3, 4))[i % 3]
            moment_product *= acceptance * v / (1 - (1 - acceptance) * v)
        assert gamma ** fallback * moment_product <= F(16, 15)
        conditional_copy_products += 1

    cdf_laws = defective_laws = transforms = 0
    for denominator in range(2, 7):
        for numerator in range(1, denominator + 1):
            a = F(numerator, denominator)
            z = 1 / (1 + 4 * a)
            h = int(1 / (2 * a))
            assert (1 + 4 * a) ** (h + 1) > 3
            # The last atom is tau=infinity. The finite atoms are 1..m.
            for counts in compositions(denominator, denominator + 1):
                if any(sum(counts[:s]) > numerator * s
                       for s in range(1, denominator + 1)):
                    continue
                probabilities = [F(x, denominator) for x in counts[:-1]]
                transform = sum(p * z ** (i + 1) for i, p in enumerate(probabilities))
                q = sum(probabilities[:h])
                assert q <= F(1, 2)
                assert transform <= q + (1 - q) * z ** (h + 1)
                assert transform < F(2, 3)
                assert F(16, 15) * transform < 1
                cdf_laws += 1
                defective_laws += bool(counts[-1])
                transforms += 1

    # Marginal moment bounds do NOT imply a joint contraction if reset work
    # depends on the trigger. Here K*beta<1 but E[gamma^W*z^tau]>1.
    marginal_cost = F(5, 2)
    epoch_transform = F(3, 8)
    joint_transform = F(9, 8)
    assert marginal_cost * epoch_transform == F(15, 16) < 1
    assert joint_transform > 1
    assert sum(F(1, 2) * 2 ** cost * F(1, 2) ** tau
               for tau, cost in ((1, 2), (2, 0))) == joint_transform

    # With a bound conditional on each trigger, dependence is allowed.
    # Two possible mark laws: P(W=2)=1/4 or 3/4, otherwise W=1.
    # Gamma=2, z=1/2, tau uniform on 1..4: K=7/2, beta=15/64.
    gamma, z, reset_moment, beta = F(2), F(1, 2), F(7, 2), F(15, 64)
    assert reset_moment * beta == F(105, 128) < 1
    dependent_joint_laws = 0
    for choices in product((F(1, 4), F(3, 4)), repeat=4):
        joint = sum(z ** tau * ((1 - p) * gamma + p * gamma ** 2) / 4
                    for tau, p in enumerate(choices, 1))
        assert joint <= reset_moment * beta
        dependent_joint_laws += 1

    # Optimize a finite timing observer's probability of crossing the
    # compound boundary. It may select any singleton bad-seed query,
    # repeat queries, and adapt after observing all reset durations.
    # Four uniform installed seeds give P(tau<=s)<=s/4 for every policy.
    # The rebuild mark depends on tau parity, so marks and epochs are not
    # independent. Fresh accepted seeds remain uniform for either mark.
    initial_moment = F(3)  # Initial W=1 or 2 with probability 1/2 each.
    optimized_games = []
    posterior_states = 0
    positive_crossing_games = 0
    for delta in (F(1, 4), F(1, 2), F(3, 4)):
        @lru_cache(None)
        def crossing(remaining, belief, elapsed, cost, epoch_age):
            if delta * gamma ** cost * z ** elapsed >= initial_moment:
                return F(1)
            if remaining == 0:
                return F(0)
            size = belief.bit_count()
            best = F(0)
            for query in range(4):
                hit = bool(belief & (1 << query))
                probability = F(int(hit), size)
                value = F(0)
                if hit:
                    p_high = F(3, 4) if (epoch_age + 1) % 2 else F(1, 4)
                    value += probability * sum(
                        p * crossing(remaining - 1, 15, elapsed + 1,
                                     cost + work, 0)
                        for work, p in ((1, 1 - p_high), (2, p_high)))
                survivors = belief & ~(1 << query)
                if probability < 1:
                    value += (1 - probability) * crossing(
                        remaining - 1, survivors, elapsed + 1, cost, epoch_age + 1)
                best = max(best, value)
            return best

        for horizon in (2, 4, 6, 8):
            probability = sum(crossing(horizon, 15, 0, work, 0) / 2
                              for work in (1, 2))
            assert probability <= delta
            positive_crossing_games += probability > 0
            optimized_games.append({"horizon": horizon, "delta": str(delta),
                                    "maximum_crossing_probability": str(probability)})
        posterior_states += crossing.cache_info().currsize
    assert positive_crossing_games > 0

    # A build count has an unbounded tail even with a fixed known mean.
    geometric_tail_witnesses = []
    for threshold in (3, 5, 9, 17):
        tail = F(1, 2 ** (threshold - 1))
        assert tail > 0 and threshold > 2
        geometric_tail_witnesses.append({"mean_trials": 2, "threshold": threshold,
                                         "probability_at_least_threshold": str(tail)})
    return {"exact_geometric_moments": geometric_moments,
            "conditional_copy_products": conditional_copy_products,
            "admissible_epoch_cdf_laws": cdf_laws,
            "laws_including_infinite_epoch_mass": defective_laws,
            "exact_epoch_transforms": transforms,
            "dependent_duration_cost_laws": dependent_joint_laws,
            "marginal_moments_counterexample": {
                "marginal_product": str(marginal_cost * epoch_transform),
                "actual_joint_moment": str(joint_transform)},
            "optimized_adaptive_crossing_games": optimized_games,
            "adaptive_posterior_states": posterior_states,
            "positive_crossing_games": positive_crossing_games,
            "unbounded_geometric_build_tail_witnesses": geometric_tail_witnesses,
            "scope": "Exact finite controls for the written uniform compound-supermartingale proof. Finite timing games use bounded two-point rebuild marks, with trigger-dependent laws and fresh uniform seeds; actual unbounded geometric trial counts are checked by their closed infinite generating function. No novelty, architecture advantage or exhaustive infinite-adversary computation is claimed."}


def spectral_omv_controls():
    """Check the geometric reduction, not the conditional OMv conjecture."""
    pair_distances = kernel_rows = rectangle_queries = exact_reads = 0
    rounded_decisions = independent_payload_assignments = 0
    occupancy_decisions = active_principal_rows = 0
    pure_rows = basin_queries = rejected_basin_certificates = 0
    configurations = []

    def pack(bits):
        return sum(bit << i for i, bit in enumerate(bits))

    for m, repeats, slack_bits in product(range(2, 7), (1, 2, 3), (1, 3)):
        capacity = m * m
        label_bits = (capacity - 1).bit_length()
        tag_repeats = label_bits + slack_bits
        tag_weight = label_bits * tag_repeats
        dimension = 4 * repeats * m + 2 * tag_weight
        inner_distance = 2 * repeats * (m - 1) + tag_weight
        radius = inner_distance + repeats - 1
        assert radius == dimension // 2 - repeats - 1
        eta = repeats - 1

        def repeated(mask):
            return [((mask >> i) & 1) for i in range(m) for _ in range(repeats)]

        def query(left, right):
            u, v = repeated(left), repeated(right)
            return pack(u + [1 - bit for bit in u] + v + [1 - bit for bit in v]
                        + [0] * (2 * tag_weight))

        keys = []
        for i, j in product(range(m), repeat=2):
            label = i * m + j
            tag = [value for bit in range(label_bits)
                   for value in ([((label >> bit) & 1)] * tag_repeats
                                 + [1 - ((label >> bit) & 1)] * tag_repeats)]
            assert len(tag) == 2 * tag_weight and sum(tag) == tag_weight
            bits = (repeated(1 << i) + [0] * (repeats * m)
                    + repeated(1 << j) + [0] * (repeats * m) + tag)
            assert len(bits) == dimension
            key = pack(bits)
            assert key.bit_count() == 2 * repeats + tag_weight
            keys.append(key)
        off_diagonal_kernel = []
        for index, key in enumerate(keys):
            row_sum = F(0)
            kernel_row = []
            i, j = divmod(index, m)
            for other, other_key in enumerate(keys):
                k, ell = divmod(other, m)
                distance = (key ^ other_key).bit_count()
                expected = (2 * repeats * ((i != k) + (j != ell))
                            + 2 * tag_repeats * (index ^ other).bit_count())
                assert distance == expected
                pair_distances += 1
                if other != index:
                    assert distance >= 2 * tag_repeats and distance % 2 == 0
                    entry = F(1, 2 ** (distance // 2))
                    row_sum += entry
                    kernel_row.append(entry)
                else:
                    kernel_row.append(F(0))
            assert row_sum <= F(capacity - 1, 2 ** tag_repeats)
            assert row_sum < F(1, 2 ** slack_bits)
            kernel_rows += 1
            off_diagonal_kernel.append(kernel_row)

        if m <= 3:
            payloads = [tuple((mask >> k) & 1 for k in range(capacity))
                        for mask in range(1 << capacity)]
        else:
            payloads = [tuple(0 for _ in keys), tuple(1 for _ in keys),
                        tuple(int(i == j) for i, j in product(range(m), repeat=2)),
                        tuple((i + j) % 2 for i, j in product(range(m), repeat=2)),
                        tuple((i * 3 + j * 5) % 2 for i, j in product(range(m), repeat=2))]
        payloads.append(tuple((k % 3) - 1 for k in range(capacity)))
        independent_payload_assignments += len(payloads)
        for occupancy in payloads:
            if not all(value in (0, 1) for value in occupancy):
                continue
            for index, active in enumerate(occupancy):
                if active:
                    assert sum(entry for entry, occupied
                               in zip(off_diagonal_kernel[index], occupancy)
                               if occupied) < F(1, 2 ** slack_bits)
                    active_principal_rows += 1
        grid = F(1, 1 << (4 * capacity - 1).bit_length())
        assert grid <= F(1, 4 * capacity)
        for left, right in product(range(1 << m), repeat=2):
            q = query(left, right)
            assert q.bit_count() == 2 * repeats * m
            distances = [(q ^ key).bit_count() for key in keys]
            weights = [max(0, radius + 1 - distance) for distance in distances]
            selected = [i * m + j for i, j in product(range(m), repeat=2)
                        if (left & (1 << i)) and (right & (1 << j))]
            assert [k for k, weight in enumerate(weights) if weight] == selected
            for index, distance in enumerate(distances):
                i, j = divmod(index, m)
                expected = (2 * repeats * m + 2 * repeats + tag_weight
                            - 2 * repeats * (((left >> i) & 1) + ((right >> j) & 1)))
                assert distance == expected
                assert weights[index] == (repeats if index in selected else 0)
            denominator = sum(weights)
            assert len(selected) == left.bit_count() * right.bit_count()
            assert denominator == repeats * len(selected)
            rectangle_queries += 1
            for payload in payloads:
                numerator = sum(weight * value for weight, value in zip(weights, payload))
                direct = sum(payload[i * m + j]
                             for i in range(m) if left & (1 << i)
                             for j in range(m) if right & (1 << j))
                actual = F(numerator, denominator) if denominator else None
                expected = F(direct, len(selected)) if selected else None
                assert actual == expected
                exact_reads += 1
                if all(value in (0, 1) for value in payload):
                    # The same Boolean assignment can encode occupancy
                    # instead of payload. All active payloads are then 1.
                    active_weights = [weight for weight, occupied
                                      in zip(weights, payload) if occupied]
                    active_total = sum(active_weights)
                    occupied_read = F(active_total, active_total) if active_total else None
                    assert occupied_read == (F(1) if direct > 0 else None)
                    occupancy_decisions += 1
                    if actual is None:
                        assert direct == 0
                    else:
                        scaled = actual / grid
                        rounded = ((scaled + F(1, 2)).numerator
                                   // (scaled + F(1, 2)).denominator) * grid
                        assert abs(rounded - actual) <= grid / 2
                        assert (rounded > F(1, 2 * capacity)) == (direct > 0)
                    rounded_decisions += 1

        # Certify a pure prototype's full raw Hamming basin using a
        # classical repetition decoder and an explicit final distance test.
        def certify_basin(raw):
            masks = []
            for offset in (0, 2 * repeats * m):
                mask = 0
                for i in range(m):
                    ones = sum((raw >> (offset + i * repeats + k)) & 1
                               for k in range(repeats))
                    ones += sum(1 - ((raw >> (offset + repeats * m + i * repeats + k)) & 1)
                                for k in range(repeats))
                    if ones > repeats:
                        mask |= 1 << i
                masks.append(mask)
            if any(mask.bit_count() != 1 for mask in masks):
                return None
            prototype = query(*masks)
            if (raw ^ prototype).bit_count() > eta:
                return None
            return ((masks[0].bit_length() - 1) * m
                    + masks[1].bit_length() - 1)

        for index in range(capacity):
            i, j = divmod(index, m)
            q = query(1 << i, 1 << j)
            assert certify_basin(q) == index
            assert [k for k, key in enumerate(keys)
                    if (q ^ key).bit_count() <= radius] == [index]
            pure_rows += 1
            if m <= 3 and slack_bits == 1:
                for size in range(eta + 1):
                    for changed in combinations(range(dimension), size):
                        raw = q ^ sum(1 << bit for bit in changed)
                        assert certify_basin(raw) == index
                        neighbors = [k for k, key in enumerate(keys)
                                     if (raw ^ key).bit_count() <= radius]
                        assert neighbors == [index]
                        basin_queries += 1
        for raw in (0, query(0, 0), query(3, 1), query(1, 3)):
            assert certify_basin(raw) is None
            rejected_basin_certificates += 1
        configurations.append({"side": m, "capacity": capacity,
                               "repetition": repeats, "spectral_slack_bits": slack_bits,
                               "raw_dimension": dimension, "fixed_radius": radius,
                               "pure_basin_radius": eta})

    return {"parameter_configurations": configurations,
            "raw_key_pair_distance_identities": pair_distances,
            "exact_kernel_off_diagonal_row_bounds": kernel_rows,
            "complete_rectangle_prototype_queries": rectangle_queries,
            "independent_payload_assignments": independent_payload_assignments,
            "complete_exact_scalar_reads": exact_reads,
            "finite_grid_boolean_decisions": rounded_decisions,
            "occupied_subset_empty_decisions": occupancy_decisions,
            "active_principal_kernel_row_bounds": active_principal_rows,
            "independent_pure_recall_rows": pure_rows,
            "exhaustive_small_raw_basin_queries": basin_queries,
            "invalid_basin_certificates_rejected": rejected_basin_certificates,
            "scope": "Exact bounded certificates for an explicit near-identity-kernel encoding of Boolean rectangle means and a raw-bit repetition-basin decoder. All Boolean payload assignments are exhausted for side 2 and 3; larger cases use stated deterministic payload profiles. Every Boolean rectangle query is exhausted for each listed configuration. All raw bit perturbations up to the basin radius are exhausted only for side <=3 and spectral_slack_bits=1. The uniform spectral/retention/reduction proofs are analytical. The OMv conjecture, arbitrary raw-query performance, small-radius hardness, novelty and an architectural advantage are not verified or claimed."}


def main():
    sources = ["scripts/run_timing_restart_checks.py", "docs/associative_memory_timing.md"]
    report = {"schema_version": 1,
              "title": "Bounded exact controls for the timing-adaptive conversion candidate",
              "status": "passed",
              "renewal": renewal_controls(),
              "adaptive_policies": adaptive_policy_controls(),
              "rebuilding": rebuilding_controls(),
              "backlog": backlog_controls(),
              "canonical_lookup": canonical_lookup_controls(),
              "covering_memory": covering_memory_controls(),
              "all_prefix_restart": all_prefix_restart_controls(),
              "covering_spectral": covering_spectral_controls(),
              "covering_fingerprints": covering_fingerprint_controls(),
              "constant_families": constant_family_controls(),
              "xor_covering": xor_covering_controls(),
              "direct_buckets": direct_bucket_memory_controls(),
              "sparse_shell": sparse_shell_controls(),
              "compact_timing_attack": compact_timing_attack_controls(),
              "full_range_baseline": full_range_baseline_controls(),
              "sparse_signature_baseline": sparse_signature_controls(),
              "product_spectral_baseline": product_spectral_controls(),
              "compound_rebuild_cost": compound_rebuild_controls(),
              "spectral_mixture_hardness": spectral_omv_controls(),
              "novelty_established": False,
              "strongest_baseline_advantage_established": False,
              "publication_ready": False,
              "source_sha256": {f: sha256((ROOT / f).read_bytes()).hexdigest() for f in sources}}
    target = ROOT / "results/research/timing_restart_checks.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"],
                      "adaptive_games": len(report["adaptive_policies"]["optimized_games"]),
                      "cdf_distributions": report["renewal"]["cdf_distributions"],
                      "output": target.relative_to(ROOT).as_posix()}))


if __name__ == "__main__":
    main()
