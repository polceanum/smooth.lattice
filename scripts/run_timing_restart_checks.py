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


def subset_aggregate_controls():
    """Independent Hamming reads versus subset sums and owner corrections."""
    from math import comb, isqrt

    class AggregateIndex:
        # Exact tuple-map control, not a word-RAM benchmark. The written
        # bound uses preallocated tries; the previous signature controls
        # independently exercise that trie/free-list primitive.
        def __init__(self, capacity, dimension, weight, level):
            self.n, self.d, self.w, self.k = capacity, dimension, weight, level
            self.keys = [None] * capacity
            self.values = [(0, 0)] * capacity
            self.handles = [()] * capacity
            self.sums, self.owners = {}, {}
            self.marks = [False] * capacity
            self.entry_limit = capacity * sum(comb(weight, j) for j in range(level + 1))
            self.queries = self.gates = self.candidate_occurrences = 0
            self.candidate_visits = self.corrections = self.low_only = 0
            self.bulk_low_only = self.subset_lookups = self.max_entries = 0
            self.admitted = self.rejected = self.deletions = self.overwrites = 0

        def support(self, key):
            return tuple(i for i in range(self.d) if key & (1 << i))

        def snapshot(self):
            return (tuple(self.keys), tuple(self.values), tuple(self.handles),
                    tuple(sorted((a, count, tuple(total))
                                 for a, (count, total) in self.sums.items())),
                    tuple(sorted(self.owners.items())))

        def remove(self, slot):
            key = self.keys[slot]
            if key is None:
                return
            value = self.values[slot]
            for a in self.handles[slot]:
                count, total = self.sums[a]
                count -= 1
                total = [x - y for x, y in zip(total, value)]
                if count:
                    self.sums[a] = (count, total)
                else:
                    assert total == [0, 0]
                    del self.sums[a]
            for a in combinations(self.support(key), self.k):
                assert self.owners.pop(a) == slot
            self.keys[slot], self.values[slot], self.handles[slot] = None, (0, 0), ()
            self.deletions += 1

        def replace(self, slot, key, value):
            before = self.snapshot()
            valid = 0 <= key < (1 << self.d) and key.bit_count() == self.w
            support = self.support(key) if valid else ()
            if valid:
                valid = all(a not in self.owners or self.owners[a] == slot
                            for a in combinations(support, self.k))
            brute_valid = (0 <= key < (1 << self.d) and key.bit_count() == self.w
                           and all(other is None or i == slot
                                   or (key & other).bit_count() < self.k
                                   for i, other in enumerate(self.keys)))
            assert valid == brute_valid
            if not valid:
                assert self.snapshot() == before
                self.rejected += 1
                return False
            self.remove(slot)
            handles = tuple(a for j in range(self.k) for a in combinations(support, j))
            for a in handles:
                count, total = self.sums.get(a, (0, [0, 0]))
                self.sums[a] = (count + 1, [x + y for x, y in zip(total, value)])
            for a in combinations(support, self.k):
                assert a not in self.owners
                self.owners[a] = slot
            self.keys[slot], self.values[slot], self.handles[slot] = key, tuple(value), handles
            self.admitted += 1
            self.check_state()
            return True

        def overwrite(self, slot, value):
            assert self.keys[slot] is not None
            delta = [x - y for x, y in zip(value, self.values[slot])]
            for a in self.handles[slot]:
                count, total = self.sums[a]
                self.sums[a] = (count, [x + y for x, y in zip(total, delta)])
            self.values[slot] = tuple(value)
            self.overwrites += 1
            self.check_state()

        def check_state(self):
            expected, owners = {}, {}
            for slot, key in enumerate(self.keys):
                if key is None:
                    assert self.handles[slot] == ()
                    continue
                support = self.support(key)
                for j in range(self.k):
                    for a in combinations(support, j):
                        count, total = expected.get(a, (0, [0, 0]))
                        expected[a] = (count + 1, [x + y for x, y in zip(total, self.values[slot])])
                for a in combinations(support, self.k):
                    assert a not in owners
                    owners[a] = slot
            assert self.sums == expected and self.owners == owners
            assert len(self.sums) + len(self.owners) <= self.entry_limit
            assert not any(self.marks)
            self.max_entries = max(self.max_entries, len(self.sums) + len(self.owners))

        def read(self, raw, radius):
            self.queries += 1
            self.last_totals = ((0, 0), 0)
            support = self.support(raw)
            size = len(support)
            if abs(size - self.w) > radius:
                self.gates += 1
                return None
            values = [max(0, radius + 1 - self.w - size + 2 * s)
                      for s in range(self.k)]
            coefficients = []
            while values:
                coefficients.append(values[0])
                values = [y - x for x, y in zip(values, values[1:])]
            polynomial = [sum(coefficient * comb(s, j)
                              for j, coefficient in enumerate(coefficients))
                          for s in range(self.w + 1)]
            assert all(polynomial[s] == max(0, radius + 1 - self.w - size + 2 * s)
                       for s in range(self.k))
            numerator, denominator = [0, 0], 0
            touched = []
            for j, coefficient in enumerate(coefficients):
                for a in combinations(support, j):
                    self.subset_lookups += 1
                    count, total = self.sums.get(a, (0, [0, 0]))
                    denominator += coefficient * count
                    numerator = [x + coefficient * y for x, y in zip(numerator, total)]
            for a in combinations(support, self.k):
                self.subset_lookups += 1
                slot = self.owners.get(a)
                if slot is not None:
                    self.candidate_occurrences += 1
                    if not self.marks[slot]:
                        self.marks[slot] = True
                        touched.append(slot)
            for slot in touched:
                key = self.keys[slot]
                common = (raw & key).bit_count()
                correction = max(0, radius + 1 - self.w - size + 2 * common) - polynomial[common]
                denominator += correction
                numerator = [x + correction * y for x, y in zip(numerator, self.values[slot])]
                self.candidate_visits += 1
                self.corrections += correction != 0
                self.marks[slot] = False
            assert all(not self.marks[slot] for slot in touched) and denominator >= 0
            if not touched:
                self.low_only += 1
            self.last_totals = (tuple(numerator), denominator)
            if denominator == 0:
                assert numerator == [0, 0]
                return None
            return tuple(F(x, denominator) for x in numerator)

        def oracle(self, raw, radius):
            weights = [(slot, max(0, radius + 1 - (raw ^ key).bit_count()))
                       for slot, key in enumerate(self.keys) if key is not None]
            denominator = sum(weight for _, weight in weights)
            numerator = tuple(sum(weight * self.values[slot][coordinate]
                                  for slot, weight in weights) for coordinate in range(2))
            self.oracle_totals = (numerator, denominator)
            if not denominator:
                return None
            return tuple(F(value, denominator) for value in numerator)

    cases = []
    # Fano triples are not a one-hot polynomial codebook.
    families = [(6, 2, 1, [3, 12, 48]),
                (9, 4, 2, [sum(1 << i for i in a) for a in
                           ((0, 1, 2, 3), (0, 4, 5, 6), (1, 4, 7, 8))]),
                (7, 3, 2, [sum(1 << i for i in a) for a in
                           ((0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5),
                            (1, 4, 6), (2, 3, 6), (2, 4, 5))])]
    rs = {}
    for prime, level in ((3, 2), (5, 2), (5, 3), (7, 2), (13, 2)):
        keys = []
        for coefficients in product(range(prime), repeat=level):
            symbols = [sum(c * pow(a, j, prime) for j, c in enumerate(coefficients)) % prime
                       for a in range(prime)]
            keys.append(sum(1 << (a * prime + y) for a, y in enumerate(symbols)))
        families.append((prime * prime, prime, level, keys))
        rs[prime, level] = keys

    complete_reads = dynamic_states = zero_sum_live_entries = 0
    arbitrary_admitted_keys = parameter_bounds = spectral_rows = 0
    pure_noise_queries = admitted_noise_queries = mixture_witnesses = 0
    for dimension, weight, level, keys in families:
        index = AggregateIndex(len(keys), dimension, weight, level)
        for slot, key in enumerate(keys):
            assert index.replace(slot, key, (2 * slot - 11, 17 - 3 * slot))
        if dimension <= 9:
            queries = list(range(1 << dimension))
        else:
            query_set = {0, (1 << dimension) - 1, *keys}
            for size in range(dimension + 1):
                query_set.add((1 << size) - 1)
                query_set.add(((1 << size) - 1) << (dimension - size))
            for key in keys:
                for bit in range(min(dimension, 31)):
                    query_set.add(key ^ (1 << bit))
            if weight == 5:
                query_set.update(sum(1 << (a * weight + symbol)
                                     for a, symbol in enumerate(word))
                                 for word in product(range(weight), repeat=weight))
            queries = sorted(query_set)

        def check(radii):
            nonlocal complete_reads, dynamic_states, zero_sum_live_entries
            index.check_state()
            dynamic_states += 1
            zero_sum_live_entries += sum(count > 0 and total == [0, 0]
                                         for count, total in index.sums.values())
            for radius in radii:
                for raw in queries:
                    before = index.candidate_visits
                    assert index.read(raw, radius) == index.oracle(raw, radius)
                    assert index.last_totals == index.oracle_totals
                    # Oracle-only diagnostic, outside the candidate reader.
                    if index.candidate_visits == before:
                        index.bulk_low_only += sum(
                            key is not None and (raw ^ key).bit_count() <= radius
                            for key in index.keys) >= 3
                    complete_reads += 1

        check(range(2 * weight + 1) if dimension <= 9 else (0, weight, 2 * weight))
        for slot in range(0, len(keys), 3):
            index.overwrite(slot, (0, 0))
        check((weight,))
        # Nonzero independent values can cancel, leaving live zero sums.
        for slot in range(len(keys)):
            value = 0 if slot == len(keys) - 1 and len(keys) % 2 else (-1) ** slot
            index.overwrite(slot, (value, -value))
        assert index.sums[()][1] == [0, 0]
        check((weight,))
        for slot in range(1, len(keys), 2):
            index.remove(slot)
        check((weight,))
        # Invalid weight, duplicate owner and a distinct high-overlap key
        # must all reject without altering sums or current values.
        assert not index.replace(0, 0, (99, -99))
        if len(keys) > 2:
            assert not index.replace(0, keys[2], (99, -99))
        base_support = index.support(keys[2])
        changed = tuple(base_support[:level]) + tuple(
            i for i in range(dimension) if i not in base_support)[:weight - level]
        if len(changed) == weight:
            before = index.snapshot()
            assert not index.replace(0, sum(1 << i for i in changed), (99, -99))
            assert index.snapshot() == before
        # The generic admission rule may accept keys outside the original
        # codebook. Exhaust every weight-w proposal for one small slot.
        if dimension <= 9:
            original = set(keys)
            for support in combinations(range(dimension), weight):
                key = sum(1 << i for i in support)
                accepted = index.replace(1, key, (23, -31))
                arbitrary_admitted_keys += accepted and key not in original
            check((weight,))
            admitted_eta = min(1, weight - 2 * level + 1)
            for slot, key in enumerate(index.keys):
                if key is None or admitted_eta < 0:
                    continue
                for size in range(admitted_eta + 1):
                    for bits in combinations(range(dimension), size):
                        raw = key ^ sum(1 << bit for bit in bits)
                        expected = tuple(F(x) for x in index.values[slot])
                        assert index.read(raw, weight) == expected == index.oracle(raw, weight)
                        admitted_noise_queries += 1
        # Repeated arena-entry deletion/reuse, and cancellation with
        # nonzero support, have no accumulated history parameter.
        for iteration in range(12):
            index.remove(0)
            assert index.replace(0, keys[0], (iteration - 6, 6 - iteration))
        for slot in range(len(keys)):
            index.remove(slot)
        assert index.sums == {} and index.owners == {}
        check((weight,))
        for slot, key in enumerate(keys):
            assert index.replace(slot, key, ((-1) ** slot, -((-1) ** slot)))
        check((weight,))
        if dimension <= 25:
            minimum = min((a ^ b).bit_count() for i, a in enumerate(keys) for b in keys[i + 1:])
            radius = weight
            eta = min(radius, minimum - radius - 1, 2)
            if eta >= 0:
                for slot, key in enumerate(keys):
                    for size in range(eta + 1):
                        for bits in combinations(range(dimension), size):
                            raw = key ^ sum(1 << bit for bit in bits)
                            assert index.read(raw, radius) == tuple(F(x) for x in index.values[slot])
                            assert index.oracle(raw, radius) == tuple(F(x) for x in index.values[slot])
                            pure_noise_queries += 1
        cases.append({"capacity": len(keys), "raw_dimension": dimension,
                      "weight": weight, "intersection_level": level,
                      "queries": index.queries, "weight_gate_returns": index.gates,
                      "subset_lookups": index.subset_lookups,
                      "candidate_occurrences": index.candidate_occurrences,
                      "deduplicated_candidate_visits": index.candidate_visits,
                      "nonzero_owner_corrections": index.corrections,
                      "reads_without_candidate_visits": index.low_only,
                      "bulk_reads_without_candidate_visits": index.bulk_low_only,
                      "maximum_live_subset_entries": index.max_entries,
                      "subset_entry_capacity": index.entry_limit,
                      "admitted_key_replacements": index.admitted,
                      "atomic_invalid_rejections": index.rejected,
                      "key_deletions": index.deletions, "payload_overwrites": index.overwrites})

    # Exact full-code spectral row sums, independently computed from
    # raw Hamming distances rather than a large Gram eigensolver.
    for (prime, level), keys in rs.items():
        predicted = 1 + F(sum(comb(prime, j) * (prime ** (level - j) - 1)
                              for j in range(level)), 2 ** prime)
        assert predicted < 1 + F(3 * len(keys), 2 ** prime)
        for key in keys:
            row = sum(F(1, 2 ** ((key ^ other).bit_count() // 2)) for other in keys)
            assert row == predicted
            spectral_rows += 1

    # Exact 2- and 3-value blends among the linear polynomial keys.
    for prime, count in ((5, 2), (7, 2), (13, 2), (13, 3)):
        assert (prime - 1) % count == 0
        agreement = 1 + (prime - 1) // count
        assert agreement > count  # level=2, every other polynomial hits <=count rows.
        keys = rs[prime, 2]
        index = AggregateIndex(len(keys), prime * prime, prime, 2)
        coefficients = list(product(range(prime), repeat=2))
        chosen = [coefficients.index((0, slope)) for slope in range(count)]
        for slot, key in enumerate(keys):
            assert index.replace(slot, key, (slot * slot + 1, -3 * slot - 2))
        symbols = [0] + [((a - 1) // ((prime - 1) // count)) * a % prime
                         for a in range(1, prime)]
        raw = sum(1 << (a * prime + symbol) for a, symbol in enumerate(symbols))
        radius = 2 * (prime - agreement)
        neighbors = [slot for slot, key in enumerate(keys) if (raw ^ key).bit_count() <= radius]
        assert neighbors == chosen
        expected = tuple(F(sum(index.values[slot][j] for slot in chosen), count) for j in range(2))
        assert index.read(raw, radius) == expected == index.oracle(raw, radius)
        if count == 2:
            assert index.read(raw, prime) == expected == index.oracle(raw, prime)
        index.overwrite(chosen[0], (-37, 41))
        expected = tuple(F(sum(index.values[slot][j] for slot in chosen), count) for j in range(2))
        assert index.read(raw, radius) == expected == index.oracle(raw, radius)
        index.remove(chosen[-1])
        expected = tuple(F(sum(index.values[slot][j] for slot in chosen[:-1]), count - 1)
                         for j in range(2))
        assert index.read(raw, radius) == expected == index.oracle(raw, radius)
        mixture_witnesses += 1

    parameter_examples = []
    for prime in (127, 257, 521, 1031, 2053, 4093):
        assert all(prime % divisor for divisor in range(2, isqrt(prime) + 1))
        level = prime // (4 * (prime - 1).bit_length())
        assert 2 <= level <= prime // 2
        capacity = prime ** level
        total = sum(comb(prime, j) for j in range(level + 1))
        query = sum(comb(3 * prime, j) for j in range(level + 1))
        assert total * level ** level <= (3 * prime) ** level
        assert query * level ** level <= (9 * prime) ** level
        assert capacity ** 4 <= 2 ** prime
        parameter_examples.append({"prime": prime, "level": level,
                                   "capacity_bits": capacity.bit_length(),
                                   "key_subset_count_bits": total.bit_length(),
                                   "query_subset_count_bits": query.bit_length()})
        parameter_bounds += 1
    assert arbitrary_admitted_keys > 0
    assert sum(c['bulk_reads_without_candidate_visits'] for c in cases) > 0
    assert sum(c['nonzero_owner_corrections'] for c in cases) > 0
    return {"cases": cases, "complete_exact_reads": complete_reads,
            "dynamic_snapshots": dynamic_states,
            "zero_vector_sums_retained_with_positive_support": zero_sum_live_entries,
            "admitted_keys_outside_initial_codebook": arbitrary_admitted_keys,
            "exhaustive_small_pure_noise_reads": pure_noise_queries,
            "post_admission_pure_noise_reads": admitted_noise_queries,
            "exact_full_code_spectral_row_sums": spectral_rows,
            "independent_two_or_three_value_mixture_witnesses": mixture_witnesses,
            "integer_asymptotic_parameter_controls": parameter_bounds,
            "parameter_bit_length_examples": parameter_examples,
            "scope": "Exact tuple-map algorithm versus independent complete raw Hamming reads. Small dimensions exhaust all raw queries; larger cases exhaust stated one-hot queries for weight 5 and use systematic raw support/bit-flip controls. Dynamic admission is checked against all pair intersections, including proposals outside the initial codebook. Signed cancellations, atomic rejection, deletion/reuse, pure-noise reads, nontrivial independently overwritten/deleted mixtures, spectral row sums and exact parameter bounds are checked. The code is a bounded correctness control, not a finite-word arena implementation, runtime benchmark, learned encoder or novelty/strongest-baseline determination."}


def affine_coded_controls():
    """Classical soft interpolation plus partial evaluation sums, not novelty."""
    from math import comb

    stats = dict(reads=0, low_reads=0, high_reads=0, empty_gates=0,
                 fiber_addresses_examined=0, interpolation_updates=0,
                 verified_polynomial_roots=0, nonzero_corrections=0,
                 overwrites=0, deletions=0, admissions=0, rejections=0,
                 live_zero_sums=0, high_two_neighbor_queries=0)

    class Field:
        def __init__(self, prime):
            self.p = prime
            self.inverse = [0] + [pow(a, -1, prime) for a in range(1, prime)]

        def trim(self, a):
            a = [x % self.p for x in a]
            while a and a[-1] == 0:
                a.pop()
            return tuple(a)

        def add(self, a, b):
            return self.trim([(a[i] if i < len(a) else 0)
                              + (b[i] if i < len(b) else 0)
                              for i in range(max(len(a), len(b)))])

        def scale(self, a, c):
            return self.trim([x * c for x in a])

        def mul(self, a, b):
            if not a or not b:
                return ()
            out = [0] * (len(a) + len(b) - 1)
            for i, x in enumerate(a):
                for j, y in enumerate(b):
                    out[i + j] += x * y
            return self.trim(out)

        def divide(self, a, b):
            assert b
            remainder = list(a)
            out = [0] * max(0, len(a) - len(b) + 1)
            while len(remainder) >= len(b):
                shift = len(remainder) - len(b)
                c = remainder[-1] * self.inverse[b[-1]] % self.p
                out[shift] = c
                for j, y in enumerate(b):
                    remainder[shift + j] = (remainder[shift + j] - c * y) % self.p
                while remainder and remainder[-1] == 0:
                    remainder.pop()
            return self.trim(out), tuple(remainder)

        def evaluate(self, a, x):
            total = 0
            for c in reversed(a):
                total = (total * x + c) % self.p
            return total

        def substitute(self, q, f):
            total = ()
            for a in reversed(q):
                total = self.add(self.mul(total, f), a)
            return total

        def square_root(self, a):
            if not a:
                return ()
            if (len(a) - 1) % 2:
                return None
            alpha = next((x for x in range(1, self.p)
                          if x * x % self.p == a[-1]), None)
            if alpha is None:
                return None
            m = (len(a) - 1) // 2
            out = [0] * (m + 1)
            out[m] = alpha
            for j in range(1, m + 1):
                exponent = 2 * m - j
                known = sum(out[u] * out[exponent - u]
                            for u in range(m - j + 1, m + 1)
                            if m - j < exponent - u <= m)
                out[m - j] = (a[exponent] - known) * self.inverse[2 * alpha % self.p] % self.p
            out = self.trim(out)
            return out if self.mul(out, out) == a else None

        def roots(self, q, level):
            c, b, a = q
            answers = set()
            if a:
                disc = self.add(self.mul(b, b), self.scale(self.mul(a, c), -4))
                root = self.square_root(disc)
                numerators = [] if root is None else [self.add(self.scale(b, -1), root),
                                                       self.add(self.scale(b, -1), self.scale(root, -1))]
                divisor = self.scale(a, 2)
            elif b:
                numerators, divisor = [self.scale(c, -1)], b
            else:
                assert c  # Interpolation must never return the zero polynomial.
                numerators, divisor = [], (1,)
            for numerator in numerators:
                f, remainder = self.divide(numerator, divisor)
                if not remainder and len(f) <= level:
                    assert not self.substitute(q, f)
                    answers.add(f + (0,) * (level - len(f)))
            assert len(answers) <= 2
            stats['verified_polynomial_roots'] += len(answers)
            return answers

    def interpolate(points, level, field):
        # Multiplicity one, Y degree <=2. Distinct leading Y positions
        # make the weighted-degree argument and its quadratic cost explicit.
        basis = [[(1,), (), ()], [(), (1,), ()], [(), (), (1,)]]

        def leading(q):
            return max((i + (level - 1) * j, j, i)
                       for j, a in enumerate(q) for i, c in enumerate(a) if c)

        for seen, (a, b) in enumerate(points, 1):
            delta = [sum(field.evaluate(poly, a) * pow(b, j, field.p)
                         for j, poly in enumerate(q)) % field.p for q in basis]
            nonzero = [j for j in range(3) if delta[j]]
            if nonzero:
                pivot = min(nonzero, key=lambda j: leading(basis[j]))
                old = basis[pivot]
                old_leads = [leading(q) for q in basis]
                for j in nonzero:
                    if j != pivot:
                        ratio = -delta[j] * field.inverse[delta[pivot]]
                        basis[j] = [field.add(x, field.scale(y, ratio))
                                    for x, y in zip(basis[j], old)]
                        assert leading(basis[j]) == old_leads[j]
                basis[pivot] = [field.mul(poly, (-a, 1)) for poly in old]
                assert leading(basis[pivot]) == (old_leads[pivot][0] + 1,
                                                 pivot, old_leads[pivot][2] + 1)
                stats['interpolation_updates'] += 1
            leads = [leading(q) for q in basis]
            assert [item[1] for item in leads] == [0, 1, 2]
            assert sum(item[2] for item in leads) <= seen
        q = min(basis, key=leading)
        assert all(sum(field.evaluate(poly, a) * pow(b, j, field.p)
                       for j, poly in enumerate(q)) % field.p == 0 for a, b in points)
        assert len(points) > 6 * (level - 1)
        assert leading(q)[0] < (len(points) + 1) // 2
        return q

    def affine_fiber(points, level, field):
        # RREF on the actual evaluation equations, including repeated rows.
        matrix = [[pow(a, j, field.p) for j in range(level)] + [b] for a, b in points]
        pivot_columns, row = [], 0
        for col in range(level):
            pivot = next((j for j in range(row, len(matrix)) if matrix[j][col]), None)
            if pivot is None:
                continue
            matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
            matrix[row] = [x * field.inverse[matrix[row][col]] % field.p for x in matrix[row]]
            for j in range(len(matrix)):
                if j != row:
                    c = matrix[j][col]
                    matrix[j] = [(x - c * y) % field.p for x, y in zip(matrix[j], matrix[row])]
            pivot_columns.append(col)
            row += 1
        if any(not any(line[:level]) and line[-1] for line in matrix):
            return
        free = [j for j in range(level) if j not in pivot_columns]
        for assignment in product(range(field.p), repeat=len(free)):
            f = [0] * level
            for j, x in zip(free, assignment):
                f[j] = x
            for j, line in zip(pivot_columns, matrix):
                f[j] = (line[-1] - sum(line[u] * f[u] for u in free)) % field.p
            yield tuple(f)

    class CodedIndex:
        # Sparse dictionaries model exact occupied cells only. No timings or
        # storage benchmark: the theorem allocates all stated direct arrays.
        def __init__(self, prime, level, cached):
            assert 1 <= cached < level <= prime
            self.p, self.k, self.h = prime, level, cached
            self.n, self.d = prime ** level, prime * prime
            self.field = Field(prime)
            self.slots, self.owners, self.cells = {}, {}, {}
            self.last_totals = self.oracle_totals = None
            self.last_neighbors = set()

        def symbols(self, f):
            return tuple(self.field.evaluate(f, a) for a in range(self.p))

        def raw_key(self, f):
            return sum(1 << (a * self.p + b) for a, b in enumerate(self.symbols(f)))

        def support(self, raw):
            assert 0 <= raw < 1 << self.d
            return tuple((i // self.p, i % self.p) for i in range(self.d) if raw >> i & 1)

        def address(self, points):
            j = len(points)
            if len({a for a, _ in points}) < j:
                return None
            offset = sum(comb(self.p, u) * self.p ** u for u in range(j))
            rank = sum(comb(a, u + 1) for u, (a, _) in enumerate(points))
            return offset + rank * self.p ** j + sum(b * self.p ** u
                                                     for u, (_, b) in enumerate(points))

        def descriptors(self, f):
            points = tuple(enumerate(self.symbols(f)))
            return (self.address(a) for j in range(self.h + 1) for a in combinations(points, j))

        def snapshot(self):
            return (tuple(sorted(self.slots.items())), tuple(sorted(self.owners.items())),
                    tuple(sorted((a, count, tuple(total)) for a, (count, total) in self.cells.items())))

        def remove(self, slot):
            if slot not in self.slots:
                return
            f, value = self.slots.pop(slot)
            assert self.owners.pop(f) == slot
            for a in self.descriptors(f):
                count, total = self.cells[a]
                count -= 1
                total = [x - y for x, y in zip(total, value)]
                if count:
                    self.cells[a] = count, total
                else:
                    assert total == [0, 0]
                    del self.cells[a]
            stats['deletions'] += 1

        def replace(self, slot, raw, value):
            assert 0 <= slot < self.n
            before = self.snapshot()
            points = self.support(raw) if 0 <= raw < 1 << self.d else ()
            valid = len(points) == self.p and len({a for a, _ in points}) == self.p
            f = next(affine_fiber(points[:self.k], self.k, self.field), None) if valid else None
            valid = f is not None and self.raw_key(f) == raw
            valid = valid and (f not in self.owners or self.owners[f] == slot)
            brute_valid = (raw in {self.raw_key(poly) for poly in product(range(self.p), repeat=self.k)}
                           if self.n <= 125 else f is not None and self.raw_key(f) == raw)
            brute_valid = brute_valid and all(i == slot or self.raw_key(poly) != raw
                                              for i, (poly, _) in self.slots.items())
            assert bool(valid) == bool(brute_valid)
            if not valid:
                assert self.snapshot() == before
                stats['rejections'] += 1
                return False
            self.remove(slot)
            self.slots[slot], self.owners[f] = (f, tuple(value)), slot
            for a in self.descriptors(f):
                count, total = self.cells.get(a, (0, [0, 0]))
                self.cells[a] = count + 1, [x + y for x, y in zip(total, value)]
            stats['admissions'] += 1
            return True

        def overwrite(self, slot, value):
            f, old = self.slots[slot]
            delta = [x - y for x, y in zip(value, old)]
            for a in self.descriptors(f):
                count, total = self.cells[a]
                self.cells[a] = count, [x + y for x, y in zip(total, delta)]
            self.slots[slot] = f, tuple(value)
            stats['overwrites'] += 1

        def check_state(self):
            expected = {}
            for f, value in self.slots.values():
                for a in self.descriptors(f):
                    count, total = expected.get(a, (0, [0, 0]))
                    expected[a] = count + 1, [x + y for x, y in zip(total, value)]
            assert expected == self.cells
            assert {f: slot for slot, (f, _) in self.slots.items()} == self.owners
            stats['live_zero_sums'] += sum(count > 0 and total == [0, 0]
                                          for count, total in self.cells.values())

        def read(self, raw):
            stats['reads'] += 1
            points = self.support(raw)
            size = len(points)
            numerator, denominator = [0, 0], 0
            if size > 2 * self.p:
                stats['empty_gates'] += 1
            elif size > 6 * (self.k - 1):
                stats['high_reads'] += 1
                q = interpolate(points, self.k, self.field)
                candidates = self.field.roots(q, self.k)
                self.last_neighbors = set()
                for f in candidates:
                    if f not in self.owners:
                        continue
                    weight = max(0, self.p + 1 - (raw ^ self.raw_key(f)).bit_count())
                    if weight:
                        self.last_neighbors.add(f)
                    value = self.slots[self.owners[f]][1]
                    denominator += weight
                    numerator = [x + weight * y for x, y in zip(numerator, value)]
                stats['high_two_neighbor_queries'] += len(self.last_neighbors) == 2
            else:
                stats['low_reads'] += 1
                differences = [max(0, 1 - size + 2 * s) for s in range(self.k)]
                coefficients = []
                while differences:
                    coefficients.append(differences[0])
                    differences = [y - x for x, y in zip(differences, differences[1:])]
                for j, c in enumerate(coefficients):
                    if c == 0:
                        continue
                    for subset in combinations(points, j):
                        if j <= self.h:
                            address = self.address(subset)
                            count, total = self.cells.get(address, (0, [0, 0]))
                        else:
                            count, total = 0, [0, 0]
                            for f in affine_fiber(subset, self.k, self.field):
                                stats['fiber_addresses_examined'] += 1
                                if f in self.owners:
                                    value = self.slots[self.owners[f]][1]
                                    count += 1
                                    total = [x + y for x, y in zip(total, value)]
                        denominator += c * count
                        numerator = [x + c * y for x, y in zip(numerator, total)]
                candidates = set()
                for subset in combinations(points, self.k):
                    # Inconsistent repeated rows are rejected by the solver.
                    for f in affine_fiber(subset, self.k, self.field):
                        assert len({a for a, _ in subset}) == self.k
                        candidates.add(f)
                for f in candidates & self.owners.keys():
                    overlap = (raw & self.raw_key(f)).bit_count()
                    poly = sum(c * comb(overlap, j) for j, c in enumerate(coefficients))
                    correction = max(0, 1 - size + 2 * overlap) - poly
                    value = self.slots[self.owners[f]][1]
                    denominator += correction
                    numerator = [x + correction * y for x, y in zip(numerator, value)]
                    stats['nonzero_corrections'] += correction != 0
            assert denominator >= 0
            self.last_totals = tuple(numerator), denominator
            assert denominator or numerator == [0, 0]
            return tuple(F(x, denominator) for x in numerator) if denominator else None

        def oracle(self, raw):
            weighted = [(value, max(0, self.p + 1 - (raw ^ self.raw_key(f)).bit_count()))
                        for f, value in self.slots.values()]
            denominator = sum(w for _, w in weighted)
            numerator = tuple(sum(value[j] * w for value, w in weighted) for j in range(2))
            self.oracle_totals = numerator, denominator
            return tuple(F(x, denominator) for x in numerator) if denominator else None

    primitive_roots = primitive_fibers = direct_addresses = 0
    for prime in (3, 5, 7):
        field = Field(prime)
        all_f = list(product(range(prime), repeat=2))
        # Nonconstant A and common factors test polynomial rather than
        # pointwise division. Include nonsquares, repeated roots and linear Q.
        for first, second in ((all_f[i], all_f[j]) for i in range(0, len(all_f), 3)
                              for j in range(0, len(all_f), 5)):
            factor = (1, 0, 1)
            q = [field.mul(field.mul(first, second), factor),
                 field.scale(field.mul(field.add(first, second), factor), -1), factor]
            found = field.roots(q, 2)
            brute = {f for f in all_f if not field.substitute(q, f)}
            assert found == brute == {first, second}
            primitive_roots += 1
        for q in ([(), (1, 1), ()], [(1,), (), (1,)], [(0, 1), (), (1,)],
                  [(1, 0, 1), (0, 1), (1,)], [(1,), (), ()]):
            q = [field.trim(a) for a in q]
            assert field.roots(q, 2) == {f for f in all_f if not field.substitute(q, f)}
            primitive_roots += 1
        for points in ((), ((0, 1),), ((0, 1), (0, 2)),
                       ((0, 1), (1, 2)), ((0, 1), (1, 2), (2, 0))):
            found = set(affine_fiber(points, 2, field))
            assert found == {f for f in all_f if all(field.evaluate(f, a) == b for a, b in points)}
            primitive_fibers += 1
        sample = CodedIndex(prime, 2, 1)
        addresses = [sample.address(tuple(zip(rows, values)))
                     for j in range(2) for rows in combinations(range(prime), j)
                     for values in product(range(prime), repeat=j)]
        assert sorted(addresses) == list(range(1 + prime * prime))
        direct_addresses += len(addresses)

    # Two nonzero contributions cancel in an actually shared cached cell.
    cancellation = CodedIndex(5, 2, 1)
    for slot, f, value in ((0, (0, 0), (11, -17)), (1, (0, 1), (-11, 17))):
        assert cancellation.replace(slot, cancellation.raw_key(f), value)
    cancellation.check_state()
    assert cancellation.cells[cancellation.address(((0, 0),))] == (2, [0, 0])
    cancellation.remove(1)
    cancellation.check_state()
    assert cancellation.cells[cancellation.address(((0, 0),))] == (1, [11, -17])
    for raw in (0, 1, cancellation.raw_key((0, 0))):
        assert cancellation.read(raw) == cancellation.oracle(raw)
        assert cancellation.last_totals == cancellation.oracle_totals
    cancellation.overwrite(0, ((1 << 257) + 17, -(1 << 263) + 19))
    for raw in (0, 1, cancellation.raw_key((0, 0))):
        assert cancellation.read(raw) == cancellation.oracle(raw)
        assert cancellation.last_totals == cancellation.oracle_totals

    empty = CodedIndex(7, 2, 1)
    for raw in (0, 1, empty.raw_key((0, 0)), (1 << empty.d) - 1):
        assert empty.read(raw) is None and empty.oracle(raw) is None
        assert empty.last_totals == empty.oracle_totals == ((0, 0), 0)

    address_index = CodedIndex(5, 4, 2)
    addresses = [address_index.address(tuple(zip(rows, values)))
                 for j in range(3) for rows in combinations(range(5), j)
                 for values in product(range(5), repeat=j)]
    assert sorted(addresses) == list(range(1 + 25 + 250))
    direct_addresses += len(addresses)

    address_index = CodedIndex(7, 6, 3)
    addresses = [address_index.address(tuple(zip(rows, values)))
                 for j in range(4) for rows in combinations(range(7), j)
                 for values in product(range(7), repeat=j)]
    assert sorted(addresses) == list(range(sum(comb(7, j) * 7 ** j for j in range(4))))
    direct_addresses += len(addresses)

    cases, sensitivity_reads, pure_reads, atomic_rejections = [], 0, 0, 0
    for prime, level, cached, complete in ((3, 2, 1, True), (5, 2, 1, True),
                                           (5, 3, 1, True), (5, 3, 2, True), (7, 2, 1, True),
                                           (7, 4, 3, False),
                                           (13, 3, 1, True), (7, 6, 3, False), (17, 5, 2, False),
                                           (29, 5, 2, False)):
        index = CodedIndex(prime, level, cached)
        field = index.field
        polynomials = list(product(range(prime), repeat=level)) if complete else [
            tuple([a, b] + [0] * (level - 2)) for a in range(4) for b in range(min(8, prime))]
        initial_stats = stats.copy()
        for slot, f in enumerate(polynomials):
            assert index.replace(slot, index.raw_key(f), (slot * slot + 1, -3 * slot - 2))
        index.check_state()

        queries = {0, (1 << index.d) - 1}
        if prime == 3:
            queries.update(range(1 << index.d))
        else:
            # Raw multi-hot sets, not a predecoded sequence of symbols.
            for size in range(1, min(2 * prime + 2, 16)):
                for offset in range(6):
                    coords = {(offset + 3 * j + j * j) % index.d for j in range(size)}
                    queries.add(sum(1 << i for i in coords))
            for f in polynomials[:4]:
                raw = index.raw_key(f)
                queries.add(raw)
                queries.update(raw ^ (1 << j) for j in range(0, index.d, max(1, index.d // 7)))
            first, second = polynomials[0], polynomials[1]
            union = index.raw_key(first) | index.raw_key(second)
            queries.add(union)
            # Test both sides of the strict interpolation threshold, even
            # when the low branch has a large fixed-k constant.
            support = index.support(union)
            for size in (6 * (level - 1), 6 * (level - 1) + 1, 2 * prime, 2 * prime + 1):
                if size <= index.d:
                    pairs = list(support[:size])
                    for i in range(index.d):
                        point = divmod(i, prime)
                        if len(pairs) >= size:
                            break
                        if point not in pairs:
                            pairs.append(point)
                    queries.add(sum(1 << (a * prime + b) for a, b in pairs))
        for raw in sorted(queries):
            assert index.read(raw) == index.oracle(raw)
            assert index.last_totals == index.oracle_totals

        # Current arbitrary payloads, signed cancellation, latest overwrite,
        # deletion/reuse, duplicate rejection and malformed raw admission.
        index.overwrite(0, (11, -17))
        index.overwrite(1, (-11, 17))
        index.check_state()
        invalid = [0, index.raw_key(polynomials[0]) ^ 1,
                   index.raw_key(polynomials[1]), 1 << index.d]
        # One-hot but not a degree-<k polynomial: change one evaluation.
        if level < prime:
            invalid.append(index.raw_key(polynomials[0]) ^ 1 ^ 2)
        for raw in invalid:
            before = index.snapshot()
            assert not index.replace(0, raw, (23, 31))
            assert index.snapshot() == before
            atomic_rejections += 1
        index.remove(1)
        assert index.replace(1, index.raw_key(polynomials[1]), (0, 0))
        index.check_state()
        dynamic_queries = sorted(queries)[:12] + [index.raw_key(polynomials[0]),
                                                index.raw_key(polynomials[0]) | index.raw_key(polynomials[1])]
        for raw in dynamic_queries:
            assert index.read(raw) == index.oracle(raw)
            assert index.last_totals == index.oracle_totals
        if prime > 6 * (level - 1):
            raw = index.raw_key(polynomials[0]) | index.raw_key(polynomials[1])
            assert index.read(raw) == index.oracle(raw)
            assert index.last_totals == index.oracle_totals
            assert len(index.last_neighbors) == 2
            index.overwrite(1, (-37, 41))
            assert index.read(raw) == index.oracle(raw)
            assert index.last_totals == index.oracle_totals
            index.remove(1)
            assert index.read(raw) == index.oracle(raw)
            assert index.last_totals == index.oracle_totals
            assert len(index.last_neighbors) == 1
            assert index.replace(1, index.raw_key(polynomials[1]), (-37, 41))
        if prime > 2 * level - 2:
            for slot in range(min(4, len(polynomials))):
                f, value = index.slots[slot]
                key = index.raw_key(f)
                for radius in range(min(2, prime - 2 * level + 1) + 1):
                    coords = range(index.d) if prime <= 5 else range(0, index.d, max(1, index.d // 8))
                    for flips in combinations(coords, radius):
                        raw = key ^ sum(1 << i for i in flips)
                        assert index.read(raw) == tuple(F(x) for x in value) == index.oracle(raw)
                        assert index.last_totals == index.oracle_totals
                        pure_reads += 1
        # Exact every-bit sensitivity at the zero query. This witness is
        # complete-code only; pairwise-independent evaluations prove uniformity.
        if complete and prime <= 7:
            def label(y):
                return 1 if y < prime - 1 else -(prime - 1)
            for slot, (f, _) in list(index.slots.items()):
                value = sum(label(y) for y in index.symbols(f))
                index.overwrite(slot, (value, -value))
            assert index.read(0) == (F(0), F(0)) == index.oracle(0)
            assert index.last_totals == index.oracle_totals
            sensitivity_reads += 1
            for a in range(prime):
                for y in range(prime):
                    raw = 1 << (a * prime + y)
                    assert index.read(raw) == (F(label(y)), F(-label(y))) == index.oracle(raw)
                    assert index.last_totals == index.oracle_totals
                    sensitivity_reads += 1
        index.check_state()
        # Replacing a key by a different valid raw codeword updates both
        # inverse ownership and every cached aggregate, including row changes.
        old_first, old_second = index.slots[0], index.slots[2]
        index.remove(2)
        assert index.replace(0, index.raw_key(old_second[0]), (101, -103))
        assert index.replace(2, index.raw_key(old_first[0]), (-107, 109))
        index.check_state()
        for raw in (0, 1, index.raw_key(old_first[0]), index.raw_key(old_second[0])):
            assert index.read(raw) == index.oracle(raw)
            assert index.last_totals == index.oracle_totals
        cases.append({"prime": prime, "polynomial_dimension": level, "cached_orders_through": cached,
                      "capacity": prime ** level, "raw_dimension": index.d,
                      "linear_space_parameter_condition": 2 * cached <= level,
                      "complete_codebook_values_installed": complete,
                      "occupied_values": len(index.slots),
                      "cache_cell_capacity": sum(comb(prime, j) * prime ** j for j in range(cached + 1)),
                      "cells_per_value_edit": sum(comb(prime, j) for j in range(cached + 1)),
                      "query_controls": {name: stats[name] - initial_stats[name]
                                         for name in ('reads', 'low_reads', 'high_reads', 'empty_gates',
                                                      'fiber_addresses_examined', 'nonzero_corrections',
                                                      'high_two_neighbor_queries')}})
    assert stats['low_reads'] and stats['high_reads'] and stats['empty_gates']
    assert stats['nonzero_corrections'] and stats['fiber_addresses_examined']
    assert stats['high_two_neighbor_queries'] and stats['live_zero_sums']
    growing_parameters = []
    for prime in (127, 257, 521, 1031, 2053, 4093):
        level = prime // (4 * (prime - 1).bit_length())
        capacity = prime ** level
        edit_cells = sum(comb(prime, j) for j in range(level))
        cache_cells = sum(comb(prime, j) * prime ** j for j in range(level))
        query_subsets = sum(comb(6 * (level - 1), j) for j in range(level + 1))
        assert cache_cells * prime <= capacity * edit_cells
        assert query_subsets <= 2 ** (6 * (level - 1))
        assert capacity ** 4 <= 2 ** prime
        growing_parameters.append({"prime": prime, "polynomial_dimension": level,
                                   "capacity_bits": capacity.bit_length(),
                                   "cache_cell_capacity_bits": cache_cells.bit_length(),
                                   "cells_per_edit_bits": edit_cells.bit_length(),
                                   "low_query_subset_bound_bits": query_subsets.bit_length()})

    # Published ranked up/down zeta + Pascal inversion, independently
    # compared to our truncated moments and unique-owner tail compression.
    transform_reads = 0
    for dimension, weight, level, supports in (
            (6, 2, 1, [(0, 1), (2, 3), (4, 5)]),
            (7, 3, 2, [(0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5),
                       (1, 4, 6), (2, 3, 6), (2, 4, 5)]),
            (9, 4, 2, [(0, 1, 2, 3), (0, 4, 5, 6), (1, 4, 7, 8)])):
        assert all(len(set(x) & set(y)) < level for x, y in combinations(supports, 2))
        values = [(1, i * i - 3, (-1) ** i * (i + 2)) for i in range(len(supports))]
        up, owners = {}, {}
        for owner, support in enumerate(supports):
            for j in range(weight + 1):
                for subset in combinations(support, j):
                    previous = up.get(subset, (0, 0, 0))
                    up[subset] = tuple(x + y for x, y in zip(previous, values[owner]))
                    if j == level:
                        assert subset not in owners
                        owners[subset] = owner
        for raw in range(1 << dimension):
            query = tuple(i for i in range(dimension) if raw >> i & 1)
            overlaps = [len(set(query) & set(support)) for support in supports]
            moments = [tuple(sum(up.get(subset, (0, 0, 0))[v]
                                 for subset in combinations(query, j)) for v in range(3))
                       for j in range(weight + 1)]
            distribution = [tuple(sum((-1) ** (i + j) * comb(i, j) * moments[i][v]
                                      for i in range(j, weight + 1)) for v in range(3))
                            for j in range(weight + 1)]
            assert distribution == [tuple(sum(value[v] for value, s in zip(values, overlaps) if s == j)
                                          for v in range(3)) for j in range(weight + 1)]
            candidates = {owners[a] for a in combinations(query, level) if a in owners}
            assert candidates == {i for i, s in enumerate(overlaps) if s >= level}
            for j in range(level, weight + 1):
                assert moments[j] == tuple(sum(comb(overlaps[i], j) * values[i][v] for i in candidates)
                                           for v in range(3))
            for radius in (0, weight, 2 * weight):
                g = [max(0, radius + 1 - weight - len(query) + 2 * s) for s in range(weight + 1)]
                differences, coefficients = g[:], []
                while differences:
                    coefficients.append(differences[0])
                    differences = [y - x for x, y in zip(differences, differences[1:])]
                full = tuple(sum(g[j] * distribution[j][v] for j in range(weight + 1)) for v in range(3))
                assert full == tuple(sum(coefficients[j] * moments[j][v] for j in range(weight + 1))
                                     for v in range(3))
                compressed = tuple(
                    sum(coefficients[j] * moments[j][v] for j in range(level))
                    + sum((g[overlaps[i]] - sum(coefficients[j] * comb(overlaps[i], j)
                                               for j in range(level))) * values[i][v] for i in candidates)
                    for v in range(3))
                assert full == compressed
                transform_reads += 1
    return {"cases": cases, "exact_controls": stats,
            "polynomial_root_primitive_cases": primitive_roots,
            "affine_fiber_primitive_cases": primitive_fibers,
            "distinct_direct_cell_addresses": direct_addresses,
            "all_bit_sensitivity_reads": sensitivity_reads,
            "pure_current_value_noise_reads": pure_reads,
            "atomic_invalid_or_duplicate_rejections": atomic_rejections,
            "large_exact_payload_guard_control_bits": 264,
            "empty_state_raw_queries": 4,
            "published_intersection_transform_equivalence_reads": transform_reads,
            "intersection_transform_scope": "On three non-polynomial low-intersection families, exhaust all raw queries and three radii, checking the published full ranked up/down zeta moments and inverse Pascal distribution against direct intersection classes, then our low-order truncation plus exact unique-owner tail. This proves the algebraic equivalence, not that the published full-down-closure theorem already states the compressed dynamic bound or that the elementary compression is a novel contribution.",
            "growing_dimension_parameter_bit_lengths": growing_parameters,
            "general_space_bound": "For any 1<=h<k, cache C_h=sum_{j<=h}binom(p,j)*p^j cells and keep n=p^k independent value slots. State O_k,h((n+C_h)*(v+1)+p^2). The linear-space simplification additionally needs 2h<=k at fixed k. The h=k-1 controls also exercise the general, potentially superlinear fixed-k tradeoff.",
            "uniform_bound": "For fixed k,h with 1<=h<k and 2h<=k, odd prime p>=k, n=p^k, raw dimension p^2 and radius p: state O_k,h(n*(v+1)), exact reads O_k,h(p^2+p^(k-h-1)*(v+1)), overwrites/deletes O_k,h(p^h*(v+1)), raw admitted rekeys O_k,h(p^2+p^h*(v+1)). Field/precision bit costs are additional; guard bits and initialization are charged in the note. For k=5,h=2 scalar reads/rekeys are n^(2/5) up to word/precision factors, with linear state.",
            "growing_dimension_scope": "With k=floor(p/(4ceil(log2 p))) and h=k-1, keep exact C_h cells rather than the loose p^(2h) bound. This yields n^(1+o(1)) state/preprocessing, n^o(1) updates, and 2^O(k)*poly(p,k)*(v+1) complete raw reads. No implicit constant is treated as uniform in growing k; the note expands it. This is a codebook-specific comparator in the preceding generic candidate's asymptotic regime, not a novelty determination.",
            "scope": "Exact sparse-map correctness control, not a direct-array implementation, hardware benchmark or novelty proof. All raw-query numerators/denominators are checked against independent Hamming scans. Small complete codebooks and sparsely occupied larger fixed codebooks are distinguished; larger capacities are conceptual. The construction uses established soft interpolation and affine evaluation constraints with partial cached sums. It restricts key admission to a known polynomial codebook, so does not settle the priority or usefulness of the generic low-intersection candidate. No neural training or strongest-baseline separation is claimed."}


def weighted_exponential_controls():
    """Full signed-integer-query softmax versus independent dyadic sums."""
    from math import comb

    class Index:
        # Correctness map, not an arena/word-RAM benchmark. The written
        # bound uses the preserved finite-capacity trie and explicit handles.
        def __init__(self, capacity, dimension, weight, level, budget, value_bits=32, load_limit=None):
            self.n, self.d, self.w, self.k = capacity, dimension, weight, level
            self.budget, self.value_bits = budget, value_bits
            self.keys, self.values = {}, {}
            self.sums, self.owners, self.handles = {}, {}, {}
            self.postings, self.total_value = {}, [0, 0]
            self.load_limit = load_limit
            self.f = sum(comb(weight, j) for j in range(level + 1))
            self.qmax = sum(comb(min(dimension, budget), j) for j in range(level + 1))
            self.guard = value_bits + budget + (capacity + 1).bit_length() + self.f.bit_length() + self.qmax.bit_length() + 4
            self.stats = dict(reads=0, subset_lookups=0, candidate_visits=0,
                              negative_coefficients=0, negative_corrections=0,
                              no_candidate_reads=0, bulk_no_candidate_reads=0,
                              overwrites=0, deletes=0, admitted=0, rejected=0,
                              live_zero_sum_observations=0, max_intermediate_bits=0)
            self.stats.update(posting_reads=0, posting_hits=0, posting_value_visits=0,
                              max_coordinate_load=0, load_only_rejections=0)

        def support(self, key):
            return tuple(i for i in range(self.d) if key >> i & 1)

        def snapshot(self):
            return (tuple(sorted(self.keys.items())), tuple(sorted(self.values.items())),
                    tuple(sorted((a, count, tuple(total)) for a, (count, total) in self.sums.items())),
                    tuple(sorted(self.owners.items())), tuple(sorted(self.handles.items())),
                    tuple(sorted((a, tuple(sorted(slots))) for a, slots in self.postings.items())),
                    tuple(self.total_value))

        def remove(self, slot):
            if slot not in self.keys:
                return
            support = self.support(self.keys[slot])
            for a in self.handles.pop(slot):
                count, total = self.sums[a]
                count -= 1
                total = [x - y for x, y in zip(total, self.values[slot])]
                if count:
                    self.sums[a] = count, total
                else:
                    assert total == [0, 0]
                    del self.sums[a]
            for a in combinations(support, self.k):
                assert self.owners.pop(a) == slot
            del self.keys[slot]
            for a in support:
                self.postings[a].remove(slot)
                if not self.postings[a]:
                    del self.postings[a]
            self.total_value = [x - y for x, y in zip(self.total_value, self.values[slot])]
            del self.values[slot]
            self.stats['deletes'] += 1

        def replace(self, slot, key, value):
            assert 0 <= slot < self.n and all(abs(x) < 1 << self.value_bits for x in value)
            before = self.snapshot()
            support = self.support(key) if 0 <= key < 1 << self.d else ()
            valid = len(support) == self.w
            valid = valid and all(a not in self.owners or self.owners[a] == slot
                                  for a in combinations(support, self.k))
            load_valid = (self.load_limit is None or all(
                len(self.postings.get(a, ())) - (slot in self.postings.get(a, ())) < self.load_limit
                for a in support))
            self.stats['load_only_rejections'] += valid and not load_valid
            valid = valid and load_valid
            brute = (0 <= key < 1 << self.d and key.bit_count() == self.w
                     and all(i == slot or (other & key).bit_count() < self.k for i, other in self.keys.items()))
            brute = brute and (self.load_limit is None or all(
                sum(i != slot and (other >> a & 1) for i, other in self.keys.items()) < self.load_limit
                for a in support))
            assert valid == brute
            if not valid:
                assert self.snapshot() == before
                self.stats['rejected'] += 1
                return False
            self.remove(slot)
            self.keys[slot], self.values[slot] = key, tuple(value)
            self.total_value = [x + y for x, y in zip(self.total_value, value)]
            for a in support:
                self.postings.setdefault(a, set()).add(slot)
                self.stats['max_coordinate_load'] = max(self.stats['max_coordinate_load'], len(self.postings[a]))
            handles = tuple(a for j in range(self.k) for a in combinations(support, j))
            self.handles[slot] = handles
            for a in handles:
                count, total = self.sums.get(a, (0, [0, 0]))
                self.sums[a] = count + 1, [x + y for x, y in zip(total, value)]
            for a in combinations(support, self.k):
                assert a not in self.owners
                self.owners[a] = slot
            self.stats['admitted'] += 1
            return True

        def overwrite(self, slot, value):
            assert slot in self.keys and all(abs(x) < 1 << self.value_bits for x in value)
            delta = [x - y for x, y in zip(value, self.values[slot])]
            for a in self.handles[slot]:
                count, total = self.sums[a]
                self.sums[a] = count, [x + y for x, y in zip(total, delta)]
            self.values[slot] = tuple(value)
            self.total_value = [x + y for x, y in zip(self.total_value, delta)]
            self.stats['overwrites'] += 1

        def check_state(self):
            expected, owners = {}, {}
            for slot, key in self.keys.items():
                support = self.support(key)
                assert self.handles[slot] == tuple(a for j in range(self.k) for a in combinations(support, j))
                for a in self.handles[slot]:
                    count, total = expected.get(a, (0, [0, 0]))
                    expected[a] = count + 1, [x + y for x, y in zip(total, self.values[slot])]
                for a in combinations(support, self.k):
                    assert a not in owners
                    owners[a] = slot
            assert self.sums == expected and self.owners == owners
            postings = {}
            for slot, key in self.keys.items():
                for a in self.support(key):
                    postings.setdefault(a, set()).add(slot)
            assert postings == self.postings
            assert self.total_value == [sum(value[j] for value in self.values.values()) for j in range(2)]
            assert self.load_limit is None or all(len(slots) <= self.load_limit for slots in self.postings.values())
            self.stats['live_zero_sum_observations'] += sum(count > 0 and total == [0, 0]
                                                          for count, total in self.sums.values())
            assert len(self.sums) + len(self.owners) <= self.n * self.f

        def bounded(self, *numbers):
            bits = max((abs(x).bit_length() + 1 for x in numbers), default=1)
            self.stats['max_intermediate_bits'] = max(bits, self.stats['max_intermediate_bits'])
            assert bits <= self.guard

        def read(self, raw):
            if len(raw) != self.d or any(type(x) is not int for x in raw) or sum(abs(x) for x in raw) > self.budget:
                raise ValueError('query outside the stated integer norm budget')
            self.stats['reads'] += 1
            support = tuple(i for i in range(self.d) if raw[i])
            negative_sum = sum(-x for x in raw if x < 0)
            factor = {a: (1 << abs(raw[a])) - 1 for a in support}
            numerator, denominator = [0, 0], 0
            for j in range(self.k):
                for a in combinations(support, j):
                    negative = sum(-raw[i] for i in a if raw[i] < 0)
                    coefficient = 1 << (negative_sum - negative)
                    for i in a:
                        coefficient *= factor[i] * (-1 if raw[i] < 0 else 1)
                    self.stats['negative_coefficients'] += coefficient < 0
                    self.stats['subset_lookups'] += 1
                    count, total = self.sums.get(a, (0, [0, 0]))
                    denominator += coefficient * count
                    numerator = [x + coefficient * y for x, y in zip(numerator, total)]
                    self.bounded(coefficient, denominator, *numerator)
            candidates = {self.owners[a] for a in combinations(support, self.k) if a in self.owners}
            self.stats['subset_lookups'] += comb(len(support), self.k)
            for slot in candidates:
                overlap = [a for a in self.support(self.keys[slot]) if raw[a]]
                # E_j = 2^negative_sum times the j-th elementary symmetric
                # polynomial. Dividing by a new negative factor's power
                # of two is exact before it is introduced.
                elementary = [1 << negative_sum] + [0] * (self.k - 1)
                for a in overlap:
                    shift = max(0, -raw[a])
                    scale = 1 << shift
                    for j in range(self.k - 1, 0, -1):
                        assert elementary[j - 1] % scale == 0
                        elementary[j] += (elementary[j - 1] // scale) * factor[a] * (-1 if raw[a] < 0 else 1)
                        self.bounded(elementary[j])
                exponent = negative_sum + sum(raw[a] for a in overlap)
                assert 0 <= exponent <= self.budget
                correction = (1 << exponent) - sum(elementary)
                denominator += correction
                numerator = [x + correction * y for x, y in zip(numerator, self.values[slot])]
                self.stats['candidate_visits'] += 1
                self.stats['negative_corrections'] += correction < 0
                self.bounded(correction, denominator, *numerator)
            assert denominator > 0 if self.keys else denominator == 0 and numerator == [0, 0]
            self.last_totals = tuple(numerator), denominator, negative_sum
            self.last_candidates = candidates
            if not candidates:
                self.stats['no_candidate_reads'] += 1
                self.stats['bulk_no_candidate_reads'] += len(self.keys) > 1
            return tuple(F(x, denominator) for x in numerator) if denominator else None

        def oracle(self, raw):
            # Independent dot products and rational exponentiation, with
            # no subset expansion, truncation, owners or elementary DP.
            weighted = [(slot, F(2) ** sum(raw[a] for a in self.support(key)))
                        for slot, key in self.keys.items()]
            denominator = sum((weight for _, weight in weighted), F(0))
            numerator = tuple(sum((weight * self.values[slot][j] for slot, weight in weighted), F(0))
                              for j in range(2))
            scale = 1 << sum(-x for x in raw if x < 0)
            scaled = tuple(x * scale for x in numerator), denominator * scale
            assert all(x.denominator == 1 for x in scaled[0]) and scaled[1].denominator == 1
            self.oracle_totals = tuple(int(x) for x in scaled[0]), int(scaled[1]), scale.bit_length() - 1
            return tuple(x / denominator for x in numerator) if denominator else None

        def posting_read(self, raw):
            # Standalone h=0 baseline: total plus exact corrections at
            # keys sharing >=1 coordinate. This control shares values with
            # the subset index; no combined-state performance is claimed.
            self.stats['posting_reads'] += 1
            scores = {}
            for a, value in enumerate(raw):
                if value:
                    for slot in self.postings.get(a, ()):
                        scores[slot] = scores.get(slot, 0) + value
                        self.stats['posting_hits'] += 1
            negative = sum(-value for value in raw if value < 0)
            base = 1 << negative
            numerator = [base * x for x in self.total_value]
            denominator = base * len(self.keys)
            for slot, score in scores.items():
                assert 0 <= negative + score <= self.budget
                correction = (1 << (negative + score)) - base
                numerator = [x + correction * y for x, y in zip(numerator, self.values[slot])]
                denominator += correction
                assert max(abs(x).bit_length() + 1 for x in [denominator, *numerator]) <= self.value_bits + self.budget + (self.n + 1).bit_length() + 4
                self.stats['posting_value_visits'] += 1
            self.posting_totals = tuple(numerator), denominator, negative
            return tuple(F(x, denominator) for x in numerator) if denominator else None

    def verify(index, raw):
        actual = index.read(raw)
        assert actual == index.oracle(raw)
        assert index.last_totals == index.oracle_totals
        assert index.posting_read(raw) == actual
        assert index.posting_totals == index.oracle_totals
        return actual

    cases, outside_admissions, rejected_queries = [], 0, 0
    families = [(6, 2, 1, [(0, 1), (2, 3), (4, 5)]),
                (7, 3, 2, [(0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5),
                           (1, 4, 6), (2, 3, 6), (2, 4, 5)]),
                (9, 4, 2, [(0, 1, 2, 3), (0, 4, 5, 6), (1, 4, 7, 8)])]
    for prime, level in ((5, 2), (5, 3)):
        supports = [tuple(a * prime + sum(c * pow(a, j, prime) for j, c in enumerate(f)) % prime
                          for a in range(prime)) for f in product(range(prime), repeat=level)]
        families.append((prime * prime, prime, level, supports))
    for dimension, weight, level, supports in families:
        index = Index(len(supports) + 3, dimension, weight, level, 2 * weight)
        initial = [sum(1 << i for i in support) for support in supports]
        for slot, key in enumerate(initial):
            assert index.replace(slot, key, (slot * slot - 3, (-1) ** slot * (slot + 2)))
        index.check_state()
        queries = []
        if dimension <= 9:
            queries = [q for q in product((-1, 0, 1), repeat=dimension) if sum(abs(x) for x in q) <= index.budget]
        else:
            queries.append((0,) * dimension)
            for size in range(1, 2 * weight + 1):
                for offset in range(5):
                    q = [0] * dimension
                    for j in range(size):
                        q[(offset + 7 * j + j * j) % dimension] = (-1) ** j
                    queries.append(tuple(q))
        # Nonunit amplitudes, unequal signed coordinates, and full supports.
        for a in range(min(dimension, 7)):
            for b in range(a + 1, min(dimension, 7)):
                for first, second in ((2, 1), (1, 2), (-2, 1), (2, -1), (-2, -1)):
                    q = [0] * dimension
                    q[a], q[b] = first, second
                    queries.append(tuple(q))
        for q in queries:
            verify(index, q)
        for raw in ((), (index.budget + 1,) + (0,) * (dimension - 1), (0.5,) + (0,) * (dimension - 1)):
            before = index.snapshot()
            try:
                index.read(raw)
                assert False
            except ValueError:
                pass
            assert index.snapshot() == before
            rejected_queries += 1
        index.overwrite(0, (11, -17))
        index.overwrite(1, (-11, 17))
        for q in queries[:12] + queries[-12:] + [(0,) * dimension]:
            verify(index, q)
        for proposal in (0, initial[0] ^ 1, initial[1], 1 << dimension):
            before = index.snapshot()
            assert not index.replace(0, proposal, (23, 31))
            assert index.snapshot() == before
        index.remove(0)
        index.remove(1)
        # Generic certified admission, not just a polynomial codebook.
        for support in combinations(range(dimension), weight):
            proposal = sum(1 << a for a in support)
            if proposal in initial:
                continue
            if index.replace(0, proposal, (19, -23)):
                outside_admissions += 1
                break
            if index.stats['rejected'] > 80:
                break
        index.check_state()
        for q in queries[:12] + queries[-12:] + [(0,) * dimension]:
            verify(index, q)
        for slot in list(index.keys):
            index.remove(slot)
        index.check_state()
        verify(index, (0,) * dimension)
        assert index.replace(0, initial[0], (0, 0))
        index.check_state()
        for q in queries[-12:]:
            verify(index, q)
        cases.append({'dimension': dimension, 'weight': weight, 'intersection_threshold': level,
                      'capacity': index.n, 'query_l1_budget': index.budget,
                      'preselected_guard_bits': index.guard, 'controls': index.stats.copy()})

    # Same selected coordinates, different amplitudes: no unweighted
    # intersection-count score table can implement both rows of this witness.
    witness = Index(3, 6, 2, 1, 4)
    for slot, key in enumerate((3, 12, 48)):
        assert witness.replace(slot, key, (slot, -slot))
    first, second = (2, 0, 1, 0, 0, 0), (1, 0, 2, 0, 0, 0)
    first_result = verify(witness, first)
    second_result = verify(witness, second)
    assert first_result != second_result

    guard = Index(3, 6, 2, 2, 512, 264)
    for slot, key in enumerate((3, 12, 48)):
        assert guard.replace(slot, key, (((-1) ** slot) * ((1 << 257) + slot), -(1 << 263) + 19 + slot))
    for q in ((128, -127, 0, 0, 0, 0), (-128, -127, 0, 0, 0, 0),
              (128, -1, 0, 0, 0, 0), (-1, 0, 2, 0, -3, 0)):
        verify(guard, q)
    guard.check_state()

    # Uniform grid retention and dense influence, with every active value
    # included rather than masked out. Radius one is exhaustive at four keys.
    prime, level = 13, 2
    retention = Index(prime ** level, prime * prime, prime, level, 2 * prime)
    polys = list(product(range(prime), repeat=level))
    for slot, f in enumerate(polys):
        key = sum(1 << (a * prime + (f[0] + f[1] * a) % prime) for a in range(prime))
        assert retention.replace(slot, key, (slot % 7 - 3, (3 * slot) % 7 - 3))
    retention.check_state()
    noise_reads = 0
    ratio_bound = F(retention.n - 1, 2 ** (prime - level + 1 - 1))
    assert 2 * 3 * ratio_bound < F(1, 2)
    for slot in (0, 1, 47, 168):
        q = tuple(int(retention.keys[slot] >> a & 1) for a in range(retention.d))
        for a in range(-1, retention.d):
            for delta in ((0,) if a == -1 else (-1, 1)):
                raw = list(q)
                if a >= 0:
                    raw[a] += delta
                result = verify(retention, tuple(raw))
                rounded = tuple((2 * x.numerator + x.denominator) // (2 * x.denominator)
                                for x in result)
                assert rounded == retention.values[slot]
                noise_reads += 1
        retention.overwrite(slot, tuple(-x for x in retention.values[slot]))
        assert tuple(round(x) for x in verify(retention, q)) == retention.values[slot]
    q = tuple(int(retention.keys[0] >> a & 1) for a in range(retention.d))
    old_output, old_value = verify(retention, q), retention.values[43]
    new_value = old_value[0] + 1, old_value[1] - 1
    retention.overwrite(43, new_value)
    new_output = verify(retention, q)
    weight = F(2) ** sum(q[a] for a in retention.support(retention.keys[43]))
    unscaled_den = F(retention.last_totals[1], 1 << retention.last_totals[2])
    assert tuple(x - y for x, y in zip(new_output, old_output)) == (weight / unscaled_den, -weight / unscaled_den)
    assert 43 not in retention.last_candidates  # A dense low-score value still changes the answer.
    row_errors = []
    for key in retention.keys.values():
        weights = [1 << (key & other).bit_count() for other in retention.keys.values()]
        total = sum(weights)
        error = F(2 * (total - (1 << prime)), total)
        assert error <= F(2 * (retention.n - 1), 2 ** (prime - level + 1)) < 1
        row_errors.append(error)
    rho = max(row_errors)
    load_control = Index(3, 6, 2, 2, 4, load_limit=1)
    for slot, key in enumerate((3, 12, 48)):
        assert load_control.replace(slot, key, (slot + 1, -slot - 2))
    before = load_control.snapshot()
    assert not load_control.replace(0, 5, (17, -19))
    assert load_control.snapshot() == before
    assert load_control.stats['load_only_rejections'] == 1
    load_control.remove(1)
    assert load_control.replace(0, 5, (17, -19))
    load_control.check_state()
    for raw in ((0,) * 6, (2, -1, 1, 0, 0, 0), (0, 0, -2, 0, 1, 0)):
        verify(load_control, raw)
    random_admission_states = random_admission_proposals = 0
    # Exhaustive hypergeometric union bounds in genuinely generic states.
    for dimension, weight, level, supports in families[:3]:
        keys = [sum(1 << a for a in support) for support in supports]
        proposals = [sum(1 << a for a in support) for support in combinations(range(dimension), weight)]
        for bits in range(1 << len(keys)):
            active = [key for i, key in enumerate(keys) if bits >> i & 1]
            rejected = sum(any((proposal & key).bit_count() >= level for key in active) for proposal in proposals)
            union_bound = F(len(active) * comb(weight, level) ** 2, comb(dimension, level))
            assert F(rejected, len(proposals)) <= min(F(1), union_bound)
            random_admission_states += 1
            random_admission_proposals += len(proposals)
    parameter_examples = []
    for m in (16, 32, 64, 128, 256, 512):
        capacity, weight = 1 << m, 8 * m
        dimension, logarithm = weight * weight, m.bit_length()
        level = (6 * m + logarithm - 1) // logarithm
        assert 1 <= level <= weight // 2
        pair_num, pair_den = comb(weight, level) ** 2, comb(dimension, level)
        assert capacity ** 3 * pair_num <= pair_den  # rejection <= n^-2.
        noise = weight // 4
        gap = weight - level + 1 - noise
        assert 4 * (capacity - 1) * 2 ** (2 * m) < 2 ** gap
        features = sum(comb(weight, j) for j in range(level + 1))
        query_subsets = sum(comb(2 * weight, j) for j in range(level + 1))
        parameter_examples.append({'log2_capacity': m, 'weight': weight, 'dimension': dimension,
                                   'threshold': level, 'query_l1_budget': 2 * weight,
                                   'grid_precision_bits': 2 * m, 'retained_l1_noise_radius': noise,
                                   'per_proposal_rejection_at_most_inverse_capacity_squared': True,
                                   'grid_retention_strict_integer_bound_passed': True,
                                   'key_subset_count_bits': features.bit_length(),
                                   'max_query_subset_count_bits': query_subsets.bit_length(),
                                   'key_subset_count_exceeds_capacity': features > capacity})
    for index in (witness, guard, retention):
        index.check_state()
    all_stats = [case['controls'] for case in cases] + [witness.stats, guard.stats, retention.stats, load_control.stats]
    assert outside_admissions > 0
    assert sum(s['negative_corrections'] for s in all_stats) > 0
    assert sum(s['negative_coefficients'] for s in all_stats) > 0
    assert sum(s['live_zero_sum_observations'] for s in all_stats) > 0
    return {'cases': cases, 'total_exact_raw_reads': sum(s['reads'] for s in all_stats),
            'total_candidate_visits': sum(s['candidate_visits'] for s in all_stats),
            'negative_subset_coefficients': sum(s['negative_coefficients'] for s in all_stats),
            'negative_owner_corrections': sum(s['negative_corrections'] for s in all_stats),
            'bulk_reads_without_owner_candidates': sum(s['bulk_no_candidate_reads'] for s in all_stats),
            'generic_outside_initial_family_admissions': outside_admissions,
            'atomic_invalid_queries_rejected': rejected_queries,
            'live_zero_sum_observations': sum(s['live_zero_sum_observations'] for s in all_stats),
            'exact_inverted_index_comparator': {'reads': sum(s['posting_reads'] for s in all_stats),
                                                'coordinate_posting_hits': sum(s['posting_hits'] for s in all_stats),
                                                'corrected_value_visits': sum(s['posting_value_visits'] for s in all_stats),
                                                'atomic_load_only_rejections': load_control.stats['load_only_rejections'],
                                                'scope': 'Maintained coordinate postings plus total value sum give every full raw-query numerator/denominator independently of subset truncation. The verifier stores both comparison states and does not benchmark either. A standalone linked-posting implementation and its load-capped costs are accounted in the note; the prototype uses exact sets.'},
            'same_support_amplitude_witness': {'first_query': first, 'second_query': second,
                                              'first_output': [str(x) for x in first_result],
                                              'second_output': [str(x) for x in second_result]},
            'large_dyadic_precision_control': {'budget': guard.budget, 'input_value_bits': guard.value_bits,
                                               'preselected_guard_bits': guard.guard,
                                               'observed_max_intermediate_bits': guard.stats['max_intermediate_bits']},
            'grid_retention': {'capacity': retention.n, 'weight': retention.w, 'threshold': retention.k,
                               'query_l1_noise_radius': 1, 'value_integer_bound': 3,
                               'off_target_ratio_upper_bound': str(ratio_bound),
                               'conservative_value_error_bound': str(6 * ratio_bound),
                               'exhaustive_noise_queries_at_four_keys': noise_reads,
                               'independently_overwritten_target_checks': 4,
                               'noncandidate_value_exactly_changes_dense_output': True,
                               'all_exact_query_at_key_row_bounds': len(row_errors),
                               'normalized_operator_distance_from_identity': str(rho),
                               'rank_deficient_uniform_value_error_lower_bound': str(1 - rho)},
            'random_admission_controls': {'exhaustive_current_states': random_admission_states,
                                           'uniform_proposals_examined': random_admission_proposals,
                                           'bound': 'At any current state with <=n keys, a fresh independent uniform weight-w proposal is rejected with probability <=min(1,n*binom(w,k)^2/binom(d,k)). This is conditional per proposal, not an event of no future rejection or a worst-case resampling bound.'},
            'uniform_parameter_integer_examples': parameter_examples,
            'finite_size_warning': 'Every shown key-subset count exceeds capacity, despite the eventual subpolynomial asymptotic bound. These examples certify inequalities and expose expensive finite constants; they do not certify useful ML latency.',
            'scope': 'Exact tuple-map controls for full base-two exponential attention under a signed-integer raw-query l1 budget and certified low-intersection binary key admission. Independent direct dyadic dot-product sums check normalized output and scaled numerator/denominator, including negative coefficients/corrections, updates, generic admission, cancellation, deletion/reuse, malformed inputs, large precision and noisy grid retrieval. Python arbitrary integers are bounded against preselected guard bits; this is not a finite-word arena implementation or performance/novelty certificate. The 169-value retention example is fully installed. Other capacities are explicit. No learned encoder, training, strongest-baseline advantage or publication priority is established.'}


def capacity_frontier_controls():
    """Capacity/admission audit and a classical unrestricted subset baseline."""
    from math import comb

    class FullSubsets:
        # Independent current values live at the full (length, support)
        # cells themselves. A separate reference is used only by the oracle.
        def __init__(self, capacity, dimension, max_weight, noise=1, grid_bits=8, query_support=None):
            self.n, self.d, self.w = capacity, dimension, max_weight
            self.eta, self.g = noise, grid_bits
            self.bias = noise + grid_bits + (8 * (capacity + 1) - 1).bit_length()
            self.tmax = min(dimension, max_weight + noise if query_support is None else query_support)
            self.budget = 2 * self.bias * max_weight + noise
            self.qmax = sum(sum(comb(self.tmax, j) for j in range(min(t, self.tmax) + 1))
                            for t in range(max_weight + 1))
            self.guard = (2 * grid_bits + 2 + self.budget + self.bias * max_weight
                          + (capacity + 1).bit_length() + self.qmax.bit_length() + 4)
            self.keys, self.handles, self.cells = {}, {}, {}
            self.stats = dict(reads=0, lookups=0, overwrites=0, deletes=0,
                              admitted=0, rejected=0, negative_coefficients=0,
                              max_cells=0, max_intermediate_bits=0, live_zero_sums=0)

        def support(self, raw):
            return tuple(i for i in range(self.d) if raw >> i & 1)

        def full(self, raw):
            support = self.support(raw)
            return len(support), support

        def value(self, slot):
            count, value = self.cells[self.full(self.keys[slot])]
            assert count == 1
            return tuple(value)

        def snapshot(self):
            return (tuple(sorted(self.keys.items())), tuple(sorted(self.handles.items())),
                    tuple(sorted((a, count, tuple(total)) for a, (count, total) in self.cells.items())))

        def remove(self, slot):
            if slot not in self.keys:
                return
            old = self.value(slot)
            for a in self.handles.pop(slot):
                count, total = self.cells[a]
                count -= 1
                total = [x - y for x, y in zip(total, old)]
                if count:
                    self.cells[a] = count, total
                else:
                    assert total == [0, 0]
                    del self.cells[a]
            del self.keys[slot]
            self.stats['deletes'] += 1

        def replace(self, slot, raw, value):
            assert 0 <= slot < self.n and all(abs(x) <= 1 << self.g for x in value)
            valid = 0 <= raw < 1 << self.d and raw.bit_count() <= self.w
            if valid:
                valid = self.full(raw) not in self.cells or self.keys.get(slot) == raw
            if not valid:
                self.stats['rejected'] += 1
                return False
            self.remove(slot)
            support = self.support(raw)
            handles = tuple((len(support), a) for j in range(len(support) + 1)
                            for a in combinations(support, j))
            self.keys[slot], self.handles[slot] = raw, handles
            for a in handles:
                count, total = self.cells.get(a, (0, [0, 0]))
                self.cells[a] = count + 1, [x + y for x, y in zip(total, value)]
            assert self.value(slot) == tuple(value)
            self.stats['admitted'] += 1
            self.stats['max_cells'] = max(self.stats['max_cells'], len(self.cells))
            return True

        def overwrite(self, slot, value):
            assert all(abs(x) <= 1 << self.g for x in value)
            delta = [x - y for x, y in zip(value, self.value(slot))]
            for a in self.handles[slot]:
                count, total = self.cells[a]
                self.cells[a] = count, [x + y for x, y in zip(total, delta)]
            assert self.value(slot) == tuple(value)
            self.stats['overwrites'] += 1

        def check(self, reference):
            assert self.keys == {slot: raw for slot, (raw, _) in reference.items()}
            expected = {}
            for slot, (raw, value) in reference.items():
                assert self.value(slot) == value
                support = self.support(raw)
                for j in range(len(support) + 1):
                    for subset in combinations(support, j):
                        a = len(support), subset
                        count, total = expected.get(a, (0, [0, 0]))
                        expected[a] = count + 1, [x + y for x, y in zip(total, value)]
            assert self.cells == expected
            assert len(self.cells) <= len(reference) * 2 ** self.w
            self.stats['live_zero_sums'] += sum(count > 0 and total == [0, 0]
                                               for count, total in self.cells.values())

        def read(self, query):
            if (len(query) != self.d or any(type(x) is not int for x in query)
                    or sum(x != 0 for x in query) > self.tmax
                    or sum(abs(x) for x in query) > self.budget):
                raise ValueError('raw query exceeds the declared support/amplitude budgets')
            self.stats['reads'] += 1
            support = tuple(a for a, x in enumerate(query) if x)
            negative = sum(-x for x in query if x < 0)
            numerator, denominator = [0, 0], 0
            for t in range(self.w + 1):
                for j in range(min(t, len(support)) + 1):
                    for subset in combinations(support, j):
                        coefficient = 1 << (negative + self.bias * (self.w - t)
                                            - sum(-query[a] for a in subset if query[a] < 0))
                        for a in subset:
                            coefficient *= ((1 << abs(query[a])) - 1) * (-1 if query[a] < 0 else 1)
                        count, total = self.cells.get((t, subset), (0, [0, 0]))
                        denominator += coefficient * count
                        numerator = [x + coefficient * y for x, y in zip(numerator, total)]
                        self.stats['lookups'] += 1
                        self.stats['negative_coefficients'] += coefficient < 0
                        bits = max(abs(x).bit_length() + 1 for x in [coefficient, denominator, *numerator])
                        self.stats['max_intermediate_bits'] = max(self.stats['max_intermediate_bits'], bits)
                        assert bits <= self.guard
            assert denominator > 0 if self.keys else denominator == 0 and numerator == [0, 0]
            self.last_totals = tuple(numerator), denominator, negative + self.bias * self.w
            return tuple(F(x, denominator) for x in numerator) if denominator else None

    def replace(index, reference, slot, raw, value):
        expected = (0 <= raw < 1 << index.d and raw.bit_count() <= index.w
                    and all(i == slot or key != raw for i, (key, _) in reference.items()))
        before = index.snapshot()
        assert index.replace(slot, raw, value) == expected
        if expected:
            reference[slot] = raw, tuple(value)
        else:
            assert index.snapshot() == before
        index.check(reference)
        return expected

    def verify(index, reference, query):
        result = index.read(query)
        scale = sum(-x for x in query if x < 0) + index.bias * index.w
        numerator, denominator = [0, 0], 0
        for raw, value in reference.values():
            # Direct biased dot-product weights, independent of subset cells.
            score = sum(query[a] for a in range(index.d) if raw >> a & 1) - index.bias * raw.bit_count()
            assert 0 <= scale + score <= index.budget + index.bias * index.w
            weight = 1 << (scale + score)
            denominator += weight
            numerator = [x + weight * y for x, y in zip(numerator, value)]
        assert index.last_totals == (tuple(numerator), denominator, scale)
        expected = tuple(F(x, denominator) for x in numerator) if denominator else None
        assert result == expected
        return result

    cases, noise_reads, raw_rejections = [], 0, 0
    for dimension, weight in ((5, 2), (6, 3), (7, 3), (8, 4)):
        keys = [sum(1 << a for a in subset) for t in range(weight + 1)
                for subset in combinations(range(dimension), t)]
        index, reference = FullSubsets(len(keys), dimension, weight, query_support=2 * weight), {}
        for slot, raw in enumerate(keys):
            assert replace(index, reference, slot, raw, (slot % 513 - 256, (3 * slot) % 513 - 256))
        assert len(index.cells) == sum(comb(dimension, j) * (weight - j + 1) for j in range(weight + 1))
        queries = [q for q in product((-1, 0, 1), repeat=dimension)
                   if sum(x != 0 for x in q) <= index.tmax]
        for a in range(dimension):
            for b in range(a + 1, dimension):
                q = [0] * dimension
                q[a], q[b] = 3, -2
                queries.append(tuple(q))
        for q in queries:
            verify(index, reference, q)
        for slot in range(3):
            for phase in range(2):
                if phase:
                    value = tuple(-x for x in reference[slot][1])
                    index.overwrite(slot, value)
                    reference[slot] = reference[slot][0], value
                    index.check(reference)
                center = tuple(2 * index.bias * int(reference[slot][0] >> a & 1) for a in range(dimension))
                for a in range(-1, dimension):
                    for delta in ((0,) if a < 0 else (-1, 1)):
                        q = list(center)
                        if a >= 0:
                            q[a] += delta
                        result = verify(index, reference, tuple(q))
                        assert tuple(round(x) for x in result) == reference[slot][1]
                        noise_reads += 1
        for q in ((), (index.budget + 1,) + (0,) * (dimension - 1), (0.5,) + (0,) * (dimension - 1)):
            before = index.snapshot()
            try:
                index.read(q)
                assert False
            except ValueError:
                pass
            assert index.snapshot() == before
            raw_rejections += 1
        assert not replace(index, reference, 0, reference[1][0], (17, -19))
        assert not replace(index, reference, 0, 1 << dimension, (17, -19))
        assert not replace(index, reference, 0, (1 << (weight + 1)) - 1, (17, -19))
        first, second = reference[0], reference[1]
        index.remove(0); reference.pop(0)
        index.remove(1); reference.pop(1)
        assert replace(index, reference, 0, second[0], (23, -31))
        assert replace(index, reference, 1, first[0], (-37, 41))
        for q in queries[:12] + queries[-12:]:
            verify(index, reference, q)
        for slot in list(reference):
            index.remove(slot); reference.pop(slot)
        index.check(reference)
        verify(index, reference, (0,) * dimension)
        assert replace(index, reference, 0, 0, (0, 0))
        verify(index, reference, (0,) * dimension)
        cases.append({'dimension': dimension, 'maximum_key_weight': weight, 'capacity': len(keys),
                      'query_support_budget': index.tmax, 'query_l1_budget': index.budget,
                      'length_bias': index.bias, 'preselected_guard_bits': index.guard,
                      'controls': index.stats.copy()})

    def johnson(dimension, weight, threshold):
        bound = 1
        for j in range(threshold - 1, -1, -1):
            bound = (dimension - j) * bound // (weight - j)
        return bound

    packing_controls = 0
    for dimension in range(2, 13):
        for weight in range(1, dimension + 1):
            for threshold in range(1, weight + 1):
                upper = johnson(dimension, weight, threshold)
                assert upper <= comb(dimension, threshold) // comb(weight, threshold)
                packing_controls += 1
    finite_frontier = []
    for dimension in (128, 512, 4096, 65536):
        for weight in (4, 8, 32, 64):
            capacity = 1 << 20
            threshold = next((k for k in range(1, weight + 1)
                              if johnson(dimension, weight, k) >= capacity), None)
            if threshold is None:
                finite_frontier.append({'dimension': dimension, 'weight': weight,
                                        'capacity': capacity, 'necessary_packing_feasible': False})
                continue
            finite_frontier.append({'dimension': dimension, 'weight': weight, 'capacity': capacity,
                                    'necessary_packing_feasible': True,
                                    'minimum_threshold_not_ruled_out_by_johnson': threshold,
                                    'materialized_owner_entries': capacity * comb(weight, threshold),
                                    'payload_cells_changed_per_overwrite': sum(comb(weight, j) for j in range(threshold)),
                                    'binary_query_subset_visits_at_support_twice_weight': sum(comb(min(dimension, 2 * weight), j)
                                                                                          for j in range(threshold + 1)),
                                    'existence_or_retention_not_certified': True})

    jam_cases, exhaustive_insertions, exhaustive_replacements = [], 0, 0
    repair_example = None
    for prime, exponent in ((2, 2), (2, 3), (3, 2), (3, 3), (5, 2), (7, 2)):
        points = list(product(range(prime), repeat=exponent))
        point_id = {point: i for i, point in enumerate(points)}
        directions = [p for p in points if any(p) and next(x for x in p if x) == 1]
        local_lines = []
        for direction in directions:
            pivot = next(i for i, x in enumerate(direction) if x)
            for base in points:
                if base[pivot] == 0:
                    line = tuple(sorted(point_id[tuple((x + t * y) % prime for x, y in zip(base, direction))]
                                        for t in range(prime)))
                    local_lines.append(line)
        v = prime ** exponent
        dimension = (prime - 1) * v
        supports = [tuple(group * v + a for a in line)
                    for group in range(prime - 1) for line in local_lines]
        owners, loads = {}, [0] * dimension
        for slot, support in enumerate(supports):
            assert len(set(support)) == prime
            for a in support:
                loads[a] += 1
            for pair in combinations(support, 2):
                assert pair not in owners
                owners[pair] = slot
        active = len(supports)
        assert active == v * (v - 1) // prime
        assert set(loads) == {(v - 1) // (prime - 1)}
        assert len(owners) == (prime - 1) * comb(v, 2)
        assert all((a, b) in owners for group in range(prime - 1)
                   for a, b in combinations(range(group * v, (group + 1) * v), 2))
        # The elementary Turan/Caro--Wei lower bound is attained exactly.
        assert active * prime * (prime - 1) ** 2 == dimension * (dimension - prime + 1)
        proposal_count = comb(dimension, prime)
        replacement_count = sum(comb(prime, j) * comb(prime - 2, prime - j) * v ** (prime - j)
                                for j in range(2, prime + 1))
        if proposal_count <= 30000:
            actual_replacements = 0
            for proposal in combinations(range(dimension), prime):
                blockers = {owners[pair] for pair in combinations(proposal, 2) if pair in owners}
                assert blockers
                actual_replacements += blockers == {0}
                exhaustive_insertions += 1
                exhaustive_replacements += 1
            assert actual_replacements == replacement_count
        acceptance = F(replacement_count, proposal_count)
        assert acceptance * active * 2 ** (prime // 2) <= 64 * prime
        intersecting = prime * ((v - 1) // (prime - 1) - 1)
        off_numerator = intersecting * (2 ** 8 - 1) + active - 1
        off = F(off_numerator, 2 ** (8 * prime))
        # Each affine line meets q*(r-1) other lines at one point.
        raw_keys = [sum(1 << a for a in support) for support in supports]
        for key in raw_keys:
            assert sum(1 << (8 * (key & other).bit_count()) for other in raw_keys if key != other) == off_numerator
        assert 0 < off < 1
        jam_cases.append({'prime': prime, 'affine_dimension': exponent, 'key_weight': prime,
                          'raw_dimension': dimension, 'active_keys': active, 'slot_capacity': 2 * active,
                          'actual_coordinate_load': loads[0], 'declared_load_limit': 2 * loads[0],
                          'owned_pair_fraction': str(F(active * comb(prime, 2), comb(dimension, 2))),
                          'all_within_group_pairs_certified': len(owners),
                          'all_fresh_insertions_blocked': True, 'minimum_maximal_size_bound_attained': True,
                          'single_deleted_slot_uniform_replacement_acceptance': str(acceptance),
                          'near_identity_kernel_at_inverse_temperature_8_off_row_sum': str(off),
                          'proposal_family_exhausted': proposal_count <= 30000})
        if prime == 3 and exponent == 2:
            index, reference = FullSubsets(2 * active, dimension, prime), {}
            for slot, support in enumerate(supports):
                raw = sum(1 << a for a in support)
                assert replace(index, reference, slot, raw, (slot - 12, 12 - slot))
            proposal = (0, 1, v)
            assert any(pair in owners for pair in combinations(proposal, 2))
            raw = sum(1 << a for a in proposal)
            assert replace(index, reference, active, raw, (201, -197))
            for slot in (0, active):
                query = tuple(2 * index.bias * int(reference[slot][0] >> a & 1) for a in range(dimension))
                assert tuple(round(x) for x in verify(index, reference, query)) == reference[slot][1]
            repair_example = {'old_active_keys': active, 'new_active_keys': len(reference),
                              'previously_blocked_raw_support': proposal,
                              'generic_full_subset_baseline_admits_it': True,
                              'old_and_new_values_retained_by_grid_rounding': True}

    n, d, w, tmax = 1 << 20, 128, 4, 8
    low_cells = sum(sum(comb(d, j) for j in range(t)) for t in range(1, w + 1))
    resource = {'capacity': n, 'dimension': d, 'max_key_weight': w, 'query_support_budget': tmax,
                'subsets_changed_per_largest_key': 2 ** w,
                'query_cell_lookups': sum(sum(comb(tmax, j) for j in range(t + 1)) for t in range(w + 1)),
                'proper_subset_payload_cell_upper_bound': low_cells,
                'all_payload_cells_including_independent_values_upper_bound': n + low_cells,
                'input_and_state_precision_and_all_dictionary_metadata_are_additional': True,
                'no_million_slot_allocation_or_runtime_benchmark_performed': True}
    return {'generic_baseline_cases': cases, 'complete_raw_reads': sum(c['controls']['reads'] for c in cases) + 2,
            'exhaustive_integer_noise_grid_reads': noise_reads, 'atomic_invalid_queries': raw_rejections,
            'johnson_vs_simple_packing_integer_controls': packing_controls,
            'finite_frontier': finite_frontier,
            'admission_jam_families': jam_cases, 'exhaustive_blocked_insertions': exhaustive_insertions,
            'exhaustive_single_slot_replacement_proposals': exhaustive_replacements,
            'relaxed_admission_witness': repair_example, 'finite_full_subset_resource_example': resource,
            'scope': 'Exact oracle checks for the classical length-stratified all-subsets kernel with a fixed length bias and separate raw-query support/amplitude budgets. Arbitrary distinct keys of weight <=w are admitted; no low-intersection rule or independent value-cache duplicate is used. Full cells isolate each current value. Standard packing/Johnson bounds, affine-geometry minimal-maximal families and exact rekey acceptance counts audit the older admission promise. These are known kernel/design mechanisms, not a novelty or ML-advantage certificate. Fixed-grid amplitude jitter is distinct from binary-key erasure correction. Large finite resource counts are arithmetic bounds, not allocated memory or measured runtime.'}


def continuous_amplitude_controls():
    """Finite-grid sparse Taylor baseline, with exact integer certificates."""
    from math import comb, factorial

    def degree(radius, accuracy_bits):
        m = 1
        while m < radius or 9 ** radius * radius ** m * 2 ** accuracy_bits > factorial(m):
            m += 1
        return m - 1

    class SparsePolynomial:
        def __init__(self, capacity, dimension, weight, query_support, radius, grid=2, accuracy=8, value_bits=4):
            self.n, self.d, self.w, self.t = capacity, dimension, weight, query_support
            self.h, self.r, self.p, self.g = radius, grid, accuracy, value_bits
            self.l = degree(radius, accuracy)
            self.factorials = [factorial(j) for j in range(self.l + 1)]
            self.scale = self.factorials[-1] * grid ** (2 * self.l)
            self.fq = sum(comb(query_support, j) * comb(self.l, j)
                          for j in range(min(weight, query_support, self.l) + 1))
            self.guard = (value_bits + 2 * self.l * grid.bit_length() + self.factorials[-1].bit_length()
                          + self.l * max(1, radius).bit_length() + (capacity + 1).bit_length()
                          + self.fq.bit_length() + 8)
            self.keys, self.values, self.handles, self.cells = {}, {}, {}, {}
            self.stats = dict(reads=0, lookups=0, admitted=0, rejected=0, overwrites=0, deletes=0,
                              signed_features=0, live_cancelled_moments=0, max_cells=0, max_bits=0,
                              pruned_query_monomials=0)

        def guarded(self, *values):
            bits = max(abs(v).bit_length() + 1 for v in values)
            self.stats['max_bits'] = max(bits, self.stats['max_bits'])
            assert bits <= self.guard

        def features(self, raw, query=False):
            active = [(a, x) for a, x in enumerate(raw) if x]
            powers = [[x ** j for j in range(self.l + 1)] for _, x in active]
            limit = min(len(active), self.w if query else len(active), self.l)
            if query:
                self.stats['pruned_query_monomials'] += (comb(self.l + len(active), len(active))
                                                        - sum(comb(len(active), j) * comb(self.l, j)
                                                              for j in range(limit + 1)))
            constant = self.r ** self.l * (self.factorials[-1] if query else 1)
            self.guarded(constant)
            yield (), constant
            for size in range(1, limit + 1):
                for selected in combinations(range(len(active)), size):
                    for total in range(size, self.l + 1):
                        for extra in compositions(total - size, size):
                            exponents = tuple(x + 1 for x in extra)
                            alpha = tuple((active[j][0], exponent) for j, exponent in zip(selected, exponents))
                            value = self.r ** (self.l - total)
                            divisor = 1
                            for j, exponent in zip(selected, exponents):
                                value *= powers[j][exponent]
                                divisor *= self.factorials[exponent]
                            if query:
                                assert self.factorials[-1] % divisor == 0
                                value *= self.factorials[-1] // divisor
                            self.guarded(value)
                            self.stats['signed_features'] += value < 0
                            yield alpha, value

        def snapshot(self):
            return (dict(self.keys), dict(self.values), dict(self.handles), dict(self.cells))

        def remove(self, slot):
            if slot not in self.keys:
                return
            old = self.values.pop(slot)
            for alpha, feature in self.handles.pop(slot):
                count, moment, sums = self.cells[alpha]
                count -= 1
                moment -= feature
                sums = tuple(x - feature * y for x, y in zip(sums, old))
                self.guarded(moment, *sums)
                if count:
                    self.cells[alpha] = count, moment, sums
                else:
                    assert moment == 0 and sums == (0, 0)
                    del self.cells[alpha]
            del self.keys[slot]
            self.stats['deletes'] += 1

        def replace(self, slot, raw, value):
            assert 0 <= slot < self.n and all(type(x) is int and abs(x) <= 2 ** self.g for x in value)
            valid = (len(raw) == self.d and all(type(x) is int for x in raw)
                     and sum(x != 0 for x in raw) <= self.w and sum(x * x for x in raw) <= self.r ** 2
                     and all(i == slot or key != raw for i, key in self.keys.items()))
            if not valid:
                self.stats['rejected'] += 1
                return False
            handles = tuple(self.features(raw))
            assert len(handles) == comb(self.l + sum(x != 0 for x in raw), sum(x != 0 for x in raw))
            self.remove(slot)
            self.keys[slot], self.values[slot], self.handles[slot] = tuple(raw), tuple(value), handles
            for alpha, feature in handles:
                count, moment, sums = self.cells.get(alpha, (0, 0, (0, 0)))
                moment += feature
                sums = tuple(x + feature * y for x, y in zip(sums, value))
                self.guarded(moment, *sums)
                self.cells[alpha] = count + 1, moment, sums
            self.stats['admitted'] += 1
            self.stats['max_cells'] = max(self.stats['max_cells'], len(self.cells))
            return True

        def overwrite(self, slot, value):
            assert all(type(x) is int and abs(x) <= 2 ** self.g for x in value)
            delta = tuple(x - y for x, y in zip(value, self.values[slot]))
            for alpha, feature in self.handles[slot]:
                count, moment, sums = self.cells[alpha]
                sums = tuple(x + feature * y for x, y in zip(sums, delta))
                self.guarded(moment, *sums)
                self.cells[alpha] = count, moment, sums
            self.values[slot] = tuple(value)
            self.stats['overwrites'] += 1

        def check(self, reference):
            assert self.keys == {slot: raw for slot, (raw, _) in reference.items()}
            assert self.values == {slot: value for slot, (_, value) in reference.items()}
            expected = {}
            for slot, (raw, value) in reference.items():
                for alpha, feature in self.features(raw):
                    count, moment, sums = expected.get(alpha, (0, 0, (0, 0)))
                    expected[alpha] = (count + 1, moment + feature,
                                       tuple(x + feature * y for x, y in zip(sums, value)))
            assert self.cells == expected
            assert len(self.cells) <= len(self.keys) * comb(self.l + self.w, self.w)
            self.stats['live_cancelled_moments'] += sum(c > 0 and m == 0 for c, m, _ in self.cells.values())

        def read(self, raw):
            if (len(raw) != self.d or any(type(x) is not int for x in raw)
                    or sum(x != 0 for x in raw) > self.t
                    or sum(x * x for x in raw) > (self.h * self.r) ** 2):
                raise ValueError('raw query exceeds declared grid/support/norm bounds')
            numerator, denominator = [0, 0], 0
            for alpha, feature in self.features(raw, query=True):
                _, moment, sums = self.cells.get(alpha, (0, 0, (0, 0)))
                denominator += feature * moment
                numerator = [x + feature * y for x, y in zip(numerator, sums)]
                self.guarded(denominator, *numerator)
                self.stats['lookups'] += 1
            assert denominator > 0 if self.keys else denominator == 0 and numerator == [0, 0]
            self.stats['reads'] += 1
            self.last = tuple(numerator), denominator
            return tuple(F(x, denominator) for x in numerator) if denominator else None

    def direct_weight(index, raw, query):
        dot = sum(x * y for x, y in zip(raw, query))
        assert abs(dot) <= index.h * index.r ** 2
        # Powers of the scalar dot product: no monomial cells are consulted.
        weight = sum(index.factorials[-1] // index.factorials[j]
                     * dot ** j * index.r ** (2 * (index.l - j)) for j in range(index.l + 1))
        assert weight > 0
        return weight

    def verify(index, reference, query):
        result = index.read(query)
        numerator, denominator = [0, 0], 0
        for raw, value in reference.values():
            weight = direct_weight(index, raw, query)
            denominator += weight
            numerator = [x + weight * y for x, y in zip(numerator, value)]
        assert index.last == (tuple(numerator), denominator)
        assert result == (tuple(F(x, denominator) for x in numerator) if denominator else None)
        return result

    cases, retention_reads, invalid_queries = [], 0, 0
    for dimension, weight, support, radius, grid, value_bits in (
            (3, 2, 3, 1, 2, 4), (4, 2, 3, 2, 2, 4), (1, 1, 1, 9, 2, 4), (2, 2, 2, 2, 2 ** 40, 80)):
        index = SparsePolynomial(12, dimension, weight, support, radius, grid=grid, value_bits=value_bits)
        reference = {}
        if dimension >= 3:
            cores = [(0, 0, 0), (1, 0, 0), (2, 0, 0), (-1, 0, 0),
                     (0, 1, -1), (0, -1, 1), (1, 1, 0), (-1, 1, 0)]
            keys = [key + (0,) * (dimension - 3) for key in cores]
        elif dimension == 2:
            keys = [(grid - 1, 1), (grid // 2, grid // 3), (0, 0), (-grid // 2, grid // 3)]
        else:
            keys = [(2,), (-2,)]
        for slot, key in enumerate(keys):
            value = ((2 ** 78 + slot, -2 ** 77 - slot) if value_bits == 80 else (slot - 4, 4 - slot))
            assert index.replace(slot, key, value)
            reference[slot] = tuple(key), tuple(value)
        index.check(reference)
        queries = [q for q in product((-1, 0, 1), repeat=dimension)
                   if sum(x != 0 for x in q) <= support and sum(x * x for x in q) <= (grid * radius) ** 2]
        queries.extend(tuple((grid if j == a else 0) for j in range(dimension)) for a in range(dimension))
        if dimension == 2:
            queries.extend(((grid // 3, -grid // 5), (-grid, 0), (0, 2 * grid)))
        for q in queries:
            verify(index, reference, q)
        for slot in range(min(3, len(keys))):
            value = ((2 ** 79 - 1 - slot, -2 ** 78 + slot) if value_bits == 80 else (7 - slot, -7 + slot))
            index.overwrite(slot, value)
            reference[slot] = reference[slot][0], value
            index.check(reference)
            for q in queries[:5]:
                verify(index, reference, q)
        if dimension == 1:
            epsilon = F(1, 2 ** index.p)
            for slot in (0, 1):
                for phase in range(2):
                    if phase:
                        index.overwrite(slot, (-11 + slot, 11 - slot))
                        reference[slot] = reference[slot][0], (-11 + slot, 11 - slot)
                    for noise in (-2, -1, 0, 1, 2):
                        q = (8 * index.keys[slot][0] + noise,)
                        result = verify(index, reference, q)
                        # Each true exp weight is in [P/(1+eps),P/(1-eps)].
                        low, high = [F(0), F(0)], [F(0), F(0)]
                        zlow = zhigh = F(0)
                        for other, key in index.keys.items():
                            pweight = direct_weight(index, key, q)
                            wl, wu = F(pweight, 1 + epsilon), F(pweight, 1 - epsilon)
                            zlow += wl; zhigh += wu
                            for j, value in enumerate(index.values[other]):
                                low[j] += min(wl * value, wu * value)
                                high[j] += max(wl * value, wu * value)
                        for j, target in enumerate(index.values[slot]):
                            endpoints = [a / b for a in (low[j], high[j]) for b in (zlow, zhigh)]
                            assert min(endpoints) > target - F(1, 2)
                            assert max(endpoints) < target + F(1, 2)
                        assert tuple(round(x) for x in result) == index.values[slot]
                        retention_reads += 1
            # Distinct amplitudes with identical supports are valid admissions.
            assert index.replace(2, (1,), (13, -13))
            reference[2] = (1,), (13, -13)
            verify(index, reference, (3,))
            index.check(reference)
        before = index.snapshot()
        assert not index.replace(0, index.keys[1], (1, -1))
        assert not index.replace(0, (grid + 1,) + (0,) * (dimension - 1), (1, -1))
        if dimension > weight:
            assert not index.replace(0, (1,) * (weight + 1) + (0,) * (dimension - weight - 1), (1, -1))
        assert index.snapshot() == before
        invalid = [(), (grid * radius + 1,) + (0,) * (dimension - 1), (0.5,) + (0,) * (dimension - 1)]
        if dimension > support:
            invalid.append((1,) * (support + 1) + (0,) * (dimension - support - 1))
        for raw in invalid:
            before = index.snapshot()
            try:
                index.read(raw)
                assert False
            except ValueError:
                pass
            assert index.snapshot() == before
            invalid_queries += 1
        for slot in list(index.keys):
            index.remove(slot)
            reference.pop(slot)
        index.check(reference)
        assert verify(index, reference, (0,) * dimension) is None
        assert index.replace(0, (0,) * dimension, (0, 0))
        reference[0] = (0,) * dimension, (0, 0)
        assert verify(index, reference, (0,) * dimension) == (F(0), F(0))
        cases.append({'dimension': dimension, 'key_support_bound': weight, 'query_support_bound': support,
                      'query_l2_radius': radius, 'grid_denominator': index.r, 'relative_error_bits': index.p,
                      'value_unit_bits': index.g,
                      'polynomial_degree': index.l, 'preselected_guard_bits': index.guard,
                      'controls': dict(index.stats)})
    assert sum(c['controls']['signed_features'] for c in cases) > 0
    assert sum(c['controls']['live_cancelled_moments'] for c in cases) > 0
    assert sum(c['controls']['pruned_query_monomials'] for c in cases) > 0
    resource = []
    for radius in (1, 2, 4, 8, 16, 32):
        level = degree(radius, 16)
        resource.append({'score_radius': radius, 'relative_error_bits': 16, 'degree': level,
                         'features_per_weight_4_key': comb(level + 4, 4),
                         'unpruned_features_per_support_8_query': comb(level + 8, 8),
                         'query_features_after_weight_4_key_support_pruning':
                             sum(comb(8, j) * comb(level, j) for j in range(min(4, level) + 1))})
    # Finite-difference witness for any uniform scalar polynomial, not
    # just this Taylor selection. Nodes -64+j*log(10) lie in [-64,64].
    steps, ratio, radius, precision = 55, 10, 64, 16
    exp_lower = sum(F((2 * radius) ** j, factorial(j)) for j in range(2 * radius + 1))
    assert exp_lower > ratio ** steps
    assert (ratio - 1) ** steps * 2 ** precision > (ratio + 1) ** steps
    for power in range(steps):
        assert sum((-1) ** (steps - j) * comb(steps, j) * j ** power for j in range(steps + 1)) == 0
    assert sum((-1) ** (steps - j) * comb(steps, j) * ratio ** j
               for j in range(steps + 1)) == (ratio - 1) ** steps
    leading_features = sum(comb(8, j) * comb(steps - 1, j - 1) for j in range(1, 5))
    assert leading_features == 1817936 > 2 ** 20
    catalog = list(combinations(range(8), 4))
    assert all(any(set(a) <= set(key) for key in catalog)
               for j in range(1, 5) for a in combinations(range(8), j))
    obstruction = {'score_radius': radius, 'relative_error_bits': precision,
                   'minimum_possible_degree_lower_bound': steps, 'geometric_node_ratio': ratio,
                   'rational_positive_exp_series_terms': 2 * radius + 1,
                   'annihilated_monomial_degrees': steps,
                   'leading_degree_query_monomials_with_weight_4_keys_support_8_query': leading_features,
                   'small_positive_catalog_activating_all_these_support_patterns': len(catalog),
                   'scope': 'Any single scalar polynomial uniformly approximating exp on the entire real score interval [-64,64] with relative error <=2^-16 has degree >=55. Its nonzero leading coefficient yields at least 1817936 potentially active top-degree monomials even after key-support pruning. This is a feature-recipe cost obstruction, not a general query lower bound, a finite-grid-only obstruction, or a bound against local/query-adaptive approximation, alternate representations, direct scanning or sketches.'}
    existence = []
    for m in (16, 32, 64, 128, 256, 512, 1024):
        logarithm = (m + 2).bit_length() - 1
        weight = (64 * m + logarithm - 1) // logarithm
        threshold, dimension, capacity = weight // 16, weight ** 4, 2 ** m
        assert capacity ** 3 * comb(weight, threshold) ** 2 <= comb(dimension, threshold)
        # R is a power of two with 2*sqrt(w)<=R<4*sqrt(w).
        grid = 1 << (1 + ((weight - 1).bit_length() + 1) // 2)
        assert 4 * weight <= grid ** 2 < 16 * weight
        gap = F(weight - 4 * (threshold - 1), grid ** 2)
        assert gap >= F(1, 32)
        existence.append({'log2_capacity': m, 'key_weight': weight, 'raw_dimension': dimension,
                          'intersection_threshold': threshold, 'amplitude_grid_scale': grid,
                          'self_score_gap_lower_bound': str(gap),
                          'uniform_fresh_proposal_rejection_at_most_inverse_capacity_squared': True})
    # Sparse feature scores do not make finite, unmasked softmax weights sparse.
    zero_score = {'queries': [[1, 0], [1, 0]], 'keys': [[0, 1], [0, -1]],
                  'values': [2, 6], 'coordinate_overlap_events': 0,
                  'exact_outputs': ['4', '4'], 'positive_attention_entries': 4}
    assert all(sum(x * y for x, y in zip(q, k)) == 0 for q in zero_score['queries'] for k in zero_score['keys'])
    assert sum(zero_score['values']) == 2 * 4
    return {'finite_cases': cases, 'complete_polynomial_raw_reads': sum(c['controls']['reads'] for c in cases),
            'exp_interval_certified_retention_reads': retention_reads,
            'atomic_invalid_queries': invalid_queries,
            'finite_degree_resource_counts': resource, 'uniform_parameter_integer_controls': existence,
            'uniform_scalar_polynomial_obstruction': obstruction,
            'zero_score_softmax_witness': zero_score,
            'scope': 'Classical sparse monomial/Taylor feature attention with dyadic raw amplitudes, exact canonical current-value deltas, all independent values separately stored, explicit word guards and a uniform relative-exp remainder bound. Independent scalar-dot polynomial sums certify every full numerator/denominator, and rational exp-weight enclosures certify a small retention example. Seven conceptual packing/gap controls support the asymptotic nontrivial-noise family but do not allocate those states or certify a learned encoder. The feature mechanisms are prior art; no strongest-method, novelty, ML latency or publication certificate is claimed.'}


def hybrid_moment_controls():
    """Bounded fixed-point moments with complete overlap corrections."""
    from math import comb, factorial

    def ceilf(x):
        return -(-x.numerator // x.denominator)

    def clog(x):
        return (max(1, x) - 1).bit_length()

    def nearest(n, d):
        q, r = divmod(n, d)
        return q + (2 * r >= d)

    def quant(x, bits):
        return nearest(x.numerator << bits, x.denominator)

    def mul(x, y, bits):
        return nearest(x * y, 1 << bits)

    def parameters(n, d, w, k, t, amplitude, query_bound, accuracy, value_bits, input_bits):
        small = (k - 1) * amplitude * query_bound
        h = max(1, ceilf(small))
        cap = 1 << ceilf(F(3, 2) * small)
        m = 1
        while m < h or cap * h ** m * 2 ** (accuracy + 3) > factorial(m):
            m += 1
        level = m - 1
        fk = sum(comb(w, j) * comb(level, j) for j in range(1, min(k - 1, level) + 1))
        fq = sum(comb(t, j) * comb(level, j) for j in range(1, min(k - 1, t, level) + 1))
        candidate_bound = comb(t, k) if t >= k else 0
        b = clog(32 * level * cap * 2 ** accuracy)
        c = max(clog(16 * comb(level + k - 1, k - 1) * 2 ** accuracy),
                accuracy + clog(8 * max(1, candidate_bound)) - b)
        work = c + clog(64 * k * (2 * cap + 2) ** k * (level + 1))
        loads = {j: min(n, comb(d - j, k - j) // comb(w - j, k - j)) for j in range(1, k)}
        largest_load = max(loads.values(), default=n)
        big = w * amplitude * query_bound
        expcap = 1 << ceilf(F(3, 2) * big)
        bound = n + (fq * largest_load + candidate_bound * (fk + 1)) * (cap + 1)
        output_scale = accuracy + 4
        final_bound = n + fq * largest_load * (cap + 1) + candidate_bound * (expcap + 2)
        read_bits = max(b + c + value_bits + clog(bound + 1) + 2,
                        output_scale + value_bits + clog(final_bound + 2) + 2,
                        accuracy + clog(8 * max(1, n)) + value_bits + clog(n + 1) + 3)
        coefficient_bits = max(2 * work + 2 * clog(cap + 2) + 4,
                               work + clog(cap + 2) + input_bits + clog(ceilf(query_bound) + 2)
                               + amplitude.numerator.bit_length() + 4,
                               input_bits + amplitude.denominator.bit_length() + clog(level + 1) + 4)
        key_bits = 2 * b + 2 * input_bits + amplitude.numerator.bit_length() + amplitude.denominator.bit_length() + 16
        cold_scale = accuracy + clog(8 * max(1, n))
        def exp_guard(scale, radius):
            h_exp = max(1, ceilf(radius))
            exp_bound = 1 << ceilf(F(3, 2) * h_exp)
            terms = max(6 * h_exp, 2 * (scale + 8))
            frac = scale + clog(32 * (terms + 2) * exp_bound)
            score_bits = 2 * input_bits + clog(2 * ceilf(big) + 2) + 3
            return max(frac + clog(exp_bound + 2) + score_bits + 5, scale + frac + 3)
        exponential_bits = max(exp_guard(output_scale, big),
                               exp_guard(cold_scale, min(2 * big, F(cold_scale + 2))))
        return dict(degree=level, key_fraction_bits=b, coefficient_fraction_bits=c,
                    coefficient_work_bits=work, key_feature_count=fk + 1, query_feature_bound=fq + 1,
                    moment_product_fraction_bits=b + c, output_fraction_bits=output_scale,
                    candidate_bound=candidate_bound, load_bounds=loads, small_exp_bound=cap,
                    small_score_bound=str(small), full_score_bound=str(big),
                    state_payload_bits={j: b + value_bits + clog(v + 1) + 2 for j, v in loads.items()},
                    total_value_bits=value_bits + clog(n + 1) + 2, read_accumulator_bits=read_bits,
                    key_temporary_bits=key_bits, coefficient_temporary_bits=coefficient_bits,
                    exponential_temporary_bits=exponential_bits)

    interval_cache = {}

    def exp_interval(score, precision):
        cache_key = score, precision
        if cache_key in interval_cache:
            return interval_cache[cache_key]
        u, term, total, j = abs(score), F(1), F(1), 0
        while True:
            next_term = term * u / (j + 1)
            if j + 2 > u:
                tail = next_term / (1 - u / (j + 2))
                lo, hi = (total, total + tail) if score >= 0 else (1 / (total + tail), 1 / total)
                if hi - lo <= F(1, 2 ** precision):
                    interval_cache[cache_key] = lo, hi
                    return lo, hi
            j += 1
            term = next_term
            total += term

    class Hybrid:
        def __init__(self, n, d, w, k, t, amplitude=F(1, 2), query_bound=F(1), accuracy=8, value_bits=7, input_bits=4,
                     minimum_amplitude=F(0)):
            self.n, self.d, self.w, self.k, self.t = n, d, w, k, t
            self.a, self.q, self.p, self.g, self.s = amplitude, query_bound, accuracy, value_bits, input_bits
            self.minimum_amplitude = minimum_amplitude
            self.r = 1 << input_bits
            self.par = parameters(n, d, w, k, t, amplitude, query_bound, accuracy, value_bits, input_bits)
            self.l, self.b, self.c = self.par['degree'], self.par['key_fraction_bits'], self.par['coefficient_fraction_bits']
            self.work = self.par['coefficient_work_bits']
            self.keys, self.values, self.cells, self.owners = {}, {}, {}, {}
            self.loads, self.total = [0] * d, (0, 0)
            self.stats = dict(reads=0, warm_reads=0, cold_reads=0, candidate_visits=0, negative_corrections=0,
                              overwrites=0, admitted=0, rejected=0, deletes=0, live_zero_moments=0,
                              max_state_bits=0, max_read_bits=0, max_key_bits=0, max_coefficient_bits=0,
                              max_exponential_bits=0, coefficient_error_checks=0, feature_error_checks=0,
                              true_exp_enclosure_checks=0)

        def guard(self, kind, limit, *values):
            bits = max(abs(v).bit_length() + 1 for v in values)
            self.stats['max_' + kind + '_bits'] = max(self.stats['max_' + kind + '_bits'], bits)
            assert bits <= limit, (kind, bits, limit)

        def support(self, raw):
            return tuple(a for a, value in enumerate(raw) if value)

        def powers(self, raw):
            result = {}
            for a in self.support(raw):
                exact = F(raw[a], self.r) / self.a
                self.guard('key', self.par['key_temporary_bits'], exact.numerator << self.b, exact.denominator)
                root = quant(exact, self.b)
                row = [1 << self.b]
                for j in range(1, self.l + 1):
                    self.guard('key', self.par['key_temporary_bits'], row[-1] * root)
                    row.append(mul(row[-1], root, self.b))
                assert all(abs(x) <= 1 << self.b for x in row)
                result[a] = row
            return result

        def alphas(self, support):
            for size in range(1, min(self.k - 1, len(support), self.l) + 1):
                for subset in combinations(support, size):
                    for total in range(size, self.l + 1):
                        for extra in compositions(total - size, size):
                            yield tuple((a, x + 1) for a, x in zip(subset, extra))

        def feature(self, alpha, powers):
            result = 1 << self.b
            for a, exponent in alpha:
                self.guard('key', self.par['key_temporary_bits'], result * powers[a][exponent])
                result = mul(result, powers[a][exponent], self.b)
            return result

        def coefficients(self, raw):
            result = {}
            for a in self.support(raw):
                t = self.a * F(raw[a], self.r)
                self.guard('coefficient', self.par['coefficient_temporary_bits'], t.numerator, t.denominator)
                row = [1 << self.work]
                for j in range(1, self.l + 1):
                    self.guard('coefficient', self.par['coefficient_temporary_bits'], row[-1] * t.numerator)
                    row.append(nearest(row[-1] * t.numerator, t.denominator * j))
                result[a] = row
            return result

        def coefficient(self, alpha, powers):
            result = 1 << self.work
            for a, exponent in alpha:
                self.guard('coefficient', self.par['coefficient_temporary_bits'], result * powers[a][exponent])
                result = mul(result, powers[a][exponent], self.work)
            return nearest(result, 1 << (self.work - self.c))

        def exponential(self, score, scale, radius):
            # All operations are fixed-point except parsing the exact raw score.
            if score == 0:
                return 1 << scale
            if score < -(scale + 2):
                return 0
            h = max(1, ceilf(radius))
            assert abs(score) <= h
            bound = 1 << ceilf(F(3, 2) * h)
            terms = max(6 * h, 2 * (scale + 8))
            frac = scale + clog(32 * (terms + 2) * bound)
            term, total, u = 1 << frac, 1 << frac, abs(score)
            for j in range(1, terms + 1):
                self.guard('exponential', self.par['exponential_temporary_bits'], term * u.numerator)
                term = nearest(term * u.numerator, u.denominator * j)
                total += term
                self.guard('exponential', self.par['exponential_temporary_bits'], total)
            if score > 0:
                value = nearest(total, 1 << (frac - scale))
            else:
                self.guard('exponential', self.par['exponential_temporary_bits'], 1 << (scale + frac))
                value = nearest(1 << (scale + frac), total)
            lo, hi = exp_interval(score, scale + 6)  # Verification-only enclosure.
            assert F(value, 1 << scale) - F(1, 1 << scale) <= lo
            assert hi <= F(value, 1 << scale) + F(1, 1 << scale)
            return value

        def snapshot(self):
            return dict(self.keys), dict(self.values), dict(self.cells), dict(self.owners), tuple(self.loads), self.total

        def remove(self, slot):
            if slot not in self.keys:
                return
            raw, old = self.keys[slot], self.values[slot]
            support, powers = self.support(raw), self.powers(raw)
            for alpha in self.alphas(support):
                feature = self.feature(alpha, powers)
                count, moment, sums = self.cells[alpha]
                count -= 1; moment -= feature
                sums = tuple(x - feature * y for x, y in zip(sums, old))
                self.guard('state', self.par['state_payload_bits'][len(alpha)], moment, *sums)
                if count:
                    self.cells[alpha] = count, moment, sums
                else:
                    assert moment == 0 and sums == (0, 0)
                    del self.cells[alpha]
            for subset in combinations(support, self.k):
                assert self.owners.pop(subset) == slot
            for a in support:
                self.loads[a] -= 1
            self.total = tuple(x - y for x, y in zip(self.total, old))
            self.guard('state', self.par['total_value_bits'], *self.total)
            del self.keys[slot]; del self.values[slot]
            self.stats['deletes'] += 1

        def replace(self, slot, raw, value):
            assert 0 <= slot < self.n and len(value) == 2 and all(type(x) is int and abs(x) <= 2 ** self.g for x in value)
            valid = (len(raw) == self.d and all(type(x) is int for x in raw)
                     and sum(x != 0 for x in raw) == self.w
                     and all(abs(x) <= self.a * self.r and (not x or abs(x) >= self.minimum_amplitude * self.r)
                             for x in raw))
            support = self.support(raw) if valid else ()
            if valid:
                valid = all(subset not in self.owners or self.owners[subset] == slot
                            for subset in combinations(support, self.k))
            if not valid:
                self.stats['rejected'] += 1
                return False
            self.remove(slot)
            powers = self.powers(raw)
            for alpha in self.alphas(support):
                feature = self.feature(alpha, powers)
                count, moment, sums = self.cells.get(alpha, (0, 0, (0, 0)))
                sums = tuple(x + feature * y for x, y in zip(sums, value))
                self.guard('state', self.par['state_payload_bits'][len(alpha)], moment + feature, *sums)
                self.cells[alpha] = count + 1, moment + feature, sums
            for subset in combinations(support, self.k):
                self.owners[subset] = slot
            for a in support:
                self.loads[a] += 1
            self.keys[slot], self.values[slot] = tuple(raw), tuple(value)
            self.total = tuple(x + y for x, y in zip(self.total, value))
            self.guard('state', self.par['total_value_bits'], *self.total)
            self.stats['admitted'] += 1
            return True

        def overwrite(self, slot, value):
            assert len(value) == 2 and all(type(x) is int and abs(x) <= 2 ** self.g for x in value)
            delta = tuple(x - y for x, y in zip(value, self.values[slot]))
            powers = self.powers(self.keys[slot])
            for alpha in self.alphas(self.support(self.keys[slot])):
                feature = self.feature(alpha, powers)
                count, moment, sums = self.cells[alpha]
                sums = tuple(x + feature * y for x, y in zip(sums, delta))
                self.guard('state', self.par['state_payload_bits'][len(alpha)], moment, *sums)
                self.cells[alpha] = count, moment, sums
            self.total = tuple(x + y for x, y in zip(self.total, delta))
            self.guard('state', self.par['total_value_bits'], *self.total)
            self.values[slot] = tuple(value)
            self.stats['overwrites'] += 1

        def check(self, reference):
            assert self.keys == {i: raw for i, (raw, _) in reference.items()}
            assert self.values == {i: value for i, (_, value) in reference.items()}
            expected, owners, loads, total = {}, {}, [0] * self.d, (0, 0)
            for slot, (raw, value) in reference.items():
                support, powers = self.support(raw), self.powers(raw)
                total = tuple(x + y for x, y in zip(total, value))
                for a in support:
                    loads[a] += 1
                for subset in combinations(support, self.k):
                    assert subset not in owners
                    owners[subset] = slot
                for alpha in self.alphas(support):
                    feature = self.feature(alpha, powers)
                    exact = F(1)
                    for a, exponent in alpha:
                        exact *= (F(raw[a], self.r) / self.a) ** exponent
                    assert abs(F(feature, 1 << self.b) - exact) <= F(2 * self.l, 1 << self.b)
                    self.stats['feature_error_checks'] += 1
                    count, moment, sums = expected.get(alpha, (0, 0, (0, 0)))
                    expected[alpha] = count + 1, moment + feature, tuple(x + feature * y for x, y in zip(sums, value))
            assert self.cells == expected and self.owners == owners and self.loads == loads and self.total == total
            for alpha, (count, moment, sums) in self.cells.items():
                assert count <= self.par['load_bounds'][len(alpha)]
                self.guard('state', self.par['state_payload_bits'][len(alpha)], moment, *sums)
                self.stats['live_zero_moments'] += moment == 0
            self.guard('state', self.par['total_value_bits'], *self.total)

        def score(self, raw, query):
            return F(sum(x * y for x, y in zip(raw, query)), self.r ** 2)

        def read(self, raw):
            if (len(raw) != self.d or any(type(x) is not int for x in raw)
                    or sum(x != 0 for x in raw) > self.t or any(abs(x) > self.q * self.r for x in raw)):
                raise ValueError('query violates declared raw grid/support/amplitude bounds')
            self.stats['reads'] += 1
            if not self.keys:
                self.last = dict(mode='empty', numerator=(0, 0), denominator=0, scale=0, candidates=(), shift=F(0))
                return None
            support = self.support(raw)
            hits = sum(self.loads[a] for a in support)
            if len(self.keys) < 2 * hits:
                scores = {i: self.score(key, raw) for i, key in self.keys.items()}
                shift = max(scores.values())
                scale = self.p + clog(8 * len(self.keys))
                radius = min(2 * self.w * self.a * self.q, F(scale + 2))
                numerator, denominator = [0, 0], 0
                for i, score in scores.items():
                    weight = self.exponential(score - shift, scale, radius)
                    denominator += weight
                    numerator = [x + weight * y for x, y in zip(numerator, self.values[i])]
                    self.guard('read', self.par['read_accumulator_bits'], denominator, *numerator)
                self.last = dict(mode='cold', numerator=tuple(numerator), denominator=denominator,
                                 scale=scale, candidates=(), shift=shift)
                self.stats['cold_reads'] += 1
            else:
                scale, powers = self.b + self.c, self.coefficients(raw)
                numerator, denominator = [x << scale for x in self.total], len(self.keys) << scale
                for alpha in self.alphas(support):
                    coefficient = self.coefficient(alpha, powers)
                    exact = F(1)
                    for a, exponent in alpha:
                        exact *= (self.a * F(raw[a], self.r)) ** exponent / factorial(exponent)
                    assert abs(F(coefficient, 1 << self.c) - exact) <= F(1, 1 << self.c)
                    self.stats['coefficient_error_checks'] += 1
                    _, moment, sums = self.cells.get(alpha, (0, 0, (0, 0)))
                    denominator += coefficient * moment
                    numerator = [x + coefficient * y for x, y in zip(numerator, sums)]
                    self.guard('read', self.par['read_accumulator_bits'], denominator, *numerator)
                candidates = sorted({self.owners[a] for a in combinations(support, self.k) if a in self.owners})
                outscale = self.par['output_fraction_bits']
                assert outscale <= scale
                divisor = 1 << (scale - outscale)
                candidate_numerator, candidate_denominator = [0, 0], 0
                for slot in candidates:
                    key = self.keys[slot]
                    overlap = tuple(a for a in support if key[a])
                    keypowers = self.powers(key)
                    base = 1 << scale
                    for alpha in self.alphas(overlap):
                        base += self.coefficient(alpha, powers) * self.feature(alpha, keypowers)
                    weight = self.exponential(self.score(key, raw), outscale, self.w * self.a * self.q)
                    self.stats['negative_corrections'] += weight < nearest(base, divisor)
                    denominator -= base
                    numerator = [x - base * y for x, y in zip(numerator, self.values[slot])]
                    candidate_denominator += weight
                    candidate_numerator = [x + weight * y for x, y in zip(candidate_numerator, self.values[slot])]
                    self.guard('read', self.par['read_accumulator_bits'], denominator, *numerator)
                    self.guard('read', self.par['read_accumulator_bits'], candidate_denominator, *candidate_numerator)
                denominator = nearest(denominator, divisor) + candidate_denominator
                numerator = [nearest(x, divisor) + y for x, y in zip(numerator, candidate_numerator)]
                self.guard('read', self.par['read_accumulator_bits'], denominator, *numerator)
                self.last = dict(mode='warm', numerator=tuple(numerator), denominator=denominator,
                                 scale=outscale, candidates=tuple(candidates), shift=F(0))
                self.stats['warm_reads'] += 1
                self.stats['candidate_visits'] += len(candidates)
            assert denominator > 0
            return tuple(F(x, denominator) for x in numerator)

    def verify(index, reference, query, true_enclosure=True):
        result = index.read(query)
        if not reference:
            assert result is None
            return result
        mode, scale = index.last['mode'], index.last['scale']
        numerator, denominator = [0, 0], 0
        low_numerator, low_denominator = [0, 0], 0
        qsupport = tuple(a for a, x in enumerate(query) if x)
        powers = index.coefficients(query) if mode == 'warm' else None
        candidates = []
        lower, upper, zlo, zhi = [F(0), F(0)], [F(0), F(0)], F(0), F(0)
        for slot, (key, value) in reference.items():
            score = index.score(key, query)
            overlap = tuple(a for a in qsupport if key[a])
            if mode == 'cold':
                weight = index.exponential(score - index.last['shift'], scale,
                                           min(2 * index.w * index.a * index.q, F(scale + 2)))
            elif len(overlap) >= index.k:
                candidates.append(slot)
                weight = index.exponential(score, scale, index.w * index.a * index.q)
            else:
                weight = 1 << (index.b + index.c)
                kpowers = index.powers(key)
                # Per-key expansion independent of aggregate cells and owners.
                for total in range(1, index.l + 1):
                    if overlap:
                        for exponents in compositions(total, len(overlap)):
                            alpha = tuple((a, e) for a, e in zip(overlap, exponents) if e)
                            weight += index.coefficient(alpha, powers) * index.feature(alpha, kpowers)
            if mode == 'warm' and len(overlap) < index.k:
                low_denominator += weight
                low_numerator = [x + weight * y for x, y in zip(low_numerator, value)]
            else:
                denominator += weight
                numerator = [x + weight * y for x, y in zip(numerator, value)]
            if true_enclosure:
                lo, hi = exp_interval(score - index.last['shift'], index.p + clog(len(reference) + 1) + 12)
                zlo += lo; zhi += hi
                for j, v in enumerate(value):
                    lower[j] += min(lo * v, hi * v)
                    upper[j] += max(lo * v, hi * v)
        if mode == 'warm':
            divisor = 1 << (index.b + index.c - scale)
            denominator += nearest(low_denominator, divisor)
            numerator = [x + nearest(y, divisor) for x, y in zip(numerator, low_numerator)]
        assert index.last['numerator'] == tuple(numerator) and index.last['denominator'] == denominator
        assert mode != 'warm' or index.last['candidates'] == tuple(sorted(candidates))
        if true_enclosure:
            tolerance = F(2 ** (index.g + 1), 2 ** index.p - 1)
            for j, answer in enumerate(result):
                endpoints = [x / y for x in (lower[j], upper[j]) for y in (zlo, zhi)]
                assert answer - tolerance <= min(endpoints) and max(endpoints) <= answer + tolerance
            index.stats['true_exp_enclosure_checks'] += 1
        return result

    cases, outside_admissions, noncandidate_influence, invalid_queries = [], 0, 0, 0
    for threshold, input_bits in ((2, 4), (3, 4), (2, 80)):
        groups = 8
        if threshold == 2:
            local = [(0, 1, 3), (0, 2, 6), (0, 4, 5), (1, 2, 4), (1, 5, 6), (2, 3, 5), (3, 4, 6)]
            base_d, weight = 7, 3
        else:
            local = [tuple(x for x in range(8) if (x & normal).bit_count() % 2 == parity)
                     for normal in range(1, 8) for parity in range(2)]
            base_d, weight = 8, 4
        supports = [tuple(g * base_d + a for a in key) for g in range(groups) for key in local]
        index = Hybrid(2 * len(supports), groups * base_d, weight, threshold, 8, input_bits=input_bits)
        reference = {}
        for slot, support in enumerate(supports):
            raw = [0] * index.d
            for j, a in enumerate(support):
                raw[a] = (index.r // 2 - (slot + j) % 3) * (-1 if (slot + j) % 4 == 0 else 1)
            value = slot % 31 - 15, 15 - slot % 31
            assert index.replace(slot, tuple(raw), value)
            reference[slot] = tuple(raw), value
            if slot in (0, 3):
                verify(index, reference, tuple(index.r if a in support else 0 for a in range(index.d)))
        index.check(reference)
        queries = [tuple(x * index.r for x in small) + (0,) * (index.d - 4)
                   for small in product((-1, 0, 1), repeat=4)]
        for slot in range(4):
            q = [0] * index.d
            for j, a in enumerate(supports[slot]):
                q[a] = index.r * (-1 if j % 2 else 1)
            queries.append(tuple(q))
        fractional = [0] * index.d
        fractional[0], fractional[1] = index.r // 3 + 1, -index.r // 5 - 1
        queries.append(tuple(fractional))
        for q in queries:
            verify(index, reference, q)
        for slot in range(3):
            value = 31 - slot, -31 + slot
            index.overwrite(slot, value); reference[slot] = reference[slot][0], value
            index.check(reference)
            verify(index, reference, queries[-1])
        q = queries[-1]
        outside = len(local) * (groups - 1)
        assert not any(reference[outside][0][a] and q[a] for a in range(index.d))
        before = verify(index, reference, q)
        old = reference[outside][1]
        value = old[0] + 1, old[1] - 1
        index.overwrite(outside, value); reference[outside] = reference[outside][0], value
        after = verify(index, reference, q)
        assert before != after
        noncandidate_influence += 1
        raw = [0] * index.d
        for g in range(weight):
            raw[g * base_d] = index.r // 2
        slot = len(supports)
        assert index.replace(slot, tuple(raw), (37, -37))
        reference[slot] = tuple(raw), (37, -37)
        outside_admissions += 1
        index.check(reference)
        verify(index, reference, tuple(index.r if x else 0 for x in raw))
        changed = list(reference[0][0])
        changed[index.support(changed)[0]] *= -1
        assert index.replace(0, tuple(changed), reference[0][1])
        reference[0] = tuple(changed), reference[0][1]
        index.check(reference)
        changed = [0] * index.d
        for g in range(weight):
            changed[g * base_d + 1] = index.r // 2
        assert index.replace(0, tuple(changed), reference[0][1])
        reference[0] = tuple(changed), reference[0][1]
        outside_admissions += 1
        index.check(reference)
        verify(index, reference, queries[-1])
        before = index.snapshot()
        assert not index.replace(0, reference[1][0], (1, -1))
        assert not index.replace(0, (index.r,) * weight + (0,) * (index.d - weight), (1, -1))
        assert index.snapshot() == before
        for invalid in ((), (index.r + 1,) + (0,) * (index.d - 1), (0.5,) + (0,) * (index.d - 1),
                        (1,) * 9 + (0,) * (index.d - 9)):
            before = index.snapshot()
            try:
                index.read(invalid)
                assert False
            except ValueError:
                pass
            assert index.snapshot() == before
            invalid_queries += 1
        for slot in list(reference):
            index.remove(slot); reference.pop(slot)
        index.check(reference)
        assert verify(index, reference, (0,) * index.d) is None
        assert index.replace(0, tuple(raw), (0, 0))
        reference[0] = tuple(raw), (0, 0)
        verify(index, reference, tuple(index.r if x else 0 for x in raw))
        index.check(reference)
        cases.append({'threshold': threshold, 'dimension': index.d, 'weight': weight,
                      'initial_keys': len(supports), 'input_fractional_bits': input_bits,
                      'parameters': index.par, 'controls': dict(index.stats)})
    assert sum(c['controls']['negative_corrections'] for c in cases) > 0
    assert sum(c['controls']['cold_reads'] for c in cases) > 0
    assert sum(c['controls']['warm_reads'] for c in cases) > 0
    assert sum(c['controls']['live_zero_moments'] for c in cases) > 0
    # Exercise the advertised scalar widths on a small, actually installed
    # subfamily at the large capacity's fixed parameters, not a large model.
    declared_n, declared_d = 1031 ** 2, 8 * 1031
    retention = Hybrid(declared_n, declared_d, 8, 2, 16, amplitude=F(3, 8), query_bound=F(12),
                       accuracy=14, value_bits=8, input_bits=16, minimum_amplitude=F(11, 32))
    reference = {}
    for slope in range(16):
        for intercept in range(4):
            slot = 4 * slope + intercept
            raw = [0] * declared_d
            for j in range(8):
                a = j * 1031 + (slope * j + intercept) % 1031
                raw[a] = 11 * retention.r // 32 + ((17 * slot + 13 * j) % (retention.r // 32 + 1))
            value = 3 * slot - 128, 128 - 3 * slot
            assert retention.replace(slot, tuple(raw), value)
            reference[slot] = tuple(raw), value
    retention.check(reference)
    too_small = list(reference[0][0])
    too_small[retention.support(too_small)[0]] = 11 * retention.r // 32 - 1
    before = retention.snapshot()
    assert not retention.replace(0, tuple(too_small), reference[0][1])
    assert retention.snapshot() == before
    retention_reads = 0
    for slot in (0, 4):
        for phase in range(2):
            if phase:
                value = 243 - slot, -197 + slot
                retention.overwrite(slot, value)
                reference[slot] = reference[slot][0], value
                retention.check(reference)
            center = tuple(29 * x for x in reference[slot][0])
            coordinates = retention.support(reference[slot][0]) + (declared_d - 1,)
            for a in (-1,) + coordinates:
                for delta in ((0,) if a < 0 else (-retention.r // 4, retention.r // 4)):
                    query = list(center)
                    if a >= 0:
                        query[a] += delta
                    answer = verify(retention, reference, tuple(query))
                    assert retention.last['mode'] == 'warm'
                    assert tuple(round(x) for x in answer) == reference[slot][1]
                    retention_reads += 1
    assert retention.stats['max_state_bits'] <= 64 and retention.stats['max_read_bits'] <= 128
    assert retention.stats['max_key_bits'] <= 128 and retention.stats['max_coefficient_bits'] <= 128
    assert retention.stats['max_exponential_bits'] <= 256
    finite_retention = {'actually_installed_keys': len(reference), 'declared_capacity': declared_n,
                        'dimension': declared_d, 'noise_queries_before_and_after_target_overwrites': retention_reads,
                        'query_axis_noise_radius': '1/4', 'full_ball_guarantee_comes_from_uniform_proof': True,
                        'minimum_key_amplitude_enforced_on_raw_admission': '11/32',
                        'parameters': retention.par, 'controls': dict(retention.stats)}
    # Explicit finite regime; no million-value state is allocated.
    prime = 1031
    assert all(prime % j for j in range(2, 33))
    n, d, w, t, g, accuracy = prime ** 2, 8 * prime, 8, 16, 8, 14
    par = parameters(n, d, w, 2, t, F(3, 8), F(12), accuracy, g, 16)
    assert par['state_payload_bits'][1] <= 64 and par['read_accumulator_bits'] <= 128
    assert par['exponential_temporary_bits'] <= 256
    assert par['key_temporary_bits'] <= 128 and par['coefficient_temporary_bits'] <= 128
    assert n >= 2 * t * par['load_bounds'][1]
    score_gap = F(29) * (w * F(11, 32) ** 2 - F(3, 8) ** 2) - 2 * F(9, 8) * F(1, 4)
    exp_lower = sum(score_gap ** j / factorial(j) for j in range(80))
    assert exp_lower > 16 * n * 2 ** g
    assert F(29) * F(3, 8) + F(1, 4) <= 12
    rows, kmax = par['query_feature_bound'], par['candidate_bound']
    posting_rows = t * prime - (w - 1) * kmax
    hybrid_payload = rows * 8 + kmax * 2
    v = 64
    state_bytes = n * v * 2 + d * par['degree'] * (v + 1) * 8
    state_bytes += comb(d, 2) * 4 + n * w * 4 + n * 4 + (n + 7) // 8 + d * 4 + v * 4
    workspace_bytes = d * 4 + t * (par['degree'] + 1) * 8 + w * (par['degree'] + 1) * 8
    workspace_bytes += t * 2 + kmax * 4 + 4 * w * 4 + 8 * v + 3 * (v + 1) * 16
    workspace_bytes += 16 * 32 + 16 * 16 + 128 * 8 + 2048
    posting_persistent = n * (v * 2 + w * 4) + n * w * 8 + d * 4 + n * 4 + (n + 7) // 8 + v * 4
    posting_workspace = d * 4 + t * 2 + t * par['load_bounds'][1] * 4 + 2 * (v + 1) * 16
    posting_workspace += 16 * 32 + 16 * 16 + 128 * 8 + 2048
    finite = {'capacity': n, 'dimension': d, 'weight': w, 'query_support_bound': t,
              'parameters': par, 'persistent_array_bytes_at_value_width_64': state_bytes,
              'reserved_operation_input_output_and_arithmetic_workspace_bytes': workspace_bytes,
              'packed_arena_byte_upper_bound_at_value_width_64': state_bytes + workspace_bytes,
              'sparse_key_value_bytes_at_value_width_64': n * (v * 2 + w * 4),
              'posting_arena_bytes_with_doubly_linked_edges_and_workspace': posting_persistent + posting_workspace,
              'posting_arena_bytes_if_same_pair_admission_certification_is_required':
                  posting_persistent + posting_workspace + comb(d, 2) * 4,
              'hybrid_payload_bytes_per_value_coordinate_upper_bound': hybrid_payload,
              'posting_payload_bytes_per_value_coordinate_lower_bound_on_full_line_family': posting_rows * 2,
              'retention_grid_fractional_bits': g, 'retention_scale': 29,
              'retention_euclidean_noise_radius': '1/4', 'retention_score_gap': str(score_gap),
              'retention_query_support_is_separately_required': True,
              'geometry': 'All degree-one graphs at eight evaluation points over F_1031. Pairwise overlap<=1; every coordinate has load1031. The usual polynomial construction is prior art. Generic pair-owner admission need not restrict future keys to that codebook.',
              'write_payload_moment_updates_per_value_overwrite': w * par['degree'],
              'maximum_cold_scan_active_keys': 2 * t * par['load_bounds'][1] - 1,
              'scope': 'Arithmetic accounting for one packed arena with 16-bit independent values/coordinates/amplitudes, 32-bit owners/free slots, 64-bit moment cells, full raw operation input, powers, candidate IDs, separate base/candidate/output buffers, scalar temporaries, parameter storage and 2048 bytes of alignment allowance for at most 32 arrays. The posting comparator includes 32-bit next/previous edge links, coordinate heads, free/active metadata, raw input, sorting IDs, arithmetic and output workspace; pair certification is also shown separately. A caller retaining past outputs, executable/runtime memory and a concurrent/batched implementation are outside these sequential data-structure arena counts. Construction, encoder and arithmetic operation costs remain additional to byte capacity. No large state or arena implementation, hardware benchmark, encoder or neural training was run. Posting figures count logical payload loads for one faithful method, not physical bandwidth, a speedup, an equal-memory result or a lower bound for all published hybrids or sketches.'}
    return {'cases': cases, 'total_complete_reads': sum(c['controls']['reads'] for c in cases) + retention.stats['reads'],
            'atomic_invalid_queries': invalid_queries,
            'finite_width_installed_retention': finite_retention,
            'generic_admissions_outside_initial_family': outside_admissions,
            'noncandidate_payload_influence_witnesses': noncandidate_influence,
            'finite_packed_layout_accounting': finite,
            'scope': 'Quantized support-truncated polynomial moments plus complete high-overlap residuals in the existing Scatterbrain sparse-plus-low-rank template. The warm path certifies a large exact zero-score background; the cold path scans the explicitly charged smaller/current state. Independent per-key quantized formulas and rational true-exp enclosures verify full outputs. Cached integer widths are independent of input-precision times degree; scalar exponential evaluation is finite fixed-point and separately guarded. This is a candidate cost/precision theorem, not established novelty, strongest-method advantage, measured speed or publication readiness.'}


def exponential_comparator_controls():
    """Classical coordinate sums and a stronger full-noise recall comparator."""
    from math import comb, factorial, isqrt

    def nearest(n, d):
        q, r = divmod(n, d)
        return q + (2 * r >= d)

    def clog(n):
        return (max(1, n) - 1).bit_length()

    @lru_cache(maxsize=None)
    def exp_interval(score, precision=64):
        u, term, total, j = abs(score), F(1), F(1), 0
        while True:
            next_term = term * u / (j + 1)
            if j + 2 > u:
                tail = next_term / (1 - u / (j + 2))
                lo, hi = (total, total + tail) if score >= 0 else (1 / (total + tail), 1 / total)
                if hi - lo <= F(1, 2 ** precision):
                    return lo, hi
            j += 1
            term = next_term
            total += term

    coordinate_cases = []
    for input_bits in (4, 80):
        r, d, level, b, c, work = 1 << input_bits, 56, 6, 17, 15, 30
        reference, groups, moments = {}, {}, {}
        stats = dict(reads=0, snapshots=0, group_deltas=0, moment_deltas=0,
                     candidate_visits=0, cancelled_live_groups=0, noncandidate_influence=0)

        def features(u):
            root = nearest(2 * u * (1 << b), r)
            values = [1 << b]
            for j in range(level):
                values.append(nearest(values[-1] * root, 1 << b))
            return values[1:]

        def change(slot, new):
            old = reference.pop(slot, None)
            for record, sign in ((old, -1), (new, 1)):
                if record is None:
                    continue
                raw, value = record
                for a, u in enumerate(raw):
                    if not u:
                        continue
                    key = a, u
                    previous = groups.get(key, (0, 0, 0))
                    current = tuple(x + sign * y for x, y in zip(previous, (1,) + value))
                    if current[0]:
                        groups[key] = current
                    else:
                        assert current == (0, 0, 0)
                        groups.pop(key)
                    stats['group_deltas'] += 1
                    for j, phi in enumerate(features(u), 1):
                        key = a, j
                        previous = moments.get(key, (0, 0, 0, 0))
                        current = tuple(x + sign * y for x, y in zip(previous, (1, phi) + tuple(phi * v for v in value)))
                        if current[0]:
                            moments[key] = current
                        else:
                            assert current == (0, 0, 0, 0)
                            moments.pop(key)
                        stats['moment_deltas'] += 1
            if new is not None:
                reference[slot] = new

        def check():
            expected_groups, expected_moments = {}, {}
            for raw, value in reference.values():
                for a, u in enumerate(raw):
                    if not u:
                        continue
                    key = a, u
                    previous = expected_groups.get(key, (0, 0, 0))
                    expected_groups[key] = tuple(x + y for x, y in zip(previous, (1,) + value))
            for (a, u), (count, v0, v1) in expected_groups.items():
                for j, phi in enumerate(features(u), 1):
                    key = a, j
                    previous = expected_moments.get(key, (0, 0, 0, 0))
                    expected_moments[key] = tuple(x + y for x, y in zip(previous, (count, count * phi, phi * v0, phi * v1)))
            assert groups == expected_groups and moments == expected_moments
            stats['cancelled_live_groups'] += sum(v[0] > 0 and v[1:] == (0, 0) for v in groups.values())
            stats['snapshots'] += 1

        def read(query):
            coefficients = {}
            for a, q in enumerate(query):
                if not q:
                    continue
                row = [1 << work]
                for j in range(1, level + 1):
                    row.append(nearest(row[-1] * q, 2 * r * j))
                    coefficients[a, j] = nearest(row[-1], 1 << (work - c))
            total = (len(reference),) + tuple(sum(v[j] for _, v in reference.values()) for j in range(2))
            cached = [x << (b + c) for x in total]
            grouped = list(cached)
            direct_low, symbolic = [0, 0, 0], {}
            for (a, j), coefficient in coefficients.items():
                cell = moments.get((a, j), (0, 0, 0, 0))
                cached = [x + coefficient * y for x, y in zip(cached, cell[1:])]
            for (a, u), (count, v0, v1) in groups.items():
                univariate = sum(coefficients.get((a, j), 0) * phi for j, phi in enumerate(features(u), 1))
                grouped = [x + univariate * y for x, y in zip(grouped, (count, v0, v1))]
            assert cached == grouped
            for raw, value in reference.values():
                overlap = [a for a, u in enumerate(raw) if u and query[a]]
                base = (1 << (b + c)) + sum(coefficients[a, j] * phi for a in overlap
                                                        for j, phi in enumerate(features(raw[a]), 1))
                vector = (base,) + tuple(base * v for v in value)
                if len(overlap) >= 2:
                    cached = [x - y for x, y in zip(cached, vector)]
                    grouped = [x - y for x, y in zip(grouped, vector)]
                    score = F(sum(x * q for x, q in zip(raw, query)), r * r)
                    previous = symbolic.get(score, (0, 0, 0))
                    symbolic[score] = tuple(x + y for x, y in zip(previous, (1,) + value))
                    stats['candidate_visits'] += 1
                else:
                    direct_low = [x + y for x, y in zip(direct_low, vector)]
            assert cached == grouped == direct_low
            # With the same exponential primitive all three complete expressions
            # are identical, including final low-subtotal rounding and candidates.
            rounded = tuple(nearest(x, 1 << (b + c - 12)) for x in cached)
            stats['reads'] += 1
            return rounded, tuple(sorted((str(s), v) for s, v in symbolic.items()))

        local = [(0, 1, 3), (0, 2, 6), (0, 4, 5), (1, 2, 4), (1, 5, 6), (2, 3, 5), (3, 4, 6)]
        for slot in range(56):
            raw = [0] * d
            for a in local[slot % 7]:
                raw[(slot // 7) * 7 + a] = r // 2
            value = (-1 if slot % 2 else 1, 1 if slot % 2 else -1)
            change(slot, (tuple(raw), value))
        check()
        queries = [tuple(q * r for q in small) + (0,) * (d - 4) for small in product((-1, 0, 1), repeat=4)]
        queries.append((r // 3 + 1, -r // 5 - 1) + (0,) * (d - 2))
        for query in queries:
            read(query)
        for phase in range(4):
            raw, value = reference[phase]
            changed = list(raw)
            a = next(a for a, u in enumerate(raw) if u)
            changed[a] = -r // 2 + phase + 1
            change(phase, (tuple(changed), (11 + phase, -11 - phase)))
            check()
            for query in queries[-2:]:
                read(query)
        query = queries[-1]
        before = read(query)
        raw, value = reference[49]
        assert all(not (u and q) for u, q in zip(raw, query))
        change(49, (raw, (value[0] + 1, value[1] - 1)))
        assert read(query) != before
        stats['noncandidate_influence'] += 1
        for slot in list(reference):
            change(slot, None)
        check()
        assert not groups and not moments
        assert read((0,) * d) == ((0, 0, 0), ())
        coordinate_cases.append(dict(input_fractional_bits=input_bits, degree=level, key_fraction_bits=b,
                                     coefficient_fraction_bits=c, coefficient_work_bits=work, controls=stats))

    capacity, dimension, weight, r = 1031 ** 2, 8248, 8, 1 << 16
    load, beta, eta, minimum, amplitude = (dimension - 1) // (weight - 1), 29, F(1, 4), F(11, 32), F(3, 8)
    high_distractors = weight * (load - 1)
    self_exponent = beta * weight * minimum ** 2
    self_exp_lower = sum(self_exponent ** j / factorial(j) for j in range(80))
    # e^(2 eta * 9/8)<2 and e^(beta*A^2)<128, certified by e^2<8.
    mass_bracket = capacity - 1 + 127 * high_distractors
    output_error_bound = 4 * mass_bracket / self_exp_lower
    assert 2 * eta * F(9, 8) < F(2, 3)
    assert beta * amplitude ** 2 < F(14, 3)
    assert output_error_bound < F(1, 2 ** 16)
    assert beta * minimum - eta > 1 and eta < 1
    coefficient_mass = sum(F(9, 2) ** j / factorial(j) for j in range(1, 25))
    _, mistaken_bound = exp_interval(F(33, 8))
    assert coefficient_mass > mistaken_bound and coefficient_mass < 128

    class RecallDictionary:
        def __init__(self):
            self.keys, self.values, self.owners, self.loads = {}, {}, {}, [0] * dimension
            self.stats = dict(reads=0, fast_reads=0, full_scans=0, empty_reads=0, value_overwrites=0,
                              admitted=0, rejected=0, deletes=0, dense_noise_fast_reads=0,
                              true_exp_enclosure_checks=0, max_distance_bits=0, max_exponential_bits=0)

        def snapshot(self):
            return dict(self.keys), dict(self.values), dict(self.owners), tuple(self.loads)

        def remove(self, slot):
            if slot not in self.keys:
                return
            coordinates = tuple(a for a, _ in self.keys[slot])
            for pair in combinations(coordinates, 2):
                assert self.owners.pop(pair) == slot
            for a in coordinates:
                self.loads[a] -= 1
            del self.keys[slot]; del self.values[slot]
            self.stats['deletes'] += 1

        def replace(self, slot, raw, value):
            assert 0 <= slot < capacity and len(value) == 2 and all(type(v) is int and abs(v) <= 256 for v in value)
            valid = len(raw) == dimension and all(type(u) is int and (not u or minimum * r <= abs(u) <= amplitude * r) for u in raw)
            sparse = tuple((a, u) for a, u in enumerate(raw) if u) if valid else ()
            valid = valid and len(sparse) == weight
            coordinates = tuple(a for a, _ in sparse)
            if valid:
                valid = all(pair not in self.owners or self.owners[pair] == slot for pair in combinations(coordinates, 2))
            if not valid:
                self.stats['rejected'] += 1
                return False
            self.remove(slot)
            self.keys[slot], self.values[slot] = sparse, tuple(value)
            for pair in combinations(coordinates, 2):
                self.owners[pair] = slot
            for a in coordinates:
                self.loads[a] += 1
            self.stats['admitted'] += 1
            return True

        def overwrite(self, slot, value):
            assert len(value) == 2 and all(type(v) is int and abs(v) <= 256 for v in value)
            self.values[slot] = tuple(value)
            self.stats['value_overwrites'] += 1

        def score(self, slot, query):
            return F(sum(u * query[a] for a, u in self.keys[slot]), r * r)

        def exponential(self, score, scale):
            assert score <= 0
            if score == 0:
                return 1 << scale
            if score < -(scale + 2):
                return 0
            h = min(72, scale + 2)
            cap = 1 << ((3 * h + 1) // 2)
            terms = max(6 * h, 2 * (scale + 8))
            frac = scale + clog(32 * (terms + 2) * cap)
            u, term, total = -score, 1 << frac, 1 << frac
            for j in range(1, terms + 1):
                product = term * u.numerator
                self.stats['max_exponential_bits'] = max(self.stats['max_exponential_bits'], product.bit_length() + 1)
                term = nearest(product, u.denominator * j)
                total += term
            reciprocal_numerator = 1 << (scale + frac)
            self.stats['max_exponential_bits'] = max(self.stats['max_exponential_bits'], reciprocal_numerator.bit_length() + 1, total.bit_length() + 1)
            assert self.stats['max_exponential_bits'] <= 219
            return nearest(reciprocal_numerator, total)

        def read(self, query):
            if len(query) != dimension or any(type(q) is not int or abs(q) > 12 * r for q in query):
                raise ValueError('raw query violates declared grid, dimension or coordinate bound')
            self.stats['reads'] += 1
            self.last = 'empty', None
            if not self.keys:
                self.stats['empty_reads'] += 1
                return None
            heavy, squared = [], 0
            for a, q in enumerate(query):
                squared += q * q
                if abs(q) > r and len(heavy) <= weight:
                    heavy.append(a)
            if len(heavy) == weight:
                slot = self.owners.get(tuple(heavy[:2]))
                if slot is not None and tuple(a for a, _ in self.keys[slot]) == tuple(heavy):
                    cross = sum(u * query[a] for a, u in self.keys[slot])
                    key_squared = sum(u * u for _, u in self.keys[slot])
                    distance = squared - 2 * beta * cross + beta * beta * key_squared
                    self.stats['max_distance_bits'] = max(self.stats['max_distance_bits'], squared.bit_length() + 1,
                                                         (2 * beta * cross).bit_length() + 1, (beta * beta * key_squared).bit_length() + 1)
                    assert 0 <= distance and self.stats['max_distance_bits'] <= 64
                    if distance <= (r // 4) ** 2:
                        self.last = 'fast', slot
                        self.stats['fast_reads'] += 1
                        self.stats['dense_noise_fast_reads'] += sum(q != 0 for q in query) > 16
                        return tuple(F(v) for v in self.values[slot])
            self.last = 'scan', None
            self.stats['full_scans'] += 1
            # Two score passes; no n-length score array is required.
            shift = max(self.score(slot, query) for slot in self.keys)
            scale = 14 + clog(8 * len(self.keys))
            denominator, numerator = 0, [0, 0]
            for slot in self.keys:
                kernel = self.exponential(self.score(slot, query) - shift, scale)
                denominator += kernel
                numerator = [x + kernel * v for x, v in zip(numerator, self.values[slot])]
                assert max(abs(x).bit_length() + 1 for x in [denominator] + numerator) <= 128
            return tuple(F(x, denominator) for x in numerator)

    index, reference = RecallDictionary(), {}

    def verify(query, expect=None):
        answer = index.read(query)
        if expect is not None:
            assert index.last[0] == expect
        if not reference:
            assert answer is None
            return
        scores = {i: F(sum(x * q for x, q in zip(raw, query)), r * r) for i, (raw, _) in reference.items()}
        shift = max(scores.values())
        lo_num, hi_num, lo_den, hi_den = [F(0), F(0)], [F(0), F(0)], F(0), F(0)
        for slot, score in scores.items():
            lo, hi = exp_interval(score - shift)
            lo_den += lo; hi_den += hi
            for j, v in enumerate(reference[slot][1]):
                lo_num[j] += min(lo * v, hi * v)
                hi_num[j] += max(lo * v, hi * v)
        tolerance = F(1, 2 ** (16 if index.last[0] == 'fast' else 13))
        for j, value in enumerate(answer):
            endpoints = [x / y / 256 for x in (lo_num[j], hi_num[j]) for y in (lo_den, hi_den)]
            assert F(value, 256) - tolerance <= min(endpoints) and max(endpoints) <= F(value, 256) + tolerance
        index.stats['true_exp_enclosure_checks'] += 1

    for slope in range(16):
        for intercept in range(4):
            slot, raw = 4 * slope + intercept, [0] * dimension
            for j in range(weight):
                raw[j * 1031 + (slope * j + intercept) % 1031] = 11 * r // 32 + (17 * slot + 13 * j) % (r // 32 + 1)
            value = 3 * slot - 128, 128 - 3 * slot
            assert index.replace(slot, tuple(raw), value)
            reference[slot] = tuple(raw), value
    outside = [0] * dimension
    for a in range(80, 88):
        outside[a] = (1 if a % 2 else -1) * (3 * r // 8)
    assert index.replace(64, tuple(outside), (251, -239))
    reference[64] = tuple(outside), (251, -239)
    step = isqrt((r // 4) ** 2 // dimension)
    for slot in (0, 4, 64):
        for phase in range(2):
            if phase:
                index.overwrite(slot, (243 - slot, -181 - slot))
                reference[slot] = reference[slot][0], index.values[slot]
            center = tuple(beta * u for u in reference[slot][0])
            coordinates = tuple(a for a, u in enumerate(reference[slot][0]) if u) + (dimension - 1,)
            verify(center, 'fast')
            for a in coordinates:
                for delta in (-r // 4, r // 4):
                    query = list(center); query[a] += delta
                    verify(tuple(query), 'fast')
            for pattern in (0, 1, 2):
                noise = tuple(step * (1 if pattern == 0 or pattern == 2 and a % 2 else -1) for a in range(dimension))
                assert sum(u * u for u in noise) <= (r // 4) ** 2
                verify(tuple(x + y for x, y in zip(center, noise)), 'fast')
    # Own amplitude edits and outside-codebook keys remain within the same API.
    changed = list(reference[0][0]); changed[next(a for a, u in enumerate(changed) if u)] *= -1
    assert index.replace(0, tuple(changed), reference[0][1])
    reference[0] = tuple(changed), reference[0][1]
    verify(tuple(beta * u + step for u in changed), 'fast')
    center = tuple(beta * u for u in reference[4][0])
    query = list(center); query[next(a for a, u in enumerate(reference[4][0]) if u)] += r // 4 + 1
    verify(tuple(query), 'scan')
    verify((0,) * dimension, 'scan')
    verify((step,) * dimension, 'scan')
    verify((2 * r,) + (0,) * (dimension - 1), 'scan')
    before = index.snapshot()
    assert not index.replace(0, reference[4][0], (0, 0))
    invalid = list(reference[0][0]); invalid[next(a for a, u in enumerate(invalid) if u)] = 11 * r // 32 - 1
    assert not index.replace(0, tuple(invalid), (0, 0))
    assert before == index.snapshot()
    rejected_queries = 0
    for query in ((), (12 * r + 1,) + (0,) * (dimension - 1), (0.5,) + (0,) * (dimension - 1)):
        try:
            index.read(query)
            assert False
        except ValueError:
            rejected_queries += 1
        assert before == index.snapshot()
    for slot in list(reference):
        index.remove(slot); reference.pop(slot)
    verify((0,) * dimension, 'empty')
    assert index.replace(0, tuple(outside), (0, 0))
    reference[0] = tuple(outside), (0, 0)
    verify(tuple(beta * u for u in outside), 'fast')

    lam_lo, lam_hi = exp_interval(F(3, 2))
    assert 4 < lam_lo < lam_hi < 5
    # A four-key dyadic witness with six private padding coordinates per key.
    supports = [pair + tuple(range(4 + 6 * j, 10 + 6 * j)) for j, pair in enumerate(((0, 1), (0, 2), (1, 3), (2, 3)))]
    assert all(len(set(x) & set(y)) <= 1 for x, y in combinations(supports, 2))
    contrast = (1, -1, -1, 1)
    unary_cancellations = []
    for query_support in ((), (0,), (1,)):
        symbolic = {}
        for support, v in zip(supports, contrast):
            score = F(3, 2) * len(set(support) & set(query_support))
            symbolic[score] = symbolic.get(score, 0) + v
        assert all(v == 0 for v in symbolic.values())
        unary_cancellations.append(tuple((str(k), v) for k, v in sorted(symbolic.items())))
    interaction_lower = ((lam_lo - 1) / (lam_lo + 1)) ** 2
    top_one_error_lower = 2 * (1 - (lam_hi / (lam_hi + 1)) ** 2)
    assert interaction_lower > F(9, 25) and top_one_error_lower > F(11, 18)
    probabilities = (F(16, 25), F(4, 25), F(4, 25), F(1, 25))
    duality_checks = 0
    for approximation in ((F(1), F(0), F(0), F(0)), (F(1, 4),) * 4, (F(1, 2), F(1, 2), F(0), F(0))):
        error = max(abs(sum((p - q) * v for p, q, v in zip(probabilities, approximation, values)))
                    for values in product((-1, 1), repeat=4))
        assert error == sum(abs(p - q) for p, q in zip(probabilities, approximation))
        duality_checks += 1

    def laplacian(weights):
        m = len(weights)
        matrix = [[F(0) for _ in range(2 * m)] for _ in range(2 * m)]
        for i, row in enumerate(weights):
            for j, edge in enumerate(row):
                matrix[i][i] += edge; matrix[m + j][m + j] += edge
                matrix[i][m + j] -= edge; matrix[m + j][i] -= edge
        return matrix

    def psd(matrix):
        matrix = [list(row) for row in matrix]
        for a in range(len(matrix)):
            pivot = matrix[a][a]
            if pivot < 0:
                return False
            if not pivot:
                if any(matrix[a][j] for j in range(a + 1, len(matrix))):
                    return False
                continue
            for i in range(a + 1, len(matrix)):
                for j in range(a + 1, len(matrix)):
                    matrix[i][j] -= matrix[i][a] * matrix[a][j] / pivot
        return True

    # An actually smaller rational spectral surrogate, independently PSD-checked.
    dense_graph = [[F(1) for _ in range(4)] for _ in range(4)]
    sparse_graph = [[F(4, 3) if i != j else F(0) for j in range(4)] for i in range(4)]
    original_lap, sparse_lap, delta = laplacian(dense_graph), laplacian(sparse_graph), F(1, 3)
    for sign in (-1, 1):
        certificate = [[delta * original_lap[i][j] + sign * (sparse_lap[i][j] - original_lap[i][j])
                        for j in range(8)] for i in range(8)]
        assert psd(certificate)
    def simplex(m):
        return [tuple(F(u, sum(raw)) for u in raw) for raw in product((0, 1, 2), repeat=m) if sum(raw)]

    spectral_checks = 0
    for left, right in product(simplex(4), repeat=2):
        error = sum(left[i] * (sparse_graph[i][j] - 1) * right[j] for i in range(4) for j in range(4))
        assert error ** 2 <= 16 * delta ** 2 * sum(x * x for x in left) * sum(x * x for x in right)
        spectral_checks += 1
    heavy_light_checks = 0
    m, h = 3, F(3, 2)
    matrices = [tuple(tuple(F(((3 * i + j + phase) % 5) - 2, 2) for j in range(m)) for i in range(m))
                for phase in range(3)]
    for matrix in matrices:
        # These signed controls use valid scaled spectral approximants, without
        # claiming compression for every small matrix or a dynamic algorithm.
        surrogate = [[(1 + delta) * max(F(0), matrix[i][j]) - (1 - delta) * max(F(0), -matrix[i][j])
                      for j in range(m)] for i in range(m)]
        for left, right in product(simplex(m), repeat=2):
            left_heavy, right_heavy = {i for i, x in enumerate(left) if x > h / m}, {j for j, x in enumerate(right) if x > h / m}
            exact, approximation = F(0), F(0)
            for i, j in product(range(m), repeat=2):
                exact += left[i] * matrix[i][j] * right[j]
                approximation += left[i] * (matrix[i][j] if i in left_heavy or j in right_heavy else surrogate[i][j]) * right[j]
            left_light = [x if i not in left_heavy else F(0) for i, x in enumerate(left)]
            right_light = [x if j not in right_heavy else F(0) for j, x in enumerate(right)]
            assert abs(exact - approximation) <= delta * h
            assert (exact - approximation) ** 2 <= delta ** 2 * m ** 2 * sum(x * x for x in left_light) * sum(x * x for x in right_light)
            assert len(left_heavy) < m / h and len(right_heavy) < m / h
            heavy_light_checks += 1
    asymptotic_balances = []
    for j in range(1, 5):
        m, h_integer = 2 ** (6 * j), 2 ** (2 * j - 2)
        assert m * m // h_integer == 64 * m * h_integer * h_integer
        asymptotic_balances.append(dict(side=m, capacity=m * m, threshold=h_integer,
                                         exact_heavy_term=m * m // h_integer,
                                         sparsifier_leading_term_without_polylogs=64 * m * h_integer * h_integer))
    spectral_transfer = {'actual_compression_control_vertices': 8, 'original_edges': 16, 'surrogate_edges': 12,
                         'relative_spectral_error': '1/3', 'exact_psd_certificates': 2,
                         'complete_simplex_pair_controls': spectral_checks, 'signed_heavy_light_controls': heavy_light_checks,
                         'symbolic_balance_controls_at_error_one_eighth': asymptotic_balances,
                         'scope': 'Uniform algebraic transfer: supplied valid relative spectral approximants of the positive and negative payload graphs imply bilinear-read error <=delta*m*||p||2*||t||2 for independent bounded values. Exact heavy row/column scans at probability threshold h/m reduce the remaining error to delta*h and cost to O(m^2/h+S). This is a conditional classical-sparsification consequence, not a maintained dynamic or finite-arena algorithm, a finite runtime advantage, a novelty claim or the goal result. All payloads remain stored; fixed full Cartesian addresses are assumed. Adaptive updates, horizon, sparse active sets, finite grid weights and strongest method comparison require further work.'}
    fixed_sample_attacks = []
    for side in (4, 8, 16):
        sampled = {(2 * a, 2 * b) for a in range(side // 2) for b in range(side // 2) if (a + b) % 2 == 0}
        payload = [[0 for _ in range(side)] for _ in range(side)]
        untouched = 0
        for a, b in product(range(side // 2), repeat=2):
            block = {(2 * a + i, 2 * b + j) for i, j in product(range(2), repeat=2)}
            if block & sampled:
                continue
            untouched += 1
            for i, j in product(range(2), repeat=2):
                payload[2 * a + i][2 * b + j] = 1 if i == j else -1
        assert len(sampled) == side * side // 8 and untouched == side * side // 8
        assert all(payload[i][j] == 0 for i, j in sampled)
        assert all(sum(row) == 0 for row in payload)
        assert all(sum(payload[i][j] for i in range(side)) == 0 for j in range(side))
        output_lower = F(4 * untouched, side * side) * ((lam_lo - 1) / (lam_lo + 1)) ** 2
        assert output_lower > F(9, 50)
        fixed_sample_attacks.append({'side': side, 'sampled_entries': len(sampled), 'untouched_two_by_two_blocks': untouched,
                                     'sampled_values_and_all_row_column_sums_match_zero_matrix': True,
                                     'dense_dyadic_logit_query_complete_output_above': '9/50',
                                     'any_indistinguishable_reader_worst_case_error_above': '9/100'})
    from bisect import bisect_right
    cdf_cases, enumerated_draws, laws = 0, 0, []
    for side in (2, 3, 5):
        for raw in product((0, 1, 2), repeat=side):
            total = sum(raw)
            if not total:
                continue
            truth = tuple(F(x, total) for x in raw)
            for bits in (3, 5, 7):
                cumulative, cuts = 0, []
                for value in raw:
                    cumulative += value
                    cuts.append((cumulative << bits) // total)
                assert cuts[-1] == 1 << bits
                counts = [0] * side
                for draw in range(1 << bits):
                    counts[bisect_right(cuts, draw)] += 1
                    enumerated_draws += 1
                law = tuple(F(x, 1 << bits) for x in counts)
                assert sum(law) == 1 and sum(abs(x - y) for x, y in zip(law, truth)) <= F(2 * side, 1 << bits)
                if len(laws) < 20:
                    laws.append((truth, law))
                cdf_cases += 1
    product_law_checks = 0
    for (left, approx_left), (right, approx_right) in product(laws, repeat=2):
        error = sum(abs(x * y - a * b) for x, a in zip(left, approx_left) for y, b in zip(right, approx_right))
        assert error <= sum(abs(x - a) for x, a in zip(left, approx_left)) + sum(abs(y - b) for y, b in zip(right, approx_right))
        product_law_checks += 1
    binomial_tail = F(2 * sum(comb(32, j) for j in range(8)), 2 ** 32)
    _, exp4_hi = exp_interval(F(4))
    assert binomial_tail <= 2 / exp4_hi
    fresh_sampling = {'dyadic_cdf_cases': cdf_cases, 'all_fair_bit_draws_enumerated': enumerated_draws,
                       'complete_product_law_checks': product_law_checks,
                       'exact_binomial_tail_at_32_samples_and_error_half': str(binomial_tail),
                       'finite_horizon_samples_suffice': 'K>=8*epsilon^-2*log(2*v*H/rho)',
                       'scope': 'Classical fresh product-pair sampling gives the complete bounded-value mixture after charged query setup, finite CDF rounding and exact integer averaging. Conditional Hoeffding and a declared finite-horizon union bound allow adaptive values/queries chosen before each read draws new independent bits. It is stronger than the provisional spectral route in that probability model, with O(v) overwrites; no fixed-K infinite-horizon guarantee, changed active-set sampler, new theorem or ML advantage is claimed.'}
    value_width = 64
    persistent = capacity * value_width * 2 + capacity * weight * 4 + comb(dimension, 2) * 4
    persistent += capacity * 4 + (capacity + 7) // 8 + dimension * 4
    workspace = dimension * 4 + (weight + 1) * 2 + 4 * weight * 4 + 8 * value_width + (value_width + 1) * 16
    workspace += 16 * 32 + 16 * 16 + 128 * 8 + 2048
    return {'coordinate_aggregation_equivalence': coordinate_cases,
            'coefficient_mass_proof_clarification': {'support_size': 1, 'degree': 24, 'amplitude_bound': '3/8',
                                                    'query_coordinate': 12, 'actual_key_coordinate': '11/32',
                                                    'correct_mass_upper_bound': 128,
                                                    'scope': 'Coefficients use A*q, so their unweighted absolute sum is bounded by exp(A*sum|q|), not necessarily exp(sum|q*x|). The displayed intermediate estimate in the preceding note needs this clarification; its final C bound and all selected precisions remain valid.'},
            'recall_dictionary': {'declared_capacity': capacity, 'dimension': dimension, 'maximum_coordinate_load': load,
                                  'maximum_high_score_distractors': high_distractors, 'self_score_before_noise': str(self_exponent),
                                  'positive_exp_series_terms': 80, 'positive_exp_lower_numerator': self_exp_lower.numerator,
                                  'positive_exp_lower_denominator': self_exp_lower.denominator,
                                  'uniform_complete_softmax_error_below': '1/65536', 'noise_radius': '1/4',
                                  'query_support_restriction': None, 'heavy_coordinate_threshold': 1,
                                  'dense_noise_grid_units': step, 'controls': dict(index.stats), 'atomic_invalid_queries': rejected_queries,
                                  'persistent_array_bytes_at_value_width_64': persistent, 'operation_workspace_bytes': workspace,
                                  'sequential_arena_bytes_at_value_width_64': persistent + workspace,
                                  'scope': 'A checked classical dictionary shortcut is uniformly more accurate than required on the entire exhibited recall noise ball, including dense noise. It keeps every independent value, raw pair-admission state and all metadata. Outside the certified ball it pays a complete finite-precision scan, so this is not a dominance theorem over the hybrid for all mixture queries. The Python oracle installs only 65 keys, not the declared capacity or a packed arena.'},
            'spectral_value_transfer': spectral_transfer,
            'fresh_product_sampling_baseline': fresh_sampling,
            'fixed_sample_adaptive_attack': {'cases': fixed_sample_attacks,
                                             'scope': 'After a fixed entry sample is exposed, disjoint unsampled checkerboard blocks keep every sampled entry and every row/column sum zero while a dense product query with even logits 4, odd logits 0 and A=3/8 has output >9/50. The probabilities are below 2/m, so a heavy rule h>=2 reads no exact rows or columns. This invalidates uniform guarantees for readers using only that exposed sample and those sums under adaptive value edits; it does not refute oblivious-sequence matrix concentration or value-dependent/adaptively refreshed sketches, and is not a general data-structure lower bound.'},
            'mixture_gate': {'dimension': 28, 'key_weight': 8, 'key_nonzero_amplitude': '3/8', 'query_nonzero_amplitude': 4,
                             'key_supports': supports, 'independent_scalar_value_contrast': contrast,
                             'zero_and_unary_query_exact_numerators': unary_cancellations,
                             'joint_query_complete_output_above': '9/25',
                             'worst_case_top_one_all_value_error_above': '11/18',
                             'finite_l1_duality_checks': duality_checks,
                             'scope': 'Elementary exact witness against mean/unary-output reconstruction and a one-value answer at the same accuracy. Standard univariate sums plus one residual candidate already handle it; it is not an algorithmic separation or a new hard family.'},
            'scope': 'The k=2 cache is ordinary univariate kernel aggregation plus the published sparse-residual template. The displayed recall family is also solved more cheaply by a checked dictionary, with full Euclidean rather than sparse-only noise. A conditional classical spectral-transfer lemma opens a distinct value-dependent mixture route, with exact controls but without a maintained dynamic architecture claim. These are stronger baseline evidence and workload gates, not a new contribution or completion of the active goal.'}


def active_mask_sampling_controls():
    """Finite conditioning controls and a sparse-key constant-error barrier."""
    from math import comb, factorial
    from bisect import bisect_right

    def clog(n):
        return (max(1, n) - 1).bit_length()

    def nearest(n, d):
        q, r = divmod(n, d)
        return q + (2 * r >= d)

    @lru_cache(maxsize=None)
    def exp_interval(score, precision=72):
        u, term, total, j = abs(score), F(1), F(1), 0
        while True:
            next_term = term * u / (j + 1)
            if j + 2 > u:
                tail = next_term / (1 - u / (j + 2))
                lo, hi = (total, total + tail) if score >= 0 else (1 / (total + tail), 1 / total)
                if hi - lo <= F(1, 2 ** precision):
                    return lo, hi
            j += 1; term = next_term; total += term

    hardness_cases, total_boolean_queries, raw_dot_checks = [], 0, 0
    for m in (1, 2, 3):
        side, dimension = 2 * m + 1, 4 * m + 2
        capacity = side * side
        level = 2 * clog(16 * (2 * m * m + 1))
        threshold = 1 << (level // 2)
        tlo, thi = exp_interval(F(level, 2), 96)
        powers_lo, powers_hi = [tlo ** j for j in range(5)], [thi ** j for j in range(5)]
        queries = negative = positive = 0
        for mask in range(1 << (m * m)):
            edges = [(a, b) for a, b in product(range(m), repeat=2) if mask & (1 << (a * m + b))]
            active = edges + [(a, b) for a, b in product(range(m, 2 * m), repeat=2)] + [(2 * m, 2 * m)]
            n = len(active)
            assert F(n, capacity) >= F(1, 9)
            assert (n - 1) * 16 < threshold
            for left, right in product(product((0, 1), repeat=m), repeat=2):
                counts = [m * m, 0, 0]
                for a, b in edges:
                    counts[left[a] + right[b]] += 1
                truth = any(left[a] and right[b] for a, b in edges)
                assert truth == (counts[2] > 0)
                low_num = sum(counts[j] * powers_lo[2 * j] for j in range(3))
                high_num = sum(counts[j] * powers_hi[2 * j] for j in range(3))
                lo, hi = low_num / (high_num + powers_hi[3]), high_num / (low_num + powers_lo[3])
                if truth:
                    assert lo > F(15, 16) and lo - F(1, 4) > F(1, 2)
                    positive += 1
                else:
                    assert hi < F(1, 16) and hi + F(1, 4) < F(1, 2)
                    negative += 1
                # Independent raw vectors for a representative matrix catalog.
                if mask in (0, 1, (1 << (m * m)) - 1):
                    query = tuple(2 * level * x for x in left) + (0,) * m + (3 * level // 2,)
                    query += tuple(2 * level * x for x in right) + (0,) * m + (3 * level // 2,)
                    assert len(query) == dimension and max(query) <= 2 * level
                    actual_counts, sentinel = [0, 0, 0], 0
                    supports = []
                    for a, b in active:
                        raw = [F(0)] * dimension
                        raw[a] = raw[side + b] = F(1, 2)
                        assert sum(x * x for x in raw) == F(1, 2)
                        score = sum(x * q for x, q in zip(raw, query))
                        if (a, b) == (2 * m, 2 * m):
                            assert score == F(3 * level, 2)
                            sentinel += 1
                        else:
                            assert score.denominator == 1 and score % level == 0
                            actual_counts[int(score // level)] += 1
                        supports.append({a, side + b})
                        raw_dot_checks += 1
                    assert actual_counts == counts and sentinel == 1
                    assert all(len(a & b) <= 1 for a, b in combinations(supports, 2))
                queries += 1
        total_boolean_queries += queries
        hardness_cases.append(dict(boolean_side=m, raw_dimension=dimension, full_cartesian_capacity=capacity,
                                   matrix_masks=1 << (m * m), complete_queries=queries, yes_queries=positive,
                                   no_queries=negative, integer_score_step=level, maximum_raw_query_coordinate=2 * level,
                                   minimum_active_density='1/9', all_keys_support=2, squared_key_norm='1/2',
                                   scalar_error_one_quarter_still_decides=True))

    class MaskedProduct:
        def __init__(self, m=2, query_bound=2):
            self.m, self.r, self.p, self.g, self.q = m, 16, 4, 4, query_bound
            self.active, self.values, self.n = [False] * (m * m), [(0, 0)] * (m * m), 0
            self.max_precision = self.p + clog(m * m) + 3 * query_bound + clog(4096 * m)
            self.word_bits = 6 * self.max_precision + 8 + clog(m * m) + clog(query_bound + 1) + self.g + 32
            self.stats = dict(plans=0, empty_plans=0, nonempty_law_checks=0, proposal_cell_checks=0,
                              max_precision=0, max_intermediate_bits=0, max_zeta_exponent=0,
                              admitted=0, rejected_keys=0, overwrites=0, deletes=0)

        def guard(self, *values):
            bits = max(abs(x).bit_length() + 1 for x in values)
            self.stats['max_intermediate_bits'] = max(self.stats['max_intermediate_bits'], bits)
            assert bits <= self.word_bits

        def snapshot(self):
            return tuple(self.active), tuple(self.values), self.n

        def raw_key(self, a, b):
            return tuple(self.r // 2 if i in (a, self.m + b) else 0 for i in range(2 * self.m))

        def admit(self, raw, value):
            assert len(value) == 2 and all(type(x) is int and abs(x) <= 1 << self.g for x in value)
            valid = len(raw) == 2 * self.m and all(type(x) is int and x in (0, self.r // 2) for x in raw)
            support = [i for i, x in enumerate(raw) if x] if valid else []
            valid = valid and len(support) == 2 and support[0] < self.m <= support[1]
            slot = support[0] * self.m + support[1] - self.m if valid else None
            if not valid or self.active[slot]:
                self.stats['rejected_keys'] += 1
                return None
            self.active[slot], self.values[slot] = True, tuple(value)
            self.n += 1; self.stats['admitted'] += 1
            return slot

        def overwrite(self, slot, value):
            assert self.active[slot] and len(value) == 2 and all(type(x) is int and abs(x) <= 1 << self.g for x in value)
            self.values[slot] = tuple(value); self.stats['overwrites'] += 1

        def remove(self, slot):
            assert self.active[slot]
            self.active[slot], self.values[slot] = False, (0, 0)
            self.n -= 1; self.stats['deletes'] += 1

        def exponential(self, score, precision):
            assert score <= 0
            if score == 0:
                return 1 << precision
            if score < -(precision + 2):
                return 0
            h = max(1, min(self.q, precision + 2))
            cap = 1 << ((3 * h + 1) // 2)
            terms = max(6 * h, 2 * (precision + 8))
            frac = precision + clog(32 * (terms + 2) * cap)
            term, total, u = 1 << frac, 1 << frac, -score
            for j in range(1, terms + 1):
                self.guard(term * u.numerator, u.denominator * j)
                term = nearest(term * u.numerator, u.denominator * j)
                total += term; self.guard(total)
            self.guard(1 << (precision + frac))
            value = nearest(1 << (precision + frac), total)
            lo, hi = exp_interval(score)  # Verification-only; not sampler state.
            assert F(value - 1, 1 << precision) <= lo and hi <= F(value + 1, 1 << precision)
            return value

        def plan(self, query):
            if len(query) != 2 * self.m or any(type(x) is not int or abs(x) > self.q * self.r for x in query):
                raise ValueError('raw query violates dimension/grid/coordinate bound')
            self.stats['plans'] += 1
            if not self.n:
                self.stats['empty_plans'] += 1
                return None
            left, right = query[:self.m], query[self.m:]
            spread = F(max(left) - min(left) + max(right) - min(right), 2 * self.r)
            exponent = clog((self.m * self.m + self.n - 1) // self.n)
            exponent += (3 * spread.numerator + 2 * spread.denominator - 1) // (2 * spread.denominator)
            precision = self.p + exponent + clog(4096 * self.m)
            assert precision <= self.max_precision
            tables = []
            for half in (left, right):
                shift, total, prefix = max(half), 0, []
                for raw in half:
                    total += self.exponential(F(raw - shift, 2 * self.r), precision)
                    prefix.append(total)
                assert total >= 1 << precision
                for value in prefix:
                    self.guard(value << precision)
                tables.append(tuple((value << precision) // total for value in prefix))
                assert tables[-1][-1] == 1 << precision
            self.stats['max_precision'] = max(self.stats['max_precision'], precision)
            self.stats['max_zeta_exponent'] = max(self.stats['max_zeta_exponent'], exponent)
            return precision, exponent, tuple(tables)

        def proposal(self, plan, left_bits, right_bits):
            precision, _, tables = plan
            assert 0 <= left_bits < 1 << precision and 0 <= right_bits < 1 << precision
            a, b = bisect_right(tables[0], left_bits), bisect_right(tables[1], right_bits)
            slot = a * self.m + b
            return slot if self.active[slot] else None

    index, reference = MaskedProduct(), {}

    def verify(query):
        plan = index.plan(query)
        if not reference:
            assert plan is None
            return
        precision, exponent, tables = plan
        widths = [tuple(x - y for x, y in zip(table, (0,) + table[:-1])) for table in tables]
        denominator, numerator = 0, [0, 0]
        lo_den, hi_den, lo_num, hi_num = F(0), F(0), [F(0), F(0)], [F(0), F(0)]
        for a, b in product(range(index.m), repeat=2):
            slot, mass = a * index.m + b, widths[0][a] * widths[1][b]
            if mass:
                sample = index.proposal(plan, 0 if a == 0 else tables[0][a - 1], 0 if b == 0 else tables[1][b - 1])
                assert sample == (slot if slot in reference else None)
                index.stats['proposal_cell_checks'] += 1
            if slot not in reference:
                continue
            raw, values = reference[slot]
            score = F(sum(x * y for x, y in zip(raw, query)), index.r * index.r)
            lo, hi = exp_interval(score)
            lo_den += lo; hi_den += hi
            denominator += mass
            numerator = [x + mass * v for x, v in zip(numerator, values)]
            for j, value in enumerate(values):
                lo_num[j] += min(lo * value, hi * value)
                hi_num[j] += max(lo * value, hi * value)
        assert denominator > 0 and F(denominator, 1 << (2 * precision)) >= F(1, 1 << (exponent + 1))
        index.guard(denominator, *numerator)
        for j, value in enumerate(numerator):
            answer = F(value, denominator * (1 << index.g))
            endpoints = [x / y / (1 << index.g) for x in (lo_num[j], hi_num[j]) for y in (lo_den, hi_den)]
            tolerance = F(1, 1 << (index.p + 4))
            assert answer - tolerance <= min(endpoints) and max(endpoints) <= answer + tolerance
        assert tuple(index.active) == tuple(i in reference for i in range(index.m * index.m))
        assert all(index.values[i] == reference[i][1] for i in reference)
        index.stats['nonempty_law_checks'] += 1

    queries = [tuple(x * index.r for x in q) for q in product((-2, 0, 2), repeat=4)]
    for mask in range(16):
        for slot in list(reference):
            index.remove(slot); reference.pop(slot)
        for slot in range(4):
            if mask & (1 << slot):
                raw, value = index.raw_key(slot // 2, slot % 2), (5 * slot - 8, 8 - 5 * slot)
                assert index.admit(raw, value) == slot
                reference[slot] = raw, value
        for query in queries:
            verify(query)
    assert len(reference) == 4
    for slot in range(4):
        index.overwrite(slot, (16 - slot, slot - 16))
        reference[slot] = reference[slot][0], index.values[slot]
        verify(queries[-1])
        index.remove(slot); reference.pop(slot); verify(queries[0])
        raw = index.raw_key(slot // 2, slot % 2)
        assert index.admit(raw, (0, 0)) == slot
        reference[slot] = raw, (0, 0); verify(queries[-1])
    before = index.snapshot()
    assert index.admit(reference[0][0], (0, 0)) is None
    assert index.admit((8, 8, 0, 0), (0, 0)) is None
    assert index.admit((7, 0, 8, 0), (0, 0)) is None
    assert index.snapshot() == before
    invalid_queries = 0
    for query in ((), (33, 0, 0, 0), (0.5, 0, 0, 0)):
        try:
            index.plan(query)
            assert False
        except ValueError:
            invalid_queries += 1
        assert index.snapshot() == before
    saved_index, saved_reference = index, reference
    index, reference = MaskedProduct(query_bound=16), {}
    raw, value = index.raw_key(1, 0), (16, -16)
    assert index.admit(raw, value) == 2
    reference[2] = raw, value
    rare_query = (16 * index.r, 0, 0, 16 * index.r)
    verify(rare_query)
    exp8_lo, _ = exp_interval(F(8))
    assert exp8_lo > 256
    # The right half's first probability is below 2^-8, so an eight-bit
    # cumulative floor deletes the only active edge's mass entirely.
    assert (256 / (1 + exp8_lo)).numerator // (256 / (1 + exp8_lo)).denominator == 0
    rare_conditioning = {'side': 2, 'active_edges': [[1, 0]], 'active_density': '1/4',
                         'raw_physical_query': [16, 0, 0, 16], 'key_nonzero_amplitude': '1/2',
                         'naive_fractional_cdf_bits': 8, 'naive_conditional_normalizer_is_zero': True,
                         'exact_conditional_output': [1, -1],
                         'corrected_maximum_precision': index.stats['max_precision'],
                         'dyadic_acceptance_lower_bound_exponent': index.stats['max_zeta_exponent'],
                         'scope': 'A constant-density mask can select an extremely rare product event. Marginal CDF rounding that is adequate without conditioning can erase all active mass. The corrected precision preserves the conditional law, but its charged rejection bound is expensive; a complete scan is the appropriate finite fallback.'}
    index, reference = saved_index, saved_reference
    retention_certificates = []
    for log_capacity in (4, 16, 32, 64):
        capacity, grid_bits = 1 << log_capacity, 8
        beta = 4 * clog(32 * capacity * (1 << grid_bits))
        # e^(beta/4)>2^(beta/4), e^(1/2)<2. Noise radius1/4,
        # key norm<=1 and pair overlap<=1 give gap>=beta/4-1/2.
        assert F(4 * capacity, 1 << (beta // 4)) < F(1, 1 << (grid_bits + 1))
        retention_certificates.append(dict(log2_capacity=log_capacity, value_fraction_bits=grid_bits,
                                           query_scale=beta, euclidean_noise_radius='1/4',
                                           independent_nearest_grid_recall_follows=True))
    return {'sparse_key_mask_hardness': {'cases': hardness_cases, 'total_boolean_queries': total_boolean_queries,
                                        'independent_raw_dot_product_checks': raw_dot_checks,
                                        'scope': 'The elementary sentinel reduction gives a constant-error bounded-word OuMv corollary, with weight-two norm-bounded keys, integer raw queries of logarithmic amplitude and constant-density nonempty masks. It verifies the construction and numeric gap, not the conjecture or novelty. No key edits or online value changes are needed in the reduction.'},
            'finite_conditioned_product_laws': {'side': index.m, 'raw_fraction_bits': 4, 'value_fraction_bits': index.g,
                                               'accuracy_bits': index.p, 'declared_query_coordinate_bound': index.q,
                                               'preselected_maximum_fractional_precision': index.max_precision,
                                               'preselected_intermediate_guard_bits': index.word_bits,
                                               'controls': index.stats, 'atomic_invalid_queries': invalid_queries,
                                               'scope': 'Plans form only O(m) finite CDF tables. Verification independently enumerates all active cells to check the conditional law and true exponential mean; that enumeration is not part of the sampling algorithm. Density and full-score-spread determine a dyadic acceptance lower bound, charged proposal count and extra precision. Raw admission, canonical overwrite, deletion and zero-value reuse are checked. Conditional Hoeffding/horizon analysis is reused from the preceding fresh-sampling proof; no stochastic accuracy simulation or native arena claim.'},
            'independent_retention_controls': retention_certificates,
            'rare_conditioning_precision_witness': rare_conditioning,
            'scope': 'Classical conditioning, concentration and an OuMv sentinel transfer define a useful baseline boundary: bounded density alone does not keep masked product attention easy when raw logit spread grows logarithmically. They do not establish a novel theorem beyond published methods or complete the research goal.'}


def finite_radon_attention_controls():
    """Exact characteristic-two Radon/transpose and complete-mask controls.

    Classical Fourier-slice machinery; not a new transform, retention
    compression, strongest-method advantage, or trained model.
    """
    class Radon:
        def __init__(self, s, dimension, polynomial):
            self.s, self.dimension, self.q = s, dimension, 1 << s
            assert s >= 1 and dimension >= 1 and polynomial.bit_length() == s + 1

            def remainder(a, b):
                while a and a.bit_length() >= b.bit_length():
                    a ^= b << (a.bit_length() - b.bit_length())
                return a

            for degree in range(1, s // 2 + 1):
                assert all(remainder(polynomial, f) for f in range(1 << degree, 1 << (degree + 1))), 'reducible field polynomial'

            def multiply(a, b):
                result = 0
                while b:
                    if b & 1:
                        result ^= a
                    b >>= 1; a <<= 1
                    if a & self.q:
                        a ^= polynomial
                return result

            self.mul = [[multiply(a, b) for b in range(self.q)] for a in range(self.q)]
            self.points = [tuple((i >> (s * j)) & (self.q - 1) for j in range(dimension))
                           for i in range(self.q ** dimension)]
            self.directions = [a for a in self.points if any(a) and next(x for x in a if x) == 1]
            assert len(self.directions) == (len(self.points) - 1) // (self.q - 1)

            def trace(a):
                total, x = 0, a
                for _ in range(s):
                    total ^= x; x = self.mul[x][x]
                assert total in (0, 1)
                return total

            self.dual = [sum(trace(self.mul[a][1 << j]) << j for j in range(s)) for a in range(self.q)]
            assert sorted(self.dual) == list(range(self.q))
            self.permutation = [sum(self.dual[x] << (s * j) for j, x in enumerate(k)) for k in self.points]
            self.slices = [tuple(sum(self.mul[t][a[j]] << (s * j) for j in range(dimension))
                                 for t in range(self.q)) for a in self.directions]
            visits = [0] * len(self.points)
            for line in self.slices:
                for k in line:
                    visits[k] += 1
            assert visits[0] == len(self.directions) and all(v == 1 for v in visits[1:])
            self.types = len(self.directions) * self.q
            self.max_bits = 0

        def transform(self, values, permutation, guard, stats):
            data, stride = list(values), 1
            while stride < len(data):
                for start in range(0, len(data), 2 * stride):
                    for j in range(start, start + stride):
                        a, b = data[j], data[j + stride]
                        data[j], data[j + stride] = a + b, a - b
                        stats['add_subtracts'] += 2
                        self.max_bits = max(self.max_bits, abs(a + b).bit_length() + 1, abs(a - b).bit_length() + 1)
                        assert self.max_bits <= guard
                stride *= 2
            return [data[i] for i in permutation]

        def forward(self, values, payload_bits=12):
            assert len(values) == len(self.points) and all(type(x) is int and abs(x) < 1 << payload_bits for x in values)
            guard = payload_bits + self.s * self.dimension + self.s + 5
            stats = {'add_subtracts': 0, 'exact_divisions': 0}
            spectrum = self.transform(values, self.permutation, guard, stats)
            answer = []
            for line in self.slices:
                bins = self.transform([spectrum[i] for i in line], self.dual, guard, stats)
                assert all(x % self.q == 0 for x in bins)
                answer.extend(x // self.q for x in bins)
                stats['exact_divisions'] += self.q
            expected = len(self.points) * self.s * self.dimension + self.types * self.s
            assert stats['add_subtracts'] == expected
            return answer, stats

        def transpose(self, values, payload_bits=12):
            assert len(values) == self.types and all(type(x) is int and abs(x) < 1 << payload_bits for x in values)
            guard = payload_bits + self.s * self.dimension + self.s + 5
            stats = {'add_subtracts': 0, 'scatter_adds': 0, 'exact_divisions': 0}
            spectrum = [0] * len(self.points)
            for j, line in enumerate(self.slices):
                values_hat = self.transform(values[j * self.q:(j + 1) * self.q], self.dual, guard, stats)
                for i, value in zip(line, values_hat):
                    spectrum[i] += value; stats['scatter_adds'] += 1
                    self.max_bits = max(self.max_bits, abs(spectrum[i]).bit_length() + 1)
                    assert self.max_bits <= guard
            answer = self.transform(spectrum, self.permutation, guard, stats)
            assert all(x % self.q == 0 for x in answer)
            stats['exact_divisions'] = len(answer)
            return [x // self.q for x in answer], stats

    def explicit_rows(operator):
        # Independent dense incidence oracle, never part of fast operator state.
        rows = []
        for a in operator.directions:
            bins = [[] for _ in range(operator.q)]
            for i, x in enumerate(operator.points):
                value = 0
                for u, v in zip(a, x):
                    value ^= operator.mul[u][v]
                bins[value].append(i)
            rows.extend(bins)
        return rows

    cases, basis_columns, basis_rows, integer_outputs = [], 0, 0, 0
    choices = [(1, d, 3) for d in range(1, 6)] + [(2, d, 7) for d in range(1, 4)]
    choices += [(3, 1, 11), (3, 2, 11), (3, 3, 11), (4, 2, 19)]
    for s, dimension, polynomial in choices:
        operator = Radon(s, dimension, polynomial)
        rows = explicit_rows(operator)
        n, b = len(operator.points), operator.types
        values = [(37 * i + 11) % 2001 - 1000 for i in range(n)]
        covector = [(29 * i + 19) % 2001 - 1000 for i in range(b)]
        actual, forward_stats = operator.forward(values)
        expected = [sum(values[i] for i in row) for row in rows]
        assert actual == expected
        back, transpose_stats = operator.transpose(covector)
        expected_back = [sum(covector[j] for j, row in enumerate(rows) if i in row) for i in range(n)]
        assert back == expected_back and sum(a * c for a, c in zip(actual, covector)) == sum(v * z for v, z in zip(values, back))
        integer_outputs += n + b
        if n <= 64:
            for i in range(n):
                unit = [int(j == i) for j in range(n)]
                assert operator.forward(unit)[0] == [int(i in row) for row in rows]
                basis_columns += 1
            for j in range(b):
                unit = [int(k == j) for k in range(b)]
                assert operator.transpose(unit)[0] == [int(i in rows[j]) for i in range(n)]
                basis_rows += 1
        a_count = len(operator.directions)
        off_diagonal = (operator.q ** (dimension - 1) - 1) // (operator.q - 1)
        gram = operator.transpose(actual, payload_bits=20)[0]
        assert gram == [operator.q ** (dimension - 1) * v + off_diagonal * sum(values) for v in values]
        for i in range(n):
            assert sum(int(i in row) for row in rows) == a_count
        cases.append({'q': operator.q, 'ambient_dimension': dimension, 'irreducible_polynomial_bits': polynomial,
                      'profiles': n, 'types': b, 'direct_state_additions': a_count * n,
                      'forward_add_subtracts': forward_stats['add_subtracts'],
                      'forward_exact_divisions': forward_stats['exact_divisions'],
                      'transpose_add_subtracts': transpose_stats['add_subtracts'],
                      'transpose_scatter_adds': transpose_stats['scatter_adds'],
                      'observed_maximum_signed_bits': operator.max_bits,
                      'gram_and_independent_profile_invertibility_checked': True,
                      'basis_operator_checks': n <= 64,
                      'nontrivial_trace_dual_permutation': operator.dual != list(range(operator.q))})

    malformed_polynomials = 0
    for s, polynomial in [(2, 5), (3, 9), (4, 21)]:
        try:
            Radon(s, 2, polynomial)
            assert False
        except AssertionError as error:
            assert str(error) == 'reducible field polynomial'
            malformed_polynomials += 1

    attention_cases = []
    for s, dimension, polynomial in [(1, 2, 3), (2, 2, 7), (2, 3, 7)]:
        operator = Radon(s, dimension, polynomial)
        rows = explicit_rows(operator)
        capacity, types, r, v, scale = len(operator.points), operator.types, 3, 2, 8
        distant, recent = capacity + 1, 2 * types
        length = distant + recent
        profiles = [j % capacity for j in range(distant)]
        routing = [j % types for j in range(recent)]
        psi = [[1 + ((profiles[j] if j < distant else j) + 2 * t) % 5 for t in range(r)] for j in range(length)]
        phi = [[1 + (3 * j + t) % 7 for t in range(r)] for j in range(length)]
        payload = [[(j * (c + 3) + 2 * c) % 17 - 8 for c in range(v)] for j in range(length)]
        states = [[psi[j][t] * x for t in range(r) for x in payload[j] + [scale]] for j in range(length)]
        profiles_table = [[sum(states[j][c] for j in range(distant) if profiles[j] == x)
                           for c in range(r * (v + 1))] for x in range(capacity)]
        type_table = [[0] * (r * (v + 1)) for _ in range(types)]
        for c in range(r * (v + 1)):
            transformed, _ = operator.forward([row[c] for row in profiles_table])
            for i, x in enumerate(transformed):
                type_table[i][c] = x
        prefix = [0] * (r * (v + 1))
        mask_entries = 0
        for i in range(length):
            if i == distant:
                prefix = [0] * len(prefix)
            prefix = [a + b for a, b in zip(prefix, states[i])]
            combined = prefix if i < distant else [a + b for a, b in zip(prefix, type_table[routing[i - distant]])]
            output = [sum(phi[i][t] * combined[t * (v + 1) + c] for t in range(r)) for c in range(v + 1)]
            selected = [j for j in range(length) if (j <= i if i < distant else
                        (profiles[j] in rows[routing[i - distant]] if j < distant else j <= i))]
            weights = [sum(a * b for a, b in zip(phi[i], psi[j])) for j in selected]
            denominator = sum(weights)
            assert denominator > 0
            direct = [F(sum(w * payload[j][c] for j, w in zip(selected, weights)), scale * denominator) for c in range(v)]
            assert [F(x, output[-1]) for x in output[:-1]] == direct
            mask_entries += len(selected)
        # Equal feature vectors and opposite independent payload changes at
        # duplicate profiles cancel in every pooled/type summary.
        assert profiles[0] == profiles[-1] and psi[0] == psi[distant - 1]
        assert all(psi[0][t] - psi[distant - 1][t] == 0 for t in range(r))
        opposite_changes = [[0] * (r * (v + 1)) for _ in range(capacity)]
        for j, sign in [(0, 1), (distant - 1, -1)]:
            for t in range(r):
                opposite_changes[profiles[j]][t * (v + 1)] += sign * psi[j][t]
        assert all(x == 0 for row in opposite_changes for x in row)
        # Actual current-value overwrite: update its raw payload and profile,
        # touch every incident cached type, then independently replay all types.
        slot = distant // 2
        previous = states[slot][:]
        payload[slot] = [scale, -scale]
        states[slot] = [psi[slot][t] * x for t in range(r) for x in payload[slot] + [scale]]
        delta = [a - b for a, b in zip(states[slot], previous)]
        profiles_table[profiles[slot]] = [a + b for a, b in zip(profiles_table[profiles[slot]], delta)]
        touched = 0
        for i, row in enumerate(rows):
            if profiles[slot] in row:
                type_table[i] = [a + b for a, b in zip(type_table[i], delta)]
                touched += 1
        assert touched == len(operator.directions)
        for c in range(r * (v + 1)):
            replayed, _ = operator.forward([row[c] for row in profiles_table])
            assert replayed == [row[c] for row in type_table]
        attention_cases.append({'q': operator.q, 'ambient_dimension': dimension, 'tokens': length,
                                'query_feature_dimension': r, 'value_width': v,
                                'complete_normalized_rows_checked': length,
                                'explicit_mask_entries_checked': mask_entries,
                                'original_independent_payloads_retained': True,
                                'overwrite_incident_type_updates': touched,
                                'duplicate_profile_summary_collision_checked': True})
    prime_slice_cases = []
    for q, dimension in [(3, 2), (3, 3), (5, 2)]:
        # Verification-only exact arithmetic in Z[z]/(1+...+z^(q-1)).
        # This quadratic DFT oracle does not implement the fast prime FFT.
        points = list(product(range(q), repeat=dimension))
        directions = [a for a in points if any(a) and next(x for x in a if x) == 1]
        values = [(11 * i + 7) % 31 - 15 for i in range(len(points))]
        def reduce(coefficients):
            return tuple(c - coefficients[-1] for c in coefficients[:-1])
        spectrum = {}
        for k in points:
            coefficients = [0] * q
            for x, value in zip(points, values):
                coefficients[sum(a * b for a, b in zip(k, x)) % q] += value
            spectrum[k] = reduce(coefficients)
        checked = 0
        for a in directions:
            for b in range(q):
                coefficients = [0] * q
                for frequency in range(q):
                    k = tuple(frequency * x % q for x in a)
                    for exponent, value in enumerate(spectrum[k]):
                        coefficients[(exponent - frequency * b) % q] += value
                answer = reduce(coefficients)
                expected = sum(value for x, value in zip(points, values) if sum(u * v for u, v in zip(a, x)) % q == b)
                assert answer == (q * expected,) + (0,) * (q - 2)
                checked += 1
        prime_slice_cases.append({'prime': q, 'ambient_dimension': dimension, 'profiles': len(points),
                                  'exact_cyclotomic_output_checks': checked,
                                  'oracle_is_quadratic_not_a_fast_fft_implementation': True})
    return {'operator_cases': cases, 'complete_basis_columns': basis_columns, 'complete_basis_rows': basis_rows,
            'signed_integer_outputs': integer_outputs, 'reducible_field_polynomial_rejections': malformed_polynomials,
            'complete_attention_cases': attention_cases, 'prime_field_slice_controls': prime_slice_cases,
            'scope': 'Classical trace Fourier-slice/Walsh specialization exactly matches finite-field incidence multiplication and its transpose. Finite basis and signed controls, Gram identities and complete causal normalized finite-feature reads are checked. This is not a new Radon transform, softmax kernel approximation guarantee, learned-selector speedup, native benchmark or strongest-published-method separation. Characteristic-two geometry is an explicit field extension of the prime-field SMat presentation. All independent original payloads must remain charged; pooled summaries alone collide on duplicate profiles.'}


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
              "subset_aggregate_candidate": subset_aggregate_controls(),
              "affine_coded_baseline": affine_coded_controls(),
              "weighted_exponential_candidate": weighted_exponential_controls(),
              "capacity_frontier_baselines": capacity_frontier_controls(),
              "continuous_amplitude_baseline": continuous_amplitude_controls(),
              "hybrid_fixed_point_candidate": hybrid_moment_controls(),
              "exponential_aggregation_comparators": exponential_comparator_controls(),
              "active_mask_sampling_boundary": active_mask_sampling_controls(),
              "finite_radon_attention_boundary": finite_radon_attention_controls(),
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
