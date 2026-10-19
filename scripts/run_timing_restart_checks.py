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
