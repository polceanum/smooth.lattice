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
