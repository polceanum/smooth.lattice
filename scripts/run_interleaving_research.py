#!/usr/bin/env python3
"""Exact exponent-stream checks and the bounded rounded-occurrence study."""
from __future__ import annotations

import argparse
from bisect import bisect_right
from collections import Counter
from fractions import Fraction
import hashlib
import heapq
from itertools import combinations, product
import json
from math import factorial, gcd, isqrt, lcm, prod
from pathlib import Path
import shlex
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
MODELS = {
    "primes_2_3": (2, 3),
    "primes_2_3_5": (2, 3, 5),
    "primes_2_3_5_7": (2, 3, 5, 7),
    "squared_2_3": (4, 9),
    "squared_2_3_5": (4, 9, 25),
    "deformed_2_3_5": (4, 9, 125),
}


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def exact_prefix(generators, count):
    """All configured generators are powers of distinct primes: no value ties."""
    zero = (0,) * len(generators)
    heap = [(1, zero)]
    seen = {zero}
    out = []
    while len(out) < count:
        value, exps = heapq.heappop(heap)
        out.append((value, exps))
        for i, base in enumerate(generators):
            child = list(exps)
            child[i] += 1
            child = tuple(child)
            if child not in seen:
                seen.add(child)
                heapq.heappush(heap, (value * base, child))
    return out


def pointer_values(generators, count):
    """Independent exact pointer oracle for the enumeration, on a short prefix."""
    values = [1]
    pointers = [0] * len(generators)
    while len(values) < count:
        candidates = [p * values[j] for p, j in zip(generators, pointers)]
        value = min(candidates)
        values.append(value)
        for i, candidate in enumerate(candidates):
            if candidate == value:
                pointers[i] += 1
    return values


def check_trimming(word):
    counts = []
    for i in range(len(word[0])):
        trimmed = []
        for exps in word:
            if exps[i] > 0:
                parent = list(exps)
                parent[i] -= 1
                trimmed.append(tuple(parent))
        require(trimmed == word[:len(trimmed)], ("coordinate trimming", i))
        counts.append(len(trimmed))
    return counts


def check_embeddings(word):
    positions = {exps: n for n, exps in enumerate(word, 1)}
    k = len(word[0])
    checks = 0
    for exps in word:
        for i in range(k):
            for j in range(i + 1, k):
                child = list(exps)
                child[i] += 1
                child[j] += 1
                if tuple(child) not in positions:
                    continue
                left = list(exps)
                left[i] += 1
                right = list(exps)
                right[j] += 1
                # Both intermediate ranks must exist, not just the final vector.
                ni = positions[tuple(left)]
                nj = positions[tuple(right)]
                left = list(word[ni - 1])
                left[j] += 1
                right = list(word[nj - 1])
                right[i] += 1
                require(positions[tuple(left)] == positions[tuple(right)],
                        ("embedding commutation", exps, i, j))
                checks += 1
    return checks


def occurrence_transform(word):
    """Replace each entry by its number of strictly earlier occurrences."""
    counts = {}
    out = []
    for entry in word:
        previous = counts.get(entry, 0)
        out.append(previous)
        counts[entry] = previous + 1
    return out


def is_lattice_word(word):
    """Every prefix has count(0) >= count(1) >= count(2) >= ... ."""
    counts = {}
    for entry in word:
        counts[entry] = counts.get(entry, 0) + 1
        if entry > 0 and counts.get(entry - 1, 0) < counts[entry]:
            return False
    return True


def equality_patterns(length, prefix=(0,)):
    """One word per set partition, labelling symbols by first appearance."""
    if len(prefix) == length:
        yield prefix
        return
    for entry in range(max(prefix) + 2):
        yield from equality_patterns(length, prefix + (entry,))


def check_occurrence_dynamics(max_length=10):
    require(occurrence_transform([0, 0, 0, 0]) == [0, 1, 2, 3],
            "constant to occurrence indices")
    require(occurrence_transform([0, 1, 2, 3]) == [0, 0, 0, 0],
            "distinct entries to zeroes")
    checked = 0
    periods = set()
    for length in range(1, max_length + 1):
        for pattern in equality_patterns(length):
            original = list(pattern)
            first = occurrence_transform(original)
            second = occurrence_transform(first)
            third = occurrence_transform(second)
            require(third == first, ("T cubed equals T", original))
            require((second == original) == is_lattice_word(original),
                    ("two-pass return criterion", original))
            periods.add(1 if first == second else 2)
            checked += 1
    return {
        "max_exhaustive_length": max_length,
        "equality_patterns_checked": checked,
        "observed_cycle_lengths": sorted(periods),
        "universal_identity": "T(T(T(word))) == T(word)",
        "nonreturning_input_example": [0, 1, 1],
        "nonreturning_input_orbit": [[0, 1, 1], [0, 0, 1], [0, 1, 0], [0, 0, 1]],
    }


def check_restricted_occurrence_cycles():
    """Reconstructed candidate rules, not identified as the user's original rule."""
    examples = [
        ("round_counts_up_to_even", [0, 2, 0, 2, 2, 2],
         lambda count: 2 * ((count + 1) // 2)),
        ("round_counts_down_to_squares", [0, 0, 1, 1, 1, 1, 1],
         lambda count: isqrt(count) ** 2),
    ]
    results = []
    for name, seed, project in examples:
        cycle = [seed]
        current = seed
        for step in range(1, 7):
            current = [project(count) for count in occurrence_transform(current)]
            if step < 6:
                require(current != seed and current not in cycle,
                        ("premature return in restricted-count candidate", name, step))
                cycle.append(current)
            else:
                require(current == seed, ("six-pass return", name, current))
        results.append({"rule": name, "minimal_period": 6,
                        "cycle": cycle, "scope": "finite word only"})
    return {
        "status": "candidate reconstructions; original operation still unknown",
        "examples": results,
    }


def floor_to_base_power(count, base):
    if base < 2 or count < 0:
        raise ValueError("base must be >= 2 and count must be nonnegative")
    if count == 0:
        return 0
    power = 1
    while power * base <= count:
        power *= base
    return power


def rounded_occurrence_transform(word, base):
    return tuple(floor_to_base_power(count, base)
                 for count in occurrence_transform(word))


def frequency_partition(word):
    return tuple(sorted(Counter(word).values(), reverse=True))


def grouped_column_histogram(shape, project):
    """Independent Ferrers-cell oracle, retaining output labels (O(sum(shape)))."""
    counts = Counter()
    for row_length in shape:
        for column in range(row_length):
            counts[project(column)] += 1
    return counts


def prefix_partition_chain(word):
    counts = Counter()
    chain = []
    for symbol in word:
        counts[symbol] += 1
        chain.append(tuple(sorted(counts.values(), reverse=True)))
    return chain


def reconstruct_counted_word(chain, project):
    """Recover q(T(word)) from its INPUT prefix partitions, with original labels."""
    previous = Counter()
    word = []
    for shape in chain:
        histogram = grouped_column_histogram(shape, project)
        difference = {label: histogram[label] - previous[label]
                      for label in histogram.keys() | previous.keys()
                      if histogram[label] != previous[label]}
        if len(difference) != 1 or next(iter(difference.values())) != 1:
            raise ValueError("prefix partitions do not give a valid one-entry histogram increment")
        word.append(next(iter(difference)))
        previous = histogram
    return tuple(word)


def base_power_partition_step(shape, base):
    """Transpose, aggregate geometric column blocks, and sort; any finite mass."""
    if base < 2:
        raise ValueError("base must be >= 2")
    if not shape:
        return ()
    counts = [len(shape)]
    lower = 1
    while lower < max(shape):
        upper = base * lower
        counts.append(sum(max(0, min(part, upper) - lower) for part in shape))
        lower = upper
    return tuple(sorted((count for count in counts if count), reverse=True))


def restricted_partition_step(shape, base):
    """Three-group formula, valid only when sum(shape) <= base**2."""
    mass = sum(shape)
    if base < 3 or not 1 <= mass <= base * base:
        raise ValueError("restricted study requires base >= 3 and 1 <= mass <= base**2")
    tail = sum(max(0, part - base) for part in shape)
    return tuple(sorted((x for x in (len(shape), mass - len(shape) - tail, tail)
                         if x), reverse=True))


def canonical_cycle(cycle):
    """Rotation only: reversing a directed cycle is not an equivalence."""
    cycle = tuple(cycle)
    start = min(range(len(cycle)), key=cycle.__getitem__)
    return cycle[start:] + cycle[:start]


def finite_orbit(seed, step, max_steps=100000):
    """Return the first repeat, separating the transient from the minimal cycle."""
    positions = {}
    states = []
    current = tuple(seed)
    while current not in positions:
        if len(states) >= max_steps:
            raise RuntimeError("orbit limit reached; no cycle has been established")
        positions[current] = len(states)
        states.append(current)
        current = tuple(step(current))
    transient = positions[current]
    return {
        "transient_length": transient,
        "minimal_period": len(states) - transient,
        "seed_is_periodic": transient == 0,
        "cycle": tuple(states[transient:]),
    }


def partitions_at_most_three(mass):
    yield (mass,)
    for second in range(1, mass // 2 + 1):
        yield (mass - second, second)
    for third in range(1, mass // 3 + 1):
        for second in range(third, (mass - third) // 2 + 1):
            yield (mass - second - third, second, third)


def restricted_partition_cycles(mass, base):
    """Closed classification from Proposition 12; does not iterate either map."""
    if base < 3 or not 1 <= mass <= base * base:
        raise ValueError("restricted study requires base >= 3 and 1 <= mass <= base**2")

    def shape(*parts):
        return tuple(sorted(parts, reverse=True))

    if mass == 1:
        cycles = [[(1,)]]
    elif mass == 2:
        cycles = [[(2,), (1, 1)]]
    elif mass <= base + 2:
        cycles = [[shape(mass - 2, 2)]]
    elif mass == base + 3:
        cycles = [[shape(base + 1, 2), shape(base, 2, 1), shape(base, 3)]]
    elif mass == base + 4:
        cycles = [[shape(base + 1, 3), shape(base + 1, 2, 1), shape(base, 3, 1)]]
    else:
        remainder = mass - 3
        if remainder <= 2 * base:
            reflection = remainder - base
            lower, upper = 1, reflection - 1
        elif remainder < 3 * base:
            excess = remainder - 2 * base
            reflection = base + excess
            lower, upper = excess, base
        else:
            return (canonical_cycle([shape(3, 2 * base, remainder - 2 * base)]),)
        cycles = []
        for parameter in range(lower, upper + 1):
            partner = reflection - parameter
            if parameter > partner:
                continue
            cycle = [shape(3, parameter, remainder - parameter)]
            if parameter != partner:
                cycle.append(shape(3, partner, remainder - partner))
            cycles.append(cycle)
    return tuple(sorted(canonical_cycle(cycle) for cycle in cycles))


def check_restricted_partition_graph(mass, base):
    """Visit EVERY reduced state, without restricting its basin or seed shape."""
    done = {}
    cycles = set()
    enumerated = 0
    for initial in partitions_at_most_three(mass):
        enumerated += 1
        if initial in done:
            continue
        trail, positions = [], {}
        current = initial
        while current not in done and current not in positions:
            positions[current] = len(trail)
            trail.append(current)
            following = base_power_partition_step(current, base)
            require(following == restricted_partition_step(current, base),
                    ("geometric groups versus three-group formula", base, current))
            current = following
        if current in positions:
            start = positions[current]
            cycle = canonical_cycle(trail[start:])
            cycles.add(cycle)
            for state in trail[start:]:
                done[state] = (0, cycle)
            trail = trail[:start]
        for state in reversed(trail):
            distance, cycle = done[current]
            done[state] = (distance + 1, cycle)
            current = state
    require(enumerated == len(done), ("complete reduced state enumeration", base, mass))
    expected = set(restricted_partition_cycles(mass, base))
    require(cycles == expected, ("closed cycle classification", base, mass, cycles, expected))
    for cycle in cycles:
        require(len(set(cycle)) == len(cycle), ("cycle minimality", base, mass))
        for index, state in enumerate(cycle):
            require(base_power_partition_step(state, base) == cycle[(index + 1) % len(cycle)],
                    ("directed cycle transition", base, mass, state))
    return {
        "mass": mass,
        "reduced_states_checked": enumerated,
        "max_reduced_transient": max(distance for distance, _ in done.values()),
        "basin_states_by_period": dict(sorted(Counter(len(cycle)
                                                       for _, cycle in done.values()).items())),
        "cycles": [{"minimal_period": len(cycle), "states": cycle}
                   for cycle in sorted(cycles)],
    }


def check_restricted_partition_classification():
    bases = []
    for base in range(3, 17):
        rows = [check_restricted_partition_graph(mass, base)
                for mass in range(1, base * base + 1)]
        bases.append({"base": base, "mass_limit": base * base,
                      "reduced_states_checked": sum(row["reduced_states_checked"] for row in rows),
                      "masses": rows})
    return {
        "scope": "all partitions with <= 3 parts, for every mass 1..base**2; "
                 "every unrestricted partition enters this reduced set in one step",
        "theorem": "Proposition 12 in docs/interleaving_research.md",
        "bases": bases,
        "total_reduced_states_checked": sum(row["reduced_states_checked"] for row in bases),
    }


def check_word_partition_link(max_length=7):
    """Word oracles test the factor formula and prefix-period lifting separately."""
    checked = transitions = 0
    witnesses = []
    for base in range(3, 17):
        project = lambda count: floor_to_base_power(count, base)
        step = lambda shape: base_power_partition_step(shape, base)
        orbit_cache = {}

        def prefix_orbits(chain):
            for shape in chain:
                if shape not in orbit_cache:
                    orbit_cache[shape] = finite_orbit(shape, step)
            return [orbit_cache[shape] for shape in chain]

        def check(word):
            nonlocal checked, transitions
            shape = frequency_partition(word)
            transformed = rounded_occurrence_transform(word, base)
            require(Counter(transformed) == grouped_column_histogram(shape, project),
                    ("labelled Ferrers-cell oracle", base, word))
            require(frequency_partition(transformed) == step(shape),
                    ("word/partition factor", base, word))
            chain = prefix_partition_chain(word)
            require(reconstruct_counted_word(chain, project) == transformed,
                    ("general prefix reconstruction", base, word))
            orbit = finite_orbit(word, lambda state: rounded_occurrence_transform(state, base))
            for index, state in enumerate(orbit["cycle"]):
                following = orbit["cycle"][(index + 1) % orbit["minimal_period"]]
                require(frequency_partition(following) == step(frequency_partition(state)),
                        ("cycle word/partition factor", base, state))
                transitions += 1
            shape_orbits = prefix_orbits(chain)
            periods = [row["minimal_period"] for row in shape_orbits]
            require(orbit["minimal_period"] == lcm(*periods),
                    ("all-prefix period lifting", base, word, periods))
            max_transient = max(row["transient_length"] for row in shape_orbits)
            require(max_transient <= orbit["transient_length"] <= max_transient + 1,
                    ("general one-step transient bound", base, word, orbit))
            orbit["max_prefix_partition_transient"] = max_transient
            if 1 <= len(word) <= base * base:
                expected = 1 if len(word) == 1 else (2 if len(word) <= base + 2 else 6)
                require(orbit["minimal_period"] == expected,
                        ("uniform eventual word period", base, word, orbit))
            checked += 1
            return orbit, periods

        for length in range(1, max_length + 1):
            for pattern in equality_patterns(length):
                check(pattern)
        # Every geometric threshold and its neighbour, including the upper endpoint.
        for count in (0, 1, base - 1, base, base + 1, base * base - 1, base * base):
            word = (0,) * (count + 1)
            project_counts = occurrence_transform(word)
            require(frequency_partition(rounded_occurrence_transform(word, base)) ==
                    step(frequency_partition(word)), ("rounding threshold", base, count))
            require(project_counts[-1] == count, ("strictly earlier count", count))
        for length in sorted({1, 2, base + 2, base + 3, base * base}):
            seed = (0,) * length
            orbit, periods = check(seed)
            witnesses.append({"base": base, "length": length, "seed": seed,
                              **orbit, "prefix_partition_periods": periods,
                              "scope": "within the proved range"})
        if base == 3:
            seed = (0,) * 64
            orbit, periods = check(seed)
            require((orbit["transient_length"], orbit["minimal_period"]) == (10, 12),
                    "known outside-range twelve-cycle after transient")
            witnesses.append({"base": base, "length": len(seed), "seed": seed,
                              **orbit, "prefix_partition_periods": periods,
                              "scope": "outside proved classification; exact finite witness only"})
        require(base_power_partition_step((), base) == () and
                rounded_occurrence_transform((), base) == (), "empty-word convention")
    # The lifting lemma needs no monotonicity or injectivity of q.
    general_cases = 0
    for name, project in (("identity", lambda c: c), ("constant", lambda c: 0),
                          ("parity", lambda c: c % 2)):
        for word in product(range(3), repeat=5):
            def general_step(shape):
                return tuple(sorted(grouped_column_histogram(shape, project).values(), reverse=True))
            orbit = finite_orbit(word, lambda state: tuple(project(c)
                                                           for c in occurrence_transform(state)))
            chain = prefix_partition_chain(word)
            shape_orbits = [finite_orbit(shape, general_step) for shape in chain]
            periods = [row["minimal_period"] for row in shape_orbits]
            require(orbit["minimal_period"] == lcm(*periods), ("general q lifting", name, word))
            max_transient = max(row["transient_length"] for row in shape_orbits)
            require(max_transient <= orbit["transient_length"] <= max_transient + 1,
                    ("general q transient bound", name, word))
            require(reconstruct_counted_word(chain, project) ==
                    tuple(project(c) for c in occurrence_transform(word)),
                    ("general q reconstruction", name, word))
            relabelled = finite_orbit(word, lambda state: tuple(10 + 7 * project(c)
                                                               for c in occurrence_transform(state)))
            expected_cycle = [tuple(10 + 7 * entry for entry in state) for state in orbit["cycle"]]
            require(canonical_cycle(relabelled["cycle"]) == canonical_cycle(expected_cycle),
                    ("injective output labels preserve the entire eventual cycle", name, word))
            general_cases += 1
    return {"max_exhaustive_equality_pattern_length": max_length,
            "base_power_words_checked": checked, "cycle_transitions_checked": transitions,
            "rounding_threshold_cases_checked": 14 * 7,
            "general_count_map_words_checked": general_cases, "witnesses": witnesses}


def check_pairwise_order_obstruction(prefix_length=10000):
    """Prime triples with long identical PAIR exponent words but different rank 5."""
    lower, upper = (101, 683, 68963), (101, 683, 68993)
    for value in set(lower + upper):
        require(value >= 2 and all(value % divisor for divisor in range(2, isqrt(value) + 1)),
                ("independent trial-division prime check", value))
    require(lower[2] < lower[0] * lower[1] < upper[2], "third primes straddle product")
    pairs = []
    for i, j in ((0, 1), (0, 2), (1, 2)):
        left = exact_prefix((lower[i], lower[j]), prefix_length)
        right = exact_prefix((upper[i], upper[j]), prefix_length)
        left_word = [exponents for _, exponents in left]
        right_word = [exponents for _, exponents in right]
        require(left_word == right_word, ("pair prefix equality", i, j))
        oracle_length = min(1000, prefix_length)
        for generators, prefix in (((lower[i], lower[j]), left), ((upper[i], upper[j]), right)):
            require([value for value, _ in prefix[:oracle_length]] ==
                    pointer_values(generators, oracle_length), ("pair pointer oracle", generators))
        pairs.append({"coordinates": [i, j], "identical_exponent_prefix_length": prefix_length,
                      "exponent_word_sha256": hashlib.sha256(
                          json.dumps(left_word, separators=(",", ":")).encode()).hexdigest(),
                      "pointer_oracle_ranks_per_triple": oracle_length})
    left, right = exact_prefix(lower, 6), exact_prefix(upper, 6)
    difference = next(rank for rank, (a, b) in enumerate(zip(left, right), 1) if a[1] != b[1])
    require(difference == 5 and left[4][1] == (0, 0, 1) and right[4][1] == (1, 1, 0),
            ("triple order differs at rank five", left, right))
    return {"lower_prime_triple": lower, "upper_prime_triple": upper,
            "pair_comparisons": pairs, "first_triple_exponent_difference": difference,
            "lower_first_six": left, "upper_first_six": right,
            "scope": "equality of exponent patterns, not of numerical sequence values; "
                     "a finite witness for the general obstruction in Proposition 15"}


def check_base_power_cycle_family():
    examples = []
    for base in range(3, 32):
        def a(last):
            return [0, 1, 0] + [1] * (base - 1) + [last]

        def b(last):
            return [0, 0] + [1] * base + [last]

        cycle = [a(1), b(base), a(0), b(1), a(base), b(0)]
        require(len({tuple(word) for word in cycle}) == 6,
                ("six distinct words in base-power family", base))
        for step, word in enumerate(cycle):
            transformed = [floor_to_base_power(count, base)
                           for count in occurrence_transform(word)]
            require(transformed == cycle[(step + 1) % 6],
                    ("base-power family transition", base, step))
        if base in (3, 5, 7):
            examples.append({"base": base, "minimal_period": 6, "cycle": cycle})
    return {
        "status": "leading reconstruction; user reports base-3 example looks familiar",
        "rule": "count earlier occurrences, then round down to 0 or a power of the base",
        "bases_checked": list(range(3, 32)),
        "general_claim": "The proved six-cycle templates work for every integer base >= 3.",
        "scope": "finite words; the original infinite fractal construction is not recovered",
        "examples": examples,
    }


def check_model(generators, prefix):
    values = [v for v, _ in prefix]
    word = [e for _, e in prefix]
    oracle_n = min(len(prefix), 1000)
    require(values[:oracle_n] == pointer_values(generators, oracle_n),
            ("independent pointer oracle", generators))
    require(all(a < b for a, b in zip(values, values[1:])), "strict ordering")
    ranks = {value: n for n, value in enumerate(values, 1)}
    displacement_checks = 0
    occurrence_checks = 0
    coordinates = []
    for i, base in enumerate(generators):
        face = [v for v, e in prefix if e[i] == 0]
        face_ranks = {e[:i] + e[i + 1:]: rank
                      for rank, (_, e) in enumerate(
                          (entry for entry in prefix if entry[1][i] == 0))}
        exponents = [e[i] for e in word]
        partner = occurrence_transform(exponents)
        expected = [face_ranks[e[:i] + e[i + 1:]] for e in word]
        require(partner == expected, ("occurrences give remaining-face rank", generators, i))
        require(occurrence_transform(partner) == exponents,
                ("exponent stream two-pass return", generators, i))
        occurrence_checks += len(word)
        sample = []
        for n, value in enumerate(values, 1):
            target = base * value
            if target > values[-1]:
                break
            rank = ranks[target]
            face_count = bisect_right(face, target)
            require(rank == n + face_count,
                    ("rank displacement", generators, i, n, rank, face_count))
            displacement_checks += 1
            if n in (1, 10, 100, 1000):
                sample.append({"n": n, "R_i_n": rank, "displacement": rank - n})
        coordinates.append({
            "generator": base,
            "zero_count": len(face),
            "exponent_word_first_terms": exponents[:16],
            "occurrence_partner_first_terms": partner[:16],
            "sample_embeddings": sample,
        })
    zero_counts = []
    for i in range(len(generators)):
        counts = [0]
        for exps in word:
            counts.append(counts[-1] + (exps[i] == 0))
        zero_counts.append(counts)
    coupled_checks = 0
    pair_summaries = []
    for i in range(len(generators)):
        for j in range(i + 1, len(generators)):
            face_count = 0
            samples = []
            zi, zj = zero_counts[i], zero_counts[j]
            for n, exps in enumerate(word, 1):
                face_count += exps[i] == 0 and exps[j] == 0
                qi, qj = n - zi[n], n - zj[n]
                require(zi[n] - zi[qj] == face_count,
                        ("coupled zero masks", generators, i, j, n))
                require(zj[n] - zj[qi] == face_count,
                        ("symmetric coupled zero masks", generators, i, j, n))
                require(qj - zi[qj] == qi - zj[qi],
                        ("count map commutation", generators, i, j, n))
                if len(generators) == 2:
                    require(face_count == 1, ("two-generator correction", n))
                if n in (10, 100, 1000, len(word)):
                    samples.append({"n": n, "remaining_face_count": face_count})
                coupled_checks += 1
            pair_summaries.append({"coordinates": [i, j], "samples": samples})
    return {
        "generators": generators,
        "ranks": len(prefix),
        "pointer_oracle_ranks": oracle_n,
        "trimmed_prefix_lengths": check_trimming(word),
        "commutation_checks": check_embeddings(word),
        "displacement_checks": displacement_checks,
        "occurrence_return_and_face_rank_checks": occurrence_checks,
        "coupled_zero_mask_checks": coupled_checks,
        "face_intersections": pair_summaries,
        "coordinates": coordinates,
    }


def check_prime_bridge():
    """Use an independent sieve as the oracle for first gaps and prime counts."""
    flags = [True] * 257
    flags[0] = flags[1] = False
    for p in range(2, 17):
        if flags[p]:
            for multiple in range(p * p, len(flags), p):
                flags[multiple] = False
    primes = [n for n, flag in enumerate(flags) if flag]
    rows = []
    thresholds = 0
    for k in range(1, 17):
        generators = primes[:k]
        y = generators[-1]
        # The first 2*y values include all smooth values at most 2*y.
        prefix = exact_prefix(generators, 2 * y)
        values = [v for v, _ in prefix]
        present = set(values)
        gap = next(n for n in range(2, 2 * y + 1) if n not in present)
        require(gap == primes[k], ("first missing integer", k, gap))
        first_rank_gap = next(n for n, v in enumerate(values, 1) if v > n)
        require(first_rank_gap == gap, ("first rank gap", k, gap))
        for x in range(y, 2 * y + 1):
            smooth_count = bisect_right(values, x)
            prime_count = sum(flags[:x + 1])
            require(prime_count == k + x - smooth_count,
                    ("prime count bridge", k, x))
            thresholds += 1
        rows.append({"k": k, "y": y, "first_missing": gap})
    return {"families": rows, "integer_threshold_checks": thresholds}


def check_arbitrary_zero_mask(count):
    """An independently chosen mask: zero at every square position."""
    stream = []
    zero_count = 0
    for n in range(1, count + 1):
        if isqrt(n) ** 2 == n:
            zero_count += 1
            stream.append(0)
        else:
            stream.append(1 + stream[n - zero_count - 1])
    trimmed = [e - 1 for e in stream if e > 0]
    require(trimmed == stream[:len(trimmed)], "arbitrary square zero mask")
    return {"zero_positions": "positive squares", "ranks": count,
            "trimmed_prefix_length": len(trimmed), "first_terms": stream[:24]}


def integer_determinant(rows):
    """Fraction-free elimination with exact divisions (Bareiss)."""
    matrix = [list(row) for row in rows]
    size = len(matrix)
    sign, previous = 1, 1
    for column in range(size - 1):
        pivot = next((i for i in range(column, size) if matrix[i][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            matrix[pivot], matrix[column] = matrix[column], matrix[pivot]
            sign *= -1
        value = matrix[column][column]
        for i in range(column + 1, size):
            for j in range(column + 1, size):
                numerator = matrix[i][j] * value - matrix[i][column] * matrix[column][j]
                quotient, remainder = divmod(numerator, previous)
                require(remainder == 0, "exact Bareiss division")
                matrix[i][j] = quotient
            matrix[i][column] = 0
        previous = value
    return sign * matrix[-1][-1]


def solve_unit_system(rows):
    """M x=1 as integer numerators and positive denominator; None if singular."""
    determinant = integer_determinant(rows)
    if not determinant:
        return None
    sign = 1 if determinant > 0 else -1
    numerators = tuple(sign * integer_determinant(
        [[1 if j == column else value for j, value in enumerate(row)] for row in rows])
        for column in range(len(rows)))
    denominator = abs(determinant)
    require(all(sum(a * x for a, x in zip(row, numerators)) == denominator for row in rows),
            "independent substitution of Cramer solution")
    return numerators, denominator


def minimal_excluded_border(word):
    """Minimal complement of an explicitly supplied finite nonnegative downset."""
    points = set(word)
    require(bool(points), "nonempty prefix")
    dimension = len(word[0])
    require(all(len(point) == dimension and all(isinstance(x, int) and x >= 0 for x in point)
                for point in points), "nonnegative exponent vectors of equal dimension")
    require(len(points) == len(word), "distinct prefix vectors")
    candidates = set()
    for point in points:
        for i in range(dimension):
            if point[i]:
                parent = list(point)
                parent[i] -= 1
                require(tuple(parent) in points, ("prefix downset", point, i))
            child = list(point)
            child[i] += 1
            child = tuple(child)
            if child not in points:
                candidates.add(child)
    return sorted(point for point in candidates if all(
        tuple(x - (j == i) for j, x in enumerate(point)) in points
        for i in range(dimension) if point[i]))


def ordered_prefix_integer_weights(word, max_bases=2000000):
    """Exact small-instance oracle for Proposition 17, not a fast unrank kernel.

    Enumerates bases of the ordered corner-cut LP. The cap is a failure to
    finish the search, never evidence of infeasibility. Numerical logs are absent.
    """
    border = minimal_excluded_border(word)
    dimension = len(word[0])
    require(dimension > 0 and word[0] == (0,) * dimension, "prefix begins at zero")
    bounds = tuple(max(point[i] for point in word) for i in range(dimension))
    active = tuple(i for i, bound in enumerate(bounds) if bound)
    if not active:
        return {"weights": [1] * dimension, "cutoff_score": 0,
                "active_coordinates": [], "axis_bounds": list(bounds),
                "border": border, "constraints": [], "vertex_basis": [],
                "vertex_denominator": 1, "bases_examined": 0}
    reduced_word = [tuple(point[i] for i in active) for point in word]
    reduced_border = [tuple(point[i] for i in active) for point in border
                      if all(point[i] == 0 for i in range(dimension) if i not in active)]
    size = len(active)
    constraints = {tuple(int(i == j) for j in range(size)) for i in range(size)}
    constraints.update(tuple(y - x for x, y in zip(left, right))
                       for left, right in zip(reduced_word, reduced_word[1:]))
    constraints.update(tuple(y - x for x, y in zip(reduced_word[-1], point))
                       for point in reduced_border)
    constraints = sorted(constraints)
    for examined, basis in enumerate(combinations(constraints, size), 1):
        if examined > max_bases:
            raise RuntimeError("prefix compression basis-search cap reached")
        solved = solve_unit_system(basis)
        if solved is None:
            continue
        integers, denominator = solved
        if not all(x >= denominator for x in integers):
            continue
        if not all(sum(c * x for c, x in zip(row, integers)) >= denominator for row in constraints):
            continue
        common = gcd(*integers)
        integers = [x // common for x in integers]
        cutoff = sum(x * a for x, a in zip(reduced_word[-1], integers))
        weights = [cutoff + 1] * dimension
        for i, value in zip(active, integers):
            weights[i] = value
        return {"weights": weights, "cutoff_score": cutoff,
                "active_coordinates": list(active), "axis_bounds": list(bounds),
                "border": border, "constraints": constraints, "vertex_basis": basis,
                "vertex_denominator": denominator, "bases_examined": examined}
    raise ValueError("ordered prefix LP infeasible")


def integer_score_prefix(weights, count):
    """Independent additive heap: returns scores/vectors plus the next score."""
    zero = (0,) * len(weights)
    heap = [(0, zero)]
    seen = {zero}
    result = []
    for _ in range(count):
        score, point = heapq.heappop(heap)
        result.append((score, point))
        for i, weight in enumerate(weights):
            child = list(point)
            child[i] += 1
            child = tuple(child)
            if child not in seen:
                seen.add(child)
                heapq.heappush(heap, (score + weight, child))
    return result, heap[0][0]


def check_prefix_compression():
    """Check full prefix order and exclusion, Cramer bounds, and short witnesses."""
    fixtures = ((2,), (2, 3), (101, 103), (2, 3, 5), (2, 3, 101),
                (101, 683, 68963), (101, 683, 68993), (2, 3, 5, 7),
                (2, 3, 101, 103), (4, 9, 125))
    ranks = list(range(1, 21)) + [32, 64, 100]
    cases = []
    for generators in fixtures:
        prefix = exact_prefix(generators, max(ranks))
        require([value for value, _ in prefix] == pointer_values(generators, max(ranks)),
                ("compression input pointer oracle", generators))
        for count in ranks:
            word = [point for _, point in prefix[:count]]
            witness = ordered_prefix_integer_weights(word)
            weights = witness["weights"]
            score = witness["cutoff_score"]
            additive, next_score = integer_score_prefix(weights, count)
            require([point for _, point in additive] == word,
                    ("integer scores preserve entire ordered prefix", generators, count))
            require(all(a < b for (a, _), (b, _) in zip(additive, additive[1:])),
                    ("integer scores have no prefix ties", generators, count))
            require(next_score > score and all(sum(x * a for x, a in zip(point, weights)) > score
                                               for point in witness["border"]),
                    ("all excluded points have larger score", generators, count))
            require(all(a * b // gcd(a, b) > score for a, b in combinations(weights, 2)),
                    ("an axis collision cannot occur inside a strict prefix", generators, count))
            active = witness["active_coordinates"]
            size = len(active)
            bounds = [witness["axis_bounds"][i] for i in active]
            constant = size * factorial(size) ** 2 * 2 ** (size - 1) if size else 0
            require(count - 1 <= score <= constant * count,
                    ("linear score bound", generators, count, score))
            require(prod(bounds) <= factorial(size) * count,
                    ("simplex product bound", generators, count))
            if size:
                for j, i in enumerate(active):
                    bound = factorial(size) * 2 ** (size - 1) * prod(
                        value for k, value in enumerate(bounds) if k != j)
                    require(weights[i] <= bound, ("Cramer coordinate bound", generators, count, i))
                require(all(abs(row[j]) <= 2 * bounds[j]
                            for row in witness["constraints"] for j in range(size)),
                        ("constraint column bound", generators, count))
            require(all(weights[i] == score + 1 for i in range(len(weights)) if i not in active),
                    ("inactive coordinate exclusion", generators, count))
            cases.append({"generators": list(generators), "rank": count, **witness,
                          "next_integer_score": next_score,
                          "linear_score_bound": constant * count,
                          "last_exponent": list(word[-1])})
    # A downset whose supplied order cannot come from positive linear weights:
    # x < y but 2y < 2x. A failed LP must not silently produce a certificate.
    incompatible = [(0, 0), (1, 0), (0, 1), (0, 2), (1, 1), (2, 0)]
    try:
        ordered_prefix_integer_weights(incompatible)
    except ValueError:
        pass
    else:
        raise AssertionError("incompatible prefix order accepted")
    return {"status": "proved existence bound; novelty and mathematical usefulness unresolved",
            "method": "explicit ordered-prefix LP, exact rational vertex enumeration",
            "scope": "small finite consistency checks; input prefix already enumerated",
            "not_claimed": "fast random access, a new general weight-compression method, "
                           "or a certified replacement for logarithm comparisons",
            "cases_checked": len(cases), "infeasible_order_controls": 1,
            "cases": cases}


def check_monomial_prefix_realization():
    """Check the known finite-weight bridge for matrix-defined orders.

    Matrix heaps provide the inputs; the independent scalar heap checks every
    resulting LP certificate. This does not establish novelty of the bridge.
    """
    cases, effective_dimensions = [], set()
    for dimension in range(2, 5):
        units = [tuple(int(i == j) for j in range(dimension))
                 for i in range(dimension)]
        matrices = {
            "lex": units,
            "reverse_variable_lex": list(reversed(units)),
            "graded_lex": [(1,) * dimension] + units[:-1],
            "graded_reverse_lex": [(1,) * dimension]
                                  + [tuple(-v for v in row) for row in units[:0:-1]],
            "elimination": [units[0], (0,) + (1,) * (dimension - 1)]
                           + units[:1:-1],
        }
        for name, matrix in matrices.items():
            require(integer_determinant(matrix) != 0, "matrix order is total")
            require(all(next(x for x in column if x) > 0 for column in zip(*matrix)),
                    "each variable is positive in the matrix lexicographic order")
            zero = (0,) * dimension
            heap, seen, word = [((0,) * dimension, zero)], {zero}, []
            for _ in range(13):
                _, point = heapq.heappop(heap)
                word.append(point)
                for i in range(dimension):
                    child = tuple(x + (j == i) for j, x in enumerate(point))
                    if child not in seen:
                        seen.add(child)
                        key = tuple(sum(a * b for a, b in zip(row, child)) for row in matrix)
                        heapq.heappush(heap, (key, child))
            for count in (1, 2, 3, 5, 8, 13):
                prefix = word[:count]
                witness = ordered_prefix_integer_weights(prefix)
                weights, score = witness["weights"], witness["cutoff_score"]
                scalar, next_score = integer_score_prefix(weights, count)
                require([point for _, point in scalar] == prefix,
                        ("matrix prefix and scalar certificate agree", name, dimension, count))
                require(all(a < b for (a, _), (b, _) in zip(scalar, scalar[1:])),
                        "all included scalar scores are distinct")
                require(next_score > score, "strict exclusion of the entire complement")
                constant = dimension * factorial(dimension) ** 2 * 2 ** (dimension - 1)
                require(count - 1 <= score <= constant * count,
                        "the existing linear cutoff bound covers matrix-order prefixes")
                require(max(weights) <= constant * count + 1,
                        "positive integer certificate has linearly bounded coordinates")
                axis_levels = [a * b // gcd(a, b) for a, b in combinations(weights, 2)]
                require(all(level > score for level in axis_levels),
                        "every pure-axis collision lies beyond the strict cutoff")
                require(sorted(weights)[1] ** 2 > count - 1,
                        "all but the smallest coordinate weight exceed sqrt(n-1)")
                if count > 1:
                    anchor_scale = max(30, dimension)
                    relative_width = Fraction(1, 25 * anchor_scale ** 2 * (count - 1))
                    require((dimension - 1) * (count - 1) * relative_width < anchor_scale,
                            "short-prime-interval errors remain below one scaled score gap")
                effective_dimensions.add(len(witness["active_coordinates"]))
                cases.append({"order": name, "dimension": dimension, "matrix": matrix,
                              "length": count, "prefix": prefix, "weights": weights,
                              "cutoff_score": score, "next_integer_score": next_score,
                              "axis_collision_levels": axis_levels})
    require(effective_dimensions == set(range(5)), "singleton and every active face exercised")
    return {"cases_checked": len(cases), "dimensions": [2, 4],
            "effective_dimensions": sorted(effective_dimensions), "cases": cases,
            "attribution": "Finite-weight realization follows from Robbiano's ordered-group "
                           "classification; the known correspondence is not introduced here.",
            "scope": "Exact finite falsification checks, not a proof or a novelty assessment."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks", type=int, default=20000)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results/research/exponent_interleaving.json")
    args = parser.parse_args()
    if args.ranks < 2:
        parser.error("--ranks must be at least 2")
    prefixes = {}
    models = {}
    for name, generators in MODELS.items():
        prefixes[name] = exact_prefix(generators, args.ranks)
        models[name] = check_model(generators, prefixes[name])
    scaling = []
    for original, squared in (("primes_2_3", "squared_2_3"),
                              ("primes_2_3_5", "squared_2_3_5")):
        require(all(v * v == w and e == f for (v, e), (w, f)
                    in zip(prefixes[original], prefixes[squared])),
                ("common scaling", original, squared))
        scaling.append({"original": original, "squared": squared,
                        "identical_exponent_ranks": args.ranks})
    original = prefixes["primes_2_3_5"]
    deformed = prefixes["deformed_2_3_5"]
    mismatch = next((n for n, ((_, e), (_, f))
                     in enumerate(zip(original, deformed), 1) if e != f), None)
    # Degree then first coordinate: translation invariant, but not globally
    # representable by one positive real weight per coordinate (see the note).
    graded = []
    degree = 0
    while len(graded) < args.ranks:
        graded.extend((a, degree - a) for a in range(degree + 1))
        degree += 1
    graded = graded[:args.ranks]
    graded_results = {
        "ranks": len(graded),
        "trimmed_prefix_lengths": check_trimming(graded),
        "commutation_checks": check_embeddings(graded),
    }
    result = {
        "task": "exact finite exponent interleaving exploration",
        "scope": "Finite checks of identities proved in docs/interleaving_research.md; "
                 "not an audit of fast kernels or evidence of a new prime bound.",
        "command": shlex.join(["python3", str(Path(__file__).resolve().relative_to(ROOT)), *sys.argv[1:]]),
        "git_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "git_status": subprocess.check_output(
            ["git", "status", "--short", "--", "scripts/run_interleaving_research.py",
             "docs/interleaving_research.md"], cwd=ROOT, text=True).splitlines(),
        "git_status_scope": "Only the two hashed source files; unrelated workspace paths are private.",
        "source_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in (Path(__file__).resolve(), ROOT / "docs/interleaving_research.md")},
        "python_version": sys.version,
        "research_decision": {
            "feasibility": "restricted classification and universal eventual word periods proved",
            "novelty": "unresolved; different from specific cited moves is not proof of novelty",
            "primary_direction": "general structure and certification of fixed-prime rank/unrank",
            "active_goal_candidate": "uniform matching corner-cut count in every fixed dimension; "
                                     "ordered-prefix and sharp prime-height candidates retained for "
                                     "proof and priority audit; structural note remains separate",
            "sequence_direction": "general count-map reconstruction, period and transient laws",
            "expansion": "hold specialized cycle catalogues and paper claims; require a new "
                         "structural or quantitative contribution beyond known general methods",
            "prime_application": "finite-prefix prime-size realization recorded separately; no "
                                 "distribution or fixed-prime unranking consequence. The cycle "
                                 "theorem also applies to composite bases",
        },
        "models": models,
        "common_scaling": scaling,
        "deformation_first_word_difference": mismatch,
        "graded_lex_control": graded_results,
        "arbitrary_scalar_zero_mask": check_arbitrary_zero_mask(args.ranks),
        "occurrence_transform": check_occurrence_dynamics(),
        "restricted_occurrence_candidates": check_restricted_occurrence_cycles(),
        "base_power_cycle_family": check_base_power_cycle_family(),
        "restricted_partition_classification": check_restricted_partition_classification(),
        "word_partition_link": check_word_partition_link(),
        "pairwise_order_obstruction": check_pairwise_order_obstruction(),
        "prefix_compression": check_prefix_compression(),
        "monomial_prefix_realization": check_monomial_prefix_realization(),
        "prime_bridge": check_prime_bridge(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "models": len(models),
        "ranks_per_model": args.ranks,
        "displacement_checks": sum(m["displacement_checks"] for m in models.values()),
        "commutation_checks": sum(m["commutation_checks"] for m in models.values())
                              + graded_results["commutation_checks"],
        "coupled_zero_mask_checks": sum(m["coupled_zero_mask_checks"]
                                        for m in models.values()),
        "occurrence_return_and_face_rank_checks": sum(
            m["occurrence_return_and_face_rank_checks"] for m in models.values()),
        "occurrence_equality_patterns_checked": result["occurrence_transform"][
            "equality_patterns_checked"],
        "restricted_count_candidate_periods": [
            row["minimal_period"]
            for row in result["restricted_occurrence_candidates"]["examples"]],
        "base_power_six_cycle_bases_checked": len(
            result["base_power_cycle_family"]["bases_checked"]),
        "restricted_partition_states_checked": result["restricted_partition_classification"][
            "total_reduced_states_checked"],
        "word_partition_link_words_checked": result["word_partition_link"]["base_power_words_checked"],
        "general_count_map_words_checked": result["word_partition_link"]["general_count_map_words_checked"],
        "identical_pair_prefix_lengths": [row["identical_exponent_prefix_length"]
                                          for row in result["pairwise_order_obstruction"]["pair_comparisons"]],
        "triple_obstruction_first_difference": result["pairwise_order_obstruction"][
            "first_triple_exponent_difference"],
        "prime_bridge_thresholds": result["prime_bridge"]["integer_threshold_checks"],
        "prefix_compression_cases_checked": result["prefix_compression"]["cases_checked"],
        "monomial_prefix_realizations_checked": result["monomial_prefix_realization"]["cases_checked"],
        "all_checks_passed": True,
        "output": str(args.output.resolve()),
    }, indent=2))


if __name__ == "__main__":
    main()
