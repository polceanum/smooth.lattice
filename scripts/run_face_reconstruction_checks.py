#!/usr/bin/env python3
"""Exact bounded controls for occurrence streams and face reconstruction.

The uniform arguments are in docs/occurrence_reconstruction.md. These tests
do not establish novelty, optimal constants, or a faster unranking algorithm.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
import hashlib
import heapq
from itertools import combinations, product
import json
from math import comb, factorial, gcd, isqrt, prod
from pathlib import Path

from run_interleaving_research import (exact_prefix, minimal_excluded_border,
                                       occurrence_transform, pointer_values)


ROOT = Path(__file__).resolve().parents[1]


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def weighted_prefix(weights, count):
    """Heap oracle for positive integer weights, retaining exact exponent words."""
    zero = (0,) * len(weights)
    heap, seen, result, scores = [(0, zero)], {zero}, [], []
    while len(result) < count:
        score, point = heapq.heappop(heap)
        result.append(point)
        scores.append(score)
        for i, weight in enumerate(weights):
            child = tuple(x + (j == i) for j, x in enumerate(point))
            if child not in seen:
                seen.add(child)
                heapq.heappush(heap, (score + weight, child))
    return result, scores


def pointer_scores(weights, count):
    """Independent merge/pointer oracle, not the heap's lattice frontier."""
    values, indices = [0], [0] * len(weights)
    while len(values) < count:
        candidates = [w + values[i] for w, i in zip(weights, indices)]
        value = min(candidates)
        values.append(value)
        for i, candidate in enumerate(candidates):
            if candidate == value:
                indices[i] += 1
    return values


def sparse_difference_controls():
    pairs = occurrence_entries = 0
    fixtures = []
    for generators, count in [((2, 3), 256), ((2, 3, 5), 256),
                              ((2, 3, 5, 7), 128),
                              ((2, 3, 5, 7, 11, 13), 64),
                              ((101, 113, 127, 139, 149, 157, 163, 167), 128),
                              ((4, 9, 25), 128)]:
        prefix = exact_prefix(generators, count)
        require([x for x, _ in prefix] == pointer_values(generators, count),
                ("independent integer-product oracle", generators))
        largest = largest_support = 0
        for j, (_, b) in enumerate(prefix):
            for _, a in prefix[:j]:
                volume = prod(abs(x - y) + 1 for x, y in zip(a, b))
                support = sum(x != y for x, y in zip(a, b))
                require(volume <= 2 * j, ("difference box", generators, j + 1, a, b))
                require(support <= (2 * j).bit_length() - 1, "support consequence")
                largest, largest_support = max(largest, volume), max(largest_support, support)
                pairs += 1
        for i in range(len(generators)):
            stream = [point[i] for _, point in prefix]
            partner = occurrence_transform(stream)
            # Rank in the remaining face, independently ordered by integer values.
            face_values = sorted({value // generators[i] ** point[i]
                                  for value, point in prefix})
            face_rank = {value: rank for rank, value in enumerate(face_values)}
            expected = [face_rank[value // generators[i] ** point[i]]
                        for value, point in prefix]
            require(partner == expected, ("occurrence/face-rank identity", generators, i))
            require(occurrence_transform(partner) == stream, "two-pass return")
            occurrence_entries += count
        fixtures.append({"generators": generators, "prefix_length": count,
                         "maximum_difference_box": largest,
                         "maximum_difference_support": largest_support})
    return {"exponent_pairs": pairs, "occurrence_entries": occurrence_entries,
            "fixtures": fixtures}


def wall_box_controls():
    cases = entries = 0
    by_dimension = defaultdict(int)
    for dimension in range(2, 6):
        for z in product(range(-2, 3), repeat=dimension):
            if not any(x > 0 for x in z) or not any(x < 0 for x in z):
                continue
            if next(x for x in z if x) < 0 or gcd(*z) != 1:
                continue
            j = max(i for i, x in enumerate(z) if x < 0)
            t = -z[j]
            base, scale = 31, 31 ** (dimension + 1)
            boost = 1 + sum(-x for i, x in enumerate(z) if i != j and x < 0)
            values = {i: scale * (boost if z[i] > 0 else 1) + base ** k
                      for k, i in enumerate(i for i in range(dimension) if i != j)}
            weights = tuple(t * values[i] if i != j else
                            sum(z[k] * values[k] for k in values)
                            for i in range(dimension))
            require(min(weights) > 0 and dot(z, weights) == 0, "positive wall witness")
            u, v = tuple(max(x, 0) for x in z), tuple(max(-x, 0) for x in z)
            height, volume = dot(u, weights), prod(abs(x) + 1 for x in z)
            strict, equal = 0, []
            for point in product(*(range(abs(x) + 1) for x in z)):
                score = dot(point, weights)
                strict += score < height
                if score == height:
                    equal.append(point)
                entries += 1
            require(set(equal) == {u, v}, ("only primitive endpoints tie", z, weights))
            require(2 * (strict + 1) == volume, ("first-tie rank lower bound", z))
            cases += 1
            by_dimension[dimension] += 1
    planar = []
    for a in range(1, 21):
        for b in range(1, 21):
            if gcd(a, b) != 1:
                continue
            height = a * b
            below = sum(b * x + a * y < height
                        for x in range(a + 1) for y in range(b + 1))
            require(2 * (below + 1) == (a + 1) * (b + 1), "known planar threshold")
            planar.append((a, b, below + 1))
    return {"primitive_wall_boxes": cases, "box_points": entries,
            "cases_by_dimension": dict(by_dimension),
            "known_planar_threshold_cases": len(planar),
            "planar_attribution": "Adams-Watters, OEIS A163823; not a new threshold formula"}


def bounded_face_reconstruction_controls():
    groups, checked, nontrivial_collisions = 0, 0, 0
    dimensions = []
    for dimension, side in ((3, 5), (4, 4), (5, 3), (6, 2)):
        lengths = [n for n in (2, 3, 4, 7, 8, 15, 16, 31) if n.bit_length() < dimension]
        maximum = max(lengths)
        base = 2 * maximum + 3
        scale = base ** (dimension + 1)
        faces = [s for k in range(1, max(n.bit_length() for n in lengths) + 1)
                 for s in combinations(range(dimension), k)]
        buckets = {n: {} for n in lengths}
        for coarse in product(range(1, side + 1), repeat=dimension):
            weights = tuple(scale * x + base ** i for i, x in enumerate(coarse))
            full, scores = weighted_prefix(weights, maximum)
            require(scores == pointer_scores(weights, maximum), "independent weighted pointer")
            require(len(set(scores)) == maximum, "strict checked prefix")
            face_words = {s: weighted_prefix(tuple(weights[i] for i in s), maximum)[0]
                          for s in faces}
            for n in lengths:
                signature = tuple(tuple(face_words[s][:n]) for s in faces if len(s) <= n.bit_length())
                word = tuple(full[:n])
                previous = buckets[n].get(signature)
                if previous is not None:
                    require(previous == word, ("face reconstruction collision", dimension, n, coarse))
                    nontrivial_collisions += 1
                else:
                    buckets[n][signature] = word
                checked += 1
        groups += sum(len(bucket) for bucket in buckets.values())
        dimensions.append({"dimension": dimension, "coarse_weight_vectors": side ** dimension,
                           "prefix_lengths": lengths})
    return {"weight_prefix_cases": checked, "distinct_face_data_classes": groups,
            "repeated_face_data_controls": nontrivial_collisions,
            "families": dimensions,
            "scope": "Exhaustive listed integer grids with explicit strict tie refinement; not all chambers"}


def degree_points(dimension, maximum, prefix=()):
    if dimension == 1:
        for x in range(maximum + 1):
            yield prefix + (x,)
    else:
        for x in range(maximum + 1):
            yield from degree_points(dimension - 1, maximum - x, prefix + (x,))


def proper_face_obstruction_controls():
    witnesses, faces_checked, face_entries = [], 0, 0
    for r in range(1, 6):
        dimension, face_length = 2 * r, 64
        upper = comb(3 * r, r)
        limit = max(face_length, upper)
        base = 8 * limit + 9
        scale = base ** (dimension + 1)
        z = (1,) * r + (-1,) * r
        central = [scale + base ** i for i in range(dimension - 1)]
        central.append(scale + sum(z[i] * base ** i for i in range(dimension - 1)))
        denominator = 4 * limit + 1
        minus, plus = [[denominator * x for x in central] for _ in range(2)]
        minus[-1] -= 1
        plus[-1] += 1
        u, v = (1,) * r + (0,) * r, (0,) * r + (1,) * r
        height = dot(u, central)
        require(height == dot(v, central), "full-support wall")
        require((r + 1) * min(central) > height, "degree truncation is rigorous")
        candidates = list(degree_points(dimension, r))
        require(len(candidates) == upper, "stars-and-bars bound")
        equal = [p for p in candidates if dot(p, central) == height]
        require(set(equal) == {u, v}, "only the target pair ties below this cutoff")
        first_tie = 1 + sum(dot(p, central) < height for p in candidates)
        left, left_scores = weighted_prefix(minus, upper)
        right, right_scores = weighted_prefix(plus, upper)
        require(left_scores == pointer_scores(minus, upper), "left pointer oracle")
        require(right_scores == pointer_scores(plus, upper), "right pointer oracle")
        require(len(set(left_scores)) == upper and len(set(right_scores)) == upper,
                "strict full fixture prefixes")
        difference = next(i for i, (a, b) in enumerate(zip(left, right), 1) if a != b)
        require(difference == first_tie and difference <= upper, "first full difference")
        require({left[difference - 1], right[difference - 1]} == {u, v}, "target pair switches")
        require(2 ** (dimension - 1) <= difference, "wall support lower bound")
        for size in range(1, dimension):
            for s in combinations(range(dimension), size):
                a, sa = weighted_prefix(tuple(minus[i] for i in s), face_length)
                b, sb = weighted_prefix(tuple(plus[i] for i in s), face_length)
                require(a == b, ("all proper face prefixes agree", dimension, s))
                require(len(set(sa)) == face_length and len(set(sb)) == face_length,
                        "strict face fixture prefixes")
                faces_checked += 1
                face_entries += face_length
        witnesses.append({"dimension": dimension, "proper_face_prefix_length": face_length,
                          "central_integer_weights": central,
                          "minus_integer_weights": minus, "plus_integer_weights": plus,
                          "first_full_difference": difference,
                          "upper_bound_on_first_difference": upper,
                          "switched_vectors": [u, v]})
    return {"proper_faces_checked": faces_checked, "proper_face_entries": face_entries,
            "witnesses": witnesses,
            "scope": "Exact rational-weight fixtures. The uniform prime realization uses the separate proof, not these as prime data."}


def prime_obstruction_control():
    minus, plus = (168323, 106937, 130069, 138373), (168323, 106937, 130069, 138389)
    length, faces = 1000, []
    for p in set(minus + plus):
        require(p >= 2 and all(p % k for k in range(2, isqrt(p) + 1)),
                ("independent trial division", p))
    for size in range(1, 4):
        for s in combinations(range(4), size):
            a = exact_prefix(tuple(minus[i] for i in s), length)
            b = exact_prefix(tuple(plus[i] for i in s), length)
            word = [p for _, p in a]
            require(word == [p for _, p in b], ("prime proper face words", s))
            for bases, prefix in ((tuple(minus[i] for i in s), a),
                                  (tuple(plus[i] for i in s), b)):
                require([x for x, _ in prefix] == pointer_values(bases, length),
                        ("independent proper-face product oracle", s))
            faces.append({"coordinates": s, "prefix_length": length,
                          "exponent_word_sha256": hashlib.sha256(
                              json.dumps(word, separators=(",", ":")).encode()).hexdigest()})
    a, b = exact_prefix(minus, 16), exact_prefix(plus, 16)
    first = next(i for i, (x, y) in enumerate(zip(a, b), 1) if x[1] != y[1])
    require(first == 10, "four-prime first difference")
    require(minus[2] * minus[3] < minus[0] * minus[1] < plus[2] * plus[3],
            "opposite sides of the full-support product comparison")
    for bases, prefix in ((minus, a), (plus, b)):
        require([x for x, _ in prefix] == pointer_values(bases, 16), "full prime pointer oracle")
    return {"minus_prime_bases": minus, "plus_prime_bases": plus,
            "first_full_exponent_difference": first, "proper_faces": faces,
            "minus_first_eleven": a[:11], "plus_first_eleven": b[:11],
            "scope": "Exponent patterns agree; the numerical products are different inputs. A finite witness, not the uniform proof."}


def single_stream_recovery_controls():
    """Exact finite controls; limits and transcendence inputs are proved/cited separately."""
    def ceiling(x):
        return -((-x.numerator) // x.denominator)

    def log_bounds(n):
        def unit(x):
            z = (x - 1) / (x + 1)
            total = sum((z ** (2 * k + 1) / (2 * k + 1) for k in range(48)), Fraction(0)) * 2
            tail = 2 * z ** 97 / (97 * (1 - z * z))
            return total, total + tail
        k = n.bit_length() - 1
        a, b = unit(Fraction(n, 1 << k))
        c, d = unit(Fraction(2))
        return a + k * c, b + k * d

    def ratio_bounds(n, base):
        a, b = log_bounds(n)
        c, d = log_bounds(base)
        return a / d, b / c

    def serial(x):
        return [x.numerator, x.denominator]

    def mask_reconstruct(stream):
        out, zeroes = [0], 0
        for n, value in enumerate(stream, 1):
            zeroes += value == 0
            out.append(0 if value == 0 else 1 + out[n - zeroes])
        require(out[1:] == stream, "zero-mask reconstruction")

    # The grid perturbation argument behind the finite rotation bound, including
    # wrapped intervals, right-closed endpoints and nonzero row offsets.
    rotation_cases = 0
    for denominator in (*range(1, 11), 13, 29, 70):
        numerators = [a for a in range(1, 2 * denominator + 1) if gcd(a, denominator) == 1][:3]
        for numerator in numerators:
            for delta in (Fraction(-1, 2), Fraction(0), Fraction(1, 2)):
                alpha = Fraction(numerator, denominator) + delta / denominator ** 2
                for phase in (Fraction(0), Fraction(1, 7), Fraction(5, 7)):
                    for length in sorted({1, denominator, denominator + 1,
                                          3 * denominator + 2, 5 * denominator}):
                        values = [(phase + j * alpha) % 1 for j in range(length)]
                        for lower, upper in product(range(5), repeat=2):
                            if lower == upper:
                                continue
                            left, right = Fraction(lower, 5), Fraction(upper, 5)
                            width = (right - left) % 1
                            count = sum((left < x <= right) if left < right else
                                        (x > left or x <= right) for x in values)
                            require(abs(count - length * width) <=
                                    8 * (length // denominator) + length % denominator,
                                    "finite rotation discrepancy")
                            rotation_cases += 1

    # An independent exact lattice count, allowing rational test weights and
    # counting all tuples at a tied height. These cases test the finite counting
    # inequality only, not occurrence-array identities under rational weights.
    rational_cases = 0
    for weights in ((Fraction(2, 3),), (Fraction(2, 3), Fraction(5, 4)),
                    (Fraction(2, 3), Fraction(5, 4), Fraction(7, 5))):
        m = len(weights)
        for clock in (2, 4, 8):
            for target in (Fraction(0), weights[0], weights[-1], sum(weights)):
                cap = max(1, ceiling(target))
                points = list(product(*(range((clock + cap) // w + 1) for w in weights)))
                heights = [dot(p, weights) for p in points]
                count = lambda x: sum(y <= x for y in heights)
                full = lambda x: sum((x - y) // 1 + 1 for y in heights if y <= x)
                denominator = count(clock + 1)
                require(full(clock + 1) - full(clock) == denominator, "exact clock increment")
                estimate = Fraction(full(clock + target) - full(clock), denominator)
                bound = (Fraction(8, 3) + Fraction(2 ** m * 3, clock) +
                         (2 * cap + 1) * Fraction(count(clock + cap) - count(clock), denominator))
                require(abs(estimate - target) <= bound, "finite recovery on rational lattice")
                half_count = count(Fraction(clock, 2))
                row_count = sum(dot(p, weights[1:]) <= clock
                                for p in product(*(range(clock // w + 1) for w in weights[1:])))
                row_half = sum(dot(p, weights[1:]) <= Fraction(clock, 2)
                               for p in product(*(range(clock // w + 1) for w in weights[1:])))
                require(row_count <= 2 ** (m - 1) * row_half, "halving fibers")
                require(count(clock) >= (clock // (2 * weights[0]) + 1) * row_half,
                        "product subset in simplex")
                require(row_count * clock <= 2 ** m * weights[0] * count(clock),
                        "row-to-volume bound")
                require(half_count >= 1, "origin convention")
                rational_cases += 1

    entries = rank_checks = certificate_checks = 0
    stream_examples = []
    prefixes = {}
    for bases, count in (((2, 3), 20000), ((2, 3, 5), 20000),
                         ((2, 3, 5, 7), 20000), ((4, 9, 25), 20000)):
        prefix = exact_prefix(bases, count)
        prefixes[bases] = prefix
        require([x for x, _ in prefix] == pointer_values(bases, count), "stream product oracle")
        for anchor in range(len(bases)):
            stream = [p[anchor] for _, p in prefix]
            mask_reconstruct(stream)
            entries += count
            positions = defaultdict(list)
            zero_counts = [0]
            for n, value in enumerate(stream, 1):
                positions[value].append(n)
                zero_counts.append(zero_counts[-1] + (value == 0))
            faces = [value for value, point in prefix if point[anchor] == 0]
            for a, row in positions.items():
                for r, n in enumerate(row):
                    require(prefix[n - 1][0] == bases[anchor] ** a * faces[r],
                            "placement array identifies translated face point")
                    rank_checks += 1
                if a + 1 in positions:
                    require(positions[a + 1][0] - row[0] == zero_counts[positions[a + 1][0]],
                            "clock increment equals face count")
            other = min(x for j, x in enumerate(bases) if j != anchor)
            alpha_lo, alpha_hi = ratio_bounds(other, bases[anchor])
            calibration = next((a, q) for q in range(1, 50) for a in
                               range(max(1, (alpha_lo * q).numerator // (alpha_lo * q).denominator),
                                     ceiling(alpha_hi * q) + 1)
                               if gcd(a, q) == 1 and Fraction(a, q) - Fraction(1, q*q) < alpha_lo
                               and alpha_hi < Fraction(a, q) + Fraction(1, q*q))
            # A short calibration is enough to test correctness. It does not
            # certify a useful numerical precision at these small clock heights.
            p, q = calibration
            height_bound = ceiling(alpha_hi)
            witnesses = []
            for clock in (4, 8, 12, 20, 32, 48):
                if clock + 1 not in positions:
                    continue
                denominator = positions[clock + 1][0] - positions[clock][0]
                for r in range(min(8, len(positions[clock]))):
                    lo, hi = ratio_bounds(faces[r], bases[anchor])
                    cap = max(1, ceiling(hi))
                    if clock + cap not in positions:
                        continue
                    before = zero_counts[positions[clock][0]]
                    after = zero_counts[positions[clock + cap][0]]
                    estimate = Fraction(positions[clock][r] - positions[clock][0], denominator)
                    bound = (Fraction(8, q) + Fraction(2 ** (len(bases) - 1) * height_bound * q, clock) +
                             (2 * cap + 1) * Fraction(after - before, denominator))
                    require(estimate - bound <= lo <= hi <= estimate + bound,
                            "certified prime logarithm lies in recovery interval")
                    certificate_checks += 1
                    if r == 1:
                        witnesses.append({"clock": clock, "estimate": serial(estimate),
                                          "target_log_ratio_interval": [serial(lo), serial(hi)]})
            stream_examples.append({"bases": bases, "anchor_coordinate": anchor,
                                    "checked_stream_length": count,
                                    "first_face_height_estimates": witnesses})
    require([p[0] for _, p in prefixes[(2, 3, 5)]] ==
            [p[0] for _, p in prefixes[(4, 9, 25)]], "composite common-scale control")

    # Non-vacuous finite certificates at large clock heights, without pretending
    # to have enumerated the quadratically long stream. The exact two-weight
    # rank-difference formula is checked separately by direct finite lattice sums.
    root_scale = 1 << 120
    root_floor = isqrt(2 * root_scale ** 2)
    root_lo, root_hi = Fraction(root_floor, root_scale), Fraction(root_floor + 1, root_scale)
    require(root_lo ** 2 < 2 < root_hi ** 2, "rigorous square-root interval")
    # Position-array order alone need not determine numerical heights. A
    # nontrivial monotone map commuting with integer translations preserves it.
    original, deformed = [], []
    for a, r in product(range(25), repeat=2):
        k = isqrt(2 * r*r)
        original.append((Fraction(a) + r * root_lo, Fraction(a) + r * root_hi, (a, r)))
        u, v = k + 2 * r*r + k*k, -2 * r*k
        deformed.append((a + u + v * root_hi, a + u + v * root_lo, (a, r)))
    original.sort()
    deformed.sort()
    for entries_to_check in (original, deformed):
        require(all(x[1] < y[0] for x, y in zip(entries_to_check, entries_to_check[1:])),
                "interval-certified strict rectangle order")
    require([x[2] for x in original] == [x[2] for x in deformed],
            "periodic monotone deformation preserves placement order")
    require(4 - 2 * root_hi > 0 and 2 * (4 - 2 * root_lo) < 14 - 8 * root_hi,
            "deformed heights fail additive closure")
    approximants = [(1, 1), (3, 2)]
    for _ in range(11):
        a, b = approximants[-1]
        c, d = approximants[-2]
        approximants.append((2 * a + c, 2 * b + d))
    for a, b in approximants:
        require(gcd(a, b) == 1 and Fraction(a, b) - Fraction(1, b*b) < root_lo and
                root_hi < Fraction(a, b) + Fraction(1, b*b), "certified quadratic calibration")
    large = []
    direct_cases = 0
    for clock in (4, 16, 64):
        for r in range(8):
            limit = isqrt(clock * clock // 2)
            floor_root = lambda j: isqrt(2 * j*j) if j >= 0 else -isqrt(2 * j*j) - 1
            full_target = sum(clock + floor_root(r - b) + 1 for b in range(limit + r + 1))
            full_zero = sum(clock + floor_root(-b) + 1 for b in range(limit + 1))
            difference = r * (clock + 1) + sum(isqrt(2 * j*j) for j in range(1, r + 1))
            require(full_target - full_zero == difference, "independent quadratic rank sums")
            direct_cases += 1
    for clock in (100, 1000, 10000, 1000000, 1000000000):
        a, q = min(approximants, key=lambda x: Fraction(8, x[1]) + Fraction(4 * x[1], clock))
        count = lambda t: isqrt(t*t // 2) + 1
        for r in (1, 2, 5, 10):
            cap = isqrt(2 * r*r) + 1
            difference = r * (clock + 1) + sum(isqrt(2 * j*j) for j in range(1, r + 1))
            estimate = Fraction(difference, count(clock + 1))
            bound = (Fraction(8, q) + Fraction(4 * q, clock) +
                     (2 * cap + 1) * Fraction(count(clock + cap) - count(clock), count(clock + 1)))
            require(estimate - bound <= r * root_lo < r * root_hi <= estimate + bound,
                    "non-vacuous quadratic finite certificate")
            large.append({"clock": clock, "face_rank": r, "calibration": [a, q],
                          "estimate": serial(estimate), "certified_error_bound": serial(bound)})
    require(large[-4]["certified_error_bound"][0] * 1000 < large[-4]["certified_error_bound"][1],
            "large-clock first-height certificate is narrower than one thousandth")

    minus, plus, count = (101, 683, 466483), (101, 683, 466517), 20000
    for p in set(minus + plus):
        require(all(p % k for k in range(2, isqrt(p) + 1)), "hidden-wall primality")
    left, right = exact_prefix(minus, count), exact_prefix(plus, count)
    for bases, word in ((minus, left), (plus, right)):
        require([x for x, _ in word] == pointer_values(bases, count), "hidden-wall pointer oracle")
    stream = [p[0] for _, p in left]
    require(stream == [p[0] for _, p in right], "long identical anchor stream")
    require(next(n for n, (a, b) in enumerate(zip(left, right), 1) if a[1] != b[1]) == 6,
            "full order changes at rank six")
    return {"zero_mask_entries_reconstructed": entries,
            "placement_rank_identities": rank_checks,
            "finite_rotation_cases": rotation_cases,
            "finite_rational_lattice_bound_cases": rational_cases,
            "finite_prime_stream_certificates": certificate_checks,
            "prime_stream_examples": stream_examples,
            "quadratic_direct_rank_checks": direct_cases,
            "periodic_deformation_array_entries": len(original),
            "periodic_deformation_not_additively_closed": True,
            "large_clock_quadratic_certificates": large,
            "hidden_comparison_witness": {
                "minus_prime_bases": minus, "plus_prime_bases": plus,
                "identical_anchor_coordinate": 0, "identical_stream_prefix_length": count,
                "anchor_word_sha256": hashlib.sha256(json.dumps(stream, separators=(",", ":")).encode()).hexdigest(),
                "first_full_word_difference": 6,
                "minus_first_eight": left[:8], "plus_first_eight": right[:8],
                "scope": "The infinite streams differ; no first mismatch beyond the checked prefix is located."},
            "transcendence_scope": "The six-exponentials input is a cited theorem used in the written proof; these finite tests do not verify it or establish novelty.",
            "scope": "Exact finite identities and interval controls. Large-clock quadratic cases use a proved formula, not a generated huge stream. No unranking speedup or uniform prefix-recovery modulus is claimed."}


def finite_recovery_parameters(n, m, h, height_bound, c, mu, cutoff=None):
    """Explicit bounds under the Diophantine promise stated in Theorem 10.

    This computes a sufficient observation length, not a practical estimate.
    It does not verify the global promise for a proposed unknown source.
    """
    h, height_bound, c = map(Fraction, (h, height_bound, c))
    if n < 2 or m < 1 or mu < 1 or not isinstance(mu, int):
        raise ValueError("n>=2, m>=1 and integer mu>=1 are required")
    if not 0 < h <= height_bound or not 0 < c <= 1:
        raise ValueError("positive weight bounds and 0<c<=1 are required")
    ceil = lambda x: -((-x.numerator) // x.denominator)
    k = max(1, ceil(height_bound)) if cutoff is None else cutoff
    if not isinstance(k, int) or k < 1 or k < height_bound:
        raise ValueError("integer cutoff must cover all remaining weights")
    atom_coefficient_bound = max(1, 2 * ceil(Fraction(k, 1) / h))
    eta = c / atom_coefficient_bound ** mu
    delta = c / n ** mu
    epsilon = min(Fraction(1, 4), h / 4, eta / 16, delta / (8 * m * n))
    pair_c = c / (ceil(height_bound) + 1) ** mu
    b0 = ceil(32 / epsilon)
    budget = ceil(Fraction(b0 ** mu, 1) / pair_c) + 1
    clock = ceil(max(Fraction(1), 4 * 2 ** m * height_bound * budget / epsilon,
                     4 * m * (2 * k + 1) * (k + m * height_bound) / epsilon))
    observations = ceil((clock + k + 1 + m * height_bound) ** (m + 1) /
                        (factorial(m + 1) * h ** m))
    face_bound = ceil((k + m * height_bound) ** m / (factorial(m) * h ** m))
    return {"requested_rank": n, "remaining_dimension": m, "h": h,
            "height_bound": height_bound, "c": c, "mu": mu, "cutoff": k,
            "atom_coefficient_bound": atom_coefficient_bound, "eta": eta,
            "delta": delta, "epsilon": epsilon, "pair_c": pair_c,
            "denominator_minimum": b0, "dirichlet_budget": budget,
            "clock": clock, "observation_bound": observations,
            "face_count_bound": face_bound}


def decode_height_snapshot(n, m, numerators, denominator, radius):
    """Decode accurate relative placement ranks; accuracy is an input promise.

    Theorem 10 supplies that promise from an appropriate exact zero-mask prefix.
    This function alone is not a recognizer for arbitrary numerical snapshots.
    """
    radius = Fraction(radius)
    if denominator <= 0 or not numerators or numerators[0] != 0 or radius <= 0:
        raise ValueError("positive denominator/radius and origin are required")
    if any(a >= b for a, b in zip(numerators, numerators[1:])):
        raise ValueError("estimated heights must be strictly increasing")
    atoms, decompositions = [], {}
    for r in range(1, len(numerators)):
        witness = next(((s, t) for s in range(1, r) for t in range(1, r)
                        if abs(numerators[r] - numerators[s] - numerators[t]) <=
                        3 * radius * denominator), None)
        if witness is None:
            atoms.append(r)
        else:
            decompositions[r] = witness
    if len(atoms) != m:
        raise ValueError("snapshot does not meet the promised atom count/accuracy")
    weights = (denominator,) + tuple(numerators[r] for r in atoms)
    if min(weights) <= 0:
        raise ValueError("recovered weights must be positive")
    word, scores = weighted_prefix(weights, n)
    if len(set(scores)) != n:
        raise ValueError("recovered prefix has a tie; separation was not certified")
    return {"canonical_prefix": word, "integer_weights": weights,
            "atom_face_indices": atoms, "decomposition_witnesses": decompositions,
            "conditional_on_accuracy_and_separation_promises": True}


def recover_from_zero_mask(mask, parameters):
    """Read exact mask bits (1 means exponent zero), charging every entry read.

    Source validity and the numerical class are promises, not checked here.
    The routine rejects a prefix that ends before the required placement data.
    """
    a, k = parameters["clock"], parameters["cutoff"]
    history, positions, zeroes, face_count, reached = [0], defaultdict(list), 0, None, False
    for n, flag in enumerate(mask, 1):
        if flag not in (0, 1) or (n == 1 and not flag):
            raise ValueError("an exact binary mask must start with a zero position")
        zeroes += flag
        value = 0 if flag else 1 + history[n - zeroes]
        history.append(value)
        positions[value].append(n)
        if value == k and len(positions[value]) == 1:
            face_count = zeroes
        if value == a + k and len(positions[value]) == 1:
            reached = True
            break
    if not reached or face_count is None or len(positions[a]) < face_count or not positions[a + 1]:
        raise ValueError("insufficient zero-mask prefix for the required clock")
    origin = positions[a][0]
    result = decode_height_snapshot(parameters["requested_rank"], parameters["remaining_dimension"],
                                    [x - origin for x in positions[a][:face_count]],
                                    positions[a + 1][0] - origin, parameters["epsilon"])
    result["mask_entries_read"] = len(history) - 1
    result["observation_bound"] = parameters["observation_bound"]
    return result


def finite_recovery_controls():
    arithmetic_cases = 0
    parameter_examples = []
    for m in range(1, 5):
        for mu in range(m, m + 3):
            for h, upper in ((Fraction(1), Fraction(2)), (Fraction(1, 2), Fraction(3, 2))):
                for c in (Fraction(1, 3), Fraction(1, 125)):
                    for n in (2, 3, 8, 32):
                        p = finite_recovery_parameters(n, m, h, upper, c, mu)
                        e, a, q, k = p["epsilon"], p["clock"], p["dirichlet_budget"], p["cutoff"]
                        b0 = p["denominator_minimum"]
                        require(q * p["pair_c"] > b0 ** mu, "Dirichlet denominator is forced large")
                        require(Fraction(8, b0) <= e / 4, "rotation error budget")
                        require(Fraction(2 ** m * upper * q, a) <= e / 4, "row error budget")
                        require(2 * m * (2 * k + 1) * (k + m * upper) / a <= e / 2,
                                "boundary error budget")
                        require(a >= m * (k + m * upper), "binomial boundary estimate regime")
                        require(6 * e < p["eta"] and m * n * e < p["delta"] / 2,
                                "atom and word comparison margins")
                        volume_bound = (a + k + 1 + m * upper) ** (m + 1) / (factorial(m + 1) * h ** m)
                        require(p["observation_bound"] - 1 < volume_bound <= p["observation_bound"],
                                "exact rounded observation budget")
                        arithmetic_cases += 1
                        if m == 2 and mu == 3 and h == 1 and c == Fraction(1, 125):
                            parameter_examples.append({"requested_rank": n,
                                                       "clock_bit_length": a.bit_length(),
                                                       "observation_bound_bit_length": p["observation_bound"].bit_length()})

    scale = 1 << 160
    roots = {d: (Fraction(isqrt(d * scale ** 2), scale),
                 Fraction(isqrt(d * scale ** 2) + 1, scale)) for d in (2, 3)}
    for d, (lo, hi) in roots.items():
        require(lo * lo < d < hi * hi, "exact quadratic-root intervals")
    norm_checks = 0
    for z in product(range(-5, 6), repeat=3):
        if not any(z):
            continue
        a, b, c = z
        norm = (a*a + 2*b*b - 3*c*c) ** 2 - 8*a*a*b*b
        require(norm != 0, "nonzero algebraic norm")
        norm_checks += 1

    def interval(a, b, c, divisor=1):
        return (a + (b * roots[2][0] + c * roots[3][0]) / divisor,
                a + (b * roots[2][1] + c * roots[3][1]) / divisor)

    def box_prefix(n, divisor):
        candidates = [(interval(a, b, c, divisor), (a, b, c))
                      for a, b, c in product(range(n), repeat=3)]
        candidates.sort()
        require(all(x[0][1] < y[0][0] for x, y in zip(candidates[:n], candidates[1:n + 1])),
                "independent strict full-box oracle")
        return [point for _, point in candidates[:n]]

    decoder_cases = composite_witnesses = 0
    for divisor, cutoff in ((1, 4), (2, 2), (5, 1)):
        h, upper = Fraction(1, divisor), Fraction(2, divisor)
        c = Fraction(1, divisor * (divisor + 4) ** 3)
        bound = int(Fraction(cutoff, 1) / h)
        face = []
        for b, d in product(range(bound + 1), repeat=2):
            lo, hi = interval(0, b, d, divisor)
            if hi <= cutoff:
                face.append((lo, hi, (b, d)))
            else:
                require(lo > cutoff, "face cutoff decided by rigorous intervals")
        face.sort()
        expected_atoms = [i for i, (_, _, point) in enumerate(face) if sum(point) == 1]
        for n in (2, 5, 9, 16):
            parameters = finite_recovery_parameters(n, 2, h, upper, c, 3, cutoff)
            epsilon = parameters["epsilon"]
            denominator = -((-Fraction(64, 1) / epsilon).numerator // (Fraction(64, 1) / epsilon).denominator)
            oracle = box_prefix(n, divisor)
            for pattern in range(4):
                numerators = [0]
                for r, (lo, hi, _) in enumerate(face[1:], 1):
                    middle = (lo + hi) / 2
                    rounded = (middle * denominator).numerator // (middle * denominator).denominator
                    sign = 1 if (pattern == 0 or (r + pattern) % 3) else -1
                    jitter = (epsilon * denominator * Fraction(4, 5)).numerator // (epsilon * denominator * Fraction(4, 5)).denominator
                    value = rounded + sign * jitter
                    require(Fraction(value, denominator) - epsilon <= lo <= hi <=
                            Fraction(value, denominator) + epsilon, "adversarial snapshot accuracy")
                    numerators.append(value)
                result = decode_height_snapshot(n, 2, numerators, denominator, epsilon)
                require(result["atom_face_indices"] == expected_atoms, "all generators, not just first terms")
                require(result["canonical_prefix"] == oracle, "end-to-end snapshot prefix decoding")
                for r, (s, t) in result["decomposition_witnesses"].items():
                    require(tuple(x + y for x, y in zip(face[s][2], face[t][2])) == face[r][2],
                            "independent additive witness")
                    composite_witnesses += 1
                decoder_cases += 1

    # A complete placement-statistic test for (1,sqrt(2)) at the theorem's
    # actual sufficient clock. No enormous raw mask is generated or read.
    quadratic = []
    for n in (2, 5, 16, 32):
        p = finite_recovery_parameters(n, 1, 1, 2, Fraction(1, 3), 1, cutoff=4)
        a, k = p["clock"], p["cutoff"]
        denominator = isqrt((a + 1) ** 2 // 2) + 1
        face_count = isqrt(k*k // 2) + 1
        numerators = [r * (a + 1) + sum(isqrt(2 * j*j) for j in range(1, r + 1))
                      for r in range(face_count)]
        for r, value in enumerate(numerators):
            require(Fraction(value, denominator) - p["epsilon"] <= r * roots[2][0] <=
                    r * roots[2][1] <= Fraction(value, denominator) + p["epsilon"],
                    "theorem-clock estimates contain the actual heights")
        result = decode_height_snapshot(n, 1, numerators, denominator, p["epsilon"])
        candidates = [(a + r * roots[2][0], a + r * roots[2][1], (a, r))
                      for a, r in product(range(n), repeat=2)]
        candidates.sort()
        require(all(x[1] < y[0] for x, y in zip(candidates[:n], candidates[1:n + 1])),
                "quadratic independent full-box order")
        require(result["canonical_prefix"] == [x[2] for x in candidates[:n]],
                "full theorem-clock snapshot result")
        quadratic.append({"requested_rank": n, "clock": a,
                          "observation_bound": p["observation_bound"],
                          "canonical_prefix": result["canonical_prefix"],
                          "integer_weights": result["integer_weights"],
                          "atom_face_indices": result["atom_face_indices"]})

    # Exercise the complete raw-mask path at smaller clocks whose actual
    # height errors are independently certified. These are integration
    # controls, not executions of the conservative clock bound in (FR1).
    raw_controls = []
    for n, a in ((2, 32), (5, 39), (16, 237), (32, 237)):
        p = finite_recovery_parameters(n, 1, 1, 2, Fraction(1, 3), 1, cutoff=2)
        p["clock"] = a
        p["observation_bound"] = ((a + 5) ** 2 + 1) // 2
        denominator = isqrt((a + 1) ** 2 // 2) + 1
        estimate = Fraction(a + 2, denominator)
        require(estimate - p["epsilon"] <= roots[2][0] <= roots[2][1] <=
                estimate + p["epsilon"], "independent reduced-clock height certificate")
        height = a + p["cutoff"]
        # Exact lattice sum, independent of the stream reader and heap.
        mask_length = height + 1 + sum(height - isqrt(2 * r*r)
                                      for r in range(1, isqrt(height*height // 2) + 1))
        require(3 * mask_length ** 2 < scale, "quadratic ordering survives weight rounding")
        source, _ = weighted_prefix((scale, isqrt(2 * scale ** 2)), mask_length)
        result = recover_from_zero_mask((int(point[0] == 0) for point in source), p)
        candidates = sorted((a0 + r * roots[2][0], a0 + r * roots[2][1], (a0, r))
                            for a0, r in product(range(n), repeat=2))
        require(all(x[1] < y[0] for x, y in zip(candidates[:n], candidates[1:n + 1])),
                "independent raw-control target order")
        require(result["canonical_prefix"] == [x[2] for x in candidates[:n]],
                "complete raw-mask reconstruction")
        require(result["mask_entries_read"] == mask_length <= p["observation_bound"],
                "raw input cost agrees with independent full lattice count")
        require(result["integer_weights"] == (denominator, a + 2),
                "raw placements agree with independent quadratic identity")
        raw_controls.append({"requested_rank": n, "clock": a,
                             "mask_entries_read": mask_length,
                             "height_error_certified_directly": True,
                             "uses_conservative_theorem_clock": False})

    # A real finite mask cannot stand in for the larger input budget.
    count = 512
    require(2 * count * 125 * count ** 3 < scale, "rounding is below the algebraic separation")
    word, _ = weighted_prefix((scale, isqrt(2 * scale ** 2), isqrt(3 * scale ** 2)), count)
    mask = [int(point[0] == 0) for point in word]
    p = finite_recovery_parameters(2, 2, 1, 2, Fraction(1, 125), 3)
    try:
        recover_from_zero_mask(mask, p)
    except ValueError as error:
        require("insufficient" in str(error), "short input rejected for its actual reason")
    else:
        raise AssertionError("short raw prefix was silently accepted")
    require(roots[2][0] > 1 and roots[2][1] < 2, "restricted solvability counterexample")
    return {"parameter_budget_cases": arithmetic_cases,
            "parameter_examples": parameter_examples,
            "algebraic_norm_controls": norm_checks,
            "adversarial_multi_generator_snapshot_cases": decoder_cases,
            "independent_composite_witness_checks": composite_witnesses,
            "theorem_clock_quadratic_snapshots": quadratic,
            "complete_raw_mask_integration_controls": raw_controls,
            "short_raw_mask_entries": count,
            "insufficient_raw_prefix_rejected": True,
            "scope": "Uniform proof establishes the observation bound. Quadratic tests execute the full placement-statistic decoder at that clock using an exact formula. Multigenerator tests exercise adversarial accurate snapshots against independent full-box orders. Separate successful raw-mask controls use smaller clocks with directly certified actual height errors, not the conservative uniform clock. No enormous raw mask is read, and no arbitrary-input promise recognizer or faster unranking claim is made."}


def support_rank_threshold(s):
    """Sufficient lower bound on the first tied rank of an s-support wall.

    The box, strict cyclic and fourth-moment bounds are combined. The
    uniform arguments, including constant weights, are in Theorems 11-14.
    """
    if not isinstance(s, int) or s < 2:
        raise ValueError("mixed comparisons require at least two coordinates")
    q, strict_degree = s // 2, (s - 1) // 2
    degree_count = comb(s + q - 1, q)
    cyclic_count = comb(s + strict_degree, strict_degree)
    return max(2 ** (s - 1), 1 + (cyclic_count + s - 1) // s,
               1 + (degree_count + 47) // 48)


def sharp_face_dimension_bound(n):
    """Uniform sufficient face dimension; this is not the exact finite R(n)."""
    if not isinstance(n, int) or n < 1:
        raise ValueError("a positive integer rank is required")
    if n == 1:
        return 0
    return max(s for s in range(2, n.bit_length() + 1) if support_rank_threshold(s) <= n)


def odd_support_obstruction_controls():
    witnesses, checked_faces = [], 0
    for s in (3, 5, 7, 9):
        a, b = s // 2, (s + 1) // 2
        limit = max(64, comb(s + b, b))
        base, face_length = 8 * limit + 9, 64
        scale = base ** (s + 1)
        z = (1,) * a + (-1,) * b
        central = [scale * (b if i < a else a) + base ** i for i in range(s - 1)]
        central.append(scale * a + sum(z[i] * base ** i for i in range(s - 1)))
        u, v = (1,) * a + (0,) * b, (0,) * a + (1,) * b
        height = dot(u, central)
        require(height == dot(v, central) and (b + 1) * min(central) > height,
                "odd-support wall and rigorous degree cutoff")
        candidates = list(degree_points(s, b))
        require(len(candidates) == comb(s + b, b), "odd-support stars and bars")
        require({p for p in candidates if dot(p, central) == height} == {u, v},
                "only the target pair ties at the relevant height")
        first = 1 + sum(dot(p, central) < height for p in candidates)
        denominator = 4 * limit + 1
        left, right = [denominator * w for w in central], [denominator * w for w in central]
        left[-1] -= 1
        right[-1] += 1
        wl, sl = weighted_prefix(left, len(candidates))
        wr, sr = weighted_prefix(right, len(candidates))
        require(sl == pointer_scores(left, len(candidates)) and
                sr == pointer_scores(right, len(candidates)), "independent odd-support pointer")
        require(len(set(sl)) == len(sl) and len(set(sr)) == len(sr), "strict odd-support prefixes")
        difference = next(i for i, (x, y) in enumerate(zip(wl, wr), 1) if x != y)
        require(difference == first and {wl[first - 1], wr[first - 1]} == {u, v},
                "exact first odd-support change")
        require(support_rank_threshold(s) <= first <= len(candidates), "odd-support rank bracket")
        for size in range(1, s):
            for face in combinations(range(s), size):
                x, sx = weighted_prefix(tuple(left[i] for i in face), face_length)
                y, sy = weighted_prefix(tuple(right[i] for i in face), face_length)
                require(x == y and len(set(sx)) == face_length and len(set(sy)) == face_length,
                        "all proper odd-support faces agree strictly")
                checked_faces += 1
        witnesses.append({"dimension": s, "central_integer_weights": central,
                          "minus_integer_weights": left, "plus_integer_weights": right,
                          "proper_face_prefix_length": face_length,
                          "first_full_difference": first,
                          "upper_bound_on_first_difference": len(candidates),
                          "switched_vectors": [u, v]})
    return {"proper_faces_checked": checked_faces, "witnesses": witnesses,
            "scope": "Exact rational fixtures for odd support. Arbitrary face lengths and prime realizations use the separate uniform argument, not an extrapolation from these cases."}


def sharp_interaction_controls(existing_witnesses):
    moment_cases = constant_cases = composition_scores = orbit_cases = 0
    dimensions = []
    for s in range(2, 11):
        q = s // 2
        points = [x for x in degree_points(s, q) if sum(x) == q]
        all_points = list(degree_points(s, q))
        strict_points = set(degree_points(s, (s - 1) // 2))
        require(len(points) == comb(s + q - 1, q), "exact-degree stars and bars")
        if s <= 6:
            fixtures = list(product(range(1, 4), repeat=s))
        else:
            fixtures = [(1,) * s, tuple(range(1, s + 1)),
                        (s ** 3,) + (1,) * (s - 1),
                        (s ** 2, s ** 2) + (1,) * (s - 2),
                        tuple(1 + i % 2 for i in range(s)),
                        tuple(2 ** i for i in range(s))]
            fixtures += [tuple((i - j) ** 2 + 1 for i in range(s)) for j in range(s)]
        u = [Fraction(1)]
        for k in range(1, 5):
            u.append(u[-1] * Fraction(q - k + 1, s + k - 1))
        r = u[1]
        for weights in fixtures:
            total = sum(weights)
            beta = tuple(s * w - total for w in weights)
            values = [dot(x, beta) for x in points]
            composition_scores += len(points)
            a2, a4 = sum(b*b for b in beta), sum(b**4 for b in beta)
            second = (r + u[2]) * a2
            fourth = ((r + 7*u[2] + 12*u[3] + 6*u[4]) * a4 +
                      (3*u[2] + 6*u[3] + 3*u[4]) * a2*a2)
            require(sum(values) == 0, "centroid exact, without sampling")
            require(Fraction(sum(v*v for v in values), len(points)) == second,
                    "second moment versus complete enumeration")
            require(Fraction(sum(v**4 for v in values), len(points)) == fourth,
                    "fourth moment versus complete enumeration")
            require(fourth <= 12 * second*second, "uniform kurtosis bound")
            if a2:
                require(48 * sum(v < 0 for v in values) >= len(points) and
                        48 * sum(v > 0 for v in values) >= len(points),
                        "both strict centroid halfspaces have positive depth")
                moment_cases += 1
            else:
                constant_cases += 1
            below = sum(2 * dot(x, weights) < total for x in all_points)
            require(48 * below >= len(points), "strict simplex count, including equal weights")
            require(s * below >= len(strict_points), "strict cyclic count")
            # Check every distinct cyclic orbit once for selected contrasting
            # fixtures, separately from the moment and halfspace counts.
            if weights in (fixtures[0], fixtures[1], fixtures[-1]):
                unseen = set(strict_points)
                while unseen:
                    x = min(unseen)
                    orbit = {x[j:] + x[:j] for j in range(s)}
                    scores = [dot(y, weights) for y in orbit]
                    require(s * sum(scores) == len(orbit) * sum(x) * total,
                            "cyclic orbit average")
                    require(2 * min(scores) < total, "a strict point in each cyclic orbit")
                    unseen.difference_update(orbit)
                    orbit_cases += 1
        dimensions.append({"dimension": s, "exact_degree": q,
                           "exact_degree_points": len(points), "weight_cases": len(fixtures)})

    envelope_cases = 0
    for s in range(2, 2049):
        q = s // 2
        r = Fraction(q, s)
        t = Fraction(q * (q - 1), s * (s + 1))
        kappa = 10 + 18*r + 9*r*r
        require(Fraction(1, 3) <= r <= Fraction(1, 2) and 0 <= t <= r*r,
                "uniform moment-parameter regime")
        require(r + kappa*t <= 12 * (r + t) ** 2, "exact rational moment envelope")
        require(6*r + 15*r*r - 10 < 0 and 2*r + 6*r*r + 3*r**3 - 1 > 0,
                "two polynomial steps of uniform envelope proof")
        envelope_cases += 1

    for witness in existing_witnesses:
        require(support_rank_threshold(witness["dimension"]) <= witness["first_full_difference"],
                "new threshold versus independently counted old wall fixture")

    brackets = []
    for exponent in (3, 6, 12, 30, 100, 1000):
        n = 10 ** exponent
        upper = sharp_face_dimension_bound(n)
        k = 1
        while comb(3 * (k + 1), k + 1) <= n:
            k += 1
        lower = 2 * k
        require(comb(3*k, k) <= n < comb(3*(k+1), k+1), "exact lower witness budget")
        require(lower <= upper <= n.bit_length(), "rigorous finite face-dimension bracket")
        require(all(support_rank_threshold(s) > n
                    for s in range(upper + 1, n.bit_length() + 1)), "all omitted supports excluded")
        brackets.append({"rank": "10^" + str(exponent), "prime_obstruction_lower_bound": lower,
                         "new_sufficient_face_dimension": upper,
                         "previous_sufficient_face_dimension": n.bit_length()})
    return {"nonconstant_weight_moment_cases": moment_cases,
            "constant_weight_boundary_cases": constant_cases,
            "composition_scores_checked": composition_scores,
            "distinct_cyclic_orbits_checked": orbit_cases,
            "exact_rational_envelope_cases": envelope_cases,
            "old_wall_fixtures_checked": len(existing_witnesses),
            "dimension_families": dimensions, "finite_dimension_brackets": brackets,
            "odd_support_obstructions": odd_support_obstruction_controls(),
            "uniform_result": "R(n)=(log(n)+0.5*log(log(n)))/log(3*sqrt(3)/2)+O(1), for both independent positive weights and actual prime bases.",
            "scope": "Finite exact moment, strict-boundary, orbit and parameter controls. The uniform asymptotic theorem rests on the written proof and Stirling estimates, not a numerical fit. The finite upper bound is sufficient, not the exact minimum R(n). No novelty or unranking speedup is certified."}


def refined_support_cutoff(n):
    """Use the proved small thresholds J(2)=2, J(3)=5, J(4)=10 as well."""
    if n == 1:
        return 0
    exact_small = {2: 2, 3: 5, 4: 10}
    return max(s for s in range(2, n.bit_length() + 1)
               if max(support_rank_threshold(s), exact_small.get(s, 0)) <= n)


def prefix_cone_constraints(word):
    """Standard strict prefix inequalities, with oriented primitive normals.

    Sparse deletion requires a previously checked strict realization. Use
    prefix_interval_certificate to enforce that prerequisite for a new word.
    """
    word = [tuple(x) for x in word]
    border = minimal_excluded_border(word)
    d = len(word[0])
    if word[0] != (0,) * d:
        raise ValueError("the word must start with the zero vector")
    def primitive(z):
        g = 0
        for x in z:
            g = gcd(g, abs(x))
        if not g:
            raise ValueError("a strict comparison cannot have zero normal")
        return tuple(x // g for x in z)
    positivity = {tuple(int(i == j) for j in range(d)) for i in range(d)}
    normals = set(positivity)
    normals.update(primitive(tuple(b - a for a, b in zip(x, y)))
                   for x, y in zip(word, word[1:]))
    normals.update(primitive(tuple(b - a for a, b in zip(word[-1], x))) for x in border)
    cutoff = min(d, refined_support_cutoff(len(word)))
    retained = sorted(z for z in normals if z in positivity or sum(x != 0 for x in z) <= cutoff)
    return {"all_normals": sorted(normals), "retained_normals": retained,
            "removed_normals": sorted(normals.difference(retained)),
            "support_cutoff": cutoff, "minimal_border": border}


def prefix_interval_certificate(word, center, bounds):
    """Research-only interval certificate, separate from the CLI rank auditor.

    An exact integer center is replayed through n+1 points to establish that
    the proposed word has a nonempty strict realization cone. Arbitrary
    unrealizable words cannot use the sparse-constraint deletion theorem.
    """
    word = [tuple(x) for x in word]
    center = tuple(center)
    if not word or len(center) != len(word[0]) or len(bounds) != len(center):
        raise ValueError("word, center and box dimensions must match")
    if any(not isinstance(x, int) or x <= 0 for x in center):
        raise ValueError("an exact positive integer center is required")
    n = len(word)
    replay, scores = weighted_prefix(center, n + 1)
    if replay[:n] != word or len(set(scores[:n])) != n or scores[n - 1] >= scores[n]:
        raise ValueError("the supplied center does not realize the word strictly")
    bounds = [tuple(map(Fraction, pair)) for pair in bounds]
    if any(lo <= 0 or lo > hi for lo, hi in bounds):
        raise ValueError("positive rational interval endpoints are required")
    result = prefix_cone_constraints(word)
    margins = [sum(z[i] * bounds[i][0 if z[i] >= 0 else 1] for i in range(len(center)))
               for z in result["retained_normals"]]
    result["box_preserves_prefix"] = all(x > 0 for x in margins)
    result["minimum_retained_margin"] = min(margins)
    result["center_replayed_exactly"] = True
    return result


def quotient_fan_boundary_controls():
    """Buchberger controls for the socle-binomial boundary example."""
    witnesses, reductions = [], 0
    for n in (2, 3, 4, 6, 7, 14, 15, 30, 31, 62, 70, 126, 127, 254):
        k = (n + 2).bit_length() - 2
        base_length = 2 ** (k + 1) - 2
        padding, d = n - base_length, 2 * k + 1
        zero = (0,) * d
        u = (1,) * k + (0,) * (k + 1)
        v = (0,) * k + (1,) * k + (0,)
        staircase = set()
        for bits in product((0, 1), repeat=k):
            staircase.add(tuple(bits) + (0,) * (k + 1))
            staircase.add((0,) * k + tuple(bits) + (0,))
        staircase.update((0,) * (d - 1) + (j,) for j in range(1, padding + 1))
        require(len(staircase) == n + 1 and u in staircase and v in staircase,
                "quotient length before the one socle relation")
        generators = minimal_excluded_border(sorted(staircase))
        for top in (u, v):
            for i in range(d):
                child = tuple(x + (i == j) for j, x in enumerate(top))
                require(child not in staircase, "both relation monomials are in the socle")
        for lead, tail in ((u, v), (v, u)):
            standard = staircase.difference({lead})
            require(len(standard) == n, "standard staircase length")
            minimal_excluded_border(sorted(standard))  # validates divisor closure
            for g in generators:
                # S(g,lead-tail) is, up to sign, this single monomial.
                lcm = tuple(max(a, b) for a, b in zip(g, lead))
                remainder = tuple(a - b + c for a, b, c in zip(lcm, lead, tail))
                require(remainder not in staircase and
                        any(all(a <= b for a, b in zip(h, remainder)) for h in generators),
                        "every nontrivial Buchberger pair reduces by the monomial ideal")
                reductions += 1
        normal = tuple(a - b for a, b in zip(u, v))
        left, right = [2] * d, [2] * d
        left[2 * k - 1] -= 1
        right[2 * k - 1] += 1
        require(dot(normal, left) > 0 > dot(normal, right), "positive weights cross the sole fan wall")
        require(sum(x != 0 for x in normal) == 2 * k <= 2 * (n.bit_length() - 1),
                "general quotient support bracket")
        witnesses.append({"quotient_length": n, "half_support": k, "padding": padding,
                          "facet_support": 2*k, "normal": normal,
                          "monomial_generator_count": len(generators),
                          "monomial_generators_sha256": hashlib.sha256(
                              json.dumps(generators, separators=(",", ":")).encode()).hexdigest(),
                          "explicit_small_generators": generators if n == 6 else None,
                          "socle_relation_exponents": [u, v]})
    return {"buchberger_reductions_checked": reductions, "witnesses": witnesses,
            "scope": "The exact monomial-plus-one-binomial presentations have two initial ideals. These are a boundary on extending the free-lattice theorem to arbitrary quotients, not a new universal Groebner-basis algorithm."}


def sparse_cone_controls():
    small_thresholds = []
    for weights, u, v, target in (((1, 1), (1, 0), (0, 1), 2),
                                 ((100, 49, 51), (1, 0, 0), (0, 1, 1), 5),
                                 ((100, 103, 101, 102), (1, 1, 0, 0), (0, 0, 1, 1), 10)):
        height = dot(weights, u)
        degree_limit = height // min(weights)
        points = list(degree_points(len(weights), degree_limit))
        require(dot(weights, v) == height, "small-support wall")
        require({x for x in points if dot(weights, x) == height} == {u, v},
                "only the target pair ties at the relevant height")
        rank = 1 + sum(dot(weights, x) < height for x in points)
        require(rank == target, "exact small-support threshold witness")
        small_thresholds.append({"support": len(weights), "minimum_first_tied_rank": target,
                                 "central_integer_weights": weights, "height": height})

    grid_cases = removed = 0
    for d in (3, 4, 5):
        word = [(0,) * d] + [tuple(int(i == j) for i in range(d)) for j in range(3)]
        center = (10, 12, 14) + (100,) * (d - 3)
        result = prefix_interval_certificate(word, center,
                                             [(Fraction(x) - Fraction(1, 4), Fraction(x) + Fraction(1, 4))
                                              for x in center])
        require(result["box_preserves_prefix"] and result["support_cutoff"] == 2,
                "successful interval certificate with a removed dense constraint")
        target = (1, 1, -1) + (0,) * (d - 3)
        a = (2, 0, -1) + (0,) * (d - 3)
        b = (-1, 1, 0) + (0,) * (d - 3)
        require(result["removed_normals"] == [target] and
                a in result["retained_normals"] and b in result["retained_normals"] and
                tuple(x + y for x, y in zip(a, b)) == target, "independent exact conic identity")
        removed += 1
        accepted = rejected = 0
        for weights in product(range(1, 8), repeat=d):
            full = all(dot(z, weights) > 0 for z in result["all_normals"])
            sparse = all(dot(z, weights) > 0 for z in result["retained_normals"])
            require(full == sparse, "full and sparse cones agree on complete listed integer grid")
            accepted += full
            rejected += not full
            grid_cases += 1
        require(accepted > 0 and rejected > 0, "both outcomes occur in every cone grid")
        bad = [(Fraction(x) - Fraction(1, 4), Fraction(x) + Fraction(1, 4)) for x in center]
        bad[2] = (Fraction(14), Fraction(25))
        require(not prefix_interval_certificate(word, center, bad)["box_preserves_prefix"],
                "box crossing a real boundary is rejected")
    invalid = [(0, 0), (1, 0), (0, 1), (1, 1)]
    try:
        prefix_interval_certificate(invalid, (10, 12), [(10, 10), (12, 12)])
    except ValueError:
        pass
    else:
        raise AssertionError("an unrealizable prefix cannot invoke sparse deletion")

    antichain_controls = []
    for primes in ((168323, 106937, 130069, 138373), (168323, 106937, 130069, 138389)):
        prefix = exact_prefix(primes, 11)
        chain = [x for _, x in prefix[:10]]
        previous, candidate = chain[-1], prefix[10][1]
        maximal = [x for x in chain if not any(x != y and all(a <= b for a, b in zip(x, y))
                                             for y in chain)]
        normal = tuple(a - b for a, b in zip(candidate, previous))
        require(previous in maximal and sum(x != 0 for x in normal) == 4 and
                all(min(a, b) == 0 for a, b in zip(previous, candidate)),
                "published antichain reduction retains the coprime full-support comparison")
        antichain_controls.append({"primes": primes, "known_chain_length": 10,
                                   "retained_previous_monomial": previous,
                                   "next_candidate": candidate, "normal": normal})

    padding_controls = []
    for padding in (1, 3, 7):
        original, height = (100, 103, 101, 102), 203
        scale = 4 * (padding + 1)
        center = tuple(scale * w for w in original) + tuple(
            scale * height - 100 * j for j in range(1, padding + 1))
        d = len(center)
        u, v = (1, 1, 0, 0) + (0,) * padding, (0, 0, 1, 1) + (0,) * padding
        points = list(degree_points(d, 2))
        require(3 * min(center) > scale * height, "padded corner-cut degree cutoff")
        below = {x for x in points if dot(x, center) < scale * height}
        require({x for x in points if dot(x, center) == scale * height} == {u, v} and
                len(below) + 1 == 10 + padding, "one isolated point per padding coordinate")
        left, right = list(center), list(center)
        left[3] -= 1
        right[3] += 1
        n = 10 + padding
        a, sa = weighted_prefix(left, n + 1)
        b, sb = weighted_prefix(right, n + 1)
        require(set(a[:n]) == below.union({v}) and set(b[:n]) == below.union({u}) and
                sa[n - 1] < sa[n] and sb[n - 1] < sb[n], "opposite strict corner-cut sets")
        sums_a = tuple(sum(x[j] for x in a[:n]) for j in range(d))
        sums_b = tuple(sum(x[j] for x in b[:n]) for j in range(d))
        require(tuple(x - y for x, y in zip(sums_a, sums_b)) ==
                tuple(x - y for x, y in zip(v, u)), "state-polyhedron edge normal retains support")
        padding_controls.append({"padding_coordinates": padding, "corner_cut_size": n,
                                 "wall_support": 4, "central_integer_weights": center,
                                 "normal": tuple(x - y for x, y in zip(u, v))})

    def log_interval(n):
        def unit(x):
            z = (x - 1) / (x + 1)
            lower = 2 * sum((z ** (2*j + 1) / (2*j + 1) for j in range(48)), Fraction(0))
            return lower, lower + 2 * z ** 97 / (97 * (1 - z*z))
        k = n.bit_length() - 1
        a, b = unit(Fraction(n, 1 << k))
        c, e = unit(Fraction(2))
        return a + k*c, b + k*e

    prime_certificates = []
    for primes in ((2, 3, 5), (101, 103, 107), (2, 3, 5, 7), (167, 173, 179, 181, 191)):
        bounds = [log_interval(p) for p in primes]
        center = tuple(int((lo + hi) * 10 ** 40 / 2) for lo, hi in bounds)
        for n in (4, 8, 16, 32):
            word, _ = weighted_prefix(center, n)
            result = prefix_interval_certificate(word, center, bounds)
            require(result["box_preserves_prefix"], "rigorous prime-log interval certificate")
            exact = exact_prefix(primes, n)
            require(word == [x for _, x in exact] and [v for v, _ in exact] == pointer_values(primes, n),
                    "independent exact prime product validation")
            prime_certificates.append({"primes": primes, "prefix_length": n,
                                       "support_cutoff": result["support_cutoff"],
                                       "all_constraints": len(result["all_normals"]),
                                       "retained_constraints": len(result["retained_normals"]),
                                       "removed_constraints": len(result["removed_normals"]),
                                       "box_preserves_prefix": True})
    return {"small_exact_thresholds": small_thresholds,
            "integer_grid_cone_comparisons": grid_cases,
            "independent_dense_constraint_identities": removed,
            "invalid_realization_rejected": True, "boundary_crossing_boxes_rejected": 3,
            "prime_interval_certificates": prime_certificates,
            "published_antichain_full_support_controls": antichain_controls,
            "prescribed_size_corner_cut_padding": padding_controls,
            "arbitrary_quotient_boundary": quotient_fan_boundary_controls(),
            "scope": "Exact prefix-cone certificates replay a strict rational center and use rigorous interval signs, validated against exact prime products. This is separate from the CLI rank auditor and adds no certified fast-unranking range. Uniform sharp arity and exact inversion rest on the written proof; novelty remains unresolved."}


def primary_toolkit_comparison_controls():
    """Exact specializations of the cited moment and Borel-fixed frameworks."""
    factorial_cases = masses = 0
    families = []
    def falling(x, order):
        return prod(range(x - order + 1, x + 1)) if order <= x else 0
    for s in range(2, 7):
        for q in (s // 2, s // 2 + 1):
            points = [x for x in degree_points(s, q) if sum(x) == q]
            count = comb(s + q - 1, q)
            require(len(points) == count, "complete composition population")
            # One multinomial row, alpha_i=1, in the primary source's
            # Dirichlet-multinomial formula: every count vector is uniform.
            for x in points:
                cell_factorials = prod(factorial(a) for a in x)
                mass = Fraction(factorial(q) * cell_factorials,
                                cell_factorials * prod(range(s, s + q)))
                require(mass == Fraction(1, count), "published distribution specializes to uniform compositions")
                masses += 1
            for orders in degree_points(s, 4):
                total = sum(orders)
                actual = Fraction(sum(prod(falling(xi, ai) for xi, ai in zip(x, orders))
                                      for x in points), count)
                expected = Fraction(falling(q, total) * prod(factorial(a) for a in orders),
                                    prod(range(s, s + total)))
                require(actual == expected, "source equation (3) versus complete mixed-factorial enumeration")
                factorial_cases += 1
            families.append({"dimension": s, "total_degree": q, "composition_count": count})

    borel_controls = []
    for primes in ((168323, 106937, 130069, 138373), (168323, 106937, 130069, 138389)):
        word = [x for _, x in exact_prefix(primes, 10)]
        permutation = tuple(sorted(range(len(primes)), key=lambda i: primes[i], reverse=True))
        staircase = {tuple(x[i] for i in permutation) for x in word}
        transfers = 0
        for x in staircase:
            for i in range(len(primes)):
                if x[i]:
                    for j in range(i + 1, len(primes)):
                        shifted = tuple(a - (k == i) + (k == j) for k, a in enumerate(x))
                        require(shifted in staircase, "Borel-fixed transfer after descending weight permutation")
                        transfers += 1
        borel_controls.append({"primes": primes, "prefix_length": 10,
                               "descending_weight_permutation": permutation,
                               "borel_transfers_checked": transfers,
                               "last_vector": word[-1]})
    require(borel_controls[0]["descending_weight_permutation"] ==
            borel_controls[1]["descending_weight_permutation"] and
            borel_controls[0]["last_vector"] != borel_controls[1]["last_vector"],
            "same Borel coordinate order does not decide the full interleaving")
    return {"uniform_probability_cases": masses,
            "mixed_factorial_moment_cases": factorial_cases,
            "composition_families": families,
            "borel_fixed_counterpart_controls": borel_controls,
            "scope": "Positive primary-source identification: the auxiliary distribution and factorial moments are established Dirichlet-multinomial specializations; both differing corner cuts obey the published Borel-fixed condition in the same coordinate order. These tests do not establish novelty of the sharp support law."}


def prefix_boundary_convention_controls():
    """Check the empty positive-facet cases separately from global maxima."""
    singleton_cases = one_coordinate_cases = 0
    for d in range(1, 7):
        word = [(0,) * d]
        center = tuple(range(2, d + 2))
        result = prefix_interval_certificate(word, center,
                                             [(Fraction(1, 100), Fraction(10**6))] * d)
        positivity = {tuple(int(i == j) for j in range(d)) for i in range(d)}
        require(result["support_cutoff"] == 0 and result["box_preserves_prefix"] and
                set(result["all_normals"]) == set(result["retained_normals"]) == positivity,
                "a singleton prefix has only coordinate boundaries, no positive facet")
        singleton_cases += 1
    for n in range(2, 17):
        word = [(j,) for j in range(n)]
        result = prefix_interval_certificate(word, (7,), [(Fraction(1, 100), Fraction(10**6))])
        require(result["box_preserves_prefix"] and result["all_normals"] == [(1,)] and
                result["retained_normals"] == [(1,)],
                "one ambient coordinate has no positive facet at any listed prefix length")
        one_coordinate_cases += 1
    return {"singleton_prefix_cases": singleton_cases,
            "one_coordinate_prefix_cases": one_coordinate_cases,
            "global_facet_maximum_statements_require_n_at_least": 2,
            "scope": "The n=1 reconstruction parameter is zero, but no positive-weight facet exists. Global facet maxima are stated for n>=2; a one-dimensional ambient cone can still have no positive facet."}


def main():
    report = {"schema_version": 1, "status": "passed",
              "command": "python3 -B scripts/run_face_reconstruction_checks.py",
              "sparse_difference": sparse_difference_controls(),
              "wall_boxes": wall_box_controls(),
              "face_reconstruction": bounded_face_reconstruction_controls(),
              "proper_face_obstructions": proper_face_obstruction_controls(),
              "four_prime_obstruction": prime_obstruction_control(),
              "single_stream_recovery": single_stream_recovery_controls(),
              "finite_stream_reconstruction": finite_recovery_controls(),
              "novelty_established": False, "publication_ready": False,
              "scope": "Finite exact controls support the uniform proofs; no prime-distribution or unranking speedup claim.",
              "source_sha256": {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
                                for p in ("scripts/run_face_reconstruction_checks.py",
                                          "scripts/run_interleaving_research.py",
                                          "docs/occurrence_reconstruction.md")}}
    report["sharp_interaction_hierarchy"] = sharp_interaction_controls(
        report["proper_face_obstructions"]["witnesses"])
    report["sparse_prefix_cones"] = sparse_cone_controls()
    report["primary_toolkit_comparison"] = primary_toolkit_comparison_controls()
    report["prefix_boundary_conventions"] = prefix_boundary_convention_controls()
    path = ROOT / "results/research/face_reconstruction_checks.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"],
                      "exponent_pairs": report["sparse_difference"]["exponent_pairs"],
                      "primitive_wall_boxes": report["wall_boxes"]["primitive_wall_boxes"],
                      "weight_prefix_cases": report["face_reconstruction"]["weight_prefix_cases"],
                      "proper_faces_checked": report["proper_face_obstructions"]["proper_faces_checked"],
                      "single_stream_entries": report["single_stream_recovery"]["zero_mask_entries_reconstructed"],
                      "single_stream_height_certificates": report["single_stream_recovery"]["finite_prime_stream_certificates"],
                      "output": str(path.relative_to(ROOT))}))


if __name__ == "__main__":
    main()
