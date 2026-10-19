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
from math import comb, gcd, isqrt, prod
from pathlib import Path

from run_interleaving_research import exact_prefix, occurrence_transform, pointer_values


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


def main():
    report = {"schema_version": 1, "status": "passed",
              "command": "python3 -B scripts/run_face_reconstruction_checks.py",
              "sparse_difference": sparse_difference_controls(),
              "wall_boxes": wall_box_controls(),
              "face_reconstruction": bounded_face_reconstruction_controls(),
              "proper_face_obstructions": proper_face_obstruction_controls(),
              "four_prime_obstruction": prime_obstruction_control(),
              "single_stream_recovery": single_stream_recovery_controls(),
              "novelty_established": False, "publication_ready": False,
              "scope": "Finite exact controls support the uniform proofs; no prime-distribution or unranking speedup claim.",
              "source_sha256": {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
                                for p in ("scripts/run_face_reconstruction_checks.py",
                                          "scripts/run_interleaving_research.py",
                                          "docs/occurrence_reconstruction.md")}}
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
