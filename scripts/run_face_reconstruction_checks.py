#!/usr/bin/env python3
"""Exact bounded controls for occurrence streams and face reconstruction.

The uniform arguments are in docs/occurrence_reconstruction.md. These tests
do not establish novelty, optimal constants, or a faster unranking algorithm.
"""
from __future__ import annotations

from collections import defaultdict
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


def main():
    report = {"schema_version": 1, "status": "passed",
              "command": "python3 -B scripts/run_face_reconstruction_checks.py",
              "sparse_difference": sparse_difference_controls(),
              "wall_boxes": wall_box_controls(),
              "face_reconstruction": bounded_face_reconstruction_controls(),
              "proper_face_obstructions": proper_face_obstruction_controls(),
              "four_prime_obstruction": prime_obstruction_control(),
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
                      "output": str(path.relative_to(ROOT))}))


if __name__ == "__main__":
    main()
