#!/usr/bin/env python3
"""Exact, bounded controls for the general initial-ideal research goal.

This checks proof-route obstructions and a conditional rank lemma, not a new
unrestricted bound.
Uses only Python's standard library. Run from the repository root:
    python3 -B scripts/run_structural_complexity_checks.py
"""

from fractions import Fraction
from heapq import heappop, heappush
from hashlib import sha256
from itertools import combinations, product
from math import comb, factorial, gcd, lcm
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def determinant(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    result = Fraction(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            result = -result
        value = a[i][i]
        result *= value
        a[i] = [x / value for x in a[i]]
        for j in range(i + 1, len(a)):
            value = a[j][i]
            a[j] = [x - value * y for x, y in zip(a[j], a[i])]
    return result


def volume_controls():
    results = []
    for d in range(2, 7):
        for r in range(1, 7):
            staircase = [x for x in product(range(r + 1), repeat=d)
                         if sum(x) <= r]
            n = len(staircase)
            require(n == comb(r + d, d), (d, r, "simplex cardinality"))
            moment = [sum(x[i] for x in staircase) for i in range(d)]
            b = Fraction(n * r, d + 1)
            require(moment == [b] * d, (d, r, "simplex moment"))
            axis_moment = n * (n - 1) // 2
            matrix = [[axis_moment * (i == j) - moment[i]
                       for j in range(d)] for i in range(d)]
            actual = determinant(matrix)
            formula = axis_moment ** (d - 1) * (axis_moment - d * b)
            require(actual == formula > 0, (d, r, "volume determinant"))
            results.append({
                "dimension": d, "degree": r, "n": n,
                "axis_moment": axis_moment, "simplex_moment": moment,
                "simplex_volume": str(actual / factorial(d)),
                "volume_divided_by_n_to_2d":
                    str(actual / (factorial(d) * n ** (2 * d))),
            })
    return {
        "scope": "Exact determinant controls for an infinite symbolic construction; the uniform proof is in the research note.",
        "conclusion": "The full moment-polytope volume has order n^(2d) even in the generic point-ideal case. A smaller universal full-volume bound cannot improve the existing proof.",
        "fixtures": results,
    }


def primitive(normal):
    factor = gcd(*normal)
    normal = tuple(x // factor for x in normal)
    return normal if next(x for x in normal if x) > 0 else tuple(-x for x in normal)


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def area2(polygon):
    return abs(sum(p[0] * q[1] - p[1] * q[0]
                   for p, q in zip(polygon, polygon[1:] + polygon[:1])))


def clip(polygon, normal, sign):
    output = []
    for p, q in zip(polygon, polygon[1:] + polygon[:1]):
        fp, fq = sign * dot(p, normal), sign * dot(q, normal)
        if fp >= 0:
            output.append(p)
        if (fp < 0 < fq) or (fq < 0 < fp):
            t = fp / (fp - fq)
            output.append(tuple(x + t * (y - x) for x, y in zip(p, q)))
    unique = []
    for point in output:
        if not unique or point != unique[-1]:
            unique.append(point)
    if len(unique) > 1 and unique[0] == unique[-1]:
        unique.pop()
    return unique if len(unique) >= 3 and area2(unique) else []


def positive_chambers(normals):
    triangle = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
    polygons = [triangle]
    independent_region_count = 1
    for index, normal in enumerate(normals):
        require(min(normal) < 0 < max(normal), "line crosses the open weight simplex")
        intersections = set()
        for previous in normals[:index]:
            point = cross(normal, previous)
            total = sum(point)
            if total:
                point = tuple(Fraction(x, total) for x in point)
                if min(point) > 0:
                    intersections.add(point)
        independent_region_count += 1 + len(intersections)
        new = []
        for polygon in polygons:
            for sign in (-1, 1):
                piece = clip(polygon, normal, sign)
                if piece:
                    new.append(piece)
        polygons = new
    require(len(polygons) == independent_region_count,
            "polygon clipping agrees with exact line-intersection region recurrence")
    weights = []
    for polygon in polygons:
        centroid = tuple(sum(p[i] for p in polygon) / len(polygon) for i in range(3))
        require(min(centroid) > 0, "strictly positive chamber witness")
        require(all(dot(centroid, normal) != 0 for normal in normals), "strict chamber witness")
        denominator = lcm(*(x.denominator for x in centroid))
        weight = tuple(int(x * denominator) for x in centroid)
        factor = gcd(*weight)
        weights.append(tuple(x // factor for x in weight))
    return weights


def matvec(columns, vector):
    result = [0] * len(vector)
    for coefficient, column in zip(vector, columns):
        if coefficient:
            for i, value in enumerate(column):
                result[i] += coefficient * value
    return tuple(result)


def greedy(points, columns, weight, characteristic):
    pivots, selected = {}, []
    for point in sorted(points, key=lambda p: (dot(weight, p), p)):
        vector = ([x % characteristic for x in columns[point]] if characteristic
                  else [Fraction(x) for x in columns[point]])
        for index, row in sorted(pivots.items()):
            coefficient = vector[index]
            vector = [x - coefficient * y for x, y in zip(vector, row)]
            if characteristic:
                vector = [x % characteristic for x in vector]
        index = next((i for i, value in enumerate(vector) if value), None)
        if index is not None:
            inverse = (pow(vector[index], -1, characteristic) if characteristic
                       else 1 / vector[index])
            row = [x * inverse for x in vector]
            if characteristic:
                row = [x % characteristic for x in row]
            pivots[index] = row
            selected.append(point)
    return set(selected)


def socle_fixture(weights, cutoff, block_size):
    domain = {p for p in product(*(range(cutoff // w + 1) for w in weights))
              if dot(weights, p) <= cutoff}
    socle = sorted(p for p in domain if all(
        tuple(x + (i == j) for j, x in enumerate(p)) not in domain for i in range(3)))
    blocks = [socle[i:i + block_size] for i in range(0, len(socle), block_size)]
    blocks = [block for block in blocks if len(block) >= 2]
    require(blocks, "nontrivial relation fixture")
    pivots = {block[-1]: block[:-1] for block in blocks}
    basis = sorted(domain - set(pivots))
    n, size = len(basis), len(domain)
    coordinates = {point: i for i, point in enumerate(basis)}
    zero = (0,) * n

    def residue(point):
        values = [0] * n
        if point in coordinates:
            values[coordinates[point]] = 1
        elif point in pivots:
            for other in pivots[point]:
                values[coordinates[other]] = -1
        return tuple(values)

    matrices = [tuple(residue(tuple(x + (i == j) for j, x in enumerate(point)))
                      for point in basis) for i in range(3)]
    for i, j in combinations(range(3), 2):
        for a in range(n):
            require(matvec(matrices[i], matrices[j][a]) ==
                    matvec(matrices[j], matrices[i][a]), "commuting multiplication matrices")
    boundary = {tuple(x + (i == j) for j, x in enumerate(point))
                for point in domain for i in range(3)} - domain
    columns = {}
    for point in sorted(domain | boundary):
        vector = residue((0, 0, 0))
        for matrix, power in zip(matrices, point):
            for _ in range(power):
                vector = matvec(matrix, vector)
        require(vector == residue(point), "cyclic monomial evaluation including the border")
        columns[point] = vector
    normals = sorted({primitive(tuple(x - y for x, y in zip(p, q)))
                      for block in blocks for p, q in combinations(block, 2)})
    chamber_weights = positive_chambers(normals)
    witnesses = {}
    for weight in chamber_weights:
        removed = tuple(sorted(max(block, key=lambda p: dot(weight, p)) for block in blocks))
        expected = domain - set(removed)
        for characteristic in (0, 2, 3, 29):
            actual = greedy(domain | boundary, columns, weight, characteristic)
            require(actual == expected and len(actual) == n,
                    (weights, cutoff, block_size, weight, characteristic, "independent greedy basis"))
        witnesses.setdefault(removed, weight)
    h = 1
    while h ** 3 < size:
        h += 1
    require(len(socle) * h <= h ** 3 + 3 * size, "socle chain-packing bound")
    require(size <= 2 * n, "disjoint nontrivial relations preserve at least half the length")
    require(len(witnesses) <= 1 + (3 * len(socle)) * (3 * len(socle) + 1) // 2,
            "three-dimensional convex-hull edge bound")
    if len(blocks) == 1:
        require(len(witnesses) <= len(blocks[0]), "one relation deletes one exposed monomial")
    return {
        "cut_weights": weights, "cutoff": cutoff, "block_size": block_size,
        "monomial_order_ideal": sorted(domain), "socle": socle,
        "relation_supports_all_coefficients_one": blocks,
        "original_length": size, "quotient_length": n,
        "quotient_basis_for_multiplication_matrices": basis,
        "multiplication_matrices_commute": True, "cyclicity_checked": True,
        "all_pair_hyperplanes": normals,
        "exact_arrangement_chambers": len(chamber_weights),
        "distinct_initial_staircases": len(witnesses),
        "fields_checked": ["Q", "F_2", "F_3", "F_29"],
        "greedy_checks": 4 * len(chamber_weights),
        "witnesses": [{"deleted_monomials": removed, "positive_integer_weights": weight}
                      for removed, weight in sorted(witnesses.items())],
    }


def leading_row_pivots(rows, points, weight, characteristic):
    """Initial monomials of a finite relation space, via descending row reduction."""
    order = sorted(range(len(points)), key=lambda i: (dot(weight, points[i]), points[i]),
                   reverse=True)
    matrix = [[row[i] % characteristic if characteristic else Fraction(row[i])
               for i in order] for row in rows]
    rank, pivots = 0, []
    for column in range(len(points)):
        pivot = next((r for r in range(rank, len(matrix)) if matrix[r][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = (pow(matrix[rank][column], -1, characteristic) if characteristic
                   else 1 / matrix[rank][column])
        matrix[rank] = [x * inverse for x in matrix[rank]]
        if characteristic:
            matrix[rank] = [x % characteristic for x in matrix[rank]]
        for r in range(rank + 1, len(matrix)):
            coefficient = matrix[r][column]
            matrix[r] = [x - coefficient * y for x, y in zip(matrix[r], matrix[rank])]
            if characteristic:
                matrix[r] = [x % characteristic for x in matrix[r]]
        pivots.append(points[order[column]])
        rank += 1
        if rank == len(matrix):
            break
    require(rank == len(rows), "relation space has the claimed dimension")
    return set(pivots)


def overlapping_socle_fixture(weights, cutoff, socle_rank):
    domain = {p for p in product(*(range(cutoff // w + 1) for w in weights))
              if dot(weights, p) <= cutoff}
    socle = sorted(p for p in domain if all(
        tuple(x + (i == j) for j, x in enumerate(p)) not in domain for i in range(3)))
    require(1 <= socle_rank < len(socle), "proper nonzero quotient of the socle")
    free = socle[:socle_rank]
    vectors = {p: tuple(int(i == j) for i in range(socle_rank))
               for j, p in enumerate(free)}
    for j, point in enumerate(socle[socle_rank:], start=2):
        vectors[point] = tuple(j ** i for i in range(socle_rank))
    rows = []
    for point in socle[socle_rank:]:
        row = [0] * len(socle)
        row[socle.index(point)] = 1
        for i, coefficient in enumerate(vectors[point]):
            row[i] = -coefficient
        rows.append(row)
    fixed = domain - set(socle)
    basis = sorted(fixed | set(free))
    n = len(basis)
    indices = {p: i for i, p in enumerate(basis)}

    def residue(point):
        vector = [0] * n
        if point in indices:
            vector[indices[point]] = 1
        elif point in vectors:
            for i, coefficient in enumerate(vectors[point]):
                vector[indices[free[i]]] = coefficient
        return tuple(vector)

    matrices = [tuple(residue(tuple(x + (i == j) for j, x in enumerate(p)))
                      for p in basis) for i in range(3)]
    for i, j in combinations(range(3), 2):
        for a in range(n):
            require(matvec(matrices[i], matrices[j][a]) ==
                    matvec(matrices[j], matrices[i][a]), "overlapping relations commute")
    border = {tuple(x + (i == j) for j, x in enumerate(p))
              for p in domain for i in range(3)} - domain
    columns = {}
    for point in sorted(domain | border):
        vector = residue((0, 0, 0))
        for matrix, power in zip(matrices, point):
            for _ in range(power):
                vector = matvec(matrix, vector)
        require(vector == residue(point), "overlapping quotient is cyclic on D and its border")
        columns[point] = vector
    budget = [sum(p[i] for p in socle) for i in range(3)]
    require(all(x + len(socle) <= len(domain) for x in budget), "antichain chain budget")
    require(all(x <= n for x in budget), "all relative state moments fit in the length-n cube")
    normals = sorted({primitive(tuple(x - y for x, y in zip(p, q)))
                      for p, q in combinations(socle, 2)})
    witnesses = positive_chambers(normals)
    fields = []
    full_checks = 0
    for characteristic in (0, 2, 3, 29):
        states = {}
        for weight in witnesses:
            excluded = leading_row_pivots(rows, socle, weight, characteristic)
            selected = greedy(set(socle), vectors, weight, characteristic)
            require(selected == set(socle) - excluded, "relation reduction agrees with residue greediness")
            require(len(selected) == socle_rank, "socle rank")
            moment = tuple(sum(p[i] for p in selected) for i in range(3))
            require(all(0 <= moment[i] <= budget[i] <= n for i in range(3)), "relative state moment budget")
            states.setdefault(tuple(sorted(selected)), weight)
        for selected, weight in states.items():
            actual = greedy(domain | border, columns, weight, characteristic)
            require(actual == fixed | set(selected), "full cyclic-matrix greedy basis")
            full_checks += 1
        fields.append({
            "characteristic": characteristic, "distinct_initial_staircases": len(states),
            "all_chambers_checked": len(witnesses),
            "witnesses": [{"selected_socle_monomials": selected, "positive_integer_weights": weight}
                          for selected, weight in sorted(states.items())],
        })
    return {
        "cut_weights": weights, "cutoff": cutoff, "monomial_order_ideal": sorted(domain),
        "socle": socle, "socle_quotient_rank": socle_rank, "quotient_length": n,
        "relation_matrix_columns_indexed_by_socle": rows,
        "fixed_standard_monomials": sorted(fixed), "socle_coordinate_budgets": budget,
        "multiplication_matrices_commute": True, "cyclicity_checked": True,
        "exact_arrangement_chambers": len(witnesses), "fields": fields,
        "row_relation_and_socle_greedy_checks": 4 * len(witnesses),
        "full_cyclic_matrix_greedy_checks": full_checks,
    }


def integer_prefix(weights, length):
    zero = (0,) * len(weights)
    heap, seen, prefix = [(0, zero)], {zero}, []
    while len(prefix) < length:
        score, point = heappop(heap)
        prefix.append((score, point))
        for i, coefficient in enumerate(weights):
            other = tuple(x + (i == j) for j, x in enumerate(point))
            if other not in seen:
                seen.add(other)
                heappush(heap, (score + coefficient, other))
    return prefix


def rectangular_spanning_fixture(points, characteristic, label):
    n, d = len(points), len(points[0])

    def column(exponent):
        return tuple(product_value(point, exponent) for point in points)

    def product_value(point, exponent):
        value = 1
        for coordinate, power in zip(point, exponent):
            value = value * pow(coordinate, power, characteristic) % characteristic
        return value

    deficient, box_count = [], 0
    for sides in product(range(1, n + 1), repeat=d):
        box = set(product(*(range(x) for x in sides)))
        columns = {p: column(p) for p in box}
        rank = len(greedy(box, columns, (1,) * d, characteristic))
        box_count += 1
        if rank < n:
            deficient.append((len(box), sides, rank))
        if max(sides) == n:
            require(rank == n, "axis cyclicity bounds the entire infinite box check")
    largest = max(deficient) if deficient else (0, (), 0)
    kappa = largest[0] // n + 1
    length = d ** d * kappa * n
    tests = []
    rough_weights = ([(1, 1), (1, 3), (3, 1)] if d == 2
                     else [(1, 1, 1), (1, 2, 4), (4, 2, 1)])
    for rough in rough_weights:
        base = length + 1
        scale = d * length * base ** (d - 1) + 1
        weight = tuple(scale * x + base ** i for i, x in enumerate(rough))
        word = integer_prefix(weight, length + 1)
        require(all(a[0] < b[0] for a, b in zip(word, word[1:])), "strict finite word")
        threshold = Fraction(word[-2][0] + word[-1][0], 2)
        sides = tuple(int(threshold // (d * x)) + 1 for x in weight)
        box = set(product(*(range(x) for x in sides)))
        require(max(dot(weight, p) for p in box) < threshold, "inscribed box is strictly inside the cut")
        require(len(box) * d ** d >= length, "geometric box-to-prefix cardinality bound")
        require(len(box) >= kappa * n, "box satisfies the spanning threshold")
        require(len(greedy(box, {p: column(p) for p in box}, weight, characteristic)) == n,
                "inscribed box spans the quotient")
        prefix = {p for _, p in word[:-1]}
        selected = greedy(prefix, {p: column(p) for p in prefix}, weight, characteristic)
        require(len(selected) == n, "length-M word determines a full standard basis")
        completion = max(i + 1 for i, (_, p) in enumerate(word[:-1]) if p in selected)
        tests.append({"positive_integer_weights": weight, "prefix_length": length,
                      "inscribed_box_sides": sides, "box_volume": len(box),
                      "rank_completion_index": completion,
                      "selected_standard_monomials": sorted(selected)})
    return {
        "label": label, "field_characteristic": characteristic, "points": points,
        "dimension": d, "n": n, "rectangular_spanning_parameter": kappa,
        "finite_boxes_checked": box_count,
        "largest_deficient_box": {"volume": largest[0], "sides": largest[1], "rank": largest[2]},
        "infinite_box_reduction": "Every box with an axis side at least n contains a verified n-element axis basis; all remaining boxes were enumerated.",
        "prefix_checks": tests,
    }


def cartesian_rank_delay_controls():
    """Point ideals with one initial staircase but no short geometric cover."""
    results = []
    d, characteristic = 3, 101
    for axis_length in (3, 6, 9):
        zero = (0,) * d
        staircase = {zero} | {tuple(k * (i == j) for j in range(d))
                              for i in range(d) for k in range(1, axis_length + 1)}
        shift = axis_length + 1
        points = [tuple(x + shift for x in p) for p in sorted(staircase)]
        require(all(0 < x < characteristic for p in points for x in p),
                "translated point coordinates are nonzero")
        for rough in ((1, 1, 1), (1, 2, 3), (3, 2, 1)):
            exponent_bound = d * axis_length * (max(rough) + 1)
            base = exponent_bound + 1
            scale = d * exponent_bound * base ** (d - 1) + 1
            weight = tuple(scale * x + base ** i for i, x in enumerate(rough))
            last = max(staircase, key=lambda p: dot(weight, p))
            threshold = dot(weight, last)
            heap, seen, word = [(0, zero)], {zero}, []
            while heap and heap[0][0] <= threshold:
                score, point = heappop(heap)
                word.append((score, point))
                for i, coefficient in enumerate(weight):
                    other = tuple(x + (i == j) for j, x in enumerate(point))
                    if other not in seen:
                        seen.add(other)
                        heappush(heap, (score + coefficient, other))
            require(all(a[0] < b[0] for a, b in zip(word, word[1:])), "strict Cartesian cover")
            prefix = {p for _, p in word}
            require(staircase <= prefix and word[-1][1] == last, "entire axis staircase is covered")
            columns = {}
            for a in prefix:
                values = []
                for point in points:
                    value = 1
                    for coordinate, power in zip(point, a):
                        value = value * pow(coordinate, power, characteristic) % characteristic
                    values.append(value)
                columns[a] = tuple(values)
            selected = greedy(prefix, columns, weight, characteristic)
            require(selected == staircase, "unique staircase recovered from point evaluation")
            require(len(greedy(prefix - {last}, columns, weight, characteristic)) == len(points) - 1,
                    "rank completes only when the last axis monomial appears")
            lower_side = (axis_length - 1) // d + 1
            require(len(prefix) >= lower_side ** d, "uniform geometric cover lower bound")
            results.append({
                "axis_length": axis_length, "n": len(points), "field_characteristic": characteristic,
                "translated_points": points, "unique_standard_staircase": sorted(staircase),
                "positive_integer_weights": weight, "rank_completion_index": len(prefix),
                "uniform_cover_lower_bound": lower_side ** d,
                "nonzero_point_coordinates": True,
            })
    return results


def main():
    fixtures = []
    for weights, cutoff in [((1, 1, 1), 2), ((1, 1, 1), 3),
                            ((2, 3, 5), 8), ((2, 3, 5), 11),
                            ((2, 3, 5), 14), ((1, 2, 4), 9)]:
        for block_size in (2, 3, 100):
            fixtures.append(socle_fixture(weights, cutoff, block_size))
    overlapping = [overlapping_socle_fixture(weights, cutoff, rank)
                   for weights, cutoff in [((1, 1, 1), 3), ((2, 3, 5), 11), ((2, 3, 5), 14)]
                   for rank in (2, 3)]
    rectangular = [
        rectangular_spanning_fixture(
            [(0, 0), (1, 1), (2, 2), (3, 4), (5, 7), (11, 13), (17, 10)],
            29, "published nongeneric seven-point control"),
        rectangular_spanning_fixture([(t, t ** 5, t ** 25) for t in range(1, 6)],
                                     1009, "five-point exponential evaluation configuration"),
        rectangular_spanning_fixture([(t, t, t) for t in range(1, 8)],
                                     29, "diagonal configuration: spanning parameter need not be one"),
    ]
    cartesian = cartesian_rank_delay_controls()
    sources = ["scripts/run_structural_complexity_checks.py", "docs/interleaving_research.md"]
    report = {
        "schema_version": 1, "command": "python3 -B scripts/run_structural_complexity_checks.py",
        "scope": "Bounded exact checks of route obstructions and a conditional counting lemma; no unrestricted asymptotic advance or novelty assertion.",
        "source_sha256": {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in sources},
        "volume_obstruction": volume_controls(),
        "disjoint_socle_relation_control": {
            "uniform_conclusion": "For this whole three-variable construction, g(I)=O(n^(4/3)); it cannot beat the existing generic n^2(log n)^2 order.",
            "proof_location": "docs/interleaving_research.md, First feasibility pass: two rigorous route obstructions",
            "arbitrary_overlapping_relation_spaces_covered": False,
            "fixtures": fixtures,
        },
        "overlapping_socle_relations": {
            "uniform_conclusion": "Arbitrary relation spaces confined to the socle give O_d(n^(d(d-1)/(d+1))) initial staircases, by an antichain coordinate budget and the established lattice-polytope bound.",
            "scope_exclusion": "Relations with surviving nonconstant multiples are not covered.",
            "fixtures": overlapping,
        },
        "rectangular_spanning_control": {
            "conditional_conclusion": "If every monomial box of volume at least kappa*n spans, every standard basis is determined by a geometric word of length d^d*kappa*n. The existing ordered-prefix theorem then gives O_d((kappa*n)^(d-1)(log(kappa*n+1))^(d-1)).",
            "scope_exclusion": "The spanning parameter is not universally bounded; the diagonal family has growing parameter and only d initial ideals.",
            "fixtures": rectangular,
        },
        "cartesian_rank_delay_control": {
            "uniform_conclusion": "Translated Cartesian axis staircases give radical point ideals with g(I)=1, no zero monomial evaluations, and minimum geometric rank-completion length Omega_d(n^d). A linear rank-delay bound cannot be imposed on all point ideals.",
            "fixtures": cartesian,
        },
        "summary": {"volume_fixtures": 30, "socle_fixtures": len(fixtures),
                    "exact_chambers": sum(f["exact_arrangement_chambers"] for f in fixtures),
                    "independent_greedy_checks": sum(f["greedy_checks"] for f in fixtures),
                    "overlapping_socle_fixtures": len(overlapping),
                    "overlapping_row_and_greedy_checks": sum(f["row_relation_and_socle_greedy_checks"] for f in overlapping),
                    "overlapping_full_matrix_checks": sum(f["full_cyclic_matrix_greedy_checks"] for f in overlapping),
                    "rectangular_boxes_checked": sum(f["finite_boxes_checked"] for f in rectangular),
                    "rectangular_prefix_checks": sum(len(f["prefix_checks"]) for f in rectangular),
                    "cartesian_rank_delay_fixtures": len(cartesian)},
        "goal_qualifying_advance": False,
    }
    target = ROOT / "results/research/structural_complexity_checks.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
