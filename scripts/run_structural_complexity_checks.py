#!/usr/bin/env python3
"""Exact, bounded controls for the general initial-ideal research goal.

This checks proof-route obstructions, a conditional rank lemma and the
binomial/parallel-class reduction, not a new unrestricted bound.
Uses only Python's standard library. Run from the repository root:
    python3 -B scripts/run_structural_complexity_checks.py
"""

from fractions import Fraction
from heapq import heappop, heappush
from hashlib import sha256
from itertools import combinations, permutations, product
from math import comb, factorial, gcd, isqrt, lcm
from pathlib import Path
from random import Random
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


def hyperbolic_cross(n, d):
    """All possible members of an n-element monomial staircase."""
    return [a for a in product(range(n), repeat=d)
            if _box_volume(a) <= n]


def _box_volume(a):
    value = 1
    for coordinate in a:
        value *= coordinate + 1
    return value


def binomial_controls():
    """Compare class minima with cyclic matrices and independent elimination."""
    fixtures = []
    for n, powers in ((4, (1, 2, 3)), (5, (1, 2, 3)),
                      (6, (1, 2, 4)), (7, (1, 2, 4)), (9, (1, 3, 5))):
        ground = hyperbolic_cross(n, 3)
        for cyclic in (False, True):
            classes = {j: [] for j in range(n)}
            loops = []
            for a in ground:
                degree = dot(a, powers)
                if not cyclic and degree >= n:
                    loops.append(a)
                else:
                    classes[degree % n].append(a)
            require(all(classes.values()), "all n residue classes occur in V_n")
            normals = sorted({primitive(tuple(x - y for x, y in zip(a, b)))
                              for group in classes.values() for a, b in combinations(group, 2)
                              if min(x - y for x, y in zip(a, b)) < 0 <
                              max(x - y for x, y in zip(a, b))})
            weights = positive_chambers(normals)
            expected = []
            witnesses = {}
            for weight in weights:
                selected = frozenset(min(group, key=lambda a: (dot(weight, a), a))
                                     for group in classes.values())
                require(len(selected) == n, "one representative per residue class")
                require(all(tuple(x - (i == j) for j, x in enumerate(a)) in selected
                            for a in selected for i in range(3) if a[i]),
                        "class minima form a monomial staircase")
                moment = tuple(sum(a[i] for a in selected) for i in range(3))
                if selected not in witnesses:
                    witnesses[selected] = {"positive_integer_weights": weight,
                                           "standard_monomials": sorted(selected),
                                           "moment_sum": moment}
                expected.append(selected)
            require(len({tuple(w["moment_sum"]) for w in witnesses.values()}) == len(witnesses),
                    "different coherent staircases have different moment vertices")
            checks = 0
            field_results = []
            for characteristic in (0, 2, 3, 29):
                scalars = ((1, 1, 1) if characteristic == 2 else
                           (1, 2, 1) if characteristic == 3 else (1, 2, 3))
                matrices = []
                for power, scalar in zip(powers, scalars):
                    matrices.append([tuple(scalar * (j == ((k + power) % n))
                                           * (cyclic or k + power < n)
                                           for j in range(n)) for k in range(n)])
                for i, j in combinations(range(3), 2):
                    for k in range(n):
                        require(matvec(matrices[i], matrices[j][k]) ==
                                matvec(matrices[j], matrices[i][k]),
                                "actual binomial quotient multiplication commutes")
                unit = (1,) + (0,) * (n - 1)
                axis = []
                vector = unit
                for _ in range(n):
                    axis.append(vector)
                    vector = matvec(matrices[0], vector)
                require(axis == [tuple(i == j for i in range(n)) for j in range(n)],
                        "x powers make the representation cyclic of dimension n")
                columns = {}
                for a in sorted(ground, key=lambda a: (sum(a), a)):
                    if not any(a):
                        columns[a] = unit
                    else:
                        i = next(i for i, value in enumerate(a) if value)
                        parent = tuple(x - (i == j) for j, x in enumerate(a))
                        columns[a] = matvec(matrices[i], columns[parent])
                    degree, scalar = dot(a, powers), 1
                    for coefficient, exponent in zip(scalars, a):
                        scalar *= coefficient ** exponent
                    formula = tuple(scalar * (j == degree % n) * (cyclic or degree < n)
                                    for j in range(n))
                    require(columns[a] == formula,
                            "recursive matrix evaluation agrees with class formula")
                for weight, selected in zip(weights, expected):
                    require(greedy(ground, columns, weight, characteristic) == selected,
                            "independent Gaussian greedy agrees with class minima")
                    checks += 1
                field_results.append({"characteristic": characteristic, "variable_scalars": scalars,
                                      "initial_ideal_count": len(witnesses)})
            fixtures.append({
                "n": n, "powers": powers,
                "model": "K[t]/(t^n-1)" if cyclic else "K[t]/(t^n)",
                "kernel_description": "(x^n-epsilon, y-c_y*x^a, z-c_z*x^b), epsilon=1 for cyclic and 0 for nilpotent; (1,a,b)=powers; scalars are recorded per field.",
                "ground_monomials": ground, "ground_size": len(ground),
                "parallel_classes": classes, "zero_residue_monomials": loops,
                "comparison_hyperplanes": normals, "exact_arrangement_chambers": len(weights),
                "initial_ideal_count": len(witnesses), "gaussian_checks": checks,
                "multiplication_matrices_commute": True, "cyclicity_checked": True,
                "relations_have_surviving_nonconstant_multiples": True,
                "field_results": field_results, "witnesses": list(witnesses.values()),
            })
    return fixtures


def nonpartition_point_control():
    """A radical point quotient where independent class selection fails."""
    characteristic = 101
    points = [(t, t ** 2, t ** 3) for t in range(1, 5)]
    monomials = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 2, 0)]
    columns = {a: tuple(pow(t, a[0] + 2 * a[1] + 3 * a[2], characteristic)
                        for t in range(1, 5)) for a in monomials}
    coefficients = (-24, 50, -35, 10, -1)
    require(all(sum(c * columns[a][j] for c, a in zip(coefficients, monomials))
                % characteristic == 0 for j in range(4)), "explicit five-term circuit")
    for omitted in monomials:
        selected = [a for a in monomials if a != omitted]
        require(determinant([columns[a] for a in selected]).denominator == 1,
                "integer determinant")
        require(int(determinant([columns[a] for a in selected])) % characteristic != 0,
                "each four-element subset is independent")
    require(all(_box_volume(a) <= 4 for a in monomials), "circuit lies inside V_4^3")
    return {"n": 4, "points": points, "characteristic": characteristic,
            "monomials": monomials, "circuit_coefficients": coefficients,
            "every_four_columns_independent": True,
            "conclusion": "Five distinct nonzero residue lines have rank four. Selecting one representative from each line does not give a basis. The binomial proof does not apply to all point ideals.",
            "scope": "This refutes only the proposed automatic partition-matroid extension, not every possible Minkowski decomposition."}


def truncated_partition_point_control():
    """A realizable nonbinomial control; count all 86 seven-box staircases."""
    n, d, powers = 7, 3, (1, 6, 9)
    states = {frozenset({(0,) * d})}
    for _ in range(1, n):
        following = set()
        for state in states:
            border = {tuple(x + (i == j) for j, x in enumerate(a))
                      for a in state for i in range(d)} - state
            for a in border:
                if all(tuple(x - (i == j) for j, x in enumerate(a)) in state
                       for i in range(d) if a[i]):
                    following.add(state | {a})
        states = following
    require(len(states) == 86, "all three-dimensional seven-box staircases")
    ground = hyperbolic_cross(n, d)

    def feasible(state, degree):
        representatives = {degree[a]: a for a in state}
        if len(representatives) != n:
            return None
        maximal = [a for a in state
                   if not any(tuple(x + (i == j) for j, x in enumerate(a)) in state
                              for i in range(d))]
        normals = set()
        for b in ground:
            if b in state:
                continue
            targets = ([representatives[degree[b]]] if degree[b] in representatives
                       else maximal)
            for a in targets:
                normal = tuple(x - y for x, y in zip(b, a))
                if max(normal) <= 0:
                    return None
                if min(normal) >= 0:
                    continue
                factor = gcd(*normal)
                normals.add(tuple(x // factor for x in normal))
        polygon = [tuple(Fraction(i == j) for i in range(d)) for j in range(d)]
        for normal in sorted(normals):
            polygon = clip(polygon, normal, 1)
            if not polygon:
                return None
        centroid = tuple(sum(p[i] for p in polygon) / len(polygon) for i in range(d))
        require(min(centroid) > 0 and all(dot(centroid, normal) > 0 for normal in normals),
                "strict feasibility of all selected-class and excluded-class comparisons")
        denominator = lcm(*(x.denominator for x in centroid))
        weight = tuple(int(x * denominator) for x in centroid)
        factor = gcd(*weight)
        return tuple(x // factor for x in weight)

    degree = {a: dot(a, powers) for a in ground}
    points = [(t, t ** 6, t ** 9) for t in range(1, n + 1)]
    columns = {a: tuple(t ** degree[a] for t in range(1, n + 1)) for a in ground}
    generic, witnesses = {}, []
    for state in sorted(states, key=lambda state: tuple(sorted(state))):
        weight = feasible(state, {a: a for a in ground})
        if weight is not None:
            generic[state] = weight
        weight = feasible(state, degree)
        if weight is None:
            continue
        require(greedy(ground, columns, weight, 0) == state,
                "independent rational evaluation agrees with truncated-class inequalities")
        witnesses.append({"positive_integer_weights": weight,
                          "standard_monomials": sorted(state),
                          "is_a_generic_corner_cut": state in generic})
    require(len(generic) == 36 and len(witnesses) == 37,
            "exact bounded control: 37 point initials versus 36 generic corner cuts")
    return {"n": n, "dimension": d, "field": "Q", "positive_sample_parameters": list(range(1, n + 1)),
            "powers": powers, "points": points, "staircases_checked": len(states),
            "generic_corner_cuts": len(generic), "point_initial_ideals": len(witnesses),
            "independent_rational_gaussian_checks": len(witnesses),
            "non_generic_corner_cut_staircases": sum(not w["is_a_generic_corner_cut"] for w in witnesses),
            "witnesses": witnesses,
            "rank_rule": "Over positive real samples, columns of distinct t-degrees form a generalized Vandermonde matrix; rank equals min(n, number of distinct degrees). This is a truncation of the degree-class partition matroid.",
            "scope": "A finite realizable test of a proposed lower-family mechanism, not an asymptotic advance or a novelty claim."}


def pure_power_sector_controls(binomial, point_fixture):
    """Check the boxes supplied by primitive pure-power relations."""
    results = []
    sources = [(f"binomial fixture {index}", f) for index, f in enumerate(binomial)]
    sources.append(("positive monomial-curve points", point_fixture))
    for label, fixture in sources:
        n, powers = fixture["n"], tuple(fixture["powers"])
        ground = hyperbolic_cross(n, 3)
        base = max(powers) + 1
        delta = (1, base, base ** 2)
        scale = 1 + (n - 1) * sum(delta) + 2 * max(powers) * max(delta)
        for index, witness in enumerate(fixture["witnesses"]):
            original = tuple(witness["positive_integer_weights"])
            weight = tuple(scale * x + y for x, y in zip(original, delta))
            for a, b in combinations(ground, 2):
                normal = tuple(x - y for x, y in zip(a, b))
                before, after = dot(original, normal), dot(weight, normal)
                require(not before or before * after > 0,
                        "integer perturbation preserves all strict ground comparisons")
            ratios = [Fraction(w, a) for w, a in zip(weight, powers)]
            require(len(set(ratios)) == 3, "unique degree-normalized cheapest variable")
            axis = min(range(3), key=lambda i: ratios[i])
            sides = [n if j == axis else min(n, powers[axis] // gcd(powers[axis], powers[j]))
                     for j in range(3)]
            selected = set(map(tuple, witness["standard_monomials"]))
            require(all(a[j] < sides[j] for a in selected for j in range(3)),
                    "standard staircase lies in the certified sector box")
            moment = [sum(a[j] for a in selected) for j in range(3)]
            require(all(2 * moment[j] <= n * (sides[j] - 1) for j in range(3)),
                    "downward-column moment bound")
            if label == "positive monomial-curve points":
                columns = {a: tuple(t ** dot(a, powers) for t in range(1, n + 1))
                           for a in ground}
                require(greedy(ground, columns, weight, 0) == selected,
                        "independent point evaluation after removing sector-boundary ties")
            else:
                chosen = {min(map(tuple, group), key=lambda a: (dot(weight, a), a))
                          for group in fixture["parallel_classes"].values()}
                require(chosen == selected, "binomial class minima survive perturbation")
            results.append({"source": label, "witness_index": index, "n": n,
                            "powers": powers, "perturbed_integer_weights": weight,
                            "cheapest_variable_index": axis, "box_sides": sides,
                            "moment_sum": moment,
                            "doubled_moment_coordinate_bounds": [n * (side - 1) for side in sides]})
    return results


def coarse_class_path_control(fixture):
    """Same connected class-selection region, different initial staircases."""
    powers, n = fixture["powers"], fixture["n"]
    ground = hyperbolic_cross(n, 3)
    classes = {}
    for a in ground:
        classes.setdefault(dot(a, powers), []).append(a)
    groups = {}
    for index, witness in enumerate(fixture["witnesses"]):
        degrees = tuple(sorted(dot(a, powers) for a in witness["standard_monomials"]))
        groups.setdefault(degrees, []).append(index)
    target, indices = next((degrees, indices) for degrees, indices in groups.items() if len(indices) > 1)
    indices = indices[:2]
    left, right = [tuple(fixture["witnesses"][i]["positive_integer_weights"]) for i in indices]
    roots = {Fraction(0), Fraction(1)}
    for a, b in combinations(ground, 2):
        normal = tuple(x - y for x, y in zip(a, b))
        start, finish = dot(left, normal), dot(right, normal)
        if start != finish:
            root = Fraction(start, start - finish)
            if 0 < root < 1:
                roots.add(root)
    roots = sorted(roots)
    samples = sorted(set(roots) | {(a + b) / 2 for a, b in zip(roots, roots[1:])})
    for parameter in samples:
        weight = tuple((1 - parameter) * a + parameter * b for a, b in zip(left, right))
        require(min(weight) > 0, "entire parameter path has positive weights")
        minima = sorted((min(dot(weight, a) for a in group), degree)
                        for degree, group in classes.items())
        require(tuple(sorted(degree for _, degree in minima[:n])) == target,
                "same selected degree classes on every order interval and boundary")
        require(minima[n - 1][0] < minima[n][0],
                "continuous halfspace threshold separates selected and unselected classes")
    require(fixture["witnesses"][indices[0]]["standard_monomials"] !=
            fixture["witnesses"][indices[1]]["standard_monomials"],
            "initial ideals change inside one connected coarse region")
    return {"point_initial_ideals": fixture["point_initial_ideals"],
            "distinct_selected_degree_sets": len(groups),
            "same_connected_class_region": {"degree_set": target, "witness_indices": indices,
                                            "endpoint_integer_weights": [left, right],
                                            "complete_order_breakpoints": [str(x) for x in roots],
                                            "open_intervals_checked": len(roots) - 1,
                                            "boundaries_checked": len(roots)},
            "conclusion": "The two different initial staircases are in one connected halfspace-deformation class with the same selected degree classes. A coarse colored k-set count cannot be used directly as the number of initial ideals."}


def normal_curve_lower_controls():
    """Exact allocation certificates for the uniform rational-normal-curve family."""
    characteristic = 1_000_003
    require(all(characteristic % p for p in range(2, 1001)), "prime modulus for independent rank checks")
    fixtures = []
    for d in range(2, 7):
        m = 16 * d
        h = m // (8 * d)
        for remainder in (0, d):
            n = (d + 1) * m + remainder
            certificates = []
            for offsets in product(range(h + 1), repeat=d - 1):
                counts = [m + offset for offset in offsets]
                counts.insert(0, n - sum(counts))
                weight = tuple(d * (counts[0] - counts[i]) + i for i in range(1, d)) + (d,)
                require(sum(counts) == n and min(counts) > 0, "positive allocation of n classes")
                require(counts[0] >= max(counts[1:]), "ray union is a staircase")
                margins = []
                for i in range(1, d):
                    for j in range(i, d):
                        quotient, residue = divmod(i + j, d)
                        tail = quotient * weight[-1] + (weight[residue - 1] if residue else 0)
                        margins.append(weight[i - 1] + weight[j - 1] - tail)
                require(min(margins) > 0, "every product of small variables is strictly reducible")
                last_selected = [d * (counts[0] - 1)] + [weight[i - 1] + d * (counts[i] - 1)
                                                         for i in range(1, d)]
                first_unselected = [d * counts[0]] + [weight[i - 1] + d * counts[i]
                                                     for i in range(1, d)]
                doubled_threshold = 2 * d * counts[0] - 1
                require(2 * max(last_selected) < doubled_threshold < 2 * min(first_unselected),
                        "exact threshold cuts the prescribed n monomials in the d rays")
                certificates.append({"ray_counts": counts, "positive_integer_weights": weight,
                                     "minimum_quadratic_reduction_margin": min(margins),
                                     "doubled_cut_threshold": doubled_threshold})
            require(len({tuple(c["ray_counts"]) for c in certificates}) == (h + 1) ** (d - 1),
                    "different allocations give distinct staircases")
            rank_check = None
            if d <= 3:
                certificate = certificates[-1]
                selected = []
                for residue, count in enumerate(certificate["ray_counts"]):
                    for k in range(count):
                        a = [0] * d
                        a[-1] = k
                        if residue:
                            a[residue - 1] = 1
                        selected.append(tuple(a))
                columns = {a: tuple(pow(t, dot(a, range(1, d + 1)), characteristic)
                                    for t in range(1, n + 1)) for a in selected}
                require(len(greedy(selected, columns, certificate["positive_integer_weights"], characteristic)) == n,
                        "independent modular Vandermonde rank certificate implies nonzero rational determinant")
                rank_check = {"characteristic": characteristic, "rank": n,
                              "allocation": certificate["ray_counts"],
                              "meaning": "The integer positive-point evaluation minor is nonzero modulo a proved prime, so it is also nonzero over Q. No uniform finite-field assertion is made."}
            fixtures.append({"dimension": d, "n": n, "m": m, "remainder": remainder,
                             "offset_range": [0, h], "allocation_certificates_checked": len(certificates),
                             "canonical_certificates_sha256": sha256(json.dumps(certificates, sort_keys=True).encode()).hexdigest(),
                             "first_certificate": certificates[0], "last_certificate": certificates[-1],
                             "independent_modular_rank_check": rank_check})
    return fixtures


def relative_socle_fixture(length):
    """Coupled relation modules killed by x,y, with nonzero z multiples."""
    plane = {(0, 0), (1, 0), (2, 0), (3, 0),
             (0, 1), (0, 2), (0, 3), (1, 1)}
    tops = {(0, 3), (1, 1), (3, 0)}
    domain = {(*a, k) for a in plane for k in range(length)}
    support = sorted((*a, k) for a in tops for k in range(length))
    fixed = domain - set(support)
    basis = sorted(fixed | {(0, 3, k) for k in range(length)})
    n = len(basis)
    require(n == 6 * length, "exact quotient length")
    indices = {p: i for i, p in enumerate(basis)}
    polynomial = [1]
    for root in range(1, length + 1):
        following = [0] * (len(polynomial) + 1)
        for k, value in enumerate(polynomial):
            following[k] -= root * value
            following[k + 1] += value
        polynomial = following
    remainders = [(1,) + (0,) * (length - 1)]
    for _ in range(3 * length):
        current = remainders[-1]
        following = [0] + list(current[:-1])
        for k in range(length):
            following[k] -= current[-1] * polynomial[k]
        remainders.append(tuple(following))
    offsets = {(0, 3): 0, (1, 1): length, (3, 0): 2 * length}

    def residue(point):
        a, k = point[:2], point[2]
        vector = [0] * n
        if a in plane:
            anchor = (0, 3) if a in tops else a
            for j, coefficient in enumerate(remainders[k + offsets.get(a, 0)]):
                vector[indices[(*anchor, j)]] = coefficient
        return tuple(vector)

    matrices = [tuple(residue(tuple(x + (i == j) for j, x in enumerate(a)))
                      for a in basis) for i in range(3)]
    for i, j in combinations(range(3), 2):
        for k in range(n):
            require(matvec(matrices[i], matrices[j][k]) ==
                    matvec(matrices[j], matrices[i][k]), "relative-socle multiplication commutes")
    border = {tuple(x + (i == j) for j, x in enumerate(a))
              for a in domain for i in range(3)} - domain
    columns = {}
    for a in sorted(domain | border):
        vector = residue((0, 0, 0))
        for matrix, exponent in zip(matrices, a):
            for _ in range(exponent):
                vector = matvec(matrix, vector)
        require(vector == residue(a), "relative-socle quotient is cyclic on the domain and border")
        columns[a] = vector
    socle_columns = {a: remainders[a[2] + offsets[a[:2]]] for a in support}
    rows = []
    for a in support:
        if a[:2] == (0, 3):
            continue
        row = [0] * len(support)
        row[support.index(a)] = 1
        for k, value in enumerate(socle_columns[a]):
            row[support.index((0, 3, k))] -= value
        rows.append(row)
    require(any(socle_columns[(1, 1, 1)]), "z*xy remains nonzero over Q")
    budgets = [sum(a[i] for a in support) for i in (0, 1)]
    require(all(b + len(support) <= len(domain) and b <= n for b in budgets),
            "relative-antichain coordinate budgets")
    normals = sorted({primitive(tuple(x - y for x, y in zip(a, b)))
                      for a, b in combinations(support, 2)
                      if min(x - y for x, y in zip(a, b)) < 0 < max(x - y for x, y in zip(a, b))})
    weights = positive_chambers(normals)
    fields, full_checks = [], 0
    evaluation = {a: tuple(t ** (a[2] + offsets[a[:2]]) for t in range(1, length + 1))
                  for a in support}
    for characteristic in (0, 2, 3, 29):
        states = {}
        for weight in weights:
            excluded = leading_row_pivots(rows, support, weight, characteristic)
            selected = greedy(set(support), socle_columns, weight, characteristic)
            require(selected == set(support) - excluded and len(selected) == length,
                    "relative module row reduction and residue greedy agree")
            if not characteristic:
                require(greedy(set(support), evaluation, weight, 0) == selected,
                        "independent positive evaluation agrees with companion residues")
                require(selected == set(sorted(support, key=lambda a: (dot(weight, a), a))[:length]),
                        "positive Vandermonde model chooses the cheapest L socle monomials")
            standard = fixed | selected
            require(all(tuple(x - (i == j) for j, x in enumerate(a)) in standard
                        for a in standard for i in range(3) if a[i]), "relative-module standard staircase")
            moment = tuple(sum(a[i] for a in selected) for i in range(3))
            require(moment[0] <= n and moment[1] <= n and moment[2] <= n * (n - 1) // 2,
                    "translated state lies in the n by n by n-squared box")
            states.setdefault(tuple(sorted(selected)), weight)
        for selected, weight in states.items():
            require(greedy(domain | border, columns, weight, characteristic) == fixed | set(selected),
                    "full commuting-matrix greedy agrees with the relation-module basis")
            full_checks += 1
        fields.append({"characteristic": characteristic, "distinct_initial_staircases": len(states),
                       "all_chambers_checked": len(weights),
                       "witnesses": [{"selected_relative_socle_monomials": s, "positive_integer_weights": w}
                                     for s, w in sorted(states.items())]})
    return {"L": length, "quotient_length": n, "reference_domain": sorted(domain),
            "relative_socle": support, "fixed_standard_monomials": sorted(fixed),
            "monic_z_polynomial_low_degree_first": polynomial,
            "socle_relation_matrix": rows,
            "relative_coordinate_budgets_xy": budgets,
            "multiplication_matrices_commute": True, "cyclicity_checked": True,
            "nonconstant_z_multiples_of_relation_terms_nonzero_over_Q": True,
            "exact_chambers": len(weights), "fields": fields,
            "row_and_greedy_checks": 4 * len(weights), "full_matrix_checks": full_checks}


def relative_socle_lower_controls():
    """A symbolic quadratic-size family of choices in the relative socle."""
    fixtures = []
    for m in (8, 16, 24):
        h = m // 8
        for remainder in range(4):
            length = 4 * m + remainder
            certificates = []
            for left, right in product(range(m, m + h + 1), repeat=2):
                middle = length - left - right
                first, second = middle - left, middle - right
                weight = (3 * first + 6 * second + 5, 6 * first + 3 * second + 4, 9)
                require(min(weight) > 0, "strict positive integer allocation weight")
                offsets = [3 * weight[1], weight[0] + weight[1], 3 * weight[0]]
                counts = [left, middle, right]
                last = max(a + 9 * (r - 1) for a, r in zip(offsets, counts))
                following = min(a + 9 * r for a, r in zip(offsets, counts))
                threshold = 2 * (offsets[1] + 9 * middle) - 1
                require(2 * last < threshold < 2 * following,
                        "the exact cutoff selects the prescribed L relative-socle monomials")
                certificates.append({"selected_ray_counts": counts, "positive_integer_weights": weight,
                                     "doubled_cutoff": threshold})
            require(len({tuple(c["selected_ray_counts"]) for c in certificates}) == (h + 1) ** 2,
                    "distinct allocations give distinct initial staircases")
            fixtures.append({"L": length, "quotient_length": 6 * length, "m": m,
                             "remainder": remainder, "outer_count_range": [m, m + h],
                             "certificates_checked": len(certificates), "certificates": certificates})
    return fixtures


def hook_deformation_fixture(r, s):
    """A flat first-order deformation with both surviving module actions."""
    base = {(i, 0) for i in range(r)} | {(0, j) for j in range(s)}
    generators = [(r, 0), (1, 1), (0, s)]
    square = {tuple(a + b for a, b in zip(p, q))
              for p in generators for q in generators}
    outer_domain = {a for a in product(range(2 * r), range(2 * s))
                    if not any(all(x >= y for x, y in zip(a, g)) for g in square)}
    m = len(base)
    require(m == r + s - 1 and len(outer_domain) == 3 * m + 1, "hook reference lengths")
    domain = {(0, *a) for a in outer_domain} | {(1, *a) for a in base}
    fixed = {(0, *a) for a in base}
    support = sorted(domain - fixed)
    basis = sorted(fixed | {(1, *a) for a in base})
    indices = {a: i for i, a in enumerate(basis)}
    base_order = sorted(base)
    f0 = {(i, 0): i + 1 for i in range(1, r)} | {(0, s - 1): r + 3}
    f1 = {(i, 0): (i + 2) ** 2 for i in range(1, r)} | {
        (0, j): (j + 3) ** 3 for j in range(1, s)}
    f2 = {(r - 1, 0): s + 5} | {(0, j): (j + 5) ** 2 for j in range(1, s)}
    tails = [f0, f1, f2]

    def multiply(poly, monomial):
        return {tuple(a + b for a, b in zip(point, monomial)): value
                for point, value in poly.items()
                if tuple(a + b for a, b in zip(point, monomial)) in base and value}

    require(multiply(f0, (0, 1)) == multiply(f1, (r - 1, 0)), "first adjacent monomial syzygy")
    require(multiply(f1, (0, s - 1)) == multiply(f2, (1, 0)), "second adjacent monomial syzygy")
    invalid = dict(f1)
    invalid[(0, 0)] = 1
    require(multiply(f0, (0, 1)) != multiply(invalid, (r - 1, 0)) and
            multiply(invalid, (0, s - 1)) != multiply(f2, (1, 0)),
            "an arbitrary constant tail violates both module compatibilities")

    def phi(point):
        candidates = [multiply(tail, tuple(a - b for a, b in zip(point, g)))
                      for g, tail in zip(generators, tails)
                      if all(a >= b for a, b in zip(point, g))]
        require(candidates and all(p == candidates[0] for p in candidates),
                "generator-independent module map")
        return candidates[0]

    def residue(point):
        power, a = point[0], point[1:]
        vector = [0] * (2 * m)
        if power < 2 and a in base:
            vector[indices[(power, *a)]] = 1
        elif not power and a in outer_domain:
            for b, value in phi(a).items():
                vector[indices[(1, *b)]] = value
        return tuple(vector)

    matrices = [tuple(residue(tuple(x + (i == j) for j, x in enumerate(a)))
                      for a in basis) for i in range(3)]
    for i, j in combinations(range(3), 2):
        for k in range(2 * m):
            require(matvec(matrices[i], matrices[j][k]) ==
                    matvec(matrices[j], matrices[i][k]), "first-order multiplication actions commute")
    border = {tuple(x + (i == j) for j, x in enumerate(a))
              for a in domain for i in range(3)} - domain
    columns = {}
    for a in sorted(domain | border):
        vector = residue((0, 0, 0))
        for matrix, exponent in zip(matrices, a):
            for _ in range(exponent):
                vector = matvec(matrix, vector)
        require(vector == residue(a), "first-order quotient is cyclic on reference and border")
        columns[a] = vector
    socle_columns = {a: tuple(columns[a][indices[(1, *b)]] for b in base_order) for a in support}
    rows = []
    for point in sorted(outer_domain - base):
        row = [0] * len(support)
        row[support.index((0, *point))] = 1
        for b, value in phi(point).items():
            row[support.index((1, *b))] -= value
        rows.append(row)
    connectors = {(0, r, 0), (0, 0, s), (0, 1, 1)}
    constant = (1, 0, 0)
    left = ({(1, i, 0) for i in range(1, r)} |
            {(0, i, 0) for i in range(r + 1, 2 * r)} |
            {(0, i, 1) for i in range(2, r + 1)})
    right = ({(1, 0, j) for j in range(1, s)} |
             {(0, 0, j) for j in range(s + 1, 2 * s)} |
             {(0, 1, j) for j in range(2, s + 1)})
    require(set(support) == left | right | connectors | {constant}, "two blocks and three connectors")
    for points, allowed in ((left, {(i, 0) for i in range(1, r)}),
                            (right, {(0, j) for j in range(1, s)})):
        require(all(not value or b in allowed for a in points
                    for b, value in zip(base_order, socle_columns[a])), "block residue subspaces are disjoint")
    normals = lambda points: {primitive(tuple(x - y for x, y in zip(a, b)))
                              for a, b in combinations(sorted(points), 2)
                              if min(x - y for x, y in zip(a, b)) < 0 < max(x - y for x, y in zip(a, b))}
    internal = sorted(normals(left) | normals(right))
    require(len(internal) <= 6 * (r + s), "linear number of internal order hyperplanes")
    weights = positive_chambers(sorted(normals(support)))
    internal_regions = len(positive_chambers(internal))
    fields, matrix_checks = [], 0
    for characteristic in (0, 2, 3, 29):
        states, by_order = {}, {}
        for weight in weights:
            excluded = leading_row_pivots(rows, support, weight, characteristic)
            selected = greedy(set(support), socle_columns, weight, characteristic)
            require(selected == set(support) - excluded and len(selected) == m,
                    "first-order relation rows agree with Gaussian residue selection")
            require(constant in selected, "the x constant column is a coloop")
            key = (tuple(sorted(left, key=lambda a: dot(weight, a))),
                   tuple(sorted(right, key=lambda a: dot(weight, a))))
            by_order.setdefault(key, set()).add(tuple(sorted(selected)))
            states.setdefault(tuple(sorted(selected)), weight)
        require(len(by_order) == internal_regions and all(len(v) <= 60 for v in by_order.values()),
                "bounded-connector counting lemma covers every full chamber")
        for selected, weight in states.items():
            require(greedy(domain | border, columns, weight, characteristic) == fixed | set(selected),
                    "full first-order algebra basis agrees with the module reduction")
            matrix_checks += 1
        fields.append({"characteristic": characteristic, "initial_ideal_count": len(states),
                       "maximum_states_per_internal_order_region": max(map(len, by_order.values())),
                       "witnesses": [{"selected_variable_monomials": a, "positive_integer_weights": w}
                                     for a, w in sorted(states.items())]})
    return {"r": r, "s": s, "base_colength": m, "quotient_length": 2 * m,
            "plane_generators": generators,
            "tails": [[{"monomial": a, "coefficient": c} for a, c in sorted(f.items())] for f in tails],
            "reference_domain": sorted(domain), "relation_support": support,
            "relation_matrix": rows, "left_block": sorted(left), "right_block": sorted(right),
            "connectors": sorted(connectors), "internal_hyperplanes": internal,
            "internal_order_regions": internal_regions, "full_comparison_chambers": len(weights),
            "row_and_greedy_checks": 4 * len(weights), "full_matrix_checks": matrix_checks,
            "both_surviving_actions_nonzero_over_Q": bool(multiply(f1, (1, 0))) and bool(multiply(f1, (0, 1))),
            "invalid_constant_tail_rejected": True, "multiplication_matrices_commute": True,
            "cyclicity_checked": True, "fields": fields}


def sparse_corner_deformation_fixture(k, r, s):
    """Bounded controls for the growing-corner class, not a fan enumeration."""
    require(k >= 2 and min(r, s) >= 2 * k + 1, "separated long arms")
    generators = [(r, 0)] + [(k + 1 - i, i) for i in range(1, k + 1)] + [(0, s)]
    divisible = lambda a, g: all(x >= y for x, y in zip(a, g))
    base = {a for a in product(range(r), range(s))
            if not any(divisible(a, g) for g in generators)}
    m = len(base)
    require(m == r + s - 1 + k * (k - 1) // 2, "sparse corner colength")
    square = {tuple(x + y for x, y in zip(a, b)) for a in generators for b in generators}
    outer = {a for a in product(range(2 * r), range(2 * s))
             if not any(divisible(a, g) for g in square)}
    socle = {(r - 1, 0), (0, s - 1)} | {(i, k - i) for i in range(1, k)}
    y_arm = {(i, 0) for i in range(k, r)}
    z_arm = {(0, j) for j in range(k, s)}
    allowed = [y_arm | socle, y_arm | socle] + [socle] * (k - 2) + [z_arm | socle, z_arm | socle]
    require(sum(map(len, allowed)) == 2 * m, "explicit first-order parameter count")
    tails = [{a: (i + 2) * (a[0] + 1) ** 2 + (a[1] + 2) ** 3 for a in points}
             for i, points in enumerate(allowed)]
    base_order = sorted(base)
    unknowns = [(i, a) for i in range(k + 2) for a in base_order]
    equations = []
    for i, (a, b) in enumerate(zip(generators, generators[1:])):
        shifts = [(0, b[1] - a[1]), (a[0] - b[0], 0)]
        for target in base_order:
            equations.append([sum(sign for j, shift, sign in ((i, shifts[0], 1), (i + 1, shifts[1], -1))
                                  if index == j and tuple(x + y for x, y in zip(point, shift)) == target)
                              for index, point in unknowns])
    free_indices = {j for j, (i, a) in enumerate(unknowns) if a in allowed[i]}
    require(all(not row[j] for row in equations for j in free_indices),
            "every displayed free tail satisfies all adjacent syzygies")

    def multiply(poly, monomial):
        return {tuple(x + y for x, y in zip(a, monomial)): c for a, c in poly.items()
                if tuple(x + y for x, y in zip(a, monomial)) in base and c}

    def phi(a):
        images = [multiply(f, tuple(x - y for x, y in zip(a, g)))
                  for g, f in zip(generators, tails) if divisible(a, g)]
        require(images and all(v == images[0] for v in images), "sparse generator-independent module map")
        return images[0]

    fixed = {(0, *a) for a in base}
    basis = sorted(fixed | {(1, *a) for a in base})
    indices = {a: i for i, a in enumerate(basis)}
    domain = {(0, *a) for a in outer} | {(1, *a) for a in base}
    support = sorted(domain - fixed)

    def residue(point):
        power, a = point[0], point[1:]
        vector = [0] * (2 * m)
        if power < 2 and a in base:
            vector[indices[(power, *a)]] = 1
        elif power == 0 and a in outer:
            for b, c in phi(a).items():
                vector[indices[(1, *b)]] = c
        return tuple(vector)

    matrices = [tuple(residue(tuple(x + (i == j) for j, x in enumerate(a)))
                      for a in basis) for i in range(3)]
    for i, j in combinations(range(3), 2):
        for a in range(2 * m):
            require(matvec(matrices[i], matrices[j][a]) == matvec(matrices[j], matrices[i][a]),
                    "sparse first-order multiplication matrices commute")
    border = {tuple(x + (i == j) for j, x in enumerate(a)) for a in domain for i in range(3)} - domain
    columns = {}
    for a in sorted(domain | border):
        vector = residue((0, 0, 0))
        for matrix, exponent in zip(matrices, a):
            for _ in range(exponent):
                vector = matvec(matrix, vector)
        require(vector == residue(a), "sparse first-order cyclicity including border")
        columns[a] = vector
    variable_columns = {a: tuple(columns[a][indices[(1, *b)]] for b in base_order) for a in support}
    rows = []
    for a in sorted(outer - base):
        row = [0] * len(support)
        row[support.index((0, *a))] = 1
        for b, c in phi(a).items():
            row[support.index((1, *b))] -= c
        rows.append(row)
    left = ({(1, *a) for a in y_arm - socle} |
            {(0, i, 0) for i in range(r + 1, 2 * r)} |
            {(0, i, 1) for i in range(k + 1, r + k)})
    right = ({(1, *a) for a in z_arm - socle} |
             {(0, 0, j) for j in range(s + 1, 2 * s)} |
             {(0, 1, j) for j in range(k + 1, s + k)})
    core = {(0, *a) for a in generators} | {(1, *a) for a in socle}
    coloops = {(1, *a) for a in base - y_arm - z_arm - socle}
    blocks = [left, right, core]
    ground = set().union(*blocks)
    require(sum(map(len, blocks)) == len(ground) and ground <= set(support), "three disjoint ground blocks")
    loops = set(support) - ground - coloops
    require(all(not any(variable_columns[a]) for a in loops), "all unassigned columns are loops")
    for points, allowed_residues in ((left, y_arm), (right, z_arm)):
        require(all(not c or b in allowed_residues for a in points
                    for b, c in zip(base_order, variable_columns[a])), "private arm residue support")
    normals = set()
    for block in blocks:
        normals.update(primitive(tuple(x - y for x, y in zip(a, b))) for a, b in combinations(sorted(block), 2)
                       if min(x - y for x, y in zip(a, b)) < 0 < max(x - y for x, y in zip(a, b)))
    require(len(normals) <= 10 * (r + s + k + 1), "linear internal order arrangement size")
    weights = [(1, 2, 3), (3, 2, 1), (2, 3, 1), (1, 3, 2), (7, 2, 3), (11, 3, 2),
               (2, 5, 7), (31, 7, 11), (5, 13, 17), (17, 5, 13), (23, 19, 7), (41, 11, 19)]
    # Preserve every strict small-weight comparison and break its ties by a
    # large-base integer encoding of the bounded reference exponents.
    scale = 4 * max(r, s) + 1
    weights = [tuple(scale ** 3 * x + scale ** i for i, x in enumerate(w)) for w in weights]
    weights += [(1, scale, scale ** 2), (scale ** 2, 1, scale), (scale, scale ** 2, 1)]
    require(len(set(weights)) >= 3 and all(len({dot(w, a) for a in support}) == len(support) for w in weights),
            "strict finite positive weight controls")

    def ordered_selection(points, vectors, characteristic):
        labeled = {(i, 0, 0): vectors[a] for i, a in enumerate(points)}
        return {points[a[0]] for a in greedy(set(labeled), labeled, (1, 1, 1), characteristic)}

    fields = []
    for characteristic in (0, 2, 3, 29):
        equation_columns = {(j, 0, 0): tuple(row[j] for row in equations) for j in range(len(unknowns))}
        tangent_dimension = len(unknowns) - len(greedy(set(equation_columns), equation_columns, (1, 1, 1), characteristic))
        require(tangent_dimension == 2 * m, "displayed tails exhaust the tangent space in each fixture")
        rank = lambda points: len(greedy(set(points), variable_columns, (1, 1, 1), characteristic))
        rank_all = rank(ground)
        interfaces = [rank(block) + rank(ground - block) - rank_all for block in blocks]
        require(all(a <= b for a, b in zip(interfaces, (3, 3, 6))), "bounded interface dimensions")
        witnesses = []
        for weight in weights:
            selected = greedy(set(support), variable_columns, weight, characteristic)
            require(selected == set(support) - leading_row_pivots(rows, support, weight, characteristic),
                    "sparse relation rows agree with greedy residues")
            require(greedy(domain | border, columns, weight, characteristic) == fixed | selected,
                    "sparse full algebra selection agrees with module")
            forced, eligible = set(), set()
            for block, interface in zip(blocks, interfaces):
                order = sorted(block, key=lambda a: dot(weight, a))
                local = greedy(block, variable_columns, weight, characteristic)
                private = ordered_selection(sorted(ground - block) + order, variable_columns, characteristic) & block
                require(private <= local and len(local - private) == interface,
                        "only interface-rank many local basis elements can be optional")
                forced.update(private)
                eligible.update(local)
            require(forced <= selected - coloops <= eligible and coloops <= selected,
                    "global basis lies between forced private columns and local greedy bases")
            require(len(eligible - forced) == sum(interfaces) <= 12, "at most twelve optional columns")
            reduced = ordered_selection(sorted(forced) + sorted(eligible - forced, key=lambda a: dot(weight, a)),
                                        variable_columns, characteristic)
            require(reduced == selected - coloops,
                    "contracting forced columns and ordering only optional columns recovers the basis")
            witnesses.append({"positive_integer_weights": weight, "selected_variable_monomials": sorted(selected),
                              "forced_private_columns": sorted(forced), "optional_columns": sorted(eligible - forced)})
        fields.append({"characteristic": characteristic, "tangent_dimension": tangent_dimension,
                       "interface_dimensions": interfaces, "witnesses": witnesses})
    return {"k": k, "r": r, "s": s, "base_colength": m, "quotient_length": 2 * m,
            "plane_generators": generators, "number_of_corners": len(generators),
            "tails": [[{"monomial": a, "coefficient": c} for a, c in sorted(f.items())] for f in tails],
            "reference_domain": sorted(domain), "relation_support": support, "relation_matrix": rows,
            "blocks": [sorted(b) for b in blocks], "coloop_columns": sorted(coloops), "loop_columns": sorted(loops),
            "internal_hyperplanes": sorted(normals), "weight_controls_per_field": len(weights),
            "row_and_greedy_checks": 4 * len(weights), "full_matrix_checks": 4 * len(weights),
            "interface_flag_checks": 12 * len(weights), "optional_recovery_checks": 4 * len(weights),
            "multiplication_matrices_commute": True,
            "cyclicity_checked": True, "exhaustive_fan_enumeration": False, "fields": fields}


def interface_pairwise_countercontrol():
    points = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    columns = dict(zip(points, ((1, 0), (0, 1), (1, 1))))
    fields = []
    for characteristic in (0, 2, 3, 29):
        rank = lambda subset: len(greedy(set(subset), columns, (1, 2, 3), characteristic))
        require(rank(points) == 2 and all(rank(pair) == 2 for pair in combinations(points, 2)),
                "three distinct lines have zero pairwise intersection but global dependencies")
        witnesses = [{"weight": w, "basis": sorted(greedy(set(points), columns, w, characteristic))}
                     for w in ((1, 2, 3), (1, 3, 2), (3, 1, 2))]
        require(len({tuple(map(tuple, x["basis"])) for x in witnesses}) == 3,
                "the false pairwise-interface count would predict one basis instead of three")
        fields.append({"characteristic": characteristic, "witnesses": witnesses})
    return {"columns": [{"monomial": a, "residue": columns[a]} for a in points],
            "sum_pairwise_intersection_dimensions": 0, "sum_intersections_with_all_other_blocks": 3,
            "joint_internal_order_regions": 1, "fields": fields}


def thick_hook_model(a, b, r, s, tails=None):
    """Construct the exact dual-number algebra used in the thick-hook gate."""
    require(a >= 1 and b >= 1 and r >= 2 * a + 1 and s >= 2 * b + 1, "thick hook parameters")
    generators = [(r, 0), (a, b), (0, s)]
    divisible = lambda p, q: all(x >= y for x, y in zip(p, q))
    base = {p for p in product(range(r), range(s)) if not any(divisible(p, g) for g in generators)}
    m = len(base)
    require(m == r * b + s * a - a * b, "thick hook colength")
    y_arm = {p for p in base if p[0] >= a}
    z_arm = {p for p in base if p[1] >= b}
    y_top = {p for p in y_arm if p[0] >= r - a}
    z_top = {p for p in z_arm if p[1] >= s - b}
    allowed = [y_arm | z_top, y_arm | z_arm, y_top | z_arm]
    require([len(p) for p in allowed] == [r * b, m - a * b, s * a], "full tail coordinate dimensions")
    if tails is None:
        coefficients = Random(271828 + 997 * a + 10007 * b + 100003 * r + 1000003 * s)
        tails = [{p: coefficients.randrange(1, 1009) for p in sorted(points)} for points in allowed]
    require(all(set(f) <= points for f, points in zip(tails, allowed)), "compatible tail support")
    square = {tuple(x + y for x, y in zip(g, h)) for g in generators for h in generators}
    outer = {p for p in product(range(2 * r), range(2 * s))
             if not any(divisible(p, g) for g in square)}
    base_order = sorted(base)

    def multiply(f, shift):
        return {tuple(x + y for x, y in zip(p, shift)): c for p, c in f.items()
                if tuple(x + y for x, y in zip(p, shift)) in base and c}

    require(multiply(tails[0], (0, b)) == multiply(tails[1], (r - a, 0)), "first thick syzygy")
    require(multiply(tails[1], (0, s - b)) == multiply(tails[2], (a, 0)), "second thick syzygy")

    def phi(p):
        images = [multiply(f, tuple(x - y for x, y in zip(p, g)))
                  for g, f in zip(generators, tails) if divisible(p, g)]
        require(images and all(v == images[0] for v in images), "thick generator-independent module map")
        return images[0]

    fixed = {(0, *p) for p in base}
    basis = sorted(fixed | {(1, *p) for p in base})
    indices = {p: i for i, p in enumerate(basis)}
    domain = {(0, *p) for p in outer} | {(1, *p) for p in base}
    support = sorted(domain - fixed)

    def residue(p):
        power, point = p[0], p[1:]
        vector = [0] * (2 * m)
        if power < 2 and point in base:
            vector[indices[(power, *point)]] = 1
        elif not power and point in outer:
            for q, c in phi(point).items():
                vector[indices[(1, *q)]] = c
        return tuple(vector)

    matrices = [tuple(residue(tuple(x + (i == j) for j, x in enumerate(p))) for p in basis) for i in range(3)]
    for i, j in combinations(range(3), 2):
        for p in range(2 * m):
            require(matvec(matrices[i], matrices[j][p]) == matvec(matrices[j], matrices[i][p]),
                    "thick multiplication matrices commute")
    border = {tuple(x + (i == j) for j, x in enumerate(p)) for p in domain for i in range(3)} - domain
    columns = {}
    for p in sorted(domain | border):
        vector = residue((0, 0, 0))
        for matrix, exponent in zip(matrices, p):
            for _ in range(exponent):
                vector = matvec(matrix, vector)
        require(vector == residue(p), "thick cyclicity including reference border")
        columns[p] = vector
    variable = {p: tuple(columns[p][indices[(1, *q)]] for q in base_order) for p in support}
    rows = []
    for p in sorted(outer - base):
        row = [0] * len(support)
        row[support.index((0, *p))] = 1
        for q, c in phi(p).items():
            row[support.index((1, *q))] -= c
        rows.append(row)
    return {"base": base, "m": m, "generators": generators, "allowed": allowed, "tails": tails,
            "y_arm": y_arm, "z_arm": z_arm, "y_top": y_top, "z_top": z_top,
            "base_order": base_order, "domain": domain, "border": border, "fixed": fixed,
            "support": support, "columns": columns, "variable": variable, "rows": rows, "phi": phi}


def thick_hook_fixture(a, b, r, s):
    model = thick_hook_model(a, b, r, s)
    m, support, variable = model["m"], model["support"], model["variable"]
    left = ({(1, *p) for p in model["y_arm"] - model["y_top"]} |
            {(0, r + i, j) for i in range(a, r) for j in range(b)} |
            {(0, a + i, b + j) for i in range(a, r) for j in range(b)})
    right = ({(1, *p) for p in model["z_arm"] - model["z_top"]} |
             {(0, a + i, b + j) for i in range(a) for j in range(b, s)} |
             {(0, i, s + j) for i in range(a) for j in range(b, s)})
    core = ({(1, *p) for p in model["y_top"] | model["z_top"]} |
            {(0, g[0] + i, g[1] + j) for g in model["generators"] for i in range(a) for j in range(b)})
    blocks = [left, right, core]
    ground = set().union(*blocks)
    coloops = {(1, *p) for p in model["base"] - model["y_arm"] - model["z_arm"]}
    loops = set(support) - ground - coloops
    require(sum(map(len, blocks)) == len(ground) and ground <= set(support), "thick block partition")
    require(len(core) == 5 * a * b and len(coloops) == a * b, "thick core and coloop sizes")
    require(all(not any(variable[p]) for p in loops), "remaining thick columns are zero")
    for block, allowed in ((left, model["y_arm"]), (right, model["z_arm"])):
        require(all(not c or q in allowed for p in block for q, c in zip(model["base_order"], variable[p])),
                "thick arm residue support")
    normals = set()
    for block in blocks:
        normals.update(primitive(tuple(x - y for x, y in zip(p, q))) for p, q in combinations(sorted(block), 2)
                       if min(x - y for x, y in zip(p, q)) < 0 < max(x - y for x, y in zip(p, q)))
    require(len(normals) <= 120 * m, "uniform linear internal hyperplane budget")
    scale = 4 * max(r, s) + 1
    weights = [tuple(scale ** 3 * x + scale ** i for i, x in enumerate(w))
               for w in ((1, 2, 3), (3, 2, 1), (2, 3, 1), (1, 3, 2), (7, 2, 3), (2, 5, 7))]
    require(all(len({dot(w, p) for p in support}) == len(support) for w in weights), "strict thick weights")
    unknowns = [(i, p) for i in range(3) for p in model["base_order"]]
    equations = []
    for i, (g, h) in enumerate(zip(model["generators"], model["generators"][1:])):
        for target in model["base_order"]:
            equations.append([sum(sign for index, shift, sign in ((i, (0, h[1] - g[1]), 1),
                                                                  (i + 1, (g[0] - h[0], 0), -1))
                                  if j == index and tuple(x + y for x, y in zip(p, shift)) == target)
                              for j, p in unknowns])
    free = {i for i, (j, p) in enumerate(unknowns) if p in model["allowed"][j]}
    require(len(free) == 2 * m and all(not row[i] for row in equations for i in free), "thick free tangent coordinates")
    fields = []
    for characteristic in (0, 2, 3, 29):
        rank = lambda points: len(greedy(set(points), variable, (1, 1, 1), characteristic))
        total_rank = rank(ground)
        interfaces = [rank(block) + rank(ground - block) - total_rank for block in blocks]
        require(all(x <= c * a * b for x, c in zip(interfaces, (3, 3, 6))), "thick interface dimensions")
        tangent_columns = {(i, 0, 0): tuple(row[i] for row in equations) for i in range(len(unknowns))}
        dimension = len(unknowns) - len(greedy(set(tangent_columns), tangent_columns, (1, 1, 1), characteristic))
        require(dimension == 2 * m, "full thick tangent dimension")
        witnesses = []
        for weight in weights:
            selected = greedy(set(support), variable, weight, characteristic)
            require(selected == set(support) - leading_row_pivots(model["rows"], support, weight, characteristic),
                    "thick relation rows agree with residue greedy")
            require(greedy(model["domain"] | model["border"], model["columns"], weight, characteristic) ==
                    model["fixed"] | selected, "thick full algebra basis agrees with module")
            require(coloops <= selected, "thick bottom identity columns are coloops")
            witnesses.append({"positive_integer_weights": weight,
                              "selected_column_indices": [support.index(p) for p in sorted(selected)]})
        fields.append({"characteristic": characteristic, "interface_dimensions": interfaces,
                       "tangent_dimension": dimension, "witnesses": witnesses})
    residue_classes = {tuple(x % modulus for x, modulus in zip(p, (a, b)))
                       for tail in model["tails"] for p in tail}
    require(len(residue_classes) == a * b, "integer tails mix all power-substitution residue classes")
    return {"a": a, "b": b, "R": r, "S": s, "base_colength": m, "quotient_length": 2 * m,
            "plane_generators": model["generators"],
            "coefficient_seed": 271828 + 997 * a + 10007 * b + 100003 * r + 1000003 * s,
            "coefficient_rule": "Single Python Random instance; randrange(1,1009), traversing each allowed tail support in sorted order. No genericity assertion.",
            "tails": [[{"monomial": p, "coefficient": c} for p, c in sorted(f.items())] for f in model["tails"]],
            "reference_domain": sorted(model["domain"]), "relation_support": support,
            "relation_matrix": model["rows"], "blocks": [sorted(v) for v in blocks],
            "internal_hyperplanes": sorted(normals), "tail_residue_classes_over_Q": sorted(residue_classes),
            "witness_encoding": "selected_column_indices index relation_support",
            "row_and_greedy_checks": 4 * len(weights), "full_matrix_checks": 4 * len(weights),
            "multiplication_matrices_commute": True, "cyclicity_checked": True,
            "exhaustive_fan_enumeration": False, "fields": fields}


def thick_power_substitution_control():
    thin_tails = [{(1, 0): 2, (2, 0): 3, (0, 2): 6},
                  {(1, 0): 9, (2, 0): 16, (0, 1): 64, (0, 2): 125},
                  {(2, 0): 8, (0, 1): 36, (0, 2): 49}]
    thin = thick_hook_model(1, 1, 3, 3, thin_tails)
    a, b = 2, 2
    big_tails = [{(a * p[0], b * p[1]): c for p, c in tail.items()} for tail in thin_tails]
    big = thick_hook_model(a, b, 6, 6, big_tails)
    require(big["m"] == a * b * thin["m"], "substitution length scaling")
    weights = [(1, 31, 961), (961, 1, 31), (31, 961, 1), (101, 17, 23), (19, 101, 37)]
    fields = []
    for characteristic in (0, 2, 3, 29):
        witnesses = []
        for weight in weights:
            small_weight = (weight[0], a * weight[1], b * weight[2])
            small_basis = greedy(thin["domain"] | thin["border"], thin["columns"], small_weight, characteristic)
            expected = {(p[0], a * p[1] + i, b * p[2] + j) for p in small_basis
                        for i in range(a) for j in range(b)}
            actual = greedy(big["domain"] | big["border"], big["columns"], weight, characteristic)
            require(actual == expected, "free power extension preserves every checked coherent basis")
            witnesses.append({"large_weights": weight, "rescaled_small_weights": small_weight,
                              "small_standard_monomials": sorted(small_basis)})
        fields.append({"characteristic": characteristic, "witnesses": witnesses})
    return {"small_quotient_length": 2 * thin["m"], "large_quotient_length": 2 * big["m"],
            "power_scaling": [a, b], "small_tails": [[{"monomial": p, "coefficient": c} for p, c in sorted(f.items())]
                                                     for f in thin_tails],
            "basis_lifting_checks": 4 * len(weights), "fields": fields,
            "scope": "Exact selected-weight controls support the uniform free-extension proof; this does not enumerate either full fan."}


def thick_koszul_countercontrol():
    """Commutativity forces a dependence despite a full support matching."""
    a, b, r, s = 2, 3, 5, 7
    f = {(0, 0): 1, (1, 0): 2, (0, 1): 3}
    g = {(0, 0): 5, (1, 0): 7, (0, 1): 11}
    shifts = [(0, 0), (1, 0), (0, 1)]
    model = thick_hook_model(a, b, r, s,
                             [{(a + p[0], p[1]): c for p, c in f.items()},
                              {(a + p[0], p[1]): c for p, c in g.items()}, {}])
    points = [(0, r + p[0], p[1]) for p in shifts] + [(0, a + p[0], b + p[1]) for p in shifts]
    outputs = sorted({(a + p[0] + q[0], p[1] + q[1]) for p in shifts for q in shifts})
    matrix = [[model["phi"](p[1:]).get(q, 0) for p in points] for q in outputs]
    columns = {p: tuple(row[i] for row in matrix) for i, p in enumerate(points)}
    matching = next(perm for perm in permutations(range(6)) if all(matrix[perm[i]][i] for i in range(6)))
    relation = [g[p] for p in shifts] + [-f[p] for p in shifts]
    require(all(sum(c * x for c, x in zip(relation, row)) == 0 for row in matrix),
            "G*F-F*G is an exact integer relation between allowed thick columns")
    basis = sorted(greedy(set(points), columns, (1, 17, 31), 0))
    require(len(basis) == 5 and determinant(matrix) == 0, "support matching rank six differs from true rank five")
    selected = [points.index(p) for p in basis]
    minor_rows = next(rows for rows in combinations(range(6), 5)
                      if determinant([[matrix[i][j] for j in selected] for i in rows]))
    minor = determinant([[matrix[i][j] for j in selected] for i in minor_rows])
    return {"parameters": [a, b, r, s], "quotient_length": 2 * model["m"],
            "F": [{"monomial": p, "coefficient": c} for p, c in sorted(f.items())],
            "G": [{"monomial": p, "coefficient": c} for p, c in sorted(g.items())],
            "input_monomials": points, "output_monomials_after_removing_x": outputs,
            "residue_matrix_rows": matrix, "full_support_matching_rows_by_column": matching,
            "integer_syzygy_coefficients": relation, "rational_rank": 5, "support_matching_size": 6,
            "nonzero_minor": {"rows": minor_rows, "columns": selected, "determinant": str(minor)},
            "uniform_constraint": "For any coefficients, the shifted columns of two polynomials satisfy G*F-F*G=0 whenever these products fit without truncation. The displayed integer specialization has rank exactly five.",
            "scope": "A classical Koszul relation used as a realizability control, not a novel syzygy theorem or an asymptotic lower bound."}


def untruncated_convolution_controls():
    """Check the principal-syzygy rank and basis formula on rectangles."""
    f = {(0, 0): 1, (1, 0): 1, (0, 1): 1}
    g = {(0, 0): 1, (1, 0): 2, (0, 1): 3}
    offsets = [(0, 2, 0), (1, 0, 0)]
    weights = [(1, 13, 169), (169, 1, 13), (13, 169, 1), (19, 5, 7), (5, 19, 7), (7, 5, 19),
               (-13, 1, 2), (1, -13, 2), (1, 2, -13), (-1, -13, -169)]
    results = []
    for sizes in (((2, 2), (2, 2)), ((4, 3), (4, 3)), ((3, 3), (4, 2))):
        domains = [set(product(range(u), range(v))) for u, v in sizes]
        inputs = [(i, p) for i, domain in enumerate(domains) for p in sorted(domain)]
        labels = [tuple(x + y for x, y in zip(offsets[i], (0, *p))) for i, p in inputs]
        outputs = sorted({tuple(x + y for x, y in zip(p, q))
                          for i, p in inputs for q in (f, g)[i]})
        fields = []
        for characteristic in (0, 2, 3, 29):
            polys = [{p: c for p, c in h.items() if not characteristic or c % characteristic}
                     for h in (f, g)]
            degrees = [tuple(max(p[j] for p in h) for j in range(2)) for h in polys]
            h_sizes = [max(0, min(sizes[0][j] - degrees[1][j], sizes[1][j] - degrees[0][j]))
                       for j in range(2)]
            h_domain = set(product(*(range(v) for v in h_sizes)))
            columns = {label: tuple(polys[i].get(tuple(x - y for x, y in zip(q, p)), 0) for q in outputs)
                       for label, (i, p) in zip(labels, inputs)}
            expected_rank = len(inputs) - len(h_domain)
            require(len(greedy(set(labels), columns, (1, 13, 169), characteristic)) == expected_rank,
                    "untruncated convolution rank equals domain count minus eroded-domain count")
            witnesses = []
            for weight in weights:
                syzygy_terms = [(tuple(x + y for x, y in zip(offsets[0], (0, *p))), 0, p) for p in polys[1]]
                syzygy_terms += [(tuple(x + y for x, y in zip(offsets[1], (0, *p))), 1, p) for p in polys[0]]
                lead_label, side, monomial = max(syzygy_terms, key=lambda item: (dot(weight, item[0]), item[0]))
                deleted = {tuple(x + y for x, y in zip(offsets[side], (0, h[0] + monomial[0], h[1] + monomial[1])))
                           for h in h_domain}
                selected = greedy(set(labels), columns, weight, characteristic)
                require(deleted <= set(labels) and selected == set(labels) - deleted,
                        "tagged Newton leading term translates the entire kernel domain")
                total_moment = tuple(sum(p[j] for p in labels) for j in range(3))
                kernel_moment = (0, sum(p[0] for p in h_domain), sum(p[1] for p in h_domain))
                predicted = tuple(x - y - len(h_domain) * z for x, y, z in zip(total_moment, kernel_moment, lead_label))
                require(tuple(sum(p[j] for p in selected) for j in range(3)) == predicted,
                        "projected basis polytope has the translated reflected Newton support function")
                witnesses.append({"weights": weight, "winning_tagged_term": lead_label,
                                  "deleted_column_indices": [labels.index(p) for p in sorted(deleted)]})
            fields.append({"characteristic": characteristic, "kernel_rectangle_sizes": h_sizes,
                           "rank": expected_rank, "polytope_scale": len(h_domain),
                           "polytope_translation": [x - y for x, y in zip(total_moment, kernel_moment)],
                           "tagged_support_points": sorted({t[0] for t in syzygy_terms}), "witnesses": witnesses})
        results.append({"multiplier_rectangle_sizes": sizes,
                        "input_labels": labels, "output_monomials": outputs, "fields": fields})
    return {"F": [{"monomial": p, "coefficient": c} for p, c in sorted(f.items())],
            "G": [{"monomial": p, "coefficient": c} for p, c in sorted(g.items())],
            "input_offsets": offsets, "fixtures": results, "exact_basis_checks": 4 * len(weights) * len(results),
            "uniform_conclusion": "For coprime nonzero F,G and convex lattice-complete multiplier domains, the kernel is (H*G,-H*F) with H in the Newton-polytope erosion domain. Its translated leading terms determine every coherent basis, and the projected basis polytope is a translate of minus the kernel dimension times the tagged Newton polytope. Common factors may be removed before applying the formula. No quotient truncation is allowed.",
            "priority": "An elementary consequence of the classical two-generator syzygy and Newton-product identity; novelty is not established."}


def factorial_boundary_fixture(length):
    """An exact finite-algebra family with an unbounded boundary fan effect."""
    constant = factorial(length - 1)
    coefficients = [constant // factorial(j) for j in range(length)]
    fixed = {(0, 0, j) for j in range(length)}
    basis = sorted(fixed | {(1, 0, j) for j in range(length)})
    indices = {p: i for i, p in enumerate(basis)}
    domain = {(0, 0, j) for j in range(2 * length)} | {(1, 0, j) for j in range(length)}
    support = sorted(domain - fixed)

    def residue(p):
        x, y, z = p
        vector = [0] * (2 * length)
        if not y and x < 2 and z < length:
            vector[indices[p]] = 1
        elif not x and not y and length <= z < 2 * length:
            for j in range(2 * length - z):
                vector[indices[(1, 0, j + z - length)]] = coefficients[j]
        return tuple(vector)

    matrices = [tuple(residue(tuple(x + (i == j) for j, x in enumerate(p))) for p in basis) for i in range(3)]
    for i, j in combinations(range(3), 2):
        for h in range(2 * length):
            require(matvec(matrices[i], matrices[j][h]) == matvec(matrices[j], matrices[i][h]),
                    "factorial boundary multiplication matrices commute")
    border = {tuple(x + (i == j) for j, x in enumerate(p)) for p in domain for i in range(3)} - domain
    columns = {}
    for p in sorted(domain | border):
        value = residue((0, 0, 0))
        for matrix, exponent in zip(matrices, p):
            for _ in range(exponent):
                value = matvec(matrix, value)
        require(value == residue(p), "factorial boundary cyclicity including border")
        columns[p] = value
    variable = {p: tuple(columns[p][indices[(1, 0, j)]] for j in range(length)) for p in support}
    rows = []
    for j in range(length):
        row = [0] * (2 * length)
        row[support.index((0, 0, length + j))] = 1
        for h in range(length - j):
            row[support.index((1, 0, h + j))] -= coefficients[h]
        rows.append(row)
    untruncated = {}
    for p in support:
        if p[0]:
            untruncated[p] = tuple(int(i == p[2]) for i in range(2 * length - 1))
        else:
            shift = p[2] - length
            untruncated[p] = tuple(coefficients[i - shift] if 0 <= i - shift < length else 0
                                   for i in range(2 * length - 1))
    minors = []
    for a in range(length + 1):
        for b in range(length - a + 1):
            matrix = [[Fraction(constant, factorial(a + i - j)) if a + i >= j else 0
                       for j in range(b)] for i in range(b)]
            formula = Fraction(constant ** b)
            for j in range(b):
                formula *= Fraction(factorial(j), factorial(a + j))
            require(determinant(matrix) == formula > 0, "factorial Toeplitz prefix determinant")
            minors.append({"identity_prefix": a, "convolution_prefix": b, "determinant": str(formula)})
    fields, full_checks = [], 0
    for characteristic in (0, 2, 3, 29):
        states, untruncated_states = {}, set()
        uniform = characteristic == 0 or characteristic > length - 1
        for k in range(2 * length):
            weight = (2 * k + 1, 1, 2)
            selected = greedy(set(support), variable, weight, characteristic)
            require(selected == set(support) - leading_row_pivots(rows, support, weight, characteristic),
                    "factorial boundary rows agree with residue selection")
            states.setdefault(tuple(sorted(selected)), weight)
            if uniform:
                b = (2 * k + 3) // 4
                expected = {(1, 0, j) for j in range(length - b)} | {(0, 0, length + j) for j in range(b)}
                require(selected == expected, "proved factorial boundary allocation")
                moment = [sum(p[i] for p in fixed | selected) for i in range(3)]
                require(moment == [length - b, 0, length * (length - 1) + b * b],
                        "factorial boundary quadratic state-moment chain")
            old = greedy(set(support), untruncated, weight, characteristic)
            require(len(old) == 2 * length - 1, "untruncated two-copy rank")
            untruncated_states.add(tuple(sorted(old)))
        require(len(untruncated_states) == 2, "untruncated positive Newton fan has two bases")
        if uniform:
            require(len(states) == length + 1, "factorial boundary exact positive fan count")
        for selected, weight in states.items():
            actual = greedy(domain | border, columns, weight, characteristic)
            require(actual == fixed | set(selected), "factorial full algebra agrees with the variable block")
            require(all(tuple(x - (i == j) for j, x in enumerate(p)) in actual
                        for p in actual for i in range(3) if p[i]), "factorial selected set is a staircase")
            full_checks += 1
        fields.append({"characteristic": characteristic, "exact_initial_ideal_count": len(states),
                       "uniform_characteristic_hypothesis": uniform, "untruncated_positive_basis_count": 2,
                       "witnesses": [{"weights": w, "selected_column_indices": [support.index(p) for p in selected]}
                                     for selected, w in sorted(states.items())]})
    same_cell = None
    if length >= 3:
        weights = [(5, 1, 2), (9, 1, 2)]
        old = [greedy(set(support), untruncated, w, 0) for w in weights]
        new = [greedy(set(support), variable, w, 0) for w in weights]
        require(old[0] == old[1] and new[0] != new[1],
                "one untruncated Newton cell contains two different boundary bases")
        same_cell = {"weights": weights, "untruncated_selected_column_indices": [support.index(p) for p in sorted(old[0])],
                     "boundary_selected_column_indices": [[support.index(p) for p in sorted(v)] for v in new]}
    return {"L": length, "quotient_length": 2 * length, "integer_exponential_coefficients": coefficients,
            "ideal": "(y,x^2,z^L-x*A_L(z)), A_L=(L-1)!*sum_{j=0}^{L-1}z^j/j!",
            "reference_domain": sorted(domain), "relation_support": support, "relation_matrix": rows,
            "prefix_minors": minors, "exact_internal_order_intervals": 2 * length,
            "row_and_greedy_checks": 8 * length, "full_matrix_checks": full_checks,
            "multiplication_matrices_commute": True, "cyclicity_checked": True,
            "exhaustive_positive_fan_covered": True, "same_newton_cell_countercontrol": same_cell, "fields": fields}


def independent_boundary_arms_fixture(t, one_arm_counts):
    constant = factorial(t - 1)
    coefficients = [constant // factorial(j) for j in range(t)]
    model = thick_hook_model(1, 1, t + 1, t + 1,
                             [{(j + 1, 0): c for j, c in enumerate(coefficients)}, {},
                              {(0, j + 1): c for j, c in enumerate(coefficients)}])
    support, variable = model["support"], model["variable"]
    left = {(1, j, 0) for j in range(1, t + 1)} | {(0, t + 1 + j, 0) for j in range(t)}
    right = {(1, 0, j) for j in range(1, t + 1)} | {(0, 0, t + 1 + j) for j in range(t)}
    coloop = (1, 0, 0)
    loops = set(support) - left - right - {coloop}
    require(all(not any(variable[p]) for p in loops), "independent arm leftover columns are loops")
    for block, allowed in ((left, {(j, 0) for j in range(1, t + 1)}),
                            (right, {(0, j) for j in range(1, t + 1)})):
        require(all(not c or q in allowed for p in block for q, c in zip(model["base_order"], variable[p])),
                "independent arm residue spaces are disjoint")
    fields, full_checks = [], 0
    for characteristic in (0, 2, 3, 29):
        states, left_states, right_states = {}, set(), set()
        uniform = characteristic == 0 or characteristic > t - 1
        for i, j in product(range(2 * t), repeat=2):
            weight = ((2 * i + 1) * (2 * j + 1), 2 * (2 * j + 1), 2 * (2 * i + 1))
            selected = greedy(set(support), variable, weight, characteristic)
            require(selected == set(support) - leading_row_pivots(model["rows"], support, weight, characteristic),
                    "independent boundary arm rows agree with residue selection")
            require(coloop in selected, "the constant x column remains a coloop")
            states.setdefault(tuple(sorted(selected)), weight)
            left_states.add(tuple(sorted(selected & left)))
            right_states.add(tuple(sorted(selected & right)))
            if uniform:
                b, c = (2 * i + 3) // 4, (2 * j + 3) // 4
                expected = ({coloop} | {(1, h, 0) for h in range(1, t - b + 1)} |
                            {(0, t + 1 + h, 0) for h in range(b)} |
                            {(1, 0, h) for h in range(1, t - c + 1)} |
                            {(0, 0, t + 1 + h) for h in range(c)})
                require(selected == expected, "independent factorial arm allocation")
                moment = [sum(p[h] for p in model["fixed"] | selected) for h in range(3)]
                require(moment == [2 * t + 1 - b - c, t * (t + 1) + b * b, t * (t + 1) + c * c],
                        "independent arm quadratic moment surface")
        require(len(left_states) == len(right_states) == one_arm_counts[characteristic] and
                len(states) == one_arm_counts[characteristic] ** 2,
                "independent arm counts equal the product of the separately computed one-arm counts")
        if uniform:
            require(len(states) == (t + 1) ** 2, "sharp quadratic hook-family count")
        for selected, weight in states.items():
            require(greedy(model["domain"] | model["border"], model["columns"], weight, characteristic) ==
                    model["fixed"] | set(selected), "independent arm full quotient agrees with module")
            full_checks += 1
        fields.append({"characteristic": characteristic, "exact_initial_ideal_count": len(states),
                       "one_arm_count": one_arm_counts[characteristic], "uniform_characteristic_hypothesis": uniform,
                       "witnesses": [{"weights": w, "selected_column_indices": [support.index(p) for p in selected]}
                                     for selected, w in sorted(states.items())]})
    return {"t": t, "quotient_length": 4 * t + 2, "integer_exponential_coefficients": coefficients,
            "ideal": "(x^2,y*z,y^(t+1)-x*y*A_t(y),z^(t+1)-x*z*A_t(z))",
            "reference_domain": sorted(model["domain"]), "relation_support": support,
            "relation_matrix": model["rows"], "arm_blocks": [sorted(left), sorted(right)],
            "exact_internal_order_cells": 4 * t * t, "row_and_greedy_checks": 16 * t * t,
            "full_matrix_checks": full_checks, "exhaustive_positive_fan_covered": True, "fields": fields}


def finite_staircases(n, d):
    """All cardinality-n monomial downsets, with no weight-feasibility filter."""
    states = {frozenset({(0,) * d})}
    for _ in range(1, n):
        following = set()
        for state in states:
            border = {tuple(x + (i == j) for j, x in enumerate(a))
                      for a in state for i in range(d)} - state
            for a in border:
                if all(tuple(x - (i == j) for j, x in enumerate(a)) in state
                       for i in range(d) if a[i]):
                    following.add(state | {a})
        states = following
    return sorted(states, key=lambda state: tuple(sorted(state)))


def univariate_product_mod(left, right, modulus):
    """Integer normal forms, independent of confluent evaluation formulas."""
    n = len(modulus) - 1
    require(modulus[-1] == 1, "monic univariate modulus")
    values = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            values[i + j] += x * y
    values.extend([0] * max(0, n - len(values)))
    for i in range(len(values) - 1, n - 1, -1):
        coefficient = values[i]
        for j in range(n):
            values[i - n + j] -= coefficient * modulus[j]
    return tuple(values[:n])


def univariate_power_mod(base, exponent, modulus):
    n = len(modulus) - 1
    result = (1,) + (0,) * (n - 1)
    while exponent:
        if exponent & 1:
            result = univariate_product_mod(result, base, modulus)
        base = univariate_product_mod(base, base, modulus)
        exponent //= 2
    return result


def cyclic_univariate_model(modulus, images, ground):
    """Construct actual commuting multiplication matrices and residue columns."""
    n, d = len(modulus) - 1, len(images)
    unit = (1,) + (0,) * (n - 1)
    coordinates = [tuple(int(i == j) for i in range(n)) for j in range(n)]
    matrices = [tuple(univariate_product_mod(image, b, modulus) for b in coordinates)
                for image in images]
    for i, j in combinations(range(d), 2):
        for b in coordinates:
            require(matvec(matrices[i], matvec(matrices[j], b)) ==
                    matvec(matrices[j], matvec(matrices[i], b)),
                    "univariate quotient multiplication matrices commute")
    columns = {}
    for a in sorted(ground, key=lambda a: (sum(a), a)):
        if not any(a):
            columns[a] = unit
        else:
            i = next(i for i, value in enumerate(a) if value)
            previous = tuple(value - (i == j) for j, value in enumerate(a))
            columns[a] = matvec(matrices[i], columns[previous])
    axis = [tuple(i if j == 0 else 0 for j in range(d)) for i in range(n)]
    require(determinant([columns[a] for a in axis]) != 0,
            "first coordinate generates the claimed length-n cyclic quotient")
    return columns


def monomial_curve_collision_controls(point_fixture):
    """Confluent evaluation, companion normal forms and all saved fan witnesses."""
    n, d, powers = point_fixture["n"], point_fixture["dimension"], point_fixture["powers"]
    ground, states = hyperbolic_cross(n, d), finite_staircases(n, d)
    profiles = [tuple((t, 1) for t in range(1, n + 1)), ((1, n),),
                ((1, 3), (2, 4)), ((1, 1), (3, 2), (5, 4))]
    results = []
    for profile in profiles:
        modulus = [1]
        for node, multiplicity in profile:
            for _ in range(multiplicity):
                following = [0] * (len(modulus) + 1)
                for i, value in enumerate(modulus):
                    following[i] -= node * value
                    following[i + 1] += value
                modulus = following
        images = [univariate_power_mod((0, 1), exponent, modulus) for exponent in powers]
        columns = cyclic_univariate_model(modulus, images, ground)
        jets = {}
        for a in ground:
            degree = dot(a, powers)
            jets[a] = tuple(comb(degree, j) * node ** (degree - j) if degree >= j else 0
                            for node, multiplicity in profile for j in range(multiplicity))
            normal_form_jets = tuple(sum(value * comb(i, j) * node ** (i - j)
                                        for i, value in enumerate(columns[a]) if i >= j)
                                     for node, multiplicity in profile for j in range(multiplicity))
            require(normal_form_jets == jets[a], "direct Hermite jets equal companion normal forms")
        ranks = []
        for state in states:
            rank = len(greedy(state, columns, (1, 2, 3), 0))
            require(rank == min(n, len({dot(a, powers) for a in state})),
                    "every staircase rank obeys the confluent degree-class formula")
            ranks.append(rank)
        for witness in point_fixture["witnesses"]:
            weight = witness["positive_integer_weights"]
            expected = {tuple(a) for a in witness["standard_monomials"]}
            require(greedy(ground, columns, weight, 0) == expected and
                    greedy(ground, jets, weight, 0) == expected,
                    "all saved initial ideals survive every checked collision profile")
        results.append({"positive_nodes_and_multiplicities": profile,
                        "monic_modulus_low_to_high": modulus,
                        "staircase_rank_checks": len(states),
                        "ranks_in_lexicographic_staircase_order": ranks,
                        "companion_and_jet_column_checks": len(ground),
                        "fan_witnesses_checked_in_both_representations": len(point_fixture["witnesses"]),
                        "exact_initial_ideal_count_by_uniform_proof": point_fixture["point_initial_ideals"]})
    return {"n": n, "powers": powers, "field": "Q", "profiles": results,
            "scope": "Finite checks use all seven-box staircases and all 37 saved coherent bases. Equality for every positive divisor of degree n is proved in the note, not inferred from these four profiles."}


def monogenic_embedding_controls():
    """Mixed-radix embeddings of fixed monogenic algebras, including nonsplit ones."""
    fixtures = []
    for n, d, modulus in [(3, 2, [0, 0, 0, 1]), (3, 3, [-1, 1, -1, 1]),
                           (4, 2, [0, 0, 0, 0, 1]), (4, 3, [1, 0, 2, 0, 1])]:
        ground, powers = sorted(hyperbolic_cross(n, d)), tuple(n ** i for i in range(d))
        degrees = {a: dot(a, powers) for a in ground}
        require(len(set(degrees.values())) == len(ground), "mixed-radix degree injection on hyperbolic cross")
        rejected, shift = [], 1
        while True:
            columns = {a: univariate_power_mod((shift, 1), degrees[a], modulus) for a in ground}
            first_zero, smallest, checked = None, None, 0
            for subset in combinations(ground, n):
                value = determinant([columns[a] for a in subset])
                checked += 1
                if not value:
                    first_zero = subset
                    break
                smallest = abs(value) if smallest is None else min(smallest, abs(value))
            if first_zero is None:
                break
            rejected.append({"shift": shift, "zero_minor_monomials": first_zero})
            shift += 1
            require(shift <= 12, "bounded fixture search; not a bound in the uniform theorem")
        images = [univariate_power_mod((shift, 1), exponent, modulus) for exponent in powers]
        actual = cyclic_univariate_model(modulus, images, ground)
        require(actual == columns, "commuting cyclic model equals direct shifted-power formula")
        fixtures.append({"n": n, "dimension": d, "powers": powers,
                         "monic_modulus_low_to_high": modulus, "working_shift": shift,
                         "rejected_shift_certificates": rejected, "ground_size": len(ground),
                         "all_rank_n_minors_checked": checked,
                         "smallest_absolute_minor": str(smallest)})
    centered = []
    for d, maximum in [(2, 7), (3, 7), (4, 5)]:
        for n in range(2, maximum + 1):
            modulus = [0] * n + [1]
            powers, ground = tuple(n ** i for i in range(d)), hyperbolic_cross(n, d)
            images = []
            for exponent in powers:
                values = list(univariate_power_mod((1, 1), exponent, modulus))
                values[0] -= 1
                images.append(tuple(values))
            columns = cyclic_univariate_model(modulus, images, ground)
            states = finite_staircases(n, d)
            for state in states:
                ordered = sorted(state)
                plain = [tuple(comb(dot(a, powers), j) for j in range(n)) for a in ordered]
                value = determinant(plain)
                degree_list = [dot(a, powers) for a in ordered]
                formula = Fraction(1, 1)
                for i, j in combinations(range(n), 2):
                    formula *= degree_list[j] - degree_list[i]
                formula /= product_factorials(n)
                require(value == formula != 0, "exact confluent Vandermonde determinant at one")
                require(determinant([columns[a] for a in ordered]) == value,
                        "centering preserves every downset basis determinant")
            centered.append({"n": n, "dimension": d, "powers": powers,
                             "all_staircase_determinants_checked": len(states)})
    # The uncentered finite-ground matroid is uniform; the centered one need not be.
    n, d, powers = 4, 3, (1, 4, 16)
    modulus = [0] * n + [1]
    images = [tuple([0] + [comb(exponent, j) for j in range(1, n)]) for exponent in powers]
    columns = cyclic_univariate_model(modulus, images, hyperbolic_cross(n, d))
    dependent = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0))
    require(determinant([columns[a] for a in dependent]) == 0,
            "translation preserves staircase bases, not the whole monomial matroid")
    small_field_basis = ((0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 3, 0))
    require(len(greedy(small_field_basis, columns, (1, 2, 3), 2)) == 1,
            "characteristic two destroys the y-axis basis for these mixed-radix powers")
    leading = []
    for modulus, degrees in [([-1, 1, -1, 1], (0, 2, 5)), ([1, 0, 2, 0, 1], (0, 1, 3, 6))]:
        n = len(modulus) - 1
        degree = sum(degrees) - n * (n - 1) // 2
        values = [determinant([univariate_power_mod((shift, 1), exponent, modulus)
                               for exponent in degrees]) for shift in range(degree + 2)]
        differences = values
        for _ in range(degree):
            differences = [b - a for a, b in zip(differences, differences[1:])]
        formula = Fraction(1, product_factorials(n))
        for i, j in combinations(range(n), 2):
            formula *= degrees[j] - degrees[i]
        require(differences == [formula * factorial(degree)] * 2,
                "exact forward differences verify determinant leading coefficient")
        leading.append({"monic_modulus_low_to_high": modulus, "degrees": degrees,
                        "polynomial_degree": degree, "leading_coefficient": str(formula),
                        "determinants_at_consecutive_shifts_starting_zero": [str(x) for x in values]})
    return {"field_scope": "Uniform statements are over characteristic zero; the characteristic-two failure is deliberate.",
            "fixed_algebra_fixtures": fixtures, "centered_curvilinear_fixtures": centered,
            "leading_coefficient_controls": leading,
            "centered_nonuniform_matroid_control": {"n": 4, "powers": powers, "dependent_monomials": dependent},
            "small_characteristic_control": {"characteristic": 2, "failed_y_axis_rank": 1},
            "scope": "Classical Vandermonde and open-minor consequences, not a new unrestricted counting theorem or an established novelty claim."}


def product_factorials(n):
    result = 1
    for j in range(n):
        result *= factorial(j)
    return result


def sparse_prefix_interpolation_control():
    """Actual quotient-module prefixes with a symbolic Hermite rank defect."""
    fixtures = []
    for k in (2, 4, 6):
        r = 5 * k // 2
        domain = sorted((i, j) for i in range(r) for j in range(r) if i < k or j < k)
        n, scale = len(domain), factorial(r - 1) ** 2
        weight = (8 * k + 1, 8 * k + 2, 8 * k + 4)
        inputs = [(tag,) + a for tag in (0, 1) for a in domain]
        prefix = sorted(inputs, key=lambda a: (dot(a, weight), a))[:n]
        identity = [a[1:] for a in prefix if a[0] == 0]
        convol = [a[1:] for a in prefix if a[0] == 1]
        require(set(identity) == {a for a in domain if sum(a) < 2 * k} and
                set(convol) == {a for a in domain if sum(a) < 2 * k - 1},
                "strict geometric prefix has the stated two total-degree cuts")
        require(n == 4 * k * k and len(identity) == k * (2 * k + 1) and
                len(convol) == k * (2 * k - 1), "uniform sparse-prefix cardinalities")
        outside = set(inputs) - set(prefix)
        require(max(dot(a, weight) for a in prefix) < min(dot(a, weight) for a in outside),
                "no tie crosses the selected prefix boundary")
        polynomial = {}
        for i in range(2 * k):
            coefficient = comb(2 * k - 1, i) * (-1) ** (2 * k - 1 - i)
            for a in ((i + 1, 2 * k - 1 - i), (i, 2 * k - i)):
                polynomial[a] = polynomial.get(a, 0) + coefficient
        polynomial = {a: value for a, value in polynomial.items() if value}
        require(set(polynomial) <= set(domain) and (k, k) not in polynomial,
                "(y-z)^(2k-1)*(y+z) has allowed sparse-staircase support")
        dual = tuple(polynomial.get(a, 0) * factorial(a[0]) * factorial(a[1]) for a in domain)
        require(any(dual) and all(polynomial.get(a, 0) == 0 for a in identity),
                "nonzero dual polynomial annihilates the identity prefix")
        columns = {}
        for a in identity:
            columns[(0,) + a] = tuple(int(a == row) for row in domain)
        for a in convol:
            columns[(1,) + a] = tuple(
                scale // (factorial(row[0] - a[0]) * factorial(row[1] - a[1]))
                if all(x >= y for x, y in zip(row, a)) else 0 for row in domain)
        require(all(dot(dual, columns[a]) == 0 for a in prefix),
                "factorial-scaled dual gives an exact nonzero left null vector")
        record = {"k": k, "R": r, "quotient_length": n,
                  "monomial_ideal": [[r, 0], [k, k], [0, r]],
                  "positive_tag_y_z_weights": weight,
                  "identity_prefix_size": len(identity), "convolution_prefix_size": len(convol),
                  "convolution_scale": scale,
                  "dual_polynomial_coefficients": [{"monomial": a, "coefficient": v} for a, v in sorted(polynomial.items())],
                  "exact_prefix_null_vector": dual,
                  "rank_upper_bound_by_nonzero_null_vector": n - 1}
        if k == 2:
            rows = [a for a in domain if a not in identity]
            matrix = [[columns[(1,) + a][domain.index(row)] for a in convol] for row in rows]
            residual = {(1,) + a: tuple(matrix[i][j] for i in range(len(rows)))
                        for j, a in enumerate(convol)}
            selected = sorted(greedy(residual, residual, weight, 0))
            require(len(selected) == 5 and len(greedy(prefix, columns, weight, 0)) == 15,
                    "six-by-six residual has rank five and full prefix has rank fifteen")
            matching = next(perm for perm in permutations(range(6))
                            if all(matrix[i][perm[i]] for i in range(6)))
            chosen = [convol.index(a[1:]) for a in selected]
            minor_rows = next(rr for rr in combinations(range(6), 5)
                              if determinant([[matrix[i][j] for j in chosen] for i in rr]))
            minor = determinant([[matrix[i][j] for j in chosen] for i in minor_rows])
            record.update({"residual_rows": rows, "residual_columns": convol,
                           "residual_matrix": matrix, "perfect_support_matching": matching,
                           "exact_residual_rank": 5, "exact_full_prefix_rank": 15,
                           "nonzero_minor": {"rows": minor_rows, "columns": chosen, "determinant": str(minor)}})
        fixtures.append(record)
    return {"field": "Q", "fixtures": fixtures,
            "scope": "Two-generator module columns inside the actual monomial quotient K[y,z]/(y^R,y^k*z^k,z^R). The full monomial configuration of a three-variable ideal has extra powers of the tag variable and is not counted here.",
            "conclusion": "Perfect support matching does not imply normal geometric-prefix rank. This subcheck certifies a nonzero null vector for every fixture and exact rank at k=2. The subsequent graded-Jordan argument determines the uniform exact defect k^2/4."}


def jordan_filtered_convolution_controls():
    """Recover graded strings from actual multiplication ranks, then check prefixes."""
    multi_checks = 0
    for length in range(1, 5):
        nodes = (0, 1, -2)
        for cuts in product(range(length + 1), repeat=3):
            keys = [(i, j) for i in range(3) for j in range(cuts[i])]
            columns = {(i, j): tuple(Fraction(nodes[i] ** (h - j), factorial(h - j)) if h >= j else 0
                                     for h in range(length)) for i, j in keys}
            require(len(greedy(keys, columns, (1, 1), 0)) == min(length, sum(cuts)),
                    "three distinct exponential nodes obey the string prefix formula")
            multi_checks += 1
    repeated = {(0, 0): (1, 0), (1, 0): (1, 0)}
    require(len(greedy(repeated, repeated, (1, 1), 0)) == 1,
            "duplicate nodes invalidate the distinct-node rank formula")
    fixtures = []
    cases = [(k, None) for k in (2, 4, 6)] + [(None, two_variables) for two_variables in (False, True)]
    for k, two_variables in cases:
        if k is not None:
            r = 5 * k // 2
            domain = sorted((i, j) for i in range(r) for j in range(r) if i < k or j < k)
            scale, two_variables = factorial(r - 1) ** 2, True
        else:
            domain, scale = list(product(range(2), repeat=2)), 1
        top = max(map(sum, domain))

        def power_entry(row, column):
            if any(x < y for x, y in zip(row, column)):
                return 0
            degree = sum(row) - sum(column)
            return comb(degree, row[0] - column[0]) if two_variables else int(row[1] == column[1])

        ranks = {}
        for a in range(top + 1):
            source = [p for p in domain if sum(p) == a]
            for b in range(a, top + 1):
                target = [p for p in domain if sum(p) == b]
                columns = {p: tuple(power_entry(row, p) for row in target) for p in source}
                ranks[a, b] = len(greedy(source, columns, (1, 1), 0))
        strings = []
        for start in range(top + 1):
            for end in range(start, top + 1):
                number = (ranks.get((start, end), 0) - ranks.get((start - 1, end), 0)
                          - ranks.get((start, end + 1), 0) + ranks.get((start - 1, end + 1), 0))
                require(number >= 0, "graded interval multiplicities from exact map ranks")
                strings.extend([(start, end - start + 1)] * number)
        require(sum(length for _, length in strings) == len(domain), "graded Jordan strings cover the quotient")
        if k is not None:
            expected = [(i - 1, r + k + 1 - (i + 1) // 2 - i) for i in range(1, 2 * k + 1)]
            require(strings == expected, "known two-variable strong Lefschetz strings have the stated endpoints")
            cuts = (list(product(range(top + 3), repeat=2)) if k == 2 else
                    [(2 * k, 2 * k - 1), (1, 2), (2 * k, 2 * k), (r, r - 1), (0, top + 1)])
        else:
            cuts = list(product(range(top + 3), repeat=2))

        def exp_entry(row, column, node=1):
            if any(x < y for x, y in zip(row, column)) or (not two_variables and row[1] != column[1]):
                return 0
            return scale * node ** (sum(row) - sum(column)) // (
                factorial(row[0] - column[0]) * factorial(row[1] - column[1]))

        prefix_checks = []
        for a, b in cuts:
            identity = [p for p in domain if sum(p) < a]
            source = [p for p in domain if sum(p) < b]
            target = [p for p in domain if sum(p) >= a]
            columns = {p: tuple(exp_entry(row, p) for row in target) for p in source}
            actual = len(identity) + len(greedy(source, columns, (1, 1), 0))
            predicted = sum(min(length, min(length, max(0, a - start)) +
                                min(length, max(0, b - start))) for start, length in strings)
            require(actual == predicted, "graded Jordan formula agrees with direct convolution rank")
            prefix_checks.append({"identity_degree_cut": a, "convolution_degree_cut": b, "rank": actual})
            if k is not None and (a, b) == (2 * k, 2 * k - 1):
                require(actual == 4 * k * k - k * k // 4, "uniform exact sparse-prefix rank defect")
        record = {"k": k, "L": "y+z" if two_variables else "y", "quotient_length": len(domain),
                  "graded_Jordan_strings_start_and_length": strings,
                  "graded_multiplication_map_ranks_checked": len(ranks), "prefix_rank_checks": prefix_checks}
        if k == 2:
            longest = max(length for _, length in strings)
            inputs = [(i,) + p for i in range(longest) for p in domain]
            columns = {a: (tuple(int(row == a[1:]) for row in domain) if not a[0] else
                           tuple(scale ** (a[0] - 1) * exp_entry(row, a[1:], a[0]) for row in domain))
                       for a in inputs}
            x_matrix = tuple(columns[(1,) + p] for p in domain)
            for a in inputs:
                if a[0]:
                    require(matvec(x_matrix, columns[(a[0] - 1,) + a[1:]]) == columns[a],
                            "all tag powers are actual powers of multiplication by the exponential")
            costs = []
            for weight in ((1, 2, 2), (2, 1, 1), (3, 2, 2), (5, 1, 1)):
                selected = greedy(inputs, columns, weight, 0)
                require(len(selected) == len(domain), "full tag-power model has the claimed quotient rank")
                actual = sum(dot(a, weight) for a in selected)
                predicted = sum(start * length * weight[1] +
                                sum(score for score, _ in integer_prefix((weight[0], weight[1]), length))
                                for start, length in strings)
                require(actual == predicted, "full diagonal state support function is the Jordan sum of planar corner-cut costs")
                costs.append({"positive_weights": weight, "minimum_basis_cost": actual,
                              "standard_monomials": sorted(selected)})
            record["full_tag_power_diagonal_cost_controls"] = costs
        fixtures.append(record)
    return {"field": "Q", "fixtures": fixtures, "three_node_string_rank_checks": multi_checks,
            "duplicate_node_rank_failure_checked": True,
            "uniform_rank_formula": "For distinct scalar nodes tau_i and total-degree cuts a_i, rank of sum_i exp(tau_i*L)*B_{<a_i} is sum_strings min(length,sum_i clamp(a_i-start,0,length)).",
            "sparse_family_exact_defect": "For every even k>=2 the selected 4k^2-column prefix has rank 4k^2-k^2/4. This follows from the known two-variable strong Lefschetz theorem and the rank formula.",
            "scope": "This controls total-degree filtrations and the full exponential-graph algebra on diagonal weights. It does not count its full positive Groebner fan or improve the unrestricted asymptotic bounds."}


def strict_linear_infeasibility_certificate(normals):
    """Small positive dependence certifying that all dot products cannot be positive."""
    vectors = sorted(set(normals) | {(1, 0, 0), (0, 1, 0), (0, 0, 1)})
    for size in (2, 3, 4):
        for selected in combinations(vectors, size):
            if size == 2:
                if cross(*selected) != (0, 0, 0) or dot(*selected) >= 0:
                    continue
                j = next(j for j in range(3) if selected[0][j])
                coefficients = (abs(selected[1][j]), abs(selected[0][j]))
            elif size == 3:
                rows = list(zip(*selected))
                coefficients = next((cross(a, b) for a, b in combinations(rows, 2)
                                     if cross(a, b) != (0, 0, 0)), (0, 0, 0))
            else:
                coefficients = tuple((-1) ** j * dot(other[0], cross(other[1], other[2]))
                                     for j in range(4)
                                     for other in [selected[:j] + selected[j + 1:]])
            if max(coefficients) <= 0:
                coefficients = tuple(-x for x in coefficients)
            if min(coefficients) < 0 or not any(coefficients):
                continue
            if any(sum(value * normal[j] for value, normal in zip(coefficients, selected))
                   for j in range(3)):
                continue
            factor = gcd(*coefficients)
            coefficients = tuple(x // factor for x in coefficients)
            return {"normals": selected, "nonnegative_coefficients": coefficients}
    raise AssertionError("missing strict-linear-system infeasibility certificate")


def unequal_weight_jordan_countercontrols():
    """Same fixed graded algebra and L; complete fans for three coordinate embeddings."""
    domain = list(product(range(2), range(3)))
    n, ground, states = len(domain), hyperbolic_cross(len(domain), 3), finite_staircases(len(domain), 3)
    coordinates = {a: i for i, a in enumerate(domain)}
    unit = tuple(int(a == (0, 0)) for a in domain)
    shift = [tuple(tuple(int(row == tuple(value + (axis == j) for j, value in enumerate(a)))
                         for row in domain) for a in domain) for axis in range(2)]
    # The tag image is 2*(exp(y+z)-1), modulo y^2,z^3.
    polynomial = {a: 2 // (factorial(a[0]) * factorial(a[1])) for a in domain if any(a)}
    x_matrix = []
    for a in domain:
        column = [0] * n
        for b, coefficient in polynomial.items():
            row = tuple(v + w for v, w in zip(a, b))
            if row in coordinates:
                column[coordinates[row]] += coefficient
        x_matrix.append(tuple(column))
    l_matrix = tuple(tuple(y + z for y, z in zip(a, b)) for a, b in zip(*shift))
    second_generator = tuple(2 * y - z for y, z in zip(matvec(shift[0], unit), matvec(shift[1], unit)))
    jordan_basis = []
    for generator, length in ((unit, 4), (second_generator, 2)):
        column = generator
        for _ in range(length):
            jordan_basis.append(column)
            column = matvec(l_matrix, column)
        require(column == (0,) * n, "the claimed homogeneous Jordan strings terminate")
    require(determinant(jordan_basis) != 0, "the two displayed homogeneous strings form a basis")

    def inverse_columns(columns):
        matrix = [[Fraction(columns[j][i]) for j in range(n)] +
                  [Fraction(i == j) for j in range(n)] for i in range(n)]
        for j in range(n):
            pivot = next((i for i in range(j, n) if matrix[i][j]), None)
            if pivot is None:
                return None
            matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
            value = matrix[j][j]
            matrix[j] = [x / value for x in matrix[j]]
            for i in range(n):
                if i != j:
                    value = matrix[i][j]
                    matrix[i] = [x - value * y for x, y in zip(matrix[i], matrix[j])]
        return [row[n:] for row in matrix]

    results = []
    for parameter, expected in ((0, 11), (1, 9), (2, 15)):
        u_matrix = tuple(tuple(x + parameter * y for x, y in zip(a, b))
                         for a, b in zip(*shift))
        matrices = (tuple(x_matrix), u_matrix, shift[1])
        for i, j in combinations(range(3), 2):
            for a in range(n):
                require(matvec(matrices[i], matrices[j][a]) ==
                        matvec(matrices[j], matrices[i][a]), "fixed-algebra coordinate matrices commute")
        for index, a in enumerate(domain):
            recovered = unit
            for axis, exponent in enumerate(a):
                for _ in range(exponent):
                    recovered = matvec(shift[axis], recovered)
            require(recovered == tuple(int(i == index) for i in range(n)), "the reference quotient is cyclic")
            v2 = matvec(matrices[2], matrices[2][index])
            uv = matvec(matrices[1], matrices[2][index])
            uv2 = matvec(matrices[1], v2)
            formula = tuple(2 * u + 2 * (1 - parameter) * v + 2 * uv_value +
                            (1 - 2 * parameter) * v2_value + uv2_value
                            for u, v, uv_value, v2_value, uv2_value in
                            zip(matrices[1][index], matrices[2][index], uv, v2, uv2))
            require(formula == matrices[0][index], "displayed exponential graph relation holds on every column")
            require(matvec(shift[0], shift[0][index]) == (0,) * n and
                    matvec(matrices[2], v2) == (0,) * n, "displayed pure-power relations hold")
        columns = {(0, 0, 0): unit}

        def residue(a):
            if a not in columns:
                i = next(i for i, value in enumerate(a) if value)
                previous = tuple(value - (i == j) for j, value in enumerate(a))
                columns[a] = matvec(matrices[i], residue(previous))
            return columns[a]

        records, coherent, independent = [], [], 0
        for state in states:
            ordered = sorted(state)
            inverse = inverse_columns([residue(a) for a in ordered])
            minor = determinant([residue(a) for a in ordered])
            require((inverse is None) == (minor == 0), "independent inverse and determinant basis checks agree")
            record = {"staircase": ordered, "basis_determinant": str(minor)}
            if inverse is None:
                record["status"] = "dependent"
                records.append(record)
                continue
            independent += 1
            border = {tuple(value + (i == j) for j, value in enumerate(a))
                      for a in state for i in range(3)} - state
            normals = set()
            for b in sorted(border):
                vector = residue(b)
                coefficients = [dot(row, vector) for row in inverse]
                require(tuple(sum(coefficients[j] * residue(ordered[j])[i] for j in range(n))
                              for i in range(n)) == vector, "exact border normal form reconstructs the residue")
                for a, coefficient in zip(ordered, coefficients):
                    if coefficient:
                        normal = tuple(x - y for x, y in zip(b, a))
                        require(max(normal) > 0, "a border point cannot divide an included staircase point")
                        if min(normal) < 0:
                            factor = gcd(*normal)
                            normals.add(tuple(x // factor for x in normal))
            polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
            for normal in sorted(normals):
                polygon = clip(polygon, normal, 1)
                if not polygon:
                    break
            record["strict_border_comparison_normals"] = sorted(normals)
            if not polygon:
                record["status"] = "independent_but_incoherent"
                record["infeasibility_certificate"] = strict_linear_infeasibility_certificate(normals)
            else:
                weight = tuple(sum(p[i] for p in polygon) / len(polygon) for i in range(3))
                denominator = lcm(*(x.denominator for x in weight))
                weight = tuple(int(x * denominator) for x in weight)
                # Break variable-weight ties without crossing any strict border constraint.
                delta = (1, 2, 3)
                scale = 4 + max((abs(dot(delta, normal)) for normal in normals), default=0)
                weight = tuple(scale * x + y for x, y in zip(weight, delta))
                factor = gcd(*weight)
                weight = tuple(x // factor for x in weight)
                require(min(weight) > 0 and len(set(weight)) == 3 and
                        all(dot(weight, normal) > 0 for normal in normals), "strict border and variable-order witness")
                require(greedy(ground, {a: residue(a) for a in ground}, weight, 0) == state,
                        "full-ground Gaussian selection independently checks the coherent staircase")
                record.update({"status": "coherent", "positive_integer_weights": weight,
                               "variable_order_high_to_low": sorted(range(3), key=lambda i: -weight[i])})
                coherent.append(record)
            records.append(record)
        require(len(states) == 48 and len(coherent) == expected, "exhaustive six-box staircase fan counts")
        costs = []
        for weight in ((1, 2, 2), (2, 1, 1), (2, 1, 3)):
            selected = greedy(ground, {a: residue(a) for a in ground}, weight, 0)
            cost = sum(dot(a, weight) for a in selected)
            expected_cost = (10 if weight == (1, 2, 2) else 8 if weight == (2, 1, 1)
                             else {0: 13, 1: 12, 2: 11}[parameter])
            require(len(selected) == n and cost == expected_cost, "diagonal equality and unequal-weight cost separation")
            ranks = []
            for threshold in range(5):
                prefix = [a for a in ground if dot(a, weight) <= threshold]
                ranks.append(len(greedy(prefix, {a: residue(a) for a in prefix}, weight, 0)))
            costs.append({"positive_weights": weight, "minimum_basis_cost": cost,
                          "standard_monomials": sorted(selected), "ranks_at_cost_thresholds_zero_through_four": ranks})
        groups = {}
        for record in coherent:
            groups.setdefault(tuple(record["variable_order_high_to_low"]), []).append(record)
        pair = next(group[:2] for group in groups.values() if len(group) >= 2)
        results.append({"lambda": parameter, "coordinate_images": "(x,u,v)=(2*(exp(y+z)-1),y+lambda*z,z)",
                        "multiplication_matrices_by_columns": matrices, "reference_basis_y_z": domain,
                        "exact_initial_ideal_count": len(coherent), "independent_staircase_count": independent,
                        "all_staircase_certificates": records, "cost_controls": costs,
                        "same_strict_variable_order_distinct_initials": [
                            {key: record[key] for key in ("staircase", "positive_integer_weights", "variable_order_high_to_low")}
                            for record in pair]})
    return {"field": "Q", "fixed_algebra": "Q[y,z]/(y^2,z^3)", "fixed_L": "y+z",
            "Hilbert_function": [1, 2, 2, 1], "graded_Jordan_strings_start_and_length": [[0, 4], [1, 2]],
            "graded_Jordan_basis_columns": jordan_basis,
            "embeddings": results,
            "conclusion": "Even the same graded algebra and the same L and tag element have different full fans under different degree-one generating coordinates. Graded Jordan data determine the proved diagonal support function but not unequal-weight costs or the full fan count.",
            "scope": "Exact falsification certificates for a transfer lemma, not an asymptotic lower bound or a novelty claim."}


def codimension_one_binomial_quotient_controls():
    """Complete small fans for a one-dimensional quotient of a binomial algebra."""
    fixtures = []
    for period, delay, expected in ((7, 0, (11, 18)), (6, 1, (8, 12)), (6, 2, (6, 8))):
        q, powers = period + delay, (1, 2, 3)
        parent_modulus = [0] * (q + 1)
        parent_modulus[delay], parent_modulus[q] = -1, 1
        child_modulus = [0] * delay + [1] * period
        # Removing one factor t-1 gives an actual ideal containment of length difference one.
        product_modulus = [0] * (q + 1)
        for i, value in enumerate(child_modulus):
            product_modulus[i] -= value
            product_modulus[i + 1] += value
        require(product_modulus == parent_modulus, "parent=(t-1)*child")
        circuit = set(range(delay, q))
        full_ground = hyperbolic_cross(q, 3)

        def residue_class(a):
            degree = dot(a, powers)
            return degree if degree < delay else delay + (degree - delay) % period

        def model(modulus, ground):
            t = (0, 1) + (0,) * (len(modulus) - 3)
            images = [univariate_power_mod(t, exponent, modulus) for exponent in powers]
            return cyclic_univariate_model(modulus, images, ground)

        parent_columns = model(parent_modulus, full_ground)
        child_columns = model(child_modulus, full_ground)
        require(all(parent_columns[a] == tuple(int(i == residue_class(a)) for i in range(q))
                    for a in full_ground), "parent columns are exactly the stated parallel classes")
        circuit_columns = [child_columns[(i, 0, 0)] for i in sorted(circuit)]
        require(tuple(sum(column[i] for column in circuit_columns) for i in range(q - 1)) ==
                (0,) * (q - 1), "the unique child circuit has all coefficients one")
        require(all(determinant([child_columns[(i, 0, 0)] for i in range(q) if i != omitted]) != 0
                    for omitted in circuit), "omitting any circuit class gives a child basis")
        if delay:
            require(determinant([child_columns[(i, 0, 0)] for i in range(q) if i == 0 or i >= 2]) == 0
                    if delay == 2 else determinant(circuit_columns) == 0,
                    "omitting a coloop does not give a child basis")

        fans, checked = [], 0
        for n, columns, removed in ((q, parent_columns, 0), (q - 1, child_columns, 1)):
            records = []
            for state in finite_staircases(n, 3):
                ordered = sorted(state)
                classes = {residue_class(a): a for a in ordered}
                missing = set(range(q)) - set(classes)
                independent = len(classes) == n and (not removed or missing <= circuit)
                minor = determinant([columns[a] for a in ordered])
                require((minor != 0) == independent, "single-circuit basis criterion agrees with exact determinant")
                require(minor in (-1, 0, 1), "these fixture basis decisions hold in every characteristic")
                record = {"staircase": ordered, "basis_determinant": str(minor)}
                checked += 1
                if not independent:
                    record["status"] = "dependent"
                    records.append(record)
                    continue
                normals = set()
                for b in full_ground:
                    if b in state:
                        continue
                    label = residue_class(b)
                    targets = ([classes[label]] if label in classes else
                               [a for key, a in classes.items() if key in circuit])
                    for a in targets:
                        normal = tuple(x - y for x, y in zip(b, a))
                        if min(normal) >= 0:
                            continue
                        divisor = gcd(*normal)
                        normals.add(tuple(x // divisor for x in normal))
                polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
                for normal in sorted(normals):
                    polygon = clip(polygon, normal, 1)
                    if not polygon:
                        break
                record["strict_class_and_circuit_comparison_normals"] = sorted(normals)
                if not polygon:
                    certificate = strict_linear_infeasibility_certificate(normals)
                    require(all(value >= 0 for value in certificate["nonnegative_coefficients"]) and
                            any(certificate["nonnegative_coefficients"]) and
                            tuple(sum(value * normal[i] for value, normal in
                                      zip(certificate["nonnegative_coefficients"], certificate["normals"]))
                                  for i in range(3)) == (0, 0, 0), "independent infeasibility certificate check")
                    record.update({"status": "independent_but_incoherent", "infeasibility_certificate": certificate})
                else:
                    weight = tuple(sum(p[i] for p in polygon) / len(polygon) for i in range(3))
                    denominator = lcm(*(value.denominator for value in weight))
                    weight = tuple(int(value * denominator) for value in weight)
                    delta = (1, q + 1, (q + 1) ** 2)
                    scale = 1 + max(dot(delta, a) for a in full_ground)
                    weight = tuple(scale * value + step for value, step in zip(weight, delta))
                    divisor = gcd(*weight)
                    weight = tuple(value // divisor for value in weight)
                    require(min(weight) > 0 and all(dot(weight, normal) > 0 for normal in normals) and
                            len({dot(weight, a) for a in full_ground}) == len(full_ground), "strict fully separated witness")
                    parent_basis = greedy(full_ground, parent_columns, weight, 0)
                    omitted = max((a for a in parent_basis if residue_class(a) in circuit),
                                  key=lambda a: dot(weight, a))
                    formula = parent_basis - {omitted} if removed else parent_basis
                    require(formula == state, "parent class minima minus largest circuit representative")
                    for characteristic in (0, 2, 3, 7):
                        require(greedy(full_ground, columns, weight, characteristic) == state,
                                "independent Gaussian selection in each recorded characteristic")
                    record.update({"status": "coherent", "positive_integer_weights": weight,
                                   "parent_standard_monomials": sorted(parent_basis),
                                   "omitted_circuit_representative": omitted if removed else None})
                records.append(record)
            fans.append({"length": n, "modulus": parent_modulus if not removed else child_modulus,
                         "initial_ideal_count": sum(record["status"] == "coherent" for record in records),
                         "all_staircase_certificates": records})
        counts = tuple(fan["initial_ideal_count"] for fan in fans)
        require(counts == expected, "complete parent and one-dimensional quotient fan counts")
        groups = {}
        for record in fans[1]["all_staircase_certificates"]:
            if record["status"] == "coherent":
                groups.setdefault(tuple(map(tuple, record["parent_standard_monomials"])), []).append(record)
        pair = next((group[:2] for group in groups.values() if len(group) >= 2), [])
        if not delay:
            require(pair, "two child initials have exactly the same parent initial")
        fixtures.append({"period": period, "nilpotent_delay": delay, "parent_length": q,
                         "coordinate_powers": powers, "field_characteristics_checked": [0, 2, 3, 7],
                         "reference_parent_basis": ["t^" + str(i) for i in range(q)],
                         "single_circuit_classes": sorted(circuit), "coloop_classes": list(range(delay)),
                         "nilpotent_parent_and_child_at_zero": delay >= 2,
                         "monomial_rank_and_staircase_checks": checked, "fans": fans,
                         "same_parent_distinct_child_initials": [
                             {key: record[key] for key in ("staircase", "positive_integer_weights", "parent_standard_monomials",
                                                          "omitted_circuit_representative")} for record in pair]})
    return {"uniform_scope": "Every ideal I containing a colength-q binomial ideal J, with dim(I/J)=1, over every field, in dimensions two and three.",
            "count_bound": "O(M*alpha(M)) in dimension two and O(M^2*alpha(M)) in dimension three, M=|V_q^d|. Thus O(q^2*log(q+1)^4*alpha(q*(1+log(q))^2)) in dimension three.",
            "geometric_input": "Pach--Sharir (1989), Theorem 1-prime, the segment-envelope bound and the overlay argument of Lemma 2.2.",
            "decision": "Deleting one point from a cyclic binomial orbit cannot produce an n^(2+epsilon) lower family in three variables for fixed epsilon>0. The unrestricted goal is not completed; logarithmic improvements are not ruled out.",
            "priority": "A consequence of established binomial normal forms and envelope bounds; novelty has not been established.",
            "fixtures": fixtures}


def uniform_deletion_and_dual_relation_controls():
    """Prime Fourier certificates and a nonuniform two-relation countercontrol."""
    characteristic, root, q, powers = 1000000009, 440359144, 7, (1, 2, 3)
    require(all(characteristic % divisor for divisor in range(2, isqrt(characteristic) + 1)),
            "the reduction modulus is prime")
    require(pow(root, q, characteristic) == 1 and
            all(pow(root, i, characteristic) != 1 for i in range(1, q)) and
            sum(pow(root, i, characteristic) for i in range(q)) % characteristic == 0,
            "a primitive seventh root defines a reduction of the cyclotomic integer ring")

    def modular_determinant(matrix):
        a = [[value % characteristic for value in row] for row in matrix]
        result = 1
        for j in range(len(a)):
            pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
            if pivot is None:
                return 0
            if pivot != j:
                a[j], a[pivot] = a[pivot], a[j]
                result = -result
            value = a[j][j]
            result = result * value % characteristic
            inverse = pow(value, characteristic - 2, characteristic)
            for i in range(j + 1, len(a)):
                scale = a[i][j] * inverse % characteristic
                for k in range(j, len(a)):
                    a[i][k] = (a[i][k] - scale * a[j][k]) % characteristic
        return result % characteristic

    ground = hyperbolic_cross(q, 3)
    classes = {j: [a for a in ground if dot(a, powers) % q == j] for j in range(q)}
    fixtures = []
    for n, expected, second_rows in ((4, 10, (0, 2, 4, 6)), (5, 14, (0, 2, 3, 5, 6))):
        minor_product, minor_count = 1, 0
        choices = list(combinations(range(q), n))
        for rows in choices:
            for columns in choices:
                value = modular_determinant([[pow(root, i * j, characteristic) for j in columns] for i in rows])
                require(value != 0, "every relevant Fourier minor is nonzero in the certified reduction")
                minor_product = minor_product * value % characteristic
                minor_count += 1
        require(minor_count == comb(q, n) ** 2 and minor_product != 0, "compact simultaneous minor certificate")
        row_sets = (tuple(range(n)), second_rows)
        evaluations = [{a: tuple(pow(root, row * dot(a, powers), characteristic) for row in rows)
                        for a in ground} for rows in row_sets]
        records = []
        for state in finite_staircases(n, 3):
            ordered = sorted(state)
            representatives = {dot(a, powers) % q: a for a in ordered}
            independent = len(representatives) == n
            minors = [modular_determinant([columns[a] for a in ordered]) for columns in evaluations]
            require(all((value != 0) == independent for value in minors), "staircase ranks agree in both point sets")
            record = {"staircase": ordered, "basis_determinants_mod_prime": minors}
            if not independent:
                record["status"] = "dependent"
                records.append(record)
                continue
            normals = set()
            for b in ground:
                if b in state:
                    continue
                label = dot(b, powers) % q
                targets = [representatives[label]] if label in representatives else state
                for a in targets:
                    normal = tuple(x - y for x, y in zip(b, a))
                    if min(normal) < 0:
                        divisor = gcd(*normal)
                        normals.add(tuple(x // divisor for x in normal))
            polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
            for normal in sorted(normals):
                polygon = clip(polygon, normal, 1)
                if not polygon:
                    break
            record["strict_truncated_class_comparison_normals"] = sorted(normals)
            if not polygon:
                certificate = strict_linear_infeasibility_certificate(normals)
                require(all(value >= 0 for value in certificate["nonnegative_coefficients"]) and
                        any(certificate["nonnegative_coefficients"]) and
                        tuple(sum(value * normal[i] for value, normal in
                                  zip(certificate["nonnegative_coefficients"], certificate["normals"]))
                              for i in range(3)) == (0, 0, 0), "independent positive-dependence verification")
                record.update({"status": "independent_but_incoherent", "infeasibility_certificate": certificate})
            else:
                weight = tuple(sum(p[i] for p in polygon) / len(polygon) for i in range(3))
                denominator = lcm(*(value.denominator for value in weight))
                weight = tuple(int(value * denominator) for value in weight)
                delta = (1, q + 1, (q + 1) ** 2)
                scale = 1 + max(dot(delta, a) for a in ground)
                weight = tuple(scale * value + step for value, step in zip(weight, delta))
                divisor = gcd(*weight)
                weight = tuple(value // divisor for value in weight)
                require(min(weight) > 0 and all(dot(weight, normal) > 0 for normal in normals) and
                        len({dot(weight, a) for a in ground}) == len(ground), "fully separated strict witness")
                minima = {j: min(group, key=lambda a: dot(weight, a)) for j, group in classes.items()}
                removed = sorted(range(q), key=lambda j: -dot(weight, minima[j]))[:q - n]
                formula = frozenset(a for j, a in minima.items() if j not in removed)
                require(formula == state, "uniform dual rule omits precisely the largest r class minima")
                for columns in evaluations:
                    require(greedy(ground, columns, weight, characteristic) == state,
                            "independent Gaussian selection agrees in consecutive and nonconsecutive point sets")
                record.update({"status": "coherent", "positive_integer_weights": weight,
                               "parent_class_minima": minima, "omitted_classes": removed,
                               "minimum_child_cost": sum(dot(weight, a) for a in state)})
            records.append(record)
        require(sum(record["status"] == "coherent" for record in records) == expected,
                "complete bounded fan counts for deletion defects two and three")
        fixtures.append({"orbit_size": q, "n": n, "deletion_defect": q - n,
                         "coordinate_powers": powers, "parameter_index_subsets": row_sets,
                         "all_row_and_column_minors_checked": minor_count,
                         "product_of_all_minors_mod_prime": minor_product,
                         "initial_ideal_count": expected, "all_staircase_certificates": records})

    # This quotient has two circuits, so its omitted classes need not be the two most expensive.
    full_ground = hyperbolic_cross(6, 3)
    parent_modulus, child_modulus = [-1, 0, 0, 0, 0, 0, 1], [1, 0, 1, 0, 1]
    models = []
    for modulus in (parent_modulus, child_modulus):
        t = (0, 1) + (0,) * (len(modulus) - 3)
        images = [univariate_power_mod(t, exponent, modulus) for exponent in powers]
        models.append(cyclic_univariate_model(modulus, images, full_ground))
    parent, child = models
    relation_rows = [tuple(int(i % 2 == parity) for i in range(6)) for parity in range(2)]
    for row in relation_rows:
        require(tuple(sum(row[j] * child[(j, 0, 0)][i] for j in range(6)) for i in range(4)) == (0,) * 4,
                "the two displayed relations vanish in the actual quotient")
    dual_checks = []
    for omitted in combinations(range(6), 2):
        kept = [j for j in range(6) if j not in omitted]
        actual = determinant([child[(j, 0, 0)] for j in kept])
        dual = determinant([[row[j] for j in omitted] for row in relation_rows])
        require((actual != 0) == (dual != 0), "complementary basis criterion for the dual relation matrix")
        dual_checks.append({"omitted_classes": omitted, "kept_classes": kept,
                            "child_basis_determinant": str(actual), "dual_basis_determinant": str(dual)})
    weight = (10, 1, 3)
    parent_basis, child_basis = [greedy(full_ground, model, weight, 0) for model in models]
    minima = {dot(a, powers) % 6: a for a in parent_basis}
    costs = {j: dot(weight, a) for j, a in minima.items()}
    admissible = [set(record["omitted_classes"]) for record in dual_checks if record["dual_basis_determinant"] != "0"]
    optimal = max(admissible, key=lambda omitted: sum(costs[j] for j in omitted))
    wrong = set(sorted(costs, key=lambda j: -costs[j])[:2])
    correct_basis = frozenset(a for j, a in minima.items() if j not in optimal)
    wrong_basis = frozenset(a for j, a in minima.items() if j not in wrong)
    require(optimal == {1, 4} and wrong == {1, 5} and correct_basis == child_basis != wrong_basis,
            "the dual matroid constraint changes the omission")
    require(determinant([child[a] for a in sorted(wrong_basis)]) == 0 and
            determinant([child[a] for a in sorted(child_basis)]) != 0,
            "the unconstrained two-largest rule produces an actually dependent staircase")
    counter = {"field": "Q", "parent_modulus": parent_modulus, "child_modulus": child_modulus,
               "kernel_description": "The six-root orbit with t=1 and t=-1 removed; (x,y,z)=(t,t^2,t^3).",
               "relation_rows_in_reference_classes": relation_rows, "complementary_basis_checks": dual_checks,
               "positive_weights": weight, "parent_class_minima": minima, "parent_class_costs": costs,
               "correct_omitted_classes": sorted(optimal), "incorrect_unconstrained_omitted_classes": sorted(wrong),
               "actual_child_standard_monomials": sorted(child_basis), "incorrect_dependent_staircase": sorted(wrong_basis),
               "dependent_circuit": {"classes": [0, 2, 4], "coefficients": [1, 1, 1]},
               "conclusion": "A general quotient of a binomial algebra requires dual-matroid basis optimization. The uniform truncation rule and level bound cannot be asserted for it without a separate hypothesis."}
    return {"uniform_geometric_bound": "For q nonempty parallel classes on M points and rank n=q-r uniform truncation, the number of positive coherent bases is O(M*(r+1)*alpha(M)) in dimension two and O(M^2*(r+1)*alpha(M)) in dimension three.",
            "algebraic_specialization": "Prime cyclic orbit subsets over C, and consecutive samples for any cyclic orbit size. With M<=q*(1+log(q))^2, the spatial bound is O(q^2*(r+1)*(1+log(q))^4*alpha(q*(1+log(q))^2)).",
            "reduction_prime": characteristic, "primitive_seventh_root_mod_prime": root,
            "cyclotomic_certificate_scope": "A nonzero product of the indicated Fourier minors after reduction certifies their simultaneous nonvanishing over C. No implication from complex nonvanishing to all finite characteristics is used.",
            "uniform_fixtures": fixtures, "nonuniform_two_relation_countercontrol": counter,
            "priority": "Uses established envelope bounds, elementary random-sampling estimates and standard matroid duality. Novelty and an unrestricted asymptotic advance are not established."}


def box_residual_duality_controls(codimension_one):
    """Box-initial colon duality, with nonreduced and nonbox negative controls."""
    domain = list(product(range(2), repeat=3))
    coordinates = {a: i for i, a in enumerate(domain)}
    top = (1, 1, 1)
    normals = sorted({primitive(tuple(x - y for x, y in zip(a, b)))
                      for a, b in combinations(domain, 2)
                      if min(x - y for x, y in zip(a, b)) < 0 < max(x - y for x, y in zip(a, b))})
    weights = positive_chambers(normals)
    fixtures = []
    specifications = [
        ("local_cube", ((0, 0, 1),) * 3, 0),
        ("mixed_repeated_grid", ((-1, 0, 1), (0, -1, 1), (0, -1, 1)), -1),
        ("reduced_binary_grid", ((0, -1, 1),) * 3, -1),
    ]
    for label, moduli, constant in specifications:
        powers = [[univariate_power_mod((0, 1), exponent, modulus) for exponent in range(3)]
                  for modulus in moduli]
        products = {(a, b): tuple(powers[0][a[0] + b[0]][c[0]] *
                                 powers[1][a[1] + b[1]][c[1]] *
                                 powers[2][a[2] + b[2]][c[2]] for c in domain)
                    for a in domain for b in domain}
        for characteristic in (0, 2, 3, 5):
            def scalar(value):
                return value % characteristic if characteristic else Fraction(value)

            def row_reduce(rows):
                matrix = [[scalar(x) for x in row] for row in rows]
                rank, pivots = 0, []
                for column in range(len(domain)):
                    pivot = next((i for i in range(rank, len(matrix)) if matrix[i][column]), None)
                    if pivot is None:
                        continue
                    matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
                    inverse = (pow(matrix[rank][column], -1, characteristic) if characteristic
                               else 1 / matrix[rank][column])
                    matrix[rank] = [scalar(x * inverse) for x in matrix[rank]]
                    for i in range(len(matrix)):
                        if i != rank:
                            coefficient = matrix[i][column]
                            matrix[i] = [scalar(x - coefficient * y) for x, y in zip(matrix[i], matrix[rank])]
                    pivots.append(column)
                    rank += 1
                    if rank == len(matrix):
                        break
                return matrix[:rank], pivots

            def reduce_vector(vector, rows, pivots):
                vector = [scalar(x) for x in vector]
                for row, pivot in zip(rows, pivots):
                    coefficient = vector[pivot]
                    vector = [scalar(x - coefficient * y) for x, y in zip(vector, row)]
                return tuple(vector)

            def nullspace(rows):
                reduced, pivots = row_reduce(rows)
                result = []
                for free in range(len(domain)):
                    if free not in pivots:
                        vector = [scalar(int(i == free)) for i in range(len(domain))]
                        for row, pivot in zip(reduced, pivots):
                            vector[pivot] = scalar(-row[free])
                        result.append(tuple(vector))
                return result

            def multiply(left, right):
                return tuple(scalar(sum(left[i] * right[j] * products[a, b][k]
                                        for i, a in enumerate(domain) for j, b in enumerate(domain)))
                             for k in range(len(domain)))

            unit_vectors = [tuple(scalar(int(i == j)) for i in range(len(domain)))
                            for j in range(len(domain))]
            variable_vectors = [unit_vectors[coordinates[e]] for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
            require(all(multiply(a, multiply(b, unit)) == multiply(b, multiply(a, unit))
                        for a, b in combinations(variable_vectors, 2) for unit in unit_vectors),
                    "the three actual parent multiplication matrices commute")
            generator = tuple(scalar(constant if a == (0, 0, 0) else int(sum(a) == 1)) for a in domain)
            multiplication_columns = [multiply(unit, generator) for unit in unit_vectors]
            relations, relation_pivots = row_reduce(multiplication_columns)
            annihilator, annihilator_pivots = row_reduce(nullspace(list(zip(*multiplication_columns))))
            require(len(relations) + len(annihilator) == len(domain), "rank-nullity in the actual grid algebra")
            require(all(not any(multiply(a, b)) for a in relations for b in annihilator),
                    "the computed colon annihilates the entire relation ideal")
            gram = [[scalar(products[a, b][coordinates[top]]) for b in domain] for a in domain]
            require(len(row_reduce(gram)[0]) == len(domain), "top-coefficient pairing is perfect in every checked field")
            orthogonal = nullspace([[scalar(sum(row[i] * gram[i][j] for i in range(len(domain))))
                                     for j in range(len(domain))] for row in relations])
            require(len(orthogonal) == len(annihilator) and
                    all(not any(reduce_vector(row, annihilator, annihilator_pivots)) for row in orthogonal),
                    "independent pairing-orthogonal and multiplication-kernel computations agree")
            for rows, pivots in ((relations, relation_pivots), (annihilator, annihilator_pivots)):
                require(all(not any(reduce_vector(multiply(unit_vectors[coordinates[e]], row), rows, pivots))
                            for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1)) for row in rows),
                        "both relation subspaces are stable under all three multiplication actions")
            quotient_columns = [{a: reduce_vector(unit, rows, pivots) for a, unit in zip(domain, unit_vectors)}
                                for rows, pivots in ((relations, relation_pivots), (annihilator, annihilator_pivots))]
            records = []
            for weight in weights:
                leading = [leading_row_pivots(rows, domain, weight, characteristic)
                           for rows in (relations, annihilator)]
                bases = [set(domain) - chosen for chosen in leading]
                reflected = {tuple(c - a for c, a in zip(top, point)) for point in leading[0]}
                require(bases[1] == reflected, "exact reflected complement for the colon standard basis")
                monomial_annihilator = {b for b in domain
                                        if all(any(x + y >= 2 for x, y in zip(a, b)) for a in leading[0])}
                require(leading[1] == monomial_annihilator, "initial colon equals colon of the two initial ideals")
                for basis, columns in zip(bases, quotient_columns):
                    require(greedy(domain, columns, weight, characteristic) == basis,
                            "ascending quotient-column selection independently checks descending relation pivots")
                records.append({"positive_integer_weights": weight, "child_standard_monomials": sorted(bases[0]),
                                "residual_standard_monomials": sorted(bases[1])})
            child_states = {tuple(record["child_standard_monomials"]) for record in records}
            residual_states = {tuple(record["residual_standard_monomials"]) for record in records}
            require(len(child_states) == len(residual_states), "pointwise bijection identifies the complete two fans")
            unit = unit_vectors[coordinates[(0, 0, 0)]]
            sum_rows, sum_pivots = row_reduce(relations + annihilator)
            comaximal = not any(reduce_vector(unit, sum_rows, sum_pivots))
            common_top = (not any(reduce_vector(unit_vectors[coordinates[top]], relations, relation_pivots)) and
                          not any(reduce_vector(unit_vectors[coordinates[top]], annihilator, annihilator_pivots)))
            if label == "local_cube":
                require(not comaximal and common_top, "colon links need not be comaximal CRT complements")
                require(all(any(reduce_vector(a, relations, relation_pivots)) for a in variable_vectors),
                        "the local child has nonzero square-zero coordinates, so this is its maximal grid")
            fixtures.append({"label": label, "univariate_moduli": moduli, "linear_generator_constant": constant,
                             "linear_generator_variable_coefficients": [1, 1, 1], "characteristic": characteristic,
                             "parent_length": len(domain), "child_length": len(annihilator),
                             "residual_length": len(relations), "number_of_initial_ideals": len(child_states),
                             "relation_rows": [[str(x) for x in row] for row in relations],
                             "annihilator_rows": [[str(x) for x in row] for row in annihilator],
                             "comaximal": comaximal, "nonzero_top_in_both_ideals": common_top,
                             "all_comparison_chamber_checks": records})

    nongrid_ground = hyperbolic_cross(8, 3)
    nongrid_models = [{a: tuple(int((a[0] + 2 * a[1], a[2]) == (i, j))
                               for i, j in product(range(length), range(2))) for a in nongrid_ground}
                      for length in (4, 3)]
    nongrid_boxes = []
    for weight in ((17, 7, 3), (5, 13, 3)):
        box = greedy(nongrid_ground, nongrid_models[0], weight, 0)
        sides = tuple(max(a[i] for a in box) + 1 for i in range(3))
        require(box == set(product(*(range(side) for side in sides))), "the nongrid parent cone has a box initial")
        child = {(0, 0, 0), (0, 0, 1)}
        residual = greedy(nongrid_ground, nongrid_models[1], weight, 0)
        corner = tuple(side - 1 for side in sides)
        require(residual == {tuple(c - a for c, a in zip(corner, point)) for point in box - child},
                "box-cone residual duality also holds for the nongrid complete intersection")
        nongrid_boxes.append({"positive_integer_weights": weight, "box_sides": sides,
                              "parent_standard_monomials": sorted(box), "child_standard_monomials": sorted(child),
                              "residual_standard_monomials": sorted(residual)})

    pair = [dict(record, positive_integer_weights=weight) for record, weight in
            zip(codimension_one["fixtures"][0]["same_parent_distinct_child_initials"], ((41, 53, 10), (23, 53, 10)))]
    require(len(pair) == 2, "existing nonbox complete-intersection parent countercontrol is present")
    parent_basis = set(map(tuple, pair[0]["parent_standard_monomials"]))
    require(parent_basis == set(map(tuple, pair[1]["parent_standard_monomials"])), "same parent initial ideal")
    socle = {a for a in parent_basis if all(tuple(x + (i == j) for j, x in enumerate(a)) not in parent_basis
                                           for i in range(3))}
    require(socle == {(0, 0, 4), (1, 0, 1)}, "the parent initial algebra has two socle monomials")
    full_ground = hyperbolic_cross(7, 3)
    cyclic_models = []
    for modulus in ((-1, 0, 0, 0, 0, 0, 0, 1), (1,) * 7):
        t = (0, 1) + (0,) * (len(modulus) - 3)
        images = [univariate_power_mod(t, power, modulus) for power in (1, 2, 3)]
        cyclic_models.append(cyclic_univariate_model(modulus, images, full_ground))
    for witness in pair:
        child = set(map(tuple, witness["staircase"]))
        for characteristic in (0, 2, 3, 7):
            require(greedy(full_ground, cyclic_models[0], witness["positive_integer_weights"], characteristic) == parent_basis and
                    greedy(full_ground, cyclic_models[1], witness["positive_integer_weights"], characteristic) == child,
                    "independent actual quotient checks at the smaller nonbox counterweights")
        removed = parent_basis - child
        require(len(removed) == 1 and removed <= socle, "different omitted socle monomials")
        monomial_annihilator = {b for b in parent_basis
                                if all(tuple(x + y for x, y in zip(a, b)) not in parent_basis for a in removed)}
        require(monomial_annihilator == parent_basis - {(0, 0, 0)}, "both monomial colons have initial ideal (x,y,z)")
    return {"uniform_statement": "For every J contained in I, on every term-order cone where in(J)=(x_1^s_1,...,x_d^s_d), in(J:I)=in(J):in(I), with reflected complementary box standard bases. This holds over every field without radicality or comaximality.",
            "grid_consequence": "Every univariate grid parent has the same box at all global orders, so I and J:I have identical positive initial-ideal fans and complementary lengths.",
            "scope_limit": "A complete-intersection parent alone is insufficient: its initial ideal must be a box for this invertible recovery. General grid covers can have residual length n^d-n and do not improve the unrestricted bound.",
            "fixtures": fixtures,
            "nongrid_box_cones": {"parent_generators": ["x^2-y", "y^2", "z^2"],
                                   "child_additional_generator": "x", "residual_additional_generator": "x*y",
                                   "checks": nongrid_boxes},
            "nonbox_countercontrol": {"parent_generators": ["y-x^2", "z-x^3", "x^7-1"],
                                       "child_generators": ["y-x^2", "z-x^3", "1+x+x^2+x^3+x^4+x^5+x^6"],
                                       "actual_residual_generators": ["x-1", "y-1", "z-1"],
                                       "parent_initial_generators": ["y", "x^2", "x*z^2", "z^5"],
                                       "socle_monomials": sorted(socle),
                                       "residual_initial_generators": ["x", "y", "z"],
                                       "distinct_child_witnesses": pair,
                                       "information_loss": "Both colons are identical, but coloning again adds both parent socle monomials rather than recovering either of the two different child initials."},
            "priority": "Standard complete-intersection duality plus the initial-colon inclusion; the grid-complement fan correspondence is established prior art. No new unrestricted theorem or novelty certification."}


def nonuniform_two_level_exchange_control():
    """A realizable defect-two vertex with two exchange heights, not one cutoff."""
    groups = (((3, 0, 0), (0, 3, 0), (1, 0, 1)),
              ((0, 0, 4), (0, 1, 3), (2, 1, 0)))
    corners = groups[0] + groups[1]
    domain = sorted({a for c in corners for a in product(*(range(value + 1) for value in c))})
    require(len(domain) == 17 and all(tuple(value + (i == j) for j, value in enumerate(a)) not in domain
                                    for a in corners for i in range(3)), "six genuine parent socle monomials")
    removed_reference = {group[0] for group in groups}
    reference_bases = [domain, [a for a in domain if a not in removed_reference]]
    column_models = []
    ground = hyperbolic_cross(len(domain), 3)
    for basis in reference_bases:
        zero = (0,) * len(basis)

        def column(a):
            if a not in domain:
                return zero
            if a in basis:
                return tuple(int(a == b) for b in basis)
            group = next(group for group in groups if a == group[0])
            return tuple(-int(b in group[1:]) for b in basis)

        matrices = [[column(tuple(value + (i == j) for j, value in enumerate(a))) for a in basis]
                    for i in range(3)]
        require(all(matvec(matrices[i], vector) == matvec(matrices[j], other)
                    for i, j in combinations(range(3), 2)
                    for vector, other in zip(matrices[j], matrices[i])), "actual two-circuit quotient matrices commute over Z")
        columns = {}
        for a in ground:
            vector = column((0, 0, 0))
            for axis in range(3):
                for _ in range(a[axis]):
                    vector = matvec(matrices[axis], vector)
            require(vector == column(a), "cyclic evaluation agrees with the independently specified socle quotient")
            columns[a] = vector
        column_models.append(columns)
    rows = [tuple(int(a in group) for a in domain) for group in groups]
    normals = sorted({primitive(tuple(x - y for x, y in zip(a, b)))
                      for group in groups for a, b in combinations(group, 2)})
    weights = positive_chambers(normals)
    records, states = [], set()
    for weight in weights:
        omitted = {max(group, key=lambda a: dot(weight, a)) for group in groups}
        basis = set(domain) - omitted
        for characteristic in (0, 2, 3, 7):
            require(leading_row_pivots(rows, domain, weight, characteristic) == omitted and
                    greedy(ground, column_models[1], weight, characteristic) == basis,
                    "independent row and actual quotient-column checks on every chamber")
        states.add(tuple(sorted(basis)))
        records.append({"positive_integer_weights": weight, "omitted_monomials": sorted(omitted),
                        "standard_monomials": sorted(basis)})
    central = (1, 1, 1)
    require([dot(central, a) for a in groups[0]] == [3, 3, 2] and
            [dot(central, a) for a in groups[1]] == [4, 4, 3], "the two simultaneous exchange heights differ")
    require(cross((3, -3, 0), (0, -1, 1)) != (0, 0, 0), "the exchange lines are independent")
    central_records = []
    for weight in ((101, 100, 102), (100, 101, 102), (101, 100, 99), (100, 101, 99)):
        omitted = {max(group, key=lambda a: dot(weight, a)) for group in groups}
        basis = set(domain) - omitted
        require(tuple(sorted(basis)) in states and
                greedy(ground, column_models[1], weight, 0) == basis, "all four local sectors are realized")
        central_records.append({"positive_integer_weights": weight, "omitted_monomials": sorted(omitted)})
    require(len({tuple(record["omitted_monomials"]) for record in central_records}) == 4,
            "a genuine four-sector vertex at two distinct costs")
    return {"parent_definition": "The monomial complement of the union of the divisor boxes of the six displayed socle corners.",
            "parent_standard_monomials": domain, "socle_circuit_supports": groups,
            "relation_coefficients": [[1, 1, 1], [1, 1, 1]], "parent_length": 17, "child_length": 15,
            "deletion_defect": 2, "characteristics_checked": [0, 2, 3, 7],
            "initial_ideal_count": len(states), "all_comparison_chambers": records,
            "central_weight": central, "two_exchange_cost_levels": [3, 4], "four_sector_witnesses": central_records,
            "conclusion": "Nonuniform rank-two relation matroids can have independent simultaneous exchanges at different costs. The single-cutoff vertex argument for uniform truncations does not extend to arbitrary quotient matroids.",
            "limit": "The example is supported on the full socle, so the already proved O(n^(3/2)) support-class bound applies. It is not an exponent-amplifying family. Nonuniform fixed defect outside that class is not excluded by the uniform-defect theorem."}


def nonuniform_grid_deletion_controls():
    """Two-point deletion on a torus graph: threshold proof and exact finite fans."""
    triangle = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]

    def integer_centroid(polygon):
        w = tuple(sum(a[i] for a in polygon) / len(polygon) for i in range(3))
        scale = lcm(*(value.denominator for value in w))
        w = tuple(int(value * scale) for value in w)
        common = gcd(*w)
        return tuple(value // common for value in w)

    def y_row_lengths(domain):
        rows = {}
        for a, b, c in domain:
            rows.setdefault((a, c), set()).add(b)
        require(all(values == set(range(len(values))) for values in rows.values()),
                "row-length witness encoding reconstructs the full staircase exactly")
        return [(a, c, len(values)) for (a, c), values in sorted(rows.items())]

    def shortest_representatives(u, v, w):
        best, queue = {(0, 0): (0, (0, 0, 0))}, [(0, (0, 0, 0), (0, 0))]
        while queue:
            cost, a, key = heappop(queue)
            if best[key] != (cost, a):
                continue
            for axis, step in enumerate(((1, -1), (1, 0), (0, 1))):
                other = ((key[0] + step[0]) % u, (key[1] + step[1]) % v)
                b = tuple(value + (axis == j) for j, value in enumerate(a))
                candidate = (cost + w[axis], b)
                if other not in best or candidate < best[other]:
                    best[other] = candidate
                    heappush(queue, (candidate[0], candidate[1], other))
        require(len(best) == u * v, "Cayley generators reach the entire actual group")
        return {key: a for key, (_, a) in best.items()}

    def regimes(u, v):
        # Each result has a fixed parent diagram and at most six candidates per color.
        box = set(product(range(u), (0,), range(v)))
        yield "y_expensive", [( -1, 1, -1)], box, [
            {(u - 1, 0, v - 1)}, {(u - 2, 0, v - 1), (u - 1, 0, v - 2)}]
        for k0 in range(v):
            for k1 in range(k0, v):
                inequalities = [(1, -1, 0)]
                for c in range(1, v):
                    for normal, cutoff in (((v - c, c - v, -c), k0), ((v - c, u - v + c, -c), k1)):
                        inequalities.append(normal if c <= cutoff else tuple(-x for x in normal))
                domain = set()
                for c in range(v):
                    for t in range(u):
                        a = (v - c) % v
                        domain.add((0, t, c) if c <= k0 or (c <= k1 and t < a)
                                   else (a, (t - a) % u, 0))
                candidates = [{(0, u - 1, 0)}, set()]
                for color in range(2):
                    values = [c for c in range(1, k0 + 1) if c % 2 == color]
                    if values:
                        candidates[color].add((0, u - 1, max(values)))
                    values = [c for c in range(k0 + 1, k1 + 1) if c % 2 == color]
                    if values:
                        c = min(values)
                        candidates[color].add((v - c, u - v + c - 1, 0))
                        candidates[color].update((0, v - c - 1, c) for c in (min(values), max(values)))
                    values = [c for c in range(k1 + 1, v) if c % 2 == color]
                    if values:
                        c = min(values)
                        candidates[color].add((v - c, u - 1, 0))
                yield "x_above_y_" + str(k0) + "_" + str(k1), inequalities, domain, candidates
        for k in range(v):
            inequalities = [(-1, 1, 0), (1, -1, 1)]
            for c in range(1, v):
                normal = (u - c, c, -c)
                inequalities.append(normal if c <= k else tuple(-x for x in normal))
            domain = set(product(range(u), range(v), (0,)))
            for c in range(1, k + 1):
                domain.difference_update((u - c, b, 0) for b in range(c, v))
                domain.update((0, b, c) for b in range(v - c))
            candidates = [set(), set()]
            for color in range(2):
                values = [a for a in range(u - k) if a % 2 == color]
                if values:
                    candidates[color].add((max(values), v - 1, 0))
                values = [c for c in range(1, k + 1) if c % 2 == color]
                if values:
                    c = max(values)
                    candidates[color].add((u - c, c - 1, 0))
                    candidates[color].update((0, v - c - 1, c) for c in (min(values), max(values)))
            yield "x_below_y_" + str(k), inequalities, domain, candidates

    fixtures = []
    for u, v, expected in ((4, 4, 14), (8, 4, 26), (12, 4, 26), (16, 4, 26), (6, 6, 38)):
        q = u * v

        def label(a):
            return ((a[0] + a[1]) % u, (a[2] - a[0]) % v)

        group = list(product(range(u), range(v)))
        relation_rows = [tuple(int(key[1] % 2 == color) for key in group) for color in range(2)]
        for color, row in enumerate(relation_rows):
            for a, b in ((1, -1), (1, 0), (0, 1)):
                shifted = tuple(int(((key[1] - b) % v) % 2 == color) for key in group)
                require(shifted == relation_rows[(color + b) % 2],
                        "the two relation rows form an invariant ideal, with invertible variable actions")
        reference = [key for key in group if key not in ((0, 0), (0, 1))]

        def quotient_column(key):
            if key in reference:
                return tuple(int(key == other) for other in reference)
            return tuple(-int(key[1] % 2 == other[1] % 2) for other in reference)

        columns = {key: quotient_column(key) for key in group}
        matrices = [[columns[((key[0] + a) % u, (key[1] + b) % v)] for key in reference]
                    for a, b in ((1, -1), (1, 0), (0, 1))]
        require(all(matvec(matrices[i], a) == matvec(matrices[j], b)
                    for i, j in combinations(range(3), 2) for a, b in zip(matrices[j], matrices[i])),
                "two-character quotient matrices commute over Z, including even characteristics")
        require(all(matvec(matrices[0], a) == b for a, b in zip(matrices[2], matrices[1])),
                "the actual quotient matrices satisfy xz=y")
        unit = columns[(0, 0)]
        for axis, period in ((1, u), (2, v)):
            vector = unit
            for _ in range(period):
                vector = matvec(matrices[axis], vector)
            require(vector == unit, "the stated grid periods hold on the cyclic vector")
        for key in group:
            vector = unit
            for _ in range(key[0]):
                vector = matvec(matrices[1], vector)
            for _ in range(key[1]):
                vector = matvec(matrices[2], vector)
            require(vector == columns[key], "the quotient representation is cyclic with the stated class columns")
        socle_parent = set(shortest_representatives(u, v, (u + 1, u, u)).values())
        socle = {a for a in socle_parent if all(tuple(value + (axis == j) for j, value in enumerate(a)) not in socle_parent
                                              for axis in range(3))}
        require(len(socle) == 2 * v - 1, "the same family has unbounded parent initial socle type")
        records, child_states, parent_states, parent_for_basis = [], {}, set(), {}
        for name, inequalities, domain, candidates in regimes(u, v):
            polygon = triangle
            for normal in inequalities:
                polygon = clip(polygon, normal, 1)
                if not polygon:
                    break
            if not polygon:
                continue
            require(len(domain) == q and len({label(a) for a in domain}) == q,
                    "each threshold regime is a full parent transversal")
            require(all(group_candidates <= domain and len(group_candidates) <= 6
                        for group_candidates in candidates), "uniformly bounded color-max candidate lists")
            parent_states.add(tuple(sorted(domain)))
            for even in sorted(candidates[0]):
                for odd in sorted(candidates[1]):
                    piece = polygon
                    for chosen, group_candidates in zip((even, odd), candidates):
                        for other in group_candidates - {chosen}:
                            piece = clip(piece, tuple(x - y for x, y in zip(chosen, other)), 1)
                            if not piece:
                                break
                        if not piece:
                            break
                    if not piece:
                        continue
                    w = integer_centroid(piece)
                    require(all(dot(w, normal) > 0 for normal in inequalities), "strict parent-regime witness")
                    actual = shortest_representatives(u, v, w)
                    require(set(actual.values()) == domain, "independent shortest paths check the canonical parent formula")
                    for color, chosen in enumerate((even, odd)):
                        require(label(chosen)[1] % 2 == color and
                                all(dot(w, chosen) > dot(w, a) for a in domain
                                    if a != chosen and label(a)[1] % 2 == color),
                                "endpoint candidates give the unique maxima among all actual group classes")
                    basis = tuple(sorted(domain - {even, odd}))
                    child_states.setdefault(basis, w)
                    parent_for_basis.setdefault(basis, domain)
                    records.append({"regime": name, "positive_integer_weights": w,
                                    "parent_y_row_lengths": y_row_lengths(domain),
                                    "omitted_even_and_odd": [even, odd], "standard_y_row_lengths": y_row_lengths(basis)})
        require(len(child_states) == expected, "complete finite counts from the proved threshold partition")
        gaussian_checks = 0
        if q <= 32:
            for basis, w in child_states.items():
                ground = parent_for_basis[basis]
                evaluated = {a: columns[label(a)] for a in ground}
                for characteristic in (0, 2, 3, 7):
                    require(greedy(ground, evaluated, w, characteristic) == set(basis),
                            "independent Gaussian selection on every small-fixture initial basis")
                    gaussian_checks += 1
        fixtures.append({"u": u, "v": v, "parent_length": q, "child_length": q - 2,
                         "parent_fan_count": len(parent_states), "child_fan_count": len(child_states),
                         "parent_initial_socle_type_at_weights_u_plus_one_u_u": len(socle),
                         "all_regime_and_maximum_cells": records, "gaussian_checks": gaussian_checks})

    local = []
    for v in (4, 6, 8):
        u = v * v
        x, y, z = 2 * v - 2, 4, 3 * v * v - 7 * v + 10
        a, b, c, e = (v - 3, u - 1, 0), (v - 2, u - v + 1, 0), (v - 1, u - v, 0), (0, v - 3, 2)
        sector_records = []
        for s, t in product((-1, 1), repeat=2):
            w = (x + s, y, z + (v - 2) * s // 2 + t)
            parent = shortest_representatives(u, v, w)
            omitted = [max((point for key, point in parent.items() if key[1] % 2 == color),
                           key=lambda point: dot(w, point)) for color in range(2)]
            require(omitted == [e if t == 1 else b, c if s == 1 else a],
                    "uniform four-sector witnesses in the nonsocle family")
            sector_records.append({"positive_integer_weights": w, "omitted_even_and_odd": omitted})
        require(len({tuple(record["omitted_even_and_odd"]) for record in sector_records}) == 4,
                "independent two-height exchanges occur for every tested parameter")
        local.append({"u": u, "v": v, "central_weight": [x, y, z],
                      "even_maximum_cost": 6 * v * v - 10 * v + 8,
                      "odd_maximum_cost": 6 * v * v - 8 * v + 2,
                      "four_sector_witnesses": sector_records})
    return {"uniform_scope": "Every even v>=4 and every positive multiple u of v with u>=v, over every field: J=(xz-y,y^u-1,z^v-1), I=J+(sum_{a=0}^{u-1} sum_{0<=b<v, b even} y^a*z^b).",
            "length_and_upper_bound": "The child has n=uv-2 and g(I)<=18*v^2+34*v+2=O(n). In characteristic zero it removes the two points (1,1,1),(-1,1,-1) from the full torus grid graph.",
            "proof_budget": {"x_above_y_parent_regimes": "at most v*(v+1)/2", "max_candidates_per_color": 6,
                             "x_below_y_parent_regimes": "at most v", "max_candidates_per_color_in_that_region": 4,
                             "y_expensive_child_states": 2},
            "nonannihilation": "On the two-dimensional relation ideal, y acts as identity and x,z swap its displayed basis. These actions are invertible, so this is outside the variable-annihilated full-socle construction.",
            "unbounded_initial_socle_type": "At weights (u+1,u,u), the parent initial algebra has 2*v-1 socle monomials. The count bound does not assume bounded initial type.",
            "fixtures": fixtures, "uniform_local_crossing_controls": local,
            "decision": "Independent exchanges persist outside socle-supported constructions, but the whole divisible-grid family has only linear growth and cannot provide the required n^(2+epsilon) lower family. No unrestricted upper advance or novelty certification."}


def coordinate_moment_range_controls(residual_duality, point_fixture):
    """Exact box moment extrema, independent height optimization and colon checks."""
    fixtures, extremizer_checks, residual_checks = [], 0, 0
    specifications = [sides for d in (2, 3)
                      for sides in product(range(1, 5), repeat=d)]
    specifications += [(5, 5), (5, 5, 5), (2, 3, 7),
                       (2, 2, 2, 2), (2, 2, 3, 3)]

    def bounds(n, side, columns):
        full, remainder = divmod(n, side)
        upper = full * comb(side, 2) + comb(remainder, 2)
        layers, remainder = divmod(n, columns)
        lower = columns * comb(layers, 2) + layers * remainder
        return lower, upper

    for sides in specifications:
        size = 1
        for side in sides:
            size *= side
        box = set(product(*(range(side) for side in sides)))
        records = []
        for axis, side in enumerate(sides):
            columns = list(product(*(range(sides[j]) for j in range(len(sides)) if j != axis)))
            count = len(columns)
            # Optimize over every column-height assignment, without imposing a downset.
            # The separate constructions below show both relaxed optima are attained by downsets.
            minimum, maximum = [0], [0]
            for _ in columns:
                next_min, next_max = [None] * (len(minimum) + side), [None] * (len(maximum) + side)
                for previous, (low, high) in enumerate(zip(minimum, maximum)):
                    for height in range(side + 1):
                        target, cost = previous + height, comb(height, 2)
                        candidate = low + cost
                        if next_min[target] is None or candidate < next_min[target]:
                            next_min[target] = candidate
                        candidate = high + cost
                        if next_max[target] is None or candidate > next_max[target]:
                            next_max[target] = candidate
                minimum, maximum = next_min, next_max

            def staircase(heights):
                result = set()
                for column, height in zip(columns, heights):
                    for value in range(height):
                        point = column[:axis] + (value,) + column[axis:]
                        result.add(point)
                return result

            for n in range(size + 1):
                low, high = bounds(n, side, count)
                require((low, high) == (minimum[n], maximum[n]),
                        (sides, axis, n, "independent column-height optimization"))
                full, rest = divmod(n, side)
                upper_heights = [side if j < full else rest if j == full else 0 for j in range(count)]
                layers, rest = divmod(n, count)
                lower_heights = [layers + int(j < rest) for j in range(count)]
                for heights, expected in ((lower_heights, low), (upper_heights, high)):
                    state = staircase(heights)
                    require(len(state) == n and state <= box and
                            all(tuple(x - int(i == j) for j, x in enumerate(a)) in state
                                for a in state for i in range(len(sides)) if a[i]),
                            (sides, axis, n, "attaining witness is a box downset"))
                    require(sum(a[axis] for a in state) == expected,
                            (sides, axis, n, "extremal moment witness"))
                    complement = {tuple(s - 1 - x for s, x in zip(sides, a)) for a in box - state}
                    shift = (size - n) * (side - 1) - sum(a[axis] for a in box)
                    require(sum(a[axis] for a in complement) == expected + shift,
                            "reflected-complement moment translation")
                    extremizer_checks += 1
                other_low, other_high = bounds(size - n, side, count)
                require(high - low == other_high - other_low,
                        "coordinate width is invariant under complementary length")
                require(0 <= high - low <= high and high + 1 <= max(1, n * side),
                        "width cap implies the simpler product-degree bound")
            records.append({"axis": axis, "column_count": count,
                            "minimum_moments_by_length": minimum,
                            "maximum_moments_by_length": maximum,
                            "moment_widths_by_length": [hi - lo for lo, hi in zip(minimum, maximum)]})
        fixtures.append({"box_sides": sides, "box_volume": size, "coordinate_extrema": records})

    for fixture in residual_duality["fixtures"]:
        size, n = fixture["parent_length"], fixture["child_length"]
        for witness in fixture["all_comparison_chamber_checks"]:
            child, residual = witness["child_standard_monomials"], witness["residual_standard_monomials"]
            for axis in range(3):
                low, high = bounds(n, 2, size // 2)
                child_moment, residual_moment = sum(a[axis] for a in child), sum(a[axis] for a in residual)
                require(low <= child_moment <= high and
                        residual_moment == child_moment + size // 2 - n,
                        "actual grid quotients obey the moment interval and translation in every checked field")
                residual_checks += 1

    finite_fields = []
    for prime in (2, 3, 5):
        for dimension in (2, 3):
            size = prime ** dimension
            comparisons = []
            for n in range(1, size):
                r = min(n, size - n)
                low, high = bounds(r, prime, size // prime)
                base = prime * prime * (r // prime) + (r % prime) ** 2
                require(base <= 4 * (high + 1) and high + 1 <= base,
                        "Theorem 8 base and positive moment-cap expression agree up to absolute constants")
                comparisons.append({"length": n, "smaller_length": r, "published_base": base,
                                    "moment_upper_plus_one": high + 1, "moment_width_plus_one": high - low + 1})
            finite_fields.append({"prime": prime, "dimension": dimension, "comparisons": comparisons})

    points, n = point_fixture["points"], point_fixture["n"]
    require(all(len({point[axis] for point in points}) == n for axis in range(3)),
            "the known nongeneric point fixture has full coordinate elimination degrees")
    require(point_fixture["point_initial_ideals"] > point_fixture["generic_corner_cuts"],
            "full coordinate degrees do not force the generic corner-cut count")
    field_equation_counters = []
    for prime in (2, 3, 5):
        length = prime + 1
        modulus = (0,) * length + (1,)
        variable = (0, 1) + (0,) * (length - 2)
        power = univariate_power_mod(variable, prime, modulus)
        remainder = [(a - b) % prime for a, b in zip(power, variable)]
        require(any(remainder), "field equations need not vanish in a nonreduced finite algebra")
        field_equation_counters.append({"prime": prime, "length": length,
                                        "quotient": "F_p[x,y,z]/(x^(p+1),y,z)",
                                        "remainder_of_x_to_p_minus_x": remainder})
    return {"scope": "Every finite-colength ideal over every field, using any univariate grid cover. The n-based consequence uses its minimal coordinate elimination degrees.",
            "uniform_bound": "g(I)<=C_d*(product_i (W_i+1))^((d-1)/(d+1)), with W_i the exact coordinate-moment width of n-element downsets in its grid box.",
            "worst_case_limit": "Minimal degrees s_i<=n recover the old unrestricted exponent. No qualifying worst-case asymptotic improvement or novelty is asserted.",
            "box_fixtures": fixtures, "attaining_downset_checks": extremizer_checks,
            "existing_quotient_coordinate_checks": residual_checks,
            "finite_field_theorem_comparisons": finite_fields,
            "finite_algebra_field_equation_counters": field_equation_counters,
            "full_coordinate_degree_countercontrol": {
                "points": points, "coordinate_degrees": [n] * 3,
                "point_initial_ideals": point_fixture["point_initial_ideals"],
                "generic_corner_cuts": point_fixture["generic_corner_cuts"]}}


def grading_and_perturbation_controls(monogenic, point_fixture):
    """Grading moment constraints and loss of the fan under small coefficient changes."""
    graded, perturbed = [], []

    def rank(rows):
        columns = {(i,): row for i, row in enumerate(rows)}
        return len(greedy(columns, columns, (1,), 0))

    for n in range(3, 8):
        modulus = (0,) * n + (1,)
        ground, states = hyperbolic_cross(n, 3), finite_staircases(n, 3)
        for degrees in ((1, 1, 1), (1, 2, 3), (1, 2, 4)):
            images = [tuple(int(j == degree) for j in range(n)) for degree in degrees]
            columns = cyclic_univariate_model(modulus, images, ground)
            bases, moments = [], []
            for state in states:
                minor = determinant([columns[a] for a in sorted(state)])
                require(bool(minor) == ({dot(a, degrees) for a in state} == set(range(n))),
                        "independent normal-form determinant and graded-degree criterion agree")
                if minor:
                    moment = tuple(sum(a[i] for a in state) for i in range(3))
                    require(dot(degrees, moment) == n * (n - 1) // 2,
                            "all homogeneous monomial bases have the same grading moment")
                    bases.append(sorted(state))
                    moments.append(moment)
            differences = [tuple(x - y for x, y in zip(moment, moments[0])) for moment in moments[1:]]
            dimension = rank(differences)
            require(dimension <= 2, "grading confines the moment polytope to a plane")
            if degrees == (1, 1, 1):
                require(len(bases) == 3 and dimension == 2,
                        "the diagonal curvilinear embedding has exactly its three axis bases")
            graded.append({"n": n, "coordinate_t_degrees": degrees,
                           "all_staircases_checked": len(states), "independent_staircases": bases,
                           "basis_moments": moments, "affine_moment_dimension": dimension,
                           "grading_moment": n * (n - 1) // 2})

        powers = (1, n, n * n)
        largest_change = max(max(power - 1, *(comb(power, j) for j in range(2, n)))
                             for power in powers[1:])
        denominator = 1000000 * (largest_change + 1)
        # These are the numerators of x=t, y=t+epsilon*((1+t)^n-1-t), etc.
        # Scaling coordinates by nonzero constants preserves every initial staircase.
        numerators = [tuple(int(j == 1) for j in range(n))]
        for power in powers[1:]:
            numerators.append(tuple(0 if j == 0 else denominator + power - 1 if j == 1
                                    else comb(power, j) for j in range(n)))
        columns = cyclic_univariate_model(modulus, numerators, ground)
        minors = [determinant([columns[a] for a in sorted(state)]) for state in states]
        require(all(minors), "every staircase remains independent at the small rational perturbation")
        endpoint = next(item for item in monogenic["centered_curvilinear_fixtures"]
                        if item["n"] == n and item["dimension"] == 3)
        require(len(states) == endpoint["all_staircase_determinants_checked"],
                "the endpoint nonvanishing theorem covers the same complete staircase set")
        differences = []
        for axis, image in enumerate(numerators[1:], 1):
            for degree, coefficient in enumerate(image):
                if coefficient:
                    differences.append(tuple(int(i == axis) - degree * int(i == 0) for i in range(3)))
        require(rank(differences) == 3,
                "the two reduced graph relations have trivial homogeneity space")
        require(all(image[1] for image in numerators),
                "all coordinates retain nilpotency index n and minimal elimination degree n")
        require(Fraction(largest_change, denominator) < Fraction(1, 1000000),
                "every original graph coefficient changes by less than one millionth")
        coefficient_bound = 2 + n * (n * n) ** (n - 1)
        determinant_height_bound = factorial(n) * coefficient_bound ** (n * (n - 1) // 2)
        perturbed.append({"n": n, "epsilon": "1/" + str(denominator),
                          "coefficient_distance_sup_norm": str(Fraction(largest_change, denominator)),
                          "coordinate_image_numerators_low_to_high": numerators,
                          "coordinate_image_denominators": [1, denominator, denominator],
                          "all_staircases_in_lexicographic_order": [sorted(state) for state in states],
                          "nonzero_scaled_basis_determinants": [str(value) for value in minors],
                          "homogeneity_dimension": 0, "coordinate_elimination_degrees": [n] * 3,
                          "uniform_integer_coefficient_bound": coefficient_bound,
                          "uniform_determinant_polynomial_height_bound": str(determinant_height_bound),
                          "uniform_safe_epsilon_denominator": str(2 * determinant_height_bound + 2)})
    require(point_fixture["n"] == 7 and point_fixture["generic_corner_cuts"] == 36,
            "the independent existing seven-box generic corner-cut count is available")
    return {"graded_scope": "Any multigraded finite quotient over every field has all monomial-basis moments on the grading affine subspace. The dimension-sensitive count uses the established relative-lattice vertex estimate.",
            "uniform_bound": "If h independent integer gradings make I homogeneous and r=d-h, then g(I)<=1 for r=0, <=2 for r=1, and O_d(n^(2*r*(r-1)/(r+1))) for r>=2.",
            "graded_fixtures": graded, "perturbation_fixtures": perturbed,
            "perturbation_scope": "Over Q, for each n>=3, a finite excluded set of epsilon values contains every possible staircase determinant zero. Epsilon=0 has three initials, whereas sufficiently small nonzero epsilon has C_3(n). The same fixed abstract algebra and all three elimination degrees n persist.",
            "seven_box_fan_counts": {"epsilon_zero": 3, "coefficient_distance_less_than_one_millionth": 36},
            "priority": "Consequences of established multigrading, lattice-polytope and open-minor arguments, not a new unrestricted asymptotic theorem or certified novelty."}


def attention_rank_sanity_control():
    """An exact Vandermonde witness, not an attention-compression algorithm."""
    parameters = (1, 2, 3, 4)
    kernel = [[t ** i for t in parameters] for i in range(len(parameters))]
    row_sums = list(map(sum, kernel))
    attention = [[Fraction(value, scale) for value in row] for row, scale in zip(kernel, row_sums)]
    product_formula = 1
    for a, b in combinations(parameters, 2):
        product_formula *= b - a
    require(determinant(kernel) == product_formula == 12, "ordinary Vandermonde determinant")
    require(determinant(attention) == Fraction(1, 10000), "positive row normalization preserves full rank")
    return {"query_scalars": list(range(4)), "key_scalars_symbolically": ["log(" + str(t) + ")" for t in parameters],
            "logit_formula": "S_ij=i*log(t_j), a rank-one outer product; d=1, with i=0,1,2,3.",
            "exponentiated_kernel_matrix": kernel, "row_sums": row_sums,
            "row_softmax_attention_matrix": [[str(value) for value in row] for row in attention],
            "kernel_determinant": str(determinant(kernel)), "attention_determinant": str(determinant(attention)),
            "scope": "For all sequence lengths, distinct positive t_j give a full-rank Vandermonde attention matrix despite rank-at-most-one logits (rank one when N>=2). This classical observation does not rule out accurate approximate low-rank computation or establish novelty."}


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
    binomial = binomial_controls()
    truncated_partition = truncated_partition_point_control()
    sectors = pure_power_sector_controls(binomial, truncated_partition)
    coarse_path = coarse_class_path_control(truncated_partition)
    normal_curve = normal_curve_lower_controls()
    relative_socle = [relative_socle_fixture(length) for length in (2, 3, 4, 5)]
    relative_lower = relative_socle_lower_controls()
    hook = [hook_deformation_fixture(r, s) for r, s in ((2, 2), (2, 3), (3, 3))]
    sparse = [sparse_corner_deformation_fixture(k, r, s)
              for k, r, s in ((2, 5, 5), (3, 7, 7), (4, 9, 9), (3, 11, 13))]
    thick = [thick_hook_fixture(a, b, r, s)
             for a, b, r, s in ((1, 2, 3, 5), (2, 2, 5, 5), (2, 3, 5, 7), (2, 2, 9, 11))]
    power_control = thick_power_substitution_control()
    untruncated = untruncated_convolution_controls()
    factorial_boundary = [factorial_boundary_fixture(length) for length in range(1, 7)]
    arm_counts = {f["L"]: {v["characteristic"]: v["exact_initial_ideal_count"] for v in f["fields"]}
                  for f in factorial_boundary}
    independent_arms = [independent_boundary_arms_fixture(t, arm_counts[t]) for t in (2, 3, 4, 5)]
    collisions = monomial_curve_collision_controls(truncated_partition)
    monogenic = monogenic_embedding_controls()
    sparse_prefix = sparse_prefix_interpolation_control()
    jordan_prefix = jordan_filtered_convolution_controls()
    unequal_jordan = unequal_weight_jordan_countercontrols()
    codimension_one = codimension_one_binomial_quotient_controls()
    uniform_deletion = uniform_deletion_and_dual_relation_controls()
    residual_duality = box_residual_duality_controls(codimension_one)
    two_level_exchange = nonuniform_two_level_exchange_control()
    grid_deletion = nonuniform_grid_deletion_controls()
    coordinate_moments = coordinate_moment_range_controls(residual_duality, truncated_partition)
    grading_controls = grading_and_perturbation_controls(monogenic, truncated_partition)
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
        "binomial_parallel_class_control": {
            "uniform_conclusion": "For every colength-n binomial ideal in fixed dimension d>=2, g(I)=O_d(n^(d-1)(log(n+1))^((d-1)^2)), by the known normal-form reduction, hyperbolic cross and Minkowski-sum vertex bound.",
            "priority": "Recorded as a consequence of existing tools, not an established novel theorem or an unrestricted advance.",
            "fixtures": binomial,
            "nonpartition_point_control": nonpartition_point_control(),
        },
        "truncated_partition_point_control": truncated_partition,
        "pure_power_sector_control": {
            "uniform_conclusion": "Primitive pure-power relations supply d sector boxes. The lattice-polytope moment bound gives O_A(n^(d-1)) for every finite quotient on a fixed monomial curve, over every field; this rules out a fixed-curve asymptotic separation from the generic order.",
            "scope": "A restricted consequence of the established state-polytope vertex estimate; not an unrestricted improvement or an established novelty claim.",
            "witness_checks": sectors,
        },
        "coarse_class_path_control": coarse_path,
        "normal_curve_lower_control": {
            "uniform_conclusion": "The ideals of n distinct positive samples (t,t^2,...,t^d) have Theta_d(n^(d-1)) initial ideals. The full proof combines the fixed-curve upper bound with an explicit open family of ray allocations.",
            "scope": "A sharp restricted-family control below the existing generic logarithmic order. It does not improve the unrestricted lower bound and novelty has not been established.",
            "point_description": "For each fixture take t=1,...,n over Q; ray residue i has monomials x_i*x_d^k, with x_0=1.",
            "fixtures": normal_curve,
        },
        "relative_socle_control": {
            "uniform_conclusion": "For the stated finite monomial reference (or the justified product reference), relations supported on the relative socle in s specified variables and closed under the remaining multiplications give O_d(n^((2*d-s)*(d-1)/(d+1))) initial ideals. In dimension three the layer argument improves s=1 to O(n^(7/3)); s=2 gives O(n^2), with a matching family of lengths 6L.",
            "scope": "A restricted mechanism with genuinely surviving z multiples; it cannot beat the existing unrestricted generic lower order. No novelty or unrestricted improvement is asserted.",
            "fixtures": relative_socle, "symbolic_lower_certificates": relative_lower,
        },
        "hook_first_order_deformation_control": {
            "uniform_conclusion": "Every flat first-order deformation of (y^r,y*z,z^s), r,s>=2, has O((r+s)^2) initial ideals as a length-2(r+s-1) algebra in three variables. Three connecting columns suffice after deleting a coloop and two blocks with O(r+s) internal comparison hyperplanes.",
            "scope": "A proved obstruction for this sparse first-order family, not an unrestricted advance or established novelty claim.",
            "fixtures": hook,
        },
        "sparse_growing_corner_deformation_control": {
            "uniform_conclusion": "Every flat first-order deformation of the sparse plane ideal (y^r,y^k*z,y^(k-1)*z^2,...,y*z^k,z^s), k>=2 and r,s>=2k+1, has O(n^2) initial ideals at quotient length n=2*(r+s-1+k*(k-1)/2). Three residue blocks have total interface dimension at most twelve and a linear number of internal comparison hyperplanes.",
            "scope": "A uniform growing-corner obstruction, not an unrestricted advance or established novelty claim. These fixtures check exact selected weights and do not enumerate the full fan.",
            "fixtures": sparse,
            "pairwise_intersection_negative_control": interface_pairwise_countercontrol(),
        },
        "thick_hook_deformation_control": {
            "uniform_conclusion": "For all compatible thick-hook first-order tails, q<=12*a*b and at most 120*m internal hyperplanes give O(n^2*(a*b)^4) initial ideals by optional-column refinement. Together with the prior annihilator bound this is O(min(n^(7/3),n^2*(a*b)^4)). Fixed thickness is therefore insufficient for a stronger unrestricted lower bound.",
            "scope": "A parameter-sensitive restricted upper bound, not an unrestricted advance or novelty claim. These are selected-weight controls, not fan enumeration.",
            "fixtures": thick, "power_substitution_negative_control": power_control,
            "koszul_support_matching_negative_control": thick_koszul_countercontrol(),
            "untruncated_convolution_rank_and_basis_control": untruncated,
        },
        "factorial_boundary_control": {
            "uniform_conclusion": "For I_L=(x^2,z^L-x*A_L(z)) over characteristic zero or p>L-1, length is 2L and there are exactly L+1 global initial ideals, versus two positive coherent bases in the corresponding untruncated convolution. Two disjoint arms give length 4t+2 and exactly (t+1)^2 initial ideals; their interaction dimension is zero.",
            "scope": "A proved boundary mechanism and sharp restricted quadratic control below the existing generic lower order. Classical Toeplitz/Padé ingredients; novelty and an unrestricted advance are not established.",
            "one_arm_fixtures": factorial_boundary, "independent_arm_fixtures": independent_arms,
        },
        "positive_curve_collision_control": {
            "uniform_conclusion": "Positive monomial-curve divisors of the same length and exponent tuple have identical monomial residue matroids and identical positive Groebner fans, including all collisions with multiplicity. Collision within this class cannot amplify the fan.",
            "priority": "A consequence of classical sparse positive-root bounds and confluent generalized Vandermonde theory; no new interpolation theorem is claimed.",
            "checks": collisions,
        },
        "monogenic_embedding_control": {
            "uniform_conclusion": "Over characteristic zero, a finite algebra admits an embedding with the labeled generic corner-cut fan if and only if it is monogenic. For every monogenic algebra a nonempty Zariski open set of coordinate tuples has this fan. The fixed abstract curvilinear algebra K[u]/(u^n) has both one-initial-ideal and C_d(n)-initial-ideal embeddings.",
            "scope": "Generic embeddings are not asserted to maximize the fan count; the existing seven-point countercontrol rules that out. This does not improve G_d or F_d and does not complete the goal.",
            "checks": monogenic,
        },
        "sparse_prefix_interpolation_control": sparse_prefix,
        "jordan_filtered_convolution_control": jordan_prefix,
        "unequal_weight_jordan_countercontrol": unequal_jordan,
        "codimension_one_binomial_quotient_control": codimension_one,
        "uniform_deletion_and_dual_relation_control": uniform_deletion,
        "attention_exact_rank_sanity_control": attention_rank_sanity_control(),
        "box_residual_duality_control": residual_duality,
        "nonuniform_two_level_exchange_control": two_level_exchange,
        "nonuniform_grid_deletion_control": grid_deletion,
        "coordinate_moment_range_control": coordinate_moments,
        "grading_and_perturbation_control": grading_controls,
        "summary": {"volume_fixtures": 30, "socle_fixtures": len(fixtures),
                    "exact_chambers": sum(f["exact_arrangement_chambers"] for f in fixtures),
                    "independent_greedy_checks": sum(f["greedy_checks"] for f in fixtures),
                    "overlapping_socle_fixtures": len(overlapping),
                    "overlapping_row_and_greedy_checks": sum(f["row_relation_and_socle_greedy_checks"] for f in overlapping),
                    "overlapping_full_matrix_checks": sum(f["full_cyclic_matrix_greedy_checks"] for f in overlapping),
                    "rectangular_boxes_checked": sum(f["finite_boxes_checked"] for f in rectangular),
                    "rectangular_prefix_checks": sum(len(f["prefix_checks"]) for f in rectangular),
                    "cartesian_rank_delay_fixtures": len(cartesian),
                    "binomial_fixtures": len(binomial),
                    "binomial_gaussian_checks": sum(f["gaussian_checks"] for f in binomial),
                    "binomial_exact_chambers": sum(f["exact_arrangement_chambers"] for f in binomial),
                    "nonpartition_point_circuits": 1,
                    "truncated_partition_staircases_checked": truncated_partition["staircases_checked"],
                    "truncated_partition_rational_checks": truncated_partition["independent_rational_gaussian_checks"],
                    "pure_power_sector_witnesses": len(sectors),
                    "coarse_class_path_intervals": coarse_path["same_connected_class_region"]["open_intervals_checked"],
                    "normal_curve_fixtures": len(normal_curve),
                    "normal_curve_allocations_checked": sum(f["allocation_certificates_checked"] for f in normal_curve),
                    "normal_curve_independent_rank_checks": sum(f["independent_modular_rank_check"] is not None for f in normal_curve),
                    "relative_socle_fixtures": len(relative_socle),
                    "relative_socle_row_and_greedy_checks": sum(f["row_and_greedy_checks"] for f in relative_socle),
                    "relative_socle_full_matrix_checks": sum(f["full_matrix_checks"] for f in relative_socle),
                    "relative_socle_lower_certificates": sum(f["certificates_checked"] for f in relative_lower),
                    "hook_deformation_fixtures": len(hook),
                    "hook_deformation_row_and_greedy_checks": sum(f["row_and_greedy_checks"] for f in hook),
                    "hook_deformation_full_matrix_checks": sum(f["full_matrix_checks"] for f in hook),
                    "sparse_corner_fixtures": len(sparse),
                    "sparse_corner_row_and_greedy_checks": sum(f["row_and_greedy_checks"] for f in sparse),
                    "sparse_corner_full_matrix_checks": sum(f["full_matrix_checks"] for f in sparse),
                    "sparse_corner_interface_flag_checks": sum(f["interface_flag_checks"] for f in sparse),
                    "sparse_corner_optional_recovery_checks": sum(f["optional_recovery_checks"] for f in sparse),
                    "thick_hook_fixtures": len(thick),
                    "thick_hook_row_and_greedy_checks": sum(f["row_and_greedy_checks"] for f in thick),
                    "thick_hook_full_matrix_checks": sum(f["full_matrix_checks"] for f in thick),
                    "power_substitution_basis_checks": power_control["basis_lifting_checks"],
                    "untruncated_convolution_basis_checks": untruncated["exact_basis_checks"],
                    "factorial_boundary_fixtures": len(factorial_boundary),
                    "factorial_boundary_prefix_minors": sum(len(f["prefix_minors"]) for f in factorial_boundary),
                    "factorial_boundary_row_and_greedy_checks": sum(f["row_and_greedy_checks"] for f in factorial_boundary),
                    "factorial_boundary_full_matrix_checks": sum(f["full_matrix_checks"] for f in factorial_boundary),
                    "independent_boundary_arm_fixtures": len(independent_arms),
                    "independent_boundary_arm_row_and_greedy_checks": sum(f["row_and_greedy_checks"] for f in independent_arms),
                    "independent_boundary_arm_full_matrix_checks": sum(f["full_matrix_checks"] for f in independent_arms),
                    "positive_collision_profiles": len(collisions["profiles"]),
                    "positive_collision_staircase_rank_checks": sum(f["staircase_rank_checks"] for f in collisions["profiles"]),
                    "positive_collision_fan_witness_checks": sum(2 * f["fan_witnesses_checked_in_both_representations"] for f in collisions["profiles"]),
                    "monogenic_fixed_algebra_fixtures": len(monogenic["fixed_algebra_fixtures"]),
                    "monogenic_all_minors_checked": sum(f["all_rank_n_minors_checked"] for f in monogenic["fixed_algebra_fixtures"]),
                    "centered_curvilinear_staircase_checks": sum(f["all_staircase_determinants_checked"] for f in monogenic["centered_curvilinear_fixtures"]),
                    "sparse_prefix_exact_matching_rank_controls": 1,
                    "sparse_prefix_uniform_dual_certificates": len(sparse_prefix["fixtures"]),
                    "jordan_graded_map_rank_checks": sum(f["graded_multiplication_map_ranks_checked"] for f in jordan_prefix["fixtures"]),
                    "jordan_three_node_string_rank_checks": jordan_prefix["three_node_string_rank_checks"],
                    "jordan_filtered_prefix_rank_checks": sum(len(f["prefix_rank_checks"]) for f in jordan_prefix["fixtures"]),
                    "jordan_full_tag_power_diagonal_cost_checks": sum(len(f.get("full_tag_power_diagonal_cost_controls", [])) for f in jordan_prefix["fixtures"]),
                    "unequal_jordan_staircases_checked": sum(len(f["all_staircase_certificates"]) for f in unequal_jordan["embeddings"]),
                    "unequal_jordan_coherent_basis_checks": sum(f["exact_initial_ideal_count"] for f in unequal_jordan["embeddings"]),
                    "unequal_jordan_infeasibility_certificates": sum(sum(c["status"] == "independent_but_incoherent" for c in f["all_staircase_certificates"]) for f in unequal_jordan["embeddings"]),
                    "codimension_one_binomial_quotient_fixtures": len(codimension_one["fixtures"]),
                    "codimension_one_staircase_rank_checks": sum(f["monomial_rank_and_staircase_checks"] for f in codimension_one["fixtures"]),
                    "codimension_one_coherent_basis_checks": sum(g["initial_ideal_count"] for f in codimension_one["fixtures"] for g in f["fans"]),
                    "codimension_one_infeasibility_certificates": sum(sum(c["status"] == "independent_but_incoherent" for c in g["all_staircase_certificates"]) for f in codimension_one["fixtures"] for g in f["fans"]),
                    "uniform_deletion_fourier_minors_checked": sum(f["all_row_and_column_minors_checked"] for f in uniform_deletion["uniform_fixtures"]),
                    "uniform_deletion_staircase_rank_checks": sum(len(f["all_staircase_certificates"]) for f in uniform_deletion["uniform_fixtures"]),
                    "uniform_deletion_coherent_basis_checks": sum(f["initial_ideal_count"] for f in uniform_deletion["uniform_fixtures"]),
                    "uniform_deletion_gaussian_checks": 2 * sum(f["initial_ideal_count"] for f in uniform_deletion["uniform_fixtures"]),
                    "nonuniform_complementary_basis_checks": len(uniform_deletion["nonuniform_two_relation_countercontrol"]["complementary_basis_checks"]),
                    "attention_full_rank_softmax_witnesses": 1,
                    "box_residual_duality_fixtures": len(residual_duality["fixtures"]),
                    "box_residual_duality_chamber_checks": sum(len(f["all_comparison_chamber_checks"]) for f in residual_duality["fixtures"]),
                    "box_residual_duality_gaussian_checks": 2 * sum(len(f["all_comparison_chamber_checks"]) for f in residual_duality["fixtures"]),
                    "nongrid_box_cone_checks": len(residual_duality["nongrid_box_cones"]["checks"]),
                    "nonbox_residual_transfer_counter_witnesses": 2,
                    "nonuniform_two_level_exchange_chambers": len(two_level_exchange["all_comparison_chambers"]),
                    "nonuniform_two_level_exchange_gaussian_checks": 4 * len(two_level_exchange["all_comparison_chambers"]),
                    "nonuniform_two_level_exchange_local_sectors": 4,
                    "nonuniform_grid_deletion_fixtures": len(grid_deletion["fixtures"]),
                    "nonuniform_grid_deletion_regime_cells": sum(len(f["all_regime_and_maximum_cells"]) for f in grid_deletion["fixtures"]),
                    "nonuniform_grid_deletion_gaussian_checks": sum(f["gaussian_checks"] for f in grid_deletion["fixtures"]),
                    "nonsocle_uniform_crossing_sectors": 4 * len(grid_deletion["uniform_local_crossing_controls"]),
                    "coordinate_moment_box_fixtures": len(coordinate_moments["box_fixtures"]),
                    "coordinate_moment_attaining_downset_checks": coordinate_moments["attaining_downset_checks"],
                    "coordinate_moment_existing_quotient_checks": coordinate_moments["existing_quotient_coordinate_checks"],
                    "finite_field_bound_comparisons": sum(len(f["comparisons"]) for f in coordinate_moments["finite_field_theorem_comparisons"]),
                    "graded_staircase_rank_checks": sum(f["all_staircases_checked"] for f in grading_controls["graded_fixtures"]),
                    "small_perturbation_staircase_checks": sum(len(f["nonzero_scaled_basis_determinants"]) for f in grading_controls["perturbation_fixtures"])},
        "goal_qualifying_advance": False,
    }
    target = ROOT / "results/research/structural_complexity_checks.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
