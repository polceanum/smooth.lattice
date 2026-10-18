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
from math import comb, factorial, gcd, lcm
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
                    "untruncated_convolution_basis_checks": untruncated["exact_basis_checks"]},
        "goal_qualifying_advance": False,
    }
    target = ROOT / "results/research/structural_complexity_checks.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
