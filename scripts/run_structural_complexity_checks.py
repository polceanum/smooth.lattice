#!/usr/bin/env python3
"""Exact, bounded controls for the general initial-ideal research goal.

This checks proof-route obstructions, conditional lemmas and finite
certificates supporting the lattice-shadow proof in the research note.
Finite checks do not prove its uniform bound or establish priority.
Uses only Python's standard library. Run from the repository root:
    python3 -B scripts/run_structural_complexity_checks.py
"""

from fractions import Fraction
from heapq import heappop, heappush
from hashlib import sha256
from itertools import combinations, permutations, product
from math import comb, factorial, gcd, isqrt, lcm, prod
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


def minimum_product_rank_controls():
    """Minimum bilinear rank, pole groups and the loss of fan data in box ranks."""
    fixtures, stabilizers = [], []

    def rank(rows, characteristic):
        columns = {(i,): row for i, row in enumerate(rows)}
        return len(greedy(columns, columns, (1,), characteristic))

    def reduce_scalar(value, characteristic):
        value = Fraction(value)
        return (value.numerator * pow(value.denominator, -1, characteristic) % characteristic
                if characteristic else value)

    for n in range(3, 8):
        ground, states = set(hyperbolic_cross(n, 3)), finite_staircases(n, 3)
        extended = ground | set(product(range(3), repeat=3))
        for algebra in ("nilpotent", "split_points"):
            modulus = [0] * n + [1]
            if algebra == "split_points":
                modulus = [1]
                for t in range(1, n + 1):
                    following = [0] * (len(modulus) + 1)
                    for i, value in enumerate(modulus):
                        following[i] -= t * value
                        following[i + 1] += value
                    modulus = following
            characteristics = (0, 2, 3, 11) if algebra == "nilpotent" else (0, 17)
            for characteristic in characteristics:
                # Independently compute the stabilizer of each proper polynomial prefix.
                prefix_checks = []
                for k in range(1, n):
                    powers = [univariate_power_mod((0, 1), i, modulus) for i in range(2 * k - 1)]
                    equations = [tuple(powers[i + j][q] for i in range(k))
                                 for j in range(k) for q in range(k, n)]
                    dimension = k - rank(equations, characteristic)
                    require(dimension == 1, "every proper polynomial prefix has scalar stabilizer")
                    prefix_checks.append({"prefix_dimension": k, "stabilizer_dimension": dimension})
                stabilizers.append({"n": n, "algebra": algebra, "characteristic": characteristic,
                                    "monic_modulus_low_to_high": modulus, "prefix_checks": prefix_checks})
                profiles = []
                for groups in ((0, 0, 0), (0, 0, 1), (0, 1, 2)):
                    if characteristic == 2 and groups == (0, 1, 2):
                        continue  # 1/(1+2*t) is constant in this characteristic.
                    if algebra == "nilpotent":
                        t_image = tuple(int(j == 1) for j in range(n))
                        shifted = [tuple(int(j == 1) + c * int(j == 0) for j in range(n))
                                   for c in (1, 2)]
                        inverse = [tuple((-c) ** j for j in range(n)) for c in (1, 2)]
                        images = ([t_image, *shifted] if groups == (0, 0, 0) else
                                  [t_image, shifted[0], inverse[0]] if groups == (0, 0, 1) else
                                  [t_image, *inverse])
                        columns = cyclic_univariate_model(modulus, images, extended)
                        presentation = {"coordinate_images_low_to_high": images}
                    else:
                        points = [(t, t + 1, t + 2) if groups == (0, 0, 0) else
                                  (t, t + 1, Fraction(1, 1 + t)) if groups == (0, 0, 1) else
                                  (t, Fraction(1, 1 + t), Fraction(1, 1 + 2 * t))
                                  for t in range(1, n + 1)]
                        columns = {a: tuple(point[0] ** a[0] * point[1] ** a[1] * point[2] ** a[2]
                                            for point in points)
                                   for a in extended}
                        presentation = {"rational_point_coordinates": [[str(x) for x in p] for p in points]}
                    columns = {a: tuple(reduce_scalar(x, characteristic) for x in values)
                               for a, values in columns.items()}
                    unit = columns[(0, 0, 0)]
                    axes = [tuple(int(i == j) for j in range(3)) for i in range(3)]
                    require(all(rank([columns[tuple(k * x for x in axis)] for k in range(n)], characteristic) == n
                                for axis in axes), "each coordinate generates the entire length-n algebra")
                    relations = ([(-1, -1, 1, 0), (-2, -1, 1, 0)] if groups == (0, 0, 0) else
                                 [(-1, -1, 1, 0), (-1, 0, 1, 1)] if groups == (0, 0, 1) else
                                 [(-1, 0, 1, 1), (-1, 0, 1, 2)])
                    grading_normals = []
                    for i, coefficients in enumerate(relations, 1):
                        pair = tuple(a + b for a, b in zip(axes[0], axes[i]))
                        monomials = [(0, 0, 0), axes[0], axes[i], pair]
                        rows = [unit, columns[axes[0]], columns[axes[i]], columns[pair]]
                        require(rank(rows, characteristic) == 3, "minimum bilinear product rank")
                        require(all(reduce_scalar(sum(c * row[j] for c, row in zip(coefficients, rows)),
                                                  characteristic) == 0 for j in range(n)),
                                "displayed bilinear relation holds in the actual quotient")
                        a, b, c, d = coefficients
                        require(reduce_scalar(a * d - b * c, characteristic) != 0,
                                "the fractional-linear coefficient matrix is nonsingular")
                        support = [a for a, c in zip(monomials, coefficients) if reduce_scalar(c, characteristic)]
                        require(rank([columns[a] for a in support], characteristic) == len(support) - 1 and
                                all(rank([columns[a] for a in support if a != omitted], characteristic) ==
                                    len(support) - 1 for omitted in support),
                                "each displayed support is a minimal monomial relation")
                        grading_normals.extend(tuple(a[j] - support[0][j] for j in range(3)) for a in support[1:])
                    require(rank(grading_normals, 0) == 3,
                            "minimal relations force every integer grading to be trivial")
                    affine_ranks = [rank([unit, columns[axes[i]], columns[axes[j]]], characteristic)
                                    for i, j in combinations(range(3), 2)]
                    require(affine_ranks == [2 if groups[i] == groups[j] else 3
                                            for i, j in combinations(range(3), 2)],
                            "three-monomial flags detect the pole groups")
                    records, coherent = [], 0
                    for state in states:
                        lengths = [sum(a[i] > 0 for a in state) for i in range(3)]
                        active_groups = [groups[i] for i in range(3) if lengths[i]]
                        predicted = (all(sum(x > 0 for x in a) <= 1 for a in state) and
                                     len(active_groups) == len(set(active_groups)))
                        minor = reduce_scalar(determinant([columns[a] for a in sorted(state)]), characteristic)
                        actual_rank = rank([columns[a] for a in sorted(state)], characteristic)
                        require(bool(minor) == (actual_rank == n) == predicted,
                                "independent determinants, Gaussian ranks and pole-group criterion agree")
                        record = {"staircase": sorted(state), "basis_determinant": str(minor),
                                  "rank": actual_rank}
                        if predicted:
                            rational_weights = [Fraction(2, 2 * k + 1) if k else Fraction(4) for k in lengths]
                            threshold = lcm(*(w.denominator for w in rational_weights))
                            weight = tuple(int(w * threshold) for w in rational_weights)
                            require(max(dot(weight, a) for a in state) < threshold and
                                    all((lengths[i] + 1) * weight[i] > threshold for i in range(3)),
                                    "strict rational cutoff separates the prescribed pure chains")
                            require(greedy(ground, columns, weight, characteristic) == state,
                                    "full hyperbolic-cross Gaussian selection checks each claimed initial ideal")
                            record.update({"positive_integer_weights": weight, "cutoff": threshold})
                            coherent += 1
                        records.append(record)
                    expected = 3 if groups == (0, 0, 0) else 2 * n - 1 if groups == (0, 0, 1) else comb(n + 1, 2)
                    require(coherent == expected, "complete three-variable pole-group count")
                    boxes = []
                    for sides in product(range(1, 4), repeat=3):
                        members = set(product(*(range(a) for a in sides)))
                        actual = rank([columns[a] for a in sorted(members)], characteristic)
                        require(actual == min(n, 1 + sum(a - 1 for a in sides)),
                                "all tested boxes attain the same minimum additive rank")
                        boxes.append({"side_lengths": sides, "rank": actual})
                    profiles.append(boxes)
                    fixtures.append({"n": n, "algebra": algebra, "characteristic": characteristic,
                                     "pole_group_labels": groups, "coordinate_elimination_degrees": [n] * 3,
                                     "homogeneity_dimension": 0,
                                     **presentation, "bilinear_relation_coefficients_A_B_C_D": relations,
                                     "pair_affine_flag_ranks": affine_ranks,
                                     "all_staircase_certificates": records, "initial_ideal_count": coherent,
                                     "box_rank_checks": boxes})
                require(all(profile == profiles[0] for profile in profiles),
                        "identical rectangular rank profiles coexist with different full fans")
    inverse_controls = []
    for n in (5, 7, 9):
        ground = set(hyperbolic_cross(n, 2)) | set(product(range(4), repeat=2))
        images = [tuple(int(j == 1) for j in range(n)),
                  tuple(int(j in (1, 2)) for j in range(n))]
        columns = cyclic_univariate_model((0,) * n + (1,), images, ground)
        for characteristic in (0, 2, 3):
            require(all(rank([columns[tuple(k * int(i == j) for j in range(2))]
                              for k in range(n)], characteristic) == n for i in range(2)),
                    "both coordinates of the inverse-rank negative control are cyclic")
            bilinear = rank([columns[a] for a in product(range(2), repeat=2)], characteristic)
            require(bilinear == 4, "the perturbation is not fractional-linear")
            boxes = []
            for sides in product(range(2, 5), repeat=2):
                actual = rank([columns[a] for a in product(*(range(s) for s in sides))], characteristic)
                require(actual >= min(n, sum(sides)),
                        "a non-fractional-linear cyclic pair has the inverse-rank extra unit")
                boxes.append({"side_lengths": sides, "rank": actual,
                              "additive_lower_bound": min(n, sum(sides) - 1)})
                if n == 5 and sides == (3, 3):
                    require(actual == sum(sides) - 1 == n,
                            "saturated minimum rank does not imply a bilinear relation")
            inverse_controls.append({"n": n, "characteristic": characteristic,
                                     "monic_modulus_low_to_high": [0] * n + [1],
                                     "coordinate_images_low_to_high": images,
                                     "bilinear_rank": bilinear, "box_rank_checks": boxes})
    return {"scope": "Uniform proof covers every field, n>=3, all cyclic-coordinate algebras satisfying the minimum bilinear ranks, including nonreduced quotients. One unsaturated minimum-growth box with every side at least two forces this class. Bounded fixtures test reduced and nilpotent presentations.",
            "uniform_count": "If the d fractional-linear coordinates have pole-group sizes m_j, then g(I)=[z^(n-1)] product_j(1+m_j*z/(1-z)).",
            "uniform_box_rank": "Every anchored monomial box of side lengths a_i has rank min(n,1+sum_i(a_i-1)); every proper box span has scalar stabilizer in the algebra.",
            "information_loss": "The same complete box-rank profile and stabilizer dimensions allow g=d or g=binom(n+d-2,d-1). Three-monomial flags retain pole information lost by box ranks.",
            "priority": "Uses elementary cyclic multiplication, fractional-linear algebra and partial fractions. The coding Kneser theorem and critical-pair theory are prior art. No unrestricted asymptotic improvement or publication novelty is asserted.",
            "fixtures": fixtures, "polynomial_prefix_stabilizer_controls": stabilizers,
            "unsaturated_inverse_negative_controls": inverse_controls}


def product_excess_and_extension_controls():
    """Rational reconstruction with an exceptional ideal, and bounded-border lifts."""
    reconstruction, threshold_counters, extensions, degree_boundary = [], [], [], []

    def scalar(value, characteristic):
        value = Fraction(value)
        return (value.numerator * pow(value.denominator, -1, characteristic) % characteristic
                if characteristic else value)

    def rank(rows, characteristic):
        columns = {(i,): tuple(scalar(x, characteristic) for x in row) for i, row in enumerate(rows)}
        return len(greedy(columns, columns, (1,), characteristic))

    def trim(values, characteristic):
        values = [scalar(x, characteristic) for x in values]
        while values and not values[-1]:
            values.pop()
        return values

    def add(left, right, characteristic, sign=1):
        values = [0] * max(len(left), len(right))
        for i, value in enumerate(left):
            values[i] += value
        for i, value in enumerate(right):
            values[i] += sign * value
        return trim(values, characteristic)

    def multiply(left, right, characteristic):
        values = [0] * max(0, len(left) + len(right) - 1)
        for i, x in enumerate(left):
            for j, y in enumerate(right):
                values[i + j] += x * y
        return trim(values, characteristic)

    def divide(left, right, characteristic):
        left, right = trim(left, characteristic), trim(right, characteristic)
        require(right, "polynomial division by a nonzero polynomial")
        quotient = [0] * max(0, len(left) - len(right) + 1)
        inverse = scalar(1 / Fraction(right[-1]), characteristic) if not characteristic else pow(right[-1], -1, characteristic)
        while len(left) >= len(right):
            position, coefficient = len(left) - len(right), scalar(left[-1] * inverse, characteristic)
            quotient[position] = coefficient
            left = add(left, [0] * position + [coefficient * x for x in right], characteristic, -1)
        return trim(quotient, characteristic), left

    def polynomial_gcd(left, right, characteristic):
        left, right = trim(left, characteristic), trim(right, characteristic)
        while right:
            left, right = right, divide(left, right, characteristic)[1]
        if left:
            inverse = (pow(left[-1], -1, characteristic) if characteristic else 1 / left[-1])
            left = trim([x * inverse for x in left], characteristic)
        return left

    def kernel(columns, characteristic):
        matrix = [[scalar(x, characteristic) for x in row] for row in zip(*columns)]
        width, pivot_row, pivots = len(columns), 0, []
        for j in range(width):
            pivot = next((i for i in range(pivot_row, len(matrix)) if matrix[i][j]), None)
            if pivot is None:
                continue
            matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
            inverse = (pow(matrix[pivot_row][j], -1, characteristic) if characteristic
                       else 1 / matrix[pivot_row][j])
            matrix[pivot_row] = [scalar(x * inverse, characteristic) for x in matrix[pivot_row]]
            for i in range(len(matrix)):
                if i != pivot_row:
                    coefficient = matrix[i][j]
                    matrix[i] = [scalar(x - coefficient * y, characteristic)
                                 for x, y in zip(matrix[i], matrix[pivot_row])]
            pivots.append(j)
            pivot_row += 1
            if pivot_row == len(matrix):
                break
        result = []
        for free in range(width):
            if free not in pivots:
                vector = [scalar(int(i == free), characteristic) for i in range(width)]
                for row, pivot in zip(matrix, pivots):
                    vector[pivot] = scalar(-row[free], characteristic)
                result.append(vector)
        require(all(all(not scalar(sum(c * column[j] for c, column in zip(v, columns)), characteristic)
                            for j in range(len(columns[0]))) for v in result),
                "computed kernel vectors annihilate the independent residue matrix")
        return result

    for n in (7, 9):
        t_image = [int(j == 1) for j in range(n)]
        cases = [
            ("affine", [int(j in (0, 1)) for j in range(n)], 1, 0),
            ("reciprocal", [(-1) ** j for j in range(n)], 1, 0),
            ("quadratic", [int(j in (1, 2)) for j in range(n)], 2, 0),
            ("quadratic_denominator", [(-1) ** ((j - 1) // 2) if j % 2 else 0 for j in range(n)], 2, 0),
            ("one_socle_error", [int(j in (1, n - 1)) for j in range(n)], 1, 1),
            ("two_dimensional_error", [int(j in (1, n - 2)) for j in range(n)], 1, 2),
            ("quadratic_with_socle_error", [int(j in (1, 2, n - 1)) for j in range(n)], 2, 1),
            ("cubic", [int(j in (1, 3)) for j in range(n)], 3, 0),
        ]
        split_modulus = [1]
        lagrange = [1]
        for node in range(1, n + 1):
            split_modulus = multiply(split_modulus, [-node, 1], 0)
            if node < n:
                lagrange = multiply(lagrange, [-node, 1], 0)
        outlier = add(t_image, [Fraction(n + 1, factorial(n - 1)) * x for x in lagrange], 0)
        for algebra, modulus, entries, characteristics in (
                ("nilpotent", [0] * n + [1], cases, (0, 2, 3)),
                ("split_points", split_modulus, [("one_outlier", outlier, 1, 1)], (0, 37))):
            for label, x_image, expected_degree, expected_error in entries:
                for characteristic in characteristics:
                    for b in range(expected_degree + expected_error + 1, 5):
                        require(2 * b - 2 < n, "the relation eliminants are below the modulus degree")
                        domain = set(hyperbolic_cross(n, 2)) | set(product(range(2), range(b)))
                        columns = cyclic_univariate_model(modulus, [x_image, t_image], domain)
                        rows = [columns[(0, j)] for j in range(b)] + [columns[(1, j)] for j in range(b)]
                        actual = rank(rows, characteristic)
                        require(all(rank([columns[tuple(k * int(i == j) for j in range(2))]
                                          for k in range(n)], characteristic) == n for i in range(2)),
                                "both coordinates of the reconstruction fixture are cyclic")
                        relations = kernel(rows, characteristic)
                        require(len(relations) == 2 * b - actual > 0, "independent kernel dimension")
                        p0, q0 = relations[0][:b], relations[0][b:]
                        common = polynomial_gcd(p0, q0, characteristic)
                        p0, q0 = divide(p0, common, characteristic)[0], divide(q0, common, characteristic)[0]
                        require(q0, "cyclic data admit no nonzero relation in y alone below degree n")
                        scale = pow(q0[-1], -1, characteristic) if characteristic else 1 / q0[-1]
                        p0, q0 = trim([scale * x for x in p0], characteristic), trim([scale * x for x in q0], characteristic)
                        require(all(not add(multiply(p0, vector[b:], characteristic),
                                            multiply(vector[:b], q0, characteristic), characteristic, -1)
                                    for vector in relations), "all low-degree relations have the same rational ratio")
                        product_residue = univariate_product_mod(x_image, q0, modulus)
                        residue = add(p0, product_residue, characteristic)
                        quotient_modulus = polynomial_gcd(modulus, residue, characteristic)
                        error = n - len(quotient_modulus) + 1
                        multiplication = [univariate_product_mod(residue, [0] * j + [1], modulus)
                                          for j in range(n)]
                        require(rank(multiplication, characteristic) == error,
                                "independent multiplication rank equals the exceptional ideal dimension")
                        degree = max(len(p0), len(q0)) - 1
                        require((degree, error) == (expected_degree, expected_error) and
                                actual == b + degree + error,
                                "rational degree plus exceptional dimension gives the exact two-row rank")
                        require(polynomial_gcd(q0, quotient_modulus, characteristic) == [1],
                                "the denominator is a unit in the retained quotient")
                        require(not divide(residue, quotient_modulus, characteristic)[1],
                                "the primitive rational relation holds in the retained quotient")
                        annihilator, remainder = divide(modulus, quotient_modulus, characteristic)
                        require(not remainder, "the exceptional annihilator divides the modulus")
                        module_basis = []
                        for j in range(b - degree - error):
                            multiplier = [0] * j + annihilator
                            p, q = multiply(multiplier, p0, characteristic), multiply(multiplier, q0, characteristic)
                            module_basis.append(p + [0] * (b - len(p)) + q + [0] * (b - len(q)))
                        require(rank(module_basis, characteristic) == len(relations) and
                                all(all(not scalar(sum(c * col[j] for c, col in zip(v, rows)), characteristic)
                                        for j in range(n)) for v in module_basis),
                                "the annihilator multiples reproduce the entire truncated relation module")
                        reconstruction.append({"n": n, "algebra": algebra, "case": label, "characteristic": characteristic,
                                               "b": b, "monic_modulus": [str(x) for x in modulus],
                                               "x_image": [str(x) for x in x_image], "rectangle_rank": actual,
                                               "primitive_P": [str(x) for x in p0], "primitive_Q": [str(x) for x in q0],
                                               "rational_degree": degree, "exceptional_ideal_dimension": error,
                                               "retained_quotient_modulus": [str(x) for x in quotient_modulus],
                                               "annihilator_polynomial": [str(x) for x in annihilator],
                                               "full_kernel_basis": [[str(x) for x in row] for row in relations]})

    for characteristic in (0, 2, 3):
        n, b = 6, 4
        image = [int(j in (1, 2, 3)) for j in range(n)]
        modulus = [0] * n + [1]
        columns = cyclic_univariate_model(modulus, [image, [0, 1, 0, 0, 0, 0]],
                                          set(hyperbolic_cross(n, 2)) | set(product(range(2), range(b))))
        rows = [columns[(0, j)] for j in range(b)] + [columns[(1, j)] for j in range(b)]
        relations = kernel(rows, characteristic)
        eliminant = add(multiply(relations[0][:b], relations[1][b:], characteristic),
                        multiply(relations[1][:b], relations[0][b:], characteristic), characteristic, -1)
        require(eliminant and not divide(eliminant, modulus, characteristic)[1],
                "at the degree threshold a nonzero eliminant can be a multiple of the modulus")
        threshold_counters.append({"n": n, "b": b, "characteristic": characteristic,
                                   "rank": rank(rows, characteristic), "nonzero_eliminant": [str(x) for x in eliminant],
                                   "kernel_basis": [[str(x) for x in row] for row in relations]})

    def addable(state):
        border = {tuple(x + int(i == j) for j, x in enumerate(a)) for a in state for i in range(3)} - state
        return {a for a in border if all(tuple(x - int(i == j) for j, x in enumerate(a)) in state
                                        for i, value in enumerate(a) if value)}

    tables = {n: finite_staircases(n, 3) for n in (5, 6, 7)}
    for n, states in tables.items():
        ground = hyperbolic_cross(n, 3)
        indices = {state: i for i, state in enumerate(states)}
        for error in (1, 2):
            m = n - error
            parent_bases = [state for state in finite_staircases(m, 3)
                            if all(sum(x > 0 for x in a) <= 1 for a in state)]
            candidates = set()
            for base in parent_bases:
                level = {base}
                for step in range(error):
                    following = set()
                    for state in level:
                        corners = addable(state)
                        require(len(corners) <= 6 + 2 * step, "bounded minimal-generator growth at each box addition")
                        following.update(state | {corner} for corner in corners)
                    level = following
                candidates.update(level)
            bound = comb(m + 1, 2)
            for step in range(error):
                bound *= 6 + 2 * step
            require(len(candidates) <= bound and candidates <= set(states), "complete extension candidates obey the uniform product bound")
            for algebra, characteristics in (("nilpotent", (0, 3, 11)), ("one_outlier", (0, 37))):
                if algebra == "one_outlier" and error != 1:
                    continue
                if algebra == "nilpotent":
                    images = [tuple(int(j == 1) for j in range(n)),
                              tuple((-1) ** j + int(j == m) for j in range(n)),
                              tuple((-2) ** j + 2 * int(j == m) + int(j == m + 1) for j in range(n))]
                    columns = cyclic_univariate_model([0] * n + [1], images, ground)
                    presentation = {"coordinate_images": images}
                else:
                    points = [(t, Fraction(1, 1 + t), Fraction(1, 1 + 2 * t)) for t in range(1, m + 1)]
                    points.append((n, n + 1, -n - 1))
                    columns = {a: tuple(p[0] ** a[0] * p[1] ** a[1] * p[2] ** a[2] for p in points) for a in ground}
                    presentation = {"rational_points": [[str(x) for x in p] for p in points]}
                for characteristic in characteristics:
                    reduced = {a: tuple(scalar(x, characteristic) for x in row) for a, row in columns.items()}
                    parent_reduced = {a: row[:m] for a, row in reduced.items()}
                    ranks = [rank([reduced[a] for a in sorted(state)], characteristic) for state in states]
                    independent = {states[i] for i, value in enumerate(ranks) if value == n}
                    require(independent <= candidates, "all independent staircases in these fixtures are covered by parent extensions")
                    witnesses = []
                    delta = (1, n + 1, (n + 1) ** 2)
                    scale = 1 + max(dot(delta, a) for a in ground)
                    for base in parent_bases:
                        lengths = [sum(a[i] > 0 for a in base) for i in range(3)]
                        weights = [Fraction(2, 2 * length + 1) if length else Fraction(4) for length in lengths]
                        denominator = lcm(*(w.denominator for w in weights))
                        weight = tuple(scale * int(w * denominator) + delta[i] for i, w in enumerate(weights))
                        small, large = greedy(ground, parent_reduced, weight, characteristic), greedy(ground, reduced, weight, characteristic)
                        require(small == base and small <= large and len(large - small) == error and large in candidates,
                                "coordinate-compatible quotient inclusion gives the predicted staircase extension")
                        witnesses.append({"positive_integer_weight": weight, "larger_staircase_index": indices[frozenset(large)]})
                    extensions.append({"n": n, "parent_length": m, "exceptional_dimension": error,
                                       "algebra": algebra, "characteristic": characteristic, **presentation,
                                       "parent_initial_ideal_count": len(parent_bases),
                                       "uniform_extension_bound": bound, "extension_candidate_indices": sorted(indices[s] for s in candidates),
                                       "all_staircase_ranks": ranks, "independent_staircase_indices": sorted(indices[s] for s in independent),
                                       "quotient_inclusion_witnesses": witnesses,
                                       "scope": "The independent-staircase count is an upper bound, not a claimed exact fan count for the larger quotient."})

    for m in (3, 4, 5):
        n = 2 * m - 2
        basis = [(i, 0, 0) for i in range(m)] + [(0, j, 0) for j in range(1, m - 1)]
        positions = {a: i for i, a in enumerate(basis)}

        def residue(a):
            values = [0] * n
            for j in range(a[2] + 1):
                x, y = a[0], a[1] + j
                if (x and y) or x >= m or y >= m:
                    continue
                monomial = (m - 1, 0, 0) if y == m - 1 else (x, y, 0)
                values[positions[monomial]] += comb(a[2], j)
            return tuple(values)

        matrices = [tuple(residue(tuple(a[j] + int(i == j) for j in range(3))) for a in basis)
                    for i in range(3)]
        require(all(residue(a) == tuple(int(i == j) for j in range(n)) for i, a in enumerate(basis)),
                "the connected-sum reference monomials give a cyclic length-n presentation")
        require(all(matvec(matrices[i], matrices[j][k]) == matvec(matrices[j], matrices[i][k])
                    for i, j in combinations(range(3), 2) for k in range(n)),
                "the connected-sum multiplication matrices commute")
        unit = residue((0, 0, 0))
        for characteristic in (0, 2, 3):
            degrees = [rank([residue(tuple(k * int(i == j) for j in range(3))) for k in range(n)],
                            characteristic) for i in range(3)]
            bilinear = [rank([unit, residue((1, 0, 0)), residue(axis),
                              residue(tuple(x + y for x, y in zip((1, 0, 0), axis)))], characteristic)
                        for axis in ((0, 1, 0), (0, 0, 1))]
            require(degrees == [m] * 3 and bilinear == [3, 3] and 2 * m == n + 2,
                    "the strict high-degree threshold cannot be weakened to equality")
            degree_boundary.append({"m": m, "n": n, "characteristic": characteristic,
                                    "coordinate_elimination_degrees": degrees, "bilinear_ranks": bilinear,
                                    "relations": ["x*y=0", "x^m=y^m=0", "x^(m-1)=y^(m-1)", "z=y+1"],
                                    "reference_basis": basis, "multiplication_matrices_by_columns": matrices})

    return {"scope": "All-field structural proofs are in the note. Reconstruction uses cyclic coordinates and 2*b-2<n. Extension counting requires coordinate-compatible containment and a bound on all smaller initial-ideal generator counts.",
            "rank_tradeoff": "For the nonzero two-by-b relation kernel, rank=b+rational_degree+exceptional_ideal_dimension; these last two quantities sum to at most excess+1 for a larger near-minimum unsaturated rectangle.",
            "extension_bound": "If I is contained in J, lengths differ by e and every initial ideal of J has at most q minimal generators, then g(I)<=g(J)*product_(j=0)^(e-1)(q+(d-1)*j).",
            "priority": "Rational reconstruction, relation modules and minimal monomial generators are established tools. These consequences do not prove a new unrestricted asymptotic bound or certified publication novelty.",
            "reconstruction_fixtures": reconstruction, "degree_threshold_counters": threshold_counters,
            "staircase_tables": {str(n): [sorted(s) for s in states] for n, states in tables.items()},
            "bounded_extension_fixtures": extensions, "high_degree_boundary_controls": degree_boundary}


def quadratic_pole_flag_controls():
    """Direct quotient ranks versus pole flags, with complete small 3D fans."""
    def scalar(value, characteristic):
        value = Fraction(value)
        return (value.numerator * pow(value.denominator, -1, characteristic) % characteristic
                if characteristic else value)

    def eliminate(columns, characteristic, invert=False):
        if not columns:
            return 0, []
        rows, size = len(columns[0]), len(columns)
        require(not invert or rows == size, "only invert a square residue matrix")
        matrix = [[scalar(columns[j][i], characteristic) for j in range(size)] +
                  ([scalar(i == j, characteristic) for j in range(rows)] if invert else [])
                  for i in range(rows)]
        pivot_row = 0
        for j in range(size):
            pivot = next((i for i in range(pivot_row, rows) if matrix[i][j]), None)
            if pivot is None:
                continue
            matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
            scale = (pow(matrix[pivot_row][j], -1, characteristic) if characteristic
                     else 1 / matrix[pivot_row][j])
            matrix[pivot_row] = [scalar(x * scale, characteristic) for x in matrix[pivot_row]]
            for i in range(rows):
                if i != pivot_row and matrix[i][j]:
                    scale = matrix[i][j]
                    matrix[i] = [scalar(x - scale * y, characteristic)
                                 for x, y in zip(matrix[i], matrix[pivot_row])]
            pivot_row += 1
            if pivot_row == rows:
                break
        inverse = [row[size:] for row in matrix] if invert and pivot_row == rows else None
        return pivot_row, inverse

    def multiply(left, right):
        values = [0] * (len(left) + len(right) - 1)
        for i, x in enumerate(left):
            for j, y in enumerate(right):
                values[i + j] += x * y
        return values

    def inverse_linear(modulus, node):
        quotient, carry = [0] * (len(modulus) - 1), modulus[-1]
        for i in range(len(quotient) - 1, -1, -1):
            quotient[i], carry = carry, modulus[i] + node * carry
        require(carry != 0, "the requested pole avoids the divisor")
        return tuple(-Fraction(x, carry) for x in quotient)

    def compositions(total, size):
        if not size:
            if total == 0:
                yield ()
        elif size == 1:
            yield (total,)
        else:
            for first in range(total + 1):
                for rest in compositions(total - first, size - 1):
                    yield (first,) + rest

    def patterns(n, d, chosen):
        other = [i for i in range(d) if i not in (1, chosen)]
        for a in range(1, n + 1):
            for b1 in range((n - a) // 2 + 1):
                if b1 and a < 2:
                    continue
                for gap in range(n - a - 2 * b1 + 1):
                    for chains in compositions(n - a - 2 * b1 - gap, len(other)):
                        state = {(0,) * d}
                        state.update(tuple(k * (i == chosen) for i in range(d)) for k in range(1, a))
                        state.update(tuple(k * (i == 1) for i in range(d)) for k in range(1, b1 + gap + 1))
                        state.update(tuple(k * (i == 1) + (i == chosen) for i in range(d))
                                     for k in range(1, b1 + 1))
                        for axis, length in zip(other, chains):
                            state.update(tuple(k * (i == axis) for i in range(d))
                                         for k in range(1, length + 1))
                        require(len(state) == n, "the pole pattern has the stated length")
                        yield (a, b1 + gap, b1, chains), frozenset(state)

    fixtures, tables, prefix_counters = [], {}, []
    models = [(3, 5, 0, "roots"), (3, 7, 0, "roots"),
              (3, 7, 0, "collision"), (3, 9, 0, "collision"),
              (3, 5, 2, "trinomial"), (3, 7, 2, "trinomial"),
              (3, 5, 2, "delayed_odd_gap"), (3, 7, 3, "trinomial"),
              (4, 6, 0, "four_coordinates"), (4, 6, 3, "four_coordinates"),
              (2, 7, 0, "roots"), (2, 7, 2, "trinomial")]
    for d, n, characteristic, kind in models:
        if kind == "roots":
            modulus = [-1] + [0] * (n - 1) + [1]
        elif kind == "collision":
            modulus = [comb(n, j) * (-1) ** (n - j) for j in range(n + 1)]
        elif kind == "delayed_odd_gap":
            modulus = [1, 1, 0, 1, 1, 1]
        elif kind == "four_coordinates":
            modulus = [2, 0, 1] + [0] * (n - 3) + [1]
        else:
            modulus = [1, 1] + [0] * (n - 2) + [1]
        poles = {i: (1 if characteristic == 2 else 2) if d == 3 else i - 1
                 for i in range(2, d)}
        inverse_t = inverse_linear(modulus, 0)
        images = [tuple(int(i == 1) for i in range(n)),
                  univariate_product_mod(inverse_t, inverse_t, modulus)]
        images.extend(inverse_linear(modulus, poles[i]) for i in range(2, d))
        ground = set(hyperbolic_cross(2 * n, d))
        raw_columns = cyclic_univariate_model(modulus, images, ground)
        columns = {a: tuple(scalar(x, characteristic) for x in vector)
                   for a, vector in raw_columns.items()}
        require(eliminate([columns[tuple(k * (i == 0) for i in range(d))]
                           for k in range(n)], characteristic)[0] == n,
                "the cyclic parameter still generates in the specified characteristic")
        simple = [0] + list(range(2, d))
        cached, sectors, flag_ranks = {}, [], {}
        pattern_lookup, parameter_lookup = {}, {}
        for chosen in simple:
            records = list(patterns(n, d, chosen))
            table_key = str((n, d, chosen))
            table = [{"parameters": parameters, "basis": sorted(state)} for parameters, state in records]
            if table_key in tables:
                require(tables[table_key] == table, "shared pattern tables agree")
            else:
                tables[table_key] = table
                allowed = lambda a: (sum(a[i] > 0 for i in simple) <= 1 and
                                     (not a[1] or (all(a[i] == 0 for i in simple if i != chosen)
                                                   and a[chosen] < 2)))
                exhaustive = {state for state in finite_staircases(n, d)
                              if all(allowed(a) for a in state)}
                require(exhaustive == {state for _, state in records},
                        "the pole patterns enumerate every downset in the parent normal domain")
            sector_records, minimum_gaps = [], {}
            for index, (parameters, state) in enumerate(records):
                a, b0, b1, chains = parameters
                gap = b0 - b1
                other = [i for i in range(d) if i not in (1, chosen)]
                chain = dict(zip(other, chains))
                if state not in cached:
                    ordered = sorted(state)
                    rank_value, inverse = eliminate([columns[p] for p in ordered], characteristic, True)
                    normals = set()
                    if inverse is not None:
                        border = {tuple(x + (i == j) for j, x in enumerate(p))
                                  for p in state for i in range(d)} - state
                        require(border <= ground, "the independent border columns lie in the larger cross")
                        for b in border:
                            coefficients = [scalar(dot(row, columns[b]), characteristic) for row in inverse]
                            for p, value in zip(ordered, coefficients):
                                if value:
                                    normal = tuple(x - y for x, y in zip(b, p))
                                    require(max(normal) > 0, "a border comparison has a positive component")
                                    if min(normal) < 0:
                                        divisor = gcd(*normal)
                                        normals.add(tuple(x // divisor for x in normal))
                    cached[state] = rank_value, sorted(normals)
                rank_value, normals = cached[state]
                if a >= 2:
                    flag = (chosen, gap, chains) if chosen == 0 else (
                        chosen, gap, a, tuple(chain[i] for i in range(2, d) if i != chosen))
                    if flag not in flag_ranks:
                        if chosen == 0:
                            common = [1]
                            for axis, length in chain.items():
                                for _ in range(length):
                                    common = multiply(common, [-poles[axis], 1])
                            independent_columns = [univariate_power_mod((0, 1), 2 * gap + j, modulus)
                                                   for j in range(n - gap)]
                            independent_columns += [univariate_product_mod(
                                univariate_power_mod((0, 1), 2 * k, modulus), common, modulus)
                                for k in range(gap)]
                            reduced_rank = eliminate(independent_columns, characteristic)[0]
                        else:
                            reciprocal = tuple(Fraction(x, modulus[0]) for x in reversed(modulus))
                            common = [1]
                            for _ in range(a - 1):
                                common = multiply(common, [1, -poles[chosen]])
                            for axis in range(2, d):
                                if axis != chosen:
                                    for _ in range(chain[axis]):
                                        common = multiply(common, [1, -poles[axis]])
                            exponent = n - a - gap - sum(chain[i] for i in range(2, d) if i != chosen)
                            top_columns = [univariate_product_mod(
                                univariate_power_mod((0, 1), exponent + 2 * k, reciprocal),
                                common, reciprocal)[n - gap:] for k in range(1, gap + 1)]
                            reduced_rank = n - gap + eliminate(top_columns, characteristic)[0]
                        flag_ranks[flag] = reduced_rank
                    require(flag_ranks[flag] == rank_value,
                            "the independent pole flag gives the complete direct quotient rank")
                    if rank_value == n and chosen != 0:
                        key = (a, tuple(chain[i] for i in range(2, d) if i != chosen), gap % 2)
                        minimum_gaps[key] = min(minimum_gaps.get(key, gap), gap)
                row = {"pattern_index": index, "rank": rank_value, "status": "dependent"}
                if rank_value == n and d == 3:
                    sector_normals = [tuple(int(i == axis) - int(i == chosen) for i in range(d))
                                      for axis in simple if axis != chosen]
                    polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
                    for normal in normals + sector_normals:
                        polygon = clip(polygon, normal, 1)
                        if not polygon:
                            break
                    row.update({"status": "independent_but_incoherent", "strict_border_normals": normals,
                                "strict_sector_normals": sector_normals, "closed_region_vertices": []})
                    if polygon:
                        center = [sum(p[i] for p in polygon) / len(polygon) for i in range(3)]
                        scale = lcm(*(x.denominator for x in center))
                        weight = tuple(int(x * scale) for x in center)
                        require(min(weight) > 0 and all(dot(weight, v) > 0 for v in normals + sector_normals),
                                "every complete-fan certificate has a strict positive integer weight")
                        require(greedy(hyperbolic_cross(n, d), columns, weight, characteristic) == set(state),
                                "independent cross-wide greedy elimination verifies the complete-fan witness")
                        row.update({"status": "coherent", "weight": weight,
                                    "closed_region_vertices": [[str(x) for x in p] for p in polygon]})
                elif rank_value == n:
                    row["status"] = "independent_not_fan_classified"
                sector_records.append(row)
                pattern_lookup[chosen, state] = index
                parameter_lookup[chosen, parameters] = state
            sectors.append({"chosen_axis": chosen, "pattern_table": table_key,
                            "minimum_admissible_gaps": [{"flag": key, "gap": value}
                                                        for key, value in sorted(minimum_gaps.items())],
                            "all_pattern_checks": sector_records})

        selected = []
        if d != 3:
            tie = tuple((2 * n + 1) ** i for i in range(d))
            scale = 1 + max(dot(tie, p) for p in ground)
            for ordering in permutations((1, 3) if d == 2 else (1, 3, 7, 15)):
                weight = tuple(scale * x + y for x, y in zip(ordering, tie))
                require(len({dot(weight, p) for p in ground}) == len(ground), "selected weights have no cross ties")
                chosen = min(simple, key=lambda i: weight[i])
                state = frozenset(greedy(hyperbolic_cross(n, d), columns, weight, characteristic))
                require((chosen, state) in pattern_lookup and cached[state][0] == n,
                        "selected global greedy bases lie in the predicted complete parent domains")
                selected.append({"chosen_axis": chosen, "pattern_index": pattern_lookup[chosen, state],
                                 "weight": weight})

        count_checks, monotone, mechanism_checks = 0, {}, []
        witnesses = [dict(row, chosen_axis=sector["chosen_axis"])
                     for sector in sectors for row in sector["all_pattern_checks"] if row["status"] == "coherent"] + selected
        for row in witnesses:
            chosen, weight = row["chosen_axis"], row["weight"]
            parameters, state = list(patterns(n, d, chosen))[row["pattern_index"]]
            a, b0, b1, chains = parameters
            gap, h = b0 - b1, Fraction(b0 + b1, 2)
            if a < 2:
                continue
            if chosen == 0:
                monotone.setdefault((chains, gap % 2), []).append((Fraction(weight[0], weight[1]), h, gap))
            else:
                chain = dict(zip([i for i in range(d) if i not in (1, chosen)], chains))
                key = (a, tuple(chain[i] for i in range(2, d) if i != chosen), gap % 2)
                sector = next(s for s in sectors if s["chosen_axis"] == chosen)
                least = next(r["gap"] for r in sector["minimum_admissible_gaps"] if r["flag"] == key)
                product_bound = max(chain[0] - 1, 0) * max(gap - least - 1, 0)
                require(product_bound < n, "minimum feasible parity gaps satisfy the characteristic-free hyperbolic bound")
                comparisons = []
                competitors = []
                if gap > least:
                    new_b1 = (b0 + b1 - least) // 2
                    competitors.append(("least_gap", (a, new_b1 + least, new_b1, chains)))
                if chain[0] >= 2:
                    axes = [i for i in range(d) if i not in (1, chosen)]
                    new_chains = tuple(value - 2 * (axis == 0) for axis, value in zip(axes, chains))
                    competitors.append(("common_shift", (a, b0 + 1, b1 + 1, new_chains)))
                for name, parameters_other in competitors:
                    competitor = parameter_lookup[chosen, parameters_other]
                    require(cached[competitor][0] == n, "the counting comparison is an actual independent quotient basis")
                    cost_difference = sum(dot(weight, p) for p in competitor) - sum(dot(weight, p) for p in state)
                    require(cost_difference > 0, "strict global coherence verifies the comparison cost inequality")
                    comparisons.append({"comparison": name, "pattern_index": pattern_lookup[chosen, competitor],
                                        "rank": n, "cost_difference": cost_difference})
                mechanism_checks.append({"chosen_axis": chosen, "pattern_index": row["pattern_index"],
                                         "weight": weight, "minimum_admissible_gap": least,
                                         "hyperbolic_product": product_bound, "comparison_basis_checks": comparisons})
                count_checks += 1
        for values in monotone.values():
            ordered = sorted(set(values))
            require(all(h1 <= h2 and l1 <= l2 for (_, h1, l1), (_, h2, l2) in zip(ordered, ordered[1:])),
                    "coherent first-axis minimizers form the predicted coordinatewise monotone chain")
        distinct = {tuple(map(tuple, tables[s["pattern_table"]][r["pattern_index"]]["basis"]))
                    for s in sectors for r in s["all_pattern_checks"] if r["status"] == "coherent"}
        fixtures.append({"dimension": d, "n": n, "characteristic": characteristic, "kind": kind,
                         "modulus": modulus, "poles_by_axis": poles,
                         "coordinate_images": [[str(v) for v in image] for image in images],
                         "sectors": sectors, "independent_flag_rank_checks": len(flag_ranks),
                         "complete_positive_fan": d == 3,
                         "initial_ideal_count": len(distinct) if d == 3 else None,
                         "selected_weight_checks": selected, "hyperbolic_witness_checks": count_checks,
                         "counting_mechanism_witnesses": mechanism_checks,
                         "monotone_chain_checks": [{"flag": key, "ordered_weights_h_and_gaps":
                                                    [[str(r), str(h), gap] for r, h, gap in sorted(set(values))]}
                                                   for key, values in sorted(monotone.items())]})

    for m in (2, 3, 5):
        n = 2 * m + 1
        weight = (2 * n * m + 1, 2 * n)
        parent = {(i, 0) for i in range(3 * m + 2)}
        parent |= {(i, j) for i in (0, 1) for j in range(1, 3 * m + 2)}
        prefix = sorted(parent, key=lambda a: dot(weight, a))[:3 * m + 1]
        expected_before = {(0, j) for j in range(2 * m)} | {(1, j) for j in range(m)}
        require(set(prefix[:-1]) == expected_before and prefix[-1] == (0, 2 * m),
                "the first rank completion occurs at the asserted parent-normal position")
        residues = [(a - 2 * b) % n for a, b in prefix]
        require(len(set(residues[:-1])) == n - 1 and len(set(residues)) == n,
                "independent standard-unit residues prove the rank delay over every field")
        require(len({(-2 * j) % n for j in range(n)}) == n,
                "both coordinate images are cyclic over every field")
        prefix_counters.append({"m": m, "n": n, "modulus": "t^n-1", "images": ["t", "t^(-2)"],
                                "integer_weight": weight, "parent_normal_prefix": prefix,
                                "standard_unit_residue_indices": residues,
                                "rank_before_completion": n - 1, "completion_position": 3 * m + 1,
                                "delay_beyond_n": m})

    delayed = next(f for f in fixtures if f["kind"] == "delayed_odd_gap")
    checks = next(s for s in delayed["sectors"] if s["chosen_axis"] == 2)
    least = next(r["gap"] for r in checks["minimum_admissible_gaps"] if r["flag"] == (2, (), 1))
    require(least == 3, "a real characteristic-two quotient disproves comparison with gap one alone")
    encoded_tables = {key: [row["parameters"] for row in table] for key, table in tables.items()}
    for fixture in fixtures:
        for sector in fixture["sectors"]:
            for row in sector["all_pattern_checks"]:
                for key in ("strict_border_normals", "strict_sector_normals", "closed_region_vertices"):
                    if key in row:
                        row[key] = [",".join(map(str, vector)) for vector in row[key]]
    return {"uniform_conclusion": "For every field and every fixed d>=3, all length-n finite quotients with coordinates t,t^(-2),(t-lambda_3)^(-1),...,(t-lambda_d)^(-1), at distinct nonzero poles avoiding the divisor, have O_d(n^(d-1)*log(n+1)) monomial initial ideals. For d=2 the bound is O(n).",
            "scope": "The proof uses minimum admissible parity gaps and is uniform in the divisor, poles and characteristic. It remains restricted and does not improve the unrestricted asymptotic bound.",
            "encoding": "A pattern table key is (n,d,chosen_axis). Each row is (a,b0,b1,other_axis_chain_lengths), with other axes increasing and excluding y and the chosen axis; QP2 reconstructs the exact basis. Border normals and polygon vertices are comma-separated exact coordinate strings.",
            "pattern_tables": encoded_tables, "fixtures": fixtures, "prefix_delay_countercontrols": prefix_counters}


def staircase_moment_cube_controls():
    """An exact staircase-moment relaxation; no common ideal is constructed."""
    row_tables, profiles, unions, ball_controls, prefix_counters = [], {}, [], [], []

    def row_profile(radius, value):
        key = radius, value
        if key in profiles:
            return profiles[key]
        pairs = {}
        for a, b in product(range(radius + 1), repeat=2):
            pairs.setdefault(a * a + b * b, (a, b))
        first = next((s for s in sorted(pairs) if value - s in pairs), None)
        require(first is not None, "the bounded four-square witness exists")
        squares = pairs[first] + pairs[value - first]
        require(sum(x * x for x in squares) == value <= radius * radius,
                "the exact four-square certificate has the requested norm")
        variation = [0] * 13
        for j, z in enumerate(squares):
            variation[3 * j + 1:3 * j + 4] = [z, -2 * z, z]
        require(variation[0] == sum(variation) == sum(i * x for i, x in enumerate(variation)) == 0,
                "the independent second-difference patches preserve both linear constraints")
        require(sum(x * x for x in variation) == 6 * value,
                "the disjoint patches have their claimed squared norm")
        slope = 4 * radius + 1
        intercept = 12 * slope + 2 * radius + 13
        baseline = [intercept - slope * i for i in range(13)]
        rows = [x + y for x, y in zip(baseline, variation)]
        require(rows[0] == intercept and min(rows) >= 13 and
                all(x > y for x, y in zip(rows, rows[1:])),
                "the row lengths define a partition with a fixed first row")
        points = {(x, y) for y, length in enumerate(rows) for x in range(length)}
        actual = tuple(sum(p[i] for p in points) for i in range(2))
        base_moment = (sum(comb(x, 2) for x in baseline), sum(i * x for i, x in enumerate(baseline)))
        require(len(points) == sum(baseline) == 338 * radius + 247 and
                actual == (base_moment[0] + 3 * value, base_moment[1]),
                "direct planar enumeration proves the fixed-cardinality moment interval")
        profiles[key] = tuple(rows)
        row_tables.append({"radius": radius, "value": value,
                           "four_squares": ",".join(map(str, squares)),
                           "row_lengths": ",".join(map(str, rows)),
                           "cardinality": len(points), "moments": actual})
        return profiles[key]

    def construct(d, radius, values, padding):
        state = set()
        for axis, value in enumerate(values):
            rows = row_profile(radius, value)
            for j, length in enumerate(rows):
                for k in range(length):
                    a = [0] * d
                    a[axis], a[(axis + 1) % d] = k, j
                    state.add(tuple(a))
        intercept = 50 * radius + 25
        additions = {tuple(k if j == 0 else 0 for j in range(d))
                     for k in range(intercept, intercept + padding)}
        require(not state & additions, "fixed-first-row padding adds exactly the stated number of cells")
        state |= additions
        require(all(tuple(x - (i == j) for j, x in enumerate(a)) in state
                    for a in state for i in range(d) if a[i]),
                "every expanded union and padding is a genuine monomial downset")
        slope = 4 * radius + 1
        baseline = [intercept - slope * i for i in range(13)]
        planar_moments = sum(comb(x, 2) for x in baseline) + sum(i * x for i, x in enumerate(baseline))
        correction = 13 * 78 if d == 2 else 78
        base = [planar_moments - correction] * d
        base[0] += padding * intercept + padding * (padding - 1) // 2
        expected_n = (676 * radius + 325 if d == 2 else 338 * d * radius + 234 * d + 1) + padding
        moments = tuple(sum(a[i] for a in state) for i in range(d))
        require(len(state) == expected_n and moments == tuple(x + 3 * v for x, v in zip(base, values)),
                "direct union enumeration proves the translated moment-cube identity")
        return expected_n, tuple(base), moments

    for radius in range(1, 7):
        for value in range(radius * radius + 1):
            row_profile(radius, value)
    for d in range(2, 6):
        for radius in (1, 2, 3):
            period = 676 if d == 2 else 338 * d
            offset = 325 if d == 2 else 234 * d + 1
            tuples = {(0,) * d, (radius * radius,) * d,
                      tuple((i + 1) * radius % (radius * radius + 1) for i in range(d))}
            tuples.update(tuple(radius * radius * (i == j) for j in range(d)) for i in range(d))
            for values in sorted(tuples):
                for padding in (0, 1, period - 1):
                    n, base, actual = construct(d, radius, values, padding)
                    require((n - offset) // period == radius,
                            "bounded padding covers the complete cardinality period")
                    unions.append({"dimension": d, "radius": radius, "padding": padding,
                                   "cube_coordinates": values, "n": n, "base_moment": base,
                                   "actual_moment": actual})

    for radius in (2, 3):
        ball_radius = radius * radius // 2
        points = sorted(p for p in product(range(-ball_radius, ball_radius + 1), repeat=3)
                        if sum(x * x for x in p) <= ball_radius * ball_radius)
        certificates = []
        for p in points:
            if max(p) > 0:
                continue
            polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
            for q in points:
                if q != p:
                    normal = tuple(y - x for x, y in zip(p, q))
                    polygon = clip(polygon, normal, 1)
                    if not polygon:
                        break
            if not polygon:
                continue
            center = [sum(v[i] for v in polygon) / len(polygon) for i in range(3)]
            scale = lcm(*(x.denominator for x in center))
            weight = tuple(int(x * scale) for x in center)
            divisor = gcd(*weight)
            weight = tuple(x // divisor for x in weight)
            margin = min(dot(weight, tuple(y - x for x, y in zip(p, q))) for q in points if q != p)
            require(min(weight) > 0 and margin > 0, "the ball vertex has a strict positive minimum witness")
            values = tuple(x + ball_radius for x in p)
            n, base, actual = construct(3, radius, values, 7)
            require(all(dot(weight, actual) < dot(weight, tuple(base[i] + 3 * (q[i] + ball_radius)
                                                             for i in range(3)))
                        for q in points if q != p),
                    "the expanded staircase is exposed among the entire ball-indexed family")
            certificates.append({"ball_vertex": p, "cube_coordinates": values, "positive_weight": weight,
                                 "strict_integer_margin": 3 * margin, "n": n, "staircase_moment": actual})
        require(certificates, "the finite ball control has positive exposed vertices")
        ball_controls.append({"radius_parameter": radius, "integer_ball_radius": ball_radius,
                              "all_ball_points": [",".join(map(str, p)) for p in points],
                              "positive_minimum_certificates": certificates,
                              "scope": "Exposed only within this explicit family of staircases; not certified initial ideals of a common algebra."})
        if radius == 2:
            def costs(point, weight):
                state = set()
                for axis, value in enumerate(tuple(x + ball_radius for x in point)):
                    for j, length in enumerate(row_profile(radius, value)):
                        for k in range(length):
                            a = [0] * 3
                            a[axis], a[(axis + 1) % 3] = k, j
                            state.add(tuple(a))
                intercept = 50 * radius + 25
                state.update((k, 0, 0) for k in range(intercept, intercept + 7))
                return sorted(dot(weight, a) for a in state)

            for certificate in certificates:
                point, weight = certificate["ball_vertex"], certificate["positive_weight"]
                own = costs(point, weight)
                counter = None
                for other in points:
                    if other == point:
                        continue
                    comparison = costs(other, weight)
                    failure = next((i for i, (x, y) in enumerate(zip(own, comparison)) if x > y), None)
                    if failure is not None:
                        threshold = comparison[failure]
                        own_prefix = sum(x <= threshold for x in own)
                        other_prefix = sum(x <= threshold for x in comparison)
                        require(sum(own) < sum(comparison) and own_prefix < other_prefix,
                                "positive total-moment exposure does not imply the necessary greedy prefix dominance")
                        counter = {"selected_ball_vertex": point, "other_ball_point": other,
                                   "positive_weight": weight, "rank_position_one_based": failure + 1,
                                   "selected_cost_at_position": own[failure],
                                   "other_cost_at_position": threshold,
                                   "selected_prefix_size": own_prefix, "other_prefix_size": other_prefix,
                                   "selected_total_cost": sum(own), "other_total_cost": sum(comparison),
                                   "conditional_conclusion": "If both encoded staircases are quotient bases, the selected one cannot be initial at this weight or any monomial order refining it. Other weights are not excluded."}
                        break
                require(counter is not None, "each saved small exposure has an independently checked prefix obstruction")
                prefix_counters.append(counter)
    return {"uniform_cube": "For every fixed d>=2 and all sufficiently large n, moments of actual n-element downsets contain a translated 3*[0,R^2]^d integer grid, with R=Theta_d(n). Four disjoint (1,-2,1) row patches, Lagrange's theorem and fixed-row padding prove this.",
            "relaxation_extremal_order": "The maximum positive exposed-moment count of an arbitrary family of n-element staircases is Theta_d(n^(2*d*(d-1)/(d+1))). The lower bound uses the classical Barany--Larman integer-ball theorem; the upper uses the existing lattice-polytope bound.",
            "scope": "A sharp obstruction to improving the unrestricted bound using only individual downsets and their moments. No common quotient, representable residue matroid or improved G_d/F_d bound is supplied.",
            "row_encoding": "Rows j=0,...,12 contain x=0,...,ell_j-1. Cyclic arm i uses coordinates i and (i+1) mod d. Padding extends the first pure axis from B=50*R+25 by the stated count.",
            "row_certificates": row_tables, "expanded_union_checks": unions,
            "integer_ball_relaxation_checks": ball_controls,
            "prefix_dominance_countercontrols": prefix_counters}


def exchange_cycle_coherence_controls():
    """Basis-exchange closure, not a new unrestricted initial-ideal bound."""
    contractions, patches, cycles, templates = [], [], [], []
    for prime in (2, 3, 5, 7, 11):
        masks = {"first_only": 0, "second_only": 0, "both": 0}
        for a, b in product(product(range(prime), repeat=2), repeat=2):
            if not (a[0] * b[1] - a[1] * b[0]) % prime:
                continue
            # Contract the common n-2 basis entries and use the old pair as coordinates.
            valid = (a[0] != 0, a[1] != 0)
            require(any(valid), "the first new entry permits at least one old-entry exchange")
            masks["both" if all(valid) else "first_only" if valid[0] else "second_only"] += 1
        require(sum(masks.values()) == (prime ** 2 - 1) * (prime ** 2 - prime),
                "all invertible rank-two contractions were examined")
        contractions.append({"characteristic": prime, "invertible_new_pairs": sum(masks.values()),
                             "permitted_removed_entry_masks": masks})

    for radius in range(1, 9):
        slope, intercept = 4 * radius + 1, 50 * radius + 25
        for z in range(radius):
            states = []
            for value in (z, z + 1):
                variation = [0] * 13
                variation[10:13] = [value, -2 * value, value]
                rows = [intercept - slope * j + variation[j] for j in range(13)]
                require(all(x > y for x, y in zip(rows, rows[1:])), "both patch profiles remain strict partitions")
                states.append({(x, y) for y, length in enumerate(rows) for x in range(length)})
            old, new = sorted(states[0] - states[1]), sorted(states[1] - states[0])
            middle = intercept - 11 * slope - 2 * z
            require(old == [(middle - 2, 11), (middle - 1, 11)] and
                    new == [(intercept - 12 * slope + z, 12), (intercept - 10 * slope + z, 10)],
                    "direct downset subtraction gives the two-for-two patch exchange")
            lower, upper = slope - 3 * z - 2, slope + 3 * z + 2
            tested = 0
            for ratio in sorted({Fraction(1, 2), Fraction(lower) - Fraction(1, 2), Fraction(lower),
                                 Fraction(lower) + Fraction(1, 2), Fraction(slope),
                                 Fraction(upper) - Fraction(1, 2), Fraction(upper),
                                 Fraction(upper) + Fraction(1, 2), Fraction(2 * upper)}):
                require(ratio > 0, "all tested weights are positive")
                costs_old = sorted(Fraction(a) + ratio * b for a, b in old)
                costs_new = sorted(Fraction(a) + ratio * b for a, b in new)
                # Common monomials cancel from every cumulative cost count.
                thresholds = sorted(set(costs_old + costs_new))
                dominance = all(sum(x <= t for x in costs_old) >= sum(x <= t for x in costs_new)
                                for t in thresholds)
                require(dominance == all(x <= y for x, y in zip(costs_old, costs_new)) ==
                        (lower <= ratio <= upper), "exact prefix dominance is precisely the claimed ratio interval")
                tested += 1
            patches.append({"R": radius, "old_patch_parameter": z, "removed": [",".join(map(str, a)) for a in old],
                            "added": [",".join(map(str, a)) for a in new],
                            "necessary_weight_ratio_interval": [lower, upper], "weight_checks": tested})

    for data in (((2, 3), (2, 1)), ((1, 2), (3, 2), (2, 1)),
                 ((1, 3), (2, 1), (2, 1), (2, 1)), ((7, 1),) * 3,
                 ((1, 2), (2, 1)), ((1, 3), (2, 1))):
        d = len(data)
        gains = [Fraction(a, b) for a, b in data]
        gain_product = Fraction(1)
        for value in gains:
            gain_product *= value
        if gain_product <= 1:
            weight = [Fraction(1)]
            for value in gains[:-1]:
                weight.append(weight[-1] * value)
            require(all(weight[(i + 1) % d] >= gains[i] * weight[i] for i in range(d)),
                    "the nonexpanding-cycle control satisfies all weak ratio inequalities")
            cycles.append({"integer_edge_offsets": data, "gain_product": str(gain_product),
                           "weak_ratio_weights": [str(x) for x in weight],
                           "scope": "Feasible weak ratios only; no initial-basis existence is asserted."})
            continue
        step = Fraction(2)
        while step ** (d - 1) >= gain_product:
            step = 1 + (step - 1) / 2
        scaled = [Fraction(1)]
        for value in gains[1:]:
            scaled.append(scaled[-1] * step / value)
        coefficients = [value / data[i][1] for i, value in enumerate(scaled)]
        total = sum(coefficients)
        coefficients = [x / total for x in coefficients]
        worst = tuple(coefficients[(i - 1) % d] * data[(i - 1) % d][1] - coefficients[i] * data[i][0]
                      for i in range(d))
        require(min(coefficients) > 0 and sum(coefficients) == 1 and max(worst) < 0,
                "positive rational convex coefficients give strict coordinatewise moment domination")
        for removals in product((0, 1), repeat=d):
            actual = tuple(worst[i] - coefficients[i] * removals[i] for i in range(d))
            require(all(a <= b < 0 for a, b in zip(actual, worst)),
                    "every possible forced old-entry removal preserves the domination certificate")
        cycles.append({"integer_edge_offsets": data, "gain_product": str(gain_product),
                       "convex_coefficients": [str(x) for x in coefficients],
                       "coordinate_drop_upper_bounds": [str(x) for x in worst],
                       "removal_choices_checked": 2 ** d})

    def state(radius, values, padding):
        result = set()
        for axis, value in enumerate(values):
            variation = [0] * 13
            if value >= 1:
                variation[10:13] = [1, -2, 1]
            if value == 2:
                variation[7:10] = [1, -2, 1]
            require(value in (0, 1, 2), "only the needed canonical four-square profiles are reconstructed")
            for row in range(13):
                length = 50 * radius + 25 - (4 * radius + 1) * row + variation[row]
                for column in range(length):
                    point = [0] * 3
                    point[axis], point[(axis + 1) % 3] = column, row
                    result.add(tuple(point))
        result.update((x, 0, 0) for x in range(50 * radius + 25, 50 * radius + 25 + padding))
        require(all(tuple(x - (i == j) for j, x in enumerate(a)) in result
                    for a in result for i in range(3) if a[i]), "each reconstructed template is a downset")
        require(len(result) == 1014 * radius + 703 + padding, "the padded template has its exact length")
        return result

    for radius in range(2, 7):
        for padding in (0, 7):
            witnesses = []
            lower = 4 * radius - 1
            for kind in ("central", "axis"):
                for axis in range(3):
                    values = [1, 1, 1] if kind == "central" else [0 if i == axis else 2 for i in range(3)]
                    changed = list(values)
                    changed[axis] += 1
                    original, alternative = state(radius, values, padding), state(radius, changed, padding)
                    removed, added = sorted(original - alternative), sorted(alternative - original)
                    require(len(removed) == len(added) == 2, "the whole-arm unions retain exactly the local exchange")
                    u = min(removed, key=lambda a: a[axis])
                    v = max(removed, key=lambda a: a[axis])
                    early = max(added, key=lambda a: a[(axis + 1) % 3])
                    require(tuple(v[i] - u[i] for i in range(3)) == tuple(i == axis for i in range(3)) and
                            tuple(early[i] - u[i] for i in range(3)) ==
                            tuple(-lower * (i == axis) + (i == (axis + 1) % 3) for i in range(3)),
                            "the forced-exchange moment offset is exactly -L*e_i+e_next")
                    require(all(sum(a[i] for a in alternative) - sum(a[i] for a in original) == 3 * (i == axis)
                                for i in range(3)), "the competing two-for-two basis has the higher original moment")
                    witnesses.append({"kind": kind, "axis": axis, "cube_coordinates": values,
                                      "competing_cube_coordinates": changed,
                                      "removed": [",".join(map(str, a)) for a in removed],
                                      "added": [",".join(map(str, a)) for a in added]})
            for bits in product((0, 1), repeat=3):
                drop = tuple(Fraction(1 - lower - bits[i], 3) for i in range(3))
                require(all(x <= -Fraction(lower - 1, 3) for x in drop),
                        "all possible three forced bases have a uniformly dominated average moment")
            templates.append({"R": radius, "padding": padding, "n": 1014 * radius + 703 + padding,
                              "gain": lower, "central_uniform_coordinate_drop": str(Fraction(lower - 1, 3)),
                              "conditional_basis_exchange_witnesses": witnesses})
    return {"general_lemma": "For a common basis D and two-for-two competing bases, contract the common n-2 entries. A new monomial u-alpha*e_i+beta*e_j replaces at least one of u,u+h*e_i. A directed coordinate cycle with product(alpha)>product(beta) therefore forces a convex combination of actual bases whose moment is strictly below mu(D) on every cycle coordinate. D cannot be initial under any global monomial order.",
            "scope": "A self-contained use of established basis exchange and state-polyhedron theory. It is a coherence obstruction, not a characterization, unrestricted fan count, or certified novel theorem.",
            "uniform_template_consequence": "For every R>=2 and any padding delta>=0, if all 33 radius-two-ball-indexed staircases from the prior moment construction are bases of one quotient, none is initial for that quotient. The three axial candidates fail their inward-neighbor ratio condition; the central candidate has a forced average moment at least (4R-2)/3 lower in each coordinate.",
            "limits": "This closes this small-ball encoding at all weights and all lengths in the family. It does not exclude other encodings, large-ball families, or improve G_d/F_d. The supplied templates are conditional basis comparisons, not a claim of global quotient realization from independent contraction matrices.",
            "rank_two_contractions": contractions, "patch_controls": patches,
            "rational_cycle_controls": cycles, "expanded_template_controls": templates}


def algebraic_attention_hotpath_controls():
    """Exact constrained polynomial attention, using an established feature principle."""
    fixtures, pair_checks, derivative_checks, output_checks, output_derivatives = [], 0, 0, 0, 0

    def monomial(point, exponent):
        answer = Fraction(1)
        for x, power in zip(point, exponent):
            answer *= x ** power
        return answer

    def partial(point, exponent, axis):
        if not exponent[axis]:
            return Fraction(0)
        reduced = tuple(x - (i == axis) for i, x in enumerate(exponent))
        return exponent[axis] * monomial(point, reduced)

    for dimension, degree in ((1, 2), (1, 4), (2, 1), (2, 2), (3, 1)):
        zero = (0,) * dimension
        support = [a for a in product(range(degree + 1), repeat=dimension) if a != zero]
        weights = {a: Fraction(1 + sum((i + 1) * x for i, x in enumerate(a)) % 5,
                               2 * len(support)) for a in support}
        coefficients = {zero: Fraction(1)}
        for a, value in weights.items():
            coefficients[a] = coefficients.get(a, 0) + value
        for a, first in weights.items():
            for b, second in weights.items():
                exponent = tuple(x + y for x, y in zip(a, b))
                coefficients[exponent] = coefficients.get(exponent, 0) + first * second / 2
        require(set(coefficients) == set(product(range(2 * degree + 1), repeat=dimension)),
                "the collapsed support is the full doubled box")
        queries = [tuple(Fraction(((i + 2) * (j + 1)) % 7 - 3, 5)
                         for j in range(dimension)) for i in range(5)]
        keys = [tuple(Fraction(((i + 1) * (j + 2)) % 7 - 3, 6)
                      for j in range(dimension)) for i in range(5)]
        values = [(Fraction(i - 2), Fraction((i * i) % 5 - 2), Fraction(1, i + 1))
                  for i in range(5)]

        def direct(q, k):
            score = sum(value * monomial(q, a) * monomial(k, a)
                        for a, value in weights.items())
            return 1 + score + score * score / 2

        def direct_partial(q, k, axis, query_side):
            score = sum(value * monomial(q, a) * monomial(k, a)
                        for a, value in weights.items())
            change = sum(value * (partial(q, a, axis) * monomial(k, a) if query_side
                                  else monomial(q, a) * partial(k, a, axis))
                         for a, value in weights.items())
            return (1 + score) * change

        def collapsed_partial(q, k, axis, query_side):
            return sum(value * (partial(q, b, axis) * monomial(k, b) if query_side
                                else monomial(q, b) * partial(k, b, axis))
                       for b, value in coefficients.items())

        for q in queries:
            for k in keys:
                actual = sum(c * monomial(q, b) * monomial(k, b) for b, c in coefficients.items())
                require(actual == direct(q, k) >= Fraction(1, 2),
                        "expanded dot-product attention and collapsed features agree exactly")
                pair_checks += 1
                for axis in range(dimension):
                    for query_side in (False, True):
                        require(direct_partial(q, k, axis, query_side) ==
                                collapsed_partial(q, k, axis, query_side),
                                "query and key derivatives survive the feature collapse")
                        derivative_checks += 1
                for a in support:
                    score = sum(c * monomial(q, b) * monomial(k, b) for b, c in weights.items())
                    expected = (1 + score) * monomial(q, a) * monomial(k, a)
                    actual = monomial(q, a) * monomial(k, a)
                    actual += sum(c * monomial(q, tuple(x + y for x, y in zip(a, b))) *
                                  monomial(k, tuple(x + y for x, y in zip(a, b)))
                                  for b, c in weights.items())
                    require(actual == expected, "the tied coefficient map retains learned-weight derivatives")
                    derivative_checks += 1

        samples = []
        for causal in (False, True):
            state = {b: [Fraction(0)] * 4 for b in coefficients}

            def add(index):
                for b, row in state.items():
                    feature = monomial(keys[index], b)
                    row[0] += feature
                    for channel in range(3):
                        row[channel + 1] += feature * values[index][channel]

            if not causal:
                for j in range(len(keys)):
                    add(j)
            outputs = []
            for i, q in enumerate(queries):
                if causal:
                    add(i)
                visible = list(range(i + 1 if causal else len(keys)))
                direct_kernel = [direct(q, keys[j]) for j in visible]
                denominator = sum(direct_kernel)
                require(denominator >= Fraction(len(visible), 2), "normalization has a strict exact lower bound")
                totals = [sum(c * monomial(q, b) * state[b][a] for b, c in coefficients.items())
                          for a in range(4)]
                numerators = [sum(cost * values[j][a] for j, cost in zip(visible, direct_kernel))
                              for a in range(3)]
                require(totals == [denominator] + numerators, "streaming sufficient statistics equal direct sums")
                output = [x / totals[0] for x in totals[1:]]
                require(output == [x / denominator for x in numerators], "normalized attention is exactly preserved")
                output_checks += len(output)
                outputs.append([str(x) for x in output])
                for axis in range(dimension):
                    # Both query gradients and a key masked out of early causal prefixes.
                    for query_side in (True, False):
                        changed_key = len(keys) - 2
                        changes = [direct_partial(q, keys[j], axis, query_side)
                                   if query_side or j == changed_key else Fraction(0) for j in visible]
                        dden = sum(changes)
                        dnum = [sum(change * values[j][a] for j, change in zip(visible, changes))
                                for a in range(3)]
                        dtotals = [Fraction(0)] * 4
                        for b, c in coefficients.items():
                            if query_side:
                                factor = c * partial(q, b, axis)
                                for a in range(4):
                                    dtotals[a] += factor * state[b][a]
                            elif changed_key in visible:
                                factor = c * monomial(q, b) * partial(keys[changed_key], b, axis)
                                dtotals[0] += factor
                                for a in range(3):
                                    dtotals[a + 1] += factor * values[changed_key][a]
                        require(dtotals == [dden] + dnum, "the streaming derivative respects the causal mask")
                        for a in range(3):
                            expected = (dnum[a] * denominator - numerators[a] * dden) / denominator ** 2
                            actual = (dtotals[a + 1] * totals[0] - totals[a + 1] * dtotals[0]) / totals[0] ** 2
                            require(actual == expected, "normalized query/key derivatives agree exactly")
                            output_derivatives += 1
            samples.append({"causal": causal, "exact_outputs": outputs})
        fixtures.append({"latent_dimension": dimension, "box_degree": degree,
                         "embedding_width": len(support), "collapsed_feature_width": len(coefficients),
                         "generic_symmetric_quadratic_width": 1 + len(support) + comb(len(support) + 1, 2),
                         "positive_weights": {",".join(map(str, a)): str(c) for a, c in weights.items()},
                         "queries": [[str(x) for x in q] for q in queries],
                         "keys": [[str(x) for x in k] for k in keys], "samples": samples})
    return {"scope": "Exact reference computation for deliberately constrained polynomial query/key embeddings, an instance of established polynomial-kernel compression. No new ML algorithm, GPU speedup or quality equivalence is asserted.",
            "identity": "For A={0,...,m}^s minus {0}, S(q,k)=sum_(a in A) theta_a*q^a*k^a and p2(S)=1+S+S^2/2 have support B={0,...,2m}^s. Positive theta give a diagonal separable representation with |B| features.",
            "cost": "With s fixed, D=(m+1)^s-1, M=(2m+1)^s, value width v and T tokens: O(T*M*(s+v+1)) arithmetic, O(M*(v+1)) streaming state, plus latent-encoder costs. Direct coefficient setup is O(D^2*s), outside the token loop but charged whenever weights change. This is a forward streaming bound, not constant-memory training.",
            "limits": ["The latent embedding and tied polynomial coefficients restrict the interaction class; arbitrary Q/K or dense learned mixing do not retain the diagonal collapse.",
                       "The kernel is quadratic Taylor attention, not exact exponential softmax. Exact rational checks do not establish floating-point conditioning or trainability.",
                       "The algebraic reduction is existing kernel-on-variety methodology; scientific novelty requires a stronger comparison and an additional theorem or useful measured model advantage."],
            "width_examples": [{"latent_dimension": s, "box_degree": m, "embedding_width": (m + 1) ** s - 1,
                                "generic_symmetric_quadratic_width": 1 + ((m + 1) ** s - 1) + comb((m + 1) ** s, 2),
                                "collapsed_width": (2 * m + 1) ** s} for s, m in ((1, 64), (2, 7), (3, 3))],
            "pair_kernel_checks": pair_checks, "kernel_derivative_checks": derivative_checks,
            "normalized_output_checks": output_checks, "normalized_derivative_checks": output_derivatives,
            "fixtures": fixtures}


def algebraic_power_delta_controls():
    """Established polynomial collection also preserves scalar-gated delta updates."""
    def multiply(a, b):
        result = {}
        for x, c in a.items():
            for y, d in b.items():
                z = tuple(v + w for v, w in zip(x, y))
                result[z] = result.get(z, Fraction(0)) + c * d
        return result

    def compositions(total, length):
        if length == 1:
            yield (total,)
            return
        for first in range(total + 1):
            for tail in compositions(total - first, length - 1):
                yield (first,) + tail

    fixtures = []
    for dimension, degree in ((1, 2), (2, 1)):
        zero = (0,) * dimension
        support = list(product(range(degree + 1), repeat=dimension))
        theta = {a: Fraction(1 + sum(a), 11) for a in support}
        theta[zero] = Fraction(1)
        points = [tuple(Fraction((i + 2) * (j + 1) % 7 - 3, 5)
                        for j in range(dimension)) for i in range(5)]

        def monomial(point, exponent):
            result = Fraction(1)
            for x, power in zip(point, exponent):
                result *= x ** power
            return result

        def base(point):
            return sum(c * monomial(point, a) for a, c in theta.items())

        for order in (2, 4, 6):
            collected = {zero: Fraction(1)}
            for _ in range(order):
                collected = multiply(collected, theta)
            generic, grouped = [], {}
            for counts in compositions(order, len(support)):
                exponent = tuple(sum(n * a[j] for n, a in zip(counts, support))
                                 for j in range(dimension))
                weight = Fraction(factorial(order))
                for a, n in zip(support, counts):
                    weight *= theta[a] ** n / factorial(n)
                generic.append((exponent, weight))
                grouped[exponent] = grouped.get(exponent, Fraction(0)) + weight
            require(grouped == collected, "independent symmetric multinomials collect to the convolution power")
            require(set(collected) == set(product(range(order * degree + 1), repeat=dimension)) and
                    all(c > 0 for c in collected.values()), "positive full-box power coefficients")
            require(len(generic) == comb(len(support) + order - 1, order), "generic symmetric width")

            def raw(point, exponent):
                return monomial(point, exponent) / base(tuple(x * x for x in point)) ** (order // 2)

            pair_checks = output_checks = 0
            for q in points:
                for k in points:
                    expected = base(tuple(x * y for x, y in zip(q, k))) ** order / (
                        base(tuple(x * x for x in q)) * base(tuple(x * x for x in k))) ** (order // 2)
                    require(sum(c * raw(q, a) * raw(k, a) for a, c in collected.items()) == expected,
                            "collected normalized power features equal the direct kernel")
                    pair_checks += 1
            generic_state = [Fraction(0)] * len(generic)
            collected_state = {a: Fraction(0) for a in collected}
            for t, key in enumerate(points):
                alpha = Fraction(4 + t % 2, 5)
                beta = Fraction(1) if t == 4 else Fraction(1 + t % 3, 3)
                value = Fraction(t - 2, 3)
                old_generic = alpha * sum(c * raw(key, a) * state
                                          for (a, c), state in zip(generic, generic_state))
                old_collected = alpha * sum(c * raw(key, a) * collected_state[a]
                                            for a, c in collected.items())
                require(old_generic == old_collected, "delta read before correction")
                generic_state = [alpha * state + beta * raw(key, a) * (value - old_generic)
                                 for (a, _), state in zip(generic, generic_state)]
                collected_state = {a: alpha * state + beta * raw(key, a) * (value - old_collected)
                                   for a, state in collected_state.items()}
                require(all(state == collected_state[a] for (a, _), state in zip(generic, generic_state)),
                        "every duplicate coordinate remains in the common feature subspace")
                for q in points:
                    require(sum(c * raw(q, a) * state for (a, c), state in zip(generic, generic_state)) ==
                            sum(c * raw(q, a) * collected_state[a] for a, c in collected.items()),
                            "scalar-gated delta outputs agree at every prefix")
                    output_checks += 1
                if beta == 1:
                    require(sum(c * raw(key, a) * collected_state[a] for a, c in collected.items()) == value,
                            "unit-feature delta writes the current value exactly at beta one")
            fixtures.append({"latent_dimension": dimension, "box_degree": degree, "power": order,
                             "generic_width": len(generic), "collected_width": len(collected),
                             "pair_kernel_checks": pair_checks, "prefix_query_checks": output_checks,
                             "points": [[str(x) for x in q] for q in points],
                             "theta": {str(a): str(c) for a, c in theta.items()},
                             "final_collected_state": {str(a): str(c) for a, c in collected_state.items()},
                             "beta_one_last_write": str(value)})
    return {"scope": "Exact reference for a constrained normalized even-power feature map and ordinary scalar-gated delta recurrence. Existing kernel-on-variety and DeltaNet methodology; no novelty, hardware or trained-quality claim.",
            "restriction": "Shared positive feature weights are constant across one sequence; query/key embeddings lie on the specified monomial variety. Arbitrary coordinate-wise gates or mixing need not preserve the subspace.",
            "width_examples": [{"power": p, "expanded_embedding_width": 64,
                                "generic_symmetric_width": comb(64 + p - 1, p),
                                "collected_width": (7 * p + 1) ** 2} for p in (2, 4, 6)],
            "pair_kernel_checks": sum(f["pair_kernel_checks"] for f in fixtures),
            "delta_prefix_query_checks": sum(f["prefix_query_checks"] for f in fixtures),
            "fixtures": fixtures}


def algebraic_channel_gate_controls():
    """Exact closure and a failed-projection control for feature-wise forgetting.

    Store S_nu=sqrt(w_nu)*h_nu, avoiding irrational arithmetic. Reads are
    sum w_nu*raw_nu*h_nu; collected coordinates use c_b in place of w_nu.
    """
    def compositions(total, length):
        if length == 1:
            yield (total,)
        else:
            for first in range(total + 1):
                for tail in compositions(total - first, length - 1):
                    yield (first,) + tail

    fixtures = []
    for dimension, degree, order in ((1, 2, 2), (2, 1, 2), (2, 1, 4)):
        support = list(product(range(degree + 1), repeat=dimension))
        theta = [Fraction(1 + sum(a), 3) for a in support]
        generic, collected = [], {}
        for counts in compositions(order, len(support)):
            exponent = tuple(sum(n * a[j] for n, a in zip(counts, support))
                             for j in range(dimension))
            weight = Fraction(factorial(order))
            for n, coefficient in zip(counts, theta):
                weight *= coefficient ** n / factorial(n)
            generic.append((exponent, weight))
            collected[exponent] = collected.get(exponent, Fraction(0)) + weight

        points = [tuple(Fraction((t + 2) * (j + 1) % 7 - 3, 5)
                        for j in range(dimension)) for t in range(5)]

        def raw(point, exponent):
            norm2 = sum(c * prod(x ** (2 * a) for x, a in zip(point, powers))
                        for powers, c in zip(support, theta))
            return prod(x ** a for x, a in zip(point, exponent)) / norm2 ** (order // 2)

        for gate_family in ("independent_collected_channels", "separable_exponent_gates"):
            full = [[Fraction(0)] * 2 for _ in generic]
            small = {a: [Fraction(0)] * 2 for a in collected}
            reads = closure_checks = latest_writes = 0
            for t, key in enumerate(points):
                if gate_family == "independent_collected_channels":
                    gates = {a: Fraction(1 + (t + sum((j + 1) * x for j, x in enumerate(a))) % 5, 6)
                             for a in collected}
                else:
                    axis_gates = [Fraction(2 + (t + j) % 3, 5) for j in range(dimension)]
                    gates = {a: Fraction(4, 5) * prod(g ** x for g, x in zip(axis_gates, a))
                             for a in collected}
                beta = Fraction(1) if t == 4 else Fraction(1 + t % 3, 3)
                value = [Fraction(t - 2, 3), Fraction(t * t % 5, 4)]
                full_decay = [[gates[a] * x for x in row] for (a, _), row in zip(generic, full)]
                small_decay = {a: [gates[a] * x for x in row] for a, row in small.items()}
                prior_full = [sum(c * raw(key, a) * row[j]
                                  for (a, c), row in zip(generic, full_decay)) for j in range(2)]
                prior_small = [sum(collected[a] * raw(key, a) * row[j]
                                   for a, row in small_decay.items()) for j in range(2)]
                require(prior_full == prior_small, "feature-wise forgetting preserves the correction read")
                full = [[row[j] + beta * raw(key, a) * (value[j] - prior_full[j]) for j in range(2)]
                        for (a, _), row in zip(generic, full_decay)]
                small = {a: [row[j] + beta * raw(key, a) * (value[j] - prior_small[j]) for j in range(2)]
                         for a, row in small_decay.items()}
                require(all(row == small[a] for (a, _), row in zip(generic, full)),
                        "every generic state remains in the grouped invariant subspace")
                closure_checks += 1
                for query in points:
                    expected = [sum(c * raw(query, a) * row[j]
                                    for (a, c), row in zip(generic, full)) for j in range(2)]
                    actual = [sum(collected[a] * raw(query, a) * row[j]
                                  for a, row in small.items()) for j in range(2)]
                    require(expected == actual, "channel-gated delta outputs agree at every prefix")
                    reads += 2
                if beta == 1:
                    require([sum(collected[a] * raw(key, a) * row[j]
                                 for a, row in small.items()) for j in range(2)] == value,
                            "feature-wise forgetting retains the exact current-key delta write")
                    latest_writes += 2
            fixtures.append({"latent_dimension": dimension, "box_degree": degree, "power": order,
                             "gate_family": gate_family, "generic_width": len(generic),
                             "collected_width": len(collected), "state_closure_checks": closure_checks,
                             "prefix_query_coordinate_checks": reads, "latest_write_coordinates": latest_writes,
                             "theta": [str(c) for c in theta],
                             "final_collected_state": {str(a): [str(x) for x in row] for a, row in small.items()}})

    # P(x)=1+x+x^2, order two: exponent two has weights 2 and 1.
    # Write value one at x=1 (unit features), then apply two disjoint masks
    # within that duplicate group. Projecting each gate to its weighted mean
    # matches the first read, but incorrectly retains a value after the second.
    entries = []
    for counts in compositions(2, 3):
        exponent = sum(i * n for i, n in enumerate(counts))
        weight = Fraction(2, prod(factorial(n) for n in counts))
        entries.append((counts, exponent, weight))
    grouped = {a: sum(w for _, b, w in entries if b == a) for a in range(5)}
    full = [Fraction(1, 3)] * len(entries)
    small = {a: Fraction(1, 3) for a in grouped}
    masks = [tuple(counts == keep for counts, _, _ in entries) for keep in ((1, 0, 1), (0, 2, 0))]
    counter_reads = []
    for mask in masks:
        full = [x * g for x, g in zip(full, mask)]
        projected_gate = {a: sum(w * g for (_, b, w), g in zip(entries, mask) if b == a) / c
                          for a, c in grouped.items()}
        small = {a: x * projected_gate[a] for a, x in small.items()}
        exact = sum(w * x / 3 for (_, _, w), x in zip(entries, full))
        projected = sum(grouped[a] * x / 3 for a, x in small.items())
        counter_reads.append({"generic_read": str(exact), "projected_read": str(projected)})
    require(counter_reads == [{"generic_read": "2/9", "projected_read": "2/9"},
                              {"generic_read": "0", "projected_read": "2/27"}],
            "arbitrary duplicate-coordinate gates cannot be silently projected at each step")
    return {"scope": "Elementary invariant-subspace closure for an explicitly constrained KDA-style architecture. No new general attention method, trained-quality advantage or GPU speedup is asserted.",
            "criterion": "For U^T U=I, a linear gate A preserves the compressed state exactly iff (I-U U^T) A U=0. For diagonal generic gates and positive grouping weights, this holds iff the gate is constant within each equal-exponent group.",
            "architecture": "Use any independently learned gate per collected exponent, then the ordinary KDA delta correction. This retains arbitrary collected-channel forgetting at M*v state and O(M*v) recurrent arithmetic, plus all feature, encoder and gate costs.",
            "fixtures": fixtures,
            "state_closure_checks": sum(f["state_closure_checks"] for f in fixtures),
            "prefix_query_coordinate_checks": sum(f["prefix_query_coordinate_checks"] for f in fixtures),
            "latest_write_coordinates": sum(f["latest_write_coordinates"] for f in fixtures),
            "incompatible_gate_countercontrol": {"base_coefficients": [1, 1, 1], "power": 2,
                                                 "write_key": 1, "write_value": 1,
                                                 "kept_generic_multiplicities": [[1, 0, 1], [0, 2, 0]],
                                                 "gate_only_steps_beta": 0, "reads": counter_reads}}


def algebraic_erase_write_controls():
    """Independent exact oracle for grouped Gated DeltaNet-2 recurrence.

    Generic states are stored after dividing row nu by sqrt(weight_nu).
    This makes polynomial unit features and all recurrence checks rational.
    An arbitrary erase gate can be averaged inside each duplicate group;
    the forgetting gate must still preserve the grouped state subspace.
    """
    def compositions(total, length):
        if length == 1:
            yield (total,)
        else:
            for first in range(total + 1):
                for tail in compositions(total - first, length - 1):
                    yield (first,) + tail

    fixtures = []
    for dimension, degree, order in ((1, 2, 2), (1, 2, 4), (2, 1, 2), (2, 1, 4)):
        support = list(product(range(degree + 1), repeat=dimension))
        theta = [Fraction(1 + sum(a), 3) for a in support]
        generic, collected = [], {}
        for counts in compositions(order, len(support)):
            exponent = tuple(sum(n * a[j] for n, a in zip(counts, support))
                             for j in range(dimension))
            weight = Fraction(factorial(order))
            for n, coefficient in zip(counts, theta):
                weight *= coefficient ** n / factorial(n)
            generic.append((exponent, weight))
            collected[exponent] = collected.get(exponent, Fraction(0)) + weight
        points = [tuple(Fraction((t + 2) * (j + 1) % 7 - 3, 5)
                        for j in range(dimension)) for t in range(6)]

        def raw(point, exponent):
            norm2 = sum(c * prod(x ** (2 * a) for x, a in zip(point, powers))
                        for powers, c in zip(support, theta))
            return prod(x ** a for x, a in zip(point, exponent)) / norm2 ** (order // 2)

        full = [[Fraction(0)] * 2 for _ in generic]
        small = {a: [Fraction(0)] * 2 for a in collected}
        reads = nonconstant_groups = 0
        for t, key in enumerate(points):
            decay = {a: Fraction(1 + (t + sum(a)) % 5, 6) for a in collected}
            erase = [Fraction((3 * i + t) % 7, 6) for i in range(len(generic))]
            write = [Fraction((t + j + 2) % 5, 4) for j in range(2)]
            value = [Fraction(t - 2, 3), Fraction(t * t % 5, 4)]
            averaged = {a: sum(w * b for (aa, w), b in zip(generic, erase) if aa == a) / c
                        for a, c in collected.items()}
            require(all(0 <= b <= 1 for b in averaged.values()), "erase averages remain admissible gates")
            nonconstant_groups += sum(len({b for (aa, _), b in zip(generic, erase) if aa == a}) > 1
                                      for a in collected)
            full_decay = [[decay[a] * x for x in row] for (a, _), row in zip(generic, full)]
            small_decay = {a: [decay[a] * x for x in row] for a, row in small.items()}
            erased_full = [sum(c * b * raw(key, a) * row[j]
                               for (a, c), b, row in zip(generic, erase, full_decay)) for j in range(2)]
            erased_small = [sum(collected[a] * averaged[a] * raw(key, a) * row[j]
                                for a, row in small_decay.items()) for j in range(2)]
            require(erased_full == erased_small, "untied generic erase gates give exactly the weighted group read")
            full = [[row[j] + raw(key, a) * (write[j] * value[j] - erased_full[j]) for j in range(2)]
                    for (a, _), row in zip(generic, full_decay)]
            small = {a: [row[j] + raw(key, a) * (write[j] * value[j] - erased_small[j]) for j in range(2)]
                     for a, row in small_decay.items()}
            require(all(row == small[a] for (a, _), row in zip(generic, full)),
                    "decoupled erase/write preserves grouped states at every prefix")
            for query in points:
                expected = [sum(c * raw(query, a) * row[j]
                                for (a, c), row in zip(generic, full)) for j in range(2)]
                actual = [sum(collected[a] * raw(query, a) * row[j]
                              for a, row in small.items()) for j in range(2)]
                require(expected == actual, "all grouped erase/write query coordinates agree")
                reads += 2
        require(nonconstant_groups > 0, "oracle exercises erase gates that differ inside a duplicate group")
        fixtures.append({"latent_dimension": dimension, "box_degree": degree, "power": order,
                         "generic_width": len(generic), "collected_width": len(collected),
                         "state_closure_checks": len(points), "prefix_query_coordinate_checks": reads,
                         "nonconstant_erase_group_steps": nonconstant_groups,
                         "final_collected_state": {str(a): [str(x) for x in row] for a, row in small.items()}})

    # Unit key (3/5,4/5), erase (1,0). The oblique erase I-k*(b*k)^T
    # expands this vector, despite every gate belonging to [0,1].
    key, erase = [Fraction(3, 5), Fraction(4, 5)], [Fraction(1), Fraction(0)]
    before = [Fraction(-1, 2), Fraction(1)]
    prior = sum(b * k * x for b, k, x in zip(erase, key, before))
    after = [x - k * prior for k, x in zip(key, before)]
    norm_before, norm_after = sum(x * x for x in before), sum(x * x for x in after)
    require(norm_before == Fraction(5, 4) and norm_after == Fraction(41, 25) > norm_before,
            "channel-wise erase need not be nonexpansive even for a unit key")
    return {"scope": "Exact algebraic compatibility with the published Gated DeltaNet-2 update, not a new recurrence, trained-model result or performance benchmark.",
            "criterion": "For f=U*phi and U^T*U=I, tied group decay preserves im(U); arbitrary generic erase gates reduce to U^T*diag(b)*U, a diagonal weighted group average. Value-axis write gates are unrestricted.",
            "cost_caveat": "M-row recurrence costs O(M*v), but computing generic F-row erase gates and averaging still costs at least O(F). Generate M-row gates directly for an economical architecture; this does not preserve arbitrary nonlinear gate-encoder parameterizations.",
            "fixtures": fixtures,
            "state_closure_checks": sum(f["state_closure_checks"] for f in fixtures),
            "prefix_query_coordinate_checks": sum(f["prefix_query_coordinate_checks"] for f in fixtures),
            "nonconstant_erase_group_steps": sum(f["nonconstant_erase_group_steps"] for f in fixtures),
            "nonexpansiveness_countercontrol": {"unit_key": [str(x) for x in key],
                                                "erase_gate": [str(x) for x in erase],
                                                "input": [str(x) for x in before],
                                                "output": [str(x) for x in after],
                                                "input_norm_squared": str(norm_before),
                                                "output_norm_squared": str(norm_after)}}


def algebraic_attention_mixing_controls():
    """Exact finite controls for established polynomial collection and state rank."""
    def monomial(point, exponent):
        answer = Fraction(1)
        for x, power in zip(point, exponent):
            answer *= x ** power
        return answer

    def partial(point, exponent, axis):
        if not exponent[axis]:
            return Fraction(0)
        reduced = tuple(x - (i == axis) for i, x in enumerate(exponent))
        return exponent[axis] * monomial(point, reduced)

    fixtures = []
    for dimension, degree in ((1, 3), (2, 1), (2, 2), (3, 1)):
        zero = (0,) * dimension
        support = [x for x in product(range(degree + 1), repeat=dimension) if x != zero]
        collected = list(product(range(2 * degree + 1), repeat=dimension))
        index = {x: i for i, x in enumerate(collected)}
        h = [[Fraction(((i + 1) * 7 + (j + 2) * 3) % 11 - 5, 3 * len(support))
              for j in range(len(support))] for i in range(len(support))]
        terms = [(x, y, h[i][j]) for i, x in enumerate(support) for j, y in enumerate(support)]
        c = [[Fraction(0) for _ in collected] for _ in collected]
        c[index[zero]][index[zero]] = Fraction(1)
        for x, y, value in terms:
            c[index[x]][index[y]] += value
        for x, y, value in terms:
            for xx, yy, second in terms:
                row = tuple(v + w for v, w in zip(x, xx))
                column = tuple(v + w for v, w in zip(y, yy))
                c[index[row]][index[column]] += value * second / 2
        points = [tuple(Fraction(((i + 1) * (axis + 2)) % 7 - 3, 5)
                        for axis in range(dimension)) for i in range(5)]
        values = [(Fraction(i - 2), Fraction((i * i) % 5, 3)) for i in range(5)]
        pair_checks = derivative_checks = output_checks = 0

        def score(q, k):
            return sum(value * monomial(q, x) * monomial(k, y) for x, y, value in terms)

        for q in points:
            for k in reversed(points):
                mq = [monomial(q, x) for x in collected]
                mk = [monomial(k, x) for x in collected]
                expected = 1 + score(q, k) + score(q, k) ** 2 / 2
                actual = sum(mq[i] * c[i][j] * mk[j]
                             for i in range(len(collected)) for j in range(len(collected)))
                require(expected == actual >= Fraction(1, 2), "dense mixing preserves the collected kernel")
                pair_checks += 1
                for axis in range(dimension):
                    for query_side in (True, False):
                        ds = sum(value * (partial(q, x, axis) * monomial(k, y) if query_side
                                          else monomial(q, x) * partial(k, y, axis))
                                 for x, y, value in terms)
                        actual_d = sum(c[i][j] * (partial(q, x, axis) * mk[j] if query_side
                                                  else mq[i] * partial(k, y, axis))
                                       for i, x in enumerate(collected) for j, y in enumerate(collected))
                        require(actual_d == (1 + score(q, k)) * ds, "mixed query and key derivatives agree")
                        derivative_checks += 1
        for causal in (False, True):
            state = [[Fraction(0)] * 3 for _ in collected]

            def add(j):
                for k, exponent in enumerate(collected):
                    feature = monomial(points[j], exponent)
                    for channel, value in enumerate((Fraction(1),) + values[j]):
                        state[k][channel] += feature * value

            if not causal:
                for j in range(len(points)):
                    add(j)
            for i, q in enumerate(points):
                if causal:
                    add(i)
                transformed = [sum(monomial(q, x) * c[k][j] for k, x in enumerate(collected))
                               for j in range(len(collected))]
                totals = [sum(feature * row[channel] for feature, row in zip(transformed, state))
                          for channel in range(3)]
                visible = range(i + 1 if causal else len(points))
                weights = [1 + score(q, points[j]) + score(q, points[j]) ** 2 / 2 for j in visible]
                expected = [sum(weight * ((Fraction(1),) + values[j])[channel]
                                for j, weight in zip(visible, weights)) for channel in range(3)]
                require(totals == expected, "the M-row mixed streaming state equals direct summation")
                require([x / totals[0] for x in totals[1:]] == [x / expected[0] for x in expected[1:]],
                        "dense mixing preserves normalized causal and noncausal outputs")
                output_checks += 2
        fixtures.append({"latent_dimension": dimension, "box_degree": degree,
                         "embedding_width": len(support), "collected_width": len(collected),
                         "score_coefficient_matrix": [[str(x) for x in row] for row in h],
                         "nonzero_kernel_coefficients": sum(bool(x) for row in c for x in row),
                         "off_diagonal_kernel_coefficients": sum(bool(c[i][j]) for i in range(len(collected))
                                                                 for j in range(len(collected)) if i != j),
                         "pair_checks": pair_checks, "query_key_derivative_checks": derivative_checks,
                         "normalized_output_checks": output_checks})
    witnesses = []
    points = [Fraction(i + 1, 6) for i in range(5)]
    for n in range(1, 6):
        moment = [[sum(x ** (i + j) for x in points[:n]) for j in range(5)] for i in range(5)]
        minor = determinant([row[:n] for row in moment[:n]])
        require(minor > 0, "separable token updates accumulate state of growing exact rank")
        witnesses.append({"tokens": n, "exact_rank": n, "nonzero_leading_minor": str(minor),
                          "points": [str(x) for x in points[:n]],
                          "state": [[str(x) for x in row] for row in moment]})
    return {"scope": "Elementary polynomial collection with arbitrary bilinear mixing and a classical Vandermonde state-rank witness; no new ML method, benchmark, unrestricted Q/K compression or trained-quality claim.",
            "identity": "S(q,k)=sum H[a,b]*q^a*k^b gives p2(S)=u(q)^T*C*u(k) on the doubled exponent box, even for dense nonsymmetric H. The coefficient matrix C is generally dense.",
            "cost": "M*(v+1) recurrent state, O(M*v) update/read plus O(M^2) query transformation and O(M^2) coefficient storage; naive coefficient preparation O(D^4*s), charged when H changes. Encoders and training/backward cost are additional.",
            "rank_limit": "For distinct real t_j, Z[r,c]=sum_j t_j^(r+c)=V*V^T has rank min(N,L) in an L-by-L moment state. Rank-one token updates therefore do not imply a fixed-rank state under accumulation.",
            "mixed_fixtures": fixtures, "state_rank_witnesses": witnesses,
            "pair_checks": sum(x["pair_checks"] for x in fixtures),
            "query_key_derivative_checks": sum(x["query_key_derivative_checks"] for x in fixtures),
            "normalized_output_checks": sum(x["normalized_output_checks"] for x in fixtures)}


def bounded_grid_shadow_controls(residual_duality):
    """Cover-sensitive bound inputs, residual duality and border budgets.

    Exact finite controls support the note's uniform proof. No finite-value
    comparison of asymptotic big-O constants or new duality theorem is claimed.
    """
    def heads(staircase, d):
        if not staircase:
            return {(0,) * d}
        border = {tuple(x + int(i == j) for j, x in enumerate(a))
                  for a in staircase for i in range(d)} - staircase
        return {a for a in border if all(tuple(x - int(i == j) for j, x in enumerate(a)) in staircase
                                        for i in range(d) if a[i])}

    def ground(sides, length):
        return {a for a in product(*(range(s) for s in sides)) if prod(x + 1 for x in a) <= length}

    box_records = []
    for sides in ((3, 3), (2, 2, 2), (2, 3, 2), (2, 2, 2, 2)):
        d = len(sides)
        box = list(product(*(range(s) for s in sides)))
        index = {a: j for j, a in enumerate(box)}
        predecessors = [sum(1 << index[tuple(x - int(i == j) for j, x in enumerate(a))]
                            for i in range(d) if a[i]) for a in box]
        count, max_inside, samples = 0, 0, []
        for mask in range(1 << len(box)):
            if any(mask >> j & 1 and mask & pred != pred for j, pred in enumerate(predecessors)):
                continue
            staircase = {a for j, a in enumerate(box) if mask >> j & 1}
            residual = {tuple(s - 1 - x for s, x in zip(sides, a)) for a in set(box) - staircase}
            colon_standard = {a for a in box if any(all(x + y < s for x, y, s in zip(a, b, sides))
                                                   for b in set(box) - staircase)}
            require(residual == colon_standard, "independent monomial product oracle gives reflected residual")
            generators = heads(staircase, d)
            inside = {a for a in generators if all(x < s for x, s in zip(a, sides))}
            outside = generators - inside
            pure_powers = {tuple(sides[i] if i == j else 0 for j in range(d)) for i in range(d)}
            require(outside <= pure_powers and len(inside) <= len(box) - len(staircase),
                    "minimal heads split into residual-box monomials and pure powers")
            require(len(generators) <= len(box) - len(staircase) + d, "residual generator budget")
            for state in (staircase, residual):
                require(state <= ground(sides, len(state)), "every box staircase lies in the capped universal ground")
            count += 1
            max_inside = max(max_inside, len(inside))
            if 0 < len(box) - len(staircase) <= 2 and len(samples) < 2:
                samples.append({"length": len(staircase), "residual_length": len(residual),
                                "staircase": sorted(staircase), "minimal_heads": sorted(generators)})
        box_records.append({"sides": sides, "staircases_checked": count,
                            "maximum_inside_heads": max_inside, "near_full_samples": samples})

    quotient_records, quotient_chambers = [], 0
    for fixture in residual_duality["fixtures"]:
        sides = tuple(len(p) - 1 for p in fixture["univariate_moduli"])
        child, other = fixture["child_length"], fixture["residual_length"]
        length = min(child, other)
        chosen = "child_standard_monomials" if child <= other else "residual_standard_monomials"
        domain = ground(sides, length)
        caps = tuple(min(s - 1, length - 1) for s in sides)
        widths = tuple(length * p + 1 for p in caps)
        surface_budget = sum(prod(widths[j] for j in range(len(sides)) if j != i)
                             for i in range(len(sides)))
        max_heads = 0
        for chamber in fixture["all_comparison_chamber_checks"]:
            state = set(map(tuple, chamber[chosen]))
            require(len(state) == length and state <= domain, "actual reduced or nonreduced smaller quotient uses capped ground")
            require(all(sum(a[i] for a in state) <= length * caps[i] for i in range(len(sides))),
                    "actual smaller quotient obeys the integer coordinate budget")
            original = set(map(tuple, chamber["child_standard_monomials"]))
            k = len(heads(original, len(sides)))
            require(k <= prod(sides) - child + len(sides), "actual quotient obeys the residual border bound")
            max_heads = max(max_heads, k)
            quotient_chambers += 1
        quotient_records.append({"label": fixture["label"], "characteristic": fixture["characteristic"],
                                 "degrees": sides, "original_length": child, "residual_length": other,
                                 "smaller_length": length, "capped_ground_size": len(domain),
                                 "coordinate_caps": caps, "regularized_side_lengths": widths,
                                 "half_surface_budget": surface_budget, "maximum_original_minimal_heads": max_heads})

    exponents = []
    for d in range(3, 33):
        D = d * d - d - 1
        alpha, beta = Fraction(d * (d - 2), D), Fraction(d - 1, D)
        old_power = Fraction(d * (d - 1), d + 1)
        threshold = Fraction(3 * d, d + 4)
        records = []
        for gamma in (Fraction(1), (1 + threshold) / 2, threshold, (threshold + d) / 2):
            new = (gamma + 1) * (d - 1) * alpha + Fraction(4, 3) * gamma * beta
            old = (gamma + 1) * old_power
            gap = Fraction(d - 1, 3 * (d + 1) * D) * (3 * d - gamma * (d + 4))
            require(old - new == gap, "exact crossover against the named finite-field upper-bound expression")
            require((gap > 0) == (gamma < threshold), "strict power improvement only in its stated range")
            records.append({"smaller_length_power_of_field_order": str(gamma), "new_power": str(new),
                            "published_expression_power": str(old), "power_gap": str(gap)})
        exponents.append({"dimension": d, "strict_crossover": str(threshold), "samples": records})
    return {"scope": "Uniform cover-sensitive consequence of the lattice-shadow bound and established grid residual duality, over every field including nonreduced separated parents. The original unrestricted objective is retained; no new duality or numerical finite-size dominance is asserted.",
            "parameters": "Q=product b_i, h=min(n,Q-n)>0, p_i=min(b_i-1,h-1), M=#{u:0<=u_i<b_i,product(u_i+1)<=h}, A=sum_i product_(j!=i)(h*p_j+1).",
            "count_bound": "g(I)<=O_d(A^(d(d-2)/(d^2-d-1))*(M*h^(1/3)+1)^((d-1)/(d^2-d-1))).",
            "generator_consequence": "Each reduced basis has at most Q-n+d elements, as well as the classical O_d(n^((d-1)/d)) bound. Multiply the smaller generator budget by the count bound for the union.",
            "box_controls": box_records, "actual_quotient_controls": quotient_records,
            "finite_field_crossover_controls": exponents,
            "box_staircases_checked": sum(f["staircases_checked"] for f in box_records),
            "actual_quotient_chambers_checked": quotient_chambers}


def lattice_shadow_sharpness_controls():
    """Finite certificates for the uniform primitive-zonotope sharpness proof.

    Region counts and integer support objectives are independent finite oracles.
    Neither these small cases nor the density samples prove the asymptotics.
    """
    from functools import cmp_to_key

    def ray_order(x, y):
        hx = 0 if x[1] > 0 or (x[1] == 0 and x[0] > 0) else 1
        hy = 0 if y[1] > 0 or (y[1] == 0 and y[0] > 0) else 1
        if hx != hy:
            return -1 if hx < hy else 1
        turn = x[0] * y[1] - x[1] * y[0]
        return -1 if turn > 0 else 1 if turn < 0 else 0

    zonotopes = []
    for bound in (1, 2, 3):
        generators = [g for g in product(range(-bound, bound + 1), repeat=3)
                      if any(g) and gcd(*g) == 1 and next(x for x in g if x) > 0]
        flats = {primitive(cross(g, h)) for g, h in combinations(generators, 2)}
        degrees = {u: sum(dot(u, g) == 0 for g in generators) for u in flats}
        # Great-circle Euler count: F=2+sum_points(valency/2-1).
        expected = 2 + 2 * sum(k - 1 for k in degrees.values())
        vertices = {}
        total = tuple(sum(g[i] for g in generators) for i in range(3))
        for u in sorted(flats):
            axes = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
            a = next(cross(u, e) for e in axes if any(cross(u, e)))
            b = cross(u, a)
            rays = set()
            for g in generators:
                if dot(u, g) == 0:
                    x, y = -dot(b, g), dot(a, g)
                    factor = gcd(x, y)
                    ray = (x // factor, y // factor)
                    rays.update((ray, tuple(-t for t in ray)))
            rays = sorted(rays, key=cmp_to_key(ray_order))
            for first, second in zip(rays, rays[1:] + rays[:1]):
                x, y = first[0] + second[0], first[1] + second[1]
                perturbation = tuple(x * aa + y * bb for aa, bb in zip(a, b))
                scale = 1 + max(abs(dot(perturbation, g)) for g in generators)
                objective = tuple(scale * uu + vv for uu, vv in zip(u, perturbation))
                signs = [dot(objective, g) for g in generators]
                require(all(signs), "integer objective is strictly inside a zonotope normal cone")
                vertex = tuple(sum(g[i] for g, sign in zip(generators, signs) if sign > 0)
                               for i in range(3))
                vertices[vertex] = objective
                vertices[tuple(t - v for t, v in zip(total, vertex))] = tuple(-x for x in objective)
        require(len(vertices) == expected, "support-vertex enumeration equals the great-circle Euler count")
        samples = []
        for vertex, objective in list(sorted(vertices.items()))[:3]:
            require(dot(objective, vertex) == sum(max(0, dot(objective, g)) for g in generators),
                    "saved support objective exposes its integer vertex")
            samples.append({"vertex": vertex, "strict_objective": objective})
        degree_histogram = {str(k): sum(v == k for v in degrees.values()) for k in sorted(set(degrees.values()))}
        zonotopes.append({"coordinate_bound": bound, "generators": generators,
                          "line_flats": len(flats), "flat_generator_multiplicities": degree_histogram,
                          "vertices_by_euler": expected, "vertices_by_support_objectives": len(vertices),
                          "generic_planar_shadow_vertices": 2 * len(generators),
                          "vertex_set_sha256": sha256(json.dumps(sorted(vertices), separators=(",", ":")).encode()).hexdigest(),
                          "strict_support_samples": samples})

    primitive_boxes = []
    for d in (3, 4):
        h = 100 * d
        columns = [list(product(*[(h, h + 1) if i == j else (0, 1) for i in range(d)]))
                   for j in range(d - 1)]
        counts, normals, samples = 0, {}, []
        for chosen in product(*columns):
            matrix = [[chosen[j][i] for j in range(d - 1)] for i in range(d)]
            cofactors = tuple((-1) ** i * int(determinant([row for j, row in enumerate(matrix) if j != i]))
                              for i in range(d))
            if gcd(*cofactors) != 1:
                continue
            require(all(gcd(*column) == 1 for column in chosen), "saturated columns are individually primitive")
            normal = primitive(cofactors)
            require(all(dot(normal, column) == 0 for column in chosen), "primitive normal annihilates its axis-box tuple")
            normals[normal] = normals.get(normal, 0) + 1
            counts += 1
            if len(samples) < 3:
                samples.append({"columns": chosen, "signed_maximal_minors": cofactors})
        require(counts > 0, "finite primitive axis-box witness exists")
        require(max(normals.values()) <= 5 ** ((d - 1) ** 2), "finite flat multiplicities respect the proved uniform cap")
        primitive_boxes.append({"dimension": d, "axis_scale": h,
                                "tuples_examined": prod(len(box) for box in columns),
                                "primitive_tuples": counts, "distinct_spans": len(normals),
                                "max_span_multiplicity": max(normals.values()), "samples": samples,
                                "scope": "Small density controls, below the sufficient asymptotic threshold; no finite threshold inference."})

    inverse_certificates = []
    exponents = []
    for d in range(3, 33):
        D = d * d - d - 1
        alpha, beta = Fraction(d * (d - 2), D), Fraction(d - 1, D)
        vertex_power, surface_power, shadow_power = d * (d - 1), d * d - 1, d
        require(alpha * surface_power + beta * shadow_power == vertex_power,
                "primitive zonotopes saturate the geometric exponent balance in fixed dimension")
        exponents.append({"dimension": d, "surface_power": surface_power, "shadow_power": shadow_power,
                          "vertex_power": vertex_power, "surface_exponent": str(alpha),
                          "shadow_exponent": str(beta)})
        if d <= 6:
            h = 100 * d
            matrix = [[h if i == j else 1 if i == j + 1 else 0
                       for j in range(d - 1)] for i in range(d)]
            left_inverse = [[0 if j <= i else (-h) ** (j - i - 1) for j in range(d)]
                            for i in range(d - 1)]
            require(all(sum(left_inverse[i][k] * matrix[k][j] for k in range(d)) == int(i == j)
                        for i in range(d - 1) for j in range(d - 1)), "integer left-inverse saturation certificate")
            inverse_certificates.append({"dimension": d, "axis_scale": h,
                                         "matrix": matrix, "integer_left_inverse": left_inverse})
    return {"scope": "Primitive zonotopes are an established family. The note supplies a uniform proof that their growth saturates the proposed lattice surface/shadow inequality; these finite controls do not establish novelty, the asymptotics or quotient realizability.",
            "uniform_growth": "For fixed d>=3 and primitive generators up to sign in [-p,p]^d: f0=Theta_d(p^(d*(d-1))), S=Theta_d(p^(d*d-1)), and mean shadow=Theta_d(p^d).",
            "limitation": "Neither geometric exponent can be reduced with the other held fixed in a uniform pure-power bound. This does not prove sharpness of G_d or of the looser Dey-based matroid specialization.",
            "three_dimensional_zonotopes": zonotopes,
            "axis_box_controls": primitive_boxes, "integer_left_inverse_certificates": inverse_certificates,
            "exponent_balance_controls": exponents,
            "support_vertices_checked": sum(z["vertices_by_support_objectives"] for z in zonotopes),
            "primitive_box_tuples_checked": sum(z["tuples_examined"] for z in primitive_boxes)}


def lattice_shadow_controls():
    """Exact local geometry, shadow sections and exponent-balance certificates."""
    def subtract(x, y):
        return tuple(a - b for a, b in zip(x, y))

    def cross(x, y):
        return (x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2],
                x[0] * y[1] - x[1] * y[0])

    def hull2(points):
        points = sorted(set(points))
        if len(points) < 3:
            return points

        def turn(a, b, c):
            return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

        halves = []
        for order in (points, list(reversed(points))):
            half = []
            for p in order:
                while len(half) >= 2 and turn(half[-2], half[-1], p) <= 0:
                    half.pop()
                half.append(p)
            halves.append(half[:-1])
        return halves[0] + halves[1]

    rng = Random(20261018)

    def polytope(name, supplied):
        points = sorted(set(supplied))
        require(any(determinant([subtract(p, points[0]) for p in triple])
                    for triple in combinations(points[1:], 3)), "the geometry fixture is full dimensional")
        planes = {}
        for a, b, c in combinations(points, 3):
            normal = cross(subtract(b, a), subtract(c, a))
            if not any(normal):
                continue
            factor = gcd(*normal)
            normal = tuple(x // factor for x in normal)
            level = dot(normal, a)
            sides = [dot(normal, p) - level for p in points]
            if max(sides) > 0 and min(sides) < 0:
                continue
            if max(sides) > 0:
                normal, level = tuple(-x for x in normal), -level
            planes[normal + (level,)] = [p for p in points if dot(normal, p) == level]
        facets, edges, vertices = [], set(), set()
        for plane, face_points in sorted(planes.items()):
            normal, level = plane[:3], plane[3]
            drop = max(range(3), key=lambda i: abs(normal[i]))
            axes = [i for i in range(3) if i != drop]
            projected = {tuple(p[i] for i in axes): p for p in face_points}
            polygon = [projected[p] for p in hull2(projected)]
            require(len(polygon) >= 3, "supporting plane is a facet")
            vertices.update(polygon)
            for a, b in zip(polygon, polygon[1:] + polygon[:1]):
                edges.add(tuple(sorted((a, b))))
            multiples = []
            for j in range(1, len(polygon) - 1):
                area_vector = cross(subtract(polygon[j], polygon[0]), subtract(polygon[j + 1], polygon[0]))
                pivot = next(i for i, x in enumerate(normal) if x)
                multiple = Fraction(area_vector[pivot], normal[pivot])
                require(multiple.denominator == 1 and multiple != 0, "integer simplex covolume multiple")
                require(area_vector == tuple(multiple * x for x in normal), "simplex area vector is normal")
                multiples.append(abs(int(multiple)))
            normalized_area = Fraction(sum(multiples), 2)
            require(len(multiples) == len(polygon) - 2, "facet fan triangulation uses every vertex")
            require(len(polygon) <= 6 * normalized_area, "facet incidence is paid by lattice area")
            facets.append({"normal": normal, "level": level, "vertices": polygon,
                           "simplex_normal_multiples": multiples,
                           "area_divided_by_normal_length": str(normalized_area)})
        vertices = sorted(vertices)
        require(len(vertices) - len(edges) + len(facets) == 2, "independent Euler control")
        cone_certificates = []
        for v in vertices:
            incident = [f["normal"] for f in facets if v in f["vertices"]]
            chosen = next((triple for triple in combinations(incident, 3) if determinant(triple)), None)
            require(chosen is not None, "incident primitive normals span at every vertex")
            det = abs(determinant(chosen))
            require(det.denominator == 1 and det >= 1, "the normal determinant supplies an integer lower bound")
            cone_certificates.append({"vertex": v, "incident_normals": incident,
                                      "independent_normals": chosen, "absolute_determinant": str(det),
                                      "normal_squared_lengths": [dot(z, z) for z in chosen]})
        thresholds = []
        for h in (1, 2, 3, 4, 8, 16):
            long_facets = [f for f in facets if dot(f["normal"], f["normal"]) > h * h]
            bad = {v for f in long_facets for v in f["vertices"]}
            incidences = sum(len(f["vertices"]) for f in long_facets)
            require(len(bad) <= incidences, "long-normal vertices are bounded by facet incidences")
            for v in set(vertices) - bad:
                require(all(dot(z, z) <= h * h for z in next(c for c in cone_certificates if c["vertex"] == v)["incident_normals"]),
                        "good vertices have only short facet normals")
            thresholds.append({"H": h, "bad_vertices": len(bad), "long_facet_incidences": incidences,
                               "good_vertices": len(vertices) - len(bad)})
        sections = []
        for _ in range(300):
            a, b = [tuple(rng.randrange(-17, 18) for _ in range(3)) for _ in range(2)]
            if not any(cross(a, b)):
                continue
            image = {(dot(a, p), dot(b, p)): p for p in vertices}
            if len(image) != len(vertices):
                continue
            shadow = hull2(image)
            crossing_edges = []
            nongeneric = False
            for x, y in sorted(edges):
                direction = subtract(y, x)
                ad, bd = dot(a, direction), dot(b, direction)
                if not ad and not bd:
                    nongeneric = True
                    break
                w = tuple(bd * aa - ad * bb for aa, bb in zip(a, b))
                for sign in (1, -1):
                    cost = tuple(sign * z for z in w)
                    maximum = max(dot(cost, p) for p in vertices)
                    exposed = [p for p in vertices if dot(cost, p) == maximum]
                    if x in exposed and y in exposed:
                        if len(exposed) != 2:
                            nongeneric = True
                            break
                        crossing_edges.append({"edge": (x, y), "weight": cost})
                if nongeneric:
                    break
            if nongeneric:
                continue
            require(len(crossing_edges) == len(shadow), "normal-cone crossings equal projected polygon vertices")
            exposed_pairs = {tuple(sorted(entry["edge"])) for entry in crossing_edges}
            require(exposed_pairs == {tuple(sorted((image[x], image[y])))
                                      for x, y in zip(shadow, shadow[1:] + shadow[:1])},
                    "independent primal hull and normal-cone section agree edge by edge")
            sections.append({"projection_rows": (a, b), "shadow_vertices": shadow,
                             "normal_cone_crossings": crossing_edges})
            if len(sections) == 8:
                break
        require(len(sections) == 8, "eight generic exact projection certificates")
        return {"name": name, "vertices": vertices, "facets": facets, "edges": sorted(edges),
                "vertex_normal_certificates": cone_certificates,
                "threshold_controls": thresholds, "shadow_section_controls": sections}

    cube = list(product((0, 1), repeat=3))
    fixtures = [polytope("unit cube", cube),
                polytope("integer simplex", [(0, 0, 0), (3, 0, 0), (0, 4, 0), (0, 0, 5)]),
                polytope("nonsimple octahedron", [tuple(sign * (i == j) for i in range(3))
                                                for j in range(3) for sign in (-1, 1)]),
                polytope("square pyramid", [(0, 0, 0), (3, 0, 0), (0, 3, 0), (3, 3, 0), (1, 1, 4)]),
                polytope("unimodular skew cube", [(x + 5 * y + 12 * z, y + 3 * z, z) for x, y, z in cube])]
    matroids = []
    for label, powers, characteristic in (("generic rational point quotient", (1, 3, 9), 0),
                                           ("finite-field point quotient", (1, 3, 9), 5),
                                           ("diagonal rational point quotient", (1, 1, 1), 0)):
        ground = hyperbolic_cross(3, 3)
        columns = {a: tuple(pow(t, dot(a, powers), characteristic) if characteristic
                            else t ** dot(a, powers) for t in (1, 2, 3)) for a in ground}
        bases = [tuple(ground.index(a) for a in subset) for subset in combinations(ground, 3)
                 if (int(determinant([columns[a] for a in subset])) % characteristic if characteristic
                     else determinant([columns[a] for a in subset]))]
        require(bases, "the point quotient has rank three")
        moments = sorted({tuple(sum(ground[j][i] for j in base) for i in range(3)) for base in bases})
        if powers == (1, 1, 1):
            require(all(sum(p) == 3 for p in moments), "the diagonal basis polytope is lower dimensional")
            triangle = [(3, 0, 0), (0, 3, 0), (0, 0, 3)]
            require(set(triangle) <= set(moments), "axis bases expose the whole triangle")
            regularized = [tuple(x + y for x, y in zip(p, c)) for p in triangle for c in cube]
            geometry = polytope("lower-dimensional basis polytope plus unit cube", regularized)
            require(len(geometry["vertices"]) >= len(triangle), "cube regularization preserves vertex count")
        else:
            geometry = polytope(label + " basis polytope", moments)
        fixtures.append(geometry)
        for _ in range(100):
            intercepts = [rng.randrange(-1000, 1001) for _ in ground]
            slopes = rng.sample(range(-1000, 1001), len(ground))
            crossings = [Fraction(intercepts[j] - intercepts[i], slopes[i] - slopes[j])
                         for i, j in combinations(range(len(ground)), 2)]
            if len(set(crossings)) == len(crossings):
                break
        require(len(set(crossings)) == len(crossings), "generic parametric element costs")
        crossings = sorted(crossings)
        samples = [crossings[0] - 1] + [(x + y) / 2 for x, y in zip(crossings, crossings[1:])] + [crossings[-1] + 1]
        sequence, changes, owners = [], [], {}
        for interval, t in enumerate(samples):
            costs = [x + t * y for x, y in zip(intercepts, slopes)]
            selected = set()
            for j in sorted(range(len(ground)), key=lambda j: costs[j]):
                if any(selected | {j} <= set(base) for base in bases):
                    selected.add(j)
            minima = min(sum(costs[j] for j in base) for base in bases)
            require(sum(costs[j] for j in selected) == minima and len(selected) == 3,
                    "greedy and independent exhaustive minimum-basis enumeration agree")
            require(sum(sum(costs[j] for j in base) == minima for base in bases) == 1,
                    "generic cost samples have a unique minimum basis")
            if not sequence:
                owners = {j: token for token, j in enumerate(sorted(selected))}
            elif selected != set(sequence[-1]["basis"]):
                previous = set(sequence[-1]["basis"])
                removed, added = previous - selected, selected - previous
                require(len(removed) == len(added) == 1, "every generic transition is a single basis exchange")
                x, y = next(iter(removed)), next(iter(added))
                crossing = crossings[interval - 1]
                require(intercepts[x] + crossing * slopes[x] == intercepts[y] + crossing * slopes[y],
                        "the exchanged element lines meet at the transition")
                require(slopes[y] < slopes[x], "each token follows a convex chain with decreasing slope")
                token = owners.pop(x)
                owners[y] = token
                changes.append({"at": str(crossing), "removed": x, "added": y, "token": token,
                                "old_slope": slopes[x], "new_slope": slopes[y]})
            require(set(owners) == selected and len(set(owners.values())) == 3,
                    "token chains occupy different element lines on every open interval")
            sequence.append({"sample": str(t), "basis": sorted(selected), "minimum_cost": str(minima)})
        matroids.append({"name": label, "characteristic": characteristic, "rank": 3, "ground": ground,
                         "columns": [[str(x) for x in columns[a]] for a in ground], "bases": bases,
                         "basis_moments": moments, "intercepts": intercepts, "slopes": slopes,
                         "parametric_samples": sequence, "token_exchanges": changes})
    exponents = []
    for d in range(3, 33):
        denominator = d * d - d - 1
        a = Fraction(d * (d - 2), d - 1)
        surface_power, shadow_power = a / (a + 1), 1 / (a + 1)
        exponent = 2 * (d - 1) * surface_power + Fraction(4, 3) * shadow_power
        log_power = (d - 1) * shadow_power
        old = Fraction(2 * d * (d - 1), d + 1)
        gain = Fraction(2 * (d - 1) * (d - 2), 3 * (d + 1) * denominator)
        require(exponent == Fraction((d - 1) * (6 * d * (d - 2) + 4), 3 * denominator), "exact exponent balance")
        require(old - exponent == gain > 0, "strict fixed-power gain over the historical exponent")
        require(log_power == Fraction((d - 1) ** 2, denominator), "exact logarithmic exponent")
        exponents.append({"dimension": d, "normal_cutoff_exponent": str((2 * (d - 1) - Fraction(4, 3)) / (a + 1)),
                          "bound_exponent": str(exponent), "log_exponent": str(log_power),
                          "historical_exponent": str(old), "power_gain": str(gain),
                          "universal_basis_exponent": str(exponent + Fraction(d - 1, d))})
    return {"scope": "Exact finite supporting certificates for the uniform proof in the note; they do not numerically integrate random shadows, prove isoperimetry, establish novelty or replace the dimension-uniform argument.",
            "polytope_fixtures": fixtures, "point_matroid_fixtures": matroids, "exponent_controls": exponents,
            "facet_certificates": sum(len(f["facets"]) for f in fixtures),
            "vertex_normal_certificates": sum(len(f["vertex_normal_certificates"]) for f in fixtures),
            "shadow_sections": sum(len(f["shadow_section_controls"]) for f in fixtures),
            "parametric_samples": sum(len(f["parametric_samples"]) for f in matroids),
            "token_exchanges": sum(len(f["token_exchanges"]) for f in matroids)}


def staircase_generator_controls():
    """Independent border/fiber checks for the known generator-count lemma."""
    fixtures = []

    def minimal_border(state, dimension):
        border = {tuple(x + (axis == j) for j, x in enumerate(a))
                  for a in state for axis in range(dimension)} - state
        return {u for u in border if all(tuple(x - (axis == j) for j, x in enumerate(u)) in state
                                         for axis in range(dimension) if u[axis])}

    for dimension, maximum_size in ((2, 7), (3, 7), (4, 5)):
        states = {frozenset({(0,) * dimension})}
        for n in range(1, maximum_size + 1):
            following = set()
            for state in sorted(states, key=lambda a: sorted(a)):
                generators = minimal_border(state, dimension)
                # Independent oracle: take the first missing point on each
                # coordinate fiber, then remove dominated caps.
                caps = set()
                for axis in range(dimension):
                    fibers = {}
                    for a in state:
                        base = a[:axis] + a[axis + 1:]
                        fibers[base] = max(fibers.get(base, -1), a[axis])
                    for base, highest in fibers.items():
                        caps.add(base[:axis] + (highest + 1,) + base[axis:])
                independent = {u for u in caps if not any(v != u and all(x <= y for x, y in zip(v, u))
                                                         for v in caps)}
                require(generators == independent, "minimal border and independent fiber-cap oracle agree")
                require(all(not all(x <= y for x, y in zip(u, v))
                            for u, v in permutations(generators, 2)), "minimal outside exponents form an antichain")
                h = next(x for x in range(1, n + 1) if x ** dimension >= n)
                small = {u for u in generators if max(u) < h}
                require(len({u[:-1] for u in small}) == len(small) <= h ** (dimension - 1),
                        "small antichain has an injective coordinate projection")
                chain_counts = []
                for axis in range(dimension):
                    large = [u for u in generators if u[axis] >= h]
                    chains = [tuple(tuple(x - t * (axis == j) for j, x in enumerate(u))
                                    for t in range(1, h + 1)) for u in large]
                    union = {p for chain in chains for p in chain}
                    require(len(union) == h * len(large) and union <= state,
                            "large minimal generators give disjoint inside chains")
                    chain_counts.append(len(large))
                bound = h ** (dimension - 1) + dimension * Fraction(n, h)
                require(len(generators) <= len(small) + sum(chain_counts) <= bound,
                        "exact chain-packing generator bound")
                fixtures.append({"dimension": dimension, "n": n, "staircase": sorted(state),
                                 "minimal_outside": sorted(generators), "h": h,
                                 "small_generators": len(small), "large_chain_counts": chain_counts,
                                 "rational_generator_upper_bound": str(bound)})
                if n < maximum_size:
                    following.update(state | {u} for u in generators)
            states = following
    return {"scope": "Finite supporting certificates for the established minimal-generator estimate used in the universal-basis corollary. No new generator theorem or uniform proof by enumeration is claimed.",
            "method": "Enumerate small staircases; compare one-step minimal borders against independent coordinate-fiber caps, and certify the disjoint-chain packing bound exactly.",
            "fixtures": fixtures, "staircases_checked": len(fixtures),
            "minimal_generators_checked": sum(len(f["minimal_outside"]) for f in fixtures)}


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
    minimum_products = minimum_product_rank_controls()
    product_excess = product_excess_and_extension_controls()
    quadratic_poles = quadratic_pole_flag_controls()
    moment_cubes = staircase_moment_cube_controls()
    algebraic_attention = algebraic_attention_hotpath_controls()
    attention_mixing = algebraic_attention_mixing_controls()
    power_delta = algebraic_power_delta_controls()
    channel_gates = algebraic_channel_gate_controls()
    erase_write = algebraic_erase_write_controls()
    lattice_shadows = lattice_shadow_controls()
    shadow_sharpness = lattice_shadow_sharpness_controls()
    bounded_grids = bounded_grid_shadow_controls(residual_duality)
    staircase_generators = staircase_generator_controls()
    exchange_cycles = exchange_cycle_coherence_controls()
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
        "minimum_product_rank_control": minimum_products,
        "product_excess_and_extension_control": product_excess,
        "quadratic_pole_flag_control": quadratic_poles,
        "staircase_moment_cube_control": moment_cubes,
        "algebraic_attention_hotpath_control": algebraic_attention,
        "algebraic_attention_mixing_control": attention_mixing,
        "algebraic_power_delta_control": power_delta,
        "algebraic_channel_gate_control": channel_gates,
        "algebraic_erase_write_control": erase_write,
        "lattice_shadow_control": lattice_shadows,
        "lattice_shadow_sharpness_control": shadow_sharpness,
        "bounded_grid_shadow_control": bounded_grids,
        "staircase_generator_control": staircase_generators,
        "exchange_cycle_coherence_control": exchange_cycles,
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
                    "small_perturbation_staircase_checks": sum(len(f["nonzero_scaled_basis_determinants"]) for f in grading_controls["perturbation_fixtures"]),
                    "minimum_product_rank_fixtures": len(minimum_products["fixtures"]),
                    "minimum_product_staircase_checks": sum(len(f["all_staircase_certificates"]) for f in minimum_products["fixtures"]),
                    "minimum_product_fan_witness_checks": sum(f["initial_ideal_count"] for f in minimum_products["fixtures"]),
                    "minimum_product_box_rank_checks": sum(len(f["box_rank_checks"]) for f in minimum_products["fixtures"]),
                    "minimum_product_stabilizer_checks": sum(len(f["prefix_checks"]) for f in minimum_products["polynomial_prefix_stabilizer_controls"]),
                    "minimum_product_inverse_rectangle_checks": sum(len(f["box_rank_checks"]) for f in minimum_products["unsaturated_inverse_negative_controls"]),
                    "product_excess_reconstruction_fixtures": len(product_excess["reconstruction_fixtures"]),
                    "product_excess_threshold_counters": len(product_excess["degree_threshold_counters"]),
                    "bounded_extension_staircase_checks": sum(len(f["all_staircase_ranks"]) for f in product_excess["bounded_extension_fixtures"]),
                    "bounded_extension_quotient_witnesses": sum(len(f["quotient_inclusion_witnesses"]) for f in product_excess["bounded_extension_fixtures"]),
                    "high_degree_boundary_controls": len(product_excess["high_degree_boundary_controls"]),
                    "quadratic_pole_fixtures": len(quadratic_poles["fixtures"]),
                    "quadratic_pole_pattern_rank_checks": sum(len(s["all_pattern_checks"]) for f in quadratic_poles["fixtures"] for s in f["sectors"]),
                    "quadratic_pole_independent_flag_checks": sum(f["independent_flag_rank_checks"] for f in quadratic_poles["fixtures"]),
                    "quadratic_pole_complete_fan_witnesses": sum(sum(r["status"] == "coherent" for r in s["all_pattern_checks"]) for f in quadratic_poles["fixtures"] for s in f["sectors"]),
                    "quadratic_pole_selected_weight_checks": sum(len(f["selected_weight_checks"]) for f in quadratic_poles["fixtures"]),
                    "quadratic_pole_hyperbolic_checks": sum(f["hyperbolic_witness_checks"] for f in quadratic_poles["fixtures"]),
                    "quadratic_pole_prefix_delay_counters": len(quadratic_poles["prefix_delay_countercontrols"]),
                    "staircase_moment_row_checks": len(moment_cubes["row_certificates"]),
                    "staircase_moment_expanded_union_checks": len(moment_cubes["expanded_union_checks"]),
                    "staircase_moment_positive_relaxation_witnesses": sum(len(f["positive_minimum_certificates"]) for f in moment_cubes["integer_ball_relaxation_checks"]),
                    "staircase_moment_prefix_dominance_counters": len(moment_cubes["prefix_dominance_countercontrols"]),
                    "algebraic_attention_pair_checks": algebraic_attention["pair_kernel_checks"],
                    "algebraic_attention_kernel_derivative_checks": algebraic_attention["kernel_derivative_checks"],
                    "algebraic_attention_output_checks": algebraic_attention["normalized_output_checks"],
                    "algebraic_attention_normalized_derivative_checks": algebraic_attention["normalized_derivative_checks"],
                    "mixed_attention_pair_checks": attention_mixing["pair_checks"],
                    "mixed_attention_query_key_derivative_checks": attention_mixing["query_key_derivative_checks"],
                    "mixed_attention_normalized_output_checks": attention_mixing["normalized_output_checks"],
                    "attention_accumulated_state_rank_witnesses": len(attention_mixing["state_rank_witnesses"]),
                    "algebraic_power_pair_checks": power_delta["pair_kernel_checks"],
                    "algebraic_power_delta_prefix_query_checks": power_delta["delta_prefix_query_checks"],
                    "algebraic_channel_gate_state_closure_checks": channel_gates["state_closure_checks"],
                    "algebraic_channel_gate_prefix_query_coordinates": channel_gates["prefix_query_coordinate_checks"],
                    "algebraic_erase_write_state_closure_checks": erase_write["state_closure_checks"],
                    "algebraic_erase_write_query_coordinates": erase_write["prefix_query_coordinate_checks"],
                    "lattice_shadow_facet_certificates": lattice_shadows["facet_certificates"],
                    "lattice_shadow_vertex_normal_certificates": lattice_shadows["vertex_normal_certificates"],
                    "lattice_shadow_exact_sections": lattice_shadows["shadow_sections"],
                    "lattice_shadow_parametric_samples": lattice_shadows["parametric_samples"],
                    "lattice_shadow_token_exchanges": lattice_shadows["token_exchanges"],
                    "lattice_shadow_exponent_controls": len(lattice_shadows["exponent_controls"]),
                    "lattice_shadow_sharpness_support_vertices": shadow_sharpness["support_vertices_checked"],
                    "lattice_shadow_sharpness_box_tuples": shadow_sharpness["primitive_box_tuples_checked"],
                    "bounded_grid_shadow_staircases": bounded_grids["box_staircases_checked"],
                    "bounded_grid_shadow_quotient_chambers": bounded_grids["actual_quotient_chambers_checked"],
                    "staircase_generator_staircases_checked": staircase_generators["staircases_checked"],
                    "staircase_generator_minimal_generators_checked": staircase_generators["minimal_generators_checked"],
                    "exchange_cycle_rank_two_contractions": sum(x["invertible_new_pairs"] for x in exchange_cycles["rank_two_contractions"]),
                    "exchange_cycle_patch_checks": sum(x["weight_checks"] for x in exchange_cycles["patch_controls"]),
                    "exchange_cycle_rational_controls": len(exchange_cycles["rational_cycle_controls"]),
                    "exchange_cycle_expanded_templates": len(exchange_cycles["expanded_template_controls"])},
        "goal_qualifying_advance": False,
    }
    target = ROOT / "results/research/structural_complexity_checks.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
