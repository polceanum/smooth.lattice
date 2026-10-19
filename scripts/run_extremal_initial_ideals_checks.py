#!/usr/bin/env python3
"""Bounded exact controls for the renewed extremal initial-ideal goal.

These are primary-bound controls and a discriminating complete-intersection
probe, not a proof of an asymptotic improvement. Core kernels are unchanged.
"""
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, gcd, lcm, prod
from pathlib import Path
import hashlib
import json

from run_interleaving_research import minimal_excluded_border
from run_structural_complexity_checks import (
    area2, clip, cross, dot, greedy, hyperbolic_cross, positive_chambers, primitive,
)

ROOT = Path(__file__).resolve().parents[1]


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def inverse_mod(matrix, prime):
    n = len(matrix)
    rows = [[x % prime for x in row] + [int(i == j) for j in range(n)]
            for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i][j]), None)
        if pivot is None:
            return None
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scale = pow(rows[j][j], -1, prime)
        rows[j] = [x * scale % prime for x in rows[j]]
        for i in range(n):
            if i != j:
                scale = rows[i][j]
                rows[i] = [(x - scale * y) % prime for x, y in zip(rows[i], rows[j])]
    return [row[n:] for row in rows]


def grid_controls():
    packing = []
    for p, d in ((2, 2), (2, 3), (2, 4), (3, 2), (3, 3)):
        ground = tuple(product(range(p), repeat=d))
        states = {frozenset()}
        maxima = [[-1] * d for _ in range(len(ground) + 1)]
        count = 0
        while states:
            following = set()
            for state in states:
                count += 1
                for i in range(d):
                    maxima[len(state)][i] = max(maxima[len(state)][i], sum(x[i] for x in state))
                for x in ground:
                    if x not in state and all(tuple(a - (i == j) for i, a in enumerate(x)) in state
                                              for j in range(d) if x[j]):
                        following.add(state | {x})
            states = following
        for n, values in enumerate(maxima):
            q, r = divmod(n, p)
            expected = q * p * (p - 1) // 2 + r * (r - 1) // 2
            require(values == [expected] * d, "published coordinate packing versus every grid downset")
        packing.append({"prime": p, "dimension": d, "downsets_checked": count,
                        "coordinate_maxima_by_cardinality": maxima})

    # Exhaust every point subset and every positive monomial-order chamber
    # of the Boolean three-cube. This checks the actual reflected complement,
    # rather than the later chapter's literal equality of standard sets.
    p, d = 2, 3
    ground = tuple(product(range(p), repeat=d))
    normals = sorted({primitive(tuple(a - b for a, b in zip(x, y)))
                      for x, y in combinations(ground, 2)
                      if min(a - b for a, b in zip(x, y)) < 0 < max(a - b for a, b in zip(x, y))})
    weights = positive_chambers(normals)
    families = {}
    for mask in range(1 << len(ground)):
        points = [x for i, x in enumerate(ground) if mask >> i & 1]
        columns = {a: tuple(product_value(a, x, p) for x in points) for a in ground}
        families[mask] = [frozenset(greedy(ground, columns, w, p)) for w in weights]
        require(all(len(basis) == len(points) for basis in families[mask]), "point quotient dimension")
    full = (1 << len(ground)) - 1
    comparisons = 0
    maximum_counts = [0] * (len(ground) + 1)
    for mask, bases in families.items():
        complement = families[full ^ mask]
        for a, b in zip(bases, complement):
            reflected = {tuple(p - 1 - x for x in u) for u in a}
            require(b == frozenset(ground).difference(reflected), "correct grid-complement bijection")
            comparisons += 1
        require(len(set(bases)) == len(set(complement)), "equal fans imply equal numbers of initial ideals")
        maximum_counts[mask.bit_count()] = max(maximum_counts[mask.bit_count()], len(set(bases)))
    require(maximum_counts == [1, 1, 3, 3, 3, 3, 3, 1, 1], "later chapter Table 2")
    singleton = 1
    require(len(families[singleton][0]) == 1 and len(families[full ^ singleton][0]) == 7,
            "standard monomial sets cannot literally be equal for complementary cardinalities")
    return {"packing": packing, "boolean_point_subsets": len(families),
            "positive_order_chambers": len(weights), "reflected_complement_checks": comparisons,
            "maximum_initial_ideal_counts_by_cardinality": maximum_counts,
            "scope": "Independent finite checks of established grid packing and reflected-complement/fan invariance. The literal standard-set equality in the later chapter is not used."}


def product_value(exponent, point, prime):
    value = 1
    for a, x in zip(exponent, point):
        value = value * pow(x, a, prime) % prime
    return value


def rectangle_model(r, s, coefficients, prime):
    """A=F_p[x,y,z]/(x^r,y^s,z-F(x,y)), length r*s."""
    cells = tuple(product(range(r), range(s)))
    index = {a: i for i, a in enumerate(cells)}
    n = len(cells)
    require((0, 0) not in coefficients, "the probe uses a local nilpotent model")
    powers = [[int(i == 0) for i in range(n)]]
    for _ in range(n):
        vector = [0] * n
        for a, value in zip(cells, powers[-1]):
            for b, coefficient in coefficients.items():
                target = (a[0] + b[0], a[1] + b[1])
                if target in index:
                    vector[index[target]] = (vector[index[target]] + value * coefficient) % prime
        powers.append(vector)

    def residue(exponent):
        a, b, c = exponent
        if c > n or a >= r or b >= s:
            return (0,) * n
        vector = [0] * n
        for (i, j), value in zip(cells, powers[c]):
            if i + a < r and j + b < s:
                vector[index[(i + a, j + b)]] = value
        return tuple(vector)

    ground = tuple(hyperbolic_cross(n, 3))
    columns = {a: residue(a) for a in ground}
    require([residue((a, b, 0)) for a, b in cells] ==
            [tuple(int(i == j) for i in range(n)) for j in range(n)], "exact rectangular reference basis")
    matrices = []
    for direction in range(3):
        matrix = [residue((a + (direction == 0), b + (direction == 1), int(direction == 2)))
                  for a, b in cells]
        matrices.append(matrix)
    apply = lambda matrix, vector: tuple(sum(x[i] * c for x, c in zip(matrix, vector)) % prime
                                        for i in range(n))
    for i, j in combinations(range(3), 2):
        for column in range(n):
            require(apply(matrices[i], matrices[j][column]) ==
                    apply(matrices[j], matrices[i][column]), "multiplication commutes on every reference vector")
    for exponent in ground:
        vector = tuple(int(i == 0) for i in range(n))
        for direction, power in enumerate(exponent):
            for _ in range(power):
                vector = apply(matrices[direction], vector)
        require(vector == columns[exponent], "matrix powers versus direct truncated polynomial products")
    return ground, columns, residue


def integer_weight(point):
    scale = lcm(*(x.denominator for x in point))
    result = tuple(int(x * scale) for x in point)
    factor = gcd(*result)
    return tuple(x // factor for x in result)


def cone_for_basis(basis, residue, prime):
    basis = tuple(sorted(basis))
    n = len(basis)
    inverse = inverse_mod([[residue(a)[i] for a in basis] for i in range(n)], prime)
    require(inverse is not None, "every visited standard set is a quotient basis")
    normals = set()
    for u in minimal_excluded_border(basis):
        column = residue(u)
        coefficients = [sum(x * y for x, y in zip(row, column)) % prime for row in inverse]
        require(tuple(sum(c * residue(a)[i] for c, a in zip(coefficients, basis)) % prime
                      for i in range(n)) == column, "independent border-coordinate replay")
        for v, coefficient in zip(basis, coefficients):
            if coefficient:
                z = tuple(a - b for a, b in zip(u, v))
                g = gcd(*z)
                normals.add(tuple(x // g for x in z))
    polygon = [tuple(Fraction(i == j) for i in range(3)) for j in range(3)]
    for normal in sorted(normals):
        polygon = clip(polygon, normal, 1)
        require(polygon, "the known strict weight realizes a full-dimensional cone")
    return polygon, sorted(normals)


def adjacent_weights(polygon, normals, comparison_normals):
    result = []
    for p, q in zip(polygon, polygon[1:] + polygon[:1]):
        z = next((a for a in normals if dot(a, p) == dot(a, q) == 0), None)
        if z is None or any(p[i] == q[i] == 0 for i in range(3)):
            continue
        # A relative edge point must avoid unrelated comparison walls.
        k = 1
        while True:
            midpoint = tuple((k * a + b) / (k + 1) for a, b in zip(p, q))
            if all(dot(a, midpoint) or cross(a, z) == (0, 0, 0) for a in comparison_normals):
                break
            k += 1
            require(k <= len(comparison_normals) + 2, "finite generic-edge avoidance")
        direction = tuple(Fraction(a) - Fraction(sum(z), 3) for a in z)
        limits = [Fraction(1)]
        limits += [a / b for a, b in zip(midpoint, direction) if b > 0]
        for a in comparison_normals:
            value, speed = dot(a, midpoint), dot(a, direction)
            if value and speed and value / speed > 0:
                limits.append(value / speed)
        epsilon = min(limits) / 2
        neighbor = tuple(a - epsilon * b for a, b in zip(midpoint, direction))
        require(min(neighbor) > 0 and dot(z, neighbor) < 0 and
                all(dot(a, neighbor) != 0 for a in comparison_normals), "strict adjacent comparison chamber")
        result.append(integer_weight(neighbor))
    return result


def rectangle_fan(r, s, coefficients, prime, exhaustive_control=False):
    ground, columns, residue = rectangle_model(r, s, coefficients, prime)
    active = tuple(a for a in ground if any(columns[a]))
    comparison_normals = sorted({primitive(tuple(x - y for x, y in zip(a, b)))
                                 for a, b in combinations(active, 2)
                                 if min(x - y for x, y in zip(a, b)) < 0 < max(x - y for x, y in zip(a, b))})
    n = r * s
    pending = [(1, n + 1, (n + 1)**2)]
    visited = {}
    while pending:
        weight = pending.pop()
        basis = frozenset(greedy(active, columns, weight, prime))
        require(len(basis) == n, "full quotient rank")
        if basis in visited:
            continue
        polygon, normals = cone_for_basis(basis, residue, prime)
        require(all(dot(z, weight) > 0 for z in normals), "strict reduced-border leading terms")
        require(all(tuple(x - (i == j) for i, x in enumerate(a)) in basis
                    for a in basis for j in range(3) if a[j]), "every visited basis is a downset")
        visited[basis] = (polygon, weight)
        pending.extend(adjacent_weights(polygon, normals, comparison_normals))
    require(sum(area2(poly) for poly, _ in visited.values()) == 1,
            "complete positive weight simplex covered exactly")

    independent = 0
    if exhaustive_control:
        weights = positive_chambers(comparison_normals)
        bases = set()
        # For n=4 use the permutation formula, independent of row reduction.
        all_bases = []
        for candidate in combinations(active, n):
            determinant = 0
            for perm in permutations(range(n)):
                term = (-1)**sum(perm[i] > perm[j] for i in range(n) for j in range(i+1, n))
                for i in range(n):
                    term *= columns[candidate[perm[i]]][i]
                determinant += term
            independent += 1
            if determinant % prime:
                all_bases.append(candidate)
        for w in weights:
            values = [(sum(dot(w, a) for a in b), frozenset(b)) for b in all_bases]
            best = min(x for x, _ in values)
            selected = [b for value, b in values if value == best]
            require(len(selected) == 1 and selected[0] == frozenset(greedy(active, columns, w, prime)),
                    "independent exhaustive determinant/minimum-weight basis oracle")
            bases.add(selected[0])
        require(bases == set(visited), "facet traversal versus complete comparison arrangement")

    witnesses = [{"basis": sorted(basis), "strict_weight": weight,
                  "cone_polygon": [[str(x) for x in vertex] for vertex in poly]}
                 for basis, (poly, weight) in sorted(visited.items(), key=lambda item: sorted(item[0]))]
    return {"rectangle": [r, s], "quotient_length": n, "characteristic": prime,
            "polynomial_terms": [[i, j, c] for (i, j), c in sorted(coefficients.items())],
            "initial_ideal_count": len(visited), "active_ground_monomials": len(active),
            "all_multiplication_commutators_checked": True, "reference_basis_cyclicity_checked": True,
            "positive_fan_exhaustive": True, "normalized_double_area": "1",
            "independent_determinants_checked": independent, "witnesses": witnesses}


def support_matching_obstruction_controls():
    """Replay an infinite-family identity, independently of fan traversal."""
    fixtures = []
    for q in range(3, 9):
        cells = tuple(product(range(q), repeat=2))
        f = {(i, j): (i + 1)**2 + 3*(j + 1) + 5*i*j for i, j in cells if i or j}

        def multiply(left, right):
            result = {}
            for (i, j), a in left.items():
                for (k, l), b in right.items():
                    if i + k < q and j + l < q:
                        target = (i + k, j + l)
                        result[target] = result.get(target, 0) + a*b
            return {p: a for p, a in result.items() if a}

        powers = [{(0, 0): 1}]
        for _ in range(2*q - 2):
            powers.append(multiply(powers[-1], f))
        u, v = powers[2*q - 3], powers[2*q - 2]
        w = multiply({(q-2, q-2): 1}, f)
        a, b = f[(1, 0)], f[(0, 1)]
        c = a*a*f[(0, 2)] + a*b*f[(1, 1)] + b*b*f[(2, 0)]
        coefficients = [2*a*a*b*b, -2*comb(2*q-3, q-1)*(a*b)**q, -(q-2)*c]
        columns = [u, w, v]
        require(all(sum(coefficient * column.get(p, 0) for coefficient, column
                        in zip(coefficients, columns)) == 0 for p in cells),
                "universal square-truncation relation replayed over the integers")
        assignments = [(q-1, q-2), (q-2, q-1), (q-1, q-1)]
        require(all(column.get(p, 0) % 101 for p, column in zip(assignments, columns)),
                "three distinct support rows give a perfect matching")
        exponents = [(0, 0, 2*q-3), (q-2, q-2, 1), (0, 0, 2*q-2)]
        residues = {e: tuple(column.get(p, 0) % 101 for p in cells)
                    for e, column in zip(exponents, columns)}
        rank = len(greedy(exponents, residues, (1, 17, 289), 101))
        require(rank == 2, "support matching does not determine quotient rank")
        padding = max(0, 2*(q-1)**2 - q*q)
        length = q*q + padding
        # Add the reduced points (j,0,0), 1<=j<=padding. All three
        # circuit columns vanish there. Pure x powers q,...,q+padding-1
        # have zero local part and a nonsingular Vandermonde evaluation
        # block, proving a cyclic quotient of the stated larger length.
        extra_block = [[pow(j, q+k, 101) for k in range(padding)]
                       for j in range(1, padding+1)]
        require(not padding or inverse_mod(extra_block, 101) is not None,
                "padding-point block gives the exact CRT quotient length")
        require(all(prod(a+1 for a in e) <= length for e in exponents),
                "all three counterexample columns belong to the universal ground set")
        fixtures.append({"square_side": q, "local_quotient_length": q*q,
                         "padding_points_on_x_axis": padding, "quotient_length": length,
                         "all_columns_in_universal_ground": True,
                         "characteristic_for_rank_control": 101,
                         "monomial_exponents": exponents, "integer_relation_coefficients": coefficients,
                         "sparse_integer_columns": [[[i, j, value] for (i, j), value in sorted(column.items())]
                                                    for column in columns],
                         "support_matching_rows": assignments, "support_matching_rank": 3,
                         "actual_rank": rank})
    return {"fixtures": fixtures, "uniform_identity_proved_in_note": True,
            "scope": "A coefficient-sensitive dependence persists for every square size in the written proof. Exact integer convolution and modular elimination validate the listed fixtures. This rejects a support-only rank shortcut; it is not an asymptotic initial-ideal theorem."}


def main():
    grid = grid_controls()
    rectangles, homogeneous = [], []
    for r, s in ((2, 2), (2, 3), (3, 3), (3, 4)):
        coefficients = {(i, j): (i + 1)**2 + 3*(j + 1) + 5*i*j
                        for i, j in product(range(r), range(s)) if i or j}
        rectangles.append(rectangle_fan(r, s, coefficients, 101, r == s == 2))
        control = rectangle_fan(r, s, {(1, 0): 1, (0, 1): 1}, 101)
        require(control["initial_ideal_count"] == 5, "homogeneous linear control")
        homogeneous.append(control)
    small_characteristic = rectangle_fan(2, 2, {(1, 0): 1, (0, 1): 1}, 2)
    require(small_characteristic["initial_ideal_count"] == 3,
            "small characteristic can collapse homogeneous cones")
    report = {"schema_version": 1, "status": "passed",
              "command": "python3 -B scripts/run_extremal_initial_ideals_checks.py",
              "grid_source_controls": grid, "rectangular_complete_intersections": rectangles,
              "homogeneous_controls": homogeneous,
              "small_characteristic_control": small_characteristic,
              "support_matching_obstruction": support_matching_obstruction_controls(),
              "asymptotic_improvement_established": False, "novelty_established": False,
              "goal_complete": False,
              "scope": "Exhaustive finite primary-bound and candidate-family diagnostics. No polynomial separation, sharp unrestricted asymptotic, or algorithmic speedup is inferred.",
              "source_sha256": {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
                                for p in ("scripts/run_extremal_initial_ideals_checks.py",
                                          "scripts/run_structural_complexity_checks.py",
                                          "scripts/run_interleaving_research.py")}}
    path = ROOT / "results/research/extremal_initial_ideals_checks.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": "passed", "grid_downsets": sum(x["downsets_checked"] for x in grid["packing"]),
                      "grid_complement_cases": grid["reflected_complement_checks"],
                      "complete_intersection_counts": [[x["quotient_length"], x["initial_ideal_count"]] for x in rectangles],
                      "goal_complete": False, "output": str(path.relative_to(ROOT))}))


if __name__ == "__main__":
    main()
