#!/usr/bin/env python3
"""Exact interval checks for the containment poset of pyramidal partitions.

No floating point, LP package, or precomputed Mobius values are used. Every
downset of each finite top is considered, recognized by strict separation,
and the Mobius recurrence is compared with the proposed single-layer formula.
This is a verification artifact, not a scalable enumeration algorithm.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import heapq
from itertools import combinations, product
import json
from math import comb, gcd, lcm
from pathlib import Path
import shlex
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def primitive_row(row):
    denominator = lcm(*(Fraction(x).denominator for x in row))
    integers = tuple(int(Fraction(x) * denominator) for x in row)
    divisor = gcd(*integers)
    return tuple(x // divisor for x in integers) if divisor else integers


def nullspace(equations, dimension):
    """Rational basis columns of the common kernel, and its codimension."""
    matrix = [[Fraction(x) for x in row] for row in equations]
    pivots = []
    for column in range(dimension):
        pivot = next((i for i in range(len(pivots), len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        at = len(pivots)
        matrix[at], matrix[pivot] = matrix[pivot], matrix[at]
        value = matrix[at][column]
        matrix[at] = [x / value for x in matrix[at]]
        for i in range(len(matrix)):
            if i != at and matrix[i][column]:
                factor = matrix[i][column]
                matrix[i] = [x - factor * y for x, y in zip(matrix[i], matrix[at])]
        pivots.append(column)
    basis = []
    for free in range(dimension):
        if free not in pivots:
            vector = [Fraction(int(j == free)) for j in range(dimension)]
            for i, pivot in enumerate(pivots):
                vector[pivot] = -matrix[i][free]
            basis.append(vector)
    return basis, len(pivots)


def strict_cone_witness(rows, dimension):
    """Exact Fourier--Motzkin for homogeneous strict inequalities row.x > 0.

    Returns a rational witness or None with a contradictory projected row.
    Choosing the least costly elimination coordinate limits small-fixture growth.
    There is no numerical feasibility tolerance.
    """
    rows = sorted(set(primitive_row(row) for row in rows))
    if any(not any(row) for row in rows):
        return None
    if dimension == 0:
        return ()
    column = min(range(dimension), key=lambda j:
                 sum(row[j] > 0 for row in rows) * sum(row[j] < 0 for row in rows))
    positive = [row for row in rows if row[column] > 0]
    negative = [row for row in rows if row[column] < 0]
    projected = [row[:column] + row[column + 1:] for row in rows if row[column] == 0]
    for p in positive:
        for n in negative:
            projected.append(tuple(-n[column] * p[j] + p[column] * n[j]
                                   for j in range(dimension) if j != column))
    smaller = strict_cone_witness(projected, dimension - 1)
    if smaller is None:
        return None

    def boundary(row):
        rest = row[:column] + row[column + 1:]
        return -sum(x * y for x, y in zip(rest, smaller)) / Fraction(row[column])

    lower = max((boundary(row) for row in positive), default=None)
    upper = min((boundary(row) for row in negative), default=None)
    if lower is not None and upper is not None:
        require(lower < upper, "strict projected interval")
        value = (lower + upper) / 2
    elif lower is not None:
        value = lower + 1
    elif upper is not None:
        value = upper - 1
    else:
        value = Fraction(0)
    result = smaller[:column] + (value,) + smaller[column:]
    require(all(sum(a * x for a, x in zip(row, result)) > 0 for row in rows),
            "substitution in original strict system")
    return result


def separated_weights(inequalities, equations, dimension):
    basis, rank = nullspace(equations, dimension)
    reduced = [tuple(sum(a * b for a, b in zip(row, vector)) for vector in basis)
               for row in inequalities]
    coordinates = strict_cone_witness(reduced, len(basis))
    if coordinates is None:
        return None, rank
    weights = tuple(sum(c * vector[i] for c, vector in zip(coordinates, basis))
                    for i in range(dimension))
    weights = primitive_row(weights)
    require(all(sum(a * x for a, x in zip(row, weights)) > 0 for row in inequalities),
            "lifted strict witness")
    require(all(sum(a * x for a, x in zip(row, weights)) == 0 for row in equations),
            "lifted equality witness")
    return weights, rank


def vertex_cone_witness(rows, dimension):
    """Independent small-dimensional oracle: enumerate vertices of row.x>=1.

    Used only with the coordinate positivity rows included, which make the
    feasible polyhedron pointed and guarantee a vertex whenever it is nonempty.
    Laplace determinants and substitution are independent of Fourier--Motzkin.
    """
    def determinant(matrix):
        if not matrix:
            return 1
        return sum((-1) ** j * x * determinant(
            [row[:j] + row[j + 1:] for row in matrix[1:]])
            for j, x in enumerate(matrix[0]))

    rows = sorted(set(primitive_row(row) for row in rows))
    for chosen in combinations(rows, dimension):
        matrix = [list(row) for row in chosen]
        denominator = determinant(matrix)
        if not denominator:
            continue
        sign = 1 if denominator > 0 else -1
        denominator = abs(denominator)
        numerators = tuple(sign * determinant(
            [[1 if j == column else x for j, x in enumerate(row)] for row in matrix])
            for column in range(dimension))
        if all(sum(a * x for a, x in zip(row, numerators)) >= denominator for row in rows):
            return numerators
    return None


def border(points, dimension):
    """Minimal excluded points, including the empty downset case."""
    points = frozenset(points)
    if not points:
        return ((0,) * dimension,)
    candidates = set()
    for point in points:
        for i in range(dimension):
            child = tuple(x + (i == j) for j, x in enumerate(point))
            if child not in points:
                candidates.add(child)
    return tuple(sorted(point for point in candidates if all(
        tuple(x - (i == j) for j, x in enumerate(point)) in points
        for i, x in enumerate(point) if x)))


def maximal(points, dimension):
    points = frozenset(points)
    return tuple(point for point in points if all(
        tuple(x + (i == j) for j, x in enumerate(point)) not in points
        for i in range(dimension)))


def cut_witness(points, dimension):
    if not points:
        return {"weights": [1] * dimension, "threshold_twice": -1}
    inequalities = [tuple(int(i == j) for j in range(dimension)) for i in range(dimension)]
    inside = maximal(points, dimension)
    outside = border(points, dimension)
    inequalities.extend(tuple(y - x for x, y in zip(a, b)) for a in inside for b in outside)
    weights, _ = separated_weights(inequalities, [], dimension)
    if weights is None:
        return None
    low = max(sum(a * w for a, w in zip(point, weights)) for point in inside)
    high = min(sum(a * w for a, w in zip(point, weights)) for point in outside)
    require(low < high, "cut excludes its entire minimal border")
    return {"weights": list(weights), "threshold_twice": low + high}


def single_layer(lower, upper, dimension):
    difference = sorted(upper - lower)
    if not difference:
        return {"mobius": 1, "affine_dimension": -1, "weights": None, "threshold": None}
    anchor = difference[0]
    equations = [tuple(y - x for x, y in zip(anchor, point)) for point in difference[1:]]
    inequalities = [tuple(int(i == j) for j in range(dimension)) for i in range(dimension)]
    inequalities.extend(tuple(x - y for x, y in zip(anchor, point))
                        for point in maximal(lower, dimension))
    inequalities.extend(tuple(y - x for x, y in zip(anchor, point))
                        for point in border(upper, dimension))
    weights, rank = separated_weights(inequalities, equations, dimension)
    return {"mobius": (-1) ** (rank + 1) if weights is not None else 0,
            "affine_dimension": rank, "weights": list(weights) if weights is not None else None,
            "threshold": sum(x * w for x, w in zip(anchor, weights)) if weights is not None else None}


def all_downsets(top, dimension):
    points = sorted(top, key=lambda point: (sum(point), point))
    positions = {point: i for i, point in enumerate(points)}
    require(len(points) == len(positions), "distinct top points")
    prerequisites = []
    for i, point in enumerate(points):
        needed = 0
        for j in range(dimension):
            if point[j]:
                parent = tuple(x - (j == c) for c, x in enumerate(point))
                require(parent in positions and positions[parent] < i, "top is a downset")
                needed |= 1 << positions[parent]
        prerequisites.append(needed)
    masks = [0]
    for i, needed in enumerate(prerequisites):
        masks += [mask | (1 << i) for mask in masks if mask & needed == needed]
    return points, sorted(masks, key=lambda mask: (mask.bit_count(), mask))


def check_fixture(name, top, dimension, crosscheck_vertices=False):
    require(cut_witness(top, dimension) is not None, ("fixture top is a corner cut", name))
    points, downsets = all_downsets(top, dimension)
    cuts = []
    cut_points = {}
    witnesses = {}
    vertex_checks = 0
    for mask in downsets:
        selected = frozenset(point for i, point in enumerate(points) if mask & (1 << i))
        witness = cut_witness(selected, dimension)
        if crosscheck_vertices and selected:
            rows = [tuple(int(i == j) for j in range(dimension)) for i in range(dimension)]
            rows.extend(tuple(y - x for x, y in zip(a, b))
                        for a in maximal(selected, dimension) for b in border(selected, dimension))
            independent = vertex_cone_witness(rows, dimension)
            require((witness is None) == (independent is None),
                    ("Fourier--Motzkin vs independent vertex feasibility", name, mask))
            vertex_checks += 1
        if witness is not None:
            cuts.append(mask)
            cut_points[mask] = selected
            witnesses[mask] = witness
    counts = Counter()
    examples = {}
    intervals = []
    for lower in cuts:
        nonzero = {}
        for upper in cuts:
            if upper & lower != lower:
                continue
            # The recurrence uses all previously processed nonzero terms.
            value = 1 if lower == upper else -sum(mu for middle, mu in nonzero.items()
                                                  if middle & upper == middle)
            if value:
                nonzero[upper] = value
            predicted = single_layer(cut_points[lower], cut_points[upper], dimension)
            require(value == predicted["mobius"],
                    ("Mobius recurrence vs single-layer criterion", name, lower, upper, value, predicted))
            counts[str(value)] += 1
            intervals.append([lower, upper, value])
            label = ("nonzero_rank_" + str(predicted["affine_dimension"] + 1)
                     if value and lower != upper else "zero" if value == 0 else "diagonal")
            if label not in examples:
                examples[label] = {"lower_mask": lower, "upper_mask": upper,
                                   "lower_cut_witness": witnesses[lower],
                                   "upper_cut_witness": witnesses[upper], **predicted}
    return {"name": name, "dimension": dimension, "top_points": points,
            "downsets_examined": len(downsets), "corner_cuts": len(cuts),
            "independent_vertex_feasibility_checks": vertex_checks,
            "cut_masks": cuts, "comparable_intervals_checked": sum(counts.values()),
            "mobius_value_counts": dict(counts), "examples": examples,
            "intervals": intervals}


def simplex_points(dimension, degree):
    return frozenset(point for point in product(range(degree + 1), repeat=dimension)
                     if sum(point) <= degree)


def removable_face_terms(upper, dimension):
    """Find every contributing face from affine bases, without listing subcuts.

    Every positive exposed layer is contained in the coordinatewise maximal
    points and has affine dimension at most k-1. An affine basis of at most k
    such points determines its full intersection with upper. The strict layer
    test filters these candidates. This is a small exact implementation, not
    a polynomial-time LP implementation in variable dimension.
    """
    extremal = sorted(maximal(upper, dimension))
    seen = set()
    terms = []
    for size in range(1, min(dimension, len(extremal)) + 1):
        for basis in combinations(extremal, size):
            anchor = basis[0]
            normals, rank = nullspace(
                [tuple(y - x for x, y in zip(anchor, point)) for point in basis[1:]],
                dimension)
            if rank != size - 1:
                continue
            layer = frozenset(point for point in upper if all(
                sum(n * (x - a) for n, x, a in zip(normal, point, anchor)) == 0
                for normal in normals))
            if layer in seen:
                continue
            seen.add(layer)
            witness = single_layer(upper - layer, upper, dimension)
            if witness["mobius"]:
                require(witness["affine_dimension"] == rank, "affine face basis rank")
                terms.append({"layer": sorted(layer), "basis": basis,
                              "coefficient": -witness["mobius"], **witness})
    return terms


def check_face_recurrence(fixture):
    """Compare geometric terms with every nonzero Mobius column and direct counts."""
    points = fixture["top_points"]
    positions = {tuple(point): i for i, point in enumerate(points)}
    cuts = fixture["cut_masks"]
    dimension = fixture["dimension"]
    counts = {upper: sum(lower & upper == lower for lower in cuts) for upper in cuts}
    columns = {upper: {} for upper in cuts}
    for lower, upper, value in fixture["intervals"]:
        if value and lower != upper:
            columns[upper][lower] = value
    cases = []
    for upper in cuts:
        selected = frozenset(tuple(point) for i, point in enumerate(points) if upper & (1 << i))
        terms = removable_face_terms(selected, dimension)
        recovered = {}
        for term in terms:
            layer_mask = sum(1 << positions[point] for point in term["layer"])
            lower = upper ^ layer_mask
            require(lower in counts, ("removed face leaves a corner cut", fixture["name"], upper))
            require(lower not in recovered, "distinct contributing faces")
            recovered[lower] = term["mobius"]
            term["lower_mask"] = lower
            term["lower_subpartition_count"] = counts[lower]
        require(recovered == columns[upper],
                ("affine-basis enumeration vs complete Mobius column", fixture["name"], upper))
        recurrence = 1 + sum(term["coefficient"] * term["lower_subpartition_count"] for term in terms)
        require(recurrence == counts[upper],
                ("removable-face recurrence vs direct subcut count", fixture["name"], upper))
        require(len(terms) <= sum(comb(len(selected), j) for j in range(1, min(dimension, len(selected)) + 1)),
                "general face-count bound")
        if dimension == 3:
            require(len(terms) <= 6 * len(selected), "linear three-dimensional face-count bound")
        cases.append({"upper_mask": upper, "size": len(selected),
                      "direct_subpartition_count": counts[upper], "recurrence_value": recurrence,
                      "contributing_faces": terms})
    return {"fixture": fixture["name"], "cases_checked": len(cases),
            "nonzero_proper_mobius_entries": sum(len(case["contributing_faces"]) for case in cases),
            "max_contributing_faces": max(len(case["contributing_faces"]) for case in cases),
            "cases": cases}


def single_layer_controls(fixtures):
    cases = [("partial_line_is_not_a_complete_layer", 2, simplex_points(2, 1),
              simplex_points(2, 1) | {(2, 0), (1, 1)}, 0, "plane_degree_4"),
             ("complete_line_of_three_points", 2, simplex_points(2, 1),
              simplex_points(2, 2), 1, "plane_degree_4"),
             ("complete_plane_of_six_points", 3, simplex_points(3, 1),
              simplex_points(3, 2), -1, "space_degree_2")]
    results = []
    for name, dimension, lower, upper, expected, fixture_name in cases:
        predicted = single_layer(lower, upper, dimension)
        fixture = next(row for row in fixtures if row["name"] == fixture_name)
        positions = {point: i for i, point in enumerate(fixture["top_points"])}
        a = sum(1 << positions[point] for point in lower)
        b = sum(1 << positions[point] for point in upper)
        recurrence = next(mu for x, y, mu in fixture["intervals"] if x == a and y == b)
        require(recurrence == predicted["mobius"] == expected, ("targeted layer control", name))
        results.append({"name": name, "dimension": dimension, "lower": sorted(lower),
                        "upper": sorted(upper), "recurrence_mobius": recurrence, **predicted})
    return results


def quadrilateral_layer_control():
    weights, threshold = (2, 3, 5), 11
    box = product(*(range(threshold // w + 1) for w in weights))
    upper = frozenset(point for point in box if sum(x * w for x, w in zip(point, weights)) <= threshold)
    lower = frozenset(point for point in upper if sum(x * w for x, w in zip(point, weights)) < threshold)
    layer = sorted(upper - lower)
    require(len(layer) == 4, "quadrilateral has four lattice points")
    masks = []
    cut_witnesses = {}
    for mask in sorted(range(16), key=lambda m: (m.bit_count(), m)):
        selected = lower | {point for i, point in enumerate(layer) if mask & (1 << i)}
        witness = cut_witness(selected, 3)
        if witness is not None:
            masks.append(mask)
            cut_witnesses[mask] = witness
    values = {}
    for mask in masks:
        values[mask] = 1 if mask == 0 else -sum(mu for old, mu in values.items() if old & mask == old)
        selected = lower | {point for i, point in enumerate(layer) if mask & (1 << i)}
        require(values[mask] == single_layer(lower, selected, 3)["mobius"],
                ("quadrilateral recurrence vs criterion", mask))
    require(len(masks) == 14 and values[15] == -1, "four-atom nonzero interval")
    joins = []
    for size in (2, 3):
        for atoms in combinations((1, 2, 4, 8), size):
            union = sum(atoms)
            bounds = [mask for mask in masks if mask & union == union]
            minimal = [mask for mask in bounds if not any(
                other != mask and other & mask == other for other in bounds)]
            require(minimal != [15], "top is not the join of two or three atoms")
            joins.append({"atom_masks": atoms, "minimal_upper_bound_masks": minimal})
    primes = (3, 5, 13)
    def prime_value(point):
        value = 1
        for prime, exponent in zip(primes, point):
            value *= prime ** exponent
        return value
    inside_values = [{"exponents": point, "value": prime_value(point)}
                     for point in sorted(maximal(lower, 3))]
    outside_values = [{"exponents": point, "value": prime_value(point)}
                      for point in border(upper, 3)]
    layer_values = sorted(prime_value(point) for point in layer)
    require(max(row["value"] for row in inside_values) == 243 < min(layer_values),
            "exact product certificate for the lower prefix")
    require(max(layer_values) < min(row["value"] for row in outside_values) == 507,
            "exact product certificate excludes the entire upper complement")
    heap, seen, smooth = [(1, (0, 0, 0))], {(0, 0, 0)}, []
    for _ in range(len(upper) + 1):
        value, point = heapq.heappop(heap)
        smooth.append((value, point))
        for i, prime in enumerate(primes):
            child = tuple(x + (i == j) for j, x in enumerate(point))
            if child not in seen:
                seen.add(child)
                heapq.heappush(heap, (value * prime, child))
    require(frozenset(point for _, point in smooth[:len(lower)]) == lower,
            "quadrilateral lower cut is an exact prime prefix")
    require(frozenset(point for _, point in smooth[:len(upper)]) == upper,
            "quadrilateral upper cut is an exact prime prefix")
    return {"dimension": 3, "cutting_weights": weights, "cutting_threshold": threshold,
            "lower_points": sorted(lower), "upper_points": sorted(upper), "layer": layer,
            "interval_cut_masks": masks, "cut_witnesses": cut_witnesses,
            "mobius_from_lower_by_mask": values, "top_mobius": values[15],
            "two_and_three_atom_upper_bounds": joins,
            "prime_border_certificate": {"maximal_lower": inside_values,
                                         "minimal_outside_upper": outside_values},
            "prime_prefix_realization": {"primes": primes, "lower_rank": len(lower),
                "upper_rank": len(upper), "boundary_terms": [
                    {"rank": rank, "value": value, "exponents": point}
                    for rank, (value, point) in enumerate(smooth, 1) if rank >= len(lower)]}}


def slice_obstruction_points(m):
    """Downset with z-heights [[2m+2,m],[m+1,1],[1,0]] (rows indexed by y)."""
    require(isinstance(m, int) and m >= 1, "positive slice-family parameter")
    heights = ((2 * m + 2, m), (m + 1, 1), (1, 0))
    return frozenset((x, y, z) for y, row in enumerate(heights)
                     for x, height in enumerate(row) for z in range(height))


def slice_certificate(m, axis, level):
    """Positive integer weights and twice the strict threshold for every slice."""
    if axis == 0:
        if level == 0:
            return (2 * m + 1, 2), 8 * m + 5
        if level == 1:
            return (2 * m - 1, 2), 4 * m - 1
    elif axis == 1:
        if level == 0:
            return (2 * m + 3, 2), 8 * m + 5
        if level == 1:
            return (m, 1), 2 * m + 1
        if level == 2:
            return (1, 1), 1
    elif axis == 2:
        if level == 0:
            return (3, 2), 11
        if 1 <= level < m:
            return (1, 1), 3
        if level == m:
            return (2, 1), 3
        if m < level <= 2 * m + 1:
            return (1, 1), 1
    return (1, 1), -1


def check_slice_obstruction_family(max_m=32):
    witnesses = []
    total_slices = 0
    for m in range(1, max_m + 1):
        points = slice_obstruction_points(m)
        require(len(points) == 4 * m + 5, ("slice-family size", m))
        require(all(tuple(v - (i == j) for j, v in enumerate(point)) in points
                    for point in points for i, x in enumerate(point) if x),
                ("slice-family downset", m))
        inside = ((0, 0, 2 * m + 1), (1, 1, 0))
        outside = ((1, 0, m), (0, 1, m + 1))
        require(all(point in points for point in inside) and all(point not in points for point in outside),
                ("balanced pair membership", m))
        require(tuple(a + b for a, b in zip(*inside)) == tuple(a + b for a, b in zip(*outside)),
                ("exact equal-sum obstruction", m))
        require(cut_witness(points, 3) is None, ("independent global LP infeasibility", m))
        slices = []
        for axis in range(3):
            last = max(point[axis] for point in points)
            for level in range(last + 2):
                section = frozenset(tuple(x for i, x in enumerate(point) if i != axis)
                                    for point in points if point[axis] == level)
                weights, threshold_twice = slice_certificate(m, axis, level)
                require(all(w > 0 for w in weights), ("positive slice weights", m, axis, level))
                require(all(2 * sum(x * w for x, w in zip(point, weights)) < threshold_twice
                            for point in section), ("slice included points", m, axis, level))
                require(all(2 * sum(x * w for x, w in zip(point, weights)) > threshold_twice
                            for point in border(section, 2)),
                        ("slice excludes entire complement", m, axis, level))
                independent = cut_witness(section, 2)
                require(independent is not None, ("independent slice LP feasibility", m, axis, level))
                if section:
                    total_slices += 1
                slices.append({"axis": axis, "level": level, "points": sorted(section),
                               "weights": weights, "threshold_twice": threshold_twice,
                               "independent_cut_witness": independent})
        witnesses.append({"parameter": m, "points": sorted(points), "size": len(points),
                          "inside_balanced_pair": inside, "outside_balanced_pair": outside,
                          "global_cut_lp_feasible": False, "slices": slices})
    return {"status": "uniform family proved; prior-art priority unresolved",
            "source_question": "Galvan thesis (2023), Conjecture 4.1, printed p. 51",
            "parameters_checked": list(range(1, max_m + 1)),
            "nonempty_slices_checked": total_slices,
            "scope": "all coordinate-plane slices of the listed three-dimensional downsets; "
                     "the proof covers all m>=1 and extends to all dimensions k>=3",
            "no_minimality_claim": True, "witnesses": witnesses}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results/research/pyramidal_poset.json")
    args = parser.parse_args()
    # Controls distinguish infeasibility, equations, and genuinely strict boundaries.
    require(strict_cone_witness([(1,), (-1,)], 1) is None, "opposite strict inequalities")
    require(strict_cone_witness([(0,)], 1) is None, "zero cannot be strictly positive")
    require(separated_weights([(1, 0), (0, 1)], [(1, 1)], 2)[0] is None,
            "positive weights cannot sum to zero")
    require(separated_weights([(1, 0), (0, 1)], [(1, -1)], 2)[0] is not None,
            "positive equal weights exist")
    fixtures = [("line_degree_6", 1, 6), ("plane_degree_4", 2, 4),
                ("space_degree_2", 3, 2), ("space_degree_3", 3, 3),
                ("four_dimensions_degree_2", 4, 2)]
    results = []
    for name, dimension, degree in fixtures:
        result = check_fixture(name, simplex_points(dimension, degree), dimension,
                               crosscheck_vertices=dimension <= 2 or dimension == 3 and degree == 2)
        results.append(result)
        print(json.dumps({key: result[key] for key in (
            "name", "downsets_examined", "corner_cuts", "comparable_intervals_checked",
            "mobius_value_counts")}), flush=True)
    result = {
        "task": "exact checks of the pyramidal-partition Mobius formula and coordinate-slice obstruction family",
        "status": "proofs via an interval representation and Edelman (1984), and explicit slice certificates; novelty audit pending",
        "scope": "all downsets and all comparable corner-cut intervals inside the listed tops; "
                 "not all partitions of a fixed mass and not a prime-distribution result",
        "command": shlex.join(["python3", str(Path(__file__).resolve().relative_to(ROOT)), *sys.argv[1:]]),
        "python_version": sys.version,
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in (Path(__file__).resolve(), ROOT / "docs/interleaving_research.md",
                                       ROOT / "paper/pyramidal_mobius.tex")},
        "all_checks_passed": True, "fixtures": results,
        "targeted_single_layer_controls": single_layer_controls(results),
        "quadrilateral_layer_control": quadrilateral_layer_control(),
        "slice_obstruction_family": check_slice_obstruction_family(),
        "removable_face_recurrence": [check_face_recurrence(fixture) for fixture in results],
        "total_intervals_checked": sum(row["comparable_intervals_checked"] for row in results),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"all_checks_passed": True,
                      "total_intervals_checked": result["total_intervals_checked"],
                      "output": str(args.output.resolve())}), flush=True)


if __name__ == "__main__":
    main()
