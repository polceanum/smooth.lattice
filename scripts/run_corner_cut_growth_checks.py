#!/usr/bin/env python3
"""Exact falsification checks for the multiscale corner-cut growth argument.

These finite checks do not prove the asymptotic density lemma or establish
novelty. They check primitive-simplex injectivity, rational parameter bounds,
and the vertex/region inequality on a small degenerate arrangement. They
also check the prescribed-size injection using exact affine-lattice patches.
"""
from __future__ import annotations

import hashlib
import heapq
from fractions import Fraction
from itertools import combinations, product
import json
from math import comb, factorial, gcd, isqrt, prod
from pathlib import Path
from random import Random
import shlex
import sys

from run_pyramidal_poset_checks import strict_cone_witness


ROOT = Path(__file__).resolve().parents[1]


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def plane_key(p, q, r):
    normal = cross(subtract(q, p), subtract(r, p))
    height = dot(normal, p)
    if not height:
        return None
    divisor = gcd(*normal, height)
    if height < 0:
        divisor = -divisor
    return tuple(x // divisor for x in normal) + (height // divisor,)


def in_parameter_box(plane, sides):
    return all(9 * plane[3] < 10 * length * coefficient < 11 * plane[3]
               for length, coefficient in zip(sides, plane[:3]))


def axis_boxes(sides):
    result = []
    for axis in range(3):
        ranges = [range(length if j == axis else 0,
                        (length if j == axis else 0) + length // 100 + 1)
                  for j, length in enumerate(sides)]
        result.append(tuple(product(*ranges)))
    return result


def primitive_fixture(sides):
    boxes = axis_boxes(sides)
    seen = {}
    nonprimitive_seen = {}
    collisions = 0
    collision_witness = None
    primitive_count = 0
    coprime_x_count = 0
    coprime_x_primitive_count = 0
    for points in product(*boxes):
        p, q, r = points
        v, u = subtract(q, p), subtract(r, p)
        normal = cross(v, u)
        index = gcd(*normal)
        require(index > 0 and all(x > 0 for x in normal), "positive normal")
        plane = plane_key(*points)
        require(in_parameter_box(plane, sides), "normal is inside parameter box")
        coprime_x = gcd(v[0], u[0]) == 1
        coprime_x_count += coprime_x
        if index == 1:
            primitive_count += 1
            coprime_x_primitive_count += coprime_x
            require(plane not in seen, ("primitive-plane collision", seen.get(plane), points))
            seen[plane] = points
        else:
            if plane in nonprimitive_seen:
                collisions += 1
                if collision_witness is None:
                    collision_witness = {"plane": plane,
                                         "first": nonprimitive_seen[plane], "second": points}
            else:
                nonprimitive_seen[plane] = points
    # Check the affine-basis consequence against every selected lattice point
    # for a deterministic sample of at most 100 primitive planes.
    sampled = list(seen.items())[:100]
    for plane, points in sampled:
        for box, expected in zip(boxes, points):
            incident = [point for point in box if dot(plane[:3], point) == plane[3]]
            require(incident == [expected], "one lattice point in each axis box")
    return {
        "sides": sides,
        "points_per_axis_box": len(boxes[0]),
        "ordered_triples_checked": prod(map(len, boxes)),
        "primitive_triangles": primitive_count,
        "distinct_primitive_planes": len(seen),
        "coprime_x_pairs": coprime_x_count,
        "coprime_x_and_primitive": coprime_x_primitive_count,
        "primitive_injectivity": True,
        "plane_incidence_samples_checked": len(sampled),
        "nonprimitive_plane_collisions": collisions,
        "nonprimitive_collision_witness": collision_witness,
        "primitive_witness": ({"plane": sampled[0][0], "points": sampled[0][1]}
                              if sampled else None),
    }


def arrangement_fixture():
    sides = (100, 100, 100)
    points = ((100, 0, 0), (0, 100, 0), (0, 0, 100), (100, 1, 0),
              (0, 100, 1), (1, 0, 100), (101, 0, 1), (1, 101, 0))
    vertices = set()
    for chosen in combinations(points, 3):
        plane = plane_key(*chosen)
        if plane is not None and in_parameter_box(plane, sides):
            vertices.add(plane)
    bounds = [(0, 0, 0, 1)]
    for axis, length in enumerate(sides):
        lower = [0, 0, 0, -9]
        upper = [0, 0, 0, 11]
        lower[axis], upper[axis] = 10 * length, -10 * length
        bounds.extend((tuple(lower), tuple(upper)))
    regions = 0
    witness = None
    for signs in product((-1, 1), repeat=len(points)):
        rows = bounds + [tuple(sign * x for x in (*point, -1))
                         for sign, point in zip(signs, points)]
        feasible = strict_cone_witness(rows, 4)
        if feasible is None:
            continue
        require(all(dot(row, feasible) > 0 for row in rows), "region substitution")
        regions += 1
        if witness is None:
            witness = {"signs": signs, "homogeneous_weights": list(map(str, feasible))}
    require(regions >= len(vertices), "regions dominate interior vertices")
    return {"points": points, "sides": sides,
            "all_sign_patterns_checked": 2 ** len(points),
            "interior_vertices": len(vertices), "nonempty_regions": regions,
            "regions_at_least_vertices": True, "region_witness": witness,
            "scope": "All cells of these eight planes in the open parameter box; not all corner cuts."}


def determinant(matrix):
    """Bareiss determinant; independent of the three-dimensional cross product."""
    n = len(matrix)
    if not n:
        return 1
    a = [list(row) for row in matrix]
    previous, sign = 1, 1
    for k in range(n - 1):
        pivot_row = next((i for i in range(k, n) if a[i][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = pivot * a[i][j] - a[i][k] * a[k][j]
                require(numerator % previous == 0, "exact Bareiss division")
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def frame_data(points):
    d = len(points)
    matrix = [[point[i] - points[0][i] for point in points[1:]] for i in range(d)]
    normal = tuple((-1) ** (i + d - 1) * determinant(matrix[:i] + matrix[i + 1:])
                   for i in range(d))
    return matrix, normal, gcd(*normal)


def rank_mod(matrix, prime):
    a = [[x % prime for x in row] for row in matrix]
    rank = 0
    for j in range(len(a[0]) if a else 0):
        pivot = next((i for i in range(rank, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inverse = pow(a[rank][j], -1, prime)
        a[rank] = [(x * inverse) % prime for x in a[rank]]
        for i in range(len(a)):
            if i != rank and a[i][j]:
                coefficient = a[i][j]
                a[i] = [(x - coefficient * y) % prime for x, y in zip(a[i], a[rank])]
        rank += 1
    return rank


def general_fixture(sides, sample_size=None):
    d = len(sides)
    denominator = 100 * d
    boxes = [tuple(product(*(range(length if i == axis else 0,
                                   (length if i == axis else 0) + length // denominator + 1)
                            for i, length in enumerate(sides)))) for axis in range(d)]
    if sample_size is None:
        tuples = product(*boxes)
    else:
        rng = Random(20261017 + d)
        selected = set()
        while len(selected) < sample_size:
            selected.add(tuple(rng.choice(box) for box in boxes))
        tuples = iter(sorted(selected))
    checked = primitive = projected_primitive = conditional_cases = congruences = 0
    seen = {}
    witness = None
    for points in tuples:
        checked += 1
        matrix, normal, index = frame_data(points)
        require(index > 0, "nonzero primitive-frame index")
        height = dot(normal, points[0])
        if height < 0:
            normal, height = tuple(-x for x in normal), -height
        require(height > 0 and all(x > 0 for x in normal), "general positive normal")
        require(all(9 * height < 10 * length * x < 11 * height
                    for length, x in zip(sides, normal)), "general parameter bounds")
        if d == 3:
            require(normal == cross(subtract(points[1], points[0]),
                                    subtract(points[2], points[0])), "cofactor/cross-product agreement")
        if index == 1:
            primitive += 1
            plane = (*normal, height)
            require(plane not in seen, "general primitive-simplex collision")
            seen[plane] = points
            if witness is None:
                witness = {"plane": plane, "points": points}
        projected_points = tuple(point[:-1] for point in points[:-1])
        _, _, projection_index = frame_data(projected_points)
        if projection_index != 1:
            continue
        projected_primitive += 1
        if conditional_cases >= 64:
            continue
        conditional_cases += 1
        top_minor = determinant(matrix[:-1])
        require(top_minor != 0, "nonzero controlling minor")
        for prime in (2, 3, 5, 7, 11, 13, 17, 19):
            if top_minor % prime:
                continue
            require(rank_mod([row[:-1] for row in matrix[:-1]], prime) == d - 2,
                    "projected rank at every tested prime")
            failures = 0
            for last in range(sides[-1], sides[-1] + sides[-1] // denominator + 1):
                changed = points[:-1] + (points[-1][:-1] + (last,),)
                test_matrix, _, test_index = frame_data(changed)
                bad = test_index % prime == 0
                require(bad == (rank_mod(test_matrix, prime) < d - 1), "rank/minor equivalence")
                failures += bad
                congruences += 1
            length = sides[-1] // denominator + 1
            require(failures <= (length + prime - 1) // prime, "one-residue-class bound")
    return {"dimension": d, "sides": sides, "box_width_denominator": denominator,
            "mode": "exhaustive" if sample_size is None else "deterministic_sample",
            "seed": None if sample_size is None else 20261017 + d,
            "tuples_checked": checked, "primitive_simplices": primitive,
            "distinct_primitive_planes": len(seen), "projected_primitive_frames": projected_primitive,
            "conditional_frames_checked": conditional_cases, "congruence_substitutions_checked": congruences,
            "primitive_witness": witness}


def explicit_sieve_bounds():
    result = []
    for d in range(2, 13):
        a, b = 100 * d, 4 * d - 3
        s = d * d + 4 * (d - 1)
        while True:
            small_error = Fraction(2 * a * (s + 1), 2 ** (2 * s)) + Fraction(a * a, 2 ** (3 * s))
            large_error = b * (Fraction(1, 2 ** s) + Fraction(a, 2 ** (2 * s)))
            if small_error + large_error <= Fraction(1, 8):
                break
            s += 1
        require(s >= d * d + 4 * (d - 1), "large-prime factor-count threshold")
        closed_s = 4 * d * d
        closed_error = (Fraction(2 * a * (closed_s + 1), 2 ** (2 * closed_s))
                        + Fraction(a * a, 2 ** (3 * closed_s))
                        + b * (Fraction(1, 2 ** closed_s)
                               + Fraction(a, 2 ** (2 * closed_s))))
        require(closed_s >= d * d + 4 * (d - 1), "closed threshold factor count")
        require(closed_error < Fraction(1, 1024), "closed uniform error bound")
        require(a <= 2 ** (4 * d) and b <= 2 ** (2 * d)
                and closed_s + 1 <= 2 ** (3 * d), "closed exponential majorants")
        r0 = 4 * d * d * (d + 1)
        # These exact inequalities audit the rounding at and after the explicit
        # threshold. Bound log_4(n) by R+1+d/2 rather than using floating logs.
        scale_checks = 0
        for r in range(r0, r0 + 4 * (d + 1)):
            least = (r + d) // (d + 1)
            residual = r - d * least
            scale_count = comb(residual + d - 1, d - 1)
            require(least >= 4 * d * d, "minimum side satisfies closed sieve threshold")
            require(r - (d - 1) * least <= 2 * least, "polynomial side comparability")
            require(residual >= Fraction(r, 2 * (d + 1)), "residual composition bound")
            logarithm_upper = Fraction(2 * r + 2 + d, 2)
            require(r >= logarithm_upper / 2, "conversion from R to logarithm")
            require(scale_count >= logarithm_upper ** (d - 1)
                    / (factorial(d - 1) * (4 * (d + 1)) ** (d - 1)),
                    "explicit cumulative lower constant")
            scale_checks += 1
        kappa_denominator = (8 ** (d - 1) * (100 * d) ** (d * d)
                             * factorial(d - 1) * (4 * (d + 1)) ** (d - 1)
                             * 2 ** (d * (d + 2)))
        result.append({"dimension": d, "s": s, "minimum_side_threshold": 2 ** (2 * s),
                       "small_prime_error_upper_bound": str(small_error),
                       "large_prime_error_upper_bound": str(large_error),
                       "conditional_failure_upper_bound": str(Fraction(3, 4) + small_error + large_error),
                       "primitive_fraction_lower_bound": str(Fraction(1, 8 ** (d - 1))),
                       "closed_s": closed_s,
                       "closed_minimum_side_threshold": 2 ** (8 * d * d),
                       "closed_error_upper_bound": str(closed_error),
                       "explicit_cumulative_bound": {
                           "minimum_n_power_of_two_exponent": d + 2 * r0,
                           "kappa": str(Fraction(1, kappa_denominator)),
                           "logarithm_base": 4,
                           "scale_rounding_cases_checked": scale_checks},
                       "scope": "Exact arithmetic in the displayed uniform bound; the written induction supplies the theorem."})
    return result


def projection_negative_control():
    points = tuple(tuple(400 * (i == j) for j in range(4)) for i in range(4))
    _, _, projection_index = frame_data(tuple(p[:-1] for p in points[:-1]))
    failures = []
    for last in (400, 401):
        changed = points[:-1] + (points[-1][:-1] + (last,),)
        _, _, index = frame_data(changed)
        failures.append(index % 2 == 0)
    require(projection_index > 1 and sum(failures) == 2, "nonprimitive projection control")
    return {"points": points, "projection_index": projection_index,
            "prime": 2, "last_coordinates": [400, 401], "failures": sum(failures),
            "one_residue_bound_without_projection_hypothesis_is_false": True}


def prime_realization_fixtures():
    """Check the prime bridge without floating logs or a precision tolerance.

    For D={x: a.x<=s}, put b=2^(d H)+1, where H=max_{x in D}|x|_1.
    In nondecreasing weight order choose p_i in (2^(i-1)b^a_i, 2^i b^a_i).
    Then inside products are below b^(s+1) and all outside products above it.
    Boundary minima and heap generation independently check the assertion.
    """
    def first_prime_after(lower):
        candidate = lower + 1
        while True:
            if candidate >= 2 and all(candidate % k for k in range(2, isqrt(candidate) + 1)):
                return candidate
            candidate += 1

    result = []
    for weights, score in (((1, 2), 3), ((2, 1), 3), ((1, 1, 1), 0),
                           ((2, 2, 2, 2), 2), ((1, 2, 3), 2), ((1, 3, 4), 2)):
        d = len(weights)
        cut = {x for x in product(*(range(score // a + 1) for a in weights))
               if dot(weights, x) <= score}
        height = max(map(sum, cut))
        base = 2 ** (d * height) + 1
        primes = [0] * d
        intervals = [None] * d
        for rank, axis in enumerate(sorted(range(d), key=weights.__getitem__)):
            lower = 2 ** rank * base ** weights[axis]
            primes[axis] = first_prime_after(lower)
            intervals[axis] = (lower, 2 * lower)
            require(primes[axis] < 2 * lower, "Bertrand interval")
        require(len(set(primes)) == d, "distinct primes even at equal weights")
        value = lambda x: prod(p ** a for p, a in zip(primes, x))
        boundary = set()
        for x in cut:
            for axis in range(d):
                y = tuple(a + (i == axis) for i, a in enumerate(x))
                if y not in cut:
                    boundary.add(y)
        last_inside = max(cut, key=value)
        first_outside = min(boundary, key=value)
        threshold = base ** (score + 1)
        require(value(last_inside) < threshold < value(first_outside), "exact product separator")
        zero = (0,) * d
        heap, seen, generated = [(1, zero)], {zero}, []
        for _ in range(len(cut) + 1):
            number, x = heapq.heappop(heap)
            generated.append(x)
            for axis, p in enumerate(primes):
                y = tuple(a + (i == axis) for i, a in enumerate(x))
                if y not in seen:
                    seen.add(y)
                    heapq.heappush(heap, (number * p, y))
        require(set(generated[:-1]) == cut and generated[-1] not in cut,
                "independent exact heap prefix")
        result.append({"integer_weights": weights, "maximum_included_score": score,
                       "base": base, "bertrand_intervals": intervals, "primes": primes,
                       "cut": sorted(cut), "cardinality": len(cut),
                       "integer_product_threshold": threshold,
                       "last_inside": {"exponents": last_inside, "value": value(last_inside)},
                       "first_outside": {"exponents": first_outside, "value": value(first_outside)},
                       "heap_prefix_matches": True})
    return result


def solve_rational(matrix, rhs):
    n = len(matrix)
    a = [[Fraction(x) for x in row] + [Fraction(y)] for row, y in zip(matrix, rhs)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        scale = a[j][j]
        a[j] = [v / scale for v in a[j]]
        for i in range(n):
            if i != j:
                scale = a[i][j]
                a[i] = [u - scale * v for u, v in zip(a[i], a[j])]
    return [row[-1] for row in a]


def bezout_vector(values):
    g, coefficients = 0, []
    for value in values:
        old_r, r, old_s, s, old_t, t = g, value, 1, 0, 0, 1
        while r:
            q = old_r // r
            old_r, r = r, old_r - q * r
            old_s, s = s, old_s - q * s
            old_t, t = t, old_t - q * t
        g = old_r
        coefficients = [old_s * c for c in coefficients] + [old_t]
    require(g == 1 and dot(values, coefficients) == 1, 'Bezout substitution')
    return coefficients


def reduce_lattice_coset(seed, target, matrix):
    coefficients = solve_rational(matrix[1:], [Fraction(x) - y for x, y in zip(seed[1:], target[1:])])
    nearest = [int((c + Fraction(1, 2)) // 1) for c in coefficients]
    residual = [c - k for c, k in zip(coefficients, nearest)]
    require(all(abs(c) <= Fraction(1, 2) for c in residual), 'centered coefficients')
    reduced = tuple(x - dot(row, nearest) for x, row in zip(seed, matrix))
    require(all(Fraction(x) - y == dot(row, residual)
                for x, y, row in zip(reduced, target, matrix)), 'full coset substitution')
    return reduced, residual


def primitive_normal_family(d, count=24):
    # Given primitive frames in the stated windows; this finite fixture does
    # not assert the density lemma at these small side lengths.
    sides = tuple(100 * d * (i + 1) for i in range(d))
    rng = Random(20261017 + 100 * d)
    frames = {}
    attempts = 0
    while len(frames) < count:
        attempts += 1
        require(attempts < 20000, 'enough distinct fixture normals')
        points = tuple(tuple((length if axis == j else 0) + rng.randrange(length // (100 * d) + 1)
                             for j, length in enumerate(sides)) for axis in range(d))
        matrix, a, index = frame_data(points)
        if index != 1:
            continue
        height = dot(a, points[0])
        if height < 0:
            a, height = tuple(-v for v in a), -height
        require(height == abs(determinant(points)), 'primitive normal height equals determinant')
        require(all(v > 0 for v in a), 'positive normal')
        require(Fraction(prod(sides), 2 ** d) <= height <= 2 ** d * prod(sides), 'height range')
        require(all(Fraction(9, 10) < Fraction(x * length, height) < Fraction(11, 10)
                    for x, length in zip(a, sides)), 'scaled normal')
        frames[a] = (points, matrix, height)
    return sides, frames, attempts


def fixed_size_patch_fixture(d):
    sides, frames, attempts = primitive_normal_family(d, 9 if d == 2 else 24)
    k = 4 ** (d + 1)
    ratio = 16 * d * d * (k + 1)
    pairs = 0
    witness = None
    for a, (points, matrix, height) in frames.items():
        h = ratio * height
        z, c = reduce_lattice_coset(bezout_vector(a), [Fraction(v, height) for v in points[0]], matrix)
        require(dot(a, z) == 1, 'unit level increment')
        target = [Fraction(h - 1, d * v) for v in a]
        x0, beta = reduce_lattice_coset([(h - 1) * v for v in z], target, matrix)
        for signs in product((-k, k), repeat=d - 1):
            low = tuple(x + dot(row, signs) for x, row in zip(x0, matrix))
            high = tuple(x + 2 * v for x, v in zip(low, z))
            require(min(*low, *high) >= 0, 'two nonnegative patches')
            require(dot(a, low) == h - 1 and dot(a, high) == h + 1, 'patch levels')
        for b in frames:
            if b == a:
                continue
            pairs += 1
            columns = tuple(zip(*matrix))
            values = [dot(b, column) for column in columns]
            require(any(values), 'nonparallel normal')
            signs = [k if v > 0 else -k if v < 0 else 0 for v in values]
            low = tuple(x + dot(row, signs) for x, row in zip(x0, matrix))
            high = tuple(x + 2 * v - dot(row, signs) for x, v, row in zip(x0, z, matrix))
            require(dot(b, low) > dot(b, high), 'distinct normal cannot separate patches')
            require(Fraction(dot(b, points[0]), height) < 2 * 4 ** d, 'uniform height-ratio bound')
            if witness is None:
                witness = {'a': a, 'b': b, 'height': height, 'boundary_score': h,
                           'unit_increment': z, 'low_point': low, 'high_point': high,
                           'a_low': dot(a, low), 'a_high': dot(a, high),
                           'b_low': dot(b, low), 'b_high': dot(b, high)}
    return {'dimension': d, 'sides': sides, 'attempts': attempts,
            'seed': 20261017 + 100 * d,
            'primitive_normal_count': len(frames), 'ordered_distinct_pairs': pairs,
            'patch_corner_pairs_checked': len(frames) * 2 ** (d - 1),
            'patch_half_width': k, 'boundary_height_ratio': ratio, 'witness': witness}


def exact_size_constants():
    result = []
    for d in range(2, 13):
        k = 4 ** (d + 1)
        scale_root = 64 * d * d * (k + 1)
        scale = scale_root ** d
        ratio = 16 * d * d * (k + 1)
        require(Fraction(ratio - 1, 2 * d) > 2 * k * d + 5 * d,
                'uniform patch nonnegativity margin')
        require(k - Fraction(1, 2) > 2 * 4 ** d,
                'uniform integer rigidity margin')
        require(Fraction(4, 5) * scale_root >= ratio,
                'prescribed cardinality supplies sufficient height')
        require(scale <= 4 ** (2 * d * d + 4 * d), 'closed scale logarithm bound')
        r0 = 4 * d * d * (d + 1)
        require(r0 >= 2 * d * d + 4 * d + 2, 'closed logarithm conversion threshold')
        rho = Fraction(1, 8 ** (d - 1) * (100 * d) ** (d * d))
        constant = rho / (2 ** d * (4 * scale) ** (d - 1)
                          * (4 * (d + 1)) ** (d - 1) * factorial(d - 1))
        result.append({'dimension': d, 'patch_half_width': k,
                       'cardinality_to_box_volume_ratio': scale,
                       'minimum_n': {'factor': scale, 'power_of_four': r0},
                       'pointwise_lower_constant': str(constant), 'logarithm_base': 4})
    return result


def integer_weight_prefix(weights, n):
    zero = (0,) * len(weights)
    heap, seen, result = [(0, zero)], {zero}, []
    for _ in range(n):
        score, x = heapq.heappop(heap)
        result.append((score, x))
        for axis, weight in enumerate(weights):
            y = tuple(v + (i == axis) for i, v in enumerate(x))
            if y not in seen:
                seen.add(y)
                heapq.heappush(heap, (score + weight, y))
    return result


def tie_resolution_fixture():
    checks, examples = 0, []
    for normal in ((1, 1), (2, 3), (1, 1, 1), (2, 2, 3), (1, 2, 3, 4)):
        d = len(normal)
        for n in range(1, 41):
            ordered = integer_weight_prefix(normal, n)
            u = tuple(n ** (d - i - 1) for i in range(d))
            multiplier = n ** d
            strict_weights = tuple(multiplier * a + v for a, v in zip(normal, u))
            strict_prefix = integer_weight_prefix(strict_weights, n + 1)
            require([x for _, x in strict_prefix[:-1]] == [x for _, x in ordered],
                    'lexicographic prefix realized by positive integer weights')
            require(strict_prefix[-2][0] < strict_prefix[-1][0],
                    'strict score gap at prescribed cardinality')
            checks += 1
            if n == 7:
                examples.append({'normal': normal, 'cardinality': n, 'strict_weights': strict_weights,
                                 'last_inside': strict_prefix[-2], 'first_outside': strict_prefix[-1]})
    return {'prefixes_checked': checks, 'cardinalities': [1, 40],
            'exact_heap_orders_agree': True, 'witnesses': examples}



def main():
    fixtures = [primitive_fixture(sides) for sides in
                ((100, 100, 100), (100, 200, 400), (200, 300, 400))]
    report = {
        "status": "finite_falsification_and_constant_checks_passed_priority_unresolved",
        "command": "python3 -B " + shlex.join([str(Path(__file__).relative_to(ROOT)), *sys.argv[1:]]),
        "primitive_fixtures": fixtures,
        "arrangement_fixture": arrangement_fixture(),
        "general_fixtures": [general_fixture(sides, sample_size) for sides, sample_size in (
            ((200, 400), None), ((300, 600, 900), None),
            ((400, 400, 400, 400), None), ((500, 1000, 1500, 500, 1000), 1000),
            ((600, 1200, 600, 1200, 600, 1200), 500))],
        "explicit_uniform_sieve_bounds": explicit_sieve_bounds(),
        "nonprimitive_projection_control": projection_negative_control(),
        "prime_realization_fixtures": prime_realization_fixtures(),
        "exact_size_patch_fixtures": [fixed_size_patch_fixture(d) for d in range(2, 7)],
        "exact_size_constants": exact_size_constants(),
        "tie_resolution_fixture": tie_resolution_fixture(),
        "source_sha256": {str(path): hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
                          for path in (Path("scripts/run_corner_cut_growth_checks.py"),
                                       Path("scripts/run_pyramidal_poset_checks.py"),
                                       Path("docs/interleaving_research.md"),
                                       Path("paper/corner_cut_growth/main.tex"))},
        "limitations": ["Finite checks do not establish a uniform asymptotic density bound.",
                        "The prescribed-size conclusion depends on the new lattice-patch proof, not the cumulative checks alone.",
                        "No claim of priority or publication readiness is made."],
    }
    target = ROOT / "results/research/corner_cut_growth_checks.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"report": str(target.relative_to(ROOT)),
                      "triples_checked": sum(x["ordered_triples_checked"] for x in fixtures),
                      "primitive_planes": sum(x["distinct_primitive_planes"] for x in fixtures),
                      "general_fixtures": [{k: x[k] for k in ("dimension", "mode", "tuples_checked",
                                               "primitive_simplices", "congruence_substitutions_checked")}
                                           for x in report["general_fixtures"]],
                      "arrangement": report["arrangement_fixture"]}, indent=2))


if __name__ == "__main__":
    main()
