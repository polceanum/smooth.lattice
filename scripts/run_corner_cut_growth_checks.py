#!/usr/bin/env python3
"""Exact falsification checks for the multiscale corner-cut growth argument.

These finite checks do not prove the asymptotic density lemma or establish
novelty. They check primitive-triangle injectivity, rational parameter bounds,
and the vertex/region inequality on a small degenerate arrangement.
"""
from __future__ import annotations

import hashlib
from itertools import combinations, product
import json
from math import gcd, prod
from pathlib import Path
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


def main():
    fixtures = [primitive_fixture(sides) for sides in
                ((100, 100, 100), (100, 200, 400), (200, 300, 400))]
    report = {
        "status": "finite_falsification_checks_passed_asymptotic_proof_and_priority_require_audit",
        "command": "python3 -B " + shlex.join([str(Path(__file__).relative_to(ROOT)), *sys.argv[1:]]),
        "primitive_fixtures": fixtures,
        "arrangement_fixture": arrangement_fixture(),
        "source_sha256": {str(path): hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
                          for path in (Path("scripts/run_corner_cut_growth_checks.py"),
                                       Path("scripts/run_pyramidal_poset_checks.py"),
                                       Path("docs/interleaving_research.md"))},
        "limitations": ["Finite checks do not establish a uniform asymptotic density bound.",
                        "The cumulative three-dimensional candidate does not establish the exact-size three-dimensional bound.",
                        "No claim of priority or publication readiness is made."],
    }
    target = ROOT / "results/research/corner_cut_growth_checks.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"report": str(target.relative_to(ROOT)),
                      "triples_checked": sum(x["ordered_triples_checked"] for x in fixtures),
                      "primitive_planes": sum(x["distinct_primitive_planes"] for x in fixtures),
                      "arrangement": report["arrangement_fixture"]}, indent=2))


if __name__ == "__main__":
    main()
