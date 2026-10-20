#!/usr/bin/env python3
"""Exact controls for the classical generic-multiaffine screening bound.

Finite checks support docs/generic_multiaffine_bound.md. They do not prove
the asymptotic bound, its priority, or the unrestricted extremal endpoint.
Only Python's standard library is required.
"""
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import hashlib
import json

from run_extremal_initial_ideals_checks import integer_weight
from run_interleaving_research import minimal_excluded_border
from run_structural_complexity_checks import (
    area2, clip, cross, dot, hyperbolic_cross, positive_chambers, primitive,
)

ROOT = Path(__file__).resolve().parents[1]
PRIME = 1000003


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def subtract(a, b):
    return tuple(x-y for x, y in zip(a, b))


def rank(columns, prime=PRIME):
    pivots = {}
    for column in columns:
        v = [x % prime for x in column]
        for i, row in sorted(pivots.items()):
            c = v[i]
            v = [(x-c*y) % prime for x, y in zip(v, row)]
        i = next((i for i, x in enumerate(v) if x), None)
        if i is not None:
            inverse = pow(v[i], -1, prime)
            pivots[i] = [x*inverse % prime for x in v]
    return len(pivots)


def choose(ground, n, weight):
    return frozenset(sorted(ground, key=lambda a: (dot(a, weight), a))[:n])


def canonical(a, head):
    return any(a[i] == 0 for i in head)


def uniform_fan(n, head, exhaustive=False):
    ground = tuple(a for a in hyperbolic_cross(n, 3) if canonical(a, head))
    if n == 1:
        return ground, {frozenset(ground): (
            [tuple(Q(i == j) for i in range(3)) for j in range(3)], (1, 1, 1))}, 0
    pending = [choose(ground, n, (1, n+1, (n+1)**2))]
    visited = {}
    while pending:
        basis = pending.pop()
        if basis in visited:
            continue
        maximal = [a for a in basis if not any(
            tuple(x+(i == j) for i, x in enumerate(a)) in basis for j in range(3))]
        border = [a for a in minimal_excluded_border(sorted(basis)) if canonical(a, head)]
        normals = sorted({tuple(x//gcd(*z) for x in z)
                          for a in border for b in maximal if (z := subtract(a, b))})
        polygon = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for z in normals:
            polygon = clip(polygon, z, 1)
            require(polygon, "a visited basis has an open positive cone")
        weight = integer_weight(tuple(sum(p[i] for p in polygon)/len(polygon)
                                      for i in range(3)))
        require(choose(ground, n, weight) == basis and
                max(dot(weight, a) for a in basis) <
                min(dot(weight, a) for a in ground if a not in basis),
                "centroid realizes a strict prefix")
        visited[basis] = (polygon, weight)
        for p, q in zip(polygon, polygon[1:]+polygon[:1]):
            z = next((z for z in normals if dot(z, p) == dot(z, q) == 0), None)
            if z is None or any(p[i] == q[i] == 0 for i in range(3)):
                continue
            # An edge interior point avoids all unrelated comparison walls.
            for k in range(1, len(ground)**2+2):
                mid = integer_weight(tuple((k*x+y)/(k+1) for x, y in zip(p, q)))
                groups = {}
                good = True
                for a in ground:
                    cost = dot(mid, a)
                    if cost in groups and cross(z, subtract(a, groups[cost])) != (0, 0, 0):
                        good = False
                        break
                    groups[cost] = a
                if good:
                    break
            else:
                raise AssertionError("finite generic-edge avoidance")
            # Exact symbolic infinitesimal, rather than a numerical epsilon.
            neighbor = frozenset(sorted(ground, key=lambda a:
                (dot(mid, a), -dot(z, a), a))[:n])
            require(neighbor != basis, "crossing an active facet changes the prefix")
            if neighbor not in visited:
                pending.append(neighbor)
    require(sum(area2(poly) for poly, _ in visited.values()) == 1,
            "exact complete coverage of the weight simplex")
    comparisons = 0
    if exhaustive:
        normals = sorted({primitive(subtract(a, b)) for a, b in combinations(ground, 2)
                          if min(subtract(a, b)) < 0 < max(subtract(a, b))})
        weights = positive_chambers(normals)
        require({choose(ground, n, w) for w in weights} == set(visited),
                "independent full comparison arrangement versus facet traversal")
        comparisons = len(weights)
    return ground, visited, comparisons


def boundary_perturbation(ground, fan):
    maximum = max(max(a) for a in ground)
    scale = 2+18*maximum**3+3*maximum**2*max(sum(w) for _, w in fan.values())
    lifted = {a: tuple(scale*x+dot(a, a) for x in a) for a in ground}
    margins = []
    supports = {}
    for a in ground:
        i = a.index(0)
        normal = tuple(scale*(i == j)-2*a[j] for j in range(3))
        margins.extend(dot(normal, subtract(lifted[b], lifted[a])) for b in ground if b != a)
        supports[a] = normal
    require(not margins or min(margins) > 0, "every perturbed point is strictly exposed")
    for basis, (_, w) in fan.items():
        require(choose(tuple(lifted.values()), len(basis), w) ==
                frozenset(lifted[a] for a in basis), "every original strict prefix survives")
    return lifted, supports, scale, min(margins, default=0)


def kset_controls(points):
    """Independently enumerate oriented triples and strict subset witnesses."""
    m = len(points)
    facets = [0]*(m-2)
    sets = [set() for _ in range(m+1)]
    for triple in combinations(range(m), 3):
        p, q, r = (points[i] for i in triple)
        u, v = subtract(q, p), subtract(r, p)
        normal = cross(u, v)
        require(normal != (0, 0, 0), "noncollinear triples")
        require(all(dot(normal, subtract(a, p)) for i, a in enumerate(points) if i not in triple),
                "general position independently checked on every triple")
        for sign in (-1, 1):
            z = tuple(sign*x for x in normal)
            lower = frozenset(i for i, a in enumerate(points) if dot(z, subtract(a, p)) < 0)
            facets[len(lower)] += 1
            for count in (1, 2):
                for selected in combinations(triple, count):
                    target = [int(i not in selected) for i in triple]
                    du, dv = target[1]-target[0], target[2]-target[0]
                    uu, uv, vv = dot(u, u), dot(u, v), dot(v, v)
                    det = uu*vv-uv*uv
                    d = tuple((du*vv-dv*uv)*x+(dv*uu-du*uv)*y for x, y in zip(u, v))
                    require(dot(d, u) == det*du and dot(d, v) == det*dv,
                            "independent tangent separator for the threshold triangle")
                    scale = 1+2*max(abs(dot(d, a)) for a in points)
                    w = tuple(scale*x+y for x, y in zip(z, d))
                    basis = lower | frozenset(selected)
                    ordered = sorted(range(m), key=lambda i: (dot(w, points[i]), i))
                    require(frozenset(ordered[:len(basis)]) == basis and
                            dot(w, points[ordered[len(basis)-1]]) < dot(w, points[ordered[len(basis)]]),
                            "every enumerated set has an exact strict separating weight")
                    sets[len(basis)].add(basis)
    require(facets == [2*(j+1)*(m-j-2) for j in range(m-2)],
            "classical convex-position oriented-facet formula")
    counts = [len(sets[r]) for r in range(1, m)]
    require(counts == [2*r*(m-r)-m+2 for r in range(1, m)],
            "independent strict-set enumeration versus the classical formula")
    return {"points": points, "oriented_facets_by_lower_count": facets,
            "strict_set_counts": counts, "every_set_has_strict_witness": True}


def general_position_control(ground, fan):
    lifted, supports, scale, margin = boundary_perturbation(ground, fan)
    original = tuple(ground)
    jitter = [(i, i*i, i*i*i) for i in range(1, len(ground)+1)]
    scale2 = 2+2*max(abs(dot(z, subtract(a, b)))
                    for z in list(supports.values())+[w for _, w in fan.values()]
                    for a, b in combinations(jitter, 2))
    for candidate in range(scale2, scale2+3*len(tuple(combinations(ground, 4)))+1):
        points = [tuple(candidate*x+y for x, y in zip(lifted[a], jitter[i]))
                  for i, a in enumerate(original)]
        if all(dot(cross(subtract(points[j], points[i]), subtract(points[k], points[i])),
                   subtract(points[l], points[i]))
               for i, j, k, l in combinations(range(len(points)), 4)):
            break
    else:
        raise AssertionError("nonzero moment-curve perturbation determinants")
    for a, z in supports.items():
        i = original.index(a)
        require(all(dot(z, subtract(b, points[i])) > 0 for j, b in enumerate(points) if i != j),
                "convex position survives the general-position perturbation")
    for basis, (_, w) in fan.items():
        require(choose(tuple(points), len(basis), w) ==
                frozenset(points[original.index(a)] for a in basis),
                "all strict prefixes survive both perturbations")
    result = kset_controls(points)
    result.update({"boundary_scale": scale, "general_position_scale": candidate,
                   "original_exposure_margin": margin})
    return result


def multiply_laurent(a, b):
    answer = {}
    for i, x in a.items():
        for j, y in b.items():
            answer[i+j] = answer.get(i+j, 0)+x*y
    return {i: x for i, x in answer.items() if x}


def symbolic_vandermonde(exponents):
    determinant = {}
    n = len(exponents)
    for order in permutations(range(n)):
        degree = sum(i*exponents[j] for i, j in enumerate(order))
        sign = (-1)**sum(order[i] > order[j] for i in range(n) for j in range(i+1, n))
        determinant[degree] = determinant.get(degree, 0)+sign
    determinant = {i: x for i, x in determinant.items() if x}
    factored = {0: 1}
    for i, j in combinations(range(n), 2):
        factored = multiply_laurent(factored, {exponents[j]: 1, exponents[i]: -1})
    require(determinant == factored and abs(determinant[max(determinant)]) == 1,
            "permutation determinant equals Vandermonde product with unit extreme coefficient")
    for characteristic in (2, 3):
        require(any(x % characteristic for x in determinant.values()),
                "symbolic determinant remains nonzero in small characteristic")


def root_of_unity_order():
    # Complete factorization of PRIME-1, checked without an external library.
    factors, residue, trial = [], PRIME-1, 2
    while trial*trial <= residue:
        if residue % trial == 0:
            factors.append(trial)
            while residue % trial == 0:
                residue //= trial
        trial += 1
    if residue > 1:
        factors.append(residue)
    generator = next(t for t in range(2, PRIME)
                     if all(pow(t, (PRIME-1)//p, PRIME) != 1 for p in factors))
    return generator


def polynomial_model(n, T, B):
    roots = [pow(T, i, PRIME) for i in range(n)]
    f = [1]
    for root in roots:
        answer = [0]*(len(f)+1)
        for i, c in enumerate(f):
            answer[i] = (answer[i]-root*c) % PRIME
            answer[i+1] = (answer[i+1]+c) % PRIME
        f = answer

    def mul(a, b):
        answer = [0]*(2*n-1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                answer[i+j] = (answer[i+j]+x*y) % PRIME
        for k in range(2*n-2, n-1, -1):
            c = answer[k]
            for j in range(n):
                answer[k-n+j] = (answer[k-n+j]-c*f[j]) % PRIME
        return tuple(answer[:n])

    one = (1,)+(0,)*(n-1)
    x = (0, 1)+(0,)*(n-2)
    inverse = tuple(-c*pow(f[0], -1, PRIME) % PRIME for c in f[1:])
    require(mul(x, inverse) == one, "independent inverse x modulo the point polynomial")
    powers = {0: one}
    for k in range(1, (B+1)*(n-1)+1):
        powers[k] = mul(powers[k-1], x)
        powers[-k] = mul(powers[1-k], inverse)

    def column(a):
        return mul(mul(powers[a[0]], powers[B*a[1]]), powers[-(B+1)*a[2]])

    def evaluate(v, t):
        answer = 0
        for c in reversed(v):
            answer = (answer*t+c) % PRIME
        return answer

    return column, lambda v: tuple(evaluate(v, t) for t in roots), f


def realization_controls(n, ground, fan, T):
    B = 2*n+1
    exponents = {a: a[0]-a[2]+B*(a[1]-a[2]) for a in ground}
    require(len(set(exponents.values())) == len(ground), "Kronecker exponent injectivity")
    require(PRIME-1 > 2*(B+1)*(n-1), "finite-field order exceeds the full exponent width")
    column, evaluate, polynomial = polynomial_model(n, T, B)
    columns = {a: column(a) for a in ground}
    for a in ground:
        require(evaluate(columns[a]) == tuple(pow(T, i*exponents[a], PRIME) for i in range(n)),
                "polynomial remainder evaluation versus direct geometric columns")
    symbolic = 0
    if n <= 4:
        for subset in combinations(ground, n):
            require(rank([columns[a] for a in subset]) == n, "all small quotient minors")
            symbolic_vandermonde([exponents[a] for a in subset])
            symbolic += 1
    for basis in fan:
        require(rank([columns[a] for a in basis]) == n, "every coherent prefix is a quotient basis")
    return {"length": n, "characteristic": PRIME, "T": T, "B": B,
            "point_polynomial": polynomial, "canonical_columns_checked": len(ground),
            "coherent_bases_checked": len(fan), "all_small_symbolic_minors_checked": symbolic,
            "small_characteristics_checked_symbolically": [2, 3] if symbolic else [],
            "points": [[pow(T, i, PRIME), pow(T, B*i, PRIME), pow(T, -(B+1)*i, PRIME)]
                       for i in range(n)]}


def special_point_control():
    values = (1, 2, 3, 4)
    points = [(t, t, pow(t*t, -1, PRIME)) for t in values]
    ground = tuple(a for a in hyperbolic_cross(4, 3) if min(a) == 0)
    prefix = sorted(ground, key=lambda a: (dot(a, (2, 3, 5)), a))[:4]
    columns = [tuple(pow(p[0], a[0], PRIME)*pow(p[1], a[1], PRIME)*pow(p[2], a[2], PRIME)
                     % PRIME for p in points) for a in prefix]
    require(prefix == [(0, 0, 0), (1, 0, 0), (0, 1, 0), (2, 0, 0)] and rank(columns) == 3,
            "a reduced special configuration on xyz=1 violates the generic rank rule")
    require(all(p[0]*p[1]*p[2] % PRIME == 1 for p in points), "same squarefree carrier")
    return {"points": points, "weight": [2, 3, 5], "first_four_canonical": prefix,
            "prefix_rank": 3, "quotient_length": 4,
            "conclusion": "Containing the carrier does not imply generic evaluation ranks."}


def evaluation_columns(ground, points):
    return {a: tuple(pow(p[0], a[0], PRIME)*pow(p[1], a[1], PRIME)*pow(p[2], a[2], PRIME)
                     % PRIME for p in points) for a in ground}


def greedy_evaluations(ground, columns, n, weight):
    selected = []
    for a in sorted(ground, key=lambda a: (dot(a, weight), a)):
        if rank([columns[b] for b in selected]+[columns[a]]) > len(selected):
            selected.append(a)
        if len(selected) == n:
            break
    require(len(selected) == n, "evaluation quotient has the claimed length")
    return frozenset(selected)


def permute_fan(model, head):
    order = list(head)+[i for i in range(3) if i not in head]

    def permute(a):
        output = [0]*3
        for i, j in enumerate(order):
            output[j] = a[i]
        return tuple(output)

    ground, fan = model
    return tuple(map(permute, ground)), {
        frozenset(map(permute, basis)): ([permute(p) for p in poly], permute(w))
        for basis, (poly, w) in fan.items()}


def common_list(n, fans):
    boolean = tuple(product((0, 1), repeat=3))
    normals = sorted({primitive(subtract(a, b)) for a, b in combinations(boolean, 2)
                      if min(subtract(a, b)) < 0 < max(subtract(a, b))})
    require(len(normals) == 6, "fixed squarefree comparison arrangement")
    weights, pieces = set(), 0
    for t in (1, 2, 3):
        for head in combinations(range(3), t):
            ground, fan = permute_fan(fans[n, t], head)
            for basis, (poly, _) in fan.items():
                cells = [poly]
                for z in normals:
                    cells = [part for cell in cells for sign in (-1, 1)
                             if (part := clip(cell, z, sign))]
                require(sum(area2(cell) for cell in cells) == area2(poly),
                        "Boolean refinement preserves each prefix cone exactly")
                pieces += len(cells)
                for cell in cells:
                    w = integer_weight(tuple(sum(p[i] for p in cell)/len(cell) for i in range(3)))
                    require(choose(ground, n, w) == basis and all(dot(w, z) for z in normals),
                            "common-list weight is strict in both arrangements")
                    weights.add(w)
    return sorted(weights), pieces


def carrier_controls(n, fans, T, common_weights):
    all_ground = tuple(hyperbolic_cross(n, 3))
    B = 2*n+1
    xy = [(pow(T, i, PRIME), pow(T, B*i, PRIME)) for i in range(n)]
    points = [(x, y, (x+y+1)*pow(x*y, -1, PRIME) % PRIME) for x, y in xy]
    require(all((x*y*z-x-y-1) % PRIME == 0 for x, y, z in points),
            "nonbinomial multiaffine carrier vanishes")
    columns = evaluation_columns(all_ground, points)
    expected = set(fans[n, 3][1])
    for basis in expected:
        require(rank([columns[a] for a in basis]) == n,
                "every uniform prefix is independent for the explicit nonbinomial carrier")
    observed = {greedy_evaluations(all_ground, columns, n, w) for w in common_weights}
    require(observed == expected, "same complete fan on a nonbinomial dominant-head surface")
    nonbinomial = {"carrier": "xyz-x-y-1", "points": points,
                  "all_canonical_prefix_ranks_checked": len(expected),
                  "common_list_complete_fan_count": len(observed)}

    points = [(x, y, x*y % PRIME) for x, y in xy]
    columns = evaluation_columns(all_ground, points)
    expected, area, pieces = set(), Q(0), 0
    for head, sign in (((0, 1), 1), ((2,), -1)):
        ground, fan = permute_fan(fans[n, len(head)], head)
        for basis, (poly, _) in fan.items():
            part = clip(poly, (1, 1, -1), sign)
            if not part:
                continue
            area += area2(part)
            pieces += 1
            w = integer_weight(tuple(sum(p[i] for p in part)/len(part) for i in range(3)))
            require(greedy_evaluations(all_ground, columns, n, w) == basis,
                    "direct quotient greedy agrees across both possible carrier heads")
            expected.add(basis)
    require(area == 1, "both carrier-head regions cover the simplex exactly")
    observed = {greedy_evaluations(all_ground, columns, n, w) for w in common_weights}
    require(observed == expected, "one common list captures the entire changing-head fan")
    return [nonbinomial, {"carrier": "xy-z", "points": points, "carrier_heads": [[1, 1, 0], [0, 0, 1]],
                         "complete_carrier_refinement_pieces": pieces,
                         "normalized_double_area": "1", "complete_fan_count": len(expected),
                         "common_list_complete_fan_count": len(observed)}]


def main():
    require(PRIME == 1000003 and all(PRIME % i for i in range(2, 1001)), "prime modulus")
    records, realizations, perturbations = [], [], []
    generator = root_of_unity_order()
    fans = {}
    for n in (2, 3, 4, 5, 6, 8, 10, 12):
        for t in (1, 2, 3) if n <= 8 else (3,):
            head = tuple(range(t))
            ground, fan, chambers = uniform_fan(n, head, exhaustive=(n <= 3 or n == 4 and t == 3))
            divisor_sum = sum(n//j for j in range(1, n+1))
            expected = t*divisor_sum-(t*(t-1)//2)*n+int(t == 3)
            require(len(ground) == expected, "independent inclusion-exclusion ground count")
            bound = 2*n*(len(ground)-n)-len(ground)+2 if len(ground) > n else 1
            require(len(fan) <= bound, "finite boundary bound control")
            lifted, supports, scale, margin = boundary_perturbation(ground, fan)
            perturbations.append({"n": n, "head": head, "ground_size": len(ground),
                                  "scale": scale, "minimum_exposure_margin": margin,
                                  "prefixes_preserved": len(fan)})
            records.append({"n": n, "head": head, "ground_size": len(ground),
                            "complete_fan_count": len(fan), "classical_bound": bound,
                            "independent_comparison_chambers": chambers,
                            "normalized_double_area": "1",
                            "witnesses": [{"basis": sorted(basis), "weight": w,
                                           "polygon": [[str(x) for x in p] for p in poly]}
                                          for basis, (poly, w) in sorted(fan.items(), key=lambda item: sorted(item[0]))]})
            fans[n, t] = (ground, fan)
            if t == 3:
                realizations.append(realization_controls(n, ground, fan, generator))
    convex = [kset_controls([(i, i*i, i*i*i) for i in range(1, m+1)]) for m in range(4, 11)]
    convex.append(general_position_control(*fans[3, 3]))
    special = special_point_control()
    lists, carriers = [], []
    for n in (4, 8):
        weights, pieces = common_list(n, fans)
        lists.append({"n": n, "number_of_weights": len(weights),
                      "refinement_pieces": pieces, "weights": weights})
        carriers.extend([{"n": n, **record} for record in carrier_controls(n, fans, generator, weights)])
    sources = ["docs/generic_multiaffine_bound.md", "scripts/run_generic_multiaffine_checks.py",
               "scripts/run_structural_complexity_checks.py", "scripts/run_interleaving_research.py",
               "scripts/run_extremal_initial_ideals_checks.py"]
    report = {"schema_version": 1,
              "command": "python3 -B scripts/run_generic_multiaffine_checks.py",
              "scope": "Finite exact controls for a uniform restricted theorem proved in the note. The counting input is classical convex-position geometry. No growth fit, novelty certificate or unrestricted endpoint claim.",
              "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              "fans": records, "boundary_perturbations": perturbations,
              "independent_convex_position_counts": convex,
              "explicit_point_quotients": realizations, "special_configuration_guard": special,
              "common_weight_lists": lists, "additional_carrier_controls": carriers,
              "summary": {"complete_fans": len(records), "complete_cells": sum(len(v[1]) for v in fans.values()),
                          "independent_comparison_chambers": sum(r["independent_comparison_chambers"] for r in records),
                          "explicit_point_quotients": len(realizations),
                          "symbolic_vandermonde_minors": sum(r["all_small_symbolic_minors_checked"] for r in realizations),
                          "convex_position_configurations": len(convex),
                          "common_lists": len(lists), "additional_carrier_point_quotients": len(carriers)}}
    path = ROOT/"results/research/generic_multiaffine_checks.json"
    path.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
