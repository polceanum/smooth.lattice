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
    x_budget = n * (r - 1) // 2
    y_budget = n * (s - 1) // 2
    degree_budget = x_budget + y_budget

    def check_moment_budget(basis):
        moment = [sum(a[i] for a in basis) for i in range(3)]
        require(moment[0] <= x_budget and moment[1] <= y_budget and
                sum(moment) <= degree_budget, "uniform determinant/filtration moment bounds")

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
        check_moment_budget(basis)
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
                check_moment_budget(candidate)
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
            "state_moment_bounds": {"x": x_budget, "y": y_budget, "total_degree": degree_budget,
                                    "all_visited_bases_checked": True,
                                    "every_determinant_basis_checked_when_exhaustive": exhaustive_control},
            "independent_determinants_checked": independent, "witnesses": witnesses}


def rational_probe_controls():
    """A uniform hidden relation distinguishes the fixed probe from generic F."""
    d = {(i, j): (-1)**(i+j) * comb(3, i) * comb(2, j)
         for i in range(4) for j in range(3)}
    numerator = {(0, 0): 4, (1, 0): -5, (2, 0): 3,
                 (0, 1): -1, (1, 1): 4, (2, 1): -5}
    constant_part = {p: 4*d.get(p, 0) - numerator.get(p, 0) for p in set(d) | set(numerator)}
    constant_part = {p: c for p, c in constant_part.items() if c}
    require(len(d) == 12 and len(constant_part) == 11, "fixed relation has 23 monomials")
    exponents = tuple(sorted([(i, j, 0) for i, j in constant_part] +
                             [(i, j, 1) for i, j in d]))
    relation = {**{(i, j, 0): c for (i, j), c in constant_part.items()},
                **{(i, j, 1): c for (i, j), c in d.items()}}
    fixtures = []
    for r, s in ((8, 6), (9, 7), (12, 8)):
        cells = tuple(product(range(r), range(s)))
        fixed_f = {(i, j): (i+1)**2 + 3*(j+1) + 5*i*j for i, j in cells if i or j}
        generic_witness_f = {(1, 0): 1, (0, 1): 1, (4, 3): 1}

        def column(exponent, polynomial):
            a, b, c = exponent
            require(c in (0, 1), "relation-control columns need at most one z factor")
            return tuple(int((i, j) == (a, b)) if c == 0 else polynomial.get((i-a, j-b), 0)
                         for i, j in cells)

        fixed_columns = {e: column(e, fixed_f) for e in exponents}
        alternative_columns = {e: column(e, generic_witness_f) for e in exponents}
        require(all(sum(relation[e] * fixed_columns[e][k] for e in exponents) == 0
                    for k in range(len(cells))), "fixed rational-function relation over the integers")
        require(any(sum(relation[e] * alternative_columns[e][k] for e in exponents)
                    for k in range(len(cells))), "the relation is not an identity for arbitrary F")
        require(all(prod(a+1 for a in e) <= r*s for e in exponents),
                "all diagnostic columns belong to the universal monomial ground")

        # On these 23 rows the alternative columns form [I, B; 0, I].
        # The unit determinant is checked entrywise over the integers, not
        # inferred from modular rank. Replacing x^4*y^3 by t*x^4*y^3 gives t^12.
        low = sorted(constant_part)
        high = [(i+4, j+3) for i, j in sorted(d)]
        ordered_columns = [(i, j, 0) for i, j in low] + [(i, j, 1) for i, j in sorted(d)]
        chosen_rows = low + high
        require(len(set(chosen_rows)) == 23, "disjoint low and shifted rectangular rows")
        selected = [[alternative_columns[e][cells.index(p)] for e in ordered_columns] for p in chosen_rows]
        for i, row in enumerate(selected):
            for j, entry in enumerate(row):
                if i >= 11 or j < 11:
                    require(entry == int(i == j), "explicit unit triangular minor over every field")
        prime = 101
        fixed_mod = {e: tuple(v % prime for v in c) for e, c in fixed_columns.items()}
        alternative_mod = {e: tuple(v % prime for v in c) for e, c in alternative_columns.items()}
        fixed_rank = len(greedy(exponents, fixed_mod, (1, 37, 1369), prime))
        alternative_rank = len(greedy(exponents, alternative_mod, (1, 37, 1369), prime))
        require(fixed_rank == 22 and alternative_rank == 23, "coefficient choice changes actual rank")
        fixtures.append({"rectangle": [r, s], "length": r*s, "fixed_probe_rank_mod_101": fixed_rank,
                         "alternative_rank_over_every_field": 23, "alternative_rank_mod_101": alternative_rank,
                         "explicit_minor_determinant": 1, "minor_rows": chosen_rows,
                         "minor_columns": ordered_columns, "unit_triangular_minor": selected})
    return {"relation": [[*e, c] for e, c in sorted(relation.items())],
            "denominator": [[i, j, c] for (i, j), c in sorted(d.items())],
            "numerator": [[i, j, c] for (i, j), c in sorted(numerator.items())],
            "generic_witness_F": [[1, 0, 1], [0, 1, 1], [4, 3, 1]],
            "fixtures": fixtures, "uniform_rank_separation_proved_in_note": True,
            "scope": "For every r>=8,s>=6 the fixed coefficient rule has a nontrivial relation on 23 ground columns, whereas a unit minor proves that a generic polynomial F has rank 23 there. This establishes specialization of the probe, not monotonicity of fan counts or an asymptotic lower bound."}


def coherent_specialization_controls():
    """A written all-parameter Groebner basis, checked by direct column ranks."""
    fixtures = []
    for a, b in ((4, 3), (5, 3), (4, 4)):
        r, s, prime = 2*a, 2*b, 101
        cells = tuple(product(range(r), range(s)))
        index = {p: i for i, p in enumerate(cells)}
        polynomial = {(1, 0): 1, (0, 1): 1, (a, b): 1}
        powers = [{(0, 0): 1}]
        for _ in range(r+s-2):
            following = {}
            for (i, j), c in powers[-1].items():
                for (k, l), d in polynomial.items():
                    if i+k < r and j+l < s:
                        p = (i+k, j+l)
                        following[p] = (following.get(p, 0) + c*d) % prime
            powers.append({p: c for p, c in following.items() if c})

        def residue(exponent):
            i, j, k = exponent
            column = [0] * len(cells)
            if k < len(powers):
                for (u, v), c in powers[k].items():
                    if i+u < r and j+v < s:
                        column[index[(i+u, j+v)]] = c
            return tuple(column)

        expected = frozenset([(i, j, 0) for i, j in cells if i < a or j < b] +
                             [(i, j, 1) for i in range(a) for j in range(b)])
        leading = [(r, 0, 0), (0, s, 0), (a, b, 0), (a, 0, 1), (0, b, 1), (0, 0, 2)]
        require(set(minimal_excluded_border(sorted(expected))) == set(leading), "six explicit monomial generators")
        relations = [
            {(r, 0, 0): 1}, {(0, s, 0): 1},
            {(a, b, 0): 1, (0, 0, 1): -1, (1, 0, 0): 1, (0, 1, 0): 1},
            {(a, 0, 1): 1, (a+1, 0, 0): -1, (a, 1, 0): -1},
            {(0, b, 1): 1, (1, b, 0): -1, (0, b+1, 0): -1},
            {(0, 0, 2): 1, (1, 0, 1): -2, (0, 1, 1): -2,
             (2, 0, 0): 1, (1, 1, 0): 2, (0, 2, 0): 1}
        ]
        weight = (2, 2, 3)
        for head, relation in zip(leading, relations):
            require(all(sum(c * residue(e)[i] for e, c in relation.items()) % prime == 0
                        for i in range(len(cells))), "explicit relation has zero residue")
            require(all(dot(weight, head) > dot(weight, e) for e in relation if e != head),
                    "uniform strict weight for the six heads")
        ground = tuple(hyperbolic_cross(r*s, 3))
        columns = {e: residue(e) for e in ground}
        actual = frozenset(greedy(ground, columns, weight, prime))
        require(len(expected) == r*s and actual == expected, "direct greedy rank versus proved standard set")
        require((3, 2, 1) in expected, "fixed-probe universal head becomes standard in the other embedding")
        fixtures.append({"half_sides": [a, b], "rectangle": [r, s], "length": r*s,
                         "F_terms": [[*p, c] for p, c in sorted(polynomial.items())],
                         "strict_weight": weight, "basis": sorted(expected),
                         "leading_monomials": leading, "relations": [[[*e, c] for e, c in sorted(f.items())]
                                                                        for f in relations],
                         "direct_greedy_rank_checked": True, "characteristic_for_check": prime,
                         "all_parameter_proof_valid_over_every_field": True})
    return {"fixtures": fixtures,
            "scope": "For r=2a,s=2b,F=x+y+t*x^a*y^b, t nonzero, the note proves one coherent standard set in every field whenever a,b>=1. For a>=4,b>=3 it contains a monomial universally excluded by the fixed rational probe. This distinguishes actual initial ideals, not just arbitrary column subsets; it does not compare total fan counts."}


def uniform_prefix_comparison(fan):
    """Test entire state cones against the cheaper ordinary-cut model."""
    r, s = fan["rectangle"]
    prime, n = fan["characteristic"], r*s
    coefficients = {(i, j): c for i, j, c in fan["polynomial_terms"]}
    ground, columns, _ = rectangle_model(r, s, coefficients, prime)
    active = [e for e in ground if any(columns[e])]
    counts, counterexamples = {}, {}
    for cone in fan["witnesses"]:
        basis = frozenset(map(tuple, cone["basis"]))
        three_generators = all(tuple(int(i == j) for i in range(3)) in basis for j in range(3))
        polygon = [tuple(map(Fraction, p)) for p in cone["cone_polygon"]]
        bad = next(((p, a, b) for p in polygon for a in sorted(basis) for b in active
                    if b not in basis and dot(p, b) < dot(p, a)), None)
        label = ("three_generators" if three_generators else "coordinate_face")
        label += "_has_nonuniform_region" if bad else "_equals_uniform_cut_cone"
        counts[label] = counts.get(label, 0) + 1
        if bad and label not in counterexamples:
            vertex, included, skipped = bad
            raw = cone["strict_weight"]
            interior = tuple(Fraction(x, sum(raw)) for x in raw)
            f = dot(vertex, skipped) - dot(vertex, included)
            g = dot(interior, skipped) - dot(interior, included)
            k = max(1, int(g // (-f)) + 1)
            limit = k + len(active)*(len(active)-1)//2 + 1
            while True:
                weight = integer_weight(tuple((k*x+y)/(k+1) for x, y in zip(vertex, interior)))
                scores = [dot(weight, e) for e in active]
                if len(set(scores)) == len(scores) and dot(weight, skipped) < dot(weight, included):
                    break
                k += 1
                require(k <= limit, "generic interior witness avoids finitely many comparison walls")
            prefix = sorted(active, key=lambda e: dot(weight, e))[:n]
            actual = frozenset(greedy(ground, columns, weight, prime))
            rank = len(greedy(prefix, columns, weight, prime))
            require(actual == basis and rank < n, "coherent basis rejects a dependent nonzero prefix")
            counterexamples[label] = {"strict_weight": weight, "basis": sorted(basis),
                                      "first_n_nonzero_monomials": prefix, "prefix_rank": rank,
                                      "cheaper_skipped_monomial": skipped, "more_expensive_basis_monomial": included}
    return {"rectangle": [r, s], "length": n, "cone_counts": counts,
            "counterexamples": counterexamples,
            "scope": "A nonuniform label means some strict weights in that cone reject the first n nonzero monomials. Equality labels are checked on every polygon vertex and hence on the whole cone. No total ordinary-cut count or asymptotic growth is inferred."}


def coherent_skip_controls():
    """Check the uniform socle-pair obstruction inside actual weight prefixes."""
    fixtures = []
    for q in range(3, 9):
        prime = 101
        cells = tuple(product(range(q), repeat=2))
        index = {p: i for i, p in enumerate(cells)}
        polynomial = {(i, j): (i+1)**2 + 3*(j+1) + 5*i*j for i, j in cells if i or j}
        powers = [{(0, 0): 1}]
        for _ in range(2*q-2):
            following = {}
            for (i, j), c in powers[-1].items():
                for (k, l), d in polynomial.items():
                    if i+k < q and j+l < q:
                        p = (i+k, j+l)
                        following[p] = (following.get(p, 0) + c*d) % prime
            powers.append({p: c for p, c in following.items() if c})

        def residue(exponent):
            _, b, c = exponent
            column = [0]*len(cells)
            for (i, j), value in powers[c].items():
                if j+b < q:
                    column[index[(i, j+b)]] = value
            return tuple(column)

        nonzero = [(0, b, c) for b in range(q) for c in range(2*q-1-b)]
        columns = {e: residue(e) for e in nonzero}
        require(all(any(c) for c in columns.values()), "nonzero yz monomials have the proved exact support range")
        weight = (4*q*q, 2*q-1, 1)
        ordered = sorted(nonzero, key=lambda e: dot(weight, e))
        require(len({dot(weight, e) for e in ordered}) == len(ordered), "separated yz rows")
        prefix = ordered[:q*q]
        u, v = (0, 0, 2*q-2), (0, 1, 2*q-3)
        require(u in prefix and v in prefix, "both socle columns occur in the first n nonzero terms")
        b = polynomial[(0, 1)]
        require(all((x-2*b*y) % prime == 0 for x, y in zip(columns[u], columns[v])),
                "socle-pair dependence in the actual prefix")
        rank = len(greedy(prefix, columns, weight, prime))
        basis = frozenset(greedy(nonzero, columns, weight, prime))
        require(rank < q*q and len(basis) == q*q and max(dot(weight, e) for e in basis) < weight[0],
                "full coherent basis is obtained before any x monomial")
        require(all(tuple(x-(i == j) for i, x in enumerate(e)) in basis
                    for e in basis for j in range(3) if e[j]), "actual basis is a downset")
        fixtures.append({"square_side": q, "length": q*q, "characteristic": prime,
                         "weight": weight, "first_n_nonzero_monomials": prefix,
                         "prefix_rank": rank, "actual_basis": sorted(basis),
                         "parallel_pair": [u, v], "pair_relation_coefficients": [1, -2*b]})
    return {"fixtures": fixtures,
            "uniform_scope": "Every q>=3 and every F with nonzero x and y linear coefficients over characteristic zero, or characteristic greater than 2q-2, has this coherent-prefix obstruction. Finite controls use the fixed coefficient rule modulo 101.",
            "scope": "A growing family of genuine dependent prefixes, not an inference from an arbitrary column circuit. This does not establish growing fan counts; the linear F control has at most five cones."}


def generic_diagonal_graph_controls():
    """Direct determinant polynomials, independent of Gaussian elimination."""
    fixtures = []
    for r, s, blocks, prime in ((2, 2, 3, 0), (2, 3, 2, 7)):
        n = r*s
        chars = tuple(product(range(r), range(s)))
        rows = chars
        if prime:
            zr, zs = 6, 2
            fourier = [[pow(zr, i*a, prime)*pow(zs, j*b, prime) % prime
                        for a, b in chars] for i, j in rows]
        else:
            fourier = [[(-1)**(i*a+j*b) for a, b in chars] for i, j in rows]
        candidates = [(a, b, c) for c in range(blocks) for a, b in chars]
        char_index = {p: i for i, p in enumerate(chars)}
        signed_permutations = [(perm, (-1)**sum(perm[i] > perm[j] for i in range(n) for j in range(i+1, n)))
                               for perm in permutations(range(n))]
        ledger, samples = [], []
        for selected in combinations(candidates, n):
            coefficients = {}
            for perm, sign in signed_permutations:
                exponent = tuple(selected[perm[i]][2] for i in range(n))
                coefficient = sign
                for i in range(n):
                    coefficient *= fourier[i][char_index[selected[perm[i]][:2]]]
                coefficients[exponent] = coefficients.get(exponent, 0) + coefficient
            coefficients = {e: (c % prime if prime else c) for e, c in coefficients.items()
                            if (c % prime if prime else c)}
            require(coefficients, "every chosen maximal determinant is a nonzero generic polynomial")
            monomial = min(coefficients)
            witness = {"columns": selected, "parameter_exponents": monomial,
                       "nonzero_coefficient": coefficients[monomial]}
            ledger.append(witness)
            if len(samples) < 3:
                samples.append(witness)
        require(len(ledger) == comb(len(candidates), n), "exhaust all tested maximal minors")
        zero_small_minors = sum(fourier[i][a]*fourier[j][b] == fourier[i][b]*fourier[j][a]
                                if not prime else
                                (fourier[i][a]*fourier[j][b]-fourier[i][b]*fourier[j][a]) % prime == 0
                                for i, j in combinations(range(n), 2) for a, b in combinations(range(n), 2))
        require(zero_small_minors > 0, "the Fourier matrix need not have all smaller minors nonzero")
        fixtures.append({"group_orders": [r, s], "size": n, "blocks": blocks,
                         "base_characteristic": prime, "fourier_matrix": fourier,
                         "zero_two_by_two_minors": zero_small_minors,
                         "maximal_generic_minors_checked": len(ledger),
                         "coefficient_witness_sha256": hashlib.sha256(json.dumps(ledger, sort_keys=True).encode()).hexdigest(),
                         "sample_coefficient_witnesses": samples})
    return {"fixtures": fixtures,
            "scope": "For finite abelian character tables the note proves the generic diagonal-power full-rank property by transitive row matroids and classical matroid partition. Here determinant polynomials are expanded directly in row parameters, not evaluated at random parameters. No statement about finite-field scalar specialization, conditioning or algorithmic speed is inferred."}


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


def pairwise_realizability_controls():
    """A connected graph of realizable exchanges is not one ideal's fan.

    The written proof is uniform. These small cases independently replay the
    binomial quotients by exact multiplication and modular linear algebra.
    They do not enlarge the nonlinear-rectangle fan search.
    """
    fixtures = []
    for q in range(2, 8):
        shell = tuple((i, q-i, 0) for i in range(q+1))
        lower = frozenset((i, j, 0) for i in range(q) for j in range(q-i))
        k = (q+1)//2
        choices = tuple(frozenset(s) for s in combinations(shell, k))
        index = {s: i for i, s in enumerate(choices)}
        graph = [set() for _ in choices]
        edge_digests = {p: hashlib.sha256() for p in (2, 3, 101)}
        edges = 0
        for i, s in enumerate(choices):
            for a in sorted(s):
                for b in shell:
                    if b in s:
                        continue
                    t = s.difference({a}) | {b}
                    j = index[t]
                    graph[i].add(j)
                    if j <= i:
                        continue
                    edges += 1
                    left, right = lower | s, lower | t
                    union = left | right
                    reference = tuple(sorted(left))
                    n = len(reference)
                    require(n == q*(q+1)//2+k == len(right), "equal quotient lengths")
                    require(all(tuple(x-(axis == d) for d, x in enumerate(m)) in union
                                for m in union for axis in range(3) if m[axis]),
                            "the union is a downset")
                    require(all(tuple(x+(axis == d) for d, x in enumerate(m)) not in union
                                for m in (a, b) for axis in range(3)),
                            "both exchanged monomials are in the union's socle")
                    position = {m: h for h, m in enumerate(reference)}

                    # In K[x,y,z]/(C, x^a-x^b), C the monomial complement
                    # of the union, b has exactly the same residue as a.
                    def residual_position(m):
                        return position.get(a if m == b else m)

                    actions = [[residual_position(tuple(x+(axis == d) for d, x in enumerate(m)))
                                for m in reference] for axis in range(3)]
                    for u, v in combinations(range(3), 2):
                        for h in range(n):
                            vu = None if actions[u][h] is None else actions[v][actions[u][h]]
                            uv = None if actions[v][h] is None else actions[u][actions[v][h]]
                            require(vu == uv, "independent multiplication actions commute")
                    columns = {m: tuple(int(h == residual_position(m)) for h in range(n))
                               for m in union}
                    for p in edge_digests:
                        bases = {frozenset(greedy(tuple(sorted(union)), columns, w, p))
                                 for w in ((1, 2, 3), (2, 1, 3))}
                        require(bases == {left, right}, "both predicted initial bases by elimination")
                        edge_digests[p].update(json.dumps([i, j, sorted(map(sorted, bases))]).encode())

        reached, pending = {0}, [0]
        while pending:
            for j in graph[pending.pop()]:
                if j not in reached:
                    reached.add(j)
                    pending.append(j)
        require(len(reached) == len(choices), "every candidate belongs to one connected exchange graph")
        require(edges == len(choices)*k*(q+1-k)//2, "Johnson graph degree and edge count")
        states = [tuple(sum(a[d] for a in lower | s) for d in range(3)) for s in choices]
        total = sum(sum(a) for a in lower)+k*q
        require(all(x+y == total and z == 0 for x, y, z in states),
                "all candidate state vectors lie on one affine line")
        # Any initial basis is a unique weight-minimum among all bases.
        # The collinearity check is supporting data, not a finite substitute
        # for the uniform at-most-two proof in the note.
        fixtures.append({"shell_degree": q, "shell_points": list(shell), "selected_shell_points": k,
                         "quotient_length": len(lower)+k, "candidate_vertices": len(choices),
                         "distinct_state_vectors": len(set(states)), "realizable_exchange_edges": edges,
                         "exchange_graph_connected": True, "every_exchange_is_a_binomial_quotient": True,
                         "all_state_vectors_collinear": True, "uniform_bound_for_one_ideal": 2,
                         "all_edge_controls_by_characteristic": {str(p): h.hexdigest()
                                                                  for p, h in edge_digests.items()}})

    # Published nonstandard-grading counterexample, not a new family:
    # I_c=(x^7,xy,y^4+c*x^6), c!=0, weights (2,3).  The degree-12
    # quotient columns x^6, x^3*y^2, y^4 are (1,0,-c) on the x^6 row.
    nonstandard = []
    for p in (2, 3, 101):
        for c in (1, p-1):
            exponents = ((6, 0, 0), (3, 2, 0), (0, 4, 0))
            columns = dict(zip(exponents, ((1,), (0,), ((-c) % p,))))
            require(all(2*a+3*b == 12 for a, b, _ in exponents), "published weighted degree")
            require(len(greedy(exponents, columns, (1, 2, 3), p)) == 1 and
                    sum(not any(v) for v in columns.values()) == 1,
                    "a forced loop prevents the uniform rank-one matroid")
            nonstandard.append({"characteristic": p, "coefficient": c,
                                "weighted_degree": 12, "column_rank": 1, "forced_loops": [[3, 2, 0]]})
    # Example 3.26 is stronger: removing loops and coloops still does
    # not give uniformity. Set its parameters (a,b,c)=(1,2,1).
    no_loops = []
    for p in (2, 3, 101):
        exponents = ((9, 0, 0), (6, 2, 0), (3, 4, 0), (0, 6, 0))
        columns = dict(zip(exponents, ((1, 0), (0, 1), (0, -1), (-1, -1))))
        rank = lambda s: len(greedy(s, columns, (1, 2, 3), p))
        require(all(2*a+3*b == 18 for a, b, _ in exponents), "second published weighted degree")
        circuits = [list(s) for k in range(1, 5) for s in combinations(exponents, k)
                    if rank(s) < k and all(rank(s[:i]+s[i+1:]) == k-1 for i in range(k))]
        expected = [[exponents[1], exponents[2]],
                    [exponents[0], exponents[1], exponents[3]],
                    [exponents[0], exponents[2], exponents[3]]]
        require(circuits == expected and all(rank((e,)) == 1 for e in exponents) and
                all(rank(tuple(e for e in exponents if e != deleted)) == 2 for deleted in exponents),
                "published nonuniform matroid has neither loops nor coloops")
        no_loops.append({"characteristic": p, "parameters_a_b_c": [1, 2 % p, 1],
                         "weighted_degree": 18, "column_rank": 2, "circuits": circuits,
                         "no_loops_or_coloops": True})
    return {"shell_exchange_fixtures": fixtures, "nonstandard_grading_source_controls": nonstandard,
            "nonstandard_grading_without_loops_or_coloops": no_loops,
            "uniform_proof_in_note": True, "new_asymptotic_bound": False,
            "scope": "Realizability of each pair is strictly weaker than realizability in one fan. The shell obstruction is an elementary state-polytope control, not a claimed new theorem. The nonstandard-grading example is published by Maclagan--Silversmith."}


def fan_product_controls():
    """Check the all-field planar-to-spatial transfer by two methods.

    Full-ground quotient-matrix selection is compared with Buchberger's
    criterion on the actual lifted polynomial generators. No new search
    family or fitted growth law is used.
    """
    from run_structural_complexity_checks import matvec, univariate_product_mod

    def normal_form(poly, generators, weight, p):
        poly = {a: c % p for a, c in poly.items() if c % p}
        remainder = {}
        leading = [(max(g, key=lambda a: (dot(weight, a), a)), g) for g in generators]
        while poly:
            a = max(poly, key=lambda b: (dot(weight, b), b))
            for b, g in leading:
                if all(x >= y for x, y in zip(a, b)):
                    factor = poly[a]*pow(g[b], -1, p) % p
                    shift = tuple(x-y for x, y in zip(a, b))
                    for c, value in g.items():
                        target = tuple(x+y for x, y in zip(c, shift))
                        value = (poly.get(target, 0)-factor*value) % p
                        if value:
                            poly[target] = value
                        else:
                            poly.pop(target, None)
                    break
            else:
                remainder[a] = poly.pop(a)
        return remainder

    def buchberger(generators, weight, p):
        heads = [max(g, key=lambda a: (dot(weight, a), a)) for g in generators]
        for i, j in combinations(range(len(generators)), 2):
            common = tuple(max(a, b) for a, b in zip(heads[i], heads[j]))
            poly = {}
            for h, sign in ((i, 1), (j, -1)):
                shift = tuple(a-b for a, b in zip(common, heads[h]))
                scale = sign*pow(generators[h][heads[h]], -1, p)
                for a, value in generators[h].items():
                    target = tuple(x+y for x, y in zip(a, shift))
                    poly[target] = (poly.get(target, 0)+scale*value) % p
            require(not normal_form(poly, generators, weight, p), "all polynomial S-pairs reduce to zero")
        return heads, comb(len(generators), 2)

    fixtures = []
    for p in (2, 3, 101):
        for kind, modulus, image in (
                ("local_cyclic", [0, 0, 0, 0, 1], [0, 1, 2, 3]),
                ("possibly_repeated_points", [-6, 11, -6, 1], [0, 1, 1]),
                ("noncyclic_shared_coordinate", None, None),
                ("nonseparating_reduced_points", [24, -50, 35, -10, 1], None)):
            if kind == "nonseparating_reduced_points":
                if p != 101:
                    continue
                m = 4
                branch_points = [(1, 4), (1, 5), (2, 4), (3, 6)]
                reference = [(0, 0), (1, 0), (2, 0), (0, 1)]
                evaluation = lambda a: tuple(prod(pow(c, e, p) for c, e in zip(point, a)) % p
                                             for point in branch_points)
                inv = inverse_mod([[evaluation(a)[i] for a in reference] for i in range(m)], p)
                require(inv is not None, "nonseparating point reference basis")
                matrices = [tuple(tuple(sum(x*y for x, y in zip(row, evaluation(
                    tuple(x+(axis == d) for d, x in enumerate(a))))) % p for row in inv)
                    for a in reference) for axis in range(2)]
            elif modulus is not None:
                m = len(modulus)-1
                coordinates = [tuple(int(i == j) for i in range(m)) for j in range(m)]
                matrices = [tuple(tuple(c % p for c in univariate_product_mod(g, column, modulus))
                                  for column in coordinates) for g in ([0, 1], image)]
            else:
                m, modulus = 4, [0, 0, 0, 0, 1]
                # J=(x^2,y^2-x-xy), reference basis 1,x,y,xy.
                matrices = [((0, 1, 0, 0), (0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 0, 0)),
                            ((0, 0, 1, 0), (0, 0, 0, 1), (0, 1, 0, 1), (0, 0, 0, 0))]
            apply = lambda matrix, vector: tuple(c % p for c in matvec(matrix, vector))
            unit = (1,)+(0,)*(m-1)

            def branch_residue(a):
                vector = unit
                for matrix, power in zip(matrices, a):
                    for _ in range(power):
                        vector = apply(matrix, vector)
                return vector

            require(all(apply(matrices[0], column) == apply(matrices[1], matrices[0][j])
                        for j, column in enumerate(matrices[1])), "planar multiplication commutes")
            for j in range(m):
                v = tuple(int(h == j) for h in range(m))
                total = [0]*m
                for c in modulus:
                    total = [(a+c*b) % p for a, b in zip(total, v)]
                    v = apply(matrices[0], v)
                require(not any(total), "the common monic annihilator has degree m")
            ground = tuple(hyperbolic_cross(m, 2))
            comparisons = tuple(hyperbolic_cross(m+1, 2))
            columns = {a: branch_residue(a) for a in ground}
            cuts = sorted({Fraction(-(a[1]-b[1]), a[0]-b[0])
                           for a, b in combinations(comparisons, 2)
                           if (a[0]-b[0])*(a[1]-b[1]) < 0})
            samples = [(a+b)/2 for a, b in zip([Fraction(0)]+cuts, cuts)]
            samples.append(cuts[-1]+1 if cuts else Fraction(1))
            states = {}
            translation_checks = 0
            for t in samples:
                weight = integer_weight((t, Fraction(1)))
                basis = tuple(sorted(greedy(ground, columns, weight, p)))
                require(len(basis) == m, "full rank on the planar universal ground")
                if kind == "nonseparating_reduced_points":
                    original_columns = {a: tuple(pow(x, a[0], p)*pow(y-4, a[1], p) % p
                                                for x, y in branch_points) for a in ground}
                    require(greedy(ground, original_columns, weight, p) == set(basis),
                            "translation preserves every planar initial ideal")
                    translation_checks += 1
                states.setdefault(basis, weight)
            plane_data, plane_pairs = [], 0
            for basis, weight in sorted(states.items()):
                inv = inverse_mod([[branch_residue(a)[i] for a in basis] for i in range(m)], p)
                require(inv is not None, "planar basis invertible")
                gb = []
                for a in minimal_excluded_border(basis):
                    column = branch_residue(a)
                    coefficients = [sum(x*y for x, y in zip(row, column)) % p for row in inv]
                    poly = {a: 1}
                    poly.update({b: -c % p for b, c in zip(basis, coefficients) if c})
                    require(all(dot(weight, a) > dot(weight, b) for b in poly if b != a),
                            "strict planar reduced-basis comparisons")
                    gb.append(poly)
                _, pairs = buchberger(gb, weight, p)
                plane_pairs += pairs
                plane_data.append((basis, weight, gb))

            n = 3*m
            lifted_ground = tuple(hyperbolic_cross(n, 3))
            pure_columns = [branch_residue((a, 0)) for a in range(m)]
            spatial_matrices = []
            for axis in range(3):
                matrix = []
                for block in range(3):
                    for j in range(m):
                        vector = [0]*n
                        if axis == 0:
                            if block == 0:
                                if j < m-1:
                                    vector[j+1] = 1
                                else:
                                    vector[:m] = [-c % p for c in modulus[:-1]]
                            else:
                                vector[block*m:(block+1)*m] = matrices[0][j]
                        elif block == 0:
                            vector[axis*m:(axis+1)*m] = pure_columns[j]
                        elif block == axis:
                            vector[block*m:(block+1)*m] = matrices[1][j]
                        matrix.append(tuple(vector))
                spatial_matrices.append(tuple(matrix))
            for u, v in combinations(range(3), 2):
                for j in range(n):
                    require(apply(spatial_matrices[u], spatial_matrices[v][j]) ==
                            apply(spatial_matrices[v], spatial_matrices[u][j]),
                            "spatial multiplication commutes on every reference vector")
            spatial_columns = {}
            for a in lifted_ground:
                vector = (1,)+(0,)*(n-1)
                for matrix, power in zip(spatial_matrices, a):
                    for _ in range(power):
                        vector = apply(matrix, vector)
                spatial_columns[a] = vector
            witnesses, spatial_pairs, point_checks = [], 0, 0
            for left, right in product(plane_data, repeat=2):
                d1, (x1, y1), g1 = left
                d2, (x2, y2), g2 = right
                weight = (x1*x2, y1*x2, y2*x1)
                expected = ({(a, 0, 0) for a in range(m)} |
                            {(a, b+1, 0) for a, b in d1} | {(a, 0, b+1) for a, b in d2})
                actual = greedy(lifted_ground, spatial_columns, weight, p)
                require(actual == expected and len(actual) == n, "all products of branch bases by direct elimination")
                if kind in ("possibly_repeated_points", "nonseparating_reduced_points") and p == 101:
                    branch_points = ([(t, t+t*t) for t in (1, 2, 3)] if m == 3 else
                                     [(1, 4), (1, 5), (2, 4), (3, 6)])
                    points = [(t, 0, 0) for t in range(1, m+1)]
                    points += [(x, y, 0) for x, y in branch_points]
                    points += [(x, 0, y) for x, y in branch_points]
                    point_columns = {a: tuple(prod(pow(c, e, p) for c, e in zip(point, a)) % p
                                              for point in points) for a in lifted_ground}
                    require(len(set(points)) == n and
                            greedy(lifted_ground, point_columns, weight, p) == expected,
                            "independent reduced point-evaluation replay")
                    point_checks += 1
                gb = [{(a, 0, 0): c % p for a, c in enumerate(modulus) if c % p}, {(0, 1, 1): 1}]
                gb += [{(a, b+1, 0): c for (a, b), c in g.items()} for g in g1]
                gb += [{(a, 0, b+1): c for (a, b), c in g.items()} for g in g2]
                heads, pairs = buchberger(gb, weight, p)
                spatial_pairs += pairs
                standard = {a for a in lifted_ground if not any(all(x >= y for x, y in zip(a, b)) for b in heads)}
                require(standard == expected, "Buchberger initial ideal versus quotient-matrix prediction")
                witnesses.append({"branch_bases": [d1, d2], "strict_weight": weight,
                                  "spatial_basis": sorted(actual)})
            require(len({tuple(w["spatial_basis"]) for w in witnesses}) == len(states)**2,
                    "different branch pairs give different spatial initial ideals")
            fixtures.append({"kind": kind, "characteristic": p, "plane_length": m,
                             "common_annihilator_low_to_high": modulus, "spatial_length": n,
                             "plane_initial_ideals": len(states), "spatial_initial_ideals": len(states)**2,
                             "plane_comparison_intervals_exhausted": len(samples),
                             "translation_invariance_checks": translation_checks,
                             "plane_s_pairs": plane_pairs, "spatial_s_pairs": spatial_pairs,
                             "independent_reduced_point_evaluation_checks": point_checks,
                             "all_spatial_multiplication_commutators_checked": True, "witnesses": witnesses})
    return {"fixtures": fixtures, "uniform_proof_in_note": True, "new_asymptotic_separation": False,
            "scope": "For every planar finite quotient the written construction squares its fan count at length 3m. Exact checks cover cyclic, noncyclic, local and possibly nonreduced point presentations. No planar superlinear-power lower family or novelty certificate is supplied."}


def plane_curve_carrier_controls():
    """Exact wall-by-wall checks of the irreducible-curve count proof."""
    from run_structural_complexity_checks import univariate_product_mod, univariate_power_mod

    def audit(n, p, curve, residue, label):
        degree = max(sum(a) for a in curve)
        ground, comparisons = tuple(hyperbolic_cross(n, 2)), tuple(hyperbolic_cross(n+1, 2))
        columns = {a: residue(a) for a in ground}
        require(all(sum(c*residue(a)[i] for a, c in curve.items()) % p == 0 for i in range(n)),
                "the claimed carrier belongs to the quotient ideal")
        cuts = sorted({Fraction(-(a[1]-b[1]), a[0]-b[0]) for a, b in combinations(comparisons, 2)
                       if (a[0]-b[0])*(a[1]-b[1]) < 0})
        samples = [(a+b)/2 for a, b in zip([Fraction(0)]+cuts, cuts)]
        samples.append(cuts[-1]+1 if cuts else Fraction(1))
        states, inverses = {}, {}

        def relation(basis, a):
            if basis not in inverses:
                inverses[basis] = inverse_mod([[residue(b)[i] for b in basis] for i in range(n)], p)
            inv = inverses[basis]
            require(inv is not None, "every planar standard set is a quotient basis")
            coefficients = [sum(x*y for x, y in zip(row, residue(a))) % p for row in inv]
            poly = {a: 1}
            poly.update({b: -c % p for b, c in zip(basis, coefficients) if c})
            require(all(sum(c*residue(b)[i] for b, c in poly.items()) % p == 0 for i in range(n)),
                    "the border or wall polynomial is an actual quotient relation")
            return poly

        for t in samples:
            weight = (t.numerator, t.denominator)
            basis = tuple(sorted(greedy(ground, columns, weight, p)))
            require(len(basis) == n, "full planar quotient rank")
            if basis in states:
                continue
            lower, upper = Fraction(0), None
            for a in minimal_excluded_border(basis):
                poly = relation(basis, a)
                for b in poly:
                    if b == a:
                        continue
                    dx, dy = a[0]-b[0], a[1]-b[1]
                    if dx > 0:
                        lower = max(lower, Fraction(-dy, dx))
                    elif dx < 0:
                        bound = Fraction(-dy, dx)
                        upper = min(upper, bound) if upper is not None else bound
                    else:
                        require(dy > 0, "positive constant border comparison")
            require(lower < t and (upper is None or t < upper), "strict complete cone witness")
            moment = tuple(sum(a[i] for a in basis) for i in range(2))
            require(max(moment) <= n*(n-1)//2, "universal downset moment budget")
            states[basis] = {"weight": weight, "interval": (lower, upper), "moment": moment}
        ordered = sorted(states, key=lambda b: states[b]["interval"][0])
        require(states[ordered[0]]["interval"][0] == 0 and
                states[ordered[-1]]["interval"][1] is None, "both endpoint sectors covered")
        walls, budget = [], 0
        for left, right in zip(ordered, ordered[1:]):
            t = states[left]["interval"][1]
            require(t == states[right]["interval"][0], "complete fan intervals meet without gaps")
            a, b = t.denominator, t.numerator
            before, after = states[left]["moment"], states[right]["moment"]
            k, remainder = divmod(before[0]-after[0], a)
            require(k > 0 and remainder == 0 and after[1]-before[1] == k*b,
                    "the integral state jump has the claimed primitive direction")
            require(max(a, b) <= n-1, "a wall direction comes from standard monomials")
            budget += k*(a+b)
            scores = {v: v[0]*b+v[1]*a for v in curve}
            heads = [v for v in curve if scores[v] == max(scores.values())]
            row = {"ratio_wx_wy": str(t), "primitive_direction": [a, -b], "state_jump_multiple": k}
            if len(heads) > 1:
                row["kind"] = "carrier_boundary"
            elif min(a, b) < degree:
                row["kind"] = "small_direction"
            else:
                r, s = heads[0]
                require(r > 0 and s > 0, "a large direction cannot occur in a pure-head strip")
                u = next(v for v in left if v not in right)
                poly = relation(right, u)
                height = u[0]*b+u[1]*a
                require(all(v[0]*b+v[1]*a <= height for v in poly), "all lift terms lie below the wall height")
                tied = [v for v in poly if v != u and v[0]*b+v[1]*a == height]
                require(len(tied) == 1, "a large-direction fiber has only two carrier-standard monomials")
                v = tied[0]
                ell, rem = divmod(u[0]-v[0], a)
                require(rem == 0 and 1 <= ell <= k and v[1]-u[1] == ell*b,
                        "the exhibited exchange is part of the total monotone jump")
                require(v[0] < r and u[1] < s and all(not (x >= r and y >= s) for x, y in poly),
                        "the nonzero wall lift is reduced modulo the carrier")
                poly_degree = max(sum(w) for w in poly)
                require(n <= degree*poly_degree <= 3*degree*ell*max(a, b),
                        "Bezout length and the uniform wall-degree bound")
                row.update({"kind": "large_direction", "carrier_head": [r, s],
                            "exchange_endpoints": [u, v], "exchange_multiple": ell,
                            "relation_degree": poly_degree,
                            "polynomial_relation_mod_p": [[x, y, c] for (x, y), c in sorted(poly.items())]})
            walls.append(row)
        require(budget <= n*(n-1), "all jumps share one global moment-variation budget")
        kinds = {kind: sum(w["kind"] == kind for w in walls)
                 for kind in ("carrier_boundary", "small_direction", "large_direction")}
        require(kinds["carrier_boundary"] <= degree and
                kinds["small_direction"] <= 2*(degree-1)*(n-1) and
                kinds["large_direction"] <= 3*degree*(n-1) and len(states) <= 5*degree*n,
                "the three uniform count estimates apply simultaneously")
        return {"model": label, "characteristic": p, "quotient_length": n,
                "carrier": [[i, j, c] for (i, j), c in sorted(curve.items())], "carrier_degree": degree,
                "initial_ideal_count": len(states), "wall_kinds": kinds, "total_moment_variation": budget,
                "bound_5_D_n": 5*degree*n, "complete_positive_fan": True,
                "witnesses": [{"row_lengths": [sum(a[1] == j for a in basis)
                                                for j in range(max(a[1] for a in basis)+1)],
                               "strict_weight": states[basis]["weight"],
                               "interval_wx_wy": [str(x) if x is not None else None
                                                  for x in states[basis]["interval"]]}
                              for basis in ordered], "walls": walls}

    fixtures = []
    for p in (5, 7, 101):
        for n in (7, 12, 18):
            modulus = [0]*n+[1]
            mul = lambda a, b: tuple(v % p for v in univariate_product_mod(a, b, modulus))
            x0, y0 = next((x, y) for x in range(p) for y in range(p)
                          if (x*x*y+x*y*y+1) % p == 0 and (x*x+2*x*y) % p)
            x, y = [x0, 1], [y0]+[0]*(n-1)
            evaluate = lambda: tuple((a+b+(i == 0)) % p for i, (a, b) in
                                    enumerate(zip(mul(mul(x, x), y), mul(x, mul(y, y)))))
            for j in range(1, n):
                y[j] = -evaluate()[j]*pow((x0*x0+2*x0*y0) % p, -1, p) % p
            require(not any(evaluate()), "exact formal implicit-function lift on the cubic")
            residue = lambda a: mul(univariate_power_mod(x, a[0], modulus),
                                    univariate_power_mod(y, a[1], modulus))
            f = audit(n, p, {(2, 1): 1, (1, 2): 1, (0, 0): 1}, residue, "smooth_cubic_local_jet")
            f["local_parameter_images"] = {"x": x, "y": y}
            fixtures.append(f)
    for p in (2, 3, 101):
        for n in (7, 12):
            modulus = [0]*n+[1]
            x, y = [1, 1], [(-1)**j % p for j in range(n)]
            mul = lambda a, b: tuple(v % p for v in univariate_product_mod(a, b, modulus))
            residue = lambda a: mul(univariate_power_mod(x, a[0], modulus),
                                    univariate_power_mod(y, a[1], modulus))
            f = audit(n, p, {(1, 1): 1, (0, 0): -1}, residue, "hyperbola_local_jet")
            require(f["initial_ideal_count"] == n, "sharp linear lower family in every characteristic")
            fixtures.append(f)
    for p in (2, 101):
        n, modulus = 7, [0]*7+[1]
        residue = lambda a: tuple(v % p for v in univariate_power_mod([1, 1], sum(a), modulus))
        f = audit(n, p, {(1, 0): 1, (0, 1): -1}, residue, "line_local_jet")
        require(f["initial_ideal_count"] == 2, "degree-one carrier has at most two initials")
        fixtures.append(f)
    n, p, q, root = 7, 23, 11, 2
    require(pow(root, q, p) == 1 and root != 1, "a primitive eleventh root in the verification field")
    points = [(pow(root, j, p), pow(root, 7*j, p)) for j in range(n)]
    residue = lambda a: tuple(pow(x, a[0], p)*pow(y, a[1], p) % p for x, y in points)
    f = audit(n, p, {(0, 2): 1, (3, 0): -1}, residue, "seven_points_on_cyclic_orbit")
    require(f["initial_ideal_count"] == 9, "generic-maximality hypothesis fails already in the plane")
    ground = tuple(hyperbolic_cross(n, 2))
    cuts = sorted({Fraction(-(a[1]-b[1]), a[0]-b[0]) for a, b in combinations(ground, 2)
                   if (a[0]-b[0])*(a[1]-b[1]) < 0})
    samples = [(a+b)/2 for a, b in zip([Fraction(0)]+cuts, cuts)]+[cuts[-1]+1]
    generic = {tuple(sorted(sorted(ground, key=lambda a: (a[0]*t+a[1], a))[:n])) for t in samples}
    require(len(generic) == 8, "independent generic corner-cut comparison")
    for witness in f["witnesses"]:
        t = Fraction(*witness["strict_weight"])
        seen, basis = set(), []
        for a in sorted(ground, key=lambda a: (a[0]*t+a[1], a)):
            color = (a[0]+7*a[1]) % q
            if color not in seen:
                seen.add(color)
                basis.append(a)
            if len(basis) == n:
                break
        expected = [(i, j) for j, length in enumerate(witness["row_lengths"]) for i in range(length)]
        require(sorted(basis) == sorted(expected), "cyclic color rule versus independent point elimination")
    f.update({"points": points, "orbit_size": q, "generic_count": len(generic),
              "positive_weight_fan_boundaries": [w["ratio_wx_wy"] for w in f["walls"]]})
    fixtures.append(f)
    require(sum(x["wall_kinds"]["large_direction"] for x in fixtures) > 0,
            "the Bezout charging step is exercised, not only pure-head strip cases")
    return {"fixtures": fixtures,
            "basis_encoding": "row_lengths[j]=lambda means x^i*y^j, 0<=i<lambda, are standard; this is an exact compact partition witness",
            "uniform_bound_proved_in_note": "g(I)<=5*D*n for an irreducible carrier",
            "sharp_fixed_degree_order": "Theta_D(n) for D>=2, taking the maximum over carriers of degree at most D",
            "new_unrestricted_asymptotic": False, "novelty_established": False,
            "scope": "Exact complete fans and wall lifts test the written proof. Smooth cubic jets include characteristics smaller than quotient length. A finite cyclic count gain is not an asymptotic lower family."}


def finite_grading_controls():
    """Complete planar fans for the classical relative-lattice bound."""
    fixtures = []

    def audit(n, p, residue, generators, label, q=None, h=None, degree=None):
        ground = tuple(hyperbolic_cross(n, 2))
        comparisons = tuple(hyperbolic_cross(n+1, 2))
        columns = {a: residue(a) for a in ground}
        differences = []
        for polynomial in generators:
            require(all(sum(c*residue(a)[i] for a, c in polynomial.items()) % p == 0
                        for i in range(n)), "a grading generator is an exact quotient relation")
            head = next(iter(polynomial))
            differences.extend(tuple(a[i]-head[i] for i in range(2)) for a in polynomial if a != head)
        index = 0
        for a, b in combinations(differences, 2):
            index = gcd(index, abs(a[0]*b[1]-a[1]*b[0]))
        lattice_rank = 2 if index else int(any(any(a) for a in differences))
        if h is not None:
            require(index == h, "support differences of the displayed reduced lex basis give its exact index")
        cuts = sorted({Fraction(-(a[1]-b[1]), a[0]-b[0]) for a, b in combinations(comparisons, 2)
                       if (a[0]-b[0])*(a[1]-b[1]) < 0})
        samples = [(a+b)/2 for a, b in zip([Fraction(0)]+cuts, cuts)]
        samples.append(cuts[-1]+1 if cuts else Fraction(1))
        states = {}
        target_histogram = [sum(k % h == c for k in range(n)) for c in range(h)] if h else None
        for t in samples:
            weight = (t.numerator, t.denominator)
            basis = tuple(sorted(greedy(ground, columns, weight, p)))
            require(len(basis) == n, "graded complete fan has full quotient rank")
            if basis in states:
                continue
            inv = inverse_mod([[residue(b)[i] for b in basis] for i in range(n)], p)
            require(inv is not None, "independent inverse verifies each graded basis")
            lower, upper = Fraction(0), None
            for a in minimal_excluded_border(basis):
                coefficients = [sum(x*y for x, y in zip(row, residue(a))) % p for row in inv]
                for b, coefficient in zip(basis, coefficients):
                    if not coefficient:
                        continue
                    dx, dy = a[0]-b[0], a[1]-b[1]
                    if dx > 0:
                        lower = max(lower, Fraction(-dy, dx))
                    elif dx < 0:
                        bound = Fraction(-dy, dx)
                        upper = min(upper, bound) if upper is not None else bound
                    else:
                        require(dy > 0, "valid vertical graded border comparison")
            require(lower < t and (upper is None or t < upper), "strict complete graded fan cone")
            moment = tuple(sum(a[i] for a in basis) for i in range(2))
            if h:
                histogram = [sum((a+q*b) % h == c for a, b in basis) for c in range(h)]
                require(histogram == target_histogram, "basis has the quotient's fixed dimensions in every finite grade")
                for (a, b), column in columns.items():
                    require(all(not value or (k-a-q*b) % h == 0 for k, value in enumerate(column)),
                            "each residue column lies entirely in its claimed group degree")
            if degree is not None:
                if t < Fraction(1, degree):
                    require(basis == tuple((i, 0) for i in range(n)), "the unique pure-x sector")
                else:
                    require(all(a < degree for a, b in basis), "all remaining bases lie in the degree-D strip")
                    require(moment[0] <= n*(degree-1)//2 and moment[1] <= n*(n-1)//2,
                            "the two strip moment bounds hold simultaneously")
            states[basis] = {"moment": moment, "strict_weight": weight, "interval": (lower, upper)}
        ordered = sorted(states, key=lambda b: states[b]["interval"][0])
        require(states[ordered[0]]["interval"][0] == 0 and states[ordered[-1]]["interval"][1] is None,
                "both graded fan endpoint sectors covered")
        require(all(states[a]["interval"][1] == states[b]["interval"][0]
                    for a, b in zip(ordered, ordered[1:])), "graded fan coverage has no gaps")
        moments = [states[b]["moment"] for b in ordered]
        anchor = moments[0]
        if h:
            require(all((a+q*b-anchor[0]-q*anchor[1]) % h == 0 for a, b in moments),
                    "every moment lies in the same index-h affine lattice")
        if lattice_rank < 2:
            require(len(states) <= lattice_rank+1, "free grading lowers the state dimension")
            if q is not None:
                require(len({a+q*b for a, b in moments}) == 1, "the free grade puts all moments on a line")
        polygon = hull2(moments)
        require(len(polygon) == len(states), "every distinct state moment is an exposed polygon vertex")
        twice_area = abs(sum(a[0]*b[1]-a[1]*b[0] for a, b in zip(polygon, polygon[1:]+polygon[:1])))
        if h:
            require(twice_area % h == 0, "relative twice-area is integral in the grading lattice")
        strip_twice_area = None
        if degree is not None:
            strip_polygon = hull2([states[b]["moment"] for b in ordered if all(a < degree for a, j in b)])
            strip_twice_area = abs(sum(a[0]*b[1]-a[1]*b[0]
                                      for a, b in zip(strip_polygon, strip_polygon[1:]+strip_polygon[:1])))
            require(2*strip_twice_area <= n*n*(degree-1)*(n-1), "strip polygon obeys the shared box area bound")
            if h:
                require(strip_twice_area % h == 0, "strip area uses the same grading lattice")
        return {"model": label, "characteristic": p, "quotient_length": n,
                "generators": [[[a, b, c] for (a, b), c in sorted(f.items())] for f in generators],
                "universal_grading_lattice_rank": lattice_rank,
                "finite_lattice_index": index if lattice_rank == 2 else None,
                "grading_degrees_x_y": [1, q] if q is not None else None,
                "graded_quotient_dimensions": target_histogram,
                "initial_ideal_count": len(states), "twice_state_polygon_area": twice_area,
                "twice_relative_area": twice_area//h if h else None,
                "twice_strip_polygon_area": strip_twice_area,
                "complete_positive_fan": True, "comparison_intervals": len(samples),
                "witnesses": [{"row_lengths": [sum(a[1] == j for a in basis)
                                                for j in range(max(a[1] for a in basis)+1)],
                               "moment": states[basis]["moment"],
                               "strict_weight": states[basis]["strict_weight"],
                               "interval_wx_wy": [str(t) if t is not None else None
                                                  for t in states[basis]["interval"]]}
                              for basis in ordered]}

    def hull2(points):
        points = sorted(set(points))
        def turn(a, b, c):
            return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
        chains = []
        for sequence in (points, list(reversed(points))):
            chain = []
            for a in sequence:
                while len(chain) >= 2 and turn(chain[-2], chain[-1], a) <= 0:
                    chain.pop()
                chain.append(a)
            chains.extend(chain[:-1])
        return chains if len(points) > 1 else points

    models = ((7, (1, 4)), (12, (1, 5)), (12, (1, 3, 5)),
              (12, (1, 3, 4)), (18, (1, 4, 7)), (18, (1, 4, 5)))
    for p in (2, 3, 101):
        for n, support in models+((7, (2,)),):
            powers = [[1]+[0]*(n-1)]
            for j in range(n-1):
                powers.append([sum(powers[-1][k-d] for d in support if k >= d) % p for k in range(n)])
            residue = lambda a: tuple([0]*min(a[0], n)+powers[a[1]][:max(n-a[0], 0)]) if a[1] < n else (0,)*n
            q, h = support[0], 0
            for d in support:
                h = gcd(h, d-q)
            f = {(0, 1): 1, **{(d, 0): -1 for d in support}}
            fixture = audit(n, p, residue, [{(n, 0): 1}, f], "truncated_polynomial_graph",
                            q=q, h=h, degree=max(support))
            fixture["polynomial_support"] = list(support)
            fixtures.append(fixture)
        n, cells = 6, tuple(product(range(3), range(2)))
        residue = lambda a: tuple(int(a == b) for b in cells)
        fixtures.append(audit(n, p, residue, [{(3, 0): 1}, {(0, 2): 1}], "monomial_control"))
    return {"fixtures": fixtures, "uses_classical_lattice_vertex_bound": True,
            "uniform_graph_bound_in_note": "g((x^n,y-G(x))) <= 3+C*n*((D-1)/h)^(1/3) for nonmonomial G of degree D<n and exact grading index h",
            "one_offset_consequence": "g((x^n,y-x-c*x^D))=O(n) uniformly in 2<=D<n, nonzero c and the field",
            "new_unrestricted_asymptotic": False, "novelty_established": False,
            "scope": "Exact complete fans check intrinsic lattice index, homogeneous residue columns, fixed dimensions in every finite grade, the affine moment coset and strip budgets. No fitted growth rate or novel lattice theorem is claimed."}


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
    comparisons = [uniform_prefix_comparison(fan) for fan in rectangles]
    require(any(c["cone_counts"].get("three_generators_has_nonuniform_region", 0) > 0
                for c in comparisons),
            "dependent-prefix regions are not confined to coordinate faces")
    report = {"schema_version": 7, "status": "passed",
              "command": "python3 -B scripts/run_extremal_initial_ideals_checks.py",
              "grid_source_controls": grid, "rectangular_complete_intersections": rectangles,
              "homogeneous_controls": homogeneous,
              "small_characteristic_control": small_characteristic,
              "rational_probe_specialization": rational_probe_controls(),
              "coherent_specialization": coherent_specialization_controls(),
              "ordinary_cut_comparisons": comparisons,
              "coherent_skip_obstruction": coherent_skip_controls(),
              "generic_diagonal_graph": generic_diagonal_graph_controls(),
              "support_matching_obstruction": support_matching_obstruction_controls(),
              "pairwise_realizability_boundary": pairwise_realizability_controls(),
              "planar_fan_product": fan_product_controls(),
              "plane_curve_carriers": plane_curve_carrier_controls(),
              "finite_gradings": finite_grading_controls(),
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
