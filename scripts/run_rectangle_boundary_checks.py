#!/usr/bin/env python3
"""Exact controls for the regional rectangle prefix rule and cubic walls.

The uniform proofs are in docs/rectangle_boundary_fan.md. These finite
checks establish neither the full fan count nor novelty of the results.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import comb, gcd, prod
from pathlib import Path
import hashlib
import json

from run_bilinear_rectangle_checks import rank_sequence
from run_initial_cone_height_checks import add, leading, normal_form, subtract_shift
from run_interleaving_research import minimal_excluded_border
from run_nested_jet_checks import inverse
from run_rectangle_rank_checks import rectangle
from run_structural_complexity_checks import area2, clip, cross, dot, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def coefficients(s, d):
    return {(0, 1): 1, **{(1, j): comb(d, j) for j in range(min(d, s-1)+1)}}


def weight_list(r, s):
    endpoints = sorted({Q(0)} | {Q(sign*k, h) for sign in (-1, 1)
                       for k in range(1, s) for h in range(1, r)})
    samples = [Q(1-2*s, 2)]+[Q(a.numerator+b.numerator, a.denominator+b.denominator)
               for a, b in zip(endpoints, endpoints[1:])]+[Q(2*s-1, 2)]
    weights = [(s*x.denominator+max(x.numerator, 0), x.denominator,
                s*x.denominator-min(x.numerator, 0)) for x in samples]
    require(len(weights) <= 2*(r-1)*(s-1)+2 and all(max(w) <= 2*r*s for w in weights),
            "linear common-list cardinality and integer height")
    return samples, weights


def slots(t, s, d):
    return [(0, j) for j in range(s)]+[(i, j) for i in range(1, t+1)
                                      for j in range(min(d, s))]


def predicted_basis(r, s, d, rho):
    if rho < -d:
        return frozenset((i, j, 0) for i in range(r) for j in range(s))
    basis = set()
    for t in range(r):
        choices = sorted(slots(t, s, d), key=lambda a: (rho*a[0]+a[1], a))
        if len(choices) > s:
            require(rho*choices[s-1][0]+choices[s-1][1] < rho*choices[s][0]+choices[s][1],
                    "common-list parameter has a strict cutoff")
        basis.update((i, j, t-i) for i, j in choices[:s])
    require(len(basis) == r*s and all(
        tuple(x-(i == j) for i, x in enumerate(a)) in basis for a in basis for j in range(3) if a[j]),
        "slot selections form a full monomial downset")
    return frozenset(basis)


def direct_greedy(model, weight, characteristic):
    ground = sorted(hyperbolic_cross(len(model['cells']), 3), key=lambda a: (dot(a, weight), a))
    ranks = rank_sequence([model['residue'](a) for a in ground], characteristic)
    selected = frozenset(a for i, a in enumerate(ground) if ranks[i] > (ranks[i-1] if i else 0))
    require(len(selected) == len(model['cells']), "direct quotient selection has full rank")
    return selected


def block_controls(model, r, s, d, rho, characteristic):
    index = {cell: i for i, cell in enumerate(model['cells'])}
    prefixes = 0
    for t in range(r):
        choices = sorted(slots(t, s, d), key=lambda a: (rho*a[0]+a[1], a))
        leading_columns = [tuple(model['residue']((i, j, t-i))[index[t, k]] for k in range(s))
                           for i, j in choices]
        require(rank_sequence(leading_columns, characteristic) ==
                tuple(min(s, k+1) for k in range(len(choices))),
                "direct leading-component ranks agree with the uniform prefix rule")
        # Check the disjoint-exponent mechanism independently of row elimination.
        present = {}
        for i, j in choices:
            present.setdefault(i, []).append(j)
            require(present[i] == list(range(j+1)), "each prefix has consecutive j values")
            nodes = [-d*a+k for a, values in present.items() for k in range(len(values))]
            require(len(nodes) == len(set(nodes)), "transformed exponent intervals are disjoint")
            if characteristic:
                require(len(nodes) == len({x % characteristic for x in nodes}),
                        "qualifying positive characteristic retains distinct exponents")
            prefixes += 1
    return prefixes


def region_polygon(s):
    poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
    for z in ((1, -(s-1), 0), (0, -(s-1), 1)):
        poly = clip(poly, z, 1)
    return poly


def rational_cone(model, basis, weight, s, audit=True):
    basis = tuple(sorted(basis))
    n = len(basis)
    matrix = [[model['residue'](a)[i] for a in basis] for i in range(n)]
    inv = inverse(matrix)
    polys, normals = [], set()
    for head in minimal_excluded_border(basis):
        col = model['residue'](head)
        coefficients = [sum(x*y for x, y in zip(row, col)) for row in inv]
        require(tuple(sum(c*model['residue'](a)[i] for c, a in zip(coefficients, basis))
                      for i in range(n)) == col, "border coordinates replay over Q")
        polynomial = {head: Q(1)}
        for a, c in zip(basis, coefficients):
            if c:
                polynomial[a] = -c
                z = tuple(x-y for x, y in zip(head, a))
                g = gcd(*z)
                normals.add(tuple(x//g for x in z))
        polys.append(polynomial)
    full = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
    for z in sorted(normals):
        require(dot(weight, z) > 0, "strict rational reduced leading terms")
        full = clip(full, z, 1)
    require(full, "nonempty full global cone")
    regional = full
    for z in ((1, -(s-1), 0), (0, -(s-1), 1)):
        regional = clip(regional, z, 1)
    require(regional, "the global cone meets the prescribed open region")
    pairs = forms = 0
    if audit:
        for g, h in combinations(polys, 2):
            a, b = leading(g, weight), leading(h, weight)
            common = tuple(max(x, y) for x, y in zip(a, b))
            first = {add(c, tuple(x-y for x, y in zip(common, a))): value for c, value in g.items()}
            spair = subtract_shift(first, h, tuple(x-y for x, y in zip(common, b)), Q(1), 0)
            require(not normal_form(spair, polys, weight, 0), "ordinary rational S-pair")
            pairs += 1
        # Border polynomials vanish in the length-n quotient; Buchberger
        # and its n-element standard set therefore certify the same ideal.
        # Full-ground normal-form replay is extra verification for small
        # fixtures. It is not needed again for the larger wall certificates.
        probes = hyperbolic_cross(n, 3) if n <= 16 else ((1, 0, 1), (0, 1, 1), (0, 0, 2))
        for a in probes:
            result = normal_form({a: Q(1)}, polys, weight, 0)
            replay = tuple(sum(c*model['residue'](b)[i] for b, c in result.items()) for i in range(n))
            require(replay == model['residue'](a), "ordinary normal form versus direct truncated multiplication")
            forms += 1
    return {"full": full, "regional": regional, "normals": sorted(normals),
            "S_pairs": pairs, "normal_forms": forms}


def regional_control(r, s, d):
    params, weights = weight_list(r, s)
    rational = rectangle(r, s, coefficients(s, d), 0)
    characteristic = 101
    require(characteristic > d*(r-1)+s-1, "positive-characteristic qualification")
    finite = rectangle(r, s, coefficients(s, d), characteristic)
    seen, prefixes, selections = {}, 0, 0
    for rho, weight in zip(params, weights):
        basis = predicted_basis(r, s, d, rho)
        for field, model in ((0, rational), (characteristic, finite)):
            require(direct_greedy(model, weight, field) == basis,
                    "full direct monomial greedy agrees with the regional slot model")
            prefixes += block_controls(model, r, s, d, rho, field)
            selections += 1
        if basis not in seen:
            cone = rational_cone(rational, basis, weight, s)
            seen[basis] = {"basis": sorted(basis), "strict_weight": weight,
                           "full_cone": [[str(x) for x in p] for p in cone['full']],
                           "regional_cone": [[str(x) for x in p] for p in cone['regional']],
                           "ordinary_S_pairs": cone['S_pairs'], "ordinary_normal_forms": cone['normal_forms']}
    area = sum(area2([tuple(map(Q, p)) for p in record['regional_cone']]) for record in seen.values())
    require(area == area2(region_polygon(s)), "complete disjoint coverage of the prescribed weight cone")
    if d == 1:
        require(len(seen) == (r+1 if r <= s else s+2), "exact whole-outer-region bilinear count")
    return {"rectangle": [r, s], "d": d, "length": r*s,
            "characteristics_for_direct_controls": [0, characteristic],
            "common_list_weights": weights, "common_list_parameters": list(map(str, params)),
            "regional_initial_ideals": len(seen), "full_direct_selections": selections,
            "leading_prefix_ranks_checked": prefixes,
            "region_normalized_double_area": str(area),
            "ordinary_S_pairs": sum(c['ordinary_S_pairs'] for c in seen.values()),
            "ordinary_normal_forms": sum(c['ordinary_normal_forms'] for c in seen.values()),
            "witnesses": list(seen.values())}


def cubic_wall_control(a):
    q = 2*a-1
    params, weights = weight_list(q, q)
    threshold = Q(a-1, a)
    before = max(i for i, rho in enumerate(params) if rho < threshold)
    after = min(i for i, rho in enumerate(params) if rho > threshold)
    model = rectangle(q, q, coefficients(q, 1), 0)
    sides = []
    normal = (a, -(a-1), -a)
    wall_weight = (q*a+a-1, a, q*a)
    for index in (before, after):
        weight, rho = weights[index], params[index]
        basis = predicted_basis(q, q, 1, rho)
        require(direct_greedy(model, weight, 0) == basis, "actual rational quotient at both wall sides")
        cone = rational_cone(model, basis, weight, q)
        require(all(dot(z, wall_weight) >= 0 for z in cone['normals']), "shared weight is in both full cone closures")
        zero = [z for z in cone['normals'] if dot(z, wall_weight) == 0]
        require(zero and all(cross(z, normal) == (0, 0, 0) for z in zero),
                "the shared positive weight lies in the relative interior of one actual facet")
        sides.append({"basis": basis, "weight": weight, "zero_constraints": zero,
                      "S_pairs": cone['S_pairs'], "normal_forms": cone['normal_forms']})
    left, right = (x['basis'] for x in sides)
    expected_left = {(a, 0, t-a) for t in range(a, q)}
    expected_right = {(0, a-1, t) for t in range(a, q)}
    require(left-right == expected_left and right-left == expected_right, "all and only the predicted parallel swaps occur")
    difference = tuple(sum(x[j] for x in left)-sum(x[j] for x in right) for j in range(3))
    require(difference == tuple((a-1)*x for x in normal), "exact exposed profile difference")
    direction_product = prod(abs(x)+1 for x in normal)
    require(gcd(*normal) == 1 and direction_product == a*(a+1)**2, "primitive cubic direction product")
    return {"a": a, "side": q, "length": q*q, "primitive_normal": normal,
            "direction_product": direction_product, "product_divided_by_length": str(Q(direction_product, q*q)),
            "shared_strictly_positive_wall_weight": wall_weight,
            "before_weight": sides[0]['weight'], "after_weight": sides[1]['weight'],
            "before_basis": sorted(left), "after_basis": sorted(right),
            "removed_monomials": sorted(expected_left), "added_monomials": sorted(expected_right),
            "profile_difference": difference,
            "zero_border_constraints": [x['zero_constraints'] for x in sides],
            "ordinary_S_pairs": sum(x['S_pairs'] for x in sides),
            "ordinary_normal_forms": sum(x['normal_forms'] for x in sides)}


def polygon_prefix_controls():
    records = []
    for r, s in ((2, 2), (2, 3), (3, 2), (3, 3), (4, 4)):
        model = rectangle(r, s, coefficients(s, 1), 0)
        D = r+s-2
        nodes = [(U, V) for U in range(D+1) for V in range(D+1)
                 if -(s-1) <= U-V <= r-1]
        canonical = {(U, V): (max(U-V, 0), max(V-U, 0), min(U, V)) for U, V in nodes}
        columns = {(U, V): tuple(U**i*V**j for i in range(r) for j in range(s)) for U, V in nodes}
        base = D+2
        scale = (D+1)*(1+base+base**2)+1
        comparisons = []
        for coarse in ((3, 4, 6), (13, 10, 5), (9, 7, 3), (10, 13, 5), (3, 10, 5), (10, 3, 5)):
            weight = tuple(scale*x+y for x, y in zip(coarse, (1, base, base**2)))
            alpha, beta, gamma = weight
            require(gamma < alpha+beta, "strict carrier head is XY")
            def cost(node):
                U, V = node
                value = alpha*U+beta*V-(alpha+beta-gamma)*min(U, V)
                require(value == max(alpha*U+(gamma-alpha)*V, (gamma-beta)*U+beta*V)
                        == dot(weight, canonical[node]), "exact convex polygon cost identity")
                return value
            ordered = sorted(nodes, key=cost)
            require(len({cost(node) for node in nodes}) == len(nodes), "generic finite node order")
            ranks = rank_sequence([columns[node] for node in ordered], 0)
            selected = frozenset(canonical[node] for i, node in enumerate(ordered)
                                 if ranks[i] > (ranks[i-1] if i else 0))
            require(selected == direct_greedy(model, weight, 0),
                    "whole-fan polygon model versus original local multiplication")
            if gamma > max(alpha, beta):
                grid = frozenset((max(U-V, 0), max(V-U, 0), min(U, V))
                                 for U in range(r) for V in range(s))
                require(selected == grid, "constant middle cone by tensor interpolation")
            comparisons.append({'weight': weight, 'basis': sorted(selected),
                                'first_n_node_rank': ranks[r*s-1]})
        records.append({'rectangle': [r, s], 'finite_ground_nodes': len(nodes),
                        'comparisons': comparisons})
    model = rectangle(3, 3, coefficients(3, 1), 0)
    weight = (13, 10, 5)
    actual = direct_greedy(model, weight, 0)
    separated = predicted_basis(3, 3, 1, Q(8, 10))
    require(actual-separated == {(0, 0, 3)} and separated-actual == {(0, 2, 0)},
            "separated degree blocks cannot be assumed outside the proved region")
    return {'models': records, 'comparisons': sum(len(r['comparisons']) for r in records),
            'outside_region_guard': {'rectangle': [3, 3], 'weight': weight,
                                     'actual_only': sorted(actual-separated),
                                     'separated_formula_only': sorted(separated-actual)}}


def main():
    cases = ((2, 2, 1), (3, 3, 1), (4, 4, 1), (2, 4, 1), (4, 2, 1),
             (3, 4, 2), (4, 3, 2), (3, 3, 5), (2, 3, 4), (4, 4, 3),
             (5, 3, 2), (3, 5, 2))
    regions, walls = [], []
    for case in cases:
        regions.append(regional_control(*case))
        print(json.dumps({'regional_fixture_passed': case}), flush=True)
    for a in (2, 3, 4):
        walls.append(cubic_wall_control(a))
        print(json.dumps({'rational_wall_fixture_passed': a}), flush=True)
    polygon = polygon_prefix_controls()
    sources = ["docs/rectangle_boundary_fan.md", "scripts/run_rectangle_boundary_checks.py",
               "scripts/run_bilinear_rectangle_checks.py",
               "scripts/run_initial_cone_height_checks.py", "scripts/run_interleaving_research.py",
               "scripts/run_nested_jet_checks.py", "scripts/run_rectangle_rank_checks.py",
               "scripts/run_structural_complexity_checks.py"]
    report = {"schema_version": 1, "command": "python3 -B scripts/run_rectangle_boundary_checks.py",
              "scope": "Exact finite controls for the uniform regional prefix theorem, common list and infinite cubic-direction proof. Coverage is of the stated open cone, not the whole positive fan. No full asymptotic, novelty or completed-goal claim.",
              "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              "complete_regional_fans": regions, "actual_rational_cubic_walls": walls,
              "whole_bilinear_fan_polygon_reduction": polygon,
              "summary": {"complete_regional_models": len(regions),
                          "regional_initial_ideals": sum(r['regional_initial_ideals'] for r in regions),
                          "full_direct_quotient_selections": sum(r['full_direct_selections'] for r in regions),
                          "leading_prefix_rank_checks": sum(r['leading_prefix_ranks_checked'] for r in regions),
                          "actual_rational_cubic_walls": len(walls),
                          "whole_fan_polygon_comparisons": polygon['comparisons'],
                          "ordinary_S_pairs": sum(r['ordinary_S_pairs'] for r in regions+walls),
                          "ordinary_normal_forms": sum(r['ordinary_normal_forms'] for r in regions+walls)}}
    (ROOT/'results/research/rectangle_boundary_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['summary'], sort_keys=True))


if __name__ == '__main__':
    main()
