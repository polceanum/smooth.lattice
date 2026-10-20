#!/usr/bin/env python3
"""Independent exact controls for the bilinear rectangle/torsion dictionary.

The uniform proof is in docs/bilinear_rectangle_fans.md. These controls do
not prove priority, an asymptotic characteristic-zero law, or the G_3 goal.
"""
from fractions import Fraction as Q
from heapq import heappop, heappush
from itertools import combinations, product
from math import comb, factorial, gcd
from pathlib import Path
import hashlib
import json

from run_extremal_initial_ideals_checks import rectangle_fan
from run_initial_cone_height_checks import add, leading, normal_form, subtract_shift
from run_interleaving_research import minimal_excluded_border
from run_nested_jet_checks import inverse
from run_rectangle_rank_checks import rectangle
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def labels(a, r, s):
    return ((a[0]+a[2]) % r, (a[1]+a[2]) % s)


def witnesses(r, s):
    require(s % r == 0, "divisibility is a theorem hypothesis")
    if r == 1:
        return [(1, 1, 1)] if s == 1 else [(1, 1, 2), (1, 2, 1)]
    weights = [(1, 1, 3)]
    for k in range(r):
        for ell in range(k+1):
            if r == s:
                weights.append((4*k+3, 4*r-4*ell-1, 4*k-4*ell+2))
            else:
                u, w = s*(4*k+3), r*(4*k-4*ell+2)
                weights.append((u, 4*r*s-u+w, w))
    if r < s:
        weights.extend((4*s, 2*j+1, 2*s+2*j+1) for j in range(r))
    require(all(min(w) > 0 and max(w) <= 8*r*s for w in weights),
            "complete witnesses have the asserted positive height bound")
    return weights


def two_candidates(r, s, weight):
    u, v, w = weight
    if r == 1:
        require(s == 1 or v != w, "strict one-dimensional region")
        return frozenset((0, b, 0) if v < w else (0, 0, b) for b in range(s))
    delta = u+v-w
    if delta < 0:
        return frozenset((a, b, 0) for a in range(r) for b in range(s))
    require(delta > 0, "strict mixed-product region")
    result = set()
    for alpha, beta in product(range(r), range(s)):
        if r == s or w < v:
            d = (beta-alpha) % r
            options = [((-d) % r, 0, beta)]
            if d:
                options.append((0, d, (beta-d) % s))
            costs = [dot(weight, a) for a in options]
            require(len(set(costs)) == len(costs), "two candidates have strict costs")
            result.add(options[costs.index(min(costs))])
        else:
            require(w > v, "unequal periods avoid the Y/Z wall")
            c = alpha if beta >= alpha or (alpha-beta)*delta > s*v else beta
            result.add(((alpha-c) % r, (beta-c) % s, c))
    return frozenset(result)


def exhaustive_representatives(r, s, weight):
    answer = set()
    for alpha, beta in product(range(r), range(s)):
        candidates = [((alpha-c) % r, (beta-c) % s, c) for c in range(s)]
        costs = [dot(weight, a) for a in candidates]
        require(costs.count(min(costs)) == 1, "exhaustive representative is unique")
        answer.add(candidates[costs.index(min(costs))])
    return frozenset(answer)


def shortest_paths(r, s, weight):
    best = {(0, 0): (0, (0, 0, 0))}
    heap = [(0, (0, 0, 0), (0, 0))]
    steps = ((1, 0), (0, 1), (1, 1))
    while heap:
        cost, exponent, residue = heappop(heap)
        if best[residue] != (cost, exponent):
            continue
        for i, step in enumerate(steps):
            node = ((residue[0]+step[0]) % r, (residue[1]+step[1]) % s)
            word = tuple(x+(i == j) for j, x in enumerate(exponent))
            candidate = (cost+weight[i], word)
            if node not in best or candidate < best[node]:
                best[node] = candidate
                heappush(heap, (*candidate, node))
    require(len(best) == r*s, "Cayley graph has all residue classes")
    return frozenset(word for _, word in best.values())


def cone_and_basis(basis, r, s, weight):
    by_label = {labels(a, r, s): a for a in basis}
    require(len(by_label) == r*s, "basis has one representative per class")
    require(all(tuple(x-(i == j) for i, x in enumerate(a)) in basis
                for a in basis for j in range(3) if a[j]), "standard set is a downset")
    polys, normals = [], []
    for a in minimal_excluded_border(sorted(basis)):
        b = by_label[labels(a, r, s)]
        z = tuple(x-y for x, y in zip(a, b))
        require(dot(weight, z) > 0, "all reduced binomial heads are strict")
        polys.append({a: Q(1), b: Q(-1)})
        g = gcd(*z)
        normals.append(tuple(x//g for x in z))
    poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
    for z in normals:
        poly = clip(poly, z, 1)
        require(poly, "nonempty complete binomial cone")
    pairs = 0
    for g, h in combinations(polys, 2):
        a, b = leading(g, weight), leading(h, weight)
        common = tuple(max(x, y) for x, y in zip(a, b))
        first = {add(c, tuple(x-y for x, y in zip(common, a))): value for c, value in g.items()}
        spair = subtract_shift(first, h, tuple(x-y for x, y in zip(common, b)), Q(1), 0)
        require(not normal_form(spair, polys, weight, 0), "ordinary rational binomial S-pair")
        pairs += 1
    checks = 0
    for a in hyperbolic_cross(r*s, 3):
        require(normal_form({a: Q(1)}, polys, weight, 0) == {by_label[labels(a, r, s)]: Q(1)},
                "ordinary division versus modular representative formula")
        checks += 1
    return poly, pairs, checks, len(polys)


def torsion_controls(r, s):
    weights = witnesses(r, s)
    seen, area, pairs, forms = {}, Q(0), 0, 0
    for w in weights:
        basis = two_candidates(r, s, w)
        require(basis == shortest_paths(r, s, w) == exhaustive_representatives(r, s, w),
                "formula versus shortest paths and every periodic candidate")
        require(basis not in seen, "explicit weights realize distinct states")
        poly, count, checks, generators = cone_and_basis(basis, r, s, w)
        seen[basis] = {"basis": sorted(basis), "weight": w,
                       "polygon": [[str(x) for x in p] for p in poly],
                       "reduced_generators": generators}
        area += area2(poly)
        pairs += count
        forms += checks
    require(area == 1, "complete disjoint positive fan has exact area one")
    expected = (1 if s == 1 else 2 if r == 1 else
                r*(r+1)//2+1 if r == s else (r+1)*(r+2)//2)
    require(len(seen) == expected, "proved exact count versus the complete fan")
    return {"rectangle": [r, s], "length": r*s, "all_fields": True,
            "initial_ideal_count": len(seen), "normalized_double_area": "1",
            "ordinary_S_pairs": pairs, "ordinary_ground_normal_forms": forms,
            "witnesses": list(seen.values())}, seen


def rank_sequence(columns, characteristic):
    pivots, ranks = {}, []
    for column in columns:
        v = [x % characteristic if characteristic else Q(x) for x in column]
        for i, row in sorted(pivots.items()):
            c = v[i]
            v = [x-c*y for x, y in zip(v, row)]
            if characteristic:
                v = [x % characteristic for x in v]
        i = next((i for i, x in enumerate(v) if x), None)
        if i is not None:
            inverse = pow(v[i], -1, characteristic) if characteristic else 1/v[i]
            pivots[i] = [(x*inverse) % characteristic if characteristic else x*inverse for x in v]
        ranks.append(len(pivots))
    return tuple(ranks)


def translation_controls(r, s, characteristic, weights):
    local = rectangle(r, s, {(1, 0): 1, (0, 1): 1, (1, 1): 1}, characteristic)
    ground = tuple(hyperbolic_cross(r*s, 3))
    cells = local['cells']
    old = {a: local['residue'](a) for a in ground}
    pascal = {a: tuple(comb(a[0]+a[2], i)*comb(a[1]+a[2], j) for i, j in cells) for a in ground}
    if characteristic:
        pascal = {a: tuple(x % characteristic for x in col) for a, col in pascal.items()}
    for a in ground:
        expansion = [0]*(r*s)
        for b in product(*(range(x+1) for x in a)):
            factor = comb(a[0], b[0])*comb(a[1], b[1])*comb(a[2], b[2])
            for i, c in enumerate(old[b]):
                expansion[i] += factor*c
        if characteristic:
            expansion = [x % characteristic for x in expansion]
        require(tuple(expansion) == pascal[a], "integer translation expansion versus Pascal formula")
    sequences = 0
    for weight in weights:
        ordered = sorted(ground, key=lambda a: (dot(weight, a), a))
        a = rank_sequence([old[x] for x in ordered], characteristic)
        b = rank_sequence([pascal[x] for x in ordered], characteristic)
        require(a == b, "every downset-prefix rank is translation invariant")
        selected = frozenset(x for i, x in enumerate(ordered) if a[i] > (a[i-1] if i else 0))
        if characteristic:
            require(selected == two_candidates(r, s, weight), "actual local quotient has the classified fan")
        sequences += len(ordered)
    if characteristic:
        refs = {(i, j): tuple(comb(i, a)*comb(j, b) % characteristic for a, b in cells)
                for i, j in product(range(r), range(s))}
        require(rank_sequence(refs.values(), characteristic)[-1] == r*s, "tensor Pascal reference is invertible")
        require(all(pascal[a] == refs[labels(a, r, s)] for a in ground), "Lucas/Frobenius residue collapse")
    # This is deliberately not an arbitrary-subset rank identity.
    zero = (r, 0, 0)
    require(not any(local['residue'](zero)) and any(
        comb(r, i)*comb(0, j) % characteristic if characteristic else comb(r, i)*comb(0, j)
        for i, j in cells), "translation changes some individual monomial ranks")
    return {"rectangle": [r, s], "characteristic": characteristic,
            "all_ground_translation_expansions": len(ground),
            "prefix_rank_equalities": sequences, "weights_checked": len(weights),
            "individual_column_rank_guard": {"exponent": zero, "before_rank": 0, "after_rank": 1}}


def rational_fan_control(q):
    modular = rectangle_fan(q, q, {(1, 0): 1, (0, 1): 1, (1, 1): 1}, 1000003)
    model = rectangle(q, q, {(1, 0): 1, (0, 1): 1, (1, 1): 1}, 0)
    area, pairs, forms = Q(0), 0, 0
    for witness in modular['witnesses']:
        basis = tuple(map(tuple, witness['basis']))
        weight = tuple(witness['strict_weight'])
        matrix = [[model['residue'](a)[i] for a in basis] for i in range(q*q)]
        inv = inverse(matrix)
        polys, normals = [], set()
        for head in minimal_excluded_border(basis):
            col = model['residue'](head)
            coefficients = [sum(x*y for x, y in zip(row, col)) for row in inv]
            require(tuple(sum(c*model['residue'](a)[i] for c, a in zip(coefficients, basis))
                          for i in range(q*q)) == col, "rational border coordinates replay independently")
            polynomial = {head: Q(1)}
            for a, c in zip(basis, coefficients):
                if c:
                    polynomial[a] = -c
                    z = tuple(x-y for x, y in zip(head, a))
                    g = gcd(*z)
                    normals.add(tuple(x//g for x in z))
            polys.append(polynomial)
        poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
        for z in sorted(normals):
            require(dot(weight, z) > 0, "rational reduced heads are strict")
            poly = clip(poly, z, 1)
        require(set(poly) == {tuple(map(Q, p)) for p in witness['cone_polygon']},
                "exact rational cone equals the whole modular cone")
        area += area2(poly)
        for g, h in combinations(polys, 2):
            a, b = leading(g, weight), leading(h, weight)
            common = tuple(max(x, y) for x, y in zip(a, b))
            first = {add(c, tuple(x-y for x, y in zip(common, a))): value for c, value in g.items()}
            spair = subtract_shift(first, h, tuple(x-y for x, y in zip(common, b)), Q(1), 0)
            require(not normal_form(spair, polys, weight, 0), "ordinary rational S-pair certifies each cone")
            pairs += 1
        for a in model['ground']:
            remainder = normal_form({a: Q(1)}, polys, weight, 0)
            replay = tuple(sum(c*model['residue'](b)[i] for b, c in remainder.items()) for i in range(q*q))
            require(replay == model['residue'](a), "ordinary rational normal form versus direct local products")
            forms += 1
    require(area == 1, "all rational initial cones cover the positive simplex")
    return {"side": q, "length": q*q, "characteristic": 0,
            "complete_rational_fan_count": len(modular['witnesses']),
            "ordinary_rational_S_pairs": pairs, "ordinary_rational_normal_forms": forms,
            "all_rational_cones_match_modular_discovery": True,
            "normalized_double_area": "1", "witnesses": modular['witnesses']}


def interpolation_control(r, s):
    cells = tuple(product(range(r), range(s)))
    ground = tuple(hyperbolic_cross(r*s, 3))
    univariate = []
    for degree in range(max(r, s)):
        coefficients = [Q(1)]
        for root in range(degree):
            result = [Q(0)]*(len(coefficients)+1)
            for i, c in enumerate(coefficients):
                result[i] -= root*c
                result[i+1] += c
            coefficients = result
        univariate.append([c/factorial(degree) for c in coefficients])
    change = [[(univariate[i][k] if k <= i else 0)*(univariate[j][ell] if ell <= j else 0)
               for k, ell in cells] for i, j in cells]
    for a in ground:
        u, v = a[0]+a[2], a[1]+a[2]
        ordinary = [u**i*v**j for i, j in cells]
        result = [sum(x*y for x, y in zip(row, ordinary)) for row in change]
        require(result == [comb(u, i)*comb(v, j) for i, j in cells],
                "exact falling-factorial row change gives every tensor interpolation column")
    require(rank_sequence(zip(*change), 0)[-1] == r*s, "row change is invertible in characteristic zero")
    return {"rectangle": [r, s], "integer_nodes_checked": len(ground),
            "tensor_row_change": [[str(x) for x in row] for row in change]}


def main():
    fixtures = ((1, 1), (1, 5), (2, 2), (2, 4), (2, 6), (3, 3),
                (3, 6), (4, 4), (4, 8), (5, 5), (5, 10), (6, 6))
    torsion = []
    for r, s in fixtures:
        report, _ = torsion_controls(r, s)
        torsion.append(report)
    translations = []
    for r, s, p in ((2, 2, 2), (2, 4, 2), (4, 4, 2), (4, 8, 2), (3, 3, 3), (3, 9, 3), (5, 5, 5)):
        translations.append(translation_controls(r, s, p, witnesses(r, s)))
    translations.append(translation_controls(3, 3, 0, [(1, 1, 1), (1, 7, 3), (9, 2, 5)]))
    rational = [rational_fan_control(q) for q in (2, 3, 4)]
    interpolation = [interpolation_control(r, s) for r, s in ((2, 2), (2, 3), (3, 3))]
    sources = ["docs/bilinear_rectangle_fans.md", "scripts/run_bilinear_rectangle_checks.py",
               "scripts/run_extremal_initial_ideals_checks.py", "scripts/run_initial_cone_height_checks.py",
               "scripts/run_nested_jet_checks.py", "scripts/run_rectangle_rank_checks.py",
               "scripts/run_structural_complexity_checks.py", "scripts/run_interleaving_research.py"]
    report = {"schema_version": 1, "command": "python3 -B scripts/run_bilinear_rectangle_checks.py",
              "scope": "Finite exact controls for the uniform divisible-period torsion count, its prime-power local consequence, and the translation/interpolation dictionary. Classical inputs are attributed. No unrestricted exponent, characteristic-zero asymptotic fit or priority claim.",
              "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              "torsion_fans": torsion, "local_translation_controls": translations,
              "complete_rational_fans": rational, "tensor_interpolation_controls": interpolation,
              "summary": {"complete_all_field_torsion_fans": len(torsion),
                          "torsion_initial_ideals": sum(r['initial_ideal_count'] for r in torsion),
                          "ordinary_binomial_S_pairs": sum(r['ordinary_S_pairs'] for r in torsion),
                          "ordinary_binomial_normal_forms": sum(r['ordinary_ground_normal_forms'] for r in torsion),
                          "local_translation_models": len(translations),
                          "downset_prefix_rank_equalities": sum(r['prefix_rank_equalities'] for r in translations),
                          "complete_rational_fans": len(rational),
                          "rational_initial_ideals": sum(r['complete_rational_fan_count'] for r in rational),
                          "ordinary_rational_S_pairs": sum(r['ordinary_rational_S_pairs'] for r in rational),
                          "ordinary_rational_normal_forms": sum(r['ordinary_rational_normal_forms'] for r in rational)}}
    (ROOT/'results/research/bilinear_rectangle_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['summary'], sort_keys=True))


if __name__ == '__main__':
    main()
