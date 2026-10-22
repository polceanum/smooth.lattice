#!/usr/bin/env python3
"""Exact integer-point evaluation, paired columns and whole update fans."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import comb, factorial, gcd, lcm
from pathlib import Path
import hashlib
import json
import sys

from run_projection_tag_checks import fan, point_model, compare_fans
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import area2, clip, cross, dot
from run_canonical_fan_deflation_checks import marked_basis
from run_algebraic_incidence_controls import pair_checks

ROOT = Path(__file__).resolve().parents[1]
TRIANGLE = [tuple(Q(i == j) for i in range(3)) for j in range(3)]


def scalar(v, p):
    return v.numerator*pow(v.denominator, -1, p) % p if p else v


def kraw(c, k, D):
    return Q(sum((-1)**h*comb(k, h)*comb(D-k, c-h)
                 for h in range(c+1) if h <= k and c-h <= D-k), comb(D, c))


def kraw_coefficients(c, D):
    values = [kraw(c, k, D) for k in range(D+1)]
    differences = []
    while values:
        differences.append(values[0])
        values = [b-a for a, b in zip(values, values[1:])]
    assert all(v == 0 for v in differences[c+1:])
    assert differences[c] == Q((-2)**c, comb(D, c))
    polynomial = [Q(0)]*(c+1)
    newton = [Q(1)]
    for j in range(c+1):
        for i, v in enumerate(newton):
            polynomial[i] += differences[j]*v
        nxt = [Q(0)]*(len(newton)+1)
        for i, v in enumerate(newton):
            nxt[i] -= j*v/(j+1)
            nxt[i+1] += v/(j+1)
        newton = nxt
    assert polynomial[-1] == Q((-2)**c, comb(D, c)*factorial(c))
    return polynomial


def primitive(v):
    den = lcm(*(a.denominator if isinstance(a, Q) else 1 for a in v))
    v = tuple(int(a*den) for a in v)
    d = gcd(*v)
    v = tuple(a//d for a in v)
    return v if next(a for a in v if a) > 0 else tuple(-a for a in v)


def clean(poly):
    poly = list(dict.fromkeys(poly))
    while len(poly) > 2:
        keep = [a for i, a in enumerate(poly)
                if any(cross(tuple(x-y for x, y in zip(a, poly[i-1])),
                             tuple(x-y for x, y in zip(poly[(i+1) % len(poly)], a))))]
        if len(keep) == len(poly):
            break
        poly = keep
    return poly


def inventory(cells):
    walls, rays, incidences = {}, set(), Counter()
    for poly in cells.values():
        P = clean(poly)
        rays.update(primitive(v) for v in P if min(v) > 0)
        for a, b in zip(P, P[1:]+P[:1]):
            if any(a[i] == b[i] == 0 for i in range(3)):
                continue
            edge = tuple(sorted((a, b)))
            q = primitive(cross(a, b))
            if edge in walls:
                assert walls[edge] == q
            walls[edge] = q
            incidences[edge] += 1
    assert set(incidences.values()) == {2}
    directions = Counter(walls.values())
    W = sum(k*sum(map(abs, q)) for q, k in directions.items())
    return directions, rays, W


def clip_cone(normals):
    poly = list(TRIANGLE)
    for q in normals:
        poly = clip(poly, q, 1)
    assert poly and area2(poly) > 0
    return clean(poly)


def predicted_polygons(D, added):
    if not added:
        polys = [list(TRIANGLE)]
        normals = [(1, 1, -a) for a in range(1, D+1, 2)]
        normals += [(1, -1, -b) for b in range(-D, D+1, 2)]
        for q in normals:
            polys = [clean(P) for poly in polys for sign in (-1, 1)
                     if (P := clip(poly, q, sign)) and area2(P) > 0]
        return polys
    polys = []
    for i in range(D):
        for j in range(D-i):
            polys.append(clip_cone([(1, 0, -i), (-1, 0, i+1),
                                    (0, 1, -j), (0, -1, j+1)]))
    for i in range(D+1):
        b = 2*i-D
        normals = [(1, 0, -i), (0, 1, -(D-i))]
        if i:
            normals.append((1, -1, -(b-1)))
        if i < D:
            normals.append((-1, 1, b+1))
        polys.append(clip_cone(normals))
    return polys


def paired_basis(box, w, added):
    D = max(a[2] for a in box)
    pairs = [((0, b, c), (1, 1-b, D-c))
             for b, c in product(range(2), range(D+1))]
    ordered = [sorted(pair, key=lambda a: (dot(w, a), a)) for pair in pairs]
    chosen = {pair[0] for pair in ordered}
    if added:
        chosen.add(min((pair[1] for pair in ordered), key=lambda a: (dot(w, a), a)))
    return frozenset(chosen)


def one_case(s, p, added, counts):
    D = 2*s-1
    points = [(i, j, k) for i, j, k in product(range(2), range(2), range(D+1))
              if (i+j+k) % 2]
    if added:
        points.append((0, 0, 0))
    n, ev = point_model(tuple(points), p)
    assert n == 4*s+int(added)
    box = tuple(product(range(2), range(2), range(D+1)))
    polys = {c: kraw_coefficients(c, D) for c in range(D+1)}

    @lru_cache(None)
    def qev(a):
        coefficients = ([Q(1)] if a[0] == 0 else [Q(1), Q(-2)],
                        [Q(1)] if a[1] == 0 else [Q(1), Q(-2)], polys[a[2]])
        result = defaultdict(Q)
        for powers in product(*(range(len(v)) for v in coefficients)):
            factor = coefficients[0][powers[0]]*coefficients[1][powers[1]]*coefficients[2][powers[2]]
            for j, v in ev(powers).items():
                result[j] += factor*v
        result = {j: scalar(v, p) for j, v in result.items() if scalar(v, p)}
        direct = {j: scalar(Q((1-2*x)**a[0]*(1-2*y)**a[1])*kraw(a[2], z, D), p)
                  for j, (x, y, z) in enumerate(points)}
        direct = {j: v for j, v in direct.items() if v}
        assert result == direct
        counts['independent_triangular_polynomial_evaluation_columns'] += 1
        return result

    reference = [(0, b, c) for b, c in product(range(2), range(D+1))]
    if added:
        reference.append((1, 1, D))
    _, original_pivots = field_basis(reference, qev, p)
    assert len(original_pivots) == n
    def structured(a):
        i, b, c = a
        r = (b, c) if i == 0 else (1-b, D-c)
        out = {r: 1 if i == 0 else -1}
        if added:
            out[2, 0] = 1
        return {j: scalar(Q(v), p) for j, v in out.items()}
    _, paired_pivots = field_basis(reference, structured, p)
    for a in box:
        rem, coordinates = field_reduce(qev(a), original_pivots, p)
        rem2, coordinates2 = field_reduce(structured(a), paired_pivots, p)
        assert not rem and not rem2 and coordinates == coordinates2
        counts['independent_full_paired_rank_model_columns'] += 1
    if not added:
        for a, b in product(reference, repeat=2):
            value = sum(Q(comb(D, k))*Q((1-2*j)**(a[1]+b[1]))*
                        kraw(a[2], k, D)*kraw(b[2], k, D) for i, j, k in points)
            expected = Q(2**(D+1), comb(D, a[2])) if a == b else 0
            assert value == expected
            counts['independent_restricted_weighted_Gram_entries'] += 1
    else:
        separator = [Q(1+(1-2*i)*(1-2*j)*(-1)**k, 2) for i, j, k in points]
        assert separator == [0]*(n-1)+[1]
        for axis in range(3):
            assert all(point[axis]*v == 0 for point, v in zip(points, separator))
            counts['independent_actual_corner_separator_actions'] += 1

    predicted = {}
    for poly in predicted_polygons(D, added):
        w = tuple(sum(v[j] for v in poly)/len(poly) for j in range(3))
        B = paired_basis(box, w, added)
        direct = frozenset(field_basis(sorted(box, key=lambda a: (dot(w, a), a)), ev, p)[0])
        assert direct == B and len(B) == n and B not in predicted
        predicted[B] = poly
        counts['independent_entire_grid_cells_vs_original_point_greedy'] += 1
    actual = fan(n, ev, p, counts)
    compare_fans(actual, predicted, counts)
    for B, poly in actual.items():
        G, _ = marked_basis(B, ev, p)
        assert len(G) <= 6
        w = tuple(sum(v[j] for v in poly)/len(poly) for j in range(3))
        counts['independent_original_reduced_basis_S_pairs'] += pair_checks(G, w, p)
    directions, rays, W = inventory(actual)
    if added:
        expected = Counter({(1, 0, -i): D+1-i for i in range(1, D+1)})
        expected.update({(0, 1, -i): D+1-i for i in range(1, D+1)})
        expected.update({(1, -1, -b): 1 for b in range(-D+1, D, 2)})
        expected_rays = {(i, j, 1) for i in range(1, D+1) for j in range(1, D+2-i)}
        assert len(actual) == s*(2*s+1)
        assert W == Q(8, 3)*s**3+6*s*s-Q(2, 3)*s-2
    else:
        expected = Counter({(1, 1, -a): a for a in range(1, D+1, 2)})
        expected.update({(1, -1, -b): (D-abs(b))//2+1 for b in range(-D, D+1, 2)})
        expected_rays = {(i, j, 1) for i in range(1, D+1) for j in range(1, D+1-i) if (i+j) % 2}
        assert len(actual) == (s+1)**2 and W == (2*s+1)*s*(s+2)
    assert directions == expected and rays == expected_rays
    counts['independent_complete_primitive_wall_multiplicity_inventories'] += 1
    return {'parameter': s, 'characteristic': p, 'point_added': added, 'length': n,
            'initial_cones': len(actual), 'interior_rays': len(rays), 'W': W,
            'undirected_wall_pieces': sum(directions.values())}


def main():
    if not __debug__:
        raise RuntimeError('Exact assertions must be enabled.')
    proofs = ['docs/parity_point_fan_updates.md']
    before = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT/'scripts').glob('*.py')}
    before.update({f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in proofs})
    counts = Counter()
    rows = []
    for p, limit in ((0, 3), (101, 6)):
        for s in range(1, limit+1):
            pair = [one_case(s, p, added, counts) for added in (False, True)]
            assert pair[1]['initial_cones']-pair[0]['initial_cones'] == s*s-s-1
            assert pair[1]['W']-pair[0]['W'] == (2*s**3+3*s*s-8*s-6)//3
            rows.extend(pair)
            print(pair, flush=True)
    for s in range(1, 257):
        D = 2*s-1
        wc = sum((D+2)*(a+2) for a in range(1, D+1, 2))
        wp = 2*sum((D+1-i)*(i+1) for i in range(1, D+1))
        wp += sum(abs(b)+2 for b in range(-D+1, D, 2))
        assert wc == (2*s+1)*s*(s+2)
        assert wp-wc == (2*s**3+3*s*s-8*s-6)//3
        if s >= 4:
            assert 3*(wp-wc) >= 2*s**3
        counts['independent_exact_large_parameter_wall_sum_identities'] += 1
    loaded = {Path(m.__file__).resolve() for m in tuple(sys.modules.values())
              if getattr(m, '__file__', None) and Path(m.__file__).resolve().parent == ROOT/'scripts'}
    loaded.add(Path(__file__).resolve())
    sources = proofs+sorted(str(p.relative_to(ROOT)) for p in loaded)
    assert all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest() == before[f] for f in sources)
    report = {'all_checks_passed': True, 'scope': 'Complete original integer point-ideal fans and explicit paired rank models, full cone equality, ordinary S-pairs, weighted orthogonality and exact wall inventories. Finite checks support the uniform formulas and do not establish their asymptotic scope or exact priority. No unrestricted G3 improvement or polynomial separation is claimed.',
              'point_ideal_examples': rows, 'counts': dict(counts),
              'loaded_source_stability_checked': True,
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    target = ROOT/'results/research/parity_point_fan_checks.json'
    target.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'all_checks_passed': True, 'counts': dict(counts),
                      'report': str(target.relative_to(ROOT))}), flush=True)


if __name__ == '__main__':
    main()
