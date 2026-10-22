#!/usr/bin/env python3
"""Exact section controls; zonotopes here are explicitly not quotient algebras."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, gcd
from pathlib import Path
import ast
import hashlib
import json

from run_planar_dyadic_bound_checks import convex_hull
from run_structural_complexity_checks import cross, dot
from run_projection_tag_checks import fan, point_model, compare_fans
from run_coordinate_slice_budget_checks import edges_from_polygons, profile
from run_signed_face_checks import field_basis
from run_canonical_fan_deflation_checks import marked_basis
from run_algebraic_incidence_controls import pair_checks
from run_structural_complexity_checks import clip

ROOT = Path(__file__).resolve().parents[1]


def primitive(v, nonnegative=False):
    d = gcd(*v)
    if not d:
        return None
    v = tuple(t//d for t in v)
    if nonnegative:
        if max(v) <= 0:
            v = tuple(-t for t in v)
        return v if min(v) >= 0 else None
    return v if next(t for t in v if t) > 0 else tuple(-t for t in v)


def source_generators(k):
    """2(m+1)^2 times RWZ (2), d=3, m=4k+1; no row rescaling."""
    m = 4*k+1
    L = m+1
    out = [(2*L*L*(k-i), 2*L*L, 2*L*L) for i in range(2*k+1)]
    out += [(L*L*(2*i-2*k+1), -2*L*L, 2*L*L) for i in range(2*k)]
    out += [(2*(k-i), 0, 2*L) for i in range(2*k+1)]
    out += [(2*i-2*k+1, 0, -2*L) for i in range(2*k)]
    assert len(out) == 2*m
    assert len({primitive(g) for g in out}) == 2*m
    return out


def facet_section(gs, height, counts):
    """Independently build all facet polygons by planar Minkowski hulls."""
    normals = {primitive(cross(a, b)) for a, b in combinations(gs, 2)}-{None}
    section = set()
    for q0 in normals:
        for sign in (-1, 1):
            q = tuple(sign*t for t in q0)
            tied = [g for g in gs if dot(q, g) == 0]
            base = tuple(sum(g[i] for g in gs if dot(q, g) < 0) for i in range(3))
            drop = next(i for i in range(3) if q[i])
            other = [i for i in range(3) if i != drop]
            hull = [(0, 0)]
            for g in tied:
                hull = convex_hull(hull+[(p[0]+g[other[0]], p[1]+g[other[1]]) for p in hull])
            assert len(hull) >= 3
            points = []
            for p in hull:
                v = list(base)
                for i, a in zip(other, p):
                    v[i] += a
                v[drop] = Q(dot(q, base)-sum(q[i]*v[i] for i in other), q[drop])
                assert v[drop].denominator == 1
                points.append(tuple(v))
            hit = set()
            for a, b in zip(points, points[1:]+points[:1]):
                if min(a[2], b[2]) < height < max(a[2], b[2]):
                    f = (height-a[2])/(b[2]-a[2])
                    hit.add(tuple(a[i]+f*(b[i]-a[i]) for i in (0, 1)))
            assert len(hit) in (0, 2)
            section.update(hit)
            counts['independent_full_zonotope_facets'] += 1
    hull = convex_hull(section)
    assert set(hull) == section
    counts['independent_rational_central_section_hulls'] += 1
    return section


def positive_edges(gs, counts):
    """Sweep each zone plane, using integer rays and exact signed sums."""
    qs = {primitive(g) for g in gs if min(g) < 0 < max(g)}
    records = []
    profile_points = set()
    for q in sorted(qs):
        rays = {r for p in qs if p != q
                and (r := primitive(cross(q, p), True))}
        rays.update(r for i in range(3)
                    if (r := primitive(cross(q, tuple(int(i == j) for j in range(3))), True)))
        rays = sorted(rays, key=lambda r: tuple(Q(t, sum(r)) for t in r))
        for r, s in zip(rays, rays[1:]):
            w = tuple(r[i]*sum(s)+s[i]*sum(r) for i in range(3))
            if min(w) <= 0:
                continue
            assert dot(w, q) == 0
            a = [0]*3
            b = [0]*3
            for g in gs:
                t = dot(w, g)
                if t < 0:
                    for i in range(3):
                        a[i] += g[i]
                        b[i] += g[i]
                elif t == 0:
                    target = b if next(x for x in g if x) > 0 else a
                    for i in range(3):
                        target[i] += g[i]
            a, b = tuple(a), tuple(b)
            assert a != b and primitive(tuple(y-x for x, y in zip(a, b))) == q
            assert dot(w, a) == dot(w, b)
            records.append((a, b))
            profile_points.update((a, b))
            counts['independent_positive_zone_segments'] += 1
    assert len(records) == len({frozenset(e) for e in records})
    return records, profile_points


def one_case(k, counts):
    gs = source_generators(k)
    m = 4*k+1
    L = m+1
    # All vertex z coordinates are multiples of 2L; the central level is
    # an odd multiple of L. Moving it by 1/2 preserves every crossing.
    center = sum(g[2] for g in gs)//2
    assert center == m*L*L+L and center % (2*L) == L
    height = Q(2*center+1, 2)
    section = facet_section(gs, height, counts)
    assert len(section) == 2*m*(m+1)
    all_crossings = set()
    orthants = []
    for signs in product((-1, 1), repeat=3):
        transformed = [tuple(s*t for s, t in zip(signs, g)) for g in gs]
        edges, vertices = positive_edges(transformed, counts)
        length = sum(sum(abs(a-b) for a, b in zip(u, v)) for u, v in edges)
        width_sum = sum(sum(abs(g[i]) for g in transformed) for i in range(3))
        assert length <= len(gs)*width_sum
        hits = set()
        for a, b in edges:
            a = tuple(s*t for s, t in zip(signs, a))
            b = tuple(s*t for s, t in zip(signs, b))
            if min(a[2], b[2]) < height < max(a[2], b[2]):
                f = (height-a[2])/(b[2]-a[2])
                hits.add(tuple(a[i]+f*(b[i]-a[i]) for i in (0, 1)))
        assert hits <= section
        all_crossings.update(hits)
        # Explicit approaching-axis limits, with a small integer perturbation
        # to avoid all ties, independently witness large positive x width.
        limit = 1+sum(sum(map(abs, g)) for g in transformed)
        witnesses = []
        for axis in (0, 2):
            w = tuple(limit*limit if i == axis else (1 if i == 0 else limit if i == 1 else 2)
                      for i in range(3))
            assert all(dot(w, g) for g in transformed)
            witnesses.append(tuple(sum(g[i] for g in transformed if dot(w, g) < 0)
                                   for i in range(3)))
        assert abs(witnesses[0][0]-witnesses[1][0]) >= L*L*k*(2*k+1)
        orthants.append({'signs': list(signs), 'positive_section_crossings': len(hits),
                         'positive_edges': len(edges), 'L': length,
                         'axis_limit_x_width_witness': abs(witnesses[0][0]-witnesses[1][0])})
        counts['independent_positive_octants_and_integer_width_witnesses'] += 1
    assert all_crossings == section
    assert max(r['positive_section_crossings'] for r in orthants) >= len(section)/8
    lower = L*L*k*(2*k+1)
    assert lower > comb(len(gs), 2)
    return {'parameter': k, 'source_m': m, 'abstract_pair_matroid_rank': len(gs),
            'integer_generators': [list(g) for g in gs], 'half_integer_section_height': str(height),
            'section_vertices': len(section), 'positive_x_width_lower_bound': lower,
            'rank_only_downset_moment_cap': comb(len(gs), 2),
            'all_eight_positive_cost_orthants': orthants,
            'scope': 'Abstract paired-column matroid, not a polynomial quotient; '
                     'the displayed integer profiles violate the colength-rank moment cap.'}


def actual_pair_controls(counts):
    rows = []
    for s in range(1, 5):
        D = 2*s-1
        points = tuple((i, j, k) for i, j, k in product(range(2), range(2), range(D+1))
                       if (i+j+k) % 2)
        n, ev = point_model(points, 101)
        cells = fan(n, ev, 101, counts)
        pairs = [((0, b, c), (1, 1-b, D-c)) for b, c in product(range(2), range(D+1))]
        gs = [tuple(y-x for x, y in zip(a, b)) for a, b in pairs]
        base = tuple(sum(a[i] for a, b in pairs) for i in range(3))
        abstract, _ = positive_edges(gs, counts)
        predicted = {frozenset(tuple(v[i]+base[i] for i in range(3)) for v in edge)
                     for edge in abstract}
        direct = {frozenset(profile(B) for B in edge)
                  for edge in edges_from_polygons(cells, counts)}
        assert predicted == direct
        moments = {v for edge in direct for v in edge}
        assert all(sum(v) <= comb(n, 2) for v in moments)
        length = sum(sum(abs(a-b) for a, b in zip(*tuple(edge))) for edge in direct)
        assert len(cells) == (s+1)**2 and length == (2*s+1)*s*(s+2)
        assert len(cells) <= 1+n*(n+1)//2 and length <= 2*n*n*(n-1)
        section_rows = []
        integrated = 0
        for axis in range(3):
            events = Counter()
            for edge in direct:
                a, b = sorted(v[axis] for v in edge)
                if a < b:
                    events[a] += 1
                    events[b] -= 1
            keys = sorted(events)
            value = maximum = total = 0
            arg = None
            for a, b in zip(keys, keys[1:]):
                value += events[a]
                total += value*(b-a)
                if value > maximum:
                    maximum, arg = value, a
            integrated += total
            section_rows.append({'axis': axis, 'maximum_crossings': maximum,
                                 'half_integer_witness': str(Q(2*arg+1, 2)),
                                 'integral': total})
        assert integrated == length
        counts['actual_point_fans_vs_complete_pair_zone_profile_edges'] += 1
        rows.append({'parameter': s, 'characteristic': 101, 'length': n,
                     'initial_cones': len(cells), 'L': length, 'sections': section_rows,
                     'scope': 'Established actual four-line parity-point family; '
                              'finite controls, not evidence for a uniform section bound.'})
    return rows


def full_section(points, level):
    """Intersect every profile pair, independently of basis-edge structure."""
    out = {tuple(v[1:]) for v in points if v[0] == level}
    for a, b in combinations(points, 2):
        if min(a[0], b[0]) < level < max(a[0], b[0]):
            f = (level-a[0])/Q(b[0]-a[0])
            out.add(tuple(a[i]+f*(b[i]-a[i]) for i in (1, 2)))
    return convex_hull(out) if out else []


def section_identities(bases, labels, counts):
    moments = {B: tuple(sum(labels[e][i] for e in B) for i in range(3)) for B in bases}
    points = set(moments.values())
    lo = min(v[0] for v in points)
    hi = max(v[0] for v in points)
    integer = {}
    for k in range(lo, hi+1):
        same = [B for B in bases if moments[B][0] == k]
        assert same
        hull = convex_hull([moments[B][1:] for B in same])
        assert hull == full_section(points, Q(k))
        integer[k] = hull
        counts['integer_colour_sections_vs_all_profile_segments'] += 1
        sides = list(zip(hull, hull[1:]+hull[:1])) if len(hull) > 1 else []
        for a, b in sides:
            d = tuple(y-x for x, y in zip(a, b))
            on_face = [B for B in same
                       if d[0]*(moments[B][2]-a[1]) == d[1]*(moments[B][1]-a[0])]
            assert any(len(U ^ V) <= 4 and moments[U] != moments[V]
                       for U, V in combinations(on_face, 2))
            counts['projected_slice_edge_directions_with_at_most_two_exchanges'] += 1
    for k in range(lo, hi):
        for theta in (Q(1, 3), Q(1, 2), Q(2, 3)):
            weighted = convex_hull([tuple((1-theta)*a[i]+theta*b[i] for i in range(2))
                                    for a in integer[k] for b in integer[k+1]])
            assert weighted == full_section(points, k+theta)
            if theta == Q(1, 2):
                assert all((2*t).denominator == 1 for v in weighted for t in v)
            counts['adjacent_integer_minkowski_sections_vs_all_profile_segments'] += 1
    # Value-function convexity is checked independently at signed planar costs.
    for v in ((1, 1), (2, -3), (-3, 1)):
        values = [min(dot(v, p) for p in integer[k]) for k in range(lo, hi+1)]
        assert all(a-2*b+c >= 0 for a, b, c in zip(values, values[1:], values[2:]))
        counts['integer_colour_value_function_convexity'] += 1
    return integer


def binary_matroid_controls(counts):
    """Every nonempty simple binary rank-at-most-three restriction and colouring."""
    planar = ((0, 0), (1, 2), (2, 1), (3, 5), (5, 3), (4, 0), (0, 4))
    def rank(vectors):
        pivots = {}
        for v in vectors:
            while v:
                j = v.bit_length()-1
                if j in pivots:
                    v ^= pivots[j]
                else:
                    pivots[j] = v
                    break
        return len(pivots)
    models = colourings = 0
    for mask in range(1, 1 << 7):
        ground = tuple(i for i in range(7) if mask >> i & 1)
        r = rank([i+1 for i in ground])
        bases = [frozenset(B) for B in combinations(ground, r)
                 if rank([i+1 for i in B]) == r]
        for colour_mask in range(1 << len(ground)):
            labels = {i: (int(bool(colour_mask >> j & 1)), *planar[i]) for j, i in enumerate(ground)}
            polygons = section_identities(bases, labels, counts)
            lower_rank = rank([i+1 for i in ground if not labels[i][0]])
            upper_rank = rank([i+1 for i in ground if labels[i][0]])
            assert set(polygons) == set(range(r-lower_rank, upper_rank+1))
            assert len(polygons) == lower_rank+upper_rank-r+1
            counts['independent_colour_rank_interval_and_overlap_identities'] += 1
            colourings += 1
        models += 1
    assert models == 127 and colourings == 3**7-1
    return {'represented_matroids': models, 'binary_colourings': colourings,
            'field': 2, 'ground': 'All nonempty restrictions of the seven distinct nonzero GF(2)^3 columns.',
            'profile_labels': [list(p) for p in planar],
            'scope': 'Exhaustive bounded projected-identity controls; no uniform counting bound.'}


def actual_nonmatroid_sections(counts):
    rows = []
    triangle = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
    for r, p in product((2, 3, 5, 8), (0, 2, 3, 101)):
        A, B = r, r+1
        C = {(0, 0, 0), (1, 0, 0)}
        C.update((0, j, 0) for j in range(1, B))
        C.update((0, 0, j) for j in range(1, A))
        xy, xz, yB, zA = (1, 1, 0), (1, 0, 1), (0, B, 0), (0, 0, A)
        D = tuple(sorted(C | {xy, xz, yB, zA}))
        reference = tuple(sorted(C))+(xy, xz)
        index = {a: i for i, a in enumerate(reference)}
        n = A+B+2
        assert len(reference) == n
        def ev(a):
            a = tuple(a)
            if a == zA:
                a = xy
            elif a == yB:
                a = xz
            return {index[a]: 1} if a in index else {}
        bases = [frozenset(S) for S in combinations(D, n) if len(field_basis(S, ev, p)[0]) == n]
        assert len(bases) == 4
        expected = {}
        q1, q2 = (1, 1, -A), (1, -B, 1)
        for first, second in product((False, True), repeat=2):
            chosen = frozenset(C | {xy if first else zA, xz if second else yB})
            poly = list(triangle)
            for q, selected in ((q1, first), (q2, second)):
                poly = clip(poly, q, -1 if selected else 1)
            expected[chosen] = poly
        direct = fan(n, ev, p, counts)
        compare_fans(direct, expected, counts)
        assert set(bases) == set(direct)
        # Multiplication is checked on a spanning monomial basis and all carrier successors.
        def move(v, axis):
            result = Counter()
            for i, c in v.items():
                a = tuple(t+int(j == axis) for j, t in enumerate(reference[i]))
                for h, d in ev(a).items():
                    result[h] += c*d
            return {i: c % p if p else c for i, c in result.items() if (c % p if p else c)}
        for a in D:
            for i in range(3):
                successor = tuple(t+int(j == i) for j, t in enumerate(a))
                assert move(ev(a), i) == ev(successor)
                counts['section_family_actual_monomial_successors'] += 1
        for a in reference:
            for i, j in combinations(range(3), 2):
                assert move(move(ev(a), i), j) == move(move(ev(a), j), i)
                counts['section_family_commuting_coordinate_actions'] += 1
            assert not move(move(ev(a), 0), 0)
        for chosen, poly in direct.items():
            G, _ = marked_basis(chosen, ev, p)
            w = tuple(sum(v[j] for v in poly)/len(poly) for j in range(3))
            counts['section_family_original_reduced_basis_S_pairs'] += pair_checks(G, w, p)
        layer = [S for S in bases if sum(a[0] for a in S) == 2]
        S, T = frozenset(C | {xy, yB}), frozenset(C | {zA, xz})
        assert set(layer) == {S, T}
        # Direct exchange failure; excluded replacements are not filtered by an assumed model.
        failures = []
        for candidate in sorted(T-S):
            swapped = (S-{xy}) | {candidate}
            colour = sum(a[0] for a in swapped)
            rk = len(field_basis(tuple(sorted(swapped)), ev, p)[0])
            assert colour != 2 or rk < n
            failures.append({'replacement': list(candidate), 'upper_count': colour, 'rank': rk})
        d = (0, B+1, -A-1)
        assert gcd(*d) == 1
        assert tuple(a-b for a, b in zip(profile(S), profile(T))) == d
        assert not any(any(v) and not any(cross(v, d))
                       for a, b in combinations(D, 2)
                       for v in [tuple(x-y for x, y in zip(a, b))])
        polys = section_identities(bases, {a: a for a in D}, counts)
        assert set(polys[2]) == {profile(S)[1:], profile(T)[1:]}
        lower_rank = len(field_basis(tuple(a for a in D if not a[0]), ev, p)[0])
        upper_rank = len(field_basis(tuple(a for a in D if a[0]), ev, p)[0])
        assert (lower_rank, upper_rank) == (n-1, 3)
        assert set(polys) == set(range(n-lower_rank, upper_rank+1))
        counts['actual_integer_sections_with_failed_matroid_exchange'] += 1
        rows.append({'parameter': r, 'A': A, 'B': B, 'characteristic': p, 'length': n,
                     'layer_ranks': [lower_rank, upper_rank], 'overlap_rank': lower_rank+upper_rank-n,
                     'initial_cones': len(direct), 'integer_x_level': 2,
                     'two_actual_initial_bases': [list(map(list, sorted(U))) for U in (S, T)],
                     'section_primitive_direction': list(d), 'exchange_removed': list(xy),
                     'failed_single_exchanges': failures,
                     'scope': 'Actual all-field socle-pair quotient; no new large-fan lower bound.'})
    return rows


def main():
    if not __debug__:
        raise RuntimeError('Assertions must be enabled.')
    dependencies = ['scripts/run_profile_section_checks.py',
                    'scripts/run_planar_dyadic_bound_checks.py',
                    'scripts/run_structural_complexity_checks.py',
                    'scripts/run_projection_tag_checks.py',
                    'scripts/run_coordinate_slice_budget_checks.py',
                    'docs/ideal_step_incidence.md',
                    'docs/profile_section_budget.md']
    # Record every repository-local module in the import closure, including
    # the independent border-fan implementation used by the actual controls.
    pending = [p for p in dependencies if p.endswith('.py')]
    while pending:
        path = pending.pop()
        for node in ast.walk(ast.parse((ROOT/path).read_text())):
            names = ([node.module] if isinstance(node, ast.ImportFrom) and node.module
                     else [a.name for a in node.names] if isinstance(node, ast.Import) else [])
            for name in names:
                candidate = 'scripts/'+name.replace('.', '/')+'.py'
                if (ROOT/candidate).is_file() and candidate not in dependencies:
                    dependencies.append(candidate)
                    pending.append(candidate)
    dependencies.sort()
    before = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in dependencies}
    counts = Counter()
    rows = []
    for k in range(1, 6):
        row = one_case(k, counts)
        rows.append(row)
        print(json.dumps({key: row[key] for key in ('parameter', 'abstract_pair_matroid_rank',
                                                   'section_vertices', 'positive_x_width_lower_bound')}), flush=True)
    actual = actual_pair_controls(counts)
    binary = binary_matroid_controls(counts)
    obstruction = actual_nonmatroid_sections(counts)
    assert before == {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in dependencies}
    report = {'command': 'python3 -B scripts/run_profile_section_checks.py',
              'status': 'passed', 'source_sha256': before,
              'source': {'url': 'https://arxiv.org/pdf/0710.3116v3',
                         'scope': 'Equation (2), d=3, and Theorem 4.4; uniform integer scaling.'},
              'counts': dict(counts), 'cases': rows, 'actual_quotient_controls': actual,
              'binary_integer_sections': binary, 'actual_nonmatroid_section_families': obstruction,
              'additional_sources': [{'url': 'https://arxiv.org/pdf/math/0411260v1',
                                      'scope': 'Section 2, Definition 2.1, Proposition 2.4 and the following face-matroid discussion; classical edge and face inputs.'}],
              'limits': ['The source zonotopes are abstract matroids; no realization of those integer profiles or unrestricted section bound is established.',
                         'The cubic fixed-pair budget in the note is a classical arrangement corollary.',
                         'Integer colour sections are integral but need not be matroid profile polygons; the actual socle-family examples certify that obstruction.',
                         'The overlap-sensitive slice count is a classical restatement of the existing slice bound and rank identity, not a new unrestricted exponent.',
                         'Exact priority and the sharp unrestricted G3 endpoint remain open.']}
    target = ROOT/'results/research/profile_section_checks.json'
    target.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
