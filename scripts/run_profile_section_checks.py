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
from run_projection_tag_checks import fan, point_model
from run_coordinate_slice_budget_checks import edges_from_polygons, profile

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
    assert before == {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in dependencies}
    report = {'command': 'python3 -B scripts/run_profile_section_checks.py',
              'status': 'passed', 'source_sha256': before,
              'source': {'url': 'https://arxiv.org/pdf/0710.3116v3',
                         'scope': 'Equation (2), d=3, and Theorem 4.4; uniform integer scaling.'},
              'counts': dict(counts), 'cases': rows, 'actual_quotient_controls': actual,
              'limits': ['No algebra realization or unrestricted section bound is established.',
                         'The cubic fixed-pair budget in the note is a classical arrangement corollary.',
                         'Exact priority and the sharp unrestricted G3 endpoint remain open.']}
    target = ROOT/'results/research/profile_section_checks.json'
    target.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
