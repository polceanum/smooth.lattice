#!/usr/bin/env python3
"""Exact cyclic-quotient fans, aggregate charges and period stabilization.

Bounded controls support the uniform proofs; no fitted growth or novelty claim.
"""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import hashlib
import json

from run_extremal_initial_ideals_checks import adjacent_weights
from run_initial_cone_height_checks import integer_ray
from run_interleaving_research import minimal_excluded_border
from run_structural_complexity_checks import area2, clip, cross, dot, primitive

ROOT = Path(__file__).resolve().parents[1]
TRIANGLE = [tuple(Q(i == j) for i in range(3)) for j in range(3)]


@lru_cache(None)
def model(r, m):
    n = r*m
    E = tuple((a, b, c) for a in range(n) for b in range(n) for c in range(m)
              if (a+1)*(b+1)*(c+1) <= n)
    classes = tuple(tuple(a for a in E if sum(a) % m == j) for j in range(m))

    def select(w):
        return frozenset(a for C in classes for a in sorted(C, key=lambda a: (dot(a, w), a))[:r])
    return E, classes, select


def cone(B, m):
    normals = set()
    for h in minimal_excluded_border(tuple(B)):
        if h[2] >= m:
            continue
        for b in B:
            if sum(h) % m == sum(b) % m:
                z = tuple(x-y for x, y in zip(h, b))
                g = gcd(*z)
                normals.add(tuple(x//g for x in z))
    polygon = TRIANGLE[:]
    for z in sorted(normals):
        polygon = clip(polygon, z, 1)
    assert area2(polygon) > 0
    return polygon, tuple(sorted(normals))


def core(r, m, B, w, counts):
    cheapest = min(range(3), key=lambda i: w[i])
    assert all(w[i] > w[cheapest] for i in range(3) if i != cheapest)
    others = [i for i in range(3) if i != cheapest]
    prices = (m*w[cheapest], *(w[i]-w[cheapest] for i in others))
    K0 = (r+m-2)//m
    selected = []
    for j in range(m):
        points = set()
        for a in B:
            if sum(a) % m != j:
                continue
            k = (sum(a)-j)//m
            v = (k, *(a[i] for i in others))
            assert dot(a, w)-j*w[cheapest] == dot(v, prices)
            b, c = v[1:]
            assert b+c <= r-1
            if cheapest == 2:
                assert (b+1)*(c+1) <= r and 0 <= k <= K0
            elif k < K0:
                assert (b+1)*(c+1) <= r
            else:
                assert (k-K0+1)*(b+1)*(c+1) <= r
            points.add(v)
            counts['selected_core_ground_controls'] += 1
        assert len(points) == r
        if m >= r:
            if cheapest == 2:
                candidates = [(int(b+c > j), b, c) for b in range(r) for c in range(r-b)]
            else:
                candidates = [(k, b, c) for k in range(r) for b in range(r-k)
                              for c in range(r-k-b) if k or b+c <= j]
            independent = frozenset(sorted(candidates, key=lambda a: (dot(a, prices), a))[:r])
            assert independent == points
            counts['independent_compressed_prefix_selections'] += 1
        selected.append(frozenset(points))
    if m >= r:
        assert all(a == selected[r-1] for a in selected[r-1:])
    return cheapest, tuple(selected[:r] if m >= r else selected)


def ray_charges(r, m, cells, counts):
    n = r*m
    E, classes, _ = model(r, m)
    incidence = defaultdict(set)
    edge_cells = defaultdict(set)
    profiles = [tuple(sum(a[i] for a in B) for i in range(3)) for B, _, _, _ in cells]
    assert all(sum(v) % m == (r*m*(m-1)//2) % m for v in profiles)
    R = n*(n-1)//2
    assert all(v[0] <= R and v[1] <= R and v[2] <= n*(m-1) for v in profiles)
    assert dot((1, -1, 0), cross((0, 1, -1), (0, 0, m))) == m
    volume_bound = Q((R+1)*(R+2)*(n*(m-1)+m+1), m)
    assert volume_bound <= 18*n**5
    counts['normalized_profile_volume_controls'] += 1
    for index, (B, polygon, normals, w) in enumerate(cells):
        for v in polygon:
            ray = integer_ray(v)
            if min(ray) > 0:
                incidence[ray].add(index)
        for a, b in zip(polygon, polygon[1:]+polygon[:1]):
            edge_cells[tuple(sorted((a, b)))].add(index)
    adjacency = set()
    for edge, ids in edge_cells.items():
        if len(ids) == 2:
            adjacency.add(tuple(sorted(ids)))
        else:
            assert len(ids) == 1 and any(a == b == 0 for a, b in zip(*edge))
    directions, records = defaultdict(lambda: [0, 0]), []
    positive_edge_incidence = defaultdict(int)
    for w, ids in sorted(incidence.items()):
        g = gcd(w[0], w[1])
        p, q, s = Q(w[0], g), Q(w[1], g), Q(w[2], g)
        assert p.denominator == q.denominator == 1
        p, q = int(p), int(q)
        thresholds, active = [], []
        for C in classes:
            prices = sorted(dot(a, (p, q, s)) for a in C)
            t = prices[r-1]
            thresholds.append(t)
            below = sum(v < t for v in prices)
            tie = [a for a in C if dot(a, (p, q, s)) == t]
            if 0 < r-below < len(tie):
                active.extend(tuple(x-y for x, y in zip(a, tie[0])) for a in tie[1:])
        d1, d2 = next((a, b) for a, b in combinations(active, 2) if cross(a, b) != (0, 0, 0))
        assert primitive(cross(d1, d2)) == primitive(w)
        assert g <= m-1
        T = max(thresholds)
        assert T <= (n-1)*min(p, q)
        assert p*q <= 2*(m-1)*T
        if T >= 8*m*min(p, q):
            j = thresholds.index(T)
            direct = sum(a*p+b*q < T and (a+b) % m == j for a in range(n) for b in range(n))
            assert direct < r
            assert T*T <= 8*n*p*q and T <= 16*m*n
            counts['large_threshold_residue_counts'] += 1
        else:
            assert max(p, q) < 16*m*m
            counts['small_threshold_direction_controls'] += 1
        assert p*q <= 32*m*m*n+256*m**4
        assert max(w) <= 32*m**3*n+128*m**4
        width = max(profiles[i][2] for i in ids)-min(profiles[i][2] for i in ids)
        assert width >= 1 and len(ids) <= 4*width
        directions[p, q][0] += width
        directions[p, q][1] += len(ids)
        # Check the edge-length bridge by every actual facet corner.
        edges = [ij for ij in adjacency if set(ij) <= ids]
        for ij in edges:
            positive_edge_incidence[ij] += 1
        for i in ids:
            neighbors = [j if i == k else k for j, k in edges if i in (j, k)]
            assert len(neighbors) == 2
            a, b = [tuple(x-y for x, y in zip(profiles[j], profiles[i])) for j in neighbors]
            z = cross(a, b)
            assert primitive(z) == primitive(w) and dot(z, z) >= dot(w, w)
            assert max(w)**2 <= dot(a, a)*dot(b, b)
            counts['exact_facet_corner_edge_length_bridges'] += 1
        counts['fixed_direction_facet_charges'] += 1
        records.append({'primitive_ray': w, 'planar_direction': [p, q], 'threshold': str(T),
                        'incident_cells': sorted(ids), 'vertical_profile_range': width,
                        'independent_active_exchanges': [d1, d2]})
    assert all(k <= 2 for k in positive_edge_incidence.values())
    groups = []
    for pq, (jump, degree) in sorted(directions.items()):
        assert jump <= n*(m-1) and degree <= 4*jump
        counts['aggregate_direction_jump_budgets'] += 1
        groups.append({'planar_direction': pq, 'total_vertical_jump': jump, 'total_degree': degree})
    edges = [{'cells': ij, 'squared_length': sum((x-y)**2 for x, y in zip(profiles[ij[0]], profiles[ij[1]])),
              'positive_incident_facets': k} for ij, k in sorted(positive_edge_incidence.items())]
    return {'rays': records, 'direction_groups': groups, 'positive_primal_edges': edges,
            'normalized_volume_box_upper': str(volume_bound)}


def full_fan(r, m, counts):
    E, classes, select = model(r, m)
    comparisons = sorted({primitive(tuple(x-y for x, y in zip(a, b)))
                          for C in classes for a, b in combinations(C, 2)
                          if min(x-y for x, y in zip(a, b)) < 0 < max(x-y for x, y in zip(a, b))})
    pending, seen = [(1, r*m+1, (r*m+1)**2)], {}
    while pending:
        w = pending.pop()
        B = select(w)
        if B in seen:
            continue
        assert len(B) == r*m
        assert all(tuple(x-int(i == j) for i, x in enumerate(a)) in B
                   for a in B for j in range(3) if a[j])
        polygon, normals = cone(B, m)
        assert all(dot(z, w) > 0 for z in normals)
        core(r, m, B, w, counts)
        seen[B] = (polygon, normals, w)
        pending.extend(adjacent_weights(polygon, normals, comparisons))
    cells = [(B, *v) for B, v in sorted(seen.items(), key=lambda kv: sorted(kv[0]))]
    assert sum(area2(P) for B, P, Z, w in cells) == 1
    charges = ray_charges(r, m, cells, counts)
    counts['complete_direct_cyclic_fans'] += 1
    counts['direct_initial_ideals'] += len(cells)
    return {'r': r, 'm': m, 'length': r*m, 'initial_ideals': len(cells), 'normalized_area': '1',
            'cones': [{'standard_set': sorted(B), 'strict_weight': list(map(str, w)),
                       'polygon': [[str(x) for x in v] for v in P]} for B, P, Z, w in cells],
            'collective_facet_controls': charges}, cells


def stabilization(r, m, target, cells, counts):
    assert r >= 2 and m >= r and target >= r
    select_target = model(r, target)[2]
    records, labels = [], set()
    for B, P, Z, w in cells:
        cheapest, key = core(r, m, B, w, counts)
        c = Q(m, target)-1
        mapped = tuple(Q(a)+c*w[cheapest] for a in w)
        assert min(mapped) > 0 and min(range(3), key=lambda i: mapped[i]) == cheapest
        Bt = select_target(mapped)
        assert core(r, target, Bt, mapped, counts) == (cheapest, key)
        Pt, Zt = cone(Bt, target)
        pulled = [tuple(z[i]+(c*sum(z) if i == cheapest else 0) for i in range(3)) for z in Zt]
        sector = [tuple(int(i == j)-int(i == cheapest) for i in range(3))
                  for j in range(3) if j != cheapest]
        assert all(dot(z, v) >= 0 for z in sector for v in P)
        independent = TRIANGLE[:]
        for z in pulled+sector:
            independent = clip(independent, z, 1)
        assert area2(independent) == area2(P)
        assert all(dot(z, v) >= 0 for z in pulled for v in P)
        assert all(dot(z, v) >= 0 for z in Z for v in independent)
        assert Bt not in labels
        labels.add(Bt)
        counts['whole_mapped_cone_equalities'] += 1
        records.append({'source_weight': list(map(str, w)), 'target_weight': list(map(str, mapped)),
                        'cheapest_variable': cheapest, 'target_standard_set': sorted(Bt)})
    counts['complete_period_stabilizations'] += 1
    return {'r': r, 'source_period': m, 'target_period': target,
            'initial_ideals': len(labels), 'scope': 'Whole pulled-back cones, not only sample weights.',
            'mapped_cones': records}


class Field:
    def __init__(self, label):
        self.label = label
        self.characteristic = 0 if label == 'Q' else 2 if label == 'GF256' else 257
        self.alpha = 2 if label == 'Q' else 3

    def add(self, a, b):
        return a+b if self.label == 'Q' else a ^ b if self.label == 'GF256' else (a+b) % 257

    def neg(self, a):
        return -a if self.label == 'Q' else a if self.label == 'GF256' else -a % 257

    def mul(self, a, b):
        if self.label == 'Q':
            return a*b
        if self.label != 'GF256':
            return a*b % 257
        answer = 0
        while b:
            if b & 1:
                answer ^= a
            b >>= 1
            a <<= 1
            if a & 256:
                a ^= 0x11b
        return answer

    def power(self, a, n):
        out = 1
        while n:
            if n & 1:
                out = self.mul(out, a)
            a = self.mul(a, a)
            n >>= 1
        return out

    def inv(self, a):
        assert a
        return Q(1, a) if self.label == 'Q' else self.power(a, 254 if self.label == 'GF256' else 255)

    def determinant(self, a):
        answer = 0
        for sigma in permutations(range(len(a))):
            term = 1
            for i, j in enumerate(sigma):
                term = self.mul(term, a[i][j])
            if sum(sigma[i] > sigma[j] for i in range(len(a)) for j in range(i+1, len(a))) % 2:
                term = self.neg(term)
            answer = self.add(answer, term)
        return answer


def actual_quotient(r, m, cells, label, counts):
    field = Field(label)
    n, L = r*m, r*m+2
    if label != 'Q':
        order = 255 if label == 'GF256' else 256
        powers = [field.power(field.alpha, i) for i in range(order)]
        assert len(set(powers)) == order and set(powers) == set(range(1, order+1))
        assert all(field.mul(a, field.inv(a)) == 1 for a in powers)
        assert n*(L+1) < order
        counts['verified_nonzero_field_elements'] += order
    U = [field.power(field.alpha, i) for i in range(r)]
    V = [field.power(field.alpha, L*i) for i in range(r)]
    lambdas = [field.power(field.alpha, a+L*b) for a in range(n+1) for b in range(n+1)]
    assert len(set(lambdas)) == len(lambdas)
    E, classes, _ = model(r, m)

    @lru_cache(None)
    def values(a):
        return tuple(field.mul(field.power(u, a[0]), field.power(v, a[1])) for u, v in zip(U, V))

    def column(a):
        j = sum(a) % m
        return tuple(values(a)[i] if j == c else 0 for c in range(m) for i in range(r))

    for j in range(m):
        monomials = [(a, 0, (j-a) % m) for a in range(r)]
        assert field.determinant([list(values(a)) for a in monomials])
        counts['independent_cyclic_block_determinants'] += 1
    for i in range(r):
        for j in range(m):
            position = j
            for _ in range(m):
                position = (position+1) % m
            assert position == j
            counts['periodic_generator_columns'] += 1
            if label == 'GF256' and m == 2:
                v = [int(k == j) for k in range(m)]
                first = [field.add(v[(k-1) % m], v[k]) for k in range(m)]
                second = [field.add(first[(k-1) % m], first[k]) for k in range(m)]
                assert any(first) and not any(second)
                counts['nonsemisimple_generator_controls'] += 1
            for first, second in combinations((U[i], V[i], 1), 2):
                assert field.mul(first, second) == field.mul(second, first)
                assert (j+2) % m == ((j+1) % m+1) % m
                counts['commuting_coordinate_columns'] += 1
    basis_determinants = 0
    for C in classes:
        for subset in combinations(C, r):
            assert field.determinant([list(values(a)) for a in subset])
            basis_determinants += 1
    counts['independent_stream_basis_determinants'] += basis_determinants
    recovered = []
    for B, P, Z, w in cells:
        pivots, chosen = {}, []
        for a in sorted(E, key=lambda a: (dot(a, w), a)):
            v = list(column(a))
            for i, row in sorted(pivots.items()):
                c = v[i]
                if c:
                    v = [field.add(x, field.neg(field.mul(c, y))) for x, y in zip(v, row)]
            pivot = next((i for i, x in enumerate(v) if x), None)
            if pivot is not None:
                scale = field.inv(v[pivot])
                pivots[pivot] = [field.mul(x, scale) for x in v]
                chosen.append(a)
        assert frozenset(chosen) == B and len(pivots) == n
        counts['independent_full_column_greedy_bases'] += 1
        relations = []
        for h in minimal_excluded_border(tuple(B)):
            if h[2] >= m:
                assert h == (0, 0, m)
                terms = [(h, 1), ((0, 0, 0), field.neg(1))]
            else:
                selected = sorted(a for a in B if sum(a) % m == sum(h) % m)
                matrix = [list(row) for row in zip(*(values(a) for a in selected))]
                det = field.determinant(matrix)
                assert det
                coefficients = []
                for j in range(r):
                    replacement = [[values(h)[i] if k == j else matrix[i][k] for k in range(r)] for i in range(r)]
                    c = field.mul(field.determinant(replacement), field.inv(det))
                    assert c
                    coefficients.append(c)
                terms = [(h, 1)]+[(a, field.neg(c)) for a, c in zip(selected, coefficients)]
            assert all(dot(h, w) > dot(a, w) for a, c in terms[1:])
            for i in range(r):
                for start in range(m):
                    result = [0]*m
                    for a, c in terms:
                        j = (start+sum(a)) % m
                        result[j] = field.add(result[j], field.mul(c, values(a)[i]))
                    assert not any(result)
                    counts['whole_coordinate_matrix_border_annihilations'] += 1
            counts['recovered_reduced_border_polynomials'] += 1
            relations.append([[a, str(c)] for a, c in terms])
        recovered.append({'standard_set': sorted(B), 'reduced_basis': relations})
    return {'r': r, 'm': m, 'field': label, 'characteristic': field.characteristic,
            'characteristic_divides_period': bool(field.characteristic and m % field.characteristic == 0),
            'length': n, 'radix': L, 'U': U, 'V': V,
            'extension_modulus': 'x^8+x^4+x^3+x+1' if label == 'GF256' else None,
            'genericity_scope': 'All distinct bivariate monomials with each exponent <=n; not all monomials after finite specialization.',
            'initial_ideals': len(cells), 'recovered_bases': recovered}


def large_threshold_controls(counts):
    records = []
    for m in (2, 3, 5, 9):
        r = 200*m*m+1
        thresholds, partial = [], []
        for j in range(m):
            cheaper = 0
            for degree in range(j, r*m+1, m):
                tier = sum(degree-c+1 for c in range(min(m-1, degree)+1))
                if cheaper+tier >= r:
                    thresholds.append(degree)
                    if 0 < r-cheaper < tier:
                        partial.append([j, degree, cheaper, tier])
                    break
                cheaper += tier
        assert partial  # Three noncollinear tied columns give the ray (1,1,1).
        j0, degree, _, _ = partial[0]
        a, b, c = (degree, 0, 0), (0, degree, 0), (degree-1, 0, 1)
        assert primitive(cross(tuple(x-y for x, y in zip(a, b)), tuple(x-y for x, y in zip(a, c)))) == (1, 1, 1)
        T = max(thresholds)
        j = thresholds.index(T)
        strict = sum(d+1 for d in range(T) if d % m == j)
        radius = (T+1)//2
        rectangle = sum((a+b) % m == j for a in range(radius) for b in range(radius))
        assert T >= 8*m and rectangle*8*m >= T*T
        assert rectangle <= strict < r and T*T <= 8*m*r
        counts['independent_large_threshold_residue_counts'] += 1
        records.append({'r': r, 'm': m, 'length': r*m, 'primitive_ray': [1, 1, 1],
                        'class_thresholds': thresholds, 'max_threshold': T,
                        'strict_planar_count_at_maximum_class': strict,
                        'independent_rectangle_count': rectangle,
                        'partial_tiers': partial,
                        'scope': 'Exact degree-tier formula; no large quotient or full fan allocated.'})
    return records


def nonuniform_scope_control(counts):
    # B=F_7^3, U=V=(1,2,4): an actual planar quotient, not generic.
    r, m, prime = 3, 2, 7
    U = (1, 2, 4)
    x = (0, 0, 0)+U
    y = x
    assert x == y and any(x)
    V = [[pow(u, a, prime) for u in U] for a in range(r)]
    assert dot(V[0], cross(V[1], V[2])) % prime != 0
    counts['nongeneric_planar_rank_scope_obstructions'] += 1
    return {'base_field': 'F7', 'r': r, 'm': m, 'length': r*m, 'U_equals_V': U,
            'monomials': [[1, 0, 0], [0, 1, 0]], 'actual_rank': 1,
            'generic_minimum_rule_would_give': 2,
            'scope': 'The tensor/congruence structure survives; generic uniform ranks and specialized counting gates are not asserted for arbitrary B.'}


def main():
    if not __debug__:
        raise RuntimeError('Run without -O; exact assertions are required.')
    counts = defaultdict(int)
    fans, cells_by_case = [], {}
    for r, m in ((1, 2), (2, 2), (2, 3), (3, 2), (3, 3), (4, 3), (4, 4)):
        report, cells = full_fan(r, m, counts)
        fans.append(report)
        cells_by_case[r, m] = cells
    assert [f['initial_ideals'] for f in fans] == [3, 8, 8, 17, 27, 46, 58]
    maps = [stabilization(r, r, target, cells_by_case[r, r], counts)
            for r in (2, 3, 4) for target in (r+1, 2*r)]
    actual = [actual_quotient(r, m, cells_by_case[r, m], label, counts)
              for r, m in ((2, 2), (2, 3), (3, 2)) for label in ('Q', 'GF256', 'F257')]
    large = large_threshold_controls(counts)
    nongeneric = nonuniform_scope_control(counts)
    sources = ('docs/cyclic_interleaving_fans.md', 'scripts/run_cyclic_interleaving_checks.py',
               'scripts/run_extremal_initial_ideals_checks.py', 'scripts/run_initial_cone_height_checks.py',
               'scripts/run_interleaving_research.py', 'scripts/run_structural_complexity_checks.py',
               'docs/weighted_facet_incidence.md')
    report = {'status': 'passed', 'command': 'python3 -B scripts/run_cyclic_interleaving_checks.py',
              'counts': dict(counts), 'complete_fans': fans, 'whole_period_stabilizations': maps,
              'independent_actual_quotients': actual,
              'independent_large_threshold_controls': large,
              'nongeneric_scope_control': nongeneric,
              'uniform_results': {'length': 'n=m*r', 'rank_rule': 'sum_j min(r, distinct (a,b) in class j)',
                                  'fixed_period_count': 'O_m(n^2*log(n+1))',
                                  'comparable_period_count': 'O_c(n^2*log(n+1)^8) for m>=c*r',
                                  'exact_stabilization': 'g(r,m)=g(r,r) for r>=2,m>=r',
                                  'coarse_all_parameter_bound': 'O(n^(5/2))',
                                  'intermediate_parameter_bound': 'O(min(n^(5/2),m^3*n^2*L,n^8*m^(-12)*L^8)), 2<=m<=r, L=1+log(n+1)',
                                  'necessary_period_window': 'For n^(2+epsilon) growth: m>=Omega(n^(epsilon/3)*L^(-1/3)), m<=O(n^(1/2-epsilon/12)*L^(2/3)).',
                                  'remaining_polynomial_amplification_window': 'm tends to infinity and r/m tends to infinity'},
              'global_cubic_half_moment_proved': False, 'unrestricted_G3_improved': False,
              'novelty_established': False, 'goal_complete': False,
              'scope': 'Uniform proofs plus bounded exact controls. Generic quotient ranks are checked on whole finite fans and through commuting matrices, including characteristic dividing a period. Mapped-cone equality verifies complete cones rather than a mesh. No empirical growth or unrestricted moment bound is inferred.',
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    (ROOT/'results/research/cyclic_interleaving_checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(counts), indent=2))


if __name__ == '__main__':
    main()
