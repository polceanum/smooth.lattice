#!/usr/bin/env python3
"""Uniform symbolic and independent finite checks for thin positive cones."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations
from math import comb, gcd
from pathlib import Path
import hashlib
import json

from run_interleaving_research import minimal_excluded_border
from run_primitive_wall_checks import threshold
from run_signed_face_checks import field_basis, field_reduce
from run_structural_complexity_checks import area2, clip, cross, dot, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]


class P(tuple):
    """Small exact integer polynomials in the nonnegative parameter j."""
    def __new__(cls, a=0):
        a = [a] if isinstance(a, int) else list(a)
        while a and not a[-1]:
            a.pop()
        return tuple.__new__(cls, a)

    def __add__(self, other):
        other = P(other)
        return P((self[i] if i < len(self) else 0)
                 +(other[i] if i < len(other) else 0)
                 for i in range(max(len(self), len(other))))

    __radd__ = __add__

    def __neg__(self):
        return P(-x for x in self)

    def __sub__(self, other):
        return self+-P(other)

    def __rsub__(self, other):
        return P(other)+-self

    def __mul__(self, other):
        other = P(other)
        out = [0]*(len(self)+len(other))
        for i, a in enumerate(self):
            for j, b in enumerate(other):
                out[i+j] += a*b
        return P(out)

    __rmul__ = __mul__


def terms(h):
    x = lambda a: (a, 0, 0)
    y = lambda a: (0, a, 0)
    return [
        ((1, 1, 0), [(0, 0, 0)]),
        (x(5*h+4), [(0, 3, 1), x(3*h-4), y(3*h-4)]),
        (y(3*h+1), [x(5*h-1), (0, 0, 1), x(3*h-1), (2, 0, 1),
                       x(3*h+1), y(3*h-1), x(5*h+1)]),
        ((3, 0, 1), [x(3*h+2), y(3*h), x(5*h), y(3*h-2),
                     x(5*h+2), (1, 0, 1), x(3*h)]),
        ((0, 4, 1), [x(3*h-5), x(5*h+3), y(3*h-3)]),
        ((0, 0, 2), [y(2*h+2), x(5*h+2), y(3*h-2)]),
    ]


def critical(h):
    return ((-2, 3*h+1, -1), (-(5*h-1), 0, 1), (5*h+4, -3, -1))


def center(h):
    return 3*h+4, 5*h+6, 15*h*h+17*h-3


def pascal(r):
    return sum((comb(r % 16, j) % 2) << j for j in range(min(r % 16, 9)+1))


def binary_rank(columns):
    pivots = {}
    for v in columns:
        while v:
            i = v.bit_length()-1
            if i not in pivots:
                pivots[i] = v
                break
            v ^= pivots[i]
    return len(pivots)


def symbolic_controls(counts):
    h = P((4, 16))
    L, k = 8*h, 5*h+1
    certificate = []
    for head, tails in terms(h):
        groups, residual, divisions = defaultdict(int), 0, []
        for a in [head]+tails:
            exponent = P(a[0])-P(a[1])-k*P(a[2])
            quotient = (exponent[0] if exponent else 0)//L[0]
            remainder = exponent-quotient*L
            # These coefficient tests prove the same division for every j>=0.
            assert all(c >= 0 for c in remainder)
            assert all(c >= 0 for c in L-remainder-1)
            assert all(c % 16 == 0 for c in remainder[1:])
            groups[remainder] ^= 1
            if quotient % 2:
                residual ^= pascal(remainder[0] if remainder else 0)
            divisions.append({'laurent_exponent_coefficients': list(exponent),
                              'quotient_by_L': quotient,
                              'remainder_coefficients': list(remainder)})
            counts['symbolic_Laurent_divisions_for_all_parameters'] += 1
        assert not any(groups.values()) and residual == 0
        certificate.append({'term_divisions': divisions,
                            'all_ordinary_groups_cancel': True,
                            'ten_dimensional_kernel_residue': residual})
        counts['symbolic_uniform_polynomial_memberships'] += 1
    indices = [5*h+i for i in range(4)]+[3*h+i for i in range(-4, 2)]
    residues = [a[0] % 16 for a in indices]
    assert sorted(residues) == list(range(4, 14))
    assert binary_rank([pascal(a) for a in residues]) == 10
    counts['symbolic_uniform_ten_column_bases'] += 1
    N = tuple(tuple(P(a) for a in q) for q in critical(h))
    rays = [cross(N[i], N[(i+1) % 3]) for i in range(3)]
    assert dot(N[0], cross(N[1], N[2])) == P(2)
    assert dot(rays[0], cross(rays[1], rays[2])) == P(4)
    assert all(all(c >= 0 for c in a) and a[0] > 0 for r in rays for a in r)
    w = center(h)
    assert all(dot(q, w) == P(1) for q in N)
    dominance = []
    for head, tails in terms(h):
        for a in tails:
            q = tuple(P(x)-P(y) for x, y in zip(head, a))
            dots = [dot(q, r) for r in rays]
            assert all(all(c >= 0 for c in v) for v in dots)
            assert any(dots)
            dominance.append([list(v) for v in dots])
            counts['symbolic_head_tail_ray_dominance_certificates'] += 1
    denominators = [dot(q, w) for q in
                    ((2, 3*h+1, 1), (5*h-1, 0, 1), (5*h+4, 3, 1))]
    assert denominators[2] == 30*h*h+64*h+31
    assert all(c >= 0 for c in denominators[2]-denominators[0])
    assert all(c >= 0 for c in denominators[2]-denominators[1])
    return {'parameter': 'h=16*j+4, j>=0',
            'uniform_membership_certificates': certificate,
            'residual_basis_indices_mod_16': residues,
            'critical_determinant': 2, 'cross_ray_determinant': 4,
            'head_tail_ray_dot_polynomials': dominance,
            'coefficient_interpretation': 'Ascending powers of j; nonnegative coefficients prove the asserted inequalities for every j>=0.',
            'center_slacks': [1, 1, 1],
            'maximum_relative_margin_denominator': list(denominators[2])}


def companion_model(h):
    L, k = 8*h, 5*h+1
    n, period = L+10, 2*L
    f = sum(1 << j for j in (0, 2, 8, 10, L, L+2, L+8, L+10))
    powers = [1]
    for e in range(1, period+1):
        v = powers[-1] << 1
        if v >> n:
            v ^= f
        powers.append(v)
    assert powers[-1] == 1

    def ev(a):
        return powers[(a[0]-a[1]-k*a[2]) % period]
    return n, L, k, f, powers, ev


def basis(h):
    return tuple([(i, 0, 0) for i in range(5*h+4)]
                 +[(0, j, 0) for j in range(1, 3*h+1)]
                 +[(i, 0, 1) for i in range(3)]
                 +[(0, j, 1) for j in range(1, 4)])


def polynomial(head, tails):
    out = {head}
    for a in tails:
        out.symmetric_difference_update({a})
    return out


def translate(poly, shift):
    return {tuple(x+y for x, y in zip(a, shift)) for a in poly}


def reduce_polynomial(poly, G, heads, w, counts):
    remaining, remainder = set(poly), set()
    while remaining:
        lead = max(remaining, key=lambda a: (dot(w, a), a))
        i = next((i for i, b in enumerate(heads)
                  if all(x >= y for x, y in zip(lead, b))), None)
        if i is None:
            remaining.remove(lead)
            remainder.add(lead)
        else:
            shift = tuple(x-y for x, y in zip(lead, heads[i]))
            remaining.symmetric_difference_update(translate(G[i], shift))
            counts['independent_polynomial_division_steps'] += 1
    return remainder


def core(q):
    u, v = tuple(max(t, 0) for t in q), tuple(max(-t, 0) for t in q)
    def divisors(a):
        return {(i, j, k) for i in range(a[0]+1)
                for j in range(a[1]+1) for k in range(a[2]+1)}
    return (divisors(u) | divisors(v))-{u, v}


def matrix_apply(M, v):
    out = 0
    while v:
        i = (v & -v).bit_length()-1
        out ^= M[i]
        v &= v-1
    return out


def greedy(n, ev, w):
    pivots, B = {}, []
    for a in sorted(hyperbolic_cross(n, 3), key=lambda a: (dot(w, a), a)):
        v = ev(a)
        while v:
            i = v.bit_length()-1
            if i not in pivots:
                pivots[i] = v
                B.append(a)
                break
            v ^= pivots[i]
        if len(B) == n:
            return frozenset(B)
    raise AssertionError('cyclic span incomplete')


def finite_controls(h, counts):
    n, L, k, f, powers, ev = companion_model(h)
    B = basis(h)
    assert len(B) == n and binary_rank([ev(a) for a in B]) == n
    assert all(tuple(e-(i == j) for j, e in enumerate(a)) in set(B)
               for a in B for i in range(3) if a[i])
    raw = terms(h)
    heads = tuple(a for a, _ in raw)
    assert set(minimal_excluded_border(B)) == set(heads)
    G = [polynomial(a, tails) for a, tails in raw]
    w, N = center(h), critical(h)
    assert all(dot(w, q) == 1 for q in N)
    for a, g in zip(heads, G):
        assert max(g, key=lambda t: (dot(w, t), t)) == a
        assert g-{a} <= set(B)
        total = 0
        for t in g:
            total ^= ev(t)
        assert not total
        counts['independent_companion_relation_replays'] += 1
    for i, j in combinations(range(len(G)), 2):
        top = tuple(max(x, y) for x, y in zip(heads[i], heads[j]))
        left = translate(G[i], tuple(x-y for x, y in zip(top, heads[i])))
        right = translate(G[j], tuple(x-y for x, y in zip(top, heads[j])))
        assert not reduce_polynomial(left ^ right, G, heads, w, counts)
        counts['independent_S_pair_reductions'] += 1
    originals = [{(a, 0, 0) for a in range(n+1) if f >> a & 1},
                 {(1, 1, 0), (0, 0, 0)}, {(0, 0, 1), (0, k, 0)}]
    for g in originals:
        assert not reduce_polynomial(g, G, heads, w, counts)
        counts['independent_original_generator_reductions'] += 1
    # An independent epsilon/Taylor representation for every group power.
    epsilon_columns = []
    for j in range(10):
        a = sum((comb(j, i) % 2) << i for i in range(j+1))
        epsilon_columns.append(a ^ (a << L))
    for E in range(2*L):
        residue, high = E % L, E//L
        v = 1 << residue
        if high:
            v ^= matrix_apply(epsilon_columns, pascal(residue))
        assert v == powers[E]
        counts['independent_epsilon_vs_companion_powers'] += 1
    matrices = [tuple(powers[(j+e) % (2*L)] for j in range(n))
                for e in (1, 2*L-1, (2*L-k) % (2*L))]
    for i, j in combinations(range(3), 2):
        for a in range(n):
            assert matrix_apply(matrices[i], matrices[j][a]) == matrix_apply(matrices[j], matrices[i][a])
            counts['independent_commuting_matrix_columns'] += 1
    rays = [cross(N[i], N[(i+1) % 3]) for i in range(3)]
    primitives = [tuple(t//gcd(*r) for t in r) for r in rays]
    assert all(min(r) > 0 for r in primitives)
    all_normals = [tuple(x-y for x, y in zip(a, b))
                   for a, g in zip(heads, G) for b in g-{a}]
    assert set(N) <= set(all_normals)
    assert all(dot(q, r) >= 0 for q in all_normals for r in rays)
    poly = [tuple(Q(i == j) for i in range(3)) for j in range(3)]
    for q in all_normals:
        poly = clip(poly, q, 1)
    assert len(poly) == 3
    assert set(poly) == {tuple(Q(t, sum(r)) for t in r) for r in rays}
    ratio = [tuple(Q(r[i], r[1]) for i in (0, 2)) for r in rays]
    ratio_area = abs(sum(a[0]*b[1]-a[1]*b[0]
                         for a, b in zip(ratio, ratio[1:]+ratio[:1])))/2
    expected_area = Q(2, 5*(5*h+1)*(5*h+6))
    assert ratio_area == expected_area
    assert min(a[0] for a in ratio)*min(a[1] for a in ratio) > h
    assert max(a[0] for a in ratio)*max(a[1] for a in ratio) < 3*h
    margin_lower = min(Q(dot(q, w), dot(tuple(abs(t) for t in q), w)) for q in N)
    assert margin_lower == Q(1, 30*h*h+64*h+31)
    cores = [core(q) for q in N]
    union = set().union(*cores)
    assert union == set(B)-{(2, 0, 1), (0, 3, 1)}
    assert len(union) == n-2
    assert [threshold(q) for q in N] == [3*h+6, 5*h, 5*h+11]
    box_rank = binary_rank([powers[a % (2*L)] for a in range(-(2*h-1), 2*h)])
    assert box_rank == 4*h-1 < n
    counts['independent_superlinear_nonsaturating_boxes'] += 1
    neighbors = []
    if h <= 20:
        assert greedy(n, ev, w) == frozenset(B)
        for i, q in enumerate(N):
            wall = tuple(sum(r[a] for r in rays if dot(q, r) == 0) for a in range(3))
            outside = tuple(2*n**4*t-s for t, s in zip(wall, w))
            assert min(outside) > 0 and dot(q, outside) < 0
            neighbor = greedy(n, ev, outside)
            assert neighbor != frozenset(B)
            delta = tuple(sum(a[j] for a in neighbor)-sum(a[j] for a in B) for j in range(3))
            factor = gcd(*delta)
            direction = tuple(t//factor for t in delta)
            assert direction == q or direction == tuple(-t for t in q)
            assert cores[i] <= neighbor
            if i == 1:
                assert not union <= neighbor
            neighbors.append({'wall_normal': q, 'exchange_multiple': factor,
                              'weight': outside, 'basis': sorted(neighbor),
                              'joint_core_retained': union <= neighbor})
            counts['independent_actual_wall_neighbors'] += 1
    counts['independent_marked_complete_cones'] += 1
    return {'h': h, 'length': n, 'characteristic': 2,
            'basis': B, 'reduced_basis': [{'head': a, 'tails': sorted(g-{a})}
                                         for a, g in zip(heads, G)],
            'primitive_normals': N, 'primitive_rays': primitives,
            'positive_simplex_polygon': [list(map(str, a)) for a in poly],
            'minimum_integer_weight': w, 'minimum_integer_height': max(w),
            'relative_margin_bounds': [str(margin_lower), str(Q(1, (3*h+1)*(5*h-1)))],
            'ratio_triangle_area': str(ratio_area),
            'logarithmic_area_bounds': [str(ratio_area/(3*h)), str(ratio_area/h)],
            'wall_costs': [threshold(q) for q in N],
            'joint_core_rank': len(union), 'actual_neighbors': neighbors,
            'nonsaturating_box': {'side_counts': [2*h, 2*h, 1],
                                  'cardinality': 4*h*h, 'rank': box_rank}}


def all_field_ranks(columns, p):
    """DFS maintains independent pivots, rather than using the proposed formula."""
    ranks = [0]*(1 << len(columns))
    def visit(start, mask, pivots):
        ranks[mask] = len(pivots)
        for i in range(start, len(columns)):
            v = list(columns[i])
            next_pivots = dict(pivots)
            for j, row in sorted(pivots.items()):
                t = v[j]
                if t:
                    v = [(a-t*b) % p for a, b in zip(v, row)]
            if any(v):
                j = next(j for j, t in enumerate(v) if t)
                inv = pow(v[j], -1, p)
                next_pivots[j] = tuple(t*inv % p for t in v)
            visit(i+1, mask | (1 << i), next_pivots)
    visit(0, 0, {})
    return ranks


def fourier_controls(counts):
    E = (0, 2, 4, 6, 8, 10, 12, 14, 1, 9)
    # 3 has order sixteen in F_17; all calculations are exact.
    assert pow(3, 16, 17) == 1 and pow(3, 8, 17) != 1
    binary = [tuple((pascal(a) >> j) & 1 for j in range(10)) for a in range(16)]
    fourier = [tuple(pow(3, e*a, 17) for e in E) for a in range(16)]
    left, right = all_field_ranks(binary, 2), all_field_ranks(fourier, 17)
    assert left == right
    for S, rk in enumerate(left):
        low, high = S & 255, S >> 8
        doubles = low & high
        formula = (low | high).bit_count()+bool(doubles & 0x55)+bool(doubles & 0xaa)
        assert rk == formula
        counts['independent_all_subset_Pascal_Fourier_rank_equalities'] += 1
    fixtures = []
    for h, p in ((4, 257), (20, 641)):
        n, L, k = 8*h+10, 8*h, 5*h+1
        M = 2*L
        factors = [d for d in range(2, M+1) if M % d == 0
                   and all(d % q for q in range(2, int(d**0.5)+1))]
        omega = next(pow(a, (p-1)//M, p) for a in range(2, p)
                     if all(pow(pow(a, (p-1)//M, p), M//d, p) != 1 for d in factors))
        assert pow(omega, M, p) == 1
        nodes = [pow(omega, 2*a, p) for a in range(L)]
        nodes += [pow(omega, 1+h*e, p) for e in E]
        assert len(set(nodes)) == n
        def ev(a):
            exponent = (a[0]-a[1]-k*a[2]) % M
            return {i: pow(t, exponent, p) for i, t in enumerate(nodes)}
        B = basis(h)
        selected, pivots = field_basis(B, ev, p)
        assert len(selected) == n
        constraints, lifted = [], []
        for head, tails in terms(h):
            rem, co = field_reduce(ev(head), pivots, p)
            assert not rem
            expected = polynomial(head, tails)-{head}
            assert set(co) == expected
            lifted.append({'head': head,
                           'tails': [{'exponents': a, 'coefficient': c}
                                     for a, c in sorted(co.items())]})
            constraints.extend(tuple(x-y for x, y in zip(head, a)) for a in co)
            counts['independent_radical_fundamental_circuit_supports'] += 1
        N = critical(h)
        rays = [cross(N[i], N[(i+1) % 3]) for i in range(3)]
        assert set(N) <= set(constraints)
        assert all(dot(q, r) >= 0 for q in constraints for r in rays)
        assert all(dot(q, center(h)) > 0 for q in constraints)
        fixtures.append({'h': h, 'length': n, 'characteristic': p,
                         'primitive_root': omega, 'root_order': M,
                         'nodes': nodes, 'marked_reduced_basis': lifted,
                         'full_marked_cone_agrees': True})
        counts['independent_radical_point_cone_equalities'] += 1
    return {'residual_frequency_set': E, 'exhaustive_subsets': len(left),
            'binary_characteristic': 2, 'Fourier_characteristic': 17,
            'formula': 'occupied pairs modulo8 plus parities among double pairs',
            'direct_radical_point_fixtures': fixtures,
            'scope': 'Finite exact rank controls accompany the uniform Fourier proof over characteristic zero. Reduction modulo one odd prime alone is not used to infer characteristic-zero ranks.'}


def main():
    counts = defaultdict(int)
    symbolic = symbolic_controls(counts)
    finite = [finite_controls(h, counts) for h in (4, 20, 36, 68, 132)]
    fourier = fourier_controls(counts)
    files = ['docs/thin_initial_cones.md', 'scripts/run_thin_initial_cone_checks.py',
             'scripts/run_interleaving_research.py', 'scripts/run_primitive_wall_checks.py',
             'scripts/run_structural_complexity_checks.py', 'scripts/run_signed_face_checks.py',
             'docs/initial_cone_stability.md', 'docs/prime_digit_lifts.md']
    report = {'schema_version': 1, 'all_checks_passed': True,
              'command': 'python3 -B scripts/run_thin_initial_cone_checks.py',
              'scope': 'Symbolic certificates establish the stated family for every nonnegative parameter over every characteristic-two field. A uniform Fourier proof transfers the whole labelled matroid to radical characteristic-zero points of the same length. Independent finite actual-algebra checks cover one full marked cone per quotient, not the entire fan. No unrestricted count improvement, lower separation or novelty is inferred.',
              'prior_general_obstruction': 'docs/initial_cone_stability.md already proves sharp unrestricted worst margins Theta(n^-3) and examples of logarithmic area Theta(n^-4). The present exponents are weaker; the additional scope is radical characteristic-zero points and simultaneous intermediate-cost cores. The general margin premise must not be reopened.',
              'uniform_results': {'relative_margin': 'Theta(n^(-2))',
                                  'logarithmic_cone_area': 'Theta(n^(-3))',
                                  'minimum_integer_height': '15*h^2+17*h-3',
                                  'wall_costs': ['3*h+6', '5*h', '5*h+11'],
                                  'union_of_three_compulsory_cores': 'B minus {x^2*z,y^3*z}, rank n-2',
                                  'radical_point_shortcut_also_rejected': 'Uniform Omega(1/n) relative margin or Omega(1/n^2) logarithmic area even when restricted to radical characteristic-zero point ideals.',
                                  'radical_characteristic_zero_transfer': 'Same length, labelled Laurent monomial matroid and whole positive fan via the explicit L-plus-ten Fourier row set.',
                                  'collective_wall_budget_refuted': False,
                                  'unrestricted_G3_improved': False,
                                  'polynomial_lower_separation_proved': False,
                                  'novelty_established': False},
              'counts': dict(counts), 'symbolic_certificates': symbolic,
              'finite_actual_quotients': finite,
              'Fourier_transfer_controls': fourier,
              'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files}}
    target = ROOT/'results/research/thin_initial_cone_checks.json'
    target.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'all_checks_passed': True, 'counts': dict(counts),
                      'lengths': [f['length'] for f in finite],
                      'report': str(target.relative_to(ROOT))}, indent=2))


if __name__ == '__main__':
    main()
