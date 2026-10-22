#!/usr/bin/env python3
"""Actual branch envelopes, one-ray contractions and complete quotient fans."""
import ast
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations
from math import comb, factorial, isqrt, prod
from pathlib import Path
import hashlib
import json

from run_forced_socle_fan_checks import model as forced_model
from run_ideal_step_incidence_checks import actual_edge_length, fan_budgets
from run_projection_tag_checks import compare_fans, fan, polygon_key
from run_signed_face_checks import field_basis, field_reduce, normals as border_normals
from run_structural_complexity_checks import area2, clip, determinant, hyperbolic_cross

ROOT = Path(__file__).resolve().parents[1]
TRIANGLE = [tuple(Q(i == j) for i in range(3)) for j in range(3)]


def clean(v, p):
    return {i: c % p if p else Q(c) for i, c in v.items()
            if (c % p if p else c)}


def linear(terms, p):
    out = defaultdict(int)
    for c, v in terms:
        for i, d in v.items():
            out[i] += c*d
    return clean(out, p)


def apply(M, v, p):
    return linear(((c, M[j]) for j, c in v.items()), p)


class Frame:
    """Sparse coordinate solver using the existing, separately checked code."""
    def __init__(self, vectors, p):
        self.vectors = tuple(vectors)
        self.p = p
        self.selected, self.pivots = field_basis(
            range(len(vectors)), lambda j: vectors[j], p)
        self.basis = tuple(vectors[j] for j in self.selected)
        self.index = {j: i for i, j in enumerate(self.selected)}

    def coords(self, v):
        rem, co = field_reduce(v, self.pivots, self.p)
        assert not rem, 'Vector is outside the claimed module.'
        return clean({self.index[j]: c for j, c in co.items()}, self.p)


def dense_kernel(columns, height, p):
    """Independent dense Gauss-Jordan elimination over Q or a prime field."""
    width = len(columns)
    A = [[columns[j].get(i, 0) for j in range(width)] for i in range(height)]
    A = [[c % p if p else Q(c) for c in row] for row in A]
    pivots, k = [], 0
    for j in range(width):
        i = next((i for i in range(k, height) if A[i][j]), None)
        if i is None:
            continue
        A[k], A[i] = A[i], A[k]
        inv = pow(A[k][j], -1, p) if p else 1/A[k][j]
        A[k] = [c*inv % p if p else c*inv for c in A[k]]
        for i in range(height):
            if i == k or not A[i][j]:
                continue
            c = A[i][j]
            A[i] = [(u-c*v) % p if p else u-c*v
                    for u, v in zip(A[i], A[k])]
        pivots.append(j)
        k += 1
    result = []
    for j in sorted(set(range(width))-set(pivots)):
        v = {j: 1}
        v.update({b: -A[i][j] for i, b in enumerate(pivots) if A[i][j]})
        v = clean(v, p)
        assert not apply(columns, v, p)
        result.append(v)
    return tuple(result), len(pivots)


def dense_rank(columns, height, p):
    return dense_kernel(columns, height, p)[1]


def polynomial(X, v, coefficients, p):
    out = {}
    for c in reversed(coefficients):
        out = linear(((1, apply(X, out, p)), (c, v)), p)
    return out


def module_generators(X, vectors, factors, p, counts):
    """Invariant-factor count using complete, explicitly known primary factors."""
    frame = Frame(vectors, p)
    assert len(frame.basis) == len(vectors)
    if not vectors:
        return 0
    for v in vectors:
        frame.coords(apply(X, v, p))
        # Every supplied factor is irreducible. Check the full support:
        # (product factors)^dim kills the entire module.
        killed = v
        for _ in vectors:
            for f in factors:
                assert len(f) == 2 or (p == 2 and f == (1, 1, 1))
                killed = polynomial(X, killed, f, p)
        assert not killed
    ranks = []
    for f in factors:
        images = [frame.coords(polynomial(X, v, f, p)) for v in vectors]
        rank = dense_rank(images, len(vectors), p)
        assert rank == len(Frame(images, p).basis)
        deficit = len(vectors)-rank
        assert deficit % (len(f)-1) == 0
        ranks.append(deficit//(len(f)-1))
        counts['primary_module_rank_comparisons'] += 1
    return max(ranks)


def cyclic_span(M, unit, p):
    columns = [unit]
    old = 0
    while True:
        frame = Frame(columns, p)
        if len(frame.basis) == old:
            return old
        old = len(frame.basis)
        columns = list(frame.basis)
        columns.extend(apply(X, v, p) for X in M for v in frame.basis)


def recover(name, n, M, unit, factors, p, counts, gorenstein=False, reduced=False):
    X, Yop, Zop = M
    assert cyclic_span(M, unit, p) == n
    for U, V in combinations(M, 2):
        for i in range(n):
            assert apply(U, V[i], p) == apply(V, U[i], p)
            counts['actual_commuting_columns'] += 1
    assert all(not apply(Yop, v, p) for v in Zop)

    FY, FZ = Frame(Yop, p), Frame(Zop, p)
    ybasis, zbasis = FY.basis, FZ.basis
    my, mz = len(ybasis), len(zbasis)
    joined = ybasis+zbasis
    joined_rank = dense_rank(joined, n, p)
    assert joined_rank == len(Frame(joined, p).basis)
    h, d = my+mz-joined_rank, n-joined_rank
    intersection_kernel, _ = dense_kernel(
        tuple(ybasis)+tuple(linear(((-1, v),), p) for v in zbasis), n, p)
    W = tuple(apply(ybasis, {j: c for j, c in v.items() if j < my}, p)
              for v in intersection_kernel)
    assert len(W) == h == len(Frame(W, p).basis)
    assert all(not apply(Yop, v, p) and not apply(Zop, v, p) for v in W)
    assert h <= d and my+mz <= n
    for op, frame in ((Yop, FY), (Zop, FZ)):
        image = [frame.coords(apply(op, v, p)) for v in frame.basis]
        nullity = len(frame.basis)-dense_rank(image, len(frame.basis), p)
        assert h <= nullity <= d
        counts['actual_branch_surjection_dimensions'] += 1

    core, power = [], unit
    while len(Frame(core+[power], p).basis) > len(core):
        core.append(power)
        power = apply(X, power, p)
        assert len(core) <= n
    FC = Frame(core, p)
    a = len(core)
    co = FC.coords(power)
    f = tuple(-co.get(j, 0) for j in range(a))+(1,)
    assert not polynomial(X, unit, f, p)
    assert all(not polynomial(X, v, f, p) for v in ybasis+zbasis)
    frames, offsets = (FC, FY, FZ), (0, a, a+my)
    pi = tuple(core)+ybasis+zbasis
    N = len(pi)
    assert N == a+my+mz == a+n+h-d and N <= 2*n

    def lift(v, group):
        return {offsets[group]+i: c for i, c in frames[group].coords(v).items()}

    HX = tuple(lift(apply(X, v, p), g)
               for g, frame in enumerate(frames) for v in frame.basis)
    HY = tuple(lift(apply(Yop, v, p), 1) for v in FC.basis+FY.basis)
    HY += tuple({} for _ in FZ.basis)
    HZ = tuple(lift(apply(Zop, v, p), 2) for v in FC.basis)
    HZ += tuple({} for _ in FY.basis)
    HZ += tuple(lift(apply(Zop, v, p), 2) for v in FZ.basis)
    HM = HX, HY, HZ
    assert cyclic_span(HM, {0: 1}, p) == N
    for U, V in combinations(HM, 2):
        for i in range(N):
            assert apply(U, V[i], p) == apply(V, U[i], p)
            counts['reconstructed_envelope_commuting_columns'] += 1
    for original, lifted in zip(M, HM):
        for i in range(N):
            assert apply(pi, lifted[i], p) == apply(original, pi[i], p)
            counts['independent_quotient_map_coordinate_columns'] += 1
    T, rank_pi = dense_kernel(pi, n, p)
    assert rank_pi == n and len(T) == a+h-d <= a
    FT = Frame(T, p)
    for v in T:
        FT.coords(apply(HX, v, p))
        assert not apply(HY, v, p) and not apply(HZ, v, p)
        counts['actual_kernel_annihilation_and_invariance'] += 1
    for start, frame in zip(offsets, frames):
        outside = tuple({i: c for i, c in v.items()
                         if not start <= i < start+len(frame.basis)} for v in T)
        assert dense_rank(outside, N, p) == len(T)
        counts['kernel_disjoint_summands'] += 1
    projection = [{i: c for i, c in v.items() if i < a} for v in T]
    FS = Frame(projection, p)
    assert len(FS.basis) == a-d and len(T)-len(FS.basis) == h
    core_X = tuple(FC.coords(apply(X, v, p)) for v in core)
    rw = module_generators(X, W, factors, p, counts)
    rt = module_generators(HX, T, factors, p, counts)
    rs = module_generators(core_X, FS.basis, factors, p, counts)
    assert rw*(rw+1) <= n and rt <= rw+1 and rs <= 1
    assert my >= rw*(rw+1)//2 and mz >= rw*(rw+1)//2
    duality = None
    if gorenstein:
        assert rw <= 1 and rt <= 2
        stacked_A = tuple({axis*n+i: c for axis, U in enumerate((Yop, Zop))
                           for i, c in U[j].items()} for j in range(n))
        ann_A, _ = dense_kernel(stacked_A, 2*n, p)
        assert len(ann_A) == d
        generators_A = module_generators(X, ann_A, factors, p, counts)
        assert generators_A <= 1
        stacked_H = tuple({axis*N+i: c for axis, U in enumerate((HY, HZ))
                           for i, c in U[j].items()} for j in range(N))
        ann_H, _ = dense_kernel(stacked_H, 2*N, p)
        generators_H = module_generators(HX, ann_H, factors, p, counts)
        assert generators_H <= 3
        duality = {'actual_joint_annihilator_length': len(ann_A),
                   'actual_joint_annihilator_generators': generators_A,
                   'envelope_joint_annihilator_length': len(ann_H),
                   'envelope_joint_annihilator_generators': generators_H}
        counts['gorenstein_full_annihilator_duality_controls'] += 1
    if reduced:
        assert h == rw == 0 and rt <= 1
    counts['complete_coordinate_preserving_reconstructions'] += 1
    return {'fixture': name, 'characteristic': p, 'length': n,
            'x_core_length': a, 'minimal_x_polynomial_low_to_high': f,
            'Y_length': my, 'Z_length': mz, 'intersection_length': h,
            'axis_quotient_length': d, 'envelope_length': N,
            'kernel_length': len(T), 'W_module_generators': rw,
            'T_module_generators': rt, 'S_module_generators': rs,
            'known_primary_factors_low_to_high': factors,
            'gorenstein_fixture': gorenstein, 'reduced_fixture': reduced,
            'gorenstein_duality': duality}


def sharp_model(r, p):
    D = tuple([(i, 0, 0) for i in range(r+1)] +
              [(i, j, 0) for j in range(1, r+1) for i in range(r+1-j)] +
              [(i, 0, j) for j in range(1, r+1) for i in range(r+1-j)])
    index = {v: i for i, v in enumerate(D)}
    relations = [{index[i, r-i, 0]: 1, index[i, 0, r-i]: -1}
                 for i in range(r)]
    relations.append({index[r, 0, 0]: 1, index[0, r, 0]: -1})
    F = Frame([clean(v, p) for v in relations], p)
    assert len(F.basis) == r+1
    free = tuple(i for i in range(len(D)) if i not in F.pivots)
    compact = {i: j for j, i in enumerate(free)}

    @lru_cache(None)
    def ev(alpha):
        if alpha not in index:
            return {}
        rem, _ = field_reduce({index[alpha]: 1}, F.pivots, p)
        return clean({compact[i]: c for i, c in rem.items()}, p)

    reference = tuple(D[i] for i in free)
    assert all(ev(v) == {i: 1} for i, v in enumerate(reference))
    M = tuple(tuple(ev(tuple(a+int(i == j) for j, a in enumerate(v)))
                    for v in reference) for i in range(3))
    for relation in relations:
        assert not linear(((c, ev(D[i])) for i, c in relation.items()), p)
        for axis in range(3):
            assert all(tuple(v+int(axis == j) for j, v in enumerate(D[i])) not in index
                       for i in relation)
    return len(reference), M, ev((0, 0, 0)), ev


def uniserial(m, p):
    n = 2*m
    def ev(v):
        exponent = v[0]+m*(v[1]+v[2])
        return {exponent: 1} if exponent < n else {}
    M = tuple(tuple(ev((i, int(axis == 1), int(axis == 2)))
                    if axis else ev((i+1, 0, 0)) for i in range(n))
              for axis in range(3))
    return n, M, {0: 1}


def axis_padding(data, t, p):
    n, M, unit = data
    axis = (tuple(clean({i: 1, i+1: 1} if i+1 < t else {i: 1}, p)
                  for i in range(t)),
            tuple({} for _ in range(t)), tuple({} for _ in range(t)))
    joined = tuple(tuple(U)+tuple({n+i: c for i, c in v.items()} for v in V)
                   for U, V in zip(M, axis))
    return n+t, joined, {**unit, n: 1}


def partitions(n, cap=None):
    if not n:
        yield ()
        return
    for a in range(min(n, n if cap is None else cap), 0, -1):
        for tail in partitions(n-a, a):
            yield (a,)+tail


def staircase_controls(counts):
    for n in range(1, 17):
        for shape in partitions(n):
            D = {(i, j) for j, a in enumerate(shape) for i in range(a)}
            tips = [v for v in D if (v[0]+1, v[1]) not in D
                    and (v[0], v[1]+1) not in D]
            assert len(tips)*(len(tips)+1)//2 <= n
            counts['independent_planar_tip_colength_controls'] += 1


def contracted_planar_intervals(E, evaluate, rank, p, counts):
    """All positive profiles of a vector matroid, certified by every circuit.

    No algebra structure or border reduction is assumed for the contraction.
    Coordinates on E are at most rank-1; distinct comparison fractions are
    separated by more than 1/(rank+1)^3.
    """
    epsilon = Q(1, (rank+1)**3)
    sample, previous_high, seen = epsilon, Q(0), set()
    while True:
        B, pivots = field_basis(sorted(E, key=lambda a: (a[0]*sample+a[1], a)),
                               evaluate, p)
        assert len(B) == rank and frozenset(B) not in seen
        seen.add(frozenset(B))
        low, high = Q(0), None
        for e in E:
            if e in B:
                continue
            rem, co = field_reduce(evaluate(e), pivots, p)
            assert not rem
            for b in co:
                dx, dy = e[0]-b[0], e[1]-b[1]
                assert dx*sample+dy > 0
                if dx > 0:
                    low = max(low, Q(-dy, dx))
                elif dx < 0:
                    value = Q(-dy, dx)
                    high = value if high is None else min(high, value)
                else:
                    assert dy > 0
                counts['contracted_matroid_fundamental_comparisons'] += 1
        assert low == previous_high and sample > low
        assert high is None or sample < high
        yield frozenset(B), low, high
        counts['complete_contracted_planar_intervals'] += 1
        if high is None:
            break
        previous_high, sample = high, high+epsilon


def one_ray_compile(n, ev, mode, p, counts, actual_zero_products=True):
    assert mode in ('parallel', 'axis')
    if actual_zero_products:
        assert not ev((0, 1, 1))
        assert not ev((0, 0, 2) if mode == 'parallel' else (1, 0, 1))
    plane = tuple(hyperbolic_cross(n, 2))
    ell = len(field_basis(plane, lambda a: ev((*a, 0)), p)[0])
    ray = tuple((i, 0, 1) if mode == 'parallel' else (0, 0, i+1)
                for i in range(n))
    independent, _ = field_basis(ray, ev, p)
    m = len(independent)
    assert independent == ray[:m]
    eta = ell+m-n
    assert 0 <= eta <= m
    if actual_zero_products and mode == 'axis':
        assert eta <= 1
    elif actual_zero_products:
        assert 2*m <= n
    predicted, ledger, total_bound = {}, [], 0
    for k in range(n-ell, min(m, n-1)+1):
        N = n-k
        C = ray[:k]
        FC = Frame([ev(c) for c in C], p)
        assert len(FC.basis) == k
        E = tuple(hyperbolic_cross(N, 2))
        @lru_cache(None)
        def contracted(e):
            return field_reduce(ev((*e, 0)), FC.pivots, p)[0]
        remaining = len(field_basis(E, contracted, p)[0])
        ceiling = 1+4*N*(N-1).bit_length()
        total_bound += ceiling
        if remaining < N:
            ledger.append({'ray_prefix': k, 'planar_rank': N,
                           'actual_ground_rank': remaining, 'intervals': 0,
                           'candidate_ceiling': ceiling})
            continue
        assert remaining == N
        intervals = list(contracted_planar_intervals(E, contracted, N, p, counts))
        assert len(intervals) <= ceiling
        accepted, not_downsets = 0, 0
        for P, low, high in intervals:
            B = frozenset(C+tuple((*a, 0) for a in P))
            assert len(B) == n
            if any(tuple(a[j]-int(j == axis) for j in range(3)) not in B
                   for a in B for axis in range(3) if a[axis]):
                not_downsets += 1
                continue
            normals, _ = border_normals(tuple(sorted(B)), ev, p)
            poly = list(TRIANGLE)
            for q in normals:
                poly = clip(poly, q, 1)
            if not area2(poly):
                continue
            assert B not in predicted
            predicted[B] = poly
            accepted += 1
            # The full cone must lie in the contracted planar interval.
            assert all(v[0] >= low*v[1] for v in poly)
            assert high is None or all(v[0] <= high*v[1] for v in poly)
            counts['one_ray_original_border_cones'] += 1
        ledger.append({'ray_prefix': k, 'planar_rank': N,
                       'actual_ground_rank': remaining,
                       'intervals': len(intervals), 'coherent_full_cones': accepted,
                       'non_downset_candidates': not_downsets,
                       'candidate_ceiling': ceiling})
    covered = sum(area2(poly) for poly in predicted.values())
    assert covered == 1 if actual_zero_products else covered <= 1
    assert len(predicted) <= total_bound
    return predicted, {'plane_length': ell, 'outside_ray_rank': m,
                       'plane_ray_overlap': eta, 'candidate_ceiling': total_bound,
                       'prefix_contractions': ledger,
                       'actual_zero_product_hypotheses': actual_zero_products,
                       'normalized_full_cone_area_covered': covered}


def pascal_ray_model(m, p):
    @lru_cache(None)
    def ev(alpha):
        i, j, k = alpha
        if j+k >= 2:
            return {}
        if not j and not k:
            return {i: 1} if i < 2*m else {}
        return clean({m+i+s: comb(m if j else 2*m, s)
                      for s in range(max(0, m-i))}, p)
    return 2*m, ev


def pascal_ray_predicted(m, collapse=False):
    core = {(i, 0, 0) for i in range(m)}
    starts = ((m, 0, 0), (0, 1, 0), (0, 0, 1))
    compositions = ([(m, 0, 0), (0, m, 0), (0, 0, m)] if collapse else
                    [(a, b, m-a-b) for a in range(m+1) for b in range(m-a+1)])
    out = {}
    for sizes in compositions:
        groups = [tuple((v[0]+i, v[1], v[2]) for i in range(size))
                  for v, size in zip(starts, sizes)]
        B = frozenset(core | set().union(*map(set, groups)))
        poly = list(TRIANGLE)
        if collapse:
            chosen = sizes.index(m)
            normals = [tuple(v[i]-starts[chosen][i] for i in range(3))
                       for j, v in enumerate(starts) if j != chosen]
        else:
            normals = []
            for j, size in enumerate(sizes):
                if not size:
                    continue
                last = groups[j][-1]
                for v, other_size in zip(starts, sizes):
                    if other_size == m:
                        continue
                    first = (v[0]+other_size, v[1], v[2])
                    normals.append(tuple(a-b for a, b in zip(first, last)))
        for q in normals:
            poly = clip(poly, q, 1)
        assert area2(poly) > 0
        out[B] = poly
    assert sum(area2(poly) for poly in out.values()) == 1
    return out


def point_evaluate(points, p):
    @lru_cache(None)
    def ev(alpha):
        return clean({i: prod(v**e for v, e in zip(point, alpha))
                      for i, point in enumerate(points)}, p)
    return ev


def axis_connected_model(size, ray_length, p, box=False):
    if box:
        a, b = size
        plane = tuple((i, j, 0) for j in range(b) for i in range(a))
        tip = a-1, b-1, 0
    else:
        plane = tuple((i, 0, 0) for i in range(size))
        tip = size-1, 0, 0
    reference = plane+tuple((0, 0, j) for j in range(1, ray_length))
    index = {v: i for i, v in enumerate(reference)}
    @lru_cache(None)
    def ev(alpha):
        i, j, k = alpha
        if k:
            if i or j or k > ray_length:
                return {}
            return {index[tip if k == ray_length else alpha]: 1}
        if box or not j:
            return {index[alpha]: 1} if alpha in index else {}
        return clean({index[i+j+s, 0, 0]: (-1)**(j+s)*comb(j+s-1, s)
                      for s in range(max(0, size-i-j))}, p)
    return len(reference), ev


def actual_evaluation_controls(n, ev, p, counts):
    E = tuple(hyperbolic_cross(n, 3))
    selected, _ = field_basis(E, ev, p)
    assert len(selected) == n
    F = Frame([ev(a) for a in selected], p)
    M = tuple(tuple(F.coords(ev(tuple(a+int(i == j) for j, a in enumerate(v))))
                    for v in selected) for i in range(3))
    for X, Y in combinations(M, 2):
        for j in range(n):
            assert apply(X, Y[j], p) == apply(Y, X[j], p)
            counts['one_ray_actual_commuting_columns'] += 1
    for a in E:
        for i in range(3):
            successor = tuple(b+int(i == j) for j, b in enumerate(a))
            assert apply(M[i], F.coords(ev(a)), p) == F.coords(ev(successor))
            counts['one_ray_actual_successor_columns'] += 1


def rotate_model(n, ev, p, shift=(0, 0, 0)):
    """One fixed integral invertible linear change, retaining actual labels."""
    E = tuple(hyperbolic_cross(n, 3))
    selected, _ = field_basis(E, ev, p)
    F = Frame([ev(v) for v in selected], p)
    old = tuple(tuple(F.coords(ev(tuple(a+int(i == j) for j, a in enumerate(v))))
                      for v in selected) for i in range(3))
    T = ((0, -1, 1), (1, 1, -1), (-1, 0, 1))
    M = tuple(tuple(linear(
        [(T[i][a], old[a][j]) for a in range(3)]+[(shift[i], {j: 1})], p)
        for j in range(n)) for i in range(3))
    unit = F.coords(ev((0, 0, 0)))
    @lru_cache(None)
    def rotated(alpha):
        if not any(alpha):
            return unit
        axis = next(i for i, a in enumerate(alpha) if a)
        predecessor = tuple(a-int(i == axis) for i, a in enumerate(alpha))
        return apply(M[axis], rotated(predecessor), p)
    return rotated


def common_factor_carriers(l0, l1, l2, p):
    """Initial ideals of l0*(l1,l2), enumerating six variable orders."""
    carriers = {}
    for order in permutations(range(3)):
        priority = tuple(reversed(order))
        forms = [clean({priority.index(i): c for i, c in enumerate(v)}, p)
                 for v in (l1, l2)]
        selected, pivots = field_basis((0, 1), lambda j: forms[j], p)
        assert len(selected) == 2
        heads = {priority[j] for j in pivots}
        z = next(i for i in priority if (l0[i] % p if p else l0[i]))
        if z in heads:
            y = next(iter(heads-{z}))
            x = next(i for i in range(3) if i not in (y, z))
            mode, perm = 'parallel', (x, y, z)
        else:
            mode, perm = 'axis', tuple(i for i in range(3) if i != z)+(z,)
        carriers[mode, perm] = tuple(tuple(int(i == z)+int(i == j) for i in range(3))
                                    for j in sorted(heads))
    assert 1 <= len(carriers) <= 6
    return carriers


def common_factor_controls(counts):
    rows = []
    def check(name, n, original, l0, l1, l2, p, shift=(0, 0, 0)):
        ev = rotate_model(n, original, p, shift)
        actual_evaluation_controls(n, ev, p, counts)
        def affine_terms(l):
            constant = -sum(c*s for c, s in zip(l, shift))
            return [((0, 0, 0), constant)]+[
                (tuple(int(i == j) for i in range(3)), c) for j, c in enumerate(l)]
        for l in (l1, l2):
            terms = []
            for u, c in affine_terms(l0):
                for v, d in affine_terms(l):
                    a = tuple(i+j for i, j in zip(u, v))
                    terms.append((c*d, ev(a)))
            assert not linear(terms, p)
            counts['actual_common_factor_quadratic_relations'] += 1
        carriers = common_factor_carriers(l0, l1, l2, p)
        union, ledger = {}, []
        for (mode, perm), generators in sorted(carriers.items()):
            def embed(alpha):
                out = [0]*3
                for i, a in enumerate(alpha):
                    out[perm[i]] = a
                return tuple(out)
            permuted = lru_cache(None)(lambda a: ev(embed(a)))
            part, row = one_ray_compile(n, permuted, mode, p, counts,
                                       actual_zero_products=False)
            for B, poly in part.items():
                basis = frozenset(embed(v) for v in B)
                polygon = [embed(v) for v in poly]
                if basis in union:
                    assert polygon_key(union[basis]) == polygon_key(polygon)
                    counts['common_factor_duplicate_whole_cone_checks'] += 1
                else:
                    union[basis] = polygon
            row.update({'mode': mode, 'permutation': perm,
                        'monomial_initial_carrier_generators': generators})
            ledger.append(row)
            counts['common_factor_initial_carrier_compilations'] += 1
        actual = fan(n, ev, p, counts)
        compare_fans(actual, union, counts)
        if any(shift):
            unshifted = fan(n, rotate_model(n, original, p), p, counts)
            compare_fans(actual, unshifted, counts)
            counts['entire_affine_translation_fan_equalities'] += 1
        assert sum(area2(poly) for poly in union.values()) == 1
        upper = 6*sum(1+4*j*(j-1).bit_length() for j in range(1, n+1))
        assert len(actual) <= upper
        counts['complete_common_factor_pencil_fans'] += 1
        rows.append({'fixture': name, 'length': n, 'characteristic': p,
                     'coordinate_change_new_in_old': ((0, -1, 1), (1, 1, -1), (-1, 0, 1)),
                     'coordinate_translation': shift,
                     'linear_common_factor': l0, 'linear_companions': (l1, l2),
                     'affine_constants': tuple(-sum(c*s for c, s in zip(l, shift))
                                               for l in (l0, l1, l2)),
                     'initial_carriers': ledger, 'full_cones': len(actual),
                     'uniform_six_carrier_bound': upper})

    f, g, h = (1, 1, 1), (1, 1, 0), (0, 1, 1)
    for p in (0, 2, 3, 101):
        for size, s, box in ((4, 3, False), ((2, 3), 3, True)):
            n, ev = axis_connected_model(size, s, p, box)
            check('rotated_transverse_axis_'+str(size), n, ev, f, g, h, p)
        n, ev = pascal_ray_model(3, p)
        check('rotated_embedded_square_zero_ray', n, ev, f, f, h, p)
        n, ev = axis_connected_model(4, 3, p)
        check('affine_translated_transverse_axis', n, ev, f, g, h, p,
              shift=(1, 2, 3))
    for p in (0, 101):
        points = ((1, 1, 0), (2, 1, 0), (1, 2, 0), (3, 2, 0),
                  (0, 0, 1), (0, 0, 2), (0, 0, 3))
        check('rotated_reduced_plane_and_axis', len(points),
              point_evaluate(points, p), f, g, h, p)
    return rows


def one_ray_controls(counts):
    rows, determinants = [], []
    def check(name, n, ev, mode, p):
        actual_evaluation_controls(n, ev, p, counts)
        compiled, row = one_ray_compile(n, ev, mode, p, counts)
        actual = fan(n, ev, p, counts)
        compare_fans(actual, compiled, counts)
        row.update({'fixture': name, 'length': n, 'characteristic': p,
                    'outside_ray': mode, 'full_cones': len(actual)})
        rows.append(row)
        counts['complete_one_ray_contraction_compilations'] += 1
        return actual, row

    for p in (0, 2, 3, 101):
        values = list(range(1, 7)) + ([8] if p == 2 else [9] if p == 3 else [])
        for m in values:
            n, ev = pascal_ray_model(m, p)
            actual, row = check('Pascal_three_rays_m='+str(m), n, ev, 'parallel', p)
            generic = not p or p > 3*m-1
            power, residue = bool(p), m
            if p:
                while residue > 1 and residue % p == 0:
                    residue //= p
                power = residue == 1
            if generic or power:
                formula = pascal_ray_predicted(m, collapse=power and not generic)
                compare_fans(actual, formula, counts)
                counts['entire_Pascal_ray_fans_vs_formula'] += 1
                assert len(actual) == (comb(m+2, 2) if generic else 3)
            row['generic_formula_applies'] = generic
            row['prime_power_collapse_applies'] = power
            _, walls, _, W = fan_budgets(actual)
            row['primitive_wall_cost'] = W
            if generic:
                assert len(walls) == 3*m*(m+1)//2
                expected = m*(m+1)*(m+2)+sum(abs(i-j) for i in range(m)
                                            for j in range(m-i))
                assert W == expected == actual_edge_length(actual)
                counts['exact_Pascal_triangular_graph_costs'] += 1
            if not p:
                for a in range(m+1):
                    for b in range(m-a+1):
                        c = m-a-b
                        exponents = tuple(range(a))+tuple(range(m, m+b))+tuple(range(2*m, 2*m+c))
                        numerator = prod(f-e for e, f in combinations(exponents, 2))
                        denominator = prod(factorial(j) for j in range(m))
                        assert numerator % denominator == 0
                        target = numerator//denominator
                        labels = (tuple((m+i, 0, 0) for i in range(a))+
                                  tuple((i, 1, 0) for i in range(b))+
                                  tuple((i, 0, 1) for i in range(c)))
                        matrix = [[ev(v).get(m+i, 0) for v in labels] for i in range(m)]
                        assert determinant(matrix) == target > 0
                        counts['independent_binomial_Vandermonde_determinants'] += 1
                        if m in (3, 6) and a == m//3 and b == m//3:
                            determinants.append({'m': m, 'composition': (a, b, c),
                                                 'exponents': exponents, 'matrix': matrix,
                                                 'determinant': target})
        for size, ray_length, box in ((4, 3, False), ((2, 3), 3, True)):
            n, ev = axis_connected_model(size, ray_length, p, box)
            _, row = check('socle_connected_axis_'+str(size), n, ev, 'axis', p)
            assert row['plane_ray_overlap'] == 1
        # A monomial plane plus square-zero parallel ray has no overlap.
        D = tuple([(i, j, 0) for j in range(3) for i in range(3)] +
                  [(i, 0, 1) for i in range(3)])
        index = {v: i for i, v in enumerate(D)}
        ev = lru_cache(None)(lambda a: {index[a]: 1} if a in index else {})
        _, row = check('monomial_plane_and_parallel_ray', len(D), ev, 'parallel', p)
        assert row['plane_ray_overlap'] == 0 and row['full_cones'] == 1
        print(json.dumps({'one_ray_characteristic': p, 'complete_fans': len(rows)}), flush=True)
    for p in (0, 101):
        points = ((1, 1, 0), (2, 1, 0), (1, 2, 0), (3, 2, 0),
                  (0, 0, 1), (0, 0, 2), (0, 0, 3))
        _, row = check('reduced_plane_and_transverse_axis', len(points),
                       point_evaluate(points, p), 'axis', p)
        assert row['plane_ray_overlap'] == 1
        check('planar_only_zero_ray', 4, point_evaluate(points[:4], p), 'parallel', p)
    return {'whole_fan_compilations': rows,
            'exact_integer_prefix_determinant_certificates': determinants,
            'uniform_statement': 'For z^2=yz=0, g is at most the sum over k from n-ell to min(m,n-1) of 1+4*(n-k)*ceil(log2(n-k)), hence O(n^2 log(n+1)). For xz=yz=0 the overlap is at most one, giving the sharp worst order Theta(n log(n+1)), with characteristic-zero lower examples.',
            'sharp_quadratic_fixture': 'K[t]/t^(2m), x=t, y=t^m(1+t)^m, z=t^m(1+t)^(2m). In characteristic zero or p>3m-1, the full fan has binomial(m+2,2) cones, a triangular lattice adjacency graph, and W=L=m(m+1)(m+2)+sum_(i+j<=m-1)|i-j|. At m=p^s in characteristic p it has exactly three cones.',
            'scope': 'Ordinary contraction and the saved planar matroid bound are established inputs. The contraction is not an algebra quotient. Exact priority is unassessed; the logarithm for the square-zero-ray class and the general two-planar-branch problem remain unresolved.'}


def source_hashes():
    paths = ['scripts/run_univariate_branch_gluing_checks.py',
             'docs/univariate_branch_gluing.md', 'docs/quadric_cross_fans.md',
             'docs/forced_socle_fan_steps.md', 'docs/planar_dyadic_bound.md',
             'docs/profile_section_budget.md', 'docs/prime_digit_lifts.md',
             'docs/univariate_cross_fans.md']
    pending = [p for p in paths if p.endswith('.py')]
    while pending:
        path = pending.pop()
        for node in ast.walk(ast.parse((ROOT/path).read_text())):
            names = ([node.module] if isinstance(node, ast.ImportFrom) and node.module
                     else [a.name for a in node.names] if isinstance(node, ast.Import) else [])
            for name in names:
                candidate = 'scripts/'+name.replace('.', '/')+'.py'
                if (ROOT/candidate).is_file() and candidate not in paths:
                    paths.append(candidate)
                    pending.append(candidate)
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def jsonable(value):
    if isinstance(value, Q):
        return value.numerator if value.denominator == 1 else [value.numerator, value.denominator]
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(v) for v in value]
    return value


def main():
    if not __debug__:
        raise RuntimeError('Assertions must be enabled.')
    before = source_hashes()
    counts, rows, fans = Counter(), [], []
    staircase_controls(counts)
    for p in (0, 2, 3, 101):
        for r in range(1, 6):
            n, M, unit, ev = sharp_model(r, p)
            row = recover('sharp_socle_r='+str(r), n, M, unit, ((0, 1),), p, counts,
                          gorenstein=r == 1)
            assert n == r*(r+1) and row['W_module_generators'] == r
            assert row['T_module_generators'] == r+1
            assert row['x_core_length'] == r+1
            rows.append(row)
            if r <= 4:
                actual = fan(n, ev, p, counts)
                expected = 3 if r == 1 else 4
                assert len(actual) == expected
                counts['complete_sharp_family_fans'] += 1
                fans.append({'r': r, 'characteristic': p, 'cones': len(actual),
                             'all_standard_sets': [sorted(B) for B in sorted(actual)]})
            if r <= 4:
                for t in (1, 2*r+1):
                    padded = axis_padding((n, M, unit), t, p)
                    row = recover('sharp_r='+str(r)+'_axis_padding='+str(t),
                                  *padded, ((0, 1), (-1, 1)), p, counts)
                    assert row['W_module_generators'] == r
                    assert (isqrt(4*padded[0]+1)-1)//2 == r
                    rows.append(row)
        for m in range(1, 7):
            row = recover('gorenstein_uniserial_m='+str(m), *uniserial(m, p),
                          ((0, 1),), p, counts, gorenstein=True)
            assert row['envelope_length'] == 2*row['length']
            assert row['T_module_generators'] == 2 and row['W_module_generators'] == 1
            rows.append(row)
        for N in (2, 3, 5):
            for kind in ('H', 'A', 'C'):
                reference, ev = forced_model(N, kind, p)
                M = tuple(tuple(ev(tuple(a+int(i == j) for j, a in enumerate(v)))
                                for v in reference) for i in range(3))
                rows.append(recover('forced_'+kind+'_N='+str(N), len(reference), M,
                                    ev((0, 0, 0)), ((0, 1),), p, counts,
                                    gorenstein=kind == 'A' and p != 2))
        for name, points in (
                ('single_axis_point', ((0, 0, 0),)),
                ('unit_y_branch', ((0, 1, 0),)),
                ('two_unit_branches', ((0, 1, 0), (1, 0, 1)))):
            M = tuple(tuple(clean({i: v[axis]}, p) for i, v in enumerate(points))
                      for axis in range(3))
            factors = tuple((-i, 1) for i in sorted({v[0] for v in points}))
            rows.append(recover(name, len(points), M,
                                {i: 1 for i in range(len(points))}, factors,
                                p, counts, gorenstein=True, reduced=True))
        print(json.dumps({'characteristic': p, 'completed_reconstructions': len(rows)}),
              flush=True)
    points = ((0, 0, 0), (1, 1, 0), (2, 2, 0), (1, 0, 1), (3, 0, 2), (2, 0, 0))
    for p in (0, 101):
        M = tuple(tuple(clean({i: v[axis]}, p) for i, v in enumerate(points))
                  for axis in range(3))
        rows.append(recover('six_reduced_points', len(points), M,
                            {i: 1 for i in range(len(points))},
                            tuple((-i, 1) for i in range(4)), p, counts,
                            gorenstein=True, reduced=True))
    X = ({1: 1}, {2: 1}, {3: 1}, {0: 1, 2: 1})
    Y = tuple(polynomial(X, {i: 1}, (1, 1, 1), 2) for i in range(4))
    row = recover('irreducible_quadratic_square_over_F2', 4, (X, Y, Y), {0: 1},
                  ((1, 1, 1),), 2, counts, gorenstein=True)
    assert row['W_module_generators'] == 1 and row['T_module_generators'] == 2
    assert row['intersection_length'] == 2 and row['envelope_length'] == 8
    rows.append(row)
    one_ray = one_ray_controls(counts)
    common_factors = common_factor_controls(counts)
    assert before == source_hashes(), 'A source file changed during verification.'
    report = {'command': 'python3 -B scripts/run_univariate_branch_gluing_checks.py',
              'status': 'passed', 'source_sha256': before, 'counts': dict(counts),
              'actual_quotient_reconstructions': rows, 'sharp_family_complete_fans': fans,
              'one_ray_full_fan_theorems': one_ray,
              'common_linear_factor_pencil_fans': common_factors,
              'uniform_common_factor_pencil_theorem': 'For every field and every length-n I containing l0*l1,l0*l2 for nonconstant affine linear l0 and affine linear l1,l2 with independent linear parts, g(I)<=6*sum_(j=1)^n (1+4*j*ceil(log2(j)))=O(n^2 log(n+1)). This uses actual original-coordinate matroids on at most six initial carriers; no fan invariance under a linear change is assumed.',
              'exact_fixtures': {
                  'sharp_family': 'H_r=(x^(r+1),yz,y(x,y)^r,z(x,z)^r), quotient by x^i(y^(r-i)-z^(r-i)) for 0<=i<r and x^r-y^r.',
                  'axis_padding': 'Direct product with K[x]/((x-1)^t), y=z=0; t=1 or 2r+1.',
                  'gorenstein_uniserial': 'K[t]/(t^(2m)), x=t, y=z=t^m.',
                  'reduced_points': points,
                  'irreducible_F2': {'polynomial': 'f=x^2+x+1', 'algebra': 'F2[x]/(f^2), y=z=f(x)',
                                     'x_matrix_columns': X, 'y_and_z_matrix_columns': Y}},
              'primary_sources': [
                  {'url': 'https://pi.math.cornell.edu/~frankmoore/Papers/connected.pdf',
                   'scope': 'Sections 1--2, equations (1.0.2), (1.1.3), (2.0.1), (2.1.3) and Theorem 2.8: classical fiber products, diagonal module quotients and additional Gorenstein hypotheses. No new gluing operation or exact-priority conclusion.'},
                  {'url': 'https://arxiv.org/html/1904.06297v3',
                   'scope': 'Introduction Theorem 1 and Theorem 4.6: graded connected sums require a common contraction of dual generators and an equality of annihilator ideals; no arbitrary nonhomogeneous fan-counting assertion follows.'},
                  {'url': 'https://arxiv.org/pdf/1406.7600',
                   'scope': 'Introduction and Theorem 4.3 concern connected sums over the residue field, rather than arbitrary shared-K[x] couplings. No exact quantitative-priority conclusion.'}],
              'uniform_proof': 'docs/univariate_branch_gluing.md',
              'limits': [
                  'Finite checks support implementations; all-field and all-length results rely on the self-contained proof.',
                  'The invariant-factor theorem, planar staircase estimate, fiber products, socle criterion and module quotients are established. Exact priority of the combined sharp resource inequalities is unassessed.',
                  'Module generator count is not vector-space dimension. The polynomial x-action and original monomial labels must be retained.',
                  'The sharp large-generator family has only three or four cones; it is not a large-fan witness.',
                  'The one-ray counts require z^2=yz=0 or xz=yz=0. They do not bound every yz=0 coupling. Their classical matroid ingredients and exact priority are distinguished in the proof.',
                  'The Pascal-ray example has a quadratic fan in characteristic zero and a constant fan at characteristic-p prime-power lengths. The square-zero-ray upper bound still has an unresolved logarithmic factor.',
                  'The common-factor pencil theorem covers arbitrary linear forms, but not a general pair of quadrics without a common linear factor. The two-quadratic-monomial corollary requires actual relations, not merely heads in one initial ideal.',
                  'No full yz=0 fan bound, unrestricted G3 improvement, supergeneric lower separation or completed goal follows.']}
    (ROOT/'results/research/univariate_branch_gluing_checks.json').write_text(
        json.dumps(jsonable(report), indent=2)+'\n')


if __name__ == '__main__':
    main()
