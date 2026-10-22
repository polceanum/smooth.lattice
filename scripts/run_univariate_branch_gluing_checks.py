#!/usr/bin/env python3
"""Exact reconstruction of actual yz=0 branch envelopes and module kernels."""
import ast
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import isqrt
from pathlib import Path
import hashlib
import json

from run_forced_socle_fan_checks import model as forced_model
from run_projection_tag_checks import fan
from run_signed_face_checks import field_basis, field_reduce

ROOT = Path(__file__).resolve().parents[1]


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


def source_hashes():
    paths = ['scripts/run_univariate_branch_gluing_checks.py',
             'docs/univariate_branch_gluing.md', 'docs/quadric_cross_fans.md',
             'docs/forced_socle_fan_steps.md']
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
    assert before == source_hashes(), 'A source file changed during verification.'
    report = {'command': 'python3 -B scripts/run_univariate_branch_gluing_checks.py',
              'status': 'passed', 'source_sha256': before, 'counts': dict(counts),
              'actual_quotient_reconstructions': rows, 'sharp_family_complete_fans': fans,
              'exact_fixtures': {
                  'sharp_family': 'H_r=(x^(r+1),yz,y(x,y)^r,z(x,z)^r), quotient by x^i(y^(r-i)-z^(r-i)) for 0<=i<r and x^r-y^r.',
                  'axis_padding': 'Direct product with K[x]/((x-1)^t), y=z=0; t=1 or 2r+1.',
                  'gorenstein_uniserial': 'K[t]/(t^(2m)), x=t, y=z=t^m.',
                  'reduced_points': points,
                  'irreducible_F2': {'polynomial': 'f=x^2+x+1', 'algebra': 'F2[x]/(f^2), y=z=f(x)',
                                     'x_matrix_columns': X, 'y_and_z_matrix_columns': Y}},
              'primary_sources': [
                  {'url': 'https://pi.math.cornell.edu/~frankmoore/Papers/connected.pdf',
                   'scope': 'Sections 1--2, equations (1.0.2), (1.1.3), (2.0.1), (2.1.3) and Theorem 2.8: classical fiber products, diagonal module quotients and additional Gorenstein hypotheses. No new gluing operation or exact-priority conclusion.'}],
              'uniform_proof': 'docs/univariate_branch_gluing.md',
              'limits': [
                  'Finite checks support implementations; all-field and all-length results rely on the self-contained proof.',
                  'The invariant-factor theorem, planar staircase estimate, fiber products, socle criterion and module quotients are established. Exact priority of the combined sharp resource inequalities is unassessed.',
                  'Module generator count is not vector-space dimension. The polynomial x-action and original monomial labels must be retained.',
                  'The sharp large-generator family has only three or four cones; it is not a large-fan witness.',
                  'No full yz=0 fan bound, unrestricted G3 improvement, supergeneric lower separation or completed goal follows.']}
    (ROOT/'results/research/univariate_branch_gluing_checks.json').write_text(
        json.dumps(jsonable(report), indent=2)+'\n')


if __name__ == '__main__':
    main()
