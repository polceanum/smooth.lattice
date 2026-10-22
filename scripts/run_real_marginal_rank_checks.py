#!/usr/bin/env python3
"""Exact real-algebraic equal-marginal examples; no numerical root tests."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import gcd
from pathlib import Path
import hashlib
import json
import sys

from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_interleaving_research import minimal_excluded_border
from run_projection_tag_checks import compare_fans
from run_quadratic_coordinate_fan_checks import dense_rank

ROOT = Path(__file__).resolve().parents[1]
TRIANGLE = [tuple(Q(i == j) for i in range(3)) for j in range(3)]


def multiply(a, b):
    """Negacyclic multiplication modulo X^d+1, d a power of two."""
    d = len(a)
    if not any(a[1:]):
        return tuple(a[0]*v for v in b)
    if not any(b[1:]):
        return tuple(b[0]*v for v in a)
    out = [Q(0)]*d
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    k = i+j
                    out[k % d] += x*y*(1 if k < d else -1)
    return tuple(out)


@lru_cache(maxsize=8192)
def inverse(a):
    if len(a) == 1:
        return (1/a[0],)
    e, o = a[::2], a[1::2]
    ee, oo = multiply(e, e), multiply(o, o)
    shifted = (-oo[-1],)+oo[:-1]
    norm = tuple(x-y for x, y in zip(ee, shifted))
    inv = inverse(norm)
    even, odd = multiply(e, inv), multiply(o, inv)
    return tuple(v for x, y in zip(even, odd) for v in (x, -y))


class Cyclo:
    """Q[X]/(X^d+1), with explicit quadratic-tower inversion."""
    def __init__(self, d, coefficients=0):
        self.d = d
        if isinstance(coefficients, Cyclo):
            assert d == coefficients.d
            self.v = coefficients.v
        elif isinstance(coefficients, (int, Q)):
            self.v = (Q(coefficients),)+(Q(0),)*(d-1)
        else:
            assert len(coefficients) == d
            self.v = tuple(map(Q, coefficients))

    def coerce(self, other):
        return Cyclo(self.d, other)

    def __add__(self, other):
        other = self.coerce(other)
        return Cyclo(self.d, tuple(a+b for a, b in zip(self.v, other.v)))

    __radd__ = __add__

    def __neg__(self):
        return Cyclo(self.d, tuple(-v for v in self.v))

    def __sub__(self, other):
        return self+-self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other)+-self

    def __mul__(self, other):
        if isinstance(other, (int, Q)):
            return Cyclo(self.d, tuple(v*other for v in self.v))
        other = self.coerce(other)
        return Cyclo(self.d, multiply(self.v, other.v))

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, (int, Q)):
            return self*(1/Q(other))
        other = self.coerce(other)
        return Cyclo(self.d, multiply(self.v, inverse(other.v)))

    def __rtruediv__(self, other):
        return self.coerce(other)/self

    def __pow__(self, k):
        if k < 0:
            return (1/self)**(-k)
        ans, v = Cyclo(self.d, 1), self
        while k:
            if k & 1:
                ans = ans*v
            v = v*v
            k //= 2
        return ans

    def __bool__(self):
        return any(self.v)

    def __eq__(self, other):
        return self.v == self.coerce(other).v


def reduce_column(v, pivots):
    v, coefficients = dict(v), defaultdict(int)
    for i, (row, rep) in sorted(pivots.items()):
        c = v.get(i, 0)
        if c:
            for j, a in row.items():
                v[j] = v.get(j, 0)-c*a
                if not v[j]:
                    del v[j]
            for j, a in rep.items():
                coefficients[j] += c*a
    return v, {j: c for j, c in coefficients.items() if c}


def basis(E, ev):
    chosen, pivots = [], {}
    for a in E:
        v, co = reduce_column(ev(a), pivots)
        if v:
            i = min(v)
            inv = 1/v[i]
            rep = {b: -c*inv for b, c in co.items()}
            rep[a] = inv
            pivots[i] = ({j: c*inv for j, c in v.items()}, rep)
            chosen.append(a)
    return tuple(chosen), pivots


def polynomial_add(P, R, shift, scale):
    P = dict(P)
    for a, c in R.items():
        b = tuple(x+y for x, y in zip(a, shift))
        P[b] = P.get(b, 0)+scale*c
        if not P[b]:
            del P[b]
    return P


def s_pairs(G, heads, w, counts):
    for i, j in combinations(range(len(G)), 2):
        lcm = tuple(max(a, b) for a, b in zip(heads[i], heads[j]))
        P = polynomial_add({}, G[i], tuple(a-b for a, b in zip(lcm, heads[i])), 1)
        P = polynomial_add(P, G[j], tuple(a-b for a, b in zip(lcm, heads[j])), -1)
        while P:
            a = max(P, key=lambda a: (dot(w, a), a))
            k = next((k for k, h in enumerate(heads) if all(b <= c for b, c in zip(h, a))), None)
            assert k is not None
            P = polynomial_add(P, G[k], tuple(b-c for b, c in zip(a, heads[k])), -P[a])
        counts['independent_original_number_field_S_pairs'] += 1


def actual_fan(n, ev, counts):
    E = tuple(hyperbolic_cross(n, 3))
    queue, cells, queried_edges = [((1, n+1, (n+1)**2), None)], {}, set()
    while queue:
        w, old = queue.pop()
        B = frozenset(basis(sorted(E, key=lambda a: (dot(w, a), a)), ev)[0])
        assert len(B) == n
        if old is not None:
            assert old != B
            queried_edges.add(frozenset((old, B)))
        if B in cells:
            continue
        _, pivots = basis(tuple(sorted(B)), ev)
        normals, G, heads = set(), [], tuple(minimal_excluded_border(tuple(B)))
        for a in heads:
            rem, co = reduce_column(ev(a), pivots)
            assert not rem
            G.append({a: 1, **{b: -c for b, c in co.items()}})
            for b in co:
                v = tuple(x-y for x, y in zip(a, b))
                d = gcd(*v)
                normals.add(tuple(x//d for x in v))
        assert all(dot(w, q) > 0 for q in normals)
        poly = list(TRIANGLE)
        for q in normals:
            poly = clip(poly, q, 1)
        assert poly and area2(poly) > 0
        cells[B] = poly
        s_pairs(G, heads, w, counts)
        counts['independent_original_number_field_border_relations'] += len(G)
        for a, b in zip(poly, poly[1:]+poly[:1]):
            if any(a[i] == b[i] == 0 for i in range(3)):
                continue
            q = next(q for q in normals if dot(q, a) == dot(q, b) == 0)
            nw, _ = exact_adjacent_weight(n, a, b, q)
            queue.append((nw, B))
    assert sum(map(area2, cells.values())) == 1
    counts['complete_original_number_field_fans'] += 1
    return cells


def chebyshev(d):
    A, B = [1], [0, 1]
    if not d:
        return A
    for j in range(1, d):
        C = [0]*(len(B)+1)
        for i, c in enumerate(B):
            C[i+1] += 2*c
        for i, c in enumerate(A):
            C[i] -= c
        A, B = B, C
    return B


def evaluate_poly(f, a):
    v = 0
    for c in reversed(f):
        v = v*a+c
    return v


def model(q, sign_tag, counts):
    assert q >= 2 and q & (q-1) == 0
    d = 4*q
    root = Cyclo(d, (0, 1)+(0,)*(d-2))
    assert root**d == -1 and root**(2*d) == 1
    nodes = tuple((root**(2*j+1)+root**(-(2*j+1)))/2 for j in range(2*q))
    c = (root**q+root**(-q))/2
    assert c*c == Q(1, 2)
    T = chebyshev(q)
    assert all(evaluate_poly(chebyshev(2*q), t) == 0 for t in nodes)
    composed = tuple(int(evaluate_poly(T, t) == c) for t in nodes)
    assert all(evaluate_poly(T, t) in (c, -c) for t in nodes)
    tags = tuple(int(j >= q) for j in range(2*q)) if sign_tag else composed
    assert sum(tags) == q and sum(composed) == q
    assert all(nodes[2*q-1-j] == -nodes[j] for j in range(2*q))
    positive = tuple(2+t for t in nodes)
    powers = [{0: Cyclo(d, 1)} for t in nodes]
    def power(j, k):
        P = powers[j]
        while len(P) <= k:
            h = len(P)
            P[h] = P[h-1]*positive[j]
        return P[k]
    @lru_cache(None)
    def ev(a):
        x, y, z = a
        return {j: power(j, y+z) for j in range(2*q) if not x or tags[j]}
    assert len(basis(tuple((0, j, 0) for j in range(2*q)), ev)[0]) == 2*q
    for j in range(2*q):
        assert ev((2, j, 0)) == ev((1, j, 0))
        assert ev((0, j, 1)) == ev((0, j+1, 0))
        counts['independent_original_coordinate_relation_columns'] += 2
    return tags, positive, ev


def predicted(q, sign_tag, counts):
    allowed = tuple(range(q+1)) if sign_tag else (0, q)
    moments = {b: (b, (2*q-b)*(2*q-b-1)//2+b*(b-1)//2) for b in allowed}
    cells = {}
    for axis in (1, 2):
        for b in allowed:
            B = frozenset(tuple(int(i == axis)*j for i in range(3)) for j in range(2*q-b))
            B |= frozenset(tuple(int(i == 0)+int(i == axis)*j for i in range(3)) for j in range(b))
            poly = clip(list(TRIANGLE), tuple(int(i == 3-axis)-int(i == axis) for i in range(3)), 1)
            for other in allowed:
                diff = tuple((moments[other][0]-moments[b][0]) if i == 0 else
                             (moments[other][1]-moments[b][1]) if i == axis else 0 for i in range(3))
                poly = clip(poly, diff, 1)
            assert poly and area2(poly) > 0
            cells[B] = poly
    assert len(cells) == (2*q+2 if sign_tag else 4)
    return cells


def one(q, sign_tag, counts):
    tags, nodes, ev = model(q, sign_tag, counts)
    n = 2*q
    ranks = []
    for b in range(q+1):
        B = tuple((0, j, 0) for j in range(n-b))+tuple((1, j, 0) for j in range(b))
        r = len(basis(B, ev)[0])
        assert r == (n if sign_tag else n-min(b, q-b))
        assert (r == n) == (sign_tag or b in (0, q))
        ranks.append(r)
        counts['independent_complete_two_row_basis_ranks'] += 1
    for exponents in ((0, 2), (0, 1, 3), tuple(range(q)), tuple(range(q+1))):
        for mask, limit in ((tuple(range(n)), n), (tuple(j for j in range(n) if tags[j]), q),
                            (tuple(j for j in range(n) if not tags[j]), q)):
            def restricted(a):
                return {j: v for j, v in ev((0, a, 0)).items() if j in mask}
            assert len(basis(exponents, restricted)[0]) == min(limit, len(set(exponents)))
            counts['independent_full_spark_planar_rank_oracle_controls'] += 1
    actual = actual_fan(n, ev, counts)
    expected = predicted(q, sign_tag, counts)
    compare_fans(actual, expected, counts)
    counts['independent_entire_algebraic_cone_equalities'] += len(actual)
    return {'q': q, 'length': n, 'tag': 'sign' if sign_tag else 'composed',
            'cyclotomic_field_degree': 4*q, 'two_row_ranks': ranks,
            'initial_ideals': len(actual)}


def independent_flag_controls(counts):
    """Separate rational matrices check the two uniform prefix-rank laws."""
    for q in range(2, 9):
        n = 2*q
        nodes = tuple(range(-q, 0))+tuple(range(1, q+1))
        for i, j in product(range(n+1), range(q+1)):
            columns = [tuple(Q(t**k) for t in nodes) for k in range(i)]
            columns += [tuple(Q(int(t < 0)*t**k) for t in nodes) for k in range(j)]
            rank = dense_rank(list(zip(*columns)), 0) if columns else 0
            assert rank == min(n, i+j)
            # x=(y^q+1)/2 in Q[y]/(y^(2q)-1); no reduction is needed here.
            columns = [tuple(Q(int(k == h)) for k in range(n)) for h in range(i)]
            columns += [tuple(Q(int(k == h)+int(k == q+h), 2) for k in range(n)) for h in range(j)]
            rank = dense_rank(list(zip(*columns)), 0) if columns else 0
            assert rank == i+j-min(j, max(0, i-q))
            counts['independent_rational_joint_prefix_rank_laws'] += 2
    for q in range(2, 65):
        n, b = 2*q, q//2
        w = (4*q-4*b+1, 2, 8*q*q)
        ordered = sorted(hyperbolic_cross(n+1, 3), key=lambda a: (dot(w, a), a))
        B = {(0, j, 0) for j in range(n-b)} | {(1, j, 0) for j in range(b)}
        assert set(ordered[:n]) == B
        assert dot(w, ordered[n-1]) == 4*q-2*b-1
        assert dot(w, ordered[n]) == 4*q-2*b
        counts['independent_strict_ordinary_prefix_rank_gap_witnesses'] += 1


def main():
    if not __debug__:
        raise RuntimeError('Exact assertions must remain enabled.')
    proofs = ['docs/real_marginal_rank_fans.md']
    before = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT/'scripts').glob('*.py')}
    before.update({f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in proofs})
    counts, rows = Counter(), []
    independent_flag_controls(counts)
    for d in (2, 4, 8, 16, 32):
        for trial in range(1, 5):
            a = tuple(Q(((i+1)*trial+3) % 7-3) for i in range(d))
            if any(a):
                assert multiply(a, inverse(a)) == (1,)+(0,)*(d-1)
                counts['independent_cyclotomic_tower_inverse_identities'] += 1
    for q in (2, 4, 8):
        for sign_tag in (False, True):
            row = one(q, sign_tag, counts)
            rows.append(row)
            print(row, flush=True)
    for q in range(2, 25):
        T, U = chebyshev(q), chebyshev(2*q)
        square = [0]*(2*len(T)-1)
        for i, a in enumerate(T):
            for j, b in enumerate(T):
                square[i+j] += 2*a*b
        square[0] -= 1
        assert square == U
        for b in range(q+1):
            a = 2*q-b
            blocks = (tuple(range(0, a, 2))+tuple(range(1, b, 2)),
                      tuple(range(1, a, 2))+tuple(range(0, b, 2)))
            assert all(len(S) == len(set(S)) == q for S in blocks)
            for j in range(q+1):
                assert (j*(2*q-2*b)+j*j-2*q*j)-(b*(2*q-2*b)+b*b-2*q*b) == (j-b)**2
            counts['independent_parity_block_and_exposure_identities'] += 1
        counts['independent_Chebyshev_polynomial_identities'] += 1
    loaded = {Path(m.__file__).resolve() for m in tuple(sys.modules.values())
              if getattr(m, '__file__', None) and Path(m.__file__).resolve().parent == ROOT/'scripts'}
    loaded.add(Path(__file__).resolve())
    sources = proofs+sorted(str(p.relative_to(ROOT)) for p in loaded)
    assert all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest() == before[f] for f in sources)
    report = {'all_checks_passed': True,
              'scope': 'Exact common real-algebraic point configurations, original whole border fans and S-pairs over cyclotomic number fields; independently clipped entire predicted cones and rank oracles. Uniform all-q and all-set assertions follow from the proof, not finite checks. No unrestricted G3 improvement, polynomial supergeneric family or exact novelty claim.',
              'models': rows, 'counts': dict(counts), 'loaded_source_stability_checked': True,
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    target = ROOT/'results/research/real_marginal_rank_checks.json'
    target.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'all_checks_passed': True, 'models': rows, 'counts': dict(counts)}), flush=True)


if __name__ == '__main__':
    main()
