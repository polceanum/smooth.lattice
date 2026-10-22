#!/usr/bin/env python3
"""Exact two-colour braids, independently checked by all-basis feasibility."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path
from random import Random
import hashlib
import json
import sys

from run_quadratic_coordinate_fan_checks import dense_rank, quadratic_tensor
from run_projection_tag_checks import fan, compare_fans, point_model
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import area2, clip, dot, hyperbolic_cross
from run_canonical_fan_deflation_checks import marked_basis
from run_algebraic_incidence_controls import pair_checks

ROOT = Path(__file__).resolve().parents[1]
TRIANGLE = [tuple(Q(i == j) for i in range(3)) for j in range(3)]


def sample(lo, hi):
    if lo is None:
        return hi-1 if hi is not None else Q(0)
    return lo+1 if hi is None else (lo+hi)/2


def matrix_ranker(columns, char):
    memo = {frozenset(): 0}
    def rank(S):
        S = frozenset(S)
        if S not in memo:
            memo[S] = dense_rank([[columns[j][i] for j in sorted(S)]
                                  for i in range(len(columns[0]))], char)
        return memo[S]
    return rank


def greedy(order, rank):
    S = []
    for e in order:
        if rank(S+[e]) > len(S):
            S.append(e)
    return tuple(S)


def joint_pairs(A, B, rank, counts):
    R = {(i, j): rank(A[:i]+B[:j])
         for i in range(len(A)+1) for j in range(len(B)+1)}
    pairs = {(A[i-1], B[j-1]) for i in range(1, len(A)+1)
             for j in range(1, len(B)+1)
             if R[i, j-1]+R[i-1, j]-R[i-1, j-1]-R[i, j] == 1}
    assert len({a for a, b in pairs}) == len(pairs)
    assert len({b for a, b in pairs}) == len(pairs)
    for i, j in R:
        assert R[i, j] == i+j-sum(a in A[:i] and b in B[:j] for a, b in pairs)
        counts['independent_joint_prefix_rank_identities'] += 1
    return pairs


def strict_witness(profile, profiles):
    """Eliminate t from all strict inequalities; then solve one interval."""
    inequalities, lower, upper = [], [], []
    for other in profiles:
        if other == profile:
            continue
        dt, du, dc = (b-a for a, b in zip(profile, other))
        if dt:
            (lower if dt > 0 else upper).append((Q(-du, dt), Q(-dc, dt)))
        else:
            inequalities.append((du, dc))
    inequalities += [(a-c, b-d) for a, b in upper for c, d in lower]
    lo = hi = None
    for a, b in inequalities:
        if not a:
            if b <= 0:
                return None
        elif a > 0:
            lo = max(lo, Q(-b, a)) if lo is not None else Q(-b, a)
        else:
            hi = min(hi, Q(-b, a)) if hi is not None else Q(-b, a)
    if lo is not None and hi is not None and lo >= hi:
        return None
    u = sample(lo, hi)
    t = sample(max((a*u+b for a, b in lower), default=None),
               min((a*u+b for a, b in upper), default=None))
    assert all(sum(w*(b-a) for w, a, b in zip((t, u, 1), profile, v)) > 0
               for v in profiles if v != profile)
    return t, u


def eligible_interval(E, F, pair, ab, u, rank, counts):
    e, f = pair
    c = ab[e][0]*u+ab[e][1]
    assert c == ab[f][0]*u+ab[f][1]
    C = [v for v in E if ab[v][0]*u+ab[v][1] < c]
    ordered = sorted(F, key=lambda v: (ab[v][0]*u+ab[v][1], v))
    yes = []
    for j in range(len(ordered)+1):
        S = C+ordered[:j]
        r = rank(S)
        if rank(S+[e]) == rank(S+[f]) == rank(S+[e, f]) == r+1:
            yes.append(j)
        counts['independent_contracted_pair_eligibility_states'] += 1
    assert not yes or yes == list(range(yes[0], yes[-1]+1))
    counts['independent_single_interval_exchange_events'] += 1
    endpoints = []
    if yes:
        if yes[0]:
            endpoints.append(c-(ab[ordered[yes[0]-1]][0]*u+ab[ordered[yes[0]-1]][1]))
        if yes[-1] < len(ordered):
            endpoints.append(c-(ab[ordered[yes[-1]]][0]*u+ab[ordered[yes[-1]]][1]))
    return set(endpoints)


def compile_braid(E0, E1, ab, rank, counts, generic=True):
    E0, E1 = tuple(E0), tuple(E1)
    E = E0+E1
    h = rank(E0)+rank(E1)-rank(E)
    P = sum(len(Ei)*(len(Ei)-1)//2 for Ei in (E0, E1))
    event = defaultdict(list)
    for colour, Ei in enumerate((E0, E1)):
        for a, b in combinations(Ei, 2):
            if ab[a][0] != ab[b][0]:
                event[Q(ab[b][1]-ab[a][1], ab[a][0]-ab[b][0])].append((colour, a, b))
    cuts = sorted(event)
    old, bends, crossings, found = None, 0, 0, set()
    for lo, hi in zip([None]+cuts, cuts+[None]):
        u = sample(lo, hi)
        A = greedy(sorted(E0, key=lambda i: (ab[i][0]*u+ab[i][1], i)), rank)
        B = greedy(sorted(E1, key=lambda i: (ab[i][0]*u+ab[i][1], i)), rank)
        pairs = joint_pairs(A, B, rank, counts)
        assert len(pairs) == h
        slopes = [(ab[a][0]-ab[b][0], ab[a][1]-ab[b][1]) for a, b in pairs]
        if old is not None:
            oldmap, newmap = defaultdict(Counter), defaultdict(Counter)
            for a, b in old:
                oldmap[a*lo+b][a] += 1
            for a, b in slopes:
                newmap[a*lo+b][a] += 1
            assert {v: sum(c.values()) for v, c in oldmap.items()} == {
                v: sum(c.values()) for v, c in newmap.items()}
            changed = {v: sum((oldmap[v]-newmap[v]).values()) for v in oldmap}
            bends += sum(changed.values())
            if generic:
                assert len(event[lo]) == 1
                colour, e, f = event[lo][0]
                allowed = eligible_interval((E0, E1)[colour], (E1, E0)[colour],
                                            (e, f), ab, lo, rank, counts)
                if colour:
                    allowed = {-v for v in allowed}
                assert {v for v, c in changed.items() if c} <= allowed
                assert sum(changed.values()) <= 2
                counts['independent_bends_at_exchange_interval_endpoints'] += 1
            counts['independent_threshold_multiset_continuities'] += 1
        old = slopes
        roots = []
        for (a, b), (c, d) in combinations(slopes, 2):
            if a != c:
                v = Q(d-b, a-c)
                if (lo is None or lo < v) and (hi is None or v < hi):
                    roots.append(v)
        crossings += len(roots)
        inner = sorted(set(roots))
        for l, r in zip([lo]+inner, inner+[hi]):
            u = sample(l, r)
            levels = sorted(set(a*u+b for a, b in slopes))
            for z, w in zip([None]+levels, levels+[None]):
                t = sample(z, w)
                direct = frozenset(greedy(sorted(E, key=lambda i: (
                    int(i in E1)*t+ab[i][0]*u+ab[i][1], i)), rank))
                matched = {i for pair in pairs for i in pair}
                predicted = set(A+B)-matched
                predicted.update(a if ab[a][0]*u+ab[a][1] < t+ab[b][0]*u+ab[b][1]
                                 else b for a, b in pairs)
                assert direct == predicted
                found.add(direct)
                counts['independent_braid_regions_vs_greedy'] += 1
    if generic:
        assert bends <= 2*P
        assert crossings <= max(0, h-1)*bends+h*(h-1)//2
    bound = 1+h+h*(h-1)//2+(3*h+1)*P if h else 1+P
    assert len(found) <= bound
    counts['entire_two_colour_braids'] += 1
    return found, {'elements': len(E), 'rank': rank(E), 'overlap': h,
                   'same_colour_pair_budget': P, 'bends': bends,
                   'proper_crossing_pairs': crossings, 'optimal_regions': len(found),
                   'proved_region_bound': bound}


def all_basis_check(E0, E1, ab, rank, found, counts):
    E = tuple(E0)+tuple(E1)
    bases = [frozenset(S) for S in combinations(E, rank(E)) if rank(S) == len(S)]
    profiles = defaultdict(list)
    for S in bases:
        v = (len(S & set(E1)), sum(ab[i][0] for i in S), sum(ab[i][1] for i in S))
        profiles[v].append(S)
    exposed = set()
    for v, Bs in profiles.items():
        w = strict_witness(v, profiles)
        if w is not None:
            assert len(Bs) == 1
            t, u = w
            assert frozenset(greedy(sorted(E, key=lambda i: (
                int(i in E1)*t+ab[i][0]*u+ab[i][1], i)), rank)) == Bs[0]
            exposed.add(Bs[0])
        counts['independent_all_basis_strict_halfplane_feasibilities'] += 1
    assert exposed == found
    counts['independent_complete_all_basis_fan_equalities'] += 1
    return bases


def one_model(name, columns, char, seed, counts, rank=None, E0=None, E1=None):
    rank = matrix_ranker(columns, char) if rank is None else rank
    m = len(columns)
    E0 = tuple(range(m//2)) if E0 is None else tuple(E0)
    E1 = tuple(i for i in range(m) if i not in E0) if E1 is None else tuple(E1)
    rng = Random(seed)
    ab = [(rng.randrange(-10000, 10001), rng.randrange(-10000, 10001)) for i in range(m)]
    found, row = compile_braid(E0, E1, ab, rank, counts)
    all_basis_check(E0, E1, ab, rank, found, counts)
    row.update(model=name, characteristic=char)
    return row


def actual_controls(counts):
    rows = []
    for char, constant, linear in product((0, 2, 3), (0, 1), (0,)):
        n, ev = quadratic_tensor(char, 2, constant, linear)
        ground = tuple(hyperbolic_cross(n, 2))
        raw = [tuple((j, *a) for a in ground) for j in (0, 1)]
        cuts = [Q(0)]+sorted({Q(b, a) for a, b in product(range(1, n), repeat=2)})+[None]
        chosen = next((lo, hi) for lo, hi in zip(cuts, cuts[1:])
                      if lo < Q(7, 5) and (hi is None or Q(7, 5) < hi))
        lo, hi = chosen
        u = sample(lo, hi)
        layers = [field_basis(sorted(E, key=lambda a: (a[1]*u+a[2], a)), ev, char)[0] for E in raw]
        labels = layers[0]+layers[1]
        columns = [tuple(ev(a).get(i, 0) for i in range(n)) for a in labels]
        rank = matrix_ranker(columns, char)
        E0, E1 = tuple(range(len(layers[0]))), tuple(range(len(layers[0]), len(labels)))
        ab = [a[1:] for a in labels]
        found, row = compile_braid(E0, E1, ab, rank, counts, generic=False)
        bases = all_basis_check(E0, E1, ab, rank, found, counts)
        predicted = {}
        for S in bases:
            v = tuple(sum(labels[j][i] for j in S) for i in range(3))
            poly = clip(list(TRIANGLE), (0, lo.denominator, -lo.numerator), 1)
            if hi is not None:
                poly = clip(poly, (0, -hi.denominator, hi.numerator), 1)
            for T in bases:
                q = tuple(sum(labels[j][i] for j in T)-v[i] for i in range(3))
                poly = clip(poly, q, 1)
                if not poly or area2(poly) == 0:
                    break
            if poly and area2(poly) > 0:
                predicted[frozenset(labels[j] for j in S)] = poly
        original = fan(n, ev, char, counts)
        for S, poly in original.items():
            G, _ = marked_basis(S, ev, char)
            w = tuple(sum(v[i] for v in poly)/len(poly) for i in range(3))
            counts['independent_original_quotient_S_pairs'] += pair_checks(G, w, char)
        actual = {}
        for S, poly in original.items():
            poly = clip(poly, (0, lo.denominator, -lo.numerator), 1)
            if hi is not None:
                poly = clip(poly, (0, -hi.denominator, hi.numerator), 1)
            if poly and area2(poly) > 0:
                actual[S] = poly
        compare_fans(actual, predicted, counts)
        counts['independent_original_quotient_complete_strip_equalities'] += 1
        row.update(model='original_quadratic_tensor', characteristic=char,
                   length=n, constant=constant, strip=[str(lo), str(hi)],
                   original_full_fan_cones=len(original), original_strip_cones=len(actual))
        rows.append(row)
    return rows


def main():
    if not __debug__:
        raise RuntimeError('Exact assertions must remain enabled.')
    proofs = ['docs/two_colour_threshold_braids.md', 'docs/quadratic_intersection_fans.md']
    before = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT/'scripts').glob('*.py')}
    before.update({f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in proofs})
    counts, rows, rng = Counter(), [], Random(424242)
    for h, char, trial in product(range(1, 5), (0, 2, 3, 101), range(3)):
        while True:
            M = [[rng.randrange(3 if not char else char) for j in range(h)] for i in range(h)]
            if dense_rank(M, char) == h:
                break
        cols = [tuple(int(i == j) for i in range(h)) for j in range(h)]+list(zip(*M))
        rows.append(one_model('two_represented_bases', cols, char, rng.randrange(10**8), counts))
    # An abstract Vamos matroid; exchange is checked directly, not presumed.
    blocks = [frozenset((2*i, 2*i+1)) for i in range(4)]
    excluded = {blocks[a] | blocks[b] for a, b in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3))}
    bases = {frozenset(S) for S in combinations(range(8), 4)}-excluded
    for A, B in product(bases, repeat=2):
        for a in A-B:
            assert any((A-{a}) | {b} in bases for b in B-A)
            counts['independent_nonrepresentable_basis_exchange_controls'] += 1
    def vrank(S):
        S = frozenset(S)
        return max(len(S & B) for B in bases)
    A, B, C, D = blocks
    ingleton_left = sum(vrank(X | Y) for X, Y in ((A, B), (A, C), (A, D), (B, C), (B, D)))
    ingleton_right = vrank(A)+vrank(B)+vrank(C | D)+vrank(A | B | C)+vrank(A | B | D)
    assert (ingleton_left, ingleton_right) == (15, 16)
    counts['independent_Ingleton_nonrepresentability_certificates'] += 1
    rows.append(one_model('abstract_Vamos', [None]*8, None, 91243, counts,
                          rank=vrank, E0=(0, 1, 2, 4)))
    for columns in (
            [(1, 0), (1, 0), (0, 0), (0, 1), (0, 1), (1, 1)],
            [tuple(int(i == j) for i in range(6)) for j in range(6)]):
        rows.append(one_model('loops_parallel_or_zero_overlap', columns, 0, 67891, counts))
    rows.extend(actual_controls(counts))
    loaded = {Path(m.__file__).resolve() for m in tuple(sys.modules.values())
              if getattr(m, '__file__', None) and Path(m.__file__).resolve().parent == ROOT/'scripts'}
    loaded.add(Path(__file__).resolve())
    sources = proofs+sorted(str(p.relative_to(ROOT)) for p in loaded)
    assert all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest() == before[f] for f in sources)
    report = {'all_checks_passed': True, 'scope': 'Exact abstract two-colour braids and independent all-basis strict halfplane feasibility; original actual-quotient cone patches independently reconstructed. The uniform cubic proof is analytical. No sharp cubic lower family, unrestricted G3 improvement, polynomial separation or exact novelty is claimed.',
              'models': rows, 'counts': dict(counts), 'loaded_source_stability_checked': True,
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    target = ROOT/'results/research/two_colour_braid_checks.json'
    target.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'all_checks_passed': True, 'models': len(rows),
                      'counts': dict(counts), 'report': str(target.relative_to(ROOT))}), flush=True)


if __name__ == '__main__':
    main()
