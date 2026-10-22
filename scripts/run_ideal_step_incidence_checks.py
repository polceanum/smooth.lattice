#!/usr/bin/env python3
"""Exact adjacent point-ideal fans and arithmetic incidence controls."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import comb, gcd, isqrt, lcm
from pathlib import Path
import hashlib
import json
import sys

from run_projection_tag_checks import fan, point_model
from run_signed_face_checks import field_basis
from run_structural_complexity_checks import cross, dot, hyperbolic_cross
from run_canonical_fan_deflation_checks import marked_basis
from run_algebraic_incidence_controls import pair_checks

ROOT = Path(__file__).resolve().parents[1]


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


def primitive(v, unoriented=False):
    den = lcm(*(a.denominator if isinstance(a, Q) else 1 for a in v))
    v = tuple(int(a*den) for a in v)
    d = gcd(*v)
    assert d
    v = tuple(a//d for a in v)
    if unoriented and next(a for a in v if a) < 0:
        v = tuple(-a for a in v)
    return v


def ceil_sqrt(h):
    r = isqrt(h)
    return r + (r*r < h)


def fan_budgets(cells):
    rays, walls, wall_counts = Counter(), {}, Counter()
    for poly in cells.values():
        P = clean(poly)
        for a in P:
            if min(a) > 0:
                rays[primitive(a)] += 1
        for a, b in zip(P, P[1:]+P[:1]):
            if any(a[i] == b[i] == 0 for i in range(3)):
                continue
            edge = tuple(sorted((a, b)))
            q = primitive(cross(a, b), True)
            assert min(q) < 0 < max(q)
            if edge in walls:
                assert walls[edge] == q
            walls[edge] = q
            wall_counts[edge] += 1
    assert set(wall_counts.values()) <= {2}
    U = sum(ceil_sqrt(max(w))*d for w, d in rays.items())
    W = sum(sum(map(abs, q)) for q in walls.values())
    assert U <= 2*W
    return rays, walls, U, W


def triples(N):
    for A in range(1, N-1):
        for B in range(1, N-A):
            yield A, B, N-A-B


def triple_ray(A, B, C):
    return primitive((B*C, A*C, A*B))


@lru_cache(None)
def S(N):
    return sum(ceil_sqrt(max(triple_ray(*t))) for t in triples(N))


def phi(n):
    return sum(gcd(a, n) == 1 for a in range(1, n+1))


def divisor_increment(n):
    return 3*sum(d*phi(d) for d in range(2, n+1) if n % d == 0)


def direct_W(n):
    return 3*sum((a+b)//gcd(a, b) for a in range(1, n)
                 for b in range(1, n+1-a))


def local_histogram(p, K, e):
    """Exact distribution of the two smallest capped valuations' sum."""
    if not K:
        return {0: Q(1)}
    if not e:
        result = {0: 1-Q(3, p*p)}
        for j in range(1, K):
            result[j] = 3*(1-Q(1, p*p))/p**(2*j)
        result[K] = Q(3, p**(2*K))
        return {j: v for j, v in result.items() if v}
    result = {0: 1-Q(1, p*p)}
    result.update({j+2: v/p**2 for j, v in local_histogram(p, K-1, e-1).items()})
    return result


def valuation(a, p, K):
    if not a:
        return K
    e = 0
    while e < K and a % p == 0:
        a //= p
        e += 1
    return e


def exact_local_weight(hist, p):
    # a+b/sqrt(p), represented by two rational numbers.
    out = [Q(0), Q(0)]
    for j, probability in hist.items():
        out[j % 2] += probability/p**(j//2)
    return tuple(out)


def local_controls(counts):
    rows = []
    for p, K in ((2, 5), (3, 3), (5, 2), (7, 2)):
        modulus = p**K
        values = [valuation(a, p, K) for a in range(modulus)]
        for N in range(modulus):
            hist = Counter()
            for a, b in product(range(modulus), repeat=2):
                vals = sorted((values[a], values[b], values[(N-a-b) % modulus]))
                hist[vals[0]+vals[1]] += 1
            hist = {j: Q(v, modulus**2) for j, v in hist.items()}
            expected = local_histogram(p, K, values[N])
            assert hist == expected and sum(hist.values()) == 1
            assert exact_local_weight(hist, p) == exact_local_weight(expected, p)
            counts['independent_complete_local_residue_distributions'] += 1
            counts['independent_residue_pairs'] += modulus**2
        rows.append({'prime': p, 'valuation_cap': K, 'modulus': modulus,
                     'all_total_residues_checked': modulus,
                     'weight_arithmetic': 'exact rational coefficients in Q(sqrt(p))'})
    # Rational certificates in the infinite-progression proof.
    assert Q(1, 3) < Q(3, 5)**2  # 1/sqrt(3)<3/5.
    assert Q(8, 9)/Q(6, 7)+Q(1, 27) == Q(29, 27)
    assert Q(1)/(1-Q(9, 200)) == Q(200, 191) < Q(29, 27)
    counts['exact_progression_gap_certificates'] += 3
    return rows


def arithmetic_controls(counts):
    rows = []
    previous = 0
    for n in range(1, 151):
        W = direct_W(n)
        inc = divisor_increment(n)
        assert W-previous == inc <= 3*n*(n-1)
        assert W == 3*sum(d*phi(d)*(n//d) for d in range(2, n+1))
        if n >= 2 and phi(n) == n-1:
            assert inc == 3*n*(n-1)
            counts['independent_prime_increment_equalities'] += 1
        counts['independent_totient_wall_recurrences'] += 1
        previous = W
    for N in (15, 45, 105, 315, 525, 735, 945):
        lower, center, upper = S(N-2), S(N), S(N+2)
        positive, negative = 3*(center-lower), 3*(upper-center)
        if N >= 105:
            assert positive > 0 > negative
        rows.append({'total': N, 'S_minus_two': lower, 'S_center': center,
                     'S_plus_two': upper,
                     'U_increment_at_length_N_minus_one': positive,
                     'U_increment_at_length_N_plus_one': negative})
        counts['independent_exact_integer_incidence_sums'] += 3
        print(rows[-1], flush=True)
    assert rows[2]['U_increment_at_length_N_minus_one'] == 71364
    assert rows[2]['U_increment_at_length_N_plus_one'] == -28866
    for N in range(3, 41):
        for d in range(2, N+1):
            actual = sum(a % d == b % d == 0 for a, b, c in triples(N))
            m = (N-1)//d
            assert actual == m*(m-1)//2
            counts['independent_divisible_pair_tail_identities'] += 1
    return rows


def actual_controls(counts):
    rows = []
    for p, maximum in ((0, 7), (101, 11)):
        previous = None
        for n in range(2, maximum+1):
            def inverse(a):
                return pow(a, -1, p) if p else Q(1, a)
            points = tuple((j, inverse(j), inverse(j-1)) for j in range(2, n+2))
            _, ev = point_model(points, p)
            # Separator F_{n-1}(x) divided by its last-point value.
            values = []
            for point in points:
                v = 1
                for j in range(2, n+1):
                    v *= point[0]-j
                values.append(v % p if p else v)
            assert all(v == 0 for v in values[:-1]) and values[-1]
            normalizer = inverse(values[-1])
            v = tuple(x*normalizer % p if p else x*normalizer for x in values)
            assert v == (0,)*(n-1)+(1,)
            for axis in range(3):
                assert all((points[j][axis]-points[-1][axis])*v[j] == 0 for j in range(n))
                counts['independent_actual_separator_eigenline_actions'] += 1
            for point in points:
                x, y, z = point
                for value in (x*y-1, (x-1)*z-1, y*z-z+y):
                    assert (value % p if p else value) == 0
                    counts['independent_actual_curve_relations'] += 1
            cells = fan(n, ev, p, counts)
            rays, walls, U, W = fan_budgets(cells)
            predicted = {triple_ray(*t) for N in (n, n+1) for t in triples(N)}
            assert set(rays) == predicted and set(rays.values()) == {3}
            assert len(cells) == comb(n+1, 2) and len(rays) == (n-1)**2
            assert U == 3*(S(n)+S(n+1)) and W == direct_W(n)
            assert len(walls) == 3*comb(n, 2)
            counts['independent_entire_actual_ray_and_wall_inventories'] += 1
            for B, poly in cells.items():
                G, _ = marked_basis(B, ev, p)
                w = tuple(sum(a[i] for a in poly)/len(poly) for i in range(3))
                counts['independent_original_reduced_basis_S_pairs'] += pair_checks(G, w, p)
            if previous is not None:
                oldcells, oldrays, oldev = previous
                shared = {triple_ray(*t) for t in triples(n)}
                assert set(rays) & set(oldrays) == shared
                assert len(set(rays)-set(oldrays)) == comb(n, 2)
                assert len(set(oldrays)-set(rays)) == comb(n-2, 2)
                # Same global order on both actual point quotients.
                ground = tuple(hyperbolic_cross(n, 3))
                for C, poly in oldcells.items():
                    w = tuple(sum(a[i] for a in poly)/len(poly) for i in range(3))
                    B = frozenset(field_basis(sorted(ground, key=lambda a: (dot(w, a), a)), ev, p)[0])
                    assert C < B
                    counts['independent_actual_ideal_step_standard_set_inclusions'] += 1
                counts['independent_exact_adjacent_ray_intersections'] += 1
            previous = cells, rays, ev
            rows.append({'characteristic': p, 'length': n, 'initial_cones': len(cells),
                         'interior_rays': len(rays), 'undirected_wall_pieces': len(walls),
                         'U': U, 'W': W})
            print(rows[-1], flush=True)
    return rows


def main():
    if not __debug__:
        raise RuntimeError('Exact assertions must be enabled.')
    proof_sources = ['docs/ideal_step_incidence.md', 'docs/weighted_facet_incidence.md']
    before = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT/'scripts').glob('*.py')}
    before.update({f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in proof_sources})
    counts = Counter()
    local = local_controls(counts)
    actual = actual_controls(counts)
    arithmetic = arithmetic_controls(counts)
    loaded = {Path(m.__file__).resolve() for m in tuple(sys.modules.values())
              if getattr(m, '__file__', None) and
              Path(m.__file__).resolve().parent == ROOT/'scripts'}
    loaded.add(Path(__file__).resolve())
    sources = proof_sources + sorted(str(p.relative_to(ROOT)) for p in loaded)
    assert all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest() == before[f] for f in sources)
    report = {'all_checks_passed': True,
              'scope': 'Actual nested radical point ideals and complete exact border fans; finite residue distributions, totient recurrences and larger integer sums. These support, and do not computationally prove, the uniform weighted asymptotic or the infinite cubic-jump theorem. No unrestricted budget, extremal solution or priority certificate is claimed.',
              'actual_point_ideals': actual, 'local_residue_controls': local,
              'integer_incidence_sums': arithmetic, 'counts': dict(counts),
              'loaded_source_stability_checked': True,
              'source_sha256': {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sources}}
    target = ROOT/'results/research/ideal_step_incidence_checks.json'
    target.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'all_checks_passed': True, 'counts': dict(counts),
                      'report': str(target.relative_to(ROOT))}), flush=True)


if __name__ == '__main__':
    main()
