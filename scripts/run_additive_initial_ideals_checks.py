#!/usr/bin/env python3
"""Exact controls for the classical additive-polynomial obstruction.

Frobenius reduction is checked against ordinary polynomial reduction and
independent quotient columns. No novelty or unrestricted growth is claimed.
"""
from itertools import combinations, product
from fractions import Fraction
from math import comb, prod
from pathlib import Path
from random import Random
import hashlib
import json

from run_structural_complexity_checks import hyperbolic_cross, positive_chambers, primitive

ROOT = Path(__file__).resolve().parents[1]


def require(condition, label):
    if not condition:
        raise AssertionError(label)


class Field:
    def __init__(self, p, extension=False):
        self.p, self.extension = p, extension
        self.size = 256 if extension else p
        require(not extension or p == 2, "the bounded extension implementation is binary")
        self.table = [[self.raw_mul(a, b) for b in range(self.size)] for a in range(self.size)]
        self.inverses = [None]+[next((b for b in range(1, self.size) if self.table[a][b] == 1), None)
                                for a in range(1, self.size)]
        require(all(x is not None for x in self.inverses[1:]), "the quotient arithmetic really is a field")

    def add(self, a, b):
        return a ^ b if self.extension else (a+b) % self.p

    def neg(self, a):
        return a if self.extension else -a % self.p

    def raw_mul(self, a, b):
        if not self.extension:
            return a*b % self.p
        out = 0
        while b:
            if b & 1:
                out ^= a
            a <<= 1
            if a & 256:
                a ^= 0x11b
            b >>= 1
        return out

    def mul(self, a, b):
        return self.table[a][b]

    def power(self, a, k):
        out = 1
        while k:
            if k & 1:
                out = self.mul(out, a)
            a = self.mul(a, a)
            k >>= 1
        return out

    def total(self, values):
        out = 0
        for a in values:
            out = self.add(out, a)
        return out


def inverse(matrix, field):
    n = len(matrix)
    rows = [list(row)+[int(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i][j]), None)
        if pivot is None:
            return None
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scale = field.inverses[rows[j][j]]
        rows[j] = [field.mul(x, scale) for x in rows[j]]
        for i in range(n):
            if i != j:
                scale = field.neg(rows[i][j])
                rows[i] = [field.add(x, field.mul(scale, y)) for x, y in zip(rows[i], rows[j])]
    return [row[n:] for row in rows]


def add_scaled(poly, other, scale, field, shift=None, frobenius=1):
    out = dict(poly)
    for exponent, value in other.items():
        target = tuple(a*frobenius+(shift[i] if shift else 0) for i, a in enumerate(exponent))
        out[target] = field.add(out.get(target, 0), field.mul(scale, field.power(value, frobenius)))
        if not out[target]:
            del out[target]
    return out


def head(poly, weight):
    return max(poly, key=lambda a: (sum(x*y for x, y in zip(a, weight)), a))


def power_degree(a, p):
    r = 0
    while a > 1 and a % p == 0:
        a //= p
        r += 1
    require(a == 1, "a pure exponent is a power of the characteristic")
    return r


def pure_data(a, p):
    nonzero = [i for i, x in enumerate(a) if x]
    require(len(nonzero) == 1, "additive nonconstant terms are pure powers")
    i = nonzero[0]
    return i, power_degree(a[i], p)


def frobenius_basis(generators, weight, field):
    basis = [dict(g) for g in generators]
    reductions = 0
    while True:
        found = None
        for i, j in combinations(range(len(basis)), 2):
            hi, hj = head(basis[i], weight), head(basis[j], weight)
            vi, ai = pure_data(hi, field.p)
            vj, aj = pure_data(hj, field.p)
            if vi == vj:
                found = (i, j) if ai >= aj else (j, i)
                break
        if found is None:
            break
        i, j = found
        hi, hj = head(basis[i], weight), head(basis[j], weight)
        power = hi[pure_data(hi, field.p)[0]]//hj[pure_data(hj, field.p)[0]]
        denominator = field.power(basis[j][hj], power)
        scale = field.neg(field.mul(basis[i][hi], field.inverses[denominator]))
        changed = add_scaled(basis[i], basis[j], scale, field, frobenius=power)
        require(not changed or (sum(x*y for x, y in zip(head(changed, weight), weight)), head(changed, weight)) <
                (sum(x*y for x, y in zip(hi, weight)), hi), "Frobenius elimination decreases the ordinary head")
        if changed:
            basis[i] = changed
        else:
            basis.pop(i)
        reductions += 1
    basis = [{a: field.mul(c, field.inverses[g[head(g, weight)]]) for a, c in g.items()} for g in basis]
    # Reduce tails by Frobenius substitutions. This keeps ordinary leading heads.
    for i in range(len(basis)):
        while True:
            hi = head(basis[i], weight)
            found = None
            for a in basis[i]:
                if a == hi or not any(a):
                    continue
                v, r = pure_data(a, field.p)
                for j, g in enumerate(basis):
                    hj = head(g, weight)
                    u, s = pure_data(hj, field.p)
                    if i != j and v == u and r >= s:
                        candidate = (a, j, field.p**(r-s))
                        if found is None or (sum(x*y for x, y in zip(a, weight)), a) > (
                                sum(x*y for x, y in zip(found[0], weight)), found[0]):
                            found = candidate
                        break
            if not found:
                break
            a, j, power = found
            basis[i] = add_scaled(basis[i], basis[j], field.neg(basis[i][a]), field, frobenius=power)
            require(head(basis[i], weight) == hi, "tail reduction preserves the head")
            reductions += 1
    return basis, reductions


def ordinary_remainder(poly, basis, weight, field):
    todo, remainder = dict(poly), {}
    while todo:
        a = head(todo, weight)
        for g in basis:
            b = head(g, weight)
            if all(x >= y for x, y in zip(a, b)):
                scale = field.neg(field.mul(todo[a], field.inverses[g[b]]))
                todo = add_scaled(todo, g, scale, field, shift=tuple(x-y for x, y in zip(a, b)))
                break
        else:
            remainder[a] = todo.pop(a)
    return remainder


def gaussian_basis(ground, columns, weight, field, n):
    pivots, selected = {}, []
    for a in sorted(ground, key=lambda a: (sum(x*y for x, y in zip(a, weight)), a)):
        v = list(columns[a])
        for j, row in sorted(pivots.items()):
            scale = field.neg(v[j])
            if scale:
                v = [field.add(x, field.mul(scale, y)) for x, y in zip(v, row)]
        pivot = next((i for i, x in enumerate(v) if x), None)
        if pivot is not None:
            pivots[pivot] = [field.mul(x, field.inverses[v[pivot]]) for x in v]
            selected.append(a)
            if len(selected) == n:
                return tuple(sorted(selected))
    raise AssertionError("independent quotient-column rank")


def compositions(k, d):
    if d == 1:
        yield (k,)
    else:
        for a in range(k+1):
            yield from ((a,)+b for b in compositions(k-a, d-1))


def audit(generators, field, k, d, residue, model, extra=None, full_weights=True, frobenius_step=1):
    n = field.p**k
    require(k % frobenius_step == 0, "length is a power of the declared Frobenius base")
    base, base_k = field.p**frobenius_step, k//frobenius_step
    ground = tuple(hyperbolic_cross(n, d))
    columns = {a: residue(a) for a in ground}
    canonical_weights = {tuple(base**(base_k-a) for a in allocation) for allocation in compositions(base_k, d)}
    weights = set(canonical_weights)
    input_degree = max(pure_data(a, field.p)[1] for g in generators for a in g if any(a))
    require(all(not any(a) or pure_data(a, field.p)[1] % frobenius_step == 0
                for g in generators for a in g), "input equations use the declared Frobenius base")
    step_bound = (len(generators)+d)*d*(input_degree//frobenius_step+base_k+1)
    if full_weights and d == 2:
        cuts = sorted({Fraction(base**max(r, 0), base**max(-r, 0)) for r in range(-base_k, base_k+1)})
        samples = [(a+b)/2 for a, b in zip([Fraction(0)]+cuts, cuts)]+[cuts[-1]+1]
        weights.update((t.numerator, t.denominator) for t in samples)
    if full_weights and d == 3:
        pure = [tuple((i == j)*base**r for i in range(d)) for j in range(d) for r in range(base_k+1)]
        normals = {primitive(tuple(x-y for x, y in zip(a, b))) for a, b in combinations(pure, 2)
                   if any(x > y for x, y in zip(a, b)) and any(x < y for x, y in zip(a, b))}
        weights.update(positive_chambers(sorted(normals)))
    states, reductions, s_pairs = {}, 0, 0
    for weight in sorted(weights):
        gb, steps = frobenius_basis(generators, weight, field)
        if weight in canonical_weights:
            require(steps <= step_bound,
                    "the explicit pure-term reduction-step budget at canonical weights")
        reductions += steps
        require(len(gb) == d, "finite additive quotient has one head for every variable")
        allocation = [None]*d
        for g in gb:
            a = head(g, weight)
            j, r = pure_data(a, field.p)
            require(r % frobenius_step == 0, "all heads retain the declared Frobenius spacing")
            require(allocation[j] is None, "leading variables are distinct")
            allocation[j] = r
            require(len(g) <= base_k+2, "reduced additive basis has logarithmically many terms")
            require(all(not any(a) or pure_data(a, field.p)[1] % frobenius_step == 0 for a in g),
                    "all reduced terms retain the declared Frobenius base")
            require(all(field.total(field.mul(c, residue(b)[i]) for b, c in g.items()) == 0 for i in range(n)),
                    "Frobenius basis independently vanishes in the target quotient")
        require(sum(allocation) == k, "ordinary quotient length is p^k")
        require(all(not ordinary_remainder(g, gb, weight, field) for g in generators),
                "original generators reduce to zero by ordinary polynomial division")
        for g, h in combinations(gb, 2):
            a, b = head(g, weight), head(h, weight)
            common = tuple(max(x, y) for x, y in zip(a, b))
            s = add_scaled({}, g, 1, field, shift=tuple(x-y for x, y in zip(common, a)))
            s = add_scaled(s, h, field.neg(1), field, shift=tuple(x-y for x, y in zip(common, b)))
            require(not ordinary_remainder(s, gb, weight, field), "ordinary Buchberger S-pair")
            s_pairs += 1
        expected = tuple(product(*(range(field.p**a) for a in allocation)))
        actual = gaussian_basis(ground, columns, weight, field, n)
        require(tuple(sorted(expected)) == actual, "ordinary monomial Gaussian selection agrees with the additive box")
        allocation = tuple(allocation)
        encoded = [[[list(a), c] for a, c in sorted(g.items())] for g in sorted(gb, key=lambda g: head(g, weight))]
        if allocation in states:
            require(states[allocation]['reduced_basis'] == encoded, "same initial ideal gives exactly the same reduced basis")
        else:
            states[allocation] = {'weight': weight, 'reduced_basis': encoded}
    require(len(states) <= comb(base_k+d-1, d-1), "the sharp composition count bound")
    # Every realized state must occur at its coefficient-independent canonical weight.
    for allocation in states:
        weight = tuple(field.p**(k-a) for a in allocation)
        gb, _ = frobenius_basis(generators, weight, field)
        heads = {pure_data(head(g, weight), field.p) for g in gb}
        require(heads == set(enumerate(allocation)), "canonical weights cover every realized initial ideal")
        require(all(sum(x*y for x, y in zip(head(g, weight), weight)) == n and
                    all(a == head(g, weight) or sum(x*y for x, y in zip(a, weight)) <= n//base for a in g)
                    for g in gb), "the uniform factor-q head/tail gap")
    unique_polynomials = {tuple((tuple(a), c) for a, c in g) for state in states.values() for g in state['reduced_basis']}
    total_terms = sum(len(g) for g in unique_polynomials)
    require(len(unique_polynomials) <= d*comb(base_k+d-1, d-1) and
            total_terms <= d*(base_k+2)*comb(base_k+d-1, d-1), "the sparse universal output bounds")
    return {'model': model, 'characteristic': field.p, 'field_size': field.size,
            'quotient_length': n, 'logarithmic_length': base_k, 'p_logarithmic_length': k,
            'frobenius_base': base, 'variables': d,
            'generators': [[[list(a), c] for a, c in sorted(g.items())] for g in generators],
            'initial_ideal_count': len(states), 'sharp_class_bound': comb(base_k+d-1, d-1),
            'weights_checked': len(weights), 'frobenius_reductions': reductions, 'ordinary_s_pairs': s_pairs,
            'canonical_reduction_step_budget_checked': True,
            'input_q_degree': input_degree//frobenius_step, 'canonical_reduction_step_bound': step_bound,
            'universal_basis_polynomials': len(unique_polynomials), 'universal_basis_terms': total_terms,
            'full_additive_comparison_arrangement_checked': full_weights and d in (2, 3),
            'states': [{'box_exponents': [r//frobenius_step for r in a],
                        'side_lengths': [field.p**r for r in a], **states[a]}
                       for a in sorted(states)], **(extra or {})}


def local_graph(p, k, support, frobenius_step=1):
    field, n = Field(p), p**k
    powers = [[1]+[0]*(n-1)]
    for j in range(n-1):
        powers.append([sum(powers[-1][r-s] for s in support if r >= s) % p for r in range(n)])
    residue = lambda a: tuple([0]*min(a[0], n)+powers[a[1]][:max(n-a[0], 0)]) if a[1] < n else (0,)*n
    return audit([{(n, 0): 1}, {(0, 1): 1, **{(j, 0): -1 % p for j in support}}],
                 field, k, 2, residue, 'nilpotent_additive_graph', {'support': support}, frobenius_step=frobenius_step)


def generic_point_control(field, k, d, translate=None, base=None):
    base = base or field.p
    step = power_degree(base, field.p)
    scalars = [a for a in range(field.size) if field.power(a, base) == a]
    require(len(scalars) == base, "the exact coefficient field contains the required scalar subfield")
    allocations = list(compositions(k, d))
    random = Random(101+k*17+d)
    for attempt in range(1, 101):
        vectors = [[random.randrange(1, field.size) for i in range(d)] for j in range(k)]
        matrices = []
        for allocation in allocations:
            matrix = [[field.power(vector[i], base**r) for i in range(d) for r in range(allocation[i])]
                      for vector in vectors]
            matrices.append(matrix)
        if all(inverse(matrix, field) is not None for matrix in matrices):
            break
    else:
        raise AssertionError('bounded deterministic specialization did not meet every Moore minor')
    points = [tuple(field.total(field.mul(c, v[i]) for c, v in zip(coefficients, vectors)) for i in range(d))
              for coefficients in product(scalars, repeat=k)]
    require(len(set(points)) == base**k and len({x[0] for x in points}) == len(points),
            "the rational F_q-linear point space has dimension k and separating first coordinate")
    polynomial = [1]
    for point in points:
        nxt = [0]*(len(polynomial)+1)
        for j, c in enumerate(polynomial):
            nxt[j] = field.add(nxt[j], field.mul(field.neg(point[0]), c))
            nxt[j+1] = field.add(nxt[j+1], c)
        polynomial = nxt
    generators = [{tuple([j]+[0]*(d-1)): c for j, c in enumerate(polynomial) if c}]
    matrix = [[field.power(vector[0], base**r) for r in range(k)] for vector in vectors]
    inv = inverse(matrix, field)
    for i in range(1, d):
        coefficients = [field.total(field.mul(x, v[i]) for x, v in zip(row, vectors)) for row in inv]
        g = {tuple(int(j == i) for j in range(d)): 1}
        for r, c in enumerate(coefficients):
            if c:
                g[tuple([base**r]+[0]*(d-1))] = field.neg(c)
        generators.append(g)
    if translate:
        zero = (0,)*d
        for g in generators:
            constant = field.neg(field.total(field.mul(c, prod_value(a, translate, field)) for a, c in g.items()))
            if constant:
                g[zero] = constant
        points = [tuple(field.add(a, b) for a, b in zip(point, translate)) for point in points]
    residue = lambda a: tuple(prod_value(a, x, field) for x in points)
    fixture = audit(generators, field, k*step, d, residue, 'affine_additive_point_space' if translate else 'finite_additive_point_space',
                    {'basis_vectors': vectors, 'points': points, 'specialization_attempts': attempt,
                     'translation': translate, 'all_composition_moore_minors_nonzero': True},
                    full_weights=d <= 3, frobenius_step=step)
    require(fixture['initial_ideal_count'] == len(allocations), "the general upper bound is attained")
    return fixture


def prod_value(a, point, field):
    out = 1
    for r, x in zip(a, point):
        out = field.mul(out, field.power(x, r))
    return out


def tensor_local_control():
    field = Field(2)
    cells = tuple(product(range(4), range(2)))
    index = {a: i for i, a in enumerate(cells)}
    z = {(1, 0): 1, (0, 1): 1, (2, 0): 1}
    powers = [[int(a == (0, 0)) for a in cells]]
    for r in range(8):
        nxt = [0]*8
        for a, c in zip(cells, powers[-1]):
            for b, t in z.items():
                target = tuple(x+y for x, y in zip(a, b))
                if target in index:
                    nxt[index[target]] ^= c & t
        powers.append(nxt)
    residue = lambda a: tuple(powers[a[2]][index[(i-a[0], j-a[1])]]
                              if (i-a[0], j-a[1]) in index else 0 for i, j in cells) if a[2] < 9 else (0,)*8
    generators = [{(4, 0, 0): 1}, {(0, 2, 0): 1},
                  {(0, 0, 1): 1, (1, 0, 0): 1, (0, 1, 0): 1, (2, 0, 0): 1}]
    return audit(generators, field, 3, 3, residue, 'noncyclic_nilpotent_tensor_graph')


def multiplicative_contrast(p, k):
    """Independent complete fan for (x^(p^k)-x-1, xy-1)."""
    field, n = Field(p), p**k
    def multiply(a, b):
        out = [0]*(2*n-1)
        for i, c in enumerate(a):
            for j, t in enumerate(b):
                out[i+j] = (out[i+j]+c*t) % p
        for i in range(len(out)-1, n-1, -1):
            out[i-n] = (out[i-n]+out[i]) % p
            out[i-n+1] = (out[i-n+1]+out[i]) % p
        return out[:n]
    x, y, one = [0, 1]+[0]*(n-2), [-1 % p]+[0]*(n-2)+[1], [1]+[0]*(n-1)
    require(multiply(x, y) == one, "the multiplicative contrast has exact Laurent inverse")
    powers = []
    for variable in (x, y):
        seq = [one]
        for i in range(n):
            seq.append(multiply(seq[-1], variable))
        powers.append(seq)
    residue = lambda a: tuple(multiply(powers[0][a[0]], powers[1][a[1]]))
    ground = tuple(hyperbolic_cross(n, 2))
    comparisons = tuple(hyperbolic_cross(n+1, 2))
    cuts = sorted({Fraction(-(a[1]-b[1]), a[0]-b[0]) for a, b in combinations(comparisons, 2)
                   if (a[0]-b[0])*(a[1]-b[1]) < 0})
    samples = [(a+b)/2 for a, b in zip([Fraction(0)]+cuts, cuts)]+[cuts[-1]+1]
    columns = {a: residue(a) for a in ground}
    states = {}
    for t in samples:
        weight = (t.numerator, t.denominator)
        basis = gaussian_basis(ground, columns, weight, field, n)
        require(all(not (a and b) for a, b in basis), "the mixed relation forces a hook basis")
        states.setdefault(basis, weight)
    require(len(states) == n, "one multiplicative relation escapes the logarithmic additive bound")
    return {'model': 'separable_affine_additive_polynomial_and_one_product', 'characteristic': p,
            'logarithmic_length': k, 'quotient_length': n, 'initial_ideal_count': len(states),
            'additive_class_bound_at_same_length': k+1, 'complete_comparison_intervals': len(samples),
            'generators': [[[n, 0, 1], [1, 0, -1 % p], [0, 0, -1 % p]], [[1, 1, 1], [0, 0, -1 % p]]],
            'reduced_over_algebraic_closure': True,
            'witnesses': [{'row_lengths': [sum(a[1] == j for a in basis)
                                          for j in range(max(a[1] for a in basis)+1)], 'weight': states[basis]}
                          for basis in sorted(states)]}


def main():
    fixtures = [local_graph(p, k, support) for p, k, support in
                ((2, 3, (1, 4)), (2, 4, (1, 2, 4, 8)), (2, 5, (1, 4, 8, 16)),
                 (3, 2, (1, 3)), (3, 3, (1, 3, 9)))]
    extension = Field(2, extension=True)
    fixtures.extend(generic_point_control(extension, k, d) for k, d in ((3, 2), (4, 2), (3, 3), (4, 3), (2, 4)))
    fixtures.append(generic_point_control(extension, 3, 3, translate=(17, 33, 65)))
    fixtures.append(generic_point_control(extension, 2, 3, base=4))
    fixtures.append(local_graph(2, 4, (1, 4), frobenius_step=2))
    fixtures.append(tensor_local_control())
    report = {'schema_version': 1, 'status': 'passed',
              'command': 'python3 -B scripts/run_additive_initial_ideals_checks.py',
              'extension_field': {'characteristic': 2, 'degree': 8,
                                  'modulus_low_to_high': [1, 1, 0, 1, 1, 0, 0, 0, 1],
                                  'encoding': 'integer binary digits are polynomial coefficients',
                                  'every_nonzero_element_inverse_checked': True},
              'fixtures': fixtures, 'multiplicative_contrasts': [multiplicative_contrast(2, 3), multiplicative_contrast(3, 2)],
              'new_unrestricted_asymptotic': False,
              'novelty_established': False, 'goal_complete': False,
            'scope': 'A classical Frobenius reduction gives power-of-q box initial ideals for finite ideals generated by affine q-linearized polynomials. Generic additive point spaces attain the composition bound. These are restricted structural controls, not an unrestricted extremal lower family.',
              'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                for name in ('scripts/run_additive_initial_ideals_checks.py',
                                             'scripts/run_structural_complexity_checks.py')}}
    path = ROOT/'results/research/additive_initial_ideals_checks.json'
    path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'fixtures': len(fixtures),
                      'states': sum(f['initial_ideal_count'] for f in fixtures),
                      'ordinary_s_pairs': sum(f['ordinary_s_pairs'] for f in fixtures),
                      'counts': [[f['characteristic'], f['variables'], f['quotient_length'], f['initial_ideal_count']]
                                 for f in fixtures], 'goal_complete': False, 'output': str(path.relative_to(ROOT))}))


if __name__ == '__main__':
    main()
