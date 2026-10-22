#!/usr/bin/env python3
"""Exact vertical-sweep charges, sharp products and an induced generic wall."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
from random import Random
import hashlib
import json
import sys

from run_structural_complexity_checks import area2, clip, dot
from run_planar_dyadic_bound_checks import convex_hull
from run_coordinate_slice_budget_checks import primitive, profile
from run_quadric_cross_checks import sharp_model, full_fan
from run_generic_hook_checks import (model, matching_basis, exchange_normals,
                                    minor_labels, bareiss, evaluate_minor, check_family)
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_square_zero_extension_checks import extension, matrix_controls
from run_signed_face_checks import normals as border_normals, field_basis

ROOT = Path(__file__).resolve().parents[1]
PROOF = 'docs/cross_layer_direction_budget.md'
OLD_CERT = 'results/research/generic_hook_checks.json'
PROOF_INPUTS = ('docs/layer_shadow_budget.md','docs/rational_shadow_bound.md',
                'docs/quadric_cross_fans.md')
TRIANGLE = [tuple(Q(i==j) for i in range(3)) for j in range(3)]


def cells_from_profiles(points):
    points = sorted(set(points))
    cells = {}
    for v in points:
        poly = list(TRIANGLE)
        for u in points:
            if u!=v:
                poly = clip(poly,tuple(b-a for a,b in zip(v,u)),1)
        if area2(poly)>0:
            cells[v] = poly
    assert sum(map(area2,cells.values()))==1
    return cells


def sweep_certificate(points,counts):
    points = sorted(set(points))
    cells = cells_from_profiles(points)
    shadow = convex_hull([(v[1],v[2]) for v in points])
    owners = defaultdict(list)
    for v,poly in cells.items():
        for a,b in zip(poly,poly[1:]+poly[:1]):
            if not any(a[i]==b[i]==0 for i in range(3)):
                owners[tuple(sorted((a,b)))].append(v)
    lower = defaultdict(list)
    cross = []
    for segment,pair in owners.items():
        assert len(pair)==2
        v,u = pair
        q = tuple(b-a for a,b in zip(v,u))
        assert all(dot(q,w)==0 for w in segment)
        if q[0]:
            low = min(pair,key=lambda a:a[0])
            lower[low].append(segment)
            cross.append((segment,pair))
        counts['independent_positive_polygon_adjacencies'] += 1
    charged,projected = set(),set()
    uncharged = []
    for v,poly in cells.items():
        if lower[v]:
            segment = sorted(lower[v])[0]
            assert segment not in charged
            charged.add(segment)
            counts['injective_lower_cross_edge_charges'] += 1
        else:
            yz = v[1],v[2]
            assert yz in shadow and yz not in projected
            assert v[0]==min(u[0] for u in points if u[1:]==v[1:])
            projected.add(yz)
            assert any(a[0]==b[0]==0 and a!=b
                       for a,b in zip(poly,poly[1:]+poly[:1]))
            uncharged.append(v)
            counts['injective_shadow_vertex_charges'] += 1
    assert len(cells)==len(charged)+len(uncharged)
    assert len(cells)<=len(cross)+len(shadow)
    counts['whole_polytope_sweep_inequalities'] += 1
    return {'input_profiles':len(points),'positive_vertices':len(cells),
            'cross_edge_pieces':len(cross),'projected_hull_vertices':len(shadow),
            'lower_edge_charges':len(charged),'shadow_charges':len(uncharged)}


def polytope_controls(counts):
    examples = [
        [(0,0,0)],
        [(0,0,1),(0,1,0),(0,2,2)],
        [(0,1,1),(1,1,1),(2,1,1)],
        [(0,2,0),(0,0,2),(1,0,0)],
        [(0,0,4),(0,1,1),(0,4,0),(2,0,0)],
        [(0,0,3),(1,0,2),(2,0,1),(3,0,0)],
    ]
    rng = Random(20261022)
    for _ in range(24):
        labels = [tuple(rng.randrange(0,8) for _ in range(3)) for j in range(8)]
        r = rng.randrange(1,5)
        # Full uniform-matroid basis profiles, followed by arbitrary subpolytopes.
        points = [tuple(sum(labels[j][i] for j in B) for i in range(3))
                  for B in combinations(range(len(labels)),r)]
        examples.append(points)
        examples.append(rng.sample(points,min(12,len(points))))
    return [sweep_certificate(v,counts) for v in examples]


def sharp_controls(counts):
    rows = []
    for L in (1,2,3,4):
        data = full_fan(sharp_model(L),counts)
        bases = [frozenset(map(tuple,c['standard_set'])) for c in data['cells']]
        points = [profile(B) for B in bases]
        sweep = sweep_certificate(points,counts)
        cells = {B:[tuple(Q(x) for x in v) for v in c['positive_cone']]
                 for B,c in zip(bases,data['cells'])}
        directions = Counter()
        owners = defaultdict(list)
        for B,poly in cells.items():
            for a,b in zip(poly,poly[1:]+poly[:1]):
                if not any(a[i]==b[i]==0 for i in range(3)):
                    owners[tuple(sorted((a,b)))].append(B)
        for pair in owners.values():
            assert len(pair)==2
            q = primitive(tuple(a-b for a,b in zip(profile(pair[0]),profile(pair[1]))))
            assert q[0]==1
            directions[q] += 1
        expected = {(1,-(2*L-2*b-1),0) for b in range(L)}
        expected |= {(1,0,-(2*L-2*b-1)) for b in range(L)}
        assert set(directions)==expected and set(directions.values())=={L+1}
        assert data['initial_ideals']==(L+1)**2
        rows.append({'L':L,'length':data['length'],'initial_ideals':data['initial_ideals'],
                     'distinct_cross_directions':len(directions),
                     'parallel_pieces_per_direction':L+1,'sweep':sweep})
        counts['actual_sharp_product_direction_inventories'] += 1
    return rows


def support_cone(B,columns,rank):
    N = exchange_normals(B,columns,rank)
    poly = list(TRIANGLE)
    for q in N:
        poly = clip(poly,q,1)
    assert area2(poly)>0
    return N,poly


def induced_wall_control(counts):
    original = json.loads((ROOT/OLD_CERT).read_text())
    family = next(f for f in original['families'] if f['thickness']==2)
    old_check = check_family(family,counts)
    assert old_check['generic_initial_ideals_in_every_characteristic']==89
    D,index,G,allowed,P,columns,fixed = model(2)
    target = (-1,5,-1)
    direct = {(-1,G[g][0]-i,G[g][1]-j) for g,i,j in P}
    assert target not in direct and tuple(-v for v in target) not in direct
    found = None
    for raw in family['support_fan_bases']:
        B = frozenset(map(tuple,raw))
        N,poly = support_cone(B,columns,len(D))
        for a,b in zip(poly,poly[1:]+poly[:1]):
            if a==b or any(a[i]==b[i]==0 for i in range(3)):
                continue
            q = next(q for q in N if dot(q,a)==dot(q,b)==0)
            if primitive(q)==primitive(target):
                found = B,N,poly,a,b,q
                break
        if found:
            break
    assert found
    B,N,poly,a,b,q = found
    w,_ = exact_adjacent_weight(32,a,b,q)
    C,_ = matching_basis(sorted(columns,key=lambda a:(dot(w,a),a)),columns,len(D))
    assert C!=B
    M,other = support_cone(C,columns,len(D))
    midpoint = tuple((u+v)/2 for u,v in zip(a,b))
    assert all(x>0 for x in midpoint)
    for normals in (N,M):
        assert all(dot(v,midpoint)>=0 for v in normals)
        assert all(primitive(v)==primitive(target) for v in normals if dot(v,midpoint)==0)
    Bs = [frozenset(map(tuple,x)) for x in family['support_fan_bases']]
    assert B in Bs and C in Bs
    Ls = [minor_labels(V,index,columns) for V in (B,C)]
    chosen = None
    for evaluation in family['integer_evaluations']:
        values = evaluation['coefficients']
        determinants = [bareiss(evaluate_minor(L,values)) for L in Ls]
        if all(determinants):
            chosen = values,determinants
            break
    assert chosen
    values,determinants = chosen
    tails = [{v:values[P.index((g,*v))] for v in A} for g,A in enumerate(allowed)]
    n,ev,ground = extension(D,D,G,tails,0)
    matrix_controls(n,ev,tuple((x,*v) for x in range(2) for v in D),ground,0,counts)
    for V,cone in ((B,poly),(C,other)):
        full = tuple(sorted(fixed|V))
        assert len(field_basis(full,ev,0)[0])==32
        normals,relations = border_normals(full,ev,0)
        assert all(dot(v,midpoint)>=0 for v in normals)
        assert any(dot(v,midpoint)==0 for v in normals)
        assert all(primitive(v)==primitive(target) for v in normals if dot(v,midpoint)==0)
        assert all(dot(v,u)>=0 for v in normals for u in cone)
        counts['independent_original_quotient_wall_sides'] += 1
    counts['full_parameter_induced_wall_counterexamples'] += 1
    return {'thickness':2,'length':32,'primitive_direction':target,
            'original_parameter_direction_count':len(direct),
            'absent_from_original_parameters':True,
            'residual_endpoint_bases':[sorted(B),sorted(C)],
            'wall_endpoints':[[str(v) for v in p] for p in (a,b)],
            'strict_wall_interior_weight':[str(v) for v in midpoint],
            'common_integer_coefficients':values,
            'nonzero_endpoint_determinants':determinants,
            'all_characteristic_scope':'Both endpoints are among the89 previously certified primitive coherent minors; the complete thickness-two certificate is independently replayed here. Generic simultaneous survival holds over every characteristic, not necessarily at the common rational specialization.'}


def main():
    if not __debug__:
        raise RuntimeError('Exact assertions must remain enabled.')
    before = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT/'scripts').glob('*.py')}
    for f in (PROOF,OLD_CERT,*PROOF_INPUTS):
        before[f] = hashlib.sha256((ROOT/f).read_bytes()).hexdigest()
    counts = Counter()
    polytopes = polytope_controls(counts)
    sharp = sharp_controls(counts)
    wall = induced_wall_control(counts)
    loaded = {Path(m.__file__).resolve() for m in tuple(sys.modules.values())
              if getattr(m,'__file__',None) and Path(m.__file__).resolve().parent==ROOT/'scripts'}
    loaded.add(Path(__file__).resolve())
    sources = [PROOF,OLD_CERT,*PROOF_INPUTS]+sorted(str(p.relative_to(ROOT)) for p in loaded)
    assert all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==before[f] for f in sources)
    report = {'all_checks_passed':True,'polytope_sweep_controls':polytopes,
              'actual_sharp_products':sharp,'induced_generic_wall':wall,
              'counts':dict(counts),'source_sha256':{f:before[f] for f in sources},
              'scope':'Uniform direction-inventory reduction from an elementary sweep and the existing layered shadow theorem; sharp up to a logarithm on an existing quadratic family. Original coefficient arrow weights are insufficient. No unconditional exponent improvement, generic growing-thickness survival theorem, exact priority, polynomial lower separation or goal completion.'}
    (ROOT/'results/research/cross_layer_direction_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_checks_passed':True,'counts':dict(counts),
                      'sharp_products':sharp,'induced_wall':wall['primitive_direction']}),flush=True)


if __name__=='__main__':
    main()
