"""Exact controls for the dyadic planar matroid profile bound.

The proof is in docs/planar_dyadic_bound.md. This checks its exchange path,
integer-direction counting and every dyadic budget by independent greedy replay.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations,product
from pathlib import Path
from random import Random
from math import gcd
import hashlib
import json,sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_structural_complexity_checks import hyperbolic_cross


def sweep(ground,rank,independent,H,label,signs=(1,1)):
    n=len(ground)
    assert H>=1 and all(x>=0 and y>=0 and (x+1)*(y+1)<=H for x,y in ground)
    sigma,tau=signs
    def oracle_greedy(weight):
        basis=set()
        for j in sorted(range(n),key=lambda j:(ground[j][0]*weight[0]+ground[j][1]*weight[1],ground[j],j)):
            if independent(frozenset(basis|{j})):basis.add(j)
        assert len(basis)==rank
        return basis
    events={}
    for i,j in combinations(range(n),2):
        dx,dy=(a-b for a,b in zip(ground[i],ground[j]))
        if sigma*tau*dx*dy<0:events.setdefault(Fraction(-tau*dy,sigma*dx),set()).update((i,j))
    events=sorted(events.items())
    order=sorted(range(n),key=lambda j:(tau*ground[j][1],sigma*ground[j][0],j))
    positions={j:i for i,j in enumerate(order)}
    first=events[0][0]/2 if events else Fraction(1)
    basis=oracle_greedy((sigma*first.numerator,tau*first.denominator))
    token={j:i for i,j in enumerate(sorted(basis))}
    token_paths=[[ground[j]] for j in sorted(basis)]
    local=defaultdict(lambda:defaultdict(lambda:[0,0]))
    token_moves=defaultdict(lambda:[0,0])
    L=(H-1).bit_length()
    block=lambda a:a[0].bit_length()
    initial_potential=sum(block(ground[j]) for j in basis)
    cross_moves=0;cross_drop=0;states=[];walls=[]
    def record(t):
        assert basis==oracle_greedy((sigma*t.numerator,tau*t.denominator)), 'direct independent greedy replay'
        moment=tuple(sum(ground[j][i] for j in basis) for i in range(2))
        if not states or states[-1]['basis']!=sorted(basis):
            if states:assert sigma*(moment[0]-states[-1]['moment'][0])<0 and tau*(moment[1]-states[-1]['moment'][1])>0
            states.append({'basis':sorted(basis),'moment':moment,'weight':[sigma*t.numerator,tau*t.denominator]})
    record(first)
    for k,(t,involved) in enumerate(events):
        groups={}
        for j in involved:
            c=sigma*ground[j][0]*t.numerator+tau*ground[j][1]*t.denominator
            groups.setdefault(c,[]).append(j)
        pivots=[]
        for group in groups.values():
            locations=sorted(positions[j] for j in group)
            assert locations==list(range(locations[0],locations[-1]+1))
            for loc,j in enumerate(sorted(group,key=lambda j:(sigma*ground[j][0],j)),locations[0]):
                while positions[j]>loc:
                    pos=positions[j];e,f=order[pos-1],j
                    if e in basis and f not in basis and independent(frozenset(basis-{e}|{f})):
                        a,b=ground[e],ground[f]
                        dx,dy=sigma*(a[0]-b[0]),tau*(b[1]-a[1])
                        assert dx>0 and dy>0 and Fraction(dy,dx)==t
                        old,new=block(a),block(b)
                        assert 0<=new<=L and 0<=old<=L and sigma*(old-new)>=0
                        i=token.pop(e);token[f]=i;token_paths[i].append(b)
                        if old!=new:
                            cross_moves+=1;cross_drop+=sigma*(old-new)
                        else:
                            assert old>0
                            local[old][t][0]+=dx;local[old][t][1]+=dy
                            token_moves[i,old][0]+=dx;token_moves[i,old][1]+=dy
                        basis.remove(e);basis.add(f)
                        pivots.append({'token':i,'from':e,'to':f,'same_block':old==new})
                    order[pos-1],order[pos]=f,e
                    positions[f],positions[e]=pos-1,pos
        if pivots:walls.append({'ratio':str(t),'pivots':pivots})
        nxt=events[k+1][0] if k+1<len(events) else t+2
        record((t+nxt)/2)
    assert len(states)==len(walls)+1
    assert sigma*(initial_potential-sum(block(ground[j]) for j in basis))==cross_drop
    assert cross_moves<=cross_drop<=rank*L
    controls=[]
    for label_bin,rows in sorted(local.items()):
        j=label_bin-1
        X,Y=2**j-1,H//2**j
        for i in range(rank):
            dx,dy=token_moves[i,label_bin]
            assert dx<=X and dy<=Y, 'each token traverses a dyadic block at most once'
        A=sum(v[0] for v in rows.values());B=sum(v[1] for v in rows.values());K=len(rows)
        assert A<=rank*X and B<=rank*Y and X*Y<=H
        primitive_pairs=[(a//gcd(a,b),b//gcd(a,b)) for a,b in rows.values()]
        assert len(set(primitive_pairs))==K
        small=[(a,b) for a,b in primitive_pairs if Fraction(a,A)+Fraction(b,B)<=Fraction(3,K)]
        assert 3*len(small)>=K
        assert len(small)<=(3*A//K)*(3*B//K)
        assert K**3<=27*A*B<=27*rank**2*H, 'distinct primitive directions versus two budgets'
        controls.append({'block':label_bin,'x_range':[2**j,2**(j+1)-1],
                         'walls':K,'total_dx':A,'total_dy':B,'token_x_cap':X,'token_y_cap':Y,
                         'short_primitive_pairs':small,
                         'directions':[{'ratio':str(t),'jump':v,'primitive_pair':[v[0]//gcd(*v),v[1]//gcd(*v)]}
                                       for t,v in sorted(rows.items())]})
    assert len(walls)<=cross_moves+sum(row['walls'] for row in controls)
    if rank==H:assert len(states)<=1+4*H*L
    return {'model':label,'H':H,'rank':rank,'ground':ground,'signs':signs,'states':states,'walls':walls,
            'cross_block_moves':cross_moves,'cross_block_potential_drop':cross_drop,'dyadic_blocks':controls,
            'token_paths':token_paths,'all_intervals_replayed':True}


def matrix_oracle(columns,r,p):
    @lru_cache(None)
    def independent(ids):
        if len(ids)>r:return False
        pivots={}
        for j in sorted(ids):
            v=list(columns[j])
            for i,row in sorted(pivots.items()):
                c=v[i]
                if c:v=[(a-c*b)%p for a,b in zip(v,row)]
            i=next((i for i,a in enumerate(v) if a),None)
            if i is None:return False
            pivots[i]=[a*pow(v[i],-1,p)%p for a in v]
        return True
    return independent


def basic_fixtures():
    fixtures=[]
    def add(ground,r,independent,H,label):
        rows=[sweep(ground,r,independent,H,label,signs) for signs in product((-1,1),repeat=2)]
        if len(ground)<=8:
            profiles={tuple(sum(ground[i][a] for i in basis) for a in range(2))
                      for basis in combinations(range(len(ground)),r) if independent(frozenset(basis))}
            actual={tuple(state['moment']) for f in rows for state in f['states']}
            assert actual==set(convex_hull(profiles)), 'independent full basis-profile hull'
        fixtures.extend(rows)
    for H,r in ((1,0),(1,1),(2,1),(2,2),(4,0),(4,4),(4,8),(7,7),(12,12),(16,8)):
        ground=tuple(sorted(hyperbolic_cross(H,2)))
        add(ground,r,lambda s:len(s)<=r,H,'uniform')
        p=101;random=Random(H*100+r)
        columns=[tuple(int(i==j) for i in range(r)) if j<r else tuple(random.randrange(p) for i in range(r))
                 for j in range(len(ground))]
        add(ground,r,matrix_oracle(columns,r,p),H,'linear_mod_101')
    ground=((0,0),(0,0),(1,0),(1,0),(2,1),(2,1),(3,0),(0,2))
    columns=[tuple(int(i==j) for i in range(4)) if j<4 else tuple((j+2)**i%7 for i in range(4))
             for j in range(len(ground))]
    add(ground,4,matrix_oracle(columns,4,7),6,'repeated_integer_weights')
    H=8;ground=((0,0),(1,0),(2,0),(3,0),(0,1),(1,1),(0,2),(0,3))
    pairs=[frozenset((2*i,2*i+1)) for i in range(4)]
    forbidden={a|b for a,b in combinations(pairs,2) if a|b!=pairs[2]|pairs[3]}
    bases={frozenset(s) for s in combinations(range(8),4)}-forbidden
    assert all(any(B-{e}|{f} in bases for f in C-B) for B in bases for C in bases for e in B-C)
    independent=lambda s:len(s)<4 or len(s)==4 and s not in forbidden
    add(ground,4,independent,H,'sparse_paving_basis_exchange_verified')
    for q in (4,8):
        H=2*q**3
        ground=tuple((q+i,(q-i)**2) for i in range(q))
        assert all((a+1)*(b+1)<=H for a,b in ground)
        add(ground,1,lambda s:len(s)<=1,H,'single_block_convex_chain')
    H=6;ground=tuple(sorted(hyperbolic_cross(H,2)))
    reference=((0,0),(1,0),(2,0),(0,1),(1,1),(3,0))
    columns=[]
    for a in ground:
        image=(3,0) if a==(2,1) else a
        columns.append(tuple(int(image==b) for b in reference))
    add(ground,6,matrix_oracle(columns,6,101),H,'actual_binomial_local_quotient')
    for p,points,label in ((23,[(pow(2,j,23),pow(2,7*j,23)) for j in range(7)],'seven_cyclic_points'),
                           (29,[(0,0),(1,1),(2,2),(3,4),(5,7),(11,13),(17,10)],'published_seven_point_control')):
        H=7;ground=tuple(sorted(hyperbolic_cross(H,2)))
        columns=[tuple(pow(x,a,p)*pow(y,b,p)%p for x,y in points) for a,b in ground]
        add(ground,7,matrix_oracle(columns,7,p),H,label)
    return fixtures


def convex_hull(points):
    points=sorted(set(points))
    if len(points)<=1:return points
    cross=lambda a,b,c:(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    chains=[]
    for seq in (points,points[::-1]):
        chain=[]
        for a in seq:
            while len(chain)>=2 and cross(chain[-2],chain[-1],a)<=0:chain.pop()
            chain.append(a)
        chains.extend(chain[:-1])
    return chains


def rank_two_controls():
    # A rank-two matroid is a partition into loops and parallel classes,
    # truncated at rank two. Restricted-growth labels enumerate all of them.
    ground=((1,16),(2,8),(3,5),(4,3),(5,2),(8,1))
    H=40
    def labels(prefix,maximum):
        if len(prefix)==len(ground):
            if maximum>=2:yield tuple(prefix)
        else:
            for c in range(maximum+2):yield from labels(prefix+[c],max(maximum,c))
    ledger=[];examples={};maximum=0
    for colors in labels([],0):
        independent=lambda s:len(s)<=2 and all(colors[i] for i in s) and len({colors[i] for i in s})==len(s)
        exposed=set();signed=[]
        for signs in product((-1,1),repeat=2):
            f=sweep(ground,2,independent,H,'all_rank_two_matroids',signs)
            exposed.update(tuple(s['moment']) for s in f['states'])
            signed.append(f)
        profiles={tuple(ground[i][a]+ground[j][a] for a in range(2))
                  for i,j in combinations(range(len(ground)),2) if independent(frozenset((i,j)))}
        hull=set(convex_hull(profiles))
        assert exposed==hull,'independent enumeration of every basis and its full convex hull'
        maximum=max(maximum,len(hull))
        ledger.append({'parallel_and_loop_labels':colors,'full_polygon_vertices':sorted(hull),
                       'quadrant_counts':[len(f['states']) for f in signed]})
        if any(f['dyadic_blocks'] for f in signed) and 'internal_exchange' not in examples:
            examples['internal_exchange']={'labels':colors,'quadrants':signed}
    assert len(ledger)==813,'all labelled rank-two matroids on six elements'
    return {'ground':ground,'hyperbolic_bound':H,'matroids_checked':len(ledger),'quadrant_sweeps':4*len(ledger),
            'full_basis_profile_hulls_independently_checked':True,'maximum_polygon_vertices':maximum,
            'enumeration_digest':hashlib.sha256(json.dumps(ledger,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
            'exact_examples':examples}


def finite_point_controls():
    ledger=[];examples={};maxima={};subsets=0;states=0
    for p in (2,3):
        universe=tuple(product(range(p),repeat=2))
        for mask in range(1,1<<len(universe)):
            points=[a for i,a in enumerate(universe) if mask>>i&1]
            n=len(points);ground=tuple(sorted(hyperbolic_cross(n,2)))
            columns=[tuple(pow(x,a,p)*pow(y,b,p)%p for x,y in points) for a,b in ground]
            f=sweep(ground,n,matrix_oracle(columns,n,p),n,'exhaustive_prime_field_points')
            for state in f['states']:
                standard={ground[j] for j in state['basis']}
                assert all(not a or (a-1,b) in standard for a,b in standard)
                assert all(not b or (a,b-1) in standard for a,b in standard)
            count=len(f['states']);subsets+=1;states+=count
            key=f'{p}:{n}'
            if count>maxima.get(key,0):
                maxima[key]=count;examples[key]={'characteristic':p,'points':points,'control':f}
            ledger.append([p,mask,n,count,len(f['walls']),f['cross_block_moves'],
                           [[b['block'],b['walls'],b['total_dx'],b['total_dy']] for b in f['dyadic_blocks']]])
    assert subsets==526,'every nonempty subset of F_2^2 and F_3^2'
    return {'point_subsets_checked':subsets,'initial_bases_checked':states,'maximum_counts_by_field_and_length':maxima,
            'ledger_sha256':hashlib.sha256(json.dumps(ledger,separators=(',',':')).encode()).hexdigest(),
            'exact_maximum_witnesses':examples}


def main():
    if not __debug__:raise RuntimeError('Run this verifier without -O: its exact assertions are required.')
    fixtures=basic_fixtures()
    rank_two=rank_two_controls()
    finite_points=finite_point_controls()
    report={'schema_version':1,'status':'passed',
            'command':'python3 -B scripts/run_planar_dyadic_bound_checks.py',
            'uniform_proof':'docs/planar_dyadic_bound.md',
            'positive_quadrant_bound':'1+ceil(log2(H))*(r+3*(r^2*H)^(1/3))',
            'planar_colength_bound':'g(I)<=1+4*n*ceil(log2(n)) over every field',
            'planar_asymptotic_order_proved':'Theta(n*log(n+1)), using the known generic corner-cut lower bound',
            'fixtures':fixtures,'all_rank_two_matroids':rank_two,'exhaustive_finite_point_ideals':finite_points,
            'spatial_goal_complete':False,'novelty_established':False,'goal_complete':False,
            'scope':'Exact checks of the uniform dyadic exchange proof, including independent greedy replay, primitive integer-direction counting, full profile hulls, all rank-two matroids on six elements, and all nonempty point subsets of two small prime-field planes. A proof is supplied for every finite planar ideal; current novelty and the three-variable endpoint remain unestablished.',
            'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                             for name in ('scripts/run_planar_dyadic_bound_checks.py','scripts/run_structural_complexity_checks.py')}}
    out=ROOT/'results/research/planar_dyadic_bound_checks.json'
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','basic_fixtures':len(fixtures),'basic_states':sum(len(f['states']) for f in fixtures),
                      'rank_two_matroids':rank_two['matroids_checked'],'rank_two_quadrants':rank_two['quadrant_sweeps'],
                      'point_subsets':finite_points['point_subsets_checked'],
                      'point_initial_bases':finite_points['initial_bases_checked'],
                      'goal_complete':False,'output':str(out.relative_to(ROOT))}))


if __name__=='__main__':main()
