#!/usr/bin/env python3
"""Exact bounded-shift lattices and complete matching-axis outer-hook fans."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product, permutations, combinations
from random import Random
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from run_hook_axis_lattice_checks import degree,subtract,multiply
from run_signed_face_checks import field_basis,field_reduce
from run_hook_axis_birth_checks import assignment

def price(label,theta,S):
 f,q=label;return q+(theta if f==1 else S if f==2 else 0)

def sources(F,H,A,t):
 P={(1,q):{q:1} for q in range(A)}
 P.update({(0,q):{q+j:c for j,c in F.items() if q+j<A} for q in range(t,A)})
 P.update({(2,q):{q+j:c for j,c in H.items() if q+j<A} for q in range(A)})
 heads=[{i+j:c for j,c in F.items() if i+j<A} for i in range(t)]
 return P,heads

def rank(V,p):return len(field_basis(tuple(range(len(V))),V.__getitem__,p)[0])

def births(P,heads,theta,S,p):
 prices={e:price(e,theta,S) for e in P};out=[]
 for i,h in enumerate(heads):
  prev=heads[:i]
  if rank(prev+[h],p)==rank(prev,p):out.append(None);continue
  for rho in sorted(set(prices.values())):
   V=prev+[v for e,v in P.items() if prices[e]<=rho]
   if rank(V+[h],p)==rank(V,p):out.append(rho);break
  else:out.append('infinity')
 return out

def full_lattice(F,H,A,theta,S,p):
 rows=[({A:1},{},{}),(dict(F),{0:1},{}),(dict(H),{},{0:1})]
 shifts=(theta,Q(0),S);steps=0
 def lp(R):return max(range(3),key=lambda j:(shifts[j]+degree(R[j]),j))
 while len({lp(R) for R in rows})<3:
  i,j=next((i,j) for i in range(3) for j in range(i) if lp(rows[i])==lp(rows[j]))
  k=lp(rows[i])
  if degree(rows[i][k])<degree(rows[j][k]):i,j=j,i
  d=degree(rows[i][k])-degree(rows[j][k]);a=rows[i][k][degree(rows[i][k])];b=rows[j][k][degree(rows[j][k])]
  c=a*pow(b,-1,p)%p if p else Q(a)/b
  rows[i]=tuple(subtract(x,y,d,c,p) for x,y in zip(rows[i],rows[j]));steps+=1
  assert steps<64*(A+1)
 rows.sort(key=lp);assert sum(degree(R[i]) for i,R in enumerate(rows))==A
 assert_determinant(rows,A,p)
 for u,b,c in rows:
  relation=subtract(subtract(u,multiply(F,b,p),0,1,p),multiply(H,c,p),0,1,p)
  assert all(q>=A for q in relation)
 COUNTS['exact_three_stream_lattice_determinants_and_congruences']+=1
 return rows,[shifts[i]+degree(R[i]) for i,R in enumerate(rows)]

def short_gamma(rows,costs,t,p):
 V={};prices={}
 for j,R in enumerate(rows):
  for k in range(t):
   V[j,k]={q+k:c for q,c in R[1].items() if q+k<t};prices[j,k]=costs[j]+k
 result=[]
 for i in range(t):
  prev=[{j:1} for j in range(i)]
  for rho in sorted(set(prices.values())):
   C=prev+[v for e,v in V.items() if prices[e]<=rho]
   if rank(C+[{i:1}],p)==rank(C,p):result.append(rho);break
  else:raise AssertionError('The full lattice must span every low b constraint.')
 return result

def factor_controls():
 rng=Random(20261023188);counts=0
 for p,(A,t) in product([0,2,3,101],[(4,2),(7,3),(12,4)]):
  for mode in ['dense','sparse','nonunit','zero']:
   def W():
    return {} if mode=='zero' else {i:c for i in range(t if mode=='nonunit' else 0,A) if (mode!='sparse' or rng.randrange(3)==0) and(c:=rng.randrange(p or 7))}
   F,H=W(),W();P,heads=sources(F,H,A,t)
   for theta,S in [(Q(-3,2),Q(7,3)),(Q(4,3),Q(-5,3)),(Q(A,2)+Q(1,3),Q(A)-Q(1,3))]:
    rows,costs=full_lattice(F,H,A,theta,S,p);gamma=short_gamma(rows,costs,t,p);original=births(P,heads,theta,S,p)
    for i,(g,a) in enumerate(zip(gamma,original)):
     assert g==(max(i,a) if a is not None else i),(p,A,t,F,H,theta,S,i,g,a)
     counts+=1;COUNTS['independent_original_vs_short_lattice_births']+=1
    # Independent perturbations of both shifts are allowed within E.
    E=Q(t);tt,ss=theta+Q(t,2),S-Q(t,3)
    original=births(P,heads,tt,ss,p)
    for i,a in enumerate(original):
     g=max(i,a) if a is not None else i
     assert any(d-E<=g<=d+t-1+E for d in costs),(p,A,t,F,H,theta,S,i,g,costs)
     counts+=1;COUNTS['independent_perturbed_shift_birth_band_checks']+=1
   print('three-stream lattice passed',p,A,t,mode,flush=True)
 print('TOTAL',counts,flush=True)
 return counts


from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
from itertools import product
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from run_hook_axis_birth_checks import assignment
from run_hook_axis_compilation_checks import hull,centroid,fixed_profile
from run_hook_bridge_checks import model,actual,normals,polygon
from run_hook_rank_departure_checks import complete_fan
from run_oracle_enumeration_checks import exact_adjacent_weight
from run_signed_face_checks import field_basis,field_reduce
from run_structural_complexity_checks import dot,area2

def label(f,q,t,a):
 return (1,t+q,0) if f==1 else (0,t+q,t) if f==0 else (0,a+q,0)

def prepare(A,t,operators,p):
 return {'A':A,'t':t,'a':A+t,'n':4*A*t+2*t*t,'p':p,'operators':operators,
         'sources':[sources(F,H,A,t) for F,H in operators],'cache':[{},{}]}

def factor(C,arm,theta,S):
 P,heads=C['sources'][arm];p=C['p'];order=tuple(sorted(P,key=lambda e:(price(e,theta,S),e)))
 if any(price(e,theta,S)==price(f,theta,S) for e,f in zip(order,order[1:])):
  raise ValueError('Price tie.')
 if order not in C['cache'][arm]:
  B,piv=field_basis(order,P.__getitem__,p);columns={e:{} for e in B}
  for i,h in enumerate(heads):
   rem,rep=field_reduce(h,piv,p);assert not rem
   for e,c in rep.items():
    if c:columns[e][i]=c
  birth=assignment(columns,{e:price(e,theta,S) for e in B},p)
  C['cache'][arm][order]=(B,birth)
 return order,*C['cache'][arm][order]

def decision(C,w):
 t,A,a=C['t'],C['A'],C['a'];x,y,z=map(Q,w)
 result=[];N=set();profile=list(fixed_profile(t,a));orders=[]
 for arm,theta,S in [(0,(x-t*z)/y,A-t*z/y),(1,(x-t*y)/z,A-t*y/z)]:
  order,B,birth=factor(C,arm,theta,S);orders.append(order);result.append(birth)
  for e,f in zip(order,order[1:]):
   le,lf=label(*e,t,a),label(*f,t,a)
   normal=tuple(v-u for u,v in zip(le,lf))
   if arm:normal=(normal[0],normal[2],normal[1])
   assert dot(w,normal)>0;N.add(normal)
  for f,q in B:
   for j in range(t):
    e=label(f,q,t,a);e=(e[0],e[1],e[2]+j)
    if arm:e=(e[0],e[2],e[1])
    profile=[v+u for v,u in zip(profile,e)]
 def gain(arm,i,e):
  f,q=e
  v=(1,q-i,-t) if f==1 else (0,q-i,0) if f==0 else (0,A+q-i,-t)
  return (v[0],v[2],v[1]) if arm else v
 winners=[]
 for i,j in product(range(t),repeat=2):
  options={(0,0,0)}
  if i in result[0]:options.add(gain(0,i,result[0][i]))
  if j in result[1]:options.add(gain(1,j,result[1][j]))
  winner=max(options,key=lambda q:dot(w,q));winners.append(winner)
  for option in options:
   d=tuple(u-v for u,v in zip(winner,option))
   if any(d):
    if not dot(w,d):raise ValueError('Head maximum tie.')
    N.add(d)
  profile=[v-u for v,u in zip(profile,winner)]
 assert all(dot(w,q)>0 for q in N)
 return (tuple(orders),tuple(winners)),tuple(profile),N

def original(C):
 t,a,p=C['t'],C['a'],C['p'];data=model(t,a);co=[]
 for g,y,z in data[2]:
  if g==0:co.append(C['operators'][0][1].get(y-t,0) if y>=t and z==0 else 0)
  elif g==2:co.append(C['operators'][1][1].get(z-t,0) if z>=t and y==0 else 0)
  elif y>=t:co.append(C['operators'][0][0].get(y-t,0) if z==0 else 0)
  else:co.append(C['operators'][1][0].get(z-t,0) if y==0 else 0)
 return actual(data,co,p,True)

def walk(C,coarse=None,verbose=False):
 n=C['n'];L=2*n+1;queue=[((1,L,L*L),None)];cells={};areas=defaultdict(Q);vertices=defaultdict(list)
 while queue:
  w,prev=queue.pop();key,profile,N=decision(C,w)
  if coarse:
   ckey,extra,inventory=coarse(C,w);key=(key,ckey);N|=extra
  assert key!=prev
  if key in cells:continue
  P=polygon(N);cells[key]=(profile,P);areas[profile]+=area2(P);vertices[profile].extend(P)
  rays=[integral(v) for v in P]
  def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
  triple=next(T for T in combinations(rays,3) if dot(T[0],cross(T[1],T[2])))
  sample=tuple(sum(v[k] for v in triple) for k in range(3))
  assert min(sample)>0 and max(sample)<=2048*C['a']**2 and all(dot(sample,q)>0 for q in N)
  assert decision(C,sample)[1]==profile
  COUNTS['whole_refinement_patches_and_small_integer_samples']+=1
  for v,u in zip(P,P[1:]+P[:1]):
   if any(v[k]==u[k]==0 for k in range(3)):continue
   normal=next(q for q in N if dot(q,v)==dot(q,u)==0)
   if coarse:
    from run_hook_rank_departure_checks import plane
    assert plane(normal) in inventory, (ckey,normal)
    COUNTS['complete_refinement_facets_in_local_inventory']+=1
   nw,_=exact_adjacent_weight(n,v,u,normal);queue.append((nw,key))
  if verbose and len(cells)%250==0:print('live fan',len(cells),len(vertices),flush=True)
 assert sum(areas.values())==1
 fan={v:hull(P) for v,P in vertices.items()}
 assert all(area2(fan[v])==a for v,a in areas.items())
 COUNTS['exact_area_one_covers_and_convex_profile_merges']+=1
 return fan,cells

def fan_controls():
 cases=[]
 for A,t,ops,p in [
  (4,2,[({0:1,2:2},{1:1,3:1}),({0:2,1:1},{0:1,2:2})],0),
  (5,2,[({0:1,1:1,2:1},{0:1,3:1}),({1:1,3:1},{2:1,4:1})],2),
  (5,3,[({3:1},{0:1,2:1}),({4:1},{1:1})],3),
  (6,2,[({0:1,3:2,5:1},{0:1,2:3}),({0:1,1:2,4:3},{0:1,5:4})],101),
  (5,2,[({0:1,1:2,3:1},{}),({2:1,4:1},{})],3),
  (4,2,[({},{0:1,2:1}),({},{1:1,3:1})],2),
 ]:
  C=prepare(A,t,ops,p)
  fan,cells=walk(C,coarse,verbose=True);cols=original(C)
  reference=complete_fan(cols,C['n'],p);translation=fixed_profile(t,A+t)
  ref={tuple(translation[k]+sum(e[k] for e in B) for k in range(3)):hull(P) for B,P in reference.items()}
  assert set(ref)==set(fan) and all(set(ref[v])==set(fan[v]) for v in fan)
  for profile,P in cells.values():
   B=field_basis(sorted(cols,key=lambda e:(dot(centroid(P),e),e)),cols.__getitem__,p)[0]
   assert profile==tuple(translation[k]+sum(e[k] for e in B) for k in range(3))
   assert all(dot(q,v)>=0 for q in normals(B,cols,p) for v in P)
  COUNTS['independent_complete_original_quotient_fans']+=1
  result={'A':A,'t':t,'characteristic':p,'length':C['n'],'refinement_patches':len(cells),
          'actual_initial_ideal_cones':len(fan),'original_fan_equal':True,'exact_total_area':'1',
          'coarse_windows_visited':len(C['window_cache']),
          'operators':[[{str(q):str(c) for q,c in sorted(W.items())} for W in pair] for pair in ops]}
  cases.append(result);print('full original outer scalar fan passed',result,flush=True)
 return cases


from fractions import Fraction as Q
from itertools import product
from collections import Counter
from math import gcd,lcm
from run_hook_axis_lattice_checks import lattice
from run_hook_rank_departure_checks import plane

COUNTS=Counter()

def assert_determinant(rows,B,p):
 out={}
 for perm in permutations(range(len(rows))):
  term={0:1}
  for i,j in enumerate(perm):term=multiply(term,rows[i][j],p)
  sign=(-1)**sum(perm[i]>perm[j] for i in range(len(perm)) for j in range(i+1,len(perm)))
  out=subtract(out,term,0,-sign,p)
 assert set(out)=={B} and out[B] in (1,-1,p-1)
 COUNTS['exact_polynomial_lattice_determinants']+=1

def integral(v):
 d=lcm(*(Q(c).denominator for c in v));v=tuple(int(Q(c)*d) for c in v);g=gcd(*v)
 assert g;return tuple(c//g for c in v)

def reduce_rows(rows,shifts,p,B):
 rows=list(rows)
 def lp(R):return max(range(len(shifts)),key=lambda j:(shifts[j]+degree(R[j]),j))
 steps=0
 while len({lp(R) for R in rows})<len(rows):
  i,j=next((i,j) for i in range(len(rows)) for j in range(i) if lp(rows[i])==lp(rows[j]))
  k=lp(rows[i])
  if degree(rows[i][k])<degree(rows[j][k]):i,j=j,i
  delta=degree(rows[i][k])-degree(rows[j][k]);a=rows[i][k][degree(rows[i][k])];b=rows[j][k][degree(rows[j][k])]
  c=a*pow(b,-1,p)%p if p else Q(a)/b
  rows[i]=tuple(subtract(x,y,delta,c,p) for x,y in zip(rows[i],rows[j]));steps+=1
  assert steps<128*(B+1+int(max(shifts)-min(shifts)))
 rows.sort(key=lp)
 assert sum(degree(R[i]) for i,R in enumerate(rows))==B
 return rows,[shifts[i]+degree(R[i]) for i,R in enumerate(rows)]

def inverse(W,B,p):
 c=pow(W[0],-1,p) if p else Q(1)/W[0];out={0:c}
 for k in range(1,B):
  value=-c*sum(W.get(j,0)*out.get(k-j,0) for j in range(1,k+1))
  if p:value%=p
  if value:out[k]=value
 assert {q:c for q,c in multiply(W,out,p).items() if q<B}=={0:1}
 return out

def lower_lattice(F,H,A,S,p):
 nu=min(min(F,default=A),min(H,default=A));B=A-nu
 if not B:return None,[],B
 F={q-nu:c for q,c in F.items()};H={q-nu:c for q,c in H.items()}
 if F.get(0):
  W={q:-c for q,c in multiply(inverse(F,B,p),H,p).items() if q<B}
  if p:W={q:c%p for q,c in W.items() if c%p}
  rows=[({B:1},{}),(W,{0:1})]
 else:
  W={q:-c for q,c in multiply(inverse(H,B,p),F,p).items() if q<B}
  if p:W={q:c%p for q,c in W.items() if c%p}
  rows=[({0:1},W),({},{B:1})]
 rows,costs=reduce_rows(rows,(Q(0),S),p,B)
 assert_determinant(rows,B,p)
 for b,c in rows:
  relation=subtract(multiply(F,b,p),multiply(H,c,p),0,-1,p)
  assert all(q>=B for q in relation)
 COUNTS['exact_lower_kernel_lattice_congruences']+=1
 return rows,costs,B

def bands(costs,t,A,ranges):
 out=set()
 for d in costs:
  for f,(lo,hi) in ranges.items():
   low=(d-hi-3*t).__floor__();high=(d-lo+3*t).__ceil__()
   out|={(f,e) for e in range(max(-t+1,low),min(A-1,high)+1)}
 return out

def upper_tail_labels(F,H,A,t,Slo,Shi,p):
 if not F:return set()
 _,costs,_=lower_lattice(F,H,A,Q(Slo)+Q(2,3),p)
 out=bands(costs,t,A,{0:(0,0),2:(Slo,Shi)})
 nu=min(F)
 out|={(1,e) for e in range(max(-t+1,nu-t+1),min(A-1,nu+t-1)+1)}
 return out

def slot(v,t,lower=None,upper=None):
 lo=t*int((v/t).__floor__());hi=lo+t
 if lower is not None:lo=max(lo,lower)
 if upper is not None:hi=min(hi,upper)
 assert lo<v<hi,(v,lo,hi)
 return lo,hi

def coarse(C,w):
 t,A,p=C['t'],C['A'],C['p'];x,y,z=map(Q,w);swap=z>y
 if swap:y,z=z,y
 assert y>z
 ops=list(reversed(C['operators'])) if swap else C['operators']
 rho=z/y;theta=(x-t*z)/y;phi=(x-t*y)/z;S=A-t*y/z;r0=Q(t,4*A)
 N={(0,1,-1)};key=[swap];typeY=typeZ=None;Yrange=Zrange=Srange=None
 def interval(L,D,lo,hi):
  if lo is not None:N.add(integral(tuple(a-lo*b for a,b in zip(L,D))))
  if hi is not None:N.add(integral(tuple(hi*b-a for a,b in zip(L,D))))
 def ywindow():
  nonlocal Yrange,typeY
  if theta>2*A:
   N.add((1,-2*A,-t));typeY='high';key.append('Yhigh')
  else:
   Yrange=slot(theta,t,-t,2*A);interval((1,0,-t),(0,1,0),*Yrange);typeY='full';key.extend(['Y',*Yrange])
 if rho<r0:
  N.add((0,t,-4*A))
  if phi<-A:
   N.add((-1,t,-A));key.append('inactive');typeY=typeZ='inactive'
  else:
   N.add((1,-t,A));key.append('contract');ywindow()
   if phi>2*A:
    N.add((1,-t,-2*A));typeZ='contract_high';key.append('Zhigh')
   else:
    Zrange=slot(phi,t,-A,2*A);interval((1,-t,0),(0,0,1),*Zrange);typeZ='contract_full';key.extend(['Z',*Zrange])
 else:
  N.add((0,-t,4*A));ywindow()
  Srange=slot(S,t,-3*A,A);interval((0,-t,A),(0,0,1),*Srange);key.extend(['S',*Srange])
  if phi>2*A:
   N.add((1,-t,-2*A));typeZ='high';key.append('Zhigh')
  else:
   N.add((-1,t,2*A));Zrange=slot(phi,t,-4*A,2*A);interval((1,-t,0),(0,0,1),*Zrange);typeZ='full';key.extend(['Z',*Zrange])
   lo=r0
   while 2*lo<rho:lo*=2
   hi=min(2*lo,Q(1));interval((0,0,1),(0,1,0),lo,hi);key.extend(['dyadic',str(lo),str(hi)])
 key=tuple(key);cache=C.setdefault('window_cache',{})
 if key not in cache:
  inventory={plane(q) for q in N}
  def pencils(arm,theta_range,S_range,full=True):
   D=(0,1,0) if not arm else(0,0,1)
   LT=(1,0,-t) if not arm else(1,-t,0);LS=(0,A,-t) if not arm else(0,-t,A)
   pairs=[]
   if theta_range is not None:pairs.append((LT,*theta_range))
   if S_range is not None:pairs.append((LS,*S_range))
   if full and theta_range is not None and S_range is not None:
    pairs.append((tuple(a-b for a,b in zip(LT,LS)),theta_range[0]-S_range[1],theta_range[1]-S_range[0]))
   for L,lo,hi in pairs:
    for k in range(int(Q(lo).__floor__())-1,int(Q(hi).__ceil__())+2):
     raw=tuple(a-k*b for a,b in zip(L,D))
     if any(raw):inventory.add(plane(raw))
  if typeY=='inactive':
   for arm in range(2):pencils(arm,(-2*A,2*A),(-2*A,2*A))
   labels=[set(),set()]
  else:
   labels=[]
   F,H=ops[0]
   if typeY=='full':
    _,costs=full_lattice(F,H,A,Q(Yrange[0])+Q(1,3),Q(A)-Q(1,3),p)
    labels.append(bands(costs,t,A,{1:Yrange,0:(0,0),2:(A-t,A)}))
    pencils(0,Yrange,(A-t,A))
   else:
    labels.append(upper_tail_labels(F,H,A,t,A-t,A,p));pencils(0,None,(A-t,A))
   F,H=ops[1]
   if typeZ=='full':
    _,costs=full_lattice(F,H,A,Q(Zrange[0])+Q(1,3),Q(Srange[0])+Q(2,3),p)
    labels.append(bands(costs,t,A,{1:Zrange,0:(0,0),2:Srange}));pencils(1,Zrange,Srange)
   elif typeZ=='high':
    labels.append(upper_tail_labels(F,H,A,t,*Srange,p));pencils(1,None,Srange)
   else:
    B=min(H,default=A);W={q:c for q,c in F.items() if q<B}
    if not B or not W:labels.append(set())
    elif typeZ=='contract_high':
     nu=min(W);labels.append({(1,e) for e in range(max(-t+1,nu-t+1),min(A-1,nu+t-1)+1)})
    else:
     _,r=lattice(W,B,Q(Zrange[0])+Q(1,3),p)
     costs=[Q(Zrange[0])+Q(1,3)+r,B-r]
     labels.append(bands(costs,t,A,{1:Zrange,0:(0,0)}))
    if Zrange is not None:pencils(1,Zrange,None,False)
   gradients=[]
   for arm,entries in enumerate(labels):
    V=set()
    for f,e in entries:
     v=(1,e,-t) if f==1 else(0,e,0) if f==0 else(0,A+e,-t)
     if arm:v=(v[0],v[2],v[1])
     V.add(v)
    gradients.append(V)
   for v in gradients[0]|gradients[1]:
    if any(v):inventory.add(plane(v))
   for u,v in product(*gradients):
    q=tuple(a-b for a,b in zip(u,v))
    if any(q):inventory.add(plane(q))
  cache[key]=(labels,inventory)
  if typeY!='inactive':
   assert all(len(L)<=128*t for L in labels)
   assert len(inventory)<=32768*(t+1)**2
  assert all(max(map(abs,q))<=8*A+8 for q in inventory)
  COUNTS['coarse_windows_with_complete_local_inventory']+=1
  COUNTS['maximum_positive_label_candidates_per_arm']=max(COUNTS['maximum_positive_label_candidates_per_arm'],*(len(L) for L in labels))
  COUNTS['coarse_region_'+str(typeZ)]+=1
 labels,inventory=cache[key]
 if typeY!='inactive':
  # Fresh actual birth labels check membership in the whole-window inventory.
  for arm,th,ss,scale in [(0,theta,A-t*rho,y),(1,phi,S,z)]:
   actual_arm=1-arm if swap else arm
   _,_,birth=factor(C,actual_arm,th,ss)
   for i,e in birth.items():
    if price(e,th,ss)>i:
     assert (e[0],e[1]-i) in labels[arm],(key,arm,th,ss,i,e,labels[arm])
     COUNTS['positive_actual_births_in_window_bands']+=1
 if swap:
  N={(q[0],q[2],q[1]) for q in N};inventory={plane((q[0],q[2],q[1])) for q in inventory}
 return key,N,inventory


def lower_controls():
 rng=Random(20261023189);fixtures=[]
 for p,(A,t),mode in product((0,2,3,101),((5,2),(8,3)),('dense','sparse','nonunit')):
  def W():
   return {q:c for q in range(A-t if mode=='nonunit' else 0,A)
           if (mode!='sparse' or rng.randrange(3)==0) and(c:=rng.randrange(p or 7))}
  F,H=W(),W();P,heads=sources(F,H,A,t);P={e:v for e,v in P.items() if e[0]!=1}
  for S in (Q(-3*A+2,3),Q(2,3),Q(A)-Q(1,3)):
   rows,costs,B=lower_lattice(F,H,A,S,p)
   if not B:continue
   teff=min(t,B);small=[];prices=[]
   for row,d in zip(rows,costs):
    for k in range(teff):
     small.append({q+k:c for q,c in row[0].items() if q+k<teff});prices.append(d+k)
   original=births(P,heads,Q(1,3),S,p)
   for i in range(teff):
    prev=[{j:1} for j in range(i)];gamma=None
    for rho in sorted(set(prices)):
     V=prev+[v for v,c in zip(small,prices) if c<=rho]
     if rank(V+[{i:1}],p)==rank(V,p):gamma=rho;break
    alpha=original[i]
    assert gamma==(None if alpha=='infinity' else i if alpha is None else max(i,alpha)),(p,A,t,F,H,S,i,gamma,alpha)
    COUNTS['independent_lower_only_kernel_birth_checks']+=1
  fixtures.append({'A':A,'t':t,'characteristic':p,'mode':mode})
 return fixtures


def large_window_controls():
 rng=Random(20261023190);cases=[]
 for A,t,p,spec in [(96,2,2,[(0,48,1,70),(12,60,20,80)]),
                    (64,3,1009,[(12,40,17,46),(0,35,20,52)])]:
  def W(lo,hi):return {q:c for q in range(lo,hi+1) if(c:=rng.randrange(1,p))}
  ops=[(W(f0,f1),W(h0,h1)) for f0,f1,h0,h1 in spec];C=prepare(A,t,ops,p);narrow=0;queries=0
  ss=[Q(1,3)+Q(1,971),t+Q(1,7)+Q(1,971),Q(A,3)+Q(1,11)+Q(1,971),A+Q(1,13)+Q(1,971)]
  rr=[Q(t,7*A)+Q(1,997*A),Q(t,3*A)+Q(1,991*A),Q(2,7)+Q(1,983*A),Q(5,8)+Q(1,977*A)]
  for s,rho in product(ss,rr):
   key,_,_=coarse(C,(s,Q(1),rho));labels=C['window_cache'][key][0]
   if 'inactive' not in key:narrow+=sum(len(L)<2*A for L in labels)
   queries+=1
  assert queries==16 and narrow>0
  COUNTS['larger_exact_weight_window_queries']+=queries
  COUNTS['larger_nonvacuous_window_exclusions']+=narrow
  cases.append({'A':A,'t':t,'characteristic':p,'exact_weight_queries':queries,
                'nonvacuous_window_exclusions':narrow,'coarse_windows':len(C['window_cache']),
                'scope':'Exact original factor elimination and window membership; no full original quotient fan at this size.'})
 return cases


def dyadic_budget_controls():
 out=[]
 for A,t in [(4,2),(16,2),(64,3),(128,2),(256,8),(1024,2)]:
  r=Q(t,4*A);bands_count=0;maximum=0
  def bins(lo,hi):return int((hi/t).__ceil__()-(lo/t).__floor__()+2)
  while r<1:
   top=min(2*r,Q(1));nz=bins(A-Q(t)/r,A-Q(t)/top)
   np=bins(Q(-4*A),Q(2*A));ny=bins(Q(-t),t+2*A*top)
   count=nz*np*ny
   assert count<=1024*Q(A*A,t*t)
   maximum=max(maximum,count);bands_count+=1;r=top
   COUNTS['independent_dyadic_grid_product_budgets']+=1
  out.append({'A':A,'t':t,'dyadic_bands':bands_count,
              'maximum_potential_bin_triples_in_one_band':maximum,
              'uniform_test_ceiling':'1024*(A/t)^2'})
 return out


def main():
 if not __debug__:raise RuntimeError('Assertions must remain enabled.')
 prior_path='results/research/hook_axis_lattice_checks.json'
 prior=json.loads((ROOT/prior_path).read_text());assert prior['status']=='passed'
 for f,h in prior['source_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
 factor_checks=factor_controls();lower=lower_controls();large=large_window_controls()
 budget=dyadic_budget_controls();fans=fan_controls()
 files=set(prior['source_sha256'])|{prior_path,'docs/hook_scalar_outer_lattice_fans.md',
                                  'scripts/run_hook_scalar_outer_lattice_checks.py'}
 report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_scalar_outer_lattice_checks.py',
  'uniform_result':{'setting':'A=a-t>=t>=2, n=4at-2t^2, arbitrary scalar central axes F,G and arbitrary scalar matching outer axes H,K, over every field.',
    'fan_and_complete_refinement_bound':'O(A^2*t^2*log(2+A/t))=O(n^2*log(2+A/t)).',
    'bounded_shift_lattice_lemma':'A low q^t affine constraint and a predictable-degree basis give gamma in the union of [d_j-E,d_j+t-1+E] whenever every shift changes by at most E.',
    'weight_list':'Coefficient-dependent complete positive integer list with O(n^2*log(2+A/t)) members, height O(a^2).',
    'limits':'Opposite-arm outer images and transverse coefficient dependence are outside this theorem.'},
  'counts':dict(COUNTS),'factor_birth_comparisons':factor_checks,'lower_kernel_fixtures':lower,'complete_original_fans':fans,
  'larger_window_controls':large,'dyadic_grid_budget_controls':budget,
  'limits':['Polynomial relation modules, weak Popov reduction and predictable degree are classical.',
    'Exact finite original fans support the uniform proof; they do not establish an asymptotic lower excess.',
    'No arbitrary transverse/opposite-outer theorem, unrestricted G_3 endpoint, exact current novelty, optimal bit complexity or smooth-number speed claim.',
    'No coefficient-independent common weight list or profile formula is asserted.'],
  'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
 (ROOT/'results/research/hook_scalar_outer_lattice_checks.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({'status':'passed','counts':dict(COUNTS),'source_hashes':len(files)}),flush=True)


if __name__=='__main__':main()
