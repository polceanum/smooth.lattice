#!/usr/bin/env python3
"""Exact tensor-prefix, birth-curve and original axis-hook cone checks."""
from pathlib import Path
import ast,hashlib,json
from collections import Counter

from fractions import Fraction as Q
from random import Random
from itertools import product
from run_signed_face_checks import field_basis,field_reduce
from run_hook_bridge_checks import model,actual,normals,polygon
from run_hook_rank_departure_checks import rank,complete_fan,plane
from run_structural_complexity_checks import dot,area2

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()

def assignment(vectors,prices,p):
 order=sorted(vectors,key=lambda e:(-prices[e],e));B,piv=field_basis(order,vectors.__getitem__,p);index={e:i for i,e in enumerate(B)}
 return {i:max(rep,key=index.__getitem__) for i,(row,rep) in piv.items()}

def single(W,A,t,p,theta):
 cols={(1,q):{q:1} for q in range(A)}
 cols.update({(0,q):{q+j:c for j,c in W.items() if q+j<A and c} for q in range(t,A)})
 cost=lambda e:e[1]+theta*e[0]
 B,piv=field_basis(sorted(cols,key=lambda e:(cost(e),e)),cols.__getitem__,p);C={e:{} for e in B}
 for i in range(t):
  H={i+j:c for j,c in W.items() if i+j<A and c};rem,rep=field_reduce(H,piv,p);assert not rem
  for e,c in rep.items():
   if c:C[e][i]=c
 return B,C,assignment(C,{e:cost(e) for e in C},p),cols

def encoded(t,a,Ws,p,w):
 lam,u,z=map(Q,w);lam/=z;u/=z;w=(lam,u,Q(1));A=a-t
 factors=[single(Ws[0],A,t,p,(lam-t)/u),single(Ws[1],A,t,p,lam-t*u)]
 pure=[]
 for arm,(B,C,birth,cols) in enumerate(factors):
  for x,q in B:
   for j in range(t):pure.append((x,t+q,t*(1-x)+j) if arm==0 else(x,t*(1-x)+j,t+q))
 grad=tuple(sum(e[k] for e in pure) for k in range(3));N=[]
 for arm,(B,C,birth,cols) in enumerate(factors):
  def key(e):return e[1]*u+(lam-t)*e[0] if arm==0 else e[1]+(lam-t*u)*e[0]
  E=sorted(cols,key=lambda e:(key(e),e))
  for e,f in zip(E,E[1:]):
   dx=f[0]-e[0];dq=f[1]-e[1]
   q=(dx,dq,-t*dx) if arm==0 else(dx,-t*dx,dq)
   assert dot(w,q)>0;N.append(q)
 for i,j in product(range(t),repeat=2):
  options=[(0,0,0)]
  if i in factors[0][2]:
   x,q=factors[0][2][i];options.append((x,q-i,-t*x))
  if j in factors[1][2]:
   x,q=factors[1][2][j];options.append((x,-t*x,q-j))
  best=max(options,key=lambda x:dot(w,x));grad=tuple(x-y for x,y in zip(grad,best))
  for alt in options:
   delta=tuple(x-y for x,y in zip(best,alt))
   if any(delta):assert dot(w,delta)>0;N.append(delta)
 return grad,N,factors


def rr(vectors,p):
 return len(field_basis(tuple(range(len(vectors))),vectors.__getitem__,p)[0])


def factor(cols,head,theta,p):
 cost=lambda e:Q(e[1])+theta*e[0]
 B,piv=field_basis(sorted(cols,key=lambda e:(cost(e),e)),cols.__getitem__,p)
 C={e:{} for e in B}
 for i,H in enumerate(head):
  rem,rep=field_reduce(H,piv,p);assert not rem
  for e,c in rep.items():
   if c:C[e][i]=c
 assigned=assignment(C,{e:cost(e) for e in C},p)
 return B,C,assigned


def direct_births(cols,head,theta,p):
 prices={e:Q(e[1])+theta*e[0] for e in cols};out=[]
 for i,H in enumerate(head):
  previous=head[:i]
  if rr(previous+[H],p)==rr(previous,p):out.append(None);continue
  for rho in sorted(set(prices.values())):
   vectors=previous+[v for e,v in cols.items() if prices[e]<=rho]
   if rr(vectors+[H],p)==rr(vectors,p):out.append(rho);break
  else:raise AssertionError('Head vector not in the spanning pure family.')
 return out


def birth_tables(cols,head,A,t,p,direct=False):
 tables=[[] for _ in head];previous=None
 for k in range(-A-1,A+1):
  theta=Q(2*k+1,2);B,C,assigned=factor(cols,head,theta,p)
  values=[Q(assigned[i][1])+theta*assigned[i][0] if i in assigned else None for i in range(t)]
  if direct and k in [-A-1,-A//2,0,A//2,A]:
   assert values==direct_births(cols,head,theta,p)
   COUNTS['independent_entire_pure_prefix_birth_thresholds']+=t
  if previous:
   for i,(old,new) in enumerate(zip(previous,values)):
    if old is None:assert new is None;continue
    assert new is not None and 0<=new-old<=1
  for i in range(t):
   if i not in assigned:continue
   x,q=assigned[i];assert q>=t if not x else q>=0
   piece={'lo':None if k==-A-1 else Q(k),'hi':None if k==A else Q(k+1),'x':x,'q':q}
   if tables[i]:
    old=tables[i][-1];boundary=piece['lo']
    assert old['hi']==boundary and old['x']*boundary+old['q']==x*boundary+q
    if (old['x'],old['q'])==(x,q):old['hi']=piece['hi'];continue
   tables[i].append(piece)
  previous=values
 for i,T in enumerate(tables):
  assert len(T)<=2*A+2
  if T:
   at_t=next(v for v in T if (v['lo'] is None or v['lo']<=t) and (v['hi'] is None or t<=v['hi']))
   assert at_t['x']*t+at_t['q']>=t>i
  COUNTS['exact_all_parameter_birth_tables']+=1
 return tables


def cross3(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])


def segment(q,inequalities):
 points=set()
 for h in inequalities:
  v=cross3(q,h);s=sum(v)
  if not s:continue
  v=tuple(Q(x,s) for x in v)
  if all(dot(h,v)>=0 for h in inequalities):points.add(v)
 if len(points)<2:return None
 assert len(points)==2
 return tuple(sorted(points))


def on_segment(v,endpoints):
 return all(min(a,b)<=x<=max(a,b) for x,a,b in zip(v,*endpoints))


def curves(TY,TZ,A,t):
 out={};gamma=(1,-t,-t)
 for i,Y in enumerate(TY):
  for j,Z in enumerate(TZ):
   if not Y or not Z:continue
   pieces=[]
   for Yp,Zp,branch in product(Y,Z,[-1,1]):
    x,q=Yp['x'],Yp['q'];z,r=Zp['x'],Zp['q']
    a=(x,q-i,-t*x);b=(z,-t*z,r-j);normal=tuple(s-v for s,v in zip(a,b));assert any(normal)
    I=[(1,0,0),(0,1,0),(0,0,1),a,tuple(branch*v for v in gamma)]
    for P,arm in [(Yp,0),(Zp,1)]:
     if P['lo'] is not None:I.append((1,-P['lo'],-t) if not arm else(1,-t,-P['lo']))
     if P['hi'] is not None:I.append((-1,P['hi'],t) if not arm else(-1,t,P['hi']))
    endpoints=segment(normal,I)
    if endpoints is None:continue
    midpoint=tuple((s+v)/2 for s,v in zip(*endpoints))
    if min(midpoint)<=0 or dot(a,midpoint)<=0 or branch*dot(gamma,midpoint)<=0:continue
    pieces.append({'plane':plane(normal),'endpoints':endpoints,'Y':Yp,'Z':Zp,'branch':branch,'indices':(i,j)})
   assert len(pieces)<=16*(A+1),(i,j,len(pieces),A)
   # A common ratio strictly decreases; theta and phi increase together.
   for P in pieces:
    i,j=P['indices'];x,q=P['Y']['x'],P['Y']['q'];z,r=P['Z']['x'],P['Z']['q'];ky=q-i+x*t;kz=r-j+z*t
    assert ky>0 and kz>0
    w=tuple((a+b)/2 for a,b in zip(*P['endpoints']));lam,u=w[0]/w[2],w[1]/w[2];theta=(lam-t)/u;phi=lam-t*u
    assert x+Q(ky)/(theta-t)==z+Q(kz)/(phi-t)
    assert phi-t==Q(kz)*(theta-t)/(ky+(x-z)*(theta-t))
    COUNTS['independent_projective_birth_curve_pieces']+=1
   out[i,j]=pieces
 return out


def curve_contacts(curve_map,A,t):
 keys=list(curve_map);max_contacts=0;overlaps=0;isolated=0;gamma=(1,-t,-t)
 for pos,key in enumerate(keys):
  for other in keys[:pos]:
   contacts=set()
   for P,R in product(curve_map[key],curve_map[other]):
    if P['branch']!=R['branch']:continue
    v=cross3(P['plane'],R['plane']);i,j=P['indices'];k,l=R['indices']
    yi,yj,yk,yl=P['Y'],P['Z'],R['Y'],R['Z']
    ki=yi['q']-i+yi['x']*t;kj=yj['q']-j+yj['x']*t
    kk=yk['q']-k+yk['x']*t;kl=yl['q']-l+yl['x']*t
    aa=Q(yk['x']-yl['x'])-Q(kl*(yi['x']-yj['x']),kj)
    bb=Q(kk)-Q(kl*ki,kj)
    if not any(v):
     assert aa==bb==0
     common={v for v in (*P['endpoints'],*R['endpoints']) if on_segment(v,P['endpoints']) and on_segment(v,R['endpoints'])}
     if len(common)>1:overlaps+=1
     contacts.update(v for v in common if min(v)>0 and dot(gamma,v))
     continue
    if not sum(v):continue
    w=tuple(Q(x,sum(v)) for x in v)
    if min(w)<=0 or not dot(gamma,w) or not on_segment(w,P['endpoints']) or not on_segment(w,R['endpoints']):continue
    assert aa and bb
    theta=t-bb/aa;den=ki+(yi['x']-yj['x'])*(theta-t);assert den
    u=Q(kj)/den;lam=u*theta+t;total=lam+u+1
    replay=(lam/total,u/total,1/total)
    assert replay==w
    contacts.add(w);isolated+=1;COUNTS['independent_original_line_vs_ratio_contact_equalities']+=1
   assert len(contacts)<=64*(A+1),(key,other,len(contacts),A)
   max_contacts=max(max_contacts,len(contacts));COUNTS['complete_interaction_curve_pair_contact_budgets']+=1
 return {'curves':len(keys),'pieces':sum(map(len,curve_map.values())),'largest_pair_contacts':max_contacts,
         'noncoincident_segment_contacts_with_multiplicity':isolated,'coincident_segment_overlaps_with_multiplicity':overlaps}


def tensor_controls(rng):
 results=[]
 for t,p in [(2,0),(2,2),(3,3),(3,101),(4,2)]:
  C={k:{i:rng.randrange(p or 5) for i in range(t)} for k in range(t+3)}
  D={k:{i:rng.randrange(p or 5) for i in range(t)} for k in range(t+2)}
  C={k:{i:c for i,c in v.items() if c} for k,v in C.items()};D={k:{i:c for i,c in v.items() if c} for k,v in D.items()}
  if t==4:D={k:{} for k in D}
  pricesC={k:Q(5*k-7,3) for k in C};pricesD={k:Q(7*k-4,5) for k in D};u=Q(7,3)
  ac=assignment(C,pricesC,p);ad=assignment(D,pricesD,p)
  for vectors,prices,assigned in [(C,pricesC,ac),(D,pricesD,ad)]:
   for rho in sorted(set(prices.values())):
    F=[v for e,v in vectors.items() if prices[e]>=rho]
    for j in range(t+1):
     E=[{i:1} for i in range(t-j,t)]
     intersection=rr(F,p)+j-rr(F+E,p)
     assert intersection==sum(i>=t-j and prices[e]>=rho for i,e in assigned.items())
     COUNTS['independent_two_flag_intersection_rank_tables']+=1
  cols={};prices={}
  for k,j in product(C,range(t)):
   e=('Y',k,j);cols[e]={i*t+j:c for i,c in C[k].items()};prices[e]=pricesC[k]+j
  for i,k in product(range(t),D):
   e=('Z',i,k);cols[e]={i*t+j:c for j,c in D[k].items()};prices[e]=i*u+pricesD[k]
  for i,j in product(range(t),repeat=2):
   e=('H',i,j);cols[e]={i*t+j:1};prices[e]=i*u+j
  token=[]
  for i,j in product(range(t),repeat=2):
   options=[i*u+j]
   if i in ac:options.append(pricesC[ac[i]]+j)
   if j in ad:options.append(i*u+pricesD[ad[j]])
   token.append(max(options))
  for rho in sorted(set(prices.values())):
   E=tuple(e for e in cols if prices[e]>=rho)
   assert rank(E,cols,p)==sum(x>=rho for x in token)
   COUNTS['independent_entire_tensor_prefix_ranks']+=1
  B=field_basis(sorted(cols,key=lambda e:(-prices[e],e)),cols.__getitem__,p)[0]
  assert sum(prices[e] for e in B)==sum(token)
  results.append({'short_dimension':t,'characteristic':p,'rank':t*t,'prefixes':len(set(prices.values()))})
 return results


def replicated_controls(rng):
 results=[]
 for A,t,p in [(6,2,0),(7,3,2),(8,3,101)]:
  sources=[]
  for arm in range(2):
   pure={(1,q):{q:1} for q in range(A)}
   pure.update({(0,q):{i:c for i in range(A) if (c:=rng.randrange(p or 5))} for q in range(t,A)})
   head=[{i:c for i in range(A) if (c:=rng.randrange(p or 5))} for _ in range(t)]
   if arm==0:head[0]=dict(pure[0,t])
   if arm==1:head[1]=dict(head[0])
   sources.append((pure,head))
  cols={}
  for arm,(pure,head) in enumerate(sources):
   for (x,q),vector in pure.items():
    for j in range(t):
     e=(x,t+q,t*(1-x)+j) if not arm else(x,t*(1-x)+j,t+q)
     cols[e]={i*t+j:c for i,c in vector.items()} if not arm else{A*t+j*A+i:c for i,c in vector.items()}
  for i,j in product(range(t),repeat=2):
   cols[0,t+i,t+j]={**{r*t+j:c for r,c in sources[0][1][i].items()},
                    **{A*t+i*A+r:c for r,c in sources[1][1][j].items()}}
  checked=0
  for trial in range(18):
   w=(rng.randrange(1,10000),rng.randrange(1,1000),rng.randrange(1,1000))
   if len({dot(w,e) for e in cols})!=len(cols):continue
   lam,u=Q(w[0],w[2]),Q(w[1],w[2]);F=[factor(*sources[0],(lam-t)/u,p),factor(*sources[1],lam-t*u,p)]
   labels=[]
   for arm,(B,C,birth) in enumerate(F):
    labels.extend((x,t+q,t*(1-x)+j) if not arm else(x,t*(1-x)+j,t+q) for x,q in B for j in range(t))
   profile=tuple(sum(e[k] for e in labels) for k in range(3))
   for i,j in product(range(t),repeat=2):
    options=[(0,0,0)]
    if i in F[0][2]:
     x,q=F[0][2][i];options.append((x,q-i,-t*x))
    if j in F[1][2]:
     x,q=F[1][2][j];options.append((x,-t*x,q-j))
    best=max(options,key=lambda e:dot(w,e));profile=tuple(x-y for x,y in zip(profile,best))
   B=field_basis(sorted(cols,key=lambda e:(dot(w,e),e)),cols.__getitem__,p)[0]
   assert len(B)==2*A*t and profile==tuple(sum(e[k] for e in B) for k in range(3))
   COUNTS['independent_non_Toeplitz_original_replicated_profiles']+=1;checked+=1
  assert checked>=10
  results.append({'long_dimension':A,'short_dimension':t,'characteristic':p,'dependent_head_vectors':True,'profile_checks':checked})
 return results


def axis_columns(data,t,Ws,p):
 co=[]
 for g,y,z in data[2]:
  arm=0 if y>=t else 1;q,j=(y-t,z) if not arm else(z-t,y)
  co.append(Ws[arm].get(q,0) if g==1 and not j else 0)
 return actual(data,co,p,True)


def main():
 if not __debug__:raise RuntimeError('Assertions must remain enabled.')
 previous=json.loads((ROOT/'results/research/hook_inverse_orientation_checks.json').read_text());assert previous['status']=='passed'
 for f,h in previous['source_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h
 rng=Random(2026102348);tensor=tensor_controls(rng);replicated=replicated_controls(rng);cases=[];curve_reports=[]
 configs=[(2,8,101,'dense'),(3,12,1009,'dense'),(3,16,3,'sparse'),(3,16,0,'binomial'),
          (3,10,2,'nonunit'),(2,9,101,'zero'),(4,17,101,'sparse'),(3,10,3,'scalar')]
 for t,a,p,mode in configs:
  A=a-t;data=model(t,a)
  if mode=='dense':Ws=[{q:rng.randrange(1,p) for q in range(A)} for _ in range(2)]
  elif mode=='sparse':Ws=[{0:1,min(5,A-1):2,min(9,A-1):1},{0:1,min(4,A-1):1,min(10,A-1):2}]
  elif mode=='binomial':Ws=[{0:1,5:1},{0:1,6:1}]
  elif mode=='nonunit':Ws=[{1:1,4:1},{2:1,5:1}]
  elif mode=='zero':Ws=[{}, {0:1,3:1}]
  else:Ws=[{0:1},{0:1}]
  Ws=[{q:(v%p if p else Q(v)) for q,v in W.items() if q<A and (v%p if p else v)} for W in Ws]
  cols=axis_columns(data,t,Ws,p);checked=0
  for _ in range(22):
   w=(rng.randrange(1,10000),rng.randrange(1,1000),rng.randrange(1,1000))
   if len({dot(w,e) for e in cols})!=len(cols):continue
   grad,N,F=encoded(t,a,Ws,p,w)
   B=field_basis(sorted(cols,key=lambda e:(dot(w,e),e)),cols.__getitem__,p)[0]
   expected=tuple(sum(e[k] for e in B) for k in range(3));assert grad==expected
   orderN=N[:sum(len(f[3])-1 for f in F)]
   P=polygon(N);R=polygon(normals(B,cols,p)|set(orderN));assert area2(P)==area2(R)
   assert all(dot(q,v)>=0 for q in N for v in R)
   COUNTS['independent_original_minimum_profiles']+=1;COUNTS['independent_entire_original_flag_cone_patches']+=1;checked+=1
  assert checked>=10
  TY=[]
  for W in Ws:
   pure={(1,q):{q:1} for q in range(A)}
   pure.update({(0,q):{q+j:c for j,c in W.items() if q+j<A} for q in range(t,A)})
   head=[{i+j:c for j,c in W.items() if i+j<A} for i in range(t)]
   table=birth_tables(pure,head,A,t,p,direct=True);TY.append(table)
   # Krylov multiplication gives the additional normalized birth interlacing.
   for theta in [Q(k,2) for k in range(-2*A-1,2*A+2,2)]:
    B,C,assigned=factor(pure,head,theta,p)
    values=[Q(assigned[i][1])+theta*assigned[i][0]-i if i in assigned else None for i in range(t)]
    assert all(values[i+1] is None or (values[i] is not None and values[i]>=values[i+1]) for i in range(t-1))
    COUNTS['independent_Krylov_normalized_birth_interlacings']+=1
  contact=curve_contacts(curves(TY[0],TY[1],A,t),A,t);curve_reports.append({'t':t,'a':a,'characteristic':p,'mode':mode,**contact})
  fan=complete_fan(cols,4*a*t-2*t*t,p) if (t,a,mode) in [(2,8,'dense'),(2,9,'zero')] else None
  if fan is not None:
   profiles={tuple(sum(e[k] for e in B) for k in range(3)) for B in fan};assert len(profiles)==len(fan)
   COUNTS['complete_original_axis_fans']+=1
  cases.append({'t':t,'a':a,'characteristic':p,'mode':mode,'operators':[{str(q):str(v) for q,v in W.items()} for W in Ws],
                'original_profile_and_cone_checks':checked,'complete_original_cones':len(fan) if fan is not None else None,
                'birth_pieces':[[len(T) for T in factor_tables] for factor_tables in TY]})
  print('Axis birth checks passed',t,a,p,mode,'weights',checked,'curves',contact['curves'],'pieces',contact['pieces'],flush=True)
 # A genuinely non-Toeplitz factor tests the larger replicated-matrix statement.
 A,t,p=7,3,101
 pure={(1,q):{q:1} for q in range(A)}
 pure.update({(0,q):{i:rng.randrange(p) for i in range(A)} for q in range(t,A)})
 head=[{i:rng.randrange(p) for i in range(A)} for _ in range(t)]
 head[0]=dict(pure[0,t])
 tab=birth_tables(pure,head,A,t,p,direct=True)
 assert tab[0] and tab[0][-1]['x']==0
 abstract_contacts=curve_contacts(curves(tab,tab,A,t),A,t)
 files={'docs/hook_axis_birth_curves.md','scripts/run_hook_axis_birth_checks.py','docs/hook_bridge_reduction.md',
        'docs/square_zero_extension_fans.md','results/research/hook_inverse_orientation_checks.json'}
 pending=['scripts/run_hook_axis_birth_checks.py']
 while pending:
  f=pending.pop()
  for node in ast.walk(ast.parse((ROOT/f).read_text())):
   names=([node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else [])
   for name in names:
    target='scripts/'+str(name)+'.py'
    if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:files.add(target);pending.append(target)
 report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_axis_birth_checks.py',
 'scope':'All-field collective crossing and exact profile reduction for zero-outer axis-dependent hooks, including nonunits and zeros, and the larger replicated-matrix family. No unrestricted endpoint, sharpness or exact strongest-current novelty claim.',
 'uniform_result':{'collective_count':'O(A^2+A*t^4), A=a-t>=t>=2','combined_box_count':'O(n^(17/8)), n asymptotic A*t',
 'critical_aspect':'A>=c*t^2 implies O_c(n^2), with arbitrary coefficient degenerations and no prefix-regular or inverse-degree bound',
 'profile_formula':'P(lambda,u)-sum_(i,j) max(u*(alpha_i((lambda-t)/u)-i),beta_j(lambda-t*u)-j,0)',
 'curve_mechanism':'R_i=delta+kappa/(theta-t), kappa>0; each interaction is one increasing graph on either side. Pair contacts cost O(A), including coincident-arc endpoints.',
 'oracle':'O(A^3+t^2) field operations and O(A^2) storage for a generic minimum profile; no optimal basis-output or bit-complexity claim'},
 'counts':dict(sorted(COUNTS.items())),'tensor_prefix_cases':tensor,'non_Toeplitz_replicated_profile_cases':replicated,'original_axis_cases':cases,'curve_contacts':curve_reports,
 'non_Toeplitz_replicated_matrix_curve_control':abstract_contacts,
 'limits':['Both outer images vanish; transverse dependence is not removed by a coordinate change.',
 'The generic compatible hook family was already quadratically bounded up to logarithms; finite fan counts are calibration, not asymptotic lower evidence.',
 'Classical flag splitting, matroid greedy/dual facts, arrangements and the lattice-volume theorem are credited. Exact priority of the combined bound is unestablished.',
 'Sharp unrestricted G_3, infinite actual polynomial excess, the consequential unrestricted endpoint and the final paper remain incomplete.'],
 'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
 (ROOT/'results/research/hook_axis_birth_checks.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}),flush=True)


if __name__=='__main__':main()
