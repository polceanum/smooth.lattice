#!/usr/bin/env python3
"""Exact full death curves, clipped variation and joint bridge-state budgets."""
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
import ast
import hashlib
import json

from run_hook_bridge_checks import model,actual,COUNTS as REPLAY
from run_hook_filtered_annihilator_checks import assemble
from run_hook_barrier_sector_checks import flag_states,state_frontier,complete_sector_control,vertical_boundary_control
from run_hook_weighted_partner_checks import own_sectors,add_plain,reduce_plain
from run_signed_face_checks import field_basis,field_reduce
from run_structural_complexity_checks import dot

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def system(cols,t,E,J,p):
    eps,Y,Z=E;u=(J[0]+J[1])/2 if J[1] else J[0]+1
    own=lambda v:(v[1]>=t and v[2]<t) if v[0] else (v[1]>=2*t and t<=v[2]<2*t)
    LY=tuple(v for v in cols if v[0]==eps and own(v) and v[1]*u+v[2]<Y*u+Z)
    H=tuple(v for v in cols if not v[0] and t<=v[1]<2*t and t<=v[2]<2*t)
    fixed=(*LY,*H)
    _,piv=field_basis(fixed,cols.__getitem__,p)
    values={v:field_reduce(c,piv,p)[0] for v,c in cols.items()}
    C=tuple(v for v in cols if v[0]==eps and
            ((eps and v[1]<t) or (not eps and t<=v[1]<2*t and v[2]>=2*t)))
    D=tuple(v for v in cols if v[0]!=eps and v not in H)
    return values,C,D,fixed


def merge(intervals):
    out=[]
    for a,b in sorted(intervals):
        if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(out[-1][1],b))
        else:out.append((a,b))
    return out


def is_z(v,t,eps):
    return (t<=v[1]<2*t and v[2]>=2*t) if eps else v[1]<t


def death_column(values,C,D,E,p,u):
    eps,Y,Z=E;cp={}
    for v in C:
        if v[1]*u+v[2]<Y*u+Z:add_plain(values[v],cp,p)
    if not reduce_plain(values[E],cp,p):return None,'loop'
    for v in sorted(D,key=lambda v:(v[1]*u+v[2],v)):
        if add_plain(values[v],cp,p) and not reduce_plain(values[E],cp,p):return v,'finite'
    return None,'unreachable'


def original_pair(cols,fixed,E,F,u,p):
    lam=(E[1]-F[1])*u+E[2]-F[2]
    if E[0]:lam=-lam
    w=(lam,u,Q(1));tau=dot(w,E)
    assert tau==dot(w,F) and all(dot(w,v)<tau for v in fixed)
    cheap=tuple(v for v in cols if dot(w,v)<tau)
    B,pp=field_basis(cheap,cols.__getitem__,p)
    assert field_reduce(cols[E],pp,p)[0] and field_reduce(cols[F],pp,p)[0]
    _,pf=field_basis((*B,F),cols.__getitem__,p)
    assert not field_reduce(cols[E],pf,p)[0]
    COUNTS['independent_original_cheaper_pair_rank_certificates']+=1


def curve(values,C,D,t,a,E,J,p,original=None,predicted=None):
    eps,Y,Z=E;alpha=Y-2*t+1;U=Q(a,alpha)
    L=J[0];R=min(U,J[1]) if J[1] else U
    if L>=R:return None
    entries={Q(v[2]-Z,Y-v[1]) for v in C}
    cuts={L,R}|{q for q in entries if L<q<R}
    DZ=tuple(v for v in D if is_z(v,t,eps));zcuts=set()
    for i,v in enumerate(D):
        for w in D[:i]:
            if v[1]!=w[1]:
                q=Q(w[2]-v[2],v[1]-w[1])
                if L<q<R:cuts.add(q)
    for i,v in enumerate(DZ):
        for w in DZ[:i]:
            if v[1]!=w[1]:
                q=Q(w[2]-v[2],v[1]-w[1])
                if L<=q<=R:zcuts.add(q)
    cuts=sorted(cuts);pieces=[];M=Q(0);partners=set();previous=None
    first=None;jumps=Q(0);integral=Q(0);jump_witnesses=[];negative_own=[]
    for left,right in zip(cuts,cuts[1:]):
        F,kind=death_column(values,C,D,E,p,(left+right)/2)
        if kind=='finite':
            m=Q(Y-F[1]);b=Q(Z-F[2])
            if eps:m,b=-m,-b
            sub={left,right}
            for mm in (m,m-t):
                if mm:
                    q=-b/mm
                    if left<q<right:sub.add(q)
            sub=sorted(sub)
        else:sub=[left,right]
        for aa,bb in zip(sub,sub[1:]):
            u=(aa+bb)/2
            if predicted and kind=='finite':assert F[1]*u+F[2]==predicted(u)
            inside=kind=='finite' and 0<m*u+b<t*u
            if kind=='loop':ff=(Q(0),Q(0)) if eps else (Q(t),Q(0))
            elif kind=='unreachable':ff=(Q(-t),Q(0)) if eps else (Q(0),Q(0))
            elif m*u+b<=0:ff=(Q(0),Q(0))
            elif m*u+b>=t*u:ff=(Q(-t if eps else t),Q(0))
            else:ff=(-m,-b) if eps else (m,b)
            f0=ff[0]*aa+ff[1];f1=ff[0]*bb+ff[1]
            if first is None:first=f0
            if previous is not None:
                assert f0>=previous,(E,J,aa,f0,previous,kind,F)
                if f0>previous:
                    jumps+=f0-previous;jump_witnesses.append([str(aa),str(f0-previous)])
                    COUNTS['nonzero_upward_clipped_jumps']+=1
                COUNTS['independent_clipped_jump_sign_checks']+=1
            previous=f1;integral+=f1-f0
            active=inside and is_z(F,t,eps)
            assert ff[0]>=-2*t
            if active:
                assert ff[0]>=alpha
                M+=bb-aa;partners.add(F)
                if original:original_pair(original[0],original[1],E,F,u,p)
                COUNTS['full_curve_active_cross_arm_pieces']+=1
            elif inside and ff[0]<0:
                negative_own.append([str(aa),str(bb),F,str(ff[0])])
                COUNTS['negative_own_arm_clipped_pieces']+=1
            pieces.append({'lo':aa,'hi':bb,'partner':F,'kind':kind,'active':active,'slope':ff[0]})
            COUNTS['independent_full_opposite_threshold_cells']+=1
    assert previous-first==integral+jumps
    assert previous-first<=t*U and integral>=alpha*M-2*t*(R-L)
    assert alpha*M<=3*t*U and M<=Q(27*t*a,Y*Y)
    active=[v for v in pieces if v['active']]
    union=merge((v['lo'],v['hi']) for v in active)
    assert M==sum((bb-aa for aa,bb in union),Q(0))
    relevant={q for q in entries if any(aa<=q<=bb for aa,bb in union)}
    events=sum(len([q for q in entries if aa<=q<=bb]) for aa,bb in union)
    assert events<=t*Y*M+2*t*len(union)
    assert len(zcuts)<=t*t*U
    runs=0;old=None
    for v in active:
        if old is None or old['hi']<v['lo']:runs+=1
        elif old['partner']!=v['partner']:
            assert v['lo'] in entries|zcuts
            runs+=1
        old=v
    assert len(partners)<=runs<=len(union)+len(relevant)+len(zcuts)
    COUNTS['whole_curve_variation_and_event_budgets']+=1
    return {'ratio_sector':list(map(lambda v:None if v is None else str(v),J)),
            'curve_cells':len(pieces),'active_measure':str(M),'active_components':len(union),
            'active_partners':sorted(partners),'partner_runs':runs,
            'component_intervals':[[str(aa),str(bb)] for aa,bb in union],
            'clipped_integral':str(integral),'upward_jump_total':str(jumps),
            'clipped_endpoint_difference':str(previous-first),
            'negative_own_arm_witnesses':negative_own[:8],'upward_jump_witnesses':jump_witnesses[:8],
            'same_layer_entry_events':len(relevant),'opposite_Z_order_events':len(zcuts),
            'partner_bound':str(30*Q(a*t*t,Y)+4*t**3),
            '_M':M,'_union':union,'_partners':partners}


def hook_case(t,a,E,p,data,co,cols):
    assert E[1]>=(3+E[0])*t
    rows=[];allpartners=set()
    for J in own_sectors(t):
        values,C,D,fixed=system(cols,t,E,J,p)
        row=curve(values,C,D,t,a,E,J,p,(cols,fixed))
        if row is None:continue
        states=flag_states(data,cols,t,a,E,J,p);records=[]
        for state in states:
            for branch in (['upper_wedge'] if E[0] else ['lower_low','lower_middle']):
                records+=state_frontier(state,t,a,E,J,branch,p)
        intervals=sorted((v['lo'],v['hi']) for v in records)
        for (aa,bb),(cc,dd) in zip(intervals,intervals[1:]):assert bb<=cc
        union=merge(intervals)
        assert union==row['_union']
        assert {v['partner'] for v in records}==row['_partners']
        assert len(union)<=len(states)<=t*t+2
        assert len(row['_partners'])<=30*Q(a*t*t,E[1])+4*t**3
        allpartners.update(row['_partners'])
        for k in ('_M','_union','_partners'):del row[k]
        rows.append(row);COUNTS['independent_full_curve_vs_critical_states_equalities']+=1
    assert len(allpartners)<=30*Q(a*t**4,E[1])+4*t**5
    print('Clipped hook case passed',t,a,E,p,'partners',len(allpartners),flush=True)
    return {'t':t,'a':a,'field_characteristic':p,'heavy_label':E,'coefficients':co,
            'distinct_partners':len(allpartners),'sector_checks':rows}


def abstract_controls():
    # This represented quotient exercises the analytic lemma, not the hook
    # realization. Its threshold has an independent explicit min/max formula.
    t=3;a=30;E=(0,12,3)
    D=((1,13,0),(1,0,20),(1,12,0),(1,1,22))
    values={E:{0:1,1:1},D[0]:{0:1},D[1]:{0:1},D[2]:{1:1},D[3]:{1:1}}
    out=[]
    for added in (False,True):
        C=((0,4,15),) if added else ()
        vv=dict(values)
        if added:vv[C[0]]={0:1}
        def predicted(u):
            first=min(13*u,Q(20));second=min(12*u,u+22)
            return second if added and u>Q(3,2) else max(first,second)
        row=curve(vv,C,D,t,a,E,(Q(1),None),101,predicted=predicted)
        if not added:assert row['active_components']==2 and row['negative_own_arm_witnesses']
        else:assert Q(row['upward_jump_total'])>0
        for k in ('_M','_union','_partners'):del row[k]
        out.append({'kind':'abstract rank-two represented quotient','same_layer_entry_added':added,**row})
    COUNTS['independent_abstract_minmax_threshold_controls']+=2
    return out


def main():
    if not __debug__:raise RuntimeError('Assertions must remain enabled.')
    old=json.loads((ROOT/'results/research/hook_barrier_sector_checks.json').read_text())
    assert old['status']=='passed'
    for f,sha in old['source_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==sha
    cases=[]
    for c in old['cases']:
        t,a,p,E=c['t'],c['a'],c['field_characteristic'],tuple(c['heavy_label'])
        data=model(t,a);co=c['coefficients'];cols=actual(data,co,p,True)
        cases.append(hook_case(t,a,E,p,data,co,cols))
    t=3;a=16;p=101;E=(0,9,3);units=({(0,0):1,(1,1):1},{(0,0):1,(1,1):1})
    data,co,cols=assemble(t,a,units,p);cases.append(hook_case(t,a,E,p,data,co,cols))
    abstract=abstract_controls()
    sector=complete_sector_control();vertical=vertical_boundary_control()
    assert COUNTS['negative_own_arm_clipped_pieces'] and COUNTS['nonzero_upward_clipped_jumps']
    files={'docs/hook_clipped_threshold_budget.md','scripts/run_hook_clipped_threshold_checks.py',
           'docs/hook_barrier_sectors.md','docs/hook_weighted_partners.md','docs/cross_layer_direction_budget.md',
           'results/research/hook_barrier_sector_checks.json'}
    pending=['scripts/run_hook_clipped_threshold_checks.py']
    while pending:
        f=pending.pop()
        for node in ast.walk(ast.parse((ROOT/f).read_text())):
            names=[node.module] if isinstance(node,ast.ImportFrom) else [v.name for v in node.names] if isinstance(node,ast.Import) else []
            for name in names:
                target='scripts/'+str(name)+'.py'
                if str(name).startswith('run_') and (ROOT/target).exists() and target not in files:
                    files.add(target);pending.append(target)
    report={'schema_version':1,'status':'passed','command':'python3 -B scripts/run_hook_clipped_threshold_checks.py',
            'scope':'Joint clipped-price variation and distinct-partner budgets in the constant-axis, zero-outer-image low/middle wedge. The resulting restricted count improves its thickness dependence but does not bound the high-price complement or improve unrestricted G_3.',
            'counts':dict(sorted(COUNTS.items())),'original_quotient_replay_counts':dict(REPLAY),'cases':cases,
            'abstract_analytic_controls':abstract,'complete_original_sector_control':sector,
            'independent_vertical_boundary_control':vertical,
            'uniform_bounds':{'partner':'30*a*t^4/Y+4*t^5','active_measure_per_own_order_sector':'3*t*a/(Y-2*t+1)^2',
                              'restricted_fan':'O(n^2*t^4*(t+log(n+1))*log(n+1))',
                              'global_polynomial_excess_localization':'If g>=c*n^(2+epsilon) and t<=n^(epsilon/5-delta), the low wedge has o(n^(2+epsilon)) cones.'},
            'limits':['Only constant-axis central units with zero outer images are covered.',
                      'Only w_x<t*max(w_y,w_z) is counted; high-price aggregate control remains open.',
                      'Literal pair eligibility may recur; a separate parallel-shadow bound controls actual pieces.',
                      'Abstract controls exercise the analytic lemma and are not hook quotients.',
                      'Exact novelty, unrestricted sharp growth, polynomial lower separation and the final paper remain open.'],
            'source_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in sorted(files)}}
    (ROOT/'results/research/hook_clipped_threshold_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','counts':report['counts'],'source_hashes':len(files)}))


if __name__=='__main__':main()
