#!/usr/bin/env python3
"""Fixed-prime smooth-number rank/unrank toolkit with k<=8 interval auditing.

Main operations:
  unrank(P,N): exponent vector for the 1-based N-th P-smooth number
  value(P,exps): exact integer value
  audit(P,N): fast unrank + independent interval-rank certificate for P subset {2,3,5,7,11,13,17,19}, k<=8

The fast high-k kernels use floating log search and exact final-band sorting. The audit step
certifies a returned exponent vector by proving count_le(value) == N using integer log intervals
and exact big-integer resolution of boundary ambiguities.
"""
from __future__ import annotations
import argparse, json, os, re, subprocess, sys
from pathlib import Path
from typing import Iterable, List, Optional, Dict, Any

ROOT = Path(__file__).resolve().parent
BIN = ROOT / "bin"
EXPS_RE = re.compile(r"exps=\[([^\]]*)\]")
SECONDS_RE = re.compile(r"(?:^|\s)seconds=([0-9]+(?:\.[0-9]+)?(?:[eE][-+]?[0-9]+)?)")
DIGITS_RE = re.compile(r"digits=([0-9]+)")
CERT_RE = re.compile(r"rank_certified=(true|false)")
COUNT_RE = re.compile(r"certified_count_le=([0-9]+)")
CERT_SECONDS_RE = re.compile(r"cert_seconds=([0-9]+(?:\.[0-9]+)?(?:[eE][-+]?[0-9]+)?)")
GROUPA_RE = re.compile(r"groupA=([0-9]+)")
GROUPB_RE = re.compile(r"groupB=([0-9]+)")
AMB_RE = re.compile(r"ambiguous_possible=([0-9]+)")

SUPPORTED_AUDIT_PRIMES = {2,3,5,7,11,13,17,19}

def parse_primes(x: str | Iterable[int]) -> List[int]:
    if isinstance(x, str): ps = [int(t) for t in x.replace(';', ',').split(',') if t.strip()]
    else: ps = [int(t) for t in x]
    if not ps: raise ValueError('empty prime list')
    return sorted(ps)

def parse_exps(x: str | Iterable[int]) -> List[int]:
    if isinstance(x, str):
        s=x.strip()
        if s.startswith('[') and s.endswith(']'): s=s[1:-1]
        return [int(t) for t in s.replace(';', ',').split(',') if t.strip()]
    return [int(t) for t in x]

def csv(v: Iterable[int]) -> str: return ','.join(str(t) for t in v)

def value(primes: Iterable[int], exps: Iterable[int]) -> int:
    ps=parse_primes(primes); es=parse_exps(exps)
    if len(ps)!=len(es): raise ValueError('prime/exponent length mismatch')
    out=1
    for p,e in zip(ps,es):
        if e<0: raise ValueError('negative exponent')
        out *= p**e
    return out

def ensure_built():
    req=[
        BIN/'smooth_3prime_beatty_ranker',
        BIN/'smooth_layer_compressed_general',
        BIN/'smooth_sums_only_scalable',
        BIN/'smooth_interval_audit_exps',
        BIN/'smooth_interval_audit_exps_k6',
    ]
    if all(p.exists() for p in req): return
    subprocess.run([str(ROOT/'build.sh')], cwd=str(ROOT), check=True)

def run_cmd(args: List[str], timeout: Optional[int]=None) -> str:
    ensure_built()
    p=subprocess.run(args, cwd=str(ROOT), text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
    if p.returncode!=0:
        raise RuntimeError('command failed: '+ ' '.join(args)+'\nSTDOUT:\n'+p.stdout+'\nSTDERR:\n'+p.stderr)
    return p.stdout.strip()

def parse_unrank_output(out: str) -> Dict[str,Any]:
    m=EXPS_RE.search(out)
    if not m: raise RuntimeError('could not parse exps from output: '+out)
    exps=[] if not m.group(1).strip() else [int(t) for t in m.group(1).split(',')]
    sec=SECONDS_RE.search(out); dig=DIGITS_RE.search(out)
    return {'exps':exps, 'seconds':float(sec.group(1)) if sec else None, 'digits':int(dig.group(1)) if dig else None, 'raw':out}

def unrank(primes: Iterable[int], n: int, method: str='auto', timeout: Optional[int]=None) -> Dict[str,Any]:
    ps=parse_primes(primes)
    if n < 1: raise ValueError('N is 1-based')
    if n == 1:
        return {'method':'trivial','primes':ps,'N':n,'exps':[0]*len(ps),'digits':1,'raw':''}
    if len(ps)==1:
        exps=[n-1]
        return {'method':'trivial','primes':ps,'N':n,'exps':exps,'digits':len(str(value(ps,exps))),'raw':''}
    if method=='auto':
        if len(ps)==3: method='beatty3'
        elif len(ps)==5: method='layer5'
        else: method='sums'
    if method=='beatty3':
        out=run_cmd([str(BIN/'smooth_3prime_beatty_ranker'), 'nth', csv(ps), str(n)], timeout)
    elif method=='layer5':
        out=run_cmd([str(BIN/'smooth_layer_compressed_general'), 'nth', csv(ps), str(n)], timeout)
    elif method=='sums':
        out=run_cmd([str(BIN/'smooth_sums_only_scalable'), 'nth', csv(ps), str(n)], timeout)
    else:
        raise ValueError('unknown method '+method)
    d=parse_unrank_output(out); d.update({'method':method,'primes':ps,'N':n}); return d

def audit(primes: Iterable[int], n: int, timeout: Optional[int]=None) -> Dict[str,Any]:
    ps=parse_primes(primes)
    if len(ps)>8 or any(p not in SUPPORTED_AUDIT_PRIMES for p in ps):
        raise ValueError('audit supports k<=8 and primes drawn from {2,3,5,7,11,13,17,19}')
    fast=unrank(ps,n,timeout=timeout)
    out=run_cmd([str(BIN/'smooth_interval_audit_exps'), csv(ps), csv(fast['exps']), str(n)], timeout)
    cert=CERT_RE.search(out); cnt=COUNT_RE.search(out); csec=CERT_SECONDS_RE.search(out)
    ga=GROUPA_RE.search(out); gb=GROUPB_RE.search(out); amb=AMB_RE.search(out)
    return {
        'method':'fast-unrank-plus-interval-audit-k8',
        'primes':ps,'N':n,'exps':fast['exps'],'fast_method':fast['method'],
        'fast_seconds':fast.get('seconds'),'digits':fast.get('digits'),
        'rank_certified': cert.group(1)=='true' if cert else False,
        'certified_count_le': int(cnt.group(1)) if cnt else None,
        'cert_seconds': float(csec.group(1)) if csec else None,
        'groupA': int(ga.group(1)) if ga else None,
        'groupB': int(gb.group(1)) if gb else None,
        'ambiguous_possible': int(amb.group(1)) if amb else None,
        'audit_raw': out,
        'fast_raw': fast['raw'],
    }


def _heap_enumerate_until(primes: Iterable[int], limit_value: int | None = None, limit_rank: int | None = None):
    """Exact heap enumeration helper for small/medium values or ranks.

    This is deliberately simple and independent from the fast kernels. It is used for
    smoke tests and convenience helpers, not for the large-rank headline algorithms.
    """
    import heapq
    ps = parse_primes(primes)
    if limit_value is None and limit_rank is None:
        raise ValueError('provide limit_value or limit_rank')
    heap = [(1, tuple([0] * len(ps)))]
    seen = {1}
    out = []
    while heap:
        v, exps = heapq.heappop(heap)
        if limit_value is not None and v > limit_value:
            if limit_rank is None:
                return out, (v, exps)
        out.append((v, exps))
        if limit_rank is not None and len(out) >= limit_rank:
            return out, None
        for i, prime in enumerate(ps):
            nv = v * prime
            if nv not in seen:
                seen.add(nv)
                ne = list(exps)
                ne[i] += 1
                heapq.heappush(heap, (nv, tuple(ne)))
    return out, None

def count_le(primes: Iterable[int], x: int) -> int:
    if x < 1:
        return 0
    vals, _ = _heap_enumerate_until(primes, limit_value=x)
    return len(vals)

def prev_leq(primes: Iterable[int], x: int) -> Dict[str, Any]:
    if x < 1:
        raise ValueError('no positive P-smooth number <= x')
    vals, _ = _heap_enumerate_until(primes, limit_value=x)
    if not vals:
        raise ValueError('no positive P-smooth number <= x')
    v, exps = vals[-1]
    return {'primes': parse_primes(primes), 'x': x, 'value': v, 'exps': list(exps), 'rank': len(vals)}

def next_geq(primes: Iterable[int], x: int) -> Dict[str, Any]:
    if x <= 1:
        return {'primes': parse_primes(primes), 'x': x, 'value': 1, 'exps': [0]*len(parse_primes(primes)), 'rank': 1}
    vals, first_over = _heap_enumerate_until(primes, limit_value=x-1)
    if first_over is None:
        raise RuntimeError('internal enumeration failed')
    v, exps = first_over
    return {'primes': parse_primes(primes), 'x': x, 'value': v, 'exps': list(exps), 'rank': len(vals)+1}


def main():
    ap=argparse.ArgumentParser()
    sub=ap.add_subparsers(dest='cmd', required=True)
    p=sub.add_parser('unrank'); p.add_argument('primes'); p.add_argument('N', type=int); p.add_argument('--method', default='auto'); p.add_argument('--timeout', type=int)
    p=sub.add_parser('audit'); p.add_argument('primes'); p.add_argument('N', type=int); p.add_argument('--timeout', type=int, default=300)
    p=sub.add_parser('value'); p.add_argument('primes'); p.add_argument('exps')
    p=sub.add_parser('count-le'); p.add_argument('primes'); p.add_argument('x', type=int)
    p=sub.add_parser('prev-leq'); p.add_argument('primes'); p.add_argument('x', type=int)
    p=sub.add_parser('next-geq'); p.add_argument('primes'); p.add_argument('x', type=int)
    args=ap.parse_args()
    if args.cmd=='unrank': print(json.dumps(unrank(args.primes,args.N,args.method,args.timeout), indent=2))
    elif args.cmd=='audit': print(json.dumps(audit(args.primes,args.N,args.timeout), indent=2))
    elif args.cmd=='value': print(value(args.primes,args.exps))
    elif args.cmd=='count-le': print(count_le(args.primes,args.x))
    elif args.cmd=='prev-leq': print(json.dumps(prev_leq(args.primes,args.x), indent=2))
    elif args.cmd=='next-geq': print(json.dumps(next_geq(args.primes,args.x), indent=2))
if __name__=='__main__': main()
