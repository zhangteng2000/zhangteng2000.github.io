#!/usr/bin/env python3
"""Propose a finite cover; acceptance requires verify_mid.py.

Generation uses NumPy floating point. Its output is untrusted until the
independent integer checker validates the entire partition and every bound.
"""
import argparse, json, time
from pathlib import Path
from propose_mid import build


def generate(nmin=14,nmax=100000):
    if not 14 <= nmin <= nmax:
        raise ValueError('This driver supports 14 <= nmin <= nmax.')
    blocks=[];reports=[];m=nmin-1;stop=nmax-1;t0=time.time()
    while m<=stop:
        if m<30:u=m
        elif m<100:u=min(stop,m+max(1,m//12))
        else:u=min(stop,m+max(1,m//4))
        while True:
            ref=16 if m==13 else 4
            block,report=build((m,u),maxnodes=300000,ref=ref,progress=0)
            if block is not None:break
            if u==m:raise RuntimeError(f'UNRESOLVED degree {m+1}: {report}')
            u=(m+u)//2
        blocks.append(block);reports.append({'block':[m,u],**report})
        print(f'PROPOSED n={m+1}..{u+1}: leaves={report["leaves"]}',flush=True)
        m=u+1
    result={'schema':'borcea-p2-box-cover-v2',
            'warning':'Only the listed degrees are certified. No claim for degrees 4 through 13. Validate using verify_exact.py.',
            'minimum_degree':nmin,'maximum_degree':nmax,'blocks':blocks}
    return result,{'elapsed_seconds':time.time()-t0,'blocks':reports}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--nmin',type=int,default=14)
    ap.add_argument('--nmax',type=int,default=100000)
    ap.add_argument('--output',type=Path,default=Path('regenerated_certificate.json'))
    args=ap.parse_args()
    result,report=generate(args.nmin,args.nmax)
    args.output.write_text(json.dumps(result,separators=(',',':')),encoding='utf-8')
    Path(str(args.output)+'.generation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('PROPOSAL WRITTEN:',args.output,flush=True)
