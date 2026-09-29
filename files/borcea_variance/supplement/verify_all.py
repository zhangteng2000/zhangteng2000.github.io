#!/usr/bin/env python3
"""Run both exact finite-degree checkers and assert the full finite coverage.

This wrapper does not verify the human mathematical reductions or the
infinite analytic tail. Read the accompanying borcea_variance.pdf, including Appendices A and B.
"""
from __future__ import annotations
import argparse, json, subprocess, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def main()->None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bits',type=int,choices=(100,160,256),default=100)
    ap.add_argument('--log-dir',type=Path,default=ROOT/'rerun_logs')
    args=ap.parse_args()
    args.log_dir.mkdir(parents=True,exist_ok=True)
    entries=[('low','verify_low.py','certificate_n4_to_n13.json',4,13),
             ('mid','verify_mid.py','certificate_n14_to_n100000.json',14,100000)]
    reports=[]
    for name,checker,cert,nmin,nmax in entries:
        target=(args.log_dir/f'{name}_{args.bits}bit.json').resolve()
        subprocess.run([sys.executable,str(ROOT/checker),str(ROOT/cert),
                        '--bits',str(args.bits),'--log-json',str(target)],check=True)
        data=json.loads(target.read_text(encoding='utf-8'))
        if data.get('status')!='PASS' or (data['minimum_degree'],data['maximum_degree'])!=(nmin,nmax):
            raise RuntimeError('Expected degree range was not certified.')
        reports.append(data)
    if reports[0]['maximum_degree']+1!=reports[1]['minimum_degree']:
        raise RuntimeError('Degree gap between certificates.')
    result={'status':'PASS','precision_bits':args.bits,'minimum_degree':4,'maximum_degree':100000,
            'degree_blocks':sum(x['blocks'] for x in reports),'leaf_boxes':sum(x['leaves'] for x in reports),
            'scope':'Exact finite inequalities and coverage; human lemmas and infinite tail are documented separately.'}
    (args.log_dir/f'all_{args.bits}bit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
