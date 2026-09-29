#!/usr/bin/env python3
"""Regenerate an untrusted low-degree proposal, then compare its bytes.

NumPy is required only for this optional proposal stage. The exact checkers
use the standard library alone and do not trust the proposal generator.
"""
from pathlib import Path
import argparse, hashlib, json
from generate_low import build

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=Path('regenerated_low.json'))
    args=ap.parse_args()
    blocks=[]
    for m in range(3,13):
        block,stats=build((m,m),ref=32,progress=0)
        if block is None:raise RuntimeError(f'UNRESOLVED n={m+1}: {stats}')
        blocks.append(block)
        print(f'PROPOSED n={m+1}: leaves={len(block["leaves"])}',flush=True)
    raw=json.dumps({'blocks':blocks},separators=(',',':')).encode('utf-8')
    args.output.write_bytes(raw)
    expected=Path(__file__).resolve().parent/'certificate_n4_to_n13.json'
    report={'sha256':hashlib.sha256(raw).hexdigest(),'matches_supplied_bytes':raw==expected.read_bytes()}
    print(json.dumps(report,indent=2),flush=True)
    if not report['matches_supplied_bytes']:
        print('A changed floating-point proposal is not automatically invalid. It must be checked independently.')
if __name__=='__main__':main()
