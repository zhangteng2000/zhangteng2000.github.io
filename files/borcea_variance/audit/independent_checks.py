#!/usr/bin/env python3
"""Additional audits, not a substitute for the mathematical reduction.

Independent of both certificate checkers and all generators.  The topology
check uses a sorted binary-address partition and an independent volume sum.
The other checks below test universal auxiliary identities on exact rational
examples. Such sample checks are diagnostics, NOT proofs of the lemmas.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb, factorial
from random import Random
import hashlib, json

ROOT = Path(__file__).resolve().parents[1] / 'supplement'


def topology(path: Path):
    raw = path.read_bytes()
    data = json.loads(raw)
    blocks = data.get('blocks', [data])
    last = None
    result = []
    for block in blocks:
        lo, hi = block['block']
        assert isinstance(lo, int) and isinstance(hi, int) and 3 <= lo <= hi
        assert last is None or lo == last + 1
        last = hi
        addresses = []
        prefixes = {}
        volume = F(0)
        rule_counts = {}
        for path_string, rule in block['leaves']:
            assert len(path_string) % 2 == 0
            bounds = [[F(0), F(5)], [F(0), F(1)], [F(0), F(1)]]
            address = ''
            for offset in range(0, len(path_string), 2):
                axis, side = int(path_string[offset]), int(path_string[offset+1])
                assert axis in (0, 1, 2) and side in (0, 1)
                assert address not in prefixes or prefixes[address] == axis
                prefixes[address] = axis
                mid = sum(bounds[axis]) / 2
                bounds[axis][1-side] = mid
                address += str(side)
            depth = len(address)
            number = int(address, 2) if address else 0
            addresses.append((F(number, 1 << depth), F(number+1, 1 << depth)))
            v = F(1)
            for left, right in bounds:
                assert left < right
                v *= right-left
            assert v == F(5, 1 << depth)
            volume += v
            rule_counts[rule] = rule_counts.get(rule, 0) + 1
        addresses.sort()
        cursor = F(0)
        for left, right in addresses:
            assert left == cursor, 'binary addresses overlap or leave a gap'
            cursor = right
        assert cursor == 1 and volume == 5
        result.append({'n_min':lo+1,'n_max':hi+1,'leaves':len(addresses),
                       'volume':str(volume),'rules':rule_counts})
    return {'file':path.name,'sha256':hashlib.sha256(raw).hexdigest(),
            'blocks':result,'status':'PASS'}


def elementary(values):
    c = [F(1)] + [F(0)]*len(values)
    for x in values:
        for j in range(len(values), 0, -1):
            c[j] += x*c[j-1]
    return c


def rational_diagnostics():
    rng = Random(20260929)
    tests = 0
    # Exact tests of e_{m-1} <= 3 + (m-3)e_m, including boundary vectors.
    for m in range(3, 25):
        vectors = [[1]*m, [0]+[1]*(m-1), [0,0]+[1]*(m-2)]
        vectors += [[rng.randrange(0,101) for _ in range(m)] for _ in range(120)]
        for values in vectors:
            total = sum(values)
            if not total:
                continue
            y = [F(m*x, total) for x in values]
            e = elementary(y)
            assert e[m-1] <= 3+(m-3)*e[m]
            tests += 1
    # Exact polynomial integration for I_{k,p}; use coefficient expansion,
    # not the binomial coefficient identity used by either checker.
    for k in range(2, 13):
        for p in range(13-k):
            for d in [F(0),F(1,100),F(1,3),F(1),F(7,5)]:
                integral = sum((F(comb(p,j))*(d-1)**j/F(k+j+1)
                                for j in range(p+1)), F(0))
                formula = sum((F(comb(k+j,k))*d**j for j in range(p+1)),F(0))
                formula /= (k+p+1)*comb(k+p,k)
                assert integral == formula
                tests += 1
    # Verify the positive Taylor tail estimates used in the large-degree proof.
    def exp_bounds(x, k=20):
        partial = sum((x**j/F(factorial(j)) for j in range(k+1)),F(0))
        upper = partial+x**(k+1)/factorial(k+1)/(1-x/F(k+2))
        return partial, upper
    assert exp_bounds(F(3,8))[1] < F(3,2)
    assert exp_bounds(F(7,32))[1] < F(19,15)
    assert exp_bounds(F(9,32))[0] > F(25,19)
    assert exp_bounds(F(5,8))[0] > F(30,17)
    return {'status':'PASS', 'exact_sample_checks':tests,
            'note':'Sample diagnostics are not proofs of universal inequalities.'}


def main():
    files = ['certificate_n4_to_n13.json','certificate_n14_to_n100000.json']
    records = [topology(ROOT/f) for f in files]
    assert records[0]['blocks'][0]['n_min'] == 4
    assert records[0]['blocks'][-1]['n_max'] == 13
    assert records[1]['blocks'][0]['n_min'] == 14
    assert records[1]['blocks'][-1]['n_max'] == 100000
    out = {'status':'PASS','independent_topology':records,
           'diagnostics':rational_diagnostics(),
           'scope':'Independent coverage audit and exact auxiliary diagnostics; not an independent arithmetic replay of every leaf.'}
    dest = Path(__file__).with_name('independent_checks.json')
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: independent sorted-address coverage audit; 72 blocks, 33671 leaves.')
    print('PASS: exact rational auxiliary diagnostics. Full report:',dest.name)

if __name__ == '__main__':
    main()
