#!/usr/bin/env python3
"""Sanity and rejection tests. These are not a formal proof of the checker."""
from fractions import Fraction as F
from math import factorial, prod
from pathlib import Path
import contextlib, io, json, tempfile
import verify_mid as ar
import verify_low as lo

def rejects(action):
    try: action()
    except (ValueError,KeyError,TypeError): return
    raise AssertionError('Malformed or invalid input was accepted.')

# Exact inequalities for fixed-point rounding, with no floating-point oracle.
for q in [F(0),F(1,1000003),F(1,3),F(1),F(7,3)]:
    for k in range(21):
        assert F(ar.half_pow_up(q,k),ar.S)**2>=q**k
for x in [F(0),F(1,100001),F(1,3),F(5,2)]:
    u=lo.sqrt_rat_up(x)
    assert F(u,ar.S)**2>=x
    assert F(lo.floor_scaled(x),ar.S)<=x

# Dyadic mesh coverage and exact chord weights on constant functions.
for m in [3,12,13,99,99999]:
    points,segs=ar.mesh(m,32 if m<=12 else 4)
    assert points[0]==0 and points[-1]==1
    assert sum((w10+w11 for _,_,w10,w11,_,_ in segs),F(0))==F(1,2)
    assert sum((w20+w21 for _,_,_,_,w20,w21 in segs),F(0))==F(1,3)

# Deterministic exact examples for the symmetric-product lemma.
for m in range(3,21):
    for shift in range(1,6):
        raw=[F((j*j+shift)%11+1) for j in range(m)]
        ys=[m*x/sum(raw) for x in raw]
        p=prod(ys)
        assert p*sum(1/x for x in ys)<=3+(m-3)*p

# Exact elementary exponential bounds used by the analytic high-degree proof.
def exp_lower(x,order=32):
    return sum((x**k/factorial(k) for k in range(order+1)),F(0))
def exp_upper_small(x):
    assert 0<=x<3
    return 1+x+x*x/(2*(1-x/3))
assert exp_upper_small(F(1,2))<F(5,3)
assert exp_upper_small(F(3,8))<F(3,2)
assert exp_upper_small(F(7,32))<F(19,15)
assert exp_lower(F(9,32))>F(25,19)
assert exp_lower(F(5,8))>F(30,17)
assert exp_lower(F(3,2))>F(40,9)
assert exp_lower(F(38,5))>100
assert F(7,8)**40<F(1,100)
assert F(4,100000)+F(128,300)+F(20736,10**15)<F(1,2)
m=100000
assert F(78125,9)*F((m+1)*m*m,(m-1)**3)<10000
assert F(1440*40**6*(m+1),m)<6*10**12
assert F(25*m*m,4*(m-1)**2)<8
for l,u,h in [(F(1,10),F(1,2),1665),(F(1,2),F(3,4),170),(F(5,4),F(3,2),794),(F(3,2),F(2),2500)]:
    for a in [l,u]:assert 16*(1+a*a)**4/(a*a)<=h

# Reject bad coverage, impossible test names, and false asserted exclusions.
trie={};ar.insert_path(trie,'','product')
rejects(lambda:ar.insert_path(trie,'00','direct'))
rejects(lambda:ar.insert_path(trie,'','product'))
whole=(F(0),F(5),F(0),F(1),F(0),F(1))
rejects(lambda:lo.check_leaf(3,whole,'near_uniform',32))
rejects(lambda:lo.check_leaf(3,whole,'inverse_moment',32))
rejects(lambda:lo.check_leaf(3,whole,'invented_test',32))
rejects(lambda:ar.check_leaf(13,13,whole,'inverse_moment',4))
with tempfile.TemporaryDirectory() as tmp:
    bad=Path(tmp)/'bad.json'
    bad.write_text(json.dumps({'block':[3,3],'refinements':32,'leaves':[['00','product']]}))
    with contextlib.redirect_stdout(io.StringIO()):rejects(lambda:lo.verify_file(str(bad)))
print('PASS: exact arithmetic, chord weights, analytic constants, sample symmetric products, and malformed/false-certificate rejection tests')
