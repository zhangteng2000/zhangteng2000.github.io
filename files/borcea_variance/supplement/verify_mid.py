#!/usr/bin/env python3
"""Independent exact checker for Borcea p=2 scalar exclusion certificates.

All acceptance decisions use Python integers and fractions.Fraction.  There is
no floating-point arithmetic, numerical quadrature, or root-finding here.
A PASS concerns precisely the degree blocks listed in the checked file; it
must NOT be described as a proof for missing degrees or as a Lean proof.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import isqrt
from functools import lru_cache
import argparse, hashlib, json, time

BITS=100
S=1<<BITS
ONE=F(1)

def ceildiv(p:int,q:int)->int:
    if q<=0: raise ValueError('nonpositive denominator')
    return -((-p)//q)

def rat_up(q:F)->int:
    return ceildiv(q.numerator*S,q.denominator)

def rat_mul_up(q:F,x:int)->int:
    if q<0 or x<0: raise ValueError('negative upper-bound operand')
    return ceildiv(q.numerator*x,q.denominator)

def mul_up(x:int,y:int)->int:
    if x<0 or y<0: raise ValueError('negative upper-bound operand')
    return ceildiv(x*y,S)

def sqrt_up(x:int)->int:
    if x<0: raise ValueError('negative square root')
    k=isqrt(x*S)
    return k+(k*k<x*S)

def pow_up_int(x:int,k:int, cap:int|None=None)->int|None:
    if x<0 or k<0: raise ValueError('invalid power')
    z=S
    while k:
        if k&1:
            z=mul_up(z,x)
            if cap is not None and z>cap:return None
        k//=2
        if k:
            x=mul_up(x,x)
            if cap is not None and x>cap:return None
    return z

def half_pow_up(q:F,twice_exp:int)->int:
    if q<0 or twice_exp<0:raise ValueError('invalid half power')
    x=rat_up(q)
    y=pow_up_int(x,twice_exp//2)
    if y is None:raise RuntimeError('unexpected truncated power')
    return mul_up(y,sqrt_up(x)) if twice_exp%2 else y

@lru_cache(maxsize=None)
def mesh(m:int,ref:int):
    if ref<=0 or ref>256:raise ValueError('bad mesh refinement')
    L=max(4,(m-1).bit_length()+2)
    ts=[F(0)]
    for i in range(L,0,-1):
        lo=F(1,1<<i);hi=2*lo
        if i==L:ts.extend(lo*F(j,ref) for j in range(1,ref+1))
        ts.extend(lo+(hi-lo)*F(j,ref) for j in range(1,ref+1))
    if ts[0]!=0 or ts[-1]!=1 or any(ts[i]>=ts[i+1] for i in range(len(ts)-1)):
        raise ValueError('mesh does not partition [0,1]')
    segments=[]
    for l,u in zip(ts[:-1],ts[1:]):
        d=u-l
        segments.append((l,u,d*(2*l+u)/6,d*(l+2*u)/6,
                         d*(3*l*l+2*l*u+u*u)/12,
                         d*(l*l+2*l*u+3*u*u)/12))
    return tuple(ts),tuple(segments)

def primitive_tests(ml,mu,box,test):
    al,au,rl,ru,sl,su=box
    V=F(mu-1,mu);v=su-rl*rl
    if test=='degree_bound':return al*al>mu
    if test=='variance_nonnegative':return v<0
    if test=='inverse_moment':return (au*au+V)*su<1
    if test=='moment_A':
        dist=max(F(0),al-ru,rl-au)
        return dist*dist>max(F(0),au*au-F(1,mu))*(1-sl)
    if test=='moment_B':
        dist=max(F(0),1-au*ru,al*rl-1)
        return dist*dist>V*v
    if test=='product':
        base=F(ml+1,ml)*(1+au*au)*su
        if base<=0:return True
        target=(ml+1)**2*S
        bound=pow_up_int(rat_up(base),mu,cap=target)
        return bound is not None and bound<target
    H=max(rl*rl+(1-su)/mu,1-(au*au+V)*v)
    A=al*al*sl
    if test=='beta_negative':
        if H>1:return True
        if A>0:
            tv=(H+A)/(2*A)
            return 0<tv<1 and 1-(H+A)**2/(4*A)<0
        return False
    return None

MAIN_TESTS={'direct','direct_retained','center_ratio','center_absolute'}

def check_leaf(ml:int,mu:int,box:tuple[F,...],test:str,ref:int)->int|None:
    """Return the scaled main-test bound, or None for a passed primitive test."""
    p=primitive_tests(ml,mu,box,test)
    if p is not None:
        if not p:raise ValueError(f'exact primitive test failed: {test}, {ml,mu}, {box}')
        return None
    if test not in MAIN_TESTS:raise ValueError(f'unknown exclusion test {test}')
    al,au,rl,ru,sl,su=box
    V=F(mu-1,mu);v=su-rl*rl
    if v<0:raise ValueError('invalid scalar bound: negative upper variance')
    H=max(rl*rl+(1-su)/mu,1-(au*au+V)*v)
    A=al*al*sl
    if not 0<=H<=1:raise ValueError('invalid endpoint gap')
    # This explicitly validates the nonnegative convex quadratic used by
    # the chord quadrature. No clipping of a negative majorant is allowed.
    if A>0:
        tv=(H+A)/(2*A)
        if 0<tv<1 and 1-(H+A)**2/(4*A)<0:
            raise ValueError('quadratic majorant is negative')
    kap2=min(1-H,V*su)
    if not 0<=kap2<=1:raise ValueError('invalid kappa')
    kap_up=sqrt_up(rat_up(kap2))
    if al<=1<=au:phi_up=S
    else:
        ap=au if au<1 else al
        bp=1+(1-ap*ap)/mu
        if bp<0:raise ValueError('invalid product factor')
        phi_square_up=rat_mul_up(ap*ap,half_pow_up(bp,2*mu))
        phi_up=sqrt_up(phi_square_up)
    core=rat_mul_up(ru,mul_up(phi_up,half_pow_up(su,ml)))
    ts,segments=mesh(ml,ref)
    b=[1-H*t-A*t*(1-t) for t in ts]
    if min(b)<0:raise ValueError('negative mesh majorant')
    ell_up=[rat_mul_up(1-t,S)+rat_mul_up(t,kap_up) for t in ts]
    # The separately rounded sum could exceed S by one unit at t=0 only
    # in a non-dyadic mesh. We do not truncate it: it remains an upper bound.
    pb=[half_pow_up(z,ml-1) for z in b]
    pe=[pow_up_int(z,ml-1) for z in ell_up]
    coeff_b=min(ONE,5*su/3);coeff_e=min(ONE,3*su)
    rsum=0
    for i,(_,_,w10,w11,_,_) in enumerate(segments):
        ib=rat_mul_up(w10,pb[i])+rat_mul_up(w11,pb[i+1])
        ie=rat_mul_up(w10,pe[i])+rat_mul_up(w11,pe[i+1])
        rsum+=min(rat_mul_up(coeff_b,ib),rat_mul_up(coeff_e,ie))
    if test=='direct_retained':
        retained_sum=0
        for i,(_,tu,w10,w11,_,_) in enumerate(segments):
            delta=max(F(0),1-au*tu)
            rb=[max(F(0),(ml*b[j]-delta*delta)/(ml-1)) for j in (i,i+1)]
            re=[max(F(0),(ml*F(ell_up[j],S)-delta)/(ml-1)) for j in (i,i+1)]
            rbp=[max(half_pow_up(z,ml-1),half_pow_up(z,mu-1)) for z in rb]
            rep=[max(half_pow_up(z,2*(ml-1)),half_pow_up(z,2*(mu-1))) for z in re]
            ib=rat_mul_up(w10,pb[i])+rat_mul_up(w11,pb[i+1])
            ie=rat_mul_up(w10,pe[i])+rat_mul_up(w11,pe[i+1])
            rib=rat_mul_up(w10,rbp[0])+rat_mul_up(w11,rbp[1])
            rie=rat_mul_up(w10,rep[0])+rat_mul_up(w11,rep[1])
            retained_sum+=min(rat_mul_up(coeff_b,ib),rat_mul_up(coeff_e,ie),
                              rat_mul_up(su,rib),rat_mul_up(su,rie))
        rsum=retained_sum
    if test in ('direct','direct_retained'):
        bound=half_pow_up(kap2,ml)+rat_mul_up(F(mu,mu+1),core)+rat_mul_up(mu*au*au,rsum)
    else:
        pb2=[half_pow_up(z,ml-2) for z in b]
        pe2=[pow_up_int(z,ml-2) for z in ell_up]
        csum=0
        for i,(tl,tu,_,_,w20,w21) in enumerate(segments):
            ib=rat_mul_up(w20,pb[i])+rat_mul_up(w21,pb[i+1])
            ie=rat_mul_up(w20,pe[i])+rat_mul_up(w21,pe[i+1])
            jb=rat_mul_up(w20,pb2[i])+rat_mul_up(w21,pb2[i+1])
            je=rat_mul_up(w20,pe2[i])+rat_mul_up(w21,pe2[i+1])
            bounds=[rat_mul_up(F(3*mu*(mu-1),2),jb),
                    rat_mul_up(F(9*mu*(mu-1),2),je)]
            den=max(F(0),1-au*ru*tu,al*rl*tl-1)
            if den>0:
                bounds.extend([rat_mul_up(F(5*mu,6)/den,ib),
                               rat_mul_up(F(3*mu,2)/den,ie)])
            csum+=min(bounds)
        ecoef=rat_mul_up((mu+1)*au**3*ru,csum)
        if test=='center_ratio':
            end=mul_up(half_pow_up(V,ml+1),half_pow_up(1-rl*rl,ml-1))
            bound=rat_mul_up(1+ru,end+ecoef)
        else:
            bound=core+half_pow_up(V*v,ml+1)+rat_mul_up(v,ecoef)
    if bound>=S:
        raise ValueError(f'exact main test failed: {test}, {ml,mu}, {box}, upper={bound}/{S}')
    return bound

def insert_path(trie:dict,path:str,test:str):
    if len(path)%2 or len(path)>240:raise ValueError('invalid path length')
    node=trie
    for k in range(0,len(path),2):
        axis=int(path[k]);side=int(path[k+1])
        if axis not in (0,1,2) or side not in (0,1):raise ValueError('invalid split')
        if 'test' in node:raise ValueError('leaf has descendants')
        if 'axis' in node and node['axis']!=axis:raise ValueError('inconsistent split axes')
        node['axis']=axis
        node=node.setdefault(side,{})
    if node:raise ValueError('duplicate or overlapping leaf')
    node['test']=test

def verify_block(block:dict, stats:dict):
    ml,mu=block['block'];ref=block['refinements']
    if type(ml) is not int or type(mu) is not int or not 13<=ml<=mu:
        raise ValueError('invalid degree range')
    trie={}
    for path,test in block['leaves']:insert_path(trie,path,test)
    stack=[(trie,(F(0),F(5),F(0),F(1),F(0),F(1)),0)]
    while stack:
        node,box,depth=stack.pop()
        if 'test' in node:
            if len(node)!=1:raise ValueError('malformed leaf')
            upper=check_leaf(ml,mu,box,node['test'],ref)
            stats['leaves']+=1
            stats['tests'][node['test']]=stats['tests'].get(node['test'],0)+1
            stats['maximum_depth']=max(stats['maximum_depth'],depth)
            if upper is not None:stats['maximum_main_bound_integer']=max(stats['maximum_main_bound_integer'],upper)
        else:
            if set(node)!= {'axis',0,1}:raise ValueError('coverage gap in binary partition')
            ax=node['axis'];mid=(box[2*ax]+box[2*ax+1])/2
            lower=list(box);upper=list(box);lower[2*ax+1]=mid;upper[2*ax]=mid
            stack.append((node[1],tuple(upper),depth+1))
            stack.append((node[0],tuple(lower),depth+1))
    return ml,mu

def verify_file(path:str):
    raw=open(path,'rb').read();obj=json.loads(raw)
    blocks=obj['blocks'] if 'blocks' in obj else [obj]
    if not blocks:raise ValueError('empty certificate')
    stats={'precision_bits':BITS,'blocks':0,'leaves':0,'tests':{},'maximum_depth':0,'maximum_main_bound_integer':0}
    start=time.monotonic();first=None;last=None
    for block in blocks:
        ml,mu=verify_block(block,stats)
        if last is not None and ml!=last+1:raise ValueError('gap or overlap between degree blocks')
        if first is None:first=ml
        last=mu;stats['blocks']+=1
        print(f'PASS block: n={ml+1}..{mu+1}; cumulative leaves={stats["leaves"]}',flush=True)
    stats.update({'status':'PASS','minimum_degree':first+1,'maximum_degree':last+1,
                  'scale':str(S),'maximum_main_bound_integer':str(stats['maximum_main_bound_integer']),
                  'minimum_main_margin_numerator':str(S-int(stats['maximum_main_bound_integer'])),
                  'sha256':hashlib.sha256(raw).hexdigest(),
                  'elapsed_seconds_reporting_only':time.monotonic()-start})
    print(json.dumps(stats,indent=2),flush=True)
    return stats

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('certificate');ap.add_argument('--log-json')
    ap.add_argument('--bits',type=int,choices=(100,160,256),default=100)
    args=ap.parse_args()
    BITS=args.bits
    S=1<<BITS
    result=verify_file(args.certificate)
    if args.log_json:
        with open(args.log_json,'w') as f:json.dump(result,f,indent=2)
