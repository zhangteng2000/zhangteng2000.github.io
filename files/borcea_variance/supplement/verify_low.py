#!/usr/bin/env python3
"""Exact, independent checker of the n=4,...,13 exclusion certificate.

Acceptance uses integers and Fraction only. No floating-point operations,
logarithms, root solvers, numerical integration, or imported proposed bounds
are used in any acceptance decision. Mathematical lemmas are documented in
the accompanying borcea_variance.pdf; this is not a proof-assistant formalization.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import comb
from functools import lru_cache
import argparse, hashlib, json, time
import verify_mid as ar

S=ar.S
ZERO=F(0); ONE=F(1)


def exp_lower(x:F,terms:int=32)->F:
    """Exact positive Taylor lower bound for exp(x), x >= 0."""
    if x<0:raise ValueError('negative exponential argument')
    t=ONE;y=ONE
    for k in range(1,terms+1):
        t=t*x/k;y+=t
    return y


def floor_scaled(x:F)->int:
    return x.numerator*S//x.denominator


def sqrt_rat_up(x:F)->int:
    return ar.sqrt_up(ar.rat_up(x))


def power_sum(k:int,p:int,d:int)->int:
    """Upper bound on sum_{j=0}^p binom(k+j,k) d^j at scale S."""
    if k<0 or p<0 or d<0:raise ValueError('invalid power sum')
    ans=0;power=S
    for j in range(p+1):
        ans+=comb(k+j,k)*power
        if j<p:power=ar.mul_up(power,d)
    return ans


def integral_up(k:int,p:int,d:int)->int:
    return ar.rat_mul_up(F(1,(k+p+1)*comb(k+p,k)),power_sum(k,p,d))


def integral_exact(k:int,p:int,d:F)->F:
    return sum((F(comb(k+j,k))*d**j for j in range(p+1)),ZERO)/((k+p+1)*comb(k+p,k))


def initial(m:int,box:tuple[F,...])->dict:
    al,au,rl,ru,sl,su=box
    v=su-rl*rl;V=F(m-1,m);df=1-sl;tm=1-m*df;D2=None
    if tm>0:
        D2=df*df/tm
        V=min(V,1+df/tm-al*al)
    return dict(al=al,au=au,rl=rl,ru=ru,sl=sl,su=su,v=v,V0=V,D2=D2,
                adist=max(ZERO,al-ru,rl-au),
                ardist=max(ZERO,1-au*ru,al*rl-1),
                base=F(m+1,m)*(1+au*au)*su)


def enrich(m:int,c:dict)->None:
    al,au,rl,ru,sl,su=(c[k] for k in ('al','au','rl','ru','sl','su'))
    V=c['V0'];dist=max(ZERO,al*sl-ru,rl-au*su)
    varx=dist*dist/V if V>0 else ZERO
    loss=m*varx/(2*su) if su>0 else ZERO
    B=F(m+1,m)*(1+au*au)
    gap=max(ZERO,al-(1+au*au)*ru/2,(1+al*al)*rl/2-au)
    C=F(m+1,m*m)*max(ZERO,m-al*al)
    gamma=F(m,3)*F((m+1)**2,m*m)*gap*gap/C if C>0 else ZERO
    L=F((m+1)**2)*(1+gamma)/B**m
    c.update(loss=loss,gamma=gamma,L=L,V=V,eta=ZERO,cap=None)


def pack(m:int,c:dict)->None:
    """Refine the radial range and its two convex-moment maxima."""
    L,su=c['L'],c['su']
    if not 0<L<=su**m:return
    lo=ZERO;hi=su
    for _ in range(30):
        md=(lo+hi)/2
        if md*((m*su-md)/(m-1))**(m-1)<L:lo=md
        else:hi=md
    eta=lo
    if not 0<eta<1:return
    if not eta*((m*su-eta)/(m-1))**(m-1)<L:
        raise ValueError('unproved minimum radial square')
    k=0;power=ONE
    while k<m and L<=power*eta:
        k+=1;power*=eta
    z=L/power
    if not 0<=k<m or not eta<=z<=1:raise ValueError('invalid packing')
    cap=k*sqrt_rat_up((1-eta)**2/eta)+sqrt_rat_up((1-z)**2/z)
    D2=F(cap,m*S)**2
    if c['D2'] is not None:D2=min(D2,c['D2'])
    inv=(F(k)/eta+1/z+m-k-1)/m
    c.update(eta=eta,cap=cap,D2=D2,V=min(c['V0'],inv-c['al']**2),packing_k=k,packing_z=z)


@lru_cache(maxsize=None)
def constants(m:int):
    return (ar.half_pow_up(F(m,m-1),m-1),ar.half_pow_up(F(m,m-1),2*(m-1)),
            ar.half_pow_up(F(m,m-2),m-2),ar.half_pow_up(F(m,m-2),2*(m-2)))


def check_leaf(m:int,box:tuple[F,...],test:str,ref:int)->int|None:
    if not 3<=m<=12:raise ValueError('unsupported degree')
    al,au,rl,ru,sl,su=box
    if test=='near_uniform':
        if rl<F(99999,100000):raise ValueError('near-uniform range not met')
        return None
    if test=='degree_bound':
        if not al*al>m:raise ValueError('degree bound failed')
        return None
    if test=='variance_nonnegative':
        if not su-rl*rl<0:raise ValueError('variance exclusion failed')
        return None
    c=initial(m,box);V0=c['V0'];v=c['v'];ad=c['adist'];ard=c['ardist'];base=c['base']
    primitive=None
    if test=='inverse_moment':primitive=(au*au+V0)*su<1
    if test=='moment_A':primitive=ad*ad>max(ZERO,au*au+V0-1)*(1-sl)
    if test=='moment_B':primitive=ard*ard>V0*v
    if test=='product':primitive=base**m<(m+1)**2
    if test=='finite_radial':primitive=ad*ad*(1-m*(1-sl))>(1-sl)**2
    if primitive is not None:
        if not primitive:raise ValueError('failed primitive '+test)
        return None
    if test=='critical_variance_negative' and V0<0:return None
    enrich(m,c)
    loss,gamma,L=(c[k] for k in ('loss','gamma','L'))
    if test=='product_stable':
        if not base**m<(m+1)**2*exp_lower(loss):raise ValueError('stable product failed')
        return None
    if test=='product_root_stable':
        if not base**m<(m+1)**2*(1+gamma)*exp_lower(loss):raise ValueError('root-stable product failed')
        return None
    if test=='product_radial':
        if not (L>1 or (0<L<=1 and m*m*ad*ad*L>(1-L)**2)):
            raise ValueError('radial product exclusion failed')
        return None
    pack(m,c)
    V=c['V'];eta=c['eta'];cap=c['cap']
    if test=='product_logcap':
        if cap is None or not F(m)*ad>F(cap,S):raise ValueError('packing exclusion failed')
        return None
    if test=='critical_variance_negative':
        if not V<0:raise ValueError('negative critical variance exclusion failed')
        return None
    if test=='moment_A_packed':
        if not ad*ad>max(ZERO,au*au+V-1)*(1-sl):raise ValueError('packed moment A failed')
        return None
    if test=='moment_B_packed':
        if not ard*ard>V*v:raise ValueError('packed moment B failed')
        return None
    if not (V>=0 and v>=0):raise ValueError('invalid main-bound moment domain')
    H=max(rl*rl+(1-V)*(1-su),1-(au*au+V)*v)
    if c['D2'] is not None:H=max(H,al*al*(1-su)+rl*rl-c['D2'])
    A=al*al*sl
    neg=H>1
    if A>0:
        tv=(H+A)/(2*A)
        if 0<tv<1 and 1-(H+A)**2/(4*A)<0:neg=True
    if test=='beta_negative':
        if not neg:raise ValueError('negative beta exclusion failed')
        return None
    if neg or not 0<=H<=1:raise ValueError('invalid convex majorant')
    if al<=1<=au:phi=S
    else:
        ap=au if au<1 else al
        bp=1+(1-ap*ap)/m
        if bp<0:raise ValueError('invalid product factor')
        phi=ar.rat_mul_up(ap,ar.half_pow_up(bp,m))
    core=ar.rat_mul_up(ru/exp_lower(loss/2),ar.mul_up(phi,ar.half_pow_up(su,m)))
    vl=max(ZERO,sl-ru*ru)
    d2=max(ZERO,min(V*v,1-H-al*al*vl));d=sqrt_rat_up(d2)
    hL=al*rl;hU=au*ru
    if test in ('center_coeff','center_positive'):
        err2=0;err3=0
        for k in range(2,m+1):
            term=ar.mul_up(ar.half_pow_up(au*au*v*F(k-1,m-1),k),power_sum(k,m-k,d))
            err2+=term
            if k>=3:err3+=term
        end=ar.half_pow_up(d2,m+1)
        if test=='center_coeff':bound=core+end+ar.rat_mul_up(hU,err2)
        else:
            if not (hU<=1 and rl>0 and au>0):raise ValueError('positive-quadratic domain not met')
            tmin=max(ZERO,eta,1-m*(1-sl))
            minre=max(F(-1),m*rl-(m-1),(m*m*rl*rl+tmin-(m-1)**2)/(2*m*rl))
            realvar=min((m-1)*(1-rl)**2,(1-rl)*(ru-minre))
            L2=F(m,2)*(vl-2*realvar)
            if not L2>0:raise ValueError('quadratic sign not certified')
            xL=(H+al*al*sl)/(2*au)
            eps=sqrt_rat_up(max(ZERO,2*(1-xL/ru)))
            K2L=(m+1)*hL*integral_exact(2,m-2,1-hU)
            idiff=3*integral_up(2,m-2,d)+ar.rat_mul_up((m-2)*hU,integral_up(3,m-3,d))
            kd=ar.rat_mul_up((m+1)*hU,ar.mul_up(eps,idiff))
            bound=core+end+ar.rat_mul_up(au*au*F(m,2)*v,kd)+ar.rat_mul_up(hU,err3)-floor_scaled(al*al*K2L*L2)
        if bound>=S:raise ValueError(f'main test failed {test}: {bound}/{S}')
        return bound
    if test not in {'direct','direct_retained','center_ratio','center_absolute'}:
        raise ValueError('unknown exclusion test '+test)
    kap2=min(1-H,V*su)
    if not 0<=kap2<=1:raise ValueError('invalid endpoint majorant')
    kap=sqrt_rat_up(kap2)
    csq,cabs,dsq,dabs=constants(m)
    cb=min(S,ar.rat_mul_up(su,csq));ce=min(S,ar.rat_mul_up(su,cabs))
    ts,segments=ar.mesh(m,ref)
    b=[1-H*t-A*t*(1-t) for t in ts]
    if min(b)<0:raise ValueError('negative majorant')
    ell=[ar.rat_mul_up(1-t,S)+ar.rat_mul_up(t,kap) for t in ts]
    pb=[ar.half_pow_up(x,m-1) for x in b]
    pe=[ar.pow_up_int(x,m-1) for x in ell]
    if test in ('direct','direct_retained'):
        rsum=0
        for i,(_,tu,w10,w11,_,_) in enumerate(segments):
            ib=ar.rat_mul_up(w10,pb[i])+ar.rat_mul_up(w11,pb[i+1])
            ie=ar.rat_mul_up(w10,pe[i])+ar.rat_mul_up(w11,pe[i+1])
            bounds=[ar.mul_up(cb,ib),ar.mul_up(ce,ie)]
            if test=='direct_retained':
                delta=max(ZERO,1-au*tu)
                rb=[max(ZERO,(m*b[j]-delta*delta)/(m-1)) for j in (i,i+1)]
                re=[max(ZERO,(m*F(ell[j],S)-delta)/(m-1)) for j in (i,i+1)]
                rib=ar.rat_mul_up(w10,ar.half_pow_up(rb[0],m-1))+ar.rat_mul_up(w11,ar.half_pow_up(rb[1],m-1))
                rie=ar.rat_mul_up(w10,ar.half_pow_up(re[0],2*(m-1)))+ar.rat_mul_up(w11,ar.half_pow_up(re[1],2*(m-1)))
                bounds.extend([ar.rat_mul_up(su,rib),ar.rat_mul_up(su,rie)])
            rsum+=min(bounds)
        bound=ar.half_pow_up(kap2,m)+ar.rat_mul_up(F(m,m+1),core)+ar.rat_mul_up(m*au*au,rsum)
    else:
        # For m=3 the exponent (m-2)/2 is below one. Do not silently
        # apply convex chord quadrature to that branch.
        if m<4:raise ValueError('centered integral test requires m >= 4')
        pb2=[ar.half_pow_up(x,m-2) for x in b]
        pe2=[ar.pow_up_int(x,m-2) for x in ell]
        csum=0
        for i,(tl,tu,_,_,w20,w21) in enumerate(segments):
            ib=ar.rat_mul_up(w20,pb[i])+ar.rat_mul_up(w21,pb[i+1])
            ie=ar.rat_mul_up(w20,pe[i])+ar.rat_mul_up(w21,pe[i+1])
            jb=ar.rat_mul_up(w20,pb2[i])+ar.rat_mul_up(w21,pb2[i+1])
            je=ar.rat_mul_up(w20,pe2[i])+ar.rat_mul_up(w21,pe2[i+1])
            bounds=[ar.rat_mul_up(F(m*(m-1),2),ar.mul_up(dsq,jb)),ar.rat_mul_up(F(m*(m-1),2),ar.mul_up(dabs,je))]
            den=max(ZERO,1-au*ru*tu,al*rl*tl-1)
            if den>0:
                bounds.extend([ar.rat_mul_up(F(m,2)/den,ar.mul_up(csq,ib)),ar.rat_mul_up(F(m,2)/den,ar.mul_up(cabs,ie))])
            csum+=min(bounds)
        coef=ar.rat_mul_up((m+1)*au**3*ru,csum)
        if test=='center_ratio':
            end=ar.mul_up(ar.half_pow_up(V,m+1),ar.half_pow_up(1-rl*rl,m-1))
            bound=ar.rat_mul_up(1+ru,end+coef)
        else:bound=core+ar.half_pow_up(V*v,m+1)+ar.rat_mul_up(v,coef)
    if bound>=S:raise ValueError(f'main test failed {test}: {bound}/{S}')
    return bound


def verify_file(path:str)->dict:
    raw=open(path,'rb').read();obj=json.loads(raw)
    blocks=obj['blocks'] if 'blocks' in obj else [obj]
    if not blocks:raise ValueError('empty certificate')
    stats={'precision_bits':ar.BITS,'blocks':0,'leaves':0,'tests':{},'maximum_depth':0,'maximum_main_bound_integer':0,'degree_details':[]}
    start=time.monotonic();first=None;last=None
    for block in blocks:
        ml,mu=block['block'];ref=block['refinements']
        if type(ml) is not int or ml!=mu or not 3<=ml<=12:raise ValueError('invalid singleton degree')
        if last is not None and ml!=last+1:raise ValueError('degree coverage gap')
        if first is None:first=ml
        last=ml;trie={}
        for path_,test in block['leaves']:ar.insert_path(trie,path_,test)
        stack=[(trie,(ZERO,F(5),ZERO,ONE,ZERO,ONE),0)]
        count=0
        while stack:
            node,box,depth=stack.pop()
            if 'test' in node:
                if len(node)!=1:raise ValueError('malformed leaf')
                try:upper=check_leaf(ml,box,node['test'],ref)
                except Exception as e:raise ValueError(f'n={ml+1}, test={node["test"]}, box={box}: {e}') from e
                count+=1;stats['leaves']+=1
                stats['tests'][node['test']]=stats['tests'].get(node['test'],0)+1
                stats['maximum_depth']=max(stats['maximum_depth'],depth)
                if upper is not None:stats['maximum_main_bound_integer']=max(stats['maximum_main_bound_integer'],upper)
            else:
                if set(node)!={'axis',0,1}:raise ValueError('parameter coverage gap')
                ax=node['axis'];mid=(box[2*ax]+box[2*ax+1])/2
                lo=list(box);hi=list(box);lo[2*ax+1]=mid;hi[2*ax]=mid
                stack.append((node[1],tuple(hi),depth+1));stack.append((node[0],tuple(lo),depth+1))
        stats['blocks']+=1
        stats['degree_details'].append({'n':ml+1,'leaves':count})
        print(f'PASS n={ml+1}: {count} leaf boxes; cumulative={stats["leaves"]}',flush=True)
    stats.update(status='PASS',minimum_degree=first+1,maximum_degree=last+1,scale=str(S),maximum_main_bound_integer=str(stats['maximum_main_bound_integer']),minimum_main_margin_numerator=str(S-int(stats['maximum_main_bound_integer'])),sha256=hashlib.sha256(raw).hexdigest(),elapsed_seconds_reporting_only=time.monotonic()-start)
    print(json.dumps(stats,indent=2),flush=True)
    return stats

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('certificate');ap.add_argument('--bits',type=int,choices=(100,160,256),default=100);ap.add_argument('--log-json')
    args=ap.parse_args();ar.BITS=args.bits;ar.S=1<<args.bits;S=ar.S
    report=verify_file(args.certificate)
    if args.log_json:
        with open(args.log_json,'w') as f:json.dump(report,f,indent=2)
