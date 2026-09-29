"""Floating-point proposal generator; NOT a proof checker."""
import math, json, time, sys
import numpy as np

def mesh(m, refinements=4):
    L=max(4,math.ceil(math.log2(m))+2)
    tt=[0.0]
    for i in range(L,0,-1):
        lo=2.0**(-i); hi=2.0**(-i+1)
        if i==L:
            tt += [lo*j/refinements for j in range(1,refinements+1)]
        tt += [lo+(hi-lo)*j/refinements for j in range(1,refinements+1)]
    return np.array(tt)

MESH={}
def getmesh(m, ref=4):
    key=(m,ref)
    if key not in MESH:
        tt=mesh(m,ref); l=tt[:-1];u=tt[1:];d=u-l
        w10=d*(2*l+u)/6;w11=d*(l+2*u)/6
        w20=d*(3*l*l+2*l*u+u*u)/12
        w21=d*(l*l+2*l*u+3*u*u)/12
        MESH[key]=tt,l,u,w10,w11,w20,w21
    return MESH[key]

def evaluate(block, box, ref=4):
    ml,mu=block; nlo=ml+1; nhi=mu+1
    al,au,rl,ru,sl,su=box
    v=su-rl*rl;VU=1-1/mu
    if al*al>mu: return True,'degree_bound',0.
    if v<0: return True,'variance_nonnegative',0.
    if (au*au+VU)*su<1: return True,'inverse_moment',0.
    if max(0, al-ru,rl-au)**2>max(0,au*au-1/mu)*(1-sl):return True,'moment_A',0.
    arlo=al*rl; arhi=au*ru
    if max(0,1-arhi,arlo-1)**2>VU*v:return True,'moment_B',0.
    base=(1+1/ml)*(1+au*au)*su
    if base<=0 or mu*math.log(base)<2*math.log(ml+1):return True,'product',0.
    H=max(rl*rl+(1-su)/mu,1-(au*au+VU)*v)
    A=al*al*sl
    if H>1:return True,'beta_negative',0.
    if A>0:
        tv=(H+A)/(2*A)
        if 0<tv<1 and 1-(H+A)**2/(4*A)<0:return True,'beta_negative',0.
    kap2=max(0,min(1-H,VU*su));kap=math.sqrt(kap2)
    if al<=1<=au:phi=1.
    else:
        ap=au if au<1 else al
        bphi=1+(1-ap*ap)/mu
        phi=ap*bphi**(mu/2) if bphi>=0 else 1.
    core=phi*ru*su**(ml/2)
    tt,tl,tu,w10,w11,w20,w21=getmesh(ml,ref)
    b=np.maximum(0,1-H*tt-A*tt*(1-tt));ell=1-(1-kap)*tt
    pb=b**((ml-1)/2);pe=ell**(ml-1)
    ib=w10*pb[:-1]+w11*pb[1:];ie=w10*pe[:-1]+w11*pe[1:]
    rr=np.minimum(min(1,5*su/3)*ib,min(1,3*su)*ie).sum()
    direct=kap**ml+mu/nhi*core+mu*au*au*rr
    if direct<.99999999:return True,'direct',direct
    retained=math.inf
    if ml<=13:
        # Retain a certified lower bound for the factor that was omitted.
        # Freeze this lower bound on each mesh interval so endpoint chords
        # continue to bound convex functions from above.
        delta=np.maximum(0,1-au*tu)
        rb0=np.maximum(0,(ml*b[:-1]-delta**2)/(ml-1))
        rb1=np.maximum(0,(ml*b[1:]-delta**2)/(ml-1))
        re0=np.maximum(0,(ml*ell[:-1]-delta)/(ml-1))
        re1=np.maximum(0,(ml*ell[1:]-delta)/(ml-1))
        def degree_power(z,k0,k1):return np.maximum(z**k0,z**k1)
        rib=w10*degree_power(rb0,(ml-1)/2,(mu-1)/2)+w11*degree_power(rb1,(ml-1)/2,(mu-1)/2)
        rie=w10*degree_power(re0,ml-1,mu-1)+w11*degree_power(re1,ml-1,mu-1)
        rret=np.minimum.reduce([min(1,5*su/3)*ib,min(1,3*su)*ie,su*rib,su*rie]).sum()
        retained=kap**ml+mu/nhi*core+mu*au*au*rret
        if retained<.99999999:return True,'direct_retained',retained
    den=np.maximum(0,np.maximum(1-arhi*tu,arlo*tl-1))
    ib2=w20*pb[:-1]+w21*pb[1:];ie2=w20*pe[:-1]+w21*pe[1:]
    pb2=b**((ml-2)/2);pe2=ell**(ml-2)
    jb2=w20*pb2[:-1]+w21*pb2[1:];je2=w20*pe2[:-1]+w21*pe2[1:]
    with np.errstate(divide='ignore',invalid='ignore'):
        q1=np.where(den>0,(5*mu/6)*ib2/den,np.inf)
        q2=np.where(den>0,(3*mu/2)*ie2/den,np.inf)
    q3=3*mu*(mu-1)/2*jb2;q4=9*mu*(mu-1)/2*je2
    Ecoef=nhi*au**3*ru*np.minimum.reduce([q1,q2,q3,q4]).sum()
    ratio=(1+ru)*(VU**(nlo/2)*(1-rl*rl)**((ml-1)/2)+Ecoef)
    if ratio<.99999999:return True,'center_ratio',ratio
    absolute=core+(VU*v)**(nlo/2)+Ecoef*v
    if absolute<.99999999:return True,'center_absolute',absolute
    return False,'unresolved',min(direct,retained,ratio,absolute)

def build(block, maxnodes=2000000, ref=4, progress=10000):
    root=[0.,5.,0.,1.,0.,1.];stack=[(root,0,'')];leaves=[];nodes=0;start=time.time();worst=0.;maxdepth=0
    while stack:
        box,dep,path=stack.pop();nodes+=1;maxdepth=max(maxdepth,dep)
        ok,test,bound=evaluate(block,box,ref)
        if ok:
            leaves.append((path,test))
            worst=max(worst,bound)
        else:
            widths=[box[2*j+1]-box[2*j] for j in range(3)]
            ax=int(np.argmax(widths))
            mid=(box[2*ax]+box[2*ax+1])/2
            if dep>=60 or nodes>maxnodes:
                return None,{'nodes':nodes,'failed':box,'bound':bound,'depth':dep,'elapsed':time.time()-start}
            lower=box.copy();upper=box.copy();lower[2*ax+1]=mid;upper[2*ax]=mid
            stack.append((upper,dep+1,path+str(ax)+'1'))
            stack.append((lower,dep+1,path+str(ax)+'0'))
        if progress and nodes%progress==0:print('progress',block,nodes,len(stack),dep,time.time()-start,flush=True)
    return {'block':block,'refinements':ref,'leaves':leaves},{'nodes':nodes,'leaves':len(leaves),'depth':maxdepth,'elapsed':time.time()-start,'worst':worst}

if __name__=='__main__':
    import argparse
    from pathlib import Path
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('minimum_m',type=int)
    ap.add_argument('maximum_m',type=int,nargs='?')
    ap.add_argument('--refinements',type=int,default=16)
    ap.add_argument('--output',type=Path,default=Path('proposal.json'))
    args=ap.parse_args()
    upper=args.minimum_m if args.maximum_m is None else args.maximum_m
    cert,stats=build((args.minimum_m,upper),ref=args.refinements)
    print(stats,flush=True)
    if cert is None:
        raise SystemExit('UNRESOLVED: no proof certificate was produced.')
    args.output.write_text(json.dumps(cert,separators=(',',':')),encoding='utf-8')
    print('PROPOSAL WRITTEN:',args.output,flush=True)
