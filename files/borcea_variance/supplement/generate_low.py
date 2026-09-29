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
    if 3<=ml<=mu<=12 and rl>=1-1/100000:return True,'near_uniform',0.
    v=su-rl*rl;VU=1-1/mu
    df=1-sl;tm=1-mu*df
    Dfinite2=math.inf
    if tm>0:
        Dfinite2=df*df/tm
        VU=min(VU,1+df/tm-al*al)
        if VU<0:return True,'critical_variance_negative',0.
    if al*al>mu: return True,'degree_bound',0.
    if v<0: return True,'variance_nonnegative',0.
    if (au*au+VU)*su<1: return True,'inverse_moment',0.
    if max(0, al-ru,rl-au)**2>max(0,au*au+VU-1)*(1-sl):return True,'moment_A',0.
    arlo=al*rl; arhi=au*ru
    if max(0,1-arhi,arlo-1)**2>VU*v:return True,'moment_B',0.
    base=(1+1/ml)*(1+au*au)*su
    if base<=0 or mu*math.log(base)<2*math.log(ml+1):return True,'product',0.
    # Covariance forces a radial variance, hence a strict AM--GM loss.
    dist=max(0,al*sl-ru,rl-au*su)
    variance_x=dist*dist/VU if VU>0 else 0.
    loss=ml*variance_x/(2*su) if su>0 else 0.
    if mu*math.log(base)-loss<2*math.log(ml+1):return True,'product_stable',0.
    adist=max(0,al-ru,rl-au)
    if adist**2*(1-mu*(1-sl))>(1-sl)**2:return True,'finite_radial',0.
    minq=0.
    rootloss=0.
    if ml==mu and ml<=20:
        B=(ml+1)/ml*(1+au*au)
        # A-B*mu/2 = (n/m)(a-(1+a^2)*mu/2)
        gap=max(0,al-(1+au*au)*ru/2,(1+al*al)*rl/2-au)
        UV=(ml+1)/ml**2*max(0,ml-al*al)
        rootloss=(ml/3)*((ml+1)/ml*gap)**2/UV if UV>0 else 0.
        if mu*math.log(base)-loss-math.log1p(rootloss)<2*math.log(ml+1):return True,'product_root_stable',0.
        Qlow=(ml+1)**2/B**ml*(1+rootloss)
        if Qlow>1 or ml**2*adist**2*Qlow>(1-Qlow)**2:return True,'product_radial',0.
        if Qlow<=su**ml and Qlow>0:
            lo=0.;hi=su
            for _ in range(30):
                md=(lo+hi)/2
                if md*((ml*su-md)/(ml-1))**(ml-1)<Qlow:lo=md
                else:hi=md
            eta=lo
            if eta>0 and eta<1:
                kk=min(ml,int(math.log(Qlow)/math.log(eta)))
                zz=Qlow/eta**kk
                ff=lambda z:1/math.sqrt(z)-math.sqrt(z)
                cap=kk*ff(eta)+ff(zz)
                if ml*adist>cap:return True,'product_logcap',0.
                Dfinite2=min(Dfinite2,(cap/ml)**2)
                minq=eta
                tU=(kk/eta+1/zz+ml-kk-1)/ml
                VU=min(VU,tU-al*al)
                if VU<0:return True,'critical_variance_negative',0.
                if max(0,1-arhi,arlo-1)**2>VU*v:return True,'moment_B_packed',0.
                if adist**2>max(0,au*au+VU-1)*(1-sl):return True,'moment_A_packed',0.
    H=max(rl*rl+(1-VU)*(1-su),1-(au*au+VU)*v,al*al*(1-su)+rl*rl-Dfinite2)
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
    core=phi*ru*su**(ml/2)*math.exp(-loss/2)
    Csq=(mu/(mu-1))**((mu-1)/2)
    Cabs=(mu/(mu-1))**(mu-1)
    Dsq=(mu/(mu-2))**((mu-2)/2)
    Dabs=(mu/(mu-2))**(mu-2)
    coeffbound=math.inf
    positivebound=math.inf
    if ml==mu and ml<=20:
        vl=max(0,sl-ru*ru)
        d2=max(0,min(VU*v,1-H-al*al*vl))
        d=math.sqrt(d2)
        err=0.
        for k in range(2,ml+1):
            J=sum(math.comb(k+j,k)*d**j for j in range(ml-k+1))
            err+=(au*au*v*(k-1)/(ml-1))**(k/2)*J
        coeffbound=core+d**(ml+1)+au*ru*err
        if coeffbound<.99999999:return True,'center_coeff',coeffbound
        if au*ru<=1 and rl>0:
            tmin=max(0,minq,1-ml*(1-sl))
            minre=max(-1,ml*rl-(ml-1),(ml*ml*rl*rl+tmin-(ml-1)**2)/(2*ml*rl))
            realvar=min((ml-1)*(1-rl)**2,(1-rl)*(ru-minre))
            L2=ml/2*(vl-2*realvar)
            if L2>0:
                xl=(H+al*al*sl)/(2*au)
                eps=math.sqrt(max(0,2*(1-xl/ru)))
                def integral(k,p,z):
                    return sum(math.comb(k+j,k)*z**j for j in range(p+1))/((k+p+1)*math.comb(k+p,k))
                hL=al*rl;hU=au*ru
                K2L=(ml+1)*hL*integral(2,ml-2,1-hU)
                Kdiff=(ml+1)*hU*eps*(3*integral(2,ml-2,d)+(ml-2)*hU*integral(3,ml-3,d))
                err3=0.
                for k in range(3,ml+1):
                    J=sum(math.comb(k+j,k)*d**j for j in range(ml-k+1))
                    err3+=(au*au*v*(k-1)/(ml-1))**(k/2)*J
                positivebound=core+d**(ml+1)-al*al*K2L*L2+au*au*ml*v/2*Kdiff+hU*err3
                if positivebound<.99999999:return True,'center_positive',positivebound
    tt,tl,tu,w10,w11,w20,w21=getmesh(ml,ref)
    b=np.maximum(0,1-H*tt-A*tt*(1-tt));ell=1-(1-kap)*tt
    pb=b**((ml-1)/2);pe=ell**(ml-1)
    ib=w10*pb[:-1]+w11*pb[1:];ie=w10*pe[:-1]+w11*pe[1:]
    rr=np.minimum(min(1,Csq*su)*ib,min(1,Cabs*su)*ie).sum()
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
        rret=np.minimum.reduce([min(1,Csq*su)*ib,min(1,Cabs*su)*ie,su*rib,su*rie]).sum()
        retained=kap**ml+mu/nhi*core+mu*au*au*rret
        if retained<.99999999:return True,'direct_retained',retained
    den=np.maximum(0,np.maximum(1-arhi*tu,arlo*tl-1))
    ib2=w20*pb[:-1]+w21*pb[1:];ie2=w20*pe[:-1]+w21*pe[1:]
    pb2=b**((ml-2)/2);pe2=ell**(ml-2)
    jb2=w20*pb2[:-1]+w21*pb2[1:];je2=w20*pe2[:-1]+w21*pe2[1:]
    with np.errstate(divide='ignore',invalid='ignore'):
        q1=np.where(den>0,(Csq*mu/2)*ib2/den,np.inf)
        q2=np.where(den>0,(Cabs*mu/2)*ie2/den,np.inf)
    q3=Dsq*mu*(mu-1)/2*jb2 if ml>=4 else np.full_like(jb2,np.inf)
    q4=Dabs*mu*(mu-1)/2*je2
    Ecoef=nhi*au**3*ru*np.minimum.reduce([q1,q2,q3,q4]).sum()
    ratio=(1+ru)*(VU**(nlo/2)*(1-rl*rl)**((ml-1)/2)+Ecoef)
    if ratio<.99999999:return True,'center_ratio',ratio
    absolute=core+(VU*v)**(nlo/2)+Ecoef*v
    if absolute<.99999999:return True,'center_absolute',absolute
    return False,'unresolved',min(direct,retained,ratio,absolute,coeffbound,positivebound)

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
