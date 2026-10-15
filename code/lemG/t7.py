import numpy as np, mpmath as mp
from resid import exact, primes
P=np.array(primes(20000),float); lp=np.log(P)
def Aval(w,s1,s2):
    a=np.exp(-(1+w+s1)*lp); b=np.exp(-(1+w+s2)*lp); c=np.exp((1+w)*lp)
    return np.exp(np.sum(np.log((1-a-b+a*b*c)*(1-a*b*c)/((1-a)*(1-b)))))
def zeta(s): return complex(mp.zeta(s))
def zl(s): return complex(mp.zeta(s,derivative=1)/mp.zeta(s))
def Ix(w,zi,zj,c=0.3):
    ts=np.concatenate([np.linspace(-400,-20,4000,endpoint=False),np.linspace(-20,20,8001),np.linspace(20,400,4000)[1:]])
    vals=[]
    for t in ts:
        s2=c+1j*t
        vals.append(zi**(-w-s2)*zj**s2*zeta(1+w)*Aval(w,-w-s2,s2)/(zeta(1-s2)*zeta(1+w+s2)*(w+s2)**2*s2**2))
    return np.trapezoid(np.array(vals),ts).real/(2*np.pi)
def main(w,zi,zj):
    h=1e-4
    A12=(Aval(w,h,h)-Aval(w,h,-h)-Aval(w,-h,h)+Aval(w,-h,-h))/(4*h*h)
    dzl=(zl(1+w+h)-zl(1+w-h))/(2*h)
    return (np.log(zi)*np.log(zj)+dzl+A12-zj**(-w)/w**2+Ix(w,zi,zj)).real
for zi,zj in ((60,60),(60,200),(300,300)):
    for w in (0.05,0.01):
        ex=float(exact(mp.mpf(w),zi,zj)); m=main(w,zi,zj)
        print(zi,zj,w,"exact",ex,"main",m,"diff",ex-m,flush=True)
