# K(w) = zeta(1+w) H(w);  H(w) = (1/2pi i) int_{Re s=c} A(w;-w-s,s) / (zeta(1-s) zeta(1+w+s) (w+s)^2 s^2) ds
import numpy as np, mpmath as mp, sys
from resid import primes
P=np.array(primes(int(sys.argv[1]) if len(sys.argv)>1 else 200000),float); lp=np.log(P)
c=0.3; T=float(sys.argv[2]) if len(sys.argv)>2 else 300; h=0.02
ts=np.arange(-T,T+h/2,h); S=c+1j*ts
iz1=np.array([1/complex(mp.zeta(1-s)) for s in S])
def A(w,s1,s2):
    a=np.exp(-np.outer(1+w+s1,lp)); b=np.exp(-np.outer(1+w+s2,lp)); cc=np.exp((1+w)*lp)
    f=(1-a-b+a*b*cc)*(1-a*b*cc)/((1-a)*(1-b)); return np.exp(np.sum(np.log(f),axis=1))
def H(w):
    Av=np.concatenate([A(w,-w-S[i:i+2000],S[i:i+2000]) for i in range(0,len(S),2000)])
    iz2=np.array([1/complex(mp.zeta(1+w+s)) for s in S])
    f=Av*iz1*iz2/((w+S)**2*S**2)
    return np.trapezoid(f,ts)/(2*np.pi)
e=2e-3
H0=H(0.0); Hp=H(e); Hm=H(-e)
k1=H0; k0=(Hp-Hm)/(2*e)+float(mp.euler)*H0
print("H(0)=",H0," H'(0)=",(Hp-Hm)/(2*e)); print("k1=",k1.real," k0=",k0.real, " (imag parts", k1.imag, k0.imag,")")
