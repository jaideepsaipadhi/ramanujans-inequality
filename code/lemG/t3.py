import mpmath as mp, numpy as np
mp.mp.dps=30
def mu_arr(N):
    mu=np.ones(N+1,int); isp=np.ones(N+1,bool)
    for p in range(2,N+1):
        if isp[p]: isp[2*p::p]=False; mu[p::p]*=-1; mu[p*p::p*p]=0
    return mu
for z in (60,300,3000,30000):
    mu=mu_arr(z)
    for w in (mp.mpf('0.1'),mp.mpf('0.01')):
        s=mp.fsum(mu[d]*mp.log(mp.mpf(z)/d)/mp.mpf(d)**(1+w) for d in range(1,z+1) if mu[d])
        ex=mp.zeta(1+w)*s; main=mp.log(z)-mp.zeta(1+w,derivative=1)/mp.zeta(1+w)
        print(z,float(w),float(ex),float(main),float(ex-main),float((ex-main)*w))
