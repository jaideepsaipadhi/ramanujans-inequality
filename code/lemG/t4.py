import numpy as np, math
def mu_arr(N):
    mu=np.ones(N+1,int); isp=np.ones(N+1,bool)
    for p in range(2,N+1):
        if isp[p]: isp[2*p::p]=False; mu[p::p]*=-1; mu[p*p::p*p]=0
    return mu
for z in (60,200,1000,3000):
    mu=mu_arr(z); ds=[d for d in range(1,z+1) if mu[d]]
    Q=0.0
    for d in ds:
        for e in ds:
            Q+=mu[d]*mu[e]*math.log(z/d)*math.log(z/e)*math.gcd(d,e)/(d*e)
    print(z, Q, Q-math.log(z))
