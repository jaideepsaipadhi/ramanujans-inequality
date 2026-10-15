from resid import *
mp.mp.dps=30
for z in (60,300,1000):
    for w in (mp.mpf('0.05'),mp.mpf('0.01')):
        r,m,r0,r1=R(w,z,z); ex=exact(w,z,z)
        print(z,float(w),"exact",float(ex),"main",float(r),"diff*w",float((ex-r)*w),flush=True)
