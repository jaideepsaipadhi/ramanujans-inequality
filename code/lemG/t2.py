from resid import *
mp.mp.dps=40
for zi,zj in ((60,200),(150,400),(300,1000)):
    for w in (mp.mpf('0.1'),mp.mpf('0.03'),mp.mpf('0.01')):
        r,m,r0,r1=R(w,zi,zj); ex=exact(w,zi,zj)
        print(f"zi={zi} zj={zj} w={float(w)}: exact {float(ex):.5f} main {float(r):.5f} diff {float(ex-r):+.3e}",flush=True)
