# Lemma 5 (non-negativity (21)) for general eta0, sigma0: B(y) = sum_{m=1}^{M+1} b_m B_m(y) >= 0 for all real y.
import numpy as np
from fractions import Fraction as Fr
km = [Fr(1), Fr(-851,859), Fr(780,859), Fr(-525,859), Fr(171,859), Fr(28,859), Fr(-29,859)]
b = [Fr(0)]*(len(km)+1)
for m,k in enumerate(km): b[m] += k; b[m+1] += k
def Bmin(eps0, c0, ys=np.concatenate([np.linspace(0,20,200001), np.geomspace(20,1e6,20000)])):
    tot = np.zeros_like(ys)
    for m in range(1, len(b)):
        bm = float(b[m])
        if bm <= 0: tot += bm*(m + eps0/m)/(c0*c0*m*m + ys*ys)
        else:       tot += bm*(m - eps0/m)/(m*m + ys*ys)
    i = np.argmin(tot*(1+ys*ys)); return tot[i]*(1+ys[i]**2), ys[i]
if __name__ == "__main__":
    print('b =', [str(x) for x in b])
    print('BTY  (eps0=1/2000, c0=151/153):', Bmin(1/2000, 151/153))
    for e in [1/2000, 0.001, 0.002, 0.004, 0.008]:
        for c0 in [0.9864, 0.98]:
            print(e, c0, Bmin(e, c0))
