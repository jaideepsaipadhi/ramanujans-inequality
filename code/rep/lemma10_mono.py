# BTY Lemma 10 needs dG_p/dsigma < 0 on sigma in [sigma0, 1], x in [0, 2pi]; BTY checked sigma0 = 0.99352. We need 0.9930.
# G_p(sigma, x) = sum_{k<=16} a_k sum_{m<=15} p^{-m sigma} cos(m k x);  dG/dsigma = -log p * sum a_k sum m p^{-m sigma} cos(mkx).
from flint import arb, ctx
import math
ctx.prec = 64
cc = [4,-8,2,20,-9,-34,27,91,-27,-201,32,895,1949,2389,1896,949,239]
ak = [arb(sum(cc[j]*cc[j+k] for j in range(17-k))*(1 if k == 0 else 2))/14912370 for k in range(17)]
def dG(p, sg, x):
    lp = arb(p).log(); s = arb(0)
    for m in range(1, 16):
        pm = (-m*sg*lp).exp()
        s += m*pm*sum((ak[k]*(m*k*x).cos() for k in range(17)), arb(0))
    return -lp*s
def check(p, s0=0.9930):
    stack = [(s0, 1.0, 0.0, math.pi)]          # cos even in x and 2pi-periodic: x in [0, pi] suffices
    n = 0
    while stack:
        a, b, c, d = stack.pop(); n += 1
        v = dG(p, arb((a+b)/2, (b-a)/2), arb((c+d)/2, (d-c)/2))
        if v < 0: continue
        if d - c < 1e-6: return False, n, (a, b, c, d), v
        if (b - a)/(1 - s0) > (d - c)/math.pi: m = (a+b)/2; stack += [(a, m, c, d), (m, b, c, d)]
        else: m = (c+d)/2; stack += [(a, b, c, m), (a, b, m, d)]
    return True, n
import sys
for p in [int(x) for x in sys.argv[1].split(",")]:
    print(p, check(p), flush=True)
