# Independent ball-arithmetic certification of BTY's Table 1:  G_p(1, x) = sum_{m<=15} p^{-m} P(m x) > m_p  for x in [0, pi]
# (P even, 2pi-periodic). Mean-value form f(X) in f(c) + f'(X)(X - c) with f' evaluated on the ball: quadratic convergence.
from flint import arb, ctx
import math, sys
ctx.prec = 96
cc = [4,-8,2,20,-9,-34,27,91,-27,-201,32,895,1949,2389,1896,949,239]
ak = [arb(sum(cc[j]*cc[j+k] for j in range(17-k))*(1 if k == 0 else 2))/14912370 for k in range(17)]
from bty2 import mp_tab
def f(p, x):
    return sum((arb(p)**(-m)*sum((ak[k]*(m*k*x).cos() for k in range(17)), arb(0)) for m in range(1, 16)), arb(0))
def fd(p, x):
    return sum((arb(p)**(-m)*sum((-m*k*ak[k]*(m*k*x).sin() for k in range(17)), arb(0)) for m in range(1, 16)), arb(0))
def check(p, mp):
    mp = arb(mp); stack = [(0.0, math.pi)]; n = 0
    while stack:
        a, b = stack.pop(); n += 1
        c = arb((a+b)/2); X = arb((a+b)/2, (b-a)/2)
        if f(p, c) + fd(p, X)*(X - c) - mp > 0: continue
        if b - a < 1e-12: return False, n, (a, b)
        m = (a+b)/2; stack += [(a, m), (m, b)]
    return True, n
if __name__ == "__main__":
    ok = True
    for p in sorted(mp_tab):
        r = check(p, mp_tab[p]); ok &= r[0]; print(p, mp_tab[p], r, flush=True)
    print('ALL CERTIFIED' if ok else 'FAILURE')
