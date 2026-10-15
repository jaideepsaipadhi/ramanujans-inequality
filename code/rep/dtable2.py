# d*(lk, k): two zeros with lambda_i in [lam*(lk), Lam], Lam = k e-4 * lk, heights in (e^{lk-1}, e^{lk}], |gamma0-gamma1| <= d
# are impossible.  Grid: n x n in (lam0, lam1) INCLUDING endpoints, d bisection then verified on a d-grid.
import json, math, sys
import numpy as np
from rep2 import *
lamt = {float(k): v for k, v in json.load(open('lamstar2.json')).items()}
def lam_at(l):   # lambda* valid at log-height l (grid point >= l)
    ups = [k for k in lamt if k >= l - 1e-12]; return lamt[min(ups)] if ups else 1/4.896
def ok_pair(L, sig, eta0, l0, l1, d):
    for p, q in ((l0, l1), (l1, l0)):
        if rhs1(L, p, sig) + rhs2(L, p, q, d, sig) - lhs(L, p, sig, eta0) > 0: return True
    return False
def dstar(lk, Lam, n=9, nt=5):
    """zeros in the cell: height log t in (lk-1, lk], eta = 1 - beta in [lam*(lk)/lk, Lam/lk].  The test inequality is
    evaluated at the actual height (main term (a kappa/2) log t and sigma = 1 - A/log(Kt+T0) both depend on t): on a grid
    of nt heights across the cell; primary zero at height Lt."""
    de = Lam/lk; emin = lam_at(lk)/lk
    ets = np.linspace(emin, de, n)
    tests = []
    for Lt in np.linspace(max(lk - 1, 28.7), lk, nt):
        A = lam_at(math.log(K*math.exp(Lt) + T0)); sig = sigma(Lt, A); tests.append((Lt, sig))
    def ok(d):
        for Lt, sig in tests:
            for e0 in ets:
                for e1 in ets:
                    if not ok_pair(Lt, sig, de, e0*Lt, e1*Lt, d): return False
        return True
    if not ok(0.0025): return 0.0
    lo, hi = 0.0025, 1.0
    if ok(hi): return hi
    for _ in range(10):
        mid = (lo+hi)/2
        if ok(mid): lo = mid
        else: hi = mid
    for dd in np.linspace(0.0025, lo, 25):
        if not ok(dd): return float(dd) - 0.002 if dd > 0.0045 else 0.0
    return lo
if __name__ == "__main__":
    lks = [int(x) for x in sys.argv[1].split(',')]
    out = {}
    for lk in lks:
        lmin = lam_at(lk)
        for k in range(1, 400):
            Lam = k*1e-4*lk
            if Lam <= lmin: continue
            if Lam > 0.0125*(lk - 1): break
            out['%d,%d' % (lk, k)] = d = dstar(lk, Lam)
            print(lk, k, round(Lam, 4), round(d, 4), flush=True)
    json.dump(out, open('dstar2_%s.json' % sys.argv[1].replace(',', '_'), 'w'))
