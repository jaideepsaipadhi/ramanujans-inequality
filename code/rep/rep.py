# Two-zero repulsion test built on Bellotti-Trudgian-Yang (arXiv:2603.21490): their Lemma 15 inequality with a second
# zero rho1 = 1 - eta1 + i(t+d) kept (with its Stechkin partner) in the k = 1 zero sum instead of being discarded via (21).
import numpy as np, math
th = 1.1338; sec2 = 1/math.cos(th)**2; tt = math.tan(th); ct = 1/tt; U = 2*th*ct
def w(u):
    u = np.asarray(u, float)
    return np.where(u <= U, sec2*(sec2*(th*ct - u/2)*np.cos(u*tt) + 2*th*ct - u + np.sin(2*th - u*tt)/math.sin(2*th)
                                  - 2*(1 + np.sin(th - u*tt)/math.sin(th))), 0.0)
w0 = float(w(0.0))
uu = np.linspace(0, U, 20001); wu = w(uu); du = uu[1]-uu[0]
wt = np.full(uu.size, du); wt[0] = wt[-1] = du/2            # trapezoid weights
def W(z): return complex(np.sum(np.exp(-z*uu)*wu*wt))
km = [1, -851/859, 780/859, -525/859, 171/859, 28/859, -29/859]
kap = sum(km); a0 = 1.0; a1 = 865534/497079; a = 2919857/828465
def F(z, eta, sig): return sum(k*W((z + (2*sig-1)*m)/eta) for m, k in enumerate(km))
def lhs(L, eta0, sig, kfix=True, errfac=1.0):
    f0 = eta0*w0*(kap if kfix else 1.0)
    return (eta0*w0*((a*kap/2)*L - 1.568) - 0.1185*f0 + a0*F(sig-1, eta0, sig).real + 1e-10*eta0 + 1e-7
            + errfac*(187*eta0**2 + 8013*eta0**3))
def rhs1(eta0, sig): return a1*(F(sig-1+eta0, eta0, sig).real + F(sig-eta0, eta0, sig).real)
def rhs2(eta0, eta1, d, sig): return a1*(F(sig-1+eta1-1j*d, eta0, sig).real + F(sig-eta1-1j*d, eta0, sig).real)
if __name__ == "__main__":
    print("w0 = %.9f  kappa = %.6f  a = %.6f  a1 = %.6f  a*kap*w0/2 = %.5f" % (w0, kap, a, a1, a*kap*w0/2))
    # (1) reproduce the single-zero threshold at t = exp(L): smallest lambda with lhs <= rhs1 (i.e. zero allowed)
    A = 1/4.896; K, T0 = 16, 1e10
    for L in (28.73, 30, 35, 50, 76.47):
        sig = 1 - A/math.log(K*math.exp(L) + T0)
        for kfix in (False, True):
            lo, hi = 0.15, 0.30
            for _ in range(40):
                mid = (lo+hi)/2; e0 = mid/L
                if lhs(L, e0, sig, kfix) < rhs1(e0, sig): lo = mid      # contradiction -> zero impossible at this lambda
                else: hi = mid
            print("L = %6.2f  kappa-fix %-5s : single-zero threshold lambda = %.5f  (R = %.4f)" % (L, kfix, lo, 1/lo))

def lemma11_ok(etamax, u0=59.0, sig0=0.9935164):
    """BTY Lemma 11 with eta0 replaced by etamax: w(etamax u) sum kappa_m e^{-d_m u} >= kappa w0 on [0, u0]."""
    u = np.linspace(0, u0, 20001)
    s = sum(k*np.exp(-(((2*sig0-1)*m) if k < 0 else m)*u) for m, k in enumerate(km))
    return float(np.min(w(etamax*u)*s - kap*w0))
def lam_star(L, A=1/4.896, errfac=1.5, kfix=True):
    sig = 1 - A/math.log(16*math.exp(L) + 1e10)
    lo, hi = 0.15, 0.30
    for _ in range(40):
        mid = (lo+hi)/2; e0 = mid/L
        if lhs(L, e0, sig, kfix, errfac) < rhs1(e0, sig): lo = mid
        else: hi = mid
    return lo

def lam_star_A(L, A, errfac=1.5):
    sig = 1 - A/math.log(16*math.exp(L) + 1e10)
    lo, hi = 0.15, 0.32
    for _ in range(40):
        mid = (lo+hi)/2; e0 = mid/L
        if lhs(L, e0, sig, True, errfac) < rhs1(e0, sig): lo = mid
        else: hi = mid
    return lo
def iterate(Ls, n_iter=6):
    # A(l): zero-free constant valid at height e^l; start from BTY's 1/4.896 everywhere (also used beyond the grid)
    lam = {L: 1/4.896 for L in Ls}
    def Aat(l):   # conservative: lambda at the nearest grid point ABOVE l (lambda* decreasing in l); beyond grid: BTY
        ups = [L for L in Ls if L >= l]
        return min(lam[min(ups)], 10) if ups else 1/4.896
    for it in range(n_iter):
        new = {}
        for L in Ls:
            lt = math.log(16*math.exp(L) + 1e10)
            new[L] = max(1/4.896, lam_star_A(L, Aat(lt)))
        lam = new
        print("iter", it, [(L, round(lam[L], 5)) for L in Ls[::4]], flush=True)
    return lam

def dstar(L, Lam, A=None, ngrid=7):
    """largest d (on a grid) s.t. any two zeros at heights ~e^L with lambda_i in [lam*(L), Lam] and |gamma0-gamma1| <= d
    give a contradiction (for at least one choice of primary zero)."""
    import json
    lamt = {float(k): v for k, v in json.load(open('lamstar.json')).items()}
    lmin = min(v for k, v in lamt.items() if k >= L)
    if A is None:
        lt = math.log(16*math.exp(L) + 1e10); A = min(v for k, v in lamt.items() if k >= lt - 1e-9) if any(k >= lt for k in lamt) else 1/4.896
    sig = 1 - A/math.log(16*math.exp(L) + 1e10)
    lams = np.linspace(lmin, Lam, ngrid)
    def ok(d):
        for l0 in lams:
            for l1 in lams:
                good = False
                for (p, q) in ((l0, l1), (l1, l0)):
                    e0, e1 = p/L, q/L
                    m = rhs1(e0, sig) + rhs2(e0, e1, d, sig) - lhs(L, e0, sig, True, 1.5)
                    if m > 0: good = True; break
                if not good: return False
        return True
    if not ok(0.0025): return 0.0
    lo, hi = 0.0025, 2.0
    if ok(hi): return hi
    for _ in range(9):
        mid = (lo+hi)/2
        if ok(mid): lo = mid
        else: hi = mid
    return lo          # NOTE: bisection assumes the contradiction region in d is an interval [0, d*]; verified on a grid afterwards
