# Lemma H / Lemma Rep, v2: exact W (closed form), constants re-derived for general eta0 (bty2.py), Lemma 5 per-sigma.
import math, json, functools
from Wexact import W
import bty2
w0 = float(bty2.w0); km = bty2.km; kap = bty2.kap; a0, a1, a = bty2.a0, bty2.a1, bty2.a
K, T0 = 16, 1e10
def F(z, eta, sig): return sum(k*W((z + (2*sig-1)*m)/eta) for m, k in enumerate(km))
@functools.lru_cache(maxsize=None)
def cst(eta0r, sig0r):
    eta0, sig0 = eta0r/1e5, sig0r/1e6
    E2, E3 = bty2.errconsts(eta0, sig0); c10 = bty2.lemma10(eta0, sig0)[0]; G = bty2.gamma_const(sig0)
    # Lemma 5 per sigma: need eps0 = Cw eta0^2/(w0 (2 sig0 - 1)^2) <= 0.0859
    c = 2*sig0 - 1; eps0 = bty2.Cw(c/eta0, c/eta0)*eta0**2/(w0*c*c); assert eps0 < 0.002, eps0                    # lemma5_cert.py certifies eps0 = 1/500
    return E2, E3, c10, G
def cfor(eta0, sig): return cst(math.ceil(eta0*1e4)*10, math.floor(sig*1e5)*10)   # coarse: conservative (eta0 up, sigma0 down)
def lhs(L, lam0, sig, eta0):
    e = lam0/L; E2, E3, c10, G = cfor(eta0, sig)
    return (e*w0*((a*kap/2)*L - G) - c10*e*w0 + a0*F(sig-1, e, sig).real + 1e-10*e + 2e-7 + E2*lam0*e*e + E3*e**3)
def rhs1(L, lam0, sig):
    e = lam0/L; return a1*(F(sig-1+e, e, sig).real + F(sig-e, e, sig).real)
def rhs2(L, lam0, lam1, d, sig):
    e0, e1 = lam0/L, lam1/L
    return a1*(F(sig-1+e1-1j*d, e0, sig).real + F(sig-e1-1j*d, e0, sig).real)
def sigma(L, A): return 1 - A/math.log(K*math.exp(L) + T0)
def lam_star(L, A, eta0=0.013):
    sig = sigma(L, A); lo, hi = 0.15, 0.32
    for _ in range(40):
        mid = (lo+hi)/2
        if lhs(L, mid, sig, eta0) < rhs1(L, mid, sig): lo = mid      # contradiction: no zero at this lambda
        else: hi = mid
    return lo
def iterate(Ls, n_iter=5):
    lam = {L: 1/4.896 for L in Ls}
    def Aat(l):
        ups = [L for L in Ls if L >= l]
        return lam[min(ups)] if ups else 1/4.896
    for it in range(n_iter):
        lam = {L: max(1/4.896, lam_star(L, Aat(math.log(K*math.exp(L)+T0)))) for L in Ls}
        print("iter", it, [(L, round(lam[L], 5)) for L in Ls[::4]], flush=True)
    return lam
if __name__ == "__main__":
    import sys
    Ls = [28 + 0.25*i for i in range(61)]
    lam = iterate(Ls)
    json.dump({str(L): v for L, v in lam.items()}, open('lamstar2.json', 'w'))
