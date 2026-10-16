"""The certified scan (Proposition 3.6 of the paper).
 * Covers every x in a given range [e^L1, e^L2]: for each interval [L1,L2] the kernel (c, eps) is designed
   for the smallest tolerance tol(L2), and the zero tail is evaluated at L1 (weakest damping x^(beta-1)).
 * pi vs psi: pi-li = (theta-x)/L + int (theta-t)/(t log^2 t) dt  =>  need psi error <= tol/(1+2/L1).
 * Brun-Titchmarsh small-|t| slack: smoothing cost inflated by 2%."""
import sys, math, os, ramanujan_upper as r, logan, lemmaA, lemmaD
WD = int(os.environ.get('WEIGHTED', '0'))      # Lemma W cut: zeros with 1 - beta <= WD*1e-4 handled by the weighted bound
def edge_weighted(L1, H, eps, c, nt=80):
    """sum over zeros with gamma in [H, T1], 1 - beta <= WD e-4 of x^{beta-1} f(gamma), f = w(gamma eps, c)/gamma decreasing,
    via Stieltjes in gamma with W(t) = sum_{gamma <= t} x^{beta-1} <= exp(weighted_logA(WD, ceil log t, Lx))."""
    T1 = c/eps
    if T1 <= H: return 0.0
    us = [math.log(H)+(math.log(T1)-math.log(H))*i/nt for i in range(nt+1)]
    f = [logan.w(math.exp(u)*eps, c)/math.exp(u) for u in us]
    Lk = int(math.floor(L1/10)*10)                 # x^{-delta} at L >= Lk (weaker damping: conservative)
    def Wb(u):
        v = lemmaD.weighted_logA(WD, math.ceil(u), Lk)
        return 0.0 if v is None else math.exp(v)
    acc = sum(Wb(us[i+1])*(f[i]-f[i+1]) for i in range(nt)) + Wb(us[-1])*f[-1]
    return acc
CS=list(range(40,105,5)); FR=[0.45+0.05*i for i in range(11)]
def err_interval(L1, L2, H):
    # REFEREE FIX (Oct 10): the inequality at x needs the psi-error at x AND at y = x/e below eps_max(log x).
    # Certify every interval at tol(L2 + 1): then for y in an interval [L1', L2'] with log x = log y + 1 <= L2' + 1,
    # tol(L2' + 1) <= eps_max(log x) (eps_max decreasing), and for x itself tol(L2 + 1) <= eps_max(L2).
    tol = r.eps_max(L2 + 1)/(1+2/L1)
    def one(c, frac):
        eps = frac*tol*math.sqrt(c)/(0.7946*1.02)
        for _ in range(3):
            lam = math.log(10*math.sqrt(c)/eps); kap = (L2+lam)/(L1-lam)
            eps = frac*tol*math.sqrt(c)/(0.7946*1.02*kap)
        T1 = c/eps
        tr = 0.32/math.sinh(c)*math.exp(0.71*math.sqrt(c*eps))*math.log(3*c)*math.log(T1)
        tail = 0.0
        if T1 > H:
            # RIGOROUS: on [a,b], S(s) <= S(a) (S decreasing in s), times the exact integral of x^(s-1) (Oct 13)
            smax = 1-1/(r.Rof(math.log(T1))*math.log(T1))
            s0 = 1-WD/1e4 if WD else smax          # Lemma W handles beta >= s0
            smax_int = min(smax, s0); s_mid = smax_int-0.02
            grid = [0.6+(s_mid-0.6)*k/60 for k in range(61)]
            hf = 1/(40*L1); k=1
            while s_mid+k*hf < smax_int: grid.append(s_mid+k*hf); k+=1
            grid.append(smax_int)
            acc=0.0
            for a,b in zip(grid[:-1],grid[1:]):
                acc += logan.zsum(a,H,eps,c)*(math.exp((b-1)*L1) - math.exp((a-1)*L1))/L1   # S decreasing: S(a) * exact int of x^(s-1)
            tail = 2*(L2*acc + math.exp(-0.4*L1)*(logan.zsum(0.6,H,eps,c)+math.log(T1)**2))
            if WD and s0 < smax: tail += 2*edge_weighted(L1, H, eps, c)
        e = frac*tol + tr + tail + math.exp(-L1/2)*math.log(H)**2
        return e/tol
    # any (c, frac) gives a valid bound: coarse grid, then local refinement around the best (Oct 13, speed)
    v, c0, f0 = min((one(c, f), c, f) for c in range(40, 101, 10) for f in (0.5, 0.6, 0.7, 0.8, 0.9))
    for c in range(max(30, c0-10), c0+11, 5):
        for f in [f0 + 0.05*k for k in range(-2, 3)]:
            if 0.3 <= f <= 0.98: v = min(v, one(c, f))
    return v
if __name__ == "__main__":
    H=float(sys.argv[1]); lo=int(sys.argv[2]); hi=int(sys.argv[3]); w=int(sys.argv[4])
    worst=(0,None); fails=[]
    L=lo
    while L < hi:
        q=err_interval(L, min(L+w,hi), H)
        if q>1: fails.append(L)
        if q>worst[0]: worst=(q,L)
        L+=w
    print('H=%.3e [%d,%d) width %d: worst ratio %.3f at L=%d; failing intervals start at %s'%(H,lo,hi,w,worst[0],worst[1],fails[:3]+(['...',fails[-1]] if len(fails)>3 else [])),flush=True)
