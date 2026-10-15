"""Lemma F: rigorous explicit upper bound for B = eps*log(z2/z1)*sum_n Psi(n)^2 n^{-1-eps}
(Barban-Vehov weights lam_d = mu(d) min(1, log(z2/d)/log(z2/z1))), z1 = e^{Lz}, z2 = z1^tau, eps = kappa/Lz.

Proof chain (see LEMMA_F.md):
 (1) Selberg diagonalisation:  sum Psi^2 n^{-s} = zeta(s)/Delta^2 * sum_m W(m) Yt(m)^2,
     W(m) = prod_{p|m} 1/(p^s-1) <= m^{-eps} mu^2(m)/phi(m),  Yt = (N(z2/m;m) - N(z1/m;m))/H_m,
     N(y;q) = sum_{k<=y,(k,q)=1} mu(k) k^{-s} log(y/k),  H_m = prod_{p|m} (1-p^{-s})^{-1}.
 (2) N(y;q) = sum_{d|q^inf} d^{-s} N(y/d;1);  N(u;1) = A(log u + b) + r(u) (u>=1), = 0 (u<1);
     A = 1/zeta(s), b = -zeta'/zeta(s).  Residue: main_q(y) = H A (log y + Hp + b), Hp = -sum_{p|q} log p/(p^s-1).
 (3) |r(u)| <= Rbar(u): exact on [1,U0] (N(u;1) piecewise log-linear => extrema at integer endpoints),
     2c E1(eps log u) for u >= U0 >= 96955 (|sum_{k<=t} mu(k)/k| <= c/log t, c = 0.0144, Ramare);
     for u < 1: |r(u)| = A|log u + b| <= R1 + A log(1/u),  R1 := Rbar(1) (>= A b).
 (4) averaging over d with weights d^{-s}/H (sum 1), cut at D, Rankin:
     |Yt| <= |M| + sum_i chi_i [Rbar(y_i/D_i) + R1 D_i^{-al} g(m)] + 2 A h(m),
     g(m) = prod_{p|m} (1-1/p)/(1-p^{al-1}),  h(m) = sum_{p|m} log p/(p-1).
 (5) Minkowski over the three pieces; m <= M0 summed exactly, m > M0 by log-intervals with
     sum_{a<m<=b} mu^2/phi * f(m) <= P_f log(b/a) + 2*3.95 E_f/sqrt(a)  (f = 1, g^2, h^2, via f = 1*gamma;
     |sum_{m<=x} mu^2/phi - log x - 1.3326| <= 3.95/sqrt x, x >= 1 [checked to 1e8: max 1.45]).
"""
import numpy as np, mpmath as mp
from scipy.special import exp1

C_MU = 1/69; U0 = 2*10**6; M0 = 10**6; GERR = 4.0   # Ramare-Akhilesh, Colloq. Math. 147 (2017): 4/sqrt(D), D > 1

def _sieve(N):
    mu = np.ones(N + 1, np.int8); mu[0] = 0
    isp = np.ones(N + 1, bool); isp[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if isp[p]: isp[p*p::p] = False
    primes = np.nonzero(isp)[0]
    phi = np.arange(N + 1, dtype=np.float64); h = np.zeros(N + 1); lg = np.zeros(N + 1)
    for p in primes:
        mu[p::p] *= -1; phi[p::p] *= (1 - 1.0/p); h[p::p] += np.log(p)/(p - 1)
        if p*p <= N: mu[p*p::p*p] = 0
    return mu, phi, h, primes
_MU, _PHI, _H, _PR = _sieve(max(U0, M0))
_IDX = np.arange(len(_MU), dtype=np.float64)

def _gtab(al):
    g = np.ones(M0 + 1)
    for p in _PR:
        if p > M0: break
        g[p::p] *= (1 - 1.0/p)/(1 - p**(al - 1.0))
    return g
_GC = {}

def _primesum(f_small, tail_sigma, coef):
    """sum over primes of f(p): exact for p<=M0, tail <= coef * sum_{p>M0} p^{-sigma}, via pi(t) <= 1.25506 t/log t"""
    P = float(M0); sig = tail_sigma
    tail = coef*1.25506*sig/((sig - 1)*np.log(P))*P**(1 - sig)
    return f_small + tail

def addmom(al):
    """moments for h (alpha_p = log p/(p-1)) and k (alpha_p = 1/p):  returns dict of (P, E) pairs such that
    sum_{a<m<=b} mu^2/phi * F <= P log(b/a) + E 7.9/sqrt(a) for F in {k, k^2, g k, h k}.
    additive a, b >= 0:  sum mu^2/phi a = sum_p a_p/(p-1) G(I/p);  sum mu^2/phi a b <= sum_{p != q} + sum_p (diag);
    sum mu^2/phi g a = sum_p a_p g_p/(p-1) sum_{(m',p)=1} mu^2/phi g(m')."""
    p = _PR[_PR <= M0].astype(np.float64); sp = np.sqrt(p)
    lp = np.log(p)/(p - 1); kp = 1/p; gp = (1 - 1/p)/(1 - p**(al - 1))
    def S(v, sig, c): return np.sum(v) + _primesum(0, sig, c)
    Sk = S(kp/(p - 1), 2.0, 1.01); Tk = S(kp*sp/(p - 1), 1.5, 1.01)
    Sk2 = S(kp**2/(p - 1), 3.0, 1.01); Tk2 = S(kp**2*sp/(p - 1), 2.5, 1.01)
    Sh = S(lp/(p - 1), 1.9, 10.1); Th = S(lp*sp/(p - 1), 1.4, 10.1)
    Shk = S(lp*kp/(p - 1), 2.9, 10.1); Thk = S(lp*kp*sp/(p - 1), 2.4, 10.1)
    Sgk = S(kp*gp/(p - 1), 2.0, 1.1*1.01); Tgk = S(kp*gp*sp/(p - 1), 1.5, 1.1*1.01)
    return {'k': (Sk, Tk), 'kk': (Sk**2 + Sk2, Tk**2 + Tk2), 'hk': (Sh*Sk + Shk, Th*Tk + Thk), 'gk': (Sgk, Tgk)}
def euler_h(al):
    """first moments with h(m) = sum_{p|m} l_p, l_p = log p/(p-1):
       sum_{a<m<=b} mu^2/phi h <= S1 log(b/a) + T1 7.9/sqrt a;  sum mu^2/phi g h <= Sgh Pg1 log(b/a) + Tgh Eg1 7.9/sqrt a"""
    p = _PR[_PR <= M0].astype(np.float64); lp = np.log(p)/(p - 1); gp = (1 - 1/p)/(1 - p**(al - 1))
    S1 = np.sum(lp/(p - 1)) + _primesum(0, 1.9, 10.1); T1 = np.sum(lp*np.sqrt(p)/(p - 1)) + _primesum(0, 1.4, 10.1)
    Sgh = np.sum(lp*gp/(p - 1)) + _primesum(0, 1.9, 11.2); Tgh = np.sum(lp*gp*np.sqrt(p)/(p - 1)) + _primesum(0, 1.4, 11.2)
    return S1, T1, Sgh, Tgh
def euler(al):
    p = _PR[_PR <= M0].astype(np.float64)
    gp = (1 - 1/p)/(1 - p**(al - 1)); gam = gp**2 - 1                  # gamma(p) = g_p^2 - 1 >= 0
    # gamma(p) <= 2.2 p^{al-1} for p > 1e6 (al<=0.5); log(1+x) <= x
    lPg = _primesum(np.sum(np.log1p(gam/(p - 1))), 2 - al, 2.2*1.000001)
    lEg = _primesum(np.sum(np.log1p(gam*np.sqrt(p)/(p - 1))), 1.5 - al, 2.2*1.000001)
    lp = np.log(p)/(p - 1)
    # h^2 = sum_{p,p'|m} l_p l_p' : coefficient sums (tails: l_p/(p-1) <= 1.01 log p / p^2 -> crude: log p <= p^0.1 *10)
    S1 = np.sum(lp/(p - 1)) + _primesum(0, 1.9, 10.1); S2 = np.sum(lp**2/(p - 1)) + _primesum(0, 2.8, 14.0)
    T1 = np.sum(lp*np.sqrt(p)/(p - 1)) + _primesum(0, 1.4, 10.1); T2 = np.sum(lp**2*np.sqrt(p)/(p - 1)) + _primesum(0, 2.3, 14.0)
    return np.exp(lPg), np.exp(lEg), S1**2 + S2, T1**2 + T2

def rbar_factory(eps, A, b):
    N = U0; mu = _MU[:N + 1].astype(np.float64); k = _IDX[:N + 1]; k[0] = 1
    w = mu*k**(-1 - eps); w[0] = 0
    Ms = np.cumsum(w); Ls = np.cumsum(w*np.log(k))
    kk = np.arange(1, N); lo, hi = np.log(kk), np.log(kk + 1.0)
    # on [k,k+1): r(u) = (Ms_k - A) log u - Ls_k - A b ; endpoint values
    r_lo = (Ms[1:N] - A)*lo - Ls[1:N] - A*b; r_hi = (Ms[1:N] - A)*hi - Ls[1:N] - A*b
    pm = np.maximum(np.abs(r_lo), np.abs(r_hi)) + 1e-9
    big = 2*C_MU*exp1(eps*np.log(U0))
    suf = np.maximum.accumulate(pm[::-1])[::-1]; suf = np.maximum(suf, big)  # suf[k-1] = sup_{u>=k} |r|
    R1 = max(suf[0], A*b)
    def Rbar(u):
        u = np.asarray(u, float); out = np.empty_like(u)
        m1 = u < 1; out[m1] = R1 + A*np.log(1/u[m1])
        m2 = (u >= 1) & (u < U0); out[m2] = suf[np.floor(u[m2]).astype(np.int64) - 1]
        m3 = u >= U0; out[m3] = 2*C_MU*exp1(eps*np.log(u[m3]))
        return out
    return Rbar, R1

def bound(Lz, tau, kap, al=0.35, th=0.5, step=0.01, verbose=False):
    eps = kap/Lz; s = 1 + eps; mp.mp.dps = 30
    z = mp.zeta(s); zp = mp.zeta(s, derivative=1); A = float(1/z); b = float(-zp/z); zeta = float(z)
    Rbar, R1 = rbar_factory(eps, A, b)
    if al not in _GC: _GC[al] = (_gtab(al), euler(al))
    gt, (Pg, Eg, Ph, Eh) = _GC[al]
    L1, L2 = Lz, tau*Lz; Dl = L2 - L1
    def pieces(lm, gm, hm):  # lm = log m (array); returns F1, F2/g(m), F3/h(m)
        F1 = np.where(lm <= L1, A*Dl, np.where(lm <= L2, A*(L2 - lm + b), 0.0))
        F2 = np.zeros_like(lm)
        for Li in (L1, L2):
            ly = Li - lm; chi = ly >= 0; lyc = np.maximum(ly, 0)
            lD = th*lyc
            F1 = F1 + np.where(chi, Rbar(np.exp(lyc - lD)), 0.0)
            F2 = F2 + np.where(chi, R1*np.exp(-al*lD), 0.0)
        return F1, F2, 2*A*np.ones_like(lm)
    # exact part m <= M0
    m0 = M0 if m0 is None else m0
    if m0 > 1: assert L1 - np.log(m0) >= lU, "exact part needs y1 >= U0 for all m <= m0"
    m = np.arange(1, m0 + 1); sq = _MU[1:m0 + 1] != 0
    mm = m[sq].astype(float); lm = np.log(mm); W = mm**(-eps)/_PHI[1:M0 + 1][sq]
    F1, F2, F3 = pieces(lm, None, None)
    n1 = np.sum(W*F1**2); n2 = np.sum(W*(F2*gt[1:M0 + 1][sq])**2); n3 = np.sum(W*(F3*_H[1:M0 + 1][sq])**2)
    # intervals
    if np.log(M0) < L2:
        edges = np.arange(np.log(M0), L2, step); edges = np.unique(np.concatenate([edges, [L1, L2]]))
        edges = edges[edges >= np.log(M0)]
        a, bb = edges[:-1], edges[1:]
        f1a, f2a, _ = pieces(a, None, None); f1b, f2b, _ = pieces(bb - 1e-12, None, None)
        f1 = np.maximum(f1a, f1b); f2 = np.maximum(f2a, f2b)
        wa = np.exp(-eps*a); dl = bb - a; ea = 2*GERR/np.exp(a/2)
        n1 += np.sum(wa*f1**2*(dl + ea)); n2 += np.sum(wa*f2**2*(Pg*dl + Eg*ea))
        n3 += np.sum(wa*(2*A)**2*(Ph*dl + Eh*ea))
    pref = eps*zeta/Dl
    B = pref*(np.sqrt(n1) + np.sqrt(n2) + np.sqrt(n3))**2
    if verbose: print(f"  parts: main+R {pref*n1:.4f}  rankin {pref*n2:.4f}  hp {pref*n3:.5f}  R1={R1:.3f}")
    return B

def Binf(tau, k):
    from scipy.integrate import quad
    return k*(tau - 1)*(1 - np.exp(-k)) + quad(lambda v: np.exp(-k*v)*(1 + k*(tau - v))**2, 1, tau)[0]/(tau - 1)

if __name__ == "__main__":
    import sys
    for Lz in (np.log(1000), 30, 45, 55, 65):
        for tau, k in ((1.6, 0.3), (1.8, 0.3), (2.4, 0.3), (1.8, 0.5)):
            Bb = bound(Lz, tau, k, verbose=False)
            print(f"logz1={Lz:6.2f} tau={tau} kappa={k}:  bound={Bb:.4f}  Binf={Binf(tau,k):.4f}  eta={Bb/Binf(tau,k)-1:+.3f}")


# ---------------------------------------------------------------- version 2
# (a) full square for the main+Rankin pieces (no Minkowski between them): needs sum mu^2/phi * g (first moment);
# (b) small-y regime y_i < U0: N(y;q)/H = E_d[N(y/d;1)] lies in [lo(y), hi(y)], lo/hi = running min/max of N(u;1), u<=y
#     (with 0 for u<1)  -- no residue, no Rankin;
# (c) Rankin cut D chosen per interval (constant inside each interval -> endpoint maxima stay rigorous).
def euler1(al):
    p = _PR[_PR <= M0].astype(np.float64)
    g1 = (1 - 1/p)/(1 - p**(al - 1)) - 1                     # <= 1.1 p^{al-1} for p > 1e6
    lP = _primesum(np.sum(np.log1p(g1/(p - 1))), 2 - al, 1.1*1.000001)
    lE = _primesum(np.sum(np.log1p(g1*np.sqrt(p)/(p - 1))), 1.5 - al, 1.1*1.000001)
    return np.exp(lP), np.exp(lE)

def nrange_factory(eps):
    N = U0; mu = _MU[:N + 1].astype(np.float64); k = _IDX[:N + 1].copy(); k[0] = 1
    w = mu*k**(-1 - eps); w[0] = 0
    Ms = np.cumsum(w); Ls = np.cumsum(w*np.log(k))
    kk = np.arange(1, N); lo_, hi_ = np.log(kk), np.log(kk + 1.0)
    v1 = Ms[1:N]*lo_ - Ls[1:N]; v2 = Ms[1:N]*hi_ - Ls[1:N]
    mx = np.maximum.accumulate(np.maximum(np.maximum(v1, v2), 0)) + 1e-9
    mn = np.minimum.accumulate(np.minimum(np.minimum(v1, v2), 0)) - 1e-9
    def rng(y):  # y in [1,U0): range of N(u;1) over u <= y (index floor(y)-1 covers [1,floor(y)+1))
        i = np.clip(np.floor(y).astype(np.int64) - 1, 0, N - 2); return mn[i], mx[i]
    pmin = np.minimum(v1, v2) - 1e-9          # min of N on [k,k+1)
    sp = [pmin]
    while len(sp[-1]) > 1:
        q = sp[-1]; sp.append(np.minimum(q[:-1], q[1:len(q)] if False else np.minimum(q[:-1], q[1:])) if False else None); sp.pop()
        L = 1 << (len(sp) - 1)
        if 2*L > len(pmin): break
        sp.append(np.minimum(sp[-1][:-L], sp[-1][L:]))
    def wmin(u1, u2):   # min of N(u;1) over u in [u1,u2], 1 <= u1 <= u2 < U0
        i = int(np.floor(u1)) - 1; j = int(np.floor(u2)) - 1; i = max(i, 0); j = min(max(j, i), N - 2)
        k = int(np.log2(j - i + 1)); k = min(k, len(sp) - 1)
        return min(sp[k][i], sp[k][j - (1 << k) + 1])
    rng.wmin = wmin; rng.lo0 = mn[-1]
    return rng

DEBUG = None
def bound2(Lz, tau, kap, al=0.45, step=0.01, lDgrid=None, verbose=False, v3=True, m0=None):
    eps = kap/Lz; mp.mp.dps = 30; s = 1 + eps
    assert eps <= al, 'Rankin step needs 1 + (eps/R1) log(d/D) <= (d/D)^al: eps <= al suffices'
    z = mp.zeta(s); zp = mp.zeta(s, derivative=1); A = float(1/z); b = float(-zp/z); zeta = float(z)
    if v3:
        Rbar, R1 = rbar_factory4(eps, tau*Lz); A = Rbar.A; b = (1 + Rbar.c0)/A      # main1 = 1 + c0 + A log u
    else:
        Rbar, R1 = rbar_factory(eps, A, b)
    rng = nrange_factory(eps)
    key = ('v2', al)
    if key not in _GC: _GC[key] = (_gtab(al), euler(al), euler1(al), euler_h(al), addmom(al))
    gt, (Pg, Eg, Ph, Eh), (Pg1, Eg1), (S1h, T1h, Sgh, Tgh), AM = _GC[key]
    (Sk, Tk), (Pk, Ek), (Phk, Ehk), (Sgk, Tgk) = AM['k'], AM['kk'], AM['hk'], AM['gk']
    L1, L2 = Lz, tau*Lz; Dl = L2 - L1; lU = np.log(U0)
    if lDgrid is None: lDgrid = np.linspace(0, 60, 241)
    def errpart(ly):  # large-y error for one side, ly >= lU (scalar): choose lD minimising R(y/D) + R1 D^-al * rms(g)
        cand = Rbar(np.exp(np.maximum(ly - lDgrid, 0))) + R1*np.exp(-al*lDgrid)*np.sqrt(Pg)
        j = np.argmin(np.where(lDgrid <= ly, cand, np.inf)); return lDgrid[j]
    def F(a, b_, lDs):
        """sup over log m in [a, b_] of the pieces (F1 smooth, F2 x g(m), F3 x h(m)).  The interval lies in ONE regime
        (edges include L1, L2, L1-lU, L2-lU; regime read at the midpoint), and every ingredient is monotone in m:
        Rbar nonincreasing in its argument (worst at b_), main terms A(ly+b) decreasing in m (worst at a for upper,
        b_ for lower), running max/min ranges monotone in y (worst at a), window minimum over the union window."""
        mid = 0.5*(a + b_); ly1m, ly2m = L1 - mid, L2 - mid
        def large(lyb, i):          # lyb = smallest log y on the interval
            lD = lDs[i]; return Rbar(np.exp([max(lyb - lD, 0)]))[0], R1*np.exp(-al*lD)
        if ly2m < 0: return 0.0, 0.0, 0.0, 0.0
        if ly1m >= 0:   # m <= z1: Yt = (N2 - N1)/H
            e2, c2 = large(L2 - b_, 2)
            if ly1m >= lU:
                e1, c1 = large(L1 - b_, 1)
                return A*Dl + e1 + e2, c1 + c2, 2*eps, 0.0
            lo, hi = rng(np.exp(L1 - a)); mn2hi = A*(L2 - a + b); mn2lo = A*(L2 - b_ + b)
            lDp = lDs[3]; nmin = max(rng.wmin(max(np.exp(L1 - b_ - lDp), 1.0), np.exp(L1 - a)), 0.0)
            dev = max(mn2hi - nmin - min(lo, 0), hi - mn2lo)
            optA = (dev + e2, c2 + nmin*np.exp(-al*lDp), 2*eps, 0.0)
            # option C: d = 1 term only;  P(d=1) = 1/H >= 1 - k(m), k(m) = sum_{p|m} 1/p, so
            #   N1/H >= n1 - (1 - 1/H)(n1 - lo') >= n1 - k(m)(n1 - lo'),  n1 = min N over [y1(b), y1(a)], lo' = min(lo,0)
            lo_ = min(lo, 0); n1 = max(rng.wmin(max(np.exp(L1 - b_), 1.0), np.exp(L1 - a)), lo_)
            optC = (max(mn2hi - n1, hi - mn2lo) + e2, c2, 2*eps, n1 - lo_)
            if lDs[5] == 'C': return optC
            if lDs[5] is True:   # option B (decided once per interval)
                lD = lDs[4]; e1 = Rbar(np.exp([max(L1 - b_ - lD, 0)]))[0]; c1 = R1*np.exp(-al*lD)
                return (A*Dl + e1 + e2, c1 + c2, 2*eps, 0.0)
            return optA
        # z1 < m <= z2: Yt = N2/H
        if ly2m >= lU:
            e2, c2 = large(L2 - b_, 2); return A*(L2 - a + b) + e2, c2, 2*eps, 0.0
        lo, hi = rng(np.exp(L2 - a)); return max(-lo, hi), 0.0, 0.0, 0.0
    def choose(a_, b_):
        lDs = {1: errpart(L1 - b_) if L1 - b_ >= lU else 0, 2: errpart(L2 - b_) if L2 - b_ >= lU else 0, 3: 0.0, 4: 0.0, 5: False}
        if 0 <= L1 - b_ and L1 - a_ < lU:   # small-y1 regime: choose D' (option A) and D (option B), then the option
            ly1 = L1 - b_; best = None
            for lDp in np.linspace(0, max(ly1, 0), 25):
                nm = max(rng.wmin(max(np.exp(ly1 - lDp), 1.0), np.exp(L1 - a_)), 0.0)
                v = A*(L2 - a_ + b) - nm + nm*np.exp(-al*lDp)*np.sqrt(Pg)
                if best is None or v < best[0]: best = (v, lDp)
            lDs[3] = best[1]
            g_ = np.linspace(0, max(ly1, 0), 41)
            cand = Rbar(np.exp(np.maximum(ly1 - g_, 0))) + R1*np.exp(-al*g_)*np.sqrt(Pg)
            lDs[4] = g_[np.argmin(cand)]
            sc = lambda f: f[0] + f[1]*np.sqrt(Pg) + f[3]*np.sqrt(Pk)
            fA = F(a_, b_, lDs); lDs[5] = True; fB = F(a_, b_, lDs); lDs[5] = 'C'; fC = F(a_, b_, lDs)
            lDs[5] = min(((sc(fA), False), (sc(fB), True), (sc(fC), 'C')), key=lambda t: t[0])[1]
        return lDs
    # exact part m <= M0 (y's huge): vectorised large-y formulas with lD per m from a coarse lookup
    m0 = M0 if m0 is None else m0
    if m0 > 1: assert L1 - np.log(m0) >= lU, "exact part needs y1 >= U0 for all m <= m0"
    m = np.arange(1, m0 + 1); sq = _MU[1:m0 + 1] != 0
    mm = m[sq].astype(float); lm = np.log(mm); W = mm**(-eps)/_PHI[1:m0 + 1][sq]
    g = gt[1:m0 + 1][sq]; h = _H[1:m0 + 1][sq]
    lD1 = errpart(L1 - np.log(max(m0, 1))) if L1 - np.log(max(m0,1)) >= lU else 0; lD2 = errpart(L2 - np.log(max(m0,1)))   # same D for all m <= M0 (fine: rigorous for any D)
    e1 = Rbar(np.exp(L1 - lm - lD1)); e2 = Rbar(np.exp(L2 - lm - lD2)); c = R1*(np.exp(-al*lD1) + np.exp(-al*lD2))
    s12 = np.sum(W*(A*Dl + e1 + e2 + c*g + 2*eps*h)**2); s3 = 0.0     # full square incl. the Hp piece (2 eps h >= 2 A h)
    if DEBUG is not None: DEBUG.append(('M0', s12))
    if m0 == -1:   # TEST MODE (small z2 <= M0): every squarefree m evaluated individually, same per-m logic
        tot = 0.0
        for mi in range(1, int(np.exp(L2)) + 1):
            if _MU[mi] == 0: continue
            a_ = np.log(mi); b_ = a_ + 1e-9
            lDs = choose(a_, a_)
            f1, f2, f3, f4 = F(a_, a_, lDs)
            gm = np.prod([(1 - 1/p)/(1 - p**(al - 1)) for p in _PR[:200] if mi % p == 0]) if True else 1
            # exact g, h for this m
            gg = 1.0; hh = 0.0; kk_ = 0.0; r = mi
            for p in _PR:
                if p*p > r: break
                if r % p == 0: gg *= (1 - 1/p)/(1 - p**(al - 1)); hh += np.log(p)/(p - 1); kk_ += 1/p; r //= p
            if r > 1: gg *= (1 - 1/r)/(1 - r**(al - 1)); hh += np.log(r)/(r - 1); kk_ += 1/r
            tot += mi**(-eps)/_PHI[mi]*(f1 + f2*gg + f3*hh + f4*kk_)**2
            if DEBUG is not None: DEBUG.append((mi, f1 + f2*gg + f3*hh + f4*kk_))
        return eps*zeta/Dl*tot
    if m0 <= 1: s12 = 0.0; s3 = 0.0     # (m = 1 is covered by the interval sums starting at a = 1)
    edges = np.arange(np.log(max(m0, 1)), L2, step); edges = np.unique(np.concatenate([edges, [L1, L2, L1 - lU, L2 - lU]]))
    edges = edges[(edges >= np.log(max(m0, 1))) & (edges <= L2)]
    for a_, b_ in zip(edges[:-1], edges[1:]):
        lDs = choose(a_, b_)
        f1, f2, f3, f4 = F(a_, b_, lDs)
        wa = np.exp(-eps*a_); dl = b_ - a_; ea = 2*GERR*np.exp(-a_/2)
        _c = wa*(f1**2*(dl + ea) + 2*f1*f2*(Pg1*dl + Eg1*ea) + f2**2*(Pg*dl + Eg*ea)
                 + f3**2*(Ph*dl + Eh*ea) + 2*f1*f3*(S1h*dl + T1h*ea) + 2*f2*f3*(Sgh*Pg1*dl + Tgh*Eg1*ea)
                 + f4**2*(Pk*dl + Ek*ea) + 2*f1*f4*(Sk*dl + Tk*ea) + 2*f2*f4*(Sgk*Pg1*dl + Tgk*Eg1*ea)
                 + 2*f3*f4*(Phk*dl + Ehk*ea)); s12 += _c
        if DEBUG is not None: DEBUG.append((a_, _c, wa*f1**2*dl, f1, f2))
    pref = eps*zeta/Dl
    if verbose: print(f"  main+rankin {pref*s12:.4f}  hp {pref*s3:.5f}")
    return pref*(np.sqrt(s12) + np.sqrt(s3))**2


# ---------------------------------------------------------------- version 3: main term 1 + eps log u, via m-check
# N_eps(u;1) = sum mu(k) k^{-1-eps} log(u/k).  With G(x) = x e^{eps x}, two partial summations give exactly (u >= 1)
#   N_eps(u) = 1 + eps log u + r(u),  r(u) = u^{-eps} rho(u) + eps int_1^u rho(t) t^{-eps} (2 + eps log(u/t)) dt/t,
#   rho(t) = m-check(t) - 1,  m-check(t) = sum_{n<=t} mu(n)/n log(t/n).
# Ramare (Math. Comp. 2015): |rho(t)| <= 0.00257/log t for t >= 3846.  rho is exact (eps-free) below U0.
# The convolution identity works for ANY main term linear in log u, so main_q(y)/H = 1 + eps(log y + Hp).
import os as _os
# |mcheck(t) - 1| <= RHO_C/log t.  CONFIRMED (abstract of Ramare, Math. Comp. 84 (2015)): 3/14 for all t > 1.
# The sharper 0.00257 (t >= 3846) is quoted on the TME-EMT page from Ramare's work but NOT yet located in the paper;
# select it with LEMMAF_RHO=0.00257 (conditional result).
RHO_C = float(_os.environ.get('LEMMAF_RHO', 3/14))
def rbar_factory3(eps, Lmax, shift=True):
    N = U0; mu = _MU[:N + 1].astype(np.float64); k = _IDX[:N + 1].copy(); k[0] = 1
    w = mu/k; w[0] = 0
    M = np.cumsum(w); Lw = np.cumsum(w*np.log(k))
    kk = np.arange(1, N).astype(np.float64); x0, x1 = np.log(kk), np.log(kk + 1)
    Mk, Lk = M[1:N], Lw[1:N]
    r0 = Mk*x0 - Lk - 1; r1 = Mk*x1 - Lk - 1                      # rho at piece ends (linear in x on the piece)
    sup = np.maximum(np.abs(r0), np.abs(r1)) + 1e-9
    I0 = Mk*(x1**2 - x0**2)/2 - (Lk + 1)*(x1 - x0)                   # int rho dx  (eps = 0)
    I2 = Mk*(x1**3 - x0**3)/3 - (Lk + 1)*(x1**2 - x0**2)/2           # int rho x dx
    A0 = sup*(x1**2 - x0**2)/2; A2 = sup*(x1**3 - x0**3)/3          # int |rho| x dx, int |rho| x^2 dx (for |t^-eps - 1| <= eps x)
    J0 = np.concatenate([[0], np.cumsum(I0)]); J2 = np.concatenate([[0], np.cumsum(I2)])
    E0 = np.concatenate([[0], np.cumsum(A0)]); E2 = np.concatenate([[0], np.cumsum(A2)])
    jinf = J0[-1] if shift else 0.0          # main1(u) = 1 + 2 eps jinf + eps log u  (constant shift dc = 2 eps jinf)
    # bound at the left end of piece j (u in [k_j, k_j+1)): cumulative to k_j, plus the partial piece (<= sup*dx*(2+eps lu))
    lu = x1                                                             # use the larger log u on the piece (worst)
    core = np.abs(2*(J0[:-1] - jinf) + eps*(lu*J0[:-1] - J2[:-1]))
    slack = (2*eps*E0[:-1] + eps*(lu*eps*E0[:-1] + eps*E2[:-1]) + sup*(x1 - x0)*(2 + eps*lu)
             + eps*np.abs(J0[:-1])*(x1 - x0) + eps*sup*x1*(x1 - x0) + 1e-9*(x1 + 2))   # log u vs x1; partial-piece J2; float
    bnd = sup + eps*(core + slack)                                      # |u^-eps rho(u)| <= sup on the piece
    # tail v >= U0 (any v up to e^Lmax): rho_bar = 0.00257/log t
    lU = np.log(U0); J0e, J2e, E0e, E2e = J0[-1], J2[-1], E0[-1], E2[-1]
    tailint = RHO_C*np.log(Lmax/lU)                                     # int_{U0}^{v} rho_bar dt/t
    tail = (RHO_C/lU + eps*(abs(2*(J0e - jinf)) + eps*(Lmax*abs(J0e) + abs(J2e)))
            + eps*(2*E0e + eps*(Lmax*E0e + E2e)) + eps*(2 + eps*Lmax)*tailint
            + eps*eps*RHO_C*Lmax + 1e-9*(Lmax + 2))                              # tail of the J2 integral; float
    suf = np.maximum.accumulate(bnd[::-1])[::-1]; suf = np.maximum(suf, tail)
    dc = 2*eps*jinf
    R1 = max(suf[0], 1.0 + abs(dc))
    def Rbar(u):
        u = np.asarray(u, float); out = np.empty_like(u)
        m1 = u < 1; out[m1] = R1 + eps*np.log(1/u[m1])
        m2 = (u >= 1) & (u < U0); out[m2] = suf[np.floor(u[m2]).astype(np.int64) - 1]
        m3 = u >= U0; out[m3] = RHO_C/np.log(u[m3]) + (tail - RHO_C/lU)   # sup_{v>=u} v^-eps rho_bar(v) <= RHO_C/log u
        return out
    Rbar.tail = tail; Rbar.J0 = J0e; Rbar.dc = dc
    return Rbar, R1


# ---------------------------------------------------------------- version 4 (Oct 10): own rho bound, exact eps-twist
# (a) rho(x) = mcheck(x) - 1 satisfies, for x >= 2e6 (proof in LEMMA_F.md, "Lemma R"):
#       rho = (1-gamma) m(x) - M(x)/(2x) - 1/x - T(x),   T = sum_{n<=x} mu(n)/n eta(x/n),
#       eta(y) = H_[y] - log y - gamma - (1/2 - {y})/y  (continuous, |eta'| = |[y]+1/2-y|/y^2 <= 1/(2y^2)),
#     from sum_{n<=x} mu(n)/n H_[x/n] = 1 and sum mu(n)[x/n] = 1;  T by Stieltjes against m(t) with |m(t)| <= 0.010032/log t
#     (t >= 617990), |M(x)| <= 0.006688 x/log x (x >= 1798118) [arXiv 2408.05969 Cor. 1.2, Thm 1.1], exact m(t), t < 1e6.
#     Result (rho_own.py): |rho(x)| <= 0.0121/log x for x >= 2e6.   (With Ramare's published 1/69 and 0.013: 0.0187.)
# (b) the t^{-eps} weights in the rho-integrals are kept exactly (midpoint value per piece + explicit error), and the
#     eps^2 log u drift of r is absorbed into the main term (any linear main term is allowed):
#       main1(u) = 1 + c0 + A log u,  A = eps (1 + eps j),  c0 = eps (2 j - eps j2),  j = J0eps(U0), j2 = J2eps(U0).
RHO_OWN = 0.0121
def rbar_factory4(eps, Lmax):
    N = U0; mu = _MU[:N + 1].astype(np.float64); k = _IDX[:N + 1].copy(); k[0] = 1
    w = mu/k; w[0] = 0
    M = np.cumsum(w); Lw = np.cumsum(w*np.log(k))
    kk = np.arange(1, N).astype(np.float64); x0, x1 = np.log(kk), np.log(kk + 1); dx = x1 - x0
    Mk, Lk = M[1:N], Lw[1:N]
    r0 = Mk*x0 - Lk - 1; r1 = Mk*x1 - Lk - 1
    sup = np.maximum(np.abs(r0), np.abs(r1)) + 1e-9
    wm = np.exp(-eps*(x0 + x1)/2)                                   # |e^{-eps x} - wm| <= eps dx/2 on the piece
    I0 = (Mk*(x1**2 - x0**2)/2 - (Lk + 1)*dx)*wm
    I2 = (Mk*(x1**3 - x0**3)/3 - (Lk + 1)*(x1**2 - x0**2)/2)*wm
    e0 = eps*dx/2*sup*dx; e2 = eps*dx/2*sup*(x1**2 - x0**2)/2      # errors of I0, I2 vs the exact eps-twisted integrals
    J0 = np.concatenate([[0], np.cumsum(I0)]); J2 = np.concatenate([[0], np.cumsum(I2)])
    E0 = np.concatenate([[0], np.cumsum(e0)]); E2 = np.concatenate([[0], np.cumsum(e2)])
    j, j2 = J0[-1], J2[-1]
    A = eps*(1 + eps*j); c0 = eps*(2*j - eps*j2)
    # r'(u) = u^-eps rho(u) + eps[ 2(J0(u)-j) + eps lu (J0(u)-j) - eps (J2(u)-j2) ]  (u in [1,U0), cumulative to the piece)
    lu = x1
    core = np.abs(2*(J0[:-1] - j) + eps*(lu*(J0[:-1] - j) - (J2[:-1] - j2)))
    # errors: cumulative E's (both the cumulative and the totals j, j2 carry them -> factor 2), partial piece, log u vs x1
    slack = (2*(2*E0[:-1] + eps*lu*2*E0[:-1] + eps*2*E2[:-1])
             + sup*dx*(2 + eps*lu) + eps*sup*x1*dx + eps*np.abs(J0[:-1] - j)*dx + 1e-9*(x1 + 2))
    bnd = sup + eps*(core + slack)
    lU = np.log(U0)
    # u >= U0: |rho| <= RHO_OWN/log t; |J0(u)-j| <= RHO_OWN log(lu/lU); |J2(u)-j2| <= RHO_OWN (lu - lU)
    tail_rest = (eps*((2 + eps*Lmax)*RHO_OWN*np.log(Lmax/lU) + eps*RHO_OWN*(Lmax - lU))
                 + eps*(2*2*E0[-1] + eps*Lmax*2*E0[-1] + eps*2*E2[-1]) + 1e-9*(Lmax + 2))
    tail0 = RHO_OWN/lU + tail_rest
    suf = np.maximum.accumulate(bnd[::-1])[::-1]; suf = np.maximum(suf, tail0)
    R1 = max(suf[0], 1.0 + abs(c0))
    def Rbar(u):
        u = np.asarray(u, float); out = np.empty_like(u)
        m1 = u < 1; out[m1] = R1 + eps*np.log(1/u[m1])
        m2 = (u >= 1) & (u < U0); out[m2] = suf[np.floor(u[m2]).astype(np.int64) - 1]
        m3 = u >= U0; out[m3] = RHO_OWN/np.log(u[m3]) + tail_rest
        return out
    Rbar.tail = tail0; Rbar.A = A; Rbar.c0 = c0; Rbar.dc = c0
    return Rbar, R1
