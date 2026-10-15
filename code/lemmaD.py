"""Lemma D: Bellotti's log-free density argument (arXiv:2405.12545, Secs. 3-4) re-run with
 (D1) a log-free pointwise mollifier bound on Re s = 1/2:
      |M(1/2+it)| <= sum_{v<=V, w<=W} [v,w]^{-1/2} <= 4 zeta(3/2) sqrt(VW) = 10.45 sqrt(VW)
      (write v = g a, w = g b, (a,b) = 1; |psi|,|theta| <= 1; sum_{a<=A} a^{-1/2} <= 2 sqrt A), replacing the
      d2 log^4(VW) sqrt(VW) behind her B-term, together with |zeta(1/2+iu)| <= 0.618 |u|^{1/6} log|u| (HPY);
 (D2) parameters (u, v, w, x) chosen separately for every (delta-bin, log T-bin) -- the theorem holds for any
      admissible choice (w < u, u + v < x, X > U/L^2, c1^2 - c5 > 0), so per-sigma optimisation is legitimate,
      and it also covers sigma below Bellotti's alpha_0;
 (D3) Lemma A (rho, n~ = 1) and Lemma C (c4 crowding) on top.
All other constants are Bellotti's formulas (reconstruction checked against her Tables 5-6: d4 91.9 vs 90.96,
d5 1.109 vs 1.110, c5 0.063 vs 0.066, c1 0.983 vs 0.984 at [e^34, e^35])."""
import math, functools
from scipy.optimize import minimize
import pairsum, lemmaC
R0 = 4.896      # Bellotti-Trudgian-Yang, arXiv:2603.21490 (t >= 3); was 5.558691 (MTY)
IG = 5.922           # int |Gamma(-1/2 + it)| (1+|t|) dt (worst sigma in [-1/2, 1/2 - alpha]); mpmath quad
import os
USE_E = not os.environ.get('NO_LEMMA_E')
def d4_bellotti(u, t, x, L):
    v = t*u
    F = (1.084*(t + 1) + 1.301*(1 + t*t) - 0.116)/(t - 1)
    d3 = 4.0244*(u*L + 1.333)*L/((v - u)**2*L*L)
    return 3.09*(u + v)/(v - u)*F + d3*(2/L + (x - u - v))
USE_F = not os.environ.get('NO_LEMMA_F')
F_SAFE = 1.03     # table interpolation safety factor; certified a posteriori by verify_F.py (exact lemmaF.bound2 at the chosen point)
_FT = None
def _ftab():
    global _FT
    if _FT is None:
        import pickle, numpy as np
        d = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lemmaF_table.pkl'), 'rb'))
        LZ, TAU, KAP = d['LZ'], d['TAU'], d['KAP']
        arr = np.array([[[d['B'][(a, b, c)] for c in KAP] for b in TAU] for a in LZ])
        from scipy.interpolate import RegularGridInterpolator
        f = RegularGridInterpolator((np.array(LZ, float), np.array(TAU), np.log(KAP)), arr, bounds_error=False, fill_value=None)
        _FT = (f, (LZ[0], LZ[-1]), (TAU[0], TAU[-1]), (KAP[0], KAP[-1]))
    return _FT
def B_F(Lz, tau, kap):
    """Lemma F: eps log(z2/z1) sum Psi^2 n^{-1-eps} <= B_F (table value x F_SAFE); None outside the table."""
    f, lz, tt, kk = _ftab()
    if not (lz[0] <= Lz <= lz[1] and tt[0] <= tau <= tt[1] and kk[0] <= kap <= kk[1]): return None
    import numpy as np
    return F_SAFE*float(f([[Lz, tau, np.log(kap)]])[0])
# ---- Lemma G (rigorous version, Oct 13): smoothed weights + contour argument; see LEMMA_G.md "Step 3" ----
G_K1_UP = float(os.environ.get('LEMMAG_K1', '-0.5853'))   # upper end of the certified enclosure of k1 (kcert2.py)
G_K0_UP = float(os.environ.get('LEMMAG_K0', '0.4214'))     # upper end of the certified enclosure of k0
G_CT, G_ETA = 2.5, 0.02                                  # S_t <= G_CT lX eta^2/Delta^2 (stbound.py, sup over the grid)
USE_G = os.environ.get('LEMMAG', '0') == '1'
_REM = None
def d4_bound_G(u, t, x, L):
    """sum_n Psi^2 e^{-n/X}/n - e^{-1/X} <= D4m + REM + 2 sqrt(S_s S_t) + S_t, D4m = Delta^-2[(Y2^2-Y1^2)/2 - gamma Delta
    + k1 (Y1+Y2-2 gamma) + 2 k0];  REM = explicit contour remainders (lemG/remainder.py); S_t from the smoothing (stbound.py)."""
    global _REM
    lU, lV, lX = u*L, t*u*L, x*L; Y1, Y2, D = lX - lV, lX - lU, lV - lU; g = 0.5772156649
    if not (lU >= 40 and lX > lV): return None
    if _REM is None:
        import sys as _s; _s.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lemG'))
        import remainder as _r, stbound as _st; _REM = (_r, _st)
    D4m = ((Y2*Y2 - Y1*Y1)/2 - g*D + G_K1_UP*(Y1 + Y2 - 2*g) + 2*G_K0_UP)/(D*D)
    rem = _REM[0].bounds(u, t, x, L)[0]
    Ss = math.exp(-1/math.exp(lX)) + D4m + rem; St = _REM[1].St_rig(lU, lV, lX, G_ETA)
    return D4m + rem + 2*math.sqrt(Ss*St) + St + math.exp(-1/math.exp(x*L)) - math.exp(-L*L/math.exp(u*L))
def d4_bound(u, t, x, L, ret_a=False):
    v1 = d4_bound_rza(u, t, x, L, ret_a=True)
    if USE_G:
        gv = d4_bound_G(u, t, x, L)
        if gv is not None and gv < v1[0]: v1 = (gv, 0.0, 'G')
    if not USE_F: return v1 if ret_a else v1[0]
    best = v1
    import numpy as np
    f, lz, tt, kk = _ftab(); Lz = u*L
    if lz[0] <= Lz <= lz[1] and tt[0] <= t <= tt[1]:
        a = np.arange(5, 400, 2)*0.005; eps = a/(x*L); kap = eps*u*L
        ok = (kap >= kk[0]) & (kap <= kk[1])
        if ok.any():
            a, eps, kap = a[ok], eps[ok], kap[ok]
            B = F_SAFE*f(np.column_stack([np.full_like(kap, Lz), np.full_like(kap, t), np.log(kap)]))
            S = B/(eps*(t - 1)*u*L)
            val = (eps/math.e)**eps*np.exp(a)*(S - 1) + math.exp(-1/math.exp(x*L)) - math.exp(-L*L/math.exp(u*L))
            j = int(np.argmin(val))
            if val[j] < best[0]: best = (float(val[j]), float(a[j]), 'F')
    return best if ret_a else best[0]
def d4_bound_rza(u, t, x, L, ret_a=False):
    """Lemma E: sum_n Psi(n)^2 w(n)/n <= sum_n Psi(n)^2 e^{-n/X}/n <= (eps/e)^eps X^eps sum_n Psi(n)^2 n^{-1-eps}
    (e^{-y} <= (eps/e)^eps y^{-eps}), and Ramare-Zuniga Alterman (arXiv:2405.12662, Thm 1.1; z1 = U >= 100,
    z2 = V = U^t): sum Psi^2 n^{-1-eps} <= e^{gamma eps}/(eps log(V/U)) [A(t+1) + (B - C U^-eps + (B - C V^-eps) t^2) eps log U]/(t-1),
    A = 1.084, B = 1.301, C = 0.318.  eps = a/log X optimised over a."""
    A, B, C, g = 1.084, 1.301, 0.318, 0.5772156649
    best = 1e300
    for k in range(5, 400):
        a = k*0.005; eps = a/(x*L)
        val = (eps/math.e)**eps * math.exp(a) * math.exp(g*eps)/(eps*(t - 1)*u*L) \
              * (A*(t + 1) + (B - C*math.exp(-eps*u*L) + (B - C*math.exp(-eps*t*u*L))*t*t)*eps*u*L)/(t - 1)
        val -= (eps/math.e)**eps * math.exp(a)          # n = 1 term (Psi(1)^2 = 1) is exactly 1 in sum Psi^2 n^{-1-eps};
        val += math.exp(-1/math.exp(x*L)) - math.exp(-L*L/math.exp(u*L))   # ... its true weight is w(1)
        if ret_a:
            if best == 1e300 or val < best[0]: best = (val, a, 'RZA')
        else: best = min(best, val)
    return best
# crowding budget: sum_rho Re 1/(s-rho) = Re zeta'/zeta(s) + Re 1/(s-1) + Re psi(s/2+1)/2 - log(pi)/2
#   <= -zeta'/zeta(1+eta) + (1/2) log(t/(2 pi)) + O(1/t),   and  -zeta'/zeta(1+eta) <= 1/eta - GAP for 0 < eta <= 1.5
GAP = 0.36
# positivity with P(theta) = (c + cos theta)^2 = (c^2 + 1/2) + 2c cos theta + (1/2) cos 2 theta >= 0 gives also
#   sum_rho Re 1/(s-rho) <= (1/2) log(t/2pi) + ((c^2+1/2)/(2c)) (-zeta'/zeta(1+eta)) + (1/(4c)) (1/2) log(2t/2pi)
@functools.lru_cache(maxsize=None)
@functools.lru_cache(maxsize=None)
def _mzl(eta):
    """Oct 13: exact -zeta'/zeta(1+eta) (mpmath, 30 digits, + 1e-12 margin) instead of 1/eta - GAP."""
    import mpmath as _mp
    _mp.mp.dps = 30
    v = float(-_mp.zeta(1 + _mp.mpf(eta), derivative=1)/_mp.zeta(1 + _mp.mpf(eta)))
    return min(v + 1e-12, 1/eta - GAP)
def _M(logT, eta):
    F = _mzl(round(eta, 12))
    m = 0.5*(logT - math.log(2*math.pi)) + F
    for k in range(50, 201, 2):
        c = k/100
        m = min(m, 0.5*(logT - math.log(2*math.pi)) + (c*c + 0.5)/(2*c)*F + (logT - math.log(math.pi))/(8*c))
    return m
def _cost(d, eta, de): return min(a/(a*a + d*d) for a in (eta, eta + de))
@functools.lru_cache(maxsize=None)
def _gr(d, de):
    """Oct 13: |Gamma(z)| = |Gamma(1+z)|/|z| and, for Re z = s in [0, 2 delta], |Gamma(1+s+i d)| <= |Gamma(1+2delta+i d)|/Gamma(1+2delta)
    (|Gamma(1+s+id)|^2 = Gamma(1+s)^2 prod_n (1 + d^2/(1+s+n)^2)^{-1}: the product increases with s; Gamma(1+s) <= 1).
    Capped at 1 (the overall factor 1.0067 is kept outside); decreasing in d, so evaluating at the minimal distance is an upper bound."""
    import mpmath as _mp
    v = float(abs(_mp.gamma(_mp.mpc(1 + 2*de, d)))/_mp.gamma(1 + 2*de))
    return min(1.0, v*(1 + 1e-12)) if d > 0 else 1.0
_DST = None
def _dstar(logT, de):
    """Two-zero repulsion (rep/rep.py, dtable.py): two zeros with 1-beta <= de at heights <= e^logT (>= 3e12) are at least
    dstar apart.  0 where not computed.  Capped at 1."""
    global _DST
    if os.environ.get('NO_REP'): return 0.0
    if _DST is None:
        import json
        import glob
        _DST = {}
        for f in sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rep', 'dstar2c_*.json'))):
            _DST.update(json.load(open(f)))      # v2 (Oct 14): exact W, re-derived BTY constants, height-uniform cells
        # a zero counted at T = e^lk can lie at ANY height in [3e12, e^lk]: use min over the cells l <= lk (same delta bin);
        # cells where delta <= lambda*(l)/l are zero-free and impose nothing.  (rep/dstar2c_* already hold this minimum:
        # the only cells it changes are (32..37, 121) -> 0, set in the files.)
    return min(1.0, _DST.get('%d,%d' % (math.ceil(logT - 1e-9), int(round(de*1e4))), 0.0))   # half-height cells use the integer cell above (superset)
def _lp2(h, de, logT, eta):
    """max sum_s b(s) over representatives s != r with |gamma_s - gamma_r| < ~1, subject to the crowding budget at
    s0 = 1 + eta + i gamma_r:  sum_s cost <= (1/2) log T + 1/eta - cost(r).  On each side the representative of the
    k-th rectangle away lies at distance d in ((k-1)h, (k+1)h): contribution (in diagonal units / (2 delta))
    b <= 1/max(2 delta, (k-1)h)  [|Gamma(z)| X^{Re z} <= 1.0067 X^{2 delta}/max(|Delta|, 2 delta) since X^y/y increases
    on [2 nu, 2 delta]], cost >= cost((k+1)h)."""
    M = _M(logT, eta) - _cost(0, eta, de)
    dst = _dstar(logT, de)          # Oct 14 repulsion lemma (rep/): no two counted zeros within distance dst
    items = []
    for k in range(1, int(1/h) + 2):
        if (k + 1)*h <= dst: continue
        dmin = max((k - 1)*h, dst); b = _gr(round(dmin, 9), de)/max(2*de, dmin)
        c = _cost((k + 1)*h, eta, de); items += [(b/c, b, c)]*2
    items.sort(reverse=True); tot = 0.0
    for _, b, c in items:
        if M <= 0: break
        f = min(1.0, M/c); tot += f*b; M -= f*c
    return tot
def _ntil(h, de, logT):
    best = 1e9
    for k in range(1, 800):
        eta = k*0.0005; best = min(best, _M(logT, eta)/_cost(h/2, eta, de))
    return max(1, math.floor(best))
@functools.lru_cache(maxsize=None)
def pair_coef(de, L):
    best = None
    for kap in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0):
        h = kap*de; nt = _ntil(h, de, L)
        lp = min(_lp2(h, de, L, k*0.002) for k in range(1, 200))
        val = nt*1.0067*(1 + 2*de*lp)
        if best is None or val < best[0]*best[1]: best = (val/nt, nt)
    return best[0], best[1]
def consts(p, al, L):
    u, t, w, x = p; v = t*u; de = 1 - al
    if not (0.02 < w < u and t > 1.05 and x*L > u*L - 2*math.log(L)): return None
    if not USE_E and not (u + v < x): return None
    B = 10.45*0.618*L*IG/(2*math.pi) * math.exp(L*((v + w)/2 + 1/6 + x*(0.5 - al)))
    tail = 3.8*math.exp(-L*(L - 2*u))*(L*L + 1)/(w*(v - u)*L**6)
    c1 = math.exp(-L*L/math.exp(u*L)) - B - tail
    d4 = d4_bound(u, t, x, L) if USE_E else d4_bellotti(u, t, x, L)
    d5 = 1.0061*(1 + (1.333 + 4.0*math.exp(-w*L/2))/(w*L))
    # J (Oct 14): |zeta(1/2+iu)| <= max(1.461, 0.618|u|^{1/6} log|u|) (HPY, |u| <= e^105) with |u| <= T + |t|, so the whole
    # (valid for ALL |u| >= 3: for |u| > e^105 Patel-Yang 66.7 u^{27/164} <= 0.618 u^{1/6} log u since 0.618 u^{0.00203} log u >= 80.3 > 66.7 at u = e^105, increasing)
    # integral (4.16) is <= Z(T) int |Gamma(c+it)| (1+|t|/3)^{1/6}(1 + log(1+|t|/3)/log 3) dt <= 4.72 Z(T), c <= -0.44.
    J = 4.72*0.618*math.exp(L/6)*L
    c5 = d4*d5/(2*math.pi)*(math.exp(w*L)/(w*L))*(math.exp(L*x*(1.5 - 2*al)) + math.exp((L*u - 2*math.log(L))*(1.5 - 2*al)))*J
    if c1 <= 0 or c1*c1 - c5 < 0.2: return None
    return d4*d5/w, c1*c1 - c5
@functools.lru_cache(maxsize=None)
def best(de_key, lt_key):
    """minimise log N(1 - delta, T) for delta in the bin ((de_key-1)e-4, de_key e-4], log T in (lt_key-1, lt_key]."""
    de = de_key/1e4; L = float(lt_key)                 # constants at the bin top (worst), exponent at bin top
    al = 1 - de                                         # sup over beta >= al: worst al is the bin's smallest sigma
    import ramanujan_upper as _ru
    nu = 1/(_ru.Rof(L)*L)
    if de <= nu: return None
    # pair sum over |gamma_r - gamma_s| < 1 and zeros per rectangle (Lemma A''): rectangles of height h = kappa*delta,
    # one representative each; no spacing assumption between representatives of adjacent rectangles
    coef, ntil = pair_coef(de, L)
    I = lemmaC.Iprime(de_key, lt_key)
    coef = coef + 2*math.sqrt(2*math.pi)*2.4*I
    def f(p):
        r = consts(p, al, L)
        if r is None: return 1e9
        dd, den = r
        return math.log(dd*coef*ntil/den) + 2*p[3]*de*L
    starts = [(u, t, w, x) for u in (1.6, 2.0) for t in (2.0, 2.6) for w in (0.35,) for x in (5.5, 6.8)] + [(2.3, 1.8, 0.65, 5.8), (2.0, 1.75, 0.5, 5.0)]
    r = min((minimize(f, s, method='Nelder-Mead', options={'maxiter': 4000, 'xatol': 1e-7, 'fatol': 1e-10}) for s in starts), key=lambda r: r.fun)
    return r.fun, tuple(r.x)
@functools.lru_cache(maxsize=None)
def pairsum_cache(de_key, lt_key):
    import lemmaA
    return lemmaA.rho_box(de_key, lt_key)
_TAB = None
def _load():
    global _TAB
    import pickle, glob
    _TAB = {}
    for f in glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lemmaD_table*.pkl')):
        try: _TAB.update(pickle.load(open(f, 'rb')))
        except Exception: pass
def best_cached(k, lk):
    if _TAB is None: _load()
    key = (k, lk, USE_E, 13 if USE_G else (6 if USE_F else 4))
    if key not in _TAB: _TAB[key] = best(k, lk)
    return _TAB[key]
def logN(s, T):
    """log of the Lemma D bound for N(s, T), or None where it does not apply."""
    de = 1 - s; lT = math.log(T)
    if de > 0.030 or lT <= 28 or lT > 40: return None   # Oct 13: (28,29] uses the lk = 29 cell (N(s,T) <= N(s,e^29))
    k = math.ceil(de*1e4); lk = math.ceil(lT)
    if lT <= 31 and not os.environ.get('NO_HALF'): lk = math.ceil(2*lT)/2     # Oct 14: half-unit height cells below e^31
    b = best_cached(k, lk)
    if b is None: return None
    val, p = b
    # value used exponent 2 x de L at the bin's de and top L; the bound N(s,T) <= (...) T^{2x(1-s)}/((1-s) log T):
    # recompute with the actual s (<= bin de) and actual T (<= e^lk), minus log(de log T)
    x = p[3]
    # REFEREE FIX (Oct 9): the constants are evaluated at L = lk and are NOT monotone in L, so the bound is only
    # certified at T = e^lk; use N(s,T) <= N(s, e^lk) (no rescaling in T).  The rescaling in delta is kept: the
    # constants were computed at the bin's smallest sigma (worst alpha), the exponent and 1/delta use the actual delta.
    return val - 2*x*(k/1e4)*lk + 2*x*de*lk - math.log(de*lk)
if __name__ == "__main__":
    import lemmaA
    for s in [0.986, 0.988, 0.990, 0.992, 0.993, 0.994, 0.9945]:
        for T in [1e14, 5e14]:
            o = math.exp(lemmaA.logN(s, T)); n = logN(s, T)
            print(s, '%.0e' % T, 'current %.3e' % o, ' LemmaD %s' % ('%.3e' % math.exp(n) if n is not None else '-'),
                  ' params', [round(a, 3) for a in best(math.ceil((1-s)*1e4), math.ceil(math.log(T)))[1]] if n is not None else '')

def build_table(lks, ks):
    _load()
    for lk in lks:
        for k in ks: best_cached(k, lk)
    return _TAB

# ---------------------------------------------------------------------------------------------------------------
# Weighted form (Lemma W).  For the representatives r of the zeros with delta_r = 1 - beta_r <= dmax and gamma_r <= T,
# Halasz with weights a_r = 1/D_r, D_r = X^{2 delta_r}/(delta_r log T), and the Schur test (pair kernel normalised by
# sqrt(D_r D_s) is <= min(1, 2 dmax/|Delta|), row sums <= the Lemma A''/C coefficient at dmax) give
#     sum_r 1/D_r <= A := (d4 d5/w) coef / (c1^2 - c5),
# hence sum over all such zeros of x^{-delta} <= n~ A max_{nu <= delta <= dmax} x^{-delta} X^{2 delta}/(delta log T).
@functools.lru_cache(maxsize=None)
def weighted_logA(dmax_key, lk, Lx_key):
    de = dmax_key/1e4; L = float(lk); al = 1 - de; nu = 1/(R0*L); Lx = float(Lx_key)
    if de <= nu: return None
    coef, ntil = pair_coef(de, L)
    coef = coef + 2*math.sqrt(2*math.pi)*2.4*lemmaC.Iprime(dmax_key, lk)
    grid = [nu + (de - nu)*i/40 for i in range(41)]
    def f(p):
        r = consts(p, al, L)
        if r is None: return 1e9
        dd, den = r
        return math.log(dd*coef*ntil/den) + max(2*p[3]*d*L - d*Lx - math.log(d*L) for d in grid)
    starts = [(u, t, w, x) for u in (1.6, 2.0) for t in (2.0, 2.6) for w in (0.35,) for x in (5.5, 6.8)] + [(2.3, 1.8, 0.65, 5.8), (2.0, 1.75, 0.5, 5.0)]
    r = min((minimize(f, s, method='Nelder-Mead', options={'maxiter': 3000, 'xatol': 1e-6, 'fatol': 1e-9}) for s in starts), key=lambda r: r.fun)
    return r.fun     # log of  n~ A max_delta x^{-delta} X^{2 delta}/(delta log T), with x = e^{Lx}, T = e^{lk}
