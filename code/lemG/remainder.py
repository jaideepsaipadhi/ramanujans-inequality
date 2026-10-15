# Lemma G, step 3: explicit bounds for every remainder in the smoothed contour argument (see LEMMA_G.md, "Step 3").
# Everything here is an explicit majorant; the only computed constants are the arb box bounds (box2.py):
#   sigma in [0.70, 1.30], |t| <= 450:  |zeta| <= MZ (for |s-1| >= 0.1),  |1/zeta| <= MI.
# 1-D integrals of explicit majorants: rigorous upper Riemann sums (monotone factors, cell-wise sups) or closed forms.
import mpmath as mp, math, sys, functools, os, pickle
_CF = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rem_cache.pkl')
try: _CACHE = pickle.load(open(_CF, 'rb'))
except Exception: _CACHE = {}
def _memo(f):
    def g(*a):
        k = (f.__name__,) + a
        if k not in _CACHE:
            _CACHE[k] = f(*a)
            try: pickle.dump(_CACHE, open(_CF + '.%d' % os.getpid(), 'wb')); os.replace(_CF + '.%d' % os.getpid(), _CF)
            except Exception: pass
        return _CACHE[k]
    return g
functools_lru = functools.lru_cache
mp.mp.dps = 30
MZ, MI = 11.89, 14.57      # box2.py: sigma in [0.6,1.3], |t| <= 450 (box450.log, box_060.log)
SAFE = 1.0          # all integrals are now rigorous upper sums / closed forms
P = dict(eta=0.02, T2=400.0, W0=40.0, c=0.02, cp=0.01, sig0=0.01, th=0.40, b=0.38, thp=0.30, kap=0.38)
@_memo
def zeta_r(s):  return float(mp.zeta(s))                       # real s > 1
@_memo
def izeta_r(s): return float(mp.zeta(s)/mp.zeta(2*s))           # |1/zeta(s+it)| <= zeta(s)/zeta(2s), s > 1
_pr = None
@_memo
def Abound(e1, e2, e12, P_=10**6):
    """sup |A| for Re(1+w+s1) = e1, Re(1+w+s2) = e2, Re(1+w+s1+s2) = e12 (all > 1/2, pairwise sums > 1).
    A_p - 1 = (z x1 + z x2 - z^2 - x1 x2)/((1-x1)(1-x2)),  |A| <= exp(sum_p eps_p)."""
    global _pr
    if _pr is None:
        n = P_; sv = bytearray([1])*(n+1); sv[0:2] = b'\x00\x00'
        for i in range(2, int(n**0.5)+1):
            if sv[i]: sv[i*i::i] = bytearray(len(sv[i*i::i]))
        _pr = [i for i in range(n+1) if sv[i]]
    import numpy as np
    tot = 0.0
    n = 96; h = 2*math.pi/n; ph = np.arange(n)*h
    a1, a2, aw = np.meshgrid(ph, ph, ph, indexing='ij')
    for p in _pr:
        if p > 100: break
        r1, r2, rz = p**-e1, p**-e2, p**-e12
        x1 = r1*np.exp(-1j*a1); x2 = r2*np.exp(-1j*a2); z = rz*np.exp(-1j*(a1 + a2 - aw))
        # rigorous cell bound: phases within h/2 of the grid point move x1, x2 by <= r h/2 and z by <= 3 rz h/2
        e = h/2
        N1 = np.abs(1 - x1 - x2 + z) + (r1 + r2)*e + 3*rz*e
        N2 = np.abs(1 - z) + 3*rz*e
        D1 = np.abs(1 - x1) - r1*e; D2 = np.abs(1 - x2) - r2*e
        Ap = N1*N2/(D1*D2)
        tot += math.log(float(Ap.max()))       # rigorous sup over the 3-torus (cell enclosures)
    for p in _pr:
        if p <= 100: continue
        x1, x2, z = p**-e1, p**-e2, p**-e12
        tot += (z*x1 + z*x2 + z*z + x1*x2)/((1-x1)*(1-x2))
    m = min(e12+e1, e12+e2, 2*e12, e1+e2)
    tot += 4.2*P_**(1-m)/(m-1)                     # tail: each term <= 1.05 n^{-m} for p > 1e6, four terms
    return math.exp(tot)
import functools
@_memo
def _gam(sig, W0, w1):
    # rigorous upper Riemann sum: |Gamma(sig+iv)| decreasing in v >= 0 (sig > -1); weight (1+v) increasing -> right end
    n = 40000; h = W0/n; tot = 0.0
    for i in range(n):
        g = float(abs(mp.gamma(mp.mpc(sig, i*h))))
        tot += g*((1 + (i+1)*h) if w1 else 1.0)*h
    return 2*tot*(1 + 1e-12)
def gam_int(sig, W0, wfun=None):
    return _gam(sig, W0, wfun is not None)
@_memo
def _J(sig, eta, lo, hi, w1):
    top = min(hi, 3000.0); tot = 0.0; t0 = lo
    while t0 < top:
        h = 1e-4 if t0 < 0.1 else (1e-3 if t0 < 2 else 0.01)
        t1 = min(t0 + h, top)
        tot += math.exp(eta**2*(sig**2 - t0**2)/2)*((1 + t1) if w1 else 1.0)/(sig**2 + t0**2)*(t1 - t0)
        t0 = t1
    # tail t > 3000: e^{-eta^2 t^2/2} (1+t)/t^2 <= e^{-1800}: below double precision
    return 2*tot*(1 + 1e-12) + (1e-300 if hi > 3000 else 0.0)
def hs_int(sig, lU, lV, D, eta, lo=0.0, hi=mp.inf, wfun=None):
    """int over lo <= |t| <= hi of |hat rho_s(sig+it)| wfun(t) dt (majorant (V^sig + U^sig) e^{eta^2(sig^2-t^2)/2}/(D|s|^2))."""
    return (math.exp(sig*lV) + math.exp(sig*lU))/D*_J(sig, eta, float(lo), float(hi), wfun is not None)
@_memo
def _R6(eta, T2, W0, c, cp, thp, kap):
    gauss = math.exp(-eta**2*T2**2/2)
    AK = Abound(1-kap, 1+kap-thp, 1.0); Ax = Abound(1+cp, 1+c-cp, 1.0)
    def Kc_bound(v):
        a_, b_ = kap, kap - thp          # int dt/((t^2+a^2)((t+v)^2+b^2)) = pi (a+b)/(a b ((a+b)^2 + v^2))   (exact)
        seg = math.pi*(a_ + b_)/(a_*b_*((a_ + b_)**2 + v*v))
        segc = math.exp(eta**2*((kap-thp)**2 + kap**2)/2)
        # path segment Re s = kap: |zeta(1+w)| <= MZ (Re 1-thp, |w| >= thp), 1/zeta(1-s) <= MI, 1/zeta(1+w+s) explicit (Re 1+kap-thp)
        return (MZ*AK*MI*izeta_r(1+kap-thp)*seg*segc + MZ*max(AK, Ax)*MI*MI*4*(kap+cp)*gauss/T2**4)/(2*math.pi)
    Ew = math.exp(eta**2*thp**2/2)/thp**2
    n = 40000; h = W0/n; tot = 0.0
    for i in range(n):        # |Gamma(-thp+iv)| and Kc_bound decreasing in v >= 0: left endpoints
        tot += float(abs(mp.gamma(mp.mpc(-thp, i*h))))*(Kc_bound(i*h) + Ew)*h
    return 2*tot*(1 + 1e-12), Kc_bound(W0) + math.exp(eta**2*thp**2/2)/W0**2
def bounds(u, tau, x, L, Y1min=None, verbose=False):
    eta, T2, W0, c, cp, s0, th, b, thp, kap = [P[k] for k in ('eta','T2','W0','c','cp','sig0','th','b','thp','kap')]
    lU, lV, lX = u*L, tau*u*L, x*L; D = lV - lU; Y1, Y2 = lX - lV, lX - lU
    pi3 = (2*math.pi)**3
    out = {}
    one = lambda t: 1 + abs(t)
    # ---------- E1: s1 moved to C1 (segment Re -theta, |t1|<=T2 ; connectors ; tails Re sig0) ----------
    Gw = gam_int(c, W0, one)*math.exp(c*lX)
    Zb = max(MZ, 1/(th-c+cp) + 1/(1+c-th-cp))                    # |zeta(1+w+s1+s2)| <= Zb (1+|v|)(1+|t1|)(1+|t2|) at Re 0.71
    I1 = hs_int(-th, lU, lV, D, eta, 0, T2, one); I2 = hs_int(-cp, lU, lV, D, eta, 0, mp.inf, one)
    A1 = Abound(1+c-th, 1+c-cp, 1+c-th-cp)
    out['E1a'] = SAFE*Gw*zeta_r(1+c)*I1*I2*Zb*MI*izeta_r(1+c-cp)*A1/pi3
    # connectors (|t1| = T2, Re s1 in [-th, s0]) and tails (Re s1 = s0, |t1| > T2): Gaussian factor e^{-eta^2 T2^2/2}
    gauss = math.exp(-eta**2*T2**2/2)
    # connectors: Im s1 = +-T2, Re s1 in [-th, s0]; u = w+s1+s2, Re(1+u) >= 1+c-th-cp > 0, and
    #   |zeta(1+u)| <= 1/|u| + 1 + |1+u|/Re(1+u),  |1+u| <= 2 + W0 + T2 + |t2|.
    # int over the connector x t2 of |hat rho_s(s2)| |zeta(1+u)| <= len (int |hs(s2)| (2 + (2+W0+T2+|t2|)/sm) dt2) + 2 pi R0 Hinf,
    # the 1/|u| singularity handled on a disk of radius R0 = sqrt(len^2+1) (int_disk 1/r = 2 pi R0), elsewhere 1/|u| <= 1.
    ln = th + s0; sm = 1 + c - th - cp; R0 = math.sqrt(ln*ln + 1)
    Hinf = (math.exp(-cp*lV) + math.exp(-cp*lU))*math.exp(eta**2*cp**2/2)/(D*cp**2)
    J0 = hs_int(-cp, lU, lV, D, eta, 0, mp.inf); J1 = hs_int(-cp, lU, lV, D, eta, 0, mp.inf, one)
    zconn = ln*((2 + (2 + W0 + T2)/sm)*J0 + J1/sm) + 2*math.pi*R0*Hinf
    hs1 = (math.exp(s0*lV) + math.exp(s0*lU))*math.exp(eta**2*(th**2 - T2**2)/2)/(D*T2**2)      # sup |hat rho_s| on connectors
    A1c = Abound(1+c-th, 1+c-cp, 1+c-th-cp)
    E1b = 2*gam_int(c, W0)*math.exp(c*lX)*zeta_r(1+c)*hs1*zconn*MI*izeta_r(1+c-cp)*A1c/pi3
    E1c = Gw*zeta_r(1+c)*hs_int(s0, lU, lV, D, eta, T2, mp.inf, one)*zeta_r(1+c+s0-cp)*izeta_r(1+c+s0)*I2*izeta_r(1+c-cp)*A1c/pi3
    out['E1bc'] = SAFE*(E1b + E1c)
    # ---------- R_x extension (|Im(w+s2)| >= T2) and K_tail : both carry hat rho_s / E at height >= T2 ----------
    Ax = Abound(1+cp, 1+c-cp, 1.0 + 1e-9) if False else Abound(1+cp, 1+c-cp, 1.0)
    Gx = gam_int(c, W0)*max(math.exp(c*Y1), math.exp(c*Y2))
    ext = (math.exp(-(c-cp)*lU) + math.exp(-(c-cp)*lV))*math.exp(eta**2*(c-cp)**2/2)*gauss/(D*T2**2)
    out['Rext'] = SAFE*Gx*zeta_r(1+c)*Ax*izeta_r(1+cp)*izeta_r(1+c-cp)*ext*hs_int(-cp, lU, lV, D, eta, 0, mp.inf)/(2*math.pi)**2
    Ktail = 2*math.exp(eta**2*((c-cp)**2 + cp**2)/2)*gauss*2/(T2**3*3)   # int_{|t|>T2} e^{..}/(|w+s|^2|s|^2) <= 2 gauss/(3 T2^3) (per side)
    out['Ktail'] = SAFE*2*Gx*zeta_r(1+c)*Ax*izeta_r(1+cp)*izeta_r(1+c-cp)*Ktail/(2*math.pi)**2/D**2
    # ---------- cross-term shifts: rho_V (right to Re b), rho_U (left to Re -b) ----------
    AV = Abound(1-b, 1+c+b, 1.0); AU = Abound(1+b, 1+c-b, 1.0)
    segV = math.pi/(b**2*(b+c))*math.exp(eta**2*((c+b)**2 + b**2)/2)
    segU = math.pi/(b**2*(b-c))*math.exp(eta**2*((c-b)**2 + b**2)/2)
    rV = math.exp(-D*b)*zeta_r(1+c)*AV*MI*izeta_r(1+c+b)*segV/(2*math.pi)
    rU = math.exp(-D*b)*zeta_r(1+c)*AU*MI*izeta_r(1+b)*segU/(2*math.pi)
    connV = math.exp(cp*D)*gauss*4*(b+cp)*zeta_r(1+c)*max(AV, Ax)*MI*izeta_r(1+c-cp)/(T2**4)/(2*math.pi)
    out['rhoV'] = SAFE*gam_int(c, W0)*math.exp(c*Y1)*(rV + connV)/D**2
    out['rhoU'] = SAFE*gam_int(c, W0)*math.exp(c*Y2)*(rU + connV)/D**2
    # ---------- w-shift of the main term to Re w = -thp ----------
    R6, KcW0 = _R6(eta, T2, W0, c, cp, thp, kap)
    out['wshift'] = SAFE*(math.exp(-thp*Y1) + math.exp(-thp*Y2))*R6/(2*math.pi)/D**2
    # Gamma tails |v| > W0 and w-connectors: |Gamma(sig +- i W0)| <= 3 W0^{c-1/2} e^{-pi W0/2}
    # |Im w| > W0 on the original lines: |sum Psi_s^2 n^{-1-w}| <= sum tau(n)^2 n^{-1-c} = zeta(1+c)^4/zeta(2+2c) (|rho_s| <= 1);
    # |Gamma(c+iv)| <= 2.6 |v|^{c-1/2} e^{-pi|v|/2} for |v| >= 1.  w-connectors at |v| = W0 (main-term shift): |Gamma| same bound,
    # times sup |K_c| + |E|/|w|^2 (Kc_bound at |v| = W0, valid for Re w in [-thp, c]) and e^{c Y2}.
    gt = 2*2.6*W0**(c-0.5)*math.exp(-math.pi*W0/2)*2/math.pi
    out['wtail'] = gt*math.exp(c*lX)*zeta_r(1+c)**4/zeta_r(2+2*c)/(2*math.pi) + gt*(thp+c)*2*math.exp(c*Y2)*KcW0/D**2
    tot = sum(out.values())
    if verbose:
        for k, v in out.items(): print("  %-7s %.3e" % (k, v))
    return tot, out
if __name__ == "__main__":
    u, tau, x, L = [float(a) for a in sys.argv[1:5]]
    t, o = bounds(u, tau, x, L, verbose=True)
    print("TOTAL contour remainder <= %.3e   (Y1 = %.1f)" % (t, (x - tau*u)*L))
