# Self-contained re-derivation of every eta0/sigma0-dependent constant in BTY's Lemma 15 inequality.
#   W0 bound (replaces Ford's Lemma 4): integrating by parts four times (w(d)=w'(d)=w''(d)=0, w'(0)=0),
#     W0(z) = w''(0)/z^3 + (w'''(0) - e^{-zd} w'''(d))/z^4 + z^{-4} int e^{-zu} w''''(u) du,
#     so |W0(z)| <= Cw(nu, r)|z|^{-3} for Re z >= nu, |z| >= r, with
#     Cw(nu, r) = |w''(0)| + (|w'''(0)| + e^{-nu d}|w'''(d)| + int_0^d |w''''(u)| e^{-nu u} du)/r.
import math, mpmath as mp
mp.mp.dps = 20
th = mp.mpf(11338)/10000; tt = mp.tan(th); ct = 1/tt; s2 = 1/mp.cos(th)**2; D = 2*th*ct
def w(u): return s2*(s2*(th*ct-u/2)*mp.cos(u*tt)+2*th*ct-u+mp.sin(2*th-u*tt)/mp.sin(2*th)-2*(1+mp.sin(th-u*tt)/mp.sin(th)))
def dw(u, k): return mp.diff(w, u, k)
w0 = w(0)
# AUDIT FIX: exact values from arb_series (misc_cert.py), rounded UP
W2_0 = mp.mpf('50.06970149467168733130'); W3_0 = mp.mpf('117.2647765286022616124'); W3_d = W3_0
_w4 = {}
def I4(nu):
    if nu not in _w4: _w4[nu] = mp.quad(lambda u: abs(dw(u, 4))*mp.e**(-nu*u), mp.linspace(0, D, 9))
    return _w4[nu]
import json as _json, os as _os
_I4C = {int(k): v for k, v in _json.load(open(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'I4cert.json'))).items()}
def I4c(nu):   # certified (consts_cert.py, flint) upper bound; I4 decreasing in nu -> integer grid point below
    k = math.floor(nu); assert 30 <= k <= 200, nu
    return _I4C[k]
def Cw(nu, r): return (float(W2_0) + (float(W3_0) + float(mp.e**(-nu*D))*float(W3_d) + I4c(nu))/r)*(1 + 1e-12)
km = [1, -851/859, 780/859, -525/859, 171/859, 28/859, -29/859]
kap = sum(km); skap = sum(abs(k) for k in km)
ak = [1] + [0]*16
# P(x) = |sum c_k e^{ikx}|^2 / 14912370 ; a_k = 2*sum c_j c_{j+k}/14912370 (k>=1), a_0 = sum c_j^2/14912370
cc = [4,-8,2,20,-9,-34,27,91,-27,-201,32,895,1949,2389,1896,949,239]
ak = [sum(cc[j]*cc[j+k] for j in range(17-k))*(1 if k == 0 else 2)/14912370 for k in range(17)]
a0, a1, a = ak[0], ak[1], sum(ak[1:])
_I0 = {}
def I0(x0):     # int |Re psi(1/4 + iy/2)| (x0^2+y^2)^{-3/2} dy  (the true digamma; BTY (32))
    if x0 not in _I0:
        f = lambda y: abs(mp.re(mp.digamma(mp.mpf(1)/4+1j*y/2)))/(x0**2+y**2)**1.5
        _I0[x0] = 2*mp.quad(f, [0, 0.25, 0.5, 1, 2, 5, 20, 100, mp.inf])
    return float(_I0[x0])
def errconsts(eta0, sig0):
    """(E2, E3): sum_k a_k |E_k| <= E2*lam*eta^2 + E3*eta^3, where lam = eta log t (primary zero), eta <= eta0, sigma >= sig0.
    Lemma 8 for t >= H: T1 <= Cw/(2pi) eta^3 * (2/x0^2) log t (+1e-10 pieces) = Cw/(pi x0^2) lam eta^2."""
    x0 = sig0 - 0.5
    C1 = Cw(x0/eta0, x0/eta0); C2 = Cw(sig0/eta0, sig0/eta0)
    assert x0 >= 0.4930
    t0 = C1/(2*math.pi)*23.6497 + C2/x0**3      # I0(0.4930) <= 23.64970 (consts_cert.py, flint); I0 decreasing in x0
    h2 = C1/(math.pi*x0**2)*1.000001; h3 = C2/x0**3
    return a*skap*h2, (a0*skap*t0 + a*skap*h3)*1.000001
def gamma_const(sig0):
    """Lemma 12: sum_k a_k sum_m kappa_m(-1/2 log pi + 1/2 Re psi((s_k+delta_m)/2+1)) <= (a kappa/2) log t - G."""
    t1 = a*kap/2*(-math.log(2*math.pi) + 1e-10) + kap/2*sum(ak[k]*math.log(k) for k in range(1, 17))
    t0 = 0.0
    for m, k in enumerate(km):
        lo = float(mp.digamma((sig0 + (2*sig0-1)*m)/2 + 1)); hi = float(mp.digamma((m + 3)/mp.mpf(2)))
        t0 += 0.5*k*(hi if k > 0 else lo)
    t0 += a0*(-0.5*math.log(math.pi))*kap        # the -1/2 log pi of the k = 0 term
    return -(t1 + t0)
# Lemma 10/11
mp_tab = {2:0.23416332,3:0.05006283,5:0.0073564,   # AUDIT (Oct 15): BTY Table 1 says 0.00735709 but min G_5(1,x) = 0.00735649 at x = 1.64952; certified 0.0073564
          7:0.00205251,11:0.00039742,13:0.00022273,17:0.00009071,19:0.00006314,
          23:0.00003427,29:0.00001656,31:0.00001346,37:7.8e-6,41:5.69e-6,43:4.91e-6,47:3.74e-6,53:2.59e-6,59:1.87e-6,
          61:1.69e-6,67:1.27e-6,71:1.06e-6,73:9.7e-7,79:7.7e-7,83:6.6e-7,89:5.3e-7,97:4.1e-7}
def lemma11_u0(eta0, sig0):
    """largest u0 (grid 0.05) with w(eta0 u) sum kappa_m e^{-d_m u} >= kappa w0 on [0, u0]."""
    import numpy as np
    th_ = float(th); tt_ = float(tt); ct_ = float(ct); s2_ = float(s2); U = 2*th_*ct_
    def wf(u): return np.where(u <= U, s2_*(s2_*(th_*ct_-u/2)*np.cos(u*tt_)+2*th_*ct_-u+np.sin(2*th_-u*tt_)/math.sin(2*th_)-2*(1+np.sin(th_-u*tt_)/math.sin(th_))), 0.0)
    u = np.linspace(0, 80, 160001)
    s = sum(k*np.exp(-((2*sig0-1)*m if k < 0 else m)*u) for m, k in enumerate(km))
    g = wf(eta0*u)*s - kap*float(w0)
    bad = np.nonzero(g < -1e-9)[0]
    return 80.0 if bad.size == 0 else float(u[bad[0]] - 0.05)
def lemma10(eta0, sig0):
    """constant c10 with sum Lambda(n) n^-sigma f(log n) P(t log n) >= c10 * eta w0 (f(0) = kappa eta w0)."""
    import lemma11_cert
    u0 = lemma11_u0(eta0, sig0) - 0.01; assert lemma11_cert.certify(eta0, sig0, u0); s = 0.0
    for p, m_p in mp_tab.items():
        if 15*math.log(p) > u0: continue
        s += math.log(p)*m_p       # AUDIT FIX: dropped terms (p^m > e^u0, m > 15) are >= 0 since f >= 0 (BTY Lemma 3, f_nonneg_cert.py)
    return kap*kap*s, u0, s      # kappa f(0) * s = kappa^2 eta w0 s
if __name__ == "__main__":
    print("a0 a1 a", a0, a1, a, " (BTY: 1, %.6f, %.6f)" % (865534/497079, 2919857/828465))
    print("w0 =", w0, " |w''(0)| =", W2_0, " Cw(69.4) =", Cw(69.4, 69.4), " Cw(152,152) =", Cw(152.5, 152.5), " Cw(-1,1e10)", Cw(-1, 1e10))
    print("BTY reproduce (eta0=.0071093, sig0=.9935164): errconsts", errconsts(0.0071093, 0.9935164), "* (lam=A0=.2058 -> E2*lam = %.1f)" % (errconsts(0.0071093, 0.9935164)[0]/4.8596))
    print("gamma const", gamma_const(0.9935164), "(BTY 1.568)")
    print("lemma10 BTY", lemma10(0.0071093, 0.9935164))
    for eta0 in (0.0085, 0.010, 0.0115, 0.013):
        sig0 = 1 - 0.2130/31.50
        print(eta0, "u0", lemma11_u0(eta0, sig0), "c10", lemma10(eta0, sig0)[0], "err", errconsts(eta0, sig0), "G", gamma_const(sig0))
