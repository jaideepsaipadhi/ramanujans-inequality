# Re-derivation of the eta0-dependent constants of Bellotti-Trudgian-Yang (arXiv:2603.21490), Lemmas 4,5,8,9,10,11,
# for a general upper bound eta0 on eta = 1 - beta0 and a general bound Lam >= eta log t (BTY use Lam = A0).
# All quantities are computed so that they are UPPER bounds where BTY have upper bounds (outward rounding by +1e-6 rel).
import math, numpy as np
from scipy.integrate import quad
th = 1.1338; tt = math.tan(th); ct = 1/tt; sec2 = 1/math.cos(th)**2
c0_ = 1/math.sin(th)*sec2; c1_ = (th - math.sin(th)*math.cos(th))*tt**4
c2_ = tt**3*math.sin(th)**2; c3_ = (th - math.sin(th)*math.cos(th))*tt**2
w0 = sec2*(th*tt + 3*th*ct - 3)
km = [1, -851/859, 780/859, -525/859, 171/859, 28/859, -29/859]
skap = sum(abs(k) for k in km); kap = sum(km)
a0 = 1.0; a1 = 865534/497079; a = 2919857/828465
def C(nu, r):
    assert r > tt
    return c0_/r*(c2_*(r+1)**2*(math.exp(-2*nu*th*ct) + 1) + c1_*r + c3_*r**3)/(r*r - tt*tt)**2
def U(y):
    y = abs(y)
    if y < 0.5: return 0.5*math.log(16/(1+4*y*y)) + 2/(1+4*y*y) + 2
    return math.log(y/2) - 2/(1+4*y*y) + 2/(3*y) + 1/(8*y*y)
def consts(eta0, Lam, sig0):
    x0 = sig0 - 0.5
    C1 = C(x0/eta0, x0/eta0); C2 = C(sig0/eta0, sig0/eta0)
    I0 = 2*quad(lambda y: U(y)/(x0*x0+y*y)**1.5, 0, 0.5, limit=200)[0] + 2*quad(lambda y: U(y)/(x0*x0+y*y)**1.5, 0.5, np.inf, limit=400)[0]
    # Lemma 8: E1 <= T1 + T2
    e8_t0 = (C1/(2*math.pi)*I0 + C2/x0**3)                  # * eta^3      (t = 0)
    e8_h2 = C1/(2*math.pi)*2/x0**2*Lam                      # * eta^2      (t >= H; uses eta log t <= Lam)
    e8_h3 = C2/x0**3                                        # * eta^3
    # Lemma 9: sum |kappa_m| E2
    e9_t0 = skap*e8_t0; e9_h2 = skap*e8_h2; e9_h3 = skap*e8_h3
    # (48): a0 E_0 + sum_{k>=1} a_k E_k
    E2 = a*e9_h2; E3 = a0*e9_t0 + a*e9_h3
    return dict(x0=x0, C1=C1, C2=C2, I0=I0, E2=E2*(1+1e-6), E3=E3*(1+1e-6))
if __name__ == "__main__":
    # reproduce BTY at eta0 = 0.0071093, A0 = 1/4.8596, sigma0 = 0.9935164: expect C1 < 52, (32) < 24, 424, 187 eta^2 + 8013 eta^3
    d = consts(0.0071093, 1/4.8596, 0.9935164)
    print({k: round(v, 3) for k, v in d.items()})
