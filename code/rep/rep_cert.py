# Rigorous (ball arithmetic, python-flint) certification of Lemma Rep cells:
#   for all log t in [Lt_lo, Lt_hi], eta0, eta1 in [emin, de], |d| <= dcert: one of the two orderings gives
#   margin = rhs1 + rhs2 - lhs > 0  (contradiction).  Boxes are bisected adaptively.
import math, json, sys
from flint import arb, acb, ctx
ctx.prec = 80
import bty2
th = arb(11338)/10000; kk = th.tan(); ct = 1/kk; s2 = 1/(th.cos()**2); D = 2*th*ct
W0 = (s2*(th*kk + 3*th*ct - 3))
km = [arb(1), arb(-851)/859, arb(780)/859, arb(-525)/859, arb(171)/859, arb(28)/859, arb(-29)/859]
kap = arb(433)/859
cc = [4,-8,2,20,-9,-34,27,91,-27,-201,32,895,1949,2389,1896,949,239]
ak = [arb(sum(cc[j]*cc[j+k] for j in range(17-k))*(1 if k == 0 else 2))/14912370 for k in range(17)]
a0, a1 = ak[0], ak[1]; a = sum(ak[1:], arb(0))
def J0(s): return (1 - (-s*D).exp())/s
def J1(s): return (1 - (-s*D).exp()*(1 + s*D))/(s*s)
def W(z):
    ik = acb(0, kk); al, be = th*ct, arb(-1)/2
    cosp = (al*J0(z-ik) + be*J1(z-ik) + al*J0(z+ik) + be*J1(z+ik))/2
    lin = (2*th*ct - 2)*J0(z) - J1(z)
    def sp(phi): return ((acb(0, phi)).exp()*J0(z+ik) - (acb(0, -phi)).exp()*J0(z-ik))/acb(0, 2)
    return s2*(s2*cosp + lin + sp(2*th)/(2*th).sin() - 2*sp(th)/th.sin())
# Second-order enclosure:  W(z) in W(zc) + W'(zc)(z - zc) + B(0, r^2/2 * S2),  r = |z - zc|,
#   W'(z) = -int u e^{-zu} w = closed form with J0 -> -J1, J1 -> -J2;
#   |W''(zeta)| = |int u^2 w e^{-zeta u}| <= min( P(nu),  (2 w0 + Q(nu)) / |zeta|^3 ),  nu <= Re zeta, where
#   P(nu) = int u^2 w e^{-nu u},  Q(nu) = int |h3| e^{-nu u},  h = u^2 w, h3 its third derivative
#   (h(0)=h'(0)=0, h''(0)=2w0, h(D)=h'(D)=h''(D)=0; integrate by parts three times).
#   P, Q are bounded rigorously by sum over 256 pieces of (sup over the piece, by ball evaluation) * width.
from flint import arb_series
def J2(s): return (2 - (-s*D).exp()*(s*s*D*D + 2*s*D + 2))/(s*s*s)
def Wd(z):
    ik = acb(0, kk); al, be = th*ct, arb(-1)/2
    cosp = (-al*J1(z-ik) - be*J2(z-ik) - al*J1(z+ik) - be*J2(z+ik))/2
    lin = -(2*th*ct - 2)*J1(z) + J2(z)
    def sp(phi): return -((acb(0, phi)).exp()*J1(z+ik) - (acb(0, -phi)).exp()*J1(z-ik))/acb(0, 2)
    return s2*(s2*cosp + lin + sp(2*th)/(2*th).sin() - 2*sp(th)/th.sin())
def _wser(x, n):
    ctx.cap = n + 1
    u = arb_series([x, 1])
    return s2*(s2*(th*ct - u/2)*(u*kk).cos() + 2*th*ct - u + (2*th - u*kk).sin()/(2*th).sin() - 2 - 2*(th - u*kk).sin()/th.sin())
_PQ = {}
def PQ(nu):
    k = math.floor(nu*2)/2           # AUDIT FIX (Oct 15): no lower clamp (P, Q increase as nu decreases)
    if k > 60: k = 60.0
    if k not in _PQ:
        nu_ = arb(k); n = 256; P = arb(0); Q = arb(0); Dm = float(D.upper())
        for i in range(n):
            a_ = Dm*i/n; b_ = Dm*(i+1)/n; x = arb((a_+b_)/2, (b_-a_)/2)
            c = _wser(x, 4).coeffs()
            w_, w1, w2, w3 = c[0], c[1], 2*c[2], 6*c[3]
            h3 = 6*w1 + 6*x*w2 + x*x*w3
            e = (-nu_*x).exp()
            P += arb((x*x*abs(w_)*e).upper())*(b_ - a_); Q += arb((abs(h3)*e).upper())*(b_ - a_)
        _PQ[k] = (float(P.upper()), float(Q.upper()))
    return _PQ[k]
W0up = float(W0.upper())
def Wb(z):
    zc = acb(z.real.mid(), z.imag.mid())
    r = math.hypot(float(z.real.rad().upper()), float(z.imag.rad().upper()))*(1 + 1e-12) + 1e-300
    nu = float(zc.real.mid()) - r; P, Q = PQ(nu)
    za = abs(complex(float(zc.real.mid()), float(zc.imag.mid()))) - r
    S2 = P if za <= 0 else min(P, (2*W0up + Q)/za**3)
    rr = arb(0, 0.5*r*r*S2*(1 + 1e-9))
    return W(zc) + Wd(zc)*(z - zc) + acb(rr, rr)
def F(z, eta, sig): return sum((k*Wb((z + (2*sig-1)*m)/eta) for m, k in enumerate(km)), acb(0))
lamt = {float(k): v for k, v in json.load(open('lamstar2c.json' if __import__('os').path.exists('lamstar2c.json') else 'lamstar2.json')).items()}
def lam_at(l):
    ups = [k for k in lamt if k >= l - 1e-12]; return lamt[min(ups)] if ups else 1/4.896
def ball(lo, hi): lo, hi = arb(lo), arb(hi); return (lo + hi)/2 + arb(0, ((hi - lo)/2).upper())
def consts(eta0, sig0):
    E2, E3 = bty2.errconsts(eta0, sig0); c10 = bty2.lemma10(eta0, sig0)[0]; G = bty2.gamma_const(sig0)
    r = 1 + 1e-6      # outward rounding of the (mpmath, 20-digit) constants
    return arb(E2*r), arb(E3*r), arb(c10/r), arb(G/r)
def margin(Lt, sig, e0, e1, d, C):
    E2, E3, c10, G = C
    lhs = e0*W0*((a*kap/2)*Lt - G) - c10*e0*W0 + a0*F(acb(sig-1), e0, sig).real + arb('1e-10')*e0 + arb('2e-7') \
          + E2*(e0*Lt)*e0*e0 + E3*e0**3
    r1 = a1*(F(acb(sig-1+e0), e0, sig).real + F(acb(sig-e0), e0, sig).real)
    r2 = a1*(F(acb(sig-1+e1, -d), e0, sig).real + F(acb(sig-e1, -d), e0, sig).real)
    return r1 + r2 - lhs
K, T0 = 16, 1e10
def certify(lk, Lam, dcert, maxdepth=14):
    de = Lam/lk; emin = lam_at(lk)/lk; Lt_lo = max(lk - 1, 28.7); Lt_hi = float(lk)
    # sigma uses A = lambda* at the box's highest log(K t + T0) (smaller A: weaker but valid zero-free input)
    nbox = [0]
    eta0_c = math.ceil(de*1e4)/1e4
    def sig_ball(l0, l1):
        A = lam_at(math.log(K*math.exp(l1) + T0))
        lo = 1 - A/math.log(K*math.exp(l0) + T0); hi = 1 - A/math.log(K*math.exp(l1) + T0)
        return ball(lo, hi), lo
    Cc = {}
    def rec(box, depth):
        nbox[0] += 1
        (l0, l1), (x0, x1), (y0, y1), (d0, d1) = box
        sg, siglo = sig_ball(l0, l1)
        key = math.floor(siglo*1e5)
        if key not in Cc: Cc[key] = consts(eta0_c, key/1e5)
        Lt = ball(l0, l1 + 1e-12); e0 = ball(x0, x1); e1 = ball(y0, y1); d = ball(d0, d1)
        if margin(Lt, sg, e0, e1, d, Cc[key]) > 0 or margin(Lt, sg, e1, e0, d, Cc[key]) > 0: return True
        if depth >= maxdepth: return False
        widths = [(l1-l0)/1.0, (x1-x0)/max(de-emin,1e-12)*1.0, (y1-y0)/max(de-emin,1e-12), (d1-d0)/max(dcert,1e-12)]
        i = max(range(4), key=lambda j: widths[j]); b = list(box); lo, hi = b[i]; mid = (lo+hi)/2
        for half in ((lo, mid), (mid, hi)):
            b2 = list(b); b2[i] = half
            if not rec(tuple(b2), depth+1): return False
        return True
    ok = rec(((Lt_lo, Lt_hi), (emin, de), (emin, de), (0.0, dcert)), 0)
    return ok, nbox[0]
if __name__ == "__main__":
    import time
    t = time.time()
    print(certify(29, 0.29, 0.055), time.time() - t)
