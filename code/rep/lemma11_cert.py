# Rigorous Lemma 11 (ball arithmetic): g(u) = w(eta0 u) sum_m kappa_m e^{-d_m u} - kappa w0 >= 0 on [0, u0],
# d_m = m (kappa_m > 0), (2 sigma0 - 1) m (kappa_m < 0).  g(0) = 0 and g'(0) = -w0 sum kappa_m d_m > 0; near 0 we use
# g(u) >= u (g'(0) - (u/2) sup|g''|) with sup|g''| over [0, u1] by ball evaluation of the analytic second derivative.
from flint import arb, ctx
ctx.prec = 80
th = arb(11338)/10000; kk = th.tan(); ct = 1/kk; s2 = 1/(th.cos()**2)
km = [arb(1), arb(-851)/859, arb(780)/859, arb(-525)/859, arb(171)/859, arb(28)/859, arb(-29)/859]
kap = sum(km, arb(0)); w0 = s2*(th*kk + 3*th*ct - 3)
def wfun(u):   # and derivatives (closed form) for u in [0, D]
    c = (u*kk).cos(); s = (u*kk).sin()
    A = th*ct - u/2
    w = s2*(s2*A*c + 2*th*ct - u + (2*th - u*kk).sin()/(2*th).sin() - 2 - 2*(th - u*kk).sin()/th.sin())
    w1 = s2*(s2*(-c/2 - A*kk*s) - 1 - kk*(2*th - u*kk).cos()/(2*th).sin() + 2*kk*(th - u*kk).cos()/th.sin())
    w2 = s2*(s2*(kk*s - A*kk*kk*c) - kk*kk*(2*th - u*kk).sin()/(2*th).sin() + 2*kk*kk*(th - u*kk).sin()/th.sin())
    return w, w1, w2
def g_all(u, eta0, sig0):
    dm = [arb(m) if km[m] > 0 else (2*arb(sig0) - 1)*m for m in range(7)]
    S = sum((km[m]*(-dm[m]*u).exp() for m in range(7)), arb(0))
    S1 = sum((-dm[m]*km[m]*(-dm[m]*u).exp() for m in range(7)), arb(0))
    S2 = sum((dm[m]**2*km[m]*(-dm[m]*u).exp() for m in range(7)), arb(0))
    w, w1, w2 = wfun(eta0*u)
    g = w*S - kap*w0
    g2 = eta0*eta0*w2*S + 2*eta0*w1*S1 + w*S2
    g1_0 = -w0*sum((km[m]*dm[m] for m in range(7)), arb(0))
    return g, g2, g1_0
def certify(eta0, sig0, u0, u1=0.05):
    eta0 = arb(eta0)
    U = arb(0) + arb(u1/2, u1/2)               # [0, u1]
    _, g2, g1_0 = g_all(U, eta0, sig0)
    if not (g1_0 - (arb(u1)/2)*abs(g2) > 0): return False
    stack = [(u1, u0)]
    while stack:
        a_, b_ = stack.pop()
        uu = arb((a_+b_)/2, (b_-a_)/2)
        if g_all(uu, eta0, sig0)[0] > 0: continue
        if b_ - a_ < 1e-7: return False
        m = (a_ + b_)/2; stack += [(a_, m), (m, b_)]
    return True
if __name__ == "__main__":
    import bty2
    for eta0 in (0.0071093, 0.0085, 0.0125):
        s0 = 0.9932; u0 = bty2.lemma11_u0(eta0, s0)
        print(eta0, 'u0 (grid)', u0, 'certified on [0,u0-0.01]:', certify(eta0, s0, u0 - 0.01))
