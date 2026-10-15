"""Logan-kernel (Buthe) version of the error budget.
psi vs psi_{c,eps}:   smoothing cost  = 1.601*eps/sqrt(c) * L/(L-log(sqrt(c)/eps))   [Buthe E4, Brun-Titchmarsh factor]
zeros |gamma|<=c/eps: weighted by a(rho) <= (c/sinh c) sinh(sqrt(c^2-u^2))/sqrt(c^2-u^2), u=gamma*eps
zeros |gamma|> c/eps: 0.16*2/sinh(c) * e^{0.71 sqrt(c eps)} log(3c) log(c/eps)   (Buthe Prop 3, relative)
Off-line zeros above H via Stieltjes with f(t)=w(t)/t decreasing, N at right endpoints (upper bound)."""
import math, ramanujan_upper as r, bell_patch
def w(u, c):
    if u >= c: return c/math.sinh(c)
    s = math.sqrt(c*c-u*u)
    return math.exp(s - c + math.log(c/s) + math.log((1-math.exp(-2*s))/(1-math.exp(-2*c))))
def zsum(s, H, eps, c, nt=80):
    T1 = c/eps
    if T1 <= H: return 0.0
    us = [math.log(H)+(math.log(T1)-math.log(H))*i/nt for i in range(nt+1)]
    f = [w(math.exp(u)*eps, c)/math.exp(u) for u in us]
    acc = 0.0
    for i in range(nt):
        t2 = math.exp(us[i+1])
        if s >= 1-1/(r.Rof(us[i+1])*us[i+1]): continue  # rigorous: zero-free edge at the TOP of the block
        acc += math.exp(r.logN_density(s, t2))*(f[i]-f[i+1])
    if s < 1-1/(r.Rof(us[-1])*us[-1]): acc += math.exp(r.logN_density(s, T1))*f[-1]
    return acc
def err(L, H, tol, c, frac, ns=200):
    # choose eps so smoothing cost = frac*tol
    # exact first moment of the Logan transform: E|t| = 0.794 eps/sqrt(c) (c in [30,70]);
    # Brun-Titchmarsh on the one-sided deviation gives factor (L+lam)/(L-lam), lam=log(10 sqrt(c)/eps)
    eps = frac*tol*math.sqrt(c)/0.7936
    for _ in range(3):
        lam = math.log(10*math.sqrt(c)/eps); kap = (L+lam)/(L-lam)
        eps = frac*tol*math.sqrt(c)/(0.7936*kap)
    sm = frac*tol
    tr = 0.32/math.sinh(c)*math.exp(0.71*math.sqrt(c*eps))*math.log(3*c)*math.log(c/eps)
    T1 = c/eps
    tail = 0.0
    if T1 > H:
        smax = 1-1/(r.Rof(math.log(T1))*math.log(T1)); h = (smax-0.6)/ns; acc = 0.0
        for i in range(ns):
            s = 0.6+(i+0.5)*h
            acc += math.exp((s-1)*L)*zsum(s, H, eps, c)
        tail = 2*(L*acc*h + math.exp(-0.4*L)*(zsum(0.6,H,eps,c)+math.log(T1)**2))
    return sm+tr+tail, T1
def total_err(L, H, tol, bell_scale=1.0):
    best=(float('inf'),None)
    for c in [30,40,50,60,70]:
        for frac in [0.5,0.7,0.85]:
            e=err(L,H,tol,c,frac)
            if e[0]<best[0]: best=e
    return best
r.total_err = total_err
