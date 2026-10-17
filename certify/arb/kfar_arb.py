"""Ball-arithmetic bound for the parts of the k1, k0 integrals (Theorem G) not covered by the rigorous quadrature of
code/lemG/kcert2.py: the segment 150 <= |Im s| <= T2 = 400 on Re s = kap = 0.38 and the two horizontal pieces of the path P
at Im s = +-T2, Re s in [-c', kap].  At w = 0 the integrand of H_c is
    g(s) = A(s) e^{eta^2 s^2} / (zeta(1-s) zeta(1+s) s^4),     A(s) = A(0; -s, s) = prod_p A_p,
and for H_c'(0) it is multiplied by D(s) = dlogA(s) + eta^2 s - zeta'/zeta(1+s) - 2/s.
We bound sup|g| and sup|g D| on unit pieces by evaluating everything on complex balls; the Euler product over p <= 2000 is
evaluated on the ball and the tail p > 2000 is bounded by |A/A_P - 1| <= e^tau - 1, tau = 1.1 * 4 * P^{1-m}/(m-1),
|dlogA - dlogA_P| <= 2.2 (log P/(m-1) + 1/(m-1)^2) P^{1-m}, with m = 1.62 the least pairwise sum of the real parts.
Output: upper bounds F0 >= (1/2pi) int |g|, F1 >= (1/2pi) int |g D| over these parts; compare with the 'far' term of kcert2."""
from flint import arb, acb, acb_series, ctx
import math
ctx.prec = 64
eta, T2, cp, kap = arb('0.02'), 400, arb('0.01'), arb('0.38')
P = 2000
primes = [p for p in range(2, P + 1) if all(p % d for d in range(2, int(p**0.5) + 1))]
logs = {p: arb(p).log() for p in primes}
m = arb('1.62'); Pa = arb(P)
tau = arb('1.1')*4*Pa**(1 - m)/(m - 1)
dtail = arb('2.2')*(Pa.log()/(m - 1) + 1/(m - 1)**2)*Pa**(1 - m)
def g_and_D(s):
    A = acb(1); dA = acb(0)
    for p in primes:
        lp = logs[p]
        x1 = (-(1 - s)*lp).exp(); x2 = (-(1 + s)*lp).exp(); z = acb(1)/p
        N1 = 1 - x1 - x2 + z
        A *= N1*(1 - z)/((1 - x1)*(1 - x2))
        dA += -lp*x2*(1/(1 - x2) - 1/N1)
    Aabs = arb(A.abs_upper())*(tau.exp())                 # |A| <= |A_P| e^tau
    zs = acb_series([1 + s, 1], prec=2).zeta(); co = zs.coeffs(); z0, z1 = co[0], co[1]
    E = (eta**2*s**2).exp()
    MI = arb('14.57')                                    # (B1): |1/zeta| <= 14.57 on 0.6 <= sigma <= 1.3, |t| <= 450 (iota(1.3) < 14.57 above)
    inv1 = MI                                            # 1/zeta(1-s): Re(1-s) in [0.62, 1.01], |Im| <= 400
    z0l = arb(z0.abs_lower())
    inv2 = 1/z0l if z0l > 0 else MI
    g = Aabs*arb(E.abs_upper())*inv1*inv2/arb((s**4).abs_lower())
    D = arb(dA.abs_upper()) + dtail + arb((eta**2*s).abs_upper()) + arb((z1/z0).abs_upper()) + arb((2/s).abs_upper())
    return g, g*D
F0 = arb(0); F1 = arb(0)
# vertical pieces 150 <= |t| <= 400 on Re s = kap (both signs; |g| is symmetric under conjugation)
h = 0.02; N = int(round((T2 - 150)/h))
for i in range(N):
    t = 150 + i*h
    s = acb(kap, arb(t + h/2, h/2))
    g, gD = g_and_D(s)
    F0 += 2*g*h; F1 += 2*gD*h
# horizontal pieces Im s = +-T2, Re s in [-c', kap], in pieces of width 0.05
n = 40
for j in range(n):
    lo = -cp + (kap + cp)*j/n; hi = -cp + (kap + cp)*(j + 1)/n
    s = acb(arb((lo + hi)/2, ((hi - lo)/2).upper()), T2)
    g, gD = g_and_D(s)
    F0 += 2*g*(hi - lo); F1 += 2*gD*(hi - lo)
F0 /= 2*arb.pi(); F1 /= 2*arb.pi()
print('far part of k1 integral <=', F0.upper())
print('far part of H_c\'(0) integral <=', F1.upper())
far_kcert2 = 2*30*14.57*3.4*math.exp(0.02**2*0.38**2 - 0.02**2*150*150)/(3*150**3)/(2*math.pi)
print('kcert2 used: far =', far_kcert2, ', far*(10 + eta^2 T2) =', far_kcert2*(10 + 0.02**2*400))
print('OK' if F0 < far_kcert2 and F1 < far_kcert2*(10 + 0.02**2*400) else 'CHECK')
