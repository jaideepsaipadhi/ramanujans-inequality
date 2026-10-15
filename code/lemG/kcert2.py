# Rigorous enclosures of k1, k0 for the smoothed kernel K_c(w) = zeta(1+w) H_c(w)   (Lemma G, step 3):
#   H_c(w) = (1/2 pi i) int_path A(w;-w-s,s) E(w+s) E(s) / (zeta(1-s) zeta(1+w+s) (w+s)^2 s^2) ds,  E(s) = e^{eta^2 s^2/2},
#   path -c'-iT2 -> kap-iT2 -> kap+iT2 -> -c'+iT2.   k1 = H_c(0),  k0 = H_c'(0) + gamma H_c(0).
# H_c'(0) by differentiating under the integral:  d/dw log(integrand) at w = 0 equals
#   D(s) = dlogA + eta^2 s - zeta'/zeta(1+s) - 2/s,   dlogA = sum_p (-log p) x2 (1/(1-x2) - 1/(1-x1-x2+z)),
#   with x1 = p^{-(1-s)} (w-independent since s1 = -w-s), x2 = p^{-(1+s)}, z = 1/p.
# Pieces: (i) |t| <= T on Re s = kap by acb.integral with primes <= P;  (ii) Euler tails: |A - A_P| <= |A_P|(e^tau - 1),
#   |dlogA - dlogA_P| <= tau1; int |g| <= sqrt(pi int |g|^2 (1+t^2)) computed as g(t) conj(g(conj t)) (analytic);
#   (iii) |t| in [T, T2] and the two connectors: Gaussian majorant.
from flint import acb, arb, acb_series, ctx
import math
ctx.prec = 64; ctx.cap = 3
eta, T2, cp, kap = 0.02, 400.0, 0.01, 0.38
T = 150; P = 2000
primes = [p for p in range(2, P+1) if all(p % d for d in range(2, int(p**0.5)+1))]
logs = {p: acb(p).log() for p in primes}
def parts(t):
    s = acb(kap) + acb(0, 1)*t
    A = acb(1); dA = acb(0)
    for p in primes:
        lp = logs[p]
        x1 = (-(1-s)*lp).exp(); x2 = (-(1+s)*lp).exp(); z = acb(1)/p
        N1 = 1 - x1 - x2 + z
        A *= N1*(1-z)/((1-x1)*(1-x2))
        dA += -lp*x2*(1/(1-x2) - 1/N1)
    E = ((eta**2)*(s**2)).exp()
    zs = acb_series([1+s, 1]).zeta(); co = zs.coeffs(); z0 = co[0]; z1 = co[1]
    f = A*E/((1-s).zeta()*z0*s**4)
    D = dA + eta**2*s - z1/z0 - 2/s
    return f, D
def mk(which, conj=False):
    def g(t, analytic):
        tt = t.conjugate() if conj else t
        f, D = parts(tt)
        v = f if which == 0 else f*D
        return v.conjugate() if conj else v
    return g
def up(x): return float(x.mid()) + float(x.rad())
res = {}
for which in (0, 1):
    g, gc = mk(which), mk(which, True)
    I = acb.integral(g, -T, T, rel_tol=1e-8, eval_limit=10**6)/(2*arb.pi())
    Q2 = acb.integral(lambda t, an: g(t, an)*gc(t, an)*(1 + t*t), -T, T, rel_tol=1e-6, eval_limit=10**6)
    Q = math.sqrt(math.pi*up(abs(Q2.real)))
    res[which] = (I, Q); print(which, I, "int|g| <=", Q, flush=True)
m = 2 - kap
tau = 4*1.3*1.26*P**(1-m)/((m-1)*math.log(P))
tau1 = 1.1*1.26*(P**-1.0)*(1 + 1/math.log(P))*2
I0, Q0 = res[0]; I1, Q1 = res[1]
e0 = (math.exp(tau) - 1)*Q0/(2*math.pi)
e1 = ((math.exp(tau) - 1)*Q1 + tau1*Q0*math.exp(tau))/(2*math.pi)
far = 2*30*14.57*3.4*math.exp(eta**2*kap**2 - eta**2*T*T)/(3*T**3)/(2*math.pi)
g_ = 0.5772156649015329
k1 = I0.real; k0 = I1.real + g_*I0.real
print("k1 = %s  +- %.2e" % (k1, e0 + far))
print("k0 = %s  +- %.2e" % (k0, e1 + far*(10 + eta**2*T2) + g_*(e0 + far)))
