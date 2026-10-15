# S_t = sum_n Psi_t(n)^2 e^{-n/X}/n,  Psi_t = sum_{d|n} mu(d) rho_t(d),  rho_t = rho - rho_s  (Gaussian smoothing in log d, width eta).
# |rho_t(d)| <= (eta/Delta) [psi((log d - log U)/eta) + psi((log d - log V)/eta)],  psi(x) = int_{|x|}^oo Phibar = phi(x) - |x| Phibar(|x|).
# Lemma E + Selberg form: S_t <= (eps/e)^eps X^eps zeta(1+eps) (eta/Delta)^2 sum_m mu^2(m)/phi(m) beta(m)^2,
#   beta(m) = sum_K [ eta I_K(m) + 2 mx_K(m) ]  (bound for sum_{r>=1} psi((log mr - L_K)/eta)/r; f unimodal since psi log-concave).
# m-sum: blocks [u, u+du] in log m, mass <= du + 8 e^{-u/2} (Ramare-Akhilesh), times the block max of beta^2; m <= 1000 exactly.
import math, sys
from math import erfc, exp, sqrt, log, pi
def Phibar(x): return 0.5*erfc(x/sqrt(2))
def phi(x): return exp(-x*x/2)/sqrt(2*pi)
def psi(x): x = abs(x); return phi(x) - x*Phibar(x)
def Psi1(y):  # int_y^oo psi, y real
    if y >= 0:  # int_y^oo psi = int_y^oo (phi - x Phibar) = ... closed form: ((1+y^2) Phibar(y) - y phi(y))/2
        return ((1 + y*y)*Phibar(y) - y*phi(y))/2
    return 0.5 - Psi1(-y)
def beta_K(u, LK, eta):
    A = LK - u
    if A >= 0: return eta*Psi1(-A/eta) + 2*0.3990*exp(eta*eta/2)*exp(-A)
    return eta*Psi1(-A/eta) + 2*psi(A/eta)
def St_bound(lU, lV, lX, eta):
    D = lV - lU
    beta = lambda u: beta_K(u, lU, eta) + beta_K(u, lV, eta)
    # exact small m
    N0 = 1000; mu = [1]*(N0+1); ph = list(range(N0+1)); sq = [True]*(N0+1)
    for p in range(2, N0+1):
        if all(p % q for q in range(2, int(p**0.5)+1)):
            for m in range(p, N0+1, p): mu[m] = -mu[m]; ph[m] = ph[m]//p*(p-1)
            for m in range(p*p, N0+1, p*p): sq[m] = False
    S = sum(beta(log(m))**2/ph[m] for m in range(1, N0+1) if sq[m])
    du = 0.01; u = log(N0); umax = lV + 60*eta
    while u < umax:
        b = 0.0          # each beta_K is unimodal (peak at its kink): block max = max of endpoints, or the kink if inside
        for LK in (lU, lV):
            bk = max(beta_K(u, LK, eta), beta_K(u+du, LK, eta))
            if u <= LK <= u+du: bk = max(bk, beta_K(LK, LK, eta))
            b += bk
        S += b*b*(du + 8*exp(-u/2)); u += du
    S += (beta(umax)**2)*(1 + 8)          # tail beyond umax: beta decays like psi(60) ~ 0 ; generous
    best = 1e30
    for a in [0.25*k for k in range(1, 40)]:
        eps = a/lX
        z = 1/eps + 1                       # zeta(1+eps) <= 1/eps + 1
        best = min(best, (eps/math.e)**eps*exp(a)*z*(eta/D)**2*S)
    return best
if __name__ == "__main__":
    u, tau, x, L, eta = [float(a) for a in sys.argv[1:6]]
    print("S_t <=", St_bound(u*L, tau*u*L, x*L, eta))

# ---------------- closed-form rigorous version used by lemmaD (Oct 13) ----------------
# Sigma := sum_m mu^2(m)/phi(m) beta(m)^2 with beta = b(lU - log m) + b(lV - log m), where (A = L_K - log m)
#   b(A) <= eta/2 + c3 e^{-A} (A >= 0, c3 = 0.8 e^{eta^2/2}),   b(A) = eta Psi1(|A|/eta) + 2 psi(|A|/eta) (A < 0).
# Split m <= M0 = 1e6 (A >= lU - u0 >= 26): beta <= beta0 = eta + c3 (e^{-(lU-u0)} + e^{-(lV-u0)}), mass <= u0 + 1.3326 + 4/sqrt(M0).
# m > M0, h = beta^2 as a function of u = log m: partial summation with G(x) = log x + 1.3326 + R(x), |R| <= 4/sqrt(x):
#   sum <= int_{u0}^oo h du + |R(M0)| h(u0) + sup|R| TV(h);  int h <= 2 (J_U + J_V);  TV(h) <= 6 max h (<= 3 maxima),
#   J_K <= (eta^2/4)(L_K - u0) + eta c3 + c3^2/2 + eta Q0,  Q0 = int_0^oo (eta Psi1(y) + 2 psi(y))^2 dy,  max h <= (eta + 2 c3)^2.
def _Q0(eta):
    tot = 0.0; y = 0.0; h = 1e-3
    while y < 12:
        tot += (eta*Psi1(y) + 2*psi(y))**2*h; y += h        # integrand decreasing in y >= 0: left endpoints
    return tot + 1e-30
_Q0c = {}
def St_rig(lU, lV, lX, eta):
    if eta not in _Q0c: _Q0c[eta] = _Q0(eta)
    Q0 = _Q0c[eta]; D = lV - lU
    u0 = math.log(1e6); c3 = 0.8*math.exp(eta*eta/2)
    beta0 = eta + c3*(math.exp(-(lU - u0)) + math.exp(-(lV - u0)))
    small = (u0 + 1.3326 + 4/1000.0)*beta0**2
    J = lambda LK: (eta*eta/4)*(LK - u0) + eta*c3 + c3*c3/2 + eta*Q0
    mh = (eta + 2*c3)**2
    Sig = small + 2*(J(lU) + J(lV)) + (4/1000.0)*(mh + 6*mh)
    best = 1e30
    for k in range(1, 40):
        a = 0.25*k; eps = a/lX
        best = min(best, (eps/math.e)**eps*math.exp(a)*(1/eps + 1)*(eta/D)**2*Sig)
    rho_t1 = 0.5*math.erfc(lU/eta/math.sqrt(2))            # Psi_t(1) = 1 - rho_s(1) <= P(Z < -lU/eta)
    return best + rho_t1**2
