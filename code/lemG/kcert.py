# Rigorous enclosures (arb) of k1 = H(0) and k0 = H'(0) + gamma H(0), where
#   H(w) = (1/2 pi i) int_{Re s = c} A(w;-w-s,s) / (zeta(1-s) zeta(1+w+s) (w+s)^2 s^2) ds ,   c = 0.3.
# A = prod_p f_p, f_p = (1-a-b+q)(1-q)/((1-a)(1-b)), a = p^{-(1-s)}, b = p^{-(1+w+s)}, q = p^{-1}  (abc = 1/p exactly).
# Tail p > P: |log f_p| <= 1.2 (p^{-1.7+|w|} + 2 p^{-2+|w|}) (|a| <= p^{-0.7}, |b| <= p^{-1.3+|w|}); sum_{n>P} n^{-sig} <= P^{1-sig}/(sig-1).
# |t| > T: integrand bounded by Ctail/|s|^4 with Ctail = sup over the cut-off of |A||1/zeta||1/zeta| * |s|^2/|w+s|^2 checked on a grid
#   and bounded via convexity estimates (see TAIL below).
from flint import acb, arb, ctx, acb_mat
import sys
ctx.prec = 60
P = 2000
primes = [p for p in range(2, P+1) if all(p % d for d in range(2, int(p**0.5)+1))]
c = arb('0.3')
def tailfac(r):
    s1 = arb(1.7) - r; s2 = arb(2) - r
    tau = arb('1.2')*(arb(P)**(1-s1)/(s1-1) + 2*arb(P)**(1-s2)/(s2-1))
    return tau   # |log(A/A_P)| <= tau
def integrand(w):
    def f(t, analytic):
        s = acb(c, 0) + acb(0, 1)*t
        lp = None; A = acb(1)
        for p in primes:
            pp = acb(p)
            a = pp**(-(1-s)); b = pp**(-(1+w+s)); q = acb(1)/p
            A *= (1-a-b+q)*(1-q)/((1-a)*(1-b))
        return A/((1-s).zeta()*(1+w+s).zeta()*(w+s)**2*s**2)
    return f
def H(w, T):
    f = integrand(w)
    I = acb.integral(lambda t, an: f(t, an), -T, T, rel_tol=1e-6, eval_limit=200000)
    return I/(2*arb.pi())
if __name__ == "__main__":
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    print(H(acb(0), T))

def run(T=60, r=0.01, R=0.1):
    tau0 = tailfac(arb(0)); tauR = tailfac(arb(R))
    g = arb('0.5772156649015328606')
    H0 = H(acb(0), T); Hp = H(acb(r), T); Hm = H(acb(-r), T)
    # sup |H_T| on the disk |w| <= R via a ball in w
    HR = H(acb(0, 0).__class__(arb(0, R), arb(0, R)), T)
    M = abs(HR).upper()
    # Euler tail: H_T(w) computed with A_P; |A - A_P| <= |A_P| (e^tau - 1) -> relative factor on the integral bound
    k1 = H0; d = (Hp - Hm)/(2*r)
    k0 = d + g*H0
    fd_err = arb(r)**2*M/(arb(R)**3 - arb(r)**3)      # Cauchy bound on (H(r)-H(-r))/2r - H'(0)
    print("T=%d: H0=%s  H'(0)~%s  M=%s  fd_err=%s  euler tau0=%s" % (T, H0, d, M, fd_err, tau0))
    return H0, k0, fd_err, M, tauR
if __name__ == "__main__" and len(sys.argv) > 2:
    run(int(sys.argv[1]))
