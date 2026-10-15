# Ball-arithmetic upper bounds for the two integrals in bty2 that were done with mpmath:
#   I4(nu) = int_0^D |w''''(u)| e^{-nu u} du      (Cw constant; monotone decreasing in nu -> certify on an integer nu grid)
#   I0(x0) = int_R |Re psi(1/4 + iy/2)| (x0^2 + y^2)^{-3/2} dy    (Lemma 8 (32))
# Method: split at sign changes (located numerically), certify sign-constancy of the integrand piece by ball bisection,
# integrate the analytic integrand with acb.integral, sum |pieces|.
from flint import arb, acb, ctx
import mpmath as mp, math, json
ctx.prec = 100
th = arb(11338)/10000; kk = th.tan(); ct = 1/kk; s2 = 1/(th.cos()**2); D = 2*th*ct
def w4(u):   # 4th derivative of w (closed form), u acb or arb
    c = (u*kk).cos(); s = (u*kk).sin(); A = th*ct - u/2
    # d^4/du^4 [A cos(ku)] = A k^4 cos(ku) + 4*(-1/2)*(k^3 sin(ku))*(-1)... derive: f = A c; f' = -c/2 - A k s; f'' = k s + ... use generic Leibniz:
    # (A c)^{(4)} = A c^{(4)} + 4 A' c^{(3)} ;  c^{(3)} = k^3 s, c^{(4)} = k^4 c, A' = -1/2
    f4 = A*kk**4*c + 4*(-arb(1)/2)*(kk**3*s)
    # sin(phi - k u)^{(4)} = k^4 sin(phi - k u)
    return s2*(s2*f4 + kk**4*(2*th - u*kk).sin()/(2*th).sin() - 2*kk**4*(th - u*kk).sin()/th.sin())
def sign_const(f, a, b, depth=0):
    v = f(arb((a+b)/2, (b-a)/2))
    if v > 0: return 1
    if v < 0: return -1
    if depth > 60: return 0
    m = (a+b)/2; s1 = sign_const(f, a, m, depth+1); s2_ = sign_const(f, m, b, depth+1)
    return s1 if s1 == s2_ and s1 != 0 else 0
EPS = 1e-7
def abs_integral(f_acb, f_arb, roots, a, b):
    pts = [a] + [r for r in roots if a < r < b] + [b]; tot = arb(0)
    for p, q in zip(pts[:-1], pts[1:]):
        eps = EPS
        sg = sign_const(f_arb, p + eps, q - eps)
        assert sg != 0, (p, q)
        val = acb.integral(lambda x, an: f_acb(x), arb(p + eps), arb(q - eps), rel_tol=1e-20).real
        tot += abs(val)
    return tot
def roots_num(f, a, b, n=20000):
    xs = [a + (b-a)*i/n for i in range(n+1)]; vs = [f(x) for x in xs]; rs = []
    for i in range(n):
        if vs[i] == 0 or vs[i]*vs[i+1] < 0: rs.append(float(mp.findroot(f, (xs[i], xs[i+1]), solver='anderson')))
    return rs
mw4 = lambda x: float(w4(arb(x)).mid())
Dm = float(D.mid())
R4 = roots_num(mw4, 0, Dm)
def I4(nu):
    nu = arb(nu)
    fa = lambda x: w4(x)*(-nu*x).exp(); fr = lambda x: w4(x)*(-nu*x).exp()
    # the tiny intervals around roots (width 2e-9) are covered by the sup bound below
    gap = arb(0)
    for r in R4:
        gap += 2*EPS*abs(w4(arb(r, EPS)))
    gap += 2*EPS*(abs(w4(arb(0, EPS))) + abs(w4(arb(Dm, EPS))))
    return abs_integral(fa, fr, R4, 0.0, Dm) + gap
if __name__ == "__main__":
    print('roots of w4:', R4)
    out = {}
    for nu in [0] + list(range(30, 201)):
        v = I4(nu); out[nu] = float(v.upper()); 
        if nu in (0, 30, 60, 100, 150, 200): print(nu, v, flush=True)
    json.dump(out, open('I4cert.json', 'w'))
def I0(x0):
    x0 = arb(x0); q = arb(1)/4
    def f(y):
        return ((q + acb(0, 1)*y/2).digamma() + (q - acb(0, 1)*y/2).digamma())/2/((x0*x0 + y*y)**(arb(3)/2))
    mf = lambda y: float(f(acb(y)).real.mid())
    rts = roots_num(lambda y: float((((q + acb(0,1)*arb(y)/2).digamma()).real).mid()), 0.0, 50.0, 5000)
    Y = 1000.0
    body = abs_integral(lambda y: f(y), lambda y: f(acb(y)).real, rts, 0.0, Y)
    gap = arb(0)
    for r in [0.0] + rts + [Y]: gap += 2*EPS*abs(f(acb(arb(r, EPS))).real)
    tail = (2*arb(Y).log() + 1)/(4*arb(Y)**2)
    return 2*(body + gap + tail)
