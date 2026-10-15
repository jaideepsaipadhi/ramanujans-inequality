# BTY Lemma 3 (f >= 0 on [0, inf)) certified: (i) p(x) = sum kappa_m x^m > 0 on [0,1] by exact root isolation;
# (ii) w >= 0 on [0, D]: ball bisection on [0, D - e]; near D, w(u) = -w3(D)/6 (D-u)^3 + w4(xi)/24 (D-u)^4 with w3(D) < 0.
import sympy as sp
from flint import arb, arb_series, ctx
x = sp.symbols('x'); p = sum(sp.Rational(*k)*x**m for m, k in enumerate([(1,1),(-851,859),(780,859),(-525,859),(171,859),(28,859),(-29,859)]))
P = sp.Poly(p, x); print('p roots in [0,1]:', P.count_roots(0, 1), ' p(0) =', p.subs(x, 0))
ctx.prec = 100; ctx.cap = 6
th = arb(11338)/10000; kk = th.tan(); ct = 1/kk; s2 = 1/(th.cos()**2); D = 2*th*ct
def w(u): return s2*(s2*(th*ct-u/2)*(u*kk).cos()+2*th*ct-u+(2*th-u*kk).sin()/(2*th).sin()-2-2*(th-u*kk).sin()/th.sin())
e = 0.02; Dl = float(D.lower())
stack = [(0.0, Dl - e)]; n = 0; ok = True
while stack:
    a_, b_ = stack.pop(); n += 1
    if w(arb((a_+b_)/2, (b_-a_)/2)) > 0: continue
    if b_ - a_ < 1e-10: ok = False; break
    m = (a_+b_)/2; stack += [(a_, m), (m, b_)]
X = arb(Dl - e/2, e/2 + 1e-15)                     # covers [D - e, D]
c = w(arb_series([X, 1])).coeffs(); w4sup = abs(24*c[4])
w3D = 6*w(arb_series([D, 1])).coeffs()[3]
tail_ok = (-w3D/6 - w4sup/24*arb(e) > 0)
print('w > 0 on [0, D-e]:', ok, n, ' w3(D) =', w3D, ' near-D Taylor positive:', tail_ok)
print('=> f >= 0 certified' if P.count_roots(0, 1) == 0 and ok and tail_ok else 'FAIL')
