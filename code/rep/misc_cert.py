# Ball versions of the remaining float constants: Lemma 12 gamma constant, and w-derivatives used in Cw.
from flint import arb, arb_series, ctx
ctx.prec = 100
import bty2, math
cc = [4,-8,2,20,-9,-34,27,91,-27,-201,32,895,1949,2389,1896,949,239]
ak = [arb(sum(cc[j]*cc[j+k] for j in range(17-k))*(1 if k == 0 else 2))/14912370 for k in range(17)]
km = [arb(1), arb(-851)/859, arb(780)/859, arb(-525)/859, arb(171)/859, arb(28)/859, arb(-29)/859]
kap = sum(km, arb(0)); a = sum(ak[1:], arb(0))
def gamma_const(sig0):
    s0 = arb(sig0)
    t1 = a*kap/2*(-(2*arb.pi()).log() + arb('1e-10')) + kap/2*sum((ak[k]*arb(k).log() for k in range(1, 17)), arb(0))
    t0 = arb(0)
    for m, k in enumerate(km):
        lo = ((s0 + (2*s0-1)*m)/2 + 1).digamma(); hi = (arb(m + 3)/2).digamma()
        t0 += k/2*(hi if k > 0 else lo)
    t0 += ak[0]*(-arb.pi().log()/2)*kap
    return -(t1 + t0)
for s0 in (0.9930, 0.99323, 0.9935164):
    g = gamma_const(s0); print('G(%.7f) ball' % s0, g, ' float', bty2.gamma_const(s0), ' float <= ball.lower:', bty2.gamma_const(s0) <= float(g.lower()) + 1e-12)
th = arb(11338)/10000; kk = th.tan(); ct = 1/kk; s2 = 1/(th.cos()**2); D = 2*th*ct
ctx.cap = 6
def w(u): return s2*(s2*(th*ct-u/2)*(u*kk).cos()+2*th*ct-u+(2*th-u*kk).sin()/(2*th).sin()-2-2*(th-u*kk).sin()/th.sin())
c0 = w(arb_series([arb(0), 1])).coeffs(); cD = w(arb_series([D, 1])).coeffs()
W2_0 = abs(2*c0[2]); W3_0 = abs(6*c0[3]); W3_D = abs(6*cD[3])
print('w2(0)', W2_0, float(bty2.W2_0), ' w3(0)', W3_0, float(bty2.W3_0), ' w3(D)', W3_D, float(bty2.W3_d))
print('floats >= ball upper - 1e-15:', all(float(x) >= float(y.upper()) - 1e-15 for x, y in [(bty2.W2_0, W2_0), (bty2.W3_0, W3_0), (bty2.W3_d, W3_D)]))
