"""Ball-arithmetic certificates for the two Gamma integrals used by Lemma D (code/lemmaD.py, consts()):
  I_Gamma = sup_{sigma in [-1/2, -0.47]} int_R |Gamma(sigma+it)| (1+|t|) dt                                  <= 5.922
  I_J     = sup_{c in [-1/2, -0.44]}     int_R |Gamma(c+it)| (1+|t|/3)^{1/6} (1 + log(1+|t|/3)/log 3) dt      <= 4.72
Step 1 (monotonicity): d/dsigma log|Gamma(sigma+it)| = Re psi(sigma+it) > 0 on [-1/2, -0.44] x [0, 60], checked on arb boxes
  with adaptive subdivision.  Hence on |t| <= 60 the integrand is largest at the right end of the sigma-range.
Step 2: int_0^60 at the right end by rigorous quadrature (acb.integral; |Gamma|^2 = Gamma(s+it) Gamma(s-it) is analytic in t).
Step 3: |t| > 60 for any sigma in the range: |Gamma(s+it)| <= Gamma(1+s) |Gamma(1+it)|/t (|Gamma(x+it)|/Gamma(x) increases
  in x), Gamma(1+s) <= sqrt(pi), pi t/sinh(pi t) <= 2.0001 pi t e^{-pi t}, weights <= 1 + t <= 2t; integral < 12 e^{-15 pi}."""
from flint import arb, acb, ctx
ctx.prec = 80
PI = arb.pi()
def repsi_positive(slo, shi, tlo, thi, depth=0):
    s = arb((slo + shi)/2, (shi - slo)/2); t = arb((tlo + thi)/2, (thi - tlo)/2)
    v = acb(s, t).digamma().real
    if v > 0: return 1
    if depth > 40: raise RuntimeError('Re psi not certified positive near sigma=%g t=%g' % (slo, tlo))
    sm, tm = (slo + shi)/2, (tlo + thi)/2
    if (shi - slo) > (thi - tlo): return repsi_positive(slo, sm, tlo, thi, depth+1) + repsi_positive(sm, shi, tlo, thi, depth+1)
    return repsi_positive(slo, shi, tlo, tm, depth+1) + repsi_positive(slo, shi, tm, thi, depth+1)
def boxes(slo, shi, T=60.0, n=600):
    return sum(repsi_positive(slo, shi, T*i/n, T*(i+1)/n) for i in range(n))
def quad(sig, weight):
    s = arb(sig)
    def f(t, analytic):
        g2 = acb(s, 0) + acb(0, 1)*t
        z = g2.gamma()*(acb(s, 0) - acb(0, 1)*t).gamma()
        r = z.sqrt(analytic=analytic)
        return r*weight(t, analytic)
    return acb.integral(f, 0, 60, rel_tol=arb(2)**-40).real
TAIL = 12*(-15*PI).exp()
if __name__ == "__main__":
    print('Re psi > 0 boxes (sigma in [-0.5,-0.44], t in [0,60]):', boxes(-0.5, -0.44))
    IG = 2*(quad(-0.47, lambda t, a: 1 + t) + TAIL)
    print('I_Gamma <=', IG.upper(), '(used 5.922):', 'OK' if IG < arb('5.922') else 'FAIL')
    l3 = arb(3).log()
    def wJ(t, a):
        lg = (1 + t/3).log(analytic=a)
        return (lg/6).exp()*(1 + lg/l3)
    IJ = 2*(quad(-0.44, wJ) + TAIL)
    print('I_J     <=', IJ.upper(), '(used 4.72):', 'OK' if IJ < arb('4.72') else 'FAIL')
