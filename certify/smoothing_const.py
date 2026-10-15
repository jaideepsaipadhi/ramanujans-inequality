# Lemma S constant check (ball arithmetic): for every c in the scan grid [30, 110],
#   K(c) := 2 sqrt(c) [ |nu_c(0)|/2 + I0(c)/(200 sinh c) ] * e^{eps/2}  <=  0.7946 * 1.02,   eps <= 1e-4,
# where |nu_c(0)| = int_0^1 u eta_{c,1}(u) du,  eta_{c,1}(u) = c I0(c sqrt(1-u^2))/sinh(c)  (int_{-1}^1 eta = 2).
from flint import arb, acb, ctx
ctx.prec = 80
def nu(c):
    c = arb(c)
    f = lambda u, an: acb(u)*c*(c*(1 - acb(u)**2).sqrt()).bessel_i(0)/c.sinh()
    return acb.integral(f, arb(0), arb(1), rel_tol=1e-15).real
worst = 0
for c in range(30, 111):
    K = 2*arb(c).sqrt()*(nu(c)/2 + arb(c).bessel_i(0)/(200*arb(c).sinh()))*arb('1e-4').exp()**0.5
    worst = max(worst, float(K.upper()))
    if c % 10 == 0: print(c, K)
print('max K =', worst, ' bound used 0.7946*1.02 =', 0.7946*1.02, ' OK' if worst < 0.7946*1.02 else ' FAIL')
