"""Ball-arithmetic version of certify/transfer_check.py (Proposition transfer): certified psi bounds from the arb scan
==>  |pi(y) - li(y)| <= eps_R(log y + 1) y/log y and <= eps_R(log y) y/log y  for log y in [57, 3159).

Inputs: scan_arb_*.txt, lines 'L1 L2 q' with q an upper bound for the ratio of Proposition master to
tol(L1, L2) = eps_R^-(L2 + 1)/(1 + 2/L1), eps_R^- the lower end of the epsR_cert.py ball at L2 + 1 (exactly as in scan_arb.py).
Hence |psi(t) - t|/t <= q tol(L1, L2) + (a(t) + 2)/t for log t in [L1, L2].  Other inputs as in transfer_check.py
(Buthe 0.94 sqrt t below 1e19, Rosser-Schoenfeld 1.4262 sqrt t for psi - theta).  eps_R on each cell [a, a + 1/20]
is the lower end of epsR_cert.eps_ball(a, a + 1/20) (ball over the whole cell: no monotonicity used)."""
import sys, os, glob
from flint import arb, ctx
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
import epsR_cert as ER
ctx.prec = 120
pieces = []
for f in (sys.argv[1:] or sorted(glob.glob(os.path.join(HERE, 'scan_arb_*.txt')))):
    for line in open(f):
        p = line.split()
        if len(p) >= 3:
            L1, L2, q = float(p[0]), float(p[1]), float(p[2])
            tol = arb(ER.eps_ball(L2 + 1, L2 + 1)[0].lower())/(1 + 2/arb(L1))
            pieces.append((L1, L2, arb(q)*tol))
pieces.sort(key=lambda z: (z[0], z[1]))
def ecert(a, b):
    v = None
    for L1, L2, e in pieces:
        if L2 > a and L1 < b:
            v = e if v is None or e.upper() > v.upper() else v
    return v
_EPS = {}
def epsR_low(j):          # cell [43 + j/20, 43 + (j+1)/20]
    if j not in _EPS:
        lo, hi = 43 + j/20, 43 + (j + 1)/20          # doubles; the ball below covers [lo, hi] (endpoints exact or rounded outward)
        import numpy as np
        e, g = ER.eps_ball(float(np.nextafter(lo, -np.inf)), float(np.nextafter(hi, np.inf)))
        assert g < 0
        _EPS[j] = arb(e.lower())
    return _EPS[j]
h = arb(1)/20; RS = arb('1.4262'); C0 = 2/arb(2).log() - arb(2).li()
BUTHE_TOP = (arb(10)**19).log()
X0 = arb(43).exp()
I = (arb('0.94') + RS)*2*X0.sqrt()/arb(2).log()**2 + 11*11
worst = (arb(0), None); n = int(round((3159 - 43)*20))
for j in range(n):
    a = arb(43) + arb(j)/20; b = a + h
    if b < BUTHE_TOP: e = arb('0.94')*(-a/2).exp()
    else:
        e = ecert(float(a.mid()) - 1e-9, float(b.mid()) + 1e-9)
        if e is None: raise SystemExit('uncovered at %s' % a)
    e += (2 + arb('2e-4').exp()*(b + arb('1e-4'))*(arb('2.01e-4')*(b/2).exp() + ((arb(2).log() + 2*b).log())/2))*(-a).exp()
    eth = e + RS*(-a/2).exp()
    if j >= 14*20:                                       # a >= 57
        cell = eth*(b.exp() - a.exp())/a**2
        num = eth + (b/a.exp())*(I + cell + abs(C0))
        den = epsR_low(j + 20)                            # y = x/e: eps_R at log y + 1
        if j >= 15*20 and epsR_low(j) < den: den = epsR_low(j)   # y = x (needed for log x >= 58): eps_R at log y
        r = num/den
        if r.upper() > worst[0].upper(): worst = (r, a)
    I += eth*(b.exp() - a.exp())/a**2
print('max over log y in [57, 3159) of the |pi - li| bound / (eps_R y/log y):', worst[0].upper(), 'at log y =', worst[1].mid())
print('TRANSFER OK' if worst[0] < 1 else 'TRANSFER FAILS')
