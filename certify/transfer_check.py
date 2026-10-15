"""Transfer lemma check: certified psi bounds  ==>  |pi(y) - li(y)| <= eps_max(log y + 1) * y / log y  for log y in [57, 3159).

Inputs:
 * scan logs (lines 'L1 L2 q'):  |psi(t) - t| <= q * eps_max(L2 + 1)/(1 + 2/L1) * t  for log t in [L1, L2]   (Lemma S + Buthe Props 2-3),
   plus the additive terms of the smoothing/explicit formula: 2 + e^{2eps} log(e^eps t)(2.01 eps sqrt t + log log(2t^2)/2), eps <= 1e-4;
 * Buthe 2018 Thm 2: |psi(t) - t| <= 0.94 sqrt t, 11 < t <= 1e19;   Rosser-Schoenfeld: 0 <= psi - theta <= 1.42620 sqrt t;
 * pi(y) - li(y) = (theta(y) - y)/log y + int_2^y (theta(t) - t)/(t log^2 t) dt + 2/log 2 - li(2).
Rigorous upper sums on a grid of step h in u = log t; the requirement is checked on each grid cell with the numerator at the cell
start (all bounds decrease) and eps_max at the cell end (eps_max decreasing)."""
import sys, glob, mpmath as mp
mp.mp.dps = 40
def eps_max(L):
    L = mp.mpf(L); x = mp.e**L; li_x, li_xe = mp.li(x), mp.li(x/mp.e)
    g = (li_x**2 - (mp.e*x/L)*li_xe)/(x/L)**2; b = 2*li_x/(x/L) + L/(L - 1)
    return (-b + mp.sqrt(b*b - 4*g))/2
pieces = []
for f in sys.argv[1:]:
    for line in open(f):
        p = line.split()
        if len(p) == 3:
            try: pieces.append((float(p[0]), float(p[1]), float(p[2])))
            except ValueError: pass
pieces.sort()
def ecert(a, b):     # max certified relative psi error over u in [a, b] (pieces covering it)
    v = None
    for (L1, L2, q) in pieces:
        if L2 > a and L1 < b:
            e = mp.mpf(q)*eps_max(L2 + 1)/(1 + mp.mpf(2)/L1)
            v = e if v is None else max(v, e)
    return v
h = mp.mpf('0.05'); RS = mp.mpf('1.4262'); C0 = 2/mp.log(2) - mp.li(2)
BUTHE_TOP = mp.log(mp.mpf(10)**19)
# integral from 2 to e^43: |theta - t| <= (0.94 + 1.4262) sqrt t  (t >= 11; below 11 trivially |theta - t| <= 11),
# int_2^X t^{-1/2}/log^2 t <= 2 sqrt X/log(2)^2
X0 = mp.e**43
I = (mp.mpf('0.94') + RS)*2*mp.sqrt(X0)/mp.log(2)**2 + 11*mp.mpf(11)
u = mp.mpf(43); worst = (0, None); cache = {}
while u < 3159:
    a, b = u, u + h
    if b <= BUTHE_TOP: e = mp.mpf('0.94')*mp.e**(-a/2)
    else:
        e = ecert(a, b)
        if e is None: raise SystemExit('uncovered at %s' % a)
    e += (2 + mp.e**mp.mpf('2e-4')*(b + mp.mpf('1e-4'))*(mp.mpf('2.01e-4')*mp.e**(b/2) + mp.log(mp.log(2) + 2*b)/2))*mp.e**(-a)
    eth = e + RS*mp.e**(-a/2)                       # relative bound for |theta - t|/t on the cell
    if a >= 57:
        # |pi - li| at any y in the cell, relative to y/log y:  eth + (log y / y)(I(a) + cell part + |C0|)
        cell = eth*(mp.e**b - mp.e**a)/a**2
        num = eth + (b/mp.e**a)*(I + cell + abs(C0))
        r = num/eps_max(b + 1)
        if r > worst[0]: worst = (r, a)
    I += eth*(mp.e**b - mp.e**a)/a**2
    u = b
print('max over log y in [57, 3159) of |pi - li| bound / (eps_max(log y + 1) y/log y) = %s at log y = %s' % (mp.nstr(worst[0], 8), mp.nstr(worst[1], 6)))
print('TRANSFER OK' if worst[0] < 1 else 'TRANSFER FAILS')
