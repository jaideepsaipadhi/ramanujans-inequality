# Rigorous box bounds (arb ball arithmetic): for sigma in [s0, s1], 0 <= t <= Tmax (conjugate symmetry covers t < 0):
#   Mz = sup |zeta(sigma+it)|,  Mi = sup 1/|zeta(sigma+it)|   (zeta has no zeros there: sigma > 1/2 and t << 3e12, RH verified).
# Cells of size hs x ht are enclosed by complex balls; acb.zeta returns an enclosure valid on the whole ball.
# usage: python3 box.py s0 s1 Tmax  -> prints the bounds and the worst cells.  Excludes |s-1| < 0.3 (pole) from Mz: there
#   we use |zeta(s)| <= 1/|s-1| + 1 (valid for Re s > 0 near 1? no -> we instead bound (s-1)zeta(s) on the cell).
from flint import acb, arb, ctx
import sys, math
ctx.prec = 64
s0, s1, Tmax = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
hs = 0.05; ht = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
Mz = 0.0; Mi = 0.0; Mpz = 0.0; wz = wi = None
ns = int(round((s1 - s0)/hs)); nt = int(math.ceil(Tmax/ht))
for i in range(ns):
    sc = s0 + (i + 0.5)*hs
    for j in range(nt):
        tc = (j + 0.5)*ht
        z = acb(arb(sc, hs/2), arb(tc, ht/2))
        v = (z - 1)*z.zeta()           # entire; |zeta| = |(s-1)zeta|/|s-1|
        a = abs(v); d = abs(z - 1)
        up = (a/d.lower()).upper() if d.lower() > 0 else arb('inf')
        lo = (a.lower()/d.upper())
        up = float(up.mid() + up.rad()) if up.is_finite() else float('inf')
        if lo <= 0: raise SystemExit("cell (%.3f, %.3f): no lower bound for |zeta|" % (sc, tc))
        inv = float((1/lo).upper()) if hasattr(1/lo, 'upper') else 1/float(lo)
        if up > Mz: Mz, wz = up, (sc, tc)
        if inv > Mi: Mi, wi = inv, (sc, tc)
print("box sigma in [%.2f, %.2f], t in [0, %.0f]:  sup|zeta| <= %.6g at %s   sup|1/zeta| <= %.6g at %s" % (s0, s1, Tmax, Mz, wz, Mi, wi), flush=True)
