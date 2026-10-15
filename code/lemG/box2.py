# Rigorous box bounds for zeta on sigma in [s0,s1], t in [0,Tmax] via Taylor models of the ENTIRE function f(s) = (s-1) zeta(s).
# For a cell with centre c and half-diagonal r: f(c+h) = sum_{k<=K} a_k h^k + R, |R| <= M (r/Rd)^{K+1}/(1 - r/Rd), where
# M >= sup_{|s-c| = Rd} |f(s)|; we use the elementary bound (sigma > 0): zeta(s) = s/(s-1) - s int_1^oo {u} u^{-s-1} du,
# |zeta(s)| <= |s|/|s-1| + |s|/sigma, hence |f(s)| <= |s| + |s||s-1|/sigma.   Cells failing the lower bound are split.
from flint import acb, arb, acb_series, ctx
import sys, math
ctx.prec = 80
ctx.cap = 16
s0, s1, Tmax = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
K = 14; Rd = 0.2
def cellbound(sc, tc, hs, ht, depth=0):
    r = math.hypot(hs, ht)/2*1.0000001
    c = acb(sc, tc)
    ser = (acb_series([c, 1], prec=K+1) - 1)*acb_series([c, 1], prec=K+1).zeta()
    co = [ser.coeffs()[k] for k in range(K+1)]
    smin = sc - Rd; assert smin > 0
    smax_abs = math.hypot(abs(sc)+Rd, abs(tc)+Rd); sm1 = math.hypot(abs(sc-1)+Rd, abs(tc)+Rd)
    M = smax_abs + smax_abs*sm1/smin
    q = r/Rd; rem = M*q**(K+1)/(1-q)
    a0 = float(abs(co[0]).upper()); a0l = float(abs(co[0]).lower())
    tail = sum(float(abs(co[k]).upper())*r**k for k in range(1, K+1)) + rem
    up = a0 + tail; lo = a0l - tail
    dlo = max(math.hypot(sc-1, tc) - r, 0.0); dup = math.hypot(sc-1, tc) + r
    if lo <= 0:
        if depth > 6: raise SystemExit("cannot certify near (%.4f, %.4f)" % (sc, tc))
        out = [cellbound(sc + dx*hs/4, tc + dy*ht/4, hs/2, ht/2, depth+1) for dx in (-1, 1) for dy in (-1, 1)]
        return max(o[0] for o in out), max(o[1] for o in out)
    return (up/dlo if dlo > 0 else float('inf')), dup/lo          # sup |zeta|, sup 1/|zeta| on the cell
hs = 0.05; ht = 0.05
Mz = Mi = 0.0
ns = int(round((s1-s0)/hs)); nt = int(math.ceil(Tmax/ht))
for j in range(nt):
    tc = (j+0.5)*ht
    for i in range(ns):
        sc = s0 + (i+0.5)*hs
        if math.hypot(sc-1, tc) < 0.12:     # pole region: only 1/zeta needed there, and |1/zeta| <= |s-1| bound via f
            pass
        a, b = cellbound(sc, tc, hs, ht)
        if math.hypot(sc-1, tc) >= 0.12: Mz = max(Mz, a)
        Mi = max(Mi, b)
    if j % 2000 == 0: print("t<=%.1f  sup|zeta|(away from pole) %.5g  sup|1/zeta| %.5g" % (tc, Mz, Mi), flush=True)
print("FINAL box sigma[%.2f,%.2f] t[0,%.0f]: sup|zeta| (|s-1|>=0.1) <= %.6g, sup|1/zeta| <= %.6g" % (s0, s1, Tmax, Mz, Mi), flush=True)
