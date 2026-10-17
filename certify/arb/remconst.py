"""Ball-arithmetic certification of the cell-independent constants in the Lemma G contour remainder (code/lemG/remainder.py)
and in the S_t bound (code/lemG/stbound.py).  Output: remconst.json (arb balls as strings).

Parameters are the binary doubles used by the original code (eta = 0.02 etc.), converted exactly to arb, so every contour
and every real part below is the same one used in remainder.py / kcert2.py.

Abound: the sup over the 3-torus for p <= 100 is evaluated on the same 96^3 grid with the same cell enclosure, in
interval arithmetic on doubles with outward rounding after every operation (IEEE +,-,*,/,sqrt are correctly rounded;
cos/sin of the grid angles come from arb).  Primes 100 < p <= 1e6 and the tail are done in arb."""
import json, math, sys, time
import numpy as np
from flint import arb, acb, ctx
ctx.prec = 128
F = lambda x: arb(float(x))                       # exact conversion of a double
eta, T2, W0, c, cp, s0, th, b, thp, kap = [F(v) for v in (0.02, 400.0, 40.0, 0.02, 0.01, 0.01, 0.40, 0.38, 0.30, 0.38)]
one = arb(1); PI = arb.pi()
inf = np.inf

# ---------------- outward-rounded interval arithmetic on numpy doubles ----------------
def up(x): return np.nextafter(x, inf)
def dn(x): return np.nextafter(x, -inf)
def fup(a):  # arb -> double >= a
    return float(np.nextafter(np.nextafter(float(a.upper()), inf), inf))
def fdn(a):
    return float(np.nextafter(np.nextafter(float(a.lower()), -inf), -inf))
class I:
    __slots__ = ('lo', 'hi')
    def __init__(s, lo, hi): s.lo, s.hi = lo, hi
    def __add__(s, o):
        if not isinstance(o, I): o = I(o, o)
        return I(dn(s.lo + o.lo), up(s.hi + o.hi))
    __radd__ = __add__
    def __neg__(s): return I(-s.hi, -s.lo)
    def __sub__(s, o):
        if not isinstance(o, I): o = I(o, o)
        return s + (-o)
    def __rsub__(s, o): return (-s) + o
    def __mul__(s, o):
        if not isinstance(o, I): o = I(o, o)
        ps = [s.lo*o.lo, s.lo*o.hi, s.hi*o.lo, s.hi*o.hi]
        return I(dn(np.minimum(np.minimum(ps[0], ps[1]), np.minimum(ps[2], ps[3]))),
                 up(np.maximum(np.maximum(ps[0], ps[1]), np.maximum(ps[2], ps[3]))))
    __rmul__ = __mul__
    def absup(s): return np.maximum(np.abs(s.lo), np.abs(s.hi))
    def sqlo(s):   # lower bound of x^2
        m = np.where((s.lo <= 0) & (s.hi >= 0), 0.0, np.minimum(np.abs(s.lo), np.abs(s.hi)))
        return dn(m*m)
def aI(a): return I(fdn(a), fup(a))
def sqrt_up(x): return up(np.sqrt(x))
def sqrt_dn(x): return dn(np.sqrt(np.maximum(x, 0.0)))

def primes(n):
    sv = bytearray([1])*(n+1); sv[0:2] = b'\x00\x00'
    for i in range(2, int(n**0.5)+1):
        if sv[i]: sv[i*i::i] = bytearray(len(sv[i*i::i]))
    return [i for i in range(n+1) if sv[i]]
PR = primes(10**6)

def Abound(e1, e2, e12, n=96):
    """rigorous upper bound (arb) for sup |A| in remainder.Abound (same grid n = 96, same cell enclosure)."""
    h = 2*PI/n; e = h/2
    ang = [k*h for k in range(n)]
    C = [aI(a.cos()) for a in ang]; S = [aI(a.sin()) for a in ang]
    Clo = np.array([x.lo for x in C]); Chi = np.array([x.hi for x in C])
    Slo = np.array([x.lo for x in S]); Shi = np.array([x.hi for x in S])
    i1, i2, iw = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing='ij')
    ph = (i1 + i2 - iw) % n                                   # phase of z: a1 + a2 - aw (mod 2 pi), exact on the grid
    c1, s1 = I(Clo[i1], Chi[i1]), I(Slo[i1], Shi[i1])
    c2, s2 = I(Clo[i2], Chi[i2]), I(Slo[i2], Shi[i2])
    c3, s3 = I(Clo[ph], Chi[ph]), I(Slo[ph], Shi[ph])
    tot = arb(0)
    for p in PR:
        if p > 100: break
        r1a, r2a, rza = arb(p)**(-e1), arb(p)**(-e2), arb(p)**(-e12)
        r1, r2, rz = aI(r1a), aI(r2a), aI(rza)
        E = fup(e)
        # x = r e^{-i a}:  1 - x1 - x2 + z = (1 - r1 c1 - r2 c2 + rz c3) + i (r1 s1 + r2 s2 - rz s3)
        re = 1 - r1*c1 - r2*c2 + rz*c3
        im = r1*s1 + r2*s2 - rz*s3
        N1 = sqrt_up(up(up(re.absup()**2) + up(im.absup()**2)))
        enc1 = fup((r1a + r2a)*e + 3*rza*e)
        N1 = up(N1 + enc1)
        re2 = 1 - rz*c3; im2 = rz*s3
        N2 = up(sqrt_up(up(up(re2.absup()**2) + up(im2.absup()**2))) + fup(3*rza*e))
        # |1 - x1|^2 = (1 - r1 c1)^2 + (r1 s1)^2 (lower bound)
        d1 = 1 - r1*c1; q1 = r1*s1
        D1 = dn(sqrt_dn(dn(d1.sqlo() + q1.sqlo())) - fup(r1a*e))
        d2 = 1 - r2*c2; q2 = r2*s2
        D2 = dn(sqrt_dn(dn(d2.sqlo() + q2.sqlo())) - fup(r2a*e))
        assert (D1 > 0).all() and (D2 > 0).all()
        Ap = up(up(N1*N2)/dn(D1*D2))
        m = float(Ap.max())
        tot += arb(m).log()                                   # log of an upper bound (arb log is rigorous)
    acc = arb(0)
    for p in PR:
        if p <= 100: continue
        x1, x2, z = arb(p)**(-e1), arb(p)**(-e2), arb(p)**(-e12)
        acc += (z*x1 + z*x2 + z*z + x1*x2)/((1 - x1)*(1 - x2))
    tot += acc
    m = min(float((e12+e1).lower()), float((e12+e2).lower()), float((2*e12).lower()), float((e1+e2).lower()))
    m = arb(m)                                                # a lower bound for the true minimum
    P_ = arb(10**6)
    tot += 4.2*P_**(1 - m)/(m - 1)
    return tot.exp()

def zeta(s): return s.zeta()
def izeta(s): return s.zeta()/(2*s).zeta()

def gam_int(sig, w1, n=40000):
    """2 * upper Riemann sum of |Gamma(sig+iv)| (1+v)^{w1} over [0, W0] (|Gamma| decreasing in v: left end; weight: right end)."""
    h = W0/n; tot = arb(0)
    for i in range(n):
        g = acb(sig, i*h).gamma().abs_upper()
        tot += arb(g)*((1 + (i+1)*h) if w1 else 1)
    return 2*tot*h

def J(sig, lo, hi, w1):
    """2 * int_lo^hi e^{eta^2(sig^2 - t^2)/2} (1+t)^{w1}/(sig^2+t^2) dt, upper sum (decreasing factors left, (1+t) right)."""
    top = min(hi, 3000.0); tot = arb(0); t0 = lo
    E2 = eta*eta/2
    while t0 < top:
        hh = 1e-4 if t0 < 0.1 else (1e-3 if t0 < 2 else 0.01)
        t1 = min(t0 + hh, top)
        a0, a1 = arb(t0), arb(t1)
        tot += (E2*(sig*sig - a0*a0)).exp()*((1 + a1) if w1 else 1)/(sig*sig + a0*a0)*(a1 - a0)
        t0 = t1
    if hi > 3000:
        # int_3000^oo:  (1+t)^{w1}/(sig^2+t^2) <= 2/3000 ;  int_a^oo e^{-eta^2 t^2/2} <= e^{-eta^2 a^2/2}/(eta^2 a)
        a = arb(3000)
        tot += (E2*sig*sig).exp()*arb(2)/a*(-(E2*a*a)).exp()/(eta*eta*a)
    return 2*tot

def R6(AK, Ax, MZ, MI, n=40000):
    gauss = (-(eta*eta*T2*T2)/2).exp()
    a_, b_ = kap, kap - thp
    segc = (eta*eta*((kap - thp)**2 + kap**2)/2).exp()
    iz = izeta(1 + kap - thp); AKx = AK if AK > Ax else (Ax if Ax > AK else arb(max(fup(AK), fup(Ax))))
    def Kc(v):
        seg = PI*(a_ + b_)/(a_*b_*((a_ + b_)**2 + v*v))
        return (MZ*AK*MI*iz*seg*segc + MZ*AKx*MI*MI*4*(kap + cp)*gauss/T2**4)/(2*PI)
    Ew = (eta*eta*thp*thp/2).exp()/(thp*thp)
    h = W0/n; tot = arb(0)
    for i in range(n):
        v = i*h
        tot += arb(acb(-thp, v).gamma().abs_upper())*(Kc(v) + Ew)
    return 2*tot*h, Kc(W0) + (eta*eta*thp*thp/2).exp()/(W0*W0)

# ---------------- S_t pieces ----------------
def Phibar(x): return (x/arb(2).sqrt()).erfc()/2
def phi(x): return (-(x*x)/2).exp()/(2*PI).sqrt()
def psi(x): return phi(x) - x*Phibar(x)                   # x >= 0
def Psi1(y): return ((1 + y*y)*Phibar(y) - y*phi(y))/2     # y >= 0
def Q0(n=12000):
    tot = arb(0); hh = arb(12)/n
    for i in range(n):
        y = i*hh
        tot += (eta*Psi1(y) + 2*psi(y))**2
    tot *= hh
    # tail y > 12: integrand <= (eta + 2)^2 phi(y)^2;  int_12^oo e^{-y^2}/(2 pi) <= e^{-144}/(2 pi 24)
    tot += (eta + 2)**2*(arb(-144)).exp()/(2*PI*24)
    return tot

if __name__ == "__main__":
    MZ, MI = arb(11.89), arb(14.57)
    out = {}
    t0 = time.time()
    triples = {'A1': (1 + c - th, 1 + c - cp, 1 + c - th - cp), 'Ax': (1 + cp, 1 + c - cp, one),
               'AV': (1 - b, 1 + c + b, one), 'AU': (1 + b, 1 + c - b, one), 'AK': (1 - kap, 1 + kap - thp, one)}
    which = sys.argv[1:] or list(triples) + ['rest']
    try: out = {k: arb(v) for k, v in json.load(open('remconst.json')).items()}
    except Exception: out = {}
    for k, tr in triples.items():
        if k in which:
            out[k] = Abound(*tr); print(k, out[k], 'float code: see remainder.Abound', '%.0fs' % (time.time() - t0), flush=True)
            json.dump({kk: vv.str(30, radius=True) for kk, vv in out.items()}, open('remconst.json', 'w'), indent=1)
    if 'rest' in which:
        out['gam_c_w1'] = gam_int(c, True); out['gam_c'] = gam_int(c, False); print('gam', out['gam_c_w1'], out['gam_c'], flush=True)
        out['J_mth_0_T2_w1'] = J(-th, 0.0, 400.0, True)
        out['J_mcp_0_inf_w1'] = J(-cp, 0.0, math.inf, True)
        out['J_mcp_0_inf'] = J(-cp, 0.0, math.inf, False)
        out['J_s0_T2_inf_w1'] = J(s0, 400.0, math.inf, True); print('J', [str(out[k]) for k in out if k.startswith('J')], flush=True)
        r6, kcw0 = R6(out['AK'], out['Ax'], MZ, MI); out['R6'] = r6; out['KcW0'] = kcw0; print('R6', r6, kcw0, flush=True)
        out['Q0'] = Q0(); print('Q0', out['Q0'], flush=True)
        for nm, s in [('z1c', 1 + c), ('z1c_s0_cp', 1 + c + s0 - cp), ('z2_2c', 2 + 2*c)]: out[nm] = zeta(s)
        for nm, s in [('iz1c_cp', 1 + c - cp), ('iz1c_s0', 1 + c + s0), ('iz1cp', 1 + cp), ('iz1c_b', 1 + c + b), ('iz1b', 1 + b)]:
            out[nm] = izeta(s)
        json.dump({kk: vv.str(30, radius=True) for kk, vv in out.items()}, open('remconst.json', 'w'), indent=1)
    print('done %.0fs' % (time.time() - t0))
