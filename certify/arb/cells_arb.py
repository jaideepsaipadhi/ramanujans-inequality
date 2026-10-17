"""Ball-arithmetic re-evaluation of every Lemma D density cell (Theorem density) from its stored parameters.

For each table cell (k, lk) with stored parameters p = (u, t, w, x) (code/lemmaD_tableGH0.01_*.pkl) we recompute in arb
    log N^+(1 - delta_k, e^lk) = log( (d4 d5/w) * coef * ntil / (c1^2 - c5) ) + 2 x delta_k lk - log(delta_k lk)
with every ingredient evaluated in ball arithmetic:
  * d4 by Lemma G (D4m + contour remainder + S_t terms), the remainder constants from remconst.json;
  * the pair coefficient: the LP value is bounded by weak duality (lambda M + sum max(0, b - lambda c)) instead of the
    greedy knapsack; -zeta'/zeta(1+eta) is evaluated exactly (no 1/eta - GAP);  |Gamma| ratios via acb;
  * ntil and the window counts n(w) as floor of a certified upper bound;  the window sum in I' gets an explicit tail y > 60;
  * c1, c5, d5, B, J as in consts(), with I_Gamma = 5.922 and 4.72 certified by gammaint_cert.py.
Discrete choices (kappa, eta, the Gaussian-smoothing exponent a, the c in _M) are taken from a floating-point search;
any choice gives a valid bound, so only the evaluation needs to be rigorous.
Output: cells_arb.json  {"k,lk": upper bound for log N^+ (string, rounded up)}, plus the float value for comparison."""
import sys, os, json, math, pickle, glob, time
from fractions import Fraction as Fr
from flint import arb, acb, arb_series, ctx
HERE = os.path.dirname(os.path.abspath(__file__)); CODE = os.path.join(HERE, '..', '..', 'code')
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(CODE, 'lemG'))
os.environ['LEMMAG'] = '1'
import lemmaD as D, lemmaC as C, ramanujan_upper as RU
ctx.prec = 128
PI = arb.pi(); GAM = arb.const_euler(); A = arb
def fup(a): return float(a.upper()) if float(a.upper()) == float(a.upper()) else float('inf')
RC = {k: arb(v) for k, v in json.load(open(os.path.join(HERE, 'remconst.json'))).items()}
F_ = lambda v: arb(float(v))
eta, T2, W0, c, cp, s0, th, b, thp, kap = [F_(v) for v in (0.02, 400.0, 40.0, 0.02, 0.01, 0.01, 0.40, 0.38, 0.30, 0.38)]
MZ, MI = arb(11.89), arb(14.57)
K1UP, K0UP = arb(-0.5853), arb(0.4214)       # upper ends of the certified enclosures (lemG/kcert2.log)
def upper_min(*xs):        # a ball whose upper end is the least upper end (valid upper bound of the min of the true values)
    return min(xs, key=lambda a: float(a.upper()))
def upper_max(*xs):
    return max(xs, key=lambda a: float(a.upper()))
def lower_min(*xs):        # for quantities that must be bounded below: smallest lower end
    return min(xs, key=lambda a: float(a.lower()))

# ------------------------- Lemma G remainder (remainder.bounds) in arb -------------------------
def remainder(u, tau, x, L):
    lU, lV, lX = u*L, tau*u*L, x*L; Dd = lV - lU; Y1, Y2 = lX - lV, lX - lU
    pi3 = (2*PI)**3; z1c = RC['z1c']
    ex = lambda s, l: (s*l).exp()
    Gw = RC['gam_c_w1']*ex(c, lX)
    Zb = upper_max(MZ, 1/(th - c + cp) + 1/(1 + c - th - cp))
    I1 = (ex(-th, lV) + ex(-th, lU))/Dd*RC['J_mth_0_T2_w1']
    I2 = (ex(-cp, lV) + ex(-cp, lU))/Dd*RC['J_mcp_0_inf_w1']
    A1 = RC['A1']
    out = {}
    out['E1a'] = Gw*z1c*I1*I2*Zb*MI*RC['iz1c_cp']*A1/pi3
    gauss = (-(eta*eta*T2*T2)/2).exp()
    ln = th + s0; sm = 1 + c - th - cp; R0 = (ln*ln + 1).sqrt()
    Hinf = (ex(-cp, lV) + ex(-cp, lU))*(eta*eta*cp*cp/2).exp()/(Dd*cp*cp)
    J0 = (ex(-cp, lV) + ex(-cp, lU))/Dd*RC['J_mcp_0_inf']; J1 = I2
    zconn = ln*((2 + (2 + W0 + T2)/sm)*J0 + J1/sm) + 2*PI*R0*Hinf
    hs1 = (ex(s0, lV) + ex(s0, lU))*(eta*eta*(th*th - T2*T2)/2).exp()/(Dd*T2*T2)
    E1b = 2*RC['gam_c']*ex(c, lX)*z1c*hs1*zconn*MI*RC['iz1c_cp']*RC['A1box']/pi3
    hsT = (ex(s0, lV) + ex(s0, lU))/Dd*RC['J_s0_T2_inf_w1']
    E1c = Gw*z1c*hsT*RC['z1c_s0_cp']*RC['iz1c_s0']*I2*RC['iz1c_cp']*RC['A1box']/pi3
    out['E1bc'] = E1b + E1c
    Ax = RC['Ax']
    Gx = RC['gam_c']*upper_max(ex(c, Y1), ex(c, Y2))
    ext = (ex(-(c - cp), lU) + ex(-(c - cp), lV))*(eta*eta*(c - cp)**2/2).exp()*gauss/(Dd*T2*T2)
    out['Rext'] = RC['gam_c']*ex(c, lX)*z1c*Ax*RC['iz1cp']*RC['iz1c_cp']*ext*J0/(2*PI)**2   # |X^w| = e^{c lX}
    Kt = 2*(eta*eta*((c - cp)**2 + cp*cp)/2).exp()*gauss/(3*(T2 - W0)**3)   # both sides; |w+s| >= t - W0
    out['Ktail'] = 4*RC['gam_c']*ex(c, Y2)*z1c*Ax*RC['iz1cp']*RC['iz1c_cp']*Kt/(2*PI)**2/Dd**2   # four terms, each coefficient <= e^{c Y2}
    AV, AU = RC['AV'], RC['AU']
    segV = PI/(b*b*(b + c))*(eta*eta*((c + b)**2 + b*b)/2).exp()
    segU = PI/(b*b*(b - c))*(eta*eta*((c - b)**2 + b*b)/2).exp()
    rV = (-Dd*b).exp()*z1c*AV*MI*RC['iz1c_b']*segV/(2*PI)
    rU = (-Dd*b).exp()*z1c*AU*MI*RC['iz1b']*segU/(2*PI)
    gauss2 = (-(eta*eta*(T2 - W0)**2)/2).exp()
    den4 = T2**2*(T2 - W0)**2
    connV = 2*(b + cp)*(cp*Dd).exp()*(eta*eta*((c + b)**2 + b*b)/2).exp()*gauss*gauss2/den4*z1c*RC['AVbox']*MI*RC['iz1c_cp']/(2*PI)
    connU = 2*(b - cp)*(eta*eta*(b*b + b*b)/2).exp()*gauss*gauss2/den4*z1c*RC['AUbox']*RC['iz1cp']*MI/(2*PI)
    out['rhoV'] = RC['gam_c']*ex(c, Y1)*(rV + connV)/Dd**2
    out['rhoU'] = RC['gam_c']*ex(c, Y2)*(rU + connU)/Dd**2
    out['wshift'] = (ex(-thp, Y1) + ex(-thp, Y2))*RC['R6']/(2*PI)/Dd**2
    gt = 2*arb('3.4')*W0.sqrt()*(-PI*W0/2).exp()*(2/PI)/(1 - 1/(PI*W0))   # |Gamma(s+iv)| <= 3.4 v^{1/2} e^{-pi v/2}, v >= 40, -thp <= s <= c
    out['wtail'] = gt*ex(c, lX)*z1c**4/RC['z2_2c']/(2*PI) + gt*(thp + c)*2*ex(c, Y2)*RC['KcW0']/Dd**2
    return sum(out.values(), arb(0))

def _Phibar(x): return (x/arb(2).sqrt()).erfc()/2
def _phi(x): return (-(x*x)/2).exp()/(2*PI).sqrt()
def _varpi(x): return _phi(x) - x*_Phibar(x)                 # x >= 0
def _Psi1(y): return ((1 + y*y)*_Phibar(y) - y*_phi(y))/2    # y >= 0
def _bsup(Alo, Ahi):
    """upper bound for b(A) = sum_r varpi((log r - A)/eta)/r over A in [Alo, Ahi] (Lemma A.bsum)."""
    c3 = arb('0.8')*(eta*eta/2).exp()
    if Alo >= 0:
        if Alo >= 5: return eta/2 + c3*(-arb(Alo)).exp()
        return eta*(arb(1)/2 - _Psi1(arb(Ahi)/eta)) + c3*(-arb(Alo)).exp()      # eta Psi1(-Ahi/eta) = eta(1/2 - Psi1(Ahi/eta))
    if Ahi <= 0:
        y = arb(-Ahi)/eta
        return eta*_Psi1(y) + 2*_varpi(y)
    return eta/2 + c3
def St_block(lU, lV, lX, a_star, du=0.02):
    """Lemma A.St: S_t <= (eps/e)^eps X^eps (1/eps+1) (eta/Delta)^2 Sigma + Phibar(lU/eta)^2, Sigma by a block sum."""
    lUf, lVf = float(lU.mid()), float(lV.mid())
    u0 = math.log(1e6); c3 = arb('0.8')*(eta*eta/2).exp()
    beta0 = eta + c3*((-(lU - arb(u0))).exp() + (-(lV - arb(u0))).exp())
    Sig = (arb(10**6).log() + arb('1.3326') + arb(4)/1000)*beta0**2         # m <= 1e6 (c0 = 1.33258... <= 1.3326)
    u = u0; umax = lVf + 60*0.02
    while u < umax:
        ua, ub = u, u + du
        bet = arb(0)
        for K in (lU, lV):
            Kf = float(K.mid()); rk = float(K.rad()) + 1e-12
            bet += _bsup(Kf - ub - rk, Kf - ua + rk)
        Sig += bet**2*(arb(du) + 8*(-arb(ua)/2).exp())                     # block mass <= du + 8 e^{-u/2}
        u = ub
    y0 = arb(umax - lVf)/eta                                              # tail u > umax: |A| >= y0 eta for both K
    Sig += 36*(eta + 2)**2*(-(y0*y0)).exp()/(2*PI)*(1 + eta/(2*y0))
    Dd = lV - lU
    a = arb(a_star); eps = a/lX
    val = (eps/arb(1).exp())**eps*a.exp()*(1/eps + 1)*(eta/Dd)**2*Sig
    return val + _Phibar(lU/eta)**2

def St_argmin(lU, lV, lX):
    import stbound as SB
    best = None
    for k in range(1, 40):
        a = 0.25*k; eps = a/lX
        v = (eps/math.e)**eps*math.exp(a)*(1/eps + 1)
        if best is None or v < best[0]: best = (v, a)
    return best[1]

def d4G(u, t, x, L):
    lU, lV, lX = A(u)*L, A(t)*A(u)*L, A(x)*L
    assert float(lU.lower()) >= 40 and lX > lV
    Y1, Y2, Dd = lX - lV, lX - lU, lV - lU
    D4m = ((Y2*Y2 - Y1*Y1)/2 - GAM*Dd + K1UP*(Y1 + Y2 - 2*GAM) + 2*K0UP)/(Dd*Dd)
    rem = remainder(A(u), A(t), A(x), L)
    Ss = (-1/lX.exp()).exp() + D4m + rem
    St = St_block(lU, lV, lX, St_argmin(float(lU.mid()), float(lV.mid()), float(lX.mid())))
    return D4m + rem + 2*(Ss*St).sqrt() + St + (-1/(A(x)*L).exp()).exp() - (-(L*L)/(A(u)*L).exp()).exp()

# ------------------------- pair coefficient (Lemma A''/C, lemmaD._lp2/_ntil/pair_coef) -------------------------
_MZL = {}
def mzl(e):          # -zeta'/zeta(1+e), exact ball
    if e not in _MZL:
        s = arb_series([1 + arb(e), 1], prec=2).zeta(); _MZL[e] = -s[1]/s[0]
    return _MZL[e]
def M_arb(logT, e):
    Fz = mzl(e); lT = arb(logT)
    m0 = (lT - (2*PI).log())/2 + Fz
    # float argmin over the alternatives
    Ff = float(Fz.mid()); best = (0.5*(logT - math.log(2*math.pi)) + Ff, None)
    for k in range(50, 201, 2):
        cc = k/100; v = 0.5*(logT - math.log(2*math.pi)) + (cc*cc + 0.5)/(2*cc)*Ff + (logT - math.log(math.pi))/(8*cc)
        if v < best[0]: best = (v, cc)
    if best[1] is None: return m0
    cc = arb(best[1])
    alt = (lT - (2*PI).log())/2 + (cc*cc + arb('0.5'))/(2*cc)*Fz + (lT - PI.log())/(8*cc)
    return upper_min(m0, alt)
def cost_lo(d, e, de):   # lower bound of min_{a in [e, e+de]} a/(a^2+d^2)
    a1, a2 = arb(e), arb(e) + de; d = arb(d)
    return lower_min(a1/(a1*a1 + d*d), a2/(a2*a2 + d*d))
_GR = {}
def gr(d, de):
    key = (d, de.mid().str(20))
    if key not in _GR:
        if d <= 0: _GR[key] = arb(1)
        else:
            v = acb(1 + 2*de, arb(d)).gamma().abs_upper()/(1 + 2*de).gamma()
            _GR[key] = upper_min(arb(1), arb(v))
    return _GR[key]
def lp_arb(h, de, logT, e, dst):
    """weak-duality upper bound for the fractional knapsack of lemmaD._lp2 (budget M, items (b_i, c_i))."""
    M = M_arb(logT, e) - 1/(arb(e) + de)
    items = []
    for k in range(1, int(1/h) + 2):
        if (k + 1)*h <= dst: continue
        dmin = max((k - 1)*h, dst)
        # b <= |Gamma| ratio / max(2 delta, dmin); both 2 delta and dmin are lower bounds for the max
        bb = gr(dmin, de)/(arb(dmin) if dmin >= float((2*de).upper()) else 2*de)
        cc = cost_lo((k + 1)*h, e, de)
        items += [(bb, cc)]*2
    # lambda from the float greedy
    Mf = float(M.mid()); fl = sorted(((float(x.mid())/float(y.mid()), x, y) for x, y in items), key=lambda z: -z[0])
    lam = 0.0; Mr = Mf
    for r, x, y in fl:
        if Mr <= 0: break
        lam = r; Mr -= float(y.mid())
    if Mr > 0: lam = 0.0
    lamA = arb(lam)
    ub = lamA*upper_max(M, arb(0)) if lam > 0 else arb(0)
    for x, y in items:
        z = x - lamA*y
        if float(z.upper()) > 0: ub += arb(z.upper())
    return ub
def lp_float(h, de, logT, e, dst):
    M = D._M(logT, e) - D._cost(0, e, de)
    items = []
    for k in range(1, int(1/h) + 2):
        if (k + 1)*h <= dst: continue
        dmin = max((k - 1)*h, dst); bb = D._gr(round(dmin, 9), de)/max(2*de, dmin)
        cc = D._cost((k + 1)*h, e, de); items += [(bb/cc, bb, cc)]*2
    items.sort(reverse=True); tot = 0.0
    for _, bb, cc in items:
        if M <= 0: break
        f = min(1.0, M/cc); tot += f*bb; M -= f*cc
    return tot
def pair_arb(k, L):
    """returns (coef, ntil) as (arb, int) for the cell delta_k = k/1e4, height L (= lk)."""
    def_ = k/1e4; de = arb(k)/10**4
    dst = D._dstar(L, def_)
    best = None
    for kp in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0):
        h = kp*def_
        # ntil: float argmin over eta, then arb
        en = min((D._M(L, j*0.0005)/D._cost(h/2, j*0.0005, def_), j*0.0005) for j in range(1, 800))[1]
        el = min((lp_float(h, def_, L, j*0.002, dst), j*0.002) for j in range(1, 200))
        nt_f = max(1, math.floor(min(D._M(L, j*0.0005)/D._cost(h/2, j*0.0005, def_) for j in range(1, 800))))
        val = nt_f*1.0067*(1 + 2*def_*el[0])
        if best is None or val < best[0]: best = (val, h, en, el[1])
    _, h, en, el = best
    nt = max(1, math.floor(fup(M_arb(L, en)/cost_lo(h/2, en, de))))
    lp = lp_arb(h, de, L, el, dst)
    return arb('1.0067')*(1 + 2*de*lp), nt

# ------------------------- I' (Lemma C) -------------------------
def nwin_arb(w, lt, de_num):           # de = de_num/1e4 exactly
    de = arb(de_num)/10**4; def_ = de_num/1e4; d = w/2
    ff = lambda a: a/(a*a + d*d)
    ef = min((math.floor((0.5*(lt - math.log(2*math.pi)) + 1/(j*0.0025) - 0.36)/min(ff(j*0.0025), ff(j*0.0025 + def_))), j*0.0025)
             for j in range(1, 600))[1]
    M = (arb(lt) - (2*PI).log())/2 + mzl(ef)
    dA = arb(d); f1 = arb(ef)/(arb(ef)**2 + dA*dA); a2 = arb(ef) + de; f2 = a2/(a2*a2 + dA*dA)
    n1 = math.floor(fup(M/lower_min(f1, f2)))
    n2 = math.ceil(Fr(1, 4)/Fr(de_num, 10**4)) + 1
    return min(n1, n2)
def g_arb(y, de): y = arb(y); return (y/6).exp()*y**(2*de - arb('0.5'))*(-PI*y/2).exp()
def Iprime_arb(k, lt, w=0.25):
    worst = None
    dk = arb(k)/10**4
    for dn in (k - 1, k):
        n = nwin_arb(w, float(lt), dn); acc = arb(0); y = 1.0
        while y < 60:
            acc += n*g_arb(y, dk); y += w                  # g decreasing on [1, oo): window max at the left end
        acc += n*(arb('-1.4')*60).exp()/(1 - (arb('-0.35')).exp())   # windows beyond 60: g(y) <= e^{-1.4 y}
        v = arb(dn)/10**4*acc
        worst = v if worst is None else upper_max(worst, v)
    return worst

# ------------------------- consts() -------------------------
def consts_arb(p, de, L):
    u, t, w, x = [A(float(v)) for v in p]; v = t*u; al = 1 - de
    assert 0.02 < float(w.mid()) < float(u.mid()) and float(t.mid()) > 1.05 and x*L > u*L - 2*L.log()
    IG = arb('5.922')
    B = arb('10.45')*arb('0.618')*L*IG/(2*PI)*(L*((v + w)/2 + arb(1)/6 + x*(arb('0.5') - al))).exp()
    tail = arb('3.8')*(-L*(L - 2*u)).exp()*(L*L + 1)/(w*(v - u)*L**6)
    c1 = (-(L*L)/(u*L).exp()).exp() - B - tail
    d4 = d4G(float(u.mid()), float(t.mid()), float(x.mid()), L)
    d5 = arb('1.0061')*(1 + (arb('1.333') + 4*(-w*L/2).exp())/(w*L))
    J = arb('4.72')*arb('0.618')*(L/6).exp()*L
    c5 = d4*d5/(2*PI)*((w*L).exp()/(w*L))*((L*x*(arb('1.5') - 2*al)).exp() + ((L*u - 2*L.log())*(arb('1.5') - 2*al)).exp())*J
    assert c1 > 0 and c1*c1 - c5 > arb('0.2'), (c1, c5)
    return d4*d5/w, c1*c1 - c5, d4

def cell(k, lk, p):
    de = arb(k)/10**4; L = arb(lk)
    dd, den, d4 = consts_arb(p, de, L)
    coef, nt = pair_arb(k, float(lk))
    coef = coef + 2*(2*PI).sqrt()*arb('2.4')*Iprime_arb(k, lk)
    val = (dd*coef*nt/den).log() + 2*arb(float(p[3]))*de*L
    return val - (de*L).log(), val

if __name__ == "__main__":
    tab = {}
    for f in sorted(glob.glob(os.path.join(CODE, 'lemmaD_tableGH0.01_*.pkl'))): tab.update(pickle.load(open(f, 'rb')))
    keys = sorted(k for k in tab if tab[k] is not None)
    part, nparts = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (0, 1)
    outf = os.path.join(HERE, 'cells_arb_%d.json' % part)
    try: out = json.load(open(outf))
    except Exception: out = {}
    t0 = time.time(); n = 0
    for i, key in enumerate(keys):
        if i % nparts != part: continue
        k, lk = key[0], key[1]; nm = '%d,%s' % (k, lk)
        if nm in out: continue
        fv, p = tab[key]
        if fv > 1e8: out[nm] = None; continue
        lo, val = cell(k, lk, p)
        out[nm] = {'logN_up': repr(fup(lo)), 'val_arb_up': repr(fup(val)), 'val_float': repr(float(fv))}
        n += 1
        if n % 10 == 0:
            json.dump(out, open(outf + '.tmp', 'w')); os.replace(outf + '.tmp', outf)
            print('%d cells, %.0fs; last %s: arb %.6f float %.6f' % (n, time.time() - t0, nm, fup(val), fv), flush=True)
    json.dump(out, open(outf + '.tmp', 'w')); os.replace(outf + '.tmp', outf)
    print('done', n, '%.0fs' % (time.time() - t0))
