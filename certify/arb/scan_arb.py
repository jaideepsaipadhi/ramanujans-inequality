"""Ball-arithmetic re-evaluation of the scan (Proposition master + Lemma Stieltjes) on every interval [L1, L2].

For each interval the kernel parameters (c, f) and eps are chosen by the original floating-point search (code/interval.py);
any choice is admissible, so the choice need not be rigorous.  With these parameters fixed, every term of
Proposition master is evaluated in arb:
  * tol = eps_R(L2 + 1)/(1 + 2/L1) with eps_R from the rigorous expansion of certify/epsR_cert.py (lower end);
  * the smoothing term with the actual kappa(eps) (no fixed-point assumption), the truncation term, the H-term;
  * S(s) by the Stieltjes bound with the factor e^{eps(c+1)/6}; partition points in log t are doubles, the last one >= log T1;
  * N^+(s, t) = min of: the zero count N(T), Kadiri-Lumley-Ng, Bellotti (uniform, T <= e^29, Table rows), and the Lemma D
    cell value from cells_arb.json (certified);  cell indices k = ceil((1 - s) 1e4), lk computed in exact rational arithmetic;
  * blocks are dropped only when s >= the zero-free edge at the top of the block holds in ball arithmetic;
  * the integral over s uses S(a) times the exact integral of e^{(s-1) L1} on each cell [a, b] (S non-increasing).
Output: one line 'L1 L2 q_upper' per interval (q_upper an upper bound for the ratio of Proposition master to tol)."""
import sys, os, json, math, bisect, time, glob
from fractions import Fraction as Fr
from flint import arb, ctx
HERE = os.path.dirname(os.path.abspath(__file__)); CODE = os.path.join(HERE, '..', '..', 'code')
sys.path.insert(0, CODE); sys.path.insert(0, os.path.join(HERE, '..'))
os.environ['LEMMAG'] = '1'
import interval as IV, ramanujan_upper as RU, epsR_cert as ER
ctx.prec = 100
PI = arb.pi()
def fup(a):
    return float(a.upper())
def fup_safe(a):          # a double >= every point of the ball
    import numpy as np
    return float(np.nextafter(np.nextafter(float(a.upper()), np.inf), np.inf))
def upper_min(*xs): return min(xs, key=lambda a: float(a.upper()))

# --------- certified Lemma D cells ---------
CELL = {}
for f in glob.glob(os.path.join(HERE, 'cells_arb_*.json')):
    for k, v in json.load(open(f)).items():
        if v is not None: CELL[k] = arb(v['logN_up'])
# --------- zero-free region (lambda* table, rep/lamstar2c.json) ---------
LAM = RU._LAM; LAMX = [x for x, _ in LAM]; R0 = arb('4.896')
def lam_at(lt):           # lambda used on (L_i, L_{i+1}]: value at the grid point >= lt  (lt a double)
    if not LAM or lt > LAMX[-1]: return 1/R0
    i = bisect.bisect_left(LAMX, lt)
    v = arb(LAM[i][1]); return v if v > 1/R0 else 1/R0
def edge(lt):             # zero-free edge 1 - lambda(lt)/lt at height e^lt (lt a double); non-decreasing in lt (lambda/l non-increasing)
    return 1 - lam_at(lt)/arb(lt)

KLN = [(Fr(str(a)), arb(str(A)), arb(str(B))) for a, A, B in RU.KLN]
H0K = arb('3.0610046e10')
ROWS = [(a, b, Fr(str(al)), arb(str(C)), arb(str(B))) for a, b, al, C, B in __import__('bell_patch').ROWS]
_LT = {}
def _lt(lT):
    v = _LT.get(lT)
    if v is None:
        T = arb(lT).exp(); lTa = arb(lT)
        v = {'T': T, 'lTa': lTa, 'all': (T/(2*PI)*(T/(2*PI*arb(1).exp())).log() + lTa + 1).log(), 'kln': {}}
        if len(_LT) > 20000: _LT.clear()
        _LT[lT] = v
    return v
KLN_S = [r[0] for r in KLN]
def logN(s, lT):
    """upper bound (arb) for log N(s, e^lT); s and lT are doubles (exact)."""
    sF = Fr(s); P = _lt(lT); T, lTa = P['T'], P['lTa']
    cands = [P['all']]
    i = bisect.bisect_right(KLN_S, sF) - 1                  # largest tabulated sigma0 <= s
    if i >= 0 and T > H0K:
        if i not in P['kln']:
            s0, A, B = KLN[i]; s0a = arb(s0.numerator)/s0.denominator
            P['kln'][i] = (A*lTa**(2*s0a)*lTa**(5 - 4*s0a)*T**(arb(8)/3*(1 - s0a)) + B*lTa**2).log()
        cands.append(P['kln'][i])
    one_s = arb(1) - arb(s)
    if sF >= Fr('0.9927') and lT > math.log(3e12):
        cands.append(LOG_C1 + arb('1.448')*one_s*lTa)
        if lT <= 29: cands.append(LOG_C2 + arb('14.160')*one_s*lTa)
    for j, (a, b, al, C, B) in enumerate(ROWS):
        lo = math.log(a) if j == 0 else a
        if lo < lT <= b and sF >= al:
            cands.append(ROWLOG[j] + B*one_s*lTa)
    de = 1 - sF
    if de <= Fr(3, 100) and 28 < lT <= 40:
        k = math.ceil(de*10**4)
        lTF = Fr(lT)
        lk = Fr(math.ceil(2*lTF), 2) if lT <= 31 else Fr(math.ceil(lTF))
        key = '%d,%s' % (k, (int(lk) if lk.denominator == 1 else float(lk)))
        if key in CELL: cands.append(CELL[key])
    return upper_min(*cands)
LOG_C1, LOG_C2 = arb('1.62e11').log(), arb('37341.72').log()
ROWLOG = [r[3].log() for r in ROWS]

def wfun(u, c):
    """Logan weight bound w(u, c) = c sinh(sqrt(c^2-u^2))/(sqrt(c^2-u^2) sinh c) (u < c), c/sinh c (u >= c); u an arb ball."""
    s2 = c*c - u*u
    if s2 <= 0: return c/c.sinh()
    if s2 > 0:
        s = s2.sqrt(); return c*s.sinh()/(s*c.sinh())
    # straddling: sinh(s)/s is increasing in s >= 0, between 1 and its value at sqrt(s2.upper)
    su = arb(s2.upper()).sqrt(); top = su.sinh()/su
    return c/c.sinh()*arb((1 + fup(top))/2, (fup(top) - 1)/2 + 1e-30)

def zsum(s, H, eps, c, lT1_up, nt=80):
    """Stieltjes upper bound for S(s) (times e^{eps(c+1)/6}); eps, c arbs; partition in log t: doubles from log H to lT1_up."""
    lH = math.log(H)
    us = [lH + (lT1_up - lH)*i/nt for i in range(nt + 1)]; us[0] = lH; us[-1] = lT1_up
    if lT1_up <= lH: return arb(0)
    f = [wfun(arb(u).exp()*eps, c)/arb(u).exp() for u in us]
    acc = arb(0); sa = arb(s)
    for i in range(nt):
        if sa >= edge(us[i+1]): continue          # certified zero-free up to t_{i+1}
        acc += logN(s, us[i+1]).exp()*(f[i] - f[i+1])
    if not (sa >= edge(us[-1])): acc += logN(s, us[-1]).exp()*f[-1]
    return acc*(eps*(c + 1)/6).exp()

def q_arb(L1, L2, H, cc, frac):
    """upper bound for q on [L1, L2] with kernel (cc, frac) (float search output)."""
    # eps exactly as interval.one() computes it (a double: just a parameter)
    tol_f = RU.eps_max(L2 + 1)/(1 + 2/L1)
    eps_f = frac*tol_f*math.sqrt(cc)/(0.7946*1.02)
    for _ in range(3):
        lam = math.log(10*math.sqrt(cc)/eps_f); kap = (L2 + lam)/(L1 - lam)
        eps_f = frac*tol_f*math.sqrt(cc)/(0.7946*1.02*kap)
    eps = arb(eps_f); c = arb(cc); L1a, L2a = arb(L1), arb(L2)
    assert 30 <= cc <= 110 and eps_f <= 1e-4
    epsR, g = ER.eps_ball(L2 + 1, L2 + 1)
    tol = arb(epsR.lower())/(1 + 2/L1a)
    lam = (10*c.sqrt()/eps).log()
    assert lam < L1a - 1
    kap = (L2a + lam)/(L1a - lam)
    smooth = kap*arb('0.7946')*arb('1.02')*eps/c.sqrt()
    T1 = c/eps; lT1 = T1.log()
    tr = arb('0.32')/c.sinh()*(arb('0.71')*(c*eps).sqrt()).exp()*(3*c).log()*lT1
    tail = arb(0)
    if T1 > H or not (T1 < H):
        lT1_up = fup_safe(lT1)
        smax = edge(lT1_up)                         # zero-free edge at the top: >= the edge at every height <= T1
        sm_up = fup_safe(smax)
        s_mid = sm_up - 0.02
        grid = [0.6 + (s_mid - 0.6)*k/60 for k in range(61)]
        hf = 1/(40*L1); k = 1
        while s_mid + k*hf < sm_up: grid.append(s_mid + k*hf); k += 1
        grid.append(sm_up)
        acc = arb(0)
        for a, b in zip(grid[:-1], grid[1:]):
            acc += zsum(a, H, eps, c, lT1_up)*(((arb(b) - 1)*L1a).exp() - ((arb(a) - 1)*L1a).exp())/L1a
        tail = 2*(L2a*acc + (-arb('0.4')*L1a).exp()*(zsum(0.6, H, eps, c, lT1_up) + lT1**2))
    Hterm = (-L1a/2).exp()*arb(H).log()**2
    e = smooth + tr + tail + Hterm
    return e/tol, {'smooth': fup(smooth/tol), 'tail': fup(tail/tol), 'tr': fup(tr/tol), 'H': fup(Hterm/tol)}

def choose(L1, L2, H):
    """float search of interval.err_interval, returning the argmin (c, frac)."""
    one = lambda c, f: IV_one(L1, L2, H, c, f)
    v, c0, f0 = min((one(c, f), c, f) for c in range(40, 101, 10) for f in (0.5, 0.6, 0.7, 0.8, 0.9))
    best = (v, c0, f0)
    for c in range(max(30, c0 - 10), c0 + 11, 5):
        for f in [f0 + 0.05*k for k in range(-2, 3)]:
            if 0.3 <= f <= 0.98:
                w = one(c, f)
                if w < best[0]: best = (w, c, f)
    return best

# a copy of interval.err_interval's inner function (same arithmetic), exposed for the argmin
def IV_one(L1, L2, H, c, frac):
    r = RU; logan = IV.logan
    tol = r.eps_max(L2 + 1)/(1 + 2/L1)
    eps = frac*tol*math.sqrt(c)/(0.7946*1.02)
    for _ in range(3):
        lam = math.log(10*math.sqrt(c)/eps); kap = (L2 + lam)/(L1 - lam)
        eps = frac*tol*math.sqrt(c)/(0.7946*1.02*kap)
    T1 = c/eps
    tr = 0.32/math.sinh(c)*math.exp(0.71*math.sqrt(c*eps))*math.log(3*c)*math.log(T1)
    tail = 0.0
    if T1 > H:
        smax = 1 - 1/(r.Rof(math.log(T1))*math.log(T1)); smax_int = smax; s_mid = smax_int - 0.02
        grid = [0.6 + (s_mid - 0.6)*k/60 for k in range(61)]
        hf = 1/(40*L1); k = 1
        while s_mid + k*hf < smax_int: grid.append(s_mid + k*hf); k += 1
        grid.append(smax_int)
        acc = 0.0
        for a, b in zip(grid[:-1], grid[1:]):
            acc += logan.zsum(a, H, eps, c)*(math.exp((b - 1)*L1) - math.exp((a - 1)*L1))/L1
        tail = 2*(L2*acc + math.exp(-0.4*L1)*(logan.zsum(0.6, H, eps, c) + math.log(T1)**2))
    e = frac*tol + tr + tail + math.exp(-L1/2)*math.log(H)**2
    return e/tol

def choose_fast(L1, L2, H):
    """same search and s-grid as choose(), with the Stieltjes partition coarsened to nt = 20 (only used to pick (c, f))."""
    import logan
    orig = logan.zsum
    logan.zsum = lambda s_, H_, e_, c_, nt=80: orig(s_, H_, e_, c_, nt=20)
    try: return choose(L1, L2, H)
    finally: logan.zsum = orig

RANGES = [(43, 57, 2), (57, 1000, 10), (1000, 2000, 10), (2000, 2600, 10), (2600, 3159, 10)]
if __name__ == "__main__":
    H = 3.0e12
    part, nparts = int(sys.argv[1]), int(sys.argv[2])
    tops = []
    for lo, hi, w in RANGES:
        L = lo
        while L < hi:
            tops.append((L, min(L + w, hi))); L = min(L + w, hi)
    # the binding stretch first (load balance): sort by cost proxy
    outf = os.path.join(HERE, 'scan_arb_%d.txt' % part)
    done = set()
    if os.path.exists(outf):
        for line in open(outf):
            p = line.split()
            if len(p) >= 3: done.add((float(p[0]), float(p[1])))
    t0 = time.time()
    def leaf(a, b, fo):
        if (a, b) in done: return
        qf, cc, fr = (choose if a < 1100 else choose_fast)(a, b, H)
        if qf >= 1 and b - a > 0.5:
            m = (a + b)/2; leaf(a, m, fo); leaf(m, b, fo); return
        q, parts = q_arb(a, b, H, cc, fr)
        if fup(q) >= 1 and b - a > 0.5 and a >= 57:
            m = (a + b)/2; leaf(a, m, fo); leaf(m, b, fo); return
        fo.write('%r %r %r  c=%d f=%.2f float=%.6f %s\n' % (a, b, fup_safe(q), cc, fr, qf, json.dumps(parts))); fo.flush()
        print('%g %g q_arb=%.6f float=%.6f  %.0fs' % (a, b, fup(q), qf, time.time() - t0), flush=True)
    with open(outf, 'a') as fo:
        for i, (a, b) in enumerate(tops):
            if i % nparts == part: leaf(a, b, fo)
    print('done %.0fs' % (time.time() - t0))
