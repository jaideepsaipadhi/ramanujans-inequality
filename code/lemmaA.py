"""Lemma A candidate: Bellotti's Halasz-Montgomery bound with
 (i) the c3 pair sum replaced by the LP bound from Re-sum_rho 1/(s-rho) <= 0.5 log T + 1/eta (pairsum.py),
 (ii) the local count n~(alpha)=b1*delta*L+b2 replaced by 1 when the box bound < 2.
N(sigma,T) <= (c2 + rho*c3 + c4)/(c1^2-c5) * T^(B*delta)/(delta*log T),  delta = 1-sigma.
c's from Bellotti Table 6, row [e^34,e^35], inflated 6% per unit of log T above 35 (conservative)."""
import math, functools, os, ramanujan_upper as r, bell_patch, pairsum, lemmaC
USE_C = not os.environ.get('NO_LEMMA_C')
ROWS = bell_patch.ROWS
@functools.lru_cache(maxsize=None)
def rho_box(de_key, lt_key):
    de, lt = de_key/1e4, float(lt_key)
    full = 2*sum(1/(n*de) for n in range(1, int(1/de)))  # Bellotti's assumed sum (his c3)
    best = min(pairsum.lp(de, lt, k*0.002) for k in range(1, 150))
    return best/full, pairsum.box(de, lt)
base = r.logN_density
def logN(s, T, bell_scale=1.0):
    v = base(s, T, bell_scale); lT = math.log(T); de = 1-s
    if not os.environ.get('LEGACY_A'): return v   # superseded by Lemma D's corrected pair-sum accounting
    if de > 0.0095 or lT < 29 or lT > 35: return v   # only where Table 6 c-values are verified monotone (c_i increase with T)
    row = None
    for i,(a,b,al,C,B) in enumerate(ROWS):
        lo = math.log(a) if i == 0 else a
        if lo < lT <= b: row = (al, B)
    if row is None or s < row[0]: return v
    rho, bx = rho_box(math.ceil(de*1e4), math.ceil(lT))   # round delta, logT up (conservative)
    if bx >= 2: return v
    infl = 1.06**max(0.0, lT-35)
    c2, c3, c4, c1, c5 = 264.95*infl, 3128.0*infl, 447.23*infl, 0.984, 0.066   # max(table, recomputed from Bellotti formulas)
    if USE_C: c4 = lemmaC.c4prime(math.ceil(de*1e4), math.ceil(lT), c2/1.0067)   # Lemma C (crowding), c2 already inflated
    K = (c2 + rho*c3 + c4)/(c1*c1 - c5)
    return min(v, math.log(K) + row[1]*de*lT - math.log(de*lT))
r.logN_density = logN
if __name__ == "__main__":
    for s in [0.991, 0.993, 0.995]:
        for T in [1e14, 1e15]:
            print(s, '%.0e'%T, 'old %.3e'%math.exp(base(s,T)), 'new %.3e'%math.exp(logN(s,T)))
# Lemma D (lemmaD.py): improved mollifier bound + per-(sigma, T) parameter choice; take the minimum
if not os.environ.get('NO_LEMMA_D'):
    import lemmaD
    _prevD = r.logN_density
    def logN_D(s, T, bell_scale=1.0):
        v = _prevD(s, T, bell_scale)
        d = lemmaD.logN(s, T)
        return v if d is None else min(v, d)
    r.logN_density = logN_D
