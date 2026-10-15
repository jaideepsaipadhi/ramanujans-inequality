"""Lemma C: crowding bound for Bellotti's far off-diagonal term c4 (|gamma_r - gamma_s| > 1).
Bellotti bounds the representatives with |gamma_r-gamma_s| in [y, y+dy] by dy/delta (one per delta-box), giving
  c4 = 2 sqrt(2 pi) 2.4 (d4 d5 / w) * I,   I = int_1^T g(y) dy,   g(y) = e^{y/6} y^{1.5-2 alpha} e^{-pi y/2},
inside the common factor X^{2-2alpha}/(delta log T).  Representatives are zeros with beta >= 1 - delta.  For any
window of height w centred at t, with s = 1 + eta + i t:
  sum_rho Re 1/(s - rho) = Re zeta'/zeta(s) + Re 1/(s-1) + (1/2) Re psi(s/2+1) - (1/2) log pi
                         <= 1/eta + (1/2) log T          (all terms >= 0; -zeta'/zeta(1+eta) < 1/eta)
and each such zero in the window contributes >= min(f(eta), f(eta+delta)), f(a) = a/(a^2 + (w/2)^2).  So a window
holds n(w) <= min(floor(M / min f), ceil(w/delta) + 1) representatives, and
  I  ->  I' = delta * sum_i n(w) * max_{[y_i, y_i + w]} g      (windows [y_i, y_i+w] tiling [1, 60]; g < 1e-38 beyond).
Returned: c4' for given delta range [de_lo, de_hi] and log T (worst case over the range)."""
import math, functools
def g(y, de): return math.exp(y/6) * y**(2*de - 0.5) * math.exp(-math.pi*y/2)
def nwin(w, lt, de):
    best = 10**9
    for k in range(1, 600):
        eta = k*0.0025; M = 0.5*(lt - math.log(2*math.pi)) + 1/eta - 0.36; d = w/2
        f = lambda a: a/(a*a + d*d)
        best = min(best, math.floor(M/min(f(eta), f(eta+de))))
    return min(best, math.ceil(w/de) + 1)
@functools.lru_cache(maxsize=None)
def Iprime(de_key, lt_key, w=0.25):
    lt = float(lt_key)
    worst = 0.0
    for de in (de_key/1e4 - 1e-4, de_key/1e4):        # both ends of the delta bin
        n = nwin(w, lt, de); acc = 0.0; y = 1.0
        while y < 60:
            acc += n*max(g(y + w*i/40, de_key/1e4) for i in range(41)); y += w   # g at the bin's largest delta
        worst = max(worst, de*acc)
    return worst
def c4prime(de_key, lt_key, d4d5w):
    return 2*math.sqrt(2*math.pi)*2.4*d4d5w*Iprime(de_key, lt_key)
if __name__ == "__main__":
    for de in [0.0055, 0.007, 0.0095]:
        for lt in [33, 35]:
            print(de, lt, 'c4 =', round(c4prime(math.ceil(de*1e4), lt, 264.95/1.0067), 1), '(Bellotti-form 447.23)')
