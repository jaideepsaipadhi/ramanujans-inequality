# Certified iteration for Lemma H.  Step j: hypothesis (1-beta) log g >= lam_{j-1}(grid >= log g) for g >= 3e12 (j=1: BTY's
# 1/4.896 everywhere).  Claim: (1-beta) log g >= lam_j(L) for log g in (L-0.25, L].  Proof: for every l in that cell and every
# lam0 in [lam_{j-1}(L), lam_j(L)], margin1 = rhs1 - lhs > 0 (ball arithmetic, adaptive boxes), with
# sigma = 1 - A/log(K e^l + T0), A = lam_{j-1}(grid >= log(K e^l + T0)) (Lemma 6 hypothesis; lam(l)/l decreasing).
import math, json, sys
from rep_cert import *
import rep2
Ls = [28 + 0.25*i for i in range(61)]
def mk_lookup(tab):
    def f(l):
        ups = [L for L in Ls if L >= l - 1e-12]; return tab[min(ups)] if ups else 1/4.896
    return f
def margin1(Lt, sig, e0, C):
    E2, E3, c10, G = C
    lhs = e0*W0*((a*kap/2)*Lt - G) - c10*e0*W0 + a0*F(acb(sig-1), e0, sig).real + arb('1e-10')*e0 + arb('2e-7') \
          + E2*(e0*Lt)*e0*e0 + E3*e0**3
    return a1*(F(acb(sig-1+e0), e0, sig).real + F(acb(sig-e0), e0, sig).real) - lhs
def cert_cell(L, lam_lo, lam_hi, Aof, maxdepth=40):
    l_lo = max(L - 0.25, 28.7) if L > 28.7 else L - 0.25
    eta0c = math.ceil(lam_hi/l_lo*1e4)/1e4 + 1e-4
    cache = {}
    def rec(l0, l1, x0, x1, depth):
        A = Aof(math.log(K*math.exp(l1) + T0))
        slo = 1 - A/math.log(K*math.exp(l0) + T0); shi = 1 - A/math.log(K*math.exp(l1) + T0)
        key = math.floor(slo*1e5)
        if key not in cache: cache[key] = consts(eta0c, key/1e5)
        Lt = ball(l0, l1); lam = ball(x0, x1); e0 = lam/Lt
        if margin1(Lt, ball(slo, shi), e0, cache[key]) > 0: return True
        if depth >= maxdepth: return False
        if (l1 - l0)/0.25 > (x1 - x0)/max(lam_hi - lam_lo, 1e-9):
            m = (l0 + l1)/2; return rec(l0, m, x0, x1, depth+1) and rec(m, l1, x0, x1, depth+1)
        m = (x0 + x1)/2; return rec(l0, l1, x0, m, depth+1) and rec(l0, l1, m, x1, depth+1)
    return rec(l_lo, L, lam_lo, lam_hi, 0)
def _one(args):
    L, lo, tabd = args
    Aof = mk_lookup(tabd)
    cand = rep2.lam_star(L, Aof(math.log(K*math.exp(L) + T0))) - 2e-5
    hi = max(lo, cand)
    while hi > lo and not cert_cell(L, lo, hi, Aof): hi -= 2e-4
    return L, hi
if __name__ == "__main__":
    from multiprocessing import Pool
    nsteps = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    tab = {L: 1/4.896 for L in Ls}
    with Pool(2) as P:
        for j in range(1, nsteps + 1):
            tab = dict(P.map(_one, [(L, tab[L], tab) for L in Ls]))
            print('step', j, [(L, round(tab[L], 5)) for L in Ls[::4]], flush=True)
            json.dump({str(L): v for L, v in tab.items()}, open('lamstar2c.json', 'w'))
