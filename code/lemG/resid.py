# Lemma G step 1: residue main term R(w) for F_ij(w) = sum_n Lambda_{zi}(n) Lambda_{zj}(n) n^{-1-w}  (zi <= zj),
# checked against the exact finite form zeta(1+w) sum_{d,e} mu mu log+(zi/d) log+(zj/e) / [d,e]^{1+w}.
import mpmath as mp, numpy as np, sys
mp.mp.dps = 30
PR = [int(p) for p in mp.primerange(2, 20000)] if hasattr(mp,'primerange') else None
def primes(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return [int(p) for p in np.nonzero(s)[0]]
PR = primes(20000)
def logA(w, s1, s2):
    t = mp.mpf(0)
    for p in PR:
        a = mp.mpf(p)**(-1-w-s1); b = mp.mpf(p)**(-1-w-s2); c = mp.mpf(p)**(1+w)
        t += mp.log((1-a-b+a*b*c)*(1-a*b*c)/((1-a)*(1-b)))
    return t
def A(w, s1, s2): return mp.e**logA(w, s1, s2)
def zl(s): return mp.zeta(s, derivative=1)/mp.zeta(s)
def R(w, zi, zj):
    li, lj = mp.log(zi), mp.log(zj)
    A12 = mp.diff(lambda x, y: A(w, x, y), (0, 0), (1, 1))
    dzl = mp.diff(lambda s: zl(s), 1+w)
    main = li*lj + dzl + A12 - zj**(-w)/w**2
    Rx = lambda s2: zi**(-w-s2)*zj**s2*mp.zeta(1+w)*A(w, -w-s2, s2)/(mp.zeta(1-s2)*mp.zeta(1+w+s2)*(w+s2)**2*s2**2)
    h = mp.mpf('1e-8')
    r0 = (Rx(h)*h + Rx(-h)*(-h))/2                     # simple pole at 0: lim s2 Rx
    r1 = (Rx(-w+h)*h + Rx(-w-h)*(-h))/2               # simple pole at -w
    return main + r0 + r1, main, r0, r1
def exact(w, zi, zj):
    zi_, zj_ = int(zi), int(zj); N = max(zi_, zj_)
    mu = np.ones(N+1, int); isp = np.ones(N+1, bool)
    for p in range(2, N+1):
        if isp[p]:
            isp[2*p::p] = False; mu[p::p] *= -1; mu[p*p::p*p] = 0
    ds = [d for d in range(1, N+1) if mu[d]]
    s = mp.mpf(0)
    for d in ds:
        if d > zi: break
        ld = mp.log(mp.mpf(zi)/d)
        for e in ds:
            if e > zj: break
            l = d*e//np.gcd(d, e)
            s += mu[d]*mu[e]*ld*mp.log(mp.mpf(zj)/e)/mp.mpf(l)**(1+w)
    return mp.zeta(1+w)*s
if __name__ == "__main__":
    for zi, zj in ((60, 60), (60, 200), (150, 400)):
        for w in (mp.mpf('0.3'), mp.mpf('0.6')):
            r, m, r0, r1 = R(w, zi, zj); ex = exact(w, zi, zj)
            print(f"zi={zi} zj={zj} w={w}: exact {float(ex):.6f}  residue-main {float(r):.6f}  diff {float(ex-r):+.2e}   (parts {float(m):.3f} {float(r0):.3f} {float(r1):.3f})", flush=True)
