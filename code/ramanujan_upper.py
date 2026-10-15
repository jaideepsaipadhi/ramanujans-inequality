#!/usr/bin/env python3
"""
ramanujan_upper.py -- the upper threshold for Ramanujan's inequality with
REAL published constants plugged in.

Target (exact, from li):  |pi(y)-li(y)| <= eps_max(L) * y/log y  at y=x, x/e
with eps_max(L) ~ 1/(3L^4)  (computed exactly, not approximated).

Bound for (psi(x)-x)/x via a band-limited (Buthe/Logan-type) smoothing:
  * smoothing cost              ~ eps_s            (kernel width eps_s)
  * kernel truncation           ~ L * exp(-c)      (bandwidth T1 = c/eps_s)
  * zeros with gamma <= H       on the line: x^(-1/2) log^2 H / (2 pi)
  * zeros with H < gamma <= T1  possibly off-line:
        2 * sum_{blocks [T,2T]} (1/T) * [ x^(-0.4) N_all(2T)
                                        + x^(s0-1) N(s0,2T)
                                        + L * int_{s0}^{bmax} x^(s-1) N(s,2T) ds ]
    (integration by parts, s0 = 0.6, N monotone so tables are used stepwise)

Inputs (all published):
  zero-free region  beta <= 1 - 1/(R log gamma), R = 5.558691
                    (Mossinghoff-Trudgian-Yang)
  N(sigma,T), 0.60<=sigma<=0.99, T >= 3.06e10:
       A (log T)^(2 sigma) (log T)^(5-4 sigma) T^(8/3 (1-sigma)) + B (log T)^2
       Kadiri-Lumley-Ng 2101.12263, Table 1 (k=1 rows; k<1 rows treated as k=1,
       which only enlarges log(kT) -- conservative)
  N(sigma,T), sigma >= 0.9927, 3e12 < T:  C T^(B'(1-sigma))
       Bellotti 2405.12545: uniform C=1.62e11, B'=1.448; and for T<=e^29
       C=37341.72, B'=14.160.  Min of applicable bounds is used.

The smoothing-kernel constants are MODELLED (exact Buthe constants are not
in here): the kernel term is taken as L*exp(-c) and c is chosen so it costs
1% of the budget.  Treat the output as a calibrated estimate, not a proof.

Run: python3 ramanujan_upper.py
"""
import math
import mpmath as mp

mp.mp.dps = 60
R = 4.896   # Bellotti-Trudgian-Yang 2026 (arXiv:2603.21490); was 5.558691
# Height-dependent constant (rep/rep.py, Oct 14): BTY's Lemma 15 inequality evaluated EXACTLY (no Lemma 14 fit) at each height,
# iterated; zeta(s) != 0 for sigma > 1 - lam(log t)/log t.  lam on a grid in log t; for log t in (L_i, L_{i+1}] use lam(L_{i+1}).
import json as _json, os as _os, bisect as _bis
try:
    _LAM = sorted((float(k), v) for k, v in _json.load(open(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'rep', 'lamstar2c.json'))).items())
except Exception:
    _LAM = []
def Rof(lt):
    if not _LAM or lt > _LAM[-1][0] or _os.environ.get('NO_LAMH'): return R
    i = _bis.bisect_left([x for x, _ in _LAM], lt)
    return 1.0/max(_LAM[i][1], 1/R)
H0_KLN = 3.0610046e10

KLN = [  # sigma0, A, B
    (0.60, 2.177, 5.663), (0.65, 2.963, 5.249), (0.70, 3.983, 4.824),
    (0.75, 5.277, 4.403), (0.80, 6.918, 3.997), (0.85, 8.975, 3.588),
    (0.86, 9.441, 3.514), (0.87, 9.926, 3.430), (0.88, 10.431, 3.346),
    (0.89, 10.955, 3.262), (0.90, 11.499, 3.186), (0.91, 12.063, 3.102),
    (0.92, 12.646, 3.017), (0.93, 13.250, 2.941), (0.94, 13.872, 2.856),
    (0.95, 14.513, 2.772), (0.96, 15.173, 2.694), (0.97, 15.850, 2.609),
    (0.98, 16.544, 2.532), (0.99, 17.253, 2.446),
]
BELL_ALPHA = 0.9927


def eps_max(L):
    L = mp.mpf(L)
    x = mp.e ** L
    li_x, li_xe = mp.li(x), mp.li(x / mp.e)
    g = (li_x ** 2 - (mp.e * x / L) * li_xe) / (x / L) ** 2
    b = 2 * li_x / (x / L) + L / (L - 1)
    return float((-b + mp.sqrt(b * b - 4 * g)) / 2)


def N_all(T):
    return T / (2 * math.pi) * math.log(T / (2 * math.pi * math.e)) + math.log(T) + 1


def logN_density(s, T, bell_scale=1.0):
    """log of an upper bound for N(s, T) (zeros with beta > s, 0<gamma<T)."""
    lT = math.log(T)
    cands = [math.log(N_all(T))]
    # KLN: use the largest tabulated sigma0 <= s
    row = None
    for r in KLN:
        if r[0] <= s:
            row = r
    if row is not None and T >= H0_KLN:
        s0, A, B = row
        v = A * lT ** (2 * s0) * lT ** (5 - 4 * s0) * T ** (8 / 3 * (1 - s0)) + B * lT ** 2
        cands.append(math.log(v))
    if s >= BELL_ALPHA and T > 3e12:
        cands.append(math.log(1.62e11 * bell_scale) + 1.448 * (1 - s) * lT)
        if T <= math.exp(29):
            cands.append(math.log(37341.72 * bell_scale) + 14.160 * (1 - s) * lT)
    return min(cands)


def tail(L, H, T1, bell_scale=1.0, ns=400):
    if T1 <= H:
        return 0.0
    total = 0.0
    T = H
    while T < T1:
        T2 = min(2 * T, T1)
        bmax = 1 - 1 / (Rof(math.log(T2)) * math.log(T2))
        blk = math.exp(-0.4 * L + math.log(N_all(T2)))
        blk += math.exp((0.6 - 1) * L + logN_density(0.6, T2, bell_scale))
        h = (bmax - 0.6) / ns
        acc = 0.0
        for i in range(ns):
            s = 0.6 + (i + 0.5) * h
            acc += math.exp((s - 1) * L + logN_density(s, T2, bell_scale))
        blk += L * acc * h
        total += 2 * blk / T
        T = T2
    return total


def total_err(L, H, tol, bell_scale=1.0):
    """Best split of the budget; returns (error, T1)."""
    best = (float("inf"), None)
    for frac in [0.3, 0.5, 0.7, 0.9]:
        eps_s = frac * tol
        c = math.log(100 * L / tol)          # kernel term = 1% of budget
        T1 = c / eps_s
        e = (eps_s + L * math.exp(-c)
             + math.exp(-L / 2) * math.log(H) ** 2 / (2 * math.pi)
             + tail(L, H, T1, bell_scale))
        if e < best[0]:
            best = (e, T1)
    return best


def scan(H, Ls, tol, bell_scale=1.0):
    bad = [L for L in Ls if total_err(L, H, tol[L], bell_scale)[0] > tol[L]]
    return bad


def main():
    Ls = list(range(656, 3160, 25)) + [3158]
    tol = {L: eps_max(L) for L in Ls}

    print("L     tol        T1(needed)   tail@H=3e12   tail/tol")
    for L in [656, 1000, 1500, 2000, 2500, 3000, 3158]:
        if L not in tol:
            tol[L] = eps_max(L)
        e, T1 = total_err(L, 3e12, tol[L])
        tl = tail(L, 3e12, T1)
        print(f"{L:5d} {tol[L]:.3e}  {T1:.3e}   {tl:.3e}    {tl/tol[L]:.3e}")

    print("\nRH height H  -> L values in [656,3158] still FAILING")
    for H in [3e12, 1e13, 3e13, 1e14, 1e15, 1e16]:
        bad = scan(H, Ls, tol)
        rng = f"{min(bad)}..{max(bad)}" if bad else "none -- gap closed"
        print(f"  H={H:.0e}: {rng}")

    print("\nH=3e12, Bellotti constant C scaled down -> failing L")
    for sc in [1, 1e-2, 1e-4, 1e-6]:
        bad = scan(3e12, Ls, tol, sc)
        rng = f"{min(bad)}..{max(bad)}" if bad else "none -- gap closed"
        print(f"  C x {sc:g}: {rng}")


if __name__ == "__main__":
    main()
