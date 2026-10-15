# Lemma F — sharp explicit Barban–Vehov mean square (route + status)

Target: B := eps*log(z2/z1)*sum_n Psi(n)^2 n^{-1-eps} <= (1+eta)*Binf(tau,kappa),
Binf = kappa(tau-1)(1-e^{-kappa}) + (1/(tau-1)) int_1^tau e^{-kappa v}(1+kappa(tau-v))^2 dv.

## Contour-shift route: abandoned
- To shift the double Mellin integral left you need bounds for 1/zeta on Re s < 1.
- Classical zero-free region: the shift is only theta ~ 1/(5.56 log t), which gains nothing at log z1 ~ 55.
- Platt's RH to 3e12 plus Borel–Carathéodory: |log zeta| bounds of order 50–100 at sigma = 0.6–0.75, so 1/zeta is about e^{100}. Dead.

## Elementary route (no zeros needed)
1. Selberg diagonalisation, with s = 1+eps and J_s(m) = prod_{p|m}(p^s-1):
   - sum_n Psi^2 n^{-s} = zeta(s) * sum_m mu^2(m) J_s(m) y_m^2
   - y_m = mu(m) m^{-s} [N(z2/m; m) - N(z1/m; m)] / log(z2/z1)
   - N(y; q) = sum_{k<=y, (k,q)=1} mu(k) k^{-s} log(y/k)
2. Exact convolution: N(y; q) = sum_{d | q^inf} d^{-s} N(y/d; 1).
3. N(y; 1) = a log y + b + r(y), with a = 1/zeta(s) and b = -zeta'(s)/zeta(s)^2 (the residue).
   - r(y) = -int_y^inf (m_s(t) - a) dt/t
   - |m_s(t) - a| <= 2c t^{-eps}/log t, where m(t) = sum_{k<=t} mu(k)/k and |m(t)| <= c/log t
     (Ramaré: c = 0.0144 for t >= 96955).
   - Hence |r(y)| <= 2c * E1(eps log y).
4. The main terms sum over d back exactly to the residue for q. The error picks up the same q^s/phi_s(q) factor, plus a tail over d > y (bounded by Rankin).
5. For small y, use the direct bound |N(y; q)| <= 1.00303 q/phi(q), or exact tables of N(y; 1) for y <= 1e7.
6. Sum over m: use explicit asymptotics for sum mu^2(m)/phi(m) (Ramaré–Akhilesh type), twisted by m^{-eps} via partial summation.

## Model (lfm.py, q/phi factors treated as averaging out)

| tau | kappa | eta |
|---|---|---|
| 1.6 | 0.2–0.5 | 0.066–0.081 |
| 1.8 | 0.2–0.5 | 0.061–0.072 |
| 2.4 | 0.2–0.5 | 0.049–0.063 |

So eta ~ 0.1 looks reachable. That is better than the 0.2 assumed in the payoff model (gain 6.75–8x in density, which suggests H around 1e13).

## Remaining to make it rigorous
- Step 6 with explicit error terms, including how the q/phi factor correlates with the weights. This is the main risk.
- The d > y tail.
- An exact small-y table.
- Then plug in as d4_bound, rebuild the tables, bisect H.

## Pointwise test with true arithmetic (fb.c)
Setup: tau = 1.8, kappa = 0.3, B∞ = 0.8986. Each |Y_m| is replaced by |main_m| + err_m, using the actual q/phi(q) and d|q^inf structure and the exact r(y).

| z1 | true B | bound | bound/B∞ | (bound/B∞ - 1)*log z1 |
|---|---|---|---|---|
| 300 | 0.733 | 1.317 | 1.466 | 2.66 |
| 1000 | 0.759 | 1.257 | 1.398 | 2.75 |
| 3000 | 0.777 | 1.214 | 1.351 | 2.81 |
| 10000 | 0.792 | 1.178 | 1.311 | 2.86 |

The excess over B∞ behaves like about 3/log z1, i.e. a boundary layer near y ~ 1. Extrapolated to log z1 = 55 that is about 0.055 from the boundary layer. Adding about 0.07 from the E1 tails puts the expected eta around 0.10–0.13. The correlation between q/phi(q) and the weights does not blow up the bound.

## RIGOROUS version: lemmaF.py bound2 (Oct 9)
The proof chain (also in the lemmaF.py docstrings):
1. Selberg diagonalisation. W(m) = prod 1/(p^s-1) <= m^{-eps} mu^2(m)/phi(m).
2. Exact convolution: N(y;q) = sum_{d|q^inf} d^{-s} N(y/d;1). This holds for ANY main term that is linear in log u.
   We take main1 = 1 + 2 eps J + eps log u, with J = int_1^U0 (mcheck(t) - 1) dt/t (numerically -gamma).
3. Remainder: r(u) = u^{-eps} rho(u) + eps int_1^u rho(t) t^{-eps} (2 + eps log(u/t)) dt/t, with rho = mcheck - 1.
   - This is exact; it comes from two partial summations.
   - rho is computed exactly below U0 = 1e6.
   - Above U0 we use Ramare's |mcheck(t) - 1| <= 0.00257/log t (t >= 3846).
   - Resulting sup|r| <= 0.0083, versus 0.06 with the 0.0144/log t bound for m(t).
4. d-average:
   - The weights d^{-s}/H sum to 1.
   - Cut at D, with D chosen per interval.
   - Rankin bound Pr(d > D) <= D^{-al} g(m), with g(m) = prod (1-1/p)/(1-p^{al-1}).
   - Hp term: 2 eps h(m).
5. Small y (< U0): N(y;q)/H lies in [min(0, running min of N(u;1)), running max of N(u;1)], since it is an average.
   For m <= z1 with y1 small, use the windowed min with Rankin, or the residue path, whichever is smaller.
6. Sum over m:
   - m <= 1e6 is summed exactly.
   - Beyond that, log-intervals of width 0.01. On each interval: endpoint maxima, and the multiplicative sums
     sum mu^2/phi * {1, g, g^2, h^2} <= P log(b/a) + 7.9 E/sqrt(a). The Euler products use the
     pi(t) <= 1.25506 t/log t tail.
   - Input |G(x) - log x - 1.3326| <= 3.95/sqrt(x), x >= 1. Checked to 1e8 (max 1.45).
   - The square is taken in full for main + Rankin (no Minkowski); Minkowski is used only for the tiny Hp piece.

Checks:
- Every m at z1 = 300 and z1 = 2000: true |Yt_m| <= the pointwise bound, 0 violations (max ratio 1.000, tight).

| Lz | B/B∞ |
|---|---|
| 55 | 1.17-1.40 |
| 80 | 1.06-1.47 |

(tau 1.6-2.4, kappa up to 1.4)

## Integration
- lemmaD.d4_bound = min(RZA, Lemma F). Lemma F uses lemmaF_table.pkl, a grid of (Lz, tau, log kappa) interpolated
  trilinearly, times F_SAFE = 1.03.
- verify_F.py recomputes the exact bound at every chosen point and must report 0 violations.
- At the old optimum, d4 drops about 4.7x. After re-optimisation, logN drops by 1.8 at the zero-free edge, i.e. density / 6.2.
  New optimum: u ~ 2.4, tau ~ 1.82, w ~ 0.7, x ~ 6.

## Result: RH to 1.9e13 suffices (was 4.45e13)

| L range | worst ratio |
|---|---|
| 58-656 | 0.513 |
| 656-1300 | 0.875 |
| 1300-1600 | 0.958 (binding, L ~ 1380) |
| 1600-1900 | 0.828 |
| 1900-2400 | 0.633 |
| 2400-3159 | 0.521 |

Bisection: 1.75e13 fails (1.030 at L = 1370). The true threshold is about 1.8e13.

## Referee round 1 (Oct 9-10) and version 4

### Fixes from the referee
- FATAL, lemmaD.logN: the constants were computed at L = lk but are not monotone in L.
  Now N(s,T) <= N(s, e^lk), with no rescaling in T.
- Lemma F interval logic: each interval now lies in a single regime, every quantity is taken at its worst endpoint,
  and the window minimum is taken over the union window.
- IG raised to 5.922 (supremum over sigma in [-1/2, 1/2 - alpha] for delta up to 0.03).
- GERR = 4 (Ramare-Akhilesh 2017, 4/sqrt D, D > 1).
- Euler-product tail constants corrected.
- assert eps <= al. This justifies the Rankin step: 1 + (eps/R1) log(d/D) <= (d/D)^al.
- interval.py certifies every interval at tol(L2 + 1). This covers the y = x/e requirement.
  The scan must start at L = 57.

### Lemma R (own proof; removes the unlocated 0.00257 citation)
For x >= 2e6: |mcheck(x) - 1| <= 0.0121/log x.

Proof.
1. Since sum_{d|n} mu(d) = [n=1], we have sum_{n<=x} mu(n)/n H_[x/n] = 1, with H_N = sum_{k<=N} 1/k.
2. Write H_[y] = log y + gamma + delta(y) and delta = beta + eta, where beta(y) = (1/2 - {y})/y.
3. Since sum mu(n)[x/n] = 1, we get sum mu(n){x/n} = x m(x) - 1. Hence
   sum mu(n)/n beta(x/n) = M(x)/(2x) - m(x) + 1/x,
   and so
   rho = mcheck - 1 = (1 - gamma) m(x) - M(x)/(2x) - 1/x - T(x),  where T = sum_{n<=x} mu(n)/n eta(x/n).
4. Properties of eta:
   - eta is continuous on [1, inf). The jumps of delta and beta at the integers coincide.
   - |eta'(y)| = |[y] + 1/2 - y|/y^2 <= 1/(2y^2).
   - eta(1) = 1/2 - gamma.
   - V = int_1^inf |eta'| <= 0.25844 (closed form on each unit interval).
5. Stieltjes: T = eta(1) m(x) - int_1^x m(t) d[eta(x/t)].
6. Bound the integral in three ranges:
   - t >= x/Y: |m(t)| <= 0.010032/log t (t >= 617990) [arXiv 2408.05969, Cor. 1.2].
   - K0 <= t < x/Y: the same bound with log t >= log K0.
   - t < K0 = 1e6: exact |m_k|; the variation over each [k, k+1) is <= 1/(2x); sum_{k<1e6} |m_k| = 279.68.
7. Use |M(x)| <= 0.006688 x/log x (x >= 1798118) [ibid., Thm 1.1].
8. Optimise over Y on a grid. max_{x >= 2e6} rho_bar(x) log x = 0.012038.

Checks:
- The identity was verified numerically to 1e-10 at x = 1234.5, 99999.3, 777777.7, 999999.9.
- With only Ramare's published 1/69 and 0.013 the constant is 0.0187. (The 0.013 M-bound is not checked here.)

### Version 4 of the remainder
- The t^{-eps} weights in the rho-integrals are kept exactly: midpoint value per piece plus an explicit error.
- The eps^2 log u drift of r is absorbed into the main term: main1 = 1 + c0 + A log u, with A = eps(1 + eps j)
  and c0 = eps(2j - eps j2). This is legitimate because any linear main term is allowed.
- Above U0 = 2e6: Rbar(u) = 0.0121/log u + (eps terms, ~5e-4), which decreases in u.
- Brute check of |true r| / Rbar over all pieces up to 2e6, for (kappa, Lz) in {(0.5, 77.4), (1.4, 30), (0.1, 90), (3, 20)}:
  max ratio 0.98, 0.86, 0.997, 0.60.
- Rbar tail is now 0.0015, versus 0.038 with Ramare's 3/14.

### Other improvements
- Full square including the Hp piece. The coefficient is 2 eps h, which covers A Hp plus E[eps log d].
  This uses additive first and second moments:
  - sum mu^2/phi * a <= sum_p a_p/(p-1) [log(b/a) + 7.9 sqrt(p/a)]
  - sum mu^2/phi * a b <= [(sum a/(p-1))(sum b/(p-1)) + sum ab/(p-1)] log(b/a) + (the sqrt(p) analogue) * 7.9/sqrt(a)
  - sum mu^2/phi * g a <= sum_p a_p g_p/(p-1) [Pg1 log + Eg1 7.9 sqrt(p/a)]
- Option C for small y1: P(d = 1) = 1/H >= 1 - k(m), with k(m) = sum_{p|m} 1/p, so N1/H >= n1 - k(m)(n1 - lo').
  The option is chosen once per interval among A, B and C.

### Effect
- At (Lz, tau, kappa) = (77.4, 1.823, 0.5): B/B∞ went 1.31 -> 1.13.
- At the binding cells: (60, 33) logN 5.92 -> 5.61, and (70, 32) 6.30 -> 5.99.

## Result after referee round 1 + v4 (Oct 10): RH to 1.75e13 suffices
All constants are now sourced or proved here: Lemma R replaces the unlocated 0.00257.
Tables lemmaD_tableG_* (key version 6) were fully re-optimised.
verify_F: 2210 cells, 0 violations, max exact/used 0.973. The edge cell (49, 37) is OK.

Scans with tol(L2+1), starting at L = 57:

| L range | worst ratio |
|---|---|
| 57-656 | 0.521 |
| 656-1300 | 0.918 |
| 1300-1600 | 0.976 (binding, L = 1350) |
| 1600-1900 | 0.822 |
| 1900-2400 | 0.628 |
| 2400-3159 | 0.520 |

At 1.65e13 the range 1300-1600 fails (1.030 at L = 1350).
