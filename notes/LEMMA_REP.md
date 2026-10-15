# Lemmas H and Rep (Oct 14): height-dependent zero-free region and two-zero repulsion, from BTY's method

Source machinery: Bellotti-Trudgian-Yang, arXiv:2603.21490 (Lemma 15 and Lemmas 3-13). Code: rep/rep.py, rep/dtable.py.

## What we change in BTY's proof
(i) We evaluate their final inequality EXACTLY at each height t (W(z) by quadrature of the explicit w; F = sum kappa_m W((z+delta_m)/eta)),
    instead of the uniform fit C1 + C2 and the worst case t = exp(76.47). The inequality (their (44) with (45)-(48)):
      0.1185 f(0) <= eta w0 [ (a kappa/2) log t - 1.568 ] + a0 F(sigma-1) + 1e-10 eta - Z + 1.5 (187 eta^2 + 8013 eta^3),
      Z >= a1 [F(sigma-1+eta0) + F(sigma-eta0)] (+ second-zero terms, below) - 1e-7,
    with f(0) = kappa eta w0 (kappa = 433/859).  NOTE: BTY write the Lemma-10 term as -0.1186 w(0) eta, i.e. f(0) = eta w(0);
    with f(0) = kappa eta w0 the subtracted amount halves. We use the corrected form. Their Theorem 1 (4.896) still holds:
    our exact evaluation gives lambda* = 0.2054 at t = exp(76.47) with the correction (> 1/4.896 = 0.20425).
(ii) Larger eta (zeros up to lambda = 0.24): BTY's Lemma 11 (f >= kappa f(0) on [0,59]) needs eta <= eta0 = 0.0071; for
    eta <= 0.0085 it holds on [0, 49.5] only, so Lemma 10 keeps primes p <= 23 (m <= 15): sum (log p) m_p = 0.23522, giving
    0.1185 f(0) (was 0.1186). Error terms (Lemmas 8-9) recomputed constants change slightly with eta0; we inflate by 1.5.
(iii) Iteration with a height-dependent constant: sigma = 1 - A/log(Kt+T0) with A = lambda*(log(Kt+T0)) (valid since lambda(l)/l
    is decreasing). Converges in 3 steps.

## Lemma H (height-dependent zero-free region)
zeta(beta+i gamma) = 0, gamma >= 3e12  =>  (1-beta) log gamma >= lambda*(log gamma), with (grid, rep/lamstar.json)
  lambda*(28) = 0.21294, (29) 0.21263, (30) 0.21233, (31) 0.21203, (32) 0.21175, (33) 0.21147, (35) 0.21094, (37) 0.21046, (40) 0.20974
(R = 1/lambda* = 4.70 at log t = 29 versus BTY's uniform 4.896).

## Lemma Rep (two-zero repulsion)
Keep a second zero rho1 = 1 - eta1 + i(t+d) (and its Stechkin partner 1 - conj(rho1)) in the k = 1 zero sum of Lemma 6
instead of discarding its pair via (21) (pair terms are >= 0, so keeping one only strengthens Z):
   Z >= a1 [F(sigma-1+eta0) + F(sigma-eta0) + Re F(sigma-1+eta1-id) + Re F(sigma-eta1-id)] - 1e-7.
Contradiction for all (eta0, eta1) in [lambda*/L, Lambda/L]^2 (choosing either zero as primary) and |d| <= d*(L, Lambda).
The reflected term Re F(sigma - eta1 - id) ~ w0 eta0 (sigma-eta1)/((sigma-eta1)^2+d^2) ~ 0.04 is what gives repulsion at
distances d up to ~0.05-0.1.  Table rep/dstar.json (cell (lk, k): Lambda = k e-4 * lk):
   lk = 29: k=74 0.096, 75 0.088, 76-82 0.057-0.061;  lk = 30: k=71 0.135, 72 0.088, 73-80 0.057;  lk 31-37: ~0.045-0.09.
Used in lemmaD._lp2: representatives closer than d* are excluded (ntil = 1 already).
Effect: binding cells logN drop ~0.3 (e.g. (75,29) 4.39 -> 4.07).

## Rigour TODO (before the preprint)
- W(z) by trapezoid (20001 pts) -> arb/interval quadrature; F at complex arguments likewise.
- dstar via bisection assumes the contradiction set in d is an interval; verify on a fine d-grid; lambda grids (5x5) -> monotonicity
  argument in eta0, eta1 or finer grids with Lipschitz bounds.
- Re-derive Lemma 8/9 error constants for eta0 = 0.0085 (we used x1.5).
- Lemma 10 with p <= 23 needs BTY's monotonicity check dG_p/dsigma < 0 (unchanged range of sigma) - fine; m_p table reused.


# v2 (Oct 14, later): everything re-derived and certified; Lambda extended to ~0.0125 L
Code: rep/bty2.py (constants), rep/Wexact.py (closed-form W), rep/rep2.py (float), rep/lamcert.py, rep/rep_cert.py,
rep/repcert_run.py (ball arithmetic, python-flint), rep/lemma5_cert.py, rep/lemma11_cert.py, rep/consts_cert.py.

## Why the old version was not enough
Diagnosis at H = 3.1e12, L = 900 (diag900.py): ratio 1.12 = 0.60 (smoothing) + 0.52 (zero tail), the tail coming from
lambda = (1 - beta) log T in [0.20, 0.32] at log T in {29, 30}.  Half of it was above Lambda = 0.24, where repulsion was off
because BTY's constants are tied to eta <= eta0 = 0.0071 (we had stretched to 0.0085 with a x1.5 fudge).

## Re-derived BTY constants for general eta0 (and lam = eta log t of the primary zero)
* W0 bound (replaces Ford's Lemma 4): four integrations by parts (w(D)=w'(D)=w''(D)=0, w'(0)=0):
    |W(z) - w0/z| <= Cw(nu, r)|z|^-3,  Cw = |w2(0)| + (|w3(0)| + e^{-nu D}|w3(D)| + int_0^D |w4| e^{-nu u})/r   (wk = k-th derivative)
  |w2(0)| = 50.0697, |w3(0)| = |w3(D)| = 117.2648 (exact series), the integral certified on an integer nu grid
  (consts_cert.py).  Cw(69.4, 69.4) = 51.80 (BTY: < 52).
* Lemma 8: t >= H term is Cw/(pi x0^2) * lam * eta^2 (BTY use lam <= A0 here); t = 0 term uses
  int |Re psi(1/4 + iy/2)| (x0^2+y^2)^{-3/2} dy <= 23.6497 at x0 = 0.4930 (certified).  NOTE: BTY's bound U(y) in (31)
  integrates to 28.3 (29.98 with |U|), not < 24; the true digamma gives 23.59, so their (32) is right but not via (31).
  Reproduces BTY's (48): 185.4 eta^2 + 7972 eta^3 (theirs 187, 8013).
* Lemma 5 (non-negativity (21)): BTY verify it uniformly for sigma >= sigma0 using c0 = 151/153 <= 2 sigma0 - 1, i.e.
  sigma >= 0.993464.  Our height-dependent sigma is smaller (A = lam* ~ 0.213 > 1/4.86), where their sufficient test FAILS
  (min B < 0 at c0 = 0.9864).  Fix: at fixed sigma all delta_m = (2 sigma - 1) m exactly; rescaling y by (2 sigma - 1) the
  test becomes sigma-free:  B(y) = sum b_m (m - sgn(b_m) eps0/m)/(m^2 + y^2) > 0,  eps0 = Cw eta0^2/(w0 (2 sigma - 1)^2).
  Certified for eps0 = 1/500 by exact root isolation (lemma5_cert.py); this covers eta0 <= 0.0147.
* Lemma 11 at eta0: holds on [0, u0(eta0)] (u0 = 59.16 at BTY's eta0, 49.46 at 0.0085, 33.6 at 0.0125), certified in ball
  arithmetic (lemma11_cert.py; g(0) = 0, g'(0) = 1.19 > 0).  Lemma 10 then keeps primes with 15 log p <= u0 and we bound the
  dropped m > 15 terms ((a0 + a) p^{-16 sigma}/(1 - p^{-sigma}); BTY drop them as negligible without a bound).
  Constant: 0.0598 eta w0 (= 0.1186 f(0)) at BTY's eta0, 0.0595 at eta0 = 0.0115: insensitive.
* Gamma constant (Lemma 12) recomputed at our sigma0: 1.5686 (BTY 1.568).
* W(z) in closed form (Wexact.py; agrees with quadrature to 1e-15); certification uses the same formula in acb with a
  second-order enclosure W(zc) + W'(zc)(z - zc) + r^2/2 sup|W2| (W' in closed form, sup|W2| <= (2 w0 + int|(u^2 w)_3|)/|z|^3
  bounded by piecewise ball sups).

## Certified results
* Lemma H (lamcert.py, 4-step certified chain from BTY's 1/4.896; each step: for every height in each 0.25-cell and every
  lam0 between the old and new constant, the single-zero inequality is violated, in ball arithmetic):
  lam*(28) 0.21424, (29) 0.21384, (30) 0.21345, (31) 0.21308, (32) 0.21272, (33) 0.21239, (35) 0.21175, (37) 0.21117,
  (40) 0.21033  ->  rep/lamstar2c.json.   (R = 4.676 at log t = 29.)
* Lemma Rep (rep_cert.py): for each cell (lk, k) [heights (lk-1, lk], 1 - beta <= k e-4, Lambda = k e-4 lk up to
  0.0125 (lk - 1)], certified over the full continuum of (height, eta0, eta1, |d|) boxes, with sigma at the actual
  height, either ordering of the pair, at 97% of the float d*.  d* ~ 0.058 at lk = 29 up to Lambda = 0.35; cells between
  computed k's inherit the value at the next computed k above (excluded set grows with Lambda).  -> rep/dstar2c_*.json
* J-lever (Lemma D): Bellotti's (4.16) used 66.7 a^{27/164}; with HPY (0.618 u^{1/6} log u, |u| <= e^105) the whole integral
  is <= 4.72 * 0.618 T^{1/6} log T (c = 3/2 - beta_r - beta_s <= -0.44).  ~ -0.03 in log N.
* Half-unit height cells below e^31 in Lemma D (lk = 29.5, 30.5 use the integer d* cell above).
* Lemma 10 monotonicity (BTY verify dG_p/dsigma < 0 numerically for sigma >= 0.99352): not needed numerically at all.
  G_p(sigma, x) = sum_{m<=15} p^{-m sigma} P(m x) (P includes a_0), so dG_p/dsigma = -log p sum m p^{-m sigma} P(mx) <= 0
  since P >= 0.  Holds for every p and every sigma, in particular for our smaller sigma0.  (Ball check for p <= 29 agrees.)

# Oct 15 hostile audit (independent agent) and fixes
1. SOUNDNESS BUG (fixed): rep_cert.PQ (and the earlier mean-value Mnu used by lamcert) clamped nu at -2; P, Q grow as nu
   decreases and boxes reach nu ~ -3.3, so the W'' radius was under-estimated. Clamp removed. Lemma H and Rep certificates
   (and everything downstream) must be re-run: rerun_all.sh.  Spot re-checks with the fix pass ((29,80) at d = 0.0564; Lemma H
   cell 29), and the auditor re-certified (29,74) at d = 0.1423.
2. Lemma 10 truncation: dropped primes / m > 15 terms are >= 0 because f >= 0 (BTY Lemma 3), now certified (f_nonneg_cert.py:
   exact root isolation for p(x), ball bisection + endpoint Taylor for w). The (mis-normalised) tail correction removed.
3. Second zero at height t + d: Lt ball widened by 1e-12.
4. w-derivatives now exact arb values; BTY Table 1 re-certified independently (mp_cert.py): ALL hold EXCEPT p = 5:
   min_x G_5(1, x) = 0.00735649 (x = 1.64952), below BTY's "certified" 0.00735709. Replaced by 0.0073564 (certified).
   Effect on BTY's 0.23545: -1e-6 (harmless) -- remark for the preprint.
5. HPY range: for |u| > e^105, 66.7 u^{27/164} <= 0.618 u^{1/6} log u, so the majorant holds for all |u| >= 3.
Verified by the auditor: Lemma 6 with a kept pair; Lemma 5 rescaling and eps0; min-over-lower-cells d* usage and half bins;
Lemma H induction and cell coverage; error constants; Taylor remainder; J replacement (4.72, float, 0.13% margin).
