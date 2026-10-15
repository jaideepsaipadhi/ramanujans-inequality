# Lemma G (target): sharp explicit bound for the Halasz diagonal d4 = sum_n Psi(n)^2 w(n)/n

## Where the remaining loss is (lever.py, binding cell delta = 0.006, log T = 33)
Re-optimised change in log N for each hypothetical change:

| change | delta log N |
|---|---|
| d4 / 2 | -0.71 |
| B-term / 1000 | -0.30 |
| J (c5) / 1000 | -0.19 |
| (d5 - 1) / 2 | -0.03 |

So d4 is the lever, with elasticity about 1. The others are second order.

## The loss is in Lemma E, not Lemma F
Lemma E majorises the weight: e^{-n/X} <= (eps/e)^eps (n/X)^{-eps}.
- That turns the TRUNCATED mean square sum Psi^2 e^{-n/X}/n into the FULL sum sum Psi^2 n^{-1-eps}.
- Here X = T^x lies between V = T^v and V^2: x ~ 5.8, v ~ 4.1.
- The power weight drags in the tail n > X, which contributes about 1/(eps (tau-1) log U). This costs about a factor e x/(x-v) overall.

Asymptotic truth (Laplace inversion of B∞; pay.py). Let C(Θ) = sum_{n <= U^Θ} Psi^2/n.
- C = 1 for Θ <= 1.
- On [1, tau]: C is a piecewise quadratic.
- C(tau) = 3/2.
- C(Θ) = 3/2 + (Θ - tau)/(tau - 1) for Θ >= tau. This is the Graham-type linear regime, which already holds for Θ >= tau, not only Θ >= 2 tau.

Brute force (sx.c), tau = 1.71:

| log U | x/u | true S_X | C (asymptotic) |
|---|---|---|---|
| 6 | 1.3 | 1.045 | 1.089 |
| 6 | 1.71 | 1.305 | 1.5 |
| 6 | 2.0 | 1.592 | 1.908 |
| 6 | 2.41 | 2.013 | 2.486 |

So the truth sits below the asymptotic, as it did for B.

At the binding cell: Lemma E gives about 6.8 for S_X, while the truth is about 2.5 (d4 about 5.8 vs about 1.5).

Payoff, using d4 = (1 + eta) C(x/u) - 1, re-optimised:

| eta | delta log N | density factor |
|---|---|---|
| 0.15 | -1.17 | / 3.2 |
| 0.5 | -0.78 | |
| 1.0 | -0.41 | |

## Real-axis information cannot capture this (lp.py)
- Upper bounds on F(1+eps) for real eps give at best inf over positive exponential majorants, which is about Lemma E.
- Allowing signed combinations, with two-sided Lemma F at relative error eta, only gains 10-20% (eta = 0.13: 4.65 -> 4.18; truth 1.49).
- The truncation is an oscillatory (complex-w) phenomenon.

## Proof architecture (route I: vertical-line Mellin with complex Lemma F)
1. S_X = (1/2 pi i) int_{(eps)} Gamma(w) X^w F(1+w) dw, with F(1+w) = sum Psi^2 n^{-1-w} = zeta(1+w) sum_m J_{1+w}(m) y_m(w)^2 (Selberg diagonalisation, valid for complex w).
2. Complex remainder:
   N_w(u;1) = 1 + w log u + r_w(u),  r_w(u) = u^{-w} rho(u) + w int_1^u rho(t) t^{-w} (2 + w log(u/t)) dt/t.
   - This is the same partial summation as for real w, with Lemma R.
   - Relative error |r_w|/|main| is about |w| (0.45 + 0.05 log u)/log u.
3. Main term: F_asym(w) = B∞(tau, w log U)/(w Delta), the analytic continuation.
   - The w-integral of F_asym is explicit: it is C convolved with the Gumbel kernel of e^{-n/X}.
   - The error is int |Gamma(eps+it)| X^eps |F - F_asym| dt.
   - With |F - F_asym| <= eta |F_asym| it is about eta (e^a/2a) x/(u(tau-1)) B∞, minimised at a = 1, giving about 4.2 eta.
   - So eta = 0.13 gives an error of about 0.55, and d4 is about 2.0 instead of 5.8.
4. Pieces still to build:
   - (a) complex analogues of all Lemma F pointwise bounds (Rankin, ranges, option C);
   - (b) an explicit asymptotic for the complex-weighted m-sums sum mu^2(m) f_w(m) m^{-1-w} (Dirichlet series zeta(1+w+s) E_w(s), E_w absolutely convergent for Re s > -1/2);
   - (c) |F(1+w)| <= F(1+eps) for |t| > T1 (positivity), and the Gamma tail.

Route II (contour shift into Re w < 0, using RH to 3e12 and rigorous 1/zeta bounds in a box) was considered and set aside.
The left shift picks up (V^2/X)^{c} growth because V^2 > X. That growth is exactly the Graham cancellation, and it needs the oscillatory treatment of route I.

## Oct 10, session 2: what was tested

1. The easy regime X >= V^2 T^{0.45} is NOT worth it (pay2.py).
   In that regime the Selberg error terms are trivially small, so d4 would follow from Q(1), Q'(1) alone.
   But forcing x >= 2v + 0.45 costs more than the gain: delta log N = +0.21 even at eta = 0.05.

2. The Mellin route works numerically at small scale (qc.c exact complex F, drv.py).
   Setup: z1 = 403, tau = 1.71, x/u = 2.41, a = 1.
   - true S_X = 2.012 (direct). Mellin of the exact F gives 2.0196, so the quadrature is fine.
   - Mellin of the model F_asym(w) = B∞(tau, w Lz)/(w Delta) gives 2.359.
   - Error integral (1/2pi) int |Gamma| X^eps |F - F_asym| dt = 0.597.
   - So the route gives S <= 2.96, against Lemma E = 5.11.
   - Relative error |F - F_asym|/|F| is 21% near t = 0 (boundary layers at this tiny scale) and below 2% for t >= 1.

3. The obstruction to making it rigorous.
   - For |t| >> 1/log U, F(1+w) - 1 is small, about 1/(|t| Delta), only because of cancellation in the oscillating
     m-sum sum_m W_w(m) Yt_m(w)^2.
   - Pointwise diagonal bounds (|W_w| <= W_eps, |Yt - Yhat| <= e_m) lose that cancellation.
     They give errors that grow like eta |w|^2 log z, against a true size of about 1/(|t| Delta).
   - So they are only useful for |t| <~ 1/Lz.
   - Positivity gives only |F(1+w)| <= F(1+eps), which has no decay in t.
     The resulting tail int_{t0}^inf |Gamma| is about (1/pi) log(1/t0) X^eps F(1+eps), the same size as Lemma E.
   - A Gaussian-smoothed majorant c Phi(log(X'/n)/sigma) has exact Mellin transform c X'^w e^{sigma^2 w^2/2}/w.
     It cuts |t| <= ~6/sigma, but costs about sigma/Delta extra mass.
     Making the t-cut small enough forces sigma ~ Lz, which costs more than the gain.
   - Mean-value (Montgomery-Vaughan) bounds over t do not help: m runs to e^{140}, far beyond the t-range.

## Honest status
- Lemma G (about / 3 in density, about -1.2 in log N) needs genuine control of the oscillating diagonal for
  |t| up to about 1. That is a Graham-type explicit theorem; I do not see an explicit route with the tools here.
- A smaller provable lever: signed two-power majorants with a two-sided Lemma F (needs a lower-bound version of F).
  The LP says d4 drops 10-20% at eta = 0.13, i.e. about -0.1 to -0.2 in log N.

## Oct 10, session 3: a feasible rigorous route (Graham-type contour argument)

Decompose Psi = (Lambda_V - Lambda_U)/Delta, where Lambda_z(n) = sum_{d|n} mu(d) log+(z/d). Then
  S = Delta^{-2} sum_{i,j} (+-) T(z_i, z_j),   T(z_i, z_j) = sum_n Lambda_{z_i}(n) Lambda_{z_j}(n) e^{-n/X}/n,
  T = (2 pi i)^{-3} int int int Gamma(w) X^w z_i^{s1} z_j^{s2} G(w,s1,s2) dw ds1 ds2/(s1^2 s2^2),
  G = zeta(1+w) zeta(1+w+s1+s2) / (zeta(1+w+s1) zeta(1+w+s2)) * A(w,s1,s2).
- A is an absolutely convergent Euler product near 0, with A(w,0,s2) = A(w,s1,0) = 1.
- Main term = residues at s_i = 0 (double), s1 + s2 = -w, and w = 0 (Gamma, zeta).
- Asymptotically this reproduces C(Theta).

Error control: shift each s-contour to Re(1+w+s) = 1 - theta for |t| <= T0, and keep Re >= 1 + 0 beyond T0.
Three ingredients:

(a) Box: rigorous enclosure of |1/zeta(sigma+it)| for sigma >= 1 - theta and |t| <= T0, by acb ball arithmetic on cells.
    It is analytic there because RH holds to 3e12.
    Non-rigorous scan (zbox.py), max |1/zeta| for t <= 2000:

    | sigma | max abs(1/zeta) |
    |---|---|
    | 0.6 | 20.1 |
    | 0.7 | 9.1 |
    | 0.8 | 5.5 |
    | 0.9 | 4.2 |
    | 1.0 | 3.4 |

    The scan runs at about 6.7e3 points/s, so a box to T0 = 1e5 is cheap.

(b) Tail |t| > T0, from arXiv 2405.04869 (Cor. 4):
    - |1/zeta(sigma+it)| <= 30.812 log t for sigma >= 1, t >= 13;
    - |1/zeta(1+it)| <= 29.388 log t for t >= 2.
    With the 1/s^2 kernels the tail is about 900 (log T0)^3/T0, about 14 at T0 = 1e5.
    Relative to T ~ log^2 z ~ 6000, i.e. after dividing by Delta^2 ~ 3000, that is about 0.005 in S.

(c) Shifted-line part: z^{-theta} with z = e^{80} and theta = 0.3 gives e^{-24}, times the box bounds. Negligible.

So Lemma G looks FEASIBLE with explicit, checkable inputs.

Build order:
1. Exact residue main term for T(z_i, z_j; X) in closed form. Validate against brute force at small scale (sx.c).
2. Rigorous box computation of 1/zeta and zeta'/zeta (flint) for |t| <= T0, sigma in [1 - theta, 1.05].
3. Error assembly for the triple contour.
4. Plug d4 = S_main + E - (n=1 adjustment) into lemmaD, rebuild the tables, rescan. Target H ~ 1.0e13.

## Step 1 progress (Oct 10): residue main term, validated
F_ij(w) := sum_n Lambda_{z_i}(n) Lambda_{z_j}(n) n^{-1-w}, with z_i <= z_j. Then
  F_ij(w) = log z_i log z_j + (zeta'/zeta)'(1+w) + A_12(w;0,0) - z_j^{-w}/w^2 + I_x(w) + E(w), where
  I_x(w) = (1/2 pi i) int_{Re s2 = c} z_i^{-w-s2} z_j^{s2} zeta(1+w) A(w;-w-s2,s2) / (zeta(1-s2) zeta(1+w+s2) (w+s2)^2 s2^2) ds2.

Validation (lemG/t7.py) of exact (finite double sum) minus main; the remaining diff is the left-contour error, which decays in z:

| z_i | z_j | w | exact - main |
|---|---|---|---|
| 60 | 60 | 0.05 | 0.78 |
| 300 | 300 | 0.05 | 0.25 |
| 60 | 200 | 0.05 | 0.51 |

Without I_x, the error is about -0.6/w. That is exactly Q(1) - log z = -0.607, the constant in the Selberg form.

Two cases for I_x:
- i = j: the z-powers cancel, so I_x = z^{-w} K(w) with K independent of z. It does not oscillate and is cheap to compute.
- i != j: shift left. I_x = (residues at s2 = 0 and s2 = -w) + O((z_i/z_j)^{a}). At real scale V/U = e^{57}, so the remainder is negligible.

Next steps:
1. The w-integral. Shift w left: X^w z^{-w} decays because X > V. The main term T_main is then a finite sum of residues at w = 0.
2. Compute S_main at the binding parameters and compare with C(Theta) and Lemma E.
3. Rigorous error terms.

## Step 2 (Oct 13): main term in closed form -- DONE (modulo rigorous remainder)
Combine the four T_ij with Psi = (Lambda_V - Lambda_U)/Delta. Using the validated F_ij formula:
- constants log z_i log z_j combine to Delta^2 -> exactly e^{-1/X} (since (1/2pi i) int Gamma(w) X^w dw = e^{-1/X});
- (zeta'/zeta)'(1+w) and A_12 cancel (coefficients 1 - 2 + 1);
- i != j: I_x = (V^{-w} - U^{-w})/w^2 + O((U/V)^a)  (residues at s2 = 0 and s2 = -w; A = 1 there), so F_UV = ... - U^{-w}/w^2;
- i = j: I_x = z^{-w} K(w) with the UNIVERSAL kernel K(w) = zeta(1+w) H(w),
  H(w) = (1/2 pi i) int_{Re s = 0.3} A(w;-w-s,s) / (zeta(1-s) zeta(1+w+s) (w+s)^2 s^2) ds.
Total: S - e^{-1/X} = Delta^{-2} [ J_-(Y1) + J_+(Y2) ],  J_-+(Y) = Res_{w=0} Gamma(w) e^{Yw} (K(w) -+ 1/w^2) + (exp. small, w shifted left),
Y1 = log X/V, Y2 = log X/U. K has a SIMPLE pole: K = k1/w + k0 + ..., k1 = H(0), k0 = H'(0) + gamma H(0).
  D4_main = Delta^{-2} [ (Y2^2 - Y1^2)/2 - gamma Delta + k1 (Y1 + Y2 - 2 gamma) + 2 k0 ].
Numerics (lemG/kcoef.py, primes <= 2e4, |t| <= 200): k1 = -0.60727 (= the Selberg constant Q(1) - log z), k0 = 0.35189.
Leading order reproduces C(Theta) = 3/2 + (Theta - tau)/(tau - 1) exactly.

Validation against brute force (sx.c, log U = 6, tau = 1.71):
| x/u | formula S | true S |
|---|---|---|
| 1.71 (Y1 = 0, edge) | 1.299 | 1.305 |
| 2.0 | 1.5915 | 1.592 |
| 2.41 | 2.0045 | 2.013 |
At the binding cell (L = 33, u 2.424, tau 1.712, x 5.836): S = 2.4356, d4 = 1.436 (was 5.756 from Lemma F).

## Payoff (lemmaD.d4_bound_G, env LEMMAG=1, LEMMAG_ERR = budget for the rigorous remainder; requires Y1 >= 40)
Binding cell (60, 33), base logN 5.608:
| G_ERR | logN | delta |
|---|---|---|
| 0.01 | 4.183 | -1.425 |
| 0.3 | 4.362 | -1.246 |
| 0.6 | 4.509 | -1.099 |
| 1.0 | 4.668 | -0.940 |
H ~ 1.0e13 needs delta ~ -1.17, i.e. remainder budget G_ERR ~ 0.45 on d4 -- a LOOSE target (30% of d4).

## Remaining for rigour (step 3)
|E| <= 0.45 total, i.e. Sum_ij |T_ij - T_ij,main| <= ~0.45 Delta^2 / 4 ~ 365 each, against T_ij ~ log z_i log z_j ~ 1e4 (3-4% relative).
Pieces: s-contour remainders (z^{-theta} on |Im| <= T1 via a rigorous 1/zeta box; tails at Re >= 1 via |1/zeta| <= 30.812 log t),
the w-line integral, Euler-product truncation, rigorous enclosure of k1, k0. Naive absolute bounding of the s-tails is too lossy
(needs T1 ~ 1e9); the tails must be organised so the s2-integral is done by residues first.

## Scan at H = 1.0e13 with G_ERR = 0.45 (Oct 13): ALL PASS, worst ratio 0.883 (L = 1206). See h_G_1.0e13_provisional.txt.
Tables: lemmaD_tableGG0.45_*.pkl (key version 7). Conditional on: (i) |remainder| <= 0.45 on d4 for Y1 >= 40, (ii) rigorous k1, k0.

## Step 3 (Oct 13): rigorous architecture -- the smoothing trick kills the tail problem
Problem: with the Barban-Vehov weights, hat rho(s) = (V^s - U^s)/(Delta s^2) decays only like 1/|s|^2, so s-tails beyond the
RH-certified box cost ~1/T1 with big log factors (T1 ~ 1e7 needed). Fix: SPLIT THE WEIGHTS.
  rho_s(d) = E_Z[ rho(d e^{-eta Z}) ]  (Z ~ N(0,1)),   rho_t = rho - rho_s,   Psi = Psi_s + Psi_t,
  S <= ( sqrt(S_s) + sqrt(S_t) )^2   (Minkowski in l^2(e^{-n/X}/n)).
(a) Mellin: hat rho_s(s) = hat rho(s) E(s), E(s) = e^{eta^2 s^2/2}; on vertical lines |E| = e^{eta^2(sig^2 - t^2)/2}.
    So every contour can be closed inside the box |Im| <= T2 = 8/eta (+W0), Re >= 1 - theta, where RH (verified to 3e12)
    makes 1/zeta analytic and an arb cell computation bounds it; beyond T2 all lines sit at Re >= 1 (explicit bounds) and the
    Gaussian gives e^{-eta^2 T2^2/2}.
(b) |rho_t(d)| <= (eta/Delta) [psi((log d - log U)/eta) + psi((log d - log V)/eta)], psi(x) = phi(x) - |x| Phi(-|x|) (<= 0.399,
    int psi = 1/2). S_t via Lemma E + Selberg form: estimate S_t ~ 1e-5 at eta = 0.02, so Minkowski costs ~ 0.01.
(c) S_s = (2 pi i)^{-3} int Gamma(w) X^w hat rho_s(s1) hat rho_s(s2) G(w,s1,s2).  Contours: w on Re w = c; s2 on Re s2 = -c'
    (0 < c' < c; residue at s2 = 0 gives e^{-1/X} exactly since G(w,s1,0) = 1 and rho_s(1) = 1); s1 moved to Re -theta inside
    the box (residue at s1 = 0 vanishes: int over Re -c' of hat rho_s = 0; residue at s1 = -w-s2 gives R_x).
    R_x = Delta^-2 int Gamma X^w [ V^{-w}(K_c - E/w^2) + U^{-w}(K_c + E/w^2) ] + (cross-term shift remainders, Gaussian tails),
    K_c(w) = (1/2 pi i) int over the COMPACT path (-c'-iT2 -> 0.3-iT2 -> 0.3+iT2 -> -c'+iT2) of
             zeta(1+w) A(w;-w-s,s) E(w+s) E(s) / (zeta(1-s) zeta(1+w+s) (w+s)^2 s^2) ds  (analytic for Re w >= -theta').
    BUG caught while deriving: the U-cross-term residue at s2 = -w carries (V/U)^{-w}, giving -V^{-w} E/w^2 (not -U^{-w}).
(d) Main term (unchanged form, constants now k1(eta), k0(eta) from K_c; the c2 and eta^2/2 terms cancel):
    D4_main = Delta^-2 [ (Y2^2 - Y1^2)/2 - gamma Delta + k1 (Y1 + Y2 - 2 gamma) + 2 k0 ].
(e) Remainders, all exponentially small: E1 ~ U^{-theta} (s1 line), V-cross shift right ~ (U/V)^b, U-cross shift left ~ (U/V)^b,
    w-shift ~ e^{-theta' Y1}, Gaussian tails/connectors ~ e^{-eta^2 T2^2/2}; each times box sups of |zeta|, |1/zeta|, |A|.
    Rough size at the binding cell with theta = 0.3, b = 0.4, eta = 0.02, T2 = 400: << 1e-3.
So the rigorous d4 is D4_main + ~0.015 (Minkowski + remainders), far inside the 0.45 budget.

## Step 3 implementation (Oct 13)
Parameters: eta = 0.02, T2 = 400, W0 = 40, c = 0.02, c' = 0.01, sig0 = 0.01, theta = 0.40, b = 0.38, theta' = 0.30, kap = 0.38.
- Box (lemG/box2.py, Taylor models of (s-1)zeta(s), K = 14, Cauchy remainder with |zeta(s)| <= |s|/|s-1| + |s|/sigma):
  sigma in [0.6, 1.3], |t| <= 450:  |zeta| <= 11.89 (|s-1| >= 0.1),  |1/zeta| <= 14.57.   (logs: lemG/box450.log, box_060.log)
- Euler factor A: |A_p - 1| = |z x1 + z x2 - z^2 - x1 x2| / |(1-x1)(1-x2)|; primes <= 100: rigorous sup over the 3-torus of phases
  by cell enclosures (96^3 cells); p > 100: the second-order bound; tail p > 1e6 by n^{-m}.
- Remainders (lemG/remainder.py): E1 (s1 shifted to Re -theta in the box), connectors/tails (Gaussian), R_x extension, K tail,
  cross-term shifts rho_V, rho_U (e^{-Delta b}), w-shift (e^{-theta' Y}). At the binding cell (2.192, 1.732, 5.393, L = 33):
  TOTAL 8.0e-4 (E1a 6.9e-4). Small-U cells are worse (u = 2, L = 29: 0.17) -- the bound is computed PER CELL.
- S_t (lemG/stbound.py): rigorous block sum, S_t <= 2.14 lX eta^2/Delta^2 over lU in [40,120], tau in [1.4,2.4], lX/lV in [1.2,2.5];
  lemmaD uses 2.5 and restricts G to tau in [1.3, 2.5], lX/lV in [1.1, 2.6], lU >= 40 (TODO: widen the S_t grid to cover exactly this).
- k1, k0 for the smoothed kernel: lemG/kcert2.py (derivative under the integral; Euler tails via the L^2 trick). lemmaD uses the
  provisional upper ends k1 <= -0.55, k0 <= 1.0 until the certificate finishes.
- lemmaD.d4_bound_G is now the rigorous form: D4m(k1_up, k0_up) + REM(cell) + 2 sqrt(S_s S_t) + S_t (+ n = 1 adjustment).
  Binding cell: d4 = 1.480 (was 5.756 with Lemma F). Tables: lemmaD_tableGGrig_*.pkl (key version 8).
- CERTIFIED (kcert2.log): k1 = -0.607063 +- 0.0217, k0 = 0.35026 +- 0.0711  =>  k1 <= -0.5853, k0 <= 0.4214.
  lemmaD uses k1_up = -0.55, k0_up = 1.0 (weaker, so still rigorous).

## RESULT (Oct 13): rigorous Lemma G in the chain; RH to 8.0e12 suffices (was 1.75e13). 7e12 fails. See h_Grig.txt.
Every table cell uses the G bound (exact formula, no interpolation), so verify_F is not needed for these tables.
Speed: lemG/rem_cache.pkl stores the cell-independent constants (A sups, Gamma/Gaussian integrals, K_c w-shift integral).
Known soft spots to referee: (1) pole_fac = 40 in the connector term (multiplies e^{-64}); (2) "1e8" sup|K_c| on the w-connectors
(multiplies e^{-63}); (3) mpmath quadrature of explicit majorants (x SAFE = 1.5); (4) S_t constant 2.5 from a grid sup (2.148)
of a smooth function; (5) Psi_t(1) = 1 - rho_s(1) <= Phibar(log U/eta) ignored (~e^{-10^7}).

## Closing the soft spots (Oct 13, evening)
(1) pole_fac -> derived: |zeta(1+u)| <= 1/|u| + 1 + |1+u|/Re(1+u); the 1/|u| singularity integrated over a disk (2 pi R0),
    elsewhere 1/|u| <= 1.  (2) "1e8" -> Kc_bound(W0), valid on the whole w-connector (Re w in [-thp, c]).
(3) All quadratures replaced: |Gamma(sig+iv)| is decreasing in v (sig > -1), J-integrands bounded cell-wise by factor sups,
    the inner K_c integral in closed form pi(a+b)/(ab((a+b)^2+v^2)); SAFE = 1.  (4) S_t now a closed-form rigorous bound per cell
    (stbound.St_rig: exact mass for m <= 1e6, partial summation with |R(x)| <= 4/sqrt x beyond, J_K integrals, TV <= 6 max h).
(5) Psi_t(1) = 1 - rho_s(1) <= Phibar(lU/eta) added to S_t.  k1, k0 now use the certified upper ends (-0.5853, 0.4214).
Binding cell: REM = 5.4e-4, S_t = 9.9e-5, d4 = 1.483.  Tables rebuilt: lemmaD_tableGGrig2_*.pkl (key version 9).
## FINAL (Oct 13): 8.0e12 re-scanned with the fully rigorous bounds: PASS, worst 0.937 (h_Grig_final.txt).

## Push below 8e12 (Oct 13, late)
Lever analysis with G active (binding (60,33), logN 4.199): pair off-diagonal x0.5 -> -0.34; B/10 -> -0.12; J/10 -> -0.06; d5 -> -0.03.
Improvements made:
(a) pair kernel: |Gamma(z)| = |Gamma(1+z)|/|z| with |Gamma(1+s+id)| <= |Gamma(1+2delta+id)|/Gamma(1+2delta) (decreasing in d),
    instead of the flat 1.0067 (kept as the outer factor): coef 2.755 -> 2.435.   (lemmaD._gr)
(b) crowding budget: exact -zeta'/zeta(1+eta) (mpmath, +1e-12) instead of 1/eta - 0.36.   (lemmaD._mzl)
    (a)+(b): binding logN 4.199 -> 4.080 (-0.119).
(c) interval.py zero tail: S(s) decreasing, so on [a,b] use S(a) * int_a^b x^{s-1} ds exactly (was (b-a) x^{b-1}); fine grid
    1/(40 L1) (was 1/(10 L1)); kernel grid c in 40..100 step 5, frac 0.45..0.95 step 0.05.
    (c) alone, old tables, H = 7e12: ratio at L = 1116: 1.040 -> 0.995.
Tables: lemmaD_tableGGrig3_*.pkl (key version 10).

## Oct 13 night: how far down? (tables GGrig3, interval.py v2)
Spot checks (err/tol at L1 = 1050/1100/1150): 6.5e12 -> 0.978/0.980/0.957; 6.0e12 -> 1.041 at 1100. Threshold ~6.3e12 (not yet full-scanned).
FIX: lemmaD.logN returned None for log T < 29, so for H < e^29 = 3.9e12 the zeros in [H, e^29] fell to a crude fallback
(ratio 514 at 3e12!). Now (28, 29] uses the lk = 29 cell (constants at the bin top). TO REFEREE: why was 29 the original cutoff?
3e12 diagnosis: ratio 1.73 at L = 1100; the tail is dominated by 1 - beta in [0.0056, 0.008] (just inside the MTY zero-free
region, R = 5.5587) at heights just above H. Needed: logN lower by ~0.9 in the lk = 29..32 cells (dlogN -0.9 -> 0.971).
Levers at (60, 30): B/10 -0.12, J/10 -0.06, d5 -0.03, d4 is at its true value; pair off-diag x0.5 ~ -0.3. Sum of realistic
levers ~ -0.4: not enough. 3e12 needs a structurally new ingredient (better exponent 2x in T^{2x(1-sigma)}, i.e. shorter
mollifier/Dirichlet polynomial, or a stronger zero-free region near t ~ 1e13: R = 5.0 alone cuts the 3e12 ratio by ~40%).

## Oct 13 night: new zero-free region adopted
Bellotti-Trudgian-Yang, arXiv:2603.21490 (v1, Mar 2026): zeta(s) != 0 for sigma > 1 - 1/(4.896 log t), t >= 3 (unconditional).
R0 (lemmaD) and R (ramanujan_upper) set to 4.896. Spot checks: 3e12 ratio 1.29 / 1.09 / 0.77 at L = 900/1100/1300; 4e12 0.985/0.93/0.70.
log T in (28,29]: uses the lk = 29 cell -- valid since N(s,T) <= N(s,e^29) and the cell's bin is (28,29] by construction.
Tables rebuilt as lemmaD_tableGGrig4 (key 11); run_btY.sh chains full scans at 4.0e12 and 3.5e12 (scanY_*.txt).
Scans with R = 4.896 (tables GGrig4, interval.py with coarse+local kernel search):
  3.5e12 [57,1000): 1.191 at L = 947  FAIL;  4.0e12 [57,1000): 1.063 at L = 947  FAIL;  4.5e12 [57,1000): 0.956 at L = 957  PASS.
  [1000,3159) at 4.0e12 still running (slow at large L; err is monotone decreasing in H, so a pass at 4.0e12 certifies 4.5e12 there).
  Binding now: zeros with 1 - beta just above the new edge (~0.00704) at the lowest heights (lk = 29/30 cells, logN ~ 4.3).
