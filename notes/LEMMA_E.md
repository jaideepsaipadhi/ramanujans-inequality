# Lemma E and the corrected pair sum (Lemma A'')  --  Oct 9 2026

## Lemma E: the Barban-Vehov mean square via the parametric Ramare-Zuniga Alterman bound
Bellotti bounds ||xi||^2 = sum_n Psi(n)^2 w(n)/n by Lemma 3.3 (= RZA Cor. 1.2, ~ 3.09 log N/log(V/U) F(t)/(t-1))
plus a d3-tail beyond UV, which gives d4 ~ 90.  Numerically (U, V <= 3e4, N <= 2e8) the truth is
~ (log N - log U)/log(V/U) ~ 2.5-3.7, i.e. d4 is ~30-40x too large.  [bv.c]
Lemma E.  With w(n) <= e^{-n/X} and e^{-y} <= (eps/e)^eps y^{-eps} (y > 0),
   sum_n Psi^2 w/n <= w(1) + (eps/e)^eps X^eps ( sum_n Psi(n)^2 n^{-1-eps} - 1 ),
and RZA (arXiv:2405.12662, Theorem 1.1; z1 = U >= 100, z2 = V = U^tau) gives
   sum_n Psi^2 n^{-1-eps} <= e^{gamma eps}/(eps log(V/U)) [A(tau+1) + (B - C U^{-eps} + (B - C V^{-eps}) tau^2) eps log U]/(tau-1),
A = 1.084, B = 1.301, C = 0.318.  eps = a/log X, a optimised.  Result: d4 ~ 20 instead of ~ 90 (still ~10x above the
truth: closing that gap needs an explicit Graham-type asymptotic with small relative error, not available).
Since Lemma E bounds the whole sum, Bellotti's condition u + v < x (needed only to split Lemma 4.2 at UV) is dropped.

## Lemma A'' (replaces Lemma A's pair-sum LP -- correction)
Lemma A's LP took one representative per delta-slot at distance >= delta and left |Delta| < delta to Bellotti's c2
("one representative in this window").  Representatives of ADJACENT rectangles can be closer than delta, so that
accounting is not justified.  Lemma A'' makes no spacing assumption: rectangles of height h = kappa delta, one
representative each; the representative of the k-th rectangle away on either side is at distance d in
((k-1)h, (k+1)h); since X^y/y is increasing on [2 nu, 2 delta], |Gamma(z)| X^{Re z} <= 1.0067 X^{2 delta}/max(|Delta|, 2 delta)
for every pair, so it contributes <= (diagonal) * 2 delta / max(2 delta, (k-1)h); its crowding cost at
s0 = 1 + eta + i gamma_r is >= min_{a in [eta, eta+delta]} a/(a^2 + ((k+1)h)^2); budget (1/2) log T + 1/eta - cost(r).
Zeros per rectangle: n~ = floor(min_eta M / cost(h/2)).  kappa chosen per cell.  This is ~15% weaker than the old LP.
The earlier thresholds 9.6e13 and 7.8e13 used the old LP; they are superseded.  Lemma A's own density override is now
disabled (LEGACY_A=1 restores it); only Bellotti's published rows, KLN and Lemma D (with C, E, A'') are used.

## Result
RH to ~4.6e13 (bisection binding at L ~ 1700-1800; full scans in g_*.txt).
Full scans at H = 4.6e13 (worst error/tolerance): [58,656) 0.514  [656,1300) 0.517  [1300,1600) 0.886
[1600,1900) 0.988 (L = 1730)  [1900,2400) 0.858  [2400,3159) 0.571.   => RH to 4.6e13 suffices.
History: 1.12e14 -> 1.08e14 -> 9.6e13 -> 7.8e13 (old LP) -> 4.6e13 (Lemma E + corrected A'').

## Oct 9 (later): budget refinements and what was tried
 * Crowding budget sharpened (all three uses: Lemma A'' LP, n~, Lemma C windows):
   sum_rho Re 1/(s-rho) <= (1/2) log(t/(2 pi)) + (-zeta'/zeta(1+eta)),  -zeta'/zeta(1+eta) <= 1/eta - 0.36 (0 < eta <= 1.5,
   min of 1/eta + zeta'/zeta(1+eta) on that range is 0.378, computed), and the (c + cos)^2 positivity variant
   (1/2) log(t/2pi) + ((c^2+1/2)/(2c)) (-zeta'/zeta(1+eta)) + log(t/pi)/(8c).  Pair coefficient 2.79 -> 2.67 at the edge.
 * Full scans at H = 4.45e13: [58,656) 0.514  [656,1300) 0.523  [1300,1600) 0.906  [1600,1900) 1.000- (L = 1720, passes,
   no failing interval)  [1900,2400) 0.856  [2400,3159) 0.570.   => RH to 4.45e13 (4.5e13 with margin).
 * Lemma W (weighted Halasz, weights 1/D_r, Schur test): implemented (interval.py WEIGHTED=k), rigorous, but WORSE here
   (1.110 vs 1.090 at H = 4e13, L = 1730): the single parameter set and the pair coefficient at the band top cost more
   than the sigma-integration factor it saves.  Off by default.
 * Finer kernel grids (c, frac): < 0.2% change.
 * Sharper Barban-Vehov mean square: RZA lose a factor ~(tau+1)/(tau-1) in their Step 4 ((a-b)^2 <= a^2 + b^2 for the
   inner sums mcheck_q(z2/q) - mcheck_q(z1/q), q <= z1).  With an explicit two-sided bound for the coprime Moebius
   log-sums (so that the cross term costs ~eta/(tau-1), eta ~ 0.3) the density drops ~3.4x (optimizer, [hyp.py]),
   i.e. H ~ 2.5e13.  Available explicit inputs (Granville-Ramare |sum_{d<=x,(d,q)=1} mu(d)/d| <= 1; Ramare
   <= c (q/phi(q))/log(x/q), x > q) are NOT enough: the t^{-eps} partial summation leaves an error ~ eps log(z2/z1)
   log(q) that dwarfs the main term.  Needed: |sum_{l<=t,(l,q)=1} mu(l)/l| <= c (q/phi(q))/log^2 t-type decay uniformly
   in q, or a mean-square-over-q version.  This is the single biggest remaining lever.
