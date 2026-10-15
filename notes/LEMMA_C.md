# Lemma C -- crowding bound for the far off-diagonal term (c4) in Bellotti's log-free density estimate

Setting (Bellotti, arXiv:2405.12545, Sec. 4): J representatives rho_r (one zero with beta >= alpha = 1 - delta per
delta-box), Halasz-Montgomery pair sum.  For |gamma_r - gamma_s| > 1 she bounds the representatives with
|gamma_r - gamma_s| in [y, y + dy] by dy/delta, which gives
    c4 = 2 sqrt(2 pi) 2.4 (d4 d5 / w) I,   I = int_1^T g(y) dy,   g(y) = e^{y/6} y^{1.5 - 2 alpha} e^{-pi y/2}.

Lemma C.  For a window of height w centred at t <= T and 0 < eta,
    #{rho : beta >= 1 - delta, |gamma - t| <= w/2}  <=  n(w) := min( floor(M / m), ceil(w/delta) + 1 ),
    M = (1/2) log T + 1/eta,   m = min(f(eta), f(eta + delta)),   f(a) = a / (a^2 + (w/2)^2).
Proof.  s = 1 + eta + i t.  Hadamard + Stirling: sum_rho Re 1/(s - rho) = Re zeta'/zeta(s) + Re 1/(s-1)
+ (1/2) Re psi(s/2 + 1) - (1/2) log pi <= 1/eta + (1/2) log(t/2pi) + O(1/t) <= M, using |zeta'/zeta(s)| <=
-zeta'/zeta(1 + eta) < 1/eta.  Every term is positive (Re(s - rho) = 1 + eta - beta > 0) and a zero in the window
with beta >= 1 - delta has a = 1 + eta - beta in [eta, eta + delta], |gamma - t| <= w/2, so its term is >= m
(f is unimodal in a).  The second entry: representatives lie in distinct delta-boxes.  []

Consequence.  g is decreasing on [1, oo) (d/dy log g = 1/6 + (2 delta - 1/2)/y - pi/2 < 0), so tiling
[1, oo) by windows [y_i, y_i + w] on each side,
    I  ->  I' = delta * sum_i n(w) g(y_i)            (terms beyond y = 60 are < 1e-30).
With w = 0.25, delta in [0.0055, 0.0095], log T in [32, 35]: I'/I = 0.14 - 0.28, i.e. c4 = 447 -> 61 - 126.
Implemented in lemmaC.py (worst case over each 1e-4 delta bin, eta optimised on a grid), used by lemmaA.py.

Effect on the Ramanujan threshold (interval.py, rigorous upper Riemann sums, every x in [e^58, e^3159]):
RH to height 9.6e13 now suffices (was 1.08e14 with Lemma A alone; 1.12e14 before that).
    [58,656) 0.514  [656,1300) 0.500  [1300,1900) <=0.952  [1900,2030) 0.983  [2030,2150) 0.989
    [2150,2400) <=0.950  [2400,3159) <=0.705     (worst ratio error/tolerance; last three checked at H = 9.0e13)

Things tested that do NOT help (so nobody re-tries them):
 * re-optimising Bellotti's (u, v, w, x): K and her C share the factor d4 d5 / w, so her optimum is ours.
 * one-representative-per-box LP for c3, and 3-4-1 / (c + cos)^2 budgets in place of 1/eta: rho unchanged
   (the LP optimum sits at eta ~ 0.1-0.3, where the log T term dominates).
 * kernels: super-Gaussian / flat-top Fourier profiles need larger H*m1 than Gaussian/Logan at every relevant
   damping level (p = 1: 1.10, p = 2: 1.92, p = 3: 2.88 at khat(H) = 0.39); exactly band-limited kernels need
   H*m1 >= ~2.24 (vs 1.1 now): density bounds are worth ~2x over "no off-line information".
 * dropping smoothing (truncated explicit formula): undamped tail 1.5-2.3 x tolerance at any H <= 1e15.
 * replacing Brun-Titchmarsh in the smoothing step by zeros: circular (the cost integral contains psi(x)/2),
   or worse via a smooth majorant (undamped zero sum ~ tolerance).

## Round 2 of dead ends (Oct 8): attempts at a large reduction of H
 * Known lower heights in the zero-free-region argument (trig polynomial with target frequency k*theta, frequencies
   < k at heights < H computed exactly): optimised degree <= 8 polynomials give NO gain at leading order (the optimum
   puts zero weight on the known frequencies; their terms are O(1), not O(log T)).  [trig.py]
 * Smaller Barban-Vehov lengths (u, v ~ w) to cut x (B = 2x): blocked by the d3 (x - u - v) mean-square term of Psi,
   which forces u + v ~ x; with x(beta - 1/2) > 1/6 + (v+w)/2 this forces u ~ w + 0.33 + ... -> Bellotti's u ~ 1.9.
   Removing the L^4 loss in her contour term would cut x ~6 -> ~4.7: density / ~1.6, H / ~1.3.
 * Strip test (no zeros with beta >= 0.97) instead of RH: certified evaluation of anything about zeta on [T, T + D]
   costs ~ bandwidth (log sqrt t) x D samples, the same order as RH verification, and the strip must run to the same
   height (Logan damping still needs H*m1 ~ 1.1).  No gain.
 * Weil/Guinand-type localized tests with prime sums: sources = primes up to e^U per block (>> sqrt(t) RS terms);
   more expensive than Riemann-Siegel for any useful detection power.
 * Asymmetric kernels: with Brun-Titchmarsh constant 2 the symmetric kernel is a strict local minimum of the cost.
