# Lemma D -- Bellotti's log-free density estimate, re-run with a log-free mollifier bound and per-sigma parameters

Setting and notation as in Bellotti (arXiv:2405.12545): L = log T, U = T^u, V = T^v (v = t u), W = T^w, X = T^x,
Barban-Vehov weights psi_d (U, V) and theta_d (1, W), Psi(n) = sum_{d|n} psi_d, Theta(n) = sum_{d|n} theta_d,
M(s) = sum_{a<=V, b<=W} psi_a theta_b [a,b]^{-s}, so that sum_n Psi(n)Theta(n) n^{-s} = zeta(s) M(s).

## D1. Log-free pointwise bound for the mollifier on Re s = 1/2
|psi_a|, |theta_b| <= 1 and both are supported on squarefree integers.  Writing a = g a', b = g b' with
(a', b') = 1, [a, b] = g a' b' and
    |M(1/2 + it)| <= sum_{g <= W} g^{-1/2} (sum_{a' <= V/g} a'^{-1/2}) (sum_{b' <= W/g} b'^{-1/2})
                  <= sum_g g^{-1/2} * 2 (V/g)^{1/2} * 2 (W/g)^{1/2}  <  4 zeta(3/2) (VW)^{1/2} = 10.45 (VW)^{1/2}.
In the proof of her Lemma 3.4, S(X) = sum_n Psi Theta n^{-rho} e^{-n/X} = (1/2 pi i) int zeta M(rho + s) Gamma(s) X^s ds;
shifting to Re s = 1/2 - beta picks up only s = 1 - rho (|Gamma(1 - rho)| ~ e^{-pi gamma/2}: nil for gamma >= 3e12;
the residue at s = 0 vanishes since zeta(rho) = 0), and with |zeta(1/2 + iu)| <= 0.618 |u|^{1/6} log|u| (|u| >= 3,
Hiary-Patel-Yang; <= 1.461 for |u| < 3) and (gamma + |t|)^{1/6} log(gamma + |t|) <= T^{1/6} log T (1 + |t|):
    |S(X)| <= X^{1/2 - alpha} / (2 pi) * 10.45 T^{(v+w)/2} * 0.618 T^{1/6} log T * I_Gamma,
    I_Gamma = sup_{sigma in [-1/2, 1/2 - alpha]} int |Gamma(sigma + it)| (1 + |t|) dt = 5.889.
This replaces her 44.66 d2 (v+w)^4 L^4 T^{(v+w)/2 + 27/164 + x(1/2 - alpha_0)} (a factor ~ 6e8 at L = 34 smaller,
i.e. ~ 1.1 less in x for the same c1).

## D2. Parameters per (sigma, T) cell
Her theorem holds for every admissible (u, v, w, x) (w < u, u + v < x, X > U/L^2, c1^2 > c5).  We minimise the
resulting bound separately on every cell delta in ((k-1)e-4, k e-4], log T in (l-1, l] (constants evaluated at the
cell's worst corner), using her formulas for d3, d4 (Lemma 3.3 + (3.6)), d5 (Lemma 4.3), c2, c3, c4, c5 (eqs 4.17-4.21,
J3 <= 286.077 T^{27/164}).  Reconstruction check against her Tables 5-6 at [e^34, e^35] with her parameters:
d4 91.9 (90.96), d5 1.109 (1.110), c5 0.063 (0.066), c1 0.983 (0.984), c2 267.6 (264.0).
This also extends the estimate below her alpha_0 (to sigma >= 0.984), with n~ = floor(Lemma A box count).

## D3. Combined with Lemma A (rho, n~) and Lemma C (c4 crowding).
N(sigma, T) <= n~ (d4 d5 / w)(1.0067 + rho 2.2548 log(1/nu) + 12.03 I') / (c1^2 - c5) * T^{2 x (1-sigma)} / ((1-sigma) log T).

Effect (bound on N(sigma, T)):  sigma = 0.9945: 4.0e4 -> 2.7e4;  0.993: 6.6e4 -> 4.5e4;  0.992 (T = 1e14): 1.5e6 -> 6.2e4;
0.990: 1.5e6 -> 2.3e5;  0.986: 3.6e6 -> 7.7e5.

## What blocks further gains inside this architecture (measured with the optimiser)
 * the c5 "lower" piece (U/L^2)^{3/2 - 2 alpha} W T^{27/164}: forces u >~ 2.07 (w + 0.165) + O(1/L); the 2w comes from the
   trivial bound for the Selberg-weight polynomial G on Re s = 1/2 (no provable cancellation pointwise);
 * x >= u + v (and the B-term): x ~ 1 + 5w + O(1/L); dropping u + v < x would add only ~1.25x;
 * c2 ~ 1/w (log-freeness of Theta).  Net: K e^{2 x lambda} ~ (1/w) e^{2(1 + 5w) lambda}.
A big further gain needs pointwise square-root cancellation for the Theta-mollifier on the critical line
(u ~ w + 0.4 instead of 2w + 0.4 would give x ~ 3.7 instead of ~6.5: ~4x at lambda = 0.25), i.e. a new idea.

## Result (interval.py, every x in [e^58, e^3159], rigorous upper Riemann sums): RH to 7.8e13 suffices
   [58,1300) 0.514  [1300,1800) 0.864  [1800,1950) 0.970  [1950,2100) 0.969  [2100,2150) 0.861
   [2150,2400) 0.819  [2400,3159) 0.649        (7.0e13 fails at L = 1900-2020, worst 1.067)
History: 1.12e14 (Lemma A) -> 1.08e14 -> 9.6e13 (Lemma C) -> 7.8e13 (Lemma D).
