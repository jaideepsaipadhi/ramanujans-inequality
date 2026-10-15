# Ramanujan's inequality: unconditional with H = 3.0e12 (Platt-Trudgian; no new RH computation)  -- Oct 15, 2026

Final scan (interval_adapt.py, run_final.sh): every x with log x in [57, 3159) is covered by intervals [L1, L2] certified by
interval.err_interval with ratio < 1 (start width 10, bisected where needed).
  [57, 1000):   worst 0.9935 at L = 917; [857, 917) needed width 5
  [1000, 2000): worst 0.9041 at L = 1000
  [2000, 2820): worst 0.4066 (final_3.0e12_c.txt, used up to 2820)
  [2800, 3159): worst 0.4001 (final_3.0e12_d.txt)
  No uncovered intervals.
(At H = 3.1e12 with width 10: [57,1000) worst 0.991; superseded.)

New ingredients since the 8e12 certificate (details: LEMMA_REP.md, LEMMA_G.md):
 1. Lemma H: height-dependent zero-free region from BTY's method evaluated exactly, certified 4-step chain:
    (1 - beta) log g >= 0.21384 at log g <= 29 (R = 4.676), 0.21345 at 30, ... (rep/lamstar2c.json).
 2. Lemma Rep: two near-edge zeros at the same height are >= d* apart (d* ~ 0.058 at log T = 29 down to 0.045 at 37,
    for (1 - beta) log T up to 0.0125 (log T - 1)); certified in ball arithmetic over continuous boxes (rep/dstar2c_*.json);
    used in Lemma D's pair sum, with the minimum over all lower height cells.
 3. BTY constants re-derived for general eta0 with explicit, certified pieces (W0 bound by 4 IBPs, Lemma 5 per-sigma with exact
    root isolation, Lemma 11 in ball arithmetic, Lemma 10 m > 15 tail bounded, monotonicity of G_p proved in one line).
 4. Lemma D: J-integral with HPY 0.618 t^{1/6} log t (all of (4.16) <= 4.72 Z(T)); half-unit height cells below e^31.
STATUS (Oct 16): RE-CERTIFIED after the audit fixes (rerun_all.sh): Lemma H identical, all 131 Rep cells identical d*,
tables rebuilt, scan: [57,1000) worst 0.9966 at L=907; [1000,2000) 0.8884; [2000,2810) 0.4063; [2800,3159) 0.4001; nothing uncovered.
Remarks for the preprint (BTY arXiv:2603.21490): (a) f(0) = kappa eta w(0), so their Lemma-10 term is overstated ~2x (Thm 1
still holds: lambda = 0.2054 > 1/4.896 at e^76.47); (b) their (31) U(y) integrates to 28.3, not < 24 (true digamma: 23.59,
so (32) is true); (c) Lemma 5's uniform check needs sigma >= 0.99346 (fine for them, fails for smaller sigma; per-sigma version
holds up to eta0 ~ 0.0147); (d) the numerical dG_p/dsigma < 0 check is unnecessary (P >= 0).
(e) BTY Table 1, p = 5: true min 0.00735649 < their 0.00735709 (harmless).
