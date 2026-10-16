# Ramanujan's prime-counting inequality, unconditionally for all x ≥ 38,358,837,683

**Theorem.** For every real x ≥ 38 358 837 683,

    π(x)² < (e x / log x) · π(x/e).

The bound is unconditional: the only verified-RH input is Platt–Trudgian's verification up to height 3·10¹²,
and no new zero computation is used.

The theorem combines three ranges:
- 38 358 837 683 ≤ x ≤ e¹⁰³, proved by Johnston;
- e⁵⁸ ≤ x ≤ e³¹⁵⁹, proved here (computer-assisted, certified);
- x ≥ e^3158.442, proved by Axler.

Galway and Dudek–Platt computed that the last failure is below 38 358 837 683, which gives N_R = 38 358 837 683.
Dudek and Platt had previously obtained this value under the Riemann hypothesis.

The proof is in [`paper/ramanujan.pdf`](paper/ramanujan.pdf).

## Layout
| path | contents |
|---|---|
| `paper/` | LaTeX source and PDF of the proof |
| `code/` | the certified pipeline: explicit formula with Logan kernel (`interval.py`, `logan.py`), zero-density estimates (`lemmaD.py` + `lemmaC.py`, `pairsum.py`, `bell_patch.py`), Lemma G mean square (`lemG/`), Lemma D tables (`lemmaD_tableGH0.01_*.pkl`) |
| `code/rep/` | Lemma H (zero-free region) and Lemma Rep (zero repulsion), certified in ball arithmetic (`lamcert.py`, `rep_cert.py`, `repcert_run.py`, constant certificates `*_cert.py`), results `lamstar2c.json`, `dstar2c_all.json` |
| `certify/` | the smoothing constant (`smoothing_const.py`) and the ψ → π − li transfer lemma (`transfer_check.py`) |
| `logs/` | final scan logs at H = 3·10¹² and the full re-certification log |
| `notes/` | working notes for every lemma (historical; the paper is authoritative) |

## Reproducing
Requires Python 3 with mpmath, numpy, scipy, sympy, python-flint.

```
cd code && ./rerun_all.sh          # every certificate, the tables, the scan and the final checks (~20 h on 2 cores)
cd certify && python3 coverage_check.py ../logs/final_3.0e12_*.txt && python3 transfer_check.py ../logs/final_3.0e12_*.txt
```

## Final numbers (H = 3·10¹²)
- Scan ratio q < 1 on all of [57, 3159). The worst value is 0.9947, at L = 867; the range [857, 917) needed width-5 intervals.
- Transfer to π − li: the worst ratio is 0.99383, at log y = 877.
- The final scan runs in double precision; the paper (§7) gives the rounding analysis. All other steps are certified in ball arithmetic.
- The scan logs are split across `logs/final_3.0e12_{low,a,b,c,c2,c3,c4,dpar}.txt` because a container restart interrupted the run. The logs together cover the whole range, and `coverage_check.py` verifies this.
