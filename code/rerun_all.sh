#!/bin/bash
# Full re-certification after the Oct 15 audit fixes, then tables and the final scan.  ~8-10 h on 2 cores.
set -e
cd /home/claude/rep
python3 f_nonneg_cert.py | tail -1
python3 mp_cert.py | tail -1
python3 lemma5_cert.py | tail -1
python3 misc_cert.py | tail -1
rm -f lamstar2c.json dcert_*.json
python3 lamcert.py 4                                   # Lemma H chain     -> lamstar2c.json      (~30 min)
python3 repcert_run.py 29,30,31,32,33,34,35,36,37      # Lemma Rep cells   -> dstar2c_*.json      (~3 h)
python3 min_fix.py                                     # min over lower height cells -> dstar2c_all.json
cd /home/claude
mkdir -p tabGH_preaudit && mv lemmaD_tableGH0.01_*.pkl tabGH_preaudit/ 2>/dev/null || true
export LEMMAG=1
( python3 bt4.py 29,29.5,30,30.5,31 ; python3 bt4.py 32,33,34,35 ) &
( python3 bt4.py 36,37,38,39,40 ) &
wait                                                   # Lemma D tables (~3 h)
./run_final.sh 3.0e12                                  # adaptive scan, prints the three summary lines (~4 h)
