#!/bin/bash
# Full re-certification and final scan; regenerates every log cited in the paper.  Run from code/.  ~20 h on 2 cores.
set -e
HERE=$(cd "$(dirname "$0")" && pwd); cd "$HERE/rep"
python3 f_nonneg_cert.py | tail -1; python3 mp_cert.py | tail -1; python3 lemma5_cert.py | tail -1; python3 misc_cert.py | tail -1
rm -f lamstar2c.json dcert_*.json
python3 lamcert.py 4                                   # Lemma H chain     -> lamstar2c.json
python3 repcert_run.py 29,30,31,32,33,34,35,36,37      # Lemma Rep cells   -> dcert_*.json
python3 min_fix.py                                     # min over lower height cells -> dstar2c_all.json
cd "$HERE"; export LEMMAG=1
( python3 bt4.py 29,29.5,30,30.5,31 ; python3 bt4.py 32,33,34,35 ) & ( python3 bt4.py 36,37,38,39,40 ) & wait
./run_final.sh 3.0e12
cd "$HERE/../certify" && python3 coverage_check.py ../logs/final_3.0e12_*.txt && python3 smoothing_const.py | tail -1 && python3 epsR_cert.py | tail -1 && python3 transfer_check.py ../logs/final_3.0e12_*.txt
cd "$HERE/../certify/arb" && ./run_all.sh               # ball-arithmetic density cells, scan and transfer (Section 7.3)
