#!/bin/bash
# Ball-arithmetic certification of the density tables, the scan and the transfer (paper, Section 7.3).  ~1.5 h on 2 cores.
# Needs the tables code/lemmaD_tableGH0.01_*.pkl and code/rep/{lamstar2c,dstar2c_all}.json (code/rerun_all.sh builds them).
set -e
cd "$(dirname "$0")"
python3 gammaint_cert.py                         # I_Gamma <= 5.922 and the J-integral <= 4.72
python3 remconst.py                              # Lemma G remainder constants  -> remconst.json
rm -f cells_arb_*.json scan_arb_*.txt
( python3 cells_arb.py 0 2 ) & ( python3 cells_arb.py 1 2 ) & wait      # every density cell -> cells_arb_{0,1}.json
( python3 scan_arb.py 0 2 ) & ( python3 scan_arb.py 1 2 ) & wait        # the scan          -> scan_arb_{0,1}.txt
python3 ../coverage_check.py scan_arb_0.txt scan_arb_1.txt
python3 transfer_arb.py scan_arb_0.txt scan_arb_1.txt
