#!/bin/bash
cd "$(dirname "$0")"
( python3 cells_arb.py 0 2 > cells_arb_0.log 2>&1 ) & ( python3 cells_arb.py 1 2 > cells_arb_1.log 2>&1 ) & wait
( python3 scan_arb.py 0 2 > scan_arb_0.log 2>&1 ) & ( python3 scan_arb.py 1 2 > scan_arb_1.log 2>&1 ) & wait
python3 ../coverage_check.py scan_arb_0.txt scan_arb_1.txt > coverage_arb.log 2>&1
python3 transfer_arb.py scan_arb_0.txt scan_arb_1.txt > transfer_arb.log 2>&1
echo CHAIN_DONE >> transfer_arb.log
