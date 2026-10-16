#!/bin/bash
# Final adaptive scan at H (default 3e12): writes logs/final_<H>_{low,a,b,c,d}.txt.  Run from code/.
HERE=$(cd "$(dirname "$0")" && pwd); cd "$HERE"; export LEMMAG=1; H=${1:-3.0e12}; LOG="$HERE/../logs"
( python3 interval_adapt.py $H 43 57 2 > $LOG/final_${H}_low.txt; python3 interval_adapt.py $H 57 1000 10 > $LOG/final_${H}_a.txt ) &
( python3 interval_adapt.py $H 1000 2000 10 > $LOG/final_${H}_b.txt; python3 interval_adapt.py $H 2000 2600 10 > $LOG/final_${H}_c.txt; python3 interval_adapt.py $H 2600 3159 10 > $LOG/final_${H}_d.txt ) &
wait
tail -n 1 $LOG/final_${H}_*.txt
