#!/bin/bash
cd /home/claude; export LEMMAG=1
H=${1:-3.0e12}
( python3 interval_adapt.py $H 57 1000 10 > final_${H}_a.txt 2>&1; python3 interval_adapt.py $H 2000 3159 10 > final_${H}_c.txt 2>&1 ) &
( python3 interval_adapt.py $H 1000 2000 10 > final_${H}_b.txt 2>&1 ) &
wait
tail -n 1 final_${H}_*.txt
