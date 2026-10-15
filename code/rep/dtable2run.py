# resumable driver: d* on k-steps of 5 (plus every k in the first 0.01 of Lambda above lam*); cells in between take the value
# at the next computed k ABOVE (d* is non-increasing in Lambda: the excluded set grows).  Results -> dstar2_<lk>.json
import json, os, sys, math
from dtable2 import dstar, lam_at
lk = int(sys.argv[1]); fn = 'dstar2_%d.json' % lk
done = json.load(open(fn)) if os.path.exists(fn) else {}
for line in open(sys.argv[2]) if len(sys.argv) > 2 else []:      # import already-computed cells from a log
    p = line.split()
    if len(p) == 4 and int(p[0]) == lk: done['%d,%d' % (lk, int(p[1]))] = float(p[3])
lmin = lam_at(lk); kmin = math.floor(lmin/lk*1e4) + 1; kmax = math.floor(0.0125*(lk-1)/lk*1e4)
ks = [k for k in range(kmin, kmax+1) if k*1e-4*lk <= lmin + 0.01 or (k - kmin) % 5 == 0 or k == kmax]
for k in ks:
    key = '%d,%d' % (lk, k)
    if key in done: continue
    done[key] = dstar(lk, k*1e-4*lk); print(lk, k, done[key], flush=True)
    json.dump(done, open(fn, 'w'))
comp = sorted(int(c.split(',')[1]) for c in done)
out = {}
for k in range(kmin, kmax+1):
    above = [c for c in comp if c >= k]
    out['%d,%d' % (lk, k)] = done['%d,%d' % (lk, min(above))]
json.dump(out, open('dstar2_%d.json' % lk, 'w')); print('done', lk, len(out))
