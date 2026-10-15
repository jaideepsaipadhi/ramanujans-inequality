# certify every computed d* cell (others inherit by monotonicity in Lambda); uses the CERTIFIED lambda* table.
import json, math, sys, os
import rep_cert as RC
RC.lamt = {float(k): v for k, v in json.load(open('lamstar2c.json')).items()}
from multiprocessing import Pool
def ks_for(lk):
    lmin = RC.lam_at(lk); kmin = math.floor(lmin/lk*1e4) + 1; kmax = math.floor(0.0125*(lk-1)/lk*1e4)
    return kmin, kmax, [k for k in range(kmin, kmax+1) if k*1e-4*lk <= lmin + 0.01 or (k - kmin) % 5 == 0 or k == kmax]
def job(args):
    lk, k, d = args
    d = 0.97*d                     # at the float d* the margin is ~0 at the d-edge: certify 3% inside
    while d > 0.004:
        ok, nb = RC.certify(lk, k*1e-4*lk, d, maxdepth=40)
        if ok: return lk, k, d, nb
        d *= 0.9
    return lk, k, 0.0, 0
if __name__ == "__main__":
    lks = [int(x) for x in sys.argv[1].split(',')]
    jobs = []
    for lk in lks:
        old = json.load(open('dstar2_%d.json' % lk))
        kmin, kmax, ks = ks_for(lk)
        for k in ks:
            ka = min([int(c.split(',')[1]) for c in old if int(c.split(',')[1]) >= k] or [None])
            d = old.get('%d,%d' % (lk, k), 0.0) if ka is not None else 0.0
            if d > 0: jobs.append((lk, k, min(d, 1.0)))
    res = {}
    with Pool(2) as P:
        for lk, k, d, nb in P.imap_unordered(job, jobs):
            res['%d,%d' % (lk, k)] = d; print(lk, k, d, nb, flush=True)
            json.dump(res, open('dcert_%s.json' % sys.argv[1].replace(',', '_'), 'w'))
    # expand to all k by monotonicity (value at the next certified k above; 0 if none)
    out = {}
    for lk in lks:
        kmin, kmax, ks = ks_for(lk)
        comp = sorted(int(c.split(',')[1]) for c in res if int(c.split(',')[0]) == lk)
        for k in range(kmin, kmax+1):
            above = [c for c in comp if c >= k]
            out['%d,%d' % (lk, k)] = res['%d,%d' % (lk, min(above))] if above else 0.0
    json.dump(out, open('dstar2c_%s.json' % sys.argv[1].replace(',', '_'), 'w')); print('done')
