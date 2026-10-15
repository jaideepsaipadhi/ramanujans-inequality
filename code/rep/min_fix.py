# replace d*(lk,k) by min over height cells l <= lk (zero-free cells, k e-4 l <= lambda*(l), impose nothing); writes all files
import json, glob
import repcert_run as R
D = {}
for f in sorted(glob.glob('dstar2c_*.json')): D.update(json.load(open(f)))
lam = {float(k): v for k, v in json.load(open('lamstar2c.json')).items()}
lam_at = lambda l: lam[min(k for k in lam if k >= l - 1e-12)]
out = {}
for lk in range(29, 38):
    kmin, kmax, _ = R.ks_for(lk)
    for k in range(kmin, kmax + 1):
        out['%d,%d' % (lk, k)] = min(D.get('%d,%d' % (l, k), 0.0) for l in range(29, lk + 1) if k*1e-4*l > lam_at(l))
for f in glob.glob('dstar2c_*.json'): __import__('os').remove(f)
json.dump(out, open('dstar2c_all.json', 'w')); print('min-fixed cells:', len(out))
