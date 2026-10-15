import json, sys, math
import repcert_run as R
res = json.load(open('dcert_29_30_31_32_33_34_35_36_37.json'))
for lk in [int(x) for x in sys.argv[1].split(',')]:
    kmin, kmax, ks = R.ks_for(lk)
    assert all('%d,%d' % (lk, k) in res for k in ks), 'lk %d not complete' % lk
    comp = sorted(int(c.split(',')[1]) for c in res if int(c.split(',')[0]) == lk)
    out = {}
    for k in range(kmin, kmax+1):
        above = [c for c in comp if c >= k]
        out['%d,%d' % (lk, k)] = res['%d,%d' % (lk, min(above))] if above else 0.0
    json.dump(out, open('dstar2c_%d.json' % lk, 'w')); print(lk, len(out), min(out.values()))
