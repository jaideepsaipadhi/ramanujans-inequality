import json, math, numpy as np
from rep import *
lamt = {float(k): v for k, v in json.load(open('lamstar.json')).items()}
out = {}
for lk in range(29, 38):
    lmin = min(v for k, v in lamt.items() if k >= lk)
    for k in range(1, 400):
        de = k/1e4; Lam = de*lk
        if Lam <= lmin: continue
        if Lam > 0.24: break                      # Lemma 11 validity (eta <= 0.0085)
        d = dstar(float(lk), Lam, ngrid=5)
        out['%d,%d' % (lk, k)] = d
        print(lk, k, round(Lam, 4), d, flush=True)
json.dump(out, open('dstar.json', 'w'))
