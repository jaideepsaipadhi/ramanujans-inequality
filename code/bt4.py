# resumable lemmaD table builder: python3 bt3.py 29,30,31  -> lemmaD_tableG_29_30_31.pkl (saved every 20 cells)
import sys, os, pickle; sys.path.insert(0, '/home/claude')
import lemmaD as D
lks = [(int(a) if float(a).is_integer() else float(a)) for a in sys.argv[1].split(',')]
fn = '/home/claude/lemmaD_table' + (('GH%s' % os.environ.get('LEMMAG_ERR','0.01')) if D.USE_G else 'G') + '_%s.pkl' % sys.argv[1].replace(',', '_')
tab = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
n = 0
for lk in lks:
    for k in range(50, 301):
        key = (k, lk, D.USE_E, 13 if D.USE_G else 6)
        if key in tab: continue
        tab[key] = D.best(k, lk); n += 1
        if n % 20 == 0: pickle.dump(tab, open(fn + '.tmp', 'wb')); os.replace(fn + '.tmp', fn)
pickle.dump(tab, open(fn + '.tmp', 'wb')); os.replace(fn + '.tmp', fn)
print('done', sys.argv[1], len(tab))
