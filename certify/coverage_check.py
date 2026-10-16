"""Checks that the certified intervals (q < 1) in the given scan logs cover [57, 3159), and that [43, 57) is covered by some bound."""
import sys
iv = []
for f in sys.argv[1:]:
    for line in open(f):
        p = line.split()
        if len(p) == 3:
            try: iv.append((float(p[0]), float(p[1]), float(p[2])))
            except ValueError: pass
def covered(lo, hi, need_cert):
    segs = sorted((a, b) for a, b, q in iv if (q < 1 or not need_cert) and b > lo and a < hi)
    x = lo
    for a, b in segs:
        if a > x: return False, x
        x = max(x, b)
    return x >= hi, x
print('[43,57) bounded:', covered(43, 57, False))
ok, x = covered(57, 3159, True); print('[57,3159) certified (q<1):', ok, x)
print('max q on [57,3159):', max(q for a, b, q in iv if b > 57))
