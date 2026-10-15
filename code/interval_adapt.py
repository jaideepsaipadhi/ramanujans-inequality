"""Adaptive cover of [lo, hi): each interval [L1, L2] (initial width w) is certified by interval.err_interval; if the ratio is
>= 1 it is split into two halves (down to width 0.5).  Every x with log x in [lo, hi) lies in a certified interval."""
import sys, interval
H = float(sys.argv[1]); lo = float(sys.argv[2]); hi = float(sys.argv[3]); w = float(sys.argv[4])
worst = (0, None); fails = []; log = []
def cert(a, b):
    q = interval.err_interval(a, b, H)
    if q < 1: log.append((a, b, q)); return q
    if b - a <= 0.5: fails.append((a, b, q)); return q
    m = (a + b)/2; return max(cert(a, m), cert(m, b))
L = lo
while L < hi:
    b = min(L + w, hi); q = cert(L, b)
    if q > worst[0]: worst = (q, L)
    print('%g %g %.4f' % (L, b, q), flush=True)
    L = b
print('H=%.3e [%g,%g) adaptive (start width %g): worst certified-cover ratio %.4f at L=%g; uncovered: %s; subdivided: %s' %
      (H, lo, hi, w, worst[0], worst[1], fails, sorted(set((a, b) for a, b, q in log if b - a < w))[:30]), flush=True)
