"""Constrained pair sum: zeros near the 1-line can't crowd.
For s = 1+eta+i*gamma_r:  sum_rho Re 1/(s-rho) <= M = 0.5*log(T) + 1/eta  (all terms >= 0).
A zero with beta>=1-delta at distance d costs >= min_{a in [eta,eta+delta]} a/(a^2+d^2).
LP-maximise the c3 pair sum  sum 1/(n delta)  (one representative per delta-slot, both sides)."""
import math
def cost(d, eta, de):
    return min(a/(a*a+d*d) for a in (eta, eta+de))
def lp(de, logT, eta):
    M = 0.5*logT + 1/eta - cost(0, eta, de)      # r itself uses some budget
    items = []
    for n in range(1, int(1/de)):
        b = 1/(n*de); c = cost((n+1)*de, eta, de)
        items += [(b/c, b, c)]*4                  # both sides, up to 2 reps per shifted window
    items.sort(reverse=True); tot = 0.0
    for _, b, c in items:
        if M <= 0: break
        f = min(1.0, M/c); tot += f*b; M -= f*c
    return tot
def box(de, logT):
    best = 9e9
    for k in range(1, 400):
        eta = k*0.0005
        M = 0.5*logT + 1/eta; c = cost(de, eta, de)   # Bellotti's square |t-T|<=delta
        best = min(best, M/c)
    return best
if __name__=="__main__":
  for de in [0.005, 0.006, 0.007, 0.008, 0.01]:
    for logT in [32.0, 35.0]:
        full = 2*sum(1/(n*de) for n in range(1, int(1/de)))
        best = min(lp(de, logT, k*0.001) for k in range(1, 300))
        print(f"delta={de} logT={logT}: c3 sum ratio = {best/full:.3f}   box count <= {box(de,logT):.2f}")
