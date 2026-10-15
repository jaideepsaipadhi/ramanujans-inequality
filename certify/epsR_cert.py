"""Rigorous enclosures of eps_R(L) (Definition 2.1) on every cell [a, a + h] of [57.9, 3160.1].
A(L) := L e^{-L} li(e^L).  Then g = A(L)^2 - (L/(L-1)) A(L-1),  b = 2A(L) + L/(L-1),  eps_R = -2g/(b + sqrt(b^2 - 4g)).
Expansion (K integrations by parts of int_a^x dt/log t with a = sqrt x):
  A(L) = P(u) + R(L),  P(u) = sum_{k<K} k! u^k,  u = 1/L,
  R(L) = (L/x)[ li(a) - a P(2/L)/(L/2) ] + (L/x) K! int_a^x dt/log^{K+1} t,
  |li(a) - a P(2u)/(L/2)| <= 3 a/(L/2)   (A(l) <= 1.1 and 0 < P(1/l) <= 1.1 for l >= 28.5),
  0 <= (L/x) K! int_a^x dt/log^{K+1} t <= K! u^K / (1 - 2(K+1)/L).
The cancelling part S(u) = P(u)^2 - P(u/(1-u))/(1-u) is expanded as an exact power series (rational coefficients) to order N
with tail <= (sum_{k<K} k!) (2u)^{N+1}/(1-2u); it is evaluated on u-balls as u^4 * (series / u^4) (s_0..s_3 = 0, checked)."""
from fractions import Fraction as Fr
from math import comb, factorial
from flint import arb, ctx
ctx.prec = 200
K, N = 16, 70
P = [Fr(factorial(k)) for k in range(K)]
P2 = [sum(P[i]*P[n-i] for i in range(max(0, n-K+1), min(n, K-1)+1)) for n in range(N+1)]
T = [sum(factorial(k)*comb(n, k) for k in range(0, min(n, K-1)+1)) for n in range(N+1)]
S = [P2[n] - T[n] for n in range(N+1)]
assert S[0] == S[1] == S[2] == S[3] == 0, S[:5]
print('S series: s4..s7 =', [str(x) for x in S[4:8]])
SK = sum(factorial(k) for k in range(K))
def Sball(u):
    poly = arb(0)
    for n in range(N, 3, -1): poly = poly*u + arb(S[n].numerator)/S[n].denominator
    tail = arb(SK)*(2*u)**(N+1)/(1 - 2*u)
    return u**4*poly + arb(0, tail.upper())
def Pball(u):
    s = arb(0)
    for k in range(K-1, -1, -1): s = s*u + factorial(k)
    return s
def Rball(L):
    u = 1/L
    pos = arb(factorial(K))*u**K/(1 - 2*(K + 1)*u)
    neg = 3*(-L/2).exp()*L/(L/2)          # (L/x) * 3 a/(L/2),  a = e^{L/2}
    lo = -neg; hi = pos + neg
    return (lo + hi)/2 + arb(0, ((hi - lo)/2).upper())
def eps_ball(La, Lb):
    L = arb((La + Lb)/2, (Lb - La)/2)
    u = 1/L
    R1 = Rball(L); A1 = Pball(u) + R1; A2r = Rball(L - 1)
    # g = S(u) + 2 P(u) R1 + R1^2 - R2/(1-u)   [ P(u/(1-u))/(1-u) is inside S ]
    g = Sball(u) + 2*Pball(u)*R1 + R1*R1 - A2r/(1 - u)
    b = 2*A1 + 1/(1 - u)
    return -2*g/(b + (b*b - 4*g).sqrt()), g
if __name__ == "__main__":
    import mpmath as mp, json
    mp.mp.dps = 60
    def eps_mp(L):
        L = mp.mpf(L); x = mp.e**L; li_x, li_xe = mp.li(x), mp.li(x/mp.e)
        g = (li_x**2 - (mp.e*x/L)*li_xe)/(x/L)**2; b = 2*li_x/(x/L) + L/(L - 1)
        return (-b + mp.sqrt(b*b - 4*g))/2
    for L in [58, 100, 907, 3159]:
        e, g = eps_ball(L, L)
        print('L=%d  ball %s   mpmath %s' % (L, e, mp.nstr(eps_mp(L), 15)))
    h = 0.05; a = 57.9; out = []; worst_rel = 0; gpos = False
    while a < 3160.1:
        e, g = eps_ball(a, a + h)
        if not (g < 0): gpos = True
        out.append((round(a, 4), float(e.lower())))
        worst_rel = max(worst_rel, float((e.rad()/e.mid())))
        a += h
    json.dump(out, open('epsR_lower.json', 'w'))
    print('cells:', len(out), ' max relative width of enclosure:', worst_rel, ' g<0 everywhere:', not gpos)
