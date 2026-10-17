"""Ball-arithmetic assembly (the "far" terms below, from kcert2.py, are upper bounds for the parts
150 <= |Im s| <= 400 and the horizontal pieces of P; kfar_arb.py certifies that these parts are at most 1e-11 and 1.4e-11,
far below the values used).
Ball-arithmetic assembly of the enclosures of k1, k0 (Theorem G) from the output of code/lemG/kcert2.py.
kcert2.py computes, with acb.integral, the integrals I0, I1 over |t| <= 150 on Re s = 0.38 (primes <= 2000 in the Euler
product) and upper bounds Q0, Q1 for int |g|; it then adds the Euler-tail and far-range errors in double precision.
Here those error terms are recomputed in arb.  Q0, Q1 are double-precision upper bounds of sqrt(pi * (rigorous upper bound));
they are inflated by a relative 1e-12, which exceeds the rounding error of one square root and one multiplication."""
from flint import arb, ctx
ctx.prec = 100
I0 = arb('-0.607063', '6.86e-7'); I1 = arb('0.70067', '1.72e-6')          # real parts, from kcert2.log
Q0 = arb('10.84103323353977')*(1 + arb('1e-12')); Q1 = arb('27.926256506182526')*(1 + arb('1e-12'))
eta, T2, kap, T, P = arb('0.02'), arb(400), arb('0.38'), arb(150), arb(2000)
PI = arb.pi(); g = arb.const_euler()
m = 2 - kap
tau = 4*arb('1.3')*arb('1.26')*P**(1 - m)/((m - 1)*P.log())
tau1 = arb('1.1')*arb('1.26')/P*(1 + 1/P.log())*2
e0 = (tau.exp() - 1)*Q0/(2*PI)
e1 = ((tau.exp() - 1)*Q1 + tau1*Q0*tau.exp())/(2*PI)
far = 2*30*arb('14.57')*arb('3.4')*(eta**2*kap**2 - eta**2*T*T).exp()/(3*T**3)/(2*PI)
k1 = I0 + arb(0, (e0 + far).upper())
k0 = I1 + g*I0 + arb(0, (e1 + far*(10 + eta**2*T2) + g*(e0 + far)).upper())
print('k1 in', k1, ' upper', k1.upper(), ' (used: -0.5853)', 'OK' if k1 < arb('-0.5853') else 'FAIL')
print('k0 in', k0, ' upper', k0.upper(), ' (used: 0.4214)', 'OK' if k0 < arb('0.4214') else 'FAIL')
