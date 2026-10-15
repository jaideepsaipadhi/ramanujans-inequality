# Rigorous Lemma 5 (per sigma, exact delta_m = (2 sigma - 1) m): with y' = y/(2 sigma - 1) and eps0 = Cw eta0^2/(w0 (2sigma-1)^2),
#   [(2 sigma-1)/(eta w0)] sum_{m>=1} b_m Re W((delta_m + i y)/eta) >= B(y') := sum_m b_m (m - sgn(b_m) eps0/m)/(m^2 + y'^2)
# (|W0(z)| <= Cw |z|^-3 and |z|^-3 <= eta^3/((2sigma-1) m |delta_m + iy|^2)).  Certify B > 0 on R for eps0 = 1/500.
import sympy as sp
from fractions import Fraction as Fr
from lemma5 import b
y = sp.symbols('y', real=True); eps = sp.Rational(1, 500)
B = sum(sp.Rational(b[m].numerator, b[m].denominator)*(m - (1 if b[m] > 0 else -1)*eps/m)/(m*m + y**2) for m in range(1, 8))
num, den = sp.fraction(sp.together(B))
num = sp.Poly(sp.expand(num), y); den = sp.Poly(sp.expand(den), y)
print('den real roots:', len(sp.real_roots(den)), ' num real roots:', num.count_roots(), ' B(0) =', float(B.subs(y, 0)))
print('=> B > 0 on R' if num.count_roots() == 0 and len(sp.real_roots(den)) == 0 and B.subs(y, 0) > 0 else 'FAIL')
